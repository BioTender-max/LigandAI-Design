#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ligandai_scoring_filter.py
==========================
LigandAI 多维度评分筛选与帕累托前沿分析示例
Multi-Dimensional Scoring, Filtering & Pareto Frontier Analysis

功能概述 / Feature Overview:
  1. 加载已生成的候选肽段数据（CSV 或模拟数据）
  2. 多维度加权评分函数（delta_g / ipsae / immunogenicity / solubility）
  3. 多级筛选过滤（硬阈值 + 软评分）
  4. 帕累托前沿分析（delta_g vs immunogenicity 双目标优化）
  5. matplotlib 可视化（散点图 + 帕累托前沿线 + Elite Hits 标注）
  6. 输出筛选报告（控制台 + Markdown 文件）

依赖 / Requirements:
  pip install pandas numpy matplotlib rich

作者 / Author: LigandAI Example Suite
版本 / Version: 1.0.0
"""

import os
import sys
import logging
import warnings
from typing import List, Tuple, Dict, Any, Optional
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # 无头模式，适合服务器环境 / Headless mode for server environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.rule import Rule

warnings.filterwarnings("ignore", category=FutureWarning)

# ── 日志配置 ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("ligandai.scoring_filter")

console = Console()

# ── 全局配置 ───────────────────────────────────────────────────────────────────
# Global configuration — adjust weights and thresholds to your project needs

# 加权评分权重（总和应为 1.0）/ Scoring weights (should sum to 1.0)
SCORE_WEIGHTS: Dict[str, float] = {
    "delta_g":        0.40,   # 结合自由能权重（最重要）/ Binding free energy weight
    "ipsae":          0.25,   # 结构对齐误差权重 / Structure alignment error weight
    "immunogenicity": 0.20,   # 免疫原性权重 / Immunogenicity weight
    "solubility":     0.15,   # 溶解度权重 / Solubility weight
}

# 硬阈值过滤条件 / Hard-threshold filter criteria
HARD_FILTERS: Dict[str, Any] = {
    "delta_g_max":        -5.0,   # ΔG 必须低于此值（kcal/mol）
    "ipsae_max":           5.0,   # ipSAE 必须低于此值
    "immunogenicity_max":  0.6,   # 免疫原性必须低于此值
    "solubility_min":      0.3,   # 溶解度必须高于此值
    "length_min":          8,     # 最短肽段长度
    "length_max":          20,    # 最长肽段长度
}

# Elite Hit 阈值（帕累托前沿中的精英候选）
ELITE_DELTA_G_THRESHOLD    = -15.0   # ΔG < -15 kcal/mol
ELITE_IMMUNOGENICITY_LIMIT =  0.3    # 免疫原性 < 0.3

OUTPUT_PLOT   = "ligandai_pareto_analysis.png"
OUTPUT_REPORT = "ligandai_filter_report.md"
INPUT_CSV     = "ligandai_top_candidates.csv"   # 来自 e2e_pipeline.py 的输出


# ══════════════════════════════════════════════════════════════════════════════
# 数据加载 / Data Loading
# ══════════════════════════════════════════════════════════════════════════════

def load_or_simulate_data(csv_path: str) -> pd.DataFrame:
    """
    加载候选肽段数据：优先读取 CSV，若不存在则生成模拟数据用于演示。
    Load candidate peptide data: read CSV if available, else generate mock data.

    Args:
        csv_path: CSV 文件路径（来自 e2e_pipeline.py 输出）

    Returns:
        候选肽段 DataFrame
    """
    if Path(csv_path).exists():
        df = pd.read_csv(csv_path, encoding="utf-8-sig")
        logger.info("从 CSV 加载数据: %s (%d 条)", csv_path, len(df))
        console.print(f"  [green]✓[/green] 从文件加载: [bold]{csv_path}[/bold]  ({len(df)} 条记录)")
    else:
        logger.warning("CSV 文件不存在，生成模拟数据用于演示")
        console.print(f"  [yellow]⚠[/yellow] 未找到 {csv_path}，使用模拟数据演示")
        df = _generate_mock_data(n=120)

    # 确保必要列存在 / Ensure required columns exist
    required_cols = ["sequence", "delta_g", "ipsae", "immunogenicity", "solubility"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"数据缺少必要列 / Missing required columns: {missing}")

    # 补充 length 列（若缺失）
    if "length" not in df.columns:
        df["length"] = df["sequence"].str.len()

    return df


def _generate_mock_data(n: int = 120, seed: int = 42) -> pd.DataFrame:
    """
    生成模拟候选肽段数据（仅用于演示，不代表真实预测结果）。
    Generate mock peptide candidate data for demonstration purposes only.
    """
    rng = np.random.default_rng(seed)

    # 氨基酸字母表
    aa = list("ACDEFGHIKLMNPQRSTVWY")

    sequences = []
    for _ in range(n):
        length = rng.integers(8, 21)
        seq = "".join(rng.choice(aa, size=length))
        sequences.append(seq)

    df = pd.DataFrame({
        "sequence":       sequences,
        "length":         [len(s) for s in sequences],
        # ΔG: 正态分布，均值 -10，标准差 6（kcal/mol）
        "delta_g":        rng.normal(-10.0, 6.0, n).round(3),
        # ipSAE: 对数正态分布，模拟 0-10 范围
        "ipsae":          np.clip(rng.lognormal(1.0, 0.6, n), 0.5, 9.5).round(3),
        # 免疫原性: Beta 分布，偏向低值
        "immunogenicity": rng.beta(2, 5, n).round(3),
        # 溶解度: Beta 分布，偏向高值
        "solubility":     rng.beta(5, 2, n).round(3),
        # 稳定性: 均匀分布
        "stability":      rng.uniform(0.3, 0.95, n).round(3),
    })
    return df


# ══════════════════════════════════════════════════════════════════════════════
# 多维度加权评分 / Multi-Dimensional Weighted Scoring
# ══════════════════════════════════════════════════════════════════════════════

def compute_weighted_score(df: pd.DataFrame, weights: Dict[str, float]) -> pd.DataFrame:
    """
    计算每个候选肽段的加权综合评分（归一化后加权求和）。
    Compute weighted composite score for each candidate (normalized then weighted).

    归一化策略 / Normalization strategy:
      - delta_g:        min-max 归一化后取反（越负越好 → 越高越好）
      - ipsae:          min-max 归一化后取反（越低越好 → 越高越好）
      - immunogenicity: 直接取反（越低越好 → 越高越好）
      - solubility:     直接使用（越高越好）

    Args:
        df:      候选肽段 DataFrame
        weights: 各维度权重字典

    Returns:
        附加 weighted_score 列的 DataFrame
    """
    df = df.copy()

    # ── delta_g 归一化（越负越好）────────────────────────────────────────────
    dg_min, dg_max = df["delta_g"].min(), df["delta_g"].max()
    dg_range = dg_max - dg_min if dg_max != dg_min else 1.0
    df["_norm_delta_g"] = 1.0 - (df["delta_g"] - dg_min) / dg_range  # 取反

    # ── ipsae 归一化（越低越好）──────────────────────────────────────────────
    ip_min, ip_max = df["ipsae"].min(), df["ipsae"].max()
    ip_range = ip_max - ip_min if ip_max != ip_min else 1.0
    df["_norm_ipsae"] = 1.0 - (df["ipsae"] - ip_min) / ip_range       # 取反

    # ── immunogenicity 归一化（越低越好）─────────────────────────────────────
    df["_norm_immunogenicity"] = 1.0 - df["immunogenicity"]            # 直接取反

    # ── solubility 归一化（越高越好）─────────────────────────────────────────
    df["_norm_solubility"] = df["solubility"]                          # 直接使用

    # ── 加权求和 ──────────────────────────────────────────────────────────────
    df["weighted_score"] = (
        weights["delta_g"]        * df["_norm_delta_g"]
        + weights["ipsae"]        * df["_norm_ipsae"]
        + weights["immunogenicity"] * df["_norm_immunogenicity"]
        + weights["solubility"]   * df["_norm_solubility"]
    ).round(4)

    # 清理中间列
    df.drop(columns=[c for c in df.columns if c.startswith("_norm_")], inplace=True)

    return df.sort_values("weighted_score", ascending=False).reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════════════
# 硬阈值过滤 / Hard-Threshold Filtering
# ══════════════════════════════════════════════════════════════════════════════

def apply_hard_filters(df: pd.DataFrame, filters: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    应用硬阈值过滤，返回通过和未通过的候选集合。
    Apply hard-threshold filters; return (passed, rejected) DataFrames.

    Args:
        df:      候选肽段 DataFrame
        filters: 硬阈值字典

    Returns:
        (passed_df, rejected_df) 元组
    """
    mask = (
        (df["delta_g"]        <= filters["delta_g_max"])
        & (df["ipsae"]        <= filters["ipsae_max"])
        & (df["immunogenicity"] <= filters["immunogenicity_max"])
        & (df["solubility"]   >= filters["solubility_min"])
        & (df["length"]       >= filters["length_min"])
        & (df["length"]       <= filters["length_max"])
    )
    passed   = df[mask].copy()
    rejected = df[~mask].copy()

    logger.info("硬阈值过滤: 通过 %d / 总计 %d (淘汰率 %.1f%%)",
                len(passed), len(df), 100 * len(rejected) / max(len(df), 1))
    return passed, rejected


# ══════════════════════════════════════════════════════════════════════════════
# 帕累托前沿分析 / Pareto Frontier Analysis
# ══════════════════════════════════════════════════════════════════════════════

def compute_pareto_frontier(
    df: pd.DataFrame,
    obj1_col: str = "delta_g",
    obj2_col: str = "immunogenicity",
) -> pd.DataFrame:
    """
    计算双目标帕累托前沿（最小化 delta_g，最小化 immunogenicity）。
    Compute bi-objective Pareto frontier (minimize delta_g, minimize immunogenicity).

    帕累托支配定义 / Pareto dominance definition:
      候选 A 支配候选 B，当且仅当 A 在所有目标上不劣于 B，且至少一个目标严格优于 B。
      Candidate A dominates B iff A is no worse on all objectives and strictly better on at least one.

    Args:
        df:       候选肽段 DataFrame
        obj1_col: 第一目标列名（最小化）
        obj2_col: 第二目标列名（最小化）

    Returns:
        帕累托前沿候选 DataFrame（is_pareto=True 的子集）
    """
    df = df.copy()
    n = len(df)
    is_pareto = np.ones(n, dtype=bool)  # 初始假设所有点都在前沿

    obj1 = df[obj1_col].values
    obj2 = df[obj2_col].values

    for i in range(n):
        if not is_pareto[i]:
            continue
        # 检查是否有其他点支配点 i
        # Check if any other point dominates point i
        dominated_by = (
            (obj1 <= obj1[i]) & (obj2 <= obj2[i])   # 所有目标不劣于 i
            & ((obj1 < obj1[i]) | (obj2 < obj2[i]))  # 至少一个目标严格优于 i
        )
        dominated_by[i] = False  # 排除自身
        if dominated_by.any():
            is_pareto[i] = False

    df["is_pareto"] = is_pareto

    n_pareto = is_pareto.sum()
    logger.info("帕累托前沿: %d 个候选（共 %d 个）", n_pareto, n)
    return df


def identify_elite_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    在帕累托前沿中识别 Elite Hits（同时满足 ΔG 和免疫原性双重严格阈值）。
    Identify Elite Hits on the Pareto frontier meeting strict dual thresholds.

    Args:
        df: 含 is_pareto 列的 DataFrame

    Returns:
        附加 is_elite 列的 DataFrame
    """
    df = df.copy()
    df["is_elite"] = (
        df["is_pareto"]
        & (df["delta_g"] < ELITE_DELTA_G_THRESHOLD)
        & (df["immunogenicity"] < ELITE_IMMUNOGENICITY_LIMIT)
    )
    n_elite = df["is_elite"].sum()
    logger.info("Elite Hits: %d 个", n_elite)
    return df


# ══════════════════════════════════════════════════════════════════════════════
# 可视化 / Visualization
# ══════════════════════════════════════════════════════════════════════════════

def plot_pareto_analysis(df: pd.DataFrame, output_path: str) -> None:
    """
    绘制帕累托前沿散点图，标注 Elite Hits，保存为 PNG。
    Plot Pareto frontier scatter chart with Elite Hit annotations; save as PNG.

    图层说明 / Layer description:
      - 灰色点:   普通候选（未通过帕累托）
      - 蓝色点:   帕累托前沿候选
      - 红色星形: Elite Hits（帕累托前沿 + 双重严格阈值）
      - 橙色虚线: 帕累托前沿连线（按 delta_g 排序）
      - 绿色虚线: Elite Hit 阈值参考线

    Args:
        df:          含 is_pareto / is_elite 列的 DataFrame
        output_path: 输出图片路径
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(
        "LigandAI · Pareto Frontier Analysis\n"
        "Multi-Objective Peptide Optimization: ΔG vs Immunogenicity",
        fontsize=14, fontweight="bold", y=1.01,
    )

    # ── 左图：帕累托前沿散点图 ────────────────────────────────────────────────
    ax = axes[0]

    # 普通候选（灰色）
    normal = df[~df["is_pareto"]]
    ax.scatter(
        normal["delta_g"], normal["immunogenicity"],
        c="lightgray", s=40, alpha=0.5, zorder=1, label="Non-Pareto Candidates",
    )

    # 帕累托前沿（蓝色）
    pareto = df[df["is_pareto"] & ~df["is_elite"]].sort_values("delta_g")
    ax.scatter(
        pareto["delta_g"], pareto["immunogenicity"],
        c="steelblue", s=70, alpha=0.8, zorder=2, label="Pareto Frontier",
    )

    # 帕累托前沿连线（橙色虚线）
    if len(pareto) > 1:
        pf_sorted = df[df["is_pareto"]].sort_values("delta_g")
        ax.plot(
            pf_sorted["delta_g"], pf_sorted["immunogenicity"],
            color="darkorange", linestyle="--", linewidth=1.5,
            alpha=0.7, zorder=3, label="Pareto Frontier Line",
        )

    # Elite Hits（红色星形）
    elite = df[df["is_elite"]]
    if not elite.empty:
        ax.scatter(
            elite["delta_g"], elite["immunogenicity"],
            c="crimson", s=200, marker="*", zorder=4,
            edgecolors="darkred", linewidths=0.8,
            label=f"Elite Hits (n={len(elite)})",
        )
        # 标注 Elite Hit 序列（前5个）
        for _, row in elite.head(5).iterrows():
            ax.annotate(
                row["sequence"][:10] + ("…" if len(row["sequence"]) > 10 else ""),
                xy=(row["delta_g"], row["immunogenicity"]),
                xytext=(8, 8), textcoords="offset points",
                fontsize=7, color="darkred",
                arrowprops=dict(arrowstyle="-", color="darkred", lw=0.8),
            )

    # 阈值参考线
    ax.axvline(x=ELITE_DELTA_G_THRESHOLD, color="green", linestyle=":", linewidth=1.2,
               alpha=0.7, label=f"ΔG threshold ({ELITE_DELTA_G_THRESHOLD})")
    ax.axhline(y=ELITE_IMMUNOGENICITY_LIMIT, color="purple", linestyle=":", linewidth=1.2,
               alpha=0.7, label=f"Immuno. threshold ({ELITE_IMMUNOGENICITY_LIMIT})")

    # 标注 Elite 区域
    ax.fill_betweenx(
        [0, ELITE_IMMUNOGENICITY_LIMIT],
        df["delta_g"].min() - 1, ELITE_DELTA_G_THRESHOLD,
        alpha=0.06, color="green", label="Elite Zone",
    )

    ax.set_xlabel("Binding Free Energy ΔG (kcal/mol)\n← More Favorable", fontsize=11)
    ax.set_ylabel("Immunogenicity Score\n↓ Lower Risk", fontsize=11)
    ax.set_title("Pareto Frontier: ΔG vs Immunogenicity", fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.invert_xaxis()  # ΔG 越负越好，x 轴左侧为优

    # ── 右图：加权评分分布直方图 ──────────────────────────────────────────────
    ax2 = axes[1]

    # 全部候选评分分布
    ax2.hist(
        df["weighted_score"], bins=25, color="steelblue",
        alpha=0.6, edgecolor="white", label="All Candidates",
    )
    # 帕累托前沿评分分布
    ax2.hist(
        df[df["is_pareto"]]["weighted_score"], bins=15,
        color="darkorange", alpha=0.7, edgecolor="white", label="Pareto Frontier",
    )
    # Elite Hits 评分分布
    if not elite.empty:
        ax2.hist(
            elite["weighted_score"], bins=8,
            color="crimson", alpha=0.9, edgecolor="white", label="Elite Hits",
        )

    # 均值参考线
    ax2.axvline(df["weighted_score"].mean(), color="navy", linestyle="--",
                linewidth=1.5, label=f"Mean ({df['weighted_score'].mean():.3f})")
    ax2.axvline(df[df["is_pareto"]]["weighted_score"].mean(),
                color="darkorange", linestyle="-.", linewidth=1.5,
                label=f"Pareto Mean ({df[df['is_pareto']]['weighted_score'].mean():.3f})")

    ax2.set_xlabel("Weighted Composite Score ↑", fontsize=11)
    ax2.set_ylabel("Count", fontsize=11)
    ax2.set_title("Weighted Score Distribution", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9, framealpha=0.9)
    ax2.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    logger.info("可视化图表已保存: %s", output_path)
    console.print(f"  [green]✓[/green] 可视化图表已保存: [bold]{output_path}[/bold]")


# ══════════════════════════════════════════════════════════════════════════════
# 筛选报告 / Filtering Report
# ══════════════════════════════════════════════════════════════════════════════

def display_filter_report(df_raw: pd.DataFrame, df_passed: pd.DataFrame, df_final: pd.DataFrame) -> None:
    """
    在终端展示多级筛选报告摘要。
    Display multi-stage filtering report summary in the terminal.

    Args:
        df_raw:    原始候选 DataFrame
        df_passed: 通过硬阈值过滤的 DataFrame
        df_final:  最终含帕累托/Elite 标注的 DataFrame
    """
    console.rule("[bold white]筛选报告 / Filtering Report[/bold white]")

    n_raw     = len(df_raw)
    n_passed  = len(df_passed)
    n_pareto  = int(df_final["is_pareto"].sum())
    n_elite   = int(df_final["is_elite"].sum())

    # 漏斗统计
    console.print(f"\n  [bold]筛选漏斗 / Filtering Funnel:[/bold]")
    console.print(f"  原始候选:       [white]{n_raw:>5}[/white]  (100.0%)")
    console.print(f"  硬阈值通过:     [cyan]{n_passed:>5}[/cyan]  ({100*n_passed/max(n_raw,1):.1f}%)")
    console.print(f"  帕累托前沿:     [blue]{n_pareto:>5}[/blue]  ({100*n_pareto/max(n_passed,1):.1f}% of passed)")
    console.print(f"  Elite Hits:     [bold red]{n_elite:>5}[/bold red]  ({100*n_elite/max(n_pareto,1):.1f}% of Pareto)\n")

    # Top-10 加权评分候选
    top10 = df_final.sort_values("weighted_score", ascending=False).head(10)

    table = Table(
        title="Top-10 Candidates by Weighted Score",
        header_style="bold magenta",
        border_style="bright_blue",
        show_lines=True,
    )
    table.add_column("Rank",    style="bold yellow", justify="center", width=6)
    table.add_column("Sequence",style="cyan",        justify="left",   width=22)
    table.add_column("ΔG",      style="green",       justify="right",  width=10)
    table.add_column("ipSAE",   style="blue",        justify="right",  width=8)
    table.add_column("Immuno.", style="red",         justify="right",  width=9)
    table.add_column("Solub.",  style="magenta",     justify="right",  width=8)
    table.add_column("W.Score", style="bold white",  justify="right",  width=9)
    table.add_column("Pareto",  style="yellow",      justify="center", width=8)
    table.add_column("Elite",   style="bold red",    justify="center", width=7)

    for rank, (_, row) in enumerate(top10.iterrows(), 1):
        table.add_row(
            str(rank),
            row["sequence"],
            f"{row['delta_g']:.3f}",
            f"{row['ipsae']:.3f}",
            f"{row['immunogenicity']:.3f}",
            f"{row['solubility']:.3f}",
            f"[bold]{row['weighted_score']:.4f}[/bold]",
            "✓" if row["is_pareto"] else "—",
            "[bold red]★[/bold red]" if row["is_elite"] else "—",
        )

    console.print(table)


def save_markdown_report(
    df_raw: pd.DataFrame,
    df_passed: pd.DataFrame,
    df_final: pd.DataFrame,
    output_path: str,
) -> None:
    """
    将筛选报告保存为 Markdown 文件，便于团队共享。
    Save the filtering report as a Markdown file for team sharing.

    Args:
        df_raw:      原始候选 DataFrame
        df_passed:   通过硬阈值的 DataFrame
        df_final:    最终标注 DataFrame
        output_path: 输出 Markdown 文件路径
    """
    from datetime import datetime

    n_raw    = len(df_raw)
    n_passed = len(df_passed)
    n_pareto = int(df_final["is_pareto"].sum())
    n_elite  = int(df_final["is_elite"].sum())

    top5_elite = df_final[df_final["is_elite"]].sort_values("weighted_score", ascending=False).head(5)

    lines = [
        "# LigandAI 多维度评分筛选报告",
        f"> 生成时间 / Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## 筛选漏斗 / Filtering Funnel",
        "",
        f"| 阶段 Stage              | 数量 Count | 比例 Ratio |",
        f"|------------------------|-----------|-----------|",
        f"| 原始候选 Raw Candidates | {n_raw}    | 100.0%    |",
        f"| 硬阈值通过 Hard Filter  | {n_passed} | {100*n_passed/max(n_raw,1):.1f}%  |",
        f"| 帕累托前沿 Pareto       | {n_pareto} | {100*n_pareto/max(n_passed,1):.1f}%  |",
        f"| Elite Hits             | {n_elite}  | {100*n_elite/max(n_pareto,1):.1f}%  |",
        "",
        "## 评分权重 / Score Weights",
        "",
        "| 维度 Dimension   | 权重 Weight |",
        "|-----------------|------------|",
    ]
    for dim, w in SCORE_WEIGHTS.items():
        lines.append(f"| {dim:<16} | {w:.2f}       |")

    lines += [
        "",
        "## Elite Hits Top-5",
        "",
        "| Rank | Sequence | ΔG (kcal/mol) | ipSAE | Immunogenicity | Solubility | W.Score |",
        "|------|----------|--------------|-------|----------------|------------|---------|",
    ]
    for rank, (_, row) in enumerate(top5_elite.iterrows(), 1):
        lines.append(
            f"| {rank} | `{row['sequence']}` | {row['delta_g']:.3f} | "
            f"{row['ipsae']:.3f} | {row['immunogenicity']:.3f} | "
            f"{row['solubility']:.3f} | **{row['weighted_score']:.4f}** |"
        )

    lines += [
        "",
        "## 可视化 / Visualization",
        "",
        f"![Pareto Analysis]({OUTPUT_PLOT})",
        "",
        "---",
        "*Generated by LigandAI Scoring Filter Example · ligandai==0.5.3*",
    ]

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    logger.info("Markdown 报告已保存: %s", output_path)
    console.print(f"  [green]✓[/green] 筛选报告已保存: [bold]{output_path}[/bold]")


# ══════════════════════════════════════════════════════════════════════════════
# 主流程 / Main
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """
    执行多维度评分筛选与帕累托前沿分析完整流程。
    Execute the full multi-dimensional scoring, filtering, and Pareto analysis pipeline.
    """
    console.print(Panel.fit(
        "[bold white]LigandAI Multi-Dimensional Scoring & Pareto Frontier Analysis[/bold white]\n"
        "[dim]Weights: ΔG×0.40 | ipSAE×0.25 | Immunogenicity×0.20 | Solubility×0.15[/dim]",
        border_style="bright_yellow",
        title="[bold yellow]📊 Scoring & Filter Pipeline[/bold yellow]",
    ))

    # ── Step 1: 加载数据 ──────────────────────────────────────────────────────
    console.rule("[bold cyan]Step 1 · 数据加载 / Data Loading[/bold cyan]")
    df_raw = load_or_simulate_data(INPUT_CSV)
    console.print(f"  [green]✓[/green] 原始候选数量: {len(df_raw)}")

    # ── Step 2: 加权评分 ──────────────────────────────────────────────────────
    console.rule("[bold cyan]Step 2 · 加权评分 / Weighted Scoring[/bold cyan]")
    df_scored = compute_weighted_score(df_raw, SCORE_WEIGHTS)
    console.print(f"  [green]✓[/green] 评分完成  均值={df_scored['weighted_score'].mean():.4f}  "
                  f"最高={df_scored['weighted_score'].max():.4f}")

    # ── Step 3: 硬阈值过滤 ───────────────────────────────────────────────────
    console.rule("[bold cyan]Step 3 · 硬阈值过滤 / Hard-Threshold Filtering[/bold cyan]")
    df_passed, df_rejected = apply_hard_filters(df_scored, HARD_FILTERS)
    console.print(f"  [green]✓[/green] 通过: {len(df_passed)}  淘汰: {len(df_rejected)}")

    # ── Step 4: 帕累托前沿分析 ───────────────────────────────────────────────
    console.rule("[bold cyan]Step 4 · 帕累托前沿 / Pareto Frontier Analysis[/bold cyan]")
    df_pareto = compute_pareto_frontier(df_passed, "delta_g", "immunogenicity")
    df_final  = identify_elite_hits(df_pareto)
    console.print(f"  [green]✓[/green] 帕累托前沿: {df_final['is_pareto'].sum()} 个  "
                  f"Elite Hits: [bold red]{df_final['is_elite'].sum()}[/bold red] 个")

    # ── Step 5: 可视化 ────────────────────────────────────────────────────────
    console.rule("[bold cyan]Step 5 · 可视化 / Visualization[/bold cyan]")
    plot_pareto_analysis(df_final, OUTPUT_PLOT)

    # ── Step 6: 报告输出 ──────────────────────────────────────────────────────
    console.rule("[bold cyan]Step 6 · 报告输出 / Report Export[/bold cyan]")
    display_filter_report(df_raw, df_passed, df_final)
    save_markdown_report(df_raw, df_passed, df_final, OUTPUT_REPORT)

    # 保存最终筛选结果 CSV
    out_csv = "ligandai_filtered_candidates.csv"
    df_final.to_csv(out_csv, index=False, encoding="utf-8-sig")
    console.print(f"  [green]✓[/green] 筛选结果已保存: [bold]{out_csv}[/bold]")

    console.print(Panel.fit(
        f"[bold green]✓ 筛选分析完成！[/bold green]\n"
        f"  原始候选:   {len(df_raw)}\n"
        f"  通过筛选:   {len(df_passed)}\n"
        f"  帕累托前沿: {int(df_final['is_pareto'].sum())}\n"
        f"  Elite Hits: [bold red]{int(df_final['is_elite'].sum())}[/bold red]\n"
        f"  图表:       {OUTPUT_PLOT}\n"
        f"  报告:       {OUTPUT_REPORT}",
        border_style="green",
        title="[bold green]Analysis Complete[/bold green]",
    ))


if __name__ == "__main__":
    main()
