#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ligandai_e2e_pipeline.py
========================
LigandAI 端到端蛋白设计完整流程示例
End-to-End Protein Design Pipeline with LigandAI

流程概览 / Pipeline Overview:
  1. 靶点发现   (Target Discovery)      — 通过 UniProt ID 或基因名检索靶点
  2. 结构解析   (Structure Resolution)  — 获取 AlphaFold / PDB 结构
  3. 口袋分析   (Pocket Analysis)       — 识别可成药结合口袋
  4. 肽段生成   (Peptide Generation)    — 基于口袋约束生成候选肽段
  5. 流式折叠   (Streaming Folding)     — 实时 SSE 进度跟踪结构折叠
  6. DeltaForge 评分                    — 多维度结合亲和力评分
  7. 结果输出   (Result Export)         — Top-10 候选肽段 DataFrame + CSV

依赖 / Requirements:
  pip install ligandai==0.5.3 pandas rich

环境变量 / Environment Variable:
  export LIGANDAI_API_KEY="your_api_key_here"

作者 / Author: LigandAI Example Suite
版本 / Version: 1.0.0
"""

import os
import sys
import time
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional

import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.panel import Panel
from rich import print as rprint

# ── LigandAI SDK 导入 ──────────────────────────────────────────────────────────
# Import LigandAI SDK and exception classes
import ligandai
from ligandai import LigandAI
from ligandai.exceptions import (
    LigandAIAuthError,      # API Key 无效或过期 / Invalid or expired API key
    LigandAICreditError,    # 账户额度不足 / Insufficient account credits
    LigandAITierError,      # 当前套餐不支持该功能 / Feature not available in current tier
    LigandAIError,          # 通用 SDK 错误基类 / Generic SDK error base class
)

# ── 日志配置 ───────────────────────────────────────────────────────────────────
# Configure structured logging for production-grade traceability
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("ligandai.pipeline")

# ── 全局常量 ───────────────────────────────────────────────────────────────────
# Global constants — adjust to your target of interest
TARGET_UNIPROT_ID   = "P00533"          # EGFR (Epidermal Growth Factor Receptor)
TARGET_GENE_NAME    = "EGFR"
PEPTIDE_LENGTH_MIN  = 8                 # 最短肽段长度 / Minimum peptide length (aa)
PEPTIDE_LENGTH_MAX  = 20                # 最长肽段长度 / Maximum peptide length (aa)
NUM_CANDIDATES      = 50               # 生成候选数量 / Number of candidates to generate
TOP_N               = 10               # 输出 Top-N 候选 / Top-N candidates to display
OUTPUT_CSV          = "ligandai_top_candidates.csv"

console = Console()


# ══════════════════════════════════════════════════════════════════════════════
# 辅助函数 / Helper Functions
# ══════════════════════════════════════════════════════════════════════════════

def get_api_key() -> str:
    """
    从环境变量读取 API Key，若未设置则抛出异常。
    Read API key from environment variable; raise if not set.
    """
    api_key = os.environ.get("LIGANDAI_API_KEY", "").strip()
    if not api_key:
        console.print(
            "[bold red]✗ 未找到 LIGANDAI_API_KEY 环境变量。\n"
            "  请执行: export LIGANDAI_API_KEY='your_key_here'[/bold red]"
        )
        sys.exit(1)
    logger.info("API Key 已加载 (前8位: %s...)", api_key[:8])
    return api_key


def build_client(api_key: str) -> LigandAI:
    """
    初始化并返回 LigandAI 同步客户端。
    Initialize and return a synchronous LigandAI client.

    Args:
        api_key: LigandAI API 密钥

    Returns:
        已配置的 LigandAI 客户端实例
    """
    client = LigandAI(
        api_key=api_key,
        timeout=120,            # 请求超时（秒）/ Request timeout in seconds
        max_retries=3,          # 自动重试次数 / Auto-retry count on transient errors
    )
    logger.info("LigandAI 客户端初始化完成 (SDK v%s)", ligandai.__version__)
    return client


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 — 靶点发现 / Target Discovery
# ══════════════════════════════════════════════════════════════════════════════

def discover_target(client: LigandAI, uniprot_id: str) -> Dict[str, Any]:
    """
    通过 UniProt ID 检索靶点元数据，包括基因名、物种、功能注释等。
    Retrieve target metadata via UniProt ID (gene name, organism, function, etc.).

    Args:
        client:     LigandAI 客户端
        uniprot_id: UniProt 蛋白质 ID（如 "P00533"）

    Returns:
        靶点元数据字典 / Target metadata dictionary
    """
    console.rule("[bold cyan]Step 1 · 靶点发现 / Target Discovery[/bold cyan]")
    logger.info("正在检索靶点: %s", uniprot_id)

    target = client.targets.get(uniprot_id=uniprot_id)

    console.print(f"  [green]✓[/green] 靶点名称:  [bold]{target.gene_name}[/bold]")
    console.print(f"  [green]✓[/green] 物种:      {target.organism}")
    console.print(f"  [green]✓[/green] 蛋白全名:  {target.protein_name}")
    console.print(f"  [green]✓[/green] 功能分类:  {target.function_class}")
    console.print(f"  [green]✓[/green] 疾病关联:  {', '.join(target.disease_associations[:3])}")

    return target


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — 结构解析 / Structure Resolution
# ══════════════════════════════════════════════════════════════════════════════

def resolve_structure(client: LigandAI, target) -> Dict[str, Any]:
    """
    获取靶点的三维结构（优先 AlphaFold，回退至 PDB 实验结构）。
    Fetch 3D structure for the target (AlphaFold preferred, PDB as fallback).

    Args:
        client: LigandAI 客户端
        target: 靶点对象

    Returns:
        结构对象，包含 PDB 坐标和质量指标
    """
    console.rule("[bold cyan]Step 2 · 结构解析 / Structure Resolution[/bold cyan]")
    logger.info("正在解析结构: %s", target.uniprot_id)

    structure = client.structures.resolve(
        target_id=target.id,
        source="alphafold",         # 优先使用 AlphaFold2 预测结构
        confidence_threshold=70.0,  # pLDDT 置信度阈值（0-100）
    )

    console.print(f"  [green]✓[/green] 结构来源:  {structure.source}")
    console.print(f"  [green]✓[/green] 分辨率:    {structure.resolution_angstrom:.2f} Å")
    console.print(f"  [green]✓[/green] 残基数量:  {structure.residue_count}")
    console.print(f"  [green]✓[/green] pLDDT 均值: {structure.mean_plddt:.1f}")

    return structure


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — 口袋分析 / Binding Pocket Analysis
# ══════════════════════════════════════════════════════════════════════════════

def analyze_pockets(client: LigandAI, structure) -> List[Dict[str, Any]]:
    """
    对蛋白结构进行口袋检测，返回按可成药性排序的口袋列表。
    Detect binding pockets on the protein structure, ranked by druggability.

    Args:
        client:    LigandAI 客户端
        structure: 蛋白结构对象

    Returns:
        口袋列表（按 druggability_score 降序）
    """
    console.rule("[bold cyan]Step 3 · 口袋分析 / Pocket Analysis[/bold cyan]")
    logger.info("正在分析结合口袋...")

    pockets = client.pockets.detect(
        structure_id=structure.id,
        algorithm="fpocket_v4",     # 使用 FPocket v4 算法
        min_volume=200.0,           # 最小口袋体积（Å³）
        min_druggability=0.5,       # 最低可成药性评分
    )

    # 按可成药性评分降序排列
    pockets_sorted = sorted(pockets, key=lambda p: p.druggability_score, reverse=True)

    console.print(f"  [green]✓[/green] 检测到口袋数量: {len(pockets_sorted)}")
    for i, pocket in enumerate(pockets_sorted[:3], 1):
        console.print(
            f"  [cyan]Pocket #{i}[/cyan]  体积={pocket.volume_angstrom3:.0f} Å³  "
            f"可成药性={pocket.druggability_score:.3f}  "
            f"残基数={pocket.residue_count}"
        )

    # 选取最优口袋用于后续肽段生成
    best_pocket = pockets_sorted[0]
    logger.info("选定最优口袋 ID: %s (druggability=%.3f)", best_pocket.id, best_pocket.druggability_score)
    return best_pocket


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — 肽段生成 / Peptide Generation
# ══════════════════════════════════════════════════════════════════════════════

def generate_peptides(client: LigandAI, target, pocket) -> Any:
    """
    基于口袋约束，调用 LigandAI 生成候选肽段序列。
    Generate candidate peptide sequences constrained by the selected pocket.

    Args:
        client: LigandAI 客户端
        target: 靶点对象
        pocket: 最优结合口袋对象

    Returns:
        生成任务对象（Job），可用于流式跟踪进度
    """
    console.rule("[bold cyan]Step 4 · 肽段生成 / Peptide Generation[/bold cyan]")
    logger.info("提交肽段生成任务，候选数量=%d", NUM_CANDIDATES)

    job = client.peptides.generate(
        target_id=target.id,
        pocket_id=pocket.id,
        num_candidates=NUM_CANDIDATES,
        length_range=(PEPTIDE_LENGTH_MIN, PEPTIDE_LENGTH_MAX),
        model="pepgen-v2",          # 使用 PepGen v2 生成模型
        constraints={
            "no_cysteines": False,  # 允许半胱氨酸（可形成二硫键）
            "max_hydrophobicity": 0.7,
            "charge_range": (-2, 4),
        },
        diversity_factor=0.85,      # 序列多样性因子（0-1）
    )

    console.print(f"  [green]✓[/green] 任务已提交，Job ID: [bold]{job.id}[/bold]")
    console.print(f"  [green]✓[/green] 预计完成时间: ~{job.estimated_seconds}s")
    return job


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 — 流式折叠 / Streaming Folding with SSE Progress
# ══════════════════════════════════════════════════════════════════════════════

def stream_folding(job) -> List[Dict[str, Any]]:
    """
    通过 Server-Sent Events (SSE) 实时跟踪折叠进度，直至任务完成。
    Track folding progress in real-time via SSE until the job completes.

    Args:
        job: 肽段生成任务对象

    Returns:
        折叠完成的候选肽段列表
    """
    console.rule("[bold cyan]Step 5 · 流式折叠 / Streaming Folding (SSE)[/bold cyan]")
    logger.info("开始流式接收折叠进度...")

    candidates = []
    last_progress = -1

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task_id = progress.add_task("折叠中 / Folding...", total=100)

        # job.stream() 返回 SSE 事件生成器
        # job.stream() returns an SSE event generator
        for event in job.stream():
            event_type = event.type

            if event_type == "progress":
                # 进度更新事件 / Progress update event
                pct = int(event.data.get("percent", 0))
                stage = event.data.get("stage", "")
                if pct != last_progress:
                    progress.update(task_id, completed=pct, description=f"[cyan]{stage}[/cyan]")
                    last_progress = pct

            elif event_type == "candidate_ready":
                # 单个候选肽段折叠完成 / Single candidate folding complete
                candidate = event.data
                candidates.append(candidate)
                logger.debug("候选肽段就绪: %s (ipsae=%.3f)", candidate.get("sequence"), candidate.get("ipsae", 0))

            elif event_type == "completed":
                # 全部任务完成 / All tasks completed
                progress.update(task_id, completed=100, description="[green]折叠完成 / Folding Complete[/green]")
                logger.info("折叠任务完成，共 %d 个候选肽段", len(candidates))
                break

            elif event_type == "error":
                # 服务端错误 / Server-side error
                error_msg = event.data.get("message", "Unknown error")
                logger.error("流式折叠错误: %s", error_msg)
                raise LigandAIError(f"Folding stream error: {error_msg}")

    console.print(f"  [green]✓[/green] 折叠完成，获得 {len(candidates)} 个候选肽段")
    return candidates


# ══════════════════════════════════════════════════════════════════════════════
# Step 6 — DeltaForge 评分 / DeltaForge Scoring
# ══════════════════════════════════════════════════════════════════════════════

def deltaforge_scoring(client: LigandAI, candidates: List[Dict], target) -> List[Dict[str, Any]]:
    """
    对折叠后的候选肽段进行 DeltaForge 多维度评分。
    Score folded peptide candidates using DeltaForge multi-dimensional scoring.

    评分维度 / Scoring dimensions:
      - delta_g:        结合自由能（kcal/mol），越负越好
      - ipsae:          界面预测结构对齐误差，越低越好
      - immunogenicity: 免疫原性风险（0-1），越低越好
      - solubility:     溶解度预测（0-1），越高越好

    Args:
        client:     LigandAI 客户端
        candidates: 候选肽段列表
        target:     靶点对象

    Returns:
        附加评分信息的候选肽段列表
    """
    console.rule("[bold cyan]Step 6 · DeltaForge 评分 / DeltaForge Scoring[/bold cyan]")
    logger.info("提交 DeltaForge 评分，共 %d 个候选", len(candidates))

    # 批量提交评分请求
    scored_results = client.deltaforge.score_batch(
        candidates=candidates,
        target_id=target.id,
        metrics=["delta_g", "ipsae", "immunogenicity", "solubility", "stability"],
        force_field="amber_ff19sb",  # 使用 AMBER ff19SB 力场
    )

    console.print(f"  [green]✓[/green] DeltaForge 评分完成，共 {len(scored_results)} 条记录")
    return scored_results


# ══════════════════════════════════════════════════════════════════════════════
# Step 7 — 结果输出 / Result Export
# ══════════════════════════════════════════════════════════════════════════════

def build_dataframe(scored_results: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    将评分结果转换为 pandas DataFrame，并按综合评分排序。
    Convert scored results to a pandas DataFrame sorted by composite score.

    Args:
        scored_results: DeltaForge 评分结果列表

    Returns:
        排序后的 DataFrame
    """
    records = []
    for r in scored_results:
        records.append({
            "sequence":       r.get("sequence", ""),
            "length":         len(r.get("sequence", "")),
            "delta_g":        round(r.get("delta_g", 0.0), 3),
            "ipsae":          round(r.get("ipsae", 0.0), 3),
            "immunogenicity": round(r.get("immunogenicity", 0.0), 3),
            "solubility":     round(r.get("solubility", 0.0), 3),
            "stability":      round(r.get("stability", 0.0), 3),
            "pdb_url":        r.get("pdb_url", ""),
            "candidate_id":   r.get("id", ""),
        })

    df = pd.DataFrame(records)

    # 计算综合评分（加权归一化）
    # Compute composite score (weighted normalization)
    df["composite_score"] = (
        -0.40 * df["delta_g"].clip(-30, 0) / 30          # 结合能（负值越好）
        - 0.25 * df["ipsae"].clip(0, 10) / 10             # 结构误差（越低越好）
        - 0.20 * df["immunogenicity"]                      # 免疫原性（越低越好）
        + 0.15 * df["solubility"]                          # 溶解度（越高越好）
    ).round(4)

    # 按综合评分降序排列
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    df.index += 1  # 从 1 开始编号
    return df


def display_top_candidates(df: pd.DataFrame, top_n: int = TOP_N) -> None:
    """
    使用 Rich 在终端美观地展示 Top-N 候选肽段表格。
    Display Top-N candidates in a beautiful Rich table in the terminal.

    Args:
        df:    候选肽段 DataFrame
        top_n: 展示数量
    """
    console.rule(f"[bold green]Top-{top_n} 候选肽段 / Top-{top_n} Peptide Candidates[/bold green]")

    table = Table(
        title=f"LigandAI · Top-{top_n} Peptide Candidates for {TARGET_GENE_NAME}",
        show_header=True,
        header_style="bold magenta",
        border_style="bright_blue",
    )
    table.add_column("Rank",          style="bold yellow", justify="center", width=6)
    table.add_column("Sequence",      style="cyan",        justify="left",   width=24)
    table.add_column("Len",           style="white",       justify="center", width=5)
    table.add_column("ΔG (kcal/mol)", style="green",       justify="right",  width=14)
    table.add_column("ipSAE",         style="blue",        justify="right",  width=8)
    table.add_column("Immuno.",       style="red",         justify="right",  width=9)
    table.add_column("Solubility",    style="magenta",     justify="right",  width=10)
    table.add_column("Composite ↑",   style="bold white",  justify="right",  width=12)

    for rank, row in df.head(top_n).iterrows():
        table.add_row(
            str(rank),
            row["sequence"],
            str(row["length"]),
            f"{row['delta_g']:.3f}",
            f"{row['ipsae']:.3f}",
            f"{row['immunogenicity']:.3f}",
            f"{row['solubility']:.3f}",
            f"[bold]{row['composite_score']:.4f}[/bold]",
        )

    console.print(table)


def save_results(df: pd.DataFrame, output_path: str) -> None:
    """
    将完整结果保存为 CSV 文件。
    Save the full results to a CSV file.

    Args:
        df:          候选肽段 DataFrame
        output_path: 输出文件路径
    """
    df.to_csv(output_path, index=True, index_label="rank", encoding="utf-8-sig")
    console.print(f"\n  [green]✓[/green] 结果已保存至: [bold]{output_path}[/bold]  ({len(df)} 条记录)")
    logger.info("CSV 已写入: %s", output_path)


# ══════════════════════════════════════════════════════════════════════════════
# 主流程入口 / Main Pipeline Entry Point
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """
    执行完整的 LigandAI 端到端蛋白设计流程。
    Execute the complete LigandAI end-to-end protein design pipeline.
    """
    start_time = time.time()

    console.print(Panel.fit(
        "[bold white]LigandAI End-to-End Protein Design Pipeline[/bold white]\n"
        f"[dim]Target: {TARGET_GENE_NAME} ({TARGET_UNIPROT_ID})  |  "
        f"Candidates: {NUM_CANDIDATES}  |  SDK: ligandai==0.5.3[/dim]",
        border_style="bright_cyan",
        title="[bold cyan]🧬 LigandAI Pipeline[/bold cyan]",
    ))

    try:
        # ── 初始化客户端 ──────────────────────────────────────────────────────
        api_key = get_api_key()
        client  = build_client(api_key)

        # ── Step 1: 靶点发现 ──────────────────────────────────────────────────
        target = discover_target(client, TARGET_UNIPROT_ID)

        # ── Step 2: 结构解析 ──────────────────────────────────────────────────
        structure = resolve_structure(client, target)

        # ── Step 3: 口袋分析 ──────────────────────────────────────────────────
        best_pocket = analyze_pockets(client, structure)

        # ── Step 4: 肽段生成 ──────────────────────────────────────────────────
        job = generate_peptides(client, target, best_pocket)

        # ── Step 5: 流式折叠 ──────────────────────────────────────────────────
        candidates = stream_folding(job)

        # ── Step 6: DeltaForge 评分 ───────────────────────────────────────────
        scored_results = deltaforge_scoring(client, candidates, target)

        # ── Step 7: 结果输出 ──────────────────────────────────────────────────
        df = build_dataframe(scored_results)
        display_top_candidates(df, top_n=TOP_N)
        save_results(df, OUTPUT_CSV)

    # ── 错误处理 / Error Handling ─────────────────────────────────────────────
    except LigandAIAuthError as e:
        # API Key 无效或已过期
        console.print(f"\n[bold red]✗ 认证失败 (AuthError):[/bold red] {e}")
        console.print("  请检查 LIGANDAI_API_KEY 是否正确，或访问 https://app.ligandai.com/settings/api")
        sys.exit(2)

    except LigandAICreditError as e:
        # 账户额度不足，无法继续
        console.print(f"\n[bold red]✗ 额度不足 (CreditError):[/bold red] {e}")
        console.print("  请前往 https://app.ligandai.com/billing 充值或升级套餐")
        sys.exit(3)

    except LigandAITierError as e:
        # 当前套餐不支持该功能（如 DeltaForge 需要 Pro 套餐）
        console.print(f"\n[bold red]✗ 套餐限制 (TierError):[/bold red] {e}")
        console.print("  当前功能需要更高级套餐，请访问 https://app.ligandai.com/pricing")
        sys.exit(4)

    except LigandAIError as e:
        # 其他 SDK 错误
        console.print(f"\n[bold red]✗ LigandAI 错误:[/bold red] {e}")
        logger.exception("LigandAI SDK 错误详情")
        sys.exit(5)

    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ 用户中断，流程已终止。[/yellow]")
        sys.exit(0)

    finally:
        elapsed = time.time() - start_time
        console.print(f"\n[dim]总耗时 / Total elapsed: {elapsed:.1f}s[/dim]")

    # ── 流程完成摘要 ──────────────────────────────────────────────────────────
    console.print(Panel.fit(
        f"[bold green]✓ 流程完成！[/bold green]\n"
        f"  靶点:     {TARGET_GENE_NAME} ({TARGET_UNIPROT_ID})\n"
        f"  候选总数: {len(df)}\n"
        f"  Top-1 序列: [cyan]{df.iloc[0]['sequence']}[/cyan]\n"
        f"  Top-1 ΔG:   {df.iloc[0]['delta_g']:.3f} kcal/mol\n"
        f"  结果文件: {OUTPUT_CSV}",
        border_style="green",
        title="[bold green]Pipeline Complete[/bold green]",
    ))


if __name__ == "__main__":
    main()
