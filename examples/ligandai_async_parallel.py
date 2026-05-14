#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ligandai_async_parallel.py
==========================
LigandAI 异步并行多靶点肽段设计示例
Async Parallel Multi-Target Peptide Design with LigandAI

功能概述 / Feature Overview:
  - 使用 asyncio 并发对 3 个靶点同时发起肽段生成请求
  - Semaphore 控制最大并发数，防止 API 限流
  - 异步 SSE 流式接收每个靶点的折叠进度
  - 聚合所有靶点结果，进行跨靶点比较分析
  - 输出每个靶点的 Top-5 候选及综合排名

依赖 / Requirements:
  pip install ligandai==0.5.3 pandas rich aiohttp

环境变量 / Environment Variable:
  export LIGANDAI_API_KEY="your_api_key_here"

作者 / Author: LigandAI Example Suite
版本 / Version: 1.0.0
"""

import os
import sys
import asyncio
import logging
import time
from typing import List, Dict, Any, Tuple
from dataclasses import dataclass, field

import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.columns import Columns
from rich import print as rprint

# ── LigandAI 异步客户端导入 ────────────────────────────────────────────────────
# Import LigandAI async client and exceptions
from ligandai import AsyncLigandAI
from ligandai.exceptions import (
    LigandAIAuthError,
    LigandAICreditError,
    LigandAITierError,
    LigandAIError,
)

# ── 日志配置 ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("ligandai.async_parallel")

console = Console()

# ── 并发控制参数 ───────────────────────────────────────────────────────────────
# Concurrency control parameters
MAX_CONCURRENT_TARGETS = 2      # 最大并发靶点数（Semaphore 上限）
NUM_CANDIDATES_PER_TARGET = 30  # 每个靶点生成的候选肽段数
TOP_N_PER_TARGET = 5            # 每个靶点展示的 Top-N 候选数

# ── 多靶点配置 ─────────────────────────────────────────────────────────────────
# Multi-target configuration — three oncology targets
TARGETS = [
    {
        "name":       "EGFR",
        "uniprot_id": "P00533",
        "description": "表皮生长因子受体 / Epidermal Growth Factor Receptor",
        "disease":    "非小细胞肺癌 / NSCLC",
        "color":      "cyan",
    },
    {
        "name":       "HER2",
        "uniprot_id": "P04626",
        "description": "人表皮生长因子受体2 / Human Epidermal Growth Factor Receptor 2",
        "disease":    "乳腺癌 / Breast Cancer",
        "color":      "magenta",
    },
    {
        "name":       "VEGFR2",
        "uniprot_id": "P35968",
        "description": "血管内皮生长因子受体2 / Vascular Endothelial Growth Factor Receptor 2",
        "disease":    "血管生成 / Angiogenesis",
        "color":      "yellow",
    },
]


# ══════════════════════════════════════════════════════════════════════════════
# 数据类 / Data Classes
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class TargetResult:
    """
    单个靶点的完整设计结果容器。
    Container for the complete design result of a single target.
    """
    target_name:   str
    uniprot_id:    str
    disease:       str
    candidates:    List[Dict[str, Any]] = field(default_factory=list)
    elapsed_sec:   float = 0.0
    error:         Optional[str] = None
    status:        str = "pending"   # pending | running | done | failed

    @property
    def success(self) -> bool:
        return self.status == "done" and self.error is None

    @property
    def top_candidate(self) -> Optional[Dict]:
        if not self.candidates:
            return None
        return min(self.candidates, key=lambda c: c.get("delta_g", 0))


# ══════════════════════════════════════════════════════════════════════════════
# 异步核心函数 / Async Core Functions
# ══════════════════════════════════════════════════════════════════════════════

async def design_for_target(
    client: AsyncLigandAI,
    target_cfg: Dict[str, Any],
    semaphore: asyncio.Semaphore,
) -> TargetResult:
    """
    对单个靶点执行完整的异步肽段设计流程（受 Semaphore 并发控制）。
    Execute the full async peptide design pipeline for a single target,
    gated by a Semaphore to limit concurrent API calls.

    Args:
        client:     LigandAI 异步客户端
        target_cfg: 靶点配置字典（name, uniprot_id, ...）
        semaphore:  asyncio.Semaphore，控制最大并发数

    Returns:
        TargetResult 对象，包含候选肽段及状态信息
    """
    name       = target_cfg["name"]
    uniprot_id = target_cfg["uniprot_id"]
    result     = TargetResult(
        target_name=name,
        uniprot_id=uniprot_id,
        disease=target_cfg["disease"],
    )
    t0 = time.monotonic()

    # Semaphore 限制并发：同一时刻最多 MAX_CONCURRENT_TARGETS 个靶点在运行
    # Semaphore gate: at most MAX_CONCURRENT_TARGETS targets run simultaneously
    async with semaphore:
        result.status = "running"
        logger.info("[%s] 开始异步设计流程 (UniProt: %s)", name, uniprot_id)

        try:
            # ── 2a. 异步获取靶点信息 ──────────────────────────────────────────
            target = await client.targets.get(uniprot_id=uniprot_id)
            logger.info("[%s] 靶点信息获取成功: %s", name, target.protein_name)

            # ── 2b. 异步解析结构 ──────────────────────────────────────────────
            structure = await client.structures.resolve(
                target_id=target.id,
                source="alphafold",
                confidence_threshold=65.0,
            )
            logger.info("[%s] 结构解析完成 (pLDDT=%.1f)", name, structure.mean_plddt)

            # ── 2c. 异步口袋检测 ──────────────────────────────────────────────
            pockets = await client.pockets.detect(
                structure_id=structure.id,
                min_druggability=0.45,
            )
            best_pocket = max(pockets, key=lambda p: p.druggability_score)
            logger.info("[%s] 最优口袋: druggability=%.3f", name, best_pocket.druggability_score)

            # ── 2d. 异步提交肽段生成任务 ──────────────────────────────────────
            job = await client.peptides.generate(
                target_id=target.id,
                pocket_id=best_pocket.id,
                num_candidates=NUM_CANDIDATES_PER_TARGET,
                length_range=(8, 18),
                model="pepgen-v2",
            )
            logger.info("[%s] 生成任务已提交: job_id=%s", name, job.id)

            # ── 2e. 异步 SSE 流式接收折叠进度 ────────────────────────────────
            candidates = []
            async for event in job.stream():
                if event.type == "candidate_ready":
                    candidates.append(event.data)
                elif event.type == "completed":
                    logger.info("[%s] 折叠完成，候选数=%d", name, len(candidates))
                    break
                elif event.type == "error":
                    raise LigandAIError(f"[{name}] Stream error: {event.data.get('message')}")

            # ── 2f. 异步 DeltaForge 评分 ──────────────────────────────────────
            scored = await client.deltaforge.score_batch(
                candidates=candidates,
                target_id=target.id,
                metrics=["delta_g", "ipsae", "immunogenicity", "solubility"],
            )

            # 按 delta_g 升序排列（越负越好）
            scored_sorted = sorted(scored, key=lambda c: c.get("delta_g", 0))
            result.candidates = scored_sorted
            result.status     = "done"
            logger.info("[%s] ✓ 完成，Top-1 ΔG=%.3f kcal/mol", name, scored_sorted[0].get("delta_g", 0))

        except LigandAIAuthError:
            result.status = "failed"
            result.error  = "AuthError: API Key 无效"
            logger.error("[%s] 认证失败", name)
            raise  # 认证错误应立即终止整个程序

        except LigandAICreditError:
            result.status = "failed"
            result.error  = "CreditError: 额度不足"
            logger.warning("[%s] 额度不足，跳过该靶点", name)

        except LigandAITierError:
            result.status = "failed"
            result.error  = "TierError: 套餐不支持"
            logger.warning("[%s] 套餐限制，跳过该靶点", name)

        except LigandAIError as e:
            result.status = "failed"
            result.error  = str(e)
            logger.error("[%s] LigandAI 错误: %s", name, e)

        except Exception as e:
            result.status = "failed"
            result.error  = f"未知错误: {e}"
            logger.exception("[%s] 未知异常", name)

        finally:
            result.elapsed_sec = time.monotonic() - t0

    return result


async def run_parallel_design(api_key: str) -> List[TargetResult]:
    """
    并行对所有靶点执行肽段设计，返回结果列表。
    Run peptide design for all targets in parallel; return list of results.

    Args:
        api_key: LigandAI API 密钥

    Returns:
        每个靶点的 TargetResult 列表（顺序与 TARGETS 一致）
    """
    # 初始化异步客户端（使用 async context manager 确保连接池正确关闭）
    # Initialize async client using context manager for proper connection pool cleanup
    async with AsyncLigandAI(api_key=api_key, timeout=180, max_retries=2) as client:

        # 创建 Semaphore，限制同时运行的靶点数
        # Create Semaphore to cap concurrent target processing
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_TARGETS)

        # 为每个靶点创建异步任务
        # Create async tasks for each target
        tasks = [
            asyncio.create_task(
                design_for_target(client, target_cfg, semaphore),
                name=f"design_{target_cfg['name']}",
            )
            for target_cfg in TARGETS
        ]

        console.print(
            f"  [cyan]→[/cyan] 并发任务数: {len(tasks)}  "
            f"最大并发限制: {MAX_CONCURRENT_TARGETS} (Semaphore)\n"
        )

        # 等待所有任务完成（gather 保留原始顺序）
        # Wait for all tasks; gather preserves original order
        results = await asyncio.gather(*tasks, return_exceptions=False)

    return list(results)


# ══════════════════════════════════════════════════════════════════════════════
# 结果聚合与展示 / Result Aggregation & Display
# ══════════════════════════════════════════════════════════════════════════════

def aggregate_results(results: List[TargetResult]) -> pd.DataFrame:
    """
    聚合所有靶点的候选肽段，构建跨靶点比较 DataFrame。
    Aggregate candidates from all targets into a cross-target comparison DataFrame.

    Args:
        results: 所有靶点的 TargetResult 列表

    Returns:
        聚合后的 DataFrame，包含靶点标签
    """
    all_records = []
    for res in results:
        if not res.success:
            continue
        for cand in res.candidates[:TOP_N_PER_TARGET]:
            all_records.append({
                "target":         res.target_name,
                "uniprot_id":     res.uniprot_id,
                "disease":        res.disease,
                "sequence":       cand.get("sequence", ""),
                "length":         len(cand.get("sequence", "")),
                "delta_g":        round(cand.get("delta_g", 0.0), 3),
                "ipsae":          round(cand.get("ipsae", 0.0), 3),
                "immunogenicity": round(cand.get("immunogenicity", 0.0), 3),
                "solubility":     round(cand.get("solubility", 0.0), 3),
            })

    df = pd.DataFrame(all_records)
    if df.empty:
        return df

    # 跨靶点综合评分（用于全局排名）
    df["global_score"] = (
        -0.45 * df["delta_g"].clip(-30, 0) / 30
        - 0.25 * df["ipsae"].clip(0, 10) / 10
        - 0.20 * df["immunogenicity"]
        + 0.10 * df["solubility"]
    ).round(4)

    return df.sort_values("global_score", ascending=False).reset_index(drop=True)


def display_per_target_summary(results: List[TargetResult]) -> None:
    """
    逐靶点展示设计摘要（状态、耗时、Top-1 候选）。
    Display per-target design summary (status, elapsed time, Top-1 candidate).
    """
    console.rule("[bold green]各靶点设计摘要 / Per-Target Design Summary[/bold green]")

    panels = []
    for res, cfg in zip(results, TARGETS):
        color = cfg["color"]
        if res.success and res.top_candidate:
            top = res.top_candidate
            body = (
                f"[green]✓ 成功[/green]  耗时: {res.elapsed_sec:.1f}s\n"
                f"候选数: {len(res.candidates)}\n"
                f"Top-1 序列: [{color}]{top.get('sequence', 'N/A')}[/{color}]\n"
                f"Top-1 ΔG:   {top.get('delta_g', 0):.3f} kcal/mol\n"
                f"Top-1 ipSAE: {top.get('ipsae', 0):.3f}"
            )
        else:
            body = (
                f"[red]✗ 失败[/red]  耗时: {res.elapsed_sec:.1f}s\n"
                f"错误: {res.error or '未知'}"
            )

        panels.append(Panel(
            body,
            title=f"[bold {color}]{res.target_name}[/bold {color}]  [{color}]{res.uniprot_id}[/{color}]",
            subtitle=f"[dim]{res.disease}[/dim]",
            border_style=color,
            width=40,
        ))

    console.print(Columns(panels))


def display_global_ranking(df: pd.DataFrame) -> None:
    """
    展示跨靶点全局排名表格（Top-15）。
    Display cross-target global ranking table (Top-15).
    """
    if df.empty:
        console.print("[red]无可用结果 / No results available[/red]")
        return

    console.rule("[bold white]跨靶点全局排名 / Cross-Target Global Ranking[/bold white]")

    table = Table(
        title="LigandAI · Cross-Target Global Peptide Ranking",
        header_style="bold magenta",
        border_style="bright_white",
        show_lines=True,
    )
    table.add_column("Rank",    style="bold yellow", justify="center", width=6)
    table.add_column("Target",  style="bold cyan",   justify="center", width=8)
    table.add_column("Sequence",style="white",        justify="left",   width=22)
    table.add_column("ΔG",      style="green",        justify="right",  width=10)
    table.add_column("ipSAE",   style="blue",         justify="right",  width=8)
    table.add_column("Immuno.", style="red",          justify="right",  width=9)
    table.add_column("Solub.",  style="magenta",      justify="right",  width=8)
    table.add_column("Score ↑", style="bold white",   justify="right",  width=10)

    for i, row in df.head(15).iterrows():
        table.add_row(
            str(i + 1),
            row["target"],
            row["sequence"],
            f"{row['delta_g']:.3f}",
            f"{row['ipsae']:.3f}",
            f"{row['immunogenicity']:.3f}",
            f"{row['solubility']:.3f}",
            f"[bold]{row['global_score']:.4f}[/bold]",
        )

    console.print(table)

    # 保存聚合结果
    out_path = "ligandai_parallel_results.csv"
    df.to_csv(out_path, index=False, encoding="utf-8-sig")
    console.print(f"\n  [green]✓[/green] 聚合结果已保存: [bold]{out_path}[/bold]")


# ══════════════════════════════════════════════════════════════════════════════
# 主入口 / Main Entry Point
# ══════════════════════════════════════════════════════════════════════════════

async def async_main() -> None:
    """
    异步主函数：初始化、并行设计、聚合展示。
    Async main: initialize, parallel design, aggregate and display.
    """
    console.print(Panel.fit(
        "[bold white]LigandAI Async Parallel Multi-Target Design[/bold white]\n"
        f"[dim]Targets: {', '.join(t['name'] for t in TARGETS)}  |  "
        f"Max Concurrency: {MAX_CONCURRENT_TARGETS}  |  "
        f"Candidates/Target: {NUM_CANDIDATES_PER_TARGET}[/dim]",
        border_style="bright_magenta",
        title="[bold magenta]⚡ Async Parallel Pipeline[/bold magenta]",
    ))

    # 读取 API Key
    api_key = os.environ.get("LIGANDAI_API_KEY", "").strip()
    if not api_key:
        console.print("[bold red]✗ 未设置 LIGANDAI_API_KEY 环境变量[/bold red]")
        sys.exit(1)

    t_start = time.monotonic()

    try:
        console.rule("[bold cyan]启动并行设计任务 / Launching Parallel Design Tasks[/bold cyan]")
        results = await run_parallel_design(api_key)

    except LigandAIAuthError:
        console.print("[bold red]✗ API Key 认证失败，请检查 LIGANDAI_API_KEY[/bold red]")
        sys.exit(2)

    total_elapsed = time.monotonic() - t_start
    console.print(f"\n  [dim]并行总耗时 / Total parallel elapsed: {total_elapsed:.1f}s[/dim]")

    # 展示各靶点摘要
    display_per_target_summary(results)

    # 聚合并展示全局排名
    df_global = aggregate_results(results)
    display_global_ranking(df_global)

    # 统计成功/失败数
    n_success = sum(1 for r in results if r.success)
    n_failed  = len(results) - n_success
    console.print(Panel.fit(
        f"[bold green]✓ 并行设计完成[/bold green]\n"
        f"  成功靶点: {n_success}/{len(results)}\n"
        f"  失败靶点: {n_failed}/{len(results)}\n"
        f"  总候选数: {sum(len(r.candidates) for r in results if r.success)}\n"
        f"  总耗时:   {total_elapsed:.1f}s",
        border_style="green",
        title="[bold green]Parallel Design Complete[/bold green]",
    ))


def main() -> None:
    """同步入口，启动 asyncio 事件循环 / Sync entry point launching asyncio event loop."""
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
