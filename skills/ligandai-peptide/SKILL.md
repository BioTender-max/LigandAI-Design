---
name: ligandai-peptide
description: |
  【中文】基于 LigandAI LigandForge v6.5 引擎的靶向肽段生成技能。支持口袋靶向（pocket_targeted）、
  支架模板（scaffold_templated）、去免疫原性（deimmunized）和药物偶联（drug_conjugate）四种策略，
  结合 Boltz-2 自动折叠与 SSE 流式进度监听，实现从靶点到候选肽段的端到端自动化设计。
  
  [EN] Peptide generation skill powered by LigandAI LigandForge v6.5 engine. Supports four targeting
  strategies: pocket_targeted, scaffold_templated, deimmunized, and drug_conjugate. Integrates Boltz-2
  auto-folding with SSE streaming progress, enabling end-to-end automated design from target to
  ranked peptide candidates.
license: Apache-2.0
category: design-tools
tags:
  - peptide
  - generation
  - ligandforge
  - binder
  - pocket-targeted
  - scaffold
  - deimmunization

triggers:
  # 中文触发场景
  - "为 EGFR 生成靶向肽段"
  - "用 LigandAI 设计肽段结合剂"
  - "生成去免疫原性肽段"
  - "基于支架模板设计肽段"
  - "口袋靶向肽段生成"
  - "为肝脏靶点设计递送肽"
  - "生成 300 条候选肽段并自动折叠"
  - "LigandForge 肽段设计"
  - "设计药物偶联肽段"
  - "流式监听肽段生成进度"
  # English triggers
  - "generate peptide binders for EGFR"
  - "design pocket-targeted peptides using LigandAI"
  - "create deimmunized peptide candidates"
  - "scaffold-templated peptide generation"
  - "generate 300 peptides with auto-folding"
  - "LigandForge peptide design"
  - "design drug-conjugate peptides"
  - "stream peptide generation progress"

inputs:
  - name: gene
    type: string
    required: true
    description: "目标基因符号（HGNC），如 'EGFR', 'ASGR1' / Target gene symbol (HGNC)"
  - name: num_peptides
    type: integer
    required: false
    default: 100
    description: "生成肽段数量（受 Tier 限制）/ Number of peptides to generate (Tier-limited)"
  - name: targeting_strategy
    type: string
    required: false
    default: "pocket_targeted"
    description: "靶向策略：pocket_targeted / scaffold_templated / deimmunized / drug_conjugate"
  - name: target_residues
    type: list[int]
    required: false
    description: "热点残基编号列表（口袋靶向时指定）/ Hotspot residue indices for pocket targeting"
  - name: auto_fold
    type: boolean
    required: false
    default: false
    description: "是否自动折叠 Top N 肽段 / Whether to auto-fold top N peptides"
  - name: top_n_fold
    type: integer
    required: false
    default: 10
    description: "自动折叠的肽段数量 / Number of top peptides to auto-fold"
  - name: sampling_steps
    type: integer
    required: false
    default: 50
    description: "扩散采样步数（50/100/150/200，步数越多质量越高）/ Diffusion sampling steps"
  - name: num_trajectories
    type: integer
    required: false
    default: 1
    description: "每条肽段的折叠轨迹数 / Number of folding trajectories per peptide"
  - name: deimmunize
    type: boolean
    required: false
    default: false
    description: "是否启用去免疫原性优化 / Whether to enable deimmunization optimization"
  - name: scaffold_template
    type: string
    required: false
    description: "支架模板序列（scaffold_templated 策略时必填）/ Scaffold template sequence"

outputs:
  - name: peptides
    type: list[PeptideResult]
    description: "生成的肽段列表，按综合评分排序 / Generated peptides sorted by composite score"
  - name: peptides[].sequence
    type: string
    description: "肽段氨基酸序列（单字母码）/ Peptide amino acid sequence (one-letter code)"
  - name: peptides[].ipsae
    type: float
    description: "界面预测对齐误差（<0.3 高置信）/ Interface predicted alignment error"
  - name: peptides[].delta_g
    type: float
    description: "结合自由能 kcal/mol（越负越强）/ Binding free energy in kcal/mol"
  - name: peptides[].kd
    type: float
    description: "预测解离常数 nM / Predicted dissociation constant in nM"
  - name: peptides[].immunogenicity
    type: float
    description: "免疫原性风险评分 0-1（<0.3 低风险）/ Immunogenicity risk score 0-1"
  - name: peptides[].solubility
    type: float
    description: "GRAVY 溶解度指数（负值=亲水）/ GRAVY solubility index (negative=hydrophilic)"
---

# LigandAI Peptide — 靶向肽段生成技能

## 1. 功能概述

### 中文说明

LigandAI LigandForge v6.5 是专为治疗性肽段和靶向递送肽段设计的生成式 AI 引擎，基于扩散模型（Diffusion Model）
在蛋白质序列-结构联合空间中进行采样。相比传统噬菌体展示或计算对接方法，LigandForge 可在数分钟内生成数百条
具有预测结合亲和力的候选肽段，并通过集成的 Boltz-2 折叠引擎进行结构验证。

### 能力总览表

| 策略 | 英文名 | 适用场景 | 优势 | 限制 |
|------|--------|---------|------|------|
| 口袋靶向 | `pocket_targeted` | 已知结合口袋的受体 | 高结合特异性，可指定热点残基 | 需要结构信息 |
| 支架模板 | `scaffold_templated` | 基于已知活性肽优化 | 保留核心药效团，提高成功率 | 需要模板序列 |
| 去免疫原性 | `deimmunized` | 治疗性肽段，降低免疫反应 | 减少 MHC-II 呈递风险 | 可能降低亲和力 |
| 药物偶联 | `drug_conjugate` | ADC/PDC 载体肽设计 | 保留偶联位点（Cys/Lys） | 需指定偶联化学 |

---

## 2. 快速开始

```bash
pip install ligandai>=2.1.0
```

```python
import ligandai

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# 最简单的肽段生成示例
# Simplest peptide generation example
job = client.peptides.generate(
    gene="EGFR",
    num_peptides=100,
    targeting_strategy="pocket_targeted",
)

# 等待完成并获取结果
results = job.wait()
print(f"Generated {len(results.peptides)} peptides")
print(f"Best peptide: {results.peptides[0].sequence}")
print(f"Best ΔG: {results.peptides[0].delta_g:.2f} kcal/mol")
```

---

## 3. `client.peptides.generate()` 完整参数表

```python
job = client.peptides.generate(
    # ── 必填参数 / Required ──────────────────────────────────────────────────
    gene="EGFR",                      # 目标基因（HGNC 符号）

    # ── 生成控制 / Generation control ───────────────────────────────────────
    num_peptides=300,                  # 生成数量（受 Tier 限制）
    peptide_length=(8, 20),            # 肽段长度范围（氨基酸数）
    targeting_strategy="pocket_targeted",  # 靶向策略

    # ── 口袋靶向参数 / Pocket targeting ─────────────────────────────────────
    target_residues=[696, 719, 745, 858],  # 热点残基编号（UniProt 编号）
    pocket_id="P1",                    # 指定口袋 ID（来自 structures.analyze()）

    # ── 扩散采样参数 / Diffusion sampling ───────────────────────────────────
    sampling_steps=100,                # 采样步数：50 / 100 / 150 / 200
    num_trajectories=3,                # 每条肽段折叠轨迹数

    # ── 自动折叠 / Auto-folding ──────────────────────────────────────────────
    auto_fold=True,                    # 自动折叠 Top N 肽段
    top_n_fold=10,                     # 折叠前 10 名

    # ── 去免疫原性 / Deimmunization ──────────────────────────────────────────
    deimmunize=False,                  # 启用去免疫原性优化
    immunogenicity_threshold=0.3,      # 免疫原性上限阈值

    # ── 支架模板 / Scaffold template ─────────────────────────────────────────
    scaffold_template=None,            # 模板序列（scaffold_templated 策略时必填）
    scaffold_conservation=0.6,         # 模板保守度（0-1，1=完全保守）

    # ── 药物偶联 / Drug conjugate ────────────────────────────────────────────
    conjugate_site="C_terminus",       # 偶联位点：N_terminus / C_terminus / Cys / Lys
    linker_compatible=True,            # 保留连接子兼容位点

    # ── 溶解度约束 / Solubility constraints ──────────────────────────────────
    min_solubility=-1.0,               # 最低 GRAVY 指数（负值=亲水）
    avoid_aggregation=True,            # 避免聚集倾向序列

    # ── 输出控制 / Output control ────────────────────────────────────────────
    sort_by="delta_g",                 # 排序字段：delta_g / ipsae / kd / composite
    include_structure=False,           # 是否在结果中包含 PDB 坐标
)
```

**完整参数参考表：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `gene` | `str` | 必填 | HGNC 基因符号 |
| `num_peptides` | `int` | `100` | 生成数量（受 Tier 限制） |
| `peptide_length` | `tuple[int,int]` | `(8, 20)` | 肽段长度范围（氨基酸数） |
| `targeting_strategy` | `str` | `"pocket_targeted"` | 靶向策略 |
| `target_residues` | `list[int]` | `None` | 热点残基编号列表 |
| `pocket_id` | `str` | `None` | 指定口袋 ID |
| `sampling_steps` | `int` | `50` | 扩散采样步数（50/100/150/200） |
| `num_trajectories` | `int` | `1` | 每条肽段折叠轨迹数 |
| `auto_fold` | `bool` | `False` | 自动折叠 Top N 肽段 |
| `top_n_fold` | `int` | `10` | 自动折叠数量 |
| `deimmunize` | `bool` | `False` | 启用去免疫原性优化 |
| `immunogenicity_threshold` | `float` | `0.3` | 免疫原性上限 |
| `scaffold_template` | `str` | `None` | 支架模板序列 |
| `scaffold_conservation` | `float` | `0.6` | 模板保守度 |
| `conjugate_site` | `str` | `None` | 偶联位点 |
| `min_solubility` | `float` | `None` | 最低 GRAVY 指数 |
| `avoid_aggregation` | `bool` | `True` | 避免聚集序列 |
| `sort_by` | `str` | `"delta_g"` | 结果排序字段 |
| `include_structure` | `bool` | `False` | 结果包含 PDB 坐标 |

---

## 4. `targeting_strategy` 详解

### 4.1 `pocket_targeted` — 口袋靶向策略

**适用场景：** 已知靶点三维结构，且有明确结合口袋的受体（如 GPCR、激酶、离子通道）。

**工作原理：** LigandForge 读取靶点结构的口袋几何信息，在扩散过程中施加空间约束，
使生成的肽段序列在折叠后能与指定口袋形成互补接触。

```python
# 口袋靶向示例：指定 EGFR 激酶域热点残基
job = client.peptides.generate(
    gene="EGFR",
    targeting_strategy="pocket_targeted",
    target_residues=[696, 719, 745, 858],  # ATP 结合口袋关键残基
    pocket_id="P1",                         # 来自 structures.analyze() 的口袋 ID
    num_peptides=300,
    sampling_steps=100,
)
```

| 优势 | 限制 |
|------|------|
| 高结合特异性 | 需要结构信息（PDB 或 AlphaFold） |
| 可精确控制结合位点 | 计算成本较高 |
| 适合竞争性抑制设计 | 对构象变化敏感 |

### 4.2 `scaffold_templated` — 支架模板策略

**适用场景：** 已有活性肽段或天然配体序列，希望在保留核心药效团的基础上优化性质。

```python
# 支架模板示例：基于 EGF 片段优化
job = client.peptides.generate(
    gene="EGFR",
    targeting_strategy="scaffold_templated",
    scaffold_template="CMHIESLDSYTC",  # EGF 受体结合域片段
    scaffold_conservation=0.5,          # 50% 保守度，允许较大变化
    num_peptides=200,
    deimmunize=True,                    # 同时去免疫原性
)
```

| 优势 | 限制 |
|------|------|
| 保留已知活性核心 | 需要已知活性序列 |
| 提高命中率 | 创新性受限 |
| 可与去免疫原性联用 | 可能陷入局部最优 |

### 4.3 `deimmunized` — 去免疫原性策略

**适用场景：** 治疗性肽段，需要降低 MHC-II 呈递风险，减少免疫反应。

```python
# 去免疫原性示例
job = client.peptides.generate(
    gene="ASGR1",
    targeting_strategy="deimmunized",
    deimmunize=True,
    immunogenicity_threshold=0.2,  # 严格阈值
    num_peptides=200,
    sampling_steps=150,            # 更多步数以满足约束
)
```

| 优势 | 限制 |
|------|------|
| 降低免疫原性风险 | 可能降低结合亲和力 |
| 适合慢性病治疗 | 生成多样性降低 |
| 符合 FDA 指导原则 | 需要更多采样步数 |

### 4.4 `drug_conjugate` — 药物偶联策略

**适用场景：** 抗体药物偶联物（ADC）替代方案，肽段-药物偶联物（PDC）载体设计。

```python
# 药物偶联示例：C 端偶联
job = client.peptides.generate(
    gene="HER2",
    targeting_strategy="drug_conjugate",
    conjugate_site="C_terminus",    # 偶联位点
    linker_compatible=True,         # 保留连接子兼容位点
    num_peptides=150,
)
```

---

## 5. SSE 流式进度监听

LigandForge 支持通过 Server-Sent Events (SSE) 实时监听生成进度，适合长时间任务的进度追踪。

### 事件类型说明

| 事件类型 | 触发时机 | 包含数据 |
|---------|---------|---------|
| `generation_start` | 任务开始 | `job_id`, `total_peptides`, `estimated_time` |
| `batch_complete` | 每批次生成完成 | `batch_id`, `peptides_generated`, `progress_pct` |
| `elite_hit` | 发现高分候选（ΔG < -10） | `sequence`, `delta_g`, `rank` |
| `fold_start` | 开始折叠某条肽段 | `sequence`, `peptide_rank` |
| `fold_complete` | 单条肽段折叠完成 | `sequence`, `ipsae`, `iptm`, `pdb_url` |
| `generation_complete` | 全部任务完成 | `total_generated`, `top_peptide`, `job_summary` |
| `error` | 发生错误 | `error_code`, `message`, `recoverable` |

### 流式监听代码示例

```python
import ligandai
from datetime import datetime

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# 启动生成任务
job = client.peptides.generate(
    gene="EGFR",
    num_peptides=300,
    targeting_strategy="pocket_targeted",
    auto_fold=True,
    top_n_fold=10,
    sampling_steps=100,
)

print(f"Job ID: {job.job_id}")
print(f"Estimated time: {job.estimated_time}s")
print("Streaming progress...\n")

elite_hits = []
folded_peptides = []

# 流式监听所有事件
for event in job.stream():
    ts = datetime.now().strftime("%H:%M:%S")
    
    if event.type == "generation_start":
        print(f"[{ts}] 🚀 Generation started | Total: {event.total_peptides} peptides")
    
    elif event.type == "batch_complete":
        print(f"[{ts}] 📦 Batch {event.batch_id} done | "
              f"Progress: {event.progress_pct:.1f}% "
              f"({event.peptides_generated}/{event.total_peptides})")
    
    elif event.type == "elite_hit":
        elite_hits.append(event)
        print(f"[{ts}] ⭐ ELITE HIT! {event.sequence} | "
              f"ΔG={event.delta_g:.2f} kcal/mol | Rank #{event.rank}")
    
    elif event.type == "fold_start":
        print(f"[{ts}] 🔬 Folding: {event.sequence[:15]}... (rank #{event.peptide_rank})")
    
    elif event.type == "fold_complete":
        folded_peptides.append(event)
        confidence = "HIGH" if event.ipsae < 0.3 else "MED" if event.ipsae < 0.5 else "LOW"
        print(f"[{ts}] ✅ Folded: {event.sequence[:15]}... | "
              f"iPSAE={event.ipsae:.3f} [{confidence}] | "
              f"iPTM={event.iptm:.3f}")
    
    elif event.type == "generation_complete":
        print(f"\n[{ts}] 🎉 Generation complete!")
        print(f"  Total generated: {event.total_generated}")
        print(f"  Elite hits: {len(elite_hits)}")
        print(f"  Folded: {len(folded_peptides)}")
        print(f"  Best peptide: {event.top_peptide.sequence}")
        print(f"  Best ΔG: {event.top_peptide.delta_g:.2f} kcal/mol")
    
    elif event.type == "error":
        print(f"[{ts}] ❌ Error: {event.message} (recoverable={event.recoverable})")
        if not event.recoverable:
            break
```

---

## 6. 完整端到端代码示例

```python
"""
端到端示例：靶点发现 → 结构解析 → 生成300肽段 → 流式折叠 → 输出Top10
End-to-end: Target discovery → Structure → Generate 300 peptides → Stream fold → Top 10
"""

import ligandai
import pandas as pd

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# ── Step 1: 靶点发现 / Target discovery ──────────────────────────────────────
print("="*60)
print("STEP 1: Target Discovery")
print("="*60)

targets = client.discovery.tissue_markers(
    tissue="liver",
    receptor_only=True,
    si_threshold=3.0,
    top_n=5,
)
best_gene = targets[0].gene
print(f"✓ Best target: {best_gene} (SI={targets[0].si:.2f})")

# ── Step 2: 结构解析与口袋分析 / Structure & pocket analysis ─────────────────
print("\nSTEP 2: Structure & Pocket Analysis")
print("="*60)

structure = client.structures.get(gene=best_gene)
print(f"✓ Structure source: {structure.source}")  # e.g. "ReceptorDB" / "PDB" / "AlphaFold"

pocket_analysis = client.structures.analyze(
    gene=best_gene,
    analysis_depth="standard",
)
best_pocket = pocket_analysis.recommended_pocket
print(f"✓ Best pocket: {best_pocket.pocket_id}")
print(f"  Volume: {best_pocket.volume:.1f} Å³")
print(f"  Druggability: {best_pocket.druggability_score:.2f}")
print(f"  Key residues: {best_pocket.key_residues[:5]}")

# ── Step 3: 生成300条肽段 / Generate 300 peptides ────────────────────────────
print("\nSTEP 3: Peptide Generation (300 peptides)")
print("="*60)

job = client.peptides.generate(
    gene=best_gene,
    num_peptides=300,
    targeting_strategy="pocket_targeted",
    target_residues=best_pocket.key_residues,
    pocket_id=best_pocket.pocket_id,
    sampling_steps=100,
    num_trajectories=2,
    auto_fold=True,
    top_n_fold=10,
    deimmunize=False,
    avoid_aggregation=True,
    min_solubility=-1.5,
    sort_by="delta_g",
)

print(f"✓ Job submitted: {job.job_id}")
print(f"  Estimated time: {job.estimated_time}s")
print(f"  Credit cost: {job.estimated_credits:,} credits")

# ── Step 4: 流式监听进度 / Stream progress ───────────────────────────────────
print("\nSTEP 4: Streaming Progress")
print("="*60)

results = None
for event in job.stream():
    if event.type == "elite_hit":
        print(f"  ⭐ Elite: {event.sequence} | ΔG={event.delta_g:.2f}")
    elif event.type == "fold_complete":
        print(f"  ✅ Folded: {event.sequence[:12]}... | iPSAE={event.ipsae:.3f}")
    elif event.type == "generation_complete":
        results = event
        print(f"\n✓ Complete! {event.total_generated} peptides generated")

# ── Step 5: 获取完整结果 / Get full results ───────────────────────────────────
print("\nSTEP 5: Top 10 Results")
print("="*60)

final = job.results()
df = pd.DataFrame([
    {
        "rank":           i + 1,
        "sequence":       p.sequence,
        "length":         len(p.sequence),
        "delta_g":        round(p.delta_g, 2),
        "kd_nM":          round(p.kd, 2),
        "ipsae":          round(p.ipsae, 3) if p.ipsae else None,
        "iptm":           round(p.iptm, 3) if p.iptm else None,
        "immunogenicity": round(p.immunogenicity, 3),
        "solubility":     round(p.solubility, 3),
        "folded":         p.ipsae is not None,
    }
    for i, p in enumerate(final.peptides[:10])
])

print(df.to_string(index=False))

# 保存结果
df.to_csv(f"{best_gene}_top10_peptides.csv", index=False)
print(f"\n✓ Saved to {best_gene}_top10_peptides.csv")
print(f"\n→ Next: Use ligandai-scoring for detailed thermodynamic analysis")
print(f"→ Next: Use ligandai-folding for additional trajectory validation")
```

---

## 7. 结果字段说明

| 字段名 | 类型 | 说明 | 优质阈值 |
|--------|------|------|---------|
| `sequence` | `str` | 肽段氨基酸序列（单字母码） | — |
| `length` | `int` | 肽段长度（氨基酸数） | 8-20 aa |
| `delta_g` | `float` | 结合自由能（kcal/mol），越负越强 | < -8 kcal/mol |
| `kd` | `float` | 预测解离常数（nM），越小越强 | < 100 nM |
| `ipsae` | `float` | 界面预测对齐误差（折叠后可用） | < 0.3 |
| `iptm` | `float` | 界面 TM 分数（折叠后可用） | > 0.75 |
| `immunogenicity` | `float` | 免疫原性风险（0-1） | < 0.3 |
| `solubility` | `float` | GRAVY 溶解度指数 | < 0（亲水） |
| `stability` | `float` | 热稳定性预测（Pro+） | > 0.6 |
| `half_life_h` | `float` | 血浆半衰期（小时，Pro+） | > 2 h |
| `pdb_url` | `str` | 折叠结构 PDB 下载链接 | — |
| `composite_score` | `float` | 综合评分（LigandIQ 计算） | > 0.7 |

---

## 8. Tier 限制表

| Tier | 最大 num_peptides | 自动折叠 | 去免疫原性 | 药物偶联 | 支架模板 |
|------|------------------|---------|-----------|---------|---------|
| Free | 100 | ❌ | ❌ | ❌ | ❌ |
| Basic | 300 | ❌ | ✅ | ❌ | ✅ |
| Academia | 300 | ✅ (top 5) | ✅ | ❌ | ✅ |
| Pro | 1,000 | ✅ (top 20) | ✅ | ✅ | ✅ |
| Enterprise | 5,000 | ✅ (top 50) | ✅ | ✅ | ✅ |

> **注意：** `num_trajectories > 1` 需要 Academia+ Tier。

---

## 9. Credit 成本

| 操作 | Credit 成本 | 说明 |
|------|------------|------|
| 肽段生成（无折叠） | 100 credits/peptide | 基础生成费用 |
| 自动折叠（50步，1轨迹） | +100 credits/peptide | 折叠基础费用 |
| 自动折叠（100步，1轨迹） | +200 credits/peptide | 步数翻倍 |
| 自动折叠（多轨迹） | ×N credits | N = num_trajectories |
| 去免疫原性优化 | +20 credits/peptide | 额外优化费用 |

**示例计算：**
```
300 peptides × 100 credits = 30,000 credits (generation)
+ 10 peptides × 200 credits (100-step fold, 1 trajectory) = 2,000 credits
Total: 32,000 credits ≈ $32 USD (at standard rate)
```

---

## 10. 最佳实践 Tips

### 如何选择 `num_peptides`

```
探索阶段（新靶点）：  300-500 条，宽泛采样
优化阶段（已知靶点）：100-200 条，聚焦特定口袋
验证阶段（候选确认）：50-100 条，高质量生成
```

### 如何选择 `sampling_steps`

| 步数 | 质量 | 时间 | 推荐场景 |
|------|------|------|---------|
| 50 | 基础 | ~2 min/100 | 快速探索，大批量筛选 |
| 100 | 良好 | ~4 min/100 | 标准设计流程 |
| 150 | 优秀 | ~6 min/100 | 重要靶点精细设计 |
| 200 | 最优 | ~8 min/100 | 临床前候选优化 |

### 如何指定热点残基（`target_residues`）

```python
# 方法1：从口袋分析自动获取
pocket = client.structures.analyze(gene="EGFR", analysis_depth="standard")
residues = pocket.recommended_pocket.key_residues  # 自动推荐

# 方法2：从文献/实验数据手动指定
residues = [696, 719, 745, 858]  # EGFR ATP 结合口袋

# 方法3：从 PDB 结构分析
# 使用 PyMOL 或 PLIP 分析已知配体-受体复合物
```

### 溶解度与免疫原性平衡

```python
# 治疗性肽段推荐配置
job = client.peptides.generate(
    gene="TARGET",
    num_peptides=200,
    targeting_strategy="deimmunized",
    deimmunize=True,
    immunogenicity_threshold=0.25,  # 严格免疫原性控制
    min_solubility=-1.0,            # 确保水溶性
    avoid_aggregation=True,
    sampling_steps=150,             # 更多步数满足约束
)
```

---

## 11. 下游衔接

```python
# 生成完成后的推荐下游流程
# Recommended downstream workflow after generation

# → ligandai-scoring: 热力学精评（免费，CPU-only）
scored = client.scoring.deltaforge(
    peptides=[p.sequence for p in final.peptides[:50]],
    gene="EGFR",
)

# → ligandai-folding: 额外折叠轨迹验证
for peptide in final.peptides[:5]:
    fold_result = client.structures.fold_complex(
        receptor_gene="EGFR",
        peptide_sequence=peptide.sequence,
        num_trajectories=5,
        sampling_steps=200,
    )

# → ligandai-bivalent: 双价设计（Pro+）
bivalent = client.bivalent.start(
    target_a="EGFR",
    target_b="HER2",
    linker_type="GGS",
    linker_length=3,
)
```

**推荐工作流：**
```
ligandai-peptide (生成300条)
    ↓
ligandai-scoring/LigandIQ (预筛，免费)
    ↓ Top 50
ligandai-folding (精折叠，100 credits/trajectory)
    ↓ Top 10
ligandai-scoring/DeltaForge (精评，免费)
    ↓ Top 3-5
实验验证 (SPR / ITC / 细胞实验)
```

> **技能版本：** v1.0.0 | **LigandForge 版本：** v6.5 | **SDK 最低版本：** ligandai>=2.1.0
