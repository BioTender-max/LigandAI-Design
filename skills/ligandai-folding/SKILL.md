---
name: ligandai-folding
description: |
  【中文】基于 LigandAI Boltz-2 引擎的蛋白质复合物折叠技能。支持受体结构解析（ReceptorDB/PDB/AlphaFold 三级回退）、
  结合口袋分析、肽段-受体复合物折叠，提供 iPSAE/iPTM 置信度指标，支持 MSA 缓存加速、糖基化修饰和多轨迹采样。
  
  [EN] Protein complex folding skill powered by LigandAI Boltz-2 engine. Supports receptor structure resolution
  (ReceptorDB/PDB/AlphaFold three-tier fallback), binding pocket analysis, and peptide-receptor complex folding.
  Provides iPSAE/iPTM confidence metrics with MSA caching, glycosylation support, and multi-trajectory sampling.
license: Apache-2.0
category: design-tools
tags:
  - folding
  - boltz2
  - complex
  - msa
  - glycosylation
  - structure

triggers:
  # 中文触发场景
  - "折叠肽段与受体的复合物结构"
  - "用 Boltz-2 预测蛋白质结构"
  - "分析 EGFR 的结合口袋"
  - "解析靶点三维结构"
  - "计算肽段结合的 iPSAE 置信度"
  - "多轨迹折叠采样"
  - "糖基化蛋白质结构预测"
  - "获取靶点结构并识别药物口袋"
  # English triggers
  - "fold peptide-receptor complex with Boltz-2"
  - "analyze binding pocket of EGFR"
  - "resolve target 3D structure"
  - "calculate iPSAE confidence for peptide binding"
  - "multi-trajectory folding sampling"
  - "predict glycosylated protein structure"
  - "get receptor structure and identify druggable pockets"

inputs:
  - name: gene
    type: string
    required: true
    description: "目标基因符号（HGNC）/ Target gene symbol (HGNC)"
  - name: peptide_sequence
    type: string
    required: false
    description: "待折叠的肽段序列（复合物折叠时必填）/ Peptide sequence for complex folding"
  - name: analysis_depth
    type: string
    required: false
    default: "standard"
    description: "口袋分析深度：quick / standard / full"
  - name: sampling_steps
    type: integer
    required: false
    default: 50
    description: "Boltz-2 扩散采样步数（50/100/150/200）/ Boltz-2 diffusion sampling steps"
  - name: num_trajectories
    type: integer
    required: false
    default: 1
    description: "折叠轨迹数量（多轨迹取最优）/ Number of folding trajectories"
  - name: template_mode
    type: string
    required: false
    default: "auto"
    description: "模板使用模式：auto / pdb_only / no_template / custom"
  - name: glycosylation
    type: boolean
    required: false
    default: false
    description: "是否考虑糖基化修饰 / Whether to model glycosylation"

outputs:
  - name: structure
    type: StructureResult
    description: "解析的受体结构信息 / Resolved receptor structure"
  - name: pocket_analysis
    type: PocketAnalysis
    description: "口袋分析结果 / Pocket analysis results"
  - name: complex_fold
    type: ComplexFoldResult
    description: "复合物折叠结果 / Complex folding result"
  - name: complex_fold.ipsae
    type: float
    description: "界面预测对齐误差（<0.3 高置信）/ Interface predicted alignment error"
  - name: complex_fold.iptm
    type: float
    description: "界面 TM 分数（>0.75 高置信）/ Interface TM score"
  - name: complex_fold.pdb_url
    type: string
    description: "折叠结构 PDB 文件下载链接 / PDB file download URL"
---

# LigandAI Folding — 复合物折叠技能

## 1. 功能概述

### 中文说明

LigandAI Folding 模块基于 **Boltz-2** 扩散模型，提供从受体结构解析到肽段-受体复合物折叠的完整结构生物学工作流。
Boltz-2 是目前精度最高的开源蛋白质复合物折叠模型之一，在 CASP15 复合物预测任务中超越 AlphaFold-Multimer，
特别擅长处理肽段-蛋白质相互作用（PPI）预测。

### Boltz-2 核心特性表

| 特性 | 说明 | 优势 |
|------|------|------|
| **MSA 缓存** | 多序列比对结果缓存复用 | 同一靶点重复折叠速度提升 5-10× |
| **糖基化建模** | N-/O-糖基化位点显式建模 | 提高糖蛋白结构预测精度 |
| **模板模式** | 支持 PDB 模板引导折叠 | 已知结构靶点精度更高 |
| **多轨迹采样** | 并行生成多个构象 | 捕获构象多样性，取最优 |
| **复合物折叠** | 肽段+受体联合折叠 | 直接预测结合构象 |
| **置信度指标** | iPSAE + iPTM 双指标 | 量化界面预测可靠性 |
| **GPU 加速** | A100/H100 GPU 集群 | 50步折叠约 3-5 分钟 |

---

## 2. 结构解析流程

LigandAI 采用三级回退策略解析靶点结构，确保最高质量的起始结构：

```
┌─────────────────────────────────────────────────────────────┐
│                   结构解析优先级流程                          │
│              Structure Resolution Priority                   │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  1. ReceptorDB (LigandAI 专有数据库)                        │
│     ├── 经过人工审核的高质量结构                              │
│     ├── 预计算的口袋信息和 MSA                               │
│     └── 覆盖 >8,000 个人类受体                              │
│                    ↓ (未找到时)                              │
│  2. RCSB PDB (实验结构)                                     │
│     ├── 优先选择分辨率 < 2.5 Å 的晶体结构                   │
│     ├── 自动选择最相关的 PDB 条目                            │
│     └── 支持 X-ray / Cryo-EM / NMR                         │
│                    ↓ (未找到时)                              │
│  3. AlphaFold2 预测结构                                     │
│     ├── EBI AlphaFold 数据库                                │
│     ├── pLDDT > 70 区域可信                                 │
│     └── 覆盖几乎所有人类蛋白质                               │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

```python
# 查看结构来源
structure = client.structures.get(gene="EGFR")
print(f"Source: {structure.source}")          # "ReceptorDB" / "PDB" / "AlphaFold"
print(f"PDB ID: {structure.pdb_id}")          # e.g. "1IVO"
print(f"Resolution: {structure.resolution}")  # e.g. "1.9 Å" (PDB only)
print(f"Coverage: {structure.coverage:.1%}")  # 序列覆盖率
print(f"MSA cached: {structure.msa_cached}")  # MSA 是否已缓存
```

---

## 3. API 参数详解

### 3.1 `client.structures.get()` — 获取受体结构

```python
structure = client.structures.get(
    gene="EGFR",                    # 目标基因（HGNC 符号）
    prefer_source="ReceptorDB",     # 优先数据源：ReceptorDB / PDB / AlphaFold
    chain="A",                      # 指定 PDB 链（可选）
    pdb_id="1IVO",                  # 指定 PDB ID（可选，覆盖自动选择）
    include_ligands=False,          # 是否包含共晶配体
    trim_disordered=True,           # 修剪无序区域（pLDDT < 50）
)
```

**参数详解：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `gene` | `str` | 必填 | HGNC 基因符号 |
| `prefer_source` | `str` | `"ReceptorDB"` | 优先数据源 |
| `chain` | `str` | `None` | 指定 PDB 链 ID |
| `pdb_id` | `str` | `None` | 指定 PDB 条目 ID |
| `include_ligands` | `bool` | `False` | 包含共晶配体坐标 |
| `trim_disordered` | `bool` | `True` | 修剪低置信度无序区域 |
| `domain` | `str` | `None` | 指定结构域，如 `"kinase_domain"` |

### 3.2 `client.structures.analyze()` — 口袋分析

```python
pocket_analysis = client.structures.analyze(
    gene="EGFR",
    analysis_depth="standard",      # quick / standard / full
    pocket_detection="fpocket",     # fpocket / sitemap / dogsitescorer
    min_volume=100,                 # 最小口袋体积（Å³）
    min_druggability=0.5,           # 最低可成药性评分
    include_allosteric=True,        # 包含变构口袋
    surface_probe_radius=1.4,       # 探针半径（Å）
)
```

**参数详解：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `gene` | `str` | 必填 | HGNC 基因符号 |
| `analysis_depth` | `str` | `"standard"` | 分析深度（见下节） |
| `pocket_detection` | `str` | `"fpocket"` | 口袋检测算法 |
| `min_volume` | `float` | `100` | 最小口袋体积（Å³） |
| `min_druggability` | `float` | `0.5` | 最低可成药性评分（0-1） |
| `include_allosteric` | `bool` | `True` | 包含变构口袋 |
| `surface_probe_radius` | `float` | `1.4` | 溶剂探针半径（Å） |

---

## 4. `analysis_depth` 选项详解

| 深度 | 耗时 | 功能 | 推荐场景 |
|------|------|------|---------|
| `quick` | ~30s | 口袋检测 + 体积计算 | 快速筛选，大批量靶点评估 |
| `standard` | ~2min | + 可成药性评分 + 关键残基识别 + 变构口袋 | 标准设计流程（推荐） |
| `full` | ~10min | + 分子动力学口袋稳定性 + 水分子分析 + 结合能预测 | 临床前候选精细分析 |

```python
# quick 模式：快速获取口袋列表
quick = client.structures.analyze(gene="EGFR", analysis_depth="quick")
print(f"Found {len(quick.pockets)} pockets")

# standard 模式：标准分析（推荐）
standard = client.structures.analyze(gene="EGFR", analysis_depth="standard")
for pocket in standard.pockets:
    print(f"Pocket {pocket.pocket_id}: "
          f"Vol={pocket.volume:.0f}Å³, "
          f"Drug={pocket.druggability_score:.2f}")

# full 模式：完整分析（临床前）
full = client.structures.analyze(gene="EGFR", analysis_depth="full")
print(f"MD stability: {full.recommended_pocket.md_stability:.2f}")
print(f"Water bridges: {full.recommended_pocket.water_bridges}")
```

---

## 5. `recommended_pocket` 解读

```python
pocket = pocket_analysis.recommended_pocket

# 基本几何信息
print(f"Pocket ID:          {pocket.pocket_id}")        # e.g. "P1"
print(f"Residue range:      {pocket.range}")             # e.g. "695-760, 820-870"
print(f"Volume:             {pocket.volume:.1f} Å³")    # e.g. 842.3 Å³
print(f"Surface area:       {pocket.surface_area:.1f} Å²")

# 可成药性评分
print(f"Druggability score: {pocket.druggability_score:.2f}")  # 0-1, >0.7 = druggable
print(f"Hydrophobicity:     {pocket.hydrophobicity:.2f}")      # 疏水性
print(f"Polarity:           {pocket.polarity:.2f}")            # 极性

# 关键残基
print(f"Key residues:       {pocket.key_residues}")    # [696, 719, 745, 858]
print(f"Hotspot residues:   {pocket.hotspot_residues}") # 最重要的热点残基

# 口袋类型
print(f"Pocket type:        {pocket.pocket_type}")     # "orthosteric" / "allosteric"
print(f"Known ligands:      {pocket.known_ligands}")   # 已知结合配体列表
```

**口袋体积解读指南：**

| 体积范围 | 适合肽段长度 | 说明 |
|---------|------------|------|
| < 300 Å³ | 5-8 aa | 小口袋，短肽或拟肽 |
| 300-600 Å³ | 8-12 aa | 中等口袋，标准肽段 |
| 600-1000 Å³ | 12-18 aa | 大口袋，长肽或环肽 |
| > 1000 Å³ | 15-25 aa | 超大口袋，蛋白-蛋白界面 |

---

## 6. 折叠参数详解

### `client.structures.fold_complex()` — 复合物折叠

```python
fold_result = client.structures.fold_complex(
    # ── 必填参数 / Required ──────────────────────────────────────────────────
    receptor_gene="EGFR",
    peptide_sequence="ACDEFGHIKLMNPQRSTVWY",

    # ── 采样参数 / Sampling ──────────────────────────────────────────────────
    sampling_steps=100,             # 扩散步数：50 / 100 / 150 / 200
    num_trajectories=3,             # 轨迹数（取最优 iPSAE）

    # ── 模板模式 / Template mode ─────────────────────────────────────────────
    template_mode="auto",           # auto / pdb_only / no_template / custom
    template_pdb_id="1IVO",         # 自定义模板 PDB ID（template_mode="custom" 时）

    # ── 糖基化 / Glycosylation ───────────────────────────────────────────────
    glycosylation=False,            # 是否建模糖基化
    glycan_sites=None,              # 指定糖基化位点，如 [{"residue": 128, "type": "N-linked"}]

    # ── MSA 设置 / MSA settings ──────────────────────────────────────────────
    use_msa_cache=True,             # 使用缓存 MSA（加速）
    msa_depth=512,                  # MSA 深度（序列数）

    # ── 输出设置 / Output ────────────────────────────────────────────────────
    return_all_trajectories=False,  # 返回所有轨迹（否则仅返回最优）
    include_plddt=True,             # 包含 pLDDT 逐残基置信度
)
```

**折叠参数参考表：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `receptor_gene` | `str` | 必填 | 受体基因符号 |
| `peptide_sequence` | `str` | 必填 | 肽段氨基酸序列 |
| `sampling_steps` | `int` | `50` | 扩散采样步数 |
| `num_trajectories` | `int` | `1` | 折叠轨迹数 |
| `template_mode` | `str` | `"auto"` | 模板使用模式 |
| `template_pdb_id` | `str` | `None` | 自定义模板 PDB ID |
| `glycosylation` | `bool` | `False` | 建模糖基化修饰 |
| `glycan_sites` | `list` | `None` | 指定糖基化位点 |
| `use_msa_cache` | `bool` | `True` | 使用 MSA 缓存 |
| `msa_depth` | `int` | `512` | MSA 序列深度 |
| `return_all_trajectories` | `bool` | `False` | 返回所有轨迹 |
| `include_plddt` | `bool` | `True` | 包含 pLDDT 置信度 |

---

## 7. 完整代码示例：EGFR 结构解析与复合物折叠

```python
"""
完整示例：解析 EGFR 结构 → 分析口袋 → 折叠肽段复合物 → 读取 iPSAE/iPTM
Full example: Resolve EGFR structure → Analyze pocket → Fold complex → Read iPSAE/iPTM
"""

import ligandai
import pandas as pd

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# ── Step 1: 解析 EGFR 结构 / Resolve EGFR structure ──────────────────────────
print("="*60)
print("STEP 1: Resolving EGFR Structure")
print("="*60)

structure = client.structures.get(
    gene="EGFR",
    prefer_source="ReceptorDB",
    trim_disordered=True,
)

print(f"✓ Source:      {structure.source}")
print(f"  PDB ID:      {structure.pdb_id}")
print(f"  Resolution:  {structure.resolution}")
print(f"  Chain:       {structure.chain}")
print(f"  Length:      {structure.length} residues")
print(f"  MSA cached:  {structure.msa_cached}")
print(f"  Coverage:    {structure.coverage:.1%}")

# ── Step 2: 口袋分析 / Pocket analysis ───────────────────────────────────────
print("\nSTEP 2: Pocket Analysis (standard depth)")
print("="*60)

pocket_analysis = client.structures.analyze(
    gene="EGFR",
    analysis_depth="standard",
    include_allosteric=True,
    min_druggability=0.5,
)

print(f"✓ Found {len(pocket_analysis.pockets)} pockets")
print("\nAll pockets:")
for p in pocket_analysis.pockets:
    print(f"  {p.pocket_id}: Vol={p.volume:.0f}Å³, "
          f"Drug={p.druggability_score:.2f}, "
          f"Type={p.pocket_type}")

best_pocket = pocket_analysis.recommended_pocket
print(f"\n✓ Recommended pocket: {best_pocket.pocket_id}")
print(f"  Volume:          {best_pocket.volume:.1f} Å³")
print(f"  Druggability:    {best_pocket.druggability_score:.2f}")
print(f"  Key residues:    {best_pocket.key_residues}")
print(f"  Hotspot:         {best_pocket.hotspot_residues}")
print(f"  Known ligands:   {best_pocket.known_ligands}")

# ── Step 3: 折叠肽段复合物 / Fold peptide complexes ──────────────────────────
print("\nSTEP 3: Folding Peptide-EGFR Complexes")
print("="*60)

# 待折叠的候选肽段（来自 ligandai-peptide 生成结果）
candidate_peptides = [
    "ACDEFGHIKLMN",
    "WQRSTVYACDEF",
    "MNPQRSTVWYAC",
    "GHIKLMNPQRST",
    "DEFGHIKLMNPQ",
]

fold_results = []

for i, peptide in enumerate(candidate_peptides):
    print(f"\n  Folding peptide {i+1}/{len(candidate_peptides)}: {peptide}")
    
    fold = client.structures.fold_complex(
        receptor_gene="EGFR",
        peptide_sequence=peptide,
        sampling_steps=100,
        num_trajectories=3,          # 3 轨迹取最优
        template_mode="auto",
        use_msa_cache=True,          # 利用缓存 MSA 加速
        glycosylation=False,
        include_plddt=True,
    )
    
    # 等待折叠完成
    fold.wait(timeout=600)
    
    result = {
        "peptide":       peptide,
        "ipsae":         fold.ipsae,
        "iptm":          fold.iptm,
        "plddt_mean":    fold.plddt_mean,
        "plddt_interface": fold.plddt_interface,
        "best_trajectory": fold.best_trajectory_id,
        "pdb_url":       fold.pdb_url,
        "credits_used":  fold.credits_used,
    }
    fold_results.append(result)
    
    # 置信度判断
    ipsae_level = "HIGH" if fold.ipsae < 0.3 else "MED" if fold.ipsae < 0.5 else "LOW"
    iptm_level  = "HIGH" if fold.iptm > 0.75 else "MED" if fold.iptm > 0.5 else "LOW"
    
    print(f"    iPSAE: {fold.ipsae:.3f} [{ipsae_level}]")
    print(f"    iPTM:  {fold.iptm:.3f} [{iptm_level}]")
    print(f"    pLDDT (interface): {fold.plddt_interface:.1f}")
    print(f"    PDB: {fold.pdb_url}")

# ── Step 4: 汇总结果 / Summarize results ─────────────────────────────────────
print("\nSTEP 4: Summary")
print("="*60)

df = pd.DataFrame(fold_results)
df["ipsae_confidence"] = df["ipsae"].apply(
    lambda x: "HIGH" if x < 0.3 else ("MED" if x < 0.5 else "LOW")
)
df["iptm_confidence"] = df["iptm"].apply(
    lambda x: "HIGH" if x > 0.75 else ("MED" if x > 0.5 else "LOW")
)
df_sorted = df.sort_values("ipsae").reset_index(drop=True)
df_sorted["rank"] = df_sorted.index + 1

print(df_sorted[["rank", "peptide", "ipsae", "ipsae_confidence",
                  "iptm", "iptm_confidence", "plddt_interface"]].to_string(index=False))

# 保存结果
df_sorted.to_csv("EGFR_fold_results.csv", index=False)
print(f"\n✓ Results saved to EGFR_fold_results.csv")

# 最优候选
best = df_sorted.iloc[0]
print(f"\n✓ Best candidate: {best['peptide']}")
print(f"  iPSAE={best['ipsae']:.3f}, iPTM={best['iptm']:.3f}")
print(f"  PDB: {best['pdb_url']}")
print(f"\n→ Next: Use ligandai-scoring for thermodynamic evaluation")
```

---

## 8. 置信度指标解读表

### iPSAE（界面预测对齐误差）

| iPSAE 值 | 置信度 | 解读 | 建议 |
|---------|--------|------|------|
| < 0.20 | 极高置信 | 界面结构高度可靠 | 直接用于后续分析 |
| 0.20 - 0.30 | 高置信 | 界面结构可靠 | 推荐进入实验验证 |
| 0.30 - 0.40 | 中等置信 | 界面结构基本可靠 | 建议增加轨迹数验证 |
| 0.40 - 0.50 | 低置信 | 界面结构不确定 | 需要额外验证 |
| > 0.50 | 极低置信 | 界面结构不可靠 | 不推荐，考虑重新设计 |

### iPTM（界面 TM 分数）

| iPTM 值 | 置信度 | 解读 | 建议 |
|---------|--------|------|------|
| > 0.85 | 极高置信 | 复合物拓扑结构高度可靠 | 优先候选 |
| 0.75 - 0.85 | 高置信 | 复合物拓扑结构可靠 | 推荐候选 |
| 0.60 - 0.75 | 中等置信 | 复合物拓扑基本正确 | 可接受，需验证 |
| 0.50 - 0.60 | 低置信 | 拓扑结构不确定 | 谨慎使用 |
| < 0.50 | 极低置信 | 拓扑结构不可靠 | 不推荐 |

### 综合判断标准

```python
def assess_fold_quality(ipsae: float, iptm: float) -> str:
    """综合评估折叠质量 / Comprehensive fold quality assessment"""
    if ipsae < 0.3 and iptm > 0.75:
        return "EXCELLENT - 高置信复合物，推荐进入实验验证"
    elif ipsae < 0.4 and iptm > 0.6:
        return "GOOD - 中等置信，建议增加轨迹数"
    elif ipsae < 0.5 or iptm > 0.5:
        return "MARGINAL - 低置信，需要额外验证"
    else:
        return "POOR - 不可靠，建议重新设计肽段"

# 使用示例
quality = assess_fold_quality(ipsae=0.25, iptm=0.82)
print(quality)  # "EXCELLENT - 高置信复合物，推荐进入实验验证"
```

---

## 9. Tier 限制

| Tier | 每靶点折叠次数/月 | 最大轨迹数 | 最大采样步数 | 糖基化 | 全轨迹返回 |
|------|-----------------|-----------|------------|--------|-----------|
| Free | 1 | 1 | 50 | ❌ | ❌ |
| Basic | 25/target | 2 | 100 | ❌ | ❌ |
| Academia | 50/target | 3 | 150 | ✅ | ❌ |
| Pro | 100/target | 5 | 200 | ✅ | ✅ |
| Enterprise | 无限制 | 10 | 200 | ✅ | ✅ |

> **注意：** "每靶点"指同一基因的折叠次数限制，不同靶点独立计算。

---

## 10. Credit 成本

| 操作 | Credit 成本 | 计费规则 |
|------|------------|---------|
| 50步折叠，1轨迹 | 100 credits | 基础费用 |
| 100步折叠，1轨迹 | 200 credits | 步数按50递增，每+50步+100 credits |
| 150步折叠，1轨迹 | 300 credits | |
| 200步折叠，1轨迹 | 400 credits | |
| 多轨迹 | ×N credits | N = num_trajectories |
| MSA 缓存命中 | -20% 折扣 | 同一靶点第二次起享受折扣 |

**计费公式：**
```
Credits = (sampling_steps / 50) × 100 × num_trajectories × msa_discount
msa_discount = 0.8 if msa_cached else 1.0
```

**示例：**
```
100步 × 3轨迹 × MSA缓存 = 200 × 3 × 0.8 = 480 credits
```

---

## 11. GPU 并发槽位说明

LigandAI 折叠服务基于 GPU 集群运行，并发槽位按 Tier 分配：

| Tier | GPU 并发槽位 | 队列优先级 | 预计等待时间 |
|------|------------|-----------|------------|
| Free | 1 | 低 | 5-30 分钟 |
| Basic | 2 | 中 | 2-10 分钟 |
| Academia | 3 | 中高 | 1-5 分钟 |
| Pro | 5 | 高 | < 2 分钟 |
| Enterprise | 专属 GPU | 最高 | 即时 |

```python
# 查看当前队列状态
queue_status = client.infrastructure.queue_status()
print(f"Current wait time: {queue_status.estimated_wait_s}s")
print(f"Active jobs: {queue_status.active_jobs}")
print(f"Your slots: {queue_status.your_slots}/{queue_status.max_slots}")
```

---

## 12. ReceptorDB Share-back 折扣

将折叠结果贡献回 ReceptorDB 可享受 50% Credit 折扣：

```python
# 折叠时启用 share-back
fold = client.structures.fold_complex(
    receptor_gene="NOVEL_RECEPTOR",
    peptide_sequence="ACDEFGHIKLMN",
    sampling_steps=100,
    num_trajectories=3,
    contribute_to_receptordb=True,   # 启用 share-back，享受 50% 折扣
    # 注意：贡献的结构将在 6 个月后公开
)

print(f"Credits charged: {fold.credits_used}")   # 享受 50% 折扣
print(f"Contributed: {fold.contributed_to_db}")  # True
print(f"Public date: {fold.public_date}")         # 公开日期
```

**Share-back 规则：**
- 折扣：50% off（100 credits → 50 credits）
- 公开时间：贡献后 6 个月
- 适用条件：新靶点（ReceptorDB 中尚无记录）
- 质量要求：iPSAE < 0.4 且 iPTM > 0.6

> **技能版本：** v1.0.0 | **Boltz-2 版本：** v2.3 | **SDK 最低版本：** ligandai>=2.1.0
