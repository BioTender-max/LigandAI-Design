---
name: ligandai-scoring
description: |
  【中文】基于 LigandAI DeltaForge 热力学引擎和 LigandIQ 预折叠筛选器的多维度肽段评分技能。
  提供结合自由能（ΔG）、解离常数（Kd）、免疫原性、溶解度、稳定性等全面评分，
  支持两阶段筛选策略（LigandIQ 预筛 + DeltaForge 精评），所有评分均免费（CPU-only）。
  
  [EN] Multi-dimensional peptide scoring skill powered by LigandAI DeltaForge thermodynamic engine
  and LigandIQ pre-fold screener. Provides comprehensive scoring including binding free energy (ΔG),
  dissociation constant (Kd), immunogenicity, solubility, and stability. Supports two-stage screening
  (LigandIQ pre-screen + DeltaForge precision scoring). All scoring is free (CPU-only).
license: Apache-2.0
category: design-tools
tags:
  - scoring
  - deltaforge
  - thermodynamics
  - dG
  - Kd
  - immunogenicity
  - solubility

triggers:
  # 中文触发场景
  - "评估肽段的结合自由能"
  - "计算肽段的 Kd 解离常数"
  - "筛选低免疫原性肽段"
  - "DeltaForge 热力学评分"
  - "LigandIQ 预折叠筛选"
  - "多维度肽段综合评分"
  - "评估肽段溶解度和稳定性"
  - "从候选肽段中筛选最优"
  # English triggers
  - "score peptide binding free energy with DeltaForge"
  - "calculate Kd for peptide candidates"
  - "filter low-immunogenicity peptides"
  - "LigandIQ pre-fold screening"
  - "multi-dimensional peptide scoring"
  - "rank peptide candidates by thermodynamic score"
  - "assess peptide solubility and stability"

inputs:
  - name: peptides
    type: list[str]
    required: true
    description: "待评分的肽段序列列表 / List of peptide sequences to score"
  - name: gene
    type: string
    required: true
    description: "靶点基因符号 / Target gene symbol"
  - name: scoring_mode
    type: string
    required: false
    default: "deltaforge"
    description: "评分模式：deltaforge / ligandiq / both"
  - name: weights
    type: dict
    required: false
    description: "综合评分权重字典 / Weight dictionary for composite scoring"

outputs:
  - name: scores
    type: list[ScoringResult]
    description: "评分结果列表，按综合评分排序 / Scoring results sorted by composite score"
  - name: scores[].sequence
    type: string
    description: "肽段序列 / Peptide sequence"
  - name: scores[].delta_g
    type: float
    description: "结合自由能 kcal/mol / Binding free energy in kcal/mol"
  - name: scores[].kd
    type: float
    description: "解离常数 nM / Dissociation constant in nM"
  - name: scores[].immunogenicity
    type: float
    description: "免疫原性风险 0-1 / Immunogenicity risk 0-1"
  - name: scores[].solubility
    type: float
    description: "GRAVY 溶解度指数 / GRAVY solubility index"
---

# LigandAI Scoring — 热力学评分技能

## 1. 功能概述

### 中文说明

LigandAI Scoring 模块提供两个互补的评分引擎：

- **DeltaForge**：基于物理的热力学评分引擎，通过 MM-GBSA 近似计算结合自由能（ΔG），
  与实验 SPR/ITC 数据相关系数 r=0.83，是目前计算肽段结合亲和力最准确的免费工具之一。

- **LigandIQ**：轻量级预折叠筛选器，无需折叠即可快速评估肽段的多项理化性质，
  用于在折叠前淘汰明显不合格的候选，节省昂贵的折叠 Credits。

### 评分体系总览表

| 评分指标 | 引擎 | 单位 | 计算方式 | 是否免费 |
|---------|------|------|---------|---------|
| `delta_g` | DeltaForge | kcal/mol | MM-GBSA 热力学计算 | ✅ 免费 |
| `kd` | DeltaForge | nM | 由 ΔG 换算 | ✅ 免费 |
| `ipsae` | Boltz-2 | 无量纲 | 折叠后界面误差 | 需折叠 Credits |
| `iptm` | Boltz-2 | 无量纲 | 折叠后 TM 分数 | 需折叠 Credits |
| `immunogenicity` | LigandIQ | 0-1 | MHC-II 呈递预测 | ✅ 免费 |
| `solubility` | LigandIQ | GRAVY | 亲疏水性指数 | ✅ 免费 |
| `stability` | DeltaForge | 0-1 | 热稳定性预测 | Pro+ |
| `half_life` | DeltaForge | 小时 | 血浆半衰期预测 | Pro+ |
| `cross_species` | DeltaForge | 0-1 | 跨物种脱靶风险 | Pro+ |

---

## 2. 各评分指标深度解析

### 2.1 `delta_g` — 结合自由能

**物理意义：** 肽段与受体结合过程的吉布斯自由能变化，负值表示自发结合，绝对值越大结合越强。

**计算方法：** DeltaForge 采用 MM-GBSA（分子力学-广义玻恩表面积）近似：
```
ΔG_bind = ΔH_MM + ΔG_solvation - TΔS_conf
         = ΔE_vdW + ΔE_elec + ΔG_GB + ΔG_SA - TΔS
```

**实验验证：** 与 SPR/ITC 实验数据相关系数 r=0.83（n=1,247 个肽段-蛋白质复合物）

**阈值建议：**

| ΔG 值 (kcal/mol) | 结合强度 | 对应 Kd | 推荐用途 |
|-----------------|---------|---------|---------|
| < -12 | 超强结合 | < 1 nM | 治疗性候选，优先推进 |
| -10 ~ -12 | 强结合 | 1-10 nM | 优质候选 |
| -8 ~ -10 | 中强结合 | 10-100 nM | 良好候选 |
| -6 ~ -8 | 中等结合 | 100 nM - 1 μM | 可接受，需优化 |
| > -6 | 弱结合 | > 1 μM | 不推荐 |

### 2.2 `kd` — 解离常数

**换算公式：**
```python
import math

def delta_g_to_kd(delta_g_kcal_mol: float, temperature_K: float = 310.15) -> float:
    """
    将结合自由能转换为解离常数
    Convert binding free energy to dissociation constant
    
    ΔG = RT × ln(Kd)  →  Kd = exp(ΔG / RT)
    
    R = 1.987 cal/(mol·K) = 0.001987 kcal/(mol·K)
    T = 310.15 K (37°C, 生理温度)
    """
    R = 0.001987  # kcal/(mol·K)
    kd_molar = math.exp(delta_g_kcal_mol / (R * temperature_K))
    kd_nM = kd_molar * 1e9  # 转换为 nM
    return kd_nM

# 示例
print(f"ΔG=-8  → Kd={delta_g_to_kd(-8):.1f} nM")   # ~1,100 nM
print(f"ΔG=-10 → Kd={delta_g_to_kd(-10):.2f} nM")  # ~8.5 nM
print(f"ΔG=-12 → Kd={delta_g_to_kd(-12):.3f} nM")  # ~0.065 nM
```

### 2.3 `ipsae` — 界面预测对齐误差

**来源：** Boltz-2 折叠后计算，反映界面区域的结构预测不确定性。

| iPSAE | 置信度 | 说明 |
|-------|--------|------|
| < 0.3 | 高置信 | 界面结构可靠，可信赖 |
| 0.3-0.5 | 中等置信 | 界面结构基本可靠，建议验证 |
| > 0.5 | 低置信 | 界面结构不可靠，不推荐 |

### 2.4 `iptm` — 界面 TM 分数

**来源：** Boltz-2 折叠后计算，反映复合物整体拓扑结构的可靠性。

| iPTM | 置信度 | 说明 |
|------|--------|------|
| > 0.75 | 高置信 | 复合物拓扑结构可靠 |
| 0.5-0.75 | 中等置信 | 拓扑结构基本正确 |
| < 0.5 | 低置信 | 拓扑结构不可靠 |

### 2.5 `immunogenicity` — 免疫原性风险

**计算方法：** 基于 NetMHCIIpan 预测 MHC-II 呈递概率，综合考量 HLA-DR/DP/DQ 等位基因覆盖。

| 免疫原性评分 | 风险等级 | 说明 | 建议 |
|------------|---------|------|------|
| < 0.2 | 极低风险 | 几乎无免疫原性 | 治疗性肽段首选 |
| 0.2-0.3 | 低风险 | 轻微免疫原性 | 可接受 |
| 0.3-0.5 | 中等风险 | 中等免疫原性 | 需要去免疫原性优化 |
| 0.5-0.7 | 高风险 | 高免疫原性 | 不推荐治疗用途 |
| > 0.7 | 极高风险 | 强免疫原性 | 禁止治疗用途 |

### 2.6 `solubility` — 溶解度（GRAVY 指数）

**计算方法：** GRAVY（Grand Average of Hydropathicity）指数，基于 Kyte-Doolittle 疏水性量表。

```
GRAVY = Σ(疏水性值) / 肽段长度
```

| GRAVY 值 | 溶解度 | 说明 |
|---------|--------|------|
| < -1.5 | 极高溶解度 | 高度亲水，水溶液中稳定 |
| -1.5 ~ -0.5 | 高溶解度 | 亲水，适合静脉注射 |
| -0.5 ~ 0 | 中等溶解度 | 可接受 |
| 0 ~ 0.5 | 低溶解度 | 疏水，可能聚集 |
| > 0.5 | 极低溶解度 | 高度疏水，不推荐 |

### 2.7 `stability` / `half_life` — 稳定性（Pro+）

```python
# Pro+ Tier 专属指标
if account.tier in ["pro", "enterprise"]:
    scores = client.scoring.deltaforge(
        peptides=["ACDEFGHIKLMN"],
        gene="EGFR",
        include_stability=True,    # 热稳定性
        include_half_life=True,    # 血浆半衰期
    )
    print(f"Stability: {scores[0].stability:.2f}")    # 0-1，>0.6 稳定
    print(f"Half-life: {scores[0].half_life_h:.1f}h") # 血浆半衰期（小时）
```

### 2.8 `cross_species` — 跨物种脱靶（Pro+）

```python
# 跨物种脱靶风险评估（Pro+）
if account.tier in ["pro", "enterprise"]:
    scores = client.scoring.deltaforge(
        peptides=["ACDEFGHIKLMN"],
        gene="EGFR",
        cross_species=["mouse", "rat", "monkey"],  # 评估物种
        include_cross_species=True,
    )
    print(f"Mouse off-target: {scores[0].cross_species['mouse']:.2f}")
    print(f"Monkey off-target: {scores[0].cross_species['monkey']:.2f}")
```

---

## 3. LigandIQ 预折叠筛选

LigandIQ 是一个轻量级的预折叠筛选器，**无需折叠即可快速评估**肽段的多项理化性质，
用于在昂贵的折叠步骤之前淘汰明显不合格的候选。

**LigandIQ 评估的指标（全部免费）：**
- 免疫原性（MHC-II 呈递预测）
- 溶解度（GRAVY 指数）
- 聚集倾向（TANGO/AGGRESCAN）
- 蛋白酶稳定性（蛋白酶切割位点预测）
- 序列复杂度（低复杂度区域检测）
- 初步结合亲和力估计（基于序列特征）

```python
# LigandIQ 快速预筛
ligandiq_scores = client.scoring.ligandiq(
    peptides=candidate_sequences,   # 所有候选肽段
    gene="EGFR",
    filters={
        "max_immunogenicity": 0.3,  # 免疫原性上限
        "min_solubility": -1.0,     # 溶解度下限（GRAVY）
        "max_aggregation": 0.5,     # 聚集倾向上限
    },
)

# 获取通过预筛的候选
passed = [s for s in ligandiq_scores if s.passed_filters]
print(f"Pre-screening: {len(passed)}/{len(candidate_sequences)} passed")
```

---

## 4. 完整筛选流程代码

```python
"""
完整两阶段筛选流程：
生成 → LigandIQ 预筛 → 折叠 Top 候选 → DeltaForge 精评 → 输出排名

Two-stage screening pipeline:
Generate → LigandIQ pre-screen → Fold top candidates → DeltaForge precision → Ranked output
"""

import ligandai
import pandas as pd
import numpy as np

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")
TARGET_GENE = "EGFR"

# ── Stage 0: 生成候选肽段 / Generate candidates ───────────────────────────────
print("="*60)
print("STAGE 0: Generating 300 Peptide Candidates")
print("="*60)

gen_job = client.peptides.generate(
    gene=TARGET_GENE,
    num_peptides=300,
    targeting_strategy="pocket_targeted",
    sampling_steps=50,   # 快速生成，后续精评
    auto_fold=False,     # 不自动折叠，节省 credits
)
gen_results = gen_job.wait()
all_sequences = [p.sequence for p in gen_results.peptides]
print(f"✓ Generated {len(all_sequences)} candidates")

# ── Stage 1: LigandIQ 预折叠筛选（免费）/ LigandIQ pre-screening (free) ───────
print("\nSTAGE 1: LigandIQ Pre-fold Screening (FREE)")
print("="*60)

ligandiq_scores = client.scoring.ligandiq(
    peptides=all_sequences,
    gene=TARGET_GENE,
    filters={
        "max_immunogenicity": 0.30,   # 免疫原性 < 0.3
        "min_solubility": -0.5,       # GRAVY > -0.5（亲水）
        "max_aggregation": 0.40,      # 聚集倾向 < 0.4
        "min_complexity": 0.6,        # 序列复杂度 > 0.6
    },
)

# 筛选通过的候选
passed_iq = [s for s in ligandiq_scores if s.passed_filters]
passed_sequences = [s.sequence for s in passed_iq]
print(f"✓ LigandIQ passed: {len(passed_sequences)}/{len(all_sequences)}")
print(f"  Filtered out: {len(all_sequences) - len(passed_sequences)} candidates")
print(f"  Credits saved: ~{(len(all_sequences) - len(passed_sequences)) * 100:,} credits")

# 按 LigandIQ 综合评分排序，取 Top 50 进入折叠
passed_iq_sorted = sorted(passed_iq, key=lambda x: x.composite_score, reverse=True)
top50_sequences = [s.sequence for s in passed_iq_sorted[:50]]
print(f"  Top 50 selected for folding")

# ── Stage 2: 折叠 Top 50 候选 / Fold top 50 candidates ───────────────────────
print("\nSTAGE 2: Folding Top 50 Candidates")
print("="*60)

fold_results = {}
for i, seq in enumerate(top50_sequences):
    print(f"  [{i+1:2d}/50] Folding: {seq[:12]}...", end=" ")
    
    fold = client.structures.fold_complex(
        receptor_gene=TARGET_GENE,
        peptide_sequence=seq,
        sampling_steps=100,
        num_trajectories=2,
        use_msa_cache=True,
    )
    fold.wait(timeout=600)
    
    fold_results[seq] = {
        "ipsae": fold.ipsae,
        "iptm": fold.iptm,
        "pdb_url": fold.pdb_url,
    }
    print(f"iPSAE={fold.ipsae:.3f}, iPTM={fold.iptm:.3f}")

# 筛选高置信折叠结果（iPSAE < 0.4）
high_conf_seqs = [seq for seq, r in fold_results.items() if r["ipsae"] < 0.4]
print(f"\n✓ High-confidence folds (iPSAE < 0.4): {len(high_conf_seqs)}/50")

# ── Stage 3: DeltaForge 精评（免费）/ DeltaForge precision scoring (free) ─────
print("\nSTAGE 3: DeltaForge Precision Scoring (FREE)")
print("="*60)

deltaforge_scores = client.scoring.deltaforge(
    peptides=high_conf_seqs,
    gene=TARGET_GENE,
    include_stability=True,    # Pro+ 功能
    include_half_life=True,    # Pro+ 功能
)

print(f"✓ DeltaForge scored {len(deltaforge_scores)} candidates")

# ── Stage 4: 综合排名 / Composite ranking ────────────────────────────────────
print("\nSTAGE 4: Composite Ranking")
print("="*60)

# 构建综合数据框
records = []
for score in deltaforge_scores:
    seq = score.sequence
    fold = fold_results.get(seq, {})
    iq = next((s for s in passed_iq if s.sequence == seq), None)
    
    records.append({
        "sequence":       seq,
        "length":         len(seq),
        "delta_g":        score.delta_g,
        "kd_nM":          score.kd,
        "ipsae":          fold.get("ipsae"),
        "iptm":           fold.get("iptm"),
        "immunogenicity": score.immunogenicity,
        "solubility":     score.solubility,
        "stability":      getattr(score, "stability", None),
        "half_life_h":    getattr(score, "half_life_h", None),
        "pdb_url":        fold.get("pdb_url"),
    })

df = pd.DataFrame(records)

# 多维度综合评分（见第5节）
df = compute_composite_score(df)  # 见下方函数定义
df_ranked = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
df_ranked["final_rank"] = df_ranked.index + 1

# 输出 Top 10
print("\nTOP 10 CANDIDATES:")
print("-"*80)
display_cols = ["final_rank", "sequence", "delta_g", "kd_nM",
                "ipsae", "immunogenicity", "solubility", "composite_score"]
print(df_ranked[display_cols].head(10).to_string(index=False))

# 保存结果
df_ranked.to_csv(f"{TARGET_GENE}_final_ranked.csv", index=False)
print(f"\n✓ Saved to {TARGET_GENE}_final_ranked.csv")
```

---

## 5. 多维度综合评分示例

```python
import numpy as np
import pandas as pd

def compute_composite_score(
    df: pd.DataFrame,
    weights: dict = None,
) -> pd.DataFrame:
    """
    多维度综合评分函数
    Multi-dimensional composite scoring function
    
    各指标归一化后加权求和，生成 0-1 的综合评分。
    Each metric is normalized then weighted-summed to produce a 0-1 composite score.
    """
    if weights is None:
        # 默认权重配置（可根据应用场景调整）
        # Default weights (adjust based on application)
        weights = {
            "delta_g":        0.35,   # 结合亲和力（最重要）
            "ipsae":          0.25,   # 结构置信度
            "immunogenicity": 0.20,   # 免疫原性（治疗性肽段重要）
            "solubility":     0.10,   # 溶解度
            "stability":      0.10,   # 稳定性（Pro+）
        }
    
    df = df.copy()
    
    # ── 归一化各指标 / Normalize each metric ─────────────────────────────────
    
    # delta_g：越负越好，归一化到 [0, 1]
    # delta_g: more negative = better, normalize to [0, 1]
    if "delta_g" in df.columns and df["delta_g"].notna().any():
        dg_min, dg_max = df["delta_g"].min(), df["delta_g"].max()
        if dg_min != dg_max:
            df["score_delta_g"] = (dg_max - df["delta_g"]) / (dg_max - dg_min)
        else:
            df["score_delta_g"] = 1.0
    else:
        df["score_delta_g"] = 0.5
        weights["delta_g"] = 0.0
    
    # ipsae：越小越好，归一化到 [0, 1]
    # ipsae: smaller = better, normalize to [0, 1]
    if "ipsae" in df.columns and df["ipsae"].notna().any():
        df["score_ipsae"] = 1.0 - df["ipsae"].clip(0, 1)
    else:
        df["score_ipsae"] = 0.5
        weights["ipsae"] = 0.0
    
    # immunogenicity：越小越好，归一化到 [0, 1]
    # immunogenicity: smaller = better, normalize to [0, 1]
    if "immunogenicity" in df.columns and df["immunogenicity"].notna().any():
        df["score_immunogenicity"] = 1.0 - df["immunogenicity"].clip(0, 1)
    else:
        df["score_immunogenicity"] = 0.5
        weights["immunogenicity"] = 0.0
    
    # solubility（GRAVY）：越负越好（亲水），归一化到 [0, 1]
    # solubility (GRAVY): more negative = better (hydrophilic), normalize to [0, 1]
    if "solubility" in df.columns and df["solubility"].notna().any():
        sol_min, sol_max = df["solubility"].min(), df["solubility"].max()
        if sol_min != sol_max:
            df["score_solubility"] = (sol_max - df["solubility"]) / (sol_max - sol_min)
        else:
            df["score_solubility"] = 1.0
    else:
        df["score_solubility"] = 0.5
        weights["solubility"] = 0.0
    
    # stability：越大越好，归一化到 [0, 1]
    # stability: larger = better, normalize to [0, 1]
    if "stability" in df.columns and df["stability"].notna().any():
        df["score_stability"] = df["stability"].clip(0, 1)
    else:
        df["score_stability"] = 0.5
        weights["stability"] = 0.0
    
    # ── 权重归一化 / Normalize weights ───────────────────────────────────────
    total_weight = sum(weights.values())
    norm_weights = {k: v / total_weight for k, v in weights.items()}
    
    # ── 计算综合评分 / Compute composite score ───────────────────────────────
    df["composite_score"] = (
        df["score_delta_g"]        * norm_weights.get("delta_g", 0) +
        df["score_ipsae"]          * norm_weights.get("ipsae", 0) +
        df["score_immunogenicity"] * norm_weights.get("immunogenicity", 0) +
        df["score_solubility"]     * norm_weights.get("solubility", 0) +
        df["score_stability"]      * norm_weights.get("stability", 0)
    )
    
    return df


# ── 使用示例 / Usage example ──────────────────────────────────────────────────

# 标准配置（平衡亲和力与安全性）
df_standard = compute_composite_score(df)

# 治疗性肽段配置（强调低免疫原性）
df_therapeutic = compute_composite_score(df, weights={
    "delta_g":        0.30,
    "ipsae":          0.20,
    "immunogenicity": 0.35,   # 提高免疫原性权重
    "solubility":     0.10,
    "stability":      0.05,
})

# 递送肽段配置（强调溶解度和稳定性）
df_delivery = compute_composite_score(df, weights={
    "delta_g":        0.25,
    "ipsae":          0.20,
    "immunogenicity": 0.15,
    "solubility":     0.25,   # 提高溶解度权重
    "stability":      0.15,   # 提高稳定性权重
})

print("Standard ranking top 3:")
print(df_standard.nlargest(3, "composite_score")[["sequence", "composite_score"]])

print("\nTherapeutic ranking top 3:")
print(df_therapeutic.nlargest(3, "composite_score")[["sequence", "composite_score"]])
```

---

## 6. Credit 成本

| 操作 | Credit 成本 | 说明 |
|------|------------|------|
| LigandIQ 预筛 | **0 credits（免费）** | CPU-only，无限次使用 |
| DeltaForge 精评 | **0 credits（免费）** | CPU-only，无限次使用 |
| 折叠（50步，1轨迹） | 100 credits/peptide | 唯一需要付费的步骤 |

**两阶段筛选节省的 Credits 估算：**
```
场景：300 条候选肽段，最终折叠 Top 50

无预筛方案：300 × 100 = 30,000 credits
有预筛方案：50 × 100 = 5,000 credits
节省：25,000 credits（83% 节省）
```

---

## 7. 候选肽段筛选决策树

```
所有候选肽段（N 条）
│
├─ [LigandIQ 预筛，免费]
│   ├── 免疫原性 > 0.3？ → 淘汰
│   ├── GRAVY > 0？ → 淘汰（疏水）
│   ├── 聚集倾向 > 0.4？ → 淘汰
│   └── 通过 → 进入下一步（预计保留 40-60%）
│
├─ [按 LigandIQ 综合评分排序，取 Top 50]
│
├─ [Boltz-2 折叠，100 credits/条]
│   ├── iPSAE > 0.5？ → 淘汰（低置信）
│   ├── iPTM < 0.5？ → 淘汰（低置信）
│   └── 通过 → 进入精评（预计保留 60-80%）
│
├─ [DeltaForge 精评，免费]
│   ├── ΔG > -6 kcal/mol？ → 淘汰（弱结合）
│   ├── Kd > 1000 nM？ → 淘汰
│   └── 通过 → 候选列表
│
├─ [多维度综合评分排名]
│   └── 输出 Top 10 候选
│
└─ [实验验证]
    ├── SPR / ITC（结合亲和力）
    ├── 细胞毒性测试
    └── 体内药代动力学
```

> **技能版本：** v1.0.0 | **DeltaForge 版本：** v3.2 | **LigandIQ 版本：** v2.0 | **SDK 最低版本：** ligandai>=2.1.0
