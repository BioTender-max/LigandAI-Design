---
name: ligandai-bivalent
description: |
  【中文】基于 LigandAI 平台的双价/双特异性肽段设计技能（Beta，Pro+ 专属）。
  通过诱导邻近（Induced Proximity）机制，设计同时结合两个蛋白靶点的双功能肽段，
  支持 PROTAC 类降解剂、分子胶、双特异性抗体替代肽等应用场景。
  
  [EN] Bivalent/bispecific peptide design skill powered by LigandAI platform (Beta, Pro+ exclusive).
  Designs bifunctional peptides that simultaneously engage two protein targets via induced proximity
  mechanism. Supports PROTAC-like degraders, molecular glues, and bispecific antibody-mimetic peptides.
license: Apache-2.0
category: design-tools
tags:
  - bivalent
  - bispecific
  - induced-proximity
  - protac
  - degrader
  - pro-tier

triggers:
  # 中文触发场景
  - "设计双价肽段同时结合两个靶点"
  - "PROTAC 类肽段降解剂设计"
  - "诱导邻近机制肽段设计"
  - "双特异性肽段 BRD4 CRBN"
  - "分子胶肽段设计"
  - "设计同时靶向 E3 连接酶和底物的肽段"
  - "双功能肽段 linker 设计"
  - "LigandAI bivalent 双价设计"
  # English triggers
  - "design bivalent peptide targeting two proteins"
  - "PROTAC-like peptide degrader design"
  - "induced proximity peptide design"
  - "bispecific peptide BRD4 CRBN"
  - "molecular glue peptide design"
  - "design bifunctional peptide with linker"
  - "LigandAI bivalent design"

inputs:
  - name: target_a
    type: string
    required: true
    description: "第一个靶点基因符号（底物蛋白）/ First target gene symbol (substrate protein)"
  - name: target_b
    type: string
    required: true
    description: "第二个靶点基因符号（效应蛋白，如 E3 连接酶）/ Second target gene symbol (effector, e.g. E3 ligase)"
  - name: linker_type
    type: string
    required: false
    default: "GGS"
    description: "连接子类型：GGS / PEG / rigid_helix / custom"
  - name: linker_length
    type: integer
    required: false
    default: 3
    description: "连接子重复单元数（GGS 重复数或 PEG 单元数）/ Linker repeat units"
  - name: proximity_distance
    type: float
    required: false
    default: 15.0
    description: "两靶点期望邻近距离（Å）/ Desired proximity distance between targets (Å)"
  - name: num_designs
    type: integer
    required: false
    default: 50
    description: "生成双价设计数量 / Number of bivalent designs to generate"

outputs:
  - name: designs
    type: list[BivalentDesign]
    description: "双价肽段设计列表 / Bivalent peptide design list"
  - name: designs[].full_sequence
    type: string
    description: "完整双价肽段序列（臂A + 连接子 + 臂B）/ Full bivalent sequence (arm_a + linker + arm_b)"
  - name: designs[].arm_a_sequence
    type: string
    description: "靶向 target_a 的肽段臂序列 / Peptide arm targeting target_a"
  - name: designs[].arm_b_sequence
    type: string
    description: "靶向 target_b 的肽段臂序列 / Peptide arm targeting target_b"
  - name: designs[].linker_sequence
    type: string
    description: "连接子序列 / Linker sequence"
  - name: designs[].proximity_score
    type: float
    description: "诱导邻近评分（0-1，越高越好）/ Induced proximity score (0-1)"
  - name: designs[].delta_g_a
    type: float
    description: "臂A与 target_a 的结合自由能 / Arm A binding free energy to target_a"
  - name: designs[].delta_g_b
    type: float
    description: "臂B与 target_b 的结合自由能 / Arm B binding free energy to target_b"
---

# LigandAI Bivalent — 双价/双特异性肽段设计技能

> ⚠️ **Beta 功能** | 需要 **Pro+ Tier** | 当前版本：v0.9-beta

## 1. 功能概述

### 中文说明

LigandAI Bivalent 模块实现了基于**诱导邻近（Induced Proximity）**机制的双功能肽段设计。
与传统单靶点肽段不同，双价肽段通过柔性或刚性连接子将两个功能性肽段臂连接，
使两个原本独立的蛋白质在细胞内被强制拉近，从而触发特定的生物学效应。

### 诱导邻近机制示意

```
┌─────────────────────────────────────────────────────────────────┐
│                    诱导邻近机制示意图                             │
│                 Induced Proximity Mechanism                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   Target A (底物蛋白)          Target B (效应蛋白/E3连接酶)      │
│   ┌──────────┐                    ┌──────────┐                  │
│   │  BRD4    │                    │   CRBN   │                  │
│   │  (BD1)   │                    │ (E3 lig) │                  │
│   └────┬─────┘                    └─────┬────┘                  │
│        │ 结合                           │ 结合                   │
│        │ Binding                        │ Binding               │
│   ┌────┴─────────────────────────────────┴────┐                 │
│   │  Arm A  ──── Linker ────  Arm B           │                 │
│   │  (8-15aa)  (GGS×3)    (8-15aa)           │                 │
│   └───────────────────────────────────────────┘                 │
│                    双价肽段                                      │
│                  Bivalent Peptide                               │
│                                                                 │
│   效果：BRD4 被 CRBN 泛素化 → 蛋白酶体降解                      │
│   Effect: BRD4 ubiquitinated by CRBN → Proteasomal degradation │
│                                                                 │
│   邻近距离 / Proximity distance: ~15 Å                          │
│   连接子长度 / Linker length: GGS×3 = 9 aa ≈ 30 Å              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

**诱导邻近的三种主要应用：**

1. **靶向蛋白降解（PROTAC 类）：** 将底物蛋白（如 BRD4）与 E3 泛素连接酶（如 CRBN/VHL）拉近，
   触发底物泛素化和蛋白酶体降解。

2. **分子胶（Molecular Glue）：** 稳定两个蛋白质的相互作用界面，增强天然弱相互作用。

3. **双特异性抗体替代肽：** 同时结合两个细胞表面受体，介导细胞-细胞接触（如 T 细胞与肿瘤细胞）。

---

## 2. 与标准生成的区别

### 何时使用 Bivalent vs 标准流程

| 场景 | 推荐方案 | 原因 |
|------|---------|------|
| 单靶点结合/抑制 | 标准 `ligandai-peptide` | 简单高效，成本低 |
| 靶向蛋白降解（PROTAC 类） | **Bivalent** | 需要同时结合底物和 E3 连接酶 |
| 分子胶设计 | **Bivalent** | 需要稳定两蛋白界面 |
| 双特异性细胞接触 | **Bivalent** | 需要同时结合两个细胞表面受体 |
| 变构调节（同一蛋白两个位点） | 标准（指定两个口袋） | 同一蛋白，用标准流程 |
| 同源二聚体/多聚体 | 标准（指定链） | 原生多聚体，用标准流程 |
| 受体-配体复合物 | 标准 | 天然相互作用，用标准流程 |

> **重要：** Bivalent 模块**仅适用于两个独立蛋白**之间的诱导邻近设计。
> 如果目标是同一蛋白的两个结合位点，或天然存在的蛋白复合物，请使用标准 `ligandai-peptide` 流程。

---

## 3. `client.bivalent.start()` 参数详解

```python
job = client.bivalent.start(
    # ── 必填参数 / Required ──────────────────────────────────────────────────
    target_a="BRD4",                  # 第一靶点（底物蛋白）
    target_b="CRBN",                  # 第二靶点（效应蛋白/E3连接酶）

    # ── 连接子设计 / Linker design ───────────────────────────────────────────
    linker_type="GGS",                # 连接子类型：GGS / PEG / rigid_helix / custom
    linker_length=3,                  # GGS 重复数（GGS×3 = GGSGGSGSS，9 aa）
    custom_linker=None,               # 自定义连接子序列（linker_type="custom" 时）

    # ── 邻近约束 / Proximity constraints ────────────────────────────────────
    proximity_distance=15.0,          # 期望邻近距离（Å），推荐 10-25 Å
    proximity_mode="flexible",        # flexible / rigid / optimized

    # ── 臂设计参数 / Arm design parameters ──────────────────────────────────
    arm_a_length=(8, 15),             # 臂A长度范围（氨基酸数）
    arm_b_length=(8, 15),             # 臂B长度范围（氨基酸数）
    arm_a_pocket="P1",                # 臂A靶向口袋（可选，来自 structures.analyze()）
    arm_b_pocket="P1",                # 臂B靶向口袋（可选）
    arm_a_residues=None,              # 臂A热点残基（可选）
    arm_b_residues=None,              # 臂B热点残基（可选）

    # ── 生成参数 / Generation parameters ────────────────────────────────────
    num_designs=50,                   # 生成双价设计数量
    sampling_steps=100,               # 扩散采样步数
    auto_fold=True,                   # 自动折叠三元复合物
    top_n_fold=5,                     # 折叠前 5 个设计

    # ── 约束条件 / Constraints ───────────────────────────────────────────────
    max_total_length=50,              # 最大总长度（氨基酸数）
    deimmunize=False,                 # 去免疫原性优化
    avoid_aggregation=True,           # 避免聚集
)
```

**完整参数参考表：**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `target_a` | `str` | 必填 | 第一靶点基因符号（底物蛋白） |
| `target_b` | `str` | 必填 | 第二靶点基因符号（效应蛋白） |
| `linker_type` | `str` | `"GGS"` | 连接子类型 |
| `linker_length` | `int` | `3` | 连接子重复单元数 |
| `custom_linker` | `str` | `None` | 自定义连接子序列 |
| `proximity_distance` | `float` | `15.0` | 期望邻近距离（Å） |
| `proximity_mode` | `str` | `"flexible"` | 邻近约束模式 |
| `arm_a_length` | `tuple` | `(8, 15)` | 臂A长度范围 |
| `arm_b_length` | `tuple` | `(8, 15)` | 臂B长度范围 |
| `arm_a_pocket` | `str` | `None` | 臂A靶向口袋 ID |
| `arm_b_pocket` | `str` | `None` | 臂B靶向口袋 ID |
| `arm_a_residues` | `list[int]` | `None` | 臂A热点残基 |
| `arm_b_residues` | `list[int]` | `None` | 臂B热点残基 |
| `num_designs` | `int` | `50` | 生成设计数量 |
| `sampling_steps` | `int` | `100` | 扩散采样步数 |
| `auto_fold` | `bool` | `True` | 自动折叠三元复合物 |
| `top_n_fold` | `int` | `5` | 自动折叠数量 |
| `max_total_length` | `int` | `50` | 最大总长度 |
| `deimmunize` | `bool` | `False` | 去免疫原性优化 |
| `avoid_aggregation` | `bool` | `True` | 避免聚集序列 |

---

## 4. Linker 设计指南

### 4.1 连接子类型比较

| 连接子类型 | 序列示例 | 特性 | 适用场景 |
|-----------|---------|------|---------|
| `GGS` | `(GGSGGS)×n` | 柔性，亲水，低免疫原性 | 大多数 PROTAC 类应用（推荐） |
| `PEG` | `-(CH₂CH₂O)ₙ-` | 极柔性，无免疫原性，增加溶解度 | 需要最大柔性，溶解度优先 |
| `rigid_helix` | `(EAAAK)×n` | 刚性 α-螺旋，精确控制距离 | 需要精确空间定位 |
| `custom` | 用户定义 | 完全自定义 | 特殊需求 |

### 4.2 连接子长度选择建议

**GGS 连接子长度与空间距离的关系：**

| GGS 重复数 | 氨基酸数 | 估计伸展长度 | 适用邻近距离 |
|-----------|---------|------------|------------|
| 1 | 3 aa | ~10 Å | 5-15 Å |
| 2 | 6 aa | ~20 Å | 10-25 Å |
| 3 | 9 aa | ~30 Å | 15-35 Å |
| 4 | 12 aa | ~40 Å | 25-45 Å |
| 5 | 15 aa | ~50 Å | 35-55 Å |

> **经验法则：** GGS 连接子每个重复单元（3 aa）提供约 10 Å 的空间距离。
> 建议连接子长度略大于期望邻近距离（留有柔性余量）。

### 4.3 连接子选择决策

```python
def recommend_linker(proximity_distance_A: float, application: str) -> dict:
    """
    根据邻近距离和应用场景推荐连接子
    Recommend linker based on proximity distance and application
    """
    if application == "protac":
        # PROTAC 类：优先 GGS，柔性好
        n_repeats = max(1, round(proximity_distance_A / 10))
        return {"type": "GGS", "length": n_repeats,
                "sequence": "GGS" * n_repeats}
    
    elif application == "molecular_glue":
        # 分子胶：需要精确定位，考虑刚性连接子
        if proximity_distance_A < 20:
            return {"type": "rigid_helix", "length": 1,
                    "sequence": "EAAAK"}
        else:
            n_repeats = max(1, round(proximity_distance_A / 10))
            return {"type": "GGS", "length": n_repeats,
                    "sequence": "GGS" * n_repeats}
    
    elif application == "bispecific":
        # 双特异性：细胞表面受体，需要较长连接子
        n_repeats = max(3, round(proximity_distance_A / 8))
        return {"type": "GGS", "length": n_repeats,
                "sequence": "GGS" * n_repeats}
    
    else:
        # 默认：GGS×3
        return {"type": "GGS", "length": 3, "sequence": "GGSGGSGGS"}

# 示例
linker = recommend_linker(proximity_distance_A=15.0, application="protac")
print(f"Recommended: {linker['type']}×{linker['length']} = {linker['sequence']}")
```

---

## 5. 完整代码示例：BRD4-CRBN PROTAC 类设计

```python
"""
完整示例：BRD4-CRBN PROTAC 类双价肽段设计
Full example: BRD4-CRBN PROTAC-like bivalent peptide design

目标：设计同时结合 BRD4(BD1) 和 CRBN 的双价肽段，
      诱导 BRD4 泛素化降解，用于癌症治疗。

Goal: Design bivalent peptide simultaneously binding BRD4(BD1) and CRBN,
      inducing BRD4 ubiquitination and degradation for cancer therapy.
"""

import ligandai
from ligandai.exceptions import LigandAITierError, LigandAIAPIError
import pandas as pd

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# ── Step 0: Tier 检查 / Tier check ───────────────────────────────────────────
print("="*60)
print("STEP 0: Checking Tier Requirements")
print("="*60)

account = client.account.info()
print(f"Current tier: {account.tier}")
print(f"Credits: {account.credits:,}")

if account.tier not in ["pro", "enterprise"]:
    print(f"❌ Bivalent design requires Pro+ tier (current: {account.tier})")
    print("   Please upgrade at: https://ligandai.com/pricing")
    raise SystemExit(1)

print(f"✓ Tier check passed: {account.tier}")

# ── Step 1: 解析两个靶点结构 / Resolve both target structures ─────────────────
print("\nSTEP 1: Resolving Target Structures")
print("="*60)

# BRD4 BD1 结构域
brd4_structure = client.structures.get(gene="BRD4", domain="BD1")
brd4_pocket = client.structures.analyze(gene="BRD4", analysis_depth="standard")
brd4_best_pocket = brd4_pocket.recommended_pocket
print(f"✓ BRD4 BD1: {brd4_structure.source} | Pocket: {brd4_best_pocket.pocket_id}")
print(f"  Key residues: {brd4_best_pocket.key_residues[:5]}")

# CRBN（Cereblon，E3 泛素连接酶）
crbn_structure = client.structures.get(gene="CRBN")
crbn_pocket = client.structures.analyze(gene="CRBN", analysis_depth="standard")
crbn_best_pocket = crbn_pocket.recommended_pocket
print(f"✓ CRBN: {crbn_structure.source} | Pocket: {crbn_best_pocket.pocket_id}")
print(f"  Key residues: {crbn_best_pocket.key_residues[:5]}")

# ── Step 2: 启动双价设计 / Start bivalent design ──────────────────────────────
print("\nSTEP 2: Starting Bivalent Design (BRD4-CRBN)")
print("="*60)

try:
    job = client.bivalent.start(
        # 靶点定义
        target_a="BRD4",
        target_b="CRBN",
        
        # 口袋信息（来自 Step 1）
        arm_a_pocket=brd4_best_pocket.pocket_id,
        arm_b_pocket=crbn_best_pocket.pocket_id,
        arm_a_residues=brd4_best_pocket.key_residues,
        arm_b_residues=crbn_best_pocket.key_residues,
        
        # 连接子设计（GGS×3，适合 ~15 Å 邻近距离）
        linker_type="GGS",
        linker_length=3,              # GGSGGSGSS = 9 aa ≈ 30 Å
        
        # 邻近约束
        proximity_distance=15.0,      # BRD4-CRBN 期望邻近距离
        proximity_mode="flexible",
        
        # 臂长度
        arm_a_length=(8, 14),         # BRD4 结合臂
        arm_b_length=(8, 14),         # CRBN 结合臂
        
        # 生成参数
        num_designs=50,
        sampling_steps=100,
        auto_fold=True,
        top_n_fold=5,
        
        # 约束
        max_total_length=45,          # 总长度 ≤ 45 aa
        deimmunize=True,              # 治疗性应用，去免疫原性
        avoid_aggregation=True,
    )
    
    print(f"✓ Job submitted: {job.job_id}")
    print(f"  Estimated time: {job.estimated_time}s")
    print(f"  Estimated credits: {job.estimated_credits:,}")

except LigandAITierError as e:
    # Tier 不足错误处理（见 Step 6）
    print(f"❌ Tier error: {e.message}")
    print(f"   Required: {e.required_tier}")
    print(f"   Current: {e.current_tier}")
    raise

# ── Step 3: 流式监听进度 / Stream progress ───────────────────────────────────
print("\nSTEP 3: Streaming Progress")
print("="*60)

for event in job.stream():
    if event.type == "arm_generation_complete":
        print(f"  ✓ Arm generation: {event.arm} | {event.count} sequences")
    elif event.type == "linker_optimization":
        print(f"  🔗 Linker optimization: {event.iteration}/{event.total}")
    elif event.type == "proximity_check":
        print(f"  📏 Proximity check: {event.passed}/{event.total} passed "
              f"(target: {event.target_distance:.1f}Å)")
    elif event.type == "elite_hit":
        print(f"  ⭐ Elite design: {event.full_sequence[:20]}... | "
              f"Proximity={event.proximity_score:.3f}")
    elif event.type == "fold_complete":
        print(f"  ✅ Ternary fold: iPSAE={event.ipsae:.3f}, iPTM={event.iptm:.3f}")
    elif event.type == "generation_complete":
        print(f"\n  🎉 Complete! {event.total_designs} designs generated")

# ── Step 4: 获取并分析结果 / Get and analyze results ─────────────────────────
print("\nSTEP 4: Analyzing Results")
print("="*60)

results = job.results()

records = []
for i, design in enumerate(results.designs):
    records.append({
        "rank":             i + 1,
        "full_sequence":    design.full_sequence,
        "arm_a":            design.arm_a_sequence,
        "linker":           design.linker_sequence,
        "arm_b":            design.arm_b_sequence,
        "total_length":     len(design.full_sequence),
        "proximity_score":  round(design.proximity_score, 3),
        "delta_g_a":        round(design.delta_g_a, 2),
        "delta_g_b":        round(design.delta_g_b, 2),
        "delta_g_combined": round((design.delta_g_a + design.delta_g_b) / 2, 2),
        "ipsae_ternary":    round(design.ipsae_ternary, 3) if design.ipsae_ternary else None,
        "immunogenicity":   round(design.immunogenicity, 3),
        "solubility":       round(design.solubility, 3),
        "pdb_url":          design.pdb_url,
    })

df = pd.DataFrame(records)

# 综合排名：邻近评分 × 结合亲和力
df["composite"] = (
    df["proximity_score"] * 0.4 +
    (1 - df["delta_g_combined"].clip(-15, 0) / -15) * 0.4 +
    (1 - df["immunogenicity"]) * 0.2
)
df_ranked = df.sort_values("composite", ascending=False).reset_index(drop=True)

print("\nTOP 5 BRD4-CRBN BIVALENT DESIGNS:")
print("-"*80)
display_cols = ["rank", "arm_a", "linker", "arm_b",
                "proximity_score", "delta_g_a", "delta_g_b", "immunogenicity"]
print(df_ranked[display_cols].head(5).to_string(index=False))

# 保存结果
df_ranked.to_csv("BRD4_CRBN_bivalent_designs.csv", index=False)
print(f"\n✓ Saved to BRD4_CRBN_bivalent_designs.csv")

# 最优设计
best = df_ranked.iloc[0]
print(f"\n✓ Best design:")
print(f"  Full sequence: {best['full_sequence']}")
print(f"  Arm A (BRD4):  {best['arm_a']}")
print(f"  Linker:        {best['linker']}")
print(f"  Arm B (CRBN):  {best['arm_b']}")
print(f"  Proximity:     {best['proximity_score']:.3f}")
print(f"  ΔG(BRD4):      {best['delta_g_a']:.2f} kcal/mol")
print(f"  ΔG(CRBN):      {best['delta_g_b']:.2f} kcal/mol")
```

---

## 6. Tier 检查与错误处理

```python
"""
完整的 Tier 检查与错误处理代码
Complete tier check and error handling code
"""

import ligandai
from ligandai.exceptions import (
    LigandAITierError,      # Tier 不足
    LigandAIAPIError,       # API 错误
    LigandAIQuotaError,     # 配额超限
    LigandAIBetaError,      # Beta 功能限制
    LigandAIValidationError, # 参数验证错误
)

client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

def safe_bivalent_design(target_a: str, target_b: str, **kwargs):
    """
    带完整错误处理的双价设计函数
    Bivalent design function with comprehensive error handling
    """
    
    # ── 预检查 1：Tier 验证 ───────────────────────────────────────────────────
    account = client.account.info()
    
    if account.tier not in ["pro", "enterprise"]:
        print(f"❌ Bivalent design requires Pro+ tier")
        print(f"   Current tier: {account.tier}")
        print(f"   Upgrade at: https://ligandai.com/pricing")
        print(f"\n   Alternative: Use standard peptide generation for single targets:")
        print(f"   client.peptides.generate(gene='{target_a}', ...)")
        return None
    
    # ── 预检查 2：Credits 验证 ────────────────────────────────────────────────
    num_designs = kwargs.get("num_designs", 50)
    estimated_credits = num_designs * 150  # 双价设计约 150 credits/design
    
    if account.credits < estimated_credits:
        print(f"❌ Insufficient credits")
        print(f"   Required: ~{estimated_credits:,} credits")
        print(f"   Available: {account.credits:,} credits")
        print(f"   Suggestion: Reduce num_designs to {account.credits // 150}")
        return None
    
    # ── 预检查 3：Beta 功能确认 ───────────────────────────────────────────────
    beta_features = client.account.beta_features()
    if "bivalent" not in beta_features.enabled:
        print(f"❌ Bivalent feature not enabled for your account")
        print(f"   Request access at: https://ligandai.com/beta")
        return None
    
    # ── 执行双价设计 ──────────────────────────────────────────────────────────
    try:
        job = client.bivalent.start(
            target_a=target_a,
            target_b=target_b,
            **kwargs
        )
        
        print(f"✓ Job started: {job.job_id}")
        results = job.wait(timeout=3600)
        return results
    
    except LigandAITierError as e:
        print(f"❌ Tier error: {e.message}")
        print(f"   Required tier: {e.required_tier}")
        print(f"   Current tier: {e.current_tier}")
        print(f"   Feature: {e.feature}")
        return None
    
    except LigandAIBetaError as e:
        print(f"❌ Beta feature error: {e.message}")
        print(f"   Beta feature: {e.feature}")
        print(f"   Status: {e.beta_status}")  # e.g. "waitlist" / "disabled"
        return None
    
    except LigandAIQuotaError as e:
        print(f"❌ Quota exceeded: {e.message}")
        print(f"   Quota type: {e.quota_type}")  # e.g. "monthly_designs"
        print(f"   Reset date: {e.reset_date}")
        return None
    
    except LigandAIValidationError as e:
        print(f"❌ Validation error: {e.message}")
        print(f"   Invalid parameter: {e.parameter}")
        print(f"   Provided value: {e.value}")
        print(f"   Allowed values: {e.allowed_values}")
        return None
    
    except LigandAIAPIError as e:
        print(f"❌ API error: {e.message}")
        print(f"   Status code: {e.status_code}")
        print(f"   Request ID: {e.request_id}")
        if e.recoverable:
            print(f"   This error is recoverable, retrying...")
            import time
            time.sleep(30)
            return safe_bivalent_design(target_a, target_b, **kwargs)
        return None


# 使用示例
results = safe_bivalent_design(
    target_a="BRD4",
    target_b="CRBN",
    linker_type="GGS",
    linker_length=3,
    num_designs=50,
    auto_fold=True,
    top_n_fold=5,
)

if results:
    print(f"✓ Success! {len(results.designs)} designs generated")
    print(f"  Best design: {results.designs[0].full_sequence}")
```

---

## 7. 应用场景

### 7.1 PROTAC 类靶向蛋白降解

```python
# 经典 PROTAC 类设计：底物蛋白 + E3 连接酶
protac_targets = [
    ("BRD4", "CRBN"),    # BRD4 降解（癌症）
    ("BRD4", "VHL"),     # BRD4 降解（替代 E3）
    ("KRAS", "CRBN"),    # KRAS 降解（难成药靶点）
    ("AR",   "CRBN"),    # 雄激素受体降解（前列腺癌）
    ("BCL2", "VHL"),     # BCL2 降解（淋巴瘤）
    ("CDK6", "CRBN"),    # CDK6 降解（乳腺癌）
]

for target_a, target_b in protac_targets:
    print(f"PROTAC design: {target_a} → {target_b}")
```

### 7.2 分子胶（Molecular Glue）

```python
# 分子胶：稳定两蛋白相互作用
molecular_glue = client.bivalent.start(
    target_a="IKZF1",    # Ikaros（底物）
    target_b="CRBN",     # Cereblon（E3 连接酶）
    linker_type="rigid_helix",  # 刚性连接子，精确定位
    linker_length=1,
    proximity_distance=8.0,     # 更近的邻近距离
    proximity_mode="rigid",
    num_designs=30,
)
```

### 7.3 双特异性抗体替代肽（BiTE 类）

```python
# 双特异性：T 细胞接触肿瘤细胞
bite_like = client.bivalent.start(
    target_a="CD3E",     # T 细胞表面 CD3ε
    target_b="EGFR",     # 肿瘤细胞表面 EGFR
    linker_type="GGS",
    linker_length=5,     # 较长连接子，跨越细胞间距
    proximity_distance=30.0,
    arm_a_length=(10, 18),
    arm_b_length=(10, 18),
    num_designs=50,
    deimmunize=True,     # 治疗性应用
)
```

---

## 8. 注意事项

### ⚠️ 重要限制

1. **仅适用于两个独立蛋白：**
   Bivalent 模块设计用于诱导两个**原本独立**的蛋白质相互靠近。
   如果目标蛋白在生理条件下已经形成复合物（如同源二聚体、天然异源二聚体），
   请使用标准 `ligandai-peptide` 流程，指定多个口袋或链。

2. **原生多聚体使用标准流程：**
   ```python
   # ❌ 错误：用 bivalent 处理同源二聚体
   # Wrong: Using bivalent for homodimer
   job = client.bivalent.start(target_a="TNF", target_b="TNF")  # TNF 是三聚体！
   
   # ✅ 正确：用标准流程，指定多链
   # Correct: Use standard flow with multi-chain
   job = client.peptides.generate(
       gene="TNF",
       targeting_strategy="pocket_targeted",
       pocket_id="trimer_interface",
   )
   ```

3. **连接子长度与邻近距离匹配：**
   连接子过短会导致空间张力，过长会降低诱导邻近效率。
   建议连接子伸展长度比期望邻近距离大 20-50%（留有柔性余量）。

4. **Beta 功能稳定性：**
   Bivalent 模块目前处于 Beta 阶段，API 接口可能在未来版本中变更。
   建议在生产环境中固定 SDK 版本：`pip install ligandai==2.1.x`

5. **三元复合物折叠成本：**
   双价肽段的三元复合物折叠（肽段 + 靶点A + 靶点B）比二元复合物更耗时，
   Credit 成本约为标准折叠的 2-3 倍。

### 推荐工作流

```
ligandai-discovery (发现两个靶点)
    ↓
ligandai-folding (分别解析两个靶点结构)
    ↓
ligandai-bivalent (双价设计)
    ↓
ligandai-scoring (评分筛选)
    ↓
实验验证 (三元复合物 SPR / 细胞降解实验)
```

> **技能版本：** v0.9-beta | **最低 Tier：** Pro | **SDK 最低版本：** ligandai>=2.1.0
> 
> **Beta 反馈：** https://ligandai.com/beta-feedback | **文档：** https://docs.ligandai.com/bivalent
