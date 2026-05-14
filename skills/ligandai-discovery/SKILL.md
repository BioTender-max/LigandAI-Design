---
name: ligandai-discovery
description: |
  【中文】基于 LigandAI 平台的靶点发现技能。利用 GTEx 批量 RNA-seq 数据、单细胞 RNA-seq 图谱及专有 vasculome/BBB 数据库，
  通过 Specificity Index (SI) 算法识别组织/细胞类型特异性表达的受体靶点，为后续肽段设计提供高质量候选基因列表。
  
  [EN] Target discovery skill powered by the LigandAI platform. Leverages GTEx bulk RNA-seq, single-cell RNA-seq atlases,
  and proprietary vasculome/BBB databases to identify tissue- and cell-type-specific receptor targets via the
  Specificity Index (SI) algorithm, producing a ranked candidate gene list for downstream peptide design.
license: Apache-2.0
category: design-tools
tags:
  - discovery
  - tissue-markers
  - scrna
  - gtex
  - target-identification

triggers:
  # 中文触发场景
  - "帮我找肝脏特异性受体靶点"
  - "用 LigandAI 做靶点发现"
  - "哪些基因在脑血管内皮细胞高表达"
  - "筛选肿瘤微环境中的受体"
  - "GTEx 组织特异性基因分析"
  - "单细胞数据找细胞类型标志物"
  - "发现 BBB 转运靶点"
  - "找心肌细胞特异性受体用于靶向递送"
  # English triggers
  - "find tissue-specific receptor targets using LigandAI"
  - "identify liver-specific surface receptors"
  - "run target discovery for brain endothelial cells"
  - "which genes are specifically expressed in tumor microenvironment"
  - "GTEx specificity index analysis"
  - "scRNA-seq cell type marker discovery"
  - "find BBB transport targets"
  - "discover cardiac-specific receptors for targeted delivery"

inputs:
  - name: tissue
    type: string
    required: true
    description: "目标组织名称（GTEx 标准命名），如 'liver', 'brain', 'heart', 'lung' 等 / Target tissue name (GTEx standard), e.g. 'liver', 'brain'"
  - name: cell_type
    type: string
    required: false
    description: "目标细胞类型（scRNA 模式），如 'hepatocyte', 'endothelial', 'T cell' / Target cell type for scRNA mode"
  - name: receptor_only
    type: boolean
    required: false
    default: true
    description: "是否仅返回受体类基因（推荐 True 用于肽段设计）/ Whether to return receptor-class genes only"
  - name: si_threshold
    type: float
    required: false
    default: 2.0
    description: "Specificity Index 最低阈值，建议 ≥2.0 / Minimum SI threshold, recommended ≥2.0"
  - name: top_n
    type: integer
    required: false
    default: 20
    description: "返回前 N 个候选靶点 / Number of top candidates to return"
  - name: geo_accession
    type: string
    required: false
    description: "GEO 数据集编号（可选，用于自定义数据集）/ GEO accession for custom dataset import"

outputs:
  - name: targets
    type: list[TargetResult]
    description: "按 SI 降序排列的靶点列表 / Target list sorted by SI descending"
  - name: targets[].gene
    type: string
    description: "基因符号 / Gene symbol (HGNC)"
  - name: targets[].si
    type: float
    description: "Specificity Index 值 / Specificity Index value"
  - name: targets[].tissue_specificity
    type: string
    description: "组织特异性等级：high / medium / low"
  - name: targets[].receptor_class
    type: string
    description: "受体分类：GPCR / RTK / ion_channel / transporter / other"
  - name: targets[].mean_tpm
    type: float
    description: "目标组织平均 TPM 表达量 / Mean TPM in target tissue"
  - name: targets[].rank
    type: integer
    description: "SI 排名 / Rank by SI"
  - name: targets[].uniprot_id
    type: string
    description: "UniProt 蛋白 ID / UniProt protein ID"
  - name: targets[].pdb_available
    type: boolean
    description: "是否有 PDB 结构 / Whether PDB structure is available"
---

# LigandAI Discovery — 靶点发现技能

## 1. 功能概述

### 中文说明

LigandAI Discovery 模块通过整合多层次转录组数据，系统性地识别在特定组织或细胞类型中高度特异性表达的受体蛋白，
为靶向肽段设计提供科学依据。其核心算法 **Specificity Index (SI)** 综合考量目标组织的表达量与全身其他组织的背景表达，
有效避免脱靶风险。

**主要数据源：**
- **GTEx v8**：54 种人体组织的批量 RNA-seq 数据（Free 及以上 Tier 可用）
- **单细胞 RNA-seq 图谱**：Human Cell Atlas、Tabula Sapiens 等（Academia+ Tier）
- **Vasculome 数据库**：血管内皮细胞亚型特异性表达谱（Enterprise Tier）
- **BBB 转运组**：血脑屏障转运蛋白专项数据库（Enterprise Tier）

### What It Does

| Feature | Description | Data Source | Min Tier |
|---------|-------------|-------------|----------|
| Bulk tissue markers | GTEx 54-tissue specificity ranking | GTEx v8 | Free |
| Cell-type markers | scRNA-seq cell-type resolved expression | HCA / Tabula Sapiens | Academia |
| Receptor filtering | Filter to surface-accessible receptor classes | Uniprot + GO | Free |
| GEO import | Import custom GEO datasets for SI analysis | User-provided | Academia |
| Vasculome query | Endothelial subtype-specific targets | Proprietary | Enterprise |
| BBB transport targets | Blood-brain barrier transporter discovery | Proprietary | Enterprise |
| Structure availability | Flag targets with PDB / AlphaFold structures | RCSB + AF2 | Free |
| Pathway enrichment | GO/KEGG enrichment of candidate list | MSigDB | Basic |

---

## 2. 快速开始

### 安装 / Installation

```bash
# 安装 LigandAI Python SDK
pip install ligandai>=2.1.0

# 验证安装
python -c "import ligandai; print(ligandai.__version__)"
```

### API Key 配置 / API Key Configuration

```python
import os
import ligandai

# 方式1：环境变量（推荐生产环境）
# Method 1: Environment variable (recommended for production)
os.environ["LIGANDAI_API_KEY"] = "lai-xxxxxxxxxxxxxxxxxxxxxxxx"

# 方式2：直接传入客户端
# Method 2: Pass directly to client
client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# 方式3：配置文件 ~/.ligandai/config.toml
# Method 3: Config file
# [auth]
# api_key = "lai-xxxxxxxxxxxxxxxxxxxxxxxx"

# 验证连接与 Tier 信息
# Verify connection and tier info
info = client.account.info()
print(f"Tier: {info.tier}")          # e.g. "academia"
print(f"Credits: {info.credits}")    # e.g. 45000
print(f"Org: {info.organization}")
```

---

## 3. 核心 API 详解

### 3.1 `client.discovery.tissue_markers()` — 组织标志物发现

```python
results = client.discovery.tissue_markers(
    tissue="liver",           # 目标组织（GTEx 标准名）
    receptor_only=True,       # 仅返回受体类基因
    si_threshold=2.0,         # SI 最低阈值
    top_n=50,                 # 返回前 50 个候选
    include_pdb=True,         # 标注 PDB 结构可用性
    expression_unit="tpm",    # 表达单位：tpm / rpkm / counts
)
```

**参数详解 / Parameter Reference:**

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `tissue` | `str` | 必填 | GTEx 组织名，支持 54 种标准组织名称 |
| `receptor_only` | `bool` | `True` | 过滤至受体类蛋白（GPCR/RTK/离子通道/转运体） |
| `si_threshold` | `float` | `2.0` | Specificity Index 最低阈值，建议 ≥2.0 |
| `top_n` | `int` | `20` | 返回候选数量上限 |
| `include_pdb` | `bool` | `True` | 是否查询 PDB/AlphaFold 结构可用性 |
| `expression_unit` | `str` | `"tpm"` | 表达量单位：`tpm` / `rpkm` / `counts` |
| `exclude_ubiquitous` | `bool` | `True` | 排除管家基因（如 ACTB, GAPDH） |
| `min_mean_tpm` | `float` | `1.0` | 目标组织最低平均表达量（TPM） |
| `receptor_classes` | `list` | `None` | 指定受体类别，如 `["GPCR", "RTK"]` |

### 3.2 `client.discovery.cell_type_markers()` — 细胞类型标志物（Academia+）

```python
# 需要 Academia 或更高 Tier
# Requires Academia tier or above
results = client.discovery.cell_type_markers(
    tissue="liver",
    cell_type="hepatocyte",       # 目标细胞类型
    atlas="tabula_sapiens",       # 数据集：tabula_sapiens / hca / custom
    receptor_only=True,
    si_threshold=3.0,             # 单细胞数据建议更高阈值
    top_n=30,
    resolution="leiden_1.0",      # 聚类分辨率（可选）
)
```

**支持的 Atlas 数据集 / Supported Atlas Datasets:**

| Atlas | 组织覆盖 | 细胞数 | Tier |
|-------|---------|--------|------|
| `tabula_sapiens` | 24 种人体组织 | ~500K | Academia |
| `hca` | Human Cell Atlas v3 | ~2.4M | Academia |
| `gtex_scrna` | GTEx 单细胞版 | ~200K | Academia |
| `custom` | 用户上传 GEO 数据 | 用户定义 | Academia |

### 3.3 GEO 数据集导入 / GEO Dataset Import

```python
# 导入 GEO 数据集进行自定义分析
# Import GEO dataset for custom analysis
geo_job = client.discovery.import_geo(
    accession="GSE115469",        # GEO 数据集编号
    species="human",              # 物种
    data_type="scrna",            # scrna / bulk
    cell_type_column="cell_type", # 元数据中细胞类型列名
    normalize=True,               # 自动归一化
)

# 等待导入完成（通常 2-5 分钟）
geo_job.wait(timeout=600)

# 使用导入的数据集进行分析
results = client.discovery.cell_type_markers(
    tissue="liver",
    cell_type="hepatocyte",
    atlas="custom",
    custom_dataset_id=geo_job.dataset_id,
    receptor_only=True,
)
```

---

## 4. Specificity Index (SI) 解释

### SI 算法原理

Specificity Index 是 LigandAI 专有的组织特异性评分算法，综合考量：

```
SI = log2(TPM_target / geometric_mean(TPM_all_tissues + 1)) × expression_weight
```

其中：
- `TPM_target`：目标组织的平均 TPM 表达量
- `geometric_mean(TPM_all_tissues)`：所有组织 TPM 的几何平均值
- `expression_weight`：基于表达量绝对值的权重因子（避免低表达基因的假高 SI）

### SI 值解读指南

| SI 值范围 | 特异性等级 | 推荐用途 | 脱靶风险 |
|-----------|-----------|---------|---------|
| SI ≥ 5.0 | 极高特异性 | 首选靶点，组织递送首选 | 极低 |
| 3.0 ≤ SI < 5.0 | 高特异性 | 优质靶点，适合大多数应用 | 低 |
| 2.0 ≤ SI < 3.0 | 中等特异性 | 可用靶点，需结合其他指标 | 中等 |
| 1.0 ≤ SI < 2.0 | 低特异性 | 谨慎使用，建议验证 | 较高 |
| SI < 1.0 | 无特异性 | 不推荐用于靶向设计 | 高 |

**最佳实践：** 对于靶向递送应用，建议 SI ≥ 3.0；对于治疗性肽段，SI ≥ 2.0 可接受，但需结合蛋白质表达验证。

---

## 5. 完整代码示例：肝脏靶点发现

```python
"""
完整示例：肝脏靶点发现 → 筛选受体 → SI 排名 → 导出
Full example: Liver target discovery → receptor filtering → SI ranking → export
"""

import ligandai
import pandas as pd
import json

# ── 1. 初始化客户端 / Initialize client ──────────────────────────────────────
client = ligandai.Client(api_key="lai-xxxxxxxxxxxxxxxxxxxxxxxx")

# 检查账户信息
account = client.account.info()
print(f"✓ Connected | Tier: {account.tier} | Credits: {account.credits:,}")

# ── 2. 组织标志物发现（GTEx 批量数据）/ Tissue marker discovery ──────────────
print("\n[Step 1] Discovering liver-specific receptor targets via GTEx...")

liver_targets = client.discovery.tissue_markers(
    tissue="liver",
    receptor_only=True,           # 仅受体类基因
    si_threshold=2.0,             # SI ≥ 2.0
    top_n=50,                     # 前 50 个候选
    include_pdb=True,             # 标注结构可用性
    exclude_ubiquitous=True,      # 排除管家基因
    min_mean_tpm=5.0,             # 最低表达量 5 TPM
    receptor_classes=["GPCR", "RTK", "transporter"],  # 指定受体类别
)

print(f"✓ Found {len(liver_targets)} candidates (SI ≥ 2.0, receptor-only)")

# ── 3. 转换为 DataFrame 便于分析 / Convert to DataFrame ──────────────────────
df = pd.DataFrame([
    {
        "rank":              t.rank,
        "gene":              t.gene,
        "si":                round(t.si, 3),
        "tissue_specificity": t.tissue_specificity,
        "receptor_class":    t.receptor_class,
        "mean_tpm":          round(t.mean_tpm, 2),
        "uniprot_id":        t.uniprot_id,
        "pdb_available":     t.pdb_available,
        "alphafold_available": t.alphafold_available,
    }
    for t in liver_targets
])

# ── 4. 多级筛选 / Multi-level filtering ──────────────────────────────────────
# 筛选1：高特异性（SI ≥ 3.0）
high_si = df[df["si"] >= 3.0].copy()
print(f"\n[Filter 1] High specificity (SI ≥ 3.0): {len(high_si)} targets")

# 筛选2：有结构信息（PDB 或 AlphaFold）
with_structure = df[df["pdb_available"] | df["alphafold_available"]].copy()
print(f"[Filter 2] With structure (PDB/AF2): {len(with_structure)} targets")

# 筛选3：GPCR 类受体（适合肽段设计）
gpcr_targets = df[df["receptor_class"] == "GPCR"].copy()
print(f"[Filter 3] GPCR class: {len(gpcr_targets)} targets")

# ── 5. 优先级排名 / Priority ranking ─────────────────────────────────────────
# 综合评分：SI 权重 60% + 结构可用性 20% + 表达量 20%
df["structure_score"] = (df["pdb_available"].astype(int) * 1.0 +
                          df["alphafold_available"].astype(int) * 0.5)
df["expr_score"] = df["mean_tpm"].clip(upper=100) / 100.0
df["priority_score"] = (
    df["si"] * 0.6 +
    df["structure_score"] * 0.2 * 5 +  # 归一化到 SI 量级
    df["expr_score"] * 0.2 * 5
)
df_ranked = df.sort_values("priority_score", ascending=False).reset_index(drop=True)
df_ranked["priority_rank"] = df_ranked.index + 1

# ── 6. 输出结果 / Display results ────────────────────────────────────────────
print("\n" + "="*80)
print("TOP 10 LIVER RECEPTOR TARGETS (by Priority Score)")
print("="*80)
display_cols = ["priority_rank", "gene", "si", "tissue_specificity",
                "receptor_class", "mean_tpm", "pdb_available", "uniprot_id"]
print(df_ranked[display_cols].head(10).to_string(index=False))

# ── 7. 导出结果 / Export results ─────────────────────────────────────────────
df_ranked.to_csv("liver_targets_ranked.csv", index=False)
print("\n✓ Results saved to liver_targets_ranked.csv")

# 导出 JSON 供下游使用
top_targets_json = df_ranked.head(10)[["gene", "si", "uniprot_id"]].to_dict("records")
with open("top10_liver_targets.json", "w") as f:
    json.dump(top_targets_json, f, indent=2)
print("✓ Top 10 targets saved to top10_liver_targets.json")

# ── 8. 单细胞验证（Academia+ Tier）/ scRNA validation ────────────────────────
if account.tier in ["academia", "pro", "enterprise"]:
    print("\n[Step 2] Validating top targets with scRNA-seq (hepatocyte-specific)...")
    
    top_genes = df_ranked.head(5)["gene"].tolist()
    
    scrna_results = client.discovery.cell_type_markers(
        tissue="liver",
        cell_type="hepatocyte",
        atlas="tabula_sapiens",
        receptor_only=True,
        si_threshold=2.5,
        top_n=20,
        gene_filter=top_genes,    # 仅验证 GTEx 发现的候选
    )
    
    validated_genes = [r.gene for r in scrna_results]
    print(f"✓ scRNA validation: {len(validated_genes)}/{len(top_genes)} confirmed")
    print(f"  Validated: {', '.join(validated_genes)}")
else:
    print(f"\n[Skip] scRNA validation requires Academia+ tier (current: {account.tier})")

# ── 9. 传递给肽段生成 / Pass to peptide generation ──────────────────────────
best_target = df_ranked.iloc[0]["gene"]
print(f"\n✓ Best target for peptide design: {best_target}")
print(f"  → Next step: use ligandai-peptide skill with gene='{best_target}'")
```

---

## 6. 输出字段说明表

| 字段名 | 类型 | 说明 | 示例值 |
|--------|------|------|--------|
| `gene` | `str` | HGNC 基因符号 | `"ASGR1"` |
| `si` | `float` | Specificity Index 值（越高越特异） | `6.82` |
| `tissue_specificity` | `str` | 特异性等级：`high` / `medium` / `low` | `"high"` |
| `receptor_class` | `str` | 受体分类：`GPCR` / `RTK` / `ion_channel` / `transporter` / `other` | `"transporter"` |
| `mean_tpm` | `float` | 目标组织平均 TPM 表达量 | `142.7` |
| `rank` | `int` | SI 排名（1 = 最高特异性） | `1` |
| `uniprot_id` | `str` | UniProt 蛋白 ID | `"P07306"` |
| `pdb_available` | `bool` | 是否有 PDB 实验结构 | `True` |
| `alphafold_available` | `bool` | 是否有 AlphaFold 预测结构 | `True` |
| `n_tissues_expressed` | `int` | 在多少种组织中表达（TPM > 1） | `3` |
| `max_other_tpm` | `float` | 非目标组织中最高 TPM | `8.2` |
| `go_terms` | `list[str]` | 相关 GO 生物学过程 | `["receptor activity", "endocytosis"]` |
| `pathway` | `str` | 主要信号通路 | `"Wnt signaling"` |

---

## 7. Tier 要求说明

| 功能 | Free | Basic | Academia | Pro | Enterprise |
|------|------|-------|----------|-----|------------|
| GTEx 批量 RNA-seq 分析 | ✅ | ✅ | ✅ | ✅ | ✅ |
| 受体类过滤 | ✅ | ✅ | ✅ | ✅ | ✅ |
| 结构可用性标注 | ✅ | ✅ | ✅ | ✅ | ✅ |
| 单细胞 RNA-seq 图谱 | ❌ | ❌ | ✅ | ✅ | ✅ |
| GEO 数据集导入 | ❌ | ❌ | ✅ | ✅ | ✅ |
| 通路富集分析 | ❌ | ✅ | ✅ | ✅ | ✅ |
| Vasculome 内皮亚型数据库 | ❌ | ❌ | ❌ | ❌ | ✅ |
| BBB 转运蛋白数据库 | ❌ | ❌ | ❌ | ❌ | ✅ |
| 每次查询最大候选数 | 20 | 50 | 100 | 200 | 无限制 |
| 自定义 SI 阈值 | ❌ | ✅ | ✅ | ✅ | ✅ |

> **注意：** Discovery 模块本身不消耗 Credits，但后续的结构解析和肽段生成会消耗 Credits。

---

## 8. 下游衔接

发现靶点后，推荐使用以下技能进行后续分析：

```python
# 将发现的靶点传递给肽段生成技能
# Pass discovered target to peptide generation skill

best_gene = "ASGR1"  # 来自 Discovery 结果

# → 使用 ligandai-peptide 技能生成靶向肽段
peptides = client.peptides.generate(
    gene=best_gene,
    num_peptides=300,
    targeting_strategy="pocket_targeted",
    auto_fold=True,
    top_n_fold=10,
)

# → 使用 ligandai-folding 技能进行结构验证
structure = client.structures.get(gene=best_gene)
pocket = client.structures.analyze(gene=best_gene, analysis_depth="standard")
```

**推荐工作流：**
```
ligandai-discovery
    ↓ (top targets)
ligandai-folding    ← 解析靶点结构，识别结合口袋
    ↓ (pocket info)
ligandai-peptide    ← 生成靶向肽段
    ↓ (peptide candidates)
ligandai-scoring    ← 热力学评分与筛选
    ↓ (top candidates)
ligandai-bivalent   ← 双价设计（可选，Pro+）
```

---

## 9. 参考资料

- **LigandAI 官方文档：** https://docs.ligandai.com/discovery
- **GTEx Portal：** https://gtexportal.org/home/
- **Human Cell Atlas：** https://www.humancellatlas.org/
- **Tabula Sapiens：** https://tabula-sapiens-portal.ds.czbiohub.org/
- **UniProt 受体分类：** https://www.uniprot.org/
- **SI 算法论文：** LigandAI Technical Report v2.1 (2024)

> **技能版本：** v1.0.0 | **SDK 最低版本：** ligandai>=2.1.0 | **最后更新：** 2025-01
