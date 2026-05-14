<div align="center">

<img src="https://ligandai.com/favicon.ico" width="80" alt="LigandAI" />

# LigandAI-Design

**基于 LigandAI 平台的端到端蛋白/肽段设计技能包**
**End-to-End Protein & Peptide Binder Design Skills powered by LigandAI**

---

[![Skills](https://img.shields.io/badge/Skills-5-4CAF50?style=flat-square&logo=bookstack&logoColor=white)](#skills)
[![SDK](https://img.shields.io/badge/SDK-ligandai%20v0.5.3-2196F3?style=flat-square&logo=python&logoColor=white)](https://ligandai.com/docs)
[![License](https://img.shields.io/badge/License-Apache%202.0-yellow?style=flat-square)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-LigandAI-blueviolet?style=flat-square)](https://ligandai.com)
[![Last Updated](https://img.shields.io/badge/Last%20Updated-2026--05--15-orange?style=flat-square)]()

**[English](#english) | [中文](#中文)**

</div>

---

## 中文

### 简介

本仓库收录了基于 [LigandAI](https://ligandai.com) 平台的 **5 个专业蛋白设计技能（Skills）**，覆盖从靶点发现到候选肽段筛选的完整药物发现流程。

LigandAI 是 Ligandal 公司的 Predictive Interactomics™ 平台的开发者接口，集成了：
- **ReceptorDB** — 3,500+ 精选受体复合物数据库
- **LigandForge v6.5** — 离散扩散肽段生成器（>1,000 肽段/秒）
- **Boltz-2** — 多轨迹复合物结构预测
- **DeltaForge** — 热力学 dG/Kd 预测（r=0.83）

### 技能列表

| # | 技能名称 | Skill ID | 核心功能 | Tier 要求 |
|---|---------|----------|---------|----------|
| 1 | 靶点发现 | `ligandai-discovery` | GTEx/scRNA 组织特异性受体发现，SI 排名 | Free+ |
| 2 | 肽段生成 | `ligandai-peptide` | LigandForge v6.5 口袋靶向/支架/去免疫原性设计 | Basic+ |
| 3 | 复合物折叠 | `ligandai-folding` | Boltz-2 多轨迹折叠，MSA 缓存，糖基化支持 | Free+ |
| 4 | 热力学评分 | `ligandai-scoring` | DeltaForge dG/Kd，LigandIQ 预筛，多维度评分 | Free |
| 5 | 双价设计 | `ligandai-bivalent` | 诱导邻近，PROTAC 类双特异性设计 | Pro+ |

### 完整设计流程

```
靶点发现 (ligandai-discovery)
    ↓  GTEx/scRNA 组织特异性受体 + SI 排名
结构解析 (ligandai-folding)
    ↓  ReceptorDB → Boltz-2 → RCSB → AlphaFold 四级回退
口袋分析 (ligandai-folding)
    ↓  推荐结合口袋 + 可成药性评分
肽段生成 (ligandai-peptide)
    ↓  LigandForge v6.5 口袋靶向生成 + 流式折叠
热力学评分 (ligandai-scoring)
    ↓  DeltaForge dG/Kd + 免疫原性 + 溶解度
候选优化 (ligandai-bivalent)  [可选，Pro+]
    ↓  双价/双特异性设计
Top 候选肽段输出
```

### 快速开始

```bash
# 1. 安装 SDK
pip install ligandai==0.5.3

# 2. 配置 API Key（从 https://ligandai.com/sdk 获取）
export LIGANDAI_API_KEY="lgai_pro_..."

# 3. 运行端到端示例
python examples/ligandai_e2e_pipeline.py
```

### 示例脚本

| 脚本 | 描述 | 适用场景 |
|------|------|---------|
| `examples/ligandai_e2e_pipeline.py` | 完整端到端流程（559行） | 单靶点完整设计 |
| `examples/ligandai_async_parallel.py` | 异步并行多靶点设计（481行） | 多靶点批量筛选 |
| `examples/ligandai_scoring_filter.py` | 多维度评分筛选+帕累托分析（655行） | 候选肽段精细筛选 |

---

## English

### Introduction

This repository contains **5 professional protein design skills** built on the [LigandAI](https://ligandai.com) platform, covering the complete drug discovery pipeline from target identification to candidate peptide selection.

LigandAI is the developer interface to Ligandal's Predictive Interactomics™ platform, integrating:
- **ReceptorDB** — 3,500+ curated receptor complexes
- **LigandForge v6.5** — Discrete diffusion peptide generator (>1,000 peptides/sec on B200)
- **Boltz-2** — Multi-trajectory complex structure prediction with MSA caching
- **DeltaForge** — Thermodynamic dG/Kd prediction (r=0.83 vs experimental)

### Skills Overview

| # | Skill Name | Skill ID | Core Function | Tier |
|---|-----------|----------|---------------|------|
| 1 | Target Discovery | `ligandai-discovery` | GTEx/scRNA tissue-specific receptor discovery with SI ranking | Free+ |
| 2 | Peptide Generation | `ligandai-peptide` | LigandForge v6.5 pocket-targeted / scaffold / deimmunized design | Basic+ |
| 3 | Complex Folding | `ligandai-folding` | Boltz-2 multi-trajectory folding, MSA caching, glycosylation | Free+ |
| 4 | Thermodynamic Scoring | `ligandai-scoring` | DeltaForge dG/Kd, LigandIQ pre-fold screening, multi-metric ranking | Free |
| 5 | Bivalent Design | `ligandai-bivalent` | Induced proximity, PROTAC-like bispecific design | Pro+ |

### Full Design Pipeline

```
Target Discovery (ligandai-discovery)
    ↓  GTEx/scRNA tissue-specific receptors + SI ranking
Structure Resolution (ligandai-folding)
    ↓  ReceptorDB → Boltz-2 → RCSB → AlphaFold (4-tier fallback)
Pocket Analysis (ligandai-folding)
    ↓  Recommended binding pocket + druggability score
Peptide Generation (ligandai-peptide)
    ↓  LigandForge v6.5 pocket-targeted generation + streaming fold
Thermodynamic Scoring (ligandai-scoring)
    ↓  DeltaForge dG/Kd + immunogenicity + solubility
Candidate Optimization (ligandai-bivalent)  [optional, Pro+]
    ↓  Bivalent / bispecific design
Top Candidate Peptides Output
```

### Quick Start

```bash
# 1. Install SDK
pip install ligandai==0.5.3

# 2. Set API Key (generate at https://ligandai.com/sdk)
export LIGANDAI_API_KEY="lgai_pro_..."

# 3. Run end-to-end example
python examples/ligandai_e2e_pipeline.py
```

### API Key Tiers

| Prefix | Tier | Rate Limit | Concurrent GPUs | Key Features |
|--------|------|-----------|-----------------|--------------|
| `lgai_free_*` | Free | 10/min | 2 | Browse, 1 fold, 100 peptides |
| `lgai_basic_*` | Basic | 20/min | 5 | 300 peptides, 25 folds/target |
| `lgai_edu_*` | Academia | 30/min | 12 | scRNA, GEO import, 50 folds/target |
| `lgai_pro_*` | Pro | 60/min | 24 | 1000 peptides, bivalent design |
| `lgai_ent_*` | Enterprise | 300/min | 96+ | 5000 peptides, BBB vasculome |

### Credit Costs

| Operation | Cost | Notes |
|-----------|------|-------|
| Peptide generation | 100 credits/peptide | Tier-visible count only |
| Boltz-2 fold | 100 credits/trajectory | 50-step base |
| DeltaForge / LigandIQ | **Free** | CPU-only |
| Discovery (SI) | **Free** | CPU-only |

---

## Repository Structure

```
LigandAI-Design/
├── skills/
│   ├── ligandai-discovery/
│   │   └── SKILL.md          # 靶点发现 / Target Discovery
│   ├── ligandai-peptide/
│   │   └── SKILL.md          # 肽段生成 / Peptide Generation
│   ├── ligandai-folding/
│   │   └── SKILL.md          # 复合物折叠 / Complex Folding
│   ├── ligandai-scoring/
│   │   └── SKILL.md          # 热力学评分 / Thermodynamic Scoring
│   └── ligandai-bivalent/
│       └── SKILL.md          # 双价设计 / Bivalent Design
├── examples/
│   ├── ligandai_e2e_pipeline.py       # 端到端完整流程
│   ├── ligandai_async_parallel.py     # 异步并行多靶点
│   └── ligandai_scoring_filter.py     # 多维度评分筛选
├── README.md
└── LICENSE
```

---

## References

- [LigandAI Documentation](https://ligandai.com/docs)
- [LigandAI Quickstart](https://ligandai.com/docs/quickstart)
- [LigandAI Authentication](https://ligandai.com/docs/authentication)
- [LigandAI Pricing & Tiers](https://ligandai.com/docs/pricing)
- [ReceptorDB](https://receptordb.com)

---

<div align="center">

Built with ❤️ for the protein design community · Powered by [LigandAI](https://ligandai.com)

*Skills for [Biomni](https://phylo.bio) by Phylo*

</div>
