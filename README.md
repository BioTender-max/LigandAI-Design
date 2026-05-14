<div align="center">

<a href="https://ligandai.com">
  <img src="assets/banner.png" alt="LigandAI-Design" width="100%"/>
</a>

<br/>

[![License](https://img.shields.io/badge/License-Apache_2.0-green?style=flat-square)](LICENSE)
[![SDK](https://img.shields.io/badge/ligandai-v0.5.3-blue?style=flat-square&logo=python&logoColor=white)](https://ligandai.com/docs)
[![Skills](https://img.shields.io/badge/Skills-5-blueviolet?style=flat-square)](#skills)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen?style=flat-square)](https://github.com/BioTender-max/LigandAI-Design/pulls)

</div>

---

> **LigandAI-Design** is an open-source collection of composable design skills for end-to-end peptide binder discovery — from tissue-specific target identification to thermodynamic candidate ranking — built on the [LigandAI](https://ligandai.com) Predictive Interactomics™ platform.

---

## Table of Contents

- [Overview](#overview)
- [Skills](#skills)
- [Pipeline](#pipeline)
- [Quick Start](#quick-start)
- [Examples](#examples)
- [Pricing & Tiers](#pricing--tiers)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)

---

## Overview

LigandAI exposes four core engines through a unified Python SDK:

| Engine | Description |
|--------|-------------|
| **ReceptorDB** | 3,500+ curated receptor complexes with pre-computed pockets and MSA |
| **LigandForge v6.5** | Discrete diffusion peptide generator — >1,000 peptides/sec on B200 |
| **Boltz-2** | Multi-trajectory complex structure prediction with MSA caching |
| **DeltaForge** | Thermodynamic ΔG/Kd prediction — r = 0.83 vs SPR/ITC experiments |

---

## Skills

<a name="skills"></a>

| # | Skill | Core Function | Min Tier |
|---|-------|---------------|----------|
| 1 | [`ligandai-discovery`](skills/ligandai-discovery/SKILL.md) | Tissue-specific receptor discovery via GTEx/scRNA Specificity Index | Free |
| 2 | [`ligandai-folding`](skills/ligandai-folding/SKILL.md) | Boltz-2 complex folding · 4-tier structure fallback · MSA caching | Free |
| 3 | [`ligandai-peptide`](skills/ligandai-peptide/SKILL.md) | LigandForge v6.5 · pocket-targeted / scaffold / deimmunized / drug-conjugate | Basic |
| 4 | [`ligandai-scoring`](skills/ligandai-scoring/SKILL.md) | DeltaForge ΔG/Kd · LigandIQ pre-screen · multi-metric ranking (free) | Free |
| 5 | [`ligandai-bivalent`](skills/ligandai-bivalent/SKILL.md) | Induced proximity · PROTAC-like bispecific design *(Beta)* | Pro |

---

## Pipeline

```
ligandai-discovery
  └─ GTEx / scRNA Specificity Index → ranked receptor targets
        │
ligandai-folding
  └─ ReceptorDB → Boltz-2 → RCSB → AlphaFold (4-tier fallback)
  └─ Pocket detection · druggability scoring · key residues
        │
ligandai-peptide
  └─ LigandForge v6.5 · up to 1,000 candidates · SSE streaming
  └─ Strategies: pocket_targeted · scaffold_templated · deimmunized · drug_conjugate
        │
ligandai-scoring
  └─ LigandIQ pre-screen (free) → fold top 50 → DeltaForge precision score (free)
  └─ Outputs: ΔG · Kd · iPSAE · iPTM · immunogenicity · solubility
        │
ligandai-bivalent  ← optional, Pro+
  └─ Bifunctional linker design · ternary complex folding
        │
        ▼
   Top candidate peptides → experimental validation
```

---

## Quick Start

```bash
pip install ligandai==0.5.3
```

```python
import ligandai

client = ligandai.Client(api_key="lgai_pro_...")   # ligandai.com/sdk

# 1. Discover liver-specific receptor targets
targets = client.discovery.tissue_markers(tissue="liver", receptor_only=True, top_n=10)

# 2. Resolve structure + analyze binding pocket
pocket = client.structures.analyze(gene=targets[0].gene, analysis_depth="standard")

# 3. Generate 300 peptide candidates
job = client.peptides.generate(
    gene=targets[0].gene,
    num_peptides=300,
    targeting_strategy="pocket_targeted",
    auto_fold=True,
    top_n_fold=10,
)

# 4. Stream progress and collect results
for event in job.stream():
    if event.type == "elite_hit":
        print(f"Elite: {event.sequence}  ΔG={event.delta_g:.2f}")
    elif event.type == "generation_complete":
        break

results = job.results()
print(f"Top peptide: {results.peptides[0].sequence}")
print(f"ΔG = {results.peptides[0].delta_g:.2f} kcal/mol  |  Kd = {results.peptides[0].kd:.1f} nM")
```

---

## Examples

| Script | Description |
|--------|-------------|
| [`examples/ligandai_e2e_pipeline.py`](examples/ligandai_e2e_pipeline.py) | Full 7-step end-to-end pipeline with Rich progress display |
| [`examples/ligandai_async_parallel.py`](examples/ligandai_async_parallel.py) | Async parallel design across 3 targets with `asyncio.Semaphore` |
| [`examples/ligandai_scoring_filter.py`](examples/ligandai_scoring_filter.py) | Multi-metric scoring, hard filters, Pareto frontier analysis + matplotlib |

---

## Pricing & Tiers

### API Key Tiers

| Tier | Key Prefix | Peptides | Folds / Target | Notable Features |
|------|-----------|----------|----------------|-----------------|
| Free | `lgai_free_*` | 100 | 1 | GTEx discovery, basic scoring |
| Basic | `lgai_basic_*` | 300 | 25 | Deimmunization, scaffold templates |
| Academia | `lgai_edu_*` | 300 | 50 | scRNA atlases, GEO import, auto-fold |
| Pro | `lgai_pro_*` | 1,000 | 100 | Bivalent design, stability/half-life scoring |
| Enterprise | `lgai_ent_*` | 5,000 | Unlimited | BBB vasculome, dedicated GPU slots |

### Credit Costs

| Operation | Cost | Notes |
|-----------|------|-------|
| Peptide generation | 100 credits / peptide | — |
| Boltz-2 fold (50 steps, 1 trajectory) | 100 credits | +100 per additional 50 steps |
| MSA cache hit | −20% discount | Same target, second call onward |
| DeltaForge / LigandIQ scoring | **Free** | CPU-only, unlimited |
| Discovery (SI ranking) | **Free** | CPU-only, unlimited |

> Two-stage screening example: generate 300 → LigandIQ pre-screen to 50 (free) → fold 50 → DeltaForge score (free). Total: **5,000 credits** vs 30,000 without pre-screening — **83% savings**.

---

## Repository Structure

```
LigandAI-Design/
├── assets/
│   └── banner.png
├── skills/
│   ├── ligandai-discovery/SKILL.md
│   ├── ligandai-peptide/SKILL.md
│   ├── ligandai-folding/SKILL.md
│   ├── ligandai-scoring/SKILL.md
│   └── ligandai-bivalent/SKILL.md
├── examples/
│   ├── ligandai_e2e_pipeline.py
│   ├── ligandai_async_parallel.py
│   └── ligandai_scoring_filter.py
├── README.md
└── LICENSE
```

---

## Contributing

Issues and pull requests are welcome. Please open an issue first for major changes.

- **Bug reports** — include SDK version, tier, and a minimal reproducible example
- **New skills** — follow the `SKILL.md` frontmatter schema in existing skills
- **Examples** — add a row to the Examples table and keep scripts self-contained

---

<div align="center">
<sub>Built for the protein design community · <a href="https://ligandai.com">ligandai.com</a></sub>
</div>
