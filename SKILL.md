---
name: ccc-skill
description: Domain expertise for cell-cell communication (CCC) analysis. Intent-to-tool routing, workflow patterns, code templates, and pitfall warnings for LIANA+, CellPhoneDB, CellChat, NicheNet, and COMMOT across scRNA-seq and spatial transcriptomics data.
license: MIT
homepage: https://github.com/Agents365-ai/ccc-skill
platforms: [macos, linux, windows]
metadata: {"author":"Agents365-ai","version":"0.1.0"}
---

# Cell-Cell Communication Analysis Skill

You are an expert in cell-cell communication (CCC) inference from single-cell and spatial transcriptomics. Use this skill to select the right tool, design the workflow, and generate correct code.

## Decision Tree

Use this decision tree to recommend the right tool(s). Multiple tools can be combined.

### Step 1: What is the primary analysis goal?

| Goal | Recommended Tool(s) | Guide |
|------|---------------------|-------|
| Ligand-receptor interaction inference (steady-state) | LIANA+ `rank_aggregate` or CellPhoneDB | [liana.md](guides/liana.md), [cellphonedb.md](guides/cellphonedb.md) |
| Ligand-receptor with signaling pathway hierarchy | CellChat | [cellchat.md](guides/cellchat.md) |
| Ligand → downstream target gene prediction | NicheNet / MultiNicheNet | [nichenet.md](guides/nichenet.md) |
| Spatial CCC with distance decay | COMMOT (optimal transport) or LIANA+ (bivariate) | [commot.md](guides/commot.md), [liana.md](guides/liana.md) |
| Multi-sample / multi-condition CCC comparison | LIANA+ (tensor/MOFA+) or CellChat (merged comparison) | [liana.md](guides/liana.md), [cellchat.md](guides/cellchat.md) |
| Multi-view spatial learning (what predicts what) | LIANA+ MISTy | [liana.md](guides/liana.md) |
| Spatial signaling direction / vector fields | COMMOT | [commot.md](guides/commot.md) |

### Step 2: What data type?

| Data Type | Compatible Tools |
|-----------|-----------------|
| scRNA-seq (Python/AnnData) | LIANA+, CellPhoneDB |
| scRNA-seq (R/Seurat) | CellChat, NicheNet |
| Spatial spot-based (Visium) | LIANA+ bivariate, COMMOT, CellChat v2 spatial |
| Spatial single-cell (MERFISH, Xenium) | COMMOT, LIANA+ bivariate/inflow, CellChat v3 (SpatialCellChat) |
| Multi-modal (MuData) | LIANA+ |
| Multi-sample comparison | LIANA+ tensor/MOFA+, CellChat merged, MultiNicheNet |

### Step 3: Language preference?

| Language | Tools |
|----------|-------|
| Python only | LIANA+, CellPhoneDB, COMMOT |
| R only | CellChat, NicheNet/MultiNicheNet |
| Either | Use decision tree above; recommend Python for AnnData pipelines, R for Seurat pipelines |

## Tool Comparison at a Glance

| Feature | LIANA+ | CellPhoneDB | CellChat | NicheNet | COMMOT |
|---------|--------|-------------|----------|----------|--------|
| Language | Python | Python | R | R | Python |
| Input | AnnData/MuData | AnnData/matrix | Seurat/SCE/matrix | Seurat | AnnData |
| LR database | consensus (multi-DB) | CellPhoneDB v5 (curated) | CellChatDB v2 | Zenodo model | CellChat/CellPhoneDB |
| Statistical test | per-method (permutation or analytical) | permutation | permutation (law of mass action) | AUPR enrichment | optimal transport |
| Spatial support | bivariate, inflow, MISTy | no | v2 distance-weighted, v3 single-cell | no (ligand activity only) | native (OT with distance) |
| Multi-sample | tensor, MOFA+ | no | merged object comparison | MultiNicheNet (pseudobulk) | group analysis |
| Downstream targets | no (LR only) | no | no | yes (core feature) | communication-dependent genes |
| Unique strength | unified multi-method meta-analysis | curated DB + heteromeric complexes | pathway hierarchy + rich visualization | ligand→TF→target prediction | spatial direction vector fields |

## Universal CCC Principles

### LR Database Selection
- **consensus** (LIANA+ default): union of multiple databases — broadest coverage, recommended for discovery
- **CellPhoneDB**: manually curated, strong on heteromeric complexes — best for confidence
- **CellChatDB**: pathway-annotated — best when you need pathway-level grouping
- **Custom**: always an option; filter to expressed genes for performance

### Preprocessing Requirements (all tools)
- Normalized, log-transformed expression (NOT raw counts, NOT scaled/z-scored)
- Cell type annotations in metadata
- For spatial: coordinates in the correct unit system

### Interpreting Results
- **No single tool is ground truth** — consider running 2+ tools and looking for consensus
- **Magnitude vs specificity**: high expression doesn't mean specific; highly specific doesn't mean strong
- **Permutation p-values**: affected by cell type size imbalance — small populations have less power
- **Spatial distance**: paracrine signals (~100-500um) vs juxtacrine/contact (~10-50um) need different cutoffs

### Common Mistakes to Avoid
1. Using scaled/z-scored data (destroys zeros needed for proportion calculation)
2. Using raw integer counts (most tools expect normalized data)
3. Ignoring `expr_prop` / expression threshold — lowering it inflates false positives
4. Comparing CCC across conditions without multi-sample statistics (pseudoreplication)
5. Over-interpreting low-ranked interactions from a single method

## How to Use Guides

Each guide in `guides/` contains:
1. **When to use** — scenarios where this tool is the right choice
2. **Installation** — quick setup
3. **Core workflow** — step-by-step code template
4. **Key parameters** — what matters and what to tune
5. **Output format** — where results are stored and how to access them
6. **Visualization** — plotting code
7. **Pitfalls** — tool-specific gotchas

Read the relevant guide(s) based on the decision tree above, then generate code following the templates.
