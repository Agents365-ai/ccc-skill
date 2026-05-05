---
name: ccc-skill
description: Domain expertise for cell-cell communication (CCC) analysis. Intent-to-tool routing, workflow patterns, code templates, and pitfall warnings for LIANA+, CellPhoneDB, CellChat, NicheNet, COMMOT, Squidpy, MEBOCOST, DIALOGUE, and 10+ additional tools across scRNA-seq and spatial transcriptomics data.
license: MIT
homepage: https://github.com/Agents365-ai/ccc-skill
platforms: [macos, linux, windows]
metadata: {"author":"Agents365-ai","version":"0.2.0"}
---

# Cell-Cell Communication Analysis Skill

You are an expert in cell-cell communication (CCC) inference from single-cell and spatial transcriptomics. Use this skill to select the right tool, design the workflow, and generate correct code.

### Step 0. Update check (notify, don't pull) — first use per conversation

Throttle to one check per 24 hours per installation; never mutate the skill directory without explicit user consent.

1. If `<this-skill-dir>/.last_update` exists and is less than 24 hours old, skip this step entirely.

2. Otherwise, fetch the latest tag from upstream:

   ```bash
   git -C <this-skill-dir> ls-remote --tags origin 'v*' 2>/dev/null \
     | awk '{print $2}' | sed 's|refs/tags/||' \
     | sort -V | tail -1
   ```

3. Compare with this skill's `metadata.version` from the frontmatter. If the upstream tag is strictly newer (semver), tell the user one line and ask:

   > "A newer version of this skill is available: vX.Y.Z → vA.B.C. Want me to `git pull`?"

   If they say yes, run `git -C <this-skill-dir> pull --ff-only`. Refresh `.last_update` either way so the prompt doesn't repeat for 24 hours.

4. If upstream is the same or older, refresh `.last_update` silently and continue.

5. On any failure (offline, not a git checkout — e.g. ClawHub-installed copy, read-only path, no permission), swallow the error silently and continue with the user's task. Do not mention the failure.

## Decision Tree

Use this decision tree to recommend the right tool(s). Multiple tools can be combined.

### Step 1: What is the primary analysis goal?

| Goal | Recommended Tool(s) | Guide |
|------|---------------------|-------|
| Ligand-receptor inference (steady-state) | LIANA+ `rank_aggregate` or CellPhoneDB | [liana.md](guides/liana.md), [cellphonedb.md](guides/cellphonedb.md) |
| LR with signaling pathway hierarchy | CellChat | [cellchat.md](guides/cellchat.md) |
| Ligand → downstream target gene prediction | NicheNet / MultiNicheNet | [nichenet.md](guides/nichenet.md) |
| Spatial CCC with distance decay | COMMOT or LIANA+ bivariate | [commot.md](guides/commot.md), [liana.md](guides/liana.md) |
| Spatial CCC within scanpy/scverse pipeline | Squidpy `sq.gr.ligrec()` | [squidpy.md](guides/squidpy.md) |
| Spatial CCC + H&E morphology | stLearn — [docs](https://stlearn.readthedocs.io/) |
| Spatial CCC + knowledge graph pathway cascade | SpaTalk — [GitHub](https://github.com/ZJUFanLab/SpaTalk) |
| Spatial signaling direction / vector fields | COMMOT | [commot.md](guides/commot.md) |
| Multi-sample / multi-condition comparison | LIANA+ tensor/MOFA+, CellChat merged, MultiNicheNet | [liana.md](guides/liana.md), [cellchat.md](guides/cellchat.md), [nichenet.md](guides/nichenet.md) |
| Cross-context CCC pattern discovery (tensor) | cell2cell / Tensor-cell2cell — [GitHub](https://github.com/earmingol/cell2cell) |
| Differential CCC network analysis | CrossTalkeR — [GitHub](https://github.com/CostaLab/CrossTalkeR) |
| Differential CCC with aging focus | scDiffCom — [GitHub](https://github.com/CyrilLagger/scDiffCom) |
| Multi-view spatial learning | LIANA+ MISTy | [liana.md](guides/liana.md) |
| Metabolite-mediated CCC (non-protein) | MEBOCOST | [mebocost.md](guides/mebocost.md) |
| Multicellular coordination programs | DIALOGUE | [dialogue.md](guides/dialogue.md) |
| Causal signal flow inference (post-CCC) | FlowSig — [GitHub](https://github.com/axelalmet/flowsig) |
| GRN + TF perturbation (downstream of CCC) | CellOracle — [GitHub](https://github.com/morris-lab/CellOracle) |
| Multi-omic GRN (scRNA+scATAC) | SCENIC+ — [GitHub](https://github.com/aertslab/scenicplus) |
| Single-cell resolution CCC (not aggregated) | NICHES or Scriabin — [NICHES](https://github.com/msraredon/NICHES), [Scriabin](https://github.com/BlishLab/scriabin) |
| Neural-specific CCC (brain data) | NeuronChat — [GitHub](https://github.com/Wei-BioMath/NeuronChat) |

### Step 2: What data type?

| Data Type | Compatible Tools |
|-----------|-----------------|
| scRNA-seq (Python/AnnData) | LIANA+, CellPhoneDB, MEBOCOST, Squidpy |
| scRNA-seq (R/Seurat) | CellChat, NicheNet, DIALOGUE, NICHES, Scriabin |
| Spatial spot-based (Visium) | LIANA+ bivariate, COMMOT, CellChat v2, Squidpy, stLearn, SpaTalk |
| Spatial single-cell (MERFISH, Xenium) | COMMOT, LIANA+ bivariate/inflow, CellChat v3 |
| Multi-modal (MuData) | LIANA+ |
| Multi-sample comparison | LIANA+ tensor/MOFA+, CellChat merged, MultiNicheNet, cell2cell, scDiffCom |
| Brain / neuronal | NeuronChat (specialized LR database) |

### Step 3: Language preference?

| Language | Tools |
|----------|-------|
| Python | LIANA+, CellPhoneDB, COMMOT, Squidpy, MEBOCOST, CellOracle, SCENIC+, FlowSig, stLearn, cell2cell |
| R | CellChat, NicheNet/MultiNicheNet, DIALOGUE, NICHES, Scriabin, SpaTalk, CrossTalkeR, scDiffCom, NeuronChat |

## Tools with Full Guides (8)

| Tool | Language | Guide | Unique Strength |
|------|----------|-------|-----------------|
| **LIANA+** | Python | [liana.md](guides/liana.md) | Multi-method meta-analysis + spatial + tensor/MOFA+ |
| **CellPhoneDB** | Python | [cellphonedb.md](guides/cellphonedb.md) | Curated DB + heteromeric complexes + v5 scoring |
| **CellChat** | R | [cellchat.md](guides/cellchat.md) | Pathway hierarchy + rich visualization + comparison |
| **NicheNet** | R | [nichenet.md](guides/nichenet.md) | Ligand→TF→target prediction + MultiNicheNet |
| **COMMOT** | Python | [commot.md](guides/commot.md) | Optimal transport spatial CCC + vector fields |
| **Squidpy** | Python | [squidpy.md](guides/squidpy.md) | scverse ecosystem LR analysis, zero-friction spatial |
| **MEBOCOST** | Python | [mebocost.md](guides/mebocost.md) | Metabolite-mediated CCC (non-protein signals) |
| **DIALOGUE** | R | [dialogue.md](guides/dialogue.md) | Multicellular programs (cross-cell-type coordination) |

## Additional Tools (brief reference)

These tools are referenced in the decision tree but do not have dedicated guides. Use official documentation.

| Tool | Language | Stars | When to Use | Link |
|------|----------|-------|-------------|------|
| **CellOracle** | Python | 440 | GRN + in silico TF perturbation → predict cell state after signal | [GitHub](https://github.com/morris-lab/CellOracle) |
| **SCENIC+** | Python | 251 | Multi-omic GRN (scRNA+scATAC) → regulatory context of CCC | [GitHub](https://github.com/aertslab/scenicplus) |
| **stLearn** | Python | 244 | Spatial CCC combining H&E morphology + expression | [GitHub](https://github.com/BiomedicalMachineLearning/stLearn) |
| **MultiNicheNet** | R | 185 | Multi-sample differential CCC (pseudobulk + edgeR) | [GitHub](https://github.com/saeyslab/multinichenetr) |
| **Scriabin** | R | 106 | Single-cell resolution CCC, atlas-scale | [GitHub](https://github.com/BlishLab/scriabin) |
| **FlowSig** | Python | 86 | Causal flow inference on top of CCC outputs (post-analysis) | [GitHub](https://github.com/axelalmet/flowsig) |
| **cell2cell** | Python | 79 | Tensor decomposition across contexts (time/tissue/disease) | [GitHub](https://github.com/earmingol/cell2cell) |
| **SpaTalk** | R | 76 | Spatial + knowledge graph LR→target pathway cascade | [GitHub](https://github.com/ZJUFanLab/SpaTalk) |
| **NICHES** | R | 58 | Single-cell resolution niche interactions (cell-pair objects) | [GitHub](https://github.com/msraredon/NICHES) |
| **CrossTalkeR** | R/Python | 49 | Differential CCC network visualization + centrality | [GitHub](https://github.com/CostaLab/CrossTalkeR) |
| **NeuronChat** | R | 45 | Neural-specific LR database (synaptic, gap junction, neuromodulator) | [GitHub](https://github.com/Wei-BioMath/NeuronChat) |
| **scDiffCom** | R | 25 | Differential CCC with built-in 5K LR database, aging atlas | [GitHub](https://github.com/CyrilLagger/scDiffCom) |

## Universal CCC Principles

### LR Database Selection
- **consensus** (LIANA+ default): union of multiple databases — broadest coverage, recommended for discovery
- **CellPhoneDB**: manually curated, strong on heteromeric complexes — best for confidence
- **CellChatDB**: pathway-annotated — best when you need pathway-level grouping
- **Omnipath** (Squidpy default): comprehensive multi-source integration
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

Read the relevant guide(s) based on the decision tree above, then generate code following the templates. For tools without guides, refer to the official documentation linked in the table.
