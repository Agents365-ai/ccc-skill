# DIALOGUE Guide

## When to Use

- You want to find **multicellular programs (MCPs)** — coordinated gene expression programs across multiple cell types in shared tissue microenvironments
- You want to understand **how cell types co-regulate** each other at the transcriptional level (not just LR pairs)
- You have **multiple samples/patients** and want to find programs that vary across niches
- You want to **associate MCPs with clinical phenotypes** (e.g., responder vs non-responder)
- You have **spatial data** and want to find spatially co-regulated programs across cell types

**Key difference from LR tools**: DIALOGUE does NOT identify sender→receiver signaling. It finds programs that are environmentally co-regulated — both cell types respond to (or co-create) a shared tissue state. This is complementary to LR analysis, not a replacement.

| Dimension | LR Analysis (CellChat, NicheNet) | DIALOGUE MCPs |
|-----------|----------------------------------|---------------|
| Unit | Individual ligand-receptor pairs | Whole gene programs spanning many genes |
| Evidence | Co-expression of L in sender + R in receiver | Cross-sample correlation of cell-type gene modules |
| Direction | Directional (sender → receiver) | Undirected co-regulation |
| Requires | Cell type labels | Cell type labels + sample/niche grouping |

## Installation

```r
devtools::install_github("livnatje/DIALOGUE")

# Dependencies (install if missing)
install.packages(c("lme4", "lmerTest", "PMA", "plyr", "matrixStats",
                   "psych", "nnls", "RColorBrewer", "reshape2", "ggplot2"))

library(DIALOGUE)
```

## Core Workflow

### Step 1: Build cell.type objects

**From Seurat (recommended):**

```r
library(DIALOGUE)
library(Seurat)

# Add required metadata
seurat_obj$cell.subtypes = Idents(seurat_obj)
seurat_obj$samples = seurat_obj$patient_id   # sample/patient/donor IDs

r1 = DIALOGUE_make.cell.type.seurat(
    obj = seurat_obj,
    cell.subtypes = c("CD14+ Mono", "FCGR3A+ Mono", "DC"),
    name = "Myeloid"
)

r2 = DIALOGUE_make.cell.type.seurat(
    obj = seurat_obj,
    cell.subtypes = c("CD4 T", "CD8 T", "NK"),
    name = "Lymphoid"
)

rA = list(Myeloid = r1, Lymphoid = r2)
```

**From scratch (more control):**

```r
# tpm must be genes x cells (log-normalized)
tpm_mat = as.matrix(GetAssayData(sub_obj, slot = "data"))
pca_mat = Embeddings(sub_obj, reduction = "pca")     # cells x PCs
sample_vec = sub_obj$patient_id
cellQ_vec = as.numeric(scale(log1p(Matrix::colSums(GetAssayData(sub_obj, slot = "counts")))))
names(cellQ_vec) = colnames(sub_obj)

meta_df = data.frame(cellQ = cellQ_vec, row.names = colnames(sub_obj))

r1 = make.cell.type(
    name = "Myeloid",
    tpm = tpm_mat,
    samples = sample_vec,
    X = pca_mat,
    metadata = meta_df,
    cellQ = cellQ_vec
)
```

### Step 2: Set parameters

```r
dir.create("DLG.results/", recursive = TRUE)  # MUST exist before running

param = DLG.get.param(
    k = 3,                          # number of MCPs to identify (start with 2-3)
    results.dir = "DLG.results/",
    conf = "cellQ",                 # confounders in HLM (vector OK: c("cellQ","batch"))
    covar = c("cellQ", "tme.qc"),   # covariates in HLM formula
    pheno = "response",             # binary phenotype column (optional)
    n.genes = 200,                  # genes per component in PMD step
    abn.c = 15,                     # min cells per sample to count that sample
    p.anova = 0.05,                 # ANOVA filter for PCs
)
```

### Step 3: Run DIALOGUE

```r
R = DIALOGUE.run(
    rA = rA,
    main = "MyStudy",               # prefix for output files
    param = param,
    plot.flag = TRUE,                # auto-generate PDF
)
```

This runs three internal steps:
1. **DIALOGUE1** — MultiCCA (penalized matrix decomposition) on sample-averaged PCs → find joint components
2. **DIALOGUE2** — Hierarchical linear modeling per gene across cell-type pairs → validate gene programs
3. **DIALOGUE3** — NNLS-based scoring → final per-cell MCP activities and gene signatures

## Output Format

### `R$MCPs` — Gene signatures per MCP

```r
# Named list of MCPs, each with cell-type-specific up/down gene sets
R$MCPs$MCP1$Myeloid.up      # character vector of upregulated genes
R$MCPs$MCP1$Myeloid.down
R$MCPs$MCP1$Lymphoid.up
R$MCPs$MCP1$Lymphoid.down
```

### `R$scores` — Per-cell MCP activity

```r
# Named list, one data.frame per cell type
# Columns: MCP1, MCP2, ..., samples, cells, cell.type
head(R$scores$Myeloid)
```

### `R$gene.pval` — Gene-level statistics

```r
# Named list, one data.frame per cell type
# Key columns: genes, program, up, p.up, p.down, n.up, Nf, coef
head(R$gene.pval$Myeloid[order(R$gene.pval$Myeloid$p.up), ])
```

### `R$pref` — Cross-cell-type correlations

```r
# Pearson correlation of sample-averaged MCP scores between cell types
R$pref$Myeloid.vs.Lymphoid
```

### `R$phenoZ` — Phenotype associations

```r
# Matrix: rows = cell types + "All", cols = MCPs
# Values = direction x -log10(p) from HLM
R$phenoZ
```

### `R$emp.p` — Empirical p-values

```r
# Matrix: rows = MCPs, cols = cell-type pairs
# From 100 permutations of sample labels
R$emp.p
```

## Adding MCP scores back to Seurat

```r
for (ct in names(R$scores)) {
    sc = R$scores[[ct]]
    mcp_cols = grep("^MCP", colnames(sc), value = TRUE)
    for (col in mcp_cols) {
        scores_named = setNames(sc[[col]], sc$cells)
        seurat_obj = AddMetaData(seurat_obj, scores_named, col.name = paste0(ct, "_", col))
    }
}
FeaturePlot(seurat_obj, features = "Myeloid_MCP1")
```

## Spatial Mode

```r
# For spatial data, each spot = one "sample"
param_spatial = DLG.get.param(
    k = 3,
    results.dir = "spatial_results/",
    spatial.flag = TRUE,            # skips ANOVA filter on PCs
    conf = "cellQ",
    covar = c("cellQ", "tme.qc"),
)
# The 'samples' slot in each cell.type should contain spot/niche IDs
```

## Visualization

```r
# Full PDF output
DIALOGUE.plot(R, results.dir = "DLG.results/", pheno = "response")

# Sample-level MCP correlation scatter
DIALOGUE.plot.av(R, MCPs = "MCP1")

# Gene set composition bar plot
DIALOGUE.plot.sig.comp(R, main = "MyStudy")

# Violin plots by phenotype
DIALOGUE.violin.pheno(R, pheno = "response", MCPs = c("MCP1"))

# Pairwise correlation with error bars
DLG.cor.plot(r1, r2, idx = "MCP1")
```

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `k` | 2 | Number of MCPs. Start small (2-3), check `R$emp.p` for significance |
| `abn.c` | 15 | Min cells per sample; lower for small datasets (careful: more noise) |
| `p.anova` | 0.05 | Raise to 0.1 if too many PCs filtered out |
| `n.genes` | 200 | Increase to 500 for larger datasets |
| `conf` | `"cellQ"` | Add batch, sex, etc. as confounders |
| `pheno` | `NULL` | Binary (0/1) metadata column for phenotype association |
| `spatial.flag` | `FALSE` | `TRUE` for spatial data (spots as samples) |
| `specific.pair` | `NULL` | Focus on one cell-type pair: `c("TypeA", "TypeB")` |
| `bypass.emp` | `FALSE` | Skip permutation test (only for debugging) |

## Pitfalls

1. **Too few samples.** DIALOGUE requires >=5 samples with sufficient cells. With <15-20 samples, results are unstable. This is NOT a single-sample tool.

2. **`tpm` matrix must be genes x cells.** Transposed matrix produces silent wrong results. Seurat's `GetAssayData(slot="data")` is already genes x cells — that's correct.

3. **Results directory must exist.** `dir.create("DLG.results/", recursive = TRUE)` before running.

4. **`cellQ` is mandatory.** Must be a numeric column in metadata. Use `scale(log1p(colSums(counts)))`.

5. **Sample labels must overlap across cell types.** DIALOGUE intersects samples — mismatched labels reduce your effective sample count.

6. **Choosing `k`.** No automatic selection. Check `R$emp.p` — if all MCPs have high p-values for most pairs, reduce `k`.

7. **Cell type fails ANOVA filter.** Error: "Only N features passed the ANOVA filter." Fix: raise `p.anova` to 0.1 or remove that cell type.

8. **Seurat helper uses `scale.factor = 1e5`**, not the standard 1e4. For consistency, normalize manually and use `make.cell.type()` directly.

9. **Confounders must be in metadata.** Any variable in `conf` or `covar` (except auto-computed `tme.qc`) must exist as a column in each cell type's metadata.

10. **`bypass.emp = TRUE` gives unreliable results.** The empirical permutation is what validates MCPs — only skip for debugging.
