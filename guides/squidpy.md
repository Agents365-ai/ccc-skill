# Squidpy CCC Guide

## When to Use

- Your pipeline is already **scanpy/scverse-based** and you want CCC with zero extra dependencies
- You want the **OmniPath database** (superset of CellPhoneDB, HPRD, and many others) out of the box
- You need a **fast permutation test** (Numba JIT) for large datasets
- You want results that stay in `adata.uns` and persist through `write_h5ad()`

**Important**: Squidpy's `ligrec` is a **cluster-level** permutation test (CellPhoneDB algorithm). It does NOT use the spatial graph — it tests all cluster pairs globally. For spatially-restricted CCC, use COMMOT or LIANA+ bivariate.

## Installation

```bash
pip install squidpy
# OmniPath is pulled automatically as a dependency
```

## Core Workflow

### Step 1: Preprocess with scanpy

```python
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad("data.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# CRITICAL: save raw BEFORE any further transforms
adata.raw = adata

sc.pp.highly_variable_genes(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="celltype")

# Ensure categorical
adata.obs["celltype"] = adata.obs["celltype"].astype("category")
```

### Step 2: Run LR analysis

```python
sq.gr.ligrec(
    adata,
    cluster_key="celltype",
    n_perms=1000,
    use_raw=True,                # reads from adata.raw
    threshold=0.01,              # min 1% cells expressing gene
    seed=42,
)
# Results stored in adata.uns["celltype_ligrec"]
```

### Customize the database

```python
# Use only CellPhoneDB interactions
sq.gr.ligrec(
    adata, cluster_key="celltype", n_perms=1000,
    interactions_params={"resources": "CellPhoneDB"},
)

# Restrict to ligand→receptor direction
sq.gr.ligrec(
    adata, cluster_key="celltype", n_perms=1000,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
)

# Custom LR pairs
import pandas as pd
custom_lr = pd.DataFrame({
    "source": ["TGFB1", "VEGFA", "CXCL12"],
    "target": ["TGFBR1", "KDR", "CXCR4"],
})
sq.gr.ligrec(adata, cluster_key="celltype", n_perms=1000, interactions=custom_lr)
```

### FDR correction

```python
sq.gr.ligrec(
    adata, cluster_key="celltype", n_perms=1000,
    corr_method="fdr_bh",
    corr_axis="interactions",   # correct across cluster pairs per interaction
)
```

### Restrict to specific cluster pairs

```python
sq.gr.ligrec(
    adata, cluster_key="celltype", n_perms=1000,
    clusters=[("Macrophages", "T cells"), ("Fibroblasts", "Endothelial")],
)
```

## Output Format

Results in `adata.uns["{cluster_key}_ligrec"]` — a dict with 3 DataFrames:

```python
res = adata.uns["celltype_ligrec"]

res["means"]      # MultiIndex rows: (ligand, receptor)
                  # MultiIndex cols: (cluster_A, cluster_B)
                  # Values: mean expression (ligand_in_A + receptor_in_B) / 2
                  # NaN = not tested (below threshold)

res["pvalues"]    # Same shape as means
                  # Values: permutation p-value (or FDR q-value if corrected)
                  # NaN = not tested

res["metadata"]   # OmniPath metadata columns (databases, is_directed, etc.)
```

### Accessing results

```python
# Top significant interactions for a cluster pair
pair = ("Macrophages", "T cells")
pvals = res["pvalues"][pair].dropna().sort_values()
sig = pvals[pvals < 0.05]
print(sig.index.tolist()[:10])  # list of (ligand, receptor) tuples

# Convert sparse to dense for easier manipulation
means_dense = res["means"].sparse.to_dense()
pvals_dense = res["pvalues"].sparse.to_dense()
```

## Visualization

### Dot plot

```python
# Basic
sq.pl.ligrec(adata, cluster_key="celltype")

# Filtered
sq.pl.ligrec(
    adata,
    cluster_key="celltype",
    source_groups=["Macrophages", "Fibroblasts"],
    target_groups=["T cells", "NK cells"],
    pvalue_threshold=0.05,
    means_range=(0.5, float("inf")),
    remove_empty_interactions=True,
    remove_nonsig_interactions=True,
    dendrogram="both",
    alpha=0.001,                # significance marker threshold
    figsize=(20, 15),
    save="ligrec_plot.pdf",
)

# Swap axes (cluster pairs as rows)
sq.pl.ligrec(adata, cluster_key="celltype", swap_axes=True)
```

**Plot encoding:**
- Dot color = mean expression
- Dot size = -log10(p-value), larger = more significant
- Hollow dot = significant at `alpha` threshold

## Spatial Integration (complementary analyses)

While `ligrec` itself is not spatial, combine it with Squidpy's spatial statistics:

```python
# Build spatial graph
sq.gr.spatial_neighbors(adata, coord_type="visium", n_rings=1)

# Neighborhood enrichment — which cell types co-locate?
sq.gr.nhood_enrichment(adata, cluster_key="celltype")
sq.pl.nhood_enrichment(adata, cluster_key="celltype")

# Co-occurrence — spatial co-occurrence probability
sq.gr.co_occurrence(adata, cluster_key="celltype")
sq.pl.co_occurrence(adata, cluster_key="celltype")

# Spatial autocorrelation
sq.gr.spatial_autocorr(adata, mode="moran")

# Strategy: use nhood_enrichment to find spatially co-occurring cell types,
# then focus ligrec analysis on those pairs
```

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `n_perms` | `1000` | Reduce to 100 for testing; 1000+ for publication |
| `threshold` | `0.01` | Min fraction of cells expressing gene; lower → more interactions tested |
| `use_raw` | `True` | Reads `adata.raw.X`; set False if no `.raw` |
| `complex_policy` | `"min"` | How to resolve complexes: `"min"` (lowest member) or `"all"` (all combos) |
| `corr_method` | `None` | FDR method: `"fdr_bh"`, `"bonferroni"`, etc. |
| `corr_axis` | `"clusters"` | Correction axis: `"clusters"` or `"interactions"` |
| `seed` | `None` | Set for reproducibility |
| `n_jobs` | `1` | Parallel CPU cores |

## Pitfalls

1. **Gene name case mismatch (most common).** OmniPath returns **uppercase** gene symbols. If your `adata.var_names` are lowercase or mixed case, zero interactions survive. Fix: `adata.var_names = adata.var_names.str.upper()`

2. **Missing `.raw` attribute.** `use_raw=True` is default. If you didn't save `adata.raw = adata` before scaling, you get an error. Either save raw early or use `use_raw=False`.

3. **`ligrec` does NOT use the spatial graph.** This is a common misconception. It's a cluster-level permutation test, not spatially restricted. Building `sq.gr.spatial_neighbors()` does not affect `ligrec` results.

4. **Non-categorical cluster column.** `cluster_key` column must be `pd.Categorical`. Convert with `adata.obs["col"] = adata.obs["col"].astype("category")`.

5. **Ensembl IDs instead of gene symbols.** OmniPath uses HGNC symbols. If your var_names are Ensembl IDs, use the `gene_symbols` parameter pointing to a var column with symbols.

6. **All p-values NaN.** Usually means `threshold` is too high or gene names don't match. Try lowering `threshold` or check gene name format.

7. **OmniPath requires internet.** The database is fetched at runtime (then cached). First run needs network access.

8. **Memory on large datasets.** With many clusters (>30) and default OmniPath (~40K interactions), the permutation matrix is large. Filter interactions or restrict cluster pairs.
