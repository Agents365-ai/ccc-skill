# LIANA+ Guide

## When to Use

- You want a **unified multi-method meta-analysis** (run CellPhoneDB + CellChat + NATMI + others in one call)
- You need **spatial bivariate metrics** (local/global) for spot-based or single-cell spatial data
- You need **multi-sample comparison** via tensor decomposition or MOFA+
- You need **MISTy multi-view learning** for spatial data
- Your pipeline is Python/AnnData-based

## Installation

```bash
pip install liana
```

## Workflow A: Basic scRNA-seq CCC (rank_aggregate)

This is the recommended default — runs multiple methods and aggregates rankings.

```python
import scanpy as sc
import liana as li

# Load and preprocess
adata = sc.read_h5ad("data.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Run rank_aggregate (meta-method)
li.mt.rank_aggregate(
    adata,
    groupby='cell_type',           # adata.obs column with cell type labels
    resource_name='consensus',     # LR database (use li.rs.show_resources() for options)
    expr_prop=0.1,                 # min fraction of cells expressing gene
    n_perms=1000,                  # permutations for statistical methods
    use_raw=True,                  # read from adata.raw if available
    inplace=True,                  # store in adata.uns['liana_res']
    verbose=True,
)

# Access results
liana_res = adata.uns['liana_res']
# Columns: source, target, ligand_complex, receptor_complex,
#   magnitude_rank, specificity_rank (both 0-1, lower = better)
```

### Run a single method instead

```python
# Available: cellphonedb, cellchat, natmi, connectome, logfc, singlecellsignalr, geometric_mean
li.mt.cellphonedb(adata, groupby='cell_type', n_perms=1000, inplace=True)
# Output: lr_means (magnitude, higher=stronger), cellphone_pvals (specificity, lower=better)
```

### Available methods and their scores

| Method | Magnitude Score (higher=stronger) | Specificity Score (lower=better*) |
|--------|----------------------------------|----------------------------------|
| `cellphonedb` | `lr_means` | `cellphone_pvals` |
| `geometric_mean` | `lr_gmeans` | `gmean_pvals` |
| `cellchat` | `lr_probs` | `cellchat_pvals` |
| `natmi` | `expr_prod` | `spec_weight` (higher=better) |
| `connectome` | `expr_prod` | `scaled_weight` (higher=better) |
| `logfc` | — | `lr_logfc` (higher=better) |
| `singlecellsignalr` | `LRscore` | — |

*For natmi, connectome, logfc: higher specificity score = more specific (opposite convention).

### Per-sample inference

```python
# Run CCC per sample
li.mt.rank_aggregate.by_sample(
    adata,
    groupby='cell_type',
    sample_key='sample',           # adata.obs column for sample identity
    inplace=True,
)
```

## Workflow B: Spatial Bivariate Metrics

For spatially-resolved data (Visium, MERFISH, Slide-seq, etc.).

```python
import scanpy as sc
import liana as li

adata = sc.read_h5ad("spatial.h5ad")  # needs adata.obsm['spatial']

# Step 1: Build spatial connectivity graph
li.ut.spatial_neighbors(
    adata,
    bandwidth=150,                 # kernel bandwidth (in coordinate units)
    cutoff=0.1,                    # weight cutoff (below this → 0)
    kernel='gaussian',
    set_diag=False,
    spatial_key='spatial',
    key_added='spatial',           # stores in adata.obsp['spatial_connectivities']
)

# Step 2: Run bivariate analysis
lr_adata = li.mt.bivariate(
    adata,
    local_name='cosine',           # local metric: cosine, pearson, spearman, jaccard, product, morans
    global_name='morans',          # global metric: morans, lee (or None to skip)
    resource_name='consensus',
    connectivity_key='spatial_connectivities',
    nz_prop=0.05,                  # min nonzero proportion to keep an interaction
    add_categories=True,           # assign interaction categories per spot
    n_perms=100,                   # permutations for global p-values
)

# lr_adata is an AnnData:
#   .X              → local bivariate scores (n_spots × n_interactions)
#   .var             → global stats (Moran's I, p-values) per interaction
#   .layers['cats']  → category assignments
#   .layers['pvals'] → local p-values (if applicable)
```

### Available local metrics

| Metric | Best for | Notes |
|--------|----------|-------|
| `cosine` | General spatial co-expression | Default, robust |
| `pearson` | Linear spatial correlation | Sensitive to outliers |
| `spearman` | Rank-based spatial correlation | More robust |
| `product` | Weighted product (Lee's L when z-scaled) | Simple, fast |
| `morans` | Moran's R (spatialDM-style) | Classic spatial stat |
| `masked_spearman` | scHOT-style masked correlation | Handles sparse data |

## Workflow C: Spatial Proximity-Weighted CCC

New in v1.7.0 — add spatial proximity weights to standard CCC methods.

```python
li.mt.rank_aggregate(
    adata,
    groupby='cell_type',
    spatial_key='spatial',                  # enables spatial weighting
    spatial_kwargs={'bandwidth': 150},      # kernel parameters
    inplace=True,
)
```

## Workflow D: Inflow Analysis

Trivariate spatial metric — quantifies how much signal flows into each spot from each cell type.

```python
li.ut.spatial_neighbors(adata, bandwidth=150)

inflow_adata = li.mt.inflow(
    adata,
    groupby='cell_type',
    connectivity_key='spatial_connectivities',
    resource_name='consensus',
    nz_prop=0.001,
)
# Shape: (n_spots, n_cell_types × n_interactions)
```

## Workflow E: MISTy Multi-View Learning

Identifies what predicts target gene expression from different spatial contexts.

```python
# Build views
misty = li.mt.genericMistyData(
    intra=adata,                     # intra-cellular view (expression at each spot)
    spatial_key='spatial',
    bandwidth=100,                   # paraview kernel bandwidth
    add_para=True,                   # add paracrine view
    add_juxta=False,                 # add juxtacrine view
)

# Run model
misty(
    model=li.mt.LinearModel(seed=42),    # or RobustLinearModel, RandomForestModel
    bypass_intra=False,
    k_cv=10,
)

# Results
# misty.uns['target_metrics']  → R², gain per target
# misty.uns['interactions']    → feature importance per view

# Plot
li.pl.target_metrics(misty)
li.pl.contributions(misty)
li.pl.interactions(misty)
```

### LR-specific MISTy

```python
misty = li.mt.lrMistyData(
    adata,
    resource_name='consensus',
    spatial_key='spatial',
    bandwidth=100,
)
misty(model=li.mt.LinearModel(seed=42))
```

## Workflow F: Multi-Sample Tensor Decomposition

```python
# Step 1: Run CCC per sample
li.mt.rank_aggregate(adata, groupby='cell_type', inplace=True)

# Step 2: Build tensor
tensor = li.mu.to_tensor_c2c(
    adata,
    sample_key='sample',
    score_key='specificity_rank',
    non_negative=True,
)
# Shape: (n_samples, n_senders, n_receivers, n_interactions)

# Then decompose with cell2cell or tensorly
```

## Workflow G: Multi-Sample MOFA+

```python
# Convert LR results to MOFA views
mdata = li.mu.lrs_to_views(
    adata,
    score_key='specificity_rank',
    sample_key='sample',
    lr_prop=0.5,
    lrs_per_view=20,
)

# Then run MOFA via mofapy2 or muon
```

## Visualization

### Dotplot

```python
li.pl.dotplot(
    adata,
    colour='specificity_rank',
    size='magnitude_rank',
    inverse_colour=True,           # lower rank = darker color
    inverse_size=True,             # lower rank = bigger dot
    source_labels=['TypeA', 'TypeB'],
    target_labels=['TypeC'],
    top_n=20,
    orderby='specificity_rank',
    orderby_ascending=True,
    cmap='viridis',
    figure_size=(8, 6),
)
```

### Dotplot by sample

```python
li.pl.dotplot_by_sample(
    adata,
    sample_key='sample',
    colour='specificity_rank',
    size='magnitude_rank',
    inverse_colour=True,
    inverse_size=True,
)
```

### Tileplot

```python
li.pl.tileplot(
    adata,
    fill='magnitude_rank',
    top_n=20,
    orderby='magnitude_rank',
    orderby_ascending=True,
)
```

### Circle plot

```python
li.pl.circle_plot(
    adata,
    groupby='cell_type',
    score_key='specificity_rank',
    inverse_score=True,
    top_n=30,
)
```

## Key Parameters

| Parameter | Default | When to change |
|-----------|---------|----------------|
| `resource_name` | `'consensus'` | Use specific DB if needed: `li.rs.show_resources()` |
| `expr_prop` | `0.1` | Lower (0.05) for rare cell types; raise (0.2) to reduce noise |
| `n_perms` | `1000` | Use 100 for quick testing; 1000+ for publication |
| `use_raw` | `True` | Set `False` if no `adata.raw`; or use `layer='...'` |
| `bandwidth` (spatial) | varies | Match to technology: ~150 for Visium, ~50 for MERFISH |
| `nz_prop` (bivariate) | `0.05` | Lower for sparse data; raise to reduce noise |

## Output Format

All results stored in `adata.uns['liana_res']` (DataFrame):

| Column | Meaning |
|--------|---------|
| `source` | Sending cell type |
| `target` | Receiving cell type |
| `ligand_complex` | Ligand name (complexes joined with `_`) |
| `receptor_complex` | Receptor name |
| `magnitude_rank` | Aggregate magnitude rank (0-1, lower=stronger) |
| `specificity_rank` | Aggregate specificity rank (0-1, lower=more specific) |
| Per-method columns | Varies by method (see table above) |

## Resource (LR Database) Management

```python
# List available resources
li.rs.show_resources()

# Get a specific resource as DataFrame
resource = li.rs.select_resource('consensus')  # columns: ligand, receptor

# Use a custom resource
custom_lr = pd.DataFrame({'ligand': [...], 'receptor': [...]})
li.mt.rank_aggregate(adata, groupby='cell_type', resource=custom_lr)

# Translate to another species
translated = li.rs.translate_resource(resource, target_organism='mouse')
```

## Pitfalls

1. **`use_raw=True` is the default** — LIANA reads from `adata.raw` if it exists. If you normalized `adata.X` but didn't set `adata.raw`, set `use_raw=False` or `layer='...'`.

2. **Gene coverage check** — LIANA asserts that at least 2% of resource genes are present. If your gene names don't match (e.g., Ensembl IDs vs symbols), you'll get an assertion error. Convert gene names first.

3. **Spatial prerequisite** — `adata.obsp['spatial_connectivities']` MUST exist before running `bivariate`, `inflow`, or spatial-weighted methods. Use `li.ut.spatial_neighbors()` or `squidpy.gr.spatial_neighbors()`.

4. **LR separator convention** — LIANA uses `^` to join ligand-receptor pairs (e.g., `TGFB1^TGFBR1`) and `_` for complex subunits. Don't confuse with other tools' separators.

5. **`n_perms=0`** — disables permutation testing entirely. Only Moran's I has an analytical p-value fallback; other methods return no p-values.

6. **MuData support** — most methods accept `MuData` for multi-modal. Pass `mdata_kwargs` to specify modality mappings.

7. **rank_aggregate scores are ranks, not raw scores** — `magnitude_rank=0.01` means top 1%, not a raw probability. Don't threshold on absolute values as if they were p-values.
