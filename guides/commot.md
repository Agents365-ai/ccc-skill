# COMMOT Guide

## When to Use

- You need **spatial CCC with optimal transport** — finds the optimal cell-to-cell signaling plan
- You want **signaling direction vector fields** (which direction is the signal flowing?)
- Your data is **spatially resolved** (Visium, MERFISH, Slide-seq, Xenium, etc.)
- You want to identify **communication-dependent genes** (genes whose expression correlates with signaling)
- You want distance-aware CCC that respects spatial organization

## Installation

```bash
pip install commot
# With tradeSeq integration for downstream gene detection:
pip install commot[tradeSeq]
```

## Core Workflow

### Step 1: Load and preprocess

```python
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("spatial.h5ad")
# Requires: adata.obsm['spatial'] (n_cells × 2 coordinates)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
# Do NOT sc.pp.scale() — COMMOT requires non-negative values
```

### Step 2: Load LR database

```python
df_ligrec = ct.pp.ligand_receptor_database(
    database='CellChat',            # 'CellChat' or 'CellPhoneDB_v4.0'
    species='mouse',                # 'mouse' or 'human'
    signaling_type='Secreted Signaling',  # or 'Cell-Cell Contact', 'ECM-Receptor', None for all
)
# Returns DataFrame: [ligand, receptor, pathway, signaling_type]

# Filter to expressed genes (recommended for performance)
df_ligrec = ct.pp.filter_lr_database(
    df_ligrec, adata,
    heteromeric=True,
    heteromeric_delimiter='_',
    heteromeric_rule='min',         # 'min' or 'ave' for complex subunit expression
    min_cell_pct=0.05,              # keep pairs expressed in >=5% of cells
)
```

### Step 3: Run spatial communication inference

```python
ct.tl.spatial_communication(
    adata,
    database_name='cellchat',       # key prefix for results
    df_ligrec=df_ligrec,
    pathway_sum=True,               # compute pathway-level totals
    heteromeric=True,               # handle multi-subunit complexes
    heteromeric_rule='min',
    heteromeric_delimiter='_',
    dis_thr=500,                    # distance cutoff (same units as coordinates)
    cost_type='euc',                # 'euc' or 'euc_square'
    cot_eps_p=1e-1,                 # entropy regularization
    cot_rho=1e1,                    # unmatched mass penalty
    cot_nitermax=10000,
)
```

### Step 4: Compute signaling directions

```python
ct.tl.communication_direction(
    adata,
    database_name='cellchat',
    lr_pair=('Fgf1', 'Fgfr1'),     # specific LR pair
    # pathway_name='TGFb',          # OR pathway-level
    k=5,                            # top k senders/receivers per cell
)
# Stores: adata.obsm['commot_sender_vf-cellchat-Fgf1-Fgfr1']  (n_cells × 2)
# Stores: adata.obsm['commot_receiver_vf-cellchat-Fgf1-Fgfr1']
```

### Step 5: Cluster-level analysis

```python
# Label permutation (fast)
ct.tl.cluster_communication(
    adata,
    database_name='cellchat',
    lr_pair=('Fgf1', 'Fgfr1'),
    clustering='leiden',
    n_permutations=100,
    random_seed=1,
)
# Stores: adata.uns['commot_cluster-leiden-cellchat-Fgf1-Fgfr1']
#   {'communication_matrix': DataFrame, 'communication_pvalue': DataFrame}

# Spatial permutation (more accurate for neighboring clusters)
ct.tl.cluster_communication_spatial_permutation(
    adata,
    df_ligrec=df_ligrec,
    database_name='cellchat',
    clustering='leiden',
    n_permutations=100,
)

# Compute cluster positions (needed for network plot)
ct.tl.cluster_position(adata, clustering='leiden')
```

## Output Format

```python
# Cell-by-cell signaling matrices (sparse)
adata.obsp['commot-cellchat-Fgf1-Fgfr1']      # P[i,j] = signal from cell i to cell j
adata.obsp['commot-cellchat-TGFb']             # pathway-level (if pathway_sum=True)
adata.obsp['commot-cellchat-total-total']       # total across all pairs

# Per-cell summary scores
adata.obsm['commot-cellchat-sum-sender']        # DataFrame: 's-Fgf1-Fgfr1', 's-total-total', ...
adata.obsm['commot-cellchat-sum-receiver']      # DataFrame: 'r-Fgf1-Fgfr1', 'r-total-total', ...

# Vector fields (after communication_direction)
adata.obsm['commot_sender_vf-cellchat-Fgf1-Fgfr1']    # (n_cells × 2)
adata.obsm['commot_receiver_vf-cellchat-Fgf1-Fgfr1']

# Metadata
adata.uns['commot-cellchat-info']               # {'df_ligrec': ..., 'distance_threshold': ...}
```

### Key naming convention

| Slot | Pattern | Content |
|------|---------|---------|
| `obsp` | `commot-{db}-{lig}-{rec}` | Cell×cell sparse signaling matrix |
| `obsp` | `commot-{db}-{pathway}` | Pathway-summed matrix |
| `obsp` | `commot-{db}-total-total` | All-pairs aggregate |
| `obsm` | `commot-{db}-sum-sender` | Per-cell sent signal (DataFrame) |
| `obsm` | `commot-{db}-sum-receiver` | Per-cell received signal (DataFrame) |
| `obsm` | `commot_sender_vf-{db}-{lig}-{rec}` | Sender vector field |
| `obsm` | `commot_receiver_vf-{db}-{lig}-{rec}` | Receiver vector field |

## Visualization

### Cell-level communication with directions

```python
ct.pl.plot_cell_communication(
    adata,
    database_name='cellchat',
    lr_pair=('Fgf1', 'Fgfr1'),
    plot_method='grid',             # 'cell', 'grid', 'stream'
    background='summary',           # 'summary', 'image', 'cluster'
    summary='sender',               # 'sender' or 'receiver'
    cmap='coolwarm',
    scale=0.00003,                  # arrow scale (smaller = longer arrows)
    grid_density=0.4,
    filename='ccc_direction.pdf',
)
```

### Cluster-level network

```python
ct.pl.plot_cluster_communication_network(
    adata,
    uns_names=['commot_cluster-leiden-cellchat-Fgf1-Fgfr1'],
    clustering='leiden',
    quantile_cutoff=0.99,
    p_value_cutoff=0.05,
    nx_node_pos='cluster',
    filename='cluster_network.pdf',
)
```

### Dot plot

```python
ct.pl.plot_cluster_communication_dotplot(
    adata,
    database_name='cellchat',
    pathway_name='TGFb',
    clustering='leiden',
)
```

## Downstream Analysis

### Communication-dependent genes (requires tradeSeq/R)

```python
# Requires: adata.layers['counts'] (raw), adata.raw (full normalized)
df_deg, df_yhat = ct.tl.communication_deg_detection(
    adata,
    database_name='cellchat',
    lr_pair=('total', 'total'),
    summary='receiver',
    n_var_genes=3000,
    nknots=6,
    deg_pvalue_cutoff=0.05,
)

df_deg_clus, df_yhat_clus = ct.tl.communication_deg_clustering(
    df_deg, df_yhat,
    deg_clustering_npc=10,
    n_deg_genes=200,
)

ct.pl.plot_communication_dependent_genes(df_deg_clus, df_yhat_clus, filename='deg.pdf')
```

### Communication impact (partial correlation)

```python
# Requires: adata.raw (full normalized)
df_impact = ct.tl.communication_impact(
    adata,
    database_name='cellchat',
    pathway_name='TGFb',
    method='partial_corr',          # 'partial_corr', 'semipartial_corr', 'treebased_score'
    corr_method='spearman',
    ds_genes=['Gene1', 'Gene2'],
    bg_genes=100,
)
ct.pl.plot_communication_impact(df_impact, summary='receiver', filename='impact.pdf')
```

### Group CCC patterns

```python
# Cluster LR pairs by similar CCC network
comm_id, D = ct.tl.group_cluster_communication(
    adata,
    clustering='leiden',
    keys=['cellchat-Fgf1-Fgfr1', 'cellchat-TGFb'],
    dissimilarity_method='jaccard',
)
```

## Non-Spatial Data: Infer Spatial Information

```python
# Map scRNA-seq to spatial locations via optimal transport
adata_sc_pred, adata_sp_pred, gamma = ct.pp.infer_spatial_information(
    adata_sc,                 # scRNA-seq AnnData
    adata_sp,                 # spatial AnnData with obsm['spatial']
    ot_alpha=0.2,
    ot_rho=0.05,
    ot_epsilon=0.01,
    loc_pred_k=1,
)
# adata_sc_pred.obsm['spatial'] = predicted locations for scRNA-seq cells
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `dis_thr` | — | **Most important.** Match to technology: ~500 for Visium, ~50-150 for MERFISH/Xenium. Units = `adata.obsm['spatial']` units. |
| `cot_rho` | 10 | Unmatched mass penalty. Higher = stricter mass balance. Default is robust. |
| `cot_eps_p` | 0.1 | Entropy regularization. Smaller = sparser transport plan. |
| `cot_nitermax` | 10000 | Increase if results look noisy. |
| `heteromeric` | True | **Required for CellChatDB** (complexes encoded as `Gene1_Gene2`). |
| `pathway_sum` | True | Enable pathway-level analysis. |

### `dis_thr` by technology

| Technology | Typical `dis_thr` | Notes |
|------------|-------------------|-------|
| 10x Visium | 500 | ~5 spot diameters |
| MERFISH | 50–100 | Single-cell resolution |
| Slide-seq | 100–200 | Near single-cell |
| 10x Xenium | 50–100 | Single-cell |
| STARmap | 50–100 | Single-cell |

## Pitfalls

1. **Do NOT `sc.pp.scale()` before COMMOT.** Creates negative values; COMMOT requires non-negative abundances.

2. **`dis_thr` units must match `adata.obsm['spatial']` units.** Visium `obsm['spatial']` may be in pixels, not microns. Convert or calibrate accordingly.

3. **`heteromeric=True` is required for CellChatDB.** Without it, most CellChat entries (with `Gene1_Gene2` complexes) will be missed.

4. **Filter the LR database before large runs.** Full CellChatDB (~2000 pairs) without filtering is very slow. Always use `ct.pp.filter_lr_database()`.

5. **`communication_deg_detection` requires raw counts** in `adata.layers['counts']`. Save counts before normalization.

6. **`communication_impact` requires `adata.raw`.** Save full normalized expression with `adata.raw = adata.copy()` before gene subsetting.

7. **`cluster_position()` must run before `plot_cluster_communication_network`** with `nx_node_pos='cluster'`.

8. **`pathway_sum=True` is needed** for any pathway-level downstream calls. Without it, `adata.obsp['commot-db-PathwayX']` won't exist.

9. **Memory for large datasets.** Each LR pair creates an n×n sparse matrix. For datasets >10K cells with >100 LR pairs, pre-filter aggressively.

10. **OT convergence.** If results are noisy or uniform, try increasing `cot_nitermax` or adjusting `cot_rho`. For very sparse data, try `smooth=True`.
