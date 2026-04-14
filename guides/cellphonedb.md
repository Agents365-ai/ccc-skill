# CellPhoneDB Guide

## When to Use

- You need **high-confidence curated interactions** with manual curation of heteromeric complexes
- You want **permutation-based statistical testing** for LR significance
- Your data is **human** (CellPhoneDB is human-gene-specific; mouse requires ortholog conversion)
- You want the **DEG-based method** (no permutation needed, bring your own DEGs)
- You want **interaction scoring** (v5 specificity score 0-100)

## Installation

```bash
pip install cellphonedb
```

## Database Setup

```python
from cellphonedb.utils import db_utils

# Download database
db_utils.download_database('./db', 'v5.0.0')
# Creates: ./db/cellphonedb.zip
```

## Three Analysis Methods

### Method 1: Basic Mean Analysis (no stats)

```python
from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
    cpdb_file_path='db/v5/cellphonedb.zip',
    meta_file_path='data/metadata.tsv',
    counts_file_path='data/normalised_log_counts.h5ad',
    counts_data='hgnc_symbol',          # 'ensembl' | 'gene_name' | 'hgnc_symbol'
    output_path='results/',
    score_interactions=True,            # v5 specificity scoring
    threshold=0.1,                      # min fraction of cells expressing gene
    threads=4,
)
# Returns dict: means_result, deconvoluted, deconvoluted_percents, interaction_scores
```

### Method 2: Statistical Analysis (permutation test) — Most common

```python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path='db/v5/cellphonedb.zip',
    meta_file_path='data/metadata.tsv',
    counts_file_path='data/normalised_log_counts.h5ad',
    counts_data='hgnc_symbol',
    output_path='results/',
    score_interactions=True,
    iterations=1000,                    # permutations
    threshold=0.1,
    threads=4,
    debug_seed=42,                      # >=0 for reproducibility (single-thread only)
    pvalue=0.05,
)
# Returns dict: means, pvalues, significant_means, deconvoluted,
#   deconvoluted_percents, interaction_scores
```

### Method 3: DEG-Based Analysis (no permutation)

```python
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results = cpdb_degs_analysis_method.call(
    cpdb_file_path='db/v5/cellphonedb.zip',
    meta_file_path='data/metadata.tsv',
    counts_file_path='data/normalised_log_counts.h5ad',
    degs_file_path='data/DEGs.tsv',     # two-column: cell_type, gene
    counts_data='hgnc_symbol',
    output_path='results/',
    score_interactions=True,
    threshold=0.1,
    threads=4,
)
# Returns dict: means, relevant_interactions, significant_means, deconvoluted, ...
```

## Required Input Formats

### Counts file
- `.h5ad` (recommended), `.h5`, `.txt`/`.tsv`, or 10x folder path
- **Log-normalized** counts. Do NOT z-scale.
- Can pass an in-memory AnnData object directly.

### Metadata file
- Two-column TSV: `barcode_sample`, `cell_type`
- Barcodes must match counts file.

### DEGs file (Method 3 only)
- Two-column TSV: cell_type, gene
- Pre-filter to your significance threshold before passing.

### Microenvironments file (optional)
- Two-column TSV: cell_type, microenvironment
- Restricts tested pairs to same-microenvironment cell types.

## Querying Results

```python
from cellphonedb.utils import search_utils

hits = search_utils.search_analysis_results(
    query_cell_types_1=['CellA', 'CellB'],
    query_cell_types_2=['CellC', 'CellD'],
    query_genes=['TGFBR1', 'CSF1R'],
    query_interactions=['CSF1_CSF1R'],
    query_classifications=['TGF-beta signaling'],
    query_minimum_score=50,
    significant_means=cpdb_results['significant_means'],
    deconvoluted=cpdb_results['deconvoluted'],
    interaction_scores=cpdb_results['interaction_scores'],
    long_format=True,               # melted format, drops NaN rows
)
```

## Visualization (ktplotspy)

CellPhoneDB does not bundle plots. Use `ktplotspy`:

```python
import ktplotspy as kpy
import anndata as ad

adata = ad.read_h5ad('data/normalised_log_counts.h5ad')

# Heatmap: count of significant interactions per cell-type pair
kpy.plot_cpdb_heatmap(
    pvals=cpdb_results['pvalues'],
    degs_analysis=False,
    figsize=(5, 5),
    title="Sum of significant interactions",
)

# Dot plot: specific interactions across cell-type pairs
kpy.plot_cpdb(
    adata=adata,
    cell_type1="CellA|CellB",
    cell_type2="CellC|CellD",
    means=cpdb_results['means'],
    pvals=cpdb_results['pvalues'],
    celltype_key="cell_labels",
    genes=["TGFB2", "CSF1R"],
    figsize=(10, 3),
    degs_analysis=False,
    standard_scale=True,
    interaction_scores=cpdb_results['interaction_scores'],
    scale_alpha_by_interaction_scores=True,
)
```

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `counts_data` | `'ensembl'` | **Must match your gene ID type.** Common mistake: leaving default when using symbols. |
| `threshold` | `0.1` | Min fraction of cells expressing gene |
| `iterations` | `1000` | Permutations; use 100 for testing |
| `pvalue` | `0.05` | Cutoff for significant_means |
| `score_interactions` | `False` | Enable v5 specificity scoring (0-100) |
| `debug_seed` | `-1` | Set >=0 for reproducibility (requires `threads=1`) |

## Output Format

All output tables have annotation columns (rows = interactions) + cell-type pair columns (`cellA|cellB`).

| Table | Content |
|-------|---------|
| `means` | Mean expression of interaction per cell-type pair |
| `pvalues` | Permutation p-values (Method 2) |
| `significant_means` | Means where p < threshold (NaN otherwise) |
| `relevant_interactions` | Binary 0/1 relevance (Method 3) |
| `deconvoluted` | Per-gene breakdown for complexes (mean expression) |
| `deconvoluted_percents` | Per-gene % of cells expressing |
| `interaction_scores` | Specificity score 0-100 (if enabled) |

## Pitfalls

1. **Human genes only.** Mouse data must be converted to human orthologs. Mismatched IDs → zero interactions found.

2. **`counts_data` default is `'ensembl'`**, not symbol. Explicitly set `counts_data='hgnc_symbol'` when using gene symbols. This is the #1 cause of "no interactions found."

3. **Do not z-scale counts.** Scoring requires log-normalized data where zeros remain zero.

4. **Output directory must exist** before running. CellPhoneDB will not create it.

5. **Do not use numeric cell type names** (e.g., `0`, `1`, `2`). Use string labels.

6. **No dashes in barcodes.** Dashes in cell names cause known issues.

7. **DEG method does no filtering.** All genes in your DEGs file are treated as significant — pre-filter externally.

8. **Interaction directionality is asymmetric.** `IL12_IL12R` for `A|B` (receptor in B) ≠ `B|A` (receptor in A). Both are tested.

9. **Score vs significance can disagree.** High score + non-significant = broadly expressed but not cell-type-specific. Significant + low score = specific but at minimum expression level.

10. **`debug_seed` only works single-threaded.** Set `threads=1` for reproducibility.

11. **Database v5 required for v5 package.** Don't mix old DB ZIPs with v5 code.
