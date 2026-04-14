# NicheNet Guide

## When to Use

- You want to predict **which ligands caused the observed transcriptional changes** in a receiver cell type
- You need **ligand → downstream target gene** prediction (not just LR interactions)
- You have a **case-control design** and want to know which sender ligands drive condition-specific gene expression
- You need **multi-criteria prioritization** of LR pairs (DE, expression, activity)
- For multi-sample multi-condition: use **MultiNicheNet** (pseudobulk + edgeR)

## Installation

```r
devtools::install_github("saeyslab/nichenetr")
# For multi-sample:
devtools::install_github("saeyslab/multinichenetr")
```

## Prior Knowledge Setup

```r
library(nichenetr)
library(Seurat)
library(tidyverse)

organism = "mouse"  # or "human"

# Download pre-built models from Zenodo
if (organism == "mouse") {
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
} else {
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
}

lr_network = lr_network %>% distinct(from, to)

# IMPORTANT: fix gene aliases before analysis
seuratObj = alias_to_symbol_seurat(seuratObj, organism)
```

## Step-by-Step Workflow (Case-Control)

### Step 1: Define expressed genes

```r
receiver = "CD8 T"
sender_celltypes = c("CD4 T", "Treg", "Mono", "NK", "B", "DC")

expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.05)
expressed_genes_sender = sender_celltypes %>%
  lapply(get_expressed_genes, seuratObj, 0.05) %>%
  unlist() %>% unique()
```

### Step 2: Define potential ligands

```r
all_receptors = unique(lr_network$to)
expressed_receptors = intersect(all_receptors, expressed_genes_receiver)

potential_ligands = lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>% unique()

# Sender-focused (recommended):
potential_ligands = intersect(potential_ligands, expressed_genes_sender)
```

### Step 3: Define gene set of interest (DE genes in receiver)

```r
seurat_obj_receiver = subset(seuratObj, idents = receiver)
DE_table = FindMarkers(seurat_obj_receiver,
                       ident.1 = "Disease", ident.2 = "Control",
                       group.by = "condition", min.pct = 0.05)
DE_table = DE_table %>% rownames_to_column("gene")

geneset_oi = DE_table %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
  pull(gene)

# CRITICAL: filter to genes in the model
geneset_oi = geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]

background_expressed_genes = expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)]

# Check sizes: geneset ~50-500, background ~5000-10000, ratio ~1/10 to 1/200
```

### Step 4: Predict ligand activity

```r
ligand_activities = predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Primary ranking metric in v2: aupr_corrected
ligand_activities = ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands = ligand_activities %>%
  top_n(30, aupr_corrected) %>%
  pull(test_ligand)
```

### Step 5: Infer target genes

```r
active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# For heatmap visualization
active_ligand_target_links = prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33    # quantile cutoff; lower = denser heatmap
)
```

### Step 6: Infer receptors

```r
ligand_receptor_links_df = get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig
)
```

### Step 7: Extended prioritization (v2)

```r
info_tables = generate_info_tables(
  seuratObj,
  celltype_colname = "celltype",
  senders_oi = sender_celltypes,
  receivers_oi = receiver,
  lr_network = lr_network %>% filter(from %in% expressed_genes_sender & to %in% expressed_receptors),
  condition_colname = "condition",
  condition_oi = "Disease",
  condition_reference = "Control",
  scenario = "case_control"
)

prior_table = generate_prioritization_tables(
  info_tables$sender_receiver_info,
  info_tables$sender_receiver_de,
  ligand_activities,
  info_tables$lr_condition_de,
  scenario = "case_control"
)
# Output: tibble with prioritization_score and 52 component columns
```

## Wrapper Function (Quick Version)

```r
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj,
  receiver = "CD8 T",
  condition_colname = "condition",
  condition_oi = "Disease",
  condition_reference = "Control",
  sender = c("CD4 T", "Treg", "Mono", "NK", "B", "DC"),
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)
# Returns list: top_ligands, ligand_activities, ligand_target_matrix,
#   ligand_target_heatmap, ligand_receptor_matrix, ...
```

## Visualization

```r
# Ligand-target heatmap
make_heatmap_ggplot(active_ligand_target_links,
                    y_name = "Predicted target genes",
                    x_name = "Prioritized ligands",
                    color = "purple",
                    legend_title = "Regulatory potential")

# Ligand expression dot plot
DotPlot(subset(seuratObj, idents = sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") +
  coord_flip()

# Signaling path
ligand_tf_matrix = construct_ligand_tf_matrix(
  weighted_networks, lr_network, best_upstream_ligands[1:5],
  ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5
)

signaling_path = get_ligand_signaling_path(
  ligand_tf_matrix, ligands_all = "BMP2", targets_all = "HEY1",
  top_n_regulators = 4, weighted_networks = weighted_networks
)
```

## MultiNicheNet (Multi-Sample)

For proper multi-sample statistics (pseudobulk + edgeR).

```r
library(multinichenetr)
library(SingleCellExperiment)

# Convert from Seurat
sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")

# Define contrasts
contrasts_oi = c("'Disease-(Control)/1'")
contrast_tbl = tibble(contrast = c("Disease-(Control)/1"), group = c("Disease"))

# Cell type and gene filtering
abundance_info = get_abundance_info(sce, sample_id, group_id, celltype_id,
                                     min_cells = 10, senders_oi, receivers_oi)

frq_list = get_frac_exprs(sce, sample_id, celltype_id, group_id,
                           batches = NA, min_cells = 10,
                           fraction_cutoff = 0.05, min_sample_prop = 0.50)

# Pseudobulk DE
abundance_expression_info = process_abundance_expression_info(
  sce, sample_id, group_id, celltype_id,
  min_cells = 10, senders_oi, receivers_oi,
  lr_network, batches = NA, frq_list, abundance_info
)

DE_info = get_DE_info(sce, sample_id, group_id, celltype_id,
                       batches = NA, covariates = NA, contrasts_oi,
                       min_cells = 10, expressed_df = frq_list$expressed_df)

# Ligand activity + prioritization
ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
  receiver_de = DE_info$celltype_de$de_output_tidy,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = 0.50, p_val_threshold = 0.05
)

prioritization_tables = generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes,
  contrast_tbl, sender_receiver_tbl, grouping_tbl,
  scenario = "regular"
)
```

## Key Parameters

| Parameter | Notes |
|-----------|-------|
| `pct = 0.05` (get_expressed_genes) | Min fraction of cells expressing gene; intentionally lenient |
| `aupr_corrected` | Primary ranking metric in v2 (not Pearson from v1) |
| Top N ligands | 30 is arbitrary; inspect AUPR distribution for natural cutoff |
| `cutoff = 0.33` (prepare_ligand_target_visualization) | Quantile for target filtering; lower = denser heatmap |
| `n = 250` (get_weighted_ligand_target_links) | Top targets per ligand; use 100 for small gene sets |
| MultiNicheNet: min 4 samples per group | Below this, pseudobulk benefits are unclear |

## Pitfalls

1. **Always run `alias_to_symbol_seurat()` first.** Old gene aliases silently drop genes from the model.

2. **Mouse model = converted human model.** Mouse-specific genes without human orthologs are absent.

3. **`geneset_oi` must capture ligand-driven changes, not cell-intrinsic programs.** Use condition DE (Disease vs Control within receiver), NOT cluster markers. Cluster markers conflate intrinsic + extrinsic programs.

4. **Filter `geneset_oi` to `rownames(ligand_target_matrix)`.** Genes not in the model are silently ignored.

5. **Gene set size matters.** Too large (>500) dilutes signal; too small (<20) lacks power. Ideal ratio to background: 1/200 to 1/10.

6. **Top ligands may be lowly expressed.** Ligand activity (AUPR) is independent of expression level. Always cross-check with DotPlot.

7. **Sender-agnostic vs sender-focused can give different rankings.** Run both and compare. Sender-agnostic may return IFN genes not expressed by defined senders.

8. **Cell type names must be syntactically valid R names.** Spaces and hyphens will break internal functions.

9. **MultiNicheNet `contrasts_oi` format is strict.** No spaces, single quotes per contrast, double quotes around all.

10. **The prior model is static.** It captures known biology as of the model build date. Novel interactions won't be scored.
