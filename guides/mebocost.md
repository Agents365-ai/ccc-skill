# MEBOCOST Guide

## When to Use

- You suspect **metabolite-mediated communication** (lipids, amino acids, sugars, ketone bodies) — biology invisible to protein LR tools
- You study **immunometabolism**, **tumor metabolism**, or **metabolic diseases**
- You want to go beyond protein ligand-receptor pairs to capture paracrine metabolic crosstalk
- Your data is **scRNA-seq** (human or mouse)

**Key difference from LR tools**: LIANA+/CellPhoneDB/CellChat score protein ligand→receptor pairs. MEBOCOST scores **metabolite→sensor** pairs, where metabolite abundance is inferred from enzyme gene expression and sensors include GPCRs, transporters, cell-surface enzymes, and ion channels.

## Installation

```bash
conda create -n mebocost python=3.12
conda activate mebocost
git clone https://github.com/kaifuchenlab/MEBOCOST.git
cd MEBOCOST
pip install -r requirements.txt
python -m pip install .
```

**Important**: Edit `mebocost.conf` to point to your local data directory (the bundled config has hardcoded developer paths).

## Core Workflow

### Step 1: Create object

```python
from mebocost import mebocost
import scanpy as sc

adata = sc.read_h5ad("data.h5ad")
# Use raw counts with all genes (not just HVGs)
# If adata has been filtered: adata = adata.raw.to_adata()

mebo_obj = mebocost.create_obj(
    adata=adata,
    group_col='celltype',           # column in adata.obs (must be string, not list)
    species='human',                # 'human' or 'mouse'
    config_path='/path/to/MEBOCOST/mebocost.conf',
    met_est='mebocost',            # metabolite estimation method
    cutoff_exp='auto',             # sensor expression cutoff (auto = 25th pct of non-zeros)
    cutoff_met='auto',             # metabolite enzyme cutoff
    cutoff_prop=0.15,              # min fraction of cells expressing met/sensor
    sensor_type='All',             # or ['receptor','transporter','enzyme','channel']
    thread=4,
)
```

### Step 2: Infer communication

```python
commu_res = mebo_obj.infer_commu(
    n_shuffle=1000,                 # permutations
    seed=12345,
    Return=True,
    thread=4,
    save_permuation=True,           # REQUIRED if you plan differential analysis later
    pval_method='permutation_test_fdr',
    pval_cutoff=0.05,
    min_cell_number=50,             # exclude small cell groups
)
```

### Step 3 (optional): Constrain by metabolic flux

```python
# Using COMPASS output
mebo_obj._ConstrainCompassFlux_(
    compass_folder='/path/to/compass_output/',
    efflux_cut='auto',
    influx_cut='auto',
    inplace=True,
)

# Or using any FBA tool output
mebo_obj._ConstrainFluxFromAnyTool_(
    efflux_mat=efflux_df,           # rows=metabolite names, cols=cell groups
    influx_mat=influx_df,
    efflux_cut='auto',
    influx_cut='auto',
    inplace=True,
)
# Adds 'Flux_PASS' column: 'PASS', 'UNPASS', or 'N/A'
```

### Save/load

```python
mebocost.save_obj(mebo_obj, path='mebocost_result.pk')
mebo_obj = mebocost.load_obj('mebocost_result.pk')
```

## Output Format

`mebo_obj.commu_res` — DataFrame where each row is a significant mCCC event:

| Column | Meaning |
|--------|---------|
| `Sender`, `Receiver` | Cell group names |
| `Metabolite` | HMDB ID |
| `Metabolite_Name` | Human-readable name (e.g., L-Arginine) |
| `Sensor` | Sensor gene symbol (e.g., SLC7A1) |
| `Annotation` | Sensor type: Receptor, Transporter, Enzyme, Channel |
| `Commu_Score` | avg_met × avg_sensor |
| `Norm_Commu_Score` | Normalized by mean background |
| `met_in_sender` | Mean metabolite enzyme expression in sender |
| `sensor_in_receiver` | Mean sensor expression in receiver |
| `metabolite_prop_in_sender` | Fraction of sender cells above cutoff |
| `sensor_prop_in_receiver` | Fraction of receiver cells above cutoff |
| `permutation_test_pval`, `_fdr` | Permutation p-value and FDR |
| `Flux_PASS` | Flux constraint result (after Step 3) |

## Visualization

```python
# Event count bar chart
mebo_obj.eventnum_bar(
    pval_cutoff=0.05,
    flux_pass=True,
    include=['sender-receiver', 'sensor', 'metabolite', 'metabolite-sensor'],
)

# Dot heatmap
mebo_obj.commu_dotmap(
    sender_focus=['Macrophage'],
    receiver_focus=['T_cell'],
    pval_cutoff=0.05,
    flux_pass=True,
    save='dotmap.pdf',
)

# Flow plot (Sender → Metabolite → Sensor → Receiver)
mebo_obj.FlowPlot(
    metabolite_focus=['Cholesterol'],
    save='flow.pdf',
)

# Network plot
mebo_obj.commu_network_plot(
    pval_cutoff=0.05,
    flux_pass=True,
    save='network.pdf',
)

# Violin/heatmap of sensor or metabolite distribution
mebo_obj.violin_plot(
    sensor_or_met=['CD36', 'L-Arginine'],
    cell_focus=['Macrophage', 'T_cell'],
)

# Summary dot plot (sender × receiver counts)
mebo_obj.count_dot_plot(pval_cutoff=0.05, save='count_dot.pdf')
```

## Multi-Sample / Differential Analysis

```python
# Option A: single object with condition column
mebo_obj = mebocost.create_obj(
    adata=adata_combined,
    group_col='celltype',
    condition_col='condition',      # e.g., 'Tumor' or 'Normal'
    species='human',
    config_path='/path/to/mebocost.conf',
)
commu_res = mebo_obj.infer_commu(n_shuffle=1000, save_permuation=True)

# Option B: merge two objects
merged = mebocost.concat_obj(obj_tumor, obj_normal, cond1='Tumor', cond2='Normal')

# Differential analysis
diff_result = mebo_obj.CommDiff(
    comps=['Tumor_vs_Normal'],      # cond1_vs_cond2 (cond1 = numerator in Log2FC)
    sig_mccc_only=True,
    flux_pass=True,
    thread=8,
    Return=True,
)
# Result in mebo_obj.diffcomm_res['Tumor_vs_Normal']

# Differential plots
mebo_obj.DiffSummaryPlot(comp_cond='Tumor_vs_Normal', pval_cutoff=0.05)
mebo_obj.DiffFlowPlot(comp_cond='Tumor_vs_Normal', pval_cutoff=0.05)
mebo_obj.CompScatterPlot(comp_cond='Tumor_vs_Normal')
```

## Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `met_est` | `'mebocost'` | Metabolite estimation method. Also supports scFEA and COMPASS flux |
| `cutoff_exp` | `'auto'` | Sensor expression threshold (auto = 25th percentile of non-zeros) |
| `cutoff_prop` | `0.15` | Min fraction of cells above cutoff; groups below get p=1 |
| `min_cell_number` | `50` | Cell groups with fewer cells are excluded |
| `n_shuffle` | `1000` | Permutations; use 100 for testing |
| `sensor_type` | `'All'` | Filter: `['receptor', 'transporter', 'enzyme', 'channel']` |
| `save_permuation` | `False` | **Must be True** for CommDiff |

## Pitfalls

1. **Use all genes, not just HVGs.** If `adata.X` has only highly variable genes (<5000), enzyme genes for many metabolites are missing. Use `adata = adata.raw.to_adata()` first.

2. **`group_col` must be a string, not a list (v1.2+).** Passing a list raises `KeyError`. Old v1.0 pickles are handled by backward-compatibility logic.

3. **Edit `mebocost.conf` paths.** The bundled config has hardcoded developer paths. You must update all paths to your local MEBOCOST data directory.

4. **`save_permuation=True` before CommDiff.** Without the permutation background saved, differential analysis cannot run. You'd need to re-run `infer_commu`.

5. **Avoid `~` in condition or cell-type names.** The `~` character is used internally as a separator.

6. **Flux constraint name matching.** `_ConstrainFluxFromAnyTool_` matches metabolite names via a synonym dictionary. Use names from `mebo_obj.met_ann['metabolite']` for compatibility.

7. **Cell groups <50 cells are silently excluded.** Check `commu_res` to verify expected cell types appear.

8. **Receptor sensors only check efflux, not influx.** For `Annotation == 'Receptor'`, only sender efflux is checked (GPCRs don't import the metabolite). This is intentional but can surprise users.

9. **Differential analysis requires identical cell type names across conditions.** Mismatched names (e.g., `"T cell"` vs `"T_cell"`) cause silent exclusion.

10. **Human and mouse only.** Other species are not supported as of v1.2.2.
