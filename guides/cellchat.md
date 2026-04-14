# CellChat Guide

## When to Use

- You want **pathway-level hierarchy** for CCC (not just LR pairs, but signaling pathways)
- You need **rich built-in visualization** (circle plots, chord diagrams, heatmaps, river plots)
- You want to **compare two conditions** with differential interaction analysis
- You need **communication pattern recognition** (NMF decomposition)
- You need **spatial CCC** in R (v2 for Visium, v3/SpatialCellChat for single-cell resolution)
- Your pipeline is R/Seurat-based

## Installation

```r
devtools::install_github("jinworks/CellChat")
# For single-cell resolution spatial:
# devtools::install_github("jinworks/SpatialCellChat")
```

## Single-Dataset Workflow

### Step 1: Create CellChat object

```r
library(CellChat)
library(patchwork)

# From Seurat object
data.input = seurat_obj[["RNA"]]$data    # Seurat v5+ (use @data for v4)
labels = Idents(seurat_obj)
meta = data.frame(labels = labels, row.names = names(labels))

cellchat = createCellChat(object = data.input, meta = meta, group.by = "labels")
# Or directly: createCellChat(object = seurat_obj, group.by = "ident", assay = "RNA")
```

### Step 2: Set database

```r
CellChatDB = CellChatDB.human   # or CellChatDB.mouse

# Recommended: use all protein signaling (exclude non-protein)
CellChatDB.use = subsetDB(CellChatDB)
# Or specific category:
# CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

cellchat@DB = CellChatDB.use
```

### Step 3: Preprocessing

```r
cellchat = subsetData(cellchat)                          # REQUIRED: subset to signaling genes
future::plan("multisession", workers = 4)                # parallel processing

cellchat = identifyOverExpressedGenes(cellchat)          # DEG analysis
cellchat = identifyOverExpressedInteractions(cellchat)   # map DEGs to LR pairs
```

### Step 4: Inference

```r
cellchat = computeCommunProb(cellchat, type = "triMean")
# type options: "triMean" (default, conservative), "truncatedMean" (with trim=0.1, more permissive)

cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = computeCommunProbPathway(cellchat)            # aggregate to pathway level
cellchat = aggregateNet(cellchat)                        # compute count/weight matrices
```

### Step 5: Systems analysis

```r
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
```

## Output Access

```r
# LR-level results as data.frame
df.net = subsetCommunication(cellchat)

# Pathway-level
df.netP = subsetCommunication(cellchat, slot.name = "netP")

# Filtered by source/target
df.net = subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))

# Filtered by pathway
df.net = subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

# Raw probability arrays
cellchat@net$prob    # 3D: source × target × LR-pair
cellchat@net$pval    # corresponding p-values
cellchat@net$count   # K×K count matrix
cellchat@net$weight  # K×K weight matrix
cellchat@netP$prob   # 3D: source × target × pathway
```

## Visualization

### Overview plots

```r
groupSize = as.numeric(table(cellchat@idents))

# Circle plot — number of interactions
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE)

# Circle plot — interaction strength
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE)

# Heatmap
netVisual_heatmap(cellchat, measure = "weight")
```

### Pathway-specific

```r
pathways.show = c("CXCL")

# Hierarchy plot
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = seq(1,4))

# Circle
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# LR contribution to pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
```

### Bubble plot

```r
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# Filter by pathway
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),
                 signaling = c("CCL", "CXCL"))
```

### Signaling roles

```r
# 2D scatter: outgoing vs incoming strength per cell type
netAnalysis_signalingRole_scatter(cellchat)

# Heatmap of roles
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
```

### Communication patterns (NMF)

```r
library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")  # determine optimal k
cellchat = identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 4)
netAnalysis_river(cellchat, pattern = "outgoing")   # alluvial plot
netAnalysis_dot(cellchat, pattern = "outgoing")     # dot plot
```

## Multi-Dataset Comparison

```r
# Run CellChat on each dataset separately, then merge
object.list = list(Ctrl = cellchat.ctrl, Disease = cellchat.disease)
cellchat = mergeCellChat(object.list, add.names = names(object.list))

# Compare total interactions
compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))

# Differential network (red = increased, blue = decreased in dataset 2)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)

# Information flow comparison
gg1 = rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 = rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
gg1 + gg2

# Bubble plot comparison
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),
                 comparison = c(1,2), max.dataset = 2, remove.isolate = TRUE)
```

### Identify up/down-regulated LR pairs via DEG

```r
pos.dataset = "Disease"
features.name = paste0(pos.dataset, ".merged")

cellchat = identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1,
                                       thresh.fc = 0.05, thresh.p = 0.05)
net = netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up = subsetCommunication(cellchat, net = net, datasets = "Disease", ligand.logFC = 0.05)
net.down = subsetCommunication(cellchat, net = net, datasets = "Ctrl", ligand.logFC = -0.05)
```

## Spatial Workflow (CellChat v2)

```r
# Prepare spatial info
spatial.locs = Seurat::GetTissueCoordinates(visium.obj, scale = NULL,
                                             cols = c("imagerow", "imagecol"))
scalefactors = jsonlite::fromJSON("scalefactors_json.json")
spot.size = 65  # um
conversion.factor = spot.size / scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size / 2)

# Meta must include 'samples' column
meta = data.frame(labels = Idents(visium.obj),
                  samples = "sample1",
                  row.names = names(Idents(visium.obj)))

cellchat = createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs,
                           spatial.factors = spatial.factors)

# Inference with spatial parameters
cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE,
                              interaction.range = 250,   # max diffusion range (um)
                              scale.distance = 0.01,
                              contact.dependent = TRUE,
                              contact.range = 100)       # 100 for Visium, ~10 for single-cell
```

## Key Parameters

| Parameter | Values | Notes |
|-----------|--------|-------|
| `type` (computeCommunProb) | `"triMean"`, `"truncatedMean"` | triMean is conservative (~25% trimmed); truncatedMean with `trim=0.1` is more permissive |
| `trim` | 0–0.25 | Only used with truncatedMean |
| `population.size` | TRUE/FALSE | TRUE for unsorted cells (weights by proportion); FALSE for FACS-sorted |
| `nboot` | ~100 | Permutations for p-values; 20 for quick test |
| `distance.use` | TRUE/FALSE | Enable spatial distance decay |
| `interaction.range` | 250 (um) | Max range for secreted signals |
| `contact.range` | 100 (Visium) / 10 (single-cell) | Range for cell-cell contact interactions |

## Pitfalls

1. **Never use the `integrated` assay from Seurat** — it has negative values. Use `RNA` or `SCT`.

2. **Seurat v5 API change**: Use `seurat_obj[["RNA"]]$data` (not `@data`).

3. **`subsetData()` is always required**, even when using the full database.

4. **triMean zeros out genes expressed in <25% of a cell group.** If expected pathways are missing, switch to `truncatedMean` with `trim = 0.1`.

5. **Non-protein Signaling** in CellChatDB v2 uses estimated proxy expression (not directly measured). Exclude with `subsetDB(CellChatDB)` (no search argument) for standard analysis.

6. **Objects from CellChat <1.6.0** must be updated with `updateCellChat()`.

7. **Spatial `meta` must include a `samples` column** (even for single-sample data).

8. **`variable.both = FALSE`** for spatial data — require only ligand OR receptor to be over-expressed (more permissive, recommended).

9. **Functional similarity analysis requires same cell type composition** across datasets. For different compositions, use structural similarity or `liftCellChat()`.

10. **CellChat v3 / SpatialCellChat** is a separate package (`jinworks/SpatialCellChat`) for single-cell resolution spatial inference.
