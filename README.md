<table>
  <tr>
    <td>
      <h1>scPairs: Identifying Synergistic Gene Pairs in Single-Cell and Spatial Transcriptomics</h1>
      <p>
        <a href="https://github.com/zhaoqing-wang/scPairs/releases"><img src="https://img.shields.io/github/r-package/v/zhaoqing-wang/scPairs?label=Version&color=blue" alt="GitHub Package Version" /></a>
        <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT" /></a>
        <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-%3E%3D4.1.0-blue" alt="R Version" /></a>
        <a href="https://github.com/zhaoqing-wang"><img src="https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-green" alt="GitHub Maintainer" /></a>
      </p>
    </td>
    <td width="200">
      <img src="docs/Sticker.png" alt="scPairs Logo" width="200" />
    </td>
  </tr>
</table>

## Overview

**scPairs** discovers cooperative gene pairs in scRNA-seq and spatial transcriptomics data. It integrates **14 complementary metrics** across **five evidence layers** — cell-level co-expression, neighbourhood-aware smoothing, prior biological knowledge (GO/KEGG), trans-cellular interaction, and spatial co-variation — to rank gene pairs by a composite **synergy score** with optional permutation-based significance testing. Three unified workflows (global discovery, query-centric search, single-pair assessment) share a common computation engine, and all outputs feed directly into a rich suite of publication-ready visualizations.

### Key Features

- **Multi-layer scoring** — moves beyond simple co-expression to capture neighbourhood context, pathway-level co-annotation, bridge gene connectivity, and spatial co-localization.
- **Three analysis workflows** — `FindAllPairs()` for genome-wide screening, `FindGenePairs()` for query-centric partner discovery, `AssessGenePair()` for in-depth pair evaluation.
- **Selectable computation modes** — run all 14 metrics, expression-only metrics, or prior-knowledge-only screening to balance depth and speed.
- **Unified output schema** — every workflow returns the same `pairs` data.table format; every visualization function accepts any result class interchangeably.
- **14 publication-ready visualizations** — network graphs, heatmaps, UMAP overlays, smoothed expression panels, spatial maps, synergy dashboards, bridge gene networks, and more.
- **Automatic spatial detection** — Visium, MERFISH, Slide-seq and other spatial modalities are detected automatically; spatial metrics are computed without additional configuration.

## Installation

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("zhaoqing-wang/scPairs")

# For prior knowledge integration (recommended)
BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db", "AnnotationDbi"))
```

**Required:** Seurat (≥ 4.0), data.table, ggplot2, ggraph, igraph, Matrix, patchwork, tidygraph, tidyr

**Optional:** org.Mm.eg.db / org.Hs.eg.db / AnnotationDbi (prior knowledge), RANN (fast KNN), ggExtra (marginal densities)

## Quick Start

```r
library(scPairs)
library(Seurat)

sce <- readRDS("your_data.rds")
```

### 1. Global Pair Discovery

```r
# All 14 metrics (default)
result <- FindAllPairs(sce, n_top_genes = 1000, top_n = 200)

# Expression metrics only (no annotation databases needed)
result <- FindAllPairs(sce, mode = "expression")

# Prior knowledge only (fast pathway-based screening)
result <- FindAllPairs(sce, mode = "prior_only", organism = "mouse")

# With permutation testing
result <- FindAllPairs(sce, n_top_genes = 500, n_perm = 999)

# Visualize
PlotPairNetwork(result, top_n = 50)
PlotPairHeatmap(result, top_n = 25)
```

### 2. Query-Centric Partner Search

```r
tp53_partners <- FindGenePairs(sce, gene = "TP53", top_n = 20)
PlotPairNetwork(tp53_partners)
PlotPairDimplot(sce, gene1 = "TP53", gene2 = tp53_partners$pairs$gene2[1])
```

### 3. In-Depth Pair Assessment

```r
assessment <- AssessGenePair(sce, gene1 = "CD8A", gene2 = "CD8B")
assessment

PlotPairSmoothed(sce, gene1 = "CD8A", gene2 = "CD8B")
PlotPairSummary(sce, gene1 = "CD8A", gene2 = "CD8B", result = assessment)

# All visualization functions accept any result class
PlotPairNetwork(assessment)
PlotPairHeatmap(assessment)
```

### 4. Prior Knowledge & Bridge Gene Network

```r
# Requires: BiocManager::install(c("org.Mm.eg.db", "AnnotationDbi"))
result <- AssessGenePair(sce, gene1 = "Adora2a", gene2 = "Ido1",
                          organism = "mouse")

# 6-panel synergy dashboard (expression + neighbourhood + prior evidence)
PlotPairSynergy(sce, gene1 = "Adora2a", gene2 = "Ido1")

# Standalone bridge gene network
PlotBridgeNetwork(sce, gene1 = "Adora2a", gene2 = "Ido1",
                   organism = "mouse", top_bridges = 15)

# Custom interaction databases
custom_db <- data.frame(gene1 = c("Adora2a"), gene2 = c("Ido1"))
result <- FindGenePairs(sce, gene = "Adora2a", custom_pairs = custom_db)
```

### 5. Spatial Transcriptomics

Spatial metrics are automatically detected for Visium, MERFISH, Slide-seq, etc.:

```r
result <- FindAllPairs(spatial_obj, n_top_genes = 500)
PlotPairSpatial(spatial_obj, gene1 = "EPCAM", gene2 = "KRT8")
```

## Architecture

### Computation Modes

| Mode | Metrics Computed | Use Case |
|------|-----------------|----------|
| `"all"` (default) | All 14 metrics | Full analysis |
| `"expression"` | Co-expression + neighbourhood | No annotation databases |
| `"prior_only"` | Prior knowledge scores only | Fast pathway-based screening |

### Five Evidence Layers (14 Metrics)

| Layer | Metrics | Weight |
|-------|---------|--------|
| **Cell-level** | Pearson, Spearman, Biweight midcorrelation, Mutual information, Ratio consistency | 1.0 – 1.5 |
| **Neighbourhood** | KNN-smoothed correlation, Neighbourhood score, Cluster correlation, Cross-cell-type score, Neighbourhood synergy | 1.2 – 1.5 |
| **Prior knowledge** | GO/KEGG co-annotation (Jaccard), Pathway bridge score | 1.8 – 2.0 |
| **Spatial** | Lee's L, Co-location quotient | 1.2 – 1.5 |

Metrics are rank-normalised to [0, 1] and combined via weighted summation:

$$\text{synergy\\_score} = \frac{\sum_i w_i \cdot \text{rank\\_norm}(m_i)}{\sum_i w_i}$$

### Unified Output

All three workflows produce a `pairs` data.table with consistent columns: `gene1`, `gene2`, `synergy_score`, `rank`, `confidence`, plus individual metric columns. Every visualization function accepts any result class interchangeably.

## Functions

### Discovery

| Function | Purpose |
|----------|--------|
| `FindAllPairs()` | Global discovery of synergistic gene pairs |
| `FindGenePairs()` | Find partners for a specific query gene |
| `AssessGenePair()` | Detailed assessment of a specific pair |

### Visualization

| Function | Purpose |
|----------|--------|
| `PlotPairNetwork()` | Synergy-weighted gene interaction network |
| `PlotPairHeatmap()` | Synergy score heatmap |
| `PlotPairDimplot()` | UMAP co-expression overlay (3-panel) |
| `PlotPairSmoothed()` | Raw + KNN-smoothed UMAP (6-panel) |
| `PlotPairSummary()` | Comprehensive multi-panel summary figure |
| `PlotPairSpatial()` | Spatial co-expression map (3-panel) |
| `PlotPairCrossType()` | Cross-cell-type interaction heatmap |
| `PlotPairViolin()` | Expression distributions by cluster |
| `PlotPairScatter()` | Gene–gene scatter plot with marginal densities |
| `PlotPairSynergy()` | Multi-evidence synergy dashboard (6-panel) |
| `PlotBridgeNetwork()` | Bridge gene network with pathway-weighted edges |

## Performance

Vectorised implementations yield 5–20× speedups over naive approaches.

| Scale | Time (no permutation) |
|-------|-----------------------|
| 1 000 genes × 5 000 cells | ~ 5–10 s |
| 2 000 genes × 10 000 cells | ~ 30–60 s |
| + 999 permutations | + ~ 2–5 min |

## Citation

```bibtex
@Manual{scPairs,
  title = {scPairs: Identifying Synergistic Gene Pairs in Single-Cell
           and Spatial Transcriptomics},
  author = {Zhaoqing Wang},
  year = {2026},
  url = {https://github.com/zhaoqing-wang/scPairs},
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [GitHub Issues](https://github.com/zhaoqing-wang/scPairs/issues)
