# scPairs: Identifying Synergistic Gene Pairs in Single-Cell and Spatial Transcriptomics

[![GitHub Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/scPairs?label=Version&color=blue)](https://github.com/zhaoqing-wang/scPairs/releases) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![R ≥ 4.1.0](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/) [![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-green)](https://github.com/zhaoqing-wang)

## Overview

<img src="docs/Sticker.png" alt="scPairs Logo" width="200" align="right"/>

scPairs is an R package for the systematic identification and multi-evidence evaluation of synergistic gene pairs in single-cell and spatial transcriptomics data. Unlike conventional pairwise co-expression analyses that rely on a single correlation metric, scPairs integrates **14 complementary metrics** across **five orthogonal evidence layers** — cell-level co-expression, neighbourhood-aware smoothing, prior biological knowledge (GO/KEGG), trans-cellular interaction, and spatial co-variation — to compute a composite **synergy score** with optional permutation-based significance testing. This multi-scale, multi-evidence design enables researchers to move beyond simple co-expression towards a comprehensive characterization of cooperative gene regulation at both transcriptomic and spatial resolution.

## Table of Contents

1. [Preparation](#1-preparation)
    - [1.1 Installation](#11-installation)
    - [1.2 Dependencies](#12-dependencies)
2. [Quick Start](#2-quick-start)
    - [2.1 Global Pair Discovery](#21-global-pair-discovery)
    - [2.2 Query-Centric Partner Search](#22-query-centric-partner-search)
    - [2.3 In-Depth Pair Assessment](#23-in-depth-pair-assessment)
    - [2.4 Prior Knowledge & Bridge Gene Network](#24-prior-knowledge--bridge-gene-network)
    - [2.5 Spatial Transcriptomics](#25-spatial-transcriptomics)
3. [Architecture](#3-architecture)
    - [3.1 Computation Modes](#31-computation-modes)
    - [3.2 Five Evidence Layers (14 Metrics)](#32-five-evidence-layers-14-metrics)
    - [3.3 Unified Output Schema](#33-unified-output-schema)
4. [Functions](#4-functions)
    - [4.1 Discovery (3 functions)](#41-discovery-3-functions)
    - [4.2 Visualization (11 functions)](#42-visualization-11-functions)
5. [Citation](#5-citation)
6. [License](#6-license)
7. [Contact](#7-contact)

---

## 1. Preparation

### 1.1 Installation

```r
# Install from GitHub (recommended)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("zhaoqing-wang/scPairs")

# For prior knowledge integration (recommended)
BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db", "AnnotationDbi"))
```

### 1.2 Dependencies

**Required:** Seurat (≥ 4.0), data.table, ggplot2, ggraph, ggrepel, igraph, Matrix, patchwork, tidygraph, tidyr

**Optional:** org.Mm.eg.db / org.Hs.eg.db / AnnotationDbi (prior knowledge), RANN (fast KNN), ggExtra (marginal densities)

<details>
<summary><b>Install missing dependencies manually</b></summary>

```r
install.packages(c("data.table", "ggplot2", "ggraph", "ggrepel",
                    "igraph", "Matrix", "patchwork", "Seurat",
                    "tidygraph", "tidyr"))
# Optional
install.packages(c("RANN", "ggExtra"))
BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db", "AnnotationDbi"))
```

</details>

---

## 2. Quick Start

```r
library(scPairs)
library(Seurat)

sce <- readRDS("your_data.rds")
```

### 2.1 Global Pair Discovery

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

### 2.2 Query-Centric Partner Search

```r
tp53_partners <- FindGenePairs(sce, gene = "TP53", top_n = 20)
PlotPairNetwork(tp53_partners)
PlotPairDimplot(sce, gene1 = "TP53", gene2 = tp53_partners$pairs$gene2[1])
```

### 2.3 In-Depth Pair Assessment

```r
assessment <- AssessGenePair(sce, gene1 = "CD8A", gene2 = "CD8B")
assessment

PlotPairSmoothed(sce, gene1 = "CD8A", gene2 = "CD8B")
PlotPairSummary(sce, gene1 = "CD8A", gene2 = "CD8B", result = assessment)

# All visualization functions accept any result class
PlotPairNetwork(assessment)
PlotPairHeatmap(assessment)
```

### 2.4 Prior Knowledge & Bridge Gene Network

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

### 2.5 Spatial Transcriptomics

Spatial metrics are automatically detected for Visium, MERFISH, Slide-seq, etc.:

```r
result <- FindAllPairs(spatial_obj, n_top_genes = 500)
PlotPairSpatial(spatial_obj, gene1 = "EPCAM", gene2 = "KRT8")
```

---

## 3. Architecture

### 3.1 Computation Modes

| Mode | Metrics Computed | Use Case |
|------|-----------------|----------|
| `"all"` (default) | All 14 metrics | Full analysis |
| `"expression"` | Co-expression + neighbourhood | No annotation databases needed |
| `"prior_only"` | Prior knowledge scores only | Fast pathway-based screening |

### 3.2 Five Evidence Layers (14 Metrics)

| Layer | Metrics | Weight |
|-------|---------|--------|
| **Cell-level** (5) | Pearson (`cor_pearson`), Spearman (`cor_spearman`), Biweight midcorrelation (`cor_biweight`), Mutual information (`mi_score`), Ratio consistency (`ratio_consistency`) | 1.0 – 1.5 |
| **Neighbourhood** (5) | KNN-smoothed correlation (`smoothed_cor`), Neighbourhood score (`neighbourhood_score`), Cluster correlation (`cluster_cor`), Cross-cell-type score (`cross_celltype_score`), Neighbourhood synergy (`neighbourhood_synergy`) | 1.2 – 1.5 |
| **Prior knowledge** (2) | GO/KEGG co-annotation Jaccard (`prior_score`), Pathway bridge score (`bridge_score`) | 1.8 – 2.0 |
| **Spatial** (2) | Lee's L (`spatial_lee_L`), Co-location quotient (`spatial_clq`) | 1.2 – 1.5 |

Metrics are rank-normalised to [0, 1] and combined via weighted summation:

$$\text{synergy\\_score} = \frac{\sum_i w_i \cdot \text{rank\\_norm}(m_i)}{\sum_i w_i}$$

### 3.3 Unified Output Schema

All three workflows produce a `pairs` data.table with consistent columns: `gene1`, `gene2`, `synergy_score`, `rank`, `confidence`, plus individual metric columns. Every visualization function accepts any result class interchangeably.

**Confidence assignment:**

| With permutation p-values | Without p-values (score quantiles) |
|---|---|
| High: p_adj < 0.01 | High: ≥ 95th percentile |
| Medium: p_adj < 0.05 | Medium: ≥ 80th percentile |
| Low: p_adj < 0.10 | Low: ≥ 50th percentile |

---

## 4. Functions

### 4.1 Discovery (3 functions)

| Function | Purpose |
|----------|---------|
| `FindAllPairs()` | Global discovery of synergistic gene pairs across all variable genes |
| `FindGenePairs()` | Find synergistic partners for a specific query gene |
| `AssessGenePair()` | In-depth assessment of a specific gene pair with per-cluster detail |

### 4.2 Visualization (11 functions)

| Function | Purpose |
|----------|---------|
| `PlotPairNetwork()` | Synergy-weighted gene interaction network |
| `PlotPairHeatmap()` | Synergy score heatmap across top gene pairs |
| `PlotPairDimplot()` | UMAP co-expression overlay (3-panel) |
| `PlotPairSmoothed()` | Raw + KNN-smoothed UMAP expression (6-panel) |
| `PlotPairSummary()` | Comprehensive multi-panel summary figure |
| `PlotPairSpatial()` | Spatial co-expression map (3-panel) |
| `PlotPairCrossType()` | Cross-cell-type interaction heatmap |
| `PlotPairViolin()` | Expression distributions by cluster |
| `PlotPairScatter()` | Gene–gene scatter plot with marginal densities |
| `PlotPairSynergy()` | Multi-evidence synergy dashboard (6-panel) |
| `PlotBridgeNetwork()` | Bridge gene network with pathway-weighted edges |

---

## 5. Citation

```
Wang Z (2026). scPairs: Identifying Synergistic Gene Pairs in single-Cell and Spatial Transcriptomics.
https://github.com/zhaoqing-wang/scPairs
```

## 6. License

[MIT](LICENSE.md)

## 7. Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) | **Email:** <zhaoqingwang@mail.sdu.edu.cn> | **Issues:** [scPairs Issues](https://github.com/zhaoqing-wang/scPairs/issues)
