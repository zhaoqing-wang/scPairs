# scPairs 0.1.0 (2026-02-07)

## Initial Release 🎉

First public release of **scPairs** - a comprehensive R package for discovering synergistic gene pairs in single-cell RNA-seq and spatial transcriptomics data.

### Core Features

#### **Three Complementary Workflows**

* **`FindAllPairs()`** - Global discovery of all synergistic gene pairs
  - Analyzes all pairwise combinations of selected genes
  - Returns ranked list of gene pairs by composite synergy score
  - Optional permutation testing for statistical significance
  - Supports custom metric weights

* **`FindGenePairs()`** - Query-centric partner discovery
  - Find all genes that synergize with a specific gene of interest
  - Efficient computation (only pairs involving query gene)
  - Ideal for hypothesis-driven research

* **`AssessGenePair()`** - In-depth assessment of specific pairs
  - Comprehensive metrics for a single gene pair
  - Per-cluster correlation analysis
  - Permutation-based p-values (default n_perm = 999)
  - Jaccard index for expression overlap

#### **Multi-Evidence Integration Framework**

* **Five Co-Expression Metrics:**
  - Pearson correlation (linear association)
  - Spearman correlation (rank-based, robust to outliers)
  - Biweight midcorrelation (resistant to scRNA-seq dropout)
  - Mutual information (non-linear dependencies)
  - Ratio consistency (stoichiometric stability across clusters)

* **Two Spatial Metrics** (auto-detected):
  - Lee's L statistic (bivariate spatial autocorrelation)
  - Co-location quotient (CLQ, spatial proximity enrichment)

* **Score Integration:**
  - Rank-normalization of all metrics to [0, 1]
  - Weighted summation with customizable weights
  - Default weights emphasize robust metrics (biweight: 1.5, spatial: 1.5)

* **Statistical Rigor:**
  - Permutation-based empirical p-values
  - Benjamini-Hochberg FDR correction
  - Confidence categories: High (p_adj < 0.01), Medium (< 0.05), Low (< 0.1), NS

#### **Publication-Ready Visualizations**

* **`PlotPairNetwork()`** - Gene interaction network
  - ggraph-based network layout
  - Edge width proportional to synergy score
  - Node size = degree (number of connections)
  
* **`PlotPairHeatmap()`** - Synergy score matrix
  - Symmetric heatmap of top gene pairs
  - Hierarchical clustering of genes
  - Color scale from low to high synergy

* **`PlotPairDimplot()`** - UMAP/tSNE co-expression overlay
  - 3-panel layout: gene1, gene2, co-expression product
  - Compatible with any Seurat dimensionality reduction
  - Optional marginal density plots (requires ggExtra)

* **`PlotPairSpatial()`** - Spatial tissue maps
  - 3-panel spatial co-expression visualization
  - Auto-detects spatial coordinates from Seurat object
  - Supports Visium, MERFISH, Slide-seq, etc.

* **`PlotPairViolin()`** - Cluster-level distributions
  - Side-by-side violin plots by cluster
  - Shows expression patterns across cell types

* **`PlotPairScatter()`** - Gene-gene scatter plots
  - Cell-level expression scatter
  - Optional marginal density plots (ggExtra)
  - Color by co-expression product

#### **Performance Optimizations**

* **Vectorized correlation computation** - Full correlation matrices via `cor()` on transposed expression
* **Sparse matrix support** - Preserves Seurat's sparse format until necessary
* **data.table backend** - Fast pair-level operations with minimal memory overhead
* **Fast spatial KNN** - RANN::nn2() for O(n log n) neighbor search (with fallback)
* **Efficient filtering** - Early pruning of pairs with insufficient co-expression

### Technical Details

#### **Compatibility**

* **R version:** ≥ 4.1.0
* **Seurat compatibility:** v4 and v5 (handles both slot and layer naming)
* **Data types:** scRNA-seq, spatial transcriptomics (Visium, MERFISH, Slide-seq, etc.)
* **Matrix formats:** Sparse (dgCMatrix) and dense matrices

#### **Dependencies**

* **Imports:** Seurat (≥4.0), SeuratObject, data.table, ggplot2, ggraph, igraph, Matrix, patchwork, tidygraph, tidyr, stats
* **Suggests:** RANN, ggExtra, pheatmap, crayon, testthat (≥3.0.0)

#### **Testing**

* **44 unit tests** covering:
  - All three main workflows (FindAllPairs, FindGenePairs, AssessGenePair)
  - Internal metric functions (bicor, mutual_info, ratio_consistency)
  - Spatial metric computation
  - Input validation and error handling
  - Edge cases (constant vectors, insufficient co-expression, etc.)

* **Test coverage:** Main functions, internal utilities, visualization functions, error handling

### Documentation

* **Comprehensive README.md** with:
  - Quick start guide for all three workflows
  - Detailed methodology description
  - Output structure reference
  - Performance benchmarks
  - Citation information

* **26 Rd manual pages** (roxygen2-generated):
  - 9 exported functions with examples
  - 17 internal functions with technical details
  - Cross-referenced documentation

* **Package-level help:** `?scPairs` provides overview and citation

* **Startup message:** Citation reminder when loading package

### Known Limitations

* **Computational cost:** Permutation testing (999 perms) adds 5-10 minutes for typical datasets
* **Memory usage:** Dense correlation matrices for large gene sets (>3000 genes) may require 16+ GB RAM
* **Spatial metrics:** Require at least 50 cells with spatial coordinates; k-NN computation scales O(n log n)

### Future Roadmap (v0.2.0+)

Planned features for upcoming releases:

* **Vignettes:** Long-form tutorials for scRNA-seq and spatial workflows
* **Parallel processing:** Multi-core support for permutation testing
* **Additional metrics:** Transfer entropy, distance correlation, conditional MI
* **Interactive visualizations:** Shiny app for exploratory analysis
* **Export functions:** Save results to CSV, Excel, or Cytoscape format
* **Benchmarking datasets:** Pre-computed results on public data (PBMC, mouse brain, etc.)
* **pkgdown website:** Full documentation site with searchable reference

---

## Installation

```r
# Install from GitHub
devtools::install_github("zhaoqing-wang/scPairs")
```

---

## Citation

Wang Z (2026). scPairs: Go Beyond Marker Genes – Discover Synergistic Gene Pairs in scRNA-seq and Spatial Maps. R package version 0.1.0. https://github.com/zhaoqing-wang/scPairs

---

**Full Changelog:** https://github.com/zhaoqing-wang/scPairs/commits/main
