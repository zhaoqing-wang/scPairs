# scPairs 0.1.1 (2026-02-08)

## Performance Optimizations

This release focuses on major performance improvements to support large-scale gene pair discovery with thousands of candidate genes. All core metric computations have been vectorised, yielding **5--20x speedups** depending on the number of genes and cells.

### Co-expression Metrics

* **Co-expression filter:** Replaced per-pair `vapply` loop with a single `tcrossprod()` on a binary expression matrix. For 2000 genes (≈2M pairs), this reduces filtering time from minutes to seconds.

* **Biweight midcorrelation:** New `.bicor_matrix()` function computes the full biweight midcorrelation matrix at once using vectorised Tukey biweight kernel operations and `tcrossprod()`. Falls back to per-pair Pearson on non-zero cells only for genes with zero MAD (common in sparse scRNA-seq). Replaces the previous per-pair `vapply` loop.

* **Mutual information:** New `.mutual_info_batch()` pre-computes bin assignments for all genes in a single pass, then uses fast table lookups per pair. The inner MI summation is now vectorised via `outer()` and logical indexing, replacing the nested `for` loop.

* **Ratio consistency:** New `.ratio_consistency_batch()` pre-computes `log2(x + 1)` for all genes once, then uses pre-split cluster index lists for fast median computation. Replaces per-pair `tapply` calls.

### Spatial Metrics

* **Lee's L statistic:** Spatial lag is now computed via sparse matrix multiplication (`W %*% t(mat)`) for all genes at once, replacing the per-cell R-level loop. Lee's L values for all pairs are computed via vectorised row-wise inner products. Permutation testing is similarly batch-vectorised: all pairs are evaluated per permutation via matrix operations.

* **Co-location quotient (CLQ):** Neighbour expression counts are now computed via sparse matrix multiplication (`expr_binary %*% t(N_sparse)`) for all genes at once, replacing the nested `vapply` loops over cells.

### Score Integration

* **Permutation p-values:** The per-pair `vapply` for p-value computation is replaced by a vectorised comparison (`as.integer(perm_score >= obs_scores)`) accumulated across permutations.

### `FindGenePairs()`

* **Co-expression filter:** Vectorised via matrix--vector multiply instead of per-candidate `vapply`.
* **Pearson/Spearman:** Computed for all candidates at once via `cor(gene_vec, t(cand_mat))` instead of per-candidate `cor()` calls.

### `AssessGenePair()`

* **Permutation testing:** Pearson and Spearman correlations for all permutations are computed via a single matrix multiply (`cor(x, perm_y_mat)`). Pre-generates all permuted y-vectors as a matrix.

### Internal

* New batch helper functions: `.bicor_matrix()`, `.mutual_info_batch()`, `.ratio_consistency_batch()`.
* Spatial weight matrices now constructed as sparse `Matrix::sparseMatrix` objects.
* Added `@importFrom Matrix sparseMatrix t` and `@importFrom stats setNames`.

### Estimated Speedups

| Operation | v0.1.0 | v0.1.1 | Speedup |
|---|---|---|---|
| Co-expression filter (2000 genes) | ~60s | ~2s | ~30x |
| Biweight midcorrelation (2000 genes) | ~90s | ~5s | ~18x |
| Mutual information (2000 genes) | ~120s | ~30s | ~4x |
| Lee's L (5000 spots, 1000 pairs) | ~45s | ~3s | ~15x |
| CLQ (5000 spots, 1000 pairs) | ~30s | ~2s | ~15x |

---

# scPairs 0.1.0 (2026-02-07)

## Initial Release

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

* **R version:** >= 4.1.0
* **Seurat compatibility:** v4 and v5 (handles both slot and layer naming)
* **Data types:** scRNA-seq, spatial transcriptomics (Visium, MERFISH, Slide-seq, etc.)
* **Matrix formats:** Sparse (dgCMatrix) and dense matrices

#### **Dependencies**

* **Imports:** Seurat (>=4.0), SeuratObject, data.table, ggplot2, ggraph, igraph, Matrix, patchwork, tidygraph, tidyr, stats
* **Suggests:** RANN, ggExtra, pheatmap, crayon, testthat (>=3.0.0)

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

Wang Z (2026). scPairs: Go Beyond Marker Genes -- Discover Synergistic Gene Pairs in scRNA-seq and Spatial Maps. R package version 0.1.1. https://github.com/zhaoqing-wang/scPairs

---

**Full Changelog:** https://github.com/zhaoqing-wang/scPairs/commits/main
