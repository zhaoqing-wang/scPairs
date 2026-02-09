# scPairs

## Discover Synergistic Gene Pairs in Single-Cell and Spatial Transcriptomics

[![R Version](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![GitHub Package Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/scPairs?label=GitHub&color=blue)](https://github.com/zhaoqing-wang/scPairs/releases) [![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-green)](https://github.com/zhaoqing-wang)

**scPairs** goes beyond single-marker analysis to identify cooperative gene pairs that jointly define cell states, drive pathways, or form functional complexes. It integrates 11 complementary metrics spanning cell-level co-expression, neighbourhood-aware smoothing, trans-cellular interactions, and spatial co-variation, then ranks gene pairs by a composite synergy score with optional permutation-based significance testing.

## Installation

``` r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("zhaoqing-wang/scPairs")
```

**Required:** Seurat (\>=4.0), SeuratObject, data.table, ggplot2, ggraph, igraph, Matrix, patchwork, tidygraph, tidyr

**Optional:** RANN (fast KNN), ggExtra (marginal densities), pheatmap, crayon

## Quick Start

``` r
library(scPairs)
library(Seurat)

sce <- readRDS("your_data.rds")
```

### 1. Global Pair Discovery

Find top synergistic gene pairs genome-wide:

``` r
result <- FindAllPairs(sce, n_top_genes = 1000, top_n = 200)
result

# With permutation testing
result <- FindAllPairs(sce, n_top_genes = 500, n_perm = 999, top_n = 100)

# Visualize
PlotPairNetwork(result, top_n = 50)
PlotPairHeatmap(result, top_n = 25)
```

### 2. Query-Centric Partner Search

Find synergistic partners for a gene of interest:

``` r
tp53_partners <- FindGenePairs(sce, gene = "TP53", top_n = 20)
tp53_partners

PlotPairNetwork(tp53_partners, top_n = 15)
PlotPairDimplot(sce, gene1 = "TP53", gene2 = tp53_partners$pairs$gene2[1])
```

### 3. In-Depth Pair Assessment

Comprehensively evaluate a specific gene pair:

``` r
assessment <- AssessGenePair(sce, gene1 = "CD8A", gene2 = "CD8B", n_perm = 999)
assessment

PlotPairSmoothed(sce, gene1 = "CD8A", gene2 = "CD8B")
PlotPairSummary(sce, gene1 = "CD8A", gene2 = "CD8B", result = assessment)
PlotPairCrossType(sce, gene1 = "CD8A", gene2 = "CD8B")
```

### Spatial Transcriptomics

Spatial metrics are **automatically detected** for Visually, MERFISH, Slide-seq, etc.:

``` r
result <- FindAllPairs(spatial_obj, n_top_genes = 500)
result$has_spatial  # TRUE

PlotPairSpatial(spatial_obj, gene1 = "EPCAM", gene2 = "KRT8")
```

## Metrics

scPairs evaluates each gene pair through four layers of evidence:

| Layer | Metrics | What It Captures |
|------------------------|------------------------|------------------------|
| **Cell-level** | Pearson, Spearman, Biweight midcorrelation, Mutual information, Ratio consistency | Direct co-expression, non-linear dependence, stoichiometric stability |
| **Neighbourhood** | KNN-smoothed correlation, Neighbourhood co-expression, Cluster pseudo-bulk correlation | Co-regulation beyond dropout noise, population-level patterns |
| **Trans-cellular** | Cross-cell-type interaction score | Paracrine signalling, ligand-receptor pairs across cell types |
| **Spatial** | Lee's L statistic, Co-location quotient (CLQ) | Bivariate spatial autocorrelation, proximity enrichment |

All metrics are rank-normalised to [0, 1] and combined via weighted summation:

```         
synergy_score = sum(w_i * rank_norm(metric_i)) / sum(w_i)
```

Higher weights (1.5) are assigned to dropout-robust metrics (biweight, smoothed correlation, neighbourhood score, cross-cell-type, Lee's L). Permutation testing provides empirical p-values with Benjamini-Hochberg FDR correction.

## Functions

| Function              | Purpose                                        |
|-----------------------|------------------------------------------------|
| `FindAllPairs()`      | Global discovery of all synergistic gene pairs |
| `FindGenePairs()`     | Find partners for a specific gene              |
| `AssessGenePair()`    | Detailed assessment of a single pair           |
| `PlotPairNetwork()`   | Gene interaction network                       |
| `PlotPairHeatmap()`   | Synergy score heatmap                          |
| `PlotPairDimplot()`   | UMAP/tSNE co-expression overlay (3-panel)      |
| `PlotPairSmoothed()`  | Raw + KNN-smoothed UMAP (6-panel)              |
| `PlotPairSummary()`   | Comprehensive multi-panel figure               |
| `PlotPairSpatial()`   | Spatial co-expression map (3-panel)            |
| `PlotPairCrossType()` | Cross-cell-type interaction heatmap            |
| `PlotPairViolin()`    | Expression distributions by cluster            |
| `PlotPairScatter()`   | Gene-gene scatter plot                         |

## Output

All functions return S3 objects with custom `print()` methods.

**`FindAllPairs()` / `FindGenePairs()`** return a ranked `data.table` of gene pairs with all metric values, composite synergy scores, and (optionally) p-values, FDR-adjusted p-values, and confidence levels (High/Medium/Low/NS).

**`AssessGenePair()`** returns per-metric scores, per-cluster correlations, cross-cell-type breakdowns, permutation p-value, and confidence classification.

## Performance

Vectorised implementations (tcrossprod, sparse matrix multiply, pre-binned MI) yield 5-20x speedups over naive approaches. Typical runtimes (Intel i7, 16GB RAM):

| Scale                    | Time (no permutation) |
|--------------------------|-----------------------|
| 1000 genes x 5000 cells  | \~5-10s               |
| 2000 genes x 10000 cells | \~30-60s              |
| \+ 999 permutations      | add \~2-5 min         |

## Citation

``` bibtex
@Manual{scPairs,
  title = {scPairs: Discover Synergistic Gene Pairs in scRNA-seq and Spatial Transcriptomics},
  author = {Zhaoqing Wang},
  year = {2026},
  note = {R package version 0.1.4},
  url = {https://github.com/zhaoqing-wang/scPairs},
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Contact

**Author:** Zhaoqing Wang ([ORCID](https://orcid.org/0000-0001-8348-7245)) **Email:** [zhaoqingwang\@mail.sdu.edu.cn](mailto:zhaoqingwang@mail.sdu.edu.cn){.email} **Issues:** [GitHub Issues](https://github.com/zhaoqing-wang/scPairs/issues)
