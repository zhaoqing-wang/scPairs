# scPairs

## Go Beyond Marker Genes — Discover Synergistic Gene Pairs in scRNA-seq and Spatial Transcriptomics

[![R Version](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub Package Version](https://img.shields.io/github/r-package/v/zhaoqing-wang/scPairs?label=GitHub&color=blue)](https://github.com/zhaoqing-wang/scPairs/releases) 
[![GitHub Maintainer](https://img.shields.io/badge/Maintainer-Zhaoqing_Wang-green)](https://github.com/zhaoqing-wang)

## Overview

**scPairs** identifies synergistic gene pairs in single-cell RNA-seq and spatial transcriptomics data by integrating multiple lines of evidence: co-expression correlation, mutual information, expression ratio consistency, and spatial co-variation.

---

## Why Gene Pairs?

Single-marker analysis misses **cooperative regulatory programs** where two genes jointly define cell states, drive pathways, or form functional complexes. scPairs discovers these relationships through a rigorous, multi-metric framework:

| Evidence Layer | Metric | What It Captures |
|---|---|---|
| **Linear co-expression** | Pearson / Spearman | Monotonic co-regulation |
| **Robust co-expression** | Biweight midcorrelation | Outlier-resistant association<sup>1</sup> |
| **Non-linear dependence** | Mutual information | Arbitrary statistical relationships |
| **Stoichiometric consistency** | Ratio consistency | Stable expression ratios across clusters |
| **Spatial co-variation** | Lee's L statistic<sup>2</sup> | Bivariate spatial autocorrelation |
| **Spatial co-location** | Co-location quotient (CLQ)<sup>3</sup> | Spatial proximity enrichment |

All metrics are rank-normalised and combined via weighted summation. Optional permutation testing provides empirical p-values with FDR correction.

<sup>1</sup> Langfelder & Horvath (2012) *BMC Bioinformatics*  
<sup>2</sup> Lee (2001) *Int J Geogr Inf Sci*  
<sup>3</sup> Leslie & Kronenfeld (2011) *Geogr Anal*

---

## Installation

```r
# Install from GitHub (recommended)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("zhaoqing-wang/scPairs")

# Or install from local source
devtools::install_local("path/to/scPairs")
```

### Dependencies

**Required:** `Seurat` (≥4.0), `SeuratObject`, `data.table`, `ggplot2`, `ggraph`, `igraph`, `Matrix`, `patchwork`, `tidygraph`, `tidyr`, `stats`

**Optional:** `RANN` (fast spatial KNN), `ggExtra` (marginal density plots), `pheatmap`, `crayon`, `testthat`

---

## Quick Start

### Load Package and Data

```r
library(scPairs)
library(Seurat)

# Load your Seurat object (scRNA-seq or spatial)
sce <- readRDS("your_data.rds")
```

### Workflow 1: Global Discovery of All Synergistic Pairs

Find the top synergistic gene pairs genome-wide:

```r
# Fast mode (no permutation testing)
result <- FindAllPairs(sce, n_top_genes = 1000, top_n = 200)
print(result)
head(result$pairs)

# With statistical testing (999 permutations)
result <- FindAllPairs(sce, n_top_genes = 500, n_perm = 999, top_n = 100)

# Visualizations
PlotPairNetwork(result, top_n = 50)   # Gene interaction network
PlotPairHeatmap(result, top_n = 25)   # Synergy score heatmap
```

**Output:** All gene pairs ranked by composite synergy score, with optional p-values and FDR-adjusted significance.

### Workflow 2: Query Partners for a Gene of Interest

Find all genes that synergize with your gene of interest:

```r
# What genes work synergistically with TP53?
tp53_partners <- FindGenePairs(sce, gene = "TP53", top_n = 20)
print(tp53_partners)

# Visualize the partnership network
PlotPairNetwork(tp53_partners, top_n = 15)

# UMAP co-expression with the top partner
top_partner <- tp53_partners$pairs$gene2[1]
PlotPairDimplot(sce, gene1 = "TP53", gene2 = top_partner)
```

**Output:** Ranked list of partner genes with detailed synergy metrics.

### Workflow 3: Assess a Specific Gene Pair

In-depth analysis of a known or suspected gene pair:

```r
# Evaluate the CD8A-CD8B pair with permutation testing
assessment <- AssessGenePair(sce, gene1 = "CD8A", gene2 = "CD8B", n_perm = 999)
print(assessment)

# Multiple visualization options
PlotPairDimplot(sce, gene1 = "CD8A", gene2 = "CD8B")   # UMAP co-expression
PlotPairScatter(sce, gene1 = "CD8A", gene2 = "CD8B")   # Gene-gene scatter
PlotPairViolin(sce, gene1 = "CD8A", gene2 = "CD8B")    # Cluster distributions
```

**Output:** Comprehensive metrics, per-cluster correlations, synergy score, p-value, and confidence level.

### Spatial Transcriptomics

Spatial metrics are **automatically computed** when the Seurat object contains spatial coordinates (Visium, MERFISH, Slide-seq, etc.):

```r
# Spatial data automatically detected
result <- FindAllPairs(spatial_obj, n_top_genes = 500)
result$has_spatial  # TRUE

# Spatial co-expression map (3-panel: gene1, gene2, co-expression)
PlotPairSpatial(spatial_obj, gene1 = "EPCAM", gene2 = "KRT8")
```

**Spatial metrics:** Lee's L (bivariate spatial autocorrelation) and CLQ (co-location quotient) are computed using k-nearest neighbors (default k=6).

---

## Core Functions

| Function | Purpose |
|---|---|
| **`FindAllPairs()`** | Global discovery of all synergistic gene pairs |
| **`FindGenePairs()`** | Query-centric: find partners for one gene |
| **`AssessGenePair()`** | In-depth assessment of a specific gene pair |
| **`PlotPairNetwork()`** | Gene interaction network (ggraph) |
| **`PlotPairHeatmap()`** | Synergy score heatmap (pheatmap-style) |
| **`PlotPairDimplot()`** | UMAP/tSNE co-expression overlay (3-panel) |
| **`PlotPairSpatial()`** | Spatial tissue co-expression map (3-panel) |
| **`PlotPairViolin()`** | Cluster-level expression violin plots |
| **`PlotPairScatter()`** | Cell-level gene-gene scatter (with optional marginals) |

---

## Methodology

### Multi-Evidence Integration

Each gene pair is evaluated using up to 7 metrics:

**Co-Expression Metrics:**
1. **Pearson correlation** – Linear association
2. **Spearman correlation** – Rank-based, monotonic association
3. **Biweight midcorrelation** – Robust to scRNA-seq dropout noise
4. **Mutual information** – Non-linear statistical dependence
5. **Ratio consistency** – Expression ratio stability across cell clusters

**Spatial Metrics** (when applicable):
6. **Lee's L statistic** – Bivariate spatial autocorrelation
7. **Co-location quotient (CLQ)** – Spatial proximity enrichment

### Score Integration Formula

```
synergy_score = Σ(w_i × rank_norm(metric_i)) / Σw_i
```

Each metric is rank-normalised to [0, 1] within the dataset, then combined via weighted summation.

**Default Weights:**

| Metric | Weight | Rationale |
|---|---|---|
| Pearson | 1.0 | Standard linear correlation |
| Spearman | 1.0 | Rank-robust alternative |
| Biweight midcorrelation | **1.5** | Resistant to dropout artifacts |
| Mutual information | 1.0 | Captures non-linear patterns |
| Ratio consistency | **1.2** | Biological plausibility filter |
| Lee's L (spatial) | **1.5** | Strong spatial validation |
| CLQ (spatial) | **1.2** | Spatial co-location evidence |

**Custom weights** can be specified via the `weights` parameter.

### Statistical Significance

- **Permutation testing:** Cell labels are shuffled `n_perm` times (default: 0 for speed; use 999 for publication). P-values calculated as: `(n_null ≥ observed + 1) / (n_perm + 1)`
- **Multiple testing correction:** Benjamini-Hochberg FDR adjustment
- **Confidence categories:**
  - **High:** p_adj < 0.01
  - **Medium:** p_adj < 0.05
  - **Low:** p_adj < 0.1
  - **NS:** Not significant

---

## Output Structure

All result objects are S3 classes with custom `print()` methods for clear summaries:

### `scPairs_result` (from `FindAllPairs()`)

```r
result$pairs           # data.table with columns:
                       #   gene1, gene2, cor_pearson, cor_spearman, cor_biweight,
                       #   mi_score, ratio_consistency, spatial_lee_L, spatial_clq,
                       #   synergy_score, rank, p_value, p_adj, confidence
result$parameters      # Analysis parameters for reproducibility
result$n_genes         # Number of genes analyzed
result$n_cells         # Number of cells
result$has_spatial     # Logical: spatial metrics computed?
```

### `scPairs_gene_result` (from `FindGenePairs()`)

```r
result$query_gene      # The input gene
result$pairs           # data.table of partners ranked by synergy_score
result$parameters      # Analysis parameters
result$n_candidates    # Number of candidate partners tested
result$n_cells         # Number of cells
result$has_spatial     # Logical: spatial metrics computed?
```

### `scPairs_pair_result` (from `AssessGenePair()`)

```r
result$gene1, $gene2   # The query genes
result$metrics         # Named list of all individual metrics
result$per_cluster     # data.frame with per-cluster correlations
result$synergy_score   # Composite score [0, 1]
result$p_value         # Permutation-based p-value (NA if n_perm=0)
result$p_adj           # FDR-adjusted p-value
result$confidence      # "High", "Medium", "Low", or "NS"
result$jaccard_index   # Expression overlap (Jaccard index)
result$has_spatial     # Logical: spatial metrics computed?
result$n_cells         # Number of cells
```

---

## Performance Notes

**v0.1.1** introduces major vectorisation improvements across all metric computations:

- **Co-expression filter:** Uses `tcrossprod()` on binary expression matrices — handles millions of pairs in seconds.
- **Biweight midcorrelation:** Full correlation matrix computed via vectorised Tukey biweight kernel and `tcrossprod()` — no per-pair R-level loops.
- **Mutual information:** Pre-computed bin assignments for all genes, vectorised MI summation via `outer()`.
- **Spatial lag (Lee's L):** Sparse matrix multiplication (`W %*% t(mat)`) for all genes at once — replaces per-cell loops.
- **CLQ:** Neighbour counts via sparse matrix multiply for all genes at once.
- **Permutation p-values:** Vectorised null-vs-observed comparison across all pairs per permutation.
- **Sparse matrix support:** Expression data remains sparse until dense operations are required, minimizing memory footprint.
- **`data.table` backend:** All pair-level computations use `data.table` for speed and efficiency.
- **Fast spatial KNN:** Uses `RANN::nn2()` when available for O(n log n) neighbor search; falls back to base R if unavailable.

**Typical runtimes** (Intel i7, 16GB RAM):
- 1000 genes × 5000 cells: ~5-10 seconds (no permutation)
- 2000 genes × 10000 cells: ~30-60 seconds (no permutation)
- Add ~2-5 minutes for 999 permutations

---

## Citation

If you use scPairs in published research, please cite:

```
Wang Z (2026). scPairs: Go Beyond Marker Genes – Discover Synergistic Gene
Pairs in scRNA-seq and Spatial Maps. R package version 0.1.1.
https://github.com/zhaoqing-wang/scPairs
```

BibTeX:
```bibtex
@Manual{scPairs,
  title = {scPairs: Go Beyond Marker Genes - Discover Synergistic Gene Pairs in scRNA-seq and Spatial Maps},
  author = {Zhaoqing Wang},
  year = {2026},
  note = {R package version 0.1.1},
  url = {https://github.com/zhaoqing-wang/scPairs},
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for full details.

---

## Contact & Bug Reports

**Author:** Zhaoqing Wang  
**Email:** zhaoqingwang@mail.sdu.edu.cn  
**ORCID:** [0000-0001-8348-7245](https://orcid.org/0000-0001-8348-7245)  
**Issues:** [GitHub Issues](https://github.com/zhaoqing-wang/scPairs/issues)

---

## Acknowledgments

scPairs builds on methodologies from:
- Langfelder & Horvath (2012) - Biweight midcorrelation for network analysis
- Lee (2001) - Bivariate spatial autocorrelation statistics
- Leslie & Kronenfeld (2011) - Co-location quotient methodology
