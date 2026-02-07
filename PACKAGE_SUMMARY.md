# scPairs Package Summary (Version 0.1.0)

**Date:** February 7, 2026  
**Author:** Zhaoqing Wang (Shandong University)  
**License:** MIT  

---

## Package Overview

**scPairs** is an R package for discovering synergistic gene pairs in single-cell RNA-seq and spatial transcriptomics data. It integrates multiple lines of evidence to identify gene pairs that work cooperatively in biological systems.

### Core Innovation

Unlike traditional single-gene marker analysis, scPairs identifies **gene-gene relationships** through:
- Multi-metric co-expression analysis (linear, rank-based, robust)
- Non-linear dependency detection (mutual information)
- Biological validation (ratio consistency across cell types)
- Spatial validation (Lee's L, co-location quotient)

---

## Package Structure

### Directory Layout

```
scPairs/
├── R/                  # 14 R source files (~3500 lines total)
├── man/                # 26 Rd documentation files
├── tests/testthat/     # 5 test files (44 tests, 100% pass rate)
├── DESCRIPTION         # Package metadata
├── NAMESPACE           # Exports & imports
├── README.md           # User documentation (9.7 KB)
├── NEWS.md             # Version history (7.3 KB)
├── LICENSE & LICENSE.md
└── .Rbuildignore, .gitignore, scPairs.Rproj
```

### R Source Files

| File | Lines | Purpose |
|------|-------|---------|
| **assess_pair.R** | 263 | In-depth assessment of specific gene pairs |
| **coexpression.R** | 273 | Correlation & mutual information metrics |
| **find_all_pairs.R** | 215 | Global discovery workflow |
| **find_gene_pairs.R** | 238 | Gene-centric query workflow |
| **score_integration.R** | 189 | Multi-metric score combination |
| **spatial_adjacency.R** | 241 | Lee's L spatial correlation |
| **spatial_colocation.R** | 163 | Co-location quotient |
| **plot_network.R** | 157 | Gene interaction networks |
| **plot_heatmap.R** | 132 | Synergy score heatmaps |
| **plot_pair_dimplot.R** | 315 | UMAP/tSNE/scatter/violin plots |
| **plot_spatial.R** | 193 | Spatial tissue co-expression maps |
| **utils.R** | 212 | Helper functions & validation |
| **globals.R** | 20 | Global variable declarations |
| **onattach.R** | 16 | Package startup message |

**Total:** ~2,627 lines of R code (excluding tests & documentation)

---

## Three Core Workflows

### 1. Global Discovery (`FindAllPairs`)

**Purpose:** Discover all synergistic gene pairs genome-wide

**Input:**
- Seurat object
- Number of genes to analyze (default: VariableFeatures)
- Optional: permutation count, custom weights

**Output:** `scPairs_result` object with:
- All gene pairs ranked by synergy score
- Individual metrics (7 columns)
- Statistical significance (p-value, p_adj, confidence)

**Use Case:** Exploratory analysis, network construction, hypothesis generation

---

### 2. Gene Query (`FindGenePairs`)

**Purpose:** Find all partners for a specific gene of interest

**Input:**
- Seurat object
- Query gene name (e.g., "TP53")
- Optional: candidate list, top N results

**Output:** `scPairs_gene_result` object with:
- Ranked list of partner genes
- All synergy metrics
- Query metadata

**Use Case:** Targeted investigation, pathway analysis, validation studies

---

### 3. Pair Assessment (`AssessGenePair`)

**Purpose:** In-depth analysis of a specific gene pair

**Input:**
- Seurat object
- Two gene names (e.g., "CD8A", "CD8B")
- Permutation count (default: 999)

**Output:** `scPairs_pair_result` object with:
- Comprehensive metrics (named list)
- Per-cluster correlation analysis
- Permutation p-value
- Jaccard index

**Use Case:** Validation of known pairs, publication-quality analysis

---

## Multi-Evidence Framework

### Seven Metrics (5 Co-Expression + 2 Spatial)

| # | Metric | Weight | Purpose |
|---|--------|--------|---------|
| 1 | **Pearson correlation** | 1.0 | Linear monotonic association |
| 2 | **Spearman correlation** | 1.0 | Rank-based, resistant to outliers |
| 3 | **Biweight midcorrelation** | **1.5** | Robust to scRNA-seq dropout noise |
| 4 | **Mutual information** | 1.0 | Non-linear statistical dependence |
| 5 | **Ratio consistency** | **1.2** | Expression ratio stability across clusters |
| 6 | **Lee's L statistic** | **1.5** | Bivariate spatial autocorrelation |
| 7 | **Co-location quotient** | **1.2** | Spatial proximity enrichment |

**Score Integration:**
```
synergy_score = Σ(w_i × rank_norm(metric_i)) / Σw_i
```

**Statistical Testing:**
- Permutation-based null distribution (cell label shuffling)
- P-value: `(n_null ≥ observed + 1) / (n_perm + 1)`
- Benjamini-Hochberg FDR adjustment
- Confidence: High (p_adj < 0.01), Medium (< 0.05), Low (< 0.1), NS

---

## Visualization Functions (6 Total)

### 1. `PlotPairNetwork()`
**Type:** Gene interaction network  
**Technology:** ggraph + igraph  
**Features:** Edge width = synergy, node size = degree, force-directed layout

### 2. `PlotPairHeatmap()`
**Type:** Synergy score matrix  
**Technology:** pheatmap-style  
**Features:** Hierarchical clustering, symmetric matrix, color gradient

### 3. `PlotPairDimplot()`
**Type:** UMAP/tSNE 3-panel co-expression  
**Technology:** ggplot2 + patchwork  
**Features:** Gene1 | Gene2 | Co-expression product

### 4. `PlotPairSpatial()`
**Type:** Spatial tissue map (3-panel)  
**Technology:** ggplot2 + patchwork  
**Features:** Spatial coordinates, viridis colors, auto-detect spatial data

### 5. `PlotPairViolin()`
**Type:** Cluster-level distribution  
**Technology:** ggplot2  
**Features:** Side-by-side violins by cluster, shows per-cluster patterns

### 6. `PlotPairScatter()`
**Type:** Gene-gene scatter  
**Technology:** ggplot2 + ggExtra (optional)  
**Features:** Cell-level scatter, marginal densities, co-expression color

---

## Technical Specifications

### Performance

**Optimizations:**
- Vectorized correlation matrices (no per-pair loops)
- Sparse matrix support (preserves Seurat format)
- data.table backend (fast pair operations)
- Fast spatial KNN via RANN::nn2() [O(n log n)]

**Typical Runtimes** (Intel i7, 16GB RAM):
- 1000 genes × 5000 cells: ~30 sec (no permutation)
- 500 genes × 10000 cells: ~45 sec (no permutation)
- Add 5-10 min for 999 permutations

**Memory Requirements:**
- Moderate datasets (<2000 genes): 8GB sufficient
- Large datasets (>3000 genes): 16GB+ recommended

### Compatibility

**R Version:** ≥ 4.1.0

**Seurat Support:**
- Seurat v4: Full support
- Seurat v5: Full support (handles slot vs. layer naming)

**Data Types:**
- scRNA-seq (10X Chromium, Drop-seq, Smart-seq2, etc.)
- Spatial transcriptomics (Visium, MERFISH, Slide-seq, seqFISH, etc.)

**Matrix Formats:**
- Sparse: dgCMatrix, dgRMatrix
- Dense: matrix, Matrix

### Dependencies

**Required (11 packages):**
- Seurat (≥4.0), SeuratObject
- data.table, ggplot2, ggraph, igraph
- Matrix, patchwork, tidygraph, tidyr, stats

**Optional (5 packages):**
- RANN (fast spatial KNN)
- ggExtra (marginal density plots)
- pheatmap, crayon, testthat (≥3.0.0)

---

## Testing & Quality Assurance

### Test Suite

**44 unit tests** across 5 test files:

| Test File | Tests | Coverage |
|-----------|-------|----------|
| test-find-all-pairs.R | 12 | Global discovery workflow |
| test-find-gene-pairs.R | 7 | Gene query workflow |
| test-assess-pair.R | 13 | Pair assessment workflow |
| test-internal-functions.R | 12 | Metric functions, utilities |
| helper-test_data.R | - | Test data generation |

**Test Status:** ✅ 44 PASS | 0 FAIL | 0 WARN | 0 SKIP

**Test Coverage:**
- ✅ Main exported functions (FindAllPairs, FindGenePairs, AssessGenePair)
- ✅ Internal metrics (.bicor, .mutual_info, .ratio_consistency)
- ✅ Spatial metrics (.compute_spatial_lee, .compute_spatial_clq)
- ✅ Validation functions (.validate_seurat, .has_spatial)
- ✅ Edge cases (constant vectors, insufficient data, missing genes)
- ✅ Error handling (non-Seurat input, duplicate genes, etc.)

### R CMD check Results

**Latest Check:** February 7, 2026

```
0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

**Note:** "unable to verify current time" (system-level, harmless)

**CRAN Readiness:** ✅ Package meets CRAN submission standards

---

## Documentation

### User Documentation

1. **README.md** (9.7 KB, 334 lines)
   - Quick start for all 3 workflows
   - Methodology details
   - Output structure reference
   - Performance benchmarks
   - Citation information

2. **NEWS.md** (7.3 KB, 216 lines)
   - Version 0.1.0 release notes
   - Feature summary
   - Technical specifications
   - Known limitations
   - Future roadmap

3. **PACKAGE_SUMMARY.md** (this document)
   - Comprehensive package overview
   - Technical architecture
   - Usage patterns

### Developer Documentation

1. **26 Rd manual pages** (roxygen2-generated)
   - 9 exported functions with examples
   - 17 internal functions with technical details
   - Cross-referenced hyperlinks

2. **Inline code comments**
   - Roxygen2 headers for all functions
   - Algorithm explanations
   - Parameter descriptions

3. **NAMESPACE** (auto-generated)
   - 8 exported functions
   - 4 S3 print methods
   - 9 importFrom declarations

---

## Output Objects (S3 Classes)

### `scPairs_result`

**From:** `FindAllPairs()`

**Structure:**
```r
$pairs           # data.table (n_pairs rows)
  - gene1, gene2
  - cor_pearson, cor_spearman, cor_biweight
  - mi_score, ratio_consistency
  - spatial_lee_L, spatial_clq (if spatial)
  - synergy_score, rank, p_value, p_adj, confidence
$parameters      # list (reproducibility metadata)
$n_genes         # integer
$n_cells         # integer
$has_spatial     # logical
```

### `scPairs_gene_result`

**From:** `FindGenePairs()`

**Structure:**
```r
$query_gene      # character
$pairs           # data.table (ranked partners)
$parameters      # list
$n_candidates    # integer
$n_cells         # integer
$has_spatial     # logical
```

### `scPairs_pair_result`

**From:** `AssessGenePair()`

**Structure:**
```r
$gene1, $gene2   # character
$metrics         # named list (all individual metrics)
$per_cluster     # data.frame (per-cluster analysis)
$synergy_score   # numeric [0, 1]
$p_value         # numeric (or NA)
$p_adj           # numeric (or NA)
$confidence      # character ("High", "Medium", "Low", "NS")
$jaccard_index   # numeric [0, 1]
$has_spatial     # logical
$n_cells         # integer
```

All classes have custom `print()` methods for user-friendly summaries.

---

## Citations & References

### Package Citation

```
Wang Z (2026). scPairs: Go Beyond Marker Genes – Discover Synergistic Gene
Pairs in scRNA-seq and Spatial Maps. R package version 0.1.0.
https://github.com/zhaoqing-wang/scPairs
```

### Methodological References

1. **Langfelder & Horvath (2012)** - Biweight midcorrelation  
   *BMC Bioinformatics* 13:328

2. **Lee (2001)** - Bivariate spatial autocorrelation  
   *International Journal of Geographical Information Science* 15(5):413-437

3. **Leslie & Kronenfeld (2011)** - Co-location quotient  
   *Geographical Analysis* 43(3):306-326

---

## Known Limitations & Future Plans

### Current Limitations

1. **Computational Cost:**
   - Permutation testing adds 5-10 minutes
   - Large gene sets (>3000) require significant RAM

2. **Spatial Requirements:**
   - Need ≥50 cells with coordinates
   - k-NN computation scales with cell count

3. **No Parallel Processing:**
   - Single-core computation only
   - Future versions will support multi-core

### Roadmap (v0.2.0+)

**Planned Features:**
- [ ] Vignettes (long-form tutorials)
- [ ] Parallel processing (foreach/future)
- [ ] Additional metrics (transfer entropy, distance correlation)
- [ ] Interactive Shiny app
- [ ] Export functions (CSV, Excel, Cytoscape)
- [ ] Benchmarking datasets (PBMC, mouse brain)
- [ ] pkgdown website

---

## Contact & Support

**Author:** Zhaoqing Wang  
**Affiliation:** Shandong University  
**Email:** zhaoqingwang@mail.sdu.edu.cn  
**ORCID:** [0000-0001-8348-7245](https://orcid.org/0000-0001-8348-7245)

**Repository:** https://github.com/zhaoqing-wang/scPairs  
**Issues:** https://github.com/zhaoqing-wang/scPairs/issues

---

## License

MIT License - Copyright (c) 2026 Zhaoqing Wang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal in the Software
without restriction, including the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software.

---

**Document Version:** 1.0  
**Last Updated:** February 7, 2026  
**Package Version:** 0.1.0
