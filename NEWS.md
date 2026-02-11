# scPairs 0.1.5 (2026-02-11)

## Prior Knowledge Integration & Synergy-Aware Scoring

-   **New evidence layer: Prior biological knowledge** — shifts discovery from pure co-expression towards functional synergy.
    -   `prior_score`: Jaccard similarity of GO Biological Process and KEGG pathway co-annotations between gene pairs (weight 2.0). Genes sharing many functional terms are up-weighted as synergy candidates.
    -   `bridge_score`: Identifies pathway bridge genes — intermediary genes C that share functional annotations with BOTH genes A and B AND are expressed in the current dataset (weight 1.8). Captures indirect synergy (A→C→B) through shared pathway intermediaries, critical for detecting pairs like Adora2a-Ido1 that synergise through shared immunosuppressive pathway context.
    -   Supports GO (Biological Process), KEGG, and user-supplied interaction databases (CellChatDB, CellPhoneDB, SCENIC regulon targets) via the `custom_pairs` parameter.
-   **New metric: `neighbourhood_synergy`** — directional paracrine enrichment score (weight 1.5). For cells highly expressing gene A, measures whether gene B's expression in their biological neighbourhood is enriched beyond expectation (and vice versa). Captures juxtacrine/paracrine-like interactions where synergy operates between neighbouring cells rather than within the same cell.
-   **New visualization: `PlotPairSynergy()`** — publication-ready 6-panel figure integrating:
    1.  UMAP expression maps for each gene and their neighbourhood synergy signal
    2.  Prior knowledge bridge gene network (igraph/ggraph) showing functional connections through intermediary genes
    3.  Per-cluster expression comparison
    4.  Multi-evidence metric bar chart (expression + neighbourhood + prior knowledge)
-   **14 total metrics** across 5 evidence layers (cell-level, neighbourhood, prior knowledge, trans-cellular, spatial).
-   All three main functions (`FindAllPairs()`, `FindGenePairs()`, `AssessGenePair()`) gain `use_prior`, `organism`, and `custom_pairs` parameters.
-   `AssessGenePair()` now returns bridge genes, shared GO/KEGG terms, and prior network object for downstream analysis.
-   **Optional dependencies:** `org.Mm.eg.db`, `org.Hs.eg.db`, `AnnotationDbi` (for prior knowledge; gracefully skipped if unavailable).

------------------------------------------------------------------------

# scPairs 0.1.4 (2026-02-08)

## Robustness & Validation

-   **Fixed spurious cross-cell-type correlations** with sparse genes by adding `min_pct_expressed` expression validation.
-   **Comprehensive input validation** across all analysis and plotting functions with informative error messages.
-   **Code consolidation** — unified internal dispatch functions (`.compute_*()`) reduce duplication by \~40%; new centralized schema (`schema.R`) and shared plot utilities (`plot_utils.R`).
-   **Sparse matrix optimization** — KNN smoothing defers densification; vectorised cross-cell-type binning.
-   16 new tests (115 total).

------------------------------------------------------------------------

# scPairs 0.1.3 (2026-02-08)

## Cross-Cell-Type Interaction Metric

-   **New metric: `cross_celltype_score`** — detects trans-cellular synergies (paracrine signalling, ligand-receptor pairs) by correlating gene A in one cell type with gene B in neighbouring cells of a different type. Score = geometric mean of bidirectional correlations; weight 1.5.
-   **New visualization: `PlotPairCrossType()`** — heatmap of directed cell-type pair correlations.
-   Integrated into `FindAllPairs()`, `FindGenePairs()`, `AssessGenePair()`, and permutation testing.

------------------------------------------------------------------------

# scPairs 0.1.2 (2026-02-08)

## Neighbourhood-Aware Metrics

-   **3 new metrics** to overcome scRNA-seq dropout:
    -   `smoothed_cor` — Pearson correlation on KNN-smoothed expression (weight 1.5)
    -   `neighbourhood_score` — average partner expression in KNN neighbours (weight 1.5)
    -   `cluster_cor` — cluster-level pseudo-bulk correlation (weight 1.2)
-   **New visualizations:** `PlotPairSmoothed()` (6-panel raw + smoothed UMAP), `PlotPairSummary()` (comprehensive multi-panel figure).
-   All three main functions gain `use_neighbourhood`, `neighbourhood_k`, and `smooth_alpha` parameters.

------------------------------------------------------------------------

# scPairs 0.1.1 (2026-02-08)

## Performance Optimizations

-   **5-20x speedups** through vectorisation of all core metric computations:
    -   Co-expression filter via `tcrossprod()` (\~30x)
    -   Biweight midcorrelation via matrix kernel operations (\~18x)
    -   Mutual information via pre-binned `outer()` (\~4x)
    -   Lee's L and CLQ via sparse matrix multiply (\~15x)
    -   Permutation p-values vectorised across all pairs

------------------------------------------------------------------------

# scPairs 0.1.0 (2026-02-07)

## Initial Release

-   **3 core workflows:** `FindAllPairs()` (global discovery), `FindGenePairs()` (query-centric), `AssessGenePair()` (detailed assessment).
-   **7 metrics:** Pearson, Spearman, biweight midcorrelation, mutual information, ratio consistency, Lee's L, co-location quotient (CLQ).
-   **Score integration:** rank-normalised weighted summation with permutation p-values and FDR correction.
-   **6 visualization functions:** `PlotPairNetwork()`, `PlotPairHeatmap()`, `PlotPairDimplot()`, `PlotPairSpatial()`, `PlotPairViolin()`, `PlotPairScatter()`.
-   Supports Seurat v4/v5, scRNA-seq and spatial transcriptomics (Visium, MERFISH, Slide-seq).
-   44 unit tests.
