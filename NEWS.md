# scPairs 0.1.6 (2026-02-11)

## Unified Computation Engine & API Consolidation

-   **Shared metric engine** (`compute_metrics.R`): all three discovery functions (`FindAllPairs`, `FindGenePairs`, `AssessGenePair`) now share a single `.compute_pair_metrics()` backend. Eliminates \~600 lines of duplicated neighbourhood/prior/spatial metric code that was previously copy-pasted across the three entry points.

-   **`mode` parameter** for selective computation:

    -   `"all"` (default) -- compute all 14 metrics.
    -   `"expression"` -- cell-level co-expression and neighbourhood metrics only; skips prior knowledge.
    -   `"prior_only"` -- prior knowledge scores only; no expression computation (fast pathway-based screening).

-   **Unified output format**: `AssessGenePair()` now includes a `$pairs` data.table with the same column schema as `FindAllPairs()` / `FindGenePairs()` (`gene1`, `gene2`, `synergy_score`, `rank`, `confidence`, plus all metric columns). All visualization functions accept any of the three result classes interchangeably.

-   **Unified visualization input**: `PlotPairNetwork()`, `PlotPairHeatmap()`, and all other plot functions now use a shared `.extract_pairs_df()` helper that dispatches on `scPairs_result`, `scPairs_gene_result`, `scPairs_pair_result`, or plain `data.frame`.

-   **Removed dead dependency**: `pheatmap` dropped from `Suggests` (unused).

-   **Internal refactoring**:

    -   `.resolve_cluster_ids()` extracted to eliminate repeated cluster-resolution code.
    -   `.extract_pair_vals()` for symmetric matrix lookups.
    -   `.extract_pairs_df()` for unified result-to-data.frame conversion.

------------------------------------------------------------------------

# scPairs 0.1.5 (2026-02-11)

## Prior Knowledge Integration & Synergy-Aware Scoring

-   **New evidence layer: Prior biological knowledge** -- shifts discovery from pure co-expression towards functional synergy.
    -   `prior_score`: Jaccard similarity of GO-BP and KEGG pathway co-annotations (weight 2.0).
    -   `bridge_score`: Pathway bridge genes connecting the focal pair through shared intermediaries (weight 1.8).
    -   Supports GO, KEGG, and user-supplied interaction databases via `custom_pairs`.
-   **New metric: `neighbourhood_synergy`** -- directional paracrine enrichment score (weight 1.5).
-   **New visualization: `PlotPairSynergy()`** -- 6-panel figure with neighbourhood synergy, bridge gene network, cluster expression, and multi-evidence bar chart.
-   **14 total metrics** across 5 evidence layers.
-   All three main functions gain `use_prior`, `organism`, and `custom_pairs` parameters.
-   **Optional dependencies:** `org.Mm.eg.db`, `org.Hs.eg.db`, `AnnotationDbi`.

------------------------------------------------------------------------

# scPairs 0.1.4 (2026-02-08)

## Robustness & Validation

-   Fixed spurious cross-cell-type correlations with sparse genes via `min_pct_expressed`.
-   Comprehensive input validation with informative error messages.
-   Code consolidation: unified `.compute_*()` dispatch functions; centralized schema (`schema.R`); shared plot utilities (`plot_utils.R`).
-   Sparse matrix optimization for KNN smoothing.

------------------------------------------------------------------------

# scPairs 0.1.3 (2026-02-08)

## Cross-Cell-Type Interaction Metric

-   **New metric: `cross_celltype_score`** -- trans-cellular synergy detection via micro-environment binning (weight 1.5).
-   **New visualization: `PlotPairCrossType()`** -- directed cell-type pair heatmap.

------------------------------------------------------------------------

# scPairs 0.1.2 (2026-02-08)

## Neighbourhood-Aware Metrics

-   3 new metrics: `smoothed_cor`, `neighbourhood_score`, `cluster_cor`.
-   New visualizations: `PlotPairSmoothed()`, `PlotPairSummary()`.

------------------------------------------------------------------------

# scPairs 0.1.1 (2026-02-08)

## Performance Optimizations

-   5--20x speedups through vectorisation of all core metric computations.

------------------------------------------------------------------------

# scPairs 0.1.0 (2026-02-07)

## Initial Release

-   3 core workflows: `FindAllPairs()`, `FindGenePairs()`, `AssessGenePair()`.
-   7 metrics: Pearson, Spearman, biweight, MI, ratio consistency, Lee's L, CLQ.
-   6 visualization functions.
