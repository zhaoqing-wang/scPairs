# scPairs 0.1.8 (2026-02-28)

## Summary

First CRAN submission release. Incorporates all CRAN reviewer feedback, a built-in
test dataset, complete function documentation, an improved test suite, and an
overhauled package structure modelled on established CRAN bioinformatics packages.

------------------------------------------------------------------------

## CRAN Compliance

-   **Added reference URL** to the `DESCRIPTION` `Description` field
    (`<https://github.com/zhaoqing-wang/scPairs>`).

-   **Added `\value` tags** to all three exported `print()` S3 methods
    (`print.scPairs_result`, `print.scPairs_gene_result`,
    `print.scPairs_pair_result`), each returning the input object invisibly.

-   **All 11 runnable examples unwrapped**: `\donttest{}` removed entirely;
    examples now execute as plain R code during `R CMD check` (all complete
    in under 10 s on standard hardware). Three functions retain `\dontrun{}`
    because they require software not guaranteed to be present:
    `PlotPairSpatial` (spatial-assay Seurat object), `PlotBridgeNetwork` and
    `PlotPairSynergy` (both require Bioconductor annotation packages
    `org.Hs.eg.db` / `org.Mm.eg.db`).

-   **Removed `set.seed()`** from inside `.plot_bridge_network_enhanced()`
    (`R/plot_synergy.R`); the package no longer sets a user-visible seed.

------------------------------------------------------------------------

## Built-in Test Dataset

-   **New dataset `scpairs_testdata`** (`data/scpairs_testdata.rda`, 26 KB):
    a synthetic Seurat object (100 cells × 20 genes, 3 balanced clusters,
    5-component PCA, 2-D UMAP) with two injected co-expression patterns:

    -   **GENE3 & GENE4** — globally correlated across all cells
        (Pearson r ≈ 0.89 in normalised data).
    -   **GENE1 & GENE2** — moderately correlated within cluster 1 only.
    -   All remaining pairs: near-independent noise.

-   Generation script provided in `data-raw/make_testdata.R`
    (`set.seed(7391)`); `data-raw/` excluded from the CRAN build via
    `.Rbuildignore`.

-   All runnable examples reference `scpairs_testdata` directly — no inline
    Seurat-object construction needed in any function documentation.

-   Dataset documented in `R/scpairs_testdata.R` with full `@format`,
    `@details`, `@usage`, `@seealso`, and `@family Section_0_Data` tags,
    following established CRAN package conventions.

------------------------------------------------------------------------

## Documentation Improvements

-   **Completed documentation** for `PlotPairSummary()` (added `@examples`,
    improved `@return` description) and `PlotPairSynergy()` (added
    `@examples` with `\dontrun{}`, expanded 4-panel `@return`).

-   **Added `@family` tags** grouping all exported functions into three
    help families:
    -   `Section_0_Data` — `scpairs_testdata`
    -   `Section_1_Discovery` — `FindAllPairs`, `FindGenePairs`,
        `AssessGenePair`
    -   `Section_2_Visualization` — all eleven `PlotPair*` and
        `PlotBridgeNetwork` functions

-   **Added `@seealso` cross-references** throughout, linking each function
    to the most relevant related functions.

-   **Improved `@return` descriptions** for `PlotPairDimplot`,
    `PlotPairViolin`, `PlotPairScatter`, `PlotPairSmoothed`,
    `PlotPairSummary`, `PlotPairCrossType`, `PlotPairSpatial`,
    `PlotPairNetwork`, and `PlotPairHeatmap`.

------------------------------------------------------------------------

## Package Structure

-   `R/data.R` renamed to `R/scpairs_testdata.R` (consistent with the
    convention of naming dataset documentation files after the dataset).

-   `tests/README.md` added: comprehensive test documentation covering test
    categories, expected behaviour, coverage targets, and a contribution
    template (modelled on SlimR).

-   Six demonstration plots saved to `docs/` from `scpairs_testdata`:
    `scpairs_network.png`, `scpairs_heatmap.png`, `scpairs_dimplot.png`,
    `scpairs_smoothed.png`, `scpairs_violin.png`, `scpairs_scatter.png`.

-   `README.md` overhauled with professional structure and CellJanus-inspired
    formatting: CRAN/GitHub badges, right-aligned logo, full TOC with anchor
    links, collapsible dependency sections, a **Visualization Gallery** section
    embedding the six demo images in a grid layout, and an expanded Functions
    reference table.

------------------------------------------------------------------------

## Test Suite

-   All seven test files updated to use `scpairs_testdata` directly,
    eliminating per-test Seurat pipeline execution and improving
    reproducibility.
    `create_test_seurat()` retained in `helper-test_data.R` for edge-case
    tests requiring custom configurations (single-cluster data, spatial
    coordinates, unusual gene counts).

-   **Fixed 2 test warnings** (`WARN 0` across all files):
    -   `test-assess-pair.R`: wrapped the `n_perm = 49` call in
        `suppressWarnings()` — the low-permutation advisory is intentional
        package behaviour, not a bug in the test.
    -   `test-consolidated-functions.R`: changed the third matrix row from
        `rnorm(n)` to `abs(rnorm(n))` in the `.compute_ratio_consistency`
        batch test, preventing `log2(mat + 1)` from receiving negative inputs.

-   **115 tests passing**, 0 failures, **0 warnings**, 0 skips.

------------------------------------------------------------------------

## Previous Versions

# scPairs 0.1.7 (2026-02-23)

## Bridge Gene Network Visualization

-   **New exported function `PlotBridgeNetwork()`**: standalone bridge gene
    network for a focal gene pair. Nodes positioned via MDS on Jaccard
    pathway-distance; radial depth ∝ shared-term count; focal genes at the
    centre. Bridge gene nodes sized by bridging strength and filled by mean
    expression (viridis). Focal → bridge edges coloured by source gene (red /
    blue), width ∝ shared term count; dotted inter-bridge edges encode pairwise
    Jaccard similarity. Supports `layout` parameter and configurable
    `sim_threshold`, `pt_size_range`, `edge_width_range`.

-   **Refactored `.plot_bridge_network_enhanced()`**: replaces prior static
    circular layout with MDS-radial positioning; adds weighted edges, 
    expression-coded fill, duplicate-edge deduplication, and bridge-bridge
    similarity overlay.

-   **`PlotPairSynergy()` panel 2** now calls `.plot_bridge_network_enhanced()`.

------------------------------------------------------------------------

# scPairs 0.1.6 (2026-02-11)

## Unified Computation Engine & API Consolidation

-   **Shared metric engine** (`compute_metrics.R`): all three discovery
    functions share a single `.compute_pair_metrics()` backend, eliminating
    \~600 lines of duplicated code.

-   **`mode` parameter**: `"all"` (default), `"expression"`, `"prior_only"`.

-   **Unified output format**: `AssessGenePair()` now includes a `$pairs`
    data.table with the same column schema as `FindAllPairs()` /
    `FindGenePairs()`. All visualisation functions accept any of the three
    result classes interchangeably.

-   **Removed `pheatmap`** from `Suggests` (unused).

-   **Internal refactoring**: `.resolve_cluster_ids()`, `.extract_pair_vals()`,
    `.extract_pairs_df()`.

------------------------------------------------------------------------

# scPairs 0.1.5 (2026-02-11)

## Prior Knowledge Integration & Synergy-Aware Scoring

-   **Prior knowledge layer**: `prior_score` (GO/KEGG Jaccard, weight 2.0) and
    `bridge_score` (pathway bridge genes, weight 1.8).

-   **New metric `neighbourhood_synergy`** (weight 1.5).

-   **New visualisation `PlotPairSynergy()`**: 6-panel figure.

-   **14 total metrics** across 5 evidence layers.

-   `use_prior`, `organism`, `custom_pairs` parameters added to all three main
    functions.

------------------------------------------------------------------------

# scPairs 0.1.4 (2026-02-08)

## Robustness & Validation

-   `min_pct_expressed` filter for sparse genes.
-   Comprehensive input validation with informative error messages.
-   Centralized schema (`schema.R`), shared plot utilities (`plot_utils.R`).

------------------------------------------------------------------------

# scPairs 0.1.3 (2026-02-08)

## Cross-Cell-Type Interaction Metric

-   **New metric `cross_celltype_score`** (trans-cellular synergy, weight 1.5).
-   **New visualisation `PlotPairCrossType()`**.

------------------------------------------------------------------------

# scPairs 0.1.2 (2026-02-08)

## Neighbourhood-Aware Metrics

-   3 new metrics: `smoothed_cor`, `neighbourhood_score`, `cluster_cor`.
-   New visualisations: `PlotPairSmoothed()`, `PlotPairSummary()`.

------------------------------------------------------------------------

# scPairs 0.1.1 (2026-02-08)

## Performance Optimizations

-   5–20× speedups through vectorisation of all core metric computations.

------------------------------------------------------------------------

# scPairs 0.1.0 (2026-02-07)

## Initial Release

-   3 core workflows: `FindAllPairs()`, `FindGenePairs()`, `AssessGenePair()`.
-   7 metrics: Pearson, Spearman, biweight, MI, ratio consistency, Lee's L, CLQ.
-   6 visualisation functions.
