# scPairs 0.1.8 (2026-02-28)

## Summary

First CRAN submission release. Addresses all reviewer feedback from the initial
check, introduces a built-in test dataset, completes function documentation,
and extends the test suite. Package structure has been overhauled to follow
established CRAN bioinformatics conventions.

------------------------------------------------------------------------

## CRAN Compliance

-   **Added reference URL** to the `Description` field of `DESCRIPTION`
    (`<https://github.com/zhaoqing-wang/scPairs>`), as required.

-   **Added `\value` tags** to all three exported `print()` S3 methods
    (`print.scPairs_result`, `print.scPairs_gene_result`,
    `print.scPairs_pair_result`); each returns its input object invisibly.

-   **Rationalised examples**: nine exported functions now run plain
    R examples during `R CMD check`. Two functions
    (`AssessGenePair`, `PlotPairSummary`) are wrapped in `\donttest{}`
    because their examples exceed the five-second threshold on the check
    servers. Three functions retain `\dontrun{}` because they require
    software not guaranteed to be present: `PlotPairSpatial` (spatial-assay
    Seurat object), `PlotBridgeNetwork` and `PlotPairSynergy` (Bioconductor
    annotation packages `org.Hs.eg.db` / `org.Mm.eg.db`).

-   **Removed `set.seed()`** from inside `.plot_bridge_network_enhanced()`
    (`R/plot_synergy.R`); the package no longer modifies the user's RNG
    state.

-   **Fixed Unicode in documentation**: replaced non-ASCII characters
    (U+2248 approximately-equal, U+00D7 multiplication sign, U+2013/U+2014
    dashes) in `scpairs_testdata` Rd documentation with ASCII equivalents,
    resolving the LaTeX PDF-manual build error.

-   **Fixed invalid URI in `README.md`**: the `LICENSE.md` hyperlink now
    points to the full GitHub URL rather than a bare file path, removing
    the URI validity note.

-   **Added `inst/WORDLIST`** listing domain-specific terms (Jaccard, KEGG,
    KNN, midcorrelation, transcriptomics, etc.) to suppress false-positive
    spell-check notes.

------------------------------------------------------------------------

## Built-in Test Dataset

-   **New dataset `scpairs_testdata`** (`data/scpairs_testdata.rda`, 26 KB):
    a synthetic Seurat object (100 cells, 20 genes, 3 balanced clusters,
    5-component PCA, 2-D UMAP) with two injected co-expression patterns:

    -   **GENE3 & GENE4**: globally correlated across all cells
        (Pearson r approximately 0.89 in normalised data).
    -   **GENE1 & GENE2**: moderately correlated within cluster 1 only.
    -   All remaining pairs: near-independent noise.

-   Generation script provided in `data-raw/make_testdata.R`
    (`set.seed(7391)`); `data-raw/` is excluded from the CRAN build via
    `.Rbuildignore`.

-   All runnable examples reference `scpairs_testdata` directly, eliminating
    the need for inline Seurat-object construction in any function
    documentation.

-   Dataset documented in `R/scpairs_testdata.R` with `@format`,
    `@details`, `@usage`, `@source`, `@seealso`, and
    `@family Section_0_Data` tags.

------------------------------------------------------------------------

## Documentation Improvements

-   **Completed documentation** for `PlotPairSummary()` (added `@examples`,
    improved `@return`) and `PlotPairSynergy()` (added `@examples` with
    `\dontrun{}`, expanded `@return` for all four panels).

-   **Added `@family` tags** grouping all exported functions into three
    help families:
    -   `Section_0_Data` -- `scpairs_testdata`
    -   `Section_1_Discovery` -- `FindAllPairs`, `FindGenePairs`,
        `AssessGenePair`
    -   `Section_2_Visualization` -- all eleven `PlotPair*` and
        `PlotBridgeNetwork` functions

-   **Added `@seealso` cross-references** throughout, linking each function
    to its most closely related counterparts.

-   **Improved `@return` descriptions** for `PlotPairDimplot`,
    `PlotPairViolin`, `PlotPairScatter`, `PlotPairSmoothed`,
    `PlotPairSummary`, `PlotPairCrossType`, `PlotPairSpatial`,
    `PlotPairNetwork`, and `PlotPairHeatmap`.

------------------------------------------------------------------------

## Package Structure

-   `R/data.R` renamed to `R/scpairs_testdata.R`, following the convention
    of naming dataset documentation files after the dataset itself.

-   `tests/README.md` added: documents test categories, expected behaviour,
    coverage targets, and a contribution template.

-   Six demonstration plots saved to `docs/` from `scpairs_testdata`:
    `scpairs_network.png`, `scpairs_heatmap.png`, `scpairs_dimplot.png`,
    `scpairs_smoothed.png`, `scpairs_violin.png`, `scpairs_scatter.png`.

-   `README.md` restructured: CRAN/GitHub badges, right-aligned logo,
    full table of contents, collapsible dependency section, Visualization
    Gallery embedding the six demo images, and an expanded Functions
    reference table.

------------------------------------------------------------------------

## Test Suite

-   All seven test files updated to use `scpairs_testdata` directly,
    removing per-test Seurat pipeline construction and improving
    reproducibility. `create_test_seurat()` is retained in
    `helper-test_data.R` for edge-case tests requiring custom
    configurations (single-cluster data, spatial coordinates, unusual
    gene counts).

-   **Two latent test warnings resolved**:
    -   `test-assess-pair.R`: the `n_perm = 49` call is now wrapped in
        `suppressWarnings()` -- the low-permutation advisory is intentional
        package behaviour, not a test defect.
    -   `test-consolidated-functions.R`: the third matrix row in the
        `.compute_ratio_consistency` batch test now uses `abs(rnorm(n))`
        to prevent `log2(mat + 1)` from receiving negative inputs.

-   **115 tests passing**, 0 failures, 0 warnings, 0 skips.

------------------------------------------------------------------------

# scPairs 0.1.7 (2026-02-23)

## Bridge Gene Network Visualization

-   **New exported function `PlotBridgeNetwork()`**: standalone bridge gene
    network for a focal gene pair. Nodes are positioned via MDS on Jaccard
    pathway distance; radial depth reflects shared-term count; focal genes
    appear at the centre. Bridge gene nodes are sized by bridging strength
    and filled by mean expression (viridis palette). Focal-to-bridge edges
    are coloured by source gene (red/blue) and scaled in width by shared
    term count; dotted inter-bridge edges encode pairwise Jaccard similarity.
    Supports a `layout` parameter and configurable `sim_threshold`,
    `pt_size_range`, and `edge_width_range`.

-   **Refactored `.plot_bridge_network_enhanced()`**: replaces the prior
    static circular layout with MDS-radial positioning; adds weighted edges,
    expression-coded fill, duplicate-edge deduplication, and bridge-bridge
    similarity overlay.

-   **`PlotPairSynergy()` panel 2** now calls
    `.plot_bridge_network_enhanced()`.

------------------------------------------------------------------------

# scPairs 0.1.6 (2026-02-11)

## Unified Computation Engine & API Consolidation

-   **Shared metric engine** (`compute_metrics.R`): all three discovery
    functions now share a single `.compute_pair_metrics()` backend,
    eliminating approximately 600 lines of duplicated code.

-   **`mode` parameter** added to all discovery functions:
    `"all"` (default), `"expression"`, `"prior_only"`.

-   **Unified output format**: `AssessGenePair()` now returns a `$pairs`
    data.table with the same column schema as `FindAllPairs()` and
    `FindGenePairs()`. All visualisation functions accept any of the three
    result classes interchangeably.

-   **Removed `pheatmap`** from `Suggests` (unused dependency).

-   **Internal helpers added**: `.resolve_cluster_ids()`,
    `.extract_pair_vals()`, `.extract_pairs_df()`.

------------------------------------------------------------------------

# scPairs 0.1.5 (2026-02-11)

## Prior Knowledge Integration & Synergy-Aware Scoring

-   **Prior knowledge layer**: `prior_score` (GO/KEGG Jaccard co-annotation,
    weight 2.0) and `bridge_score` (pathway bridge genes, weight 1.8).

-   **New metric `neighbourhood_synergy`** (weight 1.5).

-   **New visualisation `PlotPairSynergy()`**: 6-panel evidence dashboard.

-   Package now integrates **14 metrics** across **5 evidence layers**.

-   `use_prior`, `organism`, and `custom_pairs` parameters added to all
    three main discovery functions.

------------------------------------------------------------------------

# scPairs 0.1.4 (2026-02-08)

## Robustness & Input Validation

-   `min_pct_expressed` filter added to exclude sparse genes from analysis.
-   Comprehensive input validation with informative error messages.
-   Centralised schema (`schema.R`) and shared plot utilities
    (`plot_utils.R`).

------------------------------------------------------------------------

# scPairs 0.1.3 (2026-02-08)

## Cross-Cell-Type Interaction Metric

-   **New metric `cross_celltype_score`**: measures trans-cellular synergy
    between cell populations (weight 1.5).
-   **New visualisation `PlotPairCrossType()`**: cross-cell-type interaction
    heatmap.

------------------------------------------------------------------------

# scPairs 0.1.2 (2026-02-08)

## Neighbourhood-Aware Metrics

-   Three new metrics: `smoothed_cor`, `neighbourhood_score`, `cluster_cor`.
-   Two new visualisations: `PlotPairSmoothed()`, `PlotPairSummary()`.

------------------------------------------------------------------------

# scPairs 0.1.1 (2026-02-08)

## Performance Optimizations

-   5--20x speedups achieved through vectorisation of all core metric
    computations.

------------------------------------------------------------------------

# scPairs 0.1.0 (2026-02-07)

## Initial Release

-   Three core workflows: `FindAllPairs()`, `FindGenePairs()`,
    `AssessGenePair()`.
-   Seven metrics: Pearson, Spearman, biweight midcorrelation, mutual
    information, ratio consistency, Lee's L, and co-location quotient.
-   Six visualisation functions.
