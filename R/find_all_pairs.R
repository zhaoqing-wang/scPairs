#' Discover All Synergistic Gene Pairs in a Single-Cell or Spatial Dataset
#'
#' @description
#' The primary discovery function of **scPairs**.
#' Given a Seurat object (scRNA-seq or spatial transcriptomics), `FindAllPairs`
#' identifies statistically significant synergistic gene pairs by integrating
#' multiple lines of evidence:
#'
#' 1. **Co-expression** -- Pearson, Spearman, and biweight midcorrelation
#'    capture linear, rank-based, and robust associations.
#' 2. **Mutual information** -- detects non-linear dependencies missed by
#'    correlation.
#' 3. **Ratio consistency** -- tests whether the expression ratio of two genes
#'    is stable across cell clusters, a hallmark of genuine co-regulation.
#' 4. **Spatial co-variation** (spatial data only) -- Lee's L statistic measures
#'    bivariate spatial autocorrelation; the co-location quotient (CLQ) tests
#'    whether expressing cells are spatially attracted.
#'
#' Metrics are rank-normalised and combined via weighted summation.
#' Optional permutation testing provides empirical p-values.
#'
#' @details
#' **Performance (v0.1.1):** All core metric computations are vectorised.
#' The co-expression filter uses matrix crossproduct, biweight midcorrelation
#' is computed as a full matrix via vectorised kernel operations, and spatial
#' metrics use sparse matrix multiplication for the spatial lag.  These changes
#' yield 5--20x speedups on datasets with >500 genes.
#'
#' @param object A Seurat object (scRNA-seq or spatial).
#' @param features Character vector of gene names to consider.  NULL (default)
#'     uses Seurat `VariableFeatures`; if unavailable, selects the top
#'     `n_top_genes` by mean expression.
#' @param n_top_genes Integer; maximum number of genes to analyse when
#'     `features = NULL`.  Default 2000.
#' @param assay Character; assay to use.  Default: `DefaultAssay(object)`.
#' @param slot Character; data slot.  Default `"data"` (log-normalised).
#' @param cluster_col Character; column in `meta.data` with cluster IDs
#'     (for ratio-consistency).  NULL = use `Idents(object)`.
#' @param cor_method Character vector; correlation methods to compute.
#'     Default `c("pearson", "spearman", "biweight")`.
#' @param n_mi_bins Integer; bins for mutual information.  0 = skip MI.
#' @param min_cells_expressed Integer; minimum co-expressing cells to keep a
#'     pair. Default 10.  Set to 0 to retain all pairs (recommended when
#'     neighbourhood metrics are enabled).
#' @param use_neighbourhood Logical; compute neighbourhood-aware metrics
#'     (KNN-smoothed correlation, neighbourhood co-expression score, and
#'     cluster-level pseudo-bulk correlation).  Default TRUE.  These metrics
#'     are critical for detecting synergistic pairs where cell-level
#'     co-expression is absent due to sparsity.
#' @param neighbourhood_k Integer; number of nearest neighbours for the
#'     neighbourhood graph.  Default 20.
#' @param neighbourhood_reduction Character; reduction to use for building the
#'     KNN graph.  Default `"pca"`.
#' @param smooth_alpha Numeric in \[0,1\]; self-weight for KNN smoothing.
#'     0 = pure neighbour average, 1 = no smoothing.  Default 0.3.
#' @param use_spatial Logical; compute spatial metrics when available.
#'     Default TRUE.
#' @param spatial_k Integer; KNN neighbourhood size for spatial analyses.
#' @param n_perm Integer; permutations for p-values.  0 = skip (fast mode).
#' @param weights Named numeric vector; metric weights for score integration.
#'     NULL = use defaults.
#' @param top_n Integer or NULL; return only the top *n* pairs by synergy
#'     score.  NULL = return all.
#' @param verbose Logical; print progress.
#'
#' @return A list with class `"scPairs_result"` containing:
#' \describe{
#'   \item{`pairs`}{`data.table` of gene pairs with all metric columns,
#'     `synergy_score`, `rank`, `p_value` (if permutation), `p_adj`,
#'     `confidence`.}
#'   \item{`parameters`}{List of analysis parameters for reproducibility.}
#'   \item{`n_genes`}{Number of genes analysed.}
#'   \item{`n_cells`}{Number of cells.}
#'   \item{`has_spatial`}{Logical; whether spatial metrics were computed.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Discover top 100 synergistic pairs (fast, no permutation)
#' result <- FindAllPairs(seurat_obj, n_top_genes = 500, top_n = 100)
#' head(result$pairs)
#'
#' # Full analysis with permutation p-values
#' result <- FindAllPairs(seurat_obj, n_perm = 999)
#'
#' # Visualise the network
#' PlotPairNetwork(result)
#' }
#'
FindAllPairs <- function(object,
                         features              = NULL,
                         n_top_genes           = 2000,
                         assay                 = NULL,
                         slot                  = "data",
                         cluster_col           = NULL,
                         cor_method            = c("pearson", "spearman", "biweight"),
                         n_mi_bins             = 5,
                         min_cells_expressed   = 10,
                         use_neighbourhood     = TRUE,
                         neighbourhood_k       = 20,
                         neighbourhood_reduction = "pca",
                         smooth_alpha          = 0.3,
                         use_spatial           = TRUE,
                         spatial_k             = 6,
                         n_perm                = 0,
                         weights               = NULL,
                         top_n                 = NULL,
                         verbose               = TRUE) {

  # --- Input validation -----------------------------------------------------
  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)

  # --- Select features ------------------------------------------------------
  features <- .select_features(object, features, n_top = n_top_genes, assay = assay)
  .msg("Selected ", length(features), " genes for analysis.", verbose = verbose)

  # --- Extract data ---------------------------------------------------------
  mat <- .get_expression_matrix(object, features = features, assay = assay, slot = slot)

  # Cluster IDs
  if (!is.null(cluster_col)) {
    if (!(cluster_col %in% colnames(object@meta.data))) {
      stop(sprintf("cluster_col '%s' not found in meta.data.", cluster_col),
           call. = FALSE)
    }
    cluster_ids <- as.factor(object@meta.data[[cluster_col]])
  } else {
    cluster_ids <- as.factor(Seurat::Idents(object))
  }

  # --- Co-expression metrics ------------------------------------------------
  pair_dt <- .compute_coexpression(
    mat                 = mat,
    features            = features,
    cluster_ids         = cluster_ids,
    cor_method          = cor_method,
    n_mi_bins           = n_mi_bins,
    min_cells_expressed = min_cells_expressed,
    verbose             = verbose
  )

  if (nrow(pair_dt) == 0) {
    warning("No gene pairs passed the minimum co-expression filter.",
            call. = FALSE)
    return(.build_result(pair_dt, features, object, FALSE, list()))
  }

  # --- Neighbourhood metrics (v0.2.0) --------------------------------------
  has_neighbourhood <- FALSE
  if (use_neighbourhood) {
    W <- tryCatch(
      .build_knn_graph(object, reduction = neighbourhood_reduction,
                       k = neighbourhood_k),
      error = function(e) {
        .msg("Could not build KNN graph: ", conditionMessage(e),
             ". Skipping neighbourhood metrics.", verbose = verbose)
        NULL
      }
    )

    if (!is.null(W)) {
      has_neighbourhood <- TRUE
      .msg("Computing neighbourhood-aware metrics ...", verbose = verbose)

      .msg("  KNN-smoothed correlation ...", verbose = verbose)
      pair_dt[, smoothed_cor := .smoothed_cor_batch(
        mat, pair_dt, W, alpha = smooth_alpha, method = "pearson")]

      .msg("  Neighbourhood co-expression score ...", verbose = verbose)
      pair_dt[, neighbourhood_score := .neighbourhood_coexpr_batch(
        mat, pair_dt, W)]

      .msg("  Cluster-level correlation ...", verbose = verbose)
      pair_dt[, cluster_cor := .cluster_cor_batch(mat, cluster_ids, pair_dt)]

      .msg("  Cross-cell-type interaction score ...", verbose = verbose)
      pair_dt[, cross_celltype_score := .cross_celltype_batch(
        mat, pair_dt, W, cluster_ids)]
    }
  }

  # --- Spatial metrics (if available) ---------------------------------------
  has_spatial <- FALSE
  if (use_spatial && .has_spatial(object)) {
    has_spatial <- TRUE
    .msg("Spatial data detected. Computing spatial metrics ...", verbose = verbose)
    coords <- .get_spatial_coords(object)

    pair_dt <- .compute_spatial_lee(
      coords  = coords,
      mat     = mat,
      pair_dt = pair_dt,
      k       = spatial_k,
      n_perm  = min(n_perm, 199),
      verbose = verbose
    )

    pair_dt <- .compute_spatial_clq(
      coords         = coords,
      mat            = mat,
      pair_dt        = pair_dt,
      k              = spatial_k,
      expr_threshold = 0,
      verbose        = verbose
    )
  }

  # --- Score integration ----------------------------------------------------
  pair_dt <- .integrate_scores(
    pair_dt     = pair_dt,
    weights     = weights,
    n_perm      = n_perm,
    mat         = mat,
    cluster_ids = cluster_ids,
    coords      = if (has_spatial) coords else NULL,
    spatial_k   = spatial_k,
    verbose     = verbose
  )

  # --- Top N filtering ------------------------------------------------------
  data.table::setorder(pair_dt, rank)
  if (!is.null(top_n) && top_n < nrow(pair_dt)) {
    pair_dt <- pair_dt[seq_len(top_n), ]
  }

  # --- Return ---------------------------------------------------------------
  .build_result(pair_dt, features, object, has_spatial,
                list(cor_method = cor_method, n_mi_bins = n_mi_bins,
                     min_cells_expressed = min_cells_expressed,
                     use_neighbourhood = use_neighbourhood,
                     neighbourhood_k = neighbourhood_k,
                     smooth_alpha = smooth_alpha,
                     spatial_k = spatial_k, n_perm = n_perm,
                     weights = weights, top_n = top_n))
}


#' Build standardised scPairs result object
#' @keywords internal
.build_result <- function(pair_dt, features, object, has_spatial, params) {
  structure(
    list(
      pairs       = pair_dt,
      parameters  = params,
      n_genes     = length(features),
      n_cells     = ncol(object),
      has_spatial  = has_spatial
    ),
    class = "scPairs_result"
  )
}


#' Print method for scPairs_result
#' @param x An scPairs_result object.
#' @param ... Ignored.
#' @export
#' @method print scPairs_result
print.scPairs_result <- function(x, ...) {
  cat("scPairs result\n")
  cat(sprintf("  Genes analysed : %d\n", x$n_genes))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Pairs found    : %d\n", nrow(x$pairs)))
  cat(sprintf("  Spatial metrics: %s\n", ifelse(x$has_spatial, "Yes", "No")))
  if (nrow(x$pairs) > 0 && "confidence" %in% colnames(x$pairs)) {
    conf_tbl <- table(x$pairs$confidence)
    cat("  Confidence: ",
        paste(names(conf_tbl), "=", conf_tbl, collapse = ", "), "\n")
  }
  invisible(x)
}
