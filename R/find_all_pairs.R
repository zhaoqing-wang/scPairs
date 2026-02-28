#' Discover All Synergistic Gene Pairs
#'
#' @description
#' The primary discovery function of **scPairs**.
#' Given a Seurat object, `FindAllPairs` identifies synergistic gene pairs by
#' integrating multiple lines of evidence: co-expression, neighbourhood
#' smoothing, prior biological knowledge, and spatial co-variation.
#'
#' @details
#' Metrics are rank-normalised and combined via weighted summation.
#' Optional permutation testing provides empirical p-values.
#'
#' The `mode` parameter controls which metric layers are computed:
#' * `"all"` (default) -- compute all available metrics.
#' * `"expression"` -- expression and neighbourhood metrics only (no prior
#'   knowledge).
#' * `"prior_only"` -- prior knowledge scores only (fast).
#'
#' @param object A Seurat object (scRNA-seq or spatial).
#' @param features Character vector of gene names to consider.
#'     NULL (default) uses Seurat `VariableFeatures`; if unavailable, selects
#'     the top `n_top_genes` by mean expression.
#' @param n_top_genes Integer; maximum number of genes to analyse when
#'     `features = NULL`.  Default 2000.
#' @param assay Character; assay to use.  Default: `DefaultAssay(object)`.
#' @param slot Character; data slot.  Default `"data"` (log-normalised).
#' @param cluster_col Character; column in `meta.data` with cluster IDs.
#'     NULL = use `Idents(object)`.
#' @param mode Character; `"all"`, `"expression"`, or `"prior_only"`.
#' @param cor_method Character vector; correlation methods to compute.
#'     Default `c("pearson", "spearman", "biweight")`.
#' @param n_mi_bins Integer; bins for mutual information.  0 = skip MI.
#' @param min_cells_expressed Integer; minimum co-expressing cells to keep a
#'     pair.  Default 10.
#' @param use_prior Logical; integrate prior knowledge (GO/KEGG).  Default TRUE.
#' @param organism Character; `"mouse"` or `"human"`.
#' @param custom_pairs Optional data.frame with columns `gene1`, `gene2`.
#' @param use_neighbourhood Logical; compute neighbourhood-aware metrics.
#'     Default TRUE.
#' @param neighbourhood_k Integer; KNN k.  Default 20.
#' @param neighbourhood_reduction Character; reduction for KNN.  Default `"pca"`.
#' @param smooth_alpha Numeric in \[0,1\]; self-weight for KNN smoothing.
#' @param use_spatial Logical; compute spatial metrics when available.
#' @param spatial_k Integer; spatial KNN k.
#' @param n_perm Integer; permutations for p-values.  0 = skip.
#' @param weights Named numeric; metric weights for score integration.
#' @param top_n Integer or NULL; return only top *n* pairs.
#' @param verbose Logical.
#'
#' @return A list with class `"scPairs_result"` containing:
#' \describe{
#'   \item{`pairs`}{`data.table` of gene pairs with all metric columns,
#'     `synergy_score`, `rank`, `p_value` (if permutation), `p_adj`,
#'     `confidence`.}
#'   \item{`parameters`}{List of analysis parameters.}
#'   \item{`n_genes`}{Number of genes analysed.}
#'   \item{`n_cells`}{Number of cells.}
#'   \item{`has_spatial`}{Logical.}
#'   \item{`mode`}{Character; the mode used.}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create a minimal Seurat object for demonstration
#' counts <- matrix(rpois(600, 5), nrow = 20, ncol = 30,
#'   dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:30)))
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'
#' result <- FindAllPairs(obj, n_top_genes = 20, top_n = 10,
#'   use_neighbourhood = FALSE, verbose = FALSE)
#' }
FindAllPairs <- function(object,
                         features              = NULL,
                         n_top_genes           = 2000,
                         assay                 = NULL,
                         slot                  = "data",
                         cluster_col           = NULL,
                         mode                  = c("all", "expression", "prior_only"),
                         cor_method            = c("pearson", "spearman", "biweight"),
                         n_mi_bins             = 5,
                         min_cells_expressed   = 10,
                         use_prior             = TRUE,
                         organism              = "mouse",
                         custom_pairs          = NULL,
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

  mode <- match.arg(mode)

  # --- Input validation ---
  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)
  n_cells_total <- ncol(object)

  if (mode != "prior_only") {
    .validate_cor_method(cor_method)
    .validate_min_cells_expressed(min_cells_expressed, n_cells_total)
    .validate_percentage(smooth_alpha, "smooth_alpha")
    if (use_neighbourhood) {
      .validate_neighbourhood_params(neighbourhood_k, n_cells_total)
    }
  }
  if (n_perm > 0) .validate_n_perm(n_perm)

  # Override flags based on mode
  if (mode == "expression")  use_prior <- FALSE
  if (mode == "prior_only") {
    use_neighbourhood <- FALSE
    use_spatial       <- FALSE
  }

  # --- Select features ---
  features <- .select_features(object, features, n_top = n_top_genes, assay = assay)
  if (!is.null(features) && length(features) > 0) {
    features <- .validate_features(features, object, assay)
  }
  .msg("Selected ", length(features), " genes for analysis.", verbose = verbose)

  # --- Extract data ---
  mat <- .get_expression_matrix(object, features = features,
                                assay = assay, slot = slot)

  cluster_ids <- .resolve_cluster_ids(object, cluster_col)

  # --- Build pair table via co-expression ---
  if (mode == "prior_only") {
    # Generate all unique pairs without co-expression filter
    n_genes <- length(features)
    pair_idx <- utils::combn(n_genes, 2)
    pair_dt <- data.table::data.table(
      gene1 = features[pair_idx[1, ]],
      gene2 = features[pair_idx[2, ]]
    )
  } else {
    pair_dt <- .compute_coexpression(
      mat                 = mat,
      features            = features,
      cluster_ids         = cluster_ids,
      cor_method          = cor_method,
      n_mi_bins           = n_mi_bins,
      min_cells_expressed = min_cells_expressed,
      verbose             = verbose
    )
  }

  if (nrow(pair_dt) == 0) {
    warning("No gene pairs passed the minimum co-expression filter.",
            call. = FALSE)
    return(.build_result(pair_dt, features, object, FALSE, list(), mode))
  }

  # --- Compute remaining metrics via shared engine ---
  # Skip expression metrics already computed by .compute_coexpression
  if (mode != "prior_only") {
    # Neighbourhood, prior, spatial, integration
    # Expression metrics already computed above; skip in engine
    engine <- .compute_pair_metrics(
      mat = mat, pair_dt = pair_dt, cluster_ids = cluster_ids,
      object = object, mode = mode,
      cor_method = character(0), n_mi_bins = 0,
      min_cells_expressed = 0,
      use_neighbourhood = use_neighbourhood,
      neighbourhood_k = neighbourhood_k,
      neighbourhood_reduction = neighbourhood_reduction,
      smooth_alpha = smooth_alpha,
      use_prior = use_prior, organism = organism,
      custom_pairs = custom_pairs,
      use_spatial = use_spatial, spatial_k = spatial_k,
      n_perm = n_perm, weights = weights, verbose = verbose
    )
  } else {
    engine <- .compute_pair_metrics(
      mat = mat, pair_dt = pair_dt, cluster_ids = cluster_ids,
      object = object, mode = "prior_only",
      use_prior = use_prior, organism = organism,
      custom_pairs = custom_pairs,
      n_perm = n_perm, weights = weights, verbose = verbose
    )
  }

  pair_dt <- engine$pair_dt

  # --- Top N filtering ---
  if (!is.null(top_n) && top_n < nrow(pair_dt)) {
    pair_dt <- pair_dt[seq_len(top_n), ]
  }

  # --- Return ---
  .build_result(pair_dt, features, object, engine$has_spatial,
                list(cor_method = cor_method, n_mi_bins = n_mi_bins,
                     min_cells_expressed = min_cells_expressed,
                     use_prior = use_prior, organism = organism,
                     use_neighbourhood = use_neighbourhood,
                     neighbourhood_k = neighbourhood_k,
                     smooth_alpha = smooth_alpha,
                     spatial_k = spatial_k, n_perm = n_perm,
                     weights = weights, top_n = top_n),
                mode)
}


#' Build standardised scPairs result object
#' @keywords internal
.build_result <- function(pair_dt, features, object, has_spatial, params,
                          mode = "all") {
  structure(
    list(
      pairs       = pair_dt,
      parameters  = params,
      n_genes     = length(features),
      n_cells     = ncol(object),
      has_spatial  = has_spatial,
      mode        = mode
    ),
    class = "scPairs_result"
  )
}


#' Resolve cluster IDs from Seurat object
#' @keywords internal
.resolve_cluster_ids <- function(object, cluster_col = NULL) {
  if (!is.null(cluster_col)) {
    if (!(cluster_col %in% colnames(object@meta.data))) {
      stop(sprintf("cluster_col '%s' not found in meta.data.", cluster_col),
           call. = FALSE)
    }
    as.factor(object@meta.data[[cluster_col]])
  } else {
    as.factor(Seurat::Idents(object))
  }
}


#' Print method for scPairs_result
#' @param x An scPairs_result object.
#' @param ... Ignored.
#' @return The input object \code{x}, returned invisibly.
#' @export
#' @method print scPairs_result
print.scPairs_result <- function(x, ...) {
  cat("scPairs result\n")
  cat(sprintf("  Genes analysed : %d\n", x$n_genes))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Pairs found    : %d\n", nrow(x$pairs)))
  cat(sprintf("  Spatial metrics: %s\n", ifelse(x$has_spatial, "Yes", "No")))
  if (!is.null(x$mode)) cat(sprintf("  Mode           : %s\n", x$mode))
  if (nrow(x$pairs) > 0 && "confidence" %in% colnames(x$pairs)) {
    conf_tbl <- table(x$pairs$confidence)
    cat("  Confidence: ",
        paste(names(conf_tbl), "=", conf_tbl, collapse = ", "), "\n")
  }
  invisible(x)
}
