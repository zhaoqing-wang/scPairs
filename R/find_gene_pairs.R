#' Find Synergistic Partners for a Given Gene
#'
#' @description
#' Given a gene of interest, `FindGenePairs` identifies and ranks all genes
#' that act synergistically with it.  Uses the same multi-evidence framework
#' as [FindAllPairs()] but focuses computation on pairs involving the query
#' gene, making it much faster for targeted queries.
#'
#' @details
#' The `mode` parameter controls which metric layers are computed:
#' * `"all"` (default) -- all available metrics.
#' * `"expression"` -- expression and neighbourhood only.
#' * `"prior_only"` -- prior knowledge scores only.
#'
#' @param object A Seurat object.
#' @param gene Character; the query gene name.
#' @param candidates Character vector of candidate partner genes.
#'     NULL = auto-select.
#' @param n_top_genes Integer; max candidates when `candidates = NULL`.
#' @param assay Character; assay name.
#' @param slot Character; data slot.
#' @param cluster_col Character; cluster column in meta.data.
#' @param mode Character; `"all"`, `"expression"`, or `"prior_only"`.
#' @param cor_method Correlation methods.
#' @param n_mi_bins Bins for mutual information.
#' @param min_cells_expressed Minimum cells co-expressing both genes.
#' @param use_prior Logical; compute prior knowledge metrics.
#' @param organism Character; `"mouse"` or `"human"`.
#' @param custom_pairs Optional data.frame of custom interactions.
#' @param use_neighbourhood Logical; compute neighbourhood metrics.
#' @param neighbourhood_k Integer; KNN k.
#' @param neighbourhood_reduction Character; reduction for KNN.
#' @param smooth_alpha Numeric; self-weight for smoothing.
#' @param use_spatial Logical; compute spatial metrics.
#' @param spatial_k Integer; spatial neighbourhood k.
#' @param n_perm Integer; permutations for p-values.
#' @param weights Named numeric; metric weights.
#' @param top_n Integer; return only top partners.
#' @param verbose Logical.
#'
#' @return A list with class `"scPairs_gene_result"`:
#' \describe{
#'   \item{`query_gene`}{The input gene.}
#'   \item{`pairs`}{`data.table` of partners ranked by synergy score.}
#'   \item{`parameters`}{Analysis parameters.}
#'   \item{`n_candidates`}{Number of candidates tested.}
#'   \item{`n_cells`}{Number of cells.}
#'   \item{`has_spatial`}{Logical.}
#'   \item{`mode`}{Character.}
#' }
#'
#' @family Section_1_Discovery
#'
#' @seealso \code{\link{FindAllPairs}} for genome-wide screening,
#'   \code{\link{AssessGenePair}} for in-depth single-pair assessment,
#'   \code{\link{PlotPairNetwork}} for visualising partner networks.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Find synergistic partners of GENE3.  GENE4 is expected to rank first.
#' result <- FindGenePairs(scpairs_testdata,
#'                         gene    = "GENE3",
#'                         top_n   = 10,
#'                         mode    = "expression",
#'                         verbose = FALSE)
#' print(result)
#' }
FindGenePairs <- function(object,
                          gene,
                          candidates              = NULL,
                          n_top_genes             = 2000,
                          assay                   = NULL,
                          slot                    = "data",
                          cluster_col             = NULL,
                          mode                    = c("all", "expression", "prior_only"),
                          cor_method              = c("pearson", "spearman", "biweight"),
                          n_mi_bins               = 5,
                          min_cells_expressed     = 10,
                          use_prior               = TRUE,
                          organism                = "mouse",
                          custom_pairs            = NULL,
                          use_neighbourhood       = TRUE,
                          neighbourhood_k         = 20,
                          neighbourhood_reduction = "pca",
                          smooth_alpha            = 0.3,
                          use_spatial             = TRUE,
                          spatial_k               = 6,
                          n_perm                  = 0,
                          weights                 = NULL,
                          top_n                   = NULL,
                          verbose                 = TRUE) {

  mode <- match.arg(mode)

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

  if (mode == "expression")  use_prior <- FALSE
  if (mode == "prior_only") {
    use_neighbourhood <- FALSE
    use_spatial       <- FALSE
  }

  .validate_features(gene, object, assay)

  # --- Select candidates ---
  candidates <- .select_features(object, candidates, n_top = n_top_genes,
                                 assay = assay)
  candidates <- setdiff(candidates, gene)
  if (length(candidates) == 0) {
    stop("No candidate partner genes available.", call. = FALSE)
  }

  features <- c(gene, candidates)
  .msg("Querying partners for ", gene, " against ", length(candidates),
       " candidates.", verbose = verbose)

  mat <- .get_expression_matrix(object, features = features,
                                assay = assay, slot = slot)
  cluster_ids <- .resolve_cluster_ids(object, cluster_col)

  # --- Build pair table ---
  pair_dt <- data.table::data.table(
    gene1 = rep(gene, length(candidates)),
    gene2 = candidates
  )

  # --- Co-expression filter (vectorised) ---
  if (mode != "prior_only" && min_cells_expressed > 0) {
    if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
      mat_dense <- as.matrix(mat)
    } else {
      mat_dense <- mat
    }
    gene_vec <- mat_dense[gene, ]
    gene_nonzero <- gene_vec > 0
    cand_mat <- mat_dense[candidates, , drop = FALSE]
    cand_nonzero <- cand_mat > 0
    storage.mode(cand_nonzero) <- "double"
    storage.mode(gene_nonzero) <- "double"
    co_expr <- as.numeric(cand_nonzero %*% gene_nonzero)
    keep <- co_expr >= min_cells_expressed
    pair_dt <- pair_dt[keep, ]
    candidates <- candidates[keep]
    .msg("  Retained ", sum(keep), " candidates after co-expression filter.",
         verbose = verbose)
    if (nrow(pair_dt) == 0) {
      warning("No pairs passed the co-expression filter.", call. = FALSE)
      return(.build_gene_result(pair_dt, gene, candidates, object, FALSE,
                                list(), mode))
    }
  }

  # --- Vectorised correlations for query gene ---
  if (mode != "prior_only" && nrow(pair_dt) > 0) {
    if (!exists("mat_dense", inherits = FALSE)) {
      mat_dense <- if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix"))
        as.matrix(mat) else mat
    }
    gene_vec <- mat_dense[gene, ]
    cand_mat <- mat_dense[candidates, , drop = FALSE]
    n_cand <- length(candidates)

    if ("pearson" %in% cor_method) {
      .msg("  Pearson correlation ...", verbose = verbose)
      pair_dt[, cor_pearson := as.numeric(stats::cor(gene_vec, t(cand_mat)))]
    }
    if ("spearman" %in% cor_method) {
      .msg("  Spearman correlation ...", verbose = verbose)
      gene_rank <- rank(gene_vec)
      cand_ranks <- t(apply(cand_mat, 1, rank))
      pair_dt[, cor_spearman := as.numeric(stats::cor(gene_rank, t(cand_ranks)))]
    }
    if ("biweight" %in% cor_method) {
      .msg("  Biweight midcorrelation ...", verbose = verbose)
      pair_dt[, cor_biweight := vapply(seq_len(n_cand), function(k) {
        .bicor(gene_vec, cand_mat[k, ])
      }, numeric(1))]
    }
    if (n_mi_bins > 0) {
      .msg("  Mutual information ...", verbose = verbose)
      pair_dt[, mi_score := vapply(seq_len(n_cand), function(k) {
        .mutual_info(gene_vec, cand_mat[k, ], n_bins = n_mi_bins)
      }, numeric(1))]
    }
    if (!is.null(cluster_ids)) {
      .msg("  Ratio consistency ...", verbose = verbose)
      pair_dt[, ratio_consistency := vapply(seq_len(n_cand), function(k) {
        .ratio_consistency(gene_vec, cand_mat[k, ], cluster_ids)
      }, numeric(1))]
    }
  }

  # --- Remaining metrics via shared engine (expression already computed) ---
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

  pair_dt <- engine$pair_dt

  if (!is.null(top_n) && top_n < nrow(pair_dt)) {
    pair_dt <- pair_dt[seq_len(top_n), ]
  }

  .build_gene_result(pair_dt, gene, candidates, object, engine$has_spatial,
                     list(cor_method = cor_method, n_mi_bins = n_mi_bins,
                          min_cells_expressed = min_cells_expressed,
                          spatial_k = spatial_k, n_perm = n_perm),
                     mode)
}


#' @keywords internal
.build_gene_result <- function(pair_dt, gene, candidates, object,
                               has_spatial, params, mode = "all") {
  structure(
    list(
      query_gene   = gene,
      pairs        = pair_dt,
      parameters   = params,
      n_candidates = length(candidates),
      n_cells      = ncol(object),
      has_spatial   = has_spatial,
      mode         = mode
    ),
    class = "scPairs_gene_result"
  )
}


#' Print method for scPairs_gene_result
#' @param x An scPairs_gene_result object.
#' @param ... Ignored.
#' @return The input object \code{x}, returned invisibly.
#' @export
#' @method print scPairs_gene_result
print.scPairs_gene_result <- function(x, ...) {
  cat("scPairs gene query result\n")
  cat(sprintf("  Query gene     : %s\n", x$query_gene))
  cat(sprintf("  Candidates     : %d\n", x$n_candidates))
  cat(sprintf("  Partners found : %d\n", nrow(x$pairs)))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Spatial metrics: %s\n", ifelse(x$has_spatial, "Yes", "No")))
  if (!is.null(x$mode)) cat(sprintf("  Mode           : %s\n", x$mode))
  if (nrow(x$pairs) > 0) {
    cat("  Top 5 partners:\n")
    top5 <- utils::head(x$pairs, 5)
    for (i in seq_len(nrow(top5))) {
      cat(sprintf("    %d. %s (score = %.3f)\n",
                  i, top5$gene2[i], top5$synergy_score[i]))
    }
  }
  invisible(x)
}
