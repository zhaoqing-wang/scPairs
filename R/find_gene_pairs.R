#' Find Synergistic Partners for a Given Gene
#'
#' @description
#' Given a gene of interest, `FindGenePairs` identifies and ranks all genes
#' that act synergistically with it.  The function uses the same multi-evidence
#' framework as [FindAllPairs()] but focuses computation on pairs involving the
#' query gene, making it much faster for targeted queries.
#'
#' @details
#' **Performance (v0.1.1):** Pearson and Spearman correlations for the query
#' gene against all candidates are computed via a single vectorised matrix
#' multiply instead of per-candidate loops.  The co-expression filter uses
#' a vectorised matrix operation.  Biweight midcorrelation, mutual information,
#' and ratio consistency are computed using batch-optimised helpers where
#' possible.
#'
#' @param object A Seurat object.
#' @param gene Character; the query gene name.
#' @param candidates Character vector of candidate partner genes.  NULL =
#'     auto-select from `VariableFeatures` or top-expressed genes.
#' @param n_top_genes Integer; max candidates when `candidates = NULL`.
#' @param assay Character; assay name.
#' @param slot Character; data slot.
#' @param cluster_col Character; cluster column in meta.data (NULL = Idents).
#' @param cor_method Correlation methods to use.
#' @param n_mi_bins Bins for mutual information.
#' @param min_cells_expressed Minimum cells co-expressing both genes.
#' @param use_neighbourhood Logical; compute neighbourhood-aware metrics
#'     (KNN-smoothed correlation and neighbourhood co-expression score).
#' @param neighbourhood_k Integer; number of nearest neighbours for the
#'     neighbourhood graph. Default 20.
#' @param neighbourhood_reduction Character; reduction to use for building the
#'     neighbourhood graph. Default "pca".
#' @param smooth_alpha Numeric in \[0,1\]; self-weight for KNN smoothing.
#'     0 = pure neighbour average, 1 = no smoothing. Default 0.3.
#' @param use_spatial Logical; compute spatial metrics if available.
#' @param spatial_k Integer; neighbourhood size for spatial metrics.
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
#'   \item{`n_candidates`}{Number of candidate genes tested.}
#'   \item{`n_cells`}{Number of cells.}
#'   \item{`has_spatial`}{Logical.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Find top 20 partners for TP53
#' result <- FindGenePairs(seurat_obj, gene = "TP53", top_n = 20)
#' head(result$pairs)
#' PlotPairNetwork(result)
#' }
#'
FindGenePairs <- function(object,
                          gene,
                          candidates              = NULL,
                          n_top_genes             = 2000,
                          assay                   = NULL,
                          slot                    = "data",
                          cluster_col             = NULL,
                          cor_method              = c("pearson", "spearman", "biweight"),
                          n_mi_bins               = 5,
                          min_cells_expressed     = 10,
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

  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)

  # Validate query gene exists
  all_genes <- rownames(tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = slot)
  ))
  if (!(gene %in% all_genes)) {
    stop(sprintf("Query gene '%s' not found in the expression matrix.", gene),
         call. = FALSE)
  }

  # Select candidates
  candidates <- .select_features(object, candidates, n_top = n_top_genes, assay = assay)
  candidates <- setdiff(candidates, gene)
  if (length(candidates) == 0) {
    stop("No candidate partner genes available.", call. = FALSE)
  }

  features <- c(gene, candidates)
  .msg("Querying partners for ", gene, " against ", length(candidates),
       " candidates.", verbose = verbose)

  # Extract data
  mat <- .get_expression_matrix(object, features = features, assay = assay, slot = slot)

  # Cluster IDs
  if (!is.null(cluster_col)) {
    cluster_ids <- as.factor(object@meta.data[[cluster_col]])
  } else {
    cluster_ids <- as.factor(Seurat::Idents(object))
  }

  # ------------------------------------------------------------------
  # Build pair table: gene vs. each candidate
  # ------------------------------------------------------------------
  pair_dt <- data.table::data.table(
    gene1 = rep(gene, length(candidates)),
    gene2 = candidates
  )

  # Convert to dense for computation
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense <- as.matrix(mat)
  } else {
    mat_dense <- mat
  }

  n_cells <- ncol(mat_dense)
  gene_vec <- mat_dense[gene, ]

  # --- Filter by co-expression (vectorised) ---
  if (min_cells_expressed > 0) {
    gene_nonzero <- gene_vec > 0
    # Count co-expression for all candidates at once
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
      return(.build_gene_result(pair_dt, gene, candidates, object, FALSE, list()))
    }
    cand_mat <- cand_mat[candidates, , drop = FALSE]
  } else {
    cand_mat <- mat_dense[candidates, , drop = FALSE]
  }

  n_cand <- length(candidates)

  # --- Pearson (vectorised: cor(gene_vec, cand_mat_t)) ---
  if ("pearson" %in% cor_method) {
    .msg("  Pearson correlation ...", verbose = verbose)
    cor_vals <- as.numeric(stats::cor(gene_vec, t(cand_mat)))
    pair_dt[, cor_pearson := cor_vals]
  }

  # --- Spearman (vectorised via ranks) ---
  if ("spearman" %in% cor_method) {
    .msg("  Spearman correlation ...", verbose = verbose)
    gene_rank <- rank(gene_vec)
    cand_ranks <- t(apply(cand_mat, 1, rank))
    cor_sp <- as.numeric(stats::cor(gene_rank, t(cand_ranks)))
    pair_dt[, cor_spearman := cor_sp]
  }

  # --- Biweight (per-candidate, using fast .bicor) ---
  if ("biweight" %in% cor_method) {
    .msg("  Biweight midcorrelation ...", verbose = verbose)
    pair_dt[, cor_biweight := vapply(seq_len(n_cand), function(k) {
      .bicor(gene_vec, cand_mat[k, ])
    }, numeric(1))]
  }

  # --- Mutual information ---
  if (n_mi_bins > 0) {
    .msg("  Mutual information ...", verbose = verbose)
    pair_dt[, mi_score := vapply(seq_len(n_cand), function(k) {
      .mutual_info(gene_vec, cand_mat[k, ], n_bins = n_mi_bins)
    }, numeric(1))]
  }

  # --- Ratio consistency ---
  if (!is.null(cluster_ids)) {
    .msg("  Ratio consistency ...", verbose = verbose)
    pair_dt[, ratio_consistency := vapply(seq_len(n_cand), function(k) {
      .ratio_consistency(gene_vec, cand_mat[k, ], cluster_ids)
    }, numeric(1))]
  }

  # --- Neighbourhood metrics (v0.2.0) ---
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
      .msg("  KNN-smoothed correlation ...", verbose = verbose)
      pair_dt[, smoothed_cor := .smoothed_cor_batch(
        mat, pair_dt, W, alpha = smooth_alpha, method = "pearson")]

      .msg("  Neighbourhood co-expression score ...", verbose = verbose)
      pair_dt[, neighbourhood_score := .neighbourhood_coexpr_batch(
        mat, pair_dt, W)]

      .msg("  Cluster-level correlation ...", verbose = verbose)
      pair_dt[, cluster_cor := .cluster_cor_batch(mat, cluster_ids, pair_dt)]

      .msg("  Cross-cell-type interaction score ...", verbose = verbose)
      embed <- tryCatch(
        Seurat::Embeddings(object, reduction = neighbourhood_reduction),
        error = function(e) NULL
      )
      if (!is.null(embed)) {
        dims_use <- seq_len(min(30, ncol(embed)))
        pair_dt[, cross_celltype_score := .cross_celltype_batch(
          mat, pair_dt, cluster_ids, embed[, dims_use, drop = FALSE])]
      }
    }
  }

  # --- Spatial metrics ---
  has_spatial <- FALSE
  if (use_spatial && .has_spatial(object)) {
    has_spatial <- TRUE
    coords <- .get_spatial_coords(object)
    pair_dt <- .compute_spatial_lee(coords, mat, pair_dt, k = spatial_k,
                                    n_perm = min(n_perm, 199), verbose = verbose)
    pair_dt <- .compute_spatial_clq(coords, mat, pair_dt, k = spatial_k,
                                     verbose = verbose)
  }

  # --- Score integration ---
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

  data.table::setorder(pair_dt, rank)
  if (!is.null(top_n) && top_n < nrow(pair_dt)) {
    pair_dt <- pair_dt[seq_len(top_n), ]
  }

  .build_gene_result(pair_dt, gene, candidates, object, has_spatial,
                     list(cor_method = cor_method, n_mi_bins = n_mi_bins,
                          min_cells_expressed = min_cells_expressed,
                          spatial_k = spatial_k, n_perm = n_perm))
}


#' @keywords internal
.build_gene_result <- function(pair_dt, gene, candidates, object, has_spatial, params) {
  structure(
    list(
      query_gene   = gene,
      pairs        = pair_dt,
      parameters   = params,
      n_candidates = length(candidates),
      n_cells      = ncol(object),
      has_spatial   = has_spatial
    ),
    class = "scPairs_gene_result"
  )
}


#' @export
#' @method print scPairs_gene_result
print.scPairs_gene_result <- function(x, ...) {
  cat("scPairs gene query result\n")
  cat(sprintf("  Query gene     : %s\n", x$query_gene))
  cat(sprintf("  Candidates     : %d\n", x$n_candidates))
  cat(sprintf("  Partners found : %d\n", nrow(x$pairs)))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Spatial metrics: %s\n", ifelse(x$has_spatial, "Yes", "No")))
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
