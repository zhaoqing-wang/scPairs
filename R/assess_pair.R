#' Assess the Synergy of a Specific Gene Pair
#'
#' @description
#' Given two genes, `AssessGenePair` performs an in-depth evaluation of their
#' co-regulatory relationship.  In addition to the standard multi-metric
#' scoring, it computes:
#'
#' * **Per-cluster co-expression** -- correlation within each cell cluster,
#'   revealing whether synergy is global or cluster-specific.
#' * **Expression distribution overlap** -- Jaccard index of the sets of cells
#'   expressing each gene.
#' * **Permutation-based significance** -- by default 999 permutations are run
#'   to produce a robust p-value.
#'
#' For spatial data, spatial metrics (Lee's L, CLQ) are also computed.
#'
#' @details
#' **Performance (v0.1.1):** The permutation loop now pre-generates all
#' permuted y-vectors as a matrix and uses vectorised correlation, biweight,
#' MI, and ratio consistency computations.  Pearson and Spearman correlations
#' for all permutations are computed via a single matrix multiply.
#'
#' @param object A Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param assay Character; assay name.
#' @param slot Character; data slot.
#' @param cluster_col Character; cluster column.
#' @param use_neighbourhood Logical; compute neighbourhood-aware metrics
#'     (KNN-smoothed correlation and neighbourhood co-expression score).
#' @param neighbourhood_k Integer; number of nearest neighbours for the
#'     neighbourhood graph. Default 20.
#' @param neighbourhood_reduction Character; reduction to use for building the
#'     neighbourhood graph. Default "pca".
#' @param smooth_alpha Numeric in \[0,1\]; self-weight for KNN smoothing.
#'     0 = pure neighbour average, 1 = no smoothing. Default 0.3.
#' @param use_spatial Logical.
#' @param spatial_k Integer; KNN k.
#' @param n_perm Integer; number of permutations (default 999).
#' @param verbose Logical.
#'
#' @return A list with class `"scPairs_pair_result"`:
#' \describe{
#'   \item{`gene1`, `gene2`}{The query genes.}
#'   \item{`metrics`}{Named list of all computed metrics.}
#'   \item{`per_cluster`}{data.frame of per-cluster correlations.}
#'   \item{`synergy_score`}{Composite score.}
#'   \item{`p_value`}{Permutation p-value.}
#'   \item{`confidence`}{Categorical confidence label.}
#'   \item{`jaccard_index`}{Expression overlap Jaccard index.}
#'   \item{`has_spatial`}{Logical.}
#'   \item{`n_cells`}{Integer.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- AssessGenePair(seurat_obj, gene1 = "CD8A", gene2 = "CD8B")
#' print(result)
#'
#' # Visualise
#' PlotPairDimplot(seurat_obj, gene1 = "CD8A", gene2 = "CD8B")
#' PlotPairScatter(seurat_obj, gene1 = "CD8A", gene2 = "CD8B")
#' }
#'
AssessGenePair <- function(object,
                           gene1,
                           gene2,
                           assay                   = NULL,
                           slot                    = "data",
                           cluster_col             = NULL,
                           use_neighbourhood       = TRUE,
                           neighbourhood_k         = 20,
                           neighbourhood_reduction = "pca",
                           smooth_alpha            = 0.3,
                           use_spatial             = TRUE,
                           spatial_k               = 6,
                           n_perm                  = 999,
                           verbose                 = TRUE) {

  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)

  # Validate genes
  all_genes <- rownames(tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = slot)
  ))
  for (g in c(gene1, gene2)) {
    if (!(g %in% all_genes)) {
      stop(sprintf("Gene '%s' not found in the expression matrix.", g),
           call. = FALSE)
    }
  }
  if (gene1 == gene2) {
    stop("gene1 and gene2 must be different genes.", call. = FALSE)
  }

  mat <- .get_expression_matrix(object, features = c(gene1, gene2),
                                 assay = assay, slot = slot)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense <- as.matrix(mat)
  } else {
    mat_dense <- mat
  }

  x <- mat_dense[gene1, ]
  y <- mat_dense[gene2, ]
  n_cells <- length(x)

  .msg("Assessing pair: ", gene1, " -- ", gene2, " (", n_cells, " cells)",
       verbose = verbose)

  # Cluster info
  if (!is.null(cluster_col)) {
    cluster_ids <- as.factor(object@meta.data[[cluster_col]])
  } else {
    cluster_ids <- as.factor(Seurat::Idents(object))
  }

  # --- Global metrics -------------------------------------------------------
  metrics <- list()
  metrics$cor_pearson  <- stats::cor(x, y, method = "pearson")
  metrics$cor_spearman <- stats::cor(x, y, method = "spearman")
  metrics$cor_biweight <- .bicor(x, y)
  metrics$mi_score     <- .mutual_info(x, y, n_bins = 5)
  metrics$ratio_consistency <- .ratio_consistency(x, y, cluster_ids)

  # Correlation test p-values
  ct_pearson  <- stats::cor.test(x, y, method = "pearson")
  ct_spearman <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
  metrics$p_pearson  <- ct_pearson$p.value
  metrics$p_spearman <- ct_spearman$p.value

  # Expression overlap (Jaccard)
  expr_a <- x > 0
  expr_b <- y > 0
  intersection <- sum(expr_a & expr_b)
  union_ab     <- sum(expr_a | expr_b)
  jaccard <- if (union_ab > 0) intersection / union_ab else 0
  metrics$jaccard_index <- jaccard

  # Co-expression rate
  metrics$n_coexpressing <- intersection
  metrics$pct_coexpressing <- intersection / n_cells

  # --- Per-cluster correlations ---------------------------------------------
  clusters <- levels(cluster_ids)
  per_cluster <- data.frame(
    cluster    = clusters,
    n_cells    = integer(length(clusters)),
    cor_pearson  = numeric(length(clusters)),
    cor_spearman = numeric(length(clusters)),
    pct_coexpr   = numeric(length(clusters)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(clusters)) {
    idx <- which(cluster_ids == clusters[i])
    n_cl <- length(idx)
    per_cluster$n_cells[i] <- n_cl
    if (n_cl < 5) {
      per_cluster$cor_pearson[i]  <- NA
      per_cluster$cor_spearman[i] <- NA
      per_cluster$pct_coexpr[i]   <- NA
    } else {
      per_cluster$cor_pearson[i]  <- stats::cor(x[idx], y[idx], method = "pearson")
      per_cluster$cor_spearman[i] <- stats::cor(x[idx], y[idx], method = "spearman")
      per_cluster$pct_coexpr[i]   <- sum(x[idx] > 0 & y[idx] > 0) / n_cl
    }
  }

  # --- Neighbourhood metrics (v0.2.0) ----------------------------------------
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

      # KNN-smoothed correlation
      smoothed <- .smooth_expression(mat, W, alpha = smooth_alpha)
      if (inherits(smoothed, "sparseMatrix")) smoothed <- as.matrix(smoothed)
      sx <- smoothed[gene1, ]
      sy <- smoothed[gene2, ]
      metrics$smoothed_cor <- stats::cor(sx, sy, method = "pearson")

      # Neighbourhood co-expression score
      metrics$neighbourhood_score <- .neighbourhood_coexpr(x, y, W)

      # Cluster-level correlation
      metrics$cluster_cor <- .cluster_cor(x, y, cluster_ids)

      # Cross-cell-type interaction
      cross_ct <- .cross_celltype(x, y, W, cluster_ids)
      metrics$cross_celltype_score <- cross_ct$score
      metrics$cross_celltype_r_ab  <- cross_ct$r_ab
      metrics$cross_celltype_r_ba  <- cross_ct$r_ba
    }
  }

  # --- Spatial metrics -------------------------------------------------------
  has_spatial <- FALSE
  if (use_spatial && .has_spatial(object)) {
    has_spatial <- TRUE
    coords <- .get_spatial_coords(object)

    pair_dt <- data.table::data.table(gene1 = gene1, gene2 = gene2)
    pair_dt <- .compute_spatial_lee(coords, mat, pair_dt, k = spatial_k,
                                    n_perm = min(n_perm, 199), verbose = verbose)
    pair_dt <- .compute_spatial_clq(coords, mat, pair_dt, k = spatial_k,
                                     verbose = verbose)

    metrics$spatial_lee_L <- pair_dt$spatial_lee_L[1]
    if ("spatial_lee_p" %in% colnames(pair_dt)) {
      metrics$spatial_lee_p <- pair_dt$spatial_lee_p[1]
    }
    metrics$spatial_clq <- pair_dt$spatial_clq[1]
  }

  # --- Composite score via integration framework ----------------------------
  synergy <- mean(c(
    abs(metrics$cor_pearson),
    abs(metrics$cor_spearman),
    abs(metrics$cor_biweight),
    min(metrics$mi_score / log(5), 1),
    metrics$ratio_consistency,
    if (has_neighbourhood) abs(metrics$smoothed_cor) else NULL,
    if (has_neighbourhood) min(pmax(metrics$neighbourhood_score, 0) / 2, 1) else NULL,
    if (has_neighbourhood && !is.na(metrics$cluster_cor)) abs(metrics$cluster_cor) else NULL,
    if (has_neighbourhood && !is.na(metrics$cross_celltype_score)) metrics$cross_celltype_score else NULL,
    if (has_spatial) abs(metrics$spatial_lee_L) else NULL,
    if (has_spatial) min(metrics$spatial_clq / 2, 1) else NULL
  ), na.rm = TRUE)

  # --- Permutation p-value (vectorised) -------------------------------------
  p_value <- NA_real_
  if (n_perm > 0) {
    .msg("Permutation test (", n_perm, " permutations) ...", verbose = verbose)

    # Pre-generate all permuted y vectors as a matrix (n_cells × n_perm)
    perm_y_mat <- vapply(seq_len(n_perm), function(p) sample(y), numeric(n_cells))

    # Vectorised Pearson: cor(x, perm_y_mat) gives 1 × n_perm
    perm_pearson <- abs(as.numeric(stats::cor(x, perm_y_mat)))

    # Vectorised Spearman: rank then Pearson
    x_rank <- rank(x)
    perm_y_rank <- apply(perm_y_mat, 2, rank)
    perm_spearman <- abs(as.numeric(stats::cor(x_rank, perm_y_rank)))

    # Biweight, MI, ratio consistency: still per-perm but using fast helpers
    perm_bicor <- numeric(n_perm)
    perm_mi <- numeric(n_perm)
    perm_rc <- numeric(n_perm)
    perm_smoothed <- numeric(n_perm)
    perm_ncs <- numeric(n_perm)
    perm_clcor <- numeric(n_perm)
    perm_cross_ct <- numeric(n_perm)

    for (p in seq_len(n_perm)) {
      yp <- perm_y_mat[, p]
      perm_bicor[p] <- abs(.bicor(x, yp))
      perm_mi[p] <- min(.mutual_info(x, yp, 5) / log(5), 1)
      perm_rc[p] <- .ratio_consistency(x, yp, cluster_ids)
      if (has_neighbourhood) {
        perm_smoothed[p] <- abs(.smoothed_cor(
          rbind(x), rbind(yp), W, alpha = smooth_alpha))
        perm_ncs[p] <- min(pmax(.neighbourhood_coexpr(x, yp, W), 0) / 2, 1)
        perm_clcor[p] <- abs(.cluster_cor(x, yp, cluster_ids))
        ct_res <- .cross_celltype(x, yp, W, cluster_ids)
        perm_cross_ct[p] <- if (!is.na(ct_res$score)) ct_res$score else 0
      }
    }

    # Null composite scores (vectorised)
    n_metrics <- 5L
    null_scores <- perm_pearson + perm_spearman + perm_bicor +
                      perm_mi + perm_rc
    if (has_neighbourhood) {
      null_scores <- null_scores + perm_smoothed + perm_ncs + perm_clcor +
                        perm_cross_ct
      n_metrics <- n_metrics + 4L
    }
    null_scores <- null_scores / n_metrics

    p_value <- (sum(null_scores >= synergy) + 1) / (n_perm + 1)
  }

  confidence <- if (!is.na(p_value)) {
    if (p_value < 0.01) "High"
    else if (p_value < 0.05) "Medium"
    else if (p_value < 0.1) "Low"
    else "NS"
  } else {
    if (synergy >= 0.6) "High"
    else if (synergy >= 0.4) "Medium"
    else if (synergy >= 0.2) "Low"
    else "NS"
  }

  # --- Cross-cell-type detail ------------------------------------------------
  cross_celltype_detail <- if (has_neighbourhood && exists("cross_ct")) {
    cross_ct$per_celltype_pair
  } else {
    data.frame()
  }

  # --- Return ---------------------------------------------------------------
  structure(
    list(
      gene1                = gene1,
      gene2                = gene2,
      metrics              = metrics,
      per_cluster          = per_cluster,
      cross_celltype_detail = cross_celltype_detail,
      synergy_score        = synergy,
      p_value              = p_value,
      confidence           = confidence,
      jaccard_index        = jaccard,
      has_spatial           = has_spatial,
      n_cells              = n_cells
    ),
    class = "scPairs_pair_result"
  )
}


#' @export
#' @method print scPairs_pair_result
print.scPairs_pair_result <- function(x, ...) {
  cat("scPairs gene pair assessment\n")
  cat(sprintf("  Pair           : %s -- %s\n", x$gene1, x$gene2))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Synergy score  : %.4f\n", x$synergy_score))
  cat(sprintf("  Confidence     : %s\n", x$confidence))
  if (!is.na(x$p_value)) {
    cat(sprintf("  p-value        : %.4g\n", x$p_value))
  }
  cat(sprintf("  Jaccard overlap: %.3f\n", x$jaccard_index))
  cat("\n  Metrics:\n")
  for (nm in names(x$metrics)) {
    val <- x$metrics[[nm]]
    if (is.numeric(val)) {
      cat(sprintf("    %-20s: %.4f\n", nm, val))
    }
  }
  if (x$has_spatial) cat("  (Spatial metrics included)\n")
  cat("\n  Per-cluster correlations:\n")
  print(x$per_cluster, row.names = FALSE)
  if (!is.null(x$cross_celltype_detail) && nrow(x$cross_celltype_detail) > 0) {
    cat("\n  Cross-cell-type interactions:\n")
    print(x$cross_celltype_detail, row.names = FALSE)
  }
  invisible(x)
}
