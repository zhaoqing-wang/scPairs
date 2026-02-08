#' Cross-cell-type interaction metrics
#'
#' Standard co-expression metrics measure whether two genes are expressed
#' together in the **same** cell.  However, many biologically important
#' interactions operate **across** cell types -- gene A expressed in one cell
#' type signals to (or synergises with) gene B expressed in a neighbouring cell
#' of a different type.  Classic examples include Adora2a--Ido1 trans-cellular
#' signalling and ligand--receptor pairs.
#'
#' This module provides a **cross-cell-type interaction score** that detects
#' such trans-cellular synergies.  The key insight is that "neighbouring" in
#' this context means cells that could plausibly interact -- **not** cells that
#' are close in PCA/UMAP space (which groups cells of the *same* type
#' together).
#'
#' **Algorithm:**
#' For each ordered pair of cell types (type_i, type_j where i != j):
#'
#' \enumerate{
#'   \item Partition cells of each type into \code{n_bins} matched groups
#'         (bins) by splitting the shared PCA embedding into spatial
#'         micro-environments, so that cells in the same bin represent the
#'         same local tissue context.
#'   \item Compute pseudo-bulk expression of gene A in type_i cells per bin,
#'         and pseudo-bulk expression of gene B in type_j cells per bin.
#'   \item Correlate these paired pseudo-bulks across bins.
#' }
#'
#' This captures whether tissue regions where type_i cells express high
#' gene A tend to be regions where type_j cells express high gene B --
#' exactly the signal expected from paracrine signalling or trans-cellular
#' co-regulation.
#'
#' @name cross-celltype-metrics
#' @keywords internal
NULL


# ===========================================================================
#  Helper: Build micro-environment bins from PCA embedding
# ===========================================================================

#' Assign cells to micro-environment bins using k-means on PCA embedding
#'
#' Cells are clustered into spatial micro-environment bins so that each bin
#' represents a local neighbourhood of the tissue.  Crucially, cells of
#' *all* types within the same bin are considered to share a micro-environment
#' and can thus interact.
#'
#' @param embed Numeric matrix (n_cells x n_dims) of PCA (or other) embedding.
#' @param n_bins Integer; number of micro-environment bins.
#' @return Integer vector of bin assignments (1..n_bins).
#' @keywords internal
.assign_microenv_bins <- function(embed, n_bins) {
  n <- nrow(embed)
  n_bins <- min(n_bins, n)

  km <- tryCatch(
    stats::kmeans(embed, centers = n_bins, nstart = 3, iter.max = 50),
    error = function(e) NULL
  )

  if (is.null(km)) {
    # Fallback: assign bins by cutting each PCA dimension
    scores <- embed[, 1]
    return(as.integer(cut(scores, breaks = n_bins, labels = FALSE)))
  }

  km$cluster
}


# ===========================================================================
#  1. Cross-cell-type interaction score (batch)
# ===========================================================================

#' Compute cross-cell-type interaction scores for gene pairs
#'
#' For each gene pair (A, B), this metric measures whether expression of
#' gene A in cells of one type correlates with expression of gene B in
#' cells of a *different* type that share the same tissue micro-environment.
#'
#' @section Algorithm:
#'
#' 1. Partition all cells into micro-environment bins using k-means on the
#'    PCA embedding.  Each bin represents a local tissue context containing
#'    cells of multiple types.
#' 2. For each cell-type pair (type_i, type_j, i != j), compute pseudo-bulk
#'    expression per bin: mean(gene_A) in type_i cells of that bin and
#'    mean(gene_B) in type_j cells of that bin.
#' 3. Correlate the paired pseudo-bulk vectors across bins (Pearson r).
#'    This gives the directed score A -> B.
#' 4. Repeat for B -> A.
#' 5. Aggregate across all cell-type pairs using a weighted mean (weighted
#'    by number of bins with sufficient cells of both types).
#' 6. Final score = geometric mean of |aggregated r(A->B)| and
#'    |aggregated r(B->A)|.
#'
#' @param mat Expression matrix (genes x cells, dense or sparse).
#' @param pair_dt \code{data.table} with columns \code{gene1}, \code{gene2}.
#' @param cluster_ids Factor of cell-type / cluster assignments.
#' @param embed Numeric matrix of PCA embedding (n_cells x n_dims).
#' @param n_bins Integer; number of micro-environment bins.  Default 50.
#' @param min_cells_per_bin Integer; minimum cells of a given type in a bin
#'   for that bin to contribute to the correlation.  Default 5.
#' @param min_bins Integer; minimum bins with both types present to compute
#'   a correlation for a given cell-type pair.  Default 8.
#'
#' @return Numeric vector of cross-cell-type interaction scores (one per row
#'   of \code{pair_dt}).
#'
#' @keywords internal
.cross_celltype_batch <- function(mat, pair_dt, cluster_ids, embed,
                                  n_bins            = 50,
                                  min_cells_per_bin = 5,
                                  min_bins          = 8) {

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  # Align cells across mat, embed, and cluster_ids
  if (!is.null(rownames(embed))) {
    common <- intersect(colnames(mat), rownames(embed))
    if (!is.null(names(cluster_ids))) {
      common <- intersect(common, names(cluster_ids))
    }
    if (length(common) >= 2) {
      mat <- mat[, common, drop = FALSE]
      embed <- embed[common, , drop = FALSE]
      if (!is.null(names(cluster_ids))) {
        cluster_ids <- cluster_ids[common]
      }
    }
  }

  n <- ncol(mat)
  cluster_ids <- as.factor(cluster_ids)
  cl_levels <- levels(cluster_ids)
  n_types <- length(cl_levels)

  if (n_types < 2) return(rep(NA_real_, nrow(pair_dt)))

  # --- Assign micro-environment bins ---
  bins <- .assign_microenv_bins(embed, n_bins)

  # --- Pre-compute pseudo-bulk: mean expression per gene per type per bin ---
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))
  gene_idx <- stats::setNames(seq_along(unique_genes), unique_genes)

  bin_levels <- sort(unique(bins))
  n_actual_bins <- length(bin_levels)

  # 3D array: gene x type x bin (pseudo-bulk means)
  pb <- array(NA_real_, dim = c(length(unique_genes), n_types, n_actual_bins),
              dimnames = list(unique_genes, cl_levels, as.character(bin_levels)))
  # Count of cells per type per bin
  n_cells_tb <- matrix(0L, nrow = n_types, ncol = n_actual_bins,
                        dimnames = list(cl_levels, as.character(bin_levels)))

  for (bi in seq_along(bin_levels)) {
    b <- bin_levels[bi]
    cells_in_bin <- which(bins == b)
    for (ti in seq_along(cl_levels)) {
      ct <- cl_levels[ti]
      cells_ct <- cells_in_bin[cluster_ids[cells_in_bin] == ct]
      n_ct <- length(cells_ct)
      n_cells_tb[ti, bi] <- n_ct
      if (n_ct >= min_cells_per_bin) {
        pb[, ti, bi] <- rowMeans(mat[unique_genes, cells_ct, drop = FALSE])
      }
    }
  }

  # --- Compute per-pair scores ---
  n_pairs <- nrow(pair_dt)
  scores <- numeric(n_pairs)

  for (pi in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[pi]
    g2 <- pair_dt$gene2[pi]

    if (!(g1 %in% unique_genes) || !(g2 %in% unique_genes)) {
      scores[pi] <- NA_real_
      next
    }

    # Aggregate r(A->B) and r(B->A) across all type pairs
    r_ab_vals <- numeric(0)
    r_ba_vals <- numeric(0)
    w_ab <- numeric(0)
    w_ba <- numeric(0)

    for (ti in seq_along(cl_levels)) {
      for (tj in seq_along(cl_levels)) {
        if (ti == tj) next

        # Bins where both types have enough cells
        valid_bins <- which(n_cells_tb[ti, ] >= min_cells_per_bin &
                            n_cells_tb[tj, ] >= min_cells_per_bin)
        if (length(valid_bins) < min_bins) next

        # A->B: gene1 in type_i, gene2 in type_j
        x_ab <- pb[g1, ti, valid_bins]
        y_ab <- pb[g2, tj, valid_bins]

        # B->A: gene2 in type_i, gene1 in type_j
        x_ba <- pb[g2, ti, valid_bins]
        y_ba <- pb[g1, tj, valid_bins]

        # Skip if either vector is constant or all NA
        sd_x_ab <- stats::sd(x_ab, na.rm = TRUE)
        sd_y_ab <- stats::sd(y_ab, na.rm = TRUE)
        if (is.na(sd_x_ab) || is.na(sd_y_ab) ||
            sd_x_ab < 1e-10 || sd_y_ab < 1e-10) next

        r_ab <- stats::cor(x_ab, y_ab, use = "pairwise.complete.obs")
        if (!is.na(r_ab)) {
          r_ab_vals <- c(r_ab_vals, r_ab)
          w_ab <- c(w_ab, length(valid_bins))
        }

        sd_x_ba <- stats::sd(x_ba, na.rm = TRUE)
        sd_y_ba <- stats::sd(y_ba, na.rm = TRUE)
        if (is.na(sd_x_ba) || is.na(sd_y_ba) ||
            sd_x_ba < 1e-10 || sd_y_ba < 1e-10) next

        r_ba <- stats::cor(x_ba, y_ba, use = "pairwise.complete.obs")
        if (!is.na(r_ba)) {
          r_ba_vals <- c(r_ba_vals, r_ba)
          w_ba <- c(w_ba, length(valid_bins))
        }
      }
    }

    if (length(r_ab_vals) == 0 || length(r_ba_vals) == 0) {
      scores[pi] <- NA_real_
      next
    }

    # Weighted mean across type pairs
    agg_ab <- stats::weighted.mean(r_ab_vals, w_ab)
    agg_ba <- stats::weighted.mean(r_ba_vals, w_ba)

    abs_ab <- abs(agg_ab)
    abs_ba <- abs(agg_ba)

    scores[pi] <- if (abs_ab > 0 && abs_ba > 0) {
      sqrt(abs_ab * abs_ba)
    } else {
      0
    }
  }

  scores
}


# ===========================================================================
#  2. Cross-cell-type interaction score (single pair)
# ===========================================================================

#' Cross-cell-type interaction score for a single gene pair
#'
#' Measures whether gene A expressed in one cell type correlates with gene B
#' expressed in a different cell type sharing the same tissue micro-environment.
#'
#' Unlike the KNN-graph-based approach, this method does **not** require
#' cells of different types to be neighbours in PCA/UMAP space.  Instead,
#' it partitions cells into micro-environment bins (using k-means on the
#' embedding) and computes pseudo-bulk correlations across bins.
#'
#' @param x Numeric vector; expression of gene1 (length n_cells).
#' @param y Numeric vector; expression of gene2 (length n_cells).
#' @param cluster_ids Factor of cluster / cell-type assignments.
#' @param embed Numeric matrix of PCA embedding (n_cells x n_dims).
#' @param n_bins Integer; number of micro-environment bins.  Default 50.
#' @param min_cells_per_bin Integer; minimum cells of a type per bin.
#'   Default 5.
#' @param min_bins Integer; minimum valid bins per type pair.  Default 8.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{score}}{Geometric mean of |agg_r(A->B)| and |agg_r(B->A)|.}
#'     \item{\code{r_ab}}{Weighted mean r for gene1-in-source, gene2-in-neighbour.}
#'     \item{\code{r_ba}}{Weighted mean r for gene2-in-source, gene1-in-neighbour.}
#'     \item{\code{n_type_pairs}}{Number of cell-type pairs contributing.}
#'     \item{\code{per_celltype_pair}}{data.frame with per-type-pair breakdown.}
#'   }
#'
#' @keywords internal
.cross_celltype <- function(x, y, cluster_ids, embed,
                            n_bins            = 50,
                            min_cells_per_bin = 5,
                            min_bins          = 8) {

  # Align cells across x/y, embed, and cluster_ids
  if (!is.null(names(x)) && !is.null(rownames(embed))) {
    common <- intersect(names(x), rownames(embed))
    if (!is.null(names(cluster_ids))) {
      common <- intersect(common, names(cluster_ids))
    }
    if (length(common) >= 2) {
      x <- x[common]
      y <- y[common]
      embed <- embed[common, , drop = FALSE]
      if (!is.null(names(cluster_ids))) {
        cluster_ids <- cluster_ids[common]
      }
    } else {
      return(empty_result)
    }
  }

  n <- length(x)
  cluster_ids <- as.factor(cluster_ids)
  cl_levels <- levels(cluster_ids)
  n_types <- length(cl_levels)

  empty_result <- list(
    score = NA_real_, r_ab = NA_real_, r_ba = NA_real_,
    n_type_pairs = 0L, per_celltype_pair = data.frame()
  )

  if (n_types < 2) return(empty_result)

  # --- Assign micro-environment bins ---
  bins <- .assign_microenv_bins(embed, n_bins)

  bin_levels <- sort(unique(bins))
  n_actual_bins <- length(bin_levels)

  # --- Per-type per-bin pseudo-bulk ---
  # For a single pair we only need 2 genes; store as type x bin matrices
  pb_x <- matrix(NA_real_, nrow = n_types, ncol = n_actual_bins,
                  dimnames = list(cl_levels, as.character(bin_levels)))
  pb_y <- matrix(NA_real_, nrow = n_types, ncol = n_actual_bins,
                  dimnames = list(cl_levels, as.character(bin_levels)))
  n_cells_tb <- matrix(0L, nrow = n_types, ncol = n_actual_bins,
                        dimnames = list(cl_levels, as.character(bin_levels)))

  for (bi in seq_along(bin_levels)) {
    b <- bin_levels[bi]
    cells_in_bin <- which(bins == b)
    for (ti in seq_along(cl_levels)) {
      ct <- cl_levels[ti]
      cells_ct <- cells_in_bin[cluster_ids[cells_in_bin] == ct]
      n_ct <- length(cells_ct)
      n_cells_tb[ti, bi] <- n_ct
      if (n_ct >= min_cells_per_bin) {
        pb_x[ti, bi] <- mean(x[cells_ct])
        pb_y[ti, bi] <- mean(y[cells_ct])
      }
    }
  }

  # --- Per cell-type pair results ---
  ct_results <- data.frame(
    source_type     = character(0),
    neighbour_type  = character(0),
    n_bins_valid    = integer(0),
    r_g1_to_g2      = numeric(0),
    r_g2_to_g1      = numeric(0),
    stringsAsFactors = FALSE
  )

  r_ab_vals <- numeric(0)
  r_ba_vals <- numeric(0)
  w_ab <- numeric(0)
  w_ba <- numeric(0)

  for (ti in seq_along(cl_levels)) {
    for (tj in seq_along(cl_levels)) {
      if (ti == tj) next

      valid_bins <- which(n_cells_tb[ti, ] >= min_cells_per_bin &
                          n_cells_tb[tj, ] >= min_cells_per_bin)
      n_valid <- length(valid_bins)
      if (n_valid < min_bins) next

      # A->B: gene1 in type_i, gene2 in type_j
      xi <- pb_x[ti, valid_bins]
      yj <- pb_y[tj, valid_bins]

      # B->A: gene2 in type_i, gene1 in type_j
      yi <- pb_y[ti, valid_bins]
      xj <- pb_x[tj, valid_bins]

      sd_xi <- stats::sd(xi, na.rm = TRUE)
      sd_yj <- stats::sd(yj, na.rm = TRUE)
      r_ab_k <- if (!is.na(sd_xi) && !is.na(sd_yj) &&
                     sd_xi > 1e-10 && sd_yj > 1e-10) {
        stats::cor(xi, yj, use = "pairwise.complete.obs")
      } else {
        NA_real_
      }

      sd_yi <- stats::sd(yi, na.rm = TRUE)
      sd_xj <- stats::sd(xj, na.rm = TRUE)
      r_ba_k <- if (!is.na(sd_yi) && !is.na(sd_xj) &&
                     sd_yi > 1e-10 && sd_xj > 1e-10) {
        stats::cor(yi, xj, use = "pairwise.complete.obs")
      } else {
        NA_real_
      }

      ct_results <- rbind(ct_results, data.frame(
        source_type    = cl_levels[ti],
        neighbour_type = cl_levels[tj],
        n_bins_valid   = n_valid,
        r_g1_to_g2     = r_ab_k,
        r_g2_to_g1     = r_ba_k,
        stringsAsFactors = FALSE
      ))

      if (!is.na(r_ab_k)) {
        r_ab_vals <- c(r_ab_vals, r_ab_k)
        w_ab <- c(w_ab, n_valid)
      }
      if (!is.na(r_ba_k)) {
        r_ba_vals <- c(r_ba_vals, r_ba_k)
        w_ba <- c(w_ba, n_valid)
      }
    }
  }

  if (length(r_ab_vals) == 0 || length(r_ba_vals) == 0) return(empty_result)

  agg_ab <- stats::weighted.mean(r_ab_vals, w_ab)
  agg_ba <- stats::weighted.mean(r_ba_vals, w_ba)

  score <- if (abs(agg_ab) > 0 && abs(agg_ba) > 0) {
    sqrt(abs(agg_ab) * abs(agg_ba))
  } else {
    0
  }

  list(
    score             = score,
    r_ab              = agg_ab,
    r_ba              = agg_ba,
    n_type_pairs      = nrow(ct_results),
    per_celltype_pair = ct_results
  )
}
