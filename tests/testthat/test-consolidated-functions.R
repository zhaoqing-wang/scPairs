# Tests for consolidated (unified interface) functions

test_that(".compute_bicor dispatches correctly", {
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100, sd = 0.1)

  # Single-pair mode
  bc_single <- scPairs:::.compute_bicor(x, y)
  expect_true(is.numeric(bc_single))
  expect_length(bc_single, 1)
  expect_true(bc_single > 0.8)

  # Should match .bicor directly
  bc_direct <- scPairs:::.bicor(x, y)
  expect_equal(bc_single, bc_direct)

  # Matrix mode
  mat <- rbind(x, y, rnorm(100))
  rownames(mat) <- c("A", "B", "C")
  bc_mat <- scPairs:::.compute_bicor(mat)
  expect_true(is.matrix(bc_mat))
  expect_equal(nrow(bc_mat), 3)
  expect_equal(ncol(bc_mat), 3)
})


test_that(".compute_mi dispatches correctly", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  # Single-pair mode
  mi_single <- scPairs:::.compute_mi(x, y, n_bins = 5)
  expect_true(mi_single >= 0)

  # Should match .mutual_info directly
  mi_direct <- scPairs:::.mutual_info(x, y, n_bins = 5)
  expect_equal(mi_single, mi_direct)

  # Batch mode
  mat <- rbind(x, y, rnorm(n))
  pair_idx <- matrix(c(1, 2, 1, 3), nrow = 2)
  mi_batch <- scPairs:::.compute_mi(mat, pair_idx = pair_idx, n_bins = 5)
  expect_length(mi_batch, 2)
  expect_true(all(mi_batch >= 0))
})


test_that(".compute_ratio_consistency dispatches correctly", {
  set.seed(42)
  n <- 100
  x <- rnorm(n, 2, 0.1)
  y <- x * 1.5 + rnorm(n, 0, 0.05)
  cl <- factor(rep(1:3, length.out = n))

  # Single-pair mode
  rc_single <- scPairs:::.compute_ratio_consistency(x, y, cluster_ids = cl)
  expect_true(rc_single >= 0 && rc_single <= 1)

  # Should match .ratio_consistency directly
  rc_direct <- scPairs:::.ratio_consistency(x, y, cl)
  expect_equal(rc_single, rc_direct)

  # Batch mode
  mat <- rbind(x, y, rnorm(n))
  pair_idx <- matrix(c(1, 2, 1, 3), nrow = 2)
  rc_batch <- scPairs:::.compute_ratio_consistency(mat,
                                                    pair_idx = pair_idx,
                                                    cluster_ids = cl)
  expect_length(rc_batch, 2)
  expect_true(all(rc_batch >= 0 & rc_batch <= 1))
})


test_that(".compute_cross_celltype dispatches correctly", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 200, n_genes = 20, n_clusters = 3)
  mat <- as.matrix(Seurat::GetAssayData(sce, layer = "data"))
  embed <- Seurat::Embeddings(sce, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(sce))

  # Single-pair mode
  result <- scPairs:::.compute_cross_celltype(
    mat["GENE3", ], mat["GENE4", ],
    cluster_ids = cluster_ids,
    embed = embed,
    n_bins = 15, min_cells_per_bin = 2, min_bins = 3
  )
  expect_type(result, "list")
  expect_true("score" %in% names(result))

  # Batch mode
  pair_dt <- data.table::data.table(
    gene1 = c("GENE1", "GENE3"),
    gene2 = c("GENE2", "GENE4")
  )
  scores <- scPairs:::.compute_cross_celltype(
    mat, pair_dt = pair_dt,
    cluster_ids = cluster_ids,
    embed = embed,
    n_bins = 15, min_cells_per_bin = 2, min_bins = 3
  )
  expect_length(scores, 2)
  expect_true(is.numeric(scores))
})


test_that("schema constants are consistent", {
  # All weight names should be in either ABS_METRICS or RAW_METRICS
  all_metrics <- c(scPairs:::ABS_METRICS, scPairs:::RAW_METRICS)
  weight_names <- names(scPairs:::DEFAULT_WEIGHTS)
  expect_true(all(weight_names %in% all_metrics))

  # All result columns should be character vectors
  for (nm in names(scPairs:::RESULT_COLUMNS)) {
    expect_true(is.character(scPairs:::RESULT_COLUMNS[[nm]]))
  }
})
