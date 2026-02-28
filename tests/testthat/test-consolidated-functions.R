# Tests for consolidated (unified interface) functions

test_that(".compute_bicor dispatches correctly in single-pair and matrix modes", {
  x <- rnorm(100)
  y <- x + rnorm(100, sd = 0.1)

  bc_single <- scPairs:::.compute_bicor(x, y)
  expect_true(is.numeric(bc_single) && length(bc_single) == 1)
  expect_true(bc_single > 0.8)
  expect_equal(bc_single, scPairs:::.bicor(x, y))

  mat <- rbind(x, y, rnorm(100))
  rownames(mat) <- c("A", "B", "C")
  bc_mat <- scPairs:::.compute_bicor(mat)
  expect_true(is.matrix(bc_mat))
  expect_equal(dim(bc_mat), c(3L, 3L))
})


test_that(".compute_mi dispatches correctly in single-pair and batch modes", {
  n <- 200
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  mi_single <- scPairs:::.compute_mi(x, y, n_bins = 5)
  expect_true(mi_single >= 0)
  expect_equal(mi_single, scPairs:::.mutual_info(x, y, n_bins = 5))

  mat      <- rbind(x, y, rnorm(n))
  pair_idx <- matrix(c(1, 2, 1, 3), nrow = 2)
  mi_batch <- scPairs:::.compute_mi(mat, pair_idx = pair_idx, n_bins = 5)
  expect_length(mi_batch, 2)
  expect_true(all(mi_batch >= 0))
})


test_that(".compute_ratio_consistency dispatches correctly", {
  n  <- 100
  x  <- rnorm(n, 2, 0.1)
  y  <- x * 1.5 + rnorm(n, 0, 0.05)
  cl <- factor(rep(1:3, length.out = n))

  rc_single <- scPairs:::.compute_ratio_consistency(x, y, cluster_ids = cl)
  expect_true(rc_single >= 0 && rc_single <= 1)
  expect_equal(rc_single, scPairs:::.ratio_consistency(x, y, cl))

  mat      <- rbind(x, y, abs(rnorm(n)))
  pair_idx <- matrix(c(1, 2, 1, 3), nrow = 2)
  rc_batch <- scPairs:::.compute_ratio_consistency(mat,
                                                   pair_idx   = pair_idx,
                                                   cluster_ids = cl)
  expect_length(rc_batch, 2)
  expect_true(all(rc_batch >= 0 & rc_batch <= 1))
})


test_that(".compute_cross_celltype dispatches correctly", {
  mat         <- as.matrix(Seurat::GetAssayData(scpairs_testdata, layer = "data"))
  embed       <- Seurat::Embeddings(scpairs_testdata, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(scpairs_testdata))

  # Single-pair mode
  result <- scPairs:::.compute_cross_celltype(
    mat["GENE3", ], mat["GENE4", ],
    cluster_ids = cluster_ids,
    embed       = embed,
    n_bins = 15, min_cells_per_bin = 2, min_bins = 3
  )
  expect_type(result, "list")
  expect_true("score" %in% names(result))

  # Batch mode
  pair_dt <- data.table::data.table(gene1 = c("GENE1", "GENE3"),
                                    gene2 = c("GENE2", "GENE4"))
  scores <- scPairs:::.compute_cross_celltype(
    mat, pair_dt = pair_dt,
    cluster_ids = cluster_ids,
    embed       = embed,
    n_bins = 15, min_cells_per_bin = 2, min_bins = 3
  )
  expect_length(scores, 2)
  expect_true(is.numeric(scores))
})


test_that("schema constants are self-consistent", {
  all_metrics  <- c(scPairs:::ABS_METRICS, scPairs:::RAW_METRICS)
  weight_names <- names(scPairs:::DEFAULT_WEIGHTS)
  expect_true(all(weight_names %in% all_metrics))

  for (nm in names(scPairs:::RESULT_COLUMNS)) {
    expect_true(is.character(scPairs:::RESULT_COLUMNS[[nm]]))
  }
})
