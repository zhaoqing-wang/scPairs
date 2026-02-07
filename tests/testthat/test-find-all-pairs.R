test_that("FindAllPairs returns correct structure", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 100, n_genes = 20)

  result <- FindAllPairs(
    sce,
    n_top_genes = 20,
    cor_method = c("pearson", "spearman"),
    n_mi_bins = 3,
    min_cells_expressed = 5,
    n_perm = 0,
    verbose = FALSE
  )

  expect_s3_class(result, "scPairs_result")
  expect_true(is.data.frame(result$pairs) || data.table::is.data.table(result$pairs))
  expect_true("gene1" %in% colnames(result$pairs))
  expect_true("gene2" %in% colnames(result$pairs))
  expect_true("synergy_score" %in% colnames(result$pairs))
  expect_true("rank" %in% colnames(result$pairs))
  expect_true("confidence" %in% colnames(result$pairs))
  expect_true(result$n_genes > 0)
  expect_true(result$n_cells == 100)
})


test_that("FindAllPairs detects known co-expression", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 200, n_genes = 20)

  result <- FindAllPairs(
    sce,
    features = paste0("GENE", 1:10),
    cor_method = c("pearson"),
    n_mi_bins = 0,
    min_cells_expressed = 5,
    n_perm = 0,
    verbose = FALSE
  )

  # GENE3/GENE4 should be among top pairs
  pairs <- result$pairs
  top10 <- pairs[pairs$rank <= 10, ]
  gene34_present <- any(
    (top10$gene1 == "GENE3" & top10$gene2 == "GENE4") |
    (top10$gene1 == "GENE4" & top10$gene2 == "GENE3")
  )
  expect_true(gene34_present)
})


test_that("FindAllPairs top_n works", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 100, n_genes = 15)

  result <- FindAllPairs(
    sce,
    n_top_genes = 15,
    cor_method = "pearson",
    n_mi_bins = 0,
    min_cells_expressed = 3,
    top_n = 5,
    n_perm = 0,
    verbose = FALSE
  )

  expect_true(nrow(result$pairs) <= 5)
})


test_that("FindAllPairs handles too few features gracefully", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 50, n_genes = 10)

  expect_error(
    FindAllPairs(sce, features = "GENE1", verbose = FALSE),
    "At least 2"
  )
})
