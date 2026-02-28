test_that("FindAllPairs returns correct structure", {
  result <- FindAllPairs(
    scpairs_testdata,
    n_top_genes         = 20,
    cor_method          = c("pearson", "spearman"),
    n_mi_bins           = 3,
    min_cells_expressed = 5,
    n_perm              = 0,
    verbose             = FALSE
  )

  expect_s3_class(result, "scPairs_result")
  expect_true(is.data.frame(result$pairs) || data.table::is.data.table(result$pairs))
  expect_true("gene1"        %in% colnames(result$pairs))
  expect_true("gene2"        %in% colnames(result$pairs))
  expect_true("synergy_score" %in% colnames(result$pairs))
  expect_true("rank"         %in% colnames(result$pairs))
  expect_true("confidence"   %in% colnames(result$pairs))
  expect_true(result$n_genes > 0)
  expect_equal(result$n_cells, ncol(scpairs_testdata))
})


test_that("FindAllPairs detects known co-expression (GENE3 / GENE4)", {
  result <- FindAllPairs(
    scpairs_testdata,
    features            = paste0("GENE", 1:10),
    cor_method          = "pearson",
    n_mi_bins           = 0,
    min_cells_expressed = 5,
    n_perm              = 0,
    verbose             = FALSE
  )

  pairs  <- result$pairs
  top10  <- pairs[pairs$rank <= 10, ]
  gene34 <- any(
    (top10$gene1 == "GENE3" & top10$gene2 == "GENE4") |
    (top10$gene1 == "GENE4" & top10$gene2 == "GENE3")
  )
  expect_true(gene34)
})


test_that("FindAllPairs top_n limits output", {
  result <- FindAllPairs(
    scpairs_testdata,
    n_top_genes         = 20,
    cor_method          = "pearson",
    n_mi_bins           = 0,
    min_cells_expressed = 3,
    top_n               = 5,
    n_perm              = 0,
    verbose             = FALSE
  )

  expect_true(nrow(result$pairs) <= 5)
})


test_that("FindAllPairs errors when fewer than 2 features are provided", {
  expect_error(
    FindAllPairs(scpairs_testdata, features = "GENE1", verbose = FALSE),
    "At least 2"
  )
})
