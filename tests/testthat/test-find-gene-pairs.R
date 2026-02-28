test_that("FindGenePairs returns correct structure", {
  result <- FindGenePairs(
    scpairs_testdata,
    gene                = "GENE3",
    candidates          = paste0("GENE", 1:10),
    cor_method          = "pearson",
    n_mi_bins           = 0,
    min_cells_expressed = 5,
    n_perm              = 0,
    verbose             = FALSE
  )

  expect_s3_class(result, "scPairs_gene_result")
  expect_equal(result$query_gene, "GENE3")
  expect_true("gene2"        %in% colnames(result$pairs))
  expect_true("synergy_score" %in% colnames(result$pairs))
  expect_true(all(result$pairs$gene1 == "GENE3"))
})


test_that("FindGenePairs identifies GENE4 as top partner of GENE3", {
  result <- FindGenePairs(
    scpairs_testdata,
    gene                = "GENE3",
    candidates          = paste0("GENE", c(1, 2, 4, 5, 6, 7)),
    cor_method          = c("pearson", "spearman"),
    n_mi_bins           = 0,
    min_cells_expressed = 5,
    n_perm              = 0,
    verbose             = FALSE
  )

  expect_equal(result$pairs$gene2[1], "GENE4")
})


test_that("FindGenePairs errors for a missing gene", {
  expect_error(
    FindGenePairs(scpairs_testdata, gene = "NONEXISTENT", verbose = FALSE),
    "not found|None of the requested features"
  )
})
