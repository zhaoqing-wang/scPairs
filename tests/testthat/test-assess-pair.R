test_that("AssessGenePair returns correct structure", {
  result <- AssessGenePair(
    scpairs_testdata,
    gene1   = "GENE3",
    gene2   = "GENE4",
    n_perm  = 0,
    verbose = FALSE
  )

  expect_s3_class(result, "scPairs_pair_result")
  expect_equal(result$gene1, "GENE3")
  expect_equal(result$gene2, "GENE4")
  expect_true(is.numeric(result$synergy_score))
  expect_true(result$synergy_score >= 0 && result$synergy_score <= 1)
  expect_true(is.data.frame(result$per_cluster))
  expect_true("cor_pearson"  %in% names(result$metrics))
  expect_true("jaccard_index" %in% names(result$metrics))
})


test_that("AssessGenePair gives higher score for known co-expressed pair", {
  result_good <- AssessGenePair(
    scpairs_testdata,
    gene1   = "GENE3",
    gene2   = "GENE4",
    n_perm  = 0,
    verbose = FALSE
  )

  result_bad <- AssessGenePair(
    scpairs_testdata,
    gene1   = "GENE3",
    gene2   = "GENE10",
    n_perm  = 0,
    verbose = FALSE
  )

  expect_true(result_good$synergy_score > result_bad$synergy_score)
})


test_that("AssessGenePair permutation p-value is in [0, 1]", {
  result <- AssessGenePair(
    scpairs_testdata,
    gene1   = "GENE3",
    gene2   = "GENE4",
    n_perm  = 49,
    verbose = FALSE
  )

  expect_true(!is.na(result$p_value))
  expect_true(result$p_value > 0 && result$p_value <= 1)
  expect_true(result$confidence %in% c("High", "Medium", "Low", "NS"))
})


test_that("AssessGenePair errors when gene1 == gene2", {
  expect_error(
    AssessGenePair(scpairs_testdata,
                   gene1 = "GENE1", gene2 = "GENE1", verbose = FALSE),
    "must be different"
  )
})
