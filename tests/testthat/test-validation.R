# Tests for R/validation.R helper functions

test_that(".validate_min_cells_expressed rejects invalid input", {
  expect_error(
    scPairs:::.validate_min_cells_expressed("abc", 100),
    "single numeric"
  )
  expect_error(
    scPairs:::.validate_min_cells_expressed(-5, 100),
    "non-negative"
  )
})

test_that(".validate_min_cells_expressed warns for unusually large values", {
  expect_warning(
    scPairs:::.validate_min_cells_expressed(50, 100),
    "very high"
  )
  expect_silent(
    scPairs:::.validate_min_cells_expressed(5, 100)
  )
})


test_that(".validate_neighbourhood_params rejects bad k values", {
  expect_error(
    scPairs:::.validate_neighbourhood_params("not_a_number", 100),
    "single numeric"
  )
  expect_error(
    scPairs:::.validate_neighbourhood_params(1, 100),
    ">= 2"
  )
})

test_that(".validate_neighbourhood_params warns for large k", {
  expect_warning(
    scPairs:::.validate_neighbourhood_params(60, 100),
    "more than half"
  )
  expect_silent(
    scPairs:::.validate_neighbourhood_params(20, 100)
  )
})


test_that(".validate_percentage rejects values outside [0, 1]", {
  expect_error(scPairs:::.validate_percentage(-0.1, "alpha"), "between 0 and 1")
  expect_error(scPairs:::.validate_percentage( 1.5, "alpha"), "between 0 and 1")
  expect_silent(scPairs:::.validate_percentage(0.5,  "alpha"))
})


test_that(".validate_cor_method rejects invalid correlation methods", {
  expect_error(
    scPairs:::.validate_cor_method(c("pearson", "invalid")),
    "Invalid"
  )
  expect_error(
    scPairs:::.validate_cor_method(character(0)),
    "specified"
  )
  expect_silent(scPairs:::.validate_cor_method(c("pearson", "spearman")))
})


test_that(".validate_n_perm rejects negative values and warns for small counts", {
  expect_error(
    scPairs:::.validate_n_perm(-1),
    "non-negative"
  )
  expect_warning(
    scPairs:::.validate_n_perm(50),
    "< 100"
  )
  expect_silent(scPairs:::.validate_n_perm(999))
})


test_that(".validate_features handles missing genes correctly", {
  # All missing -> error
  expect_error(
    scPairs:::.validate_features(c("NONEXISTENT1", "NONEXISTENT2"), scpairs_testdata),
    "None of the requested features"
  )

  # Some missing -> warning, returns the intersection
  expect_warning(
    result <- scPairs:::.validate_features(c("GENE1", "NONEXISTENT"), scpairs_testdata),
    "not found"
  )
  expect_equal(result, "GENE1")

  # All present -> silent pass
  expect_silent(
    result2 <- scPairs:::.validate_features(c("GENE1", "GENE2"), scpairs_testdata)
  )
  expect_equal(sort(result2), c("GENE1", "GENE2"))
})


test_that(".validate_binning_params warns for infeasible parameters", {
  expect_warning(
    scPairs:::.validate_binning_params(
      n_bins = 100, min_cells_per_bin = 5,
      min_bins = 8, n_cells = 50, n_types = 3
    ),
    "too strict"
  )
  expect_silent(
    scPairs:::.validate_binning_params(
      n_bins = 10, min_cells_per_bin = 3,
      min_bins = 5, n_cells = 500, n_types = 3
    )
  )
})
