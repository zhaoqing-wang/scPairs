#' Centralized Schema Definitions for scPairs
#'
#' @description
#' Single source of truth for column names, default parameters, metric weights,
#' and classification thresholds used throughout the package.
#'
#' @name schema
#' @keywords internal
NULL


# ===========================================================================
#  Result column names
# ===========================================================================

#' Standard column names for result data.tables
#' @keywords internal
RESULT_COLUMNS <- list(
  # Core pair identifiers
  pair_ids = c("gene1", "gene2"),

  # Co-expression metrics
  coexpr = c("cor_pearson", "cor_spearman", "cor_biweight",
             "mi_score", "ratio_consistency"),

  # Neighbourhood metrics
  neighbourhood = c("smoothed_cor", "neighbourhood_score",
                     "cluster_cor", "cross_celltype_score"),

  # Spatial metrics
  spatial = c("spatial_lee_L", "spatial_lee_p", "spatial_clq"),

  # Score integration output
  integration = c("synergy_score", "p_value", "p_adj", "rank", "confidence")
)


# ===========================================================================
#  Default parameters
# ===========================================================================

#' Default analysis parameters
#' @keywords internal
DEFAULT_PARAMS <- list(
  # Feature selection
  n_top_genes         = 2000,


  # Correlation
  cor_method          = c("pearson", "spearman", "biweight"),
  n_mi_bins           = 5,

  # Co-expression filtering
  min_cells_expressed = 10,

  # Neighbourhood
  neighbourhood_k     = 20,
  smooth_alpha        = 0.3,

  # Cross-cell-type
  n_bins              = 50,
  min_cells_per_bin   = 5,
  min_bins            = 8,
  min_pct_expressed   = 0.01,

  # Spatial
  spatial_k           = 6,

  # Permutation
  n_perm              = 0,
  n_perm_assess       = 999
)


# ===========================================================================
#  Metric weights
# ===========================================================================

#' Default metric weights for score integration
#' @keywords internal
DEFAULT_WEIGHTS <- c(
  cor_pearson          = 1.0,
  cor_spearman         = 1.0,

  cor_biweight         = 1.5,
  mi_score             = 1.0,
  ratio_consistency    = 1.2,
  smoothed_cor         = 1.5,
  neighbourhood_score  = 1.5,
  cluster_cor          = 1.2,
  cross_celltype_score = 1.5,
  spatial_lee_L        = 1.5,
  spatial_clq          = 1.2
)


# ===========================================================================
#  Metric classification
# ===========================================================================

#' Metrics that should be absolute-valued before rank normalisation
#' @keywords internal
ABS_METRICS <- c("cor_pearson", "cor_spearman", "cor_biweight",
                 "smoothed_cor", "cluster_cor", "spatial_lee_L")

#' Metrics used as raw values (not absolute-valued)
#' @keywords internal
RAW_METRICS <- c("mi_score", "ratio_consistency",
                 "neighbourhood_score", "cross_celltype_score",
                 "spatial_clq")


# ===========================================================================
#  Confidence thresholds
# ===========================================================================

#' Confidence classification thresholds (based on adjusted p-values)
#' @keywords internal
CONFIDENCE_THRESHOLDS <- list(
  High   = 0.01,
  Medium = 0.05,
  Low    = 0.10
)

#' Score-based confidence thresholds (when no p-values available)
#' @keywords internal
SCORE_CONFIDENCE_QUANTILES <- list(
  High   = 0.95,
  Medium = 0.80,
  Low    = 0.50
)
