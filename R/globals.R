#' @importFrom data.table := data.table frank setorder setnames
#' @importFrom stats cor median sd var quantile setNames
#' @importFrom Matrix sparseMatrix t
NULL

utils::globalVariables(
  c("gene1", "gene2", "cor_pearson", "cor_spearman", "cor_biweight",
    "mean_cor", "mi_score", "ratio_consistency", "spatial_lee_L",
    "spatial_clq", "synergy_score", "p_value", "p_adj", "rank",
    "confidence", "category", "weight", "from", "to", "name",
    "x_coord", "y_coord", "expr1", "expr2", "coexpr_product",
    "cell_id", "cluster", "mean_expr", "pct_expressed",
    ".", ".N", ".SD", "V1", "value", "variable",
    "node_size", "edge_width", "edge_color",
    "pair_label", "score_type", "score_value",
    "dim1", "dim2", "x", "y", "gene", "spatial_lee_p")
)
