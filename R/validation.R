#' Input Validation Helpers for scPairs
#'
#' @description
#' Internal functions for validating user inputs to main scPairs functions.
#' These ensure parameters are within reasonable ranges and provide helpful
#' error messages.
#'
#' @name validation
#' @keywords internal
NULL


#' Validate minimum cells expressed parameter
#' @keywords internal
.validate_min_cells_expressed <- function(min_cells_expressed, n_cells) {
  if (!is.numeric(min_cells_expressed) || length(min_cells_expressed) != 1) {
    stop("min_cells_expressed must be a single numeric value", call. = FALSE)
  }
  
  if (min_cells_expressed < 0) {
    stop("min_cells_expressed must be non-negative", call. = FALSE)
  }
  
  max_reasonable <- n_cells / 10
  if (min_cells_expressed > max_reasonable) {
    warning(sprintf(
      paste0("min_cells_expressed (%d) is very high relative to n_cells (%d). ",
             "This may filter out most gene pairs. ",
             "Consider using a value < %d."),
      min_cells_expressed, n_cells, floor(max_reasonable)
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Validate neighbourhood parameters
#' @keywords internal
.validate_neighbourhood_params <- function(neighbourhood_k, n_cells) {
  if (!is.numeric(neighbourhood_k) || length(neighbourhood_k) != 1) {
    stop("neighbourhood_k must be a single numeric value", call. = FALSE)
  }
  
  if (neighbourhood_k < 2) {
    stop("neighbourhood_k must be >= 2", call. = FALSE)
  }
  
  if (neighbourhood_k > n_cells / 2) {
    warning(sprintf(
      paste0("neighbourhood_k (%d) is more than half of n_cells (%d). ",
             "This may produce overly smoothed results."),
      neighbourhood_k, n_cells
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Validate binning parameters for cross-cell-type analysis
#' @keywords internal
.validate_binning_params <- function(n_bins, min_cells_per_bin, min_bins, 
                                     n_cells, n_types) {
  if (!is.numeric(n_bins) || length(n_bins) != 1 || n_bins < 1) {
    stop("n_bins must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(min_cells_per_bin) || length(min_cells_per_bin) != 1 || 
      min_cells_per_bin < 1) {
    stop("min_cells_per_bin must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(min_bins) || length(min_bins) != 1 || min_bins < 1) {
    stop("min_bins must be a positive integer", call. = FALSE)
  }
  
  # Check if binning is feasible
  min_cells_needed <- min_cells_per_bin * min_bins * n_types
  if (min_cells_needed > n_cells) {
    warning(sprintf(
      paste0("Binning parameters may be too strict: require %d cells but only %d available. ",
             "Consider reducing n_bins, min_cells_per_bin, or min_bins."),
      min_cells_needed, n_cells
    ), call. = FALSE)
  }
  
  max_reasonable_bins <- n_cells / (min_cells_per_bin * 5)
  if (n_bins > max_reasonable_bins) {
    warning(sprintf(
      paste0("n_bins (%d) may be too large for dataset size. ",
             "Consider using <= %d bins."),
      n_bins, floor(max_reasonable_bins)
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Validate percentage parameter
#' @keywords internal
.validate_percentage <- function(value, param_name) {
  if (!is.numeric(value) || length(value) != 1) {
    stop(sprintf("%s must be a single numeric value", param_name), call. = FALSE)
  }
  
  if (value < 0 || value > 1) {
    stop(sprintf("%s must be between 0 and 1", param_name), call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Validate feature names exist in Seurat object
#' @keywords internal
.validate_features <- function(features, object, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(object)
  
  available_features <- rownames(object[[assay]])
  missing <- setdiff(features, available_features)
  
  if (length(missing) > 0) {
    if (length(missing) == length(features)) {
      stop(sprintf(
        "None of the requested features found in assay '%s'",
        assay
      ), call. = FALSE)
    } else {
      warning(sprintf(
        "%d feature(s) not found in assay '%s': %s",
        length(missing), assay,
        paste(head(missing, 5), collapse = ", ")
      ), call. = FALSE)
    }
  }
  
  invisible(intersect(features, available_features))
}


#' Validate correlation method
#' @keywords internal
.validate_cor_method <- function(cor_method) {
  valid_methods <- c("pearson", "spearman", "biweight")
  
  if (length(cor_method) == 0) {
    stop("cor_method must be specified", call. = FALSE)
  }
  
  invalid <- setdiff(cor_method, valid_methods)
  if (length(invalid) > 0) {
    stop(sprintf(
      "Invalid correlation method(s): %s. Valid options: %s",
      paste(invalid, collapse = ", "),
      paste(valid_methods, collapse = ", ")
    ), call. = FALSE)
  }
  
  invisible(cor_method)
}


#' Validate permutation number
#' @keywords internal
.validate_n_perm <- function(n_perm) {
  if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm < 0) {
    stop("n_perm must be a non-negative integer", call. = FALSE)
  }
  
  if (n_perm > 0 && n_perm < 100) {
    warning(
      "n_perm < 100 may not provide reliable p-values. ",
      "Consider using n_perm >= 999 for publication-quality results.",
      call. = FALSE
    )
  }
  
  if (n_perm > 10000) {
    message(sprintf(
      paste0("n_perm = %d will be computationally intensive. ",
             "This may take considerable time."), n_perm
    ))
  }
  
  invisible(TRUE)
}
