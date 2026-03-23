#' Filter repertoire features based on imputation error
#'
#' Removes features that cannot be reliably imputed given their observed
#' level of missingness, based on RMSE benchmarking.
#'
#' @param feature_matrix Numeric matrix or data.frame (samples × features).
#' @param benchmark_results Output from benchmark_imputation().
#' @param method Character. Imputation method used for filtering (default "missForest").
#' @param error_threshold Numeric. Error metric cutoff above which features are removed.
#' @param error_metric Character. Metric for filtering: "nrmse_mad" (default),
#'   "nrmse_sd", "nrmse_range", or "rmse".
#'
#' @return A list with:
#' \item{filtered_matrix}{Feature matrix after removing poorly imputed features.}
#' \item{removed_features}{Character vector of removed feature names.}
#' \item{summary}{Data frame summarising filtering decisions.}
#'
#' @details
#' A feature is removed if its error metric exceeds the threshold at any NA-seeding
#' frequency less than or equal to its observed missingness.
#'
#' @export
filter_features_by_imputation_error_variant <- function(
    feature_matrix,
    benchmark_results,
    method = "missForest",
    error_threshold = 0.5,
    error_metric = c("nrmse_mad", "nrmse_sd", "nrmse_range", "rmse")
) {
  stopifnot(is.matrix(feature_matrix) || is.data.frame(feature_matrix))
  stopifnot(is.data.frame(benchmark_results))

  error_metric <- match.arg(error_metric)

  mat <- as.matrix(feature_matrix)

  # Observed missingness per feature
  observed_na <- colMeans(is.na(mat))

  # Restrict to selected imputation method
  bench_method <- benchmark_results[
    benchmark_results$method == method, , drop = FALSE
  ]

  if (nrow(bench_method) == 0) {
    stop("No benchmark results found for method: ", method)
  }
  
  # Determine which error metric to use for filtering
  if (!error_metric %in% colnames(benchmark_results)) {
    stop("Metric column not found in benchmark_results: ", error_metric)
  }

  # Identify problematic features
  remove_features <- c()
  summary_list <- list()

  for (feat in unique(bench_method$feature)) {

    obs_na <- observed_na[feat]

    feat_rows <- bench_method[
      bench_method$feature == feat &
        bench_method$na_fraction <= obs_na,
      , drop = FALSE
    ]

    remove_flag <- any(feat_rows[[error_metric]] > error_threshold, na.rm = TRUE)

    if (remove_flag) {
      remove_features <- c(remove_features, feat)
    }

    summary_list[[feat]] <- data.frame(
      feature = feat,
      observed_na = obs_na,
      # Keeping "max_rmse" as other functions might expect this column name
      max_rmse = ifelse(nrow(feat_rows) > 0,
                        max(feat_rows[[error_metric]], na.rm = TRUE),
                        NA_real_),
      removed = remove_flag
    )
  }

  summary_df <- do.call(rbind, summary_list)

  filtered_mat <- mat[, !colnames(mat) %in% remove_features, drop = FALSE]

  list(
    filtered_matrix = filtered_mat,
    removed_features = unique(remove_features),
    summary = summary_df
  )
}