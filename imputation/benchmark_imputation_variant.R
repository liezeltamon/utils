#' Benchmark imputation methods for repertoire feature matrices
#'
#' Benchmarks multiple imputation strategies by seeding missing values
#' into complete data and quantifying per-feature RMSE.
#'
#' @param feature_matrix Numeric matrix or data.frame (samples × features).
#' @param na_frequencies Numeric vector of NA fractions to seed (e.g. c(0.1, 0.2)).
#' @param seed_strategy Character: "feature" (column-wise) or "sample" (row-wise).
#' @param methods Character vector of imputation methods to benchmark.
#' @param random_seed Integer for reproducibility.
#'
#' @return A data.frame with columns:
#'   method, feature, na_fraction, rmse, nrmse_sd, nrmse_mad, nrmse_range
#'
#' @details
#' RMSE is computed by comparing imputed values against original values
#' at artificially masked entries only.
#'
#' @export
benchmark_imputation <- function(
    feature_matrix,
    na_frequencies = c(0.1, 0.2, 0.3),
    seed_strategy = c("feature", "sample"),
    methods = c("mean", "median", "random", "knn", "missMDA_reg", "missMDA_EM", "missForest"),
    random_seed = 1
) {
  seed_strategy <- match.arg(seed_strategy)

  stopifnot(is.matrix(feature_matrix) || is.data.frame(feature_matrix))
  mat <- as.matrix(feature_matrix)
  mat <- apply(mat, 2, as.numeric)

  # Use only complete cases as ground truth
  complete_mat <- mat[complete.cases(mat), , drop = FALSE]

  if (nrow(complete_mat) < 5) {
    stop("Not enough complete samples to benchmark imputation.")
  }

  set.seed(random_seed)

  results <- list()

  for (na_frac in na_frequencies) {

    # Seed missingness
    if (seed_strategy == "feature") {
      mat_na <- apply(
        complete_mat,
        2,
        function(x) {
          missForest::prodNA(data.frame(x), noNA = na_frac)[, 1]
        }
      )
    } else {
      mat_na <- missForest::prodNA(complete_mat, noNA = na_frac)
    }

    mat_na <- as.matrix(mat_na)
    rownames(mat_na) <- rownames(complete_mat)
    colnames(mat_na) <- colnames(complete_mat)

    # Identify masked positions
    na_mask <- is.na(mat_na)

    for (method in methods) {

      imputed <- NULL

      # ---- Imputation methods ----
      if (method == "mean") {
        imputed <- apply(mat_na, 2, function(x) Hmisc::impute(x, fun = mean))
      }

      if (method == "median") {
        imputed <- apply(mat_na, 2, function(x) Hmisc::impute(x, fun = median))
      }

      if (method == "random") {
        imputed <- apply(mat_na, 2, function(x) Hmisc::impute(x, fun = "random"))
      }

      if (method == "knn") {
        tmp <- try(VIM::kNN(mat_na, k = 5), silent = TRUE)
        if (!inherits(tmp, "try-error")) {
          imputed <- as.matrix(tmp[, colnames(mat_na)])
        }
      }

      if (method == "missMDA_reg") {
        tmp <- try(missMDA::estim_ncpPCA(mat_na), silent = TRUE)
        if (!inherits(tmp, "try-error")) {
          imp <- missMDA::imputePCA(mat_na, ncp = tmp$ncp)$completeObs
          imputed <- as.matrix(imp)
        }
      }

      if (method == "missMDA_EM") {
        tmp <- try(missMDA::estim_ncpPCA(mat_na), silent = TRUE)
        if (!inherits(tmp, "try-error")) {
          imp <- missMDA::imputePCA(mat_na, ncp = tmp$ncp, method = "EM")$completeObs
          imputed <- as.matrix(imp)
        }
      }

      if (method == "missForest") {
        tmp <- try(missForest::missForest(mat_na, xtrue = complete_mat), silent = TRUE)
        if (!inherits(tmp, "try-error")) {
          imputed <- as.matrix(tmp$ximp)
        }
      }

      if (is.null(imputed)) next

      # ---- RMSE per feature ----
      for (j in seq_len(ncol(mat_na))) {
        idx <- na_mask[, j]
        if (!any(idx)) next

        rmse <- sqrt(mean((imputed[idx, j] - complete_mat[idx, j])^2))
        
        # Factor to scale RMSE below
        complete_vals <- complete_mat[, j]
        sd_j <- sd(complete_vals)
        mad_j <- mad(complete_vals)
        range_j <- diff(range(complete_vals))

        results[[length(results) + 1]] <- data.frame(
          method = method,
          feature = colnames(mat_na)[j],
          na_fraction = na_frac,
          rmse = rmse,
          nrmse_sd = if (sd_j > 0) rmse / sd_j else NA_real_,
          nrmse_mad = if (mad_j > 0) rmse / mad_j else NA_real_,
          nrmse_range = if (range_j > 0) rmse / range_j else NA_real_
        )
      }
    }
  }

  do.call(rbind, results)
}