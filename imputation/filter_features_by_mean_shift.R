#' Filter features with significantly altered means after imputation
#'
#' Removes features whose post-imputation mean deviates substantially
#' from the pre-imputation observed-value distribution, using a bootstrap
#' z-score approach.
#'
#' @param pre_imputation_matrix Numeric matrix or data.frame (samples × features) with NAs.
#' @param post_imputation_matrix Numeric matrix or data.frame (samples × features), fully imputed.
#' @param z_threshold Numeric. Absolute z-score cutoff (default 3).
#' @param n_boot Integer. Number of bootstrap resamples per feature (default 1000).
#' @param random_seed Integer for reproducibility.
#'
#' @return A list with:
#' \item{filtered_matrix}{Post-imputation matrix after removing flagged features.}
#' \item{removed_features}{Character vector of removed feature names.}
#' \item{summary}{Data frame with per-feature z-scores and filtering decisions.}
#'
#' @details
#' For each feature, the observed (non-NA) values are bootstrapped to build a
#' null distribution of means. The post-imputation mean is then expressed as a
#' z-score relative to this distribution. Features with \code{|z| > z_threshold}
#' are removed. Features with no missing values are always retained.
#'
#' @export
filter_features_by_mean_shift <- function(
    pre_imputation_matrix,
    post_imputation_matrix,
    z_threshold = 3,
    n_boot = 1000,
    random_seed = 1
) {
  stopifnot(is.matrix(pre_imputation_matrix) || is.data.frame(pre_imputation_matrix))
  stopifnot(is.matrix(post_imputation_matrix) || is.data.frame(post_imputation_matrix))

  pre <- as.matrix(pre_imputation_matrix)
  post <- as.matrix(post_imputation_matrix)

  stopifnot(identical(dim(pre), dim(post)))
  stopifnot(identical(colnames(pre), colnames(post)))

  set.seed(random_seed)

  # Identify features with shifted means
  remove_features <- c()
  summary_list <- list()

  for (feat in colnames(pre)) {

    obs <- pre[!is.na(pre[, feat]), feat]
    n_obs <- length(obs)
    n_na <- sum(is.na(pre[, feat]))

    if (n_na == 0) {
      summary_list[[feat]] <- data.frame(
        feature = feat,
        n_observed = n_obs,
        n_missing = 0L,
        pre_mean = mean(obs),
        post_mean = mean(post[, feat]),
        boot_mean = mean(obs),
        boot_sd = NA_real_,
        z_score = 0,
        removed = FALSE
      )
      next
    }

    # Bootstrap the observed-value mean
    boot_means <- vapply(seq_len(n_boot), function(b) {
      mean(sample(obs, n_obs, replace = TRUE))
    }, numeric(1))

    boot_mu <- mean(boot_means)
    boot_sd <- sd(boot_means)

    post_mean <- mean(post[, feat])

    z <- if (boot_sd > 0) (post_mean - boot_mu) / boot_sd else 0

    remove_flag <- abs(z) > z_threshold

    if (remove_flag) {
      remove_features <- c(remove_features, feat)
    }

    summary_list[[feat]] <- data.frame(
      feature = feat,
      n_observed = n_obs,
      n_missing = n_na,
      pre_mean = mean(obs),
      post_mean = post_mean,
      boot_mean = boot_mu,
      boot_sd = boot_sd,
      z_score = z,
      removed = remove_flag
    )
  }

  summary_df <- do.call(rbind, summary_list)

  filtered_mat <- post[, !colnames(post) %in% remove_features, drop = FALSE]

  list(
    filtered_matrix = filtered_mat,
    removed_features = remove_features,
    summary = summary_df
  )
}
