#' Filter features with significantly altered means after imputation
#'
#' Removes features whose post-imputation mean deviates substantially
#' from the pre-imputation observed-value distribution, using a bootstrap
#' z-score approach.
#'
#' @param pre_imputation_matrix Numeric matrix or data.frame (samples × features) with NAs.
#' @param post_imputation_matrix Numeric matrix or data.frame (samples × features), fully imputed.
#' @param grouping_vector Optional character/factor vector aligned to matrix rows.
#'   If provided, mean-shift filtering is performed separately within each group.
#' @param z_threshold Numeric. Absolute z-score cutoff (default 3).
#' @param n_boot Integer. Number of bootstrap resamples per feature (default 1000).
#' @param random_seed Integer for reproducibility.
#'
#' @return A list with:
#' \item{filtered_matrix}{Post-imputation matrix after removing flagged features.}
#' \item{removed_features}{Character vector of removed feature names.}
#' \item{summary}{Data frame with per-feature z-scores and filtering decisions.
#'   In grouped mode this contains one row per feature-group combination.}
#' \item{feature_summary}{Feature-level summary in grouped mode, with counts of
#'   groups tested and groups that failed.}
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
    grouping_vector = NULL,
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

  is_grouped <- !is.null(grouping_vector)
  if (is_grouped) {
    stopifnot(length(grouping_vector) == nrow(pre))
    if (any(is.na(grouping_vector))) {
      stop("grouping_vector cannot contain missing values")
    }
    grouping_vector <- as.character(grouping_vector)
  }

  set.seed(random_seed)

  # Identify features with shifted means
  remove_features <- c()
  summary_list <- list()
  feature_summary_list <- list()

  for (feat in colnames(pre)) {
    if (!is_grouped) {
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
      next
    }

    group_levels <- unique(grouping_vector)
    feat_summary_rows <- vector("list", length(group_levels))
    names(feat_summary_rows) <- group_levels

    for (grp in group_levels) {
      grp_mask <- grouping_vector == grp
      pre_grp <- pre[grp_mask, feat]
      post_grp <- post[grp_mask, feat]

      obs <- pre_grp[!is.na(pre_grp)]
      n_obs <- length(obs)
      n_na <- sum(is.na(pre_grp))

      if (n_obs == 0) {
        feat_summary_rows[[grp]] <- data.frame(
          feature = feat,
          grouping_level = grp,
          n_observed = 0L,
          n_missing = n_na,
          pre_mean = NA_real_,
          post_mean = mean(post_grp),
          boot_mean = NA_real_,
          boot_sd = NA_real_,
          z_score = NA_real_,
          group_removed = TRUE
        )
        next
      }

      if (n_na == 0) {
        feat_summary_rows[[grp]] <- data.frame(
          feature = feat,
          grouping_level = grp,
          n_observed = n_obs,
          n_missing = 0L,
          pre_mean = mean(obs),
          post_mean = mean(post_grp),
          boot_mean = mean(obs),
          boot_sd = NA_real_,
          z_score = 0,
          group_removed = FALSE
        )
        next
      }

      boot_means <- vapply(seq_len(n_boot), function(b) {
        mean(sample(obs, n_obs, replace = TRUE))
      }, numeric(1))

      boot_mu <- mean(boot_means)
      boot_sd <- sd(boot_means)

      post_mean <- mean(post_grp)

      z <- if (boot_sd > 0) (post_mean - boot_mu) / boot_sd else 0
      group_remove_flag <- abs(z) > z_threshold

      feat_summary_rows[[grp]] <- data.frame(
        feature = feat,
        grouping_level = grp,
        n_observed = n_obs,
        n_missing = n_na,
        pre_mean = mean(obs),
        post_mean = post_mean,
        boot_mean = boot_mu,
        boot_sd = boot_sd,
        z_score = z,
        group_removed = group_remove_flag
      )
    }

    feat_summary_df <- do.call(rbind, feat_summary_rows)
    n_groups_total <- nrow(feat_summary_df)
    n_groups_tested <- sum(feat_summary_df$n_missing > 0)
    n_groups_failed <- sum(feat_summary_df$group_removed)
    pct_groups_failed <- if (n_groups_total > 0) {
      100 * n_groups_failed / n_groups_total
    } else {
      NA_real_
    }
    feature_remove_flag <- n_groups_failed > 0

    feat_summary_df$feature_removed <- feature_remove_flag
    feat_summary_df$n_groups_total <- n_groups_total
    feat_summary_df$n_groups_tested <- n_groups_tested
    feat_summary_df$n_groups_failed <- n_groups_failed
    feat_summary_df$pct_groups_failed <- pct_groups_failed

    if (feature_remove_flag) {
      remove_features <- c(remove_features, feat)
    }

    summary_list[[feat]] <- feat_summary_df
    feature_summary_list[[feat]] <- data.frame(
      feature = feat,
      n_groups_total = n_groups_total,
      n_groups_tested = n_groups_tested,
      n_groups_failed = n_groups_failed,
      pct_groups_failed = pct_groups_failed,
      feature_removed = feature_remove_flag
    )
  }

  summary_df <- do.call(rbind, summary_list)
  feature_summary_df <- if (is_grouped) do.call(rbind, feature_summary_list) else NULL

  filtered_mat <- post[, !colnames(post) %in% remove_features, drop = FALSE]

  out <- list(
    filtered_matrix = filtered_mat,
    removed_features = remove_features,
    summary = summary_df
  )
  if (is_grouped) {
    out$feature_summary <- feature_summary_df
  }

  out
}
