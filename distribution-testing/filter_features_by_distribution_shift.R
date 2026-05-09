#' Compute distribution-shift statistics for a feature
#'
#' @param pre_vec Numeric vector containing pre-imputation values with NAs.
#' @param post_vec Numeric vector containing post-imputation values.
#' @param z_threshold Numeric. Absolute z-score cutoff.
#' @param ks_alpha Numeric. Two-sample KS test p-value cutoff.
#' @param min_observed Integer. Minimum observed values required for testing.
#' @param sd_floor Numeric. Small SD used when observed values have no variation.
#'
#' @return One-row data.frame containing distribution-shift statistics.
#'
#' @keywords internal
flag_features <- function(
    pre_vec,
    post_vec,
    z_threshold,
    ks_alpha,
    min_observed = 3,
    sd_floor = 1e-6
) {

  obs <- pre_vec[!is.na(pre_vec)]

  n_obs <- length(obs)
  n_na <- sum(is.na(pre_vec))
  pre_mean <- if (n_obs > 0) mean(obs) else NA_real_
  post_mean <- mean(post_vec)
  pre_sd <- if (n_obs > 1) sd(obs) else NA_real_

  out <- data.frame(
    n_observed = n_obs,
    n_missing = n_na,
    pre_mean = pre_mean,
    post_mean = post_mean,
    pre_sd = pre_sd,
    sd_used = NA_real_,
    z_score = NA_real_,
    z_shift_flag = NA,
    ks_statistic = NA_real_,
    ks_p_value = NA_real_,
    ks_flag = NA,
    removed = NA,
    test_performed = FALSE,
    reason = NA_character_
  )

  # 1. Features without imputed values are always retained.
  if (n_na == 0) {
    out$sd_used <- if (!is.na(pre_sd) && pre_sd == 0) sd_floor else pre_sd
    out$z_score <- 0
    out$z_shift_flag <- FALSE
    out$ks_flag <- FALSE
    out$removed <- FALSE
    out$reason <- "no_missing_values"
    return(out)
  }

  # 2. Features with missingness need enough observed values for comparison.
  if (n_obs < min_observed) {
    out$reason <- "insufficient_observed_values"
    out$removed <- TRUE
    return(out)
  }

  # 3. Features eligible for testing
  # Use a tiny SD for constant observed values so z-score shifts remain testable.
  out$sd_used <- if (pre_sd == 0) sd_floor else pre_sd
  out$z_score <- (post_mean - pre_mean) / out$sd_used
  out$z_shift_flag <- abs(out$z_score) > z_threshold

  # Compare the observed pre-imputation distribution to the post-imputation one.
  ks_res <- suppressWarnings(stats::ks.test(obs, post_vec))
  out$ks_statistic <- unname(ks_res$statistic)
  out$ks_p_value <- ks_res$p.value
  out$ks_flag <- out$ks_p_value < ks_alpha
  out$removed <- isTRUE(out$z_shift_flag) || isTRUE(out$ks_flag)
  out$test_performed <- TRUE

  if (out$removed && out$z_shift_flag && out$ks_flag) {
    out$reason <- "ks_and_z_shift"
  } else if (out$removed && out$ks_flag) {
    out$reason <- "ks_distribution_shift"
  } else if (out$removed && out$z_shift_flag) {
    out$reason <- "z_shift"
  } else {
    out$reason <- "passed"
  }

  stopifnot(all(!is.na(out$removed)))
  out
}


#' Filter features with distribution shifts after imputation
#'
#' Removes features whose post-imputation distribution differs substantially
#' from the pre-imputation observed values using a two-sample KS test or
#' z-score shift.
#'
#' @param pre_imputation_matrix Numeric matrix or data.frame (samples × features) with NAs.
#' @param post_imputation_matrix Numeric matrix or data.frame (samples × features), fully imputed.
#' @param grouping_vector Optional character/factor vector aligned to matrix rows.
#'   If provided, filtering is performed separately within each group.
#' @param z_threshold Numeric. Absolute z-score shift cutoff (default 0.1).
#' @param ks_alpha Numeric. Two-sample KS test p-value cutoff (default 0.05).
#' @param min_observed Integer. Minimum observed values required for testing.
#' @param sd_floor Numeric. Small SD used when observed values have no variation.
#'
#' @return A list containing:
#' \item{filtered_matrix}{Post-imputation matrix after feature removal.}
#' \item{removed_features}{Character vector of removed feature names.}
#' \item{summary}{Per-feature statistics table.}
#' \item{group_summary}{Feature-level grouped summary (grouped mode only).}
#'
#' @details
#' For each feature with missing values, observed pre-imputation values are
#' compared against post-imputation values. Features are removed when the
#' two-sample KS test has \code{p < ks_alpha} or the absolute z-score shift is
#' greater than \code{z_threshold}.
#'
#' Features with no missing values are always retained.
#'
#' @export
filter_features_by_distribution_shift <- function(
    pre_imputation_matrix,
    post_imputation_matrix,
    grouping_vector = NULL,
    z_threshold = 0.1,
    ks_alpha = 0.05,
    min_observed = 3,
    sd_floor = 1e-6
) {

  # Validate inputs
  stopifnot(
    is.matrix(pre_imputation_matrix) || is.data.frame(pre_imputation_matrix),
    is.matrix(post_imputation_matrix) || is.data.frame(post_imputation_matrix)
  )

  pre <- as.matrix(pre_imputation_matrix)
  post <- as.matrix(post_imputation_matrix)
  stopifnot(is.numeric(pre), is.numeric(post))

  stopifnot(
    identical(dim(pre), dim(post)),
    identical(colnames(pre), colnames(post))
  )

  # Validate grouping structure
  is_grouped <- !is.null(grouping_vector)
  if (is_grouped) {
    stopifnot(length(grouping_vector) == nrow(pre))
    if (any(is.na(grouping_vector))) {
      stop("grouping_vector cannot contain missing values")
    }
    grouping_vector <- as.character(grouping_vector)
  }

  summary_list <- list()
  group_summary_list <- list()
  remove_features <- character()

  for (feat in colnames(pre)) {

    # Non-grouped analysis
    if (!is_grouped) {
      res <- flag_features(
        pre_vec = pre[, feat],
        post_vec = post[, feat],
        z_threshold = z_threshold,
        ks_alpha = ks_alpha,
        min_observed = min_observed,
        sd_floor = sd_floor
      )
      res$feature <- feat
      if (isTRUE(res$removed)) {
        remove_features <- c(remove_features, feat)
      }
      summary_list[[feat]] <- res
      next
    }

    # Grouped analysis
    feat_summary_df <- do.call(
      rbind,
      lapply(unique(grouping_vector), function(grp) {

        grp_mask <- grouping_vector == grp

        res <- flag_features(
          pre_vec = pre[grp_mask, feat],
          post_vec = post[grp_mask, feat],
          z_threshold = z_threshold,
          ks_alpha = ks_alpha,
          min_observed = min_observed,
          sd_floor = sd_floor
        )

        res$feature <- feat
        res$grouping_level <- grp

        res
      })
    )

    # Summarise grouped failures
    n_groups_total <- nrow(feat_summary_df)
    n_groups_tested <- sum(feat_summary_df$test_performed)
    n_groups_failed <- sum(
      feat_summary_df$removed, na.rm = TRUE
    )
    pct_groups_failed <- 100 * n_groups_failed / n_groups_total
    feature_removed <- any(
      feat_summary_df$removed, na.rm = TRUE
    )
    feat_summary_df$feature_removed <- feature_removed
    feat_summary_df$n_groups_total <- n_groups_total
    feat_summary_df$n_groups_tested <- n_groups_tested
    feat_summary_df$n_groups_failed <- n_groups_failed
    feat_summary_df$pct_groups_failed <- pct_groups_failed

    if (feature_removed) {
      remove_features <- c(remove_features, feat)
    }

    summary_list[[feat]] <- feat_summary_df

    group_summary_list[[feat]] <- data.frame(
      feature = feat,
      n_groups_total = n_groups_total,
      n_groups_tested = n_groups_tested,
      n_groups_failed = n_groups_failed,
      pct_groups_failed = pct_groups_failed,
      feature_removed = feature_removed
    )
  }

  # Combine outputs
  summary_df <- do.call(rbind, summary_list)

  group_summary_df <- if (is_grouped) {
    do.call(rbind, group_summary_list)
  } else {
    NULL
  }

  remove_features <- unique(remove_features)
  filtered_matrix <- post[
    ,
    !colnames(post) %in% remove_features,
    drop = FALSE
  ]

  out <- list(
    filtered_matrix = filtered_matrix,
    removed_features = remove_features,
    summary = summary_df
  )

  if (is_grouped) {
    out$group_summary <- group_summary_df
  }

  out
}
