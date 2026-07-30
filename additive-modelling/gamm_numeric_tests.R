# Generic inference and trend helpers for models returned by fit_gamm().
#
# The contrast construction follows tradeSeq's associationTest(),
# startVsEndTest(), and predictSmooth(): predictions are represented with the
# GAM linear-predictor matrix and differences are tested with the fitted
# coefficient covariance. Generic clustering of fitted patterns is implemented
# separately in utils/pattern_clustering.R. These functions operate on one
# continuous predictor and continuous or otherwise mgcv-compatible responses,
# without any SLE-specific column names.

.validate_feature_gamm <- function(fit) {
  if (!inherits(fit, "feature_gamm") || is.null(fit$full_model)) {
    stop("fit must be an object returned by fit_gamm().")
  }
  invisible(TRUE)
}

.gamm_reference_value <- function(x) {
  if (is.factor(x)) {
    observed <- as.character(x[!is.na(x)])
    value <- names(sort(table(observed), decreasing = TRUE))[[1L]]
    return(factor(value, levels = levels(x), ordered = is.ordered(x)))
  }
  if (is.character(x)) {
    return(names(sort(table(x), decreasing = TRUE))[[1L]])
  }
  if (is.logical(x)) return(names(sort(table(x), decreasing = TRUE))[[1L]] == "TRUE")
  if (is.numeric(x)) return(stats::median(x))
  x[[1L]]
}

# Construct prediction data while holding nuisance effects constant. Additive
# fixed effects cancel from all contrasts; the modal/median values simply give
# population-level fitted curves a reproducible vertical position.
gamm_prediction_data <- function(fit, predictor_values) {
  .validate_feature_gamm(fit)
  model_frame <- fit$full_model$model
  predictor_values <- as.numeric(predictor_values)
  newdata <- model_frame[rep(1L, length(predictor_values)), , drop = FALSE]

  for (variable in unique(c(fit$fixed_effects, fit$random_effects))) {
    newdata[[variable]] <- .gamm_reference_value(model_frame[[variable]])
  }
  newdata$.gamm_predictor <- predictor_values
  newdata$.gamm_response <- 0
  newdata
}

# Return an lpmatrix with random-intercept columns set to zero. This gives a
# population-level trajectory rather than the curve for an arbitrary subject.
gamm_population_lpmatrix <- function(fit, predictor_values) {
  .validate_feature_gamm(fit)
  newdata <- gamm_prediction_data(fit, predictor_values)
  design <- stats::predict(
    fit$full_model, newdata = newdata, type = "lpmatrix"
  )

  if (length(fit$random_effects) > 0L) {
    for (smooth in fit$full_model$smooth) {
      if (any(smooth$term %in% fit$random_effects)) {
        design[, smooth$first.para:smooth$last.para] <- 0
      }
    }
  }
  design
}

.gamm_wald_test <- function(estimate, covariance, tolerance = 1e-8) {
  estimate <- as.numeric(estimate)
  covariance <- (covariance + t(covariance)) / 2
  decomposition <- eigen(covariance, symmetric = TRUE)
  largest <- max(abs(decomposition$values), 0)
  retained <- if (largest > 0) {
    decomposition$values > tolerance * largest
  } else {
    rep(FALSE, length(decomposition$values))
  }
  degrees_freedom <- sum(retained)

  if (degrees_freedom == 0L) {
    return(c(wald_statistic = NA_real_, df = 0, p_value = NA_real_))
  }
  projected <- crossprod(decomposition$vectors[, retained, drop = FALSE], estimate)
  statistic <- sum(projected^2 / decomposition$values[retained])
  c(
    wald_statistic = statistic,
    df = degrees_freedom,
    p_value = stats::pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  )
}

# Omnibus test of whether the fitted trajectory changes over the predictor.
# As in tradeSeq associationTest(), several locations are contrasted with the
# start and tested jointly; this is not merely the p-value of one coefficient.
gamm_association_test <- function(fit, n_points = NULL, tolerance = 1e-8) {
  .validate_feature_gamm(fit)
  if (is.null(n_points)) n_points <- max(2L * fit$effective_k, 10L)
  n_points <- as.integer(n_points)
  if (n_points < 3L) stop("n_points must be at least 3.")

  predictor_values <- seq(
    fit$predictor_range[[1L]], fit$predictor_range[[2L]],
    length.out = n_points
  )
  design <- gamm_population_lpmatrix(fit, predictor_values)
  contrasts <- design[-1L, , drop = FALSE] -
    design[rep(1L, n_points - 1L), , drop = FALSE]
  coefficients <- stats::coef(fit$full_model)
  coefficient_covariance <- stats::vcov(fit$full_model)
  differences <- contrasts %*% coefficients
  contrast_covariance <- contrasts %*% coefficient_covariance %*% t(contrasts)
  wald <- .gamm_wald_test(differences, contrast_covariance, tolerance)

  smooth_table <- summary(fit$full_model)$s.table
  predictor_row <- grep(
    "^s\\(.gamm_predictor\\)", rownames(smooth_table), value = TRUE
  )
  if (length(predictor_row) == 1L) {
    smooth_values <- smooth_table[predictor_row, , drop = TRUE]
    smooth_edf <- unname(smooth_values[["edf"]])
    smooth_reference_df <- unname(smooth_values[["Ref.df"]])
    smooth_statistic <- unname(smooth_values[[ncol(smooth_table) - 1L]])
    smooth_p_value <- unname(smooth_values[[ncol(smooth_table)]])
  } else {
    smooth_edf <- smooth_reference_df <- smooth_statistic <- smooth_p_value <- NA_real_
  }

  data.frame(
    association_wald_statistic = unname(wald[["wald_statistic"]]),
    association_df = unname(wald[["df"]]),
    association_p_value = unname(wald[["p_value"]]),
    association_mean_abs_change = mean(abs(differences)),
    smooth_edf = smooth_edf,
    smooth_reference_df = smooth_reference_df,
    smooth_statistic = smooth_statistic,
    smooth_p_value = smooth_p_value,
    stringsAsFactors = FALSE
  )
}

# Wald contrast between two predictor locations. By default these are the
# observed start and end for that feature, matching tradeSeq startVsEndTest().
gamm_start_vs_end_test <- function(
    fit,
    predictor_values = NULL,
    tolerance = 1e-12
) {
  .validate_feature_gamm(fit)
  if (is.null(predictor_values)) predictor_values <- fit$predictor_range
  if (length(predictor_values) != 2L || any(!is.finite(predictor_values))) {
    stop("predictor_values must contain two finite values.")
  }
  if (predictor_values[[1L]] >= predictor_values[[2L]]) {
    stop("The start predictor value must be smaller than the end value.")
  }

  design <- gamm_population_lpmatrix(fit, predictor_values)
  contrast <- design[2L, , drop = FALSE] - design[1L, , drop = FALSE]
  estimate <- as.numeric(contrast %*% stats::coef(fit$full_model))
  covariance <- stats::vcov(fit$full_model)
  variance <- as.numeric(contrast %*% covariance %*% t(contrast))
  standard_error <- if (is.finite(variance) && variance > tolerance) {
    sqrt(variance)
  } else {
    NA_real_
  }
  z_value <- estimate / standard_error
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)

  data.frame(
    start_predictor = predictor_values[[1L]],
    end_predictor = predictor_values[[2L]],
    start_vs_end_change = estimate,
    start_vs_end_se = standard_error,
    start_vs_end_z = z_value,
    start_vs_end_p_value = p_value,
    stringsAsFactors = FALSE
  )
}

# Population-level fitted trajectory. Random intercepts are excluded. When
# predictor_values is supplied, callers can use a shared predictor grid across
# features; values outside a feature's observed range are model extrapolations.
gamm_predict_smooth <- function(
    fit,
    n_points = 100L,
    predictor_values = NULL,
    original_response_scale = FALSE
) {
  .validate_feature_gamm(fit)
  if (is.null(predictor_values)) {
    predictor_position <- seq(0, 1, length.out = n_points)
    predictor_values <- fit$predictor_range[[1L]] + predictor_position *
      diff(fit$predictor_range)
  } else {
    predictor_values <- as.numeric(predictor_values)
    predictor_position <- (predictor_values - fit$predictor_range[[1L]]) /
      diff(fit$predictor_range)
  }

  design <- gamm_population_lpmatrix(fit, predictor_values)
  fitted <- as.numeric(design %*% stats::coef(fit$full_model))
  standard_error <- sqrt(rowSums((design %*% stats::vcov(fit$full_model)) * design))
  if (original_response_scale) {
    fitted <- fitted * fit$response_scale + fit$response_center
    standard_error <- standard_error * fit$response_scale
  }

  data.frame(
    grid_index = seq_along(predictor_values),
    predictor_position = predictor_position,
    predictor = predictor_values,
    fitted = fitted,
    standard_error = standard_error,
    lower_95 = fitted - stats::qnorm(0.975) * standard_error,
    upper_95 = fitted + stats::qnorm(0.975) * standard_error,
    stringsAsFactors = FALSE
  )
}
