# Generic inference and trend helpers for models returned by fit_gamm().
#
# The contrast construction follows tradeSeq's associationTest(),
# startVsEndTest(), predictSmooth(), and clusterExpressionPatterns(): predictions
# are represented with the GAM linear-predictor matrix; differences are tested
# with the fitted coefficient covariance; fitted trajectories are row-scaled
# before clustering. These functions operate on one continuous predictor and
# continuous or otherwise mgcv-compatible responses, without any SLE-specific
# column names.

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

# Population-level fitted trajectory. Random intercepts are excluded. Relative
# positions make curves with feature-specific available ranges clusterable
# without extrapolating beyond a feature's observations.
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

# Hierarchical clustering of row-standardized fitted trajectories, mirroring
# tradeSeq's central idea while avoiding a dependency on clusterExperiment.
cluster_gamm_trends <- function(
    trend_data,
    feature_col = "feature",
    value_col = "fitted",
    grid_col = "grid_index",
    n_clusters = 6L,
    method = "ward.D2"
) {
  required <- c(feature_col, value_col, grid_col)
  missing_cols <- setdiff(required, colnames(trend_data))
  if (length(missing_cols) > 0L) {
    stop("Missing trend columns: ", paste(missing_cols, collapse = ", "))
  }

  features <- unique(as.character(trend_data[[feature_col]]))
  grid <- sort(unique(trend_data[[grid_col]]))
  trend_matrix <- matrix(
    NA_real_, nrow = length(features), ncol = length(grid),
    dimnames = list(features, as.character(grid))
  )
  for (feature in features) {
    rows <- trend_data[as.character(trend_data[[feature_col]]) == feature, ]
    trend_matrix[feature, match(rows[[grid_col]], grid)] <- rows[[value_col]]
  }
  if (anyNA(trend_matrix)) stop("Each feature must have a value at every grid point.")

  scaled_matrix <- t(scale(t(trend_matrix)))
  usable <- rowSums(is.finite(scaled_matrix)) == ncol(scaled_matrix)
  if (!all(usable)) {
    warning(sum(!usable), " constant fitted trajectories were excluded from clustering.")
    scaled_matrix <- scaled_matrix[usable, , drop = FALSE]
  }
  if (nrow(scaled_matrix) == 0L) stop("No non-constant trajectories to cluster.")

  if (nrow(scaled_matrix) == 1L) {
    clustering <- NULL
    cluster <- stats::setNames(1L, rownames(scaled_matrix))
    order <- rownames(scaled_matrix)
  } else {
    clustering <- stats::hclust(stats::dist(scaled_matrix), method = method)
    effective_clusters <- min(as.integer(n_clusters), nrow(scaled_matrix))
    if (effective_clusters < 1L) stop("n_clusters must be at least one.")
    cluster <- stats::cutree(clustering, k = effective_clusters)
    order <- clustering$labels[clustering$order]
  }

  assignments <- data.frame(
    feature = names(cluster),
    trend_cluster = as.integer(cluster),
    heatmap_order = match(names(cluster), order),
    stringsAsFactors = FALSE
  )
  scaled_long <- data.frame(
    feature = rep(rownames(scaled_matrix), each = ncol(scaled_matrix)),
    grid_index = rep(grid, times = nrow(scaled_matrix)),
    scaled_trend = as.vector(t(scaled_matrix)),
    stringsAsFactors = FALSE
  )
  scaled_long <- merge(scaled_long, assignments, by = "feature", sort = FALSE)

  list(
    assignments = assignments,
    scaled_trends = scaled_long,
    scaled_matrix = scaled_matrix,
    hclust = clustering
  )
}

# Evaluate candidate cuts of the same hierarchical tree used for trend
# clustering. Within-cluster sum of squares supplies an elbow diagnostic, while
# average silhouette width measures separation and cohesion. The elbow is the
# candidate with greatest perpendicular distance below the line joining the
# first and last points of the normalized WSS curve.
evaluate_gamm_trend_clusters <- function(
    trend_data,
    feature_col = "feature",
    value_col = "fitted",
    grid_col = "grid_index",
    max_clusters = 10L,
    method = "ward.D2"
) {
  prepared <- cluster_gamm_trends(
    trend_data = trend_data,
    feature_col = feature_col,
    value_col = value_col,
    grid_col = grid_col,
    n_clusters = 1L,
    method = method
  )
  trend_matrix <- prepared$scaled_matrix
  n_features <- nrow(trend_matrix)
  if (n_features < 2L) {
    metrics <- data.frame(
      k = 1L,
      within_cluster_ss = 0,
      explained_between_fraction = NA_real_,
      average_silhouette = NA_real_
    )
    return(list(
      metrics = metrics,
      elbow_recommended_k = 1L,
      silhouette_recommended_k = 1L
    ))
  }

  maximum_k <- min(as.integer(max_clusters), n_features)
  if (!is.finite(maximum_k) || maximum_k < 2L) {
    stop("max_clusters must be an integer of at least two.")
  }
  clustering <- prepared$hclust
  distance_matrix <- as.matrix(stats::dist(trend_matrix))
  total_ss <- sum(sweep(
    trend_matrix, 2L, colMeans(trend_matrix), FUN = "-"
  )^2)

  calculate_wss <- function(cluster) {
    sum(vapply(split(seq_len(n_features), cluster), function(indices) {
      cluster_matrix <- trend_matrix[indices, , drop = FALSE]
      sum(sweep(cluster_matrix, 2L, colMeans(cluster_matrix), FUN = "-")^2)
    }, numeric(1)))
  }

  calculate_silhouette <- function(cluster) {
    if (length(unique(cluster)) < 2L) return(NA_real_)
    silhouette <- vapply(seq_len(n_features), function(i) {
      same_cluster <- which(cluster == cluster[[i]])
      same_cluster <- setdiff(same_cluster, i)
      if (length(same_cluster) == 0L) return(0)
      within_distance <- mean(distance_matrix[i, same_cluster])
      other_clusters <- setdiff(unique(cluster), cluster[[i]])
      nearest_other_distance <- min(vapply(other_clusters, function(other) {
        mean(distance_matrix[i, cluster == other])
      }, numeric(1)))
      (nearest_other_distance - within_distance) /
        max(within_distance, nearest_other_distance)
    }, numeric(1))
    mean(silhouette)
  }

  metrics <- do.call(rbind, lapply(seq_len(maximum_k), function(k) {
    cluster <- stats::cutree(clustering, k = k)
    wss <- calculate_wss(cluster)
    data.frame(
      k = k,
      within_cluster_ss = wss,
      explained_between_fraction = if (total_ss > 0) 1 - wss / total_ss else NA_real_,
      average_silhouette = calculate_silhouette(cluster)
    )
  }))

  if (nrow(metrics) <= 2L) {
    elbow_k <- metrics$k[[nrow(metrics)]]
  } else {
    normalized_k <- (metrics$k - min(metrics$k)) /
      (max(metrics$k) - min(metrics$k))
    wss_range <- range(metrics$within_cluster_ss)
    normalized_wss <- if (diff(wss_range) > 0) {
      (metrics$within_cluster_ss - wss_range[[1L]]) / diff(wss_range)
    } else {
      rep(0, nrow(metrics))
    }
    x1 <- normalized_k[[1L]]
    y1 <- normalized_wss[[1L]]
    x2 <- normalized_k[[nrow(metrics)]]
    y2 <- normalized_wss[[nrow(metrics)]]
    distance_from_line <- abs(
      (y2 - y1) * normalized_k - (x2 - x1) * normalized_wss +
        x2 * y1 - y2 * x1
    ) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
    distance_from_line[c(1L, length(distance_from_line))] <- -Inf
    elbow_k <- metrics$k[[which.max(distance_from_line)]]
  }

  silhouette_candidates <- which(is.finite(metrics$average_silhouette))
  silhouette_k <- if (length(silhouette_candidates) > 0L) {
    candidate <- silhouette_candidates[
      which.max(metrics$average_silhouette[silhouette_candidates])
    ]
    metrics$k[[candidate]]
  } else {
    elbow_k
  }

  list(
    metrics = metrics,
    elbow_recommended_k = as.integer(elbow_k),
    silhouette_recommended_k = as.integer(silhouette_k)
  )
}
