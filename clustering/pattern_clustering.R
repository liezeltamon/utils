# Generic clustering helpers for ordered numeric patterns.
#
# Input matrices always have items to cluster in rows and ordered positions or
# measurements in columns. The functions contain no trajectory-, GAMM-, or
# project-specific assumptions, so they can be reused for other time-course or
# profile-clustering tasks.

validate_pattern_matrix <- function(pattern_matrix) {
  pattern_matrix <- as.matrix(pattern_matrix)
  if (!is.numeric(pattern_matrix)) {
    stop("pattern_matrix must be numeric.")
  }
  if (length(dim(pattern_matrix)) != 2L ||
      nrow(pattern_matrix) < 1L ||
      ncol(pattern_matrix) < 2L) {
    stop("pattern_matrix must contain at least one row and two columns.")
  }
  if (is.null(rownames(pattern_matrix)) ||
      anyNA(rownames(pattern_matrix)) ||
      any(rownames(pattern_matrix) == "") ||
      anyDuplicated(rownames(pattern_matrix))) {
    stop("pattern_matrix must have unique, non-empty row names.")
  }
  if (any(!is.finite(pattern_matrix))) {
    stop("pattern_matrix must contain only finite values.")
  }
  if (is.null(colnames(pattern_matrix))) {
    colnames(pattern_matrix) <- as.character(seq_len(ncol(pattern_matrix)))
  }
  pattern_matrix
}

long_patterns_to_matrix <- function(
    pattern_data,
    item_col = "item",
    value_col = "value",
    position_col = "position"
) {
  required <- c(item_col, value_col, position_col)
  missing_cols <- setdiff(required, colnames(pattern_data))
  if (length(missing_cols) > 0L) {
    stop("Missing pattern column(s): ", paste(missing_cols, collapse = ", "))
  }

  item <- as.character(pattern_data[[item_col]])
  position <- pattern_data[[position_col]]
  value <- suppressWarnings(as.numeric(pattern_data[[value_col]]))
  if (anyNA(item) || any(item == "")) {
    stop("Pattern item identifiers must be non-empty.")
  }
  if (anyNA(position)) {
    stop("Pattern positions must not be missing.")
  }
  if (any(!is.finite(value))) {
    stop("Pattern values must be finite.")
  }

  cell_id <- paste(item, as.character(position), sep = "\r")
  if (anyDuplicated(cell_id)) {
    stop("Pattern data contain duplicate item-position combinations.")
  }

  items <- unique(item)
  positions <- sort(unique(position))
  pattern_matrix <- matrix(
    NA_real_,
    nrow = length(items),
    ncol = length(positions),
    dimnames = list(items, as.character(positions))
  )
  pattern_matrix[
    cbind(match(item, items), match(position, positions))
  ] <- value
  if (anyNA(pattern_matrix)) {
    stop("Each item must have a value at every pattern position.")
  }
  validate_pattern_matrix(pattern_matrix)
}

scale_pattern_rows <- function(pattern_matrix, drop_constant = TRUE) {
  pattern_matrix <- validate_pattern_matrix(pattern_matrix)
  row_centres <- rowMeans(pattern_matrix)
  row_scales <- apply(pattern_matrix, 1L, stats::sd)
  usable <- is.finite(row_scales) & row_scales > 0
  excluded_items <- rownames(pattern_matrix)[!usable]

  if (length(excluded_items) > 0L && !isTRUE(drop_constant)) {
    stop(
      "Constant pattern row(s) cannot be standardized: ",
      paste(utils::head(excluded_items, 10L), collapse = ", ")
    )
  }
  if (length(excluded_items) > 0L) {
    warning(
      length(excluded_items),
      " constant pattern(s) were excluded from clustering."
    )
  }

  pattern_matrix <- pattern_matrix[usable, , drop = FALSE]
  if (nrow(pattern_matrix) == 0L) {
    stop("No non-constant patterns remain for clustering.")
  }
  scaled_matrix <- sweep(pattern_matrix, 1L, row_centres[usable], FUN = "-")
  scaled_matrix <- sweep(scaled_matrix, 1L, row_scales[usable], FUN = "/")

  list(
    scaled_matrix = scaled_matrix,
    excluded_items = excluded_items,
    row_centres = row_centres[usable],
    row_scales = row_scales[usable]
  )
}

pattern_matrix_to_long <- function(
    pattern_matrix,
    item_col = "item",
    position_col = "position",
    value_col = "value"
) {
  pattern_matrix <- validate_pattern_matrix(pattern_matrix)
  out <- data.frame(
    item = rep(rownames(pattern_matrix), each = ncol(pattern_matrix)),
    position = rep(colnames(pattern_matrix), times = nrow(pattern_matrix)),
    value = as.vector(t(pattern_matrix)),
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(item_col, position_col, value_col)
  out
}

calculate_pattern_silhouette <- function(
    pattern_matrix,
    assignments,
    unassigned_values = -1L
) {
  pattern_matrix <- validate_pattern_matrix(pattern_matrix)
  if (is.null(names(assignments))) {
    if (length(assignments) != nrow(pattern_matrix)) {
      stop("Unnamed assignments must match the number of pattern rows.")
    }
    names(assignments) <- rownames(pattern_matrix)
  }
  missing_assignments <- setdiff(rownames(pattern_matrix), names(assignments))
  if (length(missing_assignments) > 0L) {
    stop(
      "Missing cluster assignment(s) for: ",
      paste(utils::head(missing_assignments, 10L), collapse = ", ")
    )
  }

  assignments <- assignments[rownames(pattern_matrix)]
  keep <- !is.na(assignments) & !(assignments %in% unassigned_values)
  kept_matrix <- pattern_matrix[keep, , drop = FALSE]
  kept_assignments <- assignments[keep]
  result <- data.frame(
    item = rownames(pattern_matrix),
    cluster = as.integer(assignments),
    silhouette = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(kept_matrix) < 2L ||
      length(unique(kept_assignments)) < 2L) {
    return(list(values = result, average = NA_real_))
  }

  distance_matrix <- as.matrix(stats::dist(kept_matrix))
  silhouettes <- vapply(seq_len(nrow(kept_matrix)), function(index) {
    same_cluster <- which(kept_assignments == kept_assignments[[index]])
    same_cluster <- setdiff(same_cluster, index)
    if (length(same_cluster) == 0L) {
      return(0)
    }
    within_distance <- mean(distance_matrix[index, same_cluster])
    other_clusters <- setdiff(
      unique(kept_assignments),
      kept_assignments[[index]]
    )
    nearest_other_distance <- min(vapply(other_clusters, function(cluster_id) {
      mean(distance_matrix[index, kept_assignments == cluster_id])
    }, numeric(1)))
    denominator <- max(within_distance, nearest_other_distance)
    if (denominator == 0) {
      0
    } else {
      (nearest_other_distance - within_distance) / denominator
    }
  }, numeric(1))

  result$silhouette[match(rownames(kept_matrix), result$item)] <- silhouettes
  list(values = result, average = mean(silhouettes))
}

cluster_patterns_hierarchical <- function(
    pattern_matrix,
    n_clusters = 6L,
    linkage = "ward.D2",
    scale_rows = TRUE
) {
  pattern_matrix <- validate_pattern_matrix(pattern_matrix)
  prepared <- if (isTRUE(scale_rows)) {
    scale_pattern_rows(pattern_matrix)
  } else {
    list(scaled_matrix = pattern_matrix, excluded_items = character())
  }
  scaled_matrix <- prepared$scaled_matrix

  n_clusters <- as.integer(n_clusters)
  if (length(n_clusters) != 1L ||
      !is.finite(n_clusters) ||
      n_clusters < 1L) {
    stop("n_clusters must be one positive integer.")
  }
  effective_clusters <- min(n_clusters, nrow(scaled_matrix))

  if (nrow(scaled_matrix) == 1L) {
    clustering <- NULL
    cluster <- stats::setNames(1L, rownames(scaled_matrix))
    item_order <- rownames(scaled_matrix)
  } else {
    clustering <- stats::hclust(
      stats::dist(scaled_matrix),
      method = linkage
    )
    cluster <- stats::cutree(clustering, k = effective_clusters)
    item_order <- clustering$labels[clustering$order]
  }

  assignments <- data.frame(
    item = names(cluster),
    cluster = as.integer(cluster),
    item_order = match(names(cluster), item_order),
    clustering_method = "hierarchical",
    stringsAsFactors = FALSE
  )

  list(
    assignments = assignments,
    scaled_matrix = scaled_matrix,
    item_order = item_order,
    excluded_items = prepared$excluded_items,
    model = clustering,
    co_clustering = NULL,
    silhouette = calculate_pattern_silhouette(
      scaled_matrix,
      stats::setNames(assignments$cluster, assignments$item)
    )
  )
}

evaluate_hierarchical_pattern_clusters <- function(
    pattern_matrix,
    max_clusters = 10L,
    linkage = "ward.D2",
    scale_rows = TRUE
) {
  prepared <- cluster_patterns_hierarchical(
    pattern_matrix = pattern_matrix,
    n_clusters = 1L,
    linkage = linkage,
    scale_rows = scale_rows
  )
  scaled_matrix <- prepared$scaled_matrix
  n_items <- nrow(scaled_matrix)
  if (n_items < 2L) {
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

  maximum_k <- min(as.integer(max_clusters), n_items)
  if (length(maximum_k) != 1L ||
      !is.finite(maximum_k) ||
      maximum_k < 2L) {
    stop("max_clusters must be one integer of at least two.")
  }
  clustering <- prepared$model
  total_ss <- sum(sweep(
    scaled_matrix,
    2L,
    colMeans(scaled_matrix),
    FUN = "-"
  )^2)

  metrics <- do.call(rbind, lapply(seq_len(maximum_k), function(k) {
    cluster <- stats::cutree(clustering, k = k)
    within_cluster_ss <- sum(vapply(
      split(seq_len(n_items), cluster),
      function(indices) {
        cluster_matrix <- scaled_matrix[indices, , drop = FALSE]
        sum(sweep(
          cluster_matrix,
          2L,
          colMeans(cluster_matrix),
          FUN = "-"
        )^2)
      },
      numeric(1)
    ))
    silhouette <- calculate_pattern_silhouette(
      scaled_matrix,
      cluster
    )
    data.frame(
      k = k,
      within_cluster_ss = within_cluster_ss,
      explained_between_fraction = if (total_ss > 0) {
        1 - within_cluster_ss / total_ss
      } else {
        NA_real_
      },
      average_silhouette = silhouette$average
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
      (y2 - y1) * normalized_k -
        (x2 - x1) * normalized_wss +
        x2 * y1 -
        y2 * x1
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

complete_co_clustering_matrix <- function(co_clustering, item_names) {
  if (is.null(co_clustering) ||
      (length(co_clustering) == 1L && is.na(co_clustering))) {
    return(NULL)
  }
  co_clustering <- as.matrix(co_clustering)
  if (!all(dim(co_clustering) == length(item_names))) {
    stop("RSEC co-clustering matrix dimensions do not match the items.")
  }
  if (is.null(rownames(co_clustering))) {
    rownames(co_clustering) <- item_names
  }
  if (is.null(colnames(co_clustering))) {
    colnames(co_clustering) <- item_names
  }
  co_clustering <- co_clustering[item_names, item_names, drop = FALSE]

  upper_values <- co_clustering[upper.tri(co_clustering)]
  lower_values <- co_clustering[lower.tri(co_clustering)]
  if (all(is.na(upper_values) | upper_values == 0) &&
      any(is.finite(lower_values) & lower_values != 0)) {
    co_clustering[upper.tri(co_clustering)] <-
      t(co_clustering)[upper.tri(co_clustering)]
  } else if (all(is.na(lower_values) | lower_values == 0) &&
             any(is.finite(upper_values) & upper_values != 0)) {
    co_clustering[lower.tri(co_clustering)] <-
      t(co_clustering)[lower.tri(co_clustering)]
  } else {
    co_clustering <- (co_clustering + t(co_clustering)) / 2
  }
  diag(co_clustering) <- 1
  co_clustering
}

cluster_patterns_rsec <- function(
    pattern_matrix,
    scale_rows = TRUE,
    reduce_method = "PCA",
    n_reduced_dims = 10L,
    k0s = 4:15,
    resamples = 100L,
    subsample_fraction = 0.7,
    subsample_cluster_function = "pam",
    main_cluster_function = "hierarchical01",
    alphas = c(0.1, 0.2, 0.3),
    beta = 0.9,
    consensus_proportion = 0.7,
    min_cluster_size = 6L,
    remaining_items_to_stop = 30L,
    random_seed = 176201L,
    n_cores = 1L,
    verbose = FALSE
) {
  if (!requireNamespace("clusterExperiment", quietly = TRUE)) {
    stop(
      "Package 'clusterExperiment' is required for RSEC clustering."
    )
  }

  pattern_matrix <- validate_pattern_matrix(pattern_matrix)
  prepared <- if (isTRUE(scale_rows)) {
    scale_pattern_rows(pattern_matrix)
  } else {
    list(scaled_matrix = pattern_matrix, excluded_items = character())
  }
  scaled_matrix <- prepared$scaled_matrix
  remaining_items_to_stop <- as.integer(remaining_items_to_stop)
  if (length(remaining_items_to_stop) != 1L ||
      !is.finite(remaining_items_to_stop) ||
      remaining_items_to_stop < 0L) {
    stop("remaining_items_to_stop must be one non-negative integer.")
  }
  if (nrow(scaled_matrix) <= remaining_items_to_stop) {
    stop(
      "RSEC sequential clustering requires more than ",
      remaining_items_to_stop,
      " non-constant patterns with the current remaining_items_to_stop ",
      "setting; received ",
      nrow(scaled_matrix),
      "."
    )
  }
  effective_reduced_dims <- min(
    as.integer(n_reduced_dims),
    nrow(scaled_matrix) - 1L,
    ncol(scaled_matrix)
  )
  if (!is.finite(effective_reduced_dims) || effective_reduced_dims < 1L) {
    stop("RSEC requires at least two non-constant patterns.")
  }

  rsec <- clusterExperiment::RSEC(
    t(scaled_matrix),
    isCount = FALSE,
    transFun = NULL,
    reduceMethod = reduce_method,
    nReducedDims = effective_reduced_dims,
    k0s = as.integer(k0s),
    subsample = TRUE,
    sequential = TRUE,
    clusterFunction = main_cluster_function,
    alphas = alphas,
    betas = beta,
    minSizes = as.integer(min_cluster_size),
    consensusProportion = consensus_proportion,
    verbose = verbose,
    parameterWarnings = TRUE,
    subsampleArgs = list(
      clusterFunction = subsample_cluster_function,
      resamp.num = as.integer(resamples),
      samp.p = subsample_fraction
    ),
    seqArgs = list(
      remain.n = remaining_items_to_stop
    ),
    ncores = as.integer(n_cores),
    random.seed = as.integer(random_seed)
  )

  cluster_raw <- clusterExperiment::primaryCluster(rsec)
  if (is.factor(cluster_raw)) {
    cluster_raw <- as.character(cluster_raw)
  }
  cluster <- suppressWarnings(as.integer(cluster_raw))
  if (length(cluster) != nrow(scaled_matrix) || anyNA(cluster)) {
    stop("RSEC returned invalid primary cluster assignments.")
  }
  names(cluster) <- rownames(scaled_matrix)

  rsec_order <- clusterExperiment::orderSamples(rsec)
  if (is.numeric(rsec_order) &&
      length(rsec_order) == nrow(scaled_matrix) &&
      setequal(rsec_order, seq_len(nrow(scaled_matrix)))) {
    item_order <- rownames(scaled_matrix)[rsec_order]
  } else if (is.character(rsec_order) &&
             setequal(rsec_order, rownames(scaled_matrix))) {
    item_order <- rsec_order
  } else {
    item_order <- names(cluster)[order(cluster, names(cluster))]
  }

  co_clustering <- complete_co_clustering_matrix(
    clusterExperiment::coClustering(rsec),
    rownames(scaled_matrix)
  )
  assignments <- data.frame(
    item = names(cluster),
    cluster = as.integer(cluster),
    item_order = match(names(cluster), item_order),
    clustering_method = "rsec",
    stringsAsFactors = FALSE
  )

  list(
    assignments = assignments,
    scaled_matrix = scaled_matrix,
    item_order = item_order,
    excluded_items = prepared$excluded_items,
    model = rsec,
    co_clustering = co_clustering,
    silhouette = calculate_pattern_silhouette(
      scaled_matrix,
      cluster,
      unassigned_values = -1L
    ),
    parameters = list(
      effective_reduced_dims = effective_reduced_dims,
      k0s = as.integer(k0s),
      resamples = as.integer(resamples),
      subsample_fraction = subsample_fraction,
      min_cluster_size = as.integer(min_cluster_size),
      remaining_items_to_stop = remaining_items_to_stop
    )
  )
}
