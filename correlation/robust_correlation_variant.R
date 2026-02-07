robust_correlation_variant <- function(
    feature_matrix,
    features = colnames(feature_matrix),
    subsample_depth,
    n_repeats = 10000,
    method = "pearson",
    n_cores = parallel::detectCores() - 1,
    random_seed = 1,
    use_crossprod = FALSE
) {
  stopifnot(is.matrix(feature_matrix) || is.data.frame(feature_matrix))
  stopifnot(subsample_depth > 1)
  
  # Validate crossprod + method combination
  if (use_crossprod && method != "pearson") {
    stop("use_crossprod = TRUE only works with method = 'pearson'")
  }
  
  set.seed(random_seed)
  
  mat <- as.matrix(feature_matrix[, features])
  n_samples <- nrow(mat)
  
  # Setup parallel backend
  n_cores <- max(1, n_cores)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  })
  
  `%dopar%` <- foreach::`%dopar%`
  
  # Parallelize over subsamples (not feature pairs)
  corr_sum <- foreach::foreach(
    rep = seq_len(n_repeats),
    .combine = "+",
    .packages = c("stats")
  ) %dopar% {
    # Subsample rows
    sel <- sample(n_samples, subsample_depth, replace = FALSE)
    subsample_mat <- mat[sel, , drop = FALSE]
    
    # Compute full correlation matrix for this subsample
    suppressWarnings({
      if (use_crossprod) {
        mat_scaled <- scale(subsample_mat, center = TRUE, scale = TRUE)
        crossprod(mat_scaled) / (subsample_depth - 1)
      } else {
        cor(subsample_mat, method = method)
      }
    })
  }
  
  # Average across all subsamples
  corr_mat <- corr_sum / n_repeats
  dimnames(corr_mat) <- list(features, features)
  
  corr_mat
}

# Example usage:

# # Simple test data
# set.seed(123)

# n_samples <- 50
# n_features <- 10

# # Create feature matrix
# feature_matrix <- matrix(rnorm(n_samples * n_features), 
#                         nrow = n_samples, 
#                         ncol = n_features)

# # Add correlation: features 1-10 are correlated
# common_signal <- rnorm(n_samples)
# feature_matrix[, 1:10] <- feature_matrix[, 1:10] + 
#   outer(common_signal, rep(1, 10)) * 2

# colnames(feature_matrix) <- paste0("Feature_", 1:n_features)

# # Test it
# cor_result <- robust_correlation_variant(
#   feature_matrix,
#   subsample_depth = 30,
#   n_repeats = 10,
#   n_cores = 1,
#   use_crossprod = FALSE
# )

# # Check results
# cat("Correlation matrix dimensions:", dim(cor_result), "\n")
# cat("\nFeatures 1-5 correlations (should be high):\n")
# print(round(cor_result[1:5, 1:5], 3))


# cor_result_old <- vdjremix::robust_correlation(
#   feature_matrix,
#   subsample_depth = 30,
#   n_repeats = 10,
#   n_cores = 1
# )

# # 1. Overall similarity (should be > 0.99 if same method)
# cor(as.vector(cor_result), as.vector(cor_result_old))

# # 2. Biggest difference
# max(abs(cor_result - cor_result_old), na.rm = TRUE)

# # 3. Average difference
# mean(abs(cor_result - cor_result_old), na.rm = TRUE)
