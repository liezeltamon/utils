#' Test normality of numeric features
#'
#' Assesses the normality of every numeric column in a data frame or matrix
#' using visual diagnostics (QQ plots, density plots) and a statistical test.
#' Non-numeric columns are dropped automatically. Plots are written to PDF files in
#' \code{output_dir}, and a summary data frame is returned.
#'
#' @param x A data frame or matrix. Non-numeric columns are dropped
#'   automatically.
#' @param output_dir Character. Path to the directory where PDF plots are
#'   saved. Created recursively if it does not exist. Default
#'   \code{"normality_output"}.
#' @param alpha Numeric in (0, 1). Significance threshold applied to the
#'   multiple-testing-adjusted p-values when determining \code{is_normal}.
#'   Default \code{0.05}.
#' @param correct Character. Multiple-testing correction method passed to
#'   \code{\link[stats]{p.adjust}}. One of \code{"BH"} (Benjamini-Hochberg,
#'   default), \code{"bonferroni"}, \code{"holm"}, \code{"none"}, or any other
#'   method accepted by \code{p.adjust}.
#' @param test Character. Statistical test to use: \code{"shapiro"} (Shapiro-
#'   Wilk, requires n \eqn{\le} 5,000), \code{"ks"}, or \code{"auto"}
#'   (default), which selects Shapiro-Wilk when n \eqn{\le} 5,000 and otherwise
#'   selects \code{"ks"}. Note: the \code{"ks"} branch is not yet implemented
#'   and currently stops with an error.
#' @param verbose Logical. Whether to emit progress messages. Default
#'   \code{TRUE}.
#'
#' @details
#' When testing many features simultaneously, the probability of spurious
#' rejections grows with the number of tests. \code{test_normality} therefore
#' adjusts raw p-values via \code{p.adjust} before comparing against
#' \code{alpha}. The default Benjamini-Hochberg correction controls the false
#' discovery rate, which is appropriate for exploratory feature screening.
#' Use \code{correct = "bonferroni"} for stricter family-wise error rate
#' control, or \code{correct = "none"} to work with unadjusted p-values.
#'
#' For features with fewer than 3 non-missing observations, no statistical test
#' is performed and \code{NA} values are returned for the test result fields.
#'
#' @return A data frame with one row per numeric feature and the following
#'   columns:
#'   \describe{
#'     \item{feature}{Column name.}
#'     \item{n}{Number of non-missing observations.}
#'     \item{n_missing}{Number of missing values.}
#'     \item{test_used}{Name of the statistical test applied.}
#'     \item{statistic}{Test statistic (W for Shapiro-Wilk).}
#'     \item{p_value}{Raw p-value.}
#'     \item{p_adj}{Adjusted p-value after multiple-testing correction.}
#'     \item{is_normal}{Numeric indicator; \code{1} when \code{p_adj >= alpha}
#'       and \code{0} otherwise.}
#'   }
#'
#' @importFrom ggplot2 ggplot aes stat_qq stat_qq_line geom_density geom_line
#'   labs theme_bw theme element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices pdf dev.off
#' @importFrom stats shapiro.test ks.test dnorm p.adjust p.adjust.methods
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' demo <- data.frame(
#'   id = paste0("s", 1:200),
#'   normal_A = rnorm(200, mean = 5, sd = 1),
#'   normal_B = rnorm(200, mean = 0, sd = 3),
#'   skewed_C = rexp(200, rate = 0.5),
#'   bimodal_D = c(rnorm(100, -3), rnorm(100, 3))
#' )
#'
#' results <- test_normality(demo, output_dir = "normality_output")
#' print(results)
#' }
#'
#' @export
test_normality <- function(
    x,
    output_dir = "normality_output",
    alpha = 0.05,
    correct = "BH",
    test = c("auto", "shapiro", "ks"),
    verbose = TRUE
) {
  test <- match.arg(test)
  correct <- match.arg(correct, p.adjust.methods)

  # Check inputs
  stopifnot(is.numeric(alpha) && alpha > 0 && alpha < 1)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }
  if (!is.data.frame(x)) x <- as.data.frame(x)

  # Identify numeric columns and drop non-numeric ones
  numeric_cols <- names(x)[sapply(x, is.numeric)]
  dropped <- setdiff(names(x), numeric_cols)
  if (length(dropped) > 0) {
    data_num <- x[, numeric_cols, drop = FALSE]
    if (verbose) message(sprintf(
      "Dropping %d non-numeric column(s): %s", length(dropped), paste(dropped, collapse = ", ")
    ))
  } else {
    data_num <- x
  }

  # Run tests for each numeric column
  if (verbose) message("Running statistical tests...")
  test_results <- lapply(data_num, .run_test, test = test)

  # Collect results into a data frame
  results_df <- data.frame(
    feature = colnames(data_num),
    n_nonmissing = sapply(data_num, function(v) sum(!is.na(v))),
    n_missing = sapply(data_num, function(v) sum(is.na(v))),
    test_used = sapply(test_results, `[[`, "method"),
    statistic = sapply(test_results, `[[`, "statistic"),
    p_value = sapply(test_results, `[[`, "p_value"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  results_df$p_adj <- p.adjust(results_df$p_value, method = correct)
  results_df$is_normal <- as.numeric(!is.na(results_df$p_adj) & (results_df$p_adj >= alpha)) 

  if (verbose) message("Save normality results to CSV...")
  write.csv(
    results_df, file.path(output_dir, "normality_test_results.csv"), row.names = FALSE
  )

  if (verbose) {
    n_normal <- sum(results_df$is_normal, na.rm = TRUE)
    message(sprintf(
      "Correction: %s | alpha: %.3f | Normal features: %d / %d", correct, alpha, n_normal, nrow(results_df)
    ))
  }

  if (verbose) message("Generating QQ plots...")
  .save_plot_grid(
    Map(.make_qq_plot, colnames(data_num), data_num),
    file.path(output_dir, "plots_qq.pdf"),
    verbose = verbose
  )

  if (verbose) message("Generating density plots...")
  .save_plot_grid(
    Map(.make_density_plot, colnames(data_num), data_num),
    file.path(output_dir, "plots_density.pdf"),
    verbose = verbose
  )

  if (verbose) {
    message("Done.")
    print(results_df)
  }

  return(results_df)
}

#' @noRd
.choose_test <- function(col_vals, test) {
  n_complete <- sum(!is.na(col_vals))
  if (test == "auto") return(if (n_complete <= 5000) "shapiro" else "ks")
  if (test == "shapiro" && n_complete > 5000)
    warning("Shapiro-Wilk is appropriate for n <= 5000. Consider using ks.")
  test
}

#' @noRd
.run_test <- function(col_vals, test) {
  col_vals <- col_vals[!is.na(col_vals)]
  if (length(col_vals) < 3)
    return(list(statistic = NA_real_, p_value = NA_real_, method = NA_character_))

  chosen <- .choose_test(col_vals, test)
  if (chosen == "shapiro") {
    res <- tryCatch(
      shapiro.test(col_vals),
      error = function(e) list(statistic = NA_real_, p.value = NA_real_)
    )
    list(statistic = unname(res$statistic), p_value = res$p.value, method = "Shapiro-Wilk")
  } else if (chosen == "ks") {
    # res <- tryCatch(
    #   ks.test(col_vals, "pnorm"), #, mean = mean(col_vals), sd = sd(col_vals)),
    #   error = function(e) list(statistic = NA_real_, p.value = NA_real_)
    # )
    stop(
      "Need to implement Lilliefors-corrected KS test for normality which does not require pre-specifying mean and sd. Consider using the 'nortest' package or similar. Check if Lilliefors test applicable for large samples."
    )
    #list(statistic = unname(res$statistic), p_value = res$p.value, method = "Kolmogorov-Smirnov")
  } else {
    stop(".run_test(): Invalid test argument.")
  }
}

#' @noRd
.make_qq_plot <- function(feature_name, vals) {
  df_plot <- data.frame(v = vals[!is.na(vals)])
  ggplot2::ggplot(df_plot, ggplot2::aes(sample = v)) +
    ggplot2::stat_qq(colour = "#2166ac", alpha = 0.6, size = 1) +
    ggplot2::stat_qq_line(colour = "#d6604d", linewidth = 0.8) +
    ggplot2::labs(
      title = feature_name,
      subtitle = sprintf("n = %d", nrow(df_plot)),
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' @noRd
.make_density_plot <- function(feature_name, vals) {
  vals_clean <- vals[!is.na(vals)]
  df_plot <- data.frame(v = vals_clean)
  x_seq <- seq(min(vals_clean), max(vals_clean), length.out = 300)
  normal_df <- data.frame(
    x = x_seq,
    y = dnorm(x_seq, mean = mean(vals_clean), sd = sd(vals_clean))
  )
  ggplot2::ggplot(df_plot, ggplot2::aes(x = v)) +
    ggplot2::geom_density(fill = "#2166ac", colour = "#2166ac", alpha = 0.3) +
    ggplot2::geom_line(
      data = normal_df,
      ggplot2::aes(x = x, y = y),
      colour = "#d6604d", linewidth = 0.9, linetype = "dashed"
    ) +
    ggplot2::labs(
      title = feature_name,
      subtitle = sprintf("n = %d  |  blue = empirical  |  red dashed = normal fit", nrow(df_plot)),
      x = feature_name,
      y = "Density"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' @noRd
.save_plot_grid <- function(plot_list, filepath, ncol = 2, verbose = TRUE) {
  nrow_est <- ceiling(length(plot_list) / ncol)
  grDevices::pdf(filepath, width = 7 * ncol, height = min(4 * nrow_est, 50))
  gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
  grDevices::dev.off()
  if (verbose) message("  Saved: ", filepath)
}
