#' Fit a generalized additive mixed model for a continuous predictor
#'
#' This is the nonlinear counterpart of `fit_lmm()`. It fits one smooth for a
#' numeric predictor, optional ordinary fixed effects, and optional random
#' intercepts represented by mgcv random-effect smooths. The returned design
#' information is deliberately retained so generic Wald contrasts and smooth
#' predictions can be calculated without knowing anything about the analysis
#' that called this function.
#'
#' @param dataset Data frame containing the response and model variables.
#' @param response Name of the numeric response column.
#' @param predictor Name of the numeric continuous predictor column.
#' @param random_effects Character vector of grouping variables for random
#'   intercepts, for example `"individual_id"` for longitudinal samples.
#' @param fixed_effects Character vector of additional additive fixed effects.
#' @param k Basis dimension for the predictor smooth. It is reduced when a
#'   feature has fewer unique predictor values.
#' @param smooth_basis Basis passed to `mgcv::s()`; cubic regression splines are
#'   the default, matching the basis used by tradeSeq.
#' @param family An mgcv family object. Gaussian is appropriate for the current
#'   continuous, transformed feature values.
#' @param method Smoothing-parameter estimation method passed to `mgcv::gam()`.
#' @param scale_response If TRUE, standardize the response before fitting and
#'   retain the centre/scale required to return predictions on its original
#'   scale.
#' @param verbose Print the fitted formula.
#' @param save_models Save the returned object as `gamm.modelobj.rds`.
#' @param save_model_dir Destination used when `save_models` is TRUE.
#' @param ... Additional arguments passed to `mgcv::gam()`.
#'
#' @return A `feature_gamm` object containing `full_model`, its formula, the
#'   complete model data, variable names, predictor range, and response scaling.
fit_gamm <- function(
    dataset,
    response,
    predictor,
    random_effects = NULL,
    fixed_effects = NULL,
    k = 6L,
    smooth_basis = "cr",
    family = stats::gaussian(),
    method = "REML",
    scale_response = FALSE,
    verbose = FALSE,
    save_models = FALSE,
    save_model_dir = ".",
    ...
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required by fit_gamm(). Install it before running ",
      "the GAMM analysis."
    )
  }

  if (!is.data.frame(dataset)) stop("dataset must be a data frame.")
  if (length(response) != 1L || length(predictor) != 1L) {
    stop("response and predictor must each name exactly one column.")
  }

  fixed_effects <- unique(fixed_effects %||% character())
  random_effects <- unique(random_effects %||% character())
  required_cols <- unique(c(response, predictor, fixed_effects, random_effects))
  missing_cols <- setdiff(required_cols, colnames(dataset))
  if (length(missing_cols) > 0L) {
    stop("Missing model columns: ", paste(missing_cols, collapse = ", "))
  }
  if (anyNA(dataset[, required_cols, drop = FALSE])) {
    stop("Model columns must not contain missing values; subset complete rows first.")
  }
  if (!is.numeric(dataset[[response]])) stop("response must be numeric.")
  if (!is.numeric(dataset[[predictor]])) stop("predictor must be numeric.")
  if (any(!is.finite(dataset[[response]]))) stop("response must be finite.")
  if (any(!is.finite(dataset[[predictor]]))) stop("predictor must be finite.")

  n_unique_predictor <- length(unique(dataset[[predictor]]))
  if (n_unique_predictor < 4L) {
    stop("At least four unique predictor values are required for a smooth.")
  }
  response_sd <- stats::sd(dataset[[response]])
  if (!is.finite(response_sd) || response_sd == 0) {
    stop("Response has zero or undefined variance.")
  }

  for (variable in c(fixed_effects, random_effects)) {
    if (length(unique(dataset[[variable]])) < 2L) {
      stop("Model variable has fewer than two observed values: ", variable)
    }
  }

  model_data <- dataset
  response_center <- if (scale_response) mean(model_data[[response]]) else 0
  response_scale <- if (scale_response) response_sd else 1
  model_data$.gamm_response <-
    (model_data[[response]] - response_center) / response_scale
  model_data$.gamm_predictor <- model_data[[predictor]]

  requested_k <- as.integer(k)
  if (!is.finite(requested_k) || requested_k < 3L) {
    stop("k must be an integer of at least 3.")
  }
  effective_k <- min(requested_k, n_unique_predictor)

  quote_name <- function(x) paste0("`", gsub("`", "\\`", x), "`")
  model_terms <- c(
    paste0(
      "s(.gamm_predictor, k = ", effective_k,
      ", bs = '", smooth_basis, "')"
    ),
    vapply(fixed_effects, quote_name, character(1)),
    vapply(
      random_effects,
      function(x) paste0("s(", quote_name(x), ", bs = 're')"),
      character(1)
    )
  )
  model_formula <- stats::as.formula(paste(
    ".gamm_response ~", paste(model_terms, collapse = " + ")
  ))

  if (verbose) message("Fitting GAMM with formula: ", deparse(model_formula))
  full_model <- mgcv::gam(
    formula = model_formula,
    data = model_data,
    family = family,
    method = method,
    ...
  )

  output <- list(
    full_model = full_model,
    model_formula = model_formula,
    model_data = model_data,
    response = response,
    predictor = predictor,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    requested_k = requested_k,
    effective_k = effective_k,
    smooth_basis = smooth_basis,
    family = family$family,
    method = method,
    response_center = response_center,
    response_scale = response_scale,
    predictor_range = range(model_data$.gamm_predictor),
    n_observations = nrow(model_data),
    n_individuals = if (length(random_effects) > 0L) {
      length(unique(model_data[[random_effects[[1L]]]]))
    } else {
      NA_integer_
    },
    model_type = if (length(random_effects) > 0L) "gamm" else "gam"
  )
  class(output) <- c("feature_gamm", "list")

  if (save_models) {
    dir.create(save_model_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(output, file.path(save_model_dir, "gamm.modelobj.rds"))
  }

  output
}

# Internal null-coalescing helper kept local to this standalone utility.
`%||%` <- function(x, y) if (is.null(x)) y else x

