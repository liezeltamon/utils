#' Fit a Linear Mixed-Effects Model (LMM) or Linear Model (LM)
#'
#' Fits a linear mixed-effects model using `lmer()` from `lmerTest` package,
#' with user-specified fixed effects, random effects, and contrasts.
#' Automatically falls back to `lm()` if no random effects specified or if
#' singular fit is detected.
#'
#' @param dataset A data frame containing the contrast factor, random, and fixed effects for the model
#' @param response A `string` specifying the name of the response variable
#' @param contrast A character vector specifying the variables (main effects or interactions)
#'   for which contrasts should be computed (e.g., `c("group", "expanded", "group:expanded")`). The first variable in the vector will be used for Levene's test of homogeneity of variance and make sure it is a main effect rather than an interaction for this reason.
#' @param fixed_effects A character vector of additional fixed effects
#' @param random_effects A character vector of grouping variables for random intercepts
#' @param verbose Logical; if `TRUE`, prints model-fitting progress and messages
#' @param save_models Logical; if `TRUE`, saves the fitted model object(s) to disk
#' @param save_model_dir Character; directory path where models should be saved (default: current directory)
#' @param method Character; method for handling non-parametric data. Options are "lmer" (default) for linear mixed model, "rank_transform" for rank-transforming the response variable and fitting a linear model, and "art" for using the aligned rank transform method (to be implemented).
#' @param additional_contrasts Optional named list of user-defined emmeans contrasts to compute in addition to pairwise contrasts.
#' @param equivalence_deltas Optional numeric vector of equivalence margins for emmeans TOST tests.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{full_model}{The fitted model object returned by `lmer()` or `lm()`.}
#'   \item{model_terms_test}{ANOVA results for fixed effects terms.}
#'   \item{model_df}{A tidy `data.frame` containing the ANOVA table and contrasts.}
#' }
#'
#' @examples
#' \dontrun{
#' # With random effects (mixed model)
#' output <- fit_lmm(
#'   dataset = cell_data,
#'   response = "score1",
#'   contrast = c("group", "expanded", "group:expanded"),
#'   random_effects = "individual",
#'   fixed_effects = "pool_id",
#'   save_models = FALSE
#' )
#'
#' # Without random effects (linear model)
#' output <- fit_lmm(
#'   dataset = cell_data,
#'   response = "score1",
#'   contrast = "group",
#'   random_effects = NULL,
#'   save_models = FALSE
#' )
#' }
#'
#' @import tidyverse
#' @import lme4
#' @import lmerTest
source(file.path("utils", "lm.R"))

fit_lmm <- function (
    dataset,
    response,
    contrast,
    random_effects = NULL,
    fixed_effects = NULL,
    verbose = FALSE,
    save_models = FALSE,
    save_model_dir = ".",
    method = c("lmer", "rank_transform", "art"), # To implement "art" for non-parametric alternative involving interactions (see https://depts.washington.edu/acelab/proj/art), check how to interpret main effects from ART models
    additional_contrasts = NULL,
    equivalence_deltas = numeric()
) {

  stopifnot(
    "Response should not have missing values." = !anyNA(dataset[[response]])
  )

  method <- match.arg(method)

  # Rank-transform the response variable if requested
  if (method == "rank_transform") {
    if (verbose) message("Rank-transforming response variable: ", response)
    dataset[[response]] <- rank(dataset[[response]])
  }

  # Determine if we should try mixed model
  use_lmer <- !is.null(random_effects) && length(random_effects) > 0

  model_rhs <- create_model_formula(
    fixed_effects = fixed_effects,
    random_effects = if (use_lmer) random_effects else NULL,
    verbose = verbose
  )

  full_fm <- as.formula(paste0(
    c(paste0(response, " ~ ", paste(contrast, collapse = " + "), " + "), model_rhs),
    collapse = ""
  ))
  if (verbose) {
    message("Fitting model with formula: ", deparse(full_fm))
  }

  # Try mixed model first if random effects specified
  if (use_lmer) {
    full_model <- lmerTest::lmer(formula = full_fm, data = dataset)

    # Check for singular fit
    if (lme4::isSingular(full_model)) {
      message(
        "Singular fit detected. Refitting as linear model (lm) without random effects."
      )
      # Refit as simple linear model
      full_fm_lm <- as.formula(paste0(
        response, " ~ ", paste(c(contrast, fixed_effects), collapse = " + ")
      ))
      full_model <- lm(formula = full_fm_lm, data = dataset)
      use_lmer <- FALSE  # Update flag for appropriate downstream processing
    }
  } else {
    # No random effects specified - use lm directly
    if (verbose) {
      message("No random effects specified. Fitting linear model (lm).")
    }
    full_model <- lm(formula = full_fm, data = dataset)
  }

  # Extract results based on model type
  if (use_lmer) {
    # Mixed model results
    model_terms_test <- lmerTest:::anova.lmerModLmerTest(
      full_model,
      type = "III",
      ddf = "Satterthwaite"
    )
    fixed_effect_names <- names(lme4::fixef(full_model))
    ci_mx <- confint.merMod(full_model, method = "Wald", parm = fixed_effect_names)
    coeffs <- summary(full_model)$coefficients
  } else {
    # Linear model results
    model_terms_test <- anova(full_model)
    ci_mx <- confint(full_model)
    coeffs <- summary(full_model)$coefficients
  }

  colnames(ci_mx) <- c("coef.95pct.ci.lower", "coef.95pct.ci.upper")

  # Collect model info including coefficients and CIs
  model_info <- ci_mx
  colnames(coeffs)[1] <- "coef"
  model_info <- cbind(coeffs[rownames(model_info), ], model_info)

  # Tidy model_info
  model_df <- model_info %>%
    as.data.frame() %>%
    tibble::rownames_to_column("var") %>%
    pivot_longer(-var, names_to = "statistic") %>%
    unite("statistic", var, statistic, sep = ".", remove = FALSE)

  # Tidy model_terms_test
  model_terms_test_df <- model_terms_test %>%
    as.data.frame() %>%
    tibble::rownames_to_column("var") %>%
    pivot_longer(-var, names_to = "statistic") %>%
    unite("statistic", var, statistic, sep = ".", remove = FALSE)

  if (!is.null(additional_contrasts) && length(contrast) > 1) {
    if (!is.list(additional_contrasts) || is.null(names(additional_contrasts))) {
      stop(
        "When multiple contrast terms are supplied, additional_contrasts ",
        "must be named by contrast term."
      )
    }
    extra_contrast_names <- setdiff(names(additional_contrasts), contrast)
    if (length(extra_contrast_names) > 0) {
      stop(
        "additional_contrasts contains names that are not in contrast: ",
        paste(extra_contrast_names, collapse = ", ")
      )
    }
  }

  # Pairwise contrasts with emmeans()
  contrasts <- lapply(contrast, function(x) {
    additional_contrasts_for_x <- NULL
    if (!is.null(additional_contrasts)) {
      if (!is.list(additional_contrasts)) {
        stop("additional_contrasts must be a list.")
      }
      if (length(contrast) > 1) {
        if (x %in% names(additional_contrasts)) {
          additional_contrasts_for_x <- additional_contrasts[[x]]
        }
      } else if (
        !is.null(names(additional_contrasts)) &&
        x %in% names(additional_contrasts) &&
        is.list(additional_contrasts[[x]]) &&
        !is.numeric(additional_contrasts[[x]])
      ) {
        additional_contrasts_for_x <- additional_contrasts[[x]]
      } else {
        additional_contrasts_for_x <- additional_contrasts
      }
    }

    pairwise_contrasts(
      full_model,
      contrast = x,
      adjust_contrast = "none", #"tukey",
      additional_contrasts = additional_contrasts_for_x,
      equivalence_deltas = equivalence_deltas,
      save_objects = TRUE,
      save_dir = save_model_dir
    )
  })
  contrasts_df <- do.call("rbind", contrasts)

  # Combine
  model_df <- bind_rows(model_df, model_terms_test_df, contrasts_df) %>%
    select(colnames(contrasts_df))

  # Save model objects to list
  models <- list()
  models$full_model <- full_model
  models$model_terms_test <- model_terms_test
  models$model_df <- model_df
  models$model_type <- if(use_lmer) "lmer" else "lm"  # Track which type was used
  models$method <- method

  # Assumption checks
  # performance::check_model(full_model) # Visual checks of assumptions
  shapiro <- shapiro.test(residuals(full_model))
  levene_group <- contrast[1]
  if (grepl(":", levene_group)) {
    warning(
      "Levene's test skipped: '", levene_group, "' is an interaction term. Reorder contrast so a main effect is first."
    )
    levene <- NULL
  } else {
    model_df <- model.frame(full_model)
    levene <- car::leveneTest(residuals(full_model) ~ model_df[[levene_group]])
  }
  models$assumption_tests <- list(
    shapiro = shapiro,
    levene = levene,
    assumptions_violated = shapiro$p.value < 0.05 | (!is.null(levene) && levene[1, "Pr(>F)"] < 0.05)
  )

  if (save_models) {
    saveRDS(models, file.path(save_model_dir, "lmm.modelobj.rds"))
  }

  return(models)
}
