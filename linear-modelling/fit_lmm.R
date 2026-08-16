#' Fit a Linear Mixed-Effects Model (LMM) or Linear Model (LM)
#'
#' Fits a linear mixed-effects model using `lmer()` from `lmerTest` package,
#' with user-specified fixed effects, random effects, and contrasts.
#' Retries non-converged mixed models with a pre-specified optimizer sequence.
#' For singular random-intercept models, removes only the near-zero variance
#' component, refits, and repeats. Falls back to `lm()` only when no supported
#' random effects remain.
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
#' @details Before fitting, all model-required columns are checked for missing
#' values, response and contrast estimability is validated, constant adjustment
#' terms are removed, random-intercept replication is diagnosed, and the fixed
#' design matrix is checked for full rank. Invalid fits stop with data
#' diagnostics; conditions that remain estimable are retained and flagged.
#'
#' Simple continuous numeric variables supplied through `fixed_effects` are
#' centered and scaled internally before fitting. This is an invertible change
#' of units that improves numerical conditioning; the applied centers and
#' scales are returned in `lmm_diagnostics$scaled_fixed_effects`. Contrast terms,
#' binary numeric fixed effects, and transformed terms are not scaled.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{full_model}{The fitted model object returned by `lmer()` or `lm()`.}
#'   \item{model_terms_test}{ANOVA results for fixed effects terms.}
#'   \item{model_df}{A tidy `data.frame` containing the ANOVA table and contrasts.}
#'   \item{lmm_diagnostics}{Pre-fit data summaries and warning flags, optimizer
#'   attempts, fixed- and random-effect removal decisions, retained model terms,
#'   final variance components, and continuous-covariate scaling metadata.}
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
# Resolve dependencies beside this script so sourcing does not depend on the
# caller's working directory. R_UTILS_DIR supports execution contexts that do
# not expose the current source-file path.
fit_lmm_source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(condition) NULL
)
linear_modelling_dir <- if (
  length(fit_lmm_source_file) == 1L &&
    !is.null(fit_lmm_source_file) &&
    !is.na(fit_lmm_source_file) &&
    nzchar(fit_lmm_source_file)
) {
  dirname(normalizePath(fit_lmm_source_file, mustWork = TRUE))
} else {
  file.path(Sys.getenv("R_UTILS_DIR", unset = "."), "linear-modelling")
}

# Load the contrast and model-diagnostic helpers into the same environment as
# fit_lmm(), including when this file is sourced into a test environment.
source(file.path(linear_modelling_dir, "lm.R"), local = environment())
source(
  file.path(linear_modelling_dir, "fit_lmm_helpers.R"),
  local = environment()
)

# Avoid exposing path-resolution implementation details to downstream callers.
rm(fit_lmm_source_file, linear_modelling_dir)

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
    equivalence_deltas = numeric(),
    compute_contrasts = TRUE,
    check_assumptions = TRUE
) {
  method <- match.arg(method)
  requested_fixed_effects <- unique(fixed_effects)
  requested_random_effects <- unique(random_effects)
  prepared_inputs <- prepare_lmm_inputs(
    dataset = dataset,
    response = response,
    contrast = contrast,
    fixed_effects = requested_fixed_effects,
    random_effects = requested_random_effects
  )
  dataset <- prepared_inputs$dataset
  fixed_effects <- prepared_inputs$fixed_effects
  current_random_effects <- prepared_inputs$random_effects
  data_diagnostics <- prepared_inputs$diagnostics

  # Rank-transform the response variable if requested
  if (method == "rank_transform") {
    if (verbose) message("Rank-transforming response variable: ", response)
    dataset[[response]] <- rank(dataset[[response]])
  }

  # Determine if any structurally estimable random intercepts remain.
  use_lmer <- length(current_random_effects) > 0
  full_fm <- build_lmm_formula(
    response = response,
    contrast = contrast,
    fixed_effects = fixed_effects,
    random_effects =
    if (use_lmer) current_random_effects else NULL
  )
  if (verbose) {
    message("Fitting model with formula: ", deparse(full_fm))
  }

  lmm_diagnostics <- list(
    data_diagnostics = data_diagnostics,
    attempts = data.frame(),
    simplification_steps = data.frame(),
    requested_fixed_effects = requested_fixed_effects,
    retained_fixed_effects = fixed_effects,
    removed_fixed_effects = data_diagnostics$fixed_effect_removals,
    requested_random_effects = requested_random_effects,
    retained_random_effects = current_random_effects,
    removed_random_effects =
      data_diagnostics$structural_random_effect_removals$random_effect,
    final_variance_components = NULL,
    final_singular = FALSE,
    final_formula = paste(deparse(full_fm), collapse = " "),
    scaled_fixed_effects = data_diagnostics$scaled_fixed_effects
  )

  # Try a mixed model first if random effects are specified. A converged but
  # singular fit is simplified one near-zero random intercept at a time. Each
  # reduced model starts the convergence checks again.
  if (use_lmer) {
    repeat {
      retry_fit <- fit_lmm_with_retries(full_fm, dataset)
      retry_fit$attempts$formula <- paste(deparse(full_fm), collapse = " ")
      retry_fit$attempts$random_effects <- paste(
        current_random_effects,
        collapse = ","
      )
      lmm_diagnostics$attempts <- dplyr::bind_rows(
        lmm_diagnostics$attempts,
        retry_fit$attempts
      )
      fit_converged_before_removal <- !is.null(retry_fit$model)

      if (fit_converged_before_removal) {
        full_model <- retry_fit$model
        if (!lme4::isSingular(full_model)) {
          break
        }

        variance_components <- get_random_intercept_variances(full_model)
        zero_variance_components <- variance_components %>%
          dplyr::filter(.data$sd_relative_to_residual < 1e-4)

        if (nrow(zero_variance_components) == 0) {
          stop_lmm_with_diagnostics(
            paste0(
              "Mixed model is singular, but no near-zero random-intercept ",
              "variance component could be identified safely."
            ),
            data_diagnostics
          )
        }

        component_to_remove <- zero_variance_components$grp[
          which.min(zero_variance_components$sd_relative_to_residual)
        ]
        component_variance <- zero_variance_components$vcov[
          zero_variance_components$grp == component_to_remove
        ][1]
        component_relative_sd <-
          zero_variance_components$sd_relative_to_residual[
            zero_variance_components$grp == component_to_remove
          ][1]
        simplification_evidence <- "converged singular fit"
      } else {
        if (length(retry_fit$zero_random_effect_candidates) == 0) {
          stop_lmm_with_diagnostics(
            paste0(
              "Mixed model failed to converge after all optimizer attempts. ",
              retry_fit$failure_details
            ),
            data_diagnostics
          )
        }

        candidate_components <- retry_fit$candidate_variance_components %>%
          dplyr::group_by(.data$grp) %>%
          dplyr::summarise(
            optimizer_support = dplyr::n_distinct(.data$attempt),
            max_relative_sd = max(.data$sd_relative_to_residual),
            median_variance = stats::median(.data$vcov),
            .groups = "drop"
          ) %>%
          dplyr::arrange(
            dplyr::desc(.data$optimizer_support),
            .data$max_relative_sd,
            .data$grp
          )
        component_to_remove <- candidate_components$grp[1]
        component_variance <- candidate_components$median_variance[1]
        component_relative_sd <- candidate_components$max_relative_sd[1]
        simplification_evidence <- paste0(
          candidate_components$optimizer_support[1],
          " non-converged singular optimizer fit(s) identified this near-zero component"
        )
      }

      lmm_diagnostics$simplification_steps <- dplyr::bind_rows(
        lmm_diagnostics$simplification_steps,
        data.frame(
          removed_random_effect = component_to_remove,
          variance = component_variance,
          sd_relative_to_residual = component_relative_sd,
          formula_before_removal = paste(deparse(full_fm), collapse = " "),
          fit_converged_before_removal = fit_converged_before_removal,
          evidence = simplification_evidence,
          stringsAsFactors = FALSE
        )
      )
      lmm_diagnostics$removed_random_effects <- c(
        lmm_diagnostics$removed_random_effects,
        component_to_remove
      )

      if (verbose) {
        message(
          "Singular fit: removing near-zero random intercept '",
          component_to_remove,
          "' and refitting."
        )
      }

      current_random_effects <- setdiff(
        current_random_effects,
        component_to_remove
      )

      if (length(current_random_effects) == 0) {
        full_fm <- build_lmm_formula(
          response = response,
          contrast = contrast,
          fixed_effects = fixed_effects,
          random_effects = NULL
        )
        full_model <- tryCatch(
          stats::lm(formula = full_fm, data = dataset),
          error = function(e) {
            stop_lmm_with_diagnostics(
              paste0("Linear-model fallback failed: ", conditionMessage(e)),
              data_diagnostics
            )
          }
        )
        use_lmer <- FALSE
        break
      }

      full_fm <- build_lmm_formula(
        response = response,
        contrast = contrast,
        fixed_effects = fixed_effects,
        random_effects = current_random_effects
      )
    }
  } else {
    # No random effects specified - use lm directly
    if (verbose) {
      message("No random effects specified. Fitting linear model (lm).")
    }
    full_model <- tryCatch(
      stats::lm(formula = full_fm, data = dataset),
      error = function(e) {
        stop_lmm_with_diagnostics(
          paste0("Linear model failed: ", conditionMessage(e)),
          data_diagnostics
        )
      }
    )
  }

  lmm_diagnostics$retained_random_effects <- current_random_effects
  lmm_diagnostics$final_formula <- paste(deparse(full_fm), collapse = " ")
  if (use_lmer) {
    lmm_diagnostics$final_variance_components <-
      get_random_intercept_variances(full_model)
    lmm_diagnostics$final_singular <- lme4::isSingular(full_model)
  } else {
    lmm_diagnostics$retained_random_effects <- character()
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
    ci_mx <- lme4::confint.merMod(full_model, method = "Wald", parm = fixed_effect_names)
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

  if (compute_contrasts) {
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
  } else {
    model_df <- bind_rows(model_df, model_terms_test_df)
  }

  # Save model objects to list
  models <- list()
  models$full_model <- full_model
  models$model_terms_test <- model_terms_test
  models$model_df <- model_df
  models$model_type <- if(use_lmer) "lmer" else "lm"  # Track which type was used
  models$method <- method
  models$lmm_diagnostics <- lmm_diagnostics

  if (check_assumptions) {
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
  } else {
    models$assumption_tests <- NULL
  }

  if (save_models) {
    saveRDS(models, file.path(save_model_dir, "lmm.modelobj.rds"))
  }

  return(models)
}
