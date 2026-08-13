# Helpers for fit_lmm().
#
# The public fitting function remains in utils/fit_lmm.R. This file contains
# formula-aware validation, diagnostic collection, optimizer retries, and
# random-intercept variance helpers. The implementation is generic: variables
# are identified from the function arguments and formula, never by project-
# specific names.

lmm_empty_data_frame <- function(columns) {
  values <- lapply(columns, function(type) {
    switch(
      type,
      character = character(),
      integer = integer(),
      numeric = numeric(),
      logical = logical(),
      stop("Unsupported empty-column type: ", type, call. = FALSE)
    )
  })
  names(values) <- names(columns)
  as.data.frame(values, stringsAsFactors = FALSE)
}

lmm_term_variables <- function(term) {
  all.vars(stats::as.formula(paste("~", term)))
}

build_lmm_formula <- function(
    response,
    contrast,
    fixed_effects = NULL,
    random_effects = NULL
) {
  fixed_terms <- unique(c(contrast, fixed_effects))
  fixed_terms <- fixed_terms[!is.na(fixed_terms) & nzchar(fixed_terms)]
  random_terms <- unique(random_effects)
  random_terms <- random_terms[!is.na(random_terms) & nzchar(random_terms)]

  random_formula_terms <- if (length(random_terms) == 0) {
    character()
  } else {
    paste0("(1|", random_terms, ")")
  }
  rhs_terms <- c(fixed_terms, random_formula_terms)
  if (length(rhs_terms) == 0) {
    rhs_terms <- "1"
  }
  stats::as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
}

lmm_unique_count <- function(x) {
  length(unique(x))
}

lmm_value_label <- function(x) {
  if (length(x) == 0) {
    return("")
  }
  paste(as.character(unique(x)), collapse = ",")
}

lmm_is_categorical <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

collect_lmm_response_diagnostics <- function(dataset, response) {
  values <- dataset[[response]]
  counts <- table(values, useNA = "no")
  modal_count <- if (length(counts) == 0) 0L else max(counts)
  n <- length(values)

  data.frame(
    response = response,
    n = n,
    unique_values = lmm_unique_count(values),
    unique_fraction = if (n == 0) NA_real_ else lmm_unique_count(values) / n,
    modal_count = as.integer(modal_count),
    modal_fraction = if (n == 0) NA_real_ else modal_count / n,
    stringsAsFactors = FALSE
  )
}

collect_lmm_variable_diagnostics <- function(dataset, variables, role) {
  variables <- unique(variables)
  if (length(variables) == 0) {
    return(lmm_empty_data_frame(c(
      variable = "character",
      role = "character",
      class = "character",
      n_missing = "integer",
      n_unique_observed = "integer",
      observed_values = "character"
    )))
  }

  do.call(rbind, lapply(variables, function(variable) {
    values <- dataset[[variable]]
    observed <- values[!is.na(values)]
    data.frame(
      variable = variable,
      role = role,
      class = paste(class(values), collapse = ","),
      n_missing = sum(is.na(values)),
      n_unique_observed = lmm_unique_count(observed),
      observed_values = if (lmm_unique_count(observed) <= 10) {
        lmm_value_label(observed)
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  }))
}

collect_lmm_random_effect_diagnostics <- function(dataset, random_effects) {
  random_effects <- unique(random_effects)
  if (length(random_effects) == 0) {
    return(lmm_empty_data_frame(c(
      random_effect = "character",
      n_levels = "integer",
      singleton_levels = "integer",
      repeated_levels = "integer",
      singleton_fraction = "numeric",
      min_observations = "integer",
      median_observations = "numeric",
      max_observations = "integer",
      structurally_estimable = "logical",
      structural_reason = "character"
    )))
  }

  do.call(rbind, lapply(random_effects, function(random_effect) {
    counts <- table(dataset[[random_effect]], useNA = "no")
    n_levels <- length(counts)
    singleton_levels <- sum(counts == 1)
    repeated_levels <- sum(counts >= 2)
    reason <- if (n_levels < 2) {
      "fewer than two observed levels"
    } else if (repeated_levels == 0) {
      "every observed level occurs once"
    } else {
      ""
    }

    data.frame(
      random_effect = random_effect,
      n_levels = n_levels,
      singleton_levels = singleton_levels,
      repeated_levels = repeated_levels,
      singleton_fraction = if (n_levels == 0) NA_real_ else {
        singleton_levels / n_levels
      },
      min_observations = if (n_levels == 0) NA_integer_ else min(counts),
      median_observations = if (n_levels == 0) NA_real_ else {
        stats::median(as.numeric(counts))
      },
      max_observations = if (n_levels == 0) NA_integer_ else max(counts),
      structurally_estimable = !nzchar(reason),
      structural_reason = reason,
      stringsAsFactors = FALSE
    )
  }))
}

collect_lmm_overlap_diagnostics <- function(
    dataset,
    contrast_variables,
    random_effects
) {
  contrast_variables <- unique(contrast_variables)
  random_effects <- unique(random_effects)
  rows <- list()
  row_index <- 0L

  for (contrast_variable in contrast_variables) {
    contrast_values <- dataset[[contrast_variable]]
    if (!lmm_is_categorical(contrast_values)) {
      next
    }
    observed_contrast_levels <- unique(as.character(contrast_values))
    if (length(observed_contrast_levels) < 2) {
      next
    }

    for (random_effect in random_effects) {
      random_values <- as.character(dataset[[random_effect]])
      random_levels <- unique(random_values)
      groups_per_random_level <- vapply(
        random_levels,
        function(level) {
          length(unique(as.character(contrast_values[random_values == level])))
        },
        integer(1)
      )
      mixed_levels <- sum(groups_per_random_level >= 2)
      contrast_pairs <- utils::combn(observed_contrast_levels, 2, simplify = FALSE)

      for (contrast_pair in contrast_pairs) {
        shared_levels <- sum(vapply(
          random_levels,
          function(level) {
            values_in_level <- unique(as.character(
              contrast_values[random_values == level]
            ))
            all(contrast_pair %in% values_in_level)
          },
          logical(1)
        ))
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          contrast_variable = contrast_variable,
          random_effect = random_effect,
          contrast_pair = paste(contrast_pair, collapse = "_vs_"),
          random_levels = length(random_levels),
          mixed_random_levels = mixed_levels,
          mixed_random_level_fraction = mixed_levels / length(random_levels),
          shared_random_levels = shared_levels,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(lmm_empty_data_frame(c(
      contrast_variable = "character",
      random_effect = "character",
      contrast_pair = "character",
      random_levels = "integer",
      mixed_random_levels = "integer",
      mixed_random_level_fraction = "numeric",
      shared_random_levels = "integer"
    )))
  }
  do.call(rbind, rows)
}

format_lmm_data_diagnostics <- function(diagnostics) {
  parts <- character()

  if (!is.null(diagnostics$response) && nrow(diagnostics$response) > 0) {
    response <- diagnostics$response[1, ]
    parts <- c(parts, paste0(
      "n=", response$n,
      ", response_unique=", response$unique_values,
      ", response_unique_fraction=", format(response$unique_fraction, digits = 4),
      ", modal_count=", response$modal_count,
      ", modal_fraction=", format(response$modal_fraction, digits = 4)
    ))
  }

  if (!is.null(diagnostics$variable_summary) &&
      nrow(diagnostics$variable_summary) > 0) {
    missing_rows <- diagnostics$variable_summary[
      diagnostics$variable_summary$n_missing > 0,
      ,
      drop = FALSE
    ]
    if (nrow(missing_rows) > 0) {
      parts <- c(parts, paste0(
        "missing={",
        paste0(missing_rows$variable, ":", missing_rows$n_missing,
               collapse = ","),
        "}"
      ))
    }
  }

  if (!is.null(diagnostics$random_effects) &&
      nrow(diagnostics$random_effects) > 0) {
    random_parts <- apply(diagnostics$random_effects, 1, function(row) {
      paste0(
        row[["random_effect"]],
        "(levels=", row[["n_levels"]],
        ",singletons=", row[["singleton_levels"]],
        ",repeated=", row[["repeated_levels"]], ")"
      )
    })
    parts <- c(parts, paste0(
      "random_effect_replication={",
      paste(random_parts, collapse = ","),
      "}"
    ))
  }

  if (!is.null(diagnostics$fixed_design) &&
      nrow(diagnostics$fixed_design) > 0) {
    fixed_design <- diagnostics$fixed_design[1, ]
    parts <- c(parts, paste0(
      "fixed_design_rank=", fixed_design$rank,
      "/", fixed_design$columns
    ))
  }

  if (!is.null(diagnostics$overlap) && nrow(diagnostics$overlap) > 0) {
    overlap <- diagnostics$overlap
    minimum_shared <- min(overlap$shared_random_levels)
    parts <- c(parts, paste0(
      "minimum_shared_random_levels=", minimum_shared
    ))
  }

  if (length(parts) == 0) {
    return("Data diagnostics unavailable.")
  }
  paste0("Data diagnostics: ", paste(parts, collapse = "; "), ".")
}

stop_lmm_with_diagnostics <- function(message, diagnostics = NULL) {
  if (!is.null(diagnostics)) {
    message <- paste(message, format_lmm_data_diagnostics(diagnostics))
  }
  stop(message, call. = FALSE)
}

prepare_lmm_inputs <- function(
    dataset,
    response,
    contrast,
    fixed_effects = NULL,
    random_effects = NULL,
    modal_warning_threshold = 0.90,
    few_random_levels_threshold = 4L
) {
  requested_fixed_effects <- unique(fixed_effects)
  requested_random_effects <- unique(random_effects)
  unsupported_random_effects <- setdiff(requested_random_effects, names(dataset))
  if (length(unsupported_random_effects) > 0) {
    stop(
      "random_effects must name dataset columns used as separate random ",
      "intercepts. Unsupported entries: ",
      paste(unsupported_random_effects, collapse = ", "),
      call. = FALSE
    )
  }
  requested_formula <- build_lmm_formula(
    response = response,
    contrast = contrast,
    fixed_effects = requested_fixed_effects,
    random_effects = requested_random_effects
  )
  required_variables <- all.vars(requested_formula)
  missing_columns <- setdiff(required_variables, names(dataset))
  if (length(missing_columns) > 0) {
    stop(
      "Model-required columns are absent from dataset: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  contrast_variables <- unique(unlist(
    lapply(contrast, lmm_term_variables),
    use.names = FALSE
  ))
  fixed_variables <- unique(unlist(
    lapply(requested_fixed_effects, lmm_term_variables),
    use.names = FALSE
  ))
  random_variables <- unique(unlist(
    lapply(requested_random_effects, lmm_term_variables),
    use.names = FALSE
  ))

  variable_summary <- dplyr::bind_rows(
    collect_lmm_variable_diagnostics(dataset, response, "response"),
    collect_lmm_variable_diagnostics(dataset, contrast_variables, "contrast"),
    collect_lmm_variable_diagnostics(dataset, fixed_variables, "fixed_effect"),
    collect_lmm_variable_diagnostics(dataset, random_variables, "random_effect")
  )
  diagnostics <- list(
    requested_formula = paste(deparse(requested_formula), collapse = " "),
    required_variables = required_variables,
    variable_summary = variable_summary,
    response = collect_lmm_response_diagnostics(dataset, response),
    fixed_effect_removals = lmm_empty_data_frame(c(
      fixed_effect = "character",
      constant_variables = "character",
      observed_values = "character",
      reason = "character"
    )),
    random_effects = collect_lmm_random_effect_diagnostics(
      dataset,
      requested_random_effects
    ),
    structural_random_effect_removals = lmm_empty_data_frame(c(
      random_effect = "character",
      reason = "character"
    )),
    fixed_design = lmm_empty_data_frame(c(
      rank = "integer",
      columns = "integer",
      rank_deficient = "logical",
      non_estimable_columns = "character"
    )),
    overlap = collect_lmm_overlap_diagnostics(
      dataset,
      contrast_variables,
      requested_random_effects
    ),
    warning_flags = lmm_empty_data_frame(c(
      code = "character",
      variable = "character",
      message = "character"
    )),
    scaled_fixed_effects = lmm_empty_data_frame(c(
      variable = "character",
      center = "numeric",
      scale = "numeric"
    ))
  )

  missing_counts <- vapply(
    required_variables,
    function(variable) sum(is.na(dataset[[variable]])),
    integer(1)
  )
  missing_required <- missing_counts[missing_counts > 0]
  if (length(missing_required) > 0) {
    stop_lmm_with_diagnostics(
      paste0(
        "Model-required variables contain missing values: ",
        paste0(
          names(missing_required),
          "=",
          missing_required,
          collapse = ", "
        ),
        ". Upstream code must resolve missingness explicitly; fit_lmm() will not omit rows."
      ),
      diagnostics
    )
  }

  if (nrow(dataset) == 0) {
    stop_lmm_with_diagnostics("Model dataset has zero rows.", diagnostics)
  }
  if (diagnostics$response$unique_values < 2) {
    stop_lmm_with_diagnostics(
      "Response has fewer than two unique values and cannot be modelled.",
      diagnostics
    )
  }
  if (!is.numeric(dataset[[response]])) {
    stop_lmm_with_diagnostics(
      paste0(
        "Response '", response,
        "' is not numeric; fit_lmm() fits Gaussian linear models."
      ),
      diagnostics
    )
  }

  for (contrast_variable in contrast_variables) {
    values <- dataset[[contrast_variable]]
    if (lmm_unique_count(values) < 2) {
      stop_lmm_with_diagnostics(
        paste0(
          "Contrast variable '", contrast_variable,
          "' has fewer than two observed values."
        ),
        diagnostics
      )
    }
    if (is.factor(values)) {
      absent_declared_levels <- setdiff(levels(values), unique(as.character(values)))
      if (length(absent_declared_levels) > 0) {
        stop_lmm_with_diagnostics(
          paste0(
            "Contrast variable '", contrast_variable,
            "' has declared levels with no observations: ",
            paste(absent_declared_levels, collapse = ", "),
            "."
          ),
          diagnostics
        )
      }
    }
  }

  # Drop unused factor levels only after verifying that no declared contrast
  # level is absent. Unused adjustment/grouping levels contain no information.
  dataset <- droplevels(dataset)

  retained_fixed_effects <- requested_fixed_effects
  for (fixed_effect in requested_fixed_effects) {
    dependencies <- lmm_term_variables(fixed_effect)
    constant_dependencies <- dependencies[vapply(
      dependencies,
      function(variable) lmm_unique_count(dataset[[variable]]) < 2,
      logical(1)
    )]
    if (length(constant_dependencies) == 0) {
      next
    }

    retained_fixed_effects <- setdiff(retained_fixed_effects, fixed_effect)
    diagnostics$fixed_effect_removals <- rbind(
      diagnostics$fixed_effect_removals,
      data.frame(
        fixed_effect = fixed_effect,
        constant_variables = paste(constant_dependencies, collapse = ","),
        observed_values = paste(vapply(
          constant_dependencies,
          function(variable) lmm_value_label(dataset[[variable]]),
          character(1)
        ), collapse = ";"),
        reason = "constant adjustment variable",
        stringsAsFactors = FALSE
      )
    )
  }

  retained_random_effects <- requested_random_effects
  if (nrow(diagnostics$random_effects) > 0) {
    structural_removals <- diagnostics$random_effects[
      !diagnostics$random_effects$structurally_estimable,
      c("random_effect", "structural_reason"),
      drop = FALSE
    ]
    if (nrow(structural_removals) > 0) {
      names(structural_removals)[2] <- "reason"
      diagnostics$structural_random_effect_removals <- structural_removals
      retained_random_effects <- setdiff(
        retained_random_effects,
        structural_removals$random_effect
      )
    }
  }

  response_diagnostics <- diagnostics$response
  if (response_diagnostics$unique_values == 2) {
    diagnostics$warning_flags <- rbind(
      diagnostics$warning_flags,
      data.frame(
        code = "response_two_unique_values",
        variable = response,
        message = "Response has only two unique values; Gaussian-model suitability should be reviewed.",
        stringsAsFactors = FALSE
      )
    )
  }
  if (is.finite(response_diagnostics$modal_fraction) &&
      response_diagnostics$modal_fraction >= modal_warning_threshold) {
    diagnostics$warning_flags <- rbind(
      diagnostics$warning_flags,
      data.frame(
        code = "response_dominant_value",
        variable = response,
        message = paste0(
          "One response value accounts for ",
          format(100 * response_diagnostics$modal_fraction, digits = 4),
          "% of observations."
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  if (nrow(diagnostics$random_effects) > 0) {
    partly_singleton <- diagnostics$random_effects[
      diagnostics$random_effects$structurally_estimable &
        diagnostics$random_effects$singleton_levels > 0,
      ,
      drop = FALSE
    ]
    for (row_index in seq_len(nrow(partly_singleton))) {
      diagnostics$warning_flags <- rbind(
        diagnostics$warning_flags,
        data.frame(
          code = "random_effect_singleton_levels",
          variable = partly_singleton$random_effect[row_index],
          message = paste0(
            partly_singleton$singleton_levels[row_index], " of ",
            partly_singleton$n_levels[row_index],
            " random-effect levels are represented by one observation."
          ),
          stringsAsFactors = FALSE
        )
      )
    }

    few_levels <- diagnostics$random_effects[
      diagnostics$random_effects$structurally_estimable &
        diagnostics$random_effects$n_levels <= few_random_levels_threshold,
      ,
      drop = FALSE
    ]
    for (row_index in seq_len(nrow(few_levels))) {
      diagnostics$warning_flags <- rbind(
        diagnostics$warning_flags,
        data.frame(
          code = "random_effect_few_levels",
          variable = few_levels$random_effect[row_index],
          message = paste0(
            "Random effect has only ", few_levels$n_levels[row_index],
            " observed levels; its variance estimate may be imprecise."
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  if (nrow(diagnostics$overlap) > 0) {
    no_overlap <- diagnostics$overlap[
      diagnostics$overlap$shared_random_levels == 0,
      ,
      drop = FALSE
    ]
    if (nrow(no_overlap) > 0) {
      diagnostics$warning_flags <- rbind(
        diagnostics$warning_flags,
        data.frame(
          code = "no_within_random_level_contrast_overlap",
          variable = paste(unique(no_overlap$random_effect), collapse = ","),
          message = paste0(
            nrow(no_overlap),
            " contrast/random-effect combinations have no shared random-effect levels."
          ),
          stringsAsFactors = FALSE
        )
      )
    }
    limited_overlap <- diagnostics$overlap[
      diagnostics$overlap$shared_random_levels == 1,
      ,
      drop = FALSE
    ]
    if (nrow(limited_overlap) > 0) {
      diagnostics$warning_flags <- rbind(
        diagnostics$warning_flags,
        data.frame(
          code = "limited_within_random_level_contrast_overlap",
          variable = paste(unique(limited_overlap$random_effect), collapse = ","),
          message = paste0(
            nrow(limited_overlap),
            " contrast/random-effect combinations share only one random-effect level."
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Scale only simple, continuous numeric adjustment terms. Scaling the
  # underlying variable of a transformed term such as log(x) would change that
  # model rather than merely changing coefficient units.
  simple_fixed_variables <- intersect(retained_fixed_effects, names(dataset))
  protected_variables <- unique(c(contrast_variables, retained_random_effects))
  simple_fixed_variables <- setdiff(simple_fixed_variables, protected_variables)
  for (variable in simple_fixed_variables) {
    values <- dataset[[variable]]
    if (!is.numeric(values) || lmm_unique_count(values) <= 2) {
      next
    }
    center <- mean(values)
    scale <- stats::sd(values)
    if (!is.finite(scale) || scale == 0) {
      next
    }
    dataset[[variable]] <- (values - center) / scale
    diagnostics$scaled_fixed_effects <- rbind(
      diagnostics$scaled_fixed_effects,
      data.frame(
        variable = variable,
        center = center,
        scale = scale,
        stringsAsFactors = FALSE
      )
    )
  }

  fixed_formula <- build_lmm_formula(
    response = response,
    contrast = contrast,
    fixed_effects = retained_fixed_effects,
    random_effects = NULL
  )
  fixed_design <- tryCatch(
    stats::model.matrix(fixed_formula, data = dataset),
    error = function(e) e
  )
  if (inherits(fixed_design, "error")) {
    stop_lmm_with_diagnostics(
      paste0(
        "Could not construct the fixed-effect design matrix: ",
        conditionMessage(fixed_design)
      ),
      diagnostics
    )
  }
  fixed_qr <- qr(fixed_design)
  non_estimable_columns <- if (fixed_qr$rank < ncol(fixed_design)) {
    colnames(fixed_design)[fixed_qr$pivot[seq.int(fixed_qr$rank + 1, ncol(fixed_design))]]
  } else {
    character()
  }
  diagnostics$fixed_design <- data.frame(
    rank = fixed_qr$rank,
    columns = ncol(fixed_design),
    rank_deficient = fixed_qr$rank < ncol(fixed_design),
    non_estimable_columns = paste(non_estimable_columns, collapse = ","),
    stringsAsFactors = FALSE
  )
  if (diagnostics$fixed_design$rank_deficient) {
    stop_lmm_with_diagnostics(
      paste0(
        "Fixed-effect design is rank deficient (rank ", fixed_qr$rank,
        " of ", ncol(fixed_design), "). Non-estimable columns: ",
        paste(non_estimable_columns, collapse = ", "),
        ". Varying fixed effects will not be dropped automatically."
      ),
      diagnostics
    )
  }

  list(
    dataset = dataset,
    fixed_effects = retained_fixed_effects,
    random_effects = retained_random_effects,
    diagnostics = diagnostics
  )
}

# Return the convergence information that lme4 stores on a fitted model. A
# mixed model is accepted only when both the optimizer return code and lme4's
# derivative-based checks indicate convergence.
get_lmm_convergence_status <- function(model) {
  opt_code <- model@optinfo$conv$opt
  convergence_messages <- model@optinfo$conv$lme4$messages
  optimizer_message <- model@optinfo$message

  optimizer_converged <-
    is.null(opt_code) || length(opt_code) == 0 || all(opt_code == 0)
  checks_converged <-
    is.null(convergence_messages) || length(convergence_messages) == 0

  gradient <- model@optinfo$derivs$gradient
  max_abs_gradient <- if (is.null(gradient) || length(gradient) == 0) {
    NA_real_
  } else {
    max(abs(gradient), na.rm = TRUE)
  }

  list(
    converged = optimizer_converged && checks_converged,
    optimizer_code = opt_code,
    convergence_messages = convergence_messages,
    optimizer_message = optimizer_message,
    max_abs_gradient = max_abs_gradient
  )
}

# Extract random-intercept variances and express their standard deviations
# relative to the residual SD. fit_lmm() supports separate random intercepts.
get_random_intercept_variances <- function(model) {
  variance_components <- as.data.frame(lme4::VarCorr(model))
  residual_sd <- stats::sigma(model)

  variance_components <- variance_components[
    is.na(variance_components$var2) &
      variance_components$grp != "Residual",
    c("grp", "vcov", "sdcor"),
    drop = FALSE
  ]
  variance_components$sd_relative_to_residual <-
    variance_components$sdcor / residual_sd
  variance_components
}

# Fit the same mixed model with a fixed retry sequence. Only numerical controls
# change; formula, data, likelihood, and the default REML setting remain fixed.
fit_lmm_with_retries <- function(formula, dataset) {
  collapse_diagnostic <- function(value) {
    if (is.null(value) || length(value) == 0) "" else {
      paste(value, collapse = " | ")
    }
  }

  attempt_specs <- list(
    list(name = "default", control = NULL),
    list(
      name = "nloptwrap_maxeval_200000",
      control = lme4::lmerControl(
        optimizer = "nloptwrap",
        optCtrl = list(maxeval = 200000)
      )
    ),
    list(
      name = "bobyqa_maxfun_200000",
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 200000)
      )
    )
  )

  attempt_diagnostics <- list()
  attempt_models <- list()
  attempt_variance_components <- list()

  for (attempt_index in seq_along(attempt_specs)) {
    attempt_spec <- attempt_specs[[attempt_index]]
    captured_warnings <- character()
    captured_messages <- character()
    fit <- tryCatch(
      withCallingHandlers(
        {
          if (is.null(attempt_spec$control)) {
            lmerTest::lmer(formula = formula, data = dataset)
          } else {
            lmerTest::lmer(
              formula = formula,
              data = dataset,
              control = attempt_spec$control
            )
          }
        },
        warning = function(w) {
          captured_warnings <<- c(captured_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          captured_messages <<- c(captured_messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      attempt_models[attempt_index] <- list(NULL)
      attempt_variance_components[attempt_index] <- list(NULL)
      attempt_diagnostics[[attempt_index]] <- data.frame(
        attempt = attempt_index,
        requested_optimizer = attempt_spec$name,
        fitted_optimizer = NA_character_,
        converged = FALSE,
        optimizer_code = NA_character_,
        convergence_messages = conditionMessage(fit),
        optimizer_message = NA_character_,
        warnings = paste(unique(captured_warnings), collapse = " | "),
        messages = paste(unique(captured_messages), collapse = " | "),
        max_abs_gradient = NA_real_,
        log_likelihood = NA_real_,
        singular = NA,
        near_zero_random_effects = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    convergence <- get_lmm_convergence_status(fit)
    singular <- lme4::isSingular(fit)
    variance_components <- get_random_intercept_variances(fit)
    near_zero_random_effects <- variance_components$grp[
      variance_components$sd_relative_to_residual < 1e-4
    ]
    variance_components$attempt <- attempt_index
    attempt_models[[attempt_index]] <- fit
    attempt_variance_components[[attempt_index]] <- variance_components
    attempt_diagnostics[[attempt_index]] <- data.frame(
      attempt = attempt_index,
      requested_optimizer = attempt_spec$name,
      fitted_optimizer = paste(fit@optinfo$optimizer, collapse = ","),
      converged = convergence$converged,
      optimizer_code = collapse_diagnostic(convergence$optimizer_code),
      convergence_messages = collapse_diagnostic(
        convergence$convergence_messages
      ),
      optimizer_message = collapse_diagnostic(convergence$optimizer_message),
      warnings = paste(unique(captured_warnings), collapse = " | "),
      messages = paste(unique(captured_messages), collapse = " | "),
      max_abs_gradient = convergence$max_abs_gradient,
      log_likelihood = as.numeric(stats::logLik(fit)),
      singular = singular,
      near_zero_random_effects = paste(near_zero_random_effects, collapse = ","),
      stringsAsFactors = FALSE
    )

    if (convergence$converged) {
      return(list(
        model = fit,
        attempts = do.call(rbind, attempt_diagnostics),
        zero_random_effect_candidates = character(),
        candidate_variance_components = data.frame(),
        variance_components = dplyr::bind_rows(attempt_variance_components),
        failure_details = ""
      ))
    }
  }

  attempt_diagnostics <- do.call(rbind, attempt_diagnostics)
  failure_details <- paste0(
    attempt_diagnostics$requested_optimizer,
    ": ",
    ifelse(
      nzchar(attempt_diagnostics$convergence_messages),
      attempt_diagnostics$convergence_messages,
      ifelse(
        nzchar(attempt_diagnostics$optimizer_message),
        attempt_diagnostics$optimizer_message,
        "did not satisfy convergence checks"
      )
    ),
    collapse = "; "
  )
  fitted_attempt_indices <- which(!vapply(attempt_models, is.null, logical(1)))
  singular_attempt_indices <- fitted_attempt_indices[
    attempt_diagnostics$singular[fitted_attempt_indices] %in% TRUE
  ]
  candidate_variance_components <- dplyr::bind_rows(
    attempt_variance_components[singular_attempt_indices]
  )
  if (nrow(candidate_variance_components) > 0) {
    candidate_variance_components <- candidate_variance_components[
      candidate_variance_components$sd_relative_to_residual < 1e-4,
      ,
      drop = FALSE
    ]
  } else {
    candidate_variance_components <- data.frame(
      grp = character(),
      vcov = numeric(),
      sdcor = numeric(),
      sd_relative_to_residual = numeric(),
      attempt = integer(),
      stringsAsFactors = FALSE
    )
  }

  list(
    model = NULL,
    attempts = attempt_diagnostics,
    zero_random_effect_candidates = unique(candidate_variance_components$grp),
    candidate_variance_components = candidate_variance_components,
    variance_components = dplyr::bind_rows(attempt_variance_components),
    failure_details = failure_details
  )
}
