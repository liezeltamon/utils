#' Compute pairwise contrasts with emmeans and return tidied results
#'
#' This function calculates estimated marginal means (EMMs) and pairwise contrasts 
#' for a specified grouping variable from a fitted model using the `emmeans` package.
#' It returns a tidy data frame with odds ratios, confidence intervals, and p-values, 
#' and optionally saves the raw emmeans output.
#'
#'#' @import tidyverse
#'#' @import emmeans
#'
#' @param model A fitted model object compatible with `emmeans`.
#' @param contrast A string specifying the grouping variable for which to compute EMMs and contrasts.
#' @param adjust_contrast A character string specifying the method to adjust p-values for multiple comparisons (e.g., "bonferroni", "holm", "none").
#' @param save_objects Logical; if TRUE, saves the emmeans output as an RDS file in `save_dir`.
#' @param save_dir Directory path to save the emmeans output; defaults to current working directory.
#' @param ... Additional arguments passed to `emmeans` or `contrast` functions.
#'
#' @return A data frame containing pairwise comparison results including odds ratios (OR), 
#' confidence intervals, and adjusted p-values for each contrast.
#'
#' @details
#' The function first computes EMMs for \code{contrast} without adjustment, then calculates pairwise contrasts with the specified p-value adjustment.
#' Confidence intervals and significance tests are computed on the response scale (e.g., odds ratios for logistic models).
#' The results are combined, tidied, and returned as a long-format data frame.
#' 

pairwise_contrasts <- function (
    model,
    contrast,
    adjust_contrast,
    save_objects = TRUE,
    save_dir = ".",
    ...
) {

  emm_emmeans <- emmeans(
     model, as.formula(paste("~ ", contrast)), adjust = "none"
  )
  emm_contrasts <- contrast(emm_emmeans, "revpairwise", adjust = adjust_contrast) %>% 
    # Add confidence intervals and test for significance of estimates
    summary(level = 0.95, infer = c(TRUE, TRUE), type = "response")
  
  emm_emmeans <- emm_emmeans %>% 
    summary(level = 0.95, infer = c(TRUE, TRUE), type = "response")
  emm <- list(emm_emmeans, emm_contrasts)
  
  if (save_objects) {
    saveRDS(emm, file.path(save_dir, "emmeans_output.rds"))
  }
  
  # Store information into a data frame
  lst <- lapply(emm, function(x) {
     tibble(x) %>% 
      pivot_longer(cols = where(is.numeric), names_to = "statistic") %>% 
      unite(contrast, -statistic, -value, sep = "..", remove = TRUE) %>% 
      rename(contrast = 1)
  })
  df <- do.call("rbind", lst) %>% as.data.frame() %>% 
    mutate(
      statistic = recode(
        statistic,
        odds.ratio = "coef",
        p.value = "coef.pvalue",
        asymp.UCL = "coef.95pct.ci.upper",
        asymp.LCL = "coef.95pct.ci.lower"
      )
    ) %>% 
    separate(
      contrast,
      c("var", "contrast"),
      sep = "\\s*/\\s*|\\s*-\\s*",
      fill = "right",
      remove = FALSE
    ) %>% 
    mutate(
      var = str_replace_all(var, "\\s", ".."),
      contrast = if_else(is.na(contrast), NA_character_, str_replace_all(contrast, "\\s", "..")),
    )
  
  return(df)
  
}

#' Create model formulas - code from original MASC implementation

create_model_formula <- function (fixed_effects, random_effects, verbose = TRUE) {
  
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(
      c(paste0(fixed_effects, collapse = " + "),
        paste0("(1|", random_effects, ")", collapse = " + ")),
      collapse = " + "
    )
    if (verbose == TRUE) {
      message(paste("Using null model:", "response ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "response ~", model_rhs))
      #stop("No random effects specified")
      model_rhs <- paste0(fixed_effects, collapse = " + ")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "response ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "response ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  return(model_rhs)
}

