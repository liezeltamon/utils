#' Fit a Linear Mixed-Effects Model (LMM)
#'
#' Fits a linear mixed-effects model using `lmer()` from `lmerTest` package,
#' with user-specified fixed effects, random effects, and contrasts.
#'
#' @param dataset A data frame containing the contrast factor, random, and fixed effects for the model
#' @param response A `string` specifying the name of the response variable
#' @param contrast A character vector specifying the variables (main effects or interactions)
#'   for which contrasts should be computed (e.g., `c("group", "expanded", "group:expanded")`).
#' @param fixed_effects A character vector of additional fixed effects
#' @param random_effects A character vector of grouping variables for random intercepts
#' @param verbose Logical; if `TRUE`, prints model-fitting progress and messages
#' @param save_models Logical; if `TRUE`, saves the fitted model object(s) to disk
#' @param save_model_dir Character; directory path where models should be saved (default: current directory)
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{full_model}{The fitted model object returned by `lmer()`.}
#'   \item{model_terms_test}{ANOVA results for fixed effects terms, computed using `anova()` from `lmerTest`.}
#'   \item{model_df}{A tidy `data.frame` containing the Type III ANOVA table.}
#' }
#'
#' @examples
#' \dontrun{
#' output <- fit_lmm(
#'   dataset = cell_data,
#'   response = "score1",
#'   contrast = c("group", "expanded", "group:expanded"),
#'   random_effects = "individual",
#'   fixed_effects = "pool_id",
#'   save_models = FALSE
#' )
#' }
#'
#' @import tidyverse
#' @import lme4
#' @import lmerTest
source(file.path("utils", "lm.R"))

# library(lme4)
# library(lmerTest)
# source(file.path("utils", "lm.R"))
# 
# dataset = cell_data
# response = "score1"
# contrast = c("group", "expanded", "group:expanded")
# fixed_effects = c("pool_id")
# random_effects = "individual"
# 
# output <- fit_lmm(
#   dataset = cell_data,
#   response = "score1",
#   contrast = c("group", "expanded", "group:expanded"),
#   random_effects = "individual",
#   fixed_effects = c("pool_id"),
#   save_models = FALSE
# )

fit_lmm <- function (
    dataset,
    response,
    contrast,
    random_effects = NULL,
    fixed_effects = NULL,
    verbose = FALSE,
    save_models = FALSE,
    save_model_dir = "."
    #n_cores = 1,
    #p.adjust_method = "BH"
) {
  
  model_rhs <- create_model_formula(
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    verbose = verbose
  )
  full_fm <- as.formula(paste0(
    c(paste0(response, " ~ ", paste(contrast, collapse = " + "), " + "), model_rhs), collapse = ""
  ))
  full_model <- lmerTest::lmer(
    formula = full_fm,
    data = dataset
  )
  # Equivalent to anova(full_model, type = "III", ddf = "Satterthwaite")
  model_terms_test <- lmerTest:::anova.lmerModLmerTest(
    full_model,
    type = "III", ddf = "Satterthwaite"
  )
  
  fixed_effect_names <- names(fixef(full_model))
  ci_mx <- confint.merMod(full_model, method = "Wald", parm = fixed_effect_names)
  colnames(ci_mx) <- c("coef.95pct.ci.lower", "coef.95pct.ci.upper")
  # Collect other info including odds ratio
  model_info <- ci_mx 
  coeffs <- summary(full_model)$coefficients
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
  
  # Pairwise contrasts with emmeans()
  
  contrasts <- lapply(contrast, function(x) {
    pairwise_contrasts(
      full_model,
      contrast = x,
      adjust_contrast = "tukey",
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
  
  if (save_models) {
    saveRDS(models, file.path(save_model_dir, "lmm.modelobj.rds"))
  }
  return(models)
  
}
