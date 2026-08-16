# Functions to assess batch effects, integrating with different methods and regressing out different potential batch effect sources, and assessing integration

library(cowplot)
library(ggplot2)
library(harmony)
library(lisi)
library(pheatmap)
library(RColorBrewer)
library(scater)
library(Seurat)
library(BiocParallel)

# Diagnose batch effects from covariates

# Plot contribution of each covariate in variance of each gene and data as a whole
# Concerning if technical covariates have higher contribution than variables likely driving biology
#library(cowplot)
#library(scater)
plot_variableContribution <- function(sce,
                                      covars,
                                      subset_row, nCPU,
                                      dimred, n_dimred, # For scater::getExplanatoryPCs()
                                      plot_title = NULL
                                      ){
  
  # Remove covars with <= 2 levels
  remove_covars <- names(which(
    sapply(covars, USE.NAMES = TRUE, simplify = TRUE, function(var_nme) length(unique(sce[[var_nme]]))) <= 2
  ))
  if (length(remove_covars > 0)) {
    covars <- setdiff(covars, remove_covars)
    warning("plot_variableContribution(): Removing covars with <= 2 levels")
  }
  
  # Per-gene variance explained by a variable
  covars_varexp <- scater::getVarianceExplained(sce, variables = covars, subset_row = subset_row, BPPARAM = MulticoreParam(nCPU))
  p_varexp <- scater::plotExplanatoryVariables(covars_varexp, nvars_to_plot = Inf) + ggtitle(plot_title)
  exp_pcs <- scater::getExplanatoryPCs(sce, dimred = dimred, n_dimred = n_dimred, variables = covars)
  p_exp_pcs <- plotExplanatoryPCs(exp_pcs) # From source code -> scale_y_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100))
  p <- plot_grid(p_varexp, p_exp_pcs, ncol = 1)
  print(p)
  return(p)

}

# Integrate using different methods

#library(harmony)
#library(plyr)
#library(Seurat)
integrate_data <- function(obj, # Seurat object
                           dimred_name,
                           dims_use,
                           get_umap = TRUE,
                           get_tsne = FALSE,
                           seed_val = 290,
                           # Harmony parameters
                           batch_vars, # List of vectors of variables to regress together. If NULL, skip harmony::RunHamorny()
                           max_iter = 50, # Higher than default, increase if not converging based on Harmony plots generated
                           # Seurat::IntegrateData() parameters
                           split_vars = NULL, # If NULL, skip Seurat::IntegrateLayers()
                           seurat_assay = "RNA",
                           norm_method = "LogNormalize",
                           scaledata_covars = NULL) {

  DefaultAssay(obj) <- seurat_assay

  # Remove batch_vars and split vars with 1 level
  vars_unique <- unique(unlist(batch_vars, split_vars))
  is_vars_1lvl <- unlist(lapply(obj[[vars_unique]], function(col) {
    length(unique(col)) == 1
  }))
  vars_1lvl <- names(which(is_vars_1lvl))

  if (length(vars_1lvl) > 0) {

    batch_vars <- unique(lapply(batch_vars, function(vars_set){
      vars_trimmed <- setdiff(vars_set, vars_1lvl)
      if (length(vars_trimmed) == 0){
        return(NULL)
      } else {
        return(vars_trimmed)
      }
    }))
    batch_vars <- plyr::compact(batch_vars)

    split_vars <- setdiff(split_vars, vars_1lvl)

    message("integrate_data(): Removing batch_vars and split_vars with 1 level")
    if (length(batch_vars) == 0) {
      warning("integrate_data(): No batch_vars left. Skipping.")
      batch_vars <- NULL
    } else if (length(split_vars) == 0) {
      warning("integrate_data(): No split_vars left. Skipping.")
      split_vars <- NULL
    }

  }
  
  # Method 1: RNA assay + `RunHarmony()`
  
  batch_vars_integspace_ids <- NULL
  if (!is.null(batch_vars)) {
    
    # Using this separator because some symbols may be automatically converted when stored in the Seurat object
    batch_vars_ids <- lapply(batch_vars, function(x) paste(x, collapse = "zz"))
    batch_vars_integspace_ids <- lapply(batch_vars_ids, function(x) paste0(x, ".harmony"))
    
    batch_vars_len <- length(batch_vars)
    for (i in 1:batch_vars_len) {
      
      # Lower theta for each variable when regressing out multiple variables
      # Default is theta = 2, # `theta` is a positive scalar vector that determines the coefficient of harmony's diversity penalty for each corrected experimental covariate. In challenging experimental conditions, increasing theta may result in better integration results. Theta is an expontential parameter of the diversity penalty, thus setting `theta=0` disables this penalty while increasing it to greater values than 1 will perform more aggressive corrections in an expontential manner. By default, it will set `theta=2` for each experimental covariate.
      
      set.seed(seed_val * i)
      vars_count = length(batch_vars[[i]])
      if (vars_count == 1) {theta_value <- 2} else {theta_value <- rep(0.5, times = vars_count)}
      message("integrate_data(): ", batch_vars_integspace_ids[[i]], " theta = ", paste(theta_value, collapse = "zz"), ". Can change code how to control theta for every group of variables used.")
      obj <- RunHarmony(obj, group.by.vars = batch_vars[[i]], reduction.use = dimred_name, dims.use = dims_use, max_iter = max_iter, reduction.save = batch_vars_integspace_ids[[i]], plot_convergence = TRUE, project.dim = TRUE, verbose = TRUE, theta = theta_value)
      
    }
    
  }
  
  # Method 2: RNA assay + Seurat::IntegrateData(method = CCAIntegration) - only on single variable split_var
  
  seurat_integspace_ids <- NULL
  if (!is.null(split_vars)) {
    
    warning("integrate_data: Review code for Seurat::IntegrateLayers().")
    
    seurat_integspace_ids <- lapply(split_vars, FUN = function(x) paste0(x, ".", tolower(seurat_assay), ".", dimred_name, ".integrated.cca"))
    
    split_vars_len <- length(split_vars)
    for(i in 1:split_vars_len){
      
      split_var <- split_vars[[i]]
      obj_split <- obj
      obj_split[[seurat_assay]] <- split(obj[[seurat_assay]], f = factor(as.character(obj[[split_var]][,1])))
      
      obj_split <- FindVariableFeatures(obj_split, nfeatures = 2000)
      obj_split <- ScaleData(obj_split, vars.to.regress = scaledata_covars, verbose = TRUE)
      obj_split <- RunPCA(obj_split, npcs = 50, verbose = FALSE, seed.use = 290, reduction.name = dimred_name)
      
      obj_split <- IntegrateLayers(object = obj_split, method = CCAIntegration, orig.reduction = dimred_name, new.reduction = seurat_integspace_ids[[i]], assay = seurat_assay, normalization.method = norm_method, dims = dims_use, verbose = TRUE)
      
      # Store integspace for this split_var in the original obj, no need to join layers
      obj[[ seurat_integspace_ids[[i]] ]] <- obj_split[[ seurat_integspace_ids[[i]] ]]
      rm(obj_split)
      
    }
    
  }

  # Can add other complementary integration methods in the future
  
  # (Optional) Add corresponding UMAPs and TSNEs for each integrated space
  
  integspace_names <- unlist(c(batch_vars_integspace_ids, seurat_integspace_ids))
    
  if (get_umap) {
    
    umap_names <- paste0(integspace_names, ".umap")
    for (i in 1:length(umap_names)) {
      obj <- RunUMAP(obj, reduction = integspace_names[i], dims = dims_use, verbose = TRUE, reduction.name = umap_names[i], seed.use = seed_val)
    }
    
  }
  
  if (get_tsne) {
    
    tsne_names <- paste0(integspace_names, ".tsne")
    for (i in 1:length(tsne_names)) {
      obj <- RunUMAP(obj, reduction = integspace_names[i], dims = dims_use, verbose = TRUE, reduction.name = tsne_names[i], seed.use = seed_val)
    }
    
  }
  
  return(obj)
  
}

# Assess before and after integration e.g. use cell metadata (e.g. reference-based mapping) to help assess over-integration

# Plot to show group diversity/mixing/imbalance per label/cluster and arranging based on metric that show imbalance e.g. lisi

# Measure label imbalance given group, can be used to rank label-group based on diversity/imbalance
measure_label_groupdiversity <- function(cellbyfeature_mx, group_df, method = "lisi", summarise_group_fun = median) {
  
  #message("measure_label_groupdiversity(): Only taking group with > 1 unique levels.")
  #num_levels <- apply(group_df, MARGIN = 2, function(x) length(unique(x)))
  #group_df <- group_df[,num_levels > 1]
  
  if (method == "lisi") { # See https://github.com/immunogenomics/LISI
    
    warning("measure_label_groupdiversity(): Scaling of lisi to be comparable across variables likely wrong. Check with developers.")
    # "If the cells are well-mixed, then we might expect the LISI score to be near 2 for a categorical variable with 2 categories."
    lisi_df <- compute_lisi(cellbyfeature_mx, as.data.frame(group_df), colnames(group_df))
    # The higher, the more diverse e.g. e.g. lower imbalance
    summary_measures <- unname(apply(lisi_df, MARGIN = 2, FUN = summarise_group_fun))
    # Scale relative to number of group levels because range of lisi depends on this number
    scale_factors <- vapply(group_df, function(x) length(unique(x)), integer(1))
    summary_measures_scaled <- summary_measures / scale_factors
    output <- cbind(summary_measure = summary_measures, scale_factor = scale_factors, summary_measure_scaled = summary_measures_scaled)
    
  } else {
    message("measure_label_groupdiversity(): Can add more method.")
  }
  
  return(output)
  
}

#library(cowplot)
#library(pheatmap)
#library(RColorBrewer)
plot_label_groupdiversity <- function(label, group_df, num_col = 3, summarise_group_fun = "median", orderby = NULL, same_scale = FALSE) {
  
  group_len <- length(group_df[1,])
  p_lst <- vector("list", group_len)
  summary_groupvar <- vector("numeric", group_len)
  
  # Heatmap showing proportion of each label per group category
  
  for (group_ind in 1:group_len) {
    
    group_name <- colnames(group_df)[[group_ind]]
    group <- group_df[[group_ind]]
    message("plot_label_groupdiversity(): Plotting ", group_name, "...")
    
    tbl <- table(label, group)
    prop_tbl <- proportions(tbl, margin = 1)
    summary_groupvar[group_ind] <- do.call(summarise_group_fun, list(apply(prop_tbl, MARGIN = 1, FUN = function(x) sd(x) / mean(x))))
    # Coefficient of variation because the proportion given equal distribution of labels across groups depends on number of unique group levels. If group has 2 levels, equal distribution means 0.5 proportion, but if group has 10 levels, equal distribution means 0.1 proportion.
    plot_title = paste0(group_name, ": ", summarise_group_fun, " cv = ", round(summary_groupvar[group_ind], 4), "; Given equal distribution, proportion per level = ", 1 / length(unique(group)))
    
    pheatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100) # Default
    if (length(unique(prop_tbl)) == 1) {prop_tbl <- tbl}

    if (same_scale) {
      
      extendPaletteFUN = colorRampPalette(brewer.pal(11, "Spectral"))
      breakvals = seq(0, 1, by = 0.05)
    
      p_lst[[group_ind]] <- pheatmap::pheatmap(prop_tbl, display_numbers = tbl, cluster_rows = FALSE, cluster_cols = FALSE, main = plot_title, fontsize = 5, silent = TRUE, breaks = breakvals, color = rev(extendPaletteFUN(length(breakvals) - 1)))[[4]] # Access for plot_grid to work
      
    } else {
      p_lst[[group_ind]] <- pheatmap::pheatmap(prop_tbl, display_numbers = tbl, cluster_rows = FALSE, cluster_cols = FALSE, main = plot_title, fontsize = 5, silent = TRUE, color = pheatmap_color)[[4]] # Access for plot_grid to work
    }
    
    
  }
  
  # (Optional) Sort plots based on vector of numeric values
  
  if (!is.null(orderby) & is.numeric(orderby)) {
    p_lst <- p_lst[order(orderby)]
    message("plot_label_groupdiversity(): Sorted by orderby argument.")
  } else {
    p_lst <- p_lst[order(summary_groupvar, decreasing = TRUE)]
    message("plot_label_groupdiversity(): Sorted by summary statistic value on plot.")
  }
  plots <- plot_grid(plotlist = p_lst, ncol = num_col)
  return(plots)

}

# Method 2: Plot pca/harmony, umaps, tsne colouring by categorical and continuous covariates, and by clustering and other informative labels e.g. from reference-based mapping

#library(cowplot)
#library(scater)
plot_on_dimreds <- function(sce, label_name, dimred_names = NULL, dimred_name_pattern = NULL, dims_use = 1:2, num_col = 3) {
  
  if (!is.null(dimred_name_pattern)) {
    dimred_names <- reducedDimNames(sce)
    dimred_names <- dimred_names[grepl(dimred_name_pattern, dimred_names, fixed = FALSE)]
  }
  
  p_lst <- lapply(dimred_names, function(nme) {
    plotReducedDim(sce, dimred = nme, ncomponents = dims_use, colour_by = label_name)
  })
  plot_grid(plotlist = p_lst, ncol = num_col)
  
}
