# Trajectory based differential expression analysis - tradeSeq wrapper functions

library(tradeSeq)

# Run evaluate N times and summarise result to robustly identify k as recommended by developers
run_evaluatek_ntimes <- function(sce, k = 3:10, n = 3, seed_val = 290, ncpu = 1, 
                                 plot_path = NULL, # Size set assumming this is pdf
                                 # Or vector of condition colData from sce
                                 conditions = NULL) {
  
  if (class(sce) != "SingleCellExperiment") {
    stop("run_evaluatek_ntimes(): Currently accepts only sce output from Slingshot. Yet to be modified to be more flexible.")
  } else {
    
    parallel <- ifelse(ncpu > 1, TRUE, FALSE)
    bpparam = BiocParallel::MulticoreParam(ncpu, RNGseed = seed_val * 4)
    
    if (!is.null(plot_path)) {pdf(plot_path, width = 10, height = 5)}
    
    runs <- lapply(1:n, function(run_ind){
      set.seed(seed_val * run_ind)
      eval_k_mx <- evaluateK(counts = as.matrix(assays(sce)$counts),
                             k = k,
                             nGenes = 500, # Default 500, chosen randomly
                             pseudotime = pathStats(colData(sce)$slingshot)$pseudotime,
                             cellWeights = pathStats(colData(sce)$slingshot)$weights,
                             conditions = conditions,
                             plot = TRUE,
                             parallel = parallel,
                             BPPARAM = bpparam)
    })
    names(runs) <- paste0("run", 1:n, "_seed", seed_val * 1:n)
    
    if (!is.null(plot_path)) {dev.off()}
    
  }
    
  return(runs)
  
}

# Run all possible differential expression tests from tradeSeq depending on data (e.g. number of lineages/paths) then return list of outputss
run_de_tests_tradeseq <- function(
    models,
    num_paths, # To decide whether to run between-lineage comparisons, can be determined from SCE slingshot output but removing that code makes function applicable to any output of fitGAM()
    l2fc = log2(1), # Default
    lineages = FALSE, # Default, but preferable to set to TRUE
    p_adj_method = NULL, # See p.adjust() methods. If NULL, do not adjust
    
    within_pairwise = FALSE, # Default, set to TRUE For conditionTest() to return pairwise comparisons between conditions
    
    # Specific for between-lineage tests
    between_pairwise = FALSE # Default
  ) {
  
  is_fitgam_withconditions = any(grepl("conditions", colnames(models$tradeSeq)))
  #num_paths = try(ncol(pathStats(models$slingshot)$pseudotime))
  num_tests = 6
  outputs <- setNames(as.list(rep(NA, times = num_tests)),
                     c("association", "condition", "startVsEnd", "diffEnd", "pattern", "earlyDE"))
  
  # Within-lineage comparisons
  
  # Association of gene expression with pseudotime
  outputs[["association"]] <- associationTest(models, lineages = lineages, l2fc = l2fc
                                             # contrastType = "start", # Default
                                             # nPoints arguments for test
                                             )
  # Gene expression across conditions within lineage
  if (is_fitgam_withconditions) {
    message("run_de_tests_tradeseq(): fitGAM() specified with conditions, doing conditionTest()...")
    outputs[["condition"]] <- conditionTest(models, lineages = lineages, l2fc = l2fc,
                                            # To return pairwise comparisons between conditions
                                            pairwise = within_pairwise)
  }
  
  # Discovering progenitor marker genes
  outputs[["startVsEnd"]] <- startVsEndTest(models, lineages = lineages, l2fc = l2fc
                                            # pseudotimeValues, # To compare two pseudotime points
                                            )
  
  # Between-lineage comparisons
  
  if (num_paths > 1) {
    
    # Discovering differentiated cell type markers
    outputs[["diffEnd"]]  <- diffEndTest(models, l2fc = l2fc,
                                         pairwise = between_pairwise, # Can set to TRUE if > 2 lineages
                                         )
    # Discovering genes with different expression patterns
    outputs[["pattern"]] <- patternTest(models, pairwise = between_pairwise, l2fc = l2fc)
    outputs[["earlyDE"]] <- earlyDETest(models, pairwise = between_pairwise, l2fc = l2fc)
    
  }
  
  # (Optional) Adjust p-values - look for all columns that contain p-values then create extra columns containing corresponding p-adjusted values
  
  outputs <- lapply(outputs, function(test_df) {
    
    if (is.data.frame(test_df)) {
      pval_col_inds <- grep("pvalue", colnames(test_df))
      padjust_mx <- apply(test_df[,pval_col_inds, drop = FALSE], MARGIN = 2, function(col) {
        p.adjust(col, method = p_adj_method)
      })
      colnames(padjust_mx) <- gsub("pvalue", "padjust", colnames(test_df)[pval_col_inds])
      return(cbind(test_df, padjust_mx))
    } else {
      return(NA)
    }
    
  })
  
  # Order genes based on main (global, not from each lineage) p-values or p adjusted (if available)
  
  outputs <- lapply(outputs, function(test_df) {
    
    if (is.data.frame(test_df)) {
      orderby_colnme <- ifelse("padjust" %in% colnames(test_df), "padjust", "pvalue")
      return(test_df[order(-test_df[[orderby_colnme]]),])
    } else {
      return(NA)
    }
    
  })
    
  return(outputs)
  
}

# Bin cells (e.g. based on pseudotime) and calculate a value per bin
get_value_per_bin <- function (mx, bins, fun_per_bin = mean) {
  
  message('function1(): Set fun_per_bin = "percent_expressing" to use built-in function')
  
  if (!is.factor(bins)) {
    warning("function1(): bins should an ordered factor. Converting...")
    bins <- factor(as.character(bins))
  }
  
  if (identical(fun_per_bin, "percent_expressing")) {
    fun_per_bin <- function (x) {
      sum(x > 0) / length(x) * 100
    }
  }
  
  mx_transposed <- t(mx)
  mx_transposed_agg <- aggregate(as.matrix(mx_transposed), by = list(bins = bins), FUN = fun_per_bin)
  value_per_bin_mx <- t(mx_transposed_agg[,colnames(mx_transposed_agg) != "bins"])
  
}

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyr)

# Plot
plot_value_per_bin <- function (fill_mx, colour_mx, row_features_order, row_group, 
                                scale_fill = TRUE, scale_colour = TRUE) {
  
  if (!identical(dimnames(fill_mx), dimnames(colour_mx))) {
    stop("function2(): Dimnames of fill_mx and colour_mx not identical")
    rm(fill_mx)
  } else {
    dimnames_frommx <- dimnames(fill_mx)
  }
  
  # Scale
  
  scale_fun <- function(data){
    (data - mean(data)) / sd(data)
    #min_value <- min(data)
    #max_value <- max(data)
    #scaled_data <- (data - min_value) / (max_value - min_value)
  }
  
  if (scale_fill) {
    fill_mx <- t(apply(fill_mx, MARGIN = 1, scale_fun))
    fill_mx[!is.finite(fill_mx)] <- NA
  }
  
  if (scale_colour) {
    colour_mx <- t(apply(colour_mx, MARGIN = 1, scale_fun))
    colour_mx[!is.finite(colour_mx)] <- NA
  }
  
  dimnames(fill_mx) <- dimnames_frommx
  dimnames(colour_mx) <- dimnames_frommx
  
  # Prepare plot data
  
  fill_df <- as.data.frame(fill_mx)
  fill_df$group <- factor(row_group, levels = as.character(unique(row_group)))
  tidy_df <- fill_df %>% 
    rownames_to_column("feature") %>% 
    pivot_longer(-c(feature, group), names_to = "bin", values_to = "value_fill") %>% 
    unite("feature_bin", feature, bin, sep = "_", remove = FALSE)
  tidy_df$bin <- factor(tidy_df$bin, levels = setdiff(colnames(fill_df), "group"))
  tidy_df$feature <- factor(tidy_df$feature, levels = row_features_order)
  
  colour_df <- colour_mx %>% 
    as.data.frame() %>% 
    rownames_to_column("feature") %>% 
    pivot_longer(-feature, names_to = "bin") %>% 
    unite("feature_bin", feature, bin, sep = "_", remove = FALSE)
  
  if (identical(tidy_df$feature_bin, colour_df$feature_bin)) {
    tidy_df$value_colour <- colour_df$value
    rm(colour_df)
  }
  
  # Heatmap
  
  p_params <- list(fill = list(pal_colours = brewer.pal(11, "RdYlBu")),
                   colour = list(pal_colours = rev(brewer.pal(9, "Greys"))))
  
  for (aest in names(p_params)) {
    
    p_params[[aest]]$extendPaletteFUN <- colorRampPalette(p_params[[aest]]$pal_colours)
    p_params[[aest]]$vals_range <- range(boxplot.stats(tidy_df[[paste0("value_", aest)]])$stats)
    p_params[[aest]]$breakvals <- unique(seq(p_params[[aest]]$vals_range[1],
                                             max(p_params[[aest]]$vals_range[2]), length.out = 50))
    p_params[[aest]]$plot_colors <- rev(p_params[[aest]]$extendPaletteFUN(
      length(p_params[[aest]]$breakvals) - 1
    ))
    
  }
  
  p <- ggplot(tidy_df, aes(x = bin, y = feature, fill = value_fill, color = value_colour)) +
    geom_tile(size = 0.5) +  # size controls the thickness of the lines
    scale_fill_gradientn(colors = p_params$fill$plot_colors, breaks = p_params$fill$breakvals, 
                         limits = p_params$fill$vals_range, oob = scales::squish) +
    scale_colour_gradientn(colors = p_params$colour$plot_colors, breaks = p_params$colour$breakvals, 
                           limits = p_params$colour$vals_range, oob = scales::squish) +
    #facet_grid(group ~ ., scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90), legend.text = element_text(size = 1)) +
    theme_minimal()
  
  return(p)
  
}
