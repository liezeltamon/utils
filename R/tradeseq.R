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
