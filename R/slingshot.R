# Trajectory inference - Slingshot wrapper functions

library(scater)
library(SingleCellExperiment)
library(slingshot)

# Reproducibility note - No need to set seed but object will be different because of mst igraph object from slingshot output, reproducibility comfirmed by comparing outputs after removing mst igraph object
run_slingshot <- function(
    data,
    labels = NULL,
    dimred_name = NULL,
    # Specify start and end clusters for semi-supervised approach
    start.clus = NULL, # Developers "generally recommend the specification of an initial cluster based on prior knowledge"
    end.clus = NULL, # "Clusters which are specified as terminal cell states will be constrained to have only one connection when the MST is constructed
    approx_points = NULL, # Default i.e. min(150, number of cells). If numeric N, will be min(150, N). If FALSE (or 0), will use all cells. Consider setting if working with large datasets.
    ... # Pass arguments to slingshot(), For advanced control, see vignette how to change parameters like dist.method and omega
    ) {
  
  if (is.matrix(data)) {
    dimred_name <- NULL
    message("run_slingshot(): dimred_name not applicable when data is a matrix.")
  }
  
  if (is.null(labels)) {
    output <- slingshot(data, reducedDim = dimred_name, start.clus = start.clus, end.clus = end.clus, approx_points = approx_points, ...)
  } else {
    output <- slingshot(data, clusterLabels = labels, reducedDim = dimred_name, start.clus = start.clus, end.clus = end.clus, approx_points = approx_points, ...)
  }
  
  return(output)
  
}

# table(slingshot::slingBranchID()) for tabulating how many cells belong to specific or multiple lineages

plot_slingshot_curves <- function(
    data,
    dimred_name = NULL,
    new_dimred_name = NULL,
    colour_by = "slingPseudotime_1"
    ) {
  
  if (class(data) == "SingleCellExperiment") {
    
    if (!is.null(new_dimred_name)) {
      embedded <- embedCurves(data, new_dimred_name)
      dimred_name <- new_dimred_name
    } else {
      embedded <- data
    }
    curves <- slingCurves(embedded, as.df = FALSE)
    
  } else {
    stop("plot_slingshot_curves(): Code for handling non-SCE input data not done yet. Use base R as in vignette.")
  }
  
  # Modified from OSCA book
  p_lst <- list()
  
  curves_len <- length(curves) # i.e. number of lineages
  if (grepl("slingPseudotime", colour_by) & curves_len > 1) {
    data$slingPseudotime_avelineages <- rowMeans(slingPseudotime(sce), na.rm = TRUE)
    colour_by = "slingPseudotime_avelineages"
  }
  p_combined <- plotReducedDim(data, colour_by = colour_by, dimred = dimred_name)
  
  for (i in 1:curves_len) {
    
    path <- curves[[i]]
    curve_path <- data.frame(path$s[path$ord,])
    colnmes <- colnames(curve_path)
    
    if (grepl("slingPseudotime", colour_by)) {
      colour_by = paste0("slingPseudotime_", i)
    }
    
    p_lst[[i]] <- plotReducedDim(data, colour_by = colour_by, dimred = dimred_name) +
      geom_path(data = curve_path, aes(x = .data[[colnmes[1]]], y = .data[[colnmes[2]]]), size = 1.2)
    
    p_combined <- p_combined + 
      geom_path(data = curve_path, aes(x = .data[[colnmes[1]]], y = .data[[colnmes[2]]]), size = 1.2)
    
  }
  
  p_lst[[curves_len + 1]] <- p_combined
  p_grid <- plot_grid(plotlist = p_lst, nrow = 1)
  
  return(p_grid)
  
}

pairwise_ks_test <- function (value, group, 
                              alternative = "two.sided", 
                              p_adj_method = NULL,
                              # Pass arguments to ks.test()
                              ...) {
  
  stop("pairwise_ks_tes(): Yet to add code!")

}
