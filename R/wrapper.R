suppressPackageStartupMessages({
  library(assertthat)
  library(BiocParallel)
  library(bluster)
  library(DelayedMatrixStats)
  library(DropletUtils)
  library(clustree)
  library(irlba)
  library(scater)
  library(scran)
  library(Seurat)
  library(tidyverse)
  # For using leidenalg python module for Seurat::FindClusters(), need RETICULATE_PYTHON to be set to python environment path in .Renviron
  #library(reticulate)
  # Ignore error in case no env with these packages yet
  try(reticulate::import("numpy"))
  try(reticulate::import("pandas"))
  try(reticulate::import("leidenalg"))
})

# Check

## Copied from https://github.com/Bioconductor/SummarizedExperiment/compare/6011d3e62f4fe9c41e4485478eb5cc0a881a0032..4aa706a2124209fee9074a029bd0829d97ddb995
splitByCol <- function(x, f, drop=FALSE) {
  by.col <- split(seq_len(ncol(x)), f, drop=drop)
  out <- lapply(by.col, function(i) x[,i])
  List(out)
}

#

#create_seuFromSce <- function(sce, min.cells = 0, min.features = 0){
create_seuFromSce <- function(sce, ...){
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)), ...)
  return(seu)
}

make_uniqueAsInSeu <- function(sce){
  gsymbols_uniq <- make.unique(rowData(sce)$Symbol)
  rownames(sce) <- gsymbols_uniq
  return(sce)
}

create_sceFromSeu <- function(seu){
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(seu, slot = "counts"),
                                            logcounts= GetAssayData(seu, slot = "data")),
                              reducedDims = list(PCA = Embeddings(seu, reduction = "pca"),
                                                 UMAP = Embeddings(seu, reduction = "umap")),
                              colData = seu[[]])
  return(sce)
}

# Based on OSCA book
plot_barcodeRanks <- function(counts_mat, sample_id = NULL){
  bcrank <- DropletUtils::barcodeRanks(counts_mat)
  # Only showing unique points for plotting speed
  uniq <- !duplicated(bcrank$rank)
  p <- ggplot(as.data.frame(bcrank[uniq,]), aes(rank, total)) +
    geom_point() +
    geom_hline(yintercept = c(metadata(bcrank)$knee, metadata(bcrank)$inflection), 
               colour = c("dodgerblue", "darkgreen")) + 
    scale_x_log10() +
    scale_y_log10() + 
    labs(lab = "log10(rank)", ylab = "log10(UMI count)", 
         title = paste0(sample_id, " knee-blue-", metadata(bcrank)$knee, " \n inflection-green-", metadata(bcrank)$inflection)) +
    theme(plot.title = element_text(size = 10))
  return(p)
}

plot_ambientPval <- function(emptyDrops_out, emptyDrops_lower = 100, seed = 290){
  set.seed(seed)
  p <- ggplot(as.data.frame(emptyDrops_out[emptyDrops_out$Total > 0 & emptyDrops_out$Total <= emptyDrops_lower,]), aes(PValue)) + 
    geom_histogram(fill = "gray70", col = "black") +
    labs(x = "Pvalue")
  return(p)
}
  
do_basicQC <- function(counts_mat = NULL, 
                       sce = NULL,
                       feature_names = NULL, # If NULL, uses row names
                       ...){
  
  assert_that(!(!is.null(counts_mat) & !is.null(sce)),
              msg = "do_basicQC(): Supply either counts_mat or sce")
  
  if (!is.null(sce)) {
    mainExpName(sce) <- "Gene Expression"
    if (is.null(feature_names)) {feature_names <- rownames(sce)}
    input_x <- sce
    rownames(input_x) <- feature_names
  }
  
  if (!is.null(counts_mat)) {
    if (is.null(feature_names)) {feature_names <- rownames(counts_mat)}
    input_x <- counts_mat
  }
  
  assert_that(!any(duplicated(feature_names)),
              msg = "do_basicQC(): Feature names must be unique")
  mito_genes <- feature_names[grep("^MT-", feature_names)]
  ribo_genes <- feature_names[grep("^RP[SL]", feature_names)]
  
  if (is.null(counts_mat)) { counts_mat <- counts(sce) }
  
  qc_tbl <- perCellQCMetrics(input_x,
                             subsets = list(mito_genes = mito_genes, 
                                            ribo_genes = ribo_genes),
                             ...)
  
  if (!is.null(sce)) {
    colData(sce) <- cbind(colData(sce), qc_tbl)
    return(sce)
  } else {
    return(qc_tbl)
  }
  
}

plot_basicQC <- function(qc_tbl, group = "Sample", metrics = c("sum", "detected"), plot_type = "violin", transform_log10 = TRUE){
  qc_tbl <- qc_tbl[,c(group, metrics)] %>% 
    as.data.frame %>% 
    rownames_to_column(var = "Barcode_uniq")
  
  colnames(qc_tbl)[colnames(qc_tbl) == group] <- "group"
  
  qc_tbl <- qc_tbl %>% 
    tidyr::pivot_longer(cols = -c(Barcode_uniq, group))
  
  if(transform_log10){ qc_tbl$value <- log10(qc_tbl$value) }
  if(plot_type == "histogram"){
    p <- ggplot(qc_tbl, aes(value)) +
      geom_histogram(bins = 100, fill = "gray70", col = "black") + 
      facet_grid(group~name)
  } else if(plot_type == "violin"){
    p <- ggplot(qc_tbl, aes(group, value)) +
      geom_violin() +
      facet_grid(.~name)
  } else {
    stop("plot_basicQC(): plot_type not recognised.")
  }
  
  p <- p +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  
  return(p)
}

get_outliers <- function (x,
                          # If NULL, treat x as one sample.
                          # Rule of thumb is get outliers per sample so we assume that
                          # the sample are mostly good cells. If data hashed, dehash then
                          # identify outliers per sample, this way each sample assumed to be
                          # mostly good cells.
                          sample_ids = NULL,
                          is_outlier_fields = NULL,
                          add_sum_fields = "altexps_Antibody Capture_sum",
                          add_detected_fields = "altexps_Antibody Capture_detected",
                          is_outlier_nmads = 3
                          ) {
  
  
  if (is.null(is_outlier_fields)) {
    is_outlier_fields <- c("low_lib_size",
                           "low_n_features",
                           "high_subsets_mito_genes_percent")
  }
  
  # Initiate output matrix
  
  assert_that(!any(duplicated(rownames(x))),
              msg = "get_outliers(): Duplicated barcodes")
  outlier_mx <- matrix(NA, nrow = nrow(x), ncol = length(is_outlier_fields) + 1,
                       dimnames = list(rownames(x),
                                       c(is_outlier_fields, "outlier")))
  
  # Identify outlier per level of sample_id column
  
  for (sample_id_lvl in unique(sample_ids)) {
    
    message("get_outliers(): sample_id = ", sample_id_lvl, "...")
    
    lvl_inds <- which(sample_ids == sample_id_lvl)
    x_sub <- x[lvl_inds, ]
    out_a <- scuttle::perCellQCFilters(x_sub, sub.fields = TRUE, nmads = is_outlier_nmads)
    out_a$discard <- NULL
    
    # Apply on additional sum and detected fields
    
    #add_sum_fields <- sort(grep("_sum", colnames(x), value = TRUE))
    #add_detected_fields <- sort(grep("_detected", colnames(x), value = TRUE))
    len <- unique(length(add_sum_fields), length(add_detected_fields))
    out_b_lst <- lapply(1:len, function (i) {
      
      out_tmp <- scuttle::perCellQCFilters(x = x_sub,
                                           sum.field = add_sum_fields[i],
                                           detected.field = add_detected_fields[i],
                                           sub.fields = NULL,
                                           nmads = is_outlier_nmads
                                           )
      out_tmp$discard <- NULL
      colnames(out_tmp) <- paste0(c(add_sum_fields[i], add_detected_fields[i]), ".",
                                  colnames(out_tmp))
      return(out_tmp)
      
    })
    out_b <- do.call("cbind", out_b_lst)
    
    out_sub_mx <- as.matrix(cbind(out_a, out_b))
    out_sub_mx <- out_sub_mx[, is_outlier_fields]
    out_sub_mx <- cbind(out_sub_mx,
                        outlier = rowSums(as.matrix(out_sub_mx)) > 0)
    outlier_mx[lvl_inds, ] <- as.matrix(out_sub_mx)
    
    rm(out_sub_mx, lvl_inds)
    
  }
  
  return(outlier_mx)
  
}

# Return list of genes based on criteria
filter_features <- function(mx, min_num_samples = NULL, min_perc_samples = NULL, drop_mito = FALSE, drop_ribo = FALSE){
  
  samples_len <- ncol(mx)
  
  if (is.null(min_perc_samples)) {
    
    # Roughly get features expressed (non-zero) in 1 per thousand of cells
    min_num_samples <- ifelse(is.null(min_num_samples), ceiling(samples_len * 0.001), min_num_samples)
    num_samples_withfeat <- scuttle::numDetectedAcrossCells(mx, rep("group", times = samples_len), threshold = 0)
    num_samples_withfeat <- assay(num_samples_withfeat, "sum")[,"group"]
    features_left <- names(num_samples_withfeat[num_samples_withfeat >= min_num_samples])
    message("filter_features(): Filtering based on min_num_samples = ", min_num_samples)
    
  } else {
    feature_metrics_dframe <- perFeatureQCMetrics(mx, threshold = 0)
    features_left <- rownames(feature_metrics_dframe)[feature_metrics_dframe$detected >= min_perc_samples]
    message("filter_features(): Filtering based on min_perc_samples = ", min_perc_samples)
  }
  
  if (drop_mito) {
    features_left <- features_left[!grepl("^MT-", features_left)]
  }
  
  if (drop_ribo) {
    features_left <- features_left[!grepl("^RP[SL]", features_left)]
  }
  
  return(features_left)
  
}

getModuleScores <- function(obj){
  obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
  return(obj)
}

get_denoisedPCs <- function(expr_mat, subset.row, block = NULL, seed = 290){
  set.seed(seed)
  var_dec <- scran::modelGeneVar(expr_mat, block = block)
  out_lst <- scran::getDenoisedPCs(expr_mat, technical = var_dec, subset.row = subset.row)
  return(out_lst)
}

# Modified from quadbio/organoid_regulomes
get_clusters <- function(expr_mat, # observations x variables
                         do_pca,
                         dims_use,
                         # If applying pca on input expression matrix
                         scale_pca,
                         n_pcs,
                         # FindClusters() parameters
                         resolution = 0.8,
                         algorithm = "louvain",
                         # FindClusters(method = "matrix") by default, change to "igraph" if having memory issues, see https://github.com/satijalab/seurat/issues/2294
                         method_leiden = "matrix", 
                         method = "leiden" # Deprecated, method confused with algorithm argument in FindClusters()
                         ){
  
  warning("Need to fix when inputting expression matrix instead.")
  
  resolution <- sort(unique(resolution))
  
  if(do_pca){
    message("get_clusters(): Applying pca on input matrix...")
    pca_mat <- irlba::prcomp_irlba(expr_mat, n = n_pcs, scale. = scale_pca)$x
    rownames(pca_mat) <- rownames(expr_mat)
  } else {
    pca_mat <- expr_mat
  }
  message("get_clusters(): Using first ", max(dims_use), " dimensions for FindNeighbors()..")
  pca_mat <- pca_mat[,dims_use]
  knn <- FindNeighbors(pca_mat, verbose = TRUE)
  
  if (algorithm %in% c("leiden", "louvain")) {
    
    if(algorithm == "leiden"){
      alg_num <- 4
    } else if (algorithm == "louvain"){
      alg_num <- 1
    }
    clusters <- FindClusters(knn$snn, verbose = TRUE, algorithm = alg_num, resolution = resolution, method = method_leiden)
    
  } else if (algorithm == "HclustParam") {
    
    clusters <- as.data.frame(matrix(NA, nrow = nrow(pca_mat), ncol = length(resolution), 
                                     dimnames = list(rownames(pca_mat), as.character(resolution))))
    hclust_tree <- bluster::clusterRows(pca_mat, BLUSPARAM = HclustParam(method = "ward.D2"), full = TRUE)$objects$hclust
    for (num_cluster in resolution) {
      
      x <- cutree(hclust_tree, k = num_cluster)
      clusters[, as.character(num_cluster)] <- factor(as.character(x), levels = as.character(sort(unique(x))))
      
    }
   
  } else {
    stop("get_clusters(): Invalid algorithm")
  }
  
  colnames(clusters) <- paste0(algorithm, "_", as.character(resolution))
  return(clusters)
  
}

get_clusters_v1 <- function(expr_mat,
                            do_pca = TRUE,
                            scale_pca,
                            n_pcs = 50,
                            method = "leiden",
                            resolution = 0.8){
                         
  if(do_pca){
    pca_mat <- irlba::prcomp_irlba(expr_mat, n = n_pcs, scale. = scale_pca)$x
    rownames(pca_mat) <- rownames(expr_mat)
  } else {
    pca_mat <- expr_mat
  }
  if(method == "leiden"){
    alg_num <- 4
  } else if (method == "louvain"){
    alg_num <- 1
  }
  knn <- FindNeighbors(pca_mat, verbose = FALSE)
  clusters <- FindClusters(knn$snn, verbose = FALSE, method = alg_num, resolution = resolution)
  colnames(clusters) <- paste0(method, "_", as.character(resolution))
  return(clusters)
}

get_clustree <- function(clusters_mat, prefix, ggsave_path = NULL, height = 30, width = 30, ...){
  ctree_obj <- clustree(clusters_mat, prefix = prefix, ...)
  # Reverse tree so at the top is the lowest resolution i.e. larger clusters
  ctree_obj$data[[prefix]] <- factor(as.character(ctree_obj$data[[prefix]]),
                                     levels = rev(levels(ctree_obj$data[[prefix]])))
  if(!is.null(ggsave_path)){
    ggsave(ggsave_path, plot = ctree_obj, height = height, width = width)
  }
  return(ctree_obj)
}

# quadbio/organoid_regulomes - https://github.com/quadbio/organoid_regulomes/blob/1e3b32f067d5d190e8ba4709e34a3ae243b5f9db/utils/scripts/wrapper.R

# Modified to return rownames
UpSetRfromList_rowsnamed <- function(input){
  # UpSetR::fromList code
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){x <- as.vector(match(elements, x))}))
  data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  rownames(data) <- elements # Added
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  return(data)
}
