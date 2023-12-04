suppressPackageStartupMessages({
  require(BiocParallel)
  require(DelayedMatrixStats)
  require(DropletUtils)
  require(clustree)
  require(irlba)
  require(scater)
  require(scran)
  require(Seurat)
  require(tidyverse)
})

# Check

## Copied from https://github.com/Bioconductor/SummarizedExperiment/compare/6011d3e62f4fe9c41e4485478eb5cc0a881a0032..4aa706a2124209fee9074a029bd0829d97ddb995
splitByCol <- function(x, f, drop=FALSE) {
  by.col <- split(seq_len(ncol(x)), f, drop=drop)
  out <- lapply(by.col, function(i) x[,i])
  List(out)
}

#

create_seuFromSce <- function(sce, min.cells = 0, min.features = 0){
  seu <- CreateSeuratObject(counts = counts(sce),
                            meta.data = as.data.frame(colData(sce)),
                            min.cells = min.cells,
                            min.features = min.features)
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
  
do_basicQC <- function(counts_mat = NULL, sce = NULL){ # Expects rownames to be gene symbols
  if(is.null(counts_mat)){ counts_mat <- counts(sce) }
  mito_genes <- rownames(counts_mat)[grep("^MT-", rownames(counts_mat))]
  ribo_genes <- rownames(counts_mat)[grep("^RP[SL]", rownames(counts_mat))]
  qc_tbl <- perCellQCMetrics(counts_mat, subsets = list(mito_genes = mito_genes, ribo_genes = ribo_genes))
  
  if(!is.null(sce)){
    colData(sce) <- cbind(colData(sce), qc_tbl); return(sce)
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
  
  return(p)
}

getModuleScores <- function(obj){
  obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
  return(obj)
}


plot_variableContribution <- function(sce,
                                      covars = c("Sample", "condition", "nCount_RNA", "nFeature_RNA", "subsets_mito_genes_percent", "subsets_ribo_genes_percent", "S.Score","G2M.Score", "CC.Difference", "Phase"),
                                      subset_row, nCPU,
                                      dimred, n_dimred # For scater::getExplanatoryPCs()
                                      ){
  covars_varexp <- scater::getVarianceExplained(sce, variables = covars, subset_row = subset_row, BPPARAM = MulticoreParam(nCPU))
  p_varexp <- scater::plotExplanatoryVariables(covars_varexp, nvars_to_plot = Inf)
  exp_pcs <- scater::getExplanatoryPCs(sce, dimred = dimred, n_dimred = n_dimred, variables = covars)
  p_exp_pcs <- plotExplanatoryPCs(exp_pcs)
  return(list(p_varexp, p_exp_pcs))
}

get_denoisedPCs <- function(expr_mat, subset.row, block = NULL, seed = 290){
  set.seed(seed)
  var_dec <- scran::modelGeneVar(expr_mat, block = block)
  out_lst <- scran::getDenoisedPCs(expr_mat, technical = var_dec, subset.row = subset.row)
  return(out_lst)
}

# Modified from quadbio/organoid_regulomes
get_clusters <- function(expr_mat,
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

get_clustree <- function(clusters_mat, prefix, ggsave_path = NULL, height = 30, width = 30){
  ctree_obj <- clustree(clusters_mat, prefix = prefix)
  # Reverse tree so at the top is the lowest resolution i.e. larger clusters
  ctree_obj$data[[prefix]] <- factor(as.character(ctree_obj$data[[prefix]]),
                                     levels = rev(levels(ctree_obj$data[[prefix]])))
  if(!is.null(ggsave_path)){
    ggsave(ggsave_path, plot = ctree_obj, height = height, width = width)
  }
  return(ctree_obj)
}

# quadbio/organoid_regulomes - https://github.com/quadbio/organoid_regulomes/blob/1e3b32f067d5d190e8ba4709e34a3ae243b5f9db/utils/scripts/wrapper.R