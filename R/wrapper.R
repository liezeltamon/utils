suppressPackageStartupMessages({
  require(tidyverse)
  require(scran)
  require(Seurat)
  require(irlba)
  require(clustree)
})

create_seuFromSce <- function(sce){
  seu <- CreateSeuratObject(counts = counts(sce),
                            meta.data = as.data.frame(colData(sce)),
                            min.cells = 0,
                            min.features = 0)
  return(seu)
}

create_sceFromSeu <- function(seu){
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(seu, slot = "counts"),
                                            logcounts= GetAssayData(seu, slot = "data")),
                              reducedDims = list(PCA = Embeddings(seu, reduction = "pca"),
                                                 UMAP = Embeddings(seu, reduction = "umap")),
                              colData = seu[[]])
  return(sce)
}

getModuleScores <- function(obj){
  obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
  return(obj)
}

get_denoisedPCs <- function(expr_mat, block, subset.row, seed = 290){
  set.seed(seed)
  var_dec <- scran::modelGeneVar(expr_mat, block = block)
  out_lst <- scran::getDenoisedPCs(mx, technical = var_dec, subset.row = subset.row)
  
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
  clusters <- FindClusters(knn$snn, verbose = FALSE, method = alg_num, resolution = resolution) %>%
    as_tibble(rownames = "cell")
  colnames(clusters) <- c("cell", paste0(method, "_", as.character(resolution)))
  
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