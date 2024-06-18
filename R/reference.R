# Functions to call different reference datasets, returning as SCE object

suppressPackageStartupMessages({
  require(BiocParallel)
  require(data.table)
  require(dplyr)
  require(DelayedMatrixStats)
  require(LoomExperiment)
  require(readr)
  require(scuttle)
  require(Seurat)
  require(Signac)
  require(SingleCellExperiment)
  require(tidyr)
})

check_specificityLabels <- function(mdta_tbl = colData(dta_raw), 
                                    grouping_label = "Clusters",
                                    toCheck_labels = c("TopLevelCluster", "CellClass")){
  tmp_dt <- as.data.table(mdta_tbl[,c(grouping_label, toCheck_labels)])
  tmp_uniq_dt <- tmp_dt[!duplicated(tmp_dt),]
  tmp_uniq_perGroupingLabel <- split(tmp_uniq_dt, f = as.character(tmp_uniq_dt[[grouping_label]]))
  class_len_perGroupingLabel <- unname(unlist(lapply(tmp_uniq_perGroupingLabel, FUN = nrow)))
  if(!all(class_len_perGroupingLabel == 1)){
    warning(paste0("Some in ", paste(toCheck_labels, collapse = ","), " are NOT unique per ", grouping_label))
  } else{
    message(paste0(paste(toCheck_labels, collapse = ","), " are unique per ", grouping_label))
  }
}
  
# Fleck, ..., Treutlein, Nature 2022 (https://doi.org/10.5281/zenodo.5242913)

FleckBrainOrganoidData <- function(dataset = "atlas_reference/fleck_data/RNA_ATAC_metacells_srt.rds",
                                   possible_labels = c("nowakowski_prediction", "RNA_snn_res.1", "seurat_clusters", "RNA_snn_res.0.8", "pred_class", "pred_region", "RNA_snn_res.2", "patterning_region", "stage", "neuron_type", "is_mesoderm", "neurogen_region", "RNA_snn_res.10", "RNA_snn_res.20", "lineage", "is_root", "lineage_coarse"),
                                   check_labels = FALSE){
  
  dta_path <- file.path("data-raw", dataset)
  dta_raw <- readRDS(dta_path)
  
  dta_raw <- DietSeurat(dta_raw, layers = c("counts", "data"), assays = "RNA", dimreducs = c("pca", "cellrank"))
  
  # QC: > 0 library size
  
  is_libsize_0 <- dta_raw$nCount_RNA <= 0
  dta <- dta_raw[,!is_libsize_0]
  message(paste0("FleckBrainOrganoidData(): ", sum(is_libsize_0), " cells with 0 library size removed."))
  
  # Determine specificity and hierarchy of labels as described in paper
  
  if(check_labels){
    message(paste0("FleckBrainOrganoidData(): Possible labels are ", paste(possible_labels, collapse = ",")))
    toCheck_labels_src <- setdiff(possible_labels, "seurat_clusters")
    for(i in 1:length(toCheck_labels_src)){
      check_specificityLabels(
        mdta_tbl = dta[[]], grouping_label = "seurat_clusters", toCheck_labels = toCheck_labels_src[i]
      )
    }
  }
  
  return(dta)
  
}

# Braun, ..., Linnarsson, bioRxiv 2022 (https://github.com/linnarsson-lab/developing-human-brain)

BraunBrainData <- function(dataset = file.path("atlas_reference", "braun_data", "human_dev_GRCh38-3.0.0.loom"),
                           braun_donor_sex_mapping_path = file.path("data", "atlas_reference", "braun_data", "donor_sex_mapping.tsv"),
                           possible_labels = c("Clusters", "TopLevelCluster", "CellClass", "Age", 
                                                "Tissue", "Subdivision", "Subregion", "Region", "Sex"),
                           celllabel_nmes = possible_labels, 
                           n_per_celllabel = 500, 
                           seed_sampling = 290, 
                           check_labels = FALSE,
                           do_qc = FALSE, 
                           do_lognormalise = TRUE, logNormCounts_cpu = 1){
  
  dta_path <- file.path("data-raw", dataset)
  dta_raw <- LoomExperiment::import(dta_path, type = "SingleCellLoomExperiment")
  
  # Add sex column based on Xist expression analysis done by authors (see paper)
  braun_donor_sex_mapping <- read_tsv(braun_donor_sex_mapping_path)
  braun_donor_sex_mapping <- setNames(braun_donor_sex_mapping$Sex_manual,
                                      braun_donor_sex_mapping$Donor)
  colData(dta_raw)$Sex <- unname(braun_donor_sex_mapping[colData(dta_raw)$Donor])
  
  # QC: > 0 library size
  
  is_libsize_0 <- dta_raw$TotalUMIs <= 0
  dta <- dta_raw[,!is_libsize_0]
  message(paste0("BraunBrainData(): ", sum(is_libsize_0), " cells with 0 library size removed."))
  
  # Determine specificity and hierarchy of labels as described in paper
  
  if(check_labels){
    
    message(paste0("BraunBrainData(): Possible labels are ", paste(possible_labels, collapse = ",")))
 
    toCheck_labels_src <- setdiff(possible_labels, "Clusters")
    for(i in 1:length(toCheck_labels_src)){
      check_specificityLabels(
        mdta_tbl = colData(dta), grouping_label = "Clusters", toCheck_labels = toCheck_labels_src[i]
      )
    }
    
  }
  
  # Rename assay to counts
  
  counts(dta) <- assay(dta, "matrix")
  assay(dta, "matrix") <- NULL
  
  # QC: Additional (Not done here because accdg to GitHub linnarsson-lab/developing-human-brain, the object (human_dev_GRCh38-3.0.0.loom) should contain high-quality cells from complete processed dataset (HumanFetalBrainPool.h5 in GitHub) except for ~8000 cells with 0 library size
  
  if(do_qc){
    # dta <- addPerCellQC(dta)
    # qc <- quickPerCellQC(colData(dta), 
    #                      percent_subsets = "altexps_ERCC_percent",
    #                      batch = dta$donor,
    #                      subset = dta$donor %in% c("D17", "D7", "D2"))
    # dta <- dta[,!qc$discard]
    message(paste0("BraunBrainData(): Add do_qc code"))
  }
 
  # Change rownames to uniqified gene symbols
  
  gsymbols_uniq <- make.unique(rowData(dta)$Gene)
  rownames(dta) <- gsymbols_uniq
  
  # Remove cells without definite label (based on celllabel_nmes)
  
  if(!all(celllabel_nmes %in% colnames(colData(dta)))){
    stop("BraunBrainData(): Not all celllabel_nmes in cell metadata!"); rm(dta)
  }
  
  #is_cellInvalidLabel <- colData(dta) %>% as_tibble() %>% select(celllabel_nmes) %>% is.na() 
  is_cellInvalidLabel <- colData(dta) %>% as_tibble() %>% select(all_of(celllabel_nmes)) %>% is.na() 
  is_cellInvalidLabel <- rowSums(is_cellInvalidLabel) > 0
  
  dta <- dta[,!is_cellInvalidLabel]
  cell_mdta <- colData(dta) %>% as_tibble()
  
  message(paste0("BraunBrainData(): ", sum(is_cellInvalidLabel), " cells with missing values for any of celllabel_nmes removed."))
  
  # Subsample 

  celllabel_tbl <- cell_mdta[,celllabel_nmes, drop = FALSE]
  celllabel_tbl$united <- celllabel_tbl %>% 
    tidyr::unite(col = "group", all_of(celllabel_nmes), sep = ";")
  
  sampling_dt <- data.table(cell_ind = 1:ncol(dta), celllabel_tbl$united)
  set.seed(seed_sampling)
  sample_inds <- sampling_dt[, .SD[sample(x = .N, size = min(n_per_celllabel, .N), replace = FALSE)], by = group]$cell_ind
  
  dta[["sampling_group"]] <- sampling_dt$group
  dta_subsample <- dta[,sample_inds]
  rm(dta)
  
  # Normalise
  
  if(do_lognormalise){
    dta_subsample <- scuttle::logNormCounts(dta_subsample, BPPARAM = MulticoreParam(logNormCounts_cpu))
  }
  
  return(dta_subsample)
  
}


subsample_data <- function(dta, feature_colnme = "Gene", celllabel_nmes){
  
  # Change rownames to uniqified gene symbols
  
  gsymbols_uniq <- make.unique(rowData(dta)[[feature_colnme]])
  rownames(dta) <- gsymbols_uniq
  
  # Remove cells without definite label (based on celllabel_nmes)
  
  if(!all(celllabel_nmes %in% colnames(colData(dta)))){
    stop("subsample_data(): Not all celllabel_nmes in cell metadata!"); rm(dta)
  }
  
  #is_cellInvalidLabel <- colData(dta) %>% as_tibble() %>% select(celllabel_nmes) %>% is.na()
  is_cellInvalidLabel <- colData(dta) %>% as_tibble() %>% select(all_of(celllabel_nmes)) %>% is.na()
  is_cellInvalidLabel <- rowSums(is_cellInvalidLabel) > 0
  
  dta <- dta[,!is_cellInvalidLabel]
  cell_mdta <- colData(dta) %>% as_tibble()
  
  message(paste0("subsample_data(): ", sum(is_cellInvalidLabel), " cells with missing values for any of celllabel_nmes removed."))
  
  # Subsample 
  
  celllabel_tbl <- cell_mdta[,celllabel_nmes, drop = FALSE]
  celllabel_tbl$united <- celllabel_tbl %>% 
    tidyr::unite(col = "group", all_of(celllabel_nmes), sep = ";")
  
  sampling_dt <- data.table(cell_ind = 1:ncol(dta), celllabel_tbl$united)
  set.seed(seed_sampling)
  sample_inds <- sampling_dt[, .SD[sample(x = .N, size = min(n_per_celllabel, .N), replace = FALSE)], by = group]$cell_ind
  
  dta[["sampling_group"]] <- sampling_dt$group
  dta_subsample <- dta[,sample_inds]
  rm(dta)
  
  # Normalise
  
  if(do_lognormalise){
    dta_subsample <- scuttle::logNormCounts(dta_subsample, BPPARAM = MulticoreParam(logNormCounts_cpu))
  }
  
  return(dta_subsample)
  
}

# Bhaduri...Kriegstein Nature 2020, https://doi.org/10.1038/s41586-020-1962-0
# Dataset is processed and contains normalised counts only

bhaduri_humanprimarycortical_data <- function(
    dataset = file.path("data", "atlas_reference", "bhaduri_data"),
    # Area is not unique per Cluster (recommended label for reference mapping)
    # Each Cluster have unique identity based on each default possible_labels
    chosenlabel_for_mapping = "Cluster",
    possible_labels = c("Class", "State", "Type", "Subtype","Type_Subtype"),
    remove_outlier_cluster = TRUE
    ) {
  
  #library(scRNAseq)
  # BhaduriOrganoidData(ensembl = FALSE) available in scRNAseq package but no cell labels
  
  #####
  
  # Load normalised cell count matrix (189,409 cells)
      
  counts_mat <- fread(file.path(dataset, "exprMatrix.tsv.gz"))
  # counts_mat first column contains genes, remove this column and use genes as row names
  genes <- counts_mat[,1][[1]]
  genes <- gsub(".+[|]", "", genes)
  counts_mat <- counts_mat[,-1]
  rownames(counts_mat) <- genes
  
  cell_meta_df <- read.table(file.path(dataset, "meta.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
  if (any(is.na(cell_meta_df))) {
    stop("bhaduri_humanprimarycortical_data: Missing values in downloaded cell metadata.")
  }
  
  # Create seu then convert to sce object
  # Did this way instead of directly creating sce object because output of function returns index as rownames instead of genes (don't know why), when doing step-by-step all good but when using as function it does not work.
  
  if (identical(rownames(cell_meta_df), colnames(counts_mat))) { # Check that cell barcode order identical
    #seu <- CreateSeuratObject(counts = counts_mat, project = "bhaduri_humanprimarycortical", meta.data = cell_meta_df)
    dta_seu <- CreateSeuratObject(counts = counts_mat, data = counts_mat, project = "bhaduri_humanprimarycortical_data", meta.data = cell_meta_df)
    #dta <- SingleCellExperiment(assays = list(logcounts = counts_mat), colData = cell_meta_df)
    dta <- as.SingleCellExperiment(dta_seu)
    # Only have normalised data, which was added to counts layer for CreateSeuratObject() to work
    counts(dta) <- NULL
  }
  
  if (any(duplicated(rownames(dta)))) {
    stop("Duplicated row/feature names in sce object.")
  }
  
  # Check whether Type-Subtype combination level unique per Cluster
  
  dta$Type_Subtype <- colData(dta) %>% as.data.frame %>% 
    tidyr::unite(col = "Type_Subtype", c("Type", "Subtype"), sep = "_") %>% 
    pull(Type_Subtype)
  
  dta$Area_Cluster <- colData(dta) %>% as.data.frame %>% 
    tidyr::unite(col = "Area_Cluster", c("Area", "Cluster"), sep = "_") %>% 
    pull(Area_Cluster)
  
  check_specificityLabels(mdta_tbl = colData(dta),
                          grouping_label = chosenlabel_for_mapping,
                          toCheck_labels = possible_labels)
  
  ## Extra check, each cluster (rows) should have 1 unique Type_Subtype (columns) identity
  tbl <- table(dta$Cluster, dta$Type_Subtype)
  if (any(rowSums(tbl > 0) > 1)) {
    stop("bhaduri_humanprimarycortical_data: Each cluster has more than 1 Type_Subtype identity.")
  }
  
  # Remove outlier cluster
  if (remove_outlier_cluster) {
    dta <- dta[,!dta$Subtype %in% c("Outlier", "Low Quality", "Microglia low quality")]
  } 
  
  # Number of unique levels per cell label column
  cat("Number of unique levels per cell label column:")
  print(
    lapply(as.data.frame(colData(dta)[,c(chosenlabel_for_mapping, possible_labels, "Area_Cluster")]), FUN = function(col_label){
      length(unique(col_label))
    })
  )
  
  cat("Note that some clusters have the same Type_Subtype identity.")
  
  return(dta)
  
}
  
# heorgan_ref <- HeOrganAtlasData()
# heorgan_ref <- scuttle::logNormCounts(heorgan_ref)
# ref <- DarmanisBrainData()
# lamanno_esc <- LaMannoBrainData('human-es')
# lamanno_midbrain <- LaMannoBrainData('human-embryo')
# lamanno_midbrain <- scuttle::logNormCounts(lamanno_midbrain)
# lamanno_ipsc <- LaMannoBrainData('human-ips')
# nowak_cortex <- NowakowskiCortexData()
# pollen_cortex <- ReprocessedFluidigmData()
# pollen_rglia <- PollenGliaData()
# pollen_rglia <- scuttle::logNormCounts(pollen_rglia)