# Functions to call different reference datasets, returning as SCE object

suppressPackageStartupMessages({
  require(tidyverse)
  require(LoomExperiment)
  require(DelayedMatrixStats)
  require(data.table)
  require(scuttle)
  require(BiocParallel)
})

# Braun, ..., Linnarsson, bioRxiv 2022 (https://github.com/linnarsson-lab/developing-human-brain)

check_specificityLabels <- function(mdta_tbl = colData(dta_raw), 
                                    grouping_label = "Clusters",
                                    toCheck_labels = c("TopLevelCluster", "CellClass")){
  tmp_dt <- as.data.table(mdta_tbl[,c(grouping_label, toCheck_labels)])
  tmp_uniq_dt <- tmp_dt[!duplicated(tmp_dt),]
  tmp_uniq_perGroupingLabel <- split(tmp_uniq_dt, f = as.character(tmp_uniq_dt[[grouping_label]]))
  class_len_perGroupingLabel <- unname(unlist(lapply(tmp_uniq_perGroupingLabel, FUN = nrow)))
  if(!all(class_len_perGroupingLabel == 1)){
    warning(paste0(paste(toCheck_labels, collapse = ","), " are NOT unique per ", grouping_label))
  } else{
    message(paste0(paste(toCheck_labels, collapse = ","), " are unique per ", grouping_label))
  }
}
  
BraunBrainData <- function(dataset = "braun_data/human_dev_GRCh38-3.0.0.loom", 
                           celllabel_nmes = possible_labels, 
                           n_per_celllabel = 500, 
                           seed_sampling = 290, 
                           check_labels = FALSE,
                           do_qc = FALSE, 
                           do_lognormalise = TRUE, logNormCounts_cpu = 1){
  
  dta_path <- file.path("data-raw/reference", dataset)
  dta_raw <- import(dta_path, type = "SingleCellLoomExperiment")
  
  braun_donor_sex_mapping <- read_tsv("data/reference/braun_data/donor_sex_mapping.tsv")
  braun_donor_sex_mapping <- setNames(braun_donor_sex_mapping$Sex_manual,
                                      braun_donor_sex_mapping$Donor)
  colData(dta_raw)$Sex <- unname(braun_donor_sex_mapping[colData(dta_raw)$Donor])
  
  possible_labels <- c("Clusters", "TopLevelCluster", "CellClass", "Age", 
                       "Tissue", "Subdivision", "Subregion", "Region", "Sex")
    
  # QC: > 0 library size
  
  is_libsize_0 <- dta_raw$TotalUMIs <= 0
  dta <- dta_raw[,!is_libsize_0]
  message(paste0("BraunBrainData(): ", sum(is_libsize_0), " cells with 0 library size removed."))
  
  # Determine specificity and hierarchy of labels as described in paper
  
  if(check_labels){
    
    message(paste0("BraunBrainData(): Possible labels are ", paste(possible_labels, collapse = ",")))
 
    toCheck_labels_src <- setdiff(possible_labels, "Clusters")
    for(i in 1:length(toCheck_labels_raw)){
      check_specificityLabels(
        mdta_tbl = colData(dta_raw), grouping_label = "Clusters", toCheck_labels = toCheck_labels_src[i]
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
  
  is_cellInvalidLabel <- colData(dta) %>% as_tibble() %>% select(celllabel_nmes) %>% is.na() 
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

#library(scRNAseq)
#bhaduriorganoid_ref <- BhaduriOrganoidData(ensembl = FALSE) # No cell labels
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