require(data.table)
require(RColorBrewer)
require(R.utils)
require(tidyverse)

##### Stash for reference
# colData(do.call("cbind", sce_list)) %>% 
#   as_tibble() %>%
#   select(Sample, sum) %>%
#   distinct() %>% 
#   group_by(Sample) %>%
#   summarise_all(list(max = max, median = median, 
#                      q0.75 = function(x) quantile(x, probs = 0.75), 
#                      q0.99 = function(x) quantile(x, probs = 0.99)))

# SCE gene ID -> unique gene symbols
# gsymbols_uniq <- make.unique(rowData(sce)$Symbol)
# rownames(sce) <- gsymbols_uniq

extendPaletteFUN <- colorRampPalette(brewer.pal(9, "Set1"))
#####

##### Make generic
check_specificityLabels <- function(mdta_tbl = colData(dta_raw), 
                                    grouping_label = "Clusters",
                                    toCheck_labels = c("TopLevelCluster", "CellClass"),
                                    return_uniq_dt = FALSE){
  tmp_dt <- as.data.table(mdta_tbl[,c(grouping_label, toCheck_labels)])
  tmp_uniq_dt <- tmp_dt[!duplicated(tmp_dt),]
  tmp_uniq_perGroupingLabel <- split(tmp_uniq_dt, f = as.character(tmp_uniq_dt[[grouping_label]]))
  class_len_perGroupingLabel <- unname(unlist(lapply(tmp_uniq_perGroupingLabel, FUN = nrow)))
  if(!all(class_len_perGroupingLabel == 1)){
    warning(paste0(paste(toCheck_labels, collapse = ","), " are NOT unique per ", grouping_label))
  } else{
    message(paste0(paste(toCheck_labels, collapse = ","), " are unique per ", grouping_label))
  }
  
  if(return_uniq_dt){
    return(tmp_uniq_dt)
  }
}
#####

#####
# No fix for levels of same size
as_factor_bySize <- function(x, decreasing = TRUE){
  x <- as.character(x)
  ordered_levels <- sort(table(x), decreasing = decreasing)
  ordered_levels <- names(ordered_levels)
  x <- factor(x, levels = ordered_levels)
  return(x)
}

# Rename fastq for cellranger count and friends (e.g. SRR14710616_S1_L001_R2_001.fastq.gz)
## Illumina naming convention - https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
## - The lane number in the function is set to 1 (L001) as this is not often shown in fastq original file name 
# src_dir = file.path("data-raw", "2023-10-06_scRNAseq")
# out_dir = create_dir(file.path("data", "2023-10-06_scRNAseq_renamed"))
# samp = "CO-NSC1"
cellranger_renamefastq <- function(src_dir, out_dir, samples, patt = paste0("^", samp, "_R.+_001.fastq.gz$"), rename = TRUE){
  for(samp in samples){
    
    samp_paths <- list.files(src_dir, pattern = patt)
    samp_paths <- samp_paths[grepl("_001.fastq.gz$", samp_paths)]
    
    for(pth in samp_paths){
      if(rename){
        pth_renamed <- gsub(samp, paste0(samp, "_S1_L001"), pth)
      } else {
        pth_renamed <- pth
      }
      R.utils::createLink(link = file.path(out_dir, pth_renamed),
                          target = file.path(src_dir, pth))
    }
    
  }
}
#rm(list=ls())