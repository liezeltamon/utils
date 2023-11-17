require(R.utils)

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
cellranger_renamefastq <- function(src_dir, out_dir, samples){
  for(samp in samples){
    
    samp_paths <- list.files(src_dir, pattern = paste0("^", samp, "_R"))
    samp_paths <- samp_paths[grepl("_001.fastq.gz$", samp_paths)]
    
    for(pth in samp_paths){
      pth_renamed <- gsub(samp, paste0(samp, "_S1_L001"), pth)
      R.utils::createLink(link = file.path(out_dir, pth_renamed),
                          target = file.path(src_dir, pth))
    }
    
  }
}

#rm(list=ls())
