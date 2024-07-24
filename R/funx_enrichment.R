################################################################################
# Functional analyses of genes via GO or KEGG (enrichGO and enrichKEGG is more 
# reliable with more updated data than many other tools. - Guangchuang Yu,
# clusterProfiler developer)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(bitr)
# library(DOSE) # for setReadable()
# library(org.Hs.eg.db)
# For help:
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
# Citation:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
#  clusterProfiler: an R package for comparing biological themes among gene clusters. 
#  OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
### FUNCTION ###################################################################
funxAnno <- function(input = list(foreground, background),
                     org_db = org.Hs.eg.db,
                     org = org.Hs.eg.db, # org.Hs.eg.db (GO, human) | "hsa" (KEGG, human)
                     # c("SYMBOL", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                     inputKey = "SYMBOL",
                     useKey = "ENTREZID",
                     approach = "GO_BP", # KEGG, MKEGG, GO_BP, GO_CC, GO_MF, GO_ALL
                     # For saving GO or KEGG table; if NULL, table not saved
                     filePath = NULL,
                     ... # Common parameters between enrichG0() and enrichKEGG()
){
  
  # enrichKEGG() accepts one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
  # Use "ENTREZID" by default so downstream checks applicable and equivalent with when using enrichGO()
  if (approach == "KEGG"){
    message("funxAnno(): Approach is ", approach, ". Using ENTREZID key type by default.")
    useKey <- "ENTREZID"
  }
  
  # Preprocessing
  
  orgdb_keytypes <- columns(org.Hs.eg.db)
  if (!any(c(inputKey, useKey) %in% orgdb_keytypes)) {
    
    message("funxAnno(): Available OrgDb key types: \n")
    cat(orgdb_keytypes)
    stop("funxAnno(): Invalid inputKey and/or useKey")
    
  }
  
  # Clean input gene lists
  input_universe <- keys(org.Hs.eg.db, keytype = inputKey)
  input_cleaned <- lapply(input, FUN = function (genes) {
    unique(genes[genes %in% input_universe])
  })
  
  message("Number of genes after removing duplicates and inputs not in OrgDb: \n")
  print(
    rbind(orig = lengths(input),
          cleaned = lengths(input_cleaned),
          `orig - cleaned` = lengths(input) - lengths(input_cleaned))
  )
  
  # (If required) Convert key type
  if (!is.null(useKey) & !identical(inputKey, useKey)) {
    
    message("funxAnno(): Converting key type")
    conversion_df <- bitr(unique(unlist(input)), fromType = inputKey, toType = useKey, 
                          OrgDb = org.Hs.eg.db, drop = TRUE)
    useInput <- lapply(input, FUN = function (genes) {
      unique(conversion_df[conversion_df[[inputKey]] %in% genes, useKey])
    })
    
    message("Number of genes after conversion: \n")
    print(
      rbind(orig = lengths(input),
            cleaned = lengths(input_cleaned),
            converted = lengths(useInput),
            `cleaned - converted` = lengths(input_cleaned) - lengths(useInput))
    )
      
  } else {
    useKey <- inputKey
    useInput <- input
  }
  
  # ifelse(inputKey == "ncbi-geneid", "ENTREZID", inputKey)
  
  # ORA
  
  if (approach %in% c("KEGG", "MKEGG")) {
    
    if (approach == "KEGG") {
      
      message("funxAnno(): Running enrichKEGG()...")
      e_out <- enrichKEGG(gene = useInput[[1]], organism = org, keyType = "ncbi-geneid",
                          #pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          universe = useInput[[2]],
                          #minGSSize = 10, maxGSSize = 500,
                          #If FALSE, download the latest KEGG data; If TRUE, use KEGG.db
                          use_internal_data = FALSE,
                          ...)
      
    } else if (approach == "MKEGG") {
      
      message("funxAnno(): Running enrichMKEGG()...")
      e_out <- enrichMKEGG(gene = useInput[[1]], organism = org, keyType = "ncbi-geneid",
                           universe = useInput[[2]], ...)
      
    }
    
    e_out <- DOSE::setReadable(e_out, org_db, keyType = "ENTREZID")
    
  } else if (approach %in% c("GO_BP", "GO_CC", "GO_MF", "GO_ALL")) {
    
    appro <- gsub("GO_", "", approach, fixed = TRUE)
    e_out <- enrichGO(gene = useInput[[1]], OrgDb = org_db, keyType = useKey, 
                      ont = appro,
                      # clusterProfiler uses the same adjustment methods 
                      # as the R function stats::p.adjust 
                      # The "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, 
                      # and Yekutieli control the false discovery rate, the expected 
                      # proportion of false discoveries amongst the rejected hypotheses. 
                      # The false discovery rate is a less stringent condition than 
                      # the family-wise error rate, so these methods are more powerful 
                      # than the others.
                      #pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      universe = useInput[[2]],
                      # Gene set size limits for the GO terms to be considered (DEF).
                      # This was added so that parent terms like "biological process"
                      # won't be tested.
                      #minGSSize = 10, maxGSSize = 500,
                      # Whether mapping gene ID to gene Name
                      readable = ifelse(useKey != "SYMBOL", TRUE, FALSE),
                      # If ont='ALL', whether pool 3 GO sub-ontologies
                      pool = ifelse(approach == "GO_ALL", TRUE, FALSE),
                      ...)
    
  } else {
    stop("funxAnnot(): Invalid approach supplied. Choices are GO_BP, GO_CC, GO_MF, GO_ALL and KEGG.
         Can extend code to allow other databases supported by clusterProfiler.")
  }
  
  e_out <- data.frame(e_out, stringsAsFactors = FALSE, row.names = NULL)
  # To avoid separation of description names when opening the csv file.
  e_out$Description <- gsub(",", " ", e_out$Description, fixed = TRUE)
  if (approach != "GO_ALL" & nrow(e_out) > 0) {
    e_out <- cbind.data.frame(ONTOLOGY = approach, e_out, stringsAsFactors = FALSE)
  } else if (approach == "GO_ALL" & nrow(e_out) > 0) {
    e_out$ONTOLOGY <- paste("GO_", e_out$ONTOLOGY, sep = "")
  }
  
  # Save table
  if( !is.null(filePath) ){
    write.csv(e_out, file = filePath, row.names = FALSE, quote = FALSE)
  }
  
  return(e_out)
  
}

# Combine GO_ALL and KEGG enrichment analyses for a set of genes in one 
# dataframe.
#bind_rows(kegg = kegg_out, go = go_out, .id = "source")
