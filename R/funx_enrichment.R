################################################################################
# Functional analyses of genes via GO or KEGG (enrichGO and enrichKEGG is more 
# reliable with more updated data than many other tools. - Guangchuang Yu,
# clusterProfiler developer)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(org.Hs.eg.db)
# BiocManager::install("clusterProfiler"); library(clusterProfiler)
# clusterProfiler v3.8.1  For help: 
# https://yulab-smu.github.io/clusterProfiler-book/chapter1.html
# https://guangchuangyu.github.io/software/clusterProfiler
# Citation:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
#  clusterProfiler: an R package for comparing biological themes among gene clusters. 
#  OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
### FUNCTION ###################################################################
funxAnno <- function(input = list(foreground, background),
                     org = org.Hs.eg.db, # org.Hs.eg.db (GO, human) | "hsa" (KEGG, human)
                     # c("SYMBOL", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                     inputKey = "SYMBOL",
                     approach = "GO_BP",
                     # For saving GO or KEGG table; if NULL, table not saved
                     filePath = "lib/ekegg.txt"
){
  
  appro <- gsub(x=approach, pattern="GO_", replacement="", fixed=TRUE)
  accepted.KEGG <- c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")
  accepted.BP <- accepted.CC <- accepted.MF <- accepted.ALL <- c("SYMBOL", "ncbi-geneid") 
  eval(parse(text=paste0(
    "acceptedKeys <- accepted.", appro
  )))
  
  if(!inputKey%in%acceptedKeys){
    stop("Key type not accepted by", approach, ". Choose among ",
         paste(acceptedKeys, collapse=", "), ".")
  }
  
  # Clean gene lists
  input <- lapply(X=input, FUN=function(genes){
    unique(genes[!is.na(genes) & genes!="NA" & genes!="" & genes!=" "])
  })
  
  if(approach=="KEGG"){
    
    e.out <- enrichKEGG(gene=input[[1]], organism=org, keyType=inputKey,
                        # clusterProfiler uses the same adjustment methods 
                        # as the R function stats::p.adjust 
                        # The "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, 
                        # and Yekutieli control the false discovery rate, the expected 
                        # proportion of false discoveries amongst the rejected hypotheses. 
                        # The false discovery rate is a less stringent condition than 
                        # the family-wise error rate, so these methods are more powerful 
                        # than the others.
                        pvalueCutoff=0.05, pAdjustMethod="BH", universe=input[[2]],
                        # Gene set size limits for the GO terms to be considered (DEF).
                        # This was added so that parent terms like "biological process"
                        # won't be tested. 
                        minGSSize=10, maxGSSize=500, qvalueCutoff=0.2,
                        #If FALSE, download the latest KEGG data; If TRUE, use KEGG.db
                        use_internal_data=FALSE)
    
  } else if( approach%in%c("GO_BP", "GO_CC", "GO_MF", "GO_ALL") ){
    
    e.out <- enrichGO(gene=input[[1]], OrgDb=org, 
                      keyType=ifelse(inputKey=="ncbi-geneid", "ENTREZID", inputKey), 
                      ont=appro, pvalueCutoff=0.05, pAdjustMethod="BH", 
                      universe=input[[2]], minGSSize=10, maxGSSize=500,
                      # Whether mapping gene ID to gene Name
                      readable=ifelse(inputKey=="ncbi-geneid", TRUE, FALSE),
                      # If ont='ALL', whether pool 3 GO sub-ontologies
                      pool=ifelse(approach=="GO_ALL", TRUE, FALSE))
    
  } else {
    stop("Invalid approach supplied. Choices are GO_BP, GO_CC, GO_MF, GO_ALL and KEGG.")
  }
  
  e.out <- data.frame(e.out, stringsAsFactors=FALSE, row.names=NULL)
  # To avoid separation of description names when opening the csv file.
  e.out$Description <- gsub(x=e.out$Description, pattern=",", replacement=" ",  
                            fixed=TRUE)
  if(approach!="GO_ALL" & nrow(e.out)>0){
    e.out <- cbind.data.frame(ONTOLOGY=approach, e.out, stringsAsFactors=FALSE)
  } else if(approach=="GO_ALL" & nrow(e.out)>0){
    e.out$ONTOLOGY <- paste("GO_", e.out$ONTOLOGY, sep="")
  }
  
  # Save table
  if( !is.null(filePath) ){
    write.csv(x=e.out, file=filePath, row.names=FALSE, quote=FALSE)
  }
  
  return(e.out)
  
}
################################################################################

################################################################################
# Combine GO_ALL and KEGG enrichment analyses for a set of genes in one 
# dataframe.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/funxAnno.R"))
# library(org.Hs.eg.db)
# library(clusterProfiler)
### FUNCTION ###################################################################
funxAnnoWrapper <- function(input = 'list(HUGO foreground genes, HUGO background genes)',
                            filePath = 'dir for saving individual GO and KEGG tables; 
                                        if NULL, tables not saved',
                            hugoEntrez.file = 'Conversion table of HUGO gene 
                                               symbols to ncbi-geneid for KEGG analysis'
){
  
  # Read conversion table for ncbi-geneid to HUGO symbols for output
  HE.conv <- read.delim(file=hugoEntrez.file, header=T, row.names=NULL, stringsAsFactors=F)
  HE.conv <- HE.conv[!is.na(HE.conv$SYMBOL) & !is.na(HE.conv$ENTREZID),]
  HE.conv$ENTREZID <- as.character(HE.conv$ENTREZID)
  
  # Gene sets
  GENE <- list()
  
  GENE$GO_ALL <- input
  GENE$GO_ALL <- lapply(X=GENE$GO_ALL, FUN=function(set){
    return( set[!is.na(set)] )
  }) 
  
  GENE$KEGG <- lapply(X=GENE$GO_ALL, FUN=function(set){
    set <- HE.conv$ENTREZID[ HE.conv$SYMBOL%in%set ]
    return( set[!is.na(set)] )
  })
  
  # Enrichment analysis
  funxAnnoOut <- sapply(X=names(GENE), simplify=F, FUN=function(app){
    
    chunk <- funxAnno(input=GENE[[app]],
                      org=ifelse(app=="GO_ALL", "org.Hs.eg.db", "hsa"), 
                      inputKey=ifelse(app=="GO_ALL", "SYMBOL", "ncbi-geneid"), 
                      approach=app, filePath=filePath)
    return(chunk)
    
  })
  
  funxAnnoOut <- do.call("rbind", funxAnnoOut)
  rownames(funxAnnoOut) <- NULL
  
  rm(GENE)
  
  # Convert ncbi geneids to HUGO
  KEGG.TF <- funxAnnoOut$ONTOLOGY=="KEGG"
  if( sum(KEGG.TF)>0 ){
    funxAnnoOut[KEGG.TF,"geneID"] <- sapply(X=funxAnnoOut[KEGG.TF, "geneID"], 
                                            simplify=TRUE, FUN=function(id){
                                              id <- unique(strsplit(x=id, split="/", fixed=TRUE)[[1]])
                                              paste(x=HE.conv[HE.conv$ENTREZID%in%id,"SYMBOL"],
                                                    collapse="/")
                                            })
  }
  
  return(funxAnnoOut)
  
} 
################################################################################
