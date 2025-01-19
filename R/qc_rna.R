# Bulk RNA qc (Based on https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html)

# Setup

#renv::load()
setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  library(dplyr)

  library(AnnotationHub)
  library(DESeq2)
  library(RNAseqQC)
  
})

# Parameters

src_dir = file.path("results", "process_droplets", "pseudobulk", mainexpname)
counts_path = file.path(src_dir, "counts.csv")
coldata_path = file.path(src_dir, "coldata.csv")
design_formula = ~ Sample + condition

# If NULL, will retrieve latest record for ahub_species
ahub_record = "AH116291" #2023-10-23
ahub_species = "Homo sapiens"
seed_val = 2871

out_dir = create_dir(file.path("results", "process_droplets", "pseudobulk", mainexpname, "qc"))

#-----MAIN-----

# Get annotation hub record to use
if (is.null(ahub_record)) {
  
  ahub_record <- mcols(query(AnnotationHub(), ahub_species)) %>% 
    tibble::as_tibble(rownames = "record_id") %>%
    dplyr::filter(rdataclass == "EnsDb") %>% 
    dplyr::arrange(desc(rdatadateadded)) %>% 
    #dplyr::select(record_id, rdatadateadded) %>% 
    dplyr::slice_head(n = 1) %>% 
    pull(record_id)
  
}

# **Note annotation hub record used**
write.csv()

# Load data

counts_se <- readRDS(counts_se_path)
mainExpName(counts_se) <- mainexpname

count_mat <- read.csv(counts_path, row.names = 1)
meta <- read.csv(coldata_path)

# Create DESeqDataSet

dds <- make_dds(counts = count_mat, 
                metadata = meta, 
                # To add rowData required for gene biotypes plot
                ah_record = ahub_record,
                design = design_formula # The design formula specified in DESeqDataSet()
                )
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = meta,
                              design = design_formula
                              )

                              
pdf(out_dir, height = 5, width = 10)

# Basic metrics - Expect similar trends of these metrics across samples

# Total sample counts
plot_total_counts(dds)

# Library complexity
plot_library_complexity(dds)

# Gene detection
plot_gene_detection(dds)

# Gene biotypes
try(plot_biotypes(dds))

# Filter low-count genes before futher qc

dds <- filter_genes(dds, min_count = 5, min_rep = 2)

# Stabilise variance

vsd <- vst(dds)
mean_sd_plot(vsd)

# Post-filter qc

# Chromosomal expression
try(map(c("1", "5", "14"), ~plot_chromosome(vsd, .x)))

# Replicate variability
colData(vsd)$trt_mut <- paste0(colData(vsd)$condition, "_", colData(vsd)$Sample)
ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

# **Clustering** - important, enable identification of samples that seem to miscluster
set.seed(seed_val)
plot_sample_clustering(vsd, anno_vars = c("Sample", "condition"), distance = "euclidean")

# PCA
plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "condition", shape_by = "Sample")
plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "ENSG00000223972")
pca_res <- plot_pca(vsd, show_plot = FALSE)
plot_loadings(pca_res, PC = 1, annotate_top_n = 5)
plot_pca_scatters(vsd, n_PCs = 5, color_by = "condition", shape_by = "Sample")

# **MDS** - see DeSeq vignette

dev.off()
