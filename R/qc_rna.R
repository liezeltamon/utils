# Bulk RNA qc (Based on https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html)

# Setup

renv::load()
setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  library(dplyr)
  library(assertthat)
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  
  library(AnnotationHub)
  library(DESeq2)
  library(RNAseqQC)
  library(glmpca)
  
})

# Parameters

seed_val = 2871

counts_path = file.path("results", "create_counts_matrix", "counts.csv")
coldata_path = file.path("results", "create_metadata", "metadata.csv")

# If NULL, will retrieve latest record for ahub_species
use_ahub_record = TRUE
ahub_record = "AH116291" #2023-10-23
ahub_species = "Homo sapiens"

# (Optional) Subset samples
# If subset_vars = NULL, no subsetting. Columns united if more than 1 input
# subset_vars = c("purpose", "treatment")
# subset_vars_lvls = c("compare_condition_01untreated",
#                      "compare_condition_treatment_01untreated")

subset_vars = "purpose"
subset_vars_lvls = "compare_condition" #"compare_condition_treatment" #"compare_condition"

# (Optional) Drop samples
drop_samples = NULL #c("ABA", "JA", "SMA")

# Add metadata column
combine_metacols = NULL

anno_vars = c("condition", "age_group") #c("condition", "treatment")
design_formula = ~ inferred_sex + age_group + condition #~ inferred_sex + age_group + condition_treatment
# Minimum number of samples in a group to be used for comparison, used by filter_genes()
min_n_sample_group = 5 #3 #5

out_dir = create_dir(file.path(
  "results", "qc", "qc_rna",
  paste0(paste(subset_vars, collapse = "_"), 
         "_", 
         paste(subset_vars_lvls, collapse = "_"),
         "_drop_samples_",
         paste(drop_samples, collapse = "_")
  )
))

save.image(file.path(out_dir, "parameters.RData"))

# Plotting parameters

plot_height = 12
plot_width = 16

# Scatter plots
point_rel_size = 2.5
  
# plot_pca()
point_label_var = "inferred_sex" # person_id
colour_var = anno_vars[1]
shape_var = anno_vars[2]

#----- MAIN -----

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
#write.csv()

# Load data

# **If SummarizedExperiment** class
# counts_se <- readRDS(counts_path)
# mainExpName(counts_se) <- mainexpname

counts_mat <- read.csv(counts_path, row.names = 1, check.names = FALSE)
# meta <- read_csv(coldata_path, na = "") %>% 
#   as.data.frame() %>% 
#   column_to_rownames("...1")
meta <- read.delim(coldata_path, row.names = 1, check.names = FALSE, na.strings = "")
  
# (Optional) Drop samples
if (!is.null(drop_samples)) {
  counts_mat <- counts_mat[, setdiff(colnames(counts_mat), drop_samples)]
  meta <- meta[setdiff(colnames(counts_mat), drop_samples), ]
}

# (Optional) Add new metadata columns by combining existing columns
if (!is.null(combine_metacols) & length(combine_metacols) > 1) {
  
  meta <- meta %>% 
    tidyr::unite(!!paste(combine_metacols, collapse = "_"),
                 !!!syms(combine_metacols),
                 remove = FALSE)

}

# (Optional) Subset samples
if (!is.null(subset_vars)) {
  
  if (length(subset_vars) > 1) {
    meta <- meta %>% 
      tidyr::unite(!!paste(subset_vars, collapse = "_"),
                   !!!syms(subset_vars),
                   remove = FALSE)
    subset_vars <- paste(subset_vars, collapse = "_")
  }
  
  is_subset <- meta[[subset_vars]] %in% subset_vars_lvls
  counts_mat <- counts_mat[ ,is_subset]
  meta <- meta[is_subset, ]
  
}

meta <- meta %>% 
  dplyr::mutate_if(is.character, as.factor)

assert_that(identical(colnames(counts_mat), rownames(meta)),
            msg = "Column names of counts matrix not identical to row names of column metadata")

# Create DESeqDataSet

if (use_ahub_record) {
  
  # Needs internet connection
  dds_ahub <- make_dds(counts = counts_mat,
                       metadata = meta,
                       # To add rowData required for gene biotypes plot
                       ah_record = ahub_record,
                       design = design_formula # The design formula specified in DESeqDataSet()
                  
  )
  
}

dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = meta,
                              design = design_formula)

saveRDS(dds, file.path(out_dir, "dds.rds"))

pdf(file.path(out_dir, "plots.pdf"), height = plot_height, width = plot_width)

# Basic metrics - Expect similar trends of these metrics across samples

# Total sample counts
plot_total_counts(dds) + 
  ggtitle("Library size: Good if samples in the same order of magnitude")

# Library complexity
plot_library_complexity(dds) + 
  ggtitle("Library complexity: Good if consistent across samples")

# Gene detection
plot_gene_detection(dds) + 
  ggtitle("Gene detection: Good if consistent across samples")

# Gene biotypes
try(plot_biotypes(dds_ahub))

# Filter low-count genes before futher qc

# "A good strategy is to determine the size of a biological group, 
# i.e. all samples under the same biological condition, and filter (keep) genes 
# with a certain count at least as often as the size of the biological group. 
# In our case, each group is a combination of treatment and mutation status with 
# 4 samples and we choose 5 counts as threshold." (RNASeqQC vignette)
dds <- filter_genes(dds, min_count = min_n_sample_group + 1, min_rep = min_n_sample_group)

# Stabilise variance

vsd <- vst(dds)
mean_sd_plot(vsd) +
  ggtitle("Variance stabilisation: red trend line should be relatively flat with respect to the scale on the y-axis")

# Post-filter qc

# Chromosomal expression
#try(map(c("1", "5", "14"), ~plot_chromosome(vsd, .x)))

# Replicate variability
# colData(vsd)$trt_mut <- paste0(colData(vsd)$condition, "_", colData(vsd)$Sample)
# ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
# cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

# For human samples, extrapolate or confirm sex (if available already)
# data.frame with chromosome information for each gene
data("geneInfo", package = "DEGreport")
DEGreport::degSignature(counts(dds),
                        signature = geneInfo,
                        group = "person_id",
                        metadata = as.data.frame(colData(dds)),
                        # Defaults to note
                        scale = FALSE) +
  labs(y = "median counts(dds)") +
  ggtitle("Check presence or absence of reads, variation across samples expected because no library size normalisation yet")

DEGreport::degSignature(assay(vsd),
                        signature = geneInfo,
                        group = "person_id",
                        metadata = as.data.frame(colData(dds)),
                        # Defaults to note
                        scale = FALSE) +
  labs(y = "median vst") +
  ggtitle("Used vst")

### Correlation of covariates

# Covariates correlation with PCs
dds <- DESeq2::estimateSizeFactors(dds)
try({
  DEGreport::degCovariates(log2(counts(dds) + 0.5), 
                           colData(dds),
                           # Defaults to note
                           fdr = 0.1,
                           correlation = "kendall")
  
  title(main = "Used log2(counts(dds) + 0.5); kendall cor")
})

# Transformation attenuates influence of library size
try({
  DEGreport::degCovariates(assay(vsd),
                           colData(dds),
                           # Defaults to note
                           fdr = 0.1,
                           correlation = "kendall")
  title(main = "Used vst; kendall cor")
})

# Covariates correlation with metrics
DEGreport::degCorCov(as.data.frame(colData(dds)))

### PCA

# Replicates or samples expected to be similar cluster together?
plot_pca(vsd,
         PC_x = 1,
         PC_y = 2,
         n_feats = 500, # Default
         scale_feats = TRUE, # Default is FALSE
         color_by = anno_vars[1],
         shape_by = anno_vars[2],
         point_rel_size = point_rel_size,
         show_plot = TRUE
)

# Same plot but points are labelled
pca_res <- plot_pca(vsd,
                    PC_x = 1,
                    PC_y = 2,
                    n_feats = 500, # Default
                    scale_feats = TRUE, # Default is FALSE
                    color_by = anno_vars[1],
                    shape_by = anno_vars[2],
                    point_rel_size = point_rel_size,
                    show_plot = FALSE
)

p_pca <- ggplot(pca_res$data, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .data[[colour_var]],
                 shape = .data[[shape_var]]),
             size = point_rel_size) +
  geom_text_repel(aes(label = .data[[point_label_var]])) +
  labs(
    x = paste0("PC1 (", round(pca_res$var_exp[1], digits = 4), "%)"),
    y = paste0("PC2 (", round(pca_res$var_exp[2], digits = 4), "%)")
  ) +
  theme_classic()

plot_pca_scatters(vsd,
                  n_PCs = min(4, ncol(dds)),
                  n_feats = 500, # Default
                  scale_feats = TRUE, # Default is FALSE
                  color_by = anno_vars[1],
                  shape_by = anno_vars[2],
                  point_rel_size = 3
)

### **Clustering** - important, enable identification of samples that seem to miscluster

set.seed(seed_val)
plot_sample_clustering(vsd, 
                       anno_vars = anno_vars,
                       distance = "euclidean")

set.seed(seed_val)
plot_sample_clustering(vsd, 
                       anno_vars = anno_vars,
                       distance = "pearson")

set.seed(seed_val)
plot_sample_clustering(vsd, 
                       anno_vars = anno_vars,
                       distance = "spearman")

### **MDS** - see DeSeq vignette

# Modified from https://github.com/liezeltamon/neuron-microglia_coculture/blob/0e21f0555e9996b05f8b5f6ce6c17db0307c553c/code/merge_de/compare_samples.R#L189
sample_dists <- dist(t(assay(vsd)))
sample_dist_mx <- as.matrix(sample_dists)
n_dim = min(3, (ncol(vsd) - 1))

mds_mx <- cmdscale(sample_dist_mx, k = n_dim)
colnames(mds_mx) <- paste0("MDS", 1:ncol(mds_mx))
mds_df <- cbind.data.frame(as.data.frame(colData(vsd)), mds_mx)

# First two dimensions
p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = .data[[colour_var]],
                 shape = .data[[shape_var]]),
             size = point_rel_size) +
  geom_text_repel(aes(label = .data[[point_label_var]])) +
  theme_classic()

ggpairs(mds_df, columns = paste0("MDS", as.character(1:n_dim)),
        mapping = ggplot2::aes(color = .data[[anno_vars[1]]]),
        upper = list(continuous = wrap("points", size = point_rel_size)),
        lower = list(continuous = wrap("points", size = point_rel_size)), legend = 2) +
  #scale_shape_manual(values = shape_ints) +
  ggtitle("MDS with VST data")

### Generalised PCA for non-normal data (refer to https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)

#gpca <- glmpca(counts(dds), L = 2)
gpca <- glmpca(assay(vsd), L = 2)
gpca_df <- gpca$factors
gpca_df <- cbind(
  gpca_df,
  colData(dds)[, c(anno_vars, point_label_var)]
)
p_glmpca <- ggplot(gpca_df, aes(x = dim1, y = dim2)) +
  geom_point(aes(colour = .data[[colour_var]],
                 shape = .data[[shape_var]]),
             size = point_rel_size) +
  labs(x = "GLMPCA1", y = "GLMPCA2") +
  geom_text_repel(aes(label = .data[[point_label_var]])) +
  theme_classic()

### MDS + PCA

cowplot::plot_grid(p_pca, p_mds, p_glmpca, ncol = 2)

dev.off()

# rm(list = ls()); gc()
