library(withr)

library(reticulate)

conda_path <- "/var/scratch/exet4759/mamba_installation/conda/condabin/conda"
reticulate::conda_create("./envs",
                         conda = conda_path)
reticulate::use_condaenv(condaenv = "reticulate-env", conda = conda_path, required = FALSE)

reticulate::import(module = 'leidenalg')
library(Seurat)
FindClusters(example_seurat, algorithm='leiden')