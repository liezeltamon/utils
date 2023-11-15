message("env.R: Setting working directory to project directory...")

suppressPackageStartupMessages(require(devtools))
set_installPath <- function(path = "package/R"){
  devtools::dev_mode(on = TRUE, path = path)
  .libPaths(path, include.site = TRUE)
  message(paste0("env.R: Setting .libPaths() to: ", paste(.libPaths(), collapse = "  ")))
}
set_installPath("package/R")

###

suppressPackageStartupMessages({
  require(reticulate)
  require(sessioninfo)
  require(ggplot2)
})

options(warnPartialMatchDollar = TRUE)
options(warn = 1)
options(bitmapType = "cairo")

use_condaenv <- function(envname,
                         packages = NULL,
                         channel = c("defaults", "conda-forge", "bioconda"),
                         conda_install_dir = "/var/scratch/exet4759/mamba_installation"){
  
  if(file.exists(envname)){
    reticulate::use_condaenv(envname)
    message(paste0("use_condaenv(): Using ", envname, " environment..."))
  } else {
    reticulate::conda_create(envname = envname, packages = packages, channel = channel,
                             conda = paste0(conda_install_dir, "/conda/condabin/conda")
                             )
    message(paste0("use_condaenv(): Creating ", envname, " environment..."))
  }
 
}

create_dir <- function(path){
  if(!dir.exists(path)){ dir.create(path) }
  return(path)
}

save_sessioninfo <- function(out_dir, out_id){
  sessioninfo::session_info(to_file = paste0(out_dir, "/", out_id, "_session-info.txt"))
  sessioninfo::session_info()
}

theme_set(theme_classic())