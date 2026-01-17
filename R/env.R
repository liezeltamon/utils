#message("env.R: Setting working directory to project directory...")
#suppressPackageStartupMessages(library(rprojroot))
#rprojroot::find_rstudio_root_file() # Has to be in r setup on markdown

suppressPackageStartupMessages(library(devtools))
set_installPath <- function(path = file.path("package", "R")){
  devtools::dev_mode(on = TRUE, path = path)
  .libPaths(path, include.site = TRUE)
  message(paste0("env.R: Setting .libPaths() to: ", paste(.libPaths(), collapse = "  ")))
}
#set_installPath("package/R")

###

# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(svglite) # For saving SVG with ggplot2
#   library(reticulate)
#   library(sessioninfo)
#   library(foreach)
#   library(doParallel)
#   library(itertools)
# })

options(warnPartialMatchDollar = TRUE)
options(warn = 1)
options(bitmapType = "cairo")

# library(reticulate)
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
  if(!dir.exists(path)){ dir.create(path, recursive = TRUE) }
  return(path)
}

#library(sessioninfo)
save_sessioninfo <- function(out_dir, out_id){
  sessioninfo::session_info(to_file = paste0(out_dir, "/", out_id, "_session-info.txt"))
  sessioninfo::session_info()
}

#library(ggplot2)

theme_default <- function(base_size = 12, base_family = "sans") { # "DejaVu Sans"
  theme_linedraw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(), #element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_blank(), #element_line(color = "grey90", size = 0.25),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      text = element_text(family = base_family),
      axis.title   = element_text(family = base_family),
      axis.text    = element_text(family = base_family),
      legend.title = element_text(family = base_family),
      legend.text  = element_text(family = base_family),
      strip.text   = element_text(family = base_family, colour = "black"),
      strip.background = element_blank()
      )
}

theme_set(theme_default())

##### Stash for reference

#library(foreach)
#library(doParallel)
#library(itertools)
setup_parallel <- function(n_cpu) {
  if(n_cpu > 1){
    registerDoParallel(cores = n_cpu)
    return(`%dopar%`)
    print(paste0("Running with ", n_cpu, " cores."), quote = FALSE)
  } else {
    return(`%do%`)
  }
}
