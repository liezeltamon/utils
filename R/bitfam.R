library(BITFAM)

run_bitfam <- function(counts_mat, species, ncpu, out_dir, convergence_tolerance = 0.005) {
  
  # Using BITFAM_preprocess() will select 5000 variable genes and hence subset the normalised counts matrix
  norm_mx <- BITFAM_preprocess(counts_mat)
  
  res <- BITFAM(data = norm_mx, species = species, interseted_TF = NA, scATAC_obj = NA, ncores = ncpu, tol_rel_obj = convergence_tolerance)
  Z <- BITFAM_activities(res)
  save(norm_mx, res, Z, file = file.path(out_dir, "BITFAM_run.RData"))
  write.csv(Z, file = file.path(out_dir, "BITFAM_activities.csv"))
  return(Z)
  
}