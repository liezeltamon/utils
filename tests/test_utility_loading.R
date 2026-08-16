# Locate the repository from this test file so the check is independent of the
# directory from which Rscript is invoked.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this check as an Rscript file.")
}
test_file <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_root <- dirname(dirname(test_file))

# Parse every R file first to catch syntax errors without loading dependencies
# or executing exploratory scripts.
utility_files <- list.files(
  repository_root,
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE
)
parse_errors <- vapply(
  utility_files,
  function(path) {
    tryCatch(
      {
        parse(file = path)
        ""
      },
      error = function(condition) conditionMessage(condition)
    )
  },
  character(1)
)
if (any(nzchar(parse_errors))) {
  failed_files <- utility_files[nzchar(parse_errors)]
  failed_paths <- substring(failed_files, nchar(repository_root) + 2L)
  stop(
    "R files failed to parse:\n",
    paste0("- ", failed_paths, ": ", parse_errors[nzchar(parse_errors)], collapse = "\n")
  )
}

# Source fit_lmm from an unrelated working directory and an isolated
# environment to verify that its sibling helpers are located portably.
original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(tempdir())
loading_environment <- new.env(parent = baseenv())
source(
  file.path(repository_root, "linear-modelling", "fit_lmm.R"),
  local = loading_environment
)

# Confirm that the public fitter and representative helper functions were
# loaded together without leaking temporary path-resolution variables.
expected_functions <- c("fit_lmm", "pairwise_contrasts", "prepare_lmm_inputs")
missing_functions <- expected_functions[!vapply(
  expected_functions,
  exists,
  logical(1),
  envir = loading_environment,
  inherits = FALSE,
  mode = "function"
)]
if (length(missing_functions) > 0L) {
  stop("Missing loaded functions: ", paste(missing_functions, collapse = ", "))
}

path_variables <- c("fit_lmm_source_file", "linear_modelling_dir")
leaked_variables <- intersect(
  path_variables,
  ls(loading_environment, all.names = TRUE)
)
if (length(leaked_variables) > 0L) {
  stop("Temporary loading variables leaked: ", paste(leaked_variables, collapse = ", "))
}

message(
  "PASS: parsed ", length(utility_files),
  " R files and loaded fit_lmm with its sibling helpers"
)
