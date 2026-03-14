# R/paths.R
# Centralized path definitions for the project.
# Detects project root automatically so scripts run correctly from any
# working directory (project root or a subdirectory like glmm_analysis/).

.find_project_root <- function() {
  path <- normalizePath(getwd(), winslash = "/")
  for (i in seq_len(10)) {
    if (dir.exists(file.path(path, "R")) && dir.exists(file.path(path, "data"))) {
      return(path)
    }
    parent <- dirname(path)
    if (parent == path) break  # reached filesystem root
    path <- parent
  }
  stop("Could not find project root (expected R/ and data/ folders).")
}

.PROJECT_ROOT <- .find_project_root()

DATA_DIR   <- file.path(.PROJECT_ROOT, "data")
OUTPUT_DIR <- file.path(.PROJECT_ROOT, "output")
PLOTS_DIR  <- file.path(OUTPUT_DIR,   "plots")

# Helper functions to build full paths
data_path   <- function(...) file.path(DATA_DIR,   ...)
output_path <- function(...) file.path(OUTPUT_DIR, ...)
plot_path   <- function(...) file.path(PLOTS_DIR,  ...)

# Ensure output folders exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
if (!dir.exists(PLOTS_DIR))  dir.create(PLOTS_DIR, recursive = TRUE)