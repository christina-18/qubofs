# Run extract_pseudobulk.R with the Ramesh cohort
# (PRJNA549712_MS_PBMC_UCSF) as the external hold-out test cohort.
#
# Note: the cohort identifier follows the integrated Seurat object metadata.
# The identifier string contains "PBMC", but the manuscript analysis and this
# extraction are restricted to the CSF compartment.
#
# Usage from the repository root:
#   Rscript scripts/01_pipeline/extract_holdout_Ramesh.R
#
# Usage from R within the scripts/01_pipeline/ directory:
#   source("extract_holdout_Ramesh.R")
#
# The variable `holdout_cohort_arg` is read by extract_pseudobulk.R when sourced.
holdout_cohort_arg <- "PRJNA549712_MS_PBMC_UCSF"

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  script_dir <- getwd()
}

source(file.path(script_dir, "extract_pseudobulk.R"))
