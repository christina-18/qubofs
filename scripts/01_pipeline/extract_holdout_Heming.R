# Run extract_pseudobulk.R with the Heming cohort (osmzhlab_MS_ence_cov)
# as the external hold-out test cohort.
#
# Usage from the repository root:
#   Rscript scripts/01_pipeline/extract_holdout_Heming.R
#
# Usage from R within the scripts/01_pipeline/ directory:
#   source("extract_holdout_Heming.R")
#
# The variable `holdout_cohort_arg` is read by extract_pseudobulk.R when sourced.
holdout_cohort_arg <- "osmzhlab_MS_ence_cov"

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  script_dir <- getwd()
}

source(file.path(script_dir, "extract_pseudobulk.R"))
