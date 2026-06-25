# 00_annotate_compartment.R
# ---------------------------------------------------------------------------
# Add the `compartment` column (CSF vs PBMC) that extract_pseudobulk.R requires.
#
# The integrated object encodes the sample/tissue in the `sid` field: CSF
# samples contain the substring "CSF" across all four cohorts (Heming,
# Pappalardo, Ramesh, Touil); the remaining samples are treated as PBMC / blood. This step
# reconstructs `compartment` so the pipeline can be reproduced from the raw
# integrated object (`so.GEX.share.Asada.rds`).
#
# Expected per-cohort split (sanity check printed below):
#   osmzhlab_MS_ence_cov     (Heming)     114432 CSF /      0 PBMC
#   PRJNA549712_MS_PBMC_UCSF (Ramesh)      60310 CSF /  96411 PBMC
#   PRJNA671484_MS_Tcell     (Pappalardo)  39269 CSF /  62871 PBMC
#   PRJNA979258_cryoCSF      (Touil)       11823 CSF /      0 PBMC
#   Total CSF cells: 225834  (-> ~221066 after the 8-cell-type filter)
#
# Usage:
#   export QUBOFS_SEURAT_RDS_RAW=/path/to/so.GEX.share.Asada.rds
#   export QUBOFS_SEURAT_RDS=/path/to/so.GEX.share.Asada_with_compartment.rds
#   Rscript 01_pipeline/00_annotate_compartment.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(Seurat))

in_rds  <- Sys.getenv("QUBOFS_SEURAT_RDS_RAW", unset = "")
out_rds <- Sys.getenv("QUBOFS_SEURAT_RDS",     unset = "")
if (in_rds == "" || out_rds == "") {
  stop("Set QUBOFS_SEURAT_RDS_RAW (raw input .rds) and QUBOFS_SEURAT_RDS (annotated output .rds).")
}

so <- readRDS(in_rds)
required_cols <- c("sid", "prj")
missing_cols  <- setdiff(required_cols, colnames(so@meta.data))
if (length(missing_cols) > 0) {
  stop("Required column(s) not found in meta.data: ", paste(missing_cols, collapse = ", "))
}
if ("compartment" %in% colnames(so@meta.data)) {
  message("Existing `compartment` column found; it will be overwritten.")
}

# CSF samples carry "CSF" in `sid`; everything else is treated as PBMC / blood.
so$compartment <- ifelse(grepl("CSF", so@meta.data$sid, ignore.case = TRUE),
                         "CSF", "PBMC")

cat("compartment by cohort:\n")
print(table(so@meta.data$prj, so@meta.data$compartment, useNA = "ifany"))
cat(sprintf("Total CSF cells: %d\n", sum(so$compartment == "CSF")))

out_dir <- dirname(out_rds)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(so, out_rds)
cat("Wrote annotated object: ", out_rds, "\n")
