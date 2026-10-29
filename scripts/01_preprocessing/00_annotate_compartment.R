# 00_annotate_compartment.R

# Assign CSF or PBMC labels based on sample IDs (sid).

# Samples containing "CSF" in sid are classified as CSF; all others as PBMC.
# The resulting compartment column is used by the downstream preprocessing scripts.
# Expected cell counts by cohort:
# Heming       114432 CSF /     0 PBMC
# Ramesh        60310 CSF / 96411 PBMC
# Pappalardo    39269 CSF / 62871 PBMC
# Touil         11823 CSF /     0 PBMC
# Total: 225834 CSF cells (221066 after filtering to eight immune cell types).

# Input:  so.GEX.share.Asada.rds
# Output: Seurat object with a compartment column
# Run: Rscript scripts/01_preprocessing/00_annotate_compartment.R

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
