#!/usr/bin/env Rscript
# Audit Azimuth l3 subtype cell counts per donor.
#
# Goal: assess feasibility of finer-grained pseudobulk for CD4_T / CD8_T / DC
# in CSF, by counting how many donors have >= 20 cells in each l3 subtype.
#
# Output:
#   /Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/results/l3_audit/
#       l3_cell_counts_per_donor.csv
#       l3_donors_passing_20.csv
#       l3_audit_summary.txt
#
# Run from Mac (R + Seurat installed). Expected runtime ~15-30 min depending
# on whether l3 annotation already exists on the object.
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

INPUT_RDS  <- "/Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/so.GEX.share.Asada.rds"
OUT_DIR    <- "/Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/results/l3_audit"
DONOR_COL  <- "patient_id"   # adjust if needed; alternatives: "sid", "donor_id"
# TISSUE_COL is auto-detected below by scanning columns for "CSF" / "PBMC" values

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("Azimuth l3 cell-count audit\n")
cat("==============================================================\n\n")
cat("Loading:", INPUT_RDS, "\n")
so <- readRDS(INPUT_RDS)
cat("Object loaded. Class:", class(so), "\n")
cat("Cells:", ncol(so), "  Genes:", nrow(so), "\n\n")

md <- so@meta.data
cat("Metadata columns available:\n")
print(colnames(md))

# --- Identify donor / tissue / l2 / l3 columns ---
stopifnot(DONOR_COL %in% colnames(md))

# Auto-detect tissue column: scan all character/factor columns for CSF / PB[MC] values
TISSUE_COL <- NA_character_
for (cn in colnames(md)) {
  vals <- unique(as.character(md[[cn]]))
  has_csf  <- any(c("CSF","csf") %in% vals)
  has_pbmc <- any(c("PBMC","pbmc","PB","pb","Blood","blood") %in% vals)
  if (has_csf && has_pbmc) {
    TISSUE_COL <- cn
    cat("\nAuto-detected tissue column:", cn, "\n")
    cat("  values:", paste(head(vals, 8), collapse=", "), "\n")
    break
  }
}
# Normalize PB -> PBMC for downstream consistency
if (!is.na(TISSUE_COL)) {
  v <- as.character(md[[TISSUE_COL]])
  v[v %in% c("PB","pb","Blood","blood")] <- "PBMC"
  v[v %in% c("csf")] <- "CSF"
  md[[TISSUE_COL]] <- v
  so[[TISSUE_COL]] <- v
  cat("  normalized values: CSF =", sum(v=="CSF"), ", PBMC =", sum(v=="PBMC"), "\n")
}
if (is.na(TISSUE_COL)) {
  # Try to derive from "sid" / "prj" via regex
  cat("\n*** No direct CSF/PBMC column found. Trying to derive from 'sid' or 'prj' ***\n")
  if ("sid" %in% colnames(md)) {
    derived <- ifelse(grepl("CSF|csf", md$sid), "CSF",
              ifelse(grepl("PBMC|pbmc|Blood|blood", md$sid), "PBMC", NA))
    cat("From sid: derived", sum(!is.na(derived)), "/", length(derived), "cells\n")
    if (sum(!is.na(derived)) / length(derived) > 0.8) {
      md$tissue_derived <- derived
      so$tissue_derived <- derived
      TISSUE_COL <- "tissue_derived"
      cat("Using derived 'tissue_derived' column.\n")
    }
  }
  if (is.na(TISSUE_COL) && "prj" %in% colnames(md)) {
    derived <- ifelse(grepl("CSF|csf", md$prj), "CSF",
              ifelse(grepl("PBMC|pbmc|Blood|blood|Tcell", md$prj), "PBMC", NA))
    cat("From prj: derived", sum(!is.na(derived)), "/", length(derived), "cells\n")
    if (sum(!is.na(derived)) / length(derived) > 0.8) {
      md$tissue_derived <- derived
      so$tissue_derived <- derived
      TISSUE_COL <- "tissue_derived"
      cat("Using derived 'tissue_derived' column from prj.\n")
    }
  }
}
if (is.na(TISSUE_COL)) {
  cat("\nFAILED to identify tissue column. Available columns and their unique values (first 8):\n")
  for (cn in colnames(md)) {
    vals <- unique(as.character(md[[cn]]))
    if (length(vals) <= 20) {
      cat(sprintf("  %s: %s\n", cn, paste(head(vals, 8), collapse=", ")))
    }
  }
  stop("Please set TISSUE_COL manually based on the column listing above.")
}
cat("Using TISSUE_COL =", TISSUE_COL, "\n\n")

# Look for l3-style annotation column
l3_candidates <- c("predicted.celltype.l3", "celltype_l3", "azimuth_l3", "celltype.l3")
l3_col <- l3_candidates[l3_candidates %in% colnames(md)][1]
l2_candidates <- c("predicted.celltype.l2", "celltype_l2", "azimuth_l2", "celltype.l2", "celltype")
l2_col <- l2_candidates[l2_candidates %in% colnames(md)][1]

cat("\nDetected l2 column:", l2_col, "\n")
cat("Detected l3 column:", l3_col, "\n")

if (is.na(l3_col)) {
  cat("\n*** l3 annotation NOT found. Running Azimuth l3 annotation now. ***\n")
  cat("Note: this requires SeuratData and Azimuth packages.\n")
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Please install Azimuth: remotes::install_github('satijalab/azimuth')")
  }
  library(Azimuth)
  so <- RunAzimuth(so, reference = "pbmcref", annotation.levels = c("celltype.l2", "celltype.l3"))
  md <- so@meta.data
  l3_col <- "predicted.celltype.l3"
  l2_col <- "predicted.celltype.l2"
} else {
  cat("\n*** l3 annotation already present. Skipping Azimuth re-run. ***\n")
}

# --- Per-donor / per-l3 cell counts ---
cat("\n--- Building (donor x l3-subtype x tissue) cell count table ---\n")
counts <- md %>%
  count(.data[[DONOR_COL]], .data[[TISSUE_COL]], .data[[l2_col]], .data[[l3_col]],
        name = "n_cells") %>%
  rename(donor_id = !!DONOR_COL, tissue = !!TISSUE_COL,
         l2 = !!l2_col, l3 = !!l3_col)

write_csv(counts, file.path(OUT_DIR, "l3_cell_counts_per_donor.csv"))
cat("Wrote:", file.path(OUT_DIR, "l3_cell_counts_per_donor.csv"), "\n")

# --- Summary: how many donors have >= 20 cells per (l3, tissue) ---
THRESH <- 20
summary_tab <- counts %>%
  group_by(tissue, l2, l3) %>%
  summarise(
    n_donors_total   = n(),
    n_donors_pass    = sum(n_cells >= THRESH),
    median_cells     = median(n_cells),
    min_cells        = min(n_cells),
    max_cells        = max(n_cells),
    .groups = "drop"
  ) %>%
  arrange(tissue, l2, desc(n_donors_pass))

write_csv(summary_tab, file.path(OUT_DIR, "l3_donors_passing_20.csv"))
cat("Wrote:", file.path(OUT_DIR, "l3_donors_passing_20.csv"), "\n")

# --- Plain-text summary, focused on CD4_T / CD8_T / DC ---
sink(file.path(OUT_DIR, "l3_audit_summary.txt"))
cat("================================================================\n")
cat("Azimuth l3 subtype cell-count audit\n")
cat(sprintf("Threshold for inclusion in pseudobulk: n_cells >= %d per donor\n", THRESH))
cat("================================================================\n\n")

for (tis in c("CSF", "PBMC")) {
  cat("##", tis, "------------------------------------------------\n\n")
  sub <- summary_tab %>% filter(tissue == tis)
  for (parent in c("CD4_T", "CD8_T", "DC", "B", "Mono", "NK", "dnT", "gdT")) {
    rows <- sub %>% filter(grepl(parent, l2, ignore.case = TRUE) | l2 == parent)
    if (nrow(rows) == 0) next
    cat("###", parent, "(parent l2)\n")
    print(rows %>% select(l3, n_donors_pass, n_donors_total,
                          median_cells, min_cells, max_cells))
    cat("\n")
  }
}
sink()
cat("Wrote:", file.path(OUT_DIR, "l3_audit_summary.txt"), "\n\n")

cat("Done. Open the summary file to see which l3 subtypes are pseudobulk-ready.\n")
cat("Recommended next step: rerun the QUBO pipeline with l3 subtypes that have\n")
cat(">= 30 donors passing the 20-cell threshold (gives enough statistical power).\n")
