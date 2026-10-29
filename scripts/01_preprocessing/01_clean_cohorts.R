## =============================================================================
## 01_clean_cohorts.R
##
## Run once after 00_annotate_compartment.R and before the pseudobulk scripts. Saves a
## corrected copy of the Seurat object with two changes:
##
##   (1) Exclude the 2 CIS donors of the Ramesh cohort.
##       Ramesh et al. (2020) group 12 RRMS and 2 CIS donors under Dx = MS. The CIS
##       donors are removed so that the MS group contains MS only.
##
##   (2) Remove re-sequenced duplicate cells.
##       In the Heming cohort the same 10x library, sequenced again on another date, was
##       deposited as a separate sample (sid), so the same cells were counted 2-3 times
##       (95-99.8% shared cell barcodes). Each cell is kept once (the copy with the highest
##       UMI count; dedup_resequenced_libraries.R). No donor is removed. The rule is applied
##       to all cohorts; only Heming is affected.
##
## Input  : QUBOFS_SEURAT_RDS        (output of 00_annotate_compartment.R)
## Output : QUBOFS_SEURAT_RDS_CLEAN  (corrected .rds), plus logs in the same folder:
##          <name>_summary_by_cohort_dx.csv, dedup_log.csv, excluded_donors.csv
## Usage  : Rscript scripts/01_preprocessing/01_clean_cohorts.R
## =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
## Script folder: works with Rscript (--file=), source() from RStudio, or QUBOFS_PIPELINE_DIR
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (nzchar(Sys.getenv("QUBOFS_PIPELINE_DIR"))) Sys.getenv("QUBOFS_PIPELINE_DIR") else if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f))) else
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
source(file.path(script_dir, "dedup_resequenced_libraries.R"))

## Paths (environment variables)
IN_RDS  <- Sys.getenv("QUBOFS_SEURAT_RDS", unset = "")
OUT_RDS <- Sys.getenv("QUBOFS_SEURAT_RDS_CLEAN", unset = "")
if (IN_RDS == "" || OUT_RDS == "") stop("Set QUBOFS_SEURAT_RDS and QUBOFS_SEURAT_RDS_CLEAN.")
OUT_DIR <- dirname(OUT_RDS); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## --- meta.data column names  ------------------
DX_COL <- "Dx"; DONOR_COL <- "patient_id"; COHORT_COL <- "prj"
SID_COL <- "sid"; TISSUE_COL <- "compartment"

## --- (1) donors to exclude --------------------------------------------------
## The 2 Ramesh CIS donors (MS5_31 and MS12_22).
EXCLUDE_DONORS <- data.frame(
  donor  = c("Ramesh_MS5_31", "Ramesh_MS12_22"),
  reason = "CIS (not MS) in Ramesh 2020; excluded from MS group"
)

## =============================================================================
so <- readRDS(IN_RDS)
## Save memory: downstream scripts use only RNA counts/data, so drop everything else
so <- tryCatch(DietSeurat(so, layers = c("counts", "data"), assays = "RNA", dimreducs = NULL),
               error = function(e) tryCatch(DietSeurat(so, slot = c("counts", "data"), assays = "RNA", dimreducs = NULL),
                                            error = function(e2) so))
invisible(gc())
summarise_obj <- function(md, stage) {
  md %>% filter(.data[[TISSUE_COL]] == "CSF") %>%
    group_by(cohort = .data[[COHORT_COL]], Dx = .data[[DX_COL]]) %>%
    summarise(donors = n_distinct(.data[[DONOR_COL]]), cells = n(), .groups = "drop") %>%
    mutate(stage = stage)
}
md <- so@meta.data
before <- summarise_obj(md, "before")

## (1) CIS exclusion: flag the cells to drop (no copy of the object yet)
hit <- md[[DONOR_COL]] %in% EXCLUDE_DONORS$donor
found <- unique(md[[DONOR_COL]][hit])
if (length(found) != nrow(EXCLUDE_DONORS))
  stop("Exclusion list mismatch: found only ", paste(found, collapse = ", "),
       ". Check donor IDs in this object before proceeding.")
write.csv(EXCLUDE_DONORS, file.path(OUT_DIR, "excluded_donors.csv"), row.names = FALSE)
message(sprintf("[1] excluding %d donors (%d cells)", length(found), sum(hit)))

## (2) Re-sequenced duplicates: decide which cells to keep (after CIS exclusion).
##     Set QUBOFS_DO_DEDUP=0 to skip this step (sensitivity run: CIS exclusion only).
DO_DEDUP <- Sys.getenv("QUBOFS_DO_DEDUP", "1") == "1"
if (DO_DEDUP) {
  keep_dedup <- dedup_keep_cells(md[!hit, ], cohort_col = COHORT_COL, donor_col = DONOR_COL,
                                 sid_col = SID_COL, tissue_col = TISSUE_COL,
                                 log_file = file.path(OUT_DIR, "dedup_log.csv"))
} else {
  message("[2] QUBOFS_DO_DEDUP=0: re-sequenced duplicates are NOT removed")
  keep_dedup <- setNames(rep(TRUE, sum(!hit)), rownames(md)[!hit])
}
keep_cells <- names(keep_dedup)[keep_dedup]
rm(md, keep_dedup); invisible(gc())

## Apply (1) + (2) in a single subset
so <- subset(so, cells = keep_cells); invisible(gc())
message(sprintf("[2] cells kept: %d", ncol(so)))

after <- summarise_obj(so@meta.data, "after")
summ <- full_join(before, after, by = c("cohort", "Dx"), suffix = c("_before", "_after"))
write.csv(summ, file.path(OUT_DIR, sub("\\.rds$", "_summary_by_cohort_dx.csv", basename(OUT_RDS))), row.names = FALSE)
print(summ)

saveRDS(so, OUT_RDS)
message("Wrote cleaned object: ", OUT_RDS)
