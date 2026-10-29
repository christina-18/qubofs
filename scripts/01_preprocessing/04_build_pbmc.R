# ==============================================================================
# 04_build_pbmc.R
# ------------------------------------------------------------------------------
# PBMC inputs (two cohorts with MS and control donors: Ramesh PRJNA549712 and
# Pappalardo PRJNA671484).
#   1. Donor x cell-type pseudobulk over all genes, PBMC compartment (same aggregation
#      and cell-type groups as the CSF analysis)                  -> pbmc_allgenes/
#   2. Evaluation units                                          -> pbmc_units/<unit>/
#        cv_r<r>_t<t> : nested CV, test fold t of repeat r (5 folds x 3 repeats,
#                       stratified by cohort x diagnosis, seeds 20261020-20261022)
#        loco_<test>  : train on the other cohort, test on <test> (reference analysis)
#      per unit: split.csv (role train/test), inner_folds.csv (training donors only),
#      final/<ct>/  edgeR on all training donors; <ct>/inner_<j>/  edgeR on inner-training donors.
# Input : QUBOFS_SEURAT_RDS_CLEAN (output of 01_clean_cohorts.R)
# Output: QUBOFS_DATA_DIR/pbmc_allgenes/, QUBOFS_DATA_DIR/pbmc_units/
# Usage : Rscript scripts/01_preprocessing/04_build_pbmc.R
#         QUBOFS_PBMC_SKIP_AGG=1 reuses step 1; existing edgeR fits are skipped (QUBOFS_OVERWRITE=1 to redo).
# ==============================================================================
suppressPackageStartupMessages({ library(Matrix); library(edgeR) })
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(file.path(if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else getwd(), "common.R"))

DATA <- need_env("QUBOFS_DATA_DIR")
RAWD <- file.path(DATA, "pbmc_allgenes"); UNITS <- file.path(DATA, "pbmc_units")
dir.create(UNITS, recursive = TRUE, showWarnings = FALSE)
COHORTS <- c(Ramesh = "PRJNA549712_MS_PBMC_UCSF", Pappalardo = "PRJNA671484_MS_Tcell")
OVERWRITE <- Sys.getenv("QUBOFS_OVERWRITE", "0") == "1"

# ---- 1. aggregation ----
if (Sys.getenv("QUBOFS_PBMC_SKIP_AGG", "0") == "1" && file.exists(file.path(RAWD, "pb_allgenes.rds"))) {
  RAW <- readRDS(file.path(RAWD, "pb_allgenes.rds")); PB <- RAW$PB; donor_info <- RAW$donor_info
  logmsg("reusing %s", RAWD)
} else {
  suppressPackageStartupMessages(library(Seurat))
  CLEAN_RDS <- need_env("QUBOFS_SEURAT_RDS_CLEAN")
  logmsg("loading %s", CLEAN_RDS)
  so <- readRDS(CLEAN_RDS); DefaultAssay(so) <- "RNA"
  md0 <- so@meta.data
  keep <- md0[[COMPARTMENT_COL]] %in% "PBMC" & md0[[COHORT_COL]] %in% COHORTS & md0[[DX_COL]] %in% c("MS", "HD")
  so <- subset(so, cells = colnames(so)[keep])
  A <- aggregate_pseudobulk(so, RAWD, "PBMC"); PB <- A$PB; donor_info <- A$donor_info
  rm(so); invisible(gc())
  print(table(donor_info$cohort, donor_info$diagnosis))
}

# ---- 2. units and edgeR ----
fit_edger <- function(ct, tr_d, fd) {
  if (file.exists(file.path(fd, "tstats_edger_counts.csv")) && !OVERWRITE) return("exists")
  res <- run_edger(PB[[ct]], tr_d, donor_info)
  if (is.null(res)) return("skip")
  write_edger(res, fd); paste(res$design, nrow(res$ts), sep = "|")
}
D <- donor_info[, c("donor_id", "cohort", "diagnosis")]; names(D)[1] <- "donor"
units <- list()
for (r in 1:3) {
  sp <- strat_split(D, 5, 20261019L + r)                    # seeds 20261020, 20261021, 20261022
  for (t in 1:5) units[[sprintf("cv_r%d_t%d", r, t)]] <- list(test = sp$donor[sp$fold == t], inner_k = 5L,
                                                             inner_seed = 20261020L + 100L * r + t, by = c("cohort", "diagnosis"))
}
for (te in names(COHORTS)) units[[paste0("loco_", te)]] <- list(test = D$donor[D$cohort == COHORTS[[te]]],
  inner_k = if (te == "Pappalardo") 3L else 5L,              # training = the other cohort: Ramesh (3 controls) -> 3 folds
  inner_seed = 20261020L + 900L + match(te, names(COHORTS)), by = "diagnosis")
log <- list()
for (u in names(units)) {
  U <- units[[u]]; ud <- file.path(UNITS, u); dir.create(ud, recursive = TRUE, showWarnings = FALSE)
  tr <- setdiff(D$donor, U$test); stopifnot(!length(intersect(tr, U$test)))
  s <- D; s$role <- ifelse(s$donor %in% U$test, "test", "train"); write.csv(s, file.path(ud, "split.csv"), row.names = FALSE)
  inn <- strat_split(D[D$donor %in% tr, ], U$inner_k, U$inner_seed, U$by); names(inn)[names(inn) == "fold"] <- "inner_fold"
  inn$seed <- U$inner_seed; write.csv(inn, file.path(ud, "inner_folds.csv"), row.names = FALSE)
  for (ct in names(PB)) {
    log[[length(log) + 1]] <- c(u, ct, "final", fit_edger(ct, tr, file.path(ud, "final", ct)))
    for (j in seq_len(U$inner_k)) {
      itr <- inn$donor[inn$inner_fold != j]; ival <- inn$donor[inn$inner_fold == j]
      stopifnot(!length(intersect(itr, ival)), !length(intersect(itr, U$test)))      # leak guards
      log[[length(log) + 1]] <- c(u, ct, paste0("inner_", j), fit_edger(ct, itr, file.path(ud, ct, paste0("inner_", j))))
    }
  }
  logmsg("%s: %d train / %d test donors, inner %d folds", u, length(tr), length(U$test), U$inner_k)
}
L <- as.data.frame(do.call(rbind, log)); names(L) <- c("unit", "cell_type", "fit", "result")
write.csv(L, file.path(UNITS, "fit_log.csv"), row.names = FALSE)
logmsg("DONE: %d new fits, %d existing, %d skipped", sum(grepl("\\|", L$result)), sum(L$result == "exists"), sum(L$result == "skip"))
