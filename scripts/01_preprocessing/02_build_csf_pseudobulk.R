# ==============================================================================
# 02_build_csf_pseudobulk.R
# ------------------------------------------------------------------------------
# CSF donor x cell-type pseudobulk over ALL genes, and edgeR statistics for every
# outer fold of the three leave-one-cohort-out settings.
#
#   1. Aggregation (Seurat object loaded once): per donor x cell type, the sum of the
#      RNA counts layer (input to edgeR) and the mean of the log-normalised data layer
#      (used for redundancy and classification)                -> csf_allgenes/
#   2. Outer folds: fixed donor-level assignment in data_release/fold_assignments.csv
#      (Pappalardo, Heming or Ramesh held out; the remaining donors in 5 folds stratified
#      by cohort x diagnosis). For each held-out setting x cell type x fold k, edgeR is
#      fitted on the outer-training donors only (all CV donors not in fold k):
#      filterByExpr + TMM + quasi-likelihood F-test, design ~ Dx + log10_n_cells + age +
#      sex + batch (~ Dx if rank-deficient). The filterByExpr universe of that fold is the
#      candidate gene universe.                                  -> csf_outer/<holdout>/
#
# Input : QUBOFS_SEURAT_RDS_CLEAN (output of 01_clean_cohorts.R)
# Output: QUBOFS_DATA_DIR/csf_allgenes/, QUBOFS_DATA_DIR/csf_outer/<holdout>/
#           folds_assignment.csv, <cell type>/fold_<k>/{tstats_edger_counts.csv, train_donors.csv},
#         QUBOFS_DATA_DIR/csf_outer/build_summary.csv
# Usage : Rscript scripts/01_preprocessing/02_build_csf_pseudobulk.R
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(edgeR) })
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(file.path(if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else getwd(), "common.R"))

CLEAN_RDS <- need_env("QUBOFS_SEURAT_RDS_CLEAN")
DATA      <- need_env("QUBOFS_DATA_DIR")
FOLDS     <- Sys.getenv("QUBOFS_FOLDS", file.path(REPO_DIR, "data_release", "fold_assignments.csv"))
TISSUE    <- "CSF"

# ---- 1. aggregation -----------------------------------------------------------
logmsg("Loading %s", CLEAN_RDS)
so <- readRDS(CLEAN_RDS)
DefaultAssay(so) <- "RNA"
so <- subset(so, cells = colnames(so)[so@meta.data[[COMPARTMENT_COL]] %in% TISSUE])
A <- aggregate_pseudobulk(so, file.path(DATA, "csf_allgenes"), TISSUE)
PB <- A$PB; donor_info <- A$donor_info
rm(so); invisible(gc())

# ---- 2. outer folds -------------------------------------------------------------
FA <- read.csv(FOLDS, stringsAsFactors = FALSE, colClasses = "character")
if (!setequal(donor_info$donor_id, unique(FA$donor)))
  stop("Donors in the object differ from ", FOLDS, ": extra ", paste(setdiff(donor_info$donor_id, FA$donor), collapse = ","),
       "; missing ", paste(setdiff(FA$donor, donor_info$donor_id), collapse = ","))
summary_rows <- list()
for (ho in names(HOLDOUTS)) {
  out <- file.path(DATA, "csf_outer", ho); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  folds <- FA[FA$holdout == ho, c("donor", "label", "cohort", "fold", "set")]
  stopifnot(all(folds$cohort[folds$set == "heldout"] == HOLDOUTS[[ho]]))
  write.csv(folds, file.path(out, "folds_assignment.csv"), row.names = FALSE)
  logmsg("=== held-out %s -> %s", ho, out)
  for (ct in names(PB)) {
    P <- PB[[ct]]
    for (k in sort(unique(folds$fold[folds$set == "cv"]))) {
      tr_d <- folds$donor[folds$set == "cv" & folds$fold != k]
      res <- run_edger(P, tr_d, donor_info)
      if (is.null(res)) { logmsg("[skip] %s fold %s: training donors lack both diagnoses", ct, k); next }
      write_edger(res, file.path(out, ct, TISSUE, paste0("fold_", k)))
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        holdout = ho, cell_type = ct, fold = k, n_train = length(res$train_donors),
        n_genes_total = nrow(P$counts), n_genes_filterByExpr = nrow(res$ts), design = res$design)
    }
    logmsg("%s done", ct)
  }
}
S <- do.call(rbind, summary_rows)
write.csv(S, file.path(DATA, "csf_outer", "build_summary.csv"), row.names = FALSE)
logmsg("DONE: %d folds written", nrow(S))
print(aggregate(n_genes_filterByExpr ~ cell_type, S, function(x) round(c(min = min(x), median = median(x), max = max(x)))))
