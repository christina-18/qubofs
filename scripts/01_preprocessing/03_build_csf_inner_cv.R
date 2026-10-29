# ==============================================================================
# 03_build_csf_inner_cv.R
# ------------------------------------------------------------------------------
# Inner cross-validation inputs without information leakage. For each held-out
# setting and outer fold k:
#   outer-training donors (CV donors not in fold k) -> 5 inner folds stratified by
#   cohort x diagnosis (shuffle within stratum, round robin; fixed seed);
#   per cell type and inner fold j: edgeR (same model as 02_build_csf_pseudobulk.R) on
#   the inner-training donors only. Inner-validation, outer-validation and held-out
#   donors are never passed to edgeR.
# Input : QUBOFS_DATA_DIR/csf_allgenes/pb_allgenes.rds, csf_outer/<holdout>/folds_assignment.csv
# Output: QUBOFS_DATA_DIR/csf_inner/<holdout>/outer_<k>/rep_1/
#           inner_folds.csv, <cell type>/inner_<j>/{tstats_edger_counts.csv, train_donors.csv}
#         QUBOFS_DATA_DIR/csf_inner/design_log.csv
# Usage : Rscript scripts/01_preprocessing/03_build_csf_inner_cv.R
#         (subset for testing: QUBOFS_INNER_HOLDOUTS=Pappalardo QUBOFS_INNER_OUTER=1 QUBOFS_INNER_CTS=B)
# ==============================================================================
suppressPackageStartupMessages({ library(Matrix); library(edgeR) })
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(file.path(if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else getwd(), "common.R"))

DATA <- need_env("QUBOFS_DATA_DIR")
RAW  <- readRDS(file.path(DATA, "csf_allgenes", "pb_allgenes.rds"))
PB <- RAW$PB; donor_info <- RAW$donor_info
OUTB <- file.path(DATA, "csf_inner"); dir.create(OUTB, showWarnings = FALSE)
split_env <- function(v, default) { x <- Sys.getenv(v, ""); if (nzchar(x)) strsplit(x, ",")[[1]] else default }
HOS     <- split_env("QUBOFS_INNER_HOLDOUTS", names(HOLDOUTS))
OUTERS  <- as.integer(split_env("QUBOFS_INNER_OUTER", as.character(1:5)))
CTS     <- split_env("QUBOFS_INNER_CTS", names(PB))
INNER_K <- 5L
REP     <- 1L
OVERWRITE <- Sys.getenv("QUBOFS_OVERWRITE", "0") == "1"

inner_split <- function(d, k, seed) {           # stratified by cohort x diagnosis, round robin after shuffle
  set.seed(seed)
  d <- d[order(d$cohort, d$diagnosis, sample(nrow(d))), ]
  d$inner_fold <- ave(seq_len(nrow(d)), d$cohort, d$diagnosis, FUN = function(i) ((seq_along(i) - 1) %% k) + 1)
  d[order(d$inner_fold, d$cohort, d$donor), ]
}

design_rows <- list()
for (ho in HOS) {
  folds <- read.csv(file.path(DATA, "csf_outer", ho, "folds_assignment.csv"), stringsAsFactors = FALSE)
  for (k in OUTERS) {
    od <- file.path(OUTB, ho, paste0("outer_", k), paste0("rep_", REP)); dir.create(od, recursive = TRUE, showWarnings = FALSE)
    outer_tr <- folds$donor[folds$set == "cv" & folds$fold != k]
    held <- folds$donor[folds$set == "heldout"]; outer_val <- folds$donor[folds$set == "cv" & folds$fold == k]
    stopifnot(!length(intersect(outer_tr, c(held, outer_val))))
    d <- donor_info[match(outer_tr, donor_info$donor_id), c("donor_id", "cohort", "diagnosis")]; names(d)[1] <- "donor"
    seed <- 20261019L + 1000L * match(ho, names(HOLDOUTS)) + 10L * k + REP
    sp <- inner_split(d, INNER_K, seed); sp$seed <- seed
    write.csv(sp, file.path(od, "inner_folds.csv"), row.names = FALSE)
    logmsg("=== %s outer %d: %d outer-training donors -> %d inner folds (seed %d)", ho, k, nrow(sp), INNER_K, seed)
    for (ct in CTS) {
      for (j in seq_len(INNER_K)) {
        fd <- file.path(od, ct, paste0("inner_", j))
        if (file.exists(file.path(fd, "tstats_edger_counts.csv")) && !OVERWRITE) next
        tr_d <- sp$donor[sp$inner_fold != j]; va_d <- sp$donor[sp$inner_fold == j]
        stopifnot(!length(intersect(tr_d, va_d)), all(tr_d %in% outer_tr))       # leak guards
        res <- run_edger(PB[[ct]], tr_d, donor_info)
        if (is.null(res)) { logmsg("[skip] %s inner %d: training donors lack both diagnoses", ct, j); next }
        write_edger(res, fd)
        design_rows[[length(design_rows) + 1]] <- data.frame(holdout = ho, outer = k, cell_type = ct, inner = j,
          n_train = length(res$train_donors), n_val = length(va_d), n_genes = nrow(res$ts), design = res$design)
      }
    }
  }
}
if (length(design_rows)) {
  D <- do.call(rbind, design_rows); f <- file.path(OUTB, "design_log.csv")
  write.table(D, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
  logmsg("edgeR fits: %d", nrow(D))
}
logmsg("DONE")
