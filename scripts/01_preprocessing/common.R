# ==============================================================================
# common.R -- shared settings and helpers for the R preprocessing scripts
# (sourced by 02_build_csf_pseudobulk.R, 03_build_csf_inner_cv.R, 04_build_pbmc.R).
# ==============================================================================

# ---- Paths (environment variables) --------------------------------------------
need_env <- function(v) {
  x <- Sys.getenv(v, "")
  if (!nzchar(x)) stop("Set the environment variable ", v, " (see scripts/01_preprocessing/README.md).")
  x
}
script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
REPO_DIR <- normalizePath(file.path(script_dir(), "..", ".."), mustWork = FALSE)

# Output layout under QUBOFS_DATA_DIR (read by scripts/02_benchmark):
#   csf_allgenes/                 donor x cell-type pseudobulk, all genes (CSF)
#   csf_outer/<holdout>/          outer folds: folds_assignment.csv, <ct>/fold_<k>/tstats_edger_counts.csv
#   csf_inner/<holdout>/outer_<k>/rep_1/   inner folds and inner-training edgeR statistics
#   pbmc_allgenes/, pbmc_units/<unit>/      PBMC equivalents
HOLDOUTS <- c(Pappalardo = "PRJNA671484_MS_Tcell", Heming = "osmzhlab_MS_ence_cov", Ramesh = "PRJNA549712_MS_PBMC_UCSF")

# ---- Seurat meta.data columns ---------------------------------------------------
DX_COL <- "Dx"; DONOR_COL <- "patient_id"; COHORT_COL <- "prj"; SID_COL <- "sid"
COMPARTMENT_COL <- "compartment"; CELLTYPE_COL <- "predicted.celltype.l2"
AGE_COL <- "age"; SEX_COL <- "sex"; BATCH_COL <- "prj"

# Eight broad immune cell types (Azimuth level-2 labels grouped)
CELL_TYPES <- list(
  B      = c("B naive", "B memory", "B intermediate", "Plasmablast"),
  Mono   = c("CD14 Mono", "CD16 Mono"),
  CD4_T  = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg"),
  CD8_T  = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "MAIT"),
  NK     = c("NK", "NK_CD56bright"),
  DC     = c("cDC1", "cDC2", "pDC"),
  dnT    = c("dnT"),
  gdT    = c("gdT"))

logmsg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ---- Donor x cell-type pseudobulk -----------------------------------------------
# counts : sum of the RNA "counts" layer over the donor's cells of the cell type (all genes)
# mean   : donor mean of the RNA "data" layer (log-normalised)  -> used for redundancy and classification
# n_cells: cells per donor and cell type
aggregate_pseudobulk <- function(so, out_dir, compartment_label) {
  md <- so@meta.data
  donor_info <- unique(data.frame(donor_id = md[[DONOR_COL]], diagnosis = md[[DX_COL]], cohort = md[[COHORT_COL]],
                                  compartment = compartment_label, age = md[[AGE_COL]], sex = as.character(md[[SEX_COL]]),
                                  batch = as.character(md[[BATCH_COL]]), stringsAsFactors = FALSE))
  stopifnot(!anyDuplicated(donor_info$donor_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cnt_all <- LayerData(so, assay = "RNA", layer = "counts"); dat_all <- LayerData(so, assay = "RNA", layer = "data")
  PB <- list()
  for (ct in names(CELL_TYPES)) {
    cells <- colnames(so)[md[[CELLTYPE_COL]] %in% CELL_TYPES[[ct]]]
    if (!length(cells)) { logmsg("no cells: %s", ct); next }
    fac <- factor(md[cells, DONOR_COL])
    ind <- sparseMatrix(i = seq_along(fac), j = as.integer(fac), x = 1, dims = c(length(fac), nlevels(fac)))
    n <- as.numeric(table(fac))
    counts <- cnt_all[, cells, drop = FALSE] %*% ind; colnames(counts) <- levels(fac)
    meanlg <- (dat_all[, cells, drop = FALSE] %*% ind) %*% Diagonal(x = 1 / n); colnames(meanlg) <- levels(fac)
    PB[[ct]] <- list(counts = counts, mean = meanlg, n_cells = setNames(n, levels(fac)))
    writeMM(as(counts, "CsparseMatrix"), file.path(out_dir, paste0(ct, "_counts.mtx")))
    writeMM(as(meanlg, "CsparseMatrix"), file.path(out_dir, paste0(ct, "_mean.mtx")))
    write.csv(data.frame(gene = rownames(counts)), file.path(out_dir, paste0(ct, "_genes.csv")), row.names = FALSE)
    write.csv(data.frame(donor = colnames(counts)), file.path(out_dir, paste0(ct, "_donors.csv")), row.names = FALSE)
    write.csv(data.frame(donor = levels(fac), n_cells = n), file.path(out_dir, paste0(ct, "_n_cells.csv")), row.names = FALSE)
    logmsg("%s: %d cells, %d donors, %d genes", ct, length(cells), nlevels(fac), nrow(counts))
  }
  write.csv(donor_info, file.path(out_dir, "donor_info.csv"), row.names = FALSE)
  saveRDS(list(PB = PB, donor_info = donor_info, date = Sys.time()), file.path(out_dir, "pb_allgenes.rds"))
  list(PB = PB, donor_info = donor_info)
}

# ---- edgeR on training donors only ------------------------------------------------
# design ~ Dx + log10_n_cells + age + sex + batch (terms dropped when constant or missing; ~ Dx if rank-deficient)
is_full_rank <- function(design) qr(design)$rank == ncol(design)
build_design_terms <- function(meta) {
  terms <- "~ Dx"
  if (all(!is.na(meta$log10_n_cells)) && length(unique(meta$log10_n_cells)) > 1) terms <- paste(terms, "+ log10_n_cells")
  if (all(!is.na(meta$age)) && length(unique(meta$age)) > 1) terms <- paste(terms, "+ age")
  if (all(!is.na(meta$sex)) && length(unique(meta$sex)) > 1) terms <- paste(terms, "+ sex")
  if (all(!is.na(meta$batch)) && length(unique(meta$batch)) > 1) terms <- paste(terms, "+ batch")
  terms
}
safe_design <- function(meta) {
  ds <- build_design_terms(meta)
  mm <- tryCatch(model.matrix(as.formula(ds), data = meta), error = function(e) NULL)
  if (is.null(mm) || !is_full_rank(mm)) "~ Dx" else ds
}
# Returns NULL if the training donors lack one diagnosis or are fewer than 4; otherwise
# list(ts = gene statistics on the filterByExpr universe, design = design string).
# Statistic t = sign(logFC) * sqrt(F) of the quasi-likelihood F-test for Dx (MS vs HD).
run_edger <- function(P, tr_d, donor_info) {
  tr_d <- intersect(tr_d, colnames(P$counts))
  mtr <- donor_info[match(tr_d, donor_info$donor_id), ]; mtr$n_cells <- as.integer(P$n_cells[tr_d])
  if (length(unique(mtr$diagnosis)) < 2 || length(tr_d) < 4) return(NULL)
  mtr$Dx <- factor(mtr$diagnosis, levels = c("HD", "MS")); mtr$log10_n_cells <- log10(pmax(mtr$n_cells, 1))
  ds <- safe_design(mtr); design <- model.matrix(as.formula(ds), data = mtr)
  cnt_tr <- as.matrix(P$counts[, tr_d, drop = FALSE]); cnt_tr <- round(cnt_tr); storage.mode(cnt_tr) <- "integer"
  y <- DGEList(counts = cnt_tr)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  qlf <- glmQLFTest(glmQLFit(y, design), coef = "DxMS")
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  ts <- data.frame(gene = rownames(tt), t = sign(tt$logFC) * sqrt(pmax(tt$F, 0)), pval = tt$PValue, padj = tt$FDR,
                   log2FC = tt$logFC, beta = tt$logFC * log(2), se = NA_real_)
  ts$t[is.na(ts$t)] <- 0
  ts$rank <- rank(-abs(ts$t), ties.method = "first")
  list(ts = ts, design = ds, train_donors = tr_d)
}
write_edger <- function(res, fd) {
  dir.create(fd, recursive = TRUE, showWarnings = FALSE)
  write.csv(res$ts, file.path(fd, "tstats_edger_counts.csv"), row.names = FALSE)
  write.csv(data.frame(donor = res$train_donors), file.path(fd, "train_donors.csv"), row.names = FALSE)
}

# Shuffle within each stratum, then assign folds round robin (seeded)
strat_split <- function(d, k, seed, by = c("cohort", "diagnosis")) {
  set.seed(seed)
  key <- do.call(paste, d[by]); d <- d[order(key, sample(nrow(d))), ]; key <- do.call(paste, d[by])
  d$fold <- ave(seq_len(nrow(d)), key, FUN = function(i) ((seq_along(i) - 1) %% k) + 1)
  d[order(d$fold, d$cohort, d$donor), ]
}
