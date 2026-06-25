# ==============================================================================
# extend_DEG_methods.R
# ------------------------------------------------------------------------------
# Adds edgeR (primary) plus limma-voom / DESeq2 (optional sensitivity) DEG
# statistics to every fold under each hold-out's pseudobulk directory. Operates
# on donor-wise sums of SoupX-corrected count-scale values in train_pb_counts.mtx;
# values are rounded to integers inside each count-based DE method. The upstream
# Seurat .rds is NOT re-read.
#
# Outputs are written with a "_counts" suffix so they do not overwrite any
# legacy tstats_<method>.csv that were computed from the (incorrect) log1p(CPM)
# *_pb_sum.mtx files:
#   tstats_<method>_counts.csv      full ranking
#   topN_genes_<method>_counts.csv  top-N candidates
# Downstream QUBO selects these via deg_source = "edger_counts" (canonical; deseq2/limmavoom available for sensitivity).
#
# Requires that extract_pseudobulk.R has been re-run so that *_pb_counts.mtx
# exist. Folds lacking train_pb_counts.mtx are skipped with a warning.
#
# Typical runtime: 3 holdouts x 60 folds x 3 methods = ~10-20 min (DESeq2 dominates).
#
# Usage:
#   Set QUBOFS_DATA_BASE to the parent directory of the pseudobulk_* subdirs
#   (defaults to ./data). Then from R within 02_deg/:
#     source("extend_DEG_methods.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
  library(edgeR)
  library(limma)
})

# Data location: prefer QUBOFS_DATA_BASE; otherwise <QUBOFS_PROJECT_ROOT>/data so
# the script works regardless of the current working directory (avoids resolving a
# bare "data" relative to wherever the script is sourced from).
DATA_BASE <- Sys.getenv("QUBOFS_DATA_BASE",
                        unset = file.path(Sys.getenv("QUBOFS_PROJECT_ROOT", "."), "data"))
# Re-compute DE even if outputs already exist (set for a clean canonical re-run).
OVERWRITE <- tolower(Sys.getenv("QUBOFS_OVERWRITE_DEG", "false")) %in% c("true", "1", "yes")
PSEUDOBULK_SUBDIR <- Sys.getenv("QUBOFS_PSEUDOBULK_SUBDIR", unset = "pseudobulk_v5_compartment")
HOLDOUT_DIRS <- c(
  PSEUDOBULK_SUBDIR,                                            # Pappalardo (default)
  paste0(PSEUDOBULK_SUBDIR, "_holdout_osmzhlab_MS_ence_cov"),   # Heming
  paste0(PSEUDOBULK_SUBDIR, "_holdout_PRJNA549712_MS_PBMC_UCSF") # Ramesh
)
CELL_TYPES <- c("B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT")
TISSUES    <- strsplit(Sys.getenv("QUBOFS_TISSUES","CSF"), ",")[[1]]  # set QUBOFS_TISSUES=PBMC for the PBMC benchmark
FOLDS      <- 1:5
TOPN       <- 100

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
load_fold <- function(dir) {
  pb_path <- file.path(dir, "train_pb_counts.mtx")
  if (!file.exists(pb_path)) return(NULL)   # count-scale pseudobulk required; re-run extract_pseudobulk.R
  pb_counts <- as.matrix(readMM(pb_path))
  rows <- read.csv(file.path(dir, "train_pb_counts_rows.csv"), stringsAsFactors = FALSE)$gene
  cols <- read.csv(file.path(dir, "train_pb_counts_cols.csv"), stringsAsFactors = FALSE)$donor
  meta <- read.csv(file.path(dir, "train_meta.csv"),           stringsAsFactors = FALSE)

  # safety: ensure donor columns are unique
  if (anyDuplicated(cols)) {
    keep_idx  <- which(!duplicated(cols))
    pb_counts <- pb_counts[, keep_idx, drop = FALSE]
    cols      <- cols[keep_idx]
  }
  rownames(pb_counts) <- rows
  colnames(pb_counts) <- cols
  meta <- meta[!duplicated(meta$donor_id), ]
  meta <- meta[match(cols, meta$donor_id), ]
  # log10-transformed per-donor cell count covariate (manuscript design term)
  if ("n_cells" %in% names(meta)) {
    meta$log10_n_cells <- log10(pmax(as.numeric(meta$n_cells), 1))
  }
  list(pb_counts = pb_counts, meta = meta)
}

# Full-rank check for a model matrix (small-cohort folds can be rank-deficient).
is_full_rank <- function(design) {
  qr(design)$rank == ncol(design)
}

# Manuscript design: ~ Dx + log10(n_cells) + age + sex + batch, including each
# covariate only when present, non-missing and estimable (>1 level / non-constant).
build_design_terms <- function(meta) {
  terms <- "~ Dx"
  if ("log10_n_cells" %in% names(meta) && all(!is.na(meta$log10_n_cells)) &&
      length(unique(meta$log10_n_cells)) > 1)
    terms <- paste(terms, "+ log10_n_cells")
  if ("age" %in% names(meta) && all(!is.na(meta$age)) && length(unique(meta$age)) > 1)
    terms <- paste(terms, "+ age")
  if ("sex" %in% names(meta) && all(!is.na(meta$sex)) && length(unique(meta$sex)) > 1)
    terms <- paste(terms, "+ sex")
  if ("batch" %in% names(meta) && all(!is.na(meta$batch)) && length(unique(meta$batch)) > 1)
    terms <- paste(terms, "+ batch")
  terms
}

# Build the design formula, falling back to the minimal "~ Dx" if the full
# covariate-adjusted design is rank-deficient on a given (small) fold.
safe_design <- function(meta) {
  ds <- build_design_terms(meta)
  mm <- tryCatch(model.matrix(as.formula(ds), data = meta), error = function(e) NULL)
  if (is.null(mm) || !is_full_rank(mm)) return("~ Dx")
  ds
}

# ------------------------------------------------------------------------------
# DESeq2: NB GLM + Wald test
# ------------------------------------------------------------------------------
# Note: meta$diagnosis uses the internal factor levels "HD" (= control donor)
# and "MS" inherited from the upstream Seurat object. The manuscript reports
# them as "control" and "MS".
run_deseq2 <- function(pb_sum, meta) {
  cnt <- round(pb_sum); storage.mode(cnt) <- "integer"
  meta$Dx <- factor(meta$diagnosis, levels = c("HD", "MS"))
  if (length(unique(meta$Dx)) < 2) return(NULL)

  design_str <- safe_design(meta)
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = meta,
                                design = as.formula(design_str))
  # low-count filter
  keep <- rowSums(counts(dds) >= 1) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "Dx_MS_vs_HD")

  data.frame(
    gene   = rownames(dds),
    t      = res$stat,           # Wald statistic
    pval   = res$pvalue,
    padj   = res$padj,
    log2FC = res$log2FoldChange,
    beta   = res$log2FoldChange * log(2),
    se     = res$lfcSE,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# edgeR: NB GLM + quasi-likelihood F test
# ------------------------------------------------------------------------------
run_edger <- function(pb_sum, meta) {
  cnt <- round(pb_sum); storage.mode(cnt) <- "integer"
  meta$Dx <- factor(meta$diagnosis, levels = c("HD", "MS"))
  if (length(unique(meta$Dx)) < 2) return(NULL)

  design_str <- safe_design(meta)
  design <- model.matrix(as.formula(design_str), data = meta)

  y <- DGEList(counts = cnt)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "DxMS")
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table

  data.frame(
    gene   = rownames(tt),
    t      = sign(tt$logFC) * sqrt(tt$F),  # signed sqrt(F) as a t-like statistic
    pval   = tt$PValue,
    padj   = tt$FDR,
    log2FC = tt$logFC,
    beta   = tt$logFC * log(2),
    se     = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# limma-voom
# ------------------------------------------------------------------------------
run_limmavoom <- function(pb_sum, meta) {
  cnt <- round(pb_sum); storage.mode(cnt) <- "integer"
  meta$Dx <- factor(meta$diagnosis, levels = c("HD", "MS"))
  if (length(unique(meta$Dx)) < 2) return(NULL)

  design_str <- safe_design(meta)
  design <- model.matrix(as.formula(design_str), data = meta)

  y <- DGEList(counts = cnt)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "DxMS", number = Inf, sort.by = "none")

  data.frame(
    gene   = rownames(tt),
    t      = tt$t,
    pval   = tt$P.Value,
    padj   = tt$adj.P.Val,
    log2FC = tt$logFC,
    beta   = tt$logFC * log(2),
    se     = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------
n_done <- 0; n_skip <- 0; n_err <- 0
design_records <- list()
for (hd in HOLDOUT_DIRS) {
  base <- file.path(DATA_BASE, hd)
  if (!dir.exists(base)) {
    logmsg(sprintf("[skip] %s not found", base)); next
  }
  logmsg(sprintf("=== %s ===", hd))
  for (ct in CELL_TYPES) {
    for (tissue in TISSUES) {
      for (fold in FOLDS) {
        fdir <- file.path(base, ct, tissue, sprintf("fold_%d", fold))
        if (!dir.exists(fdir)) next
        d <- load_fold(fdir)
        if (is.null(d)) {
          logmsg(sprintf("[skip] %s/%s/fold_%d: train_pb_counts.mtx missing (re-run extract_pseudobulk.R)",
                         ct, tissue, fold))
          n_skip <- n_skip + 1; next
        }
        if (length(unique(d$meta$diagnosis)) < 2 || ncol(d$pb_counts) < 4) next

        for (m in c("deseq2", "edger", "limmavoom")) {
          out_path <- file.path(fdir, sprintf("tstats_%s_counts.csv", m))
          if (file.exists(out_path) && !OVERWRITE) { n_skip <- n_skip + 1; next }
          res <- tryCatch(
            switch(m,
              deseq2    = run_deseq2(d$pb_counts, d$meta),
              edger     = run_edger(d$pb_counts, d$meta),
              limmavoom = run_limmavoom(d$pb_counts, d$meta)
            ),
            error = function(e) {
              logmsg(sprintf("[err] %s/%s/fold_%d %s: %s",
                             ct, tissue, fold, m, conditionMessage(e)))
              n_err <<- n_err + 1
              NULL
            }
          )
          if (is.null(res) || nrow(res) == 0) next
          # drop genes with non-finite t-statistic or missing p-value before ranking
          res <- res[is.finite(res$t) & !is.na(res$pval), ]
          if (nrow(res) == 0) next
          # record the design formula actually used (reproducibility)
          design_records[[length(design_records) + 1]] <- data.frame(
            holdout = hd, cell_type = ct, tissue = tissue, fold = fold, method = m,
            design = safe_design(d$meta), n_donors = ncol(d$pb_counts),
            n_genes = nrow(res), stringsAsFactors = FALSE)
          # rank by descending |t|
          res <- res[order(-abs(res$t)), ]
          res$rank     <- seq_len(nrow(res))
          res$top_topn <- res$rank <= TOPN
          write.csv(res, out_path, row.names = FALSE)
          # topN
          top <- head(res, TOPN)[, c("gene", "t", "pval", "padj", "log2FC", "rank")]
          write.csv(top, file.path(fdir, sprintf("topN_genes_%s_counts.csv", m)),
                    row.names = FALSE)
          n_done <- n_done + 1
          logmsg(sprintf("[ok]  %s %s/%s/fold_%d %s n=%d",
                         hd, ct, tissue, fold, m, nrow(res)))
        }
      }
    }
  }
}
if (length(design_records) > 0) {
  design_df <- do.call(rbind, design_records)
  ds_path <- file.path(DATA_BASE, "DEG_design_summary.csv")
  write.csv(design_df, ds_path, row.names = FALSE)
  logmsg(sprintf("Wrote design summary: %s (%d rows)", ds_path, nrow(design_df)))
}
logmsg(sprintf("DONE.  ok=%d  skip=%d  err=%d", n_done, n_skip, n_err))
