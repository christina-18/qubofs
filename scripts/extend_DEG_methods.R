# ==============================================================================
# extend_DEG_methods.R
# ------------------------------------------------------------------------------
# 既存の3つの holdout 抽出フォルダを走査し、各 fold に
# DESeq2 / edgeR / limma-voom の DEG 統計量 (tstats_<method>.csv) を追加する。
# 入力は train_pb_sum.mtx (donor-level 生 count 合計) を使うため、.rds の再読込は不要。
#
# 想定ランタイム: 全 3 holdout × 60 fold × 3 方法 = ~10-20 分 (DESeq2 が最重)
#
# Usage:
#   source("/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/scripts/extend_DEG_methods.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
  library(edgeR)
  library(limma)
})

DATA_BASE <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/data"
HOLDOUT_DIRS <- c(
  "pseudobulk_v5_compartment",                              # Pappalardo (default)
  "pseudobulk_v5_compartment_holdout_osmzhlab_MS_ence_cov", # Heming
  "pseudobulk_v5_compartment_holdout_PRJNA549712_MS_PBMC_UCSF" # Ramesh
)
CELL_TYPES <- c("B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT")
TISSUES    <- c("CSF", "PBMC", "ALL")
FOLDS      <- 1:5
TOPN       <- 100

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
load_fold <- function(dir) {
  pb_path <- file.path(dir, "train_pb_sum.mtx")
  if (!file.exists(pb_path)) return(NULL)
  pb_sum <- as.matrix(readMM(pb_path))
  rows <- read.csv(file.path(dir, "train_pb_sum_rows.csv"), stringsAsFactors = FALSE)$gene
  cols <- read.csv(file.path(dir, "train_pb_sum_cols.csv"), stringsAsFactors = FALSE)$donor
  meta <- read.csv(file.path(dir, "train_meta.csv"),         stringsAsFactors = FALSE)

  # ALL tissue では cols (donor) が CSF/PBMC で 2回出ることがある → unique にそろえる
  if (anyDuplicated(cols)) {
    keep_idx <- which(!duplicated(cols))
    pb_sum   <- pb_sum[, keep_idx, drop = FALSE]
    cols     <- cols[keep_idx]
  }
  rownames(pb_sum) <- rows
  colnames(pb_sum) <- cols
  meta <- meta[!duplicated(meta$donor_id), ]
  meta <- meta[match(cols, meta$donor_id), ]
  list(pb_sum = pb_sum, meta = meta)
}

build_design_terms <- function(meta) {
  terms <- "~ Dx"
  if ("age" %in% names(meta) && all(!is.na(meta$age)))
    terms <- paste(terms, "+ age")
  if ("sex" %in% names(meta) && all(!is.na(meta$sex)) && length(unique(meta$sex)) > 1)
    terms <- paste(terms, "+ sex")
  if ("batch" %in% names(meta) && all(!is.na(meta$batch)) && length(unique(meta$batch)) > 1)
    terms <- paste(terms, "+ batch")
  terms
}

# ------------------------------------------------------------------------------
# DESeq2: NB GLM + Wald test
# ------------------------------------------------------------------------------
run_deseq2 <- function(pb_sum, meta) {
  cnt <- round(pb_sum); storage.mode(cnt) <- "integer"
  meta$Dx <- factor(meta$diagnosis, levels = c("HD", "MS"))
  if (length(unique(meta$Dx)) < 2) return(NULL)

  design_str <- build_design_terms(meta)
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = meta,
                                design = as.formula(design_str))
  # 低 count フィルタ
  keep <- rowSums(counts(dds) >= 1) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "Dx_MS_vs_HD")

  data.frame(
    gene   = rownames(dds),
    t      = res$stat,           # Wald 統計量
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

  design_str <- build_design_terms(meta)
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
    t      = sign(tt$logFC) * sqrt(tt$F),  # signed sqrt(F) を t-様統計量として
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

  design_str <- build_design_terms(meta)
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
        if (is.null(d)) next
        if (length(unique(d$meta$diagnosis)) < 2 || ncol(d$pb_sum) < 4) next

        for (m in c("deseq2", "edger", "limmavoom")) {
          out_path <- file.path(fdir, sprintf("tstats_%s.csv", m))
          if (file.exists(out_path)) { n_skip <- n_skip + 1; next }
          res <- tryCatch(
            switch(m,
              deseq2    = run_deseq2(d$pb_sum, d$meta),
              edger     = run_edger(d$pb_sum, d$meta),
              limmavoom = run_limmavoom(d$pb_sum, d$meta)
            ),
            error = function(e) {
              logmsg(sprintf("[err] %s/%s/fold_%d %s: %s",
                             ct, tissue, fold, m, conditionMessage(e)))
              n_err <<- n_err + 1
              NULL
            }
          )
          if (is.null(res) || nrow(res) == 0) next
          # |t| 降順で rank
          res <- res[order(-abs(res$t)), ]
          res$rank     <- seq_len(nrow(res))
          res$top_topn <- res$rank <= TOPN
          write.csv(res, out_path, row.names = FALSE)
          # topN
          top <- head(res, TOPN)[, c("gene", "t", "pval", "padj", "log2FC", "rank")]
          write.csv(top, file.path(fdir, sprintf("topN_genes_%s.csv", m)),
                    row.names = FALSE)
          n_done <- n_done + 1
          logmsg(sprintf("[ok]  %s %s/%s/fold_%d %s n=%d",
                         hd, ct, tissue, fold, m, nrow(res)))
        }
      }
    }
  }
}
logmsg(sprintf("DONE.  ok=%d  skip=%d  err=%d", n_done, n_skip, n_err))
