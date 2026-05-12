# ==============================================================================
# extract_pseudobulk_v7_l3.R
# ------------------------------------------------------------------------------
# v7 (2026-05-07): Azimuth predicted.celltype.l3 を使った fine-grained pseudobulk.
#
# v5 との違い:
#   * CELLTYPE_COL: predicted.celltype.l2  ->  predicted.celltype.l3
#   * COMPARTMENT_COL: "compartment" -> "Tx" (値 "PB" を "PBMC" に正規化)
#   * RDS_PATH: master (l3 を含む) を使用
#   * CELL_TYPES: l3 audit (l3_donors_passing_20.csv) で
#       n_donors_pass >= 30 / 50 (CSF) または >= 20 / 28 (PBMC) を満たす subtype のみ
#   * 出力先: pseudobulk_v7_l3/
#
# audit summary (CSF, n_donors_pass >= 30):
#   CD4 TCM_1 (50), CD4 TCM_2 (49), CD4 TCM_3 (50),
#   CD4 TEM_1 (42), CD4 TEM_2 (41),
#   CD8 TEM_1 (49), CD8 TEM_2 (47), CD8 TCM_1 (48),
#   cDC2_1 (37), cDC2_2 (43),
#   NK_CD56bright (34), dnT_2 (32),
#   CD14 Mono (29), CD16 Mono (26) ← borderline
#   B intermediate kappa/lambda (11-14) ← borderline (use anyway, was empty before)
#
# 想定ランタイム: 60-120 min (l3 は cell-type 数が l2 の約 2 倍)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# --- パラメータ -----------------------------------------------------------------
RDS_PATH <- "/Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/so.GEX.share.Asada.rds"
OUT_BASE <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/data/pseudobulk_v7_l3"

# rotate held 用
if (!exists("holdout_cohort_arg")) {
  HOLDOUT_COHORTS <- c("PRJNA671484_MS_Tcell")  # = Pappalardo (default)
} else {
  HOLDOUT_COHORTS <- holdout_cohort_arg
}
HOLDOUT_TAG <- if (identical(HOLDOUT_COHORTS, c("PRJNA671484_MS_Tcell")))
                  "" else paste0("_holdout_", paste(HOLDOUT_COHORTS, collapse = "_"))
OUT_ROOT <- if (HOLDOUT_TAG == "") OUT_BASE else paste0(OUT_BASE, HOLDOUT_TAG)

N_FOLDS <- 5
SEED <- 42

# meta.data 列名
DX_COL          <- "Dx"
DONOR_COL       <- "patient_id"
COHORT_COL      <- "prj"
COMPARTMENT_COL <- "Tx"            # v7: master object uses "Tx" (CSF/PB)
CELLTYPE_COL    <- "predicted.celltype.l3"  # v7: l3 instead of l2

# 共変量候補
COV_AGE_COLS   <- c("age", "Age", "age_years")
COV_SEX_COLS   <- c("sex", "Sex", "gender", "Gender")
COV_BATCH_COLS <- c("batch", "Batch", "library", "donor_batch", "tenx_run", "prj")

# --- v7 細胞型グループ (l3 audit に基づく) -------------------------------------
# 各エントリ: name -> vector of l3 labels to aggregate (usually a single l3 label,
#   but we can collapse some near-equivalents if needed).
# Sanitization: spaces and special chars in keys avoided so directory names are safe.
CELL_TYPES <- list(
  # B
  B_intermediate   = c("B intermediate kappa", "B intermediate lambda"),
  B_memory         = c("B memory kappa", "B memory lambda"),
  B_naive          = c("B naive kappa", "B naive lambda"),
  Plasmablast      = c("Plasmablast"),
  # Mono
  CD14_Mono        = c("CD14 Mono"),
  CD16_Mono        = c("CD16 Mono"),
  # CD4 ← v5 では空だった
  CD4_TCM_1        = c("CD4 TCM_1"),
  CD4_TCM_2        = c("CD4 TCM_2"),
  CD4_TCM_3        = c("CD4 TCM_3"),
  CD4_TEM_1        = c("CD4 TEM_1"),
  CD4_TEM_2        = c("CD4 TEM_2"),
  CD4_TEM_3        = c("CD4 TEM_3"),
  CD4_Naive        = c("CD4 Naive"),
  CD4_CTL          = c("CD4 CTL"),
  Treg_Memory      = c("Treg Memory", "Treg Naive"),
  # CD8 ← v5 では空だった
  CD8_TEM_1        = c("CD8 TEM_1"),
  CD8_TEM_2        = c("CD8 TEM_2"),
  CD8_TEM_other    = c("CD8 TEM_3", "CD8 TEM_4", "CD8 TEM_5", "CD8 TEM_6"),
  CD8_TCM_1        = c("CD8 TCM_1"),
  CD8_TCM_other    = c("CD8 TCM_2", "CD8 TCM_3"),
  CD8_Naive        = c("CD8 Naive", "CD8 Naive_2"),
  MAIT             = c("MAIT"),
  # NK
  NK_CD56bright    = c("NK_CD56bright"),
  NK_other         = c("NK_1", "NK_2", "NK_3", "NK_4"),
  # DC ← v5 では空だった
  cDC1             = c("cDC1"),
  cDC2_1           = c("cDC2_1"),
  cDC2_2           = c("cDC2_2"),
  pDC              = c("pDC"),
  # dnT / gdT
  dnT_1            = c("dnT_1"),
  dnT_2            = c("dnT_2"),
  gdT_1            = c("gdT_1"),
  gdT_other        = c("gdT_2", "gdT_3", "gdT_4")
)
TISSUES <- c("CSF", "PBMC")

N_HVG              <- 3000
TOPN_PER_CELLTYPE  <- 100
SKIP_EXISTING      <- TRUE
MIN_CELLS_PER_DONOR <- 20  # consistent with audit threshold

set.seed(SEED)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ==============================================================================
# step1: 読み込み & meta.data 監査
# ==============================================================================
logmsg("Loading Seurat object: ", RDS_PATH)
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"

logmsg("DietSeurat: drop scale.data, prediction.score.* assays, all reductions")
so <- tryCatch(
  DietSeurat(so, layers = c("counts", "data"), assays = "RNA", dimreducs = NULL),
  error = function(e1) tryCatch(
    DietSeurat(so, slot   = c("counts", "data"), assays = "RNA", dimreducs = NULL),
    error = function(e2) {
      warning("DietSeurat failed; proceeding without diet")
      so
    }
  )
)
gc()

# Tx -> compartment 正規化 (PB -> PBMC)
md <- so@meta.data
v <- as.character(md[[COMPARTMENT_COL]])
v[v %in% c("PB","pb","Blood","blood")] <- "PBMC"
v[v %in% c("csf")] <- "CSF"
md[[COMPARTMENT_COL]] <- v
so@meta.data <- md
logmsg("Tissue distribution after normalization:")
print(table(so@meta.data[[COMPARTMENT_COL]], useNA = "ifany"))

needed <- c(DX_COL, DONOR_COL, COHORT_COL, COMPARTMENT_COL, CELLTYPE_COL)
missing_cols <- setdiff(needed, colnames(so@meta.data))
if (length(missing_cols)) {
  stop("meta.data に必要な列がありません: ", paste(missing_cols, collapse = ", "))
}

logmsg("predicted.celltype.l3 (top 30):")
print(head(sort(table(so@meta.data[[CELLTYPE_COL]]), decreasing = TRUE), 30))

donor_tbl <- so@meta.data %>%
  distinct(donor   = .data[[DONOR_COL]],
           label   = .data[[DX_COL]],
           cohort  = .data[[COHORT_COL]]) %>%
  filter(!is.na(donor), !is.na(label))
logmsg("Donors per cohort × Dx:")
print(donor_tbl %>% count(cohort, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))

# --- 共変量自動判定 -----------------------------------------------------------
pick_first_present <- function(candidates, df) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NA_character_
}
age_col   <- pick_first_present(COV_AGE_COLS,   so@meta.data)
sex_col   <- pick_first_present(COV_SEX_COLS,   so@meta.data)
batch_col <- pick_first_present(COV_BATCH_COLS, so@meta.data)
logmsg("Detected covariates: age=", age_col, ", sex=", sex_col, ", batch=", batch_col)

cov_avail <- tibble(
  variable  = c("age", "sex", "batch"),
  used_col  = c(age_col, sex_col, batch_col)
) %>% mutate(present = !is.na(used_col))
write.csv(cov_avail, file.path(OUT_ROOT, "covariate_availability.csv"), row.names = FALSE)

# ==============================================================================
# step2: cell-level filter (HD/MS only) and fold assignment
# ==============================================================================
logmsg("Filtering to HD/MS donors only (Touil HD-only stays in train)")
keep_dx <- so@meta.data[[DX_COL]] %in% c("HD", "MS")
so <- subset(so, cells = colnames(so)[keep_dx])
gc()
logmsg("Cells after Dx filter: ", ncol(so))

donor_meta <- so@meta.data %>%
  distinct(donor_id = .data[[DONOR_COL]],
           diagnosis = .data[[DX_COL]],
           cohort    = .data[[COHORT_COL]])

# fold assignment (donor-level, stratified by cohort × Dx)
non_holdout <- donor_meta %>% filter(!cohort %in% HOLDOUT_COHORTS)
holdout_meta <- donor_meta %>% filter(cohort %in% HOLDOUT_COHORTS) %>%
  mutate(fold = "heldout", set = "heldout")

set.seed(SEED)
non_holdout <- non_holdout %>%
  group_by(cohort, diagnosis) %>%
  mutate(fold = as.character(((sample(seq_len(n())) - 1) %% N_FOLDS) + 1)) %>%
  ungroup() %>%
  mutate(set = "train_or_val")

folds_df <- bind_rows(non_holdout, holdout_meta)
write.csv(folds_df, file.path(OUT_ROOT, "folds_assignment.csv"), row.names = FALSE)
logmsg("Fold assignment written. Distribution:")
print(folds_df %>% count(fold, diagnosis) %>%
        pivot_wider(names_from = diagnosis, values_from = n, values_fill = 0))

# ==============================================================================
# step3: per (cell_type x tissue x fold) pseudobulk
# ==============================================================================
all_l3_in_data <- unique(as.character(so@meta.data[[CELLTYPE_COL]]))

# Validate CELL_TYPES against actual l3 labels
for (ct_name in names(CELL_TYPES)) {
  defined <- CELL_TYPES[[ct_name]]
  found   <- intersect(defined, all_l3_in_data)
  missing <- setdiff(defined, all_l3_in_data)
  if (length(missing) > 0) {
    logmsg("WARN: ", ct_name, " - missing l3 labels: ", paste(missing, collapse = ", "))
  }
  if (length(found) == 0) {
    logmsg("SKIP: ", ct_name, " - none of its l3 labels found in data")
  }
}

# Helper: aggregate pseudobulk for a cell-subset
make_pb <- function(so_sub, donors, mode = c("mean","sum")) {
  mode <- match.arg(mode)
  layer_name <- if (mode == "mean") "data" else "counts"
  m <- GetAssayData(so_sub, assay = "RNA", layer = layer_name)
  if (is.null(m) || ncol(m) == 0) return(NULL)
  donor_ids <- so_sub@meta.data[[DONOR_COL]]
  result <- matrix(0, nrow = nrow(m), ncol = length(donors),
                   dimnames = list(rownames(m), donors))
  for (d in donors) {
    cells_d <- which(donor_ids == d)
    if (length(cells_d) == 0) next
    if (mode == "mean") {
      result[, d] <- as.numeric(rowMeans(m[, cells_d, drop = FALSE]))
    } else {
      result[, d] <- as.numeric(rowSums(m[, cells_d, drop = FALSE]))
    }
  }
  as(result, "dgCMatrix")
}

write_mtx_with_names <- function(m, base) {
  Matrix::writeMM(m, paste0(base, ".mtx"))
  # Match v5 format: header "gene" / "donor" with quoted values
  write.csv(data.frame(gene = rownames(m)),
            paste0(base, "_rows.csv"), row.names = FALSE)
  write.csv(data.frame(donor = colnames(m)),
            paste0(base, "_cols.csv"), row.names = FALSE)
}

# Helper: build a non-singular design matrix.
# Drops factor levels with no observations and removes linearly dependent columns
# detected via QR rank.
build_design <- function(Dx, cov_df) {
  # Drop factor columns with <2 unique non-NA values (constant) or perfectly
  # collinear with Dx
  keep_cols <- character(0)
  for (cn in colnames(cov_df)) {
    v <- cov_df[[cn]]
    if (is.factor(v)) {
      v <- droplevels(v)
      cov_df[[cn]] <- v
      if (nlevels(v) < 2) next
      # Skip factor columns perfectly collinear with Dx (e.g. cohort that has
      # only one Dx level)
      tab <- table(Dx, v)
      if (all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)) {
        next
      }
    } else {
      if (sd(v, na.rm = TRUE) < 1e-10) next
    }
    keep_cols <- c(keep_cols, cn)
  }
  cov_df <- cov_df[, keep_cols, drop = FALSE]

  # Try full design, then progressively drop covariates if rank-deficient
  designs_to_try <- list(
    if (ncol(cov_df) > 0) model.matrix(~ Dx + ., data = cov_df) else NULL,
    model.matrix(~ Dx)
  )
  designs_to_try <- Filter(Negate(is.null), designs_to_try)

  for (design in designs_to_try) {
    qr_check <- qr(design)
    if (qr_check$rank == ncol(design)) {
      return(list(design = design, kept_covs = colnames(cov_df)))
    }
    # Try dropping rank-deficient columns identified by QR
    if (qr_check$rank > 1 && qr_check$rank < ncol(design)) {
      keep_idx <- qr_check$pivot[seq_len(qr_check$rank)]
      design_red <- design[, keep_idx, drop = FALSE]
      # Make sure DxMS column is still present (column 2 in original; check by name)
      dx_col <- grep("^Dx", colnames(design))[1]
      if (!is.na(dx_col) && !any(colnames(design_red) %in% colnames(design)[dx_col])) next
      return(list(design = design_red, kept_covs = "(reduced)"))
    }
  }
  # Last resort: ~Dx only
  return(list(design = model.matrix(~ Dx), kept_covs = "(none)"))
}

# DEG with covariate-adjusted lm() per gene; robust to singular design matrices
fit_tstats <- function(pb_mat, donor_meta_sub, age_col, sex_col, batch_col) {
  donors <- colnames(pb_mat)
  meta <- donor_meta_sub[match(donors, donor_meta_sub$donor_id), ]
  Dx <- factor(meta$diagnosis, levels = c("HD","MS"))

  covs <- list()
  if (!is.na(age_col) && age_col %in% colnames(meta))   covs$age   <- as.numeric(meta[[age_col]])
  if (!is.na(sex_col) && sex_col %in% colnames(meta))   covs$sex   <- factor(meta[[sex_col]])
  if (!is.na(batch_col) && batch_col %in% colnames(meta)) covs$batch <- factor(meta[[batch_col]])
  if (!is.null(meta$n_cells)) covs$logn <- log10(meta$n_cells + 1)
  cov_df <- as.data.frame(covs)

  des <- build_design(Dx, cov_df)
  design <- des$design
  if (length(des$kept_covs) > 0) {
    logmsg("    design covariates kept: ", paste(des$kept_covs, collapse = ","),
           " (", ncol(design), " columns, n=", nrow(design), ")")
  }

  # Pre-compute (X'X)^-1 via QR for numerical stability
  qr_design <- qr(design)
  if (qr_design$rank < ncol(design)) {
    logmsg("    WARN: design still rank-deficient after pruning; falling back to ~Dx")
    design <- model.matrix(~ Dx)
    qr_design <- qr(design)
  }
  XtX_inv <- chol2inv(qr.R(qr_design))

  # Locate Dx coefficient column
  dx_col_idx <- grep("DxMS", colnames(design))
  if (length(dx_col_idx) == 0) dx_col_idx <- 2
  dx_col_idx <- dx_col_idx[1]

  expr <- as.matrix(pb_mat)
  results <- t(apply(expr, 1, function(y) {
    if (sd(y, na.rm = TRUE) < 1e-10) return(c(NA, NA, NA))
    fit <- try(qr.coef(qr_design, y), silent = TRUE)
    if (inherits(fit, "try-error") || any(is.na(fit))) return(c(NA, NA, NA))
    resid <- y - design %*% fit
    df_resid <- nrow(design) - ncol(design)
    if (df_resid <= 0) return(c(NA, NA, NA))
    rss <- sum(resid^2)
    se <- sqrt(rss / df_resid * XtX_inv[dx_col_idx, dx_col_idx])
    coef <- fit[dx_col_idx]
    if (is.na(se) || se < 1e-12) return(c(NA, NA, NA))
    t <- coef / se
    p <- 2 * pt(-abs(t), df = df_resid)
    c(t = t, pval = p, log2FC = coef)
  }))
  colnames(results) <- c("t", "pval", "log2FC")

  out <- data.frame(gene = rownames(expr), results, stringsAsFactors = FALSE)
  out$padj <- p.adjust(out$pval, method = "BH")
  out <- out[order(-abs(out$t), na.last = TRUE), ]
  out$rank <- seq_len(nrow(out))
  out$top_topn <- out$rank <= TOPN_PER_CELLTYPE
  out
}

# Loop: cell_type × tissue × fold
for (ct_name in names(CELL_TYPES)) {
  l3_labels <- CELL_TYPES[[ct_name]]
  matching_l3 <- intersect(l3_labels, all_l3_in_data)
  if (length(matching_l3) == 0) next

  for (tissue in TISSUES) {
    cells_keep <- which(
      as.character(so@meta.data[[CELLTYPE_COL]]) %in% matching_l3 &
      as.character(so@meta.data[[COMPARTMENT_COL]]) == tissue
    )
    if (length(cells_keep) < 10) {
      logmsg("SKIP: ", ct_name, " / ", tissue, " - only ", length(cells_keep), " cells total")
      next
    }

    so_ct <- subset(so, cells = colnames(so)[cells_keep])

    # cells/donor (for QC)
    cells_per_donor <- so_ct@meta.data %>%
      count(.data[[DONOR_COL]], name = "n_cells")
    qc_donors <- cells_per_donor %>%
      filter(n_cells >= MIN_CELLS_PER_DONOR) %>%
      pull(1)

    if (length(qc_donors) < 10) {
      logmsg("SKIP: ", ct_name, " / ", tissue, " - only ", length(qc_donors),
             " donors pass >=", MIN_CELLS_PER_DONOR, " cells")
      next
    }

    so_ct <- subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% qc_donors])

    for (fold in as.character(seq_len(N_FOLDS))) {
      out_dir <- file.path(OUT_ROOT, ct_name, tissue, paste0("fold_", fold))
      if (SKIP_EXISTING && file.exists(file.path(out_dir, "tstats.csv"))) {
        next
      }
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      train_donors <- folds_df %>%
        filter(set == "train_or_val", fold != !!fold) %>% pull(donor_id) %>%
        intersect(qc_donors)
      val_donors <- folds_df %>%
        filter(set == "train_or_val", fold == !!fold) %>% pull(donor_id) %>%
        intersect(qc_donors)
      heldout_donors <- folds_df %>%
        filter(set == "heldout") %>% pull(donor_id) %>%
        intersect(qc_donors)

      if (length(train_donors) < 5 || length(val_donors) < 1) {
        logmsg("SKIP: ", ct_name, "/", tissue, "/fold ", fold,
               " - insufficient donors (train=", length(train_donors),
               " val=", length(val_donors), ")")
        next
      }

      logmsg(ct_name, "/", tissue, "/fold ", fold,
             " - train=", length(train_donors),
             " val=", length(val_donors),
             " heldout=", length(heldout_donors))

      build_meta <- function(donors, set_label) {
        cpd <- cells_per_donor %>% filter(.data[[DONOR_COL]] %in% donors) %>%
          rename(donor_id = !!DONOR_COL, n_cells = n_cells)
        donor_meta %>%
          filter(donor_id %in% donors) %>%
          left_join(cpd, by = "donor_id") %>%
          mutate(set = set_label,
                 age = if (!is.na(age_col)) so_ct@meta.data[[age_col]][match(donor_id, so_ct@meta.data[[DONOR_COL]])] else NA,
                 sex = if (!is.na(sex_col)) so_ct@meta.data[[sex_col]][match(donor_id, so_ct@meta.data[[DONOR_COL]])] else NA,
                 batch = if (!is.na(batch_col)) so_ct@meta.data[[batch_col]][match(donor_id, so_ct@meta.data[[DONOR_COL]])] else NA,
                 compartment = tissue)
      }
      tm  <- build_meta(train_donors, "train")
      vm  <- build_meta(val_donors,   "val")
      hm  <- build_meta(heldout_donors, "heldout")
      write.csv(tm, file.path(out_dir, "train_meta.csv"),  row.names = FALSE)
      write.csv(vm, file.path(out_dir, "val_meta.csv"),    row.names = FALSE)
      write.csv(hm, file.path(out_dir, "heldout_meta.csv"), row.names = FALSE)

      for (mode in c("mean", "sum")) {
        pb_t <- make_pb(subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% train_donors]),  train_donors, mode)
        pb_v <- make_pb(subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% val_donors]),    val_donors,   mode)
        pb_h <- make_pb(subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% heldout_donors]), heldout_donors, mode)
        if (!is.null(pb_t)) write_mtx_with_names(pb_t, file.path(out_dir, paste0("train_pb_",   mode)))
        if (!is.null(pb_v)) write_mtx_with_names(pb_v, file.path(out_dir, paste0("val_pb_",     mode)))
        if (!is.null(pb_h) && length(heldout_donors) > 0) write_mtx_with_names(pb_h, file.path(out_dir, paste0("heldout_pb_", mode)))
      }

      # HVG (train cells のみで FindVariableFeatures)
      so_train <- subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% train_donors])
      so_train <- FindVariableFeatures(so_train, selection.method = "vst",
                                        nfeatures = N_HVG, verbose = FALSE)
      hvg <- VariableFeatures(so_train)
      write.csv(data.frame(gene = hvg), file.path(out_dir, "HVG.csv"), row.names = FALSE)

      # tstats (train mean pseudobulk + lm())
      pb_train_mean <- make_pb(so_train, train_donors, "mean")
      tm_for_fit <- tm %>% select(donor_id, diagnosis, age, sex, batch, n_cells)
      ts <- tryCatch(
        fit_tstats(pb_train_mean[hvg, , drop = FALSE], tm_for_fit, "age", "sex", "batch"),
        error = function(e) {
          logmsg("    fit_tstats FAILED: ", conditionMessage(e), " - retrying with ~Dx only")
          # Bare-bones fallback: ~Dx only
          tryCatch({
            simple_meta <- tm_for_fit
            simple_meta$age <- NA; simple_meta$sex <- NA; simple_meta$batch <- NA
            fit_tstats(pb_train_mean[hvg, , drop = FALSE], simple_meta, NA, NA, NA)
          }, error = function(e2) {
            logmsg("    Even ~Dx fallback failed: ", conditionMessage(e2), " - skipping")
            NULL
          })
        }
      )
      if (is.null(ts)) {
        logmsg("    SKIP tstats: ", ct_name, "/", tissue, "/fold ", fold)
        next
      }
      write.csv(ts, file.path(out_dir, "tstats.csv"), row.names = FALSE)

      topn <- ts %>% filter(top_topn) %>% select(gene, t, pval, padj, log2FC, rank)
      write.csv(topn, file.path(out_dir, "topN_genes.csv"), row.names = FALSE)

      gc()
    }
  }
}

logmsg("Done. v7 l3 pseudobulk written to: ", OUT_ROOT)
logmsg("Next steps:")
logmsg("  1. Run QUBO pipeline pointing to v7 directory")
logmsg("  2. Compare AUC vs v6 (l2)")
