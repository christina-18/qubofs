# ==============================================================================
# extract_pseudobulk_v8_mixed.R
# ------------------------------------------------------------------------------
# v8 mixed (2026-05-07): l2 を base にしつつ、CD4_T / CD8_T / DC のみ l3-based
# な aggregated subtype に分ける「中間版」.
#
# v6 (l2 only, 8 types) と v7 (l3 only, 32 types) の中間.
# 目的: l2 で空だった CD4/CD8/DC を 2-3 subtype に分けて recover しつつ,
#       LOCO で 1 cohort 抜けても全 subtype が donor 数を維持できるようにする.
#
# Cell types (n = 14):
#   既存 l2 ベース (8):
#     B, Mono, NK, dnT, gdT, B_naive, B_memory, Plasmablast (B 細分化のみ optional)
#   l2 で空だったところを selective aggregation (6):
#     CD4_TCM    = CD4 TCM_1 + TCM_2 + TCM_3
#     CD4_TEM    = CD4 TEM_1 + TEM_2 + TEM_3
#     CD4_other  = CD4 Naive + CD4 CTL + Treg Memory + Treg Naive + CD4 Proliferating
#     CD8_TEM    = CD8 TEM_1 + TEM_2 + TEM_3 + TEM_4 + TEM_5 + TEM_6
#     CD8_TCM    = CD8 TCM_1 + TCM_2 + TCM_3
#     CD8_other  = CD8 Naive + CD8 Naive_2 + CD8 Proliferating + MAIT
#     cDC1, cDC2 (= cDC2_1 + cDC2_2), pDC
#
# 想定ランタイム: ~30-60 min (v7 より cell type 数が少ない)
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
OUT_BASE <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/data/pseudobulk_v8_mixed"

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
COMPARTMENT_COL <- "Tx"
# v8 では l3 を使うが aggregate 単位は CELL_TYPES list で決める
CELLTYPE_COL    <- "predicted.celltype.l3"

COV_AGE_COLS   <- c("age", "Age", "age_years")
COV_SEX_COLS   <- c("sex", "Sex", "gender", "Gender")
COV_BATCH_COLS <- c("batch", "Batch", "library", "donor_batch", "tenx_run", "prj")

# --- v8 mixed 細胞型グループ (selective aggregation) ----------------------------
# 14 cell types: l2-equivalent broad + selective l3-aggregated splits
CELL_TYPES <- list(
  # Standard l2 broad groupings (already worked in v6)
  B          = c("B naive kappa", "B naive lambda",
                 "B memory kappa", "B memory lambda",
                 "B intermediate kappa", "B intermediate lambda",
                 "Plasmablast"),
  Mono       = c("CD14 Mono", "CD16 Mono"),
  NK         = c("NK_1", "NK_2", "NK_3", "NK_4", "NK_CD56bright"),
  dnT        = c("dnT_1", "dnT_2"),
  gdT        = c("gdT_1", "gdT_2", "gdT_3", "gdT_4"),
  # CD4_T split into 3 subgroups (was empty in v6 l2)
  CD4_TCM    = c("CD4 TCM_1", "CD4 TCM_2", "CD4 TCM_3"),
  CD4_TEM    = c("CD4 TEM_1", "CD4 TEM_2", "CD4 TEM_3", "CD4 TEM_4"),
  CD4_other  = c("CD4 Naive", "CD4 CTL", "Treg Memory", "Treg Naive",
                 "CD4 Proliferating"),
  # CD8_T split into 3 subgroups (was empty in v6 l2)
  CD8_TEM    = c("CD8 TEM_1", "CD8 TEM_2", "CD8 TEM_3",
                 "CD8 TEM_4", "CD8 TEM_5", "CD8 TEM_6"),
  CD8_TCM    = c("CD8 TCM_1", "CD8 TCM_2", "CD8 TCM_3"),
  CD8_other  = c("CD8 Naive", "CD8 Naive_2", "CD8 Proliferating", "MAIT"),
  # DC split into 3 subgroups (was empty in v6 l2)
  cDC1       = c("cDC1"),
  cDC2       = c("cDC2_1", "cDC2_2"),
  pDC        = c("pDC")
)
TISSUES <- c("CSF", "PBMC")

N_HVG               <- 3000
TOPN_PER_CELLTYPE   <- 100
SKIP_EXISTING       <- TRUE
MIN_CELLS_PER_DONOR <- 20

set.seed(SEED)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ==============================================================================
# step1: 読み込み & meta.data 監査
# ==============================================================================
logmsg("Loading Seurat object: ", RDS_PATH)
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"

logmsg("DietSeurat")
so <- tryCatch(
  DietSeurat(so, layers = c("counts", "data"), assays = "RNA", dimreducs = NULL),
  error = function(e1) tryCatch(
    DietSeurat(so, slot   = c("counts", "data"), assays = "RNA", dimreducs = NULL),
    error = function(e2) so
  )
)
gc()

# Tx -> compartment 正規化
md <- so@meta.data
v <- as.character(md[[COMPARTMENT_COL]])
v[v %in% c("PB","pb","Blood","blood")] <- "PBMC"
v[v %in% c("csf")] <- "CSF"
md[[COMPARTMENT_COL]] <- v
so@meta.data <- md
logmsg("Tissue distribution:")
print(table(so@meta.data[[COMPARTMENT_COL]], useNA = "ifany"))

needed <- c(DX_COL, DONOR_COL, COHORT_COL, COMPARTMENT_COL, CELLTYPE_COL)
missing_cols <- setdiff(needed, colnames(so@meta.data))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

donor_tbl <- so@meta.data %>%
  distinct(donor   = .data[[DONOR_COL]],
           label   = .data[[DX_COL]],
           cohort  = .data[[COHORT_COL]]) %>%
  filter(!is.na(donor), !is.na(label))
logmsg("Donors per cohort × Dx:")
print(donor_tbl %>% count(cohort, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))

pick_first_present <- function(candidates, df) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NA_character_
}
age_col   <- pick_first_present(COV_AGE_COLS,   so@meta.data)
sex_col   <- pick_first_present(COV_SEX_COLS,   so@meta.data)
batch_col <- pick_first_present(COV_BATCH_COLS, so@meta.data)

cov_avail <- tibble(
  variable  = c("age", "sex", "batch"),
  used_col  = c(age_col, sex_col, batch_col)
) %>% mutate(present = !is.na(used_col))
write.csv(cov_avail, file.path(OUT_ROOT, "covariate_availability.csv"), row.names = FALSE)

# ==============================================================================
# step2: cell-level filter (HD/MS only) and fold assignment
# ==============================================================================
logmsg("Filtering to HD/MS donors only")
keep_dx <- so@meta.data[[DX_COL]] %in% c("HD", "MS")
so <- subset(so, cells = colnames(so)[keep_dx])
gc()

donor_meta <- so@meta.data %>%
  distinct(donor_id = .data[[DONOR_COL]],
           diagnosis = .data[[DX_COL]],
           cohort    = .data[[COHORT_COL]])

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
logmsg("Fold assignment:")
print(folds_df %>% count(fold, diagnosis) %>%
        pivot_wider(names_from = diagnosis, values_from = n, values_fill = 0))

# ==============================================================================
# step3: per (cell_type x tissue x fold) pseudobulk
# ==============================================================================
all_l3 <- unique(as.character(so@meta.data[[CELLTYPE_COL]]))

# Helper: design builder (rank-deficiency safe, from v7)
build_design <- function(Dx, cov_df) {
  keep_cols <- character(0)
  for (cn in colnames(cov_df)) {
    v <- cov_df[[cn]]
    if (is.factor(v)) {
      v <- droplevels(v)
      cov_df[[cn]] <- v
      if (nlevels(v) < 2) next
      tab <- table(Dx, v)
      if (all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)) next
    } else {
      if (sd(v, na.rm = TRUE) < 1e-10) next
    }
    keep_cols <- c(keep_cols, cn)
  }
  cov_df <- cov_df[, keep_cols, drop = FALSE]
  designs_to_try <- list(
    if (ncol(cov_df) > 0) model.matrix(~ Dx + ., data = cov_df) else NULL,
    model.matrix(~ Dx)
  )
  designs_to_try <- Filter(Negate(is.null), designs_to_try)
  for (design in designs_to_try) {
    qr_check <- qr(design)
    if (qr_check$rank == ncol(design)) return(list(design = design, kept_covs = colnames(cov_df)))
    if (qr_check$rank > 1 && qr_check$rank < ncol(design)) {
      keep_idx <- qr_check$pivot[seq_len(qr_check$rank)]
      design_red <- design[, keep_idx, drop = FALSE]
      dx_col <- grep("^Dx", colnames(design))[1]
      if (!is.na(dx_col) && !any(colnames(design_red) %in% colnames(design)[dx_col])) next
      return(list(design = design_red, kept_covs = "(reduced)"))
    }
  }
  return(list(design = model.matrix(~ Dx), kept_covs = "(none)"))
}

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
  qr_design <- qr(design)
  if (qr_design$rank < ncol(design)) {
    design <- model.matrix(~ Dx); qr_design <- qr(design)
  }
  XtX_inv <- chol2inv(qr.R(qr_design))
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
    if (mode == "mean") result[, d] <- as.numeric(rowMeans(m[, cells_d, drop = FALSE]))
    else                result[, d] <- as.numeric(rowSums(m[, cells_d, drop = FALSE]))
  }
  as(result, "dgCMatrix")
}

write_mtx_with_names <- function(m, base) {
  Matrix::writeMM(m, paste0(base, ".mtx"))
  write.csv(data.frame(gene = rownames(m)), paste0(base, "_rows.csv"), row.names = FALSE)
  write.csv(data.frame(donor = colnames(m)), paste0(base, "_cols.csv"), row.names = FALSE)
}

# Loop: cell_type × tissue × fold
for (ct_name in names(CELL_TYPES)) {
  l3_labels <- CELL_TYPES[[ct_name]]
  matching_l3 <- intersect(l3_labels, all_l3)
  if (length(matching_l3) == 0) next

  for (tissue in TISSUES) {
    cells_keep <- which(
      as.character(so@meta.data[[CELLTYPE_COL]]) %in% matching_l3 &
      as.character(so@meta.data[[COMPARTMENT_COL]]) == tissue
    )
    if (length(cells_keep) < 10) {
      logmsg("SKIP: ", ct_name, " / ", tissue, " - ", length(cells_keep), " cells")
      next
    }
    so_ct <- subset(so, cells = colnames(so)[cells_keep])
    cells_per_donor <- so_ct@meta.data %>% count(.data[[DONOR_COL]], name = "n_cells")
    qc_donors <- cells_per_donor %>% filter(n_cells >= MIN_CELLS_PER_DONOR) %>% pull(1)
    if (length(qc_donors) < 10) {
      logmsg("SKIP: ", ct_name, " / ", tissue, " - ", length(qc_donors), " donors pass QC")
      next
    }
    so_ct <- subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% qc_donors])

    for (fold in as.character(seq_len(N_FOLDS))) {
      out_dir <- file.path(OUT_ROOT, ct_name, tissue, paste0("fold_", fold))
      if (SKIP_EXISTING && file.exists(file.path(out_dir, "tstats.csv"))) next
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      train_donors <- folds_df %>% filter(set == "train_or_val", fold != !!fold) %>%
        pull(donor_id) %>% intersect(qc_donors)
      val_donors <- folds_df %>% filter(set == "train_or_val", fold == !!fold) %>%
        pull(donor_id) %>% intersect(qc_donors)
      heldout_donors <- folds_df %>% filter(set == "heldout") %>%
        pull(donor_id) %>% intersect(qc_donors)

      if (length(train_donors) < 5 || length(val_donors) < 1) {
        logmsg("SKIP: ", ct_name, "/", tissue, "/fold ", fold,
               " - train=", length(train_donors), " val=", length(val_donors))
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

      # Helper: safe subset+pseudobulk that handles empty donor lists
      safe_pb <- function(donor_list, mode) {
        if (length(donor_list) == 0) return(NULL)
        cells_in <- colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% donor_list]
        if (length(cells_in) == 0) return(NULL)
        so_sub <- subset(so_ct, cells = cells_in)
        make_pb(so_sub, donor_list, mode)
      }
      for (mode in c("mean", "sum")) {
        pb_t <- safe_pb(train_donors,   mode)
        pb_v <- safe_pb(val_donors,     mode)
        pb_h <- safe_pb(heldout_donors, mode)
        if (!is.null(pb_t)) write_mtx_with_names(pb_t, file.path(out_dir, paste0("train_pb_",   mode)))
        if (!is.null(pb_v)) write_mtx_with_names(pb_v, file.path(out_dir, paste0("val_pb_",     mode)))
        if (!is.null(pb_h)) write_mtx_with_names(pb_h, file.path(out_dir, paste0("heldout_pb_", mode)))
      }

      so_train <- subset(so_ct, cells = colnames(so_ct)[so_ct@meta.data[[DONOR_COL]] %in% train_donors])
      so_train <- FindVariableFeatures(so_train, selection.method = "vst",
                                        nfeatures = N_HVG, verbose = FALSE)
      hvg <- VariableFeatures(so_train)
      write.csv(data.frame(gene = hvg), file.path(out_dir, "HVG.csv"), row.names = FALSE)

      pb_train_mean <- make_pb(so_train, train_donors, "mean")
      tm_for_fit <- tm %>% select(donor_id, diagnosis, age, sex, batch, n_cells)
      ts <- tryCatch(
        fit_tstats(pb_train_mean[hvg, , drop = FALSE], tm_for_fit, "age", "sex", "batch"),
        error = function(e) {
          logmsg("    fit_tstats FAILED: ", conditionMessage(e), " - retrying ~Dx only")
          tryCatch({
            simple_meta <- tm_for_fit
            simple_meta$age <- NA; simple_meta$sex <- NA; simple_meta$batch <- NA
            fit_tstats(pb_train_mean[hvg, , drop = FALSE], simple_meta, NA, NA, NA)
          }, error = function(e2) NULL)
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

logmsg("Done. v8 mixed pseudobulk written to: ", OUT_ROOT)
