# ==============================================================================
# extract_pseudobulk_v5_compartment.R
# ------------------------------------------------------------------------------
# Step 0–2 仕様 (2026-04-28 改訂版)
#
# Step 0  - 評価設計
#   * Pappalardo cohort をまるごと外部 hold-out テストとして取り置き
#   * 残り cohort (Heming / Ramesh / Touil / etc.) で donor-level GroupKFold (5)
#   * cohort-stratified: 各 fold に各 cohort・各 Dx をなるべく均等に割る
#   * 由来組織 (CSF/PBMC) と cell type を全 cell に付与
#
# Step 1  - cell type 別 pseudobulk (donor 単位)
#   * primary: 平均集約 (data slot, normalized log-expression を donor 平均)
#   * sensitivity: 合計集約 (counts を donor 合計 → log1p(CPM))
#   * 初期特徴量: HVG ~3000 (train donor の細胞のみで計算)
#
# Step 2  - cell type 別 differential score (train donor only, donor-level)
#   * 共変量: Dx + age + sex + batch + log10(n_cells_per_donor_celltype)
#     ※ meta.data に存在する列のみ使用 (NA 比率高い列は自動で外す)
#   * lm() で各遺伝子に対し ~ Dx + 共変量、Dx の t統計量と p値を取得
#   * (cell type, tissue, fold) ごとに t-stat 表を保存
#   * candidate 集合の作成は Python 側 (Step 3) でやる方が一括管理しやすいので、
#     ここでは t-stat 表 + 上位 N (TOPN_PER_CELLTYPE) を出すまでにする
#
# 出力ルート:
#   MS_scRNA_GeneSelection_QUBO/data/pseudobulk_v5_compartment/
#     folds_assignment.csv                        # donor x cohort x label x fold ("heldout"含む)
#     covariate_availability.csv                  # 自動判定した利用可能列
#     <cell_type>/<tissue>/fold_<k>/
#       train_pb_mean.{mtx, _rows.csv, _cols.csv} # genes x donors (train, mean集約, primary)
#       val_pb_mean.{mtx, _rows.csv, _cols.csv}   # genes x donors (val/test fold)
#       heldout_pb_mean.{mtx, _rows.csv, _cols.csv}  # genes x Pappalardo donors
#       train_pb_sum.{mtx, ...}                   # 合計集約 (sensitivity)
#       val_pb_sum.{mtx, ...}
#       heldout_pb_sum.{mtx, ...}
#       train_meta.csv / val_meta.csv / heldout_meta.csv
#         # donor_id, diagnosis, cohort, compartment, age, sex, batch, n_cells, set
#       HVG.csv                                   # train cell から計算した HVG (上位 ~3000)
#       tstats.csv                                # gene, t, pval, padj, log2FC, rank, top_topn
#       topN_genes.csv                            # tstats.csv 上位 TOPN_PER_CELLTYPE 抜粋
#
# 想定ランタイム: 30-90 分 (12GB Seurat, M3 Pro 想定)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# --- パラメータ -----------------------------------------------------------------
RDS_PATH <- "/Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/data/raw_data/so.GEX.share.Asada_with_compartment.rds"
OUT_BASE <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/data/pseudobulk_v5_compartment"

# rotate held 用: 呼び出し側で holdout_cohort_arg を上書きすれば別 cohort を holdout に出来る
if (!exists("holdout_cohort_arg")) {
  HOLDOUT_COHORTS <- c("PRJNA671484_MS_Tcell")  # = Pappalardo (デフォルト)
} else {
  HOLDOUT_COHORTS <- holdout_cohort_arg
}
# 出力先は holdout 名で分岐 (デフォルト Pappalardo は無印で現行と互換)
HOLDOUT_TAG <- if (identical(HOLDOUT_COHORTS, c("PRJNA671484_MS_Tcell")))
                  "" else paste0("_holdout_", paste(HOLDOUT_COHORTS, collapse = "_"))
OUT_ROOT <- if (HOLDOUT_TAG == "") OUT_BASE else paste0(OUT_BASE, HOLDOUT_TAG)
# 参考: prj 値 → cohort
#   osmzhlab_MS_ence_cov     = Heming      (9 HD / 9 MS)
#   PRJNA671484_MS_Tcell     = Pappalardo  (6 HD / 5 MS)  ← hold-out
#   PRJNA549712_MS_PBMC_UCSF = Ramesh      (3 HD /14 MS)
#   PRJNA979258_cryoCSF      = Touil       (4 HD / 0 MS)
N_FOLDS <- 5
SEED <- 42

# meta.data の列名 (so のラベリングに合わせる; 違っていたら .rds 確認後に修正)
DX_COL          <- "Dx"
DONOR_COL       <- "patient_id"
COHORT_COL      <- "prj"
COMPARTMENT_COL <- "compartment"
CELLTYPE_COL    <- "predicted.celltype.l2"
# 共変量候補 (存在チェック後に自動採用)
COV_AGE_COLS   <- c("age", "Age", "age_years")
COV_SEX_COLS   <- c("sex", "Sex", "gender", "Gender")
COV_BATCH_COLS <- c("batch", "Batch", "library", "donor_batch", "tenx_run", "prj")
# prj が無い batch 列のフォールバックとして最後に入っている (= cohort 効果を吸収)

# 細胞型グループ (predicted.celltype.l2 のラベルに合わせる; 実ラベルは step1 で表示)
# v5.1 (2026-04-29): 4 → 8 cell type に拡張 (NK, DC, dnT, gdT 追加)
# 除外したもの:
#   Eryth/Platelet (contamination), HSPC/ILC (希少), NK Proliferating/ASDC/CD8 Proliferating (希少)
CELL_TYPES <- list(
  B      = c("B naive", "B memory", "B intermediate", "Plasmablast"),
  Mono   = c("CD14 Mono", "CD16 Mono"),
  CD4_T  = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg"),
  CD8_T  = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "MAIT"),
  NK     = c("NK", "NK_CD56bright"),
  DC     = c("cDC1", "cDC2", "pDC"),
  dnT    = c("dnT"),
  gdT    = c("gdT")
)
TISSUES <- c("CSF", "PBMC", "ALL")

N_HVG              <- 3000   # train cell から計算する HVG 数 (初期プール)
TOPN_PER_CELLTYPE  <- 100    # cell type ごとに t-stat 上位何遺伝子を出すか
# (Step 3 で union を取って候補集合とする)

# v5.1: 既存の (cell_type, tissue, fold) 出力をスキップする (新規 cell type のみ追加で計算したい時用)
SKIP_EXISTING <- TRUE

set.seed(SEED)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ==============================================================================
# step1: 読み込み & meta.data 監査
# ==============================================================================
logmsg("Loading Seurat object: ", RDS_PATH)
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"
logmsg("Object summary (raw):")
print(so)

# --- メモリ削減: 不要な assay/layer/reduction を削除 ---------------------------
# Seurat v5 (layers=) と v4 (slot=) の両対応
logmsg("DietSeurat: drop scale.data, prediction.score.* assays, all reductions")
so <- tryCatch(
  DietSeurat(so, layers = c("counts", "data"), assays = "RNA", dimreducs = NULL),
  error = function(e1) tryCatch(
    DietSeurat(so, slot   = c("counts", "data"), assays = "RNA", dimreducs = NULL),
    error = function(e2) {
      warning("DietSeurat failed (", conditionMessage(e1), " / ", conditionMessage(e2),
              "), proceeding without diet")
      so
    }
  )
)
gc()
logmsg("Object summary (after diet):")
print(so)

needed <- c(DX_COL, DONOR_COL, COHORT_COL, COMPARTMENT_COL, CELLTYPE_COL)
missing_cols <- setdiff(needed, colnames(so@meta.data))
if (length(missing_cols)) {
  stop("meta.data に必要な列がありません: ", paste(missing_cols, collapse = ", "),
       "\n  現状の列: ", paste(colnames(so@meta.data), collapse = ", "))
}

logmsg("meta.data columns (", ncol(so@meta.data), "):")
print(colnames(so@meta.data))

logmsg("Dx (cell-level):")
print(table(so@meta.data[[DX_COL]], useNA = "ifany"))
logmsg("compartment (cell-level):")
print(table(so@meta.data[[COMPARTMENT_COL]], useNA = "ifany"))
logmsg("predicted.celltype.l2 (top 30):")
print(head(sort(table(so@meta.data[[CELLTYPE_COL]]), decreasing = TRUE), 30))

donor_tbl <- so@meta.data %>%
  distinct(donor   = .data[[DONOR_COL]],
           label   = .data[[DX_COL]],
           cohort  = .data[[COHORT_COL]]) %>%
  filter(!is.na(donor), !is.na(label))
logmsg("Donors per cohort × Dx:")
print(donor_tbl %>% count(cohort, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))

# --- 共変量列の自動判定 -------------------------------------------------------
pick_first_present <- function(candidates, df) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NA_character_
}
age_col   <- pick_first_present(COV_AGE_COLS,   so@meta.data)
sex_col   <- pick_first_present(COV_SEX_COLS,   so@meta.data)
batch_col <- pick_first_present(COV_BATCH_COLS, so@meta.data)

# NA 比率を計算（donor 単位で見る）
donor_meta_full <- so@meta.data %>%
  group_by(.data[[DONOR_COL]]) %>%
  summarise(
    age   = if (!is.na(age_col))   first(.data[[age_col]])   else NA,
    sex   = if (!is.na(sex_col))   first(.data[[sex_col]])   else NA,
    batch = if (!is.na(batch_col)) first(.data[[batch_col]]) else NA,
    .groups = "drop"
  ) %>%
  rename(donor = !!DONOR_COL)

cov_avail <- tibble(
  covariate = c("age", "sex", "batch"),
  src_col   = c(age_col, sex_col, batch_col),
  na_rate_donor = c(
    mean(is.na(donor_meta_full$age)),
    mean(is.na(donor_meta_full$sex)),
    mean(is.na(donor_meta_full$batch))
  )
) %>%
  mutate(usable = !is.na(src_col) & na_rate_donor < 0.5)

logmsg("Covariate availability (donor-level NA rate):")
print(cov_avail)
write.csv(cov_avail, file.path(OUT_ROOT, "covariate_availability.csv"), row.names = FALSE)

USE_AGE   <- cov_avail$usable[cov_avail$covariate == "age"]
USE_SEX   <- cov_avail$usable[cov_avail$covariate == "sex"]
USE_BATCH <- cov_avail$usable[cov_avail$covariate == "batch"]

# ==============================================================================
# step2: fold 割当 (Pappalardo hold-out + 残 cohort で cohort-stratified 5-fold)
# ==============================================================================
heldout_donors <- donor_tbl$donor[donor_tbl$cohort %in% HOLDOUT_COHORTS]
cv_donors <- donor_tbl %>% filter(!cohort %in% HOLDOUT_COHORTS)

logmsg(sprintf("Hold-out cohorts: %s -> %d donors",
               paste(HOLDOUT_COHORTS, collapse=","), length(heldout_donors)))
logmsg(sprintf("CV donors (excluding hold-out): %d", nrow(cv_donors)))
print(cv_donors %>% count(cohort, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))

# cohort × label の各層内で donor をシャッフルして fold を割る (cohort-stratified)
assign_fold_balanced <- function(df, n_folds = 5, seed = 42) {
  set.seed(seed)
  df %>%
    group_by(cohort, label) %>%
    mutate(.shuffle = sample(seq_len(n()))) %>%
    arrange(cohort, label, .shuffle) %>%
    mutate(fold = ((row_number() - 1) %% n_folds) + 1) %>%
    ungroup() %>%
    select(-.shuffle)
}
cv_assign <- assign_fold_balanced(cv_donors, n_folds = N_FOLDS, seed = SEED)

folds_tbl <- bind_rows(
  cv_assign %>% mutate(set = "cv", fold = as.character(fold)),
  donor_tbl %>% filter(donor %in% heldout_donors) %>%
    mutate(set = "heldout", fold = "heldout")
) %>%
  arrange(set, fold, cohort, donor)

write.csv(folds_tbl, file.path(OUT_ROOT, "folds_assignment.csv"), row.names = FALSE)
logmsg("Fold assignment summary:")
print(folds_tbl %>% count(set, fold, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))
print(folds_tbl %>% count(set, fold, cohort) %>%
        pivot_wider(names_from = cohort, values_from = n, values_fill = 0))

cv_fold_ids <- as.character(seq_len(N_FOLDS))

# ==============================================================================
# step3: pseudobulk ヘルパ
# ==============================================================================
make_donor_indicator <- function(donors) {
  fac <- factor(donors)
  list(
    fac = fac,
    ind = sparseMatrix(i = as.integer(fac), j = seq_along(fac), x = 1,
                       dims = c(nlevels(fac), length(fac)))
  )
}

# 平均集約 (data slot, log-normalized expression を donor 平均) — primary
pseudobulk_mean <- function(seu_subset, donor_col) {
  expr <- GetAssayData(seu_subset, slot = "data")  # genes x cells
  donors <- as.character(seu_subset@meta.data[[donor_col]])
  d <- make_donor_indicator(donors)
  pb_sum <- expr %*% Matrix::t(d$ind)              # genes x donors (sum)
  n_cells <- as.numeric(table(d$fac))
  pb_mean <- t(t(pb_sum) / n_cells)                # divide each donor column by n
  colnames(pb_mean) <- levels(d$fac)
  pb_mean
}

# 合計集約 (counts を donor 合計 → log1p(CPM)) — sensitivity
pseudobulk_sum <- function(seu_subset, donor_col) {
  cnt <- GetAssayData(seu_subset, slot = "counts") # genes x cells
  donors <- as.character(seu_subset@meta.data[[donor_col]])
  d <- make_donor_indicator(donors)
  pb_sum <- cnt %*% Matrix::t(d$ind)               # genes x donors (raw counts sum)
  libsize <- Matrix::colSums(pb_sum)
  libsize[libsize == 0] <- 1
  pb_norm <- t(t(pb_sum) * (1e6 / libsize))        # CPM
  pb_norm@x <- log1p(pb_norm@x)
  colnames(pb_norm) <- levels(d$fac)
  pb_norm
}

# n_cells per donor (集約前の細胞数)
n_cells_per_donor <- function(seu_subset, donor_col) {
  table(as.character(seu_subset@meta.data[[donor_col]]))
}

# 共変量 + Dx で donor-level lm → Dx の t統計量
compute_donor_tstats <- function(pb_train, train_meta, use_age, use_sex, use_batch) {
  X <- as.matrix(pb_train)                         # genes x donors
  if (ncol(X) < 4 || nrow(X) == 0) return(NULL)
  meta <- train_meta %>% mutate(log10_n = log10(pmax(n_cells, 1)))
  meta$Dx_bin <- ifelse(meta$diagnosis == "MS", 1L, 0L)

  rhs <- c("Dx_bin", "log10_n")
  if (use_age   && all(!is.na(meta$age)))   rhs <- c(rhs, "age")
  if (use_sex   && all(!is.na(meta$sex))   && length(unique(meta$sex))   > 1) rhs <- c(rhs, "sex")
  if (use_batch && all(!is.na(meta$batch)) && length(unique(meta$batch)) > 1) rhs <- c(rhs, "batch")
  fmla_rhs <- paste(rhs, collapse = " + ")

  # gene 1 つずつ lm するが、行列演算で一括化: y ~ X 1モデルを使う
  # (gene 数が多いので apply で回す。1 万遺伝子 × ~30 donor で数十秒のはず)
  res <- vapply(seq_len(nrow(X)), function(i) {
    y <- X[i, ]
    if (sd(y) == 0) return(c(t = 0, p = 1, beta = 0, se = 0))
    df <- meta
    df$y <- y
    fit <- tryCatch(lm(as.formula(paste("y ~", fmla_rhs)), data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(c(t = 0, p = 1, beta = 0, se = 0))
    s <- summary(fit)$coefficients
    if (!"Dx_bin" %in% rownames(s)) return(c(t = 0, p = 1, beta = 0, se = 0))
    c(t = s["Dx_bin", "t value"],
      p = s["Dx_bin", "Pr(>|t|)"],
      beta = s["Dx_bin", "Estimate"],
      se = s["Dx_bin", "Std. Error"])
  }, numeric(4))

  data.frame(
    gene = rownames(X),
    t     = res["t",],
    pval  = res["p",],
    beta  = res["beta",],
    se    = res["se",],
    log2FC = res["beta",] * log2(exp(1)),  # data slot は ln-normalized なので log2 換算
    stringsAsFactors = FALSE
  ) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(desc(abs(t))) %>%
    mutate(rank = row_number(),
           top_topn = rank <= TOPN_PER_CELLTYPE)
}

# ==============================================================================
# step4: メインループ (cell_type × tissue × fold)
# ==============================================================================
all_l2 <- unique(as.character(so@meta.data[[CELLTYPE_COL]]))
n_done <- 0; n_skip <- 0

write_pb <- function(pb, prefix, dir) {
  writeMM(pb, file.path(dir, paste0(prefix, ".mtx")))
  write.csv(data.frame(gene  = rownames(pb)), file.path(dir, paste0(prefix, "_rows.csv")), row.names = FALSE)
  write.csv(data.frame(donor = colnames(pb)), file.path(dir, paste0(prefix, "_cols.csv")), row.names = FALSE)
}

build_meta <- function(seu_subset, set_label) {
  m <- seu_subset@meta.data %>%
    distinct(donor_id   = .data[[DONOR_COL]],
             diagnosis  = .data[[DX_COL]],
             cohort     = .data[[COHORT_COL]],
             compartment = .data[[COMPARTMENT_COL]])
  if (!is.na(age_col))   m$age   <- seu_subset@meta.data[[age_col]][match(m$donor_id, seu_subset@meta.data[[DONOR_COL]])]   else m$age   <- NA
  if (!is.na(sex_col))   m$sex   <- as.character(seu_subset@meta.data[[sex_col]][match(m$donor_id, seu_subset@meta.data[[DONOR_COL]])])   else m$sex   <- NA
  if (!is.na(batch_col)) m$batch <- as.character(seu_subset@meta.data[[batch_col]][match(m$donor_id, seu_subset@meta.data[[DONOR_COL]])]) else m$batch <- NA
  nc <- n_cells_per_donor(seu_subset, DONOR_COL)
  m$n_cells <- as.integer(nc[m$donor_id])
  m$set <- set_label
  m %>% arrange(donor_id)
}

for (ct_name in names(CELL_TYPES)) {
  ct_labels <- intersect(CELL_TYPES[[ct_name]], all_l2)
  if (length(ct_labels) == 0) {
    logmsg(sprintf("[skip] %s: l2 ラベル一致なし (希望=%s)",
                   ct_name, paste(CELL_TYPES[[ct_name]], collapse=",")))
    next
  }
  logmsg(sprintf("=== cell_type=%s (l2=%s) ===", ct_name, paste(ct_labels, collapse="/")))
  ct_cells <- rownames(so@meta.data)[so@meta.data[[CELLTYPE_COL]] %in% ct_labels]
  so_ct <- subset(so, cells = ct_cells)

  for (tissue in TISSUES) {
    if (tissue == "ALL") {
      use_cells <- colnames(so_ct)
    } else {
      use_cells <- colnames(so_ct)[so_ct@meta.data[[COMPARTMENT_COL]] == tissue &
                                     !is.na(so_ct@meta.data[[COMPARTMENT_COL]])]
    }
    if (length(use_cells) < 100) {
      logmsg(sprintf("[skip] %s/%s: cells=%d (<100)", ct_name, tissue, length(use_cells)))
      n_skip <- n_skip + 1; next
    }
    # SKIP_EXISTING: この (cell_type, tissue) の全 fold が既に揃っているなら丸ごとスキップ
    if (SKIP_EXISTING) {
      done_dirs <- vapply(cv_fold_ids, function(k) {
        file.exists(file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k), "tstats.csv"))
      }, logical(1))
      if (all(done_dirs)) {
        logmsg(sprintf("[skip exist] %s/%s: all folds already extracted", ct_name, tissue))
        n_skip <- n_skip + length(cv_fold_ids); next
      }
    }

    so_t <- subset(so_ct, cells = use_cells)
    logmsg(sprintf("--- %s / %s : %d cells, %d donors ---",
                   ct_name, tissue, ncol(so_t),
                   length(unique(so_t@meta.data[[DONOR_COL]]))))

    # tissue 単位で hold-out 部分は全 fold 共通なので 1 度だけ計算
    heldout_cells <- colnames(so_t)[so_t@meta.data[[DONOR_COL]] %in% heldout_donors]
    so_held <- if (length(heldout_cells) >= 10) subset(so_t, cells = heldout_cells) else NULL
    pb_held_mean <- if (!is.null(so_held)) pseudobulk_mean(so_held, DONOR_COL) else NULL
    pb_held_sum  <- if (!is.null(so_held)) pseudobulk_sum(so_held,  DONOR_COL) else NULL
    meta_held <- if (!is.null(so_held)) build_meta(so_held, "heldout") else NULL

    for (k in cv_fold_ids) {
      # fold 単位でも既存スキップ
      fold_out_dir <- file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k))
      if (SKIP_EXISTING && file.exists(file.path(fold_out_dir, "tstats.csv"))) {
        logmsg(sprintf("[skip exist] %s/%s/fold_%s", ct_name, tissue, k))
        n_skip <- n_skip + 1; next
      }
      train_donors <- folds_tbl$donor[folds_tbl$set == "cv" & folds_tbl$fold != k]
      val_donors   <- folds_tbl$donor[folds_tbl$set == "cv" & folds_tbl$fold == k]

      train_cells <- colnames(so_t)[so_t@meta.data[[DONOR_COL]] %in% train_donors]
      val_cells   <- colnames(so_t)[so_t@meta.data[[DONOR_COL]] %in% val_donors]
      if (length(train_cells) < 50 || length(val_cells) < 10) {
        logmsg(sprintf("[skip] %s/%s/fold_%s: train_cells=%d val_cells=%d",
                       ct_name, tissue, k, length(train_cells), length(val_cells)))
        n_skip <- n_skip + 1; next
      }
      val_dx <- unique(so_t@meta.data[[DX_COL]][so_t@meta.data[[DONOR_COL]] %in% val_donors])
      if (length(val_dx) < 2) {
        logmsg(sprintf("[skip] %s/%s/fold_%s: val に %s しか無い",
                       ct_name, tissue, k, paste(val_dx, collapse="/")))
        n_skip <- n_skip + 1; next
      }

      so_train <- subset(so_t, cells = train_cells)
      so_val   <- subset(so_t, cells = val_cells)

      # ---- HVG (train cell only) ----
      so_train <- FindVariableFeatures(so_train, selection.method = "vst",
                                       nfeatures = N_HVG, verbose = FALSE)
      hvg <- VariableFeatures(so_train)
      if (length(hvg) < 500) {
        logmsg(sprintf("[skip] %s/%s/fold_%s: HVG=%d (<500)",
                       ct_name, tissue, k, length(hvg))); n_skip <- n_skip + 1; next
      }

      # ---- pseudobulk (train, val) ----
      pb_train_mean <- pseudobulk_mean(so_train, DONOR_COL)[hvg, , drop = FALSE]
      pb_val_mean   <- pseudobulk_mean(so_val,   DONOR_COL)
      pb_val_mean   <- pb_val_mean[intersect(hvg, rownames(pb_val_mean)), , drop = FALSE]
      pb_train_sum <- pseudobulk_sum(so_train, DONOR_COL)[hvg, , drop = FALSE]
      pb_val_sum   <- pseudobulk_sum(so_val,   DONOR_COL)
      pb_val_sum   <- pb_val_sum[intersect(hvg, rownames(pb_val_sum)), , drop = FALSE]

      # heldout は tissue 単位で計算済み → HVG で部分集合
      pb_held_mean_k <- if (!is.null(pb_held_mean)) pb_held_mean[intersect(hvg, rownames(pb_held_mean)), , drop = FALSE] else NULL
      pb_held_sum_k  <- if (!is.null(pb_held_sum))  pb_held_sum[intersect(hvg, rownames(pb_held_sum)),  , drop = FALSE]  else NULL

      # ---- meta ----
      meta_train <- build_meta(so_train, "train")
      meta_val   <- build_meta(so_val,   "val")

      # 列順整合
      pb_train_mean <- pb_train_mean[, meta_train$donor_id, drop = FALSE]
      pb_val_mean   <- pb_val_mean[,   meta_val$donor_id,   drop = FALSE]
      pb_train_sum  <- pb_train_sum[,  meta_train$donor_id, drop = FALSE]
      pb_val_sum    <- pb_val_sum[,    meta_val$donor_id,   drop = FALSE]
      if (!is.null(pb_held_mean_k)) {
        meta_held_aligned <- meta_held %>% filter(donor_id %in% colnames(pb_held_mean_k))
        pb_held_mean_k <- pb_held_mean_k[, meta_held_aligned$donor_id, drop = FALSE]
        pb_held_sum_k  <- pb_held_sum_k[,  meta_held_aligned$donor_id, drop = FALSE]
      } else {
        meta_held_aligned <- NULL
      }

      # ---- t-stat (train donor-level pseudobulk, mean primary) ----
      tstats_df <- compute_donor_tstats(pb_train_mean, meta_train,
                                        USE_AGE, USE_SEX, USE_BATCH)
      if (is.null(tstats_df) || nrow(tstats_df) == 0) {
        logmsg(sprintf("[skip] %s/%s/fold_%s: tstats 失敗", ct_name, tissue, k))
        n_skip <- n_skip + 1; next
      }
      topN_df <- tstats_df %>% filter(top_topn) %>% select(gene, t, pval, padj, log2FC, rank)

      # ---- 書き出し ----
      out_dir <- file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      write_pb(pb_train_mean, "train_pb_mean", out_dir)
      write_pb(pb_val_mean,   "val_pb_mean",   out_dir)
      write_pb(pb_train_sum,  "train_pb_sum",  out_dir)
      write_pb(pb_val_sum,    "val_pb_sum",    out_dir)
      if (!is.null(pb_held_mean_k) && ncol(pb_held_mean_k) > 0) {
        write_pb(pb_held_mean_k, "heldout_pb_mean", out_dir)
        write_pb(pb_held_sum_k,  "heldout_pb_sum",  out_dir)
        write.csv(meta_held_aligned, file.path(out_dir, "heldout_meta.csv"), row.names = FALSE)
      }
      write.csv(meta_train, file.path(out_dir, "train_meta.csv"), row.names = FALSE)
      write.csv(meta_val,   file.path(out_dir, "val_meta.csv"),   row.names = FALSE)
      write.csv(data.frame(gene = hvg), file.path(out_dir, "HVG.csv"), row.names = FALSE)
      write.csv(tstats_df, file.path(out_dir, "tstats.csv"), row.names = FALSE)
      write.csv(topN_df,   file.path(out_dir, "topN_genes.csv"), row.names = FALSE)

      n_done <- n_done + 1
      logmsg(sprintf("[ok]   %s/%s/fold_%s: train=%dd val=%dd held=%dd  HVG=%d  topN=%d",
                     ct_name, tissue, k,
                     ncol(pb_train_mean), ncol(pb_val_mean),
                     ifelse(is.null(pb_held_mean_k), 0, ncol(pb_held_mean_k)),
                     length(hvg), nrow(topN_df)))

      # 中間オブジェクトを明示的に削除して GC
      rm(so_train, so_val, pb_train_mean, pb_val_mean, pb_train_sum, pb_val_sum,
         pb_held_mean_k, pb_held_sum_k, tstats_df, topN_df, hvg, meta_train, meta_val)
      gc(verbose = FALSE)
    }
    rm(so_t, so_held, pb_held_mean, pb_held_sum, meta_held)
    gc(verbose = FALSE)
  }
  rm(so_ct)
  gc(verbose = FALSE)
}

logmsg(sprintf("DONE. ok=%d skip=%d  out=%s", n_done, n_skip, OUT_ROOT))
