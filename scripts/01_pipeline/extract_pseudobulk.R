# ==============================================================================
# extract_pseudobulk.R
# ------------------------------------------------------------------------------
# Per-cell-type donor-level pseudobulk extraction (Steps 0-2 of the qubofs
# pipeline). Produces train/val/heldout pseudobulk matrices (log-norm mean,
# raw/ambient-corrected count sum, and legacy log1p-CPM sum), organised by
# (cell_type, tissue, fold). It also writes optional legacy/QC t-statistic tables
# (tstats.csv / topN_genes.csv) that are NOT the manuscript feature ranking
# (see Step 2). The canonical DE ranking is computed downstream in 02_deg.
#
# Internal donor labels use "MS" and "HD" because the upstream Seurat object
# uses these factor levels; the manuscript reports them as "MS" and "control".
#
# Step 0  - Evaluation design
#   * Pappalardo cohort held out as external test cohort.
#   * Remaining cohorts (Heming / Ramesh / Touil) are assigned to five
#     donor-level folds using a deterministic cohort-by-diagnosis balanced split.
#   * Tissue/compartment and cell type annotated for every cell (CSF used here).
#
# Step 1  - Per-cell-type donor-level pseudobulk
#   * primary aggregation (mean): donor-wise mean of normalised log-expression
#       (data slot) -> *_pb_mean.mtx; used for classification, redundancy, stability
#   * count-scale aggregation: donor-wise SUM of SoupX ambient-corrected
#       count-scale values (no normalisation) -> *_pb_counts.mtx; rounded to
#       integers downstream as input for count-based pseudobulk DE
#       (edgeR primary; limma/DESeq2 sensitivity) in 02_deg
#   * normalised-sum aggregation: donor-wise sum of counts -> log1p(CPM)
#       -> *_pb_sum.mtx; retained for backward compatibility / sensitivity only.
#       NOTE: *_pb_sum.mtx is NOT count-scale; do not feed it to count-based DE.
#   * initial features: HVG ~3000 (computed from train-donor cells only)
#
# Step 2  - OPTIONAL legacy / QC ranking (NOT the manuscript ranking)
#   * lm() per gene on the log-normalised pseudobulk mean: ~ Dx + covariates;
#     Dx t-statistic and p-value retained, saved as tstats.csv / topN_genes.csv.
#   * These are retained for backward compatibility / QC ONLY. The canonical
#     differential-expression ranking used in the manuscript is the edgeR test statistic
#     |z| on *_pb_counts.mtx (weighted by cohort-consistency C_i), computed
#     downstream in 02_deg / 03_selection — NOT from these files.
#
# Output layout (relative to OUT_BASE):
#   <pseudobulk_subdir>/
#     folds_assignment.csv
#     covariate_availability.csv
#     <cell_type>/<tissue>/fold_<k>/
#       train_pb_mean.{mtx, _rows.csv, _cols.csv}
#       val_pb_mean.{mtx, _rows.csv, _cols.csv}
#       heldout_pb_mean.{mtx, _rows.csv, _cols.csv}
#       train_pb_counts.{mtx, ...}  val_pb_counts.{...}  heldout_pb_counts.{...}  (count-scale SoupX-corrected sums)
#       train_pb_sum.{mtx, ...}  val_pb_sum.{...}  heldout_pb_sum.{...}  (log1p CPM, legacy)
#       train_meta.csv / val_meta.csv / heldout_meta.csv
#       HVG.csv
#       tstats.csv       (legacy/QC only — see Step 2; not the manuscript ranking)
#       topN_genes.csv   (legacy/QC only — derived from tstats.csv)
#
# Typical runtime: 30-90 minutes (12 GB Seurat object, M3 Pro).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# --- Parameters ----------------------------------------------------------------
# Paths are taken from environment variables if set, with placeholders that the
# user must override before running. Set QUBOFS_SEURAT_RDS and
# QUBOFS_OUT_BASE before sourcing this script.
RDS_PATH <- Sys.getenv("QUBOFS_SEURAT_RDS",
                       unset = "<path/to>/integrated_scRNAseq.rds")
OUT_BASE <- Sys.getenv("QUBOFS_OUT_BASE",
                       unset = file.path("data", "pseudobulk"))

# rotate hold-out: override holdout_cohort_arg from the caller to choose a
# different cohort as the external test set.
if (!exists("holdout_cohort_arg")) {
  HOLDOUT_COHORTS <- c("PRJNA671484_MS_Tcell")  # = Pappalardo (default)
} else {
  HOLDOUT_COHORTS <- holdout_cohort_arg
}
# Output directory suffix encodes the hold-out cohort (default Pappalardo:
# no suffix, for backward compatibility).
HOLDOUT_TAG <- if (identical(HOLDOUT_COHORTS, c("PRJNA671484_MS_Tcell")))
                  "" else paste0("_holdout_", paste(HOLDOUT_COHORTS, collapse = "_"))
OUT_ROOT <- if (HOLDOUT_TAG == "") OUT_BASE else paste0(OUT_BASE, HOLDOUT_TAG)
# Reference: prj value -> cohort (donor counts; internal labels MS / HD).
#   osmzhlab_MS_ence_cov     = Heming      (9 control / 9 MS)
#   PRJNA671484_MS_Tcell     = Pappalardo  (6 control / 5 MS)   <- default hold-out
#   PRJNA549712_MS_PBMC_UCSF = Ramesh      (3 control / 14 MS)
#   PRJNA979258_cryoCSF      = Touil       (4 control / 0 MS)
N_FOLDS <- 5
SEED <- 42

# meta.data column names (match the labelling in 'so'; if different, fix after checking the .rds)
DX_COL          <- "Dx"
DONOR_COL       <- "patient_id"
COHORT_COL      <- "prj"
COMPARTMENT_COL <- "compartment"
CELLTYPE_COL    <- "predicted.celltype.l2"
# candidate covariates (auto-adopted after an existence check)
COV_AGE_COLS   <- c("age", "Age", "age_years")
COV_SEX_COLS   <- c("sex", "Sex", "gender", "Gender")
COV_BATCH_COLS <- c("batch", "Batch", "library", "donor_batch", "tenx_run", "prj")
# placed last as a fallback batch column when 'prj' is absent (absorbs cohort effects)

# cell-type groups (match predicted.celltype.l2 labels; actual labels are printed in step1)
# v5.1 (2026-04-29): expanded from 4 to 8 cell types (added NK, DC, dnT, gdT)
# excluded:
#   Eryth/Platelet (contamination), HSPC/ILC (rare), NK Proliferating/ASDC/CD8 Proliferating (rare)
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
# Compartment(s) to build pseudobulk for. Default CSF (primary analysis); set
# QUBOFS_TISSUES=PBMC (or "CSF,PBMC") to build the peripheral-blood/PBMC
# benchmark. PBMC is available only in Pappalardo and Ramesh (Heming/Touil are
# CSF-only), so the blood benchmark covers those two cohorts.
TISSUES <- strsplit(Sys.getenv("QUBOFS_TISSUES", "CSF"), ",")[[1]]

N_HVG              <- 3000   # number of HVGs computed from training cells (initial pool)

# ---------------------------------------------------------------------------
# Pure-technical gene removal BEFORE highly variable gene selection.
#
# Non-informative high-abundance genes (mitochondrial / ribosomal-protein /
# heat-shock / classical housekeeping / nuclear lncRNA / small RNA) otherwise
# occupy a large fraction of the HVG slots without carrying disease signal, so
# they are removed from the feature universe before vst HVG selection. This is
# standard scRNA-seq feature-universe practice (Heumos et al. Nat Rev Genet
# 2023; Luecken & Theis Mol Syst Biol 2019).
#
# PURE_TECH_PATTERN is kept character-for-character identical to HK_PATTERN in
# 03_selection/qubo_pipeline.py. It does NOT include V(D)J segments: clonotype-
# driven immunoglobulin/TCR V(D)J genes, together with the cell-type-aware
# detection/specificity criteria, are applied later at the selection stage
# (03_selection), because those filters define the candidate biomarker space
# rather than the pseudobulk/HVG representation itself. ER chaperones HSPA5 /
# HSP90B1 are retained because they are secretory-pathway-associated features and
# should not be removed as generic cytosolic stress genes.
PURE_TECH_PATTERN <- paste0(
  "^(MT-|MTRNR|MTATP|MTND|",                 # mitochondrial
  "RPL[0-9]|RPLP|RPS[0-9]|RPSA|MRPL|MRPS|",  # ribosomal proteins (incl. RPLP0/1/2 stalk, RPSA)
  "HSP[A0-9]|HSPB|HSPA|HSPD|",               # heat shock
  "FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|",      # classical housekeeping
  "MALAT1$|NEAT1$|XIST$|TSIX$|",             # nuclear lncRNA + X-inactivation
  "AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|",        # uncharacterized / pseudogene loci
  "MIR[0-9]|RNU[0-9]|SNORA|SNORD)"           # small RNAs (poly-A unreliable)
)
BIOLOGY_RETAIN <- c("HSPA5", "HSP90B1")      # ER chaperones kept despite HSP* pattern

# NOTE: the cell-type-aware (detection / specificity) and V(D)J filters are NOT
# applied here; they are applied once, downstream, in 03_selection (single
# source of truth), identically to all feature-selection methods. Only the
# pure-technical removal above is done at this stage.
TOPN_PER_CELLTYPE  <- 100    # number of top t-stat genes to output per cell type
# (a union is taken in Step 3 to form the candidate set)

# v5.1: skip existing (cell_type, tissue, fold) outputs (for adding only new cell types)
# set env QUBOFS_SKIP_EXISTING=0 to force recomputation of all folds (overwrite). After a
# filter change, specify 0 so the cache need not be cleared by hand.
SKIP_EXISTING <- (Sys.getenv("QUBOFS_SKIP_EXISTING", "1") == "1")

set.seed(SEED)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ==============================================================================
# step1: load & audit meta.data
# ==============================================================================
logmsg("Loading Seurat object: ", RDS_PATH)
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"
logmsg("Object summary (raw):")
print(so)

# --- memory reduction: drop unused assays/layers/reductions ---------------------------
# supports both Seurat v5 (layers=) and v4 (slot=)
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
  stop("required columns missing from meta.data: ", paste(missing_cols, collapse = ", "),
       "\n  current columns: ", paste(colnames(so@meta.data), collapse = ", "))
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

# --- automatic determination of covariate columns -------------------------------------------------------
pick_first_present <- function(candidates, df) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  NA_character_
}
age_col   <- pick_first_present(COV_AGE_COLS,   so@meta.data)
sex_col   <- pick_first_present(COV_SEX_COLS,   so@meta.data)
batch_col <- pick_first_present(COV_BATCH_COLS, so@meta.data)

# compute NA fraction (assessed per donor)
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
# step2: fold assignment (Pappalardo hold-out + remaining cohorts, cohort-stratified 5-fold)
# ==============================================================================
# Pooled-CV mode (QUBOFS_POOLED=1): no external cohort is held out. All donors of
# the selected tissue(s) are pooled and split into N_FOLDS donor-stratified CV
# folds; for each fold the held-out test set IS that CV fold (standard k-fold CV).
# Used for the PBMC generalisation benchmark (Pappalardo + Ramesh pooled), where a
# 2-cohort leave-one-cohort-out split would train on a single cohort.
POOLED <- (Sys.getenv("QUBOFS_POOLED", "0") == "1")
if (POOLED) {
  # Restrict the pooled CV to donors that actually have cells in the requested
  # tissue(s); otherwise donors lacking the tissue (e.g. the CSF-only cohorts
  # Heming/Touil under QUBOFS_TISSUES=PBMC) would be assigned to folds with no
  # evaluable cells, distorting the per-fold splits.
  donors_with_tissue <- so@meta.data %>%
    filter(.data[[COMPARTMENT_COL]] %in% TISSUES) %>%
    distinct(donor = .data[[DONOR_COL]]) %>%
    pull(donor)
  heldout_donors <- character(0)
  cv_donors <- donor_tbl %>% filter(donor %in% donors_with_tissue)
  logmsg(sprintf("[pooled] QUBOFS_POOLED=1: no external cohort held out; N-fold CV over %d donors with tissue(s) %s; per-fold held-out test = the CV fold.",
                 nrow(cv_donors), paste(TISSUES, collapse = ",")))
} else {
  heldout_donors <- donor_tbl$donor[donor_tbl$cohort %in% HOLDOUT_COHORTS]
  cv_donors <- donor_tbl %>% filter(!cohort %in% HOLDOUT_COHORTS)
}

logmsg(sprintf("Hold-out cohorts: %s -> %d donors",
               paste(HOLDOUT_COHORTS, collapse=","), length(heldout_donors)))
logmsg(sprintf("CV donors (excluding hold-out): %d", nrow(cv_donors)))
print(cv_donors %>% count(cohort, label) %>%
        pivot_wider(names_from = label, values_from = n, values_fill = 0))

# shuffle donors within each cohort x label stratum and assign folds (cohort-stratified)
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
# step3: pseudobulk helpers
# ==============================================================================
make_donor_indicator <- function(donors) {
  fac <- factor(donors)
  list(
    fac = fac,
    ind = sparseMatrix(i = as.integer(fac), j = seq_along(fac), x = 1,
                       dims = c(nlevels(fac), length(fac)))
  )
}

# Seurat v5 deprecates GetAssayData(..., slot=) as a hard error; use layer=.
# This helper works on both Seurat v5 (layer=) and v4 (slot=).
get_assay_matrix <- function(seu, name) {
  tryCatch(
    GetAssayData(seu, layer = name),
    error = function(e) GetAssayData(seu, slot = name)
  )
}

# mean aggregation (donor mean of log-normalized expression, data slot) -- primary
pseudobulk_mean <- function(seu_subset, donor_col) {
  expr <- get_assay_matrix(seu_subset, "data")  # genes x cells (log-normalized)
  donors <- as.character(seu_subset@meta.data[[donor_col]])
  d <- make_donor_indicator(donors)
  pb_sum <- expr %*% Matrix::t(d$ind)              # genes x donors (sum)
  n_cells <- as.numeric(table(d$fac))
  pb_mean <- t(t(pb_sum) / n_cells)                # divide each donor column by n
  colnames(pb_mean) <- levels(d$fac)
  pb_mean
}

# sum aggregation (donor sum of counts -> log1p(CPM)) -- sensitivity
pseudobulk_sum <- function(seu_subset, donor_col) {
  cnt <- get_assay_matrix(seu_subset, "counts") # genes x cells
  donors <- as.character(seu_subset@meta.data[[donor_col]])
  d <- make_donor_indicator(donors)
  pb_sum <- cnt %*% Matrix::t(d$ind)               # genes x donors (count-scale sum)
  libsize <- Matrix::colSums(pb_sum)
  libsize[libsize == 0] <- 1
  pb_norm <- t(t(pb_sum) * (1e6 / libsize))        # CPM
  pb_norm@x <- log1p(pb_norm@x)
  colnames(pb_norm) <- levels(d$fac)
  pb_norm
}

# count-scale sum (donor sum of SoupX-corrected counts, no normalisation) -- for edgeR (primary), limma/DESeq2 (sensitivity)
# Raw integer count sums per (donor x gene). This is the statistically correct
# input for count-based pseudobulk DE (negative-binomial GLMs). Saved as
# *_pb_counts.mtx, separate from the log1p(CPM) *_pb_sum.mtx above.
pseudobulk_count_sum <- function(seu_subset, donor_col) {
  cnt <- get_assay_matrix(seu_subset, "counts") # genes x cells (SoupX-corrected counts)
  donors <- as.character(seu_subset@meta.data[[donor_col]])
  d <- make_donor_indicator(donors)
  pb_counts <- cnt %*% Matrix::t(d$ind)            # genes x donors (count-scale sum)
  colnames(pb_counts) <- levels(d$fac)
  pb_counts
}

# n_cells per donor (cell count before aggregation)
n_cells_per_donor <- function(seu_subset, donor_col) {
  table(as.character(seu_subset@meta.data[[donor_col]]))
}

# donor-level lm with covariates + Dx -> t statistic for Dx
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

  # would lm each gene separately, but vectorise via matrix algebra: use one y ~ X model
  # (many genes, so loop with apply; ~10k genes x ~30 donors takes tens of seconds)
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
    log2FC = res["beta",] * log2(exp(1)),  # data slot is ln-normalized, so convert to log2
    stringsAsFactors = FALSE
  ) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(desc(abs(t))) %>%
    mutate(rank = row_number(),
           top_topn = rank <= TOPN_PER_CELLTYPE)
}

# ==============================================================================
# step4: main loop (cell_type x tissue x fold)
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
    logmsg(sprintf("[skip] %s: no matching l2 label (wanted=%s)",
                   ct_name, paste(CELL_TYPES[[ct_name]], collapse=",")))
    next
  }
  logmsg(sprintf("=== cell_type=%s (l2=%s) ===", ct_name, paste(ct_labels, collapse="/")))
  ct_cells <- rownames(so@meta.data)[so@meta.data[[CELLTYPE_COL]] %in% ct_labels]
  so_ct <- subset(so, cells = ct_cells)

  for (tissue in TISSUES) {
    use_cells <- colnames(so_ct)[so_ct@meta.data[[COMPARTMENT_COL]] == tissue &
                                   !is.na(so_ct@meta.data[[COMPARTMENT_COL]])]
    if (length(use_cells) < 100) {
      logmsg(sprintf("[skip] %s/%s: cells=%d (<100)", ct_name, tissue, length(use_cells)))
      n_skip <- n_skip + 1; next
    }
    # SKIP_EXISTING: skip the whole (cell_type, tissue) if all its folds already exist
    # note: the presence of count-scale counts (train_pb_counts.mtx) is part of the condition, so
    # legacy v5 directories without counts will have counts generated on re-run.
    # A fold counts as "done" only if the key outputs are all present, so a run
    # interrupted mid-fold is not skipped on restart. (tstats.csv is written last,
    # but we also require the count matrices and meta files to be safe.)
    fold_is_done <- function(fd) {
      all(file.exists(file.path(fd, c("tstats.csv",
                                      "train_pb_counts.mtx", "val_pb_counts.mtx",
                                      "train_meta.csv", "val_meta.csv"))))
    }
    if (SKIP_EXISTING) {
      done_dirs <- vapply(cv_fold_ids, function(k) {
        fold_is_done(file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k)))
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

    # the hold-out part is common to all folds per tissue, so compute it once
    heldout_cells <- colnames(so_t)[so_t@meta.data[[DONOR_COL]] %in% heldout_donors]
    so_held <- if (length(heldout_cells) >= 10) subset(so_t, cells = heldout_cells) else NULL
    pb_held_mean <- if (!is.null(so_held)) pseudobulk_mean(so_held, DONOR_COL) else NULL
    pb_held_sum  <- if (!is.null(so_held)) pseudobulk_sum(so_held,  DONOR_COL) else NULL
    pb_held_counts <- if (!is.null(so_held)) pseudobulk_count_sum(so_held, DONOR_COL) else NULL
    meta_held <- if (!is.null(so_held)) build_meta(so_held, "heldout") else NULL

    for (k in cv_fold_ids) {
      # also skip existing at the fold level
      fold_out_dir <- file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k))
      if (SKIP_EXISTING && fold_is_done(fold_out_dir)) {
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
        logmsg(sprintf("[skip] %s/%s/fold_%s: only %s in val",
                       ct_name, tissue, k, paste(val_dx, collapse="/")))
        n_skip <- n_skip + 1; next
      }

      so_train <- subset(so_t, cells = train_cells)
      so_val   <- subset(so_t, cells = val_cells)

      # ---- pure-technical gene removal BEFORE HVG (train cell only) ----
      # Remove non-informative high-abundance genes (mito / ribosomal / HSP /
      # housekeeping / nuclear lncRNA / small RNA) from the feature universe so
      # they cannot occupy HVG slots. V(D)J segments and the cell-type-aware
      # detection/specificity filter are applied downstream in 03_selection, not
      # here, since they define the candidate biomarker space rather than the
      # pseudobulk/HVG representation. ER chaperones HSPA5/HSP90B1 are retained.
      all_feats  <- rownames(so_train)
      tech_mask  <- grepl(PURE_TECH_PATTERN, all_feats, perl = TRUE) &
                    !(all_feats %in% BIOLOGY_RETAIN)
      keep_feats <- all_feats[!tech_mask]
      logmsg(sprintf("[filter] %s/%s/fold_%s: pure-technical removed %d/%d genes -> %d retained",
                     ct_name, tissue, k, sum(tech_mask), length(all_feats), length(keep_feats)))
      so_train <- subset(so_train, features = keep_feats)

      # ---- HVG (train cell only, computed on the pure-technical-filtered universe) ----
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
      pb_train_counts <- pseudobulk_count_sum(so_train, DONOR_COL)[hvg, , drop = FALSE]
      pb_val_counts   <- pseudobulk_count_sum(so_val,   DONOR_COL)
      pb_val_counts   <- pb_val_counts[intersect(hvg, rownames(pb_val_counts)), , drop = FALSE]

      # heldout already computed per tissue -> subset by HVG
      pb_held_mean_k <- if (!is.null(pb_held_mean)) pb_held_mean[intersect(hvg, rownames(pb_held_mean)), , drop = FALSE] else NULL
      pb_held_sum_k  <- if (!is.null(pb_held_sum))  pb_held_sum[intersect(hvg, rownames(pb_held_sum)),  , drop = FALSE]  else NULL
      pb_held_counts_k <- if (!is.null(pb_held_counts)) pb_held_counts[intersect(hvg, rownames(pb_held_counts)), , drop = FALSE] else NULL

      # ---- meta ----
      meta_train <- build_meta(so_train, "train")
      meta_val   <- build_meta(so_val,   "val")

      # align column order
      pb_train_mean <- pb_train_mean[, meta_train$donor_id, drop = FALSE]
      pb_val_mean   <- pb_val_mean[,   meta_val$donor_id,   drop = FALSE]
      pb_train_sum  <- pb_train_sum[,  meta_train$donor_id, drop = FALSE]
      pb_val_sum    <- pb_val_sum[,    meta_val$donor_id,   drop = FALSE]
      pb_train_counts <- pb_train_counts[, meta_train$donor_id, drop = FALSE]
      pb_val_counts   <- pb_val_counts[,   meta_val$donor_id,   drop = FALSE]
      if (POOLED) {
        # pooled CV: this fold's held-out test set IS the CV val fold
        pb_held_mean_k   <- pb_val_mean
        pb_held_sum_k    <- pb_val_sum
        pb_held_counts_k <- pb_val_counts
        meta_held        <- meta_val %>% mutate(set = "heldout")
      }
      if (!is.null(pb_held_mean_k)) {
        meta_held_aligned <- meta_held %>% filter(donor_id %in% colnames(pb_held_mean_k))
        pb_held_mean_k <- pb_held_mean_k[, meta_held_aligned$donor_id, drop = FALSE]
        pb_held_sum_k  <- pb_held_sum_k[,  meta_held_aligned$donor_id, drop = FALSE]
        pb_held_counts_k <- pb_held_counts_k[, meta_held_aligned$donor_id, drop = FALSE]
      } else {
        meta_held_aligned <- NULL
      }

      # ---- t-stat (train donor-level pseudobulk, mean primary) ----
      tstats_df <- compute_donor_tstats(pb_train_mean, meta_train,
                                        USE_AGE, USE_SEX, USE_BATCH)
      if (is.null(tstats_df) || nrow(tstats_df) == 0) {
        logmsg(sprintf("[skip] %s/%s/fold_%s: tstats failed", ct_name, tissue, k))
        n_skip <- n_skip + 1; next
      }
      topN_df <- tstats_df %>% filter(top_topn) %>% select(gene, t, pval, padj, log2FC, rank)

      # ---- write out ----
      out_dir <- file.path(OUT_ROOT, ct_name, tissue, sprintf("fold_%s", k))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      write_pb(pb_train_mean, "train_pb_mean", out_dir)
      write_pb(pb_val_mean,   "val_pb_mean",   out_dir)
      write_pb(pb_train_sum,  "train_pb_sum",  out_dir)
      write_pb(pb_val_sum,    "val_pb_sum",    out_dir)
      write_pb(pb_train_counts, "train_pb_counts", out_dir)
      write_pb(pb_val_counts,   "val_pb_counts",   out_dir)
      if (!is.null(pb_held_mean_k) && ncol(pb_held_mean_k) > 0) {
        write_pb(pb_held_mean_k, "heldout_pb_mean", out_dir)
        write_pb(pb_held_sum_k,  "heldout_pb_sum",  out_dir)
        write_pb(pb_held_counts_k, "heldout_pb_counts", out_dir)
        write.csv(meta_held_aligned, file.path(out_dir, "heldout_meta.csv"), row.names = FALSE)
      }
      write.csv(meta_train, file.path(out_dir, "train_meta.csv"), row.names = FALSE)
      write.csv(meta_val,   file.path(out_dir, "val_meta.csv"),   row.names = FALSE)
      write.csv(data.frame(gene = hvg), file.path(out_dir, "HVG.csv"), row.names = FALSE)
      # ---- Optional legacy / QC ranking — NOT used in the manuscript ----------
      # tstats.csv / topN_genes.csv are computed by compute_donor_tstats() from the
      # log-normalised pseudobulk mean via a simple donor-level linear model. They
      # are retained only for backward compatibility / QC and are NOT the canonical
      # feature ranking. The canonical differential-expression ranking (edgeR test statistic
      # |z| on *_pb_counts.mtx, weighted by cohort-consistency C_i) is computed
      # downstream in 02_deg and 03_selection.
      write.csv(tstats_df, file.path(out_dir, "tstats.csv"), row.names = FALSE)
      write.csv(topN_df,   file.path(out_dir, "topN_genes.csv"), row.names = FALSE)

      n_done <- n_done + 1
      logmsg(sprintf("[ok]   %s/%s/fold_%s: train=%dd val=%dd held=%dd  HVG=%d  topN=%d",
                     ct_name, tissue, k,
                     ncol(pb_train_mean), ncol(pb_val_mean),
                     ifelse(is.null(pb_held_mean_k), 0, ncol(pb_held_mean_k)),
                     length(hvg), nrow(topN_df)))

      # explicitly delete intermediate objects and GC (saves memory, incl. count objects)
      rm(so_train, so_val,
         pb_train_mean, pb_val_mean,
         pb_train_sum, pb_val_sum,
         pb_train_counts, pb_val_counts,
         pb_held_mean_k, pb_held_sum_k, pb_held_counts_k,
         tstats_df, topN_df, hvg, meta_train, meta_val)
      gc(verbose = FALSE)
    }
    rm(so_t, so_held, pb_held_mean, pb_held_sum, pb_held_counts, meta_held)
    gc(verbose = FALSE)
  }
  rm(so_ct)
  gc(verbose = FALSE)
}

logmsg(sprintf("DONE. ok=%d skip=%d  out=%s", n_done, n_skip, OUT_ROOT))
