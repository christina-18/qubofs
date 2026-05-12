# ==============================================================================
# run_aucell_analysis.R
# ------------------------------------------------------------------------------
# Real AUCell (Aibar et al. 2017, Nat Methods) per-cell pathway scoring.
#
# Computes AUCell AUC per cell for TWO families of gene sets in one pass:
#   (A) 7 literature-curated MS gene sets (MHC II, Iron, Cytotoxic, Type I IFN,
#       MS GWAS, B-cell DMT panel, MS DMT direct targets)
#   (B) QUBO-selected panels per (tissue × cell_type), loaded from
#       qubo_run_v6/aucell_results/qubo_panels.json (union mode)
#       — this answers: "Does the QUBO panel discriminate MS at the cell level?"
#
# Outputs (qubo_run_v6/aucell_results/):
#   cell_aucell_scores.csv.gz       cells × all_gene_sets matrix (A + B)
#   cell_metadata.csv               cell_id, donor, cell_type, tissue, Dx, cohort
#   summary_curated.csv             (ct × tissue × Dx × set) median/mean for (A)
#   summary_qubo.csv                same for (B)
#   ms_vs_hd_diff_curated.csv       MS-HD diff + Wilcoxon p + BH-FDR for (A)
#   ms_vs_hd_diff_qubo.csv          same for (B)
#
# Usage from project root:
#   1. python3 scripts/export_qubo_panels_for_aucell.py    # writes qubo_panels.json
#   2. Rscript scripts/run_aucell_analysis.R               # this script
#   3. python3 scripts/make_aucell_figure.py               # generates figures
#
# Estimated runtime: 15-40 min (cell-ranking step is dominant)
# Required R packages: Seurat, AUCell, GSEABase, dplyr, jsonlite
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(GSEABase)
  library(dplyr)
  library(jsonlite)
  library(msigdbr)   # bundled MSigDB v2024.1 — install: install.packages("msigdbr")
})

# ==============================================================================
# Paths
# ==============================================================================
PROJ <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO"
SEURAT_RDS <- "/Users/mizuhoasada/Dropbox/Research/_scRNAseq/MS/data/raw_data/so.GEX.share.Asada_with_compartment.rds"
OUT <- file.path(PROJ, "qubo_run_v6", "aucell_results")
QUBO_JSON <- file.path(OUT, "qubo_panels.json")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Project root : %s\n", PROJ))
cat(sprintf("Seurat .rds  : %s\n", SEURAT_RDS))
cat(sprintf("QUBO panels  : %s\n", QUBO_JSON))
cat(sprintf("Output dir   : %s\n", OUT))

if (!file.exists(QUBO_JSON)) {
  stop(sprintf("\nMissing %s.\nRun first:\n  python3 scripts/export_qubo_panels_for_aucell.py\n",
               QUBO_JSON))
}

# ==============================================================================
# (A) MS-relevant gene sets — HYBRID:
#     • 4 official MSigDB v2024.1 sets (Hallmark, Reactome, KEGG, GOBP)
#       pulled via msigdbr() — fully reproducible by ID
#     • 3 MS-specific sets (Iron rim, MS GWAS, MS DMT targets) curated from
#       primary literature (cited in Methods); MSigDB does not cover these.
# ==============================================================================

# --- (A1) Pull official MSigDB sets ---
cat("\nLoading MSigDB v2024.1 official gene sets via msigdbr...\n")
mdb_h  <- msigdbr(species = "Homo sapiens", category = "H")          # Hallmark
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")         # curated (Reactome/KEGG/etc)
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

pull_set <- function(df, name) {
  g <- unique(df$gene_symbol[df$gs_name == name])
  if (length(g) == 0) stop("MSigDB set not found: ", name)
  cat(sprintf("  %s: %d genes\n", name, length(g)))
  g
}

OFFICIAL_SETS <- list(
  HALLMARK_IFN_alpha_response   = pull_set(mdb_h,  "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  REACTOME_MHC_class_II_present = pull_set(mdb_c2, "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"),
  GOBP_NK_mediated_cytotoxicity = pull_set(mdb_c5, "GOBP_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"),
  KEGG_B_cell_receptor_signaling = pull_set(mdb_c2, "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY")
)

# --- (A2) MS-specific curated sets (NOT in MSigDB) ---
MS_CURATED_SETS <- list(
  Iron_metabolism_MS_lesion_Hametner2013 = c(
    "FTL","FTH1","TFRC","SLC11A2","SLC40A1",
    "ACO1","ACO2","HAMP","HFE","CP","TF","FTMT"
  ),
  MS_GWAS_topHits_IMSGC2019 = c(
    "HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DRB5","HLA-DRA","HLA-DPA1",
    "IL7R","IL2RA","CXCR4","TNFRSF1A","STAT4","EVI5",
    "CYP27B1","MERTK","RGS1","BACH2","CD58","CLEC16A",
    "TNFSF14","ZMIZ1","TYK2","CD86","MMEL1","EOMES"
  ),
  MS_DMT_direct_targets = c(
    "MS4A1","ITGA4","CD52","TYK2",
    "S1PR1","S1PR5","DHFR","GLB1"
  )
)

CURATED_SETS <- c(OFFICIAL_SETS, MS_CURATED_SETS)
cat(sprintf("\nTotal gene sets: %d (%d MSigDB official + %d MS-specific curated)\n",
            length(CURATED_SETS), length(OFFICIAL_SETS), length(MS_CURATED_SETS)))

# ==============================================================================
# (B) QUBO panels per (tissue × cell type) from JSON
# Default 'top30' = top 30 genes by selection frequency across the 15 CV runs.
# This is the standard convention for biomarker panel validation via AUCell
# (Aibar 2017 sweet spot 30-200 genes).
# ==============================================================================
qubo_json <- fromJSON(QUBO_JSON, simplifyVector = FALSE)
qubo_mode <- "top30"   # 'top30' (default, 30 genes/ct), 'union' (~30-160), or 'stable' (~3-15)
qubo_lists <- qubo_json[[qubo_mode]]
cat(sprintf("\nQUBO panel mode: %s\n", qubo_mode))

QUBO_SETS <- list()
for (tissue in names(qubo_lists)) {
  for (ct in names(qubo_lists[[tissue]])) {
    nm <- paste0("QUBO_", tissue, "_", ct)
    QUBO_SETS[[nm]] <- unlist(qubo_lists[[tissue]][[ct]])
  }
}

# Combine
ALL_SETS <- c(CURATED_SETS, QUBO_SETS)
cat(sprintf("\nGene sets:\n"))
cat(sprintf("  (A) curated  : %d sets\n", length(CURATED_SETS)))
cat(sprintf("  (B) QUBO     : %d panels (%d CSF + %d PBMC)\n",
            length(QUBO_SETS),
            sum(grepl("^QUBO_CSF",  names(QUBO_SETS))),
            sum(grepl("^QUBO_PBMC", names(QUBO_SETS)))))
cat(sprintf("  TOTAL        : %d gene sets to score\n", length(ALL_SETS)))
for (n in names(ALL_SETS)) cat(sprintf("    %-40s %d genes\n", n, length(ALL_SETS[[n]])))

gs_collection <- GeneSetCollection(
  lapply(names(ALL_SETS), function(n) {
    GeneSet(setName = n, geneIds = ALL_SETS[[n]],
            geneIdType = SymbolIdentifier())
  })
)

# ==============================================================================
# Load integrated Seurat object
# ==============================================================================
cat("\n[1/5] Loading Seurat object (1-3 min for 12 GB)...\n")
t0 <- Sys.time()
so <- readRDS(SEURAT_RDS)
cat(sprintf("  Loaded in %.1f sec, class=%s, %d cells × %d features\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            class(so)[1], ncol(so), nrow(so)))

md <- so@meta.data
cat("Available metadata columns (first 40):\n")
print(head(colnames(md), 40))

# ==============================================================================
# Cell-type / tissue / Dx column names — matched to Asada Seurat object
# ==============================================================================
CT_L2_COL  <- "predicted.celltype.l2"   # Azimuth fine annotation (~30 subtypes)
TISSUE_COL <- "compartment"             # CSF / PBMC
DX_COL     <- "Dx"                      # MS / HD
DONOR_COL  <- "patient_id"              # patient (subject) — matches pseudobulk pipeline
SID_COL    <- "sid"                     # sample (donor × tissue), kept for traceability
COHORT_COL <- "prj"                     # project / cohort identifier

required <- c(CT_L2_COL, TISSUE_COL, DX_COL, DONOR_COL, SID_COL, COHORT_COL)
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  cat("\n!!! Missing metadata columns:", paste(missing, collapse=", "), "\n")
  cat("Edit the *_COL constants above to match your Seurat object.\n")
  cat("Available columns:\n"); print(colnames(md))
  stop("Column name mismatch. Edit the script and re-run.")
}

# ------------------------------------------------------------------
# Collapse 30+ Azimuth l2 subtypes → 8 broad immune cell types
# (mirrors scripts/extract_pseudobulk_v5_compartment.R lines 88-100)
# ------------------------------------------------------------------
CELL_TYPE_MAP <- list(
  B      = c("B naive", "B memory", "B intermediate", "Plasmablast"),
  Mono   = c("CD14 Mono", "CD16 Mono"),
  CD4_T  = c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg"),
  CD8_T  = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "MAIT"),
  NK     = c("NK", "NK_CD56bright"),
  DC     = c("cDC1", "cDC2", "pDC"),
  dnT    = c("dnT"),
  gdT    = c("gdT")
)
# Build l2 → broad lookup
l2_to_broad <- unlist(lapply(names(CELL_TYPE_MAP), function(broad) {
  setNames(rep(broad, length(CELL_TYPE_MAP[[broad]])), CELL_TYPE_MAP[[broad]])
}))
md$cell_type_broad <- as.character(l2_to_broad[as.character(md[[CT_L2_COL]])])
n_mapped <- sum(!is.na(md$cell_type_broad))
n_total <- nrow(md)
cat(sprintf("\nCell-type collapse: %d / %d cells mapped to 8 broad types (%.1f%%)\n",
            n_mapped, n_total, 100 * n_mapped / n_total))
cat("Broad cell type counts:\n")
print(table(md$cell_type_broad, useNA = "ifany"))

# Use the broad cell type for all downstream aggregation
CT_COL <- "cell_type_broad"

# ============================================================================
# Combine broad-cell-type filter + subsampling into ONE subset call
# (Seurat::subset is expensive on 12GB objects; doing it once saves OOM.)
# ============================================================================
MAX_CELLS_PER_GROUP <- 1500   # max cells per (cell_type × tissue × Dx)
                              # → ~50k cells max, fits in 16 GB RAM
                              # increase if you have >32 GB RAM

cat(sprintf("\n[Filter + subsample in one pass] max %d cells per group\n",
            MAX_CELLS_PER_GROUP))
set.seed(42)
md$.cell_id <- rownames(md)
keep_ids <- md %>%
  dplyr::filter(!is.na(cell_type_broad)) %>%   # drop unmapped cell types
  dplyr::group_by(.data[[CT_COL]], .data[[TISSUE_COL]], .data[[DX_COL]]) %>%
  dplyr::slice_sample(n = MAX_CELLS_PER_GROUP, replace = FALSE) %>%
  dplyr::pull(.cell_id)
cat(sprintf("  Will keep %d cells (down from %d)\n", length(keep_ids), nrow(md)))

# Drop unused assays/layers BEFORE subsetting to reduce memory footprint
if ("DietSeurat" %in% getNamespaceExports("Seurat")) {
  cat("  Calling DietSeurat to drop unused layers...\n")
  so <- tryCatch(
    DietSeurat(so, assays = DefaultAssay(so),
               layers = c("counts"),
               dimreducs = NULL, graphs = NULL),
    error = function(e) {
      # Old Seurat (v4) uses different argument
      DietSeurat(so, assays = DefaultAssay(so),
                 dimreducs = NULL, graphs = NULL)
    }
  )
}

cat("  Subsetting Seurat object...\n")
so <- subset(so, cells = keep_ids)
gc(verbose = FALSE)   # immediate garbage collection
md <- so@meta.data
md$cell_type_broad <- as.character(l2_to_broad[as.character(md[[CT_L2_COL]])])
cat(sprintf("  After subset: %d cells × %d genes\n", ncol(so), nrow(so)))
cat("  Per-group counts:\n")
print(table(md$cell_type_broad, md[[TISSUE_COL]], md[[DX_COL]]))

# ==============================================================================
# AUCell — step 1: build cell rankings (slowest step)
# ==============================================================================
cat("\n[2/5] Building cell rankings...\n")
t0 <- Sys.time()
# Seurat 5.0+ uses 'layer' instead of 'slot'; fall back to old API for v4
expr_mat <- tryCatch(
  GetAssayData(so, assay = DefaultAssay(so), layer = "counts"),
  error = function(e) GetAssayData(so, assay = DefaultAssay(so), slot = "counts")
)
cat(sprintf("  Matrix: %d genes × %d cells (%.1f GB sparse)\n",
            nrow(expr_mat), ncol(expr_mat),
            as.numeric(object.size(expr_mat)) / 1e9))

cells_rankings <- AUCell_buildRankings(expr_mat,
                                       plotStats = FALSE, verbose = TRUE)
cat(sprintf("  Rankings built in %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

# ==============================================================================
# AUCell — step 2: AUC for each gene set per cell (cheap once rankings are built)
# ==============================================================================
cat("\n[3/5] Computing AUCell scores per (cell × gene set)...\n")
t0 <- Sys.time()
cells_AUC <- AUCell_calcAUC(gs_collection, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),
                            nCores = 1, verbose = TRUE)
auc_mat <- t(getAUC(cells_AUC))   # cells × gene_sets
cat(sprintf("  AUC computed in %.1f min, matrix: %d cells × %d sets\n",
            as.numeric(Sys.time() - t0, units = "mins"),
            nrow(auc_mat), ncol(auc_mat)))

# ==============================================================================
# Save per-cell AUC + metadata
# ==============================================================================
cat("\n[4/5] Saving cell-level outputs...\n")
auc_df <- as.data.frame(as.matrix(auc_mat))
auc_df$cell_id <- rownames(auc_df)

cell_md <- data.frame(
  cell_id    = rownames(md),
  patient_id = md[[DONOR_COL]],   # subject-level (matches pseudobulk pipeline)
  sid        = md[[SID_COL]],     # sample-level (donor × tissue)
  cell_type  = md[[CT_COL]],
  tissue     = md[[TISSUE_COL]],
  Dx         = md[[DX_COL]],
  cohort     = md[[COHORT_COL]],
  stringsAsFactors = FALSE
)
cell_md <- cell_md[cell_md$cell_id %in% auc_df$cell_id, ]

write.csv(auc_df, gzfile(file.path(OUT, "cell_aucell_scores.csv.gz")),
          row.names = FALSE)
write.csv(cell_md, file.path(OUT, "cell_metadata.csv"), row.names = FALSE)
cat(sprintf("  Wrote: cell_aucell_scores.csv.gz (%d × %d)\n",
            nrow(auc_df), ncol(auc_df)))
cat(sprintf("  Wrote: cell_metadata.csv (%d rows)\n", nrow(cell_md)))

# ==============================================================================
# Aggregate per (cell_type × tissue × Dx × gene_set)
# Two output families: curated (A) and QUBO (B)
# ==============================================================================
cat("\n[5/5] Aggregating to (cell_type × tissue × Dx × gene_set)...\n")
all_set_names <- setdiff(colnames(auc_df), "cell_id")
merged <- merge(cell_md, auc_df, by = "cell_id")

aggregate_family <- function(set_names, label) {
  summary_rows <- list(); diff_rows <- list()
  for (ct in unique(merged$cell_type)) {
    for (ts in unique(merged$tissue)) {
      sub <- merged[merged$cell_type == ct & merged$tissue == ts, , drop = FALSE]
      if (nrow(sub) < 20) next
      for (gs in set_names) {
        vals_ms <- sub[sub$Dx == "MS", gs]
        vals_hd <- sub[sub$Dx == "HD", gs]
        if (length(vals_ms) >= 10 && length(vals_hd) >= 10) {
          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            cell_type = ct, tissue = ts, gene_set = gs,
            n_MS = length(vals_ms), n_HD = length(vals_hd),
            median_MS = median(vals_ms), median_HD = median(vals_hd),
            mean_MS = mean(vals_ms), mean_HD = mean(vals_hd)
          )
          wt <- tryCatch(wilcox.test(vals_ms, vals_hd, exact = FALSE),
                         error = function(e) NULL)
          if (!is.null(wt)) {
            diff_rows[[length(diff_rows) + 1]] <- data.frame(
              cell_type = ct, tissue = ts, gene_set = gs,
              n_MS = length(vals_ms), n_HD = length(vals_hd),
              mean_MS = mean(vals_ms), mean_HD = mean(vals_hd),
              mean_diff_MS_minus_HD = mean(vals_ms) - mean(vals_hd),
              wilcox_p = wt$p.value
            )
          }
        }
      }
    }
  }
  list(summary = do.call(rbind, summary_rows),
       diff    = do.call(rbind, diff_rows))
}

curated_names <- intersect(names(CURATED_SETS), all_set_names)
qubo_names    <- intersect(names(QUBO_SETS),    all_set_names)

cat(sprintf("  Aggregating %d curated sets...\n", length(curated_names)))
A <- aggregate_family(curated_names, "curated")
A$diff$wilcox_q_BH <- p.adjust(A$diff$wilcox_p, method = "BH")
write.csv(A$summary, file.path(OUT, "summary_curated.csv"), row.names = FALSE)
write.csv(A$diff,    file.path(OUT, "ms_vs_hd_diff_curated.csv"), row.names = FALSE)
cat(sprintf("    Wrote summary_curated.csv (%d rows), ms_vs_hd_diff_curated.csv (%d rows)\n",
            nrow(A$summary), nrow(A$diff)))

cat(sprintf("  Aggregating %d QUBO panels...\n", length(qubo_names)))
B <- aggregate_family(qubo_names, "qubo")
B$diff$wilcox_q_BH <- p.adjust(B$diff$wilcox_p, method = "BH")
write.csv(B$summary, file.path(OUT, "summary_qubo.csv"), row.names = FALSE)
write.csv(B$diff,    file.path(OUT, "ms_vs_hd_diff_qubo.csv"), row.names = FALSE)
cat(sprintf("    Wrote summary_qubo.csv (%d rows), ms_vs_hd_diff_qubo.csv (%d rows)\n",
            nrow(B$summary), nrow(B$diff)))

cat("\n=== DONE ===\n")
cat(sprintf("All outputs in %s\n", OUT))
cat("Next: tell Claude 'AUCell finished' so it can build the Slide 4 figure.\n")
