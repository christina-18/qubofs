# ==============================================================================
# run_GO_enrichment.R
# ------------------------------------------------------------------------------
# Phase 2: GO / Reactome / KEGG 濃縮解析
# QUBO で選ばれた高頻度遺伝子の生物学的意義を評価する。
#
# 入力:
#   qubo_run_v6/gene_analysis/top_genes_QUBO_per_celltype.csv
#   qubo_run_v6/gene_analysis/top_genes_summary.csv
#
# 出力:
#   qubo_run_v6/gene_analysis/GO_results/
#     - GO_BP_per_celltype.csv        (各 cell type の GO Biological Process)
#     - GO_BP_overall.csv             (全 cell type union の GO BP)
#     - Reactome_overall.csv
#     - KEGG_overall.csv
#     - figures/GO_dotplot_<celltype>.pdf
#     - figures/Reactome_dotplot.pdf
#
# Usage:
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot"))
#   source("/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/scripts/run_GO_enrichment.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# --- パス ---
ROOT <- "/Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/qubo_run_v6/gene_analysis"
OUT  <- file.path(ROOT, "GO_results")
FIG  <- file.path(OUT, "figures")
dir.create(OUT, showWarnings = FALSE)
dir.create(FIG, showWarnings = FALSE)

# --- 入力読み込み ---
top_per_ct <- read.csv(file.path(ROOT, "top_genes_QUBO_per_celltype.csv"),
                       stringsAsFactors = FALSE)
top_overall <- read.csv(file.path(ROOT, "top_genes_summary.csv"),
                        stringsAsFactors = FALSE)
cat(sprintf("Loaded: %d per-celltype rows, %d overall rows\n",
            nrow(top_per_ct), nrow(top_overall)))

# Universe (background): ここでは HVG ~3000 を使うのが理想だが、
# 簡略化して "selected_genes_combined.csv" 全体を universe とする
all_selected <- read.csv(file.path(ROOT, "selected_genes_combined.csv"),
                         stringsAsFactors = FALSE)
universe_symbols <- unique(all_selected$gene)
cat(sprintf("Universe (background): %d genes\n", length(universe_symbols)))

# Symbol → Entrez 変換
sym_to_entrez <- function(symbols) {
  res <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  res$ENTREZID
}
universe_entrez <- sym_to_entrez(universe_symbols)

# ==============================================================================
# 1. Per cell type GO BP
# ==============================================================================
cell_types <- unique(top_per_ct$cell_type)
go_per_ct <- list()
for (ct in cell_types) {
  genes_ct <- top_per_ct$gene[top_per_ct$cell_type == ct]
  # 頻度 >= 0.2 のみに絞る (= 5 fold で 1 回以上)
  freq_ct <- top_per_ct$frequency[top_per_ct$cell_type == ct]
  genes_ct <- genes_ct[freq_ct >= 0.2]
  if (length(genes_ct) < 5) {
    cat(sprintf("  [skip] %s: only %d genes\n", ct, length(genes_ct)))
    next
  }
  entrez_ct <- sym_to_entrez(genes_ct)
  ego <- enrichGO(gene          = entrez_ct,
                  universe      = universe_entrez,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  if (!is.null(ego) && nrow(ego) > 0) {
    df <- as.data.frame(ego)
    df$cell_type <- ct
    go_per_ct[[ct]] <- df
    # plot top 10
    n_show <- min(10, nrow(ego))
    p <- dotplot(ego, showCategory = n_show, font.size = 9) +
      ggtitle(sprintf("GO BP enrichment — %s (n=%d genes)", ct, length(genes_ct))) +
      theme(plot.title = element_text(size = 11, face = "bold", color = "#1f3a5f"))
    ggsave(file.path(FIG, sprintf("GO_dotplot_%s.pdf", ct)),
           p, width = 8, height = 5.5)
    cat(sprintf("  [ok] %s: %d genes, %d GO terms\n", ct, length(genes_ct), nrow(ego)))
  } else {
    cat(sprintf("  [no enrichment] %s\n", ct))
  }
}

if (length(go_per_ct) > 0) {
  go_per_ct_df <- do.call(rbind, go_per_ct)
  write.csv(go_per_ct_df, file.path(OUT, "GO_BP_per_celltype.csv"), row.names = FALSE)
}

# ==============================================================================
# 2. Overall (union of top 50 across cell types) — GO + Reactome + KEGG
# ==============================================================================
cat("\n=== Overall (top 50 union) ===\n")
overall_genes <- unique(top_overall$gene[!is.na(top_overall$gene)])
overall_entrez <- sym_to_entrez(overall_genes)
cat(sprintf("Top overall genes: %d (entrez: %d)\n",
            length(overall_genes), length(overall_entrez)))

# GO BP
ego_all <- enrichGO(gene = overall_entrez, universe = universe_entrez,
                    OrgDb = org.Hs.eg.db, ont = "BP",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego_all) && nrow(ego_all) > 0) {
  write.csv(as.data.frame(ego_all), file.path(OUT, "GO_BP_overall.csv"), row.names = FALSE)
  p <- dotplot(ego_all, showCategory = 15, font.size = 10) +
    ggtitle("GO Biological Process — top 50 QUBO genes (CSF, all cell types)") +
    theme(plot.title = element_text(size = 11, face = "bold", color = "#1f3a5f"))
  ggsave(file.path(FIG, "GO_BP_overall_dotplot.pdf"), p, width = 9, height = 6.5)
  cat(sprintf("  GO BP: %d significant terms\n", nrow(ego_all)))
}

# Reactome
ereact <- enrichPathway(gene = overall_entrez, universe = universe_entrez,
                        organism = "human", pvalueCutoff = 0.05, readable = TRUE)
if (!is.null(ereact) && nrow(ereact) > 0) {
  write.csv(as.data.frame(ereact), file.path(OUT, "Reactome_overall.csv"), row.names = FALSE)
  p <- dotplot(ereact, showCategory = 15, font.size = 10) +
    ggtitle("Reactome pathway — top 50 QUBO genes") +
    theme(plot.title = element_text(size = 11, face = "bold", color = "#1f3a5f"))
  ggsave(file.path(FIG, "Reactome_overall_dotplot.pdf"), p, width = 9, height = 6.5)
  cat(sprintf("  Reactome: %d significant pathways\n", nrow(ereact)))
}

# KEGG
ekegg <- enrichKEGG(gene = overall_entrez, organism = "hsa",
                    pvalueCutoff = 0.05)
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg_readable), file.path(OUT, "KEGG_overall.csv"), row.names = FALSE)
  p <- dotplot(ekegg_readable, showCategory = 15, font.size = 10) +
    ggtitle("KEGG pathway — top 50 QUBO genes") +
    theme(plot.title = element_text(size = 11, face = "bold", color = "#1f3a5f"))
  ggsave(file.path(FIG, "KEGG_overall_dotplot.pdf"), p, width = 9, height = 6.5)
  cat(sprintf("  KEGG: %d significant pathways\n", nrow(ekegg)))
}

cat("\n=== DONE ===\n")
cat(sprintf("Output: %s\n", OUT))
