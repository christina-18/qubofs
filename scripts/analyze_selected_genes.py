"""
analyze_selected_genes.py
==========================
Aggregate v6full ensemble selected genes across holdouts × folds × methods,
generate publication-ready figures and tables for the gene-selection results.

Inputs:
  qubo_run_v6/v6full_edger*/  (3 holdouts × 2 tissues × selected_genes_*.csv)

Outputs (qubo_run_v6/gene_analysis/):
  - selected_genes_combined.csv       — all selections, long format
  - selection_frequency.csv           — gene × cell_type × method × tissue
  - top_genes_QUBO_per_celltype.csv   — top 20 per cell type for QUBO
  - top_genes_summary.csv             — top 50 across cell types (QUBO)
  - figures/
      - selection_frequency_heatmap_<tissue>.png
      - method_overlap_venn_<tissue>.png
      - top_genes_barplot.png
      - per_celltype_K_distribution.png
"""
from __future__ import annotations
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PROJECT_ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO")
RUN_ROOT = PROJECT_ROOT / "qubo_run_v6"
OUT = RUN_ROOT / "gene_analysis"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

HOLDOUTS = {
    "Pappalardo": "v6full_edger",
    "Heming":     "v6full_edger_holdout_Heming",
    "Ramesh":     "v6full_edger_holdout_Ramesh",
}
METHODS = ["QUBO", "DE_top", "HVG", "LASSO", "ElasticNet"]
TISSUES = ["CSF", "ALL"]
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]

# Housekeeping / metabolism / pseudogene patterns to flag (not "biology-relevant" for MS)
import re
HK_PATTERN = re.compile(
    r"^(MT-|MTRNR|MTATP|MTND|"           # mitochondrial
    r"RPL[0-9]|RPS[0-9]|MRPL|MRPS|"      # ribosomal
    r"HSP[A0-9]|HSPB|HSPA|HSPD|"         # heat shock
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|"  # housekeeping
    r"MALAT1$|NEAT1$|XIST$|TSIX$|"       # housekeeping lncRNA
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|" # uncharacterized / pseudogene
    r"MIR[0-9]|RNU[0-9]|SNORA|SNORD)"   # small RNAs
)
def is_biology_gene(g):
    return not bool(HK_PATTERN.match(g))


# ============================================================
# 1. Combine all selected_genes CSVs
# ============================================================
all_rows = []
for holdout, tag in HOLDOUTS.items():
    for tissue in TISSUES:
        files = sorted((RUN_ROOT / tag / tissue).glob("selected_genes_folds_*.csv"))
        for f in files:
            df = pd.read_csv(f)
            df["holdout"] = holdout
            all_rows.append(df)
combined = pd.concat(all_rows, ignore_index=True)
combined.to_csv(OUT / "selected_genes_combined.csv", index=False)
print(f"Combined: {len(combined)} rows")
print(f"  unique genes: {combined.gene.nunique()}")
print(f"  cell types: {combined.cell_type.unique().tolist()}")
print(f"  methods: {combined.method.unique().tolist()}")


# ============================================================
# 2. Selection frequency
# ============================================================
# For each (tissue, method, cell_type, gene): count how many fold-runs selected it
# Total runs per (tissue, method, cell_type) = 3 holdouts × 5 folds = 15 (max)
def compute_frequency(df, tissue, method):
    sub = df[(df.tissue == tissue) & (df.method == method)]
    # n_runs per cell type
    runs = (sub.groupby("cell_type").apply(
        lambda g: g.drop_duplicates(["holdout", "fold"]).shape[0])
        ).rename("n_runs").reset_index()
    # gene counts
    gc = (sub.groupby(["cell_type", "gene"]).apply(
        lambda g: g.drop_duplicates(["holdout", "fold"]).shape[0])
        ).rename("n_selected").reset_index()
    out = gc.merge(runs, on="cell_type")
    out["frequency"] = out["n_selected"] / out["n_runs"]
    out["tissue"] = tissue
    out["method"] = method
    return out


freq_rows = []
for tissue in TISSUES:
    for method in METHODS:
        f = compute_frequency(combined, tissue, method)
        if len(f) > 0:
            freq_rows.append(f)
freq_df = pd.concat(freq_rows, ignore_index=True)
freq_df["is_biology"] = freq_df["gene"].apply(is_biology_gene)
freq_df.to_csv(OUT / "selection_frequency.csv", index=False)
print(f"\nSelection frequency: {len(freq_df)} (gene × cell_type × method × tissue) entries")
print(f"  biology genes:    {freq_df.is_biology.sum()}")
print(f"  housekeeping etc: {(~freq_df.is_biology).sum()}")

# Biology-filtered subset (used for top-N tables and figures)
freq_bio = freq_df[freq_df.is_biology].copy()


# ============================================================
# 3. Top genes per cell type (QUBO, CSF) with t-stat info
# ============================================================
# Load tstats from one representative fold (fold_1, Pappalardo CSF) to get directionality
def get_tstats(tissue, holdout="Pappalardo"):
    """Load tstats_edger.csv from each (cell_type, fold_1) and combine."""
    tag = HOLDOUTS[holdout]
    data_root_map = {
        "Pappalardo": PROJECT_ROOT / "data" / "pseudobulk_v5_compartment",
        "Heming":     PROJECT_ROOT / "data" / "pseudobulk_v5_compartment_holdout_osmzhlab_MS_ence_cov",
        "Ramesh":     PROJECT_ROOT / "data" / "pseudobulk_v5_compartment_holdout_PRJNA549712_MS_PBMC_UCSF",
    }
    rows = []
    for ct in CELL_TYPES:
        for fold in range(1, 6):
            f = data_root_map[holdout] / ct / tissue / f"fold_{fold}" / "tstats_edger.csv"
            if f.exists():
                df = pd.read_csv(f)
                df["cell_type"] = ct
                df["fold"] = fold
                df["holdout"] = holdout
                rows.append(df[["gene", "t", "padj", "log2FC", "cell_type", "fold", "holdout"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


print("\nLoading t-stats for top-gene analysis...")
all_tstats = pd.concat([get_tstats(t) for t in TISSUES], ignore_index=True)
print(f"  loaded {len(all_tstats)} tstats rows")

# Top QUBO genes per cell type (CSF) — biology-filtered
csf_qubo = freq_bio[(freq_bio.tissue == "CSF") & (freq_bio.method == "QUBO")]
top_per_ct = []
for ct in CELL_TYPES:
    sub = csf_qubo[csf_qubo.cell_type == ct].sort_values("frequency", ascending=False).head(20)
    if len(sub) == 0: continue
    # mean t-stat across folds
    ts_ct = all_tstats[all_tstats.cell_type == ct].groupby("gene")[["t", "log2FC", "padj"]].mean()
    sub = sub.merge(ts_ct, left_on="gene", right_index=True, how="left")
    sub["direction"] = sub["t"].apply(lambda x: "MS↑" if (pd.notna(x) and x > 0) else ("HD↑" if pd.notna(x) else "?"))
    top_per_ct.append(sub)
top_per_ct_df = pd.concat(top_per_ct, ignore_index=True)
top_per_ct_df.to_csv(OUT / "top_genes_QUBO_per_celltype.csv", index=False)
print(f"  top QUBO genes per cell type: {len(top_per_ct_df)}")

# Top genes overall (frequency × cell_types involved)
overall = (csf_qubo.groupby("gene")
           .agg(total_freq=("frequency", "sum"),
                n_celltypes=("cell_type", "nunique"),
                celltypes=("cell_type", lambda s: ",".join(sorted(s.unique())))))
overall = overall.sort_values("total_freq", ascending=False).head(50).reset_index()
ts_avg = all_tstats.groupby("gene")[["t", "log2FC", "padj"]].mean()
overall = overall.merge(ts_avg, left_on="gene", right_index=True, how="left")
overall["direction"] = overall["t"].apply(lambda x: "MS↑" if (pd.notna(x) and x > 0) else ("HD↑" if pd.notna(x) else "?"))
overall.to_csv(OUT / "top_genes_summary.csv", index=False)
print(f"  top 50 overall: top 5 = {overall.head(5)['gene'].tolist()}")


# ============================================================
# 4. Selection frequency heatmap (top genes × cell_type)
# ============================================================
def plot_heatmap(tissue, method, top_n=30):
    sub = freq_bio[(freq_bio.tissue == tissue) & (freq_bio.method == method)]
    if len(sub) == 0: return
    # Top genes by total selection frequency across cell types
    gene_totals = sub.groupby("gene")["frequency"].sum().sort_values(ascending=False)
    top_genes = gene_totals.head(top_n).index.tolist()
    # pivot: rows = gene, columns = cell_type, values = frequency
    M = (sub[sub.gene.isin(top_genes)]
         .pivot_table(index="gene", columns="cell_type", values="frequency", fill_value=0)
         .reindex(index=top_genes, columns=CELL_TYPES, fill_value=0))

    fig, ax = plt.subplots(figsize=(7, 9), dpi=130)
    im = ax.imshow(M.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels(M.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(M.index)))
    ax.set_yticklabels(M.index, fontsize=8)
    ax.set_title(f"Selection frequency — {method} / {tissue} / top {top_n}",
                 fontsize=11, color="#1f3a5f", fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Frequency (over 3 holdouts × 5 folds)", fontsize=9)
    # annotate values
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if v > 0.05:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v > 0.5 else "black", fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG / f"selection_frequency_heatmap_{method}_{tissue}.png",
                dpi=130, bbox_inches="tight")
    plt.close(fig)


for tissue in TISSUES:
    for method in ["QUBO", "DE_top", "LASSO"]:
        plot_heatmap(tissue, method, top_n=30)
print("\nHeatmaps saved.")


# ============================================================
# 5. Method overlap (Venn-like for QUBO/LASSO/DE_top, CSF)
# ============================================================
def plot_venn(tissue, methods=("QUBO", "LASSO", "DE_top"), min_freq=0.3):
    """3-way set overlap counts (Venn-like; we draw circles manually)."""
    sets = {}
    for m in methods:
        sub = freq_bio[(freq_bio.tissue == tissue) & (freq_bio.method == m) &
                      (freq_bio.frequency >= min_freq)]
        sets[m] = set(sub["gene"].unique())

    s1, s2, s3 = sets[methods[0]], sets[methods[1]], sets[methods[2]]
    only1 = s1 - s2 - s3
    only2 = s2 - s1 - s3
    only3 = s3 - s1 - s2
    a12 = (s1 & s2) - s3
    a13 = (s1 & s3) - s2
    a23 = (s2 & s3) - s1
    a123 = s1 & s2 & s3

    fig, ax = plt.subplots(figsize=(7, 6), dpi=130)
    # 3 circles approx
    from matplotlib.patches import Circle
    c1 = Circle((-0.5, 0.3), 1.0, alpha=0.4, color="#d62728", label=methods[0])
    c2 = Circle(( 0.5, 0.3), 1.0, alpha=0.4, color="#ff7f0e", label=methods[1])
    c3 = Circle(( 0.0,-0.5), 1.0, alpha=0.4, color="#2ca02c", label=methods[2])
    for c in (c1, c2, c3): ax.add_patch(c)
    # Labels (counts)
    ax.text(-1.2, 0.3, f"{len(only1)}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text( 1.2, 0.3, f"{len(only2)}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text( 0.0,-1.2, f"{len(only3)}", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text( 0.0, 0.6, f"{len(a12)}", ha="center", va="center", fontsize=12)
    ax.text(-0.7,-0.4, f"{len(a13)}", ha="center", va="center", fontsize=12)
    ax.text( 0.7,-0.4, f"{len(a23)}", ha="center", va="center", fontsize=12)
    ax.text( 0.0, 0.0, f"{len(a123)}", ha="center", va="center", fontsize=14, fontweight="bold", color="#1f3a5f")
    # Method labels
    ax.text(-1.5, 1.4, methods[0], color="#d62728", fontweight="bold", fontsize=12)
    ax.text( 1.5, 1.4, methods[1], color="#ff7f0e", fontweight="bold", fontsize=12)
    ax.text( 0.0,-1.7, methods[2], color="#2ca02c", fontweight="bold", fontsize=12, ha="center")

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.0, 2.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Gene-set overlap — {tissue} / freq ≥ {min_freq}",
                 fontsize=12, color="#1f3a5f", fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG / f"method_overlap_venn_{tissue}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return {"only1": only1, "only2": only2, "only3": only3,
            "a12": a12, "a13": a13, "a23": a23, "a123": a123}


for tissue in TISSUES:
    overlap = plot_venn(tissue)
    print(f"\n{tissue} method overlap (freq ≥ 0.3):")
    print(f"  QUBO only: {len(overlap['only1'])}")
    print(f"  LASSO only: {len(overlap['only2'])}")
    print(f"  DE_top only: {len(overlap['only3'])}")
    print(f"  All 3 shared: {len(overlap['a123'])}")


# ============================================================
# 6. Top genes barplot (overall)
# ============================================================
def plot_top_barplot(tissue="CSF", n=20):
    sub = freq_bio[(freq_bio.tissue == tissue) & (freq_bio.method == "QUBO")]
    gene_totals = sub.groupby("gene")["frequency"].sum().sort_values(ascending=False).head(n)

    fig, ax = plt.subplots(figsize=(8, max(4, n * 0.3)), dpi=130)
    colors = ["#1f3a5f"] * n
    ax.barh(range(n), gene_totals.values[::-1], color=colors[::-1])
    ax.set_yticks(range(n))
    ax.set_yticklabels(gene_totals.index[::-1], fontsize=10)
    ax.set_xlabel("Total selection frequency (sum across 8 cell types)", fontsize=10)
    ax.set_title(f"Top {n} QUBO-selected genes — {tissue}",
                 fontsize=12, color="#1f3a5f", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(FIG / f"top_genes_barplot_{tissue}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


for tissue in TISSUES:
    plot_top_barplot(tissue, n=20)


# ============================================================
# 7. K distribution per cell type (QUBO)
# ============================================================
qubo_k = (combined[combined.method == "QUBO"][["holdout", "tissue", "fold", "cell_type", "K_chosen"]]
          .drop_duplicates()
          .copy())

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=130, sharey=True)
for ax, tissue in zip(axes, TISSUES):
    sub = qubo_k[qubo_k.tissue == tissue]
    counts = (sub.groupby(["cell_type", "K_chosen"]).size()
              .unstack(fill_value=0).reindex(CELL_TYPES, fill_value=0))
    bottom = np.zeros(len(counts))
    colors = {10: "#a8d5a8", 20: "#5b9bd5", 30: "#1f3a5f"}
    for k in [10, 20, 30]:
        if k not in counts.columns: continue
        ax.bar(range(len(counts)), counts[k].values, bottom=bottom,
               color=colors[k], label=f"K={k}", edgecolor="white", linewidth=0.5)
        bottom += counts[k].values
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("# folds (3 holdouts × 5 folds = 15 max)" if tissue == TISSUES[0] else "")
    ax.set_title(f"K distribution — {tissue}", fontweight="bold", color="#1f3a5f")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
axes[-1].legend(loc="upper right", frameon=False)
plt.suptitle("QUBO-selected K (gene panel size) per cell type",
             fontsize=13, color="#1f3a5f", fontweight="bold")
plt.tight_layout()
fig.savefig(FIG / "K_distribution_per_celltype.png", dpi=130, bbox_inches="tight")
plt.close(fig)
print("\nK distribution plot saved.")


print(f"\n=== DONE ===")
print(f"Outputs in: {OUT}")
