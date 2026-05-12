"""Generate Slide 4 right-panel figures from AUCell outputs.

Reads:
  qubo_run_v6/aucell_results/summary_curated.csv
  qubo_run_v6/aucell_results/ms_vs_hd_diff_curated.csv
  qubo_run_v6/aucell_results/summary_qubo.csv
  qubo_run_v6/aucell_results/ms_vs_hd_diff_qubo.csv

Writes:
  qubo_run_v6/figures_genes/fig4a_aucell_qubo_panels_csf.png
      ★ MAIN — Slide 4 right panel: per-cell-type QUBO panel applied to its
        matched cell type. Two heatmaps: median AUCell (MS cells) + MS-HD diff.
  qubo_run_v6/figures_genes/fig4b_aucell_curated_csf.png
      ★ SUPP S6 / S8 — biological-axis check across 7 curated MS gene sets:
        cell type × curated gene set, MS-HD diff with significance stars.
  PBMC versions (same layout, _pbmc.png).
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
RES = PROJ / "qubo_run_v6" / "aucell_results"
FIG_DIR = PROJ / "qubo_run_v6" / "figures_genes"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CT_ORDER = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]

CURATED_GS_ORDER = [
    "MHC_II_antigen_presentation",
    "Iron_metabolism_MSlesion",
    "Cytotoxic_NK_T_effectors",
    "Type_I_IFN_signature",
    "MS_GWAS_topHits",
    "B_cell_panel_DMT_relevant",
    "MS_DMT_direct_targets",
]
CURATED_GS_LABELS = {
    "MHC_II_antigen_presentation": "MHC II antigen presentation\n[Reactome]",
    "Iron_metabolism_MSlesion": "Iron metabolism (MS lesion)\n[Hametner 2013]",
    "Cytotoxic_NK_T_effectors": "Cytotoxic NK / T effectors\n[curated]",
    "Type_I_IFN_signature": "Type I IFN signature\n[van Langelaar 2020]",
    "MS_GWAS_topHits": "MS GWAS top hits\n[IMSGC 2019]",
    "B_cell_panel_DMT_relevant": "B cell / DMT panel\n[curated]",
    "MS_DMT_direct_targets": "MS DMT direct targets\n[curated]",
}


# ---------------------------------------------------------------------------
# (B) QUBO panels figure — diagonal: matched cell-type panel
#     Shows all 8 cell types; cell types without a QUBO panel are explicitly
#     marked "No QUBO panel (QC skipped, see S8)" rather than hidden.
# ---------------------------------------------------------------------------
def make_qubo_diagonal(tissue):
    s = pd.read_csv(RES / "summary_qubo.csv")
    d = pd.read_csv(RES / "ms_vs_hd_diff_qubo.csv")
    s = s[s["tissue"] == tissue]
    d = d[d["tissue"] == tissue]

    # gene_set name like "QUBO_CSF_B" → cell_type "B"
    s["panel_ct"] = s["gene_set"].str.extract(rf"QUBO_{tissue}_(.+)$")[0]
    d["panel_ct"] = d["gene_set"].str.extract(rf"QUBO_{tissue}_(.+)$")[0]

    # Diagonal: panel_ct == cell_type (apply panel to its own cell type)
    s_diag = s[s["cell_type"] == s["panel_ct"]].copy().set_index("cell_type")
    d_diag = d[d["cell_type"] == d["panel_ct"]].copy().set_index("cell_type")

    # Reindex over ALL 8 cell types — missing rows render as N/A
    s_diag = s_diag.reindex(CT_ORDER)
    d_diag = d_diag.reindex(CT_ORDER)
    has_panel = ~s_diag["median_MS"].isna()

    fig, axes = plt.subplots(1, 2, figsize=(13, max(4.2, 0.65 * len(CT_ORDER) + 1.5)),
                             gridspec_kw={"width_ratios": [1, 1]})

    y_labels = []
    for ct in CT_ORDER:
        if has_panel.loc[ct]:
            n_ms = int(s_diag.loc[ct, "n_MS"])
            n_hd = int(s_diag.loc[ct, "n_HD"])
            y_labels.append(f"{ct}\n(MS={n_ms}, HD={n_hd})")
        else:
            y_labels.append(f"{ct}\n(no QUBO panel)")

    # ---- Left: median AUCell of MS cells ----
    ax = axes[0]
    vals = s_diag["median_MS"].values
    colors = ["#2874a6" if hp else "#cccccc" for hp in has_panel]
    ax.barh(range(len(CT_ORDER)), [v if not np.isnan(v) else 0 for v in vals],
            color=colors)
    ax.set_yticks(range(len(CT_ORDER)))
    ax.set_yticklabels(y_labels, fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    for i, (v, hp) in enumerate(zip(vals, has_panel)):
        if hp and not np.isnan(v):
            ax.text(v + 0.005, i, f"{v:.3f}", va="center",
                    fontsize=11, fontweight="bold")
        else:
            ax.text(0.002, i, "  No QUBO panel — see S8 (CD4/CD8 dilution)",
                    va="center", fontsize=10, color="#666666", style="italic")
    ax.set_xlabel("AUCell AUC (median, MS cells)", fontsize=12, fontweight="bold")
    ax.set_title(f"QUBO panel activity — {tissue} MS cells",
                 fontsize=13, weight="bold", color="#1f3a5f", pad=10)
    ax.tick_params(axis="x", labelsize=10)
    ax.grid(axis="x", alpha=0.3)
    valid_max = np.nanmax(vals) if has_panel.any() else 0.05
    ax.set_xlim(0, max(0.05, valid_max * 1.3))

    # ---- Right: MS-HD diff with significance ----
    ax = axes[1]
    diffs = d_diag["mean_diff_MS_minus_HD"].values
    qs = d_diag["wilcox_q_BH"].values
    bar_colors = []
    for v, hp in zip(diffs, has_panel):
        if not hp or np.isnan(v):
            bar_colors.append("#cccccc")
        elif v > 0:
            bar_colors.append("#c0392b")
        else:
            bar_colors.append("#2874a6")
    ax.barh(range(len(CT_ORDER)),
            [v if not np.isnan(v) else 0 for v in diffs],
            color=bar_colors)
    ax.set_yticks(range(len(CT_ORDER))); ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.5)
    for i, (v, q, hp) in enumerate(zip(diffs, qs, has_panel)):
        if not hp or pd.isna(v):
            ax.text(0.0, i, "   N/A", va="center", ha="left",
                    fontsize=10, color="#666666", style="italic")
            continue
        star = ""
        if not pd.isna(q):
            if q < 0.001: star = "***"
            elif q < 0.01: star = "**"
            elif q < 0.05: star = "*"
        offset = 0.0005 if v >= 0 else -0.0005
        ha = "left" if v >= 0 else "right"
        ax.text(v + offset, i, f"{v:+.4f}{star}",
                va="center", ha=ha, fontsize=11, fontweight="bold")
    ax.set_xlabel("Δ AUCell (MS − HD)", fontsize=12, fontweight="bold")
    ax.set_title(f"QUBO panel: MS − HD effect ({tissue})\n* q<0.05  ** q<0.01  *** q<0.001",
                 fontsize=13, weight="bold", color="#1f3a5f", pad=10)
    ax.tick_params(axis="x", labelsize=10)
    ax.grid(axis="x", alpha=0.3)

    fig.suptitle(
        f"AUCell — QUBO panels per cell type ({tissue})  ·  Aibar et al. 2017, Nat Methods",
        fontsize=14, weight="bold", y=0.998, color="#1f3a5f")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG_DIR / f"fig4a_aucell_qubo_panels_{tissue.lower()}.png"
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# (A) Curated MS gene sets — cell_type × gene_set heatmap, MS-HD diff
# ---------------------------------------------------------------------------
def make_curated_heatmap(tissue):
    s = pd.read_csv(RES / "summary_curated.csv")
    d = pd.read_csv(RES / "ms_vs_hd_diff_curated.csv")
    s = s[s["tissue"] == tissue]
    d = d[d["tissue"] == tissue]
    if s.empty or d.empty:
        print(f"[curated / {tissue}] no data — skipping")
        return

    pivot_ms = s.pivot(index="gene_set", columns="cell_type", values="median_MS") \
                .reindex(index=CURATED_GS_ORDER, columns=CT_ORDER)
    pivot_diff = d.pivot(index="gene_set", columns="cell_type",
                         values="mean_diff_MS_minus_HD") \
                  .reindex(index=CURATED_GS_ORDER, columns=CT_ORDER)
    pivot_q = d.pivot(index="gene_set", columns="cell_type", values="wilcox_q_BH") \
                .reindex(index=CURATED_GS_ORDER, columns=CT_ORDER)

    # Cell counts per cell type (any Dx)
    n_per_ct = (s.groupby("cell_type")[["n_MS", "n_HD"]].first()
                  .reindex(CT_ORDER).fillna(0).astype(int))
    # Compact labels (cell type only, n shown separately on subtitle)
    ct_labels = [f"{ct}" for ct in CT_ORDER]
    gs_labels = [CURATED_GS_LABELS[g] for g in CURATED_GS_ORDER]
    n_total = sum(n_per_ct.loc[ct, "n_MS"] + n_per_ct.loc[ct, "n_HD"]
                  for ct in CT_ORDER)

    # Larger figure + bigger fonts so labels stay readable when shrunk on slides
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5),
                             gridspec_kw={"width_ratios": [1, 1]})

    # Left: median AUCell in MS cells
    ax = axes[0]
    im = ax.imshow(pivot_ms.values, cmap="Blues", aspect="auto",
                   vmin=np.nanpercentile(pivot_ms.values, 2),
                   vmax=np.nanpercentile(pivot_ms.values, 98))
    ax.set_xticks(range(len(CT_ORDER)))
    ax.set_xticklabels(ct_labels, fontsize=13, fontweight="bold",
                       rotation=30, ha="right")
    ax.set_yticks(range(len(CURATED_GS_ORDER)))
    ax.set_yticklabels(gs_labels, fontsize=11)
    ax.set_title(f"AUCell activity (median, {tissue} MS cells)",
                 fontsize=14, weight="bold", color="#1f3a5f", pad=10)
    for i in range(pivot_ms.shape[0]):
        for j in range(pivot_ms.shape[1]):
            v = pivot_ms.values[i, j]
            if not np.isnan(v):
                txtcol = "white" if v > np.nanmedian(pivot_ms.values) else "#303030"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=11, fontweight="bold", color=txtcol)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("AUCell AUC", fontsize=11, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    # Right: MS-HD effect
    ax = axes[1]
    vmax = np.nanmax(np.abs(pivot_diff.values))
    im = ax.imshow(pivot_diff.values, cmap="RdBu_r", aspect="auto",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CT_ORDER)))
    ax.set_xticklabels(ct_labels, fontsize=13, fontweight="bold",
                       rotation=30, ha="right")
    ax.set_yticks(range(len(CURATED_GS_ORDER)))
    ax.set_yticklabels(gs_labels, fontsize=11)
    ax.set_title(f"MS − HD effect ({tissue}; * q<0.05  ** q<0.01)",
                 fontsize=14, weight="bold", color="#1f3a5f", pad=10)
    for i in range(pivot_diff.shape[0]):
        for j in range(pivot_diff.shape[1]):
            v = pivot_diff.values[i, j]
            q = pivot_q.values[i, j]
            if not np.isnan(v):
                star = ""
                if not np.isnan(q):
                    if q < 0.01: star = "**"
                    elif q < 0.05: star = "*"
                ax.text(j, i, f"{v:+.3f}{star}", ha="center", va="center",
                        fontsize=10, fontweight="bold", color="#303030")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Δ AUCell (MS − HD)", fontsize=11, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    fig.suptitle(
        f"AUCell — cell-type × MS gene set ({tissue})  ·  Aibar et al. 2017, Nat Methods",
        fontsize=14, weight="bold", y=0.995, color="#1f3a5f")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / f"fig4b_aucell_curated_{tissue.lower()}.png"
    plt.savefig(out, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  wrote {out}")


def main():
    for tissue in ["CSF", "PBMC"]:
        print(f"\n=== {tissue} ===")
        if (RES / "summary_qubo.csv").exists():
            make_qubo_diagonal(tissue)
        if (RES / "summary_curated.csv").exists():
            make_curated_heatmap(tissue)


if __name__ == "__main__":
    main()
