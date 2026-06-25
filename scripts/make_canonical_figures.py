"""Regenerate the data-driven main and supplementary figures from canonical outputs.

Reads only canonical CSV outputs from `qubo_run/` (no step1c_*, no qubo_run_v6,
no synthetic data):
  qubo_run/primary_summary_per_holdout.csv       -> Figure 2 (cohort robustness)
  qubo_run/table1_<run_tag>.csv                  -> Figure 3 (perf vs redundancy)
  qubo_run/within_panel_redundancy_perpanel.csv  -> Figure 3 (per-panel |rho|)

Main figures generated:
  figures_oup/figure2_cohort_robustness.{png,pdf}
  figures_oup/figure3_perf_redundancy.{png,pdf}
  figures_oup/figure4_biology.{png,pdf}

Supplementary figures generated when the required input summaries are present:
  figures_oup/figureS1_dataset_composition.{png,pdf}
  figures_oup/figureS2_K_sweep.{png,pdf}
  figures_oup/figureS3_solver_sensitivity.{png,pdf}
  figures_oup/figureS4_literature_concordance.{png,pdf}

Figure 1 is a curated workflow graphic and is not overwritten by this script.

Style: clean OUP/Bioinformatics-compatible, colour-blind-friendly palette — the
Okabe–Ito qualitative palette for method comparisons (Figures 2, 3, S2, S3) and a
pale per-cell-type palette (hue = cell type, intensity = selection frequency) for
the Figure 4 selection-frequency heatmap. White background, grey axes, 300 dpi.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(os.environ.get("QUBOFS_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
RUN = PROJ / "qubo_run"
OUT = PROJ / "figures_oup"
OUT.mkdir(exist_ok=True)
TAG = os.environ.get("QUBOFS_RUN_TAG", "primary_bio_edger_counts")

# method order (QUBO first) + colourblind-friendly palette (Okabe-Ito)
METHODS = ["QUBO", "mRMR", "DE_top", "ElasticNet", "LASSO", "HVG"]
LABEL = {"QUBO": "quboFS", "DE_top": "DE-top", "ElasticNet": "Elastic Net",
         "mRMR": "mRMR", "LASSO": "LASSO", "HVG": "HVG"}
COL = {"QUBO": "#D55E00", "mRMR": "#0072B2", "DE_top": "#009E73",
       "ElasticNet": "#56B4E9", "LASSO": "#CC79A7", "HVG": "#999999"}
COHORTS = ["Pappalardo", "Heming", "Ramesh"]


def _read_tagged(stem):
    """Read an aggregation CSV, preferring the tag-suffixed file (per DE source)
    over the legacy fixed-name file. This makes the figures pull the correct DE
    source (QUBOFS_RUN_TAG) regardless of which run wrote the legacy file last."""
    tagged = RUN / f"{stem}_{TAG}.csv"
    return pd.read_csv(tagged if tagged.exists() else RUN / f"{stem}.csv")

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans",
                     "savefig.dpi": 300, "savefig.bbox": "tight",
                     "axes.spines.top": False, "axes.spines.right": False})


def figure2():
    """Cohort-level performance heatmaps: rows = methods, cols = held-out cohorts.
    (A) ROC-AUC, (B) MCC. Blues, annotated values."""
    ph = _read_tagged("primary_summary_per_holdout")
    ph = ph[ph["tissue"] == "CSF"]
    order = ["QUBO", "LASSO", "DE_top", "mRMR", "ElasticNet", "HVG"]
    rows = [LABEL[m] for m in order]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    panels = [("held_auc_mean", "A. ROC-AUC by held-out cohort", "ROC-AUC", 0.4, 1.0),
              ("held_mcc_mean", "B. MCC by held-out cohort", "MCC", -0.1, 0.6)]
    for ax, (col, title, cbar, vmin, vmax) in zip(axes, panels):
        piv = ph.pivot_table(index="method", columns="holdout", values=col).reindex(
            index=order, columns=COHORTS)
        im = ax.imshow(piv.values, aspect="auto", cmap="Blues", vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(COHORTS))); ax.set_xticklabels(COHORTS, rotation=20, ha="right")
        ax.set_yticks(range(len(order))); ax.set_yticklabels(rows)
        ax.set_title(title, fontsize=11, fontweight="bold")
        rng = vmax - vmin
        for i in range(len(order)):
            for j in range(len(COHORTS)):
                v = piv.values[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if (v - vmin) / rng > 0.6 else "#222")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbar)
    fig.tight_layout()
    fig.savefig(OUT / "figure2_cohort_robustness.png")
    fig.savefig(OUT / "figure2_cohort_robustness.pdf")
    plt.close(fig)
    print("wrote figure2_cohort_robustness")


MARK = {"QUBO": "o", "mRMR": "D", "DE_top": "s", "LASSO": "v",
        "ElasticNet": "^", "HVG": "X"}
# per-method label offset (points) to avoid overlap
LAB_OFF = {"QUBO": (4, 9), "mRMR": (7, -14), "DE_top": (-2, -16), "LASSO": (-10, -15),
           "ElasticNet": (6, 9), "HVG": (6, -14)}


def figure3():
    """Redundancy-focused figure (2 panels; no ROC-AUC, to avoid duplicating Table 2):
    A. Excess redundancy of each baseline relative to quboFS — paired
       Δ within-panel |rho| = baseline − quboFS (mean ± 95% CI over the 119
       cohort×fold×cell-type panels); positive = more redundant than quboFS;
    B. within-panel |rho| distribution across the 119 panels (ascending median)."""
    t = pd.read_csv(RUN / f"table1_{TAG}.csv").set_index("method")
    pp = _read_tagged("within_panel_redundancy_perpanel")

    # ---- A: paired per-panel difference (baseline − quboFS) ----
    wide = pp.pivot_table(index=["holdout", "fold", "cell_type"],
                          columns="method", values="rho")
    baselines = ["mRMR", "LASSO", "DE_top", "ElasticNet", "HVG"]
    stats = {}
    for m in baselines:
        d = (wide[m] - wide["QUBO"]).dropna().values
        n = len(d)
        half = 1.96 * d.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
        stats[m] = (float(d.mean()), float(half))
    order_a = sorted(baselines, key=lambda m: stats[m][0])  # ascending Δ

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    xs = np.arange(len(order_a))
    means = [stats[m][0] for m in order_a]
    halves = [stats[m][1] for m in order_a]
    ax.bar(xs, means, 0.62, yerr=halves, capsize=3,
           color=[COL[m] for m in order_a], edgecolor="black", linewidth=0.6,
           error_kw=dict(elinewidth=0.9, ecolor="#333"))
    ax.axhline(0, color="black", lw=0.9)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{LABEL[m]}\n− quboFS" for m in order_a], fontsize=8.5)
    ax.set_ylabel("Δ within-panel |rho| vs quboFS\n(higher = more redundant)")
    ax.set_title("A. Redundancy excess relative to quboFS", fontsize=11, fontweight="bold")
    for x, m in zip(xs, order_a):
        v = stats[m][0]
        ax.text(x, v + stats[m][1] + 0.004, f"+{v:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, max(mm + hh for mm, hh in zip(means, halves)) * 1.18)

    # ---- B: per-panel |rho| distribution (ascending median) ----
    order = sorted(METHODS, key=lambda m: t.loc[m, "within_panel_rho"])
    data = [pp[pp["method"] == m]["rho"].dropna().values for m in order]
    bp = axes[1].boxplot(data, vert=True, patch_artist=True, widths=0.6,
                         showfliers=False, medianprops=dict(color="black"))
    for patch, m in zip(bp["boxes"], order):
        patch.set_facecolor(COL[m]); patch.set_alpha(0.85)
    axes[1].set_xticks(range(1, len(order) + 1))
    axes[1].set_xticklabels([LABEL[m] for m in order], rotation=30, ha="right")
    axes[1].set_ylabel("within-panel |rho| per panel")
    axes[1].set_title("B. Within-panel redundancy distribution", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "figure3_perf_redundancy.png")
    fig.savefig(OUT / "figure3_perf_redundancy.pdf")
    plt.close(fig)
    print("wrote figure3_perf_redundancy (redundancy-focused: A=excess Δ|ρ|, B=distribution)")


def figureS1():
    """Dataset composition: donors per cohort by diagnosis (CSF). Fixed, known
    counts from the integrated object (50 donors, 221,066 CSF cells)."""
    comp = [("Pappalardo", 6, 5, "held-out"),
            ("Heming", 9, 9, "held-out"),
            ("Ramesh", 3, 14, "held-out"),
            ("Touil", 4, 0, "training only")]
    names = [c[0] for c in comp]
    ctrl = np.array([c[1] for c in comp]); ms = np.array([c[2] for c in comp])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(len(names))
    ax.bar(x, ctrl, 0.6, label="control", color="#56B4E9", edgecolor="black", linewidth=0.5)
    ax.bar(x, ms, 0.6, bottom=ctrl, label="MS", color="#D55E00", edgecolor="black", linewidth=0.5)
    for i, c in enumerate(comp):
        ax.text(i, c[1] + c[2] + 0.4, f"n={c[1]+c[2]}\n({c[3]})", ha="center",
                va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(names)
    ax.set_ylabel("Donors"); ax.set_ylim(0, 22)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("CSF dataset composition (50 donors, 221,066 cells)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "figureS1_dataset_composition.png")
    fig.savefig(OUT / "figureS1_dataset_composition.pdf")
    plt.close(fig)
    print("wrote figureS1_dataset_composition")


def figure1():
    """Workflow schematic with canonical pipeline wording (no data)."""
    fig, ax = plt.subplots(figsize=(8.2, 9.6)); ax.axis("off")
    ax.set_xlim(0, 10); ax.set_ylim(0, 20)
    steps = [
        ("Input — integrated CSF scRNA-seq",
         "4 cohorts (Pappalardo, Heming, Ramesh + Touil control-only)\n50 donors · 221,066 cells", "#F4D9D0"),
        ("Cell-type annotation (Azimuth)",
         "8 immune cell types: B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT", "#DDEAF6"),
        ("Per-donor pseudobulk (CSF only · neutral)",
         "log-norm mean (classification) + count-sum pseudobulk (edgeR)\ndonor = statistical unit; no gene filter here", "#DDEAF6"),
        ("Relevance + pre-specified technical/clonotype filter",
         "relevance = |z|·C  (edgeR test statistic |z| × cohort-consistency C)\nfilter: drop MT/RPL/RPS/stress-HSP/lncRNA/HK + Ig·TCR V(D)J;\nretain Ig/TCR constant + ER chaperones", "#DDEAF6"),
        ("QUBO selection (per cell type)",
         "min  −sᵀx + γ·xᵀRx + λ(Σx−K)²   ·  K fixed=10 (primary); sweep 5,10,15\nclassical simulated annealing (no quantum hardware)", "#CBD9EC"),
        ("Benchmark: 6 feature selectors",
         "QUBO vs DE-top, mRMR, LASSO, Elastic Net, HVG\n(shared candidate pool, filter, fixed K=10, L2 logistic, ensemble)", "#DDEAF6"),
        ("Evaluation — leave-one-cohort-out",
         "L2 logistic per cell type → soft-vote → donor MS/HD\n3 held-out cohorts · fixed 0.5 threshold + threshold sensitivity", "#CBD9EC"),
    ]
    n = len(steps); top = 19.2; h = 1.9; gap = (top - 0.6) / n
    for i, (title, body, col) in enumerate(steps):
        y = top - i * gap
        box = plt.Rectangle((1.0, y - h), 8.0, h, facecolor=col,
                            edgecolor="#333", linewidth=1.0, zorder=2)
        ax.add_patch(box)
        ax.text(5.0, y - 0.45, title, ha="center", va="center",
                fontsize=10.5, fontweight="bold", zorder=3)
        ax.text(5.0, y - 1.32, body, ha="center", va="center", fontsize=8.2, zorder=3)
        if i < n - 1:
            ax.annotate("", xy=(5.0, y - h - gap + h + 0.02), xytext=(5.0, y - h - 0.02),
                        arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.4))
    fig.tight_layout()
    fig.savefig(OUT / "figure1_pipeline.png")
    fig.savefig(OUT / "figure1_pipeline.pdf")
    plt.close(fig)
    print("wrote figure1_pipeline")


# Figure 4 gene annotation — curated MS-relevant biology categories (NPG palette).
FIG4_NPG = {
    "B-cell activation": "#7E6148", "B-cell secretory / ER": "#E64B35",
    "Ig constant": "#B03A20",
    "Iron metabolism": "#F39B7F", "MHC II / GWAS": "#8491B4",
    "Cytotoxic": "#00A087", "Type I IFN": "#4DBBD5", "TF / regulator": "#3C5488",
}
FIG4_BIOLOGY = {
    "LY9": "B-cell activation", "CCR7": "B-cell activation",
    "IL2RG": "B-cell activation", "CD70": "B-cell activation",
    "LTB": "B-cell activation",
    "XBP1": "B-cell secretory / ER", "CREB3L2": "B-cell secretory / ER",
    "DNAJB11": "B-cell secretory / ER", "FKBP11": "B-cell secretory / ER",
    "HSP90B1": "B-cell secretory / ER", "KDELR1": "B-cell secretory / ER",
    "LMAN2": "B-cell secretory / ER", "P4HB": "B-cell secretory / ER",
    "TMEM167A": "B-cell secretory / ER", "MZB1": "B-cell secretory / ER",
    "SEC61B": "B-cell secretory / ER", "RPN2": "B-cell secretory / ER",
    "ERLEC1": "B-cell secretory / ER",
    "IGHM": "B-cell secretory / ER", "TMED10": "B-cell secretory / ER",
    "IGHG1": "Ig constant", "IGHG3": "Ig constant", "IGHA1": "Ig constant",
    "IGHA2": "Ig constant", "IGLC1": "Ig constant", "IGKC": "Ig constant",
    "JCHAIN": "Ig constant",
    "FTL": "Iron metabolism", "FTH1": "Iron metabolism",
    "HLA-DPB1": "MHC II / GWAS", "HLA-DRA": "MHC II / GWAS",
    "ANKRD55": "MHC II / GWAS", "CD58": "MHC II / GWAS",
    "IL7R": "MHC II / GWAS", "TYK2": "MHC II / GWAS",
    "KLRB1": "Cytotoxic", "KLRC1": "Cytotoxic", "KLRD1": "Cytotoxic",
    "CCL5": "Cytotoxic", "CCL4": "Cytotoxic", "GZMA": "Cytotoxic",
    "GZMB": "Cytotoxic", "GZMH": "Cytotoxic", "GZMK": "Cytotoxic",
    "PRF1": "Cytotoxic", "NKG7": "Cytotoxic", "GNLY": "Cytotoxic",
    "FGFBP2": "Cytotoxic",
    "ISG15": "Type I IFN", "MX1": "Type I IFN", "IFI6": "Type I IFN",
    "IFIT3": "Type I IFN", "OAS1": "Type I IFN",
    "FOXP3": "TF / regulator", "EOMES": "TF / regulator", "ZNF683": "TF / regulator",
    "ETV6": "TF / regulator", "HAVCR2": "TF / regulator",
}
FIG4_CT_LABEL = {"B": "B", "Mono": "Mono", "CD4_T": "CD4 T", "CD8_T": "CD8 T",
                 "NK": "NK", "DC": "DC", "dnT": "dnT", "gdT": "γδ T"}


def figure4():
    """Figure 4: quboFS selection-frequency heatmap (sequential Blues, 0–100%).
    Rows = 8 cell types, columns = top-5 genes per cell type (union). Red boxes
    and bold NPG-coloured x-labels mark genes mapped to a curated MS-relevant
    biology category. Selection frequency = % of that cell type's evaluable
    cohort-by-fold panels in which quboFS selected the gene."""
    import glob
    from matplotlib.patches import Rectangle
    fs = sorted(glob.glob(str(RUN / TAG / "CSF" / "selected_genes_folds_*.csv"))
                + glob.glob(str(RUN / (TAG + "_holdout_*") / "CSF" / "selected_genes_folds_*.csv")))
    parts = []
    for f in fs:
        d = pd.read_csv(f); d["src"] = f
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    q = df[df["method"] == "QUBO"].copy()
    cts = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
    # panels per cell type = unique (file, fold) pairs (3 cohorts × 5 folds = up to 15)
    npan = {ct: max(q[q.cell_type == ct].drop_duplicates(["src", "fold"]).shape[0], 1)
            for ct in cts}
    cnt = q.groupby(["cell_type", "gene"]).size().rename("n").reset_index()
    cnt["pct"] = cnt.apply(lambda r: r["n"] / npan[r["cell_type"]] * 100.0, axis=1)
    top, cols = {}, []
    for ct in cts:
        g = (cnt[cnt.cell_type == ct].sort_values(["pct", "gene"], ascending=[False, True])
             ["gene"].head(5).tolist())
        top[ct] = g
        cols += g
    cols = list(dict.fromkeys(cols))
    M = cnt.pivot_table(index="cell_type", columns="gene", values="pct", fill_value=0.0)
    M = M.reindex(index=cts, columns=cols, fill_value=0.0)
    n_ct, n_g = len(cts), len(cols)

    fig, ax = plt.subplots(figsize=(max(10.5, 0.33 * n_g + 2.8), 5.0))
    im = ax.imshow(M.values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    # annotate cell values
    for i in range(n_ct):
        for j in range(n_g):
            v = M.values[i, j]
            if v > 0:
                ax.text(j, i, f"{int(round(v))}", ha="center", va="center",
                        fontsize=6.5, color="white" if v >= 50 else "#1f3a5f")
    # red boxes around biology-mapped gene columns
    for j, g in enumerate(cols):
        if g in FIG4_BIOLOGY:
            for i in range(n_ct):
                if M.values[i, j] > 0:
                    ax.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, fill=False,
                                 edgecolor="#E64B35", linewidth=1.0))
    ax.set_yticks(range(n_ct))
    ax.set_yticklabels([f"{FIG4_CT_LABEL[ct]} (n={npan[ct]})" for ct in cts], fontsize=9)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(cols, rotation=45, ha="right", rotation_mode="anchor", fontsize=7.2)
    for tlbl in ax.get_xticklabels():
        g = tlbl.get_text()
        if g in FIG4_BIOLOGY:
            tlbl.set_fontweight("bold"); tlbl.set_color(FIG4_NPG[FIG4_BIOLOGY[g]])
    ax.set_ylabel("Cell type (n = panels)", fontsize=9, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_color("#555555"); spine.set_linewidth(0.7)
    cb = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02, shrink=0.85)
    cb.set_label("Selection frequency (%)", fontsize=8, labelpad=4)
    cb.ax.tick_params(labelsize=7)
    ax.set_title("quboFS selection frequency: top 5 genes per cell type "
                 "(CSF, 3 held-out cohorts × 5 folds = up to 15 panels)",
                 fontsize=10, fontweight="bold", pad=8)
    fig.tight_layout()
    fig.savefig(OUT / "figure4_biology.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "figure4_biology.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figure4_biology (Blues transposed: {n_ct} cell types × {n_g} genes; "
          f"{sum(1 for g in cols if g in FIG4_BIOLOGY)} biology-mapped)")


def figureS2_ksweep():
    """Panel-size (K) sweep (Figure S2): ROC-AUC and within-panel |rho| vs K per
    method. Skips gracefully if the sweep summary has not been generated yet."""
    path = RUN / "sweep_all_methods_K_summary.csv"
    if not path.exists():
        print("skip figureS2_K_sweep (run sweep_collect.py first)")
        return
    s = pd.read_csv(path)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (col, lab, better) in zip(
            axes, [("roc_auc", "Held-out ROC-AUC", "higher = better"),
                   ("macro_f1", "Held-out Macro-F1", "higher = better"),
                   ("within_panel_rho", "within-panel |rho|", "lower = better")]):
        for m in METHODS:
            sm = s[s["method"] == m].sort_values("K")
            if sm.empty:
                continue
            q = (m == "QUBO")
            ax.plot(sm["K"], sm[col], "-o", color=COL[m], label=LABEL[m],
                    markersize=7 if q else 4, linewidth=2.8 if q else 1.3,
                    alpha=1.0 if q else 0.75, zorder=5 if q else 2)
        ax.set_xlabel("panel size K"); ax.set_ylabel(lab)
        ax.set_title(f"{lab}  ({better})", fontsize=11)
    axes[1].legend(ncol=6, fontsize=8, loc="lower center", frameon=False,
                   bbox_to_anchor=(0.5, -0.34))
    fig.suptitle("Panel-size sensitivity", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "figureS2_K_sweep.png")
    fig.savefig(OUT / "figureS2_K_sweep.pdf")
    plt.close(fig)
    print("wrote figureS2_K_sweep")


def figureS3_solver():
    """Solver sensitivity (Figure S3): SA vs exact QUBO optimum (Jaccard + energy
    gap). Skips gracefully if the summary has not been generated yet."""
    path = RUN / "solver_sensitivity_summary.csv"
    if not path.exists():
        print("skip figureS3_solver (run solver_sensitivity.py first)")
        return
    d = pd.read_csv(path)
    pct = 100 * d["identical"].mean()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))
    axes[0].hist(d["jaccard_sa_exact"], bins=np.linspace(0, 1, 21),
                 color="#0072B2", edgecolor="black", linewidth=0.4)
    axes[0].set_xlabel("Jaccard(SA panel, exact-optimum panel)")
    axes[0].set_ylabel("panels"); axes[0].set_title("Panel agreement", fontsize=11)
    axes[1].hist(d["energy_gap"].abs(), bins=20, color="#D55E00",
                 edgecolor="black", linewidth=0.4)
    axes[1].set_xlabel("|energy gap|  (SA − exact)")
    axes[1].set_ylabel("panels"); axes[1].set_title("Energy gap to global optimum", fontsize=11)
    fig.suptitle("Solver sensitivity", fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "figureS3_solver_sensitivity.png")
    fig.savefig(OUT / "figureS3_solver_sensitivity.pdf")
    plt.close(fig)
    print("wrote figureS3_solver_sensitivity")


RAMESH_SIG = ["IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGKC", "IGLC1", "IGLC2",
              "IGLC3", "XBP1", "MZB1", "JCHAIN", "PRDM1", "CD27", "CD80", "CD86",
              "NFKB1", "NFKB2", "REL", "RELA", "RELB", "NFKBIA", "NFKBIE",
              "TNFRSF13B", "TNFRSF17", "TNFRSF13C", "AICDA"]


def figureS4_recovery():
    """Biological concordance (Figure S4):
    A. Ramesh (2020) B-cell signature recovery per method (canonical selected genes);
    B. per-cell-type ΔAUCell (MS − control) of the quboFS panels — drawn only if the
       06_aucell output (qubo_run/aucell_results/ms_vs_hd_diff_qubo.csv) is present
       (that step needs the full Seurat .rds)."""
    import glob
    fs = sorted(glob.glob(str(RUN / TAG / "CSF" / "selected_genes_folds_*.csv"))
                + glob.glob(str(RUN / (TAG+"_holdout_*") / "CSF" / "selected_genes_folds_*.csv")))
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    b = df[df["cell_type"] == "B"]
    surv = [g for g in RAMESH_SIG if g in set(b["gene"])]
    methods = ["QUBO", "mRMR", "DE_top", "LASSO", "ElasticNet", "HVG"]
    rec = {m: sum(1 for g in surv if g in set(b[b.method == m]["gene"])) for m in methods}
    aucell = RUN / "aucell_results" / "ms_vs_hd_diff_qubo.csv"
    have_b = aucell.exists()
    fig, axes = plt.subplots(1, 2 if have_b else 1, figsize=(11 if have_b else 6.2, 4.4))
    axA = axes[0] if have_b else axes
    x = np.arange(len(methods))
    axA.bar(x, [rec[m] for m in methods], 0.6,
            color=[COL[m] for m in methods], edgecolor="black", linewidth=0.5)
    for i, m in enumerate(methods):
        axA.text(i, rec[m] + 0.1, str(rec[m]), ha="center", va="bottom", fontsize=9)
    axA.set_xticks(x); axA.set_xticklabels([LABEL[m] for m in methods], rotation=25, ha="right")
    axA.set_ylabel(f"Ramesh B-cell signature genes recovered (/{len(surv)})")
    axA.set_ylim(0, len(surv) + 1)
    axA.set_title("A. B-cell signature recovery (Ramesh 2020)", fontsize=11, fontweight="bold")
    if have_b:
        d = pd.read_csv(aucell)
        ccol = "cell_type" if "cell_type" in d.columns else d.columns[0]
        vcol = "mean_diff_MS_minus_HD" if "mean_diff_MS_minus_HD" in d.columns else d.columns[1]
        cts = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
        d = d.set_index(ccol).reindex(cts)
        vals = d[vcol].values
        colors = ["#D55E00" if ct == "B" else "#9ecae1" for ct in cts]
        axes[1].bar(range(len(cts)), vals, 0.6, color=colors, edgecolor="black", linewidth=0.5)
        axes[1].axhline(0, color="black", lw=0.8)
        axes[1].set_xticks(range(len(cts))); axes[1].set_xticklabels(cts, rotation=30, ha="right")
        axes[1].set_ylabel("ΔAUCell (MS − control)")
        axes[1].set_title("B. quboFS panel activity (AUCell)", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "figureS4_literature_concordance.png")
    fig.savefig(OUT / "figureS4_literature_concordance.pdf")
    plt.close(fig)
    print(f"wrote figureS4_literature_concordance (recovery; AUCell panel={'yes' if have_b else 'pending .rds run'})")


def _read_mm(path):
    """Minimal MatrixMarket coordinate reader -> dense ndarray (rows x cols)."""
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while lines[i].startswith("%"):
        i += 1
    nr, nc, _ = map(int, lines[i].split()); i += 1
    M = np.zeros((nr, nc))
    for ln in lines[i:]:
        r, c, v = ln.split(); M[int(r) - 1, int(c) - 1] = float(v)
    return M


def figureS6_dotplot():
    """B-cell secretory / ER-chaperone programme dot plot across cohorts (Figure S6),
    built from the donor-level B-cell log-normalised pseudobulk (no .rds needed).
    Dot size = fraction of donors expressing; colour = mean expression."""
    DATA = PROJ / "data" / os.environ.get("QUBOFS_PSEUDOBULK_SUBDIR", "pseudobulk_v5_compartment")
    base = DATA / "B" / "CSF" / "fold_1"
    if not base.exists():
        print("skip figureS6_bcell_dotplot (B-cell pseudobulk not found)")
        return
    # combine train + val + heldout (fold_1) -> all donors once
    donors, expr_cols, genes_ref = [], [], None
    meta_rows = []
    for setname in ["train", "val", "heldout"]:
        mtx = base / f"{setname}_pb_mean.mtx"
        if not mtx.exists():
            continue
        M = _read_mm(str(mtx))  # genes x donors
        gr = pd.read_csv(base / f"{setname}_pb_mean_rows.csv")["gene"].tolist()
        mc = pd.read_csv(base / f"{setname}_pb_mean_cols.csv")
        meta = pd.read_csv(base / f"{setname}_meta.csv")
        genes_ref = gr if genes_ref is None else genes_ref
        gpos = {g: i for i, g in enumerate(gr)}
        for j, don in enumerate(meta["donor_id"]):
            donors.append(don); expr_cols.append((M[:, j], gpos))
            meta_rows.append(meta.iloc[j])
    meta = pd.DataFrame(meta_rows).reset_index(drop=True)
    genes = ["FKBP11", "TMED10", "LMAN2", "DNAJB11", "XBP1", "MZB1",
             "IGHG1", "IGHA1", "IGKC", "JCHAIN"]
    cohort_name = {"PRJNA671484_MS_Tcell": "Pappalardo", "osmzhlab_MS_ence_cov": "Heming",
                   "PRJNA549712_MS_PBMC_UCSF": "Ramesh", "PRJNA979258_cryoCSF": "Touil"}
    meta["coh"] = meta["cohort"].map(lambda c: cohort_name.get(c, c))
    groups = [("Pappalardo", "HD"), ("Pappalardo", "MS"), ("Heming", "HD"), ("Heming", "MS"),
              ("Ramesh", "HD"), ("Ramesh", "MS"), ("Touil", "HD")]
    # expression matrix donors x genes
    E = np.zeros((len(donors), len(genes)))
    for di, (vec, gpos) in enumerate(expr_cols):
        for gi, g in enumerate(genes):
            E[di, gi] = vec[gpos[g]] if g in gpos else 0.0
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ylabs = []
    for yi, (coh, dx) in enumerate(groups):
        idx = meta.index[(meta["coh"] == coh) & (meta["diagnosis"] == dx)].tolist()
        ylabs.append(f"{coh} {dx} (n={len(idx)})")
        if not idx:
            continue
        sub = E[idx, :]
        for gi in range(len(genes)):
            frac = (sub[:, gi] > 0).mean()
            mean = sub[:, gi].mean()
            ax.scatter(gi, yi, s=20 + frac * 260,
                       c=[mean], cmap="Blues", vmin=0, vmax=np.percentile(E, 99),
                       edgecolor="black", linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=45, ha="right",
                                                         fontstyle="italic", fontsize=9)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels(ylabs, fontsize=9)
    ax.invert_yaxis()
    ax.set_title("B-cell secretory / ER-chaperone programme (donor pseudobulk)", fontsize=10.5)
    sm = plt.cm.ScalarMappable(cmap="Blues",
                               norm=plt.Normalize(0, np.percentile(E, 99)))
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02); cb.set_label("mean log-norm expr.")
    # size legend
    for fr, lab in [(0.25, "25%"), (0.5, "50%"), (1.0, "100%")]:
        ax.scatter([], [], s=20 + fr * 260, c="grey", edgecolor="black",
                   linewidth=0.4, label=lab)
    ax.legend(title="% donors expr.", fontsize=7, title_fontsize=7,
              loc="lower right", frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "figureS6_bcell_dotplot.png")
    fig.savefig(OUT / "figureS6_bcell_dotplot.pdf")
    plt.close(fig)
    print(f"wrote figureS6_bcell_dotplot ({len(donors)} donors)")


def figureS7_audit():
    """Manual gene-audit category distribution per cell type (Figure S7).
    Reads the curated audit file qubo_run/audit_categories.csv (cell_type, gene,
    category) — filled in by the annotators from data_release/exploratory_gene_audit.csv.
    Skips if the curated file is not present."""
    path = RUN / "audit_categories.csv"
    if not path.exists():
        print("skip figureS7_audit (provide qubo_run/audit_categories.csv from the curated audit)")
        return
    d = pd.read_csv(path)
    d = d[d["category"].astype(str).str.len() > 0]
    cts = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
    cats = sorted(d["category"].unique())
    pal = plt.cm.tab20(np.linspace(0, 1, len(cats)))
    frac = (d.groupby(["cell_type", "category"]).size()
            .groupby(level=0).apply(lambda s: s / s.sum()).unstack(fill_value=0))
    frac = frac.reindex(index=[c for c in cts if c in frac.index], columns=cats, fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 4.6))
    bottom = np.zeros(len(frac))
    for k, cat in enumerate(cats):
        ax.bar(range(len(frac)), frac[cat].values, 0.7, bottom=bottom,
               color=pal[k], edgecolor="black", linewidth=0.3, label=cat)
        bottom += frac[cat].values
    ax.set_xticks(range(len(frac))); ax.set_xticklabels(frac.index, rotation=30, ha="right")
    ax.set_ylabel("fraction of selected genes"); ax.set_ylim(0, 1)
    ax.set_title("Gene-audit category distribution per cell type", fontsize=11)
    ax.legend(fontsize=7, ncol=2, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "figureS7_audit_categories.png")
    fig.savefig(OUT / "figureS7_audit_categories.pdf")
    plt.close(fig)
    print("wrote figureS7_audit_categories")


if __name__ == "__main__":
    # NOTE: figure1() (workflow schematic) is intentionally NOT called here.
    # Figure 1 is the curated graphic in figures_oup/figure1_pipeline.png and must
    # NOT be overwritten by the plain matplotlib version. The figure1() function is
    # kept only as a fallback/reference.
    figure2()
    figure3()
    figure4()
    figureS1()           # Figure S1: dataset composition
    figureS2_ksweep()    # Figure S2: panel-size (K) sensitivity
    figureS3_solver()    # Figure S3: solver sensitivity
    figureS4_recovery()  # Figure S4: literature concordance / Ramesh B-cell signature
    # figureS6_dotplot() and figureS7_audit() are retained as optional development
    # utilities but are NOT part of the submitted manuscript supplement: the B-cell
    # dot plot and the audit-category figure were removed to keep it focused (the
    # audit is released as a CSV).
