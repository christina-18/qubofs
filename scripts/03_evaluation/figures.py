"""Figures of the CSF evaluation (reads the outputs of csf_performance.py, csf_redundancy.py, csf_genes.py and
sensitivity_noIgV_Hb.py, and the training |rho| of the benchmark models). Methods are labelled as text, never by colour alone.
Output: $QUBOFS_RESULTS_DIR/evaluation/figures/ (supporting figures; PNG 200 dpi and PDF) and the manuscript
Figures 2 and 3 in <repository>/figures_oup/ (or $QUBOFS_FIGURE_DIR): figure2_performance_redundancy and
figure3_gene_frequency, as PDF, SVG and 600-dpi PNG, with the plotted values as CSV."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIG = out("figures"); PERF = EVAL / "csf_performance"; RED = EVAL / "csf_redundancy"
SEN = pd.read_csv(EVAL / "sensitivity_noIgV_Hb" / "auc_ci.csv"); SEN = SEN[SEN.scheme == "weighted"]
COL = {"QUBO": "#2a78d6", "ElasticNet": "#eb6834", "LASSO": "#1baf7a", "mRMR": "#eda100", "DE_top": "#e87ba4"}
MK = {"QUBO": "o", "ElasticNet": "s", "LASSO": "^", "mRMR": "D", "DE_top": "v"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
COH = HOS + ["Mean"]; CNAME = {"Pappalardo": "Pappalardo", "Heming": "Heming", "Ramesh": "Ramesh (3 controls)", "Mean": "3-cohort mean"}
def save(f, n): f.savefig(FIG / f"{n}.png", dpi=200, bbox_inches="tight"); f.savefig(FIG / f"{n}.pdf", bbox_inches="tight"); plt.close(f)
def grid(ax): ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)

# Fig A1: per-cohort ROC-AUC and PR-AUC with 95% CI
L = pd.read_csv(PERF / "metrics_ci_weighted.csv")
T = L.pivot_table(index=["cohort", "method"], columns="metric", values=["value", "lo", "hi"])
T.columns = [f"{m}_auc" if s == "value" else f"{m}_{s}" for s, m in T.columns]; T = T.reset_index()
DP = pd.read_csv(PERF / "donor_predictions_weighted.csv")
prev = {h: (DP[(DP.holdout == h) & (DP.outer == 1) & (DP.method == "QUBO")].diagnosis == "MS").mean() for h in HOS}
f, axs = plt.subplots(2, 4, figsize=(9.5, 4.2), sharey=True)
for r, (met, lab) in enumerate([("roc", "ROC-AUC"), ("pr", "PR-AUC")]):
    for c, coh in enumerate(COH):
        ax = axs[r, c]; grid(ax)
        for i, m in enumerate(METHODS):
            x = T[(T.cohort == coh) & (T.method == m)].iloc[0]; y = len(METHODS) - 1 - i
            ax.plot([x[f"{met}_lo"], x[f"{met}_hi"]], [y, y], color=COL[m], lw=2, solid_capstyle="round")
            ax.plot(x[f"{met}_auc"], y, MK[m], color=COL[m], ms=6, mec="white", mew=1)
        if met == "roc": ax.axvline(0.5, color=INK2, lw=0.8, ls=":")
        elif coh in prev: ax.axvline(prev[coh], color=INK2, lw=0.8, ls=":")
        ax.set_xlim(0.45, 1.01); ax.set_yticks(range(len(METHODS))); ax.set_yticklabels([LABEL[m] for m in METHODS][::-1])
        if r == 0: ax.set_title(CNAME[coh], fontsize=9, color=INK)
        ax.set_xlabel(lab)
f.suptitle("Held-out performance per external cohort (point = mean of 5 outer-fold models; bar = 95% donor-bootstrap CI)\n"
           "dotted line: chance (ROC) or MS prevalence (PR)", fontsize=8.5, color=INK2, y=1.04)
f.tight_layout(); save(f, "figA1_auc_per_cohort")

# Fig A2: paired differences quboFS - comparator
Dd = pd.read_csv(PERF / "paired_diffs_weighted.csv")
f, axs = plt.subplots(2, 4, figsize=(9.5, 3.8), sharey=True)
for r, (met, lab) in enumerate([("roc", "Δ ROC-AUC (quboFS − comparator)"), ("pr", "Δ PR-AUC (quboFS − comparator)")]):
    for c, coh in enumerate(COH):
        ax = axs[r, c]; grid(ax); ax.axvline(0, color=INK, lw=0.8)
        for i, m in enumerate(COMPARATORS):
            x = Dd[(Dd.cohort == coh) & (Dd.comparator == m) & (Dd.metric == met)].iloc[0]; y = len(COMPARATORS) - 1 - i
            ax.plot([x.lo, x.hi], [y, y], color=COL[m], lw=2, solid_capstyle="round"); ax.plot(x["diff"], y, MK[m], color=COL[m], ms=6, mec="white", mew=1)
        ax.set_yticks(range(len(COMPARATORS))); ax.set_yticklabels([f"vs {LABEL[m]}" for m in COMPARATORS][::-1])
        ax.set_xlim(-0.2, 0.26)
        if r == 0: ax.set_title(CNAME[coh], fontsize=9, color=INK)
        ax.set_xlabel(lab, fontsize=7.5)
f.suptitle("Paired differences, same donor resamples for all methods (95% percentile CI). Ramesh CIs rest on 3 controls.",
           fontsize=8.5, color=INK2, y=1.02)
f.tight_layout(); save(f, "figA2_paired_diffs")

# Fig B: redundancy
Pc = pd.read_csv(RED / "train_rho_per_celltype.csv"); Tr = pd.read_csv(RED / "train_rho_paired.csv"); Hr = pd.read_csv(RED / "heldout_rho_ci.csv")
f, axs = plt.subplots(1, 2, figsize=(9.5, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axs[0]; grid(ax); ax.axvline(0, color=INK, lw=0.8)
for i, m in enumerate(COMPARATORS):
    y = len(COMPARATORS) - 1 - i; g = Pc[Pc.comparator == m]; t = Tr[Tr.comparator == m].iloc[0]
    ax.scatter(g.mean_d, y + np.linspace(-0.22, 0.22, len(g)), s=10, color=COL[m], alpha=0.55, lw=0)
    ax.plot([t.lo, t.hi], [y - 0.35, y - 0.35], color=INK, lw=2, solid_capstyle="round"); ax.plot(t.diff_ct_mean, y - 0.35, "|", color=INK, ms=8)
ax.set_yticks(range(len(COMPARATORS))); ax.set_yticklabels([f"vs {LABEL[m]}" for m in COMPARATORS][::-1])
ax.set_xlabel("Δ mean |ρ| on training donors (quboFS − comparator)"); ax.set_title("Training |ρ|: dots = 8 cell types; black = mean, hierarchical 95% CI", fontsize=8, color=INK)
ax = axs[1]; grid(ax); ax.axvline(0, color=INK, lw=0.8); off = {h: o for h, o in zip(COH, [0.27, 0.09, -0.09, -0.27])}
for i, m in enumerate(COMPARATORS):
    y = len(COMPARATORS) - 1 - i
    for coh in COH:
        x = Hr[(Hr.cohort == coh) & (Hr.method == f"diff_vs_{m}")].iloc[0]; yy = y + off[coh]
        ax.plot([x.lo, x.hi], [yy, yy], color=INK if coh == "Mean" else INK2, lw=2 if coh == "Mean" else 1.2)
        ax.plot(x.rho_heldout, yy, "o" if coh == "Mean" else ".", color=INK if coh == "Mean" else INK2, ms=5)
        if i == 0: ax.annotate(CNAME[coh].split(" (")[0], (x.hi, yy), xytext=(4, 0), textcoords="offset points", va="center", fontsize=6.5, color=INK2)
ax.set_yticks(range(len(COMPARATORS))); ax.set_yticklabels([f"vs {LABEL[m]}" for m in COMPARATORS][::-1])
ax.set_xlabel("Δ mean |ρ| on held-out donors (quboFS − comparator)"); ax.set_title("Held-out |ρ| of the same panels, per cohort (donor-bootstrap 95% CI)", fontsize=8, color=INK)
f.tight_layout(); save(f, "figB_redundancy")

# Fig A3: Ramesh per-control
R = pd.read_csv(PERF / "ramesh_controls.csv"); f, axs = plt.subplots(1, 2, figsize=(8, 2.6), sharey=True)
for ax, an, tt in [(axs[0], "per_control_AUC", "Share of 12 MS donors ranked above each control"), (axs[1], "leave_one_control_out", "AUC with that control left out (2 controls)")]:
    grid(ax); g = R[R.analysis == an]; ctrls = sorted(g.control.unique())
    for i, m in enumerate(METHODS):
        y = len(METHODS) - 1 - i
        for j, c in enumerate(ctrls):
            x = g[(g.method == m) & (g.control == c)].iloc[0]
            ax.plot(x.value, y + (j - 1) * 0.22, MK[m], color=COL[m], ms=5, mec="white", mew=0.8)
            if i == 0: ax.annotate(c.replace("Ramesh_", ""), (x.value, y + (j - 1) * 0.22), xytext=(4, 0), textcoords="offset points", fontsize=6, color=INK2, va="center")
    ax.set_yticks(range(len(METHODS))); ax.set_yticklabels([LABEL[m] for m in METHODS][::-1]); ax.set_xlim(0.55, 1.02)
    ax.set_title(tt, fontsize=8, color=INK); ax.set_xlabel("AUC (mean of 5 outer-fold models)")
f.tight_layout(); save(f, "figA3_ramesh_controls")

# Fig C: correlation modules
M = pd.read_csv(EVAL / "csf_genes" / "modules_panels.csv").groupby(["tau", "method"]).same_module_pairs.mean().unstack()
f, ax = plt.subplots(figsize=(4.2, 2.8)); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
for m in METHODS:
    ax.plot(M.index, M[m], marker=MK[m], color=COL[m], lw=2, ms=6, mec="white", mew=1)
    ax.annotate(LABEL[m], (M.index[-1], M[m].iloc[-1]), xytext=(6, 0), textcoords="offset points", va="center", fontsize=7, color=INK)
ax.set_xlabel("module threshold τ (|ρ| ≥ τ, complete linkage)"); ax.set_ylabel("share of panel gene pairs\nin the same module")
ax.set_title("Correlation modules (descriptive; all donors, within-cohort centred)", fontsize=8, color=INK); ax.set_xticks([0.4, 0.5, 0.6, 0.7])
f.tight_layout(); save(f, "figC_modules")

# Fig S: sensitivity analysis without Ig V/J/D and haemoglobin genes
f, axs = plt.subplots(1, 4, figsize=(9.5, 2.9), sharey=True)
for ax, coh in zip(axs, HOS + ["Mean"]):
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for i, m in enumerate(METHODS):
        y = len(METHODS) - 1 - i
        for dy, fi, mk, fill in [(0.16, CFG.PRIMARY, "o", True), (-0.16, "edger_noIgV_Hb", "o", False)]:
            x = SEN[(SEN["filter"] == fi) & (SEN.cohort == coh) & (SEN.method == m)].iloc[0]
            ax.plot([x.roc_lo, x.roc_hi], [y + dy, y + dy], color=COL[m], lw=1.6, alpha=1 if fill else 0.6)
            ax.plot(x.roc, y + dy, mk, ms=5.5, color=COL[m], mfc=COL[m] if fill else "white", mew=1.4)
    ax.set_yticks(range(len(METHODS))); ax.set_yticklabels([LABEL[m] for m in METHODS][::-1]); ax.set_xlim(0.5, 1.01)
    ax.set_title({"Mean": "3-cohort mean", "Ramesh": "Ramesh (3 controls)"}.get(coh, coh), fontsize=9); ax.set_xlabel("ROC-AUC")
axs[0].plot([], [], "o", color=INK, label="primary (edgeR only)"); axs[0].plot([], [], "o", mfc="white", color=INK, label="Ig V(D)J + Hb excluded")
f.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.06))
f.tight_layout(); save(f, "figS_sensitivity_noIgV_Hb")
print(sorted(p.name for p in FIG.glob("*.png")))


# ---- Manuscript Figure 2: gene-panel redundancy (A) and held-out ROC-AUC (B), primary analysis (edger_only) ----
import os  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

FIG2 = Path(os.environ.get("QUBOFS_FIGURE_DIR", Path(__file__).resolve().parents[2] / "figures_oup")); FIG2.mkdir(parents=True, exist_ok=True)
ORDER = ["QUBO", "mRMR", "LASSO", "ElasticNet", "DE_top"]
CT_LABEL = {"B": "B", "Mono": "Mono", "CD4_T": "CD4 T", "CD8_T": "CD8 T", "NK": "NK", "DC": "DC", "dnT": "dnT", "gdT": "γδT"}
# Panel A: mean training |rho| per method x cell type over the 15 outer models (rho.csv of each model)
RH = pd.concat([pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "rho.csv", keep_default_na=False).assign(holdout=ho, outer=k)
                for ho in HOS for k in OUTERS])
assert RH.groupby(["method", "cell_type"]).size().eq(15).all(), "every method x cell type must have 15 outer models"
A = RH.pivot_table(index="method", columns="cell_type", values="rho", aggfunc="mean").loc[ORDER, CFG.CELL_TYPES]
A["Mean"] = A[CFG.CELL_TYPES].mean(axis=1)
# Panel B: weighted-ensemble ROC-AUC per held-out cohort (mean of the 5 outer models) and 3-cohort mean
MB = pd.read_csv(PERF / "metrics_ci_weighted.csv"); MB = MB[MB.metric == "roc"]
Bm = MB.pivot_table(index="method", columns="cohort", values="value").loc[ORDER, HOS + ["Mean"]]
# agreement with the evaluation outputs
TRP = pd.read_csv(RED / "train_rho_paired.csv")
assert np.allclose(A.loc["QUBO", "Mean"], TRP.mean_quboFS.iloc[0], atol=1e-12)
for _, r in TRP.iterrows():
    assert np.isclose(A.loc[r.comparator, "Mean"], r.mean_comp, atol=1e-12)
assert np.allclose(Bm["Mean"], Bm[HOS].mean(axis=1), atol=1e-12)
print(f"Figure 2 check: quboFS ROC-AUC {Bm.loc['QUBO', 'Mean']:.3f}, mean within-panel |rho| {A.loc['QUBO', 'Mean']:.3f}")


def heat(ax, M, cmap, norm, fmt, title, xlabels, cbar_label, fig):
    """Heatmap with the 'Mean' column set apart by a gap; values printed in every cell."""
    cols = list(M.columns); n = len(cols) - 1
    xpos = list(range(n)) + [n + 0.35]
    for j, (c, x) in enumerate(zip(cols, xpos)):
        for i, m in enumerate(M.index):
            v = M.loc[m, c]; rgba = cmap(norm(v))
            ax.add_patch(plt.Rectangle((x - 0.5, i - 0.5), 1, 1, facecolor=rgba, edgecolor="white", lw=1.2))
            lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            ax.text(x, i, fmt.format(v), ha="center", va="center", fontsize=6.5, color="white" if lum < 0.5 else INK,
                    fontweight="bold" if m == "QUBO" else "normal")
    ax.set_xlim(-0.5, xpos[-1] + 0.5); ax.set_ylim(len(M) - 0.5, -0.5)
    ax.set_xticks(xpos); ax.set_xticklabels(xlabels, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(M))); ax.set_yticklabels([LABEL[m] for m in M.index])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0); ax.set_title(title, fontsize=8.5, loc="left", color=INK)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label(cbar_label, fontsize=7.5); cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)


f, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [9.4, 5.0]})
heat(axs[0], A, plt.get_cmap("YlOrRd"), Normalize(0.15, 0.50), "{:.3f}",
     "A  Within-panel redundancy (training donors)", [CT_LABEL[c] for c in CFG.CELL_TYPES] + ["Mean"], "mean |Pearson r|", f)
heat(axs[1], Bm, plt.get_cmap("Blues"), Normalize(0.70, 1.00), "{:.3f}",
     "B  Held-out ROC-AUC", ["Pappalardo", "Heming", "Ramesh", "Mean"], "ROC-AUC", f)
axs[1].set_yticklabels([])
f.tight_layout()
for ext, kw in [("pdf", {}), ("svg", {}), ("png", {"dpi": 600})]:
    f.savefig(FIG2 / f"figure2_performance_redundancy.{ext}", bbox_inches="tight", **kw)
plt.close(f)
A.to_csv(FIG2 / "figure2_panelA_redundancy.csv"); Bm.to_csv(FIG2 / "figure2_panelB_roc_auc.csv")
print("Figure 2 ->", FIG2)


# ---- Manuscript Figure 3: selection frequency of recurrent quboFS genes per cell type (CSF, 15 outer models) ----
# Columns: union over cell types of the 5 genes most often selected by quboFS for that cell type (ties broken
# alphabetically); a gene in the top 5 of several cell types appears once, in the first block. Cells: % of the 15
# outer-fold models in which quboFS selected the gene for that cell type (white = 0%, no grid lines).
SF = pd.read_csv(EVAL / "csf_genes" / "selection_frequency.csv"); SF = SF[SF.method == "QUBO"].copy()
PQ = pd.concat([pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False).assign(model=f"{ho}{k}")
                for ho in HOS for k in OUTERS]); PQ = PQ[PQ.method == "QUBO"]
N_MODELS = len(HOS) * len(OUTERS)
assert PQ.groupby("cell_type").model.nunique().reindex(CFG.CELL_TYPES).eq(N_MODELS).all(), "every cell type must have 15 quboFS panels"
chk = PQ.groupby(["cell_type", "gene"]).size().rename("n").reset_index().merge(SF, on=["cell_type", "gene"])
assert len(chk) == len(SF) and (chk.n == chk.n_folds).all(), "selection_frequency.csv disagrees with the panels"
selected, block = [], []
for ct in CFG.CELL_TYPES:
    for g in SF[SF.cell_type == ct].sort_values(["n_folds", "gene"], ascending=[False, True]).head(5).gene:
        if g not in selected:
            selected.append(g); block.append(ct)
freq = SF.set_index(["cell_type", "gene"]).n_folds
values = np.array([[100 * freq.get((ct, g), 0) / N_MODELS for g in selected] for ct in CFG.CELL_TYPES])
for g, ct, n in [("P2RY14", "DC", 5), ("FKBP5", "Mono", 7), ("ISG15", "dnT", 8)]:
    assert np.isclose(values[CFG.CELL_TYPES.index(ct), selected.index(g)], 100 * n / N_MODELS), (g, ct)
pd.DataFrame(values.round(1), index=CFG.CELL_TYPES, columns=selected).to_csv(FIG2 / "figure3_gene_frequency.csv")
fig, ax = plt.subplots(figsize=(16, 5.5))
cmap3 = plt.get_cmap("Blues").copy(); cmap3.set_bad("white")
im = ax.imshow(np.ma.masked_where(values == 0, values), cmap=cmap3, norm=Normalize(vmin=0, vmax=60), aspect="auto", interpolation="nearest")
for i in range(values.shape[0]):
    for j in range(values.shape[1]):
        if values[i, j] > 0:
            ax.text(j, i, f"{values[i, j]:.0f}", ha="center", va="center", fontsize=8, color="white" if values[i, j] >= 35 else "black")
ax.set_xticks(np.arange(len(selected))); ax.set_xticklabels(selected, rotation=55, ha="right")
ax.set_yticks(np.arange(len(CFG.CELL_TYPES))); ax.set_yticklabels([CT_LABEL[c] for c in CFG.CELL_TYPES])
ax.grid(False); ax.set_facecolor("white"); ax.tick_params(axis="both", length=0)
for x in range(1, len(selected)):                 # separators between cell-type blocks
    if block[x] != block[x - 1]:
        ax.axvline(x - 0.5, color="0.35", linewidth=1.3)
for sp in ax.spines.values(): sp.set_visible(False)
cbar = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.025); cbar.set_label("% of 15 outer-fold models")
fig.tight_layout()
for ext, kw in [("pdf", {}), ("svg", {}), ("png", {"dpi": 600})]:
    fig.savefig(FIG2 / f"figure3_gene_frequency.{ext}", bbox_inches="tight", **kw)
plt.close(fig)
print(f"Figure 3 check: {len(selected)} genes; P2RY14 DC 5/15, FKBP5 Mono 7/15, ISG15 dnT 8/15")
