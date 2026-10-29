"""How the methods choose among the SAME candidates (CSF primary model, edger_only, 3 hold-outs x 5 outer folds).
Saved panels only; candidate pools, relevance and training correlations are rebuilt from the outer-fold training donors
(deterministic; no re-selection, no re-fitting).
A. Stability per method x cell type over the 15 fold models: mean pairwise Jaccard; Nogueira stability (Nogueira et al.,
   JMLR 2018; feature space d = union of the cell type's top-20 candidate pools over the 15 folds); module-level Jaccard
   (genes mapped to correlation modules, complete linkage at |rho| >= 0.5 on all donors centred within cohort; descriptive).
B. Per fold x cell type: relevance retained = sum of |z|*C of the panel / sum of the 10 highest |z|*C in the pool
   (= DE-top, by construction 1), mean |rho| (training donors), number of gene pairs with |r| > 0.7, overlap with DE-top.
   Network example fixed in advance: B cells, Pappalardo held out, outer fold 1.
C. Key genes (FKBP5-Mono, ISG15-dnT, P2RY14-DC, IL32-NK, LEF1-gdT, XBP1-B): folds selected per method and folds in the pool.
Output: $QUBOFS_RESULTS_DIR/evaluation/csf_stability/ (tables and figures)"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

O = out("csf_stability"); FIG = O
RAW = D.Raw(CFG.CSF_RAW)
P = pd.concat([pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False).assign(holdout=ho, outer=k)
               for ho in HOS for k in OUTERS])
rowsB, pools, exB = [], {}, None
for ho in HOS:
    for k in OUTERS:
        Bd, _, _ = csf_outer_bundles(RAW, ho, k)
        for ct, b in Bd.items():
            pool = S.candidates(b); pools.setdefault(ct, []).append(set(pool))
            rel = dict(zip(pool, S.relevance(b, pool))); top10 = sum(sorted(rel.values(), reverse=True)[:10])
            gp = {g: i for i, g in enumerate(b["train"]["genes"])}; X = b["train"]["X"][[gp[g] for g in pool]]
            with np.errstate(invalid="ignore"): Rm = pd.DataFrame(np.corrcoef(X), index=pool, columns=pool)
            de = set(P.gene[(P.holdout == ho) & (P.outer == k) & (P.cell_type == ct) & (P.method == "DE_top")])
            for m in METHODS:
                g = [x for x in P.gene[(P.holdout == ho) & (P.outer == k) & (P.cell_type == ct) & (P.method == m)] if x in rel]
                if len(g) < 2: continue
                sub = Rm.loc[g, g].abs().values; iu = np.triu_indices(len(g), 1)
                rowsB.append(dict(holdout=ho, outer=k, cell_type=ct, method=m, n=len(g), relevance_retained=sum(rel[x] for x in g) / top10,
                                  mean_rel_rank=np.mean([sorted(rel, key=rel.get, reverse=True).index(x) + 1 for x in g]),
                                  mean_abs_r=np.nanmean(sub[iu]), pairs_r_gt_0_7=int(np.nansum(sub[iu] > 0.7)), n_pairs=len(iu[0]),
                                  overlap_with_DEtop=len(set(g) & de)))
            if (ho, k, ct) == ("Pappalardo", 1, "B"):
                exB = dict(pool=pool, rel=rel, R=Rm, panels={m: list(P.gene[(P.holdout == ho) & (P.outer == k) & (P.cell_type == ct) & (P.method == m)]) for m in METHODS})
TB = pd.DataFrame(rowsB); TB.to_csv(O / "relevance_redundancy_per_fold_celltype.csv", index=False)
SB = TB.groupby("method")[["relevance_retained", "mean_rel_rank", "mean_abs_r", "pairs_r_gt_0_7", "overlap_with_DEtop"]].mean().loc[METHODS]
SB.to_csv(O / "relevance_redundancy_summary.csv")
# ---- A. stability ----
info = RAW.info; f0 = csf_folds("Pappalardo")
donors = [d for d in sorted(set(f0.donor)) if info.loc[d, "cohort"] != "PRJNA979258_cryoCSF"]
rowsA = []
for ct in CFG.CELL_TYPES:
    U = sorted(set().union(*pools[ct])); d = len(U)
    blk = RAW.block(ct, U, donors, "all"); X = pd.DataFrame(blk["X"], index=U, columns=blk["donors"]); cc = info.loc[blk["donors"], "cohort"].values
    for c in np.unique(cc): X.loc[:, cc == c] = X.loc[:, cc == c].sub(X.loc[:, cc == c].mean(1), axis=0)
    ok = X.std(1) > 0; Rm = np.abs(np.corrcoef(X[ok].values)); np.fill_diagonal(Rm, 1)
    lab = dict(zip(X.index[ok], fcluster(linkage(squareform(1 - Rm, checks=False), "complete"), 0.5, "distance")))
    for i, g in enumerate(X.index[~ok]): lab[g] = -(i + 1)
    for m in METHODS:
        PS = [set(x.gene) for _, x in P[(P.cell_type == ct) & (P.method == m)].groupby(["holdout", "outer"])]
        M = len(PS); kbar = np.mean([len(s) for s in PS])
        pf = np.array([sum(g in s for s in PS) / M for g in U])
        phi = 1 - np.mean(M / (M - 1) * pf * (1 - pf)) / ((kbar / d) * (1 - kbar / d))
        jac = np.mean([len(a & b) / len(a | b) for a, b in itertools.combinations(PS, 2)])
        Sm = [{lab[g] for g in s} for s in PS]; jm = np.mean([len(a & b) / len(a | b) for a, b in itertools.combinations(Sm, 2)])
        rowsA.append(dict(cell_type=ct, method=m, d_pool_union=d, M=M, kbar=kbar, nogueira=phi, jaccard_gene=jac, jaccard_module=jm,
                          n_distinct_genes=int((pf > 0).sum()), n_modules_used=len(set().union(*Sm))))
TA = pd.DataFrame(rowsA); TA.to_csv(O / "stability.csv", index=False)
SA = TA.groupby("method")[["nogueira", "jaccard_gene", "jaccard_module", "n_distinct_genes", "n_modules_used"]].mean().loc[METHODS]
SA.to_csv(O / "stability_summary.csv")
# ---- C. key genes ----
KG = [("Mono", "FKBP5"), ("dnT", "ISG15"), ("DC", "P2RY14"), ("NK", "IL32"), ("gdT", "LEF1"), ("B", "XBP1")]
rowsC = []
for ct, g in KG:
    r = dict(cell_type=ct, gene=g, in_pool_folds=sum(g in s for s in pools[ct]))
    for m in METHODS: r[m] = int(((P.cell_type == ct) & (P.gene == g) & (P.method == m)).sum())
    rowsC.append(r)
TC = pd.DataFrame(rowsC); TC.to_csv(O / "key_genes.csv", index=False)
# quboFS genes recurrent (>=3 folds) that DE-top never or rarely picks
q = P[P.method == "QUBO"].groupby(["cell_type", "gene"]).size(); dtp = P[P.method == "DE_top"].groupby(["cell_type", "gene"]).size()
U2 = pd.DataFrame({"quboFS": q}).join(dtp.rename("DE_top"), how="left").fillna(0).astype(int).reset_index()
U2[(U2.quboFS >= 3)].sort_values(["cell_type", "quboFS"], ascending=[True, False]).to_csv(O / "recurrent_quboFS_vs_DEtop.csv", index=False)
# ---- figures ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = {"QUBO": "#2a78d6", "ElasticNet": "#eb6834", "LASSO": "#1baf7a", "mRMR": "#eda100", "DE_top": "#e87ba4"}
MK = {"QUBO": "o", "ElasticNet": "s", "LASSO": "^", "mRMR": "D", "DE_top": "v"}; INK2, GRID = "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2})
f, axs = plt.subplots(1, 3, figsize=(10, 3.0), sharey=True)
for ax, col, ttl in zip(axs, ["jaccard_gene", "nogueira", "jaccard_module"], ["Gene-level Jaccard", "Nogueira stability (gene level)", "Module-level Jaccard (|ρ| ≥ 0.5 modules)"]):
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for i, m in enumerate(METHODS):
        v = TA[TA.method == m][col].values; y = len(METHODS) - 1 - i
        ax.scatter(v, y + np.linspace(-0.2, 0.2, len(v)), s=9, color=COL[m], alpha=0.45, lw=0); ax.plot(v.mean(), y, MK[m], color=COL[m], ms=7, mec="white", mew=1)
    ax.set_title(ttl, fontsize=8.5); ax.set_yticks(range(len(METHODS))); ax.set_yticklabels([LABEL[m] for m in METHODS][::-1])
f.tight_layout(); f.savefig(FIG / "figD1_stability.png", dpi=200, bbox_inches="tight"); f.savefig(FIG / "figD1_stability.pdf", bbox_inches="tight"); plt.close(f)
f, ax = plt.subplots(figsize=(4.6, 3.4)); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
TB["rel_per_gene"] = TB.relevance_retained * 10 / TB.n   # per-gene relevance vs the mean of the 10 most relevant (quboFS panels may have < 10 genes)
cm = TB.groupby(["method", "cell_type"])[["rel_per_gene", "mean_abs_r"]].mean().reset_index().rename(columns={"rel_per_gene": "relevance_retained"})
for m in METHODS:
    x = cm[cm.method == m]; ax.scatter(x.mean_abs_r, x.relevance_retained, s=14, color=COL[m], alpha=0.45, lw=0, marker=MK[m])
    ax.plot(x.mean_abs_r.mean(), x.relevance_retained.mean(), MK[m], color=COL[m], ms=9, mec="white", mew=1.2)
    ax.annotate(LABEL[m], (x.mean_abs_r.mean(), x.relevance_retained.mean()), xytext=(7, 3), textcoords="offset points", fontsize=7.5)
ax.set_xlabel("mean |ρ| within panel (training donors)"); ax.set_ylabel("relevance per selected gene\n(mean |z|·C ÷ mean of top 10)")
ax.set_title("Same top-20 candidates: relevance vs redundancy\nsmall = cell type (mean of 15 folds); large = method mean", fontsize=8)
f.tight_layout(); f.savefig(FIG / "figD2_relevance_vs_redundancy.png", dpi=200, bbox_inches="tight"); f.savefig(FIG / "figD2_relevance_vs_redundancy.pdf", bbox_inches="tight"); plt.close(f)
# network example
pool = exB["pool"]; ang = np.linspace(0, 2 * np.pi, len(pool), endpoint=False); pos = {g: (np.cos(a), np.sin(a)) for g, a in zip(pool, ang)}
f, axs = plt.subplots(1, 3, figsize=(11, 3.9))
for ax, m in zip(axs, ["DE_top", "mRMR", "QUBO"]):
    sel = exB["panels"][m]; R = exB["R"]
    for a_, b_ in itertools.combinations(sel, 2):
        r = abs(R.loc[a_, b_])
        if r > 0.7: ax.plot([pos[a_][0], pos[b_][0]], [pos[a_][1], pos[b_][1]], color=COL[m], lw=0.6 + 2.5 * (r - 0.7) / 0.3, alpha=0.7, zorder=1)
    relmax = max(exB["rel"].values())
    for g in pool:
        on = g in sel; ax.scatter(*pos[g], s=30 + 170 * exB["rel"][g] / relmax, color=COL[m] if on else "#d8d7d2", edgecolor="white", lw=0.8, zorder=2)
        if on: ax.annotate(g, pos[g], xytext=(pos[g][0] * 9, pos[g][1] * 9), textcoords="offset points", fontsize=6.3, ha="left" if pos[g][0] >= 0 else "right", va="center")
    n_e = sum(abs(R.loc[a_, b_]) > 0.7 for a_, b_ in itertools.combinations(sel, 2))
    ax.set_title(f"{LABEL[m]}: {len(sel)} genes, {n_e} pairs with |r| > 0.7", fontsize=8.5); ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.4, 1.4); ax.axis("off"); ax.set_aspect("equal")
f.suptitle("B cells, Pappalardo held out, outer fold 1 (fixed in advance). Circle = the 20 candidates; size = relevance |z|·C; colour = selected; line = |r| > 0.7 on training donors",
           fontsize=7.5, color=INK2, y=0.02)
f.tight_layout(); f.savefig(FIG / "figD3_network_B_example.png", dpi=200, bbox_inches="tight"); f.savefig(FIG / "figD3_network_B_example.pdf", bbox_inches="tight"); plt.close(f)
pd.set_option("display.width", 230)
print("== A. stability (mean over 8 cell types)"); print(SA.round(3).to_string())
print("cell types where quboFS has the lowest gene-level Jaccard:", int(TA.pivot_table(index="cell_type", columns="method", values="jaccard_gene").idxmin(1).eq("QUBO").sum()), "/ 8")
print("\n== B. same candidates (mean over 15 folds x 8 cell types)"); print(SB.round(3).to_string())
print("\n== C. key genes: folds selected (of 15) and folds in pool"); print(TC.to_string(index=False))
print("\nB example panels:", {m: exB["panels"][m] for m in ["DE_top", "mRMR", "QUBO"]})
