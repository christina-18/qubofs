"""Sensitivity analysis: exact QUBO solution instead of simulated annealing. The exact variant (run_csf.py /
run_pbmc.py --solver exact) differs from the primary analysis only in the optimiser: same folds, candidate pools,
relevance, correlation matrices, (lambda, gamma) chosen in the primary analysis, seeds, classifier and weight rule;
the inner CV is re-predicted with the exact panels so that the cell-type weights follow the primary rule.
1. Selection (final models x cell types): genes shared by the SA and exact panels, objective x'Qx of each, panel size,
   relevance per gene, mean |rho|; the panel size of the exact optimum without the cardinality term and the smallest
   lambda that would make a 10-gene panel optimal (lambda*). Inner folds: share of identical panels.
2. Training |rho|: exact quboFS vs each comparator and vs SA (cell type = unit of inference).
3. Prediction: ROC-AUC, PR-AUC and Brier score (CSF: per held-out cohort and mean; PBMC: within-cohort stratified
   ROC-AUC, per-cohort ROC/PR/Brier, mean of 3 repeats) for SA, exact and the comparators, with paired differences.
Usage: python3 exact_solver.py csf|pbmc
Output: $QUBOFS_RESULTS_DIR/evaluation/exact_solver/<tissue>/"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

TISSUE = sys.argv[1]; assert TISSUE in ("csf", "pbmc")
EXV = CFG.variant(CFG.PRIMARY, "exact")
O = out(f"exact_solver/{TISSUE}")
RAW = D.Raw(CFG.CSF_RAW if TISSUE == "csf" else CFG.PBMC_RAW)
M2 = METHODS + ["QUBO_exact"]
if TISSUE == "csf":
    units = [(ho, k, f"{ho}_outer{k}", csf_dir() / "outer" / ho / f"outer_{k}", csf_dir(EXV) / "outer" / ho / f"outer_{k}")
             for ho in HOS for k in OUTERS]
    inner_dirs = {uid: (csf_dir() / "inner" / ho / f"outer_{k}", csf_dir(EXV) / "inner" / ho / f"outer_{k}") for ho, k, uid, _, _ in units}
else:
    units = [(None, u, u, pbmc_dir() / "final" / u, pbmc_dir(EXV) / "final" / u) for u in PBMC_CV]
    inner_dirs = {u: (pbmc_dir() / "inner" / u, pbmc_dir(EXV) / "inner" / u) for u in PBMC_CV}


def load_unit(d_sa, d_ex):
    res = model_outputs(d_sa)
    res["QUBO_exact"] = model_outputs(d_ex, ["QUBO"])["QUBO"]
    return res


def bundles_of(ho, k):
    return csf_outer_bundles(RAW, ho, k)[0] if TISSUE == "csf" else pbmc_unit_bundles(RAW, k)[0]


# ---- 1. selection ----
rows = []
for ho, k, uid, dsa, dex in units:
    Psa = pd.read_csv(dsa / "panels.csv", keep_default_na=False); Psa = Psa[Psa.method == "QUBO"]
    Pex = pd.read_csv(dex / "panels.csv", keep_default_na=False)
    setting = Psa.setting.iloc[0]; assert set(Pex.setting) == {setting}
    lam, gam = S.parse_setting(setting)
    for ct, b in bundles_of(ho, k).items():
        sa = Psa.gene[Psa.cell_type == ct].tolist(); ex = Pex.gene[Pex.cell_type == ct].tolist()
        pool = S.candidates(b)
        if len(pool) < 5 or not sa:
            continue
        Qm = S.qubo_matrix(b, pool, lam, gam)
        x_opt, E_opt, n_ties = S.solve_exact(Qm)
        assert {pool[i] for i in np.where(x_opt == 1)[0]} == set(ex), (uid, ct, "exact panel is not the optimum")
        xs = np.array([g in set(sa) for g in pool], float); xe = np.array([g in set(ex) for g in pool], float)
        rel = dict(zip(pool, S.relevance(b, pool))); top = np.mean(sorted(rel.values(), reverse=True)[:10])
        # objective without the cardinality term, by panel size (float32 enumeration) -> lambda* for a 10-gene optimum
        SZ = S.states(len(pool)).sum(1).astype(int); E = S.all_energies(Qm)
        E0m = pd.Series(E - lam * (SZ - CFG.K) ** 2).groupby(SZ).min()
        lam_star = max(0.0, max((E0m[CFG.K] - E0m[m]) / (m - CFG.K) ** 2 for m in E0m.index if m != CFG.K))
        rows.append(dict(unit=uid, cohort=ho if ho else "PBMC", cell_type=ct, setting=setting, n_pool=len(pool), n_sa=len(sa), n_exact=len(ex),
                         overlap=len(set(sa) & set(ex)), identical=set(sa) == set(ex), n_opt_ties=n_ties,
                         E_sa=float(xs @ Qm @ xs), E_exact=float(xe @ Qm @ xe),
                         rel_per_gene_sa=np.mean([rel[g] for g in sa]) / top, rel_per_gene_exact=np.mean([rel[g] for g in ex]) / top,
                         rho_sa=S.panel_rho(b, sa), rho_exact=S.panel_rho(b, ex), lambda_star_for_10=lam_star))
C = pd.DataFrame(rows); C.to_csv(O / "selection_final_per_unit_celltype.csv", index=False)
irows = []
for uid, (dsa, dex) in inner_dirs.items():
    a = pd.read_csv(dsa / "panels_qubo.csv", keep_default_na=False); e = pd.read_csv(dex / "panels_qubo.csv", keep_default_na=False)
    a = a[a.setting == e.setting.iloc[0]]
    for (j, ct), g in a.groupby(["inner_fold", "cell_type"]):
        ge = set(e.gene[(e.inner_fold == j) & (e.cell_type == ct)])
        irows.append(dict(unit=uid, inner_fold=j, cell_type=ct, n_sa=len(g), n_exact=len(ge), overlap=len(set(g.gene) & ge), identical=set(g.gene) == ge))
Ci = pd.DataFrame(irows); Ci.to_csv(O / "selection_inner_per_fold_celltype.csv", index=False)


def summary(x, level):
    return dict(level=level, n=len(x), identical=x.identical.mean(), overlap_genes=x.overlap.sum(), sa_genes=x.n_sa.sum(),
                exact_genes=x.n_exact.sum(), size_sa=x.n_sa.mean(), size_exact=x.n_exact.mean(),
                E_gap_mean=(x.E_sa - x.E_exact).mean() if "E_sa" in x else np.nan, E_gap_max=(x.E_sa - x.E_exact).max() if "E_sa" in x else np.nan,
                exact_le_sa=bool((x.E_exact <= x.E_sa + 1e-9).all()) if "E_sa" in x else np.nan,
                rel_sa=x.rel_per_gene_sa.mean() if "E_sa" in x else np.nan, rel_exact=x.rel_per_gene_exact.mean() if "E_sa" in x else np.nan,
                rho_sa=x.rho_sa.mean() if "E_sa" in x else np.nan, rho_exact=x.rho_exact.mean() if "E_sa" in x else np.nan)


SS = pd.DataFrame([summary(C, "final, all")] + [summary(x, f"final, {c}") for c, x in C.groupby("cohort")] +
                  [summary(x, f"final, {c}") for c, x in C.groupby("cell_type")] + [summary(Ci, "inner folds, all")])
SS.to_csv(O / "selection_summary.csv", index=False)
SZ_ = C.groupby("n_exact").size().rename("n_problems"); SZ_.to_csv(O / "exact_panel_sizes.csv")

# ---- 2. training |rho| ----
rr = []
for ho, k, uid, dsa, dex in units:
    R = pd.read_csv(dsa / "rho.csv", keep_default_na=False); R["unit"] = uid; rr.append(R[["cell_type", "unit", "method", "rho"]])
    c = C[C.unit == uid]; rr.append(pd.DataFrame(dict(cell_type=c.cell_type, unit=uid, method="QUBO_exact", rho=c.rho_exact)))
W = pd.concat(rr).pivot_table(index=["cell_type", "unit"], columns="method", values="rho")
rng_rho = np.random.default_rng(20261021)
RT = pd.concat([celltype_paired_test(W, "QUBO_exact", COMPARATORS + ["QUBO"], rng_rho), celltype_paired_test(W, "QUBO", COMPARATORS, rng_rho)])
RT.to_csv(O / "rho_paired_tests.csv", index=False)

# ---- 3. prediction ----
PT, PD_, CT, WT, DP = [], [], [], [], []
if TISSUE == "csf":
    rng = np.random.default_rng(SEED); U = {(ho, k): load_unit(dsa, dex) for ho, k, uid, dsa, dex in units}
    for scheme in ["weighted", "equal"]:
        if scheme == "equal": rng = np.random.default_rng(SEED)
        point, boot = {}, {}
        for ho in HOS:
            dx = U[(ho, 1)]["QUBO"][0].groupby("donor").diagnosis.first(); ms, hd = sorted(dx.index[dx == "MS"]), sorted(dx.index[dx == "HD"])
            y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]
            cm = rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B); ch = rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)
            for m in M2:
                v = {x: [] for x in ["roc", "pr", "brier"]}; vb = {x: [] for x in v}
                for k in OUTERS:
                    h, I = U[(ho, k)][m]; p, w = ensemble(h, I, scheme)
                    if scheme == "weighted": WT += [dict(cohort=ho, outer=k, method=m, cell_type=c, weight=x) for c, x in w.items()]
                    DP += [dict(scheme=scheme, cohort=ho, outer=k, method=m, donor=dn, diagnosis=dx[dn], prob_MS=float(p[dn])) for dn in p.index]
                    pm, ph = p.loc[ms].values, p.loc[hd].values; sc = np.r_[pm, ph]
                    v["roc"].append(auc_w(pm, ph, np.ones(len(ms)), np.ones(len(hd)))[0]); vb["roc"].append(auc_w(pm, ph, cm, ch))
                    v["pr"].append(ap_w(sc, y, np.ones(len(sc)))[0]); vb["pr"].append(ap_w(sc, y, np.c_[cm, ch]))
                    v["brier"].append(brier_w(sc, y, np.ones(len(sc)))[0]); vb["brier"].append(brier_w(sc, y, np.c_[cm, ch]))
                    if scheme == "weighted" and m in ("QUBO", "QUBO_exact"):
                        for c, g in h.groupby("cell_type"):
                            gg = g.set_index("donor").prob_MS; a, b_ = gg.reindex(ms).dropna().values, gg.reindex(hd).dropna().values
                            if len(a) and len(b_): CT.append(dict(cohort=ho, outer=k, method=m, cell_type=c, roc_auc=auc_w(a, b_, np.ones(len(a)), np.ones(len(b_)))[0]))
                for x in v: point[(ho, m, x)] = np.mean(v[x]); boot[(ho, m, x)] = np.mean(vb[x], 0)
        for m in M2:
            for x in ["roc", "pr", "brier"]:
                point[("Mean", m, x)] = np.mean([point[(h, m, x)] for h in HOS]); boot[("Mean", m, x)] = np.mean([boot[(h, m, x)] for h in HOS], 0)
        cohorts = HOS + ["Mean"]; metrics = ["roc", "pr", "brier"]
        PT += [dict(scheme=scheme, cohort=c, method=m, metric=x, value=point[(c, m, x)], lo=pc(boot[(c, m, x)])[0], hi=pc(boot[(c, m, x)])[1]) for c in cohorts for m in M2 for x in metrics]
        for c in cohorts:
            for x in metrics:
                for a_, b_ in [("QUBO_exact", "QUBO")] + [(q, m) for q in ("QUBO_exact", "QUBO") for m in COMPARATORS]:
                    d = boot[(c, a_, x)] - boot[(c, b_, x)]
                    PD_.append(dict(scheme=scheme, cohort=c, metric=x, model=a_, comparator=b_, diff=point[(c, a_, x)] - point[(c, b_, x)], lo=pc(d)[0], hi=pc(d)[1],
                                    lo_bonf=pc(d, 0.05 / 4)[0], hi_bonf=pc(d, 0.05 / 4)[1], frac_gt0=(d > 0).mean()))
else:
    COHN = CFG.PBMC_COHORTS
    INFO = RAW.info
    donors = {c: (sorted(INFO.index[(INFO.cohort.map(COHN) == c) & (INFO.diagnosis == "MS")]), sorted(INFO.index[(INFO.cohort.map(COHN) == c) & (INFO.diagnosis == "HD")])) for c in COHN.values()}
    U = {uid: load_unit(dsa, dex) for _, _, uid, dsa, dex in units}; rng = np.random.default_rng(20261020)
    for scheme in ["weighted", "equal"]:
        draws = {c: (rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B), rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)) for c, (ms, hd) in donors.items()}
        point, boot = {}, {}
        for m in M2:
            pr_ = {}
            for uid in U:
                h, I = U[uid][m]; p, w = ensemble(h, I, scheme)
                if scheme == "weighted": WT += [dict(cohort="PBMC", outer=uid, method=m, cell_type=c, weight=x) for c, x in w.items()]
                for dn, x in p.items(): pr_.setdefault(int(uid[4]), {})[dn] = x
                if scheme == "weighted" and m in ("QUBO", "QUBO_exact"):
                    for c, g in h.groupby("cell_type"):
                        for coh, (ms, hd) in donors.items():
                            gg = g.set_index("donor").prob_MS; a, b_ = gg.reindex(ms).dropna().values, gg.reindex(hd).dropna().values
                            if len(a) and len(b_): CT.append(dict(cohort=coh, outer=uid, method=m, cell_type=c, roc_auc=auc_w(a, b_, np.ones(len(a)), np.ones(len(b_)))[0]))
            keys = ["within"] + [f"{x}_{c}" for x in ("roc", "pr", "brier") for c in COHN.values()] + ["brier_all"]
            v = {x: [] for x in keys}; vb = {x: [] for x in keys}
            for r in (1, 2, 3):
                p = pd.Series(pr_[r]); num_p = den_p = 0; num_b = np.zeros(B); den_b = np.zeros(B); sq_p = []; sq_b = []; w_b = []
                for c, (ms, hd) in donors.items():
                    assert set(ms + hd) <= set(p.index), "missing donor prediction"
                    DP += [dict(scheme=scheme, cohort=c, outer=f"rep{r}", method=m, donor=dn, diagnosis="MS" if dn in ms else "HD", prob_MS=float(p[dn])) for dn in ms + hd]
                    pm, ph = p.loc[ms].values, p.loc[hd].values; cm, chh = draws[c]; y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]; sc = np.r_[pm, ph]
                    Mx = (pm[:, None] > ph[None, :]) + 0.5 * (pm[:, None] == ph[None, :])
                    num_p += Mx.sum(); den_p += Mx.size; num_b += np.einsum("bi,ij,bj->b", cm, Mx, chh); den_b += cm.sum(1) * chh.sum(1)
                    v[f"roc_{c}"].append(Mx.mean()); vb[f"roc_{c}"].append(auc_w(pm, ph, cm, chh))
                    v[f"pr_{c}"].append(ap_w(sc, y, np.ones(len(sc)))[0]); vb[f"pr_{c}"].append(ap_w(sc, y, np.c_[cm, chh]))
                    v[f"brier_{c}"].append(brier_w(sc, y, np.ones(len(sc)))[0]); vb[f"brier_{c}"].append(brier_w(sc, y, np.c_[cm, chh]))
                    sq_p.append((sc - y) ** 2); w_b.append(np.c_[cm, chh])
                v["within"].append(num_p / den_p); vb["within"].append(num_b / den_b)
                sq = np.concatenate(sq_p); wb = np.concatenate(w_b, 1); v["brier_all"].append(sq.mean()); vb["brier_all"].append((wb * sq).sum(1) / wb.sum(1))
            for x in keys: point[(m, x)] = np.mean(v[x]); boot[(m, x)] = np.mean(vb[x], 0); point[(m, x, "reps")] = ",".join(f"{q:.3f}" for q in v[x])
        PT += [dict(scheme=scheme, cohort="PBMC", method=m, metric=x, value=point[(m, x)], lo=pc(boot[(m, x)])[0], hi=pc(boot[(m, x)])[1], rep_values=point[(m, x, "reps")]) for m in M2 for x in keys]
        for x in keys:
            for a_, b_ in [("QUBO_exact", "QUBO")] + [(q, m) for q in ("QUBO_exact", "QUBO") for m in COMPARATORS]:
                d = boot[(a_, x)] - boot[(b_, x)]
                PD_.append(dict(scheme=scheme, cohort="PBMC", metric=x, model=a_, comparator=b_, diff=point[(a_, x)] - point[(b_, x)], lo=pc(d)[0], hi=pc(d)[1],
                                lo_bonf=pc(d, 0.05 / 4)[0], hi_bonf=pc(d, 0.05 / 4)[1], frac_gt0=(d > 0).mean()))
T = pd.DataFrame(PT); Dd = pd.DataFrame(PD_); CTd = pd.DataFrame(CT)
T.to_csv(O / "prediction_metrics.csv", index=False); Dd.to_csv(O / "prediction_paired_diffs.csv", index=False)
CTd.to_csv(O / "per_celltype_heldout_auc_per_unit.csv", index=False); pd.DataFrame(WT).to_csv(O / "ensemble_weights.csv", index=False)
pd.DataFrame(DP).to_csv(O / "donor_predictions.csv", index=False)
CTs = CTd.groupby(["cohort", "cell_type", "method"]).roc_auc.mean().unstack("method"); CTs["diff"] = CTs["QUBO_exact"] - CTs["QUBO"]; CTs.to_csv(O / "per_celltype_heldout_auc_summary.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("== selection (final models)"); print(SS[~SS.level.str.replace("final, ", "").isin(CFG.CELL_TYPES)].round(3).to_string(index=False))
print("\n== training |rho| tests"); print(RT.round(4).to_string(index=False))
print("\n== prediction (weighted)"); print(T[T.scheme == "weighted"].pivot_table(index=["method"], columns=["cohort", "metric"], values="value").round(3).to_string())
k = Dd[(Dd.scheme == "weighted") & ((Dd.comparator == "QUBO") | Dd.comparator.isin(["ElasticNet", "mRMR", "LASSO", "DE_top"]))]
print(k[k.cohort.isin(["Mean", "PBMC"])].round(3).to_string(index=False))
print("\n== exact - SA by cohort (weighted)"); print(Dd[(Dd.scheme == "weighted") & (Dd.comparator == "QUBO")].round(3).to_string(index=False))
print("\n== equal weights, exact - SA"); print(Dd[(Dd.scheme == "equal") & (Dd.comparator == "QUBO")].round(3).to_string(index=False))
print("\n== exact panel sizes:", SZ_.to_dict(), "| lambda* (exact optimum < 10 genes): median", round(C[C.n_exact < 10].lambda_star_for_10.median(), 2))
print("\n== per-cell-type held-out AUC (mean over units/cohorts)"); print(CTs.groupby(level="cell_type")[["QUBO", "QUBO_exact", "diff"]].mean().round(3).to_string())
