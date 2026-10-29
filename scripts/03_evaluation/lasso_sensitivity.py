"""Sensitivity analyses of the LASSO comparator. Each variant re-runs the whole nested CV for LASSO (inner CV, choice
of C, cell-type weights, final models); everything else is the primary analysis.
  iter20000 : 20,000 instead of 120 coordinate-descent iterations   (run_csf.py / run_pbmc.py --lasso-max-iter 20000)
  padzc     : panels with fewer than 10 non-zero coefficients completed by |z|*C instead of |z|   (--lasso-pad zc)
1. Selection (final models x cell types): chosen C, identical panels, training |rho|, number of non-zero coefficients;
   for iter20000 also the convergence of the coefficients (change from 20,000 to 40,000 iterations; relative change
   of the penalised objective).
2. Training |rho|: quboFS vs the LASSO variant, LASSO variant vs primary LASSO (cell type = unit of inference).
3. Prediction: as exact_solver.py, for LASSO (primary), LASSO_variant and quboFS.
Usage: python3 lasso_sensitivity.py csf|pbmc iter20000|padzc
Output: $QUBOFS_RESULTS_DIR/evaluation/lasso_sensitivity/<variant>/<tissue>/"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403
from qubo_core import LogRegL1  # noqa: E402

TISSUE, NAME = sys.argv[1], sys.argv[2]; assert TISSUE in ("csf", "pbmc") and NAME in ("iter20000", "padzc")
LV = CFG.variant(CFG.PRIMARY, "sa", 20000) if NAME == "iter20000" else CFG.variant(CFG.PRIMARY, "sa", lasso_pad="zc")
O = out(f"lasso_sensitivity/{NAME}/{TISSUE}")
RAW = D.Raw(CFG.CSF_RAW if TISSUE == "csf" else CFG.PBMC_RAW)
M2 = METHODS + ["LASSO_variant"]
PAIRS = [("LASSO_variant", "LASSO"), ("QUBO", "LASSO_variant"), ("QUBO", "LASSO")]
if TISSUE == "csf":
    units = [(ho, k, f"{ho}_outer{k}", csf_dir() / "outer" / ho / f"outer_{k}", csf_dir(LV) / "outer" / ho / f"outer_{k}")
             for ho in HOS for k in OUTERS]
else:
    units = [(None, u, u, pbmc_dir() / "final" / u, pbmc_dir(LV) / "final" / u) for u in PBMC_CV]


def load_unit(d_pr, d_lv):
    res = model_outputs(d_pr)
    res["LASSO_variant"] = model_outputs(d_lv, ["LASSO"])["LASSO"]
    return res


def objective(clf, X, y, C):
    p = 1 / (1 + np.exp(-np.clip(X @ clf.coef_ + clf.intercept_, -30, 30)))
    return -np.sum(y * np.log(p + 1e-300) + (1 - y) * np.log(1 - p + 1e-300)) + np.abs(clf.coef_).sum() / C


# ---- 1. selection and convergence ----
rows = []
for ho, k, uid, dpr, dlv in units:
    Bd = csf_outer_bundles(RAW, ho, k)[0] if TISSUE == "csf" else pbmc_unit_bundles(RAW, k)[0]
    Ppr = pd.read_csv(dpr / "panels.csv", keep_default_na=False); Ppr = Ppr[Ppr.method == "LASSO"]
    Plv = pd.read_csv(dlv / "panels.csv", keep_default_na=False)
    for ct, b in Bd.items():
        a, c = Ppr[Ppr.cell_type == ct], Plv[Plv.cell_type == ct]
        if a.empty or c.empty:
            continue
        r = dict(unit=uid, cohort=ho if ho else "PBMC", cell_type=ct, C_primary=a.setting.iloc[0], C_variant=c.setting.iloc[0],
                 identical=set(a.gene) == set(c.gene), overlap=len(set(a.gene) & set(c.gene)), n_primary=len(a), n_variant=len(c),
                 rho_primary=S.panel_rho(b, list(a.gene)), rho_variant=S.panel_rho(b, list(c.gene)))
        tr = S.build_X_y(b, S.candidates(b), "train"); X = tr["X"]; mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1; Xz = (X - mu) / sd
        for lab, setting in [("primary", r["C_primary"]), ("variant", r["C_variant"])]:
            Cv = float(setting[1:])
            if NAME == "padzc":
                r[f"nonzero_{lab}"] = int((np.abs(LogRegL1(C=Cv, max_iter=120).fit(Xz, tr["y"]).coef_) > 1e-6).sum())
                continue
            f120 = LogRegL1(C=Cv, max_iter=120).fit(Xz, tr["y"]); f20k = LogRegL1(C=Cv, max_iter=20000).fit(Xz, tr["y"])
            f40k = LogRegL1(C=Cv, max_iter=40000, tol=1e-12).fit(Xz, tr["y"])
            o40 = objective(f40k, Xz, tr["y"], Cv)
            r[f"nonzero_{lab}"] = int((np.abs((f120 if lab == "primary" else f20k).coef_) > 1e-6).sum())
            r[f"coefdiff_120_vs_40k_{lab}"] = float(np.abs(f120.coef_ - f40k.coef_).max())
            r[f"coefdiff_20k_vs_40k_{lab}"] = float(np.abs(f20k.coef_ - f40k.coef_).max())
            r[f"objrel_120_{lab}"] = float((objective(f120, Xz, tr["y"], Cv) - o40) / abs(o40))
            r[f"objrel_20k_{lab}"] = float((objective(f20k, Xz, tr["y"], Cv) - o40) / abs(o40))
        rows.append(r)
C = pd.DataFrame(rows); C.to_csv(O / "selection_final_per_unit_celltype.csv", index=False)
SS = dict(n=len(C), identical=C.identical.mean(), same_C=(C.C_primary == C.C_variant).mean(),
          padded_primary=int((C.nonzero_primary < CFG.K).sum()), padded_variant=int((C.nonzero_variant < CFG.K).sum()),
          rho_primary=C.rho_primary.mean(), rho_variant=C.rho_variant.mean())
if NAME == "iter20000":
    SS.update(max_coefdiff_120=C.coefdiff_120_vs_40k_primary.max(), max_coefdiff_20k=C.coefdiff_20k_vs_40k_variant.max(),
              max_objrel_120=C.objrel_120_primary.max(), max_objrel_20k=C.objrel_20k_variant.max())
SS = pd.Series(SS)
SS.to_csv(O / "selection_summary.csv")

# ---- 2. training |rho| ----
rr = []
for ho, k, uid, dpr, dlv in units:
    R = pd.read_csv(dpr / "rho.csv", keep_default_na=False); R["unit"] = uid; rr.append(R[["cell_type", "unit", "method", "rho"]])
    R = pd.read_csv(dlv / "rho.csv", keep_default_na=False); rr.append(pd.DataFrame(dict(cell_type=R.cell_type, unit=uid, method="LASSO_variant", rho=R.rho)))
W = pd.concat(rr).pivot_table(index=["cell_type", "unit"], columns="method", values="rho")
rng_rho = np.random.default_rng(20261023)
RT = pd.concat([celltype_paired_test(W, "QUBO", ["LASSO_variant", "LASSO"], rng_rho), celltype_paired_test(W, "LASSO_variant", ["LASSO"], rng_rho)])
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
                    if scheme == "weighted" and m in ("LASSO", "LASSO_variant"):
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
                for a_, b_ in PAIRS:
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
                if scheme == "weighted" and m in ("LASSO", "LASSO_variant"):
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
            for a_, b_ in PAIRS:
                d = boot[(a_, x)] - boot[(b_, x)]
                PD_.append(dict(scheme=scheme, cohort="PBMC", metric=x, model=a_, comparator=b_, diff=point[(a_, x)] - point[(b_, x)], lo=pc(d)[0], hi=pc(d)[1],
                                lo_bonf=pc(d, 0.05 / 4)[0], hi_bonf=pc(d, 0.05 / 4)[1], frac_gt0=(d > 0).mean()))
T = pd.DataFrame(PT); Dd = pd.DataFrame(PD_); CTd = pd.DataFrame(CT)
T.to_csv(O / "prediction_metrics.csv", index=False); Dd.to_csv(O / "prediction_paired_diffs.csv", index=False)
CTd.to_csv(O / "per_celltype_heldout_auc_per_unit.csv", index=False); pd.DataFrame(WT).to_csv(O / "ensemble_weights.csv", index=False)
pd.DataFrame(DP).to_csv(O / "donor_predictions.csv", index=False)
CTs = CTd.groupby(["cohort", "cell_type", "method"]).roc_auc.mean().unstack("method"); CTs["diff"] = CTs["LASSO_variant"] - CTs["LASSO"]; CTs.to_csv(O / "per_celltype_heldout_auc_summary.csv")

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("== selection / convergence (final models)"); print(SS.round(6).to_string())
print(C[~C.identical][["unit", "cell_type", "C_primary", "C_variant", "overlap", "nonzero_primary", "nonzero_variant", "rho_primary", "rho_variant"]].round(3).to_string(index=False))
print("\n== training |rho| tests"); print(RT.round(4).to_string(index=False))
print("\n== prediction (weighted)"); print(T[T.scheme == "weighted"].pivot_table(index=["method"], columns=["cohort", "metric"], values="value").round(3).to_string())
print(Dd[(Dd.scheme == "weighted")].round(3).to_string(index=False))
