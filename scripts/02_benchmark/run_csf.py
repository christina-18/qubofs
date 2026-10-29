"""CSF benchmark: nested cross-validation inside each leave-one-cohort-out setting, all five selection methods.

For one held-out cohort (Pappalardo, Heming or Ramesh) and one outer fold k:
  inner : python3 run_csf.py inner --holdout Ramesh --outer 3 --group qubo|base [--filter F] [--solver S]
          Inner 5-fold CV on the outer-training donors with inner-training edgeR statistics: every method x setting
          selects genes on the inner-training donors and predicts the inner-validation donors per cell type.
  outer : python3 run_csf.py outer --holdout Ramesh --outer 3 [--filter F] [--solver S]
          Chooses each method's setting from the inner results (within-cohort stratified AUC of the out-of-fold
          predictions; ties within 0.01 -> lower inner |rho| -> smaller setting), re-selects on the outer-training
          donors and predicts the held-out cohort per cell type. The held-out cohort is used only for prediction.
--solver exact (sensitivity) runs quboFS only, at the setting chosen in the primary (simulated annealing) analysis.
--lasso-max-iter 20000 / --lasso-pad zc (sensitivity) run LASSO only, with more coordinate-descent iterations or with
the completion of short LASSO panels by |z|*C instead of |z| (inner CV included).
Output: $QUBOFS_RESULTS_DIR/csf/<variant>/{inner,outer}/<holdout>/outer_<k>/"""
import argparse
import re

import numpy as np
import pandas as pd

import config as CFG
import data as D
import selection as S

RAW = None


def folds_of(ho, k):
    f = pd.read_csv(CFG.CSF_OUTER / ho / "folds_assignment.csv")
    tr = f.donor[(f.set == "cv") & (f.fold.astype(str) != str(k))].tolist()
    va = f.donor[(f.set == "cv") & (f.fold.astype(str) == str(k))].tolist()
    return tr, va, f.donor[f.set == "heldout"].tolist()


def out_dir(a, stage):
    return CFG.RESULTS / "csf" / CFG.variant(a.filter, a.solver, a.lasso_max_iter, a.lasso_pad) / stage / a.holdout / f"outer_{a.outer}"


def primary_setting(a, method="QUBO"):
    g = pd.read_csv(CFG.RESULTS / "csf" / a.filter / "outer" / a.holdout / f"outer_{a.outer}" / "inner_grid.csv", keep_default_na=False)
    return g[(g.method == method) & (g.chosen.astype(str) == "True")].setting.iloc[0]


def lasso_only(a):
    return a.lasso_max_iter != CFG.LASSO_MAX_ITER or a.lasso_pad != CFG.LASSO_PAD


def group_settings(a):
    if a.solver == "exact":
        return [("QUBO", primary_setting(a))] if a.group == "qubo" else []
    if lasso_only(a):
        return [(m, s) for m, s in CFG.GROUPS["base"] if m == "LASSO"] if a.group == "base" else []
    return CFG.GROUPS[a.group]


def inner(a):
    tr, va, held = folds_of(a.holdout, a.outer)
    od = CFG.CSF_INNER / a.holdout / f"outer_{a.outer}" / "rep_1"
    sp = pd.read_csv(od / "inner_folds.csv")
    assert set(sp.donor) == set(tr) and not (set(sp.donor) & (set(held) | set(va))), "inner split leak"
    out = out_dir(a, "inner"); out.mkdir(parents=True, exist_ok=True)
    pct, pan, rho, univ = [], [], [], []
    settings = group_settings(a)
    for j in sorted(sp.inner_fold.unique()):
        tr_d = sp.donor[sp.inner_fold != j].tolist(); va_d = sp.donor[sp.inner_fold == j].tolist()
        dirs = {ct: od / ct / f"inner_{j}" for ct in CFG.CELL_TYPES}
        for ct, d in dirs.items():
            if (d / "train_donors.csv").exists():
                used = set(pd.read_csv(d / "train_donors.csv").donor)
                assert not (used & set(va_d)) and used <= set(tr_d), f"{ct} inner {j}: edgeR donor leak"
        B = D.make_bundles(RAW, D.read_tstats(dirs), tr_d, va_d, j, a.filter)
        for ct, b in B.items():
            assert set(b["train"]["donors"]) <= set(tr_d) and (b["val"] is None or set(b["val"]["donors"]) <= set(va_d))
            univ.append(dict(inner_fold=j, cell_type=ct, n_edger=len(b["train"]["genes"]), n_allowed=len(b["allowed"])))
        for m, s in settings:
            sel = {}
            for ct, b in B.items():
                g = S.select(m, s, b, CFG.SEED + a.outer * 100 + 10 * j + S.ct_seed(ct), a.solver, a.lasso_max_iter, a.lasso_pad)
                if g:
                    sel[ct] = g
                    pan += [dict(method=m, setting=s, inner_fold=j, cell_type=ct, gene=x) for x in g]
                    rho.append(dict(method=m, setting=s, inner_fold=j, cell_type=ct, rho=S.panel_rho(b, g)))
            for r in S.predict_ct(B, sel, "val"):
                pct.append(dict(method=m, setting=s, inner_fold=j, cohort=RAW.info.loc[r["donor"], "cohort"], **r))
        print(f"  {a.holdout} outer {a.outer} inner fold {j}: train {len(tr_d)} / val {len(va_d)}", flush=True)
    for nm, rows in [("pred_ct", pct), ("panels", pan), ("rho", rho), ("universe", univ)]:
        pd.DataFrame(rows).to_csv(out / f"{nm}_{a.group}.csv", index=False)


def setting_key(s):
    return tuple(float(x) for x in re.findall(r"[0-9]*\.?[0-9]+", s)) if s != "NA" else (0.0,)


def choose_settings(inner_dir, groups):
    """Chosen setting per method from the inner out-of-fold predictions, plus the inner per-cell-type predictions."""
    P = pd.concat([pd.read_csv(inner_dir / f"pred_ct_{g}.csv", keep_default_na=False) for g in groups])
    R = pd.concat([pd.read_csv(inner_dir / f"rho_{g}.csv", keep_default_na=False) for g in groups])
    rows = []
    for (m, s), g in P.groupby(["method", "setting"]):
        e = g.groupby("donor").agg(prob_MS=("prob_MS", "mean"), diagnosis=("diagnosis", "first"), cohort=("cohort", "first"))
        rows.append(dict(method=m, setting=s, within_auc=D.within_auc(e), inner_rho=np.nanmean(R[(R.method == m) & (R.setting == s)].rho), n_oof=len(e)))
    G = pd.DataFrame(rows); ch = {}; out = []
    for m, g in G.groupby("method"):
        t = g[g.within_auc >= np.nanmax(g.within_auc) - CFG.TIE_TOL].copy() if g.within_auc.notna().any() else g.copy()
        t["_k"] = t.setting.map(setting_key); t = t.sort_values(["inner_rho", "_k"])
        ch[m] = t.setting.iloc[0]
        out += [dict(r, chosen=(r["setting"] == ch[m]), n_tied=len(t)) for r in g.to_dict("records")]
    return ch, P, pd.DataFrame(out)


def outer(a):
    tr, va, held = folds_of(a.holdout, a.outer)
    groups = ["qubo"] if a.solver == "exact" else ["base"] if lasso_only(a) else list(CFG.GROUPS)
    ch, P, grid = choose_settings(out_dir(a, "inner"), groups)
    if a.solver == "exact":
        assert ch["QUBO"] == primary_setting(a)
    fdir = CFG.CSF_OUTER / a.holdout
    dirs = {ct: fdir / ct / "CSF" / f"fold_{a.outer}" for ct in CFG.CELL_TYPES}
    B = D.make_bundles(RAW, D.read_tstats(dirs), tr, va, a.outer, a.filter)
    for ct, b in B.items():
        b["heldout"] = RAW.block(ct, b["train"]["genes"], held, "heldout")
        assert set(b["train"]["donors"]) <= set(tr) and not (set(b["train"]["donors"]) & set(held))
    rows, pan, rho = [], [], []
    methods = ["QUBO"] if a.solver == "exact" else ["LASSO"] if lasso_only(a) else CFG.METHODS
    for m in methods:
        s = ch[m]; sel = {}
        for ct, b in B.items():
            g = S.select(m, s, b, CFG.SEED + a.outer * 100 + S.ct_seed(ct), a.solver, a.lasso_max_iter, a.lasso_pad)
            if g:
                sel[ct] = g
                pan += [dict(method=m, setting=s, cell_type=ct, gene=x) for x in g]
                rho.append(dict(method=m, setting=s, cell_type=ct, rho=S.panel_rho(b, g)))
        rows += [dict(method=m, setting=s, **r) for r in S.predict_ct(B, sel, "heldout", require_val=True)]
    out = out_dir(a, "outer"); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "held_ct.csv", index=False)
    pd.DataFrame(pan).to_csv(out / "panels.csv", index=False)
    pd.DataFrame(rho).to_csv(out / "rho.csv", index=False)
    grid.to_csv(out / "inner_grid.csv", index=False)
    pd.concat([P[(P.method == m) & (P.setting == ch[m])] for m in methods]).to_csv(out / "inner_pred_ct_chosen.csv", index=False)
    pd.DataFrame([dict(cell_type=ct, n_edger=len(b["train"]["genes"]), n_allowed=len(b["allowed"])) for ct, b in B.items()]
                 ).to_csv(out / "universe.csv", index=False)
    print(CFG.variant(a.filter, a.solver, a.lasso_max_iter, a.lasso_pad), a.holdout, a.outer, " | ".join(f"{m}:{ch[m]}" for m in methods), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["inner", "outer"])
    ap.add_argument("--holdout", required=True, choices=list(CFG.HOLDOUTS))
    ap.add_argument("--outer", type=int, required=True, choices=CFG.OUTER_FOLDS)
    ap.add_argument("--group", choices=list(CFG.GROUPS), default="qubo")
    ap.add_argument("--filter", choices=CFG.FILTERS, default=CFG.PRIMARY)
    ap.add_argument("--solver", choices=CFG.SOLVERS, default="sa")
    ap.add_argument("--lasso-max-iter", type=int, default=CFG.LASSO_MAX_ITER,
                    help="LASSO coordinate-descent iterations; a value other than the default runs LASSO only (sensitivity)")
    ap.add_argument("--lasso-pad", choices=["z", "zc"], default=CFG.LASSO_PAD,
                    help="LASSO completion when fewer than K coefficients are non-zero: |z| (default) or |z|*C (sensitivity; LASSO only)")
    a = ap.parse_args()
    RAW = D.Raw(CFG.CSF_RAW)
    inner(a) if a.stage == "inner" else outer(a)
