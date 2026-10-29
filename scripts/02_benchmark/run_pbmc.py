"""PBMC benchmark (two cohorts with MS and control donors): the CSF procedure applied to the PBMC evaluation units.

Units ($QUBOFS_DATA_DIR/pbmc_units/<unit>/): cv_r<r>_t<t> (nested CV: test fold t of repeat r, 5 folds x 3 repeats;
primary) and loco_<cohort> (train on the other cohort; reference).
  inner : python3 run_pbmc.py inner --unit cv_r1_t1 --group qubo|base [--filter F] [--solver S]
  final : python3 run_pbmc.py final --unit cv_r1_t1 [--filter F] [--solver S]
          chooses the settings with the CSF rule, re-selects on all training donors and predicts the test donors.
Output: $QUBOFS_RESULTS_DIR/pbmc/<variant>/{inner,final}/<unit>/"""
import argparse

import pandas as pd

import config as CFG
import data as D
import selection as S
from run_csf import choose_settings, lasso_only

RAW = None
UNITS = CFG.PBMC_UNITS_DIR


def unit_info(u):
    sp = pd.read_csv(UNITS / u / "split.csv")
    return sp.donor[sp.role == "train"].tolist(), sp.donor[sp.role == "test"].tolist(), pd.read_csv(UNITS / u / "inner_folds.csv")


def useed(u):
    return CFG.SEED + 1000 * (CFG.PBMC_UNITS.index(u) + 1)


def out_dir(a, stage):
    return CFG.RESULTS / "pbmc" / CFG.variant(a.filter, a.solver, a.lasso_max_iter, a.lasso_pad) / stage / a.unit


def primary_setting(a, method="QUBO"):
    g = pd.read_csv(CFG.RESULTS / "pbmc" / a.filter / "final" / a.unit / "inner_grid.csv", keep_default_na=False)
    return g[(g.method == method) & (g.chosen.astype(str) == "True")].setting.iloc[0]


def group_settings(a):
    if a.solver == "exact":
        return [("QUBO", primary_setting(a))] if a.group == "qubo" else []
    if lasso_only(a):
        return [(m, s) for m, s in CFG.GROUPS["base"] if m == "LASSO"] if a.group == "base" else []
    return CFG.GROUPS[a.group]


def check_donors(d, train, exclude):
    if (d / "train_donors.csv").exists():
        used = set(pd.read_csv(d / "train_donors.csv").donor)
        assert used <= set(train) and not (used & set(exclude)), f"{d}: edgeR donor leak"


def inner(a):
    tr, te, inn = unit_info(a.unit)
    assert set(inn.donor) == set(tr) and not (set(inn.donor) & set(te)), "inner split leak"
    out = out_dir(a, "inner"); out.mkdir(parents=True, exist_ok=True)
    settings = group_settings(a)
    pct, pan, rho = [], [], []
    for j in sorted(inn.inner_fold.unique()):
        tr_d = inn.donor[inn.inner_fold != j].tolist(); va_d = inn.donor[inn.inner_fold == j].tolist()
        dirs = {ct: UNITS / a.unit / ct / f"inner_{j}" for ct in CFG.CELL_TYPES}
        for d in dirs.values():
            check_donors(d, tr_d, va_d + te)
        B = D.make_bundles(RAW, D.read_tstats(dirs), tr_d, va_d, j, a.filter)
        for m, s in settings:
            sel = {}
            for ct, b in B.items():
                g = S.select(m, s, b, useed(a.unit) + 10 * j + S.ct_seed(ct), a.solver, a.lasso_max_iter, a.lasso_pad)
                if g:
                    sel[ct] = g
                    pan += [dict(method=m, setting=s, inner_fold=j, cell_type=ct, gene=x) for x in g]
                    rho.append(dict(method=m, setting=s, inner_fold=j, cell_type=ct, rho=S.panel_rho(b, g)))
            for r in S.predict_ct(B, sel, "val"):
                pct.append(dict(method=m, setting=s, inner_fold=j, cohort=RAW.info.loc[r["donor"], "cohort"], **r))
        print(f"  {a.unit} inner {j}: train {len(tr_d)} / val {len(va_d)}, cell types {len(B)}", flush=True)
    for nm, rows in [("pred_ct", pct), ("panels", pan), ("rho", rho)]:
        pd.DataFrame(rows).to_csv(out / f"{nm}_{a.group}.csv", index=False)


def final(a):
    tr, te, _ = unit_info(a.unit)
    groups = ["qubo"] if a.solver == "exact" else ["base"] if lasso_only(a) else list(CFG.GROUPS)
    ch, P, grid = choose_settings(out_dir(a, "inner"), groups)
    if a.solver == "exact":
        assert ch["QUBO"] == primary_setting(a)
    dirs = {ct: UNITS / a.unit / "final" / ct for ct in CFG.CELL_TYPES}
    for d in dirs.values():
        check_donors(d, tr, te)
    B = D.make_bundles(RAW, D.read_tstats(dirs), tr, [], 0, a.filter)
    for ct, b in B.items():
        b["heldout"] = RAW.block(ct, b["train"]["genes"], te, "heldout")
        assert set(b["train"]["donors"]) <= set(tr) and not (set(b["train"]["donors"]) & set(te))
    rows, pan, rho = [], [], []
    methods = ["QUBO"] if a.solver == "exact" else ["LASSO"] if lasso_only(a) else CFG.METHODS
    for m in methods:
        s = ch[m]; sel = {}
        for ct, b in B.items():
            g = S.select(m, s, b, useed(a.unit) + S.ct_seed(ct), a.solver, a.lasso_max_iter, a.lasso_pad)
            if g:
                sel[ct] = g
                pan += [dict(method=m, setting=s, cell_type=ct, gene=x) for x in g]
                rho.append(dict(method=m, setting=s, cell_type=ct, rho=S.panel_rho(b, g)))
        rows += [dict(method=m, setting=s, **r) for r in S.predict_ct(B, sel, "heldout")]
    out = out_dir(a, "final"); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "held_ct.csv", index=False)
    pd.DataFrame(pan).to_csv(out / "panels.csv", index=False)
    pd.DataFrame(rho).to_csv(out / "rho.csv", index=False)
    grid.to_csv(out / "inner_grid.csv", index=False)
    pd.concat([P[(P.method == m) & (P.setting == ch[m])] for m in methods]).to_csv(out / "inner_pred_ct_chosen.csv", index=False)
    print(CFG.variant(a.filter, a.solver, a.lasso_max_iter, a.lasso_pad), a.unit, " | ".join(f"{m}:{ch[m]}" for m in methods), f"| test donors {len(te)}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["inner", "final"])
    ap.add_argument("--unit", required=True, choices=CFG.PBMC_UNITS)
    ap.add_argument("--group", choices=list(CFG.GROUPS), default="qubo")
    ap.add_argument("--filter", choices=CFG.FILTERS, default=CFG.PRIMARY)
    ap.add_argument("--solver", choices=CFG.SOLVERS, default="sa")
    ap.add_argument("--lasso-max-iter", type=int, default=CFG.LASSO_MAX_ITER,
                    help="LASSO coordinate-descent iterations; a value other than the default runs LASSO only (sensitivity)")
    ap.add_argument("--lasso-pad", choices=["z", "zc"], default=CFG.LASSO_PAD,
                    help="LASSO completion when fewer than K coefficients are non-zero: |z| (default) or |z|*C (sensitivity; LASSO only)")
    a = ap.parse_args()
    RAW = D.Raw(CFG.PBMC_RAW)
    inner(a) if a.stage == "inner" else final(a)
