"""Builds the tables in data_release/ from the benchmark and evaluation results.

  selected_genes.csv    every panel of every final model (CSF: held-out cohort x outer fold; PBMC: unit), all variants
  benchmark_results.csv per-model held-out metrics (ROC-AUC, PR-AUC, Brier, training |rho|, panel size) and the
                        summary rows with 95% donor-bootstrap intervals from the evaluation scripts
Usage: python3 make_data_release.py [output folder]      (default: <repository>/data_release)
data_release/fold_assignments.csv is an input of scripts/01_preprocessing and is not rebuilt here."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

DEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "data_release"
DEST.mkdir(parents=True, exist_ok=True)
VARIANTS = [CFG.PRIMARY, "edger_noIgV_Hb", CFG.variant(CFG.PRIMARY, "exact"), CFG.variant(CFG.PRIMARY, "sa", 20000),
            CFG.variant(CFG.PRIMARY, "sa", lasso_pad="zc")]


def models():
    """(tissue, variant, evaluation, model, folder) of every final model that exists."""
    for v in VARIANTS:
        for ho in HOS:
            for k in OUTERS:
                d = csf_dir(v) / "outer" / ho / f"outer_{k}"
                if (d / "panels.csv").exists():
                    yield "CSF", v, ho, f"outer_{k}", d
        for u in CFG.PBMC_UNITS:
            d = pbmc_dir(v) / "final" / u
            if (d / "panels.csv").exists():
                yield "PBMC", v, u, u, d


def metrics(p, dx, cohort):
    y = (dx == "MS").astype(float).values; sc = p.values
    ms, hd = sc[y == 1], sc[y == 0]
    out = dict(n_MS=int(y.sum()), n_control=int((1 - y).sum()))
    if len(ms) and len(hd):
        out.update(roc_auc=float(auc_w(ms, hd, np.ones(len(ms)), np.ones(len(hd)))[0]), pr_auc=float(ap_w(sc, y, np.ones(len(sc)))[0]))
    out["brier"] = float(brier_w(sc, y, np.ones(len(sc)))[0])
    w = D.within_auc(pd.DataFrame(dict(prob_MS=sc, diagnosis=dx.values, cohort=cohort.values)))
    out["within_cohort_roc_auc"] = float(w)
    return out


pbmc_info = pd.read_csv(CFG.PBMC_RAW / "donor_info.csv").set_index("donor_id")
genes, rows = [], []
for tissue, v, ev, model, d in models():
    Pn = pd.read_csv(d / "panels.csv", keep_default_na=False)
    R = pd.read_csv(d / "rho.csv", keep_default_na=False)
    genes += [dict(tissue=tissue, variant=v, evaluation=ev, model=model, method=r.method, setting=r.setting,
                   cell_type=r.cell_type, gene=r.gene) for r in Pn.itertuples()]
    for m, (h, I) in model_outputs(d, Pn.method.unique()).items():
        p, _ = ensemble(h, I, "weighted")
        dx = h.groupby("donor").diagnosis.first().loc[p.index]
        coh = (pd.Series(p.index, index=p.index).map(pbmc_info.cohort) if tissue == "PBMC"
               else pd.Series(ev, index=p.index))
        rows.append(dict(tissue=tissue, variant=v, evaluation=ev, model=model, method=m, setting=h.setting.iloc[0],
                         **metrics(p, dx, coh), train_rho=R[R.method == m].rho.mean(),
                         panel_size=Pn[Pn.method == m].groupby("cell_type").size().mean(),
                         n_cell_types=Pn[Pn.method == m].cell_type.nunique()))
G = pd.DataFrame(genes); G.to_csv(DEST / "selected_genes.csv", index=False)
M = pd.DataFrame(rows); M.insert(0, "level", "model")

# summary rows (point estimate and 95% interval) from the evaluation outputs
S = []
for scheme_file, v in [("csf_performance/metrics_ci_weighted.csv", CFG.PRIMARY)]:
    T = pd.read_csv(EVAL / scheme_file)
    S += [dict(level="summary", tissue="CSF", variant=v, evaluation=r.cohort, method=r.method, metric=r.metric, value=r.value, lo=r.lo, hi=r.hi)
          for r in T.itertuples()]
T = pd.read_csv(EVAL / "pbmc" / "cv_metrics.csv"); T = T[T.scheme == "weighted"]
S += [dict(level="summary", tissue="PBMC", variant=CFG.PRIMARY, evaluation="nested CV (3 repeats)", method=r.method, metric=r.metric,
           value=r.value, lo=r.lo, hi=r.hi) for r in T.itertuples()]
for tissue, f in [("CSF", EVAL / "csf_redundancy" / "train_rho_paired.csv"), ("PBMC", EVAL / "pbmc" / "cv_rho_paired.csv")]:
    T = pd.read_csv(f)
    S += [dict(level="summary", tissue=tissue, variant=CFG.PRIMARY, evaluation="training donors", method=f"QUBO - {r.comparator}",
               metric="train_rho_difference", value=r.diff_ct_mean, lo=r.lo, hi=r.hi, p_holm=r.p_holm) for r in T.itertuples()]
B_ = pd.concat([M, pd.DataFrame(S)], ignore_index=True)
B_.to_csv(DEST / "benchmark_results.csv", index=False)
print(f"selected_genes.csv: {len(G)} rows; benchmark_results.csv: {len(M)} model rows + {len(S)} summary rows -> {DEST}")
print(M[M.variant == CFG.PRIMARY].groupby(["tissue", "method"])[["roc_auc", "within_cohort_roc_auc", "brier", "train_rho", "panel_size"]].mean().round(3).to_string())
