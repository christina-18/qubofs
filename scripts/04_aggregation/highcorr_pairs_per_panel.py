"""Count highly-correlated gene pairs per panel, per method (Supplementary Table S8).

Assay-design view of within-panel redundancy. For every selected panel
(method x holdout x fold x cell type) it counts the within-panel gene pairs with
|Pearson r| > THRESHOLD on the training-donor pseudobulk (mean/log-normalised)
profiles, then averages per method. Shares the exact matrix-loading and
correlation logic of within_panel_redundancy.py; only the per-panel reducer
differs (count of highly-correlated pairs instead of mean |rho|).

Outputs:
  qubo_run/highcorr_pairs_summary_<tag>.csv   (method, mean_highcorr_pairs_per_panel, n_panels, threshold_abs_r)
  qubo_run/highcorr_pairs_perpanel_<tag>.csv  (method, holdout, fold, cell_type, n_highcorr_pairs)
"""
import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(os.environ.get("QUBOFS_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
RUN = PROJ / "qubo_run"
DATA = PROJ / "data"
RUN_TAG = os.environ.get("QUBOFS_RUN_TAG", "primary_bio_edger_counts")
PSB = os.environ.get("QUBOFS_PSEUDOBULK_SUBDIR", "pseudobulk_v5_compartment")
TISSUE = os.environ.get("QUBOFS_TISSUES", "CSF").split(",")[0]
THRESHOLD = float(os.environ.get("QUBOFS_HIGHCORR_THRESHOLD", "0.70"))
OUT_TAG = RUN_TAG + ("" if TISSUE == "CSF" else "_" + TISSUE)
_LEGACY = (TISSUE == "CSF")

HOLDOUTS = {
    "Pappalardo": ("", PSB),
    "Heming": ("_holdout_Heming", f"{PSB}_holdout_osmzhlab_MS_ence_cov"),
    "Ramesh": ("_holdout_Ramesh", f"{PSB}_holdout_PRJNA549712_MS_PBMC_UCSF"),
}
_HO_ENV = os.environ.get("QUBOFS_HOLDOUTS")
if _HO_ENV:
    HOLDOUTS = {k: v for k, v in HOLDOUTS.items() if k in _HO_ENV.split(",")}

_cache = {}


def load_mtx(path):
    """Minimal MatrixMarket coordinate reader -> dense ndarray (genes x donors)."""
    if path in _cache:
        return _cache[path]
    if not os.path.exists(path):
        _cache[path] = None
        return None
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while lines[i].startswith("%"):
        i += 1
    nr, nc, _ = map(int, lines[i].split())
    i += 1
    M = np.zeros((nr, nc))
    for ln in lines[i:]:
        r, c, v = ln.split()
        M[int(r) - 1, int(c) - 1] = float(v)
    _cache[path] = M
    return M


def panel_highcorr(genes, data_sub, ct, fold):
    """Count |Pearson r| > THRESHOLD pairs among panel genes (same loading as panel_rho)."""
    d = DATA / data_sub / ct / TISSUE / f"fold_{fold}"
    M = load_mtx(str(d / "train_pb_mean.mtx"))
    if M is None:
        return np.nan
    rows = pd.read_csv(d / "train_pb_mean_rows.csv")["gene"].tolist()
    pos = {g: i for i, g in enumerate(rows)}
    idx = [pos[g] for g in genes if g in pos]
    if len(idx) < 2:
        return np.nan
    X = M[idx, :]
    X = X[X.std(axis=1) > 0]
    if X.shape[0] < 2:
        return np.nan
    C = np.corrcoef(X)
    iu = np.triu_indices_from(C, k=1)
    vals = np.abs(C[iu])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float((vals > THRESHOLD).sum())


def main():
    acc = {}
    per_panel = []
    for ho, (sub, data_sub) in HOLDOUTS.items():
        for sg in glob.glob(str(RUN / f"{RUN_TAG}{sub}" / TISSUE / "selected_genes_folds_*.csv")):
            df = pd.read_csv(sg)
            for (m, fold, ct), grp in df.groupby(["method", "fold", "cell_type"]):
                n = panel_highcorr(grp["gene"].tolist(), data_sub, ct, fold)
                if not np.isnan(n):
                    acc.setdefault(m, []).append(n)
                    per_panel.append({"method": m, "holdout": ho, "fold": fold,
                                      "cell_type": ct, "n_highcorr_pairs": n})
    if not acc:
        raise FileNotFoundError(
            f"No selected_genes found under {RUN}/{RUN_TAG}*. Run 03_selection first."
        )
    out = pd.DataFrame(
        [{"method": m, "mean_highcorr_pairs_per_panel": np.mean(v),
          "n_panels": len(v), "threshold_abs_r": THRESHOLD} for m, v in acc.items()]
    ).sort_values("mean_highcorr_pairs_per_panel")
    path = RUN / f"highcorr_pairs_summary_{OUT_TAG}.csv"
    out.to_csv(path, index=False)
    if _LEGACY:
        out.to_csv(RUN / "highcorr_pairs_summary.csv", index=False)
    pp = pd.DataFrame(per_panel)
    pp.to_csv(RUN / f"highcorr_pairs_perpanel_{OUT_TAG}.csv", index=False)
    if _LEGACY:
        pp.to_csv(RUN / "highcorr_pairs_perpanel.csv", index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\nWrote {path}" + (" (+ legacy highcorr_pairs_summary.csv)" if _LEGACY else ""))
    print(f"Wrote highcorr_pairs_perpanel_{OUT_TAG}.csv ({len(pp)} panels)")


if __name__ == "__main__":
    main()
