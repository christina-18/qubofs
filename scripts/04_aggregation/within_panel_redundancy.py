"""Compute within-panel redundancy (mean |Pearson rho|) per method.

For every selected panel (method × holdout × fold × cell type), computes the mean
absolute Pearson correlation between the selected genes on the training-donor
pseudobulk (mean / log-normalised) profiles, then averages per method. This is
the within-panel |rho| reported in the main manuscript performance table
(Table 2) — generated from outputs, not hardcoded.

Outputs: qubo_run/within_panel_redundancy_summary.csv  (method, mean_abs_rho, n_panels)
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
PSB = os.environ.get("QUBOFS_PSEUDOBULK_SUBDIR", "pseudobulk_v5_compartment")  # canonical default (reproduce.sh exports this for all stages)
TISSUE = os.environ.get("QUBOFS_TISSUES", "CSF").split(",")[0]
OUT_TAG = RUN_TAG + ("" if TISSUE == "CSF" else "_" + TISSUE)
_LEGACY = (TISSUE == "CSF")

HOLDOUTS = {
    "Pappalardo": ("", PSB),
    "Heming": ("_holdout_Heming", f"{PSB}_holdout_osmzhlab_MS_ence_cov"),
    "Ramesh": ("_holdout_Ramesh", f"{PSB}_holdout_PRJNA549712_MS_PBMC_UCSF"),
}

# Restrict to a subset of holdouts (e.g. blood benchmark: QUBOFS_HOLDOUTS=Pappalardo,Ramesh)
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


def panel_rho(genes, data_sub, ct, fold):
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
    return float(np.nanmean(np.abs(C[iu])))


def main():
    acc = {}
    per_panel = []
    for ho, (sub, data_sub) in HOLDOUTS.items():
        for sg in glob.glob(str(RUN / f"{RUN_TAG}{sub}" / TISSUE / "selected_genes_folds_*.csv")):
            df = pd.read_csv(sg)
            for (m, fold, ct), grp in df.groupby(["method", "fold", "cell_type"]):
                r = panel_rho(grp["gene"].tolist(), data_sub, ct, fold)
                if not np.isnan(r):
                    acc.setdefault(m, []).append(r)
                    per_panel.append({"method": m, "holdout": ho, "fold": fold,
                                      "cell_type": ct, "rho": r})
    if not acc:
        raise FileNotFoundError(
            f"No selected_genes found under {RUN}/{RUN_TAG}*. Run 03_selection first."
        )
    out = pd.DataFrame(
        [{"method": m, "mean_abs_rho": np.mean(v), "n_panels": len(v)} for m, v in acc.items()]
    ).sort_values("mean_abs_rho")
    # Tag-suffixed (authoritative per DE source; never clobbered by another run)
    # plus the legacy fixed names (kept for backward compatibility).
    path = RUN / f"within_panel_redundancy_summary_{OUT_TAG}.csv"
    out.to_csv(path, index=False)
    if _LEGACY: out.to_csv(RUN / "within_panel_redundancy_summary.csv", index=False)
    # per-panel values (for paired statistical tests in bootstrap_stats.py)
    pp = pd.DataFrame(per_panel)
    pp.to_csv(RUN / f"within_panel_redundancy_perpanel_{OUT_TAG}.csv", index=False)
    if _LEGACY: pp.to_csv(RUN / "within_panel_redundancy_perpanel.csv", index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nWrote {path}" + (" (+ legacy within_panel_redundancy_summary.csv)" if _LEGACY else " (tissue-suffixed; legacy not overwritten)"))
    print(f"Wrote within_panel_redundancy_perpanel_{OUT_TAG}.csv" + (" (+ legacy)" if _LEGACY else "") + f" ({len(pp)} panels)")


if __name__ == "__main__":
    main()
