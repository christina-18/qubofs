"""Collect the panel-size (K) sweep into one summary CSV (Supplementary Figure S2).

Reads the per-K runs produced by the K-sweep driver, i.e. pipeline runs with
QUBOFS_PIPELINE_RUN_TAG=sweepK<K> and QUBOFS_K_SWEEP=<K>, which write to
  qubo_run/sweepK<K>_bio_edger_counts{,_holdout_Heming,_holdout_Ramesh}/CSF/
The primary benchmark uses a fixed panel size (QUBOFS_FIXED_K=10); this sweep
covers the panel sizes reported in the Supplementary, K in {5, 10, 15} by default
(override with QUBOFS_SWEEP_KS). For each K and method it computes the mean
held-out ROC-AUC, MCC and Macro-F1 (over the 15 cohort-by-fold values) and the
mean within-panel |rho| (reusing the canonical panel_rho), with no hardcoded
numbers.

Output: qubo_run/sweep_all_methods_K_summary.csv
"""
import os
import re
import glob
import sys
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(os.environ.get("QUBOFS_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
RUN = PROJ / "qubo_run"
TISSUE = "CSF"
DEG = os.environ.get("QUBOFS_DEG_SOURCE", "edger_counts")

# reuse the canonical within-panel redundancy computation
sys.path.insert(0, str(Path(__file__).resolve().parent))
from within_panel_redundancy import panel_rho, HOLDOUTS  # noqa: E402

HO_SUFFIX = {"Pappalardo": "", "Heming": "_holdout_Heming", "Ramesh": "_holdout_Ramesh"}


# Run-tag prefix for the per-K runs (default "sweepK"); override with
# QUBOFS_SWEEP_PREFIX to read an alternative family of K-sweep runs.
SWEEP_PREFIX = os.environ.get("QUBOFS_SWEEP_PREFIX", "sweepK")

# Panel sizes reported in the Supplementary K-sweep (default 5, 10, 15). Any
# other per-K runs present under qubo_run/ are ignored so the figure matches the
# manuscript.
K_ALLOWED = [int(x) for x in os.environ.get("QUBOFS_SWEEP_KS", "5,10,15").split(",") if x.strip()]


def collect():
    bases = sorted(glob.glob(str(RUN / f"{SWEEP_PREFIX}*_bio_{DEG}")))
    if not bases:
        raise FileNotFoundError(
            f"No {SWEEP_PREFIX}*_bio_{DEG} runs under {RUN}. Run the K-sweep driver first "
            f"(QUBOFS_PIPELINE_RUN_TAG={SWEEP_PREFIX}<K>, QUBOFS_FIXED_K=<K>)."
        )
    rows = []
    for base in bases:
        tag = Path(base).name                       # e.g. sweepK5_bio_edger_counts
        K = int(re.search(rf"{re.escape(SWEEP_PREFIX)}(\d+)_", tag).group(1))
        if K not in K_ALLOWED:
            continue
        # ---- performance: mean over the 15 cohort-by-fold values per method ----
        fm = []
        for ho, suf in HO_SUFFIX.items():
            for fp in glob.glob(str(RUN / f"{tag}{suf}" / TISSUE / "fold_metrics_folds_*.csv")):
                fm.append(pd.read_csv(fp))
        if not fm:
            continue
        fm = pd.concat(fm, ignore_index=True)
        agg_kwargs = dict(roc_auc=("held_auc", "mean"), mcc=("held_mcc", "mean"))
        if "held_macro_f1" in fm.columns:
            agg_kwargs["macro_f1"] = ("held_macro_f1", "mean")
        perf = fm.groupby("method").agg(**agg_kwargs)
        if "macro_f1" not in perf.columns:
            perf["macro_f1"] = np.nan
        # ---- within-panel |rho| per method (canonical panel_rho) ----
        rho_acc = {}
        for ho, suf in HO_SUFFIX.items():
            data_sub = HOLDOUTS[ho][1]
            for sg in glob.glob(str(RUN / f"{tag}{suf}" / TISSUE / "selected_genes_folds_*.csv")):
                df = pd.read_csv(sg)
                for (m, fold, ct), grp in df.groupby(["method", "fold", "cell_type"]):
                    r = panel_rho(grp["gene"].tolist(), data_sub, ct, fold)
                    if not np.isnan(r):
                        rho_acc.setdefault(m, []).append(r)
        for m in perf.index:
            rho_vals = rho_acc.get(m, [])
            rows.append({"K": K, "method": m,
                         "roc_auc": perf.loc[m, "roc_auc"],
                         "mcc": perf.loc[m, "mcc"],
                         "macro_f1": perf.loc[m, "macro_f1"],
                         "within_panel_rho": (np.mean(rho_vals) if rho_vals else np.nan),
                         "n_panels": len(rho_vals)})
    out = pd.DataFrame(rows).sort_values(["method", "K"]).reset_index(drop=True)
    path = RUN / "sweep_all_methods_K_summary.csv"
    out.to_csv(path, index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nWrote {path}")


if __name__ == "__main__":
    collect()
