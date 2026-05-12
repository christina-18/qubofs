"""Aggregate v9 mixed results with multiple soft-voting schemes.

Uses per-cell-type per-donor predictions saved by the v9 pipeline to compute:

  (a) UNIFORM:   p_donor = mean over cell types where donor has prediction
  (b) AUC_GATED: drop cell types with inner_val_auc < threshold (default 0.55),
                 then uniform mean
  (c) AUC_WEIGHTED: weight w_c = max(0, inner_val_auc_c - 0.5),
                    p_donor = sum(w_c * p_c) / sum(w_c)

The inner-val-AUC weights are derived ONLY from training cohorts (no test-set
peeking) — this is principled, similar to Stacking / Super-Learner ensembles.

Outputs:
  qubo_run_v6/v9mixed_summary_three_schemes.csv
  qubo_run_v6/v9mixed_v6_v8_v9_compare.csv
"""
import pandas as pd
import numpy as np
import glob
from pathlib import Path
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent / "qubo_run_v6"
THRESHOLD_GATED = 0.55  # for AUC_gated mode

def load_v9_dirs():
    return sorted(ROOT.glob("v9mixed_bio_edger*"))

def load_per_ct(d):
    files = sorted(d.glob("CSF/held_per_ct_predictions_*.csv"))
    if not files:
        return None
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

def load_per_ct_diag(d):
    files = sorted(d.glob("CSF/per_ct_diag_*.csv"))
    if not files:
        return None
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

def auc_safe(y, p):
    if len(set(y)) < 2 or len(y) < 2:
        return np.nan
    try:
        return roc_auc_score(y, p)
    except Exception:
        return np.nan

records = []
v9_dirs = load_v9_dirs()
print(f"Found {len(v9_dirs)} v9 result dirs")

for d in v9_dirs:
    if "_holdout_" in d.name:
        cohort = d.name.split("_holdout_")[1]
    else:
        cohort = "Pappalardo"

    per_ct = load_per_ct(d)
    diag = load_per_ct_diag(d)
    if per_ct is None or diag is None:
        print(f"  {d.name}: missing per-ct files, skipping")
        continue
    print(f"  {d.name}: {len(per_ct)} per-ct rows, {len(diag)} diag rows")

    # diag has columns: tissue, fold, cell_type, diag_method (always QUBO),
    # n_genes, K_chosen, n_val, val_auc, n_held, held_auc
    # We need val_auc per (cell_type, fold) for the QUBO method.
    # For other methods, we don't have per-cell-type val_auc saved by default.
    # As a fallback, we use the QUBO val_auc weights for all methods (simple).
    # If the user wants method-specific weights, we'd need to extend the diag.

    diag_qubo = diag[diag.diag_method == 'QUBO'][['fold','cell_type','val_auc']].copy()

    methods_in_data = sorted(per_ct.method.unique()) if 'method' in per_ct.columns else ['QUBO']

    for method in methods_in_data:
        sub = per_ct[per_ct.method == method] if 'method' in per_ct.columns else per_ct
        for fold, fg in sub.groupby('fold'):
            # weights from QUBO inner CV (proxy)
            wmap = dict(zip(diag_qubo[diag_qubo.fold==fold].cell_type,
                            diag_qubo[diag_qubo.fold==fold].val_auc))
            # for each donor, compute aggregated probability under three schemes
            for donor, dg in fg.groupby('donor'):
                y_label = 1 if dg.diagnosis.iloc[0] == 'MS' else 0
                # per-cell-type predictions
                ct_p = dict(zip(dg.cell_type, dg.prob_MS))

                # (a) UNIFORM
                p_uni = float(np.mean(list(ct_p.values())))
                # (b) AUC GATED
                kept = [(c, p) for c, p in ct_p.items()
                        if pd.notna(wmap.get(c, np.nan))
                        and wmap.get(c, 0) >= THRESHOLD_GATED]
                p_gated = float(np.mean([p for _, p in kept])) if kept else np.nan
                # (c) AUC WEIGHTED
                weighted_pairs = [(p, max(0, wmap.get(c, 0.5) - 0.5))
                                  for c, p in ct_p.items()
                                  if pd.notna(wmap.get(c, np.nan))]
                w_sum = sum(w for _, w in weighted_pairs)
                if w_sum > 0:
                    p_weighted = sum(p * w for p, w in weighted_pairs) / w_sum
                else:
                    p_weighted = np.nan

                records.append({
                    'cohort': cohort, 'fold': fold, 'method': method,
                    'donor': donor, 'label': y_label,
                    'n_ct_used': len(ct_p),
                    'p_uniform': p_uni,
                    'p_gated': p_gated,
                    'p_weighted': p_weighted,
                })

if not records:
    print("No records aggregated.")
else:
    df = pd.DataFrame(records)
    print(f"\n{len(df)} (donor x fold x method x cohort) rows")

    # AUC per cohort × fold × method × scheme
    auc_rows = []
    for (cohort, fold, method), g in df.groupby(['cohort','fold','method']):
        if len(set(g.label)) < 2:
            continue
        auc_rows.append({
            'cohort': cohort, 'fold': fold, 'method': method,
            'auc_uniform':  auc_safe(g.label, g.p_uniform),
            'auc_gated':    auc_safe(g.label, g.p_gated.fillna(g.p_uniform)),
            'auc_weighted': auc_safe(g.label, g.p_weighted.fillna(g.p_uniform)),
        })
    A = pd.DataFrame(auc_rows)

    # Per-cohort mean AUC then std across cohorts
    print("\n=== Mean held-AUC per (method × scheme), σ across 3 cohorts ===")
    for scheme in ['auc_uniform','auc_gated','auc_weighted']:
        per_cohort = A.groupby(['method','cohort'])[scheme].mean().reset_index()
        summ = per_cohort.groupby('method')[scheme].agg(['mean','std','count']).round(3)
        print(f"\n--- {scheme} ---")
        print(summ.sort_values('mean', ascending=False))

    A.to_csv(ROOT / "v9mixed_summary_three_schemes.csv", index=False)
    print(f"\nWrote {ROOT / 'v9mixed_summary_three_schemes.csv'}")

    # Final comparison table: V6 vs V8 vs V9 (uniform/gated/weighted)
    methods = ['QUBO', 'ElasticNet', 'LASSO', 'DE_top', 'HVG']
    print("\n=== V6 vs V8 vs V9 (3 schemes): held-AUC mean ± σ ===")

    def load_summary(tag):
        files = glob.glob(str(ROOT / f"{tag}*/CSF/fold_metrics*.csv"))
        rows = []
        for f in files:
            cohort = ('Pappalardo' if '_holdout_' not in f
                      else f.split('_holdout_')[1].split('/')[0])
            sub = pd.read_csv(f)
            sub['cohort'] = cohort
            rows.append(sub)
        return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    v6 = load_summary('v6entrue_bio_edger')
    v8 = load_summary('v8mixed_bio_edger')
    def auc_summary(df):
        per = df.groupby(['method','cohort']).held_auc.mean().reset_index()
        return per.groupby('method').held_auc.agg(['mean','std']).round(3)

    rows = []
    for m in methods:
        v6m = auc_summary(v6).loc[m] if m in auc_summary(v6).index else None
        v8m = auc_summary(v8).loc[m] if m in auc_summary(v8).index else None
        # v9: uniform / gated / weighted
        v9_u = A[A.method == m].groupby('cohort').auc_uniform.mean()
        v9_g = A[A.method == m].groupby('cohort').auc_gated.mean()
        v9_w = A[A.method == m].groupby('cohort').auc_weighted.mean()
        rows.append({
            'method': m,
            'v6_AUC':  v6m['mean'] if v6m is not None else np.nan,
            'v6_σ':    v6m['std']  if v6m is not None else np.nan,
            'v8_AUC':  v8m['mean'] if v8m is not None else np.nan,
            'v8_σ':    v8m['std']  if v8m is not None else np.nan,
            'v9_uniform_AUC':  v9_u.mean(),
            'v9_uniform_σ':    v9_u.std(),
            'v9_gated_AUC':    v9_g.mean(),
            'v9_gated_σ':      v9_g.std(),
            'v9_weighted_AUC': v9_w.mean(),
            'v9_weighted_σ':   v9_w.std(),
        })
    cmp = pd.DataFrame(rows).round(3)
    print(cmp.to_string(index=False))
    cmp.to_csv(ROOT / "v9mixed_v6_v8_v9_compare.csv", index=False)
    print(f"\nWrote {ROOT / 'v9mixed_v6_v8_v9_compare.csv'}")
