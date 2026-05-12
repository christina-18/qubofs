"""
qubo_pipeline_v6.py
====================
Per-cell-type QUBO + K grid search + soft-voting ensemble (8 cell types).

Major differences from v5:
  - QUBO is solved INDEPENDENTLY per cell type (Option B)
    each cell type gets its own selected gene set (size K_c)
  - K is part of the inner-CV grid search ({10, 20, 30})
  - Classifier is trained PER CELL TYPE (no concatenated features)
  - Final prediction = soft voting (mean of cell-type probabilities)

Usage:
    python3 qubo_pipeline_v6.py <holdout> <deg> <tissue> <fold1> [fold2 ...]
"""
from __future__ import annotations

import os
import re
import sys
import zlib

def _ct_seed(ct):
    """Deterministic hash for cell-type string (replaces non-reproducible hash())."""
    return int(zlib.crc32(ct.encode("utf-8")) % 10000)
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

# ============================================================
# Biology filter (housekeeping / mito / ribosomal exclusion)
# ============================================================
HK_PATTERN = re.compile(
    # =========================================================================
    #  Pre-selection biology filter (manuscript/v6entrue version).
    #  Aligned with current best practice (Heumos et al. Nat Rev Genet 2023;
    #  Luecken & Theis Mol Syst Biol 2019).
    #
    #  Genes retained for documented biological function (not pure housekeeping):
    #    * RPLP0/1/2 — acidic ribosomal stalk paralogs; specialized ribosome
    #                 programs (Genuth & Barna 2018) + immunomodulatory roles
    #                 (Wang et al. 2020).
    #    * RPSA     — encodes the 67-kDa laminin receptor (LamR/LRP); critical
    #                 for leukocyte adhesion and BBB transmigration (Nelson
    #                 et al. 2008) — directly relevant to MS pathobiology.
    #    * Tier-3 translation/RNA-binding genes (TPT1, EEF2, EIF*, TMSB4X/10,
    #                 HNRN*) retained for regulatory/signaling functions
    #                 beyond constitutive expression.
    #  See Methods_section.md for full rationale and references.
    # =========================================================================
    r"^(MT-|MTRNR|MTATP|MTND|"             # mitochondrial
    r"RPL[0-9]|RPS[0-9]|MRPL|MRPS|"        # ribosomal (RPLP*, RPSA retained)
    r"HSP[A0-9]|HSPB|HSPA|HSPD|"           # heat shock (HSPA/B/D + HSP9*)
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|"  # classical housekeeping
    r"MALAT1$|NEAT1$|XIST$|TSIX$|"         # nuclear lncRNA + X-inactivation
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|"    # uncharacterized / pseudogene loci
    r"MIR[0-9]|RNU[0-9]|SNORA|SNORD)"      # small RNAs (poly-A unreliable)
)
def is_biology_gene(g: str) -> bool:
    """Return True if gene name passes biology filter (excludes HK/mito/ribosomal)."""
    return not bool(HK_PATTERN.match(str(g)))

sys.path.insert(0, str(Path(__file__).parent))
from qubo_utils_v5 import (
    load_fold, build_score_and_redundancy, build_score_and_redundancy_MI,
    build_qubo, solve_qubo_sa,
    cohort_variance_per_gene,
    LogRegL2, LogRegL1, LogRegElasticNet, LDA, standardize,
    roc_auc, average_precision, acc_f1, mcc_score, jaccard,
)

# ============================================================
# Configuration
# ============================================================
# PROJECT_ROOT: derive from this script's location (scripts/ is under project root)
# This makes the pipeline portable between sandbox and the user's Mac.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}

# 8 cell types (extended)
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
TISSUES = ["CSF", "PBMC", "ALL"]
FOLDS = [1, 2, 3, 4, 5]

# K grid (per-cell-type cardinality)
K_GRID = [10, 20, 30]
LAMBDA_VALS = [2.0, 5.0]
GAMMA_VALS  = [1.0]
SA_READS_GRID = 8
SA_SWEEPS_GRID = 200
SA_READS = 30
SA_SWEEPS = 600
SCORE_FN = "t_squared"
ALPHA_BATCH = 1.0
INNER_CV_FOLDS = 3
N_PER_CELL_TYPE = 100   # candidate set per cell type (top-N from tstats)
SEED = 42

# Classifier per cell type (single classifier choice for ensemble simplicity)
BASE_CLF_FACTORY = lambda: LogRegL2(C=1.0, max_iter=200)
ENSEMBLE_AGG = "mean"        # "mean" | "weighted_auc" | "stacking"
META_CLF_FACTORY = lambda: LogRegL2(C=0.5, max_iter=200)   # for stacking

RUN_TAG = "v6full"
HOLDOUT_NAME = "Pappalardo"
DEG_SOURCE = "edger"     # default for v6 (edgeR is primary main-stream method)
BIOLOGY_FILTER = False   # if True, drop housekeeping/mito/ribosomal at candidate stage

# v9: MI-based QUBO (Romero et al. 2025) — replaces |t| relevance and |Pearson|
# redundancy with mutual information (relevance: I(gene_expr, MS-vs-HD label);
# redundancy: I(gene_i, gene_j)). Both use quantile binning with N_BINS bins.
USE_MI_QUBO = False
MI_N_BINS = 5

METHODS = ["QUBO", "QUBO_consensus", "QUBO_hybrid", "DE_top", "HVG", "LASSO", "ElasticNet"]
HYBRID_TOP_N = 30   # Pre-filter candidate pool to top-N by |t| for QUBO_hybrid
CONSENSUS_N_RUNS = 10   # # of independent SA runs for QUBO_consensus
CONSENSUS_SEED_STEP = 7919   # large prime to spread seeds across runs
# Lighter SA settings for the 10 consensus runs (we don't need each
# individual run to be perfectly converged — only the FREQUENCY across runs)
CONSENSUS_SA_READS = 10
CONSENSUS_SA_SWEEPS = 250


def _data_root(holdout_name: str) -> Path:
    if holdout_name == "Pappalardo":
        return PROJECT_ROOT / "data" / "pseudobulk_v5_compartment"
    return PROJECT_ROOT / "data" / f"pseudobulk_v5_compartment_holdout_{HOLDOUT_PRJ_MAP[holdout_name]}"


# ============================================================
# Per-cell-type QUBO
# ============================================================
def select_genes_for_cell_type(bundle, candidates_ct, k, lam, gamma,
                                seed, score_fn, alpha_batch,
                                sa_reads, sa_sweeps):
    """Run QUBO on a single cell type's bundle and candidate list."""
    one_bundle = {bundle["cell_type"]: bundle}
    if USE_MI_QUBO:
        s_raw, R, _ = build_score_and_redundancy_MI(
            one_bundle, candidates_ct,
            n_bins=MI_N_BINS,
            score_agg="sum", redundancy_agg="max")
    else:
        s_raw, R, _ = build_score_and_redundancy(
            one_bundle, candidates_ct,
            score_agg="sum", redundancy_agg="max", score_fn=score_fn)

    def _mm(v):
        if v.max() > v.min():
            return (v - v.min()) / (v.max() - v.min())
        return np.zeros_like(v)
    s_norm = _mm(s_raw)

    if alpha_batch > 0:
        cv = cohort_variance_per_gene(one_bundle, candidates_ct)
        cv_norm = _mm(cv)
        s = _mm(s_norm - alpha_batch * cv_norm)
    else:
        s = s_norm

    Q = build_qubo(s, R, k=k, lam=lam, gamma=gamma)
    rng = np.random.default_rng(seed)
    x, E = solve_qubo_sa(Q, k=k, n_reads=sa_reads, n_sweeps=sa_sweeps, rng=rng)
    selected_idx = np.where(x == 1)[0]
    selected = [candidates_ct[i] for i in selected_idx]
    return selected, E


def candidates_per_cell_type(bundle, n_top=N_PER_CELL_TYPE):
    """Top-N |t| genes from a single cell type's topN list.
    If BIOLOGY_FILTER is True, drop housekeeping/mito/ribosomal genes
    BEFORE taking top-N (so all methods see a biology-only candidate pool).
    """
    if bundle is None or bundle["topN"] is None:
        return []
    # use the full tstats list (more than topN_genes.csv) for re-ranking after filter
    if bundle.get("tstats") is not None and BIOLOGY_FILTER:
        ts = bundle["tstats"].copy()
        ts = ts[ts["gene"].apply(is_biology_gene)]
        ts = ts.sort_values("t", key=lambda s: -s.abs())
        return ts.head(n_top)["gene"].tolist()
    candidates = bundle["topN"]["gene"].tolist()
    if BIOLOGY_FILTER:
        candidates = [g for g in candidates if is_biology_gene(g)]
    return candidates[:n_top]


def candidates_per_cell_type_hybrid(bundle, n_top=HYBRID_TOP_N):
    """For QUBO_hybrid: pre-filter to top-N by |t| (smaller, stronger signal pool).
    Reuses candidates_per_cell_type with smaller n_top.
    """
    return candidates_per_cell_type(bundle, n_top=n_top)


# ============================================================
# Baseline selection methods (per cell type)
# ============================================================
def select_baseline_per_ct(method, bundle, candidates_ct, K):
    """Select K genes from candidates_ct using a baseline method (per cell type)."""
    if not candidates_ct or K <= 0:
        return []
    K = min(K, len(candidates_ct))

    if method == "DE_top":
        # top K |t| (already in tstats order, but candidates_ct is top-100 already)
        ts = bundle["tstats"].set_index("gene") if bundle["tstats"] is not None else None
        if ts is None:
            return list(candidates_ct[:K])
        scores = []
        for g in candidates_ct:
            t = abs(float(ts.loc[g, "t"])) if g in ts.index else 0.0
            scores.append((g, t))
        scores.sort(key=lambda x: -x[1])
        return [g for g, _ in scores[:K]]

    if method == "HVG":
        # variance across train donors for each gene
        tr = bundle.get("train")
        if tr is None:
            return list(candidates_ct[:K])
        gene_pos = {g: i for i, g in enumerate(tr["genes"])}
        var_scores = []
        for g in candidates_ct:
            if g in gene_pos:
                v = float(np.var(tr["X"][gene_pos[g], :]))
                var_scores.append((g, v))
            else:
                var_scores.append((g, 0.0))
        var_scores.sort(key=lambda x: -x[1])
        return [g for g, _ in var_scores[:K]]

    if method == "LASSO":
        # Pure L1 logistic.  Fit LogRegL1 with various C and take top K by |coef|.
        return _fit_l1_or_en_select(
            bundle, candidates_ct, K, model_kind="L1",
        )

    if method == "ElasticNet":
        # True Elastic Net (L1 + L2 mix, l1_ratio=0.5).
        # NOTE: L2 keeps weights non-zero, so we always select top K by |coef|
        #       (matching the same K=20 cardinality as LASSO).
        return _fit_l1_or_en_select(
            bundle, candidates_ct, K, model_kind="EN", l1_ratio=0.5,
        )

    raise ValueError(f"Unknown method: {method}")


def _fit_l1_or_en_select(bundle, candidates_ct, K, model_kind="L1", l1_ratio=0.5):
    """Shared selection helper for LASSO / ElasticNet.

    For L1: tune C to hit ~K nonzero, then take top K by |coef|.
    For EN: try a small C grid, pick the C whose top-K-by-|coef| has the
            most concentrated coefficient mass (proxy for "good support").
    """
    tr_info = build_X_y_per_ct(bundle, candidates_ct, "train")
    if tr_info is None or len(tr_info["X"]) < 4:
        return list(candidates_ct[:K])
    X = tr_info["X"]; y = tr_info["y"]
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    Xz = (X - mu) / sd
    C_grid = [10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    best_genes, best_n = None, 0
    best_genes_en, best_score_en = None, -1.0
    for C in C_grid:
        try:
            if model_kind == "L1":
                clf = LogRegL1(C=C, max_iter=120).fit(Xz, y)
            else:
                clf = LogRegElasticNet(C=C, l1_ratio=l1_ratio, max_iter=120).fit(Xz, y)
            coefs = np.asarray(clf.coef_)
            nz = np.where(np.abs(coefs) > 1e-6)[0]
            if model_kind == "L1":
                # LASSO-style: stop at first C where nonzero >= K
                if len(nz) >= K:
                    abs_coefs = [(tr_info["genes"][i], abs(coefs[i])) for i in nz]
                    abs_coefs.sort(key=lambda x: -x[1])
                    return [g for g, _ in abs_coefs[:K]]
                if len(nz) > best_n:
                    best_n = len(nz)
                    abs_coefs = [(tr_info["genes"][i], abs(coefs[i])) for i in nz]
                    abs_coefs.sort(key=lambda x: -x[1])
                    best_genes = [g for g, _ in abs_coefs]
            else:
                # ElasticNet: rank by |coef| over ALL features (L2 keeps them non-zero)
                # Score the C: top-K mass over total mass (concentration index).
                order = np.argsort(-np.abs(coefs))
                topK = order[:K]
                top_mass = float(np.sum(np.abs(coefs[topK])))
                tot_mass = float(np.sum(np.abs(coefs))) + 1e-12
                conc = top_mass / tot_mass
                if conc > best_score_en:
                    best_score_en = conc
                    best_genes_en = [tr_info["genes"][i] for i in topK]
        except Exception:
            continue
    if model_kind == "EN":
        if best_genes_en is None:
            return list(candidates_ct[:K])
        return best_genes_en[:K]
    # L1 fallback
    if best_genes is None:
        best_genes = []
    if len(best_genes) < K:
        ts = bundle["tstats"].set_index("gene") if bundle["tstats"] is not None else None
        if ts is not None:
            pool = [g for g in candidates_ct if g not in set(best_genes)]
            pool_scores = [(g, abs(float(ts.loc[g, "t"])) if g in ts.index else 0.0)
                            for g in pool]
            pool_scores.sort(key=lambda x: -x[1])
            best_genes = best_genes + [g for g, _ in pool_scores[: K - len(best_genes)]]
    return best_genes[:K]


def grid_search_baseline_per_ct(method, bundle, candidates_ct, fold, seed_base):
    """Inner-CV K grid search for a baseline method on a single cell type."""
    if not candidates_ct:
        return None, []
    n_cand = len(candidates_ct)
    valid_K = sorted(set(min(K, n_cand) for K in K_GRID if K <= n_cand))
    if not valid_K:
        return None, []
    best = None
    log = []
    for K in valid_K:
        sel = select_baseline_per_ct(method, bundle, candidates_ct, K)
        if not sel or len(sel) < 3:
            continue
        auc = inner_cv_auc_one_ct(bundle, sel, n_inner=INNER_CV_FOLDS, seed=seed_base + 100)
        log.append(dict(method=method, K=K, k_actual=len(sel), inner_cv_auc=auc))
        score = -auc if not np.isnan(auc) else float('inf')
        if best is None or score < best["score"]:
            best = dict(method=method, K=K, k_actual=len(sel), score=score,
                        inner_cv_auc=auc, selected=sel)
    return best, log


# ============================================================
# Per-cell-type feature build & inner CV
# ============================================================
def build_X_y_per_ct(bundle, gene_subset, split):
    """Returns (X, donor_ids, y, meta) for a single cell type & split.
    X: donors x genes (numpy float64). NaN-fill with 0.
    """
    if bundle is None or bundle.get(split) is None:
        return None
    sb = bundle[split]
    gene_pos = {g: i for i, g in enumerate(sb["genes"])}
    present = [g for g in gene_subset if g in gene_pos]
    if not present:
        return None
    idx = [gene_pos[g] for g in present]
    X = sb["X"][idx, :].T.astype(np.float64)  # donors x genes
    donors = sb["donors"]
    meta = sb["meta"]
    if meta is None:
        return None
    diag = dict(zip(meta["donor_id"], meta["diagnosis"]))
    y = np.array([1 if diag.get(d) == "MS" else 0 for d in donors])
    return dict(X=X, donors=donors, y=y, genes=present, meta=meta)


def inner_cv_auc_one_ct(bundle, gene_subset, n_inner=3, seed=0):
    """Inner CV AUC for a single cell type's classifier on its train pseudobulk."""
    info = build_X_y_per_ct(bundle, gene_subset, "train")
    if info is None or len(info["X"]) < 2 * n_inner:
        return np.nan
    X, y = info["X"], info["y"]
    n = len(X)
    rng = np.random.default_rng(seed)
    pos = [i for i in range(n) if y[i] == 1]
    neg = [i for i in range(n) if y[i] == 0]
    rng.shuffle(pos); rng.shuffle(neg)
    fold_assign = [[] for _ in range(n_inner)]
    for i, p in enumerate(pos): fold_assign[i % n_inner].append(p)
    for i, p in enumerate(neg): fold_assign[i % n_inner].append(p)
    aucs = []
    for k in range(n_inner):
        v = fold_assign[k]
        t = [i for i in range(n) if i not in v]
        if len(set(y[v])) < 2 or len(set(y[t])) < 2:
            continue
        Xt = X[t]; yt = y[t]
        Xv = X[v]; yv = y[v]
        mu = Xt.mean(0); sd = Xt.std(0); sd[sd == 0] = 1.0
        Xtz = (Xt - mu) / sd
        Xvz = (Xv - mu) / sd
        clf = BASE_CLF_FACTORY().fit(Xtz, yt)
        p = clf.predict_proba(Xvz)
        a = roc_auc(yv, p)
        if not np.isnan(a):
            aucs.append(a)
    return float(np.mean(aucs)) if aucs else np.nan


def grid_search_per_ct(bundle, candidates_ct, fold, seed_base):
    """Search over (K, λ, γ) for one cell type, picking by inner CV AUC.
    K is clamped to min(K, len(candidates_ct)-1) to avoid sampling errors.
    """
    best = None
    log = []
    n_cand = len(candidates_ct)
    valid_K = sorted(set(min(K, n_cand - 1) for K in K_GRID if K < n_cand))
    if not valid_K:
        return None, log
    for K in valid_K:
        for lam in LAMBDA_VALS:
            for gamma in GAMMA_VALS:
                seed = seed_base + K * 10 + int(lam) + int(gamma * 10)
                sel, E = select_genes_for_cell_type(
                    bundle, candidates_ct, k=K, lam=lam, gamma=gamma,
                    seed=seed, score_fn=SCORE_FN, alpha_batch=ALPHA_BATCH,
                    sa_reads=SA_READS_GRID, sa_sweeps=SA_SWEEPS_GRID)
                if len(sel) < 3:
                    continue
                auc = inner_cv_auc_one_ct(bundle, sel, n_inner=INNER_CV_FOLDS,
                                          seed=seed_base + 100)
                log.append(dict(K=K, lam=lam, gamma=gamma, k_actual=len(sel),
                                energy=E, inner_cv_auc=auc, n_genes=len(sel)))
                score = -auc if not np.isnan(auc) else float('inf')
                if best is None or score < best["score"]:
                    best = dict(K=K, lam=lam, gamma=gamma, k_actual=len(sel),
                                inner_cv_auc=auc, score=score)
    if best is None:
        return None, log
    # final SA at chosen (K, λ, γ)
    sel_final, E_final = select_genes_for_cell_type(
        bundle, candidates_ct, k=best["K"], lam=best["lam"], gamma=best["gamma"],
        seed=seed_base, score_fn=SCORE_FN, alpha_batch=ALPHA_BATCH,
        sa_reads=SA_READS, sa_sweeps=SA_SWEEPS)
    best["selected"] = sel_final
    best["energy_final"] = E_final
    return best, log


# ============================================================
# Consensus QUBO selection (stability-selection-style)
# ============================================================
def grid_search_per_ct_consensus(bundle, candidates_ct, fold, seed_base,
                                  n_runs=CONSENSUS_N_RUNS):
    """Reuses grid_search_per_ct to pick best (K, λ, γ) by inner CV AUC,
    then runs SA n_runs times with different seeds at the chosen params,
    and returns the consensus K genes (top-K by selection frequency).

    This preserves the K=K cardinality guarantee while reducing variance
    from SA's stochastic search (analogous to Meinshausen-Bühlmann
    stability selection).
    """
    best, log = grid_search_per_ct(bundle, candidates_ct, fold, seed_base)
    if best is None:
        return None, log
    K_chosen = best["K"]
    lam     = best["lam"]
    gamma   = best["gamma"]

    freq = {}            # gene -> count across runs
    energies = []
    for k in range(n_runs):
        run_seed = seed_base + 10000 + k * CONSENSUS_SEED_STEP
        sel_k, E_k = select_genes_for_cell_type(
            bundle, candidates_ct, k=K_chosen, lam=lam, gamma=gamma,
            seed=run_seed, score_fn=SCORE_FN, alpha_batch=ALPHA_BATCH,
            sa_reads=CONSENSUS_SA_READS, sa_sweeps=CONSENSUS_SA_SWEEPS)
        for g in sel_k:
            freq[g] = freq.get(g, 0) + 1
        energies.append(E_k)

    if not freq:
        return None, log

    # Rank by selection frequency, break ties by gene-level |t-stat| if available
    if bundle.get("tstats") is not None:
        ts = bundle["tstats"].set_index("gene")["t"].abs().to_dict()
    else:
        ts = {}
    ranked = sorted(freq.items(),
                    key=lambda x: (-x[1], -ts.get(x[0], 0.0)))
    consensus = [g for g, c in ranked[:K_chosen]]

    best_consensus = dict(best)
    best_consensus["selected"] = consensus
    best_consensus["n_unique_pool"] = len(freq)
    best_consensus["mean_energy_runs"] = float(np.mean(energies))
    # Store frequency annotations for downstream reporting
    best_consensus["selection_frequency"] = {g: freq[g] for g in consensus}
    return best_consensus, log


# ============================================================
# Ensemble fit & predict
# ============================================================
def fit_predict_stacking(bundles, selected_per_ct, n_inner=3, seed=42,
                         meta_factory=None):
    """
    Stacking ensemble:
      1. For each cell type, fit base classifier on full train; predict on val/held.
      2. Generate out-of-fold (OOF) predictions on TRAIN via inner CV.
      3. Stack OOF predictions into a (donors x cell_types) feature matrix.
      4. Train a meta-classifier (LR_L2 by default) on these features.
      5. Apply meta-classifier to val/held cell-type predictions.
    Returns same shape as fit_predict_ensemble.
    """
    if meta_factory is None:
        meta_factory = META_CLF_FACTORY

    p_oof_per_ct  = {}   # ct -> {donor: oof_prob}
    p_val_per_ct  = {}   # ct -> {donor: prob}
    p_held_per_ct = {}
    val_diag_map  = {}
    held_diag_map = {}
    train_diag_map = {}
    per_ct_diag    = {}

    for ct, b in bundles.items():
        if b is None:
            continue
        sel = selected_per_ct.get(ct)
        if not sel:
            continue
        tr = build_X_y_per_ct(b, sel, "train")
        va = build_X_y_per_ct(b, sel, "val")
        he = build_X_y_per_ct(b, sel, "heldout")
        if tr is None or va is None or len(tr["X"]) < 4:
            continue

        X_tr, y_tr, donors_tr = tr["X"], tr["y"], tr["donors"]
        # diag map
        for d, y in zip(donors_tr, y_tr):
            train_diag_map[d] = "MS" if y == 1 else "HD"

        # ---- inner CV for OOF ----
        rng = np.random.default_rng(seed)
        n = len(X_tr)
        pos = [i for i, y in enumerate(y_tr) if y == 1]
        neg = [i for i, y in enumerate(y_tr) if y == 0]
        rng.shuffle(pos); rng.shuffle(neg)
        fold_assign = [[] for _ in range(n_inner)]
        for i, p in enumerate(pos): fold_assign[i % n_inner].append(p)
        for i, p in enumerate(neg): fold_assign[i % n_inner].append(p)
        oof = np.full(n, np.nan)
        for k in range(n_inner):
            v_idx = fold_assign[k]
            t_idx = [i for i in range(n) if i not in v_idx]
            if len(set(y_tr[t_idx])) < 2 or len(v_idx) == 0:
                continue
            mu = X_tr[t_idx].mean(0); sd = X_tr[t_idx].std(0); sd[sd == 0] = 1.0
            Xtz = (X_tr[t_idx] - mu) / sd
            Xvz = (X_tr[v_idx] - mu) / sd
            try:
                clf = BASE_CLF_FACTORY().fit(Xtz, y_tr[t_idx])
                ps = clf.predict_proba(Xvz)
                for j, idx in enumerate(v_idx):
                    oof[idx] = float(ps[j])
            except Exception:
                continue

        oof_dict = {donors_tr[i]: float(oof[i]) for i in range(n) if not np.isnan(oof[i])}
        if not oof_dict:
            continue
        p_oof_per_ct[ct] = oof_dict

        # ---- full-train classifier for val/held ----
        mu = X_tr.mean(0); sd = X_tr.std(0); sd[sd == 0] = 1.0
        Xtz = (X_tr - mu) / sd
        Xvz = (va["X"] - mu) / sd
        clf = BASE_CLF_FACTORY().fit(Xtz, y_tr)
        p_val = clf.predict_proba(Xvz)
        for d, p, y in zip(va["donors"], p_val, va["y"]):
            p_val_per_ct.setdefault(ct, {})[d] = float(p)
            val_diag_map[d] = "MS" if y == 1 else "HD"
        if he is not None and len(he["X"]) > 0:
            Xhz = (he["X"] - mu) / sd
            p_held = clf.predict_proba(Xhz)
            for d, p, y in zip(he["donors"], p_held, he["y"]):
                p_held_per_ct.setdefault(ct, {})[d] = float(p)
                held_diag_map[d] = "MS" if y == 1 else "HD"

        per_ct_diag[ct] = {
            "n_train": len(donors_tr),
            "n_val": len(va["donors"]),
            "val_auc": roc_auc(va["y"], p_val),
            "n_oof": int(np.sum(~np.isnan(oof))),
        }

    if not p_oof_per_ct:
        return [], [], per_ct_diag

    cell_types_used = sorted(p_oof_per_ct.keys())

    # ---- build training meta-feature matrix ----
    train_donors_list = sorted({d for ct in cell_types_used for d in p_oof_per_ct[ct]})
    M_train = np.zeros((len(train_donors_list), len(cell_types_used)))
    for j, ct in enumerate(cell_types_used):
        col_vals = [p_oof_per_ct[ct].get(d, np.nan) for d in train_donors_list]
        col_vals = np.array(col_vals)
        # mean-impute NaN
        if np.any(np.isnan(col_vals)):
            valid = ~np.isnan(col_vals)
            if valid.any():
                col_vals[~valid] = col_vals[valid].mean()
            else:
                col_vals[:] = 0.5
        M_train[:, j] = col_vals
    y_train = np.array([1 if train_diag_map.get(d) == "MS" else 0
                        for d in train_donors_list])
    if len(np.unique(y_train)) < 2:
        return [], [], per_ct_diag

    # ---- train meta-classifier on OOF features ----
    mu_m = M_train.mean(0); sd_m = M_train.std(0); sd_m[sd_m == 0] = 1.0
    M_train_z = (M_train - mu_m) / sd_m
    meta_clf = meta_factory().fit(M_train_z, y_train)

    # ---- apply meta to val ----
    val_donors_list = sorted({d for ct in cell_types_used
                               for d in p_val_per_ct.get(ct, {})})
    M_val = np.zeros((len(val_donors_list), len(cell_types_used)))
    for j, ct in enumerate(cell_types_used):
        col_vals = [p_val_per_ct.get(ct, {}).get(d, np.nan) for d in val_donors_list]
        col_vals = np.array(col_vals)
        if np.any(np.isnan(col_vals)):
            valid = ~np.isnan(col_vals)
            col_vals[~valid] = col_vals[valid].mean() if valid.any() else 0.5
        M_val[:, j] = col_vals
    M_val_z = (M_val - mu_m) / sd_m
    p_val_final = meta_clf.predict_proba(M_val_z)
    val_rows = [{"donor": d, "diagnosis": val_diag_map.get(d),
                 "prob_MS": float(p), "n_cell_types_used": len(cell_types_used)}
                for d, p in zip(val_donors_list, p_val_final)]

    # ---- apply meta to held ----
    held_donors_list = sorted({d for ct in cell_types_used
                                for d in p_held_per_ct.get(ct, {})})
    if held_donors_list:
        M_held = np.zeros((len(held_donors_list), len(cell_types_used)))
        for j, ct in enumerate(cell_types_used):
            col_vals = [p_held_per_ct.get(ct, {}).get(d, np.nan) for d in held_donors_list]
            col_vals = np.array(col_vals)
            if np.any(np.isnan(col_vals)):
                valid = ~np.isnan(col_vals)
                col_vals[~valid] = col_vals[valid].mean() if valid.any() else 0.5
            M_held[:, j] = col_vals
        M_held_z = (M_held - mu_m) / sd_m
        p_held_final = meta_clf.predict_proba(M_held_z)
        held_rows = [{"donor": d, "diagnosis": held_diag_map.get(d),
                      "prob_MS": float(p), "n_cell_types_used": len(cell_types_used)}
                     for d, p in zip(held_donors_list, p_held_final)]
    else:
        held_rows = []

    return val_rows, held_rows, per_ct_diag


def fit_predict_ensemble(bundles, selected_per_ct):
    """For each cell type, fit a classifier on train pseudobulk × selected genes,
    predict on val and heldout. Then ensemble (mean) the probabilities.

    Returns:
      val_donors_y_p: list of (donor_id, diagnosis, prob_MS, used_n_ct)
      held_donors_y_p: same for heldout
      per_ct_aucs: dict[ct] -> dict(val_auc, held_auc) for diagnostics
    """
    # Per cell type predictions (mapping donor_id -> prob)
    val_probs_per_ct = {}   # ct -> {donor_id: prob}
    held_probs_per_ct = {}
    val_diag_map = {}
    held_diag_map = {}
    per_ct_diag = {}

    for ct, b in bundles.items():
        if b is None:
            continue
        sel = selected_per_ct.get(ct)
        if not sel:
            continue
        tr = build_X_y_per_ct(b, sel, "train")
        va = build_X_y_per_ct(b, sel, "val")
        he = build_X_y_per_ct(b, sel, "heldout")
        if tr is None or va is None or len(tr["X"]) < 4:
            continue
        # standardize using train
        mu = tr["X"].mean(0); sd = tr["X"].std(0); sd[sd == 0] = 1.0
        Xtz = (tr["X"] - mu) / sd
        Xvz = (va["X"] - mu) / sd
        clf = BASE_CLF_FACTORY().fit(Xtz, tr["y"])
        p_val = clf.predict_proba(Xvz)
        for d, p, y in zip(va["donors"], p_val, va["y"]):
            val_probs_per_ct.setdefault(ct, {})[d] = float(p)
            val_diag_map[d] = "MS" if y == 1 else "HD"
        per_ct_diag[ct] = {
            "n_val": len(va["donors"]),
            "val_auc": roc_auc(va["y"], p_val),
            "selected_genes": sel,
        }
        if he is not None and len(he["X"]) > 0:
            Xhz = (he["X"] - mu) / sd
            p_held = clf.predict_proba(Xhz)
            for d, p, y in zip(he["donors"], p_held, he["y"]):
                held_probs_per_ct.setdefault(ct, {})[d] = float(p)
                held_diag_map[d] = "MS" if y == 1 else "HD"
            per_ct_diag[ct]["n_held"] = len(he["donors"])
            per_ct_diag[ct]["held_auc"] = roc_auc(he["y"], p_held)

    # Aggregate: for each donor, mean of cell-type predictions where present
    def aggregate(probs_per_ct, diag_map):
        all_donors = sorted(set(d for v in probs_per_ct.values() for d in v))
        rows = []
        for d in all_donors:
            ps = [probs_per_ct[ct][d] for ct in probs_per_ct if d in probs_per_ct[ct]]
            if not ps:
                continue
            p_avg = float(np.mean(ps))
            rows.append({"donor": d, "diagnosis": diag_map.get(d),
                         "prob_MS": p_avg, "n_cell_types_used": len(ps)})
        return rows

    val_rows = aggregate(val_probs_per_ct, val_diag_map)
    held_rows = aggregate(held_probs_per_ct, held_diag_map)
    # v9: also return per-cell-type per-donor probabilities so the caller can
    # implement post-hoc AUC-weighted aggregation. Format: list of dicts.
    val_per_ct_rows = []
    for ct, dprobs in val_probs_per_ct.items():
        for donor, p in dprobs.items():
            val_per_ct_rows.append({
                "donor": donor, "cell_type": ct,
                "diagnosis": val_diag_map.get(donor),
                "prob_MS": float(p),
            })
    held_per_ct_rows = []
    for ct, dprobs in held_probs_per_ct.items():
        for donor, p in dprobs.items():
            held_per_ct_rows.append({
                "donor": donor, "cell_type": ct,
                "diagnosis": held_diag_map.get(donor),
                "prob_MS": float(p),
            })
    return val_rows, held_rows, per_ct_diag, val_per_ct_rows, held_per_ct_rows


# ============================================================
# Main per-(holdout, deg, tissue, fold) run
# ============================================================
def run_for_tissue(holdout_name, deg_source, tissue, folds):
    data_root = _data_root(holdout_name)
    tag = RUN_TAG
    if BIOLOGY_FILTER:
        tag = f"{tag}_bio"
    tag = f"{tag}_{deg_source}"
    if holdout_name != "Pappalardo":
        tag = f"{tag}_holdout_{holdout_name}"
    out_dir = PROJECT_ROOT / "qubo_run_v6" / tag / tissue
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n========== TISSUE = {tissue}  TAG = {tag}  DEG = {deg_source} ==========")

    fold_records = []
    val_pred_records = []
    held_pred_records = []
    val_per_ct_records = []   # v9: per-cell-type per-donor probabilities
    held_per_ct_records = []
    selected_log = []
    grid_log = []
    per_ct_records = []

    for fold in folds:
        # Load bundles for all cell types
        bundles = {}
        any_present = False
        for ct in CELL_TYPES:
            b = load_fold(data_root, ct, tissue, fold,
                          aggregator="mean", deg_source=deg_source)
            if b is not None and b.get("train") is not None:
                bundles[ct] = b
                any_present = True
        if not any_present:
            print(f"  [skip] {tissue}/fold_{fold}: no data")
            continue

        # Per-cell-type grid search for ALL METHODS
        # selected_per_method[method][ct] = list of genes
        selected_per_method = {m: {} for m in METHODS}
        K_per_method = {m: {} for m in METHODS}

        for ct, b in bundles.items():
            cands = candidates_per_cell_type(b, N_PER_CELL_TYPE)
            if len(cands) < 5:
                continue
            seed_ct = SEED + fold * 100 + _ct_seed(ct)

            # QUBO (uses (K, λ, γ) grid)
            if "QUBO" in METHODS:
                best, gl = grid_search_per_ct(b, cands, fold, seed_base=seed_ct)
                for entry in gl:
                    grid_log.append({"tissue": tissue, "fold": fold,
                                     "cell_type": ct, "method": "QUBO", **entry})
                if best is not None:
                    selected_per_method["QUBO"][ct] = best["selected"]
                    K_per_method["QUBO"][ct] = best["K"]
                    for g in best["selected"]:
                        selected_log.append({"tissue": tissue, "fold": fold,
                                             "cell_type": ct, "method": "QUBO",
                                             "K_chosen": best["K"], "gene": g,
                                             "inner_cv_auc": best["inner_cv_auc"]})

            # QUBO_consensus: same K-grid but final selection = top-K by
            # frequency over CONSENSUS_N_RUNS independent SA runs.
            if "QUBO_consensus" in METHODS:
                best, gl = grid_search_per_ct_consensus(b, cands, fold,
                                                         seed_base=seed_ct)
                for entry in gl:
                    grid_log.append({"tissue": tissue, "fold": fold,
                                     "cell_type": ct, "method": "QUBO_consensus",
                                     **entry})
                if best is not None:
                    selected_per_method["QUBO_consensus"][ct] = best["selected"]
                    K_per_method["QUBO_consensus"][ct] = best["K"]
                    freq_map = best.get("selection_frequency", {})
                    for g in best["selected"]:
                        selected_log.append({"tissue": tissue, "fold": fold,
                                             "cell_type": ct, "method": "QUBO_consensus",
                                             "K_chosen": best["K"], "gene": g,
                                             "inner_cv_auc": best["inner_cv_auc"],
                                             "selection_frequency": freq_map.get(g, 0)})

            # QUBO_hybrid: pre-filter to top HYBRID_TOP_N by |t|, then QUBO grid
            if "QUBO_hybrid" in METHODS:
                cands_hybrid = candidates_per_cell_type_hybrid(b, n_top=HYBRID_TOP_N)
                if len(cands_hybrid) >= 5:
                    best, gl = grid_search_per_ct(b, cands_hybrid, fold,
                                                  seed_base=seed_ct + 7)
                    for entry in gl:
                        grid_log.append({"tissue": tissue, "fold": fold,
                                         "cell_type": ct, "method": "QUBO_hybrid",
                                         "n_pool": len(cands_hybrid), **entry})
                    if best is not None:
                        selected_per_method["QUBO_hybrid"][ct] = best["selected"]
                        K_per_method["QUBO_hybrid"][ct] = best["K"]
                        for g in best["selected"]:
                            selected_log.append({"tissue": tissue, "fold": fold,
                                                 "cell_type": ct, "method": "QUBO_hybrid",
                                                 "K_chosen": best["K"], "gene": g,
                                                 "inner_cv_auc": best["inner_cv_auc"]})

            # Baselines (K grid only, no λ,γ) — NOT QUBO variants
            for m_name in [m for m in METHODS if m not in ("QUBO", "QUBO_hybrid", "QUBO_consensus")]:
                best, gl = grid_search_baseline_per_ct(m_name, b, cands, fold,
                                                        seed_base=seed_ct)
                for entry in gl:
                    grid_log.append({"tissue": tissue, "fold": fold,
                                     "cell_type": ct, **entry})
                if best is not None:
                    selected_per_method[m_name][ct] = best["selected"]
                    K_per_method[m_name][ct] = best["K"]
                    for g in best["selected"]:
                        selected_log.append({"tissue": tissue, "fold": fold,
                                             "cell_type": ct, "method": m_name,
                                             "K_chosen": best["K"], "gene": g,
                                             "inner_cv_auc": best["inner_cv_auc"]})

        # For each method, run ensemble
        for m_name in METHODS:
            sel = selected_per_method[m_name]
            if not sel:
                continue
            if ENSEMBLE_AGG == "stacking":
                ens_out = fit_predict_stacking(bundles, sel)
                # stacking still returns 3-tuple; pad with empty per-ct lists
                if len(ens_out) == 3:
                    val_rows, held_rows, per_ct = ens_out
                    val_per_ct_rows, held_per_ct_rows = [], []
                else:
                    val_rows, held_rows, per_ct, val_per_ct_rows, held_per_ct_rows = ens_out
            else:
                val_rows, held_rows, per_ct, val_per_ct_rows, held_per_ct_rows = \
                    fit_predict_ensemble(bundles, sel)
            if not val_rows:
                continue
            # save predictions for ALL methods (per-method per-cell-type predictions
            # enable post-hoc AUC-weighted aggregation)
            for r in val_rows:
                r2 = {**r, "tissue": tissue, "fold": fold,
                      "method": m_name, "set": "val"}
                val_pred_records.append(r2)
            for r in held_rows:
                r2 = {**r, "tissue": tissue, "fold": fold,
                      "method": m_name, "set": "heldout"}
                held_pred_records.append(r2)
            # v9: per-cell-type per-donor probabilities for post-hoc aggregation
            for r in val_per_ct_rows:
                r2 = {**r, "tissue": tissue, "fold": fold,
                      "method": m_name, "set": "val"}
                val_per_ct_records.append(r2)
            for r in held_per_ct_rows:
                r2 = {**r, "tissue": tissue, "fold": fold,
                      "method": m_name, "set": "heldout"}
                held_per_ct_records.append(r2)

            y_v = np.array([1 if r["diagnosis"] == "MS" else 0 for r in val_rows])
            p_v = np.array([r["prob_MS"] for r in val_rows])
            y_h = np.array([1 if r["diagnosis"] == "MS" else 0 for r in held_rows]) if held_rows else np.array([])
            p_h = np.array([r["prob_MS"] for r in held_rows]) if held_rows else np.array([])

            m = dict(tissue=tissue, fold=fold, method=m_name,
                     n_cell_types=len(sel),
                     val_auc=roc_auc(y_v, p_v),
                     val_ap=average_precision(y_v, p_v),
                     val_n=len(y_v),
                     held_auc=roc_auc(y_h, p_h) if len(y_h) else np.nan,
                     held_ap=average_precision(y_h, p_h) if len(y_h) else np.nan,
                     held_n=len(y_h))
            a, f = acc_f1(y_v, p_v); m["val_acc"] = a; m["val_f1"] = f
            m["val_mcc"] = mcc_score(y_v, p_v)
            if len(y_h):
                a, f = acc_f1(y_h, p_h); m["held_acc"] = a; m["held_f1"] = f
                m["held_mcc"] = mcc_score(y_h, p_h)
            else:
                m["held_acc"] = np.nan; m["held_f1"] = np.nan
                m["held_mcc"] = np.nan
            fold_records.append(m)

        # per-cell-type diagnostic for the primary QUBO variant in METHODS
        diag_method = "QUBO" if "QUBO" in METHODS else (
            "QUBO_consensus" if "QUBO_consensus" in METHODS else None)
        if diag_method and selected_per_method.get(diag_method):
            _, _, per_ct, _, _ = fit_predict_ensemble(bundles, selected_per_method[diag_method])
            for ct, d in per_ct.items():
                rec = {"tissue": tissue, "fold": fold, "cell_type": ct,
                       "diag_method": diag_method,
                       "n_genes": len(selected_per_method[diag_method].get(ct, [])),
                       "K_chosen": K_per_method[diag_method].get(ct)}
                rec.update({k: v for k, v in d.items() if k != "selected_genes"})
                per_ct_records.append(rec)

        # print fold summary
        line = f"  fold {fold}: " + " | ".join(
            f"{m}={(next((r for r in fold_records if r.get('method')==m and r['fold']==fold and r['tissue']==tissue), {}).get('val_auc', float('nan'))):.2f}/{(next((r for r in fold_records if r.get('method')==m and r['fold']==fold and r['tissue']==tissue), {}).get('held_auc', float('nan'))):.2f}"
            for m in METHODS if any(r['method']==m and r['fold']==fold and r['tissue']==tissue for r in fold_records)
        )
        print(line)

    # write outputs
    suf = f"_folds_{'_'.join(str(f) for f in folds)}"
    if fold_records:
        pd.DataFrame(fold_records).to_csv(out_dir / f"fold_metrics{suf}.csv", index=False)
    if val_pred_records:
        pd.DataFrame(val_pred_records).to_csv(out_dir / f"val_predictions{suf}.csv", index=False)
    if held_pred_records:
        pd.DataFrame(held_pred_records).to_csv(out_dir / f"held_predictions{suf}.csv", index=False)
    if val_per_ct_records:
        pd.DataFrame(val_per_ct_records).to_csv(out_dir / f"val_per_ct_predictions{suf}.csv", index=False)
    if held_per_ct_records:
        pd.DataFrame(held_per_ct_records).to_csv(out_dir / f"held_per_ct_predictions{suf}.csv", index=False)
    if selected_log:
        pd.DataFrame(selected_log).to_csv(out_dir / f"selected_genes{suf}.csv", index=False)
    if grid_log:
        pd.DataFrame(grid_log).to_csv(out_dir / f"grid_log{suf}.csv", index=False)
    if per_ct_records:
        pd.DataFrame(per_ct_records).to_csv(out_dir / f"per_ct_diag{suf}.csv", index=False)


def main():
    if len(sys.argv) < 5:
        print("Usage: python3 qubo_pipeline_v6.py <holdout> <deg> <tissue> <fold1> [fold2 ...]")
        sys.exit(1)
    global HOLDOUT_NAME, DEG_SOURCE
    HOLDOUT_NAME = sys.argv[1]
    DEG_SOURCE  = sys.argv[2]
    tissue = sys.argv[3]
    folds = [int(x) for x in sys.argv[4:]]
    run_for_tissue(HOLDOUT_NAME, DEG_SOURCE, tissue, folds)


if __name__ == "__main__":
    main()
