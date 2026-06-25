"""
qubo_pipeline.py
================
Per-cell-type QUBO feature selection with a fixed matched panel size for the
primary benchmark (and optional panel-size sensitivity analyses), with a
soft-voting classifier ensemble across eight broad immune cell types.

For each (leave-one-cohort-out hold-out, cross-validation fold, cell type):
  - QUBO is solved independently per cell type; each cell type receives its
    own gene panel.
  - In the primary benchmark, K is fixed by QUBOFS_FIXED_K=10 for all methods;
    the QUBO parameters lambda and gamma are selected by inner cross-validation.
  - An L2 logistic-regression classifier is trained per cell type.
  - The donor-level prediction is the unweighted soft vote (mean) over
    cell-type probabilities.

Usage:
    python3 qubo_pipeline.py <holdout> <deg> <tissue> <fold1> [fold2 ...]
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
# Pre-specified gene-level TECHNICAL filter (manuscript §2.2/§2.4)
# ------------------------------------------------------------
# A minimal, pre-specified filter applied identically to ALL methods at the
# candidate stage. It removes only genes likely to dominate technical,
# non-specific or clonotype-driven variation (mito / ribosomal / cytosolic
# stress HSP / dominant lncRNA / housekeeping; and Ig/TCR V,D,J segments via
# is_vdj_segment). It is NOT a strong biological filter: Ig/TCR constant
# regions and ER chaperone / secretion genes are retained by design.
# (Variable/function names keep the legacy "biology" label for stability.)
# ============================================================
HK_PATTERN = re.compile(
    # =========================================================================
    #  Pre-selection biology filter (manuscript primary configuration).
    #  Aligned with current best practice (Heumos et al. Nat Rev Genet 2023;
    #  Luecken & Theis Mol Syst Biol 2019).
    #
    #  All ribosomal-protein genes are excluded, including the acidic stalk
    #  paralogs (RPLP0/1/2) and the laminin-receptor RPSA, so that no
    #  ribosomal-protein-named gene can enter a biomarker panel (consistent with
    #  the "ribosomal genes excluded" statement in Supplementary Table S6).
    # =========================================================================
    r"^(MT-|MTRNR|MTATP|MTND|"             # mitochondrial
    r"RPL[0-9]|RPLP|RPS[0-9]|RPSA|MRPL|MRPS|"  # ribosomal proteins (incl. RPLP0/1/2 stalk and RPSA)
    r"HSP[A0-9]|HSPB|HSPA|HSPD|"           # heat shock (HSPA/B/D + HSP9*)
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|"  # classical housekeeping
    r"MALAT1$|NEAT1$|XIST$|TSIX$|"         # nuclear lncRNA + X-inactivation
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|"    # uncharacterized / pseudogene loci
    r"MIR[0-9]|RNU[0-9]|SNORA|SNORD)"      # small RNAs (poly-A unreliable)
)
# ER chaperone / secretion genes retained by design: the heat-shock pattern would
# otherwise drop these, but they are ER-chaperone / secretory-pathway-associated
# features and should not be removed as generic cytosolic stress genes.
#   HSPA5  = BiP / GRP78  ;  HSP90B1 = GRP94 / endoplasmin
BIOLOGY_RETAIN = {"HSPA5", "HSP90B1"}
def is_biology_gene(g: str) -> bool:
    """True if gene passes the biology filter (drops technical HK/mito/ribosomal).

    ER chaperones in BIOLOGY_RETAIN are kept despite matching the heat-shock
    pattern. Immunoglobulin/TCR constant-region genes are not matched here and are
    retained; clonotype-driven V/D/J segments are excluded separately by
    is_vdj_segment in the cell-type-aware filter.
    """
    g = str(g)
    if g in BIOLOGY_RETAIN:
        return True
    return not bool(HK_PATTERN.match(g))


_VDJ_PREFIX = ("IGHV", "IGKV", "IGLV", "TRBV", "TRAV", "TRGV", "TRDV",
               "TRBJ", "TRAJ", "TRGJ", "TRDJ", "IGHJ", "IGKJ", "IGLJ")
def is_vdj_segment(g: str) -> bool:
    """True for clonotype-driven Ig/TCR V/D/J segment genes (constant regions kept)."""
    g = str(g)
    if g.startswith(_VDJ_PREFIX):
        return True
    if re.match(r"^IGHD\d", g):   # IGHD diversity segments (keep IGHD constant)
        return True
    return False


def compute_allowed_genes(bundles):
    """Cell-type-aware candidate filter (manuscript §2.4), shared by ALL methods.

    Returns {cell_type: set(genes passing)} where a gene passes for a target cell
    type if: detected in >= DET_THR of training donors for that cell type;
    specificity ratio (target-cell-type mean / max mean across the other cell
    types) >= SPEC_THR; and it is not a V(D)J segment gene. Means/detection are
    computed on the training-donor pseudobulk.
    """
    # per-(ct) detection rate and mean expression on the training pseudobulk
    det = {}
    mean_by_ct = {}
    for ct, b in bundles.items():
        tr = b.get("train") if b else None
        if tr is None:
            continue
        X = tr["X"]
        gp = {g: i for i, g in enumerate(tr["genes"])}
        det[ct] = (gp, (X > 0).mean(axis=1))
        mean_by_ct[ct] = (gp, X.mean(axis=1))
    cts = list(mean_by_ct.keys())
    allowed = {}
    for ct in cts:
        gp_d, dvec = det[ct]
        gp_m, mvec = mean_by_ct[ct]
        passing = set()
        for g, i in gp_d.items():
            if EXCLUDE_VDJ and is_vdj_segment(g):
                continue
            if dvec[i] < DET_THR:
                continue
            tmean = mvec[gp_m[g]] if g in gp_m else 0.0
            if tmean <= 0:
                continue
            others = [mean_by_ct[oc][1][mean_by_ct[oc][0][g]]
                      for oc in cts if oc != ct and g in mean_by_ct[oc][0]]
            omax = max(others) if others else 0.0
            ratio = (tmean / omax) if omax > 0 else float("inf")
            if ratio >= SPEC_THR:
                passing.add(g)
        allowed[ct] = passing
    return allowed

sys.path.insert(0, str(Path(__file__).parent))
from qubo_utils import (
    load_fold, build_score_and_redundancy, build_qubo, solve_qubo_sa,
    cohort_variance_per_gene, cohort_consistency_per_gene,
    LogRegL2, LogRegL1, LogRegElasticNet, LDA, standardize,
    roc_auc, average_precision, acc_f1, mcc_score, jaccard,
)

# ============================================================
# Configuration
# ============================================================
# Project root is taken from the QUBOFS_PROJECT_ROOT environment variable if
# set, otherwise the script's parent-parent directory (github_release/).
PROJECT_ROOT = Path(
    os.environ.get(
        "QUBOFS_PROJECT_ROOT",
        Path(__file__).resolve().parent.parent,
    )
)
PSEUDOBULK_SUBDIR = os.environ.get("QUBOFS_PSEUDOBULK_SUBDIR", "pseudobulk_v5_compartment")

HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}

# Eight broad immune cell types (manuscript primary configuration).
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
TISSUES = ["CSF"]
FOLDS = [1, 2, 3, 4, 5]

# K grid (per-cell-type cardinality). The primary benchmark uses a single fixed
# panel size K = 10 for all methods (set via QUBOFS_FIXED_K, default 10), so
# running this script without extra environment variables reproduces the
# manuscript's matched K = 10 configuration. The panel-size sensitivity sweep
# (Supplementary Figure S2) sets QUBOFS_K_SWEEP=<K> to evaluate one alternative K.
if os.environ.get("QUBOFS_K_SWEEP"):
    K_GRID = [int(os.environ["QUBOFS_K_SWEEP"])]
else:
    K_GRID = [int(os.environ.get("QUBOFS_FIXED_K", "10"))]
LAMBDA_VALS = [2.0, 5.0]        # manuscript: λ ∈ {2, 5}
GAMMA_VALS  = [0.5, 1.0]        # manuscript: γ ∈ {0.5, 1.0} (tuned by inner CV).
                                # Expanding the grid (e.g. adding 0.25) was tested and
                                # did not improve held-out performance on this small
                                # dataset, so the manuscript grid is retained.
SA_READS_GRID = 8
SA_SWEEPS_GRID = 200
SA_READS = 30
SA_SWEEPS = 600
SCORE_FN = "abs_t"      # manuscript relevance uses |z_i| (absolute edgeR test statistic), weighted by C_i
ALPHA_BATCH = 1.0
INNER_CV_FOLDS = 5             # manuscript: inner five-fold donor-stratified CV
N_INNER_REPEATS = int(os.environ.get("QUBOFS_INNER_REPEATS", "10"))  # repeated inner CV (different shuffles) to stabilise
                              # hyperparameter selection on small cohorts; applied
                              # identically to QUBO and all baselines
SELECTION_METRIC = "auc"      # inner-CV model-selection metric, applied to ALL
                              # methods: "auc" | "mcc" | "macro_f1". Default "auc";
                              # changing this is a design choice that must be applied
                              # uniformly and reported (see CHANGELOG).
DET_THR = 0.70          # manuscript: gene detected in >= 70% of training donors for the cell type
SPEC_THR = 0.70         # manuscript: cell-type specificity ratio (target mean / max other-cell-type mean)
EXCLUDE_VDJ = (os.environ.get("QUBOFS_EXCLUDE_VDJ", "1") == "1")  # exclude clonotype-driven Ig/TCR V/D/J segment genes (constant regions kept). Configuration toggle; default True (primary).
USE_COHORT_CONSISTENCY = (os.environ.get("QUBOFS_USE_COHORT_CONSISTENCY", "1") == "1")  # weight relevance by directional cohort-consistency C_i (|z|*C). Configuration toggle; default True (primary).
N_PER_CELL_TYPE = 100   # candidate set per cell type (top-N from t-statistics); baselines use this pool
N_QUBO_SCREEN = int(os.environ.get("QUBOFS_QUBO_SCREEN", "20"))  # manuscript: QUBO uses a top-20 sure-independence-screening pre-filter
                        # (Fan and Lv 2008) before optimisation; baselines use the full top-100 pool.
                        # Env-configurable for the screen-size sensitivity analysis.
# Optional inner-CV selection of the screen size. If QUBOFS_SCREEN_GRID lists
# several sizes (e.g. "20,25,30,35,40,45,50"), the screen is chosen per
# (cell type, fold) by inner-CV AUC, alongside (K, λ, γ); held-out cohorts are
# never used. Default = the single fixed N_QUBO_SCREEN (unchanged behaviour).
_screen_grid_env = os.environ.get("QUBOFS_SCREEN_GRID", "").strip()
SCREEN_GRID = ([int(x) for x in _screen_grid_env.split(",") if x.strip()]
               if _screen_grid_env else [N_QUBO_SCREEN])
# Optional hybrid (classifier-informed) relevance for QUBO only:
#   s_i = (1-eta)*minmax(|z_i|C_i) + eta*minmax(|w_i|)
# where w_i are L2-logistic coefficients fit on the TRAINING candidate pool
# (quboFS's own downstream classifier; no held-out leakage). eta is selected by
# inner-CV from QUBOFS_ETA_GRID (default "0.0" = pure DE relevance, unchanged).
_eta_grid_env = os.environ.get("QUBOFS_ETA_GRID", "").strip()
ETA_GRID = ([float(x) for x in _eta_grid_env.split(",") if x.strip()]
            if _eta_grid_env else [0.0])
# Optional two-tier "rescue" screen for QUBO: union the top-N_QUBO_SCREEN genes by
# edgeR |z|C with the top-QUBOFS_RESCUE_N genes by multivariate (Elastic Net)
# coefficient over the full candidate pool, so classifier-informative genes that
# fall below the univariate screen can still be selected. Training only.
# Default 0 = off (univariate screen only; unchanged behaviour).
RESCUE_N = int(os.environ.get("QUBOFS_RESCUE_N", "0"))
# Optional relevance-dependent redundancy penalty: down-weight the pairwise
# redundancy |rho_ij| when the pair is highly relevant, so correlated-but-
# discriminative genes are easier to retain. R'_ij = |rho_ij| * (1 - max(s_i,s_j)),
# s in [0,1]. Default off (uniform redundancy penalty; unchanged behaviour).
RELDEP_REDUNDANCY = (os.environ.get("QUBOFS_RELDEP_REDUNDANCY", "0") == "1")
SEED = 42

# Classifier per cell type (single classifier choice for ensemble simplicity)
BASE_CLF_FACTORY = lambda: LogRegL2(C=1.0, max_iter=200)
TUNE_CLF_C = False           # tested (C tuned by inner CV, all methods) but it
                             # OVERFIT the noisy small-cohort inner CV and reduced
                             # held-out AUC; fixed C=1 generalises better here.
CLF_C_GRID = [0.01, 0.1, 1.0, 10.0]
ENSEMBLE_AGG = os.environ.get("QUBOFS_ENSEMBLE_AGG", "mean")  # unweighted soft vote (default). "weighted_mcc" (weights ∝ inner-CV
                             # MCC) was tested and also reduced held-out performance on
                             # this small dataset. "stacking" also available. Configuration toggle.
META_CLF_FACTORY = lambda: LogRegL2(C=0.5, max_iter=200)   # for stacking

RUN_TAG = os.environ.get("QUBOFS_PIPELINE_RUN_TAG", "primary")  # base tag; "_bio_<deg>" appended in run_for_tissue. Override for the K-sweep (e.g. "sweepK5").
HOLDOUT_NAME = "Pappalardo"
DEG_SOURCE = "edger_counts"  # primary relevance: edgeR empirical-Bayes test statistic on donor-level SoupX-corrected count sums (rounded for the count model)
BIOLOGY_FILTER = True    # manuscript primary config (tag "primary_bio_*"); drop housekeeping/mito/ribosomal at candidate stage

# Methods evaluated in the manuscript benchmark. Env-restrictable via QUBOFS_METHODS
# (e.g. to run a single method); default = full set.
METHODS = ["QUBO", "DE_top", "mRMR", "HVG", "LASSO", "ElasticNet"]
_M = os.environ.get("QUBOFS_METHODS")
if _M:
    METHODS = [m for m in _M.split(",") if m]


def _data_root(holdout_name: str) -> Path:
    if holdout_name == "Pappalardo":
        return PROJECT_ROOT / "data" / PSEUDOBULK_SUBDIR
    return PROJECT_ROOT / "data" / f"{PSEUDOBULK_SUBDIR}_holdout_{HOLDOUT_PRJ_MAP[holdout_name]}"


# ============================================================
# Per-cell-type QUBO
# ============================================================
def multivariate_relevance(bundle, genes):
    """Classifier-informed relevance: |coefficient| of an L2-logistic regression
    fit on the TRAINING candidate pool (the same downstream classifier quboFS
    uses). Training data only; returns a vector aligned to `genes` (0 for genes
    absent from the train split)."""
    genes = list(genes)
    out = np.zeros(len(genes), dtype=float)
    info = build_X_y_per_ct(bundle, genes, "train")
    if info is None or len(info["X"]) < 4 or len(set(info["y"].tolist())) < 2:
        return out
    Xz = standardize(info["X"])[0]
    try:
        clf = LogRegL2(C=1.0, max_iter=200).fit(Xz, info["y"])
        coef = np.abs(np.asarray(clf.coef_, dtype=float))
    except Exception:
        return out
    pos = {g: i for i, g in enumerate(info["genes"])}
    for i, g in enumerate(genes):
        if g in pos:
            out[i] = coef[pos[g]]
    return out


def rescue_candidates(bundle, cands, n_rescue):
    """Top-`n_rescue` genes from the full candidate pool by |Elastic Net coefficient|
    (multivariate, training only) — a targeted rescue of classifier-informative
    genes that fall below the univariate-relevance screen."""
    if n_rescue <= 0:
        return []
    info = build_X_y_per_ct(bundle, list(cands), "train")
    if info is None or len(info["X"]) < 4 or len(set(info["y"].tolist())) < 2:
        return []
    Xz = standardize(info["X"])[0]
    try:
        clf = LogRegElasticNet(C=1.0, l1_ratio=0.5, max_iter=200).fit(Xz, info["y"])
        coef = np.abs(np.asarray(clf.coef_, dtype=float))
    except Exception:
        return []
    order = np.argsort(-coef)
    genes = info["genes"]
    return [genes[i] for i in order[:n_rescue]]


def select_genes_for_cell_type(bundle, candidates_ct, k, lam, gamma,
                                seed, score_fn, alpha_batch,
                                sa_reads, sa_sweeps, eta=0.0):
    """Run QUBO on a single cell type's bundle and candidate list.

    `eta` (in [0,1]) blends the manuscript DE relevance with a
    classifier-informed multivariate relevance: eta=0 is the primary
    (DE-only) configuration."""
    one_bundle = {bundle["cell_type"]: bundle}
    s_raw, R, _ = build_score_and_redundancy(
        one_bundle, candidates_ct,
        score_agg="sum", redundancy_agg="max", score_fn=score_fn)

    def _mm(v):
        if v.max() > v.min():
            return (v - v.min()) / (v.max() - v.min())
        return np.zeros_like(v)

    # Manuscript relevance: |z_i| * C_i (cohort-consistency weighting), then
    # min-max rescaled to [0, 1]. C_i = directional agreement across training
    # cohorts (NOT the legacy cohort-variance penalty).
    C = cohort_consistency_per_gene(one_bundle, candidates_ct) if USE_COHORT_CONSISTENCY \
        else np.ones(len(candidates_ct))
    s = _mm(s_raw * C)
    if eta and eta > 0.0:
        m = _mm(multivariate_relevance(bundle, candidates_ct))
        s = _mm((1.0 - eta) * s + eta * m)

    if RELDEP_REDUNDANCY:
        # weaken redundancy penalty for high-relevance pairs (keep correlated-
        # but-discriminative genes): R'_ij = |rho_ij| * (1 - max(s_i, s_j)).
        R = np.asarray(R, dtype=float) * (1.0 - np.maximum.outer(np.asarray(s, float),
                                                                 np.asarray(s, float)))

    Q = build_qubo(s, R, k=k, lam=lam, gamma=gamma)
    rng = np.random.default_rng(seed)
    x, E = solve_qubo_sa(Q, k=k, n_reads=sa_reads, n_sweeps=sa_sweeps, rng=rng)
    selected_idx = np.where(x == 1)[0]
    selected = [candidates_ct[i] for i in selected_idx]
    return selected, E


def weighted_relevance(bundle, genes):
    """Manuscript relevance score |z_i| * C_i over `genes` (single cell type).

    |z_i| = absolute edgeR test statistic from the fold's tstats; C_i =
    directional cohort-consistency. This is the shared ranking quantity used by
    the candidate pool and all supervised methods (QUBO, DE-top, mRMR).
    """
    genes = list(genes)
    ts = bundle["tstats"].set_index("gene") if bundle.get("tstats") is not None else None
    abs_t = np.array([abs(float(ts.loc[g, "t"])) if (ts is not None and g in ts.index) else 0.0
                      for g in genes], dtype=float)
    if not USE_COHORT_CONSISTENCY:
        return abs_t
    C = cohort_consistency_per_gene({bundle.get("cell_type", "ct"): bundle}, genes)
    return abs_t * np.asarray(C, dtype=float)


def candidates_per_cell_type(bundle, n_top=N_PER_CELL_TYPE):
    """Top-N genes per cell type, ranked by the manuscript relevance |z_i|*C_i.
    If BIOLOGY_FILTER is True, drop housekeeping/mito/ribosomal genes BEFORE
    ranking (so all supervised methods share a biology-only candidate pool).
    """
    if bundle is None:
        return []
    allowed = bundle.get("allowed")  # cell-type-aware filter (det/spec/V(D)J)
    if bundle.get("tstats") is None:
        if bundle.get("topN") is None:
            return []
        cands = bundle["topN"]["gene"].tolist()
        if BIOLOGY_FILTER:
            cands = [g for g in cands if is_biology_gene(g)]
        if allowed:
            cands = [g for g in cands if g in allowed]
        return cands[:n_top]
    genes = bundle["tstats"]["gene"].tolist()
    if BIOLOGY_FILTER:
        genes = [g for g in genes if is_biology_gene(g)]
    if allowed:
        genes = [g for g in genes if g in allowed]
    if not genes:
        return []
    rel = weighted_relevance(bundle, genes)
    order = np.argsort(-rel)
    return [genes[i] for i in order[:n_top]]


# ============================================================
# Baseline selection methods (per cell type)
# ============================================================
def select_baseline_per_ct(method, bundle, candidates_ct, K):
    """Select K genes from candidates_ct using a baseline method (per cell type)."""
    if not candidates_ct or K <= 0:
        return []
    K = min(K, len(candidates_ct))

    if method == "DE_top":
        # top K by the manuscript relevance |z_i| * C_i
        if bundle.get("tstats") is None:
            return list(candidates_ct[:K])
        rel = weighted_relevance(bundle, candidates_ct)
        order = np.argsort(-rel)
        return [candidates_ct[i] for i in order[:K]]

    if method == "HVG":
        # Label-blind sanity-check baseline: rank by training-pseudobulk variance
        # over the BROAD HVG universe (all ~3000 training HVGs), NOT the
        # disease-informed top-100 candidate pool. This matches the manuscript
        # ("HVG used the same filtered gene universe but did not use diagnosis
        # labels") and keeps HVG a genuine label-blind baseline. Selecting from
        # the disease-informed top-100 would leak label information into HVG and
        # inflate its performance.
        tr = bundle.get("train")
        if tr is None:
            return list(candidates_ct[:K])
        universe = list(tr["genes"])
        if BIOLOGY_FILTER:
            universe = [g for g in universe if is_biology_gene(g)]
        allowed = bundle.get("allowed")  # same cell-type-aware filter as the candidate pool
        if allowed:
            universe = [g for g in universe if g in allowed]
        if len(universe) < 1:
            return list(candidates_ct[:K])
        gene_pos = {g: i for i, g in enumerate(tr["genes"])}
        var_scores = [(g, float(np.var(tr["X"][gene_pos[g], :])))
                      for g in universe if g in gene_pos]
        var_scores.sort(key=lambda x: -x[1])
        return [g for g, _ in var_scores[:K]]

    if method == "mRMR":
        # Greedy mRMR-MID (Peng et al. 2005) using the SAME relevance and
        # redundancy ingredients as the QUBO objective: relevance s_i (cohort-
        # consistency-weighted edgeR test statistic, min-max scaled) and |Pearson rho|
        # between training-donor pseudobulk profiles. The only difference from
        # QUBO is that genes are chosen greedily one at a time rather than by
        # jointly optimising the relevance–redundancy–cardinality objective over
        # the whole panel. Criterion at step t:
        #   s_g - (1/|S|) * sum_{j in S} |rho_{gj}|
        one_bundle = {bundle["cell_type"]: bundle}
        _, R, _ = build_score_and_redundancy(
            one_bundle, candidates_ct,
            score_agg="sum", redundancy_agg="max", score_fn=SCORE_FN)
        # relevance = manuscript |z_i| * C_i (same as QUBO), min-max scaled
        s = np.asarray(weighted_relevance(bundle, candidates_ct), dtype=float)
        s = (s - s.min()) / (s.max() - s.min()) if s.max() > s.min() else np.zeros_like(s)
        Cabs = np.abs(np.asarray(R, dtype=float))
        np.fill_diagonal(Cabs, 0.0)
        N = len(candidates_ct)
        selected = [int(np.argmax(s))]
        while len(selected) < K and len(selected) < N:
            remaining = [i for i in range(N) if i not in selected]
            red = Cabs[np.ix_(remaining, selected)].mean(axis=1)
            score = s[remaining] - red
            selected.append(remaining[int(np.argmax(score))])
        return [candidates_ct[i] for i in selected]

    if method == "LASSO":
        # Pure L1 logistic.  Fit LogRegL1 with various C and take top K by |coef|.
        return _fit_l1_or_en_select(
            bundle, candidates_ct, K, model_kind="L1",
        )

    if method == "ElasticNet":
        # True Elastic Net (L1 + L2 mix, l1_ratio=0.5).
        # NOTE: L2 keeps weights non-zero, so we always select top K by |coef|
        #       (matching the same fixed K as the other methods).
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
    """Inner-CV selection over the configured panel size(s) in K_GRID for a
    baseline method on a single cell type (K_GRID = [10] in the primary
    benchmark; multiple values only in the panel-size sweep)."""
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
        auc = inner_cv_metric_one_ct(bundle, sel, n_inner=INNER_CV_FOLDS,
                                     seed=seed_base + 100, n_repeats=N_INNER_REPEATS,
                                     metric=SELECTION_METRIC)
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


def inner_cv_auc_one_ct(bundle, gene_subset, n_inner=3, seed=0, n_repeats=1):
    """Repeated inner CV AUC for one cell type's classifier on its train pseudobulk.

    Runs `n_repeats` independent donor-stratified `n_inner`-fold partitions (each
    with a different shuffle) and pools all fold AUCs. Repeated CV reduces the
    variance of hyperparameter selection on small cohorts and is applied
    identically to QUBO and every baseline, so the comparison stays fair.
    """
    info = build_X_y_per_ct(bundle, gene_subset, "train")
    if info is None or len(info["X"]) < 2 * n_inner:
        return np.nan
    X, y = info["X"], info["y"]
    n = len(X)
    aucs = []
    for rep in range(max(1, n_repeats)):
        rng = np.random.default_rng(seed + rep * 1000)
        pos = [i for i in range(n) if y[i] == 1]
        neg = [i for i in range(n) if y[i] == 0]
        rng.shuffle(pos); rng.shuffle(neg)
        fold_assign = [[] for _ in range(n_inner)]
        for i, p in enumerate(pos): fold_assign[i % n_inner].append(p)
        for i, p in enumerate(neg): fold_assign[i % n_inner].append(p)
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
            a = roc_auc(yv, clf.predict_proba(Xvz))
            if not np.isnan(a):
                aucs.append(a)
    return float(np.mean(aucs)) if aucs else np.nan


def _cv_score(X, y, clf_factory, n_inner=5, seed=0, n_repeats=5, metric="auc"):
    """Repeated stratified inner-CV score for a given (raw) X, y and classifier.

    Standardises within each inner-train fold. metric: "auc" or "mcc" (at 0.5).
    Used for fair, training-only tuning of the classifier C and for the
    cell-type ensemble weights. No held-out data is touched.
    """
    n = len(X)
    scores = []
    for rep in range(max(1, n_repeats)):
        rng = np.random.default_rng(seed + rep * 1000)
        pos = [i for i in range(n) if y[i] == 1]
        neg = [i for i in range(n) if y[i] == 0]
        if len(pos) < 2 or len(neg) < 2:
            return np.nan
        rng.shuffle(pos); rng.shuffle(neg)
        fa = [[] for _ in range(n_inner)]
        for i, p in enumerate(pos): fa[i % n_inner].append(p)
        for i, p in enumerate(neg): fa[i % n_inner].append(p)
        for k in range(n_inner):
            v = fa[k]
            t = [i for i in range(n) if i not in v]
            if len(set(y[v])) < 2 or len(set(y[t])) < 2:
                continue
            mu = X[t].mean(0); sd = X[t].std(0); sd[sd == 0] = 1.0
            clf = clf_factory().fit((X[t] - mu) / sd, y[t])
            p = np.asarray(clf.predict_proba((X[v] - mu) / sd))
            if metric == "auc":
                s = roc_auc(y[v], p)
            elif metric == "mcc":
                s = mcc_score(y[v], (p >= 0.5).astype(int))
            elif metric == "macro_f1":
                s = _macro_f1(y[v], (p >= 0.5).astype(int))
            else:
                s = roc_auc(y[v], p)
            if not np.isnan(s):
                scores.append(s)
    return float(np.mean(scores)) if scores else np.nan


def _bal_acc(y_true, y_score, threshold=0.5):
    """Balanced accuracy = (sensitivity + specificity) / 2 at the given threshold."""
    y_true = np.asarray(y_true).astype(int)
    pred = (np.asarray(y_score) >= threshold).astype(int)
    tp = int(((pred == 1) & (y_true == 1)).sum()); fn = int(((pred == 0) & (y_true == 1)).sum())
    tn = int(((pred == 0) & (y_true == 0)).sum()); fp = int(((pred == 1) & (y_true == 0)).sum())
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    return float((sens + spec) / 2)


def _macro_f1(y_true, y_pred):
    """Unweighted mean of the two class-wise F1 scores (binary macro-F1)."""
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    f1s = []
    for cls in (0, 1):
        tp = int(((y_pred == cls) & (y_true == cls)).sum())
        fp = int(((y_pred == cls) & (y_true != cls)).sum())
        fn = int(((y_pred != cls) & (y_true == cls)).sum())
        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1s.append(2 * prec * rec / max(prec + rec, 1e-12))
    return float(np.mean(f1s))


def best_threshold_macro_f1(y_true, y_score):
    """Probability threshold maximising Macro-F1 on (y_true, y_score).

    Threshold-sensitivity analysis ONLY. This is computed on the inner-CV
    training-cohort validation predictions and the resulting threshold is then
    applied UNCHANGED to the held-out cohort; held-out labels are never used to
    select it. The same selection rule (maximise Macro-F1) is applied to every
    feature-selection method. Falls back to 0.5 when selection is impossible
    (empty input or a single class / single score value).
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return 0.5
    uniq = np.unique(y_score)
    if len(uniq) < 2:
        return 0.5
    mids = (uniq[:-1] + uniq[1:]) / 2.0
    cands = np.concatenate([[0.5], mids])
    best_t, best_f = 0.5, -1.0
    for t in cands:
        f = _macro_f1(y_true, (y_score >= t).astype(int))
        if f > best_f:
            best_f, best_t = f, float(t)
    return best_t


def fit_classifier_tuned(X_raw, y, seed=0):
    """Fit LogRegL2 with C tuned by inner CV on training data (fair, all methods).

    Returns (fitted_clf_on_standardised_X, mu, sd, best_C). X_raw is unstandardised.
    """
    best_C = 1.0
    if TUNE_CLF_C and len(CLF_C_GRID) > 1:
        best_s = -np.inf
        for C in CLF_C_GRID:
            s = _cv_score(X_raw, y, (lambda C=C: LogRegL2(C=C, max_iter=200)),
                          n_inner=INNER_CV_FOLDS, seed=seed, n_repeats=N_INNER_REPEATS,
                          metric="auc")
            if not np.isnan(s) and s > best_s:
                best_s, best_C = s, C
    mu = X_raw.mean(0); sd = X_raw.std(0); sd[sd == 0] = 1.0
    clf = LogRegL2(C=best_C, max_iter=200).fit((X_raw - mu) / sd, y)
    return clf, mu, sd, best_C


def inner_cv_metric_one_ct(bundle, gene_subset, n_inner=5, seed=0, n_repeats=1, metric=None):
    """Repeated inner-CV selection score for a gene panel, using SELECTION_METRIC
    (auc / mcc / macro_f1). Fixed classifier (BASE_CLF_FACTORY) for speed; the
    same metric is used for QUBO and every baseline so selection stays fair.
    """
    metric = metric or SELECTION_METRIC
    info = build_X_y_per_ct(bundle, gene_subset, "train")
    if info is None or len(info["X"]) < 2 * n_inner:
        return np.nan
    return _cv_score(info["X"], info["y"], BASE_CLF_FACTORY,
                     n_inner=n_inner, seed=seed, n_repeats=n_repeats, metric=metric)


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
                for eta in ETA_GRID:
                    seed = seed_base + K * 10 + int(lam) + int(gamma * 10) + int(eta * 7)
                    sel, E = select_genes_for_cell_type(
                        bundle, candidates_ct, k=K, lam=lam, gamma=gamma,
                        seed=seed, score_fn=SCORE_FN, alpha_batch=ALPHA_BATCH,
                        sa_reads=SA_READS_GRID, sa_sweeps=SA_SWEEPS_GRID, eta=eta)
                    if len(sel) < 3:
                        continue
                    auc = inner_cv_metric_one_ct(bundle, sel, n_inner=INNER_CV_FOLDS,
                                                 seed=seed_base + 100, n_repeats=N_INNER_REPEATS,
                                                 metric=SELECTION_METRIC)
                    log.append(dict(K=K, lam=lam, gamma=gamma, eta=eta, k_actual=len(sel),
                                    energy=E, inner_cv_auc=auc, n_genes=len(sel)))
                    score = -auc if not np.isnan(auc) else float('inf')
                    if best is None or score < best["score"]:
                        best = dict(K=K, lam=lam, gamma=gamma, eta=eta, k_actual=len(sel),
                                    inner_cv_auc=auc, score=score)
    if best is None:
        return None, log
    # final SA at chosen (K, λ, γ, η)
    sel_final, E_final = select_genes_for_cell_type(
        bundle, candidates_ct, k=best["K"], lam=best["lam"], gamma=best["gamma"],
        seed=seed_base, score_fn=SCORE_FN, alpha_batch=ALPHA_BATCH,
        sa_reads=SA_READS, sa_sweeps=SA_SWEEPS, eta=best.get("eta", 0.0))
    best["selected"] = sel_final
    best["energy_final"] = E_final
    return best, log


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
    weight_per_ct = {}      # ct -> ensemble weight (max(0, inner-CV MCC) on train)

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
        # classifier with C tuned by inner CV on train (fair, all methods)
        clf, mu, sd, best_C = fit_classifier_tuned(tr["X"], tr["y"], seed=SEED)
        Xvz = (va["X"] - mu) / sd
        # ensemble weight = max(0, inner-CV MCC) on train (training-only)
        w = _cv_score(tr["X"], tr["y"], (lambda C=best_C: LogRegL2(C=C, max_iter=200)),
                      n_inner=INNER_CV_FOLDS, seed=SEED + 7, n_repeats=N_INNER_REPEATS,
                      metric="mcc")
        weight_per_ct[ct] = max(0.0, w) if not np.isnan(w) else 0.0
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
            ps, ws = [], []
            for ct in probs_per_ct:
                if d in probs_per_ct[ct]:
                    ps.append(probs_per_ct[ct][d])
                    ws.append(weight_per_ct.get(ct, 1.0))
            if not ps:
                continue
            ws = np.asarray(ws, dtype=float)
            if ENSEMBLE_AGG == "weighted_mcc" and ws.sum() > 0:
                p_avg = float(np.average(ps, weights=ws))
            else:  # "mean" fallback (also used if all weights are zero)
                p_avg = float(np.mean(ps))
            rows.append({"donor": d, "diagnosis": diag_map.get(d),
                         "prob_MS": p_avg, "n_cell_types_used": len(ps)})
        return rows

    val_rows = aggregate(val_probs_per_ct, val_diag_map)
    held_rows = aggregate(held_probs_per_ct, held_diag_map)
    return val_rows, held_rows, per_ct_diag


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
    out_dir = PROJECT_ROOT / "qubo_run" / tag / tissue
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n========== TISSUE = {tissue}  TAG = {tag}  DEG = {deg_source} ==========")

    fold_records = []
    val_pred_records = []
    held_pred_records = []
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

        # Cell-type-aware candidate filter (detection / specificity / V(D)J),
        # shared by ALL methods (manuscript §2.4). Stored per bundle as "allowed".
        allowed = compute_allowed_genes(bundles)
        for ct, b in bundles.items():
            if b is not None:
                b["allowed"] = allowed.get(ct, set())

        # Per-cell-type grid search for ALL METHODS
        # selected_per_method[method][ct] = list of genes
        selected_per_method = {m: {} for m in METHODS}
        K_per_method = {m: {} for m in METHODS}

        for ct, b in bundles.items():
            cands = candidates_per_cell_type(b, N_PER_CELL_TYPE)
            if len(cands) < 5:
                continue
            seed_ct = SEED + fold * 100 + _ct_seed(ct)

            # QUBO (uses (K, λ, γ) grid). Screen-then-optimise: pass only the
            # top-N_QUBO_SCREEN relevance-ranked candidates to QUBO (sure
            # independence screening); cands is already ranked by |t|.
            if "QUBO" in METHODS:
                # Screen-then-optimise. When SCREEN_GRID lists several sizes, the
                # screen is selected per (cell type, fold) by inner-CV AUC
                # (training only), alongside (K, λ, γ). cands is |t|-ranked.
                best = None
                rescue = rescue_candidates(b, cands, RESCUE_N) if RESCUE_N > 0 else []
                for nscr in SCREEN_GRID:
                    cands_qubo = list(cands[:min(nscr, len(cands))])
                    if rescue:  # two-tier: union univariate screen with multivariate rescue
                        cands_qubo = list(dict.fromkeys(cands_qubo + rescue))
                    # offset is 0 for the default single screen (N_QUBO_SCREEN),
                    # preserving the canonical seed and bit-for-bit reproducibility.
                    best_s, gl = grid_search_per_ct(b, cands_qubo, fold,
                                                    seed_base=seed_ct + (nscr - N_QUBO_SCREEN))
                    for entry in gl:
                        grid_log.append({"tissue": tissue, "fold": fold,
                                         "cell_type": ct, "method": "QUBO",
                                         "screen": nscr, **entry})
                    if (best_s is not None and not np.isnan(best_s["inner_cv_auc"])
                            and (best is None or best_s["inner_cv_auc"] > best["inner_cv_auc"])):
                        best = dict(best_s); best["screen"] = nscr
                if best is not None:
                    selected_per_method["QUBO"][ct] = best["selected"]
                    K_per_method["QUBO"][ct] = best["K"]
                    for g in best["selected"]:
                        selected_log.append({"tissue": tissue, "fold": fold,
                                             "cell_type": ct, "method": "QUBO",
                                             "K_chosen": best["K"], "gene": g,
                                             "inner_cv_auc": best["inner_cv_auc"]})

            # Baselines (K grid only, no λ,γ)
            for m_name in [m for m in METHODS if m != "QUBO"]:
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
                val_rows, held_rows, per_ct = fit_predict_stacking(bundles, sel)
            else:
                val_rows, held_rows, per_ct = fit_predict_ensemble(bundles, sel)
            if not val_rows:
                continue
            # save predictions only for QUBO (to keep CSV manageable)
            if m_name == "QUBO":
                for r in val_rows:
                    r2 = {**r, "tissue": tissue, "fold": fold, "set": "val"}
                    val_pred_records.append(r2)
                for r in held_rows:
                    r2 = {**r, "tissue": tissue, "fold": fold, "set": "heldout"}
                    held_pred_records.append(r2)

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
            m["val_macro_f1"] = _macro_f1(y_v, (np.asarray(p_v) >= 0.5).astype(int))
            m["val_bal_acc"] = _bal_acc(y_v, p_v)
            if len(y_h):
                a, f = acc_f1(y_h, p_h); m["held_acc"] = a; m["held_f1"] = f
                m["held_mcc"] = mcc_score(y_h, p_h)
                m["held_macro_f1"] = _macro_f1(y_h, (np.asarray(p_h) >= 0.5).astype(int))
                m["held_bal_acc"] = _bal_acc(y_h, p_h)
            else:
                m["held_acc"] = np.nan; m["held_f1"] = np.nan
                m["held_mcc"] = np.nan
                m["held_macro_f1"] = np.nan; m["held_bal_acc"] = np.nan

            # --- Threshold sensitivity analysis (NOT the primary analysis) ----
            # Select the classification threshold by maximising Macro-F1 on the
            # inner-CV training-cohort validation predictions (p_v), then apply it
            # UNCHANGED to the held-out cohort. Done per method and per outer split
            # (fold); held-out labels are never used to choose the threshold, and
            # the same rule is applied to every method. ROC-AUC / PR-AUC are
            # ranking metrics and are unaffected; the fixed-0.5 metrics above
            # remain the primary results.
            thr = best_threshold_macro_f1(y_v, p_v)
            m["tuned_threshold"] = thr
            m["val_mcc_tuned"] = mcc_score(y_v, p_v, threshold=thr)
            m["val_macro_f1_tuned"] = _macro_f1(y_v, (np.asarray(p_v) >= thr).astype(int))
            m["val_bal_acc_tuned"] = _bal_acc(y_v, p_v, threshold=thr)
            if len(y_h):
                m["held_mcc_tuned"] = mcc_score(y_h, p_h, threshold=thr)
                m["held_macro_f1_tuned"] = _macro_f1(y_h, (np.asarray(p_h) >= thr).astype(int))
                m["held_bal_acc_tuned"] = _bal_acc(y_h, p_h, threshold=thr)
            else:
                m["held_mcc_tuned"] = np.nan
                m["held_macro_f1_tuned"] = np.nan
                m["held_bal_acc_tuned"] = np.nan
            fold_records.append(m)

        # per-cell-type diagnostic for QUBO
        diag_method = "QUBO" if "QUBO" in METHODS else None
        if diag_method and selected_per_method.get(diag_method):
            _, _, per_ct = fit_predict_ensemble(bundles, selected_per_method[diag_method])
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
    if selected_log:
        pd.DataFrame(selected_log).to_csv(out_dir / f"selected_genes{suf}.csv", index=False)
    if grid_log:
        pd.DataFrame(grid_log).to_csv(out_dir / f"grid_log{suf}.csv", index=False)
    if per_ct_records:
        pd.DataFrame(per_ct_records).to_csv(out_dir / f"per_ct_diag{suf}.csv", index=False)


def main():
    if len(sys.argv) < 5:
        print("Usage: python3 qubo_pipeline.py <holdout> <deg> <tissue> <fold1> [fold2 ...]")
        sys.exit(1)
    global HOLDOUT_NAME, DEG_SOURCE
    HOLDOUT_NAME = sys.argv[1]
    DEG_SOURCE  = sys.argv[2]
    tissue = sys.argv[3]
    folds = [int(x) for x in sys.argv[4:]]
    run_for_tissue(HOLDOUT_NAME, DEG_SOURCE, tissue, folds)


if __name__ == "__main__":
    main()
