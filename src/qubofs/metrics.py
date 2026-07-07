"""Classification and panel-quality metrics.

Lightweight implementations of the metrics reported in the accompanying
manuscript: ROC-AUC, MCC, Macro-F1, Balanced Accuracy and within-panel
redundancy (mean absolute pairwise Pearson correlation). The core metric
functions depend only on numpy (the rank-based ROC-AUC is tie-correct, matching
scikit-learn) so that ``qubofs`` can be used without scikit-learn; the
cross-cohort summary helper ``per_cohort_auc_std`` additionally uses pandas.
"""
from __future__ import annotations

import numpy as np


def _average_ranks(a: np.ndarray) -> np.ndarray:
    """1-based ranks with ties resolved by their average (like scipy rankdata).

    Pure-numpy so the metrics stay dependency-light (numpy only).
    """
    a = np.asarray(a, dtype=float)
    sorter = np.argsort(a, kind="mergesort")
    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size)
    a_sorted = a[sorter]
    obs = np.r_[True, a_sorted[1:] != a_sorted[:-1]]
    dense = obs.cumsum()[inv]                      # 1..G dense group ids
    count = np.r_[np.flatnonzero(obs), a.size]     # group start indices + end
    return 0.5 * (count[dense] + count[dense - 1] + 1)


def _binary_confusion(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[int, int, int, int]:
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tp, tn, fp, fn


def _f1(tp: int, fp: int, fn: int) -> float:
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    return (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0


def macro_f1(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Unweighted mean of per-class F1 scores (binary case)."""
    tp, tn, fp, fn = _binary_confusion(y_true, y_pred)
    f1_pos = _f1(tp, fp, fn)
    f1_neg = _f1(tn, fn, fp)
    return (f1_pos + f1_neg) / 2.0


def matthews_corr(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Matthews correlation coefficient. Chicco & Jurman 2020."""
    tp, tn, fp, fn = _binary_confusion(y_true, y_pred)
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Unweighted mean of sensitivity and specificity."""
    tp, tn, fp, fn = _binary_confusion(y_true, y_pred)
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    return (sens + spec) / 2.0


def roc_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Mann-Whitney rank-based ROC-AUC for binary labels and continuous scores.

    Returns NaN if either class is empty.
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if y_true.size == 0 or len(np.unique(y_true)) < 2:
        return float("nan")
    pos = (y_true == 1)
    n_pos = int(pos.sum())
    n_neg = int(y_true.size - n_pos)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Average ranks so that tied scores are handled exactly like scikit-learn.
    ranks = _average_ranks(y_score)
    sum_ranks_pos = float(ranks[pos].sum())
    return float((sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def within_panel_redundancy(expression: np.ndarray) -> float:
    """Mean absolute pairwise Pearson correlation among genes in a panel.

    Parameters
    ----------
    expression : ndarray, shape (n_samples, n_genes)
        Per-donor (or per-sample) expression of the K selected panel genes,
        used to compute pairwise correlations.

    Returns
    -------
    float
        Mean of the absolute upper-triangular entries of the gene-gene Pearson
        correlation matrix. Returns NaN if fewer than 2 genes or 3 samples.
    """
    X = np.asarray(expression, dtype=float)
    if X.ndim != 2 or X.shape[1] < 2 or X.shape[0] < 3:
        return float("nan")
    # A constant (zero-variance) gene makes np.corrcoef divide by zero and emit a
    # RuntimeWarning; the resulting NaN entries are dropped below by nanmean.
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.corrcoef(X.T)
    triu = R[np.triu_indices(R.shape[0], k=1)]
    vals = np.abs(triu)  # constant genes can yield NaN entries; ignore them
    return float(np.nanmean(vals)) if np.isfinite(vals).any() else float("nan")


def per_cohort_auc_std(per_fold_metrics: "list[dict] | dict") -> float:
    """σ_AUC: standard deviation of per-cohort mean ROC-AUC.

    Parameters
    ----------
    per_fold_metrics : list of dicts or dict
        Each entry must contain ``"holdout"`` (cohort label) and ``"auc"``.
        Per-cohort mean is computed first, then the std across cohorts is
        returned. Matches the cross-cohort robustness convention used in the
        manuscript.
    """
    import pandas as pd
    df = pd.DataFrame(per_fold_metrics) if isinstance(per_fold_metrics, list) else per_fold_metrics
    if "holdout" not in df.columns or "auc" not in df.columns:
        raise ValueError("per_fold_metrics must contain 'holdout' and 'auc'")
    return float(df.groupby("holdout").auc.mean().std())


def highcorr_pairs(expression: np.ndarray, threshold: float = 0.70) -> float:
    """Number of highly-correlated gene pairs within a panel.

    Counts the unordered gene pairs whose absolute Pearson correlation exceeds
    ``threshold`` on the per-donor panel expression. This is an assay-design view
    of within-panel redundancy: the count of near-duplicate marker pairs a panel
    carries (Supplementary Table S8), complementary to the mean-|rho| summary
    returned by :func:`within_panel_redundancy`.

    Parameters
    ----------
    expression : ndarray, shape (n_samples, n_genes)
        Per-donor expression of the K selected panel genes.
    threshold : float, default 0.70
        Absolute-correlation cutoff above which a pair is counted as redundant.

    Returns
    -------
    float
        Count of upper-triangular |Pearson r| entries strictly greater than
        ``threshold``. Returns NaN if fewer than 2 genes or 3 samples.
    """
    X = np.asarray(expression, dtype=float)
    if X.ndim != 2 or X.shape[1] < 2 or X.shape[0] < 3:
        return float("nan")
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.corrcoef(X.T)
    triu = np.abs(R[np.triu_indices(R.shape[0], k=1)])
    if not np.isfinite(triu).any():
        return float("nan")
    return float((triu[np.isfinite(triu)] > threshold).sum())
