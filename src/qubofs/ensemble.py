"""Cell-type-weighted soft voting (manuscript ensemble).

Each cell type has its own classifier; a donor's probability of MS is the weighted
mean of the cell-type probabilities available for that donor. The weights are
estimated from inner cross-validation only (out-of-fold predictions on training
donors), never from held-out donors:

    a_c = within-cohort stratified ROC-AUC of cell type c (inner out-of-fold predictions)
    r_c = max(a_c - 0.5, 0)
    w_c = (1 - alpha) / n + alpha * r_c / sum(r)      (uniform weights if every r_c = 0)

With ``alpha = 0`` this reduces to the unweighted soft vote.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def within_cohort_auc(prob: pd.Series, diagnosis: pd.Series, cohort: pd.Series,
                      case_label: str = "MS") -> float:
    """ROC-AUC counting only case-control pairs from the same cohort (ties count 0.5).

    Returns NaN if no cohort contains both classes.
    """
    df = pd.DataFrame({"p": np.asarray(prob, float), "y": np.asarray(diagnosis) == case_label, "c": np.asarray(cohort)})
    num = den = 0.0
    for _, g in df.groupby("c"):
        pos, neg = g.p[g.y].values, g.p[~g.y].values
        if len(pos) and len(neg):
            num += (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
            den += len(pos) * len(neg)
    return num / den if den else float("nan")


def celltype_weights(inner_predictions: pd.DataFrame, alpha: float = 0.5, case_label: str = "MS") -> dict[str, float]:
    """Cell-type weights from inner out-of-fold predictions.

    Parameters
    ----------
    inner_predictions : DataFrame
        Columns ``cell_type``, ``prob``, ``diagnosis``, ``cohort`` (one row per
        training donor and cell type, out-of-fold).
    alpha : float in [0, 1], default 0.5
        Share of the weight that follows the inner performance.
    """
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must lie in [0, 1]")
    cts = sorted(inner_predictions["cell_type"].unique())
    n = len(cts)
    score = {}
    for c in cts:
        g = inner_predictions[inner_predictions["cell_type"] == c]
        a = within_cohort_auc(g["prob"], g["diagnosis"], g["cohort"], case_label)
        score[c] = 0.5 if np.isnan(a) else a
    r = {c: max(score[c] - 0.5, 0.0) for c in cts}
    tot = sum(r.values())
    return {c: (1 - alpha) / n + alpha * (r[c] / tot if tot > 0 else 1.0 / n) for c in cts}


def weighted_soft_vote(predictions: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    """Donor-level probabilities: weighted mean over the cell types available per donor.

    ``predictions`` has columns ``donor``, ``cell_type``, ``prob``.
    """
    x = predictions.assign(w=predictions["cell_type"].map(weights))
    if x["w"].isna().any():
        raise ValueError("a cell type in predictions has no weight")
    return x.groupby("donor").apply(lambda g: float(np.average(g["prob"], weights=g["w"])), include_groups=False)
