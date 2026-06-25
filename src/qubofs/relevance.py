"""Cohort-consistency-weighted relevance scoring for cross-cohort feature selection.

For cross-cohort biomarker discovery, candidate base relevance scores (such as
:math:`|z_i|` from pooled edgeR test statistics) are multiplied by a
**cohort-consistency score** :math:`C_i` measuring the fraction of training
cohorts in which the MS-versus-control direction of differential expression
agrees with the majority direction. Genes with consistent direction across
training cohorts retain their full weight, whereas genes with inconsistent
direction are down-weighted before candidate-pool ranking (the accompanying
manuscript, §2.4).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


def per_cohort_log_fold_change(
    pseudobulk: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    donor_col: str = "donor_id",
    diagnosis_col: str = "diagnosis",
    cohort_col: str = "cohort",
    case_label: str = "MS",
    control_label: str = "HD",
    min_donors_per_class: int = 1,
) -> pd.DataFrame:
    """Compute per-training-cohort log fold change.

    Parameters
    ----------
    pseudobulk : pandas.DataFrame
        Log-normalised pseudobulk expression matrix, indexed by gene
        (rows), columns are training donor IDs.
    meta : pandas.DataFrame
        Donor metadata containing ``donor_col``, ``diagnosis_col``,
        ``cohort_col``.
    case_label, control_label : str
        Labels for the two diagnostic classes (default "MS" / "HD").
    min_donors_per_class : int, default 1
        Skip cohorts with fewer than this many donors in either class.

    Returns
    -------
    DataFrame : genes × cohort
        Values are ``mean(case) - mean(control)`` per cohort per gene.
        Cohorts with insufficient donors of either class are omitted.
    """
    meta_indexed = meta.set_index(donor_col) if donor_col in meta.columns else meta
    if not meta_indexed.index.is_unique:
        raise ValueError(f"{donor_col} must uniquely identify donors in meta")
    result: dict[str, pd.Series] = {}
    # Group donors by cohort
    grouped: dict[str, dict[str, list[str]]] = {}
    for donor in pseudobulk.columns:
        if donor not in meta_indexed.index:
            continue
        cohort = meta_indexed.loc[donor, cohort_col]
        dx = meta_indexed.loc[donor, diagnosis_col]
        grouped.setdefault(cohort, {case_label: [], control_label: []})
        if dx in grouped[cohort]:
            grouped[cohort][dx].append(donor)
    for cohort, donors in grouped.items():
        if len(donors[case_label]) < min_donors_per_class or len(donors[control_label]) < min_donors_per_class:
            continue
        case_mean = pseudobulk[donors[case_label]].mean(axis=1)
        ctrl_mean = pseudobulk[donors[control_label]].mean(axis=1)
        result[cohort] = case_mean - ctrl_mean
    return pd.DataFrame(result)


def cohort_consistency_score(
    logfc_per_cohort: pd.DataFrame,
    min_informative_cohorts: int = 2,
) -> pd.Series:
    """Compute per-gene cohort-consistency score.

    :math:`C_i = \\max(\\text{n\\_pos}, \\text{n\\_neg}) / \\text{n\\_nonzero}`,
    i.e. the fraction of informative cohorts whose effect direction agrees with
    the majority direction. Cohorts with zero log fold change for a gene are
    treated as non-informative for that gene.

    Parameters
    ----------
    logfc_per_cohort : pandas.DataFrame
        Genes × cohort log-fold-change matrix from :func:`per_cohort_log_fold_change`.
    min_informative_cohorts : int, default 2
        If fewer than this many training cohorts are available overall (i.e. the
        log-fold-change matrix has fewer columns), every gene's score defaults to
        1.0 (no down-weighting; consistency cannot be assessed).

    Returns
    -------
    pandas.Series
        Indexed by gene, values in [0, 1]. Higher = more consistent direction.
    """
    if logfc_per_cohort.shape[1] < min_informative_cohorts:
        return pd.Series(1.0, index=logfc_per_cohort.index)
    signs = np.sign(logfc_per_cohort.values)
    n_pos = (signs > 0).sum(axis=1)
    n_neg = (signs < 0).sum(axis=1)
    n_total = (signs != 0).sum(axis=1)
    majority = np.maximum(n_pos, n_neg)
    C = np.where(n_total > 0, majority / np.maximum(n_total, 1), 0.0)
    return pd.Series(C, index=logfc_per_cohort.index)


@dataclass
class CohortConsistency:
    """Cohort-consistency-weighted relevance scorer.

    Combines a pooled per-gene relevance score (e.g. ``|edgeR statistic|``) with
    a cohort-consistency score derived from per-training-cohort log fold change
    direction:

    .. math::
        s_i = |z_i| \\times C_i

    where :math:`C_i` is the fraction of training cohorts in which the
    direction of differential expression matches the majority direction.

    Parameters
    ----------
    donor_col, diagnosis_col, cohort_col : str
        Column names in the metadata table for donor identifier, diagnostic
        class label and cohort label.
    case_label, control_label : str
        Case and control class labels (default "MS" / "HD").
    min_donors_per_class : int, default 1
        Minimum donors of each class for a cohort to be informative.
    min_informative_cohorts : int, default 2
        If fewer than this many training cohorts are available overall,
        consistency is set to 1 for every gene (no down-weighting).
    """

    donor_col: str = "donor_id"
    diagnosis_col: str = "diagnosis"
    cohort_col: str = "cohort"
    case_label: str = "MS"
    control_label: str = "HD"
    min_donors_per_class: int = 1
    min_informative_cohorts: int = 2

    consistency_: pd.Series | None = field(default=None, init=False, repr=False)
    logfc_: pd.DataFrame | None = field(default=None, init=False, repr=False)

    def fit(self, pseudobulk: pd.DataFrame, meta: pd.DataFrame) -> "CohortConsistency":
        """Compute and store per-gene cohort-consistency scores.

        Parameters
        ----------
        pseudobulk : pandas.DataFrame
            Log-normalised pseudobulk expression (genes × donors) for training
            donors of one cell type.
        meta : pandas.DataFrame
            Donor metadata containing the columns specified by ``donor_col``,
            ``diagnosis_col`` and ``cohort_col``.
        """
        self.logfc_ = per_cohort_log_fold_change(
            pseudobulk, meta,
            donor_col=self.donor_col, diagnosis_col=self.diagnosis_col,
            cohort_col=self.cohort_col, case_label=self.case_label,
            control_label=self.control_label,
            min_donors_per_class=self.min_donors_per_class,
        )
        self.consistency_ = cohort_consistency_score(
            self.logfc_, min_informative_cohorts=self.min_informative_cohorts
        )
        return self

    def weight(self, base_relevance: pd.Series) -> pd.Series:
        """Multiply ``base_relevance`` (e.g. ``|Wald|``) by the cohort-consistency
        score. Genes missing from the consistency map receive weight 0."""
        if self.consistency_ is None:
            raise RuntimeError("CohortConsistency must be fitted before .weight().")
        c = base_relevance.index.map(self.consistency_).fillna(0.0).values
        return pd.Series(np.asarray(base_relevance.values) * c, index=base_relevance.index)
