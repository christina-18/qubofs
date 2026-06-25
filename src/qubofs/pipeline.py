"""End-to-end QUBO feature-selection :class:`Pipeline` for scRNA-seq pseudobulk.

Note: this lightweight package pipeline performs **feature selection only**. The
full leave-one-cohort-out classifier benchmark reported in the manuscript
(Table 2, Figures 2-4) is implemented in ``scripts/`` and is not reproduced by
running ``Pipeline`` alone.

Combines :class:`~qubofs.filter.CellTypeFilter`,
:class:`~qubofs.relevance.CohortConsistency` and
:class:`~qubofs.qubo.QUBOSelector` into a single sklearn-style estimator that
selects per-cell-type gene panels from per-cell-type pseudobulk profiles, a
donor metadata table, and per-gene pooled relevance scores.

Example
-------
>>> from qubofs import Pipeline
>>> pipe = Pipeline(K=10, det_thr=0.7, spec_thr=0.7, exclude_vdj=True,
...                 gamma=0.5, lambda_=2.0)
>>> pipe.fit(
...     pseudobulk_per_celltype={"B": pb_b, "Mono": pb_mono, ...},
...     meta=donor_meta,
...     relevance_per_celltype={"B": wald_b, "Mono": wald_mono, ...},
... )  # doctest: +SKIP
>>> pipe.selected_panels_   # doctest: +SKIP
{"B": ["XBP1", "MZB1", ...], "Mono": [...], ...}
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from qubofs.filter import CellTypeFilter
from qubofs.qubo import QUBOSelector
from qubofs.relevance import CohortConsistency


@dataclass
class Pipeline:
    """End-to-end QUBO feature-selection pipeline (sklearn-style).

    Parameters
    ----------
    K : int, default 10
        Target panel size per cell type.
    det_thr : float, default 0.7
        Stage 1 detection-rate threshold.
    spec_thr : float, default 0.7
        Stage 2 cell-type-specificity threshold.
    exclude_vdj : bool, default True
        Stage 3: exclude immunoglobulin / TCR V(D)J variable & joining segments.
    apply_cohort_consistency : bool, default True
        Weight per-gene relevance by cohort-consistency score before SA.
    n_prefilter : int, default 20
        After filtering and weighting, take the top ``n_prefilter`` genes by
        relevance × consistency before constructing the QUBO matrix
        (Sure Independence Screening; Fan & Lv, 2008). Set to None to skip.
    alpha, gamma, lambda_ : float
        QUBO term coefficients.
    sa_reads : int, default 30
    sa_sweeps : int, default 600
    seed : int or None
        RNG seed.
    donor_col, diagnosis_col, cohort_col : str
        Column names in ``meta``.
    case_label, control_label : str

    Attributes
    ----------
    selected_panels_ : dict[str, list[str]]
        Per-cell-type list of selected gene symbols.
    selector_energy_ : dict[str, float]
        Per-cell-type final QUBO energy.
    consistency_ : dict[str, pandas.Series]
        Per-cell-type cohort-consistency scores (if enabled).
    """

    K: int = 10
    det_thr: float = 0.7
    spec_thr: float = 0.7
    exclude_vdj: bool = True
    apply_cohort_consistency: bool = True
    n_prefilter: int | None = 20
    alpha: float = 1.0
    gamma: float = 0.5
    lambda_: float = 2.0
    sa_reads: int = 30
    sa_sweeps: int = 600
    seed: int | None = 42
    donor_col: str = "donor_id"
    diagnosis_col: str = "diagnosis"
    cohort_col: str = "cohort"
    case_label: str = "MS"
    control_label: str = "HD"

    selected_panels_: dict[str, list[str]] = field(default_factory=dict, init=False, repr=False)
    selector_energy_: dict[str, float] = field(default_factory=dict, init=False, repr=False)
    consistency_: dict[str, pd.Series] = field(default_factory=dict, init=False, repr=False)
    _filter: CellTypeFilter | None = field(default=None, init=False, repr=False)

    def fit(
        self,
        pseudobulk_per_celltype: dict[str, pd.DataFrame],
        meta: pd.DataFrame,
        relevance_per_celltype: dict[str, pd.Series],
    ) -> "Pipeline":
        """Run the full selection pipeline.

        Parameters
        ----------
        pseudobulk_per_celltype : dict[str, DataFrame]
            Cell-type -> log-normalised pseudobulk matrix (genes × donors)
            of training-fold donors.
        meta : DataFrame
            Donor metadata (donor_id, diagnosis, cohort).
        relevance_per_celltype : dict[str, Series]
            Cell-type -> per-gene base relevance score, such as an absolute
            Wald statistic (Series indexed by gene).

        Returns
        -------
        self
        """
        self._filter = CellTypeFilter(
            det_thr=self.det_thr, spec_thr=self.spec_thr, exclude_vdj=self.exclude_vdj
        ).fit(pseudobulk_per_celltype)

        self.selected_panels_ = {}
        self.selector_energy_ = {}
        self.consistency_ = {}

        for ct, pb in pseudobulk_per_celltype.items():
            rel = relevance_per_celltype.get(ct)
            if rel is None or rel.empty:
                continue
            # Stage 1+2+3 filter on the gene index
            keep = [g for g in rel.index if self._filter.passes(g, ct)]
            if len(keep) < self.K:
                continue
            rel_f = rel.loc[keep].copy().astype(float)
            # Cohort-consistency weighting
            if self.apply_cohort_consistency:
                cc = CohortConsistency(
                    donor_col=self.donor_col, diagnosis_col=self.diagnosis_col,
                    cohort_col=self.cohort_col,
                    case_label=self.case_label, control_label=self.control_label,
                ).fit(pb, meta)
                self.consistency_[ct] = cc.consistency_
                rel_f = cc.weight(rel_f)
            # Sort and apply SIS pre-filter
            rel_f = rel_f.sort_values(ascending=False)
            n_top = len(rel_f) if self.n_prefilter is None else self.n_prefilter
            top_idx = rel_f.index[:n_top]
            # Pull pseudobulk for top-N candidates
            cands = [g for g in top_idx if g in pb.index]
            if len(cands) < self.K:
                continue
            X = pb.loc[cands].T.values
            X_mean = X.mean(axis=0)
            X_std = X.std(axis=0) + 1e-9
            Xz = (X - X_mean) / X_std
            R = np.abs(np.corrcoef(Xz.T))
            R = np.nan_to_num(R, nan=0.0, posinf=0.0, neginf=0.0)
            np.fill_diagonal(R, 0.0)
            r_vec = rel_f.loc[cands].values.astype(float)
            selector = QUBOSelector(
                K=self.K,
                alpha=self.alpha, gamma=self.gamma, lambda_=self.lambda_,
                sa_reads=self.sa_reads, sa_sweeps=self.sa_sweeps,
                seed=self.seed,
            ).fit_select(r_vec, R)
            self.selected_panels_[ct] = [cands[i] for i in selector.selected_indices_]
            self.selector_energy_[ct] = float(selector.energy_)
        return self
