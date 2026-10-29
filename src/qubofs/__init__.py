"""qubofs: QUBO-based feature selection for compact, non-redundant biomarker panels.

This package implements the QUBO feature-selection framework described in
the accompanying manuscript by Asada et al. The core pipeline consists of four
modular components:

1. **Cohort-consistency-weighted relevance scoring** (|z_i| * C_i) — see :mod:`relevance`.
2. **QUBO matrix construction and selection** by simulated annealing or exact
   enumeration — see :mod:`qubo`.
3. **End-to-end Pipeline** for per-cell-type selection from pseudobulk — see
   :mod:`pipeline` (candidate universe: genes kept by ``edgeR::filterByExpr``; optional
   filters in :mod:`filter`).
4. **Cell-type-weighted soft voting** of per-cell-type classifiers — see :mod:`ensemble`.

Example
-------
>>> from qubofs import QUBOSelector
>>> import numpy as np
>>> rng = np.random.default_rng(0)
>>> N, n_donors = 20, 30
>>> r = rng.random(N)                    # candidate gene relevance scores
>>> X = rng.standard_normal((n_donors, N))
>>> R = np.abs(np.corrcoef(X.T))         # gene-gene |Pearson| redundancy
>>> np.fill_diagonal(R, 0)
>>> sel = QUBOSelector(K=5, gamma=0.5, lambda_=2.0).fit_select(r, R)
>>> sel.selected_indices_.shape
(5,)
"""
from __future__ import annotations

__version__ = "0.2.0"
__author__ = "Mizuho Asada"
__license__ = "MIT"

from qubofs.ensemble import celltype_weights, weighted_soft_vote, within_cohort_auc
from qubofs.filter import CellTypeFilter, is_igv_or_hb, is_technical_gene, is_vdj_segment
from qubofs.metrics import (
    balanced_accuracy,
    highcorr_pairs,
    macro_f1,
    matthews_corr,
    roc_auc,
    within_panel_redundancy,
)
from qubofs.pipeline import Pipeline
from qubofs.qubo import QUBOSelector, exact_minimum, simulated_annealing
from qubofs.relevance import CohortConsistency

__all__ = [
    "__version__",
    "CellTypeFilter",
    "is_igv_or_hb",
    "celltype_weights",
    "weighted_soft_vote",
    "within_cohort_auc",
    "exact_minimum",
    "is_technical_gene",
    "is_vdj_segment",
    "CohortConsistency",
    "QUBOSelector",
    "Pipeline",
    "simulated_annealing",
    "roc_auc",
    "macro_f1",
    "matthews_corr",
    "balanced_accuracy",
    "within_panel_redundancy",
    "highcorr_pairs",
]
