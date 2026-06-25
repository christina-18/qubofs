"""qubofs: QUBO-based feature selection for compact, non-redundant biomarker panels.

This package implements the QUBO feature-selection framework described in
the accompanying manuscript by Asada et al. The core pipeline consists of four
modular components:

1. **Cell-type-aware filtering** (Stage 1: detection rate; Stage 2: cell-type
   specificity; Stage 3: V(D)J variable-segment exclusion) — see :mod:`filter`.
2. **Cohort-consistency-weighted relevance scoring** — see :mod:`relevance`.
3. **QUBO matrix construction and simulated-annealing-based selection**
   — see :mod:`qubo`.
4. **End-to-end Pipeline** integrating the above components for scRNA-seq
   pseudobulk feature selection — see :mod:`pipeline`.

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

__version__ = "0.1.2"
__author__ = "Mizuho Asada"
__license__ = "MIT"

from qubofs.filter import CellTypeFilter
from qubofs.metrics import (
    balanced_accuracy,
    macro_f1,
    matthews_corr,
    roc_auc,
    within_panel_redundancy,
)
from qubofs.pipeline import Pipeline
from qubofs.qubo import QUBOSelector, simulated_annealing
from qubofs.relevance import CohortConsistency

__all__ = [
    "__version__",
    "CellTypeFilter",
    "CohortConsistency",
    "QUBOSelector",
    "Pipeline",
    "simulated_annealing",
    "roc_auc",
    "macro_f1",
    "matthews_corr",
    "balanced_accuracy",
    "within_panel_redundancy",
]
