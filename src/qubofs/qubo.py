"""QUBO formulation and Simulated-Annealing solver.

The QUBO objective for cell-type-specific feature selection is

.. math::
    H(\\mathbf{x}) = -\\alpha\\, \\mathbf{s}^\\top \\mathbf{x}
                  + \\gamma\\, \\mathbf{x}^\\top \\mathbf{R}\\, \\mathbf{x}
                  + \\lambda\\, (\\mathbf{1}^\\top \\mathbf{x} - K)^2

with :math:`\\mathbf{x} \\in \\{0,1\\}^N`. The three terms encode relevance,
pairwise non-redundancy and a soft cardinality constraint. The objective is
solved approximately via classical simulated annealing (default 30 reads × 600
sweeps) for reproducibility on standard hardware; no global-optimality guarantee
is claimed. The QUBO matrix is solver-agnostic and is in principle compatible
with iterated Tabu Search, classical or quantum annealing back-ends.

References
----------
Lucas A (2014) *Front. Phys.* — Ising formulations of many NP problems.
Glover F, Kochenberger G, Du Y (2018) *arXiv:1811.11538* — QUBO tutorial.
Mücke S et al. (2023) *Quantum Mach. Intell.* — QUBO feature selection.
Romero S et al. (2025) *Quantum Mach. Intell.* — scRNA-seq QUBO.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def simulated_annealing(
    relevance: np.ndarray,
    redundancy: np.ndarray,
    K: int,
    *,
    alpha: float = 1.0,
    gamma: float = 0.5,
    lambda_: float = 2.0,
    n_reads: int = 30,
    n_sweeps: int = 600,
    t_max: float = 5.0,
    t_min: float = 0.01,
    seed: int | None = None,
) -> tuple[np.ndarray, float]:
    """Solve the relevance/redundancy/cardinality QUBO by Simulated Annealing.

    Parameters
    ----------
    relevance : ndarray, shape (N,)
        Per-gene relevance scores :math:`s_i`, typically rescaled to [0, 1].
    redundancy : ndarray, shape (N, N)
        Symmetric non-negative pairwise redundancy matrix (e.g. absolute
        Pearson correlations) with zeroed diagonal.
    K : int
        Target panel size (soft constraint).
    alpha, gamma, lambda_ : float
        Coefficients of the relevance, redundancy and cardinality terms.
    n_reads : int, default 30
        Number of independent SA reads. The best solution by energy is returned.
    n_sweeps : int, default 600
        Number of sweeps over the N variables per SA read.
    t_max, t_min : float
        Initial and final annealing temperatures (geometric cooling).
    seed : int or None
        RNG seed for reproducibility.

    Returns
    -------
    x_best : ndarray of int8, shape (N,)
        Binary indicator vector for the best (lowest-energy) solution found.
    energy_best : float
        Energy of the best solution.
    """
    relevance = np.asarray(relevance, dtype=np.float64)
    redundancy = np.asarray(redundancy, dtype=np.float64)
    if relevance.ndim != 1:
        raise ValueError("relevance must be 1-D")
    if redundancy.shape != (relevance.size, relevance.size):
        raise ValueError("redundancy shape must match relevance length")
    N = relevance.size
    if K < 1 or K > N:
        raise ValueError("K must be between 1 and the number of candidate features")
    # The incremental energy update (dE_redundancy = 2*gamma*Rx[i]*dx) and the
    # energy accumulator are exact only when the redundancy matrix has a zero
    # diagonal — a gene is not redundant with itself. A natural input such as an
    # absolute Pearson correlation matrix has 1.0 on the diagonal, which would
    # silently corrupt both the reported energy and the selection. Defensively
    # copy and zero the diagonal so the public API always matches the documented
    # precondition (and never mutates the caller's array).
    redundancy = np.array(redundancy, dtype=np.float64, copy=True)
    np.fill_diagonal(redundancy, 0.0)
    rng = np.random.default_rng(seed)
    Ts = t_max * (t_min / t_max) ** (np.arange(n_sweeps) / max(1, n_sweeps - 1))

    best_x: np.ndarray | None = None
    best_E: float = np.inf

    for _ in range(n_reads):
        x = (rng.random(N) < 0.5).astype(np.float64)
        Rx = redundancy @ x
        s = x.sum()
        E = -alpha * relevance.dot(x) + gamma * x.dot(Rx) + lambda_ * (s - K) ** 2
        for sweep in range(n_sweeps):
            T = Ts[sweep]
            order = rng.permutation(N)
            rand_vec = rng.random(N)
            for k, i in enumerate(order):
                dx = 1.0 - 2.0 * x[i]
                # dx is +/-1, so (s + dx - K)^2 - (s - K)^2 = 2*dx*(s - K) + 1.
                dE = (
                    -alpha * relevance[i] * dx
                    + 2 * gamma * Rx[i] * dx
                    + lambda_ * (2 * dx * (s - K) + 1)
                )
                if dE < 0 or rand_vec[k] < np.exp(-dE / max(T, 1e-9)):
                    x[i] += dx
                    Rx += redundancy[:, i] * dx
                    s += dx
                    E += dE
        if E < best_E:
            best_E = float(E)
            best_x = x.copy()
    assert best_x is not None
    return (best_x > 0.5).astype(np.int8), best_E


@dataclass
class QUBOSelector:
    """Sklearn-style QUBO feature selector with simulated annealing.

    Parameters
    ----------
    K : int, default 10
        Target panel size (soft cardinality constraint).
    alpha, gamma, lambda_ : float
        QUBO term coefficients (relevance, redundancy, cardinality).
    sa_reads : int, default 30
        Number of independent SA runs.
    sa_sweeps : int, default 600
        Sweeps per SA run.
    t_max, t_min : float
        Annealing temperature schedule (geometric cooling).
    rescale_relevance : bool, default True
        If True, rescale ``relevance`` to [0, 1] before SA.
    seed : int or None
        RNG seed.

    Attributes
    ----------
    selected_indices_ : numpy.ndarray
        Indices (positions in the candidate set) of selected genes.
    selected_mask_ : numpy.ndarray of int8
        Binary mask of selected genes.
    energy_ : float
        Final QUBO energy.
    """

    K: int = 10
    alpha: float = 1.0
    gamma: float = 0.5
    lambda_: float = 2.0
    sa_reads: int = 30
    sa_sweeps: int = 600
    t_max: float = 5.0
    t_min: float = 0.01
    rescale_relevance: bool = True
    seed: int | None = None

    selected_indices_: np.ndarray | None = field(default=None, init=False, repr=False)
    selected_mask_: np.ndarray | None = field(default=None, init=False, repr=False)
    energy_: float | None = field(default=None, init=False, repr=False)

    def fit_select(self, relevance: np.ndarray, redundancy: np.ndarray) -> "QUBOSelector":
        """Run SA and store selected indices in ``self.selected_indices_``."""
        r = np.asarray(relevance, dtype=np.float64)
        if self.rescale_relevance and r.size > 1:
            r_min, r_max = r.min(), r.max()
            denom = max(r_max - r_min, 1e-9)
            r = (r - r_min) / denom
        mask, energy = simulated_annealing(
            r, redundancy, self.K,
            alpha=self.alpha, gamma=self.gamma, lambda_=self.lambda_,
            n_reads=self.sa_reads, n_sweeps=self.sa_sweeps,
            t_max=self.t_max, t_min=self.t_min, seed=self.seed,
        )
        self.selected_mask_ = mask
        self.selected_indices_ = np.where(mask == 1)[0]
        self.energy_ = energy
        return self
