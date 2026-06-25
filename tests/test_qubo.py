"""Tests for the QUBO selector and Simulated Annealing solver."""
import itertools

import numpy as np
import pytest

from qubofs.qubo import QUBOSelector, simulated_annealing


def _qubo_energy(x, r, R, K, alpha=1.0, gamma=0.5, lambda_=2.0):
    """Reference QUBO energy H(x) = -alpha s.x + gamma x'Rx + lambda(sum(x)-K)^2."""
    x = np.asarray(x, dtype=float)
    return float(
        -alpha * r.dot(x) + gamma * x.dot(R @ x) + lambda_ * (x.sum() - K) ** 2
    )


def _brute_force_min(r, R, K, **kw):
    N = r.size
    best_x, best_E = None, np.inf
    for bits in itertools.product([0, 1], repeat=N):
        x = np.array(bits, dtype=float)
        E = _qubo_energy(x, r, R, K, **kw)
        if E < best_E:
            best_E, best_x = E, x
    return best_x, best_E


def test_sa_energy_matches_formula_and_global_optimum():
    """At small N with generous SA settings, the solver should (a) report an
    energy consistent with the QUBO formula and (b) recover the exact optimum
    found by brute-force enumeration."""
    rng = np.random.default_rng(3)
    N, K = 8, 3
    r = rng.random(N)
    X = rng.standard_normal((40, N))
    R = np.abs(np.corrcoef(X.T))
    np.fill_diagonal(R, 0.0)

    x_sa, E_sa = simulated_annealing(
        r, R, K=K, alpha=1.0, gamma=0.5, lambda_=2.0,
        n_reads=40, n_sweeps=400, seed=0,
    )
    # (a) returned energy equals the energy of the returned vector
    assert abs(E_sa - _qubo_energy(x_sa, r, R, K)) < 1e-9
    # (b) SA reaches the exact global minimum energy
    _, E_exact = _brute_force_min(r, R, K)
    assert abs(E_sa - E_exact) < 1e-9


def test_simulated_annealing_returns_correct_shape():
    rng = np.random.default_rng(0)
    N = 15
    r = rng.random(N)
    X = rng.standard_normal((30, N))
    R = np.abs(np.corrcoef(X.T))
    np.fill_diagonal(R, 0)
    x, E = simulated_annealing(r, R, K=5, n_reads=4, n_sweeps=100, seed=0)
    assert x.shape == (N,)
    assert x.dtype == np.int8
    assert isinstance(E, float)
    assert x.sum() >= 1  # not trivially empty


def test_simulated_annealing_respects_cardinality_constraint_softly():
    """At reasonable lambda, the SA selection should be close to K."""
    rng = np.random.default_rng(42)
    N = 20
    r = rng.random(N)
    X = rng.standard_normal((30, N))
    R = np.abs(np.corrcoef(X.T))
    np.fill_diagonal(R, 0)
    K = 6
    x, _ = simulated_annealing(
        r, R, K=K, lambda_=5.0, n_reads=8, n_sweeps=300, seed=1,
    )
    # Soft constraint: should be within ±2 of K
    assert abs(int(x.sum()) - K) <= 2


def test_qubo_selector_api():
    rng = np.random.default_rng(0)
    N = 20
    r = rng.random(N)
    X = rng.standard_normal((30, N))
    R = np.abs(np.corrcoef(X.T))
    np.fill_diagonal(R, 0)
    sel = QUBOSelector(K=5, sa_reads=4, sa_sweeps=200, seed=0).fit_select(r, R)
    assert sel.selected_indices_ is not None
    assert sel.selected_mask_ is not None
    assert sel.energy_ is not None
    assert sel.selected_indices_.shape[0] == int(sel.selected_mask_.sum())


def test_simulated_annealing_rejects_invalid_K():
    """K must lie between 1 and the number of candidate features."""
    r = np.ones(5)
    R = np.zeros((5, 5))
    with pytest.raises(ValueError):
        simulated_annealing(r, R, K=0)
    with pytest.raises(ValueError):
        simulated_annealing(r, R, K=6)


def test_simulated_annealing_rejects_bad_redundancy_shape():
    """redundancy must be square and match the relevance length."""
    r = np.ones(5)
    R = np.zeros((4, 4))
    with pytest.raises(ValueError):
        simulated_annealing(r, R, K=2)


def test_simulated_annealing_seed_reproducibility():
    rng = np.random.default_rng(0)
    N = 10
    r = rng.random(N)
    R = np.abs(rng.standard_normal((N, N)))
    R = (R + R.T) / 2
    np.fill_diagonal(R, 0)
    x1, e1 = simulated_annealing(r, R, K=3, n_reads=4, n_sweeps=200, seed=7)
    x2, e2 = simulated_annealing(r, R, K=3, n_reads=4, n_sweeps=200, seed=7)
    assert np.array_equal(x1, x2)
    assert e1 == e2
