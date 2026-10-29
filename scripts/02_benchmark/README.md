# 02_benchmark (Python)

Nested cross-validation of quboFS and four comparators (DE-top, mRMR, LASSO, elastic net)
on the outputs of `01_preprocessing`. Set `QUBOFS_DATA_DIR` and `QUBOFS_RESULTS_DIR` before running.
Use `QUBOFS_JOBS` to control parallel execution. From the repository root:

```bash
bash scripts/02_benchmark/run_all.sh
```

| File | Content |
|---|---|
| `config.py` | fixed settings (K = 10, pool 20, lambda/gamma grid, C grid, SA reads/sweeps, LASSO iterations, seeds, alpha) and paths |
| `selection.py` | candidate universe, relevance, pool, quboFS (SA or exact), comparators, classifier |
| `qubo_core.py` | cohort consistency, relevance/redundancy, QUBO matrix, simulated annealing, NumPy logistic regression |
| `data.py` | pseudobulk loading and per-cell-type data bundles |
| `run_csf.py` | CSF: `inner` (every method x setting on the inner folds) and `outer` (setting choice, held-out prediction) |
| `run_pbmc.py` | PBMC: `inner` and `final` for each evaluation unit |

Variants: `--filter edger_only` (primary) or `edger_noIgV_Hb`; `--solver sa` (primary) or
`exact` (quboFS only, at the setting chosen in the primary analysis); `--lasso-max-iter 20000`
(LASSO only, more coordinate-descent iterations; primary 120); `--lasso-pad zc` (LASSO only,
short panels completed by |z|·C instead of |z|).
