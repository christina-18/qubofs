# 03_evaluation (Python)

Statistics and figures from the benchmark results. From the repository root:

```bash
bash scripts/03_evaluation/run_all.sh
```

`figures.py` also writes the manuscript Figures 2 and 3 (`figure2_performance_redundancy`,
`figure3_gene_frequency`; PDF, SVG and 600-dpi PNG, plus the plotted values as CSV) to `figures_oup/` at the repository root, or to
`QUBOFS_FIGURE_DIR`.
 The ensemble rule,
bootstrap and paired tests are in `common.py`; see `docs/reproduction.md` for the list of
scripts and `docs/method_details.md` section 7 for the methods.
