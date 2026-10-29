#!/usr/bin/env python3
"""Render quboFS supplementary figures S1–S3 from existing evaluation CSVs.

Run from any directory:
  python scripts/03_evaluation/make_supplementary_figures.py \
      --results-dir "$QUBOFS_RESULTS_DIR" \
      --output-dir figures_oup/supplementary

No benchmark models are refitted. Requires pandas, numpy, matplotlib.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METHODS = ['mRMR', 'LASSO', 'ElasticNet', 'DE_top']
DISPLAY = {'QUBO': 'quboFS', 'QUBO_exact': 'Exact', 'mRMR': 'mRMR',
           'LASSO': 'LASSO', 'ElasticNet': 'Elastic net', 'DE_top': 'DE-top'}
HOS = ['Pappalardo', 'Heming', 'Ramesh', 'Mean']


def read(path, columns):
    if not path.is_file():
        raise FileNotFoundError(f'Required existing evaluation file not found: {path}')
    d = pd.read_csv(path)
    missing = set(columns) - set(d.columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}; found {list(d.columns)}')
    return d


def write(fig, out, stem):
    out.mkdir(parents=True, exist_ok=True)
    for ext, opts in [('pdf', {}), ('svg', {}), ('png', {'dpi': 600})]:
        fig.savefig(out / f'{stem}.{ext}', bbox_inches='tight', **opts)
    plt.close(fig)
    print(f'Created {out / stem}.{{pdf,svg,png}}')


def fig_s1(root, out):
    base = root / 'evaluation/csf_redundancy'
    cell = read(base / 'train_rho_per_celltype.csv', ['comparator', 'cell_type', 'mean_d'])
    paired = read(base / 'train_rho_paired.csv', ['comparator', 'diff_ct_mean', 'lo', 'hi', 'p_holm'])
    fig, ax = plt.subplots(figsize=(7.1, 4.1))
    for row, m in enumerate(METHODS):
        p = paired.loc[paired.comparator == m]
        if len(p) != 1:
            raise ValueError(f'Expected exactly one paired row for {m}')
        p = p.iloc[0]
        c = cell.loc[cell.comparator == m].sort_values('cell_type')
        if len(c) != 8:
            raise ValueError(f'Expected eight cell-type differences for {m}, got {len(c)}')
        y = 3 - row
        ax.scatter(c.mean_d, y + np.linspace(-0.18, 0.18, 8), s=20,
                   marker='o', color='#6e98b0', alpha=.65, zorder=2)
        ax.plot([p.lo, p.hi], [y, y], color='#1a1a1a', lw=2.6, zorder=3)
        ax.plot(p.diff_ct_mean, y, 'D', color='#182c42', markersize=6, zorder=4)
    ax.axvline(0, ls='--', lw=1, color='#777777')
    ax.set_yticks(range(4), [DISPLAY[m] for m in METHODS[::-1]])
    ax.set_xlabel('Difference in mean within-panel |r| (quboFS − comparator)')
    ax.set_title('CSF: paired redundancy differences across eight cell types', loc='left')
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.grid(False)
    # Statistical description (dots, diamonds, CI, Holm-adjusted sign-flip p) is given in the figure legend.
    fig.subplots_adjust(left=0.18, right=0.97, top=0.86, bottom=0.18)
    write(fig, out, 'figureS1_redundancy_paired')


def fig_s2(root, out):
    base = root / 'evaluation/exact_solver'
    tissue_names = ['csf', 'pbmc']
    summary = []
    effects = []
    for tissue in tissue_names:
        d = read(base / tissue / 'selection_summary.csv', ['level', 'identical'])
        sel = d.loc[d.level == 'final, all']
        if len(sel) != 1:
            raise ValueError(f'No unique final, all row in {tissue} selection summary')
        frac = float(sel.identical.iloc[0]); summary.append(frac)
        p = read(base / tissue / 'prediction_metrics.csv', ['scheme','cohort','method','metric','value'])
        p = p.loc[(p.scheme == 'weighted') & (p.method.isin(['QUBO', 'QUBO_exact'])) & (p.metric.isin(['roc', 'within']))]
        metric = 'roc' if tissue == 'csf' else 'within'
        cohort = 'Mean' if tissue == 'csf' else 'PBMC'
        vals = p.loc[(p.metric == metric) & (p.cohort == cohort)]
        if set(vals.method) != {'QUBO','QUBO_exact'}:
            raise ValueError(f'Check exact prediction_metrics.csv ({tissue}): expected {metric}, cohort {cohort}; found {vals.to_dict("records")}')
        effects.append((float(vals.loc[vals.method == 'QUBO','value'].iloc[0]),
                        float(vals.loc[vals.method == 'QUBO_exact','value'].iloc[0])))
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.9))
    tissue_labels = ['CSF', 'PBMC']
    x = np.arange(2)
    axs[0].bar(x, np.array(summary)*100, color='#568daf', width=.5)
    axs[0].set_ylim(0, 100)
    axs[0].set_xticks(x, tissue_labels)
    axs[0].set_ylabel('Identical panels (%)')
    axs[0].set_title('A  Identical gene panels', loc='left')
    for i, v in enumerate(summary): axs[0].text(i, v*100 + 2, f'{v*100:.1f}%', ha='center', fontsize=9)
    w = .28
    axs[1].bar(x-w/2, [e[0] for e in effects], width=w, color='#568daf', label='Simulated annealing')
    axs[1].bar(x+w/2, [e[1] for e in effects], width=w, color='#b7b7b7', label='Exact enumeration')
    axs[1].set_xticks(x, tissue_labels)
    axs[1].set_ylim(0, 1.04)
    axs[1].set_ylabel('ROC-AUC')
    axs[1].set_title('B  Classification performance', loc='left')
    for j in range(2):
        for k in range(2):
            v = effects[j][k]
            axs[1].text(j + (-w/2 if k == 0 else w/2), v+.025, f'{v:.3f}', ha='center', fontsize=8)
    axs[1].legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=7.5, frameon=False)
    for ax in axs:
        ax.grid(False)
        ax.spines[['top','right']].set_visible(False)
    # Description of the ROC-AUC definitions and the reference to Table S3 are given in the figure legend.
    fig.subplots_adjust(left=0.10, right=0.80, top=0.86, bottom=0.14, wspace=0.32)
    write(fig, out, 'figureS2_exact_solver')


def fig_s3(root, out):
    """Paired change in CSF ROC-AUC after excluding Ig V(D)J and haemoglobin genes (excluded - primary).
    The 95% CIs come from the same donor bootstrap resamples applied to both analyses (sensitivity_noIgV_Hb.py)."""
    base = root / 'evaluation/sensitivity_noIgV_Hb'
    d = read(base / 'paired_diffs.csv', ['scheme', 'cohort', 'contrast', 'metric', 'diff', 'lo', 'hi'])
    d = d.loc[(d.scheme == 'weighted') & (d.metric == 'roc')]
    auc = read(base / 'auc_ci.csv', ['filter', 'scheme', 'cohort', 'method', 'roc'])
    auc = auc.loc[auc.scheme == 'weighted'].set_index(['filter', 'cohort', 'method']).roc
    methods = ['QUBO'] + METHODS
    cohorts = [('Pappalardo', 'o'), ('Heming', 's'), ('Ramesh', '^')]
    fig, ax = plt.subplots(figsize=(7.1, 3.6))
    for i, m in enumerate(methods):
        y = len(methods) - 1 - i
        rows = {}
        for coh in HOS:
            r = d.loc[(d.cohort == coh) & (d.contrast == f'{m}: noIgV_Hb - primary')]
            if len(r) != 1:
                raise ValueError(f'Expected one paired difference for {m}, {coh}; got {len(r)}')
            r = r.iloc[0]
            expected = auc[('edger_noIgV_Hb', coh, m)] - auc[('edger_only', coh, m)]
            if not np.isclose(r['diff'], expected, atol=1e-12):
                raise ValueError(f'Paired difference for {m}, {coh} does not match auc_ci.csv')
            rows[coh] = r
        for k, (coh, mk) in enumerate(cohorts):
            ax.plot(rows[coh]['diff'], y + 0.18 + 0.09 * k, mk, color='#6e98b0', alpha=.8, markersize=4.5,
                    markerfacecolor='white', markeredgewidth=1.1, zorder=2)
        r = rows['Mean']
        ax.plot([r.lo, r.hi], [y - 0.12, y - 0.12], color='#1a1a1a', lw=2.6, zorder=3)
        ax.plot(r['diff'], y - 0.12, 'D', color='#182c42', markersize=6, zorder=4)
    ax.axvline(0, ls='--', lw=1, color='#777777')
    ax.set_yticks(range(len(methods)), [DISPLAY[m] for m in methods[::-1]])
    ax.set_xlabel('ΔROC-AUC (Ig V(D)J/Hb excluded − primary)')
    ax.set_title('CSF: change in ROC-AUC after excluding Ig V(D)J and haemoglobin genes', loc='left')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.grid(False)
    for coh, mk in cohorts:
        ax.plot([], [], mk, color='#6e98b0', markerfacecolor='white', markeredgewidth=1.1, label=coh)
    ax.plot([], [], 'D', color='#182c42', label='3-cohort mean (95% CI)')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=7.5, frameon=False)
    # Statistical description is given in the figure legend.
    fig.subplots_adjust(left=0.15, right=0.76, top=0.88, bottom=0.16)
    write(fig, out, 'figureS3_IgV_Hb_exclusion')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--results-dir', type=Path, required=True,
                    help='Value of QUBOFS_RESULTS_DIR (existing 03_evaluation results)')
    ap.add_argument('--output-dir', type=Path, default=Path('figures_oup/supplementary'))
    args=ap.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,
                         'axes.linewidth':.6,'savefig.facecolor':'white'})
    fig_s1(args.results_dir, args.output_dir)
    fig_s2(args.results_dir, args.output_dir)
    fig_s3(args.results_dir, args.output_dir)

if __name__ == '__main__': main()
