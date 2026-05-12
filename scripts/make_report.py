"""
Build performance figures + interactive HTML report from QUBO pipeline outputs.
Reads from qubo_run/ and writes to qubo_run/figures/ and qubo_run/report.html .
"""
import json
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).parent))
from qubo_pipeline import roc_auc, average_precision, roc_curve_np, pr_curve_np, confusion_np

OUT = Path("/sessions/quirky-eloquent-clarke/mnt/outputs/qubo_run")
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)

# ---- Load artifacts ----
metrics = pd.read_csv(OUT / "fold_metrics.csv")
oof = pd.read_csv(OUT / "oof_predictions.csv")
gene_freq = pd.read_csv(OUT / "gene_selection_frequency.csv")
sel_long = pd.read_csv(OUT / "selected_genes_per_fold.csv")
energy = pd.read_csv(OUT / "qubo_energy_per_fold.csv")

models = ["logreg_l2", "logreg_l1", "lda"]
PRETTY = {"logreg_l2": "LogReg (L2)", "logreg_l1": "LogReg (L1)", "lda": "LDA"}
COLORS = {"logreg_l2": "#1f77b4", "logreg_l1": "#ff7f0e", "lda": "#2ca02c"}

# ---- ROC curves: per fold + pooled OOF ----
def fig_roc_pooled():
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="chance")
    for m in models:
        y = oof["y"].values
        p = oof[f"proba_{m}"].values
        fpr, tpr = roc_curve_np(y, p)
        auc = roc_auc(y, p)
        ax.plot(fpr, tpr, color=COLORS[m], lw=2, label=f"{PRETTY[m]}  (pooled AUC={auc:.3f})")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC — pooled out-of-fold predictions (50 donors)")
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "roc_pooled.png", dpi=160)
    plt.close(fig)


def fig_roc_per_fold():
    fig, axes = plt.subplots(1, 5, figsize=(18, 4), sharey=True)
    for i, fold in enumerate(sorted(oof["fold"].unique())):
        sub = oof[oof["fold"] == fold]
        ax = axes[i]
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        for m in models:
            fpr, tpr = roc_curve_np(sub["y"].values, sub[f"proba_{m}"].values)
            auc = roc_auc(sub["y"].values, sub[f"proba_{m}"].values)
            ax.plot(fpr, tpr, color=COLORS[m], lw=1.6,
                    label=f"{PRETTY[m]}  (AUC={auc:.2f})")
        ax.set_title(f"Fold {fold}")
        ax.set_xlabel("FPR")
        if i == 0:
            ax.set_ylabel("TPR")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        ax.grid(alpha=0.3)
    fig.suptitle("Per-fold ROC curves (patient-level CV)")
    fig.tight_layout()
    fig.savefig(FIG / "roc_per_fold.png", dpi=160)
    plt.close(fig)


def fig_pr_pooled():
    fig, ax = plt.subplots(figsize=(6, 5.5))
    base = oof["y"].mean()
    ax.axhline(base, color="k", linestyle="--", alpha=0.5, lw=1, label=f"baseline={base:.2f}")
    for m in models:
        y = oof["y"].values
        p = oof[f"proba_{m}"].values
        rec, prec = pr_curve_np(y, p)
        ap = average_precision(y, p)
        ax.plot(rec, prec, color=COLORS[m], lw=2, label=f"{PRETTY[m]}  (AP={ap:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision–Recall — pooled OOF predictions")
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "pr_pooled.png", dpi=160)
    plt.close(fig)


def fig_metric_bars():
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, metric in zip(axes, ["auc", "ap", "acc", "f1"]):
        means = metrics.groupby("model")[metric].mean().reindex(models)
        stds = metrics.groupby("model")[metric].std().reindex(models)
        bars = ax.bar(range(len(models)), means.values,
                      yerr=stds.values, capsize=5,
                      color=[COLORS[m] for m in models],
                      edgecolor="black", linewidth=0.6, alpha=0.85)
        # also overlay individual fold dots
        for i, m in enumerate(models):
            vals = metrics.loc[metrics.model == m, metric].values
            ax.scatter([i] * len(vals), vals, color="black", s=14, zorder=5, alpha=0.7)
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels([PRETTY[m] for m in models], rotation=20)
        ax.set_ylim(0, 1.05)
        ax.set_title(metric.upper())
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Patient-level 5-fold CV — performance summary")
    fig.tight_layout()
    fig.savefig(FIG / "metric_bars.png", dpi=160)
    plt.close(fig)


def fig_confusion():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, m in zip(axes, models):
        y = oof["y"].values
        pred = (oof[f"proba_{m}"].values >= 0.5).astype(int)
        cm = confusion_np(y, pred)
        im = ax.imshow(cm, cmap="Blues")
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(["HD", "MS"]); ax.set_yticklabels(["HD", "MS"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
        ax.set_title(f"{PRETTY[m]}")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black",
                        fontsize=14, fontweight="bold")
    fig.suptitle("Confusion matrices on pooled OOF predictions (threshold = 0.5)")
    fig.tight_layout()
    fig.savefig(FIG / "confusion_matrices.png", dpi=160)
    plt.close(fig)


def fig_gene_freq():
    top = gene_freq.head(25)
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.barh(range(len(top))[::-1], top["n_folds"].values, color="#4c78a8")
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels(top["gene"].values, fontsize=10)
    ax.set_xlabel("# folds selected (out of 5)")
    ax.set_title("Top 25 most frequently selected genes")
    ax.grid(axis="x", alpha=0.3)
    ax.set_xlim(0, 5.2)
    fig.tight_layout()
    fig.savefig(FIG / "gene_frequency.png", dpi=160)
    plt.close(fig)


def fig_qubo_energy():
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(energy["fold"].astype(str), energy["energy"], color="#9467bd",
           edgecolor="black", alpha=0.85)
    for i, (f, e) in enumerate(zip(energy["fold"], energy["energy"])):
        ax.text(i, e, f"{e:.1f}", ha="center", va="top", color="white", fontsize=9)
    ax.set_xlabel("Fold")
    ax.set_ylabel("QUBO energy (lower = better)")
    ax.set_title("QUBO solution energy per fold (k=20)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "qubo_energy.png", dpi=160)
    plt.close(fig)


fig_roc_pooled()
fig_roc_per_fold()
fig_pr_pooled()
fig_metric_bars()
fig_confusion()
fig_gene_freq()
fig_qubo_energy()
print("Figures written to", FIG)

# ============================================================
# HTML report
# ============================================================
summary = (metrics.groupby("model")[["auc", "ap", "acc", "f1"]]
           .agg(["mean", "std"])
           .round(3))

# pretty summary df
rows_html = []
for m in models:
    row = []
    row.append(f"<td>{PRETTY[m]}</td>")
    for metric in ["auc", "ap", "acc", "f1"]:
        mu = summary.loc[m, (metric, "mean")]
        sd = summary.loc[m, (metric, "std")]
        row.append(f"<td>{mu:.3f} ± {sd:.3f}</td>")
    rows_html.append("<tr>" + "".join(row) + "</tr>")

per_fold_html = []
for _, r in metrics.iterrows():
    per_fold_html.append(
        f"<tr><td>{int(r.fold)}</td><td>{PRETTY[r.model]}</td>"
        f"<td>{r.auc:.3f}</td><td>{r.ap:.3f}</td><td>{r.acc:.3f}</td><td>{r.f1:.3f}</td>"
        f"<td>{int(r.n_train)}</td><td>{int(r.n_test)}</td><td>{int(r.n_features)}</td></tr>"
    )

top_genes_html = []
for _, r in gene_freq.head(25).iterrows():
    top_genes_html.append(f"<tr><td>{r.gene}</td><td>{int(r.n_folds)}</td></tr>")

html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>MS scRNA QUBO – Performance Report</title>
<style>
body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
margin:32px auto;max-width:1100px;color:#1f2937;line-height:1.55}}
h1{{font-size:26px;border-bottom:2px solid #e5e7eb;padding-bottom:8px}}
h2{{margin-top:28px;color:#1f2937}}
table{{border-collapse:collapse;margin:12px 0;font-size:14px}}
th,td{{border:1px solid #e5e7eb;padding:6px 12px;text-align:center}}
th{{background:#f9fafb}}
img{{max-width:100%;border:1px solid #e5e7eb;border-radius:6px;margin:8px 0}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}
.note{{background:#f3f4f6;padding:12px 16px;border-left:3px solid #6366f1;
border-radius:4px;font-size:14px}}
code{{background:#f3f4f6;padding:1px 6px;border-radius:3px}}
</style></head><body>

<h1>MS vs HD 判別モデル — 性能レポート</h1>
<div class="note">
<b>パイプライン</b>: 患者(group)単位の 5-fold 外側 CV →
細胞種別 pseudobulk (B / Mono / CD4_T) → 統合候補遺伝子のQ行列 →
QUBO (Simulated Annealing, k=20) で遺伝子選択 →
細胞種接頭辞付き特徴量 (B__GENE / Mono__GENE / CD4_T__GENE) で分類器学習 →
ホールドアウト患者で評価。
分類器は LogReg(L2) / LogReg(L1) / LDA を比較。
</div>

<h2>1. 性能サマリ (5-fold CV)</h2>
<table>
<tr><th>モデル</th><th>AUC</th><th>AP</th><th>Accuracy</th><th>F1</th></tr>
{''.join(rows_html)}
</table>
<img src="figures/metric_bars.png">

<h2>2. ROC / Precision–Recall</h2>
<div class="grid">
<img src="figures/roc_pooled.png">
<img src="figures/pr_pooled.png">
</div>
<img src="figures/roc_per_fold.png">

<h2>3. Confusion matrices (pooled OOF, threshold = 0.5)</h2>
<img src="figures/confusion_matrices.png">

<h2>4. 選択遺伝子と QUBO エネルギー</h2>
<div class="grid">
<img src="figures/gene_frequency.png">
<img src="figures/qubo_energy.png">
</div>

<h2>5. Fold ごとの結果</h2>
<table>
<tr><th>Fold</th><th>Model</th><th>AUC</th><th>AP</th><th>Acc</th><th>F1</th>
<th>n_train</th><th>n_test</th><th>n_features</th></tr>
{''.join(per_fold_html)}
</table>

<h2>6. 全 fold で多く選ばれた遺伝子 (Top 25)</h2>
<table>
<tr><th>遺伝子</th><th># folds (max=5)</th></tr>
{''.join(top_genes_html)}
</table>

<h2>7. 再現方法</h2>
<pre><code>python3 qubo_pipeline.py        # 5-fold 学習 + 予測 + モデル保存
python3 make_report.py          # 図と HTML レポート生成</code></pre>
<p><i>生成日: 2026-04-28 — pure-numpy / pandas 実装 (sklearn / scipy / neal なしで動作可)</i></p>
</body></html>
"""

(OUT / "report.html").write_text(html, encoding="utf-8")
print("Report:", OUT / "report.html")
