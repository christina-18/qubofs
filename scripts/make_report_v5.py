"""
make_report_v5.py
==================
Generates figures (PNG) and an HTML report for the v5 pipeline.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v5")
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# load
all_metrics = pd.read_csv(ROOT / "all_fold_metrics.csv")
summary = pd.read_csv(ROOT / "all_method_summary.csv")

METHODS = ["QUBO", "DE_top", "HVG", "LASSO", "ElasticNet"]
CLFS = ["LR_L2", "LR_L1", "LDA"]
TISSUES = ["CSF", "PBMC", "ALL"]
COLORS = {"QUBO": "#d62728", "DE_top": "#1f77b4", "HVG": "#2ca02c",
          "LASSO": "#ff7f0e", "ElasticNet": "#9467bd"}


# ============================================================
# Figure 1: AUC (val & held) bar by method × tissue
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
for ti, tissue in enumerate(TISSUES):
    sub = all_metrics[all_metrics.tissue == tissue]
    for ri, metric in enumerate(["val_auc", "held_auc"]):
        ax = axes[ri, ti]
        for mi, m in enumerate(METHODS):
            for ci, c in enumerate(CLFS):
                rows = sub[(sub.method == m) & (sub.classifier == c)]
                if len(rows) == 0:
                    continue
                vals = rows[metric].dropna()
                if len(vals) == 0:
                    continue
                pos = mi * (len(CLFS) + 0.6) + ci
                ax.bar(pos, vals.mean(), yerr=vals.std(),
                       width=0.8, color=COLORS[m],
                       alpha=0.4 + 0.3 * ci, capsize=2,
                       edgecolor="k", linewidth=0.4)
        ax.axhline(0.5, ls="--", c="gray", lw=0.6)
        ax.set_xticks([mi * (len(CLFS) + 0.6) + 1 for mi in range(len(METHODS))])
        ax.set_xticklabels(METHODS, rotation=30, ha="right", fontsize=8)
        ax.set_title(f"{tissue} — {'Validation' if metric=='val_auc' else 'Hold-out (Pappalardo)'}")
        ax.set_ylim(0, 1.05)
        if ti == 0:
            ax.set_ylabel("AUROC")
fig.suptitle("MS vs HD — AUROC by feature-selection method × classifier (3 classifiers per cluster: L2 / L1 / LDA)", fontsize=11)
fig.tight_layout()
fig.savefig(FIG / "auc_by_method_tissue.png", dpi=140, bbox_inches="tight")
plt.close(fig)


# ============================================================
# Figure 2: Selection frequency (QUBO) per tissue
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ti, tissue in enumerate(TISSUES):
    p = ROOT / tissue / "gene_selection_frequency.csv"
    if not p.exists():
        axes[ti].set_visible(False); continue
    df = pd.read_csv(p)
    qubo = df[df.method == "QUBO"].head(20)
    ax = axes[ti]
    ax.barh(range(len(qubo)), qubo["freq"].values[::-1], color=COLORS["QUBO"])
    ax.set_yticks(range(len(qubo)))
    ax.set_yticklabels(qubo["gene"].tolist()[::-1], fontsize=8)
    ax.set_xlabel("Selection freq across 5 folds")
    ax.set_title(f"QUBO top-20 — {tissue}")
    ax.set_xlim(0, 1.05)
fig.tight_layout()
fig.savefig(FIG / "qubo_selection_frequency.png", dpi=140, bbox_inches="tight")
plt.close(fig)


# ============================================================
# Figure 3: Jaccard index of selected gene sets across folds (QUBO)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
for ti, tissue in enumerate(TISSUES):
    p = ROOT / tissue / "selected_genes_per_fold.csv"
    if not p.exists():
        axes[ti].set_visible(False); continue
    df = pd.read_csv(p)
    df = df[df.method == "QUBO"]
    folds = sorted(df.fold.unique().tolist())
    sets = {f: set(df[df.fold == f]["gene"].tolist()) for f in folds}
    n = len(folds)
    J = np.zeros((n, n))
    for i, fi in enumerate(folds):
        for j, fj in enumerate(folds):
            a, b = sets[fi], sets[fj]
            J[i, j] = len(a & b) / max(len(a | b), 1)
    ax = axes[ti]
    im = ax.imshow(J, cmap="viridis", vmin=0, vmax=1, aspect="equal")
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{J[i,j]:.2f}", ha="center", va="center",
                    color="white" if J[i,j] < 0.5 else "black", fontsize=8)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([f"f{f}" for f in folds]); ax.set_yticklabels([f"f{f}" for f in folds])
    ax.set_title(f"Jaccard — {tissue}")
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6)
fig.suptitle("QUBO selected-gene Jaccard across folds", y=1.02, fontsize=11)
fig.savefig(FIG / "qubo_jaccard.png", dpi=140, bbox_inches="tight")
plt.close(fig)


# ============================================================
# Figure 4: ROC curves (val pooled across folds, QUBO+LR_L1 vs DE_top+LR_L2)
# ============================================================
def pooled_roc(oof_df, label_col="diagnosis", score_col="prob_MS"):
    y = (oof_df[label_col] == "MS").astype(int).values
    s = oof_df[score_col].values
    if len(np.unique(y)) < 2:
        return None
    order = np.argsort(-s)
    y = y[order]
    tp = np.cumsum(y) / max(y.sum(), 1)
    fp = np.cumsum(1 - y) / max((1 - y).sum(), 1)
    return np.r_[0, fp, 1], np.r_[0, tp, 1]


fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
for ti, tissue in enumerate(TISSUES):
    ax = axes[ti]
    f = ROOT / tissue / "oof_predictions.csv"
    if f.exists():
        df = pd.read_csv(f)
        roc = pooled_roc(df)
        if roc:
            fpr, tpr = roc
            ax.plot(fpr, tpr, color=COLORS["QUBO"], lw=2, label="val (QUBO+LR_L2)")
    fh = ROOT / tissue / "heldout_predictions.csv"
    if fh.exists():
        df = pd.read_csv(fh)
        roc = pooled_roc(df)
        if roc:
            fpr, tpr = roc
            ax.plot(fpr, tpr, color=COLORS["DE_top"], lw=2, ls="--",
                    label="held-out (Pappalardo)")
    ax.plot([0, 1], [0, 1], "k:", lw=0.6)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
    ax.set_title(f"ROC — {tissue}")
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "roc_pooled.png", dpi=140, bbox_inches="tight")
plt.close(fig)


# ============================================================
# HTML report
# ============================================================
def df_to_html(df, classes="tab"):
    return df.to_html(index=False, classes=classes, float_format=lambda x: f"{x:.3f}",
                      border=0, escape=False)


# best (method, classifier) per tissue by val_auc and held_auc
val_best = (summary.sort_values("val_auc_mean", ascending=False)
            .groupby("tissue").head(1).reset_index(drop=True))
held_best = (summary.sort_values("held_auc_mean", ascending=False)
             .groupby("tissue").head(1).reset_index(drop=True))

# top recurring genes per tissue
gene_freq_html = []
for tissue in TISSUES:
    p = ROOT / tissue / "gene_selection_frequency.csv"
    if not p.exists():
        continue
    df = pd.read_csv(p)
    qubo = df[df.method == "QUBO"].head(15)[["gene", "freq"]]
    qubo.columns = [f"gene ({tissue})", "freq"]
    gene_freq_html.append(qubo.reset_index(drop=True))

html = f"""
<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>MS vs HD — QUBO Gene Selection v5 (compartment-aware)</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, sans-serif; max-width: 1100px;
         margin: 2rem auto; padding: 0 1rem; color: #222; line-height: 1.5; }}
  h1 {{ border-bottom: 2px solid #444; padding-bottom: 0.3rem; }}
  h2 {{ margin-top: 2rem; border-left: 4px solid #d62728; padding-left: 0.5rem; }}
  table.tab {{ border-collapse: collapse; margin: 1rem 0; font-size: 0.85rem; }}
  table.tab th, table.tab td {{ border-bottom: 1px solid #ddd; padding: 4px 8px; }}
  table.tab th {{ background: #f6f6f6; text-align: left; }}
  img {{ max-width: 100%; }}
  .caption {{ font-size: 0.85rem; color: #555; margin-top: -0.5rem; }}
  code {{ background: #f4f4f4; padding: 1px 5px; border-radius: 3px; }}
</style>
</head>
<body>

<h1>MS vs HD 判別モデル — v5 compartment-aware QUBO パイプライン</h1>

<p>
<b>方針 (Step 0):</b> donor 単位 GroupKFold (5分割, cohort-stratified)、
Pappalardo (PRJNA671484_MS_Tcell) は外部 hold-out テストとして取り置き。
全ての DEG / 候補抽出 / Q行列構築は train donor のみで実行。
由来組織 (CSF / PBMC / ALL) と細胞型 (B / Mono / CD4_T / CD8_T) を明示。
</p>

<h2>1. データ要約</h2>
<ul>
  <li>Seurat: 32,170 features × 385,116 cells. Active assay = RNA (counts/data).</li>
  <li>4 cohort, 50 donor (CV 39 + heldout 11 = Pappalardo).</li>
  <li>共変量: age, sex, batch (= prj/cohort) を t-stat に投入。</li>
  <li>細胞型 × tissue × fold = 60 想定 / 実行成功 52 (skip 8 は PBMC fold4/5 で val が MS 単一クラス)。</li>
</ul>

<h2>2. 性能サマリ — best (method, classifier) per tissue</h2>
<h3>2-1. Validation AUC ベスト</h3>
{df_to_html(val_best[['tissue','method','classifier','val_auc_mean','val_auc_std',
                      'held_auc_mean','held_auc_std']])}
<h3>2-2. Hold-out (Pappalardo) AUC ベスト</h3>
{df_to_html(held_best[['tissue','method','classifier','val_auc_mean','val_auc_std',
                       'held_auc_mean','held_auc_std']])}

<h2>3. AUROC by method × tissue × classifier</h2>
<img src="figures/auc_by_method_tissue.png">
<p class="caption">上段: Validation (5-fold, donor-level)。下段: Hold-out (Pappalardo cohort)。
各クラスタは左から L2 / L1 / LDA。</p>

<h2>4. QUBO 選択遺伝子の安定性 (5-fold 間 Jaccard)</h2>
<img src="figures/qubo_jaccard.png">

<h2>5. QUBO 高頻度遺伝子 (5 fold で繰り返し選ばれた gene)</h2>
<img src="figures/qubo_selection_frequency.png">

<h2>6. ROC (pooled across folds, QUBO + LR_L2)</h2>
<img src="figures/roc_pooled.png">

<h2>7. 全結果テーブル (method × classifier × tissue, mean ± std)</h2>
{df_to_html(summary)}

<h2>8. 出力ファイル</h2>
<ul>
  <li><code>qubo_run_v5/all_fold_metrics.csv</code> — 全 fold × method × classifier の指標</li>
  <li><code>qubo_run_v5/all_method_summary.csv</code> — 集約サマリ</li>
  <li><code>qubo_run_v5/&lt;tissue&gt;/oof_predictions.csv</code> — donor 単位 OOF 予測 (val)</li>
  <li><code>qubo_run_v5/&lt;tissue&gt;/heldout_predictions.csv</code> — Pappalardo 予測</li>
  <li><code>qubo_run_v5/&lt;tissue&gt;/selected_genes_per_fold.csv</code> — fold × method × gene</li>
  <li><code>qubo_run_v5/&lt;tissue&gt;/gene_selection_frequency.csv</code> — 選択頻度</li>
  <li><code>qubo_run_v5/&lt;tissue&gt;/Q_fold&lt;k&gt;.npy</code> — Q 行列 (debug)</li>
</ul>

<h2>9. 既知の制約 / Next steps</h2>
<ul>
  <li>QUBO の (λ, γ) はグリッドサーチ未実施。本ランは λ=2.0, γ=1.0 固定。</li>
  <li>PBMC は fold_4/5 が val=MS のみで skip → 3 fold での評価。</li>
  <li>Random Forest ベースラインは未実装 (sklearn 不在のため)。</li>
  <li>感度解析 (sum 集約) は同じスクリプトで <code>aggregator='sum'</code> 指定で再生成可能。</li>
</ul>

</body>
</html>
"""
(ROOT / "report.html").write_text(html, encoding="utf-8")
print(f"Wrote {ROOT/'report.html'}")
print(f"Figures in {FIG}")
