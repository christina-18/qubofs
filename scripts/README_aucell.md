# AUCell 解析 — 実行手順

Slide 4 右の "cell-type × gene set" ヒートマップを、本物の AUCell (Aibar et al. 2017, Nat Methods) で置き換えるための手順。R 環境がローカルにあることを前提とします。

**本解析では 2 種類の gene set を 1 回の AUCell で同時にスコアします:**
- (A) **文献キュレーション 7 セット** (MHC II / Iron / Cytotoxic / Type I IFN / MS GWAS / B-cell DMT / MS DMT targets) → 既知 MS biology の独立検証
- (B) **QUBO panel (CSF/PBMC × cell type)** → QUBO 選択遺伝子が cell-level でも MS を判別できるかの per-cell 検証

## ステップ 1: QUBO panel JSON を生成 (Python、私が既に実行済み)

```bash
cd /Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO
python3 scripts/export_qubo_panels_for_aucell.py
```

`qubo_run_v6/aucell_results/qubo_panels.json` が出力されます (既に存在)。CSF QUBO panel は B / Mono / NK / dnT / gdT の 5 cell type、PBMC は Mono / NK / DC / dnT の 4 cell type 分。

## ステップ 2: R AUCell を実行 (← ここをお願いします)

```bash
cd /Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO
Rscript scripts/run_aucell_analysis.R
```

**所要時間**: 12 GB Seurat で 15-40 分 (cell ranking 構築が律速、1 回終わればあとは速い)。

**初回のみ**:
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "GSEABase"))
install.packages("jsonlite")
```

**メタデータ列名のチェック**: スクリプト冒頭の以下を、あなたの Seurat オブジェクトの列名に合わせてください:

```r
CT_COL    <- "celltype_l2_collapsed"
TISSUE_COL <- "compartment"
DX_COL    <- "Dx"
DONOR_COL <- "donor_id"
COHORT_COL <- "cohort"
```

列名が違うと利用可能列一覧を出して停止します (それを見て修正してください)。

**サブサンプリング** (時間短縮): スクリプト内 `MAX_CELLS_PER_GROUP <- 3000` を有効化すると各 (cell_type × tissue × Dx) で max 3,000 細胞にサブサンプル。AUCell は per-cell 計算なので分布の形は変わりません。

## ステップ 3: 出力 CSV (qubo_run_v6/aucell_results/)

| ファイル | 内容 |
|---|---|
| `cell_aucell_scores.csv.gz` | cells × all gene sets (curated + QUBO) |
| `cell_metadata.csv` | cell_id / donor / cell_type / tissue / Dx / cohort |
| `summary_curated.csv` | (ct × tissue × Dx × set) median/mean for (A) |
| `summary_qubo.csv` | 同 for (B) |
| `ms_vs_hd_diff_curated.csv` | (A) MS-HD diff + Wilcoxon p + BH-FDR |
| `ms_vs_hd_diff_qubo.csv` | (B) 同 |

## ステップ 4: 図生成 + Slide 4 差し替え (Claude が引き継ぎ)

R が終わったら **「AUCell 終わったよ」** と一言ください。私が:

1. `python3 scripts/make_aucell_figure.py` で 2 図を生成
   - **fig4a_aucell_qubo_panels_csf.png** ← Slide 4 右 (★ メイン): QUBO panel が MS 細胞で点灯するか
   - **fig4b_aucell_curated_csf.png** ← サプリ S6/S8: 既知 MS axis の cell type 別活性
2. Slide 4 (JP/EN) を fig4a に差し替え + caption を **"AUCell (Aibar et al. 2017, Nat Methods)"** に更新
3. 必要ならサプリスライドにも fig4b を追加
4. ナレーションも該当箇所を AUCell の正しい説明に書き直し

PBMC 版の図も同時に出力されます (`*_pbmc.png`)。
