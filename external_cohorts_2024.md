# MS scRNA-seq 外部コホート候補 (2024 年公開)

本研究 (4 cohorts: Pappalardo / Heming / Ramesh / Touil; 50 patients, 385K cells) を
Phase 2 で外部検証するための、2024 年公開の大規模 MS CSF scRNA-seq コホート 2 件
の詳細まとめ。

---

## 1. Jacobs et al. 2024 (Cell Reports Medicine)

### 書誌情報

| 項目 | 内容 |
|---|---|
| タイトル | Single-cell analysis of cerebrospinal fluid reveals common features of neuroinflammation |
| 著者 | **Jacobs BM, Gasperi C, Kalluri SR, Al-Najjar R, McKeon MO, Else J, Pukaj A, Held F, Sawcer S, Ban M, Hemmer B** |
| 雑誌 | *Cell Reports Medicine*, 6(1), 101733 |
| 発行 | 2024 年 12 月 (online) / 2025 年 1 月 (issue) |
| DOI | 10.1016/j.xcrm.2024.101733 |
| PMID | 39708811 |
| PMC | PMC11866449 |

### コホート構成

| 群 | 患者数 | 細胞数 |
|---|---|---|
| **MS (largely untreated)** | **123** | **203,220** |
| Other inflammatory ND (OIND) | 19 | 30,796 |
| Infectious ND | 23 | 83,339 |
| Non-inflammatory ND | 36 | (~36,700 推定) |
| **合計** | **201** | **354,055 CSF cells** |

### 重要な特徴

- **未治療 MS 中心**: 123 MS 患者の大部分が未治療 → DMT 治療効果の confounder なし
- **大規模 disease control**: OIND / Infectious ND / Non-inflammatory ND 計 78 例
  → 単純な MS vs HD ではなく **MS vs 多様な ND** の比較が可能
- **マルチセンター**: Cambridge + Munich (TUM) の collaborative study
- **CSF only (PBMC なし)**: 本研究の CSF arm 拡張に向く

### 主要知見

> CSF is distinct from peripheral blood in cellular composition and gene expression.
> The cellular and transcriptional landscape of CSF is altered in neuroinflammation
> but **strikingly similar across different neuroinflammatory disorders**.
> Clonal expansion of CSF lymphocytes was found in all disorders but most pronounced
> in inflammatory diseases.

→ **MS 特異性 vs 共通的炎症シグナル**の切り分けが論点。本研究の cross-cohort
classifier が "neuroinflammation 共通シグナル" を学習しているのか、
"MS-specific" なのかを検証する材料になる。

### データアクセス

- **データ accession 番号は論文 supp / data availability セクションを要確認** (web 検索で取得しきれず)
- 多くの場合、こうした大規模 patient cohort は **EGA (European Genome-phenome Archive)** で
  controlled access になっている可能性が高い
- 本論文の preprocessed counts matrix が GEO や Zenodo に公開されているかも要確認

### 本研究との互換性

| 観点 | 評価 |
|---|---|
| Tissue (CSF) | ✅ Match |
| MS vs control | ⚠️ Control が HD ではなく ND mixed (要 stratification) |
| 治療状態 | ✅✅ Untreated 中心、本研究の Heming 等と整合 |
| Cell type annotation | 要再 annotation (Azimuth predicted.celltype.l2) |
| Patient 数 (MS) | ★ 123 例、本研究 4 cohort 合計 (28 MS) の 4.4 倍 |

---

## 2. Ban / Bredikhin / Huang et al. 2024 (Brain)

### 書誌情報

| 項目 | 内容 |
|---|---|
| タイトル | Expression profiling of cerebrospinal fluid identifies dysregulated antiviral mechanisms in multiple sclerosis |
| 著者 | **Ban M, Bredikhin D, Huang Y** ら (詳細著者リストは論文参照) |
| 雑誌 | *Brain*, 147(2), 554–565 |
| 発行 | 2024 年 2 月 |
| PMID | 38038362 |
| 所属 | Cambridge (Sawcer / Ban lab) + 共同研究機関 |

### コホート構成

| 群 | 患者数 | 細胞数 |
|---|---|---|
| MS | **33** | 48,675 |
| OND (other neurological diseases) | **48** | 48,057 |
| **合計** | **81** | **96,732 CSF cells** |

### プラットフォーム

- **10x Genomics Chromium 3' v2**
- CSF only

### 主要知見

> MS susceptibility variants **rs10271373** (chr 7q11.23) and **rs1059091** are
> **eQTLs for antiviral genes in CSF CD8+ T cells**.
> → MS の発症機序に **viral control の dysregulation** が関与している証拠。

→ EBV / 抗ウイルス応答 (Type I IFN, ISG, OAS family) と GWAS hit を直接結びつけた
研究。**本研究の ISG15 / OAS1 系の発見と高い semantic match**。

### データアクセス

- **EGA accession: `EGAS00001007478`** (controlled access、要 Data Access Committee 申請)
- **コード/ノートブック (公開)**: https://github.com/huangyh09/MSclerosisSrc

### 本研究との互換性

| 観点 | 評価 |
|---|---|
| Tissue (CSF) | ✅ Match |
| MS vs control | ⚠️ Control が OND (Other Neurological Diseases)、HD ではない |
| 患者数 (MS) | 33 (本研究 Ramesh = 17 と同等規模) |
| 細胞数 | 96,732 (本研究 Ramesh と同程度) |
| eQTL / GWAS 統合 | ★ 本研究にない genotype 情報を持つ — biology 拡張に有用 |
| Data access | ⚠️ EGA controlled access — 公開準備に時間要 |

---

## 統合計画 — 本研究 Phase 2 への組み込み

### 推奨手順

1. **Public preprocessed counts の入手可否を確認**
   - Cell Rep Med 2024: 論文 data availability + GEO/Zenodo を確認
   - Brain 2024: EGA controlled access のため、Data Access Committee 申請

2. **Cell type annotation の再実行**
   - 全コホートで Azimuth `predicted.celltype.l2` を統一適用
   - 本研究の 8 broad type (B/Mono/CD4_T/CD8_T/NK/DC/dnT/gdT) に collapse

3. **Pseudobulk + edgeR DEG**
   - Donor 単位の pseudobulk → 既存パイプライン (`extract_pseudobulk_v5_compartment.R`) を流用
   - Cohort 列を追加して LOCO 交差検証に組み込む

4. **Cross-cohort QUBO selection**
   - 既存 4 + 新規 2 = 6 cohorts での **6-fold LOCO** 検証
   - Untreated MS (Cell Rep Med) を独立 hold-out にする設計が biomarker 検証として強い

5. **Disease specificity 検証 (Cell Rep Med 2024 特有)**
   - QUBO panel の predict 確率を OIND / Infectious ND / Non-inflammatory に適用
   - MS-specific vs neuroinflammation-common の切り分け

### Slide での言及 (5/8 ラボセミナー Future Work)

> 本研究の 4 cohort (50 patients、385K cells) に加え、**2024 年に公開された
> Jacobs Cell Rep Med 2024 (123 untreated MS + 78 ND, 354K CSF cells)** および
> **Ban Brain 2024 (33 MS + 48 OND, 97K CSF cells, eQTL 統合)** を取り込んで
> cross-cohort 検証を 6 cohort 規模に拡張する予定です。特に Jacobs コホートは
> 多様な neuroinflammatory disease を含むため、QUBO panel が **MS-specific** か
> **neuroinflammation-common** かを直接検証できる枠組みになります。

---

## Caveats

- 上記のサンプル数・accession は web 検索結果に基づく要約であり、**論文本文・supplementary tables での再確認が必須**
- EGA controlled access の Brain 2024 は、Data Access Committee の審査と利用同意が必要
- Cohort 統合時は **Harmony / scVI / scANVI** などのバッチ補正が必須 (本研究の現状 4 cohort も同様)
- 本研究の "HD" 定義 (健常者) と Cell Rep Med 2024 の "Non-inflammatory ND" は厳密には別物
  → 統合時に control の階層を明示する必要あり

---

*Compiled 2026-05-05. References:*
- *Jacobs et al. 2024: https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(24)00463-4*
- *Ban / Bredikhin / Huang et al. 2024: https://academic.oup.com/brain/article/147/2/554/7457295*
- *Ban code repo: https://github.com/huangyh09/MSclerosisSrc*
