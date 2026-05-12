# Methods Section

For: *QUBO-Optimized Cell-Type-Specific Gene Panels Enable Robust Cross-Cohort
Classification of Multiple Sclerosis from Single-Cell RNA Sequencing*

Author: Mizuho Asada
Last updated: 2026-05-01

---

## 1. English (primary, for journal Methods section)

### 1.1 Construction of cell-type-specific pseudobulk profiles

To aggregate single-cell RNA-seq data at the donor level, we constructed
pseudobulk expression profiles for each cell type. For each donor in the
training set, gene expression values were averaged across all cells of the
same cell type, yielding a donor-by-gene matrix. In this matrix, rows
represent donors and columns represent genes, with the feature space
restricted to approximately 3,000 highly variable genes (HVGs). This approach
preserves cell-type-specific expression patterns while ensuring that the
donor is the unit of statistical analysis.

Mean aggregation was used in the primary analysis to facilitate comparisons
across donors. As a sensitivity analysis, aggregation by summation was also
performed to assess the robustness of the results.

> **Design notes:**
> ① The unit of statistical inference is the **donor (person)**, not the
>   individual cell, thereby avoiding pseudoreplication.
> ② Cell-type-resolved information is **preserved** by constructing
>   independent profiles per cell type rather than collapsing across types.

### 1.2 Pre-selection gene filtering

Prior to candidate-gene scoring and QUBO panel selection, transcripts known
to confound differential-expression analyses in droplet-based scRNA-seq were
excluded. The exclusion criteria followed current best-practice
recommendations¹·² and comprised:

**(i)** mitochondrial transcripts (`MT-`, `MTRNR`, `MTATP`, `MTND`), which
serve as quality-control markers and reflect cell stress or apoptosis³·⁴;

**(ii)** cytoplasmic and mitochondrial ribosomal protein genes (`RPL*`,
`RPS*`, `MRPL*`, `MRPS*`, including the acidic ribosomal stalk paralogs
`RPLP0/1/2` and the small-subunit gene `RPSA`), whose expression primarily
reflects translational activity and frequently dominates technical variance
without contributing to cell-state identity²;

**(iii)** heat-shock protein genes (`HSPA*`, `HSPB*`, `HSPD*`, `HSPE*`,
`HSP90*`), which are induced by enzymatic dissociation and represent a
documented technical artifact in single-cell preparations⁵·⁶;

**(iv)** the nuclear long non-coding RNAs `MALAT1` and `NEAT1`, which can
persist as nuclear contamination after cytoplasmic RNA depletion and exhibit
ubiquitously high expression²;

**(v)** classical housekeeping genes (`ACTB`, `ACTG1`, `GAPDH`, `B2M`, `FAU`,
`EEF1A1`)⁷;

**(vi)** the X-inactivation transcripts `XIST` and `TSIX`, to mitigate
sex-related confounding given the mixed-sex composition of the cohorts; and

**(vii)** low-confidence transcripts including unannotated loci (`AC*`,
`AL*`, `AP00*`), long intergenic non-coding RNAs (`LINC*`), and small RNAs
(`MIR*`, `RNU*`, `SNORA*`, `SNORD*`), which are unreliably captured by
poly-A-selected 10x Genomics chemistries.

Translation-related genes for which biological evidence supports a regulatory
rather than purely constitutive function were retained — namely `TPT1`
(translationally controlled tumor protein, with documented anti-apoptotic and
proliferative roles), the elongation factor `EEF2`, eIF4F-pathway components
(`EIF4A1`, `EIF4EBP1`, `EIF5A`), the actin-binding thymosins `TMSB4X` and
`TMSB10`, and the splicing factor `HNRNPH1` — to preserve potential
cell-state-specific signals. After filtering, **835 cell-type-specific
candidate genes** remained as input to the per-cell-type QUBO selection.

### References

1. Heumos L, *et al.* Best practices for single-cell analysis across
   modalities. *Nat Rev Genet* **24**, 550–572 (2023).
2. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq
   analysis: a tutorial. *Mol Syst Biol* **15**, e8746 (2019).
3. Ilicic T, *et al.* Classification of low quality cells from single-cell
   RNA-seq data. *Genome Biol* **17**, 29 (2016).
4. Lun ATL, McCarthy DJ, Marioni JC. A step-by-step workflow for low-level
   analysis of single-cell RNA-seq data with Bioconductor.
   *F1000Res* **5**, 2122 (2016).
5. O'Flanagan CH, *et al.* Dissociation of solid tumour tissues with cold
   active protease for single-cell RNA-seq minimizes conserved
   collagenase-associated stress responses. *Genome Biol* **20**, 210 (2019).
6. van den Brink SC, *et al.* Single-cell sequencing reveals
   dissociation-induced gene expression in tissue subpopulations.
   *Nat Methods* **14**, 935–936 (2017).
7. Eisenberg E, Levanon EY. Human housekeeping genes, revisited.
   *Trends Genet* **29**, 569–574 (2013).

---

## 2. 日本語版（学会抄録・国内誌投稿用）

### 2.1 細胞種特異的 pseudobulk プロファイルの構築

scRNA-seq データを donor レベルに集約するため、細胞種ごとに pseudobulk 発現
プロファイルを構築した。学習セットの各 donor について、同一細胞種の全細胞に
わたり遺伝子発現値を平均化し、donor × 遺伝子マトリクスを得た。本マトリクス
では行が donor、列が遺伝子に対応し、特徴空間は約 3,000 個の highly variable
genes (HVG) に絞り込んだ。本アプローチにより、細胞種特異的な発現パターンを
保持しつつ、**donor を統計解析の単位**とすることが保証される。

主解析には mean aggregation を採用し donor 間比較を容易にした。感度解析として
summation aggregation も実施し、結果の頑健性を確認した。

> **設計上の論点：**
> ① 統計推論の単位は **donor（人）** であり、cell ではない (pseudoreplication
>   を回避)。
> ② **cell type の情報は失わず**、cell type ごとに独立な profile を構築する
>   ことで、type-specific な signal を温存する。

### 2.2 遺伝子選択前のフィルタリング

候補遺伝子のスコアリングおよび QUBO によるパネル選択に先立ち、droplet ベース
scRNA-seq の差次的発現解析において交絡因子となることが知られている転写産物を
除外した。除外基準は現行のベストプラクティス推奨¹·² に従い、以下の 7 群とした：

**(i)** ミトコンドリア由来転写産物 (`MT-`, `MTRNR`, `MTATP`, `MTND`) —
quality control マーカーであり、細胞ストレスやアポトーシスを反映する³·⁴；

**(ii)** 細胞質およびミトコンドリアリボソームタンパク質遺伝子 (`RPL*`,
`RPS*`, `MRPL*`, `MRPS*`。酸性リボソームストーク paralog である `RPLP0/1/2`
および small-subunit gene `RPSA` を含む) — 主として翻訳活性を反映し、
cell-state identity への寄与に乏しく technical variance を支配しがちで
ある²；

**(iii)** 熱ショックタンパク質遺伝子 (`HSPA*`, `HSPB*`, `HSPD*`, `HSPE*`,
`HSP90*`) — 酵素消化に伴い誘導される dissociation artifact として確立されて
いる⁵·⁶；

**(iv)** 核内長鎖 non-coding RNA である `MALAT1` および `NEAT1` — 細胞質 RNA
抽出後も核内残存物として持続し、ubiquitous に高発現を示す²；

**(v)** 古典的 housekeeping 遺伝子 (`ACTB`, `ACTG1`, `GAPDH`, `B2M`, `FAU`,
`EEF1A1`)⁷；

**(vi)** X 染色体不活化転写産物 `XIST` および `TSIX` — 本コホートが男女混合
構成であることに鑑み、性別交絡を抑制する目的で除外；

**(vii)** 低信頼度転写産物として unannotated locus (`AC*`, `AL*`, `AP00*`)、
long intergenic non-coding RNA (`LINC*`)、small RNA (`MIR*`, `RNU*`,
`SNORA*`, `SNORD*`) — poly-A 選択型 10x Genomics chemistry では検出が不安定
である。

一方、構成的発現を超えた制御機能を支持する生物学的エビデンスのある翻訳関連
遺伝子は保持した。具体的には `TPT1` (translationally controlled tumor
protein; anti-apoptotic および増殖促進機能の報告あり)、elongation factor
`EEF2`、eIF4F 経路構成因子 (`EIF4A1`, `EIF4EBP1`, `EIF5A`)、actin 結合
thymosin である `TMSB4X` および `TMSB10`、splicing 因子 `HNRNPH1` を含む。
これは cell-state 特異的シグナルを温存する目的である。フィルタリング後、
**835 個の cell type 特異的候補遺伝子**が cell type 別 QUBO 選択の入力と
なった。

---

## 3. 実装上のフィルタ正規表現 (再現性のため)

```python
import re

# Pre-selection biology filter (manuscript / v6entrue version).
# Aligned with current best practice (Heumos et al. 2023; Luecken & Theis 2019).
HK_PATTERN = re.compile(
    r"^(MT-|MTRNR|MTATP|MTND|"             # (i)   mitochondrial
    r"RPL[0-9]|RPS[0-9]|MRPL|MRPS|"        # (ii)  ribosomal (RPLP*, RPSA retained)
    r"HSP[A0-9]|HSPB|HSPA|HSPD|"           # (iii) heat shock (HSPA/B/D + HSP9*)
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|"  # (v)   classical housekeeping
    r"MALAT1$|NEAT1$|XIST$|TSIX$|"         # (iv,vi) nuclear lncRNA + X-inactivation
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|"    # (vii) uncharacterized loci
    r"MIR[0-9]|RNU[0-9]|SNORA|SNORD)"      # (vii) small RNAs
)

def is_biology_gene(g: str) -> bool:
    """Return True if gene name passes the biology filter."""
    return not bool(HK_PATTERN.match(str(g)))
```

**Genes intentionally retained for documented biological function (not pure
housekeeping):**

- **`RPLP0`, `RPLP1`, `RPLP2`** — acidic ribosomal stalk paralogs;
  contribute to specialized-ribosome translation programs (Genuth & Barna,
  *Mol Cell* 2018) with reported immunomodulatory roles (Wang et al.,
  *J Virol* 2020).
- **`RPSA`** — encodes the 67-kDa laminin receptor (LamR/LRP); central to
  leukocyte adhesion and blood–brain-barrier transmigration (Nelson et al.,
  *J Cell Sci* 2008; DiGiacomo & Meruelo, *Biol Rev* 2016) — directly
  relevant to MS pathobiology.
- **Tier-3 translation/RNA-binding genes** (`TPT1`, `EEF2`, `EIF4A1`,
  `EIF4EBP1`, `EIF5A`, `TMSB4X`, `TMSB10`, `HNRNPH1`) — retained for
  documented regulatory or signaling functions beyond constitutive
  expression.
