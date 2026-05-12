# 日本語ナレーション — スライド 3 & 6 (アカデミック版)

**プレゼンテーション**: QUBO 最適化による細胞種特異的遺伝子パネルを用いた MS の cross-cohort 分類
**発表者**: 浅田 瑞穂 (Mizuho Asada, Ph.D)
**スタイル**: アカデミックなラボミーティング向け、です・ます調、短い文を重ねる
**ペース**: 落ち着いて、約 70-80 秒/スライド
**日付**: 2026 年 5 月 8 日

---

## Slide 3 — Data Preparation Flow  *(≈ 70 秒)*

> 続いて、データ整形のパイプラインについて説明します。
>
> 入力は、385,116 細胞、50 donor、4 cohort から構成される統合 scRNA-seq オブジェクトです。Pappalardo、Heming、Ramesh の 3 コホートを順に hold-out として cross-validation を実施します。Touil コホートは MS donor を含まないため、訓練データに固定しました。
>
> 細胞種アノテーションには Azimuth の predicted.celltype.l2 を使用しています。この annotation には CD4 Naive、CD4 TCM、CD4 TEM、B naive、B memory など 30 種以上の subtype が定義されています。本研究では donor 数が 50 と限られていることを考慮し、3 つの基準で 8 集団に集約しました。第一に、MS 病態に関与する主要 lymphoid または myeloid 群であること。第二に、各 donor 平均 20 細胞以上で安定した pseudobulk が組めること。第三に、全 cohort で再現性高く annotation されること。具体的には B、Mono、CD4_T、CD8_T、NK、DC、dnT、gdT の 8 集団です。
>
> 次に compartment を CSF、PBMC、ALL の 3 つに分割しました。本日は主軸である CSF について報告します。
>
> ここから pseudobulk の構築です。各 (donor × cell type) について、その細胞種の全細胞の発現を平均化します。これによって donor × 遺伝子の行列が cell type ごとに得られます。本設計の重要点は、統計解析の単位が donor、すなわち「人」になることです。Cell-level 解析では pseudoreplication により false positive が増加することが Squair et al. 2021 で報告されています。Pseudobulk 化はこの問題を回避するための標準的なアプローチです。
>
> 遺伝子の絞り込みは 3 段階で実施しました。最初に HVG selection を適用し、32,170 遺伝子から各 cell type 3,000 遺伝子に絞り込みます。次に biology filter を適用しました。これは Heumos et al. 2023 および Luecken & Theis 2019 のベストプラクティス推奨に準拠したものです。ミトコンドリア、リボソーム、熱ショック、核内 lncRNA、housekeeping を除外します。ただし、RPLP0/1/2 と RPSA は意図的に保持しています。RPLP family は specialized ribosome 仮説の対象です。RPSA はリボソームと注釈されますが、実体は 67-kDa laminin receptor をコードします。白血球の BBB 通過に必須の分子であり、MS 病態において保持すべき遺伝子と判断しました。
>
> その後、edgeR で MS vs HD の差次的発現解析を実施します。モデルには年齢、性別、log10 細胞数、batch を共変量として含めています。DESeq2 および limma-voom も並列で実行し、top-100 リストの overlap が 90% 以上であることを sensitivity check として確認しました。
>
> 最後に、|t| 上位 100 遺伝子を候補プールとして QUBO に渡します。各 instance での候補数は 100 個です。全訓練設定での union は CSF で約 1,090 unique 遺伝子となります。

---

## Slide 6 — Selected Genes & Biology  *(≈ 80 秒)*

> スライド 6 では、QUBO が選択した遺伝子の生物学的妥当性について報告します。
>
> 性能の数値だけでは評価の半分にとどまります。選ばれた遺伝子が生物学的にどのような意味を持つかが、本研究のもう一方の重要な評価軸です。
>
> 数値を整理します。3 cohort × 5 fold で計 69 panels が生成され、QUBO は平均 17 遺伝子 / cell type を選択しました。K は 10、20、30 から内側 CV が自動選択する設計です。全 panels の union は 448 unique 遺伝子です。このうち、各 cell type の panel 半数以上で繰り返し選択された 49 遺伝子を「stable core」と定義し、以降の enrichment 解析の foreground としました。
>
> 左の heatmap が選択頻度を示します。行が遺伝子、列が cell type、色の濃さが「その cell type の panel のうち何 % で選ばれたか」を表します。
>
> 特に注目すべきは Mono の列です。CST3 は 100%、SAT1 と FTL は 80% の panel で選ばれています。HLA-DPB1、IFI30、FTH1、CD74 が上位に並びます。これらは MHC class II 抗原提示と鉄代謝の組み合わせであり、MS の monocyte 病態の中核に対応します。dnT 列では GZMA、ISG15、CCL5 が高頻度で選ばれ、cytotoxic と Type I IFN signature を反映します。NK と gdT は KLRB1、KLRC1、GNLY、GZMA などの cytotoxic 系で固まる傾向が見られます。
>
> 右の dot plot は pathway enrichment 解析の結果です。緑が MS-curated set、紺色が GO biological process を示します。
>
> 最も顕著な結果は鉄代謝セットです。Hametner et al. 2013 の Annals of Neurology 論文では、慢性活動性 MS lesion の周縁に iron-laden microglia が rim 状に分布することが報告されています。MRI の SWI および QSM では paramagnetic rim lesion として検出されます。FTH1 と FTL はこの iron rim の中核分子です。本解析では fold enrichment 36 倍、q 値 2×10⁻³ で有意な濃縮を示しました。Cytotoxic effectors は 18 倍、MHC II 経路は 16 倍、いずれも q 値 10⁻³ レベルで有意です。
>
> 全 QUBO selection 448 遺伝子に視野を広げると、結果はさらに明瞭になります。Reactome の MHC II 経路 10 遺伝子全てが選択されています。IMSGC の MS GWAS 上位 11 hit のうち 9 個が再現されました。具体的には BACH2、CXCR4、IL7R、および HLA 一族です。
>
> 結論として、QUBO の選択は統計的最適性のみに基づくものではありません。MS 研究において長年同定されてきた病態の中核軸に収束しています。すなわち、MHC class II 抗原提示、T 細胞活性化、Type I IFN シグナリング、細胞傷害性、鉄代謝の各軸です。これらは、臨床応用を目指すバイオマーカーパネルに求められる性質と整合します。

---

## ペース・話し方のコツ

- **数字 (385,116、17、448、49、36 倍) の前後に短い間** を置く。AUC 0.788、q = 2×10⁻³ なども同様。
- **Wet-lab 用語** (oligoclonal band、iron rim、paramagnetic rim lesion、laminin receptor、specialized ribosome) は丁寧に発音する。
- スライド間の繋ぎは形式的すぎず短く：
  - 3 → 4: 「では QUBO の定式化に進みます。」
  - 5 → 6: 「次に、選択された遺伝子の生物学的解釈に入ります。」
  - 6 → 7: 「以上を踏まえ、結論をまとめます。」

## 想定される質問への返答

> Q: なぜ 8 集団に集約したのですか？

A: 「より細粒度の subdivision でも解析は可能です。ただし donor 数が 50 と限られている本研究では、CD4_Naive や CD4_TCM などに細かく分けた場合、各 subdivision の細胞数が 20 cells/donor を満たさない donor が出てきます。8 集団は、現在のサンプルサイズで全 cohort で安定した pseudobulk が組める粒度として選定したものです。Cohort を拡張すれば、より細粒度での解析も可能になります」

> Q: 鉄代謝の濃縮は生物学的にどう解釈しますか？

A: 「Hametner らが Annals of Neurology 2013 で報告したように、慢性活動性 MS lesion の周縁には iron-laden microglia が rim 状に分布します。これは MRI の SWI または QSM で paramagnetic rim lesion として検出される所見です。FTH1 と FTL はこの iron storage の中核分子です。CSF myeloid 系で発現が上がるのは病態と整合的です。本研究で QUBO がこれらを選択したことは、統計的に強い signal を持つだけでなく、実際に MS lesion で起きている現象を反映していると解釈できます」

> Q: stable 49 個と全 448 個では、どちらが最終パネルですか？

A: 「両者は補完的な見方です。Stable 49 は cohort と fold を超えて再現性の高い core panel です。全 448 は QUBO が探索した広い視野を示します。Enrichment はどちらでも有意ですが、強度は curated set で 36 倍と 2 倍程度の差があります。臨床応用を見据えた最終パネルとしては、stable core の 20 から 49 個程度を提案する方向で検討しています」
