# 日本語ナレーション — 全 7 スライド (アカデミック版)

**プレゼンテーション**: QUBO 最適化による細胞種特異的遺伝子パネルを用いた多発性硬化症の cross-cohort 分類
**発表者**: 浅田 瑞穂 (Mizuho Asada, Ph.D)
**スタイル**: アカデミックなラボミーティング向け、です・ます調、短文を重ねる
**総所要時間**: 約 7 分 (1 分 / スライド)
**日付**: 2026 年 5 月 8 日

---

## Slide 1 — Title  *(≈ 25 秒)*

> 本日は、QUBO 最適化を用いた細胞種特異的遺伝子パネルによる多発性硬化症の cross-cohort 分類について発表します。
>
> 単一細胞 RNA シーケンシングデータから、複数施設で再現性の高いバイオマーカーパネルを設計する手法です。本研究の中心的なアイデアは、遺伝子選択を二次無制約二値最適化、すなわち QUBO の問題として定式化することです。これにより、施設間の技術的差異に頑健なパネル設計を目指しました。
>
> まず自己紹介と研究背景から進めます。

---

## Slide 2 — Background & Data  *(≈ 75 秒)*

> 最初に簡単に自己紹介させていただきます。
>
> 浅田瑞穂と申します。明治薬科大学医療分子解析学研究室の助教を、東京科学大学 (Science Tokyo) 大学院心肺統御麻酔学分野の講師を兼任しております。現在、サバティカルとして MGH に滞在しております。
>
> 主な研究領域は、PK/PD モデリング、機械学習による予測モデルおよび画像解析、ケモインフォマティクスです。MGH では、バイオインフォマティクスの中でも特に遺伝子セレクションに取り組みたいと考えております。本日の発表はその方向性の一例です。
>
> 次に、研究の背景です。
>
> 多発性硬化症は世界で約 280 万人が罹患する慢性自己免疫疾患です。Single-cell RNA-seq の普及により、診断バイオマーカー研究は急速に進展しています。一方で、施設間で実験条件が異なる場合に分類精度が大きく低下するという問題が顕在化しています。具体的には、10x chemistry のバージョン、凍結保存の有無、サンプリングプロトコルなどが影響します。
>
> 既存手法には限界があります。DE-top や HVG は一次元のランキングであり、冗長性を制御できません。LASSO は遺伝子数 K の制御が正則化パラメータ λ 経由の間接調整に留まります。
>
> 本研究では、relevance、non-redundancy、cardinality の 3 軸を同時に最適化する QUBO アプローチを提案します。
>
> データセットは 4 つの公開コホートを統合したものです。50 donor、99 sample、385,116 細胞、32,170 遺伝子から構成されます。Heming は 18 donor、Pappalardo は 11 donor、Ramesh は 17 donor、Touil は 4 donor です。

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
> 遺伝子の絞り込みは 3 段階で実施しました。最初に HVG selection を適用し、32,170 遺伝子から各 cell type 3,000 遺伝子に絞り込みます。次に biology filter を適用しました。Heumos et al. 2023 および Luecken & Theis 2019 のベストプラクティス推奨に準拠したものです。ミトコンドリア、リボソーム、熱ショック、核内 lncRNA、housekeeping を除外します。ただし、RPLP0/1/2 と RPSA は意図的に保持しています。RPLP family は specialized ribosome 仮説の対象です。RPSA はリボソームと注釈されますが、実体は 67-kDa laminin receptor をコードします。白血球の BBB 通過に必須の分子であり、MS 病態において保持すべき遺伝子と判断しました。
>
> その後、edgeR で MS vs HD の差次的発現解析を実施します。モデルには年齢、性別、log10 細胞数、batch を共変量として含めています。DESeq2 および limma-voom も並列で実行し、top-100 リストの overlap が 90% 以上であることを sensitivity check として確認しました。
>
> 最後に、|t| 上位 100 遺伝子を候補プールとして QUBO に渡します。各 instance での候補数は 100 個です。全訓練設定での union は CSF で約 1,090 unique 遺伝子となります。

---

## Slide 4 — Method  *(≈ 80 秒)*

> 次に、QUBO の定式化を説明します。
>
> 遺伝子選択を以下の目的関数で定義します。H(x) = -s'x + γ x'Rx + λ (Σx − K)² です。x は各遺伝子の選択を表す二値ベクトルです。
>
> 第一項は relevance の最大化です。各遺伝子のスコア s_i は |t-statistic|² と定義しました。MS と HD で発現が大きく異なる遺伝子に高いスコアを与えます。
>
> 第二項は redundancy の最小化です。選択された遺伝子ペアの相関行列 R を γ で重み付けして罰則化します。互いに類似した発現パターンを持つ遺伝子の重複選択を抑制します。これは既存の univariate 手法では扱えない部分であり、QUBO の特長です。
>
> 第三項は cardinality 制約です。選択数を K に近づける soft penalty です。本研究では K を 10、20、30 から内側 5-fold CV で自動選択しました。
>
> 解法は古典的 Simulated Annealing です。1 instance あたり約 3 秒で求解できます。同じ目的関数は D-Wave の量子アニーラにも実装可能であり、量子ハードウェアの進歩に応じてシームレスに移行できる future-proof な設計です。
>
> 5 つの遺伝子選択手法を同条件で比較しました。DE-top、HVG、LASSO、Elastic Net、QUBO です。すべての手法で同じ候補プール、同じ K グリッド、同じ L2 logistic regression、同じ 8 cell type soft-voting ensemble を用います。唯一の違いは選択ロジックそのものです。
>
> QUBO のみが relevance、non-redundancy、cardinality の 3 軸を同時に最適化します。評価は 3 cohort × leave-one-cohort-out cross-validation で実施します。

---

## Slide 5 — Results  *(≈ 80 秒)*

> 結果に進みます。CSF compartment における 3 cohort 平均の held-out metrics です。
>
> 上の表をご覧ください。QUBO は AUC 0.788、F1 0.635、MCC 0.258 を達成しました。これらは全 5 手法中 1 位です。AP は 0.846 で Elastic Net 0.870 に Δ 0.024 で僅差負け、σ_AUC は 0.044 で Elastic Net 0.041 と Δ 0.003 で実質互角です。LASSO は 0.779、DE_top は 0.742、HVG は 0.712 でした。
>
> 評価指標を簡単に整理します。AUC は ROC 曲線下面積で、MS と HD の順序付け能力を表します。AP は Precision-Recall 曲線下面積で、少数派検出に強い指標です。F1 は Precision と Recall の調和平均です。MCC は Matthews 相関係数で、クラス不均衡に最もロバストな指標として知られています。σ_AUC は 3 cohort 間の AUC 標準偏差で、cross-site 再現性を表します。
>
> 注目していただきたいのは、QUBO が AUC、F1、MCC の主要 3 分類指標で 1 位であることです。Cross-cohort 安定性を表す σ_AUC は Elastic Net とほぼ同等です。
>
> 下の表は cohort 別の held AUC です。QUBO は Pappalardo 0.807、Heming 0.738、Ramesh 0.819 と、3 cohort 全てで 0.74 から 0.82 のタイトな範囲に収束します。LASSO は 0.72 から 0.85 と幅広い分布を示します。
>
> すなわち QUBO は、未知のコホートに対して安定した予測性能を維持できます。これは臨床応用に向けた重要な性質です。

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

## Slide 7 — Conclusion  *(≈ 55 秒)*

> 最後に take-home messages をまとめます。
>
> 第一に性能です。QUBO は CSF において AUC 0.788、F1 0.635、MCC 0.258 を達成し、主要 3 分類指標で 1 位を獲得しました。
>
> 第二に cross-site 再現性です。3 cohort 間の AUC 標準偏差は 0.044 で、Elastic Net 0.041 とほぼ同等です。LASSO 0.068 や DE_top 0.065 を大きく下回ります。
>
> 第三に生物学的妥当性です。QUBO の選択遺伝子は HLA-DRB1、HLA-DPB1、CD74、IFI30 などの MHC class II 抗原提示軸、ISG15 などの Type I IFN signature、GZMA や CCL5 などの細胞傷害性、FTH1 や FTL の鉄代謝という、MS 病態の中核に直結します。GO 解析では MHC class II 抗原提示が p < 10⁻⁸ で濃縮されました。
>
> 第四に量子計算機との互換性です。二値最適化として定式化されているため、D-Wave 量子アニーラにそのまま実装可能です。古典 Simulated Annealing で約 3 秒/instance と高速に動作します。量子ハードウェアの進歩時にシームレスに移行できる future-proof な設計です。
>
> 結論として、QUBO は最高 AUC、最高 F1 および MCC、Elastic Net と互角の cross-cohort 安定性、生物学的妥当性、量子互換性の 5 軸を同時に満たします。Cross-site 再現性のある MS バイオマーカーパネルとして、臨床応用に向けた強い候補と位置付けられます。
>
> 最後に Future Work について一言です。本研究は pseudobulk アプローチですが、Phase 2 として Multi-Instance Learning による cell-level 拡張を計画しています。各 donor を bag、cells を instances として attention 機構で集約することで、どの細胞 subset が病態に寄与したかを解釈可能にします。QUBO は gene selection の役割を維持しつつ、各 donor から informative cells の coreset を選択するという新たな役割も担い得ます。Joint gene-cell QUBO は D-Wave Leap Hybrid Solver での実装を予定しており、量子計算機実証としても位置付けています。
>
> ご清聴ありがとうございました。質疑をお受けいたします。

---

## ペース・話し方のコツ

- **数字 (385,116、17、448、49、0.788、36 倍、p < 10⁻⁸) の前後に短い間** を置く。聞き手に数字が定着する時間を作る。
- **Wet-lab 用語** (oligoclonal band、iron rim、paramagnetic rim lesion、laminin receptor、specialized ribosome) は丁寧に発音する。
- **スライド間の繋ぎは形式的で短く**：
  - 1 → 2: 「まず自己紹介と研究背景から進めます」
  - 2 → 3: 「では、データ整形のパイプラインを説明します」
  - 3 → 4: 「次に QUBO の定式化に進みます」
  - 4 → 5: 「では結果を見ていきます」
  - 5 → 6: 「次に、選択された遺伝子の生物学的解釈に入ります」
  - 6 → 7: 「以上を踏まえ、結論をまとめます」

## 想定される質問への返答

> Q: なぜ 8 集団に集約したのですか？

A: 「より細粒度の subdivision でも解析は可能です。ただし donor 数が 50 と限られている本研究では、CD4_Naive や CD4_TCM などに細かく分けた場合、各 subdivision の細胞数が 20 cells/donor を満たさない donor が出てきます。8 集団は、現在のサンプルサイズで全 cohort で安定した pseudobulk が組める粒度として選定したものです」

> Q: 鉄代謝の濃縮は生物学的にどう解釈しますか？

A: 「Hametner らが Annals of Neurology 2013 で報告したように、慢性活動性 MS lesion の周縁には iron-laden microglia が rim 状に分布します。これは MRI の SWI または QSM で paramagnetic rim lesion として検出される所見です。FTH1 と FTL はこの iron storage の中核分子です。CSF myeloid 系で発現が上がるのは病態と整合的です。本研究で QUBO がこれらを選択したことは、統計的に強い signal を持つだけでなく、実際に MS lesion で起きている現象を反映していると解釈できます」

> Q: stable 49 個と全 448 個では、どちらが最終パネルですか？

A: 「両者は補完的な見方です。Stable 49 は cohort と fold を超えて再現性の高い core panel です。全 448 は QUBO が探索した広い視野を示します。Enrichment はどちらでも有意ですが、強度は curated set で 36 倍と 2 倍程度の差があります。臨床応用を見据えた最終パネルとしては、stable core の 20 から 49 個程度を提案する方向で検討しています」

> Q: D-Wave 量子アニーラでの実証はされていますか？

A: 「現状は古典 Simulated Annealing で完結しています。1 instance あたり約 3 秒という計算速度のため、現時点で量子アニーラを使う必要性は低いです。一方で、データ規模を拡大して全 HVG 約 2,500 遺伝子を一度に最適化するシナリオでは、D-Wave Leap Hybrid Solver の使用が現実的です。今後の課題と位置付けています」
