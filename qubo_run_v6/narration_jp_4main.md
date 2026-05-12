# 日本語ナレーション — 4 main + 4 supplementary 構成

**プレゼンテーション**: QUBO 最適化による細胞種特異的遺伝子パネルを用いた多発性硬化症の cross-cohort 分類
**発表者**: 浅田 瑞穂 (Mizuho Asada, Ph.D)
**スタイル**: アカデミックなラボミーティング向け、です・ます調、短文を重ねる
**対象**: MS 専門家集団 (基礎背景は省略)
**総所要時間**: 約 5 分 (main 4 枚) + 質疑時に補助スライド

---

## Main Slide 1 — Headline + Background + Dataset  *(≈ 80 秒)*

> 本日は、QUBO 最適化を用いた細胞種特異的遺伝子パネルによる多発性硬化症の cross-cohort 分類について発表します。
>
> 結論を先にお示しします。49 個の stable QUBO 遺伝子により、3 つの hold-out コホートで AUC 0.788 ± 0.044 を達成しました。さらに、選ばれた遺伝子は MS 病態の中核軸である MHC class II 抗原提示、鉄代謝、Type I IFN、細胞傷害性を再現しています。
>
> 3 つの take-home messages があります。第一に性能。AUC、F1、MCC の主要 3 指標で 5 手法中 1 位です。第二に再現性。cross-cohort σ_AUC は 0.044、Elastic Net 0.041 とほぼ同等で、LASSO 0.068 の 6 割です。第三に生物学的妥当性。MHC II は 8.5 倍、鉄代謝は 36 倍に濃縮されました。これは Hametner 2013 の iron rim と整合的です。
>
> 背景は技術的課題に絞ります。既存の MS scRNA-seq 解析では、施設間の技術差により分類精度が大きく低下します。具体的には 10x chemistry のバージョン、凍結保存の有無、サンプリングプロトコルなどです。既存手法 DE-top、LASSO、HVG は relevance、redundancy、cardinality を同時に最適化できません。本研究は QUBO によるこの 3 軸同時最適化と cross-cohort validation でこの課題を解決します。古典的な Simulated Annealing で求解し、同じ定式化が D-Wave 量子アニーラにも実装可能です。
>
> データセットは 4 cohort 統合で 50 donor、385,000 細胞、32,000 遺伝子です。Heming 18 donor、Pappalardo 11 donor、Ramesh 17 donor、Touil 4 donor。Touil は MS donor を含まないため訓練に固定し、3 cohort × LOOCV で外部評価します。

---

## Main Slide 2 — Method  *(≈ 75 秒)*

> 次に QUBO の定式化を説明します。
>
> 遺伝子選択を 3 項からなる目的関数で定義します。
>
> 第一項は relevance の最大化です。各遺伝子のスコア s_i は edgeR の |t-statistic| の二乗です。MS と HD で発現差が大きい遺伝子に高いスコアを与えます。
>
> 第二項は redundancy の最小化です。選択された遺伝子ペアの相関行列 R をパラメータ γ で重み付けして罰則化します。互いに類似した発現パターンの遺伝子の重複選択を抑制します。これは既存の univariate 手法では扱えない部分であり、QUBO の核心的な貢献です。
>
> 第三項は cardinality 制約です。選択数を K に近づける soft penalty です。K は 10、20、30 から内側 5-fold CV で自動選択されます。
>
> 解法は古典的 Simulated Annealing で、1 instance あたり約 3 秒です。同じ定式化が D-Wave 量子アニーラにも修正なしで実装可能です。
>
> 5 つの遺伝子選択手法を同条件で比較しました。DE-top、HVG、LASSO、Elastic Net、QUBO です。重要な点として、5 手法すべてで同じ候補プール、同じ K グリッド、同じ L2 logistic 分類器、同じ 8 cell-type soft-voting ensemble を使用しています。違いは選択ロジックそのものだけです。これにより、選択方法そのものの差を切り分けて評価できます。

---

## Main Slide 3 — Results — Performance & Biological Validity  *(≈ 90 秒)*

> 結果に進みます。CSF compartment における 3 cohort 平均の held-out metrics と、選択された遺伝子の生物学的解釈を 1 枚にまとめました。
>
> 上の表をご覧ください。QUBO は AUC 0.788、F1 0.635、MCC 0.258 を達成し、これら主要 3 分類指標で 5 手法中 1 位です。AP は 0.846 で Elastic Net 0.870 に Δ 0.024 で僅差負け、σ_AUC は 0.044 で Elastic Net 0.041 と Δ 0.003 で実質互角です。Cohort 別では Pappalardo 0.807、Heming 0.738、Ramesh 0.819 と、3 cohort すべてで 0.74 から 0.82 のタイトな範囲に収束します。詳細な per-cohort 表は補助スライド S3 に用意しています。
>
> 下半分は選択遺伝子の生物学的解釈です。
>
> 左の heatmap は選択頻度を示します。行が遺伝子、列が cell type、色の濃さがその cell type の panel のうち何 % で選ばれたかを表します。特に注目すべきは Mono の列です。CST3 が 100%、SAT1 と FTL が 80%、HLA-DPB1、IFI30、FTH1、CD74 が上位に並びます。MHC class II 抗原提示と鉄代謝の組み合わせ、まさに MS の monocyte 病態の中核です。dnT 列では GZMA、ISG15、CCL5 が高頻度で、cytotoxic と Type I IFN signature を反映します。
>
> 右の dot plot は pathway enrichment です。緑が MS-curated set、紺色が GO biological process です。最も顕著なのは鉄代謝で、fold enrichment 36 倍。Hametner 2013 の chronic active lesion における iron rim と直接対応します。Cytotoxic effectors は 18 倍、MHC II 経路は 16 倍と、いずれも有意な濃縮を示します。
>
> 全 QUBO selection 448 遺伝子に視野を広げると、結果はさらに明瞭になります。Reactome MHC II 経路の 10 遺伝子全てが選ばれ、IMSGC GWAS 上位 11 hit のうち 9 個が再現されます。BACH2、CXCR4、IL7R、HLA 一族です。

---

## Main Slide 4 — Conclusion + Future Work  *(≈ 65 秒)*

> 最後に take-home messages をまとめます。
>
> 第一に性能。QUBO は CSF において主要 3 分類指標で 1 位を獲得しました。
>
> 第二に cross-site 再現性。Elastic Net と互角の 0.044 で、LASSO や DE_top よりも実質的に優れます。
>
> 第三に生物学的妥当性。MHC class II 抗原提示、Type I IFN、cytotoxicity、鉄代謝という MS 病態の中核に直結します。GO 解析では MHC class II が p < 10⁻⁸ で濃縮されました。
>
> 第四に量子計算機との互換性。二値最適化として定式化されているため、D-Wave 量子アニーラに直接実装可能です。古典 SA で約 3 秒/instance と高速で、量子ハードウェア進歩時にシームレスに移行できます。
>
> 結論として、QUBO は最高 AUC、最高 F1 および MCC、Elastic Net と互角の cross-cohort 安定性、生物学的妥当性、量子互換性の 5 軸を同時に満たします。Cross-site 再現性のある MS バイオマーカーパネルとして、臨床応用に向けた強い候補と位置付けられます。
>
> Future Work として、Multi-Instance Learning による cell-level 拡張を計画しています。各 donor を bag、cells を instances として attention 機構で集約することで、どの細胞 subset が病態に寄与したかを解釈可能にします。QUBO は gene selection の役割を維持しつつ、各 donor から informative cells の coreset を選択するという新たな役割も担い得ます。Joint gene-cell QUBO は D-Wave Leap Hybrid Solver での実装を予定しています。
>
> ご清聴ありがとうございました。質疑をお受けいたします。

---

## 補助スライドの使い方 (質疑時)

質問に応じて該当する補助スライドに切り替えて説明します。

### S1 — Presenter Background (発表者背景)

**質問例**: 「ご経歴を教えてください」「どのような研究分野ですか」「ケモインフォマティクスとどう繋がるのですか」

**説明 (約 30 秒)**:
> 簡単に経歴を補足します。明治薬科大学医療分子解析学研究室の助教を、東京科学大学 (Science Tokyo) 大学院心肺統御麻酔学分野の講師を兼任しております。現在、サバティカルとして MGH に滞在しております。専門は PK/PD モデリング、機械学習による予測モデルおよび画像解析、ケモインフォマティクスです。本研究は、ケモインフォマティクスで標準的な「relevance × diversity × cardinality を同時最適化する手法」を、scRNA-seq の遺伝子選択に応用したものです。MGH 滞在中は、バイオインフォマティクスの中でも特に gene selection に取り組みたいと考えております。

### S2 — Data Preparation Flow (データ整形パイプライン)

**質問例**: 「データはどう前処理しましたか」「8 細胞種はどう決めましたか」「biology filter は何を除外しましたか」

**説明 (約 60 秒)**:
> パイプラインを補足します。元の Azimuth annotation には 30 種以上の subtype がありますが、3 つの基準で 8 集団に集約しました。第一に MS 病態に関与する主要 lymphoid または myeloid 群、第二に各 donor 平均 20 細胞以上で安定 pseudobulk が組める集団、第三に全 cohort で再現性高く annotation される集団です。Biology filter は Heumos 2023 と Luecken & Theis 2019 のベストプラクティスに準拠し、ミトコンドリア、リボソーム、熱ショック、核内 lncRNA、housekeeping を除外しました。ただし RPLP0/1/2 と RPSA は意図的に保持しています。RPLP は specialized ribosome 仮説の対象、RPSA は実体が laminin receptor で BBB 通過に必須のため、MS 文脈では保持すべきと判断しました。pseudobulk 化により統計単位は donor となります。

### S3 — Per-cohort Detailed Results (cohort 別詳細結果)

**質問例**: 「cohort 別にどの程度ばらつきますか」「過学習の心配はないですか」「N=50 で十分ですか」

**説明 (約 60 秒)**:
> Cohort 別の詳細をお示しします。QUBO は Pappalardo 0.807、Heming 0.738、Ramesh 0.819 と、いずれも 0.74 から 0.82 のタイトな範囲に収束します。LASSO は 0.72 から 0.85 と幅広い分布です。
>
> 統計的観点では、per cell type の classifier は K=17 features × 22 events で EPV 1.3 と、classical Peduzzi 1996 の閾値である EPV 10 には届きません。しかし L2 regularization、8-cell-type ensemble、そして cross-cohort validation により過学習を抑制しています。何より、cohort 間 σ_AUC が 0.044 という empirical な安定性が、過学習が深刻でないことの直接的証拠です。Modern guidelines である van Smeden 2019 と Riley 2019 では、prediction 目的かつ regularization 付きの場合、EPV 2-5 でも実用上許容されるとされています。

### S4 — Top Genes per Cell Type & Curated Enrichment Detail (遺伝子詳細)

**質問例**: 「どんな遺伝子が選ばれましたか」「鉄代謝以外の生物学はどうですか」「MS GWAS との関係は」

**説明 (約 60 秒)**:
> Cell type 別の上位選択遺伝子と curated set 詳細をお示しします。Mono では CST3、SAT1、FTL、HLA-DPB1、LYZ、IFI30、FTH1、CD74、TPT1。NK では KLRB1、CCL5、LTB、CRIP1、KLRC1、GNLY。dnT では GZMA、IL32、ISG15、CCL5、TXK。
>
> Curated set enrichment は 49 stable 遺伝子で評価しました。鉄代謝が fold enrichment 36 倍、cytotoxic 18 倍、MHC II 16 倍。全 448 selection で見ると、Reactome MHC II 経路の 10 遺伝子全てが選ばれ、IMSGC MS GWAS の 11 hit のうち 9 個が再現されます。BACH2、CXCR4、IL7R、そして HLA 一族すべてです。B cell では IGHM、IGKC、IGLC2 が選ばれ、これは MS の oligoclonal band と整合的です。

---

## 想定 Q&A

> Q: 細胞数 ≥ 20 の閾値はどう決めましたか？

A: 「Pseudobulk の安定性の経験則です。10 cells 以下では発現平均が不安定になり、cohort 間で再現しにくくなります。20 cells 以上で各 donor の細胞種別 profile が安定する経験的閾値で、Squair 2021 や他の scRNA-seq pseudobulk 論文でも標準的に採用されています」

> Q: D-Wave 量子アニーラでの検証はされていますか？

A: 「現状は古典 Simulated Annealing で完結しています。1 instance あたり約 3 秒という計算速度のため、現時点で量子アニーラを使う必要性は低いです。一方で、データ規模を拡大して全 HVG 約 2,500 遺伝子を一度に最適化するシナリオでは、D-Wave Leap Hybrid Solver の使用が現実的です。Future Work と位置付けています」

> Q: 1 つの細胞の MS likelihood を予測することはできますか？

A: 「はい、Multi-Instance Learning による cell-level 拡張を Phase 2 として計画しています。各 donor を bag、cells を instances として attention 機構で集約することで、cell-level の解釈が可能になります。同時に、QUBO を gene selection だけでなく cell selection にも応用する余地があります。詳細は MIL_design.md にまとめています」

## ペース・話し方のコツ

- **数字 (385,000、0.788、0.044、36 倍、p < 10⁻⁸) の前後に短い間** を置く
- **Wet-lab 用語** (oligoclonal band、iron rim、paramagnetic rim lesion、laminin receptor) は丁寧に発音
- **スライド遷移は短く形式的に**: 「次に method に進みます」「結果に進みます」「最後に conclusion です」
- **質問時の補助スライド呼び出しは自然に**: 「詳細を補助スライドでお示しします」
