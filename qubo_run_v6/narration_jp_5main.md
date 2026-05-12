# 日本語ナレーション — 5 main + 5 supplementary 構成

**プレゼンテーション**: QUBO 最適化による細胞種特異的遺伝子パネルを用いた多発性硬化症の cross-cohort 分類
**発表者**: 浅田 瑞穂 (Mizuho Asada, Ph.D)
**スタイル**: アカデミックな研究発表向け、です・ます調、短文を重ねる
**対象**: MS 専門家集団 (基礎背景は省略、技術的・診断的課題に絞る)
**総所要時間**: 約 5-6 分 (main 5 枚)、質疑応答に応じて補助資料を提示
**日付**: 2026 年 5 月 7 日

---

## Main Slide 1 — Title + Self-introduction  *(≈ 50 秒)*

> 皆さんおはようございます。本日はお時間をいただきありがとうございます。「QUBO 最適化を用いた細胞種特異的遺伝子パネルによる多発性硬化症の cross-cohort 分類」というタイトルで発表させていただく、朝田瑞穂と申します。
>
> まず簡単に自己紹介させてください。明治薬科大学 医療分子解析学研究室の助教、また Institute of Science Tokyo 麻酔学分野の講師を兼任しており、現在は MGH に Visiting Researcher として滞在し、バイオインフォマティクスと最適化手法を用いた biomarker discovery に取り組んでいます。
>
> 普段は、スライド右側の **3 つのテーマ**を軸に研究をしています。
>
> **第一に、ケモインフォマティクスです**。これは化合物構造に基づく活性予測の分野で、薬剤候補スクリーニングでは「**効く** (relevance) × **重複しない** (diversity) × **適切な数** (cardinality)」 を同時最適化する手法が標準的に使われます。**本日ご紹介する QUBO、Quadratic Unconstrained Binary Optimization、二次計画によるバイナリ最適化、はまさにこの考え方を、scRNA-seq の遺伝子選択に持ち込んだもの**です。
>
> **第二に、AI 補助麻酔管理**です。PK/PD モデリングと機械学習を組み合わせた薬効予測モデルを構築しています。
>
> **第三に、臨床 AI 応用**です。医用画像解析と診断支援モデルを開発しています。
>
> 本日の研究は、これらケモインフォマティクスの最適化アプローチを **バイオインフォマティクス**、特に **MS バイオマーカー設計**に応用したもので、MGH 滞在中の主軸テーマの具体例です。それでは本論に進みます。

---

## Main Slide 2 — Introduction (Background + Aim + Dataset)  *(≈ 50 秒)*

> 背景です。
>
> 多発性硬化症の診断は、現在も MRI 画像と CSF 所見に大きく依存しており、**確立された分子・細胞レベルのバイオマーカーは依然として限られています**。Single-cell RNA-seq を用いたバイオマーカー研究は近年急速に拡大していますが、バッチ効果やドナー差により **cohort 間での再現性が不安定**です。そのため、cross-cohort で再現性の高い単一細胞バイオマーカーの枠組みが求められています。
>
> 本研究の目的は明確です。QUBO 最適化を用いて、**biomarker performance** — すなわち独立した複数 cohort 間で再現性のある分類器 — と、**biological validity** — すなわち各細胞種で MS 病態に関連する biologically relevant な遺伝子・経路の同定 — の両方を、同一フレームワーク内で同時に取り出すことです。
>
> データセットは 4 つの公開コホートを統合したもので、50 patient、99 sample、385,000 cells から構成されます。Heming 18 patient は CSF のみで患者あたり 2 sample 程度の longitudinal、Pappalardo 11 patient と Ramesh 17 patient は CSF と PBMC の両方を採取、Touil 4 patient は HD only のため学習データに固定します。次のスライドで具体的な解析手法をお話しします。

---

## Main Slide 3 — Method  *(≈ 75 秒)*

> 次に手法の中身に入ります。
>
> まず QUBO というのは、「各候補遺伝子を選ぶか選ばないか」を 0 か 1 のバイナリで表し、選んだ遺伝子集合の "コスト" を 2 次関数 H(x) で評価して、それを最小化する組合せを探索する最適化フレームワークです。本研究のコスト関数は、以下の 3 項で構成されます。
>
> 遺伝子選択を 3 項からなる目的関数で定義します。
>
> 第一項は relevance の最大化です。各遺伝子のスコアは edgeR の |t-statistic| の二乗で、MS と HD で発現差が大きい遺伝子に高いスコアを与えます。
>
> 第二項は redundancy の最小化です。選択された遺伝子ペアの相関行列を γ で重み付けして罰則化します。同じ経路に属する遺伝子の重複選択を抑制します。これは univariate な既存手法では扱えない部分であり、QUBO の核心的な貢献です。
>
> 第三項は cardinality 制約です。選択数を K に近づける soft penalty で、K は 10, 20, 30 から内側 5-fold CV で自動選択されます。K 自体もデータ駆動で決まる点が重要です。解法は古典的 Simulated Annealing で約 3 秒/instance です。
>
> パイプラインは 5 手法すべてで共通です。各細胞種ごとに DEG を計算し、その手法で K 遺伝子を選び、L2 logistic regression で MS 確率を予測します。各細胞種の分類器が出した MS 確率を単純平均し、いわゆる soft voting で patient 単位の MS / HD 判定を行います。3 cohort × Leave-One-Cohort-Out で評価します。
>
> 5 手法は DE_top, HVG, LASSO, Elastic Net, QUBO です。すべて同じ候補プール、同じ K グリッド、同じ分類器、同じ ensemble を共有しており、違いは選択ロジックそのものだけです。これにより、選択方法そのものの効果を切り分けて評価できます。

---

## Main Slide 4 — Data Preparation Flow  *(≈ 60 秒)*

> パイプラインを 7 ステップで概観します。
>
> Step 1 が入力で、385,000 細胞の統合 scRNA-seq オブジェクトです。
>
> Step 2 では Azimuth で細胞種を annotation し、30 種以上の subtype を 8 つの主要免疫サブセット — B、Mono、CD4_T、CD8_T、NK、DC、dnT、gdT — に集約します。各 donor で 20 cells 以上を確保できる粒度として選定しました。
>
> Step 3 で donor ごとに pseudobulk を構築します。CSF と PBMC の両 compartment で donor × cell type の matrix を作成します。pseudobulk により、統計の単位が donor、すなわち「人」になります。
>
> Step 4 は遺伝子の絞り込みです。HVG selection で 3,000 遺伝子に絞り、Heumos 2023 と Luecken & Theis 2019 のベストプラクティスに従って biology filter を適用します。約 7,960 遺伝子が残ります。
>
> Step 5 で edgeR DEG を per cell type × per holdout で実行し、|t-statistic| 上位 100 遺伝子を候補プールとして QUBO に渡します。
>
> Step 6 が QUBO selection で、K を内側 CV で自動選択します。平均すると cell type あたり 17 遺伝子が選ばれます。
>
> Step 7 で各細胞種ごとに L2 logistic で出した MS 確率を平均、いわゆる soft voting で patient 単位の MS / HD 判定とします。3 cohort で 1 つを held out にする LOCO 交差検証で外部評価を行います。
>

---

## Main Slide 5 — Results 1: Cross-cohort biomarker performance (Table 1)  *(≈ 50 秒)*

> 結果の前半は cross-cohort の biomarker 性能です。
>
> Table 1 は CSF と PBMC それぞれにおける 5 手法の held-out AUC とコホート間 σ_AUC です。QUBO は CSF で AUC 0.788、PBMC で 0.768 と、**いずれの compartment でも highest AUC** を達成しました。CSF の σ_AUC は 0.044 で、LASSO の 0.068 より明らかにタイトで、Elastic Net (σ 0.041) と実質互角です。PBMC は Heming に sample がないため 2 cohort 平均で、CSF との直接の数値比較は留意してください。
>
> 補足ですが、この AUC 0.788 は、各 fold の QUBO panel — cell type あたり約 17 遺伝子 × 5 cell type で計約 85 遺伝子 — を用いて算出されています。

---

## Main Slide 6 — Results 2: Selection stability of QUBO-selected genes (Figure 1)  *(≈ 60 秒)*

> 次に Figure 1 で、QUBO が各細胞種で具体的にどんな遺伝子を選んでいるかをご覧いただきます。
>
> ヒートマップは 8 細胞種すべてを列に並べ、各細胞種の上位 5 遺伝子の union を行に表示しています。色が選択頻度、右側のラベルが MS biology カテゴリです。Mono は iron rim 関連の FTL と MHC II の HLA-DPB1、NK は cytotoxic の KLRB1, KLRC1, CCL5, CRIP1、dnT は cytotoxic と Type I IFN の GZMA, ISG15、gdT は cytotoxic の GZMA, CD69, TPT1、B は plasma cell・IgM heavy chain の IGHM など、いずれも各細胞種固有の biology を反映した panel になっています。
>
> 灰色になっている 3 列、CD4_T、CD8_T、DC は **安定した選択が得られていません**。これは heterogeneous な pseudobulk 集団における signal dilution を反映している可能性があり、最後の future work でお話しする MIL (Multi-Instance Learning) の動機にもなっています。
>
> 右側の緑色の枠は独立外部検証の結果です。Ramesh ら 2020 年 PNAS の MS pathogenic B 細胞シグネチャと比較したところ、**候補プールに残った 13 遺伝子全てが QUBO によって独立に再現**されました。CSF immune dynamics シグネチャも 9/12 を捕捉しています。
>
> 重要なのは、Ramesh らの論文は QUBO に事前知識として一切与えていない、ということです。それでも独立に同じ selections に収束しており、**選択された panel の biological validity を裏付ける**結果になっています。

---

## Main Slide 7 — Conclusion + Future Work  *(≈ 55 秒)*

(参考: 補足として、QUBO 選択遺伝子そのものを panel として per-cell の AUCell でスコア化した検証では、B 細胞 panel で MS と HD の差分 +0.049、q 値 5.7×10⁻¹⁵ と圧倒的な有意差が得られています。詳細は supplementary slide を参照。)


> 最後に結論をまとめます。
>
> 本研究のメインメッセージは、QUBO は biology と biomarker を同時に提供する、ということです。
>
> 第一に biomarker 性能。CSF AUC 0.788、cross-cohort σ 0.044 を達成しました。
>
> 第二に biological validity。Mono は MHC II と iron metabolism、NK・dnT・gdT は cytotoxic axis、B 細胞は **免疫グロブリン関連の program** を捕捉し、**known MS-associated immune programs を再現**しています。
>
> 第三に cross-cohort 再現性。3 cohort 全てで AUC 0.74 から 0.82 のタイトな範囲に収束し、**translational potential** を支える臨床的に意義のある再現性を示しています。
>
> 第四に方法論的新規性です。私たちの知る限り、relevance、非冗長性、cardinality の 3 軸を同時最適化する **細胞種別 gene panel 設計の最初期の枠組みの一つ**です。冗長性を選択ロジックに **明示的に組み込んでいる**点が、DE_top や標準的 sparse selection 手法とは異なります。
>
> 結論として、本研究は MS scRNA-seq データに対して **biomarker reproducibility と biological coherence を同時に評価する unified framework を提供する**ものであり、translational potential と機構的解明の両方に資するものです。
>
> Future Work として Multi-Instance Learning による cell-level 拡張を計画しています。本研究の最大の限界は、CSF で最も細胞数の多い CD4 と CD8 が、pseudobulk dilution と biology filter のため候補プールが枯渇し、ensemble から脱落していることです。CD4 Th17 や CD8 cytotoxic は MS 病態の中核軸ですので、これを取り戻すことが Phase 2 の最優先課題です。MIL は各 donor を bag、cells を instances として、attention 機構により CD4 内の Th17 や CD8 内の Trm のような病態駆動 subset を浮上させます。pseudobulk なしで cell-level 解像度を保つので、dilution の問題が原理的に解消します。QUBO は gene selection に加え、informative cells の coreset 選択という新たな役割も担い得ると考えています。
>
> ご清聴ありがとうございました。質疑をお受けいたします。

---

## 補助スライドの使い方 (Q&A 時)

### S1 — Per-cohort Detailed Results (CSF) & EPV 議論

**質問例**: 「cohort 別にどの程度ばらつきますか」「F1, MCC, AP の詳細は」「過学習の心配は」

**説明 (約 60 秒)**:
> Cohort 別と全 metrics の詳細です。QUBO は Pappalardo 0.807、Heming 0.738、Ramesh 0.819 と、いずれも 0.74 から 0.82 のタイトな範囲に収束します。LASSO は 0.72 から 0.85 と幅広い分布です。
>
> 全 metrics で見ると、QUBO は AUC、F1、MCC の主要 3 指標で 1 位、AP は EN に Δ 0.024 で僅差負け、σ_AUC は EN に Δ 0.003 で実質互角です。
>
> 統計的観点では、per cell type の classifier では K=17 features × 22 events で EPV 1.3 と classical Peduzzi 1996 の EPV 10 ルールには届きません。しかし L2 regularization、cell type ensemble、cross-cohort validation で過学習を抑制しています。何より cohort 間 σ_AUC = 0.044 という empirical な安定性が、過学習が深刻でないことの直接的証拠です。

### S2 — Top Genes per Cell Type & Curated Enrichment

**質問例**: 「具体的にどんな遺伝子が選ばれましたか」「鉄代謝以外の生物学は」

**説明 (約 60 秒)**:
> Cell type 別の上位選択遺伝子と、curated set enrichment の詳細です。Mono では CST3、SAT1、FTL、HLA-DPB1、LYZ、IFI30、FTH1、CD74。NK では KLRB1、CCL5、LTB、CRIP1、KLRC1、GNLY。dnT では GZMA、IL32、ISG15、CCL5。
>
> Curated set enrichment では、Hametner 2013 の鉄代謝セットが fold enrichment 36 倍と q 値 2×10⁻³ で最も顕著です。Cytotoxic effectors 18 倍、MHC II 経路 16 倍と、いずれも有意な濃縮を示します。

### S3 — Figure 1 の 23 遺伝子 機能アノテーション

**質問例**: 「Figure 1 の遺伝子それぞれの機能は?」「どれが MS 関連?」

**説明 (約 60 秒)**:
> Figure 1 の y 軸 23 遺伝子それぞれの biology を細胞種別に表にまとめたサプリです。緑色マーカーが文献上 MS 関連としてキュレートされた遺伝子で、IGHM (B 細胞: plasma marker・CSF oligoclonal band の元)、FTL (Mono: iron rim)、HLA-DPB1 (Mono: MHC II・MS GWAS 最強リスク座位)、KLRB1・KLRC1・CCL5 (NK: cytotoxic axis)、GZMA (dnT/gdT: cytotoxic granule)、ISG15 (dnT: Type I IFN signature) などが該当します。MS 病態軸とのマッピングでは、Iron metabolism = FTL、MHC II/GWAS = HLA-DPB1、Cytotoxicity = GZMA・KLRB1・KLRC1・CCL5、Type I IFN = ISG15、Plasma/Ig = IGHM、Tissue residency = CD69、Th1 axis = TXK と、MS pathology の中核軸が QUBO panel に網羅されています。

### S4 — MS Biology Reference Map

**質問例**: 「文献的に重要な遺伝子と比較して、QUBO は何を捉えていますか」「どこが拾えていませんか」

**説明 (約 60 秒)**:
> 文献的に MS で重要とされる cell type × gene combinations と、QUBO 選択の対応マップです。
>
> 強く再現できているのは 5 軸: MHC II 抗原提示、iron rim biology、cytotoxic axis、Type I IFN、MS GWAS hits です。
>
> 弱いのは B cell oligoclonal — 各 donor で異なるクローンが増殖するため共通シグナルとしては抽出しにくい性質によります。CD8 cytotoxic は CSF panel 数の制約です。
>
> 拾えていないのは Th17/Treg axis (CD4) と MS DMT 標的 (CD20, ITGA4) です。CD20 は薬剤標的として優れますが、MS と HD で発現差がないため候補プールに入りません。これは「薬剤標的 ≠ 発現バイオマーカー」の典型例で、方法論の限界ではなく正常な挙動です。

### S5 — PBMC AUCell — cell-type × MS gene set 活性

**質問例**: 「PBMC でも MS biology が見えますか」「CSF と PBMC でどう違いますか」

**説明 (約 45 秒)**:
> PBMC で文献由来 MS 関連 gene set 7 種類を 8 細胞種ごとに per-cell AUCell スコア化したものです。CSF と比べると効果サイズは全般的に小さくなっており、Mono の MHC II 経路は同様に活性化していますが、Type I IFN や cytotoxic axis のシグナルは CSF より弱い。これは CSF が病巣に近接していて MS 病態シグナルをより直接捕捉していることの裏付けです。一方で PBMC は腰椎穿刺なしに採取できる利点があるので、screening 段階での層別化や経時モニタリングへの活用余地があります。

### S6 — QUBO 定式化 詳細

**質問例**: 「QUBO の数式を詳しく教えてください」「Simulated Annealing は何ですか」

**説明 (約 60 秒)**:
> QUBO の詳細定式化です。本研究の cost function H(x) は 3 項からなります。第一項 -sᵀx は relevance を最大化する項で、s は edgeR の |t-statistic| の二乗です。第二項 γxᵀRx は redundancy 最小化で、R は遺伝子ペアの相関行列、γ で罰則強度を制御します。第三項 λ(Σx − K)² は cardinality 制約で、選択数を K に近づける soft penalty です。K は {10, 20, 30}、γ と λ も同時に内側 5-fold CV で自動選択しています。解法は古典的 Simulated Annealing で約 3 秒/instance、量子アニーラへの移植可能性も担保していますが、本研究では古典 solver で完結しています。

### S7 — 細胞種カバレッジと限界 (CD4/CD8 dilution → MIL の必然)

**質問例**: 「CD4 と CD8 が空なのは?」「MIL って具体的にどう機能?」

**説明 (約 60 秒)**:
> CD4 と CD8 が CSF で QUBO 選択 0 になる原因の詳細です。両者は CSF で最も細胞数の多い細胞種で、本研究データでも CD4 は 98,000 細胞、CD8 は 26,000 細胞ありますが、機能サブタイプの混合体です。CD4 は Th1, Th17, Treg, Tfh を内包し、CD8 は Tem, Tcm, Trm, exhausted, MAIT を含みます。MS で異常を示すのは Th17 や CD8 effector など数 % のサブセットだけで、pseudobulk で donor 内全細胞を平均すると、その signal が 90% の通常 T 細胞に薄められて消えます。これが pseudobulk dilution の問題で、edgeR + biology filter の段階で候補遺伝子が CD4 で 1 個、CD8 で 0 個まで落ち、QC で除外されました。Phase 2 の MIL はまさにこの問題を解消する拡張で、各 donor を bag、cells を instances として、attention 機構で pseudobulk なしに病態駆動 subset を浮上させます。

### S8 — QUBO panel の per-cell AUCell 検証

**質問例**: 「QUBO panel は単一細胞でも MS を識別できますか」

**説明 (約 45 秒)**:
> 補足として、QUBO 選択遺伝子そのものを panel として per-cell の AUCell でスコア化した検証もしています。最大の所見は B 細胞 panel で、MS と HD の差分が +0.049、q 値 5.7×10⁻¹⁵ と圧倒的に有意でした。これは QUBO の B 細胞選択 30 遺伝子が、pseudobulk レベルだけでなく cell-level でも MS を判別できる直接証拠です。Mono panel と NK panel も同様に有意で、本研究の panel が単一細胞レベルでも biology を捉えていることを示しています。

### S9 — 各コホートの患者背景 (n=50 donors)

**質問例**: 「コホートごとの性別・年齢バランスは?」「MS と HD で年齢は揃っていますか?」「EDSS や DMT 治療歴は?」

**説明 (約 60 秒)**:
> 4 つのコホートの患者背景を表にまとめたものです。全体で 50 donors、MS 28 例、HD 22 例。Pappalardo は MS/HD バランス良く、平均 30 歳前後で最も若い集団です。Heming は MS/HD 各 9 例で完全 balance、年齢は 30 代中盤で典型的な MS 発症年齢層です。Ramesh は MS 14 例 + HD 3 例と MS-rich で、年齢は最も高めの 39 歳。Touil は HD のみ 4 例 (cryopreserved CSF) で、訓練データの HD 数を補強する役割です。
>
> 重要な点として、**MS と HD で女性比率が共に 64% で揃い、年齢も 37.2 歳 vs 34.1 歳と well-matched** です。性別・年齢の confounder は最小限に抑えられています。
>
> EDSS、disease duration、DMT 治療歴、MRI 所見など詳細な臨床情報は、4 つのコホートすべての公開メタデータから一律に取得することができなかったため非掲載となっています。これらは Phase 2 で連携施設経由での取得を計画しています。

### S10 — パイプラインの概念図

**質問例**: 「全体の流れをもう一度説明してもらえますか」「細胞種ごとに独立に学習しているの?」

**説明 (約 45 秒)**:
> 全体パイプラインの概念図です。3 つのセクションで成り立っています。第一に input、scRNA-seq の各細胞を Azimuth で B、Mono、NK、dnT、gdT などに annotation します。第二に細胞種別の QUBO 選択と分類器、ここが核心で、細胞種ごとに**独立に** QUBO が 10 から 30 個の遺伝子を選び、その panel を使って L2 logistic regression で MS 確率を出します。例えば B 細胞 panel から 0.72、Mono panel から 0.68、NK panel から 0.81 というふうに細胞種ごとに値が出ます。第三に ensemble、5 細胞種の MS 確率を単純平均する soft voting で patient 単位の MS / HD 判定とします。例えばこの図では平均が 0.65 で 0.5 の閾値を超えているので MS と判定する流れです。

---

## ペース・話し方のコツ

- **数字 (385,000、0.788、0.044、36 倍、p < 10⁻⁸) の前後に短い間** を置く
- **wet-lab 用語** (oligoclonal band、iron rim、paramagnetic rim lesion、laminin receptor) は丁寧に発音
- **スライド遷移は短く形式的に**: 「次に method に進みます」「結果に進みます」「最後に結論をまとめます」
- **質問で補助資料を出す必要があれば**: 「詳細は別の図でお示しできます」と前置き

## 想定される質問

> Q: 49 個で AUC を計算したのですか？

A: いいえ。**AUC 0.788 は per-fold の QUBO panel** — cell type あたり約 17 遺伝子 × 8 cell type で計 約 136 遺伝子 — で算出しています。stable 49 は全 panels の union 448 のうち再現性高く選ばれる core で、生物学的解釈の対象です。両者は別物です。

> Q: K の最適値はどう決めましたか。なぜ K=10 固定ではないのですか。

A: K は内側 5-fold CV で {10, 20, 30} から自動選択しています。Cell type ごとに最適 K が異なるためです。例えば B cell は IGH/IGL ファミリの多様性が高く K=20 から 30 が選ばれやすく、Mono は MHC II が支配的で K=10 で十分なケースが多い、といった具合です。K も同じデータ駆動で決めることで、cell type 別の biology に沿った panel size を実現します。実測平均は約 17 遺伝子/cell type です。

> Q: なぜ CSF で CD4 と CD8 が ensemble に入っていないのですか。これは MS と無関係という意味ですか。

A: いいえ、全く逆です。CD4 と CD8 は CSF で最も細胞数の多い細胞種で、本研究データでも CD4 が 98,000 細胞、CD8 が 26,000 細胞あります。ところが pseudobulk + edgeR + biology filter の段階で、候補遺伝子が CD4 で 1 個、CD8 で 0 個まで落ち込み、QC で除外されました。原因は、CD4 と CD8 が機能サブタイプの混合体であることです。CD4 は Th1, Th17, Treg, Tfh を内包し、CD8 は Tem, Tcm, Trm, exhausted, MAIT を内包します。MS で異常を示すのは Th17 や CD8 effector のような数 % のサブセットだけで、pseudobulk で donor 内全細胞を平均すると、その signal が 90% の通常 T 細胞に薄められて消えてしまいます。これは方法論的限界であり、CD4 Th17 や CD8 cytotoxic が MS 病態の中核軸であることは Ramesh 2020 や IL7R, IL2RA の MS GWAS hit からも明らかです。Phase 2 の MIL がまさにこの問題を直接解消するための拡張で、attention 機構で pseudobulk なしに病態駆動 subset を浮上させます。

> Q: では現状の "5 cell type ensemble" でも AUC 0.788 が出ているのは妥当ですか。

A: はい、妥当な解釈は二つあります。第一に、B 細胞の oligoclonal、Mono の MHC II と iron rim、NK の cytotoxicity という MS 病態の "壁" の biology は本 panel で十分に捕捉されており、これだけでも 3 cohort 全てで AUC 0.74-0.82 のタイトな範囲に収束します。第二に、CD4 Th17 や CD8 cytotoxic という "中心メカニズム" の signal が加われば、AUC は更に向上する余地があり、これが Phase 2 の改善余地として残っています。逆に言えば、現状でも壁の biology だけで cross-cohort 安定性を達成できているのは、QUBO の relevance × 非冗長性 × cardinality 同時最適化が効いている証拠とも解釈できます。

> Q: cell-level の予測はできますか？

A: はい、Multi-Instance Learning による cell-level 拡張を Phase 2 として計画しています。詳細は MIL_design.md にまとめております。

> Q: なぜ Touil cohort (HD only, 4 例) を train 専用に固定したのですか?

A: HD の **control diversity を増やし**、cohort-specific な技術的偽陽性に対する **robustness を向上**させるためです。Touil は cryopreserved CSF という処理プロトコルが他 3 cohort と異なるため、test に含めると "cryopreservation artifact" を MS シグナルと誤学習するリスクがあります。train に固定することで、HD 側の technical variation を学習に組み込み、MS-specific signal の robustness を高めています。実際 σ_AUC = 0.044 という cross-cohort 安定性に貢献していると考えられます。

> Q: "biomarkers remain limited" の "limited" は何を指していますか？候補は多数ありますが。

A: 「候補が多数あること」と「**clinically established / reproducible** なものがあること」は別の話です。本研究で使う表現は "established molecular or cellular biomarkers remain limited" であり、候補の存在を否定しているわけではなく、**cross-cohort で再現する形で establish されたものが限られている**ことを指しています。これが本研究の問題設定の出発点です。
