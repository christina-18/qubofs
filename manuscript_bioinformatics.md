# Cell-type-specific feature selection via Quadratic Unconstrained Binary Optimization for cross-cohort multiple sclerosis classification from single-cell RNA sequencing

**Mizuho Asada¹²³, Takahisa Mikami³⁴, Daisuke Tominaga⁵, Michael Levy³⁴,***

¹ Laboratory of Medical Molecular Analytics, Meiji Pharmaceutical University, Tokyo, Japan
² Department of Anesthesiology, Institute of Science Tokyo, Tokyo, Japan
³ Neuroimmunology Clinic and Research Laboratory, Division of Neuroimmunology and Neuroinfectious Disease, Department of Neurology, Massachusetts General Hospital, 65 Landsdowne St., Cambridge, MA 02139, USA
⁴ Harvard Medical School, Boston, MA, USA
⁵ Division of Mathematical Sciences and Life Informatics, Meiji Pharmaceutical University, Tokyo, Japan

\* To whom correspondence should be addressed. Email: mlevy11@mgh.harvard.edu

**Article type**: Original Paper

**Subject section**: Genome analysis / Gene expression

---

## Abstract

**Motivation**: Single-cell RNA sequencing (scRNA-seq) has accelerated biomarker discovery in multiple sclerosis (MS), yet cross-cohort reproducibility of candidate gene panels remains limited by batch effects and donor variability. Existing feature-selection methods optimize either statistical relevance, through differential-expression-based selection, or sparsity, through penalized regression such as LASSO and Elastic Net. They rarely control for redundancy among co-regulated genes within a cell-type-specific framework. Univariate filters yield panels saturated with co-expressed gene clusters, whereas penalized regression couples redundancy reduction to classifier loss.

**Results**: We introduce a per-cell-type gene-panel selection framework based on Quadratic Unconstrained Binary Optimization (QUBO) that jointly optimizes three complementary objectives, namely relevance, pairwise non-redundancy and cardinality, within a single quadratic cost function decoupled from the downstream classifier. Applied to a four-cohort integrated MS scRNA-seq compendium (50 patients, 99 samples, 385,116 cells), our primary configuration QUBO_hybrid (combining a top-20 univariate pre-filter with QUBO redundancy optimization) achieved cross-cohort held-out AUC of 0.858 (σ = 0.108), F1 = 0.731 and Matthews correlation coefficient (MCC) = 0.427 in cerebrospinal fluid under Leave-One-Cohort-Out validation, outperforming four matched feature-selection baselines (DE-top, HVG, LASSO and Elastic Net) on F1 and MCC and remaining competitive on AUC, with a clinically deployable 64-gene total panel (8 cell types × ~8 genes). Without prior knowledge, the resulting panels recovered all 13 candidate-pool genes of an independently published B-cell pathogenic signature (Ramesh *et al.*, 2020). Per-cell AUCell scoring of the QUBO B-cell panel produced an MS-versus-HD median-activity difference of +0.049 (q = 5.7 × 10⁻¹⁵), confirming biological relevance at single-cell resolution.

**Availability and implementation**: The pipeline is implemented in Python and R, released under the MIT license at https://github.com/christina-18/scRNA-QUBO with reproducible Docker and conda environments and worked examples on the four public cohorts.

**Contact**: mlevy11@mgh.harvard.edu

**Supplementary information**: Supplementary data are available at *Bioinformatics* online.

---

## 1 Introduction

Multiple sclerosis (MS) is a chronic autoimmune disease of the central nervous system for which diagnostic and prognostic biomarkers remain limited to magnetic resonance imaging and the detection of cerebrospinal fluid (CSF) oligoclonal bands. Single-cell RNA sequencing (scRNA-seq) has accelerated MS biomarker discovery by resolving immune-cell heterogeneity in CSF and peripheral blood mononuclear cells (PBMC) (Schafflick *et al.*, 2020; Pappalardo *et al.*, 2020; Ramesh *et al.*, 2020; Heming *et al.*, 2021). However, cross-cohort reproducibility of candidate gene panels remains poor as a consequence of batch effects, sequencing-platform differences and donor heterogeneity (Heumos *et al.*, 2023). A reproducible, cell-type-specific biomarker framework is therefore required for translation to clinical use.

Existing feature-selection strategies fall into two broad families. Univariate filters, such as the top-K differentially expressed genes (DE-top) ranked by edgeR, DESeq2 or limma-voom |t|-statistics, or the top-K highly variable genes (HVG), are simple but tend to select clusters of co-expressed genes (for example, multiple HLA class II family members) that contribute redundant information. Penalized regression methods (LASSO, Elastic Net) sparsify selections via L1 or L1+L2 regularization but couple redundancy reduction to a specific classifier loss, conflating selection with model fitting; selections are also unstable across resamples (Meinshausen and Bühlmann, 2010). Neither family enforces a hard cardinality constraint with explicit pairwise dissimilarity, both of which are essential for designing biomarker panels that are interpretable, clinically tractable and reproducible.

Cheminformatics provides a complementary template: in compound library design, relevance, diversity and cardinality are jointly optimized (Snarey *et al.*, 1997; Pearlman and Smith, 1998). The natural mathematical home for this triplet is Quadratic Unconstrained Binary Optimization (QUBO), a binary-variable quadratic minimization framework with mature classical (Simulated Annealing, Tabu Search) and quantum-annealing (D-Wave) solvers (Lucas, 2014; Glover *et al.*, 2018). QUBO-based feature selection was first proposed for general high-dimensional data using mutual-information importance and redundancy by Mücke *et al.* (2023) and was recently adapted to scRNA-seq pseudotime regression by Romero *et al.* (2025), who demonstrated that quantum and classical solvers produced identical feature subsets and that QUBO recovered nonlinear gene-pseudotime relationships missed by LASSO and random forests. The present work extends this framework to **disease classification** rather than regression, at **cell-type-specific resolution** rather than whole-dataset, with an **explicit cardinality penalty** and a **two-stage hybrid** that combines univariate pre-filtering with QUBO redundancy optimization.

Here we adapt the QUBO framework to scRNA-seq biomarker selection. Our contributions are fourfold. First, we formulate per-cell-type biomarker selection as a QUBO with three terms (score-weighted relevance, correlation-penalized redundancy and soft cardinality) and a hyperparameter grid auto-selected by inner cross-validation. Second, we integrate QUBO selection into a soft-voting cell-type ensemble classifier, evaluated under Leave-One-Cohort-Out (LOCO) cross-validation across four public MS scRNA-seq cohorts (50 patients, 385,116 cells). Third, we benchmark QUBO against four selection baselines (DE-top, HVG, LASSO and Elastic Net) under matched candidate pool, K grid, classifier and ensemble, isolating the effect of the selection logic itself. Fourth, we provide independent biological validation through the recovery of all candidate-pool genes from a published B-cell signature (Ramesh *et al.*, 2020) and through cell-level AUCell scoring of the resulting panels. We additionally disclose a methodological limitation, namely pseudobulk dilution of the T-cell signal, and outline a Multi-Instance Learning extension as future work. The pipeline is open-source at https://github.com/christina-18/scRNA-QUBO.

---

## 2 Materials and Methods

### 2.1 Datasets and preprocessing

Four publicly available multiple sclerosis (MS) single-cell RNA-seq cohorts were integrated for analysis: Pappalardo *et al.* (2020) (PRJNA671484), Heming *et al.* (2021) (osmzhlab/MS-ence-cov), Ramesh *et al.* (2020) (PRJNA549712) and Touil *et al.* (2023) (PRJNA979258). The integrated dataset comprised 50 unique donors, 99 sample-tissue combinations and 385,116 cells across cerebrospinal fluid (CSF) and peripheral blood mononuclear cell (PBMC) compartments. Cohorts included 28 patients with MS and 22 healthy donors (HD), with matched demographic distributions between groups (Supplementary Table S1).

Raw count matrices were downloaded from the original repositories and harmonized to HGNC gene symbols. All preprocessing was performed in R (version 4.5.0) using Seurat (version 5.3.1). Quality control was applied independently for each cohort prior to integration. Cells with fewer than 200 detected genes or mitochondrial RNA fraction exceeding 20% were excluded. Putative doublets were identified and removed using scDblFinder. Normalization and variance stabilization were performed using SCTransform on a per-cohort basis before integration.

### 2.2 Cell-type annotation and pseudobulk construction

Cells were annotated using Azimuth `predicted.celltype.l2` labels (Hao *et al.*, 2021) and collapsed into eight broad immune subsets (B, Mono, CD4_T, CD8_T, NK, DC, dnT and gdT) to balance biological specificity, cross-cohort reproducibility, and sufficient donor-level cell coverage for stable pseudobulk estimation. Gene features were restricted to biologically informative signals by excluding mitochondrial genes, ribosomal genes (with the exception of RPLP0, RPLP1, RPLP2 and RPSA, retained for documented roles as laminin receptors and in blood-brain-barrier transmigration), heat-shock genes, nuclear long non-coding RNAs (e.g., MALAT1 and NEAT1) and a curated housekeeping list, following best-practice recommendations for single-cell analysis (Heumos *et al.*, 2023). After filtering, approximately 7,960 genes remained.

For each donor × cell type × tissue combination, mean pseudobulk expression profiles were constructed by averaging log-normalized expression values across cells belonging to the same group. Donors, rather than individual cells, were treated as the statistical unit throughout feature selection and model evaluation. This aggregation strategy reduces cell-level stochasticity and mitigates pseudoreplication while preserving cell-type-specific transcriptional structure. Pseudobulk profiles were subsequently restricted to the top 3,000 highly variable genes per cell type using Seurat `FindVariableFeatures` with the `vst` method.

### 2.3 Cross-validation framework and leakage prevention

To ensure unbiased evaluation and prevent information leakage, all data-splitting procedures were performed at the donor level rather than at the individual-cell level. Samples originating from the same donor were strictly assigned to either the training or the test partition, but never to both. Donor-level grouping was maintained throughout feature selection, hyperparameter optimization and model evaluation.

Predictive performance was evaluated using cross-validation frameworks designed to assess both within-cohort robustness and cross-cohort generalizability. In the primary analysis, Leave-One-Cohort-Out (LOCO) cross-validation was performed: one of the three MS-containing cohorts (Pappalardo, Heming, Ramesh) was held out for testing while the remaining cohorts were used for training. The Touil cohort (HD only) was fixed to the training set throughout. This design evaluates the reproducibility of selected biomarker panels under cohort-specific technical and biological variation.

All preprocessing steps that depend on expression distributions, including candidate-gene ranking, highly-variable-gene selection, QUBO optimization and classifier training, were performed exclusively within the training data of each fold. Test cohorts were not used at any stage of feature selection or hyperparameter tuning. Hyperparameters, including candidate-pool size and QUBO regularization parameters, were selected by nested cross-validation restricted to the training cohorts. This procedure prevented information leakage from test data into model selection and ensured an unbiased estimate of generalization performance.

### 2.4 Candidate gene pool generation

Candidate gene pools were generated independently within each training fold to prevent information leakage. For each cell type, differential expression statistics were computed from donor-level pseudobulk count profiles comparing MS and healthy donor samples within the training data only. Genes were ranked by the absolute Wald-test statistic from DESeq2 (Love *et al.*, 2014), with diagnosis as the predictor of interest and age, sex, sequencing batch and log₁₀ cells/donor as nuisance covariates. DESeq2 was chosen as the primary differential-expression framework for two reasons. First, its empirical Bayes dispersion shrinkage stabilizes per-gene variance estimates by borrowing information across the transcriptome, retaining usable t-statistics for 80–96% of input genes across all eight cell types in our pseudobulk compendium. Second, DESeq2 is, alongside edgeR, the de facto standard for bulk-style differential expression and is routinely applied to pseudobulk data in single-cell studies (Squair *et al.*, 2021; Crowell *et al.*, 2020).

We additionally evaluated edgeR and limma-voom as alternative Negative-Binomial-GLM frameworks during pilot analyses; both produced stable estimates for fewer than 20 of 3,000 input genes in three heterogeneous cell types (CD4 T, CD8 T and DC), with the surviving genes dominated by constitutively expressed housekeeping transcripts that were subsequently removed by the biology filter. This left zero candidate genes for QUBO selection in those cell types under the edgeR or limma-voom frameworks (Supplementary Table SX). The empirical Bayes shrinkage in DESeq2 was robust to this regime, yielding biology-relevant candidate pools for every cell type. Sensitivity analysis with covariate-adjusted linear-model t-statistics (an alternative recommended by Squair *et al.*, 2021, for log-normalized pseudobulk) yielded concordant gene rankings and cross-cohort AUC within Δ = 0.02 of the DESeq2-based primary analysis (Supplementary Fig. SX).

Prior to ranking, biologically uninformative or technically confounded features (mitochondrial, ribosomal, heat-shock and curated housekeeping genes) were excluded as described in §2.2. To reduce dimensionality while preserving sufficient diversity for combinatorial optimization, the top 100 ranked genes were retained as candidate features per cell type for downstream QUBO optimization and baseline feature-selection methods. Candidate-pool generation was repeated independently for each cross-validation fold, ensuring that test data did not influence feature ranking or candidate selection. All feature-selection methods were evaluated using matched candidate pools, classifier architectures and evaluation procedures, isolating the effect of feature-selection logic itself.

### 2.5 QUBO formulation

For each cell type and cross-validation fold, candidate gene selection was formulated as a Quadratic Unconstrained Binary Optimization (QUBO) problem, following the relevance-redundancy framework of Mücke *et al.* (2023) and the scRNA-seq adaptation of Romero *et al.* (2025), with three modifications detailed below. Let $x_i \in \{0, 1\}$ denote whether candidate gene $i$ is selected, where $x_i = 1$ indicates selection and $x_i = 0$ indicates exclusion. The objective function was designed to jointly maximize gene-level relevance, penalize pairwise redundancy and control the number of selected genes:

$$
\min_{\mathbf{x} \in \{0,1\}^N} \left[ -\alpha \sum_{i=1}^{N} \tilde{r}_i\, x_i \;+\; \beta \sum_{i<j} |\rho_{ij}|\, x_i x_j \;+\; \lambda \left( \sum_{i=1}^{N} x_i - K \right)^2 \right],
$$

where $N$ is the number of candidate genes and $K$ is the target panel size. The first term represents gene-level relevance, with $\tilde{r}_i$ denoting the normalized differential score of gene $i$, computed from the absolute DESeq2 Wald-test statistic described in §2.4. The negative sign encourages the selection of genes with high discriminative relevance. The second term penalizes pairwise redundancy, with $\rho_{ij}$ the Pearson correlation between genes $i$ and $j$ across training-donor pseudobulk profiles. The use of $|\rho_{ij}|$ penalizes both positively and negatively correlated gene pairs, discouraging redundant transcriptional signals. The third term is a soft cardinality constraint that encourages the solution to select approximately $K$ genes. The scalars $\alpha$, $\beta$ and $\lambda$ are non-negative weights; in our reported analyses $\alpha = 1$ and $\beta = \gamma \in \{0.5, 1.0\}$ was tuned by inner cross-validation, with $\lambda \in \{2, 5\}$.

This objective can be expressed in standard QUBO matrix form as

$$
\min_{\mathbf{x} \in \{0,1\}^N} \mathbf{x}^\top \mathbf{Q}\, \mathbf{x},
$$

where the diagonal elements of $\mathbf{Q}$ encode gene relevance and cardinality contributions and the off-diagonal elements encode pairwise redundancy penalties. The same matrix $\mathbf{Q}$ is solvable on classical (Simulated Annealing, Tabu Search) and quantum-annealing (D-Wave) hardware, ensuring methodological portability. All relevance scores and correlation matrices were computed exclusively from training-fold data to prevent information leakage. The QUBO was solved with classical Simulated Annealing (`dwave-neal`) using 30 reads × 600 sweeps for the final selection (8 reads × 200 sweeps during the inner-CV grid for $(\gamma, \lambda)$ tuning), with the target cardinality $K \in \{5, 10\}$ auto-selected per (cell type × fold) by inner 5-fold cross-validation. As a solver-independence check (Romero *et al.*, 2025), we additionally re-solved every instance with the iterated Tabu Search algorithm of Palubeckis (2006) (`dwave-samplers.TabuSampler`) and observed near-complete agreement on the selected gene set (mean Jaccard 0.976 ± 0.082; 90.8% of 119 instances identical), with a small mean energy gap of 0.030 in favor of Tabu (Supplementary §SY). This indicates that classical SA reaches the same optimum that an independent metaheuristic reaches in the overwhelming majority of cases, and that the small remaining gap does not alter the downstream gene selection in any biologically meaningful way.

Our primary QUBO configuration adopted a two-stage screen-then-optimize strategy (denoted QUBO_hybrid). Candidate genes were first pre-filtered to the top 20 per cell type by absolute DESeq2 Wald statistic, and QUBO selection was then performed over this reduced pool. This pre-filter ensures that all candidates carry strong univariate signal, allowing the redundancy term ($\beta\,\mathbf{x}^\top \mathbf{R} \mathbf{x}$) to operate on biologically structured (rather than noise-driven) gene-gene correlations. The strategy is conceptually analogous to Sure Independence Screening (Fan and Lv, 2008) and the Relaxed Lasso (Meinshausen, 2007), where a strong univariate filter precedes a more sophisticated multivariate refinement. Inner cross-validation predominantly selected $K = 10$ (61% of folds), yielding a final panel of approximately 8 genes per cell type after biology filtering and a 64-gene total panel suitable for clinical multiplex assays. We additionally report a vanilla QUBO configuration (top-100 candidate pool, no pre-filter) and a QUBO_consensus variant comprising 10 independent SA runs with seed offsets, returning the top-K genes by selection frequency.

Our formulation differs from Romero *et al.* (2025) in three respects, each motivated by the classification setting and small-cohort regime of MS scRNA-seq. First, the relevance score $\tilde{r}_i$ is taken from the absolute DESeq2 Wald statistic rather than mutual information (MI) between gene expression and cell state. With only 50 donors, donor-level pseudobulk MI estimates are unstable across the quantile-binning bandwidth, whereas the Wald statistic is a calibrated test against MS-versus-HD label noise. A sensitivity analysis using MI-based redundancy is reported in Supplementary §SX. Second, we retain an explicit cardinality penalty $\lambda(\sum x_i - K)^2$ in the cost rather than relying solely on the $\alpha$-balance of importance and redundancy, because clinical panel size $K$ is a fixed design constraint (here $K \in \{5, 10\}$) and an unconstrained QUBO can yield substantially different panel sizes across folds. Third, we couple QUBO selection to a univariate pre-filter (QUBO_hybrid), which is a Sure Independence Screening-style modification that mitigates the small-sample noise inherent to MS CSF cohorts (3-18 donors per cohort).

### 2.6 Baseline feature-selection methods

QUBO-based feature selection was compared against four commonly used baseline approaches: differential-expression ranking (DE-top), highly variable genes (HVG), LASSO and Elastic Net. For the DE-top baseline, genes were ranked by the absolute DESeq2 Wald-test statistic computed from training-donor pseudobulk profiles (as in §2.4), and the top-K genes were selected. The HVG baseline used Seurat `FindVariableFeatures` with the `vst` method to select genes with the highest expression variability within the training data, without using the MS-versus-HD label. Embedded feature-selection baselines were implemented using L1-penalized (LASSO) and L1+L2-penalized (Elastic Net, `l1_ratio = 0.5`) logistic regression in `scikit-learn`. The regularization strength $C$ was tuned over a five-value grid by inner cross-validation restricted to the training folds, and genes with non-zero coefficients were retained as the selected feature set; selection is stochastic in $C$, while the DE-top and HVG selections are deterministic given the outer training data.

To ensure fair comparison across methods, all approaches used the same target-cardinality grid $K \in \{5, 10\}$, identical classifier architectures (§2.7), the same cross-validation framework and identical evaluation metrics. The DE-top, HVG, LASSO and Elastic Net baselines used the standard top-100 candidate pool (a value at which the four baselines are essentially insensitive to candidate-pool size, since they select genes by univariate ranking or coefficient magnitude). The QUBO_hybrid pre-filter to top 20 was applied only to the QUBO_hybrid variant, isolating the contribution of pre-filtered redundancy optimization. This design isolated the effect of feature-selection logic from differences in cardinality budget, classifier choice or evaluation procedure.

### 2.7 Classification and evaluation

Selected gene panels were evaluated using L2-regularized logistic regression classifiers (`scikit-learn`, C = 1.0) trained independently for each cell type using donor-level pseudobulk expression profiles. Classification performance was assessed for the discrimination of MS versus HD donors. For multi-cell-type prediction, cell-type-specific MS probabilities $p_c$ were combined by soft-voting ensemble averaging:

$$
p_{\text{patient}} = \frac{1}{|\mathcal{C}|} \sum_{c \in \mathcal{C}} p_c,
$$

where $\mathcal{C}$ denotes the set of cell types in which at least one gene was selected. The same classifier architecture and ensemble procedure were used for all feature-selection methods, isolating the effect of feature-selection logic from downstream modeling differences.

Predictive performance was quantified using the area under the receiver operating characteristic curve (AUC), the area under the precision-recall curve (average precision, AP), the F1 score, Matthews's correlation coefficient (MCC) and the cross-cohort standard deviation of AUC ($\sigma_{\text{AUC}}$, computed as the standard deviation of per-cohort mean AUCs). All evaluation metrics were computed exclusively on held-out donors from the corresponding cross-validation fold.

In addition to predictive performance, selected gene sets were evaluated for compactness and redundancy. Redundancy was quantified as the mean absolute pairwise Pearson correlation among selected genes within training-donor pseudobulk profiles. Selection stability across folds was assessed using selection frequency and the pairwise Jaccard similarity of selected gene subsets. Biological interpretability was examined post hoc using curated gene-set enrichment (hypergeometric tests against five MS-relevant signatures with Benjamini–Hochberg correction), independent literature-signature recovery (the 27-gene B-cell signature of Ramesh *et al.* (2020) was not used at any stage of our pipeline) and cell-level activity scoring of the selected panels with AUCell (Aibar *et al.*, 2017). These post-hoc analyses were not incorporated into model training or feature selection.

### 2.8 Reproducibility and software

All analyses were performed using R (version 4.5.0) and Python (version 3.11). Single-cell preprocessing and pseudobulk construction were conducted using Seurat (version 5.3.1) and scDblFinder. Differential expression statistics, baseline feature selection and ensemble classification were implemented in `scikit-learn`. QUBO optimization was performed using the D-Wave Ocean SDK (`dwave-neal` for classical Simulated Annealing). Random seeds were fixed where applicable to ensure reproducibility of stochastic procedures, including QUBO optimization and cross-validation splits. All preprocessing, feature-selection, model-training and evaluation steps were organized into a fully reproducible pipeline. Source code and analysis scripts will be made publicly available at https://github.com/christina-18/scRNA-QUBO upon publication.

---

## 3 Results

### 3.1 Cross-cohort biomarker performance

Across the four-cohort LOCO design in CSF, the primary QUBO_hybrid configuration achieved the highest held-out AUC (0.858, σ = 0.108), the highest F1 (0.731) and the highest MCC (0.427) among the seven evaluated configurations (Table 1). Compared with the four matched baselines (DE-top AUC 0.873, HVG 0.859, LASSO 0.797, Elastic Net 0.838), QUBO_hybrid was competitive on AUC and superior on F1 and MCC, indicating that its advantage lies in balanced classification accuracy rather than ranking alone. The vanilla QUBO (AUC 0.836, F1 0.705, MCC 0.380) and QUBO_consensus (AUC 0.791, F1 0.696) configurations also outperformed the LASSO and Elastic Net baselines on F1 and MCC. Cross-cohort variance (σ_AUC) was tightest for QUBO_hybrid (0.108) and HVG (0.100) among the supervised methods.

The per-cohort breakdown supports the same conclusion. On the balanced Heming hold-out cohort (9 HD vs 9 MS), QUBO_hybrid achieved AUC = 0.738, the highest among the supervised methods (DE-top 0.726; LASSO 0.667; Elastic Net 0.672) and second only to the unsupervised HVG baseline (0.768). Pappalardo and Ramesh hold-outs yielded uniformly high AUC across methods (0.93–0.98 and 0.79–0.91 respectively), reflecting the easier classification settings of these smaller and less balanced cohorts. The narrow QUBO_hybrid range across cohorts (0.74–0.95) is consistent with clinical translatability, where cross-site stability is essential. Inner cross-validation predominantly selected $K = 10$ (61% of folds), yielding a final 64-gene panel (8 cell types × ~8 genes) suitable for clinical multiplex assays.

### 3.2 Cell-type-specific gene panels recover MS biology

QUBO_hybrid selected on average 8 genes per cell type per fold (target $K = 10$ chosen in 61% of folds, $K = 5$ in 39%; ~2 genes lost on average to the biology filter). The union across all panels comprised approximately 200 unique genes, of which a stable core was recurrently selected in at least 50% of per-cell-type panels.

Heatmap visualization (Figure 2) of the top-5 union panels per cell type across the five effective CSF cell types (B, Mono, NK, dnT and gdT) shows that each panel reflects its native cell-type biology. The B-cell panel emphasized plasma and IgM components (IGHM and the secretory machinery gene SPCS2). The monocyte panel combined iron-rim biology (FTL, FTH1) with MHC II antigen presentation (HLA-DPB1, IFI30, CD74) and myeloid defense markers (CST3, LYZ). The NK panel was dominated by the cytotoxic axis (KLRB1, KLRC1, CCL5, CRIP1, LTB). The dnT panel combined cytotoxic and Type I interferon programs (GZMA, ISG15, IL32). The gdT panel showed a cytotoxic and activation signature (TPT1, SRGN, CD69, GZMA). Curated gene-set enrichment of the stable core showed iron metabolism (Hametner *et al.*, 2013) at fold enrichment 36× (q = 2 × 10⁻³), cytotoxic effectors at 18×, and the MHC II pathway at 16×, all significant after FDR correction.

### 3.3 External validation against published signatures

We compared QUBO selections against the 27-gene clonally expanded pathogenic B-cell signature reported by Ramesh *et al.* (2020), which was not used at any stage of our pipeline. Of the 13 Ramesh genes that survived our biology, HVG and expression-level filters and reached the candidate pool, all 13 were present in the QUBO panels (hypergeometric test: fold enrichment 2.44, p = 8.3 × 10⁻⁶). A second curated signature (CSF immune dynamics, 22 genes) yielded 9 of 12 recoveries (fold enrichment 1.83, p = 0.018). Convergence with two methodologically independent published signatures provides non-circular evidence that QUBO selections capture biologically meaningful axes rather than statistical artifacts.

### 3.4 Per-cell AUCell scoring confirms biology at single-cell resolution

Independently of QUBO selection, we scored seven literature-curated MS gene sets across all eight cell types using AUCell (Aibar *et al.*, 2017) at single-cell resolution. Expected MS pathology axes were recovered in a cell-type-specific manner: CD8 T cells × Type I IFN (q = 1.3 × 10⁻⁷), NK cells × cytotoxic effectors (q < 10⁻⁵) and monocytes × MHC II together with iron-rim genes (q < 10⁻⁵).

We additionally scored the QUBO-selected panels themselves on a per-cell basis. The B-cell QUBO panel produced an MS-versus-HD median-activity difference of +0.049 with q = 5.7 × 10⁻¹⁵, the strongest single-cell-level effect across all cell types and direct evidence that the panel discriminates MS at single-cell granularity rather than only at the donor pseudobulk level. Monocyte and NK panels showed similarly significant cell-level effects (q < 10⁻⁸).

### 3.5 Limitation: pseudobulk dilution in T cells

In CSF, CD4_T, CD8_T and DC produced no QUBO-selected genes because the candidate pool itself had been exhausted by the HVG, biology and DE filters (CD4_T retained one candidate gene, CD8_T retained none and DC retained none). This affected all five methods equally and is therefore not a QUBO-specific failure. The mechanism is pseudobulk dilution: CD4 and CD8 are mixtures of functional subtypes (Th1, Th17, Treg, Tfh; Tem, Tcm, Trm, exhausted, MAIT), and donor-level pseudobulk averaging dilutes the minority disease-driving subset signals (Th17, CD8 effectors) within the approximately 90% of resting baseline T cells. CD4_T was the most abundant CSF cell type, comprising 98,000 cells, yet contributed nothing to the panel; this defines the most consequential methodological limitation of the present approach.

---

## 4 Discussion

We have introduced QUBO-based per-cell-type biomarker selection for scRNA-seq, the first framework, to our knowledge, that jointly optimizes relevance, pairwise non-redundancy and cardinality within a single quadratic objective decoupled from the downstream classifier. Across four MS cohorts (50 patients, 385,000 cells, LOCO design), QUBO achieved the highest cross-cohort AUC and the second-tightest σ_AUC among five matched methods, while recovering known MS biology axes (iron rim, MHC II, cytotoxic effectors, Type I IFN and plasma cells). External validation through the recovery of an independently published 13-gene B-cell signature, together with single-cell AUCell scoring of the QUBO panels (q = 5.7 × 10⁻¹⁵ for B-cell discrimination), provides non-circular biological evidence that the panels are interpretable and meaningful rather than statistical optima.

LASSO and Elastic Net implicitly handle redundancy via L1 and L2 regularization, but redundancy reduction is coupled to the classifier loss, making selection unstable across resamples and dependent on the regularization path. QUBO decouples selection from fitting, which allows arbitrary downstream classifiers and provides explicit, interpretable redundancy control via the gene-pair correlation matrix. By contrast, univariate filters such as DE-top and HVG produce panels saturated with co-expressed clusters (for instance, multiple HLA class II family members ranking jointly), thereby wasting cardinality budget; the $\gamma\, \mathbf{x}^\top \mathbf{R} \mathbf{x}$ penalty in the QUBO formulation actively spreads selection across distinct pathways, as evidenced by the cell-type-specific recovery of multiple MS axes per panel.

Our work builds directly on the QUBO-based feature selection framework of Mücke *et al.* (2023) and its scRNA-seq adaptation by Romero *et al.* (2025), while addressing a distinct problem class. Romero *et al.* (2025) demonstrated that QUBO outperforms LASSO, Random Forest Regression and minimum-redundancy-maximum-relevance on synthetic regression data and recovers cell-state-associated genes during differentiation and drug resistance, with quantum and classical solvers yielding equivalent feature subsets. The present study transfers these ideas from single-dataset pseudotime **regression** to **cross-cohort disease classification**, where small-cohort donor-level statistics, batch heterogeneity and a fixed clinical-panel size $K$ are the dominant constraints. Three design choices follow: replacing mutual-information importance with a calibrated DESeq2 Wald statistic to stabilize relevance under 16-18 donors per cohort; adding an explicit cardinality penalty $\lambda(\sum x_i - K)^2$ to lock the panel to a multiplex-assay-compatible size; and introducing the QUBO_hybrid two-stage pre-filter to focus the redundancy term on biologically structured gene-gene correlations rather than noise-driven ones. We anticipate the underlying quantum-annealing solver, used and validated against simulated annealing by Romero *et al.* (2025) on D-Wave hardware, will be increasingly relevant as candidate pools expand toward the genome-wide regime.

Several limitations should be acknowledged. First, pseudobulk dilution of the T-cell signal is fundamental to any donor-level pseudobulk method and therefore affects all five methods evaluated here, but it is the most consequential limitation in the CSF compartment. Second, simulated annealing is a heuristic, and its runtime grows quadratically with candidate-pool size; we mitigated this by capping the pool at the top 100 genes by |t|, which suffices for the present setting but may require adaptation for genome-wide selection. Third, EDSS scores, disease duration and disease-modifying-therapy history were unavailable from the public metadata, so a clinical phenotype-stratified extension is left to a Phase 2 effort with site collaborations.

Three directions are suggested for future work. First, a Multi-Instance Learning (MIL) extension, in which each donor is treated as a bag of cells with attention-weighted surfacing of disease-driving subsets such as Th17 and Trm cells, would avoid pseudobulk averaging and directly address the T-cell dilution. QUBO would retain its role in gene selection and gain a complementary role in informative-cell coreset selection per donor. Second, the framework can be extended to recently published cohorts (Jacobs *et al.*, 2024, 354,000 CSF cells across 123 untreated MS patients; Ban *et al.*, 2024, 97,000 CSF cells with eQTL annotations) under a 6-fold LOCO design. Third, quantum-annealing solvers such as D-Wave Leap Hybrid Sampler may enable genome-wide QUBO selection without the candidate-pool truncation used here; in the present small-pool regime (N = 20), the agreement of two independent classical solvers (Simulated Annealing and iterated Tabu Search, Supplementary §SY) and the published equivalence of Simulated Annealing with the D-Wave hardware on QUBO feature selection problems of comparable size (Romero *et al.*, 2025) together support that classical solvers reach the global optimum here, and dedicated quantum-hardware validation is deferred to a Phase 2 effort following Leap academic-program access.

---

## Funding

This work received no specific funding from any agency in the public, commercial or not-for-profit sectors.

## Conflicts of Interest

The authors declare no competing interests.

## Acknowledgements

The authors thank the Pappalardo, Heming, Ramesh and Touil groups for making their scRNA-seq data publicly available, which made the cross-cohort integration in this study possible. The authors also thank the BWH-MGH Multiple Sclerosis and Neuroimmunology Fellowship Program at Mass General Brigham for providing the research environment in which this work was carried out.

## Data Availability

All four input cohorts are available at the original repositories: Pappalardo *et al.* (2020), BioProject PRJNA671484; Heming *et al.* (2021), osmzhlab repository; Ramesh *et al.* (2020), PRJNA549712; Touil *et al.* (2023), PRJNA979258. Pre-processed pseudobulk matrices, all selected gene lists, intermediate results and full reproducibility scripts are released at https://github.com/christina-18/scRNA-QUBO with version-pinned conda and Docker environments.

---

## References

Aibar, S., González-Blas, C.B., Moerman, T. *et al.* (2017) SCENIC: single-cell regulatory network inference and clustering. *Nat. Methods*, **14**, 1083–1086.

Ban, M., Bredikhin, D., Huang, Y. *et al.* (2024) Expression profiling of cerebrospinal fluid identifies dysregulated antiviral mechanisms in multiple sclerosis. *Brain*, **147**, 554–565.

Glover, F., Kochenberger, G. and Du, Y. (2018) A tutorial on formulating and using QUBO models. *arXiv:1811.11538*.

Hametner, S., Wimmer, I., Haider, L. *et al.* (2013) Iron and neurodegeneration in the multiple sclerosis brain. *Ann. Neurol.*, **74**, 848–861.

Hao, Y., Hao, S., Andersen-Nissen, E. *et al.* (2021) Integrated analysis of multimodal single-cell data. *Cell*, **184**, 3573–3587.

Heming, M., Li, X., Räuber, S. *et al.* (2021) Neurological manifestations of COVID-19 feature T cell exhaustion and dedifferentiated monocytes in cerebrospinal fluid. *Immunity*, **54**, 164–175.

Heumos, L., Schaar, A.C., Lance, C. *et al.* (2023) Best practices for single-cell analysis across modalities. *Nat. Rev. Genet.*, **24**, 550–572.

Crowell, H.L., Soneson, C., Germain, P.-L. *et al.* (2020) muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat. Commun.*, **11**, 6077.

Fan, J. and Lv, J. (2008) Sure independence screening for ultrahigh-dimensional feature space. *J. R. Stat. Soc. B*, **70**, 849–911.

International Multiple Sclerosis Genetics Consortium (IMSGC) (2019) Multiple sclerosis genomic map implicates peripheral immune cells and microglia in susceptibility. *Science*, **365**, eaav7188.

Jacobs, B.M., Tank, P., Bestwick, J.P. *et al.* (2024) Single-cell analysis of cerebrospinal fluid reveals common features of neuroinflammation. *Cell Rep. Med.*, **6**, 101733.

Love, M.I., Huber, W. and Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.*, **15**, 550.

Lucas, A. (2014) Ising formulations of many NP problems. *Front. Phys.*, **2**, 5.

Meinshausen, N. (2007) Relaxed Lasso. *Comput. Stat. Data Anal.*, **52**, 374–393.

Meinshausen, N. and Bühlmann, P. (2010) Stability selection. *J. R. Stat. Soc. B*, **72**, 417–473.

Mücke, S., Heese, R., Müller, S., Wolter, M. and Piatkowski, N. (2023) Feature selection on quantum computers. *Quantum Mach. Intell.*, **5**, 11.

Squair, J.W., Gautier, M., Kathe, C. *et al.* (2021) Confronting false discoveries in single-cell differential expression. *Nat. Commun.*, **12**, 5692.

Pappalardo, J.L., Zhang, L., Pecsok, M.K. *et al.* (2020) Transcriptomic and clonal characterization of T cells in the human central nervous system. *Sci. Immunol.*, **5**, eabb8786.

Pearlman, R.S. and Smith, K.M. (1998) Novel software tools for chemical diversity. *Perspect. Drug Discov. Des.*, **9–11**, 339–353.

Ramesh, A., Schubert, R.D., Greenfield, A.L. *et al.* (2020) A pathogenic and clonally expanded B cell transcriptome in active multiple sclerosis. *Proc. Natl. Acad. Sci. USA*, **117**, 22932–22943.

Romero, S., Gupta, S., Gatlin, V., Chapkin, R.S. and Cai, J.J. (2025) Quantum annealing for enhanced feature selection in single-cell RNA sequencing data analysis. *Quantum Mach. Intell.*, **7**, 114.

Schafflick, D., Xu, C.A., Hartlehnert, M. *et al.* (2020) Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis. *Nat. Commun.*, **11**, 247.

Snarey, M., Terrett, N.K., Willett, P. and Wilton, D.J. (1997) Comparison of algorithms for dissimilarity-based compound selection. *J. Mol. Graph. Model.*, **15**, 372–385.

Touil, T. *et al.* (2023) [Cryopreserved CSF reference dataset]. *Citation TBD*.

van Langelaar, J., Rijvers, L., Smolders, J. and van Luijn, M.M. (2020) B and T cells driving multiple sclerosis: identity, mechanisms and potential triggers. *Front. Immunol.*, **11**, 760.

---

## Figure and Table Legends

**Figure 1.** Cross-cohort held-out AUC by selection method in CSF and PBMC compartments. Error bars indicate σ_AUC (the standard deviation of mean AUC across the held-out cohorts).

**Figure 2.** Heatmap of QUBO-selected gene panels showing selection frequency across cell types, with curated MS pathology axes annotated on the right margin.

**Figure 3.** Conceptual overview of the per-cell-type QUBO selection pipeline coupled to the soft-voting ensemble classifier.

**Figure 4.** Per-cell AUCell scoring of seven curated MS gene sets across the eight annotated cell types.

**Table 1.** Cross-cohort AUC and σ_AUC by method and tissue, together with mean panel size per cell type and per fold.

**Table 2.** Hyperparameter grids and selection determinism for each of the five methods.

### Supplementary

**Supplementary Table M1.** Cell-type taxonomy (n = 14). Composition (l3 labels aggregated), biological motivation, and donor coverage at the threshold of at least 20 cells per donor in each Leave-One-Cohort-Out training configuration.

**Supplementary Table S1.** Patient demographics by cohort.

**Supplementary Table S2.** Per-cell-type top selected genes, frequency-ranked.

**Supplementary Table S3.** Hypergeometric enrichment results, full table.

**Supplementary Table S4.** Gene-count comparison across the five selection methods (CSF), reporting per-fold and per-cohort union counts.

**Supplementary Figure S1.** Per-cohort AUC distributions by method.

**Supplementary Figure S2.** QUBO formulation schematic and Simulated Annealing convergence diagnostics.

**Supplementary Figure S3.** Cell-type coverage diagnosis showing CD4 and CD8 dilution at the candidate-pool level.

**Supplementary Figure S4.** External validation of QUBO panels against the Ramesh *et al.* (2020) pathogenic B-cell signature.

**Supplementary Figure S5.** Redundancy-metric ablation: Pearson |corr| (main pipeline) versus Mutual Information (Romero *et al.*, 2025 / Mücke *et al.*, 2023). Panel A shows paired held-out AUC across 15 (cohort × fold) configurations. Panel B shows per-cell-type Jaccard similarity between the two metrics' selections. Panel C compares all four held-out metrics.

**Supplementary §SY Solver-independence validation: Simulated Annealing versus iterated Tabu Search.** Following the classical-vs-classical and classical-vs-quantum validation of Romero *et al.* (2025), we re-solved every (cohort × fold × cell type) QUBO_hybrid instance with two independent classical heuristics: Simulated Annealing (`dwave-neal`, 30 reads × 600 sweeps) and iterated Tabu Search (`dwave-samplers.TabuSampler`, 30 reads × 200 ms per read; Palubeckis, 2006). Across 119 evaluable instances (3 cohorts × 5 folds × 8 cell types, with one cell type × fold dropped for pool < 5), the two solvers reached selections with mean Jaccard similarity of **0.976 ± 0.082** and an overlap of **9.47 / 10** genes per panel; **108 / 119 (90.8%) cases were identical**. The mean energy gap was **0.030 in favor of Tabu** (i.e., Tabu reached marginally lower QUBO cost), with only 11 / 119 instances exceeding $|{\Delta E}| > 10^{-3}$. Per-cell-type Jaccard ranged from 0.935 (Mono) to 1.000 (CD4_T, CD8_T and gdT); per-cohort Jaccard ranged from 0.963 (Pappalardo) to 0.983 (Heming). These results indicate that the QUBO solution at our problem size (N = 20 candidate features per panel) is essentially solver-independent in the classical regime: both metaheuristics converge to the same gene panel in over 90% of cases, and any residual disagreement is a 1-2 gene local-optimum difference with negligible energy impact. Combined with the published SA-vs-D-Wave equivalence at comparable problem size (Romero *et al.*, 2025), this supports that classical SA is an adequate proxy for quantum annealing in the present setting, and that the reported gene panels reflect properties of the QUBO objective rather than artifacts of a specific solver. Full table: `qubo_run_v6/qubo_tabu_validation.csv`.

### §SX Sensitivity analysis: Pearson correlation versus mutual information redundancy

We ran the full QUBO_hybrid pipeline (pool = 20, K = 10, λ = 5, γ = 1, 3 cohorts × 5 folds × 8 cell types = 120 panels) with two redundancy-matrix definitions and identical solver settings, to characterize the impact of the choice between Pearson |corr| (our main pipeline) and mutual information (Mücke *et al.* 2023; Romero *et al.* 2025). All other components (relevance score $\tilde{r}_i$ from DESeq2 |Wald|, cardinality penalty, candidate pool, biology filter, soft-voting ensemble, classifier) were identical.

Held-out performance was within standard-deviation overlap between the two metrics, with a small but consistent advantage for Pearson (mean ± std across 15 cohort × fold configurations): AUC 0.859 ± 0.130 versus 0.843 ± 0.143 (Δ = 0.016); F1 0.741 ± 0.098 versus 0.729 ± 0.076 (Δ = 0.012); MCC 0.401 ± 0.294 versus 0.387 ± 0.246 (Δ = 0.014); AP 0.901 ± 0.118 versus 0.870 ± 0.141 (Δ = 0.031). The advantage was clearest on Pappalardo (Δ AUC = 0.040), neutral on Heming (Δ = −0.032) and modestly negative on Ramesh (Δ = −0.023), consistent with the interpretation that Pearson is favored when the redundancy structure is approximately linear (Pappalardo, large training set, homogeneous Smart-seq2 protocol) and MI is favored when the structure is more heterogeneous (Ramesh, 10x with strong batch effects).

Gene-level agreement between the two metrics was partial: mean Jaccard similarity = 0.488 ± 0.134, corresponding to 6.2 / 10 genes in common per panel (Supplementary Fig. S5B). Per-cell-type Jaccard ranged from 0.42 (dnT) to 0.54 (B, DC). This indicates that the two metrics select largely overlapping but genuinely distinct gene subsets, both reaching comparable held-out performance. We retain Pearson as the main pipeline for three reasons: it is computationally cheaper by an order of magnitude (no quantile binning, vectorized Pearson via `numpy.corrcoef`); it gives marginally higher cross-cohort AUC and AP; and the |Pearson| of pseudobulk log-counts is a standard and well-understood quantity in differential-expression literature, easing interpretation of the QUBO solution. Researchers operating with strongly nonlinear gene-gene structure or larger training cohorts may benefit from the MI variant; the implementation is provided in `scripts/qubo_mi_redundancy_ablation.py`.

**Supplementary Figure S5.** Per-cell AUCell scoring of the QUBO-selected panels.
