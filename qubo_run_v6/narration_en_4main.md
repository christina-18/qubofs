# English Narration — 4 main + 4 supplementary structure

**Presentation:** QUBO-Optimized Cell-Type-Specific Gene Panels for Cross-Cohort
Classification of Multiple Sclerosis from Single-Cell RNA Sequencing
**Presenter:** Mizuho Asada, Ph.D
**Style:** academic lab meeting tone. Short sentences, periods rather than dashes.
**Audience:** MS specialist group (basic background omitted)
**Total runtime:** ~5 minutes (4 main slides) + supplementary on demand

---

## Main Slide 1 — Headline + Background + Dataset  *(≈ 80 sec)*

> Today's presentation addresses cross-cohort classification of multiple sclerosis using QUBO-optimized cell-type-specific gene panels.
>
> Let me start with the bottom line. Forty-nine stable QUBO genes achieve AUC 0.788 ± 0.044 across three held-out cohorts. The selected genes recover central axes of MS pathobiology, namely MHC class II antigen presentation, iron metabolism, Type I interferon signaling, and cytotoxicity.
>
> Three take-home messages. First, performance. QUBO ranks first on the three primary classification metrics — AUC, F1, and MCC — among five methods. Second, reproducibility. The cross-cohort sigma AUC is 0.044, essentially tied with Elastic Net at 0.041 and approximately 60 percent of LASSO at 0.068. Third, biological validity. MHC class II is enriched 8.5-fold, and iron metabolism 36-fold. The latter is consistent with the iron rim reported by Hametner et al. 2013.
>
> The background focuses on the technical challenge. Existing MS scRNA-seq classifiers collapse when training and test cohorts differ in technical conditions. Specific factors include 10x chemistry version, cryopreservation status, and sampling protocol. Existing methods such as DE-top, LASSO, and HVG cannot jointly optimize relevance, redundancy, and cardinality. This work resolves the issue through QUBO-based joint optimization of all three axes, combined with cross-cohort validation. The classical solver is Simulated Annealing. The same formulation is directly deployable on D-Wave quantum annealers.
>
> The dataset integrates four cohorts, with 50 donors, approximately 385,000 cells, and 32,000 genes. Heming contributes 18 donors, Pappalardo 11, Ramesh 17, and Touil 4. Touil lacks MS donors and is therefore fixed to the training set. External evaluation uses 3-cohort leave-one-cohort-out cross-validation.

---

## Main Slide 2 — Method  *(≈ 75 sec)*

> Now to the QUBO formulation.
>
> Gene selection is defined by an objective function with three terms.
>
> The first term maximizes relevance. The score s sub i is the squared absolute t-statistic from edgeR. Genes with strong differential expression between MS and HD receive higher scores.
>
> The second term minimizes redundancy. The pairwise correlation matrix R, weighted by the parameter gamma, penalizes the selection of co-expressed gene pairs. This avoids redundant selection of genes from the same pathway. This is the component conventional univariate methods cannot handle and is the core contribution of QUBO.
>
> The third term enforces cardinality. It is a soft penalty pulling the selection count toward K. K is auto-selected from 10, 20, or 30 by inner 5-fold cross-validation.
>
> The solver is classical Simulated Annealing. Each instance takes approximately 3 seconds. The same formulation is directly deployable on D-Wave quantum annealers without modification.
>
> We compare five gene-selection methods under identical conditions. These are DE-top, HVG, LASSO, Elastic Net, and QUBO. Importantly, all five share the same candidate pool, the same K grid, the same L2 logistic classifier, and the same 8-cell-type soft-voting ensemble. The only difference is the selection logic itself. This isolates the effect of the selection method.

---

## Main Slide 3 — Results — Performance & Biological Validity  *(≈ 90 sec)*

> Now to the results. The CSF compartment held-out metrics averaged across three cohorts, together with the biological interpretation of selected genes, on a single slide.
>
> The top table shows the metrics. QUBO achieves AUC 0.788, F1 0.635, and MCC 0.258. These rank first on the three primary classification metrics. AP is 0.846, narrowly behind Elastic Net at 0.870 by 0.024. Sigma AUC is 0.044, essentially tied with Elastic Net at 0.041. Per-cohort, QUBO scores 0.807 on Pappalardo, 0.738 on Heming, and 0.819 on Ramesh — a tight 0.74 to 0.82 range across all three cohorts. Detailed per-cohort tables are available in supplementary slide S3.
>
> The lower half shows the biological interpretation of the selected genes.
>
> The heatmap on the left shows selection frequency. Rows are genes. Columns are cell types. Color intensity indicates the percentage of panels in which each gene was selected. The Mono column is particularly informative. CST3 reaches 100 percent. SAT1 and FTL reach 80 percent. HLA-DPB1, IFI30, FTH1, and CD74 follow. This combination represents MHC class II antigen presentation and iron metabolism — central to MS monocyte pathobiology. The dnT column shows GZMA, ISG15, and CCL5 at high frequencies, reflecting cytotoxic and Type I IFN signatures.
>
> The dot plot on the right shows pathway enrichment. Green dots indicate MS-curated gene sets. Navy dots indicate GO biological process terms. The most striking is iron metabolism, with fold enrichment of 36. This corresponds directly to the iron rim of chronic active lesions reported by Hametner et al. 2013. Cytotoxic effectors show 18-fold enrichment. The MHC II pathway shows 16-fold enrichment. All three are statistically significant.
>
> If we widen our view to the full 448-gene QUBO selection, the picture sharpens. All 10 Reactome MHC class II genes are selected. Of the 11 IMSGC top GWAS hits, 9 are recovered, including BACH2, CXCR4, IL7R, and the HLA cluster.

---

## Main Slide 4 — Conclusion + Future Work  *(≈ 65 sec)*

> Finally, the take-home messages.
>
> First, performance. QUBO ranks first on the three primary classification metrics on CSF.
>
> Second, cross-site reproducibility. The cohort-to-cohort standard deviation of 0.044 is essentially tied with Elastic Net at 0.041, and substantially outperforms LASSO and DE-top.
>
> Third, biological validity. The selected genes converge on the central axes of MS pathobiology — MHC class II antigen presentation, Type I IFN, cytotoxicity, and iron metabolism. GO enrichment shows MHC class II at p less than 10 to the minus 8.
>
> Fourth, quantum-annealer compatibility. Formulated as binary optimization, the same objective function is directly deployable on D-Wave quantum annealers. Classical Simulated Annealing already runs in approximately 3 seconds per instance. The design transitions seamlessly as quantum hardware matures.
>
> In conclusion, QUBO simultaneously satisfies five criteria. Highest AUC. Highest F1 and MCC. Cross-cohort stability on par with Elastic Net. Biological validity. Quantum compatibility. As a cross-site-reproducible MS biomarker panel, it is positioned as a strong candidate for clinical translation.
>
> As future work, we plan a cell-level extension based on Multi-Instance Learning. Each donor is treated as a bag of cells, and attention mechanisms aggregate cells while preserving per-cell interpretability. This identifies which cell subpopulations drive the prediction. QUBO retains its role in gene selection, with the additional novel role of selecting informative cell coresets per donor. A joint gene-cell QUBO is planned for the D-Wave Leap Hybrid Solver, also serving as direct quantum-annealer validation.
>
> Thank you for your attention. I welcome your questions.

---

## Supplementary slides — when to use them

Switch to the relevant supplementary slide in response to specific questions.

### S1 — Presenter Background

**Anticipated questions**: "Could you tell us your background?" / "What is your research field?" / "How does this connect to cheminformatics?"

**Brief talk (~30 sec)**:
> Briefly on my background. I am an Assistant Professor at the Laboratory of Medical Molecular Analytics, Meiji Pharmaceutical University, and a Lecturer at the Institute of Science Tokyo Graduate School in the Department of Anesthesiology and Pain Medicine. I am currently on sabbatical at MGH. My research focuses on PK/PD modeling, machine-learning-based predictive and image-analysis models, and cheminformatics. This work transfers a standard cheminformatics framework — joint optimization of relevance, diversity, and cardinality — into scRNA-seq gene selection. During my sabbatical at MGH, I am particularly interested in working on gene selection within bioinformatics.

### S2 — Data Preparation Flow

**Anticipated questions**: "How was the data preprocessed?" / "How were the 8 cell types chosen?" / "What does the biology filter exclude?"

**Brief talk (~60 sec)**:
> Some additional pipeline detail. The original Azimuth annotation defines more than 30 subtypes. We collapsed these into 8 populations using three criteria. First, the population must be a major MS-relevant lymphoid or myeloid group. Second, each donor must have at least 20 cells per population on average for stable pseudobulk. Third, the annotation must be reproducible across all cohorts. The biology filter follows Heumos 2023 and Luecken & Theis 2019, removing mitochondrial, ribosomal, heat-shock, nuclear lncRNA, and housekeeping genes. We deliberately retain RPLP0/1/2 and RPSA. RPLP belongs to the specialized-ribosome literature. RPSA encodes the laminin receptor and is essential for blood-brain-barrier transmigration — directly relevant to MS. Pseudobulk makes the donor the unit of statistical analysis.

### S3 — Per-cohort Detailed Results

**Anticipated questions**: "How variable is performance across cohorts?" / "Are you concerned about overfitting?" / "Is N=50 sufficient?"

**Brief talk (~60 sec)**:
> Per-cohort detail. QUBO scores Pappalardo 0.807, Heming 0.738, Ramesh 0.819. All within a tight 0.74-0.82 range. LASSO spans the wider 0.72-0.85 range.
>
> Statistically, each per-cell-type classifier has K=17 features and 22 events, giving EPV of 1.3 — below the classical Peduzzi 1996 threshold of 10. However, L2 regularization, the 8-cell-type ensemble, and cross-cohort validation mitigate overfitting. The empirical sigma AUC of 0.044 across cohorts is direct evidence that overfitting is not severe. Modern guidelines from van Smeden 2019 and Riley 2019 accept EPV in the 2-5 range when prediction is the goal and regularization is in place.

### S4 — Top Genes per Cell Type & Curated Enrichment Detail

**Anticipated questions**: "Which specific genes were selected?" / "What about biology beyond iron metabolism?" / "How does this relate to MS GWAS?"

**Brief talk (~60 sec)**:
> Per-cell-type top genes and curated enrichment detail. In monocytes: CST3, SAT1, FTL, HLA-DPB1, LYZ, IFI30, FTH1, CD74, TPT1. In NK: KLRB1, CCL5, LTB, CRIP1, KLRC1, GNLY. In dnT: GZMA, IL32, ISG15, CCL5, TXK.
>
> Curated set enrichment was computed on the 49 stable genes. Iron metabolism shows 36-fold enrichment. Cytotoxic effectors 18-fold. MHC II 16-fold. Across the broader 448-gene selection, all 10 Reactome MHC II genes are selected, and 9 of the 11 IMSGC top MS GWAS hits are recovered, including BACH2, CXCR4, IL7R, and the entire HLA cluster. In B cells, IGHM, IGKC, and IGLC2 are selected, consistent with the oligoclonal band signature characteristic of MS.

---

## Anticipated Q&A

> Q: How was the threshold of ≥20 cells per donor chosen?

A: "It is an empirical threshold for pseudobulk stability. Below 10 cells, the expression mean becomes unstable and reproducibility across cohorts deteriorates. Twenty cells is the standard threshold at which per-donor cell-type profiles stabilize, as adopted in Squair 2021 and many scRNA-seq pseudobulk publications."

> Q: Has the formulation been validated on actual D-Wave hardware?

A: "Currently the work runs entirely on classical Simulated Annealing. Each instance takes about 3 seconds, so quantum annealing is not yet required at this scale. For larger instances, for example optimizing over the full 2,500 HVGs simultaneously, the D-Wave Leap Hybrid Solver becomes appropriate. This is part of the planned future work."

> Q: Can you predict the MS likelihood of an individual cell?

A: "Yes, this is exactly the cell-level extension we plan as Phase 2 via Multi-Instance Learning. Each donor is treated as a bag, and attention mechanisms aggregate cells while preserving per-cell interpretability. QUBO can also extend beyond gene selection to selecting informative cells per donor. The full design is documented in MIL_design.md."

## Pacing & delivery tips

- **Brief pauses on key numbers**: 385,000, 0.788, 0.044, 36-fold, p less than 10 to the minus 8.
- **Wet-lab terms** to articulate carefully: oligoclonal band, iron rim, paramagnetic rim lesion, laminin receptor.
- **Slide transitions, brief and formal**:
  - 1 → 2: "Now to the QUBO formulation."
  - 2 → 3: "Now to the results."
  - 3 → 4: "Finally, the take-home messages."
- **Switching to supplementary slides should sound natural**: "Let me show the detail in the supplementary slide."
