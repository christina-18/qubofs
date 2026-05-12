# English Narration — All 7 Slides (academic register)

**Presentation:** QUBO-Optimized Cell-Type-Specific Gene Panels for Cross-Cohort
Classification of Multiple Sclerosis from Single-Cell RNA Sequencing
**Presenter:** Mizuho Asada, Ph.D
**Style:** academic lab meeting tone. Short sentences, periods rather than dashes.
**Total runtime:** ~7 minutes (~1 minute per slide)
**Date:** May 8, 2026

---

## Slide 1 — Title  *(≈ 25 sec)*

> Today's presentation addresses cross-cohort classification of multiple sclerosis using QUBO-optimized cell-type-specific gene panels.
>
> The work designs reproducible biomarker panels from single-cell RNA sequencing data, robust to technical differences across study sites. The central idea is to formulate gene selection as a quadratic unconstrained binary optimization, or QUBO, problem. This formulation allows the joint optimization of three properties that conventional methods address only partially.
>
> I will begin with a brief introduction and the motivation behind this work.

---

## Slide 2 — Background & Data  *(≈ 75 sec)*

> First, a brief introduction.
>
> My name is Mizuho Asada. I am an Assistant Professor at the Laboratory of Medical Molecular Analytics, Meiji Pharmaceutical University. I also serve as a Lecturer at the Institute of Science Tokyo Graduate School, in the Department of Anesthesiology and Pain Medicine. I am currently on sabbatical at MGH.
>
> My research focuses on PK/PD modeling, machine-learning-based predictive and image-analysis models, and cheminformatics. At MGH, I am particularly interested in gene selection within bioinformatics. Today's project is one concrete example of that direction.
>
> Now to the background.
>
> Multiple sclerosis affects approximately 2.8 million people worldwide. Single-cell RNA-seq has accelerated diagnostic biomarker research. However, classification accuracy collapses when training and test cohorts differ in technical conditions. Specific factors include 10x chemistry version, cryopreservation status, and sampling protocol.
>
> Existing methods have limitations. DE-top and HVG provide only univariate ranking and cannot control redundancy. LASSO controls the number of selected genes K only indirectly through the regularization parameter lambda.
>
> In this work, we propose a QUBO-based approach that jointly optimizes relevance, non-redundancy, and cardinality.
>
> The dataset integrates four published cohorts. It contains 50 donors, 99 samples, 385,116 cells, and 32,170 genes. Heming contributes 18 donors, Pappalardo 11, Ramesh 17, and Touil 4.

---

## Slide 3 — Data Preparation Flow  *(≈ 70 sec)*

> Next, the data preparation pipeline.
>
> The input is an integrated single-cell RNA-seq object. It contains 385,116 cells, 50 donors, and 4 cohorts. Pappalardo, Heming, and Ramesh are rotated as the held-out cohort under leave-one-cohort-out cross-validation. The Touil cohort lacks MS donors and is therefore fixed to the training set.
>
> Cell-type annotation comes from Azimuth's predicted.celltype.l2 reference. The original annotation defines more than 30 subtypes. These include CD4 Naive, CD4 TCM, CD4 TEM, B naive, B memory, and others. Given the donor count of 50, we collapsed these into 8 broader populations. Three criteria guided the choice. First, the population must be a major MS-relevant lymphoid or myeloid group. Second, each donor must have at least 20 cells per population on average for stable pseudobulk. Third, the annotation must be reproducible across all cohorts. The selected populations are B, Mono, CD4_T, CD8_T, NK, DC, dnT, and gdT.
>
> We split cells into three compartments. These are CSF, PBMC, and the integrated ALL. Today we focus on CSF as the primary compartment.
>
> The pseudobulk construction follows. For each donor and each cell type, expression values are averaged across all cells of that type. This produces a donor-by-gene matrix per cell type. The design choice is that the unit of statistical analysis is the donor, not the individual cell. Cell-level analyses are known to inflate false positives due to pseudoreplication, as reported by Squair et al. 2021. Pseudobulk avoids this issue.
>
> Gene filtering proceeds in three stages. HVG selection reduces from 32,170 genes to 3,000 per cell type. The biology filter follows current best practice from Heumos et al. 2023 and Luecken & Theis 2019. It removes mitochondrial, ribosomal, heat-shock, MALAT1, NEAT1, and classical housekeeping genes. Two specific exceptions are retained. The RPLP family contributes to specialized ribosome programs, as established by Genuth & Barna 2018. RPSA encodes the 67-kDa laminin receptor, central to leukocyte transmigration across the blood-brain barrier per Nelson et al. 2008. Both have direct relevance to MS pathobiology.
>
> edgeR then performs differential expression analysis. The model includes covariates for diagnosis, log10 cell count, age, sex, and batch. DESeq2 and limma-voom are run in parallel as sensitivity checks. Top-100 lists overlap by more than 90 percent across the three methods.
>
> The final step takes the top 100 genes by absolute t-statistic per cell type. This serves as the QUBO candidate pool. Each instance contains 100 candidates. Across all training configurations, the universe is approximately 1,090 unique genes for CSF.

---

## Slide 4 — Method  *(≈ 80 sec)*

> Now to the QUBO formulation itself.
>
> Gene selection is defined by the following objective function. H of x equals minus s prime x, plus gamma times x prime R x, plus lambda times the squared deviation of the sum of x from K. Here, x is the binary vector representing gene selection.
>
> The first term maximizes relevance. The score s sub i is defined as the squared t-statistic from the edgeR differential expression analysis. Genes with strong differential expression between MS and HD receive higher scores.
>
> The second term minimizes redundancy. The pairwise correlation matrix R, weighted by the parameter gamma, penalizes the selection of co-expressed gene pairs. This avoids redundant selection of genes carrying similar information. This is the component that conventional univariate methods cannot handle and is the key novelty of the QUBO approach.
>
> The third term enforces cardinality. It is a soft penalty pulling the selection count toward K. In this work, K is auto-selected from 10, 20, or 30 by inner 5-fold cross-validation.
>
> The solver is classical Simulated Annealing. Each instance takes approximately 3 seconds. The same objective function can also run on D-Wave quantum annealers without modification. This makes the design future-proof for emerging quantum hardware.
>
> We compare five gene-selection methods under identical conditions. These are DE-top, HVG, LASSO, Elastic Net, and QUBO. All methods share the same candidate pool, the same K grid, the same L2 logistic regression classifier, and the same 8-cell-type soft-voting ensemble. The only difference between methods is the selection logic itself.
>
> Only QUBO jointly optimizes relevance, non-redundancy, and cardinality. Evaluation uses 3-cohort leave-one-cohort-out cross-validation.

---

## Slide 5 — Results  *(≈ 80 sec)*

> Now to the results. CSF compartment, held-out metrics averaged across 3 cohorts.
>
> The top table shows the metrics. QUBO achieves AUC 0.788, F1 0.635, and MCC 0.258. These rank first among all five methods. AP is 0.846, narrowly behind Elastic Net at 0.870 by 0.024. Sigma AUC is 0.044, essentially tied with Elastic Net at 0.041. LASSO scores 0.779, DE-top 0.742, and HVG 0.712.
>
> A brief note on the metrics. AUC is the area under the ROC curve, representing the ranking capability between MS and HD. AP is the area under the Precision-Recall curve, sensitive to minority-class detection. F1 is the harmonic mean of Precision and Recall. MCC is the Matthews Correlation Coefficient, known as the most imbalance-robust metric. Sigma AUC is the standard deviation of AUC across the 3 held-out cohorts and reflects cross-site reproducibility.
>
> The key observation is that QUBO ranks first on the three primary classification metrics. These are AUC, F1, and MCC. Cross-cohort stability, measured by sigma AUC, is on par with Elastic Net.
>
> The bottom table shows per-cohort held AUC. QUBO scores 0.807 on Pappalardo, 0.738 on Heming, and 0.819 on Ramesh. All three cohorts converge to a tight range of 0.74 to 0.82. LASSO spans a wider range from 0.72 to 0.85.
>
> In other words, QUBO maintains stable performance on unseen cohorts. This is an essential property for clinical translation.

---

## Slide 6 — Selected Genes & Biology  *(≈ 80 sec)*

> Slide 6 addresses the biological validity of the QUBO-selected genes.
>
> Performance metrics tell only part of the story. The biological interpretation of the selected genes is essential, particularly for wet-lab researchers.
>
> The numbers first. Across 3 cohorts and 5 folds, 69 panels were generated. QUBO selects on average 17 genes per cell type. K is auto-selected from 10, 20, or 30 by inner cross-validation. The union across all panels is 448 unique genes. Among these, 49 genes are recurrently selected in at least half of the panels per cell type. We define this set as the stable core. It serves as the foreground for enrichment analysis.
>
> The heatmap on the left shows selection frequency. Rows are genes. Columns are cell types. Color intensity indicates the percentage of panels in which each gene was selected.
>
> The Mono column is particularly informative. CST3 reaches 100 percent. SAT1 and FTL both reach 80 percent. HLA-DPB1, IFI30, FTH1, and CD74 follow. This combination represents MHC class II antigen presentation and iron metabolism. These are central axes of MS monocyte pathobiology. The dnT column shows GZMA, ISG15, and CCL5 at high frequencies. These reflect cytotoxic and Type I interferon signatures. NK and gdT cluster around cytotoxic genes. These include KLRB1, KLRC1, GNLY, and GZMA.
>
> The dot plot on the right presents pathway enrichment. Green dots indicate MS-curated gene sets. Navy dots indicate GO biological process terms.
>
> The most striking finding is the iron-metabolism set. Hametner et al. 2013 in Annals of Neurology reported iron-laden microglia forming a rim around chronic active MS lesions. These appear as paramagnetic rim lesions on SWI and QSM MRI. FTH1 and FTL are the central iron-storage molecules. The analysis shows fold enrichment of 36, with q-value of 2 times 10 to the minus 3. Cytotoxic effectors show 18-fold enrichment. The MHC II pathway shows 16-fold enrichment. All three results are statistically significant.
>
> The picture sharpens further at the broader 448-gene selection. All 10 Reactome MHC class II genes are selected. Of the 11 IMSGC top GWAS hits, 9 are recovered. The recovered genes include BACH2, CXCR4, IL7R, and the HLA cluster.
>
> In summary, QUBO does not simply optimize statistically. It converges on the biological axes that decades of MS research have identified as central. MHC class II antigen presentation, T-cell activation, Type I IFN signaling, cytotoxicity, and iron metabolism are all represented with high reproducibility. This is the property required of a clinically translatable biomarker panel.

---

## Slide 7 — Conclusion  *(≈ 55 sec)*

> Finally, the take-home messages.
>
> First, performance. QUBO achieves AUC 0.788, F1 0.635, and MCC 0.258 on CSF. It ranks first on the three primary classification metrics.
>
> Second, cross-site reproducibility. The cohort-to-cohort standard deviation of 0.044 is essentially tied with Elastic Net at 0.041. It substantially outperforms LASSO at 0.068 and DE-top at 0.065.
>
> Third, biological validity. The QUBO-selected genes converge on the central axes of MS pathobiology. These include the MHC class II antigen presentation axis with HLA-DRB1, HLA-DPB1, CD74, and IFI30. They include the Type I interferon signature with ISG15. They include cytotoxic effectors such as GZMA and CCL5. They include iron metabolism with FTH1 and FTL. GO enrichment analysis shows MHC class II antigen presentation enriched at p less than 10 to the minus 8.
>
> Fourth, quantum-annealer compatibility. Formulated as binary optimization, the same objective function is directly deployable on D-Wave quantum annealers. Classical Simulated Annealing already runs in approximately 3 seconds per instance. The design transitions seamlessly as quantum hardware matures.
>
> In conclusion, QUBO simultaneously satisfies five criteria. Highest AUC. Highest F1 and MCC. Cross-cohort stability on par with Elastic Net. Biological validity. Quantum compatibility. As a cross-site-reproducible MS biomarker panel, it is positioned as a strong candidate for clinical translation.
>
> A brief note on future work. The current study uses a pseudobulk approach. As Phase 2, we plan a cell-level extension based on Multi-Instance Learning. Each donor is treated as a bag of cells, and attention mechanisms aggregate cells while preserving per-cell interpretability. This identifies which cell subpopulations drive the prediction. QUBO retains its role in gene selection, with the additional novel role of selecting informative cell coresets per donor. A joint gene-cell QUBO is planned for the D-Wave Leap Hybrid Solver, which will also serve as direct quantum-annealer validation.
>
> Thank you for your attention. I welcome your questions.

---

## Pacing & delivery tips

- **Brief pauses on key numbers**: 385,116, 17, 448, 49, 0.788, 36-fold, q of 2 times 10 to the minus 3, p less than 10 to the minus 8. Land each number, then breathe.
- **Wet-lab terms** to articulate carefully: oligoclonal band, iron rim, paramagnetic rim lesion, laminin receptor, specialized ribosome.
- **Slide transitions, kept brief and formal**:
  - 1 → 2: "I will begin with a brief introduction and motivation."
  - 2 → 3: "Next, the data preparation pipeline."
  - 3 → 4: "Now to the QUBO formulation itself."
  - 4 → 5: "Now to the results."
  - 5 → 6: "Next, the biological interpretation of the selected genes."
  - 6 → 7: "Pulling these results together."

## Anticipated Q&A

> Q: Why collapse the cell types to eight rather than the finer subdivision?

A: "Finer subdivision is feasible in principle. With 50 donors, however, splitting CD4 into Naive, TCM, and TEM separately leaves a number of donors with insufficient cells per subset for stable pseudobulk. We required at least 20 cells per donor per subset. Eight populations is the resolution at which all cohorts have reproducible annotation and stable pseudobulk. With more cohorts we would scale to a finer resolution."

> Q: Is the iron-metabolism enrichment biologically meaningful?

A: "We believe so. Hametner et al. 2013 in Annals of Neurology reported iron-laden microglia forming a rim around chronic active MS lesions. These appear as paramagnetic rim lesions on susceptibility-weighted MRI. FTH1 and FTL are the central iron-storage molecules in those microglia. The fact that QUBO recurrently selects them in CSF myeloid populations is consistent with the disease pathology. The signal is not merely statistical. It reflects what is actually occurring in the lesion."

> Q: Which is the final panel, the stable 49 or the full 448?

A: "Both views are useful. The stable 49 is the conservative core. These genes are reproducibly selected across cohorts and folds. The 448 is the broader exploration. It represents the genes QUBO considered important across configurations. For clinical translation we would propose the stable subset, in the range of 20 to 49 genes. The enrichment results are significant on both views, with magnitudes of 36-fold on curated sets at the stable level and approximately 2-fold across the broader selection."

> Q: Has the formulation been validated on actual D-Wave hardware?

A: "Currently the work runs entirely on classical Simulated Annealing. Each instance takes about 3 seconds, so quantum annealing is not required at the present problem size. For larger-scale runs, for example optimizing over the full set of approximately 2,500 HVGs at once, the D-Wave Leap Hybrid Solver becomes appropriate. Hardware validation is part of the planned future work."
