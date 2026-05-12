# English Narration — Slides 3 & 6 (academic register, short sentences)

**Presentation:** QUBO-Optimized Cell-Type-Specific Gene Panels for Cross-Cohort
MS Classification
**Presenter:** Mizuho Asada, Ph.D
**Style:** academic lab meeting tone. Short sentences, periods rather than dashes.
**Pace:** ~70-80 sec per slide, deliberate pauses on key numbers.
**Date:** May 8, 2026

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
> The most striking finding is the iron-metabolism set. Hametner et al. 2013 reported iron-laden microglia forming a rim around chronic active MS lesions. These appear as paramagnetic rim lesions on SWI and QSM MRI. FTH1 and FTL are the central iron-storage molecules. The analysis shows fold enrichment of 36, with q-value of 2 times 10 to the minus 3. Cytotoxic effectors show 18-fold enrichment. The MHC II pathway shows 16-fold enrichment. All three results are statistically significant.
>
> The picture sharpens further at the broader 448-gene selection. All 10 Reactome MHC class II genes are selected. Of the 11 IMSGC top GWAS hits, 9 are recovered. The recovered genes include BACH2, CXCR4, IL7R, and the HLA cluster.
>
> In summary, QUBO does not simply optimize statistically. It converges on the biological axes that decades of MS research have identified as central. MHC class II antigen presentation, T-cell activation, Type I IFN signaling, cytotoxicity, and iron metabolism are all represented with high reproducibility. This is the property required of a clinically translatable biomarker panel.

---

## Pacing & delivery tips

- **Brief pauses on key numbers**: 385,116, 17, 448, 49, 36-fold, q of 2 times 10 to the minus 3. Land each number, then breathe.
- **Wet-lab terms** to articulate carefully: oligoclonal band, iron rim, paramagnetic rim lesion, laminin receptor, specialized ribosome.
- **Slide transitions, kept brief and formal**:
  - 3 → 4: "Now to the QUBO formulation itself."
  - 5 → 6: "Next, the biological interpretation of the selected genes."
  - 6 → 7: "Pulling these results together."

## Anticipated Q&A

> Q: Why collapse the cell types to eight rather than using the finer subdivision?

A: "Finer subdivision is feasible in principle. With 50 donors, however, splitting CD4 into Naive, TCM, and TEM separately leaves a number of donors with insufficient cells per subset for stable pseudobulk. We required at least 20 cells per donor per subset. Eight populations is the resolution at which all cohorts have reproducible annotation and stable pseudobulk. With more cohorts we would scale to a finer resolution."

> Q: Is the iron-metabolism enrichment biologically meaningful?

A: "We believe so. Hametner et al. 2013 in Annals of Neurology reported iron-laden microglia forming a rim around chronic active MS lesions. These appear as paramagnetic rim lesions on susceptibility-weighted MRI. FTH1 and FTL are the central iron-storage molecules in those microglia. The fact that QUBO recurrently selects them in CSF myeloid populations is consistent with the disease pathology. The signal is not merely statistical. It reflects what is actually occurring in the lesion."

> Q: Which is the final panel, the stable 49 or the full 448?

A: "Both views are useful. The stable 49 is the conservative core. These genes are reproducibly selected across cohorts and folds. The 448 is the broader exploration. It represents the genes QUBO considered important across configurations. For clinical translation we would propose the stable subset, in the range of 20 to 49 genes. The enrichment results are significant on both views, with magnitudes of 36-fold on curated sets at the stable level and approximately 2-fold across the broader selection."
