# English Narration — 7 main slides

**Date:** May 7, 2026

---

## Main Slide 1 — Title + Self-introduction  *(≈ 50 sec)*

Thank you for your time. My name is Mizuho Asada.

I am an Assistant Professor at Meiji Pharmaceutical University, and a Lecturer at Institute of Science Tokyo.

I am currently on sabbatical here at MGH.

My research at MGH focuses on applying optimization methods to bioinformatics.

The title of my talk is

"QUBO-optimized gene panels for cross-cohort classification of MS from single-cell RNA-seq."

To give a brief overview of my background, I work on three areas.

- First, cheminformatics.
- Second, drug-effect modeling.
- Third, clinical AI.

This work reflects my main interest here at MGH:

bringing optimization methods into bioinformatics for MS biomarker design.

Now, I will move to the introduction.

---

## Main Slide 2 — Introduction (Background + Aim + Dataset) *(≈ 50 sec)*

MS diagnosis still relies on MRI and CSF findings, and established molecular or cellular biomarkers remain limited.

Single-cell RNA-seq studies have advanced rapidly, but reproducibility across cohorts remains inconsistent

due to batch effects and donor variability.

Therefore, there is a need for a framework that provides reproducible biomarkers across cohorts.

The aim of this study is simple. We use QUBO optimization to achieve both biomarker performance and biological insight within a single framework.

The dataset includes four cohorts: 50 patients, 99 samples, and about 385,000 cells.

Heming provides CSF data only. Pappalardo and Ramesh include both CSF and PBMC.

Touil includes HD only and is used for training.

---

## Main Slide 3 — Method  *(≈ 75 sec)*

First, the intuition of QUBO.

Each gene is represented as a binary variable —selected or not.

We define a cost function, H of x, and we select genes that minimize this cost.

This cost function has three components.

- First, relevance. Genes with strong differential expression receive higher scores.
- Second, non-redundancy. Highly correlated genes are penalized, so redundant genes are not selected.
- Third, cardinality. We control the number of selected genes, K, which is chosen by cross-validation.

Next, the pipeline.

For each cell type, we perform differential expression, select K genes, and train a logistic regression model.

The predictions from each cell type are averaged to obtain the final patient-level prediction.

Finally, evaluation. We use leave-one-cohort-out cross-validation across three cohorts.

We compare five methods: DE-top, HVG, LASSO, Elastic Net, and QUBO. All settings are identical — only the selection method differs.

---

## Main Slide 4 — Data Preparation Flow  *(≈ 60 sec)*

Now, the data preparation pipeline. It consists of seven steps.

- First, we start with an integrated single-cell dataset of about 385,000 cells.
- Second, we annotate cell types using Azimuth (アジマス) and group them into eight major immune cell types.
We chose this resolution to ensure at least 20 cells per donor.
（Fewer cells would lead to unstable averages and noisy estimates.）
- Third, we construct pseudobulk data per donor for both CSF and PBMC.
This makes the donor the unit of analysis.
- Fourth, we filter genes. After HVG selection and biological filtering, about 8,000 genes remain.
- Fifth, we perform differential expression analysis and select the top 100 genes as candidates.
- Sixth, QUBO selects the final gene set. On average, 17 genes per cell type are selected.
- Finally, we combine predictions across cell types to obtain the patient-level result.

Evaluation uses leave-one-cohort-out cross-validation.

---

## Main Slide 5 — Results 1: Cross-cohort biomarker performance (Table 1)  *(≈ 50 sec)*

Now, the first results: biomarker performance.

Table 1 shows AUC and variability across cohorts.

QUBO achieves the best performance in both CSF and PBMC.

In CSF, the AUC is 0.788, with a low variability of 0.044, indicating stable performance across cohorts.

For reference, LASSO shows higher variability, while Elastic Net shows similar stability.

Clinically, CSF reflects disease biology, while PBMC is non-invasive and suitable for screening.

---

## Main Slide 6 — Results 2: Selection stability of QUBO-selected genes (Figure 1)  *(≈ 60 sec)*

Next, Figure 1 shows the selected genes.

This heatmap shows which genes are consistently selected across cell types.

Each cell type reflects its own biology.

- For example, monocytes show MHC class II and iron-related genes.
- NK and DN T cells show cytotoxic markers.
- Gamma-delta T cells also show cytotoxic signals.
- B cells show the IgM signature.

**Click**

Importantly, CD4, CD8, and dendritic cells show no selection.

This is likely due to signal attenuation from pseudobulk aggregation.

This motivates our future MIL extension.

Finally, external validation.

All 13 genes from the Ramesh et al. signature are recovered in our QUBO panel.

This was done without using that paper as input,

providing independent validation of the biology.

---

## Main Slide 7 — Conclusion + Future Work  *(≈ 55 sec)*

The main message of this work is that QUBO delivers both biology and biomarker performance within a single framework.

- First, biomarker performance.
We achieved a CSF AUC of 0.788 with a cross-cohort sigma of 0.044.
- Second, biological validity.
Monocytes capture MHC class II and iron-related signals.
NK, DN T, and gamma-delta T cells capture the cytotoxic axis.
B cells capture the oligoclonal band signature — recovering central MS pathology.
- Third, cross-site reproducibility.
AUC remains stable, ranging from 0.74 to 0.82 across all three cohorts, supporting robustness for clinical translation.
- Fourth, methodological novelty.
This is the first framework that jointly optimizes relevance, non-redundancy, and panel size within a single objective function.
Importantly, redundancy control is explicitly built into the selection process — something standard methods cannot directly achieve.

In summary, this work provides a unified scRNA-seq framework that delivers both biomarker performance and biological interpretability.

It has potential for both clinical application and mechanistic insight.

As future work, we plan a cell-level extension using Multi-Instance Learning.

A key limitation of this study is that CD4 and CD8 T cells were not selected, likely due to signal attenuation from pseudobulk aggregation combined with gene filtering.

However, CD4 Th17 and CD8 cytotoxic programs are central to MS pathobiology, making their recovery a top priority.

MIL models each donor as a bag of cells and uses attention to identify disease-relevant subsets, such as Th17-like CD4 cells or Trm-like CD8 cells.

By preserving cell-level resolution, this approach addresses the dilution problem at its source.

QUBO will retain its role in gene selection and may be extended to select informative cell subsets per donor.

Thank you for your attention. I welcome your questions.

---

## Supplementary slides — when to use them

S1 = robustness, S2–4 = biology, S5 = PBMC vs CSF, S6 = QUBO math, S7–8 = limitations & future

### S1 — Detailed performance / robustness

**When to show**:
- "Are you concerned about overfitting?"
- "How stable is the model across cohorts?"
- "What about MCC or F1?"

**Main message**: Performance is stable across independent cohorts.

**Simple explanation**:
QUBO showed stable performance across cohorts, with relatively low variability.
This suggests reasonable generalization rather than dataset-specific fitting.

### S2 — Top genes & enrichment

**When to show**:
- "Which genes were selected?"
- "What biology did you capture?"

**Main message**: The selected genes recover known MS biology.

**Simple explanation**:
We identified genes related to MHC II, cytotoxicity, interferon signaling, and iron metabolism,
which are all known MS-related pathways.

### S3 — Functional annotation

**When to show**:
- "Are these genes biologically meaningful?"
- "Which genes are MS-related?"

**Main message**: Many selected genes are already linked to MS.

**Simple explanation**:
For example, HLA-DPB1 is related to MHC II, FTL to iron metabolism, and GZMA to cytotoxic activity.
These are consistent with known MS pathology.

### S4 — Literature comparison

**When to show**:
- "How does this compare with previous studies?"
- "What is missing?"

**Main message**: QUBO recovers major MS pathways, but some T-cell signals are still missing.

**Simple explanation**:
We successfully recovered several established MS pathways, especially innate immune and MHC II signals.
However, CD4 Th17 and some CD8 signals were not fully captured.

### S5 — PBMC vs CSF

**When to show**:
- "Why is CSF better?"
- "Can PBMC still be useful?"

**Main message**: CSF shows stronger disease signals, but PBMC is clinically easier to obtain.

**Simple explanation**:
CSF is closer to the lesion site, so disease-associated signals are stronger.
PBMC signals are weaker, but blood sampling is non-invasive.

### S6 — QUBO formulation

**When to show**:
- "Can you explain the equation?"
- "What does the redundancy penalty mean?"
- "What solver did you use?"

**Main message**: The objective balances relevance, non-redundancy, and panel size.

**Simple explanation**:
The first term maximizes relevance, the second term penalizes redundancy,
and the third term controls the number of selected genes.

**If asked about quantum computing**:
This study used classical simulated annealing, but the formulation is compatible with quantum annealers.

### S7 — CD4/CD8 limitation

**When to show**:
- "Why are CD4/CD8 absent?"
- "Is this biologically reasonable?"

**Main message**: The signal is diluted by pseudobulk averaging.

**Simple explanation**:
CD4 and CD8 contain many heterogeneous subtypes.
Disease-relevant subsets may represent only a small fraction, so averaging dilutes the signal.

### S8 — Single-cell validation

**When to show**:
- "Does the panel work at single-cell level?"
- "Can you distinguish MS cells directly?"

**Main message**: The QUBO panel also shows signal at single-cell resolution.

**Simple explanation**:
AUCell analysis showed that the selected gene panels remain active at the single-cell level,
especially in B cells and monocytes.

---

## Anticipated Q&A (Simplified Version)

### 1. QUBO / Method

> **Q: What is the main advantage of QUBO?**
> A: QUBO jointly optimizes three things: relevance, non-redundancy, and panel size. This is difficult with standard methods.

> **Q: Why not just use LASSO or Elastic Net?**
> A: LASSO and Elastic Net select important genes, but they do not explicitly avoid redundant genes.
> QUBO directly penalizes highly correlated features.

> **Q: What does the redundancy term do?**
> A: It prevents selecting many similar genes. This improves interpretability and robustness.

> **Q: How is K determined?**
> A: K is selected by inner cross-validation from predefined candidates. Different cell types may require different panel sizes.

> **Q: Why use logistic regression instead of deep learning?**
> A: We wanted to isolate the effect of feature selection, not model complexity.

### 2. Data Preparation / Pseudobulk

> **Q: Why did you use pseudobulk?**
> A: Pseudobulk enables donor-level analysis and reduces noise from single cells.
> It also helps avoid data leakage.

> **Q: Why require at least 20 cells per donor?**
> A: To ensure stable pseudobulk estimation.

> **Q: How did you handle batch effects?**
> A: We used integrated preprocessing and evaluated robustness by cross-cohort validation.

> **Q: What is cross-cohort validation?**
> A: We trained on some cohorts and tested on completely independent cohorts.

### 3. Performance / Robustness

> **Q: Is an AUC of 0.788 sufficient?**
> A: It is promising for biomarker discovery, although further validation is still needed.

> **Q: Why is PBMC performance lower than CSF?**
> A: CSF is closer to the disease site, so disease-related signals are stronger.

> **Q: What does sigma AUC mean?**
> A: It measures variability across cohorts. Lower sigma AUC means better robustness.

> **Q: Are you concerned about overfitting?**
> A: Cross-cohort validation and low sigma AUC suggest reasonable generalization.

> **Q: Is AUC computed using only the 49 stable genes?**
> A: No. The prediction uses fold-specific QUBO panels. The stable genes are mainly for biological interpretation.

### 4. Biological Interpretation

> **Q: How do you interpret the selected genes biologically?**
> A: The selected genes reflect known MS biology, including cytotoxicity, MHC II pathways, and interferon responses.

> **Q: Is the overlap with Ramesh et al. meaningful?**
> A: Yes. Their data were not used in training, so this serves as independent validation.

> **Q: Why do NK and innate-like T cells appear prominently?**
> A: These cell types are enriched in CSF and are important in neuroinflammation.

### 5. CD4/CD8 & MIL

> **Q: Why are CD4 and CD8 missing?**
> A: Their signals are diluted by pseudobulk aggregation because they contain many heterogeneous subtypes.

> **Q: Does this mean CD4/CD8 are not important?**
> A: No, actually the opposite. They are central to MS biology, but harder to capture with averaging.

> **Q: Why do you need MIL?**
> A: MIL can detect disease-relevant cells directly, without averaging all cells together.

> **Q: What is MIL?**
> A: MIL stands for Multi-Instance Learning. It learns donor-level labels while identifying important cells within each donor.

> **Q: How will you combine QUBO with MIL?**
> A: MIL identifies important cells, and QUBO selects informative genes and cell subsets.

### 6. Overall Contribution

> **Q: What is the main contribution of this work?**
> A: This work provides a unified framework for both biomarker prediction and biological interpretation.
