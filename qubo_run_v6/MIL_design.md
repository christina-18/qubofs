# Multi-Instance Learning (MIL) Extension — Design Document

**Project**: QUBO-Optimized Cell-Type-Specific Gene Panels for MS Classification
**Author**: Mizuho Asada, Ph.D
**Status**: Future work design (Phase 2)
**Last updated**: 2026-05-12

---

## 1. Motivation

The current implementation uses **per-donor pseudobulk** as the unit of statistical
analysis. This design has well-established advantages:

- avoids pseudoreplication (Squair et al., *Nat Comm* 2021)
- treats donor as the natural clinical unit
- yields stable, interpretable models on small N

The v6 DESeq2 tight pipeline (manuscript main result) achieves cross-cohort
**held-out AUC 0.858, F1 0.731 and MCC 0.427** for QUBO_hybrid — the
strongest values on F1 and MCC among the seven evaluated methods. However,
pseudobulk averages out **cell-level heterogeneity**, and **two specific
limitations of the current results directly motivate the MIL extension**:

### 1.1 CD4_T and CD8_T subset dilution (primary motivation)

The v6 DESeq2 tight selections in CSF show that **CD4_T and CD8_T panels are
dominated by myeloid-axis genes** (TYROBP, AIF1) rather than canonical T-helper
or cytotoxic markers (IL17/FOXP3 for CD4; GZMA/B for CD8). CD4_T (5 panels of
TYROBP/AIF1 in 15 folds) and CD8_T (CCL3L1 dominant, no canonical CD8 effector
gene) demonstrate that the **disease-relevant minority subsets — Th17, Treg,
exhausted CD8, CD8 Trm — are submerged within the 90% baseline of resting
T cells** at the donor-pseudobulk level.

This is the dilution mechanism: CD4_T was the most abundant CSF cell type
(98,000 cells across the integrated dataset) yet contributed essentially
nothing biologically interpretable to the panel. MIL with attention pooling
directly addresses this by letting the model up-weight the rare disease-driving
cells per donor.

### 1.2 Iron-laden microglia and oligoclonal B cells

A complementary motivation comes from non-CSF biology: iron-laden microglia
in chronic active MS lesions (Hametner et al., *Ann Neurol* 2013) and
clonally expanded plasma-cell B cells (Cepok 2005; Obermeier 2008) are
known to be rare subsets within their respective compartments. The CSF
v6 selections recovered the BCR repertoire diversity (IGLC1, IGHA2, IGKV4-1,
IGHV2-70) for B cells but missed canonical plasma markers (MS4A1, BTK).
A subset-aware extension is needed to surface these populations.

**Multi-Instance Learning (MIL)** offers a principled framework to retain
cell-level information while keeping the donor as the statistical unit. This
document outlines the design for an MIL extension as Phase 2 of the project.

---

## 2. Problem Formulation

### 2.1 Bag-of-instances structure

Each donor *d* is treated as a **bag** containing *N_d* individual cells:

$$\mathcal{B}_d = \{\mathbf{c}_1, \mathbf{c}_2, \ldots, \mathbf{c}_{N_d}\}, \quad y_d \in \{0, 1\}$$

where:

- $\mathbf{c}_i \in \mathbb{R}^G$ is the gene-expression vector of cell *i* (dimension *G* = number of selected genes)
- $y_d$ is the donor's diagnosis label (1 = MS, 0 = HD)
- *N_d* varies per donor (typically a few hundred to several thousand cells)

The instance labels (whether each individual cell is "MS-like") are unknown.
Only the bag-level label is observed.

### 2.2 Goal

Learn a function $f: \mathcal{B}_d \to [0, 1]$ that predicts the donor-level
MS probability while remaining permutation-invariant with respect to cell ordering.

### 2.3 Statistical unit

The donor remains the unit of cross-validation and external validation.
Cells are **not** independent samples; they are nested within donors.

---

## 3. Architecture: Attention-based MIL

We adopt the gated-attention MIL framework of Ilse et al. (*ICML* 2018).

### 3.1 Network components

**(i) Cell encoder** $\phi: \mathbb{R}^G \to \mathbb{R}^H$

Maps each cell's expression vector to an *H*-dimensional embedding.
Implementation: a small MLP (2 hidden layers, ReLU) shared across cells.

**(ii) Attention module** $\alpha: \mathbb{R}^H \to \mathbb{R}$

Computes a scalar attention weight per cell:

$$a_i = \frac{\exp\{\mathbf{w}^\top \tanh(V \phi(\mathbf{c}_i)) \odot \mathrm{sigm}(U \phi(\mathbf{c}_i))\}}{\sum_{j=1}^{N_d} \exp\{\cdot\}}$$

where *V*, *U*, *w* are learnable parameters and $\odot$ is element-wise product.

**(iii) Bag aggregator**

$$\mathbf{z}_d = \sum_{i=1}^{N_d} a_i \, \phi(\mathbf{c}_i)$$

The attention weights $\{a_i\}$ are interpretable as the relative contribution
of each cell to the donor's prediction.

**(iv) Bag classifier** $\psi: \mathbb{R}^H \to [0, 1]$

A linear layer with sigmoid output predicting MS probability.

### 3.2 Loss

Binary cross-entropy at the bag level, with optional L2 regularization on
encoder weights to mitigate overfitting given small donor count.

---

## 4. The Role of QUBO in the MIL Pipeline

QUBO retains a central role in the MIL extension. Two complementary uses are
proposed.

### 4.1 Use 1: Gene selection (current QUBO role, preserved)

The cell encoder $\phi$ takes a *G*-dimensional input. Reducing *G* via QUBO
gene selection has the same advantages as in the pseudobulk pipeline:

- noise reduction
- interpretable feature space
- compatibility with cell-type-specific biology

In v6 DESeq2 tight, QUBO_hybrid (K = 10) produced a **union of ~500 unique
genes** across 3 cohorts × 5 folds × 8 cell types, of which only 4 genes
(CCL3L1, CDKN2B, TRH, TYROBP) appeared in ≥ 50% of per-cell-type panels;
the strict DESeq2 + apeglm shrinkage favors cohort-specific hits. For the
MIL extension we propose using the **per-cohort-fold panel (~53 genes per
fold)** as the encoder input rather than the much-smaller stable set, to
retain enough diversity for the attention module to discover MS-relevant
substructure within cells.

### 4.2 Use 2: Cell selection (novel QUBO formulation)

A new QUBO formulation can select an informative **coreset of cells** per donor
rather than using all cells in the bag.

For donor *d* with *N_d* cells, define binary variables $x_i \in \{0, 1\}$:

$$H(\mathbf{x}) = -\sum_{i=1}^{N_d} s_i x_i + \gamma \sum_{i,j} R_{ij} x_i x_j + \lambda \left( \sum_i x_i - K_{cells} \right)^2$$

where:

- $s_i$ is a per-cell relevance score (e.g., distance from the healthy cell
  distribution in latent space, or pre-computed cell-level disease likelihood)
- $R_{ij}$ is a cell-cell similarity measure (e.g., cosine similarity in the
  encoder embedding space)
- $K_{cells}$ is the target number of cells to retain per donor (e.g., 100)

This produces a **sparse, diverse, disease-relevant subset of cells** per donor,
which is then passed to the MIL aggregator. Conceptually, it is a learned
coreset selection driven by binary optimization.

### 4.3 Use 3 (forward-looking): Joint gene–cell selection

If both gene selection and cell selection are formulated as a joint binary
optimization with cross-terms, the resulting QUBO has approximately
$|G_{candidates}| + N_d$ binary variables. For donor cohorts of typical size,
this exceeds the practical capacity of classical Simulated Annealing but
remains tractable on **D-Wave Leap Hybrid CQM solvers** (up to ~10⁶ variables).

This use case provides a direct demonstration of the quantum-annealer
compatibility claim made for the current pseudobulk pipeline.

---

## 5. Cross-Validation Strategy

The MIL model must be validated at the **donor level**, identical to the
pseudobulk pipeline:

- 3-cohort × Leave-One-Cohort-Out cross-validation as the primary metric
- Internal 5-fold on training cohorts for hyperparameter tuning
- AUC, AP, F1, MCC, and per-cohort σ as the comparison metrics

**Critical**: cell-level metrics (e.g., AUC computed over individual cells) are
**not interpretable** here, because cells from the same donor inherit the same
label and are not independent. Such metrics will appear inflated and must be
avoided in reporting.

---

## 6. Expected Outcomes

### 6.1 Quantitative

The MIL extension is benchmarked against the **v6 DESeq2 tight pseudobulk
baseline** (QUBO_hybrid CSF: AUC 0.858, F1 0.731, MCC 0.427, σ_AUC 0.108
across Pappalardo / Heming / Ramesh hold-outs; 3 cohorts × 5 folds × 8 cell
types = 120 panels).

The MIL extension is expected to:

- match or modestly improve donor-level AUC over the pseudobulk baseline
  (estimated +0.01 to +0.03 AUC; target ≥ 0.86)
- maintain or improve cross-cohort σ_AUC (the current 0.108 is the second-best
  among supervised methods and is the strongest current advantage over
  ElasticNet at σ 0.145)
- substantially improve F1 and MCC by surfacing minority disease-driving cells
  in T-cell compartments (current QUBO_hybrid F1 0.731, MCC 0.427 are already
  best-in-class, so improvement here is an aggressive target)
- substantially improve **interpretability** via per-cell attention weights —
  in particular, identification of the **CD4 / CD8 subset that drives the
  prediction**, which the pseudobulk pipeline cannot resolve

The principal value of MIL is **interpretability and biological discovery**
at the T-cell-subset level, not necessarily raw classification accuracy.
Concretely, the CD4_T and CD8_T panels in the v6 DESeq2 tight selections
collapse onto myeloid-axis genes (TYROBP, AIF1) rather than canonical T-cell
effector or regulatory markers; if MIL recovers Th17, Treg, exhausted CD8,
or Trm signatures as high-attention subsets, this would be a substantive
biological advance even at unchanged donor-level AUC.

### 6.2 Qualitative — biological discovery

By examining the attention weights learned by the MIL model, we expect to
identify:

- **Which CD4_T subset** drives the MS prediction — Th17 (RORC / IL17),
  exhausted Treg (FOXP3 + HAVCR2), or Tfh (CXCR5) — currently invisible in
  the pseudobulk panel which collapses onto TYROBP / AIF1
- **Which CD8_T subset** drives the prediction — CD8 effector (GZMA / GZMB /
  PRF1), exhausted (PDCD1 / LAG3 / HAVCR2), or Trm (CD69 / ITGAE) — currently
  represented only by CCL3L1
- **Whether the TYROBP / AIF1 hit in CD4_T / CD8_T fractions** reflects a
  genuine T-myeloid interaction signature or rare myeloid contamination in
  the T-cell cluster (resolvable by attention)
- **Whether the iron-laden microglia signature** (Hametner 2013) emerges in
  high-attention monocyte cells without being explicitly seeded
- **Whether plasma B-cell markers** (MS4A1, BTK) missed by the pseudobulk panel
  emerge as attention-weighted cells in the BCR-repertoire-rich B-cell bag
- **Heterogeneity within MS donors** — the proportion of "MS-like" cells per
  donor as a candidate disease-activity marker

---

## 7. Implementation Outline

### 7.1 Data preparation

Extract per-donor cell-level expression matrices for the QUBO-selected gene
set. For the v6 DESeq2 tight pipeline, the **per-cohort-fold panel
(~53 genes/fold across 8 cell types)** is the recommended input; the union
across all folds reaches approximately 500 unique genes, while the
≥ 50%-stable core is only 4 genes (CCL3L1, CDKN2B, TRH, TYROBP) and is too
small to feed an MIL encoder.

```r
# extract_cells_per_donor.R (sketch)
library(Seurat)
seurat <- readRDS("so.GEX.share.Asada.rds")
seurat <- DietSeurat(seurat, assays = "RNA",
                     dimreducs = NULL, scale.data = FALSE)
# QUBO_hybrid union from v6 DESeq2 tight (~500 unique genes across 8 cell types)
qubo_genes <- readLines("qubo_run_v6/v6deseq2tight_qubo_hybrid_union.txt")

for (donor in unique(seurat$donor_id)) {
  cells <- subset(seurat, donor_id == donor & compartment == "CSF")
  expr <- GetAssayData(cells, slot = "data")[qubo_genes, ]
  saveRDS(list(expr = expr,
               cell_types = cells$predicted.celltype.l2,
               diagnosis = unique(cells$Dx),
               cohort = unique(cells$prj)),
          file = paste0("cells_per_donor/", donor, ".rds"))
}
```

### 7.2 MIL training (PyTorch sketch)

```python
import torch, torch.nn as nn

class GatedAttentionMIL(nn.Module):
    def __init__(self, G, H=64):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(G, H), nn.ReLU(), nn.Linear(H, H), nn.ReLU())
        self.attention_V = nn.Linear(H, H)
        self.attention_U = nn.Linear(H, H)
        self.attention_w = nn.Linear(H, 1)
        self.classifier = nn.Linear(H, 1)

    def forward(self, cells):              # cells: (N_d, G)
        h = self.encoder(cells)            # (N_d, H)
        a_V = torch.tanh(self.attention_V(h))
        a_U = torch.sigmoid(self.attention_U(h))
        a = self.attention_w(a_V * a_U)    # (N_d, 1)
        a = torch.softmax(a, dim=0)
        z = (a * h).sum(dim=0)             # (H,)
        return torch.sigmoid(self.classifier(z)), a
```

Training proceeds bag by bag with binary cross-entropy at the bag level.
Donor-stratified cross-validation as in the pseudobulk pipeline.

### 7.3 Compute requirements

- Cell extraction (R): ~30–60 minutes (single pass over the 12 GB Seurat object)
- MIL training (PyTorch on CPU): ~2–4 hours per CV fold
- Total: roughly half a day on a standard workstation

GPU acceleration (single consumer-grade card) reduces training to roughly
20 minutes per fold.

---

## 8. Comparison With Current Pipeline

| Aspect | Pseudobulk + Ensemble (current) | Attention MIL (proposed) |
|---|---|---|
| Statistical unit | Donor | Donor |
| Cell-level info | Averaged out | Preserved via attention |
| Per-cell interpretability | None | Attention weights per cell |
| Number of trainable parameters | ~17 per cell-type classifier (×8) | ~5,000–20,000 (encoder + attention) |
| Suitable N | ≥30 donors | ≥50 donors (current borderline) |
| Implementation maturity | Established (current paper) | Phase 2 (this document) |
| QUBO role | Gene selection | Gene selection + (optional) cell selection |

---

## 9. Reporting Plan for the MIL Extension

The MIL extension is intended to appear in the manuscript as:

- A dedicated **"Cell-Level Extension"** section in Discussion or as a
  follow-up figure
- A **comparison table** of pseudobulk vs MIL AUC across the three cohorts
- A **case-study figure** showing attention weights for one MS donor and one
  HD donor, with high-attention cells annotated by cell type and selected
  pathway markers
- An **explicit note** that cell-level AUC is not used as a metric due to
  pseudoreplication concerns

If the MIL results meaningfully outperform pseudobulk on σ_AUC or recover
known pathology subsets (e.g., iron-laden microglia), they may motivate a
separate manuscript focused on the cell-level discovery.

---

## 10. Risks and Mitigations

| Risk | Mitigation |
|---|---|
| Overfitting on N=50 with deep MIL | Aggressive L2 regularization, dropout, encoder pretraining on cell-type prediction |
| Attention weights uninterpretable in noise | Sparsity-promoting attention (entropic regularization) |
| Cell-level extraction time on 12 GB .rds | Restrict to QUBO-selected genes only (~49 cols) and CSF compartment to keep file sizes manageable |
| Comparison favors pseudobulk on small data | Frame contribution as **interpretability + biological discovery**, not just AUC |

---

## 11. Timeline (Indicative)

| Phase | Activity | Estimated duration |
|---|---|---|
| Phase 2.1 | R cell extraction, save per-donor RDS | 1 week |
| Phase 2.2 | Implement Attention MIL in PyTorch | 1 week |
| Phase 2.3 | Train + cross-cohort validation | 1 week |
| Phase 2.4 | Cell selection QUBO formulation and integration | 2 weeks |
| Phase 2.5 | Attention-based biological discovery analysis | 2 weeks |
| **Total** | | **~7 weeks** |

This timeline assumes part-time effort and is anchored to Phase 1 (the current
pseudobulk pipeline) as the baseline.

---

## References

- Ilse M, Tomczak J, Welling M. Attention-based deep multiple instance
  learning. *ICML* 2018; 80: 2127–2136.
- Squair JW, *et al.* Confronting false discoveries in single-cell differential
  expression. *Nat Commun* 2021; 12: 5692.
- Hametner S, *et al.* Iron and neurodegeneration in the multiple sclerosis
  brain. *Ann Neurol* 2013; 74: 848–861.
- Cepok S, *et al.* Short-lived plasma blasts are the main B cell effector
  subset during the course of multiple sclerosis. *Brain* 2005; 128: 1667–1676.
- Obermeier B, *et al.* Matching of oligoclonal immunoglobulin transcriptomes
  and proteomes of cerebrospinal fluid in multiple sclerosis. *Nat Med* 2008;
  14: 688–693.
- Mücke S, Heese R, Müller S, Wolter M, Piatkowski N. Feature selection on
  quantum computers. *Quantum Mach Intell* 2023; 5: 11.
- Romero S, Gupta S, Gatlin V, Chapkin RS, Cai JJ. Quantum annealing for
  enhanced feature selection in single-cell RNA sequencing data analysis.
  *Quantum Mach Intell* 2025; 7: 114.
- Lopez R, *et al.* Deep generative modeling for single-cell transcriptomics.
  *Nat Methods* 2018; 15: 1053–1058. (scVI)
- Xu C, *et al.* Probabilistic harmonization and annotation of single-cell
  transcriptomics data with deep generative models. *Mol Syst Biol* 2021; 17:
  e9620. (scANVI)
