# Data sources

The four input CSF single-cell RNA-seq cohorts are publicly available from their
original repositories and are **not redistributed in this repository**. Download
them from the accessions below and place the integrated, SoupX-corrected Seurat
object where `QUBOFS_SEURAT_RDS` points (see `docs/reproduction.md`).

| Cohort | Reference | Accession | Internal label (in scripts) | Donors (control / MS) |
|---|---|---|---|---|
| Pappalardo | Pappalardo *et al.* (2020), *Sci. Immunol.* | BioProject PRJNA671484 | `PRJNA671484_MS_Tcell` | 6 / 5 |
| Heming | Heming *et al.* (2021), *Immunity* | GEO GSE163005 | `osmzhlab_MS_ence_cov` | 9 / 9 |
| Ramesh | Ramesh *et al.* (2020), *PNAS* | BioProject PRJNA549712 | `PRJNA549712_MS_PBMC_UCSF` | 3 / 14 |
| Touil | Touil *et al.* (2023), *Cell Rep. Methods* | BioProject PRJNA979258 | `PRJNA979258_cryoCSF` | 4 / 0 |

The "Internal label" column lists the raw dataset identifiers used in the analysis
scripts and is not meaningful to readers; in particular, the Ramesh label contains
the substring `PBMC` for historical script reasons, but the present study uses only
the **CSF compartment** of every cohort (the peripheral-blood compartment is used
solely in the optional generalisation benchmark, Supplementary Table S7).

The analysis is restricted to the CSF compartment of each cohort. The integrated
dataset comprises 50 donors (28 MS, 22 control), 71 CSF samples and 221,066 cells
across eight broad immune cell types. Heming uses idiopathic intracranial
hypertension (IIH) donors as non-inflammatory neurological controls; Touil
contributes control donors only and is retained in the training set throughout.

Upstream preprocessing before integration: ambient-RNA decontamination with SoupX
(Young & Behjati 2020) and doublet exclusion with scDblFinder (Germain *et al.*
2021). Cell-type labels use Azimuth `predicted.celltype.l2` collapsed into eight
broad immune subsets (B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT).

To reproduce the downstream analysis **without** the raw data, use the shipped
`data_release/` tables and `examples/toy_data/` (synthetic) as described in the
top-level `README.md`.
