# Data sources

## What you need, and where it comes from

**If you only want to reproduce the manuscript figures and tables**, you do
**not** need the raw data at all: the released summary tables in `data_release/`
are sufficient, and `python scripts/make_canonical_figures.py` regenerates the
canonical figures from them. See the top-level `README.md` and
`docs/reproduction.md` for that (Seurat-free) path.

**Those `data_release/` tables were produced by the authors** from the raw
single-cell data of the four published cohorts below, by running the full
pipeline (`scripts/01_pipeline` → `scripts/04_aggregation`). The raw data
themselves are publicly available but are **not redistributed in this
repository** — download them from the accessions in the table below.

- **Which files to download from each BioProject/GEO page, and how the
  downloaded raw data were turned into each `data_release/*.csv` file, are
  documented in [`data_release/README.md`](../data_release/README.md).** That
  file is the entry point for anyone wanting to regenerate the released tables
  from scratch.
- Once you have built the integrated, SoupX-corrected Seurat object, place it
  where `QUBOFS_SEURAT_RDS_RAW` points and follow `docs/reproduction.md` from
  step 0b.

## The four input cohorts

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
`data_release/` tables (see [`data_release/README.md`](../data_release/README.md)
for what each file is and how it was produced) and `examples/toy_data/`
(synthetic), as described in the top-level `README.md` and `docs/reproduction.md`.
