# Phase 0: Is the alamandine/MrgD axis a credible target in AML? A pre-registered, calibrated negative

**Question.** Is MrgD (gene *MRGPRD*), the alamandine receptor, expressed in AML blasts, haematopoietic stem and progenitor cells (HSPCs) or the marrow niche? Is there a MasR–MrgD heteromer substrate, a dependency, or a clinical signal? The study tested four pre-registered go/no-go gates, followed by specificity and robustness checks.

**Answer.** No. At transcript level, *MRGPRD* is not expressed in human AML or its niche above approximately 1 TPM, across 12 datasets with stated sensitivity.

| Gate | Prediction | Result |
|---|---|---|
| 0.1 *MRGPRD* present (bulk > 1 TPM or > 5% of cells) | FAIL | **FAIL**: 0 of 6 datasets. Bulk medians: TCGA 0.121, Beat AML 0.000 TPM (absent-gene control *ALB* 0.082 / 0.025). 6 of 289,966 cells. |
| 0.2 *MAS1* and *MRGPRD* co-detected | UNCERTAIN | **Not passed**: 0 co-detected cells; untestable (≤ 2 *MRGPRD*+ cells per compartment). |
| 0.3 Dependency in AML lines | FAIL | **FAIL**: Chronos *MRGPRD* −0.10, *MAS1* +0.05 (*MYB* −1.26). |
| 0.4 RAS score vs overall survival | WEAK PASS | **FAIL**: Beat AML HR 1.03 per SD (0.90–1.17); robust to ELN 2022. |

**Follow-up checks:**
- The same pipeline detects *MRGPRD* at 2.2–237 TPM in expressing tumours.
- *MRGPRD* reads map uniquely against its paralogues.
- AML MSCs and AML core biopsies: 0 reads where 1 TPM predicts 22–24 (p ≤ 4 × 10⁻¹⁰).
- PMA-differentiated THP-1: 0–1 reads where 1 TPM predicts 17–35.
- Secondary result: the RAS score is associated with 46 of 122 ex vivo drug responses (venetoclax ρ = 0.38), independently of *MRGPRD*. It is not interpreted as MrgD biology.

Full narrative: `manuscript/manuscript_draft_v1.md`. Gate-by-gate report: `protocol/02_PHASE0_REPORT_gate_verdicts.md`. Every deviation: `docs/DEVIATIONS_AND_ANALYST_CHOICES.md`.

---

## Folder structure

```
README.md                      this file
run_all.py                     runs steps 00–11 in order (python run_all.py [--from 05])
protocol/                      00 kickoff pack · 01 FROZEN addendum (SHA-256 466cf2ff…6844) · 02 Phase 0 report
manuscript/                    manuscript_draft_v1.md (Scientific Reports draft; [AUTHOR]/[CITATION NEEDED] placeholders)
scripts/                       config.py (all paths) · common.py · sc_common.py · eln2022.py · numbered pipeline steps
results/tables/                every table the paper cites (CSV)
results/figures/               Figures 1–4 (PNG 300 dpi + PDF)
docs/                          deviations log · figure legends · THP-1 paper appraisal · reproduction check
data/manifests/                SHA-256 manifests for open downloads and user-supplied files
env/                           package versions (Python, R) · SETUP.md · BiocSingular_stub (Windows GSVA workaround)
data/open/, data/interim/      created at run time, not distributed (see .gitignore)
```

## Pipeline

| Step | Script | Produces | Lang |
|---|---|---|---|
| 00 | `00_download_open_data.py` | downloads + `data/manifests/open_data_manifest.csv` | Py |
| 01 | `01_bulk_gate01_and_prep.py` | `gate01_bulk_detection.csv`, `gate01_bulk_per_sample_tpm.csv`, gate 0.4 inputs | Py |
| 02a–d | `02a…`, `02b…`, `02c….R`, `02d…` | single-cell cell tables → `gate01_sc_detection_summary.csv`, `gate01_sc_calibration.csv`, `gate02_codetection.csv` | Py/R |
| 03 | `03_gate03_depmap.py` | `gate03_depmap.csv`, per-line gene effects, DepMap expression (incl. THP-1) | Py |
| 04 | `04_gate04_survival.R` | `gate04_cox.csv`, `gate04_scores.csv`, `gate04_gsva_gene_log.csv` | R |
| 05 | `05_positive_control_pancancer.py` | `step1_positive_control_gdc_star.csv`, `step1_pancan_rsem_MRGPRD_ALB.csv` | Py |
| 06 | `06_mappability_paralogs.py` | `step2_mrgprd_paralog_mappability.csv` | Py |
| 07a/b | `07a_eln2022_derive.py`, `07b_eln2022_cox.R` | `step3_*` (ELN 2022 sensitivity) | Py/R |
| 08 | `08_drug_panel.py` | `step4_drug_panel_spearman.csv`, `step4_…EXPLORATORY_monocytic_adjusted.csv` | Py |
| 09 | `09_aml_stroma.py` | `step5_*` (AML MSC + core biopsies, Poisson bound) | Py |
| 10 | `10_thp1_pma.py` | `step8_pma_thp1_expression.csv` | Py |
| 11 | `11_figures.py` | `results/figures/Figure1–4` | Py |

Step 6 (proteomics) is not automated: PeptideAtlas blocks automated queries. Look up Q8TDS7 (MrgD) and P04201 (MasR) manually.

## Data

**Downloaded automatically (step 00; SHA-256 recorded):**
- Beat AML 2.0 (github.com/biodev/beataml2.0_data): counts, normalised expression, probit drug fits, clinical, sample mapping, WES mutations;
- GEO: GSE253355 (raw + authors' Seurat object), GSE107490, GSE212615, GSE348126, GSE323383, GSE318527;
- GENCODE v36 GTF.

**Supplied locally (paths in `scripts/config.py`, or set environment variables; hashes in `data/manifests/user_supplied_data_manifest.csv`):**
- `PHASE0_USER_DATA`: DepMap `CRISPRGeneEffect.csv`, `CRISPRGeneDependency.csv`, `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, downloaded 20 Nov 2026 [release name to record].
- `PHASE0_NEW_AML`: GEO RAW tars `GSE116256_RAW.tar`, `GSE241989_RAW.tar`, `GSE207104_RAW.tar` (ftp.ncbi.nlm.nih.gov/geo/series/…/suppl/).
- `PHASE0_TCGA_STAR` and `PHASE0_TCGA_FILES`: the 151 GDC TCGA-LAML "STAR - Counts" files plus their GDC file metadata JSON (`files` endpoint with `cases.submitter_id`, `cases.samples`).

**Queried live (results may change if the services update):** cBioPortal (`laml_tcga_pub` clinical and mutations; PanCancer Atlas RSEM), GDC files/data API, Ensembl REST.

## Reproduce

```
python run_all.py                 # from the repository root, in the phase0py env; set RSCRIPT to the phase0r Rscript
```
Runtime is about 20 minutes on 24 cores, excluding downloads (~4.3 GB). Peak memory is under 16 GB (step 02c loads the 1.5 GB Seurat atlas).

**Verification (21 Nov 2026).** All steps were re-run from these scripts. The regenerated tables were compared with the tables produced during the original interactive analysis (`docs/reproduction_check_vs_original_run.csv`):
- 10 tables are identical (maximum numeric difference 0).
- Three tables differ in structure, so the automated column check records `False` for them. Targeted checks are in the `targeted_check_note` column:
  - `gate02_codetection.csv` has 29 rows instead of 30, because the 10x HSPC/blast-like and GMP-like clusters are merged before tabulation. All totals are identical (289,966 cells, 6 *MRGPRD*+, 885 *MAS1*+, 0 co-detected).
  - `step3_tcga_eln2022.csv` keeps 3 of the original 9 columns; the 6 dropped helper columns are not used downstream. The ELN 2022 labels are identical for all 151 patients.
  - `step5_stroma_sensitivity_bound.csv` adds a library-size column. The shared values differ by at most 0.04%, from rounding in the original.

## Integrity notes

- The protocol was frozen locally by hash before any target value was viewed. OSF registration was waived by the study owner, and the manuscript states this.
- Bugs caught during analysis, all fixed before verdicts or re-run:
  - a DepMap column-order bug, exposed by the MYB/SPI1/ALB controls;
  - an ELN 2022 translocation-order bug, after which the Cox models were re-run;
  - a figure-only label overlap.
- The heteromer claim it tests: `docs/APPRAISAL_RukavinaMikusic2024.md`.
- AI assistance: the analysis code and drafts were produced with an AI research assistant under the owner's direction. Disclose this per journal policy.
