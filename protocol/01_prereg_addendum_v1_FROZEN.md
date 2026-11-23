# Phase 0 Pre-registration Addendum — Alamandine/MrgD axis in AML

**Status:** FROZEN v1, 20 Nov 2026. OSF waived by owner; local hash-frozen. Section F governs.
**Drafted:** 20 November 2026. **No MRGPRD/MAS1 value from any dataset had been viewed when this was written.** Only file headers, row counts, gene-presence checks and metadata were read.

This addendum makes Section 3 of the kickoff pack operational. Where a rule below conflicts with Section 3, the rule below governs, because it is more specific. Every item is marked **[DEFAULT — owner to confirm]** until the owner edits it.

## A. Guiding principle for open choices

The predicted outcome for 0.1 and 0.3 is FAIL. Where a choice is otherwise arbitrary, the **decisional** analysis uses the setting most *lenient to the hypothesis*. A FAIL under lenient settings is robust. A PASS under lenient settings counts only after a pre-specified confirmation (listed per gate). This prevents both rescuing a negative and dismissing a positive.

## B. Inputs (frozen)

| Source | File(s) | Use | Counts toward the "≥2 independent datasets" rule in 0.1? |
|---|---|---|---|
| TCGA-LAML | GDC STAR augmented counts, `tpm_unstranded` column; one primary blood-derived cancer sample per patient | 0.1 bulk (blasts), 0.4 | Yes |
| Beat AML 2.0 | biodev/beataml2.0_data: `counts`, `norm_exp`, `probit_curve_fits_v4`, `wv1to4_clinical.xlsx`, `sample_mapping.xlsx` (SHA-256 in `download_manifest.csv`) | 0.1 bulk (blasts), 0.4 | Yes |
| van Galen 2019 | GEO GSE116256 RAW (Seq-Well) | 0.1 scRNA (blasts, HSPC, immune), 0.2 | Yes |
| Bandyopadhyay 2024 | GEO GSE253355 (normal marrow incl. stromal niche) | 0.1 scRNA (stroma, normal HSPC), 0.2 | Yes |
| GSE207104, GSE241989, GSE247205 | GEO RAW | 0.1 scRNA, 0.2 | **[OWNER: include or exclude now — cannot be decided after viewing]** |
| DepMap (release: **[OWNER: fill in]**) | `CRISPRGeneEffect.csv`, `CRISPRGeneDependency.csv`, `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` (unstranded) | 0.3; cell-line expression context | No — cell lines are not blasts, HSPCs or stroma |
| HPA, Bgee, BloodSpot | as in pack | Context | No — supporting only |

SHA-256 of all user-supplied files: `user_data_manifest.csv`. Open downloads (Beat AML 2.0, GSE253355): 8 files retrieved 20 Nov 2026; SHA-256 per file in `download_manifest.csv` (attach both manifests to the OSF registration).

**Known input defects, stated before analysis:**
1. The Bgee exports are truncated at exactly 1,000 rows (an export cap). The single-cell exports contain no marrow or blood libraries in those rows. **Owner decision:** re-export fully, or drop Bgee. A truncated export must not be cited as evidence of absence.
2. BloodSpot: a MAS1 export exists (`MAS1_log2.csv`, microarray log2), but there is no MRGPRD export. `mas1.pdf` is a browser print whose title is "about:blank"; check it renders the plot. Microarray data cannot satisfy a TPM or %-cells rule, so BloodSpot is supporting only.
3. The DepMap expression file uses a newer name than the pack. It is profile-level (2,314 rows), so models are de-duplicated with `IsDefaultEntryForModel`. The unstranded file is used because it matches the file the pack names; the stranded file is not analysed.
4. The existing ILC-project Seurat objects are **not** used. Counts are rebuilt from the RAW files with no `min.cells`/`min.features` gene filtering. Before any result is produced, the script asserts that MRGPRD and MAS1 are in each feature list. If a gene is absent from a reference annotation, the result is reported as "not quantifiable", not as zero.

## C. Gate operationalisation

### Gate 0.1 — monomer presence

- **Bulk [DEFAULT]:** a cohort is "detected" if **median** TPM of MRGPRD > 1 among eligible samples. Eligible means initial-diagnosis, adult, one specimen per patient (bone marrow preferred over peripheral blood), and blast % ≥ 50 where recorded. Bulk speaks only to the *blast* compartment; stroma and HSPC are not assessable from bulk aspirates. Reported but non-decisional: fraction of samples > 1 TPM, and the all-samples median.
  - Beat AML TPM is computed from `counts` using GENCODE gene lengths **[OWNER: name the GENCODE version]**. The `norm_exp` file is not TPM and is not used for thresholding.
- **scRNA [DEFAULT]:** detection means ≥ 1 UMI in a cell. The compartment detection fraction is the **median across patients** with ≥ 50 cells in that compartment. A compartment is "detected" if that median is > 5%. No ambient-RNA correction is applied in the decisional analysis, because the correction can only lower detection (lenient setting, per §A). Compartments: blast/malignant, HSPC, stroma (MSC/osteolineage/endothelial), macrophage/monocyte, T/NK. Labels come from each dataset's published annotation; nothing is re-clustered for gate decisions.
- **Controls, reported in the same table:** positive controls CD34 (HSPC, blasts) and MPO (myeloid blasts); stromal positive control CXCL12; negative control **[DEFAULT: ALB]**, absent from marrow.
- **Sensitivity calibration (decides how "not detected" is worded):** within each scRNA dataset and compartment, bin genes by bulk expression. Reference bulk: Beat AML median TPM for blasts; for stroma, **[OWNER: specify a bulk BM-MSC reference or declare uncalibrated]**. Estimate the per-cell detection fraction expected for a gene at 1–2 TPM. If that expected fraction is < 5%, the dataset is declared **underpowered** for the 0.1 threshold in that compartment. A non-detection there is then reported as *not detected; assay insensitive at this abundance* and does not count toward FAIL or PASS.
- **PASS:** "detected" in at least one compartment (blasts, HSPC or stroma) in ≥ 2 eligible datasets.
- **Confirmation required for a PASS:** the result must survive ambient correction (SoupX, automatic contamination estimate) and doublet removal (scDblFinder) in the same datasets.

### Gate 0.2 — heteromer potential

- **Statistic:** the number of cells with ≥ 1 UMI of both MAS1 and MRGPRD, within one compartment of one dataset.
- **Null:** 10,000 permutations of MRGPRD values across cells *within each patient × compartment*. This keeps each gene's per-patient detection rate and breaks only the pairing. One-sided p = (1 + #null ≥ observed)/(1 + 10,000).
- **PASS [DEFAULT]:** in ≥ 1 dataset and compartment (blasts, HSPC, stroma or macrophage), p < 0.05 after BH correction across all tested compartment × dataset cells, **and** ≥ 20 co-detected cells from ≥ 3 patients. The minimum count is there so that a handful of cells cannot decide the gate.
- **Power statement, pre-committed:** if fewer than 20 cells detect MRGPRD in a compartment, that compartment is reported as *untestable*, not as a FAIL.
- **Protein arm [DEFAULT]:** declared **not assessable** in this dry-lab Phase 0. No AML proteomics dataset is pre-specified, and MS and antibody detection of MrgD are unreliable. The owner may instead name specific datasets here before the timestamp.

### Gate 0.3 — dependency

- Lines are those with `OncotreePrimaryDisease == "Acute Myeloid Leukemia"` (the AML set). Heme set = `OncotreeLineage ∈ {Myeloid, Lymphoid}`. Comparator = all other lines with CRISPR data.
- **PASS [DEFAULT]** for a gene (MRGPRD, MAS1) requires all of:
  1. median Chronos gene effect in AML lines ≤ −0.5;
  2. the AML median minus the comparator median ≤ −0.3;
  3. one-sided Wilcoxon AML < comparator, with BH-FDR < 0.05 across the 2 genes.
- Reported alongside: the same statistics for the heme set; the fraction of AML lines with CRISPRGeneDependency > 0.5; MRGPRD and MAS1 expression in the same lines. Gene effects in non-expressing lines are expected to be about 0 and are interpreted as such.
- Positive control: MYB and SPI1 (known AML-selective dependencies). Negative control: ALB.

### Gate 0.4 — clinical signal

- **Score:** GSVA (Bioconductor `GSVA`, `gsvaParam`, `kcdf="Gaussian"`, `minSize=1`, `maxDiff=TRUE`), fitted per cohort on log2(TPM+1) over the eight Section 3.2 genes. This is one unsigned score, as originally registered. The pack's phrase "single-sample GSVA" is corrected here: GSVA is cohort-relative. Pre-specified fallback (§3.3): the mean z-score. Genes retained after GSVA's internal filtering are logged.
- **With vs without MRGPRD:** if MRGPRD is dropped by filtering or has zero variance, the comparison is reported as **uninformative (MRGPRD absent)**. It does *not* count as "MRGPRD contributes nothing".
- **Model:** Cox regression of overall survival on the score (per 1 SD), age (continuous) and risk category.
  - **Risk category [DEFAULT]:** ELN as recorded in each cohort (Beat AML: the ELN field in `wv1to4_clinical.xlsx`; TCGA-LAML: the cytogenetic/molecular risk field). This is a declared deviation from "ELN 2022". It avoids a derived classification with its own analyst choices.
  - **Transplant [DEFAULT]:** included as a **time-dependent** covariate if transplant dates exist; otherwise **dropped**, with the omission declared. A baseline (ever/never) transplant covariate is not allowed because it causes immortal-time bias.
  - PH assumption: Schoenfeld residual test (`cox.zph`), global and per term, reported. If PH is violated for the score, the primary estimate stays as fitted and a time-split HR is reported alongside.
- **PASS [DEFAULT]:** Beat AML (primary, larger cohort) score HR with BH-FDR < 0.05, **and** the TCGA-LAML HR in the same direction (replication of direction, not significance). FDR family = the two cohort tests.
- **Drug response:** Spearman correlation of the score with ex vivo AUC for **every** inhibitor in `probit_curve_fits_v4` tested in ≥ 100 Beat AML initial-diagnosis patients. BH across the full panel; the panel size is reported. Cytarabine, daunorubicin and venetoclax are highlighted but not tested separately.

## D. Order and locking

Gates run in order 0.1 → 0.2 → 0.3 → 0.4. Each gate's outputs (tables, code version, input hashes) are saved and its verdict written before the next gate starts. No threshold, dataset, covariate or method is changed after a gate's output has been viewed. Any deviation is logged with its reason and stated in the write-up.

## F. Frozen v1 — resolutions (20 Nov 2026, before any target value was viewed)

The owner waived OSF registration ("only trying to prove our hunch"). In its place, this file is frozen locally as v1: SHA-256 recorded, saved as a versioned artifact before analysis. All [DEFAULT] items above are adopted as written, except where this section overrides them.

1. **DepMap release:** not recorded by the owner. Files were downloaded 20 Nov 2026; identity is fixed by the SHA-256 values in `user_data_manifest.csv`.
2. **scRNA eligibility (decided on study design, from GEO metadata only):**
   - GSE116256 (AML bone-marrow aspirates plus normal marrow, Seq-Well): INCLUDE. Uses the published cell-type and malignant labels from the `anno` files.
   - GSE241989 (primary AML cells from marrow/blood, 10x): INCLUDE.
   - GSE207104 (paediatric post-transplant AML marrow, 10x): INCLUDE the **non-T-enriched samples only**, because T-cell enrichment depletes blasts.
   - GSE253355 (normal human marrow including stroma, 10x): INCLUDE. Uses the authors' annotated Seurat objects (whole atlas plus MSC subset).
   - GSE247205 (expanded or blood ILC2s): EXCLUDE — not a marrow compartment.
3. **Compartment labels for datasets without published labels (GSE241989, GSE207104):** fixed scanpy pipeline.
   - Processing: normalize_total 1e4 → log1p → 2,000 HVG → PCA (30 PCs) → kNN (15) → Leiden at resolution 1.0, with seed 0.
   - Each cluster is assigned to the compartment whose marker set has the highest mean z-scored expression. Marker sets:
     - HSPC/blast-like: CD34, PROM1, KIT, SPINK2
     - promyelocyte/GMP-like: MPO, ELANE, AZU1, PRTN3
     - mono/macro: CD14, LYZ, CD68, CSF1R, FCGR3A
     - T/NK: CD3D, CD3E, NKG7, GNLY
     - B: MS4A1, CD79A
     - erythroid: HBB, GYPA, CA1
     - stroma: CXCL12, LEPR, COL1A1, VCAM1
   - "Blast compartment" for the 0.1 decision = HSPC/blast-like ∪ promyelocyte/GMP-like clusters in AML samples. Limitation, stated now: monocytic AML blasts can fall in the mono/macro compartment.
   - MRGPRD and MAS1 are excluded from HVG selection and from all marker sets, so they cannot influence the labels.
4. **Beat AML TPM:** gene length = union of exons per gene in GENCODE v36 (the annotation GDC uses for the TCGA STAR TPMs), matched on version-stripped Ensembl ID. Genes without a match are dropped from both the TPM numerator and the denominator; the number dropped is reported. Declared deviation: Beat AML was aligned to an older assembly/annotation.
5. **Beat AML cohort:** `diseaseStageAtSpecimenCollection` = Initial Diagnosis, age ≥ 18, `dxAtInclusion` AML, one RNA-seq specimen per patient (bone marrow preferred over peripheral blood, then earliest). For the 0.1 blast restriction, blast % = `%.Blasts.in.BM` for marrow specimens and `%.Blasts.in.PB` for blood specimens.
6. **Risk covariate:**
   - Beat AML: `ELN2017`.
   - TCGA-LAML: `RISK_MOLECULAR` from cBioPortal `laml_tcga_pub` (Ley 2013).
   - Declared deviation from ELN 2022.
7. **Transplant:** no transplant dates exist in either cohort. Transplant is dropped from both Cox models, and `TRANSPLANT_TYPE` (TCGA) is reported descriptively only.
8. **Stromal calibration:** no bulk stromal reference is pre-specified, so the stroma compartment is declared **uncalibrated**. A stromal non-detection is worded "not detected; sensitivity not calibrated".
9. **Bgee:** dropped from evidence; the exports are truncated. BloodSpot MAS1 export: context only.
10. **0.2 null:** 10,000 permutations; if runtime forces fewer, the number used is reported. No other change is allowed.

## E. Items the owner had to fill (resolved in Section F)

1. DepMap release name.
2. Include or exclude GSE207104 / GSE241989 / GSE247205.
3. GENCODE version for the Beat AML TPM computation.
4. Stromal bulk reference for calibration, or declare stroma uncalibrated.
5. Bgee: re-export fully or drop.
6. Accept, or edit, every [DEFAULT] above.
