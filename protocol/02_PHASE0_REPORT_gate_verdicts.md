# Phase 0 report — alamandine/MrgD axis in AML

**Run date:** 20 November 2026 · **Protocol:** kickoff pack §3 plus frozen addendum v1 ({{artifact:44077ade-b8b8-47a7-bc84-85f92500b024}}, SHA-256 `466cf2ff…6844`, frozen before any MRGPRD/MAS1 value was viewed) · **OSF:** waived by owner; local hash-freeze instead.

## 1. Verdicts

| Gate | Pre-registered prediction | Result | Verdict |
|---|---|---|---|
| 0.1 Monomer presence | FAIL | MRGPRD "detected" in 0 of 6 eligible datasets (2 bulk, 4 scRNA) | **FAIL** |
| 0.2 Heteromer potential | UNCERTAIN | 0 co-detected cells among 289,966; every compartment untestable (< 20 MRGPRD+ cells); protein arm not assessable | **NOT PASSED** (untestable at transcript level) |
| 0.3 Dependency | FAIL | AML median Chronos: MRGPRD −0.10, MAS1 +0.05 (controls behave) | **FAIL** |
| 0.4 Clinical signal | WEAK PASS | Beat AML adjusted HR 1.03 per SD (0.90–1.17), q = 0.70 | **FAIL** — the prediction was wrong |

**Branch (kickoff §6):** 0.2 not passed and 0.4 failed → **Project C plus the Phase 0 negative.** Gate 0.4's drug-response component (§5.4) has **not been run**, because a safety-check approval is pending. It cannot change the 0.4 verdict, which is defined on overall survival alone.

## 2. Gate 0.1 — monomer presence

### 2a. Bulk (blasts; blast % ≥ 50) — {{artifact:c99575ba-187f-44a0-a882-bcc21d03932c}}

| Cohort | n | MRGPRD median TPM | > 1 TPM | Any read | ALB (negative control) median | CD34 / MPO median |
|---|---|---|---|---|---|---|
| TCGA-LAML (GDC STAR TPM) | 119 | 0.121 | 4.2% | 69.7% | 0.082 | 54 / 2,470 |
| Beat AML 2.0 (TPM from counts, GENCODE v36) | 245 | 0.000 | 0.0% | 1.2% | 0.025 | 22 / 1,709 |

- MRGPRD is at the level of the absent-gene control in both cohorts.
- TCGA's higher "any read" rate applies to ALB too (84%), so it reflects that pipeline's background floor, not MRGPRD signal.
- **Wording:** in bulk, MRGPRD is *not expressed above background* in AML blasts. Deep bulk RNA-seq is sensitive at the 1-TPM threshold, as the controls show, so here "not detected" can reasonably be read as "not expressed at ≥ 1 TPM".

### 2b. Single-cell — detection fraction = cells with ≥ 1 UMI; median across patients with ≥ 50 cells

| Dataset | Platform | Cells | MRGPRD+ cells | Compartments with median detection > 5% | Blast calibration: detection of 1–2 TPM genes |
|---|---|---|---|---|---|
| GSE116256 (van Galen) | Seq-Well | 41,090 | 1 | 0 | 0.66% → **underpowered** |
| GSE241989 | 10x | 114,254 | 2 | 0 | 1.69% (HSPC/blast-like) → **underpowered** |
| GSE207104, non-T-enriched | 10x | 51,880 | 3 | 0 | 0.44% → **underpowered** |
| GSE253355, normal marrow incl. 26,150 stromal cells | 10x | 82,742 | 0 | 0 | stroma uncalibrated (F.8) |

Tables: {{artifact:69af011d-b241-425b-9281-3e7dbd18ec77}}, {{artifact:073d406f-aaf5-42b6-9eb6-4c4ccf3896c4}}, {{artifact:7a55e4d0-90bc-48b1-992f-67f5447c332d}}, {{artifact:cece8741-a63a-4341-8b7b-dff52f934db5}}, {{artifact:4bd4e736-8afc-42f5-adde-65ca166b8262}}.

**Zero-inflation:** the MRGPRD zero fraction is ≥ 0.9992 in every compartment of every dataset; the minimum is GSE207104 erythroid, 1 positive cell of 1,320. Positive controls are detected in the expected compartments:
- CD34: HSPC median 42% (van Galen) and 49% (GSE253355);
- MPO: van Galen GMP-like 98%;
- CXCL12: GSE253355 stroma 73%.

The negative control ALB has a median per-patient detection ≤ 0.24% in every compartment of GSE116256 and GSE253355, and in the blast/HSPC/stroma compartments of GSE241989 and GSE207104. The other compartments of the two 10x AML datasets were not inspected for ALB.

### 2c. "Not detected" versus "not expressed" — central limitation, stated in the main text

- **By the pre-registered rule, the scRNA blast compartments are underpowered.** At realistic depth, even a gene at 1–2 TPM in bulk blasts is detected in < 5% of droplet or Seq-Well cells (0.44–1.69%). The 5%-of-cells threshold therefore only works for genes above roughly 5–10 TPM, so the scRNA and bulk thresholds are not equivalent. This is a design flaw in kickoff §3.1, and it is recorded here rather than corrected after the fact.
- **The calibration still rules out low-level expression in blasts.** If MRGPRD were at 1–2 TPM, about 1,000 of the 59,366 GSE241989 blast-like cells (0.0169 × 59,366) should show a UMI. **0 do.** Even at 0.5–1 TPM, about 490 would be expected. This is consistent with bulk: MRGPRD in AML blasts is well below 0.5 TPM.
- **Stroma:** 0 of 26,150 stromal cells in *normal* marrow. Sensitivity is uncalibrated, and **AML-associated stroma was not assessed by any eligible dataset.** The "stroma" clusters in the AML aspirate datasets scored at or below 0 on every marker set and are not credible stroma (§5).
- **Protein:** not assessed. This is the remaining gap: transcript non-detection does not exclude low-copy GPCR protein. Given bulk expression at the negative-control level, that possibility is remote but untested.

## 3. Gate 0.2 — heteromer potential — {{artifact:511b7376-d6e9-4f59-a4e7-8a3823cfa975}}

- Across 289,966 cells: 6 are MRGPRD+, 885 are MAS1+, and **0 are co-detected**.
- The largest number of MRGPRD+ cells in any single compartment is 2, against the pre-committed testability threshold of 20. Every compartment is therefore *untestable*, and the 10,000-permutation null was not run because its precondition was never met.
- MAS1 alone is sparse; its highest compartment is perivascular smooth muscle in normal marrow, at a median 2.5% per donor. MAS1 is not part of the 0.2 pass rule without MRGPRD, and it is reported here only as context.
- The kickoff (§2.4) puts the signalling unit in the MAS1–MRGPRD heteromer. These data contain no cells in which that heteromer could form at transcript level.

## 4. Gate 0.3 — dependency — {{artifact:236a5159-e60f-4675-a97c-f11f5b96ab63}}

| Gene | AML lines | AML median Chronos | Non-heme median | Difference | q (BH) | AML lines > 1 TPM | Pass |
|---|---|---|---|---|---|---|---|
| MRGPRD | 32 | −0.099 | −0.134 | +0.035 | 0.98 | 0 / 32 | No |
| MAS1 | 32 | +0.051 | +0.072 | −0.021 | 0.13 | 1 / 32 | No |
| MYB (positive control) | 32 | −1.259 | −0.248 | −1.011 | — | 32 / 32 | Yes |
| SPI1 (positive control) | 32 | −0.693 | −0.044 | −0.649 | — | 31 / 32 | Yes |
| ALB (negative control) | 32 | +0.079 | +0.098 | −0.020 | — | — | No |

- MAS1's unadjusted heme-versus-non-heme p of 0.014 reflects a median difference of about 0.03 Chronos units around zero. It is outside the decision rule and is not a dependency.
- The gene effects are near 0 because neither gene is expressed in AML lines; this is the expected null.

## 5. Gate 0.4 — clinical signal — {{artifact:dfcf34b5-ed54-4670-ac03-aa70cbc6d2c6}}

The model is Cox regression on OS: score (per SD) + age + risk (Beat AML ELN2017; TCGA RISK_MOLECULAR). Scores are GSVA 2.4.9 (Gaussian kcdf, maxDiff); all 8 genes were retained in both cohorts.

| Cohort | n / deaths | Score | Adjusted HR per SD (95% CI) | p | q (2 cohorts) | PH, score term | PH, global |
|---|---|---|---|---|---|---|---|
| Beat AML (primary) | 407 / 227 | RAS8 | 1.03 (0.90–1.17) | 0.70 | 0.70 | 0.24 | 0.017 |
| Beat AML | 407 / 227 | RAS7 (no MRGPRD) | 1.03 (0.91–1.18) | 0.62 | 0.62 | 0.28 | 0.019 |
| TCGA-LAML | 147 / 94 | RAS8 | 1.19 (0.95–1.48) | 0.14 | 0.28 | 0.95 | 0.003 |
| TCGA-LAML | 147 / 94 | RAS7 | 1.12 (0.90–1.40) | 0.31 | 0.61 | 0.66 | 0.002 |

- **FAIL:** the primary cohort shows no association, and TCGA does not reach significance.
- PH holds for the score term. The global PH violation comes from the other covariates, not the score. It is reported, and the estimates were not modified.
- **With versus without MRGPRD:** RAS8 and RAS7 correlate at r = 0.995 (Beat AML) and 0.919 (TCGA), and their HRs are essentially identical. MRGPRD was retained (non-constant), so this is a genuine finding, not an artefact of the gene being filtered out: **MRGPRD contributes nothing beyond the classical RAS genes.**
- **Not run:** Beat AML ex vivo drug AUC versus score, BH across the full panel. A safety-check approval is pending.

## 6. Deviations and analyst choices (complete list)

1. **OSF registration waived** by the owner; replaced by a local hash-freeze before any data were viewed.
2. **Transplant dropped** from the Cox models (F.7): neither cohort has transplant dates, so a baseline term would cause immortal-time bias. TCGA recorded transplant in 68 of 151 patients (descriptive only).
3. **Risk category:** ELN2017 / RISK_MOLECULAR instead of ELN 2022 (F.6).
4. **Complete-case exclusions in 0.4:** ambiguous ELN categories (Beat AML: FavorableOrIntermediate 8, IntermediateOrAdverse 2, MissingKaryo 4); TCGA RISK_MOLECULAR "N.D." (3); OS ≤ 0 days; unknown vital status. Not fixed by the protocol; declared.
5. **Beat AML TPM:** 11,365 of 63,677 Ensembl IDs had no GENCODE v36 match and were dropped from the TPM denominator (F.4).
6. **scRNA cell QC:** ≥ 200 genes per cell only. Not fixed by the protocol; this is the most lenient common choice.
7. **van Galen:** all timepoints pooled per patient. The D0-only sensitivity analysis gives the same result (0 detection, 1 MRGPRD+ cell).
8. **GSE253355:** the authors' processed object **lacks MRGPRD** (a gene filter). Counts were rebuilt from the RAW matrices (33,538 features), and all 82,742 annotated cells were matched.
9. **Leiden clustering (GSE241989/GSE207104):** an unweighted, symmetrised k = 15 neighbour graph was built directly, because scanpy/UMAP/numba versions in the environment were incompatible. A `getnnz` shim was added for scipy. Both are implementation changes only.
10. **Marker-based "stroma" clusters in the AML aspirate datasets are not credible**: every marker set scored ≤ 0.14 and CXCL12 detection was 0.1%. They are reported, but not interpreted as stroma.
11. **GSVA runtime:** BiocSingular's compiled library crashes this conda R build (Mingw-w64 pseudo-relocation). It was replaced by a stand-in exporting only `runExactSVD`, which errors if called and is used only by PLAGE. GSVA was run in a separate Rscript process.
12. **Exploratory, not pre-specified:** HSPC calibration against the blast reference in GSE253355 (1–2 TPM genes detected in 4.6% of HSPCs). Labelled exploratory; not used for any verdict.
13. **Caught bug, corrected before any verdict:** in gate 0.3, `pd.read_csv(usecols=…)` returned columns in file order, which mislabelled the genes. The positive and negative controls exposed it; the table was regenerated with verified labels.
14. **Bgee** dropped (truncated exports). **BloodSpot:** MAS1 only; context only.

## 7. What would change the conclusion

- Protein-level MrgD in AML blasts or AML-associated niche cells, by validated antibody or targeted MS.
- Transcript detection in a targeted assay (smFISH, or deep full-length scRNA of sorted blasts and AML stroma) at a calibrated sensitivity below 0.5 TPM.
- Heteromer evidence (co-IP / BRET) in primary AML cells, rather than THP-1 alone.
