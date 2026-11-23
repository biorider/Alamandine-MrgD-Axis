# Phase 0 Quick-Start Pack — Alamandine/MrgD Axis in Leukemia

**Owner:** Uğur Toprak (ToprakLab-LHU) · **Prepared:** 20 November 2026
**Companion doc:** [Critical Appraisal and Revised Protocol](https://claude.ai/artifact/UC5ZWynSUhV4huBuzg2BKQ)

---

## 0. How to use this file

This pack is designed so a fresh Claude session can start productive work in one message, without re-deriving what has already been verified.

1. Copy **Section 1** into a new Claude session (with Python/R enabled).
2. Attach this file to that session.
3. Work through the gates in order. Do not skip to Phase 1.

**Decision already made:** run Phase 0 to completion regardless of the low prior. A rigorous negative on a hypothesis published three months ago is a publishable, citable result and protects against committing years to a dead target. Projects A/B/C in Section 6 are held as backups.

---

## 1. Kickoff prompt — copy this

```markdown
I'm running a pre-registered go/no-go analysis (Phase 0) on whether the
alamandine/MrgD axis is a credible target in acute myeloid leukemia.
The attached file (phase0-kickoff-alamandine-mrgd.md) contains the full
context pack, pre-registered decision rules, and analysis spec. Read it first.

IMPORTANT CONTEXT: this hypothesis has already failed three preliminary
checks (see Section 2). I am running Phase 0 anyway, properly, because a
rigorous negative is the deliverable. Do not try to rescue the hypothesis.
Your job is to help me test it fairly and report what is there.

Act as a skeptical computational co-researcher. Specifically:

1. Treat Section 3's decision rules as pre-registered and binding. If a
   result misses a threshold, say it missed. Do not reach for a secondary
   analysis that would have passed.
2. Flag every analyst degree of freedom you exercise (normalization choice,
   detection threshold, covariate set) and tell me when a choice could
   plausibly flip a result.
3. For single-cell work, report detection fractions and zero-inflation
   explicitly. A GPCR at 0.3 TPM is a negative, not a trend.
4. Distinguish clearly between "not detected" and "not expressed."
   Low-abundance GPCRs are undersampled by droplet methods, and this
   limitation belongs in the write-up, not in a footnote.

Start by telling me which gates you can begin immediately with no data from
me (0.1 partially, via public browser-accessible resources), and exactly
which files I need to download for the rest. Give me a download checklist
with URLs and expected file names before writing any analysis code.

Constraints: dry lab only. I can download from DepMap, cBioPortal, GDC, GEO
and Beat AML; you cannot authenticate to those. You write and debug the
pipelines, I supply the files.
```

---

## 2. Context pack — verified, do not re-derive

All figures queried 20 November 2026 and reproducible.

### 2.1 Expression (Human Protein Atlas, ENSG00000172938)

| Category | Value |
| --- | --- |
| RNA tissue specificity | Not detected |
| Immune cell specificity | Not detected; highest 0.0 nTPM (myeloid DC) |
| Single cell type specificity | Not detected |
| Blood protein (immunoassay, MS) | Not detected |
| TCGA RNA cancer | Not detected |
| Cancer cell lines | Not detected; highest 0.2 nTPM (skin) |
| Highest tissue overall | Intestine, 0.6 nTPM |

Jurkat, K-562, NB4, REH, HEL and THP-1 all appear in HPA's cell line list with no detectable value.

### 2.2 Literature volume (PubMed)

| Query | Hits |
| --- | --- |
| `alamandine AND leukemia` | **0** |
| `MRGPRD AND (hematopoietic OR leukemia OR bone marrow)` | **2** |
| `alamandine AND MrgD AND cancer AND proliferation AND tumor` | 2 |

The two hematopoietic hits:

- **Jimenez-Andrade 2010, Bone** (doi:10.1016/j.bone.2009.09.013) — bone marrow, mineralized bone and periosteum receive SP+/CGRP+ but **not** Mrgprd+ innervation. Evidence *against* marrow MrgD.
- **Lv 2026, J Pathol** (doi:10.1002/path.70125) — MrgprD present in colonic macrophages and BMDMs; agonist activation **enhanced** pro-inflammatory NF-κB/M1 polarization; myeloid-specific KO was protective. Points *opposite* to the review's protective framing.

### 2.3 Chemical tractability (ChEMBL CHEMBL4523899 / UniProt Q8TDS7)

Exactly 2 bioactivity records, both orphan GPCRome panel screens (UCSF924, 2021; CHEMBL5724552, 2024). Neither reports IC50, EC50 or pChEMBL. No ligand series, no tool antagonist, no selectivity data.

Structures do exist and are solid: **PDB 7Y12, 7Y13, 7Y14, 7Y15** (Suzuki 2022, MrgD–Gi cryo-EM) plus **9DQH, 9DQJ**.

### 2.4 The finding that reframes the target

**Rukavina Mikusic 2024**, Biochem Pharmacol 229:116480 (doi:10.1016/j.bcp.2024.116480), in human macrophages:

- MasR and MrgD form a **constitutive heteromer**, not induced by alamandine.
- Alamandine's IL-6 reduction was abolished by MasR blockade and MasR knockdown.
- Akt/ERK1/2 activation occurred **only** in heteromer-expressing cells.
- Ang-(1–7) was antiproliferative in MasR+ but **not** MrgD+ cells.
- The antiproliferative effect required the heteromer.

**Critical detail:** those antiproliferative experiments were run in **THP-1 cells — an acute monocytic leukemia (AML FAB M5) line.** Alamandine was antiproliferative in THP-1, MasR-dependent. This is the only existing leukemia-relevant functional data point in the entire literature, and no one has framed it as such.

**Consequence for Phase 0:** test *co-expression* of MAS1 and MRGPRD, not MRGPRD alone. A selective MrgD agonist would be expected to fail; the signaling unit is the heteromer.

### 2.5 Source integrity — do not inherit the review's citations

Verified misattributions in Yolcu & Haznedaroglu 2026:

| Claim | Cites | Actually is |
| --- | --- | --- |
| "Antiproliferative effect requires both receptors" | ref 57 | Zhao 2023, cardiac nuclear import of MrgD. Correct source is ref 22 |
| Ang A binds MrgD → pro-proliferative effects | ref 23 | Suzuki 2022 cryo-EM structural paper |
| MrgD directly antagonizes AT1R (rebuts their own safety concern) | ref 23 | Same structural paper |
| MrgD expressed "most importantly, in the stem and progenitor cells of the bone marrow" | **none** | Uncited; contradicted by HPA and Jimenez-Andrade 2010 |
| Gut–bone marrow axis section | refs 20–22 | Partly a bioRxiv preprint; murine; about aging |

Re-verify any reference before citing it. A sample contained this many problems; assume more exist among the 73.

---

## 3. Pre-registration block — timestamp this on OSF before running anything

Register these four items verbatim: the gene list, the scoring method, the covariates, and the decision rules.

### 3.1 Decision rules (binding)

| Gate | Pass criterion | Action if failed |
| --- | --- | --- |
| **0.1 Monomer presence** | MRGPRD detected (>1 TPM bulk, or >5% of cells in scRNA) in blasts, HSPCs, or marrow stroma in ≥2 independent datasets | Original hypothesis (project D) is dead. Proceed to 0.2 |
| **0.2 Heteromer potential** | MAS1 and MRGPRD co-detected in the same cell population in ≥1 scRNA dataset, **or** MAS1 detected in blasts with MRGPRD confirmed at protein level | Project A reduces to a MasR/Ang-(1–7) project — still viable, different paper |
| **0.3 Dependency** | DepMap CRISPR gene effect for MRGPRD or MAS1 shows selective dependency in heme lineage lines | Not fatal alone (GPCRs are often non-essential in vitro); removes a supporting argument |
| **0.4 Clinical signal** | RAS score associates with overall survival in Beat AML or TCGA-LAML at FDR < 0.05, adjusted for age and ELN risk | Project B weakens substantially; report as negative |

### 3.2 Gene list (fixed)

`AGT, REN, ACE, ACE2, AGTR1, AGTR2, MAS1, MRGPRD`

### 3.3 Scoring method (fixed, chosen before seeing outcomes)

Single-sample GSVA on log2(TPM+1). Pre-specified fallback if GSVA fails to converge: z-score mean across the eight genes. **No other method may be substituted after seeing results.**

### 3.4 Covariates (fixed)

Age at diagnosis, ELN 2022 risk category, transplant status. Cox proportional hazards; test the PH assumption and report it.

### 3.5 Recorded predictions (make before looking)

Prediction is part of the pre-registration. Stated in advance:

- **0.1 — predicted FAIL.** HPA shows non-detection across every relevant category.
- **0.2 — predicted UNCERTAIN.** This is the informative gate. MAS1 has broader hematopoietic expression than MRGPRD.
- **0.3 — predicted FAIL.** No prior reason to expect dependency.
- **0.4 — predicted WEAK PASS.** Driven by ACE/AGTR1, not by MRGPRD.

Recording these is what lets you tell afterwards whether you learned something or rationalized.

### 3.6 Multiplicity

Eight genes × several cohorts × ~120 Beat AML drugs is enough analyst freedom that something will reach p<0.05. Apply Benjamini–Hochberg across the **full** drug panel, not the three agents of interest. Report the full panel result.

---

## 4. Data acquisition checklist

Ordered by effort. Exhaust the free/instant tier before requesting anything.

- [ ] **Human Protein Atlas** — done, Section 2.1. No download.
- [ ] **BloodSpot** (servers.binf.ku.dk/bloodspot) — browser. Check MRGPRD and MAS1 side by side across normal and malignant hematopoiesis. Fastest 0.1 evidence.
- [ ] **Bgee** — browser. Independent normalization, useful as a second dataset for 0.1.
- [ ] **DepMap public 25Q2 or later** (depmap.org/portal/download) — need `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `CRISPRGeneEffect.csv`, `Model.csv`. Serves 0.1 and 0.3.
- [ ] **TCGA-LAML** via cBioPortal (`laml_tcga_pan_can_atlas_2018`) — RNA-seq + clinical. 151 samples. Serves 0.1 and 0.4.
- [ ] **Beat AML** (Vizome / dbGaP phs001657) — RNA-seq + ex vivo drug response for ~120 agents. **Submit the data request early**; this is the rate-limiting step and the highest-value dataset.
- [ ] **van Galen 2019** — GEO `GSE116256`. scRNA of AML + normal marrow. **This is where gate 0.2 is actually decided.**
- [ ] **Human Cell Atlas bone marrow** — second scRNA dataset for 0.1/0.2 replication.

---

## 5. Analysis specification

### 5.1 Gate 0.1 — detection

Report MRGPRD and MAS1 TPM distributions across blast, HSPC, stromal, macrophage and T/NK compartments. Required outputs:

- Detection fraction (% cells with ≥1 count) per compartment, not just mean expression.
- Zero-inflation explicitly stated.
- Side-by-side positive control (e.g. `CD34`, `MPO`) and a known-absent negative control, to demonstrate the pipeline can detect what is there.

The control panel matters: without it, a negative is indistinguishable from a broken pipeline, and a reviewer will say so.

### 5.2 Gate 0.2 — co-expression

Per-cell co-detection of MAS1 and MRGPRD, tested against a null of independent expression at matched detection rates.

**Power warning, stated in advance:** if MRGPRD detection is near zero, co-detection is statistically fragile. Report the power limitation rather than presenting a fraction of a percent as a finding. If the observed co-detection rate is below what the null would produce by chance, say so plainly.

### 5.3 Gate 0.3 — dependency

DepMap `CRISPRGeneEffect` for MRGPRD and MAS1. Compare heme lineage lines against all others. Report effect sizes; a gene effect near 0 across all lines is a clean negative.

### 5.4 Gate 0.4 — clinical signal

RAS score (Section 3.3) vs overall survival in Beat AML and TCGA-LAML, adjusted per Section 3.4. Then Beat AML ex vivo AUC for cytarabine, daunorubicin and venetoclax vs the score, FDR-controlled across the full panel.

Run the score both **with and without** MRGPRD. If the signal is identical either way, that is itself the finding — it means the axis contributes nothing beyond classical RAS.

---

## 6. Backup projects (held, not abandoned)

| | Project | Testable claim | First output |
| --- | --- | --- | --- |
| **A** | MasR–MrgD heteromer in AML | The alamandine response in myeloid leukemia is carried by the heteromer; heteromer abundance predicts sensitivity | Computational + THP-1 replication; heteromer-biased ligand design |
| **B** | BM RAS signature biomarker | A multi-gene RAS score stratifies ELN risk and predicts drug response in AML | Bioinformatics paper from Beat AML + TCGA-LAML |
| **C** | Scoping review / evidence map | The marrow alamandine/MrgD literature does not support the therapeutic claims made for it | PRISMA-ScR review |

Gate 0.2 is the branch point: pass → project A; fail with 0.4 passing → project B; both fail → project C plus the Phase 0 negative.

Note that gate 0.4's work **is** project B. The pipelines are shared, not duplicated.

---

## 7. Write-up template for the negative

Prepare this now, so the result is not orphaned if the gates fail.

**Working title:** *Expression-level evaluation of the alamandine/MrgD axis in acute myeloid leukemia: a pre-registered negative*

Structure:

1. The 2026 review proposes the axis as a marrow paracrine regulator and therapeutic target.
2. The central expression claim is uncited and testable.
3. Pre-registered gates, with the OSF timestamp.
4. Results across HPA, DepMap, TCGA-LAML, Beat AML and two scRNA atlases, with detection controls.
5. Conclusion scoped honestly: *not detected at transcript level in these compartments and datasets*, with droplet undersampling of low-abundance GPCRs stated as a limitation, and protein-level assessment named as the remaining gap.
6. What would change the conclusion: protein-level detection in blasts or niche, or heteromer-level evidence.

Target journals where a well-powered negative is publishable: *Scientific Reports*, *PLOS ONE*, *Frontiers in Oncology*, or a Letter to *J Renin Angiotensin Aldosterone Syst* itself.

---

## 8. Things not to delegate

- **Verify citations yourself before publishing.** A sample of the review's 73 references contained the errors in Section 2.5; assume more.
- **Read the Lv 2026 and Larrinaga 2025 figures yourself.** Whether high MrgD in ccRCC is a driver or a counter-regulatory response to aggressive disease changes the direction argument, and abstracts cannot settle it.
- **Decide the professional question around project C** — publishing an appraisal touching a Hacettepe paper is a judgment only you can make. A neutral evidence map achieves most of the scientific purpose.
- **Keep the kill criteria.** The main risk to this project is no longer scientific; it is the temptation to soften a negative after investing effort in it.
