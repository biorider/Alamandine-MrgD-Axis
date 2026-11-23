"""Step 3a: derive approximate ELN 2022 risk (scripts/eln2022.py) for Beat AML and TCGA-LAML. Sensitivity analysis only."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, re, sys, json, requests, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from config import OPEN, INTERIM, TABLES
import eln2022 as E22
# Beat AML
cl = pd.read_excel(os.path.join(OPEN, 'beataml_wv1to4_clinical.xlsx'), sheet_name='summary')
ids = pd.read_csv(os.path.join(INTERIM, 'g04_beat_clin.csv')).dbgap_rnaseq_sample
sub = cl[cl.dbgap_rnaseq_sample.isin(ids)]
bm = pd.read_csv(os.path.join(OPEN, 'beataml_wes_wv1to4_mutations_dbgap.txt'), sep='\t', usecols=['dbgap_sample_id', 'symbol', 'variant_classification', 't_vaf'])
bm = bm[bm.variant_classification.map(lambda v: any(t in E22.PROTEIN_ALTERING for t in str(v).split('&')))]
rows = []
for _, r in sub.iterrows():
    fus = set(str(r.consensusAMLFusions).upper().split('|')) if pd.notna(r.consensusAMLFusions) else set()
    wes = pd.notna(r.dbgap_dnaseq_sample); m = bm[bm.dbgap_sample_id == r.dbgap_dnaseq_sample] if wes else bm.iloc[0:0]
    lab = E22.classify(E22.parse_karyotype(r.karyotype), fus, r.NPM1 == 'positive', r['FLT3-ITD'] == 'positive', r.CEBPA_Biallelic == 'bi',
                       ((m.symbol == 'TP53') & (m.t_vaf >= 0.10)).any(), m.symbol.isin(E22.MR).any()) if wes else 'Unclassifiable (no WES)'
    rows.append(dict(id=r.dbgap_rnaseq_sample, ELN2022_derived=lab, ELN2017=r.ELN2017))
B = pd.DataFrame(rows); B.to_csv(os.path.join(TABLES, 'step3_beat_eln2022.csv'), index=False)
print(pd.crosstab(B.ELN2017, B.ELN2022_derived, margins=True))
# TCGA-LAML (cBioPortal laml_tcga_pub mutations; cytogenetics + inferred rearrangements)
API = 'https://www.cbioportal.org/api'; G = ['NPM1', 'FLT3', 'CEBPA', 'TP53'] + E22.MR
ent = {g['hugoGeneSymbol']: g['entrezGeneId'] for g in requests.post(f'{API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', json=G, timeout=60).json()}
mu = requests.post(f'{API}/molecular-profiles/laml_tcga_pub_mutations/mutations/fetch?projection=DETAILED',
                   json={'entrezGeneIds': list(ent.values()), 'sampleListId': 'laml_tcga_pub_all'}, timeout=120).json()
tm = pd.DataFrame([dict(patient=m['patientId'], gene=m['gene']['hugoGeneSymbol'], ptype=m.get('mutationType'), pchange=m.get('proteinChange'),
                        alt=m.get('tumorAltCount'), ref=m.get('tumorRefCount'), start=m.get('proteinPosStart')) for m in mu])
tm['vaf'] = tm.alt / (tm.alt + tm.ref); tm.to_csv(os.path.join(TABLES, 'step3_tcga_mutations.csv'), index=False)
tm = tm[tm.ptype.isin(E22.PROTEIN_ALTERING)]
T = pd.read_csv(os.path.join(INTERIM, 'tcga_laml_pub_clinical_cbioportal.csv'), index_col=0)
rows = []
for pid in pd.read_csv(os.path.join(INTERIM, 'g04_tcga_clin.csv')).id:
    if pid not in T.index: rows.append(dict(id=pid, ELN2022_derived='Unclassifiable (no record)')); continue
    r = T.loc[pid]; kar = str(r.CYTOGENETICS) if pd.notna(r.CYTOGENETICS) else None
    if kar is not None and pd.notna(r.INFERRED_GENOMIC_REARRANGEMENT): kar += ',' + str(r.INFERRED_GENOMIC_REARRANGEMENT)
    m = tm[tm.patient == pid]; fl = m[m.gene == 'FLT3']
    itd = ((fl.ptype == 'In_Frame_Ins') & fl.start.between(570, 650)).any() or (fl.pchange.str.contains('dup|ins', na=False) & fl.start.between(570, 650)).any()
    ceb = ((m.gene == 'CEBPA') & m.ptype.isin(['In_Frame_Ins', 'In_Frame_Del', 'Missense_Mutation']) & m.start.between(278, 345)).any()
    rows.append(dict(id=pid, ELN2022_derived=E22.classify(E22.parse_karyotype(kar), set(), (m.gene == 'NPM1').any(), itd, ceb,
                     ((m.gene == 'TP53') & (m.vaf >= 0.10)).any(), m.gene.isin(E22.MR).any()), RISK_MOLECULAR=r.RISK_MOLECULAR))
Tt = pd.DataFrame(rows); Tt.to_csv(os.path.join(TABLES, 'step3_tcga_eln2022.csv'), index=False)
print(pd.crosstab(Tt.RISK_MOLECULAR, Tt.ELN2022_derived, margins=True))
