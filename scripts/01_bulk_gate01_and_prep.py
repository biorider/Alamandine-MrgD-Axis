"""Gate 0.1 (bulk arm) + preparation of gate 0.4 inputs.
TCGA-LAML: GDC STAR tpm_unstranded (GENCODE v36). Beat AML 2.0: TPM from raw counts with GENCODE v36 union-exon lengths.
Eligibility (frozen F.5): initial diagnosis, age>=18, AML at inclusion, one RNA-seq specimen per patient (BM > PB, earliest);
decisional subset = blast% >= 50 (BM blasts for BM specimens, PB blasts for PB specimens). TCGA blast% = BM_BLAST_PERCENTAGE (laml_tcga_pub)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, json, numpy as np, pandas as pd, requests
from config import OPEN, INTERIM, TABLES, TCGA_STAR, TCGA_FILES
from common import gene_lengths, tpm
GENES = ['MRGPRD', 'MAS1', 'AGT', 'REN', 'ACE', 'ACE2', 'AGTR1', 'AGTR2', 'CD34', 'MPO', 'ALB']

# ---- TCGA-LAML
fj = json.load(open(TCGA_FILES)); mats = {}
for f in fj:
    d = pd.read_csv(os.path.join(TCGA_STAR, f['file_id'], f['file_name']), sep='\t', comment='#', usecols=['gene_id', 'gene_name', 'tpm_unstranded'])
    d = d[~d.gene_id.str.startswith('N_')]
    mats[f['cases'][0]['submitter_id']] = d.groupby('gene_name').tpm_unstranded.sum()
TT = pd.DataFrame(mats)
def cbio(kind):
    j = requests.get(f'https://www.cbioportal.org/api/studies/laml_tcga_pub/clinical-data?clinicalDataType={kind}&projection=SUMMARY', timeout=60).json()
    return pd.DataFrame(j).pivot(index='patientId', columns='clinicalAttributeId', values='value')
clin = cbio('PATIENT').join(cbio('SAMPLE'), rsuffix='_s')
clin.to_csv(os.path.join(INTERIM, 'tcga_laml_pub_clinical_cbioportal.csv'))
tcl = pd.DataFrame(index=TT.columns)
for k in ['OS_MONTHS', 'OS_STATUS', 'AGE', 'RISK_MOLECULAR', 'RISK_CYTO', 'TRANSPLANT_TYPE', 'BM_BLAST_PERCENTAGE']:
    tcl[k] = clin.reindex(tcl.index)[k]
tcl.index.name = 'id'; tcl.to_csv(os.path.join(INTERIM, 'g04_tcga_clin.csv'))
np.log2(TT + 1).astype('float32').to_csv(os.path.join(INTERIM, 'g04_tcga_log2tpm.csv.gz'))
tc = TT.reindex(GENES).T; tc['bm_blast'] = pd.to_numeric(tcl.BM_BLAST_PERCENTAGE, errors='coerce')

# ---- Beat AML
glen, gname = gene_lengths()
cnt = pd.read_csv(os.path.join(OPEN, 'beataml_waves1to4_counts_dbgap.txt'), sep='\t').set_index('stable_id')
sym = cnt.pop('display_label'); cnt = cnt.drop(columns=['description', 'biotype'])
keep = glen.reindex(cnt.index).notna()
print('Beat AML genes', len(cnt), 'dropped (no GENCODE v36 length)', int((~keep).sum()))
T = tpm(cnt, glen); symk = sym[T.index]
cl = pd.read_excel(os.path.join(OPEN, 'beataml_wv1to4_clinical.xlsx'), sheet_name='summary')
c = cl[cl.dbgap_rnaseq_sample.isin(T.columns) & (cl.diseaseStageAtSpecimenCollection == 'Initial Diagnosis')
       & (pd.to_numeric(cl.ageAtDiagnosis, errors='coerce') >= 18)
       & cl.dxAtInclusion.str.contains('ACUTE MYELOID', case=False, na=False)].copy()
c['bm'] = c.specimenType.str.contains('Marrow', case=False, na=False)
c['order'] = pd.to_numeric(c.timeOfSampleCollectionRelativeToInclusion, errors='coerce')
c = c.sort_values(['dbgap_subject_id', 'bm', 'order'], ascending=[True, False, True]).drop_duplicates('dbgap_subject_id')
c['blast'] = np.where(c.bm, pd.to_numeric(c['%.Blasts.in.BM'], errors='coerce'), pd.to_numeric(c['%.Blasts.in.PB'], errors='coerce'))
print('Beat AML eligible patients', len(c), ' blast>=50:', int((c.blast >= 50).sum()))
B = pd.DataFrame({g: T.loc[symk[symk == g].index].sum(0) for g in GENES}).loc[c.dbgap_rnaseq_sample.values]
bl = c.set_index('dbgap_rnaseq_sample').blast
BA = T[c.dbgap_rnaseq_sample.values].groupby(symk.values).sum()
np.log2(BA + 1).astype('float32').to_csv(os.path.join(INTERIM, 'g04_beat_log2tpm.csv.gz'))
c.set_index('dbgap_rnaseq_sample')[['dbgap_subject_id', 'dbgap_dnaseq_sample', 'ageAtDiagnosis', 'ELN2017', 'overallSurvival',
    'vitalStatus', 'blast', 'specimenType']].to_csv(os.path.join(INTERIM, 'g04_beat_clin.csv'))
elig = c.loc[c.blast >= 50, 'dbgap_rnaseq_sample'].values
pd.DataFrame({'symbol': symk.values, 'median_tpm': T.loc[symk.index, elig].median(1).values}).groupby('symbol').median_tpm.max() \
  .to_csv(os.path.join(INTERIM, 'beataml_blast_median_tpm.csv'))
# per-sample bulk values for figures
pd.concat([tc.assign(cohort='TCGA-LAML', blast=tc.bm_blast).drop(columns='bm_blast'),
           B.assign(cohort='Beat AML', blast=bl.reindex(B.index).values)]).to_csv(os.path.join(TABLES, 'gate01_bulk_per_sample_tpm.csv'))

def summ(df, label):
    return [dict(cohort=label, gene=g, n=len(df), median_TPM=df[g].median(), q75=df[g].quantile(.75), q95=df[g].quantile(.95),
                 max=df[g].max(), frac_gt1=(df[g] > 1).mean(), frac_gt0=(df[g] > 0).mean()) for g in GENES]
res = pd.DataFrame(summ(tc[tc.bm_blast >= 50], 'TCGA blast>=50 [DECISIONAL]') + summ(tc, 'TCGA all') +
                   summ(B[bl.reindex(B.index) >= 50], 'BeatAML blast>=50 [DECISIONAL]') + summ(B, 'BeatAML all eligible'))
res.to_csv(os.path.join(TABLES, 'gate01_bulk_detection.csv'), index=False)
print(res[res.gene.isin(['MRGPRD', 'ALB', 'CD34'])].round(3).to_string(index=False))
