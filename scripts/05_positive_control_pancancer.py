"""Step 1: MRGPRD-specific positive control. Pan-cancer MRGPRD (cBioPortal PanCancer Atlas RSEM) -> 13 top tumours ->
their GDC STAR files (same GENCODE v36 TPM pipeline as TCGA-LAML). Also LAML vs other cohorts on the RSEM scale with ALB."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import io, os, json, requests, pandas as pd
from config import TABLES
API = 'https://www.cbioportal.org/api'
studies = [s['studyId'] for s in requests.get(f'{API}/studies', timeout=60).json() if s['studyId'].endswith('_tcga_pan_can_atlas_2018')]
ent = {g['hugoGeneSymbol']: g['entrezGeneId'] for g in requests.post(f'{API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', json=['MRGPRD', 'ALB'], timeout=30).json()}
rows = []
for st in studies:
    d = requests.post(f'{API}/molecular-profiles/{st}_rna_seq_v2_mrna/molecular-data/fetch',
                      json={'entrezGeneIds': list(ent.values()), 'sampleListId': f'{st}_all'}, timeout=120).json()
    inv = {v: k for k, v in ent.items()}
    rows += [dict(study=st.replace('_tcga_pan_can_atlas_2018', ''), sample=x['sampleId'], gene=inv[x['entrezGeneId']], rsem=x['value']) for x in d]
pc = pd.DataFrame(rows); pc.to_csv(os.path.join(TABLES, 'step1_pancan_rsem_MRGPRD_ALB.csv'), index=False)
m = pc[pc.gene == 'MRGPRD']
top = m.sort_values('rsem', ascending=False).head(13)
filt = {"op": "and", "content": [{"op": "in", "content": {"field": "cases.submitter_id", "value": [s[:12] for s in top['sample']]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}]}
hits = requests.get('https://api.gdc.cancer.gov/files', params={'filters': json.dumps(filt), 'fields': 'file_id,cases.samples.submitter_id',
                    'size': '100', 'format': 'JSON'}, timeout=60).json()['data']['hits']
out = []
for h in hits:
    samp = h['cases'][0]['samples'][0]['submitter_id']
    if samp[:15] not in set(top['sample'].str[:15]): continue
    d = pd.read_csv(io.StringIO(requests.get(f"https://api.gdc.cancer.gov/data/{h['file_id']}", timeout=120).text), sep='\t', comment='#')
    d = d[~d.gene_id.str.startswith('N_')].set_index('gene_name')
    out.append(dict(sample=samp, MRGPRD_TPM=d.tpm_unstranded.get('MRGPRD'), MRGPRD_counts=d.unstranded.get('MRGPRD'),
                    ALB_TPM=d.tpm_unstranded.get('ALB'), CD34_TPM=d.tpm_unstranded.get('CD34')))
pos = pd.DataFrame(out).drop_duplicates('sample'); pos['s15'] = pos['sample'].str[:15]
pos = pos.merge(top.assign(s15=top['sample'].str[:15])[['s15', 'study', 'rsem']].rename(columns={'rsem': 'cbio_RSEM'}), on='s15').drop(columns='s15')
pos.to_csv(os.path.join(TABLES, 'step1_positive_control_gdc_star.csv'), index=False); print(pos.round(2).to_string(index=False))
