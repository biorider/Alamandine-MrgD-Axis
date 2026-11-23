"""Step 5: MRGPRD in AML-associated stroma. GSE107490 (cultured MSC: AML/healthy/MDS) and GSE212615 (AML FFPE core biopsies).
TPM with GENCODE v36 lengths; expected reads at 1 and 0.1 TPM per sample; Poisson test of observed counts."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, re, tarfile, pandas as pd
from scipy.stats import poisson
from config import OPEN, TABLES
from common import gene_lengths, tpm, expected_reads_at_1tpm
glen, gname = gene_lengths()
symlen = pd.DataFrame({'sym': gname.values, 'len': glen.values}).groupby('sym').len.max()
genes = ['MRGPRD', 'MAS1', 'CXCL12', 'LEPR', 'COL1A1', 'PTPRC', 'ALB']; Lm = glen['ENSG00000172938']
A = pd.read_csv(os.path.join(OPEN, 'stroma', 'GSE107490_all_count.txt.gz'), sep='\t').set_index('ORF')
tpm(A, symlen).reindex(genes).T.round(3).to_csv(os.path.join(TABLES, 'step5_GSE107490_tpm.csv'))
cols = {}
with tarfile.open(os.path.join(OPEN, 'stroma', 'GSE212615_RAW.tar')) as tf:
    for n in tf.getnames():
        cols[re.sub(r'^GSM\d+_|\.txt\.gz$', '', n)] = pd.read_csv(tf.extractfile(n), sep=r'\s+', compression='gzip', header=0, index_col=0).iloc[:, 0]
C = pd.DataFrame(cols); C.index = C.index.str.replace('"', '').str.split('.').str[0]
TC = tpm(C, glen); TC.groupby(gname.reindex(TC.index).values).sum().reindex(genes).T.round(3).to_csv(os.path.join(TABLES, 'step5_GSE212615_tpm.csv'))
eA = expected_reads_at_1tpm(A, symlen, Lm); eC = expected_reads_at_1tpm(C, glen, Lm)
aml = ['AML1', 'AML2', 'AML3']; rows = []
for nm, n, e, obs in [('GSE107490 AML MSC', 3, eA[aml], A.loc['MRGPRD', aml].sum()), ('GSE107490 all MSC', 12, eA, A.loc['MRGPRD'].sum()),
                      ('GSE212615 AML core biopsies', 22, eC, C.loc['ENSG00000172938'].sum())]:
    rows.append(dict(dataset=nm, n=n, obs=int(obs), exp_1TPM=round(e.sum(), 1), p_1TPM=poisson.cdf(obs, e.sum()),
                     exp_01TPM=round(e.sum() * .1, 1), p_01TPM=poisson.cdf(obs, e.sum() * .1),
                     median_library_M=round((A if 'GSE107490' in nm else C).sum().median() / 1e6, 1)))
r = pd.DataFrame(rows); r.to_csv(os.path.join(TABLES, 'step5_stroma_sensitivity_bound.csv'), index=False); print(r.to_string(index=False))
