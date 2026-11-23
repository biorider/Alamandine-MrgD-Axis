"""Step 8: MRGPRD/MAS1 in control PMA-differentiated THP-1 macrophages (GSE348126 empty vector A1-A3; GSE323383 scrambled;
GSE318527 CST3+/+, author TPM - weak macrophage markers, supportive only). Undifferentiated THP-1: DepMap (script 03)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, re, pandas as pd
from config import OPEN, TABLES
from common import gene_lengths, tpm, expected_reads_at_1tpm
glen, gname = gene_lengths(); Lm = glen['ENSG00000172938']; g8 = ['MRGPRD', 'MAS1', 'CD14', 'ITGAM', 'CD68', 'ALB']
f = lambda n: os.path.join(OPEN, 'thp1', n); rows = []
def from_counts(d, label):
    T = tpm(d, glen); Ts = T.groupby(gname.reindex(T.index).values).sum()
    rows.append(dict(dataset=label, **Ts.reindex(g8).mean(1).round(3).to_dict(), MRGPRD_reads_total=int(d.loc['ENSG00000172938'].sum()),
                     expected_reads_at_1TPM=round(expected_reads_at_1tpm(d, glen, Lm).sum(), 1)))
d = pd.read_csv(f('GSE348126_Gene_COUNT.txt.gz'), sep='\t').set_index('gene_id'); d.index = d.index.str.split('.').str[0]
from_counts(d[['A1', 'A2', 'A3']], 'GSE348126 (empty vector, n=3)')
d = pd.read_csv(f('GSE323383_THP1_bulkRNAseq_raw_counts.tsv.gz'), sep='\t').set_index('gene_id'); d.index = d.index.str.split('.').str[0]
cc = [c for c in d.columns if 'scr' in c]; from_counts(d[cc].fillna(0), f'GSE323383 (scrambled, n={len(cc)})')
d = pd.read_csv(f('GSE318527_gene_expression.txt.gz'), sep='\t').set_index('gene_symbol')
tc = [c for c in d.columns if c.startswith('tpm_control')]; rc = [c for c in d.columns if re.fullmatch(r'read_count_control_\d', c)]
rows.append(dict(dataset=f'GSE318527 (CST3+/+, n={len(tc)}; author TPM)', **d[tc].groupby(level=0).sum().reindex(g8).mean(1).round(3).to_dict(),
                 MRGPRD_reads_total=int(d.loc['MRGPRD', rc].sum()), expected_reads_at_1TPM=None))
r = pd.DataFrame(rows); r.to_csv(os.path.join(TABLES, 'step8_pma_thp1_expression.csv'), index=False); print(r.to_string(index=False))
