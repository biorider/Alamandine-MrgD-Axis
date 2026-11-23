"""GSE116256 (van Galen 2019, Seq-Well): per-cell target/control UMIs from RAW dem/anno tables; published labels."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, re, sys, tarfile, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(__file__))
from config import NEW_AML, INTERIM
import sc_common as S
tar = tarfile.open(os.path.join(NEW_AML, 'GSE116256_RAW.tar'))
names = [n for n in tar.getnames() if n.endswith('.dem.txt.gz')]
annomap = {re.sub(r'^GSM\d+_|\.anno\.txt\.gz$', '', n): n for n in tar.getnames() if n.endswith('.anno.txt.gz')}
def comp(r):
    t, p = str(r.CellType), str(r.PredictionRefined)
    if p == 'malignant': return 'blast (malignant)'
    if t in ('HSC', 'Prog'): return 'HSPC'
    if t == 'GMP': return 'promyelocyte/GMP-like'
    if t in ('ProMono', 'Mono', 'cDC'): return 'mono/macro'
    if t in ('T', 'CTL', 'NK'): return 'T/NK'
    if t in ('ProB', 'B', 'Plasma'): return 'B'
    if t in ('earlyEry', 'lateEry'): return 'erythroid'
    return 'other/unlabelled'
cts, dets = [], []
for n in names:
    samp = re.sub(r'^GSM\d+_|\.dem\.txt\.gz$', '', n)
    dem = pd.read_csv(tar.extractfile(n), sep='\t', index_col=0, compression='gzip')
    an = pd.read_csv(tar.extractfile(annomap[samp]), sep='\t', compression='gzip').set_index('Cell').reindex(dem.columns)
    ct = an[['CellType', 'PredictionRefined']].copy(); ct['compartment'] = ct.apply(comp, axis=1)
    ct['sample'] = samp; ct['patient'] = samp.split('-')[0]; ct['timepoint'] = samp.split('-')[1] if '-' in samp else ''
    X = sp.csr_matrix(dem.values.T)
    cts.append(S.cell_table(X, dem.index, ct))
    d = S.compartment_gene_detection(X, dem.index, ct.compartment.values); n_ = ct.compartment.value_counts()
    dets.append(d.mul(n_.reindex(d.columns).values, axis=1))
vg = pd.concat(cts); vg.to_parquet(os.path.join(INTERIM, 'sc_cells_GSE116256.parquet'))
num = pd.concat(dets, axis=1).T.groupby(level=0).sum().T.fillna(0)
den = vg.compartment.value_counts()
(num / den.reindex(num.columns).values).to_parquet(os.path.join(INTERIM, 'sc_det_GSE116256.parquet'))
print(vg.shape, vg.compartment.value_counts().to_dict())
