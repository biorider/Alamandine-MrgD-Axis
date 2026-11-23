"""Download all open inputs and verify/record SHA-256 (data/manifests/open_data_manifest.csv)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, hashlib, time, requests, pandas as pd
from config import OPEN, ROOT
GH = 'https://github.com/biodev/beataml2.0_data/raw/main/'
GEO = lambda acc, f: f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{f}'
FILES = [(GH + f, f) for f in ['beataml_waves1to4_counts_dbgap.txt', 'beataml_waves1to4_norm_exp_dbgap.txt',
         'beataml_probit_curve_fits_v4_dbgap.txt', 'beataml_wv1to4_clinical.xlsx', 'beataml_waves1to4_sample_mapping.xlsx',
         'beataml_wes_wv1to4_mutations_dbgap.txt']] + [
    (GEO('GSE253355', f), f) for f in ['GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz', 'GSE253355_RAW.tar']] + [
    (GEO('GSE107490', 'GSE107490_all_count.txt.gz'), 'stroma/GSE107490_all_count.txt.gz'),
    (GEO('GSE212615', 'GSE212615_RAW.tar'), 'stroma/GSE212615_RAW.tar'),
    (GEO('GSE348126', 'GSE348126_Gene_COUNT.txt.gz'), 'thp1/GSE348126_Gene_COUNT.txt.gz'),
    (GEO('GSE323383', 'GSE323383_THP1_bulkRNAseq_raw_counts.tsv.gz'), 'thp1/GSE323383_THP1_bulkRNAseq_raw_counts.tsv.gz'),
    (GEO('GSE318527', 'GSE318527_gene_expression.txt.gz'), 'thp1/GSE318527_gene_expression.txt.gz'),
    ('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz', 'gencode.v36.annotation.gtf.gz')]
man_path = os.path.join(ROOT, 'data', 'manifests', 'open_data_manifest.csv')
ref = pd.read_csv(man_path).set_index('file') if os.path.exists(man_path) else None
rows = []
for url, rel in FILES:
    fn = os.path.join(OPEN, rel); os.makedirs(os.path.dirname(fn), exist_ok=True)
    if not os.path.exists(fn):
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(fn, 'wb') as fh:
                for ch in r.iter_content(1 << 22): fh.write(ch)
    h = hashlib.sha256()
    with open(fn, 'rb') as fh:
        for ch in iter(lambda: fh.read(1 << 22), b''): h.update(ch)
    sha = h.hexdigest(); ok = None
    if ref is not None and rel in ref.index: ok = (ref.loc[rel, 'sha256'] == sha)
    rows.append(dict(file=rel, url=url, bytes=os.path.getsize(fn), sha256=sha, matches_recorded=ok))
    print(rel, os.path.getsize(fn), sha[:12], '' if ok is None else ('OK' if ok else 'HASH MISMATCH'))
if ref is None: pd.DataFrame(rows).to_csv(man_path, index=False)
