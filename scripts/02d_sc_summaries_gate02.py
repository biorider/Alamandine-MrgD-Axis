"""Gate 0.1 single-cell summaries + calibration, and gate 0.2 co-detection, for all four scRNA datasets."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from config import INTERIM, TABLES
import sc_common as S
ref = pd.read_csv(os.path.join(INTERIM, 'beataml_blast_median_tpm.csv'), index_col=0).iloc[:, 0]
cells = {'GSE116256': pd.read_parquet(os.path.join(INTERIM, 'sc_cells_GSE116256.parquet')),
         'GSE241989': pd.read_parquet(os.path.join(INTERIM, 'sc_cells_GSE241989.parquet')),
         'GSE207104': pd.read_parquet(os.path.join(INTERIM, 'sc_cells_GSE207104.parquet')),
         'GSE253355': pd.read_csv(os.path.join(INTERIM, 'sc_cells_GSE253355.csv'))}
dets = {'GSE116256': pd.read_parquet(os.path.join(INTERIM, 'sc_det_GSE116256.parquet')),
        'GSE241989': pd.read_parquet(os.path.join(INTERIM, 'sc_det_GSE241989.parquet')),
        'GSE207104': pd.read_parquet(os.path.join(INTERIM, 'sc_det_GSE207104.parquet')),
        'GSE253355': pd.read_csv(os.path.join(INTERIM, 'sc_det_GSE253355.csv'), index_col=0)}
summ, cal = [], []
vg = cells['GSE116256']
summ.append(S.summarise(vg, 'GSE116256 all timepoints'))
summ.append(S.summarise(vg[vg.timepoint.isin(['D0', '']) | ~vg.patient.str.startswith('AML')], 'GSE116256 D0+normal (sensitivity)'))
c = S.calibration(dets['GSE116256'], ref, 'blast (malignant)'); c['dataset'] = 'GSE116256'; cal.append(c)
for gse in ['GSE241989', 'GSE207104']:
    ct = cells[gse].copy(); normal = ct.patient.str.startswith('HD')
    blastlike = ct.compartment.isin(['HSPC/blast-like', 'promyelocyte/GMP-like'])
    ct['compartment'] = np.where(blastlike & ~normal, 'blast (HSPC/GMP-like, AML)', np.where(blastlike & normal, 'HSPC (normal donor)', ct.compartment))
    cells[gse] = ct; summ.append(S.summarise(ct, gse))
    for k in ['HSPC/blast-like', 'promyelocyte/GMP-like']:
        c = S.calibration(dets[gse], ref, k); c['dataset'] = gse; cal.append(c)
summ.append(S.summarise(cells['GSE253355'], 'GSE253355 (normal marrow)'))
c = S.calibration(dets['GSE253355'], ref, 'HSPC'); c['dataset'] = 'GSE253355'; c['note'] = 'EXPLORATORY (not pre-specified)'; cal.append(c)
summ = pd.concat(summ); cal = pd.concat(cal)
summ.to_csv(os.path.join(TABLES, 'gate01_sc_detection_summary.csv'), index=False)
cal.to_csv(os.path.join(TABLES, 'gate01_sc_calibration.csv'), index=False)
rows = []
for k, d in cells.items():
    for comp, s in d.groupby('compartment'):
        rows.append(dict(dataset=k, compartment=comp, n_cells=len(s), MRGPRD_pos=int((s.MRGPRD > 0).sum()), MAS1_pos=int((s.MAS1 > 0).sum()),
                         co_detected=int(((s.MRGPRD > 0) & (s.MAS1 > 0)).sum()), n_patients_MRGPRDpos=s.loc[s.MRGPRD > 0, 'patient'].nunique(),
                         expected_co_under_independence=float((s.MRGPRD > 0).sum() * (s.MAS1 > 0).mean())))
g2 = pd.DataFrame(rows); g2['testable'] = g2.MRGPRD_pos >= 20
g2['verdict'] = np.where(g2.testable, 'test', 'UNTESTABLE (pre-committed power rule)')
g2.to_csv(os.path.join(TABLES, 'gate02_codetection.csv'), index=False)
print('cells', g2.n_cells.sum(), 'MRGPRD+', g2.MRGPRD_pos.sum(), 'MAS1+', g2.MAS1_pos.sum(), 'co-detected', g2.co_detected.sum())
