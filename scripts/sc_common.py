"""Shared single-cell routines for Phase 0 gates 0.1/0.2 (frozen protocol v1, sha256 466cf2ff...6844).
No gene filtering is applied anywhere; MRGPRD/MAS1 presence in the feature list is asserted."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import numpy as np, pandas as pd, scipy.sparse as sp

TARGETS = ['MRGPRD', 'MAS1']
CONTROLS = ['CD34', 'MPO', 'CXCL12', 'ALB']
TRACK = TARGETS + CONTROLS

MARKERS = {  # frozen v1 Section F.3; MRGPRD/MAS1 deliberately absent
    'HSPC/blast-like': ['CD34', 'PROM1', 'KIT', 'SPINK2'],
    'promyelocyte/GMP-like': ['MPO', 'ELANE', 'AZU1', 'PRTN3'],
    'mono/macro': ['CD14', 'LYZ', 'CD68', 'CSF1R', 'FCGR3A'],
    'T/NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'],
    'B': ['MS4A1', 'CD79A'],
    'erythroid': ['HBB', 'GYPA', 'CA1'],
    'stroma': ['CXCL12', 'LEPR', 'COL1A1', 'VCAM1'],
}


def assert_features(genes):
    missing = [g for g in TARGETS if g not in set(genes)]
    if missing:
        raise AssertionError(f'NOT QUANTIFIABLE: {missing} absent from feature list')


def cell_table(X, genes, meta):
    """X: cells x genes sparse/dense raw UMI. Returns per-cell table of tracked genes + totals."""
    X = sp.csr_matrix(X)
    genes = pd.Index(genes)
    assert_features(genes)
    out = meta.copy()
    out['total_umi'] = np.asarray(X.sum(1)).ravel()
    out['n_genes'] = np.asarray((X > 0).sum(1)).ravel()
    for g in TRACK:
        idx = np.where(genes == g)[0]
        out[g] = np.asarray(X[:, idx].sum(1)).ravel() if len(idx) else np.nan
    return out


def compartment_gene_detection(X, genes, comp):
    """Fraction of cells with >=1 UMI per gene, per compartment (for calibration). Returns genes x compartments."""
    X = sp.csr_matrix(X)
    B = (X > 0).astype(np.float32)
    res = {}
    for k in pd.unique(comp):
        m = np.asarray(comp == k)
        if m.sum() == 0:
            continue
        res[k] = np.asarray(B[m].mean(0)).ravel()
    df = pd.DataFrame(res, index=pd.Index(genes))
    return df.groupby(level=0).max()


def patient_detection(ct, gene, min_cells=50):
    """Per patient x compartment detection fraction (>=1 UMI) for patients with >= min_cells."""
    g = ct.groupby(['compartment', 'patient'])
    d = g.agg(n_cells=(gene, 'size'), n_detected=(gene, lambda v: int((v > 0).sum())),
              umi_sum=(gene, 'sum')).reset_index()
    d['frac'] = d.n_detected / d.n_cells
    d['eligible'] = d.n_cells >= min_cells
    d['gene'] = gene
    return d


def summarise(ct, dataset, min_cells=50):
    rows = []
    for gene in TRACK:
        if ct[gene].isna().all():
            continue
        pdx = patient_detection(ct, gene, min_cells)
        for comp, sub in pdx.groupby('compartment'):
            e = sub[sub.eligible]
            allc = ct[ct.compartment == comp][gene]
            rows.append(dict(dataset=dataset, compartment=comp, gene=gene,
                             n_patients_eligible=len(e), n_cells=int(len(allc)),
                             median_patient_frac=float(e.frac.median()) if len(e) else np.nan,
                             max_patient_frac=float(e.frac.max()) if len(e) else np.nan,
                             pooled_frac=float((allc > 0).mean()), cells_detected=int((allc > 0).sum()),
                             zero_fraction=float((allc == 0).mean()),
                             mean_umi_per_cell=float(allc.mean()),
                             detected=bool(len(e) and e.frac.median() > 0.05)))
    return pd.DataFrame(rows)


def calibration(det, bulk_ref, comp, bins=((0.5, 1), (1, 2), (2, 5), (5, 10))):
    """Median per-cell detection fraction among genes in each bulk-TPM bin, for one compartment."""
    j = det[[comp]].join(bulk_ref.rename('tpm'), how='inner')
    out = []
    for lo, hi in bins:
        s = j[(j.tpm >= lo) & (j.tpm < hi)][comp]
        out.append(dict(compartment=comp, tpm_bin=f'[{lo},{hi})', n_genes=len(s),
                        median_detect=float(s.median()), q25=float(s.quantile(.25)), q75=float(s.quantile(.75))))
    return pd.DataFrame(out)
