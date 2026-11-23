"""02b: Gate 0.1/0.2 processing for GSE241989 and GSE207104 (non-T-enriched only). Frozen protocol v1 Section F.3.
Cell QC (not fixed by protocol; most lenient common choice, declared): cells with >=200 detected genes."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import sys, tarfile, re, io, gzip, os, tempfile
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import NEW_AML, INTERIM
import sc_common as S
# compatibility shim: scanpy 1.11.5 calls .getnnz on scipy csr_array, removed in recent scipy
if not hasattr(sp.csr_array, 'getnnz'):
    sp.csr_array.getnnz = lambda self, axis=None: (np.diff(self.indptr) if axis == 1 else self.nnz)

def load_241989():
    tar = tarfile.open(os.path.join(NEW_AML, 'GSE241989_RAW.tar'))
    mats = [n for n in tar.getnames() if n.endswith('matrix.mtx.gz')]
    ads = []
    for m in mats:
        pre = m.replace('-matrix.mtx.gz', '')
        samp = re.sub(r'^GSM\d+_', '', pre)
        X = sio.mmread(gzip.open(tar.extractfile(m))).T.tocsr()
        f = pd.read_csv(tar.extractfile(pre + '-features.tsv.gz'), sep='\t', header=None, compression='gzip')
        b = pd.read_csv(tar.extractfile(pre + '-barcodes.tsv.gz'), sep='\t', header=None, compression='gzip')
        if f.shape[1] >= 3:
            keep = (f[2] == 'Gene Expression').values
            X = X[:, keep]; f = f[keep]
        a = ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=samp + '_' + b[0].astype(str)),
                       var=pd.DataFrame(index=f[1].astype(str).values))
        a.var_names_make_unique(); a.obs['sample'] = samp; a.obs['patient'] = samp
        ads.append(a)
    return ad.concat(ads, join='outer', fill_value=0)

def load_207104():
    tar = tarfile.open(os.path.join(NEW_AML, 'GSE207104_RAW.tar'))
    h5s = [n for n in tar.getnames() if n.endswith('.h5') and 'NonT' in n]
    ads = []
    for n in h5s:
        samp = re.sub(r'^GSM\d+_|_GEX\.h5$', '', n)
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False, dir=INTERIM) as tf:
            tf.write(tar.extractfile(n).read()); path = tf.name
        a = sc.read_10x_h5(path); os.remove(path)
        a.var_names_make_unique()
        a.obs_names = samp + '_' + a.obs_names
        a.obs['sample'] = samp; a.obs['patient'] = samp.split('_')[0]
        ads.append(a)
    return ad.concat(ads, join='outer', fill_value=0)

def label(a):
    raw = a.X.copy()
    sc.pp.filter_cells(a, min_genes=200)
    raw = a.X.copy()
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    excl = set(S.TARGETS)
    sc.pp.highly_variable_genes(b, n_top_genes=2000)
    b.var.loc[b.var_names.isin(excl), 'highly_variable'] = False
    sc.pp.pca(b, n_comps=30, mask_var='highly_variable', random_state=0)
    # kNN (k=15) graph on 30 PCs + Leiden res 1.0 seed 0, built directly (scanpy/umap/numba incompatibility in env).
    # Declared implementation deviation: unweighted symmetrised kNN graph instead of UMAP-weighted connectivities.
    from sklearn.neighbors import NearestNeighbors
    import igraph as ig, leidenalg
    pcs = b.obsm['X_pca']
    ind = NearestNeighbors(n_neighbors=16).fit(pcs).kneighbors(return_distance=False)[:, 1:]
    src = np.repeat(np.arange(pcs.shape[0]), ind.shape[1]); dst = ind.ravel()
    g = ig.Graph(n=pcs.shape[0], edges=list(zip(src.tolist(), dst.tolist())), directed=False).simplify()
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=1.0, seed=0)
    b.obs['leiden'] = pd.Categorical([str(x) for x in part.membership])
    Z = pd.DataFrame(b.X.toarray() if sp.issparse(b.X) else b.X, columns=b.var_names) if b.n_obs < 0 else None
    # cluster mean of z-scored marker expression
    scores = {}
    for k, genes in S.MARKERS.items():
        g = [x for x in genes if x in b.var_names]
        sc.tl.score_genes(b, g, score_name='s_' + k, random_state=0, ctrl_size=50) if False else None
        Xg = b[:, g].X
        Xg = Xg.toarray() if sp.issparse(Xg) else Xg
        Xg = (Xg - Xg.mean(0)) / (Xg.std(0) + 1e-9)
        scores[k] = Xg.mean(1)
    sdf = pd.DataFrame(scores, index=b.obs_names)
    sdf['cluster'] = b.obs['leiden'].values
    cm = sdf.groupby('cluster').mean()
    assign = cm.idxmax(1)
    a.obs['cluster'] = b.obs['leiden'].values
    a.obs['compartment'] = a.obs['cluster'].map(assign).values
    return a, raw, cm

if __name__ == '__main__':
    for gse, loader in [('GSE241989', load_241989), ('GSE207104', load_207104)]:
        a = loader()
        S.assert_features(a.var_names)
        a, raw, cm = label(a)
        meta = a.obs[['sample', 'patient', 'cluster', 'compartment']].copy()
        ct = S.cell_table(raw, a.var_names, meta)
        ct.to_parquet(os.path.join(INTERIM, f'sc_cells_{gse}.parquet'))
        det = S.compartment_gene_detection(raw, a.var_names, meta.compartment.values)
        det.to_parquet(os.path.join(INTERIM, f'sc_det_{gse}.parquet'))
        cm.to_csv(os.path.join(INTERIM, f'sc_cluster_marker_means_{gse}.csv'))
        print(gse, a.shape, meta.compartment.value_counts().to_dict(), flush=True)
