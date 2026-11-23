"""Shared helpers: GENCODE v36 union-exon gene lengths and TPM."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gzip, re, os
import pandas as pd
from config import GTF, INTERIM

def gene_lengths():
    """Union-exon length per version-stripped Ensembl gene ID, and gene_name map (cached)."""
    cache = os.path.join(INTERIM, 'gencode_v36_lengths.csv')
    if os.path.exists(cache):
        d = pd.read_csv(cache, index_col=0); return d['length'], d['gene_name']
    ex, name = [], {}
    with gzip.open(GTF, 'rt') as fh:
        for l in fh:
            if l[0] == '#': continue
            a = l.split('\t', 9)
            if a[2] != 'exon': continue
            gid = re.search('gene_id "([^"]+)"', a[8]).group(1).split('.')[0]
            name[gid] = re.search('gene_name "([^"]+)"', a[8]).group(1)
            ex.append((gid, int(a[3]), int(a[4])))
    ex = pd.DataFrame(ex, columns=['g', 's', 'e']).sort_values(['g', 's'])
    def ulen(df):
        tot = 0; cs = ce = None
        for s, e in zip(df.s.values, df.e.values):
            if cs is None: cs, ce = s, e
            elif s <= ce + 1: ce = max(ce, e)
            else: tot += ce - cs + 1; cs, ce = s, e
        return tot + ce - cs + 1
    L = ex.groupby('g')[['s', 'e']].apply(ulen)
    d = pd.DataFrame({'length': L, 'gene_name': pd.Series(name).reindex(L.index)})
    d.to_csv(cache); return d['length'], d['gene_name']

def tpm(counts, lengths):
    """counts: genes x samples; lengths indexed like counts.index (bp). Genes without length are dropped."""
    L = lengths.reindex(counts.index); k = L.notna()
    rpk = counts[k].div(L[k] / 1e3, axis=0)
    return rpk.div(rpk.sum(0), axis=1) * 1e6

def expected_reads_at_1tpm(counts, lengths, gene_len_bp):
    """Per-sample reads expected for a gene of gene_len_bp at 1 TPM: 1e-6 * sum(RPK) * len_kb."""
    L = lengths.reindex(counts.index); k = L.notna()
    return counts[k].div(L[k] / 1e3, axis=0).sum(0) / 1e6 * gene_len_bp / 1e3
