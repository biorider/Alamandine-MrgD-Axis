"""Gate 0.3: DepMap CRISPR (Chronos) dependency for MRGPRD/MAS1 in AML vs non-heme lines; MYB/SPI1 positive, ALB negative controls.
Also writes MRGPRD/MAS1 expression for all lines (incl. THP-1). NB: read_csv(usecols=...) returns FILE order -> rename from headers."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from config import USER_DATA, TABLES
gs = ['MRGPRD', 'MAS1', 'MYB', 'SPI1', 'ALB']
f = lambda n: os.path.join(USER_DATA, n)
def load(fn, genes, index_col=0, extra=()):
    hdr = pd.read_csv(f(fn), nrows=0).columns; col = {c.split(' (')[0]: c for c in hdr}
    d = pd.read_csv(f(fn), usecols=list(extra) + [col[g] for g in genes] + ([hdr[0]] if index_col == 0 else []), index_col=index_col)
    return d.rename(columns={col[g]: g for g in genes})
ge = load('CRISPRGeneEffect.csv', gs)[gs]; gd = load('CRISPRGeneDependency.csv', gs)[gs]
md = pd.read_csv(f('Model.csv')).set_index('ModelID')
ex = load('OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv', gs + ['CD14', 'ITGAM'], index_col=None, extra=('ModelID', 'IsDefaultEntryForModel'))
ex = ex[ex.IsDefaultEntryForModel.astype(str).str.lower().isin(['yes', 'true', '1'])].drop_duplicates('ModelID').set_index('ModelID')
assert list(ge.columns) == gs
lin = md.reindex(ge.index)
aml = lin.OncotreePrimaryDisease.eq('Acute Myeloid Leukemia'); heme = lin.OncotreeLineage.isin(['Myeloid', 'Lymphoid']); other = ~heme
rows = []
for g in gs:
    a, o, h = ge.loc[aml, g].dropna(), ge.loc[other, g].dropna(), ge.loc[heme, g].dropna(); e = ex.reindex(ge.index[aml])[g]
    rows.append(dict(gene=g, n_AML=len(a), AML_median=a.median(), nonheme_median=o.median(), diff=a.median() - o.median(),
                     p_AML=mannwhitneyu(a, o, alternative='less').pvalue, heme_median=h.median(), p_heme=mannwhitneyu(h, o, alternative='less').pvalue,
                     frac_AML_dep_prob_gt05=(gd.loc[aml, g] > 0.5).mean(), AML_min=a.min(),
                     AML_expr_median_log2TPM1=e.median(), AML_frac_expr_TPMgt1=(e > 1).mean()))
r = pd.DataFrame(rows); t = r.gene.isin(['MRGPRD', 'MAS1'])
r.loc[t, 'q_AML_BH'] = multipletests(r.loc[t, 'p_AML'], method='fdr_bh')[1]
r['PASS_rule'] = (r.AML_median <= -0.5) & (r['diff'] <= -0.3) & (r.p_AML < 0.05)
r.loc[t, 'PASS_rule'] = r.loc[t, 'PASS_rule'] & (r.loc[t, 'q_AML_BH'] < 0.05)
r.to_csv(os.path.join(TABLES, 'gate03_depmap.csv'), index=False); print(r.round(3).to_string(index=False))
ge.join(lin[['StrippedCellLineName', 'OncotreeLineage', 'OncotreePrimaryDisease']]).to_csv(os.path.join(TABLES, 'gate03_depmap_gene_effect_per_line.csv'))
tp = (2 ** ex[gs + ['CD14', 'ITGAM']] - 1).round(4); tp = tp.join(md[['StrippedCellLineName', 'OncotreePrimaryDisease']])
tp.to_csv(os.path.join(TABLES, 'depmap_expression_TPM_per_line.csv'))
print(tp[tp.StrippedCellLineName.isin(['THP1', 'HS852T'])][['StrippedCellLineName', 'MRGPRD', 'MAS1']].to_string())
