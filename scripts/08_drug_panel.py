"""Step 4 (pre-registered secondary, kickoff 5.4): Spearman of RAS score (GSVA, from 04) with Beat AML ex vivo AUC for every
inhibitor tested in >=100 eligible patients with a score; BH across the full panel. Then an EXPLORATORY monocytic adjustment."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests
from config import OPEN, INTERIM, TABLES
dr = pd.read_csv(os.path.join(OPEN, 'beataml_probit_curve_fits_v4_dbgap.txt'), sep='\t')
ids = set(pd.read_csv(os.path.join(INTERIM, 'g04_beat_clin.csv')).dbgap_rnaseq_sample)
d = dr[dr.dbgap_rnaseq_sample.isin(ids) & (dr.paper_inclusion.astype(str).str.upper() == 'TRUE')]
d = d.groupby(['dbgap_rnaseq_sample', 'inhibitor'], as_index=False).auc.mean()
sc = pd.read_csv(os.path.join(TABLES, 'gate04_scores.csv')); sc = sc[sc.cohort == 'BeatAML'].set_index('id')
out = []
for inh, g in d.groupby('inhibitor'):
    g = g[g.dbgap_rnaseq_sample.isin(sc.index)]
    if g.dbgap_rnaseq_sample.nunique() < 100: continue
    s = sc.loc[g.dbgap_rnaseq_sample]
    for nm in ['RAS8', 'RAS7']:
        rho, p = spearmanr(s[nm].values, g.auc.values); out.append(dict(inhibitor=inh, score=nm, n=len(g), rho=rho, p=p))
res = pd.DataFrame(out)
for nm in ['RAS8', 'RAS7']:
    m = res.score == nm; res.loc[m, 'q_BH'] = multipletests(res.loc[m, 'p'], method='fdr_bh')[1]
res.to_csv(os.path.join(TABLES, 'step4_drug_panel_spearman.csv'), index=False)
r8 = res[res.score == 'RAS8']
print('panel', len(r8), 'q<0.05 RAS8', int((r8.q_BH < .05).sum()), 'RAS7', int((res[res.score == 'RAS7'].q_BH < .05).sum()))
# EXPLORATORY: monocytic score = mean z of CD14, LYZ, CSF1R, FCGR3A; rank-residualised partial Spearman
X = pd.read_csv(os.path.join(INTERIM, 'g04_beat_log2tpm.csv.gz'), index_col=0)
Z = X.loc[['CD14', 'LYZ', 'CSF1R', 'FCGR3A']].T; mono = ((Z - Z.mean()) / Z.std()).mean(1)
s8 = sc.RAS8; common = s8.index.intersection(mono.index)
def resid(y, x):
    ry, rx = rankdata(y), rankdata(x); b = np.polyfit(rx, ry, 1); return ry - np.polyval(b, rx)
ex = []
for inh, g in d.groupby('inhibitor'):
    g = g[g.dbgap_rnaseq_sample.isin(common)]
    if g.dbgap_rnaseq_sample.nunique() < 100: continue
    s, a, mo = s8[g.dbgap_rnaseq_sample].values, g.auc.values, mono[g.dbgap_rnaseq_sample].values
    pr, pp = spearmanr(resid(s, mo), resid(a, mo))
    ex.append(dict(inhibitor=inh, n=len(g), rho_mono_auc=spearmanr(mo, a)[0], partial_rho_RAS8=pr, partial_p=pp))
ex = pd.DataFrame(ex); ex['partial_q'] = multipletests(ex.partial_p, method='fdr_bh')[1]
ex = ex.merge(r8[['inhibitor', 'rho', 'q_BH']].rename(columns={'rho': 'raw_rho_RAS8', 'q_BH': 'raw_q'}), on='inhibitor')
ex.to_csv(os.path.join(TABLES, 'step4_drug_panel_EXPLORATORY_monocytic_adjusted.csv'), index=False)
print('rho(RAS8, mono) %.3f; adjusted q<0.05: %d; across-drug r %.3f' % (spearmanr(s8[common], mono[common])[0], (ex.partial_q < .05).sum(),
      np.corrcoef(ex.raw_rho_RAS8, ex.rho_mono_auc)[0, 1]))
