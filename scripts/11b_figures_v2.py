"""11b: Publication figures v2 (Fig 1 panel titles below panels; Figs 2-4 titles removed, content in legends). (PNG 300 dpi + PDF) from results/tables only. Run after steps 01-10."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, ast
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from config import TABLES, FIGS

mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7,
                     'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'xtick.direction': 'out', 'ytick.direction': 'out',
                     'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.titlelocation': 'left'})
T = lambda f: pd.read_csv(os.path.join(TABLES, f))
FOCAL, GREY, POS, NEG = '#C0392B', '#8C8C8C', '#2E6DB4', '#B8B8B8'   # MRGPRD focal; controls blue (positive) / light grey (negative)
GCOL = {'MRGPRD': FOCAL, 'MAS1': '#E08E3C', 'ALB': NEG, 'CD34': POS, 'MPO': '#5DA5DA', 'MYB': POS, 'SPI1': '#5DA5DA'}
def letter(ax, s): ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')
def logy(ax, floor=0.01, top=1e5):
    ax.set_yscale('log'); ax.set_ylim(floor * 0.6, top)
    t = [x for x in [0.01, 0.1, 1, 10, 100, 1e3, 1e4, 1e5] if floor <= x <= top]
    ax.set_yticks(t); ax.set_yticklabels([('0' if x == floor else (f'{x:g}' if x < 1e3 else f'{x/1e3:g}k')) for x in t])
def save(fig, name):
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(FIGS, f'{name}.{ext}'), dpi=300, bbox_inches='tight')
    plt.close(fig)
rng = np.random.default_rng(0)

# ---------------- Figure 1: bulk ----------------
bulk = T('gate01_bulk_per_sample_tpm.csv'); bulk = bulk[bulk.blast >= 50]
pc = T('step1_positive_control_gdc_star.csv')
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.5), gridspec_kw={'width_ratios': [1, 1, 0.9], 'wspace': 0.35})
genes = ['ALB', 'MRGPRD', 'MAS1', 'CD34', 'MPO']; FL = 0.01
for ax, coh in zip(axs[:2], ['TCGA-LAML', 'Beat AML']):
    d = bulk[bulk.cohort == coh]
    for i, g in enumerate(genes):
        v = np.clip(d[g].values, FL, None); x = i + rng.uniform(-0.28, 0.28, len(v))
        ax.scatter(x, v, s=3, color=GCOL[g], alpha=0.55, lw=0, rasterized=True)
        ax.plot([i - 0.35, i + 0.35], [max(np.median(d[g]), FL)] * 2, color='k', lw=1.2)
    ax.axhline(1, color=GREY, lw=0.6, ls='--'); ax.text(4.45, 1.15, '1 TPM', fontsize=6, color=GREY, ha='right', va='bottom')
    ax.set_xticks(range(len(genes))); ax.set_xticklabels([f'$\\it{{{g}}}$' for g in genes], rotation=40, ha='right', rotation_mode='anchor')
    logy(ax, FL, 1e5); ax.set_xlabel(f'{coh} blasts (n = {len(d)})', fontsize=8, labelpad=4)
axs[0].set_ylabel('Expression (TPM; 0 plotted at floor)')
axs[1].set_yticklabels([])
ax = axs[2]; p = pc.sort_values('MRGPRD_TPM')
ax.scatter(range(len(p)), p.MRGPRD_TPM, s=14, color=FOCAL, zorder=3)
ax.vlines(range(len(p)), FL, p.MRGPRD_TPM, color=FOCAL, lw=0.6)
tc_med = bulk[bulk.cohort == 'TCGA-LAML'].MRGPRD.median()
ax.axhline(tc_med, color='k', lw=0.8); ax.text(0, tc_med * 1.3, f'TCGA-LAML median {tc_med:.2f}', fontsize=6, va='bottom', bbox=dict(facecolor='white', edgecolor='none', pad=0.5))
ax.axhline(1, color=GREY, lw=0.6, ls='--')
ax.set_xticks(range(len(p))); ax.set_xticklabels(p.study.str.upper(), rotation=90)
logy(ax, FL, 1e3); ax.set_xlabel('Same pipeline detects $\\it{MRGPRD}$\nin expressing tumours', fontsize=8, labelpad=4)
ax.set_ylabel('$\\it{MRGPRD}$ (TPM, GDC STAR)')
for a, l in zip(axs, 'abc'): letter(a, l)
save(fig, 'Figure1_bulk_expression_and_positive_control_v2')

# ---------------- Figure 2: single cell ----------------
sm = T('gate01_sc_detection_summary.csv'); cal = T('gate01_sc_calibration.csv'); g2 = T('gate02_codetection.csv')
fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0), gridspec_kw={'width_ratios': [1.35, 1], 'wspace': 0.75})
ax = axs[0]
ds = ['GSE116256 all timepoints', 'GSE241989', 'GSE207104', 'GSE253355 (normal marrow)']
lab = {'GSE116256 all timepoints': 'van Galen\n(Seq-Well, AML)', 'GSE241989': 'GSE241989\n(10x, AML)', 'GSE207104': 'GSE207104\n(10x, AML)', 'GSE253355 (normal marrow)': 'GSE253355\n(10x, normal\nmarrow)'}
keyc = {'GSE116256 all timepoints': ['blast (malignant)', 'HSPC', 'mono/macro', 'T/NK'],
        'GSE241989': ['blast (HSPC/GMP-like, AML)', 'mono/macro', 'T/NK'],
        'GSE207104': ['blast (HSPC/GMP-like, AML)', 'mono/macro', 'T/NK'],
        'GSE253355 (normal marrow)': ['HSPC', 'stroma', 'VSMC (perivascular)', 'mono/macro', 'T/NK']}
rows = [(d, c) for d in ds for c in keyc[d]]; gl = ['MRGPRD', 'MAS1', 'CD34', 'CXCL12', 'ALB']
M = np.full((len(rows), len(gl)), np.nan); N = []
for i, (d, c) in enumerate(rows):
    s = sm[(sm.dataset == d) & (sm.compartment == c)]
    for j, g in enumerate(gl):
        v = s[s.gene == g].pooled_frac
        if len(v): M[i, j] = v.iloc[0] * 100
    N.append(int(s.n_cells.iloc[0]) if len(s) else 0)
im = ax.imshow(np.log10(np.nan_to_num(M, nan=0) + 0.01), cmap='Greys', vmin=-2, vmax=2, aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]; txt = 'n.d.' if np.isnan(v) else ('0' if v == 0 else (f'{v:.2f}' if v < 1 else f'{v:.0f}'))
        ax.text(j, i, txt, ha='center', va='center', fontsize=5.5, color='white' if (not np.isnan(v) and v > 10) else 'black')
ax.set_xticks(range(len(gl))); ax.set_xticklabels([f'$\\it{{{g}}}$' for g in gl])
ax.set_yticks(range(len(rows))); ax.set_yticklabels([f'{c.replace(" (HSPC/GMP-like, AML)", "")}  (n={n/1000:.1f}k)' for (d, c), n in zip(rows, N)])
ax.tick_params(length=0)
y0 = 0
for d in ds:
    k = len(keyc[d]); ax.axhline(y0 - 0.5, color='white', lw=2)
    ax.text(len(gl) - 0.35, y0 + k / 2 - 0.5, lab[d], rotation=0, va='center', ha='left', fontsize=6)
    y0 += k
ax = axs[1]
dsn = {'GSE116256': ('van Galen, blasts', 'blast (malignant)'), 'GSE241989': ('GSE241989, HSPC/blast-like', 'HSPC/blast-like'),
       'GSE207104': ('GSE207104, HSPC/blast-like', 'HSPC/blast-like')}
cols = {'GSE116256': '#7A7A7A', 'GSE241989': '#2E6DB4', 'GSE207104': '#5DA5DA'}
mids = {'[0.5,1)': 0.7, '[1,2)': 1.4, '[2,5)': 3.2, '[5,10)': 7.1}
for k, (name, comp) in dsn.items():
    c = cal[(cal.dataset == k) & (cal.compartment == comp)]
    x = [mids[b] for b in c.tpm_bin]
    ax.plot(x, c.median_detect * 100, '-o', ms=3, color=cols[k], lw=1)
    ax.fill_between(x, c.q25 * 100, c.q75 * 100, color=cols[k], alpha=0.12, lw=0)
    ax.text(x[-1] * 1.08, c.median_detect.iloc[-1] * 100 + {'GSE116256': -0.6, 'GSE241989': 0, 'GSE207104': 0.6}[k], name, fontsize=6, color=cols[k], va='center')
ax.axhline(5, color=FOCAL, lw=0.8, ls='--'); ax.text(0.62, 5.4, 'pre-registered 5% threshold', fontsize=6, color=FOCAL)
ax.set_xscale('log'); ax.set_xticks([0.7, 1.4, 3.2, 7.1]); ax.set_xticklabels(['0.5–1', '1–2', '2–5', '5–10'])
ax.set_xlim(0.55, 60); ax.set_ylim(0, 17)
ax.set_xlabel('Gene expression in bulk Beat AML blasts (TPM bin)'); ax.set_ylabel('Cells detecting the gene (%)\nmedian, IQR across genes')
axs[0].text(-0.55, 1.06, 'a', transform=axs[0].transAxes, fontsize=10, fontweight='bold'); letter(axs[1], 'b')
save(fig, 'Figure2_single_cell_detection_and_calibration_v2')

# ---------------- Figure 3: stroma / THP-1 bound + DepMap ----------------
st = T('step5_stroma_sensitivity_bound.csv'); th = T('step8_pma_thp1_expression.csv')
dep = pd.read_csv(os.path.join(TABLES, 'gate03_depmap_gene_effect_per_line.csv'), index_col=0)
fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7), gridspec_kw={'width_ratios': [1.1, 1], 'wspace': 0.4})
ax = axs[0]
items = [('AML MSC (GSE107490, n=3)', st.iloc[0].exp_1TPM, st.iloc[0].obs), ('AML core biopsies (GSE212615, n=22)', st.iloc[2].exp_1TPM, st.iloc[2].obs)]
for _, r in th.dropna(subset=['expected_reads_at_1TPM']).iterrows():
    items.append((f"PMA-THP-1 {r.dataset.split(' (')[0]} (n={r.dataset.split('n=')[1].split(')')[0]})", r.expected_reads_at_1TPM, r.MRGPRD_reads_total))
for i, (name, e, o) in enumerate(items):
    ax.plot([e * 0.1, e], [i, i], color=GREY, lw=1)
    ax.scatter([e], [i], marker='o', s=22, facecolor='white', edgecolor='k', zorder=3)
    ax.scatter([e * 0.1], [i], marker='o', s=14, facecolor='white', edgecolor=GREY, zorder=3)
    ax.scatter([max(o, 0.3)], [i], marker='D', s=20, color=FOCAL, zorder=4)
ax.set_xscale('log'); ax.set_xlim(0.2, 60); ax.set_xticks([0.3, 1, 10, 50]); ax.set_xticklabels(['0', '1', '10', '50'])
ax.set_yticks(range(len(items))); ax.set_yticklabels([n for n, _, _ in items]); ax.invert_yaxis()
ax.set_xlabel('$\\it{MRGPRD}$ reads (pooled over samples)')
ax.scatter([], [], marker='D', color=FOCAL, s=20, label='observed'); ax.scatter([], [], marker='o', facecolor='white', edgecolor='k', s=22, label='expected at 1 TPM')
ax.scatter([], [], marker='o', facecolor='white', edgecolor=GREY, s=14, label='expected at 0.1 TPM')
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.2), ncol=3, handletextpad=0.2, columnspacing=1.0)
ax = axs[1]
gs = ['MRGPRD', 'MAS1', 'MYB', 'SPI1', 'ALB']; aml = dep.OncotreePrimaryDisease.eq('Acute Myeloid Leukemia')
heme = dep.OncotreeLineage.isin(['Myeloid', 'Lymphoid'])
for i, g in enumerate(gs):
    o = dep.loc[~heme, g].dropna(); a = dep.loc[aml, g].dropna()
    vp = ax.violinplot([o], positions=[i - 0.18], widths=0.32, showextrema=False)
    for b in vp['bodies']: b.set_facecolor(NEG); b.set_alpha(0.7); b.set_edgecolor('none')
    ax.scatter(i + 0.18 + rng.uniform(-0.08, 0.08, len(a)), a, s=4, color=GCOL[g], lw=0, alpha=0.9)
    ax.plot([i + 0.06, i + 0.30], [a.median()] * 2, color='k', lw=1.1)
ax.axhline(-0.5, color=GREY, ls='--', lw=0.6); ax.text(-0.45, -0.53, 'dependency rule (−0.5)', fontsize=6, color=GREY, ha='left', va='top')
ax.axhline(0, color='k', lw=0.4)
ax.set_xticks(range(len(gs))); ax.set_xticklabels([f'$\\it{{{g}}}$' for g in gs])
ax.set_ylabel('CRISPR gene effect (Chronos)'); print('N_AML_LINES', int(aml.sum()), 'N_NONHEME', int((~heme).sum()))
for a, l in zip(axs, 'ab'): letter(a, l)
save(fig, 'Figure3_stroma_THP1_bounds_and_dependency_v2')

# ---------------- Figure 4: survival + drug ----------------
cx = T('gate04_cox.csv'); e22 = T('step3_eln2022_cox.csv'); dr = T('step4_drug_panel_spearman.csv'); ex = T('step4_drug_panel_EXPLORATORY_monocytic_adjusted.csv')
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), gridspec_kw={'width_ratios': [1.05, 1, 1], 'wspace': 0.45})
ax = axs[0]; rows = []
for coh in ['BeatAML', 'TCGA-LAML']:
    for mdl in ['RAS8', 'RAS7']:
        a = cx[(cx.cohort == coh) & (cx.model == mdl)].iloc[0]; b = e22[(e22.cohort == coh) & (e22.model == mdl)].iloc[0]
        rows.append((f'{coh.replace("BeatAML", "Beat AML")} {mdl}', a.HR_per_SD, a.lo, a.hi, b.HR_per_SD, b.lo, b.hi))
for i, (nm, h, l, u, h2, l2, u2) in enumerate(rows):
    ax.plot([l, u], [i - 0.12] * 2, color='k', lw=1); ax.scatter([h], [i - 0.12], color='k', s=14, zorder=3)
    ax.plot([l2, u2], [i + 0.12] * 2, color=GREY, lw=1); ax.scatter([h2], [i + 0.12], color='white', edgecolor=GREY, s=14, zorder=3)
ax.axvline(1, color='k', lw=0.5, ls=':')
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
ax.set_xlabel('HR per SD of RAS score (95% CI)')
ax.scatter([], [], color='k', s=14, label='ELN2017 / RISK_MOLECULAR'); ax.scatter([], [], color='white', edgecolor=GREY, s=14, label='ELN 2022 (derived)')
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.22), ncol=1)
ax = axs[1]; d8 = dr[dr.score == 'RAS8'].sort_values('rho').reset_index(drop=True)
sig = d8.q_BH < 0.05
ax.scatter(d8.index[~sig], d8.rho[~sig], s=6, color=NEG, lw=0)
ax.scatter(d8.index[sig], d8.rho[sig], s=6, color=POS, lw=0)
for nm in ['Venetoclax', 'Cytarabine']:
    r = d8[d8.inhibitor == nm]
    if len(r): ax.annotate(nm, (r.index[0], r.rho.iloc[0]), xytext=(-45, -28 if nm == 'Venetoclax' else -14), textcoords='offset points', fontsize=6,
                           arrowprops=dict(arrowstyle='-', lw=0.5, color='k'))
ax.axhline(0, color='k', lw=0.4); ax.set_xlabel(f'Inhibitors ranked by ρ (n = {len(d8)})'); ax.set_ylabel('Spearman ρ, RAS score vs ex vivo AUC')
ax.text(0.02, 0.97, f'{int(sig.sum())}/{len(d8)} at q < 0.05', transform=ax.transAxes, fontsize=6, color=POS, va='top')
ax.text(0.98, 0.03, 'higher AUC = more resistant', transform=ax.transAxes, fontsize=6, color=GREY, ha='right')
ax = axs[2]
ax.scatter(ex.raw_rho_RAS8, ex.partial_rho_RAS8, s=6, color=np.where(ex.partial_q < 0.05, POS, NEG), lw=0)
lim = [-0.45, 0.45]; ax.plot(lim, lim, color=GREY, lw=0.6, ls='--'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('ρ, unadjusted'); ax.set_ylabel('ρ, adjusted for monocytic score')
ax.text(0.02, 0.97, f'{int((ex.partial_q < 0.05).sum())}/{len(ex)} at q < 0.05\nafter adjustment (exploratory)', transform=ax.transAxes, fontsize=6, va='top', color=POS)
for a, l in zip(axs, 'abc'): letter(a, l)
save(fig, 'Figure4_survival_and_drug_response_v2')
print('figures written to', FIGS)
