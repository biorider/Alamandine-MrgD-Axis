"""Approximate ELN 2022 genetic risk (Döhner et al., Blood 2022) from karyotype strings, fusion calls and mutations.
Sensitivity analysis only — the frozen Phase 0 primary uses ELN2017 / RISK_MOLECULAR.
Declared approximations: ISCN parsed heuristically; complex = >=3 abnormalities in the first abnormal clone;
monosomal karyotype = >=2 autosomal monosomies, or 1 autosomal monosomy plus >=1 structural abnormality;
CEBPA bZIP in-frame approximated (Beat: CEBPA_Biallelic=='bi'; TCGA: in-frame indel/missense at aa 278-345);
TP53 adverse if protein-altering with VAF >= 0.10; APL (PML::RARA) kept Favorable to mirror the cohort ELN2017 coding;
missing/failed karyotype with no class-defining fusion -> 'Unclassifiable' (excluded, complete-case)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import re
import pandas as pd

MR = ['ASXL1', 'BCOR', 'EZH2', 'RUNX1', 'SF3B1', 'SRSF2', 'STAG2', 'U2AF1', 'ZRSR2']
PROTEIN_ALTERING = {'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins',
                    'In_Frame_Del', 'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site',
                    'missense_variant', 'stop_gained', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
                    'splice_acceptor_variant', 'splice_donor_variant', 'start_lost', 'stop_lost', 'protein_altering_variant'}


def parse_karyotype(k):
    """Return dict of cytogenetic flags, or None if karyotype is missing/failed."""
    if k is None or (isinstance(k, float) and pd.isna(k)):
        return None
    s = str(k).lower().replace(' ', '')
    # normalise translocation chromosome order: t(21;8) -> t(8;21)
    def _norm(m):
        a, b = m.group(1), m.group(2)
        ka = (99 if a in 'xy' else int(a)) if a.isdigit() or a in 'xy' else 99
        kb = (99 if b in 'xy' else int(b)) if b.isdigit() or b in 'xy' else 99
        return f't({a};{b})' if ka <= kb else f't({b};{a})'
    s = re.sub(r't\((\w{1,2});(\w{1,2})\)', _norm, s)
    if not re.search(r'\d{2}', s) or any(w in s for w in ['fail', 'nomet', 'insufficient', 'notdone', 'n/a', 'unknown']):
        return None
    clone = s.split('/')[0]
    clone = re.sub(r'\[.*?\]', '', clone)
    parts = [p for p in clone.split(',')[2:] if p]          # drop modal number and sex chromosomes
    abn = [p for p in parts if not re.fullmatch(r'(sl|sdl\d*|idem|c)', p)]
    if '/' in s and not abn:                                # first clone normal, look for any abnormal clone
        for cl in s.split('/')[1:]:
            cl = re.sub(r'\[.*?\]', '', cl); p2 = [p for p in cl.split(',')[2:] if p and p not in ('sl', 'idem')]
            if p2: abn = p2; break
    mono = [p for p in abn if re.fullmatch(r'-(\d{1,2})', p)]
    struct = [p for p in abn if re.match(r'(t|inv|del|der|add|i|dup|ins|r|ider|dic|hsr)\(', p)]
    f = dict(
        t8_21=bool(re.search(r't\(8;21\)', s)),
        inv16=bool(re.search(r'inv\(16\)|t\(16;16\)', s)),
        t15_17=bool(re.search(r't\(15;17\)', s)),
        t6_9=bool(re.search(r't\(6;9\)', s)),
        t9_11=bool(re.search(r't\(9;11\)', s)),
        kmt2a_other=bool(re.search(r'11q23|t\(\d+;11\)\([^)]*q23|t\(11;\d+\)\(q23', s)) and not re.search(r't\(9;11\)', s),
        t9_22=bool(re.search(r't\(9;22\)', s)),
        t8_16=bool(re.search(r't\(8;16\)', s)),
        inv3=bool(re.search(r'inv\(3\)|t\(3;3\)|3q26|3q21q26', s)),
        m5=bool(re.search(r'(^|,)-5(,|$)|del\(5\)\(q|5q-', s)),
        m7=bool(re.search(r'(^|,)-7(,|$)', s)),
        m17=bool(re.search(r'(^|,)-17(,|$)|17p|i\(17\)\(q10\)|del\(17\)\(p', s)),
        n_abn=len(abn),
        monosomal=(len(mono) >= 2) or (len(mono) == 1 and len(struct) >= 1),
    )
    return f


def classify(cyto, fus, npm1, itd, cebpa_bzip, tp53_10, mr_mut):
    """cyto: dict from parse_karyotype or None; fus: set of fusion tokens (upper)."""
    fus = fus or set()
    cbf = any(x in fus for x in ['RUNX1-RUNX1T1', 'CBFB-MYH11']) or (cyto and (cyto['t8_21'] or cyto['inv16']))
    apl = ('PML-RARA' in fus) or (cyto and cyto['t15_17'])
    if cbf or apl:
        return 'Favorable'
    kmt2a_other = ('KMT2A_RE' in fus) or (cyto and cyto['kmt2a_other'])
    t9_11 = ('MLLT3-KMT2A' in fus) or (cyto and cyto['t9_11'])
    adv_fus = any(x in fus for x in ['DEK-NUP214', 'BCR-ABL1', 'KAT6A-CREBBP', 'GATA2-MECOM'])
    if cyto is None and not (adv_fus or (kmt2a_other and not t9_11) or t9_11):
        return 'Unclassifiable'
    adv_cyto = adv_fus or (kmt2a_other and not t9_11) or (cyto is not None and (
        cyto['t6_9'] or cyto['t9_22'] or cyto['t8_16'] or cyto['inv3'] or cyto['m5'] or cyto['m7'] or cyto['m17']
        or cyto['n_abn'] >= 3 or cyto['monosomal']))
    if adv_cyto or tp53_10:
        return 'Adverse'
    if npm1 and not itd:
        return 'Favorable'
    if npm1 and itd:
        return 'Intermediate'
    if cebpa_bzip:
        return 'Favorable'
    if mr_mut:
        return 'Adverse'
    return 'Intermediate'   # includes t(9;11), FLT3-ITD without adverse lesions, all others
