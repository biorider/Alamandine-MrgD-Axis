"""Step 2: MRGPRD read-window uniqueness against MRG-family paralogs and MAS1/MAS1L (Ensembl cDNA, both strands)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import os, time, requests, numpy as np, pandas as pd
from config import TABLES
E = 'https://rest.ensembl.org'; H = {'Content-Type': 'application/json'}
fam = ['MRGPRD', 'MRGPRX1', 'MRGPRX2', 'MRGPRX3', 'MRGPRX4', 'MRGPRE', 'MRGPRF', 'MRGPRG', 'MAS1', 'MAS1L']
seqs = {}
for g in fam:
    j = requests.get(f'{E}/lookup/symbol/homo_sapiens/{g}?expand=1', headers=H, timeout=30).json()
    seqs[g] = [requests.get(f"{E}/sequence/id/{t['id']}?type=cdna", headers=H, timeout=30).json()['seq'].upper() for t in j.get('Transcript', [])]
    time.sleep(0.1)
rc = lambda s: s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
win = lambda s, L: np.lib.stride_tricks.sliding_window_view(np.frombuffer(s.encode(), dtype=np.uint8), L)
target = seqs['MRGPRD'][0]; pool = [s for g in fam[1:] for s in seqs[g]]
res = []
for L in [50, 76, 100]:
    T = win(target, L); P = np.vstack([win(s, L) for s in pool if len(s) >= L] + [win(rc(s), L) for s in pool if len(s) >= L])
    mind = np.full(len(T), L)
    for i in range(0, len(P), 4000): mind = np.minimum(mind, (T[:, None, :] != P[None, i:i + 4000, :]).sum(2).min(1))
    best = [(g, int((T[:, None, :] != np.vstack([win(s, L) for s in seqs[g]] + [win(rc(s), L) for s in seqs[g]])[None]).sum(2).min())) for g in fam[1:]]
    res.append(dict(read_len=L, n_windows=len(T), min_mismatch_any=int(mind.min()), median_min_mismatch=float(np.median(mind)),
                    frac_exact=float((mind == 0).mean()), frac_le4pct=float((mind <= int(L * .04)).mean()), frac_le10pct=float((mind <= int(L * .10)).mean()),
                    closest=sorted(best, key=lambda x: x[1])[:3]))
r = pd.DataFrame(res); r.to_csv(os.path.join(TABLES, 'step2_mrgprd_paralog_mappability.csv'), index=False); print(r.to_string(index=False))
