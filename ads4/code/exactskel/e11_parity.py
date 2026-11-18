"""
e11: tolerance of the near-gap thresholds to a hypothetical O(1/m) term in the heavy moments (PARITY_ODD_HEAVY.md sect. 3).
Reads cache/e10_<Dg>_Kstar_recs.pkl (per-k heavy parts P_k and forward part FW of every exact-skeleton state, D-strong 'U' and D-safe 'N').
Model: every heavy moment of the affected channel(s) multiplied by (1 + c1/m) with unknown sign -> F >= F - |c1|/m (sum_k |P_k| + |FW|).
Reports c1* = largest |c1| that leaves n unchanged (scope: '++' only = parity-odd exchanges of a parity-invariant theory; 'all' = both
channels), and n(c1) for c1 = 2, 4, 8, 16 (delta = 0 and 0.03).   Usage: python3 e11_parity.py  -> cache/e11.txt
"""
import os, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s):
    print(s, flush=True); out.append(s)
pr("# e11: O(1/m) tolerance of the thresholds (exact-skeleton states, K* kernel; |c1| bound on sum_k |P_k| + |FW|)")
for Dg in (100, 300, 1000):
    recs = pickle.load(open(os.path.join(HERE, 'cache', f'e10_{Dg}_Kstar_recs.pkl'), 'rb'))
    GRID = sorted(set(r['mfl'] for r in recs))
    def n_thr(ok):
        for n0 in GRID:
            if all(ok(r) for r in recs if r['mfl'] >= n0 - 1e-9):
                return n0
    for fx, name in (('U', 'D-strong (M02)'), ('N', 'D-safe (NF02)')):
        n0 = n_thr(lambda r: r[fx]['F'] > 0)
        A = lambda r, hs: (np.sum(np.abs(r[fx]['P'])) + abs(r[fx]['FW'])) if r['h'] in hs else 0.0
        cells = []
        for hs, lab in ((('++',), "'++' only"), (('++', '+-'), "all")):
            for delta in (0.0, 0.03):
                v = [(r[fx]['F'] - delta*r[fx]['Ns'])*r['m']/A(r, hs) for r in recs if r['mfl'] >= n0 - 1e-9 and A(r, hs) > 0]
                ns = [n_thr(lambda r: r[fx]['F'] - c1/r['m']*A(r, hs) - delta*r[fx]['Ns'] > 0) for c1 in (2, 4, 8, 16)]
                cells.append(f"{lab}, delta={delta}: c1* = {min(v):.1f}; n(c1=2,4,8,16) = {ns}")
        pr(f"Dg={Dg} {name}: n = {n0}\n   " + "\n   ".join(cells))
open(os.path.join(HERE, 'cache', 'e11.txt'), 'w').write("\n".join(out) + "\n")
