"""
n1_ratio_check (doc 24.4): normalisation-free comparison of the state/reference ratio (Delta = 100 row vs Delta_ref = 2e4, same nu,
helicity diagonal elements) from the Laplace chain at K2 and K3 (l9 records) with the direct (non-perturbative) value:
   d_K(x) = m^2 Re[ v_K/v_direct - 1 ]   (c2 units; v = exp(dl) x state/ref ratio, l9_analysis conventions),
and the RQMC error of the direct ratio in the same units.  Usage: python3 n1_ratio_check.py J [J ...] -> cache/n1_ratio_check.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import n1_analysis as N, l9_analysis as A, l8_analysis as L8

lines = []
def pr(s=""):
    print(s, flush=True); lines.append(s)
for J in [int(a) for a in sys.argv[1:]] or [10, 8, 6, 4]:
    D = 100.0; m = A.mphys(D, J)
    Rb = dict(L8.records(J)); Rb.update(A._load('C2', J, '_k2')); Rb.update(A._load('C2', J, ''))
    B = L8.basis(J, D, Rb); Bref = L8.basis(J, A.DREF, Rb)
    for kind in ('C2', 'C3H'):
        db = N.dload(kind, J)
        RCs = {K: A.Recs(kind, J, K) for K in (2, 3)}
        for j, h in enumerate(('++', '+-')):
            row = {2: [], 3: [], 'e': [], 'x': []}
            for (DD, nu), r in sorted(db.items()):
                if DD != D:
                    continue
                rr = N.dget(db, A.DREF, nu)
                if rr is None:
                    continue
                N.NU[0] = nu; Hs, ls = N.Hrep(r, B, D, J); N.NU[0] = nu; Hr, lr = N.Hrep(rr, Bref, A.DREF, J)
                lvd = ls - lr + np.log(Hs.mean(0)[j, j]/Hr.mean(0)[j, j])
                n = min(len(Hs), len(Hr))
                rel = np.sqrt(np.var((Hs[:n, j, j]/Hs.mean(0)[j, j]).real, ddof=1)/n + np.var((Hr[:n, j, j]/Hr.mean(0)[j, j]).real, ddof=1)/n)
                row['x'].append(nu/m); row['e'].append(m*m*rel)
                for K in (2, 3):
                    st, rf = RCs[K].get(D, nu), RCs[K].get(A.DREF, nu)
                    if st is None or rf is None or (K == 3 and np.all(st['M6'] == 0)):
                        row[K].append(np.nan); continue
                    Ms, Mr = A.helmat(st, B, K), A.helmat(rf, Bref, K)
                    dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                    lvK = dl + np.log(A.pred([Mo[j, j] for Mo in Ms], [Mo[j, j] for Mo in Mr], K))
                    row[K].append(m*m*(np.exp(lvK - lvd) - 1).real)
            if row['x']:
                pr(f"J={J:2d} {kind:3s} {h}: x      " + " ".join(f"{v:7.2f}" for v in row['x']))
                pr(f"            K2-dir " + " ".join(f"{v:7.0f}" for v in row[2]))
                pr(f"            K3-dir " + " ".join(f"{v:7.0f}" for v in row[3]))
                pr(f"            err    " + " ".join(f"{v:7.0f}" for v in row['e']))
open(os.path.join(HERE, 'cache', 'n1_ratio_check.txt'), 'w').write("\n".join(lines) + "\n")
