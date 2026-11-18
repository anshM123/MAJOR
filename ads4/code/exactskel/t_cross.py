"""Validation of skelexact2 (contour across the Bhat poles) against skelexact (contour left of tau) at Delta_g = 100 states."""
import sys, time
import numpy as np
sys.path.insert(0, '.')
import skelexact2
from skelexact2 import ExactState2
from skelexact import ExactState
Lam = np.sqrt(98*99)
pn = (np.polynomial.legendre.leggauss(64)[0] + 1)/2
idx = [5, 30, 45, 55, 63]
nus = Lam*pn[idx]
for (D, J, Df) in ((100.05, 0, 3.0), (138.0, 2, 3.0), (138.0, 5, 2.0)):
    ref = ExactState(D, J, Df).values(nus)
    for lc in (6.0, 2.0):
        skelexact2.LOGCANC = np.log(10.0**lc)
        t = time.time()
        es = ExactState2(D, J, Df, numax=nus.max())
        v = es.values(nus)
        print(f"D={D} J={J} Df={Df} target canc 1e{lc:.0f}: {time.time()-t:.0f}s")
        for k in (2, 3, 6, 7):
            print(f"   k={k}: c1 = {np.round(es.diag[k][1], 1)}  canc {np.array2string(es.diag[k][0], precision=1)}  rel diff vs left-of-tau "
                  f"{np.array2string(np.abs(v[k]/ref[k]-1), precision=1)}")
