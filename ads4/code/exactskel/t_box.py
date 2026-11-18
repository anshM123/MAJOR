"""Validation of the box contour (skelexact2.values_box) against skelexact (contour left of tau) at Delta_g = 100 states, and
against the CHMRS series at small nu at Delta_g = 300 / 1000."""
import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../uniform')
import skel
from skelexact2 import ExactState2
from skelexact import ExactState
pn = (np.polynomial.legendre.leggauss(64)[0] + 1)/2
Lam = np.sqrt(98*99)
idx = [3, 20, 35, 45, 52, 58, 63]
nus = Lam*pn[idx]
for (D, J, Df) in (((100.05, 0, 3.0), (138.0, 2, 3.0), (138.0, 5, 2.0)) if len(sys.argv) == 1 else ()):
    ref = ExactState(D, J, Df).values(nus)
    t = time.time(); es = ExactState2(D, J, Df, numax=nus.max()); v = es.values_box(nus)
    print(f"Dg=100 D={D} J={J} Df={Df}: {time.time()-t:.0f}s")
    for k in (2, 3, 6, 7):
        print(f"   k={k}: c1 {np.round(es.diag[k][1], 1)} canc {np.array2string(es.diag[k][0], precision=1)} |box/left-1| {np.array2string(np.abs(v[k]/ref[k]-1), precision=1)}")
for Dg in [float(a) for a in sys.argv[1:]]:
    L = np.sqrt((Dg-2)*(Dg-1))
    for D, J in ((Dg + 0.05, 0), (Dg*1.05, 4)):
        t = time.time(); es = ExactState2(D, J, 3.0, numax=L); v = es.values_box(L*pn[idx] if False else L*pn)
        print(f"Dg={Dg:g} D={D} J={J}: {time.time()-t:.0f}s")
        for k in (2, 5, 7):
            s = [v[k][i]/skel.series_value('seriesC', k, 3.0, L*pn[i], D, J) - 1 for i in (2, 6, 10)]
            print(f"   k={k}: vs series {np.array2string(np.array(s), precision=1)}  max canc {es.diag[k][0].max():.1e}  c1 max {es.diag[k][1].max():.1f}"
                  f"  x=1-: m^2k C = {v[k][-1]:.6g}")
