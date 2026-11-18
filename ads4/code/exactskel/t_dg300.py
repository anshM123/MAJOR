"""Timing and small-nu series check of skelexact2 at Delta_g = 300 near-gap states."""
import sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../uniform')
import numpy as np, skel
from skelexact2 import ExactState2
Dg = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
Lam = np.sqrt((Dg-2)*(Dg-1))
pn = (np.polynomial.legendre.leggauss(64)[0]+1)/2
for D, J in ((Dg + 0.05, 0), (Dg*1.05, 4)):
    t = time.time(); es = ExactState2(D, J, 3.0, numax=Lam); v = es.values(Lam*pn); print(D, J, 'time %.0f' % (time.time()-t), flush=True)
    for k in (2, 5, 7):
        s = [v[k][i]/skel.series_value('seriesC', k, 3.0, Lam*pn[i], D, J)-1 for i in (2, 6, 10)]
        print(k, 'vs series (small nu)', np.array2string(np.array(s), precision=1), ' max canc %.1e' % es.diag[k][0].max(), ' c1 max', es.diag[k][1].max(), flush=True)
