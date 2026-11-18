import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../uniform')
from exactmom import HeavyC
import skel
Delta, J, Df = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
nus = np.array([0.3, 2.0, 5.0, 10.0])
for kw in (dict(), dict(N0=400), dict(nq=400), dict(h=0.05), dict(jshift=2), dict(jshift=-1), dict(L=2000.0)):
    t0 = time.time()
    hc = HeavyC(Delta, J, Df, **kw)
    for k in (2, 4, 6):
        v, c, jj = hc.value(k, nus, diag=True)
        s = np.array([skel.series_value('seriesC', k, Df, nu, Delta, J) for nu in nus])
        print(kw, k, " ".join(f"{vv:.13f}({vv/ss-1:+.1e},c{cc:.0e})" for vv, ss, cc in zip(v, s, c)))
    print("  time %.1f" % (time.time()-t0))
