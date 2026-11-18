import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../uniform')
from exactmom import HeavyC
import skel
for (Delta, J, Df) in ((150.7, 2, 3.0), (150.7, 5, 3.0), (160.3, 10, 2.0), (140.2, 4, 2.5)):
    t0 = time.time()
    hc = HeavyC(Delta, J, Df)
    m = np.sqrt(hc.m2)
    nus = np.array([0.3, 3.0, 8.0, 15.0])
    for k in (2, 3, 4, 5, 6, 7):
        v, c, jj = hc.value(k, nus, diag=True)
        s = np.array([skel.series_value('seriesC', k, Df, nu, Delta, J) for nu in nus])
        print(f"D={Delta} J={J} Df={Df} m={m:.1f} k={k}: " + " ".join(f"{vv:.11f}({vv/ss-1:+.1e},c{cc:.0e})" for vv, ss, cc in zip(v, s, c)))
    print("  time %.1f" % (time.time()-t0))
