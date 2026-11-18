import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../uniform')
from exactmom import HeavyC
import skel
from scipy.special import eval_legendre
Dg = 100.0; Lam = np.sqrt((Dg-2)*(Dg-1))
def Ek(k, x): return (2-x*x)/(1-x*x)**(1+k/2)
for (Delta, J, Df) in ((100.05, 0, 3.0), (138.0, 2, 3.0)):
    m = np.sqrt((Delta-J-2)*(Delta+J-1))
    nus = Lam*np.array([0.2, 0.4, 0.6, 0.8, 0.9, 0.97, 0.999])
    for kw in (dict(), dict(jshift=-4), dict(h=0.05, N0=600, nq=320)):
        t0 = time.time()
        hc = HeavyC(Delta, J, Df, **kw)
        for k in (2, 5, 7):
            v, c, jj = hc.value(k, nus, diag=True)
            x = nus/m
            bp = Ek(k, x)*eval_legendre(J, 1-2*x*x)
            s = np.array([skel.series_value('seriesC', k, Df, nu, Delta, J) if xx < 0.45 else np.nan for nu, xx in zip(nus, x)])
            print(f"D={Delta} J={J} m={m:.2f} {kw} k={k}: " + " ".join(f"x={xx:.3f}:{vv:.10g}[/BP-1={vv/b-1:+.2e},/ser-1={vv/ss-1:+.1e},c{cc:.0e},j{j}]" for xx, vv, b, ss, cc, j in zip(x, v, bp, s, c, jj)))
        print("  time %.1f" % (time.time()-t0))
