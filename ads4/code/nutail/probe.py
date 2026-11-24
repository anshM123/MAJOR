"""
probe.py (NU_TAIL_CHECK.md): exact CHMRS per-nu action m^{2k}C_{k,nu} of ONE heavy scalar block (Df = 3), for nu from below to
above its mass m, by exactskel/exactmom.py (read-only import).  Prints value, flat bulk-point form E_k(x)P_J(1-2x^2) (x < 1),
the cancellation ratio and the contour index; repeats with a shifted contour (jshift) and a finer step to test stability.
usage: python3 -B probe.py Delta J k nu1,nu2,...  [jc offsets from j0 (<=0) or absolute (>0), comma list]  [h]
"""
import os, sys, time, warnings
sys.dont_write_bytecode = True
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import exactmom as X
from scipy.special import eval_legendre

D = float(sys.argv[1]); J = int(sys.argv[2]); k = int(sys.argv[3])
nus = np.array([float(s) for s in sys.argv[4].split(',')])
jcs = [int(a) for a in sys.argv[5].split(',')] if len(sys.argv) > 5 else [None]
h = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
m = np.sqrt((D - J - 2)*(D + J - 1))
t0 = time.time()
hc = X.HeavyC(D, J, 3.0, h=h)
res = []
for jc in jcs:
    jcc = None if jc is None else (hc.contour(k)[2] + jc if jc <= 0 else jc)
    res.append(hc.value(k, nus, diag=True, jc=jcc))
v, c, jj = res[0]
spread = np.max(np.abs(np.array([r[0] for r in res])/v - 1), axis=0) if len(res) > 1 else np.zeros(nus.size)
for nu, vv, cc, j, sp in zip(nus, v, c, jj, spread):
    x = nu/m
    bp = (2 - x*x)/(1 - x*x)**(1 + k/2)*eval_legendre(J, 1 - 2*x*x) if x < 1 else float('nan')
    print(f"D={D} J={J} k={k} m={m:.3f} nu={nu:8.2f} x={x:.4f}  exact {vv:+.8e}  flatBP {bp:+.4e}  canc {cc:.1e}  jc {j}  spread {sp:.1e}", flush=True)
print(f"time {time.time()-t0:.1f}s", flush=True)
