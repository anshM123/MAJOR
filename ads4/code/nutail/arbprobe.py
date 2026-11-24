"""
arbprobe.py (NU_TAIL_CHECK.md): exact CHMRS per-nu action m^{2k}C_{k,nu} of one heavy scalar block (Df = 3) by the contour-free
Arb Laplace/Barnes series of exactskel/skelexact3.py (read-only import), for nu from below to far above the block's mass.
Prints the value, Arb's relative error radius and (x < 1) the flat bulk-point form.
usage: python3 -B arbprobe.py Delta J k nu1,nu2,... [prec] [N0] [nq] [Mextra]
"""
import os, sys, time, warnings
sys.dont_write_bytecode = True
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import skelexact3 as S3
from scipy.special import eval_legendre

D = float(sys.argv[1]); J = int(sys.argv[2]); k = int(sys.argv[3])
nus = [float(s) for s in sys.argv[4].split(',')]
prec = int(sys.argv[5]) if len(sys.argv) > 5 else 400
N0 = int(sys.argv[6]) if len(sys.argv) > 6 else 120
nq = int(sys.argv[7]) if len(sys.argv) > 7 else 160
Mx = int(sys.argv[8]) if len(sys.argv) > 8 else 0
m = np.sqrt((D - J - 2)*(D + J - 1))
numax = max(nus)
M = int(3*numax**2/(2*D) + 6*np.sqrt(numax + 1) + 60) + Mx
t0 = time.time()
L = S3.Laplace3(D, J, 3.0, ks=(k,), prec=prec, N0=N0, nq=nq, M=M)
t1 = time.time()
for nu in nus:
    v = L.value(k, nu)
    mid = float(v.mid()); rad = float(v.rad())
    x = nu/m
    bp = (2 - x*x)/(1 - x*x)**(1 + k/2)*eval_legendre(J, 1 - 2*x*x) if x < 1 else float('nan')
    print(f"D={D} J={J} k={k} m={m:.4f} nu={nu:9.3f} x={x:.5f}  exact {mid:+.10e}  relrad {rad/abs(mid) if mid else rad:.1e}  "
          f"flatBP {bp:+.4e}  (prec {prec}, N0 {N0}, nq {nq}, M {M})", flush=True)
print(f"build {t1-t0:.1f}s, eval {time.time()-t1:.1f}s", flush=True)
