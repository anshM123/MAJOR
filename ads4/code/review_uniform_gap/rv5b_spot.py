"""rv5b_spot.py -- mpmath spot check of the double-precision chi table (breakpoints every 1/4 unit; window |v-nu|<=60)."""
import sys, numpy as np, mpmath as mp
from rv_common import *
lam2 = int(sys.argv[1]); mp.mp.dps = 25
d = np.load(f'out/chi_{lam2}.npz'); NU, chi = d['nu'], d['chi']
Lm = mp.sqrt(lam2); L = float(Lm)
def Wm(v):
    v = abs(v)
    return psi2_fact_mp(v/Lm)/mp.re(Fz(v)) if v <= Lm else mp.mpf(0)
for x in [float(a) for a in sys.argv[2:]]:
    i = np.argmin(np.abs(NU-x)); x = NU[i]; xm = mp.mpf(x)
    lo, hi = max(-Lm, xm-60), min(Lm, xm+60)
    pts = sorted(set([lo, hi] + [mp.mpf(k)/4 for k in range(int(4*lo)+1, int(4*hi)+1)] + ([xm] if xm < Lm else [])))
    if lo >= hi: print(f'nu={x:.4f} outside mp window, skipped'); continue
    pts = [p for p in pts if lo <= p <= hi]
    qm = mp.quad(lambda v: Wm(v)*mp.sech(mp.pi*(xm-v)/2), pts)/4
    cm = mp.re(Bz(xm))*qm
    print(f"nu={x:.4f} chi(double)={chi[i]:.15e} chi(mp)={mp.nstr(cm,15)} rel={float(abs(chi[i]-cm)/max(abs(cm),mp.mpf(10)**-300)):.2e}", flush=True)
