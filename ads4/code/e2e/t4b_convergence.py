"""Numerical-error controls for the tensor Regge moments: dependence on the radial truncation N (the unconverged region of
small u = sqrt(rhobar/rho) is dropped, see moments.series_converged) and on the quadrature grid.  Usage: python3 t4b_convergence.py"""
import sys, time, numpy as np, mpmath as mp
from t4_spinning_e2e import run_case
from cases import cases_fixed_beta, info
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
fine = dict(upanels=((0.0, 0.01, 10), (0.01, 0.03, 14), (0.03, 0.1, 18), (0.1, 0.3, 22), (0.3, 1.0, 32)), tpanels=((0.0, 6.0, 24), (6.0, 28.0, 36), (28.0, 85.0, 28)))
for l, mlist in ((1, (10, 20)), (2, (10, 20))):
    for m in mlist:
        Delta, J = cases_fixed_beta(1.0, [m])[0]
        ref = None
        for tag, N, gk in (("N= 800 Delta", int(800*Delta), {}), ("N= 400 Delta", int(400*Delta), {}), ("N= 200 Delta", int(200*Delta), {}),
                           ("N= 400 Delta, fine grid", int(400*Delta), fine)):
            if N > 12000: continue
            t0 = time.time()
            res = run_case(l, Delta, J, 1.75, (2, 4, 6), (0.0, 3.0), N=N, **gk)
            if ref is None: ref = res
            row = []
            for k in (2, 4, 6):
                for nu in (0.0, 3.0):
                    dev = max(abs(res[(k, nu, q1, q3)]/ref[(k, nu, q1, q3)] - 1) for (q1, q3) in ((l, -l), (l, l)))
                    row.append(f"k={k},nu={nu:.0f}: {float(dev):.1e}")
            pr(f"  l={l} Delta={Delta:.4f} J={J} [{tag:24s}] max rel. deviation from the first line (largest N) over the modes (l,-l),(l,l):  " + "  ".join(row) + f"   ({time.time()-t0:.0f} s)")
open('out_t4b_convergence.txt', 'w').write("\n".join(out) + "\n")
