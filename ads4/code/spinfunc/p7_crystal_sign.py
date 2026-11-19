"""
p7: a-independent check at the Delta = 100 rows (Delta_g = 100 binding region, m/Lam ~ 1): the crystal's own psi-integrated route-A' action
for the odd (+-) state, F_cr/a = sum psi2 w env2 (D + Gamma2/m^2) + sum psi3 w u env3 (-D + Gamma3/m^2), at the row's own m with
x = p sqrt(u), u = (Lam/m)^2 (Lam of Delta_g = 100).  a > 0 is a common positive factor of C^_2 and C^_3, so the SIGN of F_cr is independent
of the (non-converged) per-state constant; reported with N = |C2 LO part| + |C3 LO part| for K2 and K3 (where computed).
Usage: python3 p7_crystal_sign.py J [J ...]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.interpolate import CubicSpline
from spf_common import NearGap
import p4_odd_analysis as A

NG = NearGap(100); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
for J in [int(a) for a in sys.argv[1:]]:
    for K, src in ((2, False), (2, True), (3, True)):
        kinds = tuple(k for k in ('C2', 'C3H') if A.load(k, J, k3=src) and any(kk[0] == 100.0 for kk in A.load(k, J, k3=src)))
        if 'C2' not in kinds:
            continue
        try:
            g = A.grid(J, K=K, kinds=kinds, k3src=src, Ds=[100.0])[100.0]
        except Exception as ex:
            print(J, K, src, ex); continue
        m = g['m']; u = (Lm/m)**2; x = pn*np.sqrt(u)
        p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
        D = A.Dc('+-', J, x)
        G2 = CubicSpline(g['xs'], g['C2'])(np.minimum(x, g['xs'][-1]))
        f2 = (p2*A.env('+-', x, 2)*(D + G2/m**2)).sum(); n2 = abs((p2*A.env('+-', x, 2)*D).sum())
        if 'C3H' in kinds:
            G3 = CubicSpline(g['xs'], g['C3H'])(np.minimum(x, g['xs'][-1]))
            f3 = (p3*u*A.env('+-', x, 3)*(-D + G3/m**2)).sum(); n3 = abs((p3*u*A.env('+-', x, 3)*D).sum())
        else:
            f3 = (p3*u*A.env('+-', x, 3)*(-D)).sum(); n3 = abs(f3)
        print(f"J={J} K{K} ({'DEG8' if src else 'DEG6'} records, kinds {kinds}): m={m:.1f} (m/Lam={m/Lm:.3f}); F_cr/(a N) = {(f2 + f3)/(n2 + n3):+.4f} "
              f"[C2 part {f2/(n2 + n3):+.4f}, C3 part {f3/(n2 + n3):+.4f}; LO only {((p2*A.env('+-', x, 2)*D).sum() - (p3*u*A.env('+-', x, 3)*D).sum())/(n2 + n3):+.4f}]; a = {g['a']:.4g}")
