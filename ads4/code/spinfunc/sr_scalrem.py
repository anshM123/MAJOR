"""
sr_scalrem (SUBLEADING_REGGE.md §3): integrated size of the EXACT scalar-skeleton remainder beyond the hybrid (flat + leading-Regge
curvature), per state, relative to N_s:  rem2/N_s = sum p2 E2 rho_2 / N_s (k = 2 part; Df = 3; exact Mellin values of exactskel/e1).
Also its smooth (non-oscillating) part: shape fit rho = c(x)P_J + s(x) x P_J'/(2J+1) + r(x) (deg 3), and the integral of r alone.
Usage: python3 sr_scalrem.py  -> stdout (cache/sr_scalrem.txt)
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import eval_legendre
from spf_common import NearGap, load_states

lines = []
for Dg in (100, 300, 1000):
    NG = NearGap(Dg); F = NG.F; pn, pw = NG.pn, NG.pw
    seen = set()
    for h, DF, r0, ref in load_states(Dg):
        if h != '++' or DF != 3.0 or r0['J'] < 6 or r0['J'] % 2:
            continue
        J, mfl = r0['J'], r0['mfl']
        x = pn/mfl; y = x*x; E2 = (2 - y)/(1 - y)**2
        p2 = F.psi2(pn)*pw
        rho = r0['rho'][2]
        rem = (p2*E2*rho).sum()
        PJ = eval_legendre(J, 1 - 2*y); hh = 1e-6
        dP = (eval_legendre(J, 1 - 2*(x + hh)**2) - eval_legendre(J, 1 - 2*(x - hh)**2))/(2*hh)*x/(2*J + 1)
        sel = x < 0.95
        Ps = [x**k for k in range(4)]
        A = np.array([p*PJ for p in Ps] + [p*dP for p in Ps] + [p for p in Ps]).T
        c = np.linalg.lstsq(A[sel], rho[sel], rcond=None)[0]
        rs = A[:, 8:] @ c[8:]
        remS = (p2*E2*rs).sum()
        Ns = ref['N_s'] if 'N_s' in ref else ref.get('Ns', np.nan)
        lines.append((Dg, J, mfl, r0['beta'], rem/Ns, remS/Ns, np.sqrt(np.mean((r0['m']**2*rho)**2))))
for Dg in (100, 300, 1000):
    for J in sorted(set(l[1] for l in lines if l[0] == Dg)):
        L = [l for l in lines if l[0] == Dg and l[1] == J]
        w = max(L, key=lambda l: abs(l[4]))
        print(f"Dg={Dg:4d} J={J:2d} states {len(L):2d}: max|rem/N_s| {abs(w[4]):.4f} (m/L={w[2]:.2f}, beta={w[3]:.3f}); "
              f"max|smooth part/N_s| {max(abs(l[5]) for l in L):.4f}; rms m^2 rho {max(l[6] for l in L):.0f}")
