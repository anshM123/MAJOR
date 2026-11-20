"""
fw_scalfit (FOLD_WINDOW_LARGEJ.md §2): control of the smooth-part extraction on the EXACT scalar skeleton at the same state and the same
fold-window nodes: Gamma_s = m^2 (rho_2 + 2[P_{-1/2+i nu}(cosh b) - J0((2J+1)x)]/E_2) (sm_exact, rho_2 barycentric-interpolated at the
nodes, n3_exact.rho_interp), minus the kappa image (C2, ++), fitted with fw_gamma.fit (bases cs, cs+x2, cs+lin) -> r_s (no smooth part
expected beyond the fold plateau m^2 rho_2(nu -> 0) ~ 6.5 c2 of SMOOTH_BOUND §2).
Usage: python3 fw_scalfit.py J:Delta nu1,nu2,...
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import fw_gamma as FG
import fw_scalar as FS
import n3_exact as NX
import remk
from scipy.special import j0

spec, nus = sys.argv[1], np.array([float(v) for v in sys.argv[2].split(',')])
J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
r = FS.sm_rec(D, J)
m = r['m']; x = nus/m
b = 2*np.arcsinh((J + .5)/m)
rho = np.array([NX.rho_interp(v + 0j, r, NX.Lam(100)).real for v in nus])
P = remk.Pders_vec(nus, np.cosh(b), nmax=0)[0]
E2 = (2 - x*x)/(1 - x*x)**2
Gs = m*m*(rho + 2*(P - j0((2*J + 1)*x))/E2)
img = FG.CF.curv_image('C2', '++', J, m, x, FG.Q.A.env, FG.Q.A.SIG)
y = Gs - img
print(f"exact scalar J={J} Delta={D:g} ({r['src']}): nodes nu = {nus.tolist()}")
print(f"   m^2 rho_2       = {np.round(m*m*rho, 2).tolist()}")
print(f"   Gamma_s - kappa = {np.round(y, 1).tolist()}")
for bs in ('cs', 'cs+x2', 'cs+lin'):
    c, e, ir, res, rms = FG.fit(x, y, '++', J, bs)
    print(f"   fit [{bs:7s}]: r = {c[ir]:+8.2f} +- {e[ir]:6.2f} (lsq);  resid {res:.2f} of {rms:.2f}")
