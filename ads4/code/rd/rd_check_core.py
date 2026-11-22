"""
rd_check_core: validation of rd_core against the project's code (read-only imports).
 (1) flow-function gradient: rd_core.Flow.grad vs the TPS gradient of spinfunc/n1_direct.phi_w (L1_DEG = 2), at the crystal and at
     random complex points around it; also Re/Im of Phi differences (via l1 TPS constant terms) vs sum_k e_k log f_k.
 (2) factorisation Phi_EE = Phi_scal x Q on the real section D: repro_engine/integrand.point (EE k = 4, own implementation of the
     reproducer) minus [scalar sum_k Re e_k log|f_k| + log Jac_D] minus log|Q| must be a constant (nu = 0 part; A^{i nu/2} is common).
Usage: python3 rd_check_core.py  -> cache/rd_check_core.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ['L1_DEG'] = '2'
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'spinfunc'))
import numpy as np
import rd_core as C

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

import l1_tps as T
import n1_direct as N1

pr('# (1) flow-function gradient: rd_core vs n1_direct.phi_w (TPS, L1_DEG = 2)')
rng = np.random.default_rng(1)
for (nu, Dl) in [(4.0, 300.0), (12.0, 100.0), (2.0, 30.0)]:
    c, a = C.crystal(nu, Dl)
    F = C.Flow(nu, Dl, c)
    for trial in range(3):
        w = np.zeros(7, complex) if trial == 0 else 0.3*(rng.normal(size=7) + 1j*rng.normal(size=7))
        wv = [T.var(i, w[i]) for i in range(7)]
        P = N1.phi_w(wv, c, nu, Dl, 0)
        g_tps, H_tps = T.grad_hess(P)
        z, lj, dz = C.zmap(w[:, None], c)
        g0 = C.geom_plain(z)
        y0 = g0['y0']
        # TPS y0 is the principal sqrt (n0 = +1): align
        g_me, _ = F.grad(w[:, None], y0)
        H_me = F.hess(w, y0[0])
        pr(f"  nu = {nu:5.1f}, Delta = {Dl:6.1f}, w-point {trial}: max|grad diff| = {np.abs(g_me[:, 0] - g_tps).max():.2e} "
           f"(|grad| = {np.abs(g_tps).max():.2e}); max|Hess diff| = {np.abs(H_me - H_tps).max():.2e} (|H| = {np.abs(H_tps).max():.2e})")

pr('')
pr('# (2) factorisation on D: lG_EE(repro) - [scalar log|.| + log Jac_D] - log|Q| = const')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'repro_engine'))
import integrand as RE
for Dl in (30.0, 300.0, 1000.0):
    n = 400
    X1, X2, Y1, Y2 = (rng.gamma(k, size=n)*1.5 for k in (4, 3, 4, 3.5))
    s = rng.uniform(-0.95, 0.95, n); phi = rng.uniform(-1.4, 1.4, n); al = rng.uniform(0.05, 3.09, n)
    P = RE.point(X1, X2, Y1, Y2, s, phi, al, Dl)
    z, y0 = C.raw_from_D(X1, X2, s, Y1, Y2, phi, al, Dl)
    g = C.geom_plain(z, y0)
    F = C.Flow(0.0, Dl, np.ones(7))
    lsc = F.logabs_parts({k: v for k, v in g.items() if not isinstance(v, list)})
    psi = al - phi
    ljac = np.log(np.sqrt(X1*X2)*np.sqrt(Y1*Y2)*np.cos(phi)**2*np.abs(np.sin(psi)))
    Q = C.Qratio(g, Dl)
    d = P['lG'] - (lsc + ljac + np.log(np.abs(Q)))
    pr(f"  Delta = {Dl:7.1f}: spread of the difference {d.max() - d.min():.2e} (mean {d.mean():+.6f}); "
       f"spread of log|Q| itself {np.log(np.abs(Q)).max() - np.log(np.abs(Q)).min():.1f}; "
       f"max|Im Q/Q| = {np.abs(Q.imag/Q.real).max():.1e}; mean - log(1/Delta) = {d.mean() + np.log(Dl):+.2e}")
open(os.path.join(HERE, 'cache', 'rd_check_core.txt'), 'w').write('\n'.join(lines) + '\n')
