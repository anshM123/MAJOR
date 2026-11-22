"""rv2_q_largeJ_tol (REVIEW_UNIFIED): does Result Q's delta=0.2 functional need an H_LJ analogue?
At near-gap states (gap mass) of large J, compute the HYB heavy action FF/N of Q (K=0, no s_d) and the tolerance to a
UNIFORM graviton correction Gamma (c2 units) applied (a) over the whole p-range of the k=2 part, (b) k=3 part, (c) the edge
window x>0.8 (k=2+k=3), with the weights of RESULT_Q_REVISION's rv_tol (same env, same 1/m^2 convention).
Read-only use of resultq modules/caches.  Output: cache/rv2_q_largeJ_tol.txt"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); RQ = os.path.abspath(os.path.join(HERE, '..', 'resultq'))
sys.path.insert(0, RQ); sys.dont_write_bytecode = True
import numpy as np
import stage6_model as M, common_q as CQ
from stage6_scan import BMAX
from rv_smear import env
F, c2, c3 = M.load('Q3-full-delta0.2')
psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2); psi3 = lambda p: np.polynomial.polynomial.polyval(p, c3)
x_, w_ = np.polynomial.legendre.leggauss(6000); pn = 0.5*(x_ + 1); pw = 0.5*w_
out = ['# Result Q (delta=0.2) HYB at the gap mass; Gamma_tol = FF/|weight| for a uniform Gamma (c2) in the window; sign = failing sign']
for Dg in (300, 1000):
    Lam = M.lam(Dg); tab = M.Tables('Q3-full-delta0.2', Lam, BMAX[Dg])
    for J in (8, 16, 24, 32, 40, 48, 64, 100):
        u = float(M.u_gap(J, Dg)); mfl = 1/np.sqrt(u); m = mfl*Lam; x = pn*np.sqrt(u)
        bh = (2*J + 1)*np.sqrt(u)
        if bh > BMAX[Dg]:
            continue
        for h in ('++', '+-'):
            if h == '++' and J % 2:
                continue
            Js, B2, B3 = CQ.parts(h, u, max(J, 4), [c2], [c3]); k = list(Js).index(J)
            FF, N, _, parts = M.hyb(F, h, u, float(J), B2[k, 0], B3[k, 0], tab, Lam)
            FF = float(FF); N = float(N)
            w2 = pw*psi2(pn)*env(h, x, 2)/m**2; w3 = pw*psi3(pn)*u*env(h, x, 3)/m**2
            sel = x > 0.8
            def tol(w, absw=None):
                return ('-' if w > 0 else '+') + '%.1e' % (FF/abs(w)) if w != 0 else 'inf'
            out.append(f'Dg={Dg} {h} J={J:3d} m/Lam={mfl:.4f} bhat={bh:6.2f}: FF/N={FF/N:+.4f} N={N:.2e} | '
                       f'k2 all {tol(w2.sum())} | k3 all {tol(w3.sum())} | k2+k3 x>0.8 {tol(w2[sel].sum() + w3[sel].sum())} | adversarial |w| k2+k3 {FF/(np.abs(w2).sum()+np.abs(w3).sum()):.1e}')
            print(out[-1], flush=True)
open(os.path.join(HERE, 'cache', 'rv2_q_largeJ_tol.txt'), 'w').write('\n'.join(out) + '\n')
