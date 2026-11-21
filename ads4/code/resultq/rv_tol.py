"""
rv_tol (RESULT_Q_REVISION.md §1.4): tolerances of the delta = 0.2 Result Q functional at its binding states against a UNIFORM graviton
correction Gamma (c2 units, as in rv_smear) confined to a window, for the windows that no computation covers:
   W1: k = 2, x > 0.8 (edge window, held in Result D's chain)        W2: k = 3, x > 0.8
   W3: k = 3, nu < 15 (below the fold; no k = 3 engine)              W4: k = 2, nu < 41 (engine window; for scale)
Gamma_tol = FF_s / |weight|, weight = d(FF)/d(Gamma) = sum_p w psi_k env / m^2 (x u for k = 3), FF_s = FF with the measured forward s_d.
Usage: python3 rv_tol.py -> cache/rv_tol.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pickle
from rv_smear import env

if __name__ == '__main__':
    import stage6_model as M
    F, c2, c3 = M.load('Q3-full-delta0.2')
    psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
    psi3 = lambda p: np.polynomial.polynomial.polyval(p, c3)
    S = pickle.load(open(os.path.join(HERE, 'cache', 'rv_smear.pkl'), 'rb'))
    lines = ['# rv_tol: uniform-Gamma tolerances (c2) per window at the binding states (delta = 0.2; FF incl. measured s_d)']
    x_, w_ = np.polynomial.legendre.leggauss(4000); pn = 0.5*(x_ + 1); pw = 0.5*w_
    for (Dg, h, J, mfl), r in sorted(S.items()):
        Lam = M.lam(Dg); u = 1/mfl**2; m = mfl*Lam; x = pn*np.sqrt(u); nu = pn*Lam
        wins = {'W1 k2 x>0.8': (2, x > 0.8), 'W2 k3 x>0.8': (3, x > 0.8), 'W3 k3 nu<15': (3, nu < 15), 'W4 k2 nu<41': (2, nu < 41)}
        out = []
        for nm, (k, sel) in wins.items():
            if k == 2:
                wgt = np.sum((pw*psi2(pn)*env(h, x, 2))[sel])/m**2
            else:
                wgt = np.sum((pw*psi3(pn)*u*env(h, x, 3))[sel])/m**2
            tol = r['FFs']/abs(wgt) if wgt != 0 else np.inf
            out.append(f"{nm}: weight {wgt:+.2e}, Gamma_tol {'-' if wgt > 0 else '+'}{tol:.1e}")
        lines.append(f"  Dg={Dg} {h} J={J} m/Lam={mfl:.4f} (FF_s/N = {r['FFs']/r['N']:+.4f}): " + "; ".join(out))
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'rv_tol.txt'), 'w').write('\n'.join(lines) + '\n')
