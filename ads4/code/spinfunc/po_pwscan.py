"""
po_pwscan (PARITY_ODD_CLOSE.md sect. 1.4): how exact is OO = w^T pref_EE v POINT-WISE (single S^0 point n0 = +1, EE integrand, 3x3 p8/w0
prefactor), as a function of the distance |xi| from the crystal saddle (Laplace-normalised coordinates) and of Delta.
w, v from the saddle point.  Usage: python3 po_pwscan.py J  -> appends cache/po_pwscan.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import po_pwtest as PT
J = int(sys.argv[1])
for Dl, nu in ((100.0, 10.0), (100.0, 30.0), (400.0, 40.0), (1600.0, 160.0)):
    zs, L = PT.saddle_L(nu, Dl, 1.0)
    rng = np.random.default_rng(3); base = rng.standard_normal((30, 7)); base /= np.linalg.norm(base, axis=1)[:, None]
    row = []
    for sc in (0.1, 0.2, 0.4, 0.8, 1.6):
        zl = PT.PW.points(zs, L, sc*base)
        P = PT.BEE(zl, nu, Dl, J, 'EE', 1.0)[2]
        w = np.linalg.lstsq(P[:2, :, 0].T, P[2, :, 0], rcond=None)[0]; v = np.linalg.lstsq(P[:, :2, 0], P[:, 2, 0], rcond=None)[0]
        r = (P[2, 2] - np.einsum('a,c,acn->n', w, v, P[:2, :2]))/P[2, 2]
        row.append(np.median(np.abs(r[1:])))
    line = f'J={J} D={Dl} nu={nu}: median |OO - w EE v|/|OO| at |xi| = 0.1 0.2 0.4 0.8 1.6: ' + ' '.join(f'{x:.1e}' for x in row)
    print(line, flush=True)
    open(os.path.join(HERE, 'cache', 'po_pwscan.txt'), 'a').write(line + '\n')
