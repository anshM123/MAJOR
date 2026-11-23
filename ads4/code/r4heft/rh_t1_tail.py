"""Regge tail of the non-dispersive class: decay exponent p(c) of the symmetrised residues r_n(d13 = c) ~ n^{-p}, class 1243 (k = 2, D4)
and, for contrast, class s (k = 2, D8; dispersive: r_n = 0 for n > 4).  A sum over n that converges only for Re d13 > d_c makes
Phi(d13) singular at d_c, i.e. an extra y^{-d_c} in the collinear limit.  Usage: python3 rh_t1_tail.py <dphi>"""
import sys
import numpy as np
import rh_common as R
import rh_t1_collinear as T

dphi = float(sys.argv[1])
a, b, p2, p4 = T.cfg(dphi)
Sab = R.setup(2, a, b, p2, p4, classes=('1243', 's')); Sba = R.setup(2, b, a, p2, p4, classes=('1243', 's'))
L = []
for c in (1.2, 0.9, 2/3, 0.5, 0.35, 0.2, 0.05):
    r = sum((1.0 if o == 'ab' else 1.0)*R.residues(S, R.poly_D(4), c + 0j, 160)['1243'] for o, S in (('ab', Sab), ('ba', Sba)))
    n = np.arange(40, 161)
    p = -np.polyfit(np.log(n), np.log(np.abs(r[40:161]) + 1e-300), 1)[0]
    p2_ = -np.polyfit(np.log(n[-40:]), np.log(np.abs(r[121:161]) + 1e-300), 1)[0]
    L.append(f'dphi {dphi} class 1243, d13 = {c:.3f}: residue decay exponent p (fit n in [40,160]) = {p:.3f}; (n in [121,160]) = {p2_:.3f}; |r_160| = {abs(r[160]):.2e}')
    print(L[-1], flush=True)
open(f'{R.CACHE}/t1_tail_{dphi}.txt', 'w').write('\n'.join(L) + '\n')
