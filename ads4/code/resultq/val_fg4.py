"""Validation: CHLPS example functional F_g4 (arXiv:2201.06602 (3.3b)), completed with +B5|_0 (FLAT_GRAVITON_BOUND.md 1.4).
Reproduces A = 99/7 = 14.14, B (paper: -15.5), kappa = 35249/840840 = 0.042, and checks heavy positivity on a dense grid [NE]."""
import numpy as np
from fractions import Fraction as Fr
from common_q import heavy_action, pmul, mom
from rcommon import jmax_for
o3 = [1., -3., 3., -1.]
psi2 = pmul(o3, [0, 99, 282, 49]); psi3 = pmul(o3, [21, -474, 1216])
a, c = 1.0, 7.0
k = -2 * mom(psi2, 6) + 2 * mom(psi3, 4) + 4 * a
x0 = psi2[1]; Lfin = sum(psi2[i] / (i - 1) for i in range(2, len(psi2)))
print('A = x0/c = %.4f   B = Lfin/c = %.4f   kappa = k/(16c) = %.5f  (k = %.5f)' % (x0 / c, Lfin / c, k / (16 * c), k))
mn = (np.inf, None)
for u in np.unique(np.concatenate([np.logspace(-5, 0, 120), np.linspace(0.9, 1, 21)])):
    for h in ['++', '+-']:
        Js, FF, N = heavy_action(psi2, psi3, a, c, h, u, jmax_for(u, 300.0, 3000))
        i = np.argmin(FF / N)
        if FF[i] / N[i] < mn[0]:
            mn = (FF[i] / N[i], (1 / np.sqrt(u), int(Js[i]), h, FF[i]))
print('min FF/N =', mn)
