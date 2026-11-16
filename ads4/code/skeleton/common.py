"""
Shared helpers for the skeleton test (SKELETON_TEST.md).  Imports, unchanged, ads4/code/e2e (exact scalar blocks, grids, formulas),
ads4/code/fourpt (two-pair <JJJJ> partial waves, kernels, leading prediction) and ads4/code/uniform (CHMRS order-8 series).
Nothing outside ads4/code/skeleton/ is modified.  Caches: ads4/code/skeleton/cache/ (git-ignored).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
for sub in ('e2e', 'fourpt', 'uniform'):
    p = os.path.join(CODE, sub)
    if p not in sys.path:
        sys.path.insert(0, p)
CACHE = os.path.join(HERE, 'cache')
os.makedirs(CACHE, exist_ok=True)

import numpy as np, mpmath as mp

UPANELS = ((0.0, 0.01, 8), (0.01, 0.03, 10), (0.03, 0.1, 14), (0.1, 0.3, 16), (0.3, 1.0, 24))
TPANELS = ((0.0, 6.0, 16), (6.0, 28.0, 26), (28.0, 85.0, 20))

def grid_kw(fu=1.0, ft=1.0, usplit=1):
    """e2e.moments.make_grid keywords with the node counts of the u (t) panels multiplied by fu (ft); usplit > 1 splits every
    u panel into usplit equal sub-panels (each with the scaled node count / usplit, at least 8)."""
    up = []
    for (a, b, n) in UPANELS:
        nn = int(round(n*fu))
        if usplit == 1:
            up.append((a, b, nn))
        else:
            e = np.linspace(a, b, usplit + 1)
            for i in range(usplit):
                up.append((float(e[i]), float(e[i+1]), max(8, int(round(nn/usplit)))))
    tp = tuple((a, b, int(round(n*ft))) for (a, b, n) in TPANELS)
    return dict(upanels=tuple(up), tpanels=tp)

def info(Delta, J):
    """m, beta, eta_AdS (mp)"""
    D = mp.mpf(Delta)
    al, alb = D - J - 2, D + J - 1
    m = mp.sqrt(al*alb)
    return m, mp.log(alb/al), (al*al + alb*alb)/(2*al*alb)

def case(beta, J):
    """(Delta, J) with integer J at impact parameter beta (Delta rounded to a double)."""
    m = (2*J + 1)/(2*np.sinh(beta/2))
    al, alb = m*np.exp(-beta/2), m*np.exp(beta/2)
    return float((al + alb + 3)/2), J
