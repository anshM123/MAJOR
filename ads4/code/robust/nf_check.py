"""Dense floating-point check of a no-forward functional (NOFORWARD_FUNCTIONAL.md section 3) [NE].
Weights evaluated in factored form (1-p)^r p X(p) (see nf_lp.basis_values).  Grid: 150 u in [1e-6, 1] (dense near 1),
J <= min(JCAP, (BMAX/sqrt(u) - 1)/2) for u < 1/4 and J <= JHI for u >= 1/4 (large J at m ~ M: the backward region of d_44).
Usage: python3 nf_check.py NAME [JHI]"""
import sys, os, json
import numpy as np
from fractions import Fraction as Fr
from rcommon import nodes, nquad, flat_integrals, jmax_for, SIGMA, peval
from nf_lp import CDIR

def load_nf(name):
    for fn in [os.path.join(os.path.dirname(os.path.abspath(__file__)), 'functionals_noforward.json'), os.path.join(CDIR, 'functionals_nf.json')]:
        if os.path.exists(fn):
            D = json.load(open(fn))
            if name in D:
                d = D[name]
                return d['r'], np.array([float(Fr(c)) for c in d['X']]), np.array([float(Fr(c)) for c in d['Y']])
    raise KeyError(name)

def parts(r, X, Y, h, u, Jmax):
    n = nquad(Jmax) + (400 if u > 0.9 else 0)
    p, w = nodes(n)
    base = (1 - p) ** r * p
    v2 = base * peval(X, p); v3 = base * peval(Y, p)
    W = np.vstack([w * v2 * (2 - u * p * p), w * SIGMA[h] * u * v3]).T
    Js, R = flat_integrals(h, u, Jmax, W, p)
    if h == '++':
        s = Js % 2 == 0; Js, R = Js[s], R[s]
    return Js, R[:, 0], R[:, 1]

if __name__ == '__main__':
    name = sys.argv[1]; JHI = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    r, X, Y = load_nf(name)
    us = np.unique(np.concatenate([np.logspace(-6, np.log10(0.25), 70), np.linspace(0.25, 1, 61), np.linspace(0.97, 1, 13)]))
    for h in SIGMA:
        worst = (np.inf, None); thr = {C: 0.0 for C in [1, 3, 10]}
        for u in us:
            Jm = JHI if u >= 0.25 else jmax_for(u, 200.0, 3000)
            Js, B2, B3 = parts(r, X, Y, h, u, Jm)
            FF = B2 + B3; N = np.abs(B2) + np.abs(B3); rr = FF / N
            i = np.argmin(rr)
            if rr[i] < worst[0]:
                worst = (rr[i], dict(m=round(1 / np.sqrt(u), 4), J=int(Js[i]), bhat=round((2 * Js[i] + 1) * np.sqrt(u), 2)))
            for C in thr:
                if np.any(FF - C * u * N < 0): thr[C] = max(thr[C], 1 / np.sqrt(u))
        print('%s %s: min FF/N = %.5f at %s ; flat envelope thresholds n(C) (C=1,3,10): %s' % (name, h, worst[0], worst[1],
              ' '.join('%.3f' % thr[C] for C in thr)), flush=True)
