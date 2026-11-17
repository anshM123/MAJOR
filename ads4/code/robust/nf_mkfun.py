"""Rationalise a no-forward LP solution (cache/nf/pareto_r{r}_d{d}.json) into exact cofactor coefficients:
psi2 = (1-p)^r p X(p), psi3 = (1-p)^r p Y(p) (X, Y rational, scaled so that k ~ 27.485 as for F_g3).  Writes cache/nf/functionals_nf.json
(untracked; the certified one is copied by hand into the tracked functionals_noforward.json).  Usage: python3 nf_mkfun.py NAME r d delta"""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
from numpy.polynomial import legendre as Lg, polynomial as Pn
from nf_lp import CDIR

def cofactor(v, d):
    X = np.zeros(d + 1)
    for i in range(d + 1):
        li = Lg.leg2poly(np.eye(d + 1)[i]); sh = np.zeros(1)
        for k, ck in enumerate(li): sh = Pn.polyadd(sh, ck * Pn.polypow([-1., 2.], k))
        X[:len(sh)] += v[i] * sh
    return X

def mul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, s in enumerate(a):
        for j, t in enumerate(b): r[i + j] += s * t
    return r

def expand(r, X):
    omp = [Fr(1)]
    for _ in range(r): omp = mul(omp, [Fr(1), Fr(-1)])
    return mul(mul(omp, [Fr(0), Fr(1)]), X)

if __name__ == '__main__':
    name, r, d, delta = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    suffix = sys.argv[5] if len(sys.argv) > 5 else ''          # e.g. _proxy120
    rows = json.load(open(os.path.join(CDIR, 'pareto_r%d_d%d%s.json' % (r, d, suffix))))
    row = next(x for x in rows if abs(x['delta'] - delta) < 1e-12)
    v = np.array(row['v']); sc = 27.48509823509836
    X = [Fr(int(round(c * sc * 1000)), 1000) for c in cofactor(v[:d + 1], d)]
    Y = [Fr(int(round(c * sc * 1000)), 1000) for c in cofactor(v[d + 1:], d)]
    mom = lambda c, j: sum(cc / (i + j + 1) for i, cc in enumerate(c))
    psi2, psi3 = expand(r, X), expand(r, Y)
    while mom(psi2, 2) - mom(psi3, 0) > Fr(-1, 10000):          # scalar-matter coefficient strictly negative after rounding
        Y[0] += Fr(1, 1000); psi3 = expand(r, Y)
    k = -2 * mom(psi2, 6) + 2 * mom(psi3, 4)
    A = 16 * psi2[1] / k
    Lfin = sum(psi2[i] / (i - 1) for i in range(2, len(psi2)))
    print(name, 'r=%d d=%d X=%s Y=%s  k=%.6f A=%.4f B=%.4f' % (r, d, [str(c) for c in X], [str(c) for c in Y], float(k), float(A), float(16 * Lfin / k)))
    fn = os.path.join(CDIR, 'functionals_nf.json')
    D = json.load(open(fn)) if os.path.exists(fn) else {}
    D[name] = dict(r=r, X=[str(c) for c in X], Y=[str(c) for c in Y], psi2_exact=[str(c) for c in psi2], psi3_exact=[str(c) for c in psi3],
                   a_exact='0', psi2=[float(c) for c in psi2], psi3=[float(c) for c in psi3], a=0.0, source=dict(r=r, d=d, delta=delta))
    json.dump(D, open(fn, 'w'), indent=1)
