"""Rationalise LP solutions and store them in cache/functionals.json (never in the tracked functionals_certified.json, which holds
the certified M01, M02 and is read first by rcert.py, lowcheck.py, check.py) (exact rational coefficients, scaled by k(F_g3)=27.485...
so that the numbers are comparable with arXiv:2201.06602 (3.3a); the overall scale is irrelevant).
Usage: python3 mkfun.py NAME JSONFILE KEY=VALUE...   (selects the first entry of JSONFILE matching all KEY=VALUE)"""
import sys, json, os
from fractions import Fraction as Fr
import numpy as np
from rcommon import Functional, basis2, basis3
from lp import CACHE

def rationalise(v, d2, d3, scale=27.48509823509836, den=1000):
    w = [Fr(int(round(float(x) * scale * den)), den) for x in v]
    x = w[:d2 + 1]; y = w[d2 + 1:d2 + d3 + 2]; a = w[-1]
    omp = [Fr(1), Fr(-3), Fr(3), Fr(-1)]
    def mul(p, q):
        r = [Fr(0)] * (len(p) + len(q) - 1)
        for i, s in enumerate(p):
            for j, t in enumerate(q): r[i + j] += s * t
        return r
    psi2 = mul(omp, [Fr(0)] + x); psi3 = mul(omp, y)
    # scalar-matter coefficient at m_l -> 0 is mu2_2 - mu3_0 (must be < 0); the LP puts it at 0, so rounding can flip its sign:
    # raise y_0 = psi3(0) in steps of 1/1000 until mu2_2 - mu3_0 <= -1/10000 (changes psi3 by O(1e-3) out of O(10^4)).
    mom = lambda c, j: sum(cc / (i + j + 1) for i, cc in enumerate(c))
    while mom(psi2, 2) - mom(psi3, 0) > Fr(-1, 10000):
        y[0] += Fr(1, 1000); psi3 = mul(omp, y)
    return x, y, a, psi2, psi3

if __name__ == '__main__':
    name, fn = sys.argv[1], sys.argv[2]
    sel = dict(kv.split('=') for kv in sys.argv[3:])
    rows = json.load(open(os.path.join(CACHE, fn)))
    row = next(r for r in rows if all(str(r[k]) == v or (isinstance(r[k], float) and abs(r[k] - float(v)) < 1e-12) for k, v in sel.items()))
    x, y, a, psi2, psi3 = rationalise(row['v'], row['d2'], row['d3'])
    F = Functional([float(c) for c in psi2], [float(c) for c in psi3], float(a), name)
    path = os.path.join(CACHE, 'functionals.json')
    D = json.load(open(path)) if os.path.exists(path) else {}
    D[name] = dict(psi2=[float(c) for c in psi2], psi3=[float(c) for c in psi3], a=float(a),
                   psi2_exact=[str(c) for c in psi2], psi3_exact=[str(c) for c in psi3], a_exact=str(a),
                   x=[str(c) for c in x], y=[str(c) for c in y], source=dict(file=fn, **{k: row[k] for k in row if k != 'v'}))
    json.dump(D, open(path, 'w'), indent=1)
    print(name, 'x =', [str(c) for c in x], 'y =', [str(c) for c in y], 'a =', a, ' k=%.5f A=%.4f B=%.4f B_AdS=%.4f' % (F.k(), *F.AB()))
