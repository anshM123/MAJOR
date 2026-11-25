"""Step 1 (heuristic, not part of the proof): collect the canonical structure vectors of code/kernels for all computed J
and interpolate every coefficient in J.  The resulting candidate vectors are then PROVED to be conserved and
Bose-symmetric for symbolic J by ap2_conservation.py.  Output: cache/symbolic_structures.pkl"""
import json, pickle, os
from common import *
import canonical

DATA = os.path.join(HERE, '..', 'kernels', 'data')

def get(Jv):
    fn = os.path.join(DATA, f'structures_J{Jv}.json')
    out = {}
    if os.path.exists(fn):
        d = json.load(open(fn))
        for name, v in d.items():
            out[name] = {tuple(m[:6]) + (m[6],): sp.sympify(c, locals={'s': S}) for m, c in zip(v['monomials'], v['coefficients'])}
    else:
        for name, (par, monos, vec, K, P) in canonical.canonical(Jv).items():
            out[name] = {tuple(m): sp.factor(c) for m, c in zip(monos, vec)}
    return out

def rel(m, Jv):   # monomial label relative to J: m4 -> J - m4
    return m[:2] + (Jv - m[2],) + m[3:]

def fit(points):
    """points: {Jv: poly in s}; interpolate each s-coefficient in J on all but the last two points; check the rest"""
    Js = sorted(points)
    train, test = Js[:-2], Js[-2:]
    polys = {Jv: sp.Poly(sp.expand(points[Jv]), S) for Jv in Js}
    deg = max(p.degree() for p in polys.values())
    if deg < 0:
        return sp.Integer(0), all(points[Jv] == 0 for Jv in Js)
    res = 0
    for k in range(deg + 1):
        data = [(Jv, polys[Jv].coeff_monomial(S**k)) for Jv in train]
        res += sp.interpolate(data, J)*S**k
    res = sp.factor(sp.expand(res))
    ok = all(sp.expand(res.subs(J, Jv) - points[Jv]) == 0 for Jv in test)
    return res, ok

if __name__ == '__main__':
    even_J = [4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 30]
    odd_J = [5, 7, 9, 11, 15, 21]
    out = {}
    for names, Js in ((('pp', 'pm', 'odd'), even_J), (('oddJ',), odd_J)):
        data = {Jv: get(Jv) for Jv in Js}
        for name in names:
            monos = sorted({rel(m, Jv) for Jv in Js for m in data[Jv][name]})
            vec = {}
            for mr in monos:
                pts = {}
                for Jv in Js:
                    m = mr[:2] + (Jv - mr[2],) + mr[3:]
                    pts[Jv] = data[Jv][name].get(m, sp.Integer(0))
                c, ok = fit(pts)
                vec[mr] = c
                print(f"{name:5s} {mr}: {c}   [{'out-of-sample OK' if ok else 'FAIL'}]")
                assert ok
            out[name] = vec
    pickle.dump(out, open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'wb'))
