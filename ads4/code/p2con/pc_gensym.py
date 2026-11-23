"""
p2con/pc_gensym.py (P2_CONTRACT.md §1.2): exact heavy symbols of the conserved <TTO> structures to O(1/m^2) at a given rational
q = e^{beta/2}, by EXACTLY the code path of p2spin/s2_main.py (imported read-only: hel_pol, modes_from_pol, sym_modes, the helicity
basis HB of p2spin/cache/s2_main.pkl), with exact polynomial interpolation in (X, Y) (degree <= 2, 28 monomials, residual checked on
4 extra points).  Output: cache/sym_q<p>_<q>.pkl  {h: [ {mode: srepr(poly)} for orders 0..2 ]}, h in same, opp, OO, odd.
Usage: python3 pc_gensym.py <beta1> <beta2> ...   (q = rational approximation of e^{beta/2}, |error| < 1e-7)
"""
import os, sys, pickle, time, random
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'p2spin'))
import sympy as sp
import p2_series as PS
import s2_main as SM          # read-only import (its main() is not run)
import s2_symbols as S2

FAMS = [('even', 1, 0), ('even', 1, 1), ('odd', 1, 0), ('odd', -1, 0)]


def qrat(beta):
    from fractions import Fraction
    import math
    f = Fraction(math.exp(beta/2)).limit_denominator(20000)
    return sp.Rational(f.numerator, f.denominator)


def gen(qv):
    t0 = time.time()
    can = S2.CAN()
    d = pickle.load(open(os.path.join(HERE, '..', 'p2spin', 'cache', 's2_main.pkl'), 'rb'))
    ev = lambda s: sp.sympify(eval(s, {**sp.__dict__}))
    qs = sp.Symbol('q', positive=True)
    HB = {h: tuple(ev(c) for c in v) for h, v in d['HB'].items()}
    PS.ORDER = 2
    X = S2.xs; Y = S2.ys
    mons = [()] + [(i,) for i in range(6)] + [(i, j) for i in range(6) for j in range(i, 6)]
    def monval(pt, mon):
        r = sp.Integer(1)
        for i in mon:
            r *= pt[i]
        return r
    rng = random.Random(42)
    def rpt():
        while True:
            u1, v1, u2, v2 = [sp.Rational(rng.randint(1, 30), rng.randint(1, 9)) for _ in range(4)]
            y1 = sp.Rational(rng.randint(-9, 9), 10)*sp.sqrt(u1*v1); y2 = sp.Rational(rng.randint(-9, 9), 10)*sp.sqrt(u2*v2)
            y1 = sp.nsimplify(sp.floor(y1*7)/7); y2 = sp.nsimplify(sp.floor(y2*7)/7)
            if y1**2 < u1*v1 and y2**2 < u2*v2:
                return (u1, v1, y1, u2, v2, y2)
    NPT = len(mons) + 4
    pts = [rpt() for _ in range(NPT)]
    Minv = sp.Matrix([[monval(pt, mo) for mo in mons] for pt in pts[:len(mons)]]).inv()
    allv = list(X) + list(Y)
    maxres = 0
    def interp(vals):
        nonlocal maxres
        keys = set().union(*[set(v) for v in vals])
        out = {}
        for kk in keys:
            c = Minv*sp.Matrix([vals[i].get(kk, 0) for i in range(len(mons))])
            poly = sp.expand(sum(ci*monval(allv, mo) for ci, mo in zip(c, mons)))
            for i in range(len(mons), NPT):
                r = sp.simplify(poly.subs(dict(zip(allv, pts[i]))) - vals[i].get(kk, 0))
                maxres = max(maxres, abs(complex(sp.N(r))))
            if poly != 0:
                out[kk] = poly
        return out
    sm = {}
    for sec, sig, k in FAMS:
        per = [[], [], []]
        for pt in pts:
            v, ms = SM.sym_modes('L', (sec, sig), can[(sec, sig)][k], pt[:3], pt[3:], qv)
            assert v == 0
            for o in range(3):
                per[o].append(ms[o])
        sm[(sec, sig, k)] = [interp(per[o]) for o in range(3)]
    hb = {}
    for h in ('same', 'opp'):
        c0, c1 = [c.subs(qs, qv) for c in HB[h]]
        comb = []
        for o in range(3):
            dd = {}
            for kk in set(sm[('even', 1, 0)][o]) | set(sm[('even', 1, 1)][o]):
                dd[kk] = sp.expand(c0*sm[('even', 1, 0)][o].get(kk, 0) + c1*sm[('even', 1, 1)][o].get(kk, 0))
            comb.append({kk: c for kk, c in dd.items() if c != 0})
        hb[h] = comb
    hb['OO'] = sm[('odd', 1, 0)]; hb['odd'] = sm[('odd', -1, 0)]
    fn = os.path.join(HERE, 'cache', 'sym_q%s.pkl' % str(qv).replace('/', '_'))
    pickle.dump({h: [{kk: sp.srepr(c) for kk, c in dd.items()} for dd in comb] for h, comb in hb.items()}, open(fn, 'wb'))
    return time.time() - t0, maxres


if __name__ == '__main__':
    for b in sys.argv[1:]:
        qv = qrat(float(b)) if '/' not in b else sp.Rational(b)
        fn = os.path.join(HERE, 'cache', 'sym_q%s.pkl' % str(qv).replace('/', '_'))
        if os.path.exists(fn):
            print(f'beta={b}: q={qv} cached', flush=True); continue
        dt, mr = gen(qv)
        print(f'beta={b}: q={qv} (beta_q = {float(2*sp.log(qv)):.8f}) done in {dt:.0f}s, interpolation residual {mr:.1e}', flush=True)
