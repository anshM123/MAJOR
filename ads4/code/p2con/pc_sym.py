"""
p2con/pc_sym.py (P2_CONTRACT.md §1.2): the exact heavy symbols of P2_SPINNING (s2_main conventions), made drift-free and converted to
the engine's rest-frame insertion monomials.

Input per q = e^{beta/2}: for each structure h in {same, opp, OO, odd} the left symbol orders o = 0, 1, 2 as {(k1, k2): poly(X, Y)}
(light-cone components xu, xv, xy, yu, yv, yy; P2_SPINNING §2), from p2spin/cache/s2_main.pkl (q = 11/10, 3/2, 2; read-only) or from
this directory's cache/sym_q*.pkl (pc_gensym.py, same code path as s2_main).

Drift removal (congruence, W1c): with U(eps) = 1 + eps U1 + eps^2 U2 acting on the structure index, chosen so that the drift-free
series  L~ = U^-1 L  has NO X,Y-independent component along the leading class tensors at orders 1 and 2:
    L~(1) = L(1) - U1 L(0),   L~(2) = L(2) - U1 L(1) - [const class part],   U1_{ht} = const class-t coefficient of L_h(1).
Class tensors: same = (2,2) + (-2,-2), opp = (2,-2) + (-2,2) (unit amplitudes, s2_main's helicity basis); OO, odd: their own leading.
Right symbol: R_h(eps) = L_h(-eps) at (t3, t4) (reflection identity, P2 §2.3).
Rest-frame monomials: X = B(beta/2) X_rest (side A, pair 2,4), Y = B(-beta/2) Y_rest (side B, pair 1,3); the plane wave is exp(2P.X + 2Pb.Y),
so a monomial of degree d carries 2^-d (CHMRS's (1/2 d_P)^d; validated in pc_engine checks).
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
P2C = os.path.join(HERE, '..', 'p2spin', 'cache')
import numpy as np
import sympy as sp

xs = sp.symbols('xu xv xy'); ys = sp.symbols('yu yv yy')
Xr = sp.symbols('xr0:3'); Yr = sp.symbols('yr0:3')
CLASSES = {'same': [(2, 2), (-2, -2)], 'opp': [(2, -2), (-2, 2)]}
EVEN = ('same', 'opp')


def _ev(s):
    return sp.sympify(eval(s, {**sp.__dict__}))


def load_q(q):
    """{h: [ {mode: poly}, ... orders 0..2 ]} at exact rational q."""
    q = q if isinstance(q, sp.Rational) else sp.nsimplify(q)
    fn = os.path.join(HERE, 'cache', 'sym_q%s.pkl' % str(q).replace('/', '_'))
    if os.path.exists(fn):
        d = pickle.load(open(fn, 'rb'))
        return {h: [{tuple(k): _ev(c) for k, c in dd.items()} for dd in comb] for h, comb in d.items()}
    d = pickle.load(open(os.path.join(P2C, 's2_main.pkl'), 'rb'))
    for qk, hb in d['orders'].items():
        if _ev(qk) == q:
            return {h: [{tuple(kk): _ev(c) for kk, c in dd.items()} for dd in comb] for h, comb in hb.items()}
    raise KeyError(q)


def _const(p):
    return sp.expand(p).subs({**{x: 0 for x in xs}, **{y: 0 for y in ys}})


def drift_free(sym):
    """returns {h: [L~0, L~1, L~2]} and the drift matrices (U1, U2 const parts) for the report."""
    out = {}; info = {}
    # even sector: 2x2 over (same, opp)
    L = {h: sym[h] for h in EVEN}
    ref = {t: L[t][0][CLASSES[t][0]] for t in EVEN}                     # leading amplitude (= 1 by construction)
    U1 = {(h, t): _const(L[h][1].get(CLASSES[t][0], 0))/ref[t] for h in EVEN for t in EVEN}
    asym1 = {(h, t): _const(L[h][1].get(CLASSES[t][0], 0) - L[h][1].get(CLASSES[t][1], 0)) for h in EVEN for t in EVEN}
    Lt = {}
    for h in EVEN:
        o1 = dict(L[h][1])
        for t in EVEN:
            for md in set(L[t][0]):
                o1[md] = sp.expand(o1.get(md, 0) - U1[(h, t)]*L[t][0][md])
        o2 = dict(L[h][2])
        for t in EVEN:
            for md in set(L[t][1]):
                o2[md] = sp.expand(o2.get(md, 0) - U1[(h, t)]*L[t][1][md])
        Lt[h] = [dict(L[h][0]), o1, o2]
    U2 = {(h, t): _const(Lt[h][2].get(CLASSES[t][0], 0))/ref[t] for h in EVEN for t in EVEN}
    asym2 = {(h, t): _const(Lt[h][2].get(CLASSES[t][0], 0) - Lt[h][2].get(CLASSES[t][1], 0)) for h in EVEN for t in EVEN}
    for h in EVEN:
        for t in EVEN:
            for md in set(L[t][0]):
                Lt[h][2][md] = sp.expand(Lt[h][2].get(md, 0) - U2[(h, t)]*L[t][0][md])
        out[h] = [{k: v for k, v in d.items() if v != 0} for d in Lt[h]]
    info['even'] = dict(U1=U1, U2=U2, asym1=asym1, asym2=asym2)
    # single-structure sectors
    for h in ('OO', 'odd'):
        if h not in sym:
            continue
        L0 = sym[h][0]
        md0 = max(L0, key=lambda k: (abs(complex(sp.N(L0[k]))), k))
        u1 = _const(sym[h][1].get(md0, 0))/L0[md0]
        o1 = {k: sp.expand(sym[h][1].get(k, 0) - u1*L0.get(k, 0)) for k in set(sym[h][1]) | set(L0)}
        o2 = {k: sp.expand(sym[h][2].get(k, 0) - u1*sym[h][1].get(k, 0)) for k in set(sym[h][2]) | set(sym[h][1])}
        u2 = _const(o2.get(md0, 0))/L0[md0]
        o2 = {k: sp.expand(o2.get(k, 0) - u2*L0.get(k, 0)) for k in set(o2) | set(L0)}
        out[h] = [dict(L0), {k: v for k, v in o1.items() if v != 0}, {k: v for k, v in o2.items() if v != 0}]
        info[h] = dict(u1=u1, u2=u2, ref=md0)
    return out, info


def rest_subs(q):
    q = q if isinstance(q, sp.Rational) else sp.nsimplify(q)
    ch = (q + 1/q)/2; sh = (q - 1/q)/2
    XG = [ch*Xr[0] + sh*Xr[1], sh*Xr[0] + ch*Xr[1], Xr[2]]
    YG = [ch*Yr[0] - sh*Yr[1], -sh*Yr[0] + ch*Yr[1], Yr[2]]
    return {xs[0]: XG[0] + XG[1], xs[1]: XG[0] - XG[1], xs[2]: XG[2],
            ys[0]: YG[0] + YG[1], ys[1]: YG[0] - YG[1], ys[2]: YG[2]}


def to_terms(poly, q, sub=None):
    """poly in (xu..yy) -> {(aX, aY): complex coefficient * 2^-deg} in rest-frame monomials."""
    sub = sub or rest_subs(q)
    p = sp.expand(sp.sympify(poly).subs(sub, simultaneous=True))
    if p == 0:
        return {}
    P = sp.Poly(p, *Xr, *Yr)
    return {(m[:3], m[3:]): complex(sp.N(c, 30))*2.0**(-sum(m)) for m, c in P.terms()}


def seed_terms(q):
    """the scalar seed's O(eps^2) insertion s2 (CHMRS T, Delta_phi = 3; p2spin/s3_mech.seed_series), rest-frame terms."""
    sys.path.insert(0, os.path.join(HERE, '..', 'p2spin'))
    import s3_mech as M
    q = q if isinstance(q, sp.Rational) else sp.nsimplify(q)
    P, Pb, (s1, s2) = M.seed_series(q)
    assert s1 == 0
    ch = (q + 1/q)/2; sh = (q - 1/q)/2
    XG = [ch*Xr[0] + sh*Xr[1], sh*Xr[0] + ch*Xr[1], Xr[2]]
    YG = [ch*Yr[0] - sh*Yr[1], -sh*Yr[0] + ch*Yr[1], Yr[2]]
    p = sp.expand(s2.subs({**{M.X[i]: XG[i] for i in range(3)}, **{M.Y[i]: YG[i] for i in range(3)}}, simultaneous=True))
    Pl = sp.Poly(p, *Xr, *Yr)
    return {(m[:3], m[3:]): complex(sp.N(c, 30))*2.0**(-sum(m)) for m, c in Pl.terms()}


def numeric_series(Lt, q):
    """{h: [ {mode: {(aX, aY): c}} for orders ]}"""
    sub = rest_subs(q)
    return {h: [{md: to_terms(p, q, sub) for md, p in d.items()} for d in ser] for h, ser in Lt.items()}


def mul_terms(t1, t2, c=1.0):
    out = {}
    for (a1, b1), v1 in t1.items():
        for (a2, b2), v2 in t2.items():
            k = (tuple(x + y for x, y in zip(a1, a2)), tuple(x + y for x, y in zip(b1, b2)))
            out[k] = out.get(k, 0) + c*v1*v2
    return out


def block_tensor(S, h, hp, n, seed=None):
    """O(eps^n) tensor of entry (h, h'): {(k1, k2, k3, k4): {(aX, aY): c}} =
       sum_o L~_h^(o) (x) (-1)^(n-o) L~_h'^(n-o)  [+ s2 L~_h^(0) (x) L~_h'^(0) if n = 2 and seed]."""
    out = {}
    for o in range(n + 1):
        sg = (-1)**(n - o)
        for m1, t1 in S[h][o].items():
            for m2, t2 in S[hp][n - o].items():
                key = (m1[0], m1[1], m2[0], m2[1])
                pr = mul_terms(t1, t2, sg)
                d = out.setdefault(key, {})
                for k, v in pr.items():
                    d[k] = d.get(k, 0) + v
    if n == 2 and seed is not None:
        for m1, t1 in S[h][0].items():
            for m2, t2 in S[hp][0].items():
                key = (m1[0], m1[1], m2[0], m2[1])
                pr = mul_terms(mul_terms(t1, t2), seed)
                d = out.setdefault(key, {})
                for k, v in pr.items():
                    d[k] = d.get(k, 0) + v
    return {k: {kk: vv for kk, vv in v.items() if abs(vv) > 0} for k, v in out.items()}
