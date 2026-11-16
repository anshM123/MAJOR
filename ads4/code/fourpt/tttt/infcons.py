"""
Conservation of <O_1(x) O_2(0) X(inf)> at x_1 and x_2 in the frame at infinity (fock.py convention: divergence in x),
symbolic in (Delta, J).  F_w = (zeta.x)^J (x^2)^{-al} R_w;  zeta = (1, i, 0).
Returns the null space (coefficient vectors on the words), modulo the d = 3 relations (which vanish identically here).
"""
import os, sys, pickle, random
import sympy as sp
from infframe import X, Z1, Z2, W, dot, R_words, alpha_ref, CACHE, Delta, J, words
ZETA = (1, sp.I, 0)

def harm(R, Z, l):
    if l < 2: return R
    lap = sum(sp.diff(R, Z[i], 2) for i in range(3))
    return R - dot(Z, Z)*lap/6

def div_eqs(R, l, pt):
    """numerator polynomial of  G^{-1} div_pt (G R)  with G = (zeta.x)^J (x^2)^{-al}"""
    Z = Z1 if pt == 1 else Z2
    sgn = 1 if pt == 1 else -1
    al = alpha_ref(l)
    zx = dot(W, X); xx = dot(X, X)
    h = harm(R, Z, l)
    d = 0
    for m in range(3):
        dh = sp.diff(h, Z[m])
        d += (J*W[m]/zx - 2*al*X[m]/xx)*dh + sp.diff(dh, X[m])
    return sgn*d

def null_space(l, Dv=None, Jv=None, npts=None):
    R = R_words(l)
    ws = words(l)
    zs = {W[i]: ZETA[i] for i in range(3)}
    num = {}
    if Dv is not None: num[Delta] = Dv
    if Jv is not None: num[J] = Jv
    rows = []
    random.seed(3)
    cols = []
    for w in ws:
        e = []
        for pt in (1, 2):
            d = sp.together(div_eqs(R[w].subs(num), l, pt).subs(zs))
            e.append(d)
        cols.append(e)
    # common denominator per point: evaluate at random rational x, z (exact) -> linear equations
    M = []
    for k in range(6 if l == 1 else 14):
        p = {**{X[i]: sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for i in range(3)},
             **{Z1[i]: sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for i in range(3)},
             **{Z2[i]: sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for i in range(3)}}
        for pt in (0, 1):
            M.append([sp.expand(sp.together(cols[j][pt].subs(p))) for j in range(len(ws))])
    Mm = sp.Matrix(M)
    ns = Mm.nullspace(simplify=True)
    out = []
    for v in ns:
        den = sp.lcm([sp.fraction(sp.together(x))[1] for x in v])
        out.append([sp.factor(sp.cancel(x*den)) for x in v])
    return ws, out
