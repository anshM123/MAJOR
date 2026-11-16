"""
Three-point structures <O_1(x) O_2(0) X(inf, zeta)> in the frame at infinity, for the differential-basis words of fourpt/structures2
(l = 1: 5 words; l = 2: 14 words, exact expressions of fourpt/build_l2.py).  Euclidean.
Embedding at infinity: P_inf = (0,1,0), Z_inf = (0,0,zeta)  =>  s_{i5} = -1, Z_i.P_inf = 0, Z_inf.P_i = zeta.x_i;  x_1 = x, x_2 = 0:
   s12 -> -x^2, s15, s25 -> -1, z12 -> -z1.x, z21 -> z2.x, z15, z25 -> 0, z51 -> zeta.x, z52 -> 0, z1z2, z1z5 -> z1.zeta, z2z5 -> z2.zeta,
   M_ref -> v5^J |x^2|^{-al},  v5 = -zeta.x/2,  al = (2(l+1) + 2l - Delta + J)/2  (seed S3(l+1+l, l+1+l)).
Structure of word w:  F_w = (-zeta.x/2)^J (x^2)^{-al} R_w(x, z1, z2, zeta),  R_w rational (powers of 1/(zeta.x), 1/x^2).
"""
import os, sys, pickle, functools
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from embed2 import Delta, J
from structures2 import words

X = sp.symbols('x0:3'); Z1 = sp.symbols('a0:3'); Z2 = sp.symbols('b0:3'); W = sp.symbols('w0:3')   # W = zeta (symbolic)
def dot(u, v): return sum(u[i]*v[i] for i in range(3))
CACHE = os.path.join(HERE, 'cache')

@functools.lru_cache(None)
def word_exprs(l):
    if l == 1:
        from structures2 import rel_structures
        rs, S, v5 = rel_structures(1)
        return {w: Q for (w, Q) in rs}, S
    rel, S, v5 = pickle.load(open(os.path.join(HERE, '..', 'cache', 'rel_l2.pkl'), 'rb'))
    return rel, S

def inf_sub(S):
    return {S['s12']: -dot(X, X), S['s15']: -1, S['s25']: -1, S['z12']: -dot(Z1, X), S['z21']: dot(Z2, X), S['z15']: 0, S['z25']: 0,
            S['z51']: dot(W, X), S['z52']: 0, S['z1z2']: dot(Z1, Z2), S['z1z5']: dot(Z1, W), S['z2z5']: dot(Z2, W)}

def alpha_ref(l):
    return (4*l + 2 - Delta + J)/2

@functools.lru_cache(None)
def R_words(l):
    fn = os.path.join(CACHE, 'Rinf_l%d.pkl' % l)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    ex, S = word_exprs(l)
    # first in invariants (small), then in components
    AX, BX, WX, XX, AB, AW, BW = sp.symbols('AX BX WX XX AB AW BW')
    isub = {S['s12']: -XX, S['s15']: -1, S['s25']: -1, S['z12']: -AX, S['z21']: BX, S['z15']: 0, S['z25']: 0,
            S['z51']: WX, S['z52']: 0, S['z1z2']: AB, S['z1z5']: AW, S['z2z5']: BW}
    csub = {AX: dot(Z1, X), BX: dot(Z2, X), WX: dot(W, X), XX: dot(X, X), AB: dot(Z1, Z2), AW: dot(Z1, W), BW: dot(Z2, W)}
    out = {}
    for w in words(l):
        r = sp.together(ex[w].xreplace(isub))
        n, d = sp.fraction(r)
        out[w] = sp.expand(n).xreplace(csub)/d.xreplace(csub)
    pickle.dump(out, open(fn, 'wb'))
    return out
