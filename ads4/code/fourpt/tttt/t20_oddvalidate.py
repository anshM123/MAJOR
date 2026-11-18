"""
Validation of taylor4odd.Engine4Odd against a brute-force application of the same operators as explicit differential operators on
5-component embedding vectors (sympy), for a test seed  F = Kt(s) exp(c1 P + c2 Q)  (P, Q as in taylor4.seed), at a random numerical point
on the constraint surface (P_i^2 = Z_i.P_i = Z_i^2 = 0; metric of the Lorentzian frame, eta_{+-} = -1/2, diag(-1, 1, 1)).
The epsilon products are evaluated with the Gram identity of taylor4odd (det(eta_E) = -1/4) on the engine side and as explicit Levi-Civita
contractions with det(eta_L) on the brute-force side; the two must agree up to the constant det(eta_L)/det(eta_E) per epsilon PAIR
(= -1: the Lorentzian continuation of the Euclidean-normalised structure), which is checked.
Usage: python3 t20_oddvalidate.py   -> out_t20_oddvalidate.txt
"""
import os, sys, itertools, random
import numpy as np, sympy as sp
from sympy.combinatorics import Permutation
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
from taylor4 import seed, tables, SN, ZN
from taylor4odd import Engine4Odd, word_budget_odd, seed_shifts, gram_value
from embed2 import shifts as shifts_even

ETA = sp.zeros(5, 5); ETA[0, 1] = ETA[1, 0] = -sp.Rational(1, 2); ETA[2, 2] = -1; ETA[3, 3] = 1; ETA[4, 4] = 1
EIN = ETA.inv()
VN = ['P1', 'P2', 'P3', 'P4', 'Z1', 'Z2', 'Z3', 'Z4']
SYM = {n: sp.Matrix(sp.symbols(n + '_0:5')) for n in VN}
def d(u, v): return (u.T*ETA*v)[0]
LC = {}
for p in itertools.permutations(range(5)):
    LC[p] = Permutation(list(p)).signature()

def eps_vec(vs):
    """covariant-index-free: eps(v1..v4, .)_E as a contravariant vector (index raised with EIN), vs: 4 vectors"""
    out = sp.zeros(5, 1)
    for p, sg in LC.items():
        coef = sg*vs[0][p[0]]*vs[1][p[1]]*vs[2][p[2]]*vs[3][p[3]]
        if coef != 0:
            for F in range(5):
                out[F] += coef*EIN[p[4], F]
    return out

def eps5(vs):
    return sp.Matrix([list(v) for v in vs]).det()

def grad(F, name):
    return sp.Matrix([sp.diff(F, SYM[name][A]) for A in range(5)])

def apply_bf(F, kind, i, j):
    Pi, Pj, Zi, Zj = SYM['P%d' % i], SYM['P%d' % j], SYM['Z%d' % i], SYM['Z%d' % j]
    nPi, nPj, nZi, nZj = 'P%d' % i, 'P%d' % j, 'Z%d' % i, 'Z%d' % j
    dd = lambda V, W: (V.T*grad(F, W))[0]          # V.d/dW
    if kind == 'H':
        return F*(d(Zi, Zj)*d(Pi, Pj) - d(Zi, Pj)*d(Zj, Pi))
    if kind in 'AC':
        if kind == 'C': Pi, Pj, Zi, Zj, nPi, nPj, nZi, nZj = Pj, Pi, Zj, Zi, nPj, nPi, nZj, nZi
        return d(Pi, Pj)*dd(Zi, nPj) - d(Zi, Pj)*dd(Pi, nPj) - d(Zi, Zj)*dd(Pi, nZj) + d(Pi, Zj)*dd(Zi, nZj)
    if kind in 'BE':
        if kind == 'E': Pi, Pj, Zi, Zj, nPi, nPj, nZi, nZj = Pj, Pi, Zj, Zi, nPj, nPi, nZj, nZi
        return d(Pi, Pj)*dd(Zi, nPi) - d(Zi, Pj)*dd(Pi, nPi) + d(Zi, Pj)*dd(Zi, nZi)
    if kind in 'OQ':
        W = nPi if kind == 'O' else nPj
        e = eps_vec([Zi, Zj, Pi, Pj])
        return (e.T*grad(F, W))[0]
    if kind in 'UV':
        if kind == 'U': e = eps_vec([Zi, Pi, Pj, Zj]); W = nZj
        else: e = eps_vec([Zj, Pj, Pi, Zi]); W = nZi
        return (e.T*grad(F, W))[0]
    if kind in 'RS':
        Zs = Zi if kind == 'R' else Zj
        # eps(Zs, Pi, Pj, d/dPi, d/dPj) F = sum_{F,G} eps(Zs,Pi,Pj,e_F',e_G') EIN EIN dF/dPi^F' dF/dPj^G'
        out = 0
        g1 = grad(F, nPi)
        for A in range(5):
            gA = sp.diff(g1[A], SYM[nPj]) if False else None
        H = sp.Matrix(5, 5, lambda A, B: sp.diff(F, SYM[nPi][A], SYM[nPj][B]))
        for p, sg in LC.items():
            c = sg*Zs[p[0]]*Pi[p[1]]*Pj[p[2]]
            if c == 0: continue
            for A in range(5):
                for B in range(5):
                    w = EIN[p[3], A]*EIN[p[4], B]
                    if w != 0: out += c*w*H[A, B]
        return out
    raise ValueError(kind)

def test(wL, wR, cE=(sp.Rational(3, 10), sp.Rational(-1, 7)), seedd=None):
    random.seed(4)
    # numeric point on the constraint surface: P_i = (1, x^2, x), Z_i = (0, 2 x.z, z), z null (Lorentzian x.y = -x0y0 + x1y1 + x2y2)
    def mdot(a, b): return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    pts = {}
    xs = {1: [sp.Rational(0), sp.Rational(-1, 13), sp.Rational(2, 17)], 2: [sp.Rational(1), sp.Rational(1, 9), sp.Rational(-1, 11)],
          3: [sp.Rational(23, 10), sp.Rational(1, 7), sp.Rational(1, 8)], 4: [sp.Rational(37, 10), sp.Rational(-1, 5), sp.Rational(1, 6)]}
    zs = {1: [5, 3, 4], 2: [13, 5, 12], 3: [5, 4, -3], 4: [17, 8, 15]}             # null in (-,+,+)
    for i in (1, 2, 3, 4):
        x, z = xs[i], zs[i]
        pts.update({SYM['P%d' % i][A]: v for A, v in enumerate([1, mdot(x, x)] + x)})
        pts.update({SYM['Z%d' % i][A]: v for A, v in enumerate([0, 2*mdot(x, z)] + z)})
    sL = shifts_even(wL) if not any(c in wL for c in 'OQRSUV') else seed_shifts(wL)
    sR = shifts_even(wR) if not any(c in wR for c in 'OQRSUV') else seed_shifts(wR)
    D = (3 + sL[0], 3 + sL[1], 3 + sR[0], 3 + sR[1])
    s = {n: 2*d(SYM['P' + n[1]], SYM['P' + n[2]]) for n in SN}
    lu = sp.log(s['s12']) + sp.log(s['s34']) - sp.log(s['s13']) - sp.log(s['s24'])
    lv = sp.log(s['s14']) + sp.log(s['s23']) - sp.log(s['s13']) - sp.log(s['s24'])
    P = (lu + lv)/2; Q = (lv - lu)/2
    d1, d2, d3, d4 = D
    Kt = (s['s24']/s['s14'])**sp.Rational(d1 - d2, 2)*(s['s24']/s['s23'])**sp.Rational(d3 - d4, 2)/(s['s12']**sp.Rational(d1 + d2, 2)*s['s34']**sp.Rational(d3 + d4, 2))
    F = Kt*sp.exp(cE[0]*P + cE[1]*Q)
    for k in reversed(wL): F = apply_bf(F, k, 1, 2)
    for k in reversed(wR): F = apply_bf(F, k, 3, 4)
    bf = complex(sp.N(F.subs(pts), 30))
    # engine
    M = word_budget_odd(wL, wR)
    s0 = {n: np.array([float(sp.N(s[n].subs(pts)))]) for n in SN}
    gt0 = 1.0
    G = [[np.array([float(cE[0])**i*float(cE[1])**j]) for j in range(M + 1)] for i in range(M + 1)]
    # gt = exp(c1 P + c2 Q) at the point: G_ij = c1^i c2^j gt(point)
    gtv = float(sp.N(sp.exp(cE[0]*P + cE[1]*Q).subs(pts)))
    G = [[g*gtv for g in row] for row in G]
    sd = seed(M, s0, tuple(float(x) for x in D), G)[(0,)*18]
    E = Engine4Odd(M, s0)
    res = E.apply_words({(0,)*18 + ((),): sd[:tables(M)['nb'][M]]}, wL, wR)
    zval = {}
    for n in ZN:
        if n.count('z') == 2:
            zval[n] = float(sp.N(d(SYM['Z' + n[1]], SYM['Z' + n[3]]).subs(pts)))
        else:
            zval[n] = float(sp.N(d(SYM['Z' + n[1]], SYM['P' + n[2]]).subs(pts)))
    def dots(u, v):
        return float(sp.N(d(SYM[u], SYM[v]).subs(pts)))
    def embed(name):
        return np.array([float(sp.N(SYM[name][A].subs(pts))) for A in range(5)])
    tot = 0
    nps = 0
    for key, v in res.items():
        mono, eps = key[:18], key[18]
        val = v[0]
        for n, e in zip(ZN, mono):
            if e: val *= zval[n]**e
        if eps:
            val *= gram_value(eps, dots, embed); nps = len(eps)
        tot += val
    return bf, tot, nps

if __name__ == "__main__":
    out = []
    pairs = [tuple(p.split(':')) for p in sys.argv[1:]] or [('O', 'O'), ('Q', 'O'), ('HO', 'HQ'), ('AC', 'HO'), ('S', 'O'), ('R', 'Q'), ('AS', 'O'), ('S', 'S'), ('U', 'O')]
    for wL, wR in pairs:
        try:
            bf, en, nps = test(wL, wR)
            r = en/bf if bf != 0 else float('nan')
            s = f"  {wL:>5s} x {wR:<5s}: brute force {bf:.10e}, engine {en:.10e}, ratio {r:.12f}  (eps factors {nps}; expected ratio {(-1)**(nps//2) if nps == 2 else 1})"
        except Exception as ex:
            s = f"  {wL} x {wR}: error {ex!r}"
        print(s, flush=True); out.append(s)
    open(os.path.join(HERE, 'out_t20_oddvalidate.txt'), 'a').write("\n".join(out) + "\n")
