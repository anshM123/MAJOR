"""
t1a_boundary (DSKETCH_UP.md §1.5): near-boundary fall-off of the Weyl tensor on the explicit configuration
g = z^-2 (delta + a gam3 + b gam4 + a b Gam_p) of t1a_solve2 (the C1_VALUE §3.1 premise that I_+- = int sqrt(g) Tr W_+-^3 needs no
boundary counterterm: orthonormal-frame Weyl = O(z^3)).

Orthonormal-frame components W_{(abcd)} = z^4 W_{mu nu rho sigma}[g] = z^2 W~_{mu nu rho sigma}[g~] (W^a_{bcd} conformally invariant),
so the premise is W~ = O(z) at every order a, b, ab, with W~ the Weyl tensor of g~ = delta + ..., lowered with g~.
Checks: every component of W~ at orders (1,0), (0,1), (1,1) vanishes at z = 0, and its O(z) coefficient is non-zero for some
component (so O(z) is sharp: orthonormal Weyl = O(z^3), magnetic/electric data at z^3).
Usage: python3 -B t1a_boundary.py <tag of t1a_solve2>
"""
import os, sys, pickle, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from sympy import *  # noqa: F401,F403
from asd_alg import Alg, z, N, christoffel, riemann_up, ricci, inverse_conformal, mode

tag = sys.argv[1] if len(sys.argv) > 1 else '0_0_3_4_0_0'
kk = [int(v) for v in tag.split('_')]
k3, k4 = kk[:3], kk[3:]
GamP = eval(pickle.load(open(os.path.join(HERE, 'cache', f't1a_solve2_{tag}.pkl'), 'rb'))['GamP'])
A = Alg(k3, k4)
g3 = mode(k3, -1); g4 = mode(k4, -1)
H = [[{(1, 0): g3[i, j], (0, 1): g4[i, j], (1, 1): GamP[i, j]} for j in range(N)] for i in range(N)]
g = [[A.add({(0, 0): sp.Integer(1 if i == j else 0)}, H[i][j]) for j in range(N)] for i in range(N)]
gi = inverse_conformal(A, H)
G = christoffel(A, g, gi)
Rup = riemann_up(A, G)
Ric = ricci(A, Rup)
Rsc = {}
for b_, dd in itertools.product(range(N), repeat=2):
    Rsc = A.add(Rsc, A.mul(gi[b_][dd], Ric[b_][dd]))
S = [[A.scale(sp.Rational(1, 2), A.add(Ric[i][j], A.scale(-sp.Rational(1, 6), A.mul(Rsc, g[i][j])))) for j in range(N)] for i in range(N)]
lines = []
def log(x=''):
    print(x, flush=True); lines.append(str(x))
ok = True
log(f't1a_boundary: tag {tag}')
for o in [(1, 0), (0, 1), (1, 1)]:
    at0 = 0; lin = 0; nz = 0
    for a_, b_, c_, dd in itertools.product(range(N), repeat=4):
        Rl = {}
        for e in range(N):
            Rl = A.add(Rl, A.mul(g[a_][e], Rup[(e, b_, c_, dd)]))
        W = A.add(Rl, A.scale(-1, A.mul(g[a_][c_], S[b_][dd])), A.mul(g[a_][dd], S[b_][c_]),
                  A.scale(-1, A.mul(g[b_][dd], S[a_][c_])), A.mul(g[b_][c_], S[a_][dd]))
        w = W.get(o, 0)
        if w == 0:
            continue
        ser = sp.series(sp.simplify(w), z, 0, 2).removeO()
        c0 = sp.simplify(ser.coeff(z, 0)); c1 = sp.simplify(ser.coeff(z, 1))
        nz += 1
        if c0 != 0:
            at0 += 1
        if c1 != 0:
            lin += 1
    good = at0 == 0 and lin > 0
    log(f'  order {o}: components {nz:3d}; non-zero at z = 0: {at0}; with O(z) term: {lin}   => W~ = O(z) sharp: {good}')
    ok &= good
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', f't1a_boundary_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')
