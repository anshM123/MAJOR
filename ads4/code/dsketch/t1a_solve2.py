"""
t1a_solve2 (DSKETCH_UP.md §1.2-1.3): explicit regular second-order solution with vanishing second-order SD Weyl tensor.

1. Solve the O(ab) Einstein equations of t1a_einstein2 with the ansatz F_ij = P_ij(z) exp(-(|k3|+|k4|) z), P_ij polynomials
   (complex rational coefficients, degree <= DEG): a linear system, solved exactly.  Every free parameter left is set to 0
   (they are pure-gauge / homogeneous pieces of the same exponential class).
2. With g~ = delta + a gam3 + b gam4 + a b Gam_p (conformal frame; the Weyl tensor W^a_{bcd} is conformally invariant and
   the Hodge star on 2-forms is conformally invariant in 4D), compute the full Weyl tensor to O(ab), its SD projection with the
   FULL metric, P+[g] = (1 + *_g)/2, and check: W+^{(a)} = W+^{(b)} = 0 (ASD legs), W+^{(ab)} = 0 (the second-order
   solution is ASD), W-^{(ab)} != 0 (non-trivial), and the full Einstein equation at O(ab) on the solution.
3. Control: the same with the gauge-equivalent but non-ASD choice Gam_p + c * (regular label-+ mode at K): W+^{(ab)} != 0.
Writes cache/t1a_solve2_<k>.txt / .pkl.
"""
import os, sys, pickle, time, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from sympy import *  # noqa: F401,F403  (for eval of srepr strings)
from asd_alg import Alg, z, N, christoffel, riemann_up, ricci, inverse_conformal, mode

DEG = int(os.environ.get('DEG', '5'))
tag = sys.argv[1] if len(sys.argv) > 1 else '0_0_3_4_0_0'
dat = pickle.load(open(os.path.join(HERE, 'cache', f't1a_einstein2_{tag}.pkl'), 'rb'))
k3 = [sp.sympify(v) for v in dat['k3']]; k4 = [sp.sympify(v) for v in dat['k4']]
A = Alg(k3, k4)
n3 = sp.sqrt(sum(v**2 for v in k3)); n4 = sp.sqrt(sum(v**2 for v in k4))
Kv = [k3[i] + k4[i] for i in range(3)]
s = n3 + n4
lines = []
def log(x=''):
    print(x, flush=True); lines.append(str(x))

log(f't1a_solve2: k3 = {k3} (|k3| = {n3}), k4 = {k4} (|k4| = {n4}), K = {Kv} (|K| = {sp.sqrt(sum(v**2 for v in Kv))}); legs: regular label-(-) (ASD)')
eqs = {key: eval(v) for key, v in dat['eqs'].items()}       # already multiplied by exp(s z) z^2 in t1a_einstein2
# polynomial ansatz
coef = {}
Fsub = {}
for i in range(3):
    for j in range(i, 3):
        cs = sp.symbols(f'p{i}{j}_0:{DEG + 1}')
        coef[(i, j)] = cs
        Fsub[sp.Function(f'F{i}{j}')(z)] = sum(cs[n]*z**n for n in range(DEG + 1))*sp.exp(-s*z)
allc = [c for key in coef for c in coef[key]]
lin = []
for key, e in eqs.items():
    ee = sp.expand(sp.powsimp(sp.expand(e.subs(Fsub).doit())))
    P = sp.Poly(ee, z)
    lin += P.coeffs()
sol = sp.solve(lin, allc, dict=True)
log(f'polynomial ansatz degree {DEG}: {len(lin)} linear equations, {len(allc)} unknowns, solutions: {len(sol)}')
ok = len(sol) == 1
sol = sol[0]
free = [c for c in allc if c not in sol]
log(f'free parameters (set to 0): {len(free)}')
GamP = sp.zeros(4, 4)
for (i, j), cs in coef.items():
    p = sum(sp.sympify(sol.get(cs[n], cs[n])).subs({f: 0 for f in free})*z**n for n in range(DEG + 1))
    GamP[i, j] = GamP[j, i] = sp.expand(p)*sp.exp(-s*z)
    log(f'  Gam_p[{i}{j}] = ({sp.factor(sp.expand(p))}) * exp(-{s} z)')
# residual check of the Einstein equations at O(ab)
res_ok = True
for key, e in eqs.items():
    sub = {sp.Function(f'F{i}{j}')(z): GamP[i, j] for i in range(3) for j in range(i, 3)}
    r = sp.simplify(e.subs(sub).doit())
    if r != 0:
        res_ok = False; log(f'  residual {key}: {r}')
log(f'O(ab) Einstein equations solved exactly: {res_ok}')
ok &= res_ok


def weyl_sd(Gam):
    """W+ and W- of g~ = delta + a gam3 + b gam4 + a b Gam (conformal frame), projected with the full metric; returns dicts."""
    g3 = mode(k3, -1); g4 = mode(k4, -1)
    H = [[{(1, 0): g3[i, j], (0, 1): g4[i, j], (1, 1): Gam[i, j]} for j in range(N)] for i in range(N)]
    g = [[A.add({(0, 0): sp.Integer(1 if i == j else 0)}, H[i][j]) for j in range(N)] for i in range(N)]
    gi = inverse_conformal(A, H)
    G = christoffel(A, g, gi)
    Rup = riemann_up(A, G)
    Ric = ricci(A, Rup)
    # scalar curvature R = g^{bd} Ric_bd
    Rsc = {}
    for b_, dd in itertools.product(range(N), repeat=2):
        Rsc = A.add(Rsc, A.mul(gi[b_][dd], Ric[b_][dd]))
    # lower Riemann: R_{abcd} = g_{ae} R^e_{bcd}
    Rl = {}
    for a_, b_, c_, dd in itertools.product(range(N), repeat=4):
        acc = {}
        for e in range(N):
            acc = A.add(acc, A.mul(g[a_][e], Rup[(e, b_, c_, dd)]))
        Rl[(a_, b_, c_, dd)] = acc
    # Weyl (n = 4): W = R - (g_ac S_bd - g_ad S_bc + g_bd S_ac - g_bc S_ad), S = (Ric - R g/6)/2 (Schouten)
    S = [[A.scale(sp.Rational(1, 2), A.add(Ric[i][j], A.scale(-sp.Rational(1, 6), A.mul(Rsc, g[i][j])))) for j in range(N)] for i in range(N)]
    W = {}
    for a_, b_, c_, dd in itertools.product(range(N), repeat=4):
        W[(a_, b_, c_, dd)] = A.add(Rl[(a_, b_, c_, dd)],
                                     A.scale(-1, A.mul(g[a_][c_], S[b_][dd])), A.mul(g[a_][dd], S[b_][c_]),
                                     A.scale(-1, A.mul(g[b_][dd], S[a_][c_])), A.mul(g[b_][c_], S[a_][dd]))
    # Hodge dual on the first pair with the full metric: (*W)_{abcd} = 1/2 sqrt(g) eps_{abef} g^{ee'} g^{ff'} W_{e'f'cd}
    # sqrt(det g~) to O(ab): exp(tr log(1+H)/2) -> 1 + tr H/2 + [ (trH)^2/8 - tr(H.H)/4 ] at (1,1)
    trH = A.add(*[H[i][i] for i in range(N)])
    HH = {}
    for i, j in itertools.product(range(N), repeat=2):
        HH = A.add(HH, A.mul(H[i][j], H[j][i]))
    sq = A.add({(0, 0): sp.Integer(1)}, A.scale(sp.Rational(1, 2), trH), A.scale(sp.Rational(1, 8), A.mul(trH, trH)),
               A.scale(-sp.Rational(1, 4), HH))
    Wup = {}
    for e, f, c_, dd in itertools.product(range(N), repeat=4):
        acc = {}
        for e2, f2 in itertools.product(range(N), repeat=2):
            acc = A.add(acc, A.mul(A.mul(gi[e][e2], gi[f][f2]), W[(e2, f2, c_, dd)]))
        Wup[(e, f, c_, dd)] = acc
    Wp, Wm = {}, {}
    for a_, b_, c_, dd in itertools.product(range(N), repeat=4):
        acc = {}
        for e, f in itertools.product(range(N), repeat=2):
            lc = sp.LeviCivita(a_, b_, e, f)
            if lc != 0:
                acc = A.add(acc, A.scale(sp.Rational(1, 2)*lc, Wup[(e, f, c_, dd)]))
        star = A.mul(sq, acc)
        Wp[(a_, b_, c_, dd)] = A.scale(sp.Rational(1, 2), A.add(W[(a_, b_, c_, dd)], star))
        Wm[(a_, b_, c_, dd)] = A.scale(sp.Rational(1, 2), A.add(W[(a_, b_, c_, dd)], A.scale(-1, star)))
    return Wp, Wm


def report(Wp, Wm, label):
    res = {}
    for o in [(1, 0), (0, 1), (1, 1)]:
        nzp = [key for key in Wp if sp.simplify(Wp[key].get(o, 0)) != 0]
        nzm = [key for key in Wm if sp.simplify(Wm[key].get(o, 0)) != 0]
        res[o] = (len(nzp), len(nzm))
        log(f'  [{label}] order {o}: nonzero components  W+ : {len(nzp):3d}   W- : {len(nzm):3d}')
    return res


t0 = time.time()
Wp, Wm = weyl_sd(GamP)
r = report(Wp, Wm, 'Gam_p')
okp = r[(1, 0)][0] == 0 and r[(0, 1)][0] == 0 and r[(1, 1)][0] == 0 and r[(1, 1)][1] > 0
log(f'second-order solution Gam_p is ASD (W+^(ab) = 0, W-^(ab) != 0): {okp}   ({time.time() - t0:.0f}s)')
ok &= okp
# control: add the regular label-(+) mode at K (an SD linear solution): W+^(ab) must become nonzero
Kreg = mode(Kv, +1)
Wp2, Wm2 = weyl_sd(GamP + Kreg)
r2 = report(Wp2, Wm2, 'Gam_p + regular(+) mode at K')
okc = r2[(1, 1)][0] > 0
log(f'control (adding the SD regular mode at K makes W+^(ab) != 0): {okc}')
ok &= okc
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', f't1a_solve2_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')
pickle.dump({'GamP': sp.srepr(GamP), 'ok': ok}, open(os.path.join(HERE, 'cache', f't1a_solve2_{tag}.pkl'), 'wb'))
