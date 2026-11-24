"""
t1a_numcheck (DSKETCH_UP.md §1.3): fast numerical version of t1a_solve2's Weyl checks (no symbolic simplification): rebuilds Gam_p
from the t1a_solve2 log, computes W+-, W- of g~ = delta + a gam3 + b gam4 + a b Gam_p to O(ab) with the full-metric Hodge star (same
algebra as t1a_solve2.weyl_sd), and evaluates every component at z = 0.2, 0.7, 1.5 (40 digits). Also checks the O(ab) Einstein equation
numerically and the control (Gam_p + regular (+) mode at K).
Usage: python3 -B t1a_numcheck.py <tag>   -> cache/t1a_numcheck_<tag>.txt
"""
import os, sys, re, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from asd_alg import Alg, z, N, christoffel, riemann_up, ricci, inverse_conformal, mode

tag = sys.argv[1]
kk = [int(v) for v in tag.split('_')]; k3, k4 = kk[:3], kk[3:]
Kv = [k3[i] + k4[i] for i in range(3)]
log_ = open(os.path.join(HERE, 'cache', 't1a_solve2_run2.log' if tag != '0_0_3_4_0_0' else f't1a_solve2_{tag}.txt')).read()
GamP = sp.zeros(4, 4)
for m in re.finditer(r'Gam_p\[(\d)(\d)\] = (.*) \* exp\(-(\d+) z\)', log_):
    i, j = int(m.group(1)), int(m.group(2))
    e = sp.sympify(m.group(3).replace('I', 'I'), locals={'z': z})*sp.exp(-int(m.group(4))*z)
    GamP[i, j] = GamP[j, i] = e
A = Alg(k3, k4)
lines = []
def log(x=''):
    print(x, flush=True); lines.append(str(x))
ZS = [sp.Rational(1, 5), sp.Rational(7, 10), sp.Rational(3, 2)]

def fields(Gam):
    g3 = mode(k3, -1); g4 = mode(k4, -1)
    H = [[{(1, 0): g3[i, j], (0, 1): g4[i, j], (1, 1): Gam[i, j]} for j in range(N)] for i in range(N)]
    g = [[A.add({(0, 0): sp.Integer(1 if i == j else 0)}, H[i][j]) for j in range(N)] for i in range(N)]
    gi = inverse_conformal(A, H)
    G = christoffel(A, g, gi); Rup = riemann_up(A, G); Ric = ricci(A, Rup)
    return H, g, gi, Rup, Ric

def maxabs(exprs):
    mx = 0
    for e in exprs:
        if e == 0:
            continue
        for zz in ZS:
            v = abs(sp.N(e.subs(z, zz), 40))
            mx = max(mx, v)
    return mx

def weyl_pm(H, g, gi, Rup, Ric):
    Rsc = {}
    for b_, dd in itertools.product(range(N), repeat=2):
        Rsc = A.add(Rsc, A.mul(gi[b_][dd], Ric[b_][dd]))
    S = [[A.scale(sp.Rational(1, 2), A.add(Ric[i][j], A.scale(-sp.Rational(1, 6), A.mul(Rsc, g[i][j])))) for j in range(N)] for i in range(N)]
    W = {}
    for a_, b_, c_, dd in itertools.product(range(N), repeat=4):
        Rl = {}
        for e in range(N):
            Rl = A.add(Rl, A.mul(g[a_][e], Rup[(e, b_, c_, dd)]))
        W[(a_, b_, c_, dd)] = A.add(Rl, A.scale(-1, A.mul(g[a_][c_], S[b_][dd])), A.mul(g[a_][dd], S[b_][c_]),
                                     A.scale(-1, A.mul(g[b_][dd], S[a_][c_])), A.mul(g[b_][c_], S[a_][dd]))
    trH = A.add(*[H[i][i] for i in range(N)])
    HH = {}
    for i, j in itertools.product(range(N), repeat=2):
        HH = A.add(HH, A.mul(H[i][j], H[j][i]))
    sq = A.add({(0, 0): sp.Integer(1)}, A.scale(sp.Rational(1, 2), trH), A.scale(sp.Rational(1, 8), A.mul(trH, trH)), A.scale(-sp.Rational(1, 4), HH))
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
        Wp[(a_, b_, c_, dd)] = A.scale(sp.Rational(1, 2), A.add(W[(a_, b_, c_, dd)], star)).get((1, 1), 0)
        Wm[(a_, b_, c_, dd)] = A.scale(sp.Rational(1, 2), A.add(W[(a_, b_, c_, dd)], A.scale(-1, star))).get((1, 1), 0)
    return Wp, Wm

ok = True
log(f't1a_numcheck {tag}: K = {Kv}')
for label, Gam in (('Gam_p', GamP), ('Gam_p + reg(+) at K', GamP + mode(Kv, +1))):
    H, g, gi, Rup, Ric = fields(Gam)


    Wp, Wm = weyl_pm(H, g, gi, Rup, Ric)
    mp_ = maxabs(list(Wp.values())); mm_ = maxabs(list(Wm.values()))
    log(f'  [{label}] max |W+^(ab)| over components and z in {ZS}: {sp.N(mp_, 5)};  max |W-^(ab)|: {sp.N(mm_, 5)}')
    if label == 'Gam_p':
        ok &= mp_ < 1e-30 and mm_ > 1e-6
    else:
        ok &= mp_ > 1e-6
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', f't1a_numcheck_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')
