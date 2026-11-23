"""
p2spin/s1c_validate.py (P2_SPINNING.md §1.3): validations of the conserved-structure solve.
 (a) integer J = 0..9: exact solve with only the legal monomials (V5 exponent >= 0), sigma = (-1)^J -> counts per sector;
     containment of the integer-J space in the generic-J space specialised at that J (function values, random Delta);
 (b) J = 0 against repro_engine/tto_conserved.py's coefficients (read-only import of its formula);
 (c) independent mpmath check: Todorov divergence at x1 AND x2 of the canonical structures by numerical differentiation (30 digits)
     at random real points, non-integer J and Delta.
-> cache/s1c_validate.txt
"""
import os, sys, random, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'repro_engine'))
import sympy as sp, mpmath as mp
from fractions import Fraction as Fr
import p2_core as C
import s1_structures as S

OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))
Dl, Jl = S.Dl, S.Jl
CAN = {k: [[sp.sympify(eval(x, sp.__dict__)) for x in v] for v in B]
       for k, B in pickle.load(open(os.path.join(HERE, 'cache', 's1b_canonical.pkl'), 'rb')).items()}


def legal(sector, J):
    if sector == 'even':
        return [t for t in C.EVEN if t[1] + t[2] <= J]
    out = []
    for t in C.ODD:
        typ, a, b, c = t
        e5 = J - b - c - (0 if typ == '12' else 1)
        if e5 >= 0:
            out.append(t)
    return out


def phys_count(ns, rel):
    if not ns:
        return 0
    return sp.Matrix(ns).rank(simplify=True) - len(rel)


pr('# s1c_validate')
pr('(a) integer J: exact solve on legal monomials (V5 power >= 0); sigma = (-1)^J')
rng = random.Random(5)
for J in range(0, 10):
    sig = (-1)**J
    cnt = {}
    for sector in ('even', 'odd'):
        monos = legal(sector, J)
        if not monos:
            cnt[sector] = 0; continue
        ns, rel, _ = S.solve(sector, monos, sig, Jval=J, seed=10 + J)
        cnt[sector] = phys_count(ns, rel)
        # containment in the generic family specialised at J
        gen = CAN.get((sector, sig), [])
        if cnt[sector] and gen:
            full = C.EVEN if sector == 'even' else C.ODD
            idx = [full.index(t) for t in monos]
            Dv = sp.Rational(rng.randint(30, 400), 7)
            cfgs = [S.rcfg(rng) for _ in range(len(gen) + len(ns) + 4)]
            vals = [S.mono_vals(sector, cfg, full) for cfg in cfgs]
            rowsG = [[sum(sp.Rational(v[i].numerator, v[i].denominator)*g[i].subs({Dl: Dv, Jl: J}) for i in range(len(full)))
                      for v in vals] for g in gen]
            rowsI = [[sum(sp.Rational(v[i].numerator, v[i].denominator)*n[k].subs(Dl, Dv) for k, i in enumerate(idx))
                      for v in vals] for n in ns]
            rG = sp.Matrix(rowsG).rank(); rGI = sp.Matrix(rowsG + rowsI).rank()
            cnt[sector + '_contained'] = (rGI == rG, rG)
    pr(f'   J = {J}: parity-even {cnt["even"]}, parity-odd {cnt["odd"]}   '
       f'(integer-J space inside generic family at this J: even {cnt.get("even_contained")}, odd {cnt.get("odd_contained")})')

pr('(b) J = 0 vs repro_engine/tto_conserved.py')
ns, rel, _ = S.solve('even', [(0, 0, 0), (1, 0, 0), (2, 0, 0)], 1, Jval=0, seed=3)
v = [sp.factor(x/ns[0][2]) for x in ns[0]]
den = Dl**2 - 6*Dl + 6
ref = [Dl*(Dl + 2)/(2*den), 2*Dl*(Dl - 4)/den, 1]
pr('   own (V1^2V2^2, V1V2H12, H12^2):', v)
pr('   tto_conserved              :', [sp.factor(x) for x in ref], ' difference:', [sp.simplify(a - b) for a, b in zip(v, ref)])
import integrand as IG                     # read-only: its coef() is the tto_conserved result used by the J = 0 engine
for Dv in (300.0, 2.5):
    pr(f'   integrand.coef({Dv}) =', list(IG.coef(Dv)), ' own:', [float(x.subs(Dl, Dv)) for x in v])

pr('(c) mpmath: Todorov divergence at x1 and x2 of the canonical structures (numerical derivatives, 30 digits)')
mp.mp.dps = 30


def fullf(sector, coeffs, Dv, Jv, x1, z1, x2, z2, x5, z5):
    B = C.blocks(x1, z1, x2, z2, x5, z5)
    monos = C.EVEN if sector == 'even' else C.ODD
    f = C.mono_even if sector == 'even' else C.mono_odd
    N = sum(c*f(B, t) for c, t in zip(coeffs, monos))
    a3 = (10 - Dv - Jv)/2; a2 = a1 = (Dv + Jv)/2
    if sector == 'odd':
        a3 += mp.mpf(1)/2; a2 += mp.mpf(1)/2; a1 += mp.mpf(1)/2
    return N*B['V5']**Jv*B['X12']**(-a3)*B['X15']**(-a2)*B['X25']**(-a1)


G = [[0, -mp.mpf(1)/2, 0], [-mp.mpf(1)/2, 0, 0], [0, 0, 1]]; GI = [[0, -2, 0], [-2, 0, 0], [0, 0, 1]]


def todorov_div(F, x, z):
    """d_x^mu D_mu F(x, z) for F quadratic in z (h = 3/2), numerically."""
    h = mp.mpf(3)/2
    def Dlow(xx, m):
        g = lambda *zz: F(xx, zz)
        dz = lambda k: mp.diff(lambda t: g(*[z[i] + (t if i == k else 0) for i in range(3)]), 0)
        d2 = lambda a, b: mp.diff(lambda s, t: g(*[z[i] + (s if i == a else 0) + (t if i == b else 0) for i in range(3)]), (0, 0), (1, 1)) if a != b \
            else mp.diff(lambda t: g(*[z[i] + (t if i == a else 0) for i in range(3)]), 0, 2)
        dzm = dz(m)
        box = sum(GI[a][b]*d2(a, b) for a in range(3) for b in range(3) if GI[a][b] != 0)
        zl = sum(G[m][n]*z[n] for n in range(3))
        return h*dzm - zl*box/2          # (h - 1 + z.dz) dz_m F = h dz_m F  (dz_m F has degree 1)
    tot = 0
    for m in range(3):
        for r in range(3):
            if GI[m][r] == 0:
                continue
            tot += GI[m][r]*mp.diff(lambda t: Dlow([x[i] + (t if i == r else 0) for i in range(3)], m), 0)
    return tot


rr = random.Random(9)
rv = lambda: mp.mpf(rr.uniform(-2, 2))
def rn():
    u, y = rv(), rv(); return [u, y*y/u, y]
for (sector, sig), Bs in CAN.items():
    for Dv, Jv in ((mp.mpf('7.3'), mp.mpf('4.6')), (mp.mpf('12.1'), mp.mpf('9.3'))):
        for k, cf in enumerate(Bs):
            cnum = [mp.mpf(sp.N(c.subs({Dl: sp.Float(str(Dv), 30), Jl: sp.Float(str(Jv), 30)}), 30)) for c in cf]
            x1, z1, x2, z2, x5, z5 = [rv() for _ in range(3)], rn(), [rv() for _ in range(3)], rn(), [rv() for _ in range(3)], rn()
            F1 = lambda xx, zz: fullf(sector, cnum, Dv, Jv, xx, zz, x2, z2, x5, z5)
            F2 = lambda xx, zz: fullf(sector, cnum, Dv, Jv, x1, z1, xx, zz, x5, z5)
            d1 = todorov_div(F1, x1, z1); d2 = todorov_div(F2, x2, z2)
            # scale: same operator on the first monomial alone
            one = [mp.mpf(1) if i == 0 else 0 for i in range(len(cnum))]
            sc = abs(todorov_div(lambda xx, zz: fullf(sector, one, Dv, Jv, xx, zz, x2, z2, x5, z5), x1, z1))
            # Bose symmetry with the sigma sign for non-integer J: F(1<->2) vs sigma*F, using |V5|^J convention
            pr(f'   {sector} sigma={sig:+d} basis {k} (Delta, J) = ({float(Dv)}, {float(Jv)}): |div_1|/scale = {float(abs(d1)/sc):.1e}, '
               f'|div_2|/scale = {float(abs(d2)/sc):.1e}')
open(os.path.join(HERE, 'cache', 's1c_validate.txt'), 'w').write('\n'.join(OUT) + '\n')
