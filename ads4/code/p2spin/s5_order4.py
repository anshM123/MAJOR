"""
p2spin/s5_order4.py (P2_SPINNING.md §7, follow-up (a)): O(1/m^4) of the seed and of the spinning symbols; K4 in Theorem W units.

* Seed (scalar CHMRS T, Delta_phi = 3, k = 2): exact expansion to eps^4; mechanism (s3_mech, memoised derivatives, Legendre reduction
  to order 10); validated against CHMRS's ancillary 'series Pi' at order 1/m2^2 (parsed from regge_moments_order8.m, read-only,
  scratchpad copy given as argv[1]) -- all nu-dependent coefficients.
* Spinning symbols (even J, helicity classes same/opp): exact series to eps^4 at q in QS, polynomials of degree <= k in (X, Y)
  by exact interpolation (210 monomials + 6 check points); helicity-diagonal model [M] of s3_spin, direct pairing.
* K4: the nu^4 coefficients of P and P' at O(1/m^4), in Theorem W's units eta^2 nu^4 P and eta(eta^2-1) nu^4 P'.
-> cache/s5_order4.txt
"""
import os, sys, random, time, re
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import s3_mech as M

OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))
A = 7
QS = [sp.Rational(3, 2), sp.Rational(2)]


# ------------------------------------------------------------------ seed to eps^4
def seed4(qv):
    e = (1, 0, 0); n = (1, 1, 0); EPS = M.EPS; X, Y = M.X, M.Y; mdot = M.mdot
    x = [EPS*c for c in X]; y = [EPS*c for c in Y]
    D = ((qv + 1/qv)/EPS + 3)/2; J = ((qv - 1/qv)/EPS - 1)/2
    def V(sg):
        ee = [sg*c for c in e]
        a = (1 - 2*mdot(ee, x))*(1 - 2*mdot(ee, y)) - mdot(x, x)*mdot(y, y)
        ex2 = mdot([ee[i] + x[i] for i in range(3)], [ee[i] + x[i] for i in range(3)])
        ey2 = mdot([ee[i] + y[i] for i in range(3)], [ee[i] + y[i] for i in range(3)])
        ybs = [-y[i] - 2*ee[i]*mdot(y, ee) for i in range(3)] if sg == 1 else [-y[i] - 2*e[i]*mdot(y, e) for i in range(3)]
        return [a*ee[i] + ey2*x[i] + ex2*ybs[i] for i in range(3)]
    sq = lambda u: mdot(u, u); em = lambda u, s: [e[i] + s*u[i] for i in range(3)]
    terms = [(-4, 1 + 2*mdot(x, y) + sq(x)*sq(y)), (-1, sq(em(x, -1))), (-1, sq(em(y, -1))),
             ((D + J)/2, sq(em(x, -1))*sq(em(y, -1))/(sq(em(x, 1))*sq(em(y, 1)))),
             (J, -mdot(n, V(1))), (-(J + 1), mdot(n, V(-1)))]
    L = 0
    for pw, arg in terms:
        a0 = sp.expand(arg.subs(EPS, 0))
        L += pw*sp.series(sp.log(sp.expand(arg/a0)), EPS, 0, 6).removeO()
    L = sp.expand(sp.series(sp.expand(L), EPS, 0, 5).removeO())
    Lk = [sp.expand(L.coeff(EPS, k)) for k in range(5)]
    assert Lk[1] == 0 and Lk[3] == 0
    return sp.expand(Lk[2]), sp.expand(Lk[4] + Lk[2]**2/2)


# ------------------------------------------------------------------ memoised mechanism
class Mech:
    def __init__(self, P0, Pb0):
        self.Ps, self.Qs, G, F = M.deriv_table(A)
        self.G = G; self.F = F; self.memo = {(0,)*6: F}
        self.sub = dict(zip(self.Ps, P0)); self.sub.update(dict(zip(self.Qs, Pb0)))
        nP = sp.sqrt(-M.mdot(self.Ps, self.Ps)); nQ = sp.sqrt(-M.mdot(self.Qs, self.Qs))
        self.F0 = (nP**(-A)*nQ**(-A)).subs(self.sub)
        self.et0 = sp.nsimplify(-M.mdot(P0, Pb0)); self.vals = {}
        self.gs = [sp.Symbol(f'G{k}') for k in range(12)]

    def d(self, mon):
        if mon in self.memo:
            return self.memo[mon]
        i = max(j for j in range(6) if mon[j] > 0)
        prev = list(mon); prev[i] -= 1; prev = tuple(prev)
        var = self.Ps[i] if i < 3 else self.Qs[i - 3]
        r = sp.diff(self.d(prev), var)*sp.Rational(1, 2)*M.MET[i % 3, i % 3]
        self.memo[mon] = r
        return r

    def val(self, mon):
        if mon not in self.vals:
            ex = self.d(mon).subs(self.sub)
            w = sp.Symbol('w')
            ex = ex.doit()
            rep = {}
            for at in ex.atoms(sp.Subs):
                dv = at.expr
                if isinstance(dv, sp.Derivative) and dv.expr.func == self.G:
                    rep[at] = self.gs[sum(c for _, c in dv.variable_count)]
            for at in ex.atoms(sp.Derivative):
                if at.expr.func == self.G:
                    rep[at] = self.gs[sum(c for _, c in at.variable_count)]
            ex = ex.xreplace(rep)
            ex = ex.subs(self.G(self.et0), self.gs[0])
            assert not ex.has(self.G), ex
            self.vals[mon] = sp.expand(sp.nsimplify(sp.expand(ex/self.F0)))
        return self.vals[mon]

    def sol(self):
        e = self.et0
        sol = {0: (1, 0), 1: (0, 1)}
        for k in range(0, 10):
            c0 = -(2*(k + 1)*e*sol[k + 1][0] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][0])/(e**2 - 1)
            c1 = -(2*(k + 1)*e*sol[k + 1][1] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][1])/(e**2 - 1)
            sol[k + 2] = (sp.expand(c0), sp.expand(c1))
        return sol

    def apply(self, poly):
        P = sp.Poly(sp.expand(poly), *M.X, *M.Y)
        tot = sum(c*self.val(mon) for mon, c in P.terms())
        e = self.et0
        sol = {0: (1, 0), 1: (0, 1)}
        for k in range(0, 10):
            c0 = -(2*(k + 1)*e*sol[k + 1][0] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][0])/(e**2 - 1)
            c1 = -(2*(k + 1)*e*sol[k + 1][1] + (k*(k + 1) + M.nu2 + sp.Rational(1, 4))*sol[k][1])/(e**2 - 1)
            sol[k + 2] = (sp.expand(c0), sp.expand(c1))
        tot = sp.expand(tot)
        o0 = sum(tot.coeff(self.gs[k])*sol[k][0] for k in range(12)); o1 = sum(tot.coeff(self.gs[k])*sol[k][1] for k in range(12))
        return sp.Poly(sp.expand(o0), M.nu2), sp.Poly(sp.expand(o1), M.nu2)


# ------------------------------------------------------------------ CHMRS ancillary, order 1/m2^2
def chmrs_order2(path, d=3, Dphi=3, alpha=7):
    txt = open(path).read()
    s = txt[txt.index('series\\[CapitalPi] ='):txt.index('seriesC =')]
    s = s[s.index('=') + 1:].strip().rstrip(';')
    for a, b in (('\\[Eta]AdS', 'eta'), ('\\[Nu]', 'nu'), ('\\[CapitalDelta]\\[Phi]', 'Dphi'), ('\\[Alpha]', 'alpha'), ('^', '**')):
        s = s.replace(a, b)
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'pj\[\s*\(2 - d\)/2 \+ I\*nu,\s*eta,\s*(\d)\s*\]', r'G\1', s)
    assert 'pj' not in s, s[s.index('pj'):s.index('pj') + 80]
    nu, eta, m2 = sp.symbols('nu eta m2')
    ex = sp.sympify(s, locals={'nu': nu, 'eta': eta, 'm2': m2, 'd': sp.Integer(d), 'Dphi': sp.Integer(Dphi),
                                'alpha': sp.Integer(alpha), 'I': sp.I, **{f'G{k}': sp.Symbol(f'G{k}') for k in range(9)}})
    ex = sp.expand(ex)
    t2 = sp.expand(ex.coeff(m2, -2)/2)
    return t2, nu, eta


def main():
    t0 = time.time()
    anc = sys.argv[1] if len(sys.argv) > 1 else None
    seeds = {}
    for qv in QS:
        s2, s4 = seed4(qv)
        P, Pb, _ = M.seed_series(qv)
        mech = Mech(P, Pb)
        c0, c1 = mech.apply(s4)
        et = mech.et0
        seeds[qv] = (P, Pb, mech, c0, c1)
        pr(f'q = {qv} (eta = {et}): seed O(1/m^4): P: {c0.as_expr()}   P\': {c1.as_expr()}   [{time.time() - t0:.0f}s]')
        if anc:
            t2, nu, eta = chmrs_order2(anc)
            t2 = sp.expand(t2.subs(eta, et))
            sol = mech.sol()
            t2 = sp.expand(t2.subs(nu, sp.sqrt(M.nu2)))
            r0 = sp.Poly(sp.expand(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][0] for k in range(9))), M.nu2)
            r1 = sp.Poly(sp.expand(sum(t2.coeff(sp.Symbol(f'G{k}'))*sol[k][1] for k in range(9))), M.nu2)
            dP = [sp.nsimplify(c0.as_expr().coeff(M.nu2, k) - r0.as_expr().coeff(M.nu2, k)) for k in range(4)]
            dD = [sp.nsimplify(c1.as_expr().coeff(M.nu2, k) - r1.as_expr().coeff(M.nu2, k)) for k in range(4)]
            pr(f'   vs CHMRS series Pi (1/m2^2): mine - CHMRS, P: nu^0..nu^4 = {dP};  P\': {dD}')
            m2P, m2D = M.apply_insertion(s2, P, Pb, A)[0:1][0], None
            e1, et1, reps1 = M.apply_insertion(s2, P, Pb, A); a0, a1 = M.reduce_legendre(e1, reps1, et1)
            a0 = sp.expand(a0); a1 = sp.expand(a1)
            pr(f'   expected if CHMRS = (1 - 8 eta/m^2 + n2/m^4) x mine:  P\' diff = 8 eta x (first-order P\') = '
               f'{[sp.nsimplify(8*et*a1.coeff(M.nu2, k)) for k in range(3)]};  P diff (nu^2) = 8 eta x (first-order nu^2 P) = {sp.nsimplify(8*et*a0.coeff(M.nu2, 1))}')
        n4P = c0.as_expr().coeff(M.nu2, 2); n4D = c1.as_expr().coeff(M.nu2, 2)
        pr(f'   seed nu^4: P {float(n4P):.5f} = {float(n4P/et**2):.5f} eta^2 ; P\' {float(n4D):.5f} = {float(n4D/(et*(et**2 - 1))):.5f} eta(eta^2-1)')
    # ---- spinning symbols to eps^4 (even J), interpolation
    import p2_series as PS
    PS.ORDER = 4
    import s2_symbols as S2, s2_main as SM
    PS.ORDER = 4                      # s2_main sets 2 at import
    can = S2.CAN()
    allv = list(S2.xs) + list(S2.ys)
    mons = [m for deg in range(5) for m in sp.itermonomials(allv, deg) if sp.Poly(m, *allv).total_degree() == deg]
    mons = list(dict.fromkeys(mons))
    rng = random.Random(7)
    def rpt():
        while True:
            u1, v1, u2, v2 = [sp.Rational(rng.randint(1, 40), rng.randint(1, 11)) for _ in range(4)]
            y1 = sp.Rational(rng.randint(-20, 20), rng.randint(5, 13)); y2 = sp.Rational(rng.randint(-20, 20), rng.randint(5, 13))
            if y1**2 < u1*v1 and y2**2 < u2*v2:
                return (u1, v1, y1, u2, v2, y2)
    pts = [rpt() for _ in range(len(mons) + 6)]
    Vm = DomainMatrix.from_Matrix(sp.Matrix([[m.subs(dict(zip(allv, p))) for m in mons] for p in pts[:len(mons)]])).convert_to(sp.QQ)
    MODES = [(2, 2), (2, -2)]
    for qv in QS:
        P, Pb, mech, s0c, s1c = seeds[qv]
        vals = {}
        for k in (0, 1):
            cf = can[('even', 1)][k]
            per = []
            for p in pts:
                v, ms = SM.sym_modes('L', ('even', 1), cf, p[:3], p[3:], qv)
                per.append(ms)
            vals[k] = per
            assert len(per[0]) == 5
        pr(f'q = {qv}: symbols evaluated at {len(pts)} points   [{time.time() - t0:.0f}s]')
        HB = {h: [sp.sympify(c).subs(sp.Symbol('q', positive=True), qv) for c in SM_HB(h)] for h in ('same', 'opp')}
        L = {}
        keys = [(h, md, o, part) for h in ('same', 'opp') for md in MODES for o in range(5) for part in (0, 1)]
        cols = {}; yvs = {}
        for (h, md, o, part) in keys:
            if (h, md, o) not in yvs:
                yvs[(h, md, o)] = [sp.nsimplify(sp.expand(HB[h][0]*vals[0][i][o].get(md, 0) + HB[h][1]*vals[1][i][o].get(md, 0)))
                                   for i in range(len(pts))]
            cols[(h, md, o, part)] = [(sp.re(v) if part == 0 else sp.im(v)) for v in yvs[(h, md, o)]]
        Bm = DomainMatrix.from_Matrix(sp.Matrix([[cols[k][i] for k in keys] for i in range(len(mons))])).convert_to(sp.QQ)
        Cm = Vm.lu_solve(Bm).to_Matrix()                       # one factorisation, all right-hand sides
        for (h, md) in [(h, md) for h in ('same', 'opp') for md in MODES]:
            ords = []
            for o in range(5):
                poly = 0
                for part, unit in ((0, 1), (1, sp.I)):
                    j = keys.index((h, md, o, part))
                    poly += unit*sum(Cm[i, j]*m for i, m in enumerate(mons))
                for i in range(len(mons), len(pts)):
                    r = sp.expand(poly.subs(dict(zip(allv, pts[i]))) - yvs[(h, md, o)][i])
                    assert r == 0, ('interpolation check failed', h, md, o, r)
                ords.append(sp.expand(poly))
            L[(h, md)] = ords
        pr(f'   interpolation exact (6 check points, residual 0)   [{time.time() - t0:.0f}s]')
        LC = {S2.xs[0]: M.X[0] + M.X[1], S2.xs[1]: M.X[0] - M.X[1], S2.xs[2]: M.X[2],
              S2.ys[0]: M.Y[0] + M.Y[1], S2.ys[1]: M.Y[0] - M.Y[1], S2.ys[2]: M.Y[2]}
        for h in ('same', 'opp'):
            tot4 = 0
            for md in MODES:                                   # both classes (drift of the other class enters at O(eps^2))
                lead = L[(md == (2, 2) and 'same' or 'opp', md)][0]
                l = [sp.expand(c/lead) for c in L[(h, md)]]
                tot4 += sum(l[i]*l[4 - i]*(-1)**(4 - i) for i in range(5))
            tot4 = sp.expand(tot4)
            tot4 = sp.expand(((tot4 + tot4.subs(sp.I, -sp.I))/2).subs(LC))
            # full spinning O(eps^4) = seed s4 + symbol_4 + seed s2 * symbol_2 ; symbol_2 from the same series
            tot2 = 0
            for md in MODES:
                lead = L[(md == (2, 2) and 'same' or 'opp', md)][0]
                l = [sp.expand(c/lead) for c in L[(h, md)]]
                tot2 += sum(l[i]*l[2 - i]*(-1)**(2 - i) for i in range(3))
            tot2 = sp.expand(((tot2 + tot2.subs(sp.I, -sp.I))/2).subs(LC))
            s2, _ = seed4(qv)
            zero = {v: 0 for v in M.X + M.Y}
            c2 = tot2.subs(zero)                                    # normalisation drift (congruence, W1c) factored out
            ins = sp.expand(tot4 + s2*tot2 - c2*(s2 + tot2)); ins = ins - ins.subs(zero)
            c0, c1 = mech.apply(ins)
            et = mech.et0
            n4P = c0.as_expr().coeff(M.nu2, 2); n4D = c1.as_expr().coeff(M.nu2, 2)
            n2P = c0.as_expr().coeff(M.nu2, 1); n2D = c1.as_expr().coeff(M.nu2, 1)
            pr(f'   {h}: spinning symbol (incl. seed x symbol cross term) at O(1/m^4): nu^4 P {float(n4P):+.5f}, nu^4 P\' {float(n4D):+.5f}; '
               f'nu^2 P {float(n2P):+.4f}, nu^2 P\' {float(n2D):+.4f}; '
               f'K4 contribution (model): {abs(float(n4P/et**2)):.4f} (eta^2 nu^4 P), {abs(float(n4D/(et*(et**2 - 1)))):.4f} (eta sinh b nu^4 dP)')
    open(os.path.join(HERE, 'cache', 's5_order4.txt'), 'w').write('\n'.join(OUT) + '\n')


def SM_HB(h):
    q = sp.Symbol('q', positive=True)
    return {'same': (-(q**2 + 1)**2/q**2, -2*(2*q**6 + 4*q**4 - 6*q**2 - 1)/q**4),
            'opp': ((q**2 + 1)**2/q**2, 2*(2*q**6 + 4*q**4 - 6*q**2 + 1)/q**4)}[h]


if __name__ == '__main__':
    main()
