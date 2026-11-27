"""
Linear-in-coupling correction to eta/s for AdS4 black branes,  L = R + 6/L^2 + lam * O,  L_AdS = 1.  [computed]

O in { I1 = W_ab^cd W_cd^ef W_ef^ab  (Weyl^3),
       O1 = (W_abcd W^abcd)^2,  O2 = (W_abcd W~^abcd)^2,
       P3 = W W W~ (parity odd cubic) }.

Method (first principles, no effective-coupling shortcut):
  * background: ds^2 = -N(r)^2 f(r) dt^2 + dr^2/f(r) + r^2 (dx^2+dy^2); reduced action in (N,f), Euler-Lagrange
    with all higher derivatives, solved to O(lam) with horizon kept at r = 1 and boundary metric conformal to
    -dt^2 + dx^2 (N^2 f / r^2 -> 1).
  * entropy: Wald, s = -8 pi E^{trtr} g_tt g_rr r^2  (E^{abcd} = dL/dR_abcd, Riemann treated as independent).
  * shear viscosity: h^x_y = phi(r) (static) on that background; the radial quadratic Lagrangian is reduced by
    integrations by parts to alpha phi''^2 + beta phi'^2 + mu phi^2 (mu = 0 checked), whose conserved momentum is
    Pi = 2 beta chi' - (2 alpha chi'')'.  For the ingoing solution phi = e^{-i w t}e^{-i w r_*}(1 + i w g), g regular,
    the radial flux at O(w) is w * Pi[-r_*] evaluated at the horizon (alpha ~ (r-1)^2 checked so that g drops out).
    Normalisation fixed by Einstein gravity (eta/s = 1/4pi); the flux->G_R constant is coupling independent.
"""
import sympy as sp

r, lam, eps = sp.symbols('r lambda epsilon', positive=True)
t, x, y = sp.symbols('t x y')
X = [t, r, x, y]


def curvature(g):
    ginv = sp.simplify(g.inv())
    n = 4
    Gam = [[[sp.simplify(sum(ginv[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                             for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]
    # R^a_{bcd}
    Rud = sp.MutableDenseNDimArray.zeros(n, n, n, n)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(c + 1, n):
                    v = (sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
                         + sum(Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c] for e in range(n)))
                    v = sp.simplify(v)
                    Rud[a, b, c, d] = v
                    Rud[a, b, d, c] = -v
    Rl = sp.MutableDenseNDimArray.zeros(n, n, n, n)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    Rl[a, b, c, d] = sp.simplify(sum(g[a, e] * Rud[e, b, c, d] for e in range(n)))
    return ginv, Rl


def weyl_from_riem(Rl, g, ginv):
    """Algebraic map Riemann(lower) -> Weyl(lower), D = 4."""
    n = 4
    Ric = sp.Matrix(n, n, lambda b, d: sum(ginv[a, c] * Rl[a, b, c, d] for a in range(n) for c in range(n)))
    Rs = sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    W = sp.MutableDenseNDimArray.zeros(n, n, n, n)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    W[a, b, c, d] = (Rl[a, b, c, d]
                                     - sp.Rational(1, 2) * (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c]
                                                            - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c])
                                     + Rs / 6 * (g[a, c] * g[b, d] - g[a, d] * g[b, c]))
    return W


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def bivector_ops(W, g, ginv):
    """M[A,B] = W_{ab}^{cd} (A=(a<b), B=(c<d)); star[A,B] = (1/2)eps_{ab}^{cd}*2 restricted -> eps_{ab}^{cd}."""
    n = 4
    # W_{ab}^{cd} = W_{ab ef} g^{ec} g^{fd}
    M = sp.zeros(6, 6)
    for A, (a, b) in enumerate(PAIRS):
        for B, (c, d) in enumerate(PAIRS):
            M[A, B] = sum(W[a, b, e, f] * ginv[e, c] * ginv[f, d] for e in range(n) for f in range(n))
    # eps_{abcd} = sqrt(-g) [abcd]; eps_{ab}^{cd}
    sqrtmg = sp.sqrt(-g.det())
    S = sp.zeros(6, 6)
    for A, (a, b) in enumerate(PAIRS):
        for B, (c, d) in enumerate(PAIRS):
            S[A, B] = sum(sqrtmg * sp.LeviCivita(a, b, e, f) * ginv[e, c] * ginv[f, d] for e in range(n) for f in range(n))
    # On bivectors with sums over a<b: (1/2) X_{ab}^{ef} Y_{ef}^{cd} (ordered) = (X Y)[A,B] with X = M
    # W~_{ab}^{cd} = 1/2 eps_{ab}^{ef} W_{ef}^{cd} = (S M)[A,B]
    return M, S


def scalar(M, S, which):
    # traces over ordered index pairs: sum_{ab,cd,...} = 2^k tr over a<b
    if which == 'I1':
        return 8 * (M ** 3).trace()
    if which == 'W2':
        return 4 * (M ** 2).trace()
    if which == 'O1':
        return (4 * (M ** 2).trace()) ** 2
    if which == 'O2':
        return (4 * (S * M * M).trace()) ** 2
    if which == 'P3':
        return 8 * (S * M ** 3).trace()
    raise ValueError


def lagr_density(g, which, Rl=None):
    ginv, Rl0 = curvature(g)
    Rl = Rl0 if Rl is None else Rl
    W = weyl_from_riem(Rl, g, ginv)
    M, S = bivector_ops(W, g, ginv)
    return ginv, Rl0, scalar(M, S, which)


N, f, phi = sp.Function('N')(r), sp.Function('f')(r), sp.Function('phi')(r)


def EL(Lag, F):
    derivs = [F] + [sp.diff(F, r, k) for k in range(1, 5)]
    out = 0
    for k, Dk in enumerate(derivs):
        out += (-1) ** k * sp.diff(sp.diff(Lag, Dk), r, k)
    return out


def run(which):
    print(f"\n===== O = {which} =====")
    gB = sp.diag(-N ** 2 * f, 1 / f, r ** 2, r ** 2)
    ginv, Rl, Ob = lagr_density(gB, which)
    sqrtmg = N * r ** 2
    Ob = sp.simplify(Ob)
    LE = sqrtmg * sp.simplify(sum(ginv[a, c] * ginv[b, d] * Rl[a, b, c, d] for a in range(4) for b in range(4)
                                  for c in range(4) for d in range(4)) + 6)
    LO = sqrtmg * Ob
    # background to O(lam)
    f0, N0 = r ** 2 - 1 / r, sp.Integer(1)
    f1, n1 = sp.Function('f1')(r), sp.Function('n1')(r)
    subsE = {N: N0 + lam * n1, f: f0 + lam * f1}
    eqs = []
    for F in (N, f):
        e = EL(LE + lam * LO, F)
        e = e.subs(subsE).doit()
        e = sp.series(e, lam, 0, 2).removeO()
        e0 = sp.simplify(e.coeff(lam, 0))
        assert e0 == 0, e0
        eqs.append(sp.simplify(e.coeff(lam, 1)))
    sol = sp.dsolve(eqs, [f1, n1])
    sol = {s.lhs: s.rhs for s in sol}
    f1s, n1s = sp.simplify(sol[f1]), sp.simplify(sol[n1])
    consts = sorted((f1s.free_symbols | n1s.free_symbols) - {r, lam}, key=str)
    # impose f1(1) = 0 (horizon at r=1), boundary metric -dt^2+dx^2 conformally: N^2 f/r^2 -> 1
    cond = [sp.Eq(f1s.subs(r, 1), 0)]
    asym = sp.limit((2 * n1s * f0 + f1s) / r ** 2, r, sp.oo)
    cond.append(sp.Eq(asym, 0))
    # third constant (if any): fix by requiring no r^{-2}-growth beyond what EOM forces: set extra homogeneous pieces to zero
    cs = sp.solve(cond, consts[:2], dict=True)
    f1s, n1s = f1s.subs(cs[0]), n1s.subs(cs[0])
    rem = (f1s.free_symbols | n1s.free_symbols) - {r, lam}
    if rem:
        print("  remaining constants (set to 0):", rem)
        f1s, n1s = f1s.subs({c: 0 for c in rem}), n1s.subs({c: 0 for c in rem})
    print("  f1 =", sp.simplify(f1s), "   n1 =", sp.simplify(n1s))
    Nb, fb = 1 + lam * n1s, f0 + lam * f1s
    T = sp.series(sp.limit(Nb * sp.diff(fb, r), r, 1) / (4 * sp.pi), lam, 0, 2).removeO()
    print("  4 pi T =", sp.simplify(4 * sp.pi * T))

    # ---- Wald entropy: E^{trtr} via Riemann perturbation X_{trtr}
    d = sp.symbols('d')
    Xt = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    Xt[0, 1, 0, 1] = Xt[1, 0, 1, 0] = 1
    Xt[0, 1, 1, 0] = Xt[1, 0, 0, 1] = -1
    ginvB, RlB = ginv, Rl
    Rpert = sp.MutableDenseNDimArray([[[[RlB[a, b, c, e] + d * Xt[a, b, c, e] for e in range(4)] for c in range(4)]
                                       for b in range(4)] for a in range(4)])
    # Lagrangian (without sqrt-g) as function of Riemann: R + 6 + lam O
    Rsc = sum(ginvB[a, c] * ginvB[b, e] * Rpert[a, b, c, e] for a in range(4) for b in range(4)
              for c in range(4) for e in range(4))
    W = weyl_from_riem(Rpert, gB, ginvB)
    M, S = bivector_ops(W, gB, ginvB)
    Lfull = Rsc + 6 + lam * scalar(M, S, which)
    Ettrr4 = sp.diff(Lfull, d).subs(d, 0)  # = 4 E^{trtr}  (per 16 pi G)
    sW = -8 * sp.pi * (Ettrr4 / 4) * (-gB[0, 0] * gB[1, 1]) * r ** 2   # units 1/(16 pi G); expect 4 pi at lam=0
    sW = sW.subs({N: Nb, f: fb}).doit()
    sW = sp.series(sp.simplify(sW.subs(r, 1)), lam, 0, 2).removeO()
    print("  s (units 1/16piG) =", sp.simplify(sW))

    # ---- quadratic action for h^x_y = eps*phi(r)
    gP = sp.Matrix([[-N ** 2 * f, 0, 0, 0], [0, 1 / f, 0, 0], [0, 0, r ** 2, r ** 2 * eps * phi],
                    [0, 0, r ** 2 * eps * phi, r ** 2]])
    ginvP, RlP, OP = lagr_density(gP, which)
    sqrtP = sp.sqrt(-gP.det())
    LEP = sqrtP * (sum(ginvP[a, c] * ginvP[b, e] * RlP[a, b, c, e] for a in range(4) for b in range(4)
                       for c in range(4) for e in range(4)) + 6)
    LOP = sqrtP * OP
    L2E = sp.diff(LEP, eps, 2).subs(eps, 0) / 2
    L2O = sp.diff(LOP, eps, 2).subs(eps, 0) / 2
    p0, p1, p2 = sp.symbols('p0 p1 p2')
    rep = {sp.diff(phi, r, 2): p2, sp.diff(phi, r): p1, phi: p0}

    def coeffs(L2):
        L2 = sp.expand(L2.subs({sp.diff(phi, r, 2): p2}).subs({sp.diff(phi, r): p1}).subs({phi: p0}))
        P = sp.Poly(L2, p0, p1, p2)
        get = lambda m: P.coeff_monomial(m)
        return dict(a1=get(p2 ** 2), a2=get(p1 ** 2), a3=get(p0 ** 2), a4=get(p1 * p2), a5=get(p0 * p2), a6=get(p0 * p1))

    cE, cO = coeffs(L2E), coeffs(L2O)
    for k in ('a1', 'a2', 'a3', 'a4', 'a5', 'a6'):
        assert not cE[k].has(eps)
    # assemble on background to O(lam): Einstein part on corrected bg, O part on Einstein bg
    def onbg(expr, Nv, fv):
        return expr.subs({N: Nv, f: fv}).doit()

    C = {}
    for k in cE:
        C[k] = sp.simplify(sp.series(onbg(cE[k], Nb, fb) + lam * onbg(cO[k], 1, f0), lam, 0, 2).removeO())
    alpha = C['a1']
    beta = sp.simplify(C['a2'] - sp.diff(C['a4'], r) / 2 - C['a5'])
    mu = sp.simplify(C['a3'] + sp.diff(C['a5'], r, 2) / 2 - sp.diff(C['a6'], r) / 2)
    mu = sp.simplify(sp.series(mu, lam, 0, 2).removeO())
    print("  effective mass term mu (expect 0):", mu)
    print("  alpha =", sp.factor(alpha))
    rstar_p = 1 / (Nb * fb)
    Pi_h = 2 * beta * rstar_p - sp.diff(2 * alpha * sp.diff(rstar_p, r), r)
    Pi_h = sp.series(Pi_h, lam, 0, 2).removeO()
    Pi_h = sp.limit(sp.simplify(Pi_h), r, 1)
    print("  -Pi[-r*] at horizon =", sp.simplify(Pi_h))
    # Einstein calibration: eta_E = K * Pi_h(lam=0) must equal s_E/(4 pi)
    K = (sW.subs(lam, 0) / (4 * sp.pi)) / Pi_h.subs(lam, 0)
    etas = sp.series(K * Pi_h / sW, lam, 0, 2).removeO()
    res = sp.simplify(sp.expand(4 * sp.pi * etas))
    print("  4 pi eta/s =", res)
    return res


if __name__ == '__main__':
    import sys
    whichs = sys.argv[1:] or ['I1', 'O1', 'O2', 'P3']
    for w in whichs:
        run(w)
