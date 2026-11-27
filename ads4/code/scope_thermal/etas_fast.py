"""
Fast version of etas_ads4.py.  [computed]

At linear order in the coupling only the correction operator's contributions on the *Einstein* black brane are
needed: with the horizon held at r = 1, the Einstein-term contributions to both eta and s equal r_h^2/(16 pi G) x
(1, 4 pi) on any static background of the form -N^2 f dt^2 + dr^2/f + r^2 dx^2 (membrane value sqrt(-g/(g_tt g_rr))
g_xx^0 ... = r_h^2, Wald area term), so the O(lam) background correction drops out of eta/s.

  4 pi eta/s - 1 = eta_1/eta_0 - s_1/s_0,
  eta  ~ -Pi[-r_*]_h = 2 beta r_*' - (2 alpha r_*'')'  at r -> 1,   r_*' = 1/f0
  s    = Wald.

Operators (all built from the Weyl tensor W; W~ = 1/2 eps W):
  I2 = W^2, I3 = W_ab^cd W_cd^ef W_ef^ab, O2 = (W W~)^2, P3 = W W W~, I2^2, I2 I3, I2^3, I3^2.
Also prints I2, I3, dI2/dr, dI3/dr at the horizon.
"""
import sys
import sympy as sp

r, eps = sp.symbols('r epsilon', positive=True)
t, x, y = sp.symbols('t x y')
X = [t, r, x, y]
phi = sp.Function('phi')(r)
f0 = r ** 2 - 1 / r
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def riemann(g, simp=sp.cancel):
    ginv = g.inv()
    ginv = ginv.applyfunc(simp)
    n = 4
    dg = [[[sp.diff(g[a, b], X[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[simp(sum(ginv[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    Rl = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    if c == d:
                        Rl[a, b, c, d] = 0
                        continue
                    v = 0
                    for e in range(n):
                        if g[a, e] == 0:
                            continue
                        v += g[a, e] * (sp.diff(Gam[e][b][d], X[c]) - sp.diff(Gam[e][b][c], X[d])
                                        + sum(Gam[e][c][h] * Gam[h][b][d] - Gam[e][d][h] * Gam[h][b][c] for h in range(n)))
                    Rl[a, b, c, d] = simp(v)
    return ginv, Rl


def weyl(Rl, g, ginv):
    n = 4
    Ric = sp.Matrix(n, n, lambda b, d: sum(ginv[a, c] * Rl[a, b, c, d] for a in range(n) for c in range(n)))
    Rs = sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n))
    W = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    W[a, b, c, d] = (Rl[a, b, c, d] - sp.Rational(1, 2) * (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c]
                                                                          - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c])
                                     + Rs / 6 * (g[a, c] * g[b, d] - g[a, d] * g[b, c]))
    return W, Rs


def bivec(W, g, ginv):
    M = sp.zeros(6, 6)
    S = sp.zeros(6, 6)
    sq = sp.sqrt(-g.det())
    for A, (a, b) in enumerate(PAIRS):
        for B, (c, d) in enumerate(PAIRS):
            M[A, B] = sum(W[a, b, e, h] * ginv[e, c] * ginv[h, d] for e in range(4) for h in range(4))
            S[A, B] = sum(sq * sp.LeviCivita(a, b, e, h) * ginv[e, c] * ginv[h, d] for e in range(4) for h in range(4))
    return M, S


def invariants(M, S):
    I2 = 4 * (M ** 2).trace()
    I3 = 8 * (M ** 3).trace()
    J2 = 4 * (S * M * M).trace()          # W W~
    P3 = 8 * (S * M ** 3).trace()         # W W W~
    return I2, I3, J2, P3


OPS = {
    'I2': lambda I2, I3, J2, P3: I2,
    'I3': lambda I2, I3, J2, P3: I3,
    'I2sq': lambda I2, I3, J2, P3: I2 ** 2,
    'O2': lambda I2, I3, J2, P3: J2 ** 2,
    'P3': lambda I2, I3, J2, P3: P3,
    'I2I3': lambda I2, I3, J2, P3: I2 * I3,
    'I2cube': lambda I2, I3, J2, P3: I2 ** 3,
    'I3sq': lambda I2, I3, J2, P3: I3 ** 2,
    'I2J2': lambda I2, I3, J2, P3: I2 * J2,   # parity odd quartic-ish (W^2)(W W~)
}

# ---------------- background quantities
gB = sp.diag(-f0, 1 / f0, r ** 2, r ** 2)
ginvB, RlB = riemann(gB)
WB, RsB = weyl(RlB, gB, ginvB)
MB, SB = bivec(WB, gB, ginvB)
invB = [sp.simplify(v) for v in invariants(MB, SB)]
print("background: I2, I3, WW~, WWW~ =", invB)
print("  at horizon: I2 =", invB[0].subs(r, 1), " I3 =", invB[1].subs(r, 1),
      " dI2/dr =", sp.diff(invB[0], r).subs(r, 1), " dI3/dr =", sp.diff(invB[1], r).subs(r, 1))

# Wald: perturb Riemann by d * X (trtr components), Riemann treated as independent
d = sp.symbols('d')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: RlB[k] + d * Xt.get(k, 0) for k in RlB}
Wp, Rsp = weyl(Rp, gB, ginvB)
Mp, Sp_ = bivec(Wp, gB, ginvB)
invP = invariants(Mp, Sp_)
wald_fac = -8 * sp.pi * (-gB[0, 0] * gB[1, 1]) * r ** 2 / 4   # s = wald_fac * dL/dd   (units 1/16piG)
s0 = sp.simplify((wald_fac * sp.diff(Rsp, d)).subs(d, 0).subs(r, 1))
print("Einstein Wald s0 (units 1/16piG, expect 4 pi):", s0)

# ---------------- perturbed metric (static shear h^x_y = eps phi)
gP = sp.Matrix([[-f0, 0, 0, 0], [0, 1 / f0, 0, 0], [0, 0, r ** 2, r ** 2 * eps * phi], [0, 0, r ** 2 * eps * phi, r ** 2]])
ginvP, RlP = riemann(gP)
WP, RsP = weyl(RlP, gP, ginvP)
MP, SP = bivec(WP, gP, ginvP)
sqrtP = sp.sqrt(-gP.det())
invPert = invariants(MP, SP)
p0, p1, p2 = sp.symbols('p0 p1 p2')


def second_order(expr):
    """eps^2 coefficient, as polynomial in (phi, phi', phi'')."""
    e = expr.subs({sp.diff(phi, r, 2): p2}).subs({sp.diff(phi, r): p1}).subs({phi: p0})
    e2 = sp.diff(e, eps, 2).subs(eps, 0) / 2
    return sp.expand(sp.simplify(e2))


# expansions of the invariants to O(eps^2) as series: I = I^(0) + eps I^(1) + eps^2 I^(2)
def eps_series(expr):
    e = expr.subs({sp.diff(phi, r, 2): p2}).subs({sp.diff(phi, r): p1}).subs({phi: p0})
    return [sp.simplify(sp.diff(e, eps, k).subs(eps, 0) / sp.factorial(k)) for k in range(3)]


ser = [eps_series(v) for v in invPert]
sq_ser = eps_series(sqrtP)
print("first-order (eps^1) pieces of I2, I3, WW~, WWW~ (expect 0 for the shear mode):", [s_[1] for s_ in ser])
print("second-order pieces: dI2 =", sp.factor(ser[0][2]), "\n   dI3 =", sp.factor(ser[1][2]),
      "\n   d(WW~) =", sp.factor(ser[2][2]), "\n   d(WWW~) =", sp.factor(ser[3][2]))
RsSer = eps_series(RsP)


def L2_of(fun):
    """eps^2 coefficient of sqrt(-g) * fun(I2,I3,J2,P3) using the series (exact to O(eps^2))."""
    e_sym = sp.symbols('E')
    vals = [s_[0] + e_sym * s_[1] + e_sym ** 2 * s_[2] for s_ in ser]
    sqv = sq_ser[0] + e_sym * sq_ser[1] + e_sym ** 2 * sq_ser[2]
    tot = sp.expand(sqv * fun(*vals))
    return sp.expand(sp.simplify(tot.coeff(e_sym, 2)))


def coeffs(L2):
    P = sp.Poly(L2, p0, p1, p2)
    g = P.coeff_monomial
    return dict(a1=g(p2 ** 2), a2=g(p1 ** 2), a3=g(p0 ** 2), a4=g(p1 * p2), a5=g(p0 * p2), a6=g(p0 * p1))


def reduce(c):
    alpha = sp.simplify(c['a1'])
    beta = sp.simplify(c['a2'] - sp.diff(c['a4'], r) / 2 - c['a5'])
    mu = sp.simplify(c['a3'] + sp.diff(c['a5'], r, 2) / 2 - sp.diff(c['a6'], r) / 2)
    return alpha, beta, mu


rsp = 1 / f0
# Einstein
e_sym = sp.symbols('E')
LE2 = sp.expand(sp.simplify(((sq_ser[0] + e_sym * sq_ser[1] + e_sym ** 2 * sq_ser[2]) *
                             (RsSer[0] + e_sym * RsSer[1] + e_sym ** 2 * RsSer[2] + 6)).expand().coeff(e_sym, 2)))
aE, bE, mE = reduce(coeffs(LE2))
print("Einstein: alpha =", aE, " beta =", sp.factor(bE), " mu =", mE)
eta0 = sp.limit(sp.simplify(2 * bE * rsp - sp.diff(2 * aE * sp.diff(rsp, r), r)), r, 1)
print("Einstein eta0 (flux units) =", eta0, ";  s0 =", s0)

results = {}
for name in (sys.argv[1:] or list(OPS)):
    fun = OPS[name]
    L2 = L2_of(fun)
    a, b, m = reduce(coeffs(L2))
    eta1 = sp.limit(sp.simplify(2 * b * rsp - sp.diff(2 * a * sp.diff(rsp, r), r)), r, 1)
    # Wald correction: dL/dd of fun(invariants of perturbed Riemann)
    s1 = sp.simplify((wald_fac * sp.diff(fun(*invP), d)).subs(d, 0).subs(r, 1))
    res = sp.nsimplify(sp.simplify(eta1 / eta0 - s1 / s0))
    results[name] = res
    print(f"[{name}] mu (expect 0 only after adding the background-EOM pieces; info) = {sp.simplify(m)}")
    print(f"[{name}] alpha near horizon: {sp.factor(a)}")
    print(f"[{name}] eta1/eta0 = {sp.simplify(eta1/eta0)},  s1/s0 = {sp.simplify(s1/s0)},  "
          f"4 pi eta/s - 1 = {res} * coupling")
print("\nSUMMARY (coefficient of the coupling lam in 4 pi eta/s - 1, L = R + 6 + lam O, L_AdS = 1):")
for k, v in results.items():
    print(f"  {k:7s}: {v}")
