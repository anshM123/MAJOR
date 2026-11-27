"""
Thermodynamic cross-check at O(lam) for L = R + 6 + lam O on the AdS4 black brane (L_AdS = 1).  [computed]

  (a) background correction from the reduced action in (N, f), temperature shift with horizon held at r = 1;
  (b) free energy at fixed T from the first-order on-shell action:  f = f0 (1 + k lam),  k = int_1^oo r^2 O dr
      (f0 = -1/(16 pi G) at r_h = 1); W^n (n >= 3) and (W^2)^2 do not change the AdS vacuum, C_T or L_AdS;
  (c) consistency: Wald s at r_h = 1 is s0 (1 + w lam); then at fixed T, s/s(Einstein) = 1 + (w - 2 tau1) lam,
      which must equal 1 + k lam  (first law, f ~ T^3).
"""
import sympy as sp

r, lam = sp.symbols('r lambda', positive=True)
t, x, y = sp.symbols('t x y')
X = [t, r, x, y]
N, f = sp.Function('N')(r), sp.Function('f')(r)
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def curv(g):
    ginv = g.inv()
    n = 4
    Gam = [[[sp.cancel(sum(ginv[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                           for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]
    Rl = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    v = sum(g[a, e] * (sp.diff(Gam[e][b][d], X[c]) - sp.diff(Gam[e][b][c], X[d])
                                       + sum(Gam[e][c][h] * Gam[h][b][d] - Gam[e][d][h] * Gam[h][b][c] for h in range(n)))
                            for e in range(n))
                    Rl[a, b, c, d] = sp.cancel(v)
    Ric = sp.Matrix(n, n, lambda b, d: sum(ginv[a, c] * Rl[a, b, c, d] for a in range(n) for c in range(n)))
    Rs = sp.cancel(sum(ginv[b, d] * Ric[b, d] for b in range(n) for d in range(n)))
    W = {}
    for (a, b, c, d), v in Rl.items():
        W[a, b, c, d] = (v - sp.Rational(1, 2) * (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c] - g[b, c] * Ric[a, d]
                                                  + g[b, d] * Ric[a, c]) + Rs / 6 * (g[a, c] * g[b, d] - g[a, d] * g[b, c]))
    M = sp.zeros(6, 6)
    for A, (a, b) in enumerate(PAIRS):
        for B, (c, d) in enumerate(PAIRS):
            M[A, B] = sp.cancel(sum(W[a, b, e, h] * ginv[e, c] * ginv[h, d] for e in range(4) for h in range(4)))
    return Rs, sp.cancel(4 * (M ** 2).trace()), sp.cancel(8 * (M ** 3).trace())


g = sp.diag(-N ** 2 * f, 1 / f, r ** 2, r ** 2)
Rs, I2, I3 = curv(g)
f0 = r ** 2 - 1 / r


def EL(Lag, F):
    out = 0
    for k in range(0, 5):
        Dk = F if k == 0 else sp.diff(F, r, k)
        out += (-1) ** k * sp.diff(sp.diff(Lag, Dk), r, k)
    return out


f1, n1 = sp.Function('f1')(r), sp.Function('n1')(r)
for name, O, wald in (('I3', I3, 6), ('I2^2', I2 ** 2, 96)):
    Lag = N * r ** 2 * (Rs + 6 + lam * O)
    sub = {N: 1 + lam * n1, f: f0 + lam * f1}
    eqs = []
    for F in (N, f):
        e = EL(Lag, F).subs(sub).doit()
        e1 = sp.simplify(sp.diff(e, lam).subs(lam, 0))
        eqs.append(e1)
    s_n = sp.dsolve(eqs[1], n1)            # EL_f involves n1 only at this order
    n1s = s_n.rhs.subs(sp.Symbol('C1'), sp.Symbol('cn'))
    e2 = sp.simplify(eqs[0].subs(n1, n1s).doit())
    s_f = sp.dsolve(e2, f1)
    f1s = s_f.rhs
    C = sorted((f1s.free_symbols | n1s.free_symbols) - {r, lam}, key=str)
    # horizon at r = 1 and boundary conformal to -dt^2 + dx^2
    conds = [sp.Eq(f1s.subs(r, 1), 0), sp.Eq(sp.limit((2 * n1s * f0 + f1s) / r ** 2, r, sp.oo), 0)]
    sol = sp.solve(conds, C, dict=True)[0]
    f1s, n1s = sp.simplify(f1s.subs(sol)), sp.simplify(n1s.subs(sol))
    tau1 = sp.simplify(sp.limit((n1s * sp.diff(f0, r) + sp.diff(f1s, r)) / sp.diff(f0, r), r, 1))
    k = sp.integrate(r ** 2 * O.subs({N: 1, f: f0}).doit(), (r, 1, sp.oo))
    print(f"[{name}] f1 = {f1s},  n1 = {n1s}")
    print(f"[{name}] T = T0 (1 + {tau1} lam);  free-energy k = {sp.simplify(k)};  Wald w = {wald};  "
          f"first-law check w - 2 tau1 = {sp.simplify(wald - 2 * tau1)}  (should equal k)")
