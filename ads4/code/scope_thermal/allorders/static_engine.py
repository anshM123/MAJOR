"""
Static shear-channel engine (diagonal gauge), generalising ../etas_diag.py to
  * arbitrary background functions N(r), f(r) (for the O(lam^2) computation of item 5),
  * radial profiles theta(r) multiplying operators,
  * operators containing derivatives of curvature: (grad I)^2 and U = grad C . grad C (third derivatives of phi),
  * the generic membrane flux for L2 = sum_k alpha_k (phi^(k))^2 (k <= 3).   [computed]

Metric: ds^2 = -N^2 f dt^2 + dr^2/f + r^2 (e^{eps phi} dx^2 + e^{-eps phi} dy^2), horizon at r = 1.
Flux:   Pi[chi] = sum_{k>=1} (-1)^{k-1} d^{k-1}( 2 alpha_k chi^(k) ),  eta ~ -Pi[-r_*]|_{r->1},  r_*' = 1/(N f).
(The sign/normalisation is fixed by Einstein gravity: 4 pi eta/s = 1.)
"""
import itertools
import sympy as sp

r, eps = sp.symbols('r epsilon', positive=True)
t, x, y = sp.symbols('t x y')
X = [t, r, x, y]
phi = sp.Function('phi')(r)
P = sp.symbols('p0:5')   # jets phi, phi', phi'', phi''', phi''''
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def jets(e):
    for k in (4, 3, 2, 1):
        e = e.subs(sp.Derivative(phi, (r, k)), P[k])
    return e.subs(phi, P[0])


def metric(N, f, pert=True):
    e = sp.exp(eps * phi) if pert else 1
    return sp.diag(-N ** 2 * f, 1 / f, r ** 2 * e, r ** 2 / e)


def curvature(g):
    n = 4
    gi = sp.diag(*[1 / g[i, i] for i in range(n)])
    Gam = [[[sp.cancel((gi[a, a] * (sp.diff(g[a, b], X[c]) + sp.diff(g[a, c], X[b]) - sp.diff(g[b, c], X[a]))) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    Rl = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        if a == b or c == d:
            Rl[a, b, c, d] = 0
            continue
        v = g[a, a] * (sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
                       + sum(Gam[a][c][h] * Gam[h][b][d] - Gam[a][d][h] * Gam[h][b][c] for h in range(n)))
        Rl[a, b, c, d] = sp.cancel(v)
    Ric = sp.Matrix(n, n, lambda b, d: sum(gi[a, a] * Rl[a, b, a, d] for a in range(n)))
    Rs = sp.cancel(sum(gi[b, b] * Ric[b, b] for b in range(n)))
    W = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        W[a, b, c, d] = sp.cancel(Rl[a, b, c, d] - sp.Rational(1, 2) * (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c]
                                                                         - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c])
                                  + Rs / 6 * (g[a, c] * g[b, d] - g[a, d] * g[b, c]))
    return gi, Gam, Rl, Ric, Rs, W


def invariants(g, gi, W):
    """I = W.W, J = W_ab^cd W_cd^ef W_ef^ab (diagonal metric => explicit sums)."""
    n = 4
    Wm = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        Wm[a, b, c, d] = W[a, b, c, d] * gi[c, c] * gi[d, d]
    I = sum(Wm[a, b, c, d] * Wm[c, d, a, b] for a, b, c, d in itertools.product(range(n), repeat=4))
    J = 0
    for a, b, c, d in itertools.product(range(n), repeat=4):
        if Wm[a, b, c, d] == 0:
            continue
        for e, h in itertools.product(range(n), repeat=2):
            J += Wm[a, b, c, d] * Wm[c, d, e, h] * Wm[e, h, a, b]
    return I, J


def gradC2(g, gi, Gam, W):
    """U = nabla_e C_abcd nabla^e C^abcd for a static diagonal metric."""
    n = 4
    U = 0
    for e in range(n):
        for a, b, c, d in itertools.product(range(n), repeat=4):
            v = sp.diff(W[a, b, c, d], X[e])
            for h in range(n):
                v -= (Gam[h][e][a] * W[h, b, c, d] + Gam[h][e][b] * W[a, h, c, d]
                      + Gam[h][e][c] * W[a, b, h, d] + Gam[h][e][d] * W[a, b, c, h])
            if v == 0:
                continue
            U += v ** 2 * gi[e, e] * gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d]
    return U


def box(g, gi, s):
    """Scalar Laplacian of a function s(r) on a static diagonal metric."""
    sq = sp.sqrt(-g[0, 0] * g[1, 1] * g[2, 2] * g[3, 3])
    return sp.diff(sq * gi[1, 1] * sp.diff(s, r), r) / sq


def eps_series(e, order=2):
    e = jets(e)
    return [sp.simplify(sp.diff(e, eps, k).subs(eps, 0) / sp.factorial(k)) for k in range(order + 1)]


def reduce_quadratic(L2):
    """Reduce sum c_ab p_a p_b (a,b <= 4) by integrations by parts to sum_k alpha_k p_k^2; returns alpha list, mass."""
    L2 = sp.expand(L2)
    C = {}
    for a in range(5):
        for b in range(a, 5):
            if a == b:
                c = L2.coeff(P[a], 2)
                for o in range(5):
                    if o != a:
                        c = c.subs(P[o], 0)
            else:
                c = L2.coeff(P[a], 1).coeff(P[b], 1)
            c = sp.simplify(c)
            if c != 0:
                C[(a, b)] = c
    alpha = {k: 0 for k in range(5)}
    work = dict(C)
    while work:
        (a, b), c = next(iter(work.items()))
        del work[(a, b)]
        if a == b:
            alpha[a] += c
        elif b == a + 1:
            alpha[a] += -sp.diff(c, r) / 2
        else:   # c p_a p_b ~ -c' p_a p_{b-1} - c p_{a+1} p_{b-1}
            for key, val in (((a, b - 1), -sp.diff(c, r)), ((a + 1, b - 1), -c)):
                k2 = tuple(sorted(key))
                work[k2] = work.get(k2, 0) + val
    return {k: sp.simplify(v) for k, v in alpha.items() if v != 0}


def flux_at_horizon(alpha, rsp):
    """-Pi[-r_*] at r -> 1 with Pi = sum_{k>=1} (-1)^{k-1} d^{k-1}(2 alpha_k chi^(k)), chi' = -r_*'."""
    chi_der = {1: -rsp}
    for k in range(2, 5):
        chi_der[k] = sp.diff(chi_der[k - 1], r)
    Pi = 0
    for k, a in alpha.items():
        if k == 0:
            continue
        term = 2 * a * chi_der[k]
        for _ in range(k - 1):
            term = sp.diff(term, r)
        Pi += (-1) ** (k - 1) * term
    return sp.limit(sp.simplify(-Pi), r, 1)


def wald_s(Lfun_of_R, g, gi, Rl):
    """Wald entropy density (units 1/16piG) at r -> 1 for a Lagrangian given as a function of the lower Riemann dict.
    Valid for algebraic L(Riemann); for derivative operators use the first-law route (thermo_* scripts)."""
    dd = sp.Symbol('dd')
    Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
    Rp = {k: Rl[k] + dd * Xt.get(k, 0) for k in Rl}
    val = sp.diff(Lfun_of_R(Rp), dd).subs(dd, 0)
    s = -8 * sp.pi * (val / 4) * (-g[0, 0] * g[1, 1]) * r ** 2
    return sp.limit(sp.simplify(s), r, 1)


def weyl_from(Rl, g, gi):
    n = 4
    Ric = sp.Matrix(n, n, lambda b, d: sum(gi[a, a] * Rl[a, b, a, d] for a in range(n)))
    Rs = sum(gi[b, b] * Ric[b, b] for b in range(n))
    W = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        W[a, b, c, d] = (Rl[a, b, c, d] - sp.Rational(1, 2) * (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c]
                                                              - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c])
                         + Rs / 6 * (g[a, c] * g[b, d] - g[a, d] * g[b, c]))
    return W, Rs


def Dr(e):
    """total r-derivative of an expression in r and jets p0..p3."""
    return sp.diff(e, r) + sum(P[k + 1] * sp.diff(e, P[k]) for k in range(4))
