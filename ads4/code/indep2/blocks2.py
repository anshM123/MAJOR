"""Independent (0, b) continuous-spin blocks in d = 3, generic external dimension shift parameter b, generic Delta, complex l.

Definitions reused from the previous independent check (ads4/code/indep/cblocks.py docstring, itself an own derivation of the
Dolan-Osborn 1/2-Casimir in radial variables); the implementation below is new and works over any field
(sympy expressions with exact cancellation, or mpmath mpc numbers).

   (1-r^2) Den (D_r^2 + Lam - C) G - [2(1-r^2)^2((1+r^2) + 4 b r eta) + Den (1+r^2)] D_r G
        - 8 r (1-r^2) [r eta + b(1+r^2)] (eta^2-1) d_eta G = 0,
   Den = (1+r^2)^2 - 4 r^2 eta^2,  D_r = r d_r,  Lam P_mu = mu(mu+1) P_mu,  C = Delta(Delta-3) + l(l+1).
   G = r^Delta sum_N r^N sum_m B[N][m] P_{l+m}(eta),  B[0][0] = 1.
"""
import sympy as sp

_r, _e, _b = sp.symbols('r e b')


def operator_terms():
    """[(kind, k, j, coeff(b) as sympy expr)] for the r^k eta^j monomials of C1, C2, C3 (kind 1, 2, 3)."""
    Den = (1 + _r**2)**2 - 4 * _r**2 * _e**2
    C1 = sp.expand((1 - _r**2) * Den)
    C2 = sp.expand(-(2 * (1 - _r**2)**2 * ((1 + _r**2) + 4 * _b * _r * _e) + Den * (1 + _r**2)))
    C3 = sp.expand(-8 * _r * (1 - _r**2) * (_r * _e + _b * (1 + _r**2)))
    out = []
    for kind, Cx in ((1, C1), (2, C2), (3, C3)):
        P = sp.Poly(Cx, _r, _e)
        for (k, j), c in zip(P.monoms(), P.coeffs()):
            out.append((kind, int(k), int(j), sp.expand(c)))
    return out


_TERMS = operator_terms()


def recursion(l, b, Delta, Nmax, conv, simp=lambda x: x):
    """B[N][m] for N = 0..Nmax.  conv(sympy expr in b) -> field element; simp: simplifier applied to each new coefficient."""
    terms = [(kind, k, j, conv(c.subs(_b, b) if not isinstance(b, sp.Symbol) else c)) for (kind, k, j, c) in _TERMS]
    C = Delta * (Delta - 3) + l * (l + 1)
    B = [dict() for _ in range(Nmax + 1)]
    B[0][0] = conv(sp.Integer(1))

    def eta_mult(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            out[m + 1] = out.get(m + 1, 0) + c * (mu + 1) / (2 * mu + 1)
            out[m - 1] = out.get(m - 1, 0) + c * mu / (2 * mu + 1)
        return out

    def D1(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            f = mu * (mu + 1) / (2 * mu + 1)
            out[m + 1] = out.get(m + 1, 0) + c * f
            out[m - 1] = out.get(m - 1, 0) - c * f
        return out

    for N in range(1, Nmax + 1):
        src = {}
        for (kind, k, j, cf) in terms:
            if k == 0 or k > N or not B[N - k]:
                continue
            g = B[N - k]
            D = Delta + N - k
            if kind == 1:
                h = {m: c * (D * D - C + (l + m) * (l + m + 1)) for m, c in g.items()}
            elif kind == 2:
                h = {m: c * D for m, c in g.items()}
            else:
                h = D1(g)
            for _ in range(j):
                h = eta_mult(h)
            for m, c in h.items():
                src[m] = src.get(m, 0) + cf * c
        D = Delta + N
        for m, c in src.items():
            c = simp(c)
            if c == 0:
                continue
            E = D * D - 3 * D + (l + m) * (l + m + 1) - C
            B[N][m] = simp(-c / E)
    return B


def eta_mult_basis(g, l):
    """multiply sum_m g[m] P_{l+m} by eta (same basis)."""
    out = {}
    for m, c in g.items():
        mu = l + m
        out[m + 1] = out.get(m + 1, 0) + c * (mu + 1) / (2 * mu + 1)
        out[m - 1] = out.get(m - 1, 0) + c * mu / (2 * mu + 1)
    return out


def magic_residual(l, b, Nmax, conv, simp=lambda x: x, cplus=None, cminus=None, Dl=1, Dr=0):
    """Coefficients (in the basis r^N P_{l+m}) of
        2 eta (1+r^2) Ghat^{Dl}_b(l)  -  cm * Ghat^{Dr}_b(l-1)  -  cp * Ghat^{Dr}_b(l+1),
    cm = 2l/(2l+1) [= 1 + i/(2 nu)], cp = 2(l+1)/(2l+1) [= 1 - i/(2 nu)] unless overridden.
    Returns list over N of dicts m -> residual, and the list of the LHS coefficients for scale."""
    if cminus is None:
        cminus = 2 * l / (2 * l + 1)
    if cplus is None:
        cplus = 2 * (l + 1) / (2 * l + 1)
    B1 = recursion(l, b, Dl, Nmax, conv, simp)
    Bm = recursion(l - 1, b, Dr, Nmax, conv, simp)
    Bp = recursion(l + 1, b, Dr, Nmax, conv, simp)
    lhs = [dict() for _ in range(Nmax + 1)]
    for N in range(Nmax + 1):
        e1 = eta_mult_basis(B1[N], l)
        for m, c in e1.items():
            lhs[N][m] = lhs[N].get(m, 0) + 2 * c
            if N + 2 <= Nmax:
                lhs[N + 2][m] = lhs[N + 2].get(m, 0) + 2 * c
    res = [dict(d) for d in lhs]
    for N in range(Nmax + 1):
        for m, c in Bm[N].items():           # P_{(l-1)+m} = P_{l+(m-1)}
            res[N][m - 1] = res[N].get(m - 1, 0) - cminus * c
        for m, c in Bp[N].items():           # P_{(l+1)+m} = P_{l+(m+1)}
            res[N][m + 1] = res[N].get(m + 1, 0) - cplus * c
    res = [{m: simp(c) for m, c in d.items()} for d in res]
    return res, lhs, (B1, Bm, Bp)
