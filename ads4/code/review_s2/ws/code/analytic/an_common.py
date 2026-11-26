"""
analytic/an_common.py: shared exact evaluators for the ANALYTIC_SLIVER study.

State (J, a, b): n = J - b, s = 2b+1, N = n + a, C0 = n(n+s), D = N(N+s), c = C0/D, d = 1-c, h = 1/D.
Ladder U_k = prod_{j<k} (C0 - j(j+s))/(D - j(j+s)), terminating at k = n+1.
Taylor coefficients at z = 0 of F_{p,nu}(z) = (1-z)^{-p} P_{-1/2+i nu}((1+z)/(1-z)) from the ODE
    z(1-z)^2 F'' + (1 - (2+2p) z + (1+2p) z^2) F' + (nu^2 + 1/4 - p + p^2 z) F = 0:
    (k+1)^2 g_{k+1} = (2k^2 + 2pk + p - nu^2 - 1/4) g_k - (k-1+p)^2 g_{k-1},  g_0 = 1.
K_p(nu) = L(F_{p,nu}) = sum_{k<=n} U_k g_k(nu)  (exact, high precision; used only as [DV] reference).
"""
import os, sys, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, acb, ctx, arb_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
CACHE = os.path.join(HERE, 'cache')
K0 = Fr('27.48509817')
LAM2 = 998 * 999


def nf02():
    nf = json.load(open(os.path.join(ADS4, 'prl', 'data', 'nf02_psi.json')))
    return [Fr(c) / K0 for c in nf['psi2_exact']], [Fr(c) / K0 for c in nf['psi3_exact']]


def state(J, a, b):
    n = J - b
    s = 2 * b + 1
    N = n + arb(a)
    C0 = arb(n * (n + s))
    D = N * (N + s)
    c = C0 / D
    return dict(J=J, a=a, b=b, n=n, s=s, N=N, C0=C0, D=D, c=c, d=1 - c, h=1 / D)


def ladder(st):
    n, s, C0, D = st['n'], st['s'], st['C0'], st['D']
    U = [arb(1)]
    for j in range(n):
        p = j * (j + s)
        U.append(U[-1] * (C0 - p) / (D - p))
    return U  # U_0..U_n (U_{n+1} = 0)


def gcoef(p, nu, K):
    """g_0..g_K for F_{p,nu}, nu real (arb) or acb."""
    nu2 = nu * nu
    g = [arb(1) if not isinstance(nu, acb) else acb(1)]
    # k = 0: g_1 = (p - nu^2 - 1/4) g_0
    g.append((p - nu2 - arb(1) / 4) * g[0])
    for k in range(1, K):
        g.append(((2 * k * k + 2 * p * k + p - nu2 - arb(1) / 4) * g[k] - (k - 1 + p) ** 2 * g[k - 1]) / ((k + 1) ** 2))
    return g


def Kladder(p, nu, U):
    g = gcoef(p, nu, len(U) - 1)
    t = 0
    for k in range(len(U)):
        t += U[k] * g[k]
    return t


def Finf(p, nu, c):
    """F_{p,nu}(c) via Arb Legendre."""
    c = arb(c)
    X = acb((1 + c) / (1 - c))
    mu = acb(-0.5, nu)
    return ((1 - c) ** (-p) * X.legendre_p(mu, 0, 3)).real
