"""Checks of the ingredients of Lemma B (CERT_CLOSE.md §1.3) -- symbolic (sympy) and high-precision numeric (mpmath):
  (1) y(p) = Y(1 - 2 u p^2) with Y a solution of the Jacobi (0,8) equation satisfies (K y')' + 4 u lam rho y = 0,
      K = p (1-up^2)^9, rho = p (1-up^2)^8, lam = n (n+9)                                     [symbolic, generic n, u]
  (2) L f := d/dp[(f' - f/p)(1-up^2) + 16 u p f] acts on monomials as L p^m = (m-1)^2 p^(m-2) + (17-m)(m+1) u p^m   [symbolic]
  (3) (g' K) = (f' - f/p)(1-up^2) + 16 u p f  for g = f/rho                                     [symbolic]
  (4) the three-fold identity  int w_B y = -(4ul)^-3 [ beta_{f2} + int (L^3 w_B) y ]  for NF02 (pattern s2=s3=-1, delta=1/5),
      J = 60, 201 and u = 0.5, 0.9, 1, by 60-digit quadrature (Jacobi polynomials from mpmath).                [numeric]"""
import sys, os
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
p, u, n, x = sp.symbols('p u n x', positive=True)
Y = sp.Function('Y')
lam = n * (n + 9)
K = p * (1 - u * p**2) ** 9; rho = p * (1 - u * p**2) ** 8
X = 1 - 2 * u * p**2
y = Y(X)
expr = sp.diff(K * sp.diff(y, p), p) + 4 * u * lam * rho * y
Y0, Y1, Y2 = sp.symbols('Y0 Y1 Y2')
e = expr.doit()
rep = {}
for s_ in e.atoms(sp.Subs):
    order = sum(c for _, c in s_.expr.variable_count) if isinstance(s_.expr, sp.Derivative) else 0
    rep[s_] = {1: Y1, 2: Y2}[order]
e = e.xreplace(rep).xreplace({Y(X): Y0})
jac = ((1 - x**2) * Y2 + (8 - 10 * x) * Y1 + lam * Y0).subs(x, X)      # Jacobi (alpha, beta) = (0, 8), DLMF 18.8.1
print('(1) (K y\')\' + 4 u lam rho y  -  4 p u (1-up^2)^8 * [Jacobi operator at x = 1-2up^2]  =', sp.simplify(sp.expand(e - 4 * p * u * (1 - u * p**2)**8 * jac)))
m = sp.Symbol('m')
f = p**m
Lf = sp.diff((sp.diff(f, p) - f / p) * (1 - u * p**2) + 16 * u * p * f, p)
print('(2) L p^m - formula:', sp.simplify(Lf - ((m - 1)**2 * p**(m - 2) + (17 - m) * (m + 1) * u * p**m)))
F = sp.Function('f')(p)
g = F / rho
print('(3) g\'K - formula:', sp.simplify(sp.diff(g, p) * K - ((sp.diff(F, p) - F / p) * (1 - u * p**2) + 16 * u * p * F)))
# (4) numeric
import nf_cert, rcert, nf_corner
from nf_corner import P2, smoothstep
mp.mp.dps = 60
F0 = nf_cert.load('NF02')
g = rcert.gens(*F0, Fr(1, 5), '+-', 1)[3]
a, c = Fr(1, 5), Fr(3, 5)
S = smoothstep(a, c, 7)
w = P2.from_p(g.W0) + P2.from_p(g.W1, 1)
pieces0 = [(a, c, w * P2.from_p(S)), (c, Fr(1), w)]
def ev(P, pp, uu):
    return mp.fsum(mp.mpf(v.numerator) / v.denominator * pp**i * uu**j for (i, j), v in P.c.items())
for J in (60, 201):
    nn = J - 4; lamv = nn * (nn + 9)
    for uu in (mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf(1)):
        yv = lambda pp: mp.jacobi(nn, 0, 8, 1 - 2 * uu * pp**2)
        lhs = mp.fsum(mp.quad(lambda pp: ev(P, pp, uu) * yv(pp), mp.linspace(float(lo), float(hi), 40)) for lo, hi, P in pieces0)
        f3 = [(lo, hi, P.L().L().L()) for lo, hi, P in pieces0]
        f2hi = pieces0[1][2].L().L()
        dd = f2hi.dp()
        beta = -ev(dd, mp.mpf(1), uu) * (1 - uu) * yv(mp.mpf(1))
        rhs = -(4 * uu * lamv) ** -3 * (beta + mp.fsum(mp.quad(lambda pp: ev(P, pp, uu) * yv(pp), mp.linspace(float(lo), float(hi), 40)) for lo, hi, P in f3))
        print('(4) J=%d u=%s: int w_B y = %s ; IBP^3 = %s ; rel diff %.1e' % (J, mp.nstr(uu, 3), mp.nstr(lhs, 15), mp.nstr(rhs, 15), float(abs(lhs - rhs) / abs(lhs))), flush=True)
