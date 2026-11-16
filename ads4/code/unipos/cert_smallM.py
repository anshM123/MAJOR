r"""Certificate C2 (exact rational arithmetic):  for M <= M1 = 23/20 the function
        h_M(p) = phi(p) * p * coth(pi M p),      phi(p) = (1-p)^3 (65 + 155 p + 47 p^2),
is strictly decreasing on (0,1).

d/dp log h_M = phi'/phi + pi M l(pi M p),  l(x) = 1/x - 2/sinh(2x) < 2x/3  (x > 0),  and
-phi'/phi = (40 + 526 p + 235 p^2) / ((1-p)(65 + 155 p + 47 p^2)).   Hence it suffices that
        D_c(p) := 40 + 526 p + 235 p^2 - c p (1-p) (65 + 155 p + 47 p^2) > 0   on [0,1]
for one rational c >= 2 pi^2 M1^2 / 3.  We take c = 871/100 and check
  (i)  2 pi^2 M1^2/3 < c                         (interval arithmetic for pi),
  (ii) D_c has no real root in [0,1] and D_c(0) > 0   (Sturm sequence, exact, via sympy),
  (iii) independent check of (ii): interval evaluation of D_c on 4096 subintervals of [0,1].
"""
import sys
from fractions import Fraction as Fr
import sympy as sp
from common import *

c = Fr(871, 100)
ok1 = (2 * iv.pi * iv.pi * frac_iv(M1) ** 2 / 3).b < frac_iv(c).a
print(f"{'PASS' if ok1 else 'FAIL'} (i)  2 pi^2 M1^2/3 = {mp.nstr(mp.mpf((2*iv.pi*iv.pi*frac_iv(M1)**2/3).b), 12)} < c = {c}")

p = sp.symbols('p')
phi = (1 - p) ** 3 * (65 + 155 * p + 47 * p ** 2)
# exact identity  -phi' = (1-p)^2 (40 + 526 p + 235 p^2)
ok0 = sp.expand(-sp.diff(phi, p) - (1 - p) ** 2 * (40 + 526 * p + 235 * p ** 2)) == 0
print(f"{'PASS' if ok0 else 'FAIL'} (0)  phi' = -(1-p)^2 (40 + 526 p + 235 p^2)")
D = sp.Poly(sp.expand(40 + 526 * p + 235 * p ** 2 - sp.Rational(c.numerator, c.denominator) * p * (1 - p) * (65 + 155 * p + 47 * p ** 2)), p)
nroots = D.count_roots(0, 1)
ok2 = (nroots == 0) and D.eval(0) > 0
print(f"{'PASS' if ok2 else 'FAIL'} (ii) D_c = {D.as_expr()} : real roots in [0,1] = {nroots}, D_c(0) = {D.eval(0)}")

coef = [Fr(int(sp.numer(a)), int(sp.denom(a))) for a in D.all_coeffs()]      # highest degree first
n = 4096
dmin = None
ok3 = True
for i in range(n):
    x = iv.mpf([i, i + 1]) / n
    v = iv.mpf(0)
    for a in coef:
        v = v * x + frac_iv(a)
    # Horner on a subinterval of width 1/4096 (rigorous enclosure, slightly wide)
    ok3 = ok3 and v.a > 0
    dmin = lo(v) if dmin is None else min(dmin, lo(v))
print(f"{'PASS' if ok3 else 'FAIL'} (iii) interval check on {n} subintervals: D_c >= {mp.nstr(dmin, 6)} on [0,1]")
ok = ok0 and ok1 and ok2 and ok3
print(f"VERDICT C2 (h_M strictly decreasing on (0,1) for all 0 < M <= {M1}): {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
