"""Checks of ANGULAR_ROUTE §6 (author claims): generating function S(z), f_0, Y_0 = Einstein, Y_4 = 303/5120, y_m formula and sign,
min_x (x+998)^6/x^4, the 47266 constant, F(phi) ~ -(4 sqrt2/3)|phi|^{3/2}."""
import sympy as sp, mpmath as mp
from fractions import Fraction as F
z = sp.symbols('z')
f = lambda k: sp.rf(sp.Rational(1, 2), k-1)/(sp.factorial(k-1)*k**2)
S = 4*(1 - sp.sqrt(1-z) + sp.log((1 + sp.sqrt(1-z))/2))
ser = sp.series(S, z, 0, 12).removeO()
assert all(sp.simplify(ser.coeff(z, k) - f(k)) == 0 for k in range(1, 12)); print("S(z) generating function OK (m<12)")
f0 = -8*(1 - sp.log(2))
fa = lambda k: f0 if k == 0 else (f(abs(k)) if k != 0 else f0)
Y = lambda m: sp.Rational(12 - m*m, 2)*fa(m) + sp.Rational(13 - 7*m, 2)*fa(m-1) + sp.Rational(13 + 7*m, 2)*fa(m+1) + sp.Rational((m-2)**2, 4)*fa(m-2) + sp.Rational((m+2)**2, 4)*fa(m+2)
print("Y_0 =", sp.nsimplify(sp.expand(Y(0))), " claim -139/4 + 48 log2:", sp.simplify(Y(0) - (-sp.Rational(139, 4) + 48*sp.log(2))) == 0)
print("Y_4 =", Y(4), " claim 303/5120")
for m in range(3, 80):
    ym = Y(m); g = sp.Rational(3*(48*m**5 - 176*m**4 + 303*m**3 - 376*m**2 + 297*m - 480), 16*(m-1)*(m+1)**2*(2*m-5)*(2*m-3))
    assert sp.simplify(ym/f(m) - g) == 0 and ym > 0, m
mm = sp.symbols('m'); print("y_m/f_m formula and y_m>0 OK for 3<=m<80; numerator at m=n+3:", sp.expand((48*mm**5-176*mm**4+303*mm**3-376*mm**2+297*mm-480).subs(mm, mm+3)), "; y_2 =", Y(2))
x = sp.symbols('x', positive=True); g = (x + 998)**6/x**4; xs = sp.solve(sp.diff(sp.log(g), x), x); print("argmin", xs, " min == 729/16*998^2:", sp.simplify(g.subs(x, xs[0]) - sp.Rational(729, 16)*998**2) == 0)
print("constant:", float(sp.Rational(5, 4800)*sp.Rational(29999, 30000)**3*sp.Rational(729, 16)*998**2), " (claim > 47266)")
mp.mp.dps = 30
Sf = lambda w: 4*(1 - mp.sqrt(1-w) + mp.log((1 + mp.sqrt(1-w))/2))
Fphi = lambda p: 2*mp.re(Sf(mp.exp(1j*p)) - Sf(1))
for p in [1e-2, 1e-4, 1e-6]: print("F(phi)/(-(4sqrt2/3)phi^1.5) at", p, mp.nstr(Fphi(p)/(-(4*mp.sqrt(2)/3)*p**1.5), 8))
print("bound: sqrt(t4^2+aT^2) =", float(4*sp.sqrt((sp.Rational(139, 4) - 48*sp.log(2))/(47266 + sp.Rational(303, 5120)))))
