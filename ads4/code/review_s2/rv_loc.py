"""rv_loc.py [exact/Arb]: Lemma Loc'' tail ratio.  For every unit nu-cell used with Rt = 1 by loc4.log_loc4 (hypothesis
q = 4 nu_lo^2 - 2 nu_hi^2 - 3/2 > 0, J^2 q >= (1/4 + nu_hi^2)^2), check in Arb that
  sup_{j >= J, nu in cell} |c_{j+1}/c_j| <= ((x - 1/2)^2 + nu_hi^2)/(x sqrt(x^2 + 4 nu_lo^2)) < 1   (x = j + 1),
on x = J+1 ... J+20000 (and x -> oo: the bound -> 1 from below since P(x) < 0), plus the polynomial identity behind the proof."""
import math
from flint import arb, ctx, fmpq_poly, fmpq
ctx.prec = 128
# polynomial identity: ((x-1/2)^2 + h)^2 - x^2 (x^2 + 4 l) == -2x^3 - q x^2 - 2x(1/4 + h) + (1/4 + h)^2 with q = 4l - 2h - 3/2 (h = nu_hi^2, l = nu_lo^2)
import sympy as sp
x, h, l = sp.symbols('x h l')
lhs = sp.expand(((x - sp.Rational(1, 2))**2 + h)**2 - x**2*(x**2 + 4*l))
q = 4*l - 2*h - sp.Rational(3, 2)
rhs = sp.expand(-2*x**3 - q*x**2 - 2*x*(sp.Rational(1, 4) + h) + (sp.Rational(1, 4) + h)**2)
print('identity P(x):', sp.simplify(lhs - rhs) == 0)
worst = 0.0; ncell = 0
for J in (1500, 2200):
    for T in list(range(1, 1150)) + [0.4]:
        lo, hi = (0.4, 1.0) if T == 0.4 else (float(T), float(T + 1))
        cc = 4*lo*lo - 2*hi*hi - 1.5
        if not (cc > 0 and J*J*cc >= (0.25 + hi*hi)**2):
            continue
        ncell += 1
        L2, H2 = arb(lo)**2, arb(hi)**2
        for xx in list(range(J + 1, J + 200)) + list(range(J + 200, J + 20000, 97)):
            X = arb(xx)
            r = ((X - arb(0.5))**2 + H2)/(X*(X*X + 4*L2).sqrt())
            worst = max(worst, float(r.upper()))
print(f'cells using Rt = 1: {ncell}; max ratio bound over tested j >= J: {worst:.6f} (< 1 required)')
