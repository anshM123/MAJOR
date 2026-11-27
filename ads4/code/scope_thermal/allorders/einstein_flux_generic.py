"""Item 5 support: Einstein-part horizon flux for arbitrary N, f (local jets).  [computed]"""
import sympy as sp
from static_engine import *
# ---- (i) Einstein flux / Wald for arbitrary N, f near the horizon (local jets)
a = sp.symbols('a1:5'); b = sp.symbols('b0:4')
fl = sum(a[k - 1] * (r - 1) ** k for k in range(1, 5))
Nl = sum(b[k] * (r - 1) ** k for k in range(4))
gG = metric(Nl, fl)
gi, Gam, Rl, Ric, Rs, W = curvature(gG)
RsS = eps_series(Rs)
alE = reduce_quadratic(Nl * r ** 2 * RsS[2])
print("(i) Einstein part on arbitrary background: alphas =", {k: sp.factor(v) for k, v in alE.items()})
print("    -Pi[-r*] at horizon =", sp.simplify(flux_at_horizon(alE, 1 / (Nl * fl))), " (expect -1 * r_h^2 = -1)")

