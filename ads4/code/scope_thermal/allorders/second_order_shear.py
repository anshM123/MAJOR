"""
Item 5 (shear): O(lam^2) term of 4 pi eta/s for L = R + 6 + lam J  (J = C^3, lam = -alpha3/6).        [computed]
* Einstein kinetic term is exactly -N r^2 f/2 for any static N, f (einstein_alpha_generic.py), so its horizon flux is
  r_h^2 independently of the background; the Einstein Wald entropy is 4 pi r_h^2 likewise. Hence only the J part on the
  O(lam)-corrected brane is needed:  f = f0 + lam(-8/r + 24/r^4 - 16/r^7),  N = 1 - 6 lam/r^6  (../thermo_check.py).
* J-part Wald entropy on that brane: w = 6 + 312 lam (wald_J_corrected.py).
"""
import sympy as sp
from static_engine import *

lam = sp.Symbol('lam')
f0 = r ** 2 - 1 / r
N = 1 + lam * (-6 / r ** 6)
f = f0 + lam * (-8 / r + 24 / r ** 4 - 16 / r ** 7)
g = metric(N, f)
gi, Gam, Rl, Ric, Rs, W = curvature(g)
I, J = invariants(g, gi, W)
Jj = jets(J)
J2 = sp.diff(Jj, eps, 2).subs(eps, 0) / 2
L2 = N * r ** 2 * J2
L2 = sp.expand(L2.subs(lam, 0) + lam * sp.diff(L2, lam).subs(lam, 0))     # O(lam) is all that is needed
al = reduce_quadratic(L2)
rsp = 1 / f0 - lam * (n1 := -6 / r ** 6) / f0 - lam * (-8 / r + 24 / r ** 4 - 16 / r ** 7) / f0 ** 2   # 1/(N f) to O(lam)
chi_der = {1: -rsp}
for k in range(2, 5):
    chi_der[k] = sp.diff(chi_der[k - 1], r)
Pi = 0
for k, a in al.items():
    if k == 0:
        continue
    term = 2 * a * chi_der[k]
    for _ in range(k - 1):
        term = sp.diff(term, r)
    Pi += (-1) ** (k - 1) * term
mPi = sp.expand(-Pi)
c0 = sp.limit(sp.cancel(sp.together(mPi.subs(lam, 0))), r, 1)
c1 = sp.limit(sp.cancel(sp.together(sp.diff(mPi, lam).subs(lam, 0))), r, 1)
print("J-part horizon flux (Einstein = -1):", c0, "+ lam *", c1)
eta_hat = 1 - lam * (c0 + lam * c1)          # eta / eta_E at fixed r_h
s_hat = 1 + lam * (6 + 312 * lam)
ratio = sp.expand(sp.series(eta_hat / s_hat, lam, 0, 3).removeO())
a3, t4 = sp.symbols('alpha3 t4')
print("4 pi eta/s =", ratio)
print("           =", sp.expand(ratio.subs(lam, -a3 / 6)), " = ", sp.expand(ratio.subs(lam, t4 / 1440)),
      "   [Li-Lu-Mai: 1 - 3 t4/40 - 17 t4^2/3840 + ...]")
