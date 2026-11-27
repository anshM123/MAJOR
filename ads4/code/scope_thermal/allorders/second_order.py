"""
Item 5: O(coupling^2) terms in eta/s and eta_H/s for L = R + 6 + lam J + kap Q  (J = C^3, Q = C C C~),
lam = -alpha3/6, kap = -alpha3~/6, L_AdS = 1, orientation eps_{vrxy} = +r^2.                    [computed]

Structure (argued in README, checked here where computable):
  * eta: lam^2 term computed below; NO kap^2 term (Q is linear in kap, its quadratic action has only phi-chi cross
    terms, which feed back into the chi channel only at O(omega^2); Q does not change the parity-even background or
    the Wald entropy at any order in kap alone).
  * eta_H: lam*kap term computed below; no lam^2 or kap^2 terms (parity).
Method: Einstein-part contributions to the horizon flux and to the Wald entropy are background independent at fixed
r_h (checked below for arbitrary N, f); the J and Q parts are evaluated on the O(lam)-corrected brane
  f = f0 + lam (-8/r + 24/r^4 - 16/r^7),  N = 1 - 6 lam / r^6   (../thermo_check.py; r_h = 1, boundary normalised).
"""
import sympy as sp
from static_engine import *

lam = sp.Symbol('lam')
f0 = r ** 2 - 1 / r
f1 = -8 / r + 24 / r ** 4 - 16 / r ** 7
n1 = -6 / r ** 6

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

# ---- (ii) J on the corrected background, to O(lam)
N = 1 + lam * n1
f = f0 + lam * f1
g = metric(N, f)
gi, Gam, Rl, Ric, Rs, W = curvature(g)
I, J = invariants(g, gi, W)
Js = eps_series(J)
L2J = sp.series(N * r ** 2 * Js[2], lam, 0, 2).removeO()
alJ = reduce_quadratic(sp.expand(L2J))
rsp = sp.series(1 / (N * f), lam, 0, 2).removeO()


def flux_series(alpha, rsp):
    chi_der = {1: -rsp}
    for k in range(2, 5):
        chi_der[k] = sp.diff(chi_der[k - 1], r)
    Pi = 0
    for k, aa in alpha.items():
        if k == 0:
            continue
        term = 2 * aa * chi_der[k]
        for _ in range(k - 1):
            term = sp.diff(term, r)
        Pi += (-1) ** (k - 1) * term
    e = sp.series(sp.expand(-Pi), lam, 0, 2).removeO()
    return sp.expand(sum(sp.limit(sp.simplify(e.coeff(lam, j)), r, 1) * lam ** j for j in range(2)))


e_J = flux_series(alJ, rsp)          # J-part flux (units: Einstein flux = -1)
print("(ii) J-part horizon flux on corrected background (Einstein = -1):", e_J)
# Wald for J on corrected background (Riemann treated as independent; Einstein part = 4 pi r_h^2 exactly)
gB = metric(N, f, pert=False)
giB, GamB, RlB, RicB, RsB, WB = curvature(gB)
dd = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: RlB[k] + dd * Xt.get(k, 0) for k in RlB}
Wp, _ = weyl_from(Rp, gB, giB)
_, Jp = invariants(gB, giB, Wp)
sJ = -8 * sp.pi * (sp.diff(Jp, dd).subs(dd, 0) / 4) * (-gB[0, 0] * gB[1, 1]) * r ** 2
sJ = sp.series(sp.limit(sp.simplify(sp.series(sJ, lam, 0, 2).removeO().coeff(lam, 0)), r, 1)
               + lam * sp.limit(sp.simplify(sp.series(sJ, lam, 0, 2).removeO().coeff(lam, 1)), r, 1), lam, 0, 2).removeO()
print("    J-part Wald entropy (units 1/16piG; Einstein = 4 pi):", sp.expand(sJ))
eta_hat = 1 + lam * (-e_J)           # eta/eta_E at fixed r_h  (flux sign: Einstein -1 <-> +1)
s_hat = 1 + lam * sJ / (4 * sp.pi)
ratio = sp.series(eta_hat / s_hat, lam, 0, 3).removeO()
print("    4 pi eta/s = ", sp.expand(ratio))
a3 = sp.Symbol('alpha3')
rat_a3 = sp.expand(ratio.subs(lam, -a3 / 6))
print("    in alpha3 (lam = -alpha3/6):", rat_a3)
t4 = sp.Symbol('t4')
print("    in t4 (t4 = -240 alpha3, exact at tree level for pure C^3):", sp.expand(rat_a3.subs(a3, -t4 / 240)),
      "   [Li-Lu-Mai massless-cubic family: 1 - 3 t4/40 - 17 t4^2/3840 - 45/2 ((C_T - a~)/C_T)^2 + ...]")
