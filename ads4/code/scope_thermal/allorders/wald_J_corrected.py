"""Item 5 support: Wald entropy of the J term on the O(lam)-corrected brane, to O(lam) (lean evaluation).  [computed]"""
import itertools
import sympy as sp
from static_engine import *
lam, dl = sp.symbols('lam delta')
f0 = r ** 2 - 1 / r
N = 1 + lam * (-6 / r ** 6)
f = f0 + lam * (-8 / r + 24 / r ** 4 - 16 / r ** 7)
gB = metric(N, f, pert=False)
gi, Gam, Rl, Ric, Rs, W = curvature(gB)
dd = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: Rl[k] + dd * Xt.get(k, 0) for k in Rl}
Wp, _ = weyl_from(Rp, gB, gi)
# derivative of J w.r.t. dd at dd = 0: J = W_ab^cd W_cd^ef W_ef^ab -> 3 * dW . W . W
Wm = {k: Wp[k] * gi[k[2], k[2]] * gi[k[3], k[3]] for k in Wp}
dWm = {k: sp.diff(v, dd).subs(dd, 0) for k, v in Wm.items()}
W0m = {k: v.subs(dd, 0) for k, v in Wm.items()}
dJ = 0
for (a, b, c, d), x1 in dWm.items():
    if x1 == 0:
        continue
    for e, h in itertools.product(range(4), repeat=2):
        x2 = W0m[c, d, e, h]
        if x2 == 0:
            continue
        x3 = W0m[e, h, a, b]
        if x3 == 0:
            continue
        dJ += 3 * x1 * x2 * x3
sJ = -8 * sp.pi * (dJ / 4) * (-gB[0, 0] * gB[1, 1]) * r ** 2
out = []
for j in range(2):
    c = sp.diff(sJ, lam, j).subs(lam, 0) / sp.factorial(j)
    c = sp.cancel(sp.together(c))
    out.append(sp.limit(c, r, 1))
print("J-part Wald entropy at r_h = 1 (units 1/16piG; Einstein = 4 pi):", out[0], "+ lam *", out[1])
print("  => w(lam) = s_J/(4 pi) =", sp.simplify(out[0] / (4 * sp.pi)), "+ lam *", sp.simplify(out[1] / (4 * sp.pi)))
