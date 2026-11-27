"""
Item 3: Wald entropy of U = (grad C)^2 and I*U from the Iyer-Wald formula, without field redefinitions.  [computed]
For L = G(I) U (G = 1 or I): E^{abcd} = dL/dR_abcd - nabla_e dL/d(nabla_e R_abcd)
                                     = G'(I) U * 2 C^{abcd} - 2 nabla_e( G nabla^e C^{abcd} ).
At the bifurcation surface nabla I = 0 (smooth Kruskal argument) and U = 180 f/r^8 = 0, so E = -2 G(I) box C^{abcd}.
Then s1/s0 = 4 * [ -G(I) (box C)^{trtr} / C^{trtr} ]_{r=1}  (I alone has E = 2 C^{abcd}, s1/s0 = 4).
"""
import itertools
import sympy as sp
from static_engine import *

f0 = r ** 2 - 1 / r
g = metric(1, f0, pert=False)
gi, Gam, Rl, Ric, Rs, W = curvature(g)
n = 4
idx = list(itertools.product(range(n), repeat=4))


def cov1(T):
    """nabla_e T_abcd -> dict (e,a,b,c,d)."""
    out = {}
    for e in range(n):
        for (a, b, c, d) in idx:
            v = sp.diff(T[a, b, c, d], X[e])
            for h in range(n):
                v -= (Gam[h][e][a] * T[h, b, c, d] + Gam[h][e][b] * T[a, h, c, d]
                      + Gam[h][e][c] * T[a, b, h, d] + Gam[h][e][d] * T[a, b, c, h])
            out[e, a, b, c, d] = sp.cancel(v)
    return out


DC = cov1(W)
boxC = {}
for (a, b, c, d) in idx:
    s = 0
    for e in range(n):
        # nabla_e (nabla_e C_abcd) with the extra index e also connected
        v = sp.diff(DC[e, e, a, b, c, d] if False else DC[e, a, b, c, d], X[e])
        for h in range(n):
            v -= Gam[h][e][e] * DC[h, a, b, c, d]
            v -= (Gam[h][e][a] * DC[e, h, b, c, d] + Gam[h][e][b] * DC[e, a, h, c, d]
                  + Gam[h][e][c] * DC[e, a, b, h, d] + Gam[h][e][d] * DC[e, a, b, c, h])
        s += gi[e, e] * v
    boxC[a, b, c, d] = sp.cancel(s)
I = sum(W[a, b, c, d] ** 2 * gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d] for (a, b, c, d) in idx)
CboxC = sp.simplify(sum(W[a, b, c, d] * boxC[a, b, c, d] * gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d] for (a, b, c, d) in idx))
ratio = sp.simplify(boxC[0, 1, 0, 1] / W[0, 1, 0, 1])
print("C . box C =", CboxC, ";   -6 I - 3 J =", sp.simplify(-6 * 12 / r ** 6 - 3 * 12 / r ** 9))
print("(box C)_trtr / C_trtr =", ratio, " at r=1:", ratio.subs(r, 1))
I0 = 12 / r ** 6
wU = sp.simplify(4 * (-ratio)).subs(r, 1)
wIU = sp.simplify(4 * (-I0 * ratio)).subs(r, 1)
print("Wald s1/s0:  U ->", wU, "   I U ->", wIU)
print("=> dV_E[U] = e1 - w = -288 -", wU, "=", -288 - wU, "   (identity route: 3 * (-108) = -324)")
print("=> dV_E[I U] = -15768 -", wIU, "=", -15768 - wIU, "   (collaborator: -16200)")
