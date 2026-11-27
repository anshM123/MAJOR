"""Item 1/2: shear response for theta(r)-weighted I and J with GENERIC theta, and direct gradient operators.  [computed]"""
import sympy as sp
from static_engine import *

f0 = r ** 2 - 1 / r
gP = metric(1, f0)
gi, Gam, Rl, Ric, Rs, W = curvature(gP)
I, J = invariants(gP, gi, W)
sq = r ** 2                                      # sqrt(-g), exact (det e^{eps H} = 1)
Is, Js, Rss = eps_series(I), eps_series(J), eps_series(Rs)
print("background I, J:", Is[0], Js[0], "  first order:", Is[1], Js[1])
rsp = 1 / f0
aE = reduce_quadratic(sq * Rss[2])
eta0 = flux_at_horizon(aE, rsp)
print("Einstein alphas:", aE, " eta0 =", eta0)
th = sp.Function('theta')(r)
T = sp.symbols('t0:6')
thpoly = sum(T[k] * (r - 1) ** k / sp.factorial(k) for k in range(6))   # local jets of theta at the horizon
def loc(al):
    return {k: sp.simplify(v.subs(th, thpoly).doit()) for k, v in al.items()}
wI, wJ = 4, 6          # Wald s1/s0 for I and J (../etas_diag.py; theta-weighted: s1/s0 = w * theta(1))
for name, ser, w in (('I', Is, wI), ('J', Js, wJ)):
    al = reduce_quadratic(sq * th * ser[2])
    print(f"[theta {name}] alphas:", {k: sp.factor(v) for k, v in al.items()})
    e1 = flux_at_horizon(loc(al), rsp) / eta0
    dV = sp.simplify(e1 - w * T[0])
    print(f"[theta {name}] (t0 = theta(1), t1 = theta'(1), ...) eta1/eta0 = {sp.simplify(e1)}  ->  dV_E = {dV}")
