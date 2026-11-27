"""Item 5 support: Einstein-part shear kinetic coefficients on an ARBITRARY static background N(r), f(r).
Expect alpha_1 = -N r^2 f / 2 and no higher alpha  =>  horizon flux -2 alpha_1 r_*' = r^2 exactly (background independent).  [computed]"""
import sympy as sp
from static_engine import *
N, f = sp.Function('N')(r), sp.Function('f')(r)
g = metric(N, f)
gi, Gam, Rl, Ric, Rs, W = curvature(g)
e = jets(Rs)
R2 = sp.diff(e, eps, 2).subs(eps, 0) / 2
al = reduce_quadratic(sp.expand(N * r ** 2 * R2))
print("alphas:", {k: sp.simplify(v) for k, v in al.items()})
print("alpha_1 + N r^2 f/2 =", sp.simplify(al.get(1, 0) + N * r ** 2 * f / 2), "; other alphas:", {k: v for k, v in al.items() if k != 1})
