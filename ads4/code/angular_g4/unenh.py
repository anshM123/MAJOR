"""Unenhanced flip combinations: kappa with no T- (and C-) growth modulo the cross-ladder identity."""
import numpy as np, sys
from fractions import Fraction as Fr
import sympy as sp
N = int(sys.argv[1])
labs = [3] + list(range(4, N + 1)) + [-p for p in range(4, N + 1)]
idx = {l: i for i, l in enumerate(labs)}
nv = len(labs)
mu = lambda s: sp.Rational((s + 1) * (s - 3), (s + 4) * s)
def ab(p):
    a = [0] * nv; b = [0] * nv
    if p >= 4 and p <= N:
        a[idx[p]] += sp.Rational(1, p - 1); b[idx[p]] += -sp.Rational(1, p - 3)
        a[idx[-p]] += sp.Rational(1, p + 3); b[idx[-p]] += -sp.Rational(1, p + 1)
    if p == 4:
        a[idx[3]] += -sp.Rational(4, 21)
    return a, b
rows = []
for p in range(4, N + 2):
    a, _ = ab(p); _, bm = ab(p - 1) if p - 1 >= 4 else (None, [0] * nv)
    rows.append([a[i] + mu(p - 1) * bm[i] for i in range(nv)])
M = sp.Matrix(rows)
ns = M.nullspace()
rhs = sp.Matrix([[-2 if l == 3 else 1 for l in labs]])
print('N', N, 'nullspace dim', len(ns), ' normalisation values on basis:', [sp.nsimplify((rhs * v)[0]) for v in ns][:6])
# rhs in rowspace(M)?  solve y^T M = rhs
y = sp.symbols('y0:%d' % M.shape[0])
sol = sp.solve(list(sp.Matrix([y]) * M - rhs), y, dict=True)
print('rhs = y^T M solvable:', bool(sol), ' y =', [sol[0].get(v, v) for v in y] if sol else None)
