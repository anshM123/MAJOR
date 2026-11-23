"""
t1a_einstein2 (DSKETCH_UP.md §1.2): the second-order (a*b) Einstein equation on Euclidean Poincare AdS4 with two ASD legs.

g = (delta + a gam3 + b gam4 + a b Gam)/z^2,  gam3, gam4 regular helicity-(-) TT modes (ASD, c1val/v2), Gam radial gauge
(Gam_{z mu} = 0) with unknown components F_ij(z) (times exp(i K.x), K = k3 + k4).  E_{mu nu} = R_{mu nu} + 3 g_{mu nu}.
Checks the O(a), O(b) equations vanish, and writes the O(ab) equations  L[F] + S = 0  to cache/t1a_einstein2.pkl
(srepr strings) for t1a_solve2.py.
Usage: python3 -B t1a_einstein2.py k3x k3y k3z k4x k4y k4z
"""
import os, sys, pickle, time, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from asd_alg import Alg, z, N, christoffel, riemann_up, ricci, inverse_conformal, mode

t0 = time.time()
k3 = [sp.Integer(int(v)) for v in sys.argv[1:4]] if len(sys.argv) > 6 else [0, 0, 3]
k4 = [sp.Integer(int(v)) for v in sys.argv[4:7]] if len(sys.argv) > 6 else [4, 0, 0]
A = Alg(k3, k4)
SUM = sp.sqrt(sum(v**2 for v in k3)) + sp.sqrt(sum(v**2 for v in k4))
g3 = mode(k3, -1); g4 = mode(k4, -1)
F = {}
for i in range(3):
    for j in range(i, 3):
        F[(i, j)] = sp.Function(f'F{i}{j}')(z)
H = [[None]*N for _ in range(N)]
for i, j in itertools.product(range(N), repeat=2):
    t = {(1, 0): g3[i, j], (0, 1): g4[i, j]}
    if i < 3 and j < 3:
        t[(1, 1)] = F[(min(i, j), max(i, j))]
    H[i][j] = t
g = [[A.scale(1/z**2, A.add({(0, 0): sp.Integer(1 if i == j else 0)}, H[i][j])) for j in range(N)] for i in range(N)]
gic = inverse_conformal(A, H)
gi = [[A.scale(z**2, gic[i][j]) for j in range(N)] for i in range(N)]
G = christoffel(A, g, gi)
print('christoffel', round(time.time() - t0, 1), flush=True)
Rup = riemann_up(A, G)
Ric = ricci(A, Rup)
print('ricci', round(time.time() - t0, 1), flush=True)
ok = True
eqs = {}
for m, n in itertools.product(range(N), repeat=2):
    if n < m:
        continue
    E = A.add(Ric[m][n], A.scale(3, g[m][n]))
    for o in [(0, 0), (1, 0), (0, 1)]:
        v = sp.simplify(E.get(o, 0))
        if v != 0:
            print('NONZERO', (m, n), o, v); ok = False
    eqs[(m, n)] = sp.expand(E.get((1, 1), 0)*sp.exp(SUM*z)*z**2)   # strip the e^{-(|k3|+|k4|) z} of the source
print('background + linear orders vanish:', ok, ' time', round(time.time() - t0, 1), flush=True)
out = os.path.join(HERE, 'cache', 't1a_einstein2_%s.pkl' % '_'.join(str(v) for v in list(k3) + list(k4)))
pickle.dump({'k3': [str(v) for v in k3], 'k4': [str(v) for v in k4], 'ok': ok,
             'eqs': {key: sp.srepr(val) for key, val in eqs.items()}}, open(out, 'wb'))
print('wrote', out)
