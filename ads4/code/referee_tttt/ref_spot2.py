"""Referee spot check of the saved Taylor-mode wave caches: recompute the Fourier modes at a few grid points (largest r, i.e. near the
c -> 1 region of the arccosh bug, plus random ones) in float64, long double, and with the OLD (expanded) arccosh derivatives.
Usage: python3 ref_spot.py beta m npz"""
import os, sys, time, functools
import numpy as np, sympy as sp
T = '/home/user/MAJOR/ads4/code/fourpt/tttt'
sys.path.insert(0, T)
import tmoments_l2 as TM
from tmoments_l2 import TaylorWaves, COMBOS, contract, L
from moments import make_grid, frame_data, frame_positions, series_converged
from cases import cases_fixed_beta
from moments4 import pol_grid4, fourier4
from spinblock2 import invariants4

beta, mm, fn = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
(Dv, Jv), = cases_fixed_beta(beta, [mm])
N = int(min(8000, max(4000, 400*Dv)))
g = make_grid(Dv, Jv); ok = g['ok']
fd = frame_data(g['U'][ok], g['S'][ok])
conv = series_converged(fd['rs'], fd['chi'], Jv, N, cut=25.0)
convfull = np.zeros(g['S'].shape, dtype=bool); convfull[ok] = conv
fd = {k: v[conv] for k, v in fd.items()}
flat = np.flatnonzero(convfull.ravel())
d = np.load(fn)
cache = {}
for key, v in zip(d['keys'], d['modes']):
    cache[(COMBOS[int(key[0])],) + tuple(int(x) for x in key[1:])] = v.ravel()
rng = np.random.default_rng(1)
order = np.argsort(-fd['rs'])
sel = np.array([int(np.argmin(abs(fd["rs"] - t))) for t in (0.86, 0.88, 0.92, 0.95)])
print(f"# D={Dv:.6f} J={Jv} N={N}; selected r = {np.round(fd['rs'][sel], 5)} ; max r on grid {fd['rs'].max():.5f}")

def modes_at(dt, old=False):
    if old:
        @functools.lru_cache(None)
        def _old(M):
            c = sp.Symbol('c')
            fs = [sp.lambdify(c, sp.simplify(sp.diff(sp.acosh(c), c, k)), 'numpy') for k in range(1, M + 1)]
            return [(lambda f: (lambda c0, sh: f(c0)))(f) for f in fs]
        TM._acosh_np = _old
    TW = TaylorWaves(Dv, Jv, N, dt)
    x1, x2, x3, x4 = frame_positions(fd['ru'][sel], fd['xi'][sel])
    z1, z3, z2, z4 = pol_grid4(fd['xi'][sel], L)
    s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
    Cs = TW.chunk(s, fd['rs'][sel], fd['chi'][sel], fd['A0'][sel], fd['E0'][sel])
    out = {}
    for c in COMBOS:
        Mv = contract(Cs[c], z, len(sel)).reshape(len(sel), -1).reshape((len(sel),) + (5,)*4)
        for q, v in fourier4(Mv, L).items():
            out[(c,) + q] = v
    return out

t = time.time(); f64 = modes_at(np.float64); print(f"float64 {time.time()-t:.0f}s", flush=True)
t = time.time(); f80 = modes_at(np.longdouble); print(f"longdouble {time.time()-t:.0f}s", flush=True)
fold = f80
for c in COMBOS:
    keys = [k for k in f64 if k[0] == c]
    def dev(A, B):
        num = np.array([[abs(A[k][i] - B[k][i]) for k in keys] for i in range(len(sel))])
        den = np.array([max(abs(B[k][i]) for k in keys) for i in range(len(sel))])
        return num.max(axis=1)/den
    cached = {k: cache[k][flat[sel]] for k in keys}
    print(c, " |f64-cache|/max:", np.array2string(dev(f64, cached), precision=1),
          " |f64-ld|/max:", np.array2string(dev(f64, f80), precision=1),
          " |old-ld|/max:", np.array2string(dev(fold, f80), precision=1))
