"""
d2: roundoff diagnostic of the exact waves at selected grid points: the Taylor-mode word application of gwaves (float64) repeated in
long double (same float64 seed radial sums, cache/<tag>_part/seedS.npz), for W[X_s,X_s] and W[X_o,X_o].  Prints, per point, |modes|
summed over q for each wave, and the relative change float64 -> longdouble.  Small-t points vs bulk points.
Usage: python3 d2_ldcheck.py Dg mfl J
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import gwaves as GW
TM = GW.TM
from moments import make_grid, frame_data, frame_positions, series_converged
from moments4 import fourier4
from spinblock2 import invariants4
Dg, mfl, J = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
Lam, m, Dv = GW.state(Dg, mfl, J)
g = make_grid(Dv, J); ok = g['ok']
fd = frame_data(g['U'][ok], g['S'][ok])
conv = series_converged(fd['rs'], fd['chi'], J, 8000, cut=25.0)
fd = {k: v[conv] for k, v in fd.items()}
Tfull = np.broadcast_to(g['ts'][None, :], g['S'].shape)[ok][conv]
Ufull = g['U'][ok][conv]
tag = 'gw_%d_%.3f_%d' % (Dg, mfl, J)
S = np.load(os.path.join(HERE, 'cache', tag + '_part', 'seedS.npz'), allow_pickle=True)['S'].item()
rng = np.random.default_rng(1)
small = np.flatnonzero((Tfull > 0.8) & (Tfull < 2.5)); bulk = np.flatnonzero((Tfull > 8) & (Tfull < 16))
idx = np.concatenate([rng.choice(small, 12, replace=False), rng.choice(bulk, 8, replace=False)])
res = {}
for dt in (np.float64, np.longdouble):
    t0 = time.time()
    TW = TM.TaylorWaves(Dv, J, 8000, dt); TW.S = S
    x1, x2, x3, x4 = frame_positions(fd['ru'][idx], fd['xi'][idx])
    z1, z3, z2, z4 = TM.pol_grid4(fd['xi'][idx], 2)
    s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
    if dt is np.longdouble:
        s = {k: np.asarray(v, dtype=np.longdouble) for k, v in s.items()}
    Cs = TW.chunk(s, fd['rs'][idx], fd['chi'][idx], fd['A0'][idx], fd['E0'][idx], sl=idx)
    out = {}
    for c in (('same', 'same'), ('opp', 'opp')):
        Mv = TM.contract(Cs[c], z, idx.size).reshape((idx.size,) + (5,)*4)
        out[c] = fourier4(Mv, 2)
    res[dt] = out
    print(f"# {np.dtype(dt).name}: {time.time() - t0:.0f} s", flush=True)
print(f"# Dg={Dg} m/Lam={mfl} J={J}: per point t, u, |W| (sum over modes, float64) and max_q |W_ld - W_d|/max_q |W_ld|")
for j, i in enumerate(idx):
    cells = []
    for c in (('same', 'same'), ('opp', 'opp')):
        a = np.array([res[np.float64][c][q][j] for q in res[np.float64][c]]); b = np.array([res[np.longdouble][c][q][j] for q in res[np.float64][c]])
        cells.append(f"{c[0][0]}{c[1][0]}: |W| {np.abs(a).max():.3e} rel.chg {np.abs(a - b).max()/np.abs(b).max():.1e}")
    print(f"  t={Tfull[i]:6.3f} u={Ufull[i]:.4f}  " + " | ".join(cells), flush=True)
