"""hs_nodes (HIGHSPIN.md): Gamma (c2; chain conventions, NO kappa subtraction for even J = the w2/NEARGAP2 convention) at the direct
nodes of a state row, K2 next to direct +- jackknife, and summary statistics of direct - K2.
Usage: python3 hs_nodes.py J Delta [xmax]   -> cache/hs_nodes_J<J>_D<Delta>.txt"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_core as H
import numpy as np
J, D = int(sys.argv[1]), float(sys.argv[2]); xmax = float(sys.argv[3]) if len(sys.argv) > 3 else 9
gK = H.grid(J, 2)[D]
gI = H.grid(J, 2, inject=[D]); gI = gI[D]
reps = [H.grid(J, 2, inject=[D], jk=j)[D] for j in range(8)]
keys = H.keys_of(J)
L = [f"# hs_nodes J={J} row {D:g}: Gamma (c2), K2 / direct +- jk; reference zero at x = {H.ref_zero_x(J, D):.4f}"]
diffs = {k: [] for k in keys}
xd = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
for i, x in enumerate(gI['xs']):
    if round(x, 6) not in xd or x > xmax:
        continue
    s = f"x={x:.4f}"
    j0 = int(np.argmin(np.abs(gK['xs'] - x))); on = abs(gK['xs'][j0] - x) < 1e-6
    for k in keys:
        vj = np.array([r[k][int(np.argmin(np.abs(r['xs'] - x)))] for r in reps]); e = np.sqrt(7/8*((vj - vj.mean())**2).sum())
        kk = k if isinstance(k, str) else k[0] + k[1]
        s += f" | {kk}: " + (f"{gK[k][j0]:+.0f}/" if on else "") + f"{gI[k][i]:+.0f}+-{e:.0f}"
        if on:
            diffs[k].append(((gI[k][i] - gK[k][j0]), e))
    L.append(s)
for k, v in diffs.items():
    if v:
        d = np.array([a for a, b in v]); e = np.array([b for a, b in v])
        L.append(f"# {k}: direct - K2 over {len(d)} nodes: median |d| {np.median(np.abs(d)):.0f}, max |d| {np.abs(d).max():.0f}, "
                 f"median |d|/sigma {np.median(np.abs(d)/np.maximum(e, 1e-9)):.1f}, max {np.max(np.abs(d)/np.maximum(e, 1e-9)):.1f}")
print("\n".join(L))
open(os.path.join(H.HSC, 'hs_nodes_J%d_D%g.txt' % (J, D)), 'w').write("\n".join(L) + "\n")
