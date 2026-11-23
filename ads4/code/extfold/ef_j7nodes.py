"""
ef_j7nodes (EXT_FOLD.md §3.4): J = 7 direct own-row nodes (Delta = 1007, jackknife sigma) against the published chain's model of the
same row (q3 odd_model: K2 rows 800 and 1600, natural spline in x, 1/m interpolation) at the same x -- the H_M check for J = 7.
Usage: python3 ef_j7nodes.py -> cache/ef_j7nodes.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C
import numpy as np
H = C.H
J, D = 7, 1007.0
gK = H.grid(J, 2)
gI = H.grid(J, 2, inject=[D])[D]
reps = [H.grid(J, 2, inject=[D], jk=j)[D] for j in range(8)]
f, Ds = H.Q.odd_model(gK)
m = gI['m']
mod = f(gI['xs'], m)
lines = [f"# J=7 row 1007 (m={m:.3f}): Gamma (c2; no kappa subtraction, chain convention) direct +- jk vs published model (rows {Ds}) at the same x"]
dev = []
for i, x in enumerate(gI['xs']):
    s = f"   x={x:.4f}"
    for kd in ('C2', 'C3H'):
        vj = np.array([r[kd][i] for r in reps]); e = float(np.sqrt(7/8*((vj - vj.mean())**2).sum()))
        d = gI[kd][i] - mod[kd][i]; dev.append((abs(d), abs(d)/max(e, 1e-9)))
        s += f" | {kd}: dir {gI[kd][i]:+8.1f} +- {e:5.1f}  model {mod[kd][i]:+8.1f}  diff {d:+7.1f}"
    lines.append(s)
dev = np.array(dev)
lines.append(f"#  max |dir| = {max(np.abs(gI['C2']).max(), np.abs(gI['C3H']).max()):.0f} c2; |dir - model|: median {np.median(dev[:, 0]):.1f}, max {dev[:, 0].max():.1f} c2; "
             f"median {np.median(dev[:, 1]):.1f} sigma, max {dev[:, 1].max():.1f} sigma")
print("\n".join(lines))
open(os.path.join(C.EFC, 'ef_j7nodes.txt'), 'w').write("\n".join(lines) + "\n")
