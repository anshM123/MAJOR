"""
lsl/lsl_lrinf_grid.py (LARGESPIN_L.md §3): Lemma LR-inf hybrid-level relative remainder r(beta) = (2|E_2| + |E_3|)/(A1min M1/Lam^2)
(u <= 1, worse sigma) on a fine beta grid: step 0.01 on [0.9, 3], 0.05 on [3, 16], 0.5 on [16, 64]; reports monotonicity.
M1: numeric (lsl_core.M1) for beta <= 16; the certified lower bound of h13_cert_m1 region (c) beyond.  Output: cache/lrinf_grid.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import lsl_core as C
import lsl_lrinf as LI
L = C.lam(1000)
A1min = 2*C.d0(2)[1] - C.d0(3)[1]
grid = np.concatenate([np.arange(0.9, 3.0, 0.01), np.arange(3.0, 16.0, 0.05), np.arange(16.0, 64.01, 0.5)])
rs = []
for b in grid:
    if b <= 16:
        M1 = C.M1(b)
    else:
        D = np.cosh(b) - 1; K0 = D**-0.5
        M1 = (np.sqrt(2)/np.pi)*K0*(1/b - 1/(2*np.sinh(b/2)) - np.sqrt(2)*(np.exp(-b/2)/4.7 + 4/b**2) - 2*np.sqrt(2*np.log(2))/(b - np.log(2))**2)
    E2, _ = LI.bound(2, b); E3, _ = LI.bound(3, b)
    rs.append((2*E2 + E3)/(A1min*M1/L**2))
rs = np.array(rs)
lines = [__doc__.split('Output')[0]]
for b0 in (1.0, 1.2, 1.5, 2.0, 3.0, 5.0):
    sel = grid >= b0 - 1e-9
    lines.append(f"  sup_(beta >= {b0}) r(beta) on the grid = {rs[sel].max():.3e}   (at beta = {grid[sel][rs[sel].argmax()]:.2f})")
lines.append(f"  r decreasing along the whole grid: {bool(np.all(np.diff(rs) < 0))}; r(64) = {rs[-1]:.2e}")
print("\n".join(lines))
open(os.path.join(C.CACHE, 'lrinf_grid.txt'), 'w').write("\n".join(lines) + "\n")
