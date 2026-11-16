"""
Quadrature-grid check of the completed k = 2, 3 moments: the same (Delta, J) evaluated on the e2e grid (72 x 62 Gauss nodes) and on a
refined grid (more panels near u = 0, ~1.4x the nodes elsewhere), both completed in the corner with corner.TensorMomentsLC.
Usage: python3 q1_grid.py      Output: out_q1_grid.txt
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from moments import TensorMoments
from corner import TensorMomentsLC, CACHE
from t4_spinning_e2e import moments_of, default_N
from t7_corrections import measured_components
from cases import info
mp.mp.dps = 25
UP = ((0.0, 0.002, 8), (0.002, 0.01, 10), (0.01, 0.03, 12), (0.03, 0.1, 18), (0.1, 0.3, 22), (0.3, 1.0, 32))
TP = ((0.0, 6.0, 22), (6.0, 28.0, 34), (28.0, 85.0, 26))
if __name__ == "__main__":
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("# R_pres / R_flip on the e2e grid and on the refined grid (both corner-completed); difference = quadrature error estimate")
    for l in (1, 2):
        for (Delta, J) in ((24.221511, 10), (44.565202, 8)):
            from t4_spinning_e2e import main_cases
            Delta = [c[0] for c in main_cases() if c[1] == J and abs(c[0] - Delta) < 1e-3][0]
            m, b_, eta = info(Delta, J)
            A = TensorMomentsLC.build(l, Delta, J)
            fn = os.path.join(CACHE, f"lcmodes_fine_l{l}_D{Delta:.6f}_J{J}.npz")
            if os.path.exists(fn):
                B = TensorMomentsLC.load(fn)
            else:
                t0 = time.time()
                base = TensorMoments(l, Delta, J, 1.75, default_N(l, Delta), upanels=UP, tpanels=TP)
                pr(f"#   refined radial grid evaluated in {time.time()-t0:.0f} s")
                B = TensorMomentsLC.build(l, Delta, J, base=base, fn=fn)
            ra = moments_of(A, (2, 3, 4, 6), (0.0, 1.0, 3.0)); rb = moments_of(B, (2, 3, 4, 6), (0.0, 1.0, 3.0))
            for k in (2, 3, 4, 6):
                for nu in (0.0, 1.0, 3.0):
                    pa, fa, _, _ = measured_components(l, Delta, J, 1.75, k, nu, ra)
                    pb, fb, _, _ = measured_components(l, Delta, J, 1.75, k, nu, rb)
                    pr(f"  l={l} (Delta,J)=({Delta:.3f},{J}) m={float(m):.1f} beta={float(b_):.2f} k={k} nu={nu}:  pres {float(mp.re(pa)):.9f} / {float(mp.re(pb)):.9f}  "
                       f"diff {float(abs(pa-pb)):.1e}   flip {float(mp.re(fa)):.9f} / {float(mp.re(fb)):.9f}  diff {float(abs(fa-fb)):.1e}")
    open(os.path.join(HERE, 'out_q1_grid.txt'), 'w').write("\n".join(out) + "\n")
