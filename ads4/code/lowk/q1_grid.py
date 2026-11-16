"""
Quadrature-grid check of the completed moments: the same (Delta, J) on the e2e grid (72 x 62 Gauss nodes; u-panels [0,0.01] (8),
[0.01,0.03] (10), ...) and on the corner-refined grid of corner.refined_base (46 nodes on six panels below u = 0.03, from u = 3e-6;
identical nodes above), both completed with corner.TensorMomentsLC.  The difference estimates the u-quadrature error of the corner.
Usage: python3 q1_grid.py      Output: out_q1_grid.txt
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from corner import TensorMomentsLC, CACHE, refined_base, UP_REF2, moments_with_err
from t4_spinning_e2e import moments_of, main_cases
from t7_corrections import measured_components
from cases import info
mp.mp.dps = 25
CASES = [(1, 8), (1, 20), (1, 52), (1, 82), (2, 8), (2, 20)]      # (l, J): beta ~ 0.4 (m ~ 42, 102), beta ~ 1 (m ~ 101), beta ~ 2 (m ~ 70)
DOUBLE = [(1, 20), (1, 82)]                                         # also a grid C with ~1.6x the corner nodes (corner.UP_REF2)
if __name__ == "__main__":
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("# R_pres / R_flip on the e2e grid (A) and on the corner-refined grid (B; used for the results); |A-B|: quadrature error of A;")
    pr("# |B-C|: quadrature error of B (C: ~1.6x the corner nodes, two cases); [B: ...]: the error budget of B (relative errors of the two modes)")
    for (l, J) in CASES:
        Delta = [c[0] for c in main_cases() if c[1] == J][0]
        m, b_, eta = info(Delta, J)
        A = TensorMomentsLC.build(l, Delta, J)
        fn = os.path.join(CACHE, f"lcmodes_ref_l{l}_D{Delta:.6f}_J{J}.npz")
        B = TensorMomentsLC.build(l, Delta, J, base=refined_base(l, Delta, J), fn=fn)
        ra = moments_with_err(A, (2, 3, 4, 6), (0.0, 1.0, 3.0))[0]; rb = moments_with_err(B, (2, 3, 4, 6), (0.0, 1.0, 3.0))[0]
        C = None
        if (l, J) in DOUBLE:
            fn2 = os.path.join(CACHE, f"lcmodes_ref2_l{l}_D{Delta:.6f}_J{J}.npz")
            C = TensorMomentsLC.build(l, Delta, J, base=refined_base(l, Delta, J, UP_REF2), fn=fn2)
            rc = moments_with_err(C, (2, 3, 4, 6), (0.0, 1.0, 3.0))[0]
        for k in (2, 3, 4, 6):
            for nu in (0.0, 1.0, 3.0):
                pa, fa, _, _ = measured_components(l, Delta, J, 1.75, k, nu, ra)
                pb, fb, _, _ = measured_components(l, Delta, J, 1.75, k, nu, rb)
                e1 = B.moment_err(k, nu, l, -l); e2 = B.moment_err(k, nu, l, l)
                line = (f"  l={l} (Delta,J)=({Delta:.3f},{J}) m={float(m):.1f} beta={float(b_):.2f} k={k} nu={nu}:  pres {float(mp.re(pa)):.9f} / {float(mp.re(pb)):.9f} "
                        f"|A-B| {float(abs(pa-pb)):.1e}   flip {float(mp.re(fa)):.9f} / {float(mp.re(fb)):.9f} |A-B| {float(abs(fa-fb)):.1e}"
                        f"   [B: mode errs lc {e1[1]:.0e}/{e2[1]:.0e} drop {e1[2]:.0e}/{e2[2]:.0e}]")
                if C is not None:
                    pc, fc, _, _ = measured_components(l, Delta, J, 1.75, k, nu, rc)
                    line += f"   C: |B-C| {float(abs(pb-pc)):.1e} / {float(abs(fb-fc)):.1e}"
                pr(line)
    open(os.path.join(HERE, 'out_q1_grid.txt'), 'w').write("\n".join(out) + "\n")
