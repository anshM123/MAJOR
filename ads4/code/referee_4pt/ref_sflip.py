"""Referee: is the sign of S at the right pair fixed by the data beyond a relabelling of oo<->ss?  Compare fits with S_R = +S, -S:
uniform-weight residuals, per-class ratios, and the largest per-mode misfit."""
import sys, os, glob, numpy as np
sys.argv = [sys.argv[0], sys.argv[1]]
import ref_analysis as ra
from tmoments4 import TensorMoments4
import t5_analysis_l1 as t5
f = sorted(glob.glob(os.path.join(sys.argv[1], 'modes4_l1_D105*.npz')))[0]
mom = TensorMoments4.load(f)
D, J = mom.Delta, mom.J
for nu in (0.0, 1.0, 3.0):
  for name in ('v2v2', 'v1v1', 'ss', 'oo'):
    cL, cR = t5.combos(D, J)[name]
    for sR in (1.0, -1.0):
        A = ra.design(D, J, 6, nu, 1.0, sR)
        Tp = ra.Tpred(cL, cR, D, J)
        pred = {q: sum(A[q][c]*Tp[c] for c in ra.CLASSES) for q in ra.QS}
        big = max(abs(v) for v in pred.values()); logref = np.log(big)
        meas = t5.measured_combo(mom, cL, cR, 6, nu, logref)
        cl = list(ra.CLASSES)
        M = np.array([[A[q][c]*np.exp(-logref) for c in cl] for q in ra.QS]); y = np.array([meas[q] for q in ra.QS])
        Mr = np.vstack([M.real, M.imag]); yr = np.concatenate([y.real, y.imag])
        sol, *_ = np.linalg.lstsq(Mr, yr, rcond=None)
        res = np.linalg.norm(Mr@sol - yr)/np.linalg.norm(yr)
        # residual if only the leading predicted classes are free with their predicted ratio fixed (one overall scale)
        p = Mr@np.array([Tp[c] for c in cl]); s = (p@yr)/(p@p); res1 = np.linalg.norm(s*p - yr)/np.linalg.norm(yr)
        print(f"nu={nu} {name:4s} S_R={sR:+.0f}: uniform resid (8 free) {res:.2e}; one-scale resid {res1:.2e} scale {s:.4f}; "
              + " ".join(f"{t5.CNAME[c]}={sol[i]/Tp[c]:.3g}" if abs(Tp[c])>1e-9*max(abs(v) for v in Tp.values()) else f"{t5.CNAME[c]}=[{sol[i]:.2e}]" for i,c in enumerate(cl)))
