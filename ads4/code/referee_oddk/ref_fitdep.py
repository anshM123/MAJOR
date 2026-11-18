"""Referee: dependence of R_ss (W[s,s]) and R_oo (W[o,o]) on the extraction: 8-class fit with 1/|pred| weights (as in the paper),
uniform-weight 8-class fit, and the single leading mode; raw vs corner-corrected.  nu = 0.  Usage: python3 ref_fitdep.py K file..."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); FP = os.path.join(HERE, '..', 'fourpt')
sys.path.insert(0, FP); sys.path.insert(0, os.path.join(FP, 'oddk')); sys.path.insert(0, os.path.join(FP, '..', 'e2e'))
from tmoments4 import TensorMoments4
from cases import info
import t5_analysis_l1 as T5
from corner_est import moment_and_corner
K = int(sys.argv[1])
def ufit(A, meas, Tp, logref):
    cl = list(T5.CLASSES)
    M = np.array([[A[q][c]*np.exp(-logref) for c in cl] for q in T5.QS]); y = np.array([meas[q] for q in T5.QS])
    Mr = np.vstack([M.real, M.imag]); yr = np.concatenate([y.real, y.imag])
    sol = np.linalg.lstsq(Mr, yr, rcond=None)[0]
    return {c: sol[i]/Tp[c] if abs(Tp[c]) > 1e-6*max(abs(v) for v in Tp.values()) else None for i, c in enumerate(cl)}
for fn in sys.argv[2:]:
    fc = os.path.join(FP, 'oddk', 'cache', os.path.basename(fn).replace('modes4_', 'modes4c_'))
    mom = TensorMoments4.load(fc if os.path.exists(fc) else fn)
    m, beta, _ = info(mom.Delta, mom.J)
    A = T5.design(mom.Delta, mom.J, K, 0.0)
    line = f"k={K} beta={float(beta):.2f} m={float(m):6.2f}:"
    for name, cl in (('ss', (1, -1, -1, 1)), ('oo', (1, 1, -1, -1))):
        cL, cR = T5.combos(mom.Delta, mom.J)[name]
        Tp = T5.Tpred(cL, cR, mom.Delta, mom.J)
        pred = {q: sum(A[q][c]*Tp[c] for c in T5.CLASSES) for q in T5.QS}
        logref = np.log(max(abs(v) for v in pred.values()))
        w = {(a, b): cL[a]*cR[b] for a in range(2) for b in range(2) if cL[a]*cR[b] != 0}
        raw, cor = {}, {}
        for q in T5.QS:
            v, e, bd, sh, sk = moment_and_corner(mom, K, 0.0, None, q, w)
            f = np.exp(sh - logref)*(-1)**mom.J
            raw[q] = complex(v)*f; cor[q] = complex(v + e)*f
        qlead = max(T5.QS, key=lambda q: abs(pred[q]))
        for lab, meas in (('raw', raw), ('cor', cor)):
            Rw = T5.fit_components(A, meas, Tp, logref)[0][cl]; Ru = ufit(A, meas, Tp, logref)[cl]
            Rl = (meas[qlead]/(pred[qlead]*np.exp(-logref))).real
            line += f"  {name}-{lab}: wfit {Rw:.4f} ufit {Ru:.4f} lead {Rl:.4f}"
    print(line, flush=True)
