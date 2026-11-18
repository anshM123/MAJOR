"""Referee: per-mode size of the dropped tiny-t correction (estimate e and bound b relative to the accepted moment v) for the W[s,s]
combination, and structure of the dropped set (holes above the first accepted t-node; rows skipped).  Usage: python3 ref_corner.py FILE K"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); FP = os.path.join(HERE, '..', 'fourpt')
sys.path.insert(0, FP); sys.path.insert(0, os.path.join(FP, 'oddk')); sys.path.insert(0, os.path.join(FP, '..', 'e2e'))
from tmoments4 import TensorMoments4
import t5_analysis_l1 as T5
from corner_est import moment_and_corner
mom = TensorMoments4.load(sys.argv[1]); K = int(sys.argv[2])
g = mom.g; drop = g['ok'] & ~mom.conv
holes = 0; skip = 0
for iu in range(drop.shape[0]):
    acc = np.where(mom.conv[iu])[0]; dr = np.where(drop[iu])[0]
    if len(acc) < 2: skip += len(dr); continue
    holes += int(np.sum(dr > acc[0]))
print(f"Delta={mom.Delta:.3f} J={mom.J}: ok={g['ok'].sum()} conv={mom.conv.sum()} dropped={drop.sum()} holes(above first accepted)={holes} in skipped rows={skip}")
print("ts ascending?", np.all(np.diff(g['ts'][0]) > 0) if g['ts'].ndim == 2 else np.all(np.diff(g['ts']) > 0))
cL, cR = T5.combos(mom.Delta, mom.J)['ss']
w = {(a, b): cL[a]*cR[b] for a in range(2) for b in range(2)}
A = T5.design(mom.Delta, mom.J, K, 0.0); Tp = T5.Tpred(cL, cR, mom.Delta, mom.J)
pred = {q: sum(A[q][c]*Tp[c] for c in T5.CLASSES) for q in T5.QS}; pm = max(abs(v) for v in pred.values())
lst = []
for q in T5.QS:
    v, e, b, sh, sk = moment_and_corner(mom, K, 0.0, None, q, w)
    if abs(pred[q]) > 1e-3*pm: lst.append((abs(pred[q])/pm, q, abs(e/v), b/abs(v), sk))
lst.sort(reverse=True)
for r in lst[:12]: print(f"  |pred|/max={r[0]:.2e} q={r[1]} |est/v|={r[2]:.2e} bound/|v|={r[3]:.2e} skipped={r[4]}")
