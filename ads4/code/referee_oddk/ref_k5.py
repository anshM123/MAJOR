"""Referee check of SPINNING_ODDK: reproduce k=5 ratios from the cached <JJJJ> modes (independent driver, same building blocks),
and run mutations of the odd-k prediction:  b2 (b0->b2, O(1/m)),  'interp' (one-pair tables M at k=5 replaced by the mean of the
k=4 and k=6 tables: a plausible but wrong 'continuation'),  'shiftS' (S -> -S only for odd k).  Usage: python3 ref_k5.py BETA"""
import os, sys, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); FP = os.path.join(HERE, '..', 'fourpt')
sys.path.insert(0, FP); sys.path.insert(0, os.path.join(FP, 'oddk')); sys.path.insert(0, os.path.join(FP, '..', 'e2e'))
import glob
from cases import info
from tmoments4 import TensorMoments4, DATADIR
import t5_analysis_l1 as T5, predict4 as P4
from corner_est import moment_and_corner
BETA = float(sys.argv[1]); K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
orig_Mval = P4.Mval
def Mval_interp(l, q1, q3, h1, h3, nu, a):
    # tables (R, S, g1) at a-1 and a+1 averaged, gamma prefactor untouched
    R0, S0, mu = P4.Mfun(l)[(q1, q3, h1, h3)]
    def MM(aa):
        return R0(nu, aa) + P4.gamma_a(aa-1, nu)/P4.gamma_a(aa, nu)*S0(nu, aa)
    return (MM(a-1) + MM(a+1))/2, mu
def Mval_flipS(l, q1, q3, h1, h3, nu, a):
    R0, S0, mu = P4.Mfun(l)[(q1, q3, h1, h3)]
    return R0(nu, a) - P4.gamma_a(a-1, nu)/P4.gamma_a(a, nu)*S0(nu, a), mu
rows = []
for f in sorted(glob.glob(os.path.join(DATADIR, 'modes4_l1_*.npz'))):
    mom = TensorMoments4.load(f)
    m, beta, _ = info(mom.Delta, mom.J)
    if abs(float(beta) - BETA) > 0.15 or float(m) < 30: continue
    D, J = mom.Delta, mom.J
    cmb = T5.combos(D, J)
    out = {}
    for mut in (None, 'b2', 'interp', 'flipS'):
        T5.Mval = P4.Mval = {'interp': Mval_interp, 'flipS': Mval_flipS}.get(mut, orig_Mval)
        P4._Kc.clear() if False else None
        A = T5.design(D, J, K, 0.0)
        for name in ('ss', 'oo', 'v2v2'):
            cL, cR = cmb[name]
            Tp = T5.Tpred(cL, cR, D, J, mut=('b2' if mut == 'b2' else None))
            big = max(abs(sum(A[q][c]*Tp[c] for c in T5.CLASSES)) for q in T5.QS); logref = np.log(big)
            w = {(a, b): cL[a]*cR[b] for a in range(2) for b in range(2) if cL[a]*cR[b] != 0}
            raw, cor = {}, {}
            for q in T5.QS:
                v, e, bd, sh, sk = moment_and_corner(mom, K, 0.0, None, q, w)
                fac = np.exp(sh - logref)*(-1)**J
                raw[q] = complex(v)*fac; cor[q] = complex(v + e)*fac
            for lab, meas in (('raw', raw), ('cor', cor)):
                R, sol, res = T5.fit_components(A, meas, Tp, logref)
                cl = (1, -1, -1, 1) if name in ('ss', 'v2v2') else (1, 1, -1, -1)
                out[(mut, name, lab)] = R[cl]
    T5.Mval = P4.Mval = orig_Mval
    rows.append((float(m), out))
    print(f"m={float(m):.2f} beta={float(beta):.3f}", {k_: round(v, 5) for k_, v in out.items()}, flush=True)
rows.sort(key=lambda r: r[0])
for key in rows[0][1]:
    ms = [r[0] for r in rows]; Rs = [r[1][key] for r in rows]
    f = T5.fits(ms, Rs)
    print(key, " R:", " ".join(f"{x:.4f}" for x in Rs), "| c1(A) = %+.3f" % f['A'][0] if 'A' in f else "",
          "| m^2(R-1) last:", " ".join(f"{m*m*(R-1):+.1f}" for m, R in zip(ms[-2:], Rs[-2:])))
