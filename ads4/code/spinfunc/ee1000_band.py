"""
ee1000_band (EE_DG1000_LARGEJ.md §5): the uncorrected per-state value lambda0/N_s (= the leading-Regge band model M' of EE_FUNCTIONAL §4
for Omega^EE: flat NF02 + leading-Regge curvature, off-diagonal 0) at Delta_g = 1000 on the discrete states (J, m/Lam), all J up to
b^ = (2J+1)/(m/Lam) <= 2 Lam (EE1000_BJMAX), m/Lam in the ee1000_chain grid (plus each J's gap mass), both helicities (odd J: +- only).
Output: cache/ee1000_band.txt, cache/ee1000_band.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_common as E
DG = 1000
Lm = E.lam_of(DG)
MG = [1.05, 1.075, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]
Js = sorted(set(list(range(4, 130)) + list(range(130, 301, 10))))      # larger b^: EE_FUNCTIONAL §4 (continuous b^ <= 2 Lam)
FAIL = []
out, lines = [], []
for J in Js:
    mg = E.mgap(DG, J)
    for mfl in [mg] + [v for v in MG if v > mg + 1e-6]:
        if (2*J + 1)/mfl > 2*Lm:
            continue
        hs = ('+-',) if J % 2 else ('++', '+-')
        try:
            st = {h: E.hyb(DG, J, mfl, h, 128 if J < 100 else 512) for h in hs}   # 128 nodes fail at J = 300 (checked)
        except Exception as ex:
            FAIL.append((J, mfl, type(ex).__name__)); continue
        lam = min(st[h]['base']/st[h]['Ns'] for h in hs)
        out.append(dict(J=J, mfl=mfl, lam0=lam, bh=(2*J + 1)/mfl, Ns=[st[h]['Ns'] for h in hs]))
pickle.dump(out, open(os.path.join(E.C, 'ee1000_band.pkl'), 'wb'))
for lo, hi in ((4, 15), (16, 31), (32, 64), (65, 128), (129, 300)):
    sel = [r for r in out if lo <= r['J'] <= hi]
    if sel:
        w = min(sel, key=lambda r: r['lam0'])
        lines.append(f"J in [{lo},{hi}]: {len(sel)} states; min lam0/N_s {w['lam0']:+.4f} at J={w['J']}, m/Lam={w['mfl']:.4f}, b^={w['bh']:.0f}; "
                     f"near-gap (m <= 1.3 Lam) min {min(r['lam0'] for r in sel if r['mfl'] <= 1.3):+.4f}")
lines.append(f'failed evaluations (mpmath non-convergence): {FAIL}')
print("\n".join(lines))
open(os.path.join(E.C, 'ee1000_band.txt'), 'w').write("\n".join(lines) + "\n")
