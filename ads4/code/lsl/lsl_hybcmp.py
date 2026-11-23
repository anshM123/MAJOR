"""
lsl/lsl_hybcmp.py (LARGESPIN_L.md §3.3): the project's hybrid floor (flat NF02 partial wave + leading-Regge curvature,
largespin/ls_core.hyb_fast, read-only import) against the pure leading-Regge action F_LR = 2 Omega_2 + u sigma Omega_3 (lsl_core),
near-gap and 1.3 m_gap states, J = 100 ... 3000; and the spin J*(m) at which beta(J, m) >= beta_L (Lemma LR-inf threshold).
Output: cache/hybcmp.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
import numpy as np
import lsl_core as C
import ls_core as LS

L = C.lam(1000)
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

pr(__doc__.split('Output')[0])
for mf in (1.0, 1.3):
    for J in (100, 344, 1000, 2000, 3000):
        m = C.mgap(1000, J)*mf; mfl = m/L; u = 1/mfl**2
        for h in ('++', '+-'):
            hb = LS.hyb_fast(1000, J, mfl, h)
            F = C.LR_action(hb['beta'], u, C.SIGMA[h])
            pr(f"  J = {J:5d}, m/mgap = {mf}, {h}: beta = {hb['beta']:.4f}, bhat = {hb['bhat']:.1f}: hybrid base = {hb['base']:.6e}, "
               f"F_LR = {F:.6e}, (hybrid - F_LR)/F_LR = {(hb['base'] - F)/F:+.3e}, N_flat-part B2+B3 = {hb['B2'] + hb['B3']:.3e}")
pr("\n  J*(m): smallest J with beta(J, m) >= beta_L, i.e. bhat >= 2 Lam sinh(beta_L/2)")
for bL in (1.0, 1.5, 2.0):
    bh = 2*L*np.sinh(bL/2)
    row = f"  beta_L = {bL}: bhat >= {bh:.1f} ({bh/L:.3f} Lam): "
    for mf in (1.0, 1.1, 1.3):
        J = 0
        while True:
            m = C.mgap(1000, J)*mf
            if C.beta_of(J, m) >= bL:
                break
            J += 1
        row += f" m = {mf} m_gap: J* = {J};"
    for mL in (1.5, 2.0, 5.0):
        J = int(np.ceil((bh*mL - 1)/2)); row += f" m = {mL} Lam: J* = {J};"
    pr(row)
with open(os.path.join(C.CACHE, 'hybcmp.txt'), 'w') as fh:
    fh.write("\n".join(lines) + "\n")
