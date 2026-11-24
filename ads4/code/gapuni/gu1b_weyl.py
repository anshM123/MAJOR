"""
gu1b_weyl (GAP_UNIFORM.md §3.1): 2x2 (even J >= 4) chart-A tolerance with the Weyl bound
   lambda_min >= min_h rig_h - sup|Gamma| (max_h W_h + W_off)   ->   G_W = min_h rig_h / (max_h W_h + W_off)   (certified-floor form)
and the hybrid form min_h lam_h / (max_h W_h + W_off), at Delta_g = 1000, 2000, 5000 (ratio to (Lam/Lam0)^2).  Read-only large_m_map.
Usage: python3 gu1b_weyl.py -> cache/gu1b_weyl.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gu_guard  # noqa: F401
sys.path.insert(0, os.path.join(HERE, '..', 'poserr')); sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import large_m_map as LM
E = LM.E
OUT = []
L0 = E.lam_of(1000.0)
OUT.append(f"{'J':>3} {'mu':>8} | {'GW_cert1000':>11} {'r2000':>7} {'r5000':>7} | {'GW_hyb1000':>10} {'r2000':>7} {'r5000':>7}")
for J in [4, 6, 8, 10, 12, 14, 16]:
    mg = E.mgap(1000.0, J)
    for mu in sorted(set([round(mg*1.000001, 6)] + [v for v in (1.02, 1.05, 1.1, 1.3) if v > mg])):
        row = []
        for Dg in (1000.0, 2000.0, 5000.0):
            s = LM.state(Dg, J, mu)
            Wt = max(s[h]['W_lo'] + s[h]['W_hi'] for h in ('++', '+-')) + s['W_off_lo'] + s['W_off_hi']
            row.append((min(s[h]['rig'] for h in ('++', '+-'))/Wt, min(s[h]['lam'] for h in ('++', '+-'))/Wt, (E.lam_of(Dg)/L0)**2))
        (c0, h0, _), (c2, h2, r2), (c5, h5, r5) = row
        f = lambda a, b, r: (a/(b*r) if b > 0 else float('nan'))
        OUT.append(f"{J:>3} {mu:8.5f} | {c0:11.3e} {f(c2, c0, r2):7.4f} {f(c5, c0, r5):7.4f} | {h0:10.3e} {f(h2, h0, r2):7.4f} {f(h5, h0, r5):7.4f}")
        print(OUT[-1], flush=True)
open(os.path.join(HERE, 'cache', 'gu1b_weyl.txt'), 'w').write("\n".join(OUT) + "\n")
