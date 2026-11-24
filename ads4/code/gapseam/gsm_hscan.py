"""
gsm_hscan (GAP_SEAM.md): hybrid (flat + leading-Regge) near-gap value lambda_hyb/N_s of spin J on a fine m/Lam grid from the twist edge
to 1.15 Lam at Delta_g = 1000 and 2000 (ee1000_common.hyb, read-only), to locate each spin's near-gap minimum (the states evaluated
by gsm_bind are the twist edges, the 1000 twist edge at 2000, and the hybrid minima).  Even J >= 4: 2x2 (++, +-) minimum and the ++
diagonal (= OO block at LO).  Usage: python3 gsm_hscan.py J,J,...  -> cache/gsm_hscan.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gsm_guard  # noqa: F401
sys.path.insert(0, os.path.join(HERE, '..', 'neargap2'))
import numpy as np
import ng_chain as NC
E = NC.E
lines = []
for J in [int(v) for v in sys.argv[1].split(',')]:
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    for Dg in (1000, 2000):
        tw = E.mgap(Dg, J)
        grid = np.unique(np.r_[tw, np.arange(np.ceil(tw*400)/400, 1.1500001, 0.0025)])
        vals, oo = [], []
        for mfl in grid:
            st = {h: E.hyb(Dg, J, mfl, h, 256) for h in hs}
            base = np.array([st[h]['base'] for h in hs]); Ns = np.array([st[h]['Ns'] for h in hs])
            vals.append(E.lam2(base, Ns, np.zeros((len(hs), len(hs))))); oo.append(base[0]/Ns[0])
        vals = np.array(vals); oo = np.array(oo); i = int(np.argmin(vals)); j = int(np.argmin(oo))
        s = (f"J={J:2d} Dg={Dg}: tw {tw:.5f} -> {vals[0]:+.4f}; min {vals[i]:+.4f} at {grid[i]:.4f}" +
             (f"; ++ diag (OO LO) tw {oo[0]:+.4f}, min {oo[j]:+.4f} at {grid[j]:.4f}" if len(hs) == 2 else ""))
        print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 'gsm_hscan.txt'), 'w').write("\n".join(lines) + "\n")
