"""
d3: signal-to-noise of the kernel-moment densities (t19 noise floor: median |Im| of the real-by-symmetry modes, propagated through the
kernel).  Per (k, nu, wave, kernel): share of sum|d W| carried by nodes with |d| > 10 sigma, and |sum d W| / sum(sigma W) (moment over the
noise integrated over the grid).  Decides whether a moment is measured (signal) or noise-dominated (numerical).
Usage: python3 d3_snr.py cache/gw_<tag>.npz
"""
import os, sys, functools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
for p in (TT, os.path.join(TT, '..'), os.path.join(TT, '..', '..', 'e2e')):
    sys.path.insert(0, p)
import numpy as np
import t19_k2 as T19
from t13_analysis_l2 import Modes
T19.design = functools.lru_cache(None)(T19.design)
fn = sys.argv[1]; mo = Modes(fn); mo.fn = fn
W = mo.g['W']; mask = mo.conv & mo.g['ok'] & (T19.crit(mo) > 40)
print(f"# d3 {os.path.basename(fn)}: J={mo.J}; mask cut 40")
for k in (2, 4, 5, 6, 7):
    for nu in (0.0, 3.0):
        cells = []
        for wave in (('same', 'same'), ('opp', 'opp')):
            for kn in ('K*', 'A'):
                d, sig, H, Fc = T19.setup(mo, k, nu, wave, T19.KS[kn])
                a = np.abs(d*W)*mask; s = sig*W*mask
                share = a[np.abs(d) > 10*sig].sum()/a.sum()
                snr = abs((d*W*mask).sum())/s.sum()
                cells.append(f"{wave[0][0]}{wave[1][0]}-{kn}: share>10sig {share:.3f}, |M|/sum(sig) {snr:.1e}")
        print(f"  k={k} nu={nu}: " + " | ".join(cells), flush=True)
