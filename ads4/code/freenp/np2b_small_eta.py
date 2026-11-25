"""
np2b_small_eta.py -- the eta -> 0 approach of I_eta for the slowly converging cases of np2 (m = 4 on S1, S2), down to
eta = 1/1280, and the predicted onset scale eta_* ~ (typical K)/(m^2): the regulator must resolve H-frequencies ~ m^2/K.
Output: cache/np2b_small_eta.txt
"""
import numpy as np, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from lfcore import profile
from np2_onepart import lf_bump, build_state, quantities

out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

etas = tuple(0.0125/2**k for k in range(5))
for name, cq in (('S1', lf_bump(0.6, 0.06, 0.3, 0.15)), ('S2', lf_bump(0.2, 0.025, 0.4, 0.2))):
    sl, _ = build_state(cq, nK=40)
    for m in (4, 8):
        r = quantities(sl, lambda y: profile(y, m, 're'), etas)
        P(f"{name} re f_{m}: W = {r['W']:.10f}")
        for e in etas:
            P(f"   eta = {e:.6f}:  I_eta = {r['I'][e]:.10f}   I - W = {r['I'][e] - r['W']: .3e}   (I-W)/eta^2 = {(r['I'][e] - r['W'])/e**2: .2f}")
open(os.path.join(os.path.dirname(__file__), 'cache', 'np2b_small_eta.txt'), 'w').write('\n'.join(out) + '\n')
