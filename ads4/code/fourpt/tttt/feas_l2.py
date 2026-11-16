"""Feasibility of <TTTT> blocks via the e2e/fourpt symbolic route: time and size of single word pairs at jet order 8."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sympy as sp
from embed2 import word4
for (wL, wR) in (('HH', 'HH'), ('ACH', 'ACH'), ('AACC', 'HH'), ('AACC', 'AACC')):
    t = time.time()
    Phi, S, G = word4(wL, wR, 8)
    n = len(sp.Add.make_args(Phi))
    print(f"{wL} x {wR}: {n} terms, {time.time()-t:.0f} s", flush=True)
