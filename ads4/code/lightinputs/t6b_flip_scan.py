"""t6b: for every frame mode (q1,q3) of l = 1, 2: the parity-odd part of its helicity-flip response
(M_{++} - M_{--})/2 and of its preserved response (M_{+-} - M_{-+})/2 (laplace.M_table, exact).
If every odd flip part vanishes, NO parity-odd kernel at a pair can produce helicity-flip (c_same cbar_opp) pieces at leading order."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e')); os.chdir(os.path.join(HERE, '..', 'e2e'))
import sympy as sp
from laplace import M_table, nu, a, g1
for l in (1, 2):
    T = M_table(l)
    nzf = []; nzp = []
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            M = lambda h1, h3: T[(q1, q3, h1, h3)][0] + g1*T[(q1, q3, h1, h3)][1]
            fo = sp.simplify(M(1, 1) - M(-1, -1)); po = sp.simplify(M(1, -1) - M(-1, 1))
            fe = sp.simplify(M(1, 1) + M(-1, -1)); pe = sp.simplify(M(1, -1) + M(-1, 1))
            if fo != 0: nzf.append((q1, q3))
            if po != 0: nzp.append((q1, q3))
            print(f"l={l} q=({q1:+d},{q3:+d}): odd-flip {'0' if fo==0 else 'NONZERO'}  odd-pres {'0' if po==0 else 'nonzero'}"
                  f"  even-flip {'0' if fe==0 else 'nonzero'}  even-pres {'0' if pe==0 else 'nonzero'}")
    print(f"l={l}: modes with nonzero ODD FLIP response: {nzf}")
    print(f"l={l}: modes with nonzero odd pres response: {nzp}")
