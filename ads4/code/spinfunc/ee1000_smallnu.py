"""
ee1000_smallnu (EE_DG1000_LARGEJ.md §5): how large (negative) a beyond-leading-Regge graviton correction confined to the small-nu window
nu < NU_C (where neither K2 nor the direct method controls Gamma at Delta_g = 1000: §3) must be to flip a state.  Per state (J, m/Lam):
  W_<(h) = sum over p-nodes with nu = p Lam < NU_C of [psi2 w env_{2,h} + psi3 w u env_{3,h}]/(m^2 N_h)   (uniform Gamma = C in c2 units,
  same C for C^_2 and C^_3),  C*(J, m) = -lambda0/W_<  (uncorrected band value lambda0; with the computed rows' corrections where given).
Output: cache/ee1000_smallnu.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_common as E
Q = E.Q
DG = 1000; Lm = E.lam_of(DG)
NU_C = float(os.environ.get('EE1000_NUC', '50'))
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
pr(f"# small-nu window nu < {NU_C}: weight W_< per unit uniform Gamma (c2 units) / N_h, and the failure threshold C* = -lambda0/W_<")
pr("| J | m/Lam | lambda0 (min h) | W_< (++, +-) | C* (min over h) |")
pr("|---|---|---|---|---|")
for J in (16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128, 160, 200, 256, 300):
    hs = ('+-',) if J % 2 else ('++', '+-')
    for mfl in sorted(set([E.mgap(DG, J), 1.1, 1.2, 1.3, 1.5, 2.0])):
        if mfl < E.mgap(DG, J) - 1e-9:
            continue
        npn = 256 if J < 100 else 768
        Ws, lam, Cs = [], [], []
        for h in hs:
            st = E.hyb(DG, J, mfl, h, npn)
            x, m, u = st['x'], st['m'], st['u']
            sel = st['pn']*Lm < NU_C
            w = ((st['p2']*Q.A.env(h, x, 2) + st['p3']*u*Q.A.env(h, x, 3))/m**2/st['Ns'])[sel].sum()
            l0 = st['base']/st['Ns']
            Ws.append(w); lam.append(l0); Cs.append(-l0/w if w > 0 else np.inf)
        pr(f"| {J} | {mfl:.3f} | `{min(lam):+.3f}` | " + ", ".join(f"{w:.2e}" for w in Ws) + f" | `{min(Cs):.0f}` |")
open(os.path.join(E.C, 'ee1000_smallnu_%g.txt' % NU_C), 'w').write("\n".join(lines) + "\n")
