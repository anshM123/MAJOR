"""
f5: coherent (smooth, sign-definite) graviton deviations on the route-A near-gap action -- a less conservative companion of the
worst-sign tolerances of f3.  For a deviation family  F -> F_ex + eps * sum_k sum_{x_i > x0} c_k T_k(p_i),  T = dens+curv (H2-type) or
rem (H4-type), and pattern c_k in {1 (all k), (-1)^k, [k>=4]}, report the interval (eps_-, eps_+) of eps for which every state stays
positive (n = 1 kept; all exact-skeleton states, both Df).
Usage: python3 f5_coherent.py [Dg ...] -> cache/f5_coherent.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from spf_common import NearGap, load_states
from common import KS

PATS = {'all': {k: 1.0 for k in KS}, 'alt': {k: (-1.0)**k for k in KS}, 'k>=4': {k: float(k >= 4) for k in KS}}
X0S = (0.0, 0.3, 0.5, 0.7)
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def main(Dgs):
    pr("# f5: coherent deviation tolerances (eps_-, eps_+) keeping n = 1, route A, NF02; T = dens+curv or rem, x > x0")
    for Dg in Dgs:
        NG = NearGap(Dg)
        acc = {}
        for h, DF, r0, ref in load_states(Dg):
            x, P = NG.pieces(h, r0)
            for x0 in X0S:
                m = x > x0
                for pat, c in PATS.items():
                    for T in ('dens', 'rem'):
                        if T == 'dens':
                            G = sum(c[k]*(P[k]['dens'] + P[k]['curv'])[m].sum() for k in KS)
                        else:
                            G = sum(c[k]*P[k]['rem'][m].sum() for k in KS)
                        lo, hi = acc.get((x0, pat, T), (-np.inf, np.inf))
                        if G > 0:
                            lo = max(lo, -ref['Fex']/G)
                        elif G < 0:
                            hi = min(hi, -ref['Fex']/G)
                        acc[(x0, pat, T)] = (lo, hi)
        pr(f"  Delta_g = {Dg}:")
        for T in ('dens', 'rem'):
            for pat in PATS:
                pr(f"    {T:4s} {pat:5s}: " + "  ".join(f"x>{x0}: ({acc[(x0, pat, T)][0]:+.3g}, {acc[(x0, pat, T)][1]:+.3g})" for x0 in X0S))
    open(os.path.join(HERE, 'cache', 'f5_coherent.txt'), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
