"""
hs_refzero (HIGHSPIN.md §2): the reference-zero diagnostic on every K2 row of the Delta_g = 1000 chain for J = 17 ... 48.
For each row: the x at which the Delta_ref = 2e4 reference's leading partial wave has its first zero (x0 = j_h m_ref / m_row), the
K2 node nearest to it, the excess |r2/r0| of the reference record there (EE and C3H; > 0.1 = the K2 series of the ratio fails), and
whether the chain's Hampel rule (w2_chain, W2_DESPIKE=1) drops that node.
Usage: python3 hs_refzero.py   -> cache/hs_refzero.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_core as H
import numpy as np

L = []
for J in list(range(17, 41)) + [48]:
    os.environ['W2_DESPIKE'] = '1'
    H.W.DESPIKE_LOG.clear()
    try:
        g = H.Q.odd_grid(J) if J % 2 else H.Q.even_grid(J)
    except Exception as ex:
        L.append(f"J={J}: grid error {ex!r}"[:150]); continue
    desp = {(round(D, 1), round(x, 4)) for (JJ, D, xs) in [(a[0], a[1], a[2]) for a in H.W.DESPIKE_LOG] for x in xs} \
        if H.W.DESPIKE_LOG and len(H.W.DESPIKE_LOG[0]) == 3 else set()
    for D in sorted(g):
        x0 = H.ref_zero_x(J, D)
        rd = {kind: {round(x, 4): ex for x, v, ex in H.refdiag(J, D, kind)} for kind in (('EE', 'C3H'))}
        allx = sorted(set(rd['EE']) | set(rd['C3H']))
        if not allx:
            L.append(f"J={J:2d} row {D:6.0f}: no K2 reference records"); continue
        xn = min(allx, key=lambda v: abs(v - x0))
        ex = max(rd['EE'].get(xn, 0), rd['C3H'].get(xn, 0))
        kept = xn in set(np.round(g[D]['xs'], 4))
        L.append(f"J={J:2d} row {D:6.0f}: x0 = {x0:.4f}; nearest K2 node {xn:.4f} (|dx| = {abs(xn - x0):.4f}); reference excess |r2/r0| = "
                 f"{ex:.3f}{'  <-- FAILS' if ex > 0.1 else ''}; in the chain after Hampel: {'kept' if kept else 'dropped'}")
print("\n".join(L))
open(os.path.join(H.HSC, 'hs_refzero.txt'), 'w').write("\n".join(L) + "\n")
