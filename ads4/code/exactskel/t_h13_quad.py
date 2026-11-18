"""p-quadrature check for h13_exact: recompute S on one state with n and 2n Gauss nodes (n = h13_exact.nodes_for) and with 64 nodes.
Usage: python3 t_h13_quad.py Dg Df mfl J  -> printed (and appended to cache/t_h13_quad.txt)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import h13_exact as H
from common import gl01
Dg, Df, mfl, J = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
st = H.exact_state(Dg, mfl, J)
n0, _ = H.nodes_for((2*J + 1)/st['mfl'])
lines = []
for n in sorted({64, n0, 2*n0}):
    qn, qw = gl01(n)
    r = H.rho_state(Dg, Df, st, 400, pn=qn)
    for h in ('++', '+-'):
        if (h == '++' and J % 2):
            continue
        S, Sabs, _ = H.S_of(h, st['mfl'], r['rho'], qn, qw)
        s = f"Dg={Dg} Df={Df:g} m/Lam={st['mfl']:.4f} J={J} bhat={(2*J+1)/st['mfl']:.1f} {h} nodes={n:4d}: S = {S:+.10e}  sum|.| = {Sabs:.4e}"
        print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 't_h13_quad.txt'), 'a').write('\n'.join(lines) + '\n')
