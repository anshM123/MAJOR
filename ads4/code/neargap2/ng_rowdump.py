"""ng_rowdump (NEARGAP2.md): Gamma (c2) of one row at Laplace orders K = 0, 1, 2 (chain conventions, normalisation a per K).
Usage: python3 ng_rowdump.py J Delta  -> stdout (and cache/ng_rowdump_J<J>_D<Delta>.txt)"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ng_chain as NC
import numpy as np
Q = NC.Q
J, D = int(sys.argv[1]), float(sys.argv[2])
gs = {}
for K in (0, 1, 2):
    NC.KK = K
    gs[K] = (Q.odd_grid(J) if J % 2 else Q.even_grid(J))[D]
keys = ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in (('++', '+-', 'off') if J >= 4 else ('++',))]
L = [f"# J={J} Delta={D:g} m={gs[2]['m']:.3f}  a(K=0,1,2) = {gs[0]['a']:.8g} {gs[1]['a']:.8g} {gs[2]['a']:.8g};  Gamma: K2 [K1, K0]"]
xs = gs[2]['xs']
for i, x in enumerate(xs):
    s = f"x={x:.4f}"
    for k in keys:
        v = [g[k][int(np.argmin(np.abs(g['xs'] - x)))] for g in (gs[2], gs[1], gs[0])]
        s += f" | {k}: {v[0]:+9.1f} [{v[1]:+9.1f} {v[2]:+9.1f}]"
    L.append(s)
print("\n".join(L))
open(os.path.join(NC.NGC, 'ng_rowdump_J%d_D%g.txt' % (J, D)), 'w').write("\n".join(L) + "\n")
