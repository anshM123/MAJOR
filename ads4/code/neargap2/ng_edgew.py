"""ng_edgew (NEARGAP2.md §1.6): sign-free weight of the light-cone-edge windows x in [0.8, 0.9] and (0.9, Lam/m] for near-gap states
(hybrid weights, ee1000_common.hyb, 1024 nodes), and the cost of |Gamma| <= G there.  Usage: python3 ng_edgew.py Dg J:mfl[:G] ..."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ng_chain as NC
import numpy as np
E, Q = NC.E, NC.Q
Dg = float(sys.argv[1]); L = []
for spec in sys.argv[2:]:
    p = spec.split(':'); J = int(p[0]); mfl = float(p[1]); G = float(p[2]) if len(p) > 2 else 1e4
    hs = ['+-'] if J % 2 else ['++', '+-']
    for h in hs:
        st = E.hyb(Dg, J, mfl, h, 1024); x, m, u = st['x'], st['m'], st['u']
        env = Q.P4.env if J % 2 else Q.A.env
        w = (np.abs(st['p2'])*env(h, x, 2) + u*np.abs(st['p3'])*env(h, x, 3))/m**2/st['Ns']
        a = w[(x >= 0.8) & (x <= 0.9)].sum(); b = w[x > 0.9].sum()
        L.append(f"J={J} {h} m/Lam={mfl}: lam_hyb {st['base']/st['Ns']:+.4f}; W[0.8,0.9] = {a:.2e}, W(0.9, {x.max():.3f}] = {b:.2e} per c2;"
                 f" cost at G = {G:.0f}: {a*G:.3f} / {b*G:.3f} N_s")
print("\n".join(L))
open(os.path.join(NC.NGC, 'ng_edgew_%g.txt' % Dg), 'a').write("\n".join(L) + "\n")
