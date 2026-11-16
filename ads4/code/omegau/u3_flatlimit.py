"""u3: sanity check of the flat limit of the skeleton AdS action of Omega^U.

At fixed (m/Lam, J, h) (fixed bhat) and Dg -> infinity, the skeleton action F_skel = Fh + S (hybrid + signed scalar-skeleton remainder,
all k = 2..7 and forward terms) must tend to the flat heavy action of the stripped functional (F_g3 / M02 / M01), with corrections O(1/Lam^2).
Printed: flat, F_skel, Lam^2 (F_skel - flat)/Ns for Dg = 100, 1000, 10^4, 10^5 (states with all nodes controlled; else flagged).
Output: cache/u3_<Df>.log
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from common import load, KS, Transforms, state, flat_unstripped_parts, flat_parts, gl01
from adsaction import action
import remk

pn, pw = gl01(64); pf, wf = gl01(400)
out = []
def pr(s):
    print(s, flush=True); out.append(s)
W0 = Transforms({k: np.array([1.0]) for k in KS})
pr(f"# u3 skeleton {remk.TAG}: Lam^2 (F_skel - F_flat)/Ns -> const as Dg -> infinity at fixed (m/Lam, J)")
for name in ('M02', 'Fg3'):
    F, PH = load(name); TR = Transforms(PH)
    for (mfl, J, h) in ((1.2, 0, '++'), (1.5, 2, '++'), (2.0, 6, '++'), (1.2, 4, '+-'), (1.5, 6, '+-'), (3.0, 20, '+-'), (6.0, 40, '++')):
        u = 1/mfl**2
        Jsu, parts, FWu = flat_unstripped_parts(F, PH, h, u, J, pf, wf)
        Jss, B2, B3, FWs = flat_parts(F, h, u, J)
        flat = B2[-1] + B3[-1] + FWs[-1]
        s = f"  {name} {h} m/Lam={mfl} J={J}: flat = {flat:+.6f} | F_skel - flat [Lam^2 x that / Ns] at Dg ="
        for Dg in (100, 1000, 10**4, 10**5):
            Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
            st = state(Dg, Lam, mfl, J)
            if st is None:
                s += f" Dg={Dg}: unphys |"; continue
            _, W, P0 = W0.regge(Lam, st['beta'])
            sk = remk.skeleton_state(st, J, pn*np.sqrt(u), W, P0)
            r = action(F, PH, TR, h, st, J, mfl, Lam, pn, pw, sk, parts[-1], FWu[-1], (B2[-1], B3[-1], FWs[-1]))
            Fs = r['Fh'] + r['Sall']
            flag = '' if r['U'] == 0 else '*'
            s += f" {Dg}: {Fs - flat:+.2e} [{Lam**2*(Fs - flat)/r['Ns']:+7.1f}]{flag} |"
        pr(s)
open(os.path.join(HERE, 'cache', f'u3_{remk.TAG}.log'), 'w').write("\n".join(out) + "\n")
