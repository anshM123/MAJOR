"""md_decay.py -- decay rate of the spin-J EE k=2 integrand toward the B zero (t_x or t_y -> inf) on D', several base points."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import md_spinJ as SJ
np.seterr(all='ignore')
rng = np.random.default_rng(3)
Delta = float(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
for J in Js:
    out = []
    for trial in range(3):
        base = [rng.uniform(0.5, 2), rng.uniform(0.2, 0.8), rng.uniform(-0.7, 0.7), rng.uniform(0.5, 2), rng.uniform(0.2, 0.8),
                rng.uniform(-1, 1), rng.uniform(0.3, 2.8)]
        for which in (0, 3):
            vals = []
            for tt in (5*Delta, 8*Delta):
                p = [np.array([v]) for v in base]; p[which] = np.array([tt])
                o = P.point(*p, Delta)
                rr = SJ.ratio(o['g'], Delta, J)
                vals.append(o['L'][0] + np.log(np.max(np.abs(rr))))
            out.append((vals[1] - vals[0])/(3*Delta))
    print(f'J={J} Delta={Delta}: rates {np.round(out, 4)};  guess -(1-(5+J)/Delta) = {-(1-(5+J)/Delta):+.4f}; -(1-5/Delta) = {-(1-5/Delta):+.4f}', flush=True)
