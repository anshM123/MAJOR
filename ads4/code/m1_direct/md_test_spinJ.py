"""md_test_spinJ.py -- unit tests of md_spinJ: (j0) J = 0 against md_dquad.Q2 (rd_core / repro_engine EE contraction, independent
implementation); (sym) the k = 2 integrand at J = 2, 4 is real and finite on D'; boundary exponents at the B zero (t -> inf)."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import md_dquad as MD
import md_spinJ as SJ
np.seterr(all='ignore')

rng = np.random.default_rng(11)
n = 7
pts = [rng.uniform(0.2, 3, n), rng.uniform(0.05, 0.95, n), rng.uniform(-0.9, 0.9, n), rng.uniform(0.2, 3, n),
       rng.uniform(0.05, 0.95, n), rng.uniform(-1.3, 1.3, n), rng.uniform(0.1, 3.0, n)]
for Delta in (10.0, 23.7):
    o = P.point(*pts, Delta)
    q2 = MD.Q2(o['g'], Delta)
    r = SJ.ratio(o['g'], Delta, 0)[0, 0]
    print(f'j0 Delta={Delta}: ratio md_spinJ/Q2 =', np.array2string(r/q2, precision=12))
for J in (2, 4):
    for Delta in (11.0, 23.7):
        o = P.point(*pts, Delta)
        pre = SJ.prefh(o['g'], Delta, J)*SJ.jfactor(o['g'], J)[None, None]
        print(f'J={J} Delta={Delta}: nstruct {pre.shape[:2]}, max |Im|/|Re| = {np.max(np.abs(pre.imag))/np.max(np.abs(pre.real)):.1e}')
# decay toward the B zero (t -> inf), per J
base = [1.0, 0.4, 0.3, 1.2, 0.55, 0.3, 1.1]
for J in (0, 2, 4):
    for Delta in (12.0, 20.0):
        rows = []
        for which in (0, 3):
            vals = []
            for tt in (2*Delta, 4*Delta):
                p = [np.array([v]) for v in base]; p[which] = np.array([tt])
                o = P.point(*p, Delta)
                rr = SJ.ratio(o['g'], Delta, J)
                vals.append(o['L'][0] + np.log(np.max(np.abs(rr))))
            rows.append((vals[1] - vals[0])/(2*Delta))
        print(f'J={J} Delta={Delta}: d log|Phi|/dt at large t_x, t_y = {rows[0]:+.4f}, {rows[1]:+.4f};  -(1-(5+J)/Delta) = {-(1-(5+J)/Delta):+.4f}')
