"""n-sum check behind h13_laplace.py: m^4 C_{2,nu} on one large-J state at nu/Lam = 0.3, 0.7, 0.85, 0.95, 0.995 with
 (a) skelexact3.Laplace3 (double-precision Gauss tail weights) for several (N0, nq, L-factor), and
 (b) h13_laplace.Laplace3H (Arb tail weights) for several (N0, nq, L-factor).
Usage: python3 t_h13_nsum.py Dg Df mfl J  -> cache/t_h13_nsum.txt (appended)"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import h13_exact as H
import skelexact3
skelexact3.mack = H.mack
from skelexact3 import Laplace3
from h13_laplace import Laplace3H
Dg, Df, mfl, J = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
st = H.exact_state(Dg, mfl, J); Lam = H.Lam_of(Dg)
nus = Lam*np.array([0.3, 0.7, 0.85, 0.95, 0.995])
lines = [f"# t_h13_nsum Dg={Dg} Df={Df:g} m/Lam={st['mfl']:.4f} J={J}: m^4 C_2 at nu/Lam = 0.3, 0.7, 0.85, 0.95, 0.995"]
for cls, name in ((Laplace3, 'double-weights'), (Laplace3H, 'arb-weights')):
    for N0, nq, Lf in ((120, 160, 1.0), (120, 320, 0.25), (200, 320, 1.0), (320, 320, 0.5), (400, 240, 1.0), (800, 320, 1.0), (1600, 320, 1.0)):
        Ld = max(4.0*N0, 0.5*st['Delta']**2/max(2*Df - 1.5, 0.5))*Lf
        L3 = cls(st['Delta'], J, Df, ks=(2,), N0=N0, nq=nq, L=Ld, numax=Lam)
        v = [L3.value(2, nu) for nu in nus]
        s = f"  {name:15s} N0={N0:5d} nq={nq} Lf={Lf}: " + ' '.join('%+.12e' % float(a.mid()) for a in v)
        print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 't_h13_nsum.txt'), 'a').write('\n'.join(lines) + '\n')
