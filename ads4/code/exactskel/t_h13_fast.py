"""Check Laplace3H.value_fast (transposed Horner) against Laplace3H.value on one state (all even k, a few nu).
Usage: python3 t_h13_fast.py Dg Df mfl J  -> printed, appended to cache/t_h13_fast.txt"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import h13_exact as H
from h13_laplace import Laplace3H
Dg, Df, mfl, J = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
st = H.exact_state(Dg, mfl, J); Lam = H.Lam_of(Dg)
L3 = Laplace3H(st['Delta'], J, Df, ks=(2, 4, 6, 8, 10), N0=200, nq=320, numax=Lam)
lines = []
for k in (2, 4, 6, 8, 10):
    for nu in Lam*np.array([0.0, 0.3, 0.8, 0.99]):
        t = time.time(); a = L3.value(k, nu); t1 = time.time() - t
        t = time.time(); b = L3.value_fast(k, nu); t2 = time.time() - t
        s = (f"Dg={Dg} Df={Df:g} m/Lam={st['mfl']:.3f} J={J} k={k} nu/Lam={nu/Lam:.2f}: value {float(a.mid()):+.15e} (rad {float(a.rad()):.0e}, {t1:.2f}s)"
             f"  fast {float(b.mid()):+.15e} (rad {float(b.rad()):.0e}, {t2:.3f}s)  rel diff {abs(float((a - b).mid()))/max(abs(float(a.mid())), 1e-300):.1e}")
        print(s, flush=True); lines.append(s)
open(os.path.join(HERE, 'cache', 't_h13_fast.txt'), 'a').write('\n'.join(lines) + '\n')
