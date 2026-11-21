"""hc_jk: jackknife error of per-state lambda from hc_chain outputs <mode>_<prefix>jk0..7 (w2/z4 modes).  Usage: python3 hc_jk.py mode Dg prefix J"""
import os, sys, pickle
import numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
mode, Dg, pre, J = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
L = []
for i in range(8):
    f = os.path.join(C, 'hc_chain_%s_%sjk%d.pkl' % (mode, pre, i))
    if not os.path.exists(f):
        continue
    o = pickle.load(open(f, 'rb'))[Dg]
    R = sorted([r for r in o if r['J'] == J], key=lambda r: (round(r['mfl'], 6), str(r['DF']), r['lam']))
    L.append(([round(r['mfl'], 4) for r in R], [r['lam'] for r in R]))
ms = L[0][0]; V = np.array([l[1] for l in L]); n = len(V)
e = np.sqrt((n - 1)/n*((V - V.mean(0))**2).sum(0))
i = int(np.argmin(V.mean(0)))
print(f"J={J}: {n} jackknife runs; min lambda {V.mean(0)[i]:+.4f} +- {e[i]:.4f} at m/Lam={ms[i]}; max err over states {e.max():.4f}")
