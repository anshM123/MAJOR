"""Odd k from even k: cubic/quartic interpolation in k of rho_k (and of m^{2k}C_k / E_k) from k = 2,4,6,8(,10), tested against the exact odd-k
moments (Gamma-ratio continuation, skelexact.py) on Delta_g = 100 e1 states."""
import sys, pickle
sys.path.insert(0, '.')
import numpy as np
from skelexact3 import Laplace3
from scipy.special import j0, eval_legendre
db = pickle.load(open('cache/e1_100_Df3.pkl', 'rb'))
Ek = lambda k, x: (2 - x*x)/(1 - x*x)**(1 + k/2)
worst = {}
for key in [(1.0, 0), (1.04, 4), (1.08, 4), (1.1, 2), (1.2, 2), (1.3, 6), (1.5, 10), (2.0, 4)]:
    r = db[key]
    L3 = Laplace3(r['Delta'], key[1], 3.0, ks=(2, 4, 6, 8, 10), numax=r['nus'].max())
    x = r['x']
    rho = {}
    for k in (2, 4, 6, 8, 10):
        v = np.array([float(L3.value(k, nu).mid()) for nu in r['nus']])
        hyb = Ek(k, x)*eval_legendre(key[1], 1 - 2*x*x) + 2*(0)   # curvature part cancels in the interpolation test below
        rho[k] = v/Ek(k, x)
    for kind, ks in (('cubic 2-8', (2, 4, 6, 8)), ('quartic 2-10', (2, 4, 6, 8, 10))):
        for ko in (3, 5, 7):
            # Lagrange interpolation in k of v_k/E_k
            est = sum(rho[ka]*np.prod([(ko - kb)/(ka - kb) for kb in ks if kb != ka]) for ka in ks)
            exact = r['ex'][ko]/Ek(ko, x)
            d = np.abs(est - exact)
            worst[(kind, ko)] = max(worst.get((kind, ko), 0), d.max())
            print(f"{key} {kind} k={ko}: max|interp - exact| (units of E_k) = {d.max():.2e} at x = {x[np.argmax(d)]:.3f}; at x<0.9: {d[x < 0.9].max():.2e}")
print(worst)
