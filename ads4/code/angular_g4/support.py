"""Print the primal support (LP marginals) of the light-robust dual optimum."""
import sys, numpy as np
from scipy.optimize import linprog
import lp2
g = int(sys.argv[1]); N = int(sys.argv[2]); K = max(N, 16); th = float(eval(sys.argv[3]))
o = lp2.solve(g, K, N, th, ret=True)
B, J, D, Lt, A = o['B'], o['J'], o['D'], o['Lt'], o['A']
x = np.concatenate([[o['nu']], o['mu'], o['lam']])
F = A @ x
print('L', o['L'], 'L/g^2', o['L'] / g**2, 'nu', o['nu'])
# near-active constraints
act_h = np.where((Lt == 0) & (F < 1e-6 * np.abs(A).max(axis=1).clip(1)))[0]
act_l = np.where((Lt == 1) & (F + th < 1e-6 * th))[0]
print('active heavy (b,J,twist):', [(int(B[i]), int(J[i]), round(D[i] - J[i], 3)) for i in act_h][:40])
print('active light (b,J,twist):', [(int(B[i]), int(J[i]), round(D[i] - J[i], 3)) for i in act_l][:40])
print('lam*D^2:', [f'{l*((k+1)*(k+2)*(k+3)/6)**2:.3g}' for k, l in enumerate(o['lam'][:10])])
print('mu:', [f'{m:.3g}' for m in o['mu'][:10]])
