"""Light-robust dual LP at gap g with safety margin delta, rationalised output.
Constraints (sampled): heavy  W+Lam4 >= delta W, W-Lam0 >= delta W ; light |Lam| <= (1-delta) kappa W.
Saves ell_k = lambda_k D_k^2 as exact rationals (denominator 2^40) in cache/robust_g{g}_n{n}.json"""
import numpy as np, sys, json
from fractions import Fraction as Fr
from scipy.optimize import linprog
from lp_light2 import rows
g = int(sys.argv[1]); n = int(sys.argv[2]); kappa = float(sys.argv[3]); delta = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-5
A, bnd = rows(g, n, kappa, nJ=600, nth=40, nlight=40)
# heavy rows have bnd=1 -> 1-delta ; light rows bnd=kappa -> kappa(1-delta)
bnd = bnd*(1-delta)
sc = np.abs(A).max(axis=0); sc[sc == 0] = 1
c = np.zeros(n+1); c[2] = 1/sc[2]
res = linprog(c, A_ub=A/sc, b_ub=bnd, bounds=[(None, None)]*(n+1), method='highs')
assert res.status == 0, res.message
lam = res.x/sc
D2 = [Fr((k+1)*(k+2)*(k+3), 6)**2 for k in range(n+1)]
ell = [Fr(int(round(float(lam[k])*float(D2[k])*2**40)), 2**40) for k in range(n+1)]
L = -ell[2]/D2[2]
print('g', g, 'n', n, 'kappa', kappa, 'L', float(L), 'tau_heavy', 4*np.sqrt(3/float(L)), 'lightcoef', 16*(kappa-1)/float(L))
json.dump(dict(g=g, n=n, kappa=str(Fr(kappa)), delta=delta, ell=[str(x) for x in ell]), open(f'cache/robust_g{g}_n{n}.json', 'w'))
