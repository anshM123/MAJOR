import sys, math, time
import numpy as np
sys.path.insert(0, '.')
from common2 import n_of_L
from fastS import S_vec, d_of
from scipy.special import loggamma
L = float(sys.argv[1]); nu = float(sys.argv[2]); a = float(sys.argv[3]); b = int(sys.argv[4]); JJ = int(sys.argv[5])
s = 2*b + 1
n = n_of_L(L, a, b); d = d_of(n, a, s)
al = 2.5 + 1j*nu
js = np.arange(JJ)
S, lg = S_vec(n, a, s, al - js)
logSh = np.log(np.abs(S)) + lg + (2.5 - js)*math.log(d)
argSh = np.angle(S) + nu*math.log(d)   # phase of d^{i nu}
# c_j
lc = np.zeros(JJ, dtype=complex)
a0 = 0.5 - 1j*nu; b0 = 1 - 2j*nu
for j in range(1, JJ):
    lc[j] = lc[j-1] + 2*np.log(a0 + j - 1) - np.log(b0 + j - 1) - math.log(j)
lB = loggamma(2j*nu) - 2*loggamma(0.5 + 1j*nu)
# K = 2 Re B sum_j c_j S(alpha - j) = 2Re B d^{-alpha} sum c_j d^j S^_j
terms = np.exp(lc + js*math.log(d) + logSh + 1j*argSh)
part = np.cumsum(terms)
pref = np.exp(lB - al*math.log(d))
print(f'n={n} d={d:.5e} L={-math.log(d):.4f} nu={nu} a={a} b={b}')
for j in list(range(0, JJ, max(1, JJ//20))):
    print(f' j={j:3d} log|S^_j|={logSh[j]:8.3f}  log|c_j|d^j={lc[j].real + j*math.log(d):8.3f}   |partial K| = {abs(2*(pref*part[j]).real):.4e}')
print('K_true =', 2*(pref*part[-1]).real, '  |B d^-alpha| =', abs(pref))
print('termwise sum |B||c_j|d^{j-5/2}|S^_j| =', abs(pref)*np.sum(np.abs(terms)), ' log =', math.log(abs(pref)*np.sum(np.abs(terms))))
