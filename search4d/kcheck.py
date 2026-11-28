"""Print the K-condition (SSA/WM) violations of a vaidya_scan output file."""
import json, sys, itertools
import numpy as np
d = json.load(open(sys.argv[1]))
r = np.array(d['radii']); t = np.array(d['times']); B = np.array(d['tube_branches']); c = np.array(d['c'])
L = np.array([[x if x is not None else np.inf for x in row] for row in d['L']])
print('radii', np.round(r, 4)); print('times', np.round(t, 4)); print('c', np.round(c, 5))
print('one-sided', d['one_sided'])
Lp = np.minimum(L, c[:, None] + c[None, :]); n = len(r)
for i, j, k in itertools.combinations(range(n), 3):
    a = Lp[i, k] + c[j] - Lp[i, j] - c[k]; b = Lp[i, k] + c[j] - c[i] - Lp[j, k]
    if min(a, b) < -1e-7: print('K1/K2', i, j, k, a, b, 'L', L[i, k], L[i, j], L[j, k], 'branches', B[i, k], B[k, i], B[i, j], B[j, k])
for i, j, k, l in itertools.combinations(range(n), 4):
    a = Lp[i, k] + Lp[j, l] - Lp[i, j] - Lp[k, l]; b = Lp[i, k] + Lp[j, l] - Lp[i, l] - Lp[j, k]
    if min(a, b) < -1e-7: print('K4', i, j, k, l, a, b)
