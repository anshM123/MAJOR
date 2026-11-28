import numpy as np, flat, fastflat
from mpmath import mpf
flat.mp.dps = 25
for n in [0, 1, 2, 7]:
    for u in [0.0, 0.3, 3.0, 40.0, 900.0, 20000.0, 1e6]:
        f = fastflat.Pi0(n, np.array([u]))
        for d in range(3):
            ex = float(flat.Pi0(n, mpf(u), d))
            print(n, u, d, f'{float(f[d][0]):.15e}', f'{ex:.15e}', f'rel {abs(float(f[d][0])/ex-1):.1e}')
