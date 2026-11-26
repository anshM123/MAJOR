"""End-to-end sanity: on random points of every new coupled rectangle,
|F_K - R| must stay below the certified uniform residual M, and F_K >= B_R - M.
Also: true physical states (integer n, a>=998) lie in the rectangle (y<=1), and
F_17199/(1+qh)^4 at those states is positive and agrees with F_K within the
recorded physical tail.  Independent evaluator (rc.F_K), not the archived C++."""
import random, sys
from fractions import Fraction as Fr
from math import isqrt
import rc

random.seed(12345)
NEW = ['cell_5305_6041', 'cell_6041_6677', 'cell_6677_7187', 'cell_7187_7596', 'cell_7596_7929',
       'cell_7929_8204', 'cell_8204_8434', 'cell_8434_8628', 'cell_8628_8794', 'cell_8794_8936', 'cell_8936_9000']
npts = int(sys.argv[1]) if len(sys.argv) > 1 else 12
worst = 0.0
for stem in NEW:
    z, l, u, Hl, Hu = rc.cell(stem); K = z['K']
    led, suf = rc.ledger(stem)
    pts = [(Fr(0), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(0)), (Fr(1), Fr(1))]
    pts += [(Fr(random.randrange(10**6), 10**6), Fr(random.randrange(10**6), 10**6)) for _ in range(npts)]
    # bias some y near 1 (cap edge, where the prefix is least tame)
    pts += [(Fr(random.randrange(10**6), 10**6), 1 - Fr(random.randrange(10**3), 10**6)) for _ in range(4)]
    for i, b in enumerate((0, 4)):
        P = led['reference'][i]; M = Fr(led['results'][i]['uniform_residual']); BR = Fr(led['results'][i]['ref_lower'])
        mx = Fr(0); mn = None
        for x, y in pts:
            c = l + (u - l) * x; h = rc.Hsec(c, l, u, Hl, Hu) * y
            f = rc.F_K(c, h, K, b)
            d = abs(f - rc.Rval(P, x, y)); mx = max(mx, d)
            mn = f if mn is None else min(mn, f)
        r = float(mx / M); worst = max(worst, r)
        print(f'{stem} b={b} K={K} max|F-R|={float(mx):.3e} M={float(M):.3e} ratio={r:.3e} minF={float(mn):.4e} BR-M={float(BR-M):.4e}',
              'OK' if mx <= M and mn >= BR - M else 'FAIL', flush=True)
    # physical states: pick J, a giving c in cell
    for b in (0, 4):
        s = 2 * b + 1; q = b * (b + 1); cnt = 0; tries = 0
        while cnt < 3 and tries < 20000:
            tries += 1
            a = random.choice([998, 998, 999, 1000, 1003, 1100, 1500, 3000])
            # choose n so that c in [l,u]: c=n(n+s)/(N(N+s)), N=n+a
            import math
            sl, su = math.sqrt(float(l)), math.sqrt(float(u))
            n = random.randrange(max(1000, int(a*sl/(1-sl)) - 5), int(a*su/(1-su)) + 6)
            N = n + a; h = Fr(1, N * (N + s)); c = n * (n + s) * h
            if not (l <= c <= u) or n + b < 1001:
                continue
            cnt += 1
            y = h / rc.Hsec(c, l, u, Hl, Hu)
            assert 0 <= y <= 1, ('physical y>1', stem, n, a, float(y))
            fK = rc.F_K(c, h, K, b); fP = rc.F_K(c, h, min(17199, n), b)  # U_k=0 for k>n (termination)
            ph = fP / (1 + q * h) ** 4
            print(f'  phys b={b} n={n} a={a} y={float(y):.6f} F_K={float(fK):.6e} F_17199/(1+qh)^4={float(ph):.6e} |diff|={float(abs(fK-fP)):.2e} tail={float(Fr(z["physical_tail"])):.2e}',
                  'OK' if abs(fK - fP) <= Fr(z['physical_tail']) and ph > 0 else 'FAIL', flush=True)
print('worst ratio max|F-R|/M =', worst)
