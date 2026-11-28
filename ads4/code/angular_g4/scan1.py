import sys, numpy as np, time
from lp1 import solve, grid
for g in [int(x) for x in sys.argv[1].split(',')]:
    S = grid(g, 12 * g, 2000 * g, nJlog=60, nt=40)
    for N in [int(x) for x in sys.argv[2].split(',')]:
        K = max(N, 12)
        t0 = time.time()
        c2, res, _ = solve(g, K, N, S)
        if c2 is None: print(g, N, 'status', res.status, flush=True); continue
        tau = 4 * np.sqrt(max(c2, 0))
        print(f'g={g} K={K} N={N}: tau={tau:.5e} tau g^3={tau*g**3:.2f} tau g^4={tau*g**4:.1f}  ({time.time()-t0:.0f}s)', flush=True)
