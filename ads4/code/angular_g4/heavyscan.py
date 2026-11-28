import sys, time, numpy as np, lp2
g = int(sys.argv[1])
for N in [int(x) for x in sys.argv[2].split(',')]:
    for useM in (0, 1):
        t0 = time.time()
        o = lp2.solve(g, N, N if useM else 0, None, Jlin=8 * g, nJlog=150, nt=50)
        if o is None: print(g, N, useM, 'fail', flush=True); continue
        print(f"g={g} K={N} M={'N' if useM else 0}: tau={o['tau']:.5e} tau*g^3={o['tau']*g**3:.1f} nu={o['nu']:.3g} [{time.time()-t0:.0f}s]", flush=True)
