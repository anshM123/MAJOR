import sys, numpy as np
from lp_light import solve
g=int(sys.argv[1]); n=int(sys.argv[2])
for p in [2, 2.5, 3, 3.5, 4, 5, 6]:
    kappa = float(g)**p
    L, lam = solve(g, n, kappa)
    if L is None: print(g,n,p,'fail', flush=True); continue
    print(f"g={g} n={n} p={p} kappa={kappa:.3g}: L={L:.4e} log_g L={np.log(L)/np.log(g):.3f} tau_heavy*g={4*np.sqrt(3/L)*g:.3f} lightcoef={16*kappa/L:.4g}", flush=True)
