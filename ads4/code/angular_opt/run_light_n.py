import sys, numpy as np
from lp_light import solve
g=int(sys.argv[1]); p=float(sys.argv[2])
kappa=float(g)**p
for n in [int(x) for x in sys.argv[3:]]:
    L, lam = solve(g, n, kappa)
    if L is None: print(g,n,p,'fail', flush=True); continue
    print(f"g={g} p={p} n={n}: L={L:.4e} L/g^2={L/g**2:.3f} tau_heavy={4*np.sqrt(3/L):.4e} tau*g={4*np.sqrt(3/L)*g:.3f} lightcoef={16*kappa/L:.4g}", flush=True)
