import sys, numpy as np
from lp_light import solve
g=int(sys.argv[1]); n=int(sys.argv[2])
for kappa in [float(x) for x in sys.argv[3:]]:
    L, lam = solve(g, n, kappa)
    if L is None: print(g,n,kappa,'fail'); continue
    print(f"g={g} n={n} kappa={kappa:.3g}: L={L:.4e} L/kappa={L/kappa:.4g} heavy-tau={4*np.sqrt(3/L):.4e} light-coef 16kappa/L={16*kappa/L:.4g}", flush=True)
