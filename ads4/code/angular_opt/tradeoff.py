"""Trade-off curve: light-state coefficient Lambda_l = 16(kappa-1)/L versus heavy bound tau_h = 4 sqrt(3/L),
with light spin-0/2 single traces included as exact constraints.  tau^2 <= tau_h^2 (1+D/3) + Lambda_l E'_light."""
import sys, numpy as np
from lp_light2 import solve
g = int(sys.argv[1]); n = int(sys.argv[2])
K = 729/4800*(g-2)**2*(1-49/(2*g+9)**2)
print(f'# g={g} n={n} paper: tau_*={4*np.sqrt(3/K):.4e} (tau g={4*np.sqrt(3/K)*g:.2f}), light coef ~1120')
for kappa in [float(x) for x in sys.argv[3:]]:
    L, lam = solve(g, n, kappa)
    if L is None: print(g, n, kappa, 'fail', flush=True); continue
    th = 4*np.sqrt(3/L)
    print(f"g={g} n={n} kappa={kappa:.3g}: L/K={L/K:.2f} tau_h={th:.4e} tau_h*g={th*g:.3f} log_g(1/tau_h)={np.log(1/th)/np.log(g):.3f} Lambda_l={16*(kappa-1)/L:.4g}", flush=True)
