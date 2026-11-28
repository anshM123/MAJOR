import sys, numpy as np, lp2
g = int(sys.argv[1]); N = int(sys.argv[2]); K = max(N, 20)
for p in [float(x) for x in sys.argv[3].split(',')]:
    th = g ** p
    o = lp2.solve(g, K, N, th)
    if o is None: print(g, N, p, 'fail', flush=True); continue
    print(f"g={g} N={N} Theta=g^{p}: L={o['L']:.4e} log_g L={np.log(o['L'])/np.log(g):.3f} L/g^2={o['L']/g**2:.3f} tau*g={o['tau']*g:.3f} lightcoef 16Th/L={16*th/o['L']:.3g} nu={o['nu']:.3g}", flush=True)
