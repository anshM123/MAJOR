"""e(p): L(g, Theta=g^p) for the light-robust dual LP.  usage: escan.py g N K p1,p2 [Jlin_mult nJlog nt nlight]"""
import sys, time, numpy as np, lp2
g = int(sys.argv[1]); N = int(sys.argv[2]); K = int(sys.argv[3])
ps = [float(x) for x in sys.argv[4].split(',')]
jm = float(sys.argv[5]) if len(sys.argv) > 5 else 12; nJ = int(sys.argv[6]) if len(sys.argv) > 6 else 80
nt = int(sys.argv[7]) if len(sys.argv) > 7 else 50; nl = int(sys.argv[8]) if len(sys.argv) > 8 else 40
for p in ps:
    t0 = time.time(); th = g ** p
    o = lp2.solve(g, K, N, th, Jlin=int(jm * g), nJlog=nJ, nt=nt, nlight=nl)
    if o is None: print(g, N, K, p, 'fail', flush=True); continue
    print(f"g={g} N={N} K={K} p={p}: L={o['L']:.5e} e=log_g L={np.log(o['L'])/np.log(g):.4f} tau={o['tau']:.4e} tau*g={o['tau']*g:.3f} lightcoef={16*th/o['L']:.4g} nu={o['nu']:.3g} [{time.time()-t0:.0f}s]", flush=True)
