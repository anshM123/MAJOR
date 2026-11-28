import sys, time, numpy as np
from combo2 import Combo
g = int(sys.argv[1]); p = None if sys.argv[2] == 'inf' else float(sys.argv[2]); K, N, Nf = map(int, sys.argv[3:6])
th = None if p is None else g ** p
res = {}
for uf in (False, True):
    t0 = time.time()
    cm = Combo(g, K, N, Nf, th, use_flip=uf, nJlog=40, nt=20, nl=20)
    cs = cm.cstar(1e-5, 1.0, it=14, iters_s=2)
    res[uf] = cs
    print(f"g={g} p={p} K={K} N={N} Nf={Nf} flip={uf}: tau*={4*cs:.4e} tau*g={4*cs*g:.3f} tau*g^3={4*cs*g**3:.1f} [{time.time()-t0:.0f}s]", flush=True)
