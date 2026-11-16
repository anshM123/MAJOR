"""v2: independent check of the CHMRS series (and of our normalisations q, a_Df) by the exact Mellin representation for J = 0
(mellin.py), Df = 3.  Prints m^{2k} C_{k,nu}[Delta,0] from both."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp, mellin, skel
lines = []
for Delta in (30.3,):
    m2, eta = skel.m2_eta(Delta, 0)
    for k in (2, 4):
        for nu in (0.0, 2.0):
            t = time.time()
            v = mellin.heavy_C(k, nu, Delta, 3.0, dps=25)
            s = skel.series_value('seriesC', k, 3.0, nu, Delta, 0)
            l = f"Delta={Delta} m={float(mp.sqrt(m2)):.3f} k={k} nu={nu}: Mellin {mp.nstr(mp.re(v*m2**k), 10)}  series(order 8) {s:.10f}  rel.diff {float(mp.re(v*m2**k))/s-1:+.2e}  ({time.time()-t:.0f} s)"
            print(l, flush=True); lines.append(l)
open(os.path.join(HERE, 'cache', 'out_v2.txt'), 'w').write("\n".join(lines) + "\n")
