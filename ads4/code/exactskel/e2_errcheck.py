"""
e2: numerical error of the exact moments.  For a sample of e1 states the moments are recomputed with a different contour set
(other line positions), finer line step h, more direct n-terms N0 and more EM nodes; the difference in rho_k (relative to E_k) is
the error estimate.  Also prints a rounding estimate eps*canc*|ex|/E_k.
Usage: python3 e2_errcheck.py Dg Df   -> cache/e2_<Dg>_Df<Df>.txt
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import skelexact
from skelexact import ExactState

Dg = int(sys.argv[1]); DF = float(sys.argv[2])
db = pickle.load(open(os.path.join(HERE, 'cache', 'e1_%d_Df%g.pkl' % (Dg, DF)), 'rb'))
keys = sorted(db)
sample = [k for k in keys if k[0] in (1.0, 1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0)]
skelexact.JSET = (0, 1, 2, 3, 5, 7, 9, 12, 16, 20, 26, 33, 40, 48, 58, 69, 81, 95, 111, 130, 152, 178)
KS = (2, 3, 4, 5, 6, 7)
Ek = lambda k, x: (2 - x*x)/(1 - x*x)**(1 + k/2)
lines = []
worst = 0.0
for key in sample:
    r = db[key]
    es = ExactState(r['Delta'], r['J'], DF, h=0.07, N0=600, nq=320)
    v = es.values(r['nus'])
    d = max(np.max(np.abs(v[k] - r['ex'][k])/Ek(k, r['x'])) for k in KS)
    rnd = max(np.max(1.1e-16*r['canc'][k]*np.abs(r['ex'][k])/Ek(k, r['x'])) for k in KS)
    worst = max(worst, d)
    s = f"m/Lam={key[0]:.2f} J={key[1]:2d}: max_k,nodes |d rho| = {d:.2e}   rounding estimate {rnd:.2e}"
    print(s, flush=True); lines.append(s)
lines.append(f"# worst |d rho| over sample: {worst:.2e}")
print(lines[-1])
open(os.path.join(HERE, 'cache', 'e2_%d_Df%g.txt' % (Dg, DF)), 'w').write("\n".join(lines) + "\n")
