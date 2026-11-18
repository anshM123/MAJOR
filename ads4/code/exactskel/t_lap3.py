"""Validation of skelexact3 (Arb, Laplace/Barnes form, even k) against skelexact (double, contour) at Delta_g = 100."""
import sys, time
sys.path.insert(0, '.')
import numpy as np
from skelexact import ExactState
from skelexact3 import Laplace3
Lam = np.sqrt(98*99)
nus = Lam*np.array([0.05, 0.4, 0.7, 0.9, 0.999])
for (D, J, Df) in ((100.05, 0, 3.0), (138.0, 2, 3.0), (110.0, 5, 2.0)):
    ref = ExactState(D, J, Df).values(nus, ks=(2, 4, 6))
    for kw in (dict(prec=300), dict(prec=500, N0=200, nq=240)):
        t = time.time(); L3 = Laplace3(D, J, Df, ks=(2, 4, 6), numax=nus.max(), **kw); t1 = time.time() - t
        out = []
        for k in (2, 4, 6):
            vals = [L3.value(k, nu) for nu in nus]
            out.append(f"k={k}: " + " ".join(f"{float(v.mid()):.10g}(rel {float(v.mid())/r-1:+.1e}, rad {float(v.rad())/abs(float(v.mid())):.0e})" for v, r in zip(vals, ref[k])))
        print(f"D={D} J={J} Df={Df} {kw} M={L3.M} build {t1:.0f}s total {time.time()-t:.0f}s")
        for o in out:
            print("   " + o)
