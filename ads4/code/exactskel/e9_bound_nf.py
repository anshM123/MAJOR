"""e9: D-safe (NF02 / Omega^NF,SC) conditional bound at thresholds n, same formula as omegau/nf7_summary.py `bound`
(k(0) = 1 normalisation, NF02's psi2).  Usage: python3 e9_bound_nf.py Dg n [n ...]  -> printed and appended to cache/e9.txt"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from scipy.integrate import quad
from nf_common import NF
F = NF()
Dg = int(sys.argv[1]); Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
lines = []
for n in [float(a) for a in sys.argv[2:]]:
    L = Lam/n
    LA = quad(lambda t: float(F.psi2(t))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]
    X = n**8*16*LA
    s = f"NF02 Dg={Dg} n={n:.3f}: |ghat3|^2 Lam^8 <= {X:.4g};  sqrt(t4^2+aT^2) <= {240*np.sqrt(X)/Lam**4:.3g}"
    print(s); lines.append(s)
open(os.path.join(HERE, 'cache', 'e9.txt'), 'a').write("\n".join(lines) + "\n")
