"""cert_t1: test of cert_moments against skelexact3 / e1 values (Delta_g = 100): ball contains the earlier value, radius, timing."""
import sys, time, pickle
sys.path.insert(0, '.')
import numpy as np
from cert_moments import CertMoments
Lam = float(np.sqrt(98*99))
db = pickle.load(open('cache/e1_100_Df3.pkl', 'rb'))
for key in [(1.04, 4), (1.3, 2)]:
    r = db[key]
    t = time.time()
    cm = CertMoments(r['Delta'], key[1], 3.0, Lam, ks=(2, 4, 6), verbose=True)
    t1 = time.time() - t
    for k in (2, 6):
        idx = [5, 30, 50, 63]
        vals = cm.values(k, r['nus'][idx])
        for i, v in zip(idx, vals):
            ref = r['ex'][k][i]
            print(f"{key} k={k} nu={r['nus'][i]:.2f}: ball {v.mid()} +- {float(v.rad()):.2e}; e1 {ref:.10g}; contains: {v.contains(arb_(ref)) if False else abs(float(v.mid()) - ref) <= float(v.rad()) + 1e-9*abs(ref)}")
    print(f"  build {t1:.0f}s total {time.time()-t:.0f}s")
