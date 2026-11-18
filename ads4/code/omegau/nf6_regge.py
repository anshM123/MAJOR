"""nf6: leading Regge action of Omega^{NF,SC} (NF02; same construction as u6_regge.py, no forward terms):
  FR = sum_k u^{k-2} 2 eps_k sigma^k int phi_k P_{-1/2+i Lam p}(cosh beta) dp   on physical (m, beta) grids; last negative m/Lam, plain and
  with the C1 envelope C1 u N (N = sum of |terms|).  Also the minimum of FR/N over states with bhat = 2 Lam sinh(beta/2) >= 10.
  [NE; leading order only].  Output: cache/nf6.log"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from nf_common import NF, TransformsNF, KS, EPS, SIGMA
def betas(Lam):      # as u6_regge.py (copied: that module runs on import)
    b = np.concatenate([np.logspace(np.log10(0.05), np.log10(400.0), 220), np.arange(6.0, 12.0, 0.05), np.arange(12.0, 40.0, 0.1)])
    bs = np.concatenate([b/Lam, np.linspace(0.02, 8.0, 120)])
    return np.unique(bs[bs <= 8.0])
KAP = np.unique(np.concatenate([np.logspace(0, 3, 600), np.arange(1.0, 12.0, 0.01), np.arange(12, 40, 0.05)]))
out = []
def pr(s):
    print(s, flush=True); out.append(s)
F = NF(); TR = TransformsNF(F)
for Dg in (20, 100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = KAP*Lam; u = 1/KAP**2
    last = {(h, C1): 0.0 for h in SIGMA for C1 in (0, 1, 3)}; mn = {h: (np.inf, None) for h in SIGMA}
    for b in betas(Lam):
        I = TR.regge(Lam, b)
        phys = (m*np.exp(-b/2) >= Dg - 2 - 1e-9); bh = 2*Lam*np.sinh(b/2)
        for h, s in SIGMA.items():
            T = [u**(k - 2)*2*EPS[k]*s**k*I[k]*np.ones_like(u) for k in KS]
            FR = sum(T); N = sum(np.abs(t) for t in T)
            for C1 in (0, 1, 3):
                neg = np.where(phys & (FR - C1*u*N <= 0))[0]
                if neg.size:
                    last[(h, C1)] = max(last[(h, C1)], KAP[neg.max()])
            if bh >= 10:
                r = np.where(phys, FR/N, np.inf); i = int(np.argmin(r))
                if r[i] < mn[h][0]:
                    mn[h] = (r[i], (round(KAP[i], 3), round(bh, 2)))
    pr(f"  NF02 Dg={Dg}: last negative m/Lam (leading Regge; ++ / +-): C1=0: {last[('++',0)]:.2f} / {last[('+-',0)]:.2f};  C1=1: {last[('++',1)]:.2f} / {last[('+-',1)]:.2f};"
       f"  C1=3: {last[('++',3)]:.2f} / {last[('+-',3)]:.2f};  min FR/N at bhat >= 10: ++ {mn['++'][0]:+.3g} at {mn['++'][1]}, +- {mn['+-'][0]:+.3g} at {mn['+-'][1]}")
open(os.path.join(HERE, 'cache', 'nf6.log'), 'w').write("\n".join(out) + "\n")
