"""Last negative m/Lam of FR - C1 u N (nf6 definitions, bhat >= 10 only) for larger C1 (the effective values of REFEREE_UNIFORM_REMAINDER)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from nf_common import NF, TransformsNF, KS, EPS, SIGMA
F = NF(); TR = TransformsNF(F)
def betas(Lam):
    b = np.concatenate([np.logspace(np.log10(0.05), np.log10(400.0), 220), np.arange(6.0, 12.0, 0.05), np.arange(12.0, 40.0, 0.1)])
    bs = np.concatenate([b/Lam, np.linspace(0.02, 8.0, 120)]); return np.unique(bs[bs <= 8.0])
KAP = np.unique(np.concatenate([np.logspace(0, 3, 600), np.arange(1.0, 12.0, 0.01), np.arange(12, 40, 0.05)]))
C1S = (1, 3, 5)
import itertools
for (BLO, BHI), Dg in itertools.product(((10, 15), (15, 20), (20, 40)), (100, 1000)):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = KAP*Lam; u = 1/KAP**2
    last = {(h, c): 0.0 for h in SIGMA for c in C1S}
    for b in betas(Lam):
        bh = 2*Lam*np.sinh(b/2)
        if not (BLO <= bh <= BHI): continue
        I = TR.regge(Lam, b); phys = (m*np.exp(-b/2) >= Dg - 2 - 1e-9)
        for h, s in SIGMA.items():
            T = [u**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS]; FR = sum(T); N = sum(np.abs(t) for t in T)
            for c in C1S:
                neg = np.where(phys & (FR - c*u*N <= 0))[0]
                if neg.size: last[(h, c)] = max(last[(h, c)], KAP[neg.max()])
    print(f"bhat in [{BLO},{BHI}] Dg={Dg}: " + "  ".join(f"C1={c}: ++ {last[('++',c)]:.2f} / +- {last[('+-',c)]:.2f}" for c in C1S))
