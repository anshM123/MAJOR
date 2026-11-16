"""u6: leading Regge action of Omega^U vs Omega (ADS4_FUNCTIONAL.md statement R), M02 and F_g3, on physical (m, beta) grids.

Leading Regge order of each moment: eps_k m^{2k} C_k -> 2 eps_k sigma^k P_{-1/2+i nu}(cosh beta) (+ forward W, P0) [D|H*]. Then
  FR^U = sum_k u^{k-2} 2 eps_k sigma^k int phi_k P  + 2a Lam^2 u^2 W(beta) + 5 a sigma u^3 P_{-1/2}(cosh beta)
  FR   = 2 int psi2 P + sigma u int psi3 P + 2a Lam^2 u^2 W + a sigma u^3 P0          (the old Omega; robust/ads.py)
The two differ at O(u p^2) (phi_k carry p^{2(k-2)}-pieces) and in the forward B5 coefficient (5a vs a): both are O(nu^2/m^2, 1/m^2)
effects, i.e. of the size of the corrections the leading action drops.  Printed: last negative m/Lam on the grid, plain and with the
C1 envelope (C1 u N, N = sum of |terms|).  [NE; leading order only; not a positivity statement about the full action]
Output: cache/u6.log
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from common import load, KS, EPS, SIGMA, Transforms

out = []
def pr(s):
    print(s, flush=True); out.append(s)

def betas(Lam):
    b = np.concatenate([np.logspace(np.log10(0.05), np.log10(400.0), 220), np.arange(6.0, 12.0, 0.05), np.arange(12.0, 40.0, 0.1)])
    bs = np.concatenate([b/Lam, np.linspace(0.02, 8.0, 120)])
    return np.unique(bs[bs <= 8.0])

KAP = np.unique(np.concatenate([np.logspace(0, 3, 600), np.arange(1.0, 12.0, 0.01), np.arange(12, 40, 0.05)]))
for name in ('M02', 'Fg3'):
    F, PH = load(name)
    TR = Transforms(PH)
    from rcommon import CosT
    c3 = CosT(F.psi3)
    for Dg in (20, 100, 1000):
        Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = KAP*Lam; u = 1/KAP**2
        last = {(v, h, C1): 0.0 for v in ('U', 'old') for h in SIGMA for C1 in (0, 1)}
        for b in betas(Lam):
            I, W, P0 = TR.regge(Lam, b)
            # psi3 transform for the old functional
            from rcommon import glang, nth, Einvsqrt
            th, wt = glang(nth(Lam*b)); E = Einvsqrt(b, np.sin(th))
            I3 = (2/np.pi)*np.sum(wt*c3(Lam*b*np.sin(th))*E)
            phys = (m*np.exp(-b/2) >= Dg - 2 - 1e-9)
            for h, s in SIGMA.items():
                termsU = [u**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS] + [2*F.a*Lam**2*u**2*W, 5*F.a*s*u**3*P0]
                termsO = [2*I[2]*np.ones_like(u), s*u*I3, 2*F.a*Lam**2*u**2*W, F.a*s*u**3*P0]
                for v, T in (('U', termsU), ('old', termsO)):
                    FR = sum(T); N = sum(np.abs(t) for t in T)
                    for C1 in (0, 1):
                        neg = np.where(phys & (FR - C1*u*N <= 0))[0]
                        if neg.size:
                            last[(v, h, C1)] = max(last[(v, h, C1)], KAP[neg.max()])
        pr(f"  {name} Dg={Dg}: last negative m/Lam (leading Regge; ++ / +-):  Omega^U: {last[('U','++',0)]:.2f} / {last[('U','+-',0)]:.2f}"
           f" (C1=1: {last[('U','++',1)]:.2f} / {last[('U','+-',1)]:.2f});   Omega: {last[('old','++',0)]:.2f} / {last[('old','+-',0)]:.2f}"
           f" (C1=1: {last[('old','++',1)]:.2f} / {last[('old','+-',1)]:.2f})")
open(os.path.join(HERE, 'cache', 'u6.log'), 'w').write("\n".join(out) + "\n")
