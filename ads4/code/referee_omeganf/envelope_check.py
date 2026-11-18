"""(1) Where does the C1-envelope threshold of nf6 bind (bhat < 10 or >= 10)?  Restrict to bhat >= 10.
(2) Size of the Regge-envelope normalisation N = sum_k |T_k| vs the integrand norm Nabs = sum_k u^{k-2} 2 eps_k int |phi_k| |P_{-1/2+i Lam p}|
    (the natural scale of an O(1/m^2) correction density that does not share the leading action's cancellation)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np, mpmath as mp
from nf_common import NF, TransformsNF, KS, EPS, SIGMA
from common import glang, nth, Einvsqrt
from numpy.polynomial.legendre import leggauss
F = NF(); TR = TransformsNF(F)
def betas(Lam):
    b = np.concatenate([np.logspace(np.log10(0.05), np.log10(400.0), 220), np.arange(6.0, 12.0, 0.05), np.arange(12.0, 40.0, 0.1)])
    bs = np.concatenate([b/Lam, np.linspace(0.02, 8.0, 120)]); return np.unique(bs[bs <= 8.0])
KAP = np.unique(np.concatenate([np.logspace(0, 3, 600), np.arange(1.0, 12.0, 0.01), np.arange(12, 40, 0.05)]))
xg, wg = leggauss(2000); pn = (xg + 1)/2; pw = wg/2
def Pnu(nu, beta):
    th, wt = glang(nth(Lam*beta)); s = np.sin(th); E = Einvsqrt(beta, s)
    return (2/np.pi)*(np.cos(np.outer(nu, beta*s)) @ (wt*E))
for Dg in (20, 100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = KAP*Lam; u = 1/KAP**2
    res = {}
    for b in betas(Lam):
        I = TR.regge(Lam, b); phys = (m*np.exp(-b/2) >= Dg - 2 - 1e-9); bh = 2*Lam*np.sinh(b/2)
        for h, s in SIGMA.items():
            T = [u**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS]; FR = sum(T); N = sum(np.abs(t) for t in T)
            for C1 in (1, 3):
                for reg in ('all', 'bh<10', 'bh>=10'):
                    if (reg == 'bh<10' and bh >= 10) or (reg == 'bh>=10' and bh < 10):
                        continue
                    neg = np.where(phys & (FR - C1*u*N <= 0))[0]
                    if neg.size:
                        i = neg.max(); key = (h, C1, reg)
                        if KAP[i] > res.get(key, (0,))[0]:
                            res[key] = (KAP[i], bh, b, FR[i]/N[i])
    print(f"Dg={Dg}: last negative m/Lam with C1 u N [(m/Lam, bhat, beta, FR/N)]:")
    for k, v in sorted(res.items()):
        print("   ", k, tuple(round(float(x), 3) for x in v))
    # (2) integrand norm vs N at a few (m, bhat)
    for kap in (1.2, 1.5, 2.0):
        mm = kap*Lam; uu = 1/kap**2
        for bh in (5, 10, 15, 20, 40):
            b = 2*np.arcsinh(bh/(2*Lam))
            if mm*np.exp(-b/2) < Dg - 2:
                continue
            P = Pnu(Lam*pn, b); I = TR.regge(Lam, b)
            for h, s in SIGMA.items():
                T = [uu**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS]; FR = sum(T); N = sum(abs(t) for t in T)
                Tchk = [uu**(k - 2)*2*EPS[k]*s**k*np.sum(pw*F.phi(k, pn)*P) for k in KS]
                Nabs = sum(uu**(k - 2)*2*EPS[k]*np.sum(pw*np.abs(F.phi(k, pn)*P)) for k in KS)
                # coherent integrand: |sum_k ...| pointwise
                Ncoh = np.sum(pw*np.abs(sum(uu**(k - 2)*2*EPS[k]*s**k*F.phi(k, pn) for k in KS)*P))
                print(f"   Dg={Dg} m={kap}Lam bhat={bh} {h}: FR={FR:+.3e} N={N:.3e} FR/N={FR/N:+.3f} | check {abs(sum(Tchk)-FR)/N:.1e} | "
                      f"Nabs={Nabs:.3e} Ncoh={Ncoh:.3e}  FR/Ncoh={FR/Ncoh:+.2e}  N/Ncoh={N/Ncoh:.2e}")
# spot check P against mpmath
Lam = np.sqrt(98*99.)
print("P check:", float(Pnu(np.array([37.0]), 0.3)[0]), float(mp.re(mp.legenp(-0.5 + 37j, 0, mp.cosh(0.3), type=3))))
