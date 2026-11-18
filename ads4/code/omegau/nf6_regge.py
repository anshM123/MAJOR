"""nf6: leading Regge action of Omega^{NF,SC} (NF02; same construction as u6_regge.py, no forward terms), with three envelope norms
(revised after REFEREE_OMEGA_NF M1):
  FR = sum_k u^{k-2} 2 eps_k sigma^k int phi_k P_{-1/2+i Lam p}(cosh beta) dp  on physical (m, beta) grids.
  Envelope  FR - C1 u N_x > 0  with
    N_comp = sum_k |u^{k-2} 2 eps_k int phi_k P|                 (component norm; shares NF02's large-bhat cancellation)
    N_coh  = int |sum_k u^{k-2} 2 eps_k sigma^k phi_k| |P| dp    (pointwise, coherent across k)
    N_pt   = sum_k u^{k-2} 2 eps_k int |phi_k| |P| dp            (pointwise, per k)
  P at the p-nodes from the Mehler-Dirichlet representation (same quadrature as rcommon.ReggeAction).
Printed per Dg: last negative m/Lam for C1 = 1, 3, 10, 17 and each norm (all bhat; and bhat > 20 only); the tolerated C1_max(n) =
min FR/(u N_x) over physical points with m >= n Lam and bhat > 20; and the ratios FR/N_x at m = 1.5 Lam for bhat ~ 10, 20, 40.
[NE; leading order only].  Output: cache/nf6.log"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from nf_common import NF, TransformsNF, KS, EPS, SIGMA, glang, nth, Einvsqrt
from common import gl01
def betas(Lam):      # as u6_regge.py (copied: that module runs on import)
    b = np.concatenate([np.logspace(np.log10(0.05), np.log10(400.0), 220), np.arange(6.0, 12.0, 0.05), np.arange(12.0, 40.0, 0.1)])
    bs = np.concatenate([b/Lam, np.linspace(0.02, 8.0, 120)])
    return np.unique(bs[bs <= 8.0])
KAP = np.unique(np.concatenate([np.logspace(0, 3, 600), np.arange(1.0, 12.0, 0.01), np.arange(12, 40, 0.05)]))
out = []
def pr(s):
    print(s, flush=True); out.append(s)
F = NF(); TR = TransformsNF(F)
pn, pw = gl01(300)
PHI = np.array([F.phi(k, pn) for k in KS])                      # (6, nodes)
C1S = (1, 3, 10, 17); NORMS = ('comp', 'coh', 'pt'); NLIST = (1.5, 2.1, 3.0)
def Pnodes(Lam, b):
    th, wt = glang(nth(Lam*b)); s = np.sin(th); E = Einvsqrt(b, s)
    out = np.empty(len(pn))
    for i0 in range(0, len(pn), 50):
        out[i0:i0 + 50] = (2/np.pi)*(np.cos(np.outer(pn[i0:i0 + 50], Lam*b*s)) @ (wt*E))
    return out
for Dg in (20, 50, 100, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); m = KAP*Lam; u = 1/KAP**2
    last = {(nm, C1, big): 0.0 for nm in NORMS for C1 in C1S for big in (False, True)}
    c1max = {(nm, n): np.inf for nm in NORMS for n in NLIST}
    ratios = {}
    U = np.array([u**(k - 2)*2*EPS[k] for k in KS])            # (6, nu)
    for b in betas(Lam):
        I = TR.regge(Lam, b)
        phys = (m*np.exp(-b/2) >= Dg - 2 - 1e-9); bh = 2*Lam*np.sinh(b/2)
        if not phys.any():
            continue
        P = Pnodes(Lam, b); aP = np.abs(P)
        Apt = np.abs(PHI) @ (pw*aP)                             # (6,)
        Npt = Apt @ U
        for h, s in SIGMA.items():
            sg = np.array([s**k for k in KS])
            T = (U.T*np.array([sg[i]*I[k] for i, k in enumerate(KS)])).T
            FR = T.sum(axis=0); Ncomp = np.abs(T).sum(axis=0)
            Ncoh = np.abs((U*sg[:, None]).T @ PHI) @ (pw*aP)
            NN = {'comp': Ncomp, 'coh': Ncoh, 'pt': Npt}
            for nm in NORMS:
                for C1 in C1S:
                    neg = np.where(phys & (FR - C1*u*NN[nm] <= 0))[0]
                    if neg.size:
                        last[(nm, C1, False)] = max(last[(nm, C1, False)], KAP[neg.max()])
                        if bh > 20:
                            last[(nm, C1, True)] = max(last[(nm, C1, True)], KAP[neg.max()])
                if bh > 20:
                    for n in NLIST:
                        sel = phys & (KAP >= n)
                        if sel.any():
                            c1max[(nm, n)] = min(c1max[(nm, n)], np.min(FR[sel]/(u[sel]*NN[nm][sel])))
            if Dg == 100:
                i = int(np.argmin(np.abs(KAP - 1.5)))
                for bt in (10, 20, 40):
                    if abs(bh - bt) < 0.3 and phys[i]:
                        ratios[(h, bt)] = tuple(FR[i]/NN[nm][i] for nm in NORMS)
    for big in (False, True):
        pr(f"  NF02 Dg={Dg} {'bhat > 20 only' if big else 'all bhat     '}: last negative m/Lam of FR - C1 u N_x (max over ++/+-): " +
           "; ".join(f"{nm}: " + " ".join(f"C1={C1}:{last[(nm, C1, big)]:.2f}" for C1 in C1S) for nm in NORMS))
    pr(f"  NF02 Dg={Dg}: tolerated C1_max(n) = min FR/(u N_x) over bhat > 20, m >= n Lam: " +
       "; ".join(f"{nm}: " + " ".join(f"n={n}:{c1max[(nm, n)]:.3g}" for n in NLIST) for nm in NORMS))
    if ratios:
        pr("  Dg=100, m = 1.5 Lam: FR/N_comp, FR/N_coh, FR/N_pt: " + "; ".join(f"{h} bhat~{bt}: " + "/".join(f"{x:.3g}" for x in v) for (h, bt), v in sorted(ratios.items())))
open(os.path.join(HERE, 'cache', 'nf6.log'), 'w').write("\n".join(out) + "\n")
