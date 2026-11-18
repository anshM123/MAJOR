"""
o1: object identity (H14) at the flat / bulk-point level.  Bulk-point densities (large m at fixed x = nu/m, CHMRS saddle
r* = (1-s)/(1+s), s = sqrt(1-x^2), eta* = 1) of candidate objects, relative to the CHMRS dispersive C_{k,nu} (density E_k):
   C    : 1                                   (CHMRS eq. bulkpointbknuagain; what the skeleton model / B' dictionary uses for even k)
   Pi   : (1-r*^2)^k/(1+r*^2)                 (pure Regge moment Pi-hat_{k,nu}; CHMRS eq. extrafactors inverted)
   Psi  : (1-r*^2)^{k+2}/(1+r*^2) * S_k(r*)    (scalar bare light-ray commutator Phi_{-1,nu}[A G], A = r^{k+1} S_k: large-nu kernel
                                               4(1-r^2)D^{1/2} r^{-3} A P_nu from CHMRS magicblockidentity + largenukernelexpression;
                                               S_k = project shape A_k/r^{k-3}: 1 (even k, k = 3), (1+6r^2+r^4)/(1-r^2)^4 (k = 5), (1-r^2)^-4 (k = 7))
Flat heavy action of NF02's uplift with each object: parts_k = D_h @ (u^{k-2} w phi_k ek(k,y,sigma) R_X(k,x)); F/N_s with the stripped
norm N_s = |B2|+|B3| (B' identity: sum_k parts_C = B2 + B3).  Variant 'mix' = the actual Omega^{NF,SC} assignment: k = 2 -> C (C+_2 kept),
k = 3..7 -> Psi (Psi_k/N_k).  Odd k: the B' odd density times R_X (no scalar odd functional exists; R_X continued in k) [M].
Usage: python3 o1_flat.py   -> cache/o1_flat.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from common import gl01, Dh_all, ek, KS, SIGMA
from nf_common import NF, flat_parts_nf

def rstar(x):
    s = np.sqrt(np.maximum(1 - x*x, 0.0)); return (1 - s)/(1 + s)

def S_shape(k, r):
    if k == 5: return (1 + 6*r**2 + r**4)/(1 - r*r)**4
    if k == 7: return 1/(1 - r*r)**4
    return np.ones_like(r)

def R(obj, k, x):
    r = rstar(x)
    if obj == 'C': return np.ones_like(x)
    if obj == 'Pi': return (1 - r*r)**k/(1 + r*r)
    if obj == 'Psi': return (1 - r*r)**(k + 2)/(1 + r*r)*S_shape(k, r)
    if obj == 'PsiPure': return (1 - r*r)**(k + 2)/(1 + r*r)
    if obj == 'mix': return R('C', k, x) if k == 2 else R('Psi', k, x)
    raise ValueError(obj)

OBJS = ('C', 'Pi', 'Psi', 'PsiPure', 'mix')
F = NF()
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def action(h, u, Jmax, obj, npq=600):
    p, w = gl01(npq)
    s = SIGMA[h]; y = u*p*p; x = np.sqrt(y)
    Js, D = Dh_all(h, Jmax, y)
    Js2, B2, B3, partsC = flat_parts_nf(F, h, u, Jmax, p, w)
    cols = [u**(k - 2)*w*F.phi(k, p)*ek(k, y, s)*R(obj, k, x) for k in KS]
    parts = D @ np.array(cols).T
    return Js, parts.sum(1), np.abs(B2) + np.abs(B3), (B2 + B3) - partsC.sum(1)

if __name__ == "__main__":
    pr("# o1: flat (bulk-point) heavy action of the NF02 uplift with alternative object densities; min F/N_s over J <= 80 per m/M")
    pr("# density ratios R_X(k,x) at x = 0.3 / 0.6 / 0.9 / 0.99 / 1:")
    for obj in OBJS[1:4]:
        for k in KS:
            xs = np.array([0.3, 0.6, 0.9, 0.99, 1.0])
            pr(f"   {obj:7s} k={k}: " + " ".join(f"{v:.4g}" for v in R(obj, k, xs)))
    mgrid = [1.0, 1.001, 1.01, 1.03, 1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0]
    worst = {}
    for h in ('++', '+-'):
        for mfl in mgrid:
            u = 1/mfl**2
            cells = []
            for obj in OBJS:
                Js, Fv, Ns, idres = action(h, u, 80, obj)
                ok = (Js % 2 == 0) if h == '++' else (Js >= 4)
                rat = Fv[ok]/Ns[ok]; i = int(np.argmin(rat)); J0 = int(Js[ok][i])
                cells.append(f"{obj}: {rat[i]:+.4f} (J={J0})")
                if obj == 'C': cells[-1] += f" [B' id {np.max(np.abs(idres)/Ns):.0e}]"
                if rat[i] < worst.get((h, obj), (9, 0, 0))[0]: worst[(h, obj)] = (rat[i], mfl, J0)
            pr(f"  {h} m/M={mfl:6.3f}: " + " | ".join(cells))
    pr("# overall minima per object (helicity: min F/N_s at m/M, J)")
    for (h, obj), v in sorted(worst.items()):
        pr(f"   {h} {obj:7s}: {v[0]:+.4f} at m/M = {v[1]}, J = {v[2]}")
    open(os.path.join(HERE, 'cache', 'o1_flat.txt'), 'w').write("\n".join(out) + "\n")
