"""
f7 (M2 correction of f2): the bulk-point images of CLS-D-class functionals are the STRIPPED B^(1)_k (g1_crystal4d.py; H0 doc sec. 9),
so the flat limit of the B'-uplift  sum_k phi_k(p) (Lam^{2k}) X_k  with objects X_k of that class is
     F = sum_k D~_h @ (u^{k-2} w phi_k strip_k(y) R_k(x)),   strip_k = (2-y)/(1-y)^{(k-2)/2} (even), sigma/(1-y)^{(k-3)/2} (odd),
     D~ = D_h (++), D_h/(1-y)^4 (+-)   [stripped density, FLAT_GRAVITON_BOUND (2.33)/(2.35a)],
with R_k = 1 (CLS-D weights at every k) or R_k of f1 (the project's Psi_k, i.e. Omega^{NF,SC} as defined).  Neither is NF02: NF02 is
psi_2 B_2 + psi_3 B_3 in the stripped basis (= route A', object 'A2' here).  min F/N_s over J <= 80, m/M in [1, 10], N_s = |B2| + |B3|.
Usage: python3 f7_flat_stripped.py -> cache/f7_flat_stripped.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from spf_common import R_grav
from common import gl01, Dh_all, KS, SIGMA
from nf_common import NF, flat_parts_nf

F = NF()
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def strip(k, y, s):
    return (2 - y)/(1 - y)**((k - 2)/2) if k % 2 == 0 else s/(1 - y)**((k - 3)/2)

def action(h, u, Jmax, obj, npq=1200):
    p, w = gl01(npq); s = SIGMA[h]; y = u*p*p; x = np.sqrt(y)
    Js, D = Dh_all(h, Jmax, y)
    Dt = D if h == '++' else D/(1 - y)**4
    Js2, B2, B3, _ = flat_parts_nf(F, h, u, Jmax, p, w)
    Ns = np.abs(B2) + np.abs(B3)
    if obj == 'A2':
        return Js, B2 + B3, Ns
    R = (lambda k: R_grav(k, x)) if obj == 'PsiG_strip' else (lambda k: np.ones_like(x))
    cols = np.array([u**(k - 2)*w*F.phi(k, p)*strip(k, y, s)*R(k) for k in KS]).T
    return Js, (Dt @ cols).sum(1), Ns

if __name__ == "__main__":
    pr("# f7: flat limits with the stripped (CLS-D) dictionary; min F/N_s over J <= 80 (J even for ++, J >= 4 for +-)")
    objs = ('A2', 'CLS_Bp_uplift', 'PsiG_strip')
    pr("#   A2 = route A' (psi2 C2 + psi3 C3 -> NF02);  CLS_Bp_uplift = B'-uplift weights phi_k on CLS-D objects (R = 1);  "
       "PsiG_strip = Omega^{NF,SC} as defined (Psi_k, R_k of f1)")
    worst = {}
    for h in ('++', '+-'):
        for mfl in (1.0, 1.001, 1.01, 1.03, 1.05, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0, 10.0):
            cells = []
            for obj in objs:
                Js, Fv, Ns = action(h, 1/mfl**2, 80, obj)
                ok = (Js % 2 == 0) if h == '++' else (Js >= 4)
                rat = Fv[ok]/Ns[ok]; i = int(np.argmin(rat)); J0 = int(Js[ok][i])
                cells.append(f"{obj}: {rat[i]:+.4f} (J={J0})")
                if rat[i] < worst.get((h, obj), (9, 0, 0))[0]:
                    worst[(h, obj)] = (rat[i], mfl, J0)
            pr(f"  {h} m/M={mfl:6.3f}: " + " | ".join(cells))
    pr("# minima")
    for (h, obj), v in sorted(worst.items()):
        pr(f"   {h} {obj:14s}: {v[0]:+.4f} at m/M = {v[1]}, J = {v[2]}")
    open(os.path.join(HERE, 'cache', 'f7_flat_stripped.txt'), 'w').write("\n".join(out) + "\n")
