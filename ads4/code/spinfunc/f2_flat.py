"""
f2: flat (bulk-point) limit of the D-safe functional with the graviton bulk-point densities of f1, on the flat heavy states.
For object X the flat heavy action is  F_X(m, J, h) = sum_k D_h(J) @ (u^{k-2} w phi_k(p) e_k(y, sigma) R_X(k, x)),  x = p/m_flat (u = 1/m^2),
NF02 weights, J <= 80, m/M in [1, 10] (same scan as objid/o1_flat.py); normalised by the stripped norm N_s = |B2| + |B3|.
X = 'A' is NF02 itself (B' identity, certified in NOFORWARD_FUNCTIONAL.md); 'PsiG' is the functional as currently defined for gravitons
(Psi_4, Psi_5, Psi_6 with their true graviton densities); 'PsiG2', 'PsiS' see spf_common.
Also: the per-k contributions of the k = 4, 5, 6 terms near x -> 1 that the object change moves.
Usage: python3 f2_flat.py -> cache/f2_flat.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from spf_common import R, OBJS
from common import gl01, Dh_all, ek, KS, SIGMA
from nf_common import NF, flat_parts_nf

F = NF()
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def action(h, u, Jmax, obj, npq=600):
    p, w = gl01(npq)
    s = SIGMA[h]; y = u*p*p; x = np.sqrt(y)
    Js, D = Dh_all(h, Jmax, y)
    Js2, B2, B3, partsA = flat_parts_nf(F, h, u, Jmax, p, w)
    cols = np.array([u**(k - 2)*w*F.phi(k, p)*ek(k, y, s)*R(obj, k, x) for k in KS]).T
    parts = D @ cols
    return Js, parts, np.abs(B2) + np.abs(B3), (B2 + B3) - partsA.sum(1)

if __name__ == "__main__":
    pr("# f2: flat heavy action of the D-safe functional (NF02 weights) with graviton bulk-point densities; min F/N_s over J <= 80")
    mgrid = [1.0, 1.001, 1.003, 1.01, 1.03, 1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0]
    worst = {}
    for h in ('++', '+-'):
        for mfl in mgrid:
            u = 1/mfl**2
            cells = []
            for obj in OBJS:
                Js, parts, Ns, idres = action(h, u, 80, obj)
                Fv = parts.sum(1)
                ok = (Js % 2 == 0) if h == '++' else (Js >= 4)
                rat = Fv[ok]/Ns[ok]; i = int(np.argmin(rat)); J0 = int(Js[ok][i])
                cells.append(f"{obj}: {rat[i]:+.4f} (J={J0})")
                if obj == 'A':
                    cells[-1] += f" [B' id {np.max(np.abs(idres)/Ns):.0e}]"
                if rat[i] < worst.get((h, obj), (9, 0, 0))[0]:
                    worst[(h, obj)] = (rat[i], mfl, J0)
            pr(f"  {h} m/M={mfl:6.3f}: " + " | ".join(cells))
    pr("# overall minima per object (helicity: min F/N_s at m/M, J)")
    for (h, obj), v in sorted(worst.items()):
        pr(f"   {h} {obj:6s}: {v[0]:+.4f} at m/M = {v[1]}, J = {v[2]}")
    pr("# what moves: per-k parts / N_s at the worst A-state and at m = M (J = 4, 5, 6), A vs PsiG")
    for h, J in (('++', 4), ('++', 6), ('+-', 4), ('+-', 5), ('+-', 6)):
        for mfl in (1.0, 1.05, 1.3):
            u = 1/mfl**2
            rows = {}
            for obj in ('A', 'PsiG'):
                Js, parts, Ns, _ = action(h, u, 80, obj)
                i = int(np.searchsorted(Js, J)); rows[obj] = (parts[i]/Ns[i], parts[i].sum()/Ns[i])
            pr(f"   {h} J={J} m/M={mfl}: A per-k " + " ".join(f"{v:+.3f}" for v in rows['A'][0]) + f" (sum {rows['A'][1]:+.4f});  "
               f"PsiG per-k " + " ".join(f"{v:+.3f}" for v in rows['PsiG'][0]) + f" (sum {rows['PsiG'][1]:+.4f})")
    open(os.path.join(HERE, 'cache', 'f2_flat.txt'), 'w').write("\n".join(out) + "\n")
