"""
f4: cost of the eta-inversion in route A, estimated on the near-gap action [D-sketch + NE].
CLS-D's graviton weights f_{k;2} have leading Regge kernels c_k eta^j r^{k-3} P_nu (j = 1,0,3,2,1,0 for k = 2..7; f1_bulkpoint.py).
With the d = 3 Legendre recursion  eta P_nu = 1/2 (1 + 1/(2 i nu)) P_{nu-i} + 1/2 (1 - 1/(2 i nu)) P_{nu+i}  (P_nu = P_{-1/2+i nu}(eta)),
the smeared functional  int dnu phi(nu) X_nu  (X: pure-harmonic leading kernel, the object the D-safe weights phi_k are defined for) equals
int dnu phit(nu) Y_nu (Y: CLS-D object, kernel eta^j P_nu) at the level of the leading Regge kernel iff (contour shifts by +-i)
   phi(nu) = 1/2 phit(nu+i)(1 + 1/(2i(nu+i))) + 1/2 phit(nu-i)(1 - 1/(2i(nu-i)))   (j = 1; iterate for j > 1),
whose large-nu solution is  phit = phi + j [phi''/2 - phi'/(2 nu) + phi/(2 nu^2)] + O(nu^-4)  (' = d/dnu).
At the bulk point Y has the B' density (CLS-D), so route A's near-gap action differs from the R = 1 (e7) action by
   dF = sum_k sum_i (phit_k - phi_k)/phi_k (p_i) I_k(p_i),   nu = Lam p.
The large-nu symbol is not valid at nu = O(1); the contribution of nodes with nu < 3 is reported separately (it needs the exact
difference-equation solution, plan item M3).
Usage: python3 f4_etainv.py [Dg ...] -> cache/f4_etainv.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import polynomial as Pl
from spf_common import NearGap, load_states
from common import KS

JPOW = {2: 1, 3: 0, 4: 3, 5: 2, 6: 1, 7: 0}
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def main(Dgs):
    pr("# f4: eta-inversion correction to the route-A near-gap action (large-nu symbol), D-safe NF02 weights")
    for Dg in Dgs:
        NG = NearGap(Dg); L = NG.Lam; pn = NG.pn; nu = L*pn
        ratio = {}
        for k in KS:
            c = np.array([float(v) for v in NG.F.PHx[k]])
            ph = Pl.polyval(pn, c); d1 = Pl.polyval(pn, Pl.polyder(c, 1))/L; d2 = Pl.polyval(pn, Pl.polyder(c, 2))/L**2
            dphi = JPOW[k]*(0.5*d2 - d1/(2*nu) + ph/(2*nu**2))
            ratio[k] = np.where(np.abs(ph) > 0, dphi/np.where(ph == 0, 1, ph), 0.0)
        lo = nu < 3.0
        worst = (0, None); worst_hi = (0, None); recs = []
        for h, DF, r0, ref in load_states(Dg):
            x, P = NG.pieces(h, r0)
            tot = {k: P[k]['dens'] + P[k]['curv'] + P[k]['rem'] for k in KS}
            dF = sum((ratio[k]*tot[k]).sum() for k in KS)
            dF_hi = sum((ratio[k]*tot[k])[~lo].sum() for k in KS)
            recs.append((h, r0['mfl'], r0['J'], DF, ref['Fex'], ref['Ns'], dF, dF_hi))
            if abs(dF)/ref['Ns'] > worst[0]:
                worst = (abs(dF)/ref['Ns'], (h, r0['mfl'], r0['J'], DF))
            if abs(dF_hi)/ref['Ns'] > worst_hi[0]:
                worst_hi = (abs(dF_hi)/ref['Ns'], (h, r0['mfl'], r0['J'], DF))
        grid = sorted(set(r[1] for r in recs))
        n = next((n0 for n0 in grid if all(r[4] + r[6] > 0 for r in recs if r[1] >= n0 - 1e-9)), None)
        mn = min(recs, key=lambda r: (r[4] + r[6])/r[5])
        pr(f"  Delta_g = {Dg} (Lam = {L:.2f}): max |dF|/N_s = {worst[0]:.2e} {worst[1]}; nodes nu >= 3 only: {worst_hi[0]:.2e} {worst_hi[1]}; "
           f"n with dF = {n}; min (F+dF)/N_s = {(mn[4] + mn[6])/mn[5]:+.4f} [{mn[0]} m/Lam={mn[1]} J={mn[2]}]")
        pr(f"     relative weight change |phit/phi - 1| at p = 0.1, 0.5, 0.9 (k = 2, 4, 6): " +
           "; ".join(", ".join(f"{abs(np.interp(pp, pn, ratio[k])):.1e}" for pp in (0.1, 0.5, 0.9)) for k in (2, 4, 6)))
    open(os.path.join(HERE, 'cache', 'f4_etainv.txt'), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
