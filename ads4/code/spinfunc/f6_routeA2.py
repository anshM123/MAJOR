"""
f6 (M2/M3): near-gap action of the route-A' D-safe functional  Omega^{A'} = int psi_2 C^_2 + int psi_3 C^_3  (NF02 weights psi_2, psi_3; the
two UNSUBTRACTED CLS-D graviton sum rules f+_{2;2} = v'-u' (j' = 0) and f-_{3;2} = -1 (j' = 1)), whose bulk-point images are CLPS's
STRIPPED B^(1)_2, B^(1)_3 (g1_crystal4d.py + truncation argument, H0_SPINNING_FUNCTIONAL.md sec. 9), so that its flat limit is NF02
itself in the stripped basis.  Model [M], per state and node (64 GL nodes, nu = Lam p, x = p/(m/Lam), y = x^2, D~ = D_h (++), D_h/(1-y)^4 (+-)):
   flat   : B2 + B3 (stripped, 400-node, = NF02 flat action = e7's flat part)
   curv   : psi_2 * 2 (P_{-1/2+i nu}(cosh beta) - J0((2J+1)x))  +  psi_3 * u sigma (P - J0)          [2 eps_k sigma^k u^{k-2}, k = 2, 3]
   rem    : psi_2 * env_2 * rho_2  +  psi_3 * u sigma * env_3 * rho_3,  rho_k = exact scalar-skeleton remainder (e1; k = 3 by continuation [M]),
            env_k = B' envelope e_k (variant 'Bp') or the stripped one = e_k * (1-y)^{+2} (++), (1-y)^{-2} (+-) (variant 'strip').
Normalisation N_s = |B2 + curv_2| + |B3 + curv_3|.  Thresholds n and min F/N_s at Delta_g = 100, 300, 1000 (all e1 states, both Df).
M3 pieces: (i) eta-inversion only for k = 2 (j = 1): psi~_2 = psi_2 + (1/(2 Lam^2))(psi_2'' - psi_2'/p + psi_2/p^2) (large-nu symbol); its near-gap
shift; (ii) its change of k(0):  dk0 = -2 int dpsi_2 p^6 (NF normalised to k(0) = 1); A, B_AdS scale as 1/k(0) (graviton-pole term exact).
Usage: python3 f6_routeA2.py [Dg ...] -> cache/f6_routeA2.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import polynomial as Pl
from scipy.special import j0
from spf_common import NearGap, load_states
from common import gl01, Dh_one, ek, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk

out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def psi_coeffs(F, which):
    from fractions import Fraction as Fr
    return np.array([float(c) for c in (F.s2 if which == 2 else F.s3)])

def main(Dgs):
    pr("# f6: route-A' D-safe functional  int psi2 C^_2 + int psi3 C^_3  (stripped B2, B3 images), near-gap model [M]")
    for Dg in Dgs:
        NG = NearGap(Dg); F = NG.F; pn, pw, L = NG.pn, NG.pw, NG.Lam
        c2 = psi_coeffs(F, 2)
        ps2 = Pl.polyval(pn, c2); d1 = Pl.polyval(pn, Pl.polyder(c2, 1)); d2 = Pl.polyval(pn, Pl.polyder(c2, 2))
        dpsi2 = (d2 - d1/pn + ps2/pn**2)/(2*L**2)            # eta-inversion symbol, j = 1, in p-derivatives
        # dk0 with a fine quadrature
        pf, wf = gl01(400); psf = Pl.polyval(pf, c2)
        dpf = (Pl.polyval(pf, Pl.polyder(c2, 2)) - Pl.polyval(pf, Pl.polyder(c2, 1))/pf + psf/pf**2)/(2*L**2)
        dk0 = -2*np.sum(wf*dpf*pf**6)
        k0chk = -2*np.sum(wf*psf*pf**6) + 2*np.sum(wf*F.psi3(pf)*pf**4)
        recs = []
        for h, DF, r0, ref in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(L*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            cv2 = p2*2*(P - jz); cv3 = p3*u*s*(P - jz)
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            rem = {}
            for var, scale in (('Bp', 1.0), ('strip', fac)):
                rem[var] = (p2*envk(2, x)*scale*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*scale*r0['rho'][3]).sum()
            Ns = abs(B2 + cv2.sum()) + abs(B3 + cv3.sum())
            Dh = Dh_one(h, J, y); Dt = Dh if h == '++' else Dh/(1 - y)**4
            # eta-inversion shift of the k = 2 term (flat + curv + rem), symbol
            rat = np.where(ps2 != 0, dpsi2/np.where(ps2 == 0, 1, ps2), 0.0)
            t2 = p2*((2 - y)*Dt + 2*(P - jz) + envk(2, x)*fac*r0['rho'][2])
            dinv = (rat*t2).sum()
            recs.append(dict(h=h, Df=DF, mfl=mfl, J=J, Ns=Ns, base=B2 + B3 + cv2.sum() + cv3.sum(), rem=rem, dinv=dinv,
                             FexOld=ref['Fex'], NsOld=ref['Ns']))
        grid = sorted(set(r['mfl'] for r in recs))
        pr(f"  Delta_g = {Dg}: {len(recs)} records; NF k(0) check {k0chk:.12f}; eta-inversion dk0 = {dk0:+.3e} -> A, B_AdS x {1/(1 + dk0):.6f}, "
           f"bounds x {(1 + dk0)**-0.5:.6f}")
        for var in ('strip', 'Bp'):
            for inv in (False, True):
                Fv = lambda r: r['base'] + r['rem'][var] + (r['dinv'] if inv else 0.0)
                for delta in (0.0, 0.03):
                    n = next((n0 for n0 in grid if all(Fv(r) - delta*r['Ns'] > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
                    w = min(recs, key=lambda r: Fv(r)/r['Ns'])
                    bad = [r for r in recs if Fv(r) - delta*r['Ns'] <= 0]
                    pr(f"    rem-envelope {var:5s} inversion {str(inv):5s} delta={delta}: n = {n}; min F/N_s = {Fv(w)/w['Ns']:+.4f} "
                       f"[{w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}]; {len(bad)} non-positive"
                       + (f" (highest m/Lam {max(r['mfl'] for r in bad)})" if bad else ""))
        mr = max(abs(r['rem']['strip'])/r['Ns'] for r in recs); mi = max(abs(r['dinv'])/r['Ns'] for r in recs)
        pr(f"    max |remainder|/N_s (strip) = {mr:.4f}; max |inversion shift|/N_s = {mi:.2e}")
        old = min(recs, key=lambda r: r['FexOld']/r['NsOld'])
        pr(f"    (old B'-uplift Omega^NF,SC on the same states: min F/N_s = {old['FexOld']/old['NsOld']:+.4f} [{old['h']} m/Lam={old['mfl']} J={old['J']}])")
    open(os.path.join(HERE, 'cache', 'f6_routeA2_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
