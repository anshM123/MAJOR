"""
m2 (M6): H4 tolerance for route A' in the units of the calibrated expansion.  The exact scalar remainder is O(1/m^2) on the load-bearing window
(m1_scalar_scaling.py), so a graviton-specific correction to the bulk-point density is modelled as
   dens_k,h(x) -> dens_k,h(x) (1 + c_{k,h} w(x)/m^2),   w = 1 on x > x0 (else 0),
with sign patterns  'all' (c same for k = 2, 3 and both helicities), 'hel' (c_{++} = -c_{+-}), 'k' (c_2 = -c_3).
For each Delta_g, the critical |c| at which some state first becomes non-positive (route-A' functional psi2 C2 + psi3 C3, stripped densities,
f6 model, no inversion), for x0 = 0, 0.3, 0.5, and both signs.  Compare: scalar c2 ~ 10-160 (m1), measured <TTTT> Regge-moment c2 ~ 100-300
at x <~ 0.1 (SPINNING_TTTT_EXT; mostly a common factor c*).
Usage: python3 m2_tolerance_c2.py [Dg ...] -> cache/m2_tolerance_c2_<Dgs>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
from spf_common import NearGap, load_states
from common import gl01, Dh_one, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
X0S = (0.0, 0.3, 0.5)
PATS = ('all', 'hel', 'k')
def main(Dgs):
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr("# m2: critical |c| (density x (1 + c/m^2) on x > x0) breaking n = 1 for route A'; (c_-, c_+) = most negative / positive allowed")
    for Dg in Dgs:
        NG = NearGap(Dg); F = NG.F; pn, pw, L = NG.pn, NG.pw, NG.Lam
        acc = {}
        for h, DF, r0, ref in load_states(Dg):
            mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]; m = r0['m']
            x = pn*np.sqrt(u); y = x*x
            P = remk.Pders_vec(L*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
            pf2, wf2 = gl01(max(400, 2*J + 100))
            Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            fac = (1 - y)**2 if h == '++' else (1 - y)**-2
            base = B2 + B3 + (p2*2*(P - jz)).sum() + (p3*u*s*(P - jz)).sum() \
                   + (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
            Dh = Dh_one(h, J, y); Dt = Dh if h == '++' else Dh/(1 - y)**4
            d2 = p2*(2 - y)*Dt; d3 = p3*u*s*Dt
            for x0 in X0S:
                msk = x > x0
                G2, G3 = d2[msk].sum()/m**2, d3[msk].sum()/m**2
                for pat in PATS:
                    sh = 1.0 if (pat != 'hel' or h == '++') else -1.0
                    G = sh*(G2 + (G3 if pat != 'k' else -G3))
                    lo, hi = acc.get((x0, pat), (-np.inf, np.inf))
                    if G > 0:
                        lo = max(lo, -base/G)
                    elif G < 0:
                        hi = min(hi, -base/G)
                    acc[(x0, pat)] = (lo, hi)
        pr(f"  Delta_g = {Dg} (Lam = {L:.1f}):")
        for pat in PATS:
            pr(f"    {pat:3s}: " + "  ".join(f"x>{x0}: ({acc[(x0, pat)][0]:+.3g}, {acc[(x0, pat)][1]:+.3g})" for x0 in X0S))
    open(os.path.join(HERE, 'cache', 'm2_tolerance_c2_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write("\n".join(out) + "\n")
if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
