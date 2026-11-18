"""
H13' for the route-A' D-safe functional  Omega^{A'} = int psi2 C^_2 + int psi3 C^_3  (H0_SPINNING_FUNCTIONAL.md s10; model of
spinfunc/f6_routeA2.py), at large impact parameter.  Common tools (imported read-only: h13_common, omegau/nf_common).
Model [M] per state (as f6):
   flat  : B2 + B3 (stripped NF02, = h13_common.flat_rows)
   curv  : 2 int psi2 [P_{-1/2+i Lam p}(cosh beta) - J0(bhat p)] + u sigma int psi3 [P - J0]   (exact Mehler-Dirichlet / J0 transforms)
   rem   : 2? no -- int psi2 E_2 s_h rho_2 + u sigma int psi3 E_3 s_h rho_3,  E_k = nf_common.envk (B' envelope), s_h = (1-y)^{+2} (++),
           (1-y)^{-2} (+-) for the stripped envelope (variant 'strip', the headline of f6) or 1 (variant 'Bp').
           (f6 uses psi2*env2*rho2 (weight 1, not 2): the B'-normalised remainder rho_k is relative to eps_k m^{2k} C_k ~ E_k, and the
            k = 2 flat weight 2 of the curvature is the e_2(0) = 2 of E_2 -- kept identical to f6 here.)
   N_s   = |B2 + curv2| + |B3 + curv3|.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from h13_common import F, flat_rows, Lam_of, CACHE
from nf_common import ExactCosT, SIGMA, glang, nth, Einvsqrt, state, envk

CT = {2: ExactCosT(F.s2), 3: ExactCosT(F.s3)}


def transforms(Lam, beta, bh):
    b = Lam*beta; th, wt = glang(nth(b)); s = np.sin(th); E = Einvsqrt(beta, s)
    I = {k: (2/np.pi)*np.sum(wt*CT[k](b*s)*E) for k in (2, 3)}
    th2, wt2 = glang(nth(bh))
    f = {k: (2/np.pi)*np.sum(wt2*CT[k](bh*np.sin(th2))) for k in (2, 3)}
    return I, f


def hybridA(h, Dg, mfl, J, B2, B3, beta=None):
    Lam = Lam_of(Dg)
    if beta is None:
        st = state(Dg, Lam, mfl, J)
        if st is None:
            return None
        beta = st['beta']
    s = SIGMA[h]; u = 1/mfl**2; bh = (2*J + 1)/mfl
    I, f = transforms(Lam, beta, bh)
    c2 = 2*(I[2] - f[2]); c3 = u*s*(I[3] - f[3])
    return dict(h=h, Dg=Dg, mfl=mfl, J=J, bh=bh, beta=beta, B2=B2, B3=B3, curv2=c2, curv3=c3, Fh=B2 + B3 + c2 + c3,
                Ns=abs(B2 + c2) + abs(B3 + c3), Fflat=B2 + B3, Nflat=abs(B2) + abs(B3), FR=2*I[2] + u*s*I[3], FJ0=2*f[2] + u*s*f[3])


def remA(h, mfl, rho, pn, pw, variant='strip'):
    s = SIGMA[h]; u = 1/mfl**2; x = pn*np.sqrt(u); y = x*x
    fac = ((1 - y)**2 if h == '++' else (1 - y)**-2) if variant == 'strip' else 1.0
    p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
    t2 = p2*envk(2, x)*fac*rho[2]; t3 = p3*u*s*envk(3, x)*fac*rho[3]
    return float(t2.sum() + t3.sum()), float(np.abs(t2).sum() + np.abs(t3).sum())
