"""
H13' (RESULT_D.md) for the D-safe functional (NF02 / Omega^{NF,SC}): common tools.

* Conventions are those of omegau/nf_common.py and exactskel/e7_nf.py (imported read-only):
    m = mfl*Lam, Lam = sqrt((Dg-2)(Dg-1)), bhat = (2J+1)/mfl = (2J+1) Lam/m,  cosh(beta) = eta = 1 + (2J+1)^2/(2 m2)
  (so 2 sinh(beta/2) ~ bhat/Lam), u = 1/mfl^2, x = nu/m = p sqrt(u).
  Hybrid action F_hyb = B2 + B3 + sum_k curv_k  (flat stripped NF02 action + leading AdS Regge curvature
  curv_k = u^{k-2} 2 eps_k sigma^k [int phi_k P_{-1/2+i Lam p}(cosh beta) - int phi_k J0(bhat p)]),
  N_s = |B2 + curv_2| + |B3 + curv_3| (as e7/nf_common.action_nf), exact action F = F_hyb + S, S = sum_k int u^{k-2} sigma^k phi_k E_k rho_k
  with rho_k the exact scalar-skeleton remainder (e8_states_arb.py definition).
* flat_rows(): B2, B3 for many J at once by the Legendre / Jacobi(0,8) recursion on one Gauss rule that is exact for every requested J
  (integrand is a polynomial of degree 2J + deg(psi) + 2 in p).
* exact_state(): states parametrised by an exactly representable Delta (multiple of 1/64) so that the exact Mack polynomial (mack_flint,
  FLINT fmpq_mpoly; identical to mack.mack_P_fast) has small rationals.  mfl is then derived from Delta (not rounded).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
import numpy as np
from nf_common import NF, TransformsNF, KS, EPS, SIGMA, state
from common import gl01

CACHE = os.path.join(HERE, 'cache')
F = NF()
TR = TransformsNF(F)


def Lam_of(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def flat_rows(h, u, Js):
    """{J: (B2, B3)} for the stripped NF02 flat action (normal form), h in '++', '+-'."""
    Js = sorted(Js); Jmax = Js[-1]
    p, w = gl01(Jmax + 64)
    y = u*p*p; x = 1 - 2*y
    s = SIGMA[h]
    W2 = w*F.psi2(p)*(2 - y); W3 = w*s*u*F.psi3(p)
    want = set(Js); out = {}
    if h == '++':
        P0, P1 = np.ones_like(x), x.copy()
        for J in range(0, Jmax + 1):
            if J == 0:
                Pc = P0
            elif J == 1:
                Pc = P1
            else:
                Pc = ((2*J - 1)*x*P1 - (J - 1)*P0)/J
                P0, P1 = P1, Pc
            if J in want:
                out[J] = (float(W2 @ Pc), float(W3 @ Pc))
    else:
        a, b = 0.0, 8.0
        Pm, Pn = np.ones_like(x), (a + 1) + (a + b + 2)*(x - 1)/2
        for J in range(4, Jmax + 1):
            n = J - 4
            if n == 0:
                Pc = Pm
            elif n == 1:
                Pc = Pn
            else:
                c = 2*n + a + b
                Pc = (((c - 1)*(a*a - b*b) + (c - 2)*(c - 1)*c*x)*Pn - 2*(n + a - 1)*(n + b - 1)*c*Pm)/(2*n*(n + a + b)*(c - 2))
                Pm, Pn = Pn, Pc
            if J in want:
                out[J] = (float(W2 @ Pc), float(W3 @ Pc))       # dtilde44 = P^(0,8)_{J-4}
    return out


def curvature(h, Lam, mfl, J, beta):
    s = SIGMA[h]; u = 1/mfl**2; bh = (2*J + 1)/mfl
    I = TR.regge(Lam, beta); f = TR.flat(bh)
    curv = np.array([u**(k - 2)*2*EPS[k]*s**k*(I[k] - f[k]) for k in KS])
    FR = sum(u**(k - 2)*2*EPS[k]*s**k*I[k] for k in KS)
    FJ0 = sum(u**(k - 2)*2*EPS[k]*s**k*f[k] for k in KS)
    return curv, FR, FJ0


def hybrid(h, Dg, mfl, J, B2, B3, beta=None):
    Lam = Lam_of(Dg)
    if beta is None:
        st = state(Dg, Lam, mfl, J)
        if st is None:
            return None
        beta = st['beta']
    curv, FR, FJ0 = curvature(h, Lam, mfl, J, beta)
    Fh = B2 + B3 + curv.sum()
    Ns = abs(B2 + curv[0]) + abs(B3 + curv[1])
    return dict(h=h, Dg=Dg, mfl=mfl, J=J, bh=(2*J + 1)/mfl, beta=beta, B2=B2, B3=B3, curv=curv, Fh=Fh, Ns=Ns,
                Fflat=B2 + B3, Nflat=abs(B2) + abs(B3), FR=FR, FJ0=FJ0)


def exact_state(Dg, mfl_target, J, den=64):
    """state with Delta = round(Delta_target*den)/den; returns dict like omegau.common.state plus mfl (derived)."""
    Lam = Lam_of(Dg)
    st0 = state(Dg, Lam, mfl_target, J)
    if st0 is None:
        return None
    Delta = round(st0['Delta']*den)/den
    alpha = Delta - J - 2.0
    if alpha < Dg - 2 - 1e-12:
        return None
    abar = alpha + 2*J + 1
    m2 = alpha*abar; m = np.sqrt(m2)
    eta = (alpha**2 + abar**2)/(2*alpha*abar)
    return dict(m=m, m2=m2, eta=eta, beta=float(np.arccosh(eta)), Delta=Delta, mfl=m/Lam, J=J)
