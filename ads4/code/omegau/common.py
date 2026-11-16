"""Common tools for OMEGA_U.md: the corrected AdS uplift Omega^U (FLAT_UNSTRIPPED.md section 0 item 5) of M02 / M01 / F_g3.

Omega^U = int_0^1 dp sum_{k=2..7} phi_k(p) eps_k Lam^{2k} C^{(-)^k}_{k, Lam p}  - a Lam^10 d_{nu^2} C^+_{4,nu}|_0 + (5a/2) Lam^10 C^-_{5,0}
  eps_k = 1 (k even), 1/2 (k odd);  phi_2 = psi2, phi_3 = psi3 - 4p^2 psi2, phi_4 = 2p^4 psi2 - p^2 psi3, phi_5 = -3p^6 psi2 + 2p^4 psi3,
  phi_6 = (p^8 psi2 - p^6 psi3)/2, phi_7 = (-p^10 psi2 + p^8 psi3)/2.        (a = forward coefficient of the stripped functional)

Normal form (per |c^h|^2; u = Lam^2/m^2, x = nu/m = p sqrt(u), y = x^2; D_++ = P_J, D_+- = d^J_44 = (1-y)^4 P^(0,8)_{J-4}):
  part_k = u^{k-2} int phi_k(p) [eps_k m^{2k} C^h_k(Lam p)] dp,
  flat (B' dictionary [I]):  eps_k m^{2k} C^h_k -> e_k(y) D_h(1-2y),  e_k = (2-y)/(1-y)^{(k+2)/2} (even),  sigma_h/(1-y)^{(k+1)/2} (odd)
  forward flat: a u^3 (2JJ - c_h) + 5 a sigma u^3,  c_++ = 5, c_+- = 37   (= 2a JJ (++), a(2JJ-42) (+-): the stripped values).
Floating point throughout.
"""
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
ROBUST = os.path.join(HERE, '..', 'robust')
UNIFORM = os.path.join(HERE, '..', 'uniform')
sys.path.insert(0, ROBUST)
from rcommon import CosT, glang, nth, Einvsqrt, pmul, peval, gl01, SIGMA, flat_parts   # noqa: E402  (read-only imports)
from check import load as load_stripped                                                  # noqa: E402

KS = (2, 3, 4, 5, 6, 7)
EPS = {k: (1.0 if k % 2 == 0 else 0.5) for k in KS}
CH = {'++': 5.0, '+-': 37.0}

def p2n(n):
    """monomial coefficients of p^(2n)."""
    c = np.zeros(2*n + 1); c[2*n] = 1.0
    return c

def padd(*ps):
    n = max(len(p) for p in ps)
    out = np.zeros(n)
    for p in ps:
        out[:len(p)] += p
    return out

def phis(F):
    """{k: monomial coefficients of phi_k(p)} for a (normalised) stripped Functional F."""
    s2, s3 = F.psi2, F.psi3
    return {2: s2.copy(),
            3: padd(s3, -4*pmul(p2n(1), s2)),
            4: padd(2*pmul(p2n(2), s2), -pmul(p2n(1), s3)),
            5: padd(-3*pmul(p2n(3), s2), 2*pmul(p2n(2), s3)),
            6: padd(0.5*pmul(p2n(4), s2), -0.5*pmul(p2n(3), s3)),
            7: padd(-0.5*pmul(p2n(5), s2), 0.5*pmul(p2n(4), s3))}

def ek(k, y, sigma=1.0):
    """flat B' density eps_k m^{2k} C_k / D_h  (even: CHMRS bulk point E_k; odd: sigma (1-y)^{-(k+1)/2})."""
    if k % 2 == 0:
        return (2 - y)/(1 - y)**(1 + k/2)
    return sigma/(1 - y)**((k + 1)/2)

def Ek_scalar(k, y):
    """the scalar bulk-point envelope E_k (CHMRS (3.40)), continued to odd k."""
    return (2 - y)/(1 - y)**(1 + k/2)

def env(k, y):
    """envelope used for transferring the scalar relative remainder rho_k to the B' moment: eps_k * (B' envelope)."""
    return np.abs(ek(k, y))

def jacobi_all(nmax, a, b, x):
    P = np.zeros((nmax + 1,) + np.shape(x)); P[0] = 1.0
    if nmax >= 1:
        P[1] = (a + 1) + (a + b + 2)*(x - 1)/2
    for n in range(2, nmax + 1):
        c = 2*n + a + b
        P[n] = (((c - 1)*(a*a - b*b) + (c - 2)*(c - 1)*c*x)*P[n - 1] - 2*(n + a - 1)*(n + b - 1)*c*P[n - 2])/(2*n*(n + a + b)*(c - 2))
    return P

def Dh_all(h, Jmax, y):
    """rows J = J0..Jmax of D_h(1-2y) (++: P_J, J0 = 0; +-: full d^J_44, J0 = 4)."""
    x = 1 - 2*y
    if h == '++':
        return np.arange(0, Jmax + 1), jacobi_all(Jmax, 0, 0, x)
    return np.arange(4, Jmax + 1), (1 - y)**4*jacobi_all(Jmax - 4, 0, 8, x)

def Dh_one(h, J, y):
    Js, D = Dh_all(h, J, y)
    return D[-1]

def flat_unstripped_parts(F, PH, h, u, Jmax, p, w):
    """per-k flat parts (rows J, columns k = 2..7) and forward, on quadrature (p, w).  Returns Js, parts (nJ x 6), FW (nJ)."""
    s = SIGMA[h]
    y = u*p*p
    Js, D = Dh_all(h, Jmax, y)
    cols = []
    for k in KS:
        cols.append(u**(k - 2)*w*peval(PH[k], p)*ek(k, y, s))
    parts = D @ np.array(cols).T
    JJ = Js*(Js + 1.0)
    FW = F.a*u**3*(2*JJ - CH[h]) + 5*F.a*s*u**3
    return Js, parts, FW

# ---------------------------------------------------------------- AdS kinematics and Regge transforms
def state(Dg, Lam, mfl, J):
    """heavy state with m = mfl Lam, spin J: None if unphysical (tau < Dg); else dict(m, m2, eta, beta, Delta)."""
    m = mfl*Lam
    Bq = 2*J + 1.0
    alpha = (-Bq + np.sqrt(Bq*Bq + 4*m*m))/2
    if alpha < Dg - 2 - 1e-9:
        return None
    abar = alpha + 2*J + 1
    eta = (alpha**2 + abar**2)/(2*alpha*abar)
    return dict(m=m, m2=alpha*abar, eta=eta, beta=float(np.arccosh(eta)), Delta=alpha + J + 2)

class Transforms:
    """int phi_k P_{-1/2+i Lam p}(cosh beta) dp, int phi_k J_0(bhat p) dp for k = 2..7, plus W(beta), P_{-1/2}(cosh beta)."""
    def __init__(self, PH):
        self.ct = {k: CosT(PH[k]) for k in KS}
    def regge(self, Lam, beta):
        b = Lam*beta
        th, wt = glang(nth(b))
        s = np.sin(th)
        E = Einvsqrt(beta, s)
        I = {k: (2/np.pi)*np.sum(wt*self.ct[k](b*s)*E) for k in KS}
        W = (2/np.pi)*np.sum(wt*(beta*s)**2/2*E)
        P0 = (2/np.pi)*np.sum(wt*E)
        return I, W, P0
    def flat(self, bh):
        th, wt = glang(nth(bh))
        return {k: (2/np.pi)*np.sum(wt*self.ct[k](bh*np.sin(th))) for k in KS}

FUNS = ('M02', 'M01', 'Fg3')
def load(name):
    F = load_stripped(name)
    return F, phis(F)
