#!/usr/bin/env python3
"""C1_FINAL: one-point ("light-ray") matrices of E(n) = int du T(u,n) and E1(n) = int du u T(u,n) on the stress-tensor family,
computed from first principles in the two free 3d CFTs that span the parity-even <TTT> (free real boson, free Dirac fermion).

Library module (imported by c1f_onepoint.py and c1f_compose.py).  Conventions: metric (-,+,+), vectors are upper-index (t,x,y).
States |T(P),h> := eps_h^{mu nu}(P) T_{mu nu}(P)|0>, eps_h = e_h e_h, e_h(P) = Lambda(P) e_h(rest), e_pm(rest) = (0,1,+-i)/sqrt2,
Lambda(P) = pure boost (M,0,0) -> P.  Two-particle wave function Psi(k1,k2) = F(k1+k2) psi_h(k1,k2):
   boson:   psi_h = eps_h^{mn} t_{mn},   t = k1 k2 + k2 k1 - eta (k1.k2)              (symmetric pair, norm 1/2 int dPi dPi)
   fermion: psi_h = eps_h^{mn} ubar(k1) gamma_m v(k2) (k1-k2)_n                          (particle k1, antiparticle k2)
Light-ray operators at future null infinity, n = (cos phi, sin phi), nb = (1,n) (derived from the mode expansion at scri+; the
fermion kernel differs from the boson's only off the omega=omega' diagonal at second order, which int du u T does not see):
   E(phi)  = (1/8pi^2) int dw w a^dag a(w n)
   E1(phi) = (1/8pi^2) int dw w (-i/2)(a^dag d_w a - d_w a^dag a)
Matrix elements:  <Psi'|E|Psi>  = int d^3P F'* F H(P,n),
                  <Psi'|E1|Psi> = int d^3P [(-i/2)(F'* nb.dF - nb.dF'* F) H + F'* M F]
with H, M (2x2 in helicity h = +2, -2) given by a single phase-space point (k1 = w nb, k2 = P - k1 null => w = P^2/(2 P.nb)).
Normalised: Hh = H/n(P), Mh = M/n(P) with <P'h'|Ph> = n(P) delta.
"""
import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0])
def dot(a, b): return a @ ETA @ b
def lower(v): return ETA @ v

# --- gamma matrices (2+1, {g^m,g^n} = 2 eta^{mn}) and massless spinors ---
G0 = np.array([[0, 1], [-1, 0]], dtype=complex)      # i sigma_2
G1 = np.array([[0, 1], [1, 0]], dtype=complex)       # sigma_1
G2 = np.array([[1, 0], [0, -1]], dtype=complex)      # sigma_3
GAM = [G0, G1, G2]
def slash(k):  # gamma^m k_m
    kl = lower(k); return kl[0]*G0 + kl[1]*G1 + kl[2]*G2
def spinor(k):  # u(k) = sqrt(2w) (sin(th/2), cos(th/2)), solves slash(k) u = 0; v(k) := u(k); u(w n) = sqrt(w) * chi(n)
    w = k[0]; th = np.arctan2(k[2], k[1])
    return np.sqrt(2*w)*np.array([np.sin(th/2), np.cos(th/2)], dtype=complex)
def ubar(u): return u.conj() @ G0

def boost(P):
    M = np.sqrt(-dot(P, P)); g = P[0]/M; b = P[1:]/P[0]; b2 = b @ b
    L = np.eye(3); L[0, 0] = g; L[0, 1:] = g*b; L[1:, 0] = g*b
    if b2 > 0: L[1:, 1:] += (g - 1)*np.outer(b, b)/b2
    return L
ERest = {+2: np.array([0, 1, 1j])/np.sqrt(2), -2: np.array([0, 1, -1j])/np.sqrt(2)}
def evec(P, h): return boost(P) @ ERest[h]
def eps(P, h): e = evec(P, h); return np.outer(e, e)          # upper indices

def psiB(k1, k2, h):
    e = evec(k1 + k2, h)
    t = np.outer(lower(k1), lower(k2)) + np.outer(lower(k2), lower(k1)) - ETA*dot(k1, k2)   # lower indices
    return np.einsum('m,n,mn->', e, e, t)
def psiF(k1, k2, h):
    e = evec(k1 + k2, h); el = lower(e)
    ub = ubar(spinor(k1)); v = spinor(k2)
    eslash = sum(e[m]*(ETA[m, m]*GAM[m]) for m in range(3))    # e^m gamma_m
    return (ub @ eslash @ v) * (el @ (k1 - k2))
PSI = {'B': psiB, 'F': psiF}
HEL = (+2, -2)

def nbar(phi): return np.array([1.0, np.cos(phi), np.sin(phi)])

def _point(P, phi):
    nb = nbar(phi); w = dot(P, P)/(2*dot(P, nb)); k1 = w*nb; k2 = P - k1
    J = 1.0/((2*np.pi)**2*2*k2[0]*(1 - nb[1:] @ (k2[1:]/k2[0])))
    return w, k1, k2, J

def _HM_one(theory, P, phi, which, dw_rel=1e-5):
    """contribution of the particle that radiates into phi: which = 1 (particle 1 at phi) or 2 (particle 2 at phi)."""
    psi = PSI[theory]
    w, ka, kb, J = _point(P, phi)        # ka radiates, kb spectator
    nb = nbar(phi)
    def ps(wa, h):
        k = wa*nb
        return psi(k, kb, h) if which == 1 else psi(kb, k, h)
    dw = dw_rel*w
    v = np.array([ps(w, h) for h in HEL])
    dv = np.array([(ps(w + dw, h) - ps(w - dw, h))/(2*dw) for h in HEL])
    pref = w*J/(8*np.pi**2)
    H = pref*np.outer(v.conj(), v)
    M = pref*(-0.5j)*(np.outer(v.conj(), dv) - np.outer(dv.conj(), v))
    return H, M, v, dv, w

def norm_density(theory, P, nth=4000):
    M = np.sqrt(-dot(P, P)); L = boost(P); k = M/2
    th = (np.arange(nth) + 0.5)*2*np.pi/nth
    vals = []
    for t in th:
        k1 = L @ (k*np.array([1, np.cos(t), np.sin(t)])); k2 = L @ (k*np.array([1, -np.cos(t), -np.sin(t)]))
        vals.append(abs(PSI[theory](k1, k2, +2))**2)
    I = np.mean(vals)*2*np.pi/((2*np.pi)**4*4*M)
    return 0.5*I if theory == 'B' else I

_NC = {}
def HM(theory, P, phi):
    """normalised (Hh, Mh) 2x2 matrices, rows/cols ordered (h=+2, h=-2)."""
    key = (theory, round(float(np.sqrt(-dot(P, P))), 12))
    if key not in _NC: _NC[key] = norm_density(theory, P)
    n = _NC[key]
    H, M, *_ = _HM_one(theory, P, phi, 1)
    if theory == 'F':
        H2, M2, *_ = _HM_one(theory, P, phi, 2)
        H = H + H2; M = M + M2
    return H/n, M/n
