"""
largespin/ls_core.py (LARGESPIN_ANALYTIC.md): shared helpers.

  * NF02 weights psi2, psi3 (spinfunc nf_common.NF, exact rational coefficients, normalised k(0) = 1), imported read-only;
  * Phi_k(x)  = int_0^1 psi_k(p) cos(p x) dp                     (Gauss-Legendre, exact for the polynomial part up to the node count);
  * OmLR_k(t) = int_0^1 psi_k(p) P_{-1/2 + i Lam p}(cosh t) dp     (the leading-Regge action at H2 impact parameter t, per unit weight)
              = (2/pi) int_0^{pi/2} Phi_k(Lam s(phi)) / cosh(s/2) dphi,   s(phi) = 2 asinh(sinh(t/2) sin phi)
    (Mehler-Dirichlet with the substitution of Lemma C; smooth integrand, no singularity);
  * hybrid state pieces from spinfunc/ee1000_common.hyb (read-only import).
Nothing in the project outside code/largespin is written.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
CACHE = os.path.join(HERE, 'cache')
import numpy as np
import ee1000_common as E

Q = E.Q
_F = {}


def NFw(Dg=1000):
    if Dg not in _F:
        _F[Dg] = Q.NearGap(Dg).F
    return _F[Dg]


def lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


_GL = {}


def gl01(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = ((x + 1)/2, w/2)
    return _GL[n]


def Phi(k, xs, Dg=1000, n=3000):
    """Phi_k(x) = int_0^1 psi_k cos(p x) dp for an array xs (x <= ~ 4000 is exact to round-off with n = 3000)."""
    F = NFw(Dg); p, w = gl01(n)
    ps = F.psi2(p) if k == 2 else F.psi3(p)
    xs = np.atleast_1d(np.asarray(xs, float))
    out = np.empty_like(xs)
    for i0 in range(0, len(xs), 256):
        xx = xs[i0:i0 + 256]
        out[i0:i0 + 256] = (np.cos(np.outer(xx, p))*(ps*w)).sum(1)
    return out


def OmLR(k, ts, Dg=1000, nphi=4000):
    """int_0^1 psi_k(p) P_{-1/2+i Lam p}(cosh t) dp for each t in ts."""
    L = lam(Dg)
    ph, wph = gl01(nphi); ph = ph*np.pi/2; wph = wph*np.pi/2
    ts = np.atleast_1d(np.asarray(ts, float))
    out = np.empty_like(ts)
    for i, t in enumerate(ts):
        s = 2*np.arcsinh(np.sinh(t/2)*np.sin(ph))
        out[i] = (2/np.pi)*(Phi(k, L*s, Dg)/np.cosh(s/2)*wph).sum()
    return out


def OmFlat(k, bs, Dg=1000, nth=4000):
    """int_0^1 psi_k(p) J0(b p) dp (flat Regge image) via J0(z) = (2/pi) int_0^{pi/2} cos(z sin th) dth."""
    th, w = gl01(nth); th = th*np.pi/2; w = w*np.pi/2
    bs = np.atleast_1d(np.asarray(bs, float))
    return np.array([(2/np.pi)*(Phi(k, b*np.sin(th), Dg)*w).sum() for b in bs])


def beta_of(J, mfl, Dg=1000):
    """H2 impact parameter of the state: sinh(beta/2) = (J + 1/2)/m exactly."""
    return 2*np.arcsinh((J + 0.5)/(mfl*lam(Dg)))


def mgap(Dg, J):
    return E.mgap(Dg, J)


def hyb(Dg, J, mfl, h, npn=512):
    return E.hyb(Dg, J, mfl, h, npn)


def hyb_fast(Dg, J, mfl, h):
    """hybrid pieces without mpmath: flat NF02 parts (nf_common.flat_parts_nf, exact Legendre/Jacobi recursion) and the exact
    leading-Regge curvature cv_k = w_k [OmLR_k(beta) - OmFlat_k(bhat)] (Mehler-Dirichlet, ls_core)."""
    L = lam(Dg); F = NFw(Dg); u = 1/mfl**2; s = Q.SIGMA[h]
    pf2, wf2 = Q.gl01(max(400, 2*J + 100))
    Js, B2, B3, _ = Q.flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
    beta = beta_of(J, mfl, Dg); bh = (2*J + 1)/mfl
    c2 = 2*(OmLR(2, [beta], Dg)[0] - OmFlat(2, [bh], Dg)[0])
    c3 = u*s*(OmLR(3, [beta], Dg)[0] - OmFlat(3, [bh], Dg)[0])
    Ns = abs(B2 + c2) + abs(B3 + c3)
    return dict(base=B2 + B3 + c2 + c3, Ns=Ns, B2=B2, B3=B3, cv2=c2, cv3=c3, beta=beta, bhat=bh, u=u)
