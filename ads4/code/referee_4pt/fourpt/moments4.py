"""
Tensor Regge moments for two spinning Regge pairs (1,3) and (2,4), frame of CHMRS (C.2):
   (x_1, x_2, x_3, x_4) = (-y, x/x^2, y, -x/x^2),  yhat = (1,0,0),  xhat = (cosh xi, sinh xi, 0).
Polarisation circles (rule: e_par at a point a of H^2 points to the other point b; e_perp.v = det[v, a, e_par]):
   pair (1,3) at yhat:  z(th) = yhat + cos th e_par + sin th e_perp,         e_par = (0,1,0), e_perp = (0,0,1);
   pair (2,4) at xhat:  z'(th) = xhat + cos th e'_par + sin th e'_perp,     e'_par = (-sinh xi, -cosh xi, 0), e'_perp = (0,0,-1),
                        and the operators carry z = I(xhat) z'(th) = -xhat + cos th e'_par + sin th e'_perp.
Modes  G_q(r, eta),  q = (q1, q3, q2, q4):  G(z1(th1), z3(th3), z2(th2), z4(th4)) = sum_q G_q e^{i(q1 th1 + q3 th3 + q2 th2 + q4 th4)}.
nu-harmonic kernel (Fourier-Helgason; generalises e2e/moments.py, which is the case mu_R = 0):
   K_{mu_L, mu_R}(xi) = (1/2pi) int_0^{2pi} dth (cosh xi - sinh xi cos th)^{-1/2 + i nu} e^{i mu_L th} e^{i mu_R phi_R(th)},
   e^{i phi_R} = (sinh xi - cosh xi cos th - i sin th)/(cosh xi - sinh xi cos th)     (angle of n_l(xhat) in the frame at xhat),
   Pihat^q_{k,nu} = int [d eta] K_{q1+q3, q2+q4}(xi) int dr r^{k-2} G_q(r, eta)     ([d eta] = 2 d eta).
"""
import os, sys, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from moments import make_grid, frame_data, frame_positions, series_converged, series_X       # noqa: E402

def kernel_np(nu, pairs, xi, M=None):
    """K_{muL,muR}(xi) for a list of (muL, muR); xi array.  Periodic trapezoid rule (midpoint nodes)."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    out = {p: np.zeros(xi.size, dtype=complex) for p in pairs}
    for i in range(xi.size):
        Mi = M if M is not None else int(2**np.ceil(np.log2(400 + 60*np.exp(xi[i]) + 40*abs(nu)*(1 + xi[i]))))
        th = 2*np.pi*(np.arange(Mi) + 0.5)/Mi
        ch, sh = np.cosh(xi[i]), np.sinh(xi[i])
        den = ch - sh*np.cos(th)
        base = np.exp((-0.5 + 1j*nu)*np.log(den))
        eR = (sh - ch*np.cos(th) - 1j*np.sin(th))/den
        for (a, b) in pairs:
            out[(a, b)][i] = np.mean(base*np.exp(1j*a*th)*eR**b)
    return out

def kernel_mp(nu, muL, muR, xi):
    xi = mp.mpf(xi); ch, sh = mp.cosh(xi), mp.sinh(xi)
    def f(th):
        den = ch - sh*mp.cos(th)
        eR = (sh - ch*mp.cos(th) - 1j*mp.sin(th))/den
        return den**(-mp.mpf(1)/2 + 1j*nu)*mp.exp(1j*muL*th)*eR**muR
    return mp.quad(f, mp.linspace(0, 2*mp.pi, 17))/(2*mp.pi)

def pol_grid4(xi, l):
    """polarisations for the four operators on a grid of nth = 2l+1 angles each: arrays [npts, nth,nth,nth,nth, 3] (z1, z3, z2, z4 order of axes: th1, th3, th2, th4)."""
    nth = 2*l + 1
    th = 2*np.pi*np.arange(nth)/nth
    npts = xi.size
    ch, sh = np.cosh(xi), np.sinh(xi)
    zy = np.stack([np.ones(nth), np.cos(th), np.sin(th)], axis=-1)                    # [nth, 3] (point independent)
    # z = -xhat + cos th e'_par + sin th e'_perp
    zx = np.zeros((npts, nth, 3))
    zx[:, :, 0] = -ch[:, None] - np.cos(th)[None, :]*sh[:, None]
    zx[:, :, 1] = -sh[:, None] - np.cos(th)[None, :]*ch[:, None]
    zx[:, :, 2] = -np.sin(th)[None, :]
    shp = (npts, nth, nth, nth, nth, 3)
    z1 = np.broadcast_to(zy[None, :, None, None, None, :], shp).copy()
    z3 = np.broadcast_to(zy[None, None, :, None, None, :], shp).copy()
    z2 = np.broadcast_to(zx[:, None, None, :, None, :], shp).copy()
    z4 = np.broadcast_to(zx[:, None, None, None, :, :], shp).copy()
    return z1, z3, z2, z4

def fourier4(M, l):
    """M[npts, nth^4] -> dict q -> [npts]  (q = (q1, q3, q2, q4))"""
    nth = M.shape[1]
    th = 2*np.pi*np.arange(nth)/nth
    F = np.fft.fftn(M, axes=(1, 2, 3, 4))/nth**4       # F[k] = mean M e^{-2 pi i k n/nth}  -> coefficient of e^{+i k th}
    out = {}
    for q in np.ndindex(*(2*l + 1,)*4):
        qq = tuple(int(v) - l for v in q)
        idx = (slice(None),) + tuple(v % nth for v in qq)
        out[qq] = F[idx]
    return out
