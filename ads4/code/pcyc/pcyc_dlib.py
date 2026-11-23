"""
pcyc_dlib (PCYC_CLOSE.md section 3): the scalar and EE (k = 4, J = 0) CLS-D integrands on the PHYSICAL chain D' in coordinates in which
D' is a product domain and the exact B^{Delta/2} is a pure exponential:

  x-block:  t_x in (0, inf), r_x in (0, 1), s in (-1, 1):   B_x = (1 - s_x + q_x)/(1 + s_x + q_x) = exp(-2 t_x/Delta),
            u_x = sigma r_x, v_x = sigma (1 - r_x), x_p = s sqrt(u_x v_x), q_x = u_x v_x (1 - s^2) = sigma^2 k_x,
            k_x = r_x(1 - r_x)(1 - s^2), sigma = 2T/(1 + sqrt(1 - 4 k_x T^2)), T = tanh(t_x/Delta)   (exact inverse of B_x),
  y-block:  t_y, r_y, phi in (-pi/2, pi/2):  same with k_y = r_y(1 - r_y) cos^2 phi, y_p = sin(phi) sqrt(u_y v_y),
  x0:       alpha in (-pi, pi], psi = alpha - phi, x~0 = cos(phi) cos(psi), y0 = sqrt(u_y v_y) sin(psi)   (both sheets n0 = +-1).
t -> inf is the light cone of x5 = e (the B zero, boundary of D'); t, r, s, phi, alpha -> ends are the atom divisors.  The map covers D'
exactly once (pcyc_region: D' = the CLS-D causal domain 4>5>3>0>1>2).
Integrand (raw measure; constants dropped):  |Phi_scal| x Jacobian x exp((i nu/2) log A)  [A > 0 on D'];  EE = (that) x Q (rd_core.Qratio,
real on D'; the overall factor a of EE = a Q scal is dropped: it cancels in every ratio used).
Imports rd_core (RD_CLOSE) read-only.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import rd_core as C
np.seterr(all='ignore')

LOG4 = np.log(4.0)


def sig_map(t, k, Delta):
    """sigma(t) and log(sigma dsigma/dt) for B = exp(-2t/Delta) at fixed k."""
    T = np.tanh(t/Delta)
    R = np.sqrt(1.0 - 4*k*T*T)
    sig = 2*T/(1 + R)
    dsdT = 2/(1 + R) + 8*k*T*T/(R*(1 + R)**2)
    lj = np.log(sig) + np.log(dsdT) + np.log1p(-T*T) - np.log(Delta)
    return sig, lj


def point(tx, rx, s, ty, ry, phi, al, Delta, want_Q=False):
    """arrays (same shape) -> dict(L = log|integrand| incl. all Jacobians (real), logA, lat (leading-Regge log A - log 4), Q)."""
    kx = rx*(1 - rx)*(1 - s*s)
    cf = np.cos(phi); sf = np.sin(phi)
    ky = ry*(1 - ry)*cf*cf
    sx, ljx = sig_map(tx, kx, Delta)
    sy, ljy = sig_map(ty, ky, Delta)
    ux, vx = sx*rx, sx*(1 - rx)
    uy, vy = sy*ry, sy*(1 - ry)
    rxv = np.sqrt(ux*vx); ryv = np.sqrt(uy*vy)
    psi = al - phi
    z = [cf*np.cos(psi), ux, vx, uy, vy, s*rxv, sf*ryv]
    z = [np.asarray(v, complex) for v in z]
    y0 = (ryv*np.sin(psi)).astype(complex)
    xsq = (-(ux*vx)*(1 - s)*(1 + s)).astype(complex)
    g = C.geometry(z, y0val=y0, xsqval=xsq)
    F = C.Flow(0.0, Delta, np.ones(7))
    vals = {k: v for k, v in g.items() if not isinstance(v, list)}
    L = np.real(F.logabs_parts(vals))
    # raw measure -> (u_x, v_x, s, u_y, v_y, phi, psi):  sqrt(ux vx) sqrt(uy vy) cos(phi)^2 |sin psi|;  (u, v) -> (t, r): sigma dsigma dr
    L = L + np.log(rxv) + np.log(ryv) + 2*np.log(np.abs(cf)) + np.log(np.abs(np.sin(psi))) + ljx + ljy
    logA = np.real(F.logA(vals))
    lat = -np.log(rx) + np.log1p(-rx) + np.log(ry) - np.log1p(-ry) - np.log1p(-s*s) + 2*np.log(np.abs(cf)) - 2*np.log(np.abs(np.sin(al)))
    out = dict(L=L, logA=logA, lat=lat, eps=logA - LOG4 - lat, g=g)
    if want_Q:
        out['Q'] = np.real(C.Qratio(g, Delta))
    return out
