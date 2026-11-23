r"""
ee_quad (EESTEP_CLOSE.md §1.1): the EXACT per-nu EE kernel of C'_EE / C4^ds in the H7 §18.1 Regge frame, by direct 2-d quadrature.

Frame and scalar part: as lrc/lrc1_identify.py (logic re-implemented here; not imported):
  yh = (1,0,0), xh = (cosh B, sinh B, 0), |x| = 1, x2 = -xh, x4 = xh, x1 = -r yh, x3 = r yh, x0 = r n_th, n_th = (0, cos th, sin th),
  z0 = (1, cos ph, sin ph), ph = th + u;   metric (-,+,+).
  scalar integrand  S_nu = A^{i nu/2} |Nn|^{-1/2} m24^{3/4} (m20 m40)^{-1/4} c13^{5/4} (c10 c30)^{-1/4}   (CLS-D; repro_engine/integrand.py).
EE dressing (repro_engine/integrand.py header; h7i14/v3_regge_kernel.py header; lrc/lrc3_ee.py):
  * (13) side: the light-ray polarisations of T1, T3 are x10, x30, so the kernel's (13) tensor is x10^(x)2 (x) x30^(x)2; contracted with
    dual null polarisations w1, w3:  tau13 = (w1.x10)^2 (w3.x30)^2  (= r^4 x a trigonometric polynomial of degree 4 in th);
  * (24) side: EE (x) V042^{-4} = D(z2, z4) = [(H24 - 2 V_{2,04}(z2) V_{4,02}(z4))/(X20 X40)]^2, z0- and nu-independent (lrc3 (a)).
  K^EE_nu[w1, w3, z2, z4](r, B) = int dth du  S_nu * P(th),   P = tau13 * D.
Polarisations: w1, w3 = (1, cos a, sin a) (rest frame of yh, the (13) pair); z2, z4 = Lambda_B (1, cos a, sin a) (rest frame of xh, the (24)
pair), Lambda_B the boost taking (1,0,0) to xh.

Exact reduction (EESTEP_CLOSE.md §1.1).  On the circle x0 = r n_th, X20 = -(1-r^2)(1 - s c), X40 = -(1-r^2)(1 + s c), c = cos th,
s = 2 r sinh B/(1-r^2).  P is a polynomial of degree <= 8 in (cos th, sin th) divided by (X20 X40)^4; its sin-odd part integrates to zero
(S_nu is even under (th, ph) -> (-th, -ph)); with t = X20/X40, c = (1-t)/((1+t) s), (X20 X40)^{-4} = (1+t)^8/(256 t^4 (1-r^2)^8):
  P_even = sum_{b=-4}^{4} c_b(r, B; pols) t^b   exactly,   so   K^EE_nu = sum_b c_b Psi^(b)_nu,   Psi^(b)_nu := int S_nu t^b.
"""
import numpy as np
from numpy.polynomial import chebyshev as Cheb
from numpy.polynomial import polynomial as Poly


def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


def grid_u(ds=0.02, smax=60.0):
    s = np.arange(-smax, smax + ds/2, ds)
    u = np.pi/(1 + np.exp(-s))
    w = u*(1 - u/np.pi)*ds
    return u, w


def kernel_parts(r, beta, nth=128, ds=0.02):
    """log of the nu-independent scalar factor, log A on the (theta, u) grid, quadrature weights, theta array."""
    th = 2*np.pi*np.arange(nth)/nth
    u1, w1 = grid_u(ds)
    u = np.concatenate([u1, u1 + np.pi]); wu = np.concatenate([w1, w1])
    TH = th[:, None]; PH = TH + u[None, :]
    xh = np.array([np.cosh(beta), np.sinh(beta), 0.0]); yh = np.array([1.0, 0.0, 0.0])
    x2 = -xh; x4 = xh; x1 = -r*yh; x3 = r*yh
    n = np.stack([0*TH + 0*PH, np.cos(TH) + 0*PH, np.sin(TH) + 0*PH], -1)
    x0 = r*n
    z0 = np.stack([np.ones_like(PH), np.cos(PH), np.sin(PH)], -1)
    x10 = x1 - x0; x30 = x3 - x0; x20 = x2 - x0; x40 = x4 - x0
    m20 = -mdot(x20, x20); m40 = -mdot(x40, x40); m24 = -mdot(x2 - x4, x2 - x4)
    Nn = 2*(mdot(z0, x40)*(-m20) - mdot(z0, x20)*(-m40))
    c13 = 2*mdot(x10, x30)
    U = 0*TH + u[None, :]
    c10 = 4*r*np.sin(U/2)**2; c30 = 4*r*np.cos(U/2)**2
    A = Nn**2*c13/(m24*m40*m20*c10*c30)
    ok = bool(np.all(m20 > 0) and np.all(m40 > 0) and m24 > 0 and np.all(c13 > 0) and (np.all(Nn > 0) or np.all(Nn < 0)) and np.all(A > 0))
    logF = (-0.5*np.log(np.abs(Nn)) + 0.75*np.log(m24) - 0.25*np.log(m20*m40) + 1.25*np.log(c13) - 0.25*np.log(c10*c30))
    W = (2*np.pi/nth)*wu[None, :]
    return dict(logF=logF, logA=np.log(A), W=W, th=th, r=r, beta=beta, ok=ok)


def boost(beta):
    ch, sh = np.cosh(beta), np.sinh(beta)
    return np.array([[ch, sh, 0.0], [sh, ch, 0.0], [0.0, 0.0, 1.0]])


def pols(beta, a1, a3, a2, a4):
    nv = lambda a: np.array([1.0, np.cos(a), np.sin(a)])
    L = boost(beta)
    return nv(a1), nv(a3), L @ nv(a2), L @ nv(a4)


def dressing(th, r, beta, pol):
    """P(th) = tau13 D / r^4 on the x0-circle (arrays over th)."""
    w1, w3, z2, z4 = pol
    xh = np.array([np.cosh(beta), np.sinh(beta), 0.0]); yh = np.array([1.0, 0.0, 0.0])
    x2 = -xh; x4 = xh; x1 = -r*yh; x3 = r*yh
    n = np.stack([0*th, np.cos(th), np.sin(th)], -1)
    x0 = r*n
    X = lambda p, q: mdot(p - q, p - q)
    X20 = X(x2, x0); X40 = X(x4, x0); X24 = X(x2, x4)
    H24 = mdot(z2, z4)*X24 + 2*mdot(z2, x4 - x2)*mdot(z4, x2 - x4)
    V204 = (mdot(z2, x0 - x2)*X24 - mdot(z2, x4 - x2)*X20)/X40
    V402 = (mdot(z4, x0 - x4)*X24 - mdot(z4, x2 - x4)*X40)/X20
    D = ((H24 - 2*V204*V402)/(X20*X40))**2
    tau13 = (mdot(w1, x1 - x0)*mdot(w3, x3 - x0))**2/r**4
    return tau13*D


def t_of_th(th, r, beta):
    s = 2*r*np.sinh(beta)/(1 - r*r)
    c = np.cos(th)
    return (1 - s*c)/(1 + s*c), s


def laurent_cb(r, beta, pol, ncheb=33):
    """c_b, b = -4..4, with P_even = sum_b c_b t^b (P normalised by r^4 as in dressing()).  Also returns the size of the Chebyshev
    coefficients of Q(c) = P_even (X20 X40)^4 beyond degree 8 (must vanish) and the max |sin-odd part| (diagnostic only)."""
    k = np.arange(ncheb)
    cn = np.cos(np.pi*(k + 0.5)/ncheb)
    th = np.arccos(cn)
    s = 2*r*np.sinh(beta)/(1 - r*r)
    Pe = 0.5*(dressing(th, r, beta, pol) + dressing(-th, r, beta, pol))
    Po = 0.5*(dressing(th, r, beta, pol) - dressing(-th, r, beta, pol))
    Q = Pe*((1 - r*r)**2*(1 - s*s*cn*cn))**4
    ch = Cheb.chebfit(cn, Q, ncheb - 1)
    tail = np.max(np.abs(ch[9:]))/np.max(np.abs(ch[:9]))
    q = Cheb.cheb2poly(ch[:9])                      # Q(c) = sum_k q_k c^k, k = 0..8
    # P_even = sum_k q_k c^k / ((1-r^2)^8 (1 - s^2 c^2)^4),   c = (1-t)/((1+t) s)
    #        = sum_k q_k s^-k (1-t)^k (1+t)^(8-k) / (256 t^4 (1-r^2)^8)
    lau = np.zeros(9, dtype=complex if np.iscomplexobj(q) or np.iscomplexobj(s) else float)
    for kk in range(9):
        pk = Poly.polymul(Poly.polypow([1.0, -1.0], kk), Poly.polypow([1.0, 1.0], 8 - kk))   # coefficients of t^0..t^8
        lau += q[kk]*s**(-kk)*pk
    lau /= 256*(1 - r*r)**8
    cb = {b: lau[b + 4] for b in range(-4, 5)}       # t^j / t^4 = t^(j-4)
    return cb, tail, np.max(np.abs(Po))/np.max(np.abs(Pe)), q


def ukernels(r, beta, pol, ncheb=33):
    """The same reduction in the 'u-basis' that is uniform as rho -> 0:  P_even = (1-r^2)^-8/256 * sum_k q_k s^-k I_k(t),
    I_k = (1-t)^k (1+t)^(8-k) t^-4; returns the q_k s^-k and the Laurent coefficients of each I_k (for the block assembly)."""
    cb, tail, odd, q = laurent_cb(r, beta, pol, ncheb)
    s = 2*r*np.sinh(beta)/(1 - r*r)
    Ik = []
    for kk in range(9):
        pk = Poly.polymul(Poly.polypow([1.0, -1.0], kk), Poly.polypow([1.0, 1.0], 8 - kk))
        Ik.append({b: pk[b + 4] for b in range(-4, 5)})
    return [q[kk]*s**(-kk)/(256*(1 - r*r)**8) for kk in range(9)], Ik


def integrate(parts, nus, insert=None):
    """sum_grid W exp(logF + i nu logA/2) * insert(th)  for each nu (insert: array over th, broadcast along u)."""
    logF, logA, W = parts['logF'], parts['logA'], parts['W']
    base = W*np.exp(logF)
    if insert is not None:
        base = base*insert[:, None]
    return np.array([np.sum(base*np.exp(0.5j*nu*logA)) for nu in nus])
