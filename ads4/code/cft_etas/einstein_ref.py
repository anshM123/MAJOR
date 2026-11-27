"""
Einstein (AdS4 black brane) reference for the zero-momentum shear correlator, at complex frequency.

Conventions as code/scope_sumrule/einstein_ads4_sumrule.py: L=1, r_+=1, z=1/r, f=1-z^3, T=3/(4 pi), 2 pi T = 3/2,
p=1, eps=2, s=4 pi, eta=1.  Mode eq. (f phi'/z^2)' + w^2 phi/(f z^2) = 0.  Near the boundary phi = A(1+w^2 z^2/2) + B z^3 + ...
G(w) := 3 B/A  (up to a w-independent contact constant).  Flux: Im G = rho = w |phi_h/A|^2.  Vacuum: G_0 = i w^3.
Ingoing: phi = f^{-i w/3} psi, psi regular at z=1.  Works for any w in the closed upper half plane (w = i y: Euclidean).

Outputs (when run):  (i) check Im G(w) = w*eta(w) on the real axis; (ii) Euclidean values dG(iy) := G(iy) - G_0(iy);
(iii) Poisson identity  [dG(0)-dG(iy)]/y = (1/pi) Int deta(w) y/(w^2+y^2) dw  at Matsubara heights y = 1.5 n.
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

z_, w_, u_ = sp.symbols('z w u')
_f = 1 - z_**3
_B1 = 2*sp.I*w_*z_**2 - 3*z_**2 - 2*_f/z_
_B0 = w_**2*(1 + z_)*(1 + z_**2)/(1 + z_ + z_**2)
_N = 6
_cs = sp.symbols('c1:%d' % (_N + 1))
_ser = 1 + sum(_cs[k]*u_**(k + 1) for k in range(_N))
_expr = _f.subs(z_, 1 - u_)*sp.diff(_ser, u_, 2) - _B1.subs(z_, 1 - u_)*sp.diff(_ser, u_) + _B0.subs(z_, 1 - u_)*_ser
_expr = sp.series(sp.together(_expr*(1 + (1 - u_) + (1 - u_)**2)), u_, 0, _N).removeO()
_sol = {}
for _k in range(_N):
    _eq = sp.expand(_expr.coeff(u_, _k)).subs(_sol)
    _sol[_cs[_k]] = sp.simplify(sp.solve(_eq, _cs[_k])[0])
_cfun = [sp.lambdify(w_, _sol[_cs[k]], 'numpy') for k in range(_N)]


def _rhs_psi(zz, y, wv):
    ff = 1 - zz**3
    b1 = 2j*wv*zz**2 - 3*zz**2 - 2*ff/zz
    b0 = wv**2*(1 + zz)*(1 + zz**2)/(1 + zz + zz**2)
    return [y[1], -(b1*y[1] + b0*y[0])/ff]


def _rhs_phi(zz, y, wv):
    ff = 1 - zz**3
    return [y[1], -((-3*zz**2 - 2*ff/zz)*y[1] + wv**2*y[0]/ff)/ff]


def _bdy_series(wv, z0, a0, a3, M=160):
    a = np.zeros(M + 1, dtype=complex)
    a[0] = a0; a[1] = 0.0; a[2] = wv**2*a0/2; a[3] = a3
    for m in range(3, M):
        s = wv**2*a[m - 1] + (m - 2)*(7 - 2*m)*a[m - 2] + ((m - 5)**2*a[m - 5] if m >= 5 else 0.0)
        a[m + 1] = -s/((m + 1)*(m - 2))
    n = np.arange(M + 1)
    return np.sum(a*z0**n), np.sum(n[1:]*a[1:]*z0**(n[1:] - 1))


def G(wv, h=1e-3, zm=0.5, rtol=1e-12):
    """3B/A for complex wv (Im wv >= 0); psi(1)=1 normalisation.  Returns (G, A)."""
    wv = complex(wv)
    c = [complex(cf(wv)) for cf in _cfun]
    p0 = 1 + sum(c[k]*h**(k + 1) for k in range(_N))
    dp = -sum((k + 1)*c[k]*h**k for k in range(_N))
    s = solve_ivp(_rhs_psi, (1 - h, zm), [p0 + 0j, dp + 0j], args=(wv,), method='DOP853', rtol=rtol, atol=1e-300)
    psi, dpsi = s.y[0, -1], s.y[1, -1]
    ffm = 1 - zm**3
    F = ffm**(-1j*wv/3)
    g = 1j*wv*zm**2/ffm
    phi, dphi = F*psi, F*(dpsi + g*psi)
    z0 = min(0.25, 1.5/max(abs(wv), 1e-9))
    us = []
    for (a0, a3) in ((1.0, 0.0), (0.0, 1.0)):
        v, d = _bdy_series(wv, z0, a0, a3)
        r = solve_ivp(_rhs_phi, (z0, zm), [v + 0j, d + 0j], args=(wv,), method='DOP853', rtol=rtol, atol=1e-300)
        us.append((r.y[0, -1], r.y[1, -1]))
    (u1, du1), (u2, du2) = us
    W = u1*du2 - du1*u2
    A = (phi*du2 - dphi*u2)/W
    B = (u1*dphi - du1*phi)/W
    return 3*B/A, A


def deta(wr):
    """thermal-minus-vacuum eta(w) = |phi_h/A|^2 - w^2 on the real axis (Einstein, p=1 units)."""
    _, A = G(wr)
    return abs(1/A)**2 - wr**2


def dG(wv):
    return G(wv)[0] - 1j*complex(wv)**3


if __name__ == '__main__':
    print('# Einstein reference, units 2 pi T = 1.5, p = 1, eta(0) = 1')
    for wr in (0.3, 1.0, 2.5, 6.0):
        g, A = G(wr)
        print('real w=%4.1f  Im G = %.12f   w*|1/A|^2 = %.12f' % (wr, g.imag, wr*abs(1/A)**2))
    G0 = dG(1e-7).real
    print('dG(0) (contact-free part) = %.3e' % G0)
    rows = []
    for n in range(1, 9):
        y = 1.5*n
        val = (G0 - dG(1j*y).real)/y
        # Poisson average of deta, using the OPE tail deta ~ -c w^-4 beyond w=40
        f = lambda wr: deta(wr)*y/(wr**2 + y**2)
        tot = 0.0
        for (a, b) in ((0, 2), (2, 5), (5, 10), (10, 20), (20, 40)):
            tot += quad(f, a, b, limit=200, epsabs=1e-13, epsrel=1e-11)[0]
        c4 = -deta(40.0)*40.0**4
        tail = quad(lambda wr: -c4*wr**-4*y/(wr**2 + y**2), 40, np.inf)[0]
        P = 2*(tot + tail)/np.pi
        rows.append((n, y, val, P))
        print('n=%d y=%5.2f  Euclid [dG(0)-dG(iy)]/y = %.10f   Poisson avg of deta = %.10f   diff = %.1e' % (n, y, val, P, val - P))
    print('large-y: dG(0)-dG(i y) -> 3/4 (=eps*S_E):', [(y, round(G0 - dG(1j*y).real, 8)) for y in (15, 30, 60)])
