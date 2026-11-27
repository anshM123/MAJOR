"""
Independent numerical check of the Einstein (t4 = 0) value of the CDP shear sum rule in d = 3:

    S := (1/(pi eps)) Int_{-inf}^{inf} dw [rho_T(w) - rho_0(w)]/w  =  3/8 .

AdS4 planar black brane, L = 1, r_+ = 1:  ds^2 = r^2(-f dt^2 + dx^2 + dy^2) + dr^2/(r^2 f),  f = 1 - z^3, z = 1/r,
T = 3/(4 pi).  The shear mode h^x_y = phi(z) e^{-i w t} obeys
    f phi'' + (f' - 2 f/z) phi' + w^2 phi / f = 0.
Normalisation (overall 1/(2 kappa^2) stripped; it cancels against P = r_+^3/(2 kappa^2)):
    rho_T(w) = w |phi_h|^2/|phi_bdy|^2  (conserved flux; phi ~ phi_h (1-z)^{-i w/3} at the horizon),
    rho_0(w) = w^3                       (pure AdS4: phi = (1 - i w z) e^{i w z}, G_0 = -i w^3),
    P = 1, eps = 2,  eta = lim rho_T/w = 1, s = 4 pi  =>  eta/s = 1/(4 pi)  (printed check).
Hence  S = (1/pi) Int_0^inf dw [ |phi_h/phi_bdy|^2 - w^2 ].

Method: horizon side, phi = f^{-i w/3} psi, psi regular (4-term Frobenius start at z = 1 - h), integrated to z_m = 1/2.
Boundary side: the two exact Frobenius solutions u1 = 1 + w^2 z^2/2 + 0*z^3 + ..., u2 = z^3 + ... (indices 0, 3,
no logarithm in AdS4; recursion derived in the docstring of bdy_series) summed at z0 = min(0.25, 1.5/w) and
integrated to z_m.  phi_bdy = A with phi = A u1 + B u2, A = W(phi, u2)/W(u1, u2) at z_m.
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# ---------------- horizon Frobenius data (sympy, once) ----------------
z, w, u = sp.symbols('z w u')
f = 1 - z**3
Bpsi1 = 2*sp.I*w*z**2 - 3*z**2 - 2*f/z               # coefficient of psi'
Bpsi0 = w**2*(1 + z)*(1 + z**2)/(1 + z + z**2)       # coefficient of psi
N = 5
cs = sp.symbols('c1:%d' % (N + 1))
ser = 1 + sum(cs[k]*u**(k + 1) for k in range(N))    # psi in u = 1 - z
# psi_zz = psi_uu, psi_z = -psi_u
expr = f.subs(z, 1 - u)*sp.diff(ser, u, 2) - Bpsi1.subs(z, 1 - u)*sp.diff(ser, u) + Bpsi0.subs(z, 1 - u)*ser
expr = sp.series(sp.together(expr*(1 + (1 - u) + (1 - u)**2)), u, 0, N).removeO()  # clear the (1+z+z^2) denominator
sol = {}
for k in range(N):
    eq = sp.expand(expr.coeff(u, k)).subs(sol)
    sol[cs[k]] = sp.simplify(sp.solve(eq, cs[k])[0])
cfun = [sp.lambdify(w, sol[cs[k]], 'numpy') for k in range(N)]


def rhs_psi(zz, y, wv):
    ff = 1 - zz**3
    b1 = 2j*wv*zz**2 - 3*zz**2 - 2*ff/zz
    b0 = wv**2*(1 + zz)*(1 + zz**2)/(1 + zz + zz**2)
    return [y[1], -(b1*y[1] + b0*y[0])/ff]


def rhs_phi(zz, y, wv):
    ff = 1 - zz**3
    return [y[1], -((-3*zz**2 - 2*ff/zz)*y[1] + wv**2*y[0]/ff)/ff]


def bdy_series(wv, z0, a0, a3, M=120):
    """phi = sum a_n z^n; from z f^2 phi'' + f(z f' - 2 f) phi' + w^2 z phi = 0 with f = 1 - z^3:
    (m+1)(m-2) a_{m+1} + w^2 a_{m-1} + (m-2)(7-2m) a_{m-2} + (m-5)^2 a_{m-5} = 0;  a1 = 0, a2 = w^2 a0/2, a3 free."""
    a = np.zeros(M + 1, dtype=float)
    a[0] = a0; a[1] = 0.0; a[2] = wv**2*a0/2; a[3] = a3
    for m in range(3, M):
        s = wv**2*a[m - 1] + (m - 2)*(7 - 2*m)*a[m - 2] + ((m - 5)**2*a[m - 5] if m >= 5 else 0.0)
        a[m + 1] = -s/((m + 1)*(m - 2))
    n = np.arange(M + 1)
    val = np.sum(a*z0**n)
    der = np.sum(n[1:]*a[1:]*z0**(n[1:] - 1))
    return val, der


def ratio2(wv, h=1e-3, zm=0.5, rtol=1e-12):
    c = [complex(cf(wv)) for cf in cfun]
    p0 = 1 + sum(c[k]*h**(k + 1) for k in range(N))
    dp = -sum((k + 1)*c[k]*h**k for k in range(N))
    s = solve_ivp(rhs_psi, (1 - h, zm), [p0, dp], args=(wv,), method='DOP853', rtol=rtol, atol=1e-15)
    psi, dpsi = s.y[0, -1], s.y[1, -1]
    ffm = 1 - zm**3
    F = ffm**(-1j*wv/3)
    g = 1j*wv*zm**2/ffm
    phi, dphi = F*psi, F*(dpsi + g*psi)
    z0 = min(0.25, 1.5/max(wv, 1e-9))
    us = []
    for (a0, a3) in ((1.0, 0.0), (0.0, 1.0)):
        v, d = bdy_series(wv, z0, a0, a3)
        r = solve_ivp(rhs_phi, (z0, zm), [v + 0j, d + 0j], args=(wv,), method='DOP853', rtol=rtol, atol=1e-300)
        us.append((r.y[0, -1], r.y[1, -1]))
    (u1, du1), (u2, du2) = us
    A = (phi*du2 - dphi*u2)/(u1*du2 - du1*u2)
    return abs(1.0/A)**2           # |phi_h|^2 / |phi_bdy|^2 with psi(1) = 1


if __name__ == '__main__':
    print('eta check: |phi_h/phi_b|^2 at w=1e-3 =', ratio2(1e-3), '(expect 1 => eta/s = 1/4pi)')
    I = lambda wv: ratio2(wv) - wv**2
    for (hh, zm) in ((1e-3, 0.5), (5e-4, 0.6)):
        print('stability h=%g zm=%g: I(3)=%.12e I(12)=%.12e' % (hh, zm, ratio2(3.0, hh, zm) - 9, ratio2(12.0, hh, zm) - 144))
    for wv in (0.5, 1, 2, 3, 5, 8, 10, 12, 15, 20, 25):
        print('w=%5.1f  I(w)=% .12e' % (wv, I(wv)))
    from scipy.integrate import quad
    pieces = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 30)]
    tot = 0.0
    for (a, b) in pieces:
        val, err = quad(I, a, b, limit=200, epsabs=1e-12, epsrel=1e-11)
        tot += val
        print('Int_%g^%g I = % .12e (err %.1e)' % (a, b, val, err))
    # OPE tail: next thermal operator in T x T is the double trace [TT] (Delta = 6) => delta G ~ w^-3, I ~ w^-4
    ws = np.array([20.0, 25.0, 30.0, 40.0]); Is = np.array([I(x) for x in ws])
    pfit = np.polyfit(np.log(ws), np.log(-Is), 1)
    print('tail fit on w in [20,40]: I ~ -%.4f w^(%.4f)   (OPE expectation: exponent -4)' % (np.exp(pfit[1]), pfit[0]))
    c4 = -Is[-1]*ws[-1]**4
    tail = -c4/(3*30.0**3)
    print('tail Int_30^inf (c w^-4 with c from w=40) = %.3e' % tail)
    S0 = tot/np.pi; S = (tot + tail)/np.pi
    print('S_eta (Einstein) without tail = %.10f, with tail = %.10f ; prediction 3/8 ; diff = %.2e' % (S0, S, S - 0.375))
