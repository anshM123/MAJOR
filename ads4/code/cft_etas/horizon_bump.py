"""
Bulk mechanism behind Theorem C: Matsubara modes vanish at the Euclidean horizon (tip of the cigar).

Remark-8 class: shear action -(N/2) Int [A phi'^2 - B phidot^2], A = f0 a/u^2, B = b/(u^2 f0), f0 = 1-u^3, u = z.
Deformation a = b = 1 + Amp*chi(u), chi = 1 on u >= 1-w/2, 0 on u <= 1-w, C^infinity in between.
Near the horizon a = b = const, so the ingoing exponent and the Frobenius data are Einstein's; the boundary region is untouched,
so G(0) and the UV/OPE data are unchanged.  Entropy is held fixed (nu = 0 in Remark 8), so 4 pi eta/s - 1 = c(1) - 1 = Amp.
E(s) = -dG(i y)/y with y = 1.5 s (dG(0) = 0 in this contact scheme), F = E - E_Einstein.
Prints F(0+), F(n) for n = 1..3 and F(-gamma) (Theorem C's strip norm) as the bump width w shrinks.
"""
import numpy as np
from scipy.integrate import solve_ivp
import einstein_ref as er


def chi(u, w):
    t = (u - (1 - w))/(w/2)
    if t <= 0: return 0.0
    if t >= 1: return 1.0
    e = lambda x: np.exp(-1/x) if x > 0 else 0.0
    return e(t)/(e(t) + e(1 - t))


def dchi(u, w, h=1e-7):
    return (chi(u + h, w) - chi(u - h, w))/(2*h)


def Gdef(wv, Amp, w, h=1e-4, rtol=1e-11):
    """3B/A for the deformed problem at complex frequency wv."""
    wv = complex(wv)
    ua = 1 - w/4                      # inside the constant region
    c = [complex(cf(wv)) for cf in er._cfun]
    p0 = 1 + sum(c[k]*h**(k + 1) for k in range(er._N))
    dp = -sum((k + 1)*c[k]*h**k for k in range(er._N))
    s = solve_ivp(er._rhs_psi, (1 - h, ua), [p0 + 0j, dp + 0j], args=(wv,), method='DOP853', rtol=rtol, atol=1e-300)
    psi, dpsi = s.y[0, -1], s.y[1, -1]
    ff = 1 - ua**3
    F = ff**(-1j*wv/3); g = 1j*wv*ua**2/ff
    phi, dphi = F*psi, F*(dpsi + g*psi)
    # deformed equation: (a f0 phi'/u^2)' + w^2 b phi/(u^2 f0) = 0, a = b = 1 + Amp chi
    def rhs(u, y):
        f0 = 1 - u**3
        a = 1 + Amp*chi(u, w); da = Amp*dchi(u, w)
        # a f0 phi''/u^2 + (da f0/u^2 + a (f0/u^2)') phi' + w^2 a phi/(u^2 f0) = 0
        dP = (-3*u**2*u**2 - 2*u*f0)/u**4
        return [y[1], -((da*f0/u**2 + a*dP)*y[1] + wv**2*a*y[0]/(u**2*f0))/(a*f0/u**2)]
    zm = min(0.5, 1 - w)
    r = solve_ivp(rhs, (ua, zm), [phi, dphi], method='DOP853', rtol=rtol, atol=1e-300, max_step=w/40)
    phi, dphi = r.y[0, -1], r.y[1, -1]
    z0 = min(0.25, 1.5/max(abs(wv), 1e-9), zm*0.9)
    us = []
    for (a0, a3) in ((1.0, 0.0), (0.0, 1.0)):
        v, d = er._bdy_series(wv, z0, a0, a3)
        rr = solve_ivp(er._rhs_phi, (z0, zm), [v + 0j, d + 0j], args=(wv,), method='DOP853', rtol=rtol, atol=1e-300)
        us.append((rr.y[0, -1], rr.y[1, -1]))
    (u1, du1), (u2, du2) = us
    W = u1*du2 - du1*u2
    A = (phi*du2 - dphi*u2)/W
    B = (u1*dphi - du1*phi)/W
    return 3*B/A


def E(s, Amp, w):
    y = 1.5*s
    wv = 1j*y
    return -(Gdef(wv, Amp, w) - 1j*wv**3).real/y


if __name__ == '__main__':
    Amp = 0.1
    gam = 0.5
    print('# a=b=1+Amp*chi, Amp=%.2f, entropy fixed => exact 4pi eta/s - 1 = Amp' % Amp)
    print('%-7s %-12s %-12s %-12s %-12s %-12s' % ('w', 'F(1e-4)', 'F(1)', 'F(2)', 'F(3)', 'F(-0.5)'))
    for w in (0.4, 0.2, 0.1, 0.05, 0.025):
        vals = [E(s, Amp, w) - E(s, 0.0, w) for s in (1e-4, 1, 2, 3, -gam)]
        print('%-7.3f ' % w + ' '.join('%-12.4e' % v for v in vals))
    print('expected scalings: F(0)=Amp; F(n) ~ w^n ; F(-gamma) ~ w^(-gamma)')
