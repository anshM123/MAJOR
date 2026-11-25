"""kb_celestial.py -- celestial (rest-frame) kernels of the CLS-D k = 2 (j'=0) and k = 3 (j'=1) functionals and the k = 3 light
function (ads4/K_B2A4.md Sec. 2-3).

In the gauge x0 = inf, after Parseval in x24 and gauge-fixing P = (1,0,0) (rest frame of the s-channel momentum), the celestial
structure <P~_2(z1) P~_2(z3) P_delta(z0)> integrated over z0 (dual factor (z0.P)^(delta-1) is constant in this frame) is
    K^(j')(phi) = G_{1+delta/2}(phi) int_0^2pi dphi0 G_{-delta/2}(-phi0) G_{-delta/2}(phi0-phi),   phi = phi3 - phi1,
    G_s(x) = |2 sin(x/2)|^{2s}  (j'=0)   or   sgn(2 sin(x/2)) |2 sin(x/2)|^{2s}  (j'=1, antiperiodic),   delta = 1/2 + i nu.
Its Fourier coefficients k(m) = int e^{i m phi} K = (2pi)^2 sum_n b_n(s0)^2 b_{n-m}(s1), n in Z (j'=0) or Z+1/2 (j'=1),
b_n(s) = e^{-i pi n} Gamma(2s+1)/(Gamma(1+s-n)Gamma(1+s+n)) = c(s) beta_s(n), beta_s(n) = Gamma(n-s)/Gamma(n+1+s).
k = 3 light function:  -C_3|_{g3}(nu) / -C_3|_{g3}(flat) = S(4)/S(2) / lim = (nu^4 - 35/2 nu^2 + 569/16)/nu^4.
Usage: python3 -B kb_celestial.py [full]   -> cache/kb_celestial.txt   (full: adds direct double-quadrature checks, ~3 min)
"""
import sys, os
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))


def beta(s, n):
    return mp.exp(mp.loggamma(n - s) - mp.loggamma(n + 1 + s))


def Ssum(nu, m, jp, N=3000):
    d = mp.mpf(1)/2 + 1j*mp.mpf(nu)
    s0 = -d/2; s1 = 1 + d/2
    off = mp.mpf(1)/2 if jp == 1 else mp.mpf(0)
    return sum(beta(s0, j + off)**2*beta(s1, j + off - m) for j in range(-N, N))


def pref(nu, jp):
    d = mp.mpf(1)/2 + 1j*mp.mpf(nu)
    s0 = -d/2; s1 = 1 + d/2
    c = (lambda s: -1j*mp.cos(mp.pi*s)*mp.gamma(2*s+1)/mp.pi) if jp == 1 else (lambda s: -mp.sin(mp.pi*s)*mp.gamma(2*s+1)/mp.pi)
    return (2*mp.pi)**2*c(s0)**2*c(s1)


def kdirect(nu, m, jp):
    """direct double quadrature of int e^{i m phi} K^(j')(phi) dphi (independent check)."""
    old = mp.mp.dps; mp.mp.dps = 15
    d = mp.mpf(1)/2 + 1j*mp.mpf(nu); s0 = -d/2; s1 = 1 + d/2

    def G(s, x):
        v = 2*mp.sin(x/2)
        if v == 0:
            return mp.mpf(0)
        return (mp.sign(v) if jp == 1 else 1)*abs(v)**(2*s)

    def K(phi):
        return G(s1, phi)*mp.quad(lambda p0: G(s0, -p0)*G(s0, p0 - phi), [0, phi, 2*mp.pi])
    r = mp.quad(lambda ph: mp.e**(1j*m*ph)*K(ph), [0, mp.pi, 2*mp.pi])
    mp.mp.dps = old
    return r


def dixon(A, B, C):
    g = mp.gamma
    return g(1+A+B+C)*g(1+2*A)*g(1+2*B)*g(1+2*C)/(g(1+A)*g(1+B)*g(1+C)*g(1+A+B)*g(1+B+C)*g(1+C+A))


def G3poly(nu):
    nu = mp.mpf(nu)
    return (nu**4 - mp.mpf(35)/2*nu**2 + mp.mpf(569)/16)/400


def main(full=False):
    mp.mp.dps = 30
    out = []

    def pr(s):
        print(s, flush=True); out.append(s)
    pr('# kb_celestial: rest-frame celestial kernels; k = 3 light function')
    # (a) Fourier representation vs direct quadrature
    if full:
        for nu, m, jp in ((0.7, 2, 1), (1.3, 4, 1), (1.3, 0, 0), (1.3, 4, 0)):
            a = pref(nu, jp)*Ssum(nu, m, jp, N=3000); b = kdirect(nu, m, jp)
            pr('(a) nu=%.1f m=%d j\'=%d: Fourier %s | direct %s | rel diff %.1e' % (nu, m, jp, mp.nstr(a, 12), mp.nstr(b, 12), float(abs(a/b - 1))))
    else:
        pr('(a) direct-quadrature cross-checks: run with "full" (recorded results in K_B2A4.md Sec. 6)')
    # (b) Dixon closed form (j'=0, m=0) vs Fourier
    for nu in (0.9, 2.5):
        d = mp.mpf(1)/2 + 1j*mp.mpf(nu)
        a = pref(nu, 0)*Ssum(nu, 0, 0); b = (2*mp.pi)**2*dixon(1 + d/2, -d/2, -d/2)
        pr('(b) Dixon vs Fourier, nu=%.1f: %s vs %s, rel %.1e' % (nu, mp.nstr(a, 14), mp.nstr(b, 14), float(abs(a/b - 1))))
    # (c) k = 3: S(4)/S(2) (j'=1) against the exact polynomial; truncation N = 1500 vs 3000
    pr('(c) k = 3 light-function shape  G3 = S(4)/S(2) [j\'=1]  vs  (nu^4 - 35/2 nu^2 + 569/16)/400:')
    worst = 0
    for nu in ('0.1', '0.5', '1', '2', '3', '5.5', '7', '12', '20'):
        g1 = Ssum(nu, 4, 1, 1500)/Ssum(nu, 2, 1, 1500); g2 = Ssum(nu, 4, 1, 3000)/Ssum(nu, 2, 1, 3000)
        rel = abs(g2/G3poly(nu) - 1); worst = max(worst, rel)
        pr('    nu=%-4s G3 = %s  (N=1500: rel %.1e; N=3000: rel %.1e; |Im/Re| %.1e)' % (nu, mp.nstr(g2.real, 16), float(abs(g1/G3poly(nu) - 1)),
                                                                              float(rel), float(abs(g2.imag/g2.real))))
    pr('    max rel deviation (N=3000) = %.1e  -> b2 = -35/2, b0 = 569/16 (exact rational identification) [DV]' % float(worst))
    # (d) j'=0 harmonic ratios (structure; used as a cross-check of the j'=1 result)
    pr('(d) j\'=0: k(2)/k(0) vs (4nu^2+1)(4nu^2+25)/1152 ; k(4)/k(2) vs the same G3 polynomial:')
    for nu in ('0.5', '2', '4'):
        z0, z2, z4 = (Ssum(nu, m, 0) for m in (0, 2, 4))
        n_ = mp.mpf(nu)
        pr('    nu=%s: %.2e, %.2e' % (nu, float(abs(z2/z0/((4*n_**2+1)*(4*n_**2+25)/1152) - 1)), float(abs(z4/z2/G3poly(nu) - 1))))
    return out


if __name__ == '__main__':
    o = main(full=(len(sys.argv) > 1 and sys.argv[1] == 'full'))
    open(os.path.join(HERE, 'cache', 'kb_celestial.txt'), 'w').write('\n'.join(o) + '\n')
