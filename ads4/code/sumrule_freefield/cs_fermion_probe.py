"""
Feasibility probe for a Chern-Simons-fermion test of the shear/Hall sum rules.

Input: Ghosh-Mazumdar, arXiv:1912.06589 (v2), eq. (strf1a): the renormalised large-N, all-lambda thermal correlator
    G(q3) = < J_{--}(-q) J_{++}(q) >_beta ,   q^+- = 0 (zero spatial momentum), q3 = Matsubara frequency,
regular-fermion theory, with (their eqs. (propterms), (j0rf1), (j0rf2)), beta = 1:
    mu_F = lambda mu_F + (1/(pi i)) [Li2(-e^{-mu_F - i pi lambda}) - c.c.]            (lambda > 0 branch)
    F_F(q3, mu) = Int_0^inf dx x H_F,
    H_F = (4i/pi) sgn(lambda) [log cosh(E/2 - i pi |lambda|/2) - c.c.] / (E (4E^2 + q3^2)),  E = sqrt(x^2 + mu^2).
J_{--} = T1 + i T2 up to normalisation (x^+- = (x1 +- i x2)/sqrt 2), so G(q3) is a helicity correlator; its part odd
in q3 is the Hall (parity-odd) channel h, the even part the shear channel g.

What is probed (structure only; the overall normalisation is irrelevant here):
  large-|q3| expansion  G = a3 q3^3 + a2 q3^2 + a1 q3 + a0 + O(1/q3).
  a3 = vacuum (T = 0) piece; a1 = a term linear in the frequency.  At N = infinity an operator of dimension Delta with
  a thermal vev contributes <O> q3^{3 - Delta}; psibar psi has Delta = 2 + O(1/N), is parity odd, and has
  <psibar psi>_beta ~ mu_F != 0 for lambda != 0, so a1 != 0 signals the light-scalar obstruction in the Hall channel.
"""
import mpmath as mp

mp.mp.dps = 30


def mu_F(lam):
    f = lambda m: (1 - lam)*m - (1/(mp.pi*1j))*(mp.polylog(2, -mp.e**(-m - 1j*mp.pi*lam))
                                                  - mp.conj(mp.polylog(2, -mp.e**(-m - 1j*mp.pi*lam))))
    return mp.findroot(lambda m: mp.re(f(m)), 1.0)


def im_logcosh(a, b):
    # Im log cosh(a - i b) for 0 <= b < pi/2, continuous branch
    return -mp.atan(mp.tanh(a)*mp.tan(b))


def FF(q3, mu, lam):
    b = mp.pi*abs(lam)/2

    def H(x):
        E = mp.sqrt(x*x + mu*mu)
        return (4j/mp.pi)*mp.sign(lam)*(2j*im_logcosh(E/2, b))/(E*(4*E*E + q3*q3))
    return mp.quad(lambda x: x*H(x), [0, 1, 10, mp.inf])


def G(q3, lam):
    mu = mu_F(lam)
    s = mp.sign(lam)
    N = 1
    e = mp.e**(1j*q3*FF(q3, mu, lam))
    t1 = 1j*N/(128*mp.pi*lam*q3)*(1 - e)*(q3**2 + 4*mu**2)*(q3**2 - 4*mu**2 + 4j*q3*s*mu)
    t2 = 1j*N/(96*mp.pi*lam)*(2*mu**2*(3*q3 - 2j*s*mu) - 3j*q3**2*s*mu)
    t3 = N*mu**3/(12*mp.pi)
    li = mp.quad(lambda x: x*(mp.polylog(2, -mp.e**(-x + 1j*mp.pi*lam)) - mp.conj(mp.polylog(2, -mp.e**(-x + 1j*mp.pi*lam)))),
                 [mu, mu + 5, mp.inf])
    t4 = -1j*N/(4*mp.pi**2*lam)*li
    return t1 + t2 + t3 + t4, mu


if __name__ == '__main__':
    for lam in [mp.mpf('0.02'), mp.mpf('0.3'), mp.mpf('0.6')]:
        mu = mu_F(lam)
        print('lambda = %s : mu_F = %s' % (lam, mp.nstr(mu, 12)))
        # fit on q3 = 2 pi n, n large (both signs), polynomial in q3 with powers 3..-3
        # G(-q3) = conj G(q3): Re G is the even (shear) part, Im G the odd (Hall) part.  Fit q3 > 0 only,
        # separately, in powers of q3 = |q3| (vacuum: |q3|^3 in Re, q3^3 contact in Im).
        ns = [20, 30, 40, 60, 80, 120, 160, 240, 320, 480]
        rows = [(2*mp.pi*n, G(2*mp.pi*n, lam)[0]) for n in ns]
        powers = [3, 2, 1, 0, -1, -2]
        Q0 = 2*mp.pi*100
        A = mp.matrix([[(q/Q0)**p for p in powers] for q, _ in rows])
        for part, fn in [('Re (shear, even)', mp.re), ('Im (Hall, odd)', mp.im)]:
            bvec = mp.matrix([fn(g) for _, g in rows])
            coef, res = mp.qr_solve(A, bvec)
            print('   %s: ' % part + ', '.join('q^%d: %s' % (p, mp.nstr(c/Q0**p, 8)) for p, c in zip(powers, coef))
                  + '   [rel. residual %s]' % mp.nstr(res/mp.norm(bvec), 3))
        print('   compare: N mu_F^2/(16 pi lambda) =', mp.nstr(mu**2/(16*mp.pi*lam), 10),
              ';  vacuum-type a3 vs (1 - e^{+- i pi lambda}) structure')
        for n in [1, 2, 3]:
            g, _ = G(2*mp.pi*n, lam); gm, _ = G(-2*mp.pi*n, lam)
            print('   n=%d: G(+q3) = %s ; G(-q3) = %s' % (n, mp.nstr(g, 10), mp.nstr(gm, 10)))


def stability(lam, ns, powers):
    rows = [(2*mp.pi*n, G(2*mp.pi*n, lam)[0]) for n in ns]
    Q0 = 2*mp.pi*100
    A = mp.matrix([[(q/Q0)**p for p in powers] for q, _ in rows])
    out = []
    for fn in [mp.re, mp.im]:
        coef, res = mp.qr_solve(A, mp.matrix([fn(g) for _, g in rows]))
        out.append({p: c/Q0**p for p, c in zip(powers, coef)})
    return out


if __name__ == '__main__':
    print('stability of the |q3|^1 (shear) and q3^1, q3^0 (Hall) coefficients, lambda = 0.3:')
    for ns, pw in [([30, 45, 60, 90, 120, 180, 240, 360], [3, 2, 1, 0, -1, -2]),
                   ([60, 90, 120, 180, 240, 360, 480, 720, 960], [3, 2, 1, 0, -1, -2, -3]),
                   ([10, 15, 20, 25, 30, 40, 50, 60, 80, 100], [3, 2, 1, 0, -1, -2, -3, -4])]:
        re, im = stability(mp.mpf('0.3'), ns, pw)
        print('   n in [%d..%d], %d powers: Re |q|^1 = %s, Re q^0 = %s ; Im q^1 = %s, Im q^0 = %s'
              % (ns[0], ns[-1], len(pw), mp.nstr(re[1], 9), mp.nstr(re[0], 9), mp.nstr(im[1], 9), mp.nstr(im[0], 9)))
    print('small-lambda scaling of the shear |q3| coefficient a1 (expect ~ lambda^2 if it is <TT psibar psi> x <psibar psi>):')
    for lam in ['0.005', '0.01', '0.02', '0.04', '0.08']:
        lam = mp.mpf(lam)
        re, im = stability(lam, [30, 45, 60, 90, 120, 180, 240, 360], [3, 2, 1, 0, -1, -2])
        mu = mu_F(lam)
        print('   lambda = %s: mu_F = %s, a1 = %s, a1/lambda^2 = %s, a1/mu_F^2 = %s ; Hall b1 = %s, b0 = %s'
              % (mp.nstr(lam, 3), mp.nstr(mu, 8), mp.nstr(re[1], 8), mp.nstr(re[1]/lam**2, 8), mp.nstr(re[1]/mu**2, 8),
                 mp.nstr(im[1], 8), mp.nstr(im[0], 8)))


def eps_CS(lam):
    """energy density (N = 1, T = 1) of the regular CS-fermion theory, eps = 2P = -2F/V,
    Aharony-Giombi-Gur-Ari-Maldacena-Yacoby 1211.4843 eq. (F-nmass-nc-fer) with sigma~ = 0, lambda > 0 branch."""
    mu = mu_F(lam)
    li = mp.quad(lambda y: y*(mp.polylog(2, -mp.e**(-y + 1j*mp.pi*lam)) - mp.conj(mp.polylog(2, -mp.e**(-y + 1j*mp.pi*lam)))),
                 [mu, mu + 5, mp.inf])
    betaF_over_V = (1/(2*mp.pi))*(mu**3/3*(1 - 1/lam) + li/(mp.pi*1j*lam))
    return -2*mp.re(betaF_over_V)


if __name__ == '__main__':
    print('check eps(lambda -> 0) = 3 zeta(3)/(2 pi) = %s :' % mp.nstr(3*mp.zeta(3)/(2*mp.pi), 10),
          mp.nstr(eps_CS(mp.mpf('1e-4')), 10), ';  O(lambda^2) check (1211.4843):',
          mp.nstr((eps_CS(mp.mpf('0.01')) - 3*mp.zeta(3)/(2*mp.pi))/mp.mpf('0.01')**2, 8), 'vs',
          mp.nstr(-2*(mp.pi**2*mp.log(2) + 4*mp.log(2)**3)/(6*mp.pi), 8))
    print('Hall UV constant vs |alpha_T|/32 = sin(pi lambda)/8 (CDP 1707.03007: alpha_T = 4 sin(pi lambda)):')
    print('  lambda    eps         b0 = lim Im G      S_H = b0/(2 eps)   sin(pi lam)/8     ratio     a1/eps (shear |q| term)')
    for lam in ['0.01', '0.05', '0.1', '0.2', '0.3', '0.5', '0.7', '0.9']:
        lam = mp.mpf(lam)
        re, im = stability(lam, [30, 45, 60, 90, 120, 180, 240, 360], [3, 2, 1, 0, -1, -2])
        e = eps_CS(lam)
        SH = im[0]/(2*e)
        pred = mp.sin(mp.pi*lam)/8
        print('  %-6s  %-10s  %-16s  %-16s  %-14s  %-10s  %s' % (mp.nstr(lam, 3), mp.nstr(e, 8), mp.nstr(im[0], 10),
              mp.nstr(SH, 10), mp.nstr(pred, 10), mp.nstr(SH/pred, 10), mp.nstr(re[1]/e, 8)))
    print('Shear-channel UV constant a0 = lim [Re G - a3|q|^3 - a1|q|] vs the t4 term (t4 = -4 cos(pi lambda)):')
    print('  lambda   a0/(2 eps)        (cos(pi lam) - 1)/8 = -1/8 - t4/32     diff')
    for lam in ['0.01', '0.1', '0.3', '0.5', '0.7', '0.9']:
        lam = mp.mpf(lam)
        re, im = stability(lam, [30, 45, 60, 90, 120, 180, 240, 360], [3, 2, 1, 0, -1, -2])
        e = eps_CS(lam)
        v = re[0]/(2*e); pred = (mp.cos(mp.pi*lam) - 1)/8
        print('  %-6s  %-16s  %-16s  %s' % (mp.nstr(lam, 3), mp.nstr(v, 10), mp.nstr(pred, 10), mp.nstr(v - pred, 3)))
