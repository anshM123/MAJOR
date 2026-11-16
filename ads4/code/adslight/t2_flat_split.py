"""t2: flat-space low-energy sides of B_2..B_5 split by channel, and the flat-limit power counting used in ADS4_LIGHT.md section 2.

(i) [D, sympy, exact] f_low(s,u) = 8 pi G/(s t u) + 2 pi G s u |g3|^2/t + g4 + g5 t + g6 t^2 - g6' s u   (arXiv:2201.06602 (2.7a)),
    t = -s-u, u = -p^2.  Partial fractions in s at fixed u:  8 pi G/(s t u) = -(8 pi G/u^2)(1/s + 1/t).
    'pm' part = poles of f in its s-channel (heavy analogue: +- states, |c^{+-}|^2),
    'pp' part = poles of f in its t-channel (heavy analogue: ++ states, |c^{++}|^2).
    B_k^{(1)} per (2.26a,b); -B|_low = +(1/2) sum of residues of the (2.26) integrand at s = 0, p^2 (overall sign fixed by reproducing (2.31)-(2.32)).
(ii) the AdS-corrected ansatz of ADS4_LIGHT.md (2.6)-(2.9) inserted in Omega_g3 symbolically: which unknown AdS coefficients
    actually enter (moments of psi_2, psi_3 that vanish kill them).
"""
import sympy as sp

s, p, G, g3, g4, g5, g6, g6p = sp.symbols('s p G g3sq g4 g5 g6 g6p')
u = -p ** 2


def f_parts(x, uu):
    t = -x - uu
    pm = -8 * sp.pi * G / uu ** 2 / x                       # s-channel graviton pole
    pp = -8 * sp.pi * G / uu ** 2 / t + 2 * sp.pi * G * x * uu * g3 / t   # t-channel poles (Einstein + |g3|^2)
    contact = g4 + g5 * t + g6 * t ** 2 - g6p * x * uu
    return pm, pp, contact


def minusB(k, which):
    """-B_k|_low from the requested part(s) of f."""
    def F(x):
        pm, pp, c = f_parts(x, u)
        return {'pm': pm, 'pp': pp, 'c': c, 'all': pm + pp + c}[which]
    if k % 2 == 0:
        integrand = (2 * s - p ** 2) / (s * (s - p ** 2)) ** sp.Rational(k - 2, 2) * (F(s) + F(p ** 2 - s))
    else:
        integrand = (s * (s - p ** 2)) ** (-sp.Rational(k - 3, 2)) * (-F(s) + F(p ** 2 - s))
    integrand = sp.together(sp.simplify(integrand))
    res = sum(sp.residue(integrand, s, s0) for s0 in [0, p ** 2])
    return sp.simplify(sp.Rational(1, 2) * res)


def main():
    print('=== (i) flat low sides, split pm (s-channel poles) / pp (t-channel poles) / contacts ===')
    for k in [2, 3, 4, 5]:
        tot = minusB(k, 'all')
        print('k=%d  -B_k|low = %s' % (k, sp.expand(tot)))
        for w in ['pm', 'pp', 'c']:
            print('        %-2s : %s' % (w, sp.expand(minusB(k, w))))
    print('  (the overall sign convention is the one for which the totals equal FLAT_GRAVITON_BOUND.md (2.31)-(2.32))')

    print('\n=== (ii) Omega_g3 light side with the AdS-corrected forms of ADS4_LIGHT.md section 2 ===')
    q = sp.symbols('q', positive=True)     # q = nu/Lam = p
    psi2 = (1 - q) ** 3 * q * (65 + 155 * q + 47 * q ** 2)
    psi3 = 2 * (1 - q) ** 3 * (81 - 1927 * q + 4601 * q ** 2)
    mu2 = [sp.integrate(psi2 * q ** j, (q, 0, 1)) for j in range(8)]
    mu3 = [sp.integrate(psi3 * q ** j, (q, 0, 1)) for j in range(8)]
    print('moments mu2_j = int psi2 p^j :', mu2[:7])
    print('moments mu3_j = int psi3 p^j :', mu3[:7])
    k0 = 20 - 2 * mu2[6] + 2 * mu3[4]
    print('k(0) = 20 - 2 mu2_6 + 2 mu3_4 =', k0, '=', sp.N(k0, 12), '(expected 330151/12012)')
    Lam, b2, b1, b0, f1, f0, c1, t4 = sp.symbols('Lambda b2 b1 b0 f1 f0 c1 t4')
    # -C^+_2 |g3|^2-part = 2 pi G |g3|^2 (nu^6 + b2 nu^4 + b1 nu^2 + b0);  -C^-_3 = -4 pi G |g3|^2 (nu^4 + f1 nu^2 + f0);
    # -d_{nu^2} C^+_4|_0 = 4 pi G |g3|^2 (exact);  -C^-_5|_0 = 0 (|g3|^2-part).  nu = Lam q;  |g3|^2 = x^2/Lam^8.
    x2 = sp.symbols('x2')
    g3sq = x2 / Lam ** 8
    nu = Lam * q
    part2 = Lam ** 4 * sp.integrate(psi2 * 2 * sp.pi * G * g3sq * (nu ** 6 + b2 * nu ** 4 + b1 * nu ** 2 + b0), (q, 0, 1))
    part3 = sp.Rational(1, 2) * Lam ** 6 * sp.integrate(psi3 * (-4) * sp.pi * G * g3sq * (nu ** 4 + f1 * nu ** 2 + f0), (q, 0, 1))
    part4 = -5 * Lam ** 10 * 4 * sp.pi * G * g3sq
    lin = Lam ** 4 * sp.integrate(psi2 * 16 * sp.pi * G * c1 * t4, (q, 0, 1))
    tot = sp.expand((part2 + part3 + part4) / (sp.pi * G * x2 * Lam ** 2))
    print('-Omega_light |g3|^2-part / (pi G x^2 Lam^2) = -k_AdS(Lam) with k_AdS =', sp.collect(-tot, Lam))
    print('linear-t4 part of -Omega_light (u-channel L=0 contact, c1 unknown) =', sp.simplify(lin), ' [t4 = -240 Re g3 (AdS units)]')


if __name__ == '__main__':
    main()
