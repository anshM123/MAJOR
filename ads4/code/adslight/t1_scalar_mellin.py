"""t1: scalar warm-up (ADS4_LIGHT.md section 1).  d = 3, identical scalars phi of dimension Df.

Implements CHMRS (arXiv:2106.10274) eq. (CfromBMellin) for k = 2:
   C_{2,nu} = 2^{d-2-2k} pi^{3/2} Gamma((d-1)/2) a_Df / gamma_{2Df+1}(nu)^2 * int dT/(4 pi i) Gamma(T/2)^2 gamma_{2Df+1-T}(nu) Bhat_{2,T}
with Bhat_{2,T} (eq. BMellinDefinition) = (-1) * (residue at S' = infinity, counter-clockwise) of M(S', U)/[(2Df-S')(Df-S'/2)(Df-T'/2)],
U = 2Df - T fixed.  For M = A(U) S^2 + O(S) the only contribution is Bhat_{2,T} = -4 A(2Df - T)  [D, elementary].
With T = 2z, c = Df + 1/4 - i nu/2 (d = 3):
   -C_{2,nu}[A] = (pi^{3/2} a_Df / 2) / |Gamma(c)|^4 * int_{Re z = z0} dz/(2 pi i) Gamma(z)^2 Gamma(c - z) Gamma(cbar - z) A(2Df - 2z).

Checks (all [DV] in ADS4_LIGHT.md):
 (a) contact M_{1,0} = c_2 (S^2+T^2+U^2):  -C_2 = 2 for all nu   (CHMRS (C2g2)); fixes the orientation of Bhat.
 (b) contact M_{0,1} = c_3 S T U:          -C_2 = nu^2 - (2Df+1/2)(2Df-1/2)  (CHMRS (C2C4Light), d = 3).
 (c) graviton exchange, Regge form CHMRS (MGravityRegge): A(U) = sum_m R_m/(U - 1 - 2m) -- i.e. ONLY the stress-tensor
     single-trace poles (primary m = 0 and its descendants): the 'T block' in Mellin space.  Claim (CHMRS (BknuGravity)):
     -C_{2,nu} = 8 pi G/(nu^2 + 9/4) exactly.  Checked at several nu, Df; the m = 0 term alone is NOT of this form.
 (d) H_2 propagator identity CHMRS (gravitonpropagator), d = 3: int_0^oo dnu/(2pi) rho(nu) P_{-1/2-i nu}(g)/(nu^2+9/4)
     = (vol S^1/2) h_2(g-1); rho = pi nu tanh(pi nu).  Shows 9/4 <-> transverse field of dimension 2 = d-1 (graviton), and that
     1/(nu^2+1/4) would be dimension 1.  Numerically the left side equals h_{D0}(g-1) itself (ratio 1, not pi = vol S^1/2): the
     constant in our reading of CHMRS's normalisation differs by pi; only the shape (pole <-> D0) is used.
 (e) dictionary: 8 pi G d Gamma(d/2)^3/(4 pi^{d/2} Gamma(d+2)) = c_T^free/c_T  (CHMRS) => d = 3: 8 pi G = 256/c_T, i.e.
     G = 32/(pi c_T) = DICTIONARY_T4.md (c_T = 32 L^2/(pi G)).  And f^2_{phi phi T} = 3 Df^2/(8 c_T).
"""
import mpmath as mp

mp.mp.dps = 25
d = 3


def a_D(Df):
    return mp.pi ** ((d - 3) / 2.) * mp.gamma(Df) ** 2 * mp.gamma(Df - (d - 2) / 2.) ** 2 / (2 ** (d - 5) * mp.gamma((d - 1) / 2.))


def c_n(n, Df):
    return 2 ** (n - 3) * mp.gamma(2 * Df - d / 2. + n) / (mp.pi ** (d / 2.) * mp.gamma(Df) ** 2 * mp.gamma(Df - (d - 2) / 2.) ** 2)


def minusC2(Afun, Df, nu, z0):
    c = Df + mp.mpf(1) / 4 - 1j * nu / 2
    cb = mp.conj(c)
    pref = mp.pi ** 1.5 * a_D(Df) / 2 / abs(mp.gamma(c)) ** 4

    def f(y):
        z = z0 + 1j * y
        return mp.gamma(z) ** 2 * mp.gamma(c - z) * mp.gamma(cb - z) * Afun(2 * Df - 2 * z)
    val = mp.quad(f, [-mp.inf, -20, -5, 0, 5, 20, mp.inf]) / (2 * mp.pi)   # dz/(2 pi i) = dy/(2 pi)
    return pref * val


def moments_z(Df, nu, z0, jmax):
    """mu_j = int dz/(2pi i) Gamma(z)^2 Gamma(c-z) Gamma(cbar-z) z^j on Re z = z0."""
    c = Df + mp.mpf(1) / 4 - 1j * nu / 2
    cb = mp.conj(c)
    out = []
    for j in range(jmax + 1):
        def f(y, j=j):
            z = z0 + 1j * y
            return mp.gamma(z) ** 2 * mp.gamma(c - z) * mp.gamma(cb - z) * z ** j
        out.append(mp.quad(f, [-mp.inf, -20, -5, 0, 5, 20, mp.inf]) / (2 * mp.pi))
    return out


def R_m(m, Df, G):
    return (-8 * mp.pi * G * Df ** 2 * mp.sin(mp.pi / 2 * (d - 2 * Df)) ** 2 * mp.gamma((d - 2 * Df) / 2. + m) ** 2
            / (4 * mp.pi ** (d / 2. + 2) * mp.factorial(m) * mp.gamma(d / 2. + m + 1)))


def grav_minusC2(Df, nu, G, m_direct=12, m_max=20000, jmax=10, return_m0=False):
    z0 = Df - mp.mpf(1) / 4          # between the T-poles z_m = Df - 1/2 - m (left) and Gamma(c - z) poles (right)
    c = Df + mp.mpf(1) / 4 - 1j * nu / 2
    pref = mp.pi ** 1.5 * a_D(Df) / 2 / abs(mp.gamma(c)) ** 4
    tot = mp.mpf(0)
    m0 = None
    for m in range(m_direct):
        val = minusC2(lambda U, m=m: R_m(m, Df, G) / (U - (d - 2) - 2 * m), Df, nu, z0)
        if m == 0:
            m0 = val
        tot += val
    mu = moments_z(Df, nu, z0, jmax)
    # tail: 1/(2Df - 2z - 1 - 2m) = sum_j (2z)^j / w^{j+1}, w = 2Df - 1 - 2m
    tail = mp.mpf(0)
    for m in range(m_direct, m_max):
        w = 2 * Df - 1 - 2 * m
        s = sum(2 ** j * mu[j] / w ** (j + 1) for j in range(jmax + 1))
        tail += R_m(m, Df, G) * s
    tot += pref * tail
    return (tot, m0) if return_m0 else tot


def h_prop(D0, u):
    """CHMRS h_{Delta0}(u), d = 3 (AdS_{d-1} = H_2 bulk-to-bulk propagator)."""
    return (mp.gamma(D0) / (2 * mp.pi ** ((d - 2) / 2.) * mp.gamma((2 - d) / 2. + D0 + 1)) * (2 * u) ** (-D0)
            * mp.hyp2f1(D0, D0 - (d - 3) / 2., 2 * D0 - d + 3, -2 / u))


def main():
    print('=== (a),(b) contact calibrations, d = 3 ===')
    for Df in [1.2, 1.75, 2.6]:
        for nu in [0, 0.7, 3.0, 9.0]:
            z0 = mp.mpf(Df) / 2
            g2 = minusC2(lambda U: 2 * c_n(2, Df), Df, nu, z0)
            g3 = minusC2(lambda U: -U * c_n(3, Df), Df, nu, z0)
            exp3 = nu ** 2 - (2 * Df + 0.5) * (2 * Df - 0.5)
            print('Df=%.2f nu=%4.1f  -C2[M10] = %s (exp 2)   -C2[M01] = %s (exp %.10f)  diff %.1e'
                  % (Df, nu, mp.nstr(mp.re(g2), 15), mp.nstr(mp.re(g3), 15), exp3, abs(g3 - exp3)))
    print('\n=== (c) graviton exchange = stress-tensor single-trace Mellin poles; G = 1/(8 pi) so target = 1/(nu^2+9/4) ===')
    G = 1 / (8 * mp.pi)
    for Df in [1.2, 1.75, 2.6]:
        for nu in [0, 0.5, 1.5, 4.0, 12.0]:
            tot, m0 = grav_minusC2(Df, nu, G, return_m0=True)
            tgt = 1 / (mp.mpf(nu) ** 2 + mp.mpf(9) / 4)
            print('Df=%.2f nu=%5.1f  sum_m = %s  1/(nu^2+9/4) = %s  rel.diff %.1e   [m=0 term alone / target = %s]'
                  % (Df, nu, mp.nstr(mp.re(tot), 14), mp.nstr(tgt, 14), abs(tot / tgt - 1), mp.nstr(mp.re(m0 / tgt), 8)))
    print('\n=== (d) H_2 propagator: int dnu/(2pi) rho(nu) P_{-1/2-i nu}(g)/(nu^2+a^2) vs h_{D0}(g-1), D0 = 1/2 + a ===')
    for a, D0 in [(1.5, 2), (0.5, 1)]:
        for g in [1.5, 3.0, 10.0]:
            f = lambda nu: nu * mp.tanh(mp.pi * nu) * mp.legenp(-0.5 - 1j * nu, 0, g, type=3) / (nu ** 2 + a ** 2)
            lhs = mp.quad(f, [0, 5, 20, 60, 150]) / (2 * mp.pi)
            rhs = h_prop(D0, g - 1)
            print('a^2=%.2f (D0=%d) g=%5.1f  lhs=%s  rhs=%s  ratio=%s' % (a * a, D0, g, mp.nstr(mp.re(lhs), 10), mp.nstr(rhs, 10),
                                                                         mp.nstr(mp.re(lhs) / rhs, 8)))
    print('   (oscillatory tail truncated at nu = 150: agreement to ~1e-3 is the expected accuracy of this check;\n'
          '    the ratio is the same for both D0 and all g, so the pole position <-> D0 = 1/2 + a identification is what is checked)')
    print('\n=== (e) dictionary ===')
    kap = d * mp.gamma(d / 2.) ** 3 / (4 * mp.pi ** (d / 2.) * mp.gamma(d + 2))
    print('8 pi G * %s = 1/c_T  =>  8 pi G = %s / c_T ;  G = 32/(pi c_T) <=> c_T = 32 L^2/(pi G) [DICTIONARY_T4]' % (mp.nstr(kap, 12), mp.nstr(1 / kap, 12)))
    Df = mp.mpf(1.75)
    f2 = 8 * mp.pi * Df ** 2 * mp.gamma((d + 2) / 2.) ** 3 / (2 * mp.pi ** (d / 2.) * (d - 1) * d * mp.gamma(d + 2))
    print('f^2_{phi phi T}/(8 pi G) at Df=1.75: %s ; with 8piG = 256/c_T: f^2 c_T = %s ; 3 Df^2/8 = %s'
          % (mp.nstr(f2 / (8 * mp.pi), 12), mp.nstr(f2 / (8 * mp.pi) * 256, 12), mp.nstr(3 * Df ** 2 / 8, 12)))


if __name__ == '__main__':
    main()
