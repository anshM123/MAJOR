"""Flat thermal one-loop bubble of the critical O(N) model at N = infinity (beta = 1).

Pi0(n, u) = T sum_{n1} int d^2q/(2pi)^2  G(n1,q) G(n-n1,k-q),  u = k^2, G = 1/(w^2+q^2+m0^2),
m0 = 2 log(golden ratio) (critical thermal mass, 2 sinh(m0/2) = 1).

Pi0 = Pi_vac(p) + Pi_th(n,u),  Pi_vac = arctan(p/2m)/(4 pi p),  p^2 = u + w_n^2.
Pi_th (n != 0) = -(1/pi) int_0^inf q dq n_B(E)/E Re[1/S],  S = sqrt(a^2-b^2) on the branch continuous with a,
                 a = -(u+w^2) + 2 i w E,  b^2 = 4 u q^2.
Pi_th (n == 0) = (1/4pi) int_0^{pi/2} dth sin(th) n_B(E)/E,  E^2 = u sin^2(th)/4 + m0^2.
Derivatives in u are taken analytically under the integral.
"""
from mpmath import mp, mpf, mpc, sqrt, log, exp, atan, pi, quad, inf, re, sin, cos, expm1

mp.dps = 30
PHI = (1 + sqrt(5)) / 2
M0 = 2 * log(PHI)


def nB(E):
    return 1 / expm1(E)


def pi_vac(u, w, m=None, d=0):
    """vacuum 3D bubble arctan(p/2m)/(4 pi p) and its u-derivatives (d = 0,1,2)."""
    m = M0 if m is None else m
    p2 = u + w * w
    p = sqrt(p2)
    if p < m / 2:   # series in x^2 = p^2/(4 m^2): atan(x)/x = sum (-1)^k x^{2k}/(2k+1)
        X = p2 / (4 * m * m)
        tot = mpf(0)
        for k in range(d, 80):
            c = (-1) ** k / mpf(2 * k + 1)
            fac = 1
            for j in range(d):
                fac *= (k - j)
            tot += c * fac * X ** (k - d)
        return tot / (4 * m * m) ** d / (8 * pi * m)
    f0 = atan(p / (2 * m)) / (4 * pi * p)
    if d == 0:
        return f0
    # g(p) = atan(p/2m)/(4 pi p); dg/dp, d2g/dp2
    A = atan(p / (2 * m))
    D = 1 / (1 + (p / (2 * m)) ** 2) / (2 * m)          # dA/dp
    DD = -2 * (p / (2 * m)) / (1 + (p / (2 * m)) ** 2) ** 2 / (2 * m) ** 2
    g1 = (D / p - A / p2) / (4 * pi)
    g2 = (DD / p - 2 * D / p2 + 2 * A / (p2 * p)) / (4 * pi)
    # d/du = (1/2p) d/dp
    if d == 1:
        return g1 / (2 * p)
    # d2/du2 = (1/4p^2) g2 - (1/4p^3) g1
    return g2 / (4 * p2) - g1 / (4 * p2 * p)


def _branch(a, b2):
    S = a * sqrt(1 - b2 / (a * a))
    return S


def pi_th(n, u, m=None, d=0):
    m = M0 if m is None else m
    u = mpf(u)
    if n == 0:
        def f(th):
            s2 = sin(th) ** 2
            E = sqrt(u * s2 / 4 + m * m)
            g = nB(E) / E
            if d == 0:
                return sin(th) * g
            # dg/dE
            e = exp(E)
            gE = -e / (e - 1) ** 2 / E - g / E
            Eu = s2 / (8 * E)
            if d == 1:
                return sin(th) * gE * Eu
            gEE = (-(e * (e - 1) ** 2 - e * 2 * (e - 1) * e) / (e - 1) ** 4 / E
                   + e / (e - 1) ** 2 / E ** 2 - gE / E + g / E ** 2)
            Euu = -s2 / (8 * E * E) * Eu
            return sin(th) * (gEE * Eu * Eu + gE * Euu)
        sc = 2 / sqrt(u + 1)
        pts = sorted(set([mpf(0)] + [sc * 2 ** k for k in range(-1, 9) if sc * 2 ** k < pi / 2] + [pi / 4, pi / 2]))
        return quad(f, pts) / (4 * pi)
    w = 2 * pi * n

    def f(q):
        E = sqrt(q * q + m * m)
        a = -(u + w * w) + 2j * w * E
        b2 = 4 * u * q * q
        S = _branch(a, b2)
        if d == 0:
            val = 1 / S
        elif d == 1:
            val = (a + 2 * q * q) / S ** 3
        else:
            val = -1 / S ** 3 + 3 * (a + 2 * q * q) ** 2 / S ** 5
        return q * nB(E) / E * re(val)
    return -quad(f, [0, 1, 4, 12, 30, 80]) / pi


def Pi0(n, u, d=0):
    w = 2 * pi * n
    return pi_vac(mpf(u), w, d=d) + pi_th(n, u, d=d)


if __name__ == "__main__":
    import sys
    # check 1: against direct 2D integral of the closed-form Matsubara kernel
    from mpmath import quad as q2

    def K(w, E1, E2):
        n1, n2 = nB(E1), nB(E2)
        d = E1 - E2
        t2 = (n2 - n1) * 2 * d / (w * w + d * d) if (w != 0 or abs(d) > mpf(10)**-12) else 2 * exp(E1) / expm1(E1) ** 2
        return ((1 + n1 + n2) * 2 * (E1 + E2) / (w * w + (E1 + E2) ** 2) + t2) / (4 * E1 * E2)

    mp.dps = 15
    for n, u in [(0, mpf('0.7')), (1, mpf('2.3')), (2, mpf(30))]:
        w = 2 * pi * n
        k = sqrt(u)

        def integrand(q, phi):
            E1 = sqrt(q * q + M0 ** 2)
            E2 = sqrt(q * q + k * k - 2 * q * k * cos(phi) + M0 ** 2)
            return q * K(w, E1, E2)
        direct = q2(integrand, [0, 2, 8, 40, inf], [0, pi]) * 2 / (4 * pi ** 2)
        print('n,u', n, u, 'direct', direct, 'formula', Pi0(n, u))
