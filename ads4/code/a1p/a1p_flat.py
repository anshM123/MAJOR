"""
a1p_flat (A1PRIME.md section 4): flat-space action of the project's functionals on LIGHT single traces (mass m < M, spin J = 0, 2,
coupling ++ ; the parity-odd J = 0, 2 structures couple like ++ at this order), exact rational arithmetic where the weights are exact.

Light-state contributions to the low side (CHLPS 2201.06602 eq. (5.x) "Light spin-0 and spin-2 matter", App. B; f_matter in its
minimal-Regge form):   B2|m = -|g0|^2 (p^2-2m^2) - (|g2|^2/m^4)(p^2-2m^2) P(m,p),   B3|m = |g0|^2 + (|g2|^2/m^4) P(m,p),
                       B4|m = -(|g2|^2/m^4) 12 p^2,   B5|m = 0,      P = m^4 - 6 m^2 p^2 + 6 p^4  (CHLPS print m^2: typo).
A functional F = sum_k int psi_k B_k + forward terms is favourable on light matter iff F[matter] >= 0 (CHLPS convention:
"F|_matter >= 0"); then light states enter with the sign of heavy states and only strengthen the bound.

  D  (NF02, Result D):  F_D = int psi2 B2 + int psi3 B3                                     (k = 2, 3 only)
  Qu (Q3-full-delta0.2): F_Q = int psi2 B2 + int psi3 B3 + a(-d_{p^2}B4|_0 + B5|_0) - c B4|_0  (flat form of Omega^Q_up)
  Ql (Omega^Q_low):      F_L = B4|_0
Check: the CHLPS F_g3 functional reproduces their printed matter polynomials exactly (2591/210 m^2 + 1/18, ...).
"""
import os, sys, json
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'prl', 'data')
OUT = os.path.join(HERE, 'cache', 'a1p_flat.txt')

def pmul(a, b):
    r = [Fr(0)]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x*y
    return r

def mom(c, j):
    return sum(Fr(ci)/(i + j + 1) for i, ci in enumerate(c))

def matter_polys(mu2, mu3, a=Fr(0), c=Fr(0), fwd_g3=Fr(0)):
    """F[matter] per |g0|^2 and per |g2|^2/m^4, as polynomials in x = m^2 (lowest power first).
    fwd_g3: coefficient of -d_{p^2}B4|_0 used by CHLPS F_g3 (their -5 d B4)."""
    # J = 0:  int psi2 (2m^2 - p^2) + int psi3 * 1
    F0 = [mu3(0) - mu2(2), 2*mu2(0)]
    # J = 2:  int psi2 (2m^2 - p^2) P + int psi3 P ;  (2x - p^2)(x^2 - 6 x p^2 + 6 p^4) = 2x^3 - 13 x^2 p^2 + 18 x p^4 - 6 p^6
    F2 = [-6*mu2(6) + 6*mu3(4), 18*mu2(4) - 6*mu3(2), -13*mu2(2) + mu3(0), 2*mu2(0)]
    # forward pieces: -a d_{p^2} B4|0 = +12 a (per |g2|^2/m^4); B5 = 0; B4|0 = 0
    F2[0] += 12*a + 12*fwd_g3
    return F0, F2

def peval(c, x):
    return sum(float(ci)*x**i for i, ci in enumerate(c))

def scan(F, ms):
    return [(m, peval(F, m*m)) for m in ms]

def main():
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    # ---------------- CHLPS F_g3 check (their eq. F KK): psi2 = (1-p)^3 p (65+155p+47p^2), psi3 = 2(1-p)^3(81-1927p+4601p^2), -5 d B4
    om3 = [Fr(1), Fr(-3), Fr(3), Fr(-1)]
    c2 = pmul(om3, [0, 65, 155, 47]); c3 = pmul(om3, [162, -3854, 9202])
    F0, F2 = matter_polys(lambda j: mom(c2, j), lambda j: mom(c3, j), fwd_g3=Fr(5))
    pr('CHLPS F_g3 check (units M = 1): per |g0|^2: ' + ' + '.join(f'({v}) m^{2*i}' for i, v in enumerate(F0)))
    pr('   CHLPS print: 2591/210 m^2 + 1/18')
    pr('  per |g2|^2/m^4: ' + ' + '.join(f'({v}) m^{2*i}' for i, v in enumerate(F2)))
    pr('   CHLPS print: 2591/210 m^6 - 239/18 m^4 - 64573/1540 m^2 + 330151/4004')
    # ---------------- NF02 (Result D)
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    s2 = [Fr(x) for x in d['psi2_exact']]; s3 = [Fr(x) for x in d['psi3_exact']]
    k0 = -2*mom(s2, 6) + 2*mom(s3, 4)
    mu2 = lambda j: mom(s2, j)/k0; mu3 = lambda j: mom(s3, j)/k0
    pr(f'\nNF02: k(0) = -2 mu2(6) + 2 mu3(4) = {float(-2*mu2(6) + 2*mu3(4)):.12f} (normalised); k0 = {float(k0):.8f}; '
       f'L_A(flat) = int psi2/p^2 diverges -> 16 L_A is the graviton term')
    F0, F2 = matter_polys(mu2, mu3)
    pr('  F_D per |g0|^2      : ' + ', '.join(f'{float(v):+.6f} m^{2*i}' for i, v in enumerate(F0)))
    pr('  F_D per |g2|^2/m^4  : ' + ', '.join(f'{float(v):+.6f} m^{2*i}' for i, v in enumerate(F2)) + '   [m -> 0 limit = 3 k(0) = 3]')
    ms = [i/200 for i in range(0, 201)]
    for nm, P in (('J=0', F0), ('J=2', F2)):
        vals = scan(P, ms); mn = min(vals, key=lambda t: t[1])
        pr(f'  F_D {nm}: min over m in [0,1] = {mn[1]:+.6f} at m = {mn[0]:.3f}; values at m = 0, .25, .5, .75, 1: '
           + ', '.join(f'{peval(P, m*m):+.5f}' for m in (0, .25, .5, .75, 1)))
    # exact sign certificate on [0,1] for x = m^2 in [0,1]: Bernstein coefficients
    def bern(c, n=None):
        n = n or len(c) - 1; c = list(c) + [Fr(0)]*(n + 1 - len(c))
        from math import comb
        return [sum(Fr(comb(k, i), comb(n, i))*c[i] for i in range(k + 1)) for k in range(n + 1)]
    for nm, P in (('J=0', F0), ('J=2', F2)):
        b = bern(P, 12)
        pr(f'  F_D {nm}: min Bernstein coefficient on x = m^2 in [0,1] (degree-12 elevation) = {float(min(b)):+.6f} '
           f'-> {"POSITIVE on all m in [0,1] [P]" if min(b) > 0 else "not certified by Bernstein"}')
    # ---------------- Q3-full-delta0.2 (Result Q upper)
    q = json.load(open(os.path.join(DATA, 'q3_full_delta0.2.json')))
    q2 = [Fr(float(x)) for x in q['psi2']]; q3 = [Fr(float(x)) for x in q['psi3']]
    a = Fr(float(q['a'])); cc = Fr(float(q['c'])); kap = float(q['kappa'])
    mu2 = lambda j: mom(q2, j); mu3 = lambda j: mom(q3, j)
    KQ = -2*mu2(6) + 2*mu3(4) + 4*a
    pr(f'\nQ3-full-delta0.2: a = {float(a):.6f}, c = {float(cc)}, kappa (file) = {kap:.6f}; '
       f'K_Q := -2 mu2(6) + 2 mu3(4) + 4a = {float(KQ):.6f}, K_Q/8 = {float(KQ)/8:.6f} (= kappa?)')
    G0, G2 = matter_polys(mu2, mu3, a=a)
    G0k, G2k = matter_polys(mu2, mu3)            # k = 2, 3 parts alone (what the light state would give without drop-out at k = 4, 5)
    pr('  F_Q per |g0|^2      : ' + ', '.join(f'{float(v):+.6f} m^{2*i}' for i, v in enumerate(G0)) + '   (k = 4, 5 parts: 0, drop-out n = 3 < 4)')
    pr('  F_Q per |g2|^2/m^4  : ' + ', '.join(f'{float(v):+.6f} m^{2*i}' for i, v in enumerate(G2)) + f'   [m -> 0: 3 K_Q = 24 kappa = {3*float(KQ):.6f}]')
    for nm, P in (('J=0', G0), ('J=2', G2)):
        vals = scan(P, ms); mn = min(vals, key=lambda t: t[1])
        b = bern(P, 12)
        pr(f'  F_Q {nm}: min over m in [0,1] = {mn[1]:+.6f} at m = {mn[0]:.3f}; at m = 0, .25, .5, .75, 1: '
           + ', '.join(f'{peval(P, m*m):+.5f}' for m in (0, .25, .5, .75, 1))
           + f'; min Bernstein = {float(min(b)):+.5f}')
    # contrast: Q's heavy-state formula continued below the gap (i.e. NO drop-out at k = 4, 5), J = 0 ++, per-state m^4 normalisation:
    # FF(u) = 2 mu2(0) + u (mu3(0) - mu2(2)) - 2c u^2  (the a-term vanishes for J = 0 ++); root in u gives the mass below which it is < 0
    import math
    A0, A1, A2 = float(2*mu2(0)), float(mu3(0) - mu2(2)), float(-2*cc)
    ur = (-A1 - math.sqrt(A1*A1 - 4*A2*A0))/(2*A2)
    pr(f'  [contrast] without the k = 4, 5 drop-out the J = 0 light action would be the heavy formula FF(u = 1/m^2) = '
       f'{A0:.4f} + {A1:.4f} u - {-A2:.1f} u^2, negative for m < {1/math.sqrt(ur):.4f} M (and ~ -u^2 as m -> 0); '
       f'with the drop-out it is negative only for m < 0.4602 M and bounded (-0.430 per |g0|^2)')
    # heavy FF (per state, m^4 normalisation) for J = 0, 2 ++ at u = 1/m^2, m < 1: shows the forward terms would dominate
    import numpy as np
    for J in (0, 2):
        for m in (0.25, 0.5, 0.75):
            u = 1/m**2
            ps = np.linspace(0, 1, 20001); y = u*ps**2
            D = np.ones_like(y) if J == 0 else 1 - 6*y + 6*y*y
            p2 = np.polynomial.polynomial.polyval(ps, [float(x) for x in q2]); p3 = np.polynomial.polynomial.polyval(ps, [float(x) for x in q3])
            FF = np.trapezoid(p2*(2 - y)*D, ps) + u*np.trapezoid(p3*D, ps) + float(a)*(2*J*(J + 1) - 1 + 1)*u**3 - 2*float(cc)*u**2
            pr(f'  [contrast] Q heavy FF formula continued to m = {m} (J = {J}): {FF:+.4f}  (true light value per normalised coupling uses the drop-out)')
    # ---------------- Omega^Q_low
    pr('\nOmega^Q_low = B4|_0: light J = 0 (n = 3 < 4) -> 0; light J = 2 (n = 4 = k) -> arc = 12 |g2|^2 p^2/m^4 at p = 0 -> 0. '
       'Both drop out of the flat value piece; in AdS the J = 2 arc is a degree-1 polynomial in nu^2, i.e. an O(Lambda^-2) term like g5.')
    with open(OUT, 'w') as fh:
        fh.write('\n'.join(L) + '\n')

if __name__ == '__main__':
    main()
