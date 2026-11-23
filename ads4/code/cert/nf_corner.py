"""NF02 flat certificate, the backward +- corner (CERT_CLOSE.md §1.3):   +-,  J >= J1,  3/8 <= u <= 1  (M <= m <= 1.633 M).

Claim certified: for each of the four sign patterns g of nf_cert/rcert (margin DELTA, no forward term),
    FF_g(J,u) = int_0^1 w(p) y(p) dp > 0,   w = W0 + u W1,   y = P^{(0,8)}_{J-4}(1 - 2 u p^2),
for all integers J >= J1 and all u in [3/8, 1].  (FF_g >= 0 for all four g  <=>  FF >= DELTA N, N = |B2| + |B3|.)

Method (all constants in FLINT/Arb ball arithmetic or exact rationals):
  split w = w_A + w_B with a C^6 cutoff chi (chi = 1 on [0,a], 1 - I_t(7,7) on [a,c], 0 on [c,1]); w_A = w chi, w_B = w (1 - chi).
  (A) Lemma C + Lemma M (FLAT_GRAVITON_BOUND.md §2-3) for w_A (1-up^2)^-8, which vanishes on [c,1] (so it is uniformly C^4):
        b^4 I_A >= gamma(u) b/(1+ds)^3 - C_R/(1-ds)^4,   gamma = -(W0''(0) + u W1''(0))/2,  ds = 8/(2J+1),  b = (2J+1) sqrt(u).
  (B) Lemma B (integration by parts with the Jacobi operator; proof in CERT_CLOSE.md §1.3): with
        L f := d/dp[(f' - f/p)(1 - u p^2) + 16 u p f]      (L p^m = (m-1)^2 p^(m-2) + (17-m)(m+1) u p^m),
      (K y')' = -4 u lam rho y,  K = p(1-up^2)^9,  rho = p(1-up^2)^8,  lam = (J-4)(J+5),  4 u lam = b^2 - 81 u,  and
        int_0^1 f y = -(4 u lam)^-1 [ beta_f + int_0^1 (L f) y ],
        beta_f = f(1)(1-u) y'(1) - ((f'(1) - f(1))(1-u) + 16 u f(1)) y(1),     for f in C^2 with f = 0 near p = 0.
      Applied three times to f0 = w_B (f0, f1 = L f0 vanish to second order at p = 1, f2 = L^2 f0 vanishes at p = 1, all exact):
        |I_B| <= (b^2 - 81u)^-3 [ |f2'(1)| (1-u) |y(1)| + int_0^1 |L^3 w_B| |y| ],
      with |y| <= min(B_n, (1-up^2)^-4)  (B_n = binom(J+4, 8) = max|P^{(0,8)}_{J-4}|;  |d^J_44| <= 1).
  Final inequality, checked at J = J1 and monotone in J (CERT_CLOSE.md):
        gamma b/(1+ds)^3 - C_R/(1-ds)^4 - b^4 (b^2-81u)^-3 [T_0(u) + T_1(u) B_n^(1/4)] > 0,  B_n^(1/4) <= (J+4)^2 / 8!^(1/4).

Corner (u >= 1 - E0, p >= 1 - D0): L^3 w_B = p^im G(delta, eps) with G an exact polynomial in delta = 1-p, eps = 1-u whose terms all have
degree >= 2; the quadratic part is integrated against min(B_n, (eps + kappa delta)^-4) exactly by scaling (sup over eps of
G(e) in Arb), the rest by K3 (delta+eps)^3 (CERT_CLOSE.md §1.3).
Usage: python3 nf_corner.py DELTA J1 [a c]        (defaults a = 0, c = 3/5; chi = 1 - I_t(7,7) on [a, c])
Read-only imports: ads4/code/robust/{nf_cert,rcert}.py (nf_cert patches rcert to centred interval evaluation on import).
"""
import sys, os, time, math
from fractions import Fraction as Fr
import flint
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
import nf_cert, rcert                                    # noqa: E402
from rcert import A, fq, Jet, gens                       # noqa: E402

PREC = 256
flint.ctx.prec = PREC
SQ = lambda x: fq(x).sqrt()


# ------------------------------------------------------------------ exact bivariate (Laurent in p) polynomials
class P2:
    """sum c[(i, j)] p^i u^j, i in Z, j >= 0, exact Fractions."""
    def __init__(self, c=None):
        self.c = {k: v for k, v in (c or {}).items() if v != 0}

    @staticmethod
    def from_p(coeffs, uj=0):
        return P2({(i, uj): Fr(v) for i, v in enumerate(coeffs)})

    def __add__(self, o):
        r = dict(self.c)
        for k, v in o.c.items(): r[k] = r.get(k, 0) + v
        return P2(r)

    def __sub__(self, o): return self + o.scale(-1)

    def scale(self, s): return P2({k: v * s for k, v in self.c.items()})

    def __mul__(self, o):
        r = {}
        for (i, j), v in self.c.items():
            for (k, l), w in o.c.items():
                r[(i + k, j + l)] = r.get((i + k, j + l), 0) + v * w
        return P2(r)

    def dp(self):
        return P2({(i - 1, j): v * i for (i, j), v in self.c.items() if i != 0})

    def L(self):
        r = {}
        for (m, j), v in self.c.items():
            a, b = (m - 1) ** 2, (17 - m) * (m + 1)
            if a: r[(m - 2, j)] = r.get((m - 2, j), 0) + v * a
            if b: r[(m, j + 1)] = r.get((m, j + 1), 0) + v * b
        return P2(r)

    def at_p(self, p0):
        """polynomial in u (list) at p = p0 (exact)."""
        out = {}
        for (i, j), v in self.c.items(): out[j] = out.get(j, 0) + v * Fr(p0) ** i
        n = max(out) + 1 if out else 1
        return [out.get(j, Fr(0)) for j in range(n)]

    def imin(self): return min([i for i, _ in self.c] + [0])

    def ushift_ppolys(self, um):
        """G~ = p^-imin G; returns list over k of fmpq_poly in p of the coefficient of (u - um)^k."""
        im = self.imin(); dj = {}
        for (i, j), v in self.c.items():
            for k in range(j + 1):
                key = (k, i - im)
                dj[key] = dj.get(key, 0) + v * math.comb(j, k) * Fr(um) ** (j - k)
        K = max([k for k, _ in dj] + [0]); D = max([i for _, i in dj] + [0])
        out = []
        for k in range(K + 1):
            cl = [dj.get((k, i), Fr(0)) for i in range(D + 1)]
            out.append(flint.fmpq_poly([flint.fmpq(x.numerator, x.denominator) for x in cl]))
        return out, im


def box_abs(ushifted, im, p0, p1, hu):
    """upper bound for |G| on [p0,p1] x [um-hu, um+hu] (centred form in both variables); p0 > 0."""
    polys = ushifted
    pm = (Fr(p0) + Fr(p1)) / 2; hp = (Fr(p1) - Fr(p0)) / 2
    T = (fq(-hp)).union(fq(hp))
    S = (fq(-hu)).union(fq(hu)) if hu else A(0)
    sh = flint.fmpq_poly([flint.fmpq(pm.numerator, pm.denominator), 1])
    tot = A(0)
    for k in reversed(range(len(polys))):
        Q = polys[k](sh) if polys[k].degree() >= 0 else polys[k]
        r = A(0)
        for c in reversed(Q.coeffs()):
            r = r * T + A(c)
        tot = tot * S + r
    v = A(tot.abs_upper())
    if im < 0: v = v * fq(p0) ** im
    return A(v.upper())


# ------------------------------------------------------------------ exact data
def smoothstep(a, c, k=7):
    """S(p) = I_t(k,k), t = (p-a)/(c-a) (exact polynomial in p); S in C^(k-1) across a and c when glued to 0 / 1."""
    t = [Fr(-a) / (c - a), Fr(1) / (c - a)]
    one_m_t = [1 - t[0], -t[1]]
    def ppow(q, n):
        r = [Fr(1)]
        for _ in range(n): r = rcert.pmul(r, q)
        return r
    S = [Fr(0)]
    for j in range(k, 2 * k):
        S = rcert.padd(S, rcert.pscale(rcert.pmul(ppow(t, j), ppow(one_m_t, 2 * k - 1 - j)), math.comb(2 * k - 1, j)))
    return S


def peval_exact(poly, x):
    r = Fr(0)
    for c in reversed(poly): r = r * x + c
    return r


def pderiv_exact(poly, k):
    for _ in range(k): poly = rcert.pder(poly)
    return poly


# ------------------------------------------------------------------ (A): Lemma M constant of w_A (1-up^2)^-8 on [0, c]
def CR_A(gW0, gW1, chi_pieces, ulo, uhi, NP):
    """chi_pieces: list of (lo, hi, chi_poly) covering [0, c]; cells of width 1/NP (lo, hi multiples of 1/NP).
    Returns C_R of FLAT_GRAVITON_BOUND §3 (M) for w = (W0 + u W1) chi (1-up^2)^-8, sup over u in [ulo, uhi]."""
    u = fq(ulo).union(fq(uhi)); uu = fq(uhi)
    n_w4 = A(0); n_pw3 = A(0); n_w2 = A(0); s2 = A(0); s1 = A(0); width = fq(Fr(1, NP))
    for lo, hi, chi in chi_pieces:
        P0 = rcert.pmul(gW0, chi); P1 = rcert.pmul(gW1, chi)
        i0, i1 = int(Fr(lo) * NP), int(Fr(hi) * NP)
        assert Fr(i0, NP) == lo and Fr(i1, NP) == hi
        for i in range(i0, i1):
            p0 = rcert.iv(Fr(i, NP), Fr(i + 1, NP))
            pj = Jet([p0, A(1)])
            w = rcert.poly_jet(P0, pj) + rcert.poly_jet(P1, pj) * u
            q = (1 - u * pj * pj).inv(); q2 = q * q; q4 = q2 * q2; w = w * (q4 * q4)
            wd = [w.der(k) for k in range(5)]
            n_w4 += A(wd[4].abs_upper()) * width
            n_pw3 += A((p0 * wd[3]).abs_upper()) * width
            n_w2 += A(wd[2].abs_upper()) * width
            s2 = s2.max(A(wd[2].abs_upper())); s1 = s1.max(A(wd[1].abs_upper()))
    CR = 2 / A.pi() * (A(2) / 3 * n_w4 + 6 * uu * n_pw3 + 4 * uu * n_w2) + uu * uu.sqrt() * (3 * A.pi() / 4 * s2 + s1)
    return A(CR.upper())


# ------------------------------------------------------------------ (B): integrals of |L^3 w_B| against the envelope
def pcells(lo, hi, nuni, dmin, ratio):
    """cells: uniform on [lo, 1 - d0] then geometric towards 1 (widths ~ (1-ratio) * delta), last cell [1 - dmin, 1]."""
    cells = []; hi = Fr(hi)
    if hi < 1:
        for i in range(nuni): cells.append((lo + (hi - lo) * i / nuni, lo + (hi - lo) * (i + 1) / nuni))
        return cells
    d0 = Fr(1, 20)
    if 1 - d0 <= lo: d0 = 1 - lo
    else:
        for i in range(nuni): cells.append((lo + (1 - d0 - lo) * i / nuni, lo + (1 - d0 - lo) * (i + 1) / nuni))
    d = d0
    while d > dmin:
        dn = max(d * ratio, dmin); dn = Fr(dn).limit_denominator(10 ** 12) if dn != dmin else dmin
        cells.append((1 - d, 1 - dn)); d = dn
    cells.append((1 - dmin, Fr(1)))
    return cells


def T_far(f3pieces, u0, u1, pmax_corner, pc_args):
    """sup over u in [u0,u1] of int over the pieces (restricted to p <= pmax_corner) of |f3| (1-up^2)^-4 (cells)."""
    um = (Fr(u0) + Fr(u1)) / 2; hu = (Fr(u1) - Fr(u0)) / 2
    tot = A(0)
    for lo, hi, f3 in f3pieces:
        hi_eff = min(Fr(hi), Fr(pmax_corner))
        if hi_eff <= lo: continue
        sh, im = f3.ushift_ppolys(um)
        for p0, p1 in pcells(Fr(lo), hi_eff, *pc_args):
            qmin = 1 - fq(u1) * fq(p1) ** 2
            assert bool(qmin > 0), 'q must stay positive away from the corner'
            tot += box_abs(sh, im, p0, p1, hu) * qmin ** -4 * fq(p1 - p0)
    return A(tot.upper())


def corner_K(f3_hi, d0, e0):
    """f3 on [c,1] in variables delta = 1-p, eps = 1-u: G = p^-imin f3; returns K with |f3| <= K (delta+eps)^2 on
    delta <= d0, eps <= e0 (asserts that every monomial has total degree >= 2)."""
    im = f3_hi.imin()
    G = f3_hi * P2({(-im, 0): Fr(1)})
    # substitute p = 1 - d, u = 1 - e  (exact)
    co = {}
    for (i, j), v in G.c.items():
        for a in range(i + 1):
            for b in range(j + 1):
                key = (a, b)
                co[key] = co.get(key, 0) + v * math.comb(i, a) * (-1) ** a * math.comb(j, b) * (-1) ** b
    co = {k: v for k, v in co.items() if v != 0}
    low = min(a + b for a, b in co)
    assert low >= 2, 'f3 does not vanish to second order at the corner'
    s = Fr(d0) + Fr(e0)
    K3 = sum(abs(v) * s ** (a + b - 3) for (a, b), v in co.items() if a + b >= 3)
    lead = {k: co.get(k, Fr(0)) for k in [(2, 0), (1, 1), (0, 2)]}
    pf = fq(1 - Fr(d0)) ** (-im)              # |f3| = |G| p^im <= |G| / (1-d0)^(-im)
    return lead, fq(K3), pf


def G_sup(lead, kappa, emax=Fr(20), ne=400, xmax=Fr(40), nx=1600):
    """sup_{e >= 0} G(e),  G(e) = int_0^oo |c20 x^2 + c11 x e + c02 e^2| min(1, (e + kappa x)^-4) dx  (Arb, cells in e and x;
    closed-form tails).  With e = eps B^(1/4), x = delta B^(1/4):  int_0^oo |quad(delta, eps)| min(B, (eps + kappa delta)^-4) d delta
    = B^(1/4) G(e)."""
    c20, c11, c02 = (fq(lead[k]) for k in [(2, 0), (1, 1), (0, 2)])
    a20, a11, a02 = abs(c20), abs(c11), abs(c02)
    best = A(0)
    for i in range(ne):
        e0, e1 = emax * i / ne, emax * (i + 1) / ne
        E = rcert.iv(e0, e1)
        tot = A(0)
        for k in range(nx):
            x0, x1 = xmax * k / nx, xmax * (k + 1) / nx
            X = rcert.iv(x0, x1)
            qd = A((c20 * X * X + c11 * X * E + c02 * E * E).abs_upper())
            env = A(1).min((fq(e0) + kappa * fq(x0)) ** -4) if (e0 + x0) > 0 else A(1)
            tot += qd * env * fq(x1 - x0)
        X = fq(xmax); e1a = fq(e1)
        tot += (a20 / X + a11 * e1a / (2 * X * X) + a02 * e1a * e1a / (3 * X ** 3)) / kappa ** 4
        best = best.max(A(tot.upper()))
    # e >= emax (>= 1): G(e) <= (a20/(3k^3) + a11/(6k^2) + a02/(3k))/e
    tail = (a20 / (3 * kappa ** 3) + a11 / (6 * kappa ** 2) + a02 / (3 * kappa)) / fq(emax)
    return A(best.max(tail).upper())


def main():
    delta = Fr(sys.argv[1]); J1 = int(sys.argv[2])
    a = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); c = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(3, 5)
    NP = 2000; E0 = Fr(1, 100); D0 = Fr(1, 100)
    t00 = time.time()
    F = nf_cert.load('NF02')
    psi2, psi3, _ = F
    S = smoothstep(a, c, 7)
    one_m_S = rcert.padd([Fr(1)], rcert.pscale(S, -1))
    # C^6 gluing of the cutoff (exact)
    for k in range(7):
        Sk = pderiv_exact(S, k)
        assert peval_exact(Sk, a) == 0 and peval_exact(Sk, c) == (1 if k == 0 else 0), 'smoothstep gluing'
    print('cutoff: a = %s, c = %s, chi = 1 - I_t(7,7) (C^6 gluing checked exactly)' % (a, c), flush=True)
    J = J1
    ds = fq(Fr(8, 2 * J + 1))
    root8f = A(40320).root(4)
    ok_all = True
    # u cells: coarse (for C_R) and fine (for T): geometric in eps = 1-u from 5/8 to E0, then the corner cell [1-E0, 1]
    eps_cuts = [Fr(5, 8)]
    while eps_cuts[-1] > E0:
        e = eps_cuts[-1] * Fr(9, 10)
        eps_cuts.append(max(Fr(e).limit_denominator(10 ** 6), E0))
    ucells = [(1 - e1, 1 - e2) for e1, e2 in zip(eps_cuts[:-1], eps_cuts[1:])] + [(1 - E0, Fr(1))]
    coarse = [Fr(3, 8), Fr(1, 2), Fr(5, 8), Fr(3, 4), Fr(7, 8), Fr(15, 16), Fr(1)]
    for g in gens(psi2, psi3, Fr(0), delta, '+-', 1):
        w = P2.from_p(g.W0) + P2.from_p(g.W1, 1)
        wB_mid = w * P2.from_p(S); wB_hi = w
        # exact structure checks
        for k in range(7):
            d_mid = wB_mid; d_hi = wB_hi
            for _ in range(k): d_mid = d_mid.dp(); d_hi = d_hi.dp()
            assert all(x == 0 for x in d_mid.at_p(a)), 'w_B not C^6 at a'
            vm, vh = d_mid.at_p(c), d_hi.at_p(c); nmx = max(len(vm), len(vh))
            assert vm + [Fr(0)] * (nmx - len(vm)) == vh + [Fr(0)] * (nmx - len(vh)), 'w_B not C^6 at c'
        f = [wB_hi]
        for _ in range(3): f.append(f[-1].L())
        assert all(x == 0 for x in f[0].at_p(1)) and all(x == 0 for x in f[0].dp().at_p(1)), 'f0(1), f0\'(1)'
        assert all(x == 0 for x in f[1].at_p(1)) and all(x == 0 for x in f[1].dp().at_p(1)), 'f1(1), f1\'(1)'
        assert all(x == 0 for x in f[2].at_p(1)), 'f2(1)'
        f3_mid = wB_mid.L().L().L(); f3_hi = f[3]
        # boundary-term polynomial  f2'(1) = (1-u)^2 r(u)  (exact division)
        bp = f[2].dp().at_p(1)
        # divide bp(u) by (1-u)^2 exactly
        q = list(bp)
        for _ in range(2):
            assert peval_exact(q, 1) == 0, "f2'(1) not divisible by (1-u)^2"
            # synthetic division by (u - 1), then sign
            out = [Fr(0)] * (len(q) - 1); acc = Fr(0)
            for i in reversed(range(1, len(q))):
                acc = acc * 1 + q[i]; out[i - 1] = acc
            q = out
        rpoly = q        # f2'(1) = (u-1)^2 rpoly(u) = (1-u)^2 rpoly(u)
        lead, K3, pf = corner_K(f3_hi, D0, E0)
        kappa = 2 * (1 - fq(E0)) * (1 - fq(D0) / 2)          # q = eps + u delta (2 - delta) >= eps + kappa delta  (u >= 1-E0, delta <= D0)
        Gs = G_sup(lead, kappa)
        Smax = fq(E0) + 2 * fq(D0)
        print('%s: gamma(u) = %.6g + %.6g u ; corner (delta,eps)^2 part %s, sup G = %s, K3 = %s ; '
              "f2'(1) = (1-u)^2 r(u), r = %s" % (g.label, float(-g.W0[2]), float(-g.W1[2]), {k: float(v) for k, v in lead.items()},
                                                 Gs.str(5), K3.str(5), [float(x) for x in rpoly]), flush=True)
        pieces = [(a, c, f3_mid), (c, Fr(1), f3_hi)]
        chi_pieces = ([(Fr(0), a, [Fr(1)])] if a > 0 else []) + [(a, c, one_m_S)]
        CRc = {}
        for lo, hi in zip(coarse[:-1], coarse[1:]):
            CRc[(lo, hi)] = CR_A(g.W0, g.W1, chi_pieces, lo, hi, NP)
        worst = None; ok_g = True
        for u0, u1 in ucells:
            lo_c = max(x for x in coarse[:-1] if x <= u0); hi_c = min(x for x in coarse[1:] if x >= u1)
            CR = max([CRc[k] for k in CRc if k[0] < u1 and k[1] > u0], key=lambda z: float(z.upper()))
            gam = min(-g.W0[2] - g.W1[2] * u0, -g.W0[2] - g.W1[2] * u1)
            assert gam > 0
            rmax = A(0)
            for ua, ub in [(u0, u1)]:
                # |r(u)| on [u0,u1] by centred evaluation
                rr = P2({(0, j): v for j, v in enumerate(rpoly)})
                sh, im = rr.ushift_ppolys((Fr(ua) + Fr(ub)) / 2)
                rmax = box_abs(sh, im, Fr(1, 2), Fr(1, 2), (Fr(ub) - Fr(ua)) / 2)
            if u1 <= 1 - E0:
                # case A: q >= eps > 0 everywhere; boundary |f2'(1)|(1-u)|y(1)| <= |r| (1-u)^3 (1-u)^-4 = |r|/(1-u)
                T0 = T_far(pieces, u0, u1, Fr(1), (200, Fr(1, 10 ** 9), Fr(9, 10))) + rmax / (1 - fq(u1))
                T1 = A(0)
            else:
                # case B: cells up to p = 1 - D0, corner bound beyond; boundary <= |r| B^(1/4)
                T0 = T_far(pieces, u0, u1, 1 - D0, (200, Fr(1, 10 ** 9), Fr(9, 10)))
                T1 = pf * Gs + rmax
                # cubic and higher corner terms: pf K3 int (delta+eps)^3 min(B, q^-4) <= pf K3 (1/4 + log(Smax B^(1/4)))/kappa
                T0 = T0 + pf * K3 * (A(1) / 4 + (Smax * A(J + 4) ** 2 / root8f).log()) / kappa
            # J = J1 (monotone in J); b at the worst u of the cell
            bmin = (2 * J + 1) * fq(u0).sqrt()
            Fac = A((2 * J + 1) ** 4) / (fq(u0) * (A((2 * J + 1) ** 2) - 81) ** 3)       # b^4 (b^2 - 81u)^-3, max at u0
            B14 = A(J + 4) ** 2 / root8f
            main_t = fq(gam) * bmin / (1 + ds) ** 3
            val = main_t - CR / (1 - ds) ** 4 - Fac * (T0 + T1 * B14)
            ok = bool(val > 0)
            ok_all &= ok; ok_g &= ok
            rel = val / main_t
            if worst is None or float(rel.mid()) < worst[0]:
                worst = (float(rel.mid()), u0, u1, CR, Fac * T0, Fac * T1 * B14, main_t)
            if not ok:
                print('  FAIL u in [%s, %s]: gamma b/(1+ds)^3 = %s, C_R = %s, B-part = %s' % (float(u0), float(u1), main_t.str(5), CR.str(5),
                      (Fac * (T0 + T1 * B14)).str(5)), flush=True)
        print('  %s: %d u-cells %s; tightest (val / main term) = %.4f on u in [%.5f, %.5f]: main %s, C_R %s, B-far %s, B-corner %s  (%.0fs)'
              % (g.label, len(ucells), 'PASS' if ok_g else 'FAIL', worst[0], float(worst[1]), float(worst[2]), worst[6].str(5),
                 worst[3].str(5), worst[4].str(5), worst[5].str(5), time.time() - t00), flush=True)
    print('VERDICT nf_corner (delta=%s, +-, J >= %d, 3/8 <= u <= 1, a=%s, c=%s): %s  (%.0fs)' % (delta, J1, a, c, 'PASS' if ok_all else 'FAIL', time.time() - t00))


if __name__ == '__main__':
    main()
