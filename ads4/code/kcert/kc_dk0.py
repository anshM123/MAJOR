"""kc_dk0.py -- certified (Arb ball arithmetic) light-side constants of Omega^EE's k = 2 piece at Delta_g = 300, 1000.

Computes, for the DEFINED functional (continuum weight chi_re of prl/data/k2/k2_continuum_D_<Dg>.json times theta(Lam+100, 10)):
  (a) K0 = -2 int psi2 p^6 + 2 int psi3 p^4 (exact rational; NF02 with k0 = 27.48509817 exactly),
  (b) dK0 := -2 int_0^oo (Re X_theta - psi) p^6 dnu/Lam,  Re X_theta = theta chi_re (R_- + R_+), p = nu/Lam   [the published definition,
      code/cert/k2_cont.py checks(); Supplement row L2]  -> K0(1 + dK0/K0) is the |g3|^2 light coefficient up to H17,
  (c) V3 := int_0^oo theta chi_re [R_-(nu) D(nu+i) + R_+(nu) D(nu-i)] dnu, D(z) = 1/(z^2 + 9/4): the realised Einstein light value
      (16 Lam V3 replaces 16 L_A), and the exact target L_A/Lam = int_0^Lam psi(nu) D(nu) dnu (closed form),
  (d) the moments M_j := int_0^oo Re X_theta p^(2j) dnu/Lam, j = 0..3 (H17 sensitivity coefficients; M_3 = -(K0 + dK0)/2 + int psi3 p^4),
      and the same with the "shifted" light reading, Ms_j := int theta chi_re 2Re[R_-(nu)(nu+i)^(2j)] dnu/Lam^(2j+1) (diagnostic only).
All integrals by acb.integral (rigorous Gauss-Legendre with Bernstein-ellipse error bounds) on [0,Lam] and [Lam, Lam+60] in chunks,
plus a rigorous tail bound beyond Lam + 60 (tail_bound()).
Usage: python3 kc_dk0.py 300|1000 [quick]   -> cache/kc_dk0_<Dg>.txt
"""
import sys, time, os
from fractions import Fraction as Fr
from flint import arb, acb, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kc_common import *  # noqa

ctx.prec = 128


def LA_over_Lam(Dg, p2):
    """exact closed form of int_0^1 psi2(p) dp/(p^2 + a^2), a^2 = 9/(4 Lam^2); returns arb (this is L_A)."""
    a2 = Fr(9, 4 * (Dg - 2) * (Dg - 1))
    c = list(p2)                       # lowest power first
    # divide sum c_n p^n by (p^2 + a2): highest-first synthetic division
    n = len(c) - 1
    q = [Fr(0)] * max(n - 1, 1)
    r = c[:]
    for k in range(n, 1, -1):
        qk = r[k]
        q[k - 2] = qk
        r[k] -= qk
        r[k - 2] -= qk * a2
    r1, r0 = r[1], r[0]
    old = ctx.prec; ctx.prec = 512
    a = arb(fq(a2)).sqrt()
    val = arb(fq(sum(qn / (i + 1) for i, qn in enumerate(q))))
    val += arb(fq(r1)) / 2 * ((1 + arb(fq(a2))) / arb(fq(a2))).log()
    val += arb(fq(r0)) / a * (1 / a).atan()
    ctx.prec = old
    return val


def main(Dg, quick=False):
    t0 = time.time()
    out = []

    def pr(s):
        print(s, flush=True); out.append(s)
    L = Lam(Dg)
    p2, p3 = nf02()
    K0 = K0_exact()
    pr('# kc_dk0: Delta_g = %d, Lam = %s' % (Dg, L.str(20)))
    pr('K0 (exact rational) = %s ~ 1 + %.6e' % (str(K0)[:60] + '...', float(K0 - 1)))
    pr('  -2 int psi2 p^6 = %.12f ; 2 int psi3 p^4 = %.12f' % (float(-2 * poly_moment(p2, 6)), float(2 * poly_moment(p3, 4))))
    for j in range(4):
        pr('  int psi2 p^%d = %.15f' % (2 * j, float(poly_moment(p2, 2 * j))))
    for j in range(3):
        pr('  int psi3 p^%d = %.15f' % (2 * j, float(poly_moment(p3, 2 * j))))
    C = Chi(Dg)
    Rm_s, Rp_s = kernels()
    S = Rat('nu*((%s) + (%s))' % (Rm_s, Rp_s))                                   # nu (R_- + R_+), no pole at 0
    Tm = Rat('nu*(%s)' % Rm_s)                                                     # nu R_-
    Tv = Rat('nu*((%s)/((nu+I)**2 + 9/4) + (%s)/((nu-I)**2 + 9/4))' % (Rm_s, Rp_s))  # nu [R_- D(nu+i) + R_+ D(nu-i)]
    pc = [fq(c) for c in p2]

    def psi(z):
        p = z / L
        s = acb(0)
        for c in reversed(pc):
            s = s * p + arb(c)
        return s

    I = acb(0, 1)

    def make(kind, j, with_psi):
        def f(z, analytic):
            if z.rad() > 4:
                return NAN
            if kind == 'V3':
                ch = C.over_z(z, analytic)
                if not ch.is_finite():
                    return NAN
                v = theta(z, L) * ch * Tv(z)
                if with_psi:
                    v -= psi(z) / (z * z + arb(9) / 4)
                return v
            if kind == 'X':                        # Re X_theta nu^(2j) ; computed as theta chi/z * (nu S) * nu^(2j)... (chi/z)(z S)/z? use S
                ch = C(z, analytic) if j >= 1 else C.over_z(z, analytic)
                if not ch.is_finite():
                    return NAN
                if j >= 1:
                    v = theta(z, L) * ch * S(z) * z ** (2 * j - 1)
                else:
                    v = theta(z, L) * ch * S(z)
                if with_psi:
                    v -= psi(z) * z ** (2 * j)
                return v
            if kind == 'Xs':                       # theta chi 2Re[R_- (z+i)^(2j)] continued analytically: chi [R_-(z)(z+i)^2j + R_+(z)(z-i)^2j]
                ch = C.over_z(z, analytic)
                if not ch.is_finite():
                    return NAN
                Tp = Rat('nu*(%s)' % Rp_s) if not hasattr(make, 'Tp') else make.Tp
                make.Tp = Tp
                v = theta(z, L) * ch * (Tm(z) * (z + I) ** (2 * j) + Tp(z) * (z - I) ** (2 * j))
                if with_psi:
                    v -= psi(z) * z ** (2 * j)
                return v
        return f

    def integrate(f, a, b, step, tol):
        tot = acb(0); x = a
        while x < b:
            y = min(x + step, b)
            tot += acb.integral(f, arb(x) if not isinstance(x, arb) else x, arb(y) if not isinstance(y, arb) else y,
                                abs_tol=tol, eval_limit=10 ** 7, deg_limit=200)
            x = y
        return tot

    Lf = float(L)
    # breakpoint exactly at Lam (psi has compact support): integrate [0, Lam] with psi, [Lam, Lam+60] without
    def two_piece(kind, j, tol, step=10.0):
        f1 = make(kind, j, True); f2 = make(kind, j, False)
        A = acb(0); x = 0.0
        edges = [0.0]
        while edges[-1] + step < Lf:
            edges.append(edges[-1] + step)
        for i in range(len(edges)):
            lo = arb(edges[i]); hi = arb(edges[i + 1]) if i + 1 < len(edges) else L
            A += acb.integral(f1, lo, hi, abs_tol=tol, eval_limit=10 ** 7, deg_limit=200)
        Bp = acb(0); e = L
        for k in range(6):
            Bp += acb.integral(f2, e + 10 * k, e + 10 * (k + 1), abs_tol=tol, eval_limit=10 ** 7, deg_limit=200)
        return A + Bp

    # tail beyond Lam + 60: for nu >= T, |R_-+| <= 18 nu^2 (each linear factor |2nu + ai| <= 2nu(1 + 11/(2T)) <= 2.04 nu, T >= 358),
    # |Q| <= 1/nu^2, |B| <= 2 sum|x| exp(-(nu - c_max)^2/(2 0.9^2)), |theta| <= 1, sing terms <= sum|y| nu^4 (nu^2+1) log(nu^2+1) e^{-nu^2/2}.
    # So |theta chi (R_-+R_+) nu^m| <= 36 nu^m [2 sum|x| e^{-(nu-c)^2/1.62} + sum|y| 2 nu^6 log(2nu^2) e^{-nu^2/2}], m <= 6, and
    # int_T^oo nu^6 e^{-(nu-c)^2/1.62} <= (2T)^6 e^{-(T-c)^2/1.62} * 1.62/(2(T-c)) + int_{2T}^oo (2(nu-c))^6 e^{..} (negligible), bound below.
    def tail_bound(m):
        T = L + 60; c = arb(float(C.c.max()))
        assert float(T - c) >= 29
        sx = arb(float(np.abs(C.x).sum())); sy = arb(float(sum(abs(v) for v in C.y)))
        # crude: nu^m e^{-(nu-c)^2/1.62} <= (T)^m e^{6 (nu-T)/T} e^{-(nu-c)^2/1.62} and (nu-c)^2 >= (T-c)^2 + 2(T-c)(nu-T)
        rate = 2 * (T - c) / arb('1.62') - arb(m) / T
        g = T ** m * (-(T - c) ** 2 / arb('1.62')).exp() / rate
        # sing part: |sum_j y_j nu^2j l e^{-nu^2/2}| <= sy nu^4 * 2 nu^3 e^{-nu^2/2} (l <= 2nu^2 * nu), times 36 nu^2 nu^m -> 72 sy nu^(m+9) e^{..};
        # int_T^oo nu^a e^{-nu^2/2} <= T^a e^{-T^2/2}/(T - a/T)  (log-derivative a/nu - nu <= a/T - T)
        hsing = 2 * T ** (m + 9) * (-(T * T) / 2).exp() / (T - arb(m + 9) / T)
        return 36 * (2 * sx * g + sy * hsing)

    Lf7 = L ** 7
    res = {}
    tol_scale = {300: 1e-3, 1000: 1e-2}[Dg]
    # (b) dK0
    t1 = time.time()
    J6 = two_piece('X', 3, arb(tol_scale))
    tb = tail_bound(6)
    J6 = J6 + acb(arb(0, tb))
    dK0 = -2 * J6.real / Lf7
    res['dK0'] = dK0
    pr('dK0 = -2 int (Re X_theta - psi) p^6 dnu/Lam = %s   [%.0fs]' % (dK0.str(15, radius=True), time.time() - t1))
    pr('  Lam^2 dK0 = %s' % (dK0 * L * L).str(10, radius=True))
    pr('  K = K0 + dK0 (up to H17) = %s' % (arb(fq(K0)) + dK0).str(15, radius=True))
    if quick:
        return out, res
    # (c) V3
    t1 = time.time()
    V3d = two_piece('V3', 0, arb(1e-18))
    V3d = V3d + acb(arb(0, tail_bound(0)))
    LA = LA_over_Lam(Dg, p2)
    V3 = V3d.real + LA / L
    pr('V3 (realised Einstein light value) = %s ; L_A/Lam (exact target) = %s ; V3 - L_A/Lam = %s  [%.0fs]' % (
        V3.str(18, radius=True), (LA / L).str(18), V3d.real.str(5, radius=True), time.time() - t1))
    pr('  16 L_A = %s ; 16 Lam V3 = %s' % ((16 * LA).str(15, radius=True), (16 * L * V3).str(15, radius=True)))
    res['V3'] = V3; res['LA'] = LA
    # (d) moments
    for j in range(3):
        t1 = time.time()
        Jj = two_piece('X', j, arb(tol_scale) / L ** (6 - 2 * j))
        Jj = Jj + acb(arb(0, tail_bound(2 * j)))
        mj = Jj.real / L ** (2 * j + 1)
        exact = arb(fq(poly_moment(p2, 2 * j)))
        pr('M_%d - int psi2 p^%d = int (Re X_theta - psi) p^%d dnu/Lam = %s ; M_%d = %s  [%.0fs]' % (
            j, 2 * j, 2 * j, mj.str(12, radius=True), j, (mj + exact).str(15, radius=True), time.time() - t1))
        res['M%d' % j] = mj + exact
    res['M3'] = -dK0 / 2 + arb(fq(poly_moment(p2, 6)))
    pr('M_3 = %s' % res['M3'].str(15, radius=True))
    for j in (3,):
        t1 = time.time()
        Js = two_piece('Xs', j, arb(tol_scale))
        Js = Js + acb(arb(0, tail_bound(2 * j) * 2 ** 6))
        dks = -2 * Js.real / Lf7
        res['dK0s'] = dks
        pr('diagnostic: shifted reading dK0_s = -2 int (theta chi 2Re[R_-(nu+i)^6] - psi nu^6)/Lam^7 = %s ; Lam^2 dK0_s = %s [%.0fs]' % (
            dks.str(12, radius=True), (dks * L * L).str(8, radius=True), time.time() - t1))
    pr('total %.0fs' % (time.time() - t0))
    res['K0'] = arb(fq(K0)); res['Lam'] = L
    return out, res


if __name__ == '__main__':
    Dg = int(sys.argv[1]); quick = len(sys.argv) > 2
    o, _ = main(Dg, quick)
    if not quick:
        open(os.path.join(HERE, 'cache', 'kc_dk0_%d.txt' % Dg), 'w').write('\n'.join(o) + '\n')
