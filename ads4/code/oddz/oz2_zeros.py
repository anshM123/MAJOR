"""PARITY_ODD_DZ §3: exact order of the zero of every detector amplitude (all components m, source helicity +2 and -2,
zeroth moment E and first moments E^(1) [midpoint- and origin-centred]) at the N = infinity [TT] dimensions, for ALL n,
at fixed J (symbolic Delta).  M_m(s) = G_K(s) P_m(s) (pipeline.amplitude), so

   ord_{s0} M_m = #{Gamma factors of G_K with a zero at s0} + mult_{s0} P_m,

and for n beyond the largest rational root of P_m only the e13 factor 1/Gamma(e13 + Kb) vanishes (simply).  Hence the
statement 'simple zero for every n >= 0' is decided exactly by finitely many checks per (J, m).
Points: odd structures (S^odd even J >= 4, S^oddJ odd J >= 5, J = 0, 2 odd): s0 = (Delta0 + J)/2 with
        Delta0 = J+3+2n (J >= 4), 7+2n (J = 0, 2);
        controls: even structures (S^++, S^+-; J = 0, 2 even) at Delta0 = J+2+2n (J >= 4), 6+2n (J = 0, 2);
                  odd structures at the even-twist points (must NOT vanish identically) and at h = 1/2.
Also: structure regularity (the coefficient vector has no common zero at s0), and for J >= 4 the identity
   M_m = [KERNELS_D3 closed form of the ladder head] x [ladder factor]   as an exact identity in Delta (gammasimp).
usage: python3 -B oz2_zeros.py name J [J ...]   (name in odd, oddJ, pp, pm, even)   -> cache/oz2_<name>_J<J>.txt"""
import sys, time
from common import *
from oz_amps import amps

HALF = sp.Rational(1, 2)

def gamma_args(parity, K):
    Ka, Kb, KB, Kg = K
    h = HALF if parity == 'odd' else 0
    e13, e14, e34 = 5 - S + h, S + h, S + h
    return [e14 + Ka, e13 + Kb, e34 + KB, e34 - HALF + Kg]

def G_expr(parity, K):
    h = HALF if parity == 'odd' else 0
    e14, e34 = S + h, S + h
    pre = -4*sp.pi*sp.I*sp.exp(sp.I*sp.pi*e14)*sp.pi**sp.Rational(5, 2)*sp.Integer(2)**(4 - 2*e34)
    return pre/sp.Mul(*[sp.gamma(a) for a in gamma_args(parity, K)])

def n_gamma_zeros(parity, K, s0):
    return sum(1 for a in gamma_args(parity, K) if (a.subs(S, s0)).is_integer and a.subs(S, s0) <= 0)

def rational_roots(P):
    if P == 0:
        return None
    Pp = sp.Poly(P, S)
    out = {}
    for f, mlt in sp.factor_list(Pp)[1]:
        if f.degree() == 1:
            c1, c0 = f.all_coeffs()
            out[sp.Rational(-c0, c1)] = mlt
    return out

def orders(P, parity, K, s_of_n, nmax_extra=3):
    """returns (list of orders for n = 0..N, N, tail order) ; the tail order holds for every n > N."""
    roots = rational_roots(P)
    if roots is None:
        return None
    big = [r for r in roots if (r - s_of_n(0)).is_integer and r >= s_of_n(0)]
    # largest n where something special could happen: last root, or where Gamma-zero pattern stabilises
    nG = 0
    while True:          # Gamma pattern is monotone: find n after which it is constant
        if all(n_gamma_zeros(parity, K, s_of_n(q)) == n_gamma_zeros(parity, K, s_of_n(nG)) for q in range(nG, nG + 60)):
            break
        nG += 1
    N = max([nG] + [int(r - s_of_n(0)) for r in big]) + nmax_extra
    ords = [n_gamma_zeros(parity, K, s_of_n(q)) + roots.get(s_of_n(q), 0) for q in range(N + 1)]
    tail = n_gamma_zeros(parity, K, s_of_n(N + 1))
    return ords, N, tail

def closed_head(name, J, m):
    """KERNELS_D3 §5.3 / ANGULAR_PROOFS closed form x ladder factor, source +2, theta = 0 (None if not covered)."""
    Dl = 2*S - J
    if name == 'odd' and J >= 4:
        if m > 0:
            return 0
        k = -m
        M0 = 3*(-1)**(J//2)*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(4 - J)*sp.rf(Dl - J - 1, J) \
            / (sp.gamma(S + HALF)*sp.gamma(Dl + 3)*sp.gamma(J + sp.Rational(3, 2) - S))
        return M0*sp.binomial(J, k)*sp.rf(4, k)/sp.rf(Dl - 1 - k, k)
    if name == 'oddJ':
        if m < 4:
            return 0
        k = m - 4
        M4 = -3*sp.I*(-1)**((J - 5)//2)*sp.exp(sp.I*sp.pi*S)*sp.pi**3*sp.Integer(2)**(4 - J)*sp.rf(Dl - J - 1, J - 4) \
            / (sp.gamma(S + HALF)*sp.gamma(Dl - 1)*sp.gamma(J + sp.Rational(3, 2) - S))
        return M4*sp.binomial(J - 4, k)*sp.rf(4, k)/sp.rf(Dl - 5 - k, k)
    if name == 'odd' and J == 2:
        if m > 0:
            return 0
        k = -m
        M0 = sp.Rational(3, 4)*sp.pi**sp.Rational(7, 2)*sp.exp(sp.I*sp.pi*S)*sp.Integer(2)**(1 - Dl)*(Dl - 3)*(Dl - 2) \
            / (sp.gamma(2 + Dl/2)*sp.gamma((Dl + 3)/2)**2*sp.gamma((7 - Dl)/2))
        return M0*sp.binomial(J, k)*sp.rf(4, k)/sp.rf(Dl - 1 - k, k)
    return None

def run(name, J):
    parity = 'even' if name in ('pp', 'pm', 'even') else 'odd'
    out = []
    pr = lambda x: (print(x, flush=True), out.append(str(x)))
    t0 = time.time()
    r = amps(name, J)
    pr(f'=== {name}  J = {J}  ({parity}; {len(r["monos"])} monomials)  [{time.time()-t0:.0f}s]')
    # structure regularity: gcd of the coefficient vector
    import functools
    g = sp.Poly(functools.reduce(sp.gcd, [sp.sympify(c) for c in r['coeffs'] if c != 0]), S)
    pr(f'  gcd of structure coefficients (as polynomials in s): {g.as_expr()}')
    if parity == 'odd':
        s_dz = (lambda n: sp.Rational(2*J + 3, 2) + n) if J >= 4 else (lambda n: sp.Rational(7 + J, 2) + n)
        s_ctrl = (lambda n: sp.Integer(J + 1) + n) if J >= 4 else (lambda n: sp.Integer(3 + J) + n)
        lab, labc = ('Delta = J+3+2n' if J >= 4 else 'Delta = 7+2n'), ('Delta = J+2+2n' if J >= 4 else 'Delta = 6+2n')
    else:
        s_dz = (lambda n: sp.Integer(J + 1) + n) if J >= 4 else (lambda n: sp.Integer(3 + J) + n)
        s_ctrl = (lambda n: sp.Rational(2*J + 3, 2) + n) if J >= 4 else (lambda n: sp.Rational(7 + J, 2) + n)
        lab, labc = ('Delta = J+2+2n' if J >= 4 else 'Delta = 6+2n'), ('Delta = J+3+2n' if J >= 4 else 'Delta = 7+2n')
    regular = all(g.as_expr().subs(S, s_dz(n)) != 0 for n in range(0, 200))
    pr(f'  structure non-degenerate at every {lab}: {regular}  (gcd has no root of that form)')
    allok = regular
    for key in [k for k in r if isinstance(k, tuple)]:
        K, P = r[key]
        summary = []
        ok = True
        for m in sorted(P):
            res = orders(P[m], parity, K, s_dz)
            ords, N, tail = res
            if key[0] == 'zero':      # zeroth moment: exactly simple
                good = all(o == 1 for o in ords) and tail == 1
            else:                     # first moments: at least simple (higher order only helps)
                good = all(o >= 1 for o in ords) and tail >= 1
                if good and not (all(o == 1 for o in ords) and tail == 1):
                    summary.append(f'm={m}: orders {ords} tail {tail} (>= 1)')
            ok &= good
            if not good:
                summary.append(f'm={m}: FAIL orders {ords} tail {tail}')
            # control: other-twist points, generic non-vanishing (order 0 for all n beyond the roots)
            rc = orders(P[m], parity, K, s_ctrl)
            if rc[2] != 0:
                summary.append(f'm={m}: CONTROL tail order {rc[2]} at {labc}')
        allok &= ok
        pr(f'  {key[0]:4s} moment, source {key[1]:5s}, K={K}, components m={min(P)}..{max(P)} ({len(P)} nonzero):' +
           (f' simple zero at every {lab}, all n >= 0: {ok}' if key[0] == 'zero' else f' zero of order >= 1 at every {lab}, all n >= 0: {ok}') + ('' if ok else '  ' + '; '.join(summary)))
        if summary and ok:
            pr('     ' + '; '.join(summary))
    # closed-form identity (zeroth moment, source +2)
    K, P = r[('zero', 'plus')]
    if closed_head(name, J, 0) is not None:
        G = G_expr(parity, K)
        bad = []
        for m in range(-J, J + 1):
            cl = closed_head(name, J, m)
            val = G*P.get(m, 0)
            if cl == 0:
                if val != 0:
                    bad.append(m)
                continue
            rat = sp.simplify(sp.gammasimp(val/cl))
            if rat != 1:
                bad.append((m, rat))
        pr(f'  zeroth moment, source +2: M_m = closed-form head x ladder factor, exactly in Delta, all m: {not bad}' + (f' {bad}' if bad else ''))
        allok &= not bad
    if name == 'odd' and J == 0:      # |K|^2 = 12 |M_0|^2 / N_{O,0} against KERNELS_D3 §5.2 K_o^2 (numerically, 3 points)
        import mpmath as mp
        mp.mp.dps = 30
        G = G_expr(parity, K)
        rats = []
        for Dv in (sp.Rational(37, 7), sp.Rational(113, 9), sp.Rational(51, 4)):
            M = complex(sp.N((G*P[0]).subs(S, Dv/2), 30)); d = mp.mpf(Dv.p)/Dv.q
            N0 = mp.pi**2.5*mp.mpf(2)**(4 - 2*d)/(mp.gamma(d)*mp.gamma(d - 0.5))
            Ko2 = mp.mpf(27)/4*mp.pi**4.5*mp.gamma(d)*mp.gamma(d - 0.5)/(mp.gamma(2 + d/2)**2*mp.gamma((1 + d)/2)**2
                                                                     *mp.gamma((3 + d)/2)**2*mp.gamma((7 - d)/2)**2)
            rats.append(mp.nstr(12*abs(M)**2/N0/Ko2, 15))
        pr(f'  J = 0: 12|M_0|^2/N_0 / K_o^2(KERNELS_D3 §5.2) at Delta = 37/7, 113/9, 51/4: {rats}')
    # mirror h = -2
    Kp, Pp = r[('zero', 'plus')]; Km, Pm = r[('zero', 'minus')]
    sig = 1 if parity == 'even' else -1
    mir = Kp == Km and all(sp.expand(Pm.get(m, 0) - sig*Pp.get(-m, 0)) == 0 for m in range(-J - 1, J + 2))
    pr(f'  source -2 = mirror of source +2 (P^-_m = {sig:+d} P^+_-m): {mir}')
    pr(f'  ALL CHECKS PASSED: {allok}   [{time.time()-t0:.0f}s]')
    return out

if __name__ == '__main__':
    name = sys.argv[1]
    for J in [int(x) for x in sys.argv[2:]]:
        write(f'oz2_{name}_J{J}.txt', run(name, J))
