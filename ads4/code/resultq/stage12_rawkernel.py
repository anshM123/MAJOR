"""Result Q, Stage 1 (RESULT_Q_STAGE12.md §1.3): the RAW (non-de-shifted) leading-Regge heavy image of the subtracted k = 4 functional,
SCALAR PROXY.  stage12_kernel.py shows f+_{4;2} = f+_{2;2} x 16 eta^2 r^2 (1 + O(r)): the k = 4 leading kernel is r^2 eta^2 x (k = 2 kernel),
i.e. eta^3 x (f-stripped kernel).  For the scalar functional the f-stripped leading kernel is a single H_2 harmonic P_l(eta), l = -1/2 + i nu,
whose plane-wave (heavy) image is E_b(nu) P_l(cosh beta) with (v3_regge_kernel's K(c, lam) for both sides, radial weight b)
    E_b(nu) = 2^{2b-2} pi Gamma((b - 1/2 - i nu)/2)^2 Gamma((b - 1/2 + i nu)/2)^2       (b = 7 for k = 2, b = 9 for k = 4).
Multiplication by eta acts on harmonics by the Legendre recurrence eta P_mu = [(mu+1) P_{mu+1} + mu P_{mu-1}]/(2mu+1), so
    heavy image of eta^n P_l = sum_j a^{(n)}_j(l) E_b(nu - i j) P_{l+j}(cosh beta).
Check: n = 1, b = 7 reproduces CHMRS's magic-block ratio g/g~ (H7 §18-§23: h_-/h_+ = g/g~) exactly.
Then n = 3, b = 9 (k = 4): shifts j = +-1, +-3, and the forward (nu -> 0) image contains P_{5/2}(cosh beta) ~ (2J/m)^5.
EE spin factor: for k = 2 the EE image is the scalar one times c(nu') on each harmonic nu' (EE_FUNCTIONAL §2); whether the same holds at
b = 9, and with which factor, is NOT computed here [I]; c(nu') is negative at nu' = +-i, +-3i (ratio 206).
Output: cache/stage12_rawkernel.txt.  Seconds.
"""
import os, io
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()
mp.mp.dps = 40


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def recurrence_coeffs(n):
    mu = sp.Symbol('mu')
    coeffs = {0: sp.Integer(1)}
    for _ in range(n):
        new = {}
        for j, c in coeffs.items():
            m = mu + j
            new[j + 1] = sp.simplify(new.get(j + 1, 0) + c * (m + 1) / (2 * m + 1))
            new[j - 1] = sp.simplify(new.get(j - 1, 0) + c * m / (2 * m + 1))
        coeffs = {j: sp.factor(c) for j, c in new.items() if sp.simplify(c) != 0}
    return mu, coeffs


def E(b, nu):
    return 2 ** (2 * b - 2) * mp.pi * mp.gamma((b - mp.mpf(1) / 2 - 1j * nu) / 2) ** 2 * mp.gamma((b - mp.mpf(1) / 2 + 1j * nu) / 2) ** 2


def P(mu, beta):
    return mp.legenp(mu, 0, mp.cosh(beta), type=3)


def main():
    pr("# stage12_rawkernel: scalar proxy for the raw leading-Regge image of the k = 4 subtracted functional")
    mu, c1 = recurrence_coeffs(1)
    _, c3 = recurrence_coeffs(3)
    pr("eta P_mu   =", {j: str(v) for j, v in sorted(c1.items())})
    pr("eta^3 P_mu =", {j: str(v) for j, v in sorted(c3.items())})
    # check n = 1, b = 7 against CHMRS g/g~
    for nu in (mp.mpf('0.3'), mp.mpf('2.7'), mp.mpf('9.1')):
        l = -mp.mpf(1) / 2 + 1j * nu
        hm = complex(c1[-1].subs(mu, sp.Float(0)).evalf()) if False else None
        am = (l) / (2 * l + 1); ap = (l + 1) / (2 * l + 1)
        hm = am * E(7, nu + 1j); hp = ap * E(7, nu - 1j)
        lg = lambda w: mp.gamma((mp.mpf(13) / 2 - 1j * w) / 2) * mp.gamma((mp.mpf(13) / 2 + 1j * w) / 2)
        g = (-mp.mpf(1) / 2 + 1j * nu) / (1j * nu) * (lg(nu + 1j) / lg(nu)) ** 2
        gt = (-mp.mpf(1) / 2 - 1j * nu) / (-1j * nu) * (lg(-nu + 1j) / lg(-nu)) ** 2
        pr(f"   check k=2 (b=7): nu = {float(nu)}: (h_-/h_+)/(g/g~) - 1 = {complex((hm / hp) / (g / gt) - 1):.2e}")
    # k = 4 proxy: forward image near nu = 0, large beta
    pr("k = 4 scalar proxy (eta^3, b = 9): normalised raw image H(beta) = sum_j a_j E_9(nu - i j) P_{l+j}(cosh beta) / sum_j a_j E_9(nu - i j)")
    nu = mp.mpf('1e-12')
    l = -mp.mpf(1) / 2 + 1j * nu
    a = {j: complex(sp.N(c3[j].subs(mu, sp.Float(str(-0.5)) + sp.I * sp.Float('1e-12')), 40)) for j in c3}
    a = {j: mp.mpc(sp.N(sp.re(c3[j].subs(mu, sp.Rational(-1, 2) + sp.I * sp.Rational(1, 10 ** 12))), 40),
                   sp.N(sp.im(c3[j].subs(mu, sp.Rational(-1, 2) + sp.I * sp.Rational(1, 10 ** 12))), 40)) for j in c3}
    w = {j: a[j] * E(9, nu - 1j * j) for j in a}
    Nf = sum(w.values())
    pr(f"   forward normalisation N(0) = {mp.nstr(mp.re(Nf), 8)} (imag {mp.nstr(mp.im(Nf), 3)});  weights/N: " +
       ", ".join(f"j={j}: {mp.nstr(mp.re(w[j] / Nf), 6)}{'+' if mp.im(w[j] / Nf) >= 0 else ''}{mp.nstr(mp.im(w[j] / Nf), 3)}i" for j in sorted(w)))
    pr("   beta    H(beta) (raw, nu=0)     P_{-1/2}(cosh b)   P_{5/2}(cosh b)    H/P_{5/2}")
    for beta in (0.01, 0.1, 0.5, 1, 2, 4, 6, 8):
        H = sum(w[j] * P(l + j, beta) for j in w) / Nf
        pr(f"   {beta:5.2f}   {mp.nstr(mp.re(H), 10):>18s}   {mp.nstr(mp.re(P(-0.5, beta)), 8):>14s}   {mp.nstr(mp.re(P(2.5, beta)), 8):>14s}   {mp.nstr(mp.re(H) / mp.re(P(2.5, beta)), 6)}")
    cfun = lambda v: 32 * (4 * v ** 2 + 1) * (4 * v ** 2 + 9) / ((4 * v ** 2 + 25) * (4 * v ** 2 + 81))
    pr(f"   EE factor of the k = 2 functional at the shifted harmonics: c(+-i) = {complex(cfun(1j)).real:.4f}, c(+-3i) = {complex(cfun(3j)).real:.4f}"
       " (both negative, ratio 206: IF the k = 2 rule 'EE = scalar x c(nu') per harmonic' held at b = 9 with this c, the j = +-3 weight would'"
       " dominate even more and H > 0; whether it holds at b = 9 is not computed) [I]")
    open(os.path.join(CACHE, 'stage12_rawkernel.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
