"""Certificate EPS: localisation of the negative set of the IR-truncated functional, and Delta_IR under E(lambda,eta).

F^eps = F - R_eps,  m^4 R_eps = int_0^eps Psi_h D_h dp,  |m^4 R_eps| <= q(65 eps^2 + 162 u eps + 17 u eps^4) <= r0 + r1 u
(q = (1-eps^2)^-4; psi2 <= 65p, 0 < psi3 <= 162 on [0,0.04], |P_J| <= 1, |dtilde44| <= (1-y)^-4), for all eps <= EPS.
Part A (certified): for every eps <= EPS, every physical heavy state OUTSIDE
     Omega_eps := { bhat >= 120, u <= min(1/100, uN(bhat)) },
     uN(b) = [ (162/b_- + r1) + sqrt((162/b_- + r1)^2 + 10(1-eta)b^2 r0) ] / (5(1-eta) b^2),
has m^4 F > r0 + r1 u >= m^4|R_eps|, hence F^eps > 0.  (P) exact Descartes, (C) crude, (S) Cauchy cells, (M) analytic.
Part B (proved + evaluated): per state in Omega_eps,  m^4[F^eps]_- <= [m^4 R_eps]_+ <= 2(65 eps [J1(eps b)]_+/b + 14eps^3)
     + r1 u + 130(2.3 u eps^4 + 0.024 b u eps^5 + 1.6 sqrt(u) eps^3),  with u <= uN(b) <= 65.5/b^3 + 65.5eps/b^2 + 5.11eps/b;
under E(lambda,eta):  Delta_IR <= (1+eta)(16/k_lo)[65 I(lambda) + 3.5 eps/lambda^2 + X(eps,lambda)].
Usage: python3 cert_eps.py  (EPS = 1/20000)
"""
import sys, time
import flint
from fractions import Fraction as Fr
from common import poly_u, to_fmpz_poly, set_prec
import cert_P
import cert_S
import cert_C
from cert_S import fq, A
set_prec(cert_S.PREC)

EPS = Fr(1, 20000)
JP = 400; BA = Fr(120); U1 = Fr(1, 100)
q = 1 / (1 - EPS ** 2) ** 4
r0 = 65 * EPS ** 2 * q
r1 = (162 * EPS + 17 * EPS ** 4) * q
ok_all = True
t0 = time.time()
print('EPS = %s ; r0 = %.6e ; r1 = %.6e' % (EPS, float(r0), float(r1)))

# ---- (P): J < JP, all u in [0,1]: FF_J(u) - r0 - r1 u > 0
cnt = 0
for h in ['++', '+-']:
    Js = range(0, JP, 2) if h == '++' else range(4, JP)
    for J in Js:
        c = poly_u(J, h)
        c[0] -= r0; c[1] -= r1
        P, L = to_fmpz_poly(c)
        ok, n = cert_P.positive_on_01(P)
        cnt += 1
        if not ok:
            print('EPS-P %s J=%d FAIL' % (h, J))
        ok_all &= ok
print('EPS-P: %d polynomials FF_J(u) - r0 - r1 u > 0 on [0,1]: %s (%.0fs)' % (cnt, 'PASS' if ok_all else 'FAIL', time.time() - t0), flush=True)

# ---- (C): J >= JP, 1<=m<=2: m^10 F >= 10JJ-210 - 64 U(J) (cert_C), so FF = m^-6 m^10 F >= (10JJ-210-64U)/64
from math import comb
set_prec(200)
c_a = cert_C.integral_upper(lambda x: cert_C.g_hat(x) / ((1 - x) * (1 + x) ** 4), 0, Fr(1, 2), 2000)
gamma = cert_C.sup_on(lambda x: cert_C.g_hat(x) / (1 + x) ** 4, Fr(1, 2), 1, 4000)
g1 = cert_C.sup_on(cert_C.g_hat, Fr(1, 2), 1, 4000)
U = A(c_a) + A(gamma) * A(comb(JP + 4, 8)).log() / 4 + A(g1) / (4 * (A(3) / 2) ** 4)
margin = (A(10 * JP * (JP + 1) - 210) - 64 * U) / 64
okC = bool(margin > fq(r0) + fq(r1))      # (u <= 1); the margin increases with J (cert_C monotonicity)
print('EPS-C: FF >= %s > r0 + r1 : %s' % (margin.str(5), 'PASS' if okC else 'FAIL'), flush=True)
ok_all &= okC
set_prec(cert_S.PREC)

# ---- (S): J >= JP, m >= 2, bhat in (0,120]
for h in ['++', '+-']:
    okb, low = cert_S.small_b(h, Fr(1), JP, r0, r1)
    ak = cert_S.build_ak(h)
    b = Fr(1); ncell = 0; okS = okb
    while b < BA:
        r = min(Fr(2), (b + Fr(1, 10)) / 2); hh = min(r / 10, Fr(1, 10)); bc = b + hh
        while hh > min(Fr(2), bc / 2) / 10:
            hh = hh / 2; bc = b + hh
        if bc + hh > BA:
            hh = (BA - b) / 2; bc = b + hh
        okc, nbox, M, u0 = cert_S.cell(h, ak, bc, hh, JP, 3, r0, r1)
        if not okc: print('EPS-S %s cell [%s,%s] FAIL' % (h, float(b), float(b + 2 * hh)))
        okS &= okc; ncell += 1; b = b + 2 * hh
    print('EPS-S %s: small-b and %d cells, FF - r0 - r1 u > 0: %s (%.0fs)' % (h, ncell, 'PASS' if okS else 'FAIL', time.time() - t0), flush=True)
    ok_all &= okS

# ---- (M): bhat >= 120, m >= 2.  Constants from cert_M (recomputed): C_R(u<=u1) <= CR1, C_R(u1<=u<=1/4) <= CRm
import cert_M
set_prec(128)
CR1 = cert_M.norms('++', 0, U1, 4000)[0].max(cert_M.norms('+-', 0, U1, 4000)[0])
CRm = cert_M.norms('++', U1, Fr(1, 4), 4000)[0].max(cert_M.norms('+-', U1, Fr(1, 4), 4000)[0])
b = fq(BA); delta = 8 * fq(U1).sqrt() / b; eta = fq(85) / 4 / b ** 2
bp = b * (1 + delta); bm = b * (1 - delta)
# common lower bound (both channels), u <= u1:  FF >= 80/b+^3 - C/b-^4 - 162u/b- + 2.5(1-eta) b^2 u^2 ;  bracket >= 0:
br = bool(80 * b / (1 + delta) ** 3 - CR1 / (1 - delta) ** 4 > 0)      # (b^4 * bracket at b = 120, increasing in b)
# u in [u1, 1/4]: 2.5(1-eta) b^2 u^2 - (162/(b-4) + r1) u - r0 - CRm/(b-4)^4 > 0, increasing in u for u >= u1
a = fq(Fr(5, 2)) * (1 - eta) * b ** 2; cc = 162 / (b - 4) + fq(r1); uu = fq(U1)
mu1 = bool(a * uu > cc) and bool(a * uu ** 2 - cc * uu - fq(r0) - CRm / (b - 4) ** 4 > 0)
okM = br and mu1
print('EPS-M: bracket 80/b+^3 - C_R/b-^4 >= 0 for b>=120: %s ; u in [1/100,1/4] excluded: %s ; '
      'for u <= 1/100 the negative set lies in u < uN(b) (quadratic in u): PASS by construction' % (br, mu1))
ok_all &= okM

# ---- Part B: Delta_IR bound
set_prec(200)
import cert_ir
a1 = 162 / (fq(Fr(5, 2)) * (1 - eta) * (1 - delta)); a2 = fq(r1) / (fq(Fr(5, 2)) * (1 - eta)) / fq(EPS)
a3 = (fq(r0) / (fq(Fr(5, 2)) * (1 - eta))).sqrt() / fq(EPS)
print('uN(b) <= %s/b^3 + %s eps/b^2 + %s eps/b' % (a1.str(5), a2.str(5), a3.str(5)))
assert bool(a1 < fq(Fr(655, 10))) and bool(a2 < fq(Fr(655, 10))) and bool(a3 < fq(Fr(511, 100)))
def X_bound(eps, lam):
    """majorant of int_{120}^{L} e(b) b log(L/b) db, L = 1/(lam eps); e = per-state extras / 2."""
    e = fq(eps); L = 1 / (fq(lam) * e); b0 = fq(120)
    Lam = (L / b0).log()                       # log(L/b) <= log(L/120)
    r1h = fq(r1) / 2
    # int b*uN <= 65.5/120 + 65.5 eps log(L/120) + 5.11 eps L
    Ib_u = fq(Fr(655, 10)) / b0 + fq(Fr(655, 10)) * e * Lam + fq(Fr(511, 100)) * e * L
    # int b^2 uN <= 65.5 log(L/120) + 65.5 eps L + 5.11 eps L^2/2
    Ib2_u = fq(Fr(655, 10)) * Lam + fq(Fr(655, 10)) * e * L + fq(Fr(511, 100)) * e * L * L / 2
    # int b sqrt(uN) <= sqrt(65.5)*2 sqrt(L) + sqrt(65.5 eps) L + sqrt(5.11 eps)(2/3) L^(3/2)
    s655 = fq(Fr(655, 10)).sqrt()
    Ib_su = s655 * 2 * L.sqrt() + (fq(Fr(655, 10)) * e).sqrt() * L + (fq(Fr(511, 100)) * e).sqrt() * 2 * L ** fq(Fr(3, 2)) / 3
    tot = r1h * Ib_u + 65 * (fq(Fr(23, 10)) * e ** 4 * Ib_u + fq(Fr(24, 1000)) * e ** 5 * Ib2_u + fq(Fr(16, 10)) * e ** 3 * Ib_su)
    return tot * Lam, Lam
klo = cert_ir.klo
# sharper constants for eps <= EPS (same argument as cert_low, eps0 = EPS)
import sympy as sp
pp = sp.Symbol('p')
psi2s = sp.expand((1 - pp)**3*pp*(65 + 155*pp + 47*pp**2)); psi3s = sp.expand(2*(1 - pp)**3*(81 - 1927*pp + 4601*pp**2))
k0 = Fr(330151, 12012)
klo_s = k0 - 2 * Fr(str(sp.integrate(psi3s*pp**4, (pp, 0, sp.Rational(EPS.numerator, EPS.denominator)))))
khi_s = k0 + 2 * Fr(str(sp.integrate(psi2s*pp**6, (pp, 0, sp.Rational(EPS.numerator, EPS.denominator)))))
e = EPS
ell_s = 40*e + Fr(223, 2)*e**2 - Fr(259, 3)*e**3 + Fr(7, 2)*e**4 + Fr(47, 5)*e**5
Bc = fq(16 * (ell_s - Fr(1171, 15)) / khi_s); Ac = fq(Fr(1040) / klo_s)
print('for m_IR <= EPS*M:  A = 1040/k_lo = %s ,  B = %s' % (Ac.str(10), Bc.str(10)))
print('lambda   Lambda=log(1/(120 lambda eps))   X(eps*,lambda)   Delta_IR bound (eta=0)   B + Delta (valid for m_IR <= %s M)' % EPS)
for lam in [Fr(2, 7), Fr(1, 2), Fr(1), Fr(2), Fr(4)]:
    X, Lam = X_bound(EPS, lam)
    assert bool(Lam > 1)       # => every term eps^a*Lambda (a>=1) is increasing in eps on (0, EPS]
    D = 16 / fq(klo) * (65 * cert_ir.I_lambda(lam) + fq(Fr(7, 2)) * fq(EPS) / fq(lam) ** 2 + X)
    print('%-7s  %-30s  %-16s %-24s %s' % (lam, Lam.str(5), X.str(5), D.str(6), (Bc + D).str(7)))
print('VERDICT cert_eps (EPS=%s, %.0fs): %s' % (EPS, time.time() - t0, 'PASS' if ok_all else 'FAIL'))
