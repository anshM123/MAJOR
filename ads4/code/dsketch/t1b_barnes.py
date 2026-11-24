"""
t1b_barnes (DSKETCH_UP.md §2.3): the per-nu image of a top-Regge coefficient that is POLYNOMIAL in the Mellin variable T is a
polynomial in nu^2 of the same degree, with no non-polynomial remainder (this replaces input (E) for the O(t4) term).

CHMRS (2106.10274) eq. (CfromBMellin):  C_{k,nu}[Bhat] = N_k/gamma_a(nu)^2 * int dT/(4 pi i) Gamma(T/2)^2 gamma_{a-T}(nu) Bhat(T),
gamma_x(nu) = Gamma((x - 1/2 + i nu)/2) Gamma((x - 1/2 - i nu)/2),  A := a - 1/2.  With z = T/2, c, d = (A +- i nu)/2, and
Bhat(T) = (T/2)_j (rising factorial), Barnes' first lemma (Gamma(z)(z)_j = Gamma(z+j)) gives exactly
    C[(T/2)_j] / C[1] = (c)_j (d)_j / (A)_j = prod_{i<j} ((A/2 + i)^2 + nu^2/4) / (A)_j ,
a polynomial in nu^2 of degree j; j = 0 gives a nu-independent constant.  Since {(T/2)_j} is a basis of polynomials, every polynomial
Bhat of degree n maps to a polynomial in nu^2 of degree n.  In particular T -> (A^2 + nu^2)/(2A): a Mellin constant stays a constant,
and a Mellin u-linear term becomes nu^2 + A^2 (the AdS shift used in DSKETCH_UP §3.4).
Checks: (1) the Barnes integral by direct numerical contour integration (mpmath, 30 digits) for several (A, nu, j);
        (2) the closed form against sympy's product, symbolic in nu.
Writes cache/t1b_barnes.txt
"""
import os, sys
sys.dont_write_bytecode = True
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
lines = []
def log(x=''):
    print(x, flush=True); lines.append(str(x))

ok = True
log('(1) numerical Barnes integrals  I_j = int_{c0 - i oo}^{c0 + i oo} dz/(2 pi i) Gamma(z)^2 Gamma(c-z) Gamma(d-z) (z)_j')
for A in (mp.mpf('6.5'), mp.mpf('7.5'), mp.mpf('4.5')):
    for nu in (mp.mpf('0.7'), mp.mpf('3'), mp.mpf('11')):
        c = (A + 1j*nu)/2; d = (A - 1j*nu)/2
        c0 = A/4
        I0 = None
        for j in range(4):
            f = lambda t: mp.gamma(c0 + 1j*t)**2*mp.gamma(c - c0 - 1j*t)*mp.gamma(d - c0 - 1j*t)*mp.rf(c0 + 1j*t, j)/(2*mp.pi)
            num = mp.quad(f, [-mp.inf, -40, -10, 0, 10, 40, mp.inf])
            exact = mp.gamma(c + j)*mp.gamma(d + j)*mp.gamma(c)*mp.gamma(d)/mp.gamma(c + d + j)
            rel = abs(num/exact - 1)
            if j == 0:
                I0 = exact
            ratio = (exact/I0).real
            poly = mp.mpf(1)
            for i in range(j):
                poly *= (A/2 + i)**2 + nu**2/4
            poly /= mp.rf(A, j)
            log(f'  A = {mp.nstr(A, 3)}, nu = {mp.nstr(nu, 3)}, j = {j}: |num/exact - 1| = {mp.nstr(rel, 3)};  '
                f'I_j/I_0 = {mp.nstr(ratio, 12)}  closed form {mp.nstr(poly, 12)}')
            ok &= rel < 1e-20 and abs(ratio - poly) < 1e-20*abs(poly)
nu_s, A_s = sp.symbols('nu A', positive=True)
log('(2) symbolic images of the monomials T^n (T = 2z), A symbolic:')
for n in range(4):
    # T^n = 2^n z^n = 2^n sum_j S(n, j) (-1)^{n-j} (z)_j ... use falling->rising via sympy: expand z^n in rising factorials
    zz = sp.symbols('zz')
    # solve z^n = sum_j a_j (z)_j
    a = sp.symbols(f'a0:{n + 1}')
    expr = sum(a[j]*sp.rf(zz, j) for j in range(n + 1)) - zz**n
    sol = sp.solve(sp.Poly(sp.expand(expr), zz).coeffs(), a, dict=True)[0]
    img = 0
    for j in range(n + 1):
        cj = (A_s + sp.I*nu_s)/2; dj = (A_s - sp.I*nu_s)/2
        img += sol[a[j]]*sp.expand(sp.rf(cj, j)*sp.rf(dj, j))/sp.rf(A_s, j)
    img = sp.factor(sp.simplify(2**n*img))
    isp = sp.Poly(sp.expand(sp.numer(sp.together(img))), nu_s)
    even = all(m[0] % 2 == 0 for m in isp.monoms())
    deg = sp.degree(sp.numer(sp.together(img)), nu_s)
    log(f'  T^{n}  ->  {img}    (polynomial in nu^2: {even}, degree in nu: {deg})')
    ok &= even and deg == 2*n
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', 't1b_barnes.txt'), 'w').write('\n'.join(lines) + '\n')
