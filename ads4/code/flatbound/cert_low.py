"""Low-energy side of the F_{g3} functional (exact, sympy) and the resulting bound on |ghat_3|^2.

1. -B_k^{(1)}(p^2)|_low for k=2..5 from the tree-level EFT amplitude (2.7a) and from light matter (B.1), as
   (1/2) * sum of residues of the fixed-u integrand (2.26) at the low-energy poles (check of (2.31),(2.32),(3.8)).
2. Low-energy action of
      F_lit = int_eps^1 dp [psi2 B2 + psi3 B3] - 5 d_{p^2} B4|_0            (as printed, (3.3a))
      F_g5  = F_lit + 5 B5|_0                                                 (g5-free completion)
3. The bound |ghat3|^2 <= 16 L(eps)/k(eps) (pure gravity EFT), with exact L, k, and the explicit constants A, B
   for eps <= eps0 = 1/100; positivity of the matter polynomials on 0 < m_l < 1 (exact Sturm).
"""
import sympy as sp

s, p, G, g3, g4, g5, g6, g6p, ml, a0, a2, eps = sp.symbols('s p G g3 g4 g5 g6 g6p m_l a0 a2 eps', positive=True)

def flow(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U

def fmatter(S, U):
    T = -S - U
    return a0/(ml**2 - T) + a2*(1 - 6*S*U/ml**4)/(ml**2 - T)

def minus_B_low(k, f, poles):
    br_even = f(s, -p**2) + f(p**2 - s, -p**2)
    br_odd = -f(s, -p**2) + f(p**2 - s, -p**2)
    if k % 2 == 0:
        e = (2*s - p**2)/(s*(s - p**2))**sp.Rational(k - 2, 2) * br_even
    else:
        e = 1/(s*(s - p**2))**sp.Rational(k - 3, 2) * br_odd
    return sp.simplify(sum(sp.residue(e, s, z) for z in poles)/2)

ok = True
print('--- EFT (2.7a), |ghat3|^2 written g3 ---')
low = {}
for k in [2, 3, 4, 5]:
    low[k] = sp.expand(minus_B_low(k, flow, [0, p**2]))
    print('-B_%d|low =' % k, low[k])
ok &= sp.simplify(low[2] - (16*sp.pi*G/p**2 + 2*sp.pi*G*g3*p**6)) == 0
ok &= sp.simplify(low[3] - (-2*sp.pi*G*g3*p**4)) == 0
ok &= sp.simplify(low[4] - (2*g4 + (4*sp.pi*G*g3 + g5)*p**2 + (g6 + g6p)*p**4)) == 0
ok &= sp.simplify(low[5] - (g5 + (g6 - g6p)*p**2)) == 0
print('agrees with (2.31a,b),(2.32a,b):', ok)

print('--- matter (B.1): a0=|g0|^2, a2=|g2|^2 ---')
mat = {}
for k in [2, 3, 4, 5]:
    mat[k] = sp.factor(sp.expand(minus_B_low(k, fmatter, [0, p**2, ml**2, p**2 - ml**2])))
    print('-B_%d|matter =' % k, mat[k])

psi2 = (1 - p)**3*p*(65 + 155*p + 47*p**2)
psi3 = 2*(1 - p)**3*(81 - 1927*p + 4601*p**2)

def lint(expr, a):
    """int_a^1 of a Laurent polynomial in p (termwise, exact)."""
    e = sp.expand(expr); tot = 0
    for term in sp.Add.make_args(e):
        c, k = term.as_coeff_exponent(p)
        tot += c*(-sp.log(a) if k == -1 else (1 - a**(k + 1))/(k + 1))
    return sp.expand(tot)

def functional_low(expr, eps_, g5free):
    tot = lint(psi2*expr[2], eps_) + lint(psi3*expr[3], eps_)
    xx = sp.Symbol('xx', positive=True)
    d4 = sp.diff(sp.expand(expr[4]).subs(p, sp.sqrt(xx)), xx).subs(xx, 0)
    tot += -5*d4
    if g5free:
        tot += 5*sp.expand(expr[5]).subs(p, 0)
    return sp.expand(tot)

for g5free in [False, True]:
    name = 'F_g5free' if g5free else 'F_literal'
    Fl = sp.expand(functional_low(low, eps, g5free))
    print('--- %s: low-energy action (EFT part) ---' % name)
    print(sp.collect(Fl, [g3, g5]))
    L = sp.integrate(sp.expand(psi2/p**2), (p, eps, 1))
    kfun = sp.simplify(-sp.expand(Fl).coeff(g3)/(sp.pi*G))
    print('  coefficient of G pi: 16 L(eps), L =', sp.simplify(L))
    print('  k(eps) =', sp.expand(kfun), ' k(0) =', kfun.subs(eps, 0), '=', sp.N(kfun.subs(eps, 0), 12))
    print('  g5 coefficient:', sp.expand(Fl).coeff(g5))

# matter part of the action (g5-free), coefficient of a0 and a2
Fm = sp.expand(functional_low({k: sp.expand(v) for k, v in mat.items()}, 0, True))
Fm_lit = sp.expand(functional_low({k: sp.expand(v) for k, v in mat.items()}, 0, False))
c0 = sp.factor(Fm.coeff(a0)); c2 = sp.factor(Fm.coeff(a2))
print('--- matter action (g5-free) = a0*C0 + a2*C2 ;  literal: a0*%s + a2*%s' % (sp.factor(Fm_lit.coeff(a0)), sp.factor(Fm_lit.coeff(a2))))
print('C0 =', c0); print('C2 =', c2)
# The bound: action = 16 pi G L - pi G k |g3|^2 + C0 a0 + C2 a2 >= 0 ; dropping matter needs C0, C2 <= 0 on 0<ml<1
x = sp.Symbol('x')
okm = True
for nm, cc in [('C0', c0), ('C2', c2), ('C0lit', sp.factor(Fm_lit.coeff(a0))), ('C2lit', sp.factor(Fm_lit.coeff(a2)))]:
    num = sp.numer(sp.together(cc)); den = sp.denom(sp.together(cc))
    P = sp.Poly(sp.expand(num.subs(ml, x)), x); D = sp.Poly(sp.expand(den.subs(ml, x)), x)
    nroots_num = sp.count_roots(P, 0, 1) if P.degree() > 0 else 0
    # denominator must be a positive constant times a power of m_l (positive on (0,1])
    den_ok = len(D.terms()) == 1 and D.LC() > 0
    sign_mid = sp.sign(cc.subs(ml, sp.Rational(1, 2)))
    neg = (sign_mid < 0) and nroots_num == 0 and den_ok
    print('  %s < 0 on (0,1]: roots(num) in [0,1]=%d, den=c*m^j with c>0: %s, sign(1/2)=%s -> %s' % (nm, nroots_num, den_ok, sign_mid, neg))
    okm &= bool(neg)

# ---- F_g4 (3.3b): int (1-p)^3[p(99+282p+49p^2)B2 + (21-474p+1216p^2)B3] - (7 + d_{p^2}) B4|_0
psi2g4 = (1 - p)**3*p*(99 + 282*p + 49*p**2)
psi3g4 = (1 - p)**3*(21 - 474*p + 1216*p**2)
xx4 = sp.Symbol('xx4', positive=True)
Fg4 = sp.expand(sp.integrate(sp.expand(psi2g4*(low[2] - 16*sp.pi*G/p**2)), (p, 0, 1)) + sp.integrate(sp.expand(psi3g4*low[3]), (p, 0, 1))
                - 7*sp.expand(low[4]).subs(p, 0) - sp.diff(sp.expand(low[4]).subs(p, sp.sqrt(xx4)), xx4).subs(xx4, 0))
print('--- F_g4 (v1 (3.3b)) low side, non-pole part:', sp.collect(Fg4, [g3, g4, g5]))
c_g5_g4 = Fg4.coeff(g5); c_g4 = Fg4.coeff(g4)
print('  coefficient of g5 = %s (so (3.4b) also needs +B5|_0 or g5>=0); g4 coefficient = %s ; bound g4 <= 8piG(%s log(1/eps) + ...), %s = 14.14 ; |g3|^2 coefficient/(8 pi G*14) = %s' % (
    c_g5_g4, c_g4, sp.Rational(99, 7), sp.Rational(99, 7), sp.nsimplify(-Fg4.coeff(g3)/(8*sp.pi*G*14))))
ok &= (c_g5_g4 == -1)
# constants A, B for eps <= eps0
eps0 = sp.Rational(1, 100)
k0 = sp.Rational(330151, 12012)
Lfull = sp.integrate(sp.expand(psi2/p**2), (p, eps, 1))
ell = sp.simplify(Lfull - 65*sp.log(1/eps) + sp.Rational(1171, 15))
print('ell(eps) = L - 65 log(1/eps) + 1171/15 =', sp.expand(ell))
# k(eps) in [k0 - 2 int_0^eps |psi3| p^4, k0 + 2 int_0^eps psi2 p^6]; psi3>0 on [0,0.04]
klo = k0 - 2*sp.integrate(sp.expand(psi3*p**4), (p, 0, eps0))
khi = k0 + 2*sp.integrate(sp.expand(psi2*p**6), (p, 0, eps0))
assert sp.count_roots(sp.Poly(sp.expand(81 - 1927*x + 4601*x**2), x), 0, sp.Rational(1, 25)) == 0
Aconst = 16*65/klo
Bconst = 16*(-sp.Rational(1171, 15) + ell.subs(eps, eps0))/khi
print('for 0 < eps <= 1/100:  |ghat3|^2 M^8 <= A log(M/m_IR) + B,  A = %s <= %.6f,  B = %s <= %.6f' % (
    sp.nsimplify(Aconst), float(sp.N(Aconst, 20)) + 1e-6, '', float(sp.N(Bconst, 20)) + 1e-6))
print('limit eps->0 constants: A0 = 1040/k0 = %.6f, B0 = 16(-1171/15)/k0 = %.6f' % (float(1040/k0), float(16*sp.Rational(-1171, 15)/k0)))
print('VERDICT cert_low: %s' % ('PASS' if ok and okm else 'FAIL'))
