"""The functional as printed in the published version JHEP05(2023)122, eq. (3.3a):
     F^JHEP_g3 = int (1-p)^3[...B2^(1) + ...B3^(1)] - 5 d_{p^2}( B4^(1) + B4^(2) )|_{p=0}
and (3.3b) F^JHEP_g4 = int ... - (7 B4^(1) + d_{p^2}B4^(1) + d_{p^2}B4^(2))|_{p=0}.
Exact checks (sympy / Fractions):
 (a) low-energy side of d_{p^2} B4^(2)|_0 (EFT (2.7a) + matter (B.1)): identically zero, so the -5 g5 of the
     B4^(1) term survives in F^JHEP_g3 (and -g5 in F^JHEP_g4);
 (b) heavy action: B_k^(2) (fixed t = -p^2) only sees +- states; from (2.16a) the s-channel residue of f(s,u) at
     s=m^2 is |c+-|^2 dtilde44(1+2u/m^2), and at fixed t, u = p^2-m^2, so the per-state density of B4^(2) is
     (2m^2-p^2) P^{(0,8)}_{J-4}(-1+2p^2/m^2) / (m^12 (m^2-p^2)^3)   [= (2.34b) with (-1)^J dtilde_{4,-4}].
     Hence  m^10 F^JHEP(+-) = m^10 F^lit(+-) + 5 (-1)^J binom(J+4,8)(2JJ-85)/(9 m^8)  (checked exactly, J=4..40),
            m^10 F^JHEP(++) = m^10 F^lit(++);
 (c) counterexamples: m = M, J = 4 and J = 7 (+-) give negative actions (exact rationals).
"""
import sympy as sp
from fractions import Fraction as Fr
from math import comb, factorial
from common import poly_u

ok = True
s, p, G, g3, g4, g5, g6, g6p, ml, a0, a2 = sp.symbols('s p G g3 g4 g5 g6 g6p m_l a0 a2', positive=True)
def flow(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U
def fmat(S, U):
    T = -S - U
    return a0/(ml**2 - T) + a2*(1 - 6*S*U/ml**4)/(ml**2 - T)
# (a) B4^(2) = oint ds/(4 pi i) (2s-p^2) f(s,p^2-s)/[s(s-p^2)]^3 ; fixed t=-p^2
for name, f, poles in [('EFT', flow, [0, p**2]), ('matter', fmat, [0, p**2, ml**2, p**2 - ml**2])]:
    e = (2*s - p**2)*f(s, p**2 - s)/(s*(s - p**2))**3
    low = sp.simplify(sum(sp.residue(e, s, z) for z in poles)/2)
    print('(a) -B4^(2)|low [%s] =' % name, low)
    ok &= (sp.simplify(low) == 0)

# (b) heavy action of the forward term, per +- state, from the exact Jacobi polynomial
def jac_poly_z(n):
    """P^{(0,8)}_n(-1+2z) as exact coefficient list in z: 2F1(-n,n+9;1;1-z)."""
    # 2F1(-n,n+9;1;w) = sum_k (-n)_k (n+9)_k/(k!)^2 w^k, w = 1-z
    coef = [Fr(0)] * (n + 2)
    for k in range(n + 1):
        c = Fr(1)
        for i in range(k):
            c *= Fr((-n + i) * (n + 9 + i), (i + 1) ** 2)
        for j in range(k + 1):          # (1-z)^k
            coef[j] += c * comb(k, j) * (-1) ** j
    return coef
mismatch = 0
for J in range(4, 41):
    q = jac_poly_z(J - 4)                # D(z), z = p^2/m^2 ; set m = 1 (scaling: result ~ m^-18)
    D0, D1 = q[0], q[1]
    # d/dz [ (2-z)(1-z)^-3 D(z) ] at 0 = 5 D0 + 2 D1 ;  -5 x that, in m^10 F units at m=1
    val = -5 * (5 * D0 + 2 * D1)
    JJ = J * (J + 1)
    closed = Fr(5 * (-1) ** J * comb(J + 4, 8) * (2 * JJ - 85), 9)
    mismatch += (val != closed)
print('(b) forward-term closed form 5(-1)^J binom(J+4,8)(2JJ-85)/(9 m^8) verified for J=4..40:', mismatch == 0)
ok &= mismatch == 0

# (c) counterexamples at m = M (u = 1): m^10 F = FF(u=1)
def lit_pm(J):   # literal +- action at u=1 (poly_u uses Phi=10JJ-210; literal adds +5)
    return sum(poly_u(J, '+-')) + 5
for J in [4, 7, 9]:
    lit = lit_pm(J); extra = Fr(5 * (-1) ** J * comb(J + 4, 8) * (2 * J * (J + 1) - 85), 9)
    tot = lit + extra; tot_g5 = tot - 5
    print('(c) +- J=%d m=M: literal %s (%.4f), JHEP %s (%.2f), JHEP+5B5 %s : %s' % (
        J, lit, float(lit), tot, float(tot), tot_g5, 'NEGATIVE' if tot < 0 else 'positive'))
c4 = lit_pm(4) + Fr(5 * comb(8, 8) * (40 - 85), 9)
ok &= (c4 == Fr(-6281, 315))
# J = 4 as polynomial in u = 1/m^2: m^4 F^JHEP = 2591/210 - 41u/18 - 5u^3 - 25u^7
x = sp.Symbol('u')
P4 = sp.Rational(2591, 210) - sp.Rational(41, 18)*x - 5*x**3 - 25*x**7
assert sp.expand(P4 - (sum(c*x**i for i, c in enumerate([sp.Rational(c.numerator, c.denominator) for c in poly_u(4, '+-')])) + 5*x**3 - 25*x**7)) == 0
roots = [r for r in sp.Poly(P4, x).intervals(eps=sp.Rational(1, 10**8)) if r[0][0] >= 0 and r[0][1] <= 1]
assert len(roots) == 1 and P4.subs(x, 0) > 0 and P4.subs(x, 1) < 0
(a, b), mult = roots[0]
print('(c) J=4: m^4 F^JHEP = 2591/210 - 41u/18 - 5u^3 - 25u^7 has exactly one root u* in [%.8f, %.8f] on [0,1];' % (float(a), float(b)),
      'F^JHEP(J=4) < 0 for 1 <= m < m* = %.6f' % (1 / float(b) ** 0.5))
print('VERDICT cert_jhep: %s' % ('PASS' if ok else 'FAIL'))
