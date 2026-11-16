"""Exact low-energy side of a robust functional, recomputed from the EFT amplitude (2.7a) and light matter (B.1) by residues
(same method as flatbound/cert_low.py, whose -B_k|_low expressions are re-derived here, not imported).
Checks for the functional F = int psi2 B2 + psi3 B3 - a d_{p^2}B4|_0 + a B5|_0 (eps -> 0, pole term separated):
  * g4, g5 (and g6, g6') drop out;  * coefficient of |ghat3|^2 is -pi G k with k = -2 mu2_6 + 2 mu3_4 + 4a;
  * graviton pole gives 16 pi G int psi2/p^2;  * matter coefficients equal the C0, C2 used in rcert.py (negative on (0,1]).
Usage: python3 lowcheck.py NAME"""
import sys, os, json
import sympy as sp

s, p, G, g3, g4, g5, g6, g6p, ml, a0, a2 = sp.symbols('s p G g3 g4 g5 g6 g6p m_l a0 a2', positive=True)

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
    return sp.expand(sp.simplify(sum(sp.residue(e, s, z) for z in poles)/2))

def main(name):
    d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'functionals.json')))[name]
    psi2 = sum(sp.Rational(c)*p**i for i, c in enumerate(d['psi2_exact']))
    psi3 = sum(sp.Rational(c)*p**i for i, c in enumerate(d['psi3_exact']))
    a = sp.Rational(d['a_exact'])
    low = {k: minus_B_low(k, flow, [0, p**2]) for k in [2, 3, 4, 5]}
    mat = {k: minus_B_low(k, fmatter, [0, p**2, ml**2, p**2 - ml**2]) for k in [2, 3, 4, 5]}
    xx = sp.Symbol('xx', positive=True)
    def act(ex, pole_sub=True):
        e2 = sp.expand(ex[2])
        pole = e2.coeff(p, -2)
        reg = sp.expand(e2 - pole/p**2)
        tot = sp.integrate(sp.expand(psi2*reg), (p, 0, 1)) + sp.integrate(sp.expand(psi3*ex[3]), (p, 0, 1))
        tot += -a*sp.diff(sp.expand(ex[4]).subs(p, sp.sqrt(xx)), xx).subs(xx, 0) + a*sp.expand(ex[5]).subs(p, 0)
        return sp.expand(tot), pole
    tot, pole = act(low)
    k = -2*sp.integrate(psi2*p**6, (p, 0, 1)) + 2*sp.integrate(psi3*p**4, (p, 0, 1)) + 4*a
    ok = True
    for cpl in [g4, g5, g6, g6p]:
        c = tot.coeff(cpl); print('coefficient of %s: %s' % (cpl, c)); ok &= (c == 0)
    cg3 = sp.simplify(tot.coeff(g3) + sp.pi*G*k); print('coefficient of |g3|^2 + pi G k =', cg3); ok &= (cg3 == 0)
    print('pole of -B2|low: %s / p^2  (expect 16 pi G)' % pole); ok &= sp.simplify(pole - 16*sp.pi*G) == 0
    print('k =', k, '=', sp.N(k, 12))
    tm, _ = act(mat)
    C0 = sp.factor(tm.coeff(a0)); C2 = sp.factor(tm.coeff(a2))
    print('matter: C0 =', C0, ';  C2 =', C2)
    x = sp.Symbol('x')
    for nm, cc in [('C0', C0), ('C2', C2)]:
        num = sp.numer(sp.together(cc)); den = sp.denom(sp.together(cc))
        P = sp.Poly(sp.expand(num.subs(ml, x)), x); D = sp.Poly(sp.expand(den.subs(ml, x)), x)
        nr = sp.count_roots(P, 0, 1) if P.degree() > 0 else 0
        neg = nr == 0 and len(D.terms()) == 1 and D.LC() > 0 and sp.sign(cc.subs(ml, sp.Rational(1, 2))) < 0 and P.eval(0) < 0
        print('  %s < 0 on (0,1]: %s' % (nm, neg)); ok &= bool(neg)
    print('VERDICT lowcheck %s: %s' % (name, 'PASS' if ok else 'FAIL'))

if __name__ == '__main__':
    main(sys.argv[1])
