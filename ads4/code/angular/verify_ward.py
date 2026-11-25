"""ANGULAR_ROUTE §7 check: rest-frame Ward operator  T + B1^2 + B2^2 + R^2 + (4-2Delta) B1 - 2 R B2  on spin-J polynomials
e_m = a^{J+m} b^{J-m}: diagonal T + C - 2m^2, neighbour coefficients (unnormalised) and normalised -> -b_m."""
import sympy as sp
a, b, D = sp.symbols('a b Delta')
B1 = lambda f: sp.Rational(1, 2)*(b*sp.diff(f, a) + a*sp.diff(f, b))
B2 = lambda f: sp.I/2*(-b*sp.diff(f, a) + a*sp.diff(f, b))
R = lambda f: -sp.I/2*(a*sp.diff(f, a) - b*sp.diff(f, b))
T = (D-1)*(D-2)
def W(f): return sp.expand(T*f + B1(B1(f)) + B2(B2(f)) + R(R(f)) + (4-2*D)*B1(f) - 2*R(B2(f)))
for J in range(0, 7):
    C = J*(J+1)
    for m in range(-J, J+1):
        img = sp.Poly(W(a**(J+m)*b**(J-m)), a, b)
        coef = {mm: img.coeff_monomial(a**(J+mm)*b**(J-mm)) for mm in range(-J, J+1)}
        assert sp.simplify(coef[m] - (T + C - 2*m*m)) == 0, (J, m, coef[m])
        assert all(coef[mm] == 0 for mm in coef if abs(mm - m) > 1)
        # unnormalised: W e_m = sum_mm coef[mm] e_mm.  Matrix element (mm <- m).
        if m+1 <= J: up = sp.factor(coef[m+1])
        if m-1 >= -J: dn = sp.factor(coef[m-1])
        # normalised: with N_m the norm of e_m, element = coef * sqrt(N_target/N_source); symmetric result expected -b
        N = lambda k: sp.binomial(2*J, J+k)*sp.gamma(D+k-1)*sp.gamma(D-k-1)
        if m+1 <= J:
            norm_el = coef[m+1]*sp.sqrt(sp.gammasimp(N(m)/N(m+1)))   # element <m+1|W|m> between normalised states
            bm = sp.sqrt((C - m*(m+1))*(T - m*(m+1)))
            val = sp.N((norm_el/bm).subs(D, sp.Rational(1013, 3)))
            assert abs(val + 1) < 1e-12, (J, m, val)   # equals -b_m
    print(f"J={J}: diagonal T+C-2m^2 OK; tridiagonal; |normalised off-diagonal| = b_m; sample up-coefficient(m=0):",
          sp.factor(sp.Poly(W(a**J*b**J), a, b).coeff_monomial(a**(J+1)*b**(J-1))) if J > 0 else '-')
