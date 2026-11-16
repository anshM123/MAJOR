"""
Radial Casimir recursion for d = 3 scalar blocks with general external dimensions, derived symbolically.

Conventions (Dolan-Osborn):  a = -Delta_12/2, b = Delta_34/2,
   Cas g = (1/2) C g,   C = Delta(Delta-3) + J(J+1),
   Cas = D_z + D_zb + (d-2) z zb/(z-zb) ((1-z) d_z - (1-zb) d_zb),   D_z = z^2(1-z) d_z^2 - (1+a+b) z^2 d_z - a b z.
Radial variables: z = 4 rho/(1+rho)^2, rho = r t, rhob = r/t, eta = (t+1/t)/2.
We multiply (Cas - C/2) by  mult = (1-r^2)(1+2 r eta+r^2)^2 (1-2 r eta+r^2)  and write the result as
   sum_terms  r^k eta^i * coeff(E, lam; a, b, C) * {1 ('I')  or  (1-eta^2) d_eta ('D')}
acting on r^E P_j(eta):  theta = r d_r -> E,  L = -(1-eta^2) d_eta^2 + 2 eta d_eta -> lam = j(j+1).
Output: terms.json = list of [k, i, type, c_E0, c_E1, c_E2, c_lam, c_C]; coefficient = c_E0 + c_E1 E + c_E2 E^2 + c_lam lam + c_C C,
each c a polynomial in a, b (string).
"""
import sympy as sp, json

r, t, r0, e0 = sp.symbols('r t r0 e0', positive=True)
a, b, C, E, lam, eta = sp.symbols('a b C E lam eta')
rho, rhob = r*t, r/t
etat = (t + 1/t)/2

def casimir_on(F):
    th = lambda G: r*sp.diff(G, r)
    dc = lambda G: t*sp.diff(G, t)
    th1 = lambda G: (th(G) + dc(G))/2          # rho d_rho
    th2 = lambda G: (th(G) - dc(G))/2          # rhob d_rhob
    def Dz(G, q, thq):
        return thq(thq(G)) - (1+q**2)/(1-q**2)*thq(G) - (a+b)*4*q/(1-q**2)*thq(G) - a*b*4*q/(1+q)**2*G
    third = (rhob*(1-rho**2)*th1(F) - rho*(1-rhob**2)*th2(F))/((rho-rhob)*(1-r**2))
    return Dz(F, rho, th1) + Dz(F, rhob, th2) + third

X, Y = sp.symbols('X Y')
tests = {'f0': sp.S(1), 'fr': X, 'fe': Y, 'frr': X**2/2, 'fee': Y**2/2, 'fre': X*Y}
mult_t = (1-r**2)*(1+2*r*etat+r**2)**2*(1-2*r*etat+r**2)
coef = {}
for name, poly in tests.items():
    Ff = poly.subs({X: r - r0, Y: etat - e0}, simultaneous=True)
    val = casimir_on(Ff).subs({r0: r, e0: etat})
    coef[name] = sp.cancel(sp.together(val*mult_t))

def to_eta(expr, extra_den=1):
    """Laurent polynomial in t, symmetric under t -> 1/t  ->  polynomial in eta (Chebyshev)."""
    expr = sp.cancel(sp.together(expr))
    num, den = sp.fraction(expr)
    # den must be a monomial in t (and r)
    P = sp.Poly(sp.expand(num), t)
    dpoly = sp.Poly(den, t)
    assert len(dpoly.terms()) == 1, den
    (dt,), dc = dpoly.terms()[0]
    res = 0
    cf = {}
    for (n,), c in P.terms():
        cf[n - dt] = sp.cancel(c/dc)
    for n, c in cf.items():
        assert sp.simplify(c - cf.get(-n, 0)) == 0, ("not symmetric", n)
        if n > 0:
            res += 2*c*sp.chebyshevt(n, eta)
        elif n == 0:
            res += c
    return sp.expand(res)

assert sp.simplify(coef['fre']) == 0
c_f0 = to_eta(coef['f0']); c_fr = to_eta(coef['fr']); c_frr = to_eta(coef['frr'])
c_fee = to_eta(coef['fee'])
# f_eta coefficient is t-antisymmetric times sinh? No: d/d eta is even; check symmetric
c_fe = to_eta(coef['fe'])
print("f_ee coefficient / (eta^2-1):", sp.factor(c_fee))
cE2 = sp.cancel(c_frr/r**2)
cE1 = sp.cancel(c_fr/r - c_frr/r**2)
cL = sp.cancel(-c_fee/(1-eta**2))
cD = sp.cancel((c_fe + 2*eta*c_fee/(1-eta**2))/(1-eta**2))
mult = (1-r**2)*(1+2*r*eta+r**2)**2*(1-2*r*eta+r**2)
c0 = sp.expand(c_f0 - C/2*mult)
for nm, e in (('cE2', cE2), ('cE1', cE1), ('cL', cL), ('cD', cD)):
    assert sp.fraction(sp.cancel(e))[1].is_number, (nm, sp.fraction(sp.cancel(e))[1])
terms = []
for typ, expr in (('I', c0 + cE1*E + cE2*E**2 + cL*lam), ('D', cD)):
    P = sp.Poly(sp.expand(expr), r, eta)
    for (k, i), cc in P.terms():
        pc = sp.Poly(sp.expand(cc), E, lam, C)
        d = {'E0': 0, 'E1': 0, 'E2': 0, 'lam': 0, 'C': 0}
        for mon, v in pc.terms():
            d[{(0,0,0): 'E0', (1,0,0): 'E1', (2,0,0): 'E2', (0,1,0): 'lam', (0,0,1): 'C'}[mon]] = v
        terms.append([int(k), int(i), typ] + [str(d[x]) for x in ('E0', 'E1', 'E2', 'lam', 'C')])
terms.sort()
json.dump(terms, open('terms.json', 'w'), indent=0)
print(len(terms), "terms; max k =", max(x[0] for x in terms), "max i =", max(x[1] for x in terms))
for x in terms: print(x)
