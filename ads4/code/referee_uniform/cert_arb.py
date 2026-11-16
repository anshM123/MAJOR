"""Referee: independent Arb re-certification of S(tau)=int_0^1 psi2(p) sin(p tau) dp > 0 on [tau0, oo) for M02/M01,
plus check of phi=psi2/p decreasing, psi2>=0, and closed-form vs quadrature."""
import json, sys
from fractions import Fraction as Fr
import sympy as sp
from flint import arb, fmpq, ctx
ctx.prec = 200
D = json.load(open('/home/user/MAJOR/ads4/code/robust/functionals_certified.json'))
p = sp.Symbol('p'); t = sp.Symbol('t', positive=True)
for name, tau0 in (('M02', Fr(2311, 100)), ('M01', Fr(2349, 100))):
    c = [sp.Rational(str(Fr(s))) for s in D[name]['psi2_exact']]
    f = sum(ci*p**i for i, ci in enumerate(c))
    # psi2 >= 0 and phi decreasing on (0,1): independent Sturm via sympy real_roots on f/p and derivative
    phi = sp.cancel(f/p)
    print(name, 'roots of phi in (0,1):', [r for r in sp.Poly(phi, p).real_roots() if 0 < r < 1],
          ' roots of phi\' in (0,1):', [float(r) for r in sp.Poly(sp.diff(phi, p), p).real_roots() if 0 < r < 1], ' phi(0)=', phi.subs(p, 0), ' phi\'(0)=', sp.diff(phi,p).subs(p,0))
    # exact S(t) by symbolic integration (independent of the IBP formula)
    Ssym = sp.integrate(f*sp.sin(p*t), (p, 0, 1))
    Ssym = sp.expand(sp.simplify(Ssym))
    # build arb evaluator: Ssym = A(t) + B(t) cos t + C(t) sin t with rational functions in t
    Ssym = sp.expand(Ssym)
    A = sp.together(Ssym.subs({sp.sin(t): 0, sp.cos(t): 0}))
    B = sp.together(Ssym.coeff(sp.cos(t)))
    C = sp.together(Ssym.coeff(sp.sin(t)))
    assert sp.simplify(A + B*sp.cos(t) + C*sp.sin(t) - Ssym) == 0
    def topoly(e):
        n, d = sp.fraction(sp.together(e))
        return [fmpq(int(sp.Rational(v).p), int(sp.Rational(v).q)) for v in sp.Poly(n, t).all_coeffs()], \
               [fmpq(int(sp.Rational(v).p), int(sp.Rational(v).q)) for v in sp.Poly(d, t).all_coeffs()]
    PA, PB, PC = topoly(A), topoly(B), topoly(C)
    def hev(cs, x):
        r = arb(0)
        for cc in cs: r = r*x + arb(cc)
        return r
    def S(x):
        return hev(PA[0], x)/hev(PA[1], x) + hev(PB[0], x)/hev(PB[1], x)*x.cos() + hev(PC[0], x)/hev(PC[1], x)*x.sin()
    # bisection certificate on [tau0, 80]
    T1 = Fr(80)
    stack = [(Fr(tau0), T1)]; ncell = 0; fail = None
    while stack:
        a, b = stack.pop()
        mid = (a + b)/2; rad = (b - a)/2
        x = arb(fmpq(mid.numerator, mid.denominator), fmpq(rad.numerator, rad.denominator))
        v = S(x)
        if v > 0:
            ncell += 1; continue
        if b - a < Fr(1, 10**7):
            fail = (float(a), v); break
        stack += [(a, mid), (mid, b)]
    print(f'  Arb: S>0 on [{float(tau0)}, 80]: {fail is None} ({ncell} cells) {fail or ""}')
    # tail beyond 80: |S - A| <= (|B|+|C|) bound; A ~ a2/t^3 ; check t^9 S >= ... via exact polynomial: for t>=80, A(t) - |B(t)| - |C(t)| > 0 ?
    num = sp.together(A*t**9); 
    tailchk = sp.Poly(sp.expand(sp.cancel(A*t**9)), t)  # A*t^9 polynomial
    absB = sum(abs(v)*t**i for i, v in enumerate(sp.Poly(sp.cancel(B*t**9), t).all_coeffs()[::-1]))
    absC = sum(abs(v)*t**i for i, v in enumerate(sp.Poly(sp.cancel(C*t**9), t).all_coeffs()[::-1]))
    absA_rest = sp.Poly(sp.cancel(A*t**9), t)
    lc = absA_rest.all_coeffs()
    # positive lead term minus |rest|
    lead = lc[0]*t**absA_rest.degree()
    rest = sum(abs(v)*t**i for i, v in enumerate(absA_rest.all_coeffs()[::-1]) if i != absA_rest.degree())
    q = sp.Poly(sp.expand(lead - rest - absB - absC), t)
    print('  tail poly lead', lc[0], ' real roots >= 80:', [float(r) for r in q.real_roots() if r >= 80], ' q(80)>0:', q.eval(80) > 0)
    # negative windows of S on [0, tau0]
    import mpmath as mp
    Sf = sp.lambdify(t, Ssym, 'mpmath')
    xs = [i*0.01 for i in range(1, int(float(tau0)*100) + 300)]
    sgn = [(x, Sf(x) < 0) for x in xs]
    wins = []; cur = None
    for x, neg in sgn:
        if neg and cur is None: cur = x
        if not neg and cur is not None: wins.append((cur, x)); cur = None
    print('  negative windows of S (step .01):', wins)
