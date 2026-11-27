"""Proposition NQ: the heavy-scalar family F = (g2 I + g3 J)^2/(2 M^2)  (a, b, c) = (g2^2/2M^2, g2 g3/M^2, g3^2/2M^2).
Checks: Gram/causality cone b^2 = 4ac (saturated), thermal-causality K(q) = (2 g2 - 3 q g3)^2/4M^2 >= 0, planar 4pi eta/s - 1 and
hyperbolic Taylor data a_k (from THEOREM_H table) at FIXED four-point data a while g3 varies."""
import sympy as sp
g2, g3, M, q, x, n = sp.symbols('g2 g3 M q x n', real=True)
a, b, c = g2**2/(2*M**2), g2*g3/M**2, g3**2/(2*M**2)
print('cone 4ac - b^2 =', sp.simplify(4*a*c - b**2))
print('thermal causality K(q) =', sp.factor(sp.expand(2*a - 3*b*q + sp.Rational(9, 2)*c*q**2)))
planar = sp.expand(864*a - 1944*b - 6480*c)
print('planar 4pi eta/s - 1 =', sp.factor(planar), ' = (1/M^2)(432 g2^2 - 1944 g2 g3 - 3240 g3^2)')
gmax = sp.solve(sp.diff(planar, g3), g3)[0]
print('  max over g3 at g3 =', sp.simplify(gmax), ' value =', sp.simplify(planar.subs(g3, gmax)), ' ; unbounded below as |g3| grows')
# hyperbolic V(x) per operator (THEOREM_H table), Taylor coefficients at n=1
V = {'I2': 1344*x - 96*x**2 - 384*x**3, 'IJ': -2880*x**2 - 2016*x**3 + 2952*x**4, 'J2': -14976*x**3 - 2448*x**4 + 10944*x**5}
Vfam = sp.expand(a*V['I2'] + b*V['IJ'] + c*V['J2'])
xn = 1 - 9*n**2/(1 + sp.sqrt(1 + 3*n**2))**2
co = [sp.simplify(sp.diff(Vfam.subs(x, xn), n, k).subs(n, 1)/sp.factorial(k)) for k in range(4)]
print('hyperbolic Taylor a_k (k=0..3) in units 1/M^2:', [sp.factor(sp.expand(ck*M**2)) for ck in co])
print('  => a_0 = 0, a_1 = -672 g2^2 depend on g2 only (4-point); g3 enters first at a_2 (5-point) and a_3 (6-point)')
