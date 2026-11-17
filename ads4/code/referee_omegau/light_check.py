"""Referee: independent k(0), conditional bound at supported n, contact and |g3|^2 leverages (M02), pointwise identity (sympy)."""
import json, numpy as np, sympy as sp, os
from fractions import Fraction as Fr
from scipy.integrate import quad
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust', 'functionals_certified.json')))['M02']
P = np.polynomial.Polynomial
P2 = P([float(Fr(c)) for c in d['psi2_exact']]); P3 = P([float(Fr(c)) for c in d['psi3_exact']]); a = float(Fr(d['a_exact']))
p = P([0, 1]); I = lambda f: f.integ()(1) - f.integ()(0)
phi = {2: P2, 3: P3 - 4*p**2*P2, 4: 2*p**4*P2 - p**2*P3, 5: -3*p**6*P2 + 2*p**4*P3}
k0 = -(I(phi[2]*2*p**6) + I(phi[3]*2*p**4) + I(phi[4]*4*p**2) - 4*a)   # |g3|^2 coefficient / (pi G), sign of -B'
print('k(0) =', k0)
for Dg, n in ((20, 1.6), (50, 1.45), (100, 1.4), (100, 1.27), (100, 1.15), (300, 1.4), (1000, 1.35)):
    Lam = np.sqrt((Dg - 2.)*(Dg - 1.)); L = Lam/n
    X = n**8*16*quad(lambda t: P2(t)/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]/k0
    print(f'Dg={Dg} n={n}: |g3|^2 Lam^8 <= {X:.4g}, sqrt(t4^2+aT^2) <= {240*np.sqrt(X)/Lam**4:.3g}')
print('g4 lev k2,k3:', 4*I(phi[2]*p**2)/k0, 2*I(phi[3])/k0)
print('|g3|^2 lev (top/pi) k2,k3,k4:', 2*I(phi[2]*p**4)/k0, 2*I(phi[3]*p**2)/k0, 4*I(phi[4])/k0, ' (u5 prints these x pi)')
y = sp.symbols('y'); h = sp.Rational(1, 2)
for s in (1, -1):
    c2 = (2-y)/(1-y)**2 - s*4*y/(1-y)**2 + 2*y**2*(2-y)/(1-y)**3 - s*3*y**3/(1-y)**3 + h*y**4*(2-y)/(1-y)**4 - s*h*y**5/(1-y)**4
    c3 = s/(1-y)**2 - y*(2-y)/(1-y)**3 + s*2*y**2/(1-y)**3 - h*y**3*(2-y)/(1-y)**4 + s*h*y**4/(1-y)**4
    f = 1 if s == 1 else (1-y)**4
    print('sigma', s, ': psi2 coeff', sp.factor(sp.simplify(c2*f)), '; psi3 coeff / u', sp.factor(sp.simplify(c3*f)))
