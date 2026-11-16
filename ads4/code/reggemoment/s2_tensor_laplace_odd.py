"""Stage 2(b), continued: parity-odd weight-0 structure for a stress-tensor pair,
     A := eps(z1, z3, yhat) * D1 = eps(z1,z3,yhat) (z1.z3 + (z1.yhat)(z3.yhat)),
contracted with (z1, z3) = (eps_+, eps_-) and (eps_-, eps_+) (helicity preserved, helicity +2 resp. -2), and with (eps_+, eps_+).
Result:  L[A Phi_a] = C_a gamma_a |p|^{-a} * M * P,   g1 := gamma_{a-1}(nu)/gamma_a(nu)."""
import sympy as sp
from s2_tensor_laplace import *

g1 = sp.Symbol('g1')
r1 = sp.Rational(1, 2)*g1                                             # C_{a-1} gamma_{a-1} / (C_a gamma_a)
r3 = sp.Rational(1, 8)*g1*4/((a - sp.Rational(7, 2))**2 + nu**2)      # C_{a-3} gamma_{a-3} / (C_a gamma_a)
def epsD(z1, z3, f):                                                  # eps(z1, z3, d/dp) f
    return sum(sp.diff(eps3(z1, z3, v), v[m])*sum(ETA[m, n]*D(f, n) for n in range(3)) for m in range(3))
def oddTT(z1, z3):
    return epsD(z1, z3, dot(z1, z3)*r1*Psi(a-1) + zD(z1, zD(z3, r3*Psi(a-3))))
for name, (z1, z3) in {'(+2 -> +2)': (EP, EM), '(-2 -> -2)': (EM, EP), 'flip (eps_+, eps_+)': (EP, EP)}.items():
    e, P0, P1 = at_point(oddTT(z1, z3))
    print(f"  {name}:  coefficient of P: {sp.factor(sp.cancel(sp.diff(e, P0)))} ;  coefficient of P': {sp.factor(sp.cancel(sp.diff(e, P1)))}")
