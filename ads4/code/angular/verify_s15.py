"""ANGULAR_ROUTE §15 checks: finite-harmonic obstruction int_0^1 (1-x) P_w = 0 for every finitely supported w;
endpoint weights w_k = 1/D_k give P_w = 1 (as a formal series); ladder bound (C-y)/(T-y) <= C/T."""
import sympy as sp
x = sp.symbols('x'); D = lambda k: sp.Rational((k+1)*(k+2)*(k+3), 6)
for k in range(0, 30):
    term = D(k)**2*(x**k/(k+3) - x**(k+1)/(k+1))
    assert sp.integrate((1-x)*term, (x, 0, 1)) == 0
print("int (1-x) [x^k/(k+3) - x^(k+1)/(k+1)] dx = 0 for k < 30  =>  int (1-x) P_w = 0 for finite w")
K = 40; P = 3*sum(D(k)*(x**k/(k+3) - x**(k+1)/(k+1)) for k in range(K))
print("w_k = 1/D_k, truncated at K =", K, ": P_w =", sp.factor(sp.expand(P)), " (-> 1 as K -> inf for |x|<1; the truncation leaves -(K+1)(K+2)/2 x^K... )")
C, T, y = sp.symbols('C T y', positive=True)
print("C/T - (C-y)/(T-y) =", sp.factor(C/T - (C-y)/(T-y)), " (>= 0 for T >= C, 0 <= y < C)")
