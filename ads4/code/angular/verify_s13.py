"""ANGULAR_ROUTE §13 checks: Q_N(0) = -6/((N+2)(N+3)); truncated telescoping sum_{k<=N} q_{k+2} a_k (Aa)_k = (T-b(b-1))|a0|^2 - 3 B_N
with B_N = X_N |a_N|^2/((N+1) D_N); counterexample integral (1/6) int x^-3 e^{-2/x} = 1/24."""
import numpy as np, sympy as sp
from fractions import Fraction as F
for N in range(0, 60):
    M = N + 2; Q = F(-3) + 2*sum(F(6, (m-1)*m*(m+1)) for m in range(2, M+1))
    assert Q == F(-6, (N+2)*(N+3))
print("Q_N(0) = -6/((N+2)(N+3)) exact, N < 60")
rng = np.random.default_rng(11); worst = 0
for _ in range(500):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b+3, 120)); D = J + 2 + 10**rng.uniform(-1, 3)
    T = (D-1)*(D-2); C = J*(J+1); n = J-b+1; a = np.zeros(n); a[0] = 1.0
    for k in range(n-1):
        s = b+k; a[k+1] = a[k]*(k+4)/(k+1)*np.sqrt((C-s*(s+1))/(T-s*(s+1)))
    s = np.arange(b, J+1, dtype=float); A = np.diag(T+C-2*s**2)
    for i in range(n-1): A[i, i+1] = A[i+1, i] = -np.sqrt((C-s[i]*(s[i]+1))*(T-s[i]*(s[i]+1)))
    Aa = A@a; Nn = int(rng.integers(0, n-1)); Dk = lambda k: (k+1)*(k+2)*(k+3)/6
    lhs = sum(a[k]*Aa[k]/Dk(k) for k in range(Nn+1)); XN = C-(b+Nn)*(b+Nn+1)
    rhs = (T-b*(b-1))*a[0]**2 - 3*XN*a[Nn]**2/((Nn+1)*Dk(Nn))
    worst = max(worst, abs(lhs-rhs)/((T-b*(b-1))*a[0]**2))
print("truncated telescoping with B_N (factor 3 in this normalisation): max rel err", worst)
x = sp.symbols('x', positive=True); print("(1/6) int x^-3 e^{-2/x} dx =", sp.integrate(x**-3*sp.exp(-2/x), (x, 0, sp.oo))/6)
