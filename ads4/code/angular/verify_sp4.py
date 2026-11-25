"""Checks of ANGULAR_ROUTE §5 (author claims): ratio-difference formula; nu = i/2 coefficients; SP4 B = (ZA+AZ)/2 >= (sqrtT-sqrtC)^2 Z
on helicity sectors b = 0, 4 (random Delta, J); real-frequency weights z_m do NOT give arbitrary-amplitude positivity at large J."""
import numpy as np, sympy as sp
from scipy.linalg import eigh
m = sp.symbols('m')
f = lambda k: sp.rf(sp.Rational(1, 2), k-1)/(sp.factorial(k-1)*k**2)
zr = lambda k: 12*f(k) - (k-1)*f(k-1) + (k+1)*f(k+1)          # positive convention z = -Z
for k in range(2, 60):
    d = sp.simplify(zr(k+1)/zr(k) - sp.Rational(k-1, k+2))
    g = sp.Rational(12*k**3 - 17*k**2 + 43*k - 48, 2*(k+1)*(k+2)*(12*k**2 - 7*k - 24))
    assert sp.simplify(d - g) == 0 and d > 0, (k, d, g)
assert sp.expand((12*m**3 - 17*m**2 + 43*m - 48).subs(m, m + 2)) == 12*m**3 + 55*m**2 + 119*m + 66
print("ratio-difference formula OK (k<60), numerator at m=n+2 has positive coefficients")
# nu = i/2: f_m = 1/m solves (m+1)^2 f_{m+1} - 2 m^2 f_m + (m-1)^2 f_{m-1} = 0, and z_m = 12/m
for k in range(2, 50):
    assert (k+1)**2*sp.Rational(1, k+1) - 2*k**2*sp.Rational(1, k) + (k-1)**2*(sp.Rational(1, k-1) if k > 1 else 0) == (0 if k > 1 else 0)
    assert 12*sp.Rational(1, k) - (k-1)*(sp.Rational(1, k-1) if k > 1 else 0) + (k+1)*sp.Rational(1, k+1) == sp.Rational(12, k)
print("nu = i/2: f_m = 1/m, z_m = 12/m OK for m >= 2 (z_1 = 13)")
def sector(J, D, b, nmax=None):
    T = (D-1)*(D-2); C = J*(J+1); s = np.arange(b, J+1 if nmax is None else b+nmax, dtype=float)
    A = np.diag(T + C - 2*s**2)
    for i in range(len(s)-1):
        A[i, i+1] = A[i+1, i] = -np.sqrt((C - s[i]*(s[i]+1))*(T - s[i]*(s[i]+1)))
    return A, (np.sqrt(T) - np.sqrt(C))**2
rng = np.random.default_rng(1); worst = np.inf
for _ in range(400):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b, 300)); D = J + 2 + rng.random()*(10 + 3*J)
    A, g = sector(J, D, b); k = np.arange(A.shape[0]); Z = np.diag(12.0/(k+2))
    B = (Z@A + A@Z)/2; lam = eigh(B, Z, eigvals_only=True).min()
    worst = min(worst, (lam - g)/max(g, 1))
print("SP4: min over 400 random (b,J,Delta) of (lambda_min(B,Z) - (sqrtT-sqrtC)^2)/scale =", worst)
# real-frequency weights: generalized eigenvalue at large J can go negative (author's counterexample J=1e6, Delta=J+1000)
zf = np.array([float(zr(k+2)) for k in range(25)])
for J, D in [(10**4, 10**4 + 1000), (10**6, 10**6 + 1000)]:
    n = 25; A, g = sector(J, D, 0, n); Ar = A; Z = np.diag(zf[:n]); B = (Z@Ar + Ar@Z)/2
    print(f"real-freq weights, J={J}, Delta=J+1000, first {n} helicities: min eig(B,Z)/C =", eigh(B, Z, eigvals_only=True).min()/(J*(J+1)))
