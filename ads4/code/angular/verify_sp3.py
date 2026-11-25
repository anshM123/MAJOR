"""Checks of ANGULAR_ROUTE §3 (author claims, 25 Nov):
SP3 sum-of-squares identity and the Gershgorin-type bound; f_m = smeared d=0,l=0 CDH harmonics with w = 8nu^2/cosh(pi nu);
the z_m formula (true up to an overall sign), the ratio inequality, and opposite signs of z_m and y_m. Exact or high precision."""
import numpy as np, mpmath as mp
from fractions import Fraction as F
rng = np.random.default_rng(0); worst = 1e9
for _ in range(300):
    J = int(rng.integers(0, 40)); D = J + 2 + rng.random()*60; T = (D-2)*(D-1); C = J*(J+1); ms = np.arange(-J, J+1)
    A = np.diag(T + C - 2.0*ms**2)
    for i, m in enumerate(ms[:-1]):
        A[i, i+1] = A[i+1, i] = -np.sqrt((C - m*(m+1))*(T - m*(m+1)))
    worst = min(worst, np.linalg.eigvalsh(A).min() - (np.sqrt(T) - np.sqrt(C))**2 + 1e-9*T)
    v = rng.normal(size=2*J+1); x = C - ms[:-1]*(ms[:-1]+1); y = T - ms[:-1]*(ms[:-1]+1)
    rhs = np.sum((np.sqrt(x)*v[:-1] - np.sqrt(y)*v[1:])**2) + (T - C)*v[0]**2
    assert abs(v@A@v - rhs) < 1e-7*abs(rhs) + 1e-6
assert worst > 0; print("SP3 identity and (sqrtT - sqrtC)^2 bound: OK on 300 random (Delta, J)")
mp.mp.dps = 30
def Q(j, x): return sum(mp.rf(-j, k)*mp.fprod([(mp.mpf(3)/2 + i)**2 + x for i in range(k)])/(mp.rf(2, k)**2*mp.factorial(k)) for k in range(j+1))
def rf(a, n):
    r = F(1)
    for i in range(n): r *= a + i
    return r
f = {0: F(0)}
for m in range(1, 3003): f[m] = rf(F(1, 2), m-1)*rf(F(1), m-1)/rf(F(2), m-1)**2
for m in range(1, 7):
    num = mp.quad(lambda n: 8*n**2/mp.cosh(mp.pi*n)*Q(m-1, n**2), [0, 2, 6, mp.inf]); assert abs(num - mp.mpf(f[m].numerator)/f[m].denominator) < 1e-15
print("f_m = (1/2)_{m-1}(1)_{m-1}/(2)_{m-1}^2 matches quadrature, m<=6")
z = lambda m: -12*f[m] + (m-1)*f[m-1] - (m+1)*f[m+1]
y = lambda m: F(12-m*m, 2)*f[m] + F(13-7*m, 2)*f[m-1] + F(13+7*m, 2)*f[m+1] + F((m-2)**2, 4)*f.get(m-2, F(0)) + F((m+2)**2, 4)*f[m+2]
for m in range(2, 3000):
    assert z(m)/f[m] == -F(3*(12*m*m - 7*m - 24), 2*(m+1)*(2*m - 3))   # author's formula up to overall sign
    assert z(m+1)/z(m) > F(m-1, m+2); assert z(m)*y(m) < 0
print("z_m formula (overall sign -1 in this convention), ratio inequality, sign(z) = -sign(y): OK for 2 <= m < 3000")
