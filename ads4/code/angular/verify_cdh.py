"""Verify: celestial three-term recurrence (SP1) solutions = normalised continuous dual Hahn 3F2,
and reproduce the K_B2A4 k=3 polynomial. Exact sympy arithmetic. Run: python3 verify_cdh.py"""
import sympy as sp
x = sp.symbols('x')  # x = nu^2

def rec(d, l, M):
    k = {d: sp.Integer(0), d + 1: sp.Integer(1)}
    for m in range(d + 1, M):
        k[m + 1] = sp.expand(-((x + sp.Rational(1, 4) + 2*d*(d + 1) - 2*m*(m + l))*k[m]
                               + (m - 1 - d)*(m + l - 1 - d)*k[m - 1]) / ((m + 1 + d)*(m + l + 1 + d)))
    return k

def cdh(n, a, b, c):  # 3F2(-n, a+i nu, a-i nu; a+b, a+c; 1)
    return sp.expand(sum(sp.rf(-n, j)*sp.Mul(*[(a + i)**2 + x for i in range(j)])
                         / (sp.rf(a + b, j)*sp.rf(a + c, j)*sp.factorial(j)) for j in range(n + 1)))

k = rec(1, 0, 8)
assert sp.expand(k[3] + (4*x - 15)/64) == 0
assert sp.expand(400*k[4] - (x**2 - sp.Rational(35, 2)*x + sp.Rational(569, 16))) == 0   # K_B2A4 k=3 factor
assert sp.expand(921600*k[5] - (-64*x**3 + 2896*x**2 - 26956*x + 39159)) == 0
for d, l in [(1, 0), (0, 0), (1, 4), (0, 4), (1, 2), (2, 3)]:
    k = rec(d, l, d + 10)
    a, b, c = 2*d + sp.Rational(3, 2), sp.Rational(1, 2), l + sp.Rational(1, 2)
    assert all(sp.expand(k[j + d + 1] - cdh(j, a, b, c)) == 0 for j in range(9)), (d, l)
print("OK: SP1 recurrence == continuous dual Hahn (a,b,c)=(2d+3/2,1/2,l+1/2), j<=8; k4/k2 = K_B2A4 polynomial")
