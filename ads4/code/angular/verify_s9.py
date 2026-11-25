"""ANGULAR_ROUTE §9 checks: telescoping identity on ladder vectors (random weights), q_m endpoint identity, Ytilde values, 151250 constant."""
import numpy as np, sympy as sp
rng = np.random.default_rng(3)
def ladder(J, D, b):
    T = (D-1)*(D-2); C = J*(J+1); n = J - b + 1; a = np.zeros(n); a[0] = 1.0
    for k in range(n-1):
        s = b + k; a[k+1] = a[k]*(k+4)/(k+1)*np.sqrt((C - s*(s+1))/(T - s*(s+1)))
    s = np.arange(b, J+1, dtype=float); A = np.diag(T + C - 2*s**2)
    for i in range(n-1): A[i, i+1] = A[i+1, i] = -np.sqrt((C - s[i]*(s[i]+1))*(T - s[i]*(s[i]+1)))
    return a, A, T, C
worst = 0
for _ in range(300):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b+1, 60)); D = J + 3 + rng.random()*80
    a, A, T, C = ladder(J, D, b); n = len(a); w = rng.random(n) + 0.1
    lhs = np.sum(w*a*(A@a)); rho = np.array([(k+4)/(k+1) for k in range(n)]); X = np.array([C - (b+k)*(b+k+1) for k in range(n)])
    rhs = w[0]*(T - b*(b-1))*a[0]**2 + sum((rho[k]-1)*(rho[k]*w[k+1] - w[k])*X[k]*a[k]**2 for k in range(n-1))
    worst = max(worst, abs(lhs - rhs)/max(abs(lhs), 1))
print("telescoping identity on ladder vectors, max rel err:", worst)
q = lambda m: sp.Rational(6, (m-1)*m*(m+1)) if abs(m) >= 2 else (sp.Integer(-3) if m == 0 else sp.Integer(0))
qq = lambda m: q(abs(m))
worst = 0
for _ in range(200):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b+1, 60)); D = J + 3 + rng.random()*80
    a, A, T, C = ladder(J, D, b); w = np.array([float(qq(k+2)) for k in range(len(a))])
    worst = max(worst, abs(np.sum(w*a*(A@a)) - (T - b*(b-1))*a[0]**2)/((T - b*(b-1))*a[0]**2))
print("q-weighted endpoint identity, max rel err:", worst)
Yt = lambda m: qq(m) + (1 - sp.Rational(m, 2))*qq(m-1) + (1 + sp.Rational(m, 2))*qq(m+1)
print("Ytilde_0, _1, _2:", Yt(0), Yt(1), Yt(2), "; Ytilde_m = 0 for 3<=m<200:", all(Yt(m) == 0 for m in range(3, 200)))
m = sp.symbols('m'); print("general m>=3 symbolic:", sp.simplify(6/((m-1)*m*(m+1)) + (1 - m/2)*6/((m-2)*(m-1)*m) + (1 + m/2)*6/(m*(m+1)*(m+2))))
print("constant (T-12)/3/W lower bound:", float(sp.Rational(1, 300)*sp.Rational(29999, 30000)**3*sp.Rational(729, 16)*998**2), "; bound 4*sqrt(3/151250) =", float(4*sp.sqrt(sp.Rational(3, 151250))))
