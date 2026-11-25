"""ANGULAR_ROUTE §12 check: d_k = (Aa)_k/a_k on ladders equals 3[T-(b+k-1)(b+k)]/(k+3) - 3[C-(b+k)(b+k+1)]/(k+1) (interior k)."""
import numpy as np
rng = np.random.default_rng(7); worst = 0
for _ in range(1000):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b+3, 300)); D = J + 2 + 10**rng.uniform(-2, 3.5)
    T = (D-1)*(D-2); C = J*(J+1); n = J-b+1; a = np.zeros(n); a[0] = 1.0
    for k in range(n-1):
        s = b+k; a[k+1] = a[k]*(k+4)/(k+1)*np.sqrt((C-s*(s+1))/(T-s*(s+1)))
    s = np.arange(b, J+1, dtype=float); A = np.diag(T+C-2*s**2)
    for i in range(n-1): A[i, i+1] = A[i+1, i] = -np.sqrt((C-s[i]*(s[i]+1))*(T-s[i]*(s[i]+1)))
    Aa = A@a
    for k in range(0, min(n-1, 40)):
        if a[k] < 1e-250: break
        dk = Aa[k]/a[k]; f = 3*(T-(b+k-1)*(b+k))/(k+3) - 3*(C-(b+k)*(b+k+1))/(k+1)
        worst = max(worst, abs(dk-f)/(abs(T)+abs(C)))
print("d_k closed form on ladders: max rel err", worst)
