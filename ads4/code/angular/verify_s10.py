"""ANGULAR_ROUTE §10 checks: d_0, d_1 on ladders; coercivity C^2|a0|^2 <= 25/4(|(Aa)_0|^2+|(Aa)_1|^2); Q(theta) closed form;
K(g) and the bound; abstract counterexample B -> e^{-2}/6."""
import numpy as np, mpmath as mp
rng = np.random.default_rng(5)
def ladder(J, D, b, n=None):
    T = (D-1)*(D-2); C = J*(J+1); n = J-b+1 if n is None else min(n, J-b+1); a = np.zeros(n); a[0] = 1.0
    for k in range(n-1):
        s = b+k; a[k+1] = a[k]*(k+4)/(k+1)*np.sqrt((C-s*(s+1))/(T-s*(s+1)))
    s = np.arange(b, b+n, dtype=float); A = np.diag(T+C-2*s**2)
    for i in range(n-1): A[i, i+1] = A[i+1, i] = -np.sqrt((C-s[i]*(s[i]+1))*(T-s[i]*(s[i]+1)))
    return a, A, T, C
worst_d = 0; worst_c = np.inf
for _ in range(2000):
    b = int(rng.choice([0, 4])); J = int(rng.integers(b+1, 400)); tw = 2 + 10**rng.uniform(-3, 3.5); D = J + tw
    a, A, T, C = ladder(J, D, b, n=6); Aa = A@a
    if len(a) > 2:
        d0 = Aa[0]/a[0]; d1 = Aa[1]/a[1]
        worst_d = max(worst_d, abs(d0-(T-3*C+2*b*b+4*b))/abs(T+C), abs(d1-(0.75*T-1.5*C+0.75*(b+1)*(b+4)))/abs(T+C))
    worst_c = min(worst_c, 6.25*(Aa[0]**2+(Aa[1]**2 if len(a) > 1 else 0))/(C**2*a[0]**2))
print("d0,d1 formulas max rel err:", worst_d, "(rows interior: valid when ladder length >= 3)")
print("coercivity: min of 25/4(|Aa0|^2+|Aa1|^2)/(C^2|a0|^2) over 2000 random (b,J,twist>2):", worst_c, "(claim >= 1)")
mp.mp.dps = 25
q = lambda m: 6/mp.mpf((m-1)*m*(m+1)) if m >= 2 else (mp.mpf(-3) if m == 0 else mp.mpf(0))
for th in [0.3, 1.0, 2.5]:
    ser = q(0) + 2*mp.nsum(lambda m: q(int(m))*mp.cos(m*th), [2, mp.inf])
    closed = (1-mp.cos(th))*(12*mp.log(2*abs(mp.sin(th/2)))-9)
    print("Q(theta) series vs closed form at", th, mp.nstr(ser, 15), mp.nstr(closed, 15))
K = lambda g: 729/4800*(g-2)**2*(1-49/(2*g+9)**2)
print("K(1000) =", K(1000), " bound 4 sqrt(3/K) =", 4*np.sqrt(3/K(1000)), " asymptotic 160/(9g) =", 160/9/1000)
v = np.linspace(1, 1e5, 2000001); g = 1000; print("min_v (v+g-2)^6/v^4 numeric vs 729/16 (g-2)^2:", np.min((v+g-2)**6/v**4), 729/16*(g-2)**2)
