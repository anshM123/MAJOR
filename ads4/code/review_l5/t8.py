import mpmath as mp
mp.mp.dps = 30
n, a, s = 5, mp.mpf(4.5), 1
N = n + a; M = N + n + s
C0 = n*(n+s); D = N*(N+s)
U = [mp.mpf(1)]
for j in range(n): U.append(U[-1]*(C0 - j*(j+s))/(D - j*(j+s)))
def S(beta): return mp.fsum(U[k]*mp.rf(beta, k)/mp.factorial(k) for k in range(n+1))
def binom(x, y): return mp.gamma(x+1)/(mp.gamma(y+1)*mp.gamma(x-y+1))
def g(q):
    p = q/(1+q)
    return mp.fsum(binom(M, K)*p**K*(1-p)**(M-K)*binom(N-K, a)/binom(N, a) for K in range(n+1))
def kap(beta): return mp.gamma(N+s)/(mp.gamma(beta)*mp.gamma(N+s-beta))
# check Mellin rep of S
for beta in [mp.mpf(2.5), mp.mpc(2.5, 1.3)]:
    I = mp.quad(lambda q: q**(beta-1)*g(q), [0, 1, 10, mp.inf])
    print('S', S(beta), kap(beta)*I)
nu = mp.mpf(1.3)
al = mp.mpf(2.5) + 1j*nu; A_ = 0.5 - 1j*nu; B_ = 1 - 2j*nu
Bf = mp.gamma(2j*nu)/mp.gamma(0.5+1j*nu)**2
# exact K via Taylor coefficients of F_{3,nu}
g_ = [mp.mpf(1), 3 - nu**2 - 0.25]
for k in range(1, n): g_.append(((2*k*k + 6*k + 3 - nu**2 - 0.25)*g_[k] - (k + 2)**2*g_[k-1])/(k+1)**2)
K = mp.fsum(U[k]*g_[k] for k in range(n+1))
F32 = lambda q: mp.hyper([A_, A_, 1-al], [B_, N+s-al], -1/q)
Im = mp.quad(lambda q: q**(al-1)*g(q)*F32(q), [0, 0.1, 1, 10, mp.inf])
kM = Bf*kap(al)*Im
print('K exact', K, ' 2Re kappa_Mellin', 2*kM.real, kM)
# truncated j-series
tot = 0; cj = 1
for j in range(60):
    tot += cj*S(al - j); cj = cj*(A_+j)**2/((B_+j)*(j+1))
    if j in (10, 20, 30, 40, 59): print(' j-series J=', j, 2*(Bf*tot).real)
