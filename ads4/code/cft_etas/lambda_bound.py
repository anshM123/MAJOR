"""Rigorous-bound checks for Theorem C constants.
(a) beta_N <= Gamma(1+2g) exp(-2g psi(N+1)) <= Gamma(1+2g) N^(-2g)   [log-convexity of Gamma; psi(N+1) > log N]
(b) 1/|B_k(k)| <= (N+1)^(2g) C(2k-1,k) C(N+k,2k)   and   S_N := sum_k C(2k-1,k) C(N+k,2k) <= (3+2sqrt2)^N  (checked N<=400)
"""
import mpmath as mp
from math import comb
mp.mp.dps = 40
a = 3 + 2*mp.sqrt(2)
worst = 0
for N in range(1, 401):
    S = sum(comb(2*k - 1, k)*comb(N + k, 2*k) for k in range(1, N + 1))
    worst = max(worst, mp.mpf(S)/a**N)
print('max_{N<=400} S_N/(3+2sqrt2)^N =', mp.nstr(worst, 8))
ok = True
for g in (0.25, 0.5, 1.0, 1.5, 1.75):
    for N in (1, 2, 5, 10, 40, 100):
        b = mp.gamma(N + 1)*mp.gamma(1 + 2*g)/mp.gamma(N + 1 + 2*g)
        ok &= b <= mp.gamma(1 + 2*g)*mp.exp(-2*g*mp.digamma(N + 1)) <= mp.gamma(1 + 2*g)*mp.mpf(N)**(-2*g)
        for k in (1, max(1, N//2), N):
            inv = mp.fprod([(k + j + 2*g)/abs(k - j) for j in range(1, N + 1) if j != k])
            ok &= inv <= (N + 1)**(2*g)*comb(2*k - 1, k)*comb(N + k, 2*k)
print('inequalities (a),(b) hold on all sampled (g,N,k):', ok)
# (c) closed form: a_N := sum_{k=0}^N C(N+k,2k) 4^k obeys a_N = 6 a_{N-1} - a_{N-2}, a_0=1, a_1=5,
#     so a_N = A(3+2sqrt2)^N + B(3-2sqrt2)^N with A = (1+1/sqrt2)/2; hence S_N <= a_N/2 <= 0.427(3+2sqrt2)^N + 1/2.
an = [sum(comb(N + k, 2*k)*4**k for k in range(N + 1)) for N in range(60)]
print('recursion a_N = 6a_{N-1}-a_{N-2} holds for N<60:', all(an[N] == 6*an[N-1] - an[N-2] for N in range(2, 60)))
A = (1 + 1/mp.sqrt(2))/2; B = 1 - A
print('closed form matches:', all(abs(A*a**N + B*(3 - 2*mp.sqrt(2))**N - an[N]) < 1e-20*an[N] for N in range(60)), ' A =', mp.nstr(A, 6))
