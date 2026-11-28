# Term 1: -1/2 sum (2J+1) log Pi_vac(w_n,J) = -S(beta/r) (+ linear in beta), S(x) = sum_N sigma_2(N) e^{-Nx}/N
#   check S(x) = zeta(3)/x^2 + (1/12) log x + zeta'(-1) - x^2/2880 + ...
from mpmath import mp, mpf, zeta, log, exp, nsum, inf, diff
mp.dps = 30
def S(x):
    # sum_m (1/m) 1/(4 sinh^2(m x/2))
    return nsum(lambda m: 1/(m*4*mp.sinh(m*x/2)**2), [1, inf])
zp = diff(zeta, -1)
for x in [mpf('0.2'), mpf('0.1'), mpf('0.05')]:
    v = S(x) - zeta(3)/x**2 - log(x)/12 - zp
    print(x, v, v/x**2, ' expect -1/2880 =', -mpf(1)/2880)
print("zeta'(-1) =", zp)
