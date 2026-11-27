"""Theorem H numerics: analytic structure of the Einstein hyperbolic family in the Renyi index n = T0/T, and the
Cauchy/Taylor reconstruction of the planar endpoint n = 0 from vacuum Taylor data at n = 1.
Einstein AdS4, R = L = 1:  T = (3 r_h - 1/r_h)/4pi,  T0 = 1/2pi,  n = T0/T = 2/(3 r_h - 1/r_h),
r_h(n) = (1 + sqrt(1 + 3 n^2))/(3n),  x(n) = 1 - 1/r_h^2 (horizon Weyl scale; x(1) = 0, x(0) = 1).
Branch points of x(n): n = +-i/sqrt3  =>  Taylor radius about n = 1 is rho = sqrt(4/3) = 1.1547 > 1 = distance to n = 0."""
import mpmath as mp
mp.mp.dps = 40
def x_of_n(n):
    rh = (1 + mp.sqrt(1 + 3*n**2))/(3*n)
    return 1 - 1/rh**2
# 1/rh is analytic at n=0: 1/rh = 3n/(1+sqrt(1+3n^2))
def x_of_n_reg(n):
    inv = 3*n/(1 + mp.sqrt(1 + 3*n**2))
    return 1 - inv**2
print('x(1) =', x_of_n_reg(1), '  x(0) =', x_of_n_reg(0), '  x(0.5) =', x_of_n_reg(mp.mpf(0.5)))
rho = abs(1 - 1j/mp.sqrt(3))
print('Taylor radius about n=1: |1 - i/sqrt3| = %s ; ratio 1/rho = %s' % (mp.nstr(rho, 10), mp.nstr(1/rho, 10)))
# Taylor coefficients of x(n) and of a sample first-order deviation delta(x) = x^k at n=1, and partial sums at n=0
co = mp.taylor(x_of_n_reg, 1, 60)
print('|c_k| rho^k for k=10,20,40,60:', [mp.nstr(abs(co[k])*rho**k, 5) for k in (10, 20, 40, 60)])
for K in (5, 10, 20, 40, 60):
    S = sum(co[k]*(-1)**k for k in range(K + 1))
    print('K=%2d  partial sum at n=0: %s   error %s   (rho^-K = %s)' % (K, mp.nstr(S, 12), mp.nstr(abs(S - 1), 3), mp.nstr(rho**-K, 3)))
# Theorem H rate table: |V(0)| <= min_K [ eps*sum_{k<K} C^k + M rho^{1-K}/(rho-1) ]
def bound(eps, M, C, rho):
    best = mp.inf
    for K in range(1, 2000):
        b = eps*sum(mp.mpf(C)**k for k in range(K)) + M*rho**(1 - K)/(rho - 1)
        best = min(best, b)
    return best
print('\nTheorem H bound (M = 1): power law in Delta_g with eps = Dg^-4; fitted exponent = -dlog(bound)/dlog(Dg)')
for (C, rr) in ((1.0, rho), (1.1, rho), (1.5, rho), (1.0, mp.mpf(2)), (2.0, mp.mpf(2))):
    vals = []
    for Dg in (1e2, 1e4, 1e8):
        vals.append(bound(mp.mpf(Dg)**-4, 1, C, rr))
    expo = -(mp.log(vals[2]) - mp.log(vals[1]))/(mp.log(1e8) - mp.log(1e4))
    pred = 4*mp.log(rr)/mp.log(rr*max(C, 1))
    print('C=%.2f rho=%.4f : bound(Dg=1e2,1e4,1e8) = %s ; exponent %s (prediction 4 log rho/log(rho*max(C,1)) = %s)' % (
        C, rr, [mp.nstr(v, 3) for v in vals], mp.nstr(expo, 4), mp.nstr(pred, 4)))
