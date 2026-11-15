r"""Checks for the rigorous N^{3/2} bounds on the ABJM S^3 partition function (ads4/ABJM_THEOREM.md).
Z(N,k) = (1/N!) (4 pi k)^{-N} \int d^N x  prod_i 1/(2cosh(x_i/2))  prod_{i<j} tanh^2((x_i-x_j)/(2k)).
(1) constants: c = \int v = pi^2 k,  m1 = \int_0^inf u v(u) du = 7 zeta(3) k^2 / 2,  v(u) = -2 log tanh(|u|/2k);
(2) normalisation: Z(1,k) = 1/(4k), Z(2,1) = 1/(16 pi) by direct quadrature;
(3) the explicit two-sided bounds, evaluated and compared with exact values at small N (k = 1)."""
import mpmath as mp
mp.mp.dps = 30
def v(u, k): return -2*mp.log(mp.tanh(abs(u)/(2*k)))
for k in (1, 2, 3):
    c = 2*mp.quad(lambda u: v(u, k), [0, 1, 10, mp.inf]); m1 = mp.quad(lambda u: u*v(u, k), [0, 1, 10, mp.inf])
    print(f"k={k}: int v = {mp.nstr(c, 12)} (pi^2 k = {mp.nstr(mp.pi**2*k, 12)});  int_0^inf u v = {mp.nstr(m1, 12)} (7 zeta(3) k^2/2 = {mp.nstr(7*mp.zeta(3)*k*k/2, 12)})")
Z1 = mp.quad(lambda x: 1/(2*mp.cosh(x/2)), [-mp.inf, 0, mp.inf])/(4*mp.pi)
Z2 = mp.quad(lambda x, y: mp.tanh((x-y)/2)**2/(4*mp.cosh(x/2)*mp.cosh(y/2)), [-40, 0, 40], [-40, 0, 40])/(2*(4*mp.pi)**2)
print("Z(1,1) =", mp.nstr(Z1, 12), " exact 1/4;   Z(2,1) =", mp.nstr(Z2, 10), " exact 1/(16 pi) =", mp.nstr(1/(16*mp.pi), 10))

def bounds(N, k):
    """returns (lower, upper) bounds on F = -log Z from the theorem, optimising the free parameters numerically"""
    c = mp.pi**2*k; m1 = 7*mp.zeta(3)*k*k/2
    def lower(theta, a):                       # F >= ...
        if not (0 < theta < 1): return -mp.inf
        c1 = c - 2*mp.quad(lambda u: v(u, k), [0, a])
        if c1 <= 0: return -mp.inf
        kappa = theta*m1/c1
        return (mp.mpf(2)/3)*mp.sqrt(theta*c1/2)*N**mp.mpf(1.5) - N*(kappa + v(a, k)/2) + mp.loggamma(N+1) - N*mp.log(1/((1-theta)*mp.pi*k))
    best = max(lower(1-t, a) for t in (mp.mpf(1)/N, 2/mp.sqrt(N), 1/mp.sqrt(N), mp.mpf('0.5'), mp.mpf('0.2')) for a in (mp.mpf(1)/N, 1/mp.sqrt(N), mp.mpf('0.3'), mp.mpf(1), mp.mpf(2)*k))
    mu = mp.sqrt(c*N/2)
    upper = mp.loggamma(N+1) + N*mp.log(4*mp.pi*k) + (mp.mpf(2)/3)*mu*N + N*mp.log(2) - (mp.mpf(N)/2)*mp.log(2*c*N)
    return best, upper
exact = {1: mp.mpf(1)/4, 2: 1/(16*mp.pi), 3: (mp.pi-3)/(64*mp.pi), 4: (10-mp.pi**2)/(1024*mp.pi**2)}
print("\nk = 1:   N   lower bound <= F=-log Z <= upper bound      leading term (pi sqrt(2k)/3) N^1.5")
for N in (1, 2, 3, 4):
    lo, up = bounds(N, 1); F = -mp.log(exact[N])
    print(f"   {N:5d}  {mp.nstr(lo, 8):>12} <= {mp.nstr(F, 8):>10} <= {mp.nstr(up, 8):>12}    {mp.nstr(mp.pi*mp.sqrt(2)/3*N**1.5, 8)}   ok={bool(lo <= F <= up)}")
for N in (10, 100, 1000, 10**4, 10**6):
    lo, up = bounds(N, 1); lead = mp.pi*mp.sqrt(2)/3*mp.mpf(N)**1.5
    print(f"   {N:7d}  lower/lead = {mp.nstr(lo/lead, 6)}   upper/lead = {mp.nstr(up/lead, 6)}")
