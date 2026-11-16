r"""Tilted ABJM integral  Z(N,k;lam) = (1/N!)(4 pi k)^{-N} \int prod e^{lam x_i}/(2cosh(x_i/2)) prod_{i<j} tanh^2((x_i-x_j)/2k),
|lam| < 1/2  (trial R-charges (1/2+lam, 1/2-lam, 1/2, 1/2); imaginary mass m = -2 i lam).
Checks of the bounds (U_lam), (L_lam) of ads4/ABJM_THEOREM.md, section 'Extension':
  F(lam) = -log Z = f0 sqrt(1-4 lam^2) N^{3/2} + O(N log N),  and exact facts Z(1,k;lam) = 1/(4k cos(pi lam)),
  Z(lam) >= Z(0) (F maximal at lam = 0)."""
import mpmath as mp
mp.mp.dps = 20
def v(u, k): return -2*mp.log(mp.tanh(abs(u)/(2*k)))
def bounds(N, k, lam):
    c = mp.pi**2*k; m1 = 7*mp.zeta(3)*k*k/2; sp, sm = mp.mpf(1)/2-lam, mp.mpf(1)/2+lam; S = 1/sp+1/sm; smax = max(sp, sm)
    def lower(theta, a):
        c1 = c - 2*mp.quad(lambda u: v(u, k), [0, a]); kappa = 2*smax*theta*m1/c1
        return (mp.mpf(2)/3)*mp.sqrt(2*theta*c1/S)*N**mp.mpf(1.5) - N*(kappa+v(a, k)/2) + mp.loggamma(N+1) - N*mp.log(S/((1-theta)*4*mp.pi*k))
    lo = max(lower(1-t, a) for t in (mp.mpf('0.5'), mp.mpf('0.2'), 1/mp.sqrt(N+1), mp.mpf(1)/(N+1)) for a in (mp.mpf('0.3'), mp.mpf(1), 2*mp.mpf(k), 1/mp.sqrt(N+1), mp.mpf(1)/(N+1)))
    mu = mp.sqrt(2*c*N/S)
    up = mp.loggamma(N+1) + N*mp.log(4*mp.pi*k) + (mp.mpf(2)/3)*mu*N + N*mp.log(2) - N*mp.log(c*N/mu)   # entropy: f <= mu/(cN)
    return lo, up
def Z2(k, lam):
    f = lambda x, y: mp.exp(lam*(x+y))*mp.tanh((x-y)/(2*k))**2/(4*mp.cosh(x/2)*mp.cosh(y/2))
    return mp.quad(f, [-60, 0, 60], [-60, 0, 60])/(2*(4*mp.pi*k)**2)
print("N=1: Z = 1/(4k cos(pi lam)); bounds must bracket F")
for k, lam in ((1, mp.mpf('0.3')), (2, mp.mpf('-0.2')), (1, mp.mpf('0.45'))):
    F = -mp.log(1/(4*k*mp.cos(mp.pi*lam))); lo, up = bounds(1, k, lam); print(f"  k={k} lam={mp.nstr(lam,3)}: {mp.nstr(lo,6)} <= {mp.nstr(F,8)} <= {mp.nstr(up,6)}  ok={bool(lo<=F<=up)}")
mp.mp.dps = 12
print("N=2 (quadrature):")
for k, lam in ((1, mp.mpf(0)), (1, mp.mpf('0.3'))):
    Z = Z2(k, lam); F = -mp.log(Z); lo, up = bounds(2, k, lam); print(f"  k={k} lam={mp.nstr(lam,3)}: Z={mp.nstr(Z,10)}  {mp.nstr(lo,6)} <= {mp.nstr(F,8)} <= {mp.nstr(up,6)}  ok={bool(lo<=F<=up)}")
print("ratio of bounds to f0 sqrt(1-4 lam^2) N^1.5, k=1:")
for lam in (mp.mpf(0), mp.mpf('0.25'), mp.mpf('0.4')):
    for N in (10**3, 10**5, 10**7):
        lo, up = bounds(N, 1, lam); lead = mp.pi*mp.sqrt(2)/3*mp.sqrt(1-4*lam**2)*mp.mpf(N)**1.5
        print(f"  lam={mp.nstr(lam,3)} N={N:.0e}: lower/lead={mp.nstr(lo/lead,6)}  upper/lead={mp.nstr(up/lead,6)}")
