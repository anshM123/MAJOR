r"""Uniform curved-regime positivity for explicit functional weights against the scalar AdS4 kernel.
  G_M[f](beta) := \int_0^1 f(p) P_{-1/2+i M p}(cosh beta) dp
               = (sqrt2/pi) \int_beta^infty K_{f,M}(t) (cosh t - cosh beta)^{-1/2} dt,
  K_{f,M}(t)  := \int_0^1 f(p) coth(pi M p) sin(M p t) dp        (Mehler's formula with coth).
If K_{f,M}(t) >= 0 for all t > 0 then G_M[f](beta) > 0 for all beta > 0.
Weights: f_g3(p) = (1-p)^3 p (65 + 155 p + 47 p^2)   (B_2 weight of the graviton functional, arXiv:2201.06602 (3.3a), as read
by the audit), f_sc(p) = p (1-p)^2 (1310 + 1540 p - 10318 p^2 - 4697 p^3 + 15680 p^4)  (scalar functional, 2106.10274)."""
import mpmath as mp
mp.mp.dps = 20
f_g3 = lambda p: (1-p)**3*p*(65+155*p+47*p*p)
f_sc = lambda p: p*(1-p)**2*(1310+1540*p-10318*p**2-4697*p**3+15680*p**4)
def K(f, M, t):
    n = int(M*t/mp.pi)+2
    pts = [mp.mpf(i)/n for i in range(n+1)]
    return mp.quad(lambda p: f(p)*mp.coth(mp.pi*M*p)*mp.sin(M*p*t) if p > 0 else f(mp.mpf('1e-30'))*0, pts)
def fs(f, tau):                       # sine transform \int_0^1 f(p) sin(p tau) dp
    n = int(tau/mp.pi)+2; pts = [mp.mpf(i)/n for i in range(n+1)]
    return mp.quad(lambda p: f(p)*mp.sin(p*tau), pts)
for name, f in (("f_g3", f_g3), ("f_sc", f_sc)):
    print(name, " min over tau in (0,200] of tau^3 * sine-transform:", mp.nstr(min(fs(f, mp.mpf(k)/4)*(mp.mpf(k)/4)**3 for k in range(1, 801)), 6),
          "  min of p-weight f on [0,1]:", mp.nstr(min(f(mp.mpf(i)/200) for i in range(201)), 5))
    for M in (1, 2, 5, 20, 80):
        ts = [mp.mpf(k)/(4*M) for k in range(1, 160)] + [mp.mpf(k)/8 for k in range(1, 161)]
        vals = [(K(f, M, t), t) for t in ts]
        mn = min(vals, key=lambda v: v[0])
        print(f"   M={M:>3}: min_t K = {mp.nstr(mn[0], 5)} at t = {mp.nstr(mn[1], 4)}   (sign: {'>= 0 on grid' if mn[0] >= 0 else 'NEGATIVE somewhere'})")
