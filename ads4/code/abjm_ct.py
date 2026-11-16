r"""c_T of the ABJM matrix integral: exact identity and two-sided bounds (ads4/ABJM_THEOREM.md, section 'Stress-tensor
normalisation').  With E the expectation in the undeformed gas and w_ij = sech^2((x_i-x_j)/2k):
   Z(m)/Z(0) = sech^N(pi m/2) E[ prod_{i<j} (1 + sinh^2(pi m/2) w_ij)^{-1} ]      (real mass m; m = 2 i lam gives the tilt)
   c_T := (64/pi^2) d^2F/dm^2 |_0 = (16/pi^2) Var(sum x_i) = 16 N + 32 E[sum_{i<j} w_ij].
Checks: N = 2 by quadrature (three independent evaluations of c_T); then the sandwich bounds
   [F_gas^L(v_m) - F^U(0)]/sinh^2(pi m/2) <= (c_T-16N)/32 <= [F^U(0) - F^L(lam) + N log cos(pi lam)]/sin^2(pi lam)."""
import mpmath as mp
mp.mp.dps = 15
k = 1
dens = lambda x, y: mp.tanh((x-y)/(2*k))**2/(4*mp.cosh(x/2)*mp.cosh(y/2))
I = lambda f: mp.quad(lambda x, y: f(x, y)*dens(x, y), [-50, 0, 50], [-50, 0, 50])
Z0 = I(lambda x, y: 1); var = I(lambda x, y: (x+y)**2)/Z0; Ew = I(lambda x, y: mp.sech((x-y)/(2*k))**2)/Z0
m = mp.mpf('0.7'); lhs = I(lambda x, y: mp.cos(m*(x+y)/2))/Z0
rhs = mp.sech(mp.pi*m/2)**2*I(lambda x, y: 1/(1+mp.sinh(mp.pi*m/2)**2*mp.sech((x-y)/(2*k))**2))/Z0
print("N=2,k=1: identity at m=0.7:  E[e^{imX/2}] =", mp.nstr(lhs, 10), "  sech^N E[prod] =", mp.nstr(rhs, 10))
print("          c_T = (16/pi^2) Var(X) =", mp.nstr(16/mp.pi**2*var, 9), "   16N + 32 E[w] =", mp.nstr(32+32*Ew, 9), "   (literature: 112/3 = 37.3333)")
mp.mp.dps = 20
def gas_bounds(N, vfun, c, m1, sp, sm, b0, zeta):
    """two-sided bounds on -log Z for the gas with pair potential vfun (int = c, first moment m1), trap slopes sp, sm"""
    S = 1/sp+1/sm; smax = max(sp, sm)
    def lower(theta, a):
        c1 = c - 2*mp.quad(vfun, [0, a])
        return (mp.mpf(2)/3)*mp.sqrt(2*theta*c1/S)*N**mp.mpf(1.5) - N*(2*smax*theta*m1/c1+vfun(a)/2) + mp.loggamma(N+1) + N*mp.log(zeta) - N*mp.log(S/(1-theta))
    lo = max(lower(1-t, a) for t in (1/mp.sqrt(N), mp.mpf(N)**(-0.75), mp.mpf(1)/N) for a in (1/mp.sqrt(N), mp.mpf(N)**(-0.75), mp.mpf(1)/N))
    mu = mp.sqrt(2*c*N/S)
    up = mp.loggamma(N+1) + N*mp.log(zeta) + (mp.mpf(2)/3)*mu*N + N*b0 - N*mp.log(c*N/mu)
    return lo, up
def ct_bounds(N):
    zeta = 4*mp.pi*k; v0 = lambda u: -2*mp.log(mp.tanh(abs(u)/(2*k))); c0 = mp.pi**2*k; m10 = 7*mp.zeta(3)*k*k/2
    lo0, up0 = gas_bounds(N, v0, c0, m10, mp.mpf(1)/2, mp.mpf(1)/2, mp.log(2), zeta)
    best_up, best_lo = mp.inf, -mp.inf
    for lam in (mp.mpf(N)**(-0.125)/3, mp.mpf(N)**(-0.125)/6, mp.mpf('0.2'), mp.mpf('0.1')):
        if lam >= 0.5: continue
        loL, _ = gas_bounds(N, v0, c0, m10, mp.mpf(1)/2-lam, mp.mpf(1)/2+lam, mp.log(2), zeta)
        best_up = min(best_up, (up0 - loL + N*mp.log(mp.cos(mp.pi*lam)))/mp.sin(mp.pi*lam)**2)
    for mm in (mp.mpf(N)**(-0.125), mp.mpf(N)**(-0.125)/2, mp.mpf('0.5'), mp.mpf('0.25')):
        C2 = mp.cosh(mp.pi*mm/2)**2; vm = lambda u: mp.log(1+C2/mp.sinh(abs(u)/(2*k))**2)
        cm = mp.pi**2*k*(1+mm**2); m1m = mp.quad(lambda u: u*vm(u), [0, 1, 10, 60])
        loM, _ = gas_bounds(N, vm, cm, m1m, mp.mpf(1)/2, mp.mpf(1)/2, mp.log(2), zeta)
        best_lo = max(best_lo, (loM - up0)/mp.sinh(mp.pi*mm/2)**2)
    return 16*N+32*max(best_lo, 0), min(16*N+32*best_up, 16*N*N)
print("int v_m check (m=0.5): ", mp.nstr(2*mp.quad(lambda u: mp.log(1+mp.cosh(mp.pi/4)**2/mp.sinh(u/2)**2), [0, 1, 10, 80]), 12), " vs pi^2 k (1+m^2) =", mp.nstr(mp.pi**2*1.25, 12))
for N in (10**4, 10**6, 10**8, 10**10):
    lo, up = ct_bounds(N); lead = 64*mp.sqrt(2*k)/(3*mp.pi)*mp.mpf(N)**1.5
    print(f"N={N:.0e}: c_T lower/lead = {mp.nstr(lo/lead, 5)}   upper/lead = {mp.nstr(up/lead, 5)}")
