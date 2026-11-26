# r5_lightbound.py -- the rate bound of note sec 3 and sec 6 at the actual M, by direct quadrature, and the crossover
# with the bound the paper already has at every gap >= 1000 (Theorem FP at 1000 + the monotonicity Remark of sec:allgap.A).
import mpmath as mp
mp.mp.dps = 40
def run(M):
    m = M+6
    rho = [mp.mpf(j)/(100*m*m) for j in range(1, m+1)]
    Rprof = lambda z: 1/mp.fprod([1-r*z for r in rho])
    eps = mp.mpf(1)/(M+1)**8
    w = lambda z: (1+z)*(-mp.log(z))-2*(1-z)          # Q_1 in the z-variable, times 2/(1-z)^2 * (1-z)^3 / (1-z)
    pts = [0, mp.mpf(1)/(10*M), mp.mpf(1)/M, mp.mpf(3)/M, mp.mpf(10)/M, mp.mpf(30)/M, mp.mpf(100)/M, 1]
    I3 = mp.quad(lambda z: z**3*(1-z)**M*Rprof(z)*w(z), pts)
    I1 = mp.quad(lambda z: z*(1-z)**M*Rprof(z)*w(z), pts)
    L = mp.log(mp.mpf(M+4)/3)
    bound3 = mp.mpf(100)/99*(1+mp.mpf(4)/(M+5))*6*L/(M+1)**4
    exact_beta = mp.mpf(100)/99*(mp.beta(4, M+1)*(mp.harmonic(M+4)-mp.harmonic(3))+mp.beta(5, M+1)*(mp.harmonic(M+5)-mp.harmonic(4)))
    tau_ns = mp.sqrt(12800*I3)
    # strict: K Lam^6 = 100 F2 + 36 F3 exactly
    s1 = mp.fsum(rho)
    # [z^2] Rprof = sum_{i<=j} rho_i rho_j = (s1^2 + s2)/2
    s2 = mp.fsum(r*r for r in rho)
    A1 = -M+s1
    A2 = mp.mpf(M*(M-1))/2-M*s1+(s1**2+s2)/2
    KL6 = 100*eps*A1+36*(1+eps*A2)
    tau_s = mp.sqrt(460800*(I3+eps*I1)/KL6)
    claim = 280*mp.sqrt(L)/(M+1)**2
    print('M=%d  I3=%s  bound=%s (ratio %s; Beta form %s)' % (M, mp.nstr(I3, 10), mp.nstr(bound3, 10), mp.nstr(I3/bound3, 6), mp.nstr(exact_beta, 10)))
    print('      eps*I1=%s  K*Lam^6(strict)=%s' % (mp.nstr(eps*I1, 6), mp.nstr(KL6, 15)))
    print('      tau bound: non-strict sqrt(12800 I)=%s  strict=%s  claimed 280 sqrt(log)/(M+1)^2=%s  (ratios %s, %s)' %
          (mp.nstr(tau_ns, 8), mp.nstr(tau_s, 8), mp.nstr(claim, 8), mp.nstr(tau_ns/claim, 6), mp.nstr(tau_s/claim, 6)))
    return claim
for M in [998, 2000, 10**4, 10**5]:
    run(M)
# crossover with tau_0(1000)=2.8914e-8 (certified upper end), valid at every Delta_g>=1000 by the monotonicity Remark
f = lambda Dg: 280*mp.sqrt(mp.log((Dg+2)/3))/(Dg-2)**2-mp.mpf('2.8914e-8')
x = mp.findroot(f, 2e5)
print('RP bound 280 sqrt(log((Dg+2)/3))/(Dg-2)^2 equals tau_0(1000)=2.8914e-8 at Delta_g =', mp.nstr(x, 8))
for Dg in [1000, 1e4, 1e5, 2e5, 1e6]:
    Dg = mp.mpf(Dg)
    tstar = 160/(9*(Dg-2))/mp.sqrt(1-49/(2*Dg+9)**2)
    print('  Delta_g=%-8s RP=%s  tau_*=%s  tau_0(1000)=2.8914e-8' % (mp.nstr(Dg, 7), mp.nstr(280*mp.sqrt(mp.log((Dg+2)/3))/(Dg-2)**2, 6), mp.nstr(tstar, 6)))
