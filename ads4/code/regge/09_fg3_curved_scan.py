r"""Check 9 (conditional, same model as 08): scan over the heavy mass m (units of the gap) of the +- action of F_g3 with
scalar kernels, for M = 30, 100, 300 and beta = 0.1 ... 1.  Printed: the minimum over m in [1.5, 6600] of m^4 M^2 F^(+-),
its m -> infinity value 2 M^2 int psi_2 P, and M^2 int psi_3 P."""
import mpmath as mp
mp.mp.dps = 15
P = lambda nu, be: mp.re(mp.hyp2f1(0.5+1j*nu, 0.5-1j*nu, 1, -mp.sinh(be/2)**2))
psi2 = lambda p: (1-p)**3*p*(65+155*p+47*p*p)
psi3 = lambda p: 2*(1-p)**3*(81-1927*p+4601*p*p)
def parts(M, be):
    n = max(4, int(M*be/3) + 4); pts = mp.linspace(0, 1, n)
    I2a = mp.quad(lambda p: psi2(p)*P(M*p, be), pts); I2b = mp.quad(lambda p: psi2(p)*p*p*P(M*p, be), pts)
    I3 = mp.quad(lambda p: psi3(p)*P(M*p, be), pts)
    dP = -mp.diff(lambda x: P(x, be), 0, 2)/2
    return I2a, I2b, I3, dP, P(0, be)
def F(M, m, pr, sgn):
    I2a, I2b, I3, dP, P0 = pr
    return (2*m*m*I2a - I2b + sgn*I3)/m**6 + 5*(2*m*m*M*M*dP - P0)/m**10
for M in (30, 100, 300):
    for be in map(mp.mpf, ('0.1', '0.2', '0.3', '0.5', '1.0')):
        pr = parts(M, be)
        best = min(((F(M, mp.mpf(m), pr, -1)*mp.mpf(m)**4*M*M, m) for m in [1.5*1.08**k for k in range(110)]), key=lambda t: t[0])
        lim = 2*M*M*pr[0]
        print(f"M={M} beta={mp.nstr(be,2)} (M beta = {mp.nstr(M*be,3)}):  min over m>=1.5 of m^4 M^2 F^(+-) = {mp.nstr(best[0],5)} at m = {mp.nstr(best[1],4)};  m->inf value {mp.nstr(lim,5)};  M^2 I3 = {mp.nstr(M*M*pr[2],5)}")
