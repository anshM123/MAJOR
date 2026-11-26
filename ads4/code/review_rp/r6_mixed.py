# r6_mixed.py -- checks of MIXED_PROFILE_THEOREM.md at the actual size M=998 (high precision, large spins):
#  (1) H_b/(w^3 L_b(F)) >= 1/25 on heavy states incl. J up to 1500 (the coordinator's exact test reached smaller J)
#  (2) tolerance of b=4 positivity to the k=3 heavy normalisation: c* = min_states T*L[(1+D/T)F]/L(G)
#      (b=4 positivity survives if the k=3 heavy term is multiplied by any c < c*)
#  (3) K Lambda^6 = 100F_2+36F_3+8G_2 exactly for the chosen R, and its split into k=2 and k=3 parts
#  (4) I by quadrature, tau^2 <= 460800 I/(K Lambda^6) versus the claimed 420 sqrt(log Dg)/(Dg-2)^{5/2}
#  (5) robustness of K: only top light coefficients enter the k=3 part (G_0=G_1=0); the k=2 part depends on the
#      lower (AdS) coefficients of Q but is O(M) against the k=3 O(M^2)
import sys
import mpmath as mp
mp.mp.dps = 800
M = 998; D = M+6
rho = [(mp.mpf(j)+mp.mpf(1)/2)/(100*(D+1)**2) for j in range(1, D+1)]   # the t=1/2 member of the family of sec 3
eps = mp.mpf(1)/(M+1)**8; A = mp.mpf(7)*M*M/10
Jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
N = Jmax+3
Rc = [mp.mpf(1)]+[mp.mpf(0)]*N
for r in rho:
    for i in range(1, N+1): Rc[i] += r*Rc[i-1]
Bc = [(-1)**l*mp.binomial(M, l) if l <= M else mp.mpf(0) for l in range(N+1)]
A0 = [mp.fsum(Bc[i]*Rc[kk-i] for i in range(kk+1)) for kk in range(N+1)]
F = [mp.mpf(0)]*(N+1); G = [mp.mpf(0)]*(N+1)
for kk in range(N+1):
    if kk >= 2: F[kk] += A0[kk-2]; G[kk] = A*A0[kk-2]
    if kk >= 1: F[kk] += eps*A0[kk-1]
def Ulist(b, J, Delta, kmax):
    T = (Delta-1)*(Delta-2); C = J*(J+1); U = [mp.mpf(1)]
    for j in range(kmax):
        u = (b+j)*(b+j+1); U.append(U[-1]*(C-u)/(T-u) if j < J-b else mp.mpf(0))
    return U
def pieces(b, J, Delta):
    T = (Delta-1)*(Delta-2); n = J-b; U = Ulist(b, J, Delta, n+2)
    dressed = 0; LF = 0; LG = 0
    for kk in range(0, n+2):
        um1 = (b+kk-1)*(b+kk) if kk >= 1 else (b-1)*b
        v = U[kk]*(1+(mp.mpf((kk+3)**2)/2-um1)/T)
        if kk >= 1: v += mp.mpf(kk*kk)/(2*T)*U[kk-1]
        dressed += v*F[kk]; LF += U[kk]*F[kk]; LG += U[kk]*G[kk]
    return T, dressed, LF, LG
minratio = mp.inf; cstar = mp.inf; arg1 = arg2 = None
for b in (0, 4):
    for J in sorted(set([b+1, b+2, b+3, 6, 10, 30, 100, 300, 700, 1000, Jmax])):
        if J-b < 1 or J > Jmax: continue
        for a in [mp.mpf(998), mp.mpf('998.37'), mp.mpf(1005), mp.mpf(1500), mp.mpf(3000), mp.mpf(10)**4, mp.mpf(10)**6]:
            T, dr, LF, LG = pieces(b, J, J+a+2)
            sig = 1 if b == 0 else -1
            H = dr+sig*LG/T
            r = H/LF
            if r < minratio: minratio = r; arg1 = (b, J, a)
            tol = T*dr/LG if LG != 0 else mp.inf
            if b == 4 and tol < cstar: cstar = tol; arg2 = (J, a)
            print('b=%d J=%5d a=%-8s H/L(F)=%-12s T*L[(1+D/T)F]/L(G)=%s' % (b, J, mp.nstr(a, 7), mp.nstr(r, 6), mp.nstr(tol, 6)))
            sys.stdout.flush()
print('(1) min H/(w^3 L(F)) =', mp.nstr(minratio, 6), 'at', arg1, ' (claimed >= 1/25 = 0.04)')
print('(2) tolerance c* on the k=3 heavy normalisation (b=4):', mp.nstr(cstar, 6), 'at', arg2)
mp.mp.dps = 40
s1 = mp.fsum(rho); s2 = mp.fsum(r*r for r in rho)
r1 = s1; r2 = (s1**2+s2)/2
A1 = -M+r1; A2 = mp.mpf(M*(M-1))/2-M*r1+r2
F2 = 1+eps*A1; F3 = A1+eps*A2; G2 = A
k2 = 100*F2+36*F3; k3 = 8*G2
print('(3) K Lambda^6 = %s  (k=2 part %s, k=3 part %s);  K Lambda^6/M^2 = %s ; claimed lower bound 11/2 M^2 = %s' %
      (mp.nstr(k2+k3, 12), mp.nstr(k2, 8), mp.nstr(k3, 8), mp.nstr((k2+k3)/M**2, 6), mp.nstr(mp.mpf(11)/2*M*M, 8)))
Rf = lambda z: 1/mp.fprod([1-r*z for r in rho])
w = lambda z: (1+z)*(-mp.log(z))-2*(1-z)
pts = [0, mp.mpf(1)/(10*M), mp.mpf(1)/M, mp.mpf(3)/M, mp.mpf(10)/M, mp.mpf(30)/M, mp.mpf(100)/M, 1]
I = mp.quad(lambda z: (z*z+eps*z)*(1-z)**M*Rf(z)*w(z), pts)
for Dg in [mp.mpf(1000), mp.mpf('1000.99')]:
    Mg = int(mp.floor(Dg-2))
    if Mg != M: continue
    tau = mp.sqrt(460800*I/(k2+k3))
    claim = 420*mp.sqrt(mp.log(Dg))/(Dg-2)**mp.mpf(2.5)
    print('(4) Dg=%s  sqrt(460800 I/(K Lam^6)) = %s   claimed 420 sqrt(log Dg)/(Dg-2)^2.5 = %s   ratio %s' %
          (mp.nstr(Dg, 7), mp.nstr(tau, 8), mp.nstr(claim, 8), mp.nstr(tau/claim, 6)))
Ibound = mp.mpf(100)/99*mp.log(1000)/(M+1)**3*(2*(1+mp.mpf(3)/(M+4))+(1+mp.mpf(2)/(M+3))/(M+1)**7)
print('    I = %s  analytic bound (Dg=1000) = %s  ratio %s' % (mp.nstr(I, 8), mp.nstr(Ibound, 8), mp.nstr(I/Ibound, 5)))
# crossover with tau_0(1000)=2.8914e-8
f = lambda Dg: 420*mp.sqrt(mp.log(Dg))/(Dg-2)**mp.mpf(2.5)-mp.mpf('2.8914e-8')
print('    420-bound equals tau_0(1000)=2.8914e-8 at Delta_g =', mp.nstr(mp.findroot(f, (1e4, 1e5), solver='bisect'), 8))
