# r4_heavy.py -- high-precision checks of the positive-basis lemma (note 2.2) and of the full dressed action (2.3, 6)
# at the actual size M=998 (and M=1000.7 -> M=998 with fractional twists), including large spins, which the earlier
# exact-rational checks (J<=64) did not reach.  The alternating binomial sums cancel ~300 digits, so dps=800 is used.
#   (A) L_b[z^k(1-z)^m0] computed directly (alternating sum of U_{k+l}) versus U_k * sum_j p_j (a)^up_j/(N'+s')^up_j
#       (the note's E[(1-B0)^X], with E[(1-B0)^j]=(a)_j/(N'+s')_j for B0~Beta(n'+s',a)), and sign.
#   (B) V_b(F_strict) = sum_{k<=n+1} V_{b,k} F_k (fp-single, the paper's own coefficients) on a set of heavy states,
#       J up to 1500, twists a in {998 (edge, Delta_g=1000), 998.37, 1000, 2500, 10^5, 10^9}.
import sys, random
import mpmath as mp
mp.mp.dps = 800
random.seed(11)
def Ulist(b, J, Delta, kmax):
    T = (Delta-1)*(Delta-2); C = J*(J+1)
    U = [mp.mpf(1)]
    for j in range(kmax):
        if j >= J-b+0 and j >= J-b: pass
        u = (b+j)*(b+j+1)
        U.append(U[-1]*(C-u)/(T-u) if j < J-b else mp.mpf(0))
    return U
# (A) lemma identity at real size
print('(A) positive-basis lemma at M=998')
M = 998
worst = mp.mpf(1)
for trial in range(24):
    b = random.choice((0, 4))
    n = random.choice([1, 2, 5, 37, 300, 997, 998, 999, 1500, 4000])
    J = n+b
    a = mp.mpf(M) + random.choice([0, mp.mpf('0.37'), 1, 2.5, 500, 10**5])
    Delta = J+a+2
    k = random.choice([0, 1, 3, min(n, 10), n])
    m0 = random.choice([M-1, M, M+1])
    U = Ulist(b, J, Delta, k+m0+1)
    direct = mp.fsum((-1)**l*mp.binomial(m0, l)*U[k+l] for l in range(m0+1))
    np_ = n-k; Np = n+a-k; s = 2*b+1; sp_ = s+2*k
    if np_ == 0:
        pred = U[k]
    else:
        # hypergeometric pmf p_j = C(m0,j)(n')_j^down (a)^down_{m0-j}/(N')^down_{m0}
        ff = lambda x, r: mp.fprod([x-i for i in range(r)])
        den = ff(Np, m0)
        tot = 0
        for j in range(min(m0, np_)+1):
            pj = mp.binomial(m0, j)*ff(np_, j)*ff(a, m0-j)/den
            assert pj >= 0
            tot += pj*mp.rf(a, j)/mp.rf(Np+sp_, j)
        pred = U[k]*tot
    rel = abs(direct-pred)/abs(pred) if pred != 0 else abs(direct)
    print('   b=%d n=%5d a=%-8s k=%4d m0=%4d  L=%s  rel.diff=%s' % (b, n, mp.nstr(a, 7), k, m0, mp.nstr(direct, 6), mp.nstr(rel, 3)))
    if pred > 0: worst = min(worst, direct/U[k])
print('   smallest normalised value L/U_k:', mp.nstr(worst, 6))
# (B) full action
def Fstrict(N, M, eps):
    m = M+6; rho = [mp.mpf(j)/(100*m*m) for j in range(1, m+1)]
    R = [mp.mpf(1)]+[mp.mpf(0)]*N
    for r in rho:
        for i in range(1, N+1): R[i] += r*R[i-1]
    Bc = [(-1)**l*mp.binomial(M, l) if l <= M else mp.mpf(0) for l in range(N+1)]
    A = [mp.fsum(Bc[i]*R[kk-i] for i in range(kk+1)) for kk in range(N+1)]
    F = [mp.mpf(0)]*(N+1)
    for kk in range(N+1):
        if kk >= 3: F[kk] += A[kk-3]
        if kk >= 1: F[kk] += eps*A[kk-1]
    return F
Jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
print('(B) full action, building F to degree', Jmax+2); sys.stdout.flush()
eps = mp.mpf(1)/(M+1)**8
F = Fstrict(Jmax+2, M, eps)
F0 = Fstrict(Jmax+2, M, 0)
def Vb(F, b, J, Delta):
    T = (Delta-1)*(Delta-2)
    n = J-b
    U = Ulist(b, J, Delta, n+2)
    tot = 0
    for kk in range(0, n+2):
        um1 = (b+kk-1)*(b+kk) if kk >= 1 else (b-1)*b
        v = U[kk]*(1+(mp.mpf((kk+3)**2)/2-um1)/T)
        if kk >= 1: v += mp.mpf(kk*kk)/(2*T)*U[kk-1]
        tot += v*F[kk]
    return tot
res = []
for b in (0, 4):
    for J in sorted(set([b, b+1, b+2, b+3, 10, 50, 200, 500, 997, 1001, Jmax])):
        if J < b or J > Jmax: continue
        for a in [mp.mpf(998), mp.mpf('998.37'), mp.mpf(1000), mp.mpf(2500), mp.mpf(10)**5, mp.mpf(10)**9]:
            Delta = J+a+2
            v = Vb(F, b, J, Delta); v0 = Vb(F0, b, J, Delta)
            T = (Delta-1)*(Delta-2)
            res.append((b, J, a, v, v0))
            print('   b=%d J=%5d a=%-10s  T*V(F_strict)=%-14s T*V(F_nonstrict)=%s' % (b, J, mp.nstr(a, 8), mp.nstr(T*v, 6), mp.nstr(T*v0, 6)))
            sys.stdout.flush()
print('   all strict values > 0:', all(r[3] > 0 for r in res), '; all non-strict >= 0:', all(r[4] >= 0 for r in res),
      '; non-strict zero exactly for n<=1:', all((r[4] == 0) == (r[1]-r[0] <= 1) for r in res))
