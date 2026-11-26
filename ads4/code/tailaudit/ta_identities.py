"""
tailaudit/ta_identities.py [DV, mpmath 60 digits]: re-check TAIL_CERT §1/§2 identities at NON-INTEGER twist a (tc_exact.py tests
integer a only), both channels, complex beta.
 (E2) V_b on powers, (E3) three-term recurrence;
 (TP) Theorem P representation S(beta) = kappa int_0^oo q^{beta-1} g(q) dq, g(q) = sum_{K<=n} Bin(M,p)(K) C(N-K,a)/C(N,a), p = q/(1+q);
 (V)  the V-representation: for real q (p < 1)  g(q) = E_V[(1 - p/V)_+^a], V ~ Beta(N+1, n+s)   (exact for every real a > 0; the
      doc's "g = A - T, T = sum_{K>N}" is the integer-a special case, where A = E_V[(1-p/V)^a] is a polynomial);
 (LL) Lemma Loc's representation L(F) = E_B[[z^n](1-z)^{-a-1}F(Bz)]/C(N,a), B ~ Beta(n+s, a)  (checked on F = (1-z)^{-beta}).
Output: cache/ta_identities.txt
"""
import os, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60
out = []
def pr(s):
    print(s, flush=True); out.append(s)

def U(n, a, b):
    s = 2*b + 1; N = n + a
    u = [mp.mpf(1)]
    for k in range(n):
        u.append(u[-1]*(n - k)*(n + s + k)/((N - k)*(N + s + k)))
    return u

def S(beta, n, a, b):
    u = U(n, a, b); tot = 0; p = mp.mpf(1)
    for k in range(n + 1):
        tot += u[k]*p; p *= (beta + k)/(k + 1)
    return tot

def binom(x, k):
    return mp.gamma(x + 1)/(mp.gamma(k + 1)*mp.gamma(x - k + 1))

worst = {}
for (n, a, b) in ((7, mp.mpf('3.37'), 0), (11, mp.mpf('5.5'), 4), (20, mp.mpf('12.25'), 0), (15, mp.mpf('17.9'), 4)):
    s = 2*b + 1; N = n + a; T = N*(N + s) + b*(b + 1); M = N + n + s
    for beta in (mp.mpc('0.25', '3.1'), mp.mpc('2.5', '0'), mp.mpc('1.7', '-0.4'), mp.mpc('-2.3', '1.9')):
        # E2
        u = U(n, a, b) + [mp.mpf(0)]
        V = []
        for k in range(n + 2):
            um = (b + k - 1)*(b + k)
            V.append(u[k]*(1 + ((k + 3)**2/mp.mpf(2) - um)/T) + (k*k/(2*T)*u[k - 1] if k > 0 else 0))
        lhs = 0; p = mp.mpf(1)
        for k in range(n + 2):
            lhs += V[k]*p; p *= (beta + k)/(k + 1)
        A1 = beta*(beta + 1)/2 + (4 - 2*b)*beta
        A0 = -beta*(beta + 1)/2 - (mp.mpf(7)/2 - 2*b)*beta + mp.mpf(9)/2 - b*(b - 1)
        rhs = S(beta, n, a, b) + (A1*S(beta + 1, n, a, b) + A0*S(beta, n, a, b))/T
        worst['E2'] = max(worst.get('E2', 0), abs(lhs - rhs)/abs(rhs))
        # E3
        l3 = (N*(N + s) - beta*(beta - s))*S(beta - 1, n, a, b)
        r3 = (a*(N + n + s) + beta*(s - 2*beta - 1))*S(beta, n, a, b) + beta*(beta + 1)*S(beta + 1, n, a, b)
        worst['E3'] = max(worst.get('E3', 0), abs(l3 - r3)/abs(r3))
        if beta.real > 0:
            # Theorem P representation
            kap = mp.gamma(N + s)/(mp.gamma(beta)*mp.gamma(N + s - beta))
            cf = [binom(M, K)*binom(N - K, a)/binom(N, a) for K in range(n + 1)]
            g = lambda q: sum(cf[K]*(q/(1 + q))**K*(1/(1 + q))**(M - K) for K in range(n + 1))
            I = kap*mp.quad(lambda q: q**(beta - 1)*g(q), [0, 0.5, 2, 10, mp.inf])
            worst['TP'] = max(worst.get('TP', 0), abs(I - S(beta, n, a, b))/abs(S(beta, n, a, b)))
    # V-representation at real q
    Bnorm = mp.beta(N + 1, n + s)
    for q in (mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('3')):
        p = q/(1 + q)
        cf = [binom(M, K)*binom(N - K, a)/binom(N, a) for K in range(n + 1)]
        gq = sum(cf[K]*p**K*(1 - p)**(M - K) for K in range(n + 1))
        ev = mp.quad(lambda v: v**N*(1 - v)**(n + s - 1)*(1 - p/v)**a, [p, (p + 1)/2, 1])/Bnorm
        worst['V'] = max(worst.get('V', 0), abs(gq - ev)/abs(gq))
    # Lemma Loc representation on F = (1-z)^{-beta}, beta = 1.3+0.7i
    beta = mp.mpc('1.3', '0.7')
    # [z^n](1-z)^{-a-1}(1-Bz)^{-beta} = sum_k C(n-k+a, a) (beta)_k B^k / k!  ;  E B^k = (n+s)_k/(N+s)_k
    tot = 0
    for k in range(n + 1):
        tot += binom(n - k + a, a)*mp.rf(beta, k)/mp.factorial(k)*mp.rf(n + s, k)/mp.rf(N + s, k)
    worst['LL'] = max(worst.get('LL', 0), abs(tot/binom(N, a) - S(beta, n, a, b))/abs(S(beta, n, a, b)))
for k, v in worst.items():
    pr(f'{k}: max relative deviation over non-integer-a states = {mp.nstr(v, 3)}')
open(os.path.join(HERE, 'cache', 'ta_identities.txt'), 'w').write('\n'.join(out) + '\n')
