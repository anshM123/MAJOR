"""
tailcert/tc_exact.py (TAIL_CERT.md §1) [exact rationals + DV]: the identities the tail argument rests on.
 (E1) ladder: U_k = (-n)_k (n+s)_k / ((-N)_k (N+s)_k), N = n + a, s = 2b+1, T = (N(N+s) + b(b+1)), C = J(J+1) (fp_common.ladder);
      S(beta) := L((1-z)^{-beta}) = sum_k U_k (beta)_k/k!  = 3F2(beta,-n,n+s;-N,N+s;1).
 (E2) first moment: V_{b,k} (FULLPOS §32) acting on (1-z)^{-beta} equals S(beta) + tau[A1 S(beta+1) + A0 S(beta)],
      A1 = beta(beta+1)/2 + (4-2b) beta, A0 = -beta(beta+1)/2 - (7/2-2b) beta + 9/2 - b(b-1), tau = 1/T   (exact, rational beta).
 (E3) three-term recurrence (derived from the annihilator sum_k U_k[(T-u_{k-1})F_{k-1} - (C-u_k)F_k] = 0, u_k = (b+k)(b+k+1)):
      [N(N+s) - beta(beta-s)] S(beta-1) = [a(N+n+s) + beta(s-2beta-1)] S(beta) + beta(beta+1) S(beta+1)       (exact, rational beta)
 (E4) branch formula: (1-z)^{-3} P_{-1/2+i nu}((1+z)/(1-z)) = 2 Re[B(nu) sum_j c_j(nu) (1-z)^{j-5/2-i nu}],
      B = Gamma(2i nu)/Gamma(1/2+i nu)^2, c_j = ((1/2-i nu)_j)^2/((1-2i nu)_j j!)    (mpmath, |1-z| < 1).
Output: cache/tc_exact.txt
"""
import os, sys, random
from fractions import Fraction as Fr
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s):
    print(s, flush=True); out.append(s)

def poch(x, k):
    p = Fr(1)
    for i in range(k):
        p *= (x + i)
    return p

def ladder(n, a, b):
    s = 2*b + 1; N = n + a
    U = [Fr(1)]
    for k in range(n):
        U.append(U[-1]*Fr((n - k)*(n + s + k), (N - k)*(N + s + k)))
    return U

def S(beta, n, a, b):
    U = ladder(n, a, b)
    return sum(U[k]*poch(beta, k)/poch(1, k) for k in range(n + 1))

random.seed(1)
ok = True
for trial in range(30):
    b = random.choice((0, 4)); n = random.randint(1, 25); a = random.randint(1, 30); s = 2*b + 1; N = n + a
    J = n + b; C = J*(J + 1); T = N*(N + s) + b*(b + 1)
    # check U against fp_common's product form prod (C-u_j)/(T-u_j)
    U = ladder(n, a, b)
    U2 = [Fr(1)]
    for j in range(n):
        u = (b + j)*(b + j + 1); U2.append(U2[-1]*Fr(C - u, T - u))
    ok &= (U == U2)
    beta = Fr(random.randint(-40, 40), random.randint(1, 9))
    # (E2): V_{b,k} = U_k[1 + ((k+3)^2/2 - u_{k-1})/T] + (k^2/(2T)) U_{k-1}, k <= n+1 ; acting on F_k = (beta)_k/k!
    Ue = U + [Fr(0)]
    V = []
    for k in range(n + 2):
        um = (b + k - 1)*(b + k)
        V.append(Ue[k]*(1 + (Fr((k + 3)**2, 2) - um)/T) + (Fr(k*k, 2*T)*Ue[k - 1] if k > 0 else 0))
    lhs = sum(V[k]*poch(beta, k)/poch(1, k) for k in range(n + 2))
    A1 = beta*(beta + 1)/2 + (4 - 2*b)*beta
    A0 = -beta*(beta + 1)/2 - (Fr(7, 2) - 2*b)*beta + Fr(9, 2) - b*(b - 1)
    rhs = S(beta, n, a, b) + Fr(1, T)*(A1*S(beta + 1, n, a, b) + A0*S(beta, n, a, b))
    ok &= (lhs == rhs)
    # (E3)
    l3 = (N*(N + s) - beta*(beta - s))*S(beta - 1, n, a, b)
    r3 = (a*(N + n + s) + beta*(s - 2*beta - 1))*S(beta, n, a, b) + beta*(beta + 1)*S(beta + 1, n, a, b)
    ok &= (l3 == r3)
pr(f'(E1) ladder product form, (E2) first-moment operator on powers, (E3) recurrence: 30 random exact-rational states/betas: {ok}')
# (E4)
mp.mp.dps = 30
for nu, z in ((mp.mpf('0.7'), mp.mpf('0.6')), (mp.mpf(5), mp.mpf('0.9')), (mp.mpf('20.3'), mp.mpf('0.97'))):
    X = (1 + z)/(1 - z)
    lhs = (1 - z)**-3*mp.legenp(-0.5 + 1j*nu, 0, X, type=3).real
    Bn = mp.gamma(2j*nu)/mp.gamma(0.5 + 1j*nu)**2
    a_ = mp.mpf(1)/2 - 1j*nu
    rhs = 2*mp.re(Bn*(1 - z)**(-mp.mpf(5)/2 - 1j*nu)*mp.hyp2f1(a_, a_, 1 - 2j*nu, 1 - z))
    pr(f'(E4) nu={nu}, z={z}: lhs {mp.nstr(lhs, 15)}  rhs {mp.nstr(rhs, 15)}')
open(os.path.join(HERE, 'cache', 'tc_exact.txt'), 'w').write('\n'.join(out) + '\n')
