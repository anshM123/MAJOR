"""EXACT rational check of the mixture lemma (Lemma M) at random rational twists:
   m_k^(b) = E[C(Y,k)]/C(M,k), k=0..M, with P(Y=y) >= 0, for any integer M with M-1 < d = twist-2.
   Y | U ~ BetaBinomial(U; a, d), U ~ generalized Hypergeometric(N; Jb good; M draws),
   b=0: Jb=J, N=Delta-2, a=J+1 ;  b=4: Jb=J-4, N=Delta-6, a=J+5."""
from fractions import Fraction as Fr
import random
def gbin(x, m):
    r = Fr(1)
    for i in range(m): r *= (x - i)/Fr(i+1)
    return r
def poch(x, m):
    r = Fr(1)
    for i in range(m): r *= (x + i)
    return r
def binom(n, k):
    return gbin(Fr(n), k) if 0 <= k <= n else Fr(0)
def ydist(b, J, D, M):
    if b == 0: Jb, N, a = J, D-2, Fr(J+1)
    else: Jb, N, a = J-4, D-6, Fr(J+5)
    d = N - Jb
    PU = [gbin(Fr(Jb), u)*gbin(d, M-u)/gbin(N, M) for u in range(M+1)]
    P = [Fr(0)]*(M+1)
    for u in range(M+1):
        if PU[u] == 0: continue
        for y in range(u+1):
            P[y] += PU[u]*binom(u, y)*poch(a, y)*poch(d, u-y)/poch(a+d, u)
    return P, PU
def mseq(b, J, D, K):
    C = J*(J+1); T = (D-1)*(D-2); m = [Fr(1)]
    for k in range(K):
        y = (b+k)*(b+k+1)
        m.append(m[-1]*(C-y)/(T-y) if k+1 <= J-b else Fr(0))
    return m
if __name__ == '__main__':
    random.seed(7); ok = 0
    for trial in range(40):
        b = random.choice([0, 4]); J = random.randint(4, 120)
        if b == 0 and J % 2: J += 1
        h = Fr(random.randint(700, 2500), 100); D = J + h
        Mmax = int(h - 2) if (h-2).denominator != 1 else int(h-2)
        M = random.randint(1, Mmax)
        P, PU = ydist(b, J, D, M); m = mseq(b, J, D, M)
        assert min(P) >= 0 and min(PU) >= 0
        for k in range(M+1):
            assert sum(P[y]*binom(y, k) for y in range(M+1))/binom(M, k) == m[k], (b, J, h, M, k)
        ok += 1
    print('Lemma M verified exactly in', ok, 'random cases (b, J<=120, rational twist in [7,25], 1<=M<=twist-2)')
