"""Numerical check of the mixture lemma: for twist h (d=h-2 > M-1), ladder b in {0,4},
   m_k = prod_{i<k}(C-C_{b+i})/(T-C_{b+i}) = E[C(Y,k)]/C(M,k) for k<=M, with
   Y = Bin(U, p), U ~ Hyp(N real; Jb good; M draws), p ~ Beta(Jb+1+2b'..)  (b=0: Beta(J+1,d); b=4: Beta(J+5,d))."""
import mpmath as mp, random
mp.mp.dps = 50
def gbin(x, m):  # generalized binomial C(x,m), integer m>=0
    if m < 0: return mp.mpf(0)
    r = mp.mpf(1)
    for i in range(m): r *= (x - i)/(i+1)
    return r
def ydist(b, J, D, M):
    if b == 0: Jb, N, a = J, D-2, J+1
    else: Jb, N, a = J-4, D-6, J+5       # (J'+9) = J+5
    d = N - Jb
    PU = [gbin(Jb, u)*gbin(d, M-u)/gbin(N, M) for u in range(M+1)]
    P = [mp.mpf(0)]*(M+1)
    for u in range(M+1):
        if PU[u] == 0: continue
        for y in range(u+1):
            P[y] += PU[u]*mp.binomial(u, y)*mp.beta(y+a, u-y+d)/mp.beta(a, d)
    return P, PU
def mseq(b, J, D, K):
    C = J*(J+1); T = (D-1)*(D-2); m = [mp.mpf(1)]
    for k in range(K):
        y = (b+k)*(b+k+1)
        m.append(m[-1]*(C-y)/(T-y) if k+1 <= J-b else mp.mpf(0))
    return m
random.seed(1)
worst = 0
for trial in range(60):
    b = random.choice([0, 4]); J = random.randint(4, 300)
    if b == 0 and J % 2: J += 1
    h = random.uniform(8, 40); D = mp.mpf(J) + h; M = int(mp.floor(h-2))
    P, PU = ydist(b, J, D, M)
    m = mseq(b, J, D, M)
    err = max(abs(sum(P[y]*mp.binomial(y, k) for y in range(M+1))/mp.binomial(M, k) - m[k]) for k in range(M+1))
    worst = max(worst, err)
    assert min(P) >= 0 and min(PU) >= 0
print('max moment mismatch', worst, ' (all P>=0)')
