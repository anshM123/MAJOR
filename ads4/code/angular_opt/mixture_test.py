"""Test (i): m_k = prod_{i<k}(C-C_i)/(T-C_i) = E[C(Y,k)]/C(J,k) with P(Y=y)>=0 (exact inversion), for real N=Delta-2.
   Test (ii): sup_y |Q_n(y;a,b,M)| for Hahn polynomials."""
import mpmath as mp
mp.mp.dps = 80
def m_seq(b, J, D):
    C = J*(J+1); T = (D-1)*(D-2); M = J-b
    m = [mp.mpf(1)]
    for k in range(M):
        y = (b+k)*(b+k+1); m.append(m[-1]*(C-y)/(T-y))
    return m  # length M+1
def Ydist(m):
    M = len(m)-1
    P = []
    for y in range(M+1):
        s = mp.mpf(0)
        for k in range(y, M+1):
            s += (-1)**(k-y)*mp.binomial(M-y, k-y)*m[k]
        P.append(mp.binomial(M, y)*s)
    return P
worst = 1
import random
for J in list(range(0, 40)) + [50, 80, 120]:
    for b in (0, 4):
        if J - b < 0: continue
        for D in [J+2.001, J+2.5, J+3.3, J+5.7, J+10.1, J+30.5, J+100.25, J+1000.5, 2*J+7.3]:
            P = Ydist(m_seq(b, J, mp.mpf(D)))
            mn = min(P)
            if mn < worst: worst = mn; print('new min', J, b, D, mp.nstr(mn, 5), 'sum', mp.nstr(sum(P), 10))
print('overall min P(Y=y):', worst)
