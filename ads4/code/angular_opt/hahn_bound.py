"""Q_n(z;M) = sum_k c_k C(z,k)/C(M,k), c_k = (-n)_k (n+s)_k/(k! (s+1)_k).  Exact rationals.
Max/min over z in {0..M}."""
from fractions import Fraction as Fr
import sys
def ck(n, s):
    c = [Fr(1)]
    for k in range(n):
        c.append(c[-1]*Fr((-n+k)*(n+s+k), (k+1)*(s+1+k)))
    return c
def Qvals(n, M, s=0):
    c = ck(n, s)
    vals = []
    for z in range(M+1):
        # sum_k c_k C(z,k)/C(M,k) ; ratio C(z,k)/C(M,k) = prod (z-i)/(M-i)
        tot = Fr(0); r = Fr(1)
        for k in range(n+1):
            if k > z: break
            tot += c[k]*r
            r *= Fr(z-k, M-k)
        vals.append(tot)
    return vals
if __name__ == '__main__':
    M = int(sys.argv[1]); s = int(sys.argv[2])
    for n in map(int, sys.argv[3:]):
        v = Qvals(n, M, s)
        mx = max(v); mn = min(v)
        print(M, s, n, 'max', float(mx), 'at', v.index(mx), 'min', float(mn), 'at', v.index(mn), flush=True)
