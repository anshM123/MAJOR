"""Exact-rational sanity check of Lemma PB+ (positive basis for d <= a+1, real twist a): L[z^k(1-z)^d] >= 0.
Also scans the smallest d at which L[(1-z)^d] < 0 for a few states (shows the range is state dependent)."""
from fractions import Fraction as Fr
from math import comb
import random
def U(J, b, a, kmax):
    T = (a+J+1)*(a+J); C = J*(J+1); out = [Fr(1)]
    for j in range(kmax):
        uj = (b+j)*(b+j+1); out.append(out[-1]*Fr(C-uj)/(T-uj))
    return out
def Lzk1zd(Us, k, d):
    return sum((-1)**i*comb(d, i)*Us[k+i] for i in range(d+1) if k+i < len(Us))
random.seed(1)
bad = 0
for trial in range(60):
    b = random.choice((0, 4)); J = random.randint(b, 40); a = Fr(random.randint(2000, 4000), 4)  # real twist 500..1000
    n = J-b; Us = U(J, b, a, n+1) + [Fr(0)]*0
    Us = Us[:n+1]
    for d in (0, 1, 5, int(a)//2, int(a), int(a)+1):
        for k in (0, 1, n//2, n):
            v = Lzk1zd(Us, k, d)
            if v < 0: bad += 1; print('VIOLATION', J, b, a, k, d, float(v))
print('Lemma PB+ exact checks: violations =', bad)
for (J, b, a) in [(5, 4, 998), (20, 0, 998), (40, 0, 998), (40, 4, 998)]:
    Us = U(J, b, a, J-b)
    d = 0
    while Lzk1zd(Us, 0, d) >= 0 and d < 10**6: d = d+1 if d < 2000 else int(d*1.05)
    print('state J=%d b=%d a=%d: first d with L[(1-z)^d]<0 (coarse search) ~ %d  (a+1=%d)' % (J, b, a, d, a+1))
