"""Size of the certificate functional on LIGHT states (twist < g): controls the light-state term.
Reports max over sampled light states of |Lambda^(b)|/W^(b) and |Lambda|."""
import json, sys
from fractions import Fraction as Fr
import mpmath as mp
from math import comb
from cert_lemmaM import gram_rows, WB
mp.mp.dps = 200
fn = sys.argv[1]; g = int(sys.argv[2])
d = json.load(open(fn)); Mb = d['Mb']; n = d['n']
a = [Fr(x) for x in d['a']]; nus = [int(x) for x in d['nus']]; shift = Fr(d['shift'])
q = [Fr(0)]*(n+2)
for j, t in gram_rows(Mb, n):
    for z in range(n+2): q[z] += a[j]/nus[j]*t[z]
q = [v+shift for v in q]; M = g-2
diffs = []; row = q[:]
for k in range(n+1):
    diffs.append(row[0]); row = [row[i+1]-row[i] for i in range(len(row)-1)]
cm = [mp.mpf(comb(M, k)*diffs[k].numerator)/diffs[k].denominator for k in range(n+1)]
print('sum |c_k| =', mp.nstr(sum(abs(x) for x in cm), 5), ' max|c_k| =', mp.nstr(max(abs(x) for x in cm), 5))
wb = mp.mpf(WB(g).numerator)/WB(g).denominator
def F(b, J, D):
    C = J*(J+1); T = (D-1)*(D-2); m = mp.mpf(1); s = cm[0]
    for k in range(n):
        if k+1 > J-b: break
        y = (b+k)*(b+k+1); m *= (C-y)/(T-y); s += cm[k+1]*m
    return s
W = lambda b, D: ((D-1)*(D-2) - b*(b-1))/mp.mpf(3) + mp.mpf(3)/2
best = (0, None)
for h in [2.01, 2.5, 3, 4, 6, 6.5, 8, 10, 15, 20, 50, 100, 200, 400, 600, 800, 900, 950, 990, 999]:
    for J in [0, 2, 4, 5, 6, 7, 8, 10, 12, 16, 20, 30, 50, 100, 300, 1000, 3000, 10**4, 10**5]:
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            D = mp.mpf(J) + h
            v = abs(wb*F(b, J, D))/W(b, D)
            if v > best[0]: best = (v, (b, J, h))
print('max |Lambda|/W on sampled light states:', mp.nstr(best[0], 6), best[1])
print('per-twist max |Lambda|/W (over J):')
for h in [2.01, 3, 5.9, 6, 7, 8, 10, 14, 20, 30, 50, 80, 120, 200, 300, 500, 700, 900, 990]:
    best = 0
    for J in list(range(0, 41)) + [50, 70, 100, 150, 200, 300, 500, 1000, 2000, 5000, 10**4, 10**5, 10**6]:
        for b in (0, 4):
            if b == 0 and J % 2: continue
            if b == 4 and J < 4: continue
            D = mp.mpf(J) + h
            best = max(best, abs(wb*F(b, J, D))/W(b, D))
    print('  twist', h, ' max|Lambda|/W =', mp.nstr(best, 4))
