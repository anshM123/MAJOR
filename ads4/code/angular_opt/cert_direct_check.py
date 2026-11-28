"""Independent check of a Lemma-M certificate: build the multipliers lambda_k = -W_B c_k / D_k^2,
c_k = C(M,k) Delta^k q(0), and evaluate the dual inequalities DIRECTLY (no Lemma M) on many heavy states, mpmath."""
import json, sys, random
from fractions import Fraction as Fr
import mpmath as mp
from cert_lemmaM import gram_rows, WB, theta
mp.mp.dps = 120
fn = sys.argv[1]; g = int(sys.argv[2])
d = json.load(open(fn)); Mb = d['Mb']; n = d['n']
a = [Fr(x) for x in d['a']]; nus = [int(x) for x in d['nus']]; shift = Fr(d['shift'])
# q values at z = 0..n+1 exactly
q = [Fr(0)]*(n+2)
for j, t in gram_rows(Mb, n):
    for z in range(n+2): q[z] += a[j]/nus[j]*t[z]
q = [v+shift for v in q]
M = g-2
# forward differences at 0
diffs = []; row = q[:]
for k in range(n+1):
    diffs.append(row[0]); row = [row[i+1]-row[i] for i in range(len(row)-1)]
from math import comb
c = [comb(M, k)*diffs[k] for k in range(n+1)]
print('c_0..c_3 =', [float(x) for x in c[:4]], ' c_2 = C(M,2)*S:', float(c[2]), float(comb(M,2)*(q[2]-2*q[1]+1)))
cm = [mp.mpf(x.numerator)/x.denominator for x in c]
W = lambda b, D: ((D-1)*(D-2) - b*(b-1))/mp.mpf(3) + mp.mpf(3)/2
wb = mp.mpf(WB(g).numerator)/WB(g).denominator
def F(b, J, D):
    C = J*(J+1); T = (D-1)*(D-2); m = mp.mpf(1); s = cm[0]
    for k in range(n):
        if k+1 > J-b: break
        y = (b+k)*(b+k+1); m *= (C-y)/(T-y); s += cm[k+1]*m
    return s
random.seed(3); worst = (mp.inf, None)
Js = list(range(0, 60)) + [random.randint(60, 5000) for _ in range(150)] + [10**4, 10**5, 10**6, 10**8]
for J in Js:
    for b in (0, 4):
        if b == 0 and J % 2: continue
        if b == 4 and J < 4: continue
        for t in [0, 0.3, 1, 3, 10, 30, 100, 1000, 10**4, 10**6]:
            D = mp.mpf(J) + g + t
            f = F(b, J, D)
            marg = (W(b, D) - wb*f) if b == 4 else (W(b, D) + wb*f)
            marg /= W(b, D)
            if marg < worst[0]: worst = (marg, (b, J, float(t)))
print('min normalised margin over sampled heavy states:', mp.nstr(worst[0], 8), worst[1])
