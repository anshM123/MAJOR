#!/usr/bin/env python3
"""Referee check: Lemma V difference recurrence and nb_delta seed on adversarial physical states.
X(beta)=d^beta S(beta), S(beta)=sum_k (beta)_k/k! U_k, U_k=prod_{j<k}(c-h j(j+s))/(1-h j(j+s)).
"""
import mpmath as mp, sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/MAJOR/ads4/code/sliver_ext')
from contour_direct import nb_delta
from base_bounds import hcap
mp.mp.dps = 80

def state(n, a, b):
    s = 2*b+1; N = n+a; D = mp.mpf(N)*(N+s); h = 1/D; c = mp.mpf(n)*(n+s)/D
    return dict(n=n, a=a, b=b, s=s, N=N, h=h, c=c, d=1-c)

def U_list(st):
    n, s, h, c = st['n'], st['s'], st['h'], st['c']
    U = [mp.mpf(1)]
    for j in range(n):
        U.append(U[-1]*(c-h*j*(j+s))/(1-h*j*(j+s)))
    return U

def X(st, U, beta):
    tot = mp.mpc(0); t = mp.mpc(1)
    for k in range(len(U)):
        tot += t*U[k]
        t *= (beta+k)/(k+1)
    return mp.power(st['d'], beta)*tot

def lemmaV(st, U, beta):
    s, h, d, c = st['s'], st['h'], st['d'], st['c']
    g = h/d**2
    D = lambda b_: X(st, U, b_)-1
    lhs = (1-h*beta*(beta-s))*D(beta-1)
    rhs = (1+g*d*beta*(s-1-2*beta))*D(beta)+g*beta*(beta+1)*D(beta+1)+g*beta*c*(beta*c+1+s*d)
    return abs(lhs-rhs), abs(lhs)

def solve_n(a, b, L):
    # smallest integer n with d<=e^{-L}  (d decreasing in n at fixed a)
    s = 2*b+1; target = mp.e**(-L)
    lo, hi = 1, 10**7
    while hi-lo > 1:
        m = (lo+hi)//2; N = m+a
        d = mp.mpf(a)*(N+m+s)/(mp.mpf(N)*(N+s))
        if d <= target: hi = m
        else: lo = m
    return hi

if __name__ == '__main__':
    # 1. Lemma V on small and physical states
    for (n, a, b) in [(30, 40, 0), (25, 13, 4), (997, 998, 4)]:
        st = state(n, a, b); U = U_list(st)
        for beta in [mp.mpc('0.3', '0.25'), mp.mpc('-7.75', '0.25'), mp.mpc('2.5', '900')]:
            e, m = lemmaV(st, U, beta)
            print('lemmaV', n, a, b, mp.nstr(beta, 6), 'defect', mp.nstr(e, 3), 'scale', mp.nstr(m, 3))
    # 2. nb_delta vs exact |X-1| at Re alpha in [1/4, 3.5], Im alpha <= .25
    worst = 0
    for L in [2.6, 2.9, 3.65]:
        for a in [998, 1500, 4000]:
            for b in [0, 4]:
                n = solve_n(a, b, L); st = state(n, a, b); U = U_list(st)
                s = st['s']; du = F(str(mp.nstr(st['d'], 20)))*(1+F(1, 10**12))
                g = hcap(du, F(998))/du**2
                for sig in ['0.25', '0.75', '2.5', '3.5']:
                    for t in ['0', '0.25']:
                        beta = mp.mpc(sig, t)
                        ex = abs(X(st, U, beta)-1)
                        bd = F(nb_delta(g, du, s, F(sig), F(sig), F(t)), 10**40)
                        r = ex/mp.mpf(bd.numerator)/mp.mpf(bd.denominator)**-1 if False else ex/(mp.mpf(bd.numerator)/bd.denominator)
                        worst = max(worst, r)
                        if r > 0.5 or sig == '0.25' and t == '0.25' and a == 998:
                            print('nb', L, a, b, n, sig, t, 'exact', mp.nstr(ex, 4), 'bound', float(bd), 'ratio', mp.nstr(r, 4))
    print('nb worst exact/bound', mp.nstr(worst, 5), 'PASS' if worst < 1 else 'FAIL')
