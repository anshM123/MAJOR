"""rv_state.py (review): (1) adversarial scan of the state relations 0.49d <= a/N <= 0.51d, D d^2 >= 3.9 a^2 (and T = D + b(b+1)) for
L >= 4.5 in exact rationals; (2) K0 = rho^-n I_W / C(N,a) and the theta-fraction bound; (3) P(B < 1-b) vs (1.02 d/b)^998."""
from fractions import Fraction as Fr
import math, random
import mpmath as mp

def rel(n, a, b):
    s = 2*b + 1; N = n + a; C0 = n*(n + s); D = N*(N + s)
    d = Fr(D - C0, D)
    return d, Fr(a, N)/d, d*d*D/(a*a), d*d*(D + b*(b + 1))/(a*a)

def n_for_L(a, b, L):
    # smallest n with -log d >= L  (d decreasing in n)
    lo, hi = 1, 10**15
    while lo < hi:
        m = (lo + hi)//2
        d = rel(m, a, b)[0]
        if -math.log(float(d)) >= L: hi = m
        else: lo = m + 1
    return lo

def part1():
    worst = dict(rmax=(0, None), rmin=(9, None), Tmin=(99, None))
    random.seed(1)
    cases = []
    for b in (0, 4):
        for a in [998, 999, 1000, 1500, 3000, 3600, 10**4, 10**5, 10**6, 10**8, 10**10] + [random.randint(998, 10**7) for _ in range(60)]:
            for L in (4.5, 4.5000001, 4.6, 5, 6, 8, 16, 20):
                n = max(n_for_L(a, b, L), 1001 - b)
                cases.append((n, a, b))
            # also the smallest admissible n (J >= 1001) regardless of L, to see where the relations fail
    for (n, a, b) in cases:
        d, r, Dd, Td = rel(n, a, b)
        if -math.log(float(d)) < 4.5: continue
        if r > worst['rmax'][0]: worst['rmax'] = (r, (n, a, b))
        if r < worst['rmin'][0]: worst['rmin'] = (r, (n, a, b))
        if Dd < worst['Tmin'][0]: worst['Tmin'] = (Dd, (n, a, b))
    for k, v in worst.items():
        print(k, float(v[0]), v[1])
    # where would the relations fail? scan L downward at a = 998, b = 4 and a = 1e8
    for (a, b) in ((998, 4), (998, 0), (10**8, 4)):
        for L in (3.0, 2.6, 2.5, 2.0):
            n = n_for_L(a, b, L)
            d, r, Dd, Td = rel(n, a, b)
            print(f'  a={a} b={b} L={L}: a/(N d) = {float(r):.5f}, D d^2/a^2 = {float(Dd):.4f}')

def part2():
    mp.mp.dps = 40
    for L in (4.5, 5.0, 5.5, 8.0, 16.0):
        for a in (998, 1500, 3000, 10**5):
            for b in (0, 4):
                n = n_for_L(a, b, L); N = n + a
                rho = mp.mpf(n)/N
                g = mp.mpf(a + 1)/2
                IW = mp.hyp2f1(g, g, 1, rho**2)
                lC = mp.loggamma(N + 1) - mp.loggamma(a + 1) - mp.loggamma(n + 1)
                K0 = mp.exp(-n*mp.log(rho) + mp.log(IW) - lC)
                # lower bound used for I_W:  e^{-0.51} (1-rho)^{-a} / (pi sqrt a)
                lb = mp.exp(-0.51)*(1 - rho)**(-a)/(mp.pi*mp.sqrt(a))
                # exact theta-fraction at t = 0.02 vs the bound formula (with exact 1-rho)
                t = mp.mpf('0.02')
                W = lambda th: ((1 - rho)**2 + 4*rho*mp.sin(th/2)**2)**(-(a + 1)/mp.mpf(2))
                tail = mp.quad(W, [t, 2*t, 0.1, 0.5, mp.pi])/mp.pi    # (1/2pi) * 2 * int_t^pi
                fr = tail/IW
                q = mp.sqrt((1 - rho)**2 + 4*rho*mp.sin(t/2)**2)
                fb = mp.pi*mp.sqrt(a)*mp.e**0.51/(1 - rho)*((1 - rho)/q)**(a + 1)
                print(f'L={L} a={a} b={b}: K0 = {mp.nstr(K0, 8)}  IW/lowerbound = {mp.nstr(IW/lb, 5)}  log frac(0.02) exact {mp.nstr(mp.log(fr), 7)} bound {mp.nstr(mp.log(fb), 7)}', flush=True)

def part3():
    mp.mp.dps = 30
    for L in (4.5, 5.0, 5.5):
        for a in (998, 1500, 3000, 20000):
            b = 4; s = 9
            n = n_for_L(a, b, L); N = n + a
            d = float(rel(n, a, b)[0])
            for bb in (0.04, 0.0640, 0.15, 0.95):
                # B ~ Beta(n+s, a): P(B < 1 - bb)
                # 1-B ~ Beta(a, n+s) (integer params): P(1-B > bb) = P(Bin(a+n+s-1, bb) <= a-1)
                M = a + n + s - 1; lb = mp.log(bb); l1 = mp.log(1 - bb)
                terms = [mp.loggamma(M + 1) - mp.loggamma(k + 1) - mp.loggamma(M - k + 1) + k*lb + (M - k)*l1 for k in range(a)]
                mx = max(terms); P = mp.exp(mx)*mp.fsum(mp.exp(t - mx) for t in terms)
                bound = 998*mp.log(1.02*d/bb)
                print(f'L={L} a={a} b={bb}: log P(B<1-b) = {mp.nstr(mp.log(P), 8) if P > 0 else "-inf"}  bound {mp.nstr(bound, 8)}', flush=True)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1: part3()
    else: part1(); part2(); part3()
