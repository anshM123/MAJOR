"""
tailaudit/ta_lemmas_check.py [DV, mpmath]: spot checks of the remaining TAIL_CERT inequalities at adversarial tail states.
 (R)  recurrence constants (§2.3): |T/D - 1|, |beta(s-2beta-1)|/(w|D|) <= 1e-6 and w^2|D| >= 3.9 a^2, D = N(N+s) - beta(beta-s),
      for |beta| up to 6200 (tc_assemble applies the recurrence up to t = 6000 + j), L in {18, 20}, a in {998, 998.5, 1e4, 1e6}.
 (D)  dressing constants: T^{-1}/w <= w/(3.9a^2), T^{-1} <= w^2/(3.9a^2), and (k = 3) T^{-1}/w <= 0.2515 w/a^2.
 (Q)  Theorem P window: |lam - 1| <= 7 sqrt(w) on |V - N/M| <= eta, eta^2 = (0.75a + pi tau/2 + 60)/M, tau <= 1146.
 (LL) Lemma Loc, the two Cauchy integrals: I1 = rho^{-n} 2F1(g/2,g/2;1;rho^2)/C(N,a) (g = a+1-m), I2 = rho^{-n} 2F1((a+1)/2,(a+1)/2;1;
      rho^2) E(1-B)^m/C(N,a), rho = n/(n+a-m), against the stated 2.53 (a/N)^m and 2.53 ((a+m)/(a-m))^m (a/N)^m, and the stated
      (1.0006 w) >= 2a/N.
 (PB) FULLPOS3 Lemma B as used in tc_remainder: |P_{-1/2+i nu}(X)| <= (pi/2) e^{nu |Im arccosh X|} for X = (1+xi)/(1-xi), |xi| < 1;
      and the |nu| <= 0.4 Pfaff bound |P_{-1/2+i nu}(X)| <= 2^{0.9} e^{0.2 pi}/(1-|xi|) (nu complex, |nu| = 0.4).
Output: cache/ta_lemmas_check.txt
"""
import os, math, random
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s):
    print(s, flush=True); out.append(s)

mp.mp.dps = 60

def tail_state(a, b, L):
    s = 2*b + 1
    n = mp.mpf(2)*a*mp.e**L
    for _ in range(80):
        N = n + a; T = N*(N + s) + b*(b + 1); w = a*(N + n + s)/T
        n = n*mp.e**(L + mp.log(w))
    n = mp.floor(n)
    N = n + a; T = N*(N + s) + b*(b + 1); w = a*(N + n + s)/T
    return n, N, T, w, s

# (R), (D), (Q)
worst = dict(TD=0, BW=0, W2D=1e99, d1=0, d2=0, d3=0, lam=0)
for L in (18, 20, 30):
    for a in (mp.mpf(998), mp.mpf('998.5'), mp.mpf(10**4), mp.mpf(10**6), mp.mpf(10)**12):
        for b in (0, 4):
            n, N, T, w, s = tail_state(a, b, L)
            for beta in (mp.mpc(0.25, 6200), mp.mpc(-199.75, 6000), mp.mpc(3.5 - 200, 1146), mp.mpc(2.5, 0), mp.mpc(-197.5, 0)):
                D = N*(N + s) - beta*(beta - s)
                worst['TD'] = max(worst['TD'], abs(T/D - 1))
                worst['BW'] = max(worst['BW'], abs(beta*(s - 2*beta - 1))/(w*abs(D)))
                worst['W2D'] = min(worst['W2D'], w*w*abs(D)/(3.9*a*a))
            worst['d1'] = max(worst['d1'], (1/T/w)/(w/(3.9*a*a)))
            worst['d2'] = max(worst['d2'], (1/T)/(w*w/(3.9*a*a)))
            worst['d3'] = max(worst['d3'], (1/T/w)/(0.2515*w/(a*a)))
            M = N + n + s
            for tau in (0, 1146):
                eta = mp.sqrt((0.75*a + mp.pi*tau/2 + 60)/M)
                for V in (N/M - eta, N/M + eta):
                    worst['lam'] = max(worst['lam'], abs(1/V - 2)/(7*mp.sqrt(w)))
pr(f'(R) max |T/D - 1| = {mp.nstr(worst["TD"], 3)};  max |beta(s-2beta-1)|/(w|D|) = {mp.nstr(worst["BW"], 3)}  (need <= 1e-6);  '
   f'min w^2|D|/(3.9a^2) = {mp.nstr(worst["W2D"], 6)} (need >= 1)')
pr(f'(D) max (T^-1/w)/(w/(3.9a^2)) = {mp.nstr(worst["d1"], 6)}, max T^-1/(w^2/(3.9a^2)) = {mp.nstr(worst["d2"], 6)}, '
   f'max (T^-1/w)/(0.2515 w/a^2) = {mp.nstr(worst["d3"], 6)}  (need <= 1)')
pr(f'(Q) max |lam-1|/(7 sqrt w) over the window edges = {mp.nstr(worst["lam"], 6)}  (need <= 1)')

# (LL)
pr('(LL) Lemma Loc Cauchy integrals (ratio = computed / stated bound; need <= 1)')
mp.mp.dps = 40
for L in (8, 20):
    for a in (mp.mpf(998), mp.mpf('1234.5'), mp.mpf(10**4)):
        for m in (2, 50, 196):
            n, N, T, w, s = tail_state(a, 0, L)
            rho = n/(n + a - m)
            lC = mp.loggamma(N + 1) - mp.loggamma(a + 1) - mp.loggamma(n + 1)
            g = a + 1 - m
            I1 = mp.exp(-n*mp.log(rho) - lC)*mp.hyp2f1(g/2, g/2, 1, rho**2)
            EB = mp.rf(a, m)/mp.rf(N + s, m)
            I2 = mp.exp(-n*mp.log(rho) - lC)*mp.hyp2f1((a + 1)/2, (a + 1)/2, 1, rho**2)*EB
            b1 = 2.53*(a/N)**m; b2 = 2.53*((a + m)/(a - m))**m*(a/N)**m
            pr(f'   L={L} a={mp.nstr(a, 6)} m={m}: I1/bound = {mp.nstr(I1/b1, 5)}, I2/bound = {mp.nstr(I2/b2, 5)}, '
               f'2a/N / (1.0006 w) = {mp.nstr(2*a/N/(1.0006*w), 6)}')

# (PB)
mp.mp.dps = 30
random.seed(3)
wB = 0; wP = 0
for _ in range(400):
    r = random.choice((0.3, 0.7, 0.95, 0.999)); th = random.uniform(-math.pi, math.pi)
    xi = mp.mpc(r*math.cos(th), r*math.sin(th))
    X = (1 + xi)/(1 - xi)
    zeta = mp.acosh(X)
    for nu in (0.5, 3.0, 17.0, 60.0):
        P = mp.legenp(-0.5 + 1j*nu, 0, X, type=3)
        wB = max(wB, abs(P)/((mp.pi/2)*mp.e**(nu*abs(zeta.imag))))
    for k in range(8):
        nu = 0.4*mp.e**(1j*(2*math.pi*k/8 + 0.1))
        P = (1 - xi)**(0.5 - 1j*nu)*mp.hyp2f1(0.5 - 1j*nu, 0.5 - 1j*nu, 1, xi)
        wP = max(wP, abs(P)/(2**0.9*mp.e**(0.2*mp.pi)/(1 - abs(xi))))
pr(f'(PB) max |P|/((pi/2) e^(nu|Im zeta|)) over 1600 samples = {mp.nstr(wB, 5)} (need <= 1);  '
   f'|nu| = 0.4 Pfaff bound: max ratio = {mp.nstr(wP, 5)} (need <= 1)')
open(os.path.join(HERE, 'cache', 'ta_lemmas_check.txt'), 'w').write('\n'.join(out) + '\n')
