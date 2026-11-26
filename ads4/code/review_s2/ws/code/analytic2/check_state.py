"""check_state.py: [exact + DV] state relations at L >= Lmin (exact rationals), K0 and the Beta tier bound at L in {3.75, 4, 4.25}
(same formulas as code/review_l5/rv_state.py, run at the lower L)."""
import sys, math, random
sys.path.insert(0, '/home/user/MAJOR/ads4/code/review_l5')
import rv_state as RV
import mpmath as mp
Lmin = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
random.seed(2)
worst = dict(rmax=(0, None), rmin=(9, None), Dmin=(99, None), Nmin=(1e99, None))
for b in (0, 4):
    for a in [998, 999, 1000, 1500, 3000, 10**4, 10**6, 10**9] + [random.randint(998, 10**7) for _ in range(40)]:
        for L in (Lmin, Lmin + 1e-7, Lmin + 0.1, Lmin + 0.25, Lmin + 0.5, 5.0):
            n = max(RV.n_for_L(a, b, L), 1001 - b)
            d, r, Dd, Td = RV.rel(n, a, b)
            if -math.log(float(d)) < Lmin: continue
            if r > worst['rmax'][0]: worst['rmax'] = (r, (n, a, b))
            if r < worst['rmin'][0]: worst['rmin'] = (r, (n, a, b))
            if Dd < worst['Dmin'][0]: worst['Dmin'] = (Dd, (n, a, b))
            if (n + a)/a < worst['Nmin'][0]: worst['Nmin'] = ((n + a)/a, (n, a, b))
for k, v in worst.items():
    print(k, float(v[0]), v[1])
mp.mp.dps = 40
for L in (Lmin, Lmin + 0.25, Lmin + 0.5):
    for a in (998, 1500, 3000, 10**5):
        for b in (0, 4):
            n = RV.n_for_L(a, b, L); N = n + a
            rho = mp.mpf(n)/N; g = mp.mpf(a + 1)/2
            IW = mp.hyp2f1(g, g, 1, rho**2)
            lC = mp.loggamma(N + 1) - mp.loggamma(a + 1) - mp.loggamma(n + 1)
            K0 = mp.exp(-n*mp.log(rho) + mp.log(IW) - lC)
            lb = mp.exp(-0.51)*(1 - rho)**(-a)/(mp.pi*mp.sqrt(a))
            print(f'L={L} a={a} b={b}: K0 = {mp.nstr(K0, 8)}  IW/lowerbound = {mp.nstr(IW/lb, 5)}', flush=True)
