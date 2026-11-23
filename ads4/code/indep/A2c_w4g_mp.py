"""k = 4 regular weight w4g beyond the double-precision floor: S_b - S_LR in 45-digit arithmetic, rho = 0.5, b = 0, 2,
beta = 9..13 (the double-precision runs of A2 and of EESTEP/LRC hit the floor at beta ~ 9.5-10)."""
import os, sys, math, time
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A1b_negb_mp import recursion_mp
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A2c_w4g_mp.txt')
mp.mp.dps = 45
t0 = time.time()
NG, L, NMAX, rho = 160, 24, 110, mp.mpf('0.5')
xs, ws = zip(*[(x, w) for x, w in zip(*__import__('numpy').polynomial.legendre.leggauss(NG))])
# refine GL nodes
nodes, wts = [], []
for x0 in xs:
    x = mp.mpf(x0)
    for _ in range(5):
        p = mp.legendre(NG, x); dp = mp.diff(lambda u: mp.legendre(NG, u), x); x -= p / dp
    dp = mp.diff(lambda u: mp.legendre(NG, u), x)
    nodes.append(L / 2 * (x + 1)); wts.append(L / 2 * 2 / ((1 - x**2) * dp**2))
w4g = lambda v: v**2 * (v**2 + 1) * (v**2 + 4) * mp.exp(-v**2 / 9)
lines = []
for b in (0, 2):
    Bs = [recursion_mp(v, b, NMAX) for v in nodes]
    res = []
    for beta in (9, 10, 11, 12, 13):
        beta = mp.mpf(beta); r = rho * mp.exp(-beta); eta = mp.cosh(beta)
        S = 0; SLR = 0
        for v, w, B in zip(nodes, wts, Bs):
            l = mp.mpc(-0.5, v)
            Q = {0: mp.legenp(l, 0, eta, type=3), 1: r * mp.legenp(l + 1, 0, eta, type=3)}
            for m in range(1, NMAX + 1):
                mu = l + m; Q[m + 1] = ((2 * mu + 1) * (r * eta) * Q[m] - mu * r * r * Q[m - 1]) / (mu + 1)
            for m in range(1, NMAX + 2): Q[-m] = mp.conj(Q[m])
            G = 0
            for N in range(len(B)):
                for m, c in B[N].items(): G += c * Q[m] * r**(N - abs(m))
            S += w * w4g(v) * G; SLR += w * w4g(v) * Q[0]
        d = abs(eta**3 * (S - SLR))
        res.append((float(beta), d))
        print(b, float(beta), mp.nstr(d, 8), f'{time.time()-t0:.0f}s', flush=True)
    rates = [math.log(float(res[k + 1][1] / res[k][1])) for k in range(len(res) - 1)]
    s = f'b={b}, rho=0.5: |S_b - S_LR| at beta=9..13: ' + ' '.join(mp.nstr(x[1], 4) for x in res) + '; local exponents ' + ' '.join(f'{x:+.3f}' for x in rates)
    print(s); lines.append(s)
open(OUT, 'w').write('\n'.join(lines) + '\n')
