"""Target 1a (d'): pointwise check of the magic identity with direct Legendre functions of complex degree, per b, with the
truncation order varied (N = 30, 45, 60) to separate truncation from the identity.  Output: cache/t1d_pointwise.txt"""
import mpmath as mp, sympy as sp
from blocks2 import recursion
mp.mp.dps = 30
out = open('cache/t1d_pointwise.txt', 'w')
def log(s): print(s); out.write(s + '\n'); out.flush()
conv = lambda e: mp.mpc(sp.re(e).evalf(60), sp.im(e).evalf(60)) if not e.is_real else mp.mpf(e.evalf(60))
def G(B, l0, r, eta, cache):
    s = 0
    for N, d in enumerate(B):
        for m, c in d.items():
            if (l0, m) not in cache: cache[(l0, m)] = mp.legenp(l0 + m, 0, eta, type=3)
            s += c * r**N * cache[(l0, m)]
    return s
pts = [(mp.mpf('0.2'), mp.mpf('0.6')), (mp.mpf('0.3'), mp.mpf('0.4')), (mp.mpf('0.1'), mp.mpf('1.5')), (mp.mpf('0.05'), mp.mpf('1.2'))]
for nu in (mp.mpf('0.7'), mp.mpc('2.3', '0.2')):
    lv = mp.mpf(-0.5) + 1j * nu
    for b in range(-4, 5):
        row = []
        Bs = (recursion(lv, b, 1, 60, conv), recursion(lv - 1, b, 0, 60, conv), recursion(lv + 1, b, 0, 60, conv))
        for (r, beta) in pts:
            eta = mp.cosh(beta); cache = {}
            res = []
            for NM in (30, 45, 60):
                B1, Bm, Bp = (x[:NM + 1] for x in Bs)
                L = 2 * eta * (1 + r**2) * G(B1, lv, r, eta, cache)
                R = (1 + 1j / (2 * nu)) * G(Bm, lv - 1, r, eta, cache) + (1 - 1j / (2 * nu)) * G(Bp, lv + 1, r, eta, cache)
                res.append(abs(L - R) / abs(L))
            row.append('/'.join(mp.nstr(x, 2) for x in res))
        log(f'nu={mp.nstr(nu,3)} b={b:+d}: rel. residual at N=30/45/60 for (r,beta) = (0.2,0.6) (0.3,0.4) (0.1,1.5) (0.05,1.2): ' + '  '.join(row))
