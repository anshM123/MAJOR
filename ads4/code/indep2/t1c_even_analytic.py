"""Target 1 (side claims of PW2): the Dp = 0 (0, b) blocks Ghat^0_b(mu) are even in mu, and analytic at the would-be
poles mu = +-i/2, +-3i/2 (zeros of the Casimir denominators), b = -4..4, coefficients to N = 14.
Evenness: B[N][m](l) = B[N][-m](-1-l)  (P_{-1-l+m} = P_{l-m}).  Analyticity: B[N][m] at l = l0 + eps stays bounded as eps -> 0.
Also: f = v - u = eta(1+r^2)/(2r) in the H7 frame (exact, cross-ratios u = x12^2 x34^2/(x13^2 x24^2), v = x14^2 x23^2/(x13^2 x24^2)).
Output: cache/t1c_even_analytic.txt
"""
import mpmath as mp
import sympy as sp
from blocks2 import recursion

mp.mp.dps = 60
out = open('cache/t1c_even_analytic.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
conv = lambda e: mp.mpf(e.evalf(70)) if e.is_real else mp.mpc(sp.re(e).evalf(70), sp.im(e).evalf(70))

NM = 14
# evenness
worst = 0
for nu in (mp.mpf('0.7'), mp.mpc('2.3', '0.4')):
    l = mp.mpf(-0.5) + 1j * nu
    for b in range(-4, 5):
        B1 = recursion(l, b, 0, NM, conv); B2 = recursion(-1 - l, b, 0, NM, conv)
        for N in range(NM + 1):
            if not B1[N]:
                assert not B2[N]; continue
            sc = max(abs(c) for c in B1[N].values())
            for m, c in B1[N].items():
                worst = max(worst, abs(c - B2[N].get(-m, 0)) / sc)
log(f'evenness of Ghat^0_b in mu, b = -4..4, N <= {NM}, two nu: max |B(l)[N,m] - B(-1-l)[N,-m]|/max|B[N]| = {mp.nstr(worst, 3)}')

# analyticity at mu = i/2 (l = 0, -1) and mu = 3i/2 (l = 1, -2)
for l0 in (0, -1, 1, -2):
    maxes = []; worst = 0
    for e in (mp.mpf('1e-12'), mp.mpf('1e-24'), mp.mpf('1e-36')):
        mx = 0
        for b in range(-4, 5):
            B = recursion(l0 + e * (1 + 1j), b, 0, NM, conv)
            mx = max(mx, max(abs(c) for d in B for c in d.values()))
            if e == mp.mpf('1e-36'):
                Bref = recursion(l0 + mp.mpf('1e-24') * (1 + 1j), b, 0, NM, conv)
                for N in range(NM + 1):
                    if B[N]:
                        sc = max(abs(c) for c in B[N].values())
                        worst = max(worst, max(abs(c - Bref[N].get(m, 0)) for m, c in B[N].items()) / sc)
        maxes.append(mx)
    log(f'l0 = {l0:+d} (mu = {"+-i/2" if l0 in (0, -1) else "+-3i/2"}), b = -4..4, N <= {NM}: max|B| at eps = 1e-12/1e-24/1e-36: '
        f'{", ".join(mp.nstr(x, 6) for x in maxes)}; max_N |B(1e-36) - B(1e-24)|/max|B[N]| = {mp.nstr(worst, 3)}  (a pole would scale like 1/eps)')

# frame identity
r, eta = sp.symbols('r eta', positive=True)
# timelike unit vectors yhat, xhat with -yhat.xhat = eta; x1 = -r yhat, x3 = r yhat, x2 = -xhat, x4 = xhat (metric -,+,+)
def sq(a, bb, c):     # (a xhat + bb yhat)^2 with xhat^2 = yhat^2 = -1, xhat.yhat = -eta
    return -a**2 - bb**2 - 2 * a * bb * eta
x12 = sq(1, -r, 0); x34 = sq(-1, r, 0); x13 = sq(0, -2 * r, 0); x24 = sq(2, 0, 0); x14 = sq(-1, -r, 0); x23 = sq(-1, -r, 0)
u = x12 * x34 / (x13 * x24); v = x14 * x23 / (x13 * x24)
log(f'f = v - u - eta(1+r^2)/(2r) simplifies to: {sp.simplify(v - u - eta * (1 + r**2) / (2 * r))}')
