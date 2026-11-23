"""Target 2: Lemma LR-infinity (LARGESPIN_L §3).
Omega_k(beta) = int_0^1 psi_k(p) P_{-1/2+i Lam p}(cosh beta) dp   (NF02, prl/data/nf02_psi.json; Lam = sqrt(998*999)).
Independent evaluation: Mehler-Dirichlet in t with the exact Fourier transform Phi(x) = int_0^1 psi cos(px) dp (Taylor series
for x < 25, exact integration-by-parts form for x >= 25), mpmath at 40 digits, Gauss-Legendre panels of length <= pi/(2 Lam)
while Lam t <= 3000 (beyond, the oscillatory p = 1 part is < 1e-25 and only the smooth part is integrated, on coarse panels),
and t = beta - u^2 on [beta/2, beta].
Checks: (i) two-term form a1 M1/Lam^2 - a2/(2 bhat^3); (ii) the bound F_LR >= (1 - r)A1 M1/Lam^2 + 2047.1/(2 bhat^3) with the stated
r table; (iii) the Hadamard identity of step 3; (iv) the J*(m) table.   Output: cache/t2_lrinf.txt
"""
import json, time, sys
import mpmath as mp
from fractions import Fraction

mp.mp.dps = 40
DATA = '../../prl/data/'
out = open('cache/t2_lrinf.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

js = json.load(open(DATA + 'nf02_psi.json'))
k0 = mp.mpf(js['k0'])
psi = {2: [mp.mpf(Fraction(c).numerator) / Fraction(c).denominator / k0 for c in js['psi2_exact']],
       3: [mp.mpf(Fraction(c).numerator) / Fraction(c).denominator / k0 for c in js['psi3_exact']]}
Lam = mp.sqrt(mp.mpf(998) * 999)

def deriv_at(c, j, p):
    # j-th derivative of sum c_n p^n at p
    s = 0
    for n in range(j, len(c)):
        s += c[n] * mp.ff(n, j) * p**(n - j)
    return s

def phi_factory(c):
    n = len(c)
    d0 = [deriv_at(c, j, 0) for j in range(n)]
    d1 = [deriv_at(c, j, 1) for j in range(n)]
    def phi_smooth(x):          # -Re sum_j (-1)^j psi^(j)(0)/(ix)^(j+1)
        s = 0
        for j in range(1, n, 2):
            s += (-1)**((j + 1) // 2) * d0[j] / x**(j + 1)
        return s
    def phi_osc(x):
        s = 0
        e = mp.expj(x)
        for j in range(n):
            if d1[j] != 0:
                s += (-1)**j * d1[j] * e / (1j * x)**(j + 1)
        return mp.re(s)
    def phi(x):
        if x < 25:
            # sum_n c_n sum_k (-1)^k x^{2k}/((2k)! (n+2k+1))
            s = 0; k = 0; term_x = mp.mpf(1)
            while True:
                t = 0
                for m, cm in enumerate(c):
                    t += cm / (m + 2 * k + 1)
                add = (-1)**k * term_x * t
                s += add
                if k > 10 and abs(add) < mp.mpf(10)**(-38) * (1 + abs(s)):
                    break
                k += 1
                term_x = term_x * x * x / ((2 * k - 1) * (2 * k))
            return s
        return phi_smooth(x) + phi_osc(x)
    return phi, phi_smooth, d0, d1

GLN = 16
def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp; x -= dx
            if abs(dx) < mp.mpf(10)**(-45):
                break
        xs.append(x); ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws
glx, glw = gl_nodes(GLN)

def gl(f, a, b):
    h = (b - a) / 2; m = (a + b) / 2
    return h * sum(w * f(m + h * x) for x, w in zip(glx, glw))

def Omega(k, beta, phi_pack):
    phi, phi_smooth, _, _ = phi_pack
    beta = mp.mpf(beta); cb = mp.cosh(beta)
    def ff(x):
        return phi(x) if x <= 3000 else phi_smooth(x)
    K = lambda t: 1 / mp.sqrt(2 * mp.sinh((beta + t) / 2) * mp.sinh((beta - t) / 2))
    # panels on [0, beta/2]
    tot = 0
    a = mp.mpf(0); half = beta / 2
    while a < half:
        x = Lam * a
        step = (mp.pi / 2) / Lam if x < 3000 else max((mp.pi / 2) / Lam, a / 8)
        b = min(a + step, half)
        tot += gl(lambda t: ff(Lam * t) * K(t), a, b)
        a = b
    # [beta/2, beta]: t = beta - u^2, u in [0, sqrt(beta/2)];  dt K = 2u du / sqrt(2 sinh((2beta-u^2)/2) sinh(u^2/2))
    U = mp.sqrt(half)
    def g(u):
        t = beta - u * u
        if u == 0:
            return 2 * ff(Lam * t) / mp.sqrt(mp.sinh(beta))
        return ff(Lam * t) * 2 * u / mp.sqrt(2 * mp.sinh((beta + t) / 2) * mp.sinh(u * u / 2))
    # panel size in u: d(Lam t)/du = 2 Lam u ; keep phase per panel <= pi/2 when Lam t <= 3000
    a = mp.mpf(0)
    while a < U:
        if Lam * (beta - a * a) <= 3000:
            step = min(U / 64, (mp.pi / 2) / (2 * Lam * max(a, mp.mpf(10)**-8)))
        else:
            step = U / 64
        b = min(a + step, U)
        tot += gl(g, a, b)
        a = b
    return mp.sqrt(2) / mp.pi * tot

def M1(beta):
    cb = mp.cosh(beta)
    return mp.quad(lambda nu: 2 * nu * mp.legenp(-0.5 + 1j * nu, 0, cb, type=3).real / (mp.exp(2 * mp.pi * nu) + 1), [0, 1, 3, 8, 20])

if __name__ == '__main__':
    t0 = time.time()
    packs = {k: phi_factory(psi[k]) for k in (2, 3)}
    a1 = {k: packs[k][2][1] for k in (2, 3)}; a2 = {k: packs[k][2][2] for k in (2, 3)}
    log(f'Lam = {mp.nstr(Lam, 12)}; a1 = psi\'(0): {mp.nstr(a1[2], 7)} (k=2), {mp.nstr(a1[3], 7)} (k=3); a2 = psi\'\'(0): {mp.nstr(a2[2], 8)}, {mp.nstr(a2[3], 4)}')
    A1 = 2 * a1[2] - a1[3]; A2 = 2 * a2[2] + abs(a2[3])
    log(f'A1 = 2a1(2) - a1(3) = {mp.nstr(A1, 7)};  -A2 = -(2 a2(2) + |a2(3)|) = {mp.nstr(-A2, 7)};  int psi2 = {mp.nstr(sum(c/(n+1) for n,c in enumerate(psi[2])), 8)}')
    # (iii) Hadamard identity of step 3:  (sqrt2/pi)[K0/beta - int_0^beta (K - K0) t^-2 dt] = M1(beta)
    for beta in (mp.mpf('0.1'), mp.mpf(1), mp.mpf(4), mp.mpf(10)):
        cb = mp.cosh(beta); K0 = 1 / mp.sqrt(cb - 1)
        Kf = lambda t: 1 / mp.sqrt(2 * mp.sinh((beta + t) / 2) * mp.sinh((beta - t) / 2))
        # K - K0 = 2 sinh^2(t/2) / ((c - cosh t)(c - 1)(K + K0)), stable at small t
        def integrand(t):
            if t == 0:
                return K0 / (4 * (cb - 1))
            Kt = Kf(t)
            return 2 * mp.sinh(t / 2)**2 / ((cb - mp.cosh(t)) * (cb - 1) * (Kt + K0)) / t**2
        I = mp.quad(integrand, [0, beta / 2, beta * mp.mpf('0.9'), beta])
        lhs = mp.sqrt(2) / mp.pi * (K0 / beta - I)
        log(f'(iii) beta={mp.nstr(beta,3)}: Hadamard form {mp.nstr(lhs, 15)}  M1 = {mp.nstr(M1(beta), 15)}')
    # control of Omega against an independent nu-quadrature at small beta (double precision is fine there)
    for beta in (mp.mpf('0.3'), mp.mpf(1)):
        cb = mp.cosh(beta)
        mp.mp.dps = 20
        direct = mp.quad(lambda p: sum(c * p**n for n, c in enumerate(psi[2])) * mp.legenp(-0.5 + 1j * Lam * p, 0, cb, type=3).real,
                         mp.linspace(0, 1, 401))
        mp.mp.dps = 40
        log(f'control beta={mp.nstr(beta,2)}: Omega2 via Mehler-t {mp.nstr(Omega(2, beta, packs[2]), 12)}  vs direct p-quadrature (legenp) {mp.nstr(direct, 12)}  [{time.time()-t0:.0f}s]')
    # main grid
    rtab = [(1.0, 0.577), (1.2, 0.225), (1.5, 0.070), (2.0, 0.015), (3.0, 1.6e-3), (5.0, 1.0e-4), (64, 1.6e-7)]
    def r_claim(beta):
        r = None
        for b0, rv in rtab:
            if beta >= b0 - 1e-12:
                r = rv
        return r
    betas = [1.0, 1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 32.0]
    if len(sys.argv) > 1:
        betas = [float(x) for x in sys.argv[1:]]
    log('\n beta   bhat/Lam   Omega2            two-term(2)       E2/two-term    Omega3            two-term(3)        |E3|/(a1(3)M1/L^2)   r_act=(2|E2|+|E3|)/(A1 M1/L^2)  r_claim   F_LR,min=2Om2-|Om3|   bound(1-r)A1M1/L^2+2047.1/(2bh^3)  F/bound  [s]')
    worst = []
    for beta in betas:
        b = mp.mpf(beta)
        bhat = 2 * Lam * mp.sinh(b / 2)
        m1 = M1(b)
        Om = {k: Omega(k, b, packs[k]) for k in (2, 3)}
        tt = {k: a1[k] * m1 / Lam**2 - a2[k] / (2 * bhat**3) for k in (2, 3)}
        E = {k: Om[k] - tt[k] for k in (2, 3)}
        r_act = (2 * abs(E[2]) + abs(E[3])) / (A1 * m1 / Lam**2)
        F = 2 * Om[2] - abs(Om[3])
        rc = r_claim(beta)
        bound = (1 - rc) * A1 * m1 / Lam**2 + mp.mpf('2047.1') / (2 * bhat**3)
        worst.append((beta, r_act, rc, F, bound))
        log(f'{beta:5.2f}  {mp.nstr(bhat/Lam,6):>8}  {mp.nstr(Om[2],10):>16}  {mp.nstr(tt[2],10):>16}  {mp.nstr(E[2]/tt[2],3):>10}  {mp.nstr(Om[3],8):>14}  {mp.nstr(tt[3],8):>14}  {mp.nstr(abs(E[3])/(a1[3]*m1/Lam**2),3):>10}  {mp.nstr(r_act,3):>12}  {rc:8.2e}  {mp.nstr(F,8):>14}  {mp.nstr(bound,8):>14}  {mp.nstr(F/bound,6):>8}  [{time.time()-t0:.0f}]')
    ok = all(F > 0 and F >= bound and ra <= rc for (_, ra, rc, F, bound) in worst)
    log(f'\nall grid points: F_LR,min > 0, F_LR,min >= stated bound, r_act <= stated r : {ok}')
