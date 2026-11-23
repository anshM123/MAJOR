"""Target 2 (bound table): evaluate the explicit remainder bound of Lemma LR-infinity as stated in LARGESPIN_L §3,
   |E_k| <= (sqrt2/pi)[T_a + T_b + T_c],  r(beta) := (2|E_2|_bd + |E_3|_bd)/(A1 M1/Lam^2),
and compare sup_{beta' >= beta} r with the stated table (1.0: 0.577, 1.2: 0.225, 1.5: 0.070, 2: 0.015, 3: 1.6e-3, 5: 1.0e-4, 64: 1.6e-7).
Grid as stated: step 0.01 on [0.9, 3], 0.05 to 16, 0.5 to 64.   Output: cache/t2c_rbound.txt
"""
import json, sys
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 30
out = open('cache/t2c_rbound.txt' if len(sys.argv) == 1 else 'cache/t2c_rbound_extra.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

js = json.load(open('../../prl/data/nf02_psi.json'))
k0 = mp.mpf(js['k0'])
psi = {2: [mp.mpf(Fraction(c).numerator) / Fraction(c).denominator / k0 for c in js['psi2_exact']],
       3: [mp.mpf(Fraction(c).numerator) / Fraction(c).denominator / k0 for c in js['psi3_exact']]}
Lam = mp.sqrt(mp.mpf(998) * 999)

def dj(c, j, p):
    return sum(c[n] * mp.ff(n, j) * p**(n - j) for n in range(j, len(c)))

def pack(k):
    c = psi[k]; n = len(c)
    g = c[1:]                                    # g = psi/p (psi(0) = 0)
    d0 = [dj(c, j, 0) for j in range(n)]; d1 = [dj(c, j, 1) for j in range(n)]
    g0 = [dj(g, j, 0) for j in range(len(g))]; g1 = [dj(g, j, 1) for j in range(len(g))]
    absint = mp.quad(lambda p: abs(sum(cc * p**m for m, cc in enumerate(c))), mp.linspace(0, 1, 41))
    def Rinf(x):
        return sum(abs(d0[j]) / x**(j + 1) for j in range(3, n, 2)) + sum(abs(d1[j]) / x**(j + 1) for j in range(5, n))
    def R(x):
        return min(Rinf(x), abs(d0[1]) / x**2 + absint)
    def d1f(X):
        return sum(abs(g0[j]) / X**(j + 1) for j in range(2, len(g), 2)) + sum(abs(g1[j]) / X**(j + 1) for j in range(5, len(g)))
    def int_x2Rinf(X):   # int_X^oo x^2 Rinf(x) dx
        return sum(abs(d0[j]) * X**(2 - j) / (j - 2) for j in range(3, n, 2)) + sum(abs(d1[j]) * X**(2 - j) / (j - 2) for j in range(5, n))
    return dict(a1=d0[1], Rinf=Rinf, R=R, d1=d1f, ix2=int_x2Rinf, absint=absint)

P = {k: pack(k) for k in (2, 3)}
A1 = 2 * P[2]['a1'] - P[3]['a1']

def M1(beta):
    cb = mp.cosh(beta)
    return mp.quad(lambda nu: 2 * nu * mp.legenp(-0.5 + 1j * nu, 0, cb, type=3).real / (mp.exp(2 * mp.pi * nu) + 1), [0, 1, 3, 8, 20])

def Ebound(k, beta):
    p = P[k]
    D = mp.cosh(beta) - 1; K0 = 1 / mp.sqrt(D); X = Lam * beta
    Ta = K0 / Lam * p['d1'](X)
    Tb = K0 / (4 * D * Lam**3) * p['ix2'](X)
    rho_b = lambda t: K0 * (t**4 * mp.cosh(t) / (48 * D) + mp.mpf(3) / 8 * 2**mp.mpf(2.5) * ((mp.cosh(t) - 1) / D)**2)
    # R(Lam t) switches branch once; integrate with breakpoints
    Tc1 = mp.quad(lambda t: p['R'](Lam * t) * rho_b(t), [0, mp.mpf(1) / Lam, mp.mpf(10) / Lam, mp.mpf(100) / Lam, beta / 4, beta / 2])
    Tc2 = p['Rinf'](Lam * beta / 2) * (mp.pi / mp.sqrt(2)) * mp.legenp(-0.5, 0, mp.cosh(beta), type=3)
    return mp.sqrt(2) / mp.pi * (Ta + Tb + Tc1 + Tc2), (Ta, Tb, Tc1, Tc2)

grid = [mp.mpf(90 + i) / 100 for i in range(0, 211)] + [mp.mpf(3) + mp.mpf(i) / 20 for i in range(1, 261)] + [mp.mpf(16) + mp.mpf(i) / 2 for i in range(1, 97)]
if len(sys.argv) > 1:
    grid = [mp.mpf(x) for x in sys.argv[1:]]
rs = []
for b in grid:
    m1 = M1(b)
    e2, parts2 = Ebound(2, b); e3, _ = Ebound(3, b)
    r = (2 * e2 + e3) / (A1 * m1 / Lam**2)
    rs.append((b, r))
    if len(sys.argv) > 1 or (abs(b - mp.nint(b * 4) / 4) < 1e-9 and (b <= 3 or abs(b - mp.nint(b)) < 1e-9)):
        log(f'beta = {mp.nstr(b, 4):>5}: r_bound = {mp.nstr(r, 4)}   (T_a, T_b, T_c1, T_c2 for k=2: {", ".join(mp.nstr(x, 3) for x in parts2)};  M1 = {mp.nstr(m1, 6)})')
log('')
for b0, claim in (() if len(sys.argv) > 1 else ((1.0, 0.577), (1.2, 0.225), (1.5, 0.070), (2.0, 0.015), (3.0, 1.6e-3), (5.0, 1.0e-4), (64, 1.6e-7))):
    sup = max(r for (b, r) in rs if b >= b0 - 1e-12)
    log(f'sup_(beta >= {b0}) r_bound = {mp.nstr(sup, 4)}   stated {claim}')
