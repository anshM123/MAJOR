r"""NON-RIGOROUS numerical sanity checks (floating point; NOT part of any proof).

 1. the classical integral representations used (Mehler-Dirichlet, Mehler with coth, Mehler-Sonine) against the
    hypergeometric definition of the conical function / mpmath Bessel;
 2. the closed forms of S and ftil, the three representations of F_M, the kernel lower bound K_low <= K, R_M > 0;
 3. the curvature expansion (Lemma 6) against the direct definition;
 4. the ingredients (I1)-(I5) of the task statement;
 5. grid scans of F_M(beta) for many M (including small M), and of the pointwise kernel K_{f,M};
 6. the bounds |q_n| <= mu_n, mu_n <= 2/pi^(2n) on samples, and consistency of the certified lower bounds recorded in
    logs/*.log with brute-force floating-point values on randomly chosen boxes.
"""
import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from fractions import Fraction as Fr
from expansion import g_coeffs

mp.mp.dps = 30
f = lambda p: (1 - p) ** 3 * p * (65 + 155 * p + 47 * p * p)
P = lambda nu, be: mp.re(mp.hyp2f1(0.5 + 1j * nu, 0.5 - 1j * nu, 1, -mp.sinh(mp.mpf(be) / 2) ** 2))
S = lambda t: 80 / t ** 3 + 6216 / t ** 5 + 33840 / t ** 7 + (1602 / t ** 4 - 35520 / t ** 6) * mp.sin(t) + (12384 / t ** 5 - 33840 / t ** 7) * mp.cos(t)
ftil = lambda t: -65 / t ** 2 - 1338 / t ** 4 + 1680 / t ** 6 + (1602 / t ** 4 - 35520 / t ** 6) * mp.cos(t) + (-12384 / t ** 5 + 33840 / t ** 7) * mp.sin(t)

msinhc = lambda x: mp.sinh(x) / x if x != 0 else mp.mpf(1)
Efun = lambda be, t: msinhc((be + t) / 2) * msinhc((be - t) / 2)          # cosh be - cosh t = (be^2 - t^2)/2 * E



def ext_int(h, be, T, npts):
    """int_be^T h(t) (cosh t - cosh be)^(-1/2) dt; near t = be use t = be + s^2 and
    cosh t - cosh be = 2 sinh((t+be)/2) sinh((t-be)/2) = s^2 sinh(be + s^2/2) sinhc(s^2/2)."""
    near = mp.quad(lambda s: 2 * h(be + s * s) / mp.sqrt(mp.sinh(be + s * s / 2) * msinhc(s * s / 2)), [0, 0.5, 1])
    far = mp.quad(lambda t: h(t) / mp.sqrt(mp.cosh(t) - mp.cosh(be)), mp.linspace(be + 1, T, npts))
    return near + far


print("=== 1. classical integral representations (max abs deviation over samples) ===")
dev1 = dev2 = dev3 = 0
for nu, be in ((0.3, 0.5), (2.0, 1.3), (7.5, 0.4), (1.1, 4.0), (12.0, 2.0)):
    ref = P(nu, be)
    # interior, in the theta-form  P = (2/pi) int_0^{pi/2} cos(nu be sin th) E^{-1/2} dth  (t = be sin th)
    v1 = 2 / mp.pi * mp.quad(lambda th: mp.cos(mp.mpf(nu) * mp.mpf(be) * mp.sin(th)) / mp.sqrt(Efun(mp.mpf(be), mp.mpf(be) * mp.sin(th))), mp.linspace(0, mp.pi / 2, 12))
    v2 = mp.sqrt(2) / mp.pi * mp.coth(mp.pi * nu) * ext_int(lambda t: mp.sin(nu * t), mp.mpf(be), be + 150, 300)
    dev1 = max(dev1, abs(v1 - ref)); dev2 = max(dev2, abs(v2 - ref))
print("  Mehler-Dirichlet  P = (sqrt2/pi) int_0^beta cos(nu t)(cosh beta - cosh t)^(-1/2) dt :", mp.nstr(dev1, 3))
print("  Mehler (coth)     P = (sqrt2/pi) coth(pi nu) int_beta^oo sin(nu t)(cosh t - cosh beta)^(-1/2) dt :", mp.nstr(dev2, 3))
# Mehler-Sonine through its consequence fhat(b) = (2/pi) int_b^oo S(tau) (tau^2-b^2)^(-1/2) dtau
dev3 = 0
for b in (2.0, 9.9, 15.0):
    lhs = mp.quad(lambda q: f(q) * mp.besselj(0, q * b), mp.linspace(0, 1, 6))
    rhs = 2 / mp.pi * (mp.quad(lambda s: 2 * S(b + s * s) / mp.sqrt(2 * b + s * s), [0, 1]) + mp.quad(lambda t: S(t) / mp.sqrt(t * t - b * b), mp.linspace(b + 1, b + 400, 200)) + 0)
    dev3 = max(dev3, abs(lhs - rhs))
print("  flat exterior     fhat(b) = (2/pi) int_b^oo S (tau^2-b^2)^(-1/2) (tail beyond b+400 dropped, ~1e-6):", mp.nstr(dev3, 3))

print("=== 2. closed forms and representations of F_M ===")
d = max(abs(S(t) - mp.quad(lambda q: f(q) * mp.sin(q * t), mp.linspace(0, 1, 8))) for t in (3.3, 9.9, 17.0))
d2 = max(abs(ftil(t) - mp.quad(lambda q: f(q) * mp.cos(q * t), mp.linspace(0, 1, 8))) for t in (3.3, 9.9, 17.0))
print("  S closed form:", mp.nstr(d, 3), "  ftil closed form:", mp.nstr(d2, 3))
D0 = {1: 65, 2: -80, 3: -1338, 4: 6216, 5: -1680, 6: -33840}
D1 = {3: -1602, 4: -12384, 5: -35520, 6: -33840}


def Kexact(M, t):
    return mp.quad(lambda q: f(q) * mp.coth(mp.pi * M * q) * mp.sin(M * q * t) if q > 0 else mp.mpf(0), mp.linspace(0, 1, 12))


def Klow(M, tau, m):
    tot = S(tau)
    for k in range(1, m + 1):
        z = 2 * mp.pi * k * M - 1j * tau
        tot += 2 * mp.im(sum(D0[j] / z ** (j + 1) for j in D0) - mp.exp(-z) * sum(D1[j] / z ** (j + 1) for j in D1))
    return tot


for M, tau in ((1.15, 10.1), (3.0, 7.0), (6.0, 10.12), (6.0, 9.7)):
    ke = Kexact(M, tau / M)
    print(f"  M={M}, tau={tau}:  K exact = {mp.nstr(ke, 12)},  K_low^(4) = {mp.nstr(Klow(M, tau, 4), 12)},  K_low^(200) = {mp.nstr(Klow(M, tau, 200), 12)},  S = {mp.nstr(S(tau), 8)}  (K_low^(4) <= K: {Klow(M, tau, 4) <= ke}; R_M = K - S > 0: {ke - S(tau) > 0})")
_xg, _wg = leggauss(200)
_pg = (_xg + 1) / 2
_fg = f(_pg) * _wg / 2


def F_exterior_numpy(M, be):
    """(sqrt2/pi) int_be^oo K_{f,M}(t) (cosh t - cosh be)^(-1/2) dt with K by Gauss-Legendre in p (double precision)."""
    Kt = lambda t: (_fg[:, None] / np.tanh(np.pi * M * _pg[:, None]) * np.sin(M * _pg[:, None] * t[None, :])).sum(0)
    tot = 0.0
    # near part: t = be + s^2, s in [0,1]
    s = (_xg + 1) / 2
    t = be + s * s
    tot += np.sum(_wg / 2 * 2 * Kt(t) / np.sqrt(np.sinh(be + s * s / 2) * np.sinh(s * s / 2) / (s * s / 2)))
    edges = np.linspace(be + 1, be + 60, 60)
    for a0, a1 in zip(edges[:-1], edges[1:]):
        t = a0 + (a1 - a0) * (_xg + 1) / 2
        tot += np.sum(_wg * (a1 - a0) / 2 * Kt(t) / np.sqrt(np.cosh(t) - np.cosh(be)))
    return np.sqrt(2) / np.pi * tot


mp.mp.dps = 15
ftil_num = lambda t: ftil(t) if t > 2 else mp.quad(lambda q: f(q) * mp.cos(q * t), [0, 1])
for M, be in ((2.0, 1.5), (8.0, 0.9)):
    direct = mp.quad(lambda q: f(q) * P(M * q, be), mp.linspace(0, 1, 8))
    interior = 2 / mp.pi * mp.quad(lambda th: ftil_num(M * be * mp.sin(th)) / mp.sqrt(Efun(be, be * mp.sin(th))), mp.linspace(0, mp.pi / 2, 16))
    print(f"  F_M(beta), M={M}, beta={be}: direct {mp.nstr(direct, 12)}, interior(ftil) {mp.nstr(interior, 12)}, exterior(K, numpy) {F_exterior_numpy(M, be):.10f}")
mp.mp.dps = 30

print("=== 3. curvature expansion sum_n beta^(2n) G_n(M beta), N = 14, K = 36, against the direct definition ===")
N, K = 14, 36
G = g_coeffs(N, K)
for M, be in ((6, 1.8), (6, 1.0), (10, 1.0), (8, 0.3), (50, 0.2)):
    b = mp.mpf(M) * be
    ex = sum(mp.mpf(be) ** (2 * n) * sum(mp.mpf(G[n][k].numerator) / G[n][k].denominator * b ** (2 * k) for k in range(K + 1)) for n in range(N + 1))
    direct = mp.quad(lambda q: f(q) * P(M * q, be), mp.linspace(0, 1, 8))
    print(f"  M={M}, beta={be}: direct {mp.nstr(direct, 15)}  expansion {mp.nstr(ex, 15)}  diff {mp.nstr(ex - direct, 3)}")

print("=== 4. ingredients (I1)-(I5) ===")
w = 0
for nu in (0.0, 0.7, 3.0, 20.0):
    for be in (0.3, 1.0, 2.5):
        lhs = abs(P(nu, be) - mp.besselj(0, nu * be)); mid = 1 - P(0, be); rhs = mp.sinh(mp.mpf(be) / 2) ** 2 / 4
        w = max(w, lhs - mid, mid - rhs)
print("  (I1) |P - J0| <= 1 - P_{-1/2} <= sinh^2(beta/2)/4 : max violation over samples =", mp.nstr(w, 3), "(<= 0 means OK)")
print("  (I2) ftil(tau) = -65/tau^2 - 1338/tau^4 + 1680/tau^6 + (1602/tau^4 - 35520/tau^6) cos tau + (-12384/tau^5 + 33840/tau^7) sin tau : see 2.")
print("       int_0^X ftil = Psi(X) -> 0:  Psi(200) =", mp.nstr(mp.quad(lambda q: f(q) / q * mp.sin(200 * q) if q > 0 else 0, mp.linspace(0, 1, 80)), 6), " (= 65/200 + O(X^-3) = 0.325...; Psi(X) ~ 65/X)")
fh = lambda b: mp.quad(lambda q: f(q) * mp.besselj(0, q * b), mp.linspace(0, 1, int(b / 3) + 4))
for b in (40.0, 80.0, 160.0):
    nonosc = 40 / b ** 3 + 2331 / b ** 5 + 10575 / b ** 7
    print(f"  (I5) b={b}: fhat = {mp.nstr(fh(b), 10)}, 40/b^3+2331/b^5+10575/b^7 = {mp.nstr(nonosc, 10)}, (fhat - nonosc) b^4.5 = {mp.nstr((fh(b) - nonosc) * b ** 4.5, 6)}  (bounded; amplitude 1602 sqrt(2/pi) = {mp.nstr(1602 * mp.sqrt(2 / mp.pi), 6)})")

print("=== 5. grid scans (numpy, double precision) ===")
xp, wp = leggauss(100); pp = (xp + 1) / 2; fpw = f(pp) * wp / 2
xt, wt = leggauss(400); th = (xt + 1) * np.pi / 4; wth = wt * np.pi / 4


def log_sinhc_minus(a):
    """log(sinh(a)/a) - a, overflow-free:  = log((1 - e^{-2a})/(2a))  (-> -a + a^2/6 for small a)."""
    big = a > 1e-4
    asafe = np.where(big, a, 1.0)
    return np.where(big, np.log1p(-np.exp(-2 * asafe)) - np.log(2 * asafe), -a + a * a / 6)


def Fnum(M, b):
    """exp(beta/2) * F_M(beta), beta = b/M, from the theta-form of the Mehler-Dirichlet formula
    F_M = (2/pi) int_0^{pi/2} ftil(b sin th) E^{-1/2} dth,  E = sinhc((beta+t)/2) sinhc((beta-t)/2), t = beta sin th.
    Since (beta+t)/2 + (beta-t)/2 = beta, the factor exp(beta/2) cancels the exponential growth of E^{1/2} exactly."""
    beta = b / M
    t = beta * np.sin(th)
    logE_minus_beta = log_sinhc_minus((beta + t) / 2) + log_sinhc_minus((beta - t) / 2)
    ft = (fpw[:, None] * np.cos(pp[:, None] * (b * np.sin(th))[None, :])).sum(0)
    return (2 / np.pi) * np.sum(wth * ft * np.exp(-0.5 * logE_minus_beta))


bs = np.concatenate([np.linspace(0.02, 30, 750), np.linspace(30, 200, 400)[1:]])
print("  M      min over b=M*beta in (0,200] of  F_M * exp(beta/2) * max(1,b)^2      (argmin b)      [positive = OK]")
for M in (0.02, 0.1, 0.3, 0.6, 1.0, 1.15, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 45, 80, 200, 1000):
    v = np.array([Fnum(M, b) for b in bs]) * np.maximum(1, bs) ** 2
    i = int(np.argmin(v))
    print(f"  {M:<7g} {v[i]: .6e}    b = {bs[i]:.2f}")
bb = np.linspace(0, 10.85, 2171)
fl = np.array([(2 / np.pi) * np.sum(wth * (fpw[:, None] * np.cos(pp[:, None] * (b * np.sin(th))[None, :])).sum(0)) for b in bb])
print(f"  flat: min of fhat on [0,10.85] = {fl.min():.6f} at b = {bb[int(np.argmin(fl))]:.3f}   (certified lower bounds in C4/C5a are smaller, as they must be)")
taus = np.linspace(0.05, 40, 8000)


def Kgrid(M):
    return (fpw[:, None] / np.tanh(np.pi * M * pp[:, None]) * np.sin(pp[:, None] * taus[None, :])).sum(0)


Sg = (fpw[:, None] * np.sin(pp[:, None] * taus[None, :])).sum(0)
neg = taus[Sg < 0]
print(f"  S(tau) < 0 on the grid only for tau in [{neg.min():.3f}, {neg.max():.3f}]; min S = {Sg.min():.5f}")
print("  M      min_tau K_{f,M}(tau/M)   (pointwise kernel positivity holds for M <~ 7.6 and FAILS beyond: the curved theorem for M > 6 needs C4)")
for M in (1.15, 3, 5, 6, 7, 7.5, 7.7, 8, 10, 20, 80):
    k = Kgrid(M); i = int(np.argmin(k))
    print(f"  {M:<6g} {k[i]: .5f}   at tau = {taus[i]:.2f}")

print("=== 6. consistency of the certificate logs with brute-force floating-point values (if logs/ exists) ===")
import os, re, random
from expansion import q_polys, mu_coeffs
random.seed(1)
qs, mus = q_polys(14), mu_coeffs(14)
ss = np.linspace(0, 1, 401)
okq = all(max(abs(sum(float(c) * s ** (2 * j) for j, c in enumerate(qs[n]))) for s in ss) <= float(mus[n]) * (1 + 1e-12) for n in range(15))
okmu = all(float(mus[n]) * np.pi ** (2 * n) <= 2 for n in range(1, 15))
print(f"  Lemma 7: |q_n(s)| <= mu_n on a grid of s, n <= 14: {okq};   mu_n pi^(2n) <= 2, 1 <= n <= 14: {okmu}")
here = os.path.dirname(os.path.abspath(__file__))
num = r"([-+0-9.eE]+)"
p4 = os.path.join(here, "logs", "C4_expansion.log")
if os.path.exists(p4):
    lines = [l for l in open(p4) if l.startswith("PASS box")]
    worst = np.inf
    for l in random.sample(lines, 60):
        b0, b1, x0, x1, low = map(float, re.match(rf"PASS box b in \[{num}, {num}\] beta\^2 in \[{num}, {num}\]\s+FF >= {num}", l).groups())
        mn = np.inf
        for b in np.linspace(b0, b1, 4):
            for x in np.linspace(x0, x1, 9):
                be = np.sqrt(x)
                val = Fnum(b / be, b) * np.exp(-be / 2) if be > 1e-9 else (2 / np.pi) * np.sum(wth * (fpw[:, None] * np.cos(pp[:, None] * (b * np.sin(th))[None, :])).sum(0))
                mn = min(mn, val)
        worst = min(worst, mn - low)
    print(f"  C4: 60 random boxes, min over boxes of (sampled min of F) - (certified lower bound) = {worst:.3e}   (must be >= 0)")
p3 = os.path.join(here, "logs", "C3_kernel.log")
if os.path.exists(p3):
    lines = [l for l in open(p3) if l.startswith("PASS box")]
    worst = np.inf
    for l in random.sample(lines, 200):
        M0_, M1_, t0, t1, low = map(float, re.match(rf"PASS box M in \[{num}, {num}\] tau in \[{num}, {num}\]\s+K_low >= {num}", l).groups())
        mn = min(np.sum(fpw / np.tanh(np.pi * M * pp) * np.sin(pp * t)) for M in np.linspace(M0_, M1_, 3) for t in np.linspace(t0, t1, 3))
        worst = min(worst, mn - low)
    print(f"  C3: 200 random boxes, min of (sampled min of the exact kernel K) - (certified lower bound of K_low) = {worst:.3e}   (must be >= 0)")
p1 = os.path.join(here, "logs", "C1_S.log")
if os.path.exists(p1):
    lines = [l for l in open(p1) if l.startswith("PASS S>0")]
    worst = np.inf
    for l in random.sample(lines, 300):
        a0, a1 = map(float, re.search(rf"cell \[{num}, {num}\]", l).groups())
        worst = min(worst, min(np.sum(fpw * np.sin(pp * t)) for t in np.linspace(a0, a1, 5)))
    print(f"  C1: 300 random cells, min sampled S = {worst:.3e}   (must be > 0)")
