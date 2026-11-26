"""rv2_deshift_check.py -- numerical construction of q, chi = B q (note Sec.1) at a test Lambda and direct check of
   R_-(x-i)chi(x-i) + R_+(x+i)chi(x+i) = f(x) = psi2(|x|/Lambda)
by (A) limits from inside the strip y = 1-eps -> 1 (no PV formula used), (B) the stated PV form h.
Usage: python3 rv2_deshift_check.py LAMBDA2(rational, e.g. 900 or 997002) dps"""
import sys
import mpmath as mp
from fractions import Fraction as Fr
from rv_common import *

lam2 = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(900)
mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
L = mp.sqrt(mp.mpf(lam2.numerator)/lam2.denominator)
WIN = mp.mpf(70)   # sech(pi s/2) < 1e-47 beyond |s| > 70

def f(x):
    x = abs(x)
    return psi2_fact_mp(x/L) if x <= L else mp.mpf(0)

def W(v):
    return f(v)/mp.re(Fz(mp.mpf(v)))

def bps(x, extra=()):
    lo, hi = max(-L, x-WIN), min(L, x+WIN)
    pts = {lo, hi}
    # dense breakpoints (every 1/2 unit): the integrand varies on scale ~1 (sech) and steeply near +-L
    k0 = int(mp.ceil(2*lo)); k1 = int(mp.floor(2*hi))
    for k in range(k0, k1+1):
        p = mp.mpf(k)/2
        if lo < p < hi:
            pts.add(p)
    for p in [mp.mpf(0)] + list(extra):
        if lo < p < hi:
            pts.add(p)
    return sorted(pts)

def q_strip(x, y):
    """q(x+iy) = (1/4) int W(v) sech(pi(x+iy-v)/2) dv, 0<=|y|<1, by subtraction of W(x) at v=x."""
    x = mp.mpf(x); y = mp.mpf(y)
    z = mp.mpc(x, y)
    K = lambda v: mp.sech(mp.pi*(z-v)/2)
    eps = 1-abs(y)
    Wx = W(x) if abs(x) < L else mp.mpf(0)
    extra = [x] + [x+s*eps*10**j for j in range(0, 12) for s in (-1, 1) if eps*10**j < WIN]
    pts = bps(x, extra)
    I1 = mp.quad(lambda v: (W(v)-Wx)*K(v), pts, maxdegree=10)
    # int_R K = 2 exactly for |y|<1; subtract the part outside [-L,L] (only matters if window reaches +-L)
    out = mp.mpf(0)
    if x+WIN > L:
        out += mp.quad(K, [L, max(L, x)+WIN])
    if x-WIN < -L:
        out += mp.quad(K, [min(-L, x)-WIN, -L])
    return (I1 + Wx*(2-out))/4

def h_pv(x):
    """h(x) = (1/4) PV int W(v) csch(pi(x-v)/2) dv (the note's Sec.5 form), by subtraction of W(x)."""
    x = mp.mpf(x)
    Wx = W(x) if abs(x) < L else mp.mpf(0)
    K = lambda v: mp.csch(mp.pi*(x-v)/2)
    pts = bps(x, [x])
    I1 = mp.quad(lambda v: (W(v)-Wx)*K(v) if v != x else mp.mpf(0), pts, maxdegree=10)
    # PV int_{-L}^{L} csch(pi(x-v)/2) dv: symmetric window around x cancels; remaining one-sided piece
    lo, hi = max(-L, x-WIN), min(L, x+WIN)
    d1, d2 = x-lo, hi-x
    if d1 < d2:
        rest = mp.quad(K, [x+d1, x+d2])
    else:
        rest = mp.quad(K, [x-d1, x-d2])
    return (I1 + Wx*rest)/4

def lhs_strip(x, eps):
    y = 1-mp.mpf(eps)
    zm, zp = mp.mpc(x, -y), mp.mpc(x, y)
    return Rm(zm)*Bz(zm)*q_strip(x, -y) + Rp(zp)*Bz(zp)*q_strip(x, y)

def lhs_pv(x):
    """boundary values from the PV form q(x-+i) = W/2 +- i h; B(x-+i) evaluated by its *other* form to avoid circularity"""
    x = mp.mpf(x)
    h = h_pv(x); Wx = W(x) if abs(x) < L else mp.mpf(0)
    zm, zp = mp.mpc(x, -1), mp.mpc(x, 1)
    Bm = Bz_alt(zm)          # F(x-2i)/R_+(x-i)
    Bp = Bz(zp)              # F(x+2i)/R_-(x+i)
    return Rm(zm)*Bm*(Wx/2 + 1j*h) + Rp(zp)*Bp*(Wx/2 - 1j*h)

if __name__ == '__main__':
    Lf = float(L)
    xs = [0.31, 1.7, 5.0, 0.5*Lf, 0.9*Lf, Lf-0.5, Lf+3.0]
    print(f"Lambda^2 = {lam2}, Lambda = {mp.nstr(L, 12)}, dps = {mp.mp.dps}")
    print("x | f(x) | LHS(eps=1e-4) | LHS(eps=1e-7) | Richardson(1e-6,1e-7) err | PV-form LHS err | Im parts")
    worst = 0
    for x in xs:
        x = mp.mpf(x)
        fx = f(x)
        l4 = lhs_strip(x, mp.mpf('1e-4'))
        l6 = lhs_strip(x, mp.mpf('1e-6'))
        l7 = lhs_strip(x, mp.mpf('1e-7'))
        rich = l7 + (l7-l6)/9
        lp = lhs_pv(x)
        scale = max(abs(fx), mp.mpf(1))
        e_r = abs(rich-fx)/scale; e_p = abs(lp-fx)/scale
        worst = max(worst, e_p)
        print(f"{mp.nstr(x,8)} | {mp.nstr(fx,12)} | {mp.nstr(mp.re(l4),12)} | {mp.nstr(mp.re(l7),12)} | {mp.nstr(e_r,3)} | {mp.nstr(e_p,3)} | {mp.nstr(abs(mp.im(l7)),2)},{mp.nstr(abs(mp.im(lp)),2)}", flush=True)
    print("worst PV-form relative error:", mp.nstr(worst, 3))
    # relative check beyond Lambda, where f = 0 and both sides are tiny: compare the strip-limit value with the PV-form value
    for x in [Lf+0.5, Lf+3.0, Lf+25.0]:
        x = mp.mpf(x)
        l6 = lhs_strip(x, mp.mpf('1e-6')); l7 = lhs_strip(x, mp.mpf('1e-7')); rich = l7 + (l7-l6)/9
        print(f"beyond Lambda: x-L={mp.nstr(x-L,4)}: strip-limit LHS = {mp.nstr(mp.re(rich),10)} (should be 0); |LHS|/|R_-(x-i)chi(x-i)| = {mp.nstr(abs(rich)/abs(Rm(mp.mpc(x,-1))*Bz(mp.mpc(x,-1))*q_strip(x,-(1-mp.mpf('1e-7')))),3)}")
