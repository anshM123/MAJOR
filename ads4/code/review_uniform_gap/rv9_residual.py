"""rv9_residual.py -- check (4): size of the shifted cutoff defect r(x) = (Re theta(x-i)-1) f(x) - 2 Im theta(x-i) G(x),
G = F h, theta = erf cutoff (C = L+2000, width 10), computed directly (mpmath) for the new family at Lambda^2 = 997002.
Also checks the identity r = E_theta - f by evaluating E_theta = 2 Re[R_-(x-i) chi(x-i) theta(x-i)] from the strip limit."""
import sys
import mpmath as mp
sys.argv = [sys.argv[0], '997002', '40']
import rv2_deshift_check as D
from rv_common import Fz, Rm, Bz
mp.mp.dps = 40
L = D.L; C = L + 2000
def theta(z):
    # computed as 1 - (1-theta) with erfc: the erf-difference form loses Im theta to rounding
    return 1 - one_minus_theta(z)
def one_minus_theta(z):
    return (mp.erfc((z+C)/10) + mp.erfc((C-z)/10))/2
def G(x):
    x = mp.mpf(x)
    if x < L:
        return mp.re(Fz(x))*D.h_pv(x)
    # x >= L: no singularity; h = (1/4) int_{-L}^{L} W(v) csch(pi(x-v)/2) dv, dominated by v near L
    lo = max(-L, L-120)
    pts = [lo + k*(L-lo)/240 for k in range(241)]
    hv = mp.quad(lambda v: D.W(v)*mp.csch(mp.pi*(x-v)/2), pts)/4
    return mp.re(Fz(x))*hv
def r(x):
    x = mp.mpf(x); z = mp.mpc(x, -1)
    th = theta(z); omt = one_minus_theta(z)
    return -mp.re(omt)*D.f(x) + 2*mp.im(omt)*G(x)
print("x-L | log10|G(x)| | log10|Im theta(x-i)| | log10|r(x)|")
best = None
for s in [-900, -500, -100, -10, 0.5, 10, 100, 500, 1000, 1500, 1800, 1900, 1950, 2000, 2050, 2100, 2300]:
    x = L + s
    g = G(x); th = theta(mp.mpc(x, -1)); rr = r(x)
    lg = mp.log10(abs(rr)) if rr != 0 else -mp.inf
    print(f"{s} | {mp.nstr(mp.log10(abs(g)),6)} | {mp.nstr(mp.log10(abs(mp.im(one_minus_theta(mp.mpc(x,-1))))),6)} | {mp.nstr(lg,6)}", flush=True)
    if best is None or lg > best[0]:
        best = (lg, s)
print("max log10|r| on the sampled points:", mp.nstr(best[0], 6), "at x-L =", best[1], " (2^-1000 = 10^", mp.nstr(-1000*mp.log10(2), 6), ")")
# identity check at one point with a moderate cutoff (C' = L+30, width 10) where r is not astronomically small:
C = L + 30
x = L + 25
Es = []
for ep in ['1e-8', '1e-9']:
    y = 1-mp.mpf(ep); zm = mp.mpc(x, -y); zp = mp.mpc(x, y)
    Es.append(Rm(zm)*Bz(zm)*D.q_strip(x, -y)*theta(zm) + D.Rp(zp)*Bz(zp)*D.q_strip(x, y)*theta(zp))
E = mp.re(Es[1] + (Es[1]-Es[0])/9)
print("identity check with C=L+30 at x=L+25: E_theta - f =", mp.nstr(E - D.f(x), 10), "  r formula =", mp.nstr(r(x), 10))
