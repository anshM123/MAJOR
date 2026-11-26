"""rv8_bounds.py -- check (4)/(5): numerical audit of the elementary constants of the note's Secs.4-5 (NOT proofs):
 F bounds, |F(x+-i)|, |B|, Binet Gamma bounds, shifted factor, ladder coefficient bounds of Sec.4, M1 checks of Sec.5."""
import math
import mpmath as mp
from flint import fmpq_poly, fmpq, acb, arb, ctx
from rv_common import Fz, Bz, Rm
mp.mp.dps = 30
out = []
def rep(s):
    print(s, flush=True); out.append(s)

xs = [mp.mpf(0)] + [mp.mpf(10)**(e/20) for e in range(-60, 141)]   # 1e-3 .. 1e7
# (i) F(x) vs (1+x)^(23/2)
r = [mp.re(Fz(x))/(1+x)**mp.mpf(11.5) for x in xs]
rep(f"(4.1) F(x)/(1+x)^(23/2) on [0,1e7]: min {mp.nstr(min(r),6)} (claim >= 2^-11 = {2**-11:.3e}), max {mp.nstr(max(r),6)} (claim <= 2^13 = {2**13})")
r2 = [max(abs(Fz(x+1j)), abs(Fz(x-1j)))/(1+x)**mp.mpf(11.5) for x in xs]
rep(f"(4.2) |F(x+-i)|/(1+x)^(23/2): max {mp.nstr(max(r2),6)} (claim <= 2^14 = {2**14})")
r3 = [abs(Bz(x))/(1+x)**mp.mpf(9.5) for x in xs if x > 0]
rep(f"(4.3) |B(x)|/(1+x)^(19/2): max {mp.nstr(max(r3),6)} (claim <= 512); B(0+) -> {mp.nstr(abs(Bz(mp.mpf('1e-8'))),3)}")
# large-x limits
X = mp.mpf(10)**12
rep(f"      large-x limits at 1e12: F/(1+x)^11.5 = {mp.nstr(mp.re(Fz(X))/(1+X)**mp.mpf(11.5),8)}, |B|/(1+x)^9.5 = {mp.nstr(abs(Bz(X))/(1+X)**mp.mpf(9.5),8)}")
# (ii) Binet bounds: pi e^{-2a-1/(6a)} (a^2+y^2)^{a-1/2} <= cosh(pi y)|Gamma(a+iy)|^2 <= 2 pi e^{1/(6a)} (a^2+y^2)^{a-1/2}
worst_lo = worst_hi = mp.inf
for a in [mp.mpf(3)/4, mp.mpf(5)/4, mp.mpf(7)/4]:
    for y in [mp.mpf(0)] + [mp.mpf(10)**(e/10) for e in range(-30, 71)]:
        v = mp.cosh(mp.pi*y)*abs(mp.gamma(a+1j*y))**2
        base = (a*a+y*y)**(a-mp.mpf(1)/2)
        lo = mp.pi*mp.exp(-2*a-1/(6*a))*base; hi = 2*mp.pi*mp.exp(1/(6*a))*base
        worst_lo = min(worst_lo, v/lo); worst_hi = min(worst_hi, hi/v)
rep(f"(4.4) Binet bounds a in {{3/4,5/4,7/4}}, y in [0,1e7]: min value/lower = {mp.nstr(worst_lo,6)} (>=1 ok), min upper/value = {mp.nstr(worst_hi,6)} (>=1 ok)")
# (iii) shifted factor |cosh(pi(x+-i)/2) Gamma(5/4+i(x+-i)/2) Gamma(5/4-i(x+-i)/2)| <= 16 (1+x)^(3/2)
mx = 0
for x in xs:
    for s in (1, -1):
        z = x+s*1j
        v = abs(mp.cosh(mp.pi*z/2)*mp.gamma(mp.mpf(5)/4+1j*z/2)*mp.gamma(mp.mpf(5)/4-1j*z/2))/(1+x)**mp.mpf(1.5)
        mx = max(mx, v)
rep(f"(4.5) shifted Gamma/cosh factor /(1+x)^(3/2): max {mp.nstr(mx,6)} (claim <= 16)")
# (iv) (1+s)^12 e^{-s/2} <= 24^12 ; (1+s)^13 e^{-s/2} <= 26^13 < 2^62 ; csch(pi s/2) <= 3 e^{-1.5 s} for s >= 1
m12 = max((1+s)**12*math.exp(-s/2) for s in [i/100 for i in range(0, 20000)])
m13 = max((1+s)**13*math.exp(-s/2) for s in [i/100 for i in range(0, 20000)])
cs = max(1/math.sinh(math.pi*s/2)/(3*math.exp(-1.5*s)) for s in [1+i/100 for i in range(0, 20000)])
rep(f"(4.6) max (1+s)^12 e^-s/2 = {m12:.4e} <= 24^12 = {24.0**12:.4e}: {m12 <= 24.0**12};  max (1+s)^13 e^-s/2 = {m13:.4e} <= 26^13 = {26.0**13:.4e} < 2^62 = {2.0**62:.4e}: {m13 <= 26.0**13 < 2.0**62}; max csch/(3e^-1.5s) on s>=1: {cs:.4f}")

# (v) ladder-coefficient bounds of Sec.4 at nu +- i:  (k!)^2 t^k |g[k,3]((nu+-i)^2)| <= (1+y^2)^k (k<=65), y = (nu+2)/L
#     and t^k |g[k,3]((nu+-i)^2)| <= 2 exp(86 y) 1600^-k (all k), checked for k <= 250
def gpolys_flint(p, N):
    Y = fmpq_poly([0, 1]); out_ = [fmpq_poly([1])]; prev = fmpq_poly([0])
    for k in range(N):
        cur = out_[-1]
        nxt = (((8*k*k+8*p*k+4*p-1) - 4*Y)*cur - 4*(k-1+p)**2*prev)*fmpq(1, 4*(k+1)**2)
        prev = cur; out_.append(nxt)
    return out_
ctx.prec = 2000
G3 = gpolys_flint(3, 251)
def ev(poly, w):
    cs_ = poly.coeffs(); s = acb(0)
    for c in reversed(cs_):
        s = s*w + acb(arb(c))
    return s
for L2 in [997002, 10**8]:
    L = math.sqrt(L2)
    worst1 = -1e300; worst2 = -1e300; arg1 = arg2 = None; worst1b = -1e300; arg1b = None
    nus = [0, 0.5, 2, 10, 50, 200, 500, 900, 998, 1100, 1500, 2000, 3000, 5000, 1e4, 3e4, 1e5, 1e6]
    for nu in nus:
        y = (nu+2)/L
        for s in (1, -1):
            w = acb(nu, s)**2
            for k in range(0, 251):
                v = abs(ev(G3[k], w))
                lv = float((v.mid()).log()) if v > 0 else -1e300
                if k <= 65:
                    lhs = 2*math.lgamma(k+1) - k*math.log(L2) + lv
                    rhs = k*math.log1p(y*y)
                    if lhs-rhs > worst1:
                        worst1, arg1 = lhs-rhs, (nu, s, k)
                    if k >= 1 and lhs-rhs > worst1b:
                        worst1b, arg1b = lhs-rhs, (nu, s, k)
                lhs2 = -k*math.log(L2) + lv
                rhs2 = math.log(2) + 86*y - k*math.log(1600)
                if lhs2-rhs2 > worst2 or arg2 is None:
                    worst2, arg2 = lhs2-rhs2, (nu, s, k)
    rep(f"(4.7) Lambda^2={L2}: max log[(k!)^2 t^k|g_k3| / (1+y^2)^k] over k<=65, nu grid = {worst1:.4f} at {arg1} (claim <= 0); over 1<=k<=65: {worst1b:.4f} at {arg1b}")
    rep(f"(4.8) Lambda^2={L2}: max log[t^k|g_k3| / (2 e^(86y) 1600^-k)] over k<=250, nu grid = {worst2:.4f} at {arg2} (claim <= 0)")

# (vi) M1(beta) = (sqrt2/pi) int_0^beta K(t) c(t) dt, K = (cosh b - cosh t)^-1/2, c = cosh(t/2)/(4 sinh^2(t/2)) - 1/t^2
mp.mp.dps = 30
def M1(b):
    b = mp.mpf(b)
    def c(t):
        if t < mp.mpf('1e-4'):
            return mp.mpf(1)/24 - 7*t*t/5760*0 + t*t*(mp.mpf(1)/1920 - mp.mpf(1)/240*0)  # series not needed below
        return mp.cosh(t/2)/(4*mp.sinh(t/2)**2) - 1/t**2
    # substitution t = b(1-v^2) removes the endpoint singularity
    f = lambda v: c(b*(1-v*v))*2*b*v/mp.sqrt(2*mp.sinh(b*(1-v*v/2))*mp.sinh(b*v*v/2)) if v > 0 else c(b)*2*b/mp.sqrt(b*mp.sinh(b))
    return mp.sqrt(2)/mp.pi*mp.quad(f, [0, 0.5, 1])
bs = [1.75 + i*(12-1.75)/200 for i in range(201)]
vals = [(M1(b), b) for b in bs]
mn = min(vals)
rep(f"(5.1) min M1 on beta in [7/4,12] (201 pts): {mp.nstr(mn[0],6)} at beta={mn[1]:.3f} (claim > 9e-5)")
rat = []
for b in [12, 14, 16, 20, 30, 40, 60, 100]:
    K0 = 1/mp.sqrt(mp.cosh(b)-1)
    rat.append((M1(b)/(mp.sqrt(2)/mp.pi*K0/(10*b)), b))
rep(f"(5.2) M1(beta)/[(sqrt2/pi) K0/(10 beta)] for beta in 12..100: {[(b, mp.nstr(r,5)) for r, b in rat]} (claim >= 1)")
# (vii) beta*eps1 <= 1 for beta >= 12, eps1 = (sinh(b/2)-b/2)/(2(cosh b - cosh(b/2)))
e1 = max(b*(math.sinh(b/2)-b/2)/(2*(math.cosh(b)-math.cosh(b/2))) for b in [12+i/10 for i in range(0, 3000)])
rep(f"(5.3) max beta*eps1(beta) on [12,312]: {e1:.4f} (claim <= 1)")
open('out/rv8_bounds.log', 'w').write('\n'.join(out)+'\n')
