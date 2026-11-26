"""lx_bound.py: (a) independent re-implementation vs rt_pp.cross at the [3.65, 3.75] right-edge parameters;
(b) uniform-sum sanity (terms added for absent phases are >= 0); (c) monotonicity in a; (d) exact auxiliary inequalities
(Delta eps >= eps^2/2, tau*Delta(1/y) >= 1, d eps/dk <= 4/m^2, c0 <= 1 + 2|alpha|/a^2) on random states in exact rationals."""
import math, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic')
import rt_pp, crossing
from lx_common import bound_indep, LAMf
from flint import arb
LAMup = float(arb(998*999).sqrt().upper())
print('Lambda =', LAMf, 'Lam_up =', LAMup)
print('\n(a) right edge: tau in [Lam-1, Lam_up], uniform a >= 998')
print(' sigma | rt_pp.cross (c=0.7,1,1.5) | PP0 | indep (c=0.7,1,1.5) | indep phase-II-end only (c best)')
for sg in (0.0, 0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 4.5, 5.5):
    cr = [math.log(rt_pp.cross(sg, LAMf - 1, LAMup, c)) for c in (0.7, 1.0, 1.5)]
    pp0 = math.log(rt_pp.PP0(sg, LAMf - 1, LAMup))
    ind = [bound_indep(sg, LAMf - 1, LAMup, c, parts=True) for c in (0.7, 1.0, 1.5)]
    print(f' {sg:4.2f} | ' + ' '.join(f'{x:7.3f}' for x in cr) + f' | {pp0:7.3f} | ' + ' '.join(f'{x[0]:7.3f}' for x in ind)
          + ' | ' + ' '.join(f'I {x[1]["I"]:.2f} II {x[1]["II"]:.2f} III {x[1]["III"]:.2f};' for x in ind[1:2]))
print('\n(b) pieces of rt_pp.cross (c = 1, sigma = 5.5, tau = [Lam-1, Lam_up]) -- terms added for phases that are not reached must be >= 0')
s = arb(5.5); tl = arb(LAMf - 1); th = arb(LAMup); aa = arb(998)
al = (s*s + th*th).sqrt(); q = al/tl; dl = tl**(arb(-2)/3); yA = 1 - dl; yB = 1 + dl; umax = s/(aa + 1)
eta = 2*(dl + umax); dyx = 2*th/aa**2
e1 = 2*(umax + 1 - yA + dyx + ((1 - (yA - dyx)**2) + umax*umax + 2*umax).sqrt())
W = th*(1/yA - 1/yB) + 3
lII = ((2 + e1) + e1/eta.sqrt()).log() - arb('2.01').log() + W*(eta.sqrt() + eta)
rB = (yB*yB - 1 - umax*umax).sqrt(); e2 = umax + yB + dyx - 1 + ((yB + dyx)**2 - 1 + umax*umax + 2*umax*(yB + dyx)).sqrt()
lP = ((eta.sqrt()*(1 + e2) + e2)/(2*rB)).log()
ym = yB*arb('1.000001'); lmax = ym + (ym*ym - 1).sqrt() + 2*umax
print('  lII =', lII.str(5), ' lP =', lP.str(5), ' lfin =', (2*lmax*lmax).log().str(5), ' (all > 0)')
print('  y_max on the edge = Lam_up/999 =', LAMup/999, ' < yA =', float(yA.mid()), '?', LAMup/999 < float(yA.mid()), ' < 1:', LAMup/999 < 1)
print('\n(c) monotonicity in a (log cross(5.5, Lam-1, Lam_up, c=1, a)):', [round(math.log(rt_pp.cross(5.5, LAMf-1, LAMup, 1.0, a)), 4) for a in (998, 998.5, 999, 1010, 1100, 2000, 4500)])
print('\n(d) exact auxiliary inequalities on random states (Fractions)')
random.seed(1)
bad = 0; cnt = 0; worst_g = 0; worst_sp = 1e9
for t in range(400):
    a = Fr(random.randint(998*8, 4600*8), 8)
    n = random.choice([random.randint(1, 3000), random.randint(1000, 400000)])
    N = n + a
    for k in random.sample(range(0, n - 1), min(20, n - 1)) if n > 2 else []:
        Q = lambda k: (N - k)*(N + k + 1)
        e = lambda k: Fr(2*(2*k + 1))/Q(k)
        ek, ek1 = e(k), e(k + 1)
        cnt += 1
        if not (ek1 - ek >= ek*ek/2): bad += 1
        # tau Delta(1/y) = Q(k+1)/(2k+3) ... spacing of 1/y in units 1/tau: Q(k)/(2k+1) - Q(k+1)/(2k+3) >= 1
        sp = Q(k)/(2*k + 1) - Q(k + 1)/(2*k + 3)
        worst_sp = min(worst_sp, float(sp))
        if sp < 1: bad += 1
        # Delta eps <= 4/(m_{k+1})^2,  m_{k+1} = N - k - 1
        g = (ek1 - ek)*(N - k - 1)**2
        worst_g = max(worst_g, float(g))
        if g > 4: bad += 1
print(f'  {cnt} (state,k) checks: violations = {bad}; max (Delta eps)(N-k-1)^2 = {worst_g:.4f} (<= 4); min tau*Delta(1/y) = {worst_sp:.4f} (>= 1)')
# c0 at n = 0
import cmath
w0max = 0
for (a, n, sg, tau) in [(998, 1, 5.5, 1150), (998, 80000, 5.5, 998.5), (4500, 10**6, 0, 900)]:
    N = n + a; al = complex(sg, tau); w = al/(N*(N + 1))
    r = cmath.sqrt(1 + w*w); lp = w + r; lm = w - r
    c0 = max(abs(1 - lm), abs(lp - 1))/(2*abs(r))
    print(f'  c0(a={a}, n={n}, alpha={al}) = {c0:.12f}  vs 1 + 2|alpha|/a^2 = {1 + 2*abs(al)/a**2:.12f}')
