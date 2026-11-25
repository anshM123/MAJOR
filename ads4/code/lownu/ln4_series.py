"""Check of the reviewer's analytic singular de-shift term (k = 2, EE), 40-digit mpmath Taylor coefficients.
E_sing(nu) = 2 Re[R_-(nu-i) chi_sing(nu-i)], chi_sing(z) = c w log w (1 + A w) e^{-z^2/2}, w = z^2+1 = nu^2 - 2 i nu, c = -25 kappa/(36 pi sqrt e).
For nu > 0: log w = log(2 nu) - i pi/2 + log(1 + i nu/2).  So E_sing = g(nu) log(2nu) + e(nu) with g = 2Re h, e = 2Re[h (log(1+i nu/2) - i pi/2)],
h(nu) = c R_-(nu-i) w (1 + A w) e^{-(nu-i)^2/2} (analytic at 0).  Both prefactors P(z) = z^2(z^2+1)/(z^2+1/4)^q, q = 1, 2, are also checked
for E(0) = 0 (they carry w = z^2+1, a simple zero at nu = 0)."""
import mpmath as mp
from ln_common import guard
mp.mp.dps = 40
I = mp.mpc(0, 1)
def Rm(z): return -2*(2*z-I)*(2*z+I)**2*(2*z+3*I)*(2*z+5*I)*(2*z+11*I)/(z*(2*z-7*I)*(2*z-3*I)*(2*z+7*I))
out = ['R_-(-i) = %s' % mp.nstr(Rm(-I), 20), "R_-'(-i)/R_-(-i) = %s" % mp.nstr(mp.diff(Rm, -I) / Rm(-I), 20)]
kappa = mp.mpf(1); c = -25*kappa/(36*mp.pi*mp.sqrt(mp.e))
def coeffs(A):
    h = lambda n: c*Rm(n-I)*((n-I)**2+1)*(1 + A*((n-I)**2+1))*mp.exp(-(n-I)**2/2)
    e = lambda n: h(n)*(mp.log(1 + I*n/2) - I*mp.pi/2)
    th = mp.taylor(h, 0, 5); te = mp.taylor(e, 0, 5)
    return [2*mp.re(x) for x in th], [2*mp.re(x) for x in te]
g0, e0 = coeffs(0); g1, e1 = coeffs(1)
Astar = -g0[2] / (g1[2] - g0[2])
out.append('A solving nu^2 log(nu) cancellation = %s   (37/36 = %s)' % (mp.nstr(Astar, 20), mp.nstr(mp.mpf(37)/36, 20)))
g, e = coeffs(mp.mpf(37)/36)
for p in range(6):
    out.append('nu^%d: log(2nu) coeff %s ; plain coeff %s' % (p, mp.nstr(g[p], 12), mp.nstr(e[p], 12)))
out.append('=> E(0) = %s, E\'(0) = %s (kappa = 1); log terms start at nu^4 iff log coeffs of nu^1..nu^3 vanish' % (mp.nstr(e[0], 5), mp.nstr(e[1], 15)))
txt = '\n'.join(out); print(txt)
open(guard('cache/ln4_series.txt'), 'w').write(txt + '\n')
