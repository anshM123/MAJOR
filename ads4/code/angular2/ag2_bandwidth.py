"""ANGULAR_GO Sec. 2-3: does the angular (smeared-harmonic) positivity survive a wider frequency band?
For the stretched CDH density w_s(nu) = 8 nu^2 / (s^3 cosh(pi nu/s)) (s = 1 is the authors' weight), compute the smeared harmonics
f_m(s) = int w_s Q_{m-1}(nu^2), (d,l) = (0,0), and z_m^raw = -12 f_m + (m-1) f_{m-1} - (m+1) f_{m+1}  (ANGULAR_ROUTE Sec. 3 kernel identity,
applied pointwise in nu).  Single sign of z_m for all m >= 2 is what the telescoping needs.  Also the light side of
chi_s = -nu^2 (4nu^2+9)(4nu^2+49) / (32 (4nu^2+1) cosh(pi nu/s)):  Einstein I_E = int chi_s E_LR and WxW coefficient of |c|^2 = t4^2/16.
mpmath, 60 digits [NE].  Usage: python3 -B ag2_bandwidth.py -> cache/ag2_bandwidth.txt"""
import mpmath as mp
from ag_common import Qrow, write
mp.mp.dps = 60
out = []
def pr(x):
    print(x, flush=True); out.append(x)
MMAX = 40
N10 = lambda n: 1024*n**10 - 195328*n**8 + 8229504*n**6 - 90780512*n**4 + 251603732*n**2 - 114141375
for s in ('1', '1.5', '2', '3', '5'):
    s = mp.mpf(s)
    w = lambda n: 8*n**2/(s**3*mp.cosh(mp.pi*n/s))
    pts = [0, s, 3*s, 8*s, 20*s, 60*s]
    acc = [mp.mpf(0)]*(MMAX + 1)
    # simple composite Gauss-Legendre
    import numpy as np
    gx, gw = np.polynomial.legendre.leggauss(120)
    for a, b in zip(pts[:-1], pts[1:]):
        for xg, wg in zip(gx, gw):
            n = (b - a)/2*mp.mpf(xg) + (b + a)/2
            q = Qrow(0, 0, n**2, MMAX)
            wt = (b - a)/2*mp.mpf(wg)*w(n)
            for j in range(MMAX + 1):
                acc[j] += wt*q[j]
    f = [mp.mpf(0)] + acc   # f[m] = int w Q_{m-1}
    z = {m: -12*f[m] + (m - 1)*f[m - 1] - (m + 1)*f[m + 1] for m in range(2, MMAX)}
    signs = ''.join('-' if z[m] < 0 else '+' for m in range(2, MMAX))
    chi = lambda n: -n**2*(4*n**2 + 9)*(4*n**2 + 49)/(32*(4*n**2 + 1)*mp.cosh(mp.pi*n/s))
    IE = mp.quad(lambda n: n**2*(16*n**4 - 440*n**2 + 273)/((4*n**2 + 1)*mp.cosh(mp.pi*n/s)), pts)
    IW = mp.quad(lambda n: n**2*N10(n)/(117964800*mp.cosh(mp.pi*n/s)), pts)*16   # coefficient of |c|^2 = t4^2/16
    IWflat = mp.quad(lambda n: n**2*1024*n**10/(117964800*mp.cosh(mp.pi*n/s)), pts)*16
    pr('s=%-4s f_1=%s  z_m^raw signs m=2..%d: %s' % (mp.nstr(s, 3), mp.nstr(f[1], 8), MMAX - 1, signs))
    pr('       Einstein I_E = %s ; WxW coeff of |c|^2 = %s (flat nu^10 term alone: %+.4g) ; ratio |W|/I_E = %.4g' %
       (mp.nstr(IE, 12), mp.nstr(IW, 12), float(IWflat), float(abs(IW)/IE)))
write('ag2_bandwidth.txt', out)
