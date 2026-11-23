"""
lsw/lsw_weyl.py (LARGESPIN_W.md §2.2): the origin of CHMRS's "only even powers of 1/m in the symmetrical frame".

Claim [D]: with T(x, y; Delta, J) the product of three-point functions of 2106.10274 app. A (eq. 'moment appendix', d = 3),
      T(-x, -y; Delta, J) = T(x, y; 3 - Delta, -1 - J)        (Weyl reflection (Delta, J) -> (d - Delta, 2 - d - J)),
and, with m^2 = (Delta-J-2)(Delta+J-1), Delta+J-1 = m e^{beta/2}, Delta-J-2 = m e^{-beta/2}, the reflection is exactly m -> -m at fixed beta.
Hence F(eps) := T(eps X, eps Y; m = 1/eps, beta) satisfies F(-eps) = F(eps): only even powers of 1/m at fixed (X, Y, beta).

Checks (mpmath, 40 digits):
  (1) the identity at random small future-timelike x, y and random (Delta, J) (non-integer J allowed);
  (2) the Taylor coefficients of log F(eps) - log F_0 in eps (from F(eps) at eps = +-h, +-2h, ... with the reflected parameters):
      odd coefficients vanish to working precision;
  (3) the same test with UNEQUAL external dimensions (Delta_12 != 0), i.e. a CPPR seed: odd coefficients do NOT vanish
      (they are proportional to Delta_12, Delta_34) -- the reason the spinning parity statement is only [D-sketch] (§2.3).
Usage: python3 lsw_weyl.py -> cache/weyl.txt
"""
import os, sys
sys.dont_write_bytecode = True
import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
e = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
n = (mp.mpf(1), mp.mpf(0), mp.mpf(1))          # null, along e_z (here the 2-axis)


def dot(a, b):
    return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def add(*vs):
    return tuple(sum(v[i] for v in vs) for i in range(3))


def sc(c, v):
    return tuple(c*vi for vi in v)


def bar(y):                                     # space reflection ybar = -y - 2 e (y.e)  -> (y0, -y1, -y2)
    return (y[0], -y[1], -y[2])


def Vfrm(x, y, ee):
    ex, ey = dot(ee, x), dot(ee, y)
    return add(sc((1 - 2*ex)*(1 - 2*ey) - dot(x, x)*dot(y, y), ee), sc(dot(add(ee, y), add(ee, y)), x),
               sc(dot(add(ee, x), add(ee, x)), bar(y)))


def T(x, y, D, J, Dphi=3, d12=0, d34=0):
    """|T(x,y)| of CHMRS app. A (d = 3); d12, d34: extra seed factors for unequal external dimensions
    (|x15^2/x25^2|^{d12/2} and |x35^2/x45^2|^{d34/2}, x1 = -y, x2 = -x/(-x^2), x3 = y, x4 = x/(-x^2), x5 = e)."""
    d = 3
    x2, y2 = dot(x, x), dot(y, y)
    Dt = d - Dphi
    me = sc(-1, e)
    R = (dot(add(e, sc(-1, x)), add(e, sc(-1, x)))*dot(add(e, sc(-1, y)), add(e, sc(-1, y)))) / \
        (dot(add(e, x), add(e, x))*dot(add(e, y), add(e, y)))
    a = -dot(n, Vfrm(x, y, e)); b = dot(n, Vfrm(x, y, me))
    val = (abs(x2*y2))**(-Dt) * abs(1 + 2*dot(x, y) + x2*y2)**(d - 1 - 2*Dphi) / \
        abs(dot(add(e, sc(-1, x)), add(e, sc(-1, x)))*dot(add(e, sc(-1, y)), add(e, sc(-1, y)))) * \
        abs(R)**((D + J)/2) * abs(a)**J / abs(b)**(J + d - 2)
    if d12 or d34:
        X1 = sc(-1, y); X2 = sc(-1/x2, x); X3 = y; X4 = sc(1/x2, x)
        sq = lambda p, q: abs(dot(add(p, sc(-1, q)), add(p, sc(-1, q))))
        val *= (sq(X1, e)/sq(X2, e))**(mp.mpf(d12)/2) * (sq(X3, e)/sq(X4, e))**(mp.mpf(d34)/2)
    return val


def DJ(m, beta):
    tb, tt = m*mp.e**(beta/2), m*mp.e**(-beta/2)           # Delta+J-1, Delta-J-2
    J = (tb - tt - 1)/2; D = (tb + tt + 3)/2
    return D, J


out = []
pr = lambda s: (print(s, flush=True), out.append(s))
mp.mp.dps = 40
import random
random.seed(7)
pr('(1) identity T(-x,-y; D, J) = T(x, y; 3-D, -1-J):')
for t in range(4):
    x = (mp.mpf(random.uniform(0.02, 0.05)), mp.mpf(random.uniform(-0.01, 0.01)), mp.mpf(random.uniform(-0.01, 0.01)))
    y = (mp.mpf(random.uniform(0.02, 0.05)), mp.mpf(random.uniform(-0.01, 0.01)), mp.mpf(random.uniform(-0.01, 0.01)))
    D, J = mp.mpf(random.uniform(30, 60)), mp.mpf(random.uniform(5, 20))
    lhs = T(sc(-1, x), sc(-1, y), D, J); rhs = T(x, y, 3 - D, -1 - J)
    pr('   rel diff = %.2e' % float(abs(lhs/rhs - 1)))

beta = mp.mpf('0.3')
X = (mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5'))       # future timelike (X^2 = -1.69+0.16+0.25 < 0)
Y = (mp.mpf('0.9'), mp.mpf('-0.2'), mp.mpf('0.35'))


def logF(eps, d12=0, d34=0):
    m = 1/eps
    D, J = DJ(m, beta)
    return mp.log(T(sc(eps, X), sc(eps, Y), D, J, d12=d12, d34=d34))


for (d12, d34) in [(0, 0), (1, 0), (1, -1), (2, 2)]:
    h = mp.mpf('1e-4')
    # odd part of log F via symmetric differences: O(eps) coefficient = (logF(h) - logF(-h))/(2h) + O(h^2)
    c1 = (logF(h, d12, d34) - logF(-h, d12, d34))/(2*h)
    c1b = (logF(h/2, d12, d34) - logF(-h/2, d12, d34))/h
    c3 = ((logF(2*h, d12, d34) - logF(-2*h, d12, d34)) - 2*(logF(h, d12, d34) - logF(-h, d12, d34)))/(2*(8 - 2)*h**3) * 2
    ev = (logF(h, d12, d34) + logF(-h, d12, d34) - 2*logF(h/1e6, d12, d34))/(2*h**2)
    pr('(2/3) d12=%d d34=%d : odd O(1/m) coeff of log F = %.3e (h/2: %.3e), even O(1/m^2) coeff ~ %.3e' %
       (d12, d34, float(c1), float(c1b), float(ev)))
open(os.path.join(HERE, 'cache', 'weyl.txt'), 'w').write('\n'.join(out) + '\n')
