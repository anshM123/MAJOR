"""
lrtight/lrt2_m1.py (LRINF_TIGHT.md §2.2): certified LOWER bounds of M1(beta) on beta-cells covering [0.1, 24], and the analytic
lower bound for beta >= 12 (used beyond 24).

  M1(beta) = int_0^inf 2 nu P_{-1/2+i nu}(cosh beta)/(e^{2 pi nu}+1) dnu = (sqrt2/pi) I(beta),  I = int_0^beta K(t) c(t) dt,
  K = (cosh beta - cosh t)^{-1/2},  c(t) = 2 int_0^inf nu cos(nu t)/(e^{2 pi nu}+1) dnu = cosh(t/2)/(4 sinh^2(t/2)) - 1/t^2.
  (c from int_0^inf sin(a x)/(e^{2 pi x}+1) dx = 1/(2a) - 1/(4 sinh(a/2)); Fubini is absolute.)
  Method = h13_cert_m1 region (b) (code/exactskel/h13_cert_m1.py, re-implemented here): t = beta(1 - v^2),
     I = int_0^1 g(v) dv,  g = 2 sqrt(beta) c(beta(1-v^2)) / sqrt( sinh(beta(1 - v^2/2)) S(beta v^2/2) ),  S(y) = sinh(y)/y,
  with beta a BALL covering the whole cell and v in nv equal cells: I(cell) >= sum_cells |cell| inf g(v-cell, beta-cell)  [Arb].
  c on t <= 2 by the series h(x)/x^4 = sum_{n>=2} d_n x^{2n-4} + tail (x = t/2), S by monotone endpoint bounds.
  Region beta >= 12 (h13 region (c), re-derived in LRINF_TIGHT §2.2 [P]):
     I/K0 >= 1/beta - 1/(2 sinh(beta/2)) - sqrt2 (e^{-beta/2}/4.7 + 4/beta^2) - 2 sqrt(2 ln 2)/(beta - ln 2)^2.
Output: cache/lrt2_m1.json  (list of [beta_a, beta_b, M1_lower]),  cache/lrt2_m1.txt
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial
from flint import arb, ctx
from lrt_common import A, ball, CACHE

ctx.prec = 128
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

NS = 16
DN = [None, None] + [Fr(2*n*(2*n - 1) - 2**(2*n - 1), factorial(2*n)) for n in range(2, NS + 1)]


def hx4_series(x):
    """h(x)/x^4, h = x^2 cosh x - sinh^2 x = sum_{n>=2} d_n x^{2n}; valid for 0 <= x <= 1 with the tail bound."""
    s = arb(0)
    for n in range(NS, 1, -1):
        s = s*x*x + A(DN[n])
    xm = arb(abs(x).upper())
    # |d_n| <= 2^{2n-1}/(2n)! for n >= 3 and the tail sum_{n>NS} 2^{2n-1} x^{2n-4}/(2n)! <= 2 * first term (x <= 1)
    tail = 2*arb(2)**(2*NS + 1)*xm**(2*NS - 2)/factorial(2*NS + 2)
    return s + arb(0, tail.upper())


def S_mono(y):
    """sinh(y)/y for a ball y >= 0 (increasing in y >= 0): enclosure from the endpoint values."""
    yl = arb(y.lower()); yh = arb(y.upper())
    if yl <= 0:
        lo_ = arb(1)
    else:
        lo_ = yl.sinh()/yl
    if yh <= 0:
        hi_ = arb(1)
    elif yh < arb(1)/10**6:
        hi_ = 1 + yh*yh/6*(1 + yh*yh)            # sinh(y)/y <= 1 + y^2/6 cosh(y) ... crude but safe for tiny y
    else:
        hi_ = yh.sinh()/yh
    return arb(lo_.lower()).union(arb(hi_.upper()))


def c_of(t):
    """c(t) for a ball t >= 0; None if the ball straddles the two regimes."""
    if t.upper() <= 2:
        x = t/2
        return hx4_series(x)/(4*S_mono(x)**2)
    if t.lower() > 1:
        h = t/2
        return h.cosh()/(4*h.sinh()**2) - 1/t**2
    return None


def g_inf(B, v):
    t = B*(1 - v*v)
    c = c_of(t)
    if c is None:
        return None
    den = ((B*(1 - v*v/2)).sinh()*S_mono(B*v*v/2))
    g = 2*B.sqrt()*c/den.sqrt()
    return g


def lower_I(ba, bb, nv):
    B = ball(ba, bb)
    tot = arb(0)
    stack = [(Fr(j, nv), Fr(j + 1, nv)) for j in range(nv)]
    while stack:
        a, b = stack.pop()
        g = g_inf(B, ball(a, b))
        if g is None or not g.is_finite():
            m = (a + b)/2
            stack += [(a, m), (m, b)]
            continue
        tot += arb(g.lower())*A(b - a)
    return arb(tot.lower())


def M1_lower_region_c(beta):
    """analytic lower bound (beta >= 12), evaluated in Arb at a point/ball beta (each subtracted term decreasing, see doc)."""
    b = beta
    K0 = 1/((b.cosh() - 1).sqrt())
    T0 = A(Fr(47, 10))
    ln2 = arb(2).log()
    I_over_K0 = 1/b - 1/(2*(b/2).sinh()) - arb(2).sqrt()*((-b/2).exp()/T0 + 4/b**2) - 2*(2*ln2).sqrt()/(b - ln2)**2
    return arb(2).sqrt()/arb.pi()*K0*I_over_K0


def cells(spec):
    out = []
    for a, b, h in spec:
        x = Fr(a)
        while x < Fr(b) - Fr(1, 10**12):
            out.append((x, min(x + Fr(h), Fr(b)))); x += Fr(h)
    return out


SPEC = [('0.10', '0.30', '0.0025'), ('0.30', '1.50', '0.005'), ('1.50', '3', '0.01'), ('3', '6', '0.025'), ('6', '12', '0.05'), ('12', '24', '0.1')]


def nv_of(b):
    return 1600 if b <= 1.5 else (3200 if b <= 3 else (6400 if b <= 6 else (12800 if b <= 12 else 25600)))


def main():
    pr(__doc__.split('Output')[0])
    res = []
    t0 = time.time()
    for (a, b) in cells(SPEC):
        I = lower_I(a, b, nv_of(float(b)))
        M = arb(2).sqrt()/arb.pi()*I
        Ml = float(arb(M.lower()).mid()) if M.lower() > 0 else float(arb(M.lower()).mid())
        res.append([str(a), str(b), Ml])
        if (a*10).denominator == 1 or b == Fr(24):
            pr(f'  beta in [{float(a):.3f}, {float(b):.3f}]: M1 >= {Ml:.6e}   ({time.time()-t0:.0f}s)')
    ok = all(r[2] > 0 for r in res)
    pr(f'  {len(res)} cells on [0.1, 24]; all lower bounds > 0: {ok}')
    # region c at beta >= 12: values at a few points (the bound function itself is used by lrt3 for beta >= 12)
    for bb in (12, 16, 24, 32, 48, 64, 128):
        pr(f'  region (c) lower bound at beta = {bb}: M1 >= {M1_lower_region_c(arb(bb)).str(8)}')
    json.dump(res, open(os.path.join(CACHE, 'lrt2_m1.json'), 'w'))
    open(os.path.join(CACHE, 'lrt2_m1.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
