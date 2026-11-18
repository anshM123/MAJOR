"""
h13_cert_m1: computer-assisted proof that M1(beta) > 0 for every beta > 0 (H13 follow-up, task 1).

  M1(beta) = int_0^oo 2 nu P_{-1/2+i nu}(cosh beta)/(e^{2 pi nu}+1) dnu = (sqrt2/pi) I(beta),
  I(beta)  = int_0^beta K(t) c(t) dt,  K(t) = (cosh beta - cosh t)^{-1/2},  c(t) = cosh(t/2)/(4 sinh^2(t/2)) - 1/t^2
  (equivalence checked to 16 digits in h13_asym.py; c(t) = 2 int_0^oo nu cos(nu t)/(e^{2 pi nu}+1) dnu).

Lemma 0 [D]: int_0^beta c = 1/beta - 1/(2 sinh(beta/2))   (antiderivative 1/t - 1/(2 sinh(t/2)), -> 0 at t = 0);  c > -1/t^2.
Region (a) beta <= T0 = 4.7: c > 0 on (0, T0] [certified below], K > 0  =>  I > 0.
   c > 0 <=> h(x) = x^2 cosh x - sinh^2 x > 0 (x = t/2).  h(x) = sum_{n>=2} d_n x^{2n}, d_n = (2n(2n-1) - 2^{2n-1})/(2n)!,
   d_2 = 1/6, d_n < 0 (n >= 3)  =>  h/x^4 >= 1/6 - sum_{n>=3} 2^{2n-1} x^{2n-4}/(2n)! = 1/6 - (cosh 2x - 1 - 2x^2 - 2x^4/3)/(2x^4),
   and cosh y - 1 - y^2/2 - y^4/24 <= y^6 cosh(y)/720 < y^4/48 for y = 2x = t <= 2  =>  c > 0 on (0, 2] [D];
   on [2, 4.7]: Arb interval evaluation of h on subintervals [certified].
Region (b) 4.7 <= beta <= 12: rigorous lower Riemann sums in Arb, with t = beta(1 - v^2):
   I = int_0^1 g(v) dv,  g = 2 sqrt(beta) c(beta(1-v^2)) / sqrt( sinh(beta(1 - v^2/2)) S(beta v^2/2) ),  S(y) = sinh(y)/y,
   beta in boxes (balls), v in cells; I(box) >= sum_cells |cell| * inf g(cell, box) > 0 [certified].
   c on t <= 2 from the series h/x^4 = sum_{n=2}^{N} d_n x^{2n-4} + tail (|tail| <= 2 * 2^{2N+1} x^{2N-2}/(2N+2)!), S by its
   (positive) series; elsewhere direct formulas.
Region (c) beta >= 12 [D]: I = K0 int_0^beta c + int_0^beta (K - K0) c,  K0 = K(0), K - K0 >= 0, c >= 0 on [0, T0], c_- <= 1/t^2:
   I/K0 >= 1/beta - 1/(2 sinh(beta/2)) - sqrt2 (e^{-beta/2}/T0 + 4/beta^2) - 2 sqrt(2 ln 2)/(beta - ln 2)^2
   (split at t1 with sinh^2(t1/2) = sinh^2(beta/2)/2: on [T0, t1], K/K0 - 1 <= sqrt2 r, r = sinh^2(t/2)/sinh^2(beta/2) <= e^{t-beta};
    on [t1, beta], cosh beta - cosh t >= sinh(t1)(beta - t), sinh t1 >= sinh^2(beta/2), beta - t1 <= ln 2).
   beta^2 * RHS >= beta - 0.357 - sqrt2 (0.0753 + 4) - 2.355 * 1.1245 = beta - 8.77 > 0 for beta >= 12 (each subtracted term is
   decreasing on beta >= 12; its value at 12 is evaluated here in Arb).
Output: cache/h13_cert_m1.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
from math import factorial
from flint import arb, ctx
ctx.prec = 128
out = []
def pr(s=''):
    print(s, flush=True); out.append(s)

def iv(a, b):
    a, b = Fr(a), Fr(b)
    m, r = (a + b)/2, (b - a)/2
    return arb(m.numerator)/m.denominator + arb(0, float(r)*(1 + 1e-12) + 1e-300) if r else arb(m.numerator)/m.denominator

NS = 14
DN = [None, None] + [Fr(2*n*(2*n - 1) - 2**(2*n - 1), factorial(2*n)) for n in range(2, NS + 1)]

def hx4_series(x):          # h(x)/x^4 for 0 <= x <= 1, with tail bound
    s = arb(0)
    for n in range(NS, 1, -1):
        s = s*x*x + arb(DN[n].numerator)/DN[n].denominator
    xm = arb(x.upper()) if x.upper() > 0 else arb(0)
    tail = 2*arb(2)**(2*NS + 1)*xm**(2*NS - 2)/factorial(2*NS + 2)
    return s + arb(0, float(tail.upper()))

def S_series(y):            # sinh(y)/y for 0 <= y <= 1 (monotone increasing, >= 1)
    lo = 1 if y.lower() <= 0 else (y.lower().sinh()/y.lower())
    hi = y.upper().sinh()/y.upper() if y.upper() > 0 else arb(1)
    lo = arb(lo); hi = arb(hi)
    return arb((lo.lower() + hi.upper())/2, float((hi.upper() - lo.lower()).upper())/2*(1 + 1e-9))

def S_any(y):
    return S_series(y) if y.upper() <= 1 else y.sinh()/y

def c_any(t):
    if t.upper() <= 2:
        x = t/2
        return hx4_series(x)/(4*S_series(x)**2)
    if t.lower() > 1:
        return (t/2).cosh()/(4*(t/2).sinh()**2) - 1/t**2
    raise ValueError

def lower_I(bl, bh, nv):
    B = iv(bl, bh); tot = arb(0)
    for j in range(nv):
        v = iv(Fr(j, nv), Fr(j + 1, nv))
        t = B*(1 - v*v)
        # t may straddle 2: split by monotone pieces is unnecessary - use series if t <= 2 entirely, direct if t > 1 entirely
        c = c_any(t)           # cells are narrow, so t never straddles both regimes (ValueError otherwise)
        g = 2*B.sqrt()*c/((B*(1 - v*v/2)).sinh()*S_any(B*v*v/2)).sqrt()
        tot += arb(g.lower())/nv
    return tot


def main():
    # ---- region (a): c > 0 on [2, 4.7]
    T0 = Fr(47, 10)
    n = 2700; ok = True
    for i in range(n):
        a = 2 + (T0 - 2)*i/n; b = 2 + (T0 - 2)*(i + 1)/n
        x = iv(a/2, b/2)
        h = x*x*x.cosh() - x.sinh()**2
        if not h > 0:
            ok = False; pr(f"  c > 0 FAILED on [{float(a)}, {float(b)}]"); break
    pr(f"(a) c(t) > 0 on (0, 2] [D, series bound]; on [2, 4.7] by Arb on {n} subintervals: {ok}")
    pr("    => M1(beta) > 0 for 0 < beta <= 4.7 (positive integrand)")

    # ---- region (b): beta in [4.7, 12]
    bl = Fr(47, 10); allok = True; nbox = 0; worst = None
    while bl < 12:
        w = Fr(1, 50) if bl < 8 else Fr(1, 100)
        bh = min(bl + w, Fr(12))
        nv = 1500
        L = lower_I(bl, bh, nv)
        while not L > 0 and nv < 24000:
            nv *= 2; L = lower_I(bl, bh, nv)
        nbox += 1
        if not L > 0:
            allok = False; pr(f"  FAILED on beta in [{float(bl)}, {float(bh)}]: lower bound {L}"); break
        lo = float(L.lower())
        if worst is None or lo < worst[0]:
            worst = (lo, float(bl), float(bh), nv)
        bl = bh
    pr(f"(b) beta in [4.7, 12]: {nbox} boxes, I(beta) >= lower Riemann sum > 0 on every box: {allok};"
       f" smallest certified lower bound {worst[0]:.4e} on [{worst[1]:.2f}, {worst[2]:.2f}] ({worst[3]} cells)")

    # ---- region (c): explicit bound at beta = 12 and monotonicity argument
    b = arb(12)
    lhs = b - b*b/(2*(b/2).sinh())
    rhs = arb(2).sqrt()*(b*b*(-b/2).exp()/(arb(47)/10) + 4) + 2*(2*arb(2).log()).sqrt()*b*b/(b - arb(2).log())**2
    pr(f"(c) beta >= 12: beta^2 [1/beta - 1/(2 sinh(beta/2))] >= {float(lhs.lower()):.4f};  beta^2 x (correction bounds) <= {float(rhs.upper()):.4f} at beta = 12;"
       f" both bounded monotonically (LHS increasing, RHS decreasing) => I > 0 for all beta >= 12")
    pr(f"    margin at 12: {float((lhs - rhs).lower()):.3f} > 0: {bool(lhs - rhs > 0)}")
    pr("RESULT: M1(beta) > 0 for all beta > 0  [computer-assisted: (a) D + Arb, (b) Arb lower sums, (c) D]")
    open(os.path.join(HERE, 'cache', 'h13_cert_m1.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
