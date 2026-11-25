"""
contour/cb_w.py (CONTOUR_BOUNDS.md §4): the local density W(y) = -2 B(iy) Im f(iy) on 0 <= y <= 3/4, certified.

Closed form (derived in CONTOUR_BOUNDS.md §4; cross-checked here against -2 B(iy) Im f(iy) evaluated directly with cb_common):
  on nu = iy the Gaussian part and Q R_- are real; only the Log of chi(nu + i) (cut, boundary value from Re nu > 0) gives Im f:
     Im f(iy) = -pi y (2+y) theta(is) Y(is) R_+(is),   s = 1 + y,   Y(is) = (y0 - y1 s^2 + y2 s^4) e^{s^2/2};
  B(iy) = Gamma(-2y)/Gamma(1/2 - y)^2, and with the duplication formula
     W(y) = (1/2 - y)^2 V(y),
     V(y) = 4 sqrt(pi) (2+y) 2^{-2y} Gamma(1-y)/Gamma(3/2-y) * theta(is) * Y(is) * (2s-11)(2s-5)(2s-1)^2(2s+1) / (s(2s-7)(2s+3)(2s+7)),
  analytic on [0, 3/4] (the pole of B at y = 0 and its zero at y = 1/2, and the zero of R_+ at y = 1/2, are cancelled exactly).
Certificate: V is enclosed on NC cells of [0, 3/4] by interval evaluation (ball arithmetic on the whole cell); V > 0 on every cell
  => W >= 0 on [0, 3/4], with W = 0 only at y = 1/2 (double zero).
Outputs: cache/cb_w.json (cell enclosures of W), cache/cb_w.txt.
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cb_common as K
from flint import arb, acb

NC = 7500


def V_enc(W, y):
    """enclosure of V on the arb ball y."""
    s = 1 + y
    th = 1 - ((K.I*s + W.nuc)/W.acut).erfc()/2 - ((W.nuc - K.I*s)/W.acut).erfc()/2
    th = th.real                                         # theta(is) is real (theta even with real Taylor coefficients)
    Yv = (W.y[0] - W.y[1]*s**2 + W.y[2]*s**4)*(s*s/2).exp()
    rat = (2*s - 11)*(2*s - 5)*(2*s - 1)**2*(2*s + 1)/(s*(2*s - 7)*(2*s + 3)*(2*s + 7))
    pre = 4*arb.pi().sqrt()*(2 + y)*(-2*y*arb(2).log()).exp()*(1 - y).gamma()/(arb(3)/2 - y).gamma()
    return pre*th*Yv*rat


def W_direct(W, y):
    """-2 B(iy) Im f(iy) from the general evaluation (check only; y not 0 or 1/2)."""
    z = acb(0, y)
    return -2*(K.B_val(z)*K.f_point(W, z).imag).real


def main():
    t0 = time.time()
    lines = []; pr = K.C.pr_to(lines)
    pr(__doc__.split('Outputs')[0])
    W = K.load()
    # checks of the closed form against the direct evaluation
    for y in ('0.05', '0.2', '0.37', '0.49', '0.51', '0.6', '0.74'):
        ya = arb(y)
        wd = W_direct(W, ya); wc = (arb('0.5') - ya)**2*V_enc(W, ya)
        pr(f'  check y = {y}: direct {wd.str(12)}   closed form {wc.str(12)}')
        assert (wd - wc).contains(0) or abs(wd - wc) < arb('1e-25')
    cells = []
    h = arb(3)/4/NC
    Wmin = {}; Wmax = {}
    allpos = True
    for k in range(NC):
        a = h*k; b = h*(k + 1)
        y = (a + b)/2 + arb(0, 1)*0
        y = arb((a + b)/2, (b - a)/2)
        V = V_enc(W, y)
        q = (arb('0.5') - y)
        q2 = q*q
        if q.contains(0):
            q2lo = arb(0)
        else:
            q2lo = arb(q2.lower())
        Wlo = q2lo*arb(V.lower()); Whi = arb(q2.upper())*arb(V.upper())
        if not V > 0:
            allpos = False
        cells.append((float(a.mid()), float(b.mid()), float(Wlo.lower()), float(Whi.upper())))
    assert allpos
    pr(f'  V > 0 on all {NC} cells of [0, 3/4]  =>  W >= 0 on [0, 3/4]  [certified]')

    def mn(lo, hi):
        return min(c[2] for c in cells if c[0] >= lo - 1e-15 and c[1] <= hi + 1e-15)

    def mx(lo, hi):
        return max(c[3] for c in cells if c[0] >= lo - 1e-15 and c[1] <= hi + 1e-15)
    res = dict(min_0_001=mn(0, 0.01), max_0_001=mx(0, 0.01), max_0_18=mx(0, 0.125), max_18_34=mx(0.125, 0.75),
               min_0_18=mn(0, 0.125), W0=(2*(-W.y[0] + W.y[1] - W.y[2])*arb('0.5').exp()*arb(18)/25).str(15))
    pr(f'  W(0) = -2 Y(i) R_+(i) = {res["W0"]}')
    pr(f'  min W on [0, 0.01] >= {res["min_0_001"]:.6f};  max W on [0, 0.01] <= {res["max_0_001"]:.6f}')
    pr(f'  min W on [0, 1/8] >= {res["min_0_18"]:.6f};  max W on [0, 1/8] <= {res["max_0_18"]:.6f}')
    pr(f'  max W on [1/8, 3/4] <= {res["max_18_34"]:.6f}')
    for yy in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75):
        ya = arb(yy)
        pr(f'    W({yy}) = {((arb("0.5") - ya)**2*V_enc(W, ya)).str(10)}')
    res['cells'] = cells
    json.dump(res, open(os.path.join(K.CACHE, 'cb_w.json'), 'w'))
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(K.CACHE, 'cb_w.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
