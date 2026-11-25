"""
contour/cb_lines.py (CONTOUR_BOUNDS.md §§1-3, 5): certified sup / integral bounds of |f|, |B| and |f B| on the contour pieces of §38.
  short : nu = t + 3i/4, 0 <= t <= 1/4                         (5 panels, 5 subcells each)
  vert  : nu = 1/4 + iy, 3/4 <= y <= 9/4                         (75 panels of width 0.02, 5 subcells each: 375 cells of width 0.004)
  horiz Y : nu = t + iY, t >= 1/4, Y in {1, 1.25, 1.5, 1.75, 2, 2.25}: panels to T_END = c_max + 25, analytic tail beyond.
Method (cb_common.line_panel): Taylor polynomial of order N of f (and B) at the panel centre computed exactly in ball arithmetic
(power-series arithmetic: the cancellations among the RBF terms happen in the coefficients), truncation by Cauchy with a disc sup
(ring of boxes + explicit Gaussian bounds), far RBF terms bounded on the disc, subcells by exact polynomial Taylor shift,
sup|p| <= sum|c_n| hs^n, int|p| <= sum |c_n| 2 hs^(n+1)/(n+1).
Tail t >= T_END (horiz): for nu' = t + iy', |y'| <= Y + 1:
  |theta(nu')| <= e^{y'^2/a^2}   (theta(z) = (a sqrt pi)^{-1} int_{-nu_c}^{nu_c} e^{-((z+s)/a)^2} ds),
  |Q R_-(nu')| <= qbar(T_END) (product of decreasing ratios in |nu'| >= T_END),
  int_T^oo sum_k |x_k| |g_k| dt <= sum_k |x_k| e^{y'^2/(2 s_k^2)} s_k sqrt(pi/2) erfc((T - c_k)/(s_k sqrt 2)),
  Y-term <= e^{-T^2/4} (as in ew_cert.tail),
  |B(nu)| = 2^{-2Y-1} pi^{-1/2} |Gamma(z)/Gamma(z+1/2)|, z = i nu, and Gamma(z)/Gamma(z+1/2) = [(z+1/2)(z+3/2)(z+5/2)/(z(z+1)(z+2))]
     * Gamma(z+3)/Gamma(z+7/2), |Gamma(z+3)/Gamma(z+7/2)| = |B(z+3, 1/2)|/sqrt(pi) <= B(3 - Y, 1/2)/sqrt(pi)   (Re(z+3) = 3 - Y > 0).
Usage: python3 cb_lines.py short | vert | horiz Y      -> cache/cb_<piece>.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cb_common as K
from flint import arb, acb


def ser(rows):
    return [{k: (float(v.upper()) if isinstance(v, arb) else v) for k, v in r.items()} for r in rows]


def run_short(W, pr):
    rows = []
    N = 32; h = arb(1)/40; rho = 0.2
    for p in range(5):
        t0 = (2*p + 1)*h
        out = K.line_panel(W, acb(t0, arb(3)/4), acb(1), h, rho, N, 5)
        for r in out:
            r['t0'] = float((t0 + r['s0']).mid()); r['t1'] = float((t0 + r['s1']).mid())
        rows += out
    return rows


def run_vert(W, pr):
    rows = []
    N = 32; h = arb(1)/100; rho = 0.2
    for p in range(75):
        y0 = arb(3)/4 + (2*p + 1)*h
        out = K.line_panel(W, acb(arb(1)/4, y0), acb(0, 1), h, rho, N, 5)
        for r in out:
            r['y0'] = float((y0 + r['s0']).mid()); r['y1'] = float((y0 + r['s1']).mid())
        rows += out
    return rows


def tail(W, Y, T):
    """int_T^oo |f(t + iY)| dt and sup_{t >= T} |B(t + iY)| (upper bounds), see module doc."""
    Ya = arb(Y); ymax = Ya + 1
    th = (ymax*ymax/(W.acut*W.acut)).exp()
    R0 = arb(T)
    qbar = 32*R0*(2*R0 + 3)*(2*R0 + 5)*(2*R0 + 11)/((2*R0 - 1)*(2*R0 - 7)*(2*R0 - 3)*(2*R0 - 7))
    tot = arb(0)
    for (c, s, x), (ca, sa, xa, inv) in zip(W.sc, W.sc_arb):
        z = (R0 - ca)/(sa*arb(2).sqrt())
        assert z > 0
        tot += abs(xa)*(ymax*ymax*inv).exp()*sa*(arb.pi()/2).sqrt()*z.erfc()
    # Y-term: |Y l R_-| <= P(t) e^{-(t^2 - ymax^2)/2} with P polynomial-log; e^{ymax^2/2} P(t) <= e^{t^2/4} for t >= T (asserted at T,
    # and log P grows logarithmically), int_T^oo e^{-t^2/4} <= (2/T) e^{-T^2/4}
    ysum = sum((abs(y) for y in W.y), arb(0))
    m2 = R0 + ymax + 2
    Pb = ysum*m2**4*(m2*m2 + 1)*(2*m2.log() + arb.pi() + 1)*2*(2*m2 + 11)**6/(R0*(2*R0 - 7)**3)
    assert Pb*(ymax*ymax/2).exp() < (R0*R0/4).exp()
    ytail = 2/R0*(-(R0*R0)/4).exp()
    intf = 2*th*(qbar*tot + ytail)                     # two terms: G(t+iY) and conj G(t-iY), each with |y'| <= Y+1
    # |B|: z = i nu = -Y + i t, |z| >= T
    beta = (3 - Ya).gamma()*arb.pi().sqrt()/(arb(7)/2 - Ya).gamma()
    rat = (R0 + arb('0.5'))*(R0 + arb('1.5'))*(R0 + arb('2.5'))/((R0 - 3)*(R0 - 3)*(R0 - 3))   # |z+k| <= |z| + k, |z+j| >= |z| - 3
    Bsup = (-(2*Ya + 1)*arb(2).log()).exp()/arb.pi().sqrt()*rat*beta/arb.pi().sqrt()
    return intf, Bsup


def run_horiz(W, Y, pr):
    rows = []
    Ya = arb(Y)
    T_END = int(W.cmax + 25) + 1
    # region A: [1/4, 45], h = 0.02, rho = 0.2 (sigma = 0.45 terms; poles / branch points of f, B at distance >= 1/4 from the line)
    # region B: [45, T_END], h = 0.1, rho = 0.5 (sigma = 0.9 only)
    t = arb(1)/4
    edgesA = (arb(1)/4, arb(45))
    hA = arb(1)/50; NA = 40; rhoA = 0.2
    nA = int(float(((edgesA[1] - edgesA[0])/(2*hA)).upper())) + 1
    hA = (edgesA[1] - edgesA[0])/(2*nA)
    t0s = [(edgesA[0] + (2*p + 1)*hA, hA, rhoA, NA) for p in range(nA)]
    hB = arb(1)/10; NB = 40; rhoB = 0.5
    nB = int(float(((T_END - 45)/(2*hB)).upper())) + 1
    hB = (arb(T_END) - 45)/(2*nB)
    t0s += [(arb(45) + (2*p + 1)*hB, hB, rhoB, NB) for p in range(nB)]
    tt = time.time()
    for i, (t0, h, rho, N) in enumerate(t0s):
        out = K.line_panel(W, acb(t0, Ya), acb(1), h, rho, N, 4 if rho < 0.3 else 2)
        for r in out:
            r['t0'] = float((t0 + r['s0']).mid()); r['t1'] = float((t0 + r['s1']).mid())
        rows += out
        if i % 500 == 0:
            pr(f'    Y = {Y}: panel {i}/{len(t0s)} t = {float(t0.mid()):.2f}  int|f| so far {sum(r["intf"] for r in rows).upper()} '
               f'[{time.time() - tt:.0f}s]')
    intf, Bsup = tail(W, Y, T_END)
    return rows, dict(T_END=T_END, tail_intf=float(intf.upper()), tail_supB=float(Bsup.upper()),
                      tail_intfB=float((intf*Bsup).upper()))


def main():
    piece = sys.argv[1]
    t0 = time.time()
    lines = []; pr = K.C.pr_to(lines)
    W = K.load()
    extra = {}
    if piece == 'short':
        rows = run_short(W, pr); tag = 'short'
    elif piece == 'vert':
        rows = run_vert(W, pr); tag = 'vert'
    else:
        Y = float(sys.argv[2]); tag = f'horiz_{Y}'
        rows, extra = run_horiz(W, Y, pr)
    tot = {k: sum((r[k] for r in rows), arb(0)) for k in ('intf', 'intfB')}
    mx = {k: max(float(r[k].upper()) for r in rows) for k in ('supf', 'supB', 'supfB', 'ef', 'eB')}
    summ = dict(piece=tag, ncells=len(rows), int_f=float(tot['intf'].upper()), int_fB=float(tot['intfB'].upper()), **mx, **extra)
    if 'tail_intf' in extra:
        summ['int_f_total'] = float((tot['intf'] + arb(extra['tail_intf'])).upper())
        summ['int_fB_total'] = float((tot['intfB'] + arb(extra['tail_intfB'])).upper())
        summ['supB_total'] = max(mx['supB'], extra['tail_supB'])
    pr(f'RESULT {tag}: ' + ', '.join(f'{k} = {v:.6e}' if isinstance(v, float) else f'{k} = {v}' for k, v in summ.items()) +
       f'  [{time.time() - t0:.0f}s]')
    json.dump(dict(summary=summ, rows=ser(rows)), open(os.path.join(K.CACHE, f'cb_{tag}.json'), 'w'))
    open(os.path.join(K.CACHE, f'cb_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
