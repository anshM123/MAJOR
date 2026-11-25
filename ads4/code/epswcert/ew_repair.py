"""
epswcert/ew_repair.py (EPSW_CERT.md §3): minimal repairs of the published k = 2 continuum weight (prl/data/k2, read-only).

r(0) = E(0) = 2 Re[Q(-i) R_-(-i) B(-i)] (the log-basis terms vanish at mu = 0: l(nu) = 0 at nu = -i), Q(-i) R_-(-i) = -16/9 * 18/25
= -32/25, and B(-i) depends on x_0 (the c = 0 Gaussian) through g_0(-i) = exp(1/(2 sigma_0^2)).  So
    dr(0)/dx_0 = -(64/25) exp(1/(2 sigma_0^2))  (exact).
r'(0+) = E'(0+) - psi'(0) depends on y_0 (the first log-basis coefficient) through the kink of l:  d r'(0+)/d y_0 = -2 pi (18/25) e^(1/2)
(computed below in ball arithmetic from the F1 series, not assumed), and r(0) does not depend on y.

Repair 'rep'  : x_0 -> the IEEE double nearest x_0* (r(0) = 0), moved by ulps until r(0) > 0 is CERTIFIED (sign that helps T0).
Repair 'repk' : 'rep' plus y_0 -> the IEEE double nearest the value making r'(0+) = 0 (kink matched; optional, see EPSW_CERT §3).
Output: cache/k2_continuum_D_<Dg>_{rep,repk}.json (same format as prl/data/k2; only x_re[0] (and y_re[0]) changed; eps fields
removed and replaced by a pointer to the certificate), log cache/ew_repair.txt.
"""
import os, sys, json, math, copy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, acb
import ew_common as C


def r0_r1(W, N=8):
    F1, F2 = C.F12_series(W, N)
    return 2*F1[0].real, 2*F1[1].real - W.kappa


def nextafter(x, d):
    return math.nextafter(x, d)


def main():
    lines = []; pr = C.pr_to(lines)
    pr(__doc__.split('Output')[0])
    for Dg in (300, 1000):
        W = C.Weight(Dg)
        d = W.json
        assert d['centres_sigma'][0] == [0.0, 0.45]
        x0 = float(d['x_re'][0]); y0 = float(d['y_re'][0])
        r0, r1 = r0_r1(W)
        D = -arb(64)/25*(1/(2*arb(0.45)**2)).exp()
        pr(f'Dg = {Dg}: published x_0 = {x0!r}, y_0 = {y0!r};  r(0) = {r0.str(12)},  r\'(0+) = {r1.str(12)};  dr(0)/dx_0 = {D.str(15)}')
        # --- x0 repair
        xs = x0 - r0/D
        xn = float(xs.mid())
        def with_x(xv):
            d2 = copy.deepcopy(d); d2['x_re'][0] = repr(xv)
            p = os.path.join(C.CACHE, f'_tmp_{Dg}.json'); json.dump(d2, open(p, 'w'))
            W2 = C.Weight(Dg, p); os.remove(p)
            return d2, W2
        for it in range(200):
            d2, W2 = with_x(xn)
            r0n, r1n = r0_r1(W2)
            if r0n > 0:
                break
            xn = nextafter(xn, -math.inf)          # D < 0: decreasing x_0 increases r(0)
        assert r0n > 0
        dx = xn - x0
        pr(f'   rep : x_0 = {xn!r} (delta x_0 = {dx:.10e}, {abs(dx)/abs(x0):.2e} relative); r(0) = {r0n.str(6)} > 0 [certified]; '
           f'r\'(0+) = {r1n.str(12)}')
        d2['name'] = d['name'] + ' -- REPAIRED (EPSW_CERT.md): x_re[0] changed so that r(0) = E(0) - psi(0) > 0 (certified)'
        d2['repair'] = dict(changed='x_re[0]', published=d['x_re'][0], new=repr(xn), delta=dx, r0_new=r0n.str(10),
                            reason='published weight has r(0) = %s < 0: realised LR image ~ r(0)/(2 sinh(beta/2)) beats the floor '
                                   '~ e^{-beta/2}/beta at beta ~ 2.5e4 (EPSW_CERT.md)' % r0.str(8))
        for k in ('eps_W_re', 'eps_W_im_null', 'certificate', 'checks'):
            d2.pop(k, None)
        d2['certificate'] = 'see ads4/EPSW_CERT.md and code/epswcert/cache/ew_cert_%d_re_rep.txt (ball arithmetic)' % Dg
        json.dump(d2, open(os.path.join(C.CACHE, f'k2_continuum_D_{Dg}_rep.json'), 'w'), indent=0)
        # --- kink repair (on top of rep)
        F1a, _ = C.F12_series(W2, 4)
        W2.y[0] = W2.y[0] + 1
        F1b, _ = C.F12_series(W2, 4)
        W2.y[0] = W2.y[0] - 1
        dk = 2*F1b[1].real - 2*F1a[1].real                      # d r'(0+)/d y_0 (exact: linear)
        pr(f'        d r\'(0+)/d y_0 = {dk.str(15)}  (formula -2 pi (18/25) e^(1/2) = {(-2*arb.pi()*arb(18)/25*arb(0.5).exp()).str(15)})')
        yn = float((W2.y[0] - r1n/dk).mid())
        d3 = copy.deepcopy(d2); d3['y_re'][0] = repr(yn)
        p = os.path.join(C.CACHE, f'k2_continuum_D_{Dg}_repk.json')
        d3['name'] = d['name'] + ' -- REPAIRED+KINK (EPSW_CERT.md): x_re[0], y_re[0] changed (r(0) > 0, r\'(0+) ~ 0)'
        d3['repair'] = dict(d2['repair'], changed='x_re[0], y_re[0]', y0_published=d['y_re'][0], y0_new=repr(yn))
        d3['certificate'] = 'see ads4/EPSW_CERT.md and code/epswcert/cache/ew_cert_%d_re_repk.txt' % Dg
        json.dump(d3, open(p, 'w'), indent=0)
        W3 = C.Weight(Dg, p)
        r0k, r1k = r0_r1(W3)
        pr(f'   repk: y_0 = {yn!r} (delta y_0 = {yn - y0:.6e}); r(0) = {r0k.str(6)}, r\'(0+) = {r1k.str(6)}')
    open(os.path.join(C.CACHE, 'ew_repair.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
