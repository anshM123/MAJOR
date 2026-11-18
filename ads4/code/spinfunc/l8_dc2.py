"""
l8_dc2 (M6 round 4, doc 22.3): graviton-specific O(1/m^2) correction of C^_2 per helicity, and the O(1/m) helicity mixing, on the
(Delta, x) grid, with the nu-only normalisation removed EXACTLY by the heavy reference state (l8_analysis):
  D_h(x; Delta)  = L_h(nu; Delta) - L_h(nu; Delta_ref),  nu = x m(Delta)          [log graviton/scalar/W_h, cumulant form]
  Re prescription (doc 22.1): the functional is Re[C^CLS/N^]; per point Re e^{D} -> Re D + log cos(Im D).
  dplus(x)  = m^2 [Re D_++(x) - Re D_++(0)]                   = dc2^{spin,++}(x) - dc2^{spin,++}(0)        (pipeline units)
  s(x)      = m^2 [l_s(nu;Delta) - l_s(nu;Delta_ref) - log(env2(x) P_J(x))] - (same at 0)   (scalar crystal: m^2 rho_2/P_J, rel. x->0)
  Eminus(x) = m^2 [Re D_+-(x) - Re D_+-(0)] + s(x)  = graviton +- correction rel. to its own LO density (incl. scalar curvature), rel. x->0
              -> dc2^{spin,+-}(x) - dc2^{spin,+-}(0) = Eminus(x) - m^2 [rho_2(x)/d44~(x) - rho_2(0)]   (d44~ := d44, pipeline model)
  offp(x)   = m Re[M~_{+-}/M~_{++}] sqrt(c_++/c_+-) P_J(x)  : off-diagonal in pipeline units = (2-y) offp/m
              (c_h per-state constants from x -> 0; offp_im: same with Im, the part that drops under the Re prescription)
x -> 0 values: quadratic fit in x through x = 0.05, 0.1, 0.2 (checked against 0.1, 0.2 linear in x^2).
Usage: python3 l8_dc2.py [J ...] -> cache/l8_dc2.pkl, cache/l8_dc2.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import eval_legendre
import l8_analysis as L
import l8_points as P
import l7_dc2_analysis as An

DREF = 20000.0
XS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
DS = (100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0, 6400.0)


def env2(x):
    return (2 - x*x)/(1 - x*x)**2


def x0val(xs, v, w=None):
    """x -> 0 value c of a relative correction v(x), from the ABSOLUTE form w v = c w + b x^2 + d x^4 on the 3 smallest x
    (w = LO density shape; regular through zeros of w).  Error estimate: same with 2 points and no x^4 term."""
    xs = np.asarray(xs); v = np.asarray(v); w = np.ones_like(xs) if w is None else np.asarray(w)
    A = np.array([w, xs**2, xs**4]).T
    a = np.linalg.solve(A[:3], (w*v)[:3])[0]
    b = np.linalg.solve(A[:2, :2], (w*v)[:2])[0]
    return a, abs(a - b)


def analyse(J, pr=print):
    R = L.records(J)
    Bref = L.basis(J, DREF, R)
    hels = L.HEL if J >= 4 else ('++',)
    out = {}
    for D in DS:
        m = P.mphys(D, J)
        if not all((D, round(x*m, 9)) in R and (DREF, round(x*m, 9)) in R for x in XS):
            continue
        B = L.basis(J, D, R)
        rows = [(x, L.point(J, D, x*m, R, B), L.point(J, DREF, x*m, R, Bref)) for x in XS]
        rec = dict(m=m, sv3=B['sv3'], bimag=B['imag'])
        PJ = eval_legendre(J, 1 - 2*np.array(XS)**2)
        sh = lambda x: np.log(env2(x)*abs(eval_legendre(J, 1 - 2*x*x)))
        sraw = np.array([m*m*((p['ls'] - q['ls']).real - sh(x) + sh(q['x'])) for x, p, q in rows])
        s0, se = x0val(XS, sraw, PJ); rec['s'] = sraw - s0; rec['s_x0err'] = se
        for h in hels:
            Dh = np.array([p[h] - q[h] for x, p, q in rows])
            ph = np.angle(np.exp(1j*(Dh.imag - np.round(Dh.imag[0]/np.pi)*np.pi)))       # wrap; per-state sign is real
            re = Dh.real + np.log(np.abs(np.cos(ph)))
            r0, e0 = x0val(XS, re, PJ if h == '++' else An.d44(J, 1 - 2*np.array(XS)**2))
            rec['d' + h] = m*m*(re - r0); rec['x0err' + h] = m*m*e0
            rec['mIm' + h] = m*ph
        if J >= 4:
            rec['Eminus'] = rec['d+-'] + rec['s']
            lc = np.array([(p['+-'] - p['++']).real for x, p, q in rows]) - np.array([(q['+-'] - q['++']).real for x, p, q in rows])
            ratio0, _ = x0val(XS, lc, An.d44(J, 1 - 2*np.array(XS)**2))                                           # log(c_+-/c_++) (Delta vs ref)
            cref = np.mean([(q['+-'] - q['++']).real for x, p, q in rows])    # same at the reference (x_ref ~ 0)
            lcc = ratio0 + cref                                                 # log(c_+-/c_++) at Delta
            r01 = np.array([p['r01'] for x, p, q in rows])
            rec['offp'] = m*(r01*np.exp(-0.5*lcc)).real*PJ    # = m Re[M~+-/M~++] sqrt(c++/c+-) P_J
            rec['offp_im'] = m*(r01*np.exp(-0.5*lcc)).imag*PJ
            rec['asym'] = max(m*m*p['asym'] for x, p, q in rows)
        out[D] = rec
    return out


def report(J, out, pr):
    if not out:
        pr(f"J={J}: no complete rows"); return
    Ds = sorted(out)
    pr(f"J={J}: Delta rows {Ds}; x = {XS}")
    for key in (['d++', 'mIm++', 's'] + (['d+-', 'Eminus', 'mIm+-', 'offp', 'offp_im'] if J >= 4 else [])):
        pr(f"  {key}:")
        for D in Ds:
            pr(f"    Delta={D:6.0f} m={out[D]['m']:7.1f}: " + " ".join(f"{v:+8.1f}" for v in out[D][key]))
    pr("  x->0 extrapolation uncertainty (c2 units) max over Delta: " +
       ", ".join(f"{k} {max(out[D][k] for D in Ds):.1f}" for k in ('x0err++', 'x0err+-', 's_x0err') if k in out[Ds[0]]))
    if J >= 4:
        pr("  basis: sv3 " + " ".join(f"{out[D]['sv3']:.0e}" for D in Ds) + "; |asym| m^2 max " + f"{max(out[D]['asym'] for D in Ds):.1f}")


if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1:]] or [0, 2, 4, 6, 8, 10, 12]
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    fn = os.path.join(HERE, 'cache', 'l8_dc2.pkl')
    store = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    for J in Js:
        out = analyse(J, pr); report(J, out, pr); store[J] = out
    pickle.dump(store, open(fn, 'wb'))
    open(os.path.join(HERE, 'cache', 'l8_dc2_J' + '_'.join(map(str, Js)) + '.txt'), 'w').write("\n".join(lines) + "\n")


def summary(pr=print, Js=None):
    """asymptotic (a0 + a1/m fit over Delta = 400...3200) graviton-specific corrections on x in [0.3, 0.7], ABSOLUTE form
    (A+ = P_J dc2_++, A- = d44 Eminus: c2 units times the LO density shape normalised to 1 at x = 0), and the mixing (Re vs Im)."""
    store = pickle.load(open(os.path.join(HERE, 'cache', 'l8_dc2.pkl'), 'rb'))
    xs = np.array(XS); win = (xs >= 0.3) & (xs <= 0.7)
    FD = (400.0, 800.0, 1600.0, 3200.0)
    pr("# l8 summary (grid points with |P_J| < 0.05, |d44| < 0.02 excluded): asymptotic a0(x) [and value at m ~ 100 / 400] of |A_h| on x in [0.3, 0.7]; mixing max m|Re|, m|Im| (pipeline units)")
    for J in sorted(store) if Js is None else Js:
        out = store[J]
        if not set(FD) <= set(out):
            continue
        ms = np.array([out[D]['m'] for D in FD]); A = np.array([np.ones_like(ms), 1/ms]).T
        PJ = eval_legendre(J, 1 - 2*xs**2)
        keys = [('A++', lambda r: PJ*(r['m']**2 + r['s'])*np.expm1(r['d++']/r['m']**2))]
        if J >= 4:
            d44 = An.d44(J, 1 - 2*xs**2)
            keys.append(('A+-', lambda r: d44*((r['m']**2 + r['s'])*np.exp(r['d+-']/r['m']**2) - r['m']**2)))
        s = f"  J={J:2d}:"
        okp = np.abs(PJ) >= 0.05; okm = okp & (np.abs(An.d44(J, 1 - 2*xs**2)) >= 0.02) if J >= 4 else okp
        for (nm, f), ok in zip(keys, (okp, okm)):
            w = win & ok
            a = np.linalg.lstsq(A, np.array([f(out[D]) for D in FD]), rcond=None)[0]
            s += f"  {nm}: a0 max {np.abs(a[0][w]).max():6.0f} (a1/m at m=400: {np.abs(a[1][w]).max()/400:5.0f})"
            if 100.0 in out:
                s += f", m=100 direct {np.abs(f(out[100.0])[w]).max():6.0f}"
        if J >= 4:
            s += (f"  | mix m|Re| max: D=100 {np.abs(out[100.0]['offp']).max():.1f}, D=400 {np.abs(out[400.0]['offp']).max():.1f}, "
                  f"D=3200 {np.abs(out[3200.0]['offp']).max():.1f};  m|Im| max (all D>=400) {np.abs(out[3200.0]['offp_im']).max():.1f}")
        pr(s)
