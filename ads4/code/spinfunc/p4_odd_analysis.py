"""
p4 (PARITY_ODD_ROUTEA.md sect. 5-6): analysis of the parity-odd (odd J, +-) route-A' crystal action from the p3 caches.

(a) LO check: order-0 Ĉ_2 action of the odd structure / scalar k = 2 crystal action / W_{+-}(J, x) = d^J_44/((1-x^2)^2 P_J) must be
    x-independent (stripped flat +- density); Ĉ_3/Ĉ_2 at order 0 must be (2i) x^2/(2-y) (sigma_{+-} = -1 relative to ++; H0 doc 23.3).
(b) graviton-specific corrections in the LINEAR (prefactor) form with the heavy reference (Delta_ref = 2e4, same nu, same kind) removed
    order by order (= l9_analysis.grid for the 1x1 +- element), per-state constant a_{+-} from Ĉ_2 at x <= 0.25, nu >= 15, CHMRS de-shift
    (H0 doc 23.2; L9_DESHIFT semantics: P4_DESHIFT=1 default).  Output Gamma_{k,+-}(x; m) in c2 units x LO shape (l9 conventions).
Functions env, Dc, shift, pred, cum are copies of l9_analysis (so that a concurrent edit of l9_* cannot change these numbers).
Usage: python3 p4_odd_analysis.py [J ...] -> cache/p4_odd_analysis.txt, cache/p4_grid.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import legendre as LG
from scipy.special import eval_jacobi, eval_legendre

DREF = 20000.0
XFIT = 0.25
NUFIT = float(os.environ.get('P4_NUFIT', '15'))
DESHIFT = os.environ.get('P4_DESHIFT', '1') == '1'
SIG = {'++': 1.0, '+-': -1.0}


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def Dc(h, J, x):
    c = 1 - 2*x*x
    return LG.legval(c, np.eye(J + 1)[J]) if h == '++' else ((1 + c)/2)**4*eval_jacobi(J - 4, 0, 8, c)


def env(h, x, k):
    y = x*x
    e = (2 - y) if h == '++' else (2 - y)/(1 - y)**4
    return e if k == 2 else e/(2 - y)


def W(h, J, x):
    c = 1 - 2*x*x
    return (1 - x*x)**2 if h == '++' else Dc('+-', J, x)/((1 - x*x)**2*eval_legendre(J, c))


def shift(J, h, k, xs, m):
    import l2_scalar_validate as V
    sg = SIG[h] if k == 3 else 1.0
    f = lambda t: sg*env(h, t, k)*Dc(h, J, t)
    out = []
    for x in xs:
        gp, gm = V.gfac(x*m)
        out.append(((gp*f(x + 1j/m) + gm*f(x - 1j/m))/(gp + gm) - f(x)).real)
    return np.array(out)


def cum(r):
    a2, a4 = r[1]/r[0], r[2]/r[0]
    a6 = r[3]/r[0] if len(r) > 3 else 0.0
    return a2, a4 - a2*a2/2, a6 - a2*a4 + a2**3/3


def pred(a, r, K=2):
    """sum_{order<=K} of a(eps) exp(-kappa_r(eps)) / r0 (l9_analysis.pred)."""
    k2, k4, k6 = cum(r)
    e2, e4, e6 = -k2, -k4 + k2*k2/2, -k6 + k2*k4 - k2**3/6
    a = list(a) + [0.0]*(4 - len(a))
    t = [a[0], a[1] + a[0]*e2, a[2] + a[1]*e2 + a[0]*e4, a[3] + a[2]*e2 + a[1]*e4 + a[0]*e6]
    return sum(t[:K + 1])/r[0]


def load(kind, J, anti=False, k3=False):
    fn = os.path.join(HERE, 'cache', 'p3_pts_%s_J%d%s%s.pkl' % (kind, J, '_anti' if anti else '', '_k3' if k3 else ''))
    R = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    return {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}


def getrec(R, D, nu):
    k = (float(D), round(float(nu), 9))
    if k in R:
        return R[k]
    c = [kk for kk in R if kk[0] == float(D) and abs(kk[1] - nu) < 1e-7*nu]
    return R[c[0]] if c else None


def lo_check(J, D, pr):
    """order-0 shapes of the odd (last) diagonal entry."""
    R2, R3 = load('C2', J), load('C3H', J)
    m = mphys(D, J)
    xs = sorted(set(round(k[1]/m, 6) for k in R2 if k[0] == D))
    rows = []
    for x in xs:
        r = getrec(R2, D, x*m)
        o = r['M0'][-1, -1]
        val = o*np.exp(r['logC'] - r['logK'] - r['slog'])/r['sS'][0]/W('+-', J, x)
        r3 = getrec(R3, D, x*m)
        c3 = (r3['M0'][-1, -1]/o*np.exp(r3['logC'] - r['logC'])*(2 - x*x)/x**2) if r3 is not None else np.nan
        rows.append((x, val, c3))
    v = np.array([r[1] for r in rows]); x = np.array([r[0] for r in rows])
    # absolute form (regular through zeros of P_J and d44): odd/(scalar/P_J) vs c d44/(1-x^2)^2, c fitted (complex LS)
    Pj = eval_legendre(J, 1 - 2*x*x); sh = Dc('+-', J, x)/(1 - x*x)**2
    u_ = v*W('+-', J, x)*Pj
    cfit = np.sum(np.conj(sh)*u_)/np.sum(np.abs(sh)**2)
    dev = np.abs(u_ - cfit*sh)/np.abs(cfit*sh).max()
    pr(f"  J={J} Delta={D:.0f}: absolute form: max |odd/(scal/P_J) - c d44/(1-x^2)^2|/max|c d44/(1-x^2)^2| = {dev.max():.2e} "
       f"(x = {x[np.argmax(dev)]:.2f}); arg(c)/pi = {np.angle(cfit)/np.pi:+.5f}")
    sel = (x >= 0.1)
    ref = v[sel][np.argmin(np.abs(x[sel] - 0.5))]
    pr(f"  J={J} Delta={D:.0f} (m={m:.1f}): C2 order-0 odd/scalar/W_+- normalised at x=0.5: " +
       " ".join(f"{xx:.2f}:{(vv/ref).real:+.5f}{(vv/ref).imag:+.0e}i" for xx, vv in zip(x, v)) +
       f";  phase of the constant/pi {np.angle(ref)/np.pi:+.5f}")
    c3 = np.array([r[2] for r in rows])
    if np.any(np.isfinite(c3)):
        pr(f"        C3/C2 x (2-y)/x^2 at order 0: " + " ".join(f"{xx:.2f}:{cc.real:+.4f}{cc.imag:+.4f}i" for xx, cc in zip(x, c3) if np.isfinite(cc)))
    return x, v


def grid(J, K=2, Ds=None, kinds=('C2', 'C3H'), pr=None, k3src=None):
    """per Delta row: m, a (+-), Gamma_{k,+-}(x) (Re), m Im for the record."""
    h = '+-'
    RC = {kd: load(kd, J, k3=(K == 3 if k3src is None else k3src)) for kd in kinds}
    Ds = Ds or sorted(set(k[0] for k in RC['C2'] if k[0] < DREF))
    out = {}
    for D in Ds:
        m = mphys(D, J); mr = mphys(DREF, J)
        xs = sorted(set(round(k[1]/m, 6) for k in RC['C2'] if k[0] == D))
        xs = [x for x in xs if all(getrec(RC[kd], D, x*m) is not None and getrec(RC[kd], DREF, x*m) is not None for kd in kinds)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}
        c0 = None
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G = []
            for x in xs:
                nu = x*m; st = getrec(RC[kd], D, nu); rf = getrec(RC[kd], DREF, nu); xr = nu/mr
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                ords = ('M0', 'M2', 'M4', 'M6') if K == 3 else ('M0', 'M2', 'M4')
                v = np.exp(dl)*pred([st[o][-1, -1] for o in ords], [rf[o][-1, -1] for o in ords], K)*fk
                G.append((SIG[h] if k == 3 else 1.0)*v*env(h, xr, k)*Dc(h, J, xr))
            rows[kd] = (k, np.array(G))
        sel = (xs <= XFIT + 1e-9) & (xs*m >= NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]
        def fit(g):
            A_ = np.array([Dc(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A_, (g/env(h, xs, 2))[sel], rcond=None)[0][0]
        a = fit(rows['C2'][1].real)
        a_raw = a
        if DESHIFT:
            sh = {k: shift(J, h, k, xs, m) for k in (2, 3)}
            for _ in range(3):
                a = fit(rows['C2'][1].real - a*sh[2])
            rows = {kd: (k, G - a*sh[k]) for kd, (k, G) in rows.items()}
        rec = dict(m=m, xs=xs, a=a, a_raw=a_raw, c0=c0)
        for kd, (k, G) in rows.items():
            sg = SIG[h] if k == 3 else 1.0
            rec[kd] = m*m*(G.real/(a*env(h, xs, k)) - sg*Dc(h, J, xs))
            rec[(kd, 'im')] = m*G.imag/(a*env(h, xs, k))
        out[D] = rec
    return out


if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1:]] or [5, 7, 9]
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr("# p4 (a): LO check of the parity-odd Ĉ_2 action against the stripped +- density, and Ĉ_3/Ĉ_2")
    for J in Js:
        for D in (1600.0, 3200.0):
            try:
                lo_check(J, D, pr)
            except Exception as ex:
                pr(f"  J={J} Delta={D}: {ex}")
    pr(f"# p4 (b): Gamma_{{k,+-}}(x; m) (c2 units x LO shape), linear chain K2, Delta_ref = {DREF:.0f}, de-shift {DESHIFT}")
    store = {}
    KS = tuple(os.environ.get('P4_KINDS', 'C2,C3H').split(','))
    KK = int(os.environ.get('P4_K', '2'))
    for J in Js:
        g = grid(J, K=KK, kinds=KS)
        store[J] = g
        for D in sorted(g):
            r = g[D]
            pr(f"  J={J} Delta={D:.0f} m={r['m']:.1f} a={r['a']:.6g} (no de-shift {r['a_raw']:.6g})  x = " + " ".join(f"{x:.2f}" for x in r['xs']))
            for kd in ('C2', 'C3H'):
                if kd in r:
                    pr(f"    {kd:3s} Gamma: " + " ".join(f"{v:+8.1f}" for v in r[kd]))
                    pr(f"    {kd:3s} m Im : " + " ".join(f"{v:+8.2f}" for v in r[(kd, 'im')]))
    tag = ('_k3' if KK == 3 else '') + ('' if DESHIFT else '_noshift')
    pickle.dump(store, open(os.path.join(HERE, 'cache', 'p4_grid%s.pkl' % tag), 'wb'))
    open(os.path.join(HERE, 'cache', 'p4_odd_analysis%s_%s.txt' % (tag, '_'.join(KS))), 'w').write("\n".join(lines) + "\n")
