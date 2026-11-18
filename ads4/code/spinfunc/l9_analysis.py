"""
l9_analysis (M6 round 5, doc 23.2): graviton-specific absolute corrections of C^_2 AND C^_3 per helicity and the helicity mixing, in the
LINEAR (prefactor) form, with the nu-only factor removed by the heavy reference ORDER BY ORDER:

  p_red = [sum_j a_j eps^j] * exp(-kappa^ref(eps)),  kappa^ref = cumulant series of the reference element (same kind, ++ for the off-
  diagonal; the nu-only factor is helicity-independent, doc 22.1), truncated at the Laplace order K (2 or 3).
  Linear form: regular through zeros of the LO shapes D_++ = P_J, D_+- = d44 (the cumulant/log form of doc 22 is singular there, and on
  the exact scalar skeleton it is not convergent at m ~ 100: doc 23.1, l9_scalar_k3).

Normalisation (per state; per-state helicity constants a_h are COMMON to k = 2, 3 -- they come from the <TTO> structure normalisation and
beta, identical for both functionals; the relative 1/m^2 of C^_3 is the u m^2 of doc R.6, verified at LO: doc 23.3):
  value_{k,h}(x)  = e^{dlogC} p_red / ref_0 x (m^2/m_ref^2 if k = 3)
  G_{k,h}(x)      = sigma_h(k=3) value x shape_{k,h}(x_ref),   shape_{2,h} = env_h D_h, shape_{3,h} = env_h D_h/(2-y), env_++ = 2-y,
                    env_+- = (2-y)/(1-y)^4   (F units of l8_neargap / h2)
  a_h             : Re G_{2,h}(x)/env_h(x) = a_h (D_h(x) + b x^2 + d x^4) on x <= XFIT  (graviton k = 2 Regge value := 1, as l8)
  Gamma_{k,h}(x)  = m^2 [Re G_{k,h}(x)/(a_h env_{k,h}(x)) - sigma D_h(x)]     (c2 units x LO shape; regular through zeros)
  Gamma_{k,off}   = m Re[ e^{dlogC} p_red,off / ref_++ x sqrt(ref2_++/ref2_+-) ] sqrt(shape_+ shape_-)(x_ref)/sqrt(a_+ a_-)/env_{k,off}(x)
                    (env_{2,off} = sqrt(env_+ env_-), env_{3,off} = env_{2,off}/(2-y); units of l8's offp: F_+- += psi w env_off Gamma/m)
  Im parts are kept for the record (dropped by the Re prescription).
De-shift (L9_DESHIFT=1, default; doc 23.2): the CLS-D crystal computes Psi_nu ~ g C(nu+i) + g~ C(nu-i) (CHMRS C_nu = the object of
the exact-skeleton base model and of the light side); the shift of the LO density is subtracted (scalar g assumed for the graviton [I]).
Usage: python3 l9_analysis.py [J ...]  -> cache/l9_grid.pkl, cache/l9_analysis.txt  (K2 vs K3, C2 and C3H/C3V)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import eval_legendre
import l8_analysis as L
import l7_dc2_analysis as An

DREF = 20000.0
XFIT = 0.25
SIG = {'++': 1.0, '+-': -1.0}


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def _load(kind, J, sfx):
    fn = os.path.join(HERE, 'cache', 'l9_pts_%s_J%d%s.pkl' % (kind, J, sfx))
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


class Recs:
    """records for (kind, J) at Laplace order K (2: DEG-8 records truncated, else DEG-6 records, else (C2) l8/l7)."""
    def __init__(self, kind, J, K):
        self.K = K
        self.R = {}
        if K == 2:
            if kind == 'C2':
                for k, r in L.records(J).items():
                    self.R[k] = dict(r, M6=0*r['M0'])
            self.R.update(_load(kind, J, '_k2'))
        self.R.update(_load(kind, J, ''))

    def get(self, D, nu):
        k = (float(D), round(float(nu), 9))
        if k in self.R:
            return self.R[k]
        c = [kk for kk in self.R if kk[0] == float(D) and abs(kk[1] - nu) < 1e-7*nu]
        return self.R[c[0]] if c else None


def cum(r):
    """cumulant coefficients (k2, k4, k6) of a series r = [r0, r2, r4, r6]."""
    a2, a4, a6 = r[1]/r[0], r[2]/r[0], r[3]/r[0]
    return a2, a4 - a2*a2/2, a6 - a2*a4 + a2**3/3


def pred(a, r, K):
    """sum_{order<=K} of a(eps) exp(-kappa_r(eps)), divided by r0."""
    k2, k4, k6 = cum(r)
    e2, e4, e6 = -k2, -k4 + k2*k2/2, -k6 + k2*k4 - k2**3/6
    t = [a[0], a[1] + a[0]*e2, a[2] + a[1]*e2 + a[0]*e4, a[3] + a[2]*e2 + a[1]*e4 + a[0]*e6]
    return sum(t[:K + 1])/r[0]


def Dsh(h, J, x):
    return An.Dh(h, J, x)


def env(h, x, k):
    y = x*x
    e = (2 - y) if h == '++' else (2 - y)/(1 - y)**4
    return e if k == 2 else e/(2 - y)


DESHIFT = os.environ.get('L9_DESHIFT', '1') == '1'
NUFIT = float(os.environ.get('L9_NUFIT', '15'))


def Dc(h, J, x):
    from numpy.polynomial import legendre as LG
    from scipy.special import eval_jacobi
    c = 1 - 2*x*x
    return LG.legval(c, np.eye(J + 1)[J]) if h == '++' else ((1 + c)/2)**4*eval_jacobi(J - 4, 0, 8, c)


def shift(J, h, k, xs, m):
    """CLS-D vs CHMRS: X(nu) = g(nu) C(nu+i) + g(-nu) C(nu-i) (doc 20.2, scalar g of l2), applied to the LO density
    f = sigma env_{k,h} D_h (F units, k = 3 incl. sigma): returns [g f(x+i/m) + g~ f(x-i/m)]/(g+g~) - f(x) (real part)."""
    import l2_scalar_validate as V
    sg = SIG[h] if k == 3 else 1.0
    f = lambda t: sg*env(h, t, k)*Dc(h, J, t)
    out = []
    for x in xs:
        gp, gm = V.gfac(x*m)
        out.append(((gp*f(x + 1j/m) + gm*f(x - 1j/m))/(gp + gm) - f(x)).real)
    return np.array(out)


def helmat(r, B, K):
    return [B['Ti'] @ r[o] @ B['Ti'].T for o in ('M0', 'M2', 'M4', 'M6')[:K + 1]] + [0*r['M0']]*(3 - K)


def grid(J, K, kinds=('C2', 'C3H'), XS=None, Ds=None, pr=None):
    """per Delta row: m, a_h, Gamma_{k,e}(x) (Re, Im) for k in kinds; XS default: all x with complete data."""
    RC = {k: Recs(k, J, K) for k in kinds}
    Rb = dict(L.records(J)); Rb.update({kk: v for kk, v in _load('C2', J, '_k2').items()}); Rb.update(_load('C2', J, ''))
    Bref = L.basis(J, DREF, Rb) if J >= 4 else dict(Ti=np.eye(1))
    hels = ('++', '+-') if J >= 4 else ('++',)
    Ds = Ds or sorted(set(k[0] for k in RC['C2'].R if k[0] < DREF))
    out = {}
    for D in Ds:
        m = mphys(D, J); mr = mphys(DREF, J)
        B = L.basis(J, D, Rb) if J >= 4 else dict(Ti=np.eye(1))
        if B is None:
            continue
        xs = XS or sorted(set(round(k[1]/m, 6) for k in RC['C2'].R if k[0] == D))
        xs = [x for x in xs if all(RC[kd].get(D, x*m) is not None and RC[kd].get(DREF, x*m) is not None for kd in kinds)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}
        c0 = None
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G = {h: [] for h in hels}; Go = []
            for x in xs:
                nu = x*m; st = RC[kd].get(D, nu); rf = RC[kd].get(DREF, nu); xr = nu/mr
                Ms, Mr = helmat(st, B, K), helmat(rf, Bref, K)
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if 'c0' not in locals() or c0 is None:
                    c0 = dl.real                                   # common per-row scale (absorbed in a_h; same for all kinds/h)
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                for j, h in enumerate(hels):
                    v = np.exp(dl)*pred([Mo[j, j] for Mo in Ms], [Mo[j, j] for Mo in Mr], K)*fk
                    G[h].append((SIG[h] if k == 3 else 1.0)*v*env(h, xr, k)*Dsh(h, J, xr))
                if len(hels) == 2:
                    if kd == 'C2':
                        r2 = Mr
                    else:
                        r2 = helmat(RC['C2'].get(DREF, nu), Bref, K)
                    sq = np.sqrt(abs((r2[0][0, 0]/r2[0][1, 1]).real))
                    off = 0.5*np.array([Mo[0, 1] + Mo[1, 0] for Mo in Ms])
                    v = np.exp(dl)*pred(off, [Mo[0, 0] for Mo in Mr], K)*fk*sq
                    Go.append(v*np.sqrt(env('++', xr, k)*env('+-', xr, k)*Dsh('++', J, xr)*Dsh('+-', J, xr)))
            rows[kd] = (k, {h: np.array(G[h]) for h in hels}, np.array(Go))
        # per-state constants a_h from C2 (fit on x <= XFIT with nu >= NUFIT: below, the crystal misses the non-perturbative 1/nu
        # structure of g(nu), doc 23.2), iterated with the de-shift
        a = {}
        sel = (xs <= XFIT + 1e-9) & (xs*m >= NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]
        def fit(h, g):
            A = np.array([Dsh(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A, (g/env(h, xs, 2))[sel], rcond=None)[0][0]
        for h in hels:
            a[h] = fit(h, rows['C2'][1][h].real)
        if DESHIFT:
            sh = {(k, h): shift(J, h, k, xs, m) for k in (2, 3) for h in hels}
            for _ in range(3):
                for h in hels:
                    a[h] = fit(h, rows['C2'][1][h].real - a[h]*sh[(2, h)])
            for kd, (k, G, Go) in rows.items():
                for h in hels:
                    G[h] = G[h] - a[h]*sh[(k, h)]
        rec = dict(m=m, xs=xs, a=a)
        for kd, (k, G, Go) in rows.items():
            for h in hels:
                sg = SIG[h] if k == 3 else 1.0
                rec[(kd, h)] = m*m*(G[h].real/(a[h]*env(h, xs, k)) - sg*Dsh(h, J, xs))
                rec[(kd, h, 'im')] = m*G[h].imag/(a[h]*env(h, xs, k))                    # m Im (relative), for the record
            if len(hels) == 2:
                eo = np.sqrt(env('++', xs, k)*env('+-', xs, k))
                rec[(kd, 'off')] = m*Go.real/np.sqrt(a['++']*a['+-'])/eo
                rec[(kd, 'off', 'im')] = m*Go.imag/np.sqrt(a['++']*a['+-'])/eo
        out[D] = rec
    return out


def report(J, pr, kinds=('C2', 'C3H')):
    g2 = grid(J, 2, kinds)
    try:
        g3 = grid(J, 3, kinds)
    except Exception:
        g3 = {}
    hels = ('++', '+-') if J >= 4 else ('++',)
    pr(f"J={J}: rows K2 {sorted(g2)}; K3 {sorted(g3)}")
    for D in sorted(g2):
        r = g2[D]
        pr(f"  Delta={D:.0f} m={r['m']:.1f}  a = " + ", ".join(f"{h}:{r['a'][h]:.6g}" for h in hels) + "   x = " +
           " ".join(f"{x:.2f}" for x in r['xs']))
        for kd in kinds:
            for e in hels + (('off',) if J >= 4 else ()):
                s = f"    {kd:3s} {e:3s} K2: " + " ".join(f"{v:+8.1f}" for v in r[(kd, e)])
                pr(s)
                if D in g3 and len(g3[D]['xs']) == len(r['xs']):
                    pr(f"    {kd:3s} {e:3s} K3: " + " ".join(f"{v:+8.1f}" for v in g3[D][(kd, e)]))
                    pr(f"    {kd:3s} {e:3s} mIm:" + " ".join(f"{v:+8.1f}" for v in g3[D][(kd, e, 'im')]))
    return g2, g3


if __name__ == '__main__':
    Js = [int(a) for a in sys.argv[1:]] or [0, 2, 4, 6, 8, 10]
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    store = {}
    for J in Js:
        g2, g3 = report(J, pr)
        store[J] = {2: g2, 3: g3}
    fn = os.path.join(HERE, 'cache', 'l9_grid.pkl')
    old = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    old.update(store); pickle.dump(old, open(fn, 'wb'))
    open(os.path.join(HERE, 'cache', 'l9_analysis_J%s.txt' % '_'.join(map(str, Js))), 'w').write("\n".join(lines) + "\n")
