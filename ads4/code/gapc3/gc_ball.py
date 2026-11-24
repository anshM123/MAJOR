"""
gc_ball (GAP_C3.md): Arb (python-flint acb) BALL TPS backend with the l1_tps interface, so that the project's crystal integrands
(l1_crystal scalar, l5_spin / w0_ee_common EE, l9_points C2/C3H, p3_odd_action odd J) run UNCHANGED in interval arithmetic.

A "TPS" here is a numpy object array (..., 2*NM) of acb:  [0:NM]  = Taylor coefficients (Morse variables w, total degree <= DEG) at the
exact crystal point z*  ("centre half");  [NM:2NM] = Taylor coefficients at a point p that ranges over a complex polydisc P around z*
("ball half": every coefficient is an enclosure valid for EVERY p in P simultaneously).  Arithmetic is ball arithmetic throughout.

Branches: the centre half uses the principal log at z* (as l1_tps / n1_pw).  The ball half uses the analytic continuation from z*:
log a(p) = log a(z*) + Log(a(p)/a(z*)); this equals the continuation along the segment [z*, p] (and the function is analytic on P) as
soon as the enclosure of a(p)/a(z*) over P lies in Re > 0.  Every log / power checks this; failures are recorded in BAD (the caller
must reject the polydisc if BAD is non-empty).  Division by a ball containing 0 gives a non-finite enclosure (also detected).

Fresh copies of the project modules are exec'd with sys.modules['l1_tps'] -> this module (mechanism of spinfunc/n1_pw.py); the only
source changes are dtype ones (np.zeros(..., complex) -> object) and the K2-only assert of p3_odd_action; S.bm / fourier / null_circle
/ kappa are replaced by exact ball versions.  Nothing in the project is modified.
"""
import os, sys, types, itertools, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF)
import numpy as np
from flint import acb, arb, ctx

ctx.prec = int(os.environ.get('GC_PREC', '106'))
NV = 7
DEG = int(os.environ.get('GC_DEG', '3'))
MONS = []
for d_ in range(DEG + 1):
    for c in itertools.combinations_with_replacement(range(NV), d_):
        e = [0]*NV
        for i in c:
            e[i] += 1
        MONS.append(tuple(e))
NM = len(MONS)
N = NM + 1                      # Taylor coefficients + one remainder ball (Taylor model on P_rho)
EXP = np.array(MONS, dtype=np.int64)
DEGS = EXP.sum(1)
_B = DEG + 1
_CODE = (EXP*(_B**np.arange(NV))).sum(1)
_ORD = np.argsort(_CODE)
_SC = _CODE[_ORD]
_I, _J = [], []
for a_ in range(NM):
    bs = np.where(DEGS <= DEG - DEGS[a_])[0]
    _I.append(np.full(len(bs), a_)); _J.append(bs)
_I = np.concatenate(_I); _J = np.concatenate(_J)
_K = _ORD[np.searchsorted(_SC, _CODE[_I] + _CODE[_J])]
ZERO = acb(0)
BAD = []
WMAX = [0.0]
RHO = [0.1]                     # polydisc radius in the Morse variables (set by the caller before building)
_DSTART = [int(np.searchsorted(DEGS, k)) for k in range(DEG + 1)]


def A(c):
    """python / numpy number or acb -> acb (floats are exact binary numbers)."""
    if isinstance(c, acb):
        return c
    if isinstance(c, arb):
        return acb(c)
    c = complex(c)
    return acb(c.real, c.imag)


_vA = np.frompyfunc(A, 1, 1)
_vexp = np.frompyfunc(lambda z: z.exp(), 1, 1)
_vlog = np.frompyfunc(lambda z: z.log(), 1, 1)
_vabs = np.frompyfunc(lambda z: abs(z), 1, 1)


def obj(x):
    x = np.asarray(x)
    if x.dtype == object:
        return x
    return _vA(x).astype(object)


def zeros(shape):
    out = np.empty(shape, dtype=object)
    out.fill(ZERO)
    return out


def disc(X):
    """complex box enclosing the disc |w| <= X (X: arb, upper bound used)."""
    u = arb(0, X.upper()) if isinstance(X, arb) else arb(0, X)
    return acb(u, u)


_vdisc = np.frompyfunc(disc, 1, 1)


def const(c):
    a = zeros(N); a[0] = A(c)
    return a


def _rhopow():
    return [arb(RHO[0])**k for k in range(2*DEG + 2)]


def degbounds(p):
    """p: polynomial part (..., NM) -> (..., DEG+1) object array of arb: B_k = sum_{|alpha|=k} |c_alpha| rho^k."""
    ab = _vabs(p)
    s = np.add.reduceat(ab, _DSTART, axis=-1)
    rp = _rhopow()
    return s*np.array(rp[:DEG + 1], dtype=object)


def _sumlast(x):
    return np.add.reduce(x, axis=-1)


def _mulh(a, b):
    """truncated product of polynomial parts (..., NM), broadcast."""
    sh = np.broadcast_shapes(a.shape[:-1], b.shape[:-1])
    a = np.broadcast_to(a, sh + (NM,)); b = np.broadcast_to(b, sh + (NM,))
    w = a[..., _I]*b[..., _J]
    out = zeros(sh + (NM,))
    np.add.at(np.moveaxis(out, -1, 0), _K, np.moveaxis(w, -1, 0))
    return out


def _hi_pairs():
    return [(k, l) for k in range(DEG + 1) for l in range(DEG + 1) if k + l > DEG]


_HP = _hi_pairs()


def mul(a, b):
    """Taylor-model product on P_rho: (P1 + I1)(P2 + I2)."""
    a = obj(a); b = obj(b)
    if a.shape[-1] != N or b.shape[-1] != N:
        return a*b
    if a.ndim == 1 and b.ndim == 1:
        return mul(a[None], b[None])[0]
    pa, ra = a[..., :NM], a[..., NM]
    pb, rb = b[..., :NM], b[..., NM]
    prod = _mulh(pa, pb)
    Ba = degbounds(pa); Bb = degbounds(pb)
    hi = 0
    for (k, l) in _HP:
        hi = hi + Ba[..., k]*Bb[..., l]
    ta, tb = _sumlast(Ba), _sumlast(Bb)
    aa, ab_ = _vabs(ra), _vabs(rb)
    rem = _vdisc(hi + ta*ab_ + tb*aa + aa*ab_)
    sh = prod.shape[:-1]
    out = np.empty(sh + (N,), dtype=object)
    out[..., :NM] = prod
    out[..., NM] = np.broadcast_to(rem, sh)
    return out


bm = mul


def _nd(f):
    """run f on >= 2-d object arrays (1-d Taylor models are lifted to shape (1, N) and back)."""
    def g(*args, **kw):
        args = [obj(x) for x in args[:1]] + list(args[1:])
        if args[0].ndim == 1:
            r = f(args[0][None], *args[1:], **kw)
            return r[0] if isinstance(r, np.ndarray) else r
        return f(*args, **kw)
    return g


def _range_ratio(a):
    """q = sup_P |a - a0| / |a0| (upper, float) for each element; a0 = value at the crystal."""
    p, r = a[..., :NM], a[..., NM]
    a0 = p[..., 0]
    B = degbounds(p)
    dev = _sumlast(B[..., 1:]) + _vabs(r)
    q = dev/_vabs(a0)
    return q, a0


def _compose(a, coefs, remfun):
    """f(a) = sum_{n<=DEG} coefs[n] t^n + R,  t = (a - a0)/a0 (Taylor model, |t| <= q < 1 on P), |R| <= remfun(q)."""
    q, a0 = _range_ratio(a)
    for z in np.asarray(q).ravel():
        if not (z < 1):
            BAD.append(('q>=1', str(z)))
            break
        v = float(z.upper())
        if v > WMAX[0]:
            WMAX[0] = v
    t = a/a0[..., None]
    t[..., 0] = ZERO
    out = zeros(a.shape); out[..., 0] = A(coefs[0])
    p = zeros(a.shape); p[..., 0] = A(1)
    for n in range(1, DEG + 1):
        p = mul(p, t)
        out = out + A(coefs[n])*p
    out[..., NM] = out[..., NM] + np.frompyfunc(lambda z: disc(remfun(z)) if z < 1 else acb(arb('nan'), arb('nan')), 1, 1)(q)
    return out, a0


@_nd
def recip(a):
    a = obj(a)
    s, a0 = _compose(a, [(-1.0)**n for n in range(DEG + 1)], lambda q: q**(DEG + 1)/(1 - q))
    return s/a0[..., None]


def div(a, b):
    return mul(a, recip(b))


def _binom(p, n):
    c = arb(1)
    for k in range(n):
        c = c*(arb(p) - k)/(k + 1)
    return c


@_nd
def log(a, branch0=None):
    """log a = log a0 (principal, or branch0) + Log(1 + t), the continuation from the crystal (|t| < 1 on P)."""
    a = obj(a)
    s, a0 = _compose(a, [0.0] + [(-1.0)**(n + 1)/n for n in range(1, DEG + 1)],
                     lambda q: q**(DEG + 1)/((DEG + 1)*(1 - q)))
    l0 = _vlog(a0) if branch0 is None else obj(np.broadcast_to(branch0, a0.shape))
    s[..., 0] = s[..., 0] + l0
    return s


@_nd
def powr(a, p):
    """a^p = exp(p log a0) (1 + t)^p; remainder |binom(p, D+1)| q^(D+1) (1-q)^(p-D-1) (integral form, real p <= D+1)."""
    a = obj(a)
    assert p <= DEG + 1
    c = [1.0]
    for n in range(1, DEG + 1):
        c.append(c[-1]*(p - n + 1)/n)
    bnd = abs(_binom(p, DEG + 1))
    s, a0 = _compose(a, c, lambda q: bnd*q**(DEG + 1)*(1 - q)**(arb(p) - DEG - 1))
    return s*_vexp(A(p)*_vlog(a0))[..., None]


@_nd
def exp(a):
    """exp of a Taylor model: e^{a0} (sum_{n<=D} d^n/n! + R), |R| <= e^{|d|} |d|^{D+1}/(D+1)!, |d| <= sup_P |a - a0|."""
    a = obj(a)
    p, r = a[..., :NM], a[..., NM]
    a0 = p[..., 0]
    B = degbounds(p)
    dv = _sumlast(B[..., 1:]) + _vabs(r)
    d = a.copy(); d[..., 0] = ZERO
    s = zeros(a.shape); s[..., 0] = A(1)
    q = zeros(a.shape); q[..., 0] = A(1)
    for n in range(1, DEG + 1):
        q = mul(q, d)
        s = s + q/acb(math.factorial(n))
    s[..., NM] = s[..., NM] + np.frompyfunc(lambda z: disc(z.exp()*z**(DEG + 1)/math.factorial(DEG + 1)), 1, 1)(dv)
    return s*_vexp(a0)[..., None]


def bound(a):
    a = obj(a)
    if a.ndim == 1:
        a = a[None]
    """sup_P |a| (arb upper) of a Taylor model element; and per-degree bounds."""
    a = obj(a)
    B = degbounds(a[..., :NM])
    return _sumlast(B) + _vabs(a[..., NM]), B


def homog(a, d_):
    b = obj(a).copy()
    m = np.concatenate([DEGS != d_, [True]])
    b[..., m] = ZERO
    return b


def linsub_ball(zc, L):
    """Taylor-model variables z_i = zc_i + sum_j L_ij w_j on P_rho (zero remainder)."""
    out = []
    for i in range(NV):
        a = zeros(N)
        a[0] = zc[i]
        for j in range(NV):
            a[1 + j] = A(L[i, j])
        out.append(a)
    return out


assert all(MONS[1 + i][i] == 1 for i in range(NV))


# ----------------------------------------------------------------------------------------------- module loader (n1_pw mechanism)
def _load(name, swaps, patches=()):
    saved = {k: sys.modules.get(k) for k in swaps}
    sys.modules.update(swaps)
    try:
        M = types.ModuleType(name + '_ball'); M.__dict__['__file__'] = os.path.join(SF, name + '.py')
        src = open(os.path.join(SF, name + '.py')).read()
        for old, new in patches:
            assert old in src, (name, old)
            src = src.replace(old, new)
        exec(compile(src, name + '_ball', 'exec'), M.__dict__)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    return M


_self = sys.modules[__name__]
os.environ.setdefault('L1_DEG', '6')
import l1_tps as _TT                       # noqa: F401  (real engines cached before the swap)
import l1_crystal as _Cr_real              # noqa: F401
import l5_spin as _S_real
import l5_shadow_matrix as _SM_real        # noqa: F401
import l6_spin_action as _A6_real          # noqa: F401

Cr = _load('l1_crystal', {'l1_tps': _self})
S = _load('l5_spin', {'l1_tps': _self},
          patches=[("res = np.zeros((sol.shape[1], nb, n5, N), complex)", "res = np.zeros((sol.shape[1], nb, n5, N), object)")])


def _null_circle(l):
    n = 2*l + 1
    out = np.empty((n, 3), dtype=object)
    for k in range(n):
        th = arb.pi()*2*k/n
        c, s = th.cos(), th.sin()
        out[k] = [acb(1 - c), acb(1 + c), acb(s)]
    return out


def _kappa(l):
    return np.array([acb((-1)**(l + n)*math.comb(2*l, l + n))/acb(2**l) for n in range(-l, l + 1)], dtype=object)


def _fourier(vals, axis, l):
    """Fourier coefficients n = -l..l of samples on null_circle(l) along axis (exact DFT in ball arithmetic)."""
    nn = 2*l + 1
    v = np.moveaxis(obj(vals), axis, 0)
    out = []
    for n in range(-l, l + 1):
        acc = None
        for k in range(nn):
            th = -arb.pi()*2*k*n/nn
            w = acb(th.cos(), th.sin())
            t = v[k]*w
            acc = t if acc is None else acc + t
        out.append(acc/acb(nn))
    return np.moveaxis(np.stack(out, axis=0), 0, axis)


def _lin_samples(fbasis, W):
    fb = np.stack([obj(f) for f in fbasis], axis=0)
    return np.tensordot(obj(W), fb, axes=([1], [0]))


S.bm = mul
S.null_circle = _null_circle
S.kappa = _kappa
S.fourier = _fourier
S.lin_samples = _lin_samples
S.conserved = _S_real.conserved            # float coefficient vectors (real engine); see GAP_C3.md for their status
A6 = _load('l6_spin_action', {'l1_tps': _self, 'l1_crystal': Cr, 'l5_spin': S})
_W0src = [("pref = np.zeros((na, nc, T.N), complex)", "pref = np.zeros((na, nc, T.N), object)")]
W0 = _load('w0_ee_common', {}, patches=_W0src)
P9 = _load('l9_points', {'l1_tps': _self, 'l1_crystal': Cr, 'l5_spin': S, 'l6_spin_action': A6},
           patches=[("pref = np.zeros((na, nc, T.N), complex)", "pref = np.zeros((na, nc, T.N), object)")])


def load_p3():
    """odd-J builder (p3_odd_action) on the ball backend, with the EE kind of w0 installed."""
    import p1_tto_odd as _P1r  # noqa: F401
    import p2_shadow_odd as _P2r  # noqa: F401
    P3 = _load('p3_odd_action', {'l1_tps': _self, 'l1_crystal': Cr, 'l5_spin': S, 'l6_spin_action': A6},
               patches=[("assert T.DEG in (6, 8)", "assert True")])
    return P3
