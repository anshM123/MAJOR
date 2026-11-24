"""
gm_ball8 (GAP_MORSE.md §4.4): gm_ball with an 8th Taylor-model variable delta = xi - xi0 (complex disc |delta| <= RHOX), for
certificates that cover a xi-INTERVAL.  Degree bounds are bi-graded: the part of v-degree k and delta-degree j is bounded by
[p_{k,j}]_B RHO^k RHOX^j (Bombieri norm in v only), summed within each total degree (submultiplicative).  No project-module loader.

Original docstring --- gm_ball (GAP_MORSE.md): copy of gapc3/gc_ball.py.  ONLY CHANGE: with GM_NORM=bomb (default) the degree bounds are Bombieri
norms, so every Taylor model is a certified enclosure on the HERMITIAN BALL {v in C^7 : |v|_2 <= rho} instead of the polydisc
(sup_{|v|<=rho} |p_k(v)| <= [p_k]_B rho^k and [p q]_B <= [p]_B [q]_B; GAP_MORSE.md §2).  GM_NORM=poly restores the polydisc.

Original docstring --- gc_ball (GAP_C3.md): Arb (python-flint acb) BALL TPS backend with the l1_tps interface, so that the project's crystal integrands
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
NV = 8                          # 7 Morse variables v (Hermitian ball radius RHO) + delta = xi - xi0 (disc radius RHOX)
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


NORM = os.environ.get('GM_NORM', 'bomb')
# Bombieri weights alpha!/|alpha|!  (exact rationals as arb)
_BW = np.array([arb(math.prod(math.factorial(e) for e in m))/math.factorial(sum(m)) for m in MONS], dtype=object)
_vsqrt = np.frompyfunc(lambda z: arb(z.upper()).sqrt(), 1, 1)   # upper-bound sqrt (a ball [0, x] has no sqrt in Arb)


RHOX = [0.01]
_KV = np.array([sum(m[:7]) for m in MONS]); _KX = np.array([m[7] for m in MONS])
_GID = _KV*(DEG + 1) + _KX                       # group (v-degree, delta-degree)
_NG = (DEG + 1)**2


def bigraded(p):
    """p: (..., NM) -> (..., DEG+1, DEG+1) arb: [p_{k,j}]_B (Bombieri norm in v of the part of v-degree k, delta-degree j)."""
    ab = _vabs(p)
    w = ab*ab*_BW
    out = np.empty(p.shape[:-1] + (_NG,), dtype=object); out.fill(arb(0))
    np.add.at(np.moveaxis(out, -1, 0), _GID, np.moveaxis(w, -1, 0))
    return _vsqrt(out).reshape(p.shape[:-1] + (DEG + 1, DEG + 1))


def vbounds(p, rx=None):
    """(..., DEG+1) arb: sum_j [p_{k,j}]_B rx^j  (uniform-in-delta Bombieri norm of the v-degree-k part, radius-1 in v)."""
    rx = arb(RHOX[0] if rx is None else rx)
    bg = bigraded(p)
    rxp = np.array([rx**j for j in range(DEG + 1)], dtype=object)
    return np.add.reduce(bg*rxp, axis=-1)


def bombieri(p):
    """total-degree Bombieri-type bound at radius 1 in v and RHOX in delta: sum_{k+j=d} [p_{k,j}] RHOX^j."""
    bg = bigraded(p)
    rx = arb(RHOX[0])
    out = np.empty(p.shape[:-1] + (DEG + 1,), dtype=object); out.fill(arb(0))
    for k in range(DEG + 1):
        for j in range(DEG + 1 - k):
            out[..., k + j] = out[..., k + j] + bg[..., k, j]*rx**j
    return out


def degbounds(p):
    """p: polynomial part (..., NM) -> (..., DEG+1) object array of arb: an upper bound B_k for sup |p_k| on the domain
    (Bombieri: [p_k]_B rho^k on the Hermitian ball; poly: sum_{|alpha|=k} |c_alpha| rho^k on the polydisc)."""
    bg = bigraded(p)
    rx = arb(RHOX[0]); r = arb(RHO[0])
    out = np.empty(p.shape[:-1] + (DEG + 1,), dtype=object); out.fill(arb(0))
    for k in range(DEG + 1):
        for j in range(DEG + 1 - k):
            out[..., k + j] = out[..., k + j] + bg[..., k, j]*r**k*rx**j
    return out


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
    # GAP_MORSE fix: EXACT series coefficients (-1)^{n+1}/n (gc_ball used rounded floats 1/3, 1/5, ..., which the remainder
    # bound does not account for; effect ~1e-17 relative, but outside the enclosure)
    s, a0 = _compose(a, [acb(0)] + [acb((-1)**(n + 1))/n for n in range(1, DEG + 1)],
                     lambda q: q**(DEG + 1)/((DEG + 1)*(1 - q)))
    l0 = _vlog(a0) if branch0 is None else obj(np.broadcast_to(branch0, a0.shape))
    s[..., 0] = s[..., 0] + l0
    return s


@_nd
def powr(a, p):
    """a^p = exp(p log a0) (1 + t)^p; remainder |binom(p, D+1)| q^(D+1) (1-q)^(p-D-1) (integral form, real p <= D+1)."""
    a = obj(a)
    assert p <= DEG + 1
    c = [_binom(p, n) for n in range(DEG + 1)]          # GAP_MORSE: exact binomial coefficients (Arb), not floats
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
    """z_i = zc_i + sum_{j<7} L_ij v_j, zc_i a Taylor model (in delta) or a number."""
    out = []
    for i in range(7):
        a = zc[i].copy() if isinstance(zc[i], np.ndarray) else const(zc[i])
        for j in range(7):
            a[1 + j] = a[1 + j] + A(L[i, j])
        out.append(a)
    return out


def var(i):
    a = zeros(N); a[1 + i] = A(1)
    return a


assert all(MONS[1 + i][i] == 1 for i in range(NV))


