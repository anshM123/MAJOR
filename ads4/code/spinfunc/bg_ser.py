"""
bg_ser (BELOW_FOLD2.md): exact-coefficient algebra for the analytic-radial 1/Delta expansion of the crystal integrand on the real
future-cone section D of bf_core.

Variables (BELOW_FOLD §1-2): raw z = a*(scaled), a = 1/Delta; radial X1 = p1^2, X2 = p2^2, Y1 = q1^2, Y2 = q2^2 (p, q > 0),
s in (-1,1) with xp = a s p1 p2, phi in (-pi/2, pi/2), and alpha = psi + phi in (-pi, pi] (psi in (-pi, pi], sheet n0 = +1; this
covers both sheets n0 = +-1 of bf_core).  A coefficient is a Laurent polynomial in the 12 "atoms"
     0 p1, 1 p2, 2 q1, 3 q2, 4 s, 5 w = 1 - s^2, 6 cf = cos phi, 7 sf = sin phi, 8 ca = cos alpha, 9 sa = sin alpha,
     10 u = 1 - cos alpha, 11 v = 1 + cos alpha
(redundant but closed under products; canonicalised by s^2 -> 1 - w, sf^2 -> 1 - cf^2, sa^2 -> 1 - ca^2).  A quantity is an
ASer: a truncated Laurent series in a, sum_{k=0..K} a^{d+k} c_k, c_k = Poly.  Division / logs / non-integer powers need the leading
coefficient c_0 to be a single atom monomial (recognised numerically, `atomise`); this is the case for every divisor of the CLS-D
integrand (all pairwise distances, BELOW_FOLD §1.2).

Every monomial integrates in closed form against e^{-(X1+X2+Y1+Y2)} and the leading (complex-power) density (`integrate`):
  radial  int_0^oo X^c e^{-X} dX = Gamma(c+1)                                  (c = p-exponent/2, complex)
  s       int_{-1}^{1} w^G s^k ds = [k even] Gamma(1/2)Gamma(G+1)/Gamma(G+3/2)   (canonical k in {0,1})
  phi     int_{-pi/2}^{pi/2} cf^M sf^n = [n even] Gamma(1/2)Gamma((M+1)/2)/Gamma(M/2+1)
  alpha   int_{-pi}^{pi} ca^r sa^t u^A v^B = [t even] 2^{A+B+1} sum_j C(r,j)(-2)^j B(A+j+1/2, B+1/2)
"""
import numpy as np
from scipy.special import loggamma, comb

NV = 12
NAMES = ('p1', 'p2', 'q1', 'q2', 's', 'w', 'cf', 'sf', 'ca', 'sa', 'u', 'v')
OFF = 512
_W = np.array([1 << (10*i) for i in range(6)], dtype=np.int64)
TOL_PRUNE = 1e-15


def _combine(E, C):
    """merge duplicate monomials (single mixed-radix int64 key from the actual exponent ranges)."""
    if E.shape[0] == 0:
        return E, C
    E = np.asarray(E)
    mn = E.min(0).astype(np.int64); rg = E.max(0).astype(np.int64) - mn + 1
    if np.sum(np.log2(rg.astype(float))) < 62:
        mult = np.concatenate([[1], np.cumprod(rg[:-1])]).astype(np.int64)
        key = (E.astype(np.int64) - mn) @ mult
        u, inv = np.unique(key, return_inverse=True)
        Cr = np.bincount(inv, weights=C.real, minlength=len(u)); Ci = np.bincount(inv, weights=C.imag, minlength=len(u))
        first = np.full(len(u), -1, np.int64)
        first[inv[::-1]] = np.arange(len(inv))[::-1]
        return E[first], Cr + 1j*Ci
    Eo = E.astype(np.int64) - mn
    W = np.array([1 << (10*i) for i in range(6)], dtype=np.int64)
    h, l = Eo[:, :6] @ W, Eo[:, 6:] @ W
    order = np.lexsort((l, h))
    h, l, E, C = h[order], l[order], E[order], C[order]
    new = np.ones(len(h), bool)
    new[1:] = (h[1:] != h[:-1]) | (l[1:] != l[:-1])
    idx = np.flatnonzero(new)
    return E[idx], np.add.reduceat(C, idx)


MAXPAIRS = 1500000


class Poly:
    __slots__ = ('E', 'C')

    def __init__(self, E=None, C=None):
        if E is None:
            E = np.zeros((0, NV), np.int16); C = np.zeros(0, complex)
        self.E = np.asarray(E, np.int16).reshape(-1, NV); self.C = np.asarray(C, complex).reshape(-1)

    @staticmethod
    def const(c):
        if c == 0:
            return Poly()
        return Poly(np.zeros((1, NV), np.int16), [c])

    @staticmethod
    def mono(c, **ex):
        e = np.zeros((1, NV), np.int16)
        for k, v in ex.items():
            e[0, NAMES.index(k)] = v
        return Poly(e, [c])

    def __len__(self):
        return self.C.size

    def copy(self):
        return Poly(self.E.copy(), self.C.copy())

    def __add__(self, o):
        if not isinstance(o, Poly):
            o = Poly.const(o)
        E, C = _combine(np.vstack([self.E, o.E]), np.concatenate([self.C, o.C]))
        return Poly(E, C).prune()

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.E, -self.C)

    def __sub__(self, o):
        if not isinstance(o, Poly):
            o = Poly.const(o)
        return self + (-o)

    def scale(self, c):
        if c == 0:
            return Poly()
        return Poly(self.E, self.C*c)

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return self.scale(o)
        if len(self) == 0 or len(o) == 0:
            return Poly()
        a, b = (self, o) if len(self) >= len(o) else (o, self)
        step = max(1, MAXPAIRS//len(b))
        parts_E, parts_C = [], []
        for i0 in range(0, len(a), step):
            Ea, Ca = a.E[i0:i0 + step], a.C[i0:i0 + step]
            E = (Ea[:, None, :].astype(np.int32) + b.E[None, :, :]).reshape(-1, NV)
            C = (Ca[:, None]*b.C[None, :]).reshape(-1)
            E, C = _combine(E, C)
            E, C = canon(E, C)
            E, C = _combine(E, C)
            parts_E.append(E.astype(np.int16)); parts_C.append(C)
        if len(parts_E) > 1:
            E, C = _combine(np.vstack(parts_E), np.concatenate(parts_C))
        else:
            E, C = parts_E[0], parts_C[0]
        return Poly(E, C).prune()

    __rmul__ = __mul__

    def prune(self, rel=TOL_PRUNE):
        if self.C.size == 0:
            return self
        mx = np.abs(self.C).max()
        k = np.abs(self.C) > rel*mx
        if k.all():
            return self
        return Poly(self.E[k], self.C[k])

    def mul_mono(self, e, c=1.0):
        """multiply by c * atoms^e (e: length-12 int array)."""
        E = self.E.astype(np.int32) + np.asarray(e, np.int32)[None, :]
        E, C = canon(E, self.C*c)
        E, C = _combine(E, C)
        return Poly(E, C)

    def ev(self, P):
        """evaluate at points: P = (n_pts x 12) atom values."""
        if len(self) == 0:
            return np.zeros(P.shape[0], complex)
        L = np.log(P.astype(complex))                           # atoms may be negative (s, sf, ca, sa)
        return np.exp(L @ self.E.T.astype(float)) @ self.C


def canon(E, C):
    """s^2 -> 1 - w, sf^2 -> 1 - cf^2, sa^2 -> 1 - ca^2 (exponents of s, sf, sa reduced to {0, 1})."""
    E = np.asarray(E, np.int32); C = np.asarray(C, complex)
    for (i, j, sq) in ((4, 5, 'w'), (7, 6, 'c'), (9, 8, 'c')):
        while True:
            m = E[:, i] >= 2
            if not m.any():
                break
            Em, Cm = E[m].copy(), C[m]
            Em[:, i] -= 2
            E1 = Em.copy()                                     # the '1' part
            E2 = Em.copy()
            if sq == 'w':
                E2[:, j] += 1                                  # - w
            else:
                E2[:, j] += 2                                  # - c^2
            E = np.vstack([E[~m], E1, E2]); C = np.concatenate([C[~m], Cm, -Cm])
    return E, C


# ----------------------------------------------------------------------------------- random evaluation points (atoms)
def atom_points(n, seed=0, phys=True):
    rng = np.random.default_rng(seed)
    P = np.empty((n, NV))
    P[:, 0:4] = rng.uniform(0.5, 2.0, (n, 4))
    s = rng.uniform(-0.9, 0.9, n); P[:, 4] = s; P[:, 5] = 1 - s*s
    ph = rng.uniform(-1.4, 1.4, n); P[:, 6] = np.cos(ph); P[:, 7] = np.sin(ph)
    al = rng.uniform(-3.0, 3.0, n); P[:, 8] = np.cos(al); P[:, 9] = np.sin(al); P[:, 10] = 1 - P[:, 8]; P[:, 11] = 1 + P[:, 8]
    return P


_PTS = atom_points(40, 11)
_PTS2 = atom_points(20, 12)
LEAD_ATOMS = (0, 1, 2, 3, 5, 6, 10, 11)                     # atoms allowed in a leading (divisor) monomial


def atomise(P, what=''):
    """recognise a Poly as c * prod atoms^e (atoms LEAD_ATOMS, integer e); returns (c, e) or raises."""
    v = P.ev(_PTS)
    if np.max(np.abs(v)) == 0:
        raise ZeroDivisionError('atomise: zero polynomial ' + what)
    A = np.hstack([np.ones((_PTS.shape[0], 1)), np.log(np.abs(_PTS[:, LEAD_ATOMS]))])
    sol, *_ = np.linalg.lstsq(A, np.log(np.abs(v)), rcond=None)
    e = np.zeros(NV, np.int32)
    e[list(LEAD_ATOMS)] = np.round(sol[1:]).astype(np.int32)
    mono = np.exp(np.log(_PTS2.astype(complex)) @ e.astype(float))
    w = P.ev(_PTS2)/mono
    c = w[0]
    err = np.max(np.abs(w/c - 1))
    if err > 1e-9 or np.max(np.abs(sol[1:] - np.round(sol[1:]))) > 1e-6:
        raise ValueError(f'atomise failed ({what}): exps {np.round(sol[1:], 4)}, const spread {err:.2e}')
    return c, e


# ----------------------------------------------------------------------------------- a-series
class ASer:
    """sum_{k=0..K} a^{d+k} c[k]."""
    __slots__ = ('d', 'c')

    def __init__(self, d, c):
        self.d = int(d); self.c = list(c)

    @property
    def K(self):
        return len(self.c) - 1

    @staticmethod
    def const(x, K):
        p = x if isinstance(x, Poly) else Poly.const(x)
        return ASer(0, [p] + [Poly() for _ in range(K)])

    @staticmethod
    def mono(d, p, K):
        return ASer(d, [p] + [Poly() for _ in range(K)])

    def __add__(self, o):
        if not isinstance(o, ASer):
            o = ASer.const(o, self.K)
        d = min(self.d, o.d); top = min(self.d + self.K, o.d + o.K)
        c = [Poly() for _ in range(top - d + 1)]
        for S in (self, o):
            for k, p in enumerate(S.c):
                i = S.d + k - d
                if i < len(c):
                    c[i] = c[i] + p
        return ASer(d, c).norm()

    __radd__ = __add__

    def __neg__(self):
        return ASer(self.d, [-p for p in self.c])

    def __sub__(self, o):
        if not isinstance(o, ASer):
            o = ASer.const(o, self.K)
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, x):
        return ASer(self.d, [p.scale(x) for p in self.c])

    def __mul__(self, o):
        if not isinstance(o, ASer):
            if isinstance(o, Poly):
                return ASer(self.d, [p*o for p in self.c])
            return self.scale(o)
        K = min(self.K, o.K)
        c = []
        for k in range(K + 1):
            acc = Poly()
            for i in range(k + 1):
                if len(self.c[i]) and len(o.c[k - i]):
                    acc = acc + self.c[i]*o.c[k - i]
            c.append(acc)
        return ASer(self.d + o.d, c).norm()

    __rmul__ = __mul__

    def norm(self):
        """drop leading coefficients that vanish identically (checked numerically; non-canonical forms allowed)."""
        c = self.c; d = self.d
        while len(c) > 1:
            v = c[0].ev(_PTS)
            ref = max(1e-300, max(np.abs(p.C).max(initial=0) for p in c[:1]))
            if len(c[0]) == 0 or np.max(np.abs(v)) < 1e-11*ref:
                c = c[1:]; d += 1
            else:
                break
        return ASer(d, c)

    def trunc(self, K):
        return ASer(self.d, self.c[:K + 1])

    def lead(self, what=''):
        return atomise(self.c[0], what)

    def rel(self, what=''):
        """(c0, e0, eps) with self = a^d * c0 atoms^e0 * (1 + eps), eps an ASer with d = 1 (or zero)."""
        c0, e0 = self.lead(what)
        einv = -e0
        K = self.K
        eps = [Poly()] + [self.c[k].mul_mono(einv, 1.0/c0) for k in range(1, K + 1)]
        return c0, e0, ASer(0, eps)

    def recip(self, what=''):
        c0, e0, eps = self.rel(what)
        s = series_fn(eps, 'inv')
        return ASer(-self.d, [p.mul_mono(-e0, 1.0/c0) for p in s.c])

    def __truediv__(self, o):
        if not isinstance(o, ASer):
            return self.scale(1.0/o)
        return self*o.recip()

    def ev(self, P, a):
        return sum(p.ev(P)*a**(self.d + k) for k, p in enumerate(self.c))


def lmul(c1, c2, K):
    """truncated product of coefficient lists (both starting at a^0)."""
    out = []
    for k in range(K + 1):
        acc = Poly()
        for i in range(k + 1):
            if i < len(c1) and k - i < len(c2) and len(c1[i]) and len(c2[k - i]):
                acc = acc + c1[i]*c2[k - i]
        out.append(acc)
    return out


def series_fn(eps, kind, p=None):
    """f(eps) for eps = O(a) (eps.d = 0, eps.c[0] = 0): 'inv' 1/(1+eps), 'log' log(1+eps), 'pow' (1+eps)^p, 'exp' exp(eps).
    Returns an ASer with d = 0 and K = eps.K (no normalisation)."""
    K = eps.K
    assert eps.d == 0 and len(eps.c[0]) == 0
    if kind == 'inv':
        co = [(-1.0)**n for n in range(K + 1)]
    elif kind == 'log':
        co = [0.0] + [(-1.0)**(n + 1)/n for n in range(1, K + 1)]
    elif kind == 'pow':
        co = [1.0]
        for n in range(1, K + 1):
            co.append(co[-1]*(p - n + 1)/n)
    elif kind == 'exp':
        co = [1.0]
        for n in range(1, K + 1):
            co.append(co[-1]/n)
    else:
        raise ValueError(kind)
    out = [Poly.const(co[0])] + [Poly() for _ in range(K)]
    pw = [Poly.const(1.0)] + [Poly() for _ in range(K)]
    for n in range(1, K + 1):
        pw = lmul(pw, eps.c, K)
        if all(len(q) == 0 for q in pw):
            break
        out = [o + q.scale(co[n]) for o, q in zip(out, pw)]
    return ASer(0, out)


# ----------------------------------------------------------------------------------- integration of monomials
LG_HALF = loggamma(0.5).real


def _lB(x, y):
    return loggamma(x) + loggamma(y) - loggamma(x + y)


def integrate(E, C, E0, E1, nus, chunk=4000):
    """sum_terms C * int [atoms^(E + E0 + i nu E1)] over D (radial Gamma, s, phi, alpha), for each nu.
    E: (n x 12) int; C: (n,) complex; E0, E1: (12,) real (leading-density exponents; E1 multiplies i*nu).
    Returns complex array (len(nus),). Terms vanishing by parity (odd s, sf, sa) are dropped."""
    E = np.asarray(E, float); C = np.asarray(C, complex)
    keep = (np.mod(E[:, 4], 2) == 0) & (np.mod(E[:, 7], 2) == 0) & (np.mod(E[:, 9], 2) == 0)
    assert np.all(E[:, 4] <= 1) and np.all(E[:, 7] <= 1) and np.all(E[:, 9] <= 1), 'not canonical'
    E = E[keep]; C = C[keep]
    out = np.zeros(len(nus), complex)
    nus = np.asarray(nus, float)
    for c0 in range(0, len(C), chunk):
        Ec, Cc = E[c0:c0 + chunk], C[c0:c0 + chunk]
        T = Ec[:, None, :] + E0[None, None, :] + 1j*nus[None, :, None]*E1[None, None, :]    # (n, nnu, 12)
        lg = np.zeros(T.shape[:2], complex)
        for i in range(4):                                               # radial: X^{e/2}
            lg += loggamma(1 + T[..., i]/2)
        G = T[..., 5]                                                    # s (canonical: s^0)
        lg += LG_HALF + loggamma(G + 1) - loggamma(G + 1.5)
        M = T[..., 6]                                                    # phi (canonical: sf^0)
        lg += LG_HALF + loggamma((M + 1)/2) - loggamma(M/2 + 1)
        A, B, r = T[..., 10], T[..., 11], Ec[:, 8].astype(int)          # alpha
        val = np.exp(lg)
        al = np.zeros(T.shape[:2], complex)
        for rr in np.unique(r):
            m = r == rr
            acc = np.zeros((m.sum(), len(nus)), complex)
            for j in range(rr + 1):
                acc += comb(rr, j)*(-2.0)**j*np.exp(_lB(A[m] + j + 0.5, B[m] + 0.5))
            al[m] = acc
        al *= np.exp((A + B + 1)*np.log(2.0))
        out += np.sum(Cc[:, None]*val*al, axis=0)
    return out


def integrate_fast(E, C, E0, E1, nus, chunk=2000):
    """same as integrate(), with per-dimension tables over the distinct exponents (radial x4, s, phi, alpha triple)."""
    E = np.asarray(E, np.int64); C = np.asarray(C, complex)
    if len(C) == 0:
        return np.zeros(len(nus), complex)
    keep = (np.mod(E[:, 4], 2) == 0) & (np.mod(E[:, 7], 2) == 0) & (np.mod(E[:, 9], 2) == 0)
    assert np.all(E[:, 4] <= 1) and np.all(E[:, 7] <= 1) and np.all(E[:, 9] <= 1), 'not canonical'
    E = E[keep]; C = C[keep]
    nus = np.asarray(nus, float)
    z = 1j*nus
    tabs = []; idxs = []
    for i in range(4):
        u, inv = np.unique(E[:, i], return_inverse=True)
        T = (u[:, None] + E0[i] + z[None, :]*E1[i])
        tabs.append(loggamma(1 + T/2)); idxs.append(inv)
    u, inv = np.unique(E[:, 5], return_inverse=True)
    G = u[:, None] + E0[5] + z[None, :]*E1[5]
    tabs.append(LG_HALF + loggamma(G + 1) - loggamma(G + 1.5)); idxs.append(inv)
    u, inv = np.unique(E[:, 6], return_inverse=True)
    Mv = u[:, None] + E0[6] + z[None, :]*E1[6]
    tabs.append(LG_HALF + loggamma((Mv + 1)/2) - loggamma(Mv/2 + 1)); idxs.append(inv)
    trip, inv = np.unique(E[:, [8, 10, 11]], axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    al = np.zeros((len(trip), len(nus)), complex)
    for t, (rr, ea, eb) in enumerate(trip):
        A = ea + E0[10] + z*E1[10]; B = eb + E0[11] + z*E1[11]
        acc = np.zeros(len(nus), complex)
        for j in range(rr + 1):
            acc += comb(rr, j)*(-2.0)**j*np.exp(_lB(A + j + 0.5, B + 0.5))
        al[t] = acc*np.exp((A + B + 1)*np.log(2.0))
    idxs.append(inv)
    out = np.zeros(len(nus), complex)
    for c0 in range(0, len(C), chunk):
        sl = slice(c0, c0 + chunk)
        lg = sum(tabs[d][idxs[d][sl]] for d in range(6))
        out += np.sum(C[sl, None]*np.exp(lg)*al[idxs[6][sl]], axis=0)
    return out
