"""
bj_xs (BELOW_FOLD_SPINJ.md §1.3): vectorised ring arithmetic for bj_ring's closed-form pairing.
An XS is a truncated a-series whose coefficients are Laurent polynomials in the 12 bg_ser atoms AND polynomials in 4 nilpotent-truncated
polarisation variables (eps_s, eps_u, eta_s, eta_u; degree <= 4 on each side).  Stored as one table: rows (power of a, 4 eps exponents,
12 atom exponents) + complex coefficient, with ASer semantics for truncation (nominal leading power d, absolute top power):
  product: d = d1 + d2, top = min(top1 + d2, top2 + d1);  sum: d = min, top = min.
A whole ring multiplication is a few numpy operations (pairs grouped by eps-degree), instead of ~5000 separate ASer products.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S

NC = 17            # columns: 0 power, 1..4 eps (es, eu, fs, fu), 5..16 atoms
DA = 4
CHUNK = 1_500_000


def _combine(R, C):
    if R.shape[0] == 0:
        return R, C
    mn = R.min(0).astype(np.int64); rg = R.max(0).astype(np.int64) - mn + 1
    if np.sum(np.log2(rg.astype(float))) < 62:
        mult = np.concatenate([[1], np.cumprod(rg[:-1])]).astype(np.int64)
        key = (R.astype(np.int64) - mn) @ mult
        u, inv = np.unique(key, return_inverse=True)
    else:
        u, inv = np.unique(R, axis=0, return_inverse=True)
        inv = inv.reshape(-1)
        u = np.arange(len(u))
    Cr = np.bincount(inv, weights=C.real, minlength=len(u)); Ci = np.bincount(inv, weights=C.imag, minlength=len(u))
    first = np.full(len(u), -1, np.int64)
    first[inv[::-1]] = np.arange(len(inv))[::-1]
    return R[first], Cr + 1j*Ci


def _canon(R, C):
    R = R.astype(np.int32)
    for (i, j, sq) in ((5 + 4, 5 + 5, 'w'), (5 + 7, 5 + 6, 'c'), (5 + 9, 5 + 8, 'c')):
        while True:
            m = R[:, i] >= 2
            if not m.any():
                break
            Rm, Cm = R[m].copy(), C[m]
            Rm[:, i] -= 2
            R2 = Rm.copy()
            R2[:, j] += 1 if sq == 'w' else 2
            R = np.vstack([R[~m], Rm, R2]); C = np.concatenate([C[~m], Cm, -Cm])
    return R, C


def _prune(R, C, rel=None):
    """value-weighted prune (bj_ee.value_keep on the atom columns): round-off zeros only."""
    if C.size == 0:
        return R, C
    import bj_ee
    k = bj_ee.value_keep(R[:, 5:], C)
    return R[k], C[k]


class XS:
    __slots__ = ('R', 'C', 'd', 'top')

    def __init__(self, R, C, d, top):
        self.R = np.asarray(R, np.int16).reshape(-1, NC); self.C = np.asarray(C, complex).reshape(-1)
        self.d, self.top = int(d), int(top)

    @staticmethod
    def zero(d=0, top=10**6):
        return XS(np.zeros((0, NC)), np.zeros(0), d, top)

    @staticmethod
    def from_aser(A, eps=(0, 0, 0, 0)):
        rows, cs = [], []
        for k, p in enumerate(A.c):
            if len(p) == 0:
                continue
            r = np.zeros((len(p), NC), np.int16)
            r[:, 0] = A.d + k; r[:, 1:5] = eps; r[:, 5:] = p.E
            rows.append(r); cs.append(p.C)
        if not rows:
            return XS(np.zeros((0, NC)), np.zeros(0), A.d, A.d + A.K)
        return XS(np.vstack(rows), np.concatenate(cs), A.d, A.d + A.K)

    @staticmethod
    def const(x, shift, top_rel, eps=(0, 0, 0, 0)):
        """the number x counted as a^shift (stored coefficient x*Delta^shift is supplied by the caller as x)."""
        r = np.zeros((1, NC), np.int16); r[0, 0] = shift; r[0, 1:5] = eps
        return XS(r, [x], shift, shift + top_rel)

    def __len__(self):
        return self.C.size

    def scale(self, c):
        return XS(self.R, self.C*c, self.d, self.top)

    def shift(self, s, fac=1.0):
        R = self.R.copy(); R[:, 0] += s
        return XS(R, self.C*fac, self.d + s, self.top + s)

    def __add__(self, o):
        d = min(self.d, o.d); top = min(self.top, o.top)
        R = np.vstack([self.R, o.R]); C = np.concatenate([self.C, o.C])
        k = R[:, 0] <= top
        R, C = _combine(R[k], C[k])
        R, C = _prune(R, C)
        return XS(R, C, d, top)

    def __neg__(self):
        return XS(self.R, -self.C, self.d, self.top)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        if not isinstance(o, XS):
            return self.scale(o)
        d = self.d + o.d; top = min(self.top + o.d, o.top + self.d)
        if len(self) == 0 or len(o) == 0:
            return XS.zero(d, top)
        # group by eps degrees (a-side, b-side) to skip forbidden pairs
        def groups(X):
            da = X.R[:, 1].astype(int) + X.R[:, 2]; db = X.R[:, 3].astype(int) + X.R[:, 4]
            g = {}
            for key in set(zip(da.tolist(), db.tolist())):
                g[key] = np.flatnonzero((da == key[0]) & (db == key[1]))
            return g
        g1, g2 = groups(self), groups(o)
        acc_R, acc_C = np.zeros((0, NC), np.int32), np.zeros(0, complex)
        parts_R, parts_C, nrows = [], [], 0

        def flush(acc_R, acc_C, parts_R, parts_C):
            R = np.vstack(parts_R + [acc_R]); C = np.concatenate(parts_C + [acc_C])
            R, C = _combine(R, C)
            R, C = _canon(R, C)
            R, C = _combine(R, C)
            return R.astype(np.int32), C
        for k1, i1 in g1.items():
            for k2, i2 in g2.items():
                if k1[0] + k2[0] > DA or k1[1] + k2[1] > DA:
                    continue
                A, B = self.R[i1].astype(np.int32), o.R[i2].astype(np.int32)
                order = np.argsort(B[:, 0], kind='stable'); B = B[order]; Cb = o.C[i2][order]; pb = B[:, 0]
                Ca_all = self.C[i1]
                step = max(1, CHUNK//max(1, len(B)))
                for s0 in range(0, len(A), step):
                    Aa, Ca = A[s0:s0 + step], Ca_all[s0:s0 + step]
                    nb = np.searchsorted(pb, top - Aa[:, 0], side='right')
                    tot = int(nb.sum())
                    if tot == 0:
                        continue
                    rep = np.repeat(np.arange(len(Aa)), nb)
                    starts = np.concatenate([[0], np.cumsum(nb)[:-1]])
                    jdx = np.arange(tot) - np.repeat(starts, nb)
                    parts_R.append(Aa[rep] + B[jdx]); parts_C.append(Ca[rep]*Cb[jdx])
                    nrows += tot
                    if nrows > CHUNK:
                        acc_R, acc_C = flush(acc_R, acc_C, parts_R, parts_C)
                        parts_R, parts_C, nrows = [], [], 0
        if parts_R:
            acc_R, acc_C = flush(acc_R, acc_C, parts_R, parts_C)
        R, C = _prune(acc_R, acc_C)
        return XS(R, C, d, top)

    __rmul__ = __mul__

    def eps_part(self, key):
        """ASer of the rows with eps exponents == key (nominal d, top)."""
        k = np.all(self.R[:, 1:5] == np.asarray(key)[None, :], axis=1)
        R, C = self.R[k], self.C[k]
        c = []
        for p in range(self.d, self.top + 1):
            m = R[:, 0] == p
            c.append(S.Poly(R[m, 5:], C[m]) if m.any() else S.Poly())
        return S.ASer(self.d, c)

    def keys(self):
        return set(map(tuple, self.R[:, 1:5].tolist()))


# ------------------------------------------------------------------------------ batched multiply-accumulate of bg_ser.ASer
def _pcombine(E, C):
    from bg_ser import _combine as cb, canon
    E, C = cb(E, C)
    E, C = canon(E, C)
    E, C = cb(E, C)
    import bj_ee
    k = bj_ee.value_keep(E, C) if C.size else np.zeros(0, bool)
    return E[k].astype(np.int32), C[k]


def macc(pairs, maxpairs=2_000_000):
    """sum_i c_i A_i B_i for pairs [(A, B, c)] of bg_ser.ASer (c a number), with ASer truncation semantics
    (product: d = dA + dB, K = min(KA, KB); sum: d = min, top = min) and ONE combine per output order (chunked)."""
    pairs = [(A, B, c) for (A, B, c) in pairs if A is not None and B is not None]
    if not pairs:
        return None
    d = min(A.d + B.d for A, B, c in pairs)
    top = min(A.d + B.d + min(A.K, B.K) for A, B, c in pairs)
    acc = {p: (np.zeros((0, S.NV), np.int32), np.zeros(0, complex)) for p in range(d, top + 1)}
    buf = {p: ([], [], 0) for p in range(d, top + 1)}

    def flush(p):
        Es, Cs, n = buf[p]
        if not Es:
            return
        E0, C0 = acc[p]
        acc[p] = _pcombine(np.vstack([E0] + Es), np.concatenate([C0] + Cs))
        buf[p] = ([], [], 0)
    for A, B, c in pairs:
        for i, pa in enumerate(A.c):
            if len(pa) == 0:
                continue
            for j, pb in enumerate(B.c):
                pw = A.d + i + B.d + j
                if pw > top:
                    break
                if len(pb) == 0:
                    continue
                Ea, Ca, Eb, Cb = pa.E.astype(np.int32), pa.C*c, pb.E.astype(np.int32), pb.C
                step = max(1, maxpairs//max(1, len(Cb)))
                for s0 in range(0, len(Ca), step):
                    E = (Ea[s0:s0 + step, None, :] + Eb[None, :, :]).reshape(-1, S.NV)
                    C = (Ca[s0:s0 + step, None]*Cb[None, :]).reshape(-1)
                    Es, Cs, n = buf[pw]
                    Es.append(E); Cs.append(C); n += len(C)
                    buf[pw] = (Es, Cs, n)
                    if n > maxpairs:
                        flush(pw)
    for p in buf:
        flush(p)
    return S.ASer(d, [S.Poly(acc[p][0], acc[p][1]) for p in range(d, top + 1)])


# ------------------------------------------------------------------------------ dense-box multiply-accumulate (no sorting)
DENSE_MAX = int(os.environ.get('BJ_DENSE_MAX', str(12_000_000)))


def macc_dense(pairs, maxpairs=3_000_000, top=None):
    """as macc, but each output order is accumulated with np.bincount on a dense mixed-radix index over the output exponent box when the
    box has <= DENSE_MAX cells (no sort of the ~1e8-1e9 redundant pairs); falls back to macc otherwise."""
    pairs = [(A, B, c) for (A, B, c) in pairs if A is not None and B is not None]
    if not pairs:
        return None
    d = min(A.d + B.d for A, B, c in pairs)
    if top is None:
        top = min(A.d + B.d + min(A.K, B.K) for A, B, c in pairs)
    else:                     # explicit absolute top: each pair contributes the orders it has (shorter pairs = truncated terms)
        d = min(d, top)
    jobs = {p: [] for p in range(d, top + 1)}
    for A, B, c in pairs:
        for i, pa in enumerate(A.c):
            if len(pa) == 0:
                continue
            for j, pb in enumerate(B.c):
                pw = A.d + i + B.d + j
                if pw > top or pw > A.d + B.d + min(A.K, B.K):
                    break
                if len(pb):
                    jobs[pw].append((pa, pb, c))
    out = []
    for p in range(d, top + 1):
        js = jobs[p]
        if not js:
            out.append(S.Poly()); continue
        lo = np.min([pa.E.min(0).astype(np.int64) + pb.E.min(0) for pa, pb, c in js], axis=0)
        hi = np.max([pa.E.max(0).astype(np.int64) + pb.E.max(0) for pa, pb, c in js], axis=0)
        rg = hi - lo + 1
        box = float(np.prod(rg.astype(float)))
        if box > DENSE_MAX:
            r = macc([(S.ASer(0, [pa]), S.ASer(0, [pb]), c) for pa, pb, c in js], maxpairs)
            out.append(r.c[0]); continue
        mult = np.concatenate([[1], np.cumprod(rg[:-1])]).astype(np.int64)
        nbox = int(box)
        accr = np.zeros(nbox); acci = np.zeros(nbox)
        for pa, pb, c in js:
            ka = (pa.E.astype(np.int64) - lo) @ mult; kb = pb.E.astype(np.int64) @ mult
            Ca = pa.C*c; Cb = pb.C
            step = max(1, maxpairs//max(1, len(Cb)))
            for s0 in range(0, len(ka), step):
                key = (ka[s0:s0 + step, None] + kb[None, :]).ravel()
                val = (Ca[s0:s0 + step, None]*Cb[None, :]).ravel()
                accr += np.bincount(key, weights=val.real, minlength=nbox)
                acci += np.bincount(key, weights=val.imag, minlength=nbox)
        nz = np.flatnonzero((accr != 0) | (acci != 0))
        E = np.empty((len(nz), S.NV), np.int64)
        rem = nz.copy()
        for a_ in range(S.NV):
            E[:, a_] = rem % rg[a_] + lo[a_]; rem //= rg[a_]
        E, C = _pcombine(E, accr[nz] + 1j*acci[nz])
        out.append(S.Poly(E, C))
    return S.ASer(d, out)
