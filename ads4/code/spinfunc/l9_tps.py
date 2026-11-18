"""
l9_tps (M6 round 5, doc 23.1): the NEXT Laplace order K3 (relative O(1/Delta^3), "eps^6") for the l1/l6 crystal engine.

Needs the build-side TPS at DEG >= 8 (env L1_DEG=8 before importing l1_tps): (lamF)_8, G_6, q_6.
Expectations are EXACT Gaussian moments of products of polynomials whose total degree reaches 18 (h3^6), so they are taken in a separate
product space of degree 9 (every factor product formed below has true degree <= 9) and paired by the Hermite inner product
    E[a b] = sum_alpha alpha! a^_alpha b^_alpha,      a^ = coefficients of a in the multivariate probabilists' Hermite basis He_alpha,
which is exact for polynomials of any degree (no truncation).  This replaces l1_tps.expect2's O(N^2) monomial loop.

expand_pref6(lamF, G, q) -> (S0, S2, S4, S6):  E[q e^{A}] by eps-orders, A = sum_k B_k eps^k, B_k = (lamF)_{k+2} + G_k, q_j ~ eps^j;
exp(sum B_k eps^k) = sum_n E_n eps^n (complete Bell polynomials), S_{2n} = sum_j E[q_j E_{2n-j}].
S0, S2, S4 coincide with l1_tps.expand_pref (checked in _selftest).
"""
import os, sys, math, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sp
import l1_tps as T

assert T.DEG >= 8, "l9_tps needs L1_DEG=8 (set the environment variable before importing l1_tps)"
NV = T.NV
PD = 9                                                     # product-space degree


class Space:
    def __init__(self, nv, deg):
        mons = []
        for d_ in range(deg + 1):
            for c in itertools.combinations_with_replacement(range(nv), d_):
                e = [0]*nv
                for i in c:
                    e[i] += 1
                mons.append(tuple(e))
        self.N = len(mons); self.EXP = np.array(mons, dtype=np.int64); self.DEGS = self.EXP.sum(1); self.deg = deg
        B = deg + 1
        self.CODE = (self.EXP*(B**np.arange(nv))).sum(1)
        o = np.argsort(self.CODE); self._o, self._s = o, self.CODE[o]
        I, Jx = [], []
        for a in range(self.N):
            bs = np.where(self.DEGS <= deg - self.DEGS[a])[0]
            I.append(np.full(len(bs), a)); Jx.append(bs)
        self.I = np.concatenate(I); self.J = np.concatenate(Jx); self.K = self.lookup(self.CODE[self.I] + self.CODE[self.J])
        # Hermite conversion: x^n = sum_k n!/(k!(n-2k)! 2^k) He_{n-2k}  (tensor product over variables)
        rows, cols, vals = [], [], []
        for a, e in enumerate(self.EXP):
            parts = [[(e_v - 2*k, math.factorial(e_v)/(math.factorial(k)*math.factorial(e_v - 2*k)*2**k)) for k in range(e_v//2 + 1)]
                     for e_v in e]
            for combo in itertools.product(*parts):
                ee = np.array([c[0] for c in combo]); v = np.prod([c[1] for c in combo])
                rows.append(self.lookup(np.array([(ee*(B**np.arange(nv))).sum()]))[0]); cols.append(a); vals.append(v)
        self.H = sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N))
        self.fact = np.array([np.prod([math.factorial(v) for v in e]) for e in self.EXP], float)

    def lookup(self, codes):
        return self._o[np.searchsorted(self._s, codes)]

    def mul(self, a, b):
        w = a[self.I]*b[self.J]
        return np.bincount(self.K, w.real, self.N) + 1j*np.bincount(self.K, w.imag, self.N)

    def herm(self, a):
        return self.H @ a

    def E(self, a, b=None):
        """exact E[a b] (b = 1 if None), xi ~ N(0, 1_NV)."""
        ha = self.herm(a)
        if b is None:
            return ha[0]
        return np.sum(self.fact*ha*self.herm(b))


_SP = None


def space():
    global _SP
    if _SP is None:
        _SP = Space(NV, PD)
    return _SP


_LIFT = None


def lift(a):
    """T-space TPS (DEG = T.DEG) -> product space (degree PD); coefficients of degree > PD dropped (none used)."""
    global _LIFT
    S = space()
    if _LIFT is None:
        keep = np.where(T.DEGS <= PD)[0]
        codes = (T.EXP[keep]*((PD + 1)**np.arange(NV))).sum(1)
        _LIFT = (keep, S.lookup(codes))
    out = np.zeros(S.N, complex); out[_LIFT[1]] = a[_LIFT[0]]
    return out


def expand_pref6(lamF, G, q=None):
    S = space(); m = S.mul; E = S.E
    h = [lift(T.homog(lamF, k)) for k in range(9)]
    g = [lift(T.homog(G, k)) for k in range(7)]
    qq = [lift(T.homog(q if q is not None else T.const(1.0), k)) for k in range(7)]
    B = [None] + [h[k + 2] + g[k] for k in range(1, 7)]           # B1..B6
    B1, B2, B3, B4, B5, B6 = B[1:]
    B11 = m(B1, B1); B111 = m(B11, B1); B12 = m(B1, B2); B22 = m(B2, B2); B13 = m(B1, B3)
    q0, q1, q2, q3, q4, q5, q6 = qq
    c0 = q0[0]
    # E_n (only even total orders survive)
    S0 = c0
    S2 = c0*(E(B2) + 0.5*E(B11)) + E(q1, B1) + E(q2)
    e4 = E(B4) + E(B1, B3) + 0.5*E(B2, B2) + 0.5*E(B11, B2) + E(B11, B11)/24.0
    q1B1 = m(q1, B1)
    S4 = (c0*e4 + E(q1, B3) + E(q1B1, B2) + E(q1B1, B11)/6.0
          + E(q2, B2) + 0.5*E(q2, B11) + E(q3, B1) + E(q4))
    e6 = (E(B6) + E(B1, B5) + E(B2, B4) + 0.5*E(B3, B3) + 0.5*E(B11, B4) + E(B12, B3) + E(B22, B2)/6.0
          + E(B11, B13)/6.0 + 0.25*E(B11, B22) + E(B111, B12)/24.0 + E(B111, B111)/720.0)
    q1B11 = m(q1B1, B1)
    e5q1 = (E(q1, B5) + E(q1B1, B4) + E(m(q1, B2), B3) + 0.5*E(q1B11, B3) + 0.5*E(q1B1, B22)
            + E(q1B11, B12)/6.0 + E(q1B11, B111)/120.0)
    q2B1 = m(q2, B1)
    e4q2 = E(q2, B4) + E(q2B1, B3) + 0.5*E(q2, B22) + 0.5*E(m(q2B1, B1), B2) + E(q2B1, B111)/24.0
    e3q3 = E(q3, B3) + E(m(q3, B1), B2) + E(q3, B111)/6.0
    e2q4 = E(q4, B2) + 0.5*E(q4, B11)
    S6 = c0*e6 + e5q1 + e4q2 + e3q3 + e2q4 + E(q5, B1) + E(q6)
    return S0, S2, S4, S6


def _selftest():
    """(i) S0..S4 vs l1_tps.expand_pref on random TPS; (ii) exactness of E on a 1D Gaussian check: E[exp(t xi_0)] series."""
    rng = np.random.default_rng(1)
    lamF = T.const(0.3); G = T.const(0.1); q = T.const(1.2)
    for k in range(1, T.DEG + 1):
        msk = T.DEGS == k
        lamF[msk] = (rng.normal(size=msk.sum()) + 1j*rng.normal(size=msk.sum()))*0.1 if k >= 3 else 0
        G[msk] = (rng.normal(size=msk.sum()) + 1j*rng.normal(size=msk.sum()))*0.1 if k <= 6 else 0
        q[msk] = (rng.normal(size=msk.sum()) + 1j*rng.normal(size=msk.sum()))*0.1 if k <= 6 else 0
    a = expand_pref6(lamF, G, q); b = T.expand_pref(lamF, G, q)
    err = max(abs(a[i] - b[i])/max(abs(b[i]), 1e-300) for i in range(3))
    # exponent G = t xi_0 (eps-order 1): E[e^{t xi}] = e^{t^2/2} -> S2 = t^2/2, S4 = t^4/8, S6 = t^6/48
    t = 0.7; G1 = T.const(0.0); G1[1] = t
    s = expand_pref6(T.const(0.0), G1, None)
    err2 = max(abs(s[1] - t**2/2), abs(s[2] - t**4/8), abs(s[3] - t**6/48))
    # cubic exponent h3 = c xi^3 (eps-order 1): S2 = E[c xi^4 ... ] : E[(c xi^3)^2/2] = 15 c^2/2; S6 term check vs direct formula
    c = 0.05; L3 = T.const(0.0); L3[np.where([tuple(e) == (3, 0, 0, 0, 0, 0, 0) for e in T.EXP])[0][0]] = c
    s3 = expand_pref6(L3, T.const(0.0), None)
    df = lambda n: float(np.prod(np.arange(n - 1, 0, -2))) if n > 0 else 1.0
    ex = [1, c**2*df(6)/2, c**4*df(12)/24, c**6*df(18)/720]
    err3 = max(abs(s3[i] - ex[i])/abs(ex[i]) for i in range(4))
    return err, err2, err3


if __name__ == '__main__':
    print("l9_tps self-test: |S0..S4 - l1_tps.expand_pref| rel = %.1e;  E[e^{t xi}] orders err = %.1e;  E[e^{c xi^3}] orders rel err = %.1e"
          % _selftest())
