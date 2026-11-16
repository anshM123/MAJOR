"""
ctypes wrapper around cblock.c: d = 3 scalar blocks g^{a,b}_{Delta,J} with general external dimensions
(a = -Delta_12/2, b = Delta_34/2), normalised as g ~ z^{(Delta-J)/2} zbar^{(Delta+J)/2} (z << zbar << 1).

   g(r, eta) = (4 r)^Delta sum_{n,j} B_{n,j} r^n P_j(eta),      r = sqrt(rho rhobar), eta = (rho+rhobar)/(2 r).

CBlock(...).series(r, chi, sig, kmax, lmax)[p, k, l] = sum B_{n,j} sig^n (Delta+n)^k r^n P_j^{(l)}(cosh chi).
For rho, rhobar < 0 (the region w, wbar < 0 used by dispersive functionals) the real-analytic block is
   g = (-1)^J (4 r)^Delta sum B_{n,j} (-1)^n r^n P_j(|eta|)   (r > 0, eta = -|eta|),
and  d_eta^l P_j(eta) = (-1)^{j+l} P_j^{(l)}(|eta|).
"""
import os, json, ctypes
import numpy as np, mpmath as mp, sympy as sp

_here = os.path.dirname(os.path.abspath(__file__))
_libs = {}
def _lib(prec):
    if prec not in _libs:
        L = ctypes.CDLL(os.path.join(_here, 'libcblock_%s.so' % prec))
        L.block_size.restype = ctypes.c_long
        L.block_build.restype = ctypes.c_void_p
        L.block_build.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                                  ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        L.block_free.argtypes = [ctypes.c_void_p]
        L.block_get.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_void_p]
        L.block_eval.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                 ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double,
                                 ctypes.c_void_p]
        _libs[prec] = L
    return _libs[prec]

_terms = json.load(open(os.path.join(_here, 'terms.json')))
_a, _b = sp.symbols('a b')
_terms_sym = [(t[0], t[1], 0 if t[2] == 'I' else 1, [sp.sympify(c) for c in t[3:]]) for t in _terms]

def term_table(a, b):
    tk = np.array([t[0] for t in _terms_sym], dtype=np.int32)
    ti = np.array([t[1] for t in _terms_sym], dtype=np.int32)
    tt = np.array([t[2] for t in _terms_sym], dtype=np.int32)
    ar, br = sp.nsimplify(a), sp.nsimplify(b)
    tc = np.array([[float(c.subs({_a: ar, _b: br})) for c in t[3]] for t in _terms_sym], dtype=np.float64)
    return tk, ti, tt, np.ascontiguousarray(tc)

def to2(x):
    """mp/float array -> (hi, lo) double pairs"""
    x = np.atleast_1d(x)
    out = np.zeros((x.size, 2))
    for i, v in enumerate(x.flat):
        v = mp.mpf(v)
        hi = float(v); out[i, 0] = hi; out[i, 1] = float(v - mp.mpf(hi))
    return out

class CBlock:
    def __init__(self, Delta, J, a, b, N, prec='d'):
        """Delta, a, b must be exactly representable (doubles; a, b dyadic rationals)."""
        self.Delta, self.J, self.a, self.b, self.N, self.prec = float(Delta), int(J), float(a), float(b), int(N), prec
        assert float(mp.mpf(self.Delta)) == self.Delta
        self.L = _lib(prec)
        tk, ti, tt, tc = term_table(a, b)
        self.off = np.zeros(self.J + self.N + 2, dtype=np.int64)
        self.D2 = np.array([self.Delta, 0.0])
        self.ptr = self.L.block_build(self.D2.ctypes.data, self.J, self.N, len(tk), tk.ctypes.data, ti.ctypes.data, tt.ctypes.data,
                                      tc.ctypes.data, self.off.ctypes.data)

    def __del__(self):
        try:
            self.L.block_free(self.ptr)
        except Exception:
            pass

    def coeff(self, n, j):
        if abs(j - self.J) > n or (n + j - self.J) % 2 or j < 0:
            return mp.mpf(0)
        o = np.zeros(2)
        self.L.block_get(self.ptr, int(self.off[j] + (n - abs(j - self.J))//2), o.ctypes.data)
        return mp.mpf(o[0]) + mp.mpf(o[1])

    def series(self, r, chi, sig=+1, kmax=0, lmax=0, cut=None, N=None):
        """returns array [npts, kmax+1, lmax+1]; float64 for prec 'd' (inputs floats), object(mpf) for 'q'."""
        N = self.N if N is None else min(N, self.N)
        if self.prec == 'd':
            r = np.ascontiguousarray(np.atleast_1d(r), dtype=float); chi = np.ascontiguousarray(np.atleast_1d(chi), dtype=float)
            r2 = np.zeros((r.size, 2)); r2[:, 0] = r
            c2 = np.zeros((r.size, 2)); c2[:, 0] = chi
            cut = 42.0 if cut is None else cut
        else:
            r2 = to2(r); c2 = to2(chi)
            cut = 85.0 if cut is None else cut
        npts = r2.shape[0]
        out = np.zeros((npts, kmax+1, lmax+1, 2))
        self.L.block_eval(self.ptr, self.off.ctypes.data, self.D2.ctypes.data, self.J, N, npts, r2.ctypes.data, c2.ctypes.data,
                          int(sig), int(kmax), int(lmax), float(cut), out.ctypes.data)
        if self.prec == 'd':
            return out[..., 0] + out[..., 1]
        res = np.empty((npts, kmax+1, lmax+1), dtype=object)
        for idx in np.ndindex(npts, kmax+1, lmax+1):
            res[idx] = mp.mpf(out[idx + (0,)]) + mp.mpf(out[idx + (1,)])
        return res

    def g(self, rho, rhob, **kw):
        """block value for real rho, rhob of equal sign (float inputs)."""
        rho, rhob = np.atleast_1d(rho).astype(float), np.atleast_1d(rhob).astype(float)
        A, Bq = np.abs(rho), np.abs(rhob)
        r = np.sqrt(A*Bq); chi = np.abs(0.5*np.log(A/Bq))
        neg = rho < 0
        assert np.all(neg) or not np.any(neg)
        sig = -1 if np.all(neg) else 1
        S = self.series(r, chi, sig=sig, **kw)[:, 0, 0]
        if self.prec == 'q':
            return np.array([(-1)**self.J*(4*mp.mpf(rr))**mp.mpf(self.Delta)*s if sig < 0 else (4*mp.mpf(rr))**mp.mpf(self.Delta)*s
                             for rr, s in zip(r, S)], dtype=object)
        return ((-1.0)**self.J if sig < 0 else 1.0)*np.exp(self.Delta*np.log(4*r))*S


def block_mp(Delta, J, a, b, N, rho, rhob):
    """Independent slow mpmath implementation (explicit dictionaries, exact Legendre algebra done term by term)."""
    Delta, a, b = mp.mpf(Delta), mp.mpf(a), mp.mpf(b)
    C = Delta*(Delta-3) + J*(J+1)
    terms = [(t[0], t[1], t[2], [mp.mpmathify(str(sp.nsimplify(c.subs({_a: sp.nsimplify(a), _b: sp.nsimplify(b)})).evalf(40))) for c in t[3]])
             for t in _terms_sym]
    def eta_mul(v):
        o = {}
        for j, x in v.items():
            o[j+1] = o.get(j+1, 0) + x*mp.mpf(j+1)/(2*j+1)
            if j >= 1: o[j-1] = o.get(j-1, 0) + x*mp.mpf(j)/(2*j+1)
        return o
    def Dop(v):
        o = {}
        for j, x in v.items():
            c = mp.mpf(j*(j+1))/(2*j+1)
            if j >= 1: o[j-1] = o.get(j-1, 0) + c*x
            o[j+1] = o.get(j+1, 0) - c*x
        return o
    B = {0: {J: mp.factorial(J)/mp.rf(mp.mpf(1)/2, J)}}
    for n in range(1, N+1):
        rhs = {}
        for (k, i, typ, c) in terms:
            if k == 0 or k > n: continue
            E = Delta + n - k
            v = {j: (c[0] + c[1]*E + c[2]*E*E + c[3]*j*(j+1) + c[4]*C)*x for j, x in B[n-k].items()}
            if typ == 1: v = Dop(v)
            for _ in range(i): v = eta_mul(v)
            for j, x in v.items(): rhs[j] = rhs.get(j, 0) + x
        E = Delta + n
        B[n] = {j: -x/((E*(E-3) + j*(j+1) - C)/2) for j, x in rhs.items() if abs(j-J) <= n and (n+j-J) % 2 == 0}
    rho, rhob = mp.mpmathify(rho), mp.mpmathify(rhob)
    r = mp.sqrt(rho*rhob); eta = (rho+rhob)/(2*r)
    s = mp.mpf(0)
    for n, lev in B.items():
        s += r**n*sum(x*mp.legendre(j, eta) for j, x in lev.items())
    return (4*r)**Delta*s
