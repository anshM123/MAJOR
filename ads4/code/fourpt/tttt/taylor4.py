"""
Numerical Taylor-mode (jet) application of the CPPR differential operators of fourpt/embed2 (D_11, D_12, D_22, D_21, H_12 at the vertex (1,2);
the same at (3,4)) to the scalar seed  F = Kt(s) gt(P, Q)  at a numerical point, without symbolic expansion.

Representation: the operators act on functions of the 24 dot-product symbols (6 s_ij = 2 P_i.P_j, 18 z's = Z_i.P_j, Z_i.Z_j), treated as
independent (embed2: the operators preserve the constraints).  A state is a polynomial in the 18 z's (exact exponents, dict keys) whose
coefficients are truncated Taylor series in the 6 deviations ds_ij = s_ij - s0_ij (dense arrays over multi-indices of total degree <= M).
Every D is first order:  V.d/dW  = sum over symbols sigma containing W of  (d sigma/d W . V) d/d sigma ;  d/d s_ij lowers the Taylor degree by one,
so M = (number of D's in the word pair) <= 4l suffices for the value at ds = 0.

Initial series (embed2.word4 conventions):  Kt = (s24/s14)^{(d1-d2)/2} (s24/s23)^{(d3-d4)/2} s12^{-(d1+d2)/2} s34^{-(d3+d4)/2},
P = (log u + log v)/2, Q = (log v - log u)/2,  u = s12 s34/(s13 s24),  v = s14 s23/(s13 s24);  gt from its (P,Q)-jet G[i][j] = d_P^i d_Q^j gt.
Result:  value = D_wL^{(12)} D_wR^{(34)} [Kt gt]  =  Kt * Phi  of embed2.word4  (z-polynomial dict).
Number type: float64 (numpy), or any object type (mpmath mpf / flint arb) via dtype=object.
"""
import itertools, functools
import numpy as np

SN = ['s12', 's13', 's14', 's23', 's24', 's34']
ZN = ['z12', 'z13', 'z14', 'z21', 'z23', 'z24', 'z31', 'z32', 'z34', 'z41', 'z42', 'z43', 'z1z2', 'z1z3', 'z1z4', 'z2z3', 'z2z4', 'z3z4']

# symbol table: name -> (vector A, vector B)
PAIR = {}
for n in SN: PAIR[n] = ('P' + n[1], 'P' + n[2])
for n in ZN:
    PAIR[n] = ('Z' + n[1], 'Z' + n[3]) if n.count('z') == 2 else ('Z' + n[1], 'P' + n[2])
DOT = {}
for n, (a, b) in PAIR.items():
    DOT[(a, b)] = n; DOT[(b, a)] = n
SIDX = {n: i for i, n in enumerate(SN)}; ZIDX = {n: i for i, n in enumerate(ZN)}

def dot(A, B):
    """(symbol name, factor) of A.B, or None if it vanishes (P_i^2 = Z_i^2 = Z_i.P_i = 0)."""
    n = DOT.get((A, B))
    if n is None: return None
    return (n, 0.5) if n in SIDX else (n, 1.0)

# ------------------------------------------------------------------ multi-index tables
@functools.lru_cache(None)
def tables(M, nv=6):
    mons = [a for d in range(M + 1) for a in itertools.product(range(d + 1), repeat=nv) if sum(a) == d]
    mons = sorted(set(mons), key=lambda a: (sum(a), tuple(-x for x in a)))
    idx = {a: i for i, a in enumerate(mons)}
    n = len(mons)
    dn = []; dfac = []; up = []
    for k in range(nv):
        g = np.full(n, -1, dtype=np.int64); f = np.zeros(n)
        u = np.full(n, -1, dtype=np.int64)
        for a, i in idx.items():
            b = list(a); b[k] += 1; b = tuple(b)
            if b in idx: g[i] = idx[b]; f[i] = a[k] + 1
            if a[k] >= 1:
                c = list(a); c[k] -= 1; u[i] = idx[tuple(c)]
        dn.append(g); dfac.append(f); up.append(u)
    # product table
    I, Jx, K = [], [], []
    for a, i in idx.items():
        for b, j in idx.items():
            if sum(a) + sum(b) <= M:
                I.append(i); Jx.append(j); K.append(idx[tuple(x + y for x, y in zip(a, b))])
    nb = [sum(1 for a in mons if sum(a) <= b) for b in range(M + 1)]
    return dict(mons=mons, idx=idx, n=n, dn=dn, dfac=dfac, up=up, I=np.array(I), J=np.array(Jx), K=np.array(K), nb=nb)

class Series:
    """truncated Taylor series in the 6 ds variables; data shape [n] + extra (vectorised points)."""
    def __init__(self, M, data):
        self.M, self.d = M, data

def _zeros(T, extra, dtype):
    z = np.zeros((T['n'],) + extra, dtype=dtype)
    if dtype == object: z[...] = 0
    return z

def smul(T, a, b):
    out = np.zeros_like(a)
    if a.dtype == object:
        out[...] = 0
        prod = a[T['I']]*b[T['J']]
        for k, v in zip(T['K'], prod): out[k] = out[k] + v
    elif a.dtype == np.float64:
        if 'Kmat' not in T:
            import scipy.sparse as sps
            T['Kmat'] = sps.csr_matrix((np.ones(len(T['K'])), (T['K'], np.arange(len(T['K'])))), shape=(T['n'], len(T['K'])))
        out = T['Kmat'] @ (a[T['I']]*b[T['J']])
    else:
        np.add.at(out, T['K'], a[T['I']]*b[T['J']])
    return out

def sderiv(T, a, k):
    """derivative; a has length nb[b] (degree <= b), result length nb[b-1]"""
    b = T['nb'].index(a.shape[0])
    n1 = T['nb'][b - 1] if b > 0 else 0
    g, f = T['dn'][k][:n1], T['dfac'][k][:n1]
    out = np.zeros((n1,) + a.shape[1:], dtype=a.dtype)
    if a.dtype == object: out[...] = 0
    m = g >= 0
    fe = f[m].reshape((-1,) + (1,)*(a.ndim - 1))
    out[m] = a[g[m]]*fe if a.dtype != object else a[g[m]]*np.array(f[m], dtype=object).reshape(fe.shape)
    return out

def smulvar(T, a, k, s0):
    """(s0 + ds_k) * a"""
    u = T['up'][k][:a.shape[0]]
    out = a*s0
    m = u >= 0
    out[m] = out[m] + a[u[m]]
    return out

# ------------------------------------------------------------------ states: dict zexp(tuple of 18) -> series array
def st_add(A, B, c=1):
    """sum of two states; series of different length are truncated to the shorter one (a D-operator has at most one P-derivative,
    so after it at most (remaining budget) derivatives act: the Z-derivative terms may be truncated like the P-derivative ones)"""
    la = min((v.shape[0] for v in A.values()), default=None); lb = min((v.shape[0] for v in B.values()), default=None)
    L = min(x for x in (la, lb) if x is not None) if (la is not None or lb is not None) else None
    out = {k: v[:L] for k, v in A.items()}
    for k, v in B.items():
        v = v[:L]
        out[k] = out[k] + c*v if k in out else c*v
    return out

def st_scale(A, c):
    return {k: c*v for k, v in A.items()}

class Engine4:
    def __init__(self, M, s0, dtype=float):
        """s0: dict s-name -> value (scalar or array over points)"""
        self.T = tables(M); self.M = M; self.dtype = dtype
        self.s0 = [s0[n] for n in SN]

    def mul_sym(self, A, name, fac):
        out = {}
        if name in SIDX:
            k = SIDX[name]
            for key, v in A.items():
                out[key] = smulvar(self.T, v, k, self.s0[k])*fac
        else:
            m = ZIDX[name]
            for key, v in A.items():
                kk = list(key); kk[m] += 1; kk = tuple(kk)
                out[kk] = out[kk] + v*fac if kk in out else v*fac
        return out

    def mul_dot(self, A, X, Y):
        d = dot(X, Y)
        if d is None: return {}
        return self.mul_sym(A, d[0], d[1])

    def d_sym(self, A, name):
        out = {}
        if name in SIDX:
            k = SIDX[name]
            for key, v in A.items():
                out[key] = sderiv(self.T, v, k)
        else:
            m = ZIDX[name]
            for key, v in A.items():
                if key[m] > 0:
                    kk = list(key); kk[m] -= 1; kk = tuple(kk)
                    w = v*key[m]
                    out[kk] = out[kk] + w if kk in out else w
        return out

    def dirder(self, A, V, W):
        """V.d/dW acting on the state"""
        out = {}
        for name, (a, b) in PAIR.items():
            if W not in (a, b): continue
            other = b if a == W else a
            d = dot(V, other)
            if d is None: continue
            fac = 2.0 if name in SIDX else 1.0          # s = 2 P.P
            out = st_add(out, self.mul_sym(self.d_sym(A, name), d[0], d[1]*fac))
        return out

    def op(self, A, kind, i, j):
        Pi, Pj, Zi, Zj = 'P%d' % i, 'P%d' % j, 'Z%d' % i, 'Z%d' % j
        if kind == 'H':
            return st_add(self.mul_dot(self.mul_dot(A, Zi, Zj), Pi, Pj), self.mul_dot(self.mul_dot(A, Zi, Pj), Zj, Pi), -1)
        if kind in ('A', 'C'):
            if kind == 'C': Pi, Pj, Zi, Zj = Pj, Pi, Zj, Zi
            r = self.mul_dot(self.dirder(A, Zi, Pj), Pi, Pj)
            r = st_add(r, self.mul_dot(self.dirder(A, Pi, Pj), Zi, Pj), -1)
            r = st_add(r, self.mul_dot(self.dirder(A, Pi, Zj), Zi, Zj), -1)
            r = st_add(r, self.mul_dot(self.dirder(A, Zi, Zj), Pi, Zj))
            return r
        if kind in ('B', 'E'):
            if kind == 'E': Pi, Pj, Zi, Zj = Pj, Pi, Zj, Zi
            r = self.mul_dot(self.dirder(A, Zi, Pi), Pi, Pj)
            r = st_add(r, self.mul_dot(self.dirder(A, Pi, Pi), Zi, Pj), -1)
            r = st_add(r, self.mul_dot(self.dirder(A, Zi, Zi), Zi, Pj))
            return r
        raise ValueError(kind)

    def apply_words(self, A, wL, wR):
        for kind in reversed(wL):
            A = self.op(A, kind, 1, 2)
        for kind in reversed(wR):
            A = self.op(A, kind, 3, 4)
        return {k: v[0] for k, v in A.items() if v.shape[0]}             # value at ds = 0

# ------------------------------------------------------------------ initial series
def log_series(T, k, s0, extra, dtype, one=1.0):
    """log(s0 + ds_k) - log(s0)"""
    a = _zeros(T, extra, dtype)
    M = max(sum(m) for m in T['mons'])
    for n in range(1, M + 1):
        e = [0]*6; e[k] = n
        a[T['idx'][tuple(e)]] = (one*(-1)**(n + 1)/n)/s0**n
    return a

def exp_series(T, x, one=1.0):
    """exp(x) for a series with zero constant term"""
    M = max(sum(m) for m in T['mons'])
    out = _zeros(T, x.shape[1:], x.dtype); out[0] = out[0] + one
    term = out.copy()
    for n in range(1, M + 1):
        term = smul(T, term, x)/n
        out = out + term
    return out

def seed(M, s0, d, G, dtype=float, one=1.0):
    """state {0: series of Kt*gt}; d = (d1,d2,d3,d4) seed dimensions; G[i][j] = d_P^i d_Q^j gt at the point (arrays over points ok)."""
    T = tables(M)
    sv = [s0[n] for n in SN]
    extra = np.shape(sv[0])
    L = [log_series(T, k, sv[k], extra, dtype, one) for k in range(6)]
    l12, l13, l14, l23, l24, l34 = L
    U = l12 + l34 - l13 - l24; V = l14 + l23 - l13 - l24
    dP = (U + V)/2; dQ = (V - U)/2
    d1, d2, d3, d4 = d
    lK = (d1 - d2)/2*(l24 - l14) + (d3 - d4)/2*(l24 - l23) - (d1 + d2)/2*l12 - (d3 + d4)/2*l34
    K0 = (sv[4]/sv[2])**((d1 - d2)/2)*(sv[4]/sv[3])**((d3 - d4)/2)/(sv[0]**((d1 + d2)/2)*sv[5]**((d3 + d4)/2))
    Kser = exp_series(T, lK, one)*K0
    # gt series: sum_{i+j<=M} G_ij dP^i dQ^j / (i! j!)
    powP = [None]*(M + 1); powQ = [None]*(M + 1)
    e0 = _zeros(T, extra, dtype); e0[0] = e0[0] + one
    powP[0] = e0; powQ[0] = e0
    for i in range(1, M + 1):
        powP[i] = smul(T, powP[i - 1], dP); powQ[i] = smul(T, powQ[i - 1], dQ)
    g = _zeros(T, extra, dtype)
    fact = [1]
    for i in range(1, M + 1): fact.append(fact[-1]*i)
    for i in range(M + 1):
        for j in range(M + 1 - i):
            g = g + smul(T, powP[i], powQ[j])*(G[i][j]/(fact[i]*fact[j]))
    return {(0,)*18: smul(T, Kser, g)}

def word_budget(wL, wR):
    return sum(1 for c in wL + wR if c != 'H')

# ------------------------------------------------------------------ (P,Q)-jets of the scalar blocks to any order (generalises e2e/spinblock.pq_jets)
@functools.lru_cache(None)
def _acosh_derivs(M):
    """k-th derivative of acosh = P_k(c)/sh^{2k-1}, sh = sqrt(c^2-1) passed separately (no cancellation near c = 1)"""
    import sympy as sp
    c, t = sp.symbols('c t', positive=True)
    out = []
    for k in range(1, M + 1):
        Pk = sp.expand(sp.simplify(sp.diff(sp.acosh(c), c, k)*(c**2 - 1)**sp.Rational(2*k - 1, 2)))
        out.append(sp.lambdify((c, t), Pk/t**(2*k - 1), 'mpmath'))
    return out

def _pmul2(A, B, M):
    out = np.zeros_like(A)
    if A.dtype == object: out[...] = 0
    for i in range(M + 1):
        for j in range(M + 1 - i):
            a = A[..., i, j]
            for k in range(M + 1 - i - j):
                for m in range(M + 1 - i - j - k):
                    out[..., i + k, j + m] = out[..., i + k, j + m] + a*B[..., k, m]
    return out

def pq_jets_M(S, r, A0, E0, M, vpow=0.0, one=1.0):
    """single point.  S[k, l] = sum Bt (Delta+n)^k r^n P_j^{(l)}(cosh chi), k + l <= M.  Returns G[i, j] = d_P^i d_Q^j [gt/(4r)^Delta]
    (same conventions as e2e/spinblock.pq_jets, any order).  Works with floats or mpmath numbers."""
    import math
    dt = object if not isinstance(one, float) else float
    Z = lambda: (np.zeros((M + 1, M + 1), dtype=dt) if dt == float else np.full((M + 1, M + 1), 0*one, dtype=object))
    dA = Z(); dE = Z()
    for i in range(M + 1):
        for j in range(M + 1 - i):
            if i + j > 0:
                dA[i, j] = A0*(-one/2)**i*(one/2)**j/(math.factorial(i)*math.factorial(j))
        if i > 0:
            dE[0, i] = E0/math.factorial(i)
    dc = dA + dE; dx = dE - dA
    c0 = (r + one/r)/2; sh0 = (one/r - r)/2
    fk = [f(c0, sh0) for f in _acosh_derivs(M)]
    if dt == float: fk = [float(v) for v in fk]
    dL = Z(); p = Z(); p[0, 0] = one
    for k in range(1, M + 1):
        p = _pmul2(p, dc, M)
        dL = dL - p*(fk[k - 1]/math.factorial(k))
    onep = Z(); onep[0, 0] = one
    powL = [onep]; powX = [onep]
    for _ in range(M):
        powL.append(_pmul2(powL[-1], dL, M)); powX.append(_pmul2(powX[-1], dx, M))
    G = Z()
    for k in range(M + 1):
        for l in range(M + 1 - k):
            G = G + _pmul2(powL[k], powX[l], M)*(S[k, l]/(math.factorial(k)*math.factorial(l)))
    if vpow != 0:
        ex = Z()
        for i in range(M + 1):
            for j in range(M + 1 - i):
                ex[i, j] = (one*vpow)**(i + j)/(math.factorial(i)*math.factorial(j))
        G = _pmul2(G, ex, M)*((E0/A0)**(2*vpow))
    for i in range(M + 1):
        for j in range(M + 1 - i):
            G[i, j] = G[i, j]*(math.factorial(i)*math.factorial(j))
    return G
