"""
d = 3 scalar blocks near the crossed lightcone rhobar -> 1 (the corner u = sqrt(rhobar_u/rho_u) -> 0 of the Regge-moment
integration region, SPINNING_E2E.md sect. 3.2), by a lightcone reorganisation of the radial series plus an asymptotic tail.

Notation (as in ads4/code/e2e/blocks.py; the reflected block g^{c1,c2} is evaluated at positive radial coordinates):
   Gt(r, eta) = sum_{n,j} B_{n,j} r^n P_j(eta),    g = (4r)^Delta Gt,
   rho_s = r e^{-chi}, rho_b = r e^{chi}  (eta = cosh chi);   Gt = (rho_s rho_b)^{-J/2} sum_{i,q} C_{i,q} rho_s^i rho_b^q   (lcconv.c).

At fixed i, H_i(x) = sum_q C_{i,q} x^q is a function of one variable whose only singularities on |x| <= 1 are
   x = 1  (zbar = 1):  analytic + (1-x)^{-2 sigma} * analytic,     sigma = c1 + c2   (local exponents 0, -sigma of the
                                                                 Casimir equation at zbar = 1, in rho: (1-zbar) = ((1-x)/(1+x))^2)
   x = -1 (zbar = oo): (1+x)^{2 c1} * analytic + (1+x)^{2 c2} * analytic   (log if c1 = c2).
Hence for q -> oo  [derived: Darboux / transfer theorem]
   q^t C_{i,q} = sum_m s_{t,m} (2 sigma + t - m)_q / q!  +  (-1)^q sum_m [ ... q^{-2 c_j - 1 + t - m} (log q) ... ]
(the first family = Taylor coefficients of (1-x)^{m - 2 sigma - t}; 2 sigma is a half-integer for Dphi = 1.75, so no logs at x = 1).
The coefficients s_{t,m} are fitted on a window q in [Qfit_lo, Q]; the series is summed exactly up to q = Q and the
tail q > Q of the fitted main family is summed exactly:
   T(a, x) = sum_{q > Q} (a)_q/q! x^q = (a)_{Q+1}/(Q+1)! x^{Q+1} 2F1(1, Q+1+a; Q+2; x)   (mpmath).
The alternating family only enters the fit; its tail is < |last term| and is bounded in the error estimate.

The derivative sums needed by ads4/code/e2e/spinblock.pq_jets,
   S[k, l] = sum B_{n,j} (Delta+n)^k r^n P_j^{(l)}(cosh chi) = r^{-Delta} (r d_r)^k d_eta^l [r^Delta Gt],
are obtained from  E_{al,be} = r^{-Delta} th_s^al th_b^be [r^Delta Gt]  (th = rho d/drho) with r d_r = th_s + th_b,
d_eta = (1/sinh chi)(th_b - th_s).  This is free of cancellations where Gt varies on a scale (1 - rho_b) << chi (the corner);
near chi -> 0 (u -> 1) it is not, and there the ordinary radial series (blocks.CBlock) is used.
"""
import os, sys, ctypes
import numpy as np, mpmath as mp, sympy as sp
from scipy.special import gammaln, gammasgn, betainc

HERE = os.path.dirname(os.path.abspath(__file__))
E2E = os.path.join(HERE, '..', 'e2e')
if E2E not in sys.path:
    sys.path.insert(0, E2E)
from blocks import CBlock, _lib

_lc = None
def _lclib():
    global _lc
    if _lc is None:
        fn = os.path.join(HERE, 'liblcconv.so')
        if not os.path.exists(fn):
            os.system(f"gcc -O3 -march=native -fopenmp -shared -fPIC {os.path.join(HERE, 'lcconv.c')} -o {fn}")
        _lc = ctypes.CDLL(fn)
    return _lc

def lc_coefficients(blk, I, Q):
    """C_{i,q}, 0 <= i <= I, 0 <= q <= Q, from a CBlock (double precision)."""
    L = _lib('d'); size = L.block_size(blk.J, blk.N)
    Bt = np.ctypeslib.as_array((ctypes.c_double*size).from_address(blk.ptr))
    out = np.zeros((I + 1)*(Q + 1))
    _lclib().lc_convert(Bt.ctypes.data_as(ctypes.c_void_p), blk.off.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(blk.J),
                        ctypes.c_int(blk.N), ctypes.c_int(I), ctypes.c_int(Q), out.ctypes.data_as(ctypes.c_void_p))
    return out.reshape(I + 1, Q + 1)

def poch_over_fact(a, q):
    """(a)_q / q!  for real non-integer a, array q >= 0."""
    q = np.asarray(q, dtype=float)
    return gammasgn(q + a)*gammasgn(a)*np.exp(gammaln(q + a) - gammaln(q + 1) - gammaln(a))

def tail_T(a, x, Q):
    """sum_{q > Q} (a)_q/q! x^q for an array x in (-1, 1) (a non-integer).
    x > 0, a > 0:  (1-x)^{-a} I_x(Q+1, a)   (regularised incomplete beta; Taylor coefficients of (1-x)^{-a} are the negative-binomial
                   weights, so the partial sum is (1-x)^{-a} I_{1-x}(a, Q+1));
    x > 0, a < 0:  upward recursion T_{b-1} = (1-x) T_b - (b)_Q/Q! x^{Q+1} from the smallest b = a + n > 0
                   (checked against mpmath 2F1: relative error <= 1e-9 for 1-x <= 1e-3; for 1-x >~ 3e-2 the tail is < e^{-200});
    x < 0:         (a)_{Q+1}/(Q+1)! x^{Q+1} 2F1(1, Q+1+a; Q+2; x),  2F1 by the Pfaff transformation
                   (1-x)^{-1} 2F1(1, 1-a; Q+2; x/(x-1)), |x/(x-1)| <= 1/2 (60 terms)."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    pos = x > 0
    if pos.any():
        xp = x[pos]
        with np.errstate(under='ignore', over='ignore', divide='ignore', invalid='ignore'):
            if a > 0:
                T = (1 - xp)**(-a)*betainc(Q + 1, a, xp)
            else:
                n = int(np.floor(-a)) + 1
                b = a + n
                T = (1 - xp)**(-b)*betainc(Q + 1, b, xp)
                for j in range(n):
                    bb = b - j
                    T = (1 - xp)*T - poch_over_fact(bb, Q)*np.exp((Q + 1)*np.log(xp))
            T = np.where(np.isfinite(T), T, 0.0)
        out[pos] = T
    neg = ~pos
    if neg.any():
        xn = x[neg]
        z = xn/(xn - 1)
        term = np.ones_like(z); ssum = np.ones_like(z)
        for p in range(60):
            term = term*(1 - a + p)/(Q + 2 + p)*z
            ssum += term
        with np.errstate(under='ignore'):
            out[neg] = poch_over_fact(a, Q + 1)*np.sign(xn)**(Q + 1)*np.exp((Q + 1)*np.log(np.abs(xn)))/(1 - xn)*ssum
    return out

def _deta_coeffs(lmax=4):
    """d_eta^l = sum_{m=1}^{l} f_{l,m}(chi) d_chi^m  with d_eta = (1/sinh chi) d_chi;  returns numpy lambdas."""
    ch = sp.Symbol('chi'); f = sp.Function('f')
    expr = f(ch); out = {0: {0: lambda c: np.ones_like(c)}}
    for l in range(1, lmax + 1):
        expr = sp.expand(sp.diff(expr, ch)/sp.sinh(ch))
        out[l] = {}
        for m in range(1, l + 1):
            cm = sp.simplify(expr.coeff(sp.Derivative(f(ch), (ch, m))))
            out[l][m] = sp.lambdify(ch, cm, 'numpy')
    return out
_FLM = None

class LCBlock:
    """Scalar block g^{c1,c2}_{Delta,J} (c1 + c2 > 0) evaluated through the lightcone reorganisation + asymptotic tail.
    I: number of lightcone levels (power of rho_s);  Q: exact summation in rho_b up to q = Q;  fit window [qlo, Q]."""
    def __init__(self, Delta, J, c1, c2, N, I, Q, qlo=None, M=8, Malt=3, step=3, tmax=4, coeffs=None):
        assert c1 + c2 > 0
        self.Delta, self.J, self.c1, self.c2 = float(Delta), int(J), float(c1), float(c2)
        self.sig = self.c1 + self.c2
        self.N, self.I, self.Q, self.M, self.Malt, self.tmax = N, I, Q, M, Malt, tmax
        assert Q + I <= N + J, "not enough radial levels for (I, Q)"
        self.qlo = Q//2 if qlo is None else qlo
        if coeffs is None:
            blk = CBlock(Delta, J, c1, c2, N)
            self.C = lc_coefficients(blk, I, Q)
            del blk
        else:
            self.C = coeffs[:I + 1, :Q + 1]
        self._fit(step)

    # ------------------------------------------------------------------ asymptotic fit of the q-tails
    def _basis(self, q, t):
        """main family: (2 sig + t - m)_q/q!  (Taylor coefficients of (1-x)^{m-t-2sig});
        alternating family: (-1)^q (a)_q/q!, a = -2c + t - m  (Taylor coefficients of (1+x)^{2c-t+m}) for c1 != c2;
        for c1 = c2 the x = -1 singularity is logarithmic: (-1)^q q^e {1, log q}, e = -2c - 1 + t - m (tail not summed, bounded)."""
        cols, kinds = [], []
        for m in range(self.M):
            cols.append(poch_over_fact(2*self.sig + t - m, q)); kinds.append(('main', 2*self.sig + t - m))
        alt = (-1.0)**q
        cs = sorted(set([self.c1, self.c2]))
        for c in cs:
            for m in range(self.Malt):
                if len(cs) == 2:
                    a = -2*c + t - m
                    cols.append(alt*poch_over_fact(a, q)); kinds.append(('alt', a))
                else:
                    e = -2*c - 1 + t - m
                    cols.append(alt*q**e); kinds.append(('altpow', e))
                    cols.append(alt*q**e*np.log(q)); kinds.append(('altlog', e))
        return np.stack(cols, 1), kinds

    def _fit(self, step):
        q = np.arange(self.qlo, self.Q + 1, step, dtype=float)
        qi = q.astype(int)
        A0, kinds0 = self._basis(q, 0)
        nalt = A0.shape[1] - self.M
        self.coef = np.zeros((self.I + 1, self.tmax + 1, self.M))
        self.altcoef = np.zeros((self.I + 1, self.tmax + 1, nalt))
        self.altkinds = {}
        self.fitres = np.zeros((self.I + 1, self.tmax + 1))
        self.alt_last = np.zeros((self.I + 1, self.tmax + 1))     # untreated (log-case) alternating tail bound, relative
        for t in range(self.tmax + 1):
            A, kinds = self._basis(q, t)
            self.altkinds[t] = kinds[self.M:]
            sc = np.max(np.abs(A), 0); An = A/sc
            Y = self.C[:, qi]*q[None, :]**t               # (I+1, nq)
            w = 1.0/np.maximum(np.abs(Y), 1e-300)
            Aq, _ = self._basis(np.array([float(self.Q)]), t)
            for i in range(self.I + 1):
                if not np.any(Y[i]):
                    continue
                co, *_ = np.linalg.lstsq(An*w[i][:, None], Y[i]*w[i], rcond=None)
                co = co/sc
                self.coef[i, t] = co[:self.M]
                self.altcoef[i, t] = co[self.M:]
                self.fitres[i, t] = np.max(np.abs(A @ co - Y[i])*w[i])
                if kinds[self.M][0] != 'alt':
                    self.alt_last[i, t] = abs(Aq[0, self.M:] @ co[self.M:])/max(abs(Aq[0, :self.M] @ co[:self.M]), 1e-300)

    # ------------------------------------------------------------------ evaluation
    def W(self, xb):
        """W[i, t, p] = sum_q q^t C_{i,q} xb_p^q  (all q, tail included)."""
        Q = self.Q
        q = np.arange(Q + 1, dtype=float)
        npts = len(xb)
        out = np.zeros((self.I + 1, self.tmax + 1, npts))
        with np.errstate(under='ignore'):
            X = np.exp(np.outer(q, np.log(xb)))                       # (Q+1, npts)
        for t in range(self.tmax + 1):
            out[:, t, :] = (self.C*q[None, :]**t) @ X
            for m in range(self.M):
                Tm = tail_T(2*self.sig + t - m, xb, Q)
                out[:, t, :] += np.outer(self.coef[:, t, m], Tm)
            for m, kd in enumerate(self.altkinds[t]):
                if kd[0] == 'alt':
                    Tm = tail_T(kd[1], -xb, Q)
                    out[:, t, :] += np.outer(self.altcoef[:, t, m], Tm)
        return out

    def series(self, r, chi, sig=+1, kmax=4, lmax=4, cut=None, return_E=False):
        """same output as blocks.CBlock.series (sig = +1): S[p, k, l] = r^{-Delta} (r d_r)^k d_eta^l [r^Delta Gt]."""
        global _FLM
        assert sig == +1
        if _FLM is None:
            _FLM = _deta_coeffs(4)
        r = np.atleast_1d(np.asarray(r, float)); chi = np.atleast_1d(np.asarray(chi, float))
        xs, xb = r*np.exp(-chi), r*np.exp(chi)
        a = (self.Delta - self.J)/2
        Wt = self.W(xb)                                                # (I+1, 5, npts)
        nmax = min(kmax + lmax, self.tmax)        # pq_jets needs k + l <= 4
        # V[i, be, p] = sum_q (a+q)^be C_{i,q} xb^q
        from math import comb
        V = np.zeros((self.I + 1, nmax + 1, len(r)))
        for be in range(nmax + 1):
            for t in range(be + 1):
                V[:, be, :] += comb(be, t)*a**(be - t)*Wt[:, t, :]
        ii = np.arange(self.I + 1, dtype=float)
        with np.errstate(under='ignore'):
            XS = np.exp(np.outer(ii, np.log(xs)))                      # (I+1, npts)
        pref = r**(-self.J)
        E = np.zeros((nmax + 1, nmax + 1, len(r)))
        for al in range(nmax + 1):
            wi = ((a + ii)**al)[:, None]*XS
            for be in range(nmax + 1 - al):
                E[al, be] = pref*np.sum(wi*V[:, be, :], axis=0)
        self.last_iconv = np.max(np.abs(XS[-1][None, :]*V[-1, 0, :]))  # diagnostic
        self.i_tail = np.abs(XS[-1]*V[-1, 0, :])/np.maximum(np.abs(np.sum(XS*V[:, 0, :], 0)), 1e-300)
        if return_E:
            return E
        # (th_s + th_b)^k (th_b - th_s)^m  ->  sum over binomials of E_{al,be}
        from math import comb
        S = np.zeros((len(r), kmax + 1, lmax + 1))
        for k in range(kmax + 1):
            for l in range(lmax + 1 - k if lmax + 1 - k > 0 else 0):
                if k + l > nmax: continue
                tot = np.zeros(len(r))
                ms = [0] if l == 0 else range(1, l + 1)
                for m in ms:
                    acc = np.zeros(len(r))
                    for p1 in range(k + 1):              # (th_s+th_b)^k = sum C(k,p1) th_s^p1 th_b^(k-p1)
                        for p2 in range(m + 1):          # (th_b-th_s)^m = sum C(m,p2) th_b^p2 (-th_s)^(m-p2)
                            acc += comb(k, p1)*comb(m, p2)*(-1)**(m - p2)*E[p1 + m - p2, k - p1 + p2]
                    tot += _FLM[l][m](chi)*acc
                S[:, k, l] = tot
        return S
