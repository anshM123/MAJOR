"""
d = 3 scalar blocks near the crossed lightcone rho_b -> 1 (the corner u = sqrt(rhobar_u/rho_u) -> 0 of the Regge-moment
integration region, SPINNING_E2E.md sect. 3.2), by a lightcone reorganisation of the radial series plus an asymptotic tail.

Notation (as in ads4/code/e2e/blocks.py; the reflected block g^{c1,c2}, c1 + c2 > 0, is evaluated at positive radial coordinates):
   Gt(r, eta) = sum_{n,j} B_{n,j} r^n P_j(eta),    g = (4r)^Delta Gt,    rho_s = r e^{-chi}, rho_b = r e^{chi}, eta = cosh chi.
The derivative sums needed by e2e/spinblock.pq_jets are
   S[k, l] = sum_{n,j} B_{n,j} (Delta+n)^k r^n P_j^{(l)}(eta).
With P_j^{(l)} = (2l-1)!! C^{(l+1/2)}_{j-l} and C^lam_m(cos th) = sum_k a^lam_k a^lam_{m-k} e^{i(m-2k)th}, a^lam_k = (lam)_k/k! > 0:
   S[k, l] = (2l-1)!! r^{l-J} sum_{i,q >= 0} (Delta - J + l + i + q)^k C^{(l)}_{i,q} rho_s^i rho_b^q,
   C^{(l)}_{i,q} = sum_{k'} a^lam_{k'} a^lam_{|q-i|+k'} B_{i+q-J+l, |q-i|+2k'+l}        (lcconv.c; all terms >= 0, no cancellation).
At fixed i, H^{(l)}_i(x) = sum_q C^{(l)}_{i,q} x^q has on |x| <= 1 only the singularities
   x = 1  (zbar' = 1):  analytic + (1-x)^{-2 sigma - l} * analytic,   sigma = c1 + c2
            (local exponents 0, -sigma of the Casimir equation at zbar = 1; 1 - zbar = ((1-x)/(1+x))^2; each d_eta adds one power),
   x = -1 (zbar' = oo): (1+x)^{2 c1 - l} * analytic + (1+x)^{2 c2 - l} * analytic   (logarithmic if c1 = c2).
Hence for q -> oo (Darboux)
   q^t C^{(l)}_{i,q} = sum_m s_m (2 sigma + l + t - m)_q/q!  +  (-1)^q sum_{c, m} a_{c,m} (l - 2c + t - m)_q/q!     [c1 != c2]
(main family: Taylor coefficients of (1-x)^{m-t-l-2 sigma}; 2 sigma is a half-integer for Dphi = 1.75, so no logs at x = 1).
The coefficients are fitted on q in [Q/2, Q]; the series is summed exactly up to q = Q and the tails in closed form (tail_T).
"""
import os, sys, ctypes
from math import comb
import numpy as np
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

def lc_coefficients(blk, I, Q, l=0):
    """C^{(l)}_{i,q}, 0 <= i <= I, 0 <= q <= Q, from a CBlock (double precision)."""
    L = _lib('d'); size = L.block_size(blk.J, blk.N)
    Bt = np.ctypeslib.as_array((ctypes.c_double*size).from_address(blk.ptr))
    out = np.zeros((I + 1)*(Q + 1))
    _lclib().lc_convert_l(Bt.ctypes.data_as(ctypes.c_void_p), blk.off.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(blk.J),
                          ctypes.c_int(blk.N), ctypes.c_int(I), ctypes.c_int(Q), ctypes.c_int(l), out.ctypes.data_as(ctypes.c_void_p))
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

DFACT = [1, 1, 3, 15, 105]

class LCBlock:
    """Scalar block g^{c1,c2}_{Delta,J} (c1 + c2 > 0) through the lightcone reorganisation + asymptotic tail.
    I: lightcone levels (powers of rho_s);  Q: exact summation in rho_b up to q = Q;  fit window [qlo, Q];  M main + Malt alternating
    asymptotic terms;  lmax: highest eta-derivative.  coeffs: optional list of C^{(l)} arrays (at least (I+1) x (Q+1))."""
    def __init__(self, Delta, J, c1, c2, N, I, Q, qlo=None, M=8, Malt=3, step=3, lmax=4, coeffs=None):
        assert c1 + c2 > 0
        self.Delta, self.J, self.c1, self.c2 = float(Delta), int(J), float(c1), float(c2)
        self.sig = self.c1 + self.c2
        self.N, self.I, self.Q, self.M, self.Malt, self.lmax = N, I, Q, M, Malt, lmax
        assert Q + I + lmax <= N + J, "not enough radial levels for (I, Q)"
        self.qlo = Q//2 if qlo is None else qlo
        if coeffs is None:
            blk = CBlock(Delta, J, c1, c2, N)
            self.C = [lc_coefficients(blk, I, Q, l) for l in range(lmax + 1)]
            del blk
        else:
            self.C = [np.ascontiguousarray(c[:I + 1, :Q + 1]) for c in coeffs[:lmax + 1]]
        self._fit(step)

    # ------------------------------------------------------------------ asymptotic fit of the q-tails
    def _basis(self, q, l, t):
        cols, kinds = [], []
        for m in range(self.M):
            cols.append(poch_over_fact(2*self.sig + l + t - m, q)); kinds.append(('main', 2*self.sig + l + t - m))
        alt = (-1.0)**q
        cs = sorted(set([self.c1, self.c2]))
        for c in cs:
            for m in range(self.Malt):
                if len(cs) == 2:
                    a = -2*c + l + t - m
                    cols.append(alt*poch_over_fact(a, q)); kinds.append(('alt', a))
                else:
                    e = -2*c - 1 + l + t - m
                    cols.append(alt*q**e); kinds.append(('altpow', e))
                    cols.append(alt*q**e*np.log(q)); kinds.append(('altlog', e))
        return np.stack(cols, 1), kinds

    def _fit(self, step):
        q = np.arange(self.qlo, self.Q + 1, step, dtype=float)
        qi = q.astype(int)
        self.fits = {}
        self.fitres = np.zeros((self.I + 1, self.lmax + 1, self.lmax + 1))
        for l in range(self.lmax + 1):
            for t in range(self.lmax + 1 - l):
                A, kinds = self._basis(q, l, t)
                # one least-squares problem with I+1 right-hand sides; common smooth weight q^{-(2 sig + l - 1 + t)}
                w = q**(-(2*self.sig + l - 1 + t))
                Aw = A*w[:, None]
                sc = np.max(np.abs(Aw), 0); An = Aw/sc
                Y = self.C[l][:, qi]*q[None, :]**t
                co, *_ = np.linalg.lstsq(An, (Y*w[None, :]).T, rcond=None)
                co = (co/sc[:, None]).T                       # (I+1, ncols)
                self.fits[(l, t)] = (co, kinds)
                fit = co @ A.T
                scale = np.maximum(np.max(np.abs(Y), axis=1), 1e-300)
                self.fitres[:, l, t] = np.where(np.any(Y != 0, axis=1), np.max(np.abs(fit - Y), axis=1)/scale, 0.0)

    # ------------------------------------------------------------------ evaluation
    def W(self, xb, l):
        """W[i, t, p] = sum_q q^t C^{(l)}_{i,q} xb_p^q  (all q: exact up to Q, fitted asymptotic tail beyond)."""
        Q = self.Q
        q = np.arange(Q + 1, dtype=float)
        out = np.zeros((self.I + 1, self.lmax + 1 - l, len(xb)))
        with np.errstate(under='ignore'):
            X = np.exp(np.outer(q, np.log(xb)))                       # (Q+1, npts)
        for t in range(self.lmax + 1 - l):
            out[:, t, :] = (self.C[l]*q[None, :]**t) @ X
            co, kinds = self.fits[(l, t)]
            for m, kd in enumerate(kinds):
                if kd[0] == 'main':
                    out[:, t, :] += np.outer(co[:, m], tail_T(kd[1], xb, Q))
                elif kd[0] == 'alt':
                    out[:, t, :] += np.outer(co[:, m], tail_T(kd[1], -xb, Q))
                # 'altpow'/'altlog' (c1 = c2): tail not summed (alternating, < |last term|, relative <~ Q^{-2(sigma+c)} ~ 1e-8)
        return out

    def series(self, r, chi, sig=+1, kmax=4, lmax=4, cut=None):
        """same output as blocks.CBlock.series (sig = +1): S[p, k, l] = sum B (Delta+n)^k r^n P_j^{(l)}(cosh chi), k + l <= 4."""
        assert sig == +1
        r = np.atleast_1d(np.asarray(r, float)); chi = np.atleast_1d(np.asarray(chi, float))
        xs, xb = r*np.exp(-chi), r*np.exp(chi)
        ii = np.arange(self.I + 1, dtype=float)
        with np.errstate(under='ignore'):
            XS = np.exp(np.outer(ii, np.log(xs)))                      # (I+1, npts)
        S = np.zeros((len(r), kmax + 1, lmax + 1))
        for l in range(min(lmax, self.lmax) + 1):
            Wl = self.W(xb, l)
            a = self.Delta - self.J + l + ii                           # (Delta + n) = a_i + q
            pref = DFACT[l]*r**(l - self.J)
            for k in range(min(kmax, self.lmax - l) + 1):
                V = np.zeros((self.I + 1, len(r)))
                for t in range(k + 1):
                    V += comb(k, t)*(a**(k - t))[:, None]*Wl[:, t, :]
                S[:, k, l] = pref*np.sum(XS*V, axis=0)
                if k == 0 and l == 0:
                    self.i_tail = np.abs(XS[-1]*V[-1])/np.maximum(np.abs(S[:, 0, 0]/pref), 1e-300)
        return S

    def E(self, r, chi, nmax=2):
        """E[al, be, p] = r^{-Delta} th_s^al th_b^be [r^Delta Gt] (th = rho d/drho), al + be <= nmax (for the Casimir check)."""
        r = np.atleast_1d(np.asarray(r, float)); chi = np.atleast_1d(np.asarray(chi, float))
        xs, xb = r*np.exp(-chi), r*np.exp(chi)
        a = (self.Delta - self.J)/2
        W0 = self.W(xb, 0)
        ii = np.arange(self.I + 1, dtype=float)
        XS = np.exp(np.outer(ii, np.log(xs)))
        E = np.zeros((nmax + 1, nmax + 1, len(r)))
        for be in range(nmax + 1):
            V = sum(comb(be, t)*a**(be - t)*W0[:, t, :] for t in range(be + 1))
            for al in range(nmax + 1 - be):
                E[al, be] = r**(-self.J)*np.sum(((a + ii)**al)[:, None]*XS*V, axis=0)
        return E
