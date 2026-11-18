"""
Exact heavy-state action of the CHMRS functional C_{k,nu} on a d = 3 scalar block (identical external scalars Df), from the Mellin
representation of CHMRS (arXiv:2106.10274) eqs. (CfromBMellin), (Bkt action):

  C_{k,nu}[G] = 2^{1-2k} pi^{3/2} a_Df / gamma_a(nu)^2  int_{c0} dT/(4 pi i)  Gamma(T/2)^2 gamma_{a-T}(nu)  Bhat_{k,T}[G],  a = 2Df + k - 1,
  Bhat_{k,T}[G_{Delta,J}] = sum_{n>=0} (-1)^{k/2} Q^n(T') / [((2Df-tau)/2 - n)_{k/2} ((tau-T)/2 + n)_{k/2}] [1/(2Df-tau-2n) + 1/(T-tau-2n)],
  heavy density  C_{k,nu}[Delta,J] = (-1)^J C_{k,nu}[G] / (q_{Delta,J} 2 sin^2(pi(tau-2Df)/2))         (CHMRS (2.31)-(2.34)),
with Q^n from mack.py (the sin^2 cancels analytically).  0 < c0 < A := a - 1/2.

Numerics (double precision; no 1/m expansion anywhere):
 * gamma_{a-T}(nu) = Gamma((A-T+i nu)/2) Gamma((A-T-i nu)/2) makes the integrand ~ e^{-pi nu/2} on Re T = c0 while the result is
   ~ gamma_a^2 ~ e^{-pi nu}: a cancellation of up to 10^57 at nu = 100.  The contour is therefore moved to Re T = c1 (just left of
   the first Bhat pole T = tau), picking up the residues at T_j = A +- i nu + 2j (A + 2j < c1):
      I = int_{c1} dT/(4 pi i) E Bhat  -  sum_j Re[ Gamma(T_j/2)^2 (-2)(-1)^j/j! Gamma(-i nu - j) Bhat(T_j) ].
   The remaining cancellation on Re T = c1 is <= 10^9 for nu <= 98 at Delta_g = 100 (printed per state by `cancel`).
 * n-sum: direct for n < N0, the tail by the Euler--Maclaurin midpoint form  int_{N0-1/2}^oo f dn + f'(N0-1/2)/24  (all singularities of
   the summand in n lie at Re n < 0 for the T used here), the integral by Gauss--Legendre on a mapped variable.
 * line integral: trapezoid in Im T (step h), conjugate symmetry.
Odd k: Pochhammers continued as  (-1)^{k/2}/((2Df-tau)/2-n)_{k/2} -> Gamma(tau/2+n-Df+1-k/2)/Gamma(tau/2+n-Df+1),
       1/((tau-T)/2+n)_{k/2} -> Gamma((tau-T)/2+n)/Gamma((tau-T)/2+n+k/2)  (analytic in k; tested against the k-continued CHMRS series).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import mpmath as mp
from scipy.special import loggamma, gammaln
import mack

LOGPI = np.log(np.pi)


def log_q(Delta, J, Df):
    import skel
    return float(skel.log_q(Delta, J, Df))


def log_aDf(Df):
    return float(mp.log(4*mp.gamma(Df)**2*mp.gamma(Df - mp.mpf(1)/2)**2))


class BlockMellin:
    """Bhat_{k,T} of one block (Delta, J) for identical scalars Df: n-nodes, Mack coefficients, log-weights."""

    def __init__(self, Delta, J, Df, N0=400, nq=240, L=None, cheb=None):
        self.Delta, self.J, self.Df = float(Delta), int(J), float(Df)
        tau = self.Delta - J
        self.tau = tau
        # n-nodes: direct 0..N0-1 ; EM integral on [N0-1/2, oo) ; derivative stencil at N0-1/2
        nd = np.arange(N0, dtype=float)
        if L is None:
            L = max(4.0*N0, 0.5*self.Delta**2/max(2*Df - 1.5, 0.5))
        x, w = np.polynomial.legendre.leggauss(nq)
        uu = 0.5*(x + 1); wu = 0.5*w
        x0 = N0 - 0.5
        nint = x0 + L*uu/(1 - uu)
        wint = wu*L/(1 - uu)**2
        st = x0 + np.array([-2.0, -1.0, 1.0, 2.0])
        wst = np.array([1.0, -8.0, 8.0, -1.0])/12.0/24.0      # + f'(N0-1/2)/24 (midpoint Euler-Maclaurin)
        self.ns = np.concatenate([nd, nint, st])
        self.wn = np.concatenate([np.ones(N0), wint, wst])
        self.cheb = (J > 10) if cheb is None else cheb
        if self.cheb:
            self.R, self.C = mack.Qpoly_cheb(self.Delta, J, self.ns)   # stable Chebyshev form (J > 10)
        else:
            self.G = mack.Qpoly_ycoeffs(self.Delta, J, self.ns)          # [nn, J+1]
        D = self.Delta
        n = self.ns
        # log K_n without the Gamma(tau/2+n-Df+1)^2 (combined with the k-dependent ratio below)
        self.logK0 = (gammaln(D - 1 + J) - gammaln(D - 1) + gammaln(D + J) - (J - 1)*np.log(2.0) - 2*LOGPI
                      - gammaln(n + 1) - (gammaln(D - 0.5 + n) - gammaln(D - 0.5)) - 4*gammaln((D + J)/2))
        self.g1 = gammaln(tau/2 + n - Df + 1)
        self.c2 = 1.0/(2*Df - tau - 2*n)

    def Qeval(self, y):
        """Qpoly^n at y (array [nT, nn]): Horner (monomial, J <= 10) or Clenshaw (Chebyshev, scale R_n)."""
        if not self.cheb:
            Q = np.zeros(y.shape, dtype=complex) + self.G[None, :, -1]
            for i in range(self.J - 1, -1, -1):
                Q = Q*y + self.G[None, :, i]
            return Q
        z = y/self.R[None, :]
        b1 = np.zeros(y.shape, dtype=complex); b2 = np.zeros(y.shape, dtype=complex)
        for i in range(self.J, 0, -1):
            b1, b2 = self.C[None, :, i] + 2*z*b1 - b2, b1
        return self.C[None, :, 0] + z*b1 - b2

    def logw(self, k):
        """log of  K_n * Gamma(tau/2+n-Df+1-k/2)/Gamma(tau/2+n-Df+1)  (the (-1)^{k/2}/(..)_{k/2} factor), times 1/2 (heavy density)."""
        return self.logK0 + self.g1 + gammaln(self.tau/2 + self.ns - self.Df + 1 - k/2) + np.log(0.5)

    def Bhat(self, k, T, chunk=400):
        """returns (log scale Lb, array of Bhat(T)/e^Lb), T complex array.  Includes the factor -1/2 of the density (-K_n/2 Q):
        Bhat_density(T) = sum_n (-1/2) K_n R_k(n) Qpoly^n(T) / ((tau-T)/2+n)_{k/2} [c2_n + 1/(T - tau - 2n)]."""
        T = np.atleast_1d(np.asarray(T, dtype=complex))
        lw = self.logw(k)
        Lb = lw.max()
        wt = -np.exp(lw - Lb)*self.wn          # (-1/2 included in logw via log 0.5, sign here)
        out = np.zeros(T.shape, dtype=complex)
        tau = self.tau
        for i0 in range(0, T.size, chunk):
            Tc = T.flat[i0:i0 + chunk][:, None]
            y = Tc - tau/2 - self.ns[None, :]
            # Horner for Qpoly in y
            Q = self.Qeval(y)
            z = (tau - Tc)/2 + self.ns[None, :]
            if k % 2 == 0:
                poch = np.ones(z.shape, dtype=complex)
                for q in range(k//2):
                    poch = poch*(z + q)
                inv = 1.0/poch
            else:
                inv = np.exp(loggamma(z) - loggamma(z + k/2))
            br = self.c2[None, :] + 1.0/(Tc - tau - 2*self.ns[None, :])
            out.flat[i0:i0 + chunk] = np.sum(wt[None, :]*Q*inv*br, axis=1)
        return Lb, out

    def Bhat_multi(self, ks, T, chunk=300):
        """Bhat for several k of one parity at once: returns (array Lb[k], array [len(ks), len(T)])."""
        ks = list(ks)
        T = np.atleast_1d(np.asarray(T, dtype=complex)).ravel()
        lws = [self.logw(k) for k in ks]
        Lbs = np.array([lw.max() for lw in lws])
        wts = np.array([-np.exp(lw - Lb)*self.wn for lw, Lb in zip(lws, Lbs)])      # [nk, nn]
        out = np.zeros((len(ks), T.size), dtype=complex)
        tau = self.tau
        par = ks[0] % 2
        kmax = max(ks)
        for i0 in range(0, T.size, chunk):
            Tc = T[i0:i0 + chunk][:, None]
            y = Tc - tau/2 - self.ns[None, :]
            Q = self.Qeval(y)
            z = (tau - Tc)/2 + self.ns[None, :]
            br = self.c2[None, :] + 1.0/(Tc - tau - 2*self.ns[None, :])
            base = Q*br
            # 1/((z)_{k/2}) built incrementally
            if par == 0:
                inv = 1.0/z; kk = 2
            else:
                inv = np.exp(loggamma(z) - loggamma(z + 1.5)); kk = 3
            while kk <= kmax:
                if kk in ks:
                    out[ks.index(kk), i0:i0 + chunk] = (base*inv) @ wts[ks.index(kk)]
                inv = inv/(z + kk/2)
                kk += 2
        return Lbs, out


def logE(T, A, nu):
    """log of Gamma(T/2)^2 Gamma((A-T+i nu)/2) Gamma((A-T-i nu)/2) / |Gamma((A+i nu)/2)|^4."""
    return (2*loggamma(T/2) + loggamma((A - T + 1j*nu)/2) + loggamma((A - T - 1j*nu)/2)
            - 4*np.real(loggamma((A + 1j*nu)/2)))


class HeavyC:
    """m^{2k} C_{k,nu}[Delta,J] (exact, Mellin) for one block and a set of nu, k."""

    def __init__(self, Delta, J, Df, h=0.1, tmax_extra=45.0, jshift=0, **kw):
        self.B = BlockMellin(Delta, J, Df, **kw)
        self.Delta, self.J, self.Df = float(Delta), int(J), float(Df)
        self.tau = self.Delta - J
        al, alb = self.Delta - J - 2, self.Delta + J - 1
        self.m2 = al*alb
        self.logq = log_q(self.Delta, J, Df)
        self.laDf = log_aDf(Df)
        self.h = h
        self.tmax_extra = tmax_extra
        self.jshift = jshift
        self._line = {}

    def contour(self, k):
        A = 2*self.Df + k - 1.5
        # c1 midway between consecutive A + 2j, at distance >= 0.5 left of tau
        j0 = int(np.floor((self.tau - 0.5 - A - 1)/2))
        c1 = A + 2*j0 + 1
        return A, c1, j0

    def line(self, k, tmax, jc):
        key = (k, tmax, jc)
        if key not in self._line:
            A = 2*self.Df + k - 1.5
            c1 = A + 2*jc + 1
            ts = np.arange(0.0, tmax + self.h/2, self.h)
            T = c1 + 1j*ts
            Lb, Bv = self.B.Bhat(k, T)
            self._line[key] = (ts, T, Lb, Bv)
        return self._line[key]

    def jchoice(self, k, nu):
        """contour index: c1 = A + 2 jc + 1 near the saddle T* = (A^2 + nu^2)/(2A) of the Gamma factors, capped left of tau."""
        A, c1max, j0 = self.contour(k)
        Ts = (A*A + nu*nu)/(2*A)
        return int(np.clip(np.round((Ts - A - 1)/2 + self.jshift), 0, j0))

    def value(self, k, nus, tmax=None, diag=False, jc=None):
        """m^{2k} C_{k,nu}[Delta,J] for nu in nus (array).  diag: also return cancellation ratios and contour index."""
        nus = np.atleast_1d(np.asarray(nus, float))
        A = 2*self.Df + k - 1.5
        res = np.zeros(nus.size); canc = np.zeros(nus.size); jcs = np.zeros(nus.size, int)
        Lb0 = None
        for i, nu in enumerate(nus):
            jci = self.jchoice(k, nu) if jc is None else jc
            jcs[i] = jci
            tm = float(np.ceil(nu + self.tmax_extra)) if tmax is None else tmax
            tm = float(5*np.ceil(tm/5))
            ts, T, Lb, Bv = self.line(k, tm, jci)
            if Lb0 is None:
                Lb0 = Lb
            Bv = Bv*np.exp(Lb - Lb0); Lb = Lb0
            wtr = np.full(ts.size, self.h); wtr[0] = self.h/2
            js = np.arange(jci + 1)
            lE = logE(T, A, nu)
            f = np.exp(lE)*Bv
            Iline = 2*np.real(np.sum(wtr*f))/(4*np.pi)
            mag = 2*np.sum(wtr*np.abs(f))/(4*np.pi)
            # residues at T_j = A + i nu + 2j, j = 0..j0  (A + 2 j0 < c1)
            Tj = A + 1j*nu + 2*js
            Lb2, Bj = self.B.Bhat(k, Tj)
            lr = (2*loggamma(Tj/2) + loggamma(-1j*nu - js) - gammaln(js + 1) - 4*np.real(loggamma((A + 1j*nu)/2)) + (Lb2 - Lb))
            r = np.exp(lr)*(-2.0)*((-1.0)**js)*Bj
            Ires = -np.sum(np.real(r))
            tot = Iline + Ires
            res[i] = tot
            canc[i] = (mag + np.sum(np.abs(r)))/max(abs(tot), 1e-300)
        # density: (-1)^J/q * 2^{1-2k} pi^{3/2} a_Df * (I/gamma^2) * e^{Lb} ; times m^{2k}
        logpre = (1 - 2*k)*np.log(2.0) + 1.5*LOGPI + self.laDf + Lb0 - self.logq + k*np.log(self.m2)
        out = (-1.0)**self.J*np.exp(logpre)*res
        if diag:
            return out, canc, jcs
        return out


if __name__ == '__main__':
    import skel
    # quick test against the CHMRS order-8 series (v2 of uniform: Delta = 30.3, J = 0, Df = 3)
    for (Delta, J, Df) in ((30.3, 0, 3.0), (30.3, 2, 3.0)):
        hc = HeavyC(Delta, J, Df)
        for k in (2, 4):
            nus = np.array([0.05, 2.0, 5.0])
            v, c, jj = hc.value(k, nus, diag=True)
            for nu, vv, cc, j in zip(nus, v, c, jj):
                s = skel.series_value('seriesC', k, Df, nu, Delta, J)
                print(f"Delta={Delta} J={J} k={k} nu={nu}: exact {vv:.12f} series {s:.12f} rel {vv/s-1:+.2e} cancel {cc:.1e} jc {j}")
