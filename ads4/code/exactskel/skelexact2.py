"""
Exact scalar-skeleton moments with the T-contour moved across the Bhat poles (needed for Delta_g >= 300, where nu reaches ~ tau and the
contour left of tau leaves a cancellation ~ 10^{0.09 nu}).

  I = int_{c1} dT/(4 pi i) E Bhat  -  sum_{complex poles A+-i nu+2j < c1} Re[...]   (as skelexact.py)
                                    -  (1/2) sum_{real poles T = tau+2n+2i < c1} Res[E b_n]
The real-pole residues are computed per n by a 24-point circle integral of radius 1/2 (the next singularity of E b_n is at distance 2:
error ~ 4^{-24}).  For lines with c1 > tau the n-sum is direct up to N0 >= (c1-tau)/2 + 150, so that the Euler--Maclaurin tail never
meets the singularities of the summand (which lie at Re n <= (c1-tau)/2).
c1 per nu: the smallest lattice point A_kb + 2j + 1 (quantised j) at which the Gamma factor |E(c1+it)/gamma_a^2| stays below CANC
on the whole line (estimate of the residual cancellation; the realised one is returned in .diag).
Validated against skelexact.py (contour left of tau) at Delta_g = 100 (t_cross.py).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, gammaln
from exactmom import BlockMellin, log_q, log_aDf, LOGPI

LOGCANC = np.log(1e2)
NCIRC = 24


def _jq(j):
    """quantise a lattice index upward to a coarse geometric set (limits the number of distinct lines)."""
    if j <= 4:
        return int(j)
    q = int(np.ceil(np.log(j)/np.log(1.15)))
    return int(np.ceil(1.15**q))


def _panels(a, b, L, nq):
    """Gauss-Legendre nodes/weights on [a, b] (a < b or a > b), panels of length <= L with nq nodes each."""
    npan = max(1, int(np.ceil(abs(b - a)/L)))
    x, w = np.polynomial.legendre.leggauss(nq)
    e = np.linspace(a, b, npan + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5*(hi - lo)*x[None, :] + 0.5*(hi + lo)).ravel(), (0.5*(hi - lo)*w[None, :]).ravel()


class ExactState2:
    def __init__(self, Delta, J, Df, numax, h=0.1, tmax_extra=45.0, cmax=None, **kw):
        self.Delta, self.J, self.Df = float(Delta), int(J), float(Df)
        self.tau = self.Delta - J
        self.m2 = (self.Delta - J - 2)*(self.Delta + J - 1)
        self.h = h; self.tmax_extra = tmax_extra
        A2 = 2*Df + 0.5
        self.cmax = cmax if cmax is not None else max(self.tau, 6.0*numax + 60)
        N0 = int(max(400, (self.cmax - self.tau)/2 + 150))
        self.B = BlockMellin(Delta, J, Df, N0=N0, **kw)
        self.logq = log_q(self.Delta, J, Df)
        self.laDf = log_aDf(Df)
        self._lines = {}

    def choose(self, Ab, Am, nu, tmax):
        """lattice index J1 (c1 = Ab + 2 J1 + 1): smallest quantised J1 with max_t log|E/gamma^2| <= LOGCANC."""
        ts = np.linspace(0, tmax, 400)
        g4 = 4*np.real(loggamma((Am + 1j*nu)/2))
        j = 0
        jmax = int((self.cmax - Ab - 1)/2)
        while True:
            jq = min(_jq(j), jmax)
            c = Ab + 2*jq + 1
            if abs((c - self.tau)/2 - np.round((c - self.tau)/2)) < 0.2:      # keep away from the real Bhat poles
                jq = max(jq - 1, 0); c = Ab + 2*jq + 1
            T = c + 1j*ts
            lE = np.real(2*loggamma(T/2) + loggamma((Am - T + 1j*nu)/2) + loggamma((Am - T - 1j*nu)/2)) - g4
            if lE.max() <= LOGCANC or jq >= jmax:
                return jq
            j = jq + 1

    def _line(self, ks, J1, tmax):
        key = (tuple(ks), J1, tmax)
        if key not in self._lines:
            Ab = 2*self.Df + ks[0] - 1.5
            c1 = Ab + 2*J1 + 1
            ts = np.arange(0.0, tmax + self.h/2, self.h)
            Lb, Bv = self.B.Bhat_multi(ks, c1 + 1j*ts)
            self._lines[key] = (ts, c1 + 1j*ts, Lb, Bv)
        return self._lines[key]

    def _bn_poles(self, k, c1, nu, A, Lref):
        """(1/2) sum of residues of E(T) b_n(T) at the real poles T < c1 (all n), relative to e^{Lref}, divided by gamma_a^2
        (vectorised over all poles (n, i); circle integral of radius 1/2, NCIRC points)."""
        tau = self.tau
        nmax = int(np.floor((c1 - tau)/2))
        if nmax < 0:
            return 0.0
        B = self.B
        ns = np.arange(nmax + 1)
        if k % 2 == 0:
            nn, ii = np.meshgrid(ns, np.arange(k//2), indexing='ij')
        else:
            nn, ii = np.meshgrid(ns, np.arange(nmax + 1), indexing='ij')
        nn = nn.ravel(); ii = ii.ravel()
        Tp = tau + 2*nn + 2*ii
        sel = Tp < c1
        nn, ii, Tp = nn[sel], ii[sel], Tp[sel]
        if nn.size == 0:
            return 0.0
        lw = B.logw(k)[nn]
        th = 2*np.pi*(np.arange(NCIRC) + 0.5)/NCIRC
        circ = 0.5*np.exp(1j*th)
        T = Tp[:, None] + circ[None, :]
        yy = T - tau/2 - nn[:, None]
        if B.cheb:
            z = yy/B.R[nn][:, None]
            b1 = np.zeros(z.shape, dtype=complex); b2 = np.zeros(z.shape, dtype=complex)
            for q in range(B.J, 0, -1):
                b1, b2 = B.C[nn, q][:, None] + 2*z*b1 - b2, b1
            Qv = B.C[nn, 0][:, None] + z*b1 - b2
        else:
            Qv = np.zeros(T.shape, dtype=complex) + B.G[nn, -1][:, None]
            for q in range(B.J - 1, -1, -1):
                Qv = Qv*yy + B.G[nn, q][:, None]
        z = (tau - T)/2 + nn[:, None]
        if k % 2 == 0:
            poch = np.ones(T.shape, dtype=complex)
            for q in range(k//2):
                poch = poch*(z + q)
            inv = 1.0/poch
        else:
            inv = np.exp(loggamma(z) - loggamma(z + k/2))
        br = B.c2[nn][:, None] + 1.0/(T - tau - 2*nn[:, None])
        g4 = 4*np.real(loggamma((A + 1j*nu)/2))
        lE = 2*loggamma(T/2) + loggamma((A - T + 1j*nu)/2) + loggamma((A - T - 1j*nu)/2) - g4
        f = np.exp(lE + (lw - Lref)[:, None])*(-1.0)*Qv*inv*br
        return float(np.real(np.sum(f*circ[None, :]))/NCIRC/2.0)

    def values_box(self, nus, ks=(2, 3, 4, 5, 6, 7), delta=1.0, L=1.0, nq=8, ctarget=None):
        """Box contour: up the line Re T = c1 for |Im T| < nu - delta, horizontally at Im T = +-(nu - delta) to Re T = c0 = A_kb - 1,
        then up Re T = c0 to infinity.  It never crosses the rows of Gamma poles A +- i nu + 2j (no complex residues: their sum
        cancels by ~10^{0.06 nu}), only the real Bhat poles T = tau + 2n + 2i < c1 (circle residues, tiny).  c1 = the midpoint
        tau + 2n + 1 (n >= 0) or tau - 1 nearest to ctarget(nu) (default: smallest c with max|E/gamma^2| <= e^LOGCANC on the line)."""
        nus = np.atleast_1d(np.asarray(nus, float))
        out = {}; self.diag = {}
        for par in (0, 1):
            kset = [k for k in ks if k % 2 == par]
            if not kset:
                continue
            kb = kset[0]
            Ab = 2*self.Df + kb - 1.5
            Am = 2*self.Df + kset[len(kset)//2] - 1.5
            c0 = Ab - 1.0
            res = np.zeros((len(kset), nus.size)); canc = np.zeros((len(kset), nus.size)); c1s = np.zeros(nus.size)
            Lref = None
            for i, nu in enumerate(nus):
                if ctarget is None:
                    J1 = self.choose(Ab, Am, nu, float(nu + self.tmax_extra))
                    ct = Ab + 2*J1 + 1
                else:
                    ct = ctarget(nu)
                if ct < self.tau:
                    c1 = min(ct, self.tau - 1.0)
                else:
                    c1 = self.tau + 2*np.round((ct - self.tau - 1)/2) + 1
                c1 = max(c1, c0 + 2.0)
                c1s[i] = c1
                tb = nu - delta
                pts, wts = [], []
                if tb > 0.5:
                    t1, w1 = _panels(0.0, tb, L, nq)                       # up Re T = c1
                    pts.append(c1 + 1j*t1); wts.append(1j*w1)
                    x2, w2 = _panels(c1, c0, L, nq)                        # left along Im T = tb
                    pts.append(x2 + 1j*tb); wts.append(w2*np.sign(c0 - c1)*np.sign(c0 - c1))
                    t3, w3 = _panels(tb, tb + self.tmax_extra + 20, L, nq)  # up Re T = c0
                    pts.append(c0 + 1j*t3); wts.append(1j*w3)
                else:
                    t1, w1 = _panels(0.0, nu + self.tmax_extra + 20, L, nq)
                    pts.append(c1 + 1j*t1); wts.append(1j*w1)
                    c1s[i] = c1
                Tn = np.concatenate(pts); Wn = np.concatenate(wts)
                Lb, Bv = self.B.Bhat_multi(kset, Tn)
                if Lref is None:
                    Lref = Lb.copy()
                for a, k in enumerate(kset):
                    A = 2*self.Df + k - 1.5
                    g4 = 4*np.real(loggamma((A + 1j*nu)/2))
                    lE = 2*loggamma(Tn/2) + loggamma((A - Tn + 1j*nu)/2) + loggamma((A - Tn - 1j*nu)/2) - g4
                    f = np.exp(lE + (Lb[a] - Lref[a]))*Bv[a]
                    Ic = 2*np.real(np.sum(Wn*f))/(4*np.pi*1j) if False else 2*np.real(np.sum(Wn*f)/(4j*np.pi))
                    mag = 2*np.sum(np.abs(Wn*f))/(4*np.pi)
                    if tb > 0.5:
                        Ir, mr = 0.0, 0.0
                    else:                      # small nu: plain line at c1 with the few complex residues A + 2j < c1
                        jl = np.arange(max(0, int(np.floor((c1 - A)/2)) + 1)) if c1 > A else np.arange(0)
                        if jl.size:
                            Tj = A + 1j*nu + 2*jl
                            Lb2, Bj = self.B.Bhat_multi([k], Tj)
                            lr = (2*loggamma(Tj/2) + loggamma(-1j*nu - jl) - gammaln(jl + 1) - g4 + (Lb2[0] - Lref[a]))
                            r = np.exp(lr)*(-2.0)*((-1.0)**jl)*Bj[0]
                            Ir = -np.sum(np.real(r)); mr = np.sum(np.abs(r))
                        else:
                            Ir, mr = 0.0, 0.0
                    Ip = self._bn_poles(k, c1, nu, A, Lref[a])
                    tot = Ic + Ir - Ip
                    res[a, i] = tot
                    canc[a, i] = (mag + mr + abs(Ip))/max(abs(tot), 1e-300)
            for a, k in enumerate(kset):
                logpre = (1 - 2*k)*np.log(2.0) + 1.5*LOGPI + self.laDf + Lref[a] - self.logq + k*np.log(self.m2)
                out[k] = (-1.0)**self.J*np.exp(logpre)*res[a]
                self.diag[k] = (canc[a], c1s)
        return out

    def values(self, nus, ks=(2, 3, 4, 5, 6, 7)):
        nus = np.atleast_1d(np.asarray(nus, float))
        out = {}; self.diag = {}
        tmax = float(5*np.ceil((nus.max() + self.tmax_extra)/5))
        for par in (0, 1):
            kset = [k for k in ks if k % 2 == par]
            if not kset:
                continue
            kb = kset[0]
            Ab = 2*self.Df + kb - 1.5
            Am = 2*self.Df + kset[len(kset)//2] - 1.5
            res = np.zeros((len(kset), nus.size)); canc = np.zeros((len(kset), nus.size)); c1s = np.zeros(nus.size)
            Lref = None
            for i, nu in enumerate(nus):
                tmi = float(5*np.ceil((nu + self.tmax_extra)/5))
                J1 = self.choose(Ab, Am, nu, tmi)
                c1 = Ab + 2*J1 + 1
                c1s[i] = c1
                ts, T, Lb, Bv = self._line(kset, J1, tmax)
                if Lref is None:
                    Lref = Lb.copy()
                wtr = np.full(ts.size, self.h); wtr[0] = self.h/2
                jl = np.arange(J1 + 1)
                Tj = Ab + 1j*nu + 2*jl
                Lb2, Bj = self.B.Bhat_multi(kset, Tj)
                for a, k in enumerate(kset):
                    A = 2*self.Df + k - 1.5
                    g4 = 4*np.real(loggamma((A + 1j*nu)/2))
                    lE = 2*loggamma(T/2) + loggamma((A - T + 1j*nu)/2) + loggamma((A - T - 1j*nu)/2) - g4
                    f = np.exp(lE + (Lb[a] - Lref[a]))*Bv[a]
                    Il = 2*np.real(np.sum(wtr*f))/(4*np.pi)
                    mag = 2*np.sum(wtr*np.abs(f))/(4*np.pi)
                    sh = (k - kb)//2
                    jp = jl[sh:] - sh
                    Tk = Tj[sh:]
                    if jp.size:
                        lr = (2*loggamma(Tk/2) + loggamma(-1j*nu - jp) - gammaln(jp + 1) - g4 + (Lb2[a] - Lref[a]))
                        r = np.exp(lr)*(-2.0)*((-1.0)**jp)*Bj[a, sh:]
                        Ir = -np.sum(np.real(r)); mr = np.sum(np.abs(r))
                    else:
                        Ir = 0.0; mr = 0.0
                    Ip = self._bn_poles(k, c1, nu, A, Lref[a])
                    tot = Il + Ir - Ip
                    res[a, i] = tot
                    canc[a, i] = (mag + mr + abs(Ip))/max(abs(tot), 1e-300)
            for a, k in enumerate(kset):
                logpre = (1 - 2*k)*np.log(2.0) + 1.5*LOGPI + self.laDf + Lref[a] - self.logq + k*np.log(self.m2)
                out[k] = (-1.0)**self.J*np.exp(logpre)*res[a]
                self.diag[k] = (canc[a], c1s)
        return out
