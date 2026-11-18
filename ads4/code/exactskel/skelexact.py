"""
Production evaluation of exact scalar-skeleton heavy moments m^{2k} C_{k,nu}[Delta,J], k = 2..7, on many nu at once
(Mellin representation, exactmom.BlockMellin; see exactmom.py header for the method).  Lines and residue nodes are shared
between the k of one parity (A_k = 2Df + k - 3/2 lie on one lattice A_kb + 2Z per parity).

ExactState(Delta, J, Df).values(nus)  -> {k: array},  plus .diag (cancellation ratio, contour index per k, nu).
Odd k: analytic continuation in k (Gamma-ratio Pochhammers); equals the k-continued CHMRS series to 1e-10 at small nu/m (t_check2).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, gammaln
from exactmom import BlockMellin, log_q, log_aDf, LOGPI

JSET = (0, 1, 2, 3, 4, 6, 8, 11, 14, 18, 23, 29, 36, 44, 53, 63, 75, 88, 103, 120, 140, 165, 195)


class ExactState:
    def __init__(self, Delta, J, Df, h=0.1, tmax_extra=45.0, **kw):
        self.B = BlockMellin(Delta, J, Df, **kw)
        self.Delta, self.J, self.Df = float(Delta), int(J), float(Df)
        self.tau = self.Delta - J
        self.m2 = (self.Delta - J - 2)*(self.Delta + J - 1)
        self.logq = log_q(self.Delta, J, Df)
        self.laDf = log_aDf(Df)
        self.h = h
        self.tmax_extra = tmax_extra
        self._lines = {}

    def _line(self, ks, J1, tmax):
        key = (tuple(ks), J1, tmax)
        if key not in self._lines:
            Ab = 2*self.Df + ks[0] - 1.5
            c1 = Ab + 2*J1 + 1
            ts = np.arange(0.0, tmax + self.h/2, self.h)
            Lb, Bv = self.B.Bhat_multi(ks, c1 + 1j*ts)
            self._lines[key] = (ts, c1 + 1j*ts, Lb, Bv)
        return self._lines[key]

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
            j0 = int(np.floor((self.tau - 0.5 - Ab - 1)/2))
            jset = sorted(set([j for j in JSET if j <= j0] + [j0]))
            Am = 2*self.Df + kset[len(kset)//2] - 1.5
            res = np.zeros((len(kset), nus.size)); canc = np.zeros((len(kset), nus.size)); jcs = np.zeros(nus.size, int)
            Lref = None
            for i, nu in enumerate(nus):
                Ts = (Am*Am + nu*nu)/(2*Am)
                jid = (Ts - Ab - 1)/2
                J1 = max([j for j in jset if j <= max(jid, 0)] or [0])
                jcs[i] = J1
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
                    jp = jl[sh:] - sh                         # pole order j' of Gamma((A-T+i nu)/2)
                    Tk = Tj[sh:]
                    if jp.size:
                        lr = (2*loggamma(Tk/2) + loggamma(-1j*nu - jp) - gammaln(jp + 1) - g4 + (Lb2[a] - Lref[a]))
                        r = np.exp(lr)*(-2.0)*((-1.0)**jp)*Bj[a, sh:]
                        Ir = -np.sum(np.real(r)); mr = np.sum(np.abs(r))
                    else:
                        Ir = 0.0; mr = 0.0
                    res[a, i] = Il + Ir
                    canc[a, i] = (mag + mr)/max(abs(Il + Ir), 1e-300)
            for a, k in enumerate(kset):
                logpre = (1 - 2*k)*np.log(2.0) + 1.5*LOGPI + self.laDf + Lref[a] - self.logq + k*np.log(self.m2)
                out[k] = (-1.0)**self.J*np.exp(logpre)*res[a]
                self.diag[k] = (canc[a], jcs)
        return out
