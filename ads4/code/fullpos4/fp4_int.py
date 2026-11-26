"""
fullpos4/fp4_int.py (FULLPOS4.md §2) [DV]: the frequency-integrated outgoing amplitude for the target weight psi2, numpy.

  h(nu) dnu = 2 psi2(p) dp, nu = Lam p  (so G^psi(z) = 2 int psi2 F_{3,Lam p}(z) dp = (1-z)^{-3} [P+_h + P-_h]).
  P+_nu(beta) = -(i/pi) coth(pi nu) Q_{-1/2-i nu}(cosh beta)                         (checked: fp4 §2)
  => P+_h(beta) = -(i/pi) int_0^oo W_h(beta+u) (4 sinh(beta+u/2) sinh(u/2))^{-1/2} du,   (= author §58 Euler form)
     W_h(t) = 2 int_0^1 psi2(p) coth(pi Lam p) e^{i Lam p t} dp.
  Split for |t| >= T0 (exact, psi polynomial, (1-p)^5 | psi):  W_h = W_S + W_E + W_C,
     W_S(t) = -2 sum_j (-1)^j psi^(j)(0)/(i Lam t)^{j+1}           (p = 0 endpoint; non-oscillatory)
     W_E(t) =  2 e^{i Lam t} sum_{j>=5} (-1)^j psi^(j)(1)/(i Lam t)^{j+1}   (p = 1 endpoint; carrier e^{i Lam t})
     W_C(t) =  2 int psi(p) 2/(e^{2 pi Lam p} - 1) e^{i Lam p t} dp  (coth correction; |Im t| < 2 pi)
  u-integral: [0, u0] direct (W_h by Gauss-Legendre in p, |t| small), [u0, oo) split: S + C on the real axis,
  E on the ray u = u0 + r e^{i pi/4} (analytic continuation; e^{i Lam t} decays there).
"""
import sys, os, math
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fullpos3'))
import fp3_common as F3

LAM = math.sqrt(F3.LAM2)
_p2, _p3 = F3.nf02()


def derivs(c, x0):
    D = len(c) - 1; out = []
    for j in range(D + 1):
        s = 0
        for n in range(j, D + 1):
            ff = 1
            for i in range(j):
                ff *= (n - i)
            s += c[n]*ff*x0**(n - j)
        out.append(float(s))
    return out


class Weight:
    def __init__(self, k=2, lam=None):
        self.lam = LAM if lam is None else lam
        LAMl = self.lam
        self.c = _p2 if k == 2 else _p3
        self.cf = np.array([float(x) for x in self.c])
        self.d0 = derivs(self.c, 0); self.d1 = derivs(self.c, 1)
        assert all(abs(self.d1[j]) == 0 for j in range(5))
        xg, wg = np.polynomial.legendre.leggauss(200 if LAMl > 100 else 48)
        self.pg = (xg + 1)/2; self.wg = wg/2
        self.psig = np.polyval(self.cf[::-1], self.pg)
        self.cothg = 1/np.tanh(math.pi*LAMl*self.pg)
        # coth correction nodes on [0, 0.03]
        xc, wc = np.polynomial.legendre.leggauss(80)
        self.pc = 0.015*(xc + 1); self.wc = 0.015*wc
        self.psic = np.polyval(self.cf[::-1], self.pc)
        cw = min(0.015, 6.0/LAMl)
        self.pc = cw*(xc + 1); self.wc = cw*wc
        self.psic = np.polyval(self.cf[::-1], self.pc)
        self.corc = 2/np.expm1(2*math.pi*LAMl*self.pc)

    def psi(self, p):
        return np.polyval(self.cf[::-1], p)

    def W_direct(self, t):
        """W_h(t) by 200-node Gauss-Legendre in p (accurate for Lam |t| <~ 60)."""
        t = np.asarray(t, complex)
        return 2*(np.exp(1j*self.lam*np.multiply.outer(t, self.pg))*(self.psig*self.cothg*self.wg)).sum(-1)

    def W_S(self, t):
        t = np.asarray(t, complex); x = 1j*self.lam*t
        s = np.zeros_like(t)
        for j in range(len(self.d0) - 1, -1, -1):          # Horner in 1/x
            s = s/x + (-1)**j*self.d0[j]
        return -2*s/x

    def W_Eamp(self, t):
        """W_E(t) e^{-i Lam t}."""
        t = np.asarray(t, complex); x = 1j*self.lam*t
        s = np.zeros_like(t)
        for j in range(len(self.d1) - 1, -1, -1):
            s = s/x + (-1)**j*self.d1[j]
        return 2*s/x

    def W_C(self, t):
        t = np.asarray(t, complex)
        return 2*(np.exp(1j*self.lam*np.multiply.outer(t, self.pc))*(self.psic*self.corc*self.wc)).sum(-1)


def _gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return a + (b - a)*(x + 1)/2, (b - a)/2*w


def _graded(a, b, h0, r, n):
    e = [a]; h = h0
    while e[-1] + h < b:
        e.append(e[-1] + h); h *= r
    e.append(b)
    X, Wt = [], []
    for lo, hi in zip(e[:-1], e[1:]):
        x, w = _gl(lo, hi, n); X.append(x); Wt.append(w)
    return np.concatenate(X), np.concatenate(Wt)


def kern(beta, u):
    return 1/np.sqrt(4*np.sinh(beta + u/2)*np.sinh(u/2))


class Amp:
    """P+_h(beta) for arrays of complex beta (Re beta > 0, |Im beta| < pi/2 or so).
    |beta| <  B0: u in [0, u0] with W_h direct (Gauss-Legendre in p), u >= u0 split (S, C real axis; E on a pi/4 ray from u0);
    |beta| >= B0: split on the whole u-range (S, C real axis from 0; E on the ray u = r e^{i pi/4}); u = s^2 near 0."""
    B0 = 0.04

    def __init__(self, k=2, u0=0.06, lam=None, nq=16):
        self.w = Weight(k, lam); self.u0 = u0; self.lam = self.w.lam
        if lam is not None:
            self.B0 = 40.0/self.lam; u0 = 60.0/self.lam; self.u0 = u0
        s, ws = _graded(0.0, math.sqrt(u0), 1e-4, 1.6, nq)
        self.un = s*s; self.wun = 2*s*ws
        self.uf, self.wuf = _graded(u0, 80.0, 0.02, 1.25, nq)
        r, wr = _graded(0.0, 80.0/self.lam, 1e-3, 1.4, nq)
        self.ue = u0 + r*np.exp(1j*math.pi/4); self.wue = wr*np.exp(1j*math.pi/4)
        # whole-range nodes (|beta| >= B0)
        s2, ws2 = _graded(0.0, math.sqrt(80.0), 1e-3, 1.25, nq)
        self.ua = s2*s2; self.wua = 2*s2*ws2
        r2, wr2 = _graded(0.0, math.sqrt(100.0/self.lam), 1e-3, 1.3, nq)
        rot = np.exp(1j*math.pi/4)
        self.uea = r2*r2*rot; self.wuea = 2*r2*wr2*rot

    def _small(self, b):
        near = (self.w.W_direct(b + self.un[None, :])*kern(b, self.un[None, :])*self.wun).sum(1)
        tf = b + self.uf[None, :]
        sc = ((self.w.W_S(tf) + self.w.W_C(tf))*kern(b, self.uf[None, :])*self.wuf).sum(1)
        te = b + self.ue[None, :]
        e = (np.exp(1j*self.lam*self.ue[None, :])*self.w.W_Eamp(te)*kern(b, self.ue[None, :])*self.wue).sum(1)
        return near + sc, sc, e

    def _large(self, b):
        tf = b + self.ua[None, :]
        sc = ((self.w.W_S(tf) + self.w.W_C(tf))*kern(b, self.ua[None, :])*self.wua).sum(1)
        te = b + self.uea[None, :]
        e = (np.exp(1j*self.lam*self.uea[None, :])*self.w.W_Eamp(te)*kern(b, self.uea[None, :])*self.wuea).sum(1)
        return sc, sc, e


    def majorant(self, beta, factored=False):
        """Lemma-A+-type majorant (FULLPOS4 §3): (1/pi) int_0^oo |K| (|W_S| + |W_C| + e^{-Lam Im t}|W_E e^{-i Lam t}|) du on the real
        u-axis for |beta| >= B0; for |beta| < B0 the part u < u0 uses |W_h| (direct).  Not certified (quadrature), DV only."""
        beta = np.atleast_1d(np.asarray(beta, complex)); out = np.zeros(len(beta)); outE = np.zeros(len(beta))
        sm = np.abs(beta) < self.B0
        for i0 in range(0, len(beta), 64):
            b = beta[i0:i0 + 64][:, None]; m = sm[i0:i0 + 64]
            tf = b + self.ua[None, :]
            kk = np.abs(kern(b, self.ua[None, :]))*self.wua
            big = ((np.abs(self.w.W_S(tf)) + np.abs(self.w.W_C(tf)))*kk).sum(1)
            bigE = (np.exp(-self.lam*self.ua[None, :].imag)*np.abs(self.w.W_Eamp(tf))*kk).sum(1)
            if m.any():
                bs = b[m]
                near = (np.abs(self.w.W_direct(bs + self.un[None, :]))*np.abs(kern(bs, self.un[None, :]))*self.wun).sum(1)
                tf2 = bs + self.uf[None, :]
                k2 = np.abs(kern(bs, self.uf[None, :]))*self.wuf
                far = ((np.abs(self.w.W_S(tf2)) + np.abs(self.w.W_C(tf2)))*k2).sum(1)
                farE = (np.abs(self.w.W_Eamp(tf2))*k2).sum(1)
                big[m] = near + far; bigE[m] = farE
            out[i0:i0 + 64] = big/math.pi; outE[i0:i0 + 64] = bigE/math.pi
        if factored:
            return out, outE                  # majorant = out + e^{-Lam Im beta} outE
        return out + np.exp(-self.lam*beta.imag)*outE

    def __call__(self, beta, split=False, factored=False):
        """P+_h(beta).  factored=True returns (P_rest, Etil) with P+_h = P_rest + e^{i Lam beta} Etil (overflow-safe)."""
        beta = np.atleast_1d(np.asarray(beta, complex))
        rest = np.zeros_like(beta); oS = np.zeros_like(beta); oE = np.zeros_like(beta)
        sm = np.abs(beta) < self.B0
        for mask, fn in ((sm, self._small), (~sm, self._large)):
            idx = np.nonzero(mask)[0]
            for i0 in range(0, len(idx), 64):
                ii = idx[i0:i0 + 64]
                t, sc, e = fn(beta[ii][:, None])
                rest[ii] = -1j/math.pi*t; oS[ii] = -1j/math.pi*sc; oE[ii] = -1j/math.pi*e
        if factored:
            return rest, oE
        tot = rest + np.exp(1j*self.lam*beta)*oE
        return (tot, oS, np.exp(1j*self.lam*beta)*oE) if split else tot


if __name__ == '__main__':
    # self-test: (i) single-frequency Q-form vs mpmath P+; (ii) integrated amplitude 2 Re P+_h(beta) vs Mehler OmLR (real beta)
    import mpmath as mp
    mp.mp.dps = 30
    for nu, b in [(3.3, 0.4 + 0.2j), (40.0, 1.1 + 0.3j), (7.0, 0.05 + 0.01j)]:
        f = lambda u: mp.exp(1j*nu*(b + u))/mp.sqrt(4*mp.sinh(b + u/2)*mp.sinh(u/2))
        q = mp.quad(f, [0, 0.01, 0.1, 1, 5, 20, 60])
        v1 = -1j/mp.pi*mp.coth(mp.pi*nu)*q
        nn = mp.mpf(nu); bb = mp.mpc(b)
        v2 = mp.gamma(1j*nn)/(mp.sqrt(mp.pi)*mp.gamma(0.5 + 1j*nn))*mp.exp((-0.5 + 1j*nn)*bb)*mp.hyp2f1(0.5, 0.5 - 1j*nn, 1 - 1j*nn, mp.exp(-2*bb))
        print('single nu', nu, b, complex(v1), complex(v2), float(abs(v1 - v2)/abs(v2)))
    sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
    import ls_core as LSC
    A = Amp()
    bs = np.array([0.01, 0.03, 0.039, 0.041, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0, 8.0])
    v = A(bs + 0j)
    ref = 2*LSC.OmLR(2, bs)
    for b_, x, r in zip(bs, v, ref):
        print(f'beta={b_:5.2f}  2Re P+_h={2*x.real:+.10e}  Mehler={r:+.10e}  |P+_h|={abs(x):.3e}  Im={x.imag:+.3e}')
