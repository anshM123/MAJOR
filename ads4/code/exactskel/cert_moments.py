"""
cert_moments: rigorous (Arb ball) enclosures of the exact scalar-skeleton heavy moments m^{2k} C_{k,nu}[Delta, J], even k, d = 3,
identical external scalars Df, for nu in [0, Lam] -- the Laplace/Barnes form of skelexact3.py with every truncation bounded:

  I(nu) = sum_j sum_m h^(j)_m(nu) S_{j,m},  S_{j,m} = sum_{n>=0} q_{n,j} V_{n,m}   (notation of skelexact3.py; h^(j) = D^j g),
  * direct n < N (integers, no Euler-Maclaurin), all quantities in Arb (exact Mack polynomials, Gamma functions);
  * n-tail (n >= N), m < M:  |S^tail_{j,m}| <= U_{j,m}  from the Laplace kernels without partial fractions (see values in __init__),
      (exponent s = P - d_j + m + kappa), using |W_n| <= |W_N| (x_N/x_n)^P for x_n = tau/2 + n >= x_N  [certified: the polynomial
      (x+b1)(x+b2)x^P - (x+a1)(x+a2)(x+1)^P has non-negative Taylor coefficients at x = x_N, exact rationals],
      |q_{n,j}| <= Qbar_j(x) <= Qbar_j(x_N)(x/x_N)^{d_j}, beta(x+i, m+1) <= beta(x_N, m+1)(x_N/x)^{m+1}, H <= (m+1)/x;
  * m-tail (m >= M), all n:  |S_{j,m}| <= K_j beta(x_0, m+1)(1 + (m+1)/x_0)  (x_0 = tau/2), summed explicitly against |h^(j)_m(nu)| for
      M <= m < M2, and beyond M2 by  |h^(j)_m| <= (4(m+J))^J g_m(Lam)  (g_m(nu) <= g_m(Lam) for nu <= Lam; g_m(Lam) decreasing for m >= M2)
      and the certified ratio bound t_{m+1}/t_m <= (m/(m+1))^2 (polynomial check), so  sum_{m>=M2} t_m <= t_{M2}(1 + M2).
The returned value is an arb ball containing the exact moment (given the Mellin representation of CHMRS, eqs. CfromBMellin,
Bkt action, Barnes' first lemma and the Laplace representations, all [D]).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial, comb
import numpy as np
from flint import arb, ctx, arb_mat
import mack


def af(q):
    q = Fr(q)
    return arb(q.numerator)/arb(q.denominator)


def poly_nonneg_shift(coeffs_fr, x0):
    """True if the polynomial sum c_i x^i (exact Fractions) has all Taylor coefficients >= 0 at x = x0 (=> p(x) >= 0 for x >= x0)."""
    n = len(coeffs_fr)
    for r in range(n):
        s = sum(coeffs_fr[i]*comb(i, r)*x0**(i - r) for i in range(r, n))
        if s < 0:
            return False
    return True


def pmul(a, b):
    out = [Fr(0)]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x*y
    return out


def ppow(a, e):
    out = [Fr(1)]
    for _ in range(e):
        out = pmul(out, a)
    return out


def psub(a, b):
    n = max(len(a), len(b)); a = a + [Fr(0)]*(n - len(a)); b = b + [Fr(0)]*(n - len(b))
    return [x - y for x, y in zip(a, b)]


def log_q_arb(D, J, Df):
    """log q_{Delta,J} (CHMRS (2.43)), d = 3 -- skel.log_q in Arb."""
    lg = lambda v: v.lgamma()
    d = 3
    lp = (arb(2).log() + lg(D - 1) + lg(J + arb(d)/2) - 2*lg(Df) - 2*lg(Df - arb(d - 2)/2) - lg(D - arb(d)/2) - lg(arb(J) + 1)
          + 2*lg((J + D)/2) + 2*lg((D - J - d + 2)/2) - lg(J + D - 1) - lg(D - J - d + 1)
          + lg((D - J - 2*d + 2*Df + 2)/2) + lg((D + J - d + 2*Df)/2) - lg((D - J - 2*Df + 2)/2) - lg((D + J + d - 2*Df)/2))
    m2 = (D - J - d + 1)*(D + J - 1)
    nJ = 16*arb.pi()*(2*J + 1)
    return (-arb.pi().log() - lp + (2*nJ).log() - (d - 2)*m2.log() + lg(D - 1) + (2*D - d).log() - lg(D - d + 2))


class CertMoments:
    def __init__(self, Delta, J, Df, Lam, ks=(2, 4, 6, 8, 10), N=None, M=None, prec=192, verbose=False):
        ctx.prec = prec
        self.J = J; self.ks = tuple(ks)
        Dq = Fr(float(Delta)); Dfq = Fr(Df).limit_denominator(8)
        self.Dq, self.Dfq = Dq, Dfq
        D = af(Dq); Dfa = af(Dfq)
        tauq = Dq - J
        tau = D - J
        self.D, self.Dfa, self.tau = D, Dfa, tau
        self.Lam = Lam
        if N is None:
            N = int(max(3000, float(os.environ.get('CERT_NFAC', '1'))*0.6*float(Delta)**2/max(float(Df) - 1.0, 0.5)))
        if M is None:
            M = int(3*float(Lam)**2/(2*float(Delta)) + 6*np.sqrt(float(Lam) + 1) + 80)
        self.N, self.M = N, M
        P = mack.mack_P_fast(float(Delta), J)
        # q_{n,j} = (-1)^J (-1)^j sum_r P_{j r} t^r, t = tau + 2n
        Pq = {}
        for (i, r), c in P.items():
            Pq.setdefault(i, {})[r] = c
        self.dj = {j: (max(Pq[j]) if j in Pq else 0) for j in range(J + 1)}
        lg = lambda v: v.lgamma()
        base = (lg(D - 1 + J) - lg(D - 1) + lg(D + J) - (J - 1)*arb(2).log() - 2*arb.pi().log() - 4*lg((D + J)/2) + arb(0.5).log())
        x0 = tau/2
        self.x0 = x0
        Mx = M
        # q rows (arb) for n < N, and Qbar_j(x_N)
        tq = [tau + 2*n for n in range(N)]
        qrows = []
        for n in range(N):
            t = tq[n]
            tp = [arb(1)]
            for _ in range(J):
                tp.append(tp[-1]*t)
            row = []
            for j in range(J + 1):
                v = arb(0)
                for r, c in Pq.get(j, {}).items():
                    v += af(c)*tp[r]
                row.append(v*((-1)**(J + j)))
            qrows.append(row)
        xN = x0 + N
        self.xN = xN
        # |q_{n,j}| <= Qbar_j (x_n/x_N)^{d_j} for x_n >= x_N, with Qbar_j = |q_j(x_N)| prod_r (x_N + |r|)/(x_N - |r|) over the roots r of
        # q_j(x) (rigorous enclosures from flint), valid when all |r| < x_N; fallback: sum of |coefficients| (2 x_N)^r.
        from flint import fmpq_poly, fmpq
        Qbar = []
        xNq = Fr(tauq, 2) + N if False else tauq/2 + N
        for j in range(J + 1):
            cj = Pq.get(j, {})
            if not cj:
                Qbar.append(arb(0)); continue
            dmax = max(cj)
            coeffs = [Fr(0)]*(dmax + 1)
            for r, c in cj.items():
                coeffs[r] += c*Fr(2)**r*((-1)**(J + j))
            fb = sum((abs(af(c))*(2*xN)**r for r, c in cj.items()), arb(0))
            if dmax == 0:
                Qbar.append(abs(af(coeffs[0]))); continue
            poly = fmpq_poly([fmpq(c.numerator, c.denominator) for c in coeffs])
            roots = poly.complex_roots()
            qN = sum((af(c)*xN**i for i, c in enumerate(coeffs)), arb(0))
            prod = abs(qN)
            ok = True
            for rt, mult in roots:
                ra = abs(rt).upper() if hasattr(abs(rt), 'upper') else abs(rt)
                ra = arb(ra)
                if not (ra < xN):
                    ok = False; break
                prod = prod*((xN + ra)/(xN - ra))**mult
            Qbar.append(prod if ok else fb)
        self.S = {}; self.U = {}; self.K = {}; self.info = {}
        for k in self.ks:
            assert k % 2 == 0
            kap = k//2
            al = [Fr((-1)**i, factorial(i)*factorial(kap - 1 - i)) for i in range(kap)]
            a2 = Fr(1, factorial(kap - 1)); a1 = -a2*sum(Fr(1, i) for i in range(1, kap))
            bb = {}
            for i in range(1, kap):
                d = Fr(i*i)
                for l in range(1, kap):
                    if l != i:
                        d *= (l - i)
                bb[i] = 1/d
            # weights W_n (sign included), n < N
            Wl = []
            for n in range(N):
                lw = (base - lg(arb(n) + 1) - (lg(D - arb(0.5) + n) - lg(D - arb(0.5)))
                      + lg(tau/2 + n - Dfa + 1) + lg(tau/2 + n - Dfa + 1 - kap))
                Wl.append(-lw.exp())
            WN = (base - lg(arb(N) + 1) - (lg(D - arb(0.5) + N) - lg(D - arb(0.5))) + lg(tau/2 + N - Dfa + 1)
                  + lg(tau/2 + N - Dfa + 1 - kap)).exp()
            # certified power P for |W_n| <= |W_N| (x_N/x_n)^P
            xNq = tauq/2 + N
            a1q, a2q = 1 - Dfq, 1 - Dfq - kap
            b1q, b2q = 1 - tauq/2, Dq - Fr(1, 2) - tauq/2
            Pexp = int(np.floor(float(b1q + b2q - a1q - a2q) - 0.25))
            while Pexp > 1:
                f = psub(pmul(pmul([b1q, Fr(1)], [b2q, Fr(1)]), [Fr(0)]*Pexp + [Fr(1)]),
                         pmul(pmul([a1q, Fr(1)], [a2q, Fr(1)]), ppow([Fr(1), Fr(1)], Pexp)))
                if poly_nonneg_shift(f, xNq):
                    break
                Pexp -= 1
            # V rows and S (direct)
            Vm = arb_mat(N, Mx); Qm = arb_mat(N, J + 1)
            Ksum = [arb(0)]*(J + 1)
            for n in range(N):
                for j in range(J + 1):
                    Qm[n, j] = qrows[n][j]
                W = Wl[n]
                c2 = 1/(2*Dfa - tau - 2*n)
                coef = {}
                for i in range(kap):
                    coef[i] = coef.get(i, arb(0)) + c2*af(al[i])
                coef[0] = coef[0] - af(a1)/2
                for i, bv in bb.items():
                    coef[i] = coef[i] - af(bv)/2
                c_dbl = -af(a2)/2
                Ct = max(abs(c2), arb(1)/(2*kap))       # direct (non partial-fraction) Laplace bound, see header
                for j in range(J + 1):
                    Ksum[j] += abs(qrows[n][j])*abs(W)*Ct
                row = [arb(0)]*Mx
                for i, cv in coef.items():
                    ah = tau/2 + n + i
                    bm = 1/ah; H = 1/ah
                    Wc = W*cv; Wd = W*c_dbl
                    for m in range(Mx):
                        row[m] += (Wc + Wd*H)*bm if i == 0 else Wc*bm
                        bm = bm*(m + 1)/(ah + m + 1)
                        H = H + 1/(ah + m + 1)
                for m in range(Mx):
                    Vm[n, m] = row[m]
            Sm = Qm.transpose()*Vm
            self.S[k] = [[Sm[j, m] for m in range(Mx)] for j in range(J + 1)]
            # tail constants for n >= N (Laplace kernels of Gamma(z)/Gamma(z+kap) and Gamma(z)/(z Gamma(z+kap)), no partial fractions:
            #   |V_{n,m}| <= |W_n|/Gamma(kap) [ |c2_n| beta(x_n, m+kap) + beta(x_n, m+kap+1)(x_n+m+kap)/(2 kap (x_n - 1)) ]
            gk = arb(kap).gamma()
            U = []
            for j in range(J + 1):
                row = []
                b1 = (lg(xN) + lg(arb(kap)) - lg(xN + kap)).exp()        # beta(x_N, kap)
                for m in range(Mx):
                    b2 = b1*(m + kap)/(xN + m + kap)                         # beta(x_N, m+kap+1)
                    s1 = Pexp - self.dj[j] + m + kap
                    assert s1 > 1
                    row.append(WN*Qbar[j]/gk*(b1/(2*xN - 2*Dfa)*(1 + xN/(s1 - 1))
                                              + b2*(xN + m + kap)/(2*kap*(xN - 1))*(1 + xN/s1)))
                    b1 = b2
                U.append(row)
            self.U[k] = U
            # K_j: |S_{j,m}| <= K_j bx_m,  bx_m = beta(x0, m+kap)(1 + (x0+m+kap)/(x0-1))/Gamma(kap)
            Ktail = []
            for j in range(J + 1):
                s = Pexp - self.dj[j]
                assert s > 1
                Ktail.append(WN*Qbar[j]*max(1/(2*xN - 2*Dfa), arb(1)/(2*kap))*(1 + xN/(s - 1)))
            self.K[k] = [Ksum[j] + Ktail[j] for j in range(J + 1)]
            self.info[k] = dict(P=Pexp, N=N, M=Mx)
            if verbose:
                print(f"  k={k}: P={Pexp} N={N} M={Mx}", flush=True)
        self.m2 = (D - J - 2)*(D + J - 1)
        self.logq = log_q_arb(D, J, Dfa)
        self.aDf = 4*Dfa.gamma()**2*(Dfa - arb(0.5)).gamma()**2

    def _M2(self, k):
        A = 2*self.Dfa + k - arb(1.5)
        Lam = arb(self.Lam)
        M2 = int(float((A*A/4 + Lam*Lam/4 - A).mid())) + 10
        M2 = max(M2, self.M + 10)
        J = self.J
        x0q = self.Dq/2 - Fr(J, 2) if False else (self.Dq - J)/2
        kap = k//2
        # certify t_{m+1}/t_m <= (m/(m+1))^2 for m >= M2:
        #   m^2 (m+J)^J (x0+m+kap)(2x0+m+kap-1) - (m+1)^2 (m+1+J)^J (m+kap)(2x0+m+kap) >= 0
        while True:
            lhs = pmul(pmul(pmul([Fr(0), Fr(0), Fr(1)], ppow([Fr(J), Fr(1)], J)), [x0q + kap, Fr(1)]), [2*x0q + kap - 1, Fr(1)])
            rhs = pmul(pmul(pmul(ppow([Fr(1), Fr(1)], 2), ppow([Fr(1 + J), Fr(1)], J)), [Fr(kap), Fr(1)]), [2*x0q + kap, Fr(1)])
            if poly_nonneg_shift(psub(lhs, rhs), Fr(M2)):
                return M2
            M2 = int(M2*1.5) + 10

    def values(self, k, nus):
        """list of arb balls containing m^{2k} C_{k,nu} for nu in nus (floats, 0 <= nu <= Lam)."""
        ctx.prec = ctx.prec
        J = self.J
        A = 2*self.Dfa + k - arb(1.5)
        gA = A.gamma()
        M = self.M
        x0 = self.x0
        kap = k//2
        gk = arb(kap).gamma()
        S, U, K = self.S[k], self.U[k], self.K[k]
        Lam = arb(self.Lam)
        if not hasattr(self, '_M2c'):
            self._M2c = {}
        if k not in self._M2c:
            # choose M2 (>= the certified ratio point) so that the far tail is <= 1e-9 of the absolute size of the m-sum at nu = 0
            M2 = self._M2(k)
            g0 = []; gm = 1/gA
            for m in range(M + J + 2):
                g0.append(gm); gm = gm*((A/2 + m)**2)/((A + m)*(m + 1))
            h = g0; scale = arb(0)
            for j in range(J + 1):
                if j > 0:
                    h = [2*(m + 1)*h[m + 1] - 2*m*h[m] for m in range(len(h) - 1)]
                scale += sum((abs(h[m]*S[j][m]) for m in range(M)), arb(0))
            Ksum = sum((K[j] for j in range(J + 1)), arb(0))
            for _ in range(8):
                b = (x0.lgamma() + (arb(M2) + kap).lgamma() - (x0 + M2 + kap).lgamma()).exp()
                bxM2 = b*(1 + (x0 + M2 + kap)/(x0 - 1))/gk
                gL = 1/gA
                for m in range(M2):
                    gL = gL*((A/2 + m)**2 + Lam*Lam/4)/((A + m)*(m + 1))
                far = Ksum*(4*(arb(M2) + J))**J*gL*bxM2*(1 + M2)
                if far.upper() <= 1e-9*scale.lower():
                    break
                M2 *= 2
            self._M2c[k] = (M2, far)
        M2, far = self._M2c[k]
        L = M2 + J + 2
        # bx_m = beta(x0, m+kap)(1 + (x0+m+kap)/(x0-1))/Gamma(kap) for m <= M2
        bx = []; b = (x0.lgamma() + arb(kap).lgamma() - (x0 + kap).lgamma()).exp()
        for m in range(M2 + 1):
            bx.append(b*(1 + (x0 + m + kap)/(x0 - 1))/gk); b = b*(m + kap)/(x0 + m + kap)
        logpre = (arb(2).log()*(1 - 2*k) + arb.pi().log()*arb(1.5) + self.aDf.log() + k*self.m2.log() - self.logq)
        pre = ((-1)**J)*logpre.exp()
        out = []
        for nu in nus:
            nu = arb(float(nu))
            g = [arb(0)]*L; gm = 1/gA
            for m in range(L):
                g[m] = gm
                gm = gm*((A/2 + m)**2 + nu*nu/4)/((A + m)*(m + 1))
            I = arb(0); rad = arb(0)
            h = g
            for j in range(J + 1):
                if j > 0:
                    h = [2*(m + 1)*h[m + 1] - 2*m*h[m] for m in range(len(h) - 1)]
                Sj, Uj = S[j], U[j]
                for m in range(M):
                    I += h[m]*Sj[m]
                    rad += abs(h[m])*Uj[m]
                for m in range(M, M2):
                    rad += abs(h[m])*K[j]*bx[m]
            tot = pre*(I + arb(0, (rad + far).upper()))
            out.append(tot)
        return out
