"""
h13_laplace: skelexact3.Laplace3 (Arb Laplace/Barnes exact scalar-skeleton moments) with one change for large J (H13 scans):
the Euler--Maclaurin tail of the n-sum uses Gauss--Legendre nodes and weights computed in Arb (arb.legendre_p_root) instead of
numpy doubles.  At J >~ 40 the n-sum cancels strongly (summands ~ n^J W_n), and 1e-16 errors in the weights produced O(1e-2) errors
at large nu (not visible in the Arb radius, which treats the weights as exact).  Everything else is a verbatim copy of skelexact3.Laplace3.
Defaults raised: N0 = 400, nq = 240 (convergence checked in t_h13_nsum.py).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from fractions import Fraction as Fr
import flint
from flint import arb, ctx, arb_mat
import mack


EM_ORDER = int(os.environ.get('H13_EM', '1'))


def _arb_frac(q):
    return arb(q.numerator)/arb(q.denominator)


class Laplace3H:
    def __init__(self, Delta, J, Df, ks=(2, 3, 4, 5, 6, 7), prec=400, N0=400, nq=240, L=None, M=None, numax=0.0, imax=None):
        ctx.prec = prec
        self.prec = prec
        self.Delta, self.J, self.Df = float(Delta), int(J), float(Df)
        D = arb(Fr(self.Delta).numerator)/arb(Fr(self.Delta).denominator)
        Dfa = arb(Fr(self.Df).numerator)/arb(Fr(self.Df).denominator)
        tau = D - J
        self.tau = tau
        self.ks = tuple(ks)
        # n-nodes (as exactmom.BlockMellin), nodes/weights converted exactly from doubles
        if L is None:
            L = max(4.0*N0, 0.5*self.Delta**2/max(2*self.Df - 1.5, 0.5))
        # tail nodes/weights in Arb (H13 change: double-precision Gauss weights lose accuracy at large J, where the n-sum cancels)
        La = arb(Fr(float(L)).numerator)/arb(Fr(float(L)).denominator)
        x0 = arb(N0) - arb(0.5)
        ns = [arb(i) for i in range(N0)]
        wn = [arb(1)]*N0
        for i in range(nq):
            xg, wg = arb.legendre_p_root(nq, i, weight=True)
            uu = (xg + 1)/2; wu = wg/2
            ns.append(x0 + La*uu/(1 - uu)); wn.append(wu*La/(1 - uu)**2)
        if EM_ORDER == 1:          # -f'(x0)/24, 4-point stencil (skelexact3)
            st = ((-2, Fr(1, 288)), (-1, Fr(-8, 288)), (1, Fr(8, 288)), (2, Fr(-1, 288)))
        else:                      # -f'/24 + 7 f'''/5760 - 31 f^(5)/967680, 7-point central stencils (h = 1)
            d1 = {-3: Fr(-1, 60), -2: Fr(9, 60), -1: Fr(-45, 60), 1: Fr(45, 60), 2: Fr(-9, 60), 3: Fr(1, 60)}
            d3 = {-3: Fr(1, 8), -2: Fr(-1), -1: Fr(13, 8), 1: Fr(-13, 8), 2: Fr(1), 3: Fr(-1, 8)}
            d5 = {-3: Fr(-1, 2), -2: Fr(2), -1: Fr(-5, 2), 1: Fr(5, 2), 2: Fr(-2), 3: Fr(1, 2)}
            st = tuple((d, -d1[d]/24 + 7*d3[d]/5760 - 31*d5[d]/967680) for d in (-3, -2, -1, 1, 2, 3))
        for d, c in st:
            ns.append(x0 + d); wn.append(_arb_frac(Fr(c)))
        self.nn = len(ns)
        # Mack polynomial Q_n(T) = (-1)^J P(-T, tau + 2n) in powers of T
        P = mack.mack_P_fast(self.Delta, self.J)
        Pa = {k: _arb_frac(v) for k, v in P.items()}
        if M is None:
            M = int(3*numax**2/(2*self.Delta) + 6*np.sqrt(numax + 1) + 60)
        self.M = M
        Mx = M + J + 2
        self.S = {}
        lg = lambda v: v.lgamma()
        from math import factorial
        qs = []
        for a in range(self.nn):
            t = tau + 2*ns[a]
            q = [arb(0)]*(J + 1)
            tp = [arb(1)]
            for j in range(1, J + 1):
                tp.append(tp[-1]*t)
            for (i, j), pij in Pa.items():
                q[i] += pij*((-1)**i)*tp[j]
            if J % 2:
                q = [-v for v in q]
            qs.append(q)
        Ws = {}
        base = (lg(D - 1 + J) - lg(D - 1) + lg(D + J) - (J - 1)*arb(2).log() - 2*arb.pi().log() - 4*lg((D + J)/2) + arb(0.5).log())
        for k in self.ks:
            Ws[k] = []
            for a in range(self.nn):
                n = ns[a]
                logw = (base - lg(n + 1) - (lg(D - arb(0.5) + n) - lg(D - arb(0.5)))
                        + lg(tau/2 + n - Dfa + 1) + lg(tau/2 + n - Dfa + 1 - arb(k)/2))
                Ws[k].append(-logw.exp()*wn[a])
        for k in self.ks:
            assert k % 2 == 0, "even k only (finite partial fractions); odd k by interpolation in k"
            kap = k//2
            # 1/(z)_kap = sum_i al_i/(z+i);  1/(z (z)_kap) = a2/z^2 + a1/z + sum_{i>=1} b_i/(z+i)
            al = [Fr((-1)**i, factorial(i)*factorial(kap - 1 - i)) for i in range(kap)]
            a2 = Fr(1, factorial(kap - 1))
            a1 = -a2*sum(Fr(1, i) for i in range(1, kap))
            bb = {}
            for i in range(1, kap):
                d = Fr(i*i)
                for l in range(1, kap):
                    if l != i:
                        d *= (l - i)
                bb[i] = 1/d
            Vm = arb_mat(self.nn, Mx)
            Qm = arb_mat(self.nn, J + 1)
            for a in range(self.nn):
                n = ns[a]
                for j in range(J + 1):
                    Qm[a, j] = qs[a][j]
                W = Ws[k][a]
                c2 = 1/(2*Dfa - tau - 2*n)
                coef = {}
                for i in range(kap):
                    coef[i] = coef.get(i, arb(0)) + c2*_arb_frac(al[i])
                coef[0] = coef[0] - _arb_frac(a1)/2
                for i, bv in bb.items():
                    coef[i] = coef[i] - _arb_frac(bv)/2
                c_dbl = -_arb_frac(a2)/2
                row = [arb(0)]*Mx
                for i, cv in coef.items():
                    ah = tau/2 + n + i
                    bm = 1/ah                                     # beta(1)
                    H = 1/ah
                    Wc = W*cv; Wd = W*c_dbl
                    for m in range(Mx):
                        row[m] += (Wc + Wd*H)*bm if i == 0 else Wc*bm
                        bm = bm*(m + 1)/(ah + m + 1)
                        H = H + 1/(ah + m + 1)
                for m in range(Mx):
                    Vm[a, m] = row[m]
            Sm = Qm.transpose()*Vm
            S = [[Sm[j, m] for m in range(Mx)] for j in range(J + 1)]
            self.S[k] = S
        # prefactor pieces
        self.D = D; self.Dfa = Dfa
        al, alb = D - J - 2, D + J - 1
        self.m2 = al*alb

    def value(self, k, nu):
        """m^{2k} C_{k,nu}[Delta,J] (arb)."""
        ctx.prec = self.prec
        J = self.J
        A = 2*self.Dfa + k - arb(1.5)
        nu = arb(float(nu)) if not isinstance(nu, arb) else nu
        Mx = self.M + J + 2
        g = [arb(0)]*Mx
        gA = A.gamma()
        gm = 1/gA
        for m in range(Mx):
            g[m] = gm
            gm = gm*((A/2 + m)**2 + nu*nu/4)/((A + m)*(m + 1))
        S = self.S[k]
        I = arb(0)
        h = g
        for j in range(J + 1):
            if j > 0:
                h = [2*(m + 1)*h[m + 1] - 2*m*h[m] if m + 1 < len(h) else -2*m*h[m] for m in range(len(h) - 1)]
            I += sum((h[m]*S[j][m] for m in range(min(len(h), self.M))), arb(0))
        from exactmom import log_q
        lq = arb(float(log_q(self.Delta, J, self.Df)))   # log q in double (normalisation only, not cancelling)
        aDf = 4*self.Dfa.gamma()**2*(self.Dfa - arb(0.5)).gamma()**2
        pre = arb(2)**(1 - 2*k)*arb.pi()**arb(1.5)*aDf*self.m2**k*(-lq).exp()
        return ((-1)**J)*pre*I

    # ---- H13 addition: transposed Horner form of value() (identical sum, reordered) ----
    #   I(nu) = sum_j sum_{m<M} (D^j g)_m S_j[m] = sum_m g_m(nu) R_m,   R = S_0' + D^T(S_1' + D^T(S_2' + ...)),
    #   (D h)_m = 2(m+1) h_{m+1} - 2m h_m  =>  (D^T v)_m = 2m v_{m-1} - 2m v_m;  S_j' = S_j[:M] zero-padded to length Mx = M + J + 2.
    #   Per nu this costs O(Mx) instead of O(J Mx).  Checked against value() in t_h13_fast.py.
    def prepare(self, k):
        ctx.prec = self.prec
        if not hasattr(self, 'R'):
            self.R = {}
        J, M = self.J, self.M
        Mx = M + J + 2
        S = self.S[k]
        R = [arb(0)]*Mx
        for j in range(J, -1, -1):
            if j < J:
                R = [(2*m*R[m - 1] if m > 0 else arb(0)) - 2*m*R[m] for m in range(Mx)]
            for m in range(M):
                R[m] = R[m] + S[j][m]
        self.R[k] = R

    def value_fast(self, k, nu):
        ctx.prec = self.prec
        if not hasattr(self, 'R') or k not in self.R:
            self.prepare(k)
        J = self.J
        A = 2*self.Dfa + k - arb(1.5)
        nu = arb(float(nu)) if not isinstance(nu, arb) else nu
        R = self.R[k]
        Mx = len(R)
        I = arb(0)
        gm = 1/A.gamma()
        for m in range(Mx):
            I += gm*R[m]
            gm = gm*((A/2 + m)**2 + nu*nu/4)/((A + m)*(m + 1))
        if not hasattr(self, '_pre'):
            from exactmom import log_q
            self._lq = arb(float(log_q(self.Delta, J, self.Df)))
            self._aDf = 4*self.Dfa.gamma()**2*(self.Dfa - arb(0.5)).gamma()**2
            self._pre = True
        pre = arb(2)**(1 - 2*k)*arb.pi()**arb(1.5)*self._aDf*self.m2**k*(-self._lq).exp()
        return ((-1)**J)*pre*I
