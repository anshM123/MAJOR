"""
Exact scalar-skeleton heavy moments without any T-contour: Laplace/Barnes form, in Arb (python-flint) multiprecision.

Ingredients [D]:
 * Barnes' first lemma:  int dT/(4 pi i) Gamma(T/2)^2 gamma_{a-T}(nu) e^{T u} / gamma_a(nu)^2 = 2F1(c, cbar; A; 1 - e^{-2u})/Gamma(A) =: G(u),
   c = (A + i nu)/2, A = a - 1/2  (checked numerically).  So the e^{-pi nu} of gamma_a^2 cancels analytically.
   G(u) = sum_m g_m w^m,  w = 1 - e^{-2u},  g_m = |(c)_m|^2 / ((A)_m m! Gamma(A))   (real, positive).
 * T^j e^{Tu} = d_u^j e^{Tu};  in the w basis d_u w^m = 2m w^{m-1} - 2m w^m.
 * Bhat_{k,T} = sum_n W_n Q_n(T) R_n(T),  R_n = Gamma(z)/Gamma(z+kappa) [c2_n - 1/(2z)],  z = (T_n - T)/2,  T_n = tau + 2n,  kappa = k/2
   (both parities), with the Laplace forms (Re z > 0 on the original contour)
      Gamma(z)/Gamma(z+kappa)      = (2/Gamma(kappa)) int_0^oo du e^{-(T_n - T)u} w^{kappa-1},
      Gamma(z)/(z Gamma(z+kappa))  = (2/Gamma(kappa)) int_0^oo du e^{-(T_n - T)u} sum_i w^{kappa+i}/(kappa+i).
 * int_0^oo e^{-T_n u} w^{x-1} du = (1/2) Gamma(T_n/2) Gamma(x)/Gamma(T_n/2 + x) =: beta_n(x)/2.
 => I(nu) = sum_j sum_m (D^j g)_m S_{j,m},
    S_{j,m} = sum_n W_n q_{n,j} (1/Gamma(kappa)) [ c2_n beta_n(m+kappa) - (1/2) sum_i beta_n(m+kappa+i+1)/(kappa+i) ],
    Q_n(T) = sum_j q_{n,j} T^j  (exact Mack polynomial, mack.py, evaluated in Arb from the rationals).
 m^{2k} C = (-1)^J 2^{1-2k} pi^{3/2} a_Df m^{2k} / q_{Delta,J} * I   (same normalisation as exactmom.py).
The n-sum uses direct terms n < N0 and an Euler--Maclaurin tail with Gauss nodes (as exactmom.py), all weights in Arb.
Only real arithmetic.  Precision set by `prec` (bits); Arb's error radius on the result is reported.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from fractions import Fraction as Fr
import flint
from flint import arb, ctx, arb_mat
import mack


def _arb_frac(q):
    return arb(q.numerator)/arb(q.denominator)


class Laplace3:
    def __init__(self, Delta, J, Df, ks=(2, 3, 4, 5, 6, 7), prec=400, N0=120, nq=160, L=None, M=None, numax=0.0, imax=None):
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
        x, w = np.polynomial.legendre.leggauss(nq)
        uu = 0.5*(x + 1); wu = 0.5*w
        x0 = N0 - 0.5
        nint = x0 + L*uu/(1 - uu); wint = wu*L/(1 - uu)**2
        st = x0 + np.array([-2.0, -1.0, 1.0, 2.0]); wst = np.array([1.0, -8.0, 8.0, -1.0])/12.0/24.0
        ns = [arb(i) for i in range(N0)] + [arb(float(v)) for v in nint] + [arb(float(v)) for v in st]
        wn = [arb(1)]*N0 + [arb(float(v)) for v in wint] + [arb(float(v)) for v in wst]
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
