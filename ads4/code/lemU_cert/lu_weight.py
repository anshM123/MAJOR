"""
lemU_cert/lu_weight.py [Arb]: the realised (EPSW-repaired) weight of G_A, read-only from EPSW_CERT's data and code.

  f(mu) = chi(mu - i) R_-(mu - i)  (epswcert/ew_common.py, exact JSON doubles, repaired x_re[0]),
  theta(nu) = (1/2)[erf((nu + nu_c)/10) - erf((nu - nu_c)/10)],  nu_c = Lam + 100   (UNIFIED_THEOREM_v3 §2.1),
  E_theta(mu) = 2 Re[theta(mu - i) f(mu)]   (realised weight, TAIL_CERT §R.2 form),
  G_A^real(z) = (2/Lam) int_0^oo E_theta(mu) F_{3,mu}(z) dmu   (FULLPOS3 §4.1(a)).
E_theta is real-analytic on mu > 0 (continuation E(mu) = theta(mu-i) f(mu) + conj(theta(conj(mu)-i) f(conj(mu)))).
"""
import os, sys, bisect
sys.dont_write_bytecode = True
from flint import arb, acb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'epswcert'))
_p = ctx.prec
import ew_common as EW                      # sets ctx.prec = 192 and ctx.cap at import
ctx.prec = _p
REP = os.path.join(ADS4, 'code', 'epswcert', 'cache', 'k2_continuum_D_1000_rep.json')
I = acb(0, 1)


class RealWeight:
    def __init__(self, far=60.0):
        old = ctx.prec
        self.W = EW.Weight(1000, path=REP)
        assert self.W.json['repair']['changed'] == 'x_re[0]'
        self.W.CUT = far                        # instance override: near window |mu - c_k| <= far * sigma_k for the disc sups
        self.far = far
        self.sc = self.W.sc                     # (signed centre, sigma, x) floats, sorted
        self.cs = self.W.cs_sorted
        ctx.prec = old

    def _terms(self, lo, hi):
        i0 = bisect.bisect_left(self.cs, lo); i1 = bisect.bisect_right(self.cs, hi)
        return self.sc[i0:i1]

    def far_abs(self, mu, D):
        """upper bound of sum over centres with |mu - c| > D of |x| |exp(-(mu - i - c)^2/(2 s^2))| for real mu:
        each term <= |x| exp((1 - D^2)/(2 s^2))."""
        tot = arb(0)
        for s, sx in self.W.sumabs.items():
            sa = arb(s)
            tot += sx*((1 - arb(D)**2)/(2*sa*sa)).exp()
        return tot

    def f_real(self, mu):
        """f(mu) for real mu (arb ball, mu >= 0.2): all RBF terms with |mu - c| <= 60 exactly, the rest as a radius."""
        D = 60.0
        m = float(mu.mid())
        nu = acb(mu, -1)
        Bs = acb(0)
        for c, s, x in self._terms(m - D, m + D):
            Bs += arb(x)*(-((nu - arb(c))**2)/(2*arb(s)**2)).exp()
        Bs += acb(arb(0, self.far_abs(mu, D).upper()), arb(0, self.far_abs(mu, D).upper()))
        QR = EW.horner(EW.QRN, nu)/EW.horner(EW.QRD, nu)
        W = self.W
        nu2 = nu*nu
        Y = (W.y[0] + W.y[1]*nu2 + W.y[2]*nu2*nu2)*((-nu2/2).exp())
        L1 = nu2 + 1
        Rm = EW.horner(EW.RMN, nu)/EW.horner(EW.RMD, nu)
        return QR*Bs + Y*L1*L1.log()*Rm

    def theta(self, zeta):
        nuc = self.W.Lam + 100
        return ((zeta + nuc)/10).erf()/2 - ((zeta - nuc)/10).erf()/2

    def E_real(self, mu):
        """E_theta(mu) = 2 Re[theta(mu - i) f(mu)] (arb)."""
        v = self.theta(acb(mu, -1))*self.f_real(mu)
        return 2*v.real

    def E_disc_sup(self, m, rho):
        """upper bound of sup |E_theta| over the complex disc |mu - m| <= rho (m >= rho + 0.1):
        |E| <= |theta(mu-i) f(mu)| + |theta(conj mu - i) f(conj mu)| <= 2 sup_disc|theta(.-i)| sup_disc|f| (disc is conj-symmetric);
        sup|f| from ew_common.f_disc_sup (explicit RBF bounds, ring boxes, maximum modulus);
        sup|theta(mu - i)| over the disc: maximum modulus on a ring of boxes (theta entire)."""
        old = ctx.prec
        ctx.prec = 128
        try:
            M, F = EW.f_disc_sup(self.W, m, rho, nb=64)
            th = arb(0)
            for z in EW.ring_boxes(m, rho, 64):
                v = abs(self.theta(z))
                if v.upper() > th.upper():
                    th = arb(v.upper())
            return arb((2*th*M).upper())
        finally:
            ctx.prec = old
