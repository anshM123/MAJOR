"""kc_common.py -- ball-arithmetic (python-flint / Arb) evaluators for the K certificate (ads4/K_CERT.md).

Everything here is read from published data, read-only:
  * prl/data/nf02_psi.json                 NF02 psi2, psi3 (exact rationals; divided by k0 = 27.48509817, an exact decimal)
  * prl/data/kernels.json  -> D_k2_EE      R_-(nu), R_+(nu) (exact rational functions with Gaussian-rational coefficients)
  * prl/data/k2/k2_continuum_D_<Dg>.json   the continuum k = 2 weight chi = chi_re + i chi_im (IEEE doubles, exact as given)
The functional's k = 2 weight is X_theta = theta * chi * (R_- + R_+), theta = (erf((nu+nu_c)/a) - erf((nu-nu_c)/a))/2,
nu_c = Lam + 100, a = 10 (Supplement Sec. S2).  Only chi_re enters any real light value (R_- + R_+ is real on the real axis and the
functional takes Re), so chi_im is never needed.

Evaluation of chi_re on a complex ball z: the Gaussian sum is truncated to the centres within D = 10 of the ball's real range and the
omitted centres are added as a rigorous error ball (bound in gauss_tail()).  l(nu) = (nu^2+1) log(nu^2+1) is evaluated with an explicit
analyticity check (principal branch; cut where nu^2+1 <= 0).
"""
import os, json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from flint import arb, acb, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
K0_DEC = Fr('27.48509817')
NAN = acb(arb('nan'), arb('nan'))


def lam_fmpq_sq(Dg):
    return fmpq((Dg - 2) * (Dg - 1))


def Lam(Dg):
    return arb(lam_fmpq_sq(Dg)).sqrt()


# ---------------------------------------------------------------- flat weights (exact)
def nf02():
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    assert Fr(d['k0']) == K0_DEC
    ld = lambda v: v if isinstance(v, list) else eval(v)
    p2 = [Fr(s) / K0_DEC for s in ld(d['psi2_exact'])]
    p3 = [Fr(s) / K0_DEC for s in ld(d['psi3_exact'])]
    return p2, p3


def poly_moment(c, j):
    """exact int_0^1 p^j sum_n c_n p^n dp"""
    return sum(cn / (n + j + 1) for n, cn in enumerate(c))


def K0_exact():
    p2, p3 = nf02()
    return -2 * poly_moment(p2, 6) + 2 * poly_moment(p3, 4)


def fq(x):
    x = Fr(x)
    return fmpq(x.numerator, x.denominator)


# ---------------------------------------------------------------- rational kernels (exact, Gaussian rationals)
nu_s = sp.Symbol('nu')


class Rat:
    def __init__(self, expr):
        expr = sp.cancel(sp.together(sp.sympify(expr, locals={'nu': nu_s, 'I': sp.I})))
        num, den = sp.fraction(expr)
        self.expr = expr
        self.num = [sp.nsimplify(sp.expand(c)) for c in sp.Poly(sp.expand(num), nu_s).all_coeffs()]
        self.den = [sp.nsimplify(sp.expand(c)) for c in sp.Poly(sp.expand(den), nu_s).all_coeffs()]
        self._prec = None

    def _co(self):
        if self._prec != ctx.prec:
            cv = lambda c: acb(arb(fq(Fr(str(sp.re(c))))), arb(fq(Fr(str(sp.im(c))))))
            self.an = [cv(c) for c in self.num]; self.ad = [cv(c) for c in self.den]
            self._prec = ctx.prec
        return self.an, self.ad

    def __call__(self, z):
        an, ad = self._co()
        n = acb(0); d = acb(0)
        for c in an:
            n = n * z + c
        for c in ad:
            d = d * z + c
        return n / d


def kernels():
    d = json.load(open(os.path.join(DATA, 'kernels.json')))['D_k2_EE']
    return d['R_minus'], d['R_plus']


# ---------------------------------------------------------------- continuum weight chi_re
class Chi:
    D = 10.0           # truncation radius of the Gaussian sum (in nu)
    SIG_MAX = 0.9
    H_MIN = 0.15

    def __init__(self, Dg, part='re'):
        d = json.load(open(os.path.join(DATA, 'k2', 'k2_continuum_D_%d.json' % Dg)))
        assert d['Delta_g'] == Dg
        cs = np.array(d['centres_sigma'], float)
        self.c = cs[:, 0]; self.s = cs[:, 1]
        assert np.all(np.diff(self.c) > 0) and np.all(np.diff(self.c) >= self.H_MIN - 1e-12) and self.s.max() <= self.SIG_MAX
        self.x = np.array([float(v) for v in d['x_' + part]])
        self.y = [float(v) for v in d['y_' + part]]
        self.xmax = float(np.abs(self.x).max())
        self.d = d
        self._prec = None

    def _co(self):
        if self._prec != ctx.prec:
            self.ac = [arb(float(v)) for v in self.c]
            self.a2s2 = [1 / (2 * arb(float(v)) ** 2) for v in self.s]
            self.ax = [arb(float(v)) for v in self.x]
            self.ay = [arb(float(v)) for v in self.y]
            self._prec = ctx.prec
        return self.ac, self.a2s2, self.ax, self.ay

    def gauss_tail(self, Y):
        """bound on sum over omitted centres of |x_k| (|g(z-c)| + |g(z+c)|) for |Im z| <= Y, distance >= D from Re z (both sides)."""
        D, s, h = arb(self.D) - arb('1e-6'), arb(self.SIG_MAX), arb(self.H_MIN) - arb('1e-9')   # margins cover float rounding of the window
        Y = arb(Y)
        q = (-(D * h) / (s * s)).exp()
        return 4 * arb(self.xmax) * (-(D * D - Y * Y) / (2 * s * s)).exp() / (1 - q)

    @staticmethod
    def logw(w, analytic):
        """principal log of w; non-finite if the ball may touch the cut (-oo, 0] and analyticity is required."""
        if analytic or True:
            re = w.real
            if not (re > 0):
                im = w.imag
                if not (im > 0 or im < 0):
                    return NAN
        return w.log()

    def B(self, z):
        """Gaussian sum B(z) (windowed) as an acb ball including the truncation error."""
        ac, a2s2, ax, ay = self._co()
        zr = z.real; m = float(zr.mid()); r = float(zr.rad())
        lo, hi = m - r - self.D, m + r + self.D
        i0, i1 = np.searchsorted(self.c, lo), np.searchsorted(self.c, hi, side='right')
        tot = acb(0)
        for k in range(i0, i1):
            w = z - ac[k]
            tot += ax[k] * (-(w * w) * a2s2[k]).exp()
        # mirrored terms g(z + c), c > 0, centres with c <= -lo
        j1 = np.searchsorted(self.c, -lo, side='right')
        for k in range(0, j1):
            if self.c[k] > 0:
                w = z + ac[k]
                tot += ax[k] * (-(w * w) * a2s2[k]).exp()
        Y = float(z.imag.abs_upper())
        err = self.gauss_tail(Y)
        return tot + acb(arb(0, err), arb(0, err))

    def sing(self, z, analytic):
        """sum_j y_j z^(2j) l(z) e^{-z^2/2},  l = (z^2+1) log(z^2+1)."""
        ac, a2s2, ax, ay = self._co()
        w = z * z + 1
        L = self.logw(w, analytic)
        if not L.is_finite():
            return NAN
        e = (-(z * z) / 2).exp()
        z2 = z * z
        return (ay[0] + z2 * (ay[1] + z2 * ay[2])) * w * L * e

    def sing_over_z(self, z, analytic):
        """sing(z)/z, analytic at 0: uses log(1+u)/u series for |z| < 1/2 (u = z^2)."""
        ac, a2s2, ax, ay = self._co()
        if float(abs(z).upper()) < 0.5:
            u = z * z
            ua = float(abs(u).upper())
            N = 40
            s = acb(0)
            for n in range(N - 1, -1, -1):
                s = s * (-u) + acb(1) / (n + 1)
            rem = arb(ua) ** N / ((N + 1) * (1 - arb(ua)))
            s += acb(arb(0, rem), arb(0, rem))
            Lz = z * s                     # log(1+z^2)/z
            e = (-(z * z) / 2).exp()
            return (ay[0] + u * (ay[1] + u * ay[2])) * (u + 1) * Lz * e
        return self.sing(z, analytic) / z

    @staticmethod
    def Q(z):
        z2 = z * z
        return z2 / ((z2 + acb(1) / 4) ** 2)

    @staticmethod
    def Q_over_z(z):
        z2 = z * z
        return z / ((z2 + acb(1) / 4) ** 2)

    def __call__(self, z, analytic=True):
        s = self.sing(z, analytic)
        if not s.is_finite():
            return NAN
        return self.Q(z) * self.B(z) + s

    def over_z(self, z, analytic=True):
        s = self.sing_over_z(z, analytic)
        if not s.is_finite():
            return NAN
        return self.Q_over_z(z) * self.B(z) + s


def theta(z, Lam_):
    nc = Lam_ + 100; a = arb(10)
    return (((z + nc) / a).erf() - ((z - nc) / a).erf()) / 2


def psi_poly_acb(coeffs, z, Lam_):
    """sum_n c_n (z/Lam)^n  (Horner), coeffs exact Fractions."""
    p = z / Lam_
    s = acb(0)
    for c in reversed(coeffs):
        s = s * p + arb(fq(c))
    return s
