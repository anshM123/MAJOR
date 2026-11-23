"""
lrtight/lrt_common.py (LRINF_TIGHT.md): shared exact / ball-arithmetic helpers for the certification of Lemma LR-inf
(LARGESPIN_L.md §3) with the tight constant int|psi|.

Inputs: ads4/prl/data/nf02_psi.json only (exact rational coefficients; k0 = 27.48509817 is an exact decimal).
All bound quantities are evaluated in Arb ball arithmetic (python-flint), precision 128 bits; every function below returns an arb
ball that CONTAINS the exact value (or an upper bound where stated).

Notation (LARGESPIN_L §3): psi = psi_k (k = 2, 3), a1 = psi'(0), a2 = psi''(0), g = psi/p, Lam = sqrt(998*999) (Delta_g = 1000),
  R_inf(x) = sum_{j>=3 odd} |psi^(j)(0)|/x^{j+1} + sum_{j>=5} |psi^(j)(1)|/x^{j+1}            (>= |r(x)|, r = Phi + a1/x^2)
  R(x)     = min(R_inf(x), a1/x^2 + int|psi|)
  d1(X)    = sum_{j>=2 even} |g^(j)(0)|/X^{j+1} + sum_{j>=5} |g^(j)(1)|/X^{j+1}
  sumc(X)  = int_X^inf x^2 R_inf(x) dx = sum_j c_j X^{2-j}/(j-2)
No project module is imported; nothing outside code/lrtight/cache is written.
"""
import os, sys, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, ctx, fmpq, fmpq_poly

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'prl', 'data'))

TAU = 1000                       # Delta_g
LAM2 = Fr(998*999)               # Lam^2 exact


def A(q):
    """exact rational -> arb (exact if representable, else a tight ball)."""
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


LAM = A(LAM2).sqrt()


def psi_coeffs(k):
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    k0 = Fr(d['k0'])
    return [Fr(c)/k0 for c in d['psi%d_exact' % k]]


def derivs(c, p0):
    """exact c^(j)(p0) for a coefficient list c (lowest power first)."""
    D = len(c) - 1; out = []
    for j in range(D + 1):
        s = Fr(0)
        for n in range(j, D + 1):
            ff = 1
            for i in range(j):
                ff *= (n - i)
            s += c[n]*ff*Fr(p0)**(n - j)
        out.append(s)
    return out


class Weight:
    """exact data of psi_k and the ball-arithmetic bound functions of Lemma LR-inf."""

    def __init__(self, k, absint_up=None):
        self.k = k
        self.c = psi_coeffs(k)
        assert self.c[0] == 0
        self.d0 = derivs(self.c, 0); self.d1 = derivs(self.c, 1)
        assert all(self.d1[j] == 0 for j in range(5)), '(1-p)^5 | psi'
        self.g = self.c[1:]
        self.g0 = derivs(self.g, 0); self.g1 = derivs(self.g, 1)
        assert all(self.g1[j] == 0 for j in range(5))
        self.a1 = self.d0[1]; self.a2 = self.d0[2]
        n = len(self.c)
        # R_inf coefficients: list of (power p, |coef|) with term |coef|/x^p
        self.Rterms = [(j + 1, abs(self.d0[j])) for j in range(3, n, 2) if self.d0[j] != 0] + \
                      [(j + 1, abs(self.d1[j])) for j in range(5, n) if self.d1[j] != 0]
        self.dterms = [(j + 1, abs(self.g0[j])) for j in range(2, len(self.g), 2) if self.g0[j] != 0] + \
                      [(j + 1, abs(self.g1[j])) for j in range(5, len(self.g)) if self.g1[j] != 0]
        self.absint_up = absint_up        # rigorous upper bound of int_0^1 |psi| (set by lrt1)

    # ----- ball evaluations (x, X > 0 balls) -----
    def Rinf(self, x):
        s = arb(0)
        for p, cf in self.Rterms:
            s += A(cf)/x**p
        return s

    def R(self, x):
        """upper bound of |r(x)| (ball containing min(R_inf, a1/x^2 + absint_up))."""
        b1 = self.Rinf(x)
        b2 = A(self.a1)/x**2 + A(self.absint_up)
        # min of two balls: an enclosure of the min is fine for an upper bound we take .upper() later
        return b1 if b1.upper() <= b2.upper() else b2

    def d1f(self, X):
        s = arb(0)
        for p, cf in self.dterms:
            s += A(cf)/X**p
        return s

    def sumc(self, X):
        s = arb(0)
        for p, cf in self.Rterms:        # term cf x^{2-p}, p = j+1 >= 4: int_X^inf = cf X^{3-p}/(p-3)
            s += A(cf)*X**(3 - p)/(p - 3)
        return s


def up(x):
    """upper endpoint as an exact arb."""
    return arb(x.upper())


def lo(x):
    return arb(x.lower())


def ball(a, b):
    """arb ball containing the rational interval [a, b]."""
    a, b = Fr(a), Fr(b)
    return A(a).union(A(b))


def P_mhalf(beta):
    """P_{-1/2}(cosh beta) = 2F1(1/2, 1/2; 1; -sinh^2(beta/2))  (decreasing in beta: Euler integral)."""
    z = -(beta/2).sinh()**2
    return z.hypgeom_2f1(arb(1)/2, arb(1)/2, arb(1))


def mgap2(J, tau=TAU):
    """m^2 of the twist-edge state Delta = tau + J (CHMRS d = 3): (Delta-J-2)(Delta+J-1) - exact integer."""
    return (tau - 2)*(tau + 2*J - 1)
