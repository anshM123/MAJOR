"""Common tools for ADS4_FUNCTIONAL.md (AdS4 uplift of the flat graviton functional F_{g3}).

Floating point throughout: nothing here is a certificate.

Notation (AdS units R = 1):
  heavy state (tau, J), tau >= Delta_gap (twist gap), alpha = tau - 2, alphabar = tau + 2J - 1  (d = 3),
  m^2 = alpha*alphabar,  cosh(beta) = (alpha^2 + alphabar^2)/(2 alpha alphabar)   [CHMRS (2.28)-(2.30)],
  smearing scale Lam (nu in [0, Lam]);  flat variables  m_fl = m/Lam, u = 1/m_fl^2, bhat = (2J+1)/m_fl.

Weights (FLAT_GRAVITON_BOUND.md (0.1)): psi2 = (1-p)^3 p (65+155p+47p^2), psi3 = 2(1-p)^3(81-1927p+4601p^2).

Transforms of a polynomial weight w(p) = sum_i c_i p^i:
  cosT(w, tau)        = int_0^1 w(p) cos(p tau) dp                              (closed form / Taylor)
  flatT(w, b)         = int_0^1 w(p) J_0(b p) dp           = (2/pi) int_0^{pi/2} cosT(w, b sin th) dth
  curvT(w, Lam, beta) = int_0^1 w(p) P_{-1/2+i Lam p}(cosh beta) dp
                      = (2/pi) int_0^{pi/2} cosT(w, Lam beta sin th) E(beta, beta sin th)^{-1/2} dth
  (Mehler-Dirichlet in the form UNIFORM_POSITIVITY.md (2.4); E = (sinh a/a)(sinh c/c), a,c = (beta +- t)/2).
  W(beta)  := -d/d(nu^2) P_{-1/2+i nu}(cosh beta)|_{nu=0} = (2/pi) int (beta sin th)^2/2 E^{-1/2} dth
  P0(beta) := P_{-1/2}(cosh beta) = (2/pi) int E^{-1/2} dth
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from math import factorial

# ------------------------------------------------------------------ polynomials (float coefficients, index = power)
def pmul(a, b):
    r = np.zeros(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        r[i:i + len(b)] += x * np.asarray(b, float)
    return r

OMP3 = np.array([1., -3., 3., -1.])
PSI2 = pmul(OMP3, [0., 65., 155., 47.])
PSI3 = pmul(OMP3, [162., -3854., 9202.])
P2PSI2 = np.concatenate([[0., 0.], PSI2])
BASIS = {'psi2': PSI2, 'p2psi2': P2PSI2, 'psi3': PSI3}

def peval(c, x):
    r = np.zeros_like(np.asarray(x, float))
    for ci in c[::-1]:
        r = r * x + ci
    return r

def pder_vals(c, x):
    """[w(x), w'(x), ..., w^(deg)(x)]"""
    out = []
    a = np.array(c, float)
    for _ in range(len(c)):
        out.append(peval(a, x))
        a = np.array([i * a[i] for i in range(1, len(a))]) if len(a) > 1 else np.array([0.])
    return np.array(out)

def moments(c, nmax):
    return np.array([sum(ci / (i + n + 1) for i, ci in enumerate(c)) for n in range(nmax + 1)])

class CosT:
    """int_0^1 w(p) cos(p tau) dp for a fixed polynomial w, vectorised in tau."""
    TAU0 = 8.0
    K = 40
    def __init__(self, c):
        self.c = np.array(c, float)
        self.d0 = pder_vals(self.c, 0.0)
        self.d1 = pder_vals(self.c, 1.0)
        mom = moments(self.c, 2 * self.K)
        self.tay = np.array([(-1) ** k * mom[2 * k] / factorial(2 * k) for k in range(self.K)])
    def __call__(self, tau):
        tau = np.atleast_1d(np.asarray(tau, float))
        out = np.empty_like(tau)
        sm = np.abs(tau) < self.TAU0
        t2 = tau[sm] ** 2
        r = np.zeros_like(t2)
        for a in self.tay[::-1]:
            r = r * t2 + a
        out[sm] = r
        tl = tau[~sm]
        if tl.size:
            z = -1j * tl
            ez = np.exp(-z)
            acc = np.zeros_like(z)
            zp = 1.0 / z
            for j in range(len(self.c)):
                acc += (self.d0[j] - ez * self.d1[j]) * zp
                zp = zp / z
            out[~sm] = acc.real
        return out

COS = {k: CosT(v) for k, v in BASIS.items()}

_GL = {}
def gl(n):
    if n not in _GL:
        x, w = leggauss(n)
        _GL[n] = ((x + 1) * np.pi / 4, w * np.pi / 4)      # nodes/weights on [0, pi/2]
    return _GL[n]

def nth(b):
    n = int(min(12000, max(96, 1.2 * b + 96)))
    return 128 * ((n + 127) // 128)          # quantised so that the Gauss-Legendre rules are cached

def Einvsqrt(beta, s):
    """E(beta, beta s)^{-1/2}, s in [0,1]."""
    a = beta * (1 + s) / 2
    c = beta * (1 - s) / 2
    def r(x):
        x = np.asarray(x, float)
        out = np.ones_like(x)
        nz = x > 1e-8
        # log(sinh x / x) stably
        xx = x[nz]
        out[nz] = np.exp(-0.5 * (xx + np.log1p(-np.exp(-2 * xx)) - np.log(2 * xx)))
        out[~nz] = 1 - x[~nz] ** 2 / 12
        return out
    return r(a) * r(c)

def curvT_all(Lam, beta):
    """dict name -> int_0^1 w P_{-1/2+i Lam p}(cosh beta) dp, plus W(beta), P0(beta)."""
    b = Lam * beta
    th, wt = gl(nth(b))
    s = np.sin(th)
    E = Einvsqrt(beta, s)
    out = {k: (2 / np.pi) * np.sum(wt * f(b * s) * E) for k, f in COS.items()}
    out['W'] = (2 / np.pi) * np.sum(wt * (beta * s) ** 2 / 2 * E)
    out['P0'] = (2 / np.pi) * np.sum(wt * E)
    return out

def flatT_all(b):
    th, wt = gl(nth(b))
    s = np.sin(th)
    return {k: (2 / np.pi) * np.sum(wt * f(b * s)) for k, f in COS.items()}

# ------------------------------------------------------------------ flat partial waves (vector over J)
def jac_all(nmax, a, b, x):
    P = np.zeros((nmax + 1, len(x)))
    P[0] = 1
    if nmax >= 1:
        P[1] = (a + 1) + (a + b + 2) * (x - 1) / 2
    for n in range(2, nmax + 1):
        c = 2 * n + a + b
        P[n] = (((c - 1) * (a * a - b * b) + (c - 2) * (c - 1) * c * x) * P[n - 1]
                - 2 * (n + a - 1) * (n + b - 1) * c * P[n - 2]) / (2 * n * (n + a + b) * (c - 2))
    return P

SIGMA = {'++': 1, '+-': -1}
PHI_CONST = {'++': 0, '+-': 210}       # completed (g5-free) functional: Phi = 10 J(J+1) - PHI_CONST

class FlatGrid:
    """Gauss-Legendre in p for the flat heavy action at fixed u, all J <= Jmax at once."""
    def __init__(self, n=1600):
        x, w = leggauss(n)
        self.p = (x + 1) / 2
        self.w = w / 2
        p = self.p
        self.psi2 = peval(PSI2, p)
        self.psi3 = peval(PSI3, p)

    def D(self, h, Jmax, u):
        x = 1 - 2 * u * self.p ** 2
        if h == '++':
            return np.arange(Jmax + 1), jac_all(Jmax, 0, 0, x)
        return np.arange(4, Jmax + 1), jac_all(Jmax - 4, 0, 8, x)

    def parts(self, h, Jmax, u):
        """J array and the flat parts (per unit, normal form m^4 F):
        B2 = int psi2 (2 - u p^2) D,  B3 = sigma u int psi3 D,  FW = Phi u^3, and the absolute norm int|Psi||D| + |Phi|u^3."""
        Js, D = self.D(h, Jmax, u)
        s = SIGMA[h]
        p2 = self.p ** 2
        w2 = self.w * self.psi2 * (2 - u * p2)
        w3 = self.w * s * u * self.psi3
        B2 = D @ w2
        B3 = D @ w3
        Phi = 10.0 * Js * (Js + 1) - PHI_CONST[h]
        FW = Phi * u ** 3
        ab = np.abs(D) @ np.abs(w2 + w3) + np.abs(FW)
        return Js, B2, B3, FW, ab

# ------------------------------------------------------------------ AdS kinematics
def kin(tau, J):
    al = tau - 2.0
    ab = tau + 2.0 * J - 1.0
    m2 = al * ab
    ch = (al * al + ab * ab) / (2 * al * ab)
    return np.sqrt(m2), np.arccosh(ch)

def tau_of_m(m, J):
    """twist with given m and J (d=3): (tau-2)(tau+2J-1) = m^2."""
    B = 2 * J - 3.0
    C = -2 * (2 * J - 1.0) - m * m
    return (-B + np.sqrt(B * B - 4 * C)) / 2

def lam_of_gap(Dg):
    """smallest heavy m: J = 0, tau = Delta_gap."""
    return np.sqrt((Dg - 2.0) * (Dg - 1.0))

# ------------------------------------------------------------------ AdS actions (normal form: m_fl^4 x action)
def regge_parts(h, Lam, m, beta):
    """Leading Regge-limit AdS action of Omega_g3 (conditional heavy formula C^pm = 2P/m^{2k}(...)), split in components.
    B2 = 2 int psi2 P,  B3 = sigma u int psi3 P,  B4 = 10 W(beta) Lam^2 u^2,  B5 = 5 sigma u^3 P0(beta)."""
    c = curvT_all(Lam, beta)
    u = (Lam / m) ** 2
    s = SIGMA[h]
    B2 = 2 * c['psi2']
    B3 = s * u * c['psi3']
    B4 = 10 * c['W'] * Lam ** 2 * u ** 2
    B5 = 5 * s * u ** 3 * c['P0']
    return B2, B3, B4, B5, c

def hybrid_shift(h, Lam, m, beta, J):
    """Curvature shift of the hybrid model (flat + [Regge - its flat limit]) relative to the flat action, normal form.
    kernel K = D_h(1-2up^2) + P_{-1/2+i Lam p}(cosh beta) - J_0(bhat p);
    forward: Phi^hyb = Phi - 10 JJ - 5/2 + 10 m^2 W(beta) + 5(1 - sigma)(1 - P0(beta))."""
    c = curvT_all(Lam, beta)
    mfl = m / Lam
    u = 1 / mfl ** 2
    bh = (2 * J + 1) / mfl
    f = flatT_all(bh)
    s = SIGMA[h]
    d = {k: c[k] - f[k] for k in BASIS}
    dB2 = 2 * d['psi2'] - u * d['p2psi2']
    dB3 = s * u * d['psi3']
    dFW = u ** 3 * (-2.5 + 10 * m * m * c['W'] - 10 * J * (J + 1) + 5 * (1 - s) * (1 - c['P0']))
    return dB2, dB3, dFW
