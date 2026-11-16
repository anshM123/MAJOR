"""Common tools for ROBUST_FUNCTIONAL.md: flat graviton functionals for |ghat_3|^2 with general polynomial weights.

Functional (flat, units M = 1; generalises F_{g3} of arXiv:2201.06602 (3.3a), g5-free completion):

    F = int_0^1 dp [psi2(p) B_2^{(1)}(p^2) + psi3(p) B_3^{(1)}(p^2)]  - a d_{p^2} B_4^{(1)}|_0 + a B_5^{(1)}|_0 ,

    psi2(p) = (1-p)^3 p  sum_i x_i p^i ,   psi3(p) = (1-p)^3 sum_j y_j p^j ,   a >= 0 .

(F_{g3}: x = (65,155,47), y = (162,-3854,9202), a = 5.)  The coefficients a of -d B4|_0 and +B5|_0 are tied so that g5 cancels.

Low side (FLAT_GRAVITON_BOUND.md section 1.3, cert_low.py):   F|_high = 16 pi G L(eps) - pi G k |ghat3|^2 + (matter),
    k = -2 mu2_6 + 2 mu3_4 + 4 a,   L(eps) = int_eps^1 psi2/p^2 = x_0 log(1/eps) + Lfin + O(eps),
so (heavy positivity)  |ghat3|^2 M^8 <= A log(M/m_IR) + B,  A = 16 x_0/k,  B = 16 Lfin/k   (eps -> 0 constants).
AdS (ADS4_LIGHT.md): L_AdS = x_0 log(2 Lam/3) + Lfin + O(1/Lam), so B_AdS = B - A log(3/2).

Heavy action per |c^h_{J,m^2}|^2 (normal form FF = m^4 F, u = 1/m^2), FLAT_GRAVITON_BOUND.md (0.1) generalised:
    FF_h = B2 + B3 + FW,
    B2 = int psi2 (2 - u p^2) D_h(1-2up^2),  B3 = sigma_h u int psi3 D_h,  FW = a (2 JJ - c4_h + sigma_h) u^3,
    c4_{++} = 1, c4_{+-} = 41 (from d_{p^2} C_4|_0 = (c4 - 2JJ)/m^10), JJ = J(J+1).
Robustness reference: N = |B2| + |B3| + |FW| (the component norm of ADS4_FUNCTIONAL.md, N_comp).
Floating point (numpy); nothing here is a certificate.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from math import factorial, log

SIGMA = {'++': 1, '+-': -1}
C4 = {'++': 1, '+-': 41}

# ------------------------------------------------------------------ polynomials (index = power of p)
def pmul(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    r = np.zeros(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        r[i:i + len(b)] += x * b
    return r

def peval(c, x):
    x = np.asarray(x, float)
    r = np.zeros_like(x)
    for ci in np.asarray(c, float)[::-1]:
        r = r * x + ci
    return r

def mom(c, j):
    return sum(ci / (i + j + 1) for i, ci in enumerate(c))

OMP3 = np.array([1., -3., 3., -1.])

def basis2(d2, r=3):
    """psi2 basis: (1-p)^r p^{i+1}, i = 0..d2."""
    omp = np.array([1.])
    for _ in range(r):
        omp = pmul(omp, [1., -1.])
    return [pmul(omp, np.eye(i + 2)[i + 1]) for i in range(d2 + 1)]

def basis3(d3, r=3):
    """psi3 basis: (1-p)^r p^j, j = 0..d3."""
    omp = np.array([1.])
    for _ in range(r):
        omp = pmul(omp, [1., -1.])
    return [pmul(omp, np.eye(j + 1)[j]) for j in range(d3 + 1)]

class Functional:
    """psi2, psi3 (monomial coefficients) and forward coefficient a."""
    def __init__(self, psi2, psi3, a, name=''):
        self.psi2 = np.trim_zeros(np.asarray(psi2, float), 'b')
        self.psi3 = np.trim_zeros(np.asarray(psi3, float), 'b')
        self.a = float(a)
        self.name = name

    @staticmethod
    def from_vars(v, d2, d3, r=3, name=''):
        b2 = basis2(d2, r); b3 = basis3(d3, r)
        n = max(len(b) for b in b2 + b3)
        p2 = np.zeros(n); p3 = np.zeros(n)
        for i, b in enumerate(b2): p2[:len(b)] += v[i] * b
        for j, b in enumerate(b3): p3[:len(b)] += v[d2 + 1 + j] * b
        return Functional(p2, p3, v[d2 + d3 + 2], name)

    def k(self):
        return -2 * mom(self.psi2, 6) + 2 * mom(self.psi3, 4) + 4 * self.a

    def x0(self):
        return self.psi2[1]

    def Lfin(self):
        c = self.psi2.copy()
        # int_0^1 (psi2 - x0 p)/p^2
        c[1] -= self.psi2[1]
        assert abs(c[0]) < 1e-12
        return sum(c[i] / (i - 1) for i in range(2, len(c)))

    def AB(self):
        k = self.k()
        A = 16 * self.x0() / k
        B = 16 * self.Lfin() / k
        return A, B, B - A * log(1.5)

    def normalised(self):
        k = self.k()
        return Functional(self.psi2 / k, self.psi3 / k, self.a / k, self.name)

def paper_fg3():
    psi2 = pmul(OMP3, [0., 65., 155., 47.])
    psi3 = pmul(OMP3, [162., -3854., 9202.])
    return Functional(psi2, psi3, 5.0, 'F_g3 (g5-free)')

# ------------------------------------------------------------------ flat partial waves: running dot products
def gl01(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2

def flat_integrals(h, u, Jmax, W, p):
    """Matrix R[J_index, c] = sum_nodes D_h^J(1-2 u p^2) W[:, c] for J = J0..Jmax (J0 = 0 for ++, 4 for +-).
    W: (n, ncols) node-weighted functions (already multiplied by the Gauss weights).  Returns (Js, R)."""
    x = 1 - 2 * u * p * p
    if h == '++':
        a, b, J0 = 0, 0, 0
    else:
        a, b, J0 = 0, 8, 4
    nmax = Jmax - J0
    R = np.zeros((nmax + 1, W.shape[1]))
    P0 = np.ones_like(x)
    R[0] = P0 @ W
    if nmax >= 1:
        P1 = (a + 1) + (a + b + 2) * (x - 1) / 2
        R[1] = P1 @ W
        Pm, Pc = P0, P1
        for k in range(2, nmax + 1):
            c = 2 * k + a + b
            Pn = (((c - 1) * (a * a - b * b) + (c - 2) * (c - 1) * c * x) * Pc
                  - 2 * (k + a - 1) * (k + b - 1) * c * Pm) / (2 * k * (k + a + b) * (c - 2))
            R[k] = Pn @ W
            Pm, Pc = Pc, Pn
    return np.arange(J0, Jmax + 1), R

_NODES = {}
def nodes(n):
    if n not in _NODES:
        _NODES[n] = gl01(n)
    return _NODES[n]

def nquad(Jmax):
    """exact Gauss-Legendre for degree 2 Jmax + 20 in p."""
    return int(Jmax + 16)

def flat_parts_basis(h, u, Jmax, polys2, polys3):
    """For a list of psi2 basis polynomials and psi3 basis polynomials: per J,
    B2 columns (int phi (2 - u p^2) D) and B3 columns (sigma u int phi D).  Returns Js, B2 (nJ x n2), B3 (nJ x n3)."""
    n = nquad(Jmax)
    p, w = nodes(n)
    cols = []
    for c in polys2:
        cols.append(w * peval(c, p) * (2 - u * p * p))
    for c in polys3:
        cols.append(w * SIGMA[h] * u * peval(c, p))
    W = np.array(cols).T
    Js, R = flat_integrals(h, u, Jmax, W, p)
    return Js, R[:, :len(polys2)], R[:, len(polys2):]

def phi_coef(h, Js, u):
    """FW per unit a."""
    JJ = Js * (Js + 1.0)
    return (2 * JJ - C4[h] + SIGMA[h]) * u ** 3

def flat_parts(F, h, u, Jmax):
    """Flat parts (B2, B3, FW) of a given Functional for J <= Jmax at fixed u."""
    Js, B2, B3 = flat_parts_basis(h, u, Jmax, [F.psi2], [F.psi3])
    FW = F.a * phi_coef(h, Js, u)
    return Js, B2[:, 0], B3[:, 0], FW

def jmax_for(u, bmax, Jcap):
    return int(min(Jcap, max(40, np.ceil((bmax / np.sqrt(u) - 1) / 2))))

# ------------------------------------------------------------------ scaling limit (u -> 0, bhat fixed): 2 int psi2 J_0(bhat p)
def scaling_B2(polys2, bh, n=3000):
    from scipy.special import j0
    p, w = gl01(n)
    J = j0(np.outer(bh, p))
    return np.array([2 * (J @ (w * peval(c, p))) for c in polys2]).T

# ------------------------------------------------------------------ AdS: Mehler-Dirichlet transforms of arbitrary polynomial weights
class CosT:
    """int_0^1 w(p) cos(p tau) dp for a fixed polynomial w (Taylor for |tau| < 8, endpoint expansion beyond)."""
    TAU0 = 8.0
    K = 48
    def __init__(self, c):
        self.c = np.array(c, float)
        d = len(self.c)
        self.d0 = []; self.d1 = []
        a = self.c.copy()
        for _ in range(d):
            self.d0.append(peval(a, 0.0)); self.d1.append(peval(a, 1.0))
            a = np.array([i * a[i] for i in range(1, len(a))]) if len(a) > 1 else np.array([0.])
        m = [mom(self.c, n) for n in range(2 * self.K + 1)]
        self.tay = np.array([(-1) ** k * m[2 * k] / factorial(2 * k) for k in range(self.K)])
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

_GLA = {}
def glang(n):
    if n not in _GLA:
        if n <= 512:
            x, w = leggauss(n)
            _GLA[n] = ((x + 1) * np.pi / 4, w * np.pi / 4)
        else:
            K = (n + 255) // 256
            x, w = leggauss(256)
            hh = (np.pi / 2) / K
            a = np.arange(K)[:, None] * hh
            _GLA[n] = ((a + (x[None, :] + 1) * hh / 2).ravel(), np.tile(w * hh / 2, K))
    return _GLA[n]

def nth(b):
    n = int(min(12000, max(96, 1.2 * b + 96)))
    return 128 * ((n + 127) // 128)

def Einvsqrt(beta, s):
    a = beta * (1 + s) / 2
    c = beta * (1 - s) / 2
    def r(x):
        x = np.asarray(x, float)
        out = np.ones_like(x)
        nz = x > 1e-8
        xx = x[nz]
        out[nz] = np.exp(-0.5 * (xx + np.log1p(-np.exp(-2 * xx)) - np.log(2 * xx)))
        out[~nz] = 1 - x[~nz] ** 2 / 12
        return out
    return r(a) * r(c)

class ReggeAction:
    """Leading Regge AdS action of the uplifted functional (ADS4_FUNCTIONAL.md section 0 item 2, general weights):
       FR = 2 int psi2 P_{-1/2+i Lam p}(cosh b) + sigma u int psi3 P + 2a Lam^2 u^2 W(b) + a sigma u^3 P0(b)."""
    def __init__(self, F):
        self.F = F
        self.c2 = CosT(F.psi2); self.c3 = CosT(F.psi3)
    def transforms(self, Lam, beta):
        b = Lam * beta
        th, wt = glang(nth(b))
        s = np.sin(th)
        E = Einvsqrt(beta, s)
        I2 = (2 / np.pi) * np.sum(wt * self.c2(b * s) * E)
        I3 = (2 / np.pi) * np.sum(wt * self.c3(b * s) * E)
        W = (2 / np.pi) * np.sum(wt * (beta * s) ** 2 / 2 * E)
        P0 = (2 / np.pi) * np.sum(wt * E)
        return I2, I3, W, P0
