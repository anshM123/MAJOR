#!/usr/bin/env python3
"""
Squashed-sphere C_T and round-sphere F for 3d N=2 U(N) SQCD-type matrix models.

Conventions (Hama-Hosomichi-Lee; checked against free chiral, XYZ, N=4 SQED_n and N=4 SU(2) SQCD_k
closed forms of 1910.03600):
  chiral of R-charge D and weight rho:  s_b(iQ(1-D)/2 - rho(sigma))
  log s_b(x) = -i int_0^inf dw/w [ sin(2xw)/(2 sinh(bw) sinh(w/b)) - x/w ]
  vector: prod_{alpha>0} 4 sinh(pi b alpha.sigma) sinh(pi alpha.sigma/b),   Weyl 1/|W|
  F = -log|Z|,  C_T/C_T(real scalar) = (32/pi^2) d^2F/db^2 |_{b=1}.
Everything is invariant under b -> 1/b, so d/db = 0 at b=1 and d^2F/db^2 = -<d^2 log(integrand)/db^2>.

U(N) with N_f fundamental pairs at b: by Andreief,
  Z_b = det_{j,k=0}^{N-1} M_jk(b),  M_jk = int dsigma g_b(sigma) exp[2 pi (b j + k/b - (N-1)(b+1/b)/2) sigma],
  g_b = [s_b(iy - sigma) s_b(iy + sigma)]^{N_f},  y = Q(1-D)/2.
"""
import numpy as np
import mpmath as mp
from scipy.integrate import quad
from scipy.special import polygamma

mp.mp.dps = 25
PI = np.pi


# ------------------------------------------------------------ one-loop pieces
def ell(z):
    z = mp.mpc(z)
    e = mp.exp(2j * mp.pi * z)
    return -z * mp.log(1 - e) + 0.5j * (mp.pi * z ** 2 + mp.polylog(2, e) / mp.pi) - 1j * mp.pi / 12


def logpair_b1(s, D):
    """log[s_1(i(1-D)-s) s_1(i(1-D)+s)] = l(1-D+is)+l(1-D-is)  (real)."""
    v = ell(1 - D + 1j * s) + ell(1 - D - 1j * s)
    return float(mp.re(v))


def logpair_int(s, D, b):
    """Same at general b from the integral representation (checks)."""
    Q = b + 1 / b
    y = Q * (1 - D) / 2

    def f(w):
        if w < 1e-6:
            return 0.0
        r = np.exp((2 * y - Q) * w) * (-np.expm1(-4 * y * w)) / (np.expm1(-2 * b * w) * np.expm1(-2 * w / b))
        return (2 * r * np.cos(2 * s * w) - 2 * y / w) / w
    # [0,1] direct, [1,inf) split into oscillatory and non-oscillatory parts
    v1, _ = quad(f, 0, 1, limit=400, epsabs=1e-14)
    def g(w):
        return 2 * np.exp((2 * y - Q) * w) * (-np.expm1(-4 * y * w)) / (np.expm1(-2 * b * w) * np.expm1(-2 * w / b)) / w
    if s == 0:
        v2, _ = quad(g, 1, np.inf, limit=400, epsabs=1e-14)
    else:
        v2, _ = quad(g, 1, np.inf, weight="cos", wvar=2 * s, limlst=200)
    return v1 + v2 - 2 * y


def logchiral_b(D, b):
    """log s_b(iQ(1-D)/2) for a neutral chiral = half the pair at sigma=0."""
    return 0.5 * logpair_int(0.0, D, b)


def P2(s, D):
    """d^2/db^2 at b=1 of log pair(s): int_0^inf dw/w [cos(2sw) K(w) - 2a/w], a=1-D."""
    a = 1 - D

    def K(w):  # [2aw cosh(2aw) - sinh(2aw)(2w coth w - 2w^2 csch^2 w)] / sinh^2 w, stable
        em = np.exp(-2 * w)
        coth = (1 + em) / (1 - em)
        csch2 = 4 * em / (1 - em) ** 2
        # cosh(2aw)/sinh^2 w and sinh(2aw)/sinh^2 w via exponentials
        c_over = 2 * (np.exp((2 * a - 2) * w) + np.exp((-2 * a - 2) * w)) / (1 - em) ** 2
        s_over = 2 * (np.exp((2 * a - 2) * w) - np.exp((-2 * a - 2) * w)) / (1 - em) ** 2
        return 2 * a * w * c_over - s_over * (2 * w * coth - 2 * w * w * csch2)

    def f(w):
        if w < 2e-3:
            return -4 * a * s * s + 4 * a ** 3 - 10 * a / 3
        return (np.cos(2 * s * w) * K(w) - 2 * a / w) / w
    v1, _ = quad(f, 0, 1, limit=400, epsabs=1e-13)
    if s == 0:
        v2, _ = quad(lambda w: K(w) / w, 1, np.inf, limit=400, epsabs=1e-14)
    else:
        v2, _ = quad(lambda w: K(w) / w, 1, np.inf, weight="cos", wvar=2 * s, limlst=200)
    return v1 + v2 - 2 * a


def V2(x):
    """d^2/db^2 at b=1 of log[sinh(pi b x) sinh(pi x/b)]."""
    if abs(x) < 1e-8:
        return 0.0
    t = PI * x
    return 2 * t / np.tanh(t) - 2 * t * t / np.sinh(t) ** 2


# ------------------------------------------------------------ grids
class Grid:
    def __init__(self, D, L=9.0, h=0.01, need_p2=True):
        self.D = D
        s = np.arange(0, L + h / 2, h)
        lp = np.array([logpair_b1(x, D) for x in s])
        p2 = np.array([P2(x, D) for x in s]) if need_p2 else np.zeros_like(s)
        self.s = np.concatenate([-s[:0:-1], s])
        self.lp = np.concatenate([lp[:0:-1], lp])
        self.p2 = np.concatenate([p2[:0:-1], p2])
        self.h = h


def UN_logZ_and_d2(grid, N, Nf, eta=0.0):
    """Returns log Z(b=1) (up to the Weyl factor, which is exactly cancelled by Andreief's N!) and
    d log Z/db, d^2 log Z/db^2 at b=1, via the Andreief determinant (trapezoid rule, exponentially accurate)."""
    s, h = grid.s, grid.h
    lg = Nf * grid.lp
    lg2 = Nf * grid.p2
    shift = lg.max()
    g = np.exp(lg - shift)
    M = np.zeros((N, N)); M1 = np.zeros((N, N)); M2 = np.zeros((N, N))
    for j in range(N):
        for k in range(N):
            # eta: R-mixing with U(1)_T, factor exp(pi Q eta sum sigma) (R[V_pm] -> R0 -/+ eta)
            E = np.exp(2 * PI * (j + k - (N - 1) + eta) * s)
            c1 = j - k
            c2 = 2 * k - (N - 1) + eta
            M[j, k] = h * np.sum(g * E)
            M1[j, k] = h * np.sum(g * E * 2 * PI * c1 * s)
            M2[j, k] = h * np.sum(g * E * (lg2 + (2 * PI * c1 * s) ** 2 + 2 * PI * c2 * s))
    Mi = np.linalg.inv(M)
    sign, ld = np.linalg.slogdet(M)
    logZ = ld + N * shift
    d1 = np.trace(Mi @ M1)
    d2 = np.trace(Mi @ M2) - np.trace(Mi @ M1 @ Mi @ M1)
    return logZ, d1, d2, sign


def U2_direct(grid, Nf):
    """Independent check for U(2): 2D trapezoid of <V2(s1-s2) + Nf(P2(s1)+P2(s2))>."""
    s, h = grid.s, grid.h
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    x = S1 - S2
    vdm = 4 * np.sinh(PI * x) ** 2
    lw = Nf * (grid.lp[:, None] + grid.lp[None, :])
    shift = lw.max()
    w = vdm * np.exp(lw - shift)
    Z = 0.5 * h * h * w.sum()
    t = PI * x
    with np.errstate(invalid="ignore", divide="ignore"):
        v2 = np.where(np.abs(t) < 1e-8, 0.0, 2 * t / np.tanh(t) - 2 * t * t / np.sinh(t) ** 2)
    obs = v2 + Nf * (grid.p2[:, None] + grid.p2[None, :])
    d2 = (0.5 * h * h * (w * obs).sum()) / Z
    return np.log(Z) + shift, d2


def CT_from_d2logZ(d2):
    return -(32 / PI ** 2) * d2


def chiral_d2(D):
    """d^2/db^2 log s_b(iQ(1-D)/2) at b=1 (neutral chiral) = P2(0,D)/2."""
    return 0.5 * P2(0.0, D)


def chiral_logZ_b1(D):
    if abs(1 - D) < 1e-14:
        return 0.0          # l(0) = 0
    return float(mp.re(ell(1 - D)))
