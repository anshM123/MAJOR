r"""
lsd_mellin (LIGHTSIDE_DIST.md): Mellin representation of EIN's Witten-diagram sigma-pieces ON the swapped-form contour C x C, and the
three ways of applying the swapped-form functional Psi_k to them:

  Omega_pos  : t outside (the whole function, EIN's E[piece x class](a,b));
  Omega_line : Mellin outside, straight Mellin lines (Fubini partner of Omega_pos);
  Omega_OPE  : term-wise in the s-channel (12)(34) OPE of the sigma part = close the delta12 contour to the LEFT and sum the residues
               (double-trace double poles delta12 = -n, plus the single-trace pole of the s-exchange).  This is CLS-D's
               omega|_light (2311.04271 eq. omega_on_light), the object for which the sum rule is a sum over states.

Symanzik/Mellin (ein_dbar docstring):  I_Delta = (1/(2 pi i)^2) int [d delta] prod Gamma(delta_ij) X_ij^{-delta_ij}, Delta_i = 2:
  delta34 = delta12, delta24 = delta13, delta23 = delta14 = 2 - delta12 - delta13,
  prod X^{-delta} = (X14 X23)^{-2} (U/V)^{-delta12} V^{delta13},   U/V = X12 X34/(X14 X23),  V = X14 X23/(X13 X24).
Mellin amplitudes in this common measure Gamma(d12)^2 Gamma(d13)^2 Gamma(d14)^2 (exact, from ein_dbar's operator definitions):
  c0: 1;  c2: s2 = d12^2 + d13^2 + d14^2;  c3: d12 d13 d14;  c4: s2^2;
  xs: Q_s/(d12 - 1), Q_s = (d13 - d14)^2/4 + 1/12;  xt: Q_t/(d14 - 1), Q_t = (d13 - d12)^2/4 + 1/12;  xu: Q_u/(d13 - 1), Q_u = (d12-d14)^2/4+1/12.
On C every X_ij is evaluated with its PRINCIPAL branch (2t - 1 = -e^{i theta}, t = (1 - e^{i theta})/2, arg t in (-pi/2, pi/2),
Re(2t + 1) > 0); the claim that this is EIN's physical-route continuation is checked pointwise (section 'val').
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, digamma

PIECES_M = ('c0', 'c2', 'c3', 'c4', 'xs', 'xt', 'xu', 'xu4', 'r4', 'r4u')


# ------------------------------------------------------------------ Mellin amplitudes and their delta12-derivative (d14 = 2 - d12 - d13)
def M_amp(p, d12, d13):
    d14 = 2 - d12 - d13
    if p == 'c0':
        return np.ones_like(d12 + d13)
    s2 = d12**2 + d13**2 + d14**2
    if p == 'c2':
        return s2
    if p == 'c3':
        return d12*d13*d14
    if p == 'c4':
        return s2**2
    if p == 'xs':
        return ((d13 - d14)**2/4 + 1/12)/(d12 - 1)
    if p == 'xt':
        return ((d13 - d12)**2/4 + 1/12)/(d14 - 1)
    if p == 'xu':
        return ((d12 - d14)**2/4 + 1/12)/(d13 - 1)
    if p == 'r4':                   # top Regge term of c4: s2^2 = (d12 - d14)^4/4 + g (d12 - d14)^2 + g^2
        return (d12 - d14)**4/4
    if p == 'r4u':                  # top Regge term of xu4: Q_u^2 = (d12 - d14)^4/16 + ...
        return (d12 - d14)**4/(16*(d13 - 1))
    if p == 'xu4':                  # u-channel exchange with a degree-4 (polynomial, hence microcausal) residue: Regge growth 4
        return ((d12 - d14)**2/4 + 1/12)**2/(d13 - 1)
    raise KeyError(p)


def dM_amp(p, d12, d13, h=1e-4):
    """total d/d delta12 at fixed delta13 (complex-step-free central difference on a polynomial/rational: exact enough at h = 1e-4?
    No: use exact formulas)."""
    d14 = 2 - d12 - d13
    if p == 'c0':
        return np.zeros_like(d12 + d13)
    s2 = d12**2 + d13**2 + d14**2
    ds2 = 2*d12 - 2*d14
    if p == 'c2':
        return ds2
    if p == 'c3':
        return d13*d14 - d12*d13
    if p == 'c4':
        return 2*s2*ds2
    if p == 'xs':
        Q = (d13 - d14)**2/4 + 1/12; dQ = (d13 - d14)/2
        return dQ/(d12 - 1) - Q/(d12 - 1)**2
    if p == 'xt':
        Q = (d13 - d12)**2/4 + 1/12; dQ = -(d13 - d12)/2
        return dQ/(d14 - 1) + Q/(d14 - 1)**2
    if p == 'xu':
        dQ = (d12 - d14)            # d/dd12 of (d12 - d14)^2/4 with d(d12 - d14) = 2
        return dQ/(d13 - 1)
    if p == 'r4':
        return 2*(d12 - d14)**3
    if p == 'r4u':
        return (d12 - d14)**3/(2*(d13 - 1))
    if p == 'xu4':
        Q = (d12 - d14)**2/4 + 1/12
        return 2*Q*(d12 - d14)/(d13 - 1)
    raise KeyError(p)


# ------------------------------------------------------------------ logs of the cross-ratio combinations on C (principal per factor)
def logs_on_C(t1, t3, y):
    """LA = log(U/V) = [log X12 - log X14] + [log X34 - log X23];  LB = log V = log X14 + log X23 - log X13 - log X24;
    PRE = (X14 X23)^-2.  Each log principal; separable: LA = A(t1) + A(t3), LB = B(t1) + B(t3) - log(8 y)."""
    A = lambda t: np.log(2*t - 1 + 0j) - np.log(2*t + 1 + 0j)
    B = lambda t: np.log(2*t + 1 + 0j) - np.log(t + 0j)
    LA = A(t1) + A(t3)
    LB = B(t1) + B(t3) - np.log(8*y)
    PRE = ((2*t1 + 1)*(2*t3 + 1))**-2.0
    return LA, LB, PRE


def Acol(t):
    return np.log(2*t - 1 + 0j) - np.log(2*t + 1 + 0j)


def Bcol(t):
    return np.log(2*t + 1 + 0j) - np.log(t + 0j)


# ------------------------------------------------------------------ MB value of a piece at one point (validation)
def mb_point(p, X, c1=2/3, c2=2/3, L=14.0, h=0.05):
    """(1/(2 pi i)^2) int int Gamma^2 Gamma^2 Gamma^2 M prod X^{-delta}, straight lines Re d12 = c1, Re d13 = c2, trapezoid.
    X: dict pair -> complex (principal logs)."""
    lg = {k: np.log(complex(v)) for k, v in X.items()}
    la = lg[(1, 2)] + lg[(3, 4)] - lg[(1, 4)] - lg[(2, 3)]
    lb = lg[(1, 4)] + lg[(2, 3)] - lg[(1, 3)] - lg[(2, 4)]
    pre = np.exp(-2*(lg[(1, 4)] + lg[(2, 3)]))
    y = np.arange(-L, L + h/2, h)
    d12 = c1 + 1j*y[:, None]; d13 = c2 + 1j*y[None, :]; d14 = 2 - d12 - d13
    lmeas = 2*(loggamma(d12) + loggamma(d13) + loggamma(d14))
    F = np.exp(lmeas - d12*la + d13*lb)*M_amp(p, d12, d13)
    return pre*np.sum(F)*h*h/(4*np.pi**2)


# ------------------------------------------------------------------ contour on C with geometric grading towards t = 0 (both ends)
def contour_C(s0=1e-6, m=12, q=3.0):
    """t = (1 - e^{i theta})/2, theta in (0, 2 pi); composite Gauss panels [0, s0] (quadratic), [s0, q s0], ..., up to pi, mirrored.
    Returns t, w (dt included), theta."""
    x, wx = np.polynomial.legendre.leggauss(m)
    u = (x + 1)/2; wu = wx/2
    th = [s0*u**2]; wt = [wu*2*s0*u]
    a = s0
    while a < np.pi:
        b = min(a*q, np.pi)
        th.append(a + (b - a)*(x + 1)/2); wt.append(wx*(b - a)/2)
        a = b
    th = np.concatenate(th); wt = np.concatenate(wt)
    th = np.concatenate([th, 2*np.pi - th[::-1]]); wt = np.concatenate([wt, wt[::-1]])
    t = 0.5 - 0.5*np.exp(1j*th)
    dt = -0.5j*np.exp(1j*th)
    return t, wt*dt, th
