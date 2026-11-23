"""p4win/p4w_common.py (P4_WINDOW.md): shared helpers.

* Pd(beta, nu): P_{-1/2+i nu}(cosh beta) and dP/dbeta by Mehler-Dirichlet with the exact beta-derivative of the integrand
  (Gauss-Legendre, node count adaptive in nu*beta; smooth integrand).
* Series: CHMRS scalar 'series Pi' orders 1..4 (p4w_series cache), T_k(nu; beta) = r0_k P + r1_k P', flat parts g_k(b).
* NF02 weights psi2, psi3 (ls_core.NFw, read-only), the published EE grid (lsw_common, read-only: X, psi2 on the 0.04 grid, act()).
* The monopole constant D (PW_IMAGES (4.1); p2con/pc_run.D1000, read-only import).
Nothing outside code/p4win is written.
"""
import os, sys, pickle, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import sympy as sp
from scipy.special import j0, j1

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
CACHE = os.path.join(HERE, 'cache')
for sub in ('lsw', 'largespin', 'p2con', 'spinfunc', 'exactskel'):
    p = os.path.join(CODE, sub)
    if p not in sys.path:
        sys.path.append(p)
import lsw_common as W          # noqa: E402  (read-only)
import ls_core as LC            # noqa: E402  (read-only)

L = W.L
D1000 = -9.983321e-4            # PW_IMAGES (4.1), INDEP_CHECK2 1b: D(1000) = -9.9833211e-4
NF = LC.NFw(1000)


def psi2(p):
    return NF.psi2(np.asarray(p, float))


def psi3(p):
    return NF.psi3(np.asarray(p, float))


_GL = {}


def gl(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = ((x + 1)/2*np.pi/2, w/2*np.pi/2)
    return _GL[n]


def Pd(beta, nu, n=None):
    """(P, dP/dbeta) of P_{-1/2+i nu}(cosh beta) for real or complex nu array.
    P = (2/pi) int_0^{pi/2} cos(nu s)/cosh(s/2) dphi, sinh(s/2) = sinh(beta/2) sin(phi)."""
    nu = np.atleast_1d(np.asarray(nu, complex))
    numax = float(np.max(np.abs(nu))) if nu.size else 0.0
    n = n or int(max(200, 1.5*numax*beta + 120))
    ph, w = gl(n)
    sh = np.sinh(beta/2)*np.sin(ph)
    s = 2*np.arcsinh(sh); cs = np.sqrt(1 + sh*sh)          # cosh(s/2)
    ds = np.cosh(beta/2)*np.sin(ph)/cs                       # ds/dbeta
    wt = (2/np.pi)*w/cs
    dwt = -(2/np.pi)*w*np.tanh(s/2)/(2*cs)*ds              # d/dbeta of 1/cosh(s/2)
    P = np.empty(nu.shape, complex); dP = np.empty(nu.shape, complex)
    fl = nu.ravel(); Po = P.ravel(); dPo = dP.ravel()
    for i0 in range(0, len(fl), 1024):
        a = np.outer(fl[i0:i0 + 1024], s)
        c = np.cos(a); sn = np.sin(a)
        Po[i0:i0 + 1024] = (c*wt).sum(1)
        dPo[i0:i0 + 1024] = (c*dwt).sum(1) - ((sn*fl[i0:i0 + 1024, None])*(wt*ds)).sum(1)
    return P, dP


# ----------------------------------------------------------------------------------------------- CHMRS series
_SER = {}
nu2s, etas, bs, J0s, J1s = sp.symbols('nu2 eta b J0 J1')


def series():
    if not _SER:
        d = pickle.load(open(os.path.join(CACHE, 'series.pkl'), 'rb'))
        for k, v in d.items():
            r0 = sp.sympify(v['r0'], locals=sp.__dict__)
            r1 = sp.sympify(v['r1'], locals=sp.__dict__); g = sp.sympify(v['g'], locals=sp.__dict__)
            _SER[k] = dict(r0=sp.lambdify((etas, nu2s), r0, 'numpy'), r1=sp.lambdify((etas, nu2s), r1, 'numpy'),
                           g=sp.lambdify((bs, J0s, J1s), g, 'numpy'), r0s=r0, r1s=r1, gs=g)
    return _SER


def T(k, beta, nu, P=None, dPb=None):
    """k-th order CHMRS term (coefficient of 1/m^{2k}) per unit leading image, at real nu (array); also |r0 P| + |r1 P'| (condition)."""
    S = series()[k]
    if P is None:
        P, dPb = Pd(beta, nu)
    P = np.real(P); dPe = np.real(dPb)/np.sinh(beta)
    eta = np.cosh(beta); n2 = np.asarray(nu, float)**2
    a = S['r0'](eta, n2)*P; b = S['r1'](eta, n2)*dPe
    return a + b, np.abs(a) + np.abs(b)


def gflat(k, bb):
    """flat-scaling part g_k(b) of order k (coefficient of y^k = nu^{2k}/m^{2k}), at b = bhat p (array, b > 0)."""
    bb = np.asarray(bb, float)
    return series()[k]['g'](bb, j0(bb), j1(bb))


# ----------------------------------------------------------------------------------------------- states, floors, budgets
def theoremW_states():
    """Theorem W's 47 sample states with H, d2 (computed O(1/m^2) cost, block + kernel), c4 (Theorem W's K4 unit cost, per 1/m^4)."""
    d = json.load(open(os.path.join(CODE, 'p2con', 'cache', 'jstar_full.json')))
    return d['states']


def act_grid(f):
    """k = 2 de-shift action on the published 0.04 grid (n2 units, ring weight 2): 2 sum psi2 f dnu / Lam."""
    return 2*(W.PS*f).sum()*W.DN/L


def beta_of(J, mfl):
    return 2*np.arcsinh((J + 0.5)/(mfl*L))
