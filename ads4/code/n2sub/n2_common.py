"""
n2sub/n2_common.py (N2_SUBREGGE.md): shared helpers for the N2 trace (which weight multiplies which data) and Task 2.

Everything here is re-implemented from the published definitions; project modules are imported READ-ONLY (sys.path) only where noted.
  * X, Y, psi2 on the 0.04 grid: ads4/prl/data/d_abel_X_<Dg>.txt  (functional = Re sum_k dnu X_k Psi^_(EE,nu_k)/(h_- + h_+)(nu_k))
  * R_-+ (EE, k = 2): ads4/prl/data/kernels.json  'D_k2_EE'
  * the local de-shift operator of the chains (q3_analysis._shift_generic, w2 EE mode, fw_gamma):
        A[F](nu) = (R_- F(nu + i) + R_+ F(nu - i)) / (R_- + R_+)       (gp = (1+u)/2 = R_-/(R_-+R_+) on F(x + i/m))
  * P_{-1/2+i nu}(cosh b) for COMPLEX nu by Mehler-Dirichlet in the Lemma-C substitution (ls_core.OmLR form):
        P = (2/pi) int_0^{pi/2} cos(nu s)/cosh(s/2) dphi,  s = 2 asinh(sinh(b/2) sin phi)
  * D_h, env_h (l9_analysis), the kappa image (sm_curvfix.curv_image, SM_CURVIMG = kappa), all valid for complex x.
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
from scipy.special import jv

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')


def lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


_XD = {}


def xdata(Dg=1000):
    """nu grid (0.04 / 0.02), complex X, complex Y, psi2 (p-measure normalisation; functional = Re sum dnu X Psi^/Lam ...)."""
    if Dg not in _XD:
        d = np.loadtxt(os.path.join(DATA, 'd_abel_X_%d.txt' % Dg), comments='#')
        _XD[Dg] = dict(nu=d[:, 0], X=d[:, 1] + 1j*d[:, 2], Y=d[:, 3] + 1j*d[:, 4], psi2=d[:, 5], dnu=d[1, 0] - d[0, 0])
    return _XD[Dg]


_RK = {}


def Rk():
    if not _RK:
        import sympy as sp
        k = json.load(open(os.path.join(DATA, 'kernels.json')))['D_k2_EE']
        nu = sp.Symbol('nu')
        _RK['m'] = sp.lambdify(nu, sp.sympify(k['R_minus']), 'numpy')
        _RK['p'] = sp.lambdify(nu, sp.sympify(k['R_plus']), 'numpy')
    return _RK['m'], _RK['p']


def gpm(nu):
    """(g_+, g_-) = (R_-, R_+)/(R_- + R_+): weights of F(nu + i), F(nu - i) in the raw forward-normalised image."""
    Rm, Rp = Rk()
    nu = np.asarray(nu, complex)
    a, b = Rm(nu), Rp(nu)
    return a/(a + b), b/(a + b)


def Aop(F, nu):
    """local de-shift operator applied to an analytic function F(nu) (complex nu ok)."""
    gp, gm = gpm(nu)
    nu = np.asarray(nu, complex)
    return gp*F(nu + 1j) + gm*F(nu - 1j)


_GLc = {}


def gl01(n):
    if n not in _GLc:
        x, w = np.polynomial.legendre.leggauss(n)
        _GLc[n] = ((x + 1)/2, w/2)
    return _GLc[n]


def Pnu(beta, nu, nphi=600):
    """P_{-1/2 + i nu}(cosh beta), nu real or complex array (Mehler-Dirichlet, smooth integrand)."""
    ph, w = gl01(nphi); ph = ph*np.pi/2; w = w*np.pi/2
    s = 2*np.arcsinh(np.sinh(beta/2)*np.sin(ph))
    wt = (2/np.pi)*w/np.cosh(s/2)
    nu = np.atleast_1d(np.asarray(nu, complex))
    out = np.empty(nu.shape, complex)
    fl = nu.ravel(); o = out.ravel()
    for i0 in range(0, len(fl), 2048):
        o[i0:i0 + 2048] = (np.cos(np.outer(fl[i0:i0 + 2048], s))*wt).sum(1)
    return out


def jac(n, a, b, x):
    x = np.asarray(x, complex)
    p0 = np.ones_like(x); p1 = 0.5*(a - b + (a + b + 2)*x)
    if n == 0:
        return p0
    for k in range(2, n + 1):
        c1 = 2*k*(k + a + b)*(2*k + a + b - 2)
        c2 = (2*k + a + b - 1)*((2*k + a + b)*(2*k + a + b - 2)*x + a*a - b*b)
        c3 = 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
        p0, p1 = p1, (c2*p1 - c3*p0)/c1
    return p1


def Dsh(h, J, x):
    c = 1 - 2*np.asarray(x, complex)**2
    return jac(J, 0, 0, c) if h == '++' else ((1 + c)/2)**4*jac(J - 4, 0, 8, c)


def env(h, x, k):
    y = np.asarray(x, complex)**2
    e = (2 - y) if h == '++' else (2 - y)/(1 - y)**4
    return e if k == 2 else e/(2 - y)


def kappa(kind, h, J, m, x):
    """sm_curvfix.curv_image (kappa): leading-Regge curvature content of a row's Gamma, c2 units; complex x allowed."""
    x = np.asarray(x, complex)
    b = 2*np.arcsinh((J + 0.5)/m); A = np.sqrt(b/np.sinh(b))
    te = 2*np.arcsin(x)*(b/2)/np.sinh(b/2); xe = np.sin(te/2)
    osc = A*Dsh(h, J, xe) - Dsh(h, J, x)
    reg = Pnu(b, x*m) - A*jv(0, x*m*b)
    if kind == 'C2':
        return m*m*(osc + 2*reg/env(h, x, 2))
    s = 1.0 if h == '++' else -1.0
    return m*m*s*(osc + reg/env(h, x, 3))


def A_on_x(fx, nu, m):
    """A applied to a function of x = nu/m given as fx(x complex)."""
    return Aop(lambda v: fx(v/m), nu)
