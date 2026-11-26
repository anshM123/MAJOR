"""
fullpos4/fp4_pf.py (FULLPOS4.md §1) [DV]: per-frequency outgoing contour for K_p(nu) = L(F_{p,nu}), double precision (mpmath.fp).

  K_p = (C/2pi) int_{al} int_{ga} W(al, ga) dga dal,   W = e^{i S al - i a ga} xi^{s-1} F(xi^2) / sin(al),
  al = A - i eps, A in [-pi/2, pi/2]  (eps = 1/S);  base gamma path: ga = al - i v, v in [0, oo)   (= ANGULAR_ROUTE §59)
  xi = sin((ga-al)/2)/sin((ga+al)/2),  beta = log tan(ga/2) - log tan(al/2)  (= 2 atanh xi),  F = (1-xi^2)^{-p} P(beta)
  Outgoing split: P = P+ + P-,  K_p = 2 Re I+,  I+ = same integral with P+.
  Deformed path for I+ (Im beta >= 0 except on the damped tail):
     A > 0:  ga: al -> X - i eps (horizontal, X = A + 3 eps), X - i eps -> X (up), X -> pi/2 (real), pi/2 -> pi/2 - i oo
             (tail, damped by e^{-(a - nu) t}); Re beta >= 0 and Im beta >= 0 on the first three pieces
     A < 0:  base ray ga = al - i v (Im beta > 0 on it)
"""
import sys, math, cmath
sys.dont_write_bytecode = True
import numpy as np
import mpmath as mp

fp = mp.fp
GLX, GLW = np.polynomial.legendre.leggauss(24)


def gl_panels(edges, nn=None):
    x, w = (GLX, GLW) if nn is None else np.polynomial.legendre.leggauss(nn)
    X, Wt = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        X.append(a + (b - a)*(x + 1)/2); Wt.append((b - a)/2*w)
    return np.concatenate(X), np.concatenate(Wt)


def graded(a, b, h0, r=1.6):
    """edges from a to b, first panel h0, growing by r (b > a)."""
    e = [a]; h = h0
    while e[-1] + h < b:
        e.append(e[-1] + h); h *= r
    e.append(b)
    return e


from scipy.special import loggamma as _lg, digamma as _dg


def _h2f1_zb(nu, w):
    """2F1(1/2, 1/2 - i nu; 1 - i nu; w) (zero-balanced): direct series for |w| <= 0.6, the log expansion about w = 1 for
    |1 - w| <= 0.75, mpmath (30 digits) otherwise."""
    a = 0.5; b = 0.5 - 1j*nu; c = 1 - 1j*nu
    if abs(w) <= 0.6:
        t = 1.0 + 0j; s = t; k = 0
        while abs(t) > 1e-18*abs(s) or k < 3:
            t *= (a + k)*(b + k)/((c + k)*(k + 1))*w; s += t; k += 1
        return s
    u = 1 - w
    if abs(u) <= 0.75:
        L = cmath.log(u)
        pre = cmath.exp(_lg(c) - _lg(a) - _lg(b))
        t = 1.0 + 0j; s = 0j; k = 0
        da, db, d1 = _dg(a), _dg(b), _dg(1.0)
        while True:
            term = t*(2*d1 - da - db - L)
            s += term
            if k > 3 and abs(t)*(abs(L) + 40) < 1e-18*abs(s):
                break
            t *= (a + k)*(b + k)/((k + 1)**2)*u
            da += 1/(a + k); db += 1/(b + k); d1 += 1/(k + 1); k += 1
        return pre*s
    mp.mp.dps = 30
    return complex(mp.hyp2f1(mp.mpf(1)/2, mp.mpf(1)/2 - 1j*mp.mpf(nu), 1 - 1j*mp.mpf(nu), mp.mpc(w)))


def Pplus(nu, beta):
    pre = cmath.exp(_lg(1j*nu) - _lg(0.5 + 1j*nu))/math.sqrt(math.pi)
    return pre*cmath.exp((-0.5 + 1j*nu)*beta)*_h2f1_zb(nu, cmath.exp(-2*beta))


def Pfull(nu, beta):
    return fp.legenp(-0.5 + 1j*nu, 0, cmath.cosh(beta), type=3)


def lt(t):
    return cmath.log(cmath.tan(t/2))


def Cconst(n, N, s):
    return (N - n)*math.exp(math.lgamma(n + 1) + math.lgamma(N + s) - math.lgamma(n + s) - math.lgamma(N + 1))


def W(al, ga, p, nu, s, S, a, which):
    beta = lt(ga) - lt(al)
    xi = cmath.tanh(beta/2)
    P = Pplus(nu, beta) if which == '+' else Pfull(nu, beta)
    return cmath.exp(1j*S*al - 1j*a*ga)*xi**(s - 1)*(1 - xi*xi)**(-p)*P/cmath.sin(al)


def gamma_path(al, A, mode, a, eps):
    """list of (gamma nodes, d gamma weights) for the path of alpha."""
    nodes = []; wts = []
    if mode == 'base' or A < 0:
        # ga = al - i v, v in [0, oo): graded near 0, tail to 60/a
        e = graded(0.0, 60.0/a + 5*eps, eps*1e-9)
        v, wv = gl_panels(e)
        return al - 1j*v, -1j*wv
    # A > 0 deformed: horizontal al -> X - i eps (X = A + 3 eps), vertical up to X, real X -> pi/2, tail pi/2 - i t
    X = min(A + 3*eps, math.pi/2)
    xh, wxh = gl_panels(graded(A, X, eps*1e-9))
    g0 = xh - 1j*eps; w0 = wxh + 0j
    y, wy = gl_panels(graded(0.0, eps, eps*1e-2)[::-1] if False else graded(0.0, eps, eps*1e-2))
    g1 = X - 1j*(eps - y); w1 = 1j*wy
    if X < math.pi/2:
        x, wx = gl_panels(graded(X, math.pi/2, eps*1e-2))
    else:
        x, wx = np.array([]), np.array([])
    g2 = x + 0j; w2 = wx + 0j
    t, wt = gl_panels(graded(0.0, 60.0/a + 1.0, 1e-3))
    g3 = math.pi/2 - 1j*t; w3 = -1j*wt
    g1 = np.concatenate([g0, g1])
    w1 = np.concatenate([w0, w1])
    return np.concatenate([g1, g2, g3]), np.concatenate([w1, w2, w3])


def contour(p, nu, J, b, Delta, mode='def', which='+', nal=None):
    n = J - b; N = Delta - b - 2; s = 2*b + 1; a = N - n; B = 2*J + 1; S = B + a
    eps = 1.0/S
    C = Cconst(n, N, s)
    eA = graded(0.0, math.pi/2, eps*1e-3, 1.4)
    Ap, wAp = gl_panels(eA)
    As = np.concatenate([-Ap[::-1], Ap]); wA = np.concatenate([wAp[::-1], wAp])
    tot = 0j; tabs = 0.0
    for A, wa in zip(As, wA):
        al = A - 1j*eps
        g, wg = gamma_path(al, A, mode, a, eps)
        inner = 0j; iabs = 0.0
        for gg, ww in zip(g, wg):
            val = W(al, gg, p, nu, s, S, a, which)
            inner += val*ww; iabs += abs(val*ww)
        tot += inner*wa; tabs += iabs*wa
    I = C/(2*math.pi)*tot; Iabs = C/(2*math.pi)*tabs
    return I, Iabs


if __name__ == '__main__':
    sys.path.insert(0, '.')
    import fp4_kern as K
    mp.mp.dps = 30
    p, nu, J, b, Delta = (int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]))
    mode = sys.argv[6] if len(sys.argv) > 6 else 'def'; which = sys.argv[7] if len(sys.argv) > 7 else '+'
    ex = K.exact(p, nu, J - b, mp.mpf(Delta) - b - 2, 2*b + 1)
    I, Ia = contour(p, nu, J, b, Delta, mode, which)
    val = 2*I.real if which == '+' else I
    print(f'p={p} nu={nu} J={J} b={b} Delta={Delta} mode={mode} {which}: exact={mp.nstr(ex, 14)} contour={val} '
          f'(I={I}) abs={Ia:.3e} rel.err={abs(complex(val) - float(ex))/abs(float(ex)):.2e}', flush=True)
