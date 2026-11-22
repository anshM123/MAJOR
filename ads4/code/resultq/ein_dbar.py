r"""
ein_dbar (RESULT_Q_EINSTEIN_CHECK.md): Witten-diagram (D-function) correlators in the CLS-D gauge of v4c_core / og_free, analytically
continued to the swapped-form contour C x C, for the Lemma-A test.

Objects.
  I_Delta(X) := int_0^inf prod_i dt_i t_i^{Delta_i - 1} exp(-sum_{i<j} t_i t_j X_ij)       (Symanzik; = D-function / its constant)
             =  (1/(2 pi i)^2) int [d delta] prod_{i<j} Gamma(delta_ij) X_ij^{-delta_ij}  (Mellin, sum_j delta_ij = Delta_i).
  I_{1111} = Phi(z, zb) / (2 X13 X24),   Phi = [2 Li2(z) - 2 Li2(zb) + log(z zb) log((1-z)/(1-zb))]/(z - zb),
  U = X12 X34/(X13 X24) = z zb,  V = X14 X23/(X13 X24) = (1-z)(1-zb).                                   (validated in ein_run (1))
  Raising:  I_{Delta + e_i + e_j} = -d/dX_ij I_Delta  (exact).   Euler operators theta_ij = -X_ij d/dX_ij multiply the Mellin
  integrand by delta_ij, so a Mellin amplitude P(delta) (polynomial) is P(theta) I_Delta, and a pole 1/(delta_12 - 1) is
  X12^{-1} I_{Delta - e1 - e2}.  Every object is therefore a finite sum  c * prod X_ij^{e_ij} * Phi^{(m,n)},
  Phi^{(m,n)} := (d/dlogU)^m (d/dlogV)^n Phi,  and only Phi needs analytic continuation.

Gauge (og_free / v4c_core), effective squared distances (positive on the totally spacelike real domain t1, t3 > 1/2):
  X12 = 2 t1 - 1, X34 = 2 t3 - 1, X13 = 2 y t1 t3, X14 = 2 t1 + 1, X23 = 2 t3 + 1, X24 = 4,  y = 2 sin^2(dphi/2);
  with a = 1/t1, b = 1/t3:  U = (2-a)(2-b)/(8y),  V = (2+a)(2+b)/(8y)   (= og_free's u', v').
  The contour C (|t - 1/2| = 1/2) is the line Re a = 1.  The Wightman function is analytic for a, b in C \ ((-inf,0] u [2,inf))
  (light-ray analyticity + microcausality); at a = b -> 0+ (t -> +inf) the configuration is Euclidean (zb = z*), where the principal
  closed form is correct.  Phi is continued from there along the log-spiral a(s) = eps (a/eps)^s (then the same in b), tracking
  z, zb (roots of w^2 - (1+U-V) w + U), log z, log(1-z), Li2(z) (and the same for zb) by continuity.
  Phi^{(m,n)} at a node: Cauchy/FFT on a polydisc |dlogU| = |dlogV| = RHO, every polydisc point reached by continuation from the node.
"""
import numpy as np
from scipy.special import spence

TWO_PI_I = 2j*np.pi
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
PIDX = {p: i for i, p in enumerate(PAIRS)}
AU = np.array([1, -1, 0, 0, -1, 1])      # exponents of X_ij in U
BV = np.array([0, -1, 1, 1, -1, 0])      # exponents of X_ij in V


def li2p(z):
    """principal Li2 (cut [1, inf))."""
    return spence(1 - z)


# ------------------------------------------------------------------ branch tracking
class State:
    """continued values: z, zb and their blocks; arrays of a common shape."""
    __slots__ = ('z', 'zb', 'lz', 'lzb', 'l1z', 'l1zb', 'L2z', 'L2zb')

    def copy(self):
        s = State()
        for k in self.__slots__:
            setattr(s, k, getattr(self, k).copy())
        return s


def roots(U, V):
    S = 1 + U - V
    d = np.sqrt(S*S - 4*U + 0j)
    return (S + d)/2, (S - d)/2


def principal_state(U, V):
    z, zb = roots(U, V)
    s = State()
    s.z, s.zb = z, zb
    s.lz, s.lzb = np.log(z), np.log(zb)
    s.l1z, s.l1zb = np.log(1 - z), np.log(1 - zb)
    s.L2z, s.L2zb = li2p(z), li2p(zb)
    return s


def _snap_log(prev, w):
    p = np.log(w)
    return p + TWO_PI_I*np.round((prev - p).imag/(2*np.pi))


def _snap_li2(pred, w, lw):
    """nearest branch Li2p(w) + 2 pi i m Log(w) + 4 pi^2 l (m, l integers) to the prediction."""
    p = li2p(w); Lp = np.log(w)
    D = pred - p
    c = TWO_PI_I*Lp
    with np.errstate(divide='ignore', invalid='ignore'):
        m = np.where(np.abs(c.imag) > 1e-12, np.round(D.imag/np.where(np.abs(c.imag) > 1e-12, c.imag, 1.0)), 0.0)
    l = np.round((D.real - m*c.real)/(4*np.pi**2))
    out = p + m*c + 4*np.pi**2*l
    # fall back to the raw prediction if the lattice fit is poor (|z| ~ 1 degenerate case): try neighbours
    best = out; err = np.abs(out - pred)
    for dm in (-1, 0, 1):
        for dl in (-1, 0, 1):
            cand = p + (m + dm)*c + 4*np.pi**2*(l + dl)
            e = np.abs(cand - pred)
            sel = e < err
            best = np.where(sel, cand, best); err = np.where(sel, e, err)
    return best


def step(st, U, V):
    """advance state to (U, V) (assumed close)."""
    r1, r2 = roots(U, V)
    keep = np.abs(r1 - st.z) + np.abs(r2 - st.zb) <= np.abs(r2 - st.z) + np.abs(r1 - st.zb)
    z = np.where(keep, r1, r2); zb = np.where(keep, r2, r1)
    new = State()
    new.z, new.zb = z, zb
    new.lz = _snap_log(st.lz, z); new.lzb = _snap_log(st.lzb, zb)
    new.l1z = _snap_log(st.l1z, 1 - z); new.l1zb = _snap_log(st.l1zb, 1 - zb)
    # Li2' = -log(1-w)/w ; trapezoid prediction
    pz = st.L2z + 0.5*(-st.l1z/st.z - new.l1z/z)*(z - st.z)
    pzb = st.L2zb + 0.5*(-st.l1zb/st.zb - new.l1zb/zb)*(zb - st.zb)
    new.L2z = _snap_li2(pz, z, new.lz); new.L2zb = _snap_li2(pzb, zb, new.lzb)
    return new


def phi_of(st, kind='sv'):
    """Phi on the tracked branch.  kind='sv': the single-valued Dbar_1111 closed form;  'nsv': control without the log term
    (not single-valued in the Euclidean region, hence not a microcausal correlator)."""
    num = 2*st.L2z - 2*st.L2zb
    if kind == 'sv':
        num = num + (st.lz + st.lzb)*(st.l1z - st.l1zb)
    return num/(st.z - st.zb)


def UV_of(a, b, y):
    return (2 - a)*(2 - b)/(8*y), (2 + a)*(2 + b)/(8*y)


def continue_ab(a, b, y, eps=1e-3, nsteps=2500):
    """continue Phi's state from (a, b) = (eps, eps) to the targets (a, b) (broadcast arrays) along log spirals, first a then b.
    Returns the State at the targets.  Max |z - zb| collision proximity recorded in .minsep (diagnostic)."""
    a = np.asarray(a, complex); b = np.asarray(b, complex)
    a, b = np.broadcast_arrays(a, b)
    U0, V0 = UV_of(eps, eps, y)
    st = principal_state(np.full(a.shape, U0, complex), np.full(a.shape, V0, complex))
    minsep = np.abs(st.z - st.zb)
    la = np.log(a/eps); lb = np.log(b/eps)
    for s in np.linspace(0, 1, nsteps + 1)[1:]:
        U, V = UV_of(eps*np.exp(s*la), eps, y)
        st = step(st, U, V); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    for s in np.linspace(0, 1, nsteps + 1)[1:]:
        U, V = UV_of(a, eps*np.exp(s*lb), y)
        st = step(st, U, V); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    return st, minsep


def continue_path(points, y, eps=1e-3, nsub=400):
    """continue along an explicit polyline of (a, b) points (list of tuples of scalars), starting at (eps, eps).  Path-independence
    check.  Each segment is subdivided geometrically in log (spiral) if both ends have Re > 0, else linearly."""
    st = principal_state(np.array([UV_of(eps, eps, y)[0]]), np.array([UV_of(eps, eps, y)[1]]))
    cur = (eps + 0j, eps + 0j)
    for (pa, pb) in points:
        for s in np.linspace(0, 1, nsub + 1)[1:]:
            aa = cur[0] + s*(pa - cur[0]); bb = cur[1] + s*(pb - cur[1])
            U, V = UV_of(np.array([aa]), np.array([bb]), y)
            st = step(st, U, V)
        cur = (pa, pb)
    return st


def phi_derivs(st, U, V, order, rho=0.5, nc=20, nmini=12, kind='sv'):
    """table D[m][n] = (d/dlogU)^m (d/dlogV)^n Phi at the tracked node state st (arrays of shape S), m + n <= order.
    Cauchy-FFT on the polydisc, each polydisc point reached from the node by a straight path in (logU, logV) of nmini steps."""
    S = st.z.shape
    w = np.exp(2j*np.pi*np.arange(nc)/nc)
    e1 = rho*w[:, None]; e2 = rho*w[None, :]                        # (nc, nc)
    vals = np.empty(S + (nc, nc), complex)
    base = st
    # broadcast state to (S, nc, nc)
    b = State()
    for k in State.__slots__:
        setattr(b, k, np.broadcast_to(getattr(base, k)[..., None, None], S + (nc, nc)).copy())
    lU = np.log(U)[..., None, None]; lV = np.log(V)[..., None, None]
    for s in np.linspace(0, 1, nmini + 1)[1:]:
        b = step(b, np.exp(lU + s*e1), np.exp(lV + s*e2))
    vals = phi_of(b, kind)
    F = np.fft.fft2(vals, axes=(-2, -1))/(nc*nc)                  # coefficients c_{mn} rho^{m+n}
    from math import factorial
    D = {}
    for m in range(order + 1):
        for n in range(order + 1 - m):
            D[(m, n)] = F[..., m, n]*factorial(m)*factorial(n)/rho**(m + n)
    return D


# ------------------------------------------------------------------ term algebra: {(e12,e13,e14,e23,e24,e34,m,n): coef}
def I1111():
    e = [0]*6; e[PIDX[(1, 3)]] = -1; e[PIDX[(2, 4)]] = -1
    return {tuple(e) + (0, 0): 0.5}


def dX(f, pair):
    p = PIDX[pair]; out = {}
    for key, c in f.items():
        e = list(key[:6]); m, n = key[6:]
        ep = e[p]; e2 = list(e); e2[p] -= 1
        for (cc, mm, nn) in ((ep, m, n), (AU[p], m + 1, n), (BV[p], m, n + 1)):
            if cc:
                k2 = tuple(e2) + (mm, nn)
                out[k2] = out.get(k2, 0.0) + c*cc
    return {k: v for k, v in out.items() if v != 0}


def raise_pair(f, pair):
    return scale(dX(f, pair), -1.0)


def theta(f, pair):
    """theta_ij = -X_ij d/dX_ij"""
    p = PIDX[pair]; g = dX(f, pair); out = {}
    for key, c in g.items():
        e = list(key[:6]); e[p] += 1
        k2 = tuple(e) + key[6:]
        out[k2] = out.get(k2, 0.0) - c
    return {k: v for k, v in out.items() if v != 0}


def mulX(f, pair, power=1):
    p = PIDX[pair]; out = {}
    for key, c in f.items():
        e = list(key[:6]); e[p] += power
        out[tuple(e) + key[6:]] = c
    return out


def add(*fs, coefs=None):
    out = {}
    coefs = coefs or [1.0]*len(fs)
    for f, a in zip(fs, coefs):
        for k, v in f.items():
            out[k] = out.get(k, 0.0) + a*v
    return {k: v for k, v in out.items() if v != 0}


def scale(f, a):
    return {k: a*v for k, v in f.items()}


def apply_poly(f, poly):
    """poly: list of (coef, [pairs...]) = sum coef * prod theta_pair (theta's commute).  Returns P(theta) f."""
    out = {}
    for c, prs in poly:
        g = f
        for pr in prs:
            g = theta(g, pr)
        out = add(out, scale(g, c))
    return out


def max_order(f):
    return max(k[6] + k[7] for k in f) if f else 0


def evaluate(f, Xs, D):
    """Xs: dict pair -> array (effective X_ij), or callable for exact powers; D: derivative table.  Returns sum."""
    tot = 0
    for key, c in f.items():
        term = c*D[(key[6], key[7])]
        for p, ex in zip(PAIRS, key[:6]):
            if ex:
                term = term*Xs[p]**ex
        tot = tot + term
    return tot


def Xgauge(t1, t3, y):
    return {(1, 2): 2*t1 - 1, (3, 4): 2*t3 - 1, (1, 3): 2*y*t1*t3, (1, 4): 2*t1 + 1, (2, 3): 2*t3 + 1,
            (2, 4): 4.0 + 0*t1*t3}


# ------------------------------------------------------------------ the test correlators (identical scalars, Delta = 2)
def base2222():
    return raise_pair(raise_pair(I1111(), (1, 2)), (3, 4))


SYM2 = [(1.0, [(1, 2), (1, 2)]), (1.0, [(1, 3), (1, 3)]), (1.0, [(1, 4), (1, 4)])]       # sigma2 = d12^2 + d13^2 + d14^2
SYM3 = [(1.0, [(1, 2), (1, 3), (1, 4)])]                                              # sigma3 = d12 d13 d14


def sq(poly):
    out = []
    for c1, p1 in poly:
        for c2, p2 in poly:
            out.append((c1*c2, p1 + p2))
    return out


def exchange_pieces():
    """massive spin-2 exchange, Delta_e = 4 (tau = 2), between Delta = 2 scalars (truncates at the leading pole), three channels.
    s: M_s = Q(delta)/(delta12 - 1), Q = (delta13 - delta14)^2/4 + 1/12 (spin-2 continuous-Hahn polynomial at the pole), etc.
    Position space: Q(theta) [X12^{-1} I_{1122}]."""
    def Q(pa, pb):        # (theta_a - theta_b)^2/4 + 1/12
        return [(0.25, [pa, pa]), (0.25, [pb, pb]), (-0.5, [pa, pb]), (1.0/12, [])]
    Ps = mulX(raise_pair(I1111(), (3, 4)), (1, 2), -1)        # X12^-1 I_{1122}
    Pt = mulX(raise_pair(I1111(), (2, 3)), (1, 4), -1)        # X14^-1 I_{1221}
    Pu = mulX(raise_pair(I1111(), (2, 4)), (1, 3), -1)        # X13^-1 I_{1212}
    return {'s': apply_poly(Ps, Q((1, 3), (1, 4))), 't': apply_poly(Pt, Q((1, 3), (1, 2))),
            'u': apply_poly(Pu, Q((1, 2), (1, 4)))}


def pieces():
    B = base2222()
    ex = exchange_pieces()
    return {'c0': B, 'c2': apply_poly(B, SYM2), 'c3': apply_poly(B, SYM3), 'c4': apply_poly(B, sq(SYM2)),
            'xs': ex['s'], 'xt': ex['t'], 'xu': ex['u']}


# ------------------------------------------------------------------ grid evaluation on C x C
def grid_D(t, y, order, kind='sv', chunk=10, eps=1e-3, nsteps=2500, rho=0.5, nc=20):
    """Phi-derivative tables on the grid t1 = t[:,None], t3 = t[None,:]."""
    n = len(t)
    a = 1/t; b = 1/t
    D = {(m, k): np.empty((n, n), complex) for m in range(order + 1) for k in range(order + 1 - m)}
    minsep = np.empty((n, n))
    for i0 in range(0, n, chunk):
        A = a[i0:i0 + chunk, None]; Bb = b[None, :]
        st, ms = continue_ab(A, Bb, y, eps=eps, nsteps=nsteps)
        U, V = UV_of(*np.broadcast_arrays(A, Bb), y)
        Dc = phi_derivs(st, U, V, order, rho=rho, nc=nc, kind=kind)
        for key in D:
            D[key][i0:i0 + chunk] = Dc[key]
        minsep[i0:i0 + chunk] = ms
    return D, minsep


# ------------------------------------------------------------------ continuation ALONG THE CONTOUR (physical paths)
def _line_path(alpha_target, nsteps):
    """alpha(s) from 0 to alpha_target on Re a = 1, uniform in asinh(alpha)."""
    A = np.arcsinh(alpha_target)
    return [np.sinh(s*A) for s in np.linspace(0, 1, nsteps + 1)[1:]]


def continue_contour(a, b, y, order='ab', eps=1e-3, n_real=600, n_line=3000, bump=0.05):
    """Phi state at (a, b) on Re a = Re b = 1 (the contour C in 1/t), reached by the physical route:
    (i) both variables along the totally spacelike real axis from t = inf (a = eps, Euclidean) to t = 1 (a = 1), with a small imaginary
        bump i*bump*sin(pi s) to avoid the real z = zb collision;  (ii) the first variable (order[0]) along C from t = 1 to its node, the
        other held at t = 1;  (iii) the second variable along C from t = 1 to its node.
    a, b broadcast arrays with Re = 1."""
    a = np.asarray(a, complex); b = np.asarray(b, complex)
    a, b = np.broadcast_arrays(a, b)
    U0, V0 = UV_of(eps, eps, y)
    st = principal_state(np.full(a.shape, U0, complex), np.full(a.shape, V0, complex))
    minsep = np.abs(st.z - st.zb)
    for s in np.linspace(0, 1, n_real + 1)[1:]:
        r = eps + (1 - eps)*s + 1j*bump*np.sin(np.pi*s)
        st = step(st, *UV_of(r + 0*a, r + 0*b, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    al = (a - 1).imag; be = (b - 1).imag
    A_path = np.array(_line_path(al, n_line)); B_path = np.array(_line_path(be, n_line))
    if order == 'ab':
        for x in A_path:
            st = step(st, *UV_of(1 + 1j*x, 1 + 0*b, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
        for x in B_path:
            st = step(st, *UV_of(a, 1 + 1j*x, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    else:
        for x in B_path:
            st = step(st, *UV_of(1 + 0*a, 1 + 1j*x, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
        for x in A_path:
            st = step(st, *UV_of(1 + 1j*x, b, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    return st, minsep


def grid_D_contour(t, y, order, kind='sv', chunk=10, corder='ab', rho=0.5, nc=20, **kw):
    """as grid_D, but with the contour route (continue_contour).  t must lie on |t - 1/2| = 1/2."""
    n = len(t)
    a = 1/t; b = 1/t
    D = {(m, k): np.empty((n, n), complex) for m in range(order + 1) for k in range(order + 1 - m)}
    minsep = np.empty((n, n))
    for i0 in range(0, n, chunk):
        A = a[i0:i0 + chunk, None]; Bb = b[None, :]
        st, ms = continue_contour(A, Bb, y, order=corder, **kw)
        U, V = UV_of(*np.broadcast_arrays(A, Bb), y)
        Dc = phi_derivs(st, U, V, order, rho=rho, nc=nc, kind=kind)
        for key in D:
            D[key][i0:i0 + chunk] = Dc[key]
        minsep[i0:i0 + chunk] = ms
    return D, minsep


# ------------------------------------------------------------------ hairpin family: ellipses through t = 0 and t = 1
def contour_ellipse(n, kappa=1.0):
    """nodes/weights on t(theta) = (1 - cos theta)/2 - i (kappa/2) sin theta, theta = pi (1 - cos psi) (kappa = 1: og_free.contour).
    Encloses the cut (0, 1/2]; kappa -> 0 collapses onto the double-commutator segment (the hairpin)."""
    x, wx = np.polynomial.legendre.leggauss(n)
    psi = (x + 1)*np.pi/2; wpsi = wx*np.pi/2
    th = np.pi*(1 - np.cos(psi)); dth = np.pi*np.sin(psi)
    t = (1 - np.cos(th))/2 - 0.5j*kappa*np.sin(th)
    dt = (np.sin(th)/2 - 0.5j*kappa*np.cos(th))*dth
    return t, wpsi*dt, th


def _theta_path(th_n, nsteps):
    """theta(s) from pi to th_n, log-uniform in the distance to the endpoint theta = 0 (or 2 pi)."""
    d_n = np.minimum(th_n, 2*np.pi - th_n)
    lo = np.log(d_n); hi = np.log(np.pi)
    out = []
    for s in np.linspace(0, 1, nsteps + 1)[1:]:
        d = np.exp(hi + s*(lo - hi))
        out.append(np.where(th_n < np.pi, d, 2*np.pi - d))
    return out


def continue_ellipse(th1, th3, y, kappa=1.0, order='ab', eps=1e-3, n_real=600, n_line=3000, bump=0.05):
    """physical route: real spacelike axis t = inf -> 1 (both), then the first variable along the ellipse from t = 1 to its node, then
    the second.  th1, th3: broadcast arrays of node angles."""
    th1 = np.asarray(th1, float); th3 = np.asarray(th3, float)
    th1, th3 = np.broadcast_arrays(th1, th3)
    tt = lambda th: (1 - np.cos(th))/2 - 0.5j*kappa*np.sin(th)
    U0, V0 = UV_of(eps, eps, y)
    st = principal_state(np.full(th1.shape, U0, complex), np.full(th1.shape, V0, complex))
    minsep = np.abs(st.z - st.zb)
    one = np.ones(th1.shape, complex)
    for s in np.linspace(0, 1, n_real + 1)[1:]:
        r = eps + (1 - eps)*s + 1j*bump*np.sin(np.pi*s)
        st = step(st, *UV_of(r*one, r*one, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    P1 = _theta_path(th1, n_line); P3 = _theta_path(th3, n_line)
    if order == 'ab':
        for th in P1:
            st = step(st, *UV_of(1/tt(th), one, y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
        for th in P3:
            st = step(st, *UV_of(1/tt(th1), 1/tt(th), y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    else:
        for th in P3:
            st = step(st, *UV_of(one, 1/tt(th), y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
        for th in P1:
            st = step(st, *UV_of(1/tt(th), 1/tt(th3), y)); minsep = np.minimum(minsep, np.abs(st.z - st.zb))
    return st, minsep


def grid_D_ellipse(n, y, order, kappa=1.0, kind='sv', chunk=10, corder='ab', rho=0.5, nc=20, **kw):
    t, w, th = contour_ellipse(n, kappa)
    D = {(m, k): np.empty((n, n), complex) for m in range(order + 1) for k in range(order + 1 - m)}
    minsep = np.empty((n, n))
    for i0 in range(0, n, chunk):
        st, ms = continue_ellipse(th[i0:i0 + chunk, None], th[None, :], y, kappa=kappa, order=corder, **kw)
        A = 1/t[i0:i0 + chunk, None]; Bb = 1/t[None, :]
        U, V = UV_of(*np.broadcast_arrays(A, Bb), y)
        Dc = phi_derivs(st, U, V, order, rho=rho, nc=nc, kind=kind)
        for key in D:
            D[key][i0:i0 + chunk] = Dc[key]
        minsep[i0:i0 + chunk] = ms
    return D, minsep
