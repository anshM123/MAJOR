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


# ------------------------------------------------------------------ exact log-derivatives of Phi via (z, zb) Leibniz
# The (logU, logV) polydisc of phi_derivs is ill-conditioned at the light-ray endpoint (U ~ V >> 1, z ~ U - V cancels); this route is
# exact: every block depends on z or on zb only, its z-derivatives are elementary, and D_U^m D_V^n = sum_ij c_ij(z, zb) d_z^i d_zb^j with
# rational c_ij obtained once with sympy (acting on a generic f(z, zb)).
from math import comb, factorial
_OPC = {}


def _padd(p, q, c=1):
    for k, v in q.items():
        p[k] = p.get(k, 0) + c*v
        if p[k] == 0:
            del p[k]
    return p


def _pmul(p, q):
    out = {}
    for (a, b, c), u in p.items():
        for (d, e, f), v in q.items():
            k = (a + d, b + e, c + f)
            out[k] = out.get(k, 0) + u*v
    return {k: v for k, v in out.items() if v != 0}


def _pdz(p):            # d/dz on polynomials in (z, zb, w), w = 1/(z - zb): dw/dz = -w^2
    out = {}
    for (a, b, c), u in p.items():
        if a:
            k = (a - 1, b, c); out[k] = out.get(k, 0) + a*u
        if c:
            k = (a, b, c + 1); out[k] = out.get(k, 0) - c*u
    return {k: v for k, v in out.items() if v != 0}


def _pdzb(p):           # dw/dzb = +w^2
    out = {}
    for (a, b, c), u in p.items():
        if b:
            k = (a, b - 1, c); out[k] = out.get(k, 0) + b*u
        if c:
            k = (a, b, c + 1); out[k] = out.get(k, 0) + c*u
    return {k: v for k, v in out.items() if v != 0}


def _opcoefs(order):
    """D_U^m D_V^n = sum_ij c_ij d_z^i d_zb^j, c_ij integer polynomials in (z, zb, w = 1/(z - zb)) (exact)."""
    if order in _OPC:
        return _OPC[order]
    # D_U = A_z d_z + A_zb d_zb,  A_z = -z zb (1 - z) w,  A_zb = z zb (1 - zb) w ;  D_V: B_z = -z^2 (1-z)(1-zb) w... (expanded below)
    one = {(0, 0, 0): 1}
    Z = {(1, 0, 0): 1}; ZB = {(0, 1, 0): 1}; W = {(0, 0, 1): 1}
    omz = _padd(dict(one), Z, -1); omzb = _padd(dict(one), ZB, -1)
    A_z = _pmul(_pmul(_pmul(Z, ZB), omz), W); A_z = {k: -v for k, v in A_z.items()}
    A_zb = _pmul(_pmul(_pmul(Z, ZB), omzb), W)
    V = _pmul(omz, omzb)
    B_z = _pmul(_pmul(V, Z), W); B_z = {k: -v for k, v in B_z.items()}
    B_zb = _pmul(_pmul(V, ZB), W)

    def apply(op, cz, czb):
        out = {}
        for (i, j), c in op.items():
            dc = _padd(_pmul(cz, _pdz(c)), _pmul(czb, _pdzb(c)))
            if dc:
                out[(i, j)] = _padd(out.get((i, j), {}), dc)
            out[(i + 1, j)] = _padd(out.get((i + 1, j), {}), _pmul(c, cz))
            out[(i, j + 1)] = _padd(out.get((i, j + 1), {}), _pmul(c, czb))
        return {k: v for k, v in out.items() if v}
    res = {}
    colm = {(0, 0): one}
    for m in range(order + 1):
        if m > 0:
            colm = apply(colm, A_z, A_zb)
        g = colm
        for n in range(order + 1 - m):
            if n > 0:
                g = apply(g, B_z, B_zb)
            res[(m, n)] = g
    _OPC[order] = res
    return res


def _peval(p, z, zb, w, cache):
    tot = 0
    for (a, b, c), u in p.items():
        for key, base in ((('z', a), z), (('zb', b), zb), (('w', c), w)):
            if key not in cache:
                cache[key] = base**key[1]
        tot = tot + u*cache[('z', a)]*cache[('zb', b)]*cache[('w', c)]
    return tot


def _block_derivs(w, lw, l1w, L2w, K):
    """k-th derivatives (k = 0..K) of log w, log(1-w), Li2(w) on the tracked branch."""
    dl = [lw] + [(-1)**(k - 1)*factorial(k - 1)/w**k for k in range(1, K + 1)]
    dl1 = [l1w] + [-factorial(k - 1)/(1 - w)**k for k in range(1, K + 1)]
    inv = [(-1)**k*factorial(k)/w**(k + 1) for k in range(K + 1)]          # d^k (1/w)
    dL2 = [L2w]
    for k in range(1, K + 1):                                               # d^k Li2 = d^(k-1)[-l1 * (1/w)]
        q = k - 1
        dL2.append(-sum(comb(q, r)*dl1[r]*inv[q - r] for r in range(q + 1)))
    return dl, dl1, dL2


def _prod_derivs(a, b, K):
    return [sum(comb(k, r)*a[r]*b[k - r] for r in range(k + 1)) for k in range(K + 1)]


def _dphi(st, K, kind='sv'):
    dlz, dl1z, dL2z = _block_derivs(st.z, st.lz, st.l1z, st.L2z, K)
    dlzb, dl1zb, dL2zb = _block_derivs(st.zb, st.lzb, st.l1zb, st.L2zb, K)
    zz = [0*st.z]*(K + 1)
    one = [1 + 0*st.z] + zz[1:]
    pz = _prod_derivs(dlz, dl1z, K); pzb = _prod_derivs(dlzb, dl1zb, K)
    # N = 2 L2z - 2 L2zb + lz l1z - lz l1zb + lzb l1z - lzb l1zb ;  dN[i][j]
    dN = [[None]*(K + 1) for _ in range(K + 1)]
    for i in range(K + 1):
        for j in range(K + 1 - i):
            sv = 1.0 if kind == 'sv' else 0.0
            v = 0
            if j == 0:
                v = v + 2*dL2z[i] + sv*pz[i]
            if i == 0:
                v = v - 2*dL2zb[j] - sv*pzb[j]
            v = v + sv*(-dlz[i]*dl1zb[j] + dlzb[j]*dl1z[i])
            dN[i][j] = v
    r = st.z - st.zb
    dR = lambda i, j: (-1)**i*factorial(i + j)/r**(1 + i + j)
    dPhi = {}
    for i in range(K + 1):
        for j in range(K + 1 - i):
            dPhi[(i, j)] = sum(comb(i, a)*comb(j, b)*dN[a][b]*dR(i - a, j - b) for a in range(i + 1) for b in range(j + 1))
    return dPhi


def phi_derivs_exact(st, order):
    dPhi = _dphi(st, order); r = st.z - st.zb
    C = _opcoefs(order)
    w = 1/r; cache = {}
    out = {}
    for mn, co in C.items():
        out[mn] = sum(_peval(c, st.z, st.zb, w, cache)*dPhi[ij] for ij, c in co.items())
    return out


# ------------------------------------------------------------------ jet (local Taylor) arithmetic: well-conditioned D_U, D_V
def _jmul(A, B, K):
    S = A.shape[:-2]
    C = np.zeros(S + (K + 1, K + 1), complex)
    for i in range(K + 1):
        for j in range(K + 1 - i):
            acc = 0
            for a in range(i + 1):
                for b in range(j + 1):
                    acc = acc + A[..., a, b]*B[..., i - a, j - b]
            C[..., i, j] = acc
    return C


def _jdz(A, K):
    C = np.zeros_like(A)
    for i in range(K):
        for j in range(K - i):
            C[..., i, j] = (i + 1)*A[..., i + 1, j]
    return C


def _jdzb(A, K):
    C = np.zeros_like(A)
    for i in range(K):
        for j in range(K - i):
            C[..., i, j] = (j + 1)*A[..., i, j + 1]
    return C


def phi_derivs_jet(st, order, kind='sv'):
    """D_U^m D_V^n Phi, m + n <= order, by local Taylor-jet arithmetic in (z, zb) (no global polynomial expansion)."""
    K = order
    dPhi = _dphi(st, K, kind)
    S = st.z.shape
    J = np.zeros(S + (K + 1, K + 1), complex)
    for (i, j), v in dPhi.items():
        J[..., i, j] = v/(factorial(i)*factorial(j))

    def var(c0, i, j):
        A = np.zeros(S + (K + 1, K + 1), complex); A[..., 0, 0] = c0
        if i + j:
            A[..., i, j] = 1
        return A
    Z = var(st.z, 1, 0); ZB = var(st.zb, 0, 1)
    ONE = var(1.0, 0, 0)
    r0 = st.z - st.zb
    Wj = np.zeros(S + (K + 1, K + 1), complex)
    for i in range(K + 1):
        for j in range(K + 1 - i):
            k = i + j
            Wj[..., i, j] = (-1)**k*comb(k, i)*(-1)**j/r0**(k + 1)
    UZ = _jmul(Z, ZB, K); VZ = _jmul(ONE - Z, ONE - ZB, K)
    Az = -_jmul(_jmul(UZ, ONE - Z, K), Wj, K); Azb = _jmul(_jmul(UZ, ONE - ZB, K), Wj, K)
    Bz = -_jmul(_jmul(VZ, Z, K), Wj, K); Bzb = _jmul(_jmul(VZ, ZB, K), Wj, K)
    DU = lambda F: _jmul(Az, _jdz(F, K), K) + _jmul(Azb, _jdzb(F, K), K)
    DV = lambda F: _jmul(Bz, _jdz(F, K), K) + _jmul(Bzb, _jdzb(F, K), K)
    out = {}
    col = J
    for m in range(K + 1):
        if m > 0:
            col = DU(col)
        g = col
        for n in range(K + 1 - m):
            if n > 0:
                g = DV(g)
            out[(m, n)] = g[..., 0, 0].copy()
    return out


# ------------------------------------------------------------------ high-precision fallback (mpmath) for ill-conditioned nodes
def _mp_state_node(t1, t3, y, st, i):
    """mpmath (z, zb, blocks) at one node, with the branches of the float-tracked state st[i]."""
    import mpmath as mp
    a = 1/mp.mpc(complex(t1)); b = 1/mp.mpc(complex(t3)); yy = mp.mpf(float(y))
    U = (2 - a)*(2 - b)/(8*yy); V = (2 + a)*(2 + b)/(8*yy); S = 1 + U - V; d = mp.sqrt(S*S - 4*U)
    r1 = (S + d)/2; r2 = (S - d)/2
    zf = st.z[i]
    z, zb = (r1, r2) if abs(complex(r1) - zf) <= abs(complex(r2) - zf) else (r2, r1)

    def br(tr, p):
        return p + 2j*mp.pi*int(round(((tr - complex(p))/(2j*np.pi)).real))

    def l2(tr, w):
        p = mp.polylog(2, w); D = tr - complex(p); c = 2j*np.pi*complex(mp.log(w))
        m = int(round(D.imag/c.imag)) if abs(c.imag) > 1e-12 else 0
        l = int(round((D.real - m*c.real)/(4*np.pi**2)))
        return p + m*2j*mp.pi*mp.log(w) + 4*mp.pi**2*l
    return (z, zb, br(st.lz[i], mp.log(z)), br(st.lzb[i], mp.log(zb)), br(st.l1z[i], mp.log(1 - z)), br(st.l1zb[i], mp.log(1 - zb)),
            l2(st.L2z[i], z), l2(st.L2zb[i], zb))


def _mp_jet_derivs(blk, K, kind='sv'):
    """same algorithm as phi_derivs_jet, scalar mpmath."""
    import mpmath as mp
    z, zb, lz, lzb, l1z, l1zb, L2z, L2zb = blk
    fac = [mp.factorial(k) for k in range(K + 2)]

    def bd(w, lw, l1w, L2w):
        dl = [lw] + [(-1)**(k - 1)*fac[k - 1]/w**k for k in range(1, K + 1)]
        dl1 = [l1w] + [-fac[k - 1]/(1 - w)**k for k in range(1, K + 1)]
        inv = [(-1)**k*fac[k]/w**(k + 1) for k in range(K + 1)]
        dL2 = [L2w] + [-sum(mp.binomial(k - 1, r)*dl1[r]*inv[k - 1 - r] for r in range(k)) for k in range(1, K + 1)]
        return dl, dl1, dL2
    dlz, dl1z, dL2z = bd(z, lz, l1z, L2z); dlzb, dl1zb, dL2zb = bd(zb, lzb, l1zb, L2zb)
    pz = [sum(mp.binomial(k, r)*dlz[r]*dl1z[k - r] for r in range(k + 1)) for k in range(K + 1)]
    pzb = [sum(mp.binomial(k, r)*dlzb[r]*dl1zb[k - r] for r in range(k + 1)) for k in range(K + 1)]
    r0 = z - zb
    dN = {}
    for i in range(K + 1):
        for j in range(K + 1 - i):
            sv = 1 if kind == 'sv' else 0
            v = sv*(-dlz[i]*dl1zb[j] + dlzb[j]*dl1z[i])
            if j == 0:
                v += 2*dL2z[i] + sv*pz[i]
            if i == 0:
                v -= 2*dL2zb[j] + sv*pzb[j]
            dN[(i, j)] = v
    idx = [(i, j) for i in range(K + 1) for j in range(K + 1 - i)]
    J = {}
    for (i, j) in idx:
        d = sum(mp.binomial(i, a)*mp.binomial(j, b)*dN[(a, b)]*(-1)**(i - a)*fac[i - a + j - b]/r0**(1 + i - a + j - b)
                for a in range(i + 1) for b in range(j + 1))
        J[(i, j)] = d/(fac[i]*fac[j])

    def mul(A, B):
        return {(i, j): sum(A.get((a, b), 0)*B.get((i - a, j - b), 0) for a in range(i + 1) for b in range(j + 1)) for (i, j) in idx}

    def dzj(A):
        return {(i, j): (i + 1)*A.get((i + 1, j), 0) for (i, j) in idx}

    def dzbj(A):
        return {(i, j): (j + 1)*A.get((i, j + 1), 0) for (i, j) in idx}
    Z = {(0, 0): z, (1, 0): 1}; ZB = {(0, 0): zb, (0, 1): 1}
    omZ = {(0, 0): 1 - z, (1, 0): -1}; omZB = {(0, 0): 1 - zb, (0, 1): -1}
    W = {(i, j): (-1)**(i + j)*mp.binomial(i + j, i)*(-1)**j/r0**(i + j + 1) for (i, j) in idx}
    UZ = mul(Z, ZB); VZ = mul(omZ, omZB)
    neg = lambda A: {k: -v for k, v in A.items()}
    Az = neg(mul(mul(UZ, omZ), W)); Azb = mul(mul(UZ, omZB), W)
    Bz = neg(mul(mul(VZ, Z), W)); Bzb = mul(mul(VZ, ZB), W)
    add = lambda A, B: {k: A.get(k, 0) + B.get(k, 0) for k in idx}
    DU = lambda F: add(mul(Az, dzj(F)), mul(Azb, dzbj(F)))
    DV = lambda F: add(mul(Bz, dzj(F)), mul(Bzb, dzbj(F)))
    out = {}
    col = J
    for m in range(K + 1):
        if m > 0:
            col = DU(col)
        g = col
        for n in range(K + 1 - m):
            if n > 0:
                g = DV(g)
            out[(m, n)] = complex(g[(0, 0)])
    return out


def phi_derivs_mixed(st, t1, t3, y, order, zmax=0.0, sepmin=1e9, dps=32, kind="sv"):
    """D_U^m D_V^n Phi by jets.  float64 jets lose digits already at |z| ~ 5 for order >= 3 (cancellation between the d_z and d_zb parts
    of D_U, D_V; found at dphi = 0.44, RESULT_Q_EINSTEIN_CHECK.md §2.4), so by default EVERY node is recomputed in mpmath with
    dps + 4 log10 max(|z|, |zb|, 1) digits (zmax = 0).  Returns (table, number of hp nodes)."""
    import mpmath as mp
    D = phi_derivs_jet(st, order, kind)
    t1, t3 = np.broadcast_arrays(np.asarray(t1, complex), np.asarray(t3, complex))
    flag = (np.abs(st.z) > zmax) | (np.abs(st.zb) > zmax) | (np.abs(st.z - st.zb) < sepmin)
    idxs = np.argwhere(flag)
    for ix in idxs:
        ix = tuple(ix)
        big = max(abs(st.z[ix]), abs(st.zb[ix]), 1.0)
        with mp.workdps(dps + int(4*np.log10(big))):
            sub = State()
            for k in State.__slots__:
                setattr(sub, k, getattr(st, k)[ix])
            blk = _mp_state_node(t1[ix], t3[ix], y, State_wrap(sub), 0)
            r = _mp_jet_derivs(blk, order, kind)
        for key in D:
            D[key][ix] = r[key]
    return D, len(idxs)


class State_wrap:
    """index-0 view of a scalar State for _mp_state_node."""
    def __init__(self, s):
        for k in State.__slots__:
            setattr(self, k, np.atleast_1d(getattr(s, k)))


def grid_D_ellipse_hp(n, y, order, kappa=1.0, chunk=10, corder='ab', kind='sv', **kw):
    t, w, th = contour_ellipse(n, kappa)
    D = {(m, k): np.empty((n, n), complex) for m in range(order + 1) for k in range(order + 1 - m)}
    nhp = 0
    minsep = np.empty((n, n))
    for i0 in range(0, n, chunk):
        st, ms = continue_ellipse(th[i0:i0 + chunk, None], th[None, :], y, kappa=kappa, order=corder)
        T1, T3 = np.broadcast_arrays(t[i0:i0 + chunk, None], t[None, :])
        Dc, k = phi_derivs_mixed(st, T1, T3, y, order, kind=kind, **kw)
        nhp += k
        for key in D:
            D[key][i0:i0 + chunk] = Dc[key]
        minsep[i0:i0 + chunk] = ms
    return D, minsep, nhp
