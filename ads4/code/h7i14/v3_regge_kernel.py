"""v3 (H7_I14_LIGHTSIDE.md §18): the leading-Regge kernel of the CLS-D spinning functional Psi_{2,nu}[q^2] on <TTTT>, its heavy
(plane-wave) image, and the decomposition into P_{l-1}, P_l, P_{l+1} per helicity sector.

Derivation (§18.1).  CLS-D Psi^+ (eq. spinningfunctional_definition, the integrand of H0 §21.3 / l6_spin_action) in the CHMRS frame (3.23):
x1 = -y, x3 = y, x2 = x/x^2, x4 = -x/x^2 (Regge limit x, y -> 0, r = |x||y|).  The delta functions put x0 on the circle
x0 = |y| n(theta), n unit spacelike, n.y = 0; the light-ray polarisations are z1 = x_{1+0} = -|y|(yh + n), z3 = x30 = |y|(yh - n).
At leading Regge order the (2,4) dual structure is evaluated at x0 = 0 (corrections O(|x0|/|x2|) = O(r)).  With polarisations of 2, 4
referred to the origin by the inversion (z' = I(x) z, Lemma V of REGGE_MOMENT_TTTT) one finds exactly (x0 = 0, |x| = 1):
   V_{0,42} = (z0.xh)/2,  H20 = z2'.z0,  H40 = z4'.z0,  X20 = X40 = 1,  H24 = -4 z2'.I(xh).z4'   (x^2 = -1),
   q = 2 (z2'.z0)(z4'.z0) - (z0.xh)^2 (z2'.I(xh).z4'),   I(xh) = eta + 2 xh xh.
So the f-stripped leading Regge kernel factorises through the celestial point z0 (split representation):
   K0(xh, yh)[tau24, tau13] = int Dz0 A(xh, z0)[tau24] B(yh, z0)[tau13],
   A = (-z0.xh)^{-(1-delta)} * 16 q^2/(z0.xh)^4                            (scalar dual x q^2/V^4; delta = 1/2 + i nu)
   B = (-z0.yh)^{-delta} int_0^{2pi} dth |sin th|^{-delta} tau13(yh + n_th, yh - n_th)  (celestial j' = 0, delta_{1,3} = 2; th measured
                                                                                    from the direction m of z0 at yh)
up to nu-only constants and nu-only phases (all bases have a fixed sign on the whole integration region).  The scalar (Delta_phi = 3)
functional is the same with q^2/V^4 -> 1, tau13 -> 1.  Radial powers: both pairs carry a = 7 (the TTTT factors |x|^4 from q^2/V^4,
|y|^4 from tau13(x10, x30) and r^-4 from f = v'-u' instead of (v'-u')/(u'v') cancel against each other; §18.1).  f = v'-u' = eta (1+r^2)/(2r).

Heavy image (Lemma V + Lemma D3(i)): plane-wave pair e^{p.x + pb.y} with constant transverse helicity tensors tau24 at ph (z2' = eps_{h2},
z4' = eps_{-h2}) and tau13 at pbh (z1 = eps_{h1}, z3 = eps_{-h1}).  With b = a - 1 = 6 (eta K = -(x.y)/(|x||y|) K):
   L(p, pb) = - eta^{mu nu} d/dp^mu d/dpb^nu  int Dz0  Ahat_b(p, z0) Bhat_b(pb, z0),
   Ahat_b(p, z0) = int d^3x |x|^{b-3} A(xh, z0) e^{p.x}   (closed form: every term is |x|^c (-z0.x)^{-lam} (x.e)(x.f)... ->
                    (e.d_p)(f.d_p)... F_{c,lam}(p),  F_{c,lam} = K(c,lam) (-z0.p)^{-lam} (-p^2)^{(2 lam - c - 3)/2},
                    K(c,lam) = 2^{mu-1} sqrt(pi) Gamma(1 + c/2) Gamma((c+3)/2 - lam), mu = c + 3 - lam).
The f-stripped kernel (no eta, b = a = 7) is computed the same way without the gradient.
The z0 integral (celestial circle, z0 = (1, cos phi, sin phi), measure dphi) is a smooth periodic integral: trapezoid rule.

Output: functions used by v3_hn_analysis.py; run directly for the validation tests (cache/v3_regge_kernel.txt).
"""
import os, sys, io, time, functools
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.special import loggamma

ETA = np.array([-1.0, 1.0, 1.0])


def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


# ------------------------------------------------------------------------------------------------ symbolic derivative patterns
# D[i, j][mu] = (e.d)^i (f.d)^j d/dp^mu [ S^{-lam} R^{kap} ],  S = -z0.p, R = -p.p   (gradient optional: mu = None)
_p = sp.symbols('p0:3'); _z = sp.symbols('z0:3'); _e = sp.symbols('e0:3'); _f = sp.symbols('f0:3')
_lam, _kap = sp.symbols('lam kap')
_md = lambda u, v: -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
_S = -_md(_z, _p); _R = -_md(_p, _p)
_Ssym, _Rsym = sp.symbols('S R', positive=True)


def _dir(expr, vec):
    return sum(vec[k]*sp.diff(expr, _p[k]) for k in range(3))


# truncated power series in (t1, t2, t3) with orders (2, 2, 1): coefficient arrays (..., 3, 3, 2)
def _tmul(a, b):
    out = np.zeros(np.broadcast_shapes(a.shape, b.shape), complex)
    for i in range(3):
        for j in range(3):
            for k in range(2):
                ai = a[..., i, j, k]
                if not np.any(ai):
                    continue
                for i2 in range(3 - i):
                    for j2 in range(3 - j):
                        for k2 in range(2 - k):
                            out[..., i + i2, j + j2, k + k2] += ai*b[..., i2, j2, k2]
    return out


def _tpow(b, alpha):
    """b^alpha, b = b0 + h (h nilpotent, total degree <= 5)."""
    b0 = b[..., 0, 0, 0]
    h = b/b0[..., None, None, None]; h[..., 0, 0, 0] = 0
    out = np.zeros_like(h); out[..., 0, 0, 0] = 1
    term = out.copy(); coef = 1.0 + 0j
    for k in range(1, 6):
        coef = coef*(alpha - k + 1)/k
        term = _tmul(term, h)
        out = out + coef*term
    return out*(b0**alpha)[..., None, None, None]




@functools.lru_cache(maxsize=None)
def _fact():
    import math
    return np.array([[[math.factorial(i)*math.factorial(j)*math.factorial(k) for k in range(2)] for j in range(3)]
                     for i in range(3)], float)


def taylor_SR(p, z, e, f, w, lam, kap):
    """Taylor coefficients x factorials of S^-lam R^kap at p + t1 e + t2 f + t3 w: returns array (..., 3, 3, 2) of
    (e.d)^i (f.d)^j (w.d)^k [S^-lam R^kap]."""
    sh = np.broadcast_shapes(p.shape, z.shape, e.shape, f.shape, w.shape)[:-1]
    S = np.zeros(sh + (3, 3, 2), complex); R = np.zeros_like(S)
    S[..., 0, 0, 0] = -mdot(z, p); S[..., 1, 0, 0] = -mdot(z, e); S[..., 0, 1, 0] = -mdot(z, f); S[..., 0, 0, 1] = -mdot(z, w)
    R[..., 0, 0, 0] = -mdot(p, p)
    R[..., 1, 0, 0] = -2*mdot(p, e); R[..., 0, 1, 0] = -2*mdot(p, f); R[..., 0, 0, 1] = -2*mdot(p, w)
    R[..., 2, 0, 0] = -mdot(e, e); R[..., 0, 2, 0] = -mdot(f, f)
    R[..., 1, 1, 0] = -2*mdot(e, f); R[..., 1, 0, 1] = -2*mdot(e, w); R[..., 0, 1, 1] = -2*mdot(f, w)
    return _tmul(_tpow(S, -lam), _tpow(R, kap))*_fact()


def pattern(i, j, grad):
    """callable(p, z, e, f, lam, kap) -> list of (e.d)^i (f.d)^j [d/dp^mu] (S^-lam R^kap) (3 comps if grad else 1).
    d/dp^mu = derivative along the coordinate basis vector E_mu (upper-index components)."""
    def call(p, z, e, f, lam, kap):
        if not grad:
            w = np.zeros(3)
            return [taylor_SR(p, z, e, f, w, lam, kap)[..., i, j, 0]]
        return [taylor_SR(p, z, e, f, np.eye(3)[mu], lam, kap)[..., i, j, 1] for mu in range(3)]
    return call


def logKc(c, lam):
    mu = c + 3 - lam
    return (mu - 1)*np.log(2.0) + 0.5*np.log(np.pi) + loggamma(1 + c/2) + loggamma((c + 3)/2 - lam)


def Kc(c, lam):
    return np.exp(logKc(c, lam))


# ------------------------------------------------------------------------------------------------ monomial expansions (exact, sympy)
_X, _s, _ex, _fx = sp.symbols('X s ex fx')          # |x|, -z0.x, eps2.x, eps4.x
_a2, _a4, _e24, _dl = sp.symbols('a2 a4 e24 delta')


@functools.lru_cache(maxsize=None)
def A_monomials(spin=True):
    """A(x, z0) |x|^{b-3} / |x|^{b-3} as a sum of coef * X^k s^m ex^i fx^j  (b enters c below).  Returns list of
    (coef_lambda(a2, a4, e24, delta), k_offset(delta-dependent sympy), m, i, j)."""
    if not spin:
        expr = _X**(1 - _dl)*_s**(-(1 - _dl))
    else:
        br = 2*_a2*_a4 - _s**2*_X**-2*_e24 - 2*_s**2*_X**-4*_ex*_fx
        expr = 16*_X**(1 - _dl + 4)*_s**(-(1 - _dl) - 4)*sp.expand(br**2)
    expr = sp.expand(expr)
    out = []
    for term in sp.Add.make_args(expr):
        pw = term.as_powers_dict()
        k = sp.simplify(pw.get(_X, 0)); m = sp.simplify(pw.get(_s, 0)); i = int(pw.get(_ex, 0)); j = int(pw.get(_fx, 0))
        coef = sp.simplify(term/(_X**k*_s**m*_ex**i*_fx**j))
        out.append((sp.lambdify((_a2, _a4, _e24, _dl), coef, 'numpy'), sp.lambdify(_dl, k, 'numpy'),
                    sp.lambdify(_dl, m, 'numpy'), i, j))
    return out


_Y, _t, _U1, _U3, _A1, _A3, _E13 = sp.symbols('Y t U1 U3 A1 A3 E13')
_c, _sn = sp.symbols('c sn')


@functools.lru_cache(maxsize=None)
def B_monomials(spin=True):
    """B(y, z0) as sum of coef(A1, A3, E13, delta) * Y^k t^m U1^i U3^j (U = eps.y, A = eps.z0, t = -z0.y, Y = |y|), including the
    theta-average with weight |sin|^{-delta}: <c^{2i} s^{2j}> = 2 B(i+1/2, j+1/2-delta/2) (B = Euler beta; symbolic via gamma)."""
    # moments divided by the common nu-only factor <1> = 2 B(1/2, 1/2 - delta/2) (dropped; it is the same for the scalar):
    # <c^{2i} s^{2j}>/<1> = (1/2)_i (1/2 - delta/2)_j / (1 - delta/2)_{i+j}   (rational in delta; no gamma overflow at large nu)
    momv = lambda i2, j2: sp.rf(sp.Rational(1, 2), i2)*sp.rf(sp.Rational(1, 2) - _dl/2, j2)/sp.rf(1 - _dl/2, i2 + j2)
    MS = {(i2, j2): sp.Symbol('M%d%d' % (i2, j2)) for i2 in range(3) for j2 in range(3)}   # kept symbolic until the end
    mom = lambda i2, j2: MS[(i2, j2)]
    if not spin:
        W = sp.Integer(1)
    else:
        u1, u3, m1, m3, t1, t3 = sp.symbols('u1 u3 m1 m3 t1 t3')
        P = sp.expand((u1 + _c*m1 + _sn*t1)**2*(u3 - _c*m3 - _sn*t3)**2)
        W = 0
        for term in sp.Add.make_args(P):
            pw = term.as_powers_dict(); ic = int(pw.get(_c, 0)); js = int(pw.get(_sn, 0))
            if ic % 2 or js % 2:
                continue
            W += term/(_c**ic*_sn**js)*mom(ic//2, js//2)
        W = sp.expand(W)
        # t-monomials: t1^2 -> pi11 - m1^2, t3^2 -> pi33 - m3^2, t1 t3 -> pi13 - m1 m3, t1^2 t3^2 -> (pi13 - m1 m3)^2
        pi11, pi33, pi13 = u1**2, u3**2, _E13 + u1*u3
        W2 = 0
        for term in sp.Add.make_args(W):
            pw = term.as_powers_dict(); a = int(pw.get(t1, 0)); b = int(pw.get(t3, 0))
            rest = term/(t1**a*t3**b)
            if (a, b) == (0, 0):
                rep = 1
            elif (a, b) == (2, 0):
                rep = pi11 - m1**2
            elif (a, b) == (0, 2):
                rep = pi33 - m3**2
            elif (a, b) == (1, 1):
                rep = pi13 - m1*m3
            elif (a, b) == (2, 2):
                rep = (pi13 - m1*m3)**2
            else:
                raise ValueError((a, b))
            W2 += rest*rep
        W = sp.expand(W2.subs({u1: _U1/_Y, u3: _U3/_Y, m1: _A1*_Y/_t - _U1/_Y, m3: _A3*_Y/_t - _U3/_Y}))
    # W is a Laurent polynomial in Y, t (integer powers): read it off with Poly after clearing negative powers
    Wp = sp.Poly(sp.expand(W*_Y**8*_t**8), _Y, _t, _U1, _U3)
    out = []
    for (kY, kt, i, j), coef in Wp.terms():
        k0 = kY - 8; m0 = kt - 8
        coef = sp.simplify(coef.subs({MS[key]: momv(*key) for key in MS}))
        k = k0 + _dl; m = m0 - _dl
        out.append((sp.lambdify((_A1, _A3, _E13, _dl), coef, ['numpy', 'scipy']), sp.lambdify(_dl, k, 'numpy'),
                    sp.lambdify(_dl, m, 'numpy'), i, j))
    return out


# ------------------------------------------------------------------------------------------------ transforms
def side_transform(p, z, e, f, monos, b, delta, coef_args, grad, lshift=0.0):
    """sum over monomials coef * X^{k} s^{m} (e.x)^i (f.x)^j times |x|^{b-3}: Laplace transform at p (vectorised over z).
    Monomials with equal (k, m) share one Taylor evaluation.  Returns array (..., 3) if grad else (...)."""
    groups = {}
    for (cf, kf, mf, i, j) in monos:
        key = (complex(kf(delta)), complex(mf(delta)))
        groups.setdefault(key, []).append((cf, i, j))
    tot = 0
    for (k, m), lst in groups.items():
        c = (b - 3) + k; lam = -m; kap = (2*lam - c - 3)/2
        Kv = np.exp(logKc(c, lam) - lshift)      # lshift: a common nu-only log-scale (avoids overflow at large nu)
        ws = [np.eye(3)[mu] for mu in range(3)] if grad else [np.zeros(3)]
        kidx = 1 if grad else 0
        comps = []
        for w in ws:
            Tm = taylor_SR(p, z, e, f, w, lam, kap)
            acc = 0
            for (cf, i, j) in lst:
                acc = acc + cf(*coef_args, delta)*Tm[..., i, j, kidx]
            comps.append(acc*Kv)
        tot = tot + np.stack(comps, -1)
    return tot if grad else tot[..., 0]


def frames(beta):
    """ph (pair 2,4) and pbh (pair 1,3) at distance beta; parallel-transported helicity vectors (Lemma V frame)."""
    c, s = np.cosh(beta/2), np.sinh(beta/2)
    ph = np.array([c, s, 0.0]); pbh = np.array([c, -s, 0.0])
    npar_p = np.array([s, c, 0.0]); npar_pb = np.array([-s, c, 0.0]); eperp = np.array([0.0, 0.0, 1.0])
    eps = lambda npar, h: (npar + 1j*h*eperp)/np.sqrt(2)
    return ph, pbh, eps(npar_p, 1), eps(npar_p, -1), eps(npar_pb, 1), eps(npar_pb, -1)


def heavy_image(nu, beta, h1, h2, spin=True, with_eta=True, nphi=256, a=7):
    """L^{(h1,h2)}(beta): heavy (plane-wave) image of the leading Regge kernel (eta K0 if with_eta else K0), preserving helicity
    tensors: (1,3) at pbh: z1 = eps_{h1}, z3 = eps_{-h1};  (2,4) at ph: z2' = eps_{h2}, z4' = eps_{-h2}.  nu-only constants dropped."""
    delta = 0.5 + 1j*nu
    ph, pbh, ep_p, em_p, ep_pb, em_pb = frames(beta)
    E = {(+1, 'p'): ep_p, (-1, 'p'): em_p, (+1, 'pb'): ep_pb, (-1, 'pb'): em_pb}
    e2, e4 = E[(h2, 'p')], E[(-h2, 'p')]
    e1, e3 = E[(h1, 'pb')], E[(-h1, 'pb')]
    phi = 2*np.pi*np.arange(nphi)/nphi
    z = np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)
    b = a - 1 if with_eta else a
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    E2 = np.broadcast_to(e2, z.shape); E4 = np.broadcast_to(e4, z.shape)
    E1 = np.broadcast_to(e1, z.shape); E3 = np.broadcast_to(e3, z.shape)
    if spin:
        cA = (mdot(z, E2), mdot(z, E4), np.full(nphi, mdot(e2, e4)))
        cB = (mdot(z, E1), mdot(z, E3), np.full(nphi, mdot(e1, e3)))
    else:
        cA = cB = (0*phi, 0*phi, 0*phi)
    # common nu-only rescaling (the same for spin/scalar and for both b): the scalar constants at b = 6
    lsA = logKc(6 - 2 - delta, 1 - delta).real; lsB = logKc(6 - 3 + delta, delta).real
    Ah = side_transform(P, z, E2, E4, A_monomials(spin), b, delta, cA, with_eta, lsA)
    Bh = side_transform(PB, z, E1, E3, B_monomials(spin), b, delta, cB, with_eta, lsB)
    if with_eta:
        integrand = -(Ah*Bh*ETA).sum(-1)        # - eta^{mu nu} dA/dp^mu dB/dpb^nu  (eta^{-1} = eta)
    else:
        integrand = Ah*Bh
    return integrand.mean()*2*np.pi


# ------------------------------------------------------------------------------------------------ conical functions and fit
def Pcon(nu, beta):
    """P_{-1/2 + i nu}(cosh beta), nu complex (mpmath)."""
    return complex(mp.legenp(-0.5 + 1j*nu, 0, mp.cosh(beta), type=3))


def fit3(nu, betas, vals):
    """least squares vals(beta) = hm P_{l-1} + h0 P_l + hp P_{l+1} (P_{l-1} = P at nu + i, P_{l+1} at nu - i); returns coefs, rel resid."""
    M = np.array([[Pcon(nu + 1j, b), Pcon(nu, b), Pcon(nu - 1j, b)] for b in betas])
    v = np.asarray(vals)
    sol, *_ = np.linalg.lstsq(M, v, rcond=None)
    res = np.linalg.norm(M @ sol - v)/np.linalg.norm(v)
    return sol, res


def fit1(nu, betas, vals):
    M = np.array([[Pcon(nu, b)] for b in betas]); v = np.asarray(vals)
    sol, *_ = np.linalg.lstsq(M, v, rcond=None)
    return sol[0], np.linalg.norm(M @ sol - v)/np.linalg.norm(v)


def gfac(nu):
    """CHMRS magic-block coefficient (l2_scalar_validate.gfac): g = ((2-d)/2 + i nu)/(i nu) (gamma_7(nu+i)/gamma_7(nu))^2."""
    lg = lambda v: loggamma((1 + 7 - 1.5 - 1j*v)/2) + loggamma((1 + 7 - 1.5 + 1j*v)/2)
    gp = (-0.5 + 1j*nu)/(1j*nu)*np.exp(2*(lg(nu + 1j) - lg(nu)))
    gm = (-0.5 - 1j*nu)/(-1j*nu)*np.exp(2*(lg(-nu + 1j) - lg(-nu)))
    return gp, gm


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    # (V0) transform constant K(c, lam) vs direct H2 quadrature
    pr("(V0) closed-form transform int d^3x |x|^c (-z0.x)^-lam e^{p.x} = K(c,lam)(-z0.p)^-lam (-p^2)^{(2lam-c-3)/2}: check of the angular"
       " constant 2 pi int_1^oo P_{lam-1}(t) t^{-mu} dt")
    for (c, lam) in ((3.2 - 0.5j, 0.5 + 0.7j), (4.0 + 0.3j, 2.5 - 0.3j), (8.0, 4.5 + 1.1j)):
        mu = c + 3 - lam
        num = 2*mp.pi*mp.quad(lambda t: mp.legenp(lam - 1, 0, t, type=3)*t**(-mu), [1, 2, 10, mp.inf])
        clo = 2*mp.pi*2**(mu - 2)*mp.gamma((mu + lam - 1)/2)*mp.gamma((mu - lam)/2)/(mp.sqrt(mp.pi)*mp.gamma(mu))
        pr(f"   c={c}, lam={lam}: rel diff {float(abs(num/clo - 1)):.1e}")
    pr("   full transform at a boosted p (p = (cosh .8, sinh .8 cos .3, sinh .8 sin .3), z0 = (1, cos 1.1, sin 1.1)), direct H2 quadrature x"
       " Gamma(mu):")
    from numpy.polynomial.legendre import leggauss
    xg, wg = leggauss(400); rho = 20*(xg + 1); wr = 20*wg; nph = 256; ph_ = 2*np.pi*np.arange(nph)/nph
    R_, PH_ = np.meshgrid(rho, ph_, indexing='ij')
    xh = np.stack([np.cosh(R_), np.sinh(R_)*np.cos(PH_), np.sinh(R_)*np.sin(PH_)], -1)
    pv = np.array([np.cosh(.8), np.sinh(.8)*np.cos(.3), np.sinh(.8)*np.sin(.3)]); zv = np.array([1, np.cos(1.1), np.sin(1.1)])
    for (c, lam) in ((3.2 - 0.5j, 0.5 + 0.7j), (4.0 + 0.3j, 2.5 - 0.3j)):
        mu = c + 3 - lam
        dirv = np.exp(loggamma(mu))*((wr*np.sinh(rho))[:, None]*(-mdot(xh, zv))**(-lam)*(-mdot(xh, pv))**(-mu)).sum()*2*np.pi/nph
        clo = Kc(c, lam)*(-mdot(zv, pv))**(-lam)*(-mdot(pv, pv))**((2*lam - c - 3)/2)
        pr(f"   c={c}, lam={lam}: rel diff {abs(dirv/clo - 1):.1e}")
    betas = np.linspace(0.15, 3.0, 14)
    # (V1) scalar functional: eta-kernel image must be prop. to g P_{l-1} + g~ P_{l+1}, no P_l; f-stripped image prop. to P_l.
    pr("\n(V1) scalar (Delta_phi = 3) functional through the same split representation:")
    for nu in (0.3, 1.0, 2.5, 6.0):
        v0 = [heavy_image(nu, b, 1, 1, spin=False, with_eta=False) for b in betas]
        c0, r0 = fit1(nu, betas, v0)
        v = [heavy_image(nu, b, 1, 1, spin=False, with_eta=True) for b in betas]
        (hm, h0, hp), res = fit3(nu, betas, v)
        gp, gm = gfac(nu)
        pr(f"   nu={nu}: f-stripped = c P_l, resid {r0:.1e};  eta-kernel: resid {res:.1e}, |h0|/|h-| = {abs(h0)/abs(hm):.1e}, "
           f"(h-/h+)/(g/g~) - 1 = {abs((hm/hp)/(gp/gm) - 1):.1e}, h-/(c g/2) - 1 = {abs(hm/(c0*gp/2) - 1):.1e}")
    pr(f"\n[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v3_regge_kernel.txt'), 'w').write(buf.getvalue())
