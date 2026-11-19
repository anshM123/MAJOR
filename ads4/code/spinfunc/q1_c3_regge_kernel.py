"""
q1 (S4 fix, S1_S4_FIXES.md §1): the leading-Regge kernel of the CLS-D spinning functional Psi_{3,nu}[o_H q] (j' = 1, f = -1) on <TTTT>
and its heavy (plane-wave) image -- the same split-representation method as h7i14/v3_regge_kernel.py (H7_I14_LIGHTSIDE.md §18) for C^_2.
v3_regge_kernel is imported read-only (transform engine taylor_SR / logKc, frames, conical functions, fits).

Derivation (frame (3.23) of v3: x1 = -y, x3 = y, x2 = x/x^2, x4 = -x/x^2; Regge limit r = |x||y| -> 0):
 * Dual side (operators 2, 4) at x0 = 0 (corrections O(r)).  In embedding space, with X0 = (1,0,0), Z0 = (0,0,z0), X4 = (1, x2^2, -x2):
       eps(Z2,X2,X4,X0,Z0) = -2 x2^2 eps3(z2, x2, z0),    eps(Z4,X4,X2,X0,Z0) = +2 x2^2 eps3(z4, x2, z0)        (row reduction, exact)
   and with the inversion-referred polarisations z' = R(xh) z (R = eta - 2 xh xh/xh^2, det R = -1, R xh = -xh):
       eps3(z2, xh, z0) = eps3(R z2', -R xh, R R z0) = eps3(z2', xh, R z0) = eps3(z2', xh, z0),   H20 ~ z2'.z0  (v3 §18.1).
   Hence, up to a constant,  o_H(x0 = 0) = (z2'.z0) eps3(z4', xh, z0) - (z4'.z0) eps3(z2', xh, z0)   (antisymmetric under 2 <-> 4),
   and the q of v3.  A(xh, z0) = (-z0.xh)^{-(1-delta)} 16 q o/(z0.xh)^4.   eps3(e, x, z0) = det[x; z0; e] = w_e . x,
   w_e := eta.cross(z0, e) (contravariant; Minkowski dot): a z0-dependent vector, handled by the third Taylor slot of v3's engine
   ((w.x) |x|^c (-z0.x)^-lam ... -> (w.d_p) F_{c,lam}).
 * Light-ray side (operators 1, 3): x0 = |y| n_theta, z1 ~ yh + n, z3 ~ yh - n (v3).  The j' = 1 celestial structure is the j' = 0 one times
   the locally constant sign of eps(z0, z1, z3) ~ sin(theta - phi_m) (the "j' = 0 x locally constant sign" of l9_points / p3), i.e. weight
   sign(s)|s|^{-delta}: only terms ODD in s = sin(theta - phi_m) survive the theta-average, with moments
       <c^{2i} |s|^{2j+1-delta}> / <|s|^{1-delta}> = (1/2)_i (1 - delta/2)_j / (3/2 - delta/2)_{i+j}
   (the common nu-only factor <|s|^{1-delta}> is dropped).  Odd terms carry one t_k = eps_k . t_hat, t_hat = eps(yh, mh, .), and
   t_k = det[y; z0; eps_k]/(-z0.y) = (w_k . y)/t  (w_k = eta.cross(z0, eps_k)); even t-pairs reduce as in v3.
   Without spin (tau13 -> 1) the j' = 1 side vanishes identically (no scalar odd functional) -- checked.
 * f = -1 carries no eta: the kernel is (|x||y|)^{a-2} K0 (one power of r more than C^_2's eta-kernel), so the transform uses
   |x|^{b-3} with b = a + 1 = 8 and NO gradient.  b only changes the nu-only constant (b = 7 is run as a check).
 * Equivariance (H7 §18.3 argument without the gradient): with weight-0 (preserving) helicity tensors, each side's transform is an
   SO(2,1)-equivariant function of (p, z0) of z0-homogeneity -lam, hence c_A (-z0.p)^-lam; then the z0 integral gives 2 pi c_A c_B
   P_{-lam}(cosh beta) = 2 pi c_A c_B P_l: the heavy image is l-ISOTYPIC (no P_{l-1}, P_{l+1}) for ANY (2,4) and (1,3) structure,
   parity-odd ones included.  Under the reflection fixing (p, z0) both o (pseudo-scalar) and the j' = 1 sign flip, so
   c_A(-h2) = -c_A(h2), c_B(-h1) = -c_B(h1): the (h1, h2) = (+,-) [++] and (+,+) [+-] images are EXACTLY opposite (sigma_h).
Checks run here: fit of L(beta) to {P_{l-1}, P_l, P_{l+1}} (residual, |h_-+|/|h_0|), ++/+- ratio, j' = 0 / no-spin controls, the
nu-dependence and real-axis zeros of the forward normalisation c3(nu) := h_0^{++}, and b = 7 vs 8.
Output: cache/q1_c3_regge_kernel.txt, cache/q1_c3_table.pkl
"""
import os, sys, io, time, functools, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
sys.dont_write_bytecode = True
import numpy as np
import sympy as sp
import v3_regge_kernel as V

ETA = V.ETA
mdot = V.mdot
_X, _s, _ex, _fx = V._X, V._s, V._ex, V._fx
_a2, _a4, _e24, _dl = V._a2, V._a4, V._e24, V._dl
_Y, _t, _U1, _U3, _A1, _A3, _E13 = V._Y, V._t, V._U1, V._U3, V._A1, V._A3, V._E13
_c, _sn = V._c, V._sn
_WE, _WF, _W1, _W3 = sp.symbols('WE WF W1 W3')          # w_e.x, w_f.x, w_1.y, w_3.y (linear markers)


def wvec(z, e):
    """w with mdot(w, v) = det[v; z; e] for all v (z: (...,3), e: (3,) or (...,3))."""
    e = np.broadcast_to(e, z.shape)
    return ETA*np.cross(z, e)


def _split_w(expr, wsyms):
    """expr linear in at most one of wsyms -> list of (wkey or None, coefficient expr)."""
    expr = sp.expand(expr)
    out = []
    rest = expr
    for w in wsyms:
        cw = sp.expand(expr.coeff(w, 1))
        if cw != 0:
            out.append((str(w), cw))
        rest = rest.subs(w, 0)
    rest = sp.expand(rest)
    if rest != 0:
        out.append((None, rest))
    return out


@functools.lru_cache(maxsize=None)
def A3_monomials(kind='C3'):
    """A(x, z0)|x|^{b-3}/|x|^{b-3} for C^_3 (q o, kind 'C3') or C^_2 (q^2, kind 'C2', control) as a list of
    (coef(a2, a4, e24, delta), k(delta), m(delta), i, j, wkey): coef * X^k s^m ex^i fx^j (w_key . x)."""
    q = 2*_a2*_a4 - _s**2*_X**-2*_e24 - 2*_s**2*_X**-4*_ex*_fx
    if kind == 'C2':
        expr = 16*_X**(1 - _dl + 4)*_s**(-(1 - _dl) - 4)*sp.expand(q**2)
    else:
        o = (_a2*_WF - _a4*_WE)/_X                                     # (z2'.z0) eps3(z4', xh, z0) - (z4'.z0) eps3(z2', xh, z0)
        expr = 16*_X**(1 - _dl + 4)*_s**(-(1 - _dl) - 4)*sp.expand(q*o)
    out = []
    for wkey, ex in _split_w(expr, (_WE, _WF)):
        for term in sp.Add.make_args(sp.expand(ex)):
            pw = term.as_powers_dict()
            k = sp.simplify(pw.get(_X, 0)); m = sp.simplify(pw.get(_s, 0)); i = int(pw.get(_ex, 0)); j = int(pw.get(_fx, 0))
            coef = sp.simplify(term/(_X**k*_s**m*_ex**i*_fx**j))
            out.append((sp.lambdify((_a2, _a4, _e24, _dl), coef, 'numpy'), sp.lambdify(_dl, k, 'numpy'),
                        sp.lambdify(_dl, m, 'numpy'), i, j, wkey))
    return out


@functools.lru_cache(maxsize=None)
def B_monomials_j(jp=1, spin=True):
    """light-ray side with the j' celestial structure (jp = 0: v3's B_monomials re-derived here as a control; jp = 1: odd).
    Returns (coef(A1, A3, E13, delta), k(delta), m(delta), i, j, wkey): coef * Y^k t^m U1^i U3^j (w_key . y)."""
    if jp == 0:
        momv = lambda ic, js: (sp.rf(sp.Rational(1, 2), ic//2)*sp.rf(sp.Rational(1, 2) - _dl/2, js//2)/sp.rf(1 - _dl/2, ic//2 + js//2)
                               if (ic % 2 == 0 and js % 2 == 0) else 0)
    else:
        momv = lambda ic, js: (sp.rf(sp.Rational(1, 2), ic//2)*sp.rf(1 - _dl/2, js//2)/sp.rf(sp.Rational(3, 2) - _dl/2, ic//2 + js//2)
                               if (ic % 2 == 0 and js % 2 == 1) else 0)
    u1, u3, m1, m3, t1, t3 = sp.symbols('u1 u3 m1 m3 t1 t3')
    if not spin:
        P = sp.Integer(1)
    else:
        P = sp.expand((u1 + _c*m1 + _sn*t1)**2*(u3 - _c*m3 - _sn*t3)**2)
    W = 0
    MS, MV = {}, {}
    for term in sp.Add.make_args(P):
        pw = term.as_powers_dict(); ic = int(pw.get(_c, 0)); js = int(pw.get(_sn, 0))
        mv = momv(ic, js)
        if mv == 0:
            continue
        W += term/(_c**ic*_sn**js)*MS.setdefault((ic, js), sp.Symbol('MOM%d_%d' % (ic, js)))
        MV[(ic, js)] = mv
    W = sp.expand(W)
    pi11, pi33, pi13 = u1**2, u3**2, _E13 + u1*u3
    T1, T3 = _W1/_t, _W3/_t                       # t_k = (w_k . y)/t
    W2 = 0
    for term in sp.Add.make_args(W):
        pw = term.as_powers_dict(); a = int(pw.get(t1, 0)); b = int(pw.get(t3, 0))
        rest = term/(t1**a*t3**b)
        rep = {(0, 0): 1, (2, 0): pi11 - m1**2, (0, 2): pi33 - m3**2, (1, 1): pi13 - m1*m3, (2, 2): (pi13 - m1*m3)**2,
               (1, 0): T1, (0, 1): T3, (2, 1): (pi11 - m1**2)*T3, (1, 2): (pi33 - m3**2)*T1}[(a, b)]
        W2 += rest*rep
    W = sp.expand(W2.subs({u1: _U1/_Y, u3: _U3/_Y, m1: _A1*_Y/_t - _U1/_Y, m3: _A3*_Y/_t - _U3/_Y}))
    out = []
    for wkey, ex in _split_w(W, (_W1, _W3)):
        Wp = sp.Poly(sp.expand(ex*_Y**8*_t**8), _Y, _t, _U1, _U3)
        for (kY, kt, i, j), coef in Wp.terms():
            k = kY - 8 + _dl; m = kt - 8 - _dl
            coef = sp.simplify(coef.subs({MS[kk]: MV[kk] for kk in MS}))
            out.append((sp.lambdify((_A1, _A3, _E13, _dl), coef, ['numpy', 'scipy']), sp.lambdify(_dl, k, 'numpy'),
                        sp.lambdify(_dl, m, 'numpy'), i, j, wkey))
    return out


def side_transform_w(p, z, e, f, monos, b, delta, coef_args, wvecs, lshift=0.0):
    """sum of monomials coef X^k s^m (e.x)^i (f.x)^j (w.x)^{0 or 1} |x|^{b-3}: Laplace transform at p, vectorised over z (no gradient)."""
    groups = {}
    for (cf, kf, mf, i, j, wk) in monos:
        key = (complex(kf(delta)), complex(mf(delta)), wk)
        groups.setdefault(key, []).append((cf, i, j))
    tot = 0
    for (k, m, wk), lst in groups.items():
        c = (b - 3) + k; lam = -m; kap = (2*lam - c - 3)/2
        Kv = np.exp(V.logKc(c, lam) - lshift)
        w = np.zeros(3) if wk is None else wvecs[wk]
        Tm = V.taylor_SR(p, z, e, f, w, lam, kap)
        kk = 0 if wk is None else 1
        acc = 0
        for (cf, i, j) in lst:
            acc = acc + cf(*coef_args, delta)*Tm[..., i, j, kk]
        tot = tot + acc*Kv
    return tot


def heavy_image3(nu, beta, h1, h2, kindA='C3', jp=1, spin=True, nphi=256, b=8):
    """L^{(h1,h2)}(beta) = int Dz0 Ahat_b(p, z0) Bhat_b(pb, z0) (no eta: f = -1).  (h1,h2) = (+,-): ++; (+,+): +-."""
    delta = 0.5 + 1j*nu
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(beta)
    E = {(+1, 'p'): ep_p, (-1, 'p'): em_p, (+1, 'pb'): ep_pb, (-1, 'pb'): em_pb}
    e2, e4 = E[(h2, 'p')], E[(-h2, 'p')]
    e1, e3 = E[(h1, 'pb')], E[(-h1, 'pb')]
    phi = 2*np.pi*np.arange(nphi)/nphi
    z = np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    E2 = np.broadcast_to(e2, z.shape); E4 = np.broadcast_to(e4, z.shape)
    E1 = np.broadcast_to(e1, z.shape); E3 = np.broadcast_to(e3, z.shape)
    cA = (mdot(z, E2), mdot(z, E4), np.full(nphi, mdot(e2, e4)))
    cB = (mdot(z, E1), mdot(z, E3), np.full(nphi, mdot(e1, e3))) if spin else (0*phi, 0*phi, 0*phi)
    wA = {'WE': wvec(z, e2), 'WF': wvec(z, e4)}
    wB = {'W1': wvec(z, e1), 'W3': wvec(z, e3)}
    lsA = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real; lsB = V.logKc(b - 3 + delta, delta).real
    Ah = side_transform_w(P, z, E2, E4, A3_monomials(kindA), b, delta, cA, wA, lsA)
    Bh = side_transform_w(PB, z, E1, E3, B_monomials_j(jp, spin), b, delta, cB, wB, lsB)
    return (Ah*Bh).mean()*2*np.pi


def image_fit(nu, h1, h2, betas, **kw):
    v = [heavy_image3(nu, bb, h1, h2, **kw) for bb in betas]
    (hm, h0, hp), res = V.fit3(nu, betas, v)
    return dict(hm=hm, h0=h0, hp=hp, res=res, vals=np.array(v), scale=np.max(np.abs(v)))


def side_consts(nu, b=8, nz=5):
    """c_A^{(3)}(h2 = +, -), c_B^{(j'=1)}(h1 = +, -), scalar c_A, c_B, and the relative spread over nz celestial points."""
    delta = 0.5 + 1j*nu
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(0.7)
    phi = 2*np.pi*(np.arange(nz) + 0.37)/nz
    z = np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    lsA = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real; lsB = V.logKc(b - 3 + delta, delta).real
    sA = (-mdot(z, P))**(1 - delta); sB = (-mdot(z, PB))**delta
    out = []
    spA = spB = 0.0
    for h in (1, -1):
        e2, e4 = (ep_p, em_p) if h == 1 else (em_p, ep_p)
        E2 = np.broadcast_to(e2, z.shape); E4 = np.broadcast_to(e4, z.shape)
        cA = (mdot(z, E2), mdot(z, E4), np.full(nz, mdot(e2, e4)))
        A = side_transform_w(P, z, E2, E4, A3_monomials('C3'), b, delta, cA, {'WE': wvec(z, e2), 'WF': wvec(z, e4)}, lsA)*sA
        out.append(A.mean()); spA = max(spA, np.max(np.abs(A - A.mean()))/abs(A.mean()))
    for h in (1, -1):
        e1, e3 = (ep_pb, em_pb) if h == 1 else (em_pb, ep_pb)
        E1 = np.broadcast_to(e1, z.shape); E3 = np.broadcast_to(e3, z.shape)
        cB = (mdot(z, E1), mdot(z, E3), np.full(nz, mdot(e1, e3)))
        B = side_transform_w(PB, z, E1, E3, B_monomials_j(1, True), b, delta, cB, {'W1': wvec(z, e1), 'W3': wvec(z, e3)}, lsB)*sB
        out.append(B.mean()); spB = max(spB, np.max(np.abs(B - B.mean()))/abs(B.mean()))
    zero = (0*phi, 0*phi, 0*phi)
    Asc = side_transform_w(P, z, P, P, [(cf, kf, mf, i, j, None) for (cf, kf, mf, i, j) in V.A_monomials(False)], b, delta, zero, {},
                           lsA)*sA
    Bsc = side_transform_w(PB, z, PB, PB, B_monomials_j(0, False), b, delta, zero, {}, lsB)*sB
    out += [Asc.mean(), Bsc.mean(), spA, spB]
    return out


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    betas = np.linspace(0.15, 3.0, 14)
    pr("# q1: heavy image of the leading-Regge kernel of Psi_{3,nu}[o_H q] (j' = 1, f = -1), split representation (v3 engine)")
    pr("## controls")
    # (c1) the j' = 0, q^2, f-stripped image through THIS code (no gradient) must be l-isotypic (v3 §18.3 item 1)
    for nu in (1.0, 6.0):
        r = image_fit(nu, 1, -1, betas, kindA='C2', jp=0, b=7)
        pr(f"  (c1) C2-structure, j'=0, no eta, nu={nu}: fit resid {r['res']:.1e}, |h-|/|h0| {abs(r['hm'])/abs(r['h0']):.1e}, "
           f"|h+|/|h0| {abs(r['hp'])/abs(r['h0']):.1e}")
    # (c2) j' = 1 without spin vanishes identically; parity bookkeeping: under the total flip (h1, h2) -> (-h1, -h2) the o q / j' = 1 image
    # is even, while the parity-mismatched combinations (q^2 dual with j' = 1, o q dual with j' = 0) are odd, so they cancel on every
    # parity-definite block (they are not zero per helicity configuration)
    for nu in (2.0,):
        bb = 0.5
        L = lambda h1, h2, k, j, sp_=True: heavy_image3(nu, bb, h1, h2, kindA=k, jp=j, spin=sp_)
        pr(f"  (c2) nu={nu}, beta={bb}: j'=1 without spin |L| = {abs(L(1, -1, 'C3', 1, False)):.1e};  "
           f"L(-h1,-h2)/L(h1,h2): o q & j'=1 {L(-1, 1, 'C3', 1)/L(1, -1, 'C3', 1):+.6f};  q^2 & j'=1 {L(-1, 1, 'C2', 1)/L(1, -1, 'C2', 1):+.6f};"
           f"  o q & j'=0 {L(-1, 1, 'C3', 0)/L(1, -1, 'C3', 0):+.6f}")
    pr("## C^_3 image: fit to h- P_{l-1} + h0 P_l + h+ P_{l+1}, both helicity classes, b = 8 (and b = 7)")
    rows = []
    nus = [0.05, 0.3, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 5.98, 6.44, 7.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0]
    for nu in nus:
        betas = np.linspace(0.12, min(3.0, max(0.6, 30/nu)), 14)   # beta window ~ 1/nu at large nu (as v3_heavy_impact.sigma_T)
        rp = image_fit(nu, 1, -1, betas)       # ++
        rm = image_fit(nu, 1, 1, betas)        # +-
        ratio = rm['h0']/rp['h0']
        diff = np.max(np.abs(rm['vals'] + rp['vals']))/rp['scale']
        r7 = image_fit(nu, 1, -1, betas, b=7) if nu in (1.0, 6.0, 20.0) else None
        rows.append(dict(nu=nu, h0pp=rp['h0'], h0pm=rm['h0'], res=max(rp['res'], rm['res']),
                         side=max(abs(rp['hm']) + abs(rp['hp']), 0)/abs(rp['h0']), ratio=ratio, diff=diff))
        pr(f"  nu={nu:7.2f}: resid {max(rp['res'], rm['res']):.1e}; (|h-|+|h+|)/|h0| = {rows[-1]['side']:.1e}; "
           f"h0(+-)/h0(++) = {ratio.real:+.12f}{ratio.imag:+.1e}i; max|L(+-) + L(++)|/max|L| = {diff:.1e}; "
           f"c3 := h0(++) = {rp['h0']:.6e}" + (f";  b=7: resid {r7['res']:.1e}, side {(abs(r7['hm']) + abs(r7['hp']))/abs(r7['h0']):.1e}"
                                              if r7 else ""))
    # side constants (equivariance): Ahat(p, z0) = c_A (-z0.p)^{-(1-delta)}, Bhat(pb, z0) = c_B (-z0.pb)^{-delta}; image = 2 pi c_A c_B P_l.
    # Zeros of the forward normalisation on the real axis = zeros of c_A or c_B; their nu-dependence relative to the SCALAR sides
    # (A = (-z0.x)^{-(1-delta)}, j' = 0 B = (-z0.y)^{-delta}; same lshift) is a ratio of Gamma-shifted rationals (cf. v4_hclosed).
    pr("## side constants c_A (o q dual), c_B (j' = 1 light-ray side), relative to the scalar sides; equivariance spread over z0")
    tab = []
    grid = np.concatenate([np.linspace(0.02, 2, 100), np.linspace(2.02, 20, 90), np.linspace(20.5, 300, 60)])
    for nu in grid:
        tab.append(side_consts(nu))
    tab = np.array(tab)              # columns: cA3(+), cA3(-), cB1(+), cB1(-), cAsc, cBsc, spreadA, spreadB
    rA = tab[:, 0]/tab[:, 4]; rB = tab[:, 2]/tab[:, 5]
    pr(f"  equivariance: max spread of c over 5 celestial points: A {np.max(np.abs(tab[:, 6])):.1e}, B {np.max(np.abs(tab[:, 7])):.1e}")
    pr(f"  helicity: max |c_A(-h)/c_A(h) + 1| = {np.max(np.abs(tab[:, 1]/tab[:, 0] + 1)):.1e};  max |c_B(-h)/c_B(h) + 1| = "
       f"{np.max(np.abs(tab[:, 3]/tab[:, 2] + 1)):.1e}")
    for lab, r in (('c_A/c_A^sc', rA), ('c_B/c_B^sc', rB)):
        ph = np.unwrap(np.angle(r)); i = np.argmin(np.abs(r))
        pr(f"  {lab}: min |.| = {np.abs(r).min():.3e} at nu = {grid[i]:.2f}; |.| at nu = 0.02 / 1 / 6 / 20 / 300: " +
           " / ".join(f"{abs(r[np.argmin(np.abs(grid - v))]):.3e}" for v in (0.02, 1, 6, 20, 300)) +
           f"; winding of arg over (0.02, 300): {(ph[-1] - ph[0])/np.pi:+.3f} pi; large-nu power fit: "
           f"{np.polyfit(np.log(grid[-20:]), np.log(np.abs(r[-20:])), 1)[0]:+.4f}")
    vals = 2*np.pi*tab[:, 0]*tab[:, 2]
    pickle.dump(dict(rows=rows, grid=grid, c3=vals, tab=tab), open(os.path.join(HERE, 'cache', 'q1_c3_table.pkl'), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'q1_c3_regge_kernel.txt'), 'w').write(buf.getvalue())
