"""
fw_frame (FOLD_WINDOW_LARGEJ.md §1): numerically stable spin-J index contraction for the graviton crystal integrands at large J.

Diagnosis [NE]: the spin-J contraction <T, S> = sum_n T_n S_{-n}/kappa_n (l5_spin: Fourier coefficients of the samples on the FIXED real
null circle w(phi) = (1 - cos, 1 + cos, sin), kappa_n = (-1)^{J+n} 2^-J C(2J, J+n)) loses ~ J log10(1.8) digits: the spin-J factor is
~ (l.z)^J with a generic complex covector l, whose samples on the fixed circle vary over many orders of magnitude, and 1/kappa_{+-J} = 2^J
amplifies the FFT round-off.  Measured: rotating the sampling phase by 0.37 (an exact identity) changes the order-0 M0 at (1048, 48),
x = 0.3, by 2.5e-4 (relative) -- the noise floor behind the J = 48 helicity-basis defect (sv3 7e-5, imag 0.057; SMOOTH_BOUND §3) -- and
generic complex null vectors give O(1) errors from l ~ 28 on.

Fix (exact identity, in-process): sample the circle in a complex Lorentz frame adapted to the spin-J covectors,
    w(phi) = e0 + cos(phi) e1 + sin(phi) e2,   e0.e0 = -1, e1.e1 = e2.e2 = 1, mutually orthogonal (metric of l5_spin.ld),
so that w.w' = cos(phi - phi') - 1 exactly as for the standard circle and the kappa_n contraction identity holds unchanged.
e0 is chosen along l_F (the z5-covector of V_{5,12} at e, the dominant spin-J factor of <T1 T2 O>), e1 in the plane of l_F and l_G
(V_{5,34}): the <T1 T2 O> side then has (up to the |n| <= 4 H-factors) a single Fourier mode, and the <T3 T4 O> side the Laplace form
(A + B cos phi)^J (no exponential cancellation in the physical region).  The frame is taken from the constant TPS term (l9/K2 engine) or
from the median of the batch (n1_pw point-wise backend); it only affects round-off, never the value.
Spin-2 circles are unchanged.  Usage: import fw_frame; fw_frame.install_tps() / install_pw(PW) before computing points.

UPDATE (the frame idea above did NOT help -- kept only as an option, FW_FRAME=1): the measured noise comes from the EDGE modes |n| ~ J.
Both sides' true Fourier coefficients there are far below the FFT round-off floor eps max|samples|, and two round-off coefficients
multiplied and divided by kappa_{+-J} = 2^-J give a relative error ~ eps^2 4^J sqrt(J) (J = 48: ~5e-3, J = 64: ~1e5).  Cure
(install_*: default): every Fourier coefficient of a spin-J (J >= 3) sample array that is below FW_THR x eps x max|samples| (per TPS
coefficient / per point of the batch) is set to 0 before the contraction.  The dropped true terms are ~ |F_n G_-n|/kappa_n with both
coefficients at the floor, i.e. negligible wherever the spectra decay (checked: J = 48, |F_n| ~ 3e-5 at n = 24, 4e-12 at 32).
"""
import os
import numpy as np
EPS = 2.220446049250313e-16
THR = float(os.environ.get('FW_THR', '1e3'))
USE_FRAME = os.environ.get('FW_FRAME', '0') == '1'

_FR = {}          # id(S module) -> current frame (3x3: rows e0, e1, e2) or None
STATS = []


def ld(a, b):
    return -0.5*(a[0]*b[1] + a[1]*b[0]) + a[2]*b[2]


def covec(c):
    """linear form z -> sum_m c_m z_m  as a vector l with ld(l, z) = sum c_m z_m."""
    return np.array([-2*c[1], -2*c[0], c[2]], complex)


def frame(lF, lG):
    e0 = lF/np.sqrt(-ld(lF, lF) + 0j)
    g = lG + ld(lG, e0)*e0
    e1 = g/np.sqrt(ld(g, g) + 0j)
    best = None
    for r in np.eye(3):
        r = r.astype(complex)
        rp = r + ld(r, e0)*e0 - ld(r, e1)*e1
        n = ld(rp, rp)
        if best is None or abs(n) > abs(best[1]):
            best = (rp, n)
    e2 = best[0]/np.sqrt(best[1] + 0j)
    F = np.array([e0, e1, e2])
    G = np.array([[ld(a, b) for b in F] for a in F])
    assert np.abs(G - np.diag([-1, 1, 1])).max() < 1e-9, G
    return F


def _circle(F, l):
    n = 2*l + 1; ph = 2*np.pi*np.arange(n)/n
    return F[0][None, :] + np.cos(ph)[:, None]*F[1][None, :] + np.sin(ph)[:, None]*F[2][None, :]


def patch_S(S):
    if getattr(S, '_fw_frame', False):
        return
    orig = S.null_circle
    ofour = S.fourier

    def fourier(vals, axis, l):
        out = ofour(vals, axis, l)
        if l <= 2 or THR <= 0:
            return out
        a = np.abs(np.asarray(vals))
        M = a.reshape(-1, a.shape[-1]).max(axis=0)            # per TPS coefficient (or per point of a batch)
        fl = THR*EPS*np.sqrt(2*l + 1)*M
        return np.where(np.abs(out) < fl, 0.0, out)
    S.fourier = fourier

    def nc(l):
        F = _FR.get(id(S))
        if F is None or l <= 2 or not USE_FRAME:
            return orig(l)
        return _circle(F, l)
    S.null_circle = nc
    S._fw_frame = True


def set_frame(S, A6, zv, rep):
    """frame from the geometry of zv (crystal variables); rep: array (N,) -> representative complex number."""
    P = A6.geom(zv, 1.0)
    x1, x2, x3, x4, e = (P[k] for k in ('x1', 'x2', 'x3', 'x4', 'e'))
    cF = [rep(S.Vfun(e, S.BASIS[m], x1, x2)) for m in range(3)]
    cG = [rep(S.Vfun(e, S.BASIS[m], x3, x4)) for m in range(3)]
    try:
        _FR[id(S)] = frame(covec(cF), covec(cG))
    except Exception as ex:          # degenerate (null l_F): fall back to the standard circle
        _FR[id(S)] = None
        STATS.append(('fallback', str(ex)[:80]))


def wrap_build(fn, S, A6, rep, zarg=0):
    def b(*a, **k):
        set_frame(S, A6, a[zarg], rep)
        try:
            return fn(*a, **k)
        finally:
            _FR[id(S)] = None
    b._fw = True
    return b


def rep_tps(v):
    return complex(np.asarray(v).ravel()[0])


def rep_pw(v):
    v = np.asarray(v).ravel()
    return complex(np.median(v.real) + 1j*np.median(v.imag))


def install_tps(P9=None):
    """TPS (Laplace) engine: l9_points.build_spin (kinds C2/C3H/C3V and, if already wrapped, EE)."""
    import l5_spin as S, l6_spin_action as A6
    import l9_points as P
    P9 = P9 or P
    patch_S(S)
    if USE_FRAME and not getattr(P9.build_spin, '_fw', False):
        P9.build_spin = wrap_build(P9.build_spin, S, A6, rep_tps)
    return P9


def install_pw(PW):
    """point-wise backend n1_pw: PW.P9.build_spin (wraps whatever is installed, incl. z2's EE dispatcher)."""
    patch_S(PW.S)
    if USE_FRAME and not getattr(PW.P9.build_spin, '_fw', False):
        PW.P9.build_spin = wrap_build(PW.P9.build_spin, PW.S, PW.A6, rep_pw)


# ------------------------------------------------------------------------------------------------ scalar Q_J by recurrence
# l1_crystal.build evaluates Q_J = sum_p c_p (Z.W)^p (Z^2 W^2)^{(J-p)/2} with the POWER-basis Legendre coefficients c_p (|c_p| ~ 2^J,
# alternating): at J = 100 and |cos| < 1 this cancels by ~10^29 (the scalar on the direct contour had 10 % MC errors at (400, 100) that
# do not fall with N).  Replacement (even J; exact identity): Q_J = (Z^2W^2)^{J/2} P_J(t), t = Z.W (Z^2W^2)^{-1/2}, P_J by the stable
# three-term recurrence on the TPS / point arrays (P_J even: the sqrt branch is irrelevant).  Composes with z2's Q_J overflow rescaling
# QSC (Z.W -> Z.W/s, Z^2W^2 -> Z^2W^2/s^2).
QL_OLD = """        Q = T.const(0.0)
        for (p, c) in legendre_poly_coeffs(J):
            term = T.const(c)
            for _ in range(p):
                term = T.mul(term, ZW)
            for _ in range((J - p)//2):
                term = T.mul(term, ZZWW)
            Q = Q + term
"""
QL_NEW = """        if J % 2 == 0:
            tq = T.mul(ZW, T.powr(ZZWW, -0.5)); p0 = T.const(1.0); p1 = tq
            for kq in range(1, J):
                p0, p1 = p1, ((2*kq + 1)*T.mul(tq, p1) - kq*p0)/(kq + 1)
            Q = p1 if J > 1 else tq
            for _ in range(J//2):
                Q = T.mul(Q, ZZWW)
        else:
            Q = T.const(0.0)
            for (p, c) in legendre_poly_coeffs(J):
                term = T.const(c)
                for _ in range(p):
                    term = T.mul(term, ZW)
                for _ in range((J - p)//2):
                    term = T.mul(term, ZZWW)
                Q = Q + term
"""


def patch_qrec_pw(PW, here):
    """exec a copy of l1_crystal into PW.Cr with the QSC rescaling (z2) AND the recurrence; z2.patch_qscale then returns early."""
    import sys as _s
    if PW.Cr.__dict__.get('_fw_qrec'):
        return
    src = open(os.path.join(here, 'l1_crystal.py')).read()
    import z2_points as Z2
    assert src.count(Z2.QSC_SRC_OLD) == 1 and src.count(QL_OLD) == 1
    src = src.replace(Z2.QSC_SRC_OLD, Z2.QSC_SRC_NEW).replace(QL_OLD, QL_NEW)
    saved = _s.modules.get('l1_tps')
    _s.modules['l1_tps'] = PW
    try:
        PW.Cr.__dict__['QSC'] = [1.0]
        exec(compile(src, 'l1_crystal_pw_qsc_qrec', 'exec'), PW.Cr.__dict__)
    finally:
        _s.modules['l1_tps'] = saved
    PW.Cr.__dict__['_fw_qrec'] = True


def patch_qrec_tps(Cr, here):
    """the same recurrence in the TPS l1_crystal module (K2 engine), in place."""
    if Cr.__dict__.get('_fw_qrec'):
        return
    src = open(os.path.join(here, 'l1_crystal.py')).read()
    assert src.count(QL_OLD) == 1
    exec(compile(src.replace(QL_OLD, QL_NEW), 'l1_crystal_qrec', 'exec'), Cr.__dict__)
    Cr.__dict__['_fw_qrec'] = True


# ------------------------------------------------------------------------------------------------ conserved <TTO> null space at J >~ 100
# l4_tto_even.conserved_even samples random points with null z = (u, y^2/u, y), u, y ~ N(0,1): |z| is unbounded and the degree-J structure
# values overflow (J = 150: SVD 'did not converge' on inf/nan).  The structures are homogeneous in each z (degree 2, 2, J), so normalising
# every z to unit max-norm multiplies a point's structure values and conservation rows by a common constant: the null space is unchanged.
def patch_conserved(S=None):
    if S is None:
        import l5_spin as S
    L = S._l4r()
    if getattr(L, '_fw_norm', False):
        return
    orig = L.rand_pt

    def rand_pt(rng):
        x1, z1, x2, z2, x5, z5 = orig(rng)
        nz = lambda z: z/np.abs(z).max()
        return (x1, nz(z1), x2, nz(z2), x5, nz(z5))
    L.rand_pt = rand_pt
    L._fw_norm = True
    S.conserved.cache_clear()


def patch_conserved_clean(S=None):
    """conserved_even with denormal/underflow clean-up: after l4's row scaling, entries below 1e-200 (relative) are set to 0 (LAPACK SVD
    fails to converge on denormals at J >~ 150); otherwise a verbatim copy of l4_tto_even.conserved_even."""
    if S is None:
        import l5_spin as S
    L = S._l4r()
    patch_conserved(S)
    if getattr(L, '_fw_clean', False):
        return

    def conserved_even(Delta, J, npts=None, seed=3):
        rng = np.random.default_rng(seed)
        mons = L.monomials(J); K = len(mons)
        npts = npts or 3*K
        F_, C1_, C2_, S_ = [], [], [], []
        ntry = 0
        while len(F_) < npts and ntry < 20*npts:       # points whose values over/underflow (J >~ 150) are skipped
            ntry += 1
            p = L.rand_pt(rng)
            with np.errstate(all='ignore'):
                f = np.array(L.structures(*p, Delta, J), complex)
                c1 = np.array(L.conservation_row(p, Delta, J, 1), complex)
                c2 = np.array(L.conservation_row(p, Delta, J, 2), complex)
                sv_ = np.array(L.structures(p[2], p[3], p[0], p[1], p[4], p[5], Delta, J), complex) - f
            if all(np.all(np.isfinite(a)) for a in (f, c1, c2, sv_)) and np.abs(f).max() > 1e-250:
                F_.append(f); C1_.append(c1); C2_.append(c2); S_.append(sv_)
        Fv = np.array(F_); Cv = np.array(C1_ + C2_); Sv = np.array(S_)
        rn = lambda A_: A_/np.maximum(np.abs(A_).max(axis=1, keepdims=True), 1e-300)
        sc = np.maximum(np.abs(Fv).max(axis=1, keepdims=True), 1e-300)
        Fv, Cv, Sv = Fv/sc, rn(Cv), Sv/sc
        cl = lambda A_: np.where(np.abs(A_) < 1e-200, 0.0, A_)
        Fv, Cv, Sv = cl(Fv), cl(Cv), cl(Sv)
        U, sF, Vh = np.linalg.svd(Fv); rF = int(np.sum(sF > 1e-9*sF[0]))
        B = Vh[:rF].conj().T
        A = np.vstack([Cv @ B, Sv @ B])
        sA = np.linalg.svd(A, compute_uv=False)
        n0 = int(np.sum(sA < 1e-8*sA[0]))
        _, _, VhA = np.linalg.svd(A)
        sol = B @ VhA[len(VhA) - n0:].conj().T if n0 else np.zeros((K, 0))
        return dict(K=K, rank3d=rF, n_conserved_sym=n0, sv=sA[-4:]/sA[0], sol=sol, mons=mons)
    L.conserved_even = conserved_even
    L._fw_clean = True
    S.conserved.cache_clear()
