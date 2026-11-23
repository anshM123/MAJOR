"""
p2con/pc_ksh.py (P2_CONTRACT.md §1.7): the EE kernel's OWN O(r) term (Theorem W's K_sh), exactly, paired with the leading heavy block.

(a) Exact O(r) structure of the EE kernel [D, checked numerically against ee_quad's exact integrand]:
    in the H7 frame (x2 = -xh, x4 = xh, x1 = -r yh, x3 = r yh, x0 = r n_th) the scalar factor S_nu (m20 m40, Nn, c13 = 4r^2, c10, c30)
    and tau13 = (w1.x10)^2 (w3.x30)^2 / r^4 have NO O(r) term (m20 m40 = (1-r^2)^2 - 4r^2 (xh.n)^2, Nn = -4 z0.xh (1-r^2) - 8r^2(...)),
    and X20 X40 = 1 + O(r^2).  The only O(r) term is in the (24) numerator:
        D = 16 (z2'.z4')^2 + 128 r (z2'.z4') [ (z2'.n)(z4'.xh) - (z2'.xh)(z4'.n) ] + O(r^2)       (z' inversion-referred).
    Hence the O(r) kernel is of split form with ONE extra vector contraction between the two sides:
        K1 = int Dz0 { 128 A0 e24 (z4'.x)  x  B0[ |y| (z2'.n_th) ]  -  128 A0 e24 (z2'.x)  x  B0[ |y| (z4'.n_th) ] },
    A0 = (-z0.xh)^{-(1-delta)} |x|^.. (the EE A side without e24^2), B0[g] = (-z0.yh)^{-delta} int dth |sin th|^{-delta} tau13 g(th).
(b) Its plane-wave image with the leading (pure-helicity) block tensor, through pc_engine (new B monomials with the extra vector v),
    de-shifted and costed exactly as in pc_run.  This is the true 'K_sh' content of P2 (index-shifted / non-ring at O(1/m^2)).
"""
import os, sys, functools, math
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import pc_engine as E

_Y, _t, _U1, _U3, _UV = sp.symbols('Y t U1 U3 UV')
_A1, _A3, _E13, _AV, _E1V, _E3V = sp.symbols('A1 A3 E13 AV E1V E3V')
_dl = sp.Symbol('delta')
_c, _sn = sp.symbols('c sn')


@functools.lru_cache(maxsize=None)
def B1_monomials():
    """B0[|y| (v.n_th)] as sum coef(A1, A3, E13, AV, E1V, E3V, delta) Y^k t^m U1^i U3^j UV^l (v3.B_monomials + one vector v):
    v.n = c (v.m) + sn (v.t); t-component products reduced by nullness / completeness (one transverse direction)."""
    momv = lambda i2, j2: sp.rf(sp.Rational(1, 2), i2)*sp.rf(sp.Rational(1, 2) - _dl/2, j2)/sp.rf(1 - _dl/2, i2 + j2)
    u1, u3, m1, m3, t1, t3, uv, mv, tv = sp.symbols('u1 u3 m1 m3 t1 t3 uv mv tv')
    P = sp.expand((u1 + _c*m1 + _sn*t1)**2*(u3 - _c*m3 - _sn*t3)**2*(_c*mv + _sn*tv))
    W = 0
    for term in sp.Add.make_args(P):
        pw = term.as_powers_dict(); ic = int(pw.get(_c, 0)); js = int(pw.get(_sn, 0))
        if ic % 2 or js % 2:
            continue
        W += term/(_c**ic*_sn**js)*momv(ic//2, js//2)
    W = sp.expand(W)
    pi11, pi33, pi13 = u1**2, u3**2, _E13 + u1*u3
    pi1v, pi3v = _E1V + u1*uv, _E3V + u3*uv
    red = {(0, 0, 0): 1, (2, 0, 0): pi11 - m1**2, (0, 2, 0): pi33 - m3**2, (1, 1, 0): pi13 - m1*m3, (2, 2, 0): (pi13 - m1*m3)**2,
           (1, 0, 1): pi1v - m1*mv, (0, 1, 1): pi3v - m3*mv, (2, 1, 1): (pi11 - m1**2)*(pi3v - m3*mv), (1, 2, 1): (pi33 - m3**2)*(pi1v - m1*mv)}
    W2 = 0
    for term in sp.Add.make_args(W):
        pw = term.as_powers_dict(); a = int(pw.get(t1, 0)); b = int(pw.get(t3, 0)); e = int(pw.get(tv, 0))
        W2 += term/(t1**a*t3**b*tv**e)*red[(a, b, e)]
    W = sp.expand(_Y*W2.subs({u1: _U1/_Y, u3: _U3/_Y, uv: _UV/_Y, m1: _A1*_Y/_t - _U1/_Y, m3: _A3*_Y/_t - _U3/_Y,
                              mv: _AV*_Y/_t - _UV/_Y}))
    Wp = sp.Poly(sp.expand(W*_Y**10*_t**10), _Y, _t, _U1, _U3, _UV)
    out = []
    for (kY, kt, i, j, l), coef in Wp.terms():
        k = kY - 10 + _dl; m = kt - 10 - _dl
        out.append((sp.lambdify((_A1, _A3, _E13, _AV, _E1V, _E3V, _dl), sp.simplify(coef), 'numpy'), sp.lambdify(_dl, k, 'numpy'),
                    sp.lambdify(_dl, m, 'numpy'), i, j, l, coef))
    return out


def B1_direct(ypt, z0, z1, z3, v, delta, nth=4000):
    """position-space B0[|y| (v.n_th)] by direct th-quadrature (check of the monomial algebra); vectors in Cartesian (-,+,+)."""
    md = E.mdot
    Yn = np.sqrt(-md(ypt, ypt)); yh = ypt/Yn
    S = -md(z0, yh)
    mh = z0/S - yh
    # t-hat: unit vector orthogonal to yh, mh (orientation irrelevant: odd terms vanish)
    tv = np.cross(np.array([-yh[0], yh[1], yh[2]]), np.array([-mh[0], mh[1], mh[2]]))
    tv = np.array([-tv[0], tv[1], tv[2]]) if False else tv
    # orthonormalise numerically w.r.t. the Minkowski metric
    for w in (yh, mh):
        tv = tv - md(tv, w)/md(w, w)*w
    tv = tv/np.sqrt(md(tv, tv))
    th = (np.arange(nth) + 0.5)*2*np.pi/nth
    c, s = np.cos(th), np.sin(th)
    n = c[:, None]*mh[None, :] + s[:, None]*tv[None, :]
    tau = md(z1[None, :], yh[None, :] + n)**2*md(z3[None, :], yh[None, :] - n)**2
    g = Yn*md(v[None, :], n)
    w = np.abs(s)**(-delta)
    # normalised by <1> = int |s|^-delta dth (as v3: moments divided by the common nu-only factor)
    return S**(-delta)*np.sum(w*tau*g)/np.sum(w)


def B1_mono_eval(ypt, z0, z1, z3, v, delta):
    md = E.mdot
    Yn = np.sqrt(-md(ypt, ypt)); t = -md(z0, ypt)
    args = (md(z1, z0), md(z3, z0), md(z1, z3), md(v, z0), md(z1, v), md(z3, v), delta)
    U1, U3, UV = md(z1, ypt), md(z3, ypt), md(v, ypt)
    tot = 0
    for (cf, kf, mf, i, j, l, _) in B1_monomials():
        tot += cf(*args)*Yn**kf(delta)*t**mf(delta)*U1**i*U3**j*UV**l
    return tot
