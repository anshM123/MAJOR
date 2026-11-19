"""v4 (H7_I14_LIGHTSIDE.md §23): EXACT closed form of the spinning magic-block coefficients h_-+(nu) of Psi_{2,nu}[q^2] (§18) and their
analytic structure (zeros/poles in the strip |Im nu| <= 1 and on the real axis).

Method (exact computer algebra on the §18.1 kernel; no fits).
  * §18.3's proof: the gradient of the plane-wave transform of each side is an equivariant covector,
        V_mu(p, z0) = (-z0.p)^-lam [c1 p_mu + c2 z0_mu/(-z0.p)],   W_mu(pb, z0) = (-z0.pb)^-lam' [c1' pb_mu + c2' z0_mu/(-z0.pb)],
    lam = 1 - delta, lam' = delta, delta = 1/2 + i nu.  With I(a) := int dphi (-z0.p)^-a (-z0.pb)^(a-1) = 2 pi P_{-a}(cosh beta) [checked
    numerically], L = -int V.W = 2 pi [c1 c1' eta P_l + c1 c2' P_{l+1} + c2 c1' P_{l-1}], eta P_l = [(l+1) P_{l+1} + l P_{l-1}]/(2l+1):
        h_+ = 2 pi [c1 c1' (l+1)/(2l+1) + c1 c2'],    h_- = 2 pi [c1 c1' l/(2l+1) + c2 c1'],      l = -1/2 + i nu.
  * c1, c2 (c1', c2') from ONE celestial point (p = (1,0,0), z0 = (1,0,1)): the transform of each monomial |x|^c (-z0.x)^-lam (e.x)^i
    (f.x)^j is K(c,lam) (e.d)^i (f.d)^j F_{c,lam} (v3_regge_kernel); derivatives done symbolically (sympy), and all K(c,lam) of one side
    differ from a reference K by integer shifts of both Gamma arguments -> rising factorials.  Hence
        h_-+ = 2 pi K_A^ref(delta) K_B^ref(delta) R_-+(nu),    R_-+ rational in nu (exact),
    with K^ref the scalar-functional constants (same group), so the spin/scalar ratios a_-+ = h^T/h^sc (§18.4) are rational, and so is
    rho := h_+/h_- (normalisation-free).  The scalar R^sc reproduces CHMRS: R_-^sc/R_+^sc = g/g~ exactly (checked symbolically).
  * Natural normalisation (the only one fixed by the scalar): H_- := g a_-, H_+ := g~ a_+ (Psi^T with the same dual/celestial prefactors as
    the scalar Psi, divided by the scalar's nu-only factor so that the scalar has kernel (g, g~)).  Zeros/poles of H_-+, H_- + H_+ = the
    zeros/poles of these explicit meromorphic functions; g carries Gamma factors with no zeros/poles in |Im nu| <= 1 apart from those listed.
Output: cache/v4_hclosed.txt, cache/v4_hclosed.pkl (sympy srepr of R_-+^T, R_-+^sc; numpy-callable via load()).
Run time ~1-3 min (one core).
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import mpmath as mp
import v3_regge_kernel as V

d = V._dl                                     # delta symbol
nu_s = sp.Symbol('nu')
I = sp.I
OUTP = os.path.join(HERE, 'cache', 'v4_hclosed.pkl')


def _monos_sym(spin, side):
    """the monomial lists of v3_regge_kernel with SYMPY coefficients (lambdify intercepted)."""
    orig = V.sp.lambdify
    V.sp.lambdify = lambda args, expr, *a, **k: expr
    try:
        f = (V.A_monomials if side == 'A' else V.B_monomials).__wrapped__
        return f(spin)
    finally:
        V.sp.lambdify = orig


def _gratio(x, n):
    """Gamma(x + n)/Gamma(x) for integer n."""
    n = int(n)
    if n >= 0:
        return sp.Mul(*[x + k for k in range(n)])
    return 1/sp.Mul(*[x + n + k for k in range(-n)])


def _md(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


_p = sp.symbols('q0:3')
_lam, _kap = sp.symbols('Lam Kap')


def _dpat(e, f, z, i, j, mu):
    """[(e.d)^i (f.d)^j d_mu (S^-lam R^kap)] / (S^-lam R^kap) at p = (1,0,0) (S = -z.p = 1, R = -p.p = 1)."""
    S = -_md(z, _p); R = -_md(_p, _p)
    F = S**(-_lam)*R**_kap
    G = sp.diff(F, _p[mu])
    for _ in range(i):
        G = sum(e[k]*sp.diff(G, _p[k]) for k in range(3))
    for _ in range(j):
        G = sum(f[k]*sp.diff(G, _p[k]) for k in range(3))
    G = G/F
    return sp.expand(sp.simplify(G.subs({_p[0]: 1, _p[1]: 0, _p[2]: 0})))


def side_coeffs(spin, side):
    """(c1, c2) / K^ref (exact, rational in delta) for side A (pair 2,4) or B (pair 1,3); K^ref = K of the scalar group."""
    z = (1, 0, 1)
    s2 = 1/sp.sqrt(2)
    em = (0, s2, -I*s2); ep = (0, s2, I*s2)                 # eps_{-1}, eps_{+1} at p = (1,0,0) (frames(0), parallel frame)
    if side == 'A':
        e, f = em, ep                                       # z2' = eps_{h2}, z4' = eps_{-h2}, h2 = -1 (helicity-blind, §18.3)
        kref, mref = 1 - d, -(1 - d)
    else:
        e, f = ep, em                                       # z1 = eps_{h1}, z3 = eps_{-h1}, h1 = +1
        kref, mref = d, -d
    cargs = (_md(z, e), _md(z, f), _md(e, f))
    monos = _monos_sym(spin, side)
    cref = 3 + kref; lref = -mref
    V3 = [0, 0, 0]
    cache = {}
    for (coef, kf, mf, i, j) in monos:
        k = sp.sympify(kf); m = sp.sympify(mf)
        c = 3 + k; lam = -m; kap = (2*lam - c - 3)/2
        if side == 'A':
            cv = coef.subs({V._a2: cargs[0], V._a4: cargs[1], V._e24: cargs[2]})
        else:
            cv = coef.subs({V._A1: cargs[0], V._A3: cargs[1], V._E13: cargs[2]})
        # K(c, lam)/K(cref, lref): 2^{mu - muref} Gamma(1 + c/2)/Gamma(1 + cref/2) Gamma((c+3)/2 - lam)/Gamma((cref+3)/2 - lref)
        n1 = sp.simplify(c/2 - cref/2); n2 = sp.simplify((c + 3)/2 - lam - (cref + 3)/2 + lref)
        dmu = sp.simplify((c + 3 - lam) - (cref + 3 - lref))
        assert n1.is_integer and n2.is_integer and dmu.is_integer, (n1, n2, dmu)
        Kr = sp.Integer(2)**dmu*_gratio(1 + cref/2, n1)*_gratio((cref + 3)/2 - lref, n2)
        for mu in range(3):
            key = (i, j, mu)
            if key not in cache:
                cache[key] = _dpat(e, f, z, i, j, mu)
            D = cache[key].subs({_lam: lam, _kap: kap})
            V3[mu] += cv*Kr*D
    V3 = [sp.factor(sp.simplify(v)) for v in V3]
    assert sp.simplify(V3[1]) == 0, V3[1]                  # no component outside span(p, z0): equivariance
    c2 = V3[2]; c1 = sp.factor(sp.simplify(-V3[0] - c2))     # p_mu = (-1,0,0), z0_mu = (-1,0,1)
    return c1, c2


def R_pm(spin):
    c1, c2 = side_coeffs(spin, 'A'); d1, d2 = side_coeffs(spin, 'B')
    l = d - 1
    Rp = sp.factor(sp.simplify(c1*d1*(l + 1)/(2*l + 1) + c1*d2))
    Rm = sp.factor(sp.simplify(c1*d1*l/(2*l + 1) + c2*d1))
    sub = {d: sp.Rational(1, 2) + I*nu_s}
    return sp.factor(sp.simplify(Rm.subs(sub))), sp.factor(sp.simplify(Rp.subs(sub)))


def g_sym(sign=+1):
    """CHMRS g(nu) (sign=+1) or g~(nu) = g(-nu) (sign=-1), sympy."""
    v = sign*nu_s
    lg = lambda w: sp.gamma((sp.Rational(13, 2) - I*w)/2)*sp.gamma((sp.Rational(13, 2) + I*w)/2)
    return (-sp.Rational(1, 2) + I*v)/(I*v)*(lg(v + I)/lg(v))**2


def build():
    t0 = time.time()
    RmT, RpT = R_pm(True)
    RmS, RpS = R_pm(False)
    pickle.dump(dict(RmT=sp.srepr(RmT), RpT=sp.srepr(RpT), RmS=sp.srepr(RmS), RpS=sp.srepr(RpS)), open(OUTP, 'wb'))
    return RmT, RpT, RmS, RpS, time.time() - t0


def load():
    """sympy expressions R_-^T, R_+^T, R_-^sc, R_+^sc (keys RmT, RpT, RmS, RpS) of nu."""
    if not os.path.exists(OUTP):
        build()
    D = pickle.load(open(OUTP, 'rb'))
    return {k: eval(v, {**sp.__dict__}) for k, v in D.items()}


def numfuncs():
    """numpy callables of nu (complex ok) for RmT, RpT, RmS, RpS."""
    return {k: sp.lambdify(nu_s, v, 'numpy') for k, v in load().items()}


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    RmT, RpT, RmS, RpS, dt = build()
    pr(f"# v4_hclosed: exact h_-+ of Psi_2[q^2] (built in {dt:.0f}s).  h_-+ = 2 pi K_A^ref K_B^ref R_-+ ;  R rational in nu.")
    # --- scalar check: R_-^sc/R_+^sc = g/g~ (CHMRS), exactly
    gq = sp.simplify(sp.gammasimp(g_sym(+1)/g_sym(-1)))
    chk = sp.simplify(sp.gammasimp(RmS/RpS/gq))
    pr(f"scalar: R_-^sc/R_+^sc / (g/g~) = {chk}    (g/g~ = {sp.factor(gq)})")
    rho = sp.factor(sp.simplify(RpT/RmT))
    am = sp.factor(sp.simplify(RmT/RmS)); ap = sp.factor(sp.simplify(RpT/RpS))
    pr(f"rho := h_+/h_- = {rho}")
    pr(f"a_- = h_-^T/h_-^sc = {am}")
    pr(f"a_+ = h_+^T/h_+^sc = {ap}")
    pr(f"R_-^T = {RmT}")
    pr(f"R_+^T = {RpT}")
    pr(f"R_-^sc = {RmS}")
    pr(f"R_+^sc = {RpS}")
    # symmetry nu -> -nu
    pr(f"rho(nu) rho(-nu) = {sp.simplify(rho*rho.subs(nu_s, -nu_s))}")
    pr(f"a_-(-nu)/a_+(nu) = {sp.simplify(am.subs(nu_s, -nu_s)/ap)}")
    # ---- zeros and poles (CHMRS normalisation H_-+ = N_com R_-+^T, N_com = g/R_-^sc = g~/R_+^sc: even, no zeros/poles for |Im nu| < 11/2)
    Ncom = sp.simplify(sp.gammasimp(g_sym(+1)/RmS))
    pr(f"N_com := g/R_-^sc = {Ncom}   (numerical check max|g~/(R_+^sc N_com) - 1| on nu = 0.3, 2.7+0.4i, 9-0.8i: {max(abs(complex((g_sym(-1)/RpS/Ncom).subs(nu_s, v).evalf()) - 1) for v in (0.3, 2.7 + 0.4j, 9 - 0.8j)):.1e})")
    def zp(expr, lab):
        num, den = sp.fraction(sp.factor(sp.together(expr)))
        rz = np.roots([complex(c) for c in sp.Poly(sp.expand(num), nu_s).all_coeffs()]) if sp.degree(num, nu_s) > 0 else np.array([])
        rp = np.roots([complex(c) for c in sp.Poly(sp.expand(den), nu_s).all_coeffs()]) if sp.degree(den, nu_s) > 0 else np.array([])
        fm = lambda r: ", ".join(f"{z.real:+.5f}{z.imag:+.5f}i" for z in sorted(r, key=lambda z: (round(z.imag, 6), z.real)))
        ins = lambda r: [z for z in r if abs(z.imag) <= 1 + 1e-12]
        pr(f"{lab}:\n    zeros {fm(rz)}\n    poles {fm(rp)}\n    in |Im nu| <= 1: zeros [{fm(ins(rz))}], poles [{fm(ins(rp))}]")
    zp(RmT, "h_- (= N_com R_-^T)")
    zp(RpT, "h_+ (= N_com R_+^T)")
    zp(RmT + RpT, "h_- + h_+ (forward normalisation)")
    zp(rho, "rho = h_+/h_-")
    zp(RmS, "g (scalar, / N_com)"); zp(RpS, "g~ (scalar, / N_com)"); zp(RmS + RpS, "g + g~ (scalar)")
    s0 = sp.nsimplify(sp.limit(-I*nu_s*(RmT - RpT)/(RmT + RpT), nu_s, 0)); s0s = sp.limit(-I*nu_s*(RmS - RpS)/(RmS + RpS), nu_s, 0)
    pr(f"small nu: u = i s0/nu + O(1): s0 = {s0} = {float(s0):.5f} (spinning), {s0s} (scalar g), 1/2 (eta P_l); large nu: s -> 23/2")
    f = {k: sp.lambdify(nu_s, v, 'numpy') for k, v in (('m', RmT), ('p', RpT))}
    nv = np.linspace(1e-3, 300, 3000001); a, b = f['m'](nv), f['p'](nv); r = np.abs(a + b)/(np.abs(a) + np.abs(b))
    pr(f"real axis nu in (0, 300]: min |h_-| nu = {np.min(np.abs(a)*nv):.3f}, min |h_+| nu = {np.min(np.abs(b)*nv):.3f} (no zeros); "
       f"|h_- + h_+|/(|h_-| + |h_+|) local minima: " + ", ".join(f"{r[(nv > lo) & (nv < hi)].min():.4f} at {nv[(nv > lo) & (nv < hi)][np.argmin(r[(nv > lo) & (nv < hi)])]:.3f}"
                                                                 for lo, hi in ((0.5, 1.0), (1.5, 2.5), (5, 8))) + f"; max |u_T| on [0.2, 300] = {np.max(np.abs((a - b)/(a + b))[nv > 0.2]):.3f}")
    open(os.path.join(HERE, 'cache', 'v4_hclosed.txt'), 'w').write(buf.getvalue())
