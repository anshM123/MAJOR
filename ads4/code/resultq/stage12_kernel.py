"""Result Q, Stage 1 (RESULT_Q_STAGE12.md §1): the CLS-D subtraction factors f^{+}_{k;2} (even j') and f^{-}_{k;2} (odd j') for
<TTTT> (J_e = 2), CLS-D arXiv:2311.04271 eqs. (subtraction_factor_evenjprime), (subtraction_factor_oddjprime), in two limits:

 (A) the CHMRS Regge frame  x1 = -y, x3 = y, x2 = x/x^2, x4 = -x/x^2  (Lorentzian, x, y future timelike, r = |x||y|,
     eta = -x.y/(|x||y|) >= 1): exact u', v' in (r, eta), and the Laurent expansion of f_k at r -> 0 (leading power and eta-degree);
 (B) the scattering crystal (CLS-D eq. scattering_crystal): u' = m^4/nu^4, v' = (m^2-nu^2)^2/nu^4, and the CLS-D crystal values of f_k
     (checked symbolically); ratio f_4/f_2, f_5/f_3 in flat variables (s = m^2, u = -nu^2).

Everything is exact (sympy) except a random-configuration numerical check of the frame cross-ratios.
Output: cache/stage12_kernel.txt.  Run time: seconds.
"""
import os, io
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


up, vp, r, eta, m, nu = sp.symbols("u' v' r eta m nu", positive=True)


def f_plus(k, Je=2):
    """CLS-D f^+_{k;J_e}(u', v'), k even."""
    e = k - 2 * Je
    if (k - (2 * Je + 2)) % 4 == 0:
        return (vp - up) / (up * vp) ** sp.Rational(e + 2, 4)
    assert (k - 2 * Je) % 4 == 0
    return ((vp - up) ** 3 - (vp - up)) / (4 * (up * vp) ** sp.Rational(e + 4, 4))


def f_minus(k, Je=2):
    """CLS-D f^-_{k;J_e}(u', v'), k odd."""
    e = k - 2 * Je
    if (k - (2 * Je + 3)) % 4 == 0:
        return -1 / (up * vp) ** sp.Rational(e + 1, 4)
    assert (k - (2 * Je + 1)) % 4 == 0
    return -((up - vp) ** 2 - 1) / (4 * (up * vp) ** sp.Rational(e + 3, 4))


def mink(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def frame_check(ntest=2000, seed=1):
    """numerical: u' = (1 - 2 r eta + r^2)^2/(16 r^2), v' = (1 + 2 r eta + r^2)^2/(16 r^2) in the CHMRS frame."""
    rng = np.random.default_rng(seed)
    worst = 0.
    for _ in range(ntest):
        def tl():
            v = rng.normal(size=3); v[0] = np.sqrt(v[1] ** 2 + v[2] ** 2) + abs(rng.normal()) + 0.05
            return v * rng.uniform(0.05, 1.5)
        x, y = tl(), tl()
        x2 = mink(x, x); y2 = mink(y, y)
        X = [None, -y, x / x2, y, -x / x2]
        d2 = lambda i, j: mink(X[i] - X[j], X[i] - X[j])
        U = d2(1, 2) * d2(3, 4) / (d2(1, 3) * d2(2, 4)); V = d2(1, 4) * d2(2, 3) / (d2(1, 3) * d2(2, 4))
        rr = np.sqrt(x2 * y2); ee = -mink(x, y) / rr
        Uf = (1 - 2 * rr * ee + rr ** 2) ** 2 / (16 * rr ** 2); Vf = (1 + 2 * rr * ee + rr ** 2) ** 2 / (16 * rr ** 2)
        worst = max(worst, abs(U / Uf - 1), abs(V / Vf - 1))
    return worst


def main():
    pr("# stage12_kernel: CLS-D subtraction factors for <TTTT> (J_e = 2) in the Regge frame and at the scattering crystal")
    pr(f"(A0) frame cross-ratios: max relative deviation of u', v' from (1 -+ 2 r eta + r^2)^2/(16 r^2) on 2000 random configs: {frame_check():.1e}")
    Uf = (1 - 2 * r * eta + r ** 2) ** 2 / (16 * r ** 2)
    Vf = (1 + 2 * r * eta + r ** 2) ** 2 / (16 * r ** 2)
    pr("    v' - u' =", sp.factor(sp.simplify(Vf - Uf)), ";  u'v' =", sp.factor(sp.simplify(Uf * Vf)))
    pr("(A) Laurent expansion at r -> 0, fixed eta (leading term; [P], sympy).  kernel-relative factor = f_k / f_{k0} with k0 = 2 (even) or 3 (odd)")
    lead = {}
    for k in range(2, 9):
        f = f_plus(k) if k % 2 == 0 else f_minus(k)
        fr = sp.simplify(f.subs({up: Uf, vp: Vf}))
        ser = sp.series(fr, r, 0, 6).removeO()
        ser = sp.expand(ser)
        # leading term
        lt = None
        for pw in range(-4, 8):
            c = sp.simplify(ser.coeff(r, pw))
            if c != 0:
                lt = (pw, sp.factor(c)); break
        lead[k] = lt
        pr(f"    k = {k} ({'f+' if k % 2 == 0 else 'f-'}): f_k = {sp.factor(fr)}")
        pr(f"           leading: r^{lt[0]} * {lt[1]}")
    for k in (4, 6, 8):
        pr(f"    f_{k}/f_2 leading: r^{lead[k][0] - lead[2][0]} * {sp.factor(lead[k][1] / lead[2][1])}")
    for k in (5, 7):
        pr(f"    f_{k}/f_3 leading: r^{lead[k][0] - lead[3][0]} * {sp.factor(lead[k][1] / lead[3][1])}")
    pr("    => k = 4: leading kernel = (k=2 kernel) x 16 eta^2 r^2: moment shifted by 2 (Pi_4), eta^3 overall -> heavy image mixes P_{l+-1}, P_{l+-3};")
    pr("       k = 5: (k=3 kernel) x 16 eta^2 r^2: Pi_5, eta^2 -> mixes P_l, P_{l+-2}.  (CLS-D's Regge-spin shifts k-2J_e+1, k-2J_e.)")
    # (B) crystal
    pr("(B) scattering crystal (CLS-D (scattering_crystal)), lightcone coords x = (u, v, y), x^2 = -u v + y^2:")
    a = sp.Symbol('a', positive=True)
    def lc(p, q): return -(p[0] * q[1] + p[1] * q[0]) / 2 + p[2] * q[2]
    X = (-sp.I * a, sp.I * a, 0); Y = (sp.I * a, -sp.I * a, 0)
    x2 = sp.simplify(lc(X, X)); y2 = sp.simplify(lc(Y, Y)); xy = sp.simplify(lc(X, Y))
    pr(f"    x^2 = {x2}, y^2 = {y2}, x.y = {xy}")
    Uc = sp.simplify((1 + 2 * xy + x2 * y2) ** 2 / (16 * x2 * y2)); Vc = sp.simplify((1 - 2 * xy + x2 * y2) ** 2 / (16 * x2 * y2))
    # x1 = -y, x2 = x/x^2: x12^2 = (x2 y2 + 2 x.y + 1)/x2 -> u' = (1 + 2 x.y + x^2 y^2)^2 / (16 x^2 y^2)
    aval = (m - sp.sqrt(m ** 2 - nu ** 2)) / nu
    Ucs = sp.simplify(sp.radsimp(Uc.subs(a, aval))); Vcs = sp.simplify(sp.radsimp(Vc.subs(a, aval)))
    pr(f"    u' = {sp.factor(Ucs)},  v' = {sp.factor(Vcs)}  (expected m^4/nu^4, (m^2-nu^2)^2/nu^4)")
    cls = {2: -nu ** -2 * (2 * m ** 2 - nu ** 2),
           3: sp.Integer(1)}
    def clsd_plus(k): return -nu ** (2 * (k - 4 + 1)) * (2 * m ** 2 - nu ** 2) / (m ** 2 * (m ** 2 - nu ** 2)) / (m ** 2 * (m ** 2 - nu ** 2)) ** sp.Rational(k - 4, 2)
    def clsd_minus(k): return -nu ** (2 * (k - 4)) * nu ** 2 / (m ** 2 * (m ** 2 - nu ** 2)) / (m ** 2 * (m ** 2 - nu ** 2)) ** sp.Rational(k - 5, 2)
    pr("    (CLS-D's displayed unsubtracted value 'f^-_{3;2}|crystal = 1' differs in sign from their general crystal formula and from the")
    pr("     definition f^-_{3;2} = -1; the general formula is used below. Overall signs of the odd functional are a convention.)")
    for k in range(2, 8):
        f = f_plus(k) if k % 2 == 0 else f_minus(k)
        val = sp.simplify(f.subs({up: m ** 4 / nu ** 4, vp: (m ** 2 - nu ** 2) ** 2 / nu ** 4}))
        ref = clsd_plus(k) if k % 2 == 0 else clsd_minus(k)
        pr(f"    k = {k}: f_k|crystal = {sp.factor(val)};  CLS-D formula: {sp.factor(sp.simplify(ref))};  difference = {sp.simplify(val - ref)}")
    s_, u_ = sp.symbols('s u')
    rat4 = sp.factor(sp.simplify((f_plus(4) / f_plus(2)).subs({up: m ** 4 / nu ** 4, vp: (m ** 2 - nu ** 2) ** 2 / nu ** 4})))
    rat5 = sp.factor(sp.simplify((f_minus(5) / f_minus(3)).subs({up: m ** 4 / nu ** 4, vp: (m ** 2 - nu ** 2) ** 2 / nu ** 4})))
    pr(f"    f_4/f_2|crystal = {rat4} = u^2/(s(s+u)) at s = m^2, u = -nu^2: {sp.simplify(rat4.subs({m: sp.sqrt(s_), nu: sp.sqrt(-u_)}) - u_ ** 2 / (s_ * (s_ + u_))) == 0}")
    pr(f"    f_5/f_3|crystal = {rat5} = u^2/(s(s+u)): {sp.simplify(rat5.subs({m: sp.sqrt(s_), nu: sp.sqrt(-u_)}) - u_ ** 2 / (s_ * (s_ + u_))) == 0}")
    pr("    => at the crystal the k = 4 (5) functional is the k = 2 (3) one times u^2/(s(s+u)): two extra flat subtractions, as for B_4 vs B_2")
    pr("       (CHLPS kernel (2s+u)/(s(s+u))^{k/2+1}); the (2,4) dual structure and the (1,3) celestial structure are untouched, so EE's crystal")
    pr("       identity EE = q^2 (EE_FUNCTIONAL.md §1.2) and hence flip-freeness and the q^2 stripped density carry over verbatim.")
    open(os.path.join(CACHE, 'stage12_kernel.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
