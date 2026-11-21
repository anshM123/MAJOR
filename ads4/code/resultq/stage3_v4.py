"""Result Q, Stage 3, Task 1 (RESULT_Q_STAGE3.md §1): the validity input V4 for the subtracted k = 4, 5 functionals.

Exact / numerical checks (sympy + numpy; seconds):
 (1) R := f+_{4;2}/f+_{2;2} = f-_{5;2}/f-_{3;2} = ((v'-u')^2 - 1)/(4 u'v')  [one common factor for both]                      [P]
 (2) R in the CHMRS frame: 16 r^2 (eta^2 (1+r^2)^2 - 4 r^2)/(s sb)^2,  s = 1 - 2 r eta + r^2,  sb = 1 + 2 r eta + r^2              [P]
     limits: Regge (r -> 0, eta fixed): 16 eta^2 r^2;  light-cone edge (s -> 0): C(r)/s^2;  OPE corner (u'=z zb, v'=(1-z)(1-zb), z,zb -> 0):
     -(1/z + 1/zb)/2;  large-eta Regge corner (r -> 0 at fixed rho = 2 r eta < 1): 4 rho^2/(1 - rho^2)^2 (bounded, NO gain)          [P]
 (3) single light-ray endpoints (x3 -> x0 along its ray; x3 -> null infinity): R has finite limits (random configurations)         [DV]
 (4) iepsilon bookkeeping: relative order of the pairs (12), (34), (14), (23), (13) in the Wightman orderings 4312, 4132 and in the
     four orderings of each double commutator                                                                                   [P]
 (5) MFT, commutator (Wightman) form: <T(x0+b z3, z3) T(x0+a z1, z1)> with light-ray polarisations: the numerator
     (z3.z1) x^2 - 2 (z3.x)(z1.x) vanishes identically for x = b z3 - a z1, and equals eps^2 (2 (e.z1)(e.z3) + z1.z3) with the
     i eps shift x1 -> x1 - i eps e  =>  the TT commutator is O(eps^4)/(x13^2)^5; its pairing with any smooth test function, times any
     factor with a pole of finite order at x13^2 = 0 (f_4 has at most 1/x13^2), -> 0                                          [P + DV]
Output: cache/stage3_v4.txt
"""
import os, io
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


up, vp = sp.symbols("u' v'", positive=True)
r, eta, z, zb, rho = sp.symbols('r eta z zb rho', positive=True)

f2 = (vp - up)
f4 = ((vp - up) ** 3 - (vp - up)) / (4 * up * vp)
f3 = sp.Integer(-1)
f5 = -((up - vp) ** 2 - 1) / (4 * up * vp)
Rdef = ((vp - up) ** 2 - 1) / (4 * up * vp)


def part1():
    pr("(1) R = f4/f2 = f5/f3:", sp.simplify(f4 / f2 - Rdef) == 0, sp.simplify(f5 / f3 - Rdef) == 0, "  R(u',v') =", Rdef)
    U = (1 - 2 * r * eta + r ** 2) ** 2 / (16 * r ** 2); V = (1 + 2 * r * eta + r ** 2) ** 2 / (16 * r ** 2)
    s = 1 - 2 * r * eta + r ** 2; sb = 1 + 2 * r * eta + r ** 2
    Rf = sp.simplify(Rdef.subs({up: U, vp: V}))
    target = 16 * r ** 2 * (eta ** 2 * (1 + r ** 2) ** 2 - 4 * r ** 2) / (s * sb) ** 2
    pr("(2) frame form R = 16 r^2 (eta^2(1+r^2)^2 - 4r^2)/(s sb)^2:", sp.simplify(Rf - target) == 0)
    pr("    Regge r->0 (eta fixed):", sp.series(target, r, 0, 4).removeO())
    # light-cone edge: eta = (1+r^2)/(2r) - t/(2r)  =>  s = t
    t = sp.symbols('t', positive=True)
    Redge = sp.simplify(target.subs(eta, ((1 + r ** 2) - t) / (2 * r)))
    C = sp.simplify(sp.limit(Redge * t ** 2, t, 0))
    pr("    edge s -> 0: R ~ C(r)/s^2, C(r) =", sp.factor(C), "; C(0) =", C.subs(r, 0), "; C(1) =", sp.simplify(C.subs(r, 1)),
       "; C > 0 on (0,1):", all(float(C.subs(r, x)) > 0 for x in np.linspace(0.01, 0.99, 50)))
    # OPE corner
    Rz = sp.simplify(Rdef.subs({up: z * zb, vp: (1 - z) * (1 - zb)}))
    lead = sp.simplify(sp.limit(Rz * z * zb / (z + zb), z, 0))
    pr("    OPE corner u'=z zb, v'=(1-z)(1-zb): R =", sp.factor(Rz), "; z->0 then zb->0: R*z*zb/(z+zb) ->",
       sp.limit(lead, zb, 0), " i.e. R ~ -(1/z+1/zb)/2")
    # large-eta Regge corner
    Rrho = sp.simplify(target.subs(eta, rho / (2 * r)))
    pr("    large-eta Regge corner (eta = rho/(2r), r -> 0):", sp.simplify(sp.limit(Rrho, r, 0)),
       " (bounded for rho <= c < 1; diverges only at the edge rho -> 1)")
    # numerical sup of s^2 R on the physical region (bounded away from the edge => R bounded)
    rr = np.linspace(1e-3, 0.999, 400)
    worst = 0.
    for x in rr:
        em = (1 + x * x) / (2 * x)
        ee = 1 + (em - 1) * np.linspace(0, 1, 400)
        ss = 1 - 2 * x * ee + x * x; sbb = 1 + 2 * x * ee + x * x
        RR = 16 * x * x * (ee ** 2 * (1 + x * x) ** 2 - 4 * x * x) / (ss * sbb) ** 2
        ok = ss > 0
        worst = max(worst, np.max(np.abs(RR[ok] * ss[ok] ** 2)))
    pr(f"    sup over the region r in (0,1), 1 <= eta <= eta_max(r) of |s^2 R| = {worst:.4f}  (R = O(1/s^2): singular only at the edge)")
    return target


def mink(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def Rpos(x1, x2, x3, x4):
    x12, x34, x13, x24, x14, x23 = [mink(a - b, a - b) for a, b in ((x1, x2), (x3, x4), (x1, x3), (x2, x4), (x1, x4), (x2, x3))]
    U = x12 * x34 / (x13 * x24); V = x14 * x23 / (x13 * x24)
    return ((V - U) ** 2 - 1) / (4 * U * V), V - U


def part3(seed=3, ntest=300):
    rng = np.random.default_rng(seed)
    worst0, worstI, ratios = 0., 0., []
    for _ in range(ntest):
        x0 = rng.normal(size=3); x2 = rng.normal(size=3) * 2; x4 = rng.normal(size=3) * 2
        th1, th3 = rng.uniform(0, 2 * np.pi, 2)
        z1 = np.array([1, np.cos(th1), np.sin(th1)]); z3 = np.array([1, np.cos(th3), np.sin(th3)])
        a = rng.uniform(0.2, 2)
        x1 = x0 + a * z1
        # x3 -> x0 along its ray
        vals = [Rpos(x1, x2, x0 + b * z3, x4)[0] for b in (1e-3, 1e-4, 1e-5)]
        worst0 = max(worst0, abs(vals[2] - vals[1]) / max(abs(vals[2] - vals[0]), 1e-300) if abs(vals[2]) > 1e-8 else 0.)
        # x3 -> null infinity
        vI = [Rpos(x1, x2, x0 + b * z3, x4)[0] for b in (1e5, 1e6, 1e7)]
        worstI = max(worstI, abs(vI[2] - vI[1]) / max(abs(vI[2] - vI[0]), 1e-300) if abs(vI[2]) > 1e-8 else 0.)
        ratios.append(abs(vals[2]))
    pr(f"(3) single-ray endpoints (300 random configs): x3 -> x0: R converges linearly (max over configs of |R(1e-5)-R(1e-4)|/|R(1e-5)-R(1e-3)| = {worst0:.3f}, i.e. ~0.1 for O(b) approach),"
       f" median |R| = {np.median(ratios):.3g};  x3 -> null infinity: max |R(1e7)-R(1e6)|/|R(1e7)-R(1e5)| = {worstI:.3f}  => R bounded at both ends [DV]")


def part4():
    pairs = ['12', '34', '14', '23', '13']

    def rel(order):
        return {p: ('+' if order.index(p[0]) < order.index(p[1]) else '-') for p in pairs}
    pr("(4) relative order (+ = first index to the left) of each pair, Wightman orderings of the commutator <4[3,1]2>:")
    for o in ('4312', '4132'):
        pr(f"    {o}: " + "  ".join(f"({p}){s}" for p, s in rel(o).items()))
    pr("    => (12),(34),(14),(23) [the poles of f = poly(u',v')/(u'v')] carry the SAME i eps in both terms; only (13) differs, and x13^2 > 0")
    pr("       (spacelike) on the whole integration domain: f*G in both orderings are boundary values from one common domain  [P]")
    pr("    double commutator [4,3][1,2] = 4312 - 4321 - 3412 + 3421; signs of (12),(34) in each term:")
    for o in ('4312', '4321', '3412', '3421'):
        d = rel(o); pr(f"       {o}: (12){d['12']} (34){d['34']} (14){d['14']} (23){d['23']}")
    pr("    => Delta_12 Delta_34 acts on f's poles too: a Wick term without x12, x34 (e.g. <41><32>) gives W * Delta_12 Delta_34 f,")
    pr("       a distribution supported on x12^2 = x34^2 = 0 (u' = 0): the 'finitely many non-vanishing double traces' of the s-channel. [P]")


def part5():
    a, b, eps = sp.symbols('alpha beta epsilon', positive=True)
    th1, th3 = sp.symbols('theta1 theta3', real=True)
    z1 = sp.Matrix([1, sp.cos(th1), sp.sin(th1)]); z3 = sp.Matrix([1, sp.cos(th3), sp.sin(th3)]); e = sp.Matrix([1, 0, 0])
    md = lambda p, q: -p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
    x = b * z3 - a * z1 + sp.I * eps * e        # x13 with x1 -> x1 - i eps e
    num = md(z3, z1) * md(x, x) - 2 * md(z3, x) * md(z1, x)
    num_s = sp.simplify(sp.expand(num))
    target = eps ** 2 * (2 * md(e, z1) * md(e, z3) + md(z1, z3))
    pr("(5) light-ray polarised TT: numerator (z3.z1)x13^2 - 2(z3.x13)(z1.x13) at x13 = b z3 - a z1 + i eps e equals",
       "eps^2 (2(e.z1)(e.z3) + z1.z3):", sp.simplify(num_s - target) == 0, "  (eps = 0: identically 0)")
    # 1d pairing demonstration: eps^4 [(t + i eps)^-p - (t - i eps)^-p] against a smooth test function -> 0; control eps^0 -> finite
    phi = lambda t: np.exp(-(t - 0.3) ** 2) * (1 + 0.5 * t)
    tt = np.concatenate([-np.logspace(1.5, -9, 6000), np.logspace(-9, 1.5, 6000)])
    pr("    pairing I_p(eps) = int phi(t) eps^4 [(t+i eps)^-p - (t-i eps)^-p] dt  (p = 5: T 2-pt; p = 6: times f_4's 1/x13^2), control p=5 without eps^4:")
    for p in (5, 6):
        row = []
        for ep in (1e-1, 3e-2, 1e-2, 3e-3):
            g = phi(tt) * ((tt + 1j * ep) ** (-p) - (tt - 1j * ep) ** (-p))
            I = np.trapezoid(g, tt) if hasattr(np, 'trapezoid') else np.trapz(g, tt)
            row.append(f"eps={ep:.0e}: {abs(ep ** 4 * I):.2e} (control {abs(I):.2e})")
        pr(f"     p={p}: " + "; ".join(row))
    pr("    => eps^4 x (finite-order pole) -> 0 [P: all moments int s^n (s+i)^-p ds, n <= p-2, vanish; remainder O(eps^4 log eps)];")
    pr("       the light-ray-polarised <[T3,T1]> is the zero distribution, so Psi_{4,5}[f G_MFT] = 0 in commutator form, for any f   [P]")


if __name__ == '__main__':
    part1(); part3(); part4(); part5()
    open(os.path.join(CACHE, 'stage3_v4.txt'), 'w').write(buf.getvalue())
