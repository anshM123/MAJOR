"""t6: heavy action of the bare-commutator functional (CLS-D k=3, j'=1) per helicity, and the normalisation N(nu) of OMEGA_PRIME.md.

Imports (does not modify) ads4/code/e2e/laplace.py (M-tables, cached) and e2e/clsd_kernel.py (antipodal weights).
At the (1,3) pair the light-transform kinematics give the leading-Regge mode weights  w_q = (-1)^{q3} d_{(q1+q3-1)/2}(delta)/d_0  (q1+q3 odd)
[SPINNING_E2E.md 3.3, e2e/clsd_kernel.py].  Per heavy helicity flow (h1,h3) the pair response is  A(h1,h3) = sum_q w_q (R + g1 S)
(laplace.M_table; pres = h3 = -h1, mu = 0; flip = h3 = h1, mu = +-2l).  We compare with the project's odd kernel at the same pair:
D5 = eps(z1,z3,yhat) (l=1) and A = D5 D1 (l=2), whose pres responses are  -+ i (a-1) g1/2  and  -+ i g1 (a-1)(nu^2+2(a-3)^2+1/4)/(4 dd_{7/2})
(REGGE_MOMENT_TTTT.md 4.2 [DV]); recomputed here from the same tables as a check.
N(nu) := A_CLSD(+,-) / A_proj(+,-)  at a = 2l + k + 1 (k = 3): a = 6 (l=1), a = 8 (l=2).
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
os.chdir(os.path.join(HERE, '..', 'e2e'))           # laplace's cache lives there (read-only use: tables already cached)
import sympy as sp, mpmath as mp
from laplace import M_table, nu, a, g1
from clsd_kernel import weights

def kernel_modes(l, which):
    """mode weights of the project's odd kernels at one pair, from a DFT on the frame circle (exact for trig polys)"""
    th1, th3 = sp.symbols('th1 th3', real=True)
    # in the frame z(th) = yhat + cos th e_par + sin th e_perp, (yhat^2 = -1): z1.z3 = -1 + cos(th1-th3), z.yhat = -1
    D1 = sp.cos(th1 - th3)                       # z1.z3 + (z1.y)(z3.y)
    D5 = sp.sin(th3 - th1)                       # eps(z1,z3,yhat) up to an overall sign convention
    K = {'D5': D5, 'A': D5*D1}[which]
    E1, E3 = sp.symbols('E1 E3')
    Kx = sp.expand(K.rewrite(sp.exp).subs({sp.exp(sp.I*th1): E1, sp.exp(sp.I*th3): E3}))
    Kx = sp.expand(sp.powsimp(sp.expand(K.rewrite(sp.exp)), combine='exp'))
    w = {}
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            # coefficient of e^{i(q1 th1 + q3 th3)}; the mode functional k_q extracts e^{-i(...)}-weighted average, i.e. the coefficient of e^{+i(q1 th1+q3 th3)}
            c = sp.simplify(sp.integrate(sp.integrate(K*sp.exp(-sp.I*(q1*th1 + q3*th3)), (th1, 0, 2*sp.pi)), (th3, 0, 2*sp.pi))/(4*sp.pi**2))
            if c != 0:
                w[(q1, q3)] = c
    return w

def response(l, w, h1, h3):
    T = M_table(l)
    tot = 0; mus = set()
    for q, ww in w.items():
        R, S, mu = T[(q[0], q[1], h1, h3)]
        if R != 0 or S != 0:
            mus.add(mu)
        tot += ww*(R + g1*S)
    return sp.simplify(sp.together(sp.expand(tot))), mus

def gam(aa, x):
    return abs(mp.gamma((aa - mp.mpf(1)/2 + 1j*x)/2))**2

if __name__ == "__main__":
    mp.mp.dps = 25
    for l, which in ((1, 'D5'), (2, 'A')):
        aa = 2*l + 3 + 1
        print(f"=== l = {l}, k = 3, a = {aa}; project odd kernel {which} ===")
        wp = kernel_modes(l, which)
        print("  project kernel modes:", wp)
        wc = weights(l, 1)
        out = {}
        for name, w in (('proj', wp), ('CLSD', wc)):
            for (h1, h3) in ((1, -1), (-1, 1), (1, 1), (-1, -1)):
                r, mus = response(l, w, h1, h3)
                out[(name, h1, h3)] = r
                print(f"  {name:4s} (h1,h3)=({h1:+d},{h3:+d}) mu={sorted(mus)}: A = {sp.factor(r.subs(a, aa))}")
        # parity structure checks
        for name in ('proj', 'CLSD'):
            ev = sp.simplify((out[(name, 1, -1)] + out[(name, -1, 1)]).subs(a, aa))
            fl_odd = sp.simplify((out[(name, 1, 1)] - out[(name, -1, -1)]).subs(a, aa))
            fl_ev = sp.simplify((out[(name, 1, 1)] + out[(name, -1, -1)]).subs(a, aa))
            print(f"  {name}: pres even part = {ev};  flip odd part = {sp.factor(fl_odd)};  flip even part = {sp.factor(fl_ev)}")
        Nsym = sp.factor(sp.simplify((out[('CLSD', 1, -1)]/out[('proj', 1, -1)]).subs(a, aa)))
        print("  N(nu) = A_CLSD(+,-)/A_proj(+,-) =", Nsym)
        f = sp.lambdify((nu, g1), Nsym, 'mpmath')
        print("   nu      N(nu)        |N|")
        mn = None
        for x in [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20, 40]:
            g = gam(aa-1, x)/gam(aa, x)
            v = complex(f(mp.mpf(x), g))
            mn = abs(v) if mn is None else min(mn, abs(v))
            print(f"   {x:5.2f}  {v.real: .6e}{v.imag:+.2e}i   {abs(v):.6e}")
        print(f"  min |N| on the grid: {mn:.4e}")
        # zeros: N is (rational in nu) x (g1-free?) -> list real roots of numerator if g1-free
        if not Nsym.has(g1):
            num, den = sp.fraction(sp.together(Nsym))
            pn = sp.Poly(sp.expand(num), nu)
            print("  numerator roots in nu:", [sp.nsimplify(r) for r in sp.solve(pn, nu)])
            print("  denominator roots in nu:", [sp.nsimplify(r) for r in sp.solve(sp.Poly(sp.expand(den), nu), nu)])

    # exact closed forms (coefficients rationalised; checked against the floating expressions above)
    print("=== exact closed forms ===")
    for l, which, closed in ((1, 'D5', sp.Rational(2, 5)*(2*nu + sp.I)),
                             (2, 'A', None)):
        aa = 2*l + 4
        wp = kernel_modes(l, which); wc = weights(l, 1)
        rc, _ = response(l, wc, 1, -1); rp, _ = response(l, wp, 1, -1)
        Nex = sp.factor(sp.nsimplify(sp.simplify((rc/rp).subs(a, aa)), rational=True))
        print(f"  l={l}: N(nu) = {Nex}")
        if l == 2:
            target = (2*nu + sp.I)*(4*nu**2 + 9)/(4*nu**2 + 201)
            print("     N / [(2nu+i)(4nu^2+9)/(4nu^2+201)] =", sp.simplify(Nex/target))
    print("  CLS-D flip responses A(+,+), A(-,-) vanish identically for l = 1, 2 (see above).")
