"""
u2 (EE_FUNCTIONAL.md §2): EXACT spinning magic-block coefficients h_-+(nu) for any j' = 0 (2,4) structure of s1_offdiag_kernel.STRUCTS
(q2, EE, EH, Nq, ...), by the v4_hclosed method (H7 §23.1: equivariant gradient at one celestial point, rising-factorial K ratios),
with the (2,4) monomials of the structure in place of v3's q^2.  All structures share K_A^ref, K_B^ref (scalar group), so the rational
R_-+^S of different structures are directly comparable and h_-+ of a linear combination is the same combination of R_-+.
Checks: q2 reproduces v4_hclosed's R^T exactly; numerical comparison with the s3 forward ratio n_EE/n_q2 (img_vec, beta = 0).
Analysis per structure / combination: real zeros of h_-, h_+, h_- + h_+; rho = h_+/h_-; w = Z_+(rho) - P_+(rho) (H7 §24.2) -> null dim
(w-1)/2; small-nu s0 and large-nu sigma of u = (h_- - h_+)/(h_- + h_+).
Usage: python3 u2_hclosed_struct.py -> cache/u2_hclosed.txt, cache/u2_hclosed.pkl (srepr of R_-+ per structure)
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
import sympy as sp
import v3_regge_kernel as V
import v4_hclosed as H4
import s1_offdiag_kernel as S1

d = V._dl
nu_s = H4.nu_s
I = sp.I
OUTP = os.path.join(HERE, 'cache', 'u2_hclosed.pkl')


def monos_struct_sym(name):
    expr = sp.expand(16*V._X**(1 - d + 4)*V._s**(-(1 - d) - 4)*sp.expand(S1.STRUCTS[name]))
    out = []
    for term in sp.Add.make_args(expr):
        pw = term.as_powers_dict()
        k = sp.simplify(pw.get(V._X, 0)); m = sp.simplify(pw.get(V._s, 0)); i = int(pw.get(V._ex, 0)); j = int(pw.get(V._fx, 0))
        coef = sp.simplify(term/(V._X**k*V._s**m*V._ex**i*V._fx**j))
        out.append((coef, k, m, i, j))
    return out


def side_coeffs_A(name):
    """copy of v4_hclosed.side_coeffs(side='A') with the monomials of STRUCTS[name]."""
    z = (1, 0, 1)
    s2 = 1/sp.sqrt(2)
    em = (0, s2, -I*s2); ep = (0, s2, I*s2)
    e, f = em, ep
    kref, mref = 1 - d, -(1 - d)
    cargs = (H4._md(z, e), H4._md(z, f), H4._md(e, f))
    monos = monos_struct_sym(name)
    cref = 3 + kref; lref = -mref
    V3 = [0, 0, 0]
    cache = {}
    for (coef, kf, mf, i, j) in monos:
        k = sp.sympify(kf); m = sp.sympify(mf)
        c = 3 + k; lam = -m; kap = (2*lam - c - 3)/2
        cv = coef.subs({V._a2: cargs[0], V._a4: cargs[1], V._e24: cargs[2]})
        n1 = sp.simplify(c/2 - cref/2); n2 = sp.simplify((c + 3)/2 - lam - (cref + 3)/2 + lref)
        dmu = sp.simplify((c + 3 - lam) - (cref + 3 - lref))
        assert n1.is_integer and n2.is_integer and dmu.is_integer, (n1, n2, dmu)
        Kr = sp.Integer(2)**dmu*H4._gratio(1 + cref/2, n1)*H4._gratio((cref + 3)/2 - lref, n2)
        for mu in range(3):
            key = (i, j, mu)
            if key not in cache:
                cache[key] = H4._dpat(e, f, z, i, j, mu)
            D = cache[key].subs({H4._lam: lam, H4._kap: kap})
            V3[mu] += cv*Kr*D
    V3 = [sp.factor(sp.simplify(v)) for v in V3]
    assert sp.simplify(V3[1]) == 0, V3[1]
    c2 = V3[2]; c1 = sp.factor(sp.simplify(-V3[0] - c2))
    return c1, c2


def R_pm_struct(name, Bc):
    c1, c2 = side_coeffs_A(name); d1, d2 = Bc
    l = d - 1
    Rp = sp.factor(sp.simplify(c1*d1*(l + 1)/(2*l + 1) + c1*d2))
    Rm = sp.factor(sp.simplify(c1*d1*l/(2*l + 1) + c2*d1))
    sub = {d: sp.Rational(1, 2) + I*nu_s}
    return sp.factor(sp.simplify(Rm.subs(sub))), sp.factor(sp.simplify(Rp.subs(sub)))


NAMES = ['q2', 'EE', 'EH', 'Nq']


def build():
    Bc = H4.side_coeffs(True, 'B')
    D = {}
    for nm in NAMES:
        t0 = time.time()
        D[nm] = R_pm_struct(nm, Bc)
        print(f"  {nm}: built in {time.time() - t0:.0f}s", flush=True)
    pickle.dump({k: (sp.srepr(a), sp.srepr(b)) for k, (a, b) in D.items()}, open(OUTP, 'wb'))
    return D


def load():
    if not os.path.exists(OUTP):
        return build()
    D = pickle.load(open(OUTP, 'rb'))
    return {k: (eval(a, {**sp.__dict__}), eval(b, {**sp.__dict__})) for k, (a, b) in D.items()}


def numfuncs():
    """{name: (Rm(nu), Rp(nu))} numpy callables (complex nu ok). h_-+ = N_com * R (CHMRS normalisation), N_com even, positive."""
    return {k: (sp.lambdify(nu_s, a, 'numpy'), sp.lambdify(nu_s, b, 'numpy')) for k, (a, b) in load().items()}


def roots_of(expr):
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    rz = np.roots([complex(c) for c in sp.Poly(sp.expand(num), nu_s).all_coeffs()]) if sp.degree(num, nu_s) > 0 else np.array([])
    rp = np.roots([complex(c) for c in sp.Poly(sp.expand(den), nu_s).all_coeffs()]) if sp.degree(den, nu_s) > 0 else np.array([])
    return rz, rp


def winding(rho_expr):
    """w = Delta_{(0,inf)} arg rho / pi, numerically (fine log grid) and via Z_+ - P_+ (H7 §24.2)."""
    f = sp.lambdify(nu_s, rho_expr, 'numpy')
    nv = np.concatenate([np.geomspace(1e-6, 1e-2, 4000), np.linspace(1e-2, 50, 400000), np.geomspace(50, 1e7, 4000)])
    ph = np.unwrap(np.angle(f(nv)))
    rz, rp = roots_of(rho_expr)
    tol = 1e-9
    return (ph[-1] - ph[0])/np.pi, int(np.sum(rz.imag > tol)) - int(np.sum(rp.imag > tol)), rz, rp


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    D = build()
    pr("# u2: exact h_-+ (rational part R_-+, h = N_com R) for j' = 0 (2,4) structures; v4_hclosed method")
    ref = H4.load()
    pr(f"q2 check vs v4_hclosed: R_-: {sp.simplify(D['q2'][0] - ref['RmT'])}, R_+: {sp.simplify(D['q2'][1] - ref['RpT'])}")
    for nm in NAMES:
        pr(f"\n{nm}: R_- = {D[nm][0]}\n    R_+ = {D[nm][1]}")
    # numerical cross-check of the forward ratio n_EE/n_q2 with img_vec (beta = 0); forward image ~ h_- + h_+ (times nu-only common)
    pr("\nforward ratio (h_-+h_+)_S/(h_-+h_+)_q2 exact vs s1.img_vec (beta = 0):")
    F = {k: (sp.lambdify(nu_s, a, 'numpy'), sp.lambdify(nu_s, b, 'numpy')) for k, (a, b) in D.items()}
    for nm in ('EE', 'EH', 'Nq'):
        for nu in (1.0, 1.17, 3.0, 6.64, 20.0):
            ex = (F[nm][0](nu) + F[nm][1](nu))/(F['q2'][0](nu) + F['q2'][1](nu))
            nmr = S1.img_vec(nu, np.array([0.0]), S1.PRES, Amonos=S1.A_monos_struct(nm))[0]/S1.img_vec(nu, np.array([0.0]), S1.PRES)[0]
            pr(f"   {nm} nu={nu:6.2f}: exact {complex(ex):+.8f}  img_vec {complex(nmr):+.8f}")
    pr("\nzeros/poles (|Im nu| <= 2 listed) and real-axis structure:")
    for nm in ('q2', 'EE', 'EH'):
        Rm, Rp = D[nm]
        for lab, ex in (('h_-', Rm), ('h_+', Rp), ('h_- + h_+', Rm + Rp)):
            rz, rp = roots_of(ex)
            fm = lambda r: ", ".join(f"{z.real:+.5f}{z.imag:+.5f}i" for z in sorted(r, key=lambda z: (round(z.imag, 6), z.real)) if abs(z.imag) <= 2)
            pr(f"  {nm} {lab}: zeros [{fm(rz)}]; poles [{fm(rp)}]; real zeros (|Im| < 1e-8, Re > 0): "
               f"{[round(z.real, 6) for z in rz if abs(z.imag) < 1e-8 and z.real > 0]}")
        w_num, w_zp, rz, rp = winding(sp.factor(Rp/Rm))
        pr(f"  {nm}: winding of rho = h_+/h_- on (0, inf): numeric {w_num:.4f}, Z_+ - P_+ = {w_zp}  -> null dim (w-1)/2 = {(w_zp - 1)/2}")
        s0 = sp.limit(-I*nu_s*(Rm - Rp)/(Rm + Rp), nu_s, 0)
        sl = sp.limit(-I*nu_s*(Rm - Rp)/(Rm + Rp), nu_s, sp.oo)
        pr(f"  {nm}: u = i s/nu: s(0) = {sp.nsimplify(s0)}, s(inf) = {sl}")
    open(os.path.join(HERE, 'cache', 'u2_hclosed.txt'), 'w').write(buf.getvalue())
