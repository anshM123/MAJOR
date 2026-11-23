"""deshift/c5_coeffs.py -- DESHIFT_CLOSE.md §1: consistency checks of the stored spinning j' = 1 (o_H q) k = 5 coefficients R_j
(functionals_q5_deshift.json, computed in RESULT_Q_STAGE5 §1 [P + DV]).  This is NOT an independent re-derivation of the R_j (that
needs the j' = 1 light-ray monomials); it checks that the stored rationals have every structural property the de-shift and the bound use:
  (i)   real-analyticity conj R_j(x) = R_j(-x) on real x (so Re[c Psi_5/Nhat5] has a Hermitian extension; sigma = 1);
  (ii)  the 1/nu poles of R_{+-2} cancel in Nhat5 (residues opposite), Nhat5(0) finite and non-zero;
  (iii) Nhat5 has no real zero (exact: real roots of its numerator polynomial);
  (iv)  the raw nu = 0 image coefficients kappa_1..3 recomputed from the stored R_j equal the stored rationals;
  (v)   the poles of R_j in the strips crossed by the de-shift (|Im nu| < 2) and the zeros there (used by the continuum weight ab_c5).
Output: cache/c5_coeffs.txt
"""
import os
import sympy as sp
import kern

nu = kern.nu_s
out = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)


R = {j: sp.together(r.expr) for j, r in kern.load_C5().items()}
x = sp.Symbol('x', real=True)
pr("# c5_coeffs: structural checks of the stored k = 5 (o_H q, j' = 1) leading-Regge coefficients R_j, j = 0, +-2 (++ class)")
for j in sorted(R):
    d = sp.simplify(sp.conjugate(R[j].subs(nu, x)) - R[j].subs(nu, -x))
    pr(f"(i) j = {j:+d}: conj R_j(x) - R_j(-x) = {d}")
N5 = sp.cancel(sp.together(sum(R.values())))
res = {j: sp.simplify(sp.residue(R[j], nu, 0)) for j in (2, -2)}
pr(f"(ii) Res_0 R_2 = {res[2]}, Res_0 R_-2 = {res[-2]} (sum {sp.simplify(res[2] + res[-2])}); Nhat5(0) = {sp.nsimplify(N5.subs(nu, 0))}")
num, den = sp.fraction(sp.factor(N5))
P = sp.Poly(sp.expand(num), nu)
# real zeros of a complex polynomial = common real zeros of Re and Im parts
pr_ = sp.Poly(sp.expand(sp.re(P.as_expr().subs(nu, x))), x); pi_ = sp.Poly(sp.expand(sp.im(P.as_expr().subs(nu, x))), x)
g = sp.gcd(pr_, pi_)
pr(f"(iii) gcd(Re, Im) of Nhat5's numerator on the real line: {g.as_expr()}  (degree 0 => no real zero) ")
# (iv) raw nu = 0 image
h = {j: sp.cancel(R[j]/N5) for j in R}
c = sp.simplify(sp.residue(h[2], nu, 0))
dp = sp.simplify(sp.limit(h[2] - c/nu, nu, 0)); dm = sp.simplify(sp.limit(h[-2] + c/nu, nu, 0))
k1, k2, k3 = sp.nsimplify(sp.simplify(2*sp.I*c)), sp.nsimplify(sp.simplify(dp + dm)), sp.nsimplify(sp.simplify(h[0].subs(nu, 0)))
pr(f"(iv) kappa_1 = {k1}, kappa_2 = {k2}, kappa_3 = {k3}, sum = {sp.simplify(k2 + k3)} (should be 1: H5(0) = 1); "
   f"stored: -192491/1028120, 1252043/1028120, -223923/1028120")
for j in sorted(R):
    n_, d_ = sp.fraction(sp.factor(R[j]))
    zs = [z for z in sp.Poly(sp.expand(n_), nu).nroots(n=20) if abs(sp.re(z)) < 1e-15 and abs(sp.im(z)) < 2.01]
    ps = [z for z in sp.Poly(sp.expand(d_), nu).nroots(n=20) if abs(sp.im(z)) < 2.01]
    pr(f"(v) j = {j:+d}: zeros with |Im| < 2 on the imaginary axis: {[sp.N(z, 8) for z in zs]}; poles with |Im| < 2: {[sp.N(z, 8) for z in ps]}")
open(os.path.join(kern.CACHE, 'c5_coeffs.txt'), 'w').write("\n".join(out) + "\n")
