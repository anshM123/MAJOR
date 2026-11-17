"""s2: even k = 4, 6 (and k = 2 for reference): normalisation N_k of the superconvergent replacement with
(1,3) = covariant j'=0 celestial kernel (antipodal weights), (2,4) = flip-free even kernel S_ff = D1^2 + D1 D2,
relative to the project kernel D1^2 (x) D1^2 (pres-pres component).  Pair responses (preserved, F-normalised harmonics):
  j'=0 (1,3): from M-tables (s1);  D1^2: N1/(64 dd52 dd92);  D1 D2: -N2/(64 dd52 dd92);  flips D1^2: +f, D1D2: -f (REGGE_MOMENT_TTTT.md Lemma L, [DV] there;
  D1^2 weights re-checked in s0).
Also: flip/pres ratios (F-normalised) of the j'=0 (1,3) kernel and of the project's D1^2; reality/evenness (shadow factor of j'=0 weights)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e')); sys.path.insert(0, HERE); os.chdir(os.path.join(HERE, '..', 'e2e'))
import sympy as sp, mpmath as mp
from laplace import M_table, nu, a, g1
from s1_responses import w_celestial, resp, delta
dd = lambda b: nu**2 + (a - b)**2
N1 = 64*a**4 - 768*a**3 + 64*a**2*nu**2 + 3280*a**2 - 320*a*nu**2 - 5840*a + 16*nu**4 + 488*nu**2 + 3705
N2 = 64*a**3 - 32*a**2*nu**2 - 584*a**2 + 224*a*nu**2 + 1720*a - 16*nu**4 - 424*nu**2 - 1641
M_D1sq = N1/(64*dd(sp.Rational(5, 2))*dd(sp.Rational(9, 2)))
M_Sff = (N1 - N2)/(64*dd(sp.Rational(5, 2))*dd(sp.Rational(9, 2)))
# exact flip of D1^2 in F-normalisation from the tables (s0): (weights 4,4,2/9)
T = M_table(2)
W_D1sq = {(2, -2): 4, (-2, 2): 4, (0, 0): sp.Rational(2, 9)}
for k in (2, 4, 6):
    aa = k + 5
    w0 = w_celestial(0)
    pe = sp.simplify((resp(w0, 1, -1, aa) + resp(w0, -1, 1, aa))/2)
    fe = sp.simplify((resp(w0, 1, 1, aa) + resp(w0, -1, -1, aa))/2)
    pD = sp.simplify((resp(W_D1sq, 1, -1, aa) + resp(W_D1sq, -1, 1, aa))/2)
    fD = sp.simplify((resp(W_D1sq, 1, 1, aa) + resp(W_D1sq, -1, -1, aa))/2)
    assert sp.simplify(pD - M_D1sq.subs(a, aa)) == 0
    Sff = sp.factor(M_Sff.subs(a, aa))
    Nk = sp.factor(sp.simplify(pe*Sff/(M_D1sq.subs(a, aa))**2))
    print(f"=== k = {k}, a = {aa} ===")
    print(f"  (1,3) j'=0 pres = {sp.factor(pe)};  flip = {sp.factor(fe)}")
    print(f"  (2,4) S_ff = D1^2 + D1 D2: pres = {Sff}; flip = 0 (Lemma L: +1 - 1)")
    print(f"  N_{k}(nu) = {Nk}")
    num, den = sp.fraction(sp.together(Nk))
    print(f"  numerator roots (nu): {[complex(r) for r in sp.Poly(sp.expand(num), nu).nroots()]}")
    print(f"  denominator roots (nu): {[complex(r) for r in sp.Poly(sp.expand(den), nu).nroots()]}")
    f = sp.lambdify(nu, Nk, 'mpmath')
    print("  N_k on real nu:", "  ".join(f"{x}: {float(f(x)):.4g}" for x in (0, 0.5, 1, 2, 3, 5, 10, 20, 50)))
    rSC = sp.lambdify(nu, sp.simplify(fe/pe), 'mpmath'); rD = sp.lambdify(nu, sp.simplify(fD/pD), 'mpmath')
    print("  |flip/pres| (1,3) j'=0:", "  ".join(f"{x}: {abs(complex(rSC(x))):.3g}" for x in (0, 0.5, 1, 2, 3, 5, 10, 20, 50)))
    print("  |flip/pres| project D1^2:", "  ".join(f"{x}: {abs(complex(rD(x))):.3g}" for x in (0, 0.5, 1, 2, 3, 5, 10, 20, 50)))
# shadow/reality for j'=0 weights
d = delta; ok = True; okc = True
for q, wq in w_celestial(0).items():
    M = abs(q[0] + q[1])
    lhs = wq.subs(nu, -nu)*sp.rf(d, M)/sp.rf(1 - d, M)
    ok &= sp.simplify(sp.expand_func(lhs - wq)) == 0
    okc &= sp.simplify(sp.expand_complex(wq.subs(nu, -nu) - sp.conjugate(w_celestial(0)[(-q[0], -q[1])]))) == 0
print("j'=0 shadow identity w(-nu)(delta)_M/(1-delta)_M = w(nu):", ok, ";  reality w_q(-nu) = conj w_{-q}(nu):", okc)
