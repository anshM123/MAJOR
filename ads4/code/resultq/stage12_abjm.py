"""Result Q, Stage 2 (RESULT_Q_STAGE12.md §3): normalisation chain for placing the ABJM (IIA) R^4 coefficient.

Exact (sympy) checks:
 (1) BCP 1906.07195 eq. (flatSSSS): the <SSSS> flat amplitudes A^i in terms of the MHV function f_1; M_SG top <-> f_1 = 1/(stu),
     M_4 top <-> f_1 = 1 (their eqs. (SGLimit), (largeA4), read from the source in the scratchpad).
 (2) Penedones/BCP flat-space formula (FlatLimit): A = N sqrt(pi) int dalpha/(2 pi i) e^alpha alpha^{-1/2} M(L^2 s/(2 alpha), ...)
     maps a homogeneous degree-n Mellin term c (s,t)^n to N sqrt(pi) c (L^2/2)^n / Gamma(n + 1/2)  (h = sum Delta/2 - d/2 = 1/2 for Delta_S = 1).
 (3) With the IIA S-matrix A = A_SG [1 + zeta(3) l_s^6 stu/32 + ...] and L^8/l_s^8 = 4 pi^4 lambda^2 (BCP (class)):
     B^4_4/B^1_SG = 105 zeta(3)/(64 sqrt2 pi^3)  (BCP (IIARes1)) is reproduced exactly.
 (4) Result-Q variables: f = 8 pi G/(stu) (1 + g4/(8 pi G) stu + ...) => g4/(8 pi G) = zeta(3) alpha'^3/32 (l_s^2 = alpha'), M^2 = 4/alpha'
     (first massive level), Lam = M L = 2^{5/4} pi^{1/2} lambda^{1/4}  =>  g4 Lam^6/(8 pi G) = 2 zeta(3) exactly (leading order in 1/lambda),
     and directly from the CFT data: g4 Lam^6/(8 pi G) = (B^4_4/B^1_SG) lambda^{-3/2} Lam^6 Gamma(3/2)/Gamma(9/2) (L^2/2)^3 / L^6 ... = 2 zeta(3).
 (5) Placement against the Stage-0 bounds A log(2 Lam/3) + B (functionals_q.json), lambda -> Lam, with the lambda - 1/24 shift and the
     genus-one (finite mu) term as side information.
Output: cache/stage12_abjm.txt.  Seconds.
"""
import os, io, json
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def main():
    s, t, L, ls, lam, z3, al = sp.symbols('s t L l_s lambda zeta3 alpha', positive=True)
    u = -s - t
    pr("# stage12_abjm: normalisation chain for the ABJM R^4 placement (exact, sympy)")
    # (1) BCP flatSSSS
    def Avec(f):
        return [sp.simplify(-sp.Rational(1, 2) * t * u * (-s ** 2 * f(s, t) + u ** 2 * f(u, s) + t ** 2 * f(t, s))),
                sp.simplify(-sp.Rational(1, 2) * s * u * (s ** 2 * f(s, t) + u ** 2 * f(u, s) - t ** 2 * f(t, s))),
                sp.simplify(-sp.Rational(1, 2) * t * s * (s ** 2 * f(s, t) - u ** 2 * f(u, s) + t ** 2 * f(t, s))),
                sp.simplify(-sp.Rational(1, 2) * s * t * u * (u * f(u, s) + t * f(t, s))),
                sp.simplify(-sp.Rational(1, 2) * s * t * u * (u * f(u, s) + s * f(s, t))),
                sp.simplify(-sp.Rational(1, 2) * s * t * u * (s * f(s, t) + t * f(t, s)))]
    def fSG(a, b): return 1 / (a * b * (-a - b))
    def f4(a, b): return sp.Integer(1)
    ASG = Avec(fSG); A4 = Avec(f4)
    pr("(1) BCP (flatSSSS) with f1 = 1/(stu):", [sp.factor(x) for x in ASG])
    pr("    expected (SGLimit, g1 = 1): [tu/s, su/t, st/u, s/2, t/2, u/2] ->",
       all(sp.simplify(a - b) == 0 for a, b in zip(ASG, [t * u / s, s * u / t, s * t / u, s / 2, t / 2, u / 2])))
    pr("    with f1 = 1:", [sp.factor(x) for x in A4])
    pr("    expected (largeA4): [t^2u^2, s^2u^2, s^2t^2, s^2tu/2, st^2u/2, stu^2/2] ->",
       all(sp.simplify(a - b) == 0 for a, b in zip(A4, [t ** 2 * u ** 2, s ** 2 * u ** 2, s ** 2 * t ** 2, s ** 2 * t * u / 2, s * t ** 2 * u / 2, s * t * u ** 2 / 2])))
    pr("    and A4^i / A_SG^i = stu for every i:", all(sp.simplify(a / b - s * t * u) == 0 for a, b in zip(A4, ASG)))
    # (2) Penedones map of a degree-n monomial
    n = sp.Symbol('n', positive=True, integer=True)
    alpha = sp.Symbol('alpha')
    pref = lambda nn: sp.sqrt(sp.pi) * (L ** 2 / 2) ** nn / sp.gamma(nn + sp.Rational(1, 2))
    pr("(2) Penedones/BCP (FlatLimit): degree-n term c (s,t)^n -> N sqrt(pi) c (L^2/2)^n/Gamma(n+1/2)  [1/(2 pi i) int e^a a^{-1/2-n} = 1/Gamma(n+1/2)]")
    ratio_map = sp.simplify(pref(4) / pref(1))
    pr("    degree-4 / degree-1 map ratio:", sp.factor(ratio_map), " = (L^2/2)^3 Gamma(3/2)/Gamma(9/2) = L^6/105")
    # (3) BCP (IIARes1)
    # A = norm*(B_SG M_SG-map + lambda^{-3/2} B4 M4-map): A4-part/ASG-part = lambda^{-3/2} (B4/BSG) ratio_map * stu  must equal zeta3 l_s^6 stu/32
    ls_over_L6 = sp.Integer(1) / (2 * sp.sqrt(2) * sp.pi ** 3 * lam ** sp.Rational(3, 2))          # (l_s/L)^6 from L^8/l_s^8 = 4 pi^4 lambda^2
    chk = sp.simplify(((4 * sp.pi ** 4 * lam ** 2) ** sp.Rational(-3, 4)) - ls_over_L6)
    pr("    (l_s/L)^6 = (4 pi^4 lambda^2)^{-3/4} = 1/(2 sqrt2 pi^3 lambda^{3/2}):", chk == 0)
    B4overBSG = sp.simplify(z3 / 32 * ls_over_L6 * L ** 6 / ratio_map * lam ** sp.Rational(3, 2))
    pr("    => B^4_4/B^1_SG =", sp.nsimplify(B4overBSG), ";  BCP (IIARes1): 105 zeta(3)/(64 sqrt2 pi^3);  equal:",
       sp.simplify(B4overBSG - 105 * z3 / (64 * sp.sqrt(2) * sp.pi ** 3)) == 0)
    # (4) Result-Q variables
    g4_over_8piG = z3 * al ** 3 / 32                # l_s^2 = alpha'
    M2 = 4 / al
    pr("(4) g4/(8 pi G) = zeta(3) alpha'^3/32;  M^2 = 4/alpha' (first massive level; VS: s = M^2 residue = pure (+-), J = 4, RESULT_Q_STAGE0 §5)")
    pr("    g4 M^6/(8 pi G) =", sp.simplify(g4_over_8piG * M2 ** 3), "(= 2 zeta(3))")
    Lam = sp.sqrt(4 * L ** 2 / al)
    L2_over_al = sp.sqrt(2) * sp.pi * sp.sqrt(lam)    # L^2/l_s^2 = (4 pi^4 lambda^2)^{1/4}
    Lam_l = sp.simplify(sp.sqrt(4 * L2_over_al))
    pr("    Lam = M L = sqrt(4 L^2/alpha') =", Lam_l, "=", sp.N(Lam_l.subs(lam, 1), 6), "lambda^{1/4}")
    # directly from CFT data: g4 Lam^6/(8 pi G) = lambda^{-3/2} (B4/BSG) Lam^6 * [Gamma(3/2)/Gamma(9/2) (1/2)^3]  (L = 1, AdS units)
    direct = sp.simplify(lam ** sp.Rational(-3, 2) * (105 * z3 / (64 * sp.sqrt(2) * sp.pi ** 3)) * Lam_l ** 6 * sp.gamma(sp.Rational(3, 2)) / sp.gamma(sp.Rational(9, 2)) / 8)
    pr("    from CFT data (AdS units L = 1): g4 Lam^6/(8 pi G) = lambda^{-3/2} (B^4_4/B^1_SG) Lam^6 Gamma(3/2)/(8 Gamma(9/2)) =", direct)
    pr("    i.e. the R^4 / Einstein ratio of top-degree Mellin coefficients in <SSSS>, rho_M := lambda^{-3/2} B^4_4/B^1_SG, gives")
    pr("         g4/(8 pi G) [AdS units] = rho_M Gamma(h+1)/Gamma(h+4) (1/2)^3 = rho_M/105   (h = 1/2).")
    pr("    general rule (same derivation): for any component correlator with h = sum Delta_i/2 - d/2 whose SG term has top degree n_SG and R^4")
    pr("    term top degree n_SG + 3 (same flat helicity structure): g4/(8 pi G) = rho_M (1/2)^3 Gamma(h + n_SG)/Gamma(h + n_SG + 3).")
    # (5) placement
    pr("(5) placement against Stage-0 bounds g4 Lam^6/(8 pi G) <= A log(2 Lam/3) + B (flat image; the AdS bound itself is Stage 3):")
    F = json.load(open(os.path.join(HERE, 'cache', 'functionals_q.json')))
    import math
    z3n = 1.2020569031595942
    pr("    lambda       Lam(lambda)  Lam(lambda-1/24)   2zeta3   robust S_min+B3   ratio    robust S_min   ratio")
    for lv in (1e3, 1e4, 3.2e5, 1e6, 2.6e7, 1e8, 3.2e9):
        La = 2 ** 1.25 * math.sqrt(math.pi) * lv ** 0.25
        Lb = 2 ** 1.25 * math.sqrt(math.pi) * (lv - 1 / 24) ** 0.25
        b1 = F['Q0-full-robust']['A'] * math.log(2 * La / 3) + F['Q0-full-robust']['B']
        b2 = F['Q0-min-robust (no B3)']['A'] * math.log(2 * La / 3) + F['Q0-min-robust (no B3)']['B']
        pr(f"    {lv:9.2e}   {La:10.3f}   {Lb:10.3f}       {2 * z3n:.5f}   {b1:8.2f}        {2 * z3n / b1:.4f}   {b2:8.2f}      {2 * z3n / b2:.4f}")
    pr("    finite mu = N/k^5 (string loop, outside A2): 2 zeta(3) -> 2 zeta(3) + (2 pi^2/3) g_s^2 (BCP (IIAAMPAgain): zeta(3)/32 + g_s^2 pi^2/96).")
    pr("    check: 64*(zeta3/32 + g^2 pi^2/96) =", sp.expand(64 * (z3 / 32 + sp.Symbol('g_s') ** 2 * sp.pi ** 2 / 96)))
    open(os.path.join(CACHE, 'stage12_abjm.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
