# Session 7: does higher-point Regge boundedness confine g₃/g₂, and with it β₅? And what string theory says

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [Q] quoted from literature, [A] assumption.
Script: `regge5_counterexample.py` → `out_regge5_counterexample.txt`.

## 1. String theory and supersymmetry

**(a) Type II tree level has no α′⁴ terms at all [Q].**
* The four-point closed superstring amplitude is exp(Σ_k 2ζ(2k+1)(α′/4)^{2k+1}(s^{2k+1} + t^{2k+1} + u^{2k+1})/(2k+1)) times the
  supergravity amplitude. Its expansion runs ζ(3)α′³, ζ(5)α′⁵, ζ(3)²α′⁶, …, with nothing at α′⁴.
* At five points, the closed-string α′ expansion (Schlotterer–Stieberger, Stieberger: single-valued MZVs) has the same feature.
  sv(ζ(2)) = 0, and every weight-4 MZV reduces to ζ(2)², so no α′⁴ term survives.
* So the 10-derivative five-graviton contacts (R⁵, D²R⁴) are **absent at tree level**. Garousi (arXiv:0910.0180) finds the
  absence of R⁵ in the tree-level 5- and 6-graviton analysis.
* At one loop, Richards (arXiv:0807.2421) finds R⁵ and D²R⁵ absent in type II after subtracting the D^{2n}R⁴ pieces; ∇⁴R⁵ is present.

**(b) Supersymmetry.** Elvang–Freedman–Kiermaier (arXiv:1003.5018) and the Elvang–Huang notes (arXiv:1308.1697) [Q]: with 32
supercharges, the R⁵ operator **has no supersymmetrization**, and the D²R⁴ matrix element is ∝ (s+t+u) = 0 on shell. So the
10-derivative five-point graviton contact is forbidden by maximal supersymmetry in flat space.

**(c) Consequence for our β₅ [D-sketch].**
* The I₂I₃ ≅ IJ term has 10 derivatives, so its natural string order is α′⁴, i.e. Δ_g⁻⁸ in AdS units. It vanishes there for type II.
* In an AdS₄ compactification (ABJM / type IIA on CP³), a 4D C⁵ term can only come from higher superinvariants. Examples are
  α′⁵ ζ(5)D⁴R⁴, or R⁴-type invariants dressed with internal curvature, but a dimension-8 superinvariant cannot contain 10-derivative
  terms. So β₅ = O(α′⁵/L¹⁰·L²) = O(Δ_g⁻¹⁰), or smaller.
* By the THEOREM_H table, its planar η/s contribution is −1944β₅ = O(Δ_g⁻¹⁰). The R⁴ contribution is 216G₄ = O(Δ_g⁻⁶). So the
  quintic is suppressed by at least Δ_g⁻⁴ relative to the quartic.
* This is far inside what Theorem H needs (|β₅| ≲ Δ_g⁻⁴/2880).

**(d) Caveats.**
* The AdS₄ completion of the flat-space superinvariant argument (the gauged N = 8 superinvariant classification) is not checked.
* Heterotic and bosonic strings have α′⁴ terms generically, so for them β₅ ≠ 0 at natural size. That is still far below Theorem H's
  requirement, but it is not forced.

**Lesson.** In the known UV completions β₅ is at or below natural power counting, and **maximal supersymmetry removes it at the
natural order**. Supersymmetry is a genuine candidate for the "extra input". It is not available for a generic non-supersymmetric
large-gap CFT₃.

## 2. Five-point Regge boundedness gives relations, not bounds

**Structure [D].**
* At tree level, every five-graviton term linear in g₃ carries exactly one g₂, because the X leg must close. So the five-point
  amplitude is A₅ = A₅^{(0)} + A₅^{(1)}, where A₅^{(1)} is linear in the **interference measure** μ = Σ_X g₂(X)g₃(X) δ(M − M_X)
  (with spin labels), and A₅^{(0)} is independent of g₃.
* Suppose A₅ is Regge-bounded in s₁₂ (fixed other invariants), with growth ≤ s₁₂^{J₀}. A dispersion relation in s₁₂, subtracted
  ⌈J₀⌉ + 1 times, collects poles from every channel whose invariant depends on s₁₂. Matching the low-energy expansion then gives,
  for each superconvergent coefficient,

  (light coefficient: β₅-type contact + GR + g₂²-dressed GR terms) = ∫ dμ(M, J) K(M, J; …),

* with kernels K fixed by the C² and C³ vertex structures. This is a **linear relation** in a **sign-indefinite** measure. It
  has no sign and implies no inequality by itself.

**Homogeneity lemma [P].** Every tree-level constraint is homogeneous in g₃ at fixed multiplicity:
* degree 0 at 4 points;
* degree 1 in μ at 5 points;
* degree 2 at 6 points.

So the rescaling g₃ → κg₃ preserves every constraint whose g₃-independent (light) side vanishes. Confinement of κ can only come
from **inhomogeneous** sum rules, where a heavy g₃-sector moment equals a non-zero g₃-independent light quantity (GR, g₂²-terms).

**Five-point counterexample [P, linear algebra; demo in `regge5_counterexample.py`].**
* Impose any finite set of five-point relations ∫dμ K_j = L_j, together with β₅ = ∫dμ K′, where K′ is linearly independent of
  the K_j, as distinct mass powers are.
* With N ≥ J + 1 heavy species at fixed g₂(n) (so the four-point data are unchanged), solve for μ_n and set g₃(n) = μ_n/g₂(n).
* β₅ can then be any real number. The demo has two relations at residual 1e-14 and β₅ = ±10³ achieved exactly.
* The price is a large six-point diagonal: Σg₃²/M² = 4.5·10⁷ in the demo.
* Caveat [A]: this shows the *moment relations* do not confine β₅. Realizing the measure by an actual Regge-bounded spinning tower
  is not constructed.

## 3. Six-point Regge boundedness would bound β₅ (conditional) [P given the input]

The six-point diagonal is a positive measure: σ = Σ_X g₃(X)² δ(M − M_X) ≥ 0. Suppose two things hold:
* a 3 → 3 improved (superconvergent) sum rule exists whose light side is computable from GR and four-point data, L₆;
* its heavy kernel dominates M⁻², i.e. ∫dσ/M² ≤ κ₆ L₆ for M ≥ Λ (the four-point CHMRS domination principle, which applies because
  σ ≥ 0).

Then the Gram bound gives

  |β₅^{exch}| ≤ 2√(A·κ₆L₆/2),   A = Σ g₂²/(2M²) ≲ G₄ (four-point).

With natural sizes, L₆κ₆ ~ Λ⁻¹⁰ gives |β₅| ≲ Λ⁻⁸ = Δ_g⁻⁸, two-sided. Even the weak version κ₆L₆ ≲ Δ_g⁻² suffices for Theorem H.

**Status of the input.** Unproven: the 3 → 3 Regge bound in the s₁₂₃ channel, analyticity in s₁₂₃, removal of the forward
spectator terms, and multi-structure partial-wave positivity (THEOREM_Q5 §3). The 5- and 6-point positivity literature
(Cheung–Remmen 2505.05553; the six-point preprint 2609.17861) finds PSD structures of the form g₄g₆ ≥ g₅². Those are the Gram
cone, not a 3 → 3 CHMRS.

## 4. Verdict

1. **Five-point Regge boundedness does not bound β₅** [P at the level of the moment relations]. It yields linear relations in a
   sign-indefinite interference measure. Any finite set of them leaves β₅ free at fixed four-point data, via the counterexample
   above.
2. **Six-point Regge boundedness would bound β₅ two-sidedly** [P, conditional], at natural size Δ_g⁻⁸ or at the weakest useful
   size, *if* a 3 → 3 improved sum rule with positive heavy side and GR-computable light side exists. This is the single sharpened
   open problem.
3. **String theory / supersymmetry [Q + D-sketch].** Type II has no α′⁴ terms, and maximal supersymmetry forbids R⁵. So β₅ vanishes
   at its natural order, and the remainder, O(Δ_g⁻¹⁰), satisfies Theorem H with a large margin. **For maximally supersymmetric
   large-gap SCFT₃s (ABJM-type), the quintic obstruction is plausibly removed by supersymmetry**, pending the AdS₄ superinvariant
   check.

**Next steps.**
* Either construct the 3 → 3 superconvergence sum rule (flat space first).
* Or restrict the theorem to supersymmetric CFTs and verify, via the N = 8 gauged-supergravity superinvariants or the 3d N = 8
  superconformal Ward identities for the stress-tensor multiplet, that the five-point stress-tensor contact at the relevant order
  is fixed by lower-point data.
