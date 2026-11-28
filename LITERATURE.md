# LITERATURE (status audit, 15 November 2026)

Marks: **[R]** = the source text was read (by me or by a delegated reader quoting it); **[A]** = abstract page
only; **[S]** = search snippet or secondary source only. Vendor-posted AI preprints are flagged as such; none
has independent confirmation that I could find.

## 1. AI-assisted results checked so as not to duplicate them

| Result | Source | Claim | Status |
|---|---|---|---|
| Hawking–Page transition at weak coupling | Anthropic, hexagon:2611.00175v1, 7 Nov 2026 [R] | For N=4 SYM on S³ at large N the quartic Polyakov-loop coefficient is negative at order λ², so deconfinement is first order just below the Hagedorn temperature | Numerical perturbative computation, not a theorem; one diagram class disagrees with earlier literature; unrefereed |
| Spinless modular bootstrap at large c | Anthropic, hexagon:2611.00174v1, 7 Nov 2026 [R] | Optimal slope `Λ = 12π²/(π²+4 arcsinh²1) ≈ 9.127` in `Δ₁ ≤ c/Λ + o(c)` over a class of functionals | "Proved"/"computer-assisted" by its own labels; unrefereed; says nothing about spin, integrality, or pure gravity |
| Conditional `c/12` gap bound | Perlmutter, arXiv:2610.31413 [R] | `Δ₁ − (c−1)/12 ≤ o(c)` under two stated assumptions | Conditional |
| HKS mixed-temperature conjecture | Dey–Pal–Qiao, arXiv:2410.18174 [A] | Universality for all `β_Lβ_R ≠ 4π²` given sparseness | Proved (JHEP 2025) |
| OpenAI `math` release | github.com/openai/math, 6 Oct 2026 [R: contents] | 721 manuscripts; relevant families: Penrose inequalities incl. AdS, strong cosmic censorship near Kerr, BFSS bound state, scale ⇒ conformal in 4D | "Raw research data"; ~42% formalised; nothing on holography or black-hole information |
| 4D gauge-theory existence | OpenAI announcement, 26 Oct 2026 [S] | Existence of 4D gauge + fermion QFTs | Artifacts withheld; unverifiable |

No AI-assisted claim was found on the Page curve, replica wormholes, factorisation, baby universes, islands,
non-isometric codes, or the covariant entropy cone.

## 2. The covariant holographic entropy cone (the chosen target)

| Paper | What it establishes | Relation to this work |
|---|---|---|
| Bao, Nezami, Ooguri, Stoica, Sully, Walter, arXiv:1505.07839 | RT inequalities via contraction maps on a time-symmetric slice | Supplies the notion of contraction map used in Lemma A |
| Wall, arXiv:1211.3494 | Maximin; strong subadditivity and monogamy for HRT | Used only to justify Hypothesis H |
| Czech, Dong, arXiv:1905.03787 [A] | Claimed all contraction-map inequalities for HRT in 2+1 via kinematic space | Valid for simply connected bulk; error for multiply connected bulk (next row) |
| Grado-White, Grimaldi, Headrick, Hubeny, arXiv:2407.07165v2 [R, full text] | Error in Czech–Dong for multiply connected bulk ("an infinite number of surfaces enter on both sides of the inequality"); an "argument" for π₁ = ℤ; numerics for the rotating three-boundary and one-boundary-torus vacuum wormholes: 1,877 inequalities, no violation; "The validity of the RT inequalities for multiply-connected 2+1 dimensional bulk spacetimes remains an open question." Geodesic lengths `ℓ(γ) = cosh⁻¹(Tr γ_L/2) + cosh⁻¹(Tr γ_R/2)` | **States the open case proved here for vacuum.** Does not use Crofton, Liouville measure, geodesic currents or intersection numbers (searched) |
| Grado-White, Grimaldi, Headrick, Hubeny, arXiv:2502.09894 [S] | Minimax graph model; works if time-sheets are "cooperating", which may fail | Different method; conditional |
| Zhao, arXiv:2602.04888v2 [R: full PDF text searched] | Conditional no-short-cut theorem: assumes null curvature condition, maximin, AdS-hyperbolicity, and that "for every pair of HRT surfaces, the corresponding exposed region is nonempty"; "configurations in which an interaction region is covered only by a union of others" are left "to future work" | Not restricted to vacuum or 2+1, but conditional. The exposed-region condition is assumed, not proved for any class. The text contains no occurrence of "multiply", "wormhole", "BTZ", "vacuum", "kinematic"; does not cite 2407.07165 |
| Grimaldi, Headrick, Hubeny, Shteyner, arXiv:2601.09987; Czech et al., arXiv:2601.09989 [R/S] | Majorisation tests: inequalities safe under light-cone deformations | Evidence, not a proof |
| Grimaldi, Headrick, arXiv:2610.33905 [R] | Review: "The validity of higher HEIs for the HRT formula remains an open problem." | Confirms open status as of Oct 2026 |
| Bao, Furuya, Naskar, arXiv:2506.18086 [S] | Claim: an inequality is valid for RT iff it has a contraction map | Not assumed here; Theorem 2 is stated for contraction-map inequalities |

Geometry of the vacuum spacetimes: Mess (1990); Barbot, arXiv:math/0510065; Bonsante–Krasnov–Schlenker,
arXiv:math/0610429; Aminneborg–Bengtsson–Holst (1998); Maxfield, arXiv:1412.0687.

## 3. Black-hole information and factorisation (audited; not chosen as primary)

* Saad–Shenker–Stanford–Yao, arXiv:2103.16754 [R]: for the one-time-point SYK model, Gaussianity of `z` for
  `q>2` is "expect[ed] … but [we] have not proven it"; the half-wormhole identity has relative mean-square
  error ∝ `N^{2−q}`. No proof found in follow-ups (2105.08207, 2110.06221 [R]; 2411.10155, 2503.13172 [A]).
  → secondary target.
* τ-scaled spectral form factor, Saad–Stanford–Yang–Yao, arXiv:2210.11565 [R]: verified to order `τ^13`; no
  rigorous proof found.
* Non-isometric codes, Akers–Engelhardt–Harlow–Penington–Vardhan, arXiv:2207.06536 [R]: theorems for Haar
  unitaries; no concentration result for designs or pseudorandom unitaries.
* One-shot quantum focusing conjectures, arXiv:2307.13032 [R]: unproved; restricted focusing proved in JT+QFT
  with counterexamples to the original conjecture, arXiv:2510.13961 [R].
* Two-sided factorisation, Boruch–Iliesiu–Lin–Yan, arXiv:2406.04396 [R]: relies on analytic continuation of
  the resolvent; exactness conjectured.
* Closed universes: Harlow–Usatyuk–Zhao arXiv:2501.02359, Abdalla et al. arXiv:2501.02632, Gesteau
  arXiv:2509.14338 [A]: unsettled and not sharply mathematical.
* Algebraic approach (crossed products): results are conditional or qualitative; no theorem-level Page curve
  found.

## 4. A caution learned in this session

The first phase of this session (archive_qi/) targeted the ordinary-Petz recovery inequality for conditional
mutual information. My literature agent reported it open; it had in fact been settled negatively by certified
three-qubit counterexamples recorded in September 2026 (QIQCOP Zoo, entry op_87c77263c8bab523). The novelty
claims above therefore rest on (i) reading the primary sources' own statements of what is open, and (ii)
full-text searches of the two closest papers, and should still be treated as "no prior proof found", not as
"no prior proof exists".

## Definitions of the homology slice in the literature (checked 15 Nov 2026, ar5iv renderings)

| Source | Slice | Homology condition | Smooth / uniformly spacelike required? |
|---|---|---|---|
| Wall, arXiv:1211.3494, Def. 7–8, footnote 9 | "complete achronal slice", "an AdS-Cauchy surface, using the convention in which a Cauchy surface is permitted to be null" | minimal surface on the slice "homologous to A" | No. Smooth and spacelike are hypotheses of one theorem (Thm. 12) only. |
| Headrick–Hubeny–Lawrence–Rangamani, arXiv:1408.6300, §2.2, footnotes 4, 10 | Cauchy slices "required … to be acausal" | region "of a bulk Cauchy slice" with boundary `A ∪ E_A`: the "spacelike homology" condition | Not in the text read (§§2.1–3.6; §4.3 was not read verbatim). |

Consequence for this work: the slices of PROOFS §6.3 (embedded, invariant, pairwise spacelike separated) are of
the acausal kind. The pleated slices constructed there are not smooth and, at rotating closed geodesics, not
uniformly spacelike; neither source demands that. Theorems C and E in PROOFS §6.3b establish the Cauchy
identification. Proposition S+ (§6.6b, revision 8) now realizes the same cuts on smooth uniformly spacelike
completed Cauchy slices, including for endpoint regions.

The smooth extension input is Bernal–Sánchez, *Further results on the smoothability of Cauchy hypersurfaces
and Cauchy time functions*, Lett. Math. Phys. 77 (2006), 183–197, arXiv:gr-qc/0512095, Theorem 1.1 and
disconnected-case Remark 4.14. The global doubling mechanism is also used by Aké–Flores–Sánchez,
*Structure of globally hyperbolic spacetimes-with-timelike-boundary*, Rev. Mat. Iberoam. 37 (2021), 45–94,
arXiv:1808.04412v4, §5.2. The new application uses the exact smooth AdS conformal double and reflected
geodesic ends, so the prescribed set is compact and the cited extension hypotheses are verified.

Pleated surfaces and spiralling laminations: Thurston's notes, ch. 8; Bonahon, *Ann. Fac. Sci. Toulouse* 5
(1996); Casson–Bleiler, ch. 4–6. In AdS geometry pleated surfaces appear as boundaries of convex cores (Mess);
no source was found that uses a pleated slice through prescribed geodesics to verify an HRT homology
condition (search not exhaustive).
