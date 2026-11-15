# AdS4/CFT3 programme — literature and status map (audit of 15 November 2026)

Objective audited: *construct or rigorously constrain a nonperturbative quantum theory of dynamical
four-dimensional gravity using AdS4/CFT3*.

**How to read this file.**

Source marks (what I actually did in this session):

* **[R]** — full text fetched (ar5iv / arXiv HTML) and the relevant passages read or searched; quotations are verbatim.
* **[A]** — arXiv abstract page only (via the arXiv API); the description paraphrases the abstract.
* **[M]** — from memory, *not* checked in this session. Treat as unverified.

Status vocabulary (never mixed):

* **Theorem** — proved with stated hypotheses at the mathematical level of rigour.
* **Derivation** — physics-level argument (controlled approximations, exchange of limits, assumed analyticity, etc.).
* **Conjecture**, **Numerics** — as named.
* "Proved" in a physics paper's own words is recorded as *derivation* unless the paper is written at theorem level.

No semiclassical bulk reconstruction listed below is a nonperturbative definition of quantum gravity; where a
result is only valid at `N = ∞` or to fixed order in `1/N` this is said.

Searches used: arXiv API (id lookups and abstract-phrase queries sorted by date), INSPIRE citation lists for
key papers (1905.03787, 2412.06907, 2212.04893, 2311.04271, 2108.07280, 2411.04183, 2207.06536), full-text
pattern searches of ar5iv renderings, and general web search. Search coverage is *not* exhaustive: "not found"
below means "not found by these searches".

---

## (a) Conditions for a local bulk: large N, large gap (HPPS) and what is proved

| Ref. | Content | Correlator / dimension | Assumptions | Status |
|---|---|---|---|---|
| Heemskerk–Penedones–Polchinski–Sully, arXiv:0907.0151 [A] | Conjecture: any CFT with a planar expansion and a large gap in anomalous dimensions has a local bulk dual. Shown "to first nontrivial order" for an abstract CFT by counting crossing solutions against bulk contact interactions | scalar four-point; d = 2, 4 worked | large N expansion, single low-dimension scalar, truncation in spin | **Conjecture** + derivation at leading order in 1/N² |
| El-Showk–Papadodimas, arXiv:1101.4163 [A] | Factorised low-lying sector ⇒ AdS field representation; factorisation only at large c | general | large c, factorisation | Derivation |
| Fitzpatrick–Kaplan–Poland–Simmons-Duffin, arXiv:1212.3616 [A; eq. (43) R]; Komargodski–Zhiboedov, arXiv:1212.4103 [A] | Double-twist operators at large spin, anomalous dimension `γ(0,ℓ) = γ_0/ℓ^{τ_m}` with `γ_0 = −2 P_m Γ(Δ_φ)² Γ(τ_m+2ℓ_m) / (Γ(Δ_φ−τ_m/2)² Γ(τ_m/2+ℓ_m)²)` (eq. (43) of 1212.3616, verbatim structure) | scalar four-point, d > 2 | twist gap; regularity of large-spin spectrum | Derivation ("prove" in the abstract, but with the regularity assumption) |
| Fitzpatrick–Kaplan–Walters, arXiv:1403.6829 [A] | Long-distance AdS interactions universal; matches AdS gravity semiclassically | d > 2 and d = 2 | as above | Derivation |
| Maldacena–Shenker–Stanford, arXiv:1503.01409 [A] | Chaos bound `λ_L ≤ 2π k_B T/ħ` | thermal large-N systems | "plausible physical assumptions" (their words) | Derivation |
| Maldacena–Simmons-Duffin–Zhiboedov, arXiv:1509.03612 [A] | Bulk-point singularity as a diagnostic of bulk locality; absence nonperturbatively proved in 1+1 d | four-point, Lorentzian | — | Derivation; 2d statement proved by CFT methods |
| Hartman–Jain–Kundu, arXiv:1509.00014 [A] | Causality constraints from crossing in a shockwave configuration | scalar four-point, general d | reflection positivity, OPE | Derivation |
| Caron-Huot, arXiv:1703.00278 [A; intro R] | Lorentzian inversion formula: OPE data analytic in spin for `J > 1`; large-J expansions "with controlled errors"; for sparse spectrum, suppression of higher-derivative bulk interactions | scalar four-point, general d | Regge boundedness, itself derived from OPE + positivity | Derivation (contour argument); not written at theorem level |
| Simmons-Duffin–Stanford–Witten, arXiv:1711.03816 [A] | Spacetime derivation of the inversion formula | scalar, d = 1, 2, higher | Wick rotation | Derivation |
| Kravchuk–Simmons-Duffin, arXiv:1805.00098 [A]; Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov, arXiv:1905.01311 [A] | Light-ray operators; convergent light-ray OPE | general d | — | Derivation |
| Kravchuk–Qiao–Rychkov, arXiv:2001.08778, 2104.02090 [A] | Four-point functions are tempered distributions; Wightman axioms for scalar 4-pt functions from Euclidean CFT axioms. Abstract: Lorentzian properties "may represent extra assumptions" otherwise | scalar four-point, any d | Euclidean unitary CFT axioms | **Theorem** |
| Pal–Qiao–Rychkov, arXiv:2212.04893 (CMP 2023) [A; text searched R] | Twist gap ⇒ a family with `ℓ → ∞`, `τ → 2Δ_φ` | identical scalars, any unitary CFT | twist gap | **Theorem** |
| van Rees, arXiv:2412.06907 (SciPost 2025) [R] | "the spectrum and OPE coefficients at large spin and fixed twist always become that of a generalized free field theory"; all double-twist towers | identical scalars | unitarity, twist gap, discrete spectrum | **Theorem**. Outlook, verbatim: "It will be interesting to understand the size of the corrections to the results presented here." The rate (anomalous dimension) is **not** proved there |
| Mazáč–Rastelli–Zhou, arXiv:1910.12855 [A]; Carmi–Caron-Huot, arXiv:1910.12123 [A]; Penedones–Silva–Zhiboedov, arXiv:1912.11100 [A]; Caron-Huot–Mazáč–Rastelli–Simmons-Duffin, arXiv:2008.04931 [A] | Analytic functionals, conformal dispersion relation, nonperturbative Mellin amplitudes, unified dispersive CFT sum rules | scalar four-point, general d | Regge boundedness (u-channel) | Derivation |
| Caron-Huot–Mazáč–Rastelli–Simmons-Duffin, "AdS bulk locality from sharp CFT bounds", arXiv:2106.10274 [R, parts] | Bounds on bulk Wilson coefficients in terms of `Δ_gap` with the EFT scaling. **AdS4 result** (eq. (236), verbatim): `g_2 + O(1/R²_AdS) ≥ −17.64 (8πG/M²) log(0.124 M R_AdS)`. "AdS possesses a built-in IR regulator" curing the 4D flat-space IR divergence | **scalar** four-point with graviton + light scalar exchange; d = 3 treated explicitly | HPPS-type family; functional positivity established asymptotically ("We claim that ω_AdS is rigorously positive at large M") and partly numerically: "By checking numerical examples, we find that Q_n(β) is everywhere positive … for 0 ≤ n ≤ 2 in d ≥ 4, and for 0 ≤ n ≤ 1 in d = 3" | Derivation with numerically checked steps; the abstract says "We prove a sharp form of this conjecture". Not a theorem in the mathematical sense for finite `Δ_gap` |
| Caron-Huot–Van Duong, arXiv:2011.02957 [A]; Caron-Huot–Li–Parra-Martinez–Simmons-Duffin, arXiv:2201.06602, 2205.01495 [A] | Flat-space dispersive bounds; 4D graviton scattering: two-sided bounds, scaling with M "up to an infrared logarithm"; D ≥ 5 sharp bounds | S-matrix, D = 4 and D ≥ 5 | analyticity, crossing, Regge boundedness of the S-matrix (assumed) | Derivation + numerics |
| Chang–Landau–Simmons-Duffin, arXiv:2311.04271 [R, parts] | Dispersive sum rules for **spinning** operators; dictionary to flat-space graviton sum rules in the bulk-point limit. Abstract: "This dictionary is a crucial step for establishing the HPPS conjecture for stress tensor correlators." Text: "We leave an exploration of the Regge regime for four-point functions of spinning operators to future work." | stress-tensor four-point, general d | bulk-point saddle | Derivation. **HPPS for stress tensors is not established.** None of the 15 papers citing it on INSPIRE (to Nov 2026) claims to complete it (titles checked) |
| Knop–Mazáč, arXiv:2203.11170 [A] | Dispersive sum rules and EFT bounds in AdS2/CFT1 | 1D CFT | — | Derivation |
| Caron-Huot, arXiv:2211.11791 [A] | "Holographic cameras": bulk imaging from four-point functions in excited states | — | holographic state | Derivation |
| Chester–Pavarini–Piazza, arXiv:2609.31576 [A] | Dispersive sum rules + localization in the planar limit | N = 4 SYM, a 4d N = 2 theory (**d = 4**, not d = 3) | supersymmetry | Numerics |

**Summary for (a).** Rigorous: Lorentzian axioms and distributional properties of scalar four-point functions
(Kravchuk–Qiao–Rychkov); existence and GFF-limit of double-twist families (Pal–Qiao–Rychkov, van Rees).
Physics-level but quantitatively sharp: scalar-correlator locality bounds including AdS4 (2106.10274).
Open: (i) HPPS for the stress-tensor four-point function in any dimension, in particular d = 3;
(ii) a theorem-level (non-asymptotic) version of any gapped locality bound; (iii) rigorous rate of approach
of double-twist operators to GFF values.

---

## (b) Emergence of Einstein gravity: ⟨TTT⟩ in d = 3, collider bounds, entanglement

| Ref. | Content | Assumptions | Status |
|---|---|---|---|
| Giombi–Prakash–Yin, arXiv:1104.4317 [A] | In d = 3 conserved-current three-point functions: parity-even structures realised by free fields plus "at most one parity violating structure" per set of spins | conformal symmetry + conservation | Derivation (classification) |
| Maldacena–Zhiboedov, arXiv:1112.1016, 1204.3882 [A] | One higher-spin conserved current in a CFT3 ⇒ correlators of currents are those of free bosons or free fermions; slightly broken case constrained | unique stress tensor, large-N structure for the broken case | Derivation (widely treated as a theorem; not checked here at that level) |
| Hofman–Maldacena, arXiv:0803.1467 [A] | Conformal collider (energy-positivity) bounds | positivity of energy flux | Derivation |
| Hofman–Li–Meltzer–Poland–Rejon-Barrera, arXiv:1603.03771 [A] | Collider bounds hold "for any unitary parity-preserving CFT with a unique stress tensor in spacetime dimensions larger than 2" | as quoted | Derivation (bootstrap proof) |
| Faulkner–Leigh–Parrikar–Wang, arXiv:1605.08072 [A]; Hartman–Kundu–Tajdini, arXiv:1610.05308 [A] | ANEC from relative entropy / from causality; hence collider bounds | Minkowski QFT; interacting, d > 2 | Derivation |
| Chowdhury–David–Prakash, arXiv:1707.03007 [A] | **d = 3**: ⟨TTT⟩ determined by `(t_4, α_T)` with `α_T` parity-odd; collider bounds confine `(t_4, α_T)` to a disc; large-N Chern–Simons-matter theories lie on its boundary | collider bounds | Derivation |
| Meltzer, arXiv:1811.01913 [A] | Higher-spin ANEC; "prove that if the Hofman-Maldacena bounds are saturated, then the theory has a higher spin symmetry" | unitarity, Regge trajectories | Derivation |
| Cordova–Maldacena–Turiaci, arXiv:1710.03199 [A] | ANEC bounds on ⟨TTO⟩ | ANEC | Derivation |
| Dymarsky–Kos–Kravchuk–Poland–Simmons-Duffin, arXiv:1708.05718 [A] | 3d stress-tensor bootstrap (parity-preserving): bounds on `c_T` versus the ⟨TTT⟩ parameter "numerically reproduce the conformal collider bounds" | crossing of ⟨TTTT⟩ | **Numerics** |
| Chang et al., arXiv:2411.15300 [A] | 3d Ising mixed σ, ε, T bootstrap; precise data | — | Numerics (not holographic) |
| Camanho–Edelstein–Maldacena–Zhiboedov, arXiv:1407.5597 [A] | Non-Einstein graviton three-point couplings violate causality unless an infinite higher-spin tower appears | weakly coupled gravity | Derivation |
| Afkhami-Jeddi–Hartman–Kundu–Tajdini, arXiv:1610.09378 [A] | ⟨TTT⟩ approaches the Einstein structure as the higher-spin gap grows | **d = 4**, large N, sparse spectrum | Derivation |
| Meltzer–Perlmutter, arXiv:1712.04861 [A] | ⟨TTO⟩ suppressed by powers of `Δ_gap` | large N, large higher-spin gap | Derivation |
| Belin–Hofman–Mathys, arXiv:1904.05892 [A] | Vanishing ANEC commutator ⇒ `a = c` in d = 4 | large N, large gap | Derivation |
| Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov, arXiv:1904.05905 [A] | Commutativity of ANEC operators ("We prove commutativity of ANEC operators in any unitary CFT"); superconvergence sum rules; in planar 4d CFT `(a−c)/c` from single-trace data | unitarity; Regge softness | Derivation |
| Costa–Hansen–Penedones, arXiv:1707.07689 [A]; Chandorkar–Chowdhury–Kundu–Minwalla, arXiv:2102.03122 [A] | Regge-limit bounds on OPE coefficients; Regge growth of contact terms from chaos bound | large N | Derivation |
| Lashkari–McDermott–Van Raamsdonk, arXiv:1308.3716; Faulkner–Guica–Hartman–Myers–Van Raamsdonk, arXiv:1312.7856 [A] | Entanglement first law for all balls ⟺ linearised Einstein equations about AdS | **RT formula assumed**; perturbations of the vacuum | Derivation (conditional on RT) |
| Faulkner et al., arXiv:1705.03026 [A] | Second-order (nonlinear) equations from CFT entanglement for states built by Euclidean sources | large N CFT; ball regions | Derivation |
| Lewkowycz–Parrikar, arXiv:1802.10103 [A] | Shape dependence ⇒ Einstein equations; uses "the assumption of equality between bulk and boundary modular flows" | as quoted | Derivation (conditional) |

**Summary for (b).** In d = 3 the stress-tensor three-point function has two parameters beyond `c_T`
(`t_4`, parity-odd `α_T`), rigorously-in-the-physics-sense confined to a disc. That a large higher-spin gap forces
`(t_4, α_T) → 0` (Einstein) is a derivation in d = 4 (AJHKT, BHM, KKSZ); I did **not** find a d = 3 statement
with explicit bounds, and by 2311.04271 the full stress-tensor HPPS statement is open. The correspondence
"`t_4`, `α_T` ↔ parity-even/odd Riemann³ couplings in 4D" is standard lore [M]. Einstein equations from
entanglement are conditional on RT/JLMS and perturbative.

---

## (c) Nonperturbative data in AdS4/CFT3

| Ref. | Content | Status |
|---|---|---|
| Aharony–Bergman–Jafferis–Maldacena, arXiv:0806.1218 [A] | N = 6 U(N)_k × U(N)_{−k} Chern–Simons-matter theory dual to M-theory on AdS4 × S⁷/Z_k | Conjectured duality |
| Kapustin–Willett–Yaakov, arXiv:0909.4559 [A] | Localization: path integral "reduces to a non-Gaussian matrix model" | Derivation (path integral not rigorously defined) |
| Drukker–Mariño–Putrov, arXiv:1007.3837 [A]; Herzog–Klebanov–Pufu–Tesileanu, arXiv:1011.5487 [A]; Fuji–Hirano–Moriyama, arXiv:1106.4631 [A] | `N^{3/2}` free energy from the matrix model; all-genus sum → Airy | Derivation |
| Mariño–Putrov, arXiv:1110.4066 [R, parts] | Fermi gas: `Z(N) = (1/N!) Σ_σ (−1)^σ ∫ dᴺx Π_i ρ(x_i, x_σ(i))`, `ρ(x_1,x_2) = (1/2πk)(2cosh(x_1/2))^{−1/2}(2cosh(x_2/2))^{−1/2}(2cosh((x_1−x_2)/2k))^{−1}` (eqs. 3.3, 3.5); `F(N) ≈ −(π√(2k)/3) N^{3/2}` (4.27); `Z^p ∝ Ai[C^{−1/3}(N − k/24 − 1/(3k))]`, `C = 2/(π²k)` (2.14–2.15). Their Conjecture 1.1 extends Airy to other CSM theories | The Fermi-gas *identity* is exact algebra (Cauchy determinant); the asymptotics are a **derivation** (semiclassical, thermodynamic limit) |
| Hatsuda–Moriyama–Okuyama, arXiv:1207.4283 [A]; Putrov–Yamazaki, arXiv:1207.5066 [A] | Exact `Z(N)` at small N via TBA-type relations | Exact computations |
| Hatsuda–Mariño–Moriyama–Okuyama, arXiv:1306.1734 [A]; Codesido–Grassi–Mariño, arXiv:1409.1799 [A]; Hatsuda, arXiv:1503.07883 [A] | Membrane instantons from refined topological string on local P¹×P¹; closed forms for N = 8; spectral zeta function with conjectured non-perturbative poles | **Conjecture** supported by high-precision numerics |
| Bhattacharyya–Grassi–Mariño–Sen, arXiv:1210.6057 [A] | Universal log N term matches one-loop 11D supergravity | Derivation (one-loop test) |
| Dabholkar–Drukker–Gomes, arXiv:1406.0505 [A] | Supergravity localization reproduces the Airy function | Derivation with assumptions about the supergravity path integral |
| Bobev–Charles–Hristov–Reys, arXiv:2006.09390 [A]; Hristov, arXiv:2204.02992 [A]; Bobev–Hong–Reys, arXiv:2210.09318 [A] | Higher-derivative 4d supergravity; conjectured all-order Airy formula for squashed sphere with masses | **Conjecture** + tests |
| Kubo–Nosaka–Pang, arXiv:2411.07334 [A] | Mass-deformed ABJM on squashed sphere, k = 1: Fermi-gas form and all-order perturbative expansion | Derivation |
| Agmon–Chester–Pufu, arXiv:1711.07343; Chester–Pufu–Yin, arXiv:1804.00949; Binder–Chester–Pufu, arXiv:1808.10554; Binder–Chester–Pufu–Wang, arXiv:1906.07195; Alday–Chester–Raj, arXiv:2107.10274, 2207.11138; Chester–Dempsey–Pufu, arXiv:2412.14094 [A] | Stress-tensor-multiplet correlator of ABJ(M): localization constraints + bootstrap + inversion; M-theory higher-derivative terms | Derivation + numerics; relies on supersymmetry |
| Gautason–van Muiden–Kurlyand, arXiv:2304.12340 [A]; Beccaria–Giombi–Tseytlin, arXiv:2307.14112 [A] | Instanton factors from quantised strings / M2-branes | Derivation |
| Boldis–Korchemsky–Testa, arXiv:2512.02119 [A], arXiv:2602.10196 [R, parts] | Exact functional relations in the Fermi gas; "analytic derivations of several relations for the free energy that were previously known only conjecturally". Verbatim (2602.10196): "Despite this impressive amount of results, their derivation remains largely conjectural." and "The grand potential J(μ,k) has been determined at large μ and finite k by first computing it in the small-k limit and then extending the result to arbitrary k through a set of conjectures supported by high-precision numerical analysis." | Derivation; still uses an ansatz for the instanton structure |
| Hristov–Yetkin, arXiv:2611.07186 [A] | Instantons from CY4 disk partition functions | Proposal |
| Laptev–Schimmer–Takhtajan, arXiv:1510.00045 [A] | **Theorem**: Weyl law and Riesz-mean asymptotics for functional-difference operators `U+U⁻¹+V+ζV⁻¹` (mirror curves), implying trace-class inverses | Rigorous spectral theory for operators closely related to (not identical with) the ABJM Fermi-gas operator |

**Summary for (c).** There is no nonperturbative *definition* of ABJM as a mathematical QFT (no construction
satisfying Wightman/OS axioms) [my assessment; no source claims one]. The sphere partition function is a
well-defined finite-dimensional integral; its `N^{3/2}` law, Airy form, and instanton series are derivations
or conjectures. I found **no** theorem-level proof even of the leading `N^{3/2}` asymptotics of this integral
(arXiv searches "ABJM" with "rigorous / Weyl law / mathematically" returned nothing relevant).

---

## (d) Operator-algebraic emergence

| Ref. | Content | Status |
|---|---|---|
| Rehren, arXiv:hep-th/9905179 [A: title]; Dütsch–Rehren, arXiv:math-ph/0209035 [A] | Algebraic holography; explicit relation between the Klein–Gordon field on AdS and generalized free fields on the boundary | Theorem-level AQFT for free fields |
| Leutheusser–Liu, arXiv:2110.05497, 2112.12156 [A] | Emergent type III₁ algebra of single-trace operators above the Hawking–Page temperature at N = ∞; emergent times from half-sided modular inclusions | Derivation (N = ∞; uses bulk free-field input) |
| Witten, arXiv:2112.12828 [A] | 1/N corrections turn the type III₁ algebra into type II_∞ via the crossed product by the modular group | Derivation |
| Chandrasekaran–Longo–Penington–Witten, arXiv:2206.10780 [A: title]; Chandrasekaran–Penington–Witten, arXiv:2209.10454 [A] | de Sitter algebra; microcanonical type II_∞ algebra whose entropy is generalized entropy of the bifurcation surface — "a derivation of a special case of the QES prescription" | Derivation |
| Leutheusser–Liu, arXiv:2212.13266 [A] | Subregion–subalgebra duality | Proposal (N = ∞) |
| Penington–Witten, arXiv:2301.07257 [A: title]; Kolchmeyer, arXiv:2303.04701 [A: title] | Algebras in JT gravity (type II_∞) | Theorem-level in JT (2D), not AdS4 |
| Jensen–Sorce–Speranza, arXiv:2306.01837 [A]; Kudler-Flam–Leutheusser–Satishchandran, arXiv:2309.15897 [A: title]; Faulkner–Speranza, arXiv:2405.00847 [A] | Type II algebras for general subregions as `G_N → 0`; GSL from crossed products | Derivation ("we argue") |
| Furuya–Lashkari–Moosa–Ouseph, arXiv:2305.16028 [A]; Ouseph et al., arXiv:2310.13736 [A] | Mixing (decay of two-point functions) and emergent type III₁ factors; local Poincaré algebra from modular chaos | Theorems for GFF-type algebras under stated conditions (details not re-read) |
| Gesteau–Santilli, arXiv:2402.10262 [A]; Gesteau–Liu, arXiv:2408.12642 [A] | Explicit large-N type III₁ algebras from matrix models, only above the Hagedorn temperature; stringy horizons | Derivation / explicit models |
| Faulkner, arXiv:2008.04810 [A]; Faulkner–Li, arXiv:2211.12439 [R, parts] | Holographic map as conditional expectation; asymptotically isometric codes. Conjecture (verbatim): "split property for 𝒞 ⟺ local quantum field theory", with the caveat "this is a conjecture that could be wrong" | Theorems for exact codes; **Conjecture** for the split criterion |
| Herderschee–Kudler-Flam, arXiv:2510.01556 [R, parts] | Free string-field model: "the resulting Hagedorn spectrum leads to the breakdown of the split property … for regions that are within a string length of each other"; sections on Minkowski, AdS, AdS black holes. Verbatim: "we conjecture that the split property is restored for regions of O(β_H) distance from the horizon." | Derivation (AQFT techniques); black-hole case **Conjecture** |
| Buchholz–D'Antoni–Longo, arXiv:math-ph/0603083 [A] | Trace-class condition for `e^{−βL_0}` ⇒ nuclearity and (distal) split in chiral CFT; "further discussions on higher dimensional spacetimes" | **Theorem** (chiral CFT) |
| Morinelli–Tanimoto–Weiner, arXiv:1609.02196 [A] | Split property automatic for diffeomorphism-covariant nets on the circle; fails without it | **Theorem** |
| Strohmaier–Witten, arXiv:2303.16380 [A] | Timelike tube theorem in real-analytic curved spacetimes (proof sketched, "Details have appeared elsewhere") | **Theorem** |
| Bahiru et al., arXiv:2301.08753 [A] | State-dressed, approximately local diff-invariant observables commuting with charges to all orders in 1/N | Derivation (perturbative in 1/N) |
| Leutheusser–Liu, arXiv:2411.04183 [R, parts] | Superadditivity of large-N algebras; conjectured generalized connected wedge theorem (GCWT) for AdS3 | Conjecture — **now resolved**: Zhao, arXiv:2509.23119 (JHEP 2026) [A] "a geometrical or general relativity proof for a slightly modified version of their conjecture"; Lima–Pasterski–Waddell, arXiv:2509.26264 [A] prove the forward direction, "establish a modified version … and identify counter-examples to the stronger conjecture" |
| Harlow–Shao–Sorce–Srivastava, arXiv:2509.03589 [A] | "Disjoint additivity" as the locality axiom | Proposal with examples |
| Engelhardt–Liu, arXiv:2311.04281, 2509.05413 [A] | Algebraic ER = EPR; "a rigorous formulation and proof of causal wedge reconstruction for those general boundary subregions whose single trace algebras support von Neumann algebras at large-N" | Proposal; N = ∞ theorem-style statement (not re-checked) |
| Leutheusser–Liu, arXiv:2508.00056 [A] | Volume as index of a subalgebra inclusion | Proposal |
| Liu, arXiv:2510.07017 [A; text searched for stated conjectures: none flagged] | Lectures | Review |

**Summary for (d).** Everything concrete is at `N = ∞` (generalized free fields = free bulk fields) or to fixed
order in `1/N`. No operator-algebraic statement distinguishes Einstein gravity from other local bulk theories,
and none is a nonperturbative construction. Precise open items: the Faulkner–Li split conjecture; the
Herderschee–Kudler-Flam black-hole split conjecture. The LL connected-wedge conjecture is closed.

---

## (e) Quantum error correction

| Ref. | Content | Status |
|---|---|---|
| Almheiri–Dong–Harlow, arXiv:1411.7041; Jafferis–Lewkowycz–Maldacena–Suh, arXiv:1512.06431; Dong–Harlow–Wall, arXiv:1601.05416 [A: titles] | Bulk locality as QEC; bulk = boundary relative entropy (to leading orders in `G_N`); entanglement wedge reconstruction | JLMS: derivation; DHW: theorem *given* JLMS exactly on a code subspace |
| Harlow, arXiv:1607.03901 [A] | RT formula in any code with complementary recovery; "series of theorems", operator-algebra QEC | **Theorem** (finite-dimensional, exact) |
| Pastawski–Yoshida–Harlow–Preskill, arXiv:1503.06237; Hayden et al., arXiv:1601.01694 [A: titles] | HaPPY codes; random tensor networks | Theorems about toy models |
| Kang–Kolchmeyer, arXiv:1811.05482 [A]; Kamal–Penington, arXiv:1912.02240 [A]; Xu–Zhong, arXiv:2411.06361 [A] | Infinite-dimensional (von Neumann algebra) versions; algebraic RT added "for the cases of type I/II factors" | **Theorems** under exact-recovery hypotheses |
| Cotler et al., arXiv:1704.05839 [A]; Hayden–Penington, arXiv:1807.06041 [A] | Approximate recovery via universal recovery channels; state dependence, alpha-bits | Theorems (QI) + physical modelling |
| Akers–Penington, arXiv:2008.03319 [A: title], 2109.14618 [A] | Corrections near phase transitions; state-specific complementary reconstruction ⟺ a quantum minimal surface prescription | **Theorem** (2109.14618) |
| Akers–Engelhardt–Harlow–Penington–Vardhan, arXiv:2207.06536 [R, parts] | Non-isometric codes with Haar-random unitaries. Verbatim: "it is not currently possible to prove something like (3.17) for k-designs, since they have not been shown to obey any measure concentration result analogous to lemma 3.2. This is a pity, as such a result would immediately resolve many open problems in complexity theory." | Theorems for Haar unitaries; design version open and tied to complexity theory |
| Akers–Levine–Penington–Wildenhain, arXiv:2307.13032 [A] | Max/min entanglement wedges, "which we conjecture to be the bulk region reconstructible from B" | Theorems on wedge properties; reconstruction **Conjecture** |
| Gesteau–Kang, arXiv:2112.12789 [A]; Gesteau, arXiv:2302.01938 [A] | Approximate recovery with nonperturbative errors; large-N RT formulae with a bulk cutoff | Theorems in an algebraic framework |
| Akers–Hernández-Cuenca–Rath, arXiv:2108.07280 [R] | See PROBLEM_CANDIDATES.md §1.2 | Theorems + explicit conjectures |
| Czech–Shuai–Wang, arXiv:2502.12246 [A] | Holographic entropy inequalities as necessary conditions for erasure-correction schemes | Derivation |
| Chemissany et al., arXiv:2504.00096 [A]; Wang–Xu–Liu, arXiv:2606.13559 [A]; Alipour Shahmiri–Sharifian–Vardian, arXiv:2609.07866 [A] | Infinite tensor networks and type II factors; approximate-QEC theory of non-isometric codes; twirled Petz reconstruction beyond leading order | Theorems about models / derivations |

**Summary for (e).** The exact structure theorems are complete in finite dimensions and for type I/II
algebras. Nothing in QEC singles out four bulk dimensions or Einstein dynamics. Precise open statements: the
one-shot reconstruction conjectures; concentration for unitary designs (complexity-theoretic); the AHR
bulk-to-boundary cone conjectures.

---

## (f) Holographic entropy cone and HRT in higher dimensions

| Ref. | Content | Dimension / hypotheses | Status |
|---|---|---|---|
| Bao et al., arXiv:1505.07839 [A] | RT cone; contraction-map method; complete for 2, 3, 4 regions (SSA + MMI) | any d; **time-symmetric slice**; needs only a Riemannian manifold | **Theorem** (RT) |
| Hernández-Cuenca, arXiv:1903.09148 [A]; He–Hubeny–Rota, arXiv:2412.15364 [A] | N = 5 cone complete; N = 6: 208 new SSA-compatible extreme-ray orbits, 6 undecided | RT | Theorem (N = 5); partial (N = 6) |
| Czech et al., arXiv:2309.15145, 2401.13029 [A]; Czech–Wang, arXiv:2209.10547 [INSPIRE title] | Two infinite families of inequalities, shown to be facets; an N = 7 inequality | RT | Theorem |
| Wall, arXiv:1211.3494 [A] | Maximin; **SSA and MMI for HRT** | any d; null curvature condition; existence of maximin surfaces argued | Derivation at GR-theorem level, with technical existence assumptions |
| Rota–Weinberg, arXiv:1712.10004 [A] | Maximin cannot prove inequalities for ≤ 5 regions beyond SSA/MMI | — | Negative result about the method |
| Callan–He–Headrick, arXiv:1204.2309 [A]; Caceres et al., arXiv:1304.3398 [A]; Erdmenger et al., arXiv:1705.04696 [A]; Caginalp, arXiv:1905.00544 [A] | SSA / five-region inequalities in Vaidya-type spacetimes hold when NEC holds, fail otherwise | AdS3-Vaidya (and charged, with strips, in 1304.3398) | **Numerics** |
| Bao–Mezei, arXiv:1811.00019 [A; text searched R] | All static inequalities hold "for large regions and times in covariant spacetimes containing collapsing black branes" where the membrane theory applies | any d; leading extensive part only | Derivation (membrane regime) |
| Czech–Dong, arXiv:1905.03787 [R] | Contraction-map inequalities hold for HRT in 2+1 bulk dimensions. Key identity (15): LHS − RHS = ½ Σ (h_LHS − h_RHS)·CMI(i,j) with each CMI built from single-interval entropies. Their own remark: "our proof is unlikely to extend to higher-dimensional cases in its current form" | 2+1; **valid for simply connected bulk**; uses HRT SSA from maximin (NCC) | Derivation; the multiply-connected extension (their §6) is wrong (next row) |
| Grado-White–Grimaldi–Headrick–Hubeny, arXiv:2407.07165 [R] | Error in Czech–Dong for multiply connected bulk; "argument" for `π_1 = Z`; numerics for vacuum wormholes. Details in PROBLEM_CANDIDATES.md §1.4 | 2+1 | Argument with an unproved step + numerics |
| Headrick–Hubeny, arXiv:2208.10507 [A]; GGHH, arXiv:2502.09894 [A; text searched R] | Covariant bit threads, minimax; graph model if time-sheets are "cooperating", "we also find scenarios where it may fail" | any d | Conditional |
| Bousso–Kaya, arXiv:2404.11644, 2502.03516 [R] | See PROBLEM_CANDIDATES.md §1.1 | — | — |
| Grimaldi–Headrick–Hubeny, arXiv:2508.21823 [A; searched]; GHH–Shteyner, arXiv:2601.09987 [A; searched]; Czech et al., arXiv:2601.09989 [A] | Majorization test on null reductions near light-cone configurations: passed by all (contraction-provable) HEIs | any d; perturbations of light-cone configurations in vacuum | Theorem about inequalities; **evidence**, not a proof, for HRT |
| Zhao, arXiv:2602.04888 [A; searched] | Conditional no-short-cut theorem ⇒ HRT cone = RT cone "in arbitrary dimensions" **if** every pair of HRT surfaces has a non-empty exposed region | NCC, maximin, AdS-hyperbolicity + exposed-region condition (assumed, verified for no class) | Conditional theorem |
| Grimaldi–Headrick, arXiv:2610.33905 [R] | Strings 2026 review. Verbatim: "The validity of higher HEIs for the HRT formula remains an open problem." Open question 4: "Are they valid for non-static states?" | — | Confirms open status (24 Oct 2026) |
| Bao–Furuya–Naskar, arXiv:2403.13283 (Bao–Naskar), 2409.17317, 2506.18086 [A / R] | See PROBLEM_CANDIDATES.md §1.3 | RT only | — |
| Akers–Hernández-Cuenca–Rath, arXiv:2108.07280 [R] | QES entropies violate HEIs in general; bulk HEC ⇒ boundary HEC | — | Theorem |
| Czech–Shuai, arXiv:2501.17222, 2601.02472 [A]; Hubeny–Liu, arXiv:2609.22504 [A]; Hubeny–Rota, arXiv:2512.24490, 2512.18702 [A]; Bao–Ferko, arXiv:2608.25287 [A]; Ju–Zhao, arXiv:2606.15173 [A] | Structure of RT inequalities (nesting, RG reading, tripartite form, tree graph models, RG flows, multi-entropy) | RT / light-cone | Theorems or conjectures about the static cone |

**What is proved for HRT in bulk dimension ≥ 4:** SA; SSA and MMI (Wall, under NCC + maximin existence);
all inequalities in the late-time membrane regime (Bao–Mezei, derivation); invariance under light-cone
deformations to first non-trivial order (majorization results). **Open:** every higher inequality for general
time-dependent states, in every dimension ≥ 4, *including pure AdS4 with regions on a generic
non-time-symmetric boundary Cauchy slice*. I found no paper treating strip (or any symmetric) families of
regions in higher-dimensional time-dependent spacetimes beyond numerics for SSA (full-text searches of
1905.03787, 2407.07165, 2502.09894, 2508.21823, 2601.09987, 2602.04888, 2610.33905 for "strip", "planar",
"higher-dimension"; INSPIRE list of the 43 papers citing Czech–Dong).

---

## (g) 2025–2026 items, including AI-assisted

| Item | Source | Relevance | Status |
|---|---|---|---|
| RL search for graph realisations of entropy vectors (N = 3 rediscovery of MMI, application to N = 6) | He–Lee–Ooguri, arXiv:2601.19979 [A] | static cone | Machine-assisted numerics |
| Fine-tuned small language models proving Shannon-type inequalities | arXiv:2606.05729 [A] | classical entropy inequalities, not holographic | ML tooling |
| LLM benchmarks / agent frameworks in hep-th | arXiv:2610.33134, 2610.16733, 2610.14059, 2609.35958, 2609.00107 [titles only] | none claims a result on AdS4/CFT3 by title | Not read |
| Weak-coupling Hawking–Page order; spinless modular bootstrap slope; OpenAI `math` release; 4D gauge-theory existence announcement | As recorded in `/home/user/MAJOR/LITERATURE.md` §1 | 4D CFT / 2D CFT / mathematical GR | **Not re-verified by me**; the first two are not arXiv identifiers I could query |
| Bootstrapping ABJM (Boldis–Korchemsky–Testa) | arXiv:2512.02119, 2602.10196 [A/R] | (c) | Derivation |
| GCWT proof and counterexamples | arXiv:2509.23119, 2509.26264 [A] | (d) | Closes the LL conjecture |
| Stringy split property | arXiv:2510.01556 [R, parts] | (d) | Derivation |
| Majorization / combinatorics of HEIs; Zhao's conditional graph model | arXiv:2508.21823, 2601.09987, 2601.09989, 2602.04888 | (f) | Evidence / conditional |
| Rigorous large-spin results | arXiv:2412.06907, 2505.02897 (2D), 2611.11594 [A] | (a) | Theorems; none gives the anomalous-dimension rate in d ≥ 3 |
| Twirled Petz reconstruction beyond large N; AQEC theory of non-isometric codes | arXiv:2609.07866, 2606.13559 [A] | (e) | Derivations |

I found **no** AI-assisted claim on large-N CFT3, holographic bootstrap bounds, or the covariant entropy cone
(arXiv hep-th abstract searches for "AI-assisted", "language model", "autonomous", "reinforcement learning";
one general web search). This is a negative search result, not a proof of absence.

---

## What is mathematically unproved (one-page synthesis)

1. **Existence.** No holographic CFT3 (ABJM included) is constructed as a mathematical QFT. All "exact"
   data are properties of finite-dimensional integrals obtained by localization.
2. **Bulk locality.** Scalar-correlator locality bounds are physics-level with numerically checked
   positivity; stress-tensor HPPS is open; nothing is non-asymptotic in `Δ_gap`.
3. **Einstein gravity.** `(t_4, α_T) → 0` with the gap in d = 3 has no quantitative proof; Einstein equations
   from entanglement assume RT/JLMS and are perturbative.
4. **Gravitational interactions.** Universal large-spin attraction via stress-tensor exchange is a
   derivation; the rigorous results stop at the GFF limit.
5. **Nonperturbative consistency.** Airy/instanton structure of `Z_{S³}`: conjectural beyond leading order;
   QEC and algebraic results are exact only for models or at `N = ∞`.
6. **Entropy cone.** HRT validity of all inequalities beyond MMI: open in bulk dimension ≥ 4 and for
   multiply connected 2+1 spacetimes with matter.
