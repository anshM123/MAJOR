# Emergence of Einstein gravity in AdS4 from CFT3: state of the programme

Status of this file: working document, 16 Nov 2026. Result D (a theorem that local Einstein gravity emerges
from CFT3 data) is **not established**. Detailed source-by-source audit: `ads4/EINSTEIN_OBSTRUCTIONS.md`
(written by an audit agent from the primary sources; its own derivations are marked [D] there and are
unrefereed).

## 1. The target and the two known routes

**Target (P1).** For unitary CFT3s with `c_T → ∞`, large-`N` factorisation and all single-trace operators of
spin `> 2` above `Δ_gap`: `t_4² + α_T² ≤ K (log Δ_gap)^a / Δ_gap^8`-type bounds on the two non-Einstein parameters of
`⟨TTT⟩` (literature: `|ñ_s − 4ñ_f|/c_T, |ñ_odd|/c_T ≲ ln Δ_gap/Δ_gap⁴`, arXiv:1805.07393 (7.7), a Regge estimate; and
the 4D flat-space bound `|ĝ_3|² M⁸ ≤ 24.9 log(M/m_IR) − 27.6`, arXiv:2201.06602 (4.4), said to translate with
`log(M R_AdS)`). **No quantitative d = 3 theorem exists and no paper to Nov 2026 completes it.**

**Route S-A (dispersive sum rules for spinning correlators; arXiv:2311.04271, 2201.06602, 2106.10274).**

| Link | Status |
|---|---|
| Regge boundedness `J_0 ≤ 1`; spinning sum rules, general `d`; bulk-point dictionary (generic `D ≥ 8` basis) | derivation |
| Reduction to `d = 3`/`D = 4`, including the parity-odd structure; the remaining sum rules | **missing** |
| Regge regime for spinning correlators ("future work" in the source) | **missing** |
| Light-side dictionary (graviton pole in AdS, double-traces) | **missing** |
| Positivity of the flat-space graviton functional | numerics (plots) |
| Curved-regime positivity: `Q_1 > 0`, `Q_2 < 0` in `d = 3`; monotonicity of the improvement ratio | **theorem** (`POSITIVITY_LEMMA.md`: Lemma Q refereed; Lemma M under review) |

**Route S-B (commutativity of energy-flux operators; arXiv:1904.05905, 1904.05892).** Derivation in general
`d` with stated loopholes (momentum eigenstates, contact terms, convergence of the `t`-channel sum); the `d = 3`
formulas are not written (the general-`d` ones have `1/(d−3)` poles); gap suppression of the heavy side is
argued "non-rigorously".

## 2. Does the positivity formula extend to gravitons? (answer to the question asked)

* No source gives the functions that appear for external gravitons in the curved regime. By inference from
  the flat-space kernels (arXiv:2201.06602 (3.5), (A.11)), they are `ν`-moments of the associated conical
  functions `P^{−μ}_{−1/2+iν}(cosh β)`, `μ = |h − h′| ∈ {0, 4, 8}` (harmonics of SO(2) weight `μ` on `H²`).
* For the bound on the cubic coupling only the helicity-diagonal sum rules are needed, whose leading
  Regge-regime weight is `65p − 40p² + …` against the `μ = 0` kernel: the **same sign pattern** `a_1 > 0`, `a_2 < 0`
  as in the scalar case. So, *granted the missing spinning Regge formula*, the leading curved-regime
  positivity for gravitons reduces to exactly `Q_1 > 0`, `Q_2 < 0`, i.e. to Lemma Q.
* For `μ = 4, 8` the Mehler–Dirichlet representation gives closed forms with `Q_0, Q_1, Q_2 > 0` (agent's sketch
  plus numerics; note `Q_2` changes sign relative to `μ = 0`). Not written as a theorem; relevance conditional
  on the inference above.
* **The obstruction is therefore not the positivity of harmonic-function moments.** It is: (i) the spinning
  Regge-moment formula in `d = 3` (how the heavy blocks reduce to harmonic functions, with which weights, for
  each polarisation channel, uniformly); (ii) corrections of relative order `1/m` distinguishing the `++` and
  `+−` channels; (iii) the regime `m ~ M` and the threshold region, where blocks and not harmonic functions
  enter; (iv) the light side: the graviton pole with the AdS propagator, the `|c|²` term, double-traces; (v)
  the three-point tensor structures in `d = 3` including parity-odd ones; (vi) positivity of the flat-space
  functional itself, known only numerically.

## 3. A structural sum rule in `d = 3` (derivation; conditional; refereed)

**Referee verdict (hostile brief, automated, same model family):** "I could not break (★). It is correct as a
conditional derivation, but its novelty and significance are overstated in places" — corrections applied
below. **How to describe it (referee's sentence):** "In a CFT3, granted commutativity of energy detectors on
stress-tensor momentum states, the fourth angular harmonic of the energy–energy correlator in a helicity-2
stress-tensor state is at least `(t_4² + α_T²)/16` and is carried entirely by helicity-6 components of spin `≥ 6`
operators (single-trace ones in the planar limit); this is a `d = 3` specialisation of known superconvergence
relations, gives no gap dependence, and does not single out Einstein gravity."

**Two facts that limit its significance.**
* **`t_4 = α_T = 0` does not diagnose Einstein gravity.** In every CFT3 with at least `N = 1` supersymmetry the
  supercurrent three-point function has a single structure (arXiv:1503.04961), so `c = 0` — including the free
  scalar-plus-Majorana multiplet and ABJM theory at any coupling. There (★) reads `0 + 1 = 1`. The vanishing
  of the non-Einstein parameters of `⟨TTT⟩` is necessary for an Einstein dual and very far from sufficient.
* **Prior work.** Meçaj–Moult–Walters–Xin, "Energy correlator conformal blocks and positivity"
  (arXiv:2512.09986, JHEP 09 (2026) 055): scalar source, blocks for spin-`J` exchange including the `d = 3`
  parity-odd case, coefficients as sums of squares, `0 ≤ F^{(ℓ)} ≤ 1`; spinning sources and superconvergence
  listed as future work (also its sequel arXiv:2608.19322). Belin–Borak–Henriksson–Piron–Zhiboedov
  (arXiv:2608.26246, §7.1.1): polarisation-averaged stress-tensor source, general `d`, parity-even,
  `C_J = Σ_O C_J^O ≥ 0`, and `C_4^T = 2d t_4²/((d+1)²(d+3))`, which is `t_4²/16` at `d = 3`, with a double zero at
  double-twist dimensions; polarisation-resolved positivity left to future work. **So items 1, 2, 5 of
  Proposition S and the stress-tensor term in the fourth multipole are not new.** New as far as the referee
  could find: the helicity-resolved relation `a_n = a_{−n}`, i.e. (★) itself, the spin-`≥ 6` selection, and the
  inclusion of `α_T` — "a short specialisation of [Kologlu et al.], not a new method". The first audit had
  missed arXiv:2608.26246.

In three dimensions the little group of a timelike momentum is SO(2), the energy-flux operators are labelled
by one angle, and the graviton has two helicities. This makes the commutativity sum rules elementary.

**Hypotheses (E)** — standard in the physics literature, none a theorem. Unitary CFT3, unique stress tensor.
Energy-flux operators `E(θ)`, `θ ∈ S¹`, on the improper states of momentum `q = (E,0,0)`, with:
(E1) `E(θ) ≥ 0`, `U(R_φ)E(θ)U(R_φ)^{−1} = E(θ+φ)`, `[E(θ), P^μ] = 0`;
(E2) `∫ E(θ) dθ = H`, `∫ e^{∓iθ} E(θ) dθ = P^x ∓ iP^y`;
(E3) `[E(θ_1), E(θ_2)] = 0` in matrix elements between stress-tensor states for `θ_1 ≠ θ_2`. (Coincident angles
need no separate hypothesis: by (E1), (E2) and Cauchy–Schwarz for a positive form-valued measure,
`|⟨φ|E_n|ψ⟩|² ≤ ⟨φ|H|φ⟩⟨ψ|H|ψ⟩`, so `0 ≤ a_n ≤ 1` without commutativity; an antisymmetric contact term would add an
odd polynomial in `n` to `a_n − a_{−n}`, which is bounded, hence zero. Symmetric contact terms, as in free
theories, are harmless. — Referee.) In the literature (E3) is derived for `J_0 < 3`, with no restriction
excluding `d = 3`, but only on a dense set of states that "does not include important states like momentum
eigenstates"; on stress-tensor momentum eigenstates it is a physics-level hypothesis;
(E4) products `E(θ_1)E(θ_2)` in such matrix elements are computed by inserting the states `|O, m; q⟩`
(primaries `O` in `T × T`, SO(2) charge `m ∈ {−J_O, …, J_O}`; convention `U(R_φ)|m⟩ = e^{imφ}|m⟩`). This is exact: at
fixed timelike `q` a spin-`J` primary above the unitarity bound gives each charge once with positive norm,
descendants are multiples (`P^μ = q^μ`), multi-trace states are multi-trace primaries, conserved spin-`s`
currents keep only `m = ±s`; convergence of `Σ_O|â|²` is Parseval once `a_n ≤ 1`;
(E5) the one-point function in graviton states is `⟨T,h′|E(θ)|T,h⟩ ∝ M_{h′h}(θ)`, `M = [[1, c̄ e^{4iθ}],[c e^{−4iθ}, 1]]`
in the helicity basis `h = ±2`, with `c = (t_4 + iα_T)/4` (the `d = 3` conformal-collider result, arXiv:1707.03007;
phase convention-dependent).

Let `E_n := ∫ e^{−inθ} E(θ) dθ` (so `E_n† = E_{−n}` and `E_n` raises the SO(2) charge by `n`), `|+⟩` the helicity-2
graviton state, and `a_n := ‖E_n|+⟩‖² / ‖E_0|+⟩‖²`.

**Proposition S.** Under (E):
1. `a_n = Σ_O |â^O_{2+n,2}|²`, a sum of squares of the normalised amplitudes `⟨O, 2+n|E_n|T,2⟩`, to which only
   primaries of spin `J_O ≥ |2+n|` contribute.
2. `a_0 = 1` and `a_{±1} = 0`.
3. `a_n = a_{−n}` for every `n`.
4. `a_{−4} = |c|² + Σ_{O ≠ T} |â^O_{−2,2}|²` and `a_4 = Σ_{J_O ≥ 6} |â^O_{6,2}|²`. Hence

   **(★)** `(t_4² + α_T²)/16 + Σ_{O≠T} |â^O_{−2,2}|² = Σ_{J_O ≥ 6} |â^O_{6,2}|²`.
5. `(a_n)` are the Fourier coefficients of a positive measure of mass one on the circle (the
   energy–energy correlator in the helicity state), so `0 ≤ a_n ≤ 1`.

**Consequences.** (a) `t_4² + α_T² ≤ 16 Σ_{J_O≥6} |â^O_{6,2}|² ≤ 16`: the collider disc, and, more sharply, **the
helicity-6 components of operators of spin at least six in `T × T` carry at least the weight
`(t_4² + α_T²)/16`** (fourth harmonic of the detector pair). The mere existence of such operators is automatic
at finite `c_T` (double-twist `[TT]` operators of all spins always exist); only the weight statement has
content. Operators of spin four, and scalars
(which have only `m = 0`), cannot compensate: they enter (★) on the left, with the same sign as the
non-Einstein term, or not at all. (b) No assumption about light scalars is needed. (c) *Planar statement.* The squared
amplitude has a double zero at double-twist dimensions, so `[TT]_{n,ℓ}` drop out of both sides at leading
order: `|c|² + Σ_{ST ≠ T} |â_{−2,2}|² = Σ_{ST, J≥6} |â_{6,2}|²` over single-trace operators. This is meaningful only in
the strictly planar limit: light-ray integrals do not commute with the `1/N` expansion (the planar Regge
intercept is `≤ 2` but the one-loop one is `3`), so no order-by-order statement in `1/c_T` is made.

*Proof.* (1) By (E1) `E_n` maps charge `2` to charge `2+n` and preserves `q`; insert (E4); a spin-`J` primary has
charges `|m| ≤ J`. (2) (E2): `E_0 = H`, and `E_{±1} = P^x ∓ iP^y` annihilate states with `q = (E,0,0)`. (3)
`a_n − a_{−n} = ⟨+|[E_{−n}, E_n]|+⟩/‖E_0|+⟩‖² = 0` by (E3). (4) In the charge `−2` sector the stress tensor
contributes its helicity `−2` state, with amplitude `c` by (E5); charge `6` requires `J ≥ 6`; then use (3) with
`n = 4`. (5) By (E1), (E3) the products `E(f_1)E(f_2)`, `f_i ≥ 0`, are positive, so
`g ↦ ⟨+|∫∫ g(θ_1 − θ_2)E(θ_1)E(θ_2)|+⟩` is a positive functional; its Fourier coefficients are the `a_n`. ∎

**Checks.** Free scalar and free Majorana fermion: `T(q)|0⟩` is a back-to-back two-particle state, the
energy–energy correlator is `δ(θ_12) + δ(θ_12 − π)`, so `a_n = 1` for even `n`; `|c| = 1` (`t_4 = ±4`); the only
primaries are the conserved currents `J_s` with helicities `±s` (and `φ²`, resp. `ψ̄ψ`, with `m = 0`, which do not
enter (★)); (★) reads `1 + 0 = |â^{J_6}_{6,2}|² = 1`. Referee's further checks: explicit two-particle
wavefunctions give `c = +1` (scalar), `c = −1` (Majorana), and for mixtures `c = (n_s − n_f)/(n_s + n_f)` with the
`m = −2` weight splitting as `|c|²` plus `1 − |c|²` on other spin-2 currents; the one-point structure (E5) and
the disc agree in detail with arXiv:1707.03007, including `α_T`; for the 3d Ising model the published
bootstrap window gives `|c|² ≈ 0.93–0.96`, so at least 93% of the fourth harmonic sits in helicity-6 states;
large-`N` Chern–Simons-matter theories have `|c| = 1` at leading order, consistent with a conserved spin-6
current. The spin threshold is consistent with `d = 4` (there a second harmonic exists and starts at spin 4;
`d = 3` has only the fourth).

**What this is and is not.**
* It is the `d = 3`, helicity-resolved form of the "superconvergence" relations of
  Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov and Belin–Hofman–Mathys; see "Prior work" above for what is
  and is not new. It is conditional on (E), which are physics-level assumptions with known loopholes.
* It is homogeneous: the Einstein part cancels, so (★) gives **no rate** in `Δ_gap`. A rate needs the
  spin-two Regge sum rule in which the graviton pole does not cancel, i.e. route S-A.
* The amplitudes `â^O_{m,2}` are not yet expressed through OPE coefficients: the kernels
  `K^{(a)}_m(Δ, J)` with `â^O_{m,2} = c_T^{−1/2} Σ_a λ^{(a)}_{TTO} K^{(a)}_m` have not been computed. Until then (★) is a
  structural identity, not an inequality on conformal data.

## 3a. Can the sum rule give a bound that decays with the gap? (corrected after the kernel computation)

**Correction.** An earlier version of this section concluded, from a scaling lemma, that the relations of
Proposition S "cannot" give a gap-dependent bound, and a referee agreed. **That conclusion was a non
sequitur and is withdrawn.** The lemma itself is true: if `(c, x)` satisfies the relations and positivity, so
does `(tc, tx)` for `0 ≤ t ≤ 1` (the equalities for `n ≠ 0` are homogeneous quadratic, and the measure mixes with the
uniform one). But this only says the allowed set is star-shaped. Its *size* is fixed by the one
inhomogeneous constraint, `a_n ≤ 1`, together with the kernels, which tie each operator's amplitudes in all
harmonics to one or two numbers and depend on its dimension. The correct statement is: with the amplitudes
treated as free there is no gap dependence; with the actual kernels there can be.

**The kernels** (`ads4/KERNELS_D3.md`; agent computation, exact in `Δ` for all conserved `⟨TTO⟩` structures with
`J ≤ 12`, validated on free fields; general-`J` formulas are fits confirmed at `J = 9–12`, not proved; unrefereed).
* Each operator's amplitudes with `m ≤ 0` are fixed by one complex number `a^{++}_O`, and those with `m ≥ 4` by one
  number `a^{+−}_O` (`J ≥ 4`); `â_{m,2} = 0` for `m = 1, 2, 3`. The `m`-dependence is universal:
  `|â_{−k}|²/|â_0|² = W⁰_k = C(k+3,3)² (J+1−k)_{2k}/(Δ−1−k)_{2k}`,
  `|â_{4+k}|²/|â_4|² = W⁴_k = C(k+3,3)² (J−3−k)_k (J+5)_k / ((Δ−5−k)_k (Δ+3)_k)`.
* Commutativity is the family `(S_k)`, `k ≥ 0`:
  `δ_{k2}|c|² + Σ_O |a^{++}_O|² W⁰_k = Σ_{J_O ≥ 4+k} |a^{+−}_O|² W⁴_k = a_{2+k}`; `(S_2)` is (★).
* Scalar exchange: `|â_{0,2}|² = (λ_e² K_e² + λ_o² K_o²)/C_T` with explicit Gamma-function kernels, equal to the
  general-`d` formulas of Kologlu et al. and Córdova–Maldacena–Turiaci at `d = 3`; zeros at `Δ = 6+2n` (even) and
  `7+2n` (odd).
* Every amplitude has a simple zero at the `N = ∞` double-trace dimensions, so each double-trace contributes
  at first order in `1/c_T` to both sides (uniformity over the towers not shown).
* Free-field checks through the full pipeline: `c = ±1`; `|â^{J_s}_{±s,2}|² = 1`; the energy–energy correlator is
  `δ(θ) + δ(θ − π)`.

**What the constraints give.**
* A single relation is not enough: `|c|² ≤ a_2 · sup_{J ≥ 6} W⁴_2`, and `W⁴_2 → 100 (J/Δ)⁴`, which is not small for
  operators with `J` comparable to `Δ`. A gap in dimension alone gives nothing from `(S_2)`.
* The full family does appear to bound `|c|`: for heavy operators `W⁰_k` and `W⁴_k` coincide at leading order in
  `1/Δ`, so a defect at `k = 2` can only come from their mismatch. A linear programme over spins `≤ 100` finds
  maxima of `|c|²` equal to `0.475, 0.109, 0.031, 0.0054` at `Δ_gap = 10, 15, 20, 30`, about `5·10³/Δ_gap⁴`.
  **This is numerical evidence, not a proof.**
* It cannot do better: there is an exact four-operator solution of all `(S_k)` with positive energy–energy
  correlator and `|c|² ≃ 2700/Δ_gap⁴`. So this route can at best give `t_4² + α_T² ≲ 4·10⁴/Δ_gap⁴`, weaker than the
  `Δ_gap^{−8}` expected from the Regge-limit estimates in the literature.
* Open precise problem (not pursued, by decision): prove `|c|² ≤ C/Δ_gap⁴` from `(S_k)`, `0 ≤ a_n ≤ 1` and the
  explicit `W⁰_k`, `W⁴_k`, with the gap imposed on single-trace operators of spin `≥ 2` other than `T`.

All of this is conditional on hypotheses (E) and on the strictly planar reduction, concerns the three-point
function only, and by the supersymmetry remark above does not diagnose Einstein gravity.

**Corollary (infinite gap; conditional).** Assume (E), and a family of CFT3s with `c_T → ∞` such that (i) `T × T`
contains no single-trace primary of spin `≥ 6`, and (ii) the total double-trace contribution to `a_4` tends to
zero. Then `t_4² + α_T² → 0`. **This conclusion is already in the literature** (arXiv:1805.07393, (7.5)–(7.6),
from the Regge limit; the `d = 3` analogue of Belin–Hofman–Mathys); what the present route adds is only that
light scalars and single-trace operators of spin `≤ 5` need not be excluded. And by the remark above it
does not say the dual is Einstein gravity. Hypothesis (ii) is the standard planar statement (the amplitudes vanish at the double-twist
dimensions) and is **not proved here**. Note what this does *not* use: the mere existence of high-spin
double-twist operators `[TT]`, which is automatic in every CFT3, is not evidence of single-trace
higher-spin states; their contribution is what hypothesis (ii) controls. The corresponding statement in
`d = 4` (`a = c` at infinite gap) is Belin–Hofman–Mathys.

## 3b. Helicity-preserving channels: the scalar-kernel action of the graviton weight (numerical; not a theorem)

By the corrected dictionary (`SPINNING_KERNELS.md` Part B, `REGGE_FROM_BLOCKS.md`), the diagonal components
`|c^{++}|²`, `|c^{+−}|²` carry the scalar kernel at leading order. *Granted a spinning Regge-moment formula of the
scalar form* (not yet derived; in progress), the leading action of a functional with weight `f` on such heavy
states at impact parameter `β` is proportional to

`G_M[f](β) := ∫_0^1 f(p) P_{−1/2+iMp}(cosh β) dp = (√2/π) ∫_0^β f̃(Mt) (cosh β − cosh t)^{−1/2} dt`, `f̃(τ) := ∫_0^1 f(p) cos(pτ) dp`.

For the `B_2` weight of the cubic-coupling functional of arXiv:2201.06602, `f(p) = (1−p)³ p (65 + 155p + 47p²)`
(`ads4/code/uniform_positivity2.py`):

| `M` | flat regime (`Mβ ≲ 5`) | curved regime | minimum on the grid |
|---|---|---|---|
| 5 | `G` within 5% of the flat functional `f̂(Mβ)` | `G = 0.0308, 0.00300` at `β = 4, 8` against `a_1Q_1/M² + a_2Q_2/M³ = 0.0302, 0.00300` | positive, decreasing to the end of the grid |
| 20 | within 0.3% | `0.00497, 0.00185, 0.000187` at `β = 2, 4, 8` against `0.00497, 0.00185, 0.000187` | positive |
| 80 | within 0.02% | `2.924·10^{−4}` at `β = 2` against `2.924·10^{−4}`; flat functional there `9.7·10^{−6}` | positive |

So `G_M[f]` is positive on these grids and interpolates between the flat-space functional and the two-term
moment formula `a_1 Q_1(β)/M² + a_2 Q_2(β)/M³` (`a_1 = 65`, `a_2 = −40`); in the curved regime the latter exceeds the
flat functional by orders of magnitude, which is the mechanism described for scalars in arXiv:2106.10274,
and its sign is Lemma Q.

**What is not proved.** (i) Positivity of `G_M[f](β)` for all `β > 0` and all `M ≥ M_0`. The sufficient pointwise
condition in the `coth` representation fails: `K_{f,M}(t) := ∫_0^1 f(p) coth(πMp) sin(Mpt) dp` is negative near `Mt ≈ 10`
for `M ≥ 20` (`uniform_positivity.py`), so a proof has to use the integrated statement — flat-space positivity
of `f̂` (itself numerical in the source) for `Mβ ≲ B_0`, a uniform expansion of the conical function with
explicit remainder in the transition region `β ~ M^{−1/3}` where the two contributions cross, and Lemma Q
beyond. (ii) That this is the complete normalised diagonal contribution: the `1/m` corrections, one of
which has opposite sign in the two diagonal channels, are not included. (iii) Everything about the light
side.

## 3c. The four-stress-tensor Regge-moment formula: what is derived and what is not (`REGGE_MOMENT_TTTT.md`; agent derivation, **unrefereed**)

**Verified.**
* *Scalar formula reconstructed with exact blocks.* Exact `d = 3` scalar blocks (Casimir recursion, residual
  `10^{−15}`–`10^{−18}`) and the Regge moment of a block computed with no large-`Δ` approximation: ratio to the
  heavy-block formula of arXiv:2106.10274 ((3.20)/(2.59)) is `1.0385, 1.0073, 1.00072` at `m = 17.6, 36.7, 113`; with
  the paper's `1/m²` and `1/m⁴` terms, `1 + 4·10^{−6}` and `1 − 10^{−8}` at `m = 113`. A misprint was found in the paper's
  printed (3.27) (`(Δ_φ − d − 1)_3` should be `(2Δ_φ − d − 1)_3`, as in its ancillary file).
* *Heavy limit of conserved three-point structures* (`⟨JφO⟩`, `⟨TφO⟩`, their shadows, `⟨TTO⟩`) in the symmetric
  frame: a plane wave times a constant tensor, transverse at the transverse-space point of the Regge pair,
  with helicity `±ℓ_i` only and corrections `O(1/m)`. This derives the helicity-transfer dictionary.
* *Tensor Fourier–Laplace lemma:* the transform of a spinning kernel is the scalar one times explicit
  rational functions of `ν` times the harmonic of the helicity transfer (checked by direct integration to 9
  digits for currents; stress-tensor rows come from the same code and were not separately integrated).

**Answer to the key question.** In the helicity-preserving components the `β`-dependence is exactly
`P_{−1/2+iν}(cosh β)`, with no helicity-dependent factor of `β` at leading order. There is a `ν`-dependent factor,
but it is a property of the functional (the spinning analogue of the factor divided out in the scalar
case) and is identical for `++` and `+−` for parity-even kernels. Kernels antisymmetric under `1 ↔ 3` at both
pairs give the combination `|c^{++}|² − |c^{+−}|²`.

**Proposed normalised formula (conditional).**
`C^±_{k,ν}[Δ,J] = (2 P_{−1/2+iν}(η_AdS)/m^{2k}) (|c^{++}|² ± |c^{+−}|²) (1 + O(ν²/m², 1/m²))`.
The kernel shape and `ν`-factors are derived and verified piecewise; **the overall normalisation and the
helicity-diagonal pairing are inferred**, not derived.

**First subleading terms.** The `O(1/m)` terms of the two vertices cancel in the `1 ↔ 3`-symmetric part (verified
pointwise for currents, not for `⟨TTO⟩`), so corrections start at `1/m²`. Their coefficients for stress tensors
are **not derived**. The helicity-dependent part is known only in the flat corner,
`+(1/m²)[2λ² (ν/β) J_1(νβ) + λ ν² J_0(νβ)]` with `λ = 4` for `+−` (verified for flat partial waves); its completion at
`β ~ 1` is inferred.

**Blocking steps (agent's list).** (1) No end-to-end test on a spinning block (needs scalar blocks with
unequal external dimensions plus the spinning differential operators, and the spinning shadow
normalisation). (2) The spinning functionals of arXiv:2311.04271 have not been expanded in Regge moments,
and their `d = 3` reduction is absent; which functional uplifts the flat `B_2^{(1)}`, `B_3^{(1)}` is inferred.
(3) Normalisation relative to mean-field theory in the helicity basis. (4) Second-order expansion of the
conserved `⟨TTO⟩` structures. (5) Uniformity in `β`, the threshold `m ~ M`, swapping, discrete series, the light
side.

So the formula is **not established as a theorem**. What is established is the shape of the kernel in the
diagonal channels, on which the positivity statement of `UNIFORM_POSITIVITY.md` (under review) can act.

## 4. Next steps on the direct route

1. Compute the kernels `K^{(a)}_m(Δ, J)` (three-point structures of `⟨TTO⟩` in `d = 3`, light transform, Fourier
   transform), check `|â^{J_s}_{±s,2}|² = 1` in the free theories, and decide whether the spin-`≥ 6` weight is
   suppressed at large `Δ` relative to the others.
2. The spinning Regge-moment formula in `d = 3` (obstruction (i) of §2).
3. Uniform (non-asymptotic in `M`) curved-regime positivity for the explicit functionals.

## 5. Supporting results (independent; none implies Result D)

`POSITIVITY_LEMMA.md` (Lemma Q, Lemma M); `ABJM_THEOREM.md` (free energy, deformation, `c_T` and the ratio
`F/c_T`); `LIGHTCONE_THEOREM.md`; the AdS3 theorem in the repository root.
