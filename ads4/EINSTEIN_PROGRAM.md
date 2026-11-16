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

## 3. A structural sum rule in `d = 3` (derivation; conditional; **not yet refereed**)

In three dimensions the little group of a timelike momentum is SO(2), the energy-flux operators are labelled
by one angle, and the graviton has two helicities. This makes the commutativity sum rules elementary.

**Hypotheses (E)** — standard in the physics literature, none a theorem. Unitary CFT3, unique stress tensor.
Energy-flux operators `E(θ)`, `θ ∈ S¹`, on the improper states of momentum `q = (E,0,0)`, with:
(E1) `E(θ) ≥ 0`, `U(R_φ)E(θ)U(R_φ)^{−1} = E(θ+φ)`, `[E(θ), P^μ] = 0`;
(E2) `∫ E(θ) dθ = H`, `∫ e^{∓iθ} E(θ) dθ = P^x ∓ iP^y`;
(E3) `[E(θ_1), E(θ_2)] = 0` in matrix elements between stress-tensor states, as distributions on the torus
including coincident angles;
(E4) products `E(θ_1)E(θ_2)` in such matrix elements are computed by inserting the states `|O, m; q⟩`
(primaries `O` in `T × T`, SO(2) charge `m ∈ {−J_O, …, J_O}`), the sums converging as distributions;
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

**Consequences.** (a) `t_4² + α_T² ≤ 16 Σ_{J_O≥6} |â^O_{6,2}|² ≤ 16`: the collider disc, and, more sharply, **a
non-Einstein stress-tensor three-point function requires primaries of spin at least six in `T × T`, carrying
at least the weight `(t_4² + α_T²)/16` in the sixth angular harmonic.** Operators of spin four, and scalars
(which have only `m = 0`), cannot compensate: they enter (★) on the left, with the same sign as the
non-Einstein term, or not at all. (b) No assumption about light scalars is needed. (c) In a large-`N` theory
in which the double-trace contributions to the right side are `O(1/c_T)` (not shown here), `t_4² + α_T²` is
bounded by the single-trace spin-`≥ 6` weight.

*Proof.* (1) By (E1) `E_n` maps charge `2` to charge `2+n` and preserves `q`; insert (E4); a spin-`J` primary has
charges `|m| ≤ J`. (2) (E2): `E_0 = H`, and `E_{±1} = P^x ∓ iP^y` annihilate states with `q = (E,0,0)`. (3)
`a_n − a_{−n} = ⟨+|[E_{−n}, E_n]|+⟩/‖E_0|+⟩‖² = 0` by (E3). (4) In the charge `−2` sector the stress tensor
contributes its helicity `−2` state, with amplitude `c` by (E5); charge `6` requires `J ≥ 6`; then use (3) with
`n = 4`. (5) By (E1), (E3) the products `E(f_1)E(f_2)`, `f_i ≥ 0`, are positive, so
`g ↦ ⟨+|∫∫ g(θ_1 − θ_2)E(θ_1)E(θ_2)|+⟩` is a positive functional; its Fourier coefficients are the `a_n`. ∎

**Checks.** Free scalar and free Majorana fermion: `T(q)|0⟩` is a back-to-back two-particle state, the
energy–energy correlator is `δ(θ_12) + δ(θ_12 − π)`, so `a_n = 1` for even `n`; `|c| = 1` (`t_4 = ±4`); the only
primaries are the conserved currents `J_s` with helicities `±s`; (★) reads `1 + 0 = |â^{J_6}_{6,2}|² = 1`.
Einstein gravity (`c = 0`) is the unique point at which the graviton state needs no fourth harmonic from
its own helicity flip.

**What this is and is not.**
* It is the `d = 3` form of the "superconvergence" relations of Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov and
  Belin–Hofman–Mathys. According to the audit, the `d = 3` case is not written in the literature and the
  spin-six threshold with the parity-odd structure included is not stated; this has **not** been checked by
  a referee, and the statement is conditional on (E), which are physics-level assumptions with known
  loopholes.
* It is homogeneous: the Einstein part cancels, so (★) gives **no rate** in `Δ_gap`. A rate needs the
  spin-two Regge sum rule in which the graviton pole does not cancel, i.e. route S-A.
* The amplitudes `â^O_{m,2}` are not yet expressed through OPE coefficients: the kernels
  `K^{(a)}_m(Δ, J)` with `â^O_{m,2} = c_T^{−1/2} Σ_a λ^{(a)}_{TTO} K^{(a)}_m` have not been computed. Until then (★) is a
  structural identity, not an inequality on conformal data.

## 4. Next steps on the direct route

1. Compute the kernels `K^{(a)}_m(Δ, J)` (three-point structures of `⟨TTO⟩` in `d = 3`, light transform, Fourier
   transform), check `|â^{J_s}_{±s,2}|² = 1` in the free theories, and decide whether the spin-`≥ 6` weight is
   suppressed at large `Δ` relative to the others.
2. The spinning Regge-moment formula in `d = 3` (obstruction (i) of §2).
3. Uniform (non-asymptotic in `M`) curved-regime positivity for the explicit functionals.

## 5. Supporting results (independent; none implies Result D)

`POSITIVITY_LEMMA.md` (Lemma Q, Lemma M); `ABJM_THEOREM.md` (free energy, deformation, `c_T` and the ratio
`F/c_T`); `LIGHTCONE_THEOREM.md`; the AdS3 theorem in the repository root.
