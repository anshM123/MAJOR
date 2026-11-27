# m1_remainder: the light remainder R of the cubic bound (referee point M1)

27 Nov 2026. Unrefereed. Nothing outside `ads4/code/m1_remainder/` was edited, and nothing was committed.
`referee_v2/REPORT.md` does not exist in the tree. M1 was taken from the task statement: Lemma 2 of `code/no_eft/NOTE.md`,
the Supplement's "Order of limits" paragraph (supplement.tex l. 693–710), its Scope (l. 2233) and letter_v2.tex l. 136–140.

Grades: [P] proved (given the stated inputs), [D] derived, [DV] derived and checked numerically, [NE] numerical evidence,
[E] estimate (exponents computed, O(1) constants not tracked), [A] assumption.

## Verdict (short)

1. **The term-wise number is essentially right, and the threshold does not move.** We used physical tree-level Einstein data
   (γ ~ κ/(c_T J), GFF weights, exact detector kernels, the twist-2 N-structure). The 10^{691} per-state actions sit on
   twist-2 double twists at J ≈ 0.6Λ ≈ 600. Their physical weight there is about c_T|a|² ≈ 10^{-19}/c_T. That is small only
   polynomially (∝ J^{-8}), so it removes about 20 of the 690 orders.
   * Full spin sum, every J, by the exact generating function: |R₂| = c_T²|R| ≈ 10^{679}–10^{682}, with A_n = κ_n = 1.
     The central value is 5×10^{679}. The physical profile J^{−8} sits near an accidental zero crossing; generic
     profiles give 10^{681.5}.
   * Term-wise sum of |·|: 10^{687}.
   * Threshold: c_T* ≈ 10^{680}–10^{683}, against 10^{690} before. [E]
2. **There is no cancellation that saves the bound.**
   * *Identity [P given (fp-H) at light states, DV].* The spin sum of the light actions is exactly a weighted integral of
     the generating function H_b(z) = Σ_J ĥ_b(J) z^J.
   * *Smooth in J: cancels.* For spin profiles that are smooth over all J, the e^{πΛ/2} cancels exactly. H_b is modest on
     (0,1): H₀(0.995) = 4×10^{10}.
   * *Even J only: does not cancel.* Bose symmetry populates only even J on both ladders. An even-J profile also samples
     H_b on (−1,0], the back-to-back detector point, where H_b(−1) = +4.0×10^{687}.
   * *Size.* Summing over spin buys only 10^{5–7} relative to the term-wise sum. The sign of the result depends on the
     spin profile (H_b(−e^{−t}) changes sign).
   * *Heavy side.* The heavy side cannot compensate. Heavy ladders have no back-to-back component, since the edge is
     suppressed by U_{b,n} ≤ r^n.
3. **Functional modification.** The e^{πΛ/2} is the target heavy image continued to imaginary impact parameter iπ/2. So it
   is fixed by how far in ν the weight reaches. The size scales as 10^{0.681Λ} between Δ_g = 300 and 1000 (π/(2 ln 10) =
   0.682). Capping it means reaching only Λ′ ≈ (2/π) ln c_T, which costs a factor (Λ/Λ′)⁴ in the bound.
   * The cubic-type bound beats the angular bound τ_* = 0.0178 only for Δ′ ≳ 40, which needs c_T ≳ 10^{34}.
   * It would also need a new positivity certificate at Δ′. That job is easier than at 1000: the cancellation scale is
     10^{0.68Δ′}.
4. **The sum rule cannot hold with R this large [D given the framework].** Take any theory with Einstein-like light data
   and c_T below about 10^{681}.
   * If R₂ > 0, then R > 16 L_A n₀ G. This contradicts Σ_heavy ≥ 0, which is rigorous by Theorem FP.
   * If R₂ < 0, the heavy side would have to carry about 10^{681} G/c_T. No heavy state can do that.
   * So either such theories (e.g. ABJM-type M-theory duals with Δ_gap ~ 10³ and c_T ~ 10^{20–40}) violate an input, or
     the light-side bookkeeping fails at finite c_T. The bookkeeping is the continuation of (fp-target) to light states
     [D-sketch], or H_K applied block-wise to the light tower.
   * Either way, "10^{690}" is not a physical threshold. It is a symptom.
5. **The angular route τ_* has none of this.** Its remainder is 𝓛_g ≤ 70 𝓔_light, with 𝓔_light ≈ 48/c_T [E]. So τ ≤ τ_*
   holds once c_T ≳ 2×10⁸ [E].

## Ledger

| step | status | file / output |
|---|---|---|
| 0 | read NOTE.md (Lemmas 2–5, §4), supplement (Order of limits, Scope, sec_fullpos), a2prime_einstein README/scripts. REPORT.md absent | — |
| 1 | per-state ĥ_b(J) for every J ≤ 1000, twists 2, 3, 4, 6, both ladders, at Δ_g = 1000 and 300 (ball arithmetic, certified moments) | `m1a_hhat_allJ.py` → `out/m1a_*.json` [NE] |
| 2 | exact spin generating function H_b(z) of the twist-2 actions: polynomial θ-operators on Ψ_ν(z) = zP_{−½+iν}(X)/(1−z)² plus the I3 Q-part; ν-quadrature. Validated against the moment sums to 1.4×10⁻¹⁶ at z = −0.3 (H ~ 10^{437}) and 4×10⁻¹⁰ at z = 0.1 | `m1b_genfun.py`, `out/m1b_validate_*.txt`, `out/m1b_H_*.json` [D; DV] |
| 3 | physical weights (a2prime normalisation, A_n = κ_n = 1); term-wise and truncated (J ≤ 1000) sums at Δ_g = 1000 and 300 | `m1c_physical_sum.py` → `out/m1c_physical_sum*.txt` [E] |
| 4 | full even-J spin sum from H_b(±e^{−t}), physical power-law profiles; sensitivity to the profile | `m1d_spin_sum.py` → `out/m1d_spin_sum.txt`, `out/m1d_sensitivity.txt` [E] |
| 5 | T-block normalisation n₀ (graviton pole against the weight) | `m1_n0.py` [DV for the number; D-sketch for the identification] |
| 6 | route 3: gap/reach trade-off table | `m1e_gap_tradeoff.py` → `out/m1e_gap_tradeoff.txt` [E] |
| 7 | consistency test of the sum rule with the computed R | below, §5 [D given the framework] |
| 8 | proposed replacement text | below, §7 |

## 1. Setting

* **Remainder and Lemma 2.** R = Σ_{O∈𝓛′} p_O ω[O], with ω[O] = (2πGΛ⁴/9) Σ_b |a_{O,b,0}|² ĥ_b(Δ,J) (Eq. (fp-target),
  continued to light states). Lemma 2 gives τ² ≤ τ₀²(1 + X), with X = π c_T|R|/(512 n₀ L_A) and 16 L_A = 14339.59.
* **n₀.** n₀ = ∫g·16π/(ν²+9/4) dν/(16 L_A) = 4.54×10⁻⁷ = 0.453/Λ² (`m1_n0.py`; Δ_g = 300: 6.64×10⁻⁶). If n₀ = 1 is used
  instead, every log₁₀ c_T* below drops by 6.3.
* **Threshold.** X ≤ 1 needs c_T ≥ c_T* = π|R₂|/(512 n₀ L_A) ≈ 15|R₂| at Δ_g = 1000, with R₂ := c_T² R.

## 2. Per-state structure (m1a) [NE]

* **Size.** Plateau sup |ĥ/T| = 10^{685.7} (Δ_g = 1000) and 10^{209.2} (300). The ratio gives a slope of 0.6807 per unit
  Λ, against π/(2 ln 10) = 0.6822.
* **Twist 2, J ≲ Λ/2.** ĥ_b(J) alternates in sign between consecutive J, with nearly equal magnitudes. For example
  J = 64: +10^{198.6}, J = 65: −10^{200.7}.
* **Twist 2, plateau.** Beyond J ≈ 550 the signs are erratic, both in J and between ladders: ĥ₄/ĥ₀ jumps by ±1 decade
  between neighbouring even J. So no smooth ratio between the ++ and +− amplitudes can cancel them.
* **Higher twist.** Twist 2+2m is down by about J^{−2m} at the ladder edge: twist 4, 6 are 10^{−6}, 10^{−12} at J = 1000.
* **Where the giant values sit.** In S_J = Σ_{k≤J} P_k and in ĥ/T alike, the giant values are a staggered (−1)^J
  structure. This is the ladder truncation edge e^{iJφ} sampled at φ = π.

## 3. Generating function and the cancellation question (m1b, m1d)

**Identity [P given (fp-H) at light states; DV].** At twist 2, U_{b,k} = 1 and r = 1. Then
ĥ_b(J) = T_J Σ_{k≤J−b} a_k + Σ_{k≤J−b} b_k, and

H_b(z) := Σ_J ĥ_b(J) z^J = θ(θ+1)[z^b A(z)/(1−z)] + z^b B(z)/(1−z),   with θ = z d/dz.

* **Ingredients.** A and B are ν-integrals of polynomial θ-operators acting on Ψ_ν(z) = Σ m² f_m z^m =
  z P_{−½+iν}((1+z)/(1−z))/(1−z)². B also contains σ_bΛ²(1−z)^{−4}∫ψ₃(p) P_{−½+iΛp} dp.
* **Coefficients.** Every coefficient polynomial carries the factor m² that this rewriting needs (checked symbolically).
* **Checks.** The per-ν identity holds to all printed digits. The integrated identity matches the moments to 1.4×10⁻¹⁶.

For any spin profile W(J) = ∫ρ^J dμ(ρ) with μ ≥ 0:

Σ_{J even} W(J) ĥ_b(J) = ½ ∫ [H_b(ρ) + H_b(−ρ)] dμ(ρ).

**Values of H_b at z = ±e^{−t}** (`out/m1b_H_*.json`):

| t | H₀(−e^{−t}) | H₀(+e^{−t}) |
|---|---|---|
| 0 | +4.04×10^{687} | — |
| 0.005 | +1.57×10^{686} | 4.2×10^{10} |
| 0.01 | +3.5×10^{684} | 2.0×10⁹ |
| 0.0125 | +1.46×10^{683} | |
| 0.01375 | −9.7×10^{682} | |
| 0.015 | −1.27×10^{683} | |
| 0.02 | −2.1×10^{682} | 9.2×10⁷ |
| 0.03 | +9.1×10^{679} | |
| 0.05 | +6.6×10^{676} | 1.7×10⁶ |
| 0.5 | +1.0×10^{580} | 1.5×10² |

H₄ is within a few percent of H₀ on the negative side.

* **The z > 0 half is modest [DV].** It is the leading-Regge, heavy-like side, and it carries no e^{πΛ/2}. So for a
  profile smooth over all integers J the giant per-state values cancel exactly.
* **The z < 0 half is ~e^{πΛ/2} [DV].** Bose symmetry for identical T restricts both ladders to even J in a parity-even
  theory: ++ and parity-odd on b = 0, +− on b = 4 (closed_forms.py). The odd-J b = 4 structure is parity-odd. This
  restriction always brings in the back-to-back half.
* **Why H_b(−e^{−t}) changes sign.** H_b(−e^{−t}) ≈ ∫ G(ν) e^{ν(π/2 − t/2)} dν. It samples the weight χ at decreasing ν
  as t grows, and χ changes sign near the top of its support.

**Full spin sums with the physical profiles** (`out/m1d_*`; 22 t-nodes; c_T|a|² fitted to C(J+1)^{−p} with p ≈ 8.0):

* **Size.** The back-to-back part is S_pp = 4.0×10^{666}, S_pm = 1.8×10^{666}, S_N = 7.1×10^{665} (linear
  interpolation), so R₂ = 4.5×10^{679}.
* **Generic profiles.** They give |S| ≈ 2–4×10^{668} per family: p = 6, 7 give +, p = 9, 10 give −
  (`out/m1d_sensitivity.txt`).
* **Sign.** The physical p ≈ 8 sits at the zero crossing, which costs about 1.5–2 decades. The sign there is not
  resolved: the spline and linear interpolations disagree.
* **Comparison with truncated sums.** Term-wise |·| over J ≤ 1000 is 10^{687.1}. The truncated signed sum is
  −1.3×10^{684}; the sharp cut at J = 1000 is itself an e^{πΛ/2}-exciting profile.
* **Conclusion.** Summing over spin gains 10^{5–7}, not 10^{680}. A cancellation to O(1) would need the tree-level
  data of all families (pp, pm, N, odd, n = 0, 1, …) to conspire to 680 digits. Nothing in the structure produces that.

**Against the heavy side.** Heavy ladders have U_{b,k} ≤ r^k with r < 1. The edge factor U_{b,n} ≈ (τ/J)^{998} kills the
staggered/back-to-back part, so heavy actions stay ≲ 10³ (cf. twist 998–1010, no_eft §4). No heavy contribution exists
at the e^{πΛ/2} scale. The heavy side is used only through Σ_heavy ≥ 0 and is never expanded (CHMRS §1.1.3 respected).

## 4. Honest size and threshold [E]

| quantity (Δ_g = 1000, A_n = κ_n = 1) | value |
|---|---|
| term-wise Σ\|·\|, J ≤ 1000 (7 families) | R₂ = 1.3×10^{687} |
| truncated signed, J ≤ 1000 | −1.3×10^{684} |
| full even-J spin sum, twist-2 families (pp, pm, N), generating function | R₂ ≈ 5×10^{679} at the fitted p ≈ 8; \|R₂\| up to 10^{681.5} for nearby profiles; sign not determined |
| threshold c_T* = π\|R₂\|/(512 n₀ L_A) | 10^{680}–10^{683} (n₀ = 4.5×10⁻⁷); 10^{674}–10^{677} if n₀ = 1 |
| same at Δ_g = 300 (truncated signed / term-wise) | c_T* ≈ 10^{209} / 10^{211} |

* **Unknown constants.** A_n (GFF normalisation on the canonical structures) and κ_n (Einstein large-spin coefficient)
  are O(1)–O(10³). The N-structure coefficient B_J = 1/J is the a2prime estimate. None of these can move 680 orders.
* **Where the dominant states sit.** The dominant states are at J ≈ 0.6Λ (Δ_g = 1000: J ≈ 620–650; Δ_g = 300:
  J ≈ 200–210).

## 5. Consistency of the sum rule with this R [D given the framework]

* **The sum rule.** At every c_T: Σ_heavy + p_Tω_T + R = 0. Theorem FP gives Σ_heavy ≥ 0, and
  p_Tω_T ≥ −16 L_A n₀ G = −0.066/c_T.
* **Case R₂ > 0.** Then R > 16 L_A n₀ G whenever c_T < c_T* ≈ 10^{681}, and the sum rule fails. This is rigorous given the
  framework.
* **Case R₂ < 0.** Then the sum rule needs Σ_heavy ≈ |R| ~ 10^{681} G/c_T. Heavy actions are ≲ 10³ per unit detector
  weight, up to the polynomial (1−r)^{−3} growth at large impact parameter, which the shell envelope tames. So the
  heavy side would have to sit at spins J/m ~ 10^{270}, which is impossible.
* **What follows.** The block-wise light action (L0) cannot be the light side of a theory with Einstein-like light data at
  any physically realised c_T. If it were, the framework would prove c_T ≳ e^{πΔ_g/2} for every gap-1000 theory with an
  Einstein-like light sector, contradicting known holographic constructions. At least one of these fails at finite c_T:
  1. the continuation of (fp-target) to light states [D-sketch];
  2. the interchange H_K applied block-wise to the light double-twist tower (the light-cone convergence inputs LRC);
  3. the applicability of the inputs to such theories.
* **Remarks.**
  * The paper's own scalar-proxy statement ("per-ν image grows like e^{πν/2} above the mass") agrees with (1). That points
    to (2).
  * The flat-space analogue is the partial-wave expansion of the one-loop discontinuity at unphysical |t| ≫ s, which
    diverges (degenerate Lehmann ellipse). There the light side is the arc of the full amplitude, which is modest.
  * The planar statement limsup τ ≤ τ₀ is unaffected if the light action is anything bounded per unit detector weight
    (Theorem A).

**Decisive next check (not done here).** Evaluate the action of Ω^EE on one light block directly, e.g. [TT]_{0,J} at
J = 8 or 20. Do it by the position-space dispersive integral of `sec:HK`/LRC*, not by the ladder formula. Equation
(fp-target) predicts 10^{42} and 10^{80}. A modest direct value would refute the continuation and would likely dissolve
M1.

## 6. Route 3: functional modification (cost estimate) [E]

The reach of the weight sets the light action. Table (`out/m1e_gap_tradeoff.txt`, signed calibration):

| Δ′ | log₁₀ c_T* | τ₀(Δ′) | vs τ_*(1000) = 0.0178 |
|---|---|---|---|
| 30 | 27 | 2.5×10⁻² | worse |
| 40 | 34 | 7.9×10⁻³ | better |
| 50 | 40 | 3.3×10⁻³ | better |
| 100 | 74 | 2.2×10⁻⁴ | better |
| 300 | 209 | 3.1×10⁻⁶ | better |
| 1000 | 685 | 2.9×10⁻⁸ | better |

* **The trade-off.** At given c_T the method gives at best τ ≲ τ₀(Λ′) with Λ′ ≈ (2/π) ln(c_T/10^{5}), i.e. a bound
  ∝ (ln c_T)^{−4}.
* **Cost.** It needs a Theorem-FP certificate at Δ′. Only Δ′ = 1000 exists; 300 is conditional.
* **No cheaper constraint.** A constraint that caps the light action at fixed reach is not cheaper. W(iθ) ≈ Φ(iθ) is
  forced by W ≈ Φ on the real axis, and Φ(iθ) ~ e^{Λθ} comes from the target ψ₂ itself.

## 7. Proposed replacement text (not applied)

**Supplement, "Order of limits", replace from "The size of the light-side remainder is not a power of Δ_g" to "with or
without a bulk effective theory.":**

> The size of the light-side remainder is not a power of Δ_g. Light double twists of twist < Δ_g have per-state actions up
> to e^{πΛ/2}10^{4.5}T_Δ (L0), and their physical tree-level weights decay only like J^{−8}. Summed over spin by an exact
> generating-function identity, the e^{πΛ/2} would cancel for spin profiles that are smooth in J. It does not cancel for the
> even-J profiles that Bose symmetry imposes, which probe the back-to-back detector configuration. The result is
> |c_T²R| ≈ 10^{679}–10^{682} at Δ_g = 1000, with tree-level Einstein data and untracked O(1) constants. So the inequality
> constrains the exact t₄, α_T only for c_T ≳ 10^{681}. A remainder of this size is also incompatible with the sum rule
> itself at smaller c_T. The block-wise light side (L0 with H_K) is therefore not controlled at any physically realised c_T,
> and Eq. (5) is a statement about the limit c_T → ∞ at fixed Δ_g. This holds with or without a bulk effective theory. The
> angular bound, Eq. (4), carries no such constant.

**Supplement, Scope, replace "and any bound on the size of the $1/c_T$ remainder below $c_T\approx10^{690}$ (row L0)" with:**

> and any control of the light remainder of $\Omega^{EE}$ at finite $c_T$: block-wise it is $\sim e^{\pi\Lambda/2}c_T^{-2}$
> ($c_T^*\approx10^{681}$ at $\Delta_g=1000$, row L0), and at smaller $c_T$ it is inconsistent with the sum rule (open)

**Letter, replace "Bounded term by term, their remainder ... it is a statement about the planar limit." with:**

> Even summed over spin with tree-level Einstein data, their remainder is of order $e^{\pi\Lambda/2}/c_T^2\approx10^{681}/c_T^2$,
> so Eq.~\eqref{eq:cubic} is a statement about the limit $c_T\to\infty$ at fixed $\Delta_g$; Eq.~\eqref{eq:angular} has no such
> constant and holds once $c_T\gtrsim10^{8}$.

## Files

* `m1_weight.py`: the published k = 2 continuum weight χθ, R₋, in mpmath and numpy.
* `m1a_hhat_allJ.py`: per-state ĥ_b, every J; `M1_DG=300 M1_K=300` for the 300 functional.
* `m1b_genfun.py`: the generating function.
  * `python3 m1b_genfun.py validate <z>`: validation against the moment sums.
  * `python3 m1b_genfun.py point <z…>`: values; about 6 min per point.
* `m1c_physical_sum.py`: physical-weight sums and thresholds; `M1_DG=300` for the 300 functional.
* `m1d_spin_sum.py`: full spin sums from H_b. Run without arguments for the assembly, or with `sens` for the profile
  sensitivity.
* `m1_n0.py`: n₀.
* `m1e_gap_tradeoff.py`: the route-3 table.
* `out/`: all outputs.

Read-only imports: `../fullpos` (moment balls, Q_k), `../kernels/closed_forms.py`, and the weight and kernel JSON in `prl/data`.
