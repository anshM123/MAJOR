# Referee report on `ads4/OMEGA_PRIME.md` and `lightinputs/t6_normalisation.py`, `t6b_flip_scan.py`

Referee: independent and hostile. Date: 16 Nov 2026.
* Scratch code: `ads4/code/referee_omegaprime/r1_N_direct.py`. Its log is `cache/r1.log`, which is git-ignored.
* No file under review was modified.
* `t6.log` and `t6b.log` agree with §6.

## Verdicts

* **Correctness: correct, independently reproduced.**
  * `r1` recomputes the Fourier–Laplace responses by direct numerical integration on `H²`, without `laplace.py` or `clsd_kernel.py`. The CLS-D Fourier
    coefficients are obtained by quadrature, not from the Pochhammer formula. The checks are `l = 1, a = 6` and `l = 2, a = 8`, at `ν = 0.7` and `2.3`.
  * Every `A(±,∓)` agrees with the `t6` closed forms, and so does `N(ν)`:
    `l = 2`: `N(0.7) = −0.345605 − 0.246861i`, `N(2.3) = −2.854797 − 0.620608i`, with relative difference ≤ `5·10⁻⁸`.
    `l = 1`: `(2/5)(2ν+i)`, with relative difference ≤ `2·10⁻⁸`.
  * `A(+,−) = −A(−,+)` holds exactly.
  * The flip responses `A(±,±)` are ≤ `10⁻¹⁸` for both kernels.
  * The zeros of `N`, the monotonicity of `|N|` and its value `0.2047` at `ν = 0` are all correct.
  * Flip-freedom is in fact *structural*, which is stronger than "cancels in the combination":
    * at `(1,3)`, `z₁,z₃ = ŷ ± X(φ)` gives the flip tensor `(A² − (X·ε)²)^l`, which is even under `φ → φ+π`, while `sgn sin(−φ)|sin φ|^{−δ}` is odd;
    * at `(2,4)`, `𝔄` is antisymmetric under `z₂ ↔ z₄`, while the flip tensor is symmetric.

  So the flip responses vanish for every `ν`, `a` and `l`.
* **Soundness:**
  * At leading Regge order the logic holds, given its premises.
  * §7's kernel freedom is right in substance. The commutator vanishes as an operator, so it vanishes in every matrix element. The two possible sources
    of contact terms are both controlled:
    * the shared endpoints `x₀`, `x₀⁺` of the two rays, where the collision is exactly what `J₀ < 3` controls;
    * collinear `z₁ ∥ z₃`, where the celestial structure carries `|z₁₃|^{2+δ}`, with `Re = 2.5`, and kills `δ`- and `δ′`-type terms.
  * The covariant extension of `𝔄` is unique, because three points have no cross-ratio. The `(2,4)` M-tables use exactly the frame `(x̂, n_ℓ, n_ℓ^⊥)`.
  * The leading heavy action is *exactly* the old one. Nothing assumed in the heavy formula is broken:
    * it is stated for an arbitrary `K₁₃⊗K₂₄`;
    * the MFT normalisation is a property of the states;
    * the sign and phase are carried by `N`.
  * The count of 25 is right: 13 even + 12 odd.
* **Labels:** overstated (M1, M2).
* **Applicability:**
  * `k(0)` and `A` are unchanged, given (F).
  * The claim that `B_AdS` and `C₁` are unchanged is not established (M3).
* **Novelty:** low. The functional is the KKSDZ/CLS-D bare commutator. What is new: `N(ν)`, the observation about the `(2,4)` freedom, and the bookkeeping.

## Issues

**Fatal:** none.

**M1 (label / internal consistency).** After §7 the functional is no longer shown to be a CLS-D member:
* the lift of `𝔄` lives in a 25-dimensional space;
* CLS-D's duals are light transforms of `⟨TTO_{3/2+iν,3}⟩` structures, which number ≤ 23 (fewer for conserved `T`).

Despite this, §0.1, §1 and §5 row 1 still describe the CLS-D dual and cite CLS-D (200) as [Q] for the flat map.
* *Fix:* the argument of §3 (the odd fixed-`u` space is one-dimensional, plus overlap matching with `N ≠ 0`) is independent of the kernel. Relabel it as [D-sketch/I].
* The [DV]/[D] labels of §0.2, §6 and §7 also omit two premises:
  * the leading-Regge light-ray mode map (`SPINNING_E2E.md` §3.3, [Q+D]);
  * `k_r = 3`.

**M2 (label).** `t6` normalises by the `d₀`-normalised combination `Σw_qΠ̂^q`. `SPINNING_E2E` §3.3 says that the `(k,ν)` normalisation `n(ν)` relative to the true `Ψ` is not tracked.
* Hence "`N` complete" and "`N ≠ 0` on ℝ" need `n(ν)` to be finite and nonzero on ℝ. This is plausible but not shown.
* The same caveat applies to the zeros of `N`, its phase, and the "`T`-pole" remark.

**M3 (applicability).** `C̃ = Ψ̃/N` is complex, and `N(−ν) = −N̄(ν)`, so `C̃` is not manifestly `ν`-even.
* By (F) the light side `Y(ν)/N(ν)` may contain a `ν³` term (`E = −1`). That term would shift `k_AdS` at order `Λ⁻¹`, which gives `log Λ/Λ` in the bound.
  This exceeds the stated `O(log Λ/Λ²)`, so the claim that `B_AdS` and `C₁` are "unchanged" is unproven.
* *Fix:* prove a shadow identity `Ψ̃_{−ν} ∝ Ψ̃_ν`, or compute `Y(ν)`.

**m1.** Specify `Re C̃` in `Ω′`.
**m2.** The analogy "`(2,2,J≥4)`" is inapposite, since `O` has `J = 3`, for which the count is 23. The right reason for 25 is that `L[O]` has spin `−½−iν`.
**m3.** The `l = 1` row is a parallel computation, not a check of `N`.
**m4.** The "exact" closed form passes through floats: `(−1)**q3` with `q3 < 0` gives `±1.0`. This is harmless.
**m5.** "No sign-indefinite terms" means no flip interference. `|c_s|² − |c_o|²` is itself indefinite, as in `Ω`.
