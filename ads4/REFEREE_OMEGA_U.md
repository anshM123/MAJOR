# Referee report on `ads4/OMEGA_U.md` and `ads4/code/omegau/`

Referee: independent automated referee, 17 Nov 2026. Scratch code: `ads4/code/referee_omegau/` (`thr_check.py`, `light_check.py`; cache
subdir git-ignored). Nothing under review was edited. Disclosure: importing `omegau/u7_summary.py` re-ran its top level and rewrote the
git-ignored `omegau/cache/u7_summary.txt` from the same pickles; its content was unchanged.

## Reproduced

* Pointwise B′ ↔ stripped identity, both helicities: re-derived with my own sympy (coefficients of `ψ₂` give `2−y`, of `ψ₃` give `±u`). It holds by
  construction of the `φ_k`. The forward values `2aJJ` and `a(2JJ−42)` follow from `CH`.
* `N_U/N_s = 62.007` at `m = 1.01Λ`, `++`, `J = 0` (M02). The skeleton for the states that set the threshold was recomputed from scratch
  (no cache). It agrees with the cache exactly, and the `κ` values are identical.
* **M02 threshold at `Δ_g = 100` is 1.40.** At `1.35Λ` it fails on `++, J = 2` (`κ = 0.69`) and on `++, J = 0` (`κ = 0.9989`, so it misses by
  a hair). At `1.40` the binding state is `++, J = 2`, with `κ = 1.19`, `F/N_s = 0.283`, `E/N_s = 0.010`, `U/N_s = 0.229`.
* `k(0) = 27.485396` (recomputed from the unstripped low sides). The bound `n⁸·16L_AdS(Λ/n)/k(0) = 2010` at `Δ_g = 100`, `n = 1.4`, so
  `√(t₄²+α_T²) ≤ 1.14·10⁻⁴`. The other supported-`n` entries also check (2297, 1958, 2833, 2821). None of them is printed in any log.
* Leverages of `g₄` and `g₆` (k5, fw4, fw5): agree.

## Issues

**Fatal:** none.

**M1 (major, soundness): the "supported threshold" is set by a convention, not by a result.**
* **How `n` depends on the tolerance.** `n(κ_req)` at `Δ_g = 100` (M02, worst over helicities and skeletons) is:

  | `κ_req` | 0.05 | 0.1 | 0.3 | 0.5 | 1 | 2 |
  |---|---|---|---|---|---|---|
  | `n` | 1.11 | 1.15 | 1.27 | 1.35 | 1.40 | 1.45 |

  At `Δ_g = 50` the corresponding values are 1.12 … 1.50.
* **What the truncated series itself says.** On the uncontrolled nodes of the threshold states it gives `|ρ_k|` = 0.4–2 %. At `1.01Λ` it
  gives ≤ 3.5 %, and the signed sum there is `(F+S_all)/N_s = +1.0075`. "Uncontrolled" is a verdict of the `x²`-series ratio test, not an
  indicated O(1) error. The 100 % tolerance is therefore neither derived nor calibrated.
* **Amplification is not the problem at the threshold.** `U` is a per-`k` absolute envelope, so the amplification is already priced in.
  At the threshold `N_U/N_s` is only 1.3–2.5; the factor 62 sits at `1.01Λ`, where the note correctly declines `n = 1`.
* **What the criterion does not cover.** `E` and `U` contain none of the following:
  * the graviton-versus-scalar transfer (I4, I5);
  * the odd-`k` prefactor shift (I3). Interpolated to `1.35–1.4Λ` it is ≈ 1–3 % of `N_s`.

  The spare margin of the binding state is `F−E−U ≈ 0.044N_s`, which is comparable.
* **Consequence for the bound.** Read honestly, the `Δ_g = 100` bound ranges from `5.4·10⁻⁵` (`n = 1.15`) through `7.9·10⁻⁵` (`n = 1.27`) to
  `1.5·10⁻⁴`, not "1.1–1.5·10⁻⁴". The note should present `n(κ_req)` and not a single "supported" `n`.

**M2 (major, labels): "graviton pole only in `C⁺₂`" (I10).**
* **What holds.** If the Einstein light side of each `C_k` with `k ≥ 3` vanishes identically in `ν`, then the `φ_j` weights multiply zero
  and are irrelevant. So `A` and `B_AdS` are unchanged. That part is correct.
* **What is mislabelled.** The applicability of the counting rule to the project's `C⁻₃` is [I], per `LIGHT_INPUTS.md` §0 items 3–4: it
  needs a spin-3 dispersive completion, and the project's `C⁻₃` is not the CLS-D one. For `C⁻₅`, `C⁺₆` and `C⁻₇`, which are never defined
  as two-pair moments, it is weaker still. I10 labels it "[D from Q]".
* **Why it matters.** The leverages are 6.7·10³ … 36 (smeared) and 10¹³ (forward), so any failure is fatal.

**Minor**
* **m1: the `|ĝ₃|²` leverage row is too large by a factor π.** `u5` substitutes `G = 1` but keeps π, while `k(0)` is defined per `πG`.
  * Correct values: `0.043, 0.471, −1.115`, `Σ = 1.63`.
  * So `k_AdS − k(0)` is 17 % / 0.6 % / 0.006 % (`n = 1`), not 54 % / 1.9 % / 0.02 %.
  * Check: for F_g3 the corrected `k = 2` entry × `k(0)` gives 9677/13860, which matches `ADS4_LIGHT.md`.
* **m2: the residue estimates are quoted at `n = 1`, but the bound is used at `n = 1.4–1.6`.**
  * Smearing at `Λ/n` multiplies the relative residues by `n²`.
  * `|a| ≈ 36` per polynomial underestimates the next-to-top coefficient of a degree-`d` polynomial, which is ~`d·36`.
  * `g_f ≲ G_f log Λ` with a unit constant is an assumption; the flat two-sided bounds carry O(10–100) constants.
  * As a result, `Δ_g = 50` is not clearly "negligible". `Δ_g ≥ 100` survives.
* **m3: the near-gap `κ` statements are incomplete.**
  * "κ = 0.007–0.2 on `[Λ, 1.3Λ)`" contradicts the note's own `Δ_g = 1000` `++` row (0.44–0.98). That jump comes from the absolute
    `1e−4·E` control switch in `remk.residual_fixedJ`, which is an artefact of the cut.
  * Negative `κ` (F_g3 `+−`; M01 at `Δ_g = 20`) is not reported.
* **m4: the `u3` flat-limit check is close to tautological.** The model is flat plus corrections that vanish by construction. The check
  is also non-monotone at `Δ_g = 1000` in `+−`: 3.4·10⁻³ at `1.2Λ, J = 4`; 7.8·10⁻⁴ (more than at `Δ_g = 100`) at `1.5Λ, J = 6`. It is not
  evidence for I1.
* **m5: "exact scalar-skeleton moments" (item 3) is a mislabel.** What is used is the order-1/m⁸ truncated asymptotic series.
* **m6: the uncontrolled region starts lower than stated.** For `k = 7` the uncontrolled nodes start at `x ≈ 0.43–0.45`, not at
  "`x ≳ 0.55–0.6`".

## Verdicts

* **Correctness:** the code and numbers reproduce. There is one π error (m1), in the conservative direction.
* **Soundness:** the threshold, and therefore the bound, depends on the convention (M1).
* **Labels:** I10 and "exact" are overstated (M2, m5).
* **Applicability:** everything is conditional on I1–I14. The dictionary is untested for gravitons. Only `Δ_g ≳ 100` is meaningful.
* **Novelty:** low. This is a mechanical uplift and re-scan of known ingredients.
