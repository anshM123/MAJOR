# The first correction beyond leading Regge for `Ω^EE` at fixed `β`, and per-state positivity at large `J` near the gap

*20 Nov 2026, subleading-Regge agent. Unrefereed, not committed. No existing document or code file was edited. New code
(`ads4/code/spinfunc/`): `sr_scalcmp.py`, `sr_sens.py`, `sr_scalrem.py`, `sr_model.py`. Caches (the NaN `sr_scal_1040_40.pkl` was deleted, see §5):
`cache/sr_scalrem.txt`, `cache/sr_model_1000.txt`, `cache/sr_model_300.txt`. Project modules are imported
read-only. The fits of §2 were run from scratch scripts. They are reproducible from the formulas given there and the row loaders
(`ee1000_chain.trim(Q.even_grid(J))`, with `EE1000_DSTATE=1 W2_DMIN=1000` for direct rows).*

Grades: **[P]** proved, **[D]** derived, **[D-sketch]**, **[NE]** numerical evidence, **[M]** model, **[I]** open.

## 0. Summary (read this first)

1. **The task's premise needs a correction.** For per-state positivity at large `J` near the gap, the size of the beyond-Regge
   correction does not matter. What matters is its **non-oscillating content** [D mechanism + NE, §3].
   * An oscillating correction has the form `Γ = (2J+1)²·[c(x)D_h + s(x)·xD_h′/(2J+1)]` with smooth `c, s` of order 1. This covers
     every amplitude or phase correction of the partial wave, including its full resummation in `J²/m`. After the `ψ`-integration it
     moves `λ/N_s` by `≤ 6·10⁻³` (`J = 40`) to `≤ 8·10⁻²` (`J = 300`) per unit `c`, and by `≤ 5·10⁻⁴` per unit `s` or `x²c`.
   * A non-oscillating (smooth in `x`) correction is amplified by a factor `∝ b̂⁵/Λ²`.
   * The tolerance on a uniform smooth correction, at the gap state, for the worse helicity (`+−`):

     | `Δ_g` | `J` | tolerance (`c₂` units) |
     |---|---|---|
     | 1000 | 40 | 46 |
     | 1000 | 48 | 29 |
     | 1000 | 64 | 14 |
     | 1000 | 100 | 5 |
     | 1000 | 150 | 2 |
     | 1000 | 200, 300 | 1 |
     | 300 | 24 | 25 |
     | 300 | 40 | 8 |
     | 300 | 64 | 3 |
     | 300 | 100 | 2 |
     | 300 | 150 | 1 |

     Restricted to `ν < 50`, the tolerance is `5000 → 80` (1000) and `264 → 11` (300).
   * This is the precise form of "the correction weight grows like `J³`": only smooth or small-`ν` structure carries weight.
2. **First correction beyond leading Regge, exact in `β`: the oscillating part** [D for the `x → 0` limit; NE/M for `x > 0`, §1–§2].
   * The exact leading-Regge action `2P_{−½+iν}(cosh β)` fixes the β-dependence of the partial wave's phase and amplitude as `x → 0`.
     Since `sinh(β/2) = ℓ/m` and `cosh(β/2) = E/m` exactly (`ℓ = J+½`, `E = Δ−3/2`, `m² = E² − ℓ²`), the phase is
     `νβ = ℓθ·κ(β)` with `κ = (β/2)/sinh(β/2)` and `θ = 2 arcsin x`, and the amplitude is `A = √(β/sinh β)`.
   * The Regge-resummed uniform action is
     **`T_h = env_h(x)·A·D_h(θκ) + 2[P_l(cosh β) − A·J₀(νβ)]`**, plus the analogous `k = 3` term with `uσ_h env₃`.
     It is exact at leading Regge for every `β`, equals the flat partial wave at `β → 0`, and resums the phase in `J²/m`. Its first
     correction is an `O(β²x²)` change of `θ_eff` and `A`.
   * Validation against the rows `Γ = m²(G/(a env) − D_h)` [NE]:
     * On the `K₂` rows at `J = 20–36` (`x ≤ 0.6`) and `J = 40` (`x ≤ 0.45`), the κ-model removes 50–80 % of `Γ`. The residual rms is
       43–360 against row rms 160–740, and much of it sits at `β ≥ 0.15`.
     * A 4-parameter fit, `θ_eff = θ[1 − β²(1 + 1.52x² − 7.29x⁴)/24]` and `A·(1 + β²(1.48x² + 2.15x⁴))`, brings the residual to
       10–275. It fits the 24-node direct `J = 40` sub-row (`x = 0.2–0.8`) to rms 179 against a row rms of 1163.
     * It does **not** fit the direct `J = 48` sub-row (residual 1061 against 1273).
     * So agreement within the row errors holds only for `x ≲ 0.45` at `J ≤ 40`. At `x ≳ 0.5` the phase correction is not derived (the
       fitted `g(x)` vanishes near `x ≈ 0.7`).
   * The hybrid model's term `2(P − J₀)` carries `J₀`'s phase `2ℓx` instead of `ℓθ`. At `J = 40` it is out of phase with the data at
     `x ≳ 0.35` (e.g. `x = 0.5`: data −1436, `2m²(P−J₀)/env` +1462, κ-model −1051).
3. **Positivity with the uniform model** (`sr_model.py`, hybrid states, 512 `p`-nodes) [NE + M]:
   * Replacing the hybrid by the κ-model, or by the 4-parameter fit, changes `λ/N_s` by `< 10⁻³` on every state checked: `J = 40–300`
     at `Δ_g = 1000` and `J = 24–150` at `Δ_g = 300`, with `m` from the gap to `2Λ`.
   * Min `λ/N_s` per `J` at 1000:

     | `J` | min `λ/N_s` | at `m/Λ` |
     |---|---|---|
     | 40 | +0.991 | 2.0 |
     | 48 | +0.996 | 2.0 |
     | 64 | +0.993 | 1.062 (gap) |
     | 100 | +0.981 | 1.095 (gap) |
     | 150 | +0.956 | 1.140 (gap) |
     | 200 | +0.931 | 1.183 (gap) |
     | 300 | +0.902 | 1.265 (gap) |

     At `Δ_g = 300`: `J = 24` +0.961, `J = 40` +0.983, `J = 64` +0.968, `J = 100` +0.940, `J = 150` +0.922.
   * **`n = 1` holds in this model on the whole large-`J` near-gap region**, with an oscillating-truncation error `≤ 10⁻³ N_s`. The
     next-order term, estimated as the κ → fit difference, is `< 10⁻³`.
4. **What this does not settle: the non-oscillating part** **[I]**, which decides `n` (§4).
   * Evidence that it is small:
     * **Exact scalar skeleton** (Mellin, all `e1` states, `J ≤ 28`, `β ≤ 0.29`, `Δ_g = 100/300/1000`): the integrated remainder
       beyond the hybrid is `|rem/N_s| ≤ 0.0074` and **decreases with `J`** (`Δ_g = 100`: 0.0074 → 0.0025 from `J = 8` to 28) [NE,
       exact].
     * The `EE` `K₂` chain at `J = 32–40`: integrated corrections `≤ 0.1 N_s`, i.e. a smooth-equivalent level of about 4–11 `c₂` [NE].
     * Saddle structure: around the scattering-crystal saddles every order, including the `J²/m` resummation, multiplies `e^{±iℓθ_eff}`.
       Smooth terms can only come from other critical points: the fold at `ν* ≈ 11` and the light-cone edge `x → 1` [D-sketch].
   * Evidence against: the direct and `K₂` `J = 48` rows carry a smooth offset of `≈ +200–500 c₂`, giving `+2…+9 N_s`. It is positive,
     so harmless if real, but its existence is unexplained.
   * **The `J ≳ 100` verdict needs the smooth part to be `≲ 1–5 c₂`**, below anything measured so far.
5. **Side finding: curvature double-counting in the chains** [NE, code read].
   * The rows' `Γ` contains the leading-Regge curvature. Pointwise at `J = 40`, `x = 0.15 / 0.225 / 0.3`: `Γ = −699 / −930 / −1059`
     against `2m²(P−J₀)/env = −758 / −958 / −1121`.
   * The hybrid chains (`ee1000_*`, `w2_*`, `q3`) add `Σψ₂·env·Γ/m²` on top of a base that already contains `cv₂ + cv₃`. The
     exact-skeleton chains do the same: `l9_neargap` subtracts only `ρ₂ − D_hρ₂(0)` (`ρ` excludes the curvature, `e1_states`).
   * Size near the gap: `cv/N_s = 0.004 / 0.006 / 0.008 / 0.0135` at `J = 32 / 36 / 40 / 48` (1000). It is positive, so the reported
     `λ` are too high by at most this much. No verdict changes, but `dF` should subtract `cv` (or `Γ` should subtract the curvature).

**Grade of the headline:** the analytic correction is **[D]** at `x → 0` (Regge-fixed `κ(β)`, `A(β)`) and **[NE + M]** at `x > 0`, validated
on the rows only for `x ≲ 0.45`. Positivity from it is **[NE + M]**, with oscillating-truncation error `< 10⁻³ N_s`, **conditional on
[I]: the non-oscillating beyond-hybrid part being below the tolerance table of item 1.**

## 1. Variables and the exact Regge constraint [D]

* `m² = (Δ−J−2)(Δ+J−1) = E² − ℓ²` with `E = Δ − 3/2`, `ℓ = J + ½`. With `cosh β = 1 + (2J+1)²/(2m²)` this gives exactly
  `sinh(β/2) = ℓ/m` and `cosh(β/2) = E/m`.
* Flat partial wave: `P_J(cos θ)` with `sin(θ/2) = x = ν/m`. Hilb: `P_J(cos θ) ≈ (θ/sin θ)^{½}J₀(ℓθ)`.
* Leading Regge: `P_{−½+iν}(cosh β) ≈ (β/sinh β)^{½}J₀(νβ)`.
* At `x → 0`, `ℓθ → 2ℓx`, while `νβ = 2ℓx·κ(β)` with `κ = (β/2)/sinh(β/2)`. So the Regge limit fixes, exactly in `β`, the phase
  rescaling `κ` and the amplitude `√(β/sinh β)` of the partial wave at small angle.
* The uniform action `T_h` of §0.2 adds back `2[P_l − A J₀(νβ)]`, which keeps the subleading-in-`1/ν` Bessel terms of `P_l` exactly.
* At fixed `x` the first correction is the `O(β²x²)` change of `θ_eff/θ` and `A`. Its `J²/m` growth in the linearised (`K₂`) expansion is
  the phase: `δφ ≈ ℓθβ²/24 ~ J³x/m²`. That is why `K₂` breaks down at `J²/m ≳ 1.5`, while the resummed form does not.
* Tested and rejected [NE]: an "AdS Breit-frame" phase `Φ = 2m∫₀^x arcsinh(s/√(1−u²))du`. It is worse than κ at large `x`; the data want
  a phase shift that decreases at large `x`.

## 2. Validation against the rows [NE]

Rows: `EE` `K₂` (`C2`, `++`) at `J = 20 (400/800/1600)`, `24 (336/500/1600)`, `28 (336/500/1600)`, `32 (336/1040/1600)`,
`36 (1036)`, `40 (1040)`; direct sub-rows `(1040, 40)` (`x = 0.2–0.8`) and `(1048, 48)` (`x = 0.1–0.7`). A free per-state multiple
of `m²D_h` is projected out (the `a_h` convention).

| model | total rel. rms | `J = 40` `K₂` | `J = 40` direct | `J = 48` direct | `J = 20` (800) | `J = 32` (1040) | `J = 32` (336, `β = 0.2`) |
|---|---|---|---|---|---|---|---|
| rms of `Γ` itself | — | 739 | 1163 | 1273 | 190 | 488 | 667 |
| κ only (no free parameter) | 0.149 | 125 | 1764 | 1712 | 52 | 142 | 359 |
| 4-parameter `g(x)`, `a(x)` | 0.074 | 28 | 179 | 1061 | 10 | 32 | 274 |

`Γ₊₊ ≈ Γ₊₋ ≈ Γ_{C3H,++}` to about 10–30 % at `x ≤ 0.5` (`J = 40`): the large-`J` correction is kinematic and universal. Helicity
dependence appears at `x ≳ 0.6`, where `d^J₄₄ ≠ P_J`.

Not achieved: agreement within the row errors at `x ≳ 0.5`, and on `J = 48`.

## 3. Sensitivity of `λ` to correction shapes (`sr_sens.py`) [D + NE]

`dF/N_s` per unit `B`, where `Γ = (2J+1)²B(x)`, hybrid states at `Δ_g = 1000`, 512 nodes:

| `J` (gap `m`) | `x^0 D_h` | `x² D_h` | `xD_h′/(2J+1)` | smooth `1` | smooth `x²` |
|---|---|---|---|---|---|
| 40 (1.039Λ) | 6.0·10⁻³ | ≤ 3·10⁻⁵ | ≤ 2.2·10⁻⁴ | 54 (`++`) / 128 (`+−`) | 7 / 32 |
| 100 (1.096Λ) | 3.1·10⁻² | ≤ 8·10⁻⁶ | ≤ 4.6·10⁻⁴ | 3.6·10³ / 7.4·10³ | 4·10² / 1.5·10³ |
| 300 (1.265Λ) | 8.2·10⁻² | ≤ 2·10⁻⁶ | ≤ 4.1·10⁻⁴ | 1.7·10⁵ / 2.7·10⁵ | 1.5·10⁴ / 3.5·10⁴ |

* Mechanism [D]: the integrand oscillates in `p` with frequency `≈ b̂` against smooth `ψ₂`. Only the `p → 0` endpoint contributes, and
  that region is the exact-Regge region.
* A smooth term instead integrates over the whole `p`-range, against `N_s ≈ 1024/b̂³`.
* So the tolerances of §0.1 are in `c₂` units of a non-oscillating `Γ`. Equivalently, a smooth *relative* correction to the action
  must be `≲ 10⁻⁶` at `J = 300`. An `O(1/m)` smooth term would need a coefficient `≲ 1.5·10⁻³`; an `O(1/m²)` one, `≲ 2`.

## 4. The non-oscillating part [I], with the evidence

* **Exact scalar** (`sr_scalrem.py` → `cache/sr_scalrem.txt`): `max|Σψ₂E₂ρ₂|/N_s` per `J` (all `e1` states):

  | `Δ_g` | range | `β` |
  |---|---|---|
  | 100 | 0.0074 (`J = 8`) → 0.0025 (`J = 28`) | ≤ 0.29 |
  | 300 | ≤ 0.0013 | — |
  | 1000 | ≤ 0.0002 | — |

  No growth with `J`. (The decomposition-based "smooth part" printed there is a fit artefact, not used.)
* **`EE` chain at `J = 32–40`, 1000**: `|dF|/N_h ≤ 0.10`, curvature double-count included (≤ 0.008).
* **`J = 48`** (`K₂` not converged; direct row aliased, see EE_DG1000_LARGEJ §3a): smooth `+200…+500 c₂`, giving `+2…+9 N_s`. Unexplained.
* **Small `ν`**: below the fold (`ν* ≈ 11`) the crystal has two thimbles, a genuine source of non-oscillating structure. Tolerance at
  `ν < 50`: `≥ 80 c₂` (1000) and `≥ 11 c₂` (300, `J = 150`). The measured direct values at 300 (`ν = 10–17`, `J = 24, 28`) are
  `+90…+440 c₂`, i.e. positive.

## 5. Verdict on `n = 1` over the large-`J` near-gap region

| `Δ_g` | region | uniform model | grade |
|---|---|---|---|
| 1000 | `J = 40–300`, `m` from gap to `2Λ` | `n = 1`, min `+0.902` (`J = 300`, gap), oscillating truncation `< 10⁻³` | [NE + M] |
| 300 | `J = 24–150` | `n = 1`, min `+0.922` (`J = 150`, gap) | [NE + M] |

Both are conditional on the non-oscillating beyond-hybrid part being below §0.1's table **[I]**. Nowhere does `n = 1` fail in the
model. The earlier `J = 24, 28` failures at 300 came from small-`ν` extrapolation (EE_DG300_SMALLNU) and are not addressed by this
model beyond the tolerance column.

Not done:
* The exact scalar skeleton at `J = 40` (`sr_scalcmp.py 1040 40`) returned NaN: the `exactmom` Chebyshev Mack evaluation overflows or
  NaNs at `J = 40`, and the `BlockMellin` set-up takes more than 5 min. A fixed exact scalar at large `J` would directly test the
  smooth part.
* A derivation of the `x > 0` phase function.

## Hand-over

* **Reproduce:**

  | command | output |
  |---|---|
  | `python3 sr_model.py 1000 40,48,64,100,150,200,300 512` | `cache/sr_model_1000.txt` (~15 min) |
  | `python3 sr_model.py 300 24,28,32,40,48,64,100,150 512` | `cache/sr_model_300.txt` |
  | `python3 sr_sens.py 1000 40,64,100,150,200,300 512` | shape sensitivities |
  | `python3 sr_scalrem.py` | `cache/sr_scalrem.txt` |

* **Next steps (value order):**
  1. Bound the non-oscillating part: either an exact scalar skeleton at `J = 40–100` (fix the `exactmom` NaN; the scalar is the clean
     laboratory), or a direct RQMC with a *smooth-projected* estimator (integrate `Γ` against `ψ₂·env` directly instead of
     sampling nodes), which avoids aliasing.
  2. Explain the `J = 48` `+500` offset.
  3. Remove the curvature double count in the chains (subtract `cv₂ + cv₃` from `dF`, or subtract `2m²(P−J₀)/env` from `Γ`).
  4. Derive `θ_eff(x, β)` beyond `x → 0`. The fit suggests the β² phase correction vanishes near `x ≈ 0.7`.
