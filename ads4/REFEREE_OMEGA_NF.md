# Referee report: `OMEGA_NF.md` and `ads4/code/omegau/nf_*.py`

Referee: independent and hostile, 18 Nov 2026. Nothing under review was edited. Scratch code is in `ads4/code/referee_omeganf/`:
`bound_check.py`, `envelope_check.py`, `c1eff.py`, `c1scan.py`, `cutdep.py`, `cutdep2.py`, `lratio.py` and `rerun_nf4.py`. Its `cache/` is
git-ignored (the line was added to `.gitignore`).

## Reproduced

| check | result |
|---|---|
| `A`, `B`, `B_AdS`, `k(0)` (mpmath, 60 digits, from the expanded exact `ψ₂` in the JSON, not the factored form) | `3022.33769`, `−5327.20955`, `−6552.66203`, `88865643959/3233230000`: **agree** |
| Bound `240(n⁸16L_AdS(Λ/n)/k)^{1/2}/Λ⁴` (own `L_AdS`) | `Δ_g = 100`: `2.21e-4` (n=1), `1.03e-3` (1.5), `3.66e-3` (2.1). `Δ_g = 1000`: `2.89e-8`, `1.40e-7`, `5.17e-7`. **Agree.** With `Λ → Δ_g` instead of `√((Δ_g−2)(Δ_g−1))`, the bound is 6 % lower at `Δ_g = 100` (same convention as `OMEGA_U`) |
| `nf4` re-run for `Δ_g = 100`, `Δφ = 3`, output to my cache | 3733 states, **bit-identical** to the pickle |
| `nf7` tables and `nf6` thresholds against the claims in §0 | the quoted numbers match the logs: `b̂ ≤ 5/10/20` thresholds, the failure counts by `b̂`, `N_U/N_s`, the I3 shifts, the residual flip and the contact shifts |
| Regge `P_{−½+iν}` used in `nf6` vs `mpmath.legenp` | 2·10⁻¹³ |

## Issues

**Fatal: none** (for correctness). Issue M1 does, however, remove the support for the headline "central `n ≈ 1.5`".

**M1 (major, soundness): the large-`b̂` envelope `C₁uN` with `C₁ = 1` is uncalibrated, and the only available data contradict it.**
* **The envelope is not normalised like the correction it bounds.** `N = Σ_k|∫φ_kP|` inherits NF02's large-`b̂` cancellation. At `m = 1.5Λ`
  the ratio `N/∫|Σ_kφ_k-weight·P|` is 0.17 (`b̂ = 10`), 0.025 (20) and 0.005 (40), and `FR/N_coh` falls to `5·10⁻³` at `b̂ = 40`.
  * `C₁uN` therefore assumes that the `O(1/m²)` correction density has the same `1/b̂³` cancellation as the leading action. Nothing
    argues this.
  * The `(𝒥x²)ⁿ` terms grow with `b̂`.
  * A pointwise relative envelope gives no uniform `n`.
  * For M02 this did not arise, because the forward term dominated `N`.
* **The skeleton's own data, where it mostly works, imply `C₁ ≫ 1`.** Restricted to a `b̂` bin, the skeleton's `n(0.3, 0.03)` and the
  `C₁ = 1` envelope give:

  | `b̂` range | skeleton | envelope, `C₁ = 1` |
  |---|---|---|
  | (10, 12] | 1.25–1.29 | 1.50 |
  | (15, 20] | 2.00–2.10 | 1.22 (equivalent to `C₁ ≈ 4`) |
  | (20, 30] | 3.0 | 1.12 (`C₁ > 5`) |

  * The envelope binds at `b̂ ≈ 10.5`, where the two agree. It is far more lenient wherever the skeleton fails.
  * On fully controlled states with `b̂ ≥ 10`, `|S|/(uN_s)` reaches 0.7 and `E/(uN_s)` reaches 2.6–14 (`Δ_g ≥ 100`).
* **The "artefact" claim is unproven.** No resummed or independent computation is shown. If it holds, there is no evidence for any value
  of `C₁` at all.
* **This is inconsistent with `REFEREE_UNIFORM_REMAINDER`.** That report measured `C₁_eff` = 17 (`Δ_g = 20`), 3.2 (50) and ≈ 0.9 (100),
  in this same large-`J` region. My `c1scan` gives:

  | `C₁` | 3.2 | 10 | 17 |
  |---|---|---|---|
  | `n` | 2.17 | 3.4 | 4.4 |

  `C₁ = 1` is not conservative. Even "conservative" `C₁ = 3` sits below what the skeleton implies at `b̂ > 15`.
* **Verdict on the threshold.** The honest statement is that `n ≈ 2.1` at best is supported for `b̂ ≤ 20`. There is no `n ≤ 3` for
  `b̂ ≤ 40`. The large-`b̂` region is an assumption [I], not [M].

**M2 (major, labels): "combined threshold `n ≈ 1.5`" is stated without restricting `Δ_g`.** At `Δ_g = 20` the `b̂ ≤ 10` skeleton alone
gives 1.80 (κ = 0.3) and 2.15 (κ = 1). With the measured `C₁_eff`, `Δ_g ≤ 50` gives `n ≳ 2.2–4.4`.

**Minor**
* **m1: the leverage comparison is apples-to-apples only as a fraction of each functional's own bound numerator.**
  * The `k = 3` reduction of ≈ 440 at `Δ_g = 100` is ≈ 50 × 9. The factor 50 is `L_NF/L_M02` = 49 (100), 53 (1000), the same `A` that
    costs the factor 7 in the bound. Only the factor 9 is shape.
  * So "because `ψ₃(0) = 0`" over-attributes the reduction. For `k = 4` it does not apply (`φ₄(0) = 0` for both functionals).
  * The "48–200" figure holds at `Δ_g = 100` only. At `Δ_g = 1000` it is 2.6·10³–1.1·10⁴.
* **m2:** the `n = 1` row of the bound table is not a bound, because positivity fails there. The `b̂ ≤ 5` states near `m = Λ` are `++`,
  `J ≤ 2` only. Every `+−` state with `m < 1.8Λ` has `b̂ > 5`.
* **m3:** the I3 shift (0.12N_s at `+−`, `m = Λ`, which is 60 % of the flat margin) is not added to any margin. It is harmless at
  `n ≥ 1.13`.
* **m4:** "`10⁻²–10⁻⁵` of the integrand norm" depends on the norm used. With the coherent norm I find `7·10⁻²–4·10⁻³` for
  `b̂ = 10–40`.
* **m5:** contact shifts are quoted up to `n = 1.5` only. At the conservative `n = 2.1` they are roughly 2× larger, which is still
  negligible.

## Verdicts

* **Correctness:** sound. The numbers reproduce: constants, bounds, pickles and thresholds.
* **Soundness:** claim 2's combined `n` is not supported (M1). Claims 1, 3 and 4 hold as conditional arithmetic.
* **Labels:** I15 should be [I] (uncalibrated), not [M] (M1). "Central 1.5" should read "≥ 2.1, `b̂ > 20` open" (M1, M2). The leverage
  attribution is overstated (m1).
* **Applicability:** this is conditional on I1–I17 and meaningful for `Δ_g ≳ 100` only. The forward-leverage removal is real.
* **Novelty:** low. It is a mechanical re-application of the `OMEGA_U` pipeline.
