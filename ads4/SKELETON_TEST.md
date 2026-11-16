# Skeleton test: is the exact two-pair spinning moment (exact scalar skeleton) × (spin factor uniform in ν/m)?

Status: written 16 Nov 2026 by a numerics agent. **Unrefereed.** This tests an input of `UNIFORM_REMAINDER.md`: assumptions (b) "transfer of the
scalar remainder" and (c) "spin-specific `1/m²` terms within `C_s cosh β/m²`, `C_s ≤ 40`", plus the unresolved skeleton choice (`Δφ = 3` vs `2`).
The test uses currents (`⟨JJJJ⟩`, `l = 1`, `Δ_J = 2`, both Regge pairs spinning). It is not run on gravitons.

Code: `ads4/code/skeleton/` (`sh run_all.sh`; `sh run_all.sh fast` redoes only the analysis). Caches and outputs are in `ads4/code/skeleton/cache/`
(git-ignored); the main output is `cache/out_s4_summary.txt`. The code imports, unchanged, `e2e/` (exact scalar blocks, grids), `fourpt/`
(`tmoments4`, `predict4`, `t5_analysis_l1.hel_vectors`), `fourpt/tttt/t2_mft_heavy_l1.factor` (exact MFT normalisation `Q_h`) and `uniform/`
(`skel`, the CHMRS order-8 series). It reuses the `fourpt` block-tensor cache (10 cases) and adds 4 cases of its own. Nothing outside
`skeleton/` was modified.

Marks: **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence (read off exact-block numbers); **[I]** inferred or assumed;
**[ND]** not done. Fits are diagnostics only: no constant quoted here is derived.

## 0. Results

1. **The data [NE].**
   * **Spinning side.** Exact two-pair `⟨JJJJ⟩` tensor Regge moments (exact partial waves of `SPINNING_4PT.md`). They are reduced to the covariant
     even (`D_1⊗D_1`) and odd (`D_5⊗D_5`) pair kernels, on the exact helicity waves `W[ĥ_same, ĥ_same]` and `W[ĥ_opp, ĥ_opp]`. "same" is CHL
     `T^{++}` = flat `++` (`J_z = 0`); "opp" is flat `+−` (`J_z = 2`). The `±` combinations follow by linearity, `C^± ∝ |c_same|² ± |c_opp|²`.
   * **Cases.** 14 cases: `β ≈ 0.4` with `m = 12.4, 17.4, 22.4, 32.3, 42.2, 72.0, 101.8` (`J = 2 … 20`), and `β ≈ 1` with
     `m = 10.6, 16.3, 20.1, 29.7, 39.3, 70.0, 100.7` (`J = 5 … 52`). Moments `k = 2, 3, 4, 6`, at `ν = x m`, `x = 0, 0.05, …, 0.65`, `ν ≤ 16`.
   * **Skeleton side.** Exact scalar-block moments with identical external scalars `Δφ ∈ {1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5}`. Both sides are
     MFT-normalised: `Q_h` for the spinning side, the CHMRS `q_{Δ,J}` for the scalar side.
2. **Reach [NE].** Double precision is not the limit at large `ν`. At `m ≈ 12` the moments on the standard and the doubled `u` grid agree to
   `10⁻⁶…10⁻⁸` up to `ν = 12`; near zeros of `P_{−½+iν}`, the agreement is `10⁻³` of the local amplitude. The scalar moments are stable in `N` and in the
   grid, and agree with the CHMRS series to `10⁻⁶` at `m = 42`.

   The limit is the dropped small-`u` region (tail estimate `< 10⁻⁴` required). The largest `x` reached is:

   | `k` | `m ≤ 22` | `m ≈ 30–42` | `m ≈ 70` | `m ≈ 100` |
   |---|---|---|---|---|
   | 6 | 0.6–0.65 | 0.25–0.35 | 0.10–0.15 | 0.05 |
   | 4 | 0.35–0.65 | 0.15–0.25 | — | 0.05 |
   | 3 | ≤ 0.5 | — | — | — |

   At `k = 2` only `x ≤ 0.1` is reached; the `ν = 0` values are reliable for `m ≲ 32`. **`x = 0.6` is reached only for `m ≲ 22`.**
3. **Same helicity (`J_z = 0`, flat `++`): the factorisation holds, and the best skeleton is neither 2 nor 3 [NE].**
   * **Non-uniformity.** Define the spin factor `S = Q_h R_spin / (Nrat R_skel)`. Its ν-non-uniformity `m²d(x)` (§2) is small for every skeleton
     `Δφ ∈ [2, 3]`. It is smallest at `Δφ ≈ 2.5` (`k = 3, 4, 6`, odd kernel):
     * `max|m²d| ≤ 0.5–2.5` for `m ≤ 22`, `x ≤ 0.65`;
     * `≤ 3–7` for `m = 30–72`, `x ≤ 0.15–0.35`.

     The values at the `Δφ = 2`/`3` skeletons are
     * `≤ 4–10` at `β = 0.4`;
     * `≤ 12–41` at `β = 1`, oscillating near zeros of `P`.
   * **Constant offset.** What the skeleton choice mainly changes is a **ν-independent** `1/m²` offset. Within each `β` it is converged in `m` and
     nearly independent of `k`:

     | `m²(S₀−1)` | `Δφ = 2` | `Δφ = 2.5` | `Δφ = 3` |
     |---|---|---|---|
     | `β = 0.4` | `+12…+20` | `0…+3` | `−23…−33` |
     | `β = 1` | `+11…+19` | `−5…−10` | `−37…−55` |

   * **Zero offset.** The offset vanishes at `Δφ* = 2.50–2.55` (`β = 0.4`) and `2.34–2.41` (`β = 1`), for every `k = 2, 3, 4, 6`.
   * **The truth lies between the two `UNIFORM_REMAINDER` candidates.** For currents these are `Δ_J = 2` and `Δ_J + 1 = 3`. The `k = 2` offset
     gap between them, 35, has the size of the `+12.9` vs `−19.1` gap quoted there.
4. **Opposite helicity (`J_z = 2`, flat `+−`): the factorisation fails at `x ≳ 0.2`, for every scalar skeleton [NE].**
   * **The deviation is barely Δφ-dependent.** It grows with `x` *and* with `m` at fixed `x`. At `x = 0.1`, `β = 0.4`, `m²d = −3, −8, −17, −49, −79, −117`
     for `m = 12 … 72`, so it is not a `1/m²` effect.
   * **At `x ≤ 0.1–0.2` it is the curved flat-space helicity term.** This is `η_AdS → η_AdS − c/m²` with `c = 2J_z² = 8` (`SPINNING_E2E.md` §6.3),
     inferred from flat space and **not fitted**. With it, `|m²d(0.1)| ≤ 1.3` (`β = 0.4`, `m ≤ 32`) and `≤ 6.4` up to `m = 70`.
   * **At `x ≥ 0.3` nothing tried describes it.** The residual reaches `|m²d| = 50–600`.
     * Two natural flat densities with the scalar kernel `E_k`, `d̃^J_{22} = P^{(0,4)}_{J−2}` and `d^J_{22}`, overcorrect by 2–4×.
     * **The bulk-point form of this component is not identified [ND].**
   * **The offset.** It is `m²(S₀−1) = −6…−16` (`Δφ = 2`) and `−17…−40` (`Δφ = 2.5`). So the two helicities prefer different skeletons.
5. **Against `C_s ≤ 40` [NE].** Define `C_s,eff(X) = max_{x≤X}|m²(S₀−1) + m²d(x)|/cosh β` (odd kernel; `cache/out_s4_summary.txt`, TABLE 1).

   | | `X = 0.2` | `X = 0.4` | `X = 0.6` |
   |---|---|---|---|
   | same, `Δφ = 2` | 11–19 | 11–18 | 12–27 |
   | same, `Δφ = 2.5` | 0–6 | 1–6 | 1–7 |
   | same, `Δφ = 3` | 25–45 | 26–45 | 26–57 |
   | opp (`c = 0`), `Δφ = 2`, `β = 0.4` | 22–129 | 58–228 | 80–194 |
   | opp (`c = 0`), `Δφ = 2`, `β = 1` | 11–42 | 13–45 | 13–56 |
   | opp (`c = 8`) | 7–50 (outliers 77) | — | — |

   For opp with `c = 8`, the `X ≥ 0.4` values are dominated by points near zeros of `P` (up to `10³`).

   * **Same helicity:** `C_s ≤ 40` holds for `Δφ = 2` and `2.5`, but **not for `Δφ = 3`** at `β = 1` (`C_s,eff` up to 57).
   * **Opposite helicity:** `C_s ≤ 40` is **violated** at `x ≤ 0.2` already, for `m ≳ 17` at `β = 0.4`. The violation grows with `m`, which an
     envelope `C_s cosh β/m²` cannot capture.
   * **The even kernel** gives the same picture. A few points are spoiled by its leading flip leakage near zeros of `P` (§2).
6. **Graviton case (`l = 2`) — inference only [I].**
   * **(i) `++` (`J_z = 0`).** If currents are a guide, the graviton spin factor is ν/m-uniform relative to a scalar skeleton. The skeleton choice
     then mainly moves a constant `1/m²` offset, and the true effective `Δφ` sits strictly between the natural candidates. For currents that is
     `≈ 2.4–2.55`, β-dependent and not a simple rule, so no value for gravitons follows. The `Δφ = 2` vs `3` band of `UNIFORM_REMAINDER.md` is
     then a plausible bracket for `++` but not a proven one. It could equally be `[3, 4]` for `l = 2`.
   * **(ii) `+−` (`J_z = 4`).** The scalar-skeleton × uniform-factor model should fail more strongly. The Regge-regime helicity term scales as
     `2J_z²`: 32 instead of 8.
     * So assumption (b) is **unsupported for the `+−` channel**, and (c) with `C_s ≤ 40` is not met by its photon analogue.
     * The binding near-gap `+−, J = 4–6` states of `UNIFORM_REMAINDER.md` (`κ_crit` 0.29–0.66) rest on exactly this input.
     * What would replace it is the spinning bulk-point (flat helicity density) form of the AdS `+−` moment. For photons the natural candidate
       did not fit (item 4), so this is open.
   * **(iii) Odd `k`.** `k = 3` behaves like `k = 4, 6` in every respect tested (`x ≤ 0.5`). This is weak support for "odd `k` by continuation" at
     the level of these `Π̂` moments.

## 1. Objects

* **Spinning side.** `R_{kern,h}(k,ν) = (−1)^J K_kern[exact] / K_kern[leading]` (`spin.py`).
  * `K_even` and `K_odd` are the covariant pair-kernel actions of `fourpt/tttt/t4_sign_l1.py`.
  * The leading value is `fourpt/predict4` (factorised formula, no fitted constant), `→ 1 + O(1/m²)` (`SPINNING_4PT.md` §5).
  * MFT normalisation: `Q_h` (`SPINNING_4PT.md` §4.6; Weyl-even; `m²(Q_opp − 1) = 66` and `m²(Q_same − 1) = 7` at `β = 0.4`).
  * The odd kernel's leading value is purely `∝ P_{−½+iν}(cosh β)`. The even kernel contains predicted helicity-flip leakage (`≤ 5 %` for `ν ≤ 6`;
    larger near zeros of `P`).
* **Skeleton side.** `R_skel(k,ν;Δφ)` = exact scalar `Π̂_{k,ν}` / leading scalar formula (`SPINNING_E2E.md` §4.1), times
  `Nrat(Δφ) = 1 + O(1/m²)`, the CHMRS MFT normalisation (ν- and k-independent; `skeleton.py`).
  * **Only `Π̂` moments are compared.** The CHMRS `C` functional differs by `Π̂_{k+2}` subtraction terms. Their spinning analogue is not defined
    (`SPINNING_E2E.md` §3.3), so it is not tested [ND].

## 2. Measures (`s4_summary.py`)

* **Offset.** `S₀ = S(x=0)`; `m²(S₀ − 1)` is the ν-independent spin factor at `1/m²`.
* **Non-uniformity.** `d(x) = (Q_h R_spin K_p − S₀ Nrat R_skel K_p)/(A·E_k(x)/2)`.
  * `K_p` is the leading spinning value.
  * `A` is its amplitude envelope: the preserved part `|n₀₀|A_P` plus the flip leakage.
  * `A_P` is the Harish-Chandra amplitude of `P_{−½+iν}(cosh β)`, or `|P|` for `νβ < 1.5`.
  * `E_k` is the bulk-point envelope.
  * So `d` is the relative deviation from "skeleton × constant". It is finite at zeros of `P`, and normalised like the `C_s cosh β/m²·N_hyb`
    envelope of `uniform/t3_scan.py`.
* **Variants.**
  * `c = 8`: `K_p → K_p + n₀₀[P(η − c/m²) − P(η)]`, the curved helicity term.
  * `flat`, `flat2`: subtract `n₀₀(D_h − P_J)/A`, with `D_opp = d̃^J_{22}` or `d^J_{22}` at `cos θ = 1 − 2x²` [I].
* **Exclusions.** Points with a tail estimate `> 10⁻⁴` are excluded.

## 3. Checks

| check | result |
|---|---|
| two-pair moments on the standard vs the doubled-`u` grid, `m = 12.4` (`β = 0.4`), `k = 4, 6`, `ν ≤ 12` | `10⁻⁶ … 10⁻⁸` relative; `≤ 10⁻³` of the local amplitude at `ν = 6` (zero of `P`) |
| scalar moments: standard / doubled / tripled-split `u` grid / doubled `t` grid; `N` vs `4000` (`m = 42`) | `≤ 4·10⁻⁸` up to `ν = 10`; series agreement `10⁻⁶` (`cache/logs/s0*.log`) |
| `m²(S₀ − 1)` converged in `m` | yes, within `≈ 1` from `m = 17` to `100` (`k = 4, 6`); `k = 2` only to `m ≈ 32` |

## 4. What is not established

* No statement here is derived. Everything is [NE] at `β ≈ 0.4, 1`, for currents.
* The large-`x` reach shrinks with `m` (§0 item 2), so the claimed `m`-independence of `m²d` at `x ≥ 0.4` rests on `m ≤ 22`.
* The `+−` bulk-point form is not identified. `k = 2` has no `x > 0.1` data. Gravitons were not computed.
* Whether the CHMRS-`C` subtraction terms change any of this was not tested.

## 6. Follow-up: a "spinning skeleton" for the opposite-helicity component (`s5_spinning_skeleton.py`, `cache/out_s5.txt`)

**Question.** Is the right curved object for the opposite-helicity component a spinning harmonic-analysis kernel rather than a scalar skeleton?
The candidates come from `SPINNING_KERNELS.md` Theorem K and the `REGGE_FROM_BLOCKS.md` weight-shifted Legendre functions.

**Construction (no fit).**
* Exact scalar skeleton (`Δφ = 2` and `2.5`) times the MFT constant `S₀`.
* The preserved part `n₀₀P_{−½+iν}(cosh β)` of the established leading value is replaced by `n₀₀P_{−½}(cosh β)·X_ν(β)/X_0(β)`. This matches the
  leading formula at `x → 0`.
* `m²d_X(x)` is measured exactly as in §2.

The candidates `X`:

| candidate | flat limit |
|---|---|
| `Φ^ν_{1,1}` (Chebyshev–Mehler) | `J₀` |
| `Φ^ν_{2,2}` (index 2 ↔ `J_z = 2`; equals the two-pair kernel `𝒦^ν_{2,−2}` of `fourpt/moments4` in its inversion-referred frame) | `J₀` |
| `Φ^ν_{2,−2}` (`= 𝒦^ν_{2,2}`) | `J₄` |
| `(ν²+¼)P^{−1}_{−½+iν}` | — |
| `(ν²+¼)(ν²+9/4)P^{−2}_{−1/2+iν}` | — |

The control is `X = P` (= §0 item 4).

**Result: nothing works [NE].** Every candidate is worse than `X = P`, for both helicities, both kernels, `k = 4, 6` and both skeletons. Examples
(odd kernel, `k = 6`, `Δφ = 2`, `m²d_X` at `x = 0.1`):

| case | `X = P` | `Φ_{1,1}` | `Φ_{2,2}` | `Φ_{2,−2}`, `P^{−1}`, `P^{−2}` |
|---|---|---|---|---|
| opp, `β = 0.4`, `m = 12 → 72` | `−3 → −117` | `−3 → −328` | `−3 → −1030` | `10³–10⁷` |
| opp, `β = 1`, `m = 10.6 → 70` | `−1 → −6` | `−4 → +1360` | `−23 → +13100` | — |
| same (control), `Δφ = 2.5`, `β = 1`, `m = 70` | `−0.8` | `−1380` | — | — |

**Why [D, elementary].**
* **The diagonal candidates fail at fixed `ν`.** `Φ^ν_{m,m}(β) ≠ P_{−½+iν}(cosh β)` at fixed `β`, and the difference is `m`-independent. Replacing `P`
  therefore spoils the established leading formula at fixed `ν` by `O(1)`, so `m²d_X ∝ m²`. That is what the data show.
* **The off-diagonal and weight-shifted candidates have the wrong flat behaviour.** `Φ_{2,−2}` and `P^{−l}` vanish like a power of `β`, or have flat
  limit `J_{2m}`. They cannot match a component whose leading value is `∝ P` (exact for the odd kernel, `SPINNING_4PT.md` §4.5).
* **This is the corrected dictionary of `SPINNING_KERNELS.md` Part B.** The `Φ`/`P^{−l}` indices are Regge-pair helicity *transfers*. The
  opposite-helicity component measured here is helicity-preserving at both pairs, weight `(0,0)` (flat `+− → +−`), so its kernel is `P`.
  `J_z = 2` is the heavy state's `s`-channel label `h₁ − h₂`, not a Regge-pair weight. The weight-`(2,±2)` structures enter only the flip
  components; there they are already in the leading prediction (the even-kernel leakage `𝒦_{±2,·}` of §1).

**Conclusion.**
* **What does not work.** None of the `m`-independent harmonic-analysis structures replaces the scalar `β`-dependence of the opposite-helicity
  component.
* **What the data still point to** is a correction of relative size `∝ x·m/m²` with coefficient `≈ 2J_z²`. In the Regge regime this is the curved
  helicity shift `η → η − 8/m²` of §0 item 4. At `x ≥ 0.3` it is a flat-level `x`-dependence, still unidentified.
* **What would be needed** is an `m`-dependent spinning bulk-point/Regge interpolation (the spinning analogue of CHMRS's fixed-`J` series). It is
  not available [ND].

## 7. Follow-up 2: scaling of the opposite-helicity non-uniformity, and an inferred graviton projection (`s6_opp_scaling.py`, `cache/out_s6.txt`)

### 7.1 Scaling of the photon data

All fits in this subsection are diagnostic fits [NE]; none of their constants is derived. `d` is the §2 quantity at `c = 0`, i.e. relative to
scalar skeleton × constant.

**(a) `d·m/x` is not `m`-independent.** At `β = 0.4`, `x = 0.1` (odd kernel, `k = 6`, `Δφ = 2`) it is:

| `m` | 12 | 17 | 22 | 32 | 42 | 72 |
|---|---|---|---|---|---|---|
| `d·m/x` | `−2.2` | `−4.7` | `−7.7` | `−15.3` | `−18.7` | `−16.3` |

* It grows at small `m`, where `νβ < 2` and `P′/P` is not yet in its oscillatory regime.
* It saturates near the envelope `8/sinh β = 19.5`.
* At `β = 1` it oscillates in sign (`−6 … +0.6`; envelope 6.8).

So `d ∼ c x/m` holds only as an **envelope**, `|d| ≲ c x/(m sinh β)`. The signed `d` is a phase shift, not a smooth function of `x/m`.

**(b) One-parameter helicity-shift fit, `d = c·g`.** Here `g` is the exact linearised shape of `η → η − c/m²` on the preserved part.

| | `x ≤ 0.2` | `x ≤ 0.3` | residual (rms) |
|---|---|---|---|
| `β = 0.4`, all variants | `c = 5.7–9.5`, median 7.9 | `c = 5.4–9.6`, median 7.8 | `≤ 0.03` (`m ≤ 22`), `≤ 0.2` (`m ≥ 32`) |
| `β = 0.4`, `m ≥ 20` | median 8.3 | median 8.3 | |
| `β = 1`, all variants | median 8.5, range `−2.5 … 15` | median 7.8, range `−0.7 … 18` | `0.04 … 1.0` |

"All variants" means `k = 4, 6`, odd/even kernel, `Δφ = 2, 2.5`.

* At `β = 1` the fit fails at `m = 16` (`c ≈ −1`, residual 0.8–1).
* So `c ≈ 8 = 2J_z²` is supported at `β = 0.4` and only roughly at `β = 1`.

**(c) Smooth ansätze describe the data poorly.** Fitted over all `m` at fixed `β, k`, kernel, for `x ≤ 0.3`:
* `d = c₁x/m + c₂x²/m + c₃x²/m²`: rms residual 0.14–0.42 at `β = 0.4`, 0.5–0.8 at `β = 1`. The coefficients are unstable (`c₂ = ±60 … ±117`).
* `d = c₁x/(m sinh β)`: residual 0.35–0.94.
* `d = c₁x²` (flat-level): residual 0.15–1.0.

**What matters below.** "Decreases like `1/m` at fixed `x`" is true only at **fixed `β`**, since `m sinh β → ∞`. On states of **fixed `J`**,
`m sinh β = (2J+1)cosh(β/2)` stays fixed. The envelope is then `≈ c x/(2J+1)`, which does not decrease with `m`. In that case the effect is
flat-space physics at fixed impact parameter.

### 7.2 Graviton projection — inference only [I]

**Assumptions.**
* `c_g = 2J_z² = 32` for `+−` (`J_z = 4`). This assumes the photon value 8 scales as `2J_z²`, as the flat helicity term does; it was not checked for `l = 2`.
* States: `m = Δ_g m_flat`, `J = 4–6`, `m_flat ∈ [1, 1.5]`, `x ∈ (0, 1)`, `k = 2, 3`.

**Quantities.** Both are relative to the bulk-point envelope `E_k`.
* `d_A = max_x |Z(η − c/m²) − Z(η)|/|Z(η)|·2/E_k(x)`: the whole helicity shift.
* `d_B`: only its curvature part, `ν c|1/sinh β − 1/β|/m²·2/E_k`, i.e. the part beyond its flat limit.

| `Δ_g` | `d_A` (`k = 2`; J = 4 / 6) | `d_B` | `κ_crit(+−)`, M02, `[Λ, 1.5Λ)`, `Δφ = 3 / 2` |
|---|---|---|---|
| 20 | 1.59 / 0.87 | `5·10⁻²` | 0.58 / 0.62 |
| 50 | 1.59 / 0.88 | `9·10⁻³` | 0.66 / 0.53 |
| 100 | 1.59 / 0.88 | `2·10⁻³` | 0.44 / 0.41 |
| 300 | 1.59 / 0.88 | `2.5·10⁻⁴` | 0.32 / 0.32 |
| 1000 | 1.59 / 0.88 | `2·10⁻⁵` | 0.29 / 0.29 |

* `d_A` is the same at `k = 3` (1.47). It is independent of `Δ_g` and of `m_flat` (checked to `m_flat = 5`), and its maximum is at `x ≈ 0.4`.

**Which reading applies depends on the flat `+−` density in the hybrid.**

* **Reading B.** `UNIFORM_REMAINDER.md`'s hybrid uses the flat `+−` density `d̃^J_{44}` in `F_flat`. If that density correctly contains the flat
  part of the helicity shift, only `d_B` is new. Then `d_B ≤ 0.05 ≪ κ_crit` for every `Δ_g`. As far as this effect goes, the `+−` channel would stay
  positive down to `m = Λ` for all `Δ_g = 20 … 1000`. The other uncontrolled nodes behind `κ_crit` are still unbounded, so this is not a positivity claim.
* **Reading A.** If the flat helicity content of the AdS `+−` moment is *not* what `d̃^J_{44}` gives, the uncovered correction is up to `d_A`. It is
  `O(1)` and `Δ_g`-independent: it exceeds `κ_crit` for every `Δ_g`, and also the ≈`0.18N` margin of fully controlled states. With
  `d_A ≈ 0.51, 0.25, 0.12` at `J = 10, 20, 40` (`Δ_g = 100`), the projected correction stays above `κ_crit` up to `J ≈ 10–20` and above the
  0.18 margin up to `J ≈ 30`. So `+−` positivity would not be supported at **any** `n` for those `J`, not just below `1.5Λ`. This is a projection,
  not a computation.

**The photon data do not decide between A and B.** In §0 item 4 the flat candidates `d̃^J_{22}·E_k` and `d^J_{22}·E_k` did **not** reproduce the
opposite-helicity data: they overcorrect by 2–4×. The curved shift with `c = 2J_z²` did, at small `x`. So the flat density for these AdS tensor
functionals is not established to be the one `F_flat` uses, and the decisive input is open.

**Bottom line [I].**
* The opposite-helicity effect is **not** suppressed on the binding near-gap states: it is flat-level at fixed `J`, not `1/m`.
* Whether the `+−` channel of M02 survives below `1.5Λ` (or, under reading A, at all for `J ≲ 30`) hinges on identifying the bulk-point/flat density of the AdS `+−`
  moment.
* If `F_flat` is right (B), nothing changes relative to `UNIFORM_REMAINDER.md`. If not (A), the `+−` channel fails for every `Δ_g`.

## 5. Code (`ads4/code/skeleton/`)

| file | content |
|---|---|
| `common.py` | paths, grids, kinematics |
| `spin.py` | two-pair kernels on the helicity waves, leading prediction, `Q_h`; block tensors cached in `cache/modes4/` |
| `skeleton.py` | exact scalar skeletons, `Nrat`, CHMRS series ratio |
| `s0_scalar_precision.py`, `s1_spin_scan.py` | precision tests; per-case ν scans |
| `s3_compare.py` | single-case comparison and envelope |
| `s4_summary.py` | all tables (`cache/out_s4_summary.txt`) |
| `s6_opp_scaling.py` | §7: scaling fits of the opposite-helicity non-uniformity; graviton projection (`cache/out_s6.txt`) |
| `s5_spinning_skeleton.py` | §6: spinning-kernel candidates for the `β`-dependence (`cache/out_s5.txt`) |
| `run_all.sh` | everything, at most 2 cores (`OMP_NUM_THREADS=1`) |
