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

## 5. Code (`ads4/code/skeleton/`)

| file | content |
|---|---|
| `common.py` | paths, grids, kinematics |
| `spin.py` | two-pair kernels on the helicity waves, leading prediction, `Q_h`; block tensors cached in `cache/modes4/` |
| `skeleton.py` | exact scalar skeletons, `Nrat`, CHMRS series ratio |
| `s0_scalar_precision.py`, `s1_spin_scan.py` | precision tests; per-case ν scans |
| `s3_compare.py` | single-case comparison and envelope |
| `s4_summary.py` | all tables (`cache/out_s4_summary.txt`) |
| `run_all.sh` | everything, at most 2 cores (`OMP_NUM_THREADS=1`) |
