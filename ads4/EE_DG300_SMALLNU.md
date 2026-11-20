# `Ω^EE` at `Δ_g = 300`, `J ≥ 24` near the gap: the graviton diagonal at small `ν` computed (n3 contour + n4 thimbles)

*19–20 Nov 2026, EE-Δ_g=300-small-ν agent. Unrefereed, not committed by me. No existing document or project code edited. New code
`ads4/code/spinfunc/z*_*.py` (+ `z2_queue*.sh`, `z5_queueC.sh`, `z6_jk.sh`), caches `spinfunc/cache/z*` (and `w1_*_J32_k2_z5.pkl`,
`w1_shadow_k2_z5.pkl` from the `J = 32` row, §7). Project modules imported read-only (in-process wrappers only).*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

## 0. Summary

1. **The `J = 24, 28` failure of EE_NEARGAP (convention A, `−0.053 / −0.095`) is an artefact of the `K₂` node at `x = 0.05` and the
   spline below it.** Direct (all-orders) crystal values on the `Δ = 336` rows, state and `Δ_ref = 2·10⁴` on common random numbers,
   give `Γ_EE,++ = +212 / +188` at `ν = 16.7` (`J = 24 / 28`), where `K₂` has `−1331 / −1025`; at `ν = 11.7 … 13.3` the direct values are
   `+250 … +440` where the q3 spline extrapolated to `≈ −2000 … −3000` [NE; scalar checks `≤ 1.4σ`, `|Δc₂| ≤ 100` there].
2. **`n` at `Δ_g = 300` for `J = 24, 28` (min-norm gauge, non-local `X` below `ν = 25`): `n = 1` in every variant computed** [NE + M]:
   * `J = 24`: min `λ/N_s` `+0.84 … +0.86` (estimator 'best', any extrapolation below `ν = 5`; minimum at 1.5Λ); the formerly failing
     1.1Λ state is at `+1.06 ± 0.14` (jackknife); with the scalar control variate everywhere `+0.70 ± 0.10`. With all thimble nodes
     (`ν < 8.3`) dropped: `+0.82`. Failure would need a uniform `Γ_c ≲ −8900` c₂ below `ν = 8.3`; the thimble data there give
     `|Γ_c| ≲ 5000` (`ν = 5.0`: `+350 ± 4900` direct, `−4500 ± 5200` CV) [NE + M]. **Grade for `J = 24`: [NE + M].**
   * `J = 28`: min `+0.87 … +0.90` (best), `+0.89` (direct), `+0.19` (CV everywhere; minimum at 1.1Λ). But the 1.1Λ state is
     `+1.45 ± 0.58` (best) / `+0.19 ± 0.74` (CV): the `ν = 5.0, 6.67` thimble nodes have 10–60 % errors at `J = 28`
     (`Γ_c(5.0) = +10 000 ± 10 000` direct, `−8800 ± 16 000` CV), and without them (`ν ≥ 8.3` only) the 1.1Λ state is `+0.59 … +0.63`
     with failure at a uniform `Γ_c ≈ −4000` c₂ below `ν = 8.3`. **Grade for `J = 28`: passes on every central value [NE + M]; not
     statistically established at 1σ below `ν ≈ 8` [I].**
   * `J = 32` (`Δ = 336` row built, §7): see §7.
3. **Canonical gauge: not controlled.** With `X_can` (`X_can/ψ₂ ≈ −170` at `ν ≈ 5`) the same data give `−22 … −40` (`J = 24`), i.e. the
   canonical-gauge `n` at 300 is decided by `ν ≲ 6`, where the thimble errors are `≥ 5 %` [I].
4. **Reference / method findings [NE].** The fold of the crystal saddle is at `ν* ≈ 11.3–11.5` for every `(Δ, J)` used (no `J`
   dependence, §1). The heavy reference failed because `Q_J` overflows double precision at `Δ = 2·10⁴`, `J ≳ 20` (§2); with the fix,
   `Δ_ref = 2·10⁴` is computable and no change of `N̂` is needed. The n3 contour (`ν_f = 13`) works down to `ν = 8.3` at these `(Δ, J)`
   (`ν = 10`: 0.3–1 % errors; `ν = 8.3`: 2–5 %), better than the n4 thimbles (5–10 %, and C3H 11 % off n3 at `ν = 10`, 2σ) (§4).
5. **Large `J` (task 3) [M].** With direct values the small-`ν` correction is no longer a fixed `≈ −0.01` absolute: at `ν ≥ 8.3` it is
   `+0.1 … +0.5 N_s` at 1.1Λ, of the stabilising sign. What remains is the unresolved `ν < 8` window, whose failure threshold scales like
   `N_s ∝ b̂⁻³` (`−8900` → `−4000` c₂ from `J = 24` to 28). The largest near-gap `b̂` at 300 (`m ≤ 1.3Λ`, twist-allowed) is `≈ 159`
   (`J ≈ 103`); for `J = 32 … 103` (`b̂ = 59 … 159`) the threshold extrapolates to `≈ −2000 … −100` c₂, below the resolution of any
   available thimble value. **The trend can not be bounded from these data; `n = 1` at 300 for all `J ≥ 24` remains [I]** (§8).

## 1. Fold location at large `J` (`z1_fold.py`, `cache/z1_fold.txt`) [NE]

Symmetric saddle of the scalar flow function `Φ_ν` continued downwards from `ν = 16–22` (step 0.25–0.5), soft `4×4` determinant:

| block | `det_soft` at `ν = 12.5 / 12.0 / 11.5` | last symmetric saddle | `A±` pair |
|---|---|---|---|
| `Δ = 336`, `J = 24` | `1111 / 790 / 404` | `ν = 11.5` (Newton fails at 11.0) | exists at 11.2 (asym 0.19), not at 11.45 |
| `Δ = 336`, `J = 32` | identical to `J = 24` to all printed digits | 11.5 | exists at 11.2 (same `Φ₀`) |
| `Δ = 2·10⁴`, `J = 24` | `1113 / — / —` (`954` at 12.25) | scan stopped at 12.25 (same values as `Δ = 336`) | — |
| `Δ = 4000`, `J = 24` | `2329` at 14.25 (`Δ = 336`: `≈ 2330`) | stopped | — |

* `Φ_ν` (n1_direct.phi_w) came out `J`-independent along the symmetric saddle family (identical `det_soft`, `Φ₀` for `J = 24, 32`; the
  mechanism was not analysed), and its soft block is `Δ`-independent to `≲ 0.2 %` between `Δ = 336` and `2·10⁴`. **So `ν* ≈ 11.3–11.5` for every `(Δ, J)` used
  here** (H0 §25.2 had 11.28 at `Δ = 100`, `J = 4`). The fold does not move with `J`; the large-`J` difficulty is in the prefactors
  (`P_J`, the spin-`J` structures), not in the contour.
* Production rule (as NEARGAP100 §8): n3 hybrid contour for `ν ≥ 11` (`ν_f = 13` below 13), `A₊ + A₋` thimbles (n4) below.

## 2. Why the direct method failed at the heavy reference, and the fix (`z2_points.py`) [NE]

* At `Δ = 2·10⁴`, `J = 24`, `ν = 20.85`, every point of the n3 contour (4090/4096; also with the reference's own proposal) was masked.
  Cause: `l1_crystal.build` evaluates the scalar prefactor `Q_J = Σ c_p (Z·W)^p (Z²W²)^{(J−p)/2} = |Z|^J|W|^J P_J(cos)` as a polynomial;
  at the heavy reference `|Z||W| ≈ 10¹³` at the crystal, so `Q_J` overflows double precision for `J ≳ 20` (`inf`/`nan`), and the
  contour's ok-mask (scalar) also gates the graviton kinds. At `J ≤ 10` (`Q ~ 10¹³⁰`) it does not overflow, which is why the
  `J ≤ 10` references of H0 §24–25 worked. (The real-contour failure of n1 at `ν ≲ 35` reported in EE_NEARGAP §5 may be this same
  overflow, or walls; not separated here [I].)
* Fix (in-process copy of `l1_crystal` in `n1_pw` only): `Z·W → Z·W/s`, `Z²W² → Z²W²/s²` with a per-block constant `s = |Z²W²|^{1/2}` at
  the block's saddle, i.e. the scalar integrand is multiplied by the constant `s^{−J}`; the stored scalar normalisation gets `+J log s`.
  The graviton integrands do not use `Q_J` and are unchanged.
* Test (`J = 24`, `ν = 20.85`, `2⁸·16` RQMC points, proposal adapted on the state and re-used): reference `Δ = 2·10⁴` now has 0 masked
  points and relative errors 1.0 % (EE) / 1.1 % (C3H) / 1.0 % (scalar), as the state `Δ = 336` (0.95 / 1.2 / 0.9 %).
* **Consequence for task 1's alternative reference:** the heavy reference is computable as it is; the normalisation `N̂` (state/reference
  ratio at `Δ_ref = 2·10⁴`, H0 §24.1, l9_analysis) need not be re-derived.

## 3. Error structure: common random numbers and the scalar control variate [NE]

At the test point (`J = 24`, `ν = 20.85`, `2⁸·16` points):
* The `2×2` graviton matrices are nearly rank-1 with a common fluctuating factor: element ratios within a matrix have relative errors
  `4·10⁻⁵` (state) / `4·10⁻⁸` (reference), while each element has 0.85 % (state) — all four elements share the same relative error.
* State/reference ratio of an element with common random numbers: 0.85 % (vs 1.4 % if independent). In `Γ` units (`m²×` relative,
  `m = 333.6`) that is `≈ 950` c₂ — too large; it falls as `N^{-1/2}`.
* Graviton/scalar ratio on the same samples: 0.33 % (state) / 0.23 % (reference) — the common factor is largely shared with the scalar.
  Since the exact scalar state/reference ratio is known (exact Mellin skeleton, n4_exact; `J` even), the estimator
  `G = [(g/s)_state/(g/s)_ref] × (s_state/s_ref)_exact` (scalar control variate) is ≈ 2.5× more precise. Both estimators are reported.

## 4. Production and validation (`z2_points.py`, `z3_nodes.py validate` → `cache/z3_validate.txt`) [NE]

* Nodes on the `Δ = 336` rows (`m = 333.60 / 333.28` for `J = 24 / 28`; the 1.1Λ states `m = 328.3` are held at this row by the chain):
  `x = 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.0625, 0.075, 0.1` (`ν ≈ 5.0 … 33.4`), each with `Δ_ref = 2·10⁴` at the same `ν`,
  proposal adapted on the state and re-used for the reference, same Sobol seeds (8 replicas).
  n3 (`2¹¹·2·8` points): `ν ≥ 11.7` (`ν_f = max(ν, 15)` above 13, `ν_f = 13` below), and at `ν = 10.0, 8.3` (`ν_f = 13`, tag `x`); n4
  thimbles (`2⁹·2·8` per thimble, `ν_c = 10.5`): `ν = 10.0, 8.3, 6.67, 5.0`. Where both exist the chain uses n3 (`Z3_PREFER=x`).
* **Scalar check against the exact Mellin object** (`d_D − d_ref`, must vanish; exact `X` from n4_exact.record, own cache `z3_exact.pkl`;
  at `(2·10⁴, 28)` ExactState returns NaN, replaced by the heavy-block value `ρ₂ = 0`, which at `(2·10⁴, 24)` is exact to `|ρ₂| ≤ 5.5·10⁻⁸`):

| `ν` | method | `J = 24`: `Re(d − d_ref)`, `Im`; `Δc₂` | `J = 28`: same |
|---|---|---|---|
| 33.4 | n3 | `−0.0001 ± 0.0008`, `+0.0007 ± 0.0013`; `+2 ± 17` | `−0.0023 ± 0.0081`, `+0.006 ± 0.008`; `−15 ± 54` |
| 25.0 | n3 | `−0.0007 ± 0.0005`; `+30 ± 24` | `−0.0003 ± 0.0004`; `+11 ± 18` |
| 20.8 | n3 | `−0.0024 ± 0.0018`; `+74 ± 55` | `−0.0011 ± 0.0008`; `+46 ± 36` |
| 16.7 | n3 | `+0.0018 ± 0.015`; `−4 ± 35` | `−0.0024 ± 0.0038`; `+53 ± 83` |
| 13.3 | n3 | `−0.0024 ± 0.0017`; `−65 ± 48` | `−0.0036 ± 0.0046`; `−28 ± 36` |
| 11.7 | n3 | `−0.0023 ± 0.0016`; `−99 ± 71` | `−0.0024 ± 0.0027`; `−61 ± 69` |
| 10.0 | n3 (`ν_f = 13`) | `−0.0023 ± 0.0032`, `+0.0001 ± 0.004`; `−132 ± 187` | `−0.0031 ± 0.0044`; `−134 ± 194` |
| 10.0 | n4 | `+0.0018 ± 0.0155`, **`+0.099 ± 0.046`**; `+108 ± 908` | `+0.070 ± 0.039`, `−0.18 ± 0.09`; `+3044 ± 1721` |
| 8.3 | n3 (`ν_f = 13`) | `−0.0041 ± 0.014`; `−301 ± 1020` | `−0.018 ± 0.016`; `−1118 ± 959` |
| 8.3 | n4 | `+0.012 ± 0.018`; `+852 ± 1302` | `+0.0016 ± 0.012`; `+99 ± 743` |
| 6.67 | n4 | `−0.003 ± 0.028`, `+0.16 ± 0.09`; `−267 ± 2433` | `−0.078 ± 0.100`, `+0.34 ± 0.13`; `−6082 ± 7807` |
| 5.0 | n4 | `+0.058 ± 0.036`; `+5568 ± 3483` | `+0.23 ± 0.23`; `+20 939 ± 20 844` |

  All `≤ 1.8σ` in `Re`; the n4 phases are off by `0.1–0.3` rad (≤ 2.6σ, jackknife), as NEARGAP100 §4 saw at `ν = 5.45`. **The n3
  contour is validated to `|Δc₂| ≲ 100` on `ν ≥ 11.7` and `≲ 200` at 10.0; below 10 the scalar test only bounds the method at the
  10³–10⁴ c₂ level.**
* **n3 vs n4 at the same node (CV ratios, element `(0,0)`):** `ν = 10.0`: EE `−6.674` (n4, 1.0 %) vs `−6.769` (n3, 0.3 %) `×10²⁵`;
  C3H `−9.94` (n4, 5.2 %) vs `−8.93` (n3, 0.6 %) `×10²⁸` — n4 C3H 11 % high (≈ 2σ). `ν = 8.3`: EE `−9.32` (3.0 %) vs `−8.99` (2.0 %);
  C3H `−1.168` (2.4 %) vs `−1.108` (1.7 %). n3 is used where available.
* **Estimators** (`z3_nodes.records`): 'dir' = plain CRN ratio; 'cv' = scalar control variate (`G × (s_st/s_ref)_exact/(s_st/s_ref)_MC`);
  **'best'** = per node and kind the one with the smaller jackknife error (CV wins on the n3 nodes, dir on the n4 nodes `ν ≤ 6.7` and
  for C3H at `ν ≤ 8.3`; list in the chain output). Errors propagated by leave-one-replica-out chain runs (`z6_jk.sh`).

## 5. Node values (`Δ = 336` row; `cache/z4_chain_f4_besthold00.txt`, estimator 'best')

`Γ_c` = complex forward-normalised correction `m²(G/a_h − Â_h[f_LO])/env` (`Re Γ_c` = the chain's local `Γ`); `Γ_NL = a_χΓ_c` (finite
through the zero `ν = 6.643` of `h₋ + h₊`); NL density `Re[(Y/φ)Γ_NL] = Re[(X/φ)Γ_c]` (min-norm `X`); C3H local (no de-shift). `++`
shown (`+−` within 10–40 c₂ of it except where noted). Node errors (c₂, Γ units) ≈ `m²·δ_rel·|D_h|`: `≲ 30–100` for `ν ≥ 11.7`,
`≈ 200` at 10, `≈ 1500` at 8.3, `≈ 5000` at 5.0 (`J = 24`), `≈ 10⁴` at 5.0 (`J = 28`).

| `ν` | `J = 24`: `Γ_c` / NL dens / C3H | `J = 28`: `Γ_c` / NL dens / C3H | `K₂` (`J = 24 / 28`, `Γ_c`) | `Re X/φ` |
|---|---|---|---|---|
| 33.4 | `−240 / −241 / −242` | `−252 / −253 / −228` | `+39 / −93` | 1.004 |
| 25.0 | `−102 / −102 / −137` | `−233 / −235 / −253` | `0 / +407` | 1.008 |
| 20.8 | `+47 / +48 / −18` | `−56 / −57 / −111` | — | 1.013 |
| 16.7 | `+212 / +218 / +219` | `+188 / +193 / +130` | **`−1331 / −1025`** | 1.030 |
| 13.3 | `+254 / +266 / +150` | `+304 / +317 / +197` | (spline `≈ −2600`) | 1.044 |
| 11.7 | `+286 / +278 / −18` | `+337 / +328 / −33` | | 0.977 |
| 10.0 (n3) | `+132−51i / +91 / +130` | `+142−44i / +98 / +87` | | 0.710 |
| 8.3 (n3) | `−395+504i / −37 / −514` (C3H dir) | `−888+703i / −112 / +401` | | 0.186 |
| 6.67 (n4) | `≈ −1.3·10⁶` (pole) / `+9549` / `+5542` | `≈ −2·10⁶` / `+11 683` / `+5274` | | −0.007 |
| 5.0 (n4) | `+301−2788i / +844 / −6402` | `+9957+769i / +12 289 / −152` | | 1.242 |

(`K₂` columns at the shared nodes are from the run with `K₂` kept there, `i4_cvhold_k2`; `a_h` differs slightly between the runs.)
**The direct `Γ` is of moderate size (`|Γ| ≲ 450` on `ν ≥ 10`) and positive on `ν ≈ 10–17`, where `K₂` (and its spline) had
`−1000 … −3000`.** Near `ν = 6.6` `Γ_c` has the pole of the forward normalisation, as it must; the NL density is finite there (`X ≈ 0`
on `ν ≈ 6–8` makes it irrelevant).

## 6. `n` at `Δ_g = 300`, `J = 24, 28` (`z4_chain.py`; `cache/z4_chain_f*.txt`; states and conventions as EE_NEARGAP's final run:
all exact-skeleton + hybrid states `m/Λ = 1.02 … 5`, rows 336/500/1600, EE basis, Hampel; odd `J = 5` included for the binding state)

Reproduction check: `Z4_DIRECT=0 Z4_NL=0` gives `−0.0534 / −0.0950` (EE_NEARGAP convention A) exactly.

| variant (min-norm gauge unless stated) | `J = 24` min `λ/N_s` [`m/Λ`] | `J = 28` min [`m/Λ`] | `n` at 300 |
|---|---|---|---|
| EE_NEARGAP A (`K₂`, spline, local) | `−0.053` [1.1] | `−0.095` [1.1] | 1.11 |
| EE_NEARGAP C (min-norm, spline) | `+0.043` [1.1] | `+0.087` [1.1] | 1 |
| **direct nodes, 'best', `Γ_NL` held below 5.0** | **`+0.845`** [1.5] | **`+0.887`** [1.5] | **1** (min `+0.2052`, `J = 5`) |
| 'best', zero below 5.0 | `+0.841` [1.5] | `+0.869` [1.5] | 1 |
| 'best', linear below 5.0 | `+0.863` [1.5] | `+0.897` [1.75] | 1 |
| 'dir', held | `+0.846` [1.5] | `+0.887` [1.5] | 1 |
| 'cv', held | `+0.699` [1.1] | `+0.191` [1.1] | 1 |
| 'best', only `ν ≥ 8.3` direct (n4 nodes dropped), held | `+0.815` [1.3] | `+0.591` [1.1] | 1 |
| same, `Γ_c := −4000` below 8.3 | `+0.477` [1.1] | `+0.007` [1.1] | 1 |
| same, `−6000` / `−8000` / `−12000` | `+0.284 / +0.090 / −0.308` | `−0.306 / −0.618 / −1.244` | |
| canonical gauge (`X_can`), 'cv', held | `−37.5` [1.1] | `−33.1` [1.1] | 1.55 |

**Statistical errors** (leave-one-replica-out, `z6_jk.sh`; at the formerly failing states): 'best' `J = 24`: 1.1Λ `+1.06 ± 0.14`,
1.2Λ `+0.94 ± 0.07`; `J = 28`: 1.1Λ `+1.45 ± 0.58`, 1.2Λ `+1.25 ± 0.29`. 'cv': `J = 24` 1.1Λ `+0.70 ± 0.10`; `J = 28` 1.1Λ
`+0.19 ± 0.74`. The `J = 28` spread comes from the `ν = 5.0, 6.67` thimble nodes (one replica moves the 1.1Λ state by `−0.6`).

**Reading.**
* `J = 24`: positive in every variant by `≥ 5σ` (jackknife), including all thimble nodes, without them, and with the extrapolation
  below `ν = 5` held / zeroed / linear. Failure needs `Γ_c ≲ −8900` uniformly below `ν = 8.3`; the measured `ν = 5.0` value is
  `+350 ± 4900` (dir) / `−4500 ± 5200` (cv). **`n = 1` at 300 for `J = 24` [NE + M]** (M = the chain model of S1_S4/EE_NEARGAP).
* `J = 28`: every central value passes (lowest `+0.19`, CV); but the 1.1Λ state's error is `0.6–0.7 N_s`, and the decisive window
  `ν < 8.3` (failure at uniform `Γ_c ≈ −4000`) is measured only to `±10⁴` c₂. **`n = 1` at 300 for `J = 28`: [NE + M] on central values,
  [I] at the 1σ level below `ν ≈ 8`.**
* The gauge question: min-norm `n = 1`; canonical `n = 1.55` because `X_can/ψ₂ ~ 10²` at `ν ≲ 6` multiplies thimble values with 5–100 %
  errors — **canonical [I]**, not decided by these data.

## 7. `J = 32`

* No `Δ_g = 300` state-mass row existed for `J = 32`. Built here: `K₂` EE + C3H row on `Δ = 336` (`m = 332.9`, 19 nodes, `w1_rows.py`
  with tag `_z5` → `cache/w1_pts_{EE,C3H}_J32_k2_z5.pkl`; the w2 loader picks these up), then direct n3 nodes `ν = 11.65 … 33.3` and
  n3 (`ν_f = 13`) at `ν = 10.0, 8.3` (`z5_queueC2.sh`). Thimble nodes `ν < 8` were not run for `J = 32` (cost; and §6 shows they do
  not resolve the window at `J = 28` anyway).
* **Baseline with the new row, `K₂` only, EE_NEARGAP convention A** (`z4_chain_j32_base2`, rows 336/1040; the 1600 row's `+−`
  breakdown at `x ≥ 0.5` is avoided as in EE_NEARGAP): `J = 32` fails at 1.15Λ (the lowest twist-allowed chain mass), `λ/N_s = −0.062`,
  with `dF_++/N = −1.06` (the same small-`ν` pattern as `J = 24, 28`) and `dF_+−/N = +4.9`; 1.2Λ `+0.063`, 1.3Λ `+0.38`, 1.5Λ `+0.72`
  [NE + M, convention A].
* Direct-node results: **(pending — see hand-over)**.

## 8. Task 3: does every `J ≥ 24` near the gap at 300 pass? [M + I]

| `J` | twist-allowed `m_min/Λ` | `b̂_max = (2J+1)Λ/m_min` | `N_s ≈ 1024/b̂³` | failure threshold below `ν = 8.3` (uniform `Γ_c`) |
|---|---|---|---|---|
| 24 | 1.077 (chain state 1.1) | 45.5 | `1.1·10⁻²` | `≈ −8900` c₂ (computed) |
| 28 | 1.090 (1.1) | 52.3 | `7.2·10⁻³` | `≈ −4000` (computed) |
| 32 | 1.102 | 59.0 | `5.0·10⁻³` | (§7) |
| 40 | 1.126 | 71.9 | `2.8·10⁻³` | `~ −10³` [M: `∝ N_s`] |
| 64 | 1.195 | 107.9 | `8·10⁻⁴` | `~ −3·10²` [M] |
| 103 | 1.300 (last `J` with a twist-allowed state `≤ 1.3Λ`) | **159** | `2.5·10⁻⁴` | `~ −10²` [M] |

* **Largest `b̂` near the gap** (`m ≤ 1.3Λ`, twist-allowed, `Δ_g = 300`, `Λ = 298.5`): `b̂ ≈ 159` at `J ≈ 103`; above that the band
  (`m > 1.3Λ`) is covered by the leading-Regge band model and EE_NEARGAP §6, where the computed beyond-Regge corrections are small.
* **What changed relative to EE_NEARGAP's mechanism.** The "J-independent `≈ −0.01` absolute" small-`ν` correction was the `K₂`
  value at `x = 0.05` and its spline. With direct values the `ν ≥ 8.3` part is of moderate size and of the stabilising sign at
  `J = 24, 28` (the 1.1Λ states move from `−0.05/−0.10` to `+0.6 … +1.5`). The direct `Γ` on `ν = 10 … 33` is `O(10²)` c₂ and does
  not grow from `J = 24` to 28 (`ν = 13.3`: `+254 / +304`; `16.7`: `+212 / +188`; `33.4`: `−240 / −252`).
* **What can not be bounded.** The window `ν ≲ 8` carries `Σψ₂ Re X/ψ₂ ≈ 0.065` (vs `N_s ≈ 0.01`), with `Re X/ψ₂ ≈ 2–5` at `ν ≲ 4`
  in the min-norm gauge. A fixed beyond-Regge `Γ_c` of size `C` there shifts `λ/N_s` by `≈ 2·0.065 C/(m² N_s)`, i.e. the failure
  threshold scales like `N_s ∝ b̂⁻³` (`−8900 → −4000` from `J = 24` to 28, consistent with `N_s`: 0.0119 → 0.0076 and the different
  `ν ≥ 8.3` margins). For `J ≳ 40` the threshold is `≲ 10³` c₂ — below what the thimble method resolves (`±5·10³` at `J = 24`,
  `±10⁴` at `J = 28`, `ν = 5`) and below the size `K₂`-type estimates give at `ν ≈ 5–8`. **So the trend can not be bounded from these
  data: `n = 1` at 300 for all `J ≥ 24` near the gap is [I].** A bound needs either (i) the `ν < 8` beyond-Regge `Γ_c` to `≲ 10²` c₂
  (thimbles with ~10³× the statistics, or a better small-`ν` representation), or (ii) an analytic statement that the beyond-Regge
  correction vanishes as `ν → 0` fast enough to beat `X_min/ψ₂ ~ 3–5` there.

## 9. Graded statements

| statement | grade |
|---|---|
| fold `ν* ≈ 11.3–11.5`, independent of `J` (24, 32) and `Δ` (336, 4000, 2·10⁴) | [NE] |
| the n1/n3 failure at the heavy reference for `J ≳ 20` is the `Q_J` double-precision overflow; fixed by a constant rescaling; `Δ_ref = 2·10⁴` then computable, no change of `N̂` needed | [NE] (cause [D] from the code; n1's real-contour masking at `ν ≲ 35` not re-examined) |
| n3 contour (`ν_f = 13`) valid down to `ν = 8.3` at `(336, 24/28)`: scalar `|Δc₂| ≲ 200` on `ν ≥ 10`, `≲ 1100 ± 1000` at 8.3 | [NE] |
| n4 thimbles at `ν = 5–10`: scalar `Re` within 1.8σ, phases off 0.1–0.3 rad; C3H 11 % off n3 at `ν = 10` (2σ) | [NE] (weak) |
| `K₂` at `x = 0.05` on the 336 rows has the wrong sign (`−1331/−1025` vs direct `+212/+188`) | [NE] |
| `J = 24` at 300 passes (min `+0.84`, 1.1Λ `+1.06 ± 0.14`; every estimator / extrapolation; failure needs `Γ_c ≲ −8900` below 8.3) | [NE + M] |
| `J = 28` at 300 passes on every central value (`+0.19 … +0.90`) | [NE + M] |
| `J = 28` robust against the `ν < 8.3` window (needs `Γ_c > −4000`; data `±10⁴`) | [I] |
| `J = 32` | §7 |
| `n = 1` at 300 for `J ≤ 28` (min-norm gauge, non-local `X`) | [NE + M] (central), with the `J = 28`, `ν < 8` caveat [I] |
| canonical gauge at 300 | [I] (`X_can` amplifies the `ν ≲ 6` thimble region; computed values `−22 … −40`) |
| all `J ≥ 24` near the gap at 300 (`b̂ ≤ 159`) | [I] (threshold `∝ N_s ∝ b̂⁻³`; not bounded) |

## 10. Hand-over

* **Code** (`ads4/code/spinfunc/`): `z1_fold.py` (fold scan; `cache/z1_fold.txt`); `z2_points.py` (direct points: n3 contour above
  `Z2_NUSPLIT`, n4 thimbles below; CRN; the `Q_J` overflow fix `patch_qscale`/`set_qsc`; env `Z2_M`, `Z2_R`, `Z2_SEED`, `Z2_KINDS`,
  `Z2_TAG` for separate caches; skips finished nodes); `z3_nodes.py` (scalar validation `validate J…` → `cache/z3_validate.txt`;
  chain records, estimators `dir`/`cv`/`best`, jackknife; `Z3_PREFER=x` (default) lets tagged caches override); `z4_chain.py` (w2 chain
  + direct nodes + non-local `k = 2` term below `Z4_NUNL = 25`; `Z4_EST`, `Z4_EXTRAP=hold|zero|lin|const` (+`Z4_CONST`), `Z4_GAUGE=min|can`,
  `Z4_NUMIN/NUMAX`, `Z4_OVERRIDE`, `Z4_NL`, `Z4_DIRECT`, `Z4_JK`, `Z4_JS`); `z6_jk.sh EST EXTRAP TAG` (jackknife chain runs).
* **Caches:** `cache/z2_pts_{EE,C3H}_J{24,28,32}.pkl`, `cache/z2x_pts_*` (n3 at `ν = 10, 8.3`), `cache/z2s2_pts_*` (repeat of the
  `J = 28` thimble nodes, new seed), `cache/z2*_raw.pkl`, `cache/z3_exact.pkl`, `cache/z4_chain_<tag>.{txt,pkl}` (`f4_besthold00` =
  main), `cache/z6_jk_*.log`, logs `cache/z2_queue*.log`, `cache/z5_queueC*.log`.
* **Main run:** `Z4_EST=best Z4_JS=24,28 Q3_JEVEN=24,28 Q3_JODD=5 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
  Z4_TAG=x python3 z4_chain.py 300` (≈ 1 min). `Z4_DIRECT=0 Z4_NL=0` reproduces EE_NEARGAP A exactly.
* **Container restart (20 Nov ≈ 00:25)** killed the first queues; resumed with `z2_queueR{1,2}.sh`; `z2_queue{A,B,B2}.sh`, `z5_queueC.sh`
  are superseded (C's tail replaced by `z5_queueC2.sh`).
* **Next steps (value order):** (1) the `ν < 8` window at `J = 28` (and larger `J`): thimble statistics ×10–100, or a better
  representation below the fold (the n3 `ν_f = 13` contour still had 2–5 % at 8.3; try `ν_f = 12` or complex `ν_f` at 6–8); (2) `J = 32`
  thimble nodes and the 500 row's small-`ν` nodes (states at 1.15–1.3Λ use the 500 row with weight 0.1–0.6, still `K₂` there); (3) an
  analytic small-`ν` statement for the beyond-Regge `Γ` (needed for any all-`J` bound); (4) the canonical gauge needs `ν ≲ 6` at `≲ 1 %`.
