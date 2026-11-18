# `⟨TTTT⟩` beyond one point: H1 (heavy formula at more `k, ν, β, m`) and H2 (reading B for gravitons)

Status: written 18 Nov 2026 by a derivation/numerics agent; §0–§5 committed (c177cfc); §6 committed (bfdf768); **§7 added afterwards, not
refereed, not committed.** §6 corrects the even-kernel statements of §0–§4: the `t16` even kernels used mis-normalised circle weights. Code:
`ads4/code/fourpt/tttt/` (`tmoments_l2.py` production, `t16_kernels_ext.py`, `t17_readingB.py`, `t18_kstar.py`, `t19_k2.py`; logs and caches
in `tttt/cache/`, git-ignored). Builds on
`SPINNING_TTTT.md` (Taylor-mode waves, `l = 2` vertices `T^{hLhR} = 8√(λ_{hL}λ_{hR})/m⁸`, Fock data), `SPINNING_ODDK.md` (robust extraction)
and `SKELETON_TEST.md` §8 (reading A/B for currents). Targets hypotheses H1 and H2 of `RESULT_D.md`.

Marks: **[D]** derived; **[DVE]** derived and verified end to end on exact blocks; **[NE]** numerical evidence; **[ID]** closed form identified
from exact data; **[I]** inferred; **[ND]** not done.

## 0. Summary

1. **H1, heavy formula beyond one point.** Robust extraction: fixed kernel contractions and leading modes, no fit. The new production points
   are `β = 1`, `m = 100.7` and `β = 0.4`, `m = 72.0`. The data now cover `β = 1` at `m = 39, 70, 101` and `β = 0.4` at `m = 72`, for
   `k = 2…7` and `ν = 0, 1, 3`.
   * **[DVE] `k = 4…7`, `β = 1`, `ν = 0, 1, 3`, both helicities, odd kernel `𝔄⊗𝔄` and leading modes.**
     - `R = 1 + c₂/m²`, with `m²(R−1)` flat across the three masses. For example `k = 6`, `ν = 0`, same helicity: `+124, +116, +113`;
       opposite helicity: `−3, −7, −9`.
     - `|c₁| ≤ 0.04`, against `≈ 4` for the `b^{(2)}` mutation.
     - At `m = 101`: same helicity `1.007–1.018`, opposite helicity `0.970–1.011`. The `ν = 3` opposite-helicity entries are the largest
       deviations: `c₂ ≈ −140` to `−305`.
   * **[DVE] `k = 3`** for same helicity at all `ν`, with `m²(R−1) ≈ 175–250` flat. Opposite helicity under `𝔄⊗𝔄` is within `1–1.5 %` at
     `m ≥ 70` for `ν ≤ 1`. **`k = 2` is resolved only for same helicity.** Opposite helicity is dominated by the dropped tiny-`t` region
     (ratios `O(1)` off).
   * **Even kernels: superseded by §6.** The even-kernel numbers originally given here, in §2 and in §4 used the naive pointwise circle
     weights of `t16`. Those over-weight the isotropic mode `36×` per pair, and that version of "FF" is not flip-free.
     * **With the project's tensor weights, `𝔖_ff⊗𝔖_ff` passes H1 for both helicities [DVE]** (`k = 4–7`, `ν ≤ 3`, `β = 1`, `|c₁| ≤ 0.13`),
       except `k = 4`, `ν = 3`, opposite helicity (`c₁ = 1.9`).
     * The new flip-free kernel without an isotropic mode, `K* = (D₁+D₂)(D₁−½D₂)`, passes everywhere at `k ≥ 4` (`|c₁| ≤ 0.05`).
     * The opposite-helicity failure was the isotropic mode, contaminated by the dropped tiny-`t` region (an artefact).
     * The leading heavy formula holds with `c₂/m²` corrections, `|c₂| ≤ 310` at `ν ≤ 3`.
     * `k = 2` (`C⁺₂`): resolved in §7 under `K*` and `𝔄` (both helicities), after a cleaned extraction that removes two defects of the
       cached waves at the `k = 2`-weighted edges. `c₂ ≈ 140–270`. `𝔖_ff`'s opposite-helicity `k = 2` remains untestable.
   * **`β = 0.4`, `m = 72`** (single mass, `J = 14`):
     - same helicity `1.013–1.030` (all kernels, `k = 4–7`);
     - opposite helicity, odd kernel `0.93–1.05` (`0.986–1.007` at `ν ≤ 1`, `k ≥ 5`);
     - opposite helicity, even kernels (§6.4): `K*` `0.986–1.096`; `𝔖_ff` the same at `k ≥ 5`, but fails at `k = 4`.

     The `m = 42` point was not obtained (two container restarts), so no `β = 0.4` convergence is shown.
2. **The `1/4` constant** is still not tested by blocks, since the waves use unit-normalised `X_h`. Its two inputs are strengthened:
   * the Fock closed form `p[T^{E,h}] = A_hP_FK(1,1)/2304` now holds exactly on 59 entries (`N ≤ 8`, `Δ ≤ 14`; 20 new) [ID];
   * `E → 8` on a second ray (`β = 0.9`): `E = 8.678, 8.162` at `m = 44, 88`, Richardson `7.991` [NE].
3. **H2, reading B for gravitons, at `b̂ = 1.04`, `k = 6`, odd kernel, `m = 39, 70, 101` [NE].** The single ratio is limited by
   helicity-independent systematics: the same-helicity control reaches `|m²(r−1)| ≈ 12–92`. The **double ratio
   `ρ = r_opp/r_same`** cancels them, because all readings predict `D = P_J` for `J_z = 0`.
   * At `x = 0.1`: `m²d_B = −10.8, −16.7, −7.2`, which stays bounded, while `m²d_A = +18, +44, +98` grows like `m²` and `m²d_B̃ = −76, −220,
     −428` grows like `m²`.
   * The measured `ρ − 1` is `0.63, 0.72, 0.93` of the reading-B helicity term (`+0.019, +0.012, +0.010`).
   * At `x = 0.05`: `m = 39` gives `0.98` of the B term (`m²d_B = −3.9`, against `+176` for A). At `m = 70` the B term is small
     (`−0.0044`) and the data give `−0.0083` (`m²d_B = −19`, `m²d_A = −41`). At `m = 101`, `P_J(1−2x²)` is near a zero (condition `0.26`), so
     the point is not usable.
   * **Reading B (full `d^J_{44}`) is supported, and A and B̃ are disfavoured, for gravitons [NE, three masses, one `b̂`, `x ≤ 0.1`].** The
     approach to the B value is not uniform (`63 → 93 %` of the term at `x = 0.1`), and `ν = 7–10` (`x = 0.1`, `m ≥ 70`) is beyond the
     quadrature checks of the grid.

## 1. Data and method

**Production.** `⟨TTTT⟩` waves `W[X_{hL},X_{hR}]` (CHL's exact helicity structures, all 196 word pairs by Taylor mode, `SPINNING_TTTT.md` §4)
on the moment grid of `e2e/moments.py`, float64, `N = 8000`:

| `β` | `m` | `(Δ, J)` | grid points | wall time (1 core, shared machine) |
|---|---|---|---|---|
| 1 | 39.3 | (45.861, 20) | 3336 | (earlier) |
| 1 | 70.0 | (80.484, 36) | 3475 | (earlier) |
| 1 | 100.7 | (115.108, 52) | 3272 | 8820 s (new) |
| 0.4 | 42.2 | (44.565, 8) | 3463 | **not completed**: killed twice by container restarts at about 60 % (2.5 h per attempt); skipped |
| 0.4 | 72.0 | (74.964, 14) | 3397 | 8970 s (new) |

The production code now builds the 25 seed blocks one at a time and precomputes their radial sums for all points. Memory is
0.25–0.9 GB per run, against 3.8 GB before, which allows two runs side by side. The cached modes are reproduced bit for bit.

`β ≈ 0.4` is the coordinator's first choice, with masses chosen so that `A_opp` is not tiny. At fixed `β`, `J ∝ m`, so `J = 8, 14` are
forced. There `A_opp = 0.14, 0.33`, and the Weyl-odd `O(1/m)` term of the integrand is large (`λ_o/λ_s = 8.9, 2.9`; `SPINNING_TTTT.md` §3.6).

**Extraction (robust; `REFEREE_SPINNING_ODDK.md` M1).** No class fit is used. A fixed kernel contraction of the measured modes is compared
with the same contraction of the predicted modes. On the polarisation circles (`SPINNING_4PT.md` §4.5) `D_1 = cos(θ_out − θ_in)`, `D_2 = −1` and
`D_5 = sin(θ_out − θ_in)`. Per pair:

| kernel | definition | modes `(q_in, q_out)`: weight |
|---|---|---|
| E | `D_1²` (the project's even "`D₁⊗D₁`" for `l = 2`) | `(0,0)`: ½; `(±2,∓2)`: ¼ |
| FF | `𝔖_ff = D_1² + D_1D_2 = D_1(z_in·z_out)` (flip-free, `SUPERCONVERGENT_UPLIFT.md`) | E plus `(±1,∓1)`: −½ |
| O | `𝔄 = D_5D_1` | `(2,−2)`: `1/(4i)`; `(−2,2)`: `−1/(4i)` |

These are used as `E⊗E`, `FF⊗FF` and `O⊗O`, together with the leading single mode of each wave. **Correction (§6.2):** the E and FF
weights in this table are the pointwise Fourier coefficients. The project's tensor kernels need an extra factor `(1/36, 1/16, 1)` for
`|q| = 0, 1, 2`, so they are `D₁²`: `(0,0)` `1/72`, `(±2,∓2)` `¼`; `𝔖_ff`: additionally `(±1,∓1)` `−1/32`. `O` is unaffected, since it has
only `|q| = 2`. The prediction is that of
`SPINNING_TTTT.md` §4.5, with `a = k + 5` and nothing fitted.

The tiny-`t` (dropped) region is estimated only for the leading mode. The `oddk` estimator fails for modes that are small in the bulk: it
exceeds the `(0,0,0,0)` mode of `W[o,o]` thirtyfold. Kernel ratios are therefore raw.

## 2. H1: results

**Odd kernel `𝔄⊗𝔄` and leading modes, `β = 1`** (`R` at `m = 39.3 / 70.0 / 100.7`; `m²(R−1)` at the same masses; `c₁` from the fit
`c₁/m + c₂/m² + c₄/m⁴`). The leading-mode values agree with `𝔄⊗𝔄` to `≤ 10⁻³` at `k ≥ 4`.

| `k`, `ν` | same: `R` | same: `m²(R−1)`; `c₁` | opp: `R` | opp: `m²(R−1)`; `c₁` |
|---|---|---|---|---|
| 2, 0 | 1.099 / 1.030 / 1.012 | 153, 146, 123; −0.9 | 0.11 / −1.52 / 1.03 | not resolved (tiny-`t`) |
| 3, 0 | 1.129 / 1.038 / 1.017 | 200, 187, 176; −0.34 | 1.077 / 1.013 / 1.010 | 119, 63, 103; noisy |
| 3, 1 | 1.132 / 1.039 / 1.019 | 204, 193, 195; +0.19 | 1.097 / 1.013 / 1.015 | 150, 66, 150; noisy |
| 4, 0 | 1.119 / 1.035 / 1.017 | 184, 173, 170; −0.04 | 1.057 / 1.017 / 1.008 | 88, 81, 79; −0.04 |
| 4, 3 | 1.117 / 1.035 / 1.016 | 180, 170, 166; −0.09 | 0.942 / 0.984 / 0.992 | −90, −76, −78; −0.22 |
| 5, 0 | 1.103 / 1.031 / 1.015 | 160, 150, 147; −0.000 | 1.032 / 1.009 / 1.004 | 49, 42, 40; −0.005 |
| 5, 1 | 1.112 / 1.033 / 1.016 | 174, 164, 161; +0.008 | 1.056 / 1.017 / 1.008 | 87, 81, 80; +0.005 |
| 5, 3 | 1.086 / 1.026 / 1.012 | 133, 126, 124; −0.002 | 0.903 / 0.971 / 0.986 | −150, −140, −138; −0.009 |
| 6, 0 | 1.080 / 1.024 / 1.011 | 124, 116, 113; +0.003 | 1.000 / 0.999 / 0.999 | 1, −7, −9; −0.004 |
| 6, 1 | 1.093 / 1.027 / 1.013 | 144, 135, 133; +0.005 | 1.027 / 1.007 / 1.003 | 42, 36, 35; +0.000 |
| 6, 3 | 1.045 / 1.013 / 1.006 | 70, 66, 65; +0.002 | 0.857 / 0.956 / 0.979 | −222, −215, −214; −0.001 |
| 7, 0 | 1.051 / 1.014 / 1.007 | 78, 69, 67; −0.000 | 0.964 / 0.986 / 0.993 | −56, −68, −71; −0.011 |
| 7, 1 | 1.066 / 1.019 / 1.009 | 103, 94, 92; +0.002 | 0.993 / 0.996 / 0.998 | −12, −19, −21; −0.004 |
| 7, 3 | 0.995 / 0.998 / 0.999 | −7, −11, −12; +0.000 | 0.804 / 0.938 / 0.970 | −304, −305, −306; −0.002 |

**Even kernels, `β = 1`, naive weights. SUPERSEDED by §6.4** (`R` at `m = 39 / 70 / 101`). Kept for the record. The `𝔖_ff` columns
here are the naive "FF", not the project's `𝔖_ff`.

| `k`, `ν` | same, `D₁²` | same, `𝔖_ff` | opp, `D₁²` | opp, `𝔖_ff` |
|---|---|---|---|---|
| 5, 0 | 1.104 / 1.031 / 1.014 | 1.111 / 1.032 / 1.015 | 0.41 / 1.16 / 1.011 | 0.97 / 1.02 / 1.07 |
| 6, 0 | 1.082 / 1.024 / 1.011 | 1.097 / 1.026 / 1.012 | 0.785 / 1.015 / 1.003 | 0.969 / 1.004 / 1.024 |
| 6, 1 | 1.095 / 1.028 / 1.013 | 1.129 / 1.035 / 1.017 | 1.068 / 0.999 / 1.005 | 1.029 / 1.003 / 1.022 |
| 6, 3 | 1.058 / 1.017 / 1.008 | **2.18 / 1.35 / 1.17** | 0.25 / 0.96 / 1.07 | −0.17 / 0.72 / 2.06 |
| 7, 0 | 1.053 / 1.014 / 1.007 | 1.083 / 1.019 / 1.008 | 0.944 / 0.988 / 0.994 | 0.977 / 0.990 / 0.999 |
| 7, 1 | 1.068 / 1.020 / 1.009 | 1.140 / 1.035 / 1.016 | 1.000 / 0.996 / 0.998 | 0.984 / 0.989 / 0.992 |
| 7, 3 | 1.005 / 1.000 / 1.000 | **1.99 / 1.31 / 1.15** | 0.81 / 0.94 / 0.99 | 0.16 / 0.68 / 1.12 |

* At `k ≤ 4` the opposite-helicity even-kernel ratios are `O(1)` off or worse (full table in the output). Same helicity under `𝔖_ff` at
  `ν = 3` converges like `1/m²` with `m²(R−1) ≈ 1500–1800`: `c₁ = −0.07` (`k = 7`) and `−1.2` (`k = 6`).
* The opposite helicity, which the even kernels see through the isotropic `(0,0)` and `(±1,∓1)` modes, does not converge cleanly.
* Mixed wave, `|K(W[s,o])|/√|K(ss)K(oo)|`, measured against predicted at `k = 6`, `ν = 0`:
  - `m = 39`: `D₁²` `1.1·10⁻²`/`1.3·10⁻⁴`; `𝔖_ff` `1.9·10⁻³`/`5.7·10⁻⁵`.
  - `m = 101`: `D₁²` `2.7·10⁻⁴`/`1.3·10⁻⁴`; `𝔖_ff` `2.5·10⁻³`/`5.7·10⁻⁵`.

**`β = 0.4`, `m = 72.0`, `J = 14`** (one mass):
* same helicity, all kernels at `ν ≤ 1`: `1.013–1.026` (`k = 4–7`);
* `ν = 3`: `1.024–1.030`, except naive FF, which is `1.27–1.40` (the correct `𝔖_ff` gives `1.025–1.031`, §6.4);
* opposite helicity, `𝔄⊗𝔄`: `0.935–1.048` (`0.986–1.007` at `ν ≤ 1`, `k ≥ 5`);
* opposite helicity, even kernels: not resolved with naive weights; resolved for `K*` (and `𝔖_ff` at `k ≥ 5`) with the correct weights
  (§6.4).

Without the `m ≈ 42` point there is no `β = 0.4` convergence statement.

## 3. H2: reading A / B / B̃ for gravitons at fixed `b̂ = 1.04` (`t17_readingB.py`, `out_t17_readingB.txt`)

**Observable.** `R(ν)` = measured/leading prediction of the opposite-helicity wave under `𝔄⊗𝔄`, `k = 6`. `r = R(ν)/(R(0)ρ_skel(ν))`,
where `ρ_skel` is the exact scalar skeleton ratio at `Df = 3` (`Df = 2` as a check). Readings:
* A: `D = P_J`;
* B: `D = d^J_{44} = cos⁸(θ/2)P^{(0,8)}_{J−4}(cos θ)`, the full graviton Wigner function with `J_z = 4`;
* B̃: `D = d̃^J_{44} = P^{(0,8)}_{J−4}`;

all with `cos θ = 1 − 2x²`, `x = ν/m`, at fixed `β = 1` (`b̂ = 1.04`).

**Single ratio** (`Df = 3`; `m²d_X = m²(r − D_X/P_J)`). It is not conclusive. The same-helicity control (`J_z = 0`, every reading predicts
`r = 1`) gives `m²(r_same − 1) = −13.6, −5.7, −44.3` (`x = 0.05`) and `−11.9, −32.3, −91.8` (`x = 0.1`) at `m = 39, 70, 101`. These are
helicity-independent systematics of the same size as the A/B discrimination.

**Double ratio `ρ = r_opp/r_same`** (independent of `Df`, since the skeleton cancels):

| `m` | `x` | `ν` | `ρ − 1` | reading-B term `d^J_{44}/P_J − 1` | fraction of the B term | `m²d_A` | `m²d_B` | `m²d_B̃` |
|---|---|---|---|---|---|---|---|---|
| 39.3 | 0.05 | 1.97 | +0.1139 | +0.1165 | 0.98 | +176 | −3.9 | −21 |
| 70.0 | 0.05 | 3.50 | −0.0083 | −0.0044 | 1.9 | −41 | −19 | −69 |
| 100.7 | 0.05 | 5.04 | (`P_J` near a zero, condition 0.26: not usable) | | | | | |
| 39.3 | 0.1 | 3.93 | +0.0117 | +0.0187 | 0.63 | +18 | −10.8 | −76 |
| 70.0 | 0.1 | 7.00 | +0.0089 | +0.0123 | 0.72 | +44 | −16.7 | −220 |
| 100.7 | 0.1 | 10.08 | +0.0097 | +0.0104 | 0.93 | +98 | −7.2 | −428 |

* At `x = 0.1`, `m²d_B` stays bounded (`|m²d_B| ≤ 17`), while `m²d_A` grows like `m²`: `ρ − 1 ≈ 0.01` at fixed `b̂`, an `O(1)`-relative
  helicity effect that reading A lacks. `m²d_B̃` grows like `m²` (ratios `2.9`, `5.6` against `m²` ratios `3.2`, `6.6`).
* At `x = 0.05`:
  - `m = 39` is a clean B (`0.98`).
  - `m = 70` (the near-zero case asked about) has a small B term, `−0.0044`, because the Bessel-type factor `J₁(b)/J₀(b)` changes sign near
    `b = b̂ν ≈ 3.6`. Its measured value is twice that, with the right sign. A and B differ there by only `0.0044`, so it discriminates weakly
    (`m²d_A = −41`, `m²d_B = −19`).
  - `m = 101` is unusable, because `P_J(1−2x²)` is near a zero.
* **Verdict [NE]: reading B is supported for gravitons, and A and B̃ are disfavoured, in the double ratio.**
* **Caveats:**
  - one `b̂`, three masses, `x ≤ 0.1`;
  - the fraction of the B term is `0.63 → 0.93` at `x = 0.1`, not yet converged;
  - `ν = 7–10` lies beyond the quadrature checks of the e2e grid;
  - the double ratio assumes that the same-helicity moment follows `P_J`, which all readings share.

## 4. What is established, and what is not

* **[DVE]** (end to end, exact blocks, no fitted constant, robust extraction):
  - the factorised `⟨TTTT⟩` formula with the `l = 2` vertices for `β = 1`, `k = 4–7`, `ν = 0, 1, 3`, both helicities, under the odd kernel
    `𝔄⊗𝔄` and the leading modes;
  - the approach `1 + c₂/m²` with `|c₁| ≤ 0.04` over three masses (`m = 39–101`);
  - same helicity also at `k = 2, 3` (odd kernel);
  - **(§6) both helicities under the correctly normalised even kernels `𝔖_ff⊗𝔖_ff` and `K*⊗K*`, `k = 4–7`, `ν = 0, 1, 3`, with the same `c₂`
    as `𝔄⊗𝔄`** (`𝔖_ff` at `k = 4`, `ν = 3`, opposite helicity drifts, `c₁ = 1.9`);
  - the flip part of the formula through the mixed wave under `D₁²` (§6.2).
* **[NE]**
  - `β = 0.4`, `m = 72`, one mass (same helicity; opposite helicity under `𝔄⊗𝔄` and `K*`);
  - `E = 8` on a second ray;
  - the Fock closed form to `N = 8` [ID];
  - H2: reading B in the double ratio.
* **Open:**
  - ~~opposite helicity under the even kernels; `𝔖_ff` at `ν = 3` with `c₂ ≈ 1700`~~: artefacts of the naive kernel weights, resolved in
    §6;
  - `k = 2` (`C⁺₂`): resolved in §7 for `K*` and `𝔄` [NE]. Still open: `𝔖_ff` opposite helicity (isotropic-mode noise in the caches),
    `ν = 3` at `m = 39`, and `β = 0.4` opposite helicity under `K*`;
  - `β = 0.4` convergence (the `m ≈ 42` run was lost twice);
  - `ν ≳ 3` in general, and the `ν ~ m` range the functionals need;
  - the `1/4` constant end to end (blocks use unit-normalised structures);
  - H2 at other `b̂`;
  - parity-odd exchanges.
* **Consequence for `RESULT_D.md`.**
  - H1 is now [DVE] for the odd (`C⁻`) and even (`C⁺`, `𝔖_ff` or `K*`) moments, both helicities, at `k = 4–7`, `ν ≤ 3`, `m = 39–101`,
    `β = 1` (§6). It is not established at `k = 2` or at `ν ~ m`.
  - Corrections are `c₂/m²` with `|c₂| ≤ 310`. This is `3 %` at `m = Λ = 100` and `≤ 0.3 %` for `Δ_g ≥ 300` (§6.6).
  - H2 is supported for gravitons at one `b̂` [NE].

## 5. Code

| file | content | output |
|---|---|---|
| `tmoments_l2.py` | Taylor-mode `⟨TTTT⟩` waves (blocks one at a time, precomputed radial sums) | `cache/modes_l2_*.npz`, `cache/prod_*.log` |
| `t16_kernels_ext.py` | kernel moments `E⊗E`, `FF⊗FF` (**naive weights, superseded by `t18`**), `O⊗O` and leading modes vs prediction, `k = 2…7`, `ν = 0, 1, 3`; convergence in `m` | `out_t16_kernels_ext.txt` |
| `t19_k2.py` | §7: `k = 2` diagnosis (noise floor, cut scan) and cleaned extraction (strict mask, row tails, light-cone wedge), variants and envelope; `k = 4, 6` validation; `β = 0.4` | `out_t19_k2.txt`, `cache/t19b.log` |
| `t18_kstar.py` | §6: flip identity, tensor weights vs Lemma L, per-mode diagnosis, `𝔖_ff` / `K*` / `D₁²` / naive kernels both helicities, `k = 2`, mixed wave, `β = 0.4` | `out_t18_kstar.txt`, `cache/t18.log` |
| `t17_readingB.py` | reading A / B / B̃ at fixed `b̂`, single and double (opp/same) ratios | `out_t17_readingB.txt` |
| `t7_fock_l2.py l2 8 7`, `t7b` (ext) | Fock data for `[TT]` at `N = 7, 8` and the closed-form check | `out_t7_fock_l2ext.txt`, `out_t7b_closed_forms_ext.txt` |
| `t12_vertex_l2.py b0.9 1 2`, `t14_vertex_fit_l2.py` | `E` on a second ray | `out_t14_vertex_fit_l2.txt` |

## 6. Opposite helicity under even kernels: diagnosis, kernel normalisation, and `K*` (`t18_kstar.py`, `out_t18_kstar.txt`)

Added 18 Nov 2026 (same agent; not refereed; not committed). Existing caches only (`β = 1`, `m = 39.3/70.0/100.7`; `β = 0.4`, `m = 72.0`).

### 6.1 Result in one paragraph

**The even-kernel failures of §0–§2 were an error in my kernel weights, not physics.** `t16_kernels_ext.py` contracted the circle modes
with the *pointwise* Fourier coefficients of `D₁²` and `𝔖_ff` on the polarisation circles. The tensor kernels of the project instead carry
an extra factor `n_|q| = (1/36, 1/16, 1)` for `|q| = 0, 1, 2` per pair (§6.2). The naive weights over-weighted the isotropic `(0,0)` mode
by `36` per pair (`1296` for `(0,0,0,0)`), and that mode is contaminated by the dropped tiny-`t` region (§6.3). They also made "FF" not
flip-free (§6.2).
* **With the correct weights, `𝔖_ff⊗𝔖_ff` passes H1 for both helicities [DVE]** at `k = 4–7`, `ν = 0, 1, 3`, `β = 1`. The approach is
  `1 + c₂/m²` with the same `c₂` as `𝔄⊗𝔄` (§6.4).
* The one exception is `k = 4`, `ν = 3`, opposite helicity, which drifts (`c₁ = 1.9`). At `β = 0.4`, `m = 72` it fails at `k = 4`
  (`0.20–0.89`).
* A second flip-free even kernel, **`K* = D₁² + ½D₁D₂ − ½D₂² = (D₁+D₂)(D₁−½D₂)`**, has no isotropic mode. It is clean everywhere tested at
  `k ≥ 4`: `|c₁| ≤ 0.05`, and `β = 0.4` opposite helicity is `1.002–1.096`.
* **Its preserved response factorises with no real zero [D]:** `∝ (4ν²+(2a−9)²)(4ν²+(4a−7)²)`.
* Its graviton-pole residue and the absence of an Einstein term are the same as for `𝔖_ff`, by the same argument (§6.5).
* **Recommendation.** Keep `𝔖_ff` where it is clean (`k ≥ 5`, or `ν ≤ 1`). Switching to `K*` is the more robust choice. It changes only
  `N^{ff} → N^{K*}` (closed forms below).
* **`k = 2` (`C⁺₂`) is still not resolved** under any even kernel. The obstruction is the tiny-`t` region, not the kernel.

### 6.2 Kernel normalisation and the flip identity

* **[D, checked to `7·10⁻¹⁵`] Flip identity.** The one-pair flip responses of `laplace.M_table(2)` per circle mode `(q,−q)` are in the
  ratio `(0,0) : (±1,∓1) : (±2,∓2) = 36 : 16 : 1`, for every `ν` and `a` (grid `ν = 0…10`, `a = 6…12`). These numbers are
  `C(4,2+q)²`.
* **[ID → DVE] Tensor weights.** Write a pointwise kernel on the circles as `K(θ) = Σ_q w̃_q e^{iq(θ_out−θ_in)}`. The tensor kernel acting
  on the modes has weights `w_q = n_|q| w̃_q`, with `n = 1/C(4,2+q)² = (1/36, 1/16, 1)`.
  * The two ratios of `n` were fixed by the two Lemma-L flip ratios (`D₁² : D₁D₂ : D₂² = 1 : −1 : 1`, `REGGE_MOMENT_TTTT.md` §4.2).
  * They then reproduce all 18 Lemma-L preserved ratios to `10⁻⁶` (`ν = 0, 1.1, 3`; `a = 7, 9, 11`), for example
    `D₁D₂/D₁² = −0.270947` and `D₂²/D₁² = +0.129489` at `ν = 0`, `a = 7`.
  * The resulting weights per pair are:

  | kernel | `(0,0)` | `(±1,∓1)` | `(±2,∓2)` | flip |
  |---|---|---|---|---|
  | `D₁²` | `1/72` | – | `1/4` | `≠ 0` |
  | `D₁D₂` | – | `−1/32` | – | `= −D₁²` |
  | `D₂²` | `1/36` | – | – | `= D₁²` |
  | `𝔖_ff = D₁² + D₁D₂` | `1/72` | `−1/32` | `1/4` | `0` |
  | `K* = D₁² + ½D₁D₂ − ½D₂²` | `0` | `−1/64` | `1/4` | `0` |
  | t16 "E" (naive) | `1/2` | – | `1/4` | `≠ 0` |
  | t16 "FF" (naive) | `1/2` | `−1/2` | `1/4` | `≠ 0` (`|flip|/pres` up to `0.3`) |

  This is the representation caveat noted earlier (the circle representation did not reproduce Lemma L). It is now resolved.
* **[D] Flip-free criterion.** With these weights, flip-free means `36w₀ + 16(w₁+w₋₁) + (w₂+w₋₂) = 0`, which is `Σ_q w̃_q = K(θ_out = θ_in) =
  0`. **A kernel is flip-free if and only if its circle restriction vanishes at aligned polarisations.**
  * `𝔖_ff = cos θ(cos θ − 1)` and `K* = (cos θ − 1)(cos θ + ½)` both vanish there.
  * `K*` additionally has zero circle mean, so no isotropic mode. It is the unique such combination of `D₁², D₁D₂, D₂²` up to
    normalisation: the two conditions are linear, and the space is 3-dimensional.
* **Mixed wave (check of the flip prediction).** For flip-free kernels the predicted `W[X_s,X_o]` moment is zero, to `10⁻¹⁸–10⁻²¹`
  (roundoff).
  * Measured relative size: `K*` `1.6·10⁻⁵ → 1.3·10⁻⁷` (`k = 4`, `ν = 0`, `m = 39 → 101`), decreasing with `m`; `𝔖_ff` `10⁻⁵–10⁻¹⁰`,
    erratic (the isotropic mode).
  * For `D₁²`, which is not flip-free, the measured mixed wave converges to the predicted flip: `1.9·10⁻²/1.8·10⁻²` at `k = 6`, `ν = 3`, at
    all three masses. This tests the flip part of the `l = 2` formula end to end.

### 6.3 Diagnosis (route 1)

* **Per mode [NE, `t18` (a)].** Every mode with a non-negligible prediction converges as `1/m²`, with flat `m²(R−1)`.
  * Same helicity: all four main modes at `ν ≤ 1` have `m²(R−1) ≈ 70–200`.
  * Opposite helicity: the leading `(2,−2,2,−2)` mode equals the `𝔄⊗𝔄` values of §2.
* **The isotropic `(0,0,0,0)` mode of `W[X_o,X_o]` is not under control.** Its prediction is tiny (a cancellation inside the opposite-helicity
  vertex). Its measured value is dominated by the dropped tiny-`t` region (radial-series mask, `N = 8000`):
  * the boundary share is `2.6×`, `11 %` and `0.9 %` of the mode at `m = 39, 70, 101` (`k = 6`, `ν = 0`);
  * internal cancellation is `22×` at `m = 39`;
  * `m²(R−1)` is erratic and changes sign: `−2416, +582, +271` (`k = 6`, `ν = 0`) and up to `10⁷` at `ν = 3`.

  `(1,−1,1,−1)` and `(0,0,2,−2)` are partly affected at `k = 4`. **This is an extraction (completeness) artefact, not a physical
  subleading term [NE].**
* **Why the naive kernels failed.**
  * Weight `1/2` instead of `1/72` on `(0,0)` gives a `1296×` over-weighting of `(0,0,0,0)` in the opposite-helicity wave.
  * For naive FF, the prediction is also a `37–100×` cancellation, which amplifies the ordinary per-mode `c₂ ≈ 300` (`ν = 3`) to the
    `m²(R−1) ≈ 1500–1800` reported in §0/§2 for **same** helicity at `ν = 3`. With the correct weights, `𝔖_ff`'s prediction cancels only
    `1.2–2×`, and same helicity at `ν = 3` gives `m²(R−1) = 73/69/68` (`k = 6`) and `−6/−8/−9` (`k = 7`).
* **Conclusion of route 1.** The leading heavy formula is correct in every mode that the data resolve. There is no `O(1/m)` term: `|c₁| ≤ 0.05`
  for `𝔄`, `K*`, and `𝔖_ff` except `k = 4`, `ν = 3`, opposite helicity.
  * The ansatz `1 + c₁/m + c₂/m² + c₄/m⁴` fits with `c₁ ≈ 0`.
  * A `log m/m²` term of `≲ 10 %` of `c₂` cannot be excluded with three masses. It is equivalent to the `5–10 %` drift of `m²(R−1)` seen in
    some entries.

### 6.4 End-to-end tables (`β = 1`, `R` at `m = 39.3 / 70.0 / 100.7`; `m²(R−1)`; `c₁`)

Opposite helicity:

| `k`, `ν` | `𝔖_ff⊗𝔖_ff` | `K*⊗K*` |
|---|---|---|
| 4, 0 | 1.0567 / 1.0163 / 1.0080; 88, 80, 81; +0.12 | 1.0569 / 1.0166 / 1.0078; 88, 81, 79; −0.01 |
| 4, 1 | 1.0787 / 1.0231 / 1.0113; 122, 114, 115; +0.13 | 1.0787 / 1.0235 / 1.0113; 122, 116, 114; +0.01 |
| 4, 3 | 0.9390 / 0.9837 / 0.9968; −94, −80, −32; **+1.88** | 0.9408 / 0.9835 / 0.9921; −92, −81, −80; −0.05 |
| 5, 0 | 1.0312 / 1.0085 / 1.0041; 48, 42, 41; +0.04 | 1.0313 / 1.0085 / 1.0040; 48, 42, 40; −0.00 |
| 5, 1 | 1.0557 / 1.0163 / 1.0079; 86, 80, 80; +0.04 | 1.0559 / 1.0164 / 1.0078; 86, 80, 79; +0.00 |
| 5, 3 | 0.9006 / 0.9707 / 0.9862; −154, −144, −140; +0.07 | 0.9017 / 0.9710 / 0.9862; −152, −142, −140; −0.01 |
| 6, 0 | 1.0002 / 0.9985 / 0.9991; 0, −7, −9; +0.00 | 1.0003 / 0.9985 / 0.9991; 1, −7, −9; −0.00 |
| 6, 1 | 1.0269 / 1.0073 / 1.0034; 42, 36, 34; +0.00 | 1.0271 / 1.0073 / 1.0034; 42, 36, 34; +0.00 |
| 6, 3 | 0.8547 / 0.9554 / 0.9787; −225, −219, −217; +0.02 | 0.8556 / 0.9558 / 0.9788; −223, −217, −216; −0.00 |
| 7, 0 | 0.9634 / 0.9862 / 0.9930; −57, −68, −71; −0.01 | 0.9635 / 0.9862 / 0.9930; −56, −68, −71; −0.01 |
| 7, 1 | 0.9921 / 0.9960 / 0.9979; −12, −20, −22; −0.00 | 0.9923 / 0.9960 / 0.9979; −12, −19, −21; −0.00 |
| 7, 3 | 0.8019 / 0.9372 / 0.9696; −307, −308, −308; +0.00 | 0.8028 / 0.9375 / 0.9697; −305, −307, −307; −0.00 |

Same helicity:
* `𝔖_ff` and `K*` agree with each other and with `𝔄⊗𝔄` (§2) to `≤ 3·10⁻³`, at all `k = 4–7` and `ν = 0, 1, 3`. For example `k = 6`, `ν = 3`:
  `1.0469 / 1.0140 / 1.0067` (`𝔖_ff`) and `1.0462 / 1.0137 / 1.0065` (`K*`).
* `|c₁| ≤ 0.02`.
* The opposite-helicity `K*` values also agree with the `𝔄⊗𝔄` column of §2 to `≤ 2·10⁻³` (`𝔖_ff` too, except `k = 4`, `ν = 3` at
  `m = 101`). The helicity-preserving content is kernel-independent, as the factorised formula requires.

**`β = 0.4`, `m = 72.0`** (one mass; `R`):

| `k` | same (all three, `ν = 0, 1, 3`) | opp `𝔖_ff`, `ν = 0 / 1 / 3` | opp `K*`, `ν = 0 / 1 / 3` |
|---|---|---|---|
| 4 | 1.026–1.031 | 0.885 / 0.858 / **0.196** | 1.010 / 1.022 / 1.096 |
| 5 | 1.023–1.031 | 1.007 / 1.012 / 1.042 | 1.002 / 1.007 / 1.048 |
| 6 | 1.019–1.029 | 0.995 / 1.000 / 1.043 | 0.995 / 0.999 / 1.042 |
| 7 | 1.013–1.025 | 0.986 / 0.991 / 1.035 | 0.986 / 0.991 / 1.035 |

There is no convergence statement at `β = 0.4`, because the `m ≈ 42` point is missing. `K*` removes the `k = 4` failure of `𝔖_ff`.

**`k = 2` (`C⁺₂`, `a = 7`).**
* Same helicity: `1.02–1.12 / 1.05–1.09 / 1.02–1.03` at `ν ≤ 1` for all three kernels, but `c₁ ≈ −6` to `−12`. `m²(R−1)` is not flat.
* `ν = 3` and opposite helicity: not resolved (`R` is `O(1)` off, or worse, for every kernel). `K*` is the least bad (`1.94 / 1.13` at
  `m = 70 / 101`, `ν = 0`).
* The same holds for the odd kernel at `k = 2` (§2). This is the tiny-`t` region, which at low `k` holds a large share of every mode.
* **H1 at `k = 2` with opposite helicity remains open [ND].** It needs a larger `N` or a tiny-`t` completion, not a different kernel.

### 6.5 `K*` as a functional kernel (route 2)

* **Flip-free [D]:** Lemma-L flips `1 − ½ − ½ = 0`.
* **Preserved response [D, from Lemma L]:** `∝ (4ν²+(2a−9)²)(4ν²+(4a−7)²)`. Compare `𝔖_ff`, which is `∝ 2(4ν²+(2a−9)²)(4ν²+8(a−2)²+1)`. It
  is positive on ℝ for all `a` used (`a = 7, 9, 11`; no real zero for `a ≠ 9/2, 7/4`).
  * `K*/𝔖_ff` preserved = `1.097, 1.006, 0.70, → ½` at `ν = 0, 3, 10, ∞` (`a = 7`).
  * The `t18` table check gives `0.0851/0.0775 = 1.098` at `ν = 0`, `a = 7`, which cross-checks the weights.
* **Normalisations [D given Lemma L and `N^{ff}`].** These replace `N^{ff}` (the `(1,3)` ratio `j′=0 : kernel`):
  * `N₄^{K*} = 8(16ν⁴+216ν²+263)/((4ν²+81)(4ν²+841))`;
  * `N₆^{K*} = 8(16ν⁴+536ν²+2023)/((4ν²+169)(4ν²+1369))`.

  Both are positive on ℝ. `N^{K*}(0) = 0.031, 0.070` against `N^{ff}(0) = 0.033, 0.074`, so the `1/N(0)` leverage on the forward piece rises
  by `7 %`. `N^{K*}(∞) = 8`.
* **Graviton pole (I9) [D given IP, (F)], same status as for `𝔖_ff`.** The Einstein exchange enters the even moment only through its
  helicity-preserving density. For a flip-free kernel the flip parts are zero identically. The preserved part is the kernel's preserved
  response times the kernel-independent density, and dividing by `N_k` removes the kernel's response. So the normalised residue
  `16πG/(ν²+9/4)` is unchanged.
  * This needs a preserved response that is nonzero at the pole. That holds, by the factorisation above.
  * The end-to-end statement that the helicity-preserving content is kernel-independent (§6.4: `K*`, `𝔖_ff`, `𝔄` agree to `≤ 3·10⁻³`) is
    the heavy-state analogue [NE].
* **Einstein absence, `Ψ_k[G^E] = 0` for `k ≥ 3` [D given RB₂].** `SUPERCONVERGENT_UPLIFT.md` item 4 uses only microcausality of `G^E`,
  `J₀(G^E) = 2 < 3`, and the functional's convergence reach. `K*` has the same polynomial degree in `D₁, D₂` per pair as `𝔖_ff`, so the same
  Regge growth and the same reach [D-sketch]. Dividing by `N_k^{K*} > 0` creates nothing from zero.
* **Not changed by a switch:**
  * (P), the double traces and contacts of the even and forward pieces, stays [I] and load-bearing;
  * the `C⁺₂` convergence margin `J₀ < 2`.

### 6.6 Consequence for the functionals (route 3)

* **The leading heavy formula is correct for the even (`𝔖_ff`, `K*`) and odd moments, both helicities [DVE: `k = 4–7`, `ν ≤ 3`, `β = 1`,
  `m = 39–101`].**
  * Corrections are `c₂/m²` with no `1/m` term.
  * Over the tested range, `|c₂| ≤ 310`. The extremes are `−307` (`k = 7`, `ν = 3`, opposite helicity) and up to `+194` for same helicity.
  * The sign is negative for opposite helicity at `ν = 3`, so it lowers the heavy contribution.
* **Relative effect on a heavy state's contribution at `m ≥ Λ` [NE, extrapolated in `m`; `ν ≤ 3`]:**

  | `Δ_g` | `310/m²` at `m = Λ = Δ_g` | at `m = 1.4Λ` (the `OMEGA_U` binding state) |
  |---|---|---|
  | 100 | `3.1·10⁻²` | `1.6·10⁻²` |
  | 300 | `3.4·10⁻³` | `1.8·10⁻³` |
  | 1000 | `3.1·10⁻⁴` | `1.6·10⁻⁴` |

  These take `m ≈ Λ` in units of `Δ_g`, as in `OMEGA_U`. At the data points `m/Δ = 0.86–0.87`, so the effect in units of `Δ` is about 1.3×
  larger.
* **Size against the margins.**
  * At `Δ_g = 100` the correction is of the same order as the spare margin at the binding state, quoted as `0.044 N_s`
    (`SUPERCONVERGENT_UPLIFT.md` item 6), if that state's heavy action is `O(N_s)`. It should be kept explicitly, or `Δ_g ≥ 300` used.
  * At `Δ_g ≥ 300` it is `≤ 0.4 %` and negligible against margins `≥ 4 %`.
* **Caveats.**
  * This bounds only fixed `ν ≤ 3`. The functionals sample `ν ~ Λ` (`ν/m = O(1)`). There the expansion parameter is `x = ν/m`, not `1/m`,
    and the relevant check is H2 (§3, reading B at fixed `b̂`).
  * `c₂(ν)` is not monotone over `ν = 0, 1, 3`: `−56, −12, −305` (`k = 7`, opposite helicity), so `c₂` cannot be extrapolated in `ν` from
    these three points.
  * `β ≠ 1` has one mass only.
  * `k = 2` is open (§6.4).

## 7. `k = 2` (`Ĉ₂`, the graviton-pole moment): diagnosis and a cleaned extraction (`t19_k2.py`, `out_t19_k2.txt`)

Added 18 Nov 2026 (not refereed; not committed). Existing caches only. Kernel weights as in §6.2.

### 7.1 Result

**H1 at `k = 2` is now resolved under `K*⊗K*` and `𝔄⊗𝔄` for both helicities [NE, β = 1, `m = 39–101`].**
* Cleaned ratios: same helicity `m²(R−1) ≈ 196–268`; opposite helicity `≈ 140–210` (`ν ≤ 1`).
* `|c₁| ≤ 1.3` at `ν ≤ 1`, and `≈ 2` at `ν = 3`, same helicity.
* At `ν = 3`, opposite helicity, `R = 1.009 ± 0.010` and `1.002 ± 0.006` (`K*`) at `m = 70, 101`. The `m = 39` point is not resolved.

**Under `𝔖_ff⊗𝔖_ff` the opposite-helicity `k = 2` moment is still not resolved.** Its isotropic `(0,0,0,0)` mode carries noise in the
cached waves that the cleaning cannot remove (§7.2). This is a second, independent reason to use `K*` for `Ĉ₂`.

**Diagnosis, answering the three alternatives.** The failure was mostly numerical defects of the cached waves, not physics:
* The radial series is truncated next to the edge of the region where it converges.
* There is float-level noise at tiny `t`.
* Both are amplified by the `k = 2` weight and were then extrapolated by the old corner estimator.

The three alternatives:
* **The dropped region (radial series not converged) is real but small [NE].** It is a light-cone wedge `ū ≲ 0.01–0.04` (small `ρ̄`) plus
  the corner `t ≲ 1`. Its contribution, extrapolated, is `+0.6 %` to `+6 %` (`K*`, opposite helicity, `m = 101 → 39`).
* **There is no light-cone singularity [NE].** At fixed `t` the `k = 2` density vanishes towards the light cone like `u^{≈1.7}` (the
  `F`-stripped density `∝ u^{1.16}` at `m = 101`). The `k = 2` moment of a single heavy block converges there.
* **There is no slower approach [NE].** Once cleaned, `m²(R−1)` is flat to `≈ 10 %` over `m = 39–101`. That is the same `1/m²` approach
  as `k ≥ 4`, with `c₂ ≈ 200`.

### 7.2 What goes wrong in the raw `k = 2` moment

The moment weight is `(ρρ̄)^{(k−4)/2} = S^{2k−8}u^{k−4}`. Relative to `k = 4`, `k = 2` enhances:
* the tiny-`t` corner (`S → 0`) by `S^{−4}`;
* the light-cone wedge (`u → 0`) by `u^{−2}`.

Two defects of the cached waves sit exactly there. Both are invisible at `k ≥ 4`, where the cleaning changes nothing (§7.3).

* **(D1) Truncation next to the edge of convergence.** The production mask keeps nodes with `X − 10 log X > 25`. That cut was calibrated
  on the `l = 1` truncation error (`e2e/moments.py`). For the `l = 2` opposite-helicity modes, the first accepted nodes next to the edge are
  wrong by up to `10⁴×` the smooth density.
  * Example: `m = 39`, `t = 9.3`, `u = 0.013`: density `1.5·10⁵` against `≈ 5` a few nodes further in.
  * The error is real-valued and does not appear in the imaginary parts.
  * **Fix:** raise the cut to `40`. This drops `2.6–3.7 %` of the nodes (`3336/3475/3272 → 3250/3354/3150`). The results agree with cut `50`
    to `≤ 2·10⁻³` at `m ≥ 70`.
* **(D2) Noise at tiny `t`.** The modes used are real by symmetry, so their imaginary parts are pure noise.
  * The median `|Im G_q|` over the 624 complex modes gives a pointwise noise floor σ. At `t ≲ 3`, σ is `3·10⁻⁷–4·10⁻⁵` in mode units,
    comparable to the leading-mode signal there.
  * The opposite-helicity isotropic mode `(0,0,0,0)` is real by construction, so it has no such proxy. It shows sign-alternating values of
    `0.4–40` at `t ≲ 5`, `10⁵×` the leading mode at the same nodes.
  * This noise is what the earlier "tiny-`t` contamination" was.
  * **Fix:** per `u`-row, start where three consecutive nodes exceed `10σ` (propagated through the kernel) and follow a smooth power law in
    `t`. Before that start, use a power law fitted to the next three nodes. The true density there rises like `t^{4–10}`, so the replaced
    region holds `≲ 0.5 %`.
* **Light-cone wedge.** In each `t`-column (`t ≥ 3`), the dropped nodes `u < u_first` are filled with `H(u) = d/F`, extrapolated as a power
  law in `u` from the first three clean nodes and multiplied by the exact kernel `F(ξ(u))`. This keeps the `u^{±iν}` oscillation of `F`.

The old corner estimator extrapolated from the first accepted nodes, which carry the (D1) and (D2) errors. That is why it made `k = 2`
worse; at `k ≥ 4` the weight suppresses those nodes.

**Variants and envelope.**
* The model spread is the largest deviation over: fit nodes shifted by one, two-node fits, noise thresholds `3σ` and `30σ`, wedge from
  `t ≥ 5`, and cut `50`.
* The envelope adds constant (slope-free) tails and the original production cut `32`.

### 7.3 Results (`β = 1`; `R` at `m = 39.3 / 70.0 / 100.7` [model spread at `m = 101`]; `m²(R−1)`; `c₁`)

| `k = 2`, `ν` | helicity | `K*⊗K*` | `𝔄⊗𝔄` | `𝔖_ff⊗𝔖_ff` |
|---|---|---|---|---|
| 0 | same | 1.1338 / 1.0399 / 1.0208 [0.004]; 207, 196, 211; +0.8 | 1.1336 / 1.0401 / 1.0216; 207, 197, 219; +1.0 | 1.1343 / 1.0400 / 1.0204; 208, 196, 207; +0.6 |
| 0 | opp | 1.1141 / 1.0295 / 1.0141 [0.001]; 177, 145, 143; +0.3 | 1.1357 / 1.0279 / 1.0148; 210, 137, 151; +1.3 | 1.09 / 0.92 / 2.03: **not resolved** |
| 1 | same | 1.1343 / 1.0400 / 1.0197 [0.003]; 208, 196, 200; +0.3 | 1.1339 / 1.0397 / 1.0185; 207, 195, 188; −0.2 | 1.1347 / 1.0402 / 1.0195; 208, 197, 198; +0.1 |
| 1 | opp | 1.1109 / 1.0322 / 1.0157 [0.001]; 172, 158, 160; +0.2 | 1.1126 / 1.0331 / 1.0156; 174, 162, 159; −0.1 | 1.09 / 0.88 / 2.42: **not resolved** |
| 3 | same | 1.1515 / 1.0453 / 1.0264 [0.036]; 234, 222, 268; +2.1 | 1.1490 / 1.0444 / 1.0263; 231, 218, 267; +2.2 | 1.1531 / 1.0436 / 1.0249; 237, 214, 252; +1.9 |
| 3 | opp | 0.958 [0.47] / 1.0088 [0.010] / 1.0024 [0.006] | 0.882 [0.12] / 1.0019 / 1.0039 | not resolved |

Envelopes (min..max including constant tails and cut `32`), `K*`:
* same, `ν = 0`: `1.134..1.135 / 1.040..1.047 / 1.020..1.057`;
* opposite, `ν = 0`: `1.100..1.510 / 1.027..1.036 / 1.014..1.036`;
* opposite, `ν = 3`: `0.17..1.00 / 0.999..1.014 / 1.000..1.008`.

The constant-tail variant is a deliberately crude upper envelope: it ignores the measured slope `u^{1.7}`.

* **Validation at `k = 4, 6` [NE].** The same cleaning reproduces the `K*`, `𝔄` and leading-mode values of §6.4 to `≤ 10⁻³` (model spread
  `≤ 3·10⁻⁴` at `m ≥ 70`).
  * `𝔖_ff` with cleaning degrades at `k = 4`, `m = 101` (`1.127` vs raw `1.008`), because the isotropic-mode noise defeats the noise
    criterion.
  * For `𝔖_ff`, use the raw §6.4 values at `k ≥ 4`.
* **`k = 2` against the earlier §2 numbers.**
  * The leading-mode / `𝔄` same-helicity values `1.099 / 1.030 / 1.012` used the old corner estimator, biased by (D1). Cleaned:
    `1.134 / 1.040 / 1.021`.
  * The opposite-helicity entry "not resolved (tiny-`t`)" is now resolved for `𝔄` and `K*` at `ν ≤ 1` (all three masses) and at `ν = 3`
    (`m ≥ 70`).
* **`β = 0.4`, `m = 72`** (one mass; the raw values are `10²–10⁴` off):
  * same helicity: `1.028` for `K*`, `𝔄` and the leading mode, at `ν = 0, 1, 3`;
  * opposite helicity: `𝔄` `1.018 / 1.020 / 1.059` (spread `≤ 0.007`); `K*` `1.022` at `ν = 1` (spread `0.01`), but not resolved at `ν = 0`
    and `ν = 3` (spread `1.05`, `0.46`);
  * `𝔖_ff`: not resolved.

### 7.4 Status

* **H1 for `Ĉ₂` [NE].** Holds under `K*` (both helicities, `ν ≤ 1`, `m = 39–101`; `ν = 3` at `m ≥ 70`, to 1 %) and under `𝔄`. The
  approach is `1 + c₂/m²` with `c₂ ≈ 140–270`. Not established: `𝔖_ff` at opposite helicity (cache defect), `β ≠ 1` convergence, and
  `ν ~ m`.
* **The cleaning is not a proof [I].** Its envelope rests on two things:
  * smoothness: power-law behaviour beyond the clean region in `t` (towards 0) and in `u` (towards the light cone);
  * the noise proxy: `Im G` as the noise of the real parts.

  A cleaner route is a production rerun with a stricter mask and long-double sums near the edge. The radial-series edge would then need a
  larger `N`, more than 1 h per mass, so it was not done.
* **Implication for the functionals.** `Ĉ₂` should use `K*` (or the odd-kernel cross-check). `𝔖_ff`'s `k = 2` opposite-helicity moment
  cannot currently be tested end to end.
  * Under the normalised convention, the switch leaves the functionals unchanged at leading order (`EXACT_SKELETON_NEARGAP.md` §9).
  * The measured `c₂` values feed the finite-`m` sensitivity there.
