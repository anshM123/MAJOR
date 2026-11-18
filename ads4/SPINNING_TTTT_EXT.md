# `⟨TTTT⟩` beyond one point: H1 (heavy formula at more `k, ν, β, m`) and H2 (reading B for gravitons)

Status: written 18 Nov 2026 by a derivation/numerics agent; **not refereed; not committed**. Code: `ads4/code/fourpt/tttt/`
(`tmoments_l2.py` production, `t16_kernels_ext.py`, `t17_readingB.py`; logs and caches in `tttt/cache/`, git-ignored). Builds on
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
   * **Even kernels.** Same helicity behaves like `𝔄⊗𝔄` under `D₁²⊗D₁²` (`k = 4–7`, `ν ≤ 1`). Under `𝔖_ff⊗𝔖_ff` it converges like
     `1/m²` but with a large coefficient at `ν = 3`: `m²(R−1) ≈ 1500–1800`, so still `+15 %` at `m = 101`. The flip-free kernel's preserved
     response is a strong cancellation (`SUPERCONVERGENT_UPLIFT.md`: `(4ν²+201)/(2(4ν²+81))`).
   * **Opposite helicity under either even kernel is not resolved** at `m ≤ 101`, except `k = 7`, `ν ≤ 1`, which is within `1–2 %`. The
     isotropic `(0,0)` and `(±1,∓1)` modes, which the even kernels include and `𝔄` does not, carry finite-`m` contributions that the leading
     formula does not describe. They are erratic in `m`.
   * The flip-free kernel does not uniformly suppress the mixed wave `W[X_s,X_o]`: it is 6× smaller than under `D₁²` at `m = 39`, but larger at
     `m = 101`. **H1 for the even (`C⁺`) kernels at `ν ≳ 3` and for opposite helicity is therefore not established.**
   * **`β = 0.4`, `m = 72`** (single mass, `J = 14`):
     - same helicity `1.013–1.030` (all kernels, `k = 4–7`);
     - opposite helicity, odd kernel `0.93–1.05` (`0.986–1.007` at `ν ≤ 1`, `k ≥ 5`);
     - opposite helicity, even kernels not resolved (Weyl-odd terms large, `λ_o/λ_s = 2.9`).

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

These are used as `E⊗E`, `FF⊗FF` and `O⊗O`, together with the leading single mode of each wave. The prediction is that of
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

**Even kernels, `β = 1`** (`R` at `m = 39 / 70 / 101`):

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
* `ν = 3`: `1.024–1.030`, except `𝔖_ff`, which is `1.27–1.40`;
* opposite helicity, `𝔄⊗𝔄`: `0.935–1.048` (`0.986–1.007` at `ν ≤ 1`, `k ≥ 5`);
* opposite helicity, even kernels: not resolved.

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
  - same helicity also at `k = 2, 3` and under both even kernels (`D₁²`, and `𝔖_ff` at `ν ≤ 1`).
* **[NE]**
  - `β = 0.4`, `m = 72`, one mass (same helicity, opposite helicity under `𝔄⊗𝔄`);
  - `E = 8` on a second ray;
  - the Fock closed form to `N = 8` [ID];
  - H2: reading B in the double ratio.
* **Open:**
  - opposite helicity under the even kernels (`D₁²`, `𝔖_ff`) at `m ≤ 101`, except `k = 7`, `ν ≤ 1`;
  - `𝔖_ff` at `ν = 3`, where `c₂ ≈ 1700` makes the leading formula `15 %` off at `m = 101`;
  - `k = 2`, opposite helicity (dropped tiny-`t` region);
  - `β = 0.4` convergence (the `m ≈ 42` run was lost twice);
  - `ν ≳ 3` in general, and the `ν ~ m` range the functionals need;
  - the `1/4` constant end to end (blocks use unit-normalised structures);
  - H2 at other `b̂`;
  - parity-odd exchanges.
* **Consequence for `RESULT_D.md`.**
  - H1 is now [DVE] for the odd-kernel (`C⁻`) and same-helicity parts at `ν ≤ 3`, `m = 39–101`, `β = 1`, but not for the
    opposite-helicity even-kernel (`C⁺`, `𝔖_ff`) moments the functionals use, nor at `ν ~ m`.
  - H2 is supported for gravitons at one `b̂` [NE].

## 5. Code

| file | content | output |
|---|---|---|
| `tmoments_l2.py` | Taylor-mode `⟨TTTT⟩` waves (blocks one at a time, precomputed radial sums) | `cache/modes_l2_*.npz`, `cache/prod_*.log` |
| `t16_kernels_ext.py` | kernel moments `E⊗E`, `FF⊗FF`, `O⊗O` and leading modes vs prediction, `k = 2…7`, `ν = 0, 1, 3`; convergence in `m` | `out_t16_kernels_ext.txt` |
| `t17_readingB.py` | reading A / B / B̃ at fixed `b̂`, single and double (opp/same) ratios | `out_t17_readingB.txt` |
| `t7_fock_l2.py l2 8 7`, `t7b` (ext) | Fock data for `[TT]` at `N = 7, 8` and the closed-form check | `out_t7_fock_l2ext.txt`, `out_t7b_closed_forms_ext.txt` |
| `t12_vertex_l2.py b0.9 1 2`, `t14_vertex_fit_l2.py` | `E` on a second ray | `out_t14_vertex_fit_l2.txt` |
