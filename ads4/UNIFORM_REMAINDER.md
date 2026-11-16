# Uniform remainder near the gap: the AdS heavy action of the robust functionals from exact scalar-skeleton moments

Status: written 16 Nov 2026 by a derivation/numerics agent (STATUS link 9); **revised after `ads4/REFEREE_UNIFORM_REMAINDER.md`**
(changes in §8). **Not a theorem about AdS4.** Code: `ads4/code/uniform/` (`bash run_all.sh`; outputs in `uniform/cache/`, git-ignored;
it needs the CHMRS ancillary file `anc/regge_moments_order8.m` of arXiv:2106.10274 copied to `cache/`). It imports, unchanged,
`e2e/moments.py` (exact blocks), `robust/rcommon.py`, `robust/check.py`, `robust/ads.py` (functionals, flat parts, hybrid), and reads
`robust/functionals_certified.json`. Nothing in other directories is modified.

Marks: **[certified]** exact rational / interval arithmetic; **[D]** derived; **[DV]** derived and checked against exact blocks;
**[Q]** quoted from CHMRS (arXiv:2106.10274, LaTeX and ancillary file); **[NE]** floating-point numerics; **[M]** statement about a stated
model, cut or extrapolation; **[I]** inferred/assumed. CHMRS = arXiv:2106.10274.

## 0. Results

1. **Structure [D + DV].** The leading-order `⟨TTTT⟩` skeleton is the identical-scalar moment with `Δφ = 3`: `a = 2Δφ + k − 1 = k + 5`
   for both pairs (`SPINNING_E2E.md` §8). Its exact heavy action organises as
   **flat (bulk point) + leading-Regge curvature + remainder**:
   `m^{2k}C_{k,ν} = E_k(x)P_J(1−2x²) + 2[P_{−½+iν}(cosh β) − J₀((2J+1)x)] + ρ_kE_k(x)`, with `x = ν/m` and `E_k = (2−x²)/(1−x²)^{1+k/2}`.
   * The leading-ν part of the exact (CHMRS order-1/m⁸) fixed-J series **is** the bulk-point formula, at every order, symbolically in `k, Δφ, J`
     [DV, `t1`; re-derived by the referee]. This reproduces CHMRS (3.40) as the resummation of their own series, so it is not new. It means the
     `O(ν²/m²)` terms are the flat density, not an unknown correction.
   * The remainder after subtracting flat plus leading-Regge curvature has exact coefficients [Q + D]
     (`residual_coeffs.py`; the curvature part is subtracted in closed form). The subtraction lowers the `𝒥`-degree of each `x^{2n}`
     coefficient by one. The remaining `(𝒥x²)^{n}` terms are the hybrid's `O(uβ²)` cross terms.
   * The CHMRS series is validated against exact blocks: exact/series = `1 ± 10⁻⁸` at `m ≈ 67` (`Δφ = 3, 1.75`; `k = 4, 6`), and the mutation `a → a − ½` fails
     [DV, `v1`]. In the bulk-point regime, `m²(exact/BP − 1) = −4 … −80` (`J = 0`, `x ≤ 0.5`, `k = 4, 6`), while exact/Regge-leading deviates by up to +98 %
     [DV, `t2`]. The Mellin `J = 0` check agrees to `4·10⁻⁸ – 2.5·10⁻⁵` at `m = 29`. At `m = 58` it disagrees by `5·10⁻⁴`; this is an unresolved defect
     of `mellin.py` (n-sum or quadrature), so it is **not** an independent check beyond `m ≈ 30`.
2. **The skeleton is fixed only at leading order (referee issue 1).**
   * At `1/m²`, `SPINNING_MFT.md` / `REFEREE_SPINNING_MFT.md` found the per-l term clean only relative to the skeleton `(l, Δφ) → (2, Δφ)`.
     The remainder is not robust to this choice: `r̃₁(x = 0, J = 0)` is `+12.9` (`Δφ = 3`) versus `−19.1` (`Δφ = 2`), at `k = 2`. **The sign flips.**
   * Everything below is therefore given for **both** choices, `Δφ = 3` and `Δφ = 2`. Their difference is the skeleton-choice error.
3. **Convergence of the remainder series (referee issue 2) [estimate].** Only four coefficients in `x²` exist (three at `1/m⁴`).
   * A two-ratio Domb–Sykes analysis (`t4_radius.py`) is calibrated on `E_k`, whose true singularity is at `y_c = 1`; there it returns `y_c = 0.72–0.86`.
   * For `J = 0` the sign-regular ratios give `y_c ≈ 0.69–0.88` (`Δφ = 3`) and `0.59–0.66` (`Δφ = 2`). After the calibration bias this is consistent with
     a singularity at `ν = m` (`y_c ≈ 1`), or at `y_c ≈ 0.8` for `Δφ = 2`.
   * For `J ≥ 2` the ratios change sign and nothing can be inferred.
   * **"Divergent" is withdrawn.** The control criterion is a **heuristic cut** [M]: a node counts as "controlled" if the last ratio is `≤ 0.8`; the
     fixed-η fallback is used only for `x ≤ 0.6`. It is not a derived boundary.
4. **AdS action of M02 (M01, F_g3) [NE + M + I, `t3`]**: about 2800–4000 physical states per `Δ_g`, with `m ≤ 10Λ`, `b̂ ≤ 40`, both channels.
   * **Method.** The hybrid action `F_hyb` (flat plus leading-Regge curvature) minus the transferred skeleton remainder `E` (absolute-value bound on
     controlled nodes). `U` is the same weight on uncontrolled nodes. `κ_crit = min (F_hyb − E)/U`: positivity holds if the relative AdS correction on
     uncontrolled nodes is below `κ_crit`.
   * **Where the uncontrolled nodes are (referee issue 3).** Every state below `1.5Λ` has some (for `Δ_g ≤ 300`), and so do about 75 % of the states in
     `[1.5Λ, 3Λ)`, at `x < 0.6`, with `J` up to 58. These come from the unresummed `𝒥x²` terms.
   * **For `m ≥ 1.5Λ` positivity is insensitive to them.** Per-band `κ_crit` there is large (table §4).
     - `++`: `κ_crit ≥ 1.66` (`Δ_g = 20`), `≥ 5.7` (`Δ_g ≥ 50`) for `Δφ = 3`; `≥ 3.8` (`Δ_g = 20`), `≥ 7.8` (`Δ_g ≥ 50`) for `Δφ = 2`.
     - `+−`: `κ_crit ≥ 2.05` and `≥ 4.5` respectively.
     - So a relative error of order 100 % on those nodes still leaves the action positive. On fully controlled states the margin is ≥ `0.175N`
       (`Δφ = 2`, `Δ_g = 20`), `≥ 0.187N` (`Δφ = 3`, `Δ_g = 20`), `≥ 0.183N` (`Δ_g ≥ 50`, both skeletons). With the spin envelope `C_s = 40` it is ≥ `0.138N` (`Δφ = 2`, `Δ_g = 20`).
   * **For `Λ ≤ m < 1.5Λ` only `κ_crit` speaks.** For M02, `++/+−`:

     | `Δ_g` | `Δφ = 3` | `Δφ = 2` |
     |---|---|---|
     | 20 | 0.36 / 0.58 | 0.70 / 0.62 |
     | 50 | 0.79 / 0.66 | 0.78 / 0.53 |
     | 100 | 0.83 / 0.44 | 0.85 / 0.41 |
     | 300 | 0.84 / 0.32 | 0.87 / 0.32 |
     | 1000 | 7.3 / 0.29 | 8.4 / 0.29 |

     Nothing computed here bounds the relative correction on those nodes (`ν ≳ 0.6m` and/or `𝒥x² ≳ 1`), so **positivity below 1.5Λ is not established**.
   * **F_g3 for comparison.** For `+−` below 1.5Λ, `κ_crit` is 0.27 / 0.29 / 0.13 / 0.065 / 0.031 (`Δφ = 3`, `Δ_g` = 20 … 1000). For `Δφ = 2` it is
     negative at `Δ_g = 20`, i.e. already the controlled remainder exceeds `F_hyb` there.
5. **Pole model (referee issue 4) [M, illustration only].** The model `|ρ| ≤ A max(1, (1−x²)^{−4})/m²` with `A = 20` is **neither a fit nor an envelope**:
   controlled values exceed it by up to 3.4× (11× with errors).
   * At `J = 0`, `m = Λ` its integral against the `(1−p)³` weight is **logarithmically divergent**. The finite numbers come from the 48-node quadrature, so
     the earlier claim "positive at the gap state for `Δ_g ≥ 100`" was a quadrature artefact and is **withdrawn**.
   * What survives, with these caveats: with `A = 20` (and `A = 60`, the referee's run) the M02 model margin is positive on every scanned state with
     `m ≥ 1.01Λ` for `Δ_g ≥ 100`. At `A = 60` the minimum is `+0.158N` (`Δ_g = 100`). The model fails below `≈ 1.07Λ` (`Δφ = 3`) and `≈ 1.08Λ` (`Δφ = 2`) at `Δ_g = 20`, and at `m = Λ` at `Δ_g = 50`.
6. **Effective C₁ (referee issue 5).**
   * *Regge-anchored* (relative to leading Regge, as in `ROBUST_FUNCTIONAL.md` §5.1). The exact skeleton has `(exact/Regge − 1)/x² = 2.1–3.9` at
     `x = 0.17–0.5` (`k = 4, 6`, exact blocks). As `x → 0` at `J = 0` the limit is `(2k+3)/4` (the Regge term contributes `−x²/4`; corrected). It diverges as
     `x → 1`. **The assumed `C₁ = 1` is violated**, so the Regge-envelope thresholds are not statements about the true action.
   * *Flat-anchored*, scan-global `max (E+EF)/(uN)` over controlled states:

     | `Δ_g` | `Δφ = 3` | `Δφ = 2` |
     |---|---|---|
     | 20 | 17.6 | 32–42 |
     | 50 | 3.2 | 12–17 |
     | 100 | 0.9–1.3 | 3.0–4.4 |
     | 300 | 0.66–1.0 | 0.8–1.2 |
     | 1000 | 0.57–0.66 | 0.65–0.70 |

     The maxima sit at `m ≈ 10Λ`, large `J` (`u = 0.01`). There `E` is mostly the fallback's truncation-error term, not AdS physics.
   * **The earlier statement "below threshold at `n = 1` for `Δ_g ≳ 50`" is withdrawn.** As a *uniform* envelope these `C₁` give `n = (C₁/0.18)^{1/2} ≈ 2–10`.
     The envelope criterion is simply too crude near the gap. The κ-statement of item 4 is the relevant one.
7. **Rigour.** Certified [certified, `cert_regionA.py`; re-done in Arb by the referee]: **region (A) for the M02 and M01 `B₂` weights**,
   `∫ψ₂P_{−½+iMp}(cosh β)dp > 0` for all `M > 0` with `Mβ ≥ 23.11` (M02), `≥ 23.49` (M01). This is a routine application of `UNIFORM_POSITIVITY.md`
   Corollary A, given (R2) there. Nothing in the remainder is certified. A rigorous remainder needs a convergent representation with error control
   (CHMRS's Mellin formula in ball arithmetic; Mack polynomials for `J > 0`) — not done.

**Bottom line (what threshold is supported).** Assumptions:
* (a) spinning bulk-point formula — the `⟨TTTT⟩` moments at `ν ~ m` are the flat graviton densities;
* (b) transfer of the scalar skeleton's relative remainder to graviton densities, for either skeleton choice;
* (c) spin-specific `1/m²` terms within `C_s cosh β/m²`, `C_s ≤ 40`;
* (d) the `⟨TTTT⟩` heavy formula with `ε_odd = ½`;
* (e) odd `k` by continuation in `k`;
* (f) the heuristic control cut.

Under these, **M02 is positive on all scanned heavy states with `m ≥ 1.5Λ`, for `Δ_g = 20 … 1000` and both skeletons, robustly (κ_crit ≥ 1.66)** [NE/M/I].
Applying this at the gap `Δ_g/1.5` gives a conditional threshold `n ≈ 1.5`: `|ĝ₃|²Λ⁸ ≲ 3.3·10³` and `√(t₄²+α_T²) ≲ 1.5·10⁻⁴` at `Δ_g = 100`. This uses the
leading-log constants of `ROBUST_FUNCTIONAL.md` §5.4; it compares with `8.1·10⁻⁴` at `n = 2.36` and `3·10⁻⁵` at `n = 1`.
**`n = 1` is not supported.** On `[Λ, 1.5Λ)` positivity requires the uncomputed relative remainder at `ν ≳ 0.6m` (and at large `𝒥x²`) to stay below
0.29–0.85. The `+−, J = 4–6` states and the `++, J = 0` state at `m = Λ` (where `ν` reaches `m`, the saddle collision) are the binding ones.
`F_g3` would need 0.03–0.3 there.

## 1. The objects

`Ω = ∫₀¹dx[ψ₂Λ⁴C⁺_{2,Λx} + ½ψ₃Λ⁶C⁻_{3,Λx}] − aΛ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (a/2)Λ¹⁰C⁻_{5,0}` (`ADS4_FUNCTIONAL.md` §0 item 1; M02/M01 weights
from `functionals_certified.json`). The heavy action per state in normal form (`m_flat = m/Λ`, `u = Λ²/m²`) is `B2 + B3 + FW`.

**Skeleton [D at leading order].**
* For one spinning pair the exact Regge moments equal the scalar skeleton `(l+1, Δφ)` times a spin factor that is uniform to `10⁻³` for `ν/m ≤ 0.2`
  (`SPINNING_E2E.md` §7, currents; for `l = 2` up to `≈ −17 cosh β/m²`).
* In `⟨TTTT⟩` both pairs carry `a = k + 5`, so at leading order the skeleton is identical scalars with `Δφ = 3`.
* At `1/m²` the choice is not fixed. The MFT analysis favours `(2, Δφ)`, i.e. `Δφ = 2` here [I]. Both are run, selected by the environment variable `UNIFORM_DF`.
* Odd `k` (`C⁻₃`, `C⁻₅`) uses CHMRS's kernel continued in `k` [I].

**Exact scalar moments [Q].** CHMRS give no closed form valid for all `ν/m`. What is available:
* the Regge series of `m^{2k}C_{k,ν}` and `a m^{2k}Π̂_{k,ν}/γ²` to `1/m⁸` (ancillary file; `α = 2Δφ + k − 1`, identified in `v1`);
* the fixed-`J` series;
* the bulk-point formula (3.40);
* the saddle `S_ν` (app. E), which reproduces BP and Regge only when `J` and `ν` are both large and is not used;
* the Mellin representation (§3.5).

The printed CHMRS (3.29) and the ancillary file differ by `−Δφ(7Δφ² − 9Δφd + 3d² − 1)ηP/3` in the `1/m²` coefficient of `C`. The ancillary version (validated
for `Π̂`) is used.

## 2. Validation

**Bulk point = all leading-ν terms [DV, `t1`].** `Δφ = 3`, `J = 0`, `k = 2`: `R1 = (12.75, −178, −660.25, −1598)`.

**Series vs exact blocks [DV, `v1`]** (`Δφ = 3`, `β ≈ 0.4`, `m = 67.05`; `R_j` = exact/series truncated at `1/m^{2j}`):

| `k` | `ν` | `R_0` | `R_1` | `R_2` | `R_4` |
|---|---|---|---|---|---|
| 4 | 0 | 0.99764 | 1.000073 | 0.9999988 | 0.999999999 |
| 4 | 3 | 1.00542 | 1.000010 | 1.0000007 | 1.000000001 |
| 6 | 0 | 0.98809 | 1.000258 | 0.9999936 | 0.999999991 |
| 6 | 3 | 1.00045 | 0.999983 | 1.0000020 | 1.000000005 |

**Bulk-point regime, exact blocks [DV, `t2`]** (`Δφ = 3`, `J = 0`):

| `m` | `k` | `x` | ex/Regge | ex/BP | `m²(ex/BP − 1)` | ex/series₈ |
|---|---|---|---|---|---|---|
| 16 | 4 | 0.25 / 0.5 | 1.142 / 1.903 | 0.9573 / 0.8835 | −10.9 / −29.8 | 1.00021 / 1.00009 |
| 16 | 6 | 0.25 / 0.5 | 1.040 / 1.979 | 0.8175 / 0.6961 | −46.7 / −77.8 | 1.0015 / 0.9988 |
| 32 | 4 | 0.125 / 0.25 | 1.033 / 1.175 | 0.9897 / 0.9854 | −10.6 / −15.0 | 1.0000003 / 1.0000002 |
| 32 | 6 | 0.125 / 0.25 | 1.009 / 1.198 | 0.9511 / 0.9421 | −50.0 / −59.2 | 1.000002 / 0.999998 |

## 3. The remainder used in the scan (`remainder.py`)

**Nodes.** Per state and node (`ν = Λp`, `x = p√u`), the residual is `ρ_k = (exact − hyb)/E_k`. "Exact" is evaluated in this order:
1. the fixed-J series `r̃₁/m² + r̃₂/m⁴`. Its errors are a geometric tail for each truncated `x²`-series, plus a `1/m⁶` budget `|r̃₂/m⁴|·|r̃₂/r̃₁|/m²`
   (new, an estimate).
2. If that is not controlled and `x ≤ 0.6`: the order-8 fixed-η series. Its error is the last term plus **twice the exact tail of `E_k` beyond `x⁸`**,
   which replaces `3x¹⁰E_k` (that term was 1.7–4× too small for `k ≥ 3`).
3. Otherwise the node is uncontrolled.

**Forward terms (at `ν = 0`).**
* `R₄ = −m²∂_{ν²}(m⁸C₄) − (2m²W − 11/2)`
* `R₅ = m¹⁰C_{5,0} − 2P_{−½}(cosh β)`

**Transfer [I]** (relative to the envelopes): `δB2 = ∫ψ₂(2−x²)ρ₂`, `δB3 = u∫ψ₃ρ₃`, `δFW = a u³(R₄ + σR₅/2)`, bounded by absolute values.

## 4. Scan (`t3_scan.py`, `summarize.py`; `cache/out_t3_<Δ_g>[_Df2].txt`, `cache/out_summary.txt`) [NE + M + I]

M02, `C_s = 0`, per band of `m/Λ`. Each cell gives the maximum `E/N` and `κ_crit` (`++` / `+−`) for `Δφ = 3 | Δφ = 2`:

| `Δ_g` | `[1, 1.5)`: `E/N` | `[1, 1.5)`: `κ_crit` | `[1.5, 3)`: `E/N` | `[1.5, 3)`: `κ_crit` |
|---|---|---|---|---|
| 20 | 0.20/0.08 \| 0.38/0.31 | 0.36/0.58 \| 0.70/0.62 | 0.062/0.059 \| 0.096/0.15 | 1.66/2.05 \| 3.8/4.5 |
| 50 | 0.093/0.059 \| 0.11/0.12 | 0.79/0.66 \| 0.78/0.53 | 0.014/0.015 \| 0.020/0.027 | 7.9/18 \| 7.8/17.7 |
| 100 | 0.022/0.028 \| 0.026/0.047 | 0.83/0.44 \| 0.85/0.41 | 0.006/0.007 \| 0.013/0.014 | 5.8/13.9 \| 14.3/15.9 |
| 300 | 0.006/0.020 \| 0.007/0.021 | 0.84/0.32 \| 0.87/0.32 | 0.004/0.005 \| 0.004/0.005 | 5.8/16.3 \| 14.4/13.6 |
| 1000 | 0.005/0.016 \| 0.003/0.014 | 7.3/0.29 \| 8.4/0.29 | 0.003/0.004 \| 0.002/0.003 | 22/16 \| 19/16 |

Minimum `F_hyb/N` (no remainder) for M02 ranges from `0.2019 / 0.4887` (`Δ_g = 20`) to `0.1996 / 0.2131` (`Δ_g = 1000`). The binding `+−` states are `J = 4–6` at
`m ≈ 1.01–1.2Λ`, the old `F_g3` point; their flat margin is `0.2N` for M02.

## 5. Rigour (`cert_regionA.py`)

**Proposition A** [certified; given (R2) of `UNIFORM_POSITIVITY.md` and correct interval `sin`/`cos`; independently re-done in Arb by the referee].
Let `ψ₂` be the M02 (resp. M01) weight. Then for all `M > 0` and `β > 0` with `Mβ ≥ 23.11` (resp. `23.49`): `∫₀¹ψ₂P_{−½+iMp}(cosh β)dp > 0`.

*Proof.*
1. `ψ₂/p` is positive and decreasing (exact Sturm counts). Hence Lemmas 1–2 of `UNIFORM_POSITIVITY.md` apply.
2. `S(τ) = ∫ψ₂ sin pτ > 0` on `[23.11, 51]`: 32 504 interval cells with `|S′| ≤ ∫pψ₂`.
3. `S(τ) > 0` for `τ ≥ 51`: analytic tail bound.
4. Corollary A applies. ∎

This covers only the large-impact-parameter part of the `B₂` term.

## 6. What is conditional

* (a)–(f) of the bottom line in §0.
* [M] The hybrid is the anchor. The remainder is measured relative to the skeleton's own hybrid, so the hybrid's cross terms are included only on controlled nodes.
* [M] The pole model is an illustration only (item 5 of §0).
* The light side and all other inputs of `ADS4_LIGHT.md` §4 are unchanged.

## 7. Code (`ads4/code/uniform/`)

| file | content |
|---|---|
| `chmrs_series.py` | parses the CHMRS ancillary series (cached) |
| `skel.py`, `skelfast.py` | normalisations, BP, `S_ν`, numeric/vectorised series (any `Δφ`) |
| `v1_series_vs_exact.py`, `t2_bulkpoint_exact.py` | series and bulk point vs exact blocks |
| `t1_structure.py`, `residual_coeffs.py [Δφ]` | bulk-point check, remainder coefficients, curvature subtraction |
| `t4_radius.py` | Domb–Sykes ratio-trend estimate (calibrated on `E_k`) |
| `mellin.py`, `v2_mellin.py` | `J = 0` Mellin check (valid to `m ≈ 30`; defect at `m = 58`) |
| `remainder.py` (`UNIFORM_DF` = 3 or 2), `t3_scan.py`, `summarize.py` | remainder per state; scans of M02, M01, F_g3 |
| `cert_regionA.py` | Proposition A certificate |
| `run_all.sh` | everything (`bash run_all.sh cert`: certificate only) |

## 8. Revisions after `ads4/REFEREE_UNIFORM_REMAINDER.md`

1. (major) Skeleton choice. The `Δφ = 2` (MFT) skeleton was added throughout. The sign flip of the `1/m²` remainder is stated, and all tables are given for both choices.
2. (major) "Series diverges" was withdrawn. The control criterion and the `x ≤ 0.6` fallback are labelled a heuristic cut. A calibrated Domb–Sykes estimate is added (`t4_radius.py`).
3. (major) The statement about uncontrolled nodes is corrected: all states below 1.5Λ and about 75 % in `[1.5, 3)Λ` have them. Positivity is now stated per band through `κ_crit`.
4. (major) Pole model. The claim "positive at the gap for `Δ_g ≥ 100`" was a quadrature artefact (log divergence at `J = 0`, `m = Λ`) and is withdrawn. The model is labelled illustration-only, not a fit and not an envelope.
5. (major) The scan-global `C₁` values are reported. "Below threshold at `n = 1` for `Δ_g ≳ 50`" is withdrawn.
6. (minors)
   * The `x → 0` limit is `(2k+3)/4`.
   * The fallback truncation error is now twice the exact `E_k` tail.
   * A `1/m⁶` budget was added to the fixed-J remainder.
   * The Mellin defect at `m = 58` is stated, not fixed.
   * All scans were re-run with these error terms; the `Δφ = 3` numbers changed by `≤ 0.01` in `κ_crit` and margins.
