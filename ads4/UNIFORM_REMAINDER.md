# Uniform remainder near the gap: the AdS heavy action of the robust functionals from exact scalar-skeleton moments

Status: written 16 Nov 2026 by a derivation/numerics agent (STATUS link 9). **Unrefereed. Not a theorem about AdS4.**
Code: `ads4/code/uniform/` (`bash run_all.sh`; outputs in `uniform/cache/`, git-ignored; it needs the CHMRS ancillary file
`anc/regge_moments_order8.m` of arXiv:2106.10274 copied to `cache/`). It imports, unchanged, `e2e/moments.py` (exact blocks),
`robust/rcommon.py`, `robust/check.py`, `robust/ads.py` (functionals, flat parts, hybrid), and reads
`robust/functionals_certified.json`. Nothing in other directories is modified.

Marks: **[certified]** exact rational / interval arithmetic; **[D]** derived; **[DV]** derived and checked against exact blocks;
**[Q]** quoted from CHMRS (arXiv:2106.10274, LaTeX and ancillary file read); **[NE]** floating-point numerics read off from numbers;
**[M]** statement about a stated model/extrapolation; **[I]** inferred/assumed. CHMRS = arXiv:2106.10274.

## 0. Results (answers to the four questions)

1. **Structure [D + DV].** The scalar skeleton of every `⟨TTTT⟩` moment is the identical-scalar moment with `Δφ = 3`
   (`a = 2Δφ + k − 1 = k + 5`, the index of both spinning pairs, `SPINNING_E2E.md` §8). Its exact heavy action organises as
   **flat (bulk point) + leading-Regge curvature + remainder**:
   `m^{2k}C_{k,ν}[Δ,J] = E_k(x)P_J(1−2x²) + 2[P_{−½+iν}(cosh β) − J₀((2J+1)x)] + ρ_k·E_k(x)`, `x = ν/m`, `E_k = (2−x²)/(1−x²)^{1+k/2}`.
   * The **leading-ν part of the exact (CHMRS order-1/m⁸) series is exactly the bulk-point formula**, at every order `j ≤ 4`,
     symbolically in `k, Δφ, J` [DV, `t1`]. So the `O(ν²/m²)` terms are *not* an unknown correction: they are the flat density.
   * The **remainder after subtracting flat + leading-Regge curvature** has exact coefficients (CHMRS ancillary, curvature part
     subtracted in closed form, `residual_coeffs.py`): `ρ_k = [r̃₁(x,J)/m² + r̃₂(x,J)/m⁴ + …]/E_k`; e.g. `k = 2`:
     `r̃₁ = 12.875 + (−178.08 + 1.375𝒥)x² + (−660.2 + 309.1𝒥 − 5.68𝒥²)x⁴ + …` (`𝒥 = J(J+1)`). The subtraction lowers the `𝒥`-degree of each
     coefficient by one; the left-over `(𝒥x²)^{j−1}/m²` terms are the hybrid's uncontrolled `O(uβ²)` cross terms, now with exact coefficients.
   * The **CHMRS series is validated against exact blocks** (identical scalars `Δφ = 3, 1.75`, `k = 4, 6`, `ν = 0,1,3`, `β ≈ 0.4, 1`):
     exact/series = `1 ± 10⁻⁸` at `m ≈ 67` (order 8), with `m^{2(j+1)}(R_j − 1)` bounded at each order; the mutation
     `a → a − ½` fails [DV, `v1`]. Bulk-point regime (`J = 0, 2, 4`, `x ≤ 0.5`, `m = 12–32`): exact/BP − 1 = `O(1/m²)`
     (`m²(ex/BP − 1) = −4 … −80`, `k = 4, 6`), while exact/Regge-leading deviates by up to **+98 %** at `x = 0.5` [DV, `t2`]. An independent
     Mellin evaluation for `J = 0` (`mellin.py`, CHMRS §3.5 with the `J = 0` Mack coefficients) agrees with the series to `4·10⁻⁸–2.5·10⁻⁵`
     at `m = 28` [DV]; at `m = 58` it differs by `5·10⁻⁴` (unresolved numerical issue of my Mellin implementation, not used further).
   * **Spin factor.** Not computed here; taken as `1 + O(1/m²)` uniformly [I], with an envelope `C_s cosh β/m²` (`C_s = 40`) for the
     measured spin-specific `l = 2` piece (`SPINNING_E2E.md` §6.3).
2. **AdS action of M02 (M01, F_g3) on all heavy states near the gap** [NE + M + I, `t3`; ≈ 2800–4000 physical states per `Δ_g`,
   `m ≤ 10Λ`, `b̂ ≤ 40`, both channels]: hybrid action `F_hyb` (flat + leading-Regge curvature, `robust/ads.py`) minus the transferred skeleton
   remainder. Where the remainder series is controlled, **M02 stays positive with margin ≥ 0.1875·N (Δ_g = 20), ≥ 0.198·N (Δ_g ≥ 50)**;
   with the spin envelope `C_s = 40`: ≥ 0.158·N (Δ_g = 20). The **uncontrolled** nodes are `ν ≳ 0.6m` at the lowest states
   (`m ≲ 1.3Λ`, `J ≤ 6`), where the series (in `x²` and `𝒥x²`) does not converge: positivity there holds iff the relative AdS correction on those
   nodes is below `κ_crit` = **0.37 / 0.59 (Δ_g = 20, `++`/`+−`), 0.79 / 0.68 (50), 0.83 / 0.45 (100), 0.84 / 0.33 (300), 7.3 / 0.30 (1000)**.
   With a pole-model extrapolation of the remainder, `|ρ| ≤ 20 max(1,(1−x²)^{−4})/m²` [M, fitted to the controlled values], M02 is
   **positive on every scanned state for Δ_g ≥ 100, including the gap state**; for Δ_g = 50 it fails only at the single state `++, J = 0, m = Λ`;
   for Δ_g = 20 it fails for `++, J = 0, m ≤ 1.06Λ` (positive from `1.07Λ`). These failures sit exactly at the saddle collision `ν → m` (only `J = 0` at
   `τ = Δ_g` reaches `ν = m` with the smearing `Λ = min m`), where the pole model makes the integrated remainder log-divergent against the
   `(1−p)³` weight, i.e. the model, not the physics, decides them. **So: n = 1 is supported (not proved) for Δ_g ≳ 100, and n ≈ 1.01 (50),
   ≈ 1.07 (20), conditional on the ingredients of §6; the binding unknown is the remainder at `ν/m ≳ 0.6` on the lowest few states.**
   F_g3 for comparison: positive in the same scheme but fragile — `+−, J = 4, m ≈ 1.01–1.04Λ` has `κ_crit = 0.037–0.15` (Δ_g = 1000–100),
   and its `++` island margin is `0.005N`.
3. **What C₁ effectively is.** Two different envelopes were conflated in the earlier chain:
   * *Regge-anchored* (`|r| ≤ C₁Λ²/m²` relative to the leading Regge term, `ADS4_FUNCTIONAL.md` statement R, `ROBUST_FUNCTIONAL.md` §5.1):
     the exact skeleton has `(exact/Regge − 1)/x² ≈ 2.1–3.9` at `x = 0.17–0.5` (`k = 4, 6`, exact blocks, `t2`) and `→ (k+1)/2` analytically as `x → 0` at `J = 0`
     (from `Q_{k,2}`, CHMRS (3.35)), diverging like `(1−x²)^{−1−k/2}` as `x → 1`. **The assumed `C₁ = 1` is violated by a factor 2–4 at `ν ≤ m/2`
     and without bound near `ν = m`**; the Regge-envelope thresholds (2.16Λ etc.) are therefore not reliable statements about the true action.
   * *Flat-anchored* (`|r| ≤ C₁uN` relative to the flat component norm, the certified-margin threshold `n = (C₁/0.18)^{1/2}` of
     `ROBUST_FUNCTIONAL.md` §0 item 4): the transferred remainder gives, near the gap (`m < 1.5Λ`, controlled nodes),
     `E/N ≤ 0.18 (Δ_g = 20), 0.09 (50), 0.022 (100), 0.0033 (300), 0.0026 (1000)`, i.e. **`C₁^eff ≈ 0.2–0.4, 0.1–0.2, 0.02–0.05, < 0.01, < 0.01`**
     (an upper bound: absolute values inside the integrals; at large Δ_g dominated by my conservative truncation-error term). With the
     certified margin 9/50 this is below threshold already at `n = 1` for Δ_g ≳ 50 — on the controlled nodes.
4. **Rigour.** Certified: (a) **region (A) for the M02 and M01 `B₂` weights** — `∫ψ₂(p)P_{−½+iMp}(cosh β)dp > 0` for all `M > 0` with
   `Mβ ≥ 23.11` (M02), `≥ 23.49` (M01) [certified, `cert_regionA.py`: exact Sturm counts for `ψ₂/p > 0` decreasing; 32 504 / 55 696 interval
   cells for the sine transform; analytic tail; Corollary A of `UNIFORM_POSITIVITY.md` with `f → ψ₂`]. The sine transform of the M02 weight
   is negative on three windows (8.75–10.65, 15.09–17.34, 22.49–23.10), so regions (B), (C) of Theorem 1 need new certificates (not done).
   (b) The exact remainder coefficients are exact rationals [Q + D]. Everything else (series truncation, exact-block numerics, scans) is
   floating point. A rigorous remainder would need a convergent representation with error control, e.g. the CHMRS Mellin formula in ball
   arithmetic (exact for `J = 0` here; Mack polynomials for `J > 0`) — identified, not done.

## 1. The objects

`Ω = ∫₀¹dx[ψ₂Λ⁴C⁺_{2,Λx} + ½ψ₃Λ⁶C⁻_{3,Λx}] − aΛ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (a/2)Λ¹⁰C⁻_{5,0}` (`ADS4_FUNCTIONAL.md` §0 item 1; M02/M01 weights
from `functionals_certified.json`). Heavy action per state in normal form (`m_flat = m/Λ`, `u = Λ²/m²`): `B2 + B3 + FW`.

**Skeleton [D].** For one spinning pair the exact Regge moments equal the scalar skeleton with dimensions `(l+1, Δφ)` times a spin factor
(`SPINNING_E2E.md` §6.3, §7: uniform to `10⁻³` for `ν/m ≤ 0.2`, `m ≥ 40`, currents; for `l = 2` up to an extra `≈ −17 cosh β/m²`). In
`⟨TTTT⟩` both pairs carry `a = k + 5` (`SPINNING_E2E.md` §8), i.e. the skeleton is the identical-scalar moment with `Δφ = 3`. Odd `k`
(`C⁻₃`, `C⁻₅`) uses CHMRS's kernel `Ω̃_{k,ν}` continued in `k` [I]; the bulk-point check of §2 holds symbolically in `k`.

**Exact scalar moments [Q].** CHMRS give no closed form for all `ν/m`. Available exactly: (i) the Regge series of
`m^{2k}C_{k,ν}[Δ,J]` and of `a m^{2k}Π̂_{k,ν}/γ²` to `1/m⁸` (ancillary `regge_moments_order8.m`; coefficients polynomial in `ν, η`, with
`P^{(n)}_{−½+iν}(η)`, `n ≤ 4`; the symbol `α` there is `2Δφ + k − 1`, identified in `v1`: the alternative `2Δφ + k − 3/2` fails against exact
blocks); (ii) the same at fixed `J` (`Q_{k,2j}` polynomials); (iii) the bulk-point formula (3.40) (leading in `1/m` at fixed `ν/m`, `J`);
(iv) the large-(Δ,J,ν) saddle `S_ν` (app. E) — checked here: it reproduces BP and Regge only when `J`, `ν` are both large
(`S/BP = 1.11` at `J = 2`, `ν = 100`; `0.99` at `J = 1000`), so it is not used; (v) the Mellin representation (§3.5), exact.
Note: the printed CHMRS (3.29) and the ancillary differ by `−Δφ(7Δφ² − 9Δφd + 3d² − 1)η P/3` in the `1/m²` coefficient of `C`; the ancillary is the one
validated against exact blocks (for `Π̂`) here, and is used.

## 2. Structure and validation

**Bulk point = all leading-ν terms [DV, `t1`].** Splitting each `Q_{k,2j}` (`j ≤ 4`) by degree in `ν²`, the top degree equals the
`x^{2j}` coefficient of `(2−x²)(1−x²)^{−1−k/2}P_J(1−2x²)` identically in `k, Δφ, 𝒥`. The next degree gives `r₁(x) = Σ_j R1_j x^{2j−2}`, the one after
`r₂`. `Δφ = 3`, `J = 0`, `k = 2`: `R1 = (12.75, −178, −660, −1598)`, ratios `→ 2.4` (consistent with a singularity at `x → 1` of order `(1−x²)^{−5…−6}`
in `r₁`, i.e. `(1−x²)^{−3.5±1}` relative to `E_k`; extrapolation [M]). Large-`J` growth `R1_j ∝ 𝒥^j` before, `𝒥^{j−1}` after subtracting the exact Regge
curvature (closed form: `2(−e/2)^n e_r(n)/(n!)²`, `e = 2𝒥 + ½`, `e_r` elementary symmetric polynomials of `(i+½)²`, `i < n`).

**Series vs exact blocks [DV, `v1`]** (`R_j` = exact/series truncated at `1/m^{2j}`; selection, `Δφ = 3`, `β ≈ 0.4`, `m = 67.05`):

| `k` | `ν` | `R_0` | `R_1` | `R_2` | `R_4` |
|---|---|---|---|---|---|
| 4 | 0 | 0.99764 | 1.000073 | 0.9999988 | 0.999999999 |
| 4 | 3 | 1.00542 | 1.000010 | 1.0000007 | 1.000000001 |
| 6 | 0 | 0.98809 | 1.000258 | 0.9999936 | 0.999999991 |
| 6 | 3 | 1.00045 | 0.999983 | 1.0000020 | 1.000000005 |

Same quality for `β ≈ 1` and for `Δφ = 1.75`; at `m ≈ 10` the series (asymptotic) is off by up to 30 % (`Δφ = 3`).

**Bulk-point regime, exact blocks [DV, `t2`]** (`Δφ = 3`, `J = 0`; `ex` = exact normalised `Π̂`):

| `m` | `k` | `x = ν/m` | ex/Regge | ex/BP | `m²(ex/BP − 1)` | ex/series₈ |
|---|---|---|---|---|---|---|
| 16 | 4 | 0.25 / 0.5 | 1.142 / 1.903 | 0.9573 / 0.8835 | −10.9 / −29.8 | 1.00021 / 1.00009 |
| 16 | 6 | 0.25 / 0.5 | 1.040 / 1.979 | 0.8175 / 0.6961 | −46.7 / −77.8 | 1.0015 / 0.9988 |
| 32 | 4 | 0.125 / 0.25 | 1.033 / 1.175 | 0.9897 / 0.9854 | −10.6 / −15.0 | 1.0000003 / 1.0000002 |
| 32 | 6 | 0.125 / 0.25 | 1.009 / 1.198 | 0.9511 / 0.9421 | −50.0 / −59.2 | 1.000002 / 0.999998 |

`J = 2, 4` behave alike away from zeros of `P_J(1−2x²)` (where ratios are meaningless). (Normalisation: the e2e moment divided by
`q_{Δ,J}γ²/a_{Δφ}`, CHMRS (2.43), (2.57).)

## 3. The remainder used in the scan (`remainder.py`)

Per state and node `ν = Λp` (`x = p√u`): `ρ_k = (exact − hyb)/E_k` with hyb = BP + `2[P − J₀((2J+1)x)]`. Exact = fixed-J series
(`r̃₁/m² + r̃₂/m⁴`, curvature subtracted); if not convergent (last term `> 0.8×` previous and not negligible) and `x ≤ 0.6`, the fixed-η series
(order 8; error = last term + `3x¹⁰E_k` for its truncation of the bulk-point `x`-series); otherwise the node is **uncontrolled**. Forward
terms at `ν = 0` (uniform in β): `R₄ = −m²∂_{ν²}(m⁸C₄) − (2m²W − 11/2)`, `R₅ = m¹⁰C_{5,0} − 2P_{−½}(cosh β)` (hyb = flat `2𝒥 − 5`, `2` + leading
curvature); e.g. `R₄ = −35.1` against `2m²W − 5.5 = 63 320` at `m = 300, J = 200`. Transfer to the graviton [I]: relative to the envelope,
`δB2 = ∫ψ₂(2 − x²)ρ₂`, `δB3 = u∫ψ₃ρ₃`, `δFW = a u³(R₄ + σR₅/2)`; bounded by absolute values:
`E = ∫|ψ₂|(2−x²)(|ρ₂| + err) + u∫|ψ₃|(|ρ₃| + err) + a u³(|R₄| + |R₅|/2 + errs)` on controlled nodes, `U` = the same weights on uncontrolled nodes.

## 4. Scan (`t3_scan.py`; `cache/out_t3_<Δ_g>.txt`) [NE + M + I]

`min (F_hyb − E)/N` over states with no uncontrolled node; `κ_crit = min (F_hyb − E)/U` over the others; `C_s = 0` unless stated.

| `Δ_g` | M02 `++` | M02 `+−` | M02 `κ_crit` `++/+−` | M02 pole model (min over states) | M01 `++` | F_g3 `++` | F_g3 `+−` `κ_crit` |
|---|---|---|---|---|---|---|---|
| 20 | 0.1875 (0.158 with `C_s = 40`) | 0.435 | 0.37 / 0.59 | < 0 for `++, J = 0, m ≤ 1.06Λ`; ≥ +0.02 from 1.07Λ | 0.089 | 0.011 | 0.27 |
| 50 | 0.198 | 0.390 | 0.79 / 0.68 | −0.30 at `m = Λ, J = 0` only; ≥ +0.198 elsewhere | 0.098 | 0.003 | 0.32 |
| 100 | 0.199 | 0.388 | 0.83 / 0.45 | ≥ +0.199 (all states) | 0.0995 | 0.005 | 0.15 |
| 300 | 0.1995 | 0.382 | 0.84 / 0.33 | ≥ +0.199 | 0.100 | 0.005 | 0.072 |
| 1000 | 0.1996 | 0.382 | 7.3 / 0.30 | ≥ +0.199 | 0.100 | 0.006 | 0.037 |

`min F_hyb/N` (no remainder): M02 `0.2019 / 0.4887` (20) … `0.1996 / 0.2131` (1000). Near-gap band `m < 1.5Λ`, M02: `E/N ≤ 0.18, 0.088, 0.022,
0.0033, 0.0026`; `U/N ≤ 1.4, 1.15, 1.18, 1.18, 0.13` (`++`), the large values at `J = 0, m ≈ Λ`. The binding `+−` states are `J = 4–6`,
`m ≈ 1.01–1.2Λ` (the old `F_g3` point; for M02 its flat margin is 0.2N). Per-mass pole-model margins near the gap are printed at the end
of each output file.

## 5. Rigour (`cert_regionA.py`, `cache/cert_regionA.log`)

**Proposition A [certified, given (R2) of `UNIFORM_POSITIVITY.md` and the correctness of `mpmath.iv` for `+ − × ÷ sin cos`].** For the M02
(resp. M01) weight `ψ₂` and all `M > 0`, `β > 0` with `Mβ ≥ 23.11` (resp. `23.49`): `∫₀¹ψ₂(p)P_{−½+iMp}(cosh β)dp > 0`.
*Proof.* `ψ₂ = (1−p)³p g(p)`, `g` has no root on `[0,1]` and `φ = ψ₂/p` has `φ′ = (1−p)²N(p)` with `N` root-free, `N(0) < 0` (exact Sturm counts),
so Lemmas 1–2 of `UNIFORM_POSITIVITY.md` hold verbatim (kernel `S(Mt) + R_M(t)`, `R_M > 0`). `S(τ) = ∫ψ₂ sin(pτ) > 0` on `[23.11, 51]` by 32 504
interval cells (`|S′| ≤ ∫pψ₂ = 361807/100800`), and for `τ ≥ 51` from `τ⁹S ≥ (34244/125)τ⁶ − …` (exact root count). Corollary A then applies. ∎
This is the large-impact-parameter part of the `B₂` term only; it says nothing about the near-gap problem.

## 6. What is conditional

* **[I] Spinning bulk-point formula**: that the exact `⟨TTTT⟩` moments at `ν ~ m` are the flat graviton densities (the flat functional in
  `F_hyb`). Gap 2 of `ADS4_FUNCTIONAL.md`; nothing here derives it.
* **[I] Transfer** of the skeleton's relative remainder to the graviton densities ("spin factor uniform in `ν/m`"): numerical evidence only
  for currents (`l = 1`), `ν/m ≤ 0.2`; for `T` up to the `−17 cosh β/m²` piece; nothing at `ν/m ≳ 0.3`.
* **[I] Stress-tensor heavy formula** (`C^± = (2P/m^{2k})(…)`, `ε_odd = ½`, `k = 2, 3, 5`, two spinning pairs): `SPINNING_TTTT.md`, `ADS4_FUNCTIONAL.md` §4.
* **[I] Odd `k`** via analytic continuation of CHMRS's kernel in `k`.
* **[M] The hybrid** `F_hyb` is what the remainder is measured against; the remainder was computed *relative to the skeleton's own hybrid*,
  so the hybrid's cross-term error is included to the order of the series — but only on controlled nodes.
* **[M] Pole-model extrapolation** on uncontrolled nodes (`x ≳ 0.6`, small `J`, `m ≲ 1.3Λ`) — the open point; it decides the single `J = 0`
  states at `m ≈ Λ` for `Δ_g ≤ 50`.
* Light side and all other inputs of `ADS4_LIGHT.md` §4 unchanged.

## 7. Code (`ads4/code/uniform/`)

| file | content |
|---|---|
| `chmrs_series.py` | parses the CHMRS ancillary series into sympy (cached) |
| `skel.py`, `skelfast.py` | normalisations `q, a_Δφ, γ`, BP, `S_ν`, numeric/vectorised series |
| `v1_series_vs_exact.py` | series vs exact-block moments (`e2e/moments.py`) |
| `t1_structure.py`, `residual_coeffs.py` | bulk-point check, remainder coefficients, curvature subtraction |
| `t2_bulkpoint_exact.py` | exact blocks vs Regge / BP / series at `ν/m ≤ 0.5` |
| `mellin.py`, `v2_mellin.py` | exact `J = 0` action from the Mellin representation; check vs series |
| `remainder.py`, `t3_scan.py` | per-state remainder; scan of M02, M01, F_g3 |
| `cert_regionA.py` | Proposition A certificate |
| `run_all.sh` | everything (`bash run_all.sh cert`: certificate only) |
