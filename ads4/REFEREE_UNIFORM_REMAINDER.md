# Referee report: `UNIFORM_REMAINDER.md` and `ads4/code/uniform/`

Referee: independent, hostile. 16 Nov 2026. Nothing under review was edited. Scratch code: `ads4/code/referee_uniform/`
(`bp_check.py`, `cert_arb.py`, `pole_fit.py`, `trunc_m6.py`, `unc_where.py`, `t3_pole_variant.py`; outputs in its `cache/`, git-ignored).

## What I re-did

| check | result |
|---|---|
| Own parse of `seriesCfixedJ` (ancillary), d = 3; leading-ν part vs bulk point `(2−x²)(1−x²)^{−1−k/2}P_J(1−2x²)`, with `P_J` from the Legendre ODE recursion (not the 2F1 product used in `t1`) | **zero difference at orders 1/m⁰…1/m⁸, symbolically in k, Δφ, 𝒥.** `R1(k=2,Δφ=3,J=0) = (51/4, −178, −2641/4, −1598)`, as in the paper |
| Curvature subtraction (η − 1 = (2𝒥+½)/m² exactly; 1/m², 1/m⁴ coefficients via e₁, e₂) | correct; `r̃₁(0) = 12.875` reproduced |
| Prop. A certificate re-done in **Arb** (python-flint, 200 bits), S(τ) from sympy's own integral (not the IBP formula), ball bisection on [τ₀, 80], own tail polynomial for τ ≥ 80 | **PASS** for M02 (τ₀ = 23.11) and M01 (23.49); φ = ψ₂/p > 0, decreasing (exact real roots). Only sin/cos enter (no `iv.exp`). Negative windows of S: 8.76–10.66, 15.10–17.35, 22.50–23.11 (M02), as stated |
| `t3_scan` for Δ_g = 100 re-run (cached skeletons, output to my cache) | **bit-identical** to `out_t3_100.txt` (margin 0.199N; near-gap E/N ≤ 0.0221) |
| fixed-J (r̃₁/m² + r̃₂/m⁴) vs full order-8 fixed-η series near the gap (Δ_g = 20, 50; x = 0.2–0.5) | differences (1/m⁶ and beyond) lie within the code's tail errors at controlled nodes |
| Regge-anchored ratio from the `t2` table | 2.1–3.9 reproduced by hand |
| Mellin J = 0 check (`out_v2.txt`) | 4·10⁻⁸ … 2.5·10⁻⁵ at m = 28.8, as stated |

## Issues

1. **(major, applicability) The skeleton Δφ = 3 is fixed only at leading order. The remainder depends strongly on it.** `a = k+5` for both pairs
   is the leading-order statement of `SPINNING_E2E` §8. It is consistent with the (l+1, Δφ) skeleton, so (3,3) is the right *leading* choice.
   The object of this paper, however, is the 1/m² remainder. `REFEREE_SPINNING_MFT` Issue 3 found that the per-l 1/m² term is clean only
   relative to an empirically chosen skeleton (2,Δφ), not (3,Δφ). The remainder is not robust to this choice: `R1₁(k=2, J=0)` = +12.75 (Δφ = 3), −19.25 (2), −22.25 (1.75),
   so the sign flips. The C_s = 40 envelope was measured at ν/m ≤ 0.2 only. Fix: show the Δφ = 2 scan, or carry a skeleton-choice error.
2. **(major, soundness/labels) The "controlled region ν ≲ 0.6m" is a code cut, not a derived boundary.** `t3_scan` hard-codes `x ≤ 0.6` for the
   fallback. Control otherwise means a single ratio test (≤ 0.8) on 3–4 terms of the x²-series. Four terms cannot fix a radius of convergence: the ratios for
   k=2, J=0 are −14, 3.7, 2.4 (k=2, J=4: 131, −2.2, 9.0). "Does not converge / diverges" (§0 item 2) is unproven. Likely the radius is 1 and the
   real problem is non-uniformity of 1/m at x → 1. The "error = geometric tail" is an estimate, not a bound.
3. **(major, labels) Uncontrolled nodes are not confined to "m ≲ 1.3Λ, J ≤ 6, ν ≳ 0.6m".** Counts from the cached skeletons:
   - in band [1.5, 3)Λ, 738/986 (Δ_g = 20) and 976/1192 (Δ_g = 100) states have uncontrolled nodes at **x < 0.6**, with J up to 58;
   - at m ≥ 3Λ, about 100 states have them, with J up to 125;
   - in [1, 1.5)Λ, **every** state has some (Δ_g ≤ 300: `min(F−E)/N = +inf`).

   So "M02 ≥ 0.1875N where controlled" is a statement about m ≥ 1.5Λ. Near the gap, only κ_crit (0.30–0.84) speaks. The large-J
   uncontrolled nodes come from the unresummed 𝒥x² cross terms. The J-independent pole model is meaningless there, though U/N is small.
4. **(major, soundness) The pole model is neither fitted nor an envelope.**
   - On *controlled* nodes, |ρ|m²/[20 max(1,(1−x²)^{−4})] reaches 2.1 (J ≤ 2, Δ_g = 50) and 3.4 (J ≤ 6, Δ_g = 100), and up to 3.7 / 11 with err.
     The model under-bounds the data it claims to be fitted to.
   - With A = 60, Δ_g = 100 stays positive (min +0.158N) and Δ_g = 50 fails for m ≤ 1.02Λ. The Δ_g ≥ 100 conclusion is therefore robust in A.
   - **But** at J = 0, m = Λ, x = p, so the model integral `∫ψ₂(2−p²)(1−p²)^{−4}` is log-divergent for *every* Δ_g. The coefficient is 20g(1)/(16m²) with g(1) = 763.7, i.e.
     ≈ 0.098 per e-fold at Δ_g = 100. The integral is finite here only because of 48 Gauss nodes (407 at 48 nodes, 544 at 200 nodes uncapped).
   - "Positive including the gap state for Δ_g ≥ 100" is therefore a quadrature artifact, while the same divergence is blamed for the failures at Δ_g ≤ 50. This is
     inconsistent. The code comment says A = 10, but the code uses 20.
5. **(major, labels) The effective-C₁ statement is cherry-picked.**
   - The scan's own global value is `C1_eff = max (E+EF)/(uN)` = **17.3 (Δ_g = 20), 3.2 (50), 0.88–0.99 (100), 0.66–0.78 (300), 0.57–0.66 (1000)**, at m = 10Λ and large J.
   - §0 item 3 quotes only the m < 1.5Λ band (0.2–0.4, …, < 0.01). It then concludes "below threshold already at n = 1 for Δ_g ≳ 50", which is false for the
     uniform-envelope threshold n = (C₁/0.18)^{1/2} (n ≈ 4 at Δ_g = 50).
   - At large m the E term is largely the fallback's bulk-point truncation, not AdS physics. That is a further reason the number is not "what C₁ effectively is".
6. **(minor)**
   - The x → 0 limit of (exact/Regge − 1)/x² at J = 0 is **(2k+3)/4, not (k+1)/2**: the Regge term contributes −x²/4. Check: m = 32, k = 4, x = 0.125 gives
     1.033 with (2k+3)/4, versus 1.029 with (k+1)/2.
   - The fallback error 3x¹⁰E_k under-estimates the bulk-point truncation for k ≥ 3 by 1.7–4× at x = 0.45–0.6.
   - The 1/m⁶ truncation of the fixed-J remainder has no error term. It is covered numerically at the tested points.
   - The Mellin code disagrees at m = 58 (5·10⁻⁴), so the J = 0 "independent check" is weakened.

## Verdicts

* **Correctness:** claim 1 (bulk point = leading-ν part; exact coefficients; series vs blocks) is correct. Claim 4 is correct (Arb-confirmed). Issue 6 has small errors.
* **Soundness:** claim 2 is weak. Near the gap nothing is controlled; the pole model is not an envelope, and it diverges at the gap state (Issues 2–4).
  Claim 3 is misreported (Issue 5).
* **Labels:** [DV] for claim 1 and [certified] for claim 4 are fair, the latter given (R2) of `UNIFORM_POSITIVITY.md`. The phrases "series diverges", "uncontrolled = ν ≳ 0.6m at m ≲ 1.3Λ", "fitted pole model"
  and "C₁^eff < 0.01" are overstated.
* **Applicability:** this hinges on two inputs. The first is the skeleton choice at 1/m² (Issue 1). The second is the transfer of a scalar relative remainder to graviton densities at
  x up to 1. The latter is [I], with evidence only at ν/m ≤ 0.2.
* **Novelty:** low for claim 1. CHMRS's (3.40) is the resummed leading-ν part of their own fixed-J series, so this is a reproduction. Claim 4 is a routine application of
  `UNIFORM_POSITIVITY` Corollary A. The value is in the scan and the honest κ_crit tables.
