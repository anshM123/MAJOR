# Independent reproduction: the H_K EE-kernel reduction and Lemma Y

*23 Nov 2026, independent checker. Asked for by REFEREE_REPORT_2 §M3 and §6 item 3. No existing file was edited and nothing was committed.*
* *Code:* `ads4/code/indep/` (new). Outputs are in `ads4/code/indep/cache/`; its `.gitignore` ignores everything except itself.
* *Independence.* Everything was rebuilt from the statements in EESTEP_CLOSE, LRC_CLOSE, HK_CLOSE, HK_RESIDUAL, UNIFIED_THEOREM_v3,
  R4ARC_HEFT and H7 §18.1. That covers the frame, the integrand, the Casimir operator, the blocks, the D-functions, the structure
  enumeration and the explicit contacts.
  * No file from `ads4/code/` outside `indep/` was imported, read or run.
  * The original code was **not** read, even after the fact: no disagreement needed diagnosing.
  * The only inputs taken from the project are the published data in `ads4/prl/data/` (`kernels.json`,
    `q_deshift_continuum_c4{v,d}.json`, `q_deshift_xp_C4_{value,deriv}.txt`).
* *Compute:* about 1.2 CPU-hours in total, at most 2 processes, one thread each, under 3 GB.

Grades are as in the project: [P] proved; [D] derived; [D-sketch]; [DV] derived and checked numerically; [NE]; [M]; [A]; [I].
"Independently reproduced" means re-derived and re-computed by a second party with separate code.

---

## 0. Summary

| # | claim (source) | independent result | agrees? | grade after check |
|---|---|---|---|---|
| A-i.1 | Exact reduction: on the x₀-circle the EE dressing τ₁₃D is a degree-≤ 8 trigonometric polynomial over (X₂₀X₄₀)⁴; its sin-odd part integrates to 0; P_even = Σ_{b=−4}^{4} c_b t^b with ν-independent c_b, t = X₂₀/X₄₀ (EESTEP §1.2) | Chebyshev tail beyond degree 8: **5·10⁻¹⁵**. Laurent residual: 8·10⁻¹³. Odd part: up to 100 % of the even part pointwise. **K^EE = Σ_b c_bΨ^{(b)} to 2.6·10⁻¹²** (9 points × 3 ν × 3 polarisation sets), through cancellations up to 6·10⁴ | **yes** | [P + DV], independently reproduced |
| A-i.2 | Each Ψ^{(b)} = ∫S_ν t^b is the scalar continuous-spin partial wave (1, l = −½+iν) with (a, b_DO) = (0, −b), and Ψ^{(b)} = N(ν) r² v^{−b/2} Ĝ_{−b} = Ψ^{(−b)} (EESTEP §1.3–1.4) | Own d = 3 continuous-spin recursion (own derivation of the radial Casimir, checked against Dolan–Osborn to 10⁻³¹). Ratio spread over 9 points: **2.8·10⁻¹⁵ … 1.1·10⁻⁹** for every b. Fitted prefactor **r^{2.00000000} v^{−b/2}** (residual ≤ 4·10⁻¹⁰). Ψ^{(b)} = Ψ^{(−b)} to ≤ 8·10⁻¹⁴. Controls: the wrong sign of b gives 1.7–7.0; the scalar block alone gives 0.7–22 | **yes** | [D + DV], independently reproduced |
| A-i.3 | End to end: K^EE_ν(x) = N(ν) r² Σ_b c_b(x) v^{−b/2} Ĝ_{−b}(x; ν), with N(ν) fixed by the scalar kernel at one point (EESTEP §1.4: 1.6·10⁻⁷ – 1.5·10⁻⁶) | **1.1·10⁻⁹ / 6.8·10⁻⁹ / 4.7·10⁻⁸** for 3 polarisation sets, with b-sum cancellations 9·10² – 6·10⁴. Control (scalar ν-profile only): 1.2 – 4.5·10³ | **yes** (higher precision than the original) | [P + DV], independently reproduced |
| A-i.4 | Euler identity Ĝ_{−b} = v^b Ĝ_b; the direct negative-b series "loses up to 9 digits" (EESTEP §1.3) | Identity holds to **10⁻⁴⁵ / 4·10⁻³⁸** (50-digit arithmetic, ρ = 0.55). At ρ = 0.85 the direct negative-b series at N = 200 is wrong by O(1) **even in 50-digit arithmetic**: the problem is truncation (coefficients grow like N^{2\|b\|−1}), not only round-off | identity yes; diagnosis differs (immaterial) | [P + DV] |
| A-ii.1 | LRC\*: edge resummation Ĝ_b → 2Re[c_l e^{lβ}F_b(ρ)], F_b = (z/4ρ)^h ₂F₁(h, h+b; 2h; z), corrections O(e^{−2β}) (LRC §2, EESTEP §1.3) | Deviation relative to the edge envelope: 1–13 % at β = 2, then **2·10⁻⁹ – 4·10⁻⁸ at β = 10**, i.e. e^{−2β}. With F_{−b} in place of F_b: O(1) | **yes** | [D + DV], independently reproduced |
| A-ii.2 | De-shift survival exponents of (exact − LR), β ≈ 8.5–13.5, ρ = 0.2/0.5/0.85, each b-block b = 0…4: w2g −0.50; w2s ≈ −0.4 (βη^{−½}); w4g −0.50 up to the double-precision floor; w4p **+1.00**. Residue formula at 3i/2 matches to 10⁻⁶ (EESTEP §2, LRC §3.2) | w2g **−0.50** (all b, all ρ); w2s −0.36 … −0.41; w4g −0.50 for β ≲ 10, then the floor (as in the original); **in 45-digit arithmetic −0.500 out to β = 13** (b = 0, 2); w4p **+1.00** (all b, all ρ). Residue ratio **0.999992 – 1.000007**. Beyond the floor (45-digit arithmetic, ρ = 0.5): see §2.3 | **yes** | [D + DV], independently reproduced |
| A-ii.3 | Corner count c < a + 1 − J₀. Criteria: J₀ < 3/2 (k = 2, 3, 5), J₀ < 5/2 (k = 4, f_pos(±3i/2) = 0), J₀ < 2 (generic k = 4) (HK §1, LRC §3.1, EESTEP §0) | Re-derived from the Stage-3 measure. Thresholds reproduced: **3/2, 3/2, 5/2 (2 generic), 3/2**. Margins under J₀ ≤ 1: ½, ½, 3/2, ½ | **yes** | [D], with the exponent inputs [DV] |
| A-iii.1 | R₋₃^{EE3}(−3i/2) = −524288, R₋₁^{EE3}(−3i/2) = 6193152/25, R₊₁(−3i/2) = 0, R₊₃ has a pole there (EESTEP §3, HKR §1.4) | **Exactly these values** (sympy, from `kernels.json`) | **yes** | [P] |
| A-iii.2 | f_pos(±3i/2) = 0 for Q's published continuum C4 weights, since χ = Q·B has no pole at ±3i/2 (HKR §1.4: [P]) | χ(±3i/2) is finite: −2.854·10⁻⁷ (C4v), +1.636·10⁻⁷ (C4d). (1/2πi)∮χ around ±3i/2 is **≤ 10⁻⁶⁸** (numerical zero), against 5·10⁻⁹ – 10⁻² around the genuine poles ±i/2, ±5i/2, ±2.914i. 1/N_com(±3i/2) = 0.0905415, regular and non-zero. The F₋₁ decay rate tends to 2.9 (χ's nearest lower singularity, −2.914i), against **1.500** for a control with a pole inserted at ±3i/2 | **yes** | [P] (by construction) + [NE], independently reproduced |
| A-iii.3 | Stage-7 XP grid weights: F₋₁ rate 2.91 / 2.87; implied residue at −3i/2 ≤ 3.4·10⁻²⁷ (value) and 2.6·10⁻²⁵ (deriv) at t = 36 (EESTEP §3) | Rates 2.835 → **2.907** (value) and 2.555 → **2.865** (deriv); implied residue **3.39·10⁻²⁷** and **2.57·10⁻²⁵** at t = 36 | **yes**, digit for digit | [NE], independently reproduced |
| A-iii.4 | (side) the continuum weights satisfy W[χ] = Φ, with \|W − Φ\| ≤ ε_W (2.7·10⁻¹³ / 9.1·10⁻¹¹) | \|W − Φ\| ≤ 3·10⁻²⁰ (t ≤ 3), 3·10⁻¹⁷ (t = 8), 1.4·10⁻¹⁶ / 2·10⁻¹³ (t = 15) | yes | — |
| B.1 | Lemma Y: max(0, σ) ≤ m_pol for every monomial × [one ε] × D_{Δ′}; **997 836 structures, worst margin 0** (R4H §2.3) | Own enumeration. **997 836 = 35 637 × 28**: 4 539 parity-even and 31 098 non-vanishing one-ε Z-leg structures (2 118 ε's containing all of P₁, Z₁, P₃, Z₃ vanish identically), times the 28 P·P monomials of total degree ≤ 2. **Worst margin 0, no violation.** The same holds for the 3⁶ convention (2.6·10⁷ structures) and, via the closed form below, for every P·P content | **yes**, count identical | **[P]**, independently reproduced |
| B.2 | Closed form behind Lemma Y (R4H §2.3 proof) | σ = (g + e + f) − (z + e′ + f′) + (k₁₃ − k₂₄)/2, checked on all 10⁶ structures (0 mismatches). Hence σ ≤ m_pol (no ε), σ ≤ m_pol − ½ (k₁₃ ≤ 2), σ ≤ m_pol + ½ with ε = O(y^{½}) (k₁₃ = 3), ε ≡ 0 (k₁₃ = 4) | **yes** (cleaner bookkeeping than R4H's displayed line) | [P] |
| B.3 | ε y-orders: k₁₃ = 3 → y^{0.500}; k₁₃ = 4 → exactly 0; k₁₃ ≤ 2 → O(1) (R4H §2.3) | k₁₃ = 3: **0.5000** in all 24 placements; k₁₃ = 4: \|ε\| ≤ 4·10⁻¹⁸ (identically 0); k₁₃ ≤ 2: 0.000 or 0.500 (the latter when both Z₁, Z₃ are present) | **yes** | [P + DV] |
| B.4 | Conclusion at the correlator level: every term of a spin-2 contact is O(y⁰ log y) (R4H §2.3) | Five explicit bulk contacts of four spin-2, Δ = 3 fields (all index contractions of ∏(V_i·V_j)). Reduced by my own Schwinger reduction to 833 "monomial × D_{Δ′}" terms each, D-functions evaluated **without Mellin–Barnes**. The individual D's blow up like **y^{−σ}** (σ up to 4; local orders −4.0000, −2.0000, −1.0000). The products stay bounded: the largest single term grows like log y, and the full contacts behave as **A + B log y** (3 vertices) or vanish (y^{0.92}, y^{1.9}) for y = 8·10⁻² … 1.2·10⁻⁶ | **yes** | [D + DV] (correlator level) |
| B.5 | L6′: \|L_{R⁴}\| ≤ CΔφ⁻²(1 + \|log Δφ\|), p = 2 < 7/2 (R4H §2.5) | Not reproduced at the level of the light function. The Ω₄ functional (t-contours, f₄, residue sums) was not re-implemented. What is checked: the counting (B.1–B.3) and the correlator-level bound (B.4). The extra y⁻¹ of f₂, f₄ follows from X₁₃ ∝ y, which was checked: X₁₃/(2t₁t₃y) = 1 to 10⁻¹⁴ | partly | **[D]** unchanged; its counting input is now [P, independently reproduced] |
| B.6 | Non-dispersive counterexample: local slope 4.7 → 3.5 (R4H §2.6) | **Not attempted.** It needs LIGHTSIDE_DIST's product-family Mellin machinery (σ ⊗ T_free Wick classes, residue tails), which is not cheap to rebuild independently | — | [DV]/[NE] unchanged (not independently checked) |

**Net.**
* **Target A: reproduced.** Covered: the EE = Σ_b shifted-scalar identity; the edge resummation and de-shift survival that feed the
  corner count; the J₀ thresholds 3/2 (k = 2, 3, 5) and 5/2 (k = 4); and f_pos(±3i/2) = 0 for Q's published weights.
* **Target B: reproduced at the level of Lemma Y and of the correlator.** The structure count (997 836) is identical. Explicit
  contacts confirm the O(y⁰ log y) conclusion. Not covered: the step from the correlator to the light function L_{R⁴}, and the
  non-dispersive counterexample.
* **No substantive discrepancy was found.** The list of discrepancies has one diagnostic difference and one documentation
  inconsistency, neither load-bearing (§5).

---

## 1. Target A, part (i): the EE reduction (`A1_ee_identity.py`, `cblocks.py`, `A1b_negb_mp.py`)

### 1.1 What was built, from the statements only
* **Frame.** H7 §18.1 / LRC §1.1, metric (−,+,+), |x| = 1, |y| = r.
  * ŷ = (1,0,0) and x̂ = (cosh β, sinh β, 0).
  * x₁ = −rŷ, x₃ = rŷ, x₀ = r n̂_θ, x₂ = −x̂, x₄ = x̂.
  * The celestial point is z₀ = (1, cos φ, sin φ).
  * Check: X₂₀ = −(1−r²)(1 − s cos θ), with s = 2r sinh β/(1−r²), as stated.
* **Scalar integrand.** I read off LRC §1.1 that `Nn ∝ V₀,₄₂X₂₄`, `m_ij = |X_ij|`, `c₁₀ = |z₀·x₁₀|`, `c₃₀ = |z₀·x₃₀|`, and that
  `c₁₃ = |x₁₀·x₃₀| = 2r²` is θ-independent. Then

      S_ν = |V₀,₄₂|^l (|X₂₄|/|X₂₀X₄₀|)^{(1+l)/2} c₁₃^{5/4+iν/2} (r² sin²(φ−θ))^{−δ/2},   l = −½+iν,  δ = ½+iν,

  with V₀,₄₂ the CPPR structure. This identification is confirmed a posteriori: the resulting kernel is exactly N(ν) r² Ĝ_{1,l}
  (b = 0 row of §1.3).
* **Dressing.** `P = τ₁₃D` with:
  * `τ₁₃ = (w₁·x₁₀)²(w₃·x₃₀)²`;
  * `D = [(H₂₄ − 2V₂,₀₄V₄,₀₂)/(X₂₀X₄₀)]²`, built from embedding-space CPPR invariants;
  * w₁, w₃ null in the rest frame of ŷ, and z₂, z₄ null and boosted with x̂;
  * three random polarisation sets.
* **Quadrature.**
  * θ: midpoint rule with 384 nodes.
  * φ: double-exponential rule on each half-circle between the singular points φ = θ, θ + π, where the integrand goes like
    |sin(φ−θ)|^{−½−iν}.
  * Self-convergence 2·10⁻¹³ at ρ = 0.85.
* **Blocks (own derivation).**
  * I wrote the Dolan–Osborn ½-Casimir with a = 0 and general b in radial variables ρ = re^β, ρ̄ = re^{−β}, using Euler operators.
    For d = 3, after multiplying by 2(1−r²)Den with Den = (1+r²)² − 4r²η²:

        (1−r²)Den(D_r² + Λ − C)G − [2(1−r²)²((1+r²) + 4brη) + Den(1+r²)] D_rG − 8r(1−r²)[rη + b(1+r²)](η²−1)∂_ηG = 0,

    where D_r = r∂_r, Λ = ∂_η(η²−1)∂_η and C = Δ(Δ−3) + l(l+1).
  * Checks: the DO operator reproduces the d = 4 closed-form (0, b) block exactly. The radial form equals the DO operator on a test
    function to 10⁻³¹.
  * Expansion: G = r Σ_N r^N Σ_m B[N,m] P_{l+m}(η). The leading eigenvalue N(N−1) + m(m + 2iν) is the one quoted in LRC §2.
  * P_{l+m} is obtained by upward recursion, using P_{l−m} = conj P_{l+m}; it is rescaled by r^{|m|} at large η.

### 1.2 Exact reduction
* P_even(X₂₀X₄₀)⁴ is a polynomial of degree 8 in cos θ: the Chebyshev tail is ≤ 5·10⁻¹⁵.
* The Laurent form Σ_{b=−4}^{4} c_b t^b reproduces P_even to 8·10⁻¹³.
* Φ(θ) = ∫dφ S_ν is even in θ, so the sin-odd part of P drops out, although pointwise it is as large as the even part.
* **Σ_b c_bΨ^{(b)} / K^EE − 1 ≤ 2.6·10⁻¹²** on 81 cases: ρ ∈ {0.25, 0.55, 0.85} × β ∈ {0.6, 1.4, 2.2} × ν ∈ {0.4, 1.3, 2.7} × 3 polarisation
  sets. The b-sum cancels by up to 6·10⁴.

### 1.3 Identification of Ψ^{(b)} (cache/A1_ee_identity.txt)
The table gives the spread over the 9 points of Ψ^{(b)}/(r² v^{q} Ĝ_{b′}), divided by its mean.

| ν | b | claim: b′ = −b, q = −b/2 | Euler form: b′ = \|b\|, q = \|b\|/2 | control: b′ = +b, q = −b/2 | control: scalar block |
|---|---|---|---|---|---|
| 0.4 | −4 … −1 | 8.7·10⁻¹⁰ … 1.4·10⁻¹⁴ | same | 3.7 – 7.0 | 1.7 – 2.3 |
| 0.4 | 0 | 2.8·10⁻¹⁵ | — | — | — |
| 0.4 | +1 … +4 | 5·10⁻¹⁴, 6·10⁻⁹, 2·10⁻⁴, 0.73 (direct series, see §5 item 1) | **1.3·10⁻¹⁴ … 8.7·10⁻¹⁰** | 1.7 – 3.2 | 0.7 – 1.4 |
| 1.3 | all | as above, ≤ 9.4·10⁻¹⁰ (Euler form for b > 0) | | 1.7 – 7.0 | 0.9 – 8.9 |
| 2.7 | all | ≤ 1.1·10⁻⁹ (Euler form for b > 0) | | 1.7 – 7.0 | 1.0 – 22 |

* **Prefactor fit** (ν = 1.3): log\|Ψ^{(b)}/Ĝ_{−b}\| = p log r + q log v + c gives p = 2.00000000 and q = 2.0000, 1.5000, 0.5000, 0.0000
  for b = −4, −3, −1, 0. The maximum residual is 4·10⁻¹⁰.
* Ψ^{(b)} = Ψ^{(−b)} to ≤ 8·10⁻¹⁴.
* So the DO sign convention is b_DO = −b for the insertion t^b = (X₂₀/X₄₀)^b, as EESTEP states.

### 1.4 End to end
N(ν) is taken from Ψ^{(0)} at (ρ, β) = (0.25, 0.6).

| polarisation set | max \|K^EE/pred − 1\| (27 cases) | b-sum cancellation | control: b = 0 block for every b |
|---|---|---|---|
| 0 | **1.1·10⁻⁹** | 9.1·10² | 1.20 |
| 1 | **6.8·10⁻⁹** | 9.6·10³ | 4.5·10³ |
| 2 | **4.7·10⁻⁸** | 6.0·10⁴ | 34 |

The original reports 1.6·10⁻⁷ – 1.5·10⁻⁶. The difference is quadrature quality only.

---

## 2. Target A, part (ii): corner power counting (`A2_corner.py`, `A2b_edge_envelope.py`, `A2c_w4g_mp.py`)

### 2.1 Edge resummation
The table gives max_b |Ĝ_b − 2Re[c_l e^{lβ}F_b]| / (2|c_l e^{lβ}F_b|), for b = 0…4, with
F_b = (z/4ρ)^h ₂F₁(h, h+b; 2h; z) and h = (½+iν)/2.

| ν, ρ | β = 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|
| 0.7, 0.3 | 3.8·10⁻² | 2.3·10⁻⁴ | 1.1·10⁻⁵ | 1.5·10⁻⁷ | 2.8·10⁻⁹ |
| 0.7, 0.85 | 8.6·10⁻² | 1.3·10⁻³ | 1.9·10⁻⁵ | 4.9·10⁻⁷ | 2.9·10⁻⁸ |
| 3.1, 0.3 | 1.2·10⁻² | 1.8·10⁻⁴ | 4.4·10⁻⁶ | 1.0·10⁻⁷ | 2.3·10⁻⁹ |
| 3.1, 0.85 | 1.3·10⁻¹ | 2.2·10⁻³ | 4.0·10⁻⁵ | 6.5·10⁻⁷ | 4.1·10⁻⁸ |

* The decay is e^{−2β}.
* Replacing F_b by F_{−b} gives O(1). This fixes the sign convention and agrees with EESTEP.
* F_b has poles only where 2h ∈ −ℕ (Im ν = 3/2 + 2p). There c_l = Γ(l+½)/(√πΓ(l+1)) vanishes. Both facts were checked as part of the
  residue computation below.

### 2.2 De-shift survival (cache/A2_corner.txt)
* Set-up: S_b = ηⁿ∫₀³⁰ w Ĝ_b dν against S_LR = ηⁿ∫w P_l, with Gauss–Legendre (400 nodes) and N_max = 200. n = 1 for w2g and w2s;
  n = 3 for w4g and w4p.
* The entries are the mean local exponents d log|S_b − S_LR|/dβ over [8.5, 10] and [12, 13.5].

| weight | ρ | b = 0 | 1 | 2 | 3 | 4 | LR alone |
|---|---|---|---|---|---|---|---|
| w2g = ν²e^{−ν²/9} | 0.2, 0.5, 0.85 | −0.50/−0.50 | −0.50/−0.50 | −0.50/−0.50 | −0.50/−0.50 | −0.50/−0.50 | −0.50 |
| w2s = ν² sech(πν/2) | 0.2 – 0.85 | −0.22/−0.36 | −0.35/−0.41 | −0.33/−0.40 | −0.32/−0.40 | −0.31/−0.39 | −0.37/−0.41 |
| w4g (regular at 3i/2) | 0.5 | −0.52/floor | −0.52/floor | −0.50/floor | −0.50/floor | −0.50/floor | −0.52/floor |
| w4p = w4g/(ν²+9/4) | 0.2, 0.5, 0.85 | **+1.00/+1.00** | +1.00/+1.00 | +1.00/+1.00 | +1.00/+1.00 | +1.00/+1.00 | −0.49/floor |

* **"floor".** The double-precision floor, where the late-window rate equals the raw rate +2.5. At ρ = 0.2 and 0.85, single entries
  in the early window range from −0.31 to −0.86 because the floor sets in early there. This is the same limitation the original
  reports (EESTEP: "floor from β ≈ 7 at ρ = 0.2").
* **The w4p residue at 3i/2.** (S_b − S_LR)/(η³Re[2πi Res_{3i/2} w c_l e^{lβ}F_b]) at β = 10 is 0.9999919 – 1.0000069 for every b and ρ.
  For w4g the same residue is zero, since w is regular there.

### 2.3 w4g beyond the floor (45-digit arithmetic, ρ = 0.5)
The same w4g integral was run in 45-digit arithmetic: Gauss–Legendre with 160 nodes on [0, 24], N_max = 110.

| b | \|S_b − S_LR\| at β = 9, 10, 11, 12, 13 | local exponents |
|---|---|---|
| 0 | 0.1392, 0.08425, 0.05106, 0.03096, 0.01878 | **−0.502, −0.501, −0.500, −0.500** |
| 2 | 21.24, 12.89, 7.817, 4.742, 2.876 | **−0.500, −0.500, −0.500, −0.500** |

Neither this run nor the original double-precision ones could see past β ≈ 10. Here the de-shift cancellation is seen to persist to
β = 13 at the LR target rate −½. In particular there is no η^{+½} growth from the s = 1 shell at these β. So the bound c = ½ used for
J₀ < 5/2 is conservative (with c = −½ the criterion would be J₀ < 7/2, which HK calls "NE" and does not claim).

### 2.4 The count
* **Measure.** I used the Stage-3 measure (RESULT_Q_STAGE3 §1.5). The kernel obeys |K| ≲ r^a η^c on the corner {r ≤ c₀/2η}, and K1′
  gives |dDisc 𝒢| ≤ B r^{1−J₀}. Then

      ∫^∞dη η^c ∫₀^{c₀/2η} r^{a+1−J₀} dr ∝ ∫dη η^{c−a−2+J₀},

  which converges iff **c < a + 1 − J₀**.
* **Numerical check.** The closed-form r-integral, integrated in η up to 10⁴ and 10⁸, converges at J₀ = threshold − 0.1 and grows at
  threshold + 0.1 for every row below.

| piece | a | c (source of the exponent) | threshold | margin under J₀ ≤ 1 |
|---|---|---|---|---|
| k = 2, Ĉ′_EE (de-shifted, cut-off) | 0 | −½ (§2.2, w2g; LRC\*) | **J₀ < 3/2** | ½ |
| k = 3, Ĉ₃ | 0 | −½ (P_{−½}, log) | **J₀ < 3/2** | ½ |
| k = 4, Ĉ₄^{ds}, ∂Ĉ₄^{ds} with f_pos(±3i/2) = 0 | 2 | +½ (the s = 1 shell's raw η^{½}; an upper bound, the measured exponent is −½) | **J₀ < 5/2** | 3/2 |
| k = 4 with f_pos(3i/2) ≠ 0 | 2 | +1 (§2.2, w4p; residue formula) | J₀ < 2 | 1 |
| k = 5, Ĉ₅^{ds} (raw) | 2 | +3/2 (P_{3/2}) | **J₀ < 3/2** | ½ |

* All four criteria quoted in the task are reproduced.
* The k = 4 value 5/2 rests on f_pos(±3i/2) = 0 (part (iii)) and on the s = 1 shell bound c ≤ ½. That bound is conservative: the
  exponent measured with my blocks is −½ wherever precision allows (§2.2, §2.3).
* Not checked independently: K1′, K2ν, LRC-a, and the position-to-heavy normalisation (HKR §1). These are separate inputs of H_K.

---

## 3. Target A, part (iii): f_pos(±3i/2) for Q's published weights (`A3_fpos.py`, `A3b_fpos_tail.py`)

### 3.1 Exact values (sympy, `kernels.json`)
* R₋₃(−3i/2) = −524288.
* R₋₁(−3i/2) = 6193152/25.
* R₊₁(−3i/2) = 0.
* R₊₃ has a simple pole at −3i/2 (factor 2ν + 3i in the denominator).

All four are as stated.

### 3.2 Continuum weights
* **Evaluation.** χ = Q·B, with B the 401-term Gaussian sum, evaluated in 100-digit arithmetic as the README requires.

| | C4v (value) | C4d (derivative) |
|---|---|---|
| χ(±3i/2) | −2.85440177·10⁻⁷ (finite) | +1.63575176·10⁻⁷ (finite) |
| (1/2πi)∮χ, \|ν ∓ 3i/2\| = 0.1 | ≤ 2.8·10⁻⁶⁹ | ≤ 1.2·10⁻⁶⁸ |
| control: around +i/2 / +5i/2 / +2.914i (genuine poles of Q) | 1.7·10⁻⁸ / 2.8·10⁻⁵ / 7.8·10⁻⁵ | 4.8·10⁻⁹ / 1.2·10⁻⁴ / 1.2·10⁻² |
| F₋₁ local decay rate, t = 4 → 12 | 1.72, 2.64, 2.76, 2.82 (→ 2.914) | 2.78, 3.66, 1.55, 2.44 (interference, see below) |
| control χ(ν²+4)/(ν²+9/4) | 1.508, 1.497, **1.500, 1.500** | 1.646, 1.519, **1.500, 1.500** |
| \|W − Φ\| at t = 0.5 / 3 / 8 / 15 | 3·10⁻²⁰ / 1·10⁻²⁰ / 6·10⁻²¹ / 1.4·10⁻¹⁶ | 2·10⁻²¹ / 2·10⁻²⁰ / 3·10⁻¹⁷ / 2·10⁻¹³ |

* **The ν-normalisation at ±3i/2.** 1/N_com(±3i/2) = 0.0905414787, and n(3i/2) = 1024/N_com = 92.714474 (HKR: 92.714…). It is regular
  and non-zero, so a pole of the position weight w = χ·const/N_com at ±3i/2 would have to be a pole of χ. There is none.
  **f_pos(±3i/2) = 0.**
* **The C4d tail.** F₋₁ is the sum of an e^{−2.914t} term and a t²e^{−3t} term (the triple pole at −3i), so its local rate oscillates.
  A3b follows |F₋₁|e^{3t/2} out to t = 22, where a pole at −3i/2 would make it tend to a non-zero constant:

  | | t = 6 | 10 | 14 | 18 | 22 | mean rate over [18, 22] |
  |---|---|---|---|---|---|---|
  | C4v | 2.7·10⁻³ | 2.2·10⁻⁵ | 1.1·10⁻⁷ | 4.5·10⁻¹⁰ | **1.8·10⁻¹²** | 2.886 |
  | C4v control (pole inserted) | 0.258 | 0.259 | 0.259 | 0.259 | **0.259** | 1.500 |
  | C4d | 6.7·10⁻³ | 8.1·10⁻⁵ | 1.3·10⁻⁶ | 1.05·10⁻⁸ | **6.2·10⁻¹¹** | 2.783 |
  | C4d control (pole inserted) | 0.154 | 0.148 | 0.149 | 0.149 | **0.149** | 1.500 |

  The weights decay towards the 2.914 rate (the C4d approach is slower, because of the t²e^{−3t} term). The controls sit exactly on a
  constant, so they have the e^{−3t/2} pole term the weights lack. A residue of χ at −3i/2 is excluded at the level of
  \|A\| ≲ 10⁻¹² / R₋₁(−3i/2) ≈ 10⁻¹⁵ (C4v) and ≈ 10⁻¹³ (C4d) from the t = 22 values alone. The contour integrals bound it far more
  tightly (≤ 10⁻⁶⁸).

### 3.3 Stage-7 XP grid weights
The test is the midpoint sum over ν_k = (k+½)/25. Out to t = 36 it reproduces EESTEP §3's numbers digit for digit:

| weight | local rates, t = 10 → 36 | implied residue \|A(36)\| |
|---|---|---|
| C4_value | 2.835, 2.876, 2.893, 2.901, **2.907** | **3.39·10⁻²⁷** (EESTEP: 3.4·10⁻²⁷) |
| C4_deriv | 2.555, 2.743, 2.808, 2.843, **2.865** | **2.57·10⁻²⁵** (EESTEP: 2.6·10⁻²⁵) |

**Caveat.** The step "f_pos ≠ 0 ⇔ χ has a pole" uses the normalisation w = χ·const/N_com with a ν-independent const. HKR §1 checks
that normalisation through the heavy image of the EE k = 4 kernel. I did not redo that check. It is the one link of (iii) that rests
only on the original authors' work.

---

## 4. Target B: Lemma Y (`B1_lemmaY_enum.py`, `dfun.py`, `B2_contact.py`)

### 4.1 Enumeration (cache/B1_lemmaY_enum.txt)
* **What is enumerated.**
  * Monomials in Z_i·Z_j (i < j), Z_i·P_j (i ≠ j) and P_i·P_j (i < j), with each Z_i of degree exactly 2. Z_i², Z_i·P_i and P_i² are
    zero.
  * Optionally one ε of five *distinct* vectors among {P_i, Z_i}. A Z inside the ε uses one of its two legs.
* **Definitions.**
  * deg P_i counts the monomial and the ε slots; σ = (deg P₁ + deg P₃ − deg P₂ − deg P₄)/2.
  * m_pol counts the factors among {Z₁·Z₃, Z₁·P₃, Z₃·P₁, P₁·P₃}, plus ½ when the ε holds exactly three of {P₁, Z₁, P₃, Z₃}.
  * An ε holding all four of them vanishes identically, since P₁ − P₃ = −t₁Z₁ − t₃Z₃.

| convention | Z-leg structures | × P·P monomials | total | non-vanishing | worst max(0,σ) − m_pol | violations |
|---|---|---|---|---|---|---|
| no P·P | 4 539 even + 33 216 one-ε | 1 | 37 755 | 35 637 | 0 | 0 |
| **P·P total degree ≤ 2** | same | **28** | 1 057 140 | **997 836** | **0** | **0** |
| each P·P multiplicity ≤ 2 | same | 729 | 27 523 395 | 25 979 373 | 0 | 0 |

* R4H's "997 836 structures" is therefore exactly: all non-vanishing Z-leg structures (parity-even plus one-ε) times all P·P monomials
  of total degree ≤ 2.
* **Distribution.** The parity-even margins max(0,σ) − m_pol (no P·P) are: 0 for 2 082 structures, −1 for 1 862, −2 for 555, −3 for 39
  and −4 for 1.
* **Closed form** (all structures, 0 mismatches):

      σ = (g + e + f) − (z + e′ + f′) + (k₁₃ − k₂₄)/2,

  * g, e, f count Z₁·Z₃, Z₁·P₃ + Z₃·P₁, and P₁·P₃;
  * z, e′, f′ count Z₂·Z₄, Z₂·P₄ + Z₄·P₂, and P₂·P₄;
  * k₁₃ and k₂₄ count the ε slots on the 1–3 and 2–4 legs, with k₁₃ + k₂₄ = 5 when an ε is present.

  So Lemma Y holds for **any** number of P·P factors, not only for the enumerated truncation: P₁·P₃ raises σ and m_pol together, the
  other P·P factors cannot raise σ. The margin 0 is attained exactly when z = e′ = f′ = 0 (and k₁₃ = 3 for the ε terms).
* **Proof check.** The proof displayed in R4H §2.3 has a garbled intermediate line ("… − e₂₄ − …"). Its conclusion and its case split
  on k₁₃ are correct, and the closed form above is a clean replacement.

### 4.2 Explicit contacts (cache/B2_contact.txt)
* **Configuration.** 3d metric (+,−,−).
  * x₁ = −t₁z(a), x₃ = t₃z(b), x₂ = e₀, x₄ = −e₀, with t₁ = 0.3, t₃ = 0.2.
  * Polarisations z₁ = z(a), z₃ = z(b), and fixed null z₂, z₄; mean angle 0.4.
  * Then X_ij > 0 and X₁₃ = 2t₁t₃y, checked to 10⁻¹⁴. Z₁·Z₃, Z₁·P₃, Z₃·P₁ and P₁·P₃ are exactly ∝ y.
* **Contacts.** Each spin-2 Δ = 3 leg is the bulk-to-boundary propagator (V_i·W)²/W_i⁵, with V_i = W_iZ_i + 2(Z_i·X)P_i. Five quartic
  vertices: (V₁·V₂)(V₂·V₃)(V₃·V₄)(V₄·V₁), (V₁·V₃)(V₃·V₂)(V₂·V₄)(V₄·V₁), and (V_i·V_j)²(V_k·V_l) for the three pairings.
* **Reduction.** I used the Schwinger identity

      ∫dX ∏U_i^{n_i}W_i^{−b_i} ∝ Γ(Σ(b−n)/2 − h)/∏Γ(b_i) ∫∏dt t_i^{b_i−1} ∏_i(Σ_{j≠i} t_j Z_i·P_j)^{n_i} e^{−Σt_it_jX_ij},

  which gives Σ (monomial) × I_{Δ′}, with 833 distinct Δ′ per vertex and σ up to 4.
* **D-functions, without Mellin–Barnes.**
  * t₂ and t₄ are integrated in closed form (Tricomi U), and the radial integral in closed form (DLMF 13.10.7 → ₂F₁).
  * What remains is a 1-d integral, done in 30-digit arithmetic.
  * Checked against the 3-d Feynman-parameter integral to 10⁻¹¹ – 10⁻¹⁵, and against the box function (I₁₁₁₁X₁₃X₂₄ = Φ/2).
* **Single D-functions blow up as claimed.** Local y-orders for Δ′ = (7,3,7,3), (5,3,5,3), (4,3,4,3), (6,4,5,3): −4.0000, −2.0000,
  −1.0000, −2.0000, i.e. −σ.
* **The contacts do not blow up.** Over y = 7.9·10⁻² → 1.2·10⁻⁶:

| vertex | G(y) at y = 1.2·10⁻⁶ | behaviour | largest single term |
|---|---|---|---|
| (12)(23)(34)(41) | −3.69·10⁻³ | A + B log y, B = 4.9·10⁻⁴ | 0.25 → 2.18 (∝ log y) |
| (13)(32)(24)(41) | −2.9·10⁻⁸ | → 0 (local order +0.92) | 1.5·10⁻² → 2·10⁻⁶ |
| (12)²(34)² | −2.08·10⁻² | A + B log y, B = 2.13·10⁻³ | 0.83 → 8.3 (∝ log y) |
| (13)²(24)² | +3.2·10⁻¹³ | → 0 (local order +1.9) | 1.5·10⁻³ → 2.8·10⁻¹² |
| (14)²(23)² | +2.08·10⁻² | A + B log y, B = −1.55·10⁻³ | 0.19 → 1.70 (∝ log y) |

* The individual D_{Δ′} reach y⁻⁴, but every monomial × D term is O(log y), with margin 0 attained ((7,3,7,3) against (Z·P)⁴-type
  monomials). This is the correlator-level content of Lemma Y. **Confirmed.**
* **What this does not test.** The light function L_{R⁴} = Ω₄[X](a, b) needs CLS-D's og-gauge t-contours, f₄ = y⁻¹F₄ + yF₄′, and the
  finite residue sum. Those were not re-implemented, so the final step to p = 2 (R4H §2.4–2.5) is not independently checked. Neither is
  the non-dispersive counterexample (R4H §2.6 (a), (d)).

---

## 5. Discrepancies and remarks

1. **(Diagnostic, immaterial.) Negative-b radial series.**
   * EESTEP §1.3 says the radial series "loses up to 9 digits" for b = −4 at ρ = 0.85 and switches to the Euler identity.
   * I find that at N = 200 the direct series for Ĝ_{−4} at ρ = 0.85 is wrong by O(1) **also in 50-digit arithmetic**: |Ĝ_{−4}/(v⁴Ĝ₄) − 1|
     = 152 at β = 0.6 and 2.9 at β = 2.2. So the dominant problem is truncation: the coefficients grow like N^{2|b|−1}, from the
     (1+ρ)^{−4|b|} pole of v^{|b|} at ρ = −1. Double-precision round-off contributes only 5·10⁻⁹.
   * EESTEP evaluates negative b only through the Euler identity, which I confirm to 10⁻³⁸ – 10⁻⁴⁵. So no number of theirs is affected.
   * The wording should be "the negative-b series converges too slowly (and loses digits)".
2. **(Documentation, not load-bearing.) Raw k = 2 entry of the corner table.**
   * HK_CLOSE §1 lists raw Ĉ₂ per ν with c = ½ ("P_{l+1}"), giving J₀ < ½.
   * RESULT_Q_STAGE3 §1.5 lists raw Ĉ₂ with κ = ηP_{l±1} ~ η^{3/2}, giving J₀ < −½. The two differ by one power of η (the f-factor).
   * Both say "fails". Nothing uses the raw k = 2 kernel, and the de-shifted rows (c = −½, J₀ < 3/2) agree in both. Worth aligning in
     the write-up the referee asked for.
3. **(Proof text.) R4H §2.3.** The displayed chain for 2σ contains a garbled step. The statement and case analysis are right; §4.1 gives a
   clean closed form.
4. **Not independently reproduced, and still resting only on the original authors' work:**
   * the position-to-heavy normalisation of the EE k = 4 kernel (HKR §1), needed for "f_pos = 0 ⇔ χ regular";
   * K1′, K2ν, LRC-a;
   * the light-function step of L6′, i.e. the t-contours, f₄ and the residue sums (R4H §2.4–2.5);
   * the non-dispersive counterexample, slopes 4.7 → 3.5 (R4H §2.6).

**Suggested grade wording for the paper.**
* "H_K's EE reduction, LRC\* exponents, J₀ thresholds and f_pos(±3i/2) = 0: **independently reproduced** (INDEP_CHECK §1–3)."
* "Lemma Y (count and closed form) [P], **independently reproduced** including the structure count 997 836; correlator-level O(y⁰ log y)
  confirmed on explicit spin-2 contacts. L6′ remains [D] (the functional step was not re-checked)."

---

## 6. Reproduce

    cd ads4/code/indep
    ulimit -v 3000000; export OMP_NUM_THREADS=1
    python3 B1_lemmaY_enum.py        # 25 s   -> cache/B1_lemmaY_enum.txt
    python3 dfun.py                  # 60 s   self-test of the D-function evaluator
    python3 B2_contact.py            # 13 min -> cache/B2_contact.txt
    python3 cblocks.py               # 5 s    operator checks
    python3 A1_ee_identity.py        # 10 s   -> cache/A1_ee_identity.txt
    python3 A1b_negb_mp.py           # 6 min  -> cache/A1b_negb_mp.txt
    python3 A2_corner.py             # 75 s   -> cache/A2_corner.txt
    python3 A2b_edge_envelope.py     # 20 s   -> cache/A2b_edge_envelope.txt
    python3 A2c_w4g_mp.py            # 24 min  -> cache/A2c_w4g_mp.txt
    python3 A3_fpos.py               # 2.5 min -> cache/A3_fpos.txt
    python3 A3b_fpos_tail.py         # 6 min  -> cache/A3b_fpos_tail.txt

| file | content |
|---|---|
| `cblocks.py` | DO Casimir (a = 0, b) in radial variables (own derivation), checks against DO / d = 4 blocks; continuous-spin recursion (scalar and ν-vectorised); scaled Legendre functions |
| `A1_ee_identity.py` | frame, scalar integrand, EE dressing, DE/midpoint quadrature; Laurent reduction; identification; end to end |
| `A1b_negb_mp.py` | 50-digit recursion; Euler identity; diagnosis of the negative-b series |
| `A2_corner.py`, `A2b_edge_envelope.py`, `A2c_w4g_mp.py` | edge resummation; de-shift exponents per b-block; w4p residue formula; corner count; w4g beyond the floor |
| `A3_fpos.py`, `A3b_fpos_tail.py` | exact R_j at −3i/2; continuum and XP weights near ±3i/2; F₋₁ decay with controls; W = Φ |
| `B1_lemmaY_enum.py` | independent structure enumeration, σ and m_pol, closed form |
| `dfun.py` | Symanzik D-function by closed-form t₂, t₄, s integrations + 1-d quadrature (no Mellin–Barnes), self-tests |
| `B2_contact.py` | ε orders; explicit spin-2 contacts reduced to monomial × D and scanned to y = 10⁻⁶ |
