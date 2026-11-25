# FULLPOS_CERT: certified angular moments and every-mass positivity of Ω^EE at fixed spin (Δ_g = 1000)

*26 Nov 2026. Research agent, about 5 h. Unreviewed, not committed. No existing file was edited, and nothing under `prl/` was touched.*
* *New code:* `ads4/code/fullpos/`:
  * `fp_common.py` (data, CDH polynomials, ladder data);
  * `fp_moments.py` (certified moments, closed form);
  * `fp_poly.py` (Bernstein every-mass certificate per J);
  * `fp_norm.py` (fixed-r normalisation identities);
  * `fp_compare.py` and `fp_ratio.py` (H_ang vs T0's hybrid);
  * `fp_regions.py` (coefficient profile, comparison constant);
  * `fp_kernel.py` (finite-mass frequency kernel, §31 identity).
* *Outputs:* `fullpos/cache/` (its `.gitignore` ignores everything but itself).
* *Read-only inputs:*
  * data: `prl/data/k2/k2_continuum_D_1000_rep.json`, `prl/data/nf02_psi.json`, `prl/data/kernels.json`;
  * code: `largespin/ls_core.py`, which reaches `spinfunc/` (T0's hybrid). It is used only in the comparisons.
* *Compute:* about 1.5 CPU-h in total, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

**Grades.**

| grade | meaning |
|---|---|
| [certified] | Arb ball arithmetic or exact rationals, with every remainder proved |
| [DV] | numerical (mpmath or float) at high precision, or sampled |
| [D] | derived, not machine-checked |

**d = 1 bug (coordinator item 1): confirmed fixed.**
* `cdh1` divides by (n+4)² in all three files. With that, q₁ = (15−4x)/64.
* Every k = 3 quantity in this document was produced after the fix:
  * Q_k is computed afresh by `fp_poly.Q_exact` from `fp_moments.cdh1_fmpq`, in every certificate, comparison and diagnostic run;
  * the CDH identity I3 (§2) checks the fixed recurrence to 40 digits.
* The first K = 300 moment file was written by a process that had loaded the old `cdh1`. Its stale `Q` fields were deleted, and nothing ever read
  them. Z_k and Y_k (d = 0) are unaffected.

---

## 0. Verdict

**Short answer: not close to a theorem for all spins.** It is certified spin by spin up to J ≈ 240, and the routes proposed for large spin
do not work as formulated. Every-mass positivity is proved for every spin up to a computable J_max. Beyond that, no proposed asymptotic route
closes. The obstruction is structural: the angular moments are exponentially large and cancel.

| item | result | grade |
|---|---|---|
| **1. Moments** Y_k, Z_k (repaired weight, Δ_g = 1000 and 300), Q_k | Closed form (exact Gaussian moments, Faddeeva, an integration-by-parts recurrence for the log terms) plus a rigorous θ treatment. Δ_g = 1000: k ≤ 300 at relative radius ≤ 2·10⁻¹³⁴ (75 s). Δ_g = 300: k ≤ 300. Q_k exact rationals | **[certified]** |
| profile of the moments | \|Z_k\| grows to ~10⁶⁷⁰ at k ≈ Λ/2 (Δ_g = 1000). At Δ_g = 300 the certified moments show the predicted **oscillating plateau** \|Z_k\| ≈ 10^{198.5} ≈ e^{(π/2)Λ} for 160 ≤ k ≤ 300 | [certified] numbers |
| **2. Normalisation** | Fixed-r identities I2/I3 hold to 10⁻³⁹ (I3 with the fixed d = 1 recurrence). Hence H_ang → 2∫ψ₂P + σu∫ψ₃P = T0's B₂+cv₂ + B₃+cv₃ limit, **relative factor exactly 1** (9216/4608 = 2 is the factor 2 of B₂). Z₀/3 = 2∫ψ₂ = 20.77189… [certified]. The finite-T deviation scales as 1/T | identities [DV, 40 digits]; ℬ itself not re-derived |
| finite mass vs hybrid | max \|H_ang − H_hyb\|/N_s = 4.1·10⁻³ on 75 states. J ≤ 26, sampled over all masses (37 channels): \|H_ang/H_hyb − 1\| ≤ **4.7·10⁻³ < 0.01**, min H_ang/N_s = 0.2009. **This independently confirms the authors' J ≤ 26 claim**, sampled, not certified | [DV] |
| **3. Every mass at fixed J** | **H_ang,b > 0 at every mass m ≥ m_gap for every physical channel with J ≤ 239 (Δ_g = 1000) and J ≤ 285 (Δ_g = 300)**. Bernstein/bisection on 𝒫_{J,b}(w), 356 + 426 polynomials. The limit is moment precision; no sign failure was seen. POLY_EXT | **[certified]** |
| **4a. Bounded b̂** | 𝒦_b(r₀) ≥ 1.8·10⁷, 7·10¹¹, 1.8·10²⁷, 3.9·10⁸⁴ at r₀ = 10⁻⁶, 10⁻⁴, 10⁻³, 10⁻². It grows like e^{2b̂}, so it is useful only for b̂ ≲ 20 at T ≥ 10⁶–10¹³: nothing beyond T0 Part F/T territory and the fixed-J certificates | [certified] partial sums |
| **4b. r → 1 / tail envelopes** (§§29, 31–33) | The envelopes hold for **no** k ≤ 300; \|P_k\|/P*_k up to 10⁵⁵⁸ (10²¹⁰ at Δ_g = 300). The threshold is intrinsic, log K ≳ Λ (coordinator and §33 agree: K ≈ e^{10⁵}). No explicit R₀ obtainable from coefficients | [certified] numbers + [D] |
| **4c. Pivot** to §29's signed remainder in r | **Does not close the middle region.** F_j(r) is a tail of the same huge series: \|F_j\| ~ \|P_{j+1}\|. The triangle majorant of H_ang − H_LR is 6·10¹⁶ (J = 30), 2.6·10⁶² (J = 100), 1.9·10¹²¹ (J = 200), against the true remainder 10⁻⁵–10⁻⁸ (`fp_pivot.py`). At r ≳ 0.03 even H_LR cannot be summed from coefficients at K = 300. The same e^{2b̂} cancellation as in 4a | [certified] numbers |
| §31/§32 identities | 2Y_k = ((k+3)²/2)L_k + ((k+4)²/2)L_{k+1}: exact (fmpq, k ≤ 60). §32's V-representation: exact to 10⁻¹⁵³ **provided u_{−1} = (b−1)b (= 12 for b = 4) is used at k = 0**. With u_{−1} = 0 it is wrong for b = 4 by 10 % | [certified] |
| forward strip (b̂ ≤ 1, J ≥ 27) | Its inputs \|P_k\| ≤ 200Λ^{2k} and \|D_k²Q_k\| ≤ 120Λ^{2k} hold for k ≤ 300 (maxima 10^{1.32}, 10^{1.11}). The strip argument itself was not reproduced | [certified] inputs only |

**How close "full AdS positivity of Ω^EE at Δ_g = 1000" is to a certified theorem (modulo item 4).**

*Certified:*
* every state with J ≤ 239, at every mass (H_ang,b > 0);
* T0 for the hybrid (all states), with the hybrid agreeing with H_ang to ≤ 4.7·10⁻³ (sampled) where tested.

*Missing:* every J ≥ 240 at finite b̂. That is an infinite region, so no finite extension of the per-J certificate closes it.

*Cost of extending the per-J certificate:*
* moments with ≈ 0.57J + 200 correct digits (polynomial cost, Õ(J²) ball operations at O(J) digits);
* J = 10³ is a few CPU-hours;
* J = 10⁴ is impractical.

*What I recommend for large J:* bound the **frequency kernel** 𝒦_fin(ν; J, T) of §4.4 against the hybrid's kernel per ν. Measured at
J ≤ 2000, it is O(ν²/T)-close to the Legendre kernel for ν ≪ √T, with no exponential sizes. Integrate that error against |g|. Do not
manipulate moments: every moment- or coefficient-based estimate (4a, 4b, the 4c pivot) loses e^{2b̂} or e^{(π/2)Λ} to the triangle
inequality.

---

## 1. Certified moments (task 1) [certified]

### 1.1 Definitions
Everything is as in the brief.
* f_m = Q⁽⁰⁾_{m−1}(ν²) (CDH, (a,b,c) = (3/2,½,½)); z_m and y_m as given.
* q_k = Q⁽¹⁾_k(ν²) (d = 1): q₀ = 1, (n+4)²q_{n+1} = −(s+4−2(n+2)²)q_n − n²q_{n−1}.
* Y_k = ∫g y_{k+2}, Z_k = ∫g z_{k+2}, with
  g = −9216(4x+1)/((4x+9)(4x+49))·χ_re·θ/Λ, x = ν², θ = ½[erf((ν+ν_c)/10) − erf((ν−ν_c)/10)], ν_c = Λ+100.
* Q_k = ∫₀¹ψ₃(p)q_k(Λ²p²)dp, with ψ₃ = NF02 `psi3_exact`/k₀. This is the k = 3 moment whose ladder sum is σ_b u∫ψ₃P (identity I3, §2), i.e.
  T0's (B₃+cv₃) limit.
* Every double of the JSON is taken exactly.

### 1.2 Method: closed form, no quadrature of the RBF part (`fp_moments.py`)
Write g = −(1/Λ)[R(x)B(ν) + R̃(x)L(ν)]θ, where
* R = 147456x/((4x+1)(4x+9)(4x+49));
* R̃ = 9216(4x+1)/((4x+9)(4x+49));
* B is the Gaussian RBF sum;
* L = Σ_j y_j ν^{2j}(x+1)log(1+x)e^{−x/2}.

For a polynomial p (exact `fmpq_poly`), the two pieces are handled separately.

**Gaussian piece.** Exact division gives p·147456x = q·(4x+1)(4x+9)(4x+49) + rem. Then
∫₀^∞pRB = Σ_i q_i S_{2i} + Σ_a c_a F_a, where:
* S_n = ½Σ_signed x_k m_n(c_k, σ_k) are exact Gaussian moments (recurrence m_{n+1} = c m_n + nσ²m_{n−1});
* F_a = ∫₀^∞B/(4ν²+a) = ½Σ x_k (π/4b)Re w((ib−c_k)/(σ_k√2)), with b = √a/2 and w the Faddeeva function (Arb `erfc`);
* c_a = rem(−a/4)/∏_{a′≠a}(a′−a).

**Log piece.** Division by (4x+9)(4x+49) reduces it to two families.
* K_n = ∫ν^{2n}log(1+ν²)e^{−ν²/2}, by the exact integration-by-parts recurrence K_{n+1} = (2n+1)K_n + 2E_{n+1}, with
  E_{n+1} = (2n−1)!!√(π/2) − E_n and E₀ = (π/2)e^{½}erfc(1/√2).
* Three one-dimensional integrals: K₀ and ∫L/(4ν²+a), a = 9, 49. These use Arb's rigorous `acb.integral` on [0, 64], plus an e^{−2000} tail.

**θ (the cutoff).** θ is not negligible for high moments: per moment, its effect is ~10⁻²⁴ × (1.03)^{2k} of ∫|g p|. It is included as follows.
* **Edge centres** (c ∈ [Λ−100, Λ+30]). The exact integral ∫x^iR B_edge(1−θ) over [Λ−150, Λ+180] uses Gauss–Legendre on 11 panels of
  width 30, 700 nodes each (300 in the K = 300 run). The error is rigorous: (64/15)Mρ^{−2n}/(ρ²−1)(L/2), Trefethen, ATAP Thm 19.3, with
  * ρ = 1.868;
  * M = sup|h| on the Bernstein ellipse, from 256 Arb boxes covering its boundary;
  * times Xmax^i.
* **All other centres, twins and the region outside the panels.** Bounded by completing the square with 1−θ ≤ e^{−(ν_c−ν)²/100}:
  * centres c ≤ Λ−100: e^{−394};
  * twins: e^{−2750};
  * outside the panels: e^{−625}.

  Plus the log piece, e^{−3000}.
* The authors' Gaussian-convolution recurrence for M_n = ∫ν^nG_c(1−θ) would replace the panel quadrature. It was not needed.

### 1.3 Results (Δ_g = 1000; `cache/moments_1000_K*.json`, full-precision balls)

| run | k range | precision | relative radius of Z_k, Y_k | time |
|---|---|---|---|---|
| K = 300 (NGL 300) | 0–300 | 3000 bits | ≤ 1·10⁻¹⁵⁶ (k ≤ 40); 6·10⁻¹⁴² (k = 240); 1.9·10⁻¹³⁴ (k = 300) | 75 s |
| K = 1000 (NGL 700) | 0–1000 | 7000 bits | **K1000_RADII** | **K1000_TIME** |

Low moments (balls exact to the digits shown):
* Z₀ = 62.31567904599110749883…; Y₀ = −756933.09038494103655…;
* Z₁ = −756964.24822446403210…; Y₁ = 12585746515.394383818…;
* Z₀/3 = 20.7718930153 = 2∫₀¹ψ₂dp. This is the J = 0, T → ∞ value; see §2.

Q_k are exact rationals.

### 1.4 The moments are exponentially large, and the large-k asymptotics are out of reach [certified profile + D]

| k | 0 | 10 | 50 | 100 | 200 | 300 | 400 | 500 |
|---|---|---|---|---|---|---|---|---|
| log₁₀\|Z_k\| | 1.8 | 40.5 | 157.5 | 267.3 | 429.0 | 545.7 | 626.2 | 670.4 |
| Y_k/Z_k | −1.2·10⁴ | −1.2·10⁴ | −4.3·10³ | −2.2·10³ | −1.19·10³ | −710 | −411 | −116 |

(`cache/fp_regions_1000_K300.txt`; k = 400 and 500 from the K = 1000 run.)

The signs of Z_k alternate in blocks; |Z_k| grows like (Λ^k/k!)². This is not a defect of the weight. It is forced by the CDH polynomials:
* For fixed ν, the recurrence is non-oscillatory for m < ν/2. There f_m(ν) grows to about e^{(π/2)ν}: ∫₀¹acosh(1/s)ds = π/2.
* Beyond m ≈ ν/2 it oscillates with amplitude ~ e^{(π/2)ν}m^{−3/2}·cos(ν log m + φ), and it never decays further.
  * [DV] ν = 20: |f_m|m^{3/2} ≈ 10¹¹ for m = 10 … 2·10⁵;
  * [DV] ν = 60: |f_m|m^{3/2} ≈ 10³⁷ for m = 100 … 2·10⁵.
* So Z_k ~ e^{O(Λ)} up to k ≈ Λ/2 ≈ 500 (observed peak ≈ 10⁶⁷⁰). Beyond that it is an oscillatory Fourier integral of g(ν)e^{(π/2)ν} at frequency log k.
* Because g has poles at ν = ±i/2 (from R), that Fourier integral decays only like k^{−1/2}.

**Consequences.**
* The authors' asymptotic Z_k ~ z_*k^{−5/2}/log k (and Y/Z → ¼) can dominate only for log k ≳ 3Λ, i.e. **k ≳ 10¹³⁰⁰** [D, heuristic].
* For every k reachable in practice, the coefficients are ≥ 10⁵⁰⁰ times the envelope:
  * |P_k/P*_k| = 10^{5.1}, 10^{140}, 10^{278}, 10^{558} at k = 0, 40, 100, 300;
  * P*_k = ((5/2)_k/k!)(A + B/log(k+2)).

This is the key negative finding for task 4; see §4.

## 2. Normalisation (task 2)

### 2.1 Fixed impact parameter: exact identities, no relative factor [DV to 40 digits]
`fp_norm.py` uses mpmath at 40 digits, sums to k = 1500, at ν ∈ {0.3, 0.7, 2, 5.5} × r ∈ {0.05, 0.3, 0.6}.

* **(I2)** (1−r)³Σ_k D_k²r^k z_{k+2}(ν)(1/(k+3) − r/(k+1)) = −(4x+9)(4x+49)/(4608(4x+1))·[R₋P_{−½+i(ν+i)}(X) + R₊P_{−½+i(ν−i)}(X)],
  with X = (1+r)/(1−r). Maximum relative deviation: **6.5·10⁻⁴⁰**.
* **(I3)** Σ_k D_k²r^k q_k(ν²) = P_{−½+iν}(X)/(1−r)⁴. Maximum relative deviation: **2.4·10⁻⁴⁰**. This uses the fixed (n+4)² recurrence; the old
  (n+3)² fails at O(1).

**What the identities give.** As T → ∞ at fixed r (U_{b,k} → r^k; the Y, V and Q terms are O(1/T)), the k = 2 action has the limit
* H_ang^{(2)} → ∫g·(I2) = (2/Λ)∫χ_reθ[R₋P(ν+i) + R₊P(ν−i)]dν.

This is exactly **2/Λ × the realised leading-Regge image** of the k = 2 weight. By the de-shift identity it is
* = 2∫₀¹ψ₂(p)P_{−½+iΛp}(cosh β)dp + δ, with |δ| ≤ 2ε_W P_{−½}/Λ ≤ 5.4·10⁻⁹ (EPSW_CERT, repaired: ε_W ≤ 2.68·10⁻⁶).

That is T0's lim(B₂ + cv₂): u → 0 and D_h → J₀, Mehler–Heine.

**The factor 2 is 9216/4608. There is no other factor.** The k = 3 part, Λ²Q_k/T summed with (1−r)³D_k²r^k, gives
(Λ²/(T(1−r)))∫ψ₃P = σ_b u∫ψ₃P, which is T0's (B₃ + cv₃) limit. Checks:
* At J = 0, T → ∞: H_ang = Z₀/3 = 20.77189301… = 2∫ψ₂dp (§1.3) [certified].
* At fixed r, T → ∞ [DV, `fp_compare.py`]:

  | | J = 20 | J = 40 (T four times larger) |
  |---|---|---|
  | b̂ = 2 | −1.1·10⁻⁴ | −2.8·10⁻⁵ |
  | b̂ = 6.3 | −1.2·10⁻³ | −3.0·10⁻⁴ |
  | b̂ = 20 | +4.3·10⁻² | +1.3·10⁻² |

  The entries are H_ang/(2Ω_LR,2 + uΩ_LR,3) − 1. Each scales like 1/T, so the limit has unit normalisation.

**The normalisation constant I use is 1. It is exact at the level of the leading-Regge limit:** identity I2 plus the de-shift identity
(certified ε_W). This is consistent with §28's ℬ = (2πGΛ⁴/9)T⁴/(T−C)³ and the congruence lift. I did not re-derive ℬ itself (the
Einstein calibration); that remains the authors' [P].

**Impact-parameter convention.** T0 uses sinh(β/2) = (J+½)/m, m² = T − C, i.e. r_β = (C+¼)/(T+¼). Comparing with T0's β-functions at
r_β instead of r reduces the finite-T discrepancy in every row of the fixed-r table (`cache/fp_compare_1000.txt`, column "(r_beta)").

### 2.2 Finite mass: 75 states [DV]
`fp_compare.py`, `cache/fp_compare_1000.txt`.
* J ∈ {0, 2, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 60}, both channels, twists 1000, 1300, 2000.
* H_ang is exact (certified moments, finite ladder sum in Arb). H_hyb and N_s are T0's `hyb_fast` with the target ψ; the repaired weight
  moves cv₂ by ≤ 5.4·10⁻⁹.

| | value |
|---|---|
| max \|H_ang − H_hyb\|/N_s, all 75 | **4.07·10⁻³** (J = 16, `+−`, twist 1000) |
| `++`, all twists | ≤ 4.7·10⁻⁵ |
| `+−`, twist ≥ 1300 | ≤ 4.2·10⁻⁵ |
| largest values | `+−` at the twist edge (u ≈ 1): 1.5–4.1·10⁻³ for J = 8–30 |

This is consistent with the authors' 0.0026 on 39 states; it is slightly larger, because my sample includes J = 16 and 20 at the edge.

### 2.3 The authors' J ≤ 26 statement (coordinator item 2a): independent check [DV, sampled]
`fp_ratio.py`, `cache/fp_ratio_1000_J26.txt`.
* Claim checked: max over u of |H_ang/H_hyb − 1| ≤ 0.01 for J ≤ 26.
* Sampling: 24 masses per channel, uniform near the twist edge plus geometric to u = 10⁻⁴·u_max.

**Result: confirmed on every sample.**
* All 37 physical channels with J ≤ 26 satisfy max_u |H_ang/H_hyb − 1| ≤ 0.01.
* The worst case is **4.71·10⁻³** (J = 17, `+−`, at the twist edge u = 0.967). Every `+−` channel with 5 ≤ J ≤ 26 peaks at the edge with
  1–4.7·10⁻³. The `++` channels are ≤ 6·10⁻⁵.
* min_u H_ang/N_s = **0.20087** (J = 5, `+−`), consistent with the authors' H_ang ≥ 0.1979 N_s.

This is a sampled check, not a reproduction of their Bernstein certificate. I did not build the u-polynomial hybrid (exact Weyl shift
plus conical-function enclosure of cv). My own certified statement for these spins is the stronger-in-sign but weaker-in-margin
H_ang,b > 0 of §3.

## 3. Every-mass certificates at fixed J (task 3) [certified]

`fp_poly.py`. 𝒫_{J,b}(w) is assembled in Arb from the moment balls and the exact Q_k. Positivity on [0, w_max], w_max = λ/((J+Δ_g−2)(J+Δ_g−1)),
is proved from the Bernstein coefficients (all > 0), with de Casteljau bisection where needed. Channels:
* b = 0: even J (`++` and the even-J parity-odd structure);
* b = 4: J ≥ 4 (`+−`; the single parity-odd structure for odd J).

Since H_ang,b = (1−r)³·𝒫·Tⁿ/∏_{j<n}(T−u_j), and the factor is > 0 on the heavy range, 𝒫 > 0 means **H_ang,b > 0 for every mass m ≥ m_gap(J)**.

| Δ_g | moments | channels certified | first failure (cause) |
|---|---|---|---|
| 1000 | K = 300, 3000 bits | **every physical channel with J ≤ 239** (356 polynomials, ≤ 15 cells each, ≈ 2 s each at J ≈ 240) | J = 240: balls too wide, from moment precision, not a sign |
| 1000 | K = 1000, 7000 bits | **POLY1000_RESULT** | |
| 300 | K = 300, 3000 bits | **every physical channel with J ≤ 285** (426 polynomials) | J = 286, b = 0: precision |

Further details:
* **Margins.** The minimum Bernstein lower bound is smallest at the twist edge. It decreases with J: 𝒫(w_max) ≈ 0.52 (J = 4, `++`),
  2.2·10⁻³ (J = 40), 2·10⁻⁶ (J = 238). That decrease comes from the normalising factor ∏(1 − u_j w/λ); H_ang itself is ≈ H_hyb > 0.
* **Cancellation.** The terms of 𝒫 cancel by a factor ≈ e^{1.3J} (§1.4). So certifying spin J needs moments with relative radius
  ≲ 10^{−0.57J}. This, not the Bernstein step, is the cost driver. The cost is polynomial: moments to precision e^{−1.3J} cost Õ(J²) at
  O(J) digits.

## 4. Asymptotic regions (task 4)

### 4.1 Coefficient-envelope (tail) certificate of §§29, 31, 32: it cannot apply at any accessible K [certified numbers + D]
`fp_tail.py`, `cache/fp_tail_1000_K*.txt`, with P_k = D_k²L_k and P*_k = ((5/2)_k/k!)(A + B/log(k+2)).

| k | 0 | 10 | 100 | 200 | 300 | TAILK1 |
|---|---|---|---|---|---|---|
| log₁₀\|P_k\| | 1.32 | 44.28 | 275.80 | 438.96 | 556.55 | TAILK2 |
| log₁₀ P*_k | −3.76 | −2.86 | −1.70 | −1.31 | −1.08 | TAILK3 |
| log₁₀\|D_k²Q_k\| | 1.11 | 44.13 | 275.48 | 437.64 | 550.96 | TAILK4 |

* **Envelopes.** |P_k − P*_k| ≤ P*_k/4 holds for **no** k ≤ TAILKMAX. Neither does |D_k²Q_k| ≤ (k+1)P*_k.
* **Prefix constants** (k < 300, certified): log₁₀D₀ = 555.6, log₁₀D₁ = 562.6, log₁₀D₂ = 560.6, log₁₀D₃ = 550.1. With the K = 1000 data:
  TAILD.
* **Hypothetical cutoff.** If the envelopes held from K = 300 on, §32's cutoff would be w_* = 10^{−229}, i.e. β_* ≈ 528 and
  1 − r ≤ 10^{−229}. They do not hold.
* **Why.** §1.4: P_k ~ e^{O(Λ)} up to k ≈ Λ/2. After that comes an oscillating plateau of size ~ e^{(π/2)Λ}k^{3/2} × (Fourier decay
  k^{−1/2}).
  * At Δ_g = 300 the plateau is directly visible in the certified moments: |Z_k| ≈ 10^{198} for k = 200, 250, 300, with sign changes
    and Y/Z changing sign. Compare e^{(π/2)·298.5} = 10^{203.6}.
  * So the envelope can hold only once the k^{−1/2} Fourier decay has eaten e^{(π/2)Λ}/(A + B/log k), i.e. **K ≳ e^{πΛ}·10⁸
    ≈ 10^{1370}** at Δ_g = 1000 [D, heuristic]. The prefix constants up to such K are ≥ 10^{670} (certified lower bound: the k ≈ 500
    peak).
* **Verdict on the tail certificate.** As formulated (|P_k − P*_k| ≤ P*_k/4 for k ≥ K, plus absolute prefix sums), it cannot be made to
  work for the repaired weight. Asymptotically it is true, but its threshold is astronomically far away, and the prefix sums are
  ≥ 10^{670}, which forces w_* ≤ 10^{−270} even in principle. The Cauchy-majorant inputs of the forward strip do hold for k ≤ 300:
  max_k |P_k|/Λ^{2k} = 10^{1.32} ≤ 200 and max_k |D_k²Q_k|/Λ^{2k} = 10^{1.11} ≤ 120 [certified]. I did not reproduce the strip
  argument itself (b̂ ≤ 1, J ≥ 27, error < 0.06 N_s).

### 4.2 Bounded impact parameter, comparison constant 𝒦_b(r₀) [certified partial sums; a lower bound on 𝒦]
`fp_regions.py`. Define
𝒦_b(r₀) ≥ sup_{r≤r₀}(1−r)³Σ_{k≤300}D_k²[r^{k−1}S_{b,k}|Z_k||1/(k+3) − r/(k+1)| + r^k(|2Y_k + V_{b,k}Z_k| + Λ²|Q_k|)].

The truncated sum is already a lower bound. For r₀ ≤ 10⁻³ it is also essentially the full value: the terms peak at k ≈ Λ√r₀ ≤ 32.

| r₀ | b̂₀ | 𝒦₀(r₀) ≥ | H_LR(r₀) | covered only for T ≥ 𝒦/H_LR (b = 0; b = 4 similar) |
|---|---|---|---|---|
| 10⁻⁶ | 2.0 | 1.8·10⁷ | 17.9 | 9.9·10⁵ (≈ the gap: Δ ≳ 1000) |
| 10⁻⁵ | 6.3 | 8.7·10⁷ | 4.50 | 1.9·10⁷ (Δ ≳ 4400) |
| 10⁻⁴ | 20.0 | 7.0·10¹¹ | 0.120 | 5.8·10¹² (Δ ≳ 2.4·10⁶) |
| 10⁻³ | 63.2 | 1.8·10²⁷ | 4.1·10⁻³ | 4.4·10²⁹ |
| 10⁻² | 201 | 3.9·10⁸⁴ | 1.4·10⁻⁴ | 2.7·10⁸⁸ |
| 10⁻¹ | 666 | 3.5·10²⁶⁵ | 1.9·10⁻⁵ | 1.9·10²⁷⁰ |

So the absolute-value majorant 𝒦_b(r₀)/T is useful only for b̂ ≲ 20, i.e. for spins J ≲ 10⁻²√T. There it covers states that are
already inside T0's Part F/T region b̂ ≤ 16. **𝒦 grows like e^{2b̂}**, and region (a) contributes nothing beyond the fixed-J certificates.
§29's signed remainder is the same sum before taking absolute values; its terms carry the same e^{2b̂}, so it needs the same precision
as a direct certificate.

### 4.3 Large impact parameter, r → 1: no explicit threshold obtainable from coefficients [D]
* The limit 𝔑(γ) ≥ 1 − 105/(8γ²) presupposes Z_k ≈ z_*k^{−5/2}/log k for k ~ x/(1−r).
* By §§1.4 and 4.1 that holds only for k ≳ 10^{1370}. A coefficient-based threshold would therefore be R₀ = 1/(1−r₀) ≳ 10^{1370},
  i.e. β₀ ≳ 3200. Even then, the e^{(π/2)Λ}-size oscillatory part of all k below it must be controlled with its cancellation.
* That cancellation is real. It is the Mellin-type suppression of Σ_kU_k k^{iΛ}… by the smoothness of U_{b,k} in k. But no absolute
  envelope sees it.
* I did not produce a certified threshold. The Step-C identity and the 𝔑 bound are fine as calculus; the input asymptotics are what fails.

### 4.4 What does work for large J: the frequency kernel [DV]
`fp_kernel.py`. Swap the k-sum and the ν-integral: H_ang^{(2)} = ∫g(ν)𝒦_fin(ν; J, T, b)dν, where
𝒦_fin = (1−r)³Σ_k D_k²U_{b,k}[z_{k+2}(ν)(…) + (2y_{k+2}(ν) + V z_{k+2}(ν))/T]. This finite sum of polynomials is the exact finite-mass
analogue of the shifted Legendre kernel 𝒦_LR (identity I2).

Measured at J = 40 … 2000, Δ = J + 1000 (b = 0), ν = 0.5 … 1000:
* the individual terms reach 10⁵¹⁴, yet 𝒦_fin is of the same polynomial size as 𝒦_LR;
* (𝒦_fin − 𝒦_LR)·T/(1+ν²) = O(10) for ν ≤ 5 and O(10⁴–10⁵) at ν = 50;
* the deviation becomes O(1) relative only at ν ≳ 300, where ν² ~ T/10.

So **the exponential sizes live only in the moment basis**: per frequency the finite-mass kernel is tame. A large-J proof should bound
𝒦_fin − 𝒦_LR (or 𝒦_fin minus the hybrid's flat-plus-curvature kernel) per frequency, with an explicit O(ν²/T) or O(ν²/m²) error.
That error must then be integrated against |g|, exploiting (1−p)⁵ near ν = Λ. Do not try to bound moments. A closed form for 𝒦_fin
(the double sum is hypergeometric; U_{0,k} = (−J)_k(J+1)_k/((2−Δ)_k(Δ−1)_k)) would make this analytic. **This is the recommended route
for the J > J_max region.**

### 4.5 The coordinator's pivot (§29 signed remainder in r): assessment [certified numbers]
`fp_pivot.py`, `cache/fp_pivot_1000.txt`. F_j(r) = Σ_{k>j}P_k r^{k−j−1} is not a tame "fixed-r leading-Regge integral". It is the tail of the
coefficient series, so per frequency it is the tail of the CDH generating series, of size e^{(π/2)ν} until j ≳ ν/2. Hence
|F_j| ~ |P_{j+1}| ~ e^{O(Λ)}. The bracket M_j − (wu_j/(1−u_j/T))F_j is a difference of two numbers of that size, and the sum over j cancels
again.

| state | H_ang^{(2)} | truncated H_LR^{(2)} | Σ_j U_j(\|M_j\| + …\|F_j\|)·w³/T ≥ |
|---|---|---|---|
| J = 30, Δ = 1030 | 4.9485·10⁻³ | 4.9327·10⁻³ | 6.3·10¹⁶ |
| J = 100, Δ = 1100 | 1.81498·10⁻⁴ | 1.81472·10⁻⁴ | 2.6·10⁶² |
| J = 100, Δ = 3000 | 3.42607·10⁻³ | 3.42393·10⁻³ | 2.0·10²⁰ |
| J = 200, Δ = 1200 | 4.19082·10⁻⁵ | −8.5·10⁸⁹ (series not summable at K = 300) | 1.9·10¹²¹ |

* **Result.** The pivot is an exact rewriting with the same cancellation; it does not close the middle region.
* **The obstacle.** The large spins cannot be reached by any argument that treats U_{b,j} ≤ r^j and |P_j| separately. The information
  sits in the joint smoothness of U_{b,j} in j and the oscillation of P_j (a Mellin-type cancellation).
* **Cost.** Making the pivot rigorous at fixed (J, T) is the same as evaluating H_ang exactly, which §3 already does for all masses at once.
  Uniformly in J it would need a new analytic lemma, essentially §4.4's kernel bound.

### 4.6 The remaining region
* **Certified:** every physical state with J ≤ J_max (= 239 at Δ_g = 1000; 285 at Δ_g = 300), at every mass: H_ang,b > 0.
* **Not certified:** every J > J_max, at every mass. Region (a) adds only b̂ ≲ 20 at T ≥ 10⁷–10¹³, which is not a new spin range. Region (b)
  has no usable threshold. The authors' forward strip (b̂ ≤ 1, J ≥ 27) is not reproduced.
* **Cost of extending the direct certificate.** Spin J needs moments with relative radius ≲ 10^{−0.57J} for k ≤ J. In practice:
  * ≈ 2J + 500 decimal digits of working precision;
  * edge panels widened so the θ-bound is below that, e^{−(c_cut)²/101.6};
  * O(J²) ball operations for the moments and O(J²)–O(J³) per polynomial.
  * J ~ 10³ is a few CPU-hours (§3). J ~ 10⁴ would need ~2·10⁴ digits and is impractical. No finite computation reaches J = ∞, so the
    large-J region needs §4.4's kernel route (or T0's hybrid plus a per-frequency comparison).

## 5. Reproduce (`ads4/code/fullpos/`; `ulimit -v 3000000; export OMP_NUM_THREADS=1`; ≤ 2 processes)

```
python3 fp_norm.py                                  # §2.1 identities I2, I3 (mpmath, ~2 min)
python3 fp_moments.py 1000 300 3000                 # §1 certified moments, Dg = 1000, k <= 300 (~75 s) -> cache/moments_1000_K300.json
python3 fp_moments.py 300 300 3000                  # Dg = 300 (~70 s)
python3 fp_poly.py 1000 300 0-239                   # §3 certificates (~10 min) ; fp_poly.py 300 300 0-285
python3 fp_compare.py 1000 300                      # §2.1-2.2 vs T0 hybrid (largespin read-only, ~10 min)
python3 fp_ratio.py 1000 300 26 24                  # §2.3 J <= 26 check (~15 min)
python3 fp_regions.py 1000 300; python3 fp_tail.py 1000 300 300; python3 fp_pivot.py 1000 300; python3 fp_kernel.py   # §4
FP_NGL=1200 FP_CUT=300 FP_QPREC=3300 python3 fp_moments.py 1000 1000 7000      # higher-J moments (edge widened; see §1.2)
```

| file | content | grade |
|---|---|---|
| `fp_moments.py` | Y_k, Z_k in closed form with θ edge quadrature and bounds | [certified] |
| `fp_poly.py` | 𝒫_{J,b} > 0 on [0, w_max] (Bernstein + bisection); Q_k exact | [certified] given the moment balls |
| `fp_norm.py` | identities I2, I3 | [DV] 40 digits |
| `fp_compare.py`, `fp_ratio.py` | H_ang vs T0's hybrid | [DV] |
| `fp_regions.py`, `fp_tail.py`, `fp_pivot.py` | coefficient profile, 𝒦(r₀), envelopes, prefix constants, majorants, pivot | [certified] numbers |
| `fp_kernel.py` | finite-mass kernel (DV), §31 identity (exact) | [DV] / exact |

## 6. Caveats
1. **Rigour of the moments.**
   * Every ingredient is closed-form Arb, except the θ edge correction. That uses Gauss–Legendre with the Bernstein-ellipse remainder
     (Trefethen, Thm 19.3); the sup of |h| on the ellipse is enclosed by 256 Arb boxes. The box cover's radius is the arc length π/NB
     times |dz/dφ|, a deliberate overestimate.
   * The "rest" θ bound uses 1 − θ ≤ erfc((ν_c−ν)/10) ≤ e^{−(ν_c−ν)²/100} and |R| ≤ 84 on ℝ. I checked R ≤ 147456/(4·9·49) by hand.
   * Not independently reviewed.
2. **Physical identification.** "H_ang,b > 0" is positivity of the authors' finite-mass angular action (ANGULAR_ROUTE §§26–28). Its
   identification with Ω^EE on heavy states uses the retained physical identities (ladder, first-moment operator, flip cancellation), which
   are graded elsewhere (LADDER_ALLJ: [P] modulo two operator identities). The normalisation ℬ > 0 is the authors' (§28). For positivity
   only its sign matters.
3. **Hybrid side.** The hybrid comparisons use T0's float `hyb_fast` with the target ψ, not the repaired realised weight; the difference is
   ≤ 5.4·10⁻⁹. They are sampled, so [DV].
4. **K = 1000 moments.** K1000_NOTE
