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

**PLACEHOLDER**

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

**RATIO_RESULT**

This is a sampled check, not a reproduction of their Bernstein certificate. I did not build the u-polynomial hybrid (exact Weyl shift
plus conical-function enclosure of cv). My own certified statement for these spins is the stronger-in-sign but weaker-in-margin
H_ang,b > 0 of §3.
