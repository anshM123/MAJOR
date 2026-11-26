# REVIEW_AUTHOR_CHAIN: adversarial referee report on the §104–112 chain (sliver positivity, L ∈ [2.6, 3.65]; low-L prefix)

Referee agent, 27 Nov. Scope: the mathematics the coordinator did not re-derive in §105/§107/§109–111, and the coverage logic.
No existing file was edited. Scripts and logs are in `code/review_chain/`.

## 0. Verdict

| item | verdict |
|---|---|
| 1. `nb_delta` seed, downward steps, Δ-recurrence | **holds** (derived; hypotheses hold on every cell; numerics exact/bound ≤ 0.48) |
| 2. `B_box`, `inf_BH`, `pref`, `psidiv0`, contour normalisation, infinite-twist reference | **holds** |
| 3. Base-multiplier lemma (`base_bounds.py`), including negative Re α by downward steps and the b = 4 helicity reduction | **holds** |
| 4. `low_prefix.py` truncation (2^-110) and pole exclusion | **holds** (truncation only; the sign of the rational prefix is NOT certified) |
| 5. Coverage / logic (cells, d-slabs, ν-bands, channels, contour pieces, floors, Loc-Joint J₀, §73) | **holds**, with two documentation gaps (§5.3) |
| 6. Low-frequency theorem (constant 14) | **holds** |

I found **no error**. L ∈ [2.6, 3.65] positivity stands, conditional on the retained inputs listed in §7. The low-L truncation theorem stands as a truncation/pole statement only; L ∈ [0.051, 0.82] is not closed.

## 1. `nb_delta` (contour_direct.py) and the near-real recurrences

**Derivation.** X(α) = d^α Σ_{k≤n} (α)_k U_k/k!, U_k = c^k E_k, E_k = Π_{j<k} ρ_j, ρ_j = (1 − hp_j/c)/(1 − hp_j), p_j = j(j+s). At h = 0, X ≡ 1. For j < n, 0 < ρ_j ≤ 1 and 1 − ρ_j = hp_j d/(c(1 − hp_j)). ρ_n = 0, so E_k = 0 for k > n. Telescoping gives |1 − E_k| ≤ Σ_{j<k}(1 − ρ_j) for k ≤ n, and |1 − E_k| ≤ 1 always.
* k ≤ n/2: hp_j ≤ c(n+2s)/(4(n+s)) ≤ 0.2522 (n ≥ 997, s ≤ 9), so (1 − hp_j)⁻¹ ≤ 1.34 ≤ 2.
* k > n/2: f(k) := 2(hd/c)Σ_{j<k}p_j ≥ d(n²−1)/(12(n+s)), using h/c = 1/(n(n+s)). With x = n/a one has dn = n(2x+1+s/a)/((x+1)(x+1+s/a)) ≥ ~747 for n ≥ 997, a ≥ 998, so f(k) ≥ 60 > 1 ≥ |1 − E_k|. The comment "dn ≥ 497" is a valid but loose lower bound.

So |X − 1| ≤ |d^α| Σ |(α)_k|c^k f(k)/k!. Next, |(σ+it)_k| ≤ (σ)_k exp(t²(σ⁻² + σ⁻¹)/2). The negative-binomial factorial moments are E K^{(m)} = (σ)_m (c/d)^m, and Σ_{j<k}p_j = K^{(3)}/3 + (s+1)K^{(2)}/2. Together these give exactly

  |X − 1| ≤ C·2(h/d²)[σ(σ+1)(σ+2)c²/3 + (s+1)σ(σ+1)cd/2],

which is the code with c ≤ 1, cd ≤ d_hi, g = hcap(d_hi)/d_hi² ≥ h/d² (h/d² ≤ 1/(998²(1+√(1−d))²) is increasing in d), C at σ_lo and the polynomial at σ_hi. The lemma is valid for complex α with Re α = σ > 0. n ≥ 997 holds because n = J − b ≥ 997.

**Downward steps.** Lemma V was checked to 1e-80 at 8 (state, β) pairs, among them n = 997, a = 998, s = 9 at Re β = −7.75 (check_nb_lemmaV.log). One large-|β| row lost precision at 80 digits; that is an artefact of the check, and the identity is polynomial in β. step_coeff's A, C, den and force each majorise the corresponding modulus, with B ≥ |β| including the imaginary part and c ≤ 1. The weighted j-loop uses |r_j| ≤ ((j+½+y)² + t²)/((j+1)(j+1+2y)), with the tail ratio (1 + (y+½+t)/(J+1))². The p = 3 dressing term per j is (1 + τA₀)Δ_j + τA₀ + τA₁(1 + Δ_{j−1})/d with τ ≤ h. Its constants DA ≥ |A₀| and DB ≥ |A₁| match TAIL_CERT E2.

**Numerics.**
* check_nb_lemmaV.py: exact |X − 1| versus nb_delta at L ∈ {2.6, 2.9, 3.65}, a ∈ {998, 1500, 4000}, b ∈ {0, 4}, Re α ∈ {¼, ¾, 5/2, 7/2}, Im α ≤ ¼. Worst exact/bound is 0.476.
* check_near_sum.py: exact Σ_{j<1500}|c_j|d^j|dressed X_j − 1| plus the infinite-twist tail, from exact Arb ladders and the exact β-recurrence. It was run at physical states on the [2.6, 2.9] slab edges, on the bottom (ν = 0.1) and the left edge up to y = 20.25 (p = 4, Re α = −16.75). Worst exact/near_sum is 0.674 (p = 4, b = 0, y = 20.25, a = 998).

## 2. Gamma/Euler envelopes and contour normalisation

* **B_box.** Duplication gives B = 2^{2iν−1}Γ(iν)/(√π Γ(iν+½)). The shift by k > y uses factors |z+j+½|/|z+j|; (t²+A²)/(t²+δ²) is monotone in t², so endpoints suffice. Then |Γ(w)/Γ(w+½)| ≤ Γ(σ)/Γ(σ+½), from the beta integral with Re w = σ, and Gautschi gives ≤ √(σ+1)/σ. The code uses √(σ_hi+1)/σ_lo, which is looser and valid. There is exactly one 1/√π, as the duplication formula requires.
* **inf_BH.** Euler's integral with b₀ − a₀ = a₀. The prefactor |Γ(2iν)Γ(b₀)/(Γ(a₀)²Γ(½+iν)²)| = |coth πν|/(2π) was checked numerically to 12 digits at complex ν. |coth π(t+iy)| ≤ coth(πt) ≤ 1.53 < 2 for t ≥ ¼. ∫[u(1−u)]^{y−½} ≤ π for y ≥ 0 and ≤ 4^{½−y} for y ≥ ½, and (1 − du)^{−y−½} is monotone.
* **check_gamma_env.py** compared 300 complex points (t ∈ [¼, 998.5], y ∈ [0, 21], d ∈ {.026, .055, .0743}) plus 6 interval boxes with mpmath. Worst |B|/env is 0.78, worst |BH|/env 0.69 and worst box ratio 0.67.
* **pref.** d^{3 − Re α} is taken at d_hi when the exponent is ≥ 0 and at d_lo otherwise, and the exponent is monotone in y, so it is taken at y_lo. The factor 4/Λ for p = 3 comes from (2/Λ)·2Re. The factor 2τΛ for p = 4 comes from G_D = Λ∫ψ₃(ν/Λ)F₄dν with τ ≤ hcap. F_{p,ν} = (1−z)^{−p}P_{−½+iν} (fp3_common), so the outgoing exponent for p = 4 is 7/2 + iν with the same B and c_j. psidiv0 gives |ψ(ν/Λ)| ≤ ν·psidiv0, and |νB|² = ν coth(πν)/(4π) ≤ δ/(4π) + 1/(4π²) < .22².
* **Identity versus ANALYTIC_SLIVER §4.** The pieces are bottom (2/Λ)∫₀^δψ(K−K_∞), left, top, right and the extra segment [Λ, 998.5]. That segment is needed because ∫_δ^Λ = ∫_δ^{998.5} − ∫_Λ^{998.5} with the polynomial ψ continued past Λ. κ_J (J₀ = 1500) is analytic in the closed rectangle: B has poles only at ν = ik/2, the c_j poles lie in Im ν < 0, and the D_j poles are excluded by the majorant. So the normalisation matches.
* **Reference.** The main term of FULLPOS3 §4.1(c) is G_A^ψ(c) + τσ_bG_D^ψ(c). The p = 3 contour compares the dressed κ_J with the undressed κ_∞, so τ𝒟_bG_A(c) is absorbed into the error, which is correct. The p = 4 contour compares with τσ_bG_D^ψ(c) through the same B·H.
* **Transfer certificates in situ (check_transfer_points.py).** Exact P and Q (signed 1500-term recurrence, dressed for p = 3) were computed at physical states (a = 998, 1500, 40000) and compared with the bounds the [2.6, 2.9] assembly actually uses (Extend.get: real norming plus the Hadamard ellipse, or the absolute majorant). The test points are on the top line (p = 3, Y = 10.25; p = 4, Y = 20.25), the right edge, and real ν = 300 and 1000. Worst ratio is 0.96, at ν = 5 + 10.25i, where the absolute majorant is nearly sharp.

## 3. Base-multiplier lemma

I re-derived every step. The steps are:
1. The eigenframe identity ‖T⁻¹T′‖_∞ = (|λ₊′| + |λ₋′|)/(2|r|) = √(1+κ²)(|w+i|⁻¹ + |w−i|⁻¹)/2.
2. Its antiderivative A(Y, κ); the asinh(1/κ) terms cancel.
3. Branch continuity: 1 + w² stays in the upper half-plane for κ > 0.
4. |λ₋| ≤ 1 from u ≥ 0.
5. The final row norm 2cosh u = |w+i| + |w−i|, nondecreasing in y and in κ.
6. The normalisation ∫₀ⁿη = −log d and the bound η′ ≥ η², so dx ≤ ν dy/y².
7. u ≤ κy/√(1−y²) for y < 1, and u ≤ acosh y + κy/√(y²−1) for y > 1. The second uses (t² − A)² + 4t²y² ≥ A² with A = y² − 1.
8. B(Y), including acosh y ≤ √(2(y−1)).
9. c₀ ≤ 1.01 from |w₀| ≤ |α|/(a(a+1)) < .002.
10. Localisation: t′(2−t′) = d(1 − (2N+1)⁻²), y(n) = (ν/a)·2R′/(1+R′), and R′ ≤ √(1 − d_min(1 − 1/(4a²))).
11. The interval version (q_max/q_min)·A(Y, κ_min) ≤ q_max·A.

**Negative Re α.** Downward steps from exponents ≥ ¼ are valid. The denominator bound |1 − hβ(β−1)| ≥ Re(·) = 1 − h(σ² − σ) + hτ² ≥ 1 − h(|σ|+1)² is correct even though it ignores τ, because τ² enters with a plus sign.

**Helicity.** d/d₀ = 1 + 20h for s = 9. U₄^{(1)} ≥ (1 − d_hi − 12h)⁴, with h′ = 1/((N+4)(N+5)) ≤ h. The shifted d₀ lies in [d_lo/(1+20 hcap), d_hi], and a is unchanged.

**Numerics (check_base.py, Arb at 9000 bits with a radius assertion).** Exact |X(α)| and |X(α+1)| (b = 0 via the s = 1 ladder, b = 4 via the s = 9 ladder) were compared with base_bound on the [2.6, 2.9] slab. The tested set is a ∈ {998, 2500} at both slab ends, Re ν ∈ {¼, 1, 160, 600, 950, 998.5}, and Im ν = 0 and the contour heights 10.25 (p = 3) and 20.25 (p = 4; Re α = −16.75). Results are in check_base_26_p3/p4.log, which list worst ratios. All ratios are < 1, with a maximum of 0.992 at real ν = ¼, where both sides are ≈ 1.

A first run at 4000 bits produced a spurious "violation" (1.5e27 at ν = 998.5, a = 998, d = .055). Arb's radius showed it to be precision loss; the true value is 8.78. This is noted so nobody re-reports it.

## 4. low_prefix.py

All constants were derived and found correct:
* h ≤ d²/(a²(1+√c)²) = (1−√c)²/a².
* The Mehler bound on |z| ≤ .995: |1−z|^{−p} ≤ 200^p and |Im ξ| < π/2. On the segment z ∈ (−1, 0), where Re ξ = 0, Lemma B holds by continuity, with P_{−½}(cos θ) ≤ 1.18.
* Dressing ≤ 1 + (2/.005² + 5/.005 + 8)/998² < 1.082.
* 4·200⁴·120000001 < 2⁶⁰ and e^{1100π/2} < 2^{2494}, using π < 22/7.
* Cauchy |g_k| ≤ M·.99^{−k}, and the tail Σ_{k>K}|g_k|c^k ≤ M(c/.99)^{K+1}/(1 − c/.99) uses only 0 ≤ U_k ≤ c^k, which includes U_k = 0 beyond termination.
* The ν > 1100 part is 12·200⁴·gaussian_mass(1100, 11/7), where 2/Λ < 1 is dropped and ψ₃ has no content above Λ.
* The pole gap uses h_cap(c_lo) and s = 9 for both b. The prefix then reproduces the zero factor j = n for n < K, and j < K < N.

Status: truncation and poles only. **The rational numerator's sign is uncertified, so L ∈ [0.051, 0.82] is not proved.**

## 5. Coverage and logic

### 5.1 Pieces and channels

* **Channels.** Every cell charges p ∈ {3, 4} × b ∈ {0, 4}. ψ₃ is ideal, so p = 4 has no real-axis residual.
* **Frequency.** The pieces are:
  * [0, 160]: the low-frequency theorem, valid for L ∈ [0.8, 3.65];
  * [159, 1040.1]: 16 bands covered by the real boxes 159–400.1–700.1–850.1–950.1–1040.1 (choose_row's 1e-9 slack is harmless, because Extend then uses an ellipse with η > 1);
  * [1040, 1150]: the absolute majorant with Gaussian mass ∫_l^∞;
  * above 1150: the unsplit kernel ≤ 200e^{πν/2}. I checked this constant: 1.6·(π/2)·4³·(1+o(1)) ≈ 161.
  * Loc-Joint (J₀ = 1500, ν ≤ 1150) matches the J = 1500 used by major, Extend.outer (K = 1500) and near_sum.
* **Contour.** The pieces are bottom [0, ¼], left (Re ν = ¼, 0 → Y), top (¼ → 998.5 at Y), right (998.5, 0 → Y) and the extra real segment [998.4998, 998.5] ⊃ [Λ, 998.5]. The top uses t-edges 1, 10, 159, 400, 700, 850, 950, with no gap.

### 5.2 d-slabs, floors and high twist

* **d-slabs.** The slabs tile L ∈ [2.6, 3.65]. Relative margins at the top end are 2.9e-7 ([2.6, 2.9]) to 8.7e-5.
* **Floors.** Recomputed from lrt3_rcert.json with β = 2 acosh e^{L/2}. The β-ranges are tiled contiguously; Fcert is stored as an Arb lower endpoint (`lo`).

  | L-cell | rows | min F_cert | floor used |
  |---|---|---|---|
  | [2.6, 2.9] | 14 | 2.9936e-6 | 2.5e-6 |
  | [2.9, 3.2] | 13 | 2.5465e-6 | 2.5e-6 |
  | [3.2, 3.35] | 7 | 2.3458e-6 | 1.87e-6 |
  | [3.35, 3.5] | 7 | 2.1595e-6 | 1.87e-6 |
  | [3.5, 3.65] | 7 | 1.98662e-6 | 1.98e-6 |

* **High twist.** The new cells below 3.5 are uniform in a ≥ 998 (full twist), so §73 is needed only on [3.5, 3.65]. There β ∈ [4.871, 5.023] ⊂ [3.04, 8.24], where a* = 4500. ✓

### 5.3 Documentation gaps (not errors)

1. **L convention.** TAIL_CERT defines w = 1 − C/T = e^{−L}, while the slabs and floors use d = e^{−L}. For b = 4, d = w(1 + 20h). Every slab's top margin (≥ 2.9e-7 relative) exceeds 20·hcap (≤ 2.9e-8), and the floor β-ranges have slack ≥ 0.05. So both conventions are covered. The shipped checkers (assemble_full.py, assemble_26.py) only test e^{−L}; the convention should be stated once.
2. **[2.6, 2.9] floor.** assemble_26.py checks only the β-span, not the F_cert rows. The coordinator's §111 check and mine both confirm 2.9936e-6 ≥ 2.5e-6, but no shipped script asserts it.

## 6. Low-frequency theorem

I re-derived the following (constants as in check_lowfreq.verify_constants):
* 1 − r = t + v(1−t) ≥ 1.555t;
* tan|Im β| = 2√r|sin(θ/2)|/(1−r) ≤ X/1.555;
* the Cauchy phase increment .0065;
* W = t^{−a−1}(1+X²)^{−(a+1)/2};
* the Jacobian (1 − X²/X_max²)^{−½} ≤ 1.004;
* the Gaussian completion 4ν²/(7k²(a+1)) ≤ ν²/4200, and the half-line normaliser, which give ratio ≈ 2.149 ≤ 2.16;
* the lower-tail Chernoff a(−log .752 − .248) ≈ 36.9, using I ≥ −log x − (1−x);
* the upper tail I(.662; .26)·998 ≈ 1293;
* dressing ≤ 1.081;
* d ≤ 2t;
* tail monotonicity for X ≥ ½ (ν·k/(k²+X²) ≤ 160 < 399);
* K₀: a saddle computation gives K₀ → √(a/(a+1)), consistent with < 1.6.

Principal term = 1.6·1.081·(11/7)·(2/1.55)³·2.16 = 12.61 < 14.

The tail prefactor (π√a/t)e^{(a+1)/(2a)} ≈ 166/t differs from my 2π√(2(a+1))/(.999√π t) ≈ 158.6/t, but it is still an upper bound, and the tail is ≤ e^{−47} in any case. **Holds**, conditional on the retained masses and floor.

## 7. Retained (not reviewed here)

* The angular-action identity.
* EPSW residual masses, cutoff and Gaussian inventory.
* LRINF F_cert.
* Loc-Joint internals: the analytic bridge, the 208 constant and the K₀ derivation.
* The Beta–Cauchy identity.
* contour3 on [3.5, 3.65]. Its d²D ≥ 3.9a² and N ≥ 40a relations do hold there: d²D/a² ≥ (1+√c)² ≥ 3.939 at d ≤ .0302.
* §73. Its per-cell records are still only on the author's laptop.
* The helicity identity (author-checked numerically).

## 8. Files

All in code/review_chain/:

| script | log | what it checks |
|---|---|---|
| check_nb_lemmaV.py | check_nb_lemmaV.log | Lemma V identity; nb_delta |
| check_gamma_env.py | check_gamma_env.log | B, BH envelopes |
| check_near_sum.py | check_near_sum.log | bottom/left majorant vs exact |
| check_base.py | check_base_26_p3.log, check_base_26_p4.log | M0, M1 vs exact Arb |
| check_transfer_points.py | check_transfer_points.log | P, Q certificates vs exact |
