# The EE step of LRC\*, f_pos(±3i/2) for Q's weights, the a = b arc, and the g₅ deciding test

*23 Nov 2026. Unreviewed, not committed. Answers LRC_CLOSE §6.5 items 1–3 and runs F_NORM §2.6's deciding test. No existing file was
edited.*
* *New code:* `ads4/code/eestep/`: `ee_blocks.py`, `ee_quad.py` (modules); `ee1_blocks_check.py`, `ee2_identify.py`, `ee3_deshift.py`,
  `ee4_fpos.py`, `ee5_sd_run.py`, `ee5_sd_report.py`, `ee6_lrprofile.py`.
* *Outputs:* `ads4/code/eestep/cache/` (git-ignored by its own `.gitignore`).
* *Read-only imports and reuse:*
  * `lrc/lrc_blocks.py` (b = 0 comparison only);
  * the frame, scalar integrand and EE dressing are re-implemented from `lrc/lrc1_identify.py`, `lrc/lrc3_ee.py` and the header of
    `repro_engine/integrand.py`;
  * `repro_engine/direct_rho.py` and `sdchain.py` (§5; `direct_rho`'s save directory is redirected to `eestep/cache`);
  * data: `resultq/functionals_q7_xprec.json`, `resultq/functionals_q4_deshift.json`, `resultq/cache/stage6_gauge.pkl` and the cached
    `repro_engine/cache/direct_D{300,1000,20000}_n8-8-18-18-6-20-20_f64.npz`.
* *Compute:* about 25 CPU-minutes, at most 2 processes, one thread each, < 1 GB.

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]**; **[M]**; **[Q]** quoted;
**[A]** assumption; **[I]** open. A chain carries the grade of its weakest link.

---

## 0. Verdict

| item | before | after |
|---|---|---|
| 1. EE step of LRC\* ("spinning partial wave = shifted scalar ones") | [D-sketch], not checked | **[P] reduction + [D + DV] identification**: the actual EE kernel is an exact finite sum of 9 b-shifted scalar continuous-spin partial waves, with ν-independent coefficients (§1) |
| LRC\* for EE (k = 2 and k = 4) | [D-sketch] | **[D + DV]**: de-shift cancellation and strip analyticity hold for the EE kernel itself; exponents of (exact − LR) equal the scalar ones (§2) |
| LRC-a for EE at k = 4 | [D-sketch] | **[D]**: position-space LR form explicit and checked [D + DV]; the de-shift transfer to position space uses only the analytic structure of χ. The heavy image of the k = 4 EE position kernel was not recomputed (§1.4) |
| 2. f_pos(±3i/2) for Q's C4 weights | unknown | **= 0** [D given DS] + [NE]: χ is regular at ±3i/2 because R₋₃^{EE3}(−3i/2) ≠ 0. The stored weights' residue is ≤ 3·10⁻²⁷, against 9·10⁻⁷ for a scalar-type control (§3). **k = 4 criterion: J₀ < 5/2** |
| 3. a = b integrability of the borderline R⁴ arc | [I] | **Not separately load-bearing.** F_NORM's derivation does not prove it, but its input H_EFT(a) (ν-paired) already presupposes it. As a free-standing statement it stays [I] (§4) |
| 4. F_NORM's deciding test for (P) (g₅) | not run | **Run** [NE]. With an infinitely heavy reference, the ++ J = 0 derivative-piece correction is m-independent: s_d = 0.5081 / 0.5095 / 0.5096 / 0.5095 at m/Λ = 1 / 3.35 / 10 / 20. **Exact cancellation is excluded; the double-counting reading (lower slope 0.255) holds** (§5) |

**H_K: [D-sketch] → [D given DS].**
* All analytic inputs are now [D] or better: K1′, K2ν, K3 [D]; LRC\* [D + DV] (scalar and EE); LRC-a [D + DV] (scalar; EE at k = 2) and
  [D] (EE at k = 4).
* **Result D:** H_K is [D] (DS [D + NE]).
* **Result Q:** H_K is [D given DS]. DS is [NE] and is the only non-derived input left.
* **Criteria under J₀ ≤ 1:**

  | piece | criterion | margin |
  |---|---|---|
  | `Ĉ′_EE` (k = 2) | J₀ < 3/2 | ½ |
  | `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` (k = 4), Q's weights | **J₀ < 5/2** (generic weights: J₀ < 2) | **3/2** |
  | `Ĉ₅^{ds}`, `Ĉ₃` | J₀ < 3/2 | ½ |

**Letter sentence (replaces LRC_CLOSE §6.3):** see §6.3.

---

## 1. The EE kernel as a finite sum of shifted scalar partial waves

### 1.1 The object

The frame is H7 §18.1 / lrc1. The scalar per-ν integrand is `S_ν = A^{iν/2}|Nn|^{−½}m₂₄^{¾}(m₂₀m₄₀)^{−¼}c₁₃^{5/4}(c₁₀c₃₀)^{−¼}` on the
x₀-circle x₀ = r n̂_θ times the celestial circle. The EE kernel of Ĉ′_EE and Ĉ₄^{ds} carries two ν-independent factors:
* **(13) side.** The light-ray polarisations of T₁ and T₃ are x₁₀ and x₃₀ (v3_regge_kernel; repro_engine). Contracted with dual null
  polarisations: `τ₁₃ = (w₁·x₁₀)²(w₃·x₃₀)²`, i.e. r⁴ times a degree-4 trigonometric polynomial in θ.
* **(24) side.** `D = EE·V₀,₄₂⁻⁴ = [(H₂₄ − 2V₂,₀₄V₄,₀₂)/(X₂₀X₄₀)]²`. It is z₀- and ν-independent (lrc3).

So `K^EE_ν[w₁, w₃, z₂, z₄] = ∫dθ du S_ν P(θ)`, with `P = τ₁₃D`.
* Polarisations: w₁ and w₃ are taken in the rest frame of ŷ (the (13) pair); z₂ and z₄ are boosted with x̂ (the (24) pair).
* **Correction to LRC_CLOSE §4 step 1.** The (13) side is *not* literally the scalar one: τ₁₃ depends on θ.

### 1.2 Exact reduction [P + DV] (`ee_quad.py`, `ee2` (a))

1. **Everything is affine in cos θ on the circle.** On x₀ = r n̂_θ:
   * `X₂₀ = −(1−r²)(1 − s cos θ)` and `X₄₀ = −(1−r²)(1 + s cos θ)`, with `s = 2r sinh β/(1−r²)`;
   * z·x₀ is affine in (cos θ, sin θ).

   Hence P is a polynomial of degree ≤ 8 in (cos θ, sin θ) (τ₁₃: 4; numerator of D: 4), divided by (X₂₀X₄₀)⁴.
2. **The sin-odd part integrates to zero.** S_ν is even under (θ, φ) → (−θ, −φ), the reflection e₂ → −e₂ that fixes x̂ and ŷ.
3. **The even part is a Laurent polynomial.** With `t = X₂₀/X₄₀` one has `cos θ = (1−t)/((1+t)s)` and
   `(X₂₀X₄₀)⁻⁴ = (1+t)⁸/(256 t⁴(1−r²)⁸)`. So

       P_even = Σ_{b=−4}^{4} c_b(r, β; pols) t^b    exactly,      K^EE_ν = Σ_b c_b Ψ^{(b)}_ν,   Ψ^{(b)}_ν := ∫ S_ν t^b,

   with **ν-independent** c_b.
4. **Each Ψ^{(b)} is a scalar partial wave [D].** (X₂₀/X₄₀)^b is a covariant reweighting of the (24) structure. It shifts Δ₂ and Δ₄
   oppositely and leaves the x₀ weights alone, so Ψ^{(b)} is the scalar continuous-spin partial wave of the same representation
   (1, l = −½ + iν), with (a, b) = (0, ±b).

**Checks** (3 random polarisation sets × 9 points ρ ∈ {0.25, 0.55, 0.85}, β ∈ {0.6, 1.4, 2.2}):
* degree ≤ 8: the Chebyshev tail of P_even(X₂₀X₄₀)⁴ is ≤ 1.8·10⁻¹³;
* the sin-odd part is up to 100 % of the even part in size, but its integral is ≤ 6·10⁻¹⁵ of |K^EE|;
* `K^EE = Σ_b c_bΨ^{(b)}` to 1.2·10⁻⁹, through cancellations in the b-sum of up to 3.5·10⁶.

**Why this is the weight-shifting statement.** It is the CPPR/weight-shifting form, with one simplification. On the light-cone circle x₀
has a single degree of freedom, so the polarisation structures become *multiplicative*, x-dependent coefficients (no differential
operators). The (13) polarisations are traded for (24)-side b-shifts.

### 1.3 The b-shifted blocks (`ee_blocks.py`, `ee1`)

**Recursion [D].**
* Dolan–Osborn with a = 0 adds `−b(z²∂_z + z̄²∂_z̄)` to the (1/2)-Casimir. In radial variables this is
  `−(4br/den)[η(1−r²)r∂_r + (1+r²)(η²−1)∂_η]` (sympy-checked).
* Multiplied by (1−r²)den as in blocks3d, it adds odd-level terms r¹, r³, r⁵ to lrc_blocks' recursion.
* **Normalisation:** `Ĝ_b = Σ_{n≥0} rⁿ Σ_m B^{(b)}[n,m] P_{l+m}(η)` with `Ĝ_b → P_l` as r → 0, for every b.

**Checks [DV]:**

| check | result |
|---|---|
| b = 0 against lrc_blocks | identical |
| Dolan–Osborn PDE, finite differences in (z, z̄) | residual 2·10⁻¹²–5·10⁻⁸; the b = 0 operator on the same functions fails at 4·10⁻³–10¹ |
| polynomiality: (n+1)-th b-difference of B^{(b)}[n,m] | ≤ 2·10⁻¹⁵ (n-th: 10⁻²–0.2). Degree exactly n, as the recursion implies [P] |
| coefficient poles in 0 < Im ν ≤ 3 | none at Im ν = ½, 3/2, 5/2 for b = 0, 1, −2, 3. Only the integer ones from the Legendre denominators, which cancel in Ĝ as for b = 0 |
| edge structure: Ĝ_b − 2Re[c_l e^{lβ}F_b(ρ)], with F_b = (z/4ρ)^h ₂F₁(h, h+b; 2h; z) | relative to the envelope: e^{−2β} decay to 10⁻¹¹–10⁻⁸ at β = 10. One exception: b = −3 at ρ = 0.85 floors at 3·10⁻⁵ from the negative-b series (see below) |

* **Pole set of F_b.** Its poles are at Im ν = 3/2 + 2p only, for every integer b [P; lrc3 table].
* **Negative b.** The radial series loses up to 9 digits there (b = −4, ρ = 0.85). It is evaluated instead through the exact identity
  `G^{(0,−b)} = v^b G^{(0,b)}`, with v = (1−z)(1−z̄) [P: conjugation of the DO operator], checked to ≤ 10⁻¹².

### 1.4 Identification and end-to-end check (`ee2` (b)–(d), `ee6`)

**Identification (prefactor-free, as lrc1).** `[Ψ^{(b)}_ν/Ψ^{(b)}_{ν₀}]/[Ĝ_{±b}(ν)/Ĝ_{±b}(ν₀)]`:
* x-spread 2·10⁻⁷–2·10⁻⁶ for b′ = ±b (the two are equivalent by the identity above);
* 1.0–3.7 with b′ = 0.

**Prefactor.** `Ψ^{(b)}_ν = N(ν) r² v^{−b/2} Ĝ_{−b}`, with fitted exponents r^{2.000000} v^{−b/2.000000} (residual ≤ 2·10⁻⁸). Hence
`Ψ^{(b)} = Ψ^{(−b)} = N r² v^{|b|/2} Ĝ_{|b|}`, and only positive-b series are needed.

**End to end [DV].** The directly integrated EE kernel `K^EE_ν(x)` equals `N(ν) r² Σ_b c_b(x) v^{−b/2} Ĝ_{−b}(x; ν)`, where:
* N(ν) is taken from the scalar kernel at one point;
* the blocks come from the recursion.

| polarisation set | max \|K/pred − 1\| | cancellation in the b-sum | control: scalar ν-profile only |
|---|---|---|---|
| 0 | 3.1·10⁻⁷ | 8.8·10³ | 0.64 |
| 1 | 1.5·10⁻⁶ | 1.4·10⁷ | 2.7 |
| 2 | 1.6·10⁻⁷ | 8.9·10³ | 1.6 |

**Leading Regge [D + DV].**
* `L_ν(β) := lim_{r→0} K^EE/(N r²)` is extracted exactly as the r⁰ Cauchy coefficient on |r| = 0.3e^{−β}: the c_b have poles up to r⁻⁸,
  which cancel. Two radii agree to 3·10⁻¹⁰.
* `L_ν` is **not** ∝ P_l: `L/P_l` varies with ν by O(1) and changes sign (`ee6` (a)). The EE leading-Regge position kernel is
  multi-harmonic, consistent with LEMMAS' "EE3 not per-harmonic".
* **Its edge coefficients are c_l × (combinations of B^{(b)} and F_b)**, so they are analytic in 0 < Im ν < 3 except at 3i/2 (+2p),
  where c_l vanishes.
* So the position-space de-shift condition for every EE component is the scalar one: **w c_l analytic in the strip apart from the target
  singularities, and w regular at 3i/2 for the k = 4 residue.**

**Uniformity on the box [D].**
* `P_even = Σ_k q_k cos^kθ (X₂₀X₄₀)⁻⁴`, with bounded q_k.
* Each insertion is `s^{−k}` times a Laurent polynomial `(1−t)^k(1+t)^{8−k}t^{−4} = O(s^k)`.
* Because B^{(b)}[n,m] is a polynomial of degree n in b, the k-th b-difference removes all orders n < k. The combination is therefore
  O(ρ^k) uniformly, with the same ν-analytic structure.
* At fixed ρ, the c_b tend to β-independent constants (`ee6` (b): 1.89·10³ / 50.3 / 6.46 at ρ = 0.2 / 0.5 / 0.85, β = 8–12).

---

## 2. De-shift survival for the EE kernel (`ee3_deshift.py`)

**Set-up (as lrc2).**
* `S = ηⁿ∫₀⁴⁰ w E_ν`, `S_LR = ηⁿ∫ w L_ν`, with `E_ν = K^EE/(N r²) = Σ_{b=0}^{4} C_b v^{b/2}Ĝ_b`.
* n = 1 (k = 2) or n = 3 (k = 4).
* 700 Gauss–Legendre nodes, radial order 260.
* Weights as lrc2: w2g (k = 2, entire), w2s (poles at ±i), w4g (k = 4, regular at 3i/2), w4p (pole at ±3i/2).

The table gives local exponents of (exact − LR) over β ≈ 8.5–13.5, at ρ = 0.2 / 0.5 / 0.85.

| weight | each b-block alone (b = 0…4) | EE kernel, pol. set 0 | EE kernel, pol. set 1 | raw | LR target |
|---|---|---|---|---|---|
| w2g | **−0.50** (all b, all ρ) | **−0.50 / −0.50 / −0.50** | **−0.50 / −0.50 / −0.50** | +0.50 | −0.50 |
| w2s | −0.37 … −0.41 (βη^{−½} form) | −0.41 / −0.40 / −0.39 | −0.38 / −0.39 / −0.39 | +0.50 | same form |
| w4g | **−0.50** to β ≈ 9.5, then the double-precision floor | **−0.50** (ρ = 0.5, 0.85, β ≤ 9.5); floor from β ≈ 7 at ρ = 0.2 | same | +2.50 | −0.50 |
| w4p | **+1.00** | **+1.00** | **+1.00** | +2.50 | +1.00 |

* **w4p residue formula, b-blocks.** `η³Re[2πi Res_{3i/2} w c_l e^{lβ}F_b]` reproduces S_b − S_LR to 1.000000–1.000007 at β = 10, for
  all b and ρ. The same residue for w4g is ≤ 10⁻²⁴.
* **The floor.** It is `10⁻¹⁶ × Σ_b|terms| × η³`, the same limitation as lrc2's w4g, here enlarged by the b-sum cancellation (≤ 10⁴ at ρ = 0.2).

**Verdict: LRC\* holds for the actual EE kernels [D + DV].**
* k = 2 keeps the full cancellation (J₀ < 3/2).
* k = 4 has exactly one residue at 3i/2, present iff the position weight has a pole there (J₀ < 2 generic; J₀ < 5/2 otherwise, the s = 1
  shell being raw η^{½}).

---

## 3. f_pos(±3i/2) for Q's cut-off weights (`ee4_fpos.py`)

**The criterion.** By LRC_CLOSE §3.1, f_pos(3i/2) ≠ 0 iff the position weight has a pole at 3i/2, since c_l has a simple zero there.
* The position weight is `χ_c × const/N_com`.
* 1/N_com is regular and non-zero near ±3i/2.
* By §1.4, the same condition governs the EE kernel.
* χ is even. **So the question is whether χ has a pole at −3i/2.**

**Exact algebra [P].**
* For the EE3 R_j: `R₋₃^{EE3}(−3i/2) = −524288` and `R₋₁^{EE3}(−3i/2) = 6193152/25`, both non-zero.
* For the scalar-type sc3: `R₋₃^{sc3}(−3i/2) = 0`, because a₋₃(l = 1) = 0. In EE3 this zero is cancelled by the (2ν + 3i) pole of the
  ratio R^{EE3}/R^{sc3}.

**Argument [D given DS].** HK Lemma S4 makes χR₋₃ analytic in −2 < Im ν < 0. So `χ(−3i/2) = (χR₋₃)(−3i/2)/R₋₃(−3i/2)` is finite, and
**f_pos(±3i/2) = 0**. The same holds for χ_c = χ·c, since the cutoff c is entire.

**Measurement on the stored weights [NE].** A pole of χ at −3i/2 would appear as `F₋₁(t) ~ −2πi ρ_χ R₋₁(−3i/2)e^{−3t/2}`.

| weight | local decay rate of F₋₁ (t = 15 → 36) | \|A(t)\| = \|F₋₁e^{3t/2}/(2πR₋₁)\| → ρ_χ | \|F₋₃\|e^{3t} |
|---|---|---|---|
| Q C4_value (XP, × erf cutoff 150/6) | 2.85 → **2.91** | 1.8·10⁻¹⁴ (t = 15) → **3.4·10⁻²⁷** (t = 36) | 2.0 (reproduces HK) |
| Q C4_deriv | 2.63 → **2.87** | **2.6·10⁻²⁵** | t² |
| *control:* sc3 value (scalar-type, double) | → **1.500** | → **9.26·10⁻⁷** (pole present) | 2.0 |
| *control:* sc3 deriv | → **1.504** | → **4.67·10⁻⁷** | t² |

* **Where the decay rate comes from.** The rate 2.91 is χ's nearest lower singularity, the root ν = −2.914i of 4ν² + 12iν − 1 (a zero of
  R₋₃^{EE3}). It lies below 3i/2, so for k = 4 its mirror gives at most a decaying η^{−0.41}.
* **The other poles χ can have.** In the strip they are at ∓i/2 and ∓5i/2. These are cancelled by the zeros of c_l and are harmless.

**Result: f_pos(±3i/2) = 0 for Q's C4_value and C4_deriv** [D given DS + NE], so the k = 4 criterion is **J₀ < 5/2**.
* **Stated identification.** The ν-normalisation of the Stage-4 EE3 functional is the scalar-type one (const/N_com). The project's
  integrands define it so (EE = scalar × ν-independent dressing). It is exact for EE1 (LEMMAS §3.2, 2.4·10⁻¹⁶). At k = 4 it rests on the
  same code path and was not rechecked.
* **The control.** Scalar-type weights *would* have f_pos(3i/2) ≠ 0. LRC_CLOSE's generic J₀ < 2 is therefore real; Q escapes it because
  of the EE3 structure.

---

## 4. Is the borderline R⁴ arc at a = b load-bearing?

**Where integrability enters.** LIGHTSIDE_DIST §1.1 defines the light value as the per-angle term-wise function paired with C_ν. So the
R⁴ arc must be C_ν-integrable near a = b (LEMMAS §4: growth Δφ^{−p}, p < 7/2) for that pairing to exist.

**What F_NORM's calibration gives [D].**
* Lemma A holds pointwise: `Ω₄[DT(W_χ)](a,b) = −p_χ ω₄[G_χ](a,b)` for a ≠ b.
* So the *double-trace part of the χ-exchange* is integrable near a = b, because it equals minus a heavy block's function.
* It does **not** transfer integrability to the R⁴ arc. The link between them is H_EFT(a), which F_NORM B3 shows is false pointwise in
  angle and holds only in ν-paired form.
* H_EFT(a)'s right-hand side, `𝒜₄[𝒞[g₄^χ]]`, is the paired R⁴ value. Stating H_EFT(a) therefore presupposes that this value exists.

**Decision.**
* The calibration does not *prove* a = b integrability, but it does not need a separate input either. **The a = b row is subsumed by
  H_EFT(a) (ν-paired)**, which is already on Q's list (row 6a, [D-sketch] for spin).
* A non-integrable arc would falsify H_EFT(a) as stated. It cannot silently change a number: (F)'s value 2g₄ is the calibrated ν-paired one.
* **Recommendation:** drop the separate [I] row for a = b, and append to H_EFT(a): "(including C_ν-integrability of the borderline R⁴ arc
  at a = b)".

**Attempt at a free-standing argument [I]** (not completed). By LEMMAS §4.1(3), Δφ → 0 is the light-cone edge of the Regge box at the
joint corner. The arc's growth p is then fixed by:
* the large-η behaviour of the R⁴ contact's top-Regge transverse profile in light-ray polarisation;
* convolved with the per-ν k = 4 kernel, whose raw η^{5/2} is untouched by §1–§2, since the pairing is per ν.

This needs the spinning top-Regge coefficient of R⁴ (F_NORM's route I, [I]).
* LEMMAS' heuristic p = 5/2 (margin 1) is neither improved nor tested here.
* The toy family cannot fix p (LIGHTSIDE P1–P3 differ in polarisation and are non-holographic).

---

## 5. F_NORM's deciding test for (P), the g₅ cancellation (`ee5_sd_run.py`, `ee5_sd_report.py`)

**What s_d is.** The chain (sdchain) computes
`s_d(D) = 0.04 Σ χ_d Re Γ_K`, with `Γ = m²(c T_ref ρ − T_st)` and `ρ = I(D)/I(D_ref = 2·10⁴)`.
* χ_d is the single, Δ_g-independent Stage-6 H0 weight.
* So s_d(D) is the c₂-unit correction of a J = 0 state of dimension D **at fixed functional**.

**Inputs.** The direct 7-d quadrature of the exact EE k = 4 integrand (no 1/Δ expansion), on grid n8-8-18-18-6-20-20, ν < 12 (99.7 % of
s_d):
* cached D = 300, 1000 and the reference;
* new runs at D = 3000 and 6000.

**The reference matters.** With the chain's finite reference, the raw s_d falls: 0.5079, 0.5064, 0.4813, 0.3960. That is an artefact.
x = m²/m_ref² = 0.02 / 0.09 at D = 3000 / 6000, and the chain subtracts its m-independent pieces after the m² scaling.

**Infinitely heavy reference.** Take r_D = I(D;ν)/I(D;0.02), Richardson-extrapolate r_∞ from D = 6000 and 2·10⁴, and run the identical
chain:

| D | m/Λ₃₀₀ | s_d^∞ (ν < 12) | s_v^∞ | windows [0,1) [1,2) [2,3) [3,4) [4,12) |
|---|---|---|---|---|
| 300 | 1.00 | **+0.50814** | −13.474 | +0.0989 +0.5703 −0.0458 −0.1021 −0.0132 |
| 1000 | 3.35 | **+0.50952** | −13.440 | +0.0990 +0.5701 −0.0455 −0.1012 −0.0128 |
| 3000 | 10.05 | **+0.50961** | −13.437 | +0.0990 +0.5701 −0.0455 −0.1011 −0.0128 |
| 6000 | 20.10 | **+0.50945** | −13.441 | +0.0990 +0.5701 −0.0455 −0.1011 −0.0130 |

* **Fit:** `s_d^∞ = 0.50957 − 128/D²`, with spread 1.5·10⁻³ over m/Λ = 1–20.
* **Robustness:** r_∞ taken from (3000, 2·10⁴) instead gives 0.50814 / 0.50952 / 0.50946 (r_∞ agrees to 1·10⁻⁸).
* **Record note.** The finite reference shifts the production s_d at D = 300 by −2.8·10⁻⁴, which is negligible against the margin of 1.27.

**Verdict [NE].**
* The ++ J = 0 derivative-piece correction does **not** vanish as m → ∞ at fixed functional. It is m-independent to 0.3 % from m = Λ to
  20Λ.
* By F_NORM §2.6, exact (P) (λ₅ = 0) needs it to vanish. **The exact-cancellation reading is excluded**, within the engine's premises
  (real section D, i.e. premise RD / H_RS; the CLS-D integrand definition; J = 0).
* The consistent reading is RU's double counting: σ_hU is the image of the uncancelled g₅ term. Counted once, the lower slope is
  **0.255** (CHLPS 0.26).
* **Unaffected:** the δ = 0.2 / 0.05 upper bounds (s_fail > 0 in every channel, F_NORM §2.6) and every headline number.
* **(P) for g₅:** [D-sketch with tension] → **[D-sketch + NE]: the g₅ light term is not cancelled.**

---

## 6. Final statements

### 6.1 LRC\* and LRC-a (replace LRC_CLOSE §6.4 rows)

| item | before | proposed |
|---|---|---|
| EE kernel structure | [D-sketch] (CPPR, unchecked) | **exact finite sum of 9 b-shifted scalar continuous-spin partial waves, ν-independent coefficients [P + DV, 1.6·10⁻⁷–1.5·10⁻⁶ end to end]** |
| LRC\* for EE | [D-sketch] | **[D + DV]**: exponents of (exact − LR) −0.50 (k = 2), −0.50 or +1.00 (k = 4: regular / pole at 3i/2), identical to the scalar |
| LRC-a, EE k = 4 | [D-sketch] | **[D]**: position LR form explicit [D + DV]. The heavy image of the k = 4 EE position kernel was not recomputed; the ν-normalisation identification is [D] (exact at k = 2) |
| f_pos(±3i/2), Q's C4 weights | — | **0** [D given DS + NE] |
| k = 4 criterion (Q) | J₀ < 2 | **J₀ < 5/2** |

### 6.2 H_K (replaces LRC_CLOSE §6.1's grade line)

> **H_K: [D given K1′, K2ν, K3, DS, LRC-a, LRC\*]**, with:
> * K1′, K2ν, K3 [D];
> * LRC\* [D + DV] for both scalar-type and EE kernels;
> * LRC-a [D + DV] (scalar; EE at k = 2) and [D] (EE at k = 4).
>
> Net: **[D]** for Result D; **[D given DS]** for Result Q (DS [NE]).
>
> Criteria: `J₀ < 3/2` (k = 2, 5, Ĉ₃) and **`J₀ < 5/2`** (k = 4 with Q's weights).

### 6.3 Sentence for the letter

> "For the analytically cut-off de-shifted functionals the heavy OPE sum may be interchanged with the Regge-corner and ν integrations
> by Tonelli's theorem, using only the Rindler bound J₀ ≤ 1 made uniform over the Regge box. The exact per-ν kernels are
> continuous-spin conformal partial waves; for the spin-2 structures they are exact finite sums of scalar continuous-spin partial waves
> with shifted external weights. At the light-cone edge, their corrections beyond leading Regge order resum into hypergeometric factors
> whose only singularity in the de-shift strip is cancelled by the leading-Regge coefficient, so the de-shift cancellation persists to all
> orders. The heavy sum then converges absolutely for J₀ < 3/2 (cubic functional) and J₀ < 5/2 (quartic functionals), given the
> existence of the de-shift weights, which we establish numerically for the quartic case."

Optional, for (P):

> "The heavy-state correction to the derivative sum rule is independent of the state's mass (checked from m = Λ to 20Λ), so the g₅ term
> is not cancelled on the light side; counting it once gives the lower slope 0.255."

### 6.4 Proposed hypothesis-table rows (UT §3.1; not edited here)

| # | before | proposed |
|---|---|---|
| H_K | [D-sketch] | **[D given DS]**; D: [D]; Q: DS [NE] the only non-derived input |
| LRC\* | [D + DV] scalar, [D-sketch] EE | **[D + DV]** both |
| LRC-a EE k = 4 | [D-sketch] | **[D]** |
| H7 §14 corollary (Ĉ′_EE, exact kernel; Result D) | [D-sketch] (EE step) | **[D + DV]**, margin ½ |
| a = b borderline arc (LSD; LEMMAS §4) | [I] | **subsumed by H_EFT(a) (ν-paired)**; free-standing [I] |
| (P) for g₅ | [D-sketch], tension | **exact cancellation excluded [NE]**; double counting → lower slope 0.255 |

---

## 7. What remains

1. **DS for Q** [NE]: the only non-derived input of H_K for Result Q (unchanged).
2. **Heavy image of the EE k = 4 position kernel.** This is a direct check that the plane-wave map of L_ν reproduces R_j^{EE3}. It would
   raise LRC-a for EE at k = 4 and the ν-normalisation identification of §3 to [DV].
3. **The spinning K1′ normalisation.** LEMMAS §2 is [D], but the spinning version was not checked numerically. Here the polarisations of
   each pair are normalised in that pair's rest frame, where the C_b stay bounded.
4. **a = b free-standing** (§4): needs the spinning top-Regge profile of R⁴ (route I) [I].
5. **(P):** only the derivative piece and J = 0 were run; the C5 piece and J = 2 were not.

## 8. Reproduce

    cd ads4/code/eestep
    ulimit -v 3000000; export OMP_NUM_THREADS=1
    python3 ee1_blocks_check.py                 # s     -> cache/ee1_blocks_check.txt
    python3 ee2_identify.py                     # ~80 s -> cache/ee2_identify.txt
    python3 ee3_deshift.py                      # ~6 min -> cache/ee3_deshift.{txt,npz}   ('report' re-prints from the npz)
    python3 ee4_fpos.py                         # ~20 s -> cache/ee4_fpos.txt
    python3 ee5_sd_run.py 3000 8,8,18,18,6,20,20 12 & python3 ee5_sd_run.py 6000 8,8,18,18,6,20,20 12; wait   # ~4 min each
    python3 ee5_sd_report.py                    # s     -> cache/ee5_sd_report.txt
    python3 ee6_lrprofile.py                    # ~30 s -> cache/ee6_lrprofile.txt

| file | content |
|---|---|
| `ee_blocks.py` | (0, b) continuous-spin radial recursion; Euler identity for b < 0; edge factor F_b; c_μ |
| `ee_quad.py` | exact H7-frame quadrature of S_ν; EE dressing τ₁₃D; Laurent coefficients c_b; t^b insertions |
| `ee1_blocks_check.py` | b = 0 vs lrc_blocks; DO PDE; polynomiality in b; poles; edge structure |
| `ee2_identify.py` | exact reduction; identification of Ψ^{(b)}; prefactor; end-to-end EE kernel |
| `ee3_deshift.py` | smeared exact vs LR (Cauchy) for b-blocks and the EE kernel; exponents; residue formula |
| `ee4_fpos.py` | F₋₁, F₋₃ decay for Q's XP weights (× cutoff) and the sc3 control |
| `ee5_sd_run.py`, `ee5_sd_report.py` | direct runs at D = 3000, 6000; s_d with finite and infinitely heavy reference |
| `ee6_lrprofile.py` | ν-profile of L_ν; boundedness of C_b along the edge |
