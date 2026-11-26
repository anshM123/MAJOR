# ANALYTIC_SLIVER2: pushing the analytic sliver certificate below L = 5 (Δ_g = 1000)

*26–27 Nov 2026. Research agent. Running document; unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*New code: `ads4/code/analytic2/` (imports `code/analytic/` read-only; outputs in `cache/`, logs in `logs/`).*
*Grades: **[P]** proof written here; **[certified]** Arb/exact scalar constants; **[DV]** numerics (sampled); **[conj]**.*
*Continues ANALYTIC_SLIVER.md (§ numbers "AS §x" refer to it) and takes REVIEW_SLIVER_L5.md into account (E1, G1, G2, F/den).*

---

## 0. Verdict

**L ∈ [3.65, 5] CLOSED (pending review).** Together with ANALYTIC_SLIVER's L ≥ 5 (and the rest above it), H_ang,b > 0 now holds
analytically for **every sliver state with L = −log d ≥ 3.65** (β_c ≥ 5.02): all a ≥ 998, every spin J ≥ 1001, b ∈ {0, 4}.

**Theorem (sliver, L ≥ 3.65).** For Δ_g = 1000, a ≥ 998, n = J − b ≥ 1001, b ∈ {0, 4} and L ∈ [3.65, 5]:
d³ℒ(G) ≥ M − (E₂ + E₃ + E_r + E_ρ) with M ≥ F_cert(β_c) and (E₂ + E₃ + E_r + E_ρ) ≤ ρ·F_cert, ρ < 1; hence H_ang,b > 0.

| L-cell | F_cert ≥ | contour ψ₂ + ψ₃ (/F) | residual E_r/F (largest piece) | log E_ρ2, log E_ρ3 | **ρ** |
|---|---|---|---|---|---|
| [4.5, 5] | 8.557e-7 | 2.5e-4 | 0.0265 (band [40, 898]: 5.1e-3) | −421, −497 | **0.0267** |
| [4, 4.5] | 1.209e-6 | 0.0325 | 0.0569 (cell [1031, 1032]: 5.1e-3) [fine run; coarse run: 0.378] | −180, −241 | **0.0894** (coarse: 0.410) |
| [3.75, 4] | 1.609e-6 | 0.0688 (right edge) | 0.154 (cell [1031, 1032]: 8.6e-3) | −39.3, −47.4 | **0.223** |
| [3.65, 3.75] | 1.879e-6 | 0.322 (right edge 0.320) | 0.222 (cell [1031, 1032]: 9.5e-3) | −28.0, −39.4 | **0.545** |

(F_cert is ≈ e^{−13}, so E_ρ is negligible down to 3.65.) Old values: AS §7.4 had ρ = 15.9 on [4.5, 5]; the old machinery gives ≈ 1450 on
[4, 4.5].

What made the difference:
1. **Lemma F** (§3): a certified turning-point bound for the per-power multiplier Ŝ(σ + iτ), uniform in spin, sharp to e^{2–4}
   (Lemma X loses e^{10–11}): Theorem E′ stopped K = 200 steps before the end, then the last K dual-Hahn steps through the turning point
   in Riccati frames that are affine in the parameters (τ, a, 1/N, σ), with a second-order defect bound.
2. **Lemma B** (§2): the connection index j (β-chain, TAIL E3) by chained similarity transforms. It sees the true decay
   |Ŝ(α − j)| ~ (κ/2)^j beyond the turning point (κ = ν/a ≈ 1.03 gives ×0.52 per step), where the absolute-value chain grows by ×1.22 per
   step; a multi-frame minimum removes the conversion cost near the double root κ = 1.
3. **Lemma Loc″** (§5): the sharp tail ratio |c_{j+1}/c_j| ≤ 1 for j ≥ J moves the floor of the remainder lemma from L ≈ 4.4 to L ≈ 3.6;
   below ν = 999 (ψ₂ = 0 above Λ) a separate truncation J′ = 2200 is allowed.
4. EPSW per-unit-cell residual masses on [1, 898) (same machinery as `rt_rmass2`), and the top line raised to Y = 6.25 below L = 4.
5. The referee's gaps are closed in the new bookkeeping (§6): G1 (θ-cutoff), G2 (no j-splitting above ν = 1150, a proved full-kernel bound
   with an e^{πν/2}-weighted CONTOUR tail), and F/den is asserted.

**Not closed: 0.06 ≲ L < 3.65.** The obstructions are in §8. The binding one is the remainder lemma (log E_ρ2 = +6.4, log E_ρ3 = +26 on [3.58, 3.65], against log F ≈ −13),
then the contour's own j-sums (the old near-real/top-line tails stop at L ≈ 3.58) and the right edge (0.32·F at 3.65).

---

## 1. Where the old bound lost, measured [DV] (`probe3.py`, `probe4.py`)

At the binding state (a = 998, b = 0, L = 4.5, n = 178175) and ν = 1031.5:
* exact |Ŝ(5/2 + iν)| = 1546 = e^{7.34}; sup over n (L ∈ [3, 8]) at a = 998 is e^{8.04}. Lemma X (AS §5.5) gives e^{18.46}:
  a loss of e^{10.4} (phase I adiabatic 4.1, window 2.6, σ-terms 2.4, other 1.3);
* exact |Ŝ(α − j)| decays by **0.52 per step in j** (j = 0…40; ratio 0.5226 → 0.46), as predicted by the frozen β-recurrence
  (below); the old chain bound (absolute values in TAIL E3) grows by ≈ 1.22 per step;
* exact kernel K(ν) = 2Re[BΣ_j c_jS(α−j)] = 3.10·10⁶, exact termwise sum |B|Σ|c_j||S(α−j)| = 2.05·10⁷, the certified bound used
  in AS §7.4 ≈ 8.5·10¹⁴. **So the termwise structure itself is fine at L ≈ 4–5; the loss was in the per-power and chain bounds.**

The frozen β-recurrence. With Ŝ(β) = d^βS(β), g = 1/(d²D), TAIL E3 reads Ŝ(β−1) = A(β)Ŝ(β) + B(β)Ŝ(β+1),
A = [1 + β(s−2β−1)gd]/[1 − β(β−s)gd²], B = β(β+1)g/[1 − β(β−s)gd²]. For β = σ + iτ, |β| ≈ τ ≈ 10³ and g ≈ 1/(4a²): A ≈ 1 + κ²d/2,
B ≈ −κ²/4 (κ = τ/a). The characteristic roots r± = (1 ± √(1 − κ²))/2 satisfy |r±| ≤ 1 for all κ, with |r±| = κ/2 for κ > 1. The
absolute-value bound |Ŝ(β−1)| ≤ |A||Ŝ(β)| + |B||Ŝ(β+1)| replaces them by the root (1 + √(1 + κ²))/2 ≈ 1.22 of the majorant.

---

## 2. Lemma B: the β-chain by chained similarities [P + certified arithmetic] (`bchain.py`)

**Lemma B.** Let x_j := Ŝ(β₀ − j) (j ≥ −1), V_j := (x_j, x_{j−1})ᵀ, so that V_{j+1} = M_jV_j with M_j = [[A(β₀−j), B(β₀−j)], [1, 0]]
(exact, TAIL E3). Let 𝔅 be a box of parameters (τ = Im β₀ ∈ [τ_lo, τ_hi], d ∈ [d_lo, d_hi], g ∈ [g_lo, g_hi], s ∈ [1, 9]) containing every
state of interest, and 𝕄_j ∋ M_j the Arb ball matrix over 𝔅. For **any** invertible point matrices T_0, T_1, … put
c_j := T_j⁻¹V_j. Then for every state in 𝔅
  ‖c_{j+1}‖_∞ ≤ ‖T_{j+1}⁻¹𝕄_jT_j‖_∞‖c_j‖_∞,  |x_j| ≤ ‖row₁(T_j)‖₁‖c_j‖_∞,  ‖c_0‖_∞ ≤ max_r(|(T_0⁻¹)_{r1}|M₀ + |(T_0⁻¹)_{r2}|M₁),
whenever |x_0| ≤ M₀ and |x_{−1}| ≤ M₁.

*Proof.* c_{j+1} = T_{j+1}⁻¹M_jT_jc_j and M_j ∈ 𝕄_j; ball arithmetic encloses T_{j+1}⁻¹M_jT_j, and the ∞-norm of an enclosure bounds the
norm of every member. The other two inequalities are V = Tc and c_0 = T_0⁻¹V_0. ∎

*State box.* t := a/N solves t² − (2+e)t + d(1+e) = 0 (e = s/N), so t/d is increasing in d, and g = (t/d)²/(a²(1+e)). For
a ∈ [a_lo, a_hi], d ∈ [d_lo, d_hi], e ∈ [0, 9·0.51d_hi/a_lo] this gives the Arb interval of g (`bchain.state_box`); a_hi = ∞ is allowed (g_lo = 0).

*Choice of frames* (only sharpness depends on it). T_j = normalised eigenvectors of the midpoint matrix when their condition number is ≤ 60,
else a Schur basis with the second column scaled by η = 0.15 (near the double root κ = 1). Result at a ∈ [998, 999], τ ∈ [1031, 1032],
L = 4.5 (M₀ = M₁ = 1): the bound decays by 0.54 per step (exact: 0.52); at a ∈ [1031, 1032] (κ ≈ 1) by ≈ 0.6 after an initial factor 11;
at a ≥ 10⁴ it stays ≈ 1.4 (true value ≈ 1). For comparison the old chain grows by 1.22 per step at every a.

---

## 3. Lemma F: the turning point by Riccati frames [P + certified arithmetic] (`tpframes.py`, `tpf3.py`)

Setting (b = 0, s = 1; b = 4 follows from the helicity bound AS §5.4). S_n(α) = ℒ_{n,N}((1−z)^{−α}) obeys v_{k+1} = M_kv_k,
v_k = (S_k, S_{k−1}), M_k = [[2w_k, 1], [1, 0]], w_k = ε_kα/2, ε_k = 2/(N−k) − 2/(N+k+1) (AS §5.1). In the variables m = N − k and
ι = 1/N: ε = 2/m − 2ι/(2 − (m−1)ι), and d_{k+1}/d_k = 1 − 2(1 − (m−1)ι)/(m(2 − (m−1)ι)) with d_k = 1 − k(k+1)/(N(N+1)).

Fix an integer K < n and put k_A = n − K. The last K steps have m = a + K, …, a + 1: **they depend on the state only through
(a, ι)**, and they contain the whole turning-point region when K ≳ 150.

**Phase I (Theorem E′ stopped at k_A).** With the eigen-coordinates c_k = P_k⁻¹v_k of AS §5.2,
  d_{k_A}^σ|c_{k_A}|_∞ ≤ Φ_I := c₀(2/(1 + √(1−κ′²)))^σ exp(q[arcsin κ′ + (q/2)log(1/(1−κ′²))]),  κ′ = τ/(a + K − 1), q = |α|/τ.
*Proof.* AS §5.2's proof, run only over the steps k < k_A: y_k ≤ τ/(N − k) < 1 there; the σ-excess is a sum over m = N − k ≥ a + K + 1,
bounded by σ∫_{a+K}^∞; the adiabatic sum over k ≤ k_A − 1 is bounded by q∫₀^{y_{k_A+1}}f with y_{k_A+1} ≤ τ/(a + K − 1) (the convexity step
f(y_{k+1})Δy_k ≤ ∫_{y_{k+1}}^{y_{k+2}}f reaches one index further). c₀ ≤ 1 + 2|α|/a² as in AS §5.2. ∎

**Lemma F.** Let 𝔅 be a box p = (τ, a, ι, σ) ∈ Π[p⁰_l − D_l, p⁰_l + D_l] with σ ≥ 0. Let R_{c,i}(p) = R⁰_{c,i} + Σ_l g_{c,i,l}(p_l − p⁰_l)
(c = 1, 2; i = 0..K) be any affine functions with point coefficients, F_i(p) = [[R_{1,i}, R_{2,i}], [1, 1]], Q_i(p) = F_i(p)⁻¹, and
δ_{c,i}(p) := 2w_i(p) + 1/R_{c,i}(p) − R_{c,i+1}(p). Then, for every state in 𝔅 and every n > K,
  |Ŝ_n(α)| ≤ Φ_I · ‖F_0⁻¹P_{k_A}(𝔅)‖_∞ · Π_{i<K}[(d_{k+1}/d_k)^{σ_lo}·e_i] · sup_𝔅(|R_{1,K}| + |R_{2,K}|),
  e_i := max_c sup_𝔅|R_{c,i}| + max_r sup_𝔅|(Q_{i+1})_{r1}|·Σ_c sup_𝔅|R_{c,i}δ_{c,i}|
(for one index i = i_s − 1, where the frame family changes, e_i := ‖Q_{i+1}(𝔅)𝕄_iF_i(𝔅)‖_∞ instead).
*Proof.* Column c of F_i is (R_{c,i}, 1)ᵀ and M_i(R, 1)ᵀ = (2w_iR + 1, R)ᵀ = R(R′, 1)ᵀ + (Rδ, 0)ᵀ, where R′ = R_{c,i+1}. Hence the exact identity
  M_iF_i = F_{i+1}diag(R_{1,i}, R_{2,i}) + [[R_{1,i}δ_{1,i}, R_{2,i}δ_{2,i}], [0, 0]],
so E_i := Q_{i+1}M_iF_i = diag(R_{1,i}, R_{2,i}) + Q_{i+1}[[R₁δ₁, R₂δ₂], [0, 0]], whose ∞-norm is ≤ e_i. With c′_i := Q_iv_{k_A+i},
c′_{i+1} = E_ic′_i, c′_0 = F_0⁻¹P_{k_A}c_{k_A} and S_n = (F_Kc′_K)_1; the d-ratios (≤ 1, σ ≥ σ_lo ≥ 0) carry d_{k_A}^σ to d_n^σ. ∎

**How the sups are bounded** (all Arb, `tpf3.lemmaF3_log`): |R(p)|² ≤ |R⁰|² + 2Σ_l|Re(R̄⁰g_l)|D_l + (Σ_l|g_l|D_l)² (phase changes of R do
not change |R|); δ by Taylor's formula to second order around p⁰ with the exact first derivatives ∂_lδ = 2∂_lw − g_l/R² − g′_l at p⁰ and
the Hessian ∂_l∂_mδ = 2∂_l∂_mw + 2g_lg_m/R³ enclosed over 𝔅 (∂²w from ε_m = −2/m² − 2ι²/B², ε_ι = −4/B², ε_mm = 4/m³ − 4ι³/B³,
ε_mι = −8ι/B³, ε_ιι = −8(m−1)/B³, B = 2 − (m−1)ι, ∂m/∂a = 1); Q_{i+1}(𝔅) is the Arb inverse of the affine ball matrix.
*Frames (sharpness only):* the Riccati recursion R_{i+1} = 2w_i + 1/R_i at p⁰ and its p-derivative g_{i+1} = 2∂w_i − g_i/R_i². Segment 1:
two forward ratios started at λ±(M_{k_A}(p⁰)) (travelling waves). Segment 2 (i ≥ i_s): the first forward ratio and the backward ratio
R_i = 1/(R_{i+1} − 2w_i) from the small eigenvalue at the end (the recessive solution beyond the turning point). The defect δ is then of
second order in the box size, which is what makes boxes of width ~1 in a and τ affordable (zeroth-order frames cost e^{+3…+9}, `tpframes.py`).

**Numbers [certified; DV comparison].** (`tpf3.py`, `check_F.py`.) σ = 5/2, τ ∈ [1031, 1032.1], ι over the L-cell [4, 4.5]:
| a-box | log bound (K = 200 / 400) | exact log|Ŝ| (sup over the box, sampled) | Lemma X |
|---|---|---|---|
| [998, 999] | 11.9 / 11.6 | ≈ 7.6 | 18.7 |
| [1020, 1021] | 7.9 / 7.7 | ≈ 3.5 | 18.7 |
| [1028, 1029] | 7.6 / 7.4 | ≈ 3.0 | 18.7 |
Point box at (a, τ, L) = (998, 1031.5, 4.5): 9.38 (K = 400) vs exact 7.34; the remaining gap is Φ_I (≈ 1.7) plus box widths.
Random check (10 states, a ∈ [998, 1100], τ ∈ [950, 1150], σ ∈ {¼, 5/4, 5/2, 7/2, 9/2}, L ∈ [3.5, 5]): every finite bound lies above the
exact value (gaps 2.0–4.4 in the log); for τ/a ≳ 1.08 with K = 200 the frame ball can become singular (bound = ∞), and the certificate then
falls back to Lemma X — only cells with residual mass < 10⁻³⁴ are affected.

---

## 4. The per-cell kernel bound (`rt2_kernel.K_cell`) [P + certified arithmetic]

For a real frequency cell [T₀, T₁] and an L-cell [L_a, L_b] (d ∈ [d_lo, d_hi] = [e^{−L_b}, e^{−L_a}]):
  |K_J(ν)| ≤ |B(ν)| Σ_{j<J₀}|c_j(ν)| d_hi^{j−5/2} Mdiff_j,  Mdiff_j = |Ŝ_j|(1 + τ_T|A₀|) + τ_T|A₁||Ŝ_{j−1}|/d,  Ŝ_j := Ŝ(5/2 − j + iν),
τ_T = 1/T ≤ 1/D = g d² (TAIL E2 dressing; A₀, A₁ as in `rt_cert.Abounds`). The caller multiplies by d³, and d³·d^{j−5/2} is increasing in
d, so the explicit powers are taken at d_hi. (Compared with AS §7 the "+1" in Mdiff is dropped: the residual pairs r with the full K, not
with K − K_∞.) Uniformity over a ≥ 998 is by a partition of [998, ∞) into ranges (unit boxes up to T₁ + 30, width 2 up to T₁/0.85,
then ratio 1.15, last [2·10⁵, ∞)); on each range:
1. base bounds M₀ ≥ |Ŝ(5/2 + iν)|, M₁ ≥ |Ŝ(7/2 + iν)| for both channels: b = 0 is the minimum of Theorem E′ at a = a_lo (valid for all
   a ≥ a_lo, E′ is decreasing in a) when τ ≤ 0.95a_lo, of Lemma X / E′ as in `rt_pp.PP0` (uniform in a ≥ 998), and of Lemma F on the box
   (a-range, τ-cell, ι-range of the L-cell; frames are attempted when τ/a_lo > 0.8 and E′ > e³); b = 4 by the helicity bound
   (AS §5.4, `rt_pp.PP` formula with d₀ = d_hi(1 + 10⁻⁴)). The ι-range is [t_lo/a_hi·(1 − 10⁻⁴), t_hi/a_lo] with t = a/N from the exact
   quadratic (both s, and the helicity-shifted state N + 4);
2. Lemma B for j = 1 … J_max − 1 (J_max = 60, 110 below L = 3.9), s ∈ [1, 9];
3. the tail j ≥ J_max (up to J₀): |Ŝ_{j}| ≤ G^{j−J_max+1}max(x_{J_max−1}, x_{J_max−2}) with G = max(2, (A+B)/den at |β| = hypot(J₀+3, ν)+3)
   (REC, AS §7; A, B, den increase with |β|), Mdiff_j ≤ (1 + τ_T(A₀ + A₁/d))max(|Ŝ_j|, |Ŝ_{j−1}|) at the largest |β|, and |c_{j+1}/c_j| ≤
   `ratio_bound_box`; geometric sum with ratio rb = ratio·d·G < 0.9 (asserted).
K_cell returns the maximum over the a-ranges. At the binding cell ν ∈ [1031, 1032.1], L ∈ [4, 4.5]: **1.25·10¹³ → (Lemma F + B)
8.6·10¹⁰** (old AS §7 bound 6.9·10¹⁵); the maximum now sits at a ∈ [998, 999].

## 5. Lemma Loc″ (sharp connection-series tail) [P + certified arithmetic] (`loc4.py`)

Lemma Loc′ (AS §7.3) bounds the connection-series tail Σ_{j≥J}|c_j||1−ξ|^{j−e₀} on the region G = {|1 − ξ| ≤ r} with the ratio bound
R_t = (1 + ν/(J+1))², which forces r·R_t < 1, i.e. r ≲ 0.35 at ν ≈ 1030, J = 1500. That is why Loc′ failed below L ≈ 4.4.

**Lemma (tail ratio).** For real ν ∈ [ν_lo, ν_hi] and every j ≥ J: |c_{j+1}(ν)/c_j(ν)| < 1, provided
  q := 4ν_lo² − 2ν_hi² − 3/2 > 0 and J²q ≥ (¼ + ν_hi²)².
*Proof.* |c_{j+1}/c_j| = ((j+½)² + ν²)/(√((j+1)² + 4ν²)(j+1)) ≤ ((x−½)² + ν_hi²)/(x√(x² + 4ν_lo²)), x = j + 1. Squaring, the bound is < 1
iff P(x) = −2x³ − qx² − 2x(¼ + ν_hi²) + (¼ + ν_hi²)² < 0; for x ≥ J the hypothesis gives qx² ≥ (¼ + ν_hi²)², and the other terms are
negative. ∎

Hence Σ_{j≥J}|c_j|r^{j−e₀} ≤ |c_J|r^{J−e₀}/(1 − r) for every r < 1, and the regions G₁, G₂ of Loc′ may be taken with |1 − ξ| up to ≈ 0.99.
Everything else (B-tiers with weights (1.02d/b)^{998}, angular fraction, dressing 41, K₀ ≤ 1.6, the small-ν pole handling) is
unchanged; the parameter search is extended to b₁ ≤ 0.4, θ₁ ≤ 0.5. The minimum with Loc′ is taken cellwise.

| L, J = 1500 | ν = 900 | 1000 | 1029 | 1040 | 1146 |
|---|---|---|---|---|---|
| 4.0, Loc′ (AS) | — | +29.8 | +95.0 | +119.7 | +551 |
| 4.0, Loc″ | −405.9 | −220.5 | −166.8 | −146.4 | +50.2 |
| 3.75, Loc″ | −219.3 | +5.2 | +70.4 | +95.1 | +307.7 (ν = 1149) |
| 3.75, Loc″, J = 2200 | | | −135.1 | | |
(log bounds of |𝒱_bℒ(R_J)| per unit cell; the f-mass at ν ≥ 1146 is ≤ e^{−9000}.)

**State relations at L ≥ 3.75** (needed by Loc′/Loc″, REC, the dressing and Lemma B) [exact algebra + exact-rational scan]:
t = a/N solves t = d(1 + e)/(2 + e − t), so d/2 ≤ t ≤ d(1 + e)/(2 − 0.012) ≤ 0.5031d(1 + e) for d ≤ e^{−3.75}; Dd² = a²(2 − t + e)²/(1 + e) ≥
3.95a²; N/a = 1/t ≥ 84. Scan (`check_state.py`, exact rationals, 96 twists incl. 10⁹, both channels): a/(Nd) ∈ [0.50085, 0.50300],
Dd²/a² ≥ 3.9528, N/a ≥ 84.5. K₀ = 1.00033 (bound 1.6) and the I_W lower bound (ratio 2.10) at L = 3.75, 4, 4.25 [DV, mpmath].

## 6. Residual bookkeeping, with the referee's gaps closed [P + certified arithmetic] (`rt2_resid.py`)

E_r = d³(2/Λ)[2Σ_bands m_band·sup_band|K_J| + 2·sup_{ν<Λ}|K_J|·1.66·10⁻²³] + E_far, where
* bands: [0, 0.05] (maximum principle) and [0.05, 0.4] exactly as AS §7; [0.4, 40] and [40, 898] with EPSW's certified range masses
  and K_cell on sub-cells of width 2 (for L < 4: [0.4, 1] with EPSW's mass of [0, 1], then the unit cells of `rmass3.py` on [1, 898)); [898, 1059) the 161 unit cells of `rt_rmass2.json` (EPSW panel masses regrouped; reproduced by the
  referee); [1059, 1150] with EPSW's Gaussian tail mass 5.9·10⁻²⁴⁷ (K_cell without frames, cells of width 13);
* **G1 (θ-cutoff).** H is built from f_θ, the EPSW masses are for f. The difference contributes ≤ 2 sup|K_J|·∫|f_θ − f| with
  ∫|f_θ − f| ≤ 1.66·10⁻²³ (`tailcert/cache/tc_cutoff.txt`), needed only for ν < Λ: above Λ, ψ₂ = 0, so r_θ = f_θ and |f_θ| ≤ |f| = |r|
  pointwise (0 ≤ θ ≤ 1 on the real axis) and the EPSW masses already bound ∫|r_θ|; above 1150 f_θ is used directly;
* **G2 (far frequencies).** Above ν = 1150 nothing is split: the term is d³(2/Λ)∫_{1150}^∞|f_θ||𝒱_bK(ν)|dν with the full kernel bounded by
  the Loc′ representation, |𝒱_bℒ(F_{3,ν})| ≤ K₀·41·sup|F_{3,ν}| ≤ 1.6·41·(π/2)e^{πν/2}(0.235d_lo)^{−3} (FULLPOS3 Lemma B; |1 − ξ| ≥
  0.49d − 0.255d on the enlarged disc), and ∫_{1150}^∞|f_θ|e^{πν/2} ≤ `wtail(1150)` = e^{−7298}: CONTOUR's tail bound (cb_lines.tail) with the
  weight folded into every Gaussian, ∫_T^∞e^{wt−(t−c)²/2s²}dt = e^{wc+w²s²/2}s√(π/2)erfc((T−c−ws²)/(s√2)), and the Y-term
  ∫_T^∞e^{wt−t²/4} ≤ e^{w²}(2/(T−2w))e^{−(T−2w)²/4}. Correspondingly E_ρ (§5) is integrated only over ν < 1150;
* **F/den.** The geometric tails of AS §7's near-real sums need (A+B)/den ≤ 2 and F/den ≤ 1 at |β| = hypot(J+3, 1151) + 2; both are
  asserted in `rt2_loc.py` (values at L = 4: 1.97 and 0.93; the assertion uses 1.99).

**Lemma B, refinement (`bchain.chain_bounds_multi`).** Near the double root (κ ≈ 1) the eigen/Schur frames of §2 pay a conversion factor
‖T_0⁻¹‖ ≈ 10–20 at j = 0 (the phases of x_0, x_{−1} are unknown), and the termwise sum is dominated by j ≈ νd/2 ≈ 5–10, where that
factor has not yet decayed. Since Lemma B holds for *every* frame family, the elementwise minimum over several families is again a
bound. Added families: T_j = U_j diag(1, η_j), U_j the unitary Schur basis, η_0 = 1, η_{j+1} = max(η_min, qη_j), (q, η_min) ∈
{0.5, 0.7} × {0.03, 0.15}; the step norms are max(|r₁| + |t|η_j, |r₂|η_j/η_{j+1}), with no initial conversion cost. At τ ∈ [1003, 1004.1],
a ∈ [998, 999], L ∈ [4, 4.5] (M₀ = M₁ = 1) the bound at j = 6 drops from 22.7 to 0.74. (Used for L < 4.5; the [4.5, 5] certificate and the coarse [4, 4.5] run used `chain_bounds`, which is also rigorous.)

---

## 7. Certificates (`rt2_contour.py`, `rt2_loc.py`/`rt2_loc2.py`, `rt2_resid.py`, `rt2_cert.py`; `cache/rt2_cert_*.txt`)

Route T″ = Route T′ (AS §7) with: E_r by §4/§6, E_ρ by Loc″ (§5) on ν < 1150, the unsplit full-kernel bound above 1150, G1 and
F/den included. Contour pieces unchanged (`rt_cert.piece_*`, top line at Re α = ¼ − k₀).

```
Route T'' certificate, L in [4.5, 5.0] (beta_c in [5.881, 6.383]): F_cert >= 8.5569e-07
   psi2: bottom 1.35e-05  left 1.34e-05  top 2.07e-04  right 2.03e-05  (Y = 4.25)
   psi3: bottom 1.12e-08  left 1.13e-08  top 1.21e-09  right 6.42e-08
   residual E_r/F = 2.6483e-02  (small-nu 2.0e-03+5.8e-03, G1 5.0e-17, far nu>=1150 7.0e-303); worst cells: [40.0,898.0] 5.05e-03,
            [0.4,40.0] 1.26e-03, [1001.0,1002.1] 8.25e-04, [1002.0,1003.1] 8.09e-04, [1003.0,1004.1] 7.35e-04
   Loc'' remainders (J = 1500): log E_rho2 = -421.2, log E_rho3 = -496.5
   rho <= 2.6737e-02   PASS
Route T'' certificate, L in [4.0, 4.5] (beta_c in [5.377, 5.881]): F_cert >= 1.2087e-06
   psi2: bottom 1.28e-05  left 1.35e-05  top 3.00e-02  right 2.46e-03  (Y = 4.25)
   psi3: bottom 1.72e-08  left 1.86e-08  top 4.30e-07  right 1.28e-05
   residual E_r/F = 3.7767e-01  (small-nu 1.7e-03+5.6e-03, G1 3.7e-16, far nu>=1150 2.2e-302); worst cells: [40.0,898.0] 9.71e-02,
            [1002.0,1003.1] 2.19e-02, [1003.0,1004.1] 2.00e-02, [1004.0,1005.1] 1.71e-02, [1005.0,1006.1] 1.65e-02
   Loc'' remainders (J = 1500): log E_rho2 = -180.4, log E_rho3 = -241.0
   rho <= 4.1014e-01   PASS
```
("far" is the log-bound e^{−7273} floored at e^{−690} before the float conversion.)

Rerun of [4, 4.5] with the final settings (`RT2_FINE=1`, multi-frame Lemma B; `cache/rt2_cert_4.0_4.5.txt`, the coarse run kept as
`cache/rt2_cert_4.0_4.5_coarse.txt`):
```
   residual E_r/F = 5.6893e-02  (small-nu 1.7e-03+5.6e-03, G1 1.7e-16, far nu>=1150 2.2e-302); worst cells: [1031.0,1032.1] 5.14e-03,
            [1030.0,1031.1] 2.51e-03, [1025.0,1026.1] 2.45e-03, [1029.0,1030.1] 2.34e-03, [1026.0,1027.1] 2.29e-03
   rho <= 8.9364e-02   PASS
```
```
Route T'' certificate, L in [3.75, 4.0] (beta_c in [5.124, 5.377]): F_cert >= 1.6086e-06        [k0 = 4, Lemma B multi-frame,
   psi2: bottom 1.12e-05  left 1.22e-05  top 5.32e-05  right 6.83e-02  (Y = 6.25)                 EPSW per-cell masses on [1, 898),
   psi3: bottom 1.92e-08  left 2.15e-08  top 1.09e-09  right 4.56e-04                              J = 1500 below nu = 999, J' = 2200 above]
   residual E_r/F = 1.5374e-01  (small-nu 1.5e-03+5.0e-03, G1 5.7e-16, far nu>=1150 3.5e-302); worst cells: [1031.0,1032.1] 8.55e-03,
            [1029.0,1030.1] 6.22e-03, [1027.0,1028.1] 6.19e-03, [1028.0,1029.1] 6.14e-03, [1030.0,1031.1] 5.87e-03
   Loc'' remainders (J = 1500): log E_rho2 = -39.3, log E_rho3 = -47.4
   rho <= 2.2261e-01   PASS
```
```
Route T'' certificate, L in [3.65, 3.75] (beta_c in [5.023, 5.124]): F_cert >= 1.8785e-06       [same settings as [3.75, 4]]
   psi2: bottom 1.02e-05  left 1.13e-05  top 3.52e-04  right 3.20e-01  (Y = 6.25)
   psi3: bottom 1.92e-08  left 2.19e-08  top 8.49e-09  right 2.36e-03
   residual E_r/F = 2.2200e-01  (small-nu 1.3e-03+4.6e-03, G1 8.3e-16, far nu>=1150 4.1e-302); worst cells: [1031.0,1032.1] 9.53e-03,
            [1030.0,1031.1] 6.76e-03, [1029.0,1030.1] 6.70e-03, [1006.0,1007.1] 5.45e-03, [1005.0,1006.1] 4.94e-03
   Loc'' remainders (J = 1500): log E_rho2 = -28.0, log E_rho3 = -39.4
   rho <= 5.4472e-01   PASS
```
For [3.75, 4] the remainder is `rt2_loc2.py` (per-cell f-masses on [898, 1150), J′ = 2200 for ν ≥ 999; REC at J = 1500: (A+B)/den = 1.970,
F/den = 0.930), the residual is run with `RT2_FINE=1 RT2_JP=2200`, and the top line is at Re α = ¼ − 4 (k₀ = 4).
 For comparison, AS §7.4 had ρ = 15.9 on [4.5, 5]
and the old machinery gives ≈ 1450 on [4, 4.5]. On [4, 4.5] the turning band ν ∈ [998, 1035] contributes ≈ 0.27 in total (≈ 0.02 per cell),
the EPSW band [40, 898] 0.097 (one sup over the whole band), the top line 0.030.

**DV checks of the new lemmas.** `check_F.py` (30 random states, a ∈ [998, 1100], τ ∈ [950, 1150], σ ∈ {¼, 5/4, 5/2, 7/2, 9/2},
L ∈ [3.5, 5], boxes as in the certificate): every finite Lemma F bound exceeds the exact value, by 1.87–4.4 in the log (some
τ/a ≳ 1.07 boxes return ∞ and fall back to Lemma X). `check_B.py` (a ∈ {998} ∪ [998, 5000], b ∈ {0, 4}, ν ∈ [300, 1100],
L ∈ [3.75, 5], exact M₀, M₁): max_j log(exact/bound) ∈ [−1.8, −0.41] — Lemma B is within e^{0.4–1.8} of the truth.

---

## 8. Obstructions below L = 3.65 (measured)

1. **The remainder lemma (Loc′/Loc″) is now the floor.** Integrated E_ρ (rt2_loc2: J = 1500 below ν = 999 with per-cell f-masses,
   J′ = 2200 above): log E_ρ2 = −39.3 on [3.75, 4], −28.0 on [3.65, 3.75], but **+6.4 on [3.58, 3.65]** (log E_ρ3 = +26.4, the ψ₃ piece)
   and +87.4 / +107 on [3.5, 3.65];
   F ≈ e^{−13.2}. The binding cells are ν ∈ [994, 998.5] (ψ₂ still nonzero, so the contour's J = 1500 is forced there) and, for E_ρ3,
   the whole ψ₃ range near Λ. Per unit cell at L = 3.5: ν = 990: +30.7 (J = 2200), ν = 1029: +116 (J = 2200), +197 (J = 1500), +595
   (J = 3000). The loss is structural: on the Cauchy circle |z| = ρ = n/N the bound uses |(1 − ξ)^{∓iν}| ≤ e^{π|ν|/2} and
   |F_{p,ν}| ≤ (π/2)e^{π|ν|/2}|1−ξ|^{−p} (= e^{1616} at ν = 1029), which must be beaten by |c_J|r^J, (1.02d/b)^{998} and the angular fraction
   ((1−ρ)/q)^{a+1}; all three weaken as d grows. A remainder bound below L ≈ 3.6 needs a representation of ℒ(R_J) that keeps the
   phase of (1 − ξ)^{−iν} (steepest descent in θ), or no j-truncation at all (κ^M of AS §3b, which has no remainder but no bounds yet).
2. **The contour's J is pinned at ≈ 1500** by the near-real tails of AS §7 (they use (A+B)/den ≤ 2 and F/den ≤ 1 at |β| ≈ hypot(J+3, 1151);
   at J = 2200: 2.67 and 1.60). For ν ≥ 999 this is avoided here (J′ = 2200 is used only where ψ₂ = 0), but cells ν ∈ [950, 998.5] at
   L = 3.5 need J > 1500 on the contour as well; the near-real tails then have to be redone with a growth constant G > 2 and the
   inhomogeneity F/den > 1 (routine, not done).
3. **The termwise hump at moderate twist.** With Lemmas F and B the turning point is no longer special; what remains is the intrinsic
   termwise factor Σ_j|c_j|d^j r₊(κ)^j ≈ e^{r₊νd/2} (r₊ = (1 + √(1−κ²))/2 → 1 at large a, where the true kernel is O(1)). At L = 3.75,
   ν ≈ 900 it is e^{≈10}; at L = 3 it is e^{≈22}. The EPSW per-unit-cell masses on [1, 898] (computed here, `rmass3.py`, total 1.739·10⁻⁷ +
   6.085·10⁻⁷, matching EPSW's ranges) buy a factor ≈ 2·10³ against the band sup and are what make [40, 898] harmless at L = 3.75.
   Below L ≈ 3 a resummed large-a bound is needed: at large a, Ŝ_j = 1 + Δ_j with Δ_j small, and Σ_jc_jd^j = ₂F₁ is O(1) exactly,
   so K_J = Bd^{−α}[₂F₁(a₀, a₀; b₀; d)_trunc + Σc_jd^jΔ_j] would remove the hump wherever |Δ_j| ≪ 1 (κ ≲ 0.3); not done.
4. **Contour pieces (AS's `rt_cert.piece_*`, unchanged here).** With k₀ = 4 (Y = 6.25) the top line is ≤ 4·10⁻⁴·F down to L = 3.65 (k₀ = 2:
   1.0·F at L = 3.75), but (a) the right edge ν = Λ + iy (turning point, Lemma X + the old absolute-value chain) grows fast:
   0.0025·F (L = 4), 0.068·F (3.75), 0.32·F (3.65); and (b) at L ≈ 3.58 the old top-line j-sum tail (`jsum_osc`, 59 explicit terms,
   geometric ratio 2d·sup|c_{j+1}/c_j| < ½) fails its own assertion. Both are routine to fix with Lemmas F and B and a longer explicit
   j-range (as done for the residual), but were not redone here.
5. **State relations** (0.49d ≤ a/N ≤ 0.51d, Dd² ≥ 3.9a²) are verified here for L ≥ 3.75 (§5); the referee finds they first fail near
   L ≈ 2.5–3, so they are not binding before items 1–3.

---

## 9. Grades, inputs, files

* **New and unreviewed here:** Lemma B and its multi-frame refinement (§2), Lemma F with its phase-I truncation of Theorem E′ (§3),
  the per-cell kernel assembly (§4: a-range partition, helicity use, REC tail with G), the tail-ratio lemma / Loc″ (§5), the state relations
  at L ∈ [3.75, 4.5] (exact algebra + scan), the J′ ≠ J splitting for ν ≥ 999 (§6, rt2_loc2), the weighted CONTOUR tail `wtail` and the
  unsplit far-frequency bound (G2), the G1 bookkeeping.
* **Taken from AS (unreviewed there, reviewed in REVIEW_SLIVER_L5.md):** Theorem E′, the helicity bound, Lemma NB, Lemma V, Loc′'s regions and
  constants (K₀ ≤ 1.6, dressing 41), the splitting lemma, the contour pieces (`rt_cert.piece_*`), Lemma X (used only as a fallback in cells
  with residual mass ≤ 2·10⁻³⁵ or beyond, and in the contour's right edge as in AS).
* **Taken as certified from elsewhere:** LRINF_TIGHT F_cert cells, FULLPOS3 §4.1(c) and Lemma B, EPSW panel machinery and masses
  (`rt_rmass2.json`, reproduced by the referee; `rmass3.py` runs the identical computation on [1, 898)), CONTOUR tail data, TAIL_CERT
  E1–E5 and tc_cutoff.
* **Rigor notes.** All special-function / matrix / polynomial bounds are Arb balls converted to floats by `.upper()`; positive bounds are
  accumulated in floats with a final (1 + 10⁻⁹). Frames (Lemmas B, F) are floating-point *choices*; rigor does not depend on them (they
  enter as exact point matrices, inverted in Arb). The minimum over several rigorous bounds is taken cellwise.

Files (all new, `ads4/code/analytic2/`):
* lemmas: `bchain.py` (Lemma B), `tpframes.py` (phase I; zeroth-order frames, superseded), `tpf2.py` (superseded), `tpf3.py` (Lemma F),
  `loc4.py` (Loc″);
* assembly: `rt2_kernel.py` (K_cell), `rt2_resid.py` (E_r incl. G1/G2, `wtail`), `rt2_contour.py`, `rt2_loc.py`, `rt2_loc2.py`, `rt2_cert.py`,
  `rmass3.py` (EPSW per-cell masses on [1, 898));
* checks / probes: `check_F.py`, `check_B.py`, `check_state.py`, `probe1–4.py`, `probe_hump.py`, `partial.py`, `run_resid_old.py`,
  `run_loc3_probe.py`, `fastS.py`, `common2.py`;
* outputs: `cache/rt2_cert_*.txt`, `cache/resid_*.json`, `cache/loc_*.json`, `cache/contour_*.json`, `cache/rmass3_*.json`; logs in `logs/`.
