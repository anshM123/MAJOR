# P2_CONTRACT: the true three-channel EE contraction at O(1/m²), the kernel's own O(r) term, the flat tower of P4, and Theorem W's J\*

*23 Nov 2026. Unreviewed, not committed. No existing file was edited; `REGGE_TREE.md` and `code/reggetree/` were not touched.*
* *New code:* `ads4/code/p2con/`, reproduced by `run_all.sh` (≈ 4 CPU-h, at most 2 processes). Caches are in `cache/`, whose `.gitignore`
  ignores everything except itself.
  * `pc_engine.py`: an exact leading-Regge EE plane-wave image engine for arbitrary polarisation tensors and polynomial insertions.
  * `pc_sym.py`, `pc_gensym.py`: exact symbols (P2_SPINNING's code path) at every sample state, made drift-free.
  * `pc_contract.py`, `pc_run.py`: the contraction, de-shift reconstruction and exact costs.
  * `pc_ksh.py`, `pc_checks_ksh.py`, `pc_ksh_run.py`: the kernel's own O(r) term.
  * `pc_p4.py`, `pc_flat.py`, `pc_p4b.py`: the flat tower and the P4 remainder.
  * `pc_jstar.py`: Theorem W's J\*.
* *Read-only imports:*
  * `h7i14/v3_regge_kernel`, `v4_hclosed`, `spinfunc/s1_offdiag_kernel` (split representation, EE (2,4) structure);
  * `p2spin/*` (structures, symbols, CHMRS seed and mechanism, ancillary-series parser);
  * `eestep/ee_quad` (exact integrand, for checks);
  * `n2sub/n2_common`, `lsw/lsw_common` + `lsw/cache/costs_1000.json`, `pwimg/cache/pw3_cost.json`, `pw4_theoremW.json`;
  * the CHMRS ancillary `regge_moments_order8.m` (scratch copy).

Grades as in UNIFIED_THEOREM §0. Conventions as in LARGESPIN_W, PW_IMAGES and P2_SPINNING:
* f = F/c is the de-shifted image per unit leading image (f₀ = P_β);
* "c₂ units" means δI = κ S/m²;
* 'same' is the project's +−, and 'opp' is ++.

---

## 0. Summary

### 0.1 Constants: model against computed (Δ_g = 1000, all 47 sample states of Theorem W, β = 0.043–0.525)

All κ below are coefficients of the canonical ring-class (polynomial) part of f₂ = F₂/c. They are defined by the large-μ expansion of the
exactly reconstructed F₂ (§1.5) and are quoted relative to CHMRS's scalar κ^sc (κ_A = 6.46, κ_B = 11/6, κ_C1 = 1/6, κ_C2 = 1/6,
κ_C3 = −1.96).

| quantity | P2 model [M] | **computed, block side** | **computed, kernel O(r) term** | **computed, total** | grade |
|---|---|---|---|---|---|
| Λ_B (ν²P), +− (and odd J) | 1 | **1 exactly** (κ_B = 11/6; the symbol adds 0) | κ_B = +2 − O(β²) (1.75–2.00) | **κ_B = 3.58–3.83 ⇒ Λ_B = 1.95–2.09** | [DV] |
| Λ_B, ++ | 1 | **1 exactly** | κ_B = −2 + O(β²) | **κ_B = −0.17…+0.08 ⇒ \|Λ_B\| ≤ 0.09** | [DV] |
| Λ_C2 ((η²−1)ν²P′) | 1 | **1 exactly**, both helicities | 0 | **1** | [DV] |
| ν⁰ P′(η) (C1/C3), +− | Λ_C1 = 1, Λ_C3 = 2.03 | κ = −35.9…−31.9 (the flat d⁴⁴ helicity term, §2.1) | 0 | **Λ_C1 ≈ 190–215** (in C1 units) | [DV] |
| ν⁰ P′(η), ++ | as above | κ = +0.17…+0.62 | 0 | Λ_C1 ≈ 1–3.7 | [DV] |
| Λ_A (ηP), +− / ++ | 3.17 | κ_A/η = −65…−59 / +50…+59 | ±(14.5–16.5) | **Λ_A ≈ 6.6–7.8 / 5.6–6.5** (harmless: ∝ leading image, W1c) | [DV] |
| off-diagonal (++ ↔ +−) O(1/m²) | 0 | ≤ 5·10⁻¹² of the diagonal | 0 | **0** | [DV] |
| O(1/m) image (all entries, after drift removal) | antisymmetric, ∝ leading | **0** (≤ 10⁻¹¹) | — | **0** ⇒ the (S) failure cannot act at O(1/m²) | [DV] |
| monopole Q = F₂(i/2) | — | **∓12.800** (+− / ++), β-independent | **±12.800** | **0: PW-R1 holds at O(1/m²) for the full ⟨TTTT⟩** | [DV] |
| K_sh (kernel index-shift / non-ring) | 10–50 [M] | — | the kernel's O(r) term is **ring class** (ν²P, ηP), not index-shifted; its non-polynomial tail costs ≤ 10⁻⁴ of it | **K_sh premise removed** | [DV] |
| exact total O(1/m²) cost / H (max over entries) | Λ_amp·Σ\|κ^sc\|c_s/H | — | — | **≤ 0.41 (J = 26), ≤ 0.11 (J ≥ 46), ≤ 0.05 (J ≥ 64), ≤ 0.005 (J ≥ 200)** | [DV] |
| beyond-hybrid O(1/m²) cost (flat O(y) tower subtracted) | — | — | — | **≤ 0.09 × the Λ_amp = 1 scalar load** | [DV] |
| K₄ (P4) | ≤ 3.8 [M] "if the ν⁶ tower is in the hybrid" | — | — | flat tower: subtraction exact in the flat-scaling limit (§2). **Remainder over the full ν range: K₄^rem ≤ 2.9 for J ≤ 50 (one state 26.6), 10–16 at J = 64, 20–7·10³ at J ≥ 80, i.e. above tolerance** | [DV] (series costs); P4 [I] |

**Prominent negatives.**
1. **Λ_B is NOT 1 for the +− helicity once the kernel's own O(r) term is included: it is ≈ 2.0** (κ_B ≈ 23/6 instead of 11/6).
   * The block side alone gives exactly 1, as the model claimed.
   * The extra ±2 is the kernel's O(r) term, which Theorem W had filed under K_sh.
   * This is worse than the model value for +−, and better for ++.
2. **The ν⁰ P′ class for +− is ~200× its scalar value.** It is the flat d⁴⁴ helicity term (§2.1), so the model's direct pairing was wrong
   there; P2's 'conj' variant (Λ_C1 ≈ 205) was right. It is harmless in Theorem W's costs.
3. **P4 is not closed (now the only binding premise).** In the flat-scaling limit the hybrid floor subtracts the flat tower exactly.
   * At O(y) this is exact, including spin, but only after the kernel's O(r) term is added (§2.1).
   * At O(y²) it is exact for ++ up to a residual +(3/16)J₀ y².
   * But the AdS/flat kinematic mismatch of that tower (P_ν(cosh β) against J₀(b̂p)), costed over the full range ν ≤ Λ as P4 requires,
     exceeds the K₄ tolerance at every sampled J ≥ 80. **With the computed remainder there is no J\*; the passing sampled J are
     {30, 33, 38, 42, 50, 64, 301}** (§2.2, §4).

### 0.2 Bottom line
* **Task 1 (done, [DV]).** The O(1/m²) correction to every sample state was computed through the exact leading-Regge EE kernel (all three
  polarisation channels, no helicity model), with the exact drift-free spinning symbols, every leakage mode and polynomial insertion, the
  CHMRS seed, and the kernel's own O(r) term.
  * The de-shifted F₂ is reconstructed exactly: fit residual ≤ 10⁻¹², held-out nodes ≤ 10⁻¹⁰.
  * F₂ is costed exactly with PW4's formula; the X-route agrees to ≤ 10⁻⁷.
  * In Theorem W's convention (full O(1/m²), no hybrid subtraction) the computed correction uses ≤ 41 % of the 2×2 budget H at J = 26 and
    ≤ 5 % at J ≥ 64. **P2 (all of it: Λ_amp, K_sh, P2′, P3) is replaced by a computation.**
* **Validations.**
  * (i) The leading image is diagonal with equal class normalisation (10⁻¹⁰), reproducing EE_FUNCTIONAL.
  * (ii) The engine reproduces v3/s1 heavy images (ratio 1.000000) and the CHMRS ∂_p mechanism for the scalar seed (10⁻¹⁵).
  * (iii) **The flat limit of the computed spinning O(1/m²) terms reproduces exactly the O(y) terms of the flat helicity partial waves d⁰⁰
    and d⁴⁴ of the hybrid floor**, including the −31.8 P′ helicity term. A spinning flat-space fact is recovered from AdS symbols.
  * (iv) The +− value κ_B ≈ 3.8 sits inside the (1100, 100) window fits B = 3.7–5.0, and the ++ value ≈ 0 is consistent with the
    (364, 64)@300 fit B = 0.11 [NE; which helicity each fit measured is not re-checked].
* **Task 2.** The flat ν⁶ tower **is** part of the hybrid floor's flat partial wave at the flat-scaling level. Showing this needs the
  kernel's O(r) term: it converts the AdS moment's bulk-point envelope E₂(y) into the hybrid's stripped (2−y) envelope, with opposite signs
  per helicity. P4 should therefore bound only the remainder. That remainder is small (K₄ ≲ 3) up to J ≈ 64 but not beyond when costed
  over the full ν range (§2.2), so **P4 is [I]**. This is the full-range artefact LARGESPIN_W already flagged, now quantified for the
  subtracted remainder.
* **Task 3 (J = 0 engine O(a²)): not done** (§3).
* **Task 4.** Theorem W restated (§4) with premises **P0, P1′, P4 only**.
  * J\* = 26 if K₄ = 0; 30 at K₄ = 1; **38 at K₄ = 4 (model value), 50 at K₄ = 10**; 80 at K₄ = 10 with margin ×2.
  * The sampled computed P4 remainder exceeds tolerance at J ≥ 80, so **no unconditional J\* follows**.

---

## 1. Task 1: the true three-channel contraction

### 1.1 The engine (`pc_engine.py`) [P-structure + DV]
* **Kernel.** The raw EE image of a plane wave e^{p·x+p̄·y} with constant tensors on each leg and polynomial insertions Q(x)R(y) is the
  split representation of H7 §18.1 with the EE (2,4) structure V₀⁴(z₂′·z₄′)² (`s1.A_monos_struct('EE')`) and the light-ray (1,3) side
  (`v3.B_monomials`): `raw = −η^{μν}∂_μ∂̄_ν ∫Dz₀ Â₆ B̂₆`.
  * This is the exact leading-Regge kernel for every polarisation. **It contains all three channels of LARGESPIN_W §9.1**, since it has no
    helicity truncation.
* **Each side, in its own rest frame.** On the celestial circle S = 1 there, so the side is an exact trigonometric polynomial. Its harmonics
  come from an exact 32-point DFT; |k| ≤ 4 (A) and ≤ 6 (B) are populated, and the tail is ≤ 10⁻⁹.
* **Derivatives.** The exact order-8 Taylor jet of S^{−λ}R^κ in 3 variables, contracted with the operator polynomial.
* **Legs.** The null-square basis n(t_j)⊗n(t_j), t_j = 2πj/5. A mode e^{ikt} has weights e^{ikt_j}/(5c_{|k|}) with c = (3/2, −1, 1/4).
* **Celestial pairing.** M_{kk′}, an SO(2,1) matrix element, by spectral quadrature for νβ < 12 and by the closed form
  `2π√2^{n+n′}/((1−δ)_n(δ)_{n′}) Σ_j C(n,j)C(n′,j) j! η_AB^j η_A^{n−j} η_B^{n′−j} P_l^{(n+n′−j)}(η)` above. The two agree to
  ≤ 2·10⁻¹² where both are stable.
* **Checks.**
  * The engine against `v3.heavy_image` and `s1.img2(EE)`, for six polarisation configurations, ν = 0.7 and 2.3, β = 0.3–1.4: ratio
    1.000000. Helicity flips give 0, as they should for EE.
  * The engine's helicity image equals 2π(R₋P_{l−1} + R₊P_{l+1}) with the **published** R∓ (kernels.json) exactly (ratio 2π to 10⁻¹³),
    so the published normalisation is inherited.
  * **Insertions:** the scalar seed s₂ (CHMRS T, Δφ = 3) through the engine plus the de-shift reproduces the CHMRS ∂_p mechanism
    (s3_mech, a = 7) to 10⁻¹⁵ in both the P and the P′ coefficient (β = 0.81, ν = 8, 15, 30).

### 1.2 Symbols (`pc_gensym.py`, `pc_sym.py`)
* **Generation.** The exact O(1/m²) symbols of all four structure families at each of the 47 sample β, by P2_SPINNING's own code path.
  Interpolation residual 0. The q = 3/2 set is identical to `p2spin/cache/s2_main.pkl`.
* **Drift removal (W1c).** The structure basis is changed by U(ε) = 1 + εU₁ + ε²U₂, so that no X,Y-constant component remains along the
  leading class tensors at orders 1 and 2. Drift constants run from −13.6 to 69.8.
* **Rest frames and seed.** X = B(β/2)X_rest and Y = B(−β/2)Y_rest, with the weight 2^{−deg}.
* **Block tensor of entry (h, h′).** Σ_o L̃_h^{(o)} ⊗ (−1)^{n−o} L̃_{h′}^{(n−o)} + s₂ L̃⁰⊗L̃⁰.

### 1.3 Leading order and O(1/m) [DV]
* **Leading order.**
  * M₀ is helicity-diagonal (off-diagonal ≤ 8·10⁻¹²).
  * The same and opp images are both 32·(R₋P_{l−1} + R₊P_{l+1}), to ≤ 7·10⁻¹⁰.
  * The odd-J structure has normalisation 3.58.
* **O(1/m).** After drift removal the O(1/m) image vanishes identically in every entry (≤ 3·10⁻¹¹ relative), for all 47 states.
  * So **P1 holds at the level of the exact leading-Regge kernel**: P1 is [DV] at leading-Regge kernel level, and still [D-sketch] beyond the leading-Regge
    kernel. The kernel has no O(r^{1/2}) term, so any O(1/m) piece must come from the block side.
  * The O(1/m²) term g₁B₁ that the failure of (S) (P2_SPINNING §7.2) would allow is therefore 0. **(S) is no longer needed.**

### 1.4 De-shift reconstruction [DV]
* **Fit.** By PW2, raw = R₋^{sc}F(ν+i) + R₊^{sc}F(ν−i). F is fitted as a(μ)P_μ + b(μ)P′_μ, with a, b in
  span{1, μ², 1/(4μ²+(2j+1)²), j ≤ 7}, on 130 nodes ν ∈ [0.05, 1500].
* **Quality.**
  * Residual ≤ 7·10⁻¹³ (diagonal entries) and ≤ 10⁻¹¹ (symbol pieces).
  * Held-out nodes (ν = 0.37 … 777) ≤ 1.4·10⁻¹⁰.
  * The off-diagonal is numerically 0 (raw ≤ 10⁻¹¹ of the diagonal), so its relative fit is meaningless.
* **Robustness.** The polynomial (large-μ) coefficients and F on the real line are stable to ≤ 10⁻⁶ under j_max = 5 … 12. F(i/2) is
  evaluated through the analytic limit b·P′ → −r₀^{(b)}/(4(1+η)).
* **Analyticity.** a has no pole at μ = i/2 (coefficient ≤ 10⁻⁷). The would-be poles of b at i/2 and of a, b at 3i/2 cancel in F
  (res_a·η + res_b = 0, checked). This reproduces PW2's analyticity of F in |Im μ| ≤ 3/2.

### 1.5 Results: block side (`cache/run_*.log`)
* **ν² classes.**
  * F₂'s μ² coefficients are exactly 32·(11/6)η (P) and 32·(1/6)(η²−1) (P′), for both helicities, for odd J, and at all 47 states.
  * The symbol pieces have μ² coefficient ≤ 6·10⁻¹³; the whole ν² content is the seed.
  * **Λ_B = Λ_C2 = 1 exactly for the block side** [DV]. P2's model statement is confirmed, now without the model.
* **ν⁰ classes (asymptotic polynomial part).**
  * +− (and odd): ηP −65…−59 per η; P′(η): −35.9…−31.9.
  * ++: ηP +50…+59 per η; P′(η): +0.17…+0.62.
  * The seed through EE gives −13.3η (against the scalar mechanism's +14.46η). The difference is the EE multiplier tail, which is ν-flat
    at leading order.
* **Monopole Q = F₂(i/2) = −12.800 (+−) / +12.800 (++)**, constant over all β. The seed gives +88 to +100, the symbol −101 to −113.
* **Non-polynomial remainder** (multiplier tail plus monopole): its cost is ≤ 10⁻³ of the entry's total.

### 1.6 The X-route check
The total cost computed as 2Re Σ X 𝒜[F/c] dν/Λ (published X) agrees with 2[Σψ₂F/c + D·F(i/2)] to ≤ 10⁻⁷ (`pc_jstar` consistency line).
This is PW4's formula on the actual spinning corrections.

### 1.7 The kernel's own O(r) term: "K_sh" computed (`pc_ksh.py`, `pc_checks_ksh.py`, `cache/ksh_*.log`)
* **Structure [D + DV].** In the H7 frame:
  * the scalar factor S_ν has **no O(r) term** (m₂₀m₄₀ = (1−r²)² − 4r²(x̂·n)², Nn = −4z₀·x̂(1+O(r²)), c₁₃ = 4r² exactly);
  * τ₁₃ is exact and X₂₀X₄₀ = 1 + O(r²);
  * the only O(r) term is in the (2,4) numerator:

        D = 16(z₂′·z₄′)² + 128 r (z₂′·z₄′)[(z₂′·n)(z₄′·x̂) − (z₂′·x̂)(z₄′·n)] + O(r²).

  * **Checked against `ee_quad`'s exact integrand:** the O(r) coefficient of log S is 0, with the residual shrinking ×10 when r shrinks
    ×10. The D₁ formula matches with an O(r²) residual.
* **Split form.** The O(r) kernel is a split form with one extra vector contracted between the sides:
  A-side 128 e₂₄(z₄′·x) against B-side B₀[|y|(z₂′·n_θ)], minus the swap.
  * The new B monomials are checked against direct θ-quadrature to 6·10⁻⁹.
* **Image on the leading block (47 states).** It is **ring class, not index-shifted**.
  * κ_B = +2 (+−, odd) and −2 (++), with O(β²) corrections (1.75–2.00); κ_C2 = 0; ν⁰ P′ = 0; κ_A/η = ±(14.5–16.5).
  * **Q = ±12.800, cancelling the block monopole exactly**, so the total F₂(i/2) = 0.
  * The non-polynomial tail costs ≤ 7·10⁻⁹, against totals ≈ 10⁻⁵–10⁻⁴.
  * The off-diagonal is 0.
* **Interpretation (§2.1).** This ±2ν²P is exactly the difference between the bulk-point envelope E₂(y) = (2−y)/(1−y)² of the AdS moment
  and the stripped envelope of the functional, with the helicity-dependent sign of the graviton stripping.
* **Consequences.**
  * **PW-R1 holds for ⟨TTTT⟩ at O(1/m²): the monopole cancels between block and kernel** [DV]. PW_IMAGES' "[I] whether the full sum
    cancels" is answered: yes, at this order.
  * Theorem W's K_sh premise is replaced by this computation.

## 2. Task 2: the flat tower and P4

### 2.1 What the hybrid floor contains (`pc_p4.py` (a), `pc_flat.py`)
* **Method.** The flat image used by the floor (`omegau/nf_common.flat_parts_nf`) is (2−y)D̃_J(y)/2, with D̃ = P_J(1−2y) (++) and
  P^{(0,8)}_{J−4}(1−2y) (+−), and y = ν²/m². It is expanded at fixed b = (2J+1)√y (extrapolation in J; mpmath, 40 digits, at O(y²)).
  The result is compared with the flat-scaling limit of the AdS series: P → J₀(b), P′(η) → −(ν²/b)J₁(b), η²−1 → b²/ν².
* **O(y), exact [DV]** (fit residual 5·10⁻⁴ on 14 values of b):

  | object | O(y) coefficient |
  |---|---|
  | flat ++, E₂·d⁰⁰ | (11/6)J₀ − (1/6)bJ₁ − (1/6)J₁/b |
  | flat ++, hybrid (stripped) | −(1/6)J₀ − (1/6)bJ₁ − (1/6)J₁/b |
  | flat +−, E₂·d⁴⁴ | (11/6)J₀ − (1/6)bJ₁ + 31.834 J₁/b |
  | flat +−, hybrid (2−y)d⁴⁴/(1−y)⁴ | (23/6)J₀ − (1/6)bJ₁ + 31.834 J₁/b |
  | **AdS block side** (computed, §1.5) | (11/6)J₀ − (1/6)bJ₁ − κ_C J₁/b, with κ_C = +0.17 (++) and −31.87 (+−) |
  | **AdS block + kernel O(r)** | adds ∓2J₀ (++ / +−) |

  * The block side reproduces E₂·d^J exactly per helicity. Block plus kernel reproduces **the hybrid floor's flat part exactly, both
    helicities**.
  * The C1-type ν⁰P′(η) terms are flat-type (P′(η) ~ ν²/b), so the −31.8 of +− is the d⁴⁴ helicity term.
  * **The O(1/m²) flat tower is therefore exactly in the hybrid floor** [DV].
* **O(y²), scalar seed (++).**
  * The flat-scaling part of CHMRS's O(1/m⁴) series, extracted symbolically, is
    `g₂ = −b²J₀/72 + (599/240)J₀ − 7J₀/(120b²) − (11/30)bJ₁ − 7J₁/(20b) + 7J₁/(60b³)`, which contains the −(η²−1)ν⁶P/72 tower.
  * The flat E₂·d⁰⁰ partial wave equals g₂ **+ (3/16)J₀** at b = 0.8, 1.7, 2.7, 3.6: difference 0.1875J₀ to 4 digits.
  * The O(y²) kernel term (O(r²)) and the spinning O(1/m⁴) symbols were not computed. For +− the flat O(y²) coefficients are large (≈ 15–80
    in scalar units, from the d⁴⁴ helicity terms).
* **Verdict.**
  * The ν⁶ tower is part of the flat partial wave, and the hybrid subtracts it exactly at flat-scaling order: [DV] for ν⁶ and the full ++
    tower up to the (3/16)J₀ y² residual, whose K₄-equivalent is ≤ 0.3 (1.7 and 0.63 at two states).
  * At O(y) the subtraction is exact for both helicities, but only with the kernel O(r) term included.
  * At O(y²) for +−: [D-sketch] (same mechanism expected; not computed).

### 2.2 The P4 remainder (`pc_p4b.py`)
**Definition.** R₄(ν) = A₄^{CHMRS}(ν, β) − ν⁴ g₂(b̂p). This is the full scalar O(1/m⁴) series minus its flat tower in the hybrid's own
kinematics. It is costed exactly (de-shift form; both terms entire) relative to the K₄ unit 4a = η²ν⁴P.

| J | K₄ of R₄ (beyond-hybrid) | K₄ of the raw ν⁶ term (P2 §7.1) | K₄ tolerance (H − d₂)/c₄ |
|---|---|---|---|
| 26 | 0.04 / 1.07 / 0.16 | 35 / 88 / 7 | 18.8 / 0.66 / 1.03 |
| 30–42 | 0.03–2.0 | 21–86 | 1.3–30 |
| 46 | **26.6** / 0.01 / 0.42 | 242 / 85 / 75 | 4.6 / 35 / 27 |
| 50 | 2.9 / 0.74 / 1.36 | 93–135 | 12.5–17.4 |
| 64 | 16.1 / 10.5 / 14.0 | 93–224 | 16.3 / 15.9 / 15.1 |
| 80 | **163** / 35 / 20 | 224–385 | 21 / 28 / 35 |
| 100–150 | 15–5.5·10³ | 370–960 | 44–286 |
| 200–344 | **1.5·10³–7.2·10³** | 1.4·10³–3.4·10³ | 336–1090 |
| odd 33 / 51 / 101 / 201 / 301 | ≤ 0.4 / ≤ 17 / 56–949 / 726 / 3040 | 37–2650 | 5–6220 |

* The Δφ-dependent (AdS) part alone, 5189(η²−1)ν⁴P/1440, has K₄ ≤ 0.84.
* **Verdict.** Subtracting the flat tower lowers the requirement from P2's K₄ ≈ 35–3000 to ≲ 3 for J ≤ 50, except three states:
  * J = 46 at m = 1.045Λ: K₄ = 26.6 against a tolerance of 4.6;
  * J = 26 at 1.128Λ: 1.07 against 0.66;
  * J = 34 at 1.033Λ: 2.03 against 1.74.
* **At J ≥ 80 the remainder is as large as the raw tower.** The AdS/flat kinematic mismatch ν⁶[(η²−1)P_ν(cosh β) − (b̂p)²J₀/ν²], of
  relative size β², dominates once it is integrated up to ν ≈ Λ ≈ m, where the fixed-β 1/m series is not valid.
* **P4 cannot be closed with the truncated series on the full range [I].** It needs a window-limited statement: the series for ν ≪ m and an
  exact or hybrid control beyond.
* This also applies to the O(1/m²) beyond-hybrid content, but there it is tiny: its cost is ≤ 0.09 of the scalar load (mode 'bh').

## 3. Task 3: the J = 0 engine O(a²) comparison — not done
* **Why not.** The engine (`repro_engine/integrand.py`) is the k = 4 functional at J = 0, where β = O(1/m). It lies outside the fixed-β
  expansion used here (P2_SPINNING §7.3).
* **What it would take.** A k = 4 version of `pc_engine`: Lemma-B prefactor, k = 4 radial weights and the light-ray x₀ integral normalised
  as the engine. The β → 0 limit of the fixed-β symbols is not uniform (drift constants ∝ 1/β appear).
* **In its place, three independent validations of the contraction:**
  * the exact scalar mechanism (10⁻¹⁵);
  * the exact flat helicity partial waves at O(y), both helicities, block and kernel (§2.1);
  * the window B fits (§0.2 iv).

## 4. Theorem W with the computed constants (`pc_jstar.py`, `cache/jstar_{full,block,bh}.txt`)

**Theorem W (restated).** Let O be a heavy primary in H_LJ(1000) with couplings λ ∈ ℂ² (even J) or ℂ¹ (odd J). Assume:
* **(P0)** the leading-Regge hybrid floor (unchanged);
* **(P1′)** the O(1/m) part of M is a congruence. This is now **[DV] at leading-Regge kernel level** (§1.3); beyond it [D-sketch];
* **(P4)** the beyond-hybrid O(1/m⁴) remainder plus beyond-all-orders costs at most K₄·Σc₄/m² per entry.

Then M_Ω(O) ⪰ 0 for every coupling vector on every sampled state with J ≥ J\*, where the O(1/m²) correction is the **computed** one
(§1.5 + §1.7; exact costs, both forms below). **P2, P2′, P3 and K_sh are no longer premises.**
* **Grade.** [P given P0, P1′, P4] + [DV] (computed costs on sampled states). Between sampled states: [NE].

**J\*** (sign-free form: d₂ + K₄Σc₄/m² ≤ θH; d₂ = max over entries of the computed |cost|):

| O(1/m²) used | K₄ = 0 | K₄ = 1 | **K₄ = 4** | K₄ = 10 | K₄ = 35 | K₄ = 4, θ = ½ | K₄ = 10, θ = ½ |
|---|---|---|---|---|---|---|---|
| computed, block + kernel (Theorem W convention) | **26** | 30 (signed 2×2: 26) | **38** | **50** | 100 | 50 | 80 |
| computed, beyond-hybrid (flat O(y) tower subtracted) | 26 | 26 | 38 | 50 | 100 | 50 | 80 |
| block side only (kernel term omitted) | 26 | 30 | 38 | 50 | 100 | 50 | 100 |
| *P2 model (P2_SPINNING §7.5, K_sh = 10)* | — | — | *38* | *50* | — | *50* | — |
| **state-wise computed scalar P4 remainder (§2.2)** | **no J\*.** Passing sampled J: 30, 33, 38, 42, 50, 64, 301; failing: 26 (one mass), 34 (one mass), 46 (one mass), 51, and every J ≥ 80 except 301 | | | | | | |

* **K_sh = 10 or 30** (no longer a premise; kept as a margin) changes no entry at θ = 1.
* **Reading.**
  * With the true O(1/m²) the 2×2 budget is never threatened: d₂/H ≤ 0.41 at J = 26, ≤ 0.05 at J ≥ 64, and 0.0007 at J = 344.
  * Everything is decided by P4. J\* = 38 / 50 holds if K₄ ≤ 4 / 10.
  * The computed series remainder satisfies that for J ≤ 64 (with three exceptions at one mass each), but not at J ≥ 80, where the
    full-range series is not a valid representation.
  * Λ_eff (computed / scalar load at Λ_amp = 1, R1′) spans 0.02–7.7 over states. It is large only where the scalar load is small by
    cancellation (J = 38, m = 1.245Λ).

## 5. Residual list (exact)
1. **P4 [I]:** a full-range bound on the beyond-hybrid O(1/m⁴) and beyond-all-orders remainder.
   * The series remainder is costed in §2.2: small for J ≤ 64, large for J ≥ 80, where the series is not valid at ν ~ m.
   * Needs a window-limited formulation (series on ν ≪ m, exact/hybrid control above).
   * Not computed: the spinning O(1/m⁴) symbols through the true kernel and the kernel's O(r²) term. For +− the flat O(y²) tower has large
     helicity terms; they are in the hybrid if the §2.1 mechanism persists, which is [D-sketch].
2. **(3/16)J₀ y²**: an unexplained flat-tower residual between CHMRS's O(1/m⁴) series and E₂·d⁰⁰. K₄-equivalent ≤ 0.3 (1.7 at one
   state). Origin [I] (normalisation or the kernel's O(r²) term?).
3. **P1′ beyond leading-Regge kernel** [D-sketch]. The kernel has no O(r^{1/2}) term, so O(1/m) can only come from the block side, which
   is shown to vanish with the leading kernel.
4. **P0** unchanged (hybrid floor; certified flat part plus exact LR curvature).
5. **Task 3 (J = 0 engine O(a²))** not done.
6. **Between sampled states** [NE]. The constants are smooth in β: κ_B, κ_C2 exact; Q = ∓12.800 constant; the kernel κ_B = ±(2 − O(β²)).
7. **Helicity assignment of the window B fits** (+− ≈ 3.8 against ++ ≈ 0) not re-checked against DIRCV/N2's data conventions [NE].

## 6. Suggested inventory changes (not applied)
* **P2_SPINNING §3.2's [M]** is replaced:
  * the block-side Λ_B = Λ_C2 = 1 is [DV];
  * the kernel's own O(r) term adds ±2ν²ηP (+− / ++), so for +− the total is **Λ_B ≈ 2**;
  * C1 for +− is ~200× scalar (the d⁴⁴ term, flat).
* **LARGESPIN_W:**
  * K_sh is not an index-shifted remainder but a ring-class ±2ν²P term, now computed;
  * P2, P2′ and P3 are replaced by computed costs;
  * J\* = 38 (K₄ = 4) / 50 (K₄ = 10) with P4 the only magnitude premise.
* **PW_IMAGES** "[I] whether the ⟨TTTT⟩ sum cancels the monopole": it does at O(1/m²) (block −12.8, kernel +12.8; +− sign convention).
* **UNIFIED_THEOREM H_LJ(1000):** stays [D-sketch + M] → **[D given P4] + [DV]** on the sampled states with J ≥ 38 (K₄ ≤ 4) or J ≥ 50
  (K₄ ≤ 10). P4 is [I], and its computed series remainder fails at J ≥ 80. **Do not write "proved for J ≥ 38".**

## 7. Reproduce
`cd ads4/code/p2con && bash run_all.sh <path/to/regge_moments_order8.m>`. Outputs:
* `cache/run_*.log`: per-state block side;
* `cache/ksh_*.log`: kernel O(r) term;
* `cache/checks_ksh.txt`, `flat.txt`, `p4.txt`, `p4b.txt`: §1.7 and §2;
* `cache/jstar_{full,block,bh}.txt`: §4.
