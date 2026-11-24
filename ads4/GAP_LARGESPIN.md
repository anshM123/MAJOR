# GAP_LARGESPIN: the large-spin / large-impact-parameter part of Ω^EE positivity, uniformly in Δ_g ≥ 1000

*24 Nov 2026. Unreviewed. No existing file was edited, nothing was committed. Companion of `GAP_SCALING.md` (scoping; read after its
§3.1 was written here independently, and used to avoid duplication — see §0.3).*
* *New code:* `ads4/code/gaplarge/` (`run_all.sh`; ≈ 5 CPU-h in total, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`).
  Outputs in `code/gaplarge/cache/`, whose `.gitignore` ignores everything but itself.
* *Read-only imports* (no file of theirs written, no bytecode): `lrtight/` (lrt_common, lrt2_m1 + its caches), `exactskel/` (h13_common,
  h13_cert_m1), `deshift/d_abel`, `n2sub/n2_common`, `lsl/lsl_core`, `lsw/lsw_common`, `largespin/ls_core`, `p2con/` (pc_run, pc_jstar,
  run/ksh caches), `p4win/` (p4w_common, p4w_contour, p4w_states, p4w_L), `agproof/` (ag_flat, ag_chi + `chi.json`).
  Where those modules hard-code Δ_g = 1000 (module constants `L`, `NU`, `PS`, `X`, `D1000`, `Dg`), they are **re-set in my process**
  (`gl8.set_gap`, `gl6`), never edited; at Δ_g = 1000 every re-run reproduces the stored numbers (§2, §4).

Grades: [P] proved; [certified] = Arb ball arithmetic + proved monotonicity; [D] derived; [D-sketch]; [DV] derived, float numerics;
[NE]; [M] model; [A] assumption; [I] open. Notation as LRINF_TIGHT: `Λ = √((Δ_g−2)(Δ_g−1))`, `Λ₀ = Λ(1000) = 998.49987`,
`sinh(β/2) = (J+½)/m` (m in AdS units), `b̂ = (2J+1)Λ/m = 2Λ sinh(β/2)`, `μ = m/Λ`, `u = 1/μ² ≤ 1`, `F_LR = 2Ω₂ + uσΩ₃`.

---

## 0. Verdict

### 0.1 What now holds for every Δ_g ≥ 1000

| # | statement | grade | coverage in (J, m/Λ, Δ_g) |
|---|---|---|---|
| **U1** | **Lemma LR∞ for every Δ_g ≥ 1000, no upper limit.** For all Λ ≥ Λ₀, u ≤ 1, σ = ±1: `F_LR ≥ (1−r)A₁,min M₁(β)/Λ² + A₂′/(2b̂³)` with **r < 1 on β ≥ 0.37**, and **F_LR > 0 on β ≥ 0.245** (LRINF_TIGHT's own thresholds, now Δ_g-uniform). Moreover `r(β; Λ) ≤ (Λ_a/Λ)² r̂(β; Λ_a)` for Λ ≥ Λ_a (any anchor Λ_a ≥ Λ₀; r̂ = the bound with the envelope R̂, §1.2), so r falls like Λ⁻² at fixed β and the thresholds improve with Δ_g (U1′). | **[P + certified]** | every m ≥ Λ (u ≤ 1), both helicities, every **J ≥ J\*(m) = ⌈(0.37211·m − 1)/2⌉** (m in AdS units), i.e. `J + ½ ≥ 0.18606·(m/Λ)·Λ`; twist edge: J\*/Λ = 0.2233 (1000) → **0.22387** (Δ_g → ∞) |
| U1′ | Δ_g-improved certified thresholds (β₀ for r < 1 / for F_LR > 0): **0.2375 / 0.1325 for Δ_g ≥ 2606; 0.1675 / 0.10 for Δ_g ≥ 5000; 0.1175 / 0.10 for Δ_g ≥ 10⁴; 0.10 / 0.10 for Δ_g ≥ 2·10⁴** (0.10 = end of the certified M₁ cells) | [certified] | as U1 with s₀ = 2 sinh(β₀/2) |
| **U2** | **Δ_g-uniform gap-unit floor:** `Λ²·F_LR(β, u, σ; Λ) ≥ f_U(β)` for every Λ ≥ Λ₀, with f_U > 0 on every β-cell ≥ 0.245 and **f_U ≥ 0.99·A₁,min M₁(β) on β ≥ 0.34** (A₁,min M₁ = 13.34 at β = 0.37, 12.24 at 1, 9.11 at 2, 1.98 at 5) | [certified] | β ≥ 0.245, all u ≤ 1, σ |
| **U3** | **Hybrid curvature window (ROUTEA) uniform in Δ_g:** `curv ≥ A₁,min M₁(β)/Λ² − R > 0`, hence (flat NF02 ≥ N/5 [certified]) **hybrid floor > 0**, on **b̂ ∈ [114.06, 2Λ] for all m ≥ Λ; [110.74, 2Λ] for m ≥ 1.3Λ; [107.51, 2Λ] for m ≥ 1.5Λ** — for every Δ_g ≥ 1000. New Λ-uniform remainder bound for ROUTEA's region I (§3.2). | **[D + certified]** (derivation [D] as ROUTEA; all constants Arb) | every J with b̂ in the window |
| **U4** | **Lemma C′ is monotone in Λ** at fixed (J, m/Λ, h) (§3.1), so its Δ_g = 1000 conclusion (hybrid floor > 0 for b̂ ≤ 119 at m/Λ ∈ {1.3, 1.5, 2, 3, 5, 10, 20}) holds at every Δ_g ≥ 1000, with a reach in b̂ that grows with Λ | [P] (monotonicity) + LARGESPIN_ANALYTIC's [P + DV] | b̂ ≤ 119, sampled m/Λ ≥ 1.3 (sampling as at 1000) |
| **U5** | **Coverage of the leading-Regge/hybrid floor at every Δ_g ≥ 1000:** m ≥ 1.3Λ: hybrid floor > 0 on **every b̂ ≤ 2Λ** (U4 on b̂ ≤ 119 at C′'s sampled masses and b̂-samples, ∪ U3 on [110.7, 2Λ]; the old m ∈ [1.3, 1.5)Λ hole (119, 245) is closed); b̂ ≥ 0.2456Λ: F_LR > 0 (U1). m ∈ [Λ, 1.3Λ): hybrid floor > 0 on [114.1, 2Λ] (new, every m), F_LR on b̂ ≥ 0.2456Λ; b̂ < 114 only via GAP_SCALING's Lemma-C region (small J). | as U1–U4 (C′ part: sampled as at 1000) | — |
| U6 | **Theorem W / Theorem L′ premise tolerances move favourably with Δ_g in chart B** (fixed β, m/Λ): every computed load relative to the floor falls (L′: ∝ Λ^{−2.3}; W: ∝ Λ⁻¹), every tolerated premise coefficient grows (L′: Λ_amp ∝ Λ²; W's EDGE-R at fixed x₁: ∝ Λ^{3.1}). §0.2, §2, §4. | [DV] given the premises transported with Δ_g-independent AdS-unit constants [M] | 16 L′ + 47 W sampled states, images at 2000 and 5000 |

### 0.2 What moves with Δ_g (computed at 2000 and 5000; [DV] given premises transported with Δ_g-independent AdS-unit constants [M])

| quantity (fixed gap-unit state: same β, same m/Λ) | Δ_g = 1000 | 2000 | 5000 | law |
|---|---|---|---|---|
| Theorem L′ monopole tolerance Q_tol (uniform certified floor) at the J\* image (β = 0.37, 1.20Λ) | 3.35·10⁻³ | 3.34·10⁻³ | 3.33·10⁻³ | Λ-independent (Λ²D → −1000: −995.3 / −999.1 / −1001.3) |
| Theorem L′ tolerated Λ_amp, monopole (unif. floor), J\* image | 714 | 2851 | 17804 | ∝ Λ² |
| Theorem L′ tolerated Λ_amp, full scalar-shaped vector R1′, J\* image | 126 | 482 | not converged (§2.2 ‡) | ≈ Λ² (1000 → 2000) |
| Theorem L′ / Lemma P4W (uniform rate AG_U): max d/(F/2) over the 16 L′ images (unif. floor) | 0.190 (0.133 with LRINF's F_cert) | 0.035 | 0.0050 | ≈ Λ^{−2.3} |
| — min tolerated edge multiple (+− transferred scalar envelope) | 12.6 (18.9) | 228 | 9.05·10³ | grows |
| Theorem W (P4W-χ; 21 twist-edge images at fixed x₁): max d_χ/H (everything but E_χ) | 0.489 | 0.246 | 0.100 | ∝ Λ⁻¹ (floor ∝ Λ⁻³, costs ∝ Λ⁻⁴) |
| Theorem W: min tolerated E_χ / S_E,t (+− transferred), binding state J = 26 → 52 → 132 at x₁ = 0.975 | 0.406 | 3.45 | 57 | ∝ Λ^{3.1} |
| Theorem W: same, plain envelope | 158 | 1.41·10³ | 2.26·10⁴ | ∝ Λ^{3.1} |
| Theorem W, 47 chart-B images with β exact: max d_χ/H; every state's d_χ/H ratio to 1000 | 0.489 | 0.245; 0.48–0.50 | 0.100; 0.19–0.22 | ∝ Λ⁻¹ |

### 0.3 Relation to GAP_SCALING (no duplication of claims)
* GAP_SCALING §3.1 certified LR∞ (β₀ = 0.37) on Λ ∈ [998.5, 10⁶] with 454 Arb Λ-balls and listed the Λ > 10⁶ tail as [D-sketch], "not done".
  **U1 closes that tail and is an independent second certificate on [Λ₀, 2604]** by a different method: the Λ⁴-scaling lemma (§1.2,
  with an envelope that makes it exact) for Λ ≥ 2604.3, and a 13-step monotone Λ-grid below (§1.3). The two certificates agree
  (β₀ at 2000: 0.2575 there per ball vs my anchor-based 0.2375 for Δ_g ≥ 2606 and 0.2725 from the R̂-envelope at exactly 2000).
* GAP_SCALING §2 proved Lemma C monotone in Λ; U4 is the same argument for Lemma C′ (a different bound, LARGESPIN_ANALYTIC §2).
* GAP_SCALING §3.1 says Theorem L′'s tolerances and ε_R "were not recomputed (hard-coded Δ_g)": §2 here does that; §4 does its plan item 4
  (Theorem W at the chart-B images).

### 0.4 What still needs per-gap work or stays a premise
1. **Premise content (unchanged, [M]/[A]):** P_L′ magnitudes (Λ_amp, K_sh, ε_R) for the graviton, (AG_U) beyond all orders, (EDGE-R) for the
   graviton, P1′, W, the +− transfer convention. The tables of §2/§4 show their *tolerances* improve with Δ_g at fixed (β, m/Λ); whether the
   graviton's premise *values* are Δ_g-uniform in AdS units is GAP_SCALING's H_mono-B [M].
2. **The functional's realisation at each Δ_g** (the published continuum k = 2 weight and its ε-DS) exists only at 300/1000. The *target*
   weights ψ_k(p) are Δ_g-independent by definition, so U1–U4 (which use only ψ_k) are unaffected; the realisation-dependent constant D(Δ_g)
   (monopole defect) was recomputed here from re-solved Abel weights at 2000/5000.
3. **Near-gap twist-edge strip at small J (chart-A seam)** — GAP_SCALING §5 (iv); not addressed here.
4. **Between and beyond the computed gaps** the §2/§4 tolerances are known at 1000, 2000, 5000 only; their monotone trends follow from the
   explicit Λ-dependence of the computed pieces (ψ₂(ν/Λ), D(Λ) with Λ²D → −1000, 1/m², the Λ⁻³ flat floor) [D-sketch], not from a proof.
5. U1′ improvement stops at β₀ = 0.10 (end of certified M₁ cells); an M₁ lower bound below 0.1 exists in closed form (§3.2) but was not
   wired into LR∞.

### 0.5 Which large-spin statements hold for all Δ_g ≥ 1000, and under which scaling hypotheses
| statement | for all Δ_g ≥ 1000 … | hypothesis needed beyond the Δ_g = 1000 premises |
|---|---|---|
| leading-Regge floor (LR∞), J ≥ J\*(m; Δ_g); F_LR > 0 at b̂ ≥ 0.2456Λ | **unconditionally** [P + certified] (U1, U2) | none |
| hybrid floor (flat + LR curvature) on b̂ ∈ [114, 2Λ] (all m ≥ Λ) and b̂ ≤ 119 (C′ samples, m ≥ 1.3Λ) | **unconditionally** [D + certified] / [P + DV] (U3, U4) | none |
| **Theorem L′** (M_Ω ⪰ 0 for J ≥ J\*(m; Δ_g)) | [D given P_L′(Δ_g)] | **H_L:** P_L′ at Δ_g with the *same* AdS-unit constants as at 1000 (Λ_amp, K_sh relative to the scalar-sized O(1/m²) structures; e_sc coefficients; AG_U on the rectangle of height (25 + 4 ln(Λ/Λ₀))/κ₀). Under H_L every tolerance at the chart-B images is ≥ its 1000 value (×4–25 at 2000/5000) [DV, sampled] |
| **Theorem W** (sampled near-gap states, J ≥ 26 at 1000 → their chart-B images) | [P given P0, P1′, W, AG_U, EDGE-R at Δ_g] + [DV] | **H_W = GAP_SCALING's H_mono-B for the edge:** the graviton's E_χ in units of its own flat-profile scale S_E,t, at fixed (β, x₁), grows slower than Λ^{3.1}. Then the Δ_g = 1000 statement transfers with margin ×8.5 (2000) and ×140 (5000) at the binding state |
| the twist-edge seam at small J, non-sampled J, J < 26 | not addressed here | GAP_SCALING §5 (iv) |

---

## 1. Lemma LR∞ uniform in Λ (`gl_common.py`, `gl1_lrinf_uniform.py`, `gl2_lrinf_allgaps.py`, `gl3_jstar.py`)

### 1.1 Λ-dependence of every ingredient
| ingredient | Λ-dependence in gap units (β, u, σ) | source |
|---|---|---|
| ψ_k(p), its derivatives at 0 and 1, ∫\|ψ_k\|, A₁,min = 327.149, A₂′ = 2047.09 | **none** (the functional is *defined* by ψ_k(ν/Λ); LRINF_TIGHT, prl/data README §1) | [P] |
| M₁(β) and its 830 certified cell lower bounds, region-(c) bound | **none** | [certified] |
| leading term A₁M₁/Λ² | exactly Λ⁻² | [P] |
| flat tail A₂′/(2b̂³) = A₂′/(2Λ³s³), s = 2 sinh(β/2) | Λ⁻³ (relative Λ⁻¹: it helps less at larger Δ_g) | [P] |
| T_a, T_b, T_c₂ | sums of positive coefficients × Λ^{−q}, q ≥ 4 | [P] |
| T_c₁ = ∫₀^{β/2} R(Λt)ρ_b dt | R(Λt) decreasing in Λ; Λ⁴R(Λt) is *not* monotone (R's small-x branch a₁/x² + ∫\|ψ\|) → fixed by the envelope R̂ (§1.2) | [P] |
| the realised k = 2 weight (continuum fit), ε-DS, D | per gap (not used by the lemma) | — |

So `r(β; Λ) = Λ²·N(β; Λ)/(A₁,min M₁(β))` with N = (√2/π)Σ_k w_k[T_a + T_b + T_c₁ + T_c₂] carries all Λ-dependence, and it falls like Λ⁻²
at fixed β.

### 1.2 Scaling lemma [P]
Let `V_k := max(s₁⁴R_∞(s₁), a₁s₁² + ∫|ψ| s₁⁴)` for any s₁ > 0 and `R̂(s) := min(R_∞(s), V_k/s⁴)`. Since `s⁴R_∞(s) = Σ_j c_j s^{3−j}` (j ≥ 3)
is non-increasing and `a₁s² + ∫|ψ|s⁴` increasing (a₁ > 0 for k = 2, 3), `V_k ≥ sup_s s⁴R(s)`; hence `R̂ ≥ R ≥ |r|` and `s⁴R̂(s)` is
non-increasing. With R̂ in T_c₁ (and in the tail majorant), **every piece of Λ⁴·N̂(β; Λ) is non-increasing in Λ at fixed β**:
`Λ⁴T_a ∝ Λ^{3−p}` (p ≥ 3), `Λ⁴T_b ∝ Λ^{4−p}` (p ≥ 4), `Λ⁴T_c₂ ∝ Λ^{3−j}` (j ≥ 3), `Λ⁴R̂(Λt) = t⁻⁴·(s⁴R̂)(Λt)`. Therefore, for Λ ≥ Λ_a,

    r(β; Λ) ≤ (Λ_a/Λ)²·r̂(β; Λ_a),      r̃(β; Λ) ≤ (Λ_a/Λ)·r̃̂(β; Λ_a),      Λ²F_LR ≥ min(A₁,min M₁, Λ_a²F̂_cert(β; Λ_a)).

(The r̃ bound: `(a + b/Λ_a)/(a + b/Λ) ≤ Λ/Λ_a`. The floor: `a − (Λ_a/Λ)²e + b/Λ` is concave in 1/Λ, minimal at an end.)
Certified constants (`gl1`): s₁ = 15.742 / 15.686, **V₂ ≤ 684 625.18, V₃ ≤ 789 635.86** (c₃ = |ψ‴(0)| = 2258 / 6705).

**Cost of the envelope at Λ₀** (`gl1`, cell values r̂ vs LRINF's r): ×1.39 at β = 0.37 (1.383 vs 0.997), ×1.30 at 0.5, ×1.15 at 1, ×1.0 beyond 5.
With R̂ alone the uniform thresholds would be 0.40 / 0.295, which is why §1.3 is needed between Λ₀ and the first anchor.

**Direct check** (`gl1` with argv 2000 5000): N̂ evaluated directly at Λ(2000), Λ(5000) on all 910 cells: **0 violations** of
`r(Λ) ≤ (Λ₀/Λ)² r̂₀`, worst ratio 1.0000 (attained at large β, where the j = 3 term, exactly ∝ Λ⁻⁴, dominates).

### 1.3 Λ ∈ [Λ₀, 2604.28]: monotone grid with LRINF's own bound [P + certified]
Without the envelope, every piece of N(β; Λ) is still **non-increasing in Λ** (T_c₁: both branches of R are decreasing in x). So for Λ in
[Λ_i, Λ_{i+1}]: `r ≤ Λ_{i+1}²N(Λ_i)/(A₁M₁)`, `r̃ ≤ Λ_{i+1}²N(Λ_i)/(A₁M₁ + A₂′/(2Λ_{i+1}s³))`. Steps q_i = Λ_{i+1}/Λ_i chosen adaptively
(exact rationals; Λ_0 a rational lower bound of Λ₀) under the two requirements r < 1 on β ≥ 0.37, r̃ < 1 on β ≥ 0.245 (cells + tail):

| i | Λ_i | q_i | sup r (β ≥ 0.37) on the step | sup r̃ (β ≥ 0.245) |
|---|---|---|---|---|
| 0 | 998.500 | 1.0015 | 0.9998 | 0.9881 |
| 1 | 999.998 | 1.0030 | 0.9995 | 0.9900 |
| 2 | 1002.998 | 1.0050 | 0.9967 | 0.9910 |
| 3 | 1008.013 | 1.0100 | 0.9954 | 0.9973 |
| 5 | 1033.364 | 1.0200 | 0.9601 | 0.9867 |
| 8 | 1139.935 | 1.0750 | 0.8567 | 0.9833 |
| 10 | 1347.973 | 1.1500 | 0.6776 | 0.9187 |
| 12 | 1860.203 | 1.4000 | 0.5005 | 0.9564 |

13 steps reach Λ_n = 2604.284 (Δ_g = 2605.78). **Anchor Λ_n with R̂** (§1.2) closes [Λ_n, ∞): sup_{β≥0.37} r̂ = 0.147, sup_{β≥0.245} r̃̂ = 0.294,
tail (β ≥ 64) 3.0·10⁻⁸. Hence **U1 [P + certified]** (`gl2_allgaps.txt`; ≈ 1 CPU-min).

### 1.4 Δ_g-improved thresholds and the uniform floor
Anchors (R̂ at Λ_a, valid for every Λ ≥ Λ_a):

| anchor | Λ_a | sup r̂ (β ≥ 0.37) | sup r̃̂ (β ≥ 0.245) | β₀(r < 1) | β₀(F_LR > 0) |
|---|---|---|---|---|---|
| Λ_n (Δ_g ≥ 2606) | 2604.3 | 0.147 | 0.294 | 0.2375 | 0.1325 |
| Δ_g ≥ 5000 | 4998.5 | 0.035 | 0.099 | 0.1675 | 0.10 |
| Δ_g ≥ 10⁴ | 9998.5 | 0.0083 | 0.030 | 0.1175 | 0.10 |
| Δ_g ≥ 2·10⁴ | 19998.5 | 0.0020 | 0.0084 | 0.10 | 0.10 |

**Uniform floor U2** (min over all segments and the closing anchor, per cell): f_U ≥ 0.217 at [0.245, 0.247] (1.6 % of A₁M₁), 3.08 at 0.25,
11.06 at 0.30, and **f_U ≥ 0.99 A₁,min M₁ for every cell from β = 0.34** (≥ 0.9 from 0.315, ≥ 0.5 from 0.265). The binding Λ is ∞ there
(the flat tail vanishes, the remainder vanishes): in gap units the floor tends to exactly A₁,min M₁(β). At Λ₀ the numerical floor is
larger (Λ₀²F_num = 34.0 at β = 0.37, of which ≈ 20 is the flat tail): **the flat tail's share of the floor decays like Λ⁻¹**, which is
why the numerical Q_tol below drops ×1.4–2.5 to its Λ-independent limit while the certified one does not move.

### 1.5 J\*(m; Δ_g) [certified] (`gl3_jstar.txt`, Arb integer decisions)
`J ≥ J*(m) = ⌈(s₀ m − 1)/2⌉`, m in AdS units, s₀ = 0.37211 (r < 1) or 0.24561 (F_LR > 0), for **every** Δ_g ≥ 1000; in gap units
`(J+½)/Λ ≥ (s₀/2)(m/Λ)`. Twist edge: first J with `(2J+1)² ≥ s₀²(Δ_g−2)(Δ_g+2J−1)`.

| Δ_g | m_gap | 1.1 m_gap | 1.3 m_gap | 1.5Λ | 2Λ | 5Λ | J\*_edge/Λ |
|---|---|---|---|---|---|---|---|
| 1000 | 223 | 250 | 307 | 279 | 372 | 929 | 0.2233 |
| 2000 | 447 | 501 | 614 | 558 | 744 | 1859 | 0.2237 |
| 5000 | 1119 | 1253 | 1536 | 1395 | 1860 | 4650 | 0.2239 |
| 10⁴ | 2238 | 2507 | 3073 | 2790 | 3721 | 9301 | 0.2238 |
| → ∞ | | | | | | | **x\* = (s₀² + s₀√(s₀²+4))/4 = 0.22387** |

(F_LR > 0 form: 139 / 277 / 694 / 1388 at the twist edge, J\*/Λ → 0.13881.) With U1′, e.g. Δ_g ≥ 5000: J\* = 456 at the 5000 twist edge
(s₀ = 0.16770), 629 at 1.5Λ.

---

## 2. Theorem L′ uniform in Δ_g: tolerances and the remainder (`gl4_tol.py`, `gl6_p4w.py`)

### 2.1 Scaling [D] and the monopole constant
At fixed gap-unit state (β, μ):
* floor F ∝ Λ⁻² (U2; its numerical value carries an extra flat-tail part ∝ Λ⁻³);
* the O(1/m²) correction costs are `2[∫ψ₂(p) f(Λp) dp + D Q]/m²` with f, Q functions of (β, ν) in AdS units; the ν-integral is dominated by
  ν = O(1), where ψ₂(ν/Λ) ≈ a₁ν/Λ, so it is ∝ Λ⁻²; **D·Λ² → −1000** (−995.34 / −999.08 / −1001.32 at 1000 / 2000 / 5000, from Abel weights
  re-solved with `deshift/d_abel` on dν = 0.04; D(1000) reproduces PW4's −9.98332·10⁻⁴ to 10⁻¹²); and 1/m² = 1/(μΛ)².
* Hence every tolerated *coefficient* (Λ_amp, K_sh, K_tol) **grows like Λ²** and Q_tol = F/(4|D|) is Λ-independent. This presumes the premise
  is stated in AdS units with Δ_g-independent coefficients (κ^sc of CHMRS (3.20) is Δ_g-independent; the graviton's is [M]).

### 2.2 Tolerances (`gl4`; pw5/lrt5 definitions; "unif" = certified Δ_g-uniform floor f_U/Λ²) [DV]
Chart-B images m = μΛ, J = round(μΛ sinh(β/2) − ½):

| state (1000) | Δ_g | J | Q_tol num / unif | Λ_amp mono num / unif | Λ_amp vector R1′ (unif) | K_tol sh±1 (unif) | K_tol C2 (unif) |
|---|---|---|---|---|---|---|---|
| gap J = 223 (β 0.37, 1.20Λ) = J\* image | 1000 | 223 | 8.53e-3 / 3.35e-3 | 1819 / 714 | 126 | 6.5e3 | 154 |
| | 2000 | 447 | 6.00e-3 / 3.34e-3 | 5130 / 2851 | 482 | 2.6e4 | 509 |
| | 5000 | 1118 | 4.51e-3 / 3.33e-3 | 24088 / 17804 | (‡) | 1.6e5 | (‡) |
| gap J = 300 (β 0.47) | 1000 / 2000 / 5000 | 300 / 601 / 1504 | 3.33e-3 / 3.31e-3 / 3.31e-3 | 754 / 3008 / 18776 | 250 / 977 / ‡ | 8.3e3 / 3.3e4 / 2.1e5 | 252 / 926 / ‡ |
| gap J = 500 (β 0.69) | 1000 / 2000 / 5000 | 500 / 1001 / 2505 | 3.24e-3 / 3.23e-3 / 3.22e-3 | 818 / 3263 / 20365 | 712 / 2307 / ‡ | 9.5e3 / 3.8e4 / 2.4e5 | 951 / 1.3e3 / ‡ |
| gap J = 857 (β 1.0) | 1000 / 2000 / 5000 | 857 / 1716 / 4292 | 3.07e-3 / 3.06e-3 / 3.06e-3 | 853 / 3405 / 21252 | 1048 / 3539 / ‡ | 9.5e3 / 3.8e4 / 2.4e5 | 8.0e3 / 3.5e3 / ‡ |
| gap J = 10⁴ (β 3.05) | 1000 / 2000 | 10⁴ / 20016 | 1.38e-3 / 1.38e-3 | 435 / 1734 | 447 / 1735 | 5.4e3 / 2.2e4 | 5.6e3 / 7.3e3 |
| 1.5Λ, J = 279 (β 0.37) | 1000 / 2000 / 5000 | 279 / 559 / 1399 | 3.35e-3 / 3.34e-3 / 3.33e-3 | 1111 / 4434 / 27675 | 195 / 751 / ‡ | 1.0e4 / 4.0e4 / 2.5e5 | 220 / 770 / ‡ |
| 5Λ, J = 929 (β 0.37) | 1000 / 2000 / 5000 | 929 / 1860 / 4653 | 3.35e-3 / 3.34e-3 / 3.33e-3 | 12349 / 49283 / 307603 | 2186 / 8275 / ‡ | 1.1e5 / 4.5e5 / 2.8e6 | 2.7e3 / 8.4e3 / ‡ |
| twist edge J\*(Δ_g) | 1000 / 2000 / 5000 | 223 / 447 / 1119 | 3.35e-3 / 3.34e-3 / 3.33e-3 | 714 / 2853 / 17813 | 126 / 482 / ‡ | 6.5e3 / 2.6e4 / 1.6e5 | 154 / 508 / ‡ |

(‡ The C2 structure sinh β ν² ∂_βP, and through it the vector R1′ entry, is **not numerically converged at Δ_g = 5000**: at the J\* image
K_tol C2 = 2734 / 317 with 1× / 2× Mehler nodes, while at 1000 and 2000 it is stable to all digits under 1×/2×/4× (`gl4b_c2conv.py`,
`cache/gl4_c2conv.log`). The likely cause is the finite-difference ∂_β (step 10⁻⁴β) at ν ~ Λ, where the ν²-weighted terms cancel; not fixed.
Also: lrt5/pw5's node count 600 + 600β under-resolves C2 at large Λβ already at 1000 (J = 857: 1.14·10⁴ there vs 8.0·10³ here). C2 is
not the binding class for Λ_amp ≲ 100. The monopole and sh±1 columns use no β-derivative and are converged.)

Reading: Q_tol (certified floor) is **Δ_g-independent** (3.35e-3 at the J\* image); the tolerated Λ_amp grows ×4.0 (monopole) / ×3.8
(vector) at 2000 and ×25 (monopole) at 5000; K_tol sh±1 ×4 / ×25. The numerical-floor Q_tol *decreases* (8.5e-3 → 6.0e-3 → 4.5e-3) toward the certified
value, because the flat tail of F_LR is a 1/Λ effect (§1.4).

### 2.3 The remainder ε_R: Lemma P4W with the uniform rate (AG_U) on the 16 L′ images (`gl6`) [DV]
AG_PROOF §6.2's 'U' family (θ = 0.97), re-typed; Λ enters only through `p4w_common.L` (re-set per gap). Λ_amp = 3, +− transferred scalar
envelope e_sc with its Δ_g = 1000 AdS-unit coefficients [M]. At Δ_g = 1000 it reproduces AG_PROOF exactly (**max d = 0.133, edge tolerance
18.9** with LRINF's F_cert). Budget here: the certified uniform floor f_U/Λ².

| Δ_g | max d/(F_U/2) (16 images; β > 3.1 only at 1000, β > 2 excluded at 5000 — Mehler cost) | at the J\* image: window / boundary / edge+top | min edge tolerance (multiple of the transferred envelope) |
|---|---|---|---|
| 1000 | 0.190 (J = 223) | 1.7e-2 / 1.0e-1 / 7.0e-2 | 12.6 (J = 223) |
| 2000 | 0.035 | 4.2e-3 / 2.6e-2 / 4.3e-3 | 228 |
| 5000 | 0.0050 | 6.7e-4 / 4.2e-3 / 1.1e-4 | 9.05·10³ |

* **Rectangle height.** With AG_PROOF's fixed height T₀ = 25/(0.97Λβ) the top side is e⁻²⁵-suppressed but *not* Λ-suppressed, so against the
  Λ⁻² floor it grows like Λ²: with T₀ the edge+top share at β = 1 is 1.1e-4 / 3.6e-4 / 2.3e-3 (1000/2000/5000; `cache/gl6_fixedT_partial.log`).
  This is an artefact of the fixed height. The table uses `T = (25 + 4 ln(Λ/Λ₀))/(0.97Λβ)` (27.8 T₀/25 at 2000, 31.4 T₀/25 at 5000), which
  requires (AG_U) on a correspondingly taller rectangle — a slightly stronger form of the same premise.
* Hence **Theorem L′'s P4 remainder ε_R ≤ 0.19 of the certified half-floor at Δ_g = 1000 (uniform floor; 0.133 with the Δ_g = 1000 F_cert),
  0.035 at 2000, 0.005 at 5000**, at the images of all 16 sampled L′ states, under the transported premises.

---

## 3. Lemma C′ and ROUTEA uniform in Δ_g

### 3.1 Lemma C′ is monotone in Λ [P]
LARGESPIN_ANALYTIC §2: `|curv_k| ≤ |w_k|(2/π)∫₀^{π/2}[E₀(Λs)s²/8 + Λσ³/24·E₁(Λs)]dφ`, `Λσ = b̂ sin φ`, s = 2 asinh(σ/2) (code: the lower bound
σ(1 − σ²/24) for s in the arguments). At fixed (J, μ, h), b̂ = (2J+1)/μ is fixed, so as Λ grows: Λσ is fixed, σ and s decrease (s² ↓),
Λs = b̂ sin φ·(asinh y/y)|_{y = b̂ sin φ/2Λ} increases (so E_j(Λs), non-increasing envelopes, decrease), and Λσ³ = (b̂ sin φ)³/Λ² decreases.
Every term is non-increasing in Λ; N_flat (and u) are Λ-independent at fixed (J, μ, h). Hence **rig′'s numerator N_flat/5 − 2B₂ − uB₃ is
non-decreasing in Λ**, and every sampled state with rig′ > 0 at Δ_g = 1000 has it at every Δ_g ≥ 1000 (the state stays above the twist edge:
m_gap(J)/Λ decreases with Δ_g). [P] for the monotonicity; the Δ_g = 1000 inputs keep their [P + DV] grade and their sampling in b̂.
Since B_k ∝ Λ⁻² at fixed b̂ while N_flat is fixed, the reach in b̂ grows (≈ Λ^{2/3} from the flat 1/b̂³ tail) [D-sketch].

### 3.2 ROUTEA's window, uniformly in Λ (`gl7_routea_uniform.py`) [D + certified]
ROUTEA_LARGE_B §4.2/4.4: `|curv − A₁M₁/Λ²| ≤ R = [C_P + C_Q(1 − s₁²/b̂²)^{−3/2}/b̂² + C_sm/Λ² + C_osc2]/b̂` for s₁ ≤ b̂ ≤ 2Λ. C_sm and Λ²C_osc2
are Λ-independent by construction. The per-gap region-I constant C_P is not: **re-running ROUTEA's code at other gaps (`gl5_routea_spot.py`, a
verbatim copy with Δ_g as input) gives Λ²C_P(s₁ = 80) = 291 / 340 / 676 / 1878 at 1000 / 2000 / 5000 / 10⁴, and the certified window start
degrades 91.17 / 91.17 / 96.72 / 159.86** (its s < 4 crude cells and Lipschitz pads are O(1), not O(Λ⁻²)).

**New [D]:** with G = w r(τ) − r(s) (r = Φ + a₁/x², w = (1 + s²/4Λ²)^{−½}, τ = 2Λ asinh(s/2Λ)) and ε = Λ⁻²,

    Λ²G = g₁(s) + ε h(s),   g₁ = −(s³/24) r′(s) − (s²/8) r(s) = −d/ds[(s³/24) r(s)],
    |h| ≤ (3s⁵/640)|r′(s)| + (3s⁴/128)|r(s)| + (s⁶/1152) sup_{[5s/6,s]}|r″| + (s⁵/192) sup_{[5s/6,s]}|r′|

(alternating Taylor bounds of asinh x and (1+y)^{−½} on [0, 1] — b̂ ≤ 2Λ ensures x, y ≤ 1 — and τ ≥ 5s/6). Since s³r(s) → 0 at s → 0,
**∫₀^{s₁}g₁ = −s₁³r(s₁)/24 exactly**, so for every Λ ≥ Λ₀: `|∫₀^{s₁}Λ²G| ≤ s₁³R₀(s₁)/24 + ε₀∫h̄`, `∫|Λ²G|s² ≤ ∫(ḡ + ε₀h̄)s²`, with
R_n(x) ≥ |r^{(n)}(x)| the min of the exact finite integration-by-parts sum (absolute values) and the crude `∫|ψ| + c_n a₁/x^{2+n}`.
The positivity condition is checked on b̂-boxes (ratio 1.03) with `M₁` bounded below over **all** β ≤ 2 asinh(b̂_box/2Λ₀) (β ≤ 2 asinh 1
since b̂ ≤ 2Λ), using the lrt2 cells and, on (0, 0.1], the [P] bound
`M₁(β) ≥ (1/24 − 7β²/1920)/cosh(β/2) − 2e^{−2πν₁}(ν₁/2π + 1/4π²)`, ν₁ = √2/β (from cos x ≥ 1 − x²/2, |P| ≤ 1; ≥ 0.041578 there).

| m-range (u_max) | A₁,min | s₁ | Λ²C_P | Λ²C_Q | C_sm | Λ²C_osc2 | **window, every Δ_g ≥ 1000** |
|---|---|---|---|---|---|---|---|
| m ≥ 1.5Λ (4/9) | 355.28 | 80 | 17.83 | 3.66·10⁶ | 43.20 | 363.6 | **[107.51, 2Λ]** |
| m ≥ 1.3Λ (1/1.69) | 347.83 | 80 | 19.10 | 3.89·10⁶ | 49.42 | 384.4 | **[110.74, 2Λ]** |
| m ≥ Λ (1) | 327.15 | 80 | 22.62 | 4.55·10⁶ | 66.69 | 442.2 | **[114.06, 2Λ]** |

(Λ²C_P = 17.8 replaces ROUTEA's 291: the leading signed integral is the boundary term s₁³|r(s₁)|/24 = 1.49 / 4.08 for k = 2 / 3; the rest
is ε₀∫h̄ ≈ 9.5 each.) Margin at the window start (m ≥ 1.5Λ): A₁,min min M₁ = 14.60 vs U = 13.82. Tail b̂ ≥ 2Λ₀: 10.8 vs 0.21.
Grade [D + certified]: the inequalities are ROUTEA's [D] plus the expansion above [D]; all constants and decisions in Arb (`gl7_u*.log`).
Hybrid positivity then uses only the certified NF02 flat floor ≥ N/5 on every flat state (CERT_CLOSE §1), which is Λ-independent; the
old "+−, J ≥ 400, m < 1.633Λ" flat corner of ROUTEA §0 is covered by that certificate.

### 3.3 Coverage map of the hybrid / leading-Regge floor, every Δ_g ≥ 1000
| m | b̂ ≤ 107–114 | [107.5 / 110.7 / 114.1, 2Λ] | b̂ ≥ 0.2456Λ (F_LR) | b̂ > 2Λ |
|---|---|---|---|---|
| ≥ 1.5Λ | C′ (U4, b̂ ≤ 119) | **ROUTEA-U (U3)** | LR∞ (U1) | LR∞ (F_LR; hybrid-vs-F_LR [NE]) |
| [1.3, 1.5)Λ | C′ (U4) | **ROUTEA-U: closes LARGESPIN_ANALYTIC's (119, 245) hole** | LR∞ | LR∞ |
| [Λ, 1.3)Λ | Lemma C only where GAP_SCALING §2.2 finds it (J ≲ 13–35) | **ROUTEA-U [114, 2Λ]** | LR∞ | LR∞ |

---

## 4. Theorem W at the chart-B images, Δ_g = 2000 and 5000 (`gl8_thmW_chartB.py`) [DV]

This is GAP_SCALING's plan item 4 (§6.1), done without new symbol computations.

### 4.1 Why chart B is cheap: the computed O(1/m²) corrections are functions of β alone [D]
P2_CONTRACT's per-state O(1/m²) corrections (block side `pc_run.run_state`, kernel O(r) term `pc_ksh_run`) are computed from the exact
symbols at q = e^{β/2} and fitted as F(μ) = a(μ)P_μ + b(μ)P′_μ; the fit coefficients x depend **only on β and the parity of J** (code
inspection: neither J nor m nor Λ enters before `costs(x, β, m)`). Per unit leading image and in c₂ units they are therefore exactly chart-B
invariant. At the image state only the costs `2[∫ψ₂(ν/Λ′)F/c dν/Λ′ + D(Λ′)F(i/2)]/m′²`, the floors H′ (`ls_core.hyb_fast` at Δ′_g), the
window/boundary/edge contour terms (Λ′ enters through p = ν/Λ′, κ₀ = 0.97Λ′β, e_sc) and the edge scale S_E (Lemma F flat wave at (J′, m′))
are recomputed. Δ_g enters the borrowed modules only through the module constants `lsw_common.{Dg, L, NU, X, PS, DN}`,
`p4w_common.{L, D1000}`, `pc_run.D1000`, `ag_flat.L`, `ag_chi.L`, which `gl8.set_gap` re-sets in-process (the request to "make Δ_g a
parameter" in these modules is met this way; `ag_exact` is not needed — it only produced the quoted exact-scalar edge data).
**Check:** at Δ_g = 1000 the driver reproduces AG_PROOF's ag_chi rows exactly (J = 26/64/150/344: d_χ/H, ε_E tolerances
plain/transferred, and P4_WINDOW's local-rate tolerances 114 / 0.287), `cache/gl8_val.log`.

Two image families (both keep the premise forms of Theorem W with their Δ_g = 1000 AdS-unit constants [M]):
* **(B) chart-B images, β exact:** J′ the smallest integer of J's parity with J′ + ½ ≥ (J+½)Λ′/Λ₀ and above the Δ′ twist edge; then
  m′/Λ′ = (m/Λ₀)(1 + δ), 0 ≤ δ ≤ 2.75 % (largest at J = 26). Twist-edge states therefore move *off* the edge (x₁: 0.975 → 0.949 at J = 26).
* **(E) twist-edge images:** each Δ_g = 1000 twist-edge state (smallest sampled m per J; 21 states) is mapped to the **Δ′ twist edge** at
  J′ ≈ (J+½)Λ′/Λ₀ − ½ (same parity), so x₁ and the +− transfer factor (1 − x₁²)⁻² are unchanged; β′ differs from β by |δβ/β| ≤ 1.0 %
  (2000) / 0.16 % (5000), and the β-fits are re-used at β′ (an approximation of that relative size in the computed O(1/m²) costs, which are
  ≤ 0.1 H).

### 4.2 Results [DV]
**(E) twist-edge images — this is where EDGE-R binds** (`gl8_W_{2000,5000}_edge.*`):

| J (1000) | x₁ | (1−x₁²)⁻² | J′ (2000 / 5000) | d_χ/H 1000 / 2000 / 5000 | **ε_E tol / S_E,t (transferred)** 1000 / 2000 / 5000 | plain ε_E tol / S_E 1000 / 2000 / 5000 |
|---|---|---|---|---|---|---|
| **26** | 0.975 | 409 | 52 / 132 | 0.079 / 0.042 / 0.014 | **0.406 / 3.45 / 57** | 166 / 1.41e3 / 2.26e4 |
| 30 | 0.971 | 312 | 60 / 152 | 0.082 / 0.043 / 0.017 | 0.63 / 5.2 / 86 | 196 / 1.62e3 / 2.61e4 |
| 34 | 0.968 | 246 | 68 / 172 | 0.092 / 0.046 / 0.019 | 0.89 / 7.5 / 123 | 218 / 1.84e3 / 2.96e4 |
| 38 | 0.964 | 200 | 76 / 192 | 0.101 / 0.049 / 0.021 | 1.22 / 10.2 / 169 | 243 / 2.05e3 / 3.32e4 |
| 50 | 0.953 | 121 | 100 / 252 | 0.131 / 0.063 / 0.027 | 2.57 / 22.3 / 369 | 311 / 2.69e3 / 4.40e4 |
| 100 | 0.913 | 36 | 200 / 502 | 0.234 / 0.118 / 0.048 | 16.2 / 150 / 2.56e3 | 582 / 5.40e3 / 9.17e4 |
| 344 | 0.770 | 6.0 | 690 / 1724 | **0.489 / 0.246 / 0.100** | 280 / 3.34e3 / 6.23e4 | 1.69e3 / 2.0e4 / 3.74e5 |
| all 21 | | | | max **0.489 / 0.246 / 0.100** | min **0.406 / 3.45 / 57** | min 158 / 1.41e3 / 2.26e4 |

**(B) chart-B images at 2000, β exact** (47 states, `gl8_W_2000.*`): max d_χ/H = 0.245 (J = 344 → 690); every state's d_χ/H is 0.48–0.50 of
its 1000 value; min ε_E tol / S_E,t = 14.9 at (54, 1.054Λ) (x₁ = 0.949, transfer 100; the 1000 state at that x₁ would have ≈ 5), plain ≥ 1.5·10³.
(5000: `gl8_W_5000.*`, see §4.4.)

### 4.3 Laws, and what they say about EDGE-R
* **Floors:** H′ ≈ H·(Λ₀/Λ′)³ at fixed (β, x₁) (H′/H·(Λ′/Λ₀)² = 0.48–0.51 at 2000, 0.19–0.20 at 5000): near the gap the hybrid floor is the
  flat partial wave, ∝ b̂⁻³ with b̂ = Λb.
* **Everything but the edge integral:** d_χ/H ∝ Λ⁻¹ (ratios 0.50 and 0.20), i.e. the computed O(1/m²) costs + window + boundary + AG_U
  Cauchy terms fall like Λ⁻⁴ against a Λ⁻³ floor. The worst state (J = 344, β = 0.525) goes 0.489 → 0.246 → 0.100.
* **EDGE-R:** at fixed x₁ the transfer factor is fixed, and the tolerated edge integral in units of its own flat-profile scale grows like
  **Λ^{3.1}** (binding state: 0.406 → 3.45 → 57; Λ^{3.09} and Λ^{3.06}). So "a scalar-like +− edge (≤ 0.30 S_E,t, the exact scalar's
  worst at 1000)" passes at the binding state **with margin 11.5 at 2000 and 190 at 5000** (1.35 at 1000). This is the chart-B statement
  GAP_SCALING §3.3 asked for: in chart B the EDGE-R margin improves; at fixed J (chart A) the transfer factor grows as (Δ_g/J)² and it
  does not — so the near-gap strip must be read in chart B, as GAP_SCALING concluded.
* **What this does not settle:** the graviton's E_χ in S_E,t units at any gap (it is [M] at 1000 already; EDGE_R.md). The tables say: *if*
  the graviton's edge integral, in its own flat-profile units, is Δ_g-independent at fixed (β, x₁) (H_mono-B for the edge) then EDGE-R's
  margin grows like Λ³; if it is scalar-like it already passes at 2000 at every sampled image.
* **P4W with AG′'s local rate** (P4_WINDOW's form; reported for comparison, column "Λ_AG tol" in the logs): improves at small J
  (2.54 → 43.6 at the binding state) but *degrades* at J ≥ 150 at 5000 (147 → 12.1 at J = 344): its fixed rectangle height
  25/κ₀ leaves an e⁻²⁵ top term that is not Λ-suppressed, against a Λ⁻³ floor (the artefact of §2.3). The operative P4W-χ form (ag_chi,
  whose height grid T ≤ T₀ is searched) does not show it at these Λ.

### 4.4 Chart-B images at 5000
`gl8_W_5000.*` (47 states, β exact, m′/Λ′ shifted by ≤ 1.4 %): max d_χ/H = 0.100; every state's d_χ/H is 0.19–0.22 of its 1000
value; H′/H·(Λ′/Λ₀)² = 0.19–0.20; min ε_E tol/S_E,t = 133 at (134, 1.040Λ) (x₁ = 0.962, transfer 176), plain ≥ 2.3·10⁴. Consistent with
the twist-edge family (E).

---

## 5. Caveats
1. **Certified** = U1, U1′, U2, J\* tables, U3's constants. The lemma inequalities are LRINF_TIGHT's / ROUTEA's [P]/[D] plus §1.2/§3.2 here.
   M₁ cell values are read from lrt2's float-stored lower endpoints and shaved by a relative 2⁻⁴⁰ (lrt3 used them unshaved).
2. **[DV]** = §2 and §4: float quadrature, sampled states, premises transported with Δ_g = 1000 AdS-unit constants [M].
3. Chart-B images at small J are not exact images: β is kept exact (the stored O(1/m²) fits are at that β) and J′ must be an integer of
   J's parity above the Δ′ twist edge, so m′/Λ′ = m/Λ₀(1 + δ), 0 ≤ δ ≤ 2/(2J+1) (§4).
4. U1′/U2 use the end of the certified M₁ cells (β = 0.1); the improved thresholds at Δ_g ≥ 2·10⁴ are cell-limited, not bound-limited.
5. LR∞'s floor is F_LR; the hybrid-vs-F_LR identification on b̂ ≥ 2Λ stays [NE] (LRINF_TIGHT §4), as at 1000.
6. §4 family (E) re-uses the Δ_g = 1000 O(1/m²) fits at β′ = β(1 + O(1/J)) (≤ 1.0 % at 2000, ≤ 0.16 % at 5000): an approximation of that
   relative size in costs that are ≤ 0.1 H.
7. §2.3's Λ-growing rectangle height asks (AG_U) on height (25 + 4 ln(Λ/Λ₀))/κ₀; §4 uses ag_chi's own height grid (T ≤ T₀) unchanged.
8. Compute used: ≈ 6 CPU-h (incl. two aborted runs: a chart-B rounding rule that put images below the twist edge, and an under-resolved
   Mehler count in gl4; both fixed before the reported runs).

## 6. Reproduce
```
cd ads4/code/gaplarge && bash run_all.sh      # ≈ 5 CPU-h, ≤ 2 processes
```

| file | content | grade |
|---|---|---|
| `gl_common.py` | Λ-parametrised LR∞ bound pieces; envelope R̂; I-tables; tail majorant; M₁ cells | — |
| `gl1_lrinf_uniform.py` | Λ⁴-scaling form at Λ₀; envelope cost; direct checks at 2000/5000 | [certified] |
| `gl2_lrinf_allgaps.py` | U1, U1′, U2 (monotone Λ-grid + anchors) | [P + certified] |
| `gl3_jstar.py` | J\*(m; Δ_g) tables | [certified] |
| `gl4_tol.py` | Theorem L′ tolerances, D(Δ_g) (`cache/gl4.log`: 1000, 2000; `gl4_5000.log`) | [DV] |
| `gl4b_c2conv.py` | Mehler-node convergence check of the C2/vector entries | [DV] |
| `gl5_routea_spot.py` | ROUTEA per gap (verbatim copy, Δ_g input) | [certified] per gap |
| `gl6_p4w.py` | Lemma P4W (AG_U) on L′ images (`cache/gl6_final.log`; fixed-height variant `gl6_fixedT_partial.log`) | [DV] |
| `gl7_routea_uniform.py` | U3 (`GL7_UMAX` = 4/9, 100/169, 1) | [D + certified] |
| `gl8_thmW_chartB.py` | Theorem W (P4W-χ) at chart-B images and twist-edge images | [DV] |
