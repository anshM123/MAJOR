# UNCOND_INVENTORY: can per-state positivity of Ω^EE be made unconditional for all Δ_g ≥ Δ*?

*25 Nov 2026. Research-planning note. Unreviewed, not committed. Reading and analysis only: no computation was run, and no existing
file was edited. Sources: the update boxes of `UNIFIED_THEOREM.md` (to K_B2A4) and `UNIFIED_THEOREM_v3.md` §2.5, §4.7, §5; GAP_SCALING,
GAP_LARGESPIN, GAP_UNIFORM, GAP_C3, GAP_MORSE, GAP_NEARGAP, GAP_SEAM, GAP_SEAM2, LRINF_TIGHT, AG_PROOF, P4_WINDOW, P2_CONTRACT,
PW_IMAGES, LARGESPIN_W, LARGESPIN_L, EDGE_R, EDGE_R2, BATCH5_FOLD, NEARGAP_FILL, CERT_CLOSE, POS_ERRORS, RD_CLOSE, PCYC_CLOSE,
PARITY_ODD_CLOSE, SPINNING_MELLIN, CERTIFIED_NEARGAP_A, NOFORWARD_FUNCTIONAL, K_B2A4, C1_FINAL; `prl/supplement.tex` Tables
cover, coverLS, allgapU, BF, delta1 and the transfer subsection.*

Grades as in UNIFIED_THEOREM. "Unconditional" below means: given A1 (or A1′) and A2, the definition of the functional, and nothing else.
That is, no [M], [A] or [D-sketch] link and no numerical premise; every inequality is [P] or [certified].

Effort is given in **agent-hours (ah)**, calibrated on the record. GAP_SCALING estimated "3–5 days" for the transfer programme; the six
GAP_* deliverables that carried it out landed in about 9 wall-hours with two agents. One deliverable is about 3–6 ah. Risk is the
probability that the item **fails to close as stated**: L ≤ 20 %, M 20–50 %, H 50–80 %, VH > 80 %.

---

## 0. Verdict

1. **An unconditional theorem "Ω^EE ⪰ 0 on every heavy state for all Δ_g ≥ Δ\*" is not reachable with the methods in the record, at
   any finite Δ\*, on the paper's timescale.** The obstacle is not compute. It is three missing pieces of analysis, all about the one
   object that carries every non-certified grade: the **beyond-leading-Regge correction δ** of a heavy block's action. Each of the three
   is needed at every gap.
   * **(O1) A contour theorem (CT).** At bulk ξ = ν/Δ, every bound on δ needs the physical chain D′ to be homologous to the crystal
     thimble plus pieces that are lower by a margin. That is what GAP_C3/GAP_MORSE call TAIL and what GAP_UNIFORM calls Lemma R
     [D-sketch]. Prop. D′, the only real-contour bound, works only for ξ < 0.107 (GAP_UNIFORM §4.6). Every Laplace certificate stands
     on CT, and so does every evaluation of δ at ν ~ Λ.
   * **(O2) Large-b̂ cancellation (AG_U beyond all orders).** Where the state's normalisation N_s ∝ b̂⁻³ comes from cancellation, no
     per-ν error bar can be propagated sign-free (POS_ERRORS §2.2). The ν-integral has to be deformed into the complex plane. That needs
     the travelling-wave split at a uniform rate: (AG_U) [A], with no proof route except certified complex-ξ steepest descent, which
     itself needs CT.
   * **(O3) Uniformity in J across J ~ √Δ_g.** Chart A (fixed J; the Laplace expansion parameter is J²/m) and chart B (fixed β; the
     parameter is roughly m/J²) break down in the same crossover band. The supplement already records a band with no transfer route,
     33 ≤ J < J_W ∝ Λ, which grows with the gap. No certified method covers it.

   A fourth, removable obstacle:
   * **(O4) The fold at the light-cone edge.** For twist-edge states, ξ reaches 1 − (J + 3/2)/Δ_g at every gap, which is inside the Airy
     region. It can be removed by a definitional change to the functional (§4, M8b).
2. **What is reachable.**
   * **T0 (unconditional, Δ\* = 1000).** The leading-Regge (hybrid) action of Ω^EE is ≥ c·N_s > 0 on **every** heavy state at
     **every** Δ_g ≥ 1000 [certified]. Almost everything is already in hand: the flat floor N/5, Lemma C/C′, ROUTEA-U and LR∞ uniform in
     the gap. What is missing is one hole (m < 1.3Λ, b̂ ≲ 114), the sampled-mass part of C′, and two cheap edges. **≈ 10–16 ah, low
     risk.** This is a genuine unconditional theorem, but about the bulk-point (leading-Regge) action, not the full action.
   * **T1 (one named mathematical premise, Δ\* ≈ 10⁶).** For a finite set of spins (J ≤ 8, extendable), every heavy state at
     every mass is positive for Δ_g ≳ 10⁶. The only premise is CT. The other GAP_MORSE premises (Π ≤ 3, ξ-continuity, the K₂ part
     bounded, graviton kinds at J ≥ 4 and odd J) and the two windows can all be closed. **≈ 35–60 ah, medium-high risk.** It gives no
     statement for "every heavy state".
   * **T∞ (full theorem, conditional on CT only, Δ\* ≈ 10⁶–10⁷).** A uniform certified steepest-descent programme in (complex ξ,
     j = J/Δ), which removes O2 and O3 given CT. **≈ 120–200 ah, success probability ≲ 20–25 %.** It is a research programme, not a
     closing task.
3. **Infeasible in principle with available information (§5).**
   * Any certified statement about δ at Δ_g = 1000, and generally below about 2·10⁴. At Δ = 1000 the optimal-truncation floor of the
     asymptotic series is comparable to the needed tolerance, and no non-asymptotic representation of the spinning action exists
     (SPINNING_MELLIN §5).
   * Hence the range 1000 ≤ Δ_g < Δ₁ can never become better than premise B_F.
   * Below Δ_g ≈ 140–300, positivity is expected to be **false** (DG_SCAN; Δ_g = 100 is contradicted at model level).
4. **Recommendation for the headline.**
   * Keep the conditional theorem at Δ_g = 1000 as the headline. Its (Pos) evidence is [NE±] / Theorem W / Theorem L′.
   * Add T0 as the unconditional layer, worded as "leading-Regge positivity for every heavy state at every Δ_g ≥ 1000 [certified]".
   * State explicitly that the beyond-leading-Regge positivity cannot be proved at the quoted gaps by any known method.
   * Pursue T1 or T∞ only if the bounded CT probe (Phase 1, item P1-B) succeeds.

---

## 1. Partition of the heavy spectrum for Ω^EE, and current status

Coordinates:
* μ = m/Λ, with the twist edge at μ_edge(J) ≈ 1 + J/Δ_g;
* b̂ = (2J+1)Λ/m = 2Λ sinh(β/2);
* x = ν/m and ξ = ν/Δ;
* the Regge window is ν ≲ 50 (or ν < 3/β); the bulk is x ∈ [0.1, 0.85–0.9]; the edge is x > 0.85–0.9, up to 1/μ.

The per-state value is `λ·N_s = Φ_flat + curv_LR + δ` (POS_ERRORS §2.1). Rows 0.x are the leading-Regge (P0) layer; rows 1.x–9.x are
the δ layer, which is where (Pos) actually lives.

### 1.1 Sectors and light states

| sector | structure | status / missing |
|---|---|---|
| even J ≥ 4 | ++, +− (2×2), plus OO (1×1); the EO entries cancel in the S⁰ sum (PARITY_ODD_ROUTEA §3) | OO is handled through Γ_OO = Γ_++: [P given helicity lemma (A) [D]] at order 0, [NE] beyond (residual O(m⁻²) c₂). **For a proof, treat OO as its own kind in every certificate (item M17); do not use the identity.** |
| odd J ≥ 5 | +− only (parity-odd structure), 1×1 | the graviton amplitude is not certified at any odd J (M5) |
| J = 0, 2 | ++ and OO | graviton EE/C3H amplitudes certified (frozen ⟨TTO⟩) |
| heavy multi-traces (twist ≥ Δ_g) | same block shapes | included in "heavy" (REGGE_TREE route N) |
| **light, under A1** | T, plus [TT]_{n,ℓ} with twist < Δ_g | excluded by L0, the double zeros at k = 2, 3 at tree level [D] (parity-odd [TT] double zeros [D-sketch]). **No per-state positivity needed.** |
| **light, under A1′** | J = 0, 2 single traces of any Δ_L | **(Pos_light) needed.** Flat [P]. Leading-Regge form [P given the LR form] for Δ_L ≥ 2.6 (J = 0) / 4.5 (J = 2). For Δ_L = O(1) the LR form is **not** the AdS image (DSKETCH_UP); J = 0 is positive iff μ² > −0.545Λ² [D given (E′)], with (E′) [D-sketch]. KK accumulation [I]. **Recommendation:** state the unconditional theorems under A1, and treat A1′ by adding an explicit light-spectrum condition (item M18). |

### 1.2 Leading-Regge layer (P0, "hybrid floor")

| # | region | current best grade, and smallest Δ_g where it holds | missing statement |
|---|---|---|---|
| 0.0 | flat limit, every state | **[certified]** Φ_flat ≥ N/5, Λ-independent, so every Δ_g (CERT_CLOSE §1) | — |
| 0.1 | m ≥ 1.3Λ, b̂ ≤ 119 | [P + DV] Lemma C′, monotone in Λ: **Δ_g ≥ 1000**, but **only at the sampled μ ∈ {1.3, 1.5, 2, 3, 5, 10, 20} and sampled b̂** | interval-in-μ version; μ > 20 (M1) |
| 0.2 | m ≥ Λ, b̂ ∈ [114.06, 2Λ] (110.74 for μ ≥ 1.3; 107.51 for μ ≥ 1.5) | **[D + certified] for Δ_g ≥ 1000** (ROUTEA-U, GAP_LARGESPIN U3) | the derivation is ROUTEA's [D], to be written out as [P] (part of M1) |
| 0.3 | b̂ ≥ 0.2456Λ (F_LR > 0); b̂ ≥ 0.37211Λ (r < 1) | **[P + certified] for Δ_g ≥ 1000**, no upper limit (U1) | b̂ > 2Λ: hybrid vs F_LR is [NE] (M1) |
| 0.4 | m ∈ [Λ, 1.3Λ), b̂ < 114 | Lemma C sign-free near the edge only for J ≤ 13 / 17 / 23 / 35 / 61 at Δ_g = 10³ / 2·10³ / 5·10³ / 2·10⁴ / 10⁵ [certified flat + P]. Otherwise the hybrid is computed positive [NE] (≥ 0.8 N_s) | **the "intermediate-b̂ hole for m < 1.3Λ"** (b̂ from about 2J_LC up to 114, i.e. J ≈ 14–56 at μ ≈ 1 at Δ_g = 1000) (M1) |

### 1.3 Beyond-leading-Regge layer (δ): the content of (Pos)

| # | region | current best grade, and smallest Δ_g where it holds | exact missing statement |
|---|---|---|---|
| 1.1 | near gap, chart A, J ∈ {0, 2, 4, 5, 6, 7}, μ ≤ 1.3 | 1000: [NE±] (L ≥ 5.6 at the twist edge). **Every Δ_g ≥ 1000: [P given B_F + windows]** (GAP_UNIFORM U7). **Δ_g ≥ Δ₁\*(J) = 0.9·10⁵ (J = 0) … 2.5·10⁵ (J = 7): [P given TAIL, Π ≤ 3, K₂ part bounded, ξ-continuity]**, with the graviton amplitude certified only at J = 0, 2 (GAP_MORSE E″) | **B_F** on [1000, Δ₁\*) (infeasible, §5); **TAIL = CT** (M7); **Π ≤ 3** (M2); **ξ-continuity** (M3); **K₂-part bound** (M6); graviton J = 4–7 (M5); windows ν < 50 (M4) and x > 0.9 (M8); the seam images (inside M8) |
| 1.2 | near gap, J = 8 | 1000: [NE±] λ-level only. 2000: sign-free L = 5.2. B_F fails with the certified floor (F = 0.85). Δ₁\* = 4.2·10⁵ (scalar data) | as 1.1, plus the certified hybrid floor (M1) to raise the tolerance (F_hyb = 2.7) |
| 1.3 | near gap, J = 9–32, μ ≤ 1.05–1.10 | 1000: [NE±] for 25 of 31 spins. [NE+M] at the twist edges of even J = 22–32. H_hold [M] on held states; H_cont [NE] between masses. 2000: J ≤ 13 sign-free, J = 14–21 λ-level, J = 27/29 [NE+M], others not run. **No statement at any other gap** (sign-free fails for J ≥ 12–14) | everything in 1.1, plus Laplace certificates for J ≥ 9 (none exist); **H_hold**, **H_cont**; the even-J edge term (needs direct edge rows) |
| 1.4 | near-gap large spin 26 ≤ J ≤ 344, μ ≤ 1.3 (Δ_g = 1000) | **Theorem W [P given P0, P1′, W, AG_U, EDGE-R] + [DV] at 47 sampled states**, sampled J\* = 26 (or 38 if only ≤ 1.2 S_E,t is granted). Between samples [NE]. Other gaps: chart-B images under H_mono-B [M] | **P0** (row 0.4); **P1′** beyond the LR kernel [D-sketch] (M12); **W** (spinning O(1/m⁴) window coefficients, kernel O(r²)) [M] (M13); **AG_U** [A] (M10); **EDGE-R** for the graviton [M] (M8/M14); the boundary factor 2 [D-sketch]; sampled only, float only |
| 1.5 | **crossover band** 33 ≤ J < J_W(Δ_g) ≈ 26.5Λ/Λ₀, above the seam, at Δ_g > 1000 (generically J ~ √Δ_g … J_W) | **no grade at any Δ_g > 1000** (supplement, "A band with no transfer route") | a uniform-in-J control of δ (M9), plus O2 |
| 1.6 | b̂ ≥ 0.37211Λ (J ≥ J\*(m), 223 at the 1000 twist edge), J ≤ 10⁶ | **Theorem L′ [D given P_L′(i), AG_U, EDGE]**. Floor [P + certified] for every Δ_g ≥ 1000. Premises carried to other gaps with AdS-unit constants [M]; tolerances grow ∝ Λ² | **P_L′(i)** magnitudes: [DV] at J ≤ 344, [M] above (M11); **AG_U** (M10); **EDGE** (non-binding, ≥ 600× margin, but unproved) |
| 1.7 | J > 10⁶ | no theorem as stated (the tolerance falls like the floor) | large-β asymptotics of the O(1/m²) symbols (M11) |
| 1.8 | m ≥ 1.3Λ, b̂ < 5 (strip S_D) | 1000: **[NE±] L ≥ 8.7** (flat [certified] − Lemma C [P] − engine envelope [D given RD] with [NE] coefficients; the bulk ν ≥ 20 envelope is a model input, ×76 slack) | a certified 1/Δ engine truncation (M4); the bulk ν ≥ 20 for J > 0 needs CT (M7) |
| 1.9 | m ≥ 1.3Λ, 5 ≤ b̂ < 0.37Λ | H_R′: [NE+M] for b̂ < b̂_A ≈ 8.5–15, [D-sketch + M] above. **No theorem** | O2 + O3 (M9, M10). The sign-free tolerance falls to 10²–10⁻¹ c₂ (POS_ERRORS §2.2), so only signed / complex-ν methods can work |
| 1.10 | fixed J, μ → ∞ | inside 1.8 for b̂ < 5; the engine covers ν ≲ 12–20 | the matching window 50 ≲ ν ≪ m ("CHMRS matching property (ii)" [D-sketch]) (M4) |
| 1.11 | the light-cone edge x → 1/μ at twist-edge states (every Δ_g) | 1000: direct nodes to x ≈ 0.93–0.94, then a sign-free edge term (this sets [NE+M] for even J = 22–32). At large Δ_g, 1 − ξ_max ≈ (J + 3/2)/Δ_g, inside the Airy region | **fold control** (M8a) or **functional redefinition** (M8b) |

---

## 2. Missing items: route, effort, risk, dependencies

"Route" names the most promising concrete construction and the code it would start from. Items marked † are needed only in the
fixed-β-series architecture of Theorems W and L′. The uniform programme T∞ (M9–M10) supersedes them.

| id | missing statement (premise names) | most promising route (estimate, certified computation, existing code) | effort (ah) | risk | depends on |
|---|---|---|---|---|---|
| **M1** | **P0 everywhere:** the intermediate-b̂ hole for m < 1.3Λ; C′ at sampled μ only; μ > 20; hybrid vs F_LR at b̂ > 2Λ | (i) Write `curv = g(b̂, u)/Λ² + Λ⁻⁴h` with ROUTEA-U's alternating-Taylor remainder (`gaplarge/gl7_routea_uniform.py`, exact identity ∫g₁ = −s³r/24). Certify g on (b̂, u) boxes with Arb Bessel integrals (`lrtight` cells, `poserr/curv_lemma.py`), against the **cell-wise certified flat value**, not the global N/5 (`cert/nf_cert.py` machinery). (ii) C′ on μ-intervals: C′ is built from monotone envelopes (LARGESPIN_ANALYTIC §2), so take a max over the cell. (iii) μ ≥ 20: Lemma C at small β, monotone in u. (iv) b̂ > 2Λ: write hybrid − F_LR as an explicit Bessel-tail difference. Λ-uniformity: Lemma-C monotonicity plus a Λ-ball grid as in U1 (`gaplarge/gl2_lrinf_allgaps.py`) | 10–16 | L | — |
| M2 | **Π ≤ 3** (propagation of a relative error through the reference ratio and the a_h fit) | Redefine Γ against the analytic flat envelope D_h(ξ)·g_env, with no reference state and no fitted a_h (the a_h normalisation is a congruence by W1c [P], so it cannot change the sign). Then Π = 1 by construction. Alternative: certify the reference and fit nodes at small ξ (expensive, since ξ = 0.1 sets Δ₁) | 3–6 | L–M | M5 |
| M3 | **ξ-continuity** between certified nodes | δξ as an 8th Taylor-model variable (`gapmorse/gm_cert8`, `gm_ball8`, already validated on 6 intervals). About 75 intervals per spin and kind on [0.1, 0.85], ≈ 2 CPU-h per spin. Δ₁ grows ×4–8 (sets T1's Δ\* ≈ 10⁶) | 2–4 (+30–60 CPU-h) | L | M5 for graviton kinds |
| M4 | **Small-ξ window and matching** (ν < 50; 50 ≲ ν ≲ 0.1Δ; μ → ∞). Includes the **D′-chain [D-sketch] sub-step** (convergence radius of the x₀-dependent slow factors on the raw-radius-1/3 polydisc) | (a) The sub-step: TM-certify analyticity of the slow factors on the polydisc (`gapc3/gc_ball.py` |f/f*−1| < 1 checks) and give a Cauchy bound for the a-series tail. This makes Prop. D′ (RD_CLOSE §2.3) fully [D]/[P], and H_RS for Ω^Q with it. (b) Certified 1/Δ engine truncation at fixed ν (engine orders are exact rationals × one base integral, HRS_CLOSE). (c) Matching: a Morse chart rescaled to the Regge corner (crystal at a ≈ ξ/2 → 0; PCYC §2.1's "Regge-scaled system", where crystal = D is [P]), so that Δ₁(ξ) does not blow up as ξ → 0. Prop. D′'s far-region factor e^{−(0.336 − πξ)Δ} is harmless at Δ ≥ 10⁵ for ξ ≤ 0.1 | (a) 3–6; (b) 4–8; (c) 10–20 | (a) L, (b) M, (c) M–H | (c) on M7 |
| M5 | **Laplace remainder for graviton J ≥ 4 / odd J** (two parity-even structures + OO; parity-odd +− at odd J) | Use P2_SPINNING's exact ℚ(Δ, J) conserved ⟨TTO⟩ solver (no float null vectors), the exact shadow matrix (KKSD 1809.05111, already used in P2_SPINNING §7), and the helicity basis T0(Δ). Write the prefactor in closed exp–log form, as GAP_MORSE §3.2 did for Z·Z, W·W, Z·W, replacing the null-circle sums that cap the radius at 0.014–0.032. Use a Δ-ball 1/Δ ∈ [0, 1/Δ₀] (the TM backend accepts balls). Code: `gapmorse/gm_spin`, `gapc3/gc_node.py`, `p2spin` | 8–16 (+CPU) | M | — |
| M6 | **K₂-part bound** (sup over m of \|Γ_{K₂,m}\| at fixed (J, x)) | K₂ is (c₀ + c₂ε² + c₄ε⁴) with Laplace coefficients that are Gaussian expectations of exact Taylor data. Enclose them as functions of 1/m on [0, 1/m₀] (Δ enters rationally through ⟨TTO⟩ and the x ↔ ξ map) with the same TM evaluator | 3–6 | L–M | M5 |
| **M7** | **CT = TAIL = Lemma R (cycle rigidity) = "contour premise":** D′ ≃ (certified Morse ball) + chains with Re Φ_ξ ≤ Re Φ\* − δ(ξ), with explicit δ, for ξ in [0.1, 0.85] (and complex-thickened for M10) | Structural simplification already in the record: **Φ_ξ and the singular divisors are kind-independent** (RD_CLOSE §1.3: EE = a·Q·scalar, Q with integer exponents on the scalar's divisors). So CT is **one statement about the scalar phase**, and EE inherits it up to boundary integrability, which RD_CLOSE §2.2 already checks on D′. Routes: **(a) ξ-homotopy.** Start at small ξ, where Prop. D′ plus (M4a) give "D′ = crystal + exponentially small". Continue in ξ using the Δ-independence of Stokes data in chart A (GAP_UNIFORM §4.2). The obligation is that no saddle with n_σ ≠ 0 reaches the crystal's height. This needs a certified critical-point census of Φ_ξ (interval Newton / Krawczyk on a compactification with explicit divisor neighbourhoods; a BKK bound on the number of isolated critical points). PCYC's census found 0–276 points **above** the crystal per run, not saturated, so they must be shown to have n_σ = 0. **(b) Explicit deformation.** z → z + i t·v(z) along a truncated descent field, with a cell-wise certified sign of ∂_t Re Φ. 7 real dimensions: the cell count is the risk. **(c) Oracle idea (untested).** The exact scalar (Mellin/Barnes, certifiable as in CERTIFIED_NEARGAP_A) can be compared with the crystal-thimble Laplace at several Δ; a Turán/Nazarov-type lemma for exponential sums with a counted number of terms would then exclude contributing higher saddles for all Δ. Even then an explicit chain is still needed for the EE tail. | **probe 6–10; full 40–100** | **VH** | M4(a) |
| M8 | **Fold / edge window** (x → 1/μ; **EDGE-R** and EDGE at 1000; the seam images) | **(a) Certified fold chart.** A Chester–Friedman–Ursell / Thom-splitting normal form (6 Morse directions + 1 cubic), using the same Bombieri-norm TM machinery as Lemma M, with a uniform Airy remainder. The ψ_k ~ (1−p)⁵ suppression makes the Airy cap cost ~m⁻² relative (GAP_UNIFORM §4.5 [D-sketch]). **(b) Change the definition of the functional.** Use weights with support p ≤ p_max < 1 (e.g. 0.9), so that x ≤ p_max at every state and every gap and the fold is never reached. Cost: a new NF02-type LP and certificate (`robust/nf_*`, `cert/nf_corner.py`; CERT_CLOSE ran in 2.5 CPU-h); the cubic bound changes (L_A and the floor δ′); the 1000 evidence does **not** carry over to the new functional. It would also remove EDGE-R's problem region at 1000 for that functional | (a) 15–30; (b) 4–8 for the certificate + 2 for the new bound | (a) H; (b) M (flat positivity with reduced support is untested) | (a) on M7 |
| M9 | **Uniform-in-J control** (crossover band, rows 1.5 and 1.9; also replaces the chart-A/B split) | Certified Morse charts on boxes in (ξ, j = J/Δ), with P_J split into the two travelling waves W± **inside the integrand** (AG_PROOF Lemma F, exact ₂F₁ split). The j-dependence then sits in the phase (±i j θ), and the saddle moves analytically from the crystal at j = 0. No coalescence is expected away from the edge and the Regge corner (to be certified). The tolerance target uses the certified hybrid floor (M1). The amplitude is uniform in J only while J/√Δ is bounded, so it must be treated in the phase | 20–40 | H | M5, M7, M1 |
| M10 | **AG_U beyond all orders** (the large-b̂ cancellation; Theorems W and L′; row 1.9) | **(a)** With M9 at complex ξ (Im ν ≤ 25/(0.97β), i.e. Im ξ ~ 10⁻⁴ at Δ ~ 10⁶): certified G_± on the rectangles Q_±, then Lemma P4W-χ [P] with certified inputs. AG_U becomes a theorem at Δ ≥ Δ₁, given CT. **(b) J ≳ 100 and large Δ_g.** AG_PROOF §7.1's closure Lemma R [P] reduces AG_U to PW-AG for the explicit, CFT-data-free plane-wave images 𝒜_ν(p) (P2_CONTRACT's image engine). The Bernstein spread ≈ 1.5/√τ shrinks with the gap, so the near-edge obstruction (N2) disappears for β ≫ Δ_g^{−1/4}. Needs the spinning Bernstein representation [D-sketch → D] and a certified PW-AG. Neither helps at Δ_g = 1000, J ≤ 50 | (a) 15–25 after M9; (b) 15–30 | (a) H; (b) H | (a) M7, M9; (b) — |
| M11† | **P_L′(i) magnitudes** for J > 344, and the **J > 10⁶** tail | The O(1/m²) symbols are functions of β alone (GAP_LARGESPIN §4.1 [D]). Port the p2con three-channel contraction to Arb on β-cells, and add a large-β asymptotic lemma | 6–10 | M | — |
| M12† | **P1′** beyond the leading-Regge kernel (the O(1/m) part is a congruence) | Prove W1b's antisymmetry of the odd part of g(m) for the full spinning integrand: the D-operators commute with the flip [D], and the seeds are odd only ∝ (Δ₁₂ − Δ₃₄) [DV]. Then real couplings drop it exactly; the parity-mixing sector enters only through Im(λλ̄) | 4–8 | M | — |
| M13† | **W** (Regge-window bound: spinning O(1/m⁴) coefficients, kernel O(r²); boundary factor 2) | Extend `p2con` to 4th order (symbols + kernel r²). Certify convergence of the window series (successive ratio ≤ 6.5·10⁻³) with an explicit Cauchy radius; write the boundary amplitude lemma | 8–15 | M | — |
| M14 | **EDGE-R at 1000** (graviton E_χ at the 30 W states with J ≤ 51) | Needs the n3 C3H fix (another agent is on N3C3H_FIX), the batch-3 n3 jobs (≈ 90 CPU-h external), a denser 0.0125 grid and lower MC noise. **This gives [NE], never a proof.** A certified version would be M8(a) at Δ = 1000, which is infeasible (§5) | (other agent) + 2–4 fold-in | M | N3C3H_FIX |
| M15 | **B_F** on 1000 ≤ Δ_g < Δ₁\*(J) | **Cannot be proved** (§5). Mitigation only: lower Δ₁ with (i) Laplace orders beyond K₂ computed exactly in balls (≈ ×2) and (ii) a C/Arb Taylor-model port at degree 12–16 (the remainders cap the radius, not the branch checks). Floor about Δ₁ ≈ 2·10⁴ (J = 0) to 10⁵ (J = 8), set by the O(1) variation of the elementary factors at radius 0.1–0.2 | 10–20 | M (for the ×2–4 gain) | M5 |
| M16 | **H_hold, H_cont, between-samples** (Theorem W sampled; J = 11, 13, … held states; masses between budgets) | In the certified programme these become interval coverage (M3, M9). At 1000 only more samples: own rows at the held masses (external) give [NE±], not [P] | 2–4 + external CPU | L (as evidence) | — |
| M17 | **Parity:** OO as an independent kind; EO cancellation in the S⁰ sum | Certify OO directly (it is a third kind on the same phase) instead of Γ_OO = Γ_++. Check that EO cancels structurally (n₀ → −n₀ symmetry) at all orders: an algebraic lemma | 2–4 | L | M5 |
| M18 | **Light states under A1′** ((E′), the J = 0 criterion μ² > −0.545Λ², light J = 2, KK accumulation [I]) | Recommended: state T0/T1 under A1. For A1′, carry "every light J ≤ 2 single trace satisfies the explicit action-positivity condition" as a hypothesis on the spectrum. Making (E′) [D] (CHMRS's light-side formula for the spinning k = 2, 3 functional) is about 6–10 ah | 0 (state under A1) / 6–10 | M | — |
| M19 | **Implication-side [D-sketch] links** (not (Pos), but needed for a bound headline): parity-odd [TT] double zeros; H_K's EE step (spinning partial wave → shifted scalar waves, not checked end to end); K's k = 3 j′ = 1 sign and N̂₃ identification; the k = 2 α_T² part; KKSDZ loophole [Q] | End-to-end numerical checks (INDEP-style) for the EE step and the j′ = 1 sign; a written lemma for the parity-odd double zeros. None changes a number | 6–12 | L–M | — |
| M20 | **Definition of the functional at every Δ_g** (continuum k = 2 de-shift weight with an all-t bound ε_W exists only at 300 and 1000; grids at 2000 and 5000; NU_TAIL_CHECK's grid defect) | Build the continuum (strip-analytic) X(Λ) family with a Λ-uniform ε_W: use ψ₂ + S/Λ² + O(Λ⁻⁴) (GAP_SCALING §1.2) with an explicit remainder, or certified Abel solves on Λ-balls (`deshift/d_abel.py`, `cert/k2_cont.py`). **Required for any "all Δ_g ≥ Δ\*" statement**, because the object must be defined for every Δ_g | 4–8 | L–M | — |

---

## 3. Best achievable target

### 3.1 Targets ranked by achievability

| target | statement | Δ\* | conditional on | effort (ah) | success probability |
|---|---|---|---|---|---|
| **T0** | Hybrid (flat + exact leading-Regge AdS curvature) action of Ω^EE is ≥ c·N_s on every heavy state, every sector, every Δ_g ≥ Δ\* | **1000** (could be pushed toward 300 by re-running LR∞ / ROUTEA-U there, +4–6 ah) | nothing (A1/A2 are not even used: it is a statement about blocks) | 10–16 (+ M20 for the functional definition: 4–8) | ≈ 85 % |
| T0′ | T0 plus δ on the strip S_D (m ≥ 1.3Λ, b̂ < 5) and the Regge window | 1000 | the bulk ν ≥ 20 at J > 0 still needs CT, so T0′ is **not** unconditional beyond the Regge window | +10–15 | ≈ 50 % (and it stays conditional) |
| **T1** | Ω^EE ⪰ 0 on **every** heavy state of spin J ≤ J₀ (J₀ = 8; possibly 12–16 with the certified hybrid floor), all masses, every Δ_g ≥ Δ\* | **≈ 10⁶** (max over J ≤ 8 of Δ₁\* after ξ-intervals: 0.9–4.2·10⁵ × 4–8; M15 could lower it ×2–4) | **CT only** (plus A1, A2), if M8(b) is adopted; otherwise also the fold (M8a) | 35–60 | ≈ 30–40 % for "given CT"; CT itself ≲ 20 % |
| T∞ | Ω^EE ⪰ 0 on every heavy state, every Δ_g ≥ Δ\* | **≈ 10⁶–10⁷** (planning 10⁷: large-b̂ states need relative accuracy ~b̂⁻³ unless the complex-ν route M10(a) works) | CT only (plus A1, A2, and M8(b) or M8(a)) | 120–200 | ≲ 20–25 % given CT; overall ≲ 5–10 % |
| T∞-uncond | the same with CT proved | ≈ 10⁶–10⁷ | nothing | T∞ + M7 full (40–100) | ≲ 5 % on the paper's timescale |

**Answer to "smallest Δ\*":**
* For the statement the authors want (every heavy state, unconditional), **no finite Δ\* is reachable on the paper's timescale.**
  The best long-shot is about 10⁶–10⁷, after CT is proved, and CT is the item most likely to fail.
* For an unconditional statement of the type "every heavy state": **T0 at Δ\* = 1000**, about the leading-Regge action.
* For the full action with a single named mathematical (not physics) premise: **T1 at Δ\* ≈ 10⁶**, for finitely many spins.

### 3.2 Ordered work plan (at most 2 agents at a time)

**Phase 1 (≈ 1 wall-day). Decide feasibility while banking T0.**
* **P1-A (agent A): T0.** M1 (i)–(iv), then M20, then write-up "Theorem T0 [certified]". 14–24 ah. Low risk.
* **P1-B (agent B): CT probe, bounded at 8–10 ah.**
  1. Make M4(a) rigorous: certified slow-factor radius, so Prop. D′ becomes [P].
  2. On the **scalar** phase at J = 0 and one ξ-cell (ξ ∈ [0.3, 0.32]), attempt CT by route (a): certified critical-point enclosure in a compactified box, plus a homotopy from ξ = 0.1. Measure the cell count route (b) would need.

  **Go/no-go criterion:** a complete certified census on one ξ-cell in ≤ 6 CPU-h, **or** a deformation whose cell count extrapolates to ≤ 10⁷ cells for [0.1, 0.85]. If neither, stop: T1 and T∞ are infeasible in this project and the paper ships with T0 + the conditional theorem.

**Phase 2 (only if Phase 1 is go; ≈ 2–3 wall-days). T1 given CT.**
* **Agent A:** M5 (graviton J = 4–8, odd J, OO), then M17, then M6.
* **Agent B:** M2 (Π = 1 redefinition), then M3 (ξ-intervals for J ≤ 8, all kinds; CPU-heavy), then M4(b).
* **Then together:** the fold decision. Prefer M8(b), the redefined functional with p ≤ p_max, and certify its flat floor. Otherwise
  try M8(a).
* Write-up: "T1: given CT, Ω^EE (or Ω^EE_cut) ⪰ 0 on every heavy state with J ≤ 8 for Δ_g ≥ ~10⁶" [P given CT + certified].

**Phase 3 (research programme; ≈ 1–2 wall-weeks; continue only if Phase 2 lands with margin).**
* **Agent A:** M9, the uniform (ξ, j) charts. **Agent B:** M4(c), the Regge-corner matching.
* **Then:** M10(a) (complex-ξ travelling waves → certified Lemma P4W-χ) and M11 for the β → ∞ tail.
* **Optional, in parallel:** M15, to lower Δ\* toward 10⁵.
* **Last:** CT in full (M7), if Phase 1's probe indicated a route.

**Independent of the phases (any free slot).** M19 (implication-side links) and M18 (A1′ wording). Both are cheap and improve the
existing conditional headline whatever happens to (Pos).

**Not recommended.**
* Theorem-W-architecture items M11–M13 for their own sake. They improve the 1000 conditional theorem by turning [M] into [DV], but
  they cannot produce an unconditional statement, because AG_U and EDGE-R remain.
* Any attempt to certify δ at Δ_g = 1000 (§5).

---

## 4. Notes on the choices

* **Why T0 is worth stating.** CHMRS and CLS-D leave exactly this step as "sufficiently large gap". T0 makes the leading-order AdS
  (bulk-point) positivity **effective and uniform in the gap**, for every heavy block, both helicities and every coupling vector, with
  explicit constants. It is the strongest unconditional positivity statement the record supports.
* **Why a definitional change (M8b) is legitimate and useful.** The authors allow "the definition of the functional". Cutting the
  weight support below p = 1 removes the fold (O4) at every gap. Every other obstruction (O1–O3) is untouched, so it is worth doing
  only together with T1 or T∞. The headline numbers at 1000 belong to Ω^EE as now defined; a cut functional would be a second object
  with a slightly weaker bound.
* **Rejected lever:** adding a forward term to make the large-b̂ region sign-free. Ω^Q's forward term does exactly that (POS_ERRORS
  §2.5), but NOFORWARD_FUNCTIONAL removed forward terms for good reasons: a ~10¹³ light-side leverage and unanalysed double-trace and
  contact content.
* **Why the plane-wave / Bernstein positivity does not give (Pos) directly.** The representation is exact and positive, but in that
  gauge the kernel corrections cancel an O(ν²/τ) smearing (LARGESPIN_L N-L1), and 14–35 % of the measure sits below the edge at J ≤ 50
  (AG_PROOF N2). It remains useful only as a reduction tool for AG_U at large J and large Δ_g (M10b).

---

## 5. Infeasible in principle (or with the information that exists)

1. **Certified δ at Δ_g = 1000 (and anywhere below ≈ 2·10⁴), by any Laplace-type method.**
   * At Δ = 1000 the asymptotic series' optimal-truncation accuracy, e^{−Δρ_an²/2} with ρ_an ≈ 0.1–0.15, is 7·10⁻³–1.3·10⁻⁵.
     The needed tolerances are rel\* = 3.5·10⁻² (J = 0) down to 5.6·10⁻⁵ (J = 8).
   * Certifying would need a Morse chart of radius ≈ 0.24, where the integrand's elementary factors already change by O(1)
     (GAP_C3 §5, GAP_MORSE §0). The information is not in any finite-order expansion.
   * A non-asymptotic certified evaluation would need either an exact representation of the spinning action, which is blocked
     (SPINNING_MELLIN §5: the CPPR frame reduction is finite for Π̂, the exact Mellin formula exists for C, and they are different
     objects; a cancellation-free Π̂ would be "new mathematics"), or certified 7-D quadrature against an e^{−πν} cancellation (≈ 145
     digits at x = 1; 10³–10⁴ core-hours per state in float, not certified).

   **Consequence:** the conditional 1000 theorem's (Pos) can never be better than [NE±] / Theorem W / Theorem L′.
2. **Hence premise B_F on 1000 ≤ Δ_g < Δ₁\*.** It concerns a continuum of gaps where no certified evaluation is possible (item 1). Finite
   sampling (1000, 2000, 5000, the exact scalar to 10⁴) cannot replace it. An unconditional "all Δ_g ≥ 1000" for the full action is
   therefore out of reach for the same reason as item 1.
3. **Δ_g below ≈ 140–300: the statement is expected to be false.** DG_SCAN's model bracket is 140 ≲ Δ_g\* ≲ 300 [M]; Δ_g = 100
   positivity is contradicted at model level for near-gap J ≳ 18. No proof effort should target this range.
4. **Not infeasible in principle, but with no method in the record (decide before investing):**
   * CT (O1). Every object is an explicit, absolutely convergent integral over an explicit chain (PCYC_CLOSE), so it is "just"
     analysis. But a certified Picard–Lefschetz statement in 7 complex dimensions with a non-saturated census of higher critical points
     has no precedent we know of.
   * AG_U beyond all orders at Δ_g = 1000 (O2). No route there: plane waves fail at J ≤ 50 near the edge, and steepest descent needs
     the crystal thimble beyond ν ≈ 11 with a soft fibre (AG_PROOF §7.3).
   * EDGE-R for the graviton as a **proof** at 1000: only real-axis numerics exist, and those give [NE].

---

## 6. One-paragraph summary for the authors

Every **certified** layer of the positivity evidence (flat N/5, Lemma C/C′, ROUTEA-U, LR∞) already holds for every Δ_g ≥ 1000, and
closing one hole (m < 1.3Λ, b̂ ≲ 114) and three small edges turns it into an unconditional theorem about the **leading-Regge** action
of every heavy state: T0, ≈ 1–2 agent-days, low risk. Positivity of the **full** action rests, at every gap, on the
beyond-leading-Regge correction, and no certified statement about it exists anywhere. Its certification needs:
* a contour theorem (CT) that no method in the record provides;
* a uniform-in-J analysis across J ~ √Δ_g;
* a beyond-all-orders control at large impact parameter;
* fold control at the light-cone edge, or a redefined functional.

Given CT, a certified theorem for finitely many spins is plausible at Δ_g ≳ 10⁶ (T1). A full theorem is a research programme
(≈ 10⁶–10⁷, ≲ 20 % odds given CT). Nothing at Δ_g = 1000 beyond the present [NE±] / conditional grades is attainable by known
methods. Recommended:
* bank T0 now;
* spend one bounded probe (≈ 8–10 ah) on CT;
* keep the conditional Δ_g = 1000 theorem as the headline, unless the probe succeeds.
