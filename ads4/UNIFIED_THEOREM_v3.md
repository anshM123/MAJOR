# R³ and one R⁴ coupling of holographic CFT₃ from ⟨TTTT⟩ light-ray sum rules: the conditional theorem, consolidated (v3)

**Ansh Mishra and Aryan Senthilkumar**

*23 Nov 2026. Writing pass. The only new arithmetic is the `c₁²` coefficient of the ε = 0.040 repair (§3, from A1PRIME §4.3's formula).
It replaces `UNIFIED_THEOREM_v2.md` (UT2) as the source of truth for the letter. UT2's content is kept where nothing superseded it. The
following deliverables, written after UT2, were read and integrated (in this order): EESTEP_CLOSE (EES), HK_RESIDUAL (HKR), CERT_CLOSE
(CERT), POS_ERRORS (PE), A1PRIME (A1P), PCYC_CLOSE (PCYC), NEARGAP2 (NG2), C1_VALUE (C1V), R4ARC_HEFT (R4H). Where they conflict, the
later document wins unless §11 says otherwise. Every source document is unreviewed except the two marked "independent" (REPRO,
REPRO_ENGINE).*

**One item is still being computed elsewhere:** the high-spin near-gap work (`HIGHSPIN.md`): near-gap J ≥ 12 at Δ_g = 1000, the
J = 29 light-cone edge, and H_LJ(1000). Every number or grade that it may change is marked **⟦HIGHSPIN⟧** and carries the grade
current as of this writing (NG2, PE). Nothing here anticipates its result.

---

> **Update (LARGESPIN_ANALYTIC.md, unreviewed): no large-spin hypothesis is replaced by a proof. H_R′ is upgraded [M] → [D-sketch + M]; H_LJ(1000) and H_bulk keep their grades. A possible convention problem (N2) affects H_LJ's leverage.**
> * **Common missing input:** the beyond-leading-Regge correction inside the Regge window (ν ≲ 3/β). It needs the ⟨TTTT⟩ heavy-block Regge formula at subleading order (the spinning analogue of CHMRS (3.17)–(3.20)), which is not derived.
> * **Lemma I [P]:** for an image ∫P_{−½+iν}(cosh t)σ(dt), the action is ∫Ω_LR(t)σ(dt). This gives explicit tolerance profiles for H_LJ (≥ 0.062 of the ring weight at t = β) and tolerated sizes of candidate correction shapes (ring amplitude a ≥ 80/m, ν² a ≥ 37/m², index-shifted a ≥ 10.8/m, …). The remaining model input R_win (non-ring window correction): exact scalar ≈ 10 c₂ (×38); graviton J = 100 ≈ 550 c₂ (leverage 1.1–1.5); no data at J ≥ 150.
> * **H_R′ [D-sketch + M] for b̂ ≥ b̂_A(m).** Lemma C′ [P + DV], a sign-free curvature bound using the decay of Φ_k, makes the hybrid floor positive for b̂ ≤ 119 at every m/Λ from 1.3 to 20 (the old Lemma C reached b̂ < 10). With ROUTEA, the leading-Regge hybrid is rigorous for all b̂ ≤ 2Λ when m ≥ 1.5Λ. The AG saddle test passes on 5 states. 5 ≤ b̂ < b̂_A stays [NE + M].
> * **H_bulk [M], restated:** the bulk correction may reach ≥ 3.9 % of the flat envelope (9.8 % at 1.3Λ) before any state fails. There is no analytic route.
> * **N2 (open, may weaken H_LJ):** Ω^EE's k = 2 part acts on raw images with the de-shift weight X, not ψ₂, and |X−ψ₂|/ψ₂ = 0.3–0.9 at ν ≤ 10. HLJ's T_B used ψ₂ inside the window. With |X|, the leverage at (1100, 100) drops to ≈ 1.1, and the sign-free ring-amplitude tolerance (3.2) is below the measured excess (2.7–5.9). Which convention the fw/pw Γ data use must be checked. **Until then quote the H_LJ(1000) leverage as 1.1–4.8 (not 2.7–4.8).**
> * **Corrections:** ROUTEA's certified window [91.2, 2Λ] holds only for m ≥ 1.5Λ (UT3 §5.2 and PE §2.4 quote it for m ≥ 1.3Λ). HLJ §5.2's [D-sketch] (subleading terms are P_{l+n}(cosh β) of the same β) does not by itself make window corrections harmless (N3). PE's "no sign-free argument can work" is true for δ but false for the curvature.

> **Update (N2_SUBREGGE.md, unreviewed): N2's premise is right (X, not ψ₂, multiplies the window data), but the claimed leverage collapse does not happen. A dir–cv estimator systematic at (1100, 100) is now the dominant uncertainty.**
> * **What the rows are [DV].** Window rows Γ are raw images / (R₋+R₊), and X multiplies them; Re ΣX𝒜[P_β] = Σψ₂P_β holds to 2·10⁻¹². Γ has poles at the real zeros of R₋+R₊ (ν = 1.171, 6.643), which explains HLJ_CLOSE's "local Γ singular at every J". H_win should be stated on the regular aΓ_c with weight |Y|.
> * **T_B with the exact weight:** (1100, 100): 1409 → 1439 (|X|) / 1152 (|Y|). Over all sample states the change is −28 % to +3 %. LARGESPIN's "1.1" and "3.2" were artefacts (ψ₂ + |X−ψ₂| overcounts; a different floor and structure were switched in; signed ring tolerance 31 vs measured 4.5–6.6). The D_h excess is an a_h normalisation, already in de-shift form (cross-checked against k = 3, which shares a_h). The published Γ − κ contains a stripped-κ artefact of ±120–460 c₂.
> * **Corrected H_LJ(1000) leverage at (1100, 100)** (cv estimator, raw-form refit): **3.0–4.2 central; ≥ 4.1 at 2σ with the full ring class.** (364, 64): 20–36. **With the dir estimator: 0.8–1.4.** This dir–cv systematic is unresolved and dominant.
> * **New explicit hypothesis R1 [I] for EE beyond leading Regge:** ring-type window corrections enter in de-shift form. It holds exactly for the a_h mechanism and for the scalar. In stripped form they would cost 27–228× the ring action at J = 150–344.
> * **Subleading Regge:** CHMRS (3.20) is checked against the exact scalar skeleton (9 states). All its terms sit on the state's own ring, with cost 6.5·10⁻⁴ (J = 100) to ≤ 6·10⁻⁵ (J ≥ 150) of the floor, so margins are ×17 to ×10⁴. The graviton window data fit the same ring classes, with non-ring remainder consistent with zero. The **spinning formula is not derived [I]** (obstacle: the EE = Σ_b identity decomposes the kernel, not the spinning block; it needs ⟨TTO⟩ structure ratios). At each finite 1/m² order the corrections are ring-located or index-shifted at the state's β [D-sketch].
> * **Corrections:** the "leverage 1.1" in the LARGESPIN box above is superseded. Quote "3.0–4.2 (cv), ≥ 4.1 at 2σ; 0.8–1.4 with dir". HWIN_THIMBLE: refit the D_h class in de-shift form, after which the (364, 64) remainder is 7–13 c₂.

> **Update (HIGHSPIN.md, unreviewed): J = 29 is resolved and positive; near-gap J = 12, 16 move to [NE±]; no physical state came out negative.**
> * **J = 29 twist edge (1.0286Λ):** λ/N_s = **+0.972 ± 0.018** [NE + M] (was [M]). Method: a full direct n1 row (34 nodes) plus GK15 quadrature of 16 direct edge nodes, with no interpolation. Budget 0.70, so L = 1.39. The one model input is the ν < 50 window: flipping the state needs |Γ| ≥ 3.2·10³ c₂ there with the wrong sign, while the measured values at ν = 51–77 are ≤ 116 c₂. Masses 1.04–1.10Λ: +0.956 … +0.995, L = 2.8 → 9.5.
> * **Diagnosis of the −3.62 [DV]:** the K₂ method fails, not the physics. The Δ_ref = 2·10⁴ reference's leading partial wave has its first zero at x₀ = 0.8008, right at the bad node (|r₂/r₀| = 0.86 vs ≈ 0.02 elsewhere). w2's Hampel filter skips odd J, so the bad node is kept. The +1.39 is spline aliasing. The published +0.948 removed the node, and the shape chain is within 0.025 of the direct result.
> * **The same reference-zero problem in published chains:** J = 31 (+1.208 → +1.02; this one drops, though it stays above the published minimum +0.988), J = 32 (+0.904 → +0.997), and the J = 40 and J = 48 outliers explained. (Correction: not all shifts are favourable; J = 31 drops.)
> * **J = 12:** +0.8528, L = 3.6 / 3.9 / 2.5 (twist edge, 1.02Λ, 1.05Λ). **J = 16:** +0.9402, L = 1.9 / 1.7 / 1.4 / 1.2. Both now **[NE±]**, from full direct rows (38 nodes) that remove the truncation term. Odd J = 21, 23, 27 and even J = 20, 22, 28, 32 edges: [M] → [NE + M]. Not run: J = 26, 30, 48, 10, 14. **H_LJ(1000): unchanged [D-sketch + M].**
> * **Cautions:** the K₂ λ-level Laplace step cannot close a ×3 budget for J ≳ 17, although the actual K₂ error is only 4–5 % of the step and always favourable. The ν < 50 window is now the binding model input at large J. Lowest computed values: J = 21 at 1.03Λ, +0.868; J = 20 twist edge, +0.83 ± 0.14. Code bugs worked around: `p4_odd_analysis.grid` and `l9_analysis.grid` drop direct nodes whose x is not round to 6 decimals.

> **Update (DIRCV.md, unreviewed): the dir–cv discrepancy is resolved. Both estimators are biased, in opposite directions, by one per-wave systematic. H_LJ(1000) leverage at (1100, 100): 2.3–3.4 central, 1.3–2.0 at 2σ (ring classes A + C); with the full ring class (R1), 7–9 central and ≥ 2.7 at 2σ. A leverage of ≈ 1 is excluded at ≥ 2σ in every fit [NE].**
> * **Errors:** the quoted r errors were 1-dof least-squares residual errors. Replica-jackknife errors are 2–20× larger: cv −160 ± 213 / +327 ± 216, dir −1310 ± 526 / −810 ± 529, paired gap 2.4σ. The gap comes entirely from the node at ν = 13.3, where seed 4242 is a 1.4σ fluctuation (4-seed scalar bias −0.0075 ± 0.0021).
> * **Bias:** the CV mean is exact (≤ 10⁻¹⁰), but the graviton carries λ = 0.5 ± 0.25 of the scalar's per-wave bias (three tests: 0.55, 0.58, 0.43). dir assumes λ = 0 and cv assumes λ = 1, so each is disfavoured at ≈ 2σ. The unbiased P_J-cv twin agrees with the pooled result. The k = 3 "miss" disappears with jackknife errors. No heavy tails were found (Hill index 2.6–3.2).
> * **Pooled fits** (7× statistics at ν = 13.3 and 16.3): A + C, 4 nodes: r = −441 ± 180, L 2.3–2.8 (2σ 1.3). A + C, 6 nodes: +147 ± 77, L 2.7–3.4 (2σ 2.0). A + B + C: −102 ± 122, L 7.2–9.0 (2σ 2.9). J ≥ 150 is not feasible (no basis rows; windows at or below the fold), so it stays [M].
> * **Corrections:** this supersedes the leverage quoted in the N2 and LARGESPIN boxes ("3.0–4.2 (cv), ≥ 4.1 at 2σ; dir 0.8–1.4", "1.1–4.8"). Quote "2.3–3.4 (1.3–2.0 at 2σ); ≥ 2.7 at 2σ given R1". The cv window values in HWIN carry a systematic of ≈ 0.5·|Γ_dir − Γ_CV|. Next cheap step: 4× P_J-twin statistics at ν ≈ 13 (≈ 1.5 CPU-h), giving σ_r ≈ 90 with no λ assumption.

> **Update (DIRCV.md §7, unreviewed): with 4× more P_J-twin statistics the error on r did not shrink. The twin has a heavy tail at ν = 13.3 (tail index ≈ 2.3; one replica sits ≈ 8× from the typical one). Without a λ assumption, the 4-node A+C fit reaches leverage 0.8 at 2σ, so leverage ≤ 1 is no longer excluded in every fit.**
> * **P_J-cv, no λ, worse helicity:**
>   * A+C, 4 nodes: r = +534 ± 387, L = 1.8–2.3 central, 0.8 at 2σ.
>   * A+C, 6 nodes: +289 ± 124, L = 2.7–3.4, 1.7 at 2σ.
>   * A+B+C: +181 ± 231, L = 4.7–6.0, 1.6 at 2σ.
> * **λ-corrected numbers unchanged** (combined λ ≈ 0.6 ± 0.25): A+C 2.3–3.4 central (1.3–2.0 at 2σ); A+B+C 7–9 central (2.9 at 2σ). Every central leverage is ≥ 1.8.
> * **The earlier 56-replica P_J-cv cross-check was not robust** (e.g. `++` 4-node −393 ± 174 → +23 ± 380); the bracketed numbers in the DIRCV box above are superseded.
> * **Next:** add window nodes at ν = 11–15 so the fit leans less on ν = 13.3, or use a variance-bounded twin at that node. More plain samples are inefficient.

> **Update (DIRCV.md §8, unreviewed): 5 new window nodes at (1100, 100) (ν = 11.0–15.0, avoiding the D_h zero at 13.09). The A+C fit no longer hinges on ν = 13.3: dropping that node moves r by +6 to −149, inside the errors. Weighted fits (χ²/dof 0.8–1.3) are now required. λ = 0.7 ± 0.25.**
> * **Leverage, worse helicity, central / 2σ:**
>
>   | fit | λ-corrected | λ-free P_J-cv |
>   |---|---|---|
>   | A+C, 9 nodes | 2.2–2.7 / 1.4 | 2.0–2.6 / 1.0–1.3 |
>   | A+C, 11 nodes | 3.1–4.5 / 1.8 | 2.2–2.9 / 1.7 |
>   | A+B+C, 9 nodes | 2.7–3.8 / 1.4–1.7 | 1.6–2.1 / 0.5–0.6 |
>   | A+B+C, 11 nodes | 4.2–5.6 / 2.3–2.5 | 5.1–7.4 / 2.1–3.2 |
>
> * **Every central leverage is ≥ 1.6, and every A+C fit is ≥ 1.0 at 2σ.** The only 2σ entry below 1 is λ-free A+B+C on 9 nodes, where the ν²D coefficient rests on noisy P_J nodes. The 2σ floor is set by P_J noise at ν < 13 (fold side). Grade [NE].
> * The one-contour twin does not reduce variance. P_J-dir is unusable (point drops at ν < 13; only the control variate cancels them).

> **Update (SMALLNU_WINDOW.md, unreviewed): the ν < 50 window is now measured. The 10⁴ c₂ convention is kept only below the lowest measured node (6–18 % of the window weight). Nothing came out negative. J = 29 and J = 21 → [NE±] (thin); J = 12, 16 leverage roughly doubles; J = 23, 27, 31 stay [NE+M] (now limited by the K₂ bulk truncation, not the window).**
> * **Method:** direct n1 points at ν ≥ 20 and n3 points at ν = 10–17.5 on each state's own row, with jackknife errors. n1 fails below ν ≈ 20 (crystal-saddle fold); n1 and n3 agree at ν = 25 and 30. Reference zeros cannot occur in the window.
> * **J = 29 twist edge:** measured Γ ≤ 310 c₂ on ν = 25–47; ν = 10–20 is consistent with 0 (σ up to 730). Signed window: L 0.61 → **1.23**, λ ≥ +0.18. The remaining model input is |Γ| ≤ 10⁴ on ν < 10; the state tolerates 4.5·10⁴ there, 16× the largest known value. Masses 1.04–1.10Λ: L 2.2–5.7. **A sign-free ×3 envelope does not close (L = 0.94)**, because of Monte Carlo noise at ν = 10–20.
> * **J = 21:** a new 34-node direct bulk row puts the twist edge at +0.998 ± 0.020, L 0.58 → 1.20 (1.7–6.2 at 1.03–1.10Λ). **J = 16:** L → 4.4 / 3.7 / 2.6 / 1.7. **J = 12:** L → 7.4 / 7.5 / 3.6.
> * **J = 23 / 27 / 31:** window cost now 0.15 / 0.29 / 0.46, but the K₂ λ-level truncation blocks them (with calibrated truncation L ≈ 1.5 / 1.4). **J = 32:** the window costs 0.26 measured, vs 2.31 under the convention. J = 20 and J = 32 have no budget.
> * **Cautions:** the window weight grows steeply with J (1.4 → 23·10⁻⁵ N_s/c₂ from J = 12 to 32), so for J ≳ 32 it is a live input again. Lemma C/C′ cannot bound Γ (it is beyond leading Regge). HIGHSPIN's odd-J window weight omits a gauge factor (8.52 vs 8.73·10⁻⁵). EE_DG1000_LARGEJ's "direct RQMC cannot reach the precision" is false on ν = 25–50.

> **Update (LARGESPIN_W.md, unreviewed): Route W gives a conditional theorem with an explicit J\*, not an unconditional one. H_LJ(1000) stays [D-sketch + M]; do NOT read it as "proved for J ≥ 50".**
> * **Theorem W** (implication [P] given P0–P4; costs [DV] on sampled states). Premises:
>   * P0: the leading-Regge hybrid floor.
>   * P1: parity, M = 𝔫U†gU with g even in 1/m at fixed β; [D-sketch].
>   * P2: beyond-hybrid O(1/m²) corrections are ring class with |κ_s| ≤ Λ_amp·|κ_s(scalar CHMRS)|, and index-shifted harmonics have |κ| ≤ K_sh. Magnitudes are [NE] at J ≤ 100 and [M] beyond.
>   * P3: R1, i.e. de-shift form; [I] at O(1/m²) for the ν² classes.
>   * P4: O(1/m⁴) ≤ K₄, plus the AG beyond-all-orders hypothesis; [M].
>
>   If the summed Lemma-I cost d ≤ H = f₊₊f₊₋/(f₊₊+f₊₋), then M_Ω ⪰ 0 for every coupling ratio.
> * **J\*** over the samples (even J 26–344, odd J 33–301, m ≤ 1.3Λ):
>   * with (Λ_amp, K_sh, K₄) = (3, 10, 10) and R1: J\* = 50 at margin ×1, 100 at ×2;
>   * with (1, 10, 1): 30 / 38;
>   * if R1 fails (stripped form), there is no J\*: it fails at J = 300 and 344.
> * **New pieces:**
>   * W1a [D+DV]: the 1/m parity comes from a Weyl reflection T(−x,−y; Δ,J) = T(x,y; 3−Δ, −1−J), holding to 10⁻⁴⁰.
>   * W1c [P], congruence lemma: ν-independent normalisation and coupling-basis effects cannot change positivity. This disposes of the measured +21/m term and N2's a_h excess.
>   * W2 [D] / [D-sketch]: every finite order is ring class, so there is no window-smooth class.
>   * R1 at O(1/m) is no longer needed given P1.
>   * The leading term c(ν) = 32γ₃²γ₄/(γ₁γ₂γ₇) is exact.
>   * The (1100, 100) B fits correspond to Λ_amp ≈ 2–3 and sit ×10–100 inside tolerance.
> * **Decisive residuals:**
>   * the P2 magnitudes, which need the next-order EE-contracted expansion of the integrand at fixed β;
>   * R1 at O(1/m²) for the ν² classes at J ≳ 200. A sufficient, cheap check: the s = ±1 harmonics of the leading-r EE kernel share one polarisation tensor (`eestep/ee_quad.py`).

> **Update (LARGESPIN_W.md §9, unreviewed): neither premise of Theorem W closed.**
> * **R1's sufficient condition is false [DV].** The leading-r EE kernel has rank exactly 3 in polarisation (singular values 1, 0.24–0.37, 3·10⁻⁴–4·10⁻³; the next is at the 10⁻¹¹–10⁻¹⁶ floor). At leading order the heavy symbols see only the de-shift combination. At O(1/m²) the two other channels can enter in non-de-shift form; P1 forbids that at O(1/m). R1 stays [I], now pinned to these two kernel channels. If they are index-shifted harmonics with rational coefficients, they are tolerated even in stripped form (K_sh ≥ 52 at J = 344).
> * **P2 magnitudes not computed.** Three steps remain:
>   1. conserved ⟨TTO_{Δ,J}⟩ structures in d = 3 for general J;
>   2. second-order heavy symbols at fixed β;
>   3. contraction with the three kernel channels.
>
>   Λ_amp ≈ 2–3 stays [NE] (from J = 100); K_sh and K₄ stay [M].
> * **J\*:** with P3, J\* = 50 (100 at ×2). Without P3, positivity holds on 50 ≤ J ≤ 250 (and J = 301) but fails at J = 300 and 344. So R1 at O(1/m²) is load-bearing only for J ≳ 300.

> **Update (LARGESPIN_L.md, unreviewed): Route L gives no unconditional theorem either. It proves the leading-Regge floor at large impact parameter, gives a conditional Theorem L with J\* ≈ 0.52m (857 at the twist edge), and reduces the residual to two CFT-data-free statements about the EE functional's own plane-wave images.**
> * **Lemma LR∞ [P + DV]:** F_LR ≥ (1 − r)·A₁M₁/Λ² + 2047/(2b̂³) > 0 for all β ≥ 1, all m ≥ Λ and both helicities, with r ≤ 0.58 (β ≥ 1) down to 1.6·10⁻³ (β ≥ 3). The constants are float on a β grid to 64; the tail is [D]. With Lemma C′ and ROUTEA, the leading-Regge floor is rigorous on every b̂ for m ≥ 1.5Λ. This closes the old "b̂ > 2Λ: hybrid [NE]" entry and supplies Route W's P0 for J ≥ 857.
> * **Spin-Gamma representation:** [P] for the scalar heavy factor; [D-sketch] for ⟨TTO⟩⟨OTT⟩. Every heavy block is a positive average of plane waves, with J entering only through a Gamma(J+1) average.
>   * Corollary [D]: R1 holds automatically for the heavy-factor part of any correction. **R1 can fail only through the EE kernel's own plane-wave images, a property of the functional, not of the CFT.**
>   * Lemma W [P given Paley–Wiener on H₂]: a correction of exponential type T has its density within T of the ring.
> * **Theorem L [D given P_L]:** for J ≥ J\*(m) = ⌈(1.0421m − 1)/2⌉, the per-state matrix is ⪰ 0 for every coupling vector.
>   * Premise P_L, part (R1): de-shift form.
>   * Premise P_L, part (Loc): the de-shifted correction's density is ≤ tol(1)/4 at t ≥ 1 and ≤ 10⁻⁷–4·10⁻⁹ at t < 1.
>
>   P_L reduces to two open statements: **(PW-R1)** the EE functional's plane-wave images are in de-shift form, and **(PW-Loc)** they carry no density at t < 1 [I].
> * **Negative findings:**
>   * Corrections do not decay like O(1/J) at fixed twist; they tend to a J-independent lightcone amplitude (scalar ≤ 1.2·10⁻⁴ on ν ≤ 10).
>   * The de-shift tolerance has its minimum of 6.7·10⁻⁴ at J ≈ 350–500, so large J is not easier.
>   * In raw form the tolerance falls like 0.7/J.
>   * ANEC positivity holds only for the sum over states.
>   * Graviton non-ring remainder 2–5·10⁻⁴ (relative) at J = 29–100, 40× the scalar; the sign-free L ≈ 1.3–3 in the hardest band if it stays J-independent.
> * **Inventory:**
>   * J ≥ J\*(m): [D given P_L], with a rigorous floor.
>   * 150 ≤ J < 857 near the gap: the hardest band, still numeric/model.
>   * Route W covers J ≥ 50 given scalar-size coefficients and R1. Both routes stall on the same object: the EE-kernel polarisation contraction beyond leading order.

> **Update (P2_SPINNING.md, unreviewed): steps 1–2 done; step 3 done only inside a labelled model [M]. P2 is not closed, and H_LJ(1000) stays [D-sketch + M]. Do not read this as "Λ_amp = 1 derived".**
> * **Step 1 [DV + Q]: conserved ⟨TTO_{Δ,J}⟩ in d = 3, solved exactly over ℚ(Δ,J).** Generic even J: 2 parity-even + 1 parity-odd structure; odd J: 1 parity-odd. This matches Dymarsky–Kos–Kravchuk–Poland–Simmons-Duffin 1708.05718; conservation checked to 4·10⁻²⁷. The odd-J "+−" structure is parity-odd.
> * **Step 2 [D + DV]: heavy symbols to O(1/m²).** The leading symbols are pure helicity. The shadow side is exactly m → −m. O(1/m) is a basis drift plus a Wigner phase. **P1 is now [D given (S)]**, where (S) is that the shadow coefficients are diagonal in the helicity basis with positive leading entries [D-sketch]; check against 1809.05111.
> * **Step 3:** the CHMRS ∂_p mechanism reproduces CHMRS's O(1/m²) scalar series [DV] (ν⁰P off by exactly 8η, from normalisation). By-product [NE]: N2's unexplained "−8…−9·P" scalar residual is −8ηP at all 8 states. **Model [M]** (helicity-diagonal kernel): the spinning symbol adds nothing to the ν² classes, giving Λ_amp(B) = Λ_amp(C2) = 1, Λ_A = 3.17, Λ_C3 = 2.03. In the conjugate pairing, Λ_C1 ≈ 205. The (1100,100) and (364,64) fits bracket the model (inconclusive).
> * **Theorem W rerun (model constants):** J\* = 50 with R1 or in stripped form at K_sh = 10, K₄ = 10; 30 at K₄ = 1; fails at J = 300 and 344 if K_sh ≥ 30 without R1. **Now binding:** K₄ (below J = 50) and K_sh (J ≥ 300), both [M]. In the model, R1 is not needed at margin ×1. Helicity leakage into the non-de-shift channels is ≲ 10² at O(1/m²), possibly 1.8·10³ in the even basis; estimate κ ≲ 30 [M] against a tolerance of 52–59 at J = 300–344.
> * **Remaining:**
>   * the true three-channel kernel contraction at O(1/m²) (overlaps PW_IMAGES);
>   * K_sh from the kernel's O(r) terms;
>   * K₄ to O(1/m⁴) (cheap);
>   * the (S) check;
>   * the J = 0 engine O(a²) comparison.

> **Update (PW_IMAGES.md, unreviewed): R1 is replaced by an exact identity plus one bounded monopole term. PW-R1 and PW-Loc are false as literally stated, but neither is needed any more. Theorem W keeps J\* = 50 (100 at ×2) WITHOUT P3. Theorem L keeps J\* = 857 with P_L reduced to magnitudes.**
> * **Exact de-shift identity for the EE kernel.** Every b-shifted scalar block in EE = Σ_b satisfies CHMRS's magic-block identity with b-independent coefficients:
>   2η(1+r²)Ĝ¹_b(ν) = (1+i/2ν)Ĝ⁰_b(ν+i) + (1−i/2ν)Ĝ⁰_b(ν−i).
>   It is [P] to r⁸, [DV] to r⁴⁰ at 50 digits, [D] at all orders. So every plane-wave image, in any gauge and with any insertion, is exactly 𝒜_g[F] (CHMRS's scalar de-shift), with F even and analytic for |Im μ| ≤ 3/2. End to end against ee_quad: 2·10⁻⁷–1.1·10⁻⁶. **Route W's stripped (R2) worst case cannot occur.**
> * **PW-R1 is refuted** with respect to the EE operator 𝒜 = 𝒜_g∘c. It needs F(±i/2) = 0, which O(1/m²) radial insertions violate [P]: |y|² gives a ratio of 0.914, r gives 0.500. Whether the full ⟨TTTT⟩ sum cancels this is [I].
> * **The failure is exactly a monopole term, and it is bounded [DV].** The exact cost is ∫ψ₂F/c + D·F(i/2), with D = −9.98·10⁻⁴ (1000) / −1.10·10⁻² (300), ≈ −990/Λ². D is constant in impact parameter, and the equivalent cost kernel stays below 10⁻³ (no e^{t/2} growth).
> * **PW-Loc is false at O(1/m²):** index-shifted harmonics carry density 0.271 per unit coefficient at t < 1 (closed form [P]). It is no longer needed, since every finite-order structure now has an exact cost.
> * **Theorem W without P3 [P given P0, P1, P2, P2′, P4] + [DV]:** the scalar-shaped worst case is tolerated at ×24–29 for J ≥ 150 (stripped form gave ×1.8–7). J\* = 50 at (3, 10, 10); J\* = 100 at ×2. P2′ (the radial-shift family exhausts the O(1/m²) multipliers) is [D-sketch].
> * **Theorem L [D given P_L′]:** J\* = 857 at the gap. The monopole tolerance gives a tolerated Λ_amp of 1240 (J = 857), 555 (10⁴), 166 (10⁵) and 43 (10⁶).
> * **Still needed by both:** the spinning O(1/m²) magnitudes (Λ_amp ≈ 2–3 [NE] for J ≤ 100, [M] above) and the beyond-all-orders premise. Untested: polarisation-tilt O(1/m²) terms.

> **Update (P2_SPINNING.md §7, unreviewed): K₄ ≤ 3.8 in the model; the (S) check FAILS, so the P1 upgrade is withdrawn and P1 stays [D-sketch]; a new gap in P4 (the flat ν⁶ tower) is found. Model J\* = 38.**
> * **K₄:** the scalar seed matches CHMRS's O(1/m⁴) series exactly (ν⁴, ν⁶ coefficients at two β). In the helicity-diagonal model the spinning symbol adds only −(5/24)(η²−1)ν⁴P, so K₄ ≤ 3.6 + 0.21 ≤ 3.8 [M]. There is no model-independent bound.
> * **New gap in P4:** CHMRS's O(1/m⁴) term contains −(η²−1)ν⁶P/72, which is the same for every Δφ (part of the flat-space tower). Costed as a Regge term it is equivalent to K₄ ≈ 35 (J = 26) to ≈ 3000 (J = 120–200), and every sampled J fails. **P4 needs this tower's hybrid subtraction made explicit:** the identification is [DV], the subtraction is [I].
> * **(S) fails:** Karateev–Kravchuk–Simmons-Duffin 1809.05111 shadow matrices mix structures, so the two parity-even structures can mix. P1 stays [D-sketch]. What survives: the leading shadow matrix is diagonal and positive [D], and in the model the O(1/m) off-diagonal is constant, so it cannot affect positivity.
> * **−8η factor:** the first-order term of the moment's normalisation, exact at O(1/m²) and O(1/m⁴). It is harmless (a rescaling); its sign relative to N2's residual is [I].
> * **J\* (model constants, K₄ = 4):** K_sh = 10 gives 38 with R1 or stripped (θ = 1), 50 at ×2 with R1. K_sh = 30 gives 38 with R1 and none stripped. **This holds only if the flat ν⁶ tower is assigned to the hybrid floor, now a named premise.** K_sh is uncomputed [M] and binding at J ≥ 300.
> * The J = 0 engine O(a²) comparison was not done (it needs the full three-channel contraction).

> **Update (REGGE_TREE.md, unreviewed): referee-3 M2 is a presentation defect (the order of limits was never stated), not a derivation gap. Done as CHMRS do it, the Result never needs H_K at tree level; all grades stay.**
> * **H_K's J₀** is the exponent of the uniform bound on dDisc of the **exact** correlator (KKSDZ (4.70)). **Route N [D]:**
>   1. Apply the sum rule at finite c_T. H_K is [D] under J₀ ≤ 1, with margin ½ at k = 2, 3, 5 and 3/2 at k = 4.
>   2. (Pos) makes the heavy side ≥ 0 at every c_T, multi-trace states included.
>   3. Expand only the light side in 1/c_T.
>
>   The tree-level heavy sum is never expanded; CHMRS §4.4 show it need not be analytic in 1/c_T. The order-of-limits step is [D-sketch] given A2. A2 must be restated as tree-level EFT light data up to o(1/c_T), uniformly in spin.
> * **Term-by-term reading (not made):** with a uniform planar J₀ = 2, the corner count fails at k = 2, 3 **and 5**, and for generic k = 4 weights (Q's k = 4 weights pass with margin ½). "J₀ = 2 (Einstein)" is the Wightman/elastic intercept; the heavy tree dDisc is empty for pure Einstein, and for strings its fixed-η intercept is 2 − 9/(2Λ²) with an η⁻² k = 2 moment (CHMRS (4.45)).
> * **Literature [Q]:**
>   * CHMRS (§1.1.1, (1.12)–(1.14), §4.1, §4.4) and CLS-D (§3.4, §3.5.2) use J₀ ≤ 1 non-perturbatively and assume no tree-level J₀ < 2.
>   * KKSDZ §4.5 gives J₀^planar ≤ 2 (4.80); their planar sum rules need J₀ < 3.
> * **Model check:** with the uniform J₀ = 1 bound, the corner integrals converge for every k; with the uniform J₀ = 2 bound they diverge at k = 2, 3, 5; with a stringy Regge-pole model they converge for every k.
> * **Wording:** proposed replacements for letter.tex:109 (A2), :133–135 and :183–185, plus supplement paragraphs "Order of limits" and "Term-wise use (not made)", are in REGGE_TREE §5. Say "every state with Δ−J ≥ Δ_g", not "heavy single traces".
> * **Open:** the size of the light-side remainder O(Δ_g^p/c_T²) is not computed, so keep "at leading order in 1/c_T".

> **Update (INDEP_CHECK2.md, independent reproduction with new code): all five targets reproduced, nothing load-bearing refuted. Two substantive qualifications: Theorem W's J\* = 50 → 64 once the k = 3 piece is included, and Lemma LR∞'s r(β) table uses a cruder constant than the stated formula.**
> * **Magic-block identity:** exact with l and b symbolic to O(r⁶), so it is an identity in b. It is exact to O(r⁸) for every b = −4…4, and holds to 4·10⁻⁴⁸ at 50 digits to r⁴⁰ (including b = 0.37, 7, 10 and complex ν); controls fail at O(1). **The monopole constant D is reproduced to 8 digits, and its t-independence is now derived:** 1/c = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9), so D = (15/2)·defect[1/(4ν²+1)].
> * **Lemma LR∞:** reproduced; the true r(β) is ≤ 8·10⁻⁷ on β ∈ [1, 32]. **Discrepancy:** the stated formula (with ∫|ψ|) gives r(1) = 0.016, while the table (0.577) used Σ|c_n|/(n+1). Both are rigorous; the text must name the constant used. With the tighter constant the bound holds down to β = 0.5 (r = 0.28), which would give J\* = 324 at the gap instead of 857 (needs a certified ∫|ψ₃| and a denser grid).
> * **Order of limits (route N):** sound, and the quotes are verbatim or faithful. Minor fixes:
>   * the CHMRS ε quote is in §1, not §1.1.1;
>   * the "technical assumption" is eq. (5.3), not §5.3;
>   * CLS-D restates the logic in §2.2.2 but never discusses the order of limits, so do not say "the logic of CLS-D";
>   * f⁻ = −1.
> * **Theorem W:** the 2×2 criterion d ≤ H is correct and sharp (200 000 random tests). All J\* table entries are reproduced. The Gershgorin factor is ≤ 1.5 only if f₊₊ ≈ f₊₋ (otherwise up to 3; 1.50–1.53 on the samples). **Main qualification: d omits the k = 3 piece of the corrections. Including it moves J\* 50 → 64 at (3, 10, 10), 30 → 33, and 64 → 80**; the J = 50, 51 states fail at d/H ≈ 1.24. This must be stated as a named premise or used as the J\*. The cache stores the 4a/4b shapes with 1/m² instead of 1/m⁴ (documentation slip only).
> * **⟨TTO⟩ counts:** verified against DKKPS's TTTT.tex table and by an independent helicity count (ℓ ≤ 11).

> **Update (LRINF_TIGHT.md, unreviewed): Lemma LR∞ is certified down to β₀ = 0.37, so Theorem L's J\* at the twist edge goes 857 → 223. The 345–857 band is now inside Theorem L′ at every J (not only sampled J), and Theorems W (J\* = 64 with k = 3, sampled to 344) and L overlap at the gap on 223 ≤ J ≤ 344.**
> * **Constants [certified]:** ∫|ψ₂| = 10.386029907129 and ∫|ψ₃| = 12.835622283960 (Arb, ±7·10⁻²⁰; Sturm-counted sign changes). The old table used Σ|c_n|/(n+1) = 7.9·10⁵ / 1.9·10⁶, which explains the 0.577. A₁,min and A₂′ are exact rationals.
> * **r(β) [certified]:** monotone pieces evaluated in Arb on cells covering [0.1, 64], plus a new explicit tail for β ≥ 64 with r ≤ 2.03·10⁻⁷ (was [D]). sup r = 0.0162 for β ≥ 1 and 0.280 for β ≥ 0.5. **β₀ = 0.37 for r < 1.** F_LR > 0 is certified for β ≥ 0.245. Step 3 (the Hadamard identity) is now [P]; Lemma LR∞ is [P + certified].
> * **New J\*(m) = ⌈(0.37211·m − 1)/2⌉ [certified]:**
>
>   | m | old J\* | new J\* |
>   |---|---|---|
>   | twist edge | 857 | 223 |
>   | 1.1 m_gap | 987 | 250 |
>   | 1.3 m_gap | 1275 | 307 |
>   | 1.5Λ | 780 | 279 |
>   | 2Λ | 1041 | 372 |
>   | 5Λ | 2602 | 929 |
>
>   At J\* = 223 the certified floor is 0.57 of the numerical floor.
> * **Tolerances at J\* = 223 (certified floor) [DV]:**
>   * monopole Q_tol = 4.8·10⁻³;
>   * tolerated Λ_amp: 181 for the full scalar-shaped vector, 1025 for the monopole alone; ≥ 181 up to J = 10⁴, 133 to 10⁵, 32 to 10⁶;
>   * K_tol for sh±1 ≥ 449.
> * **Floor coverage:** the floor is rigorous at every b̂ for m ≥ 1.5Λ. For m ∈ [1.3, 1.5)Λ the [NE] hole shrinks to b̂ ∈ (119, 245). For m < 1.3Λ the floor is [NE] below b̂ = 245. Theorem W's P0 is now certified at its sampled gap states with J ≥ 139.
> * **Still uncovered:** non-sampled 64 < J < 223, J < 64, and each theorem's premises (P_L′ magnitudes [M], beyond-all-orders). This is not an unconditional theorem.

> **Update (P2_CONTRACT.md, unreviewed): the helicity model is replaced by an exact three-channel computation at all 47 Theorem W sample states. All of P2 (Λ_amp, K_sh, P2′, P3) is now computed [DV]. P4 (beyond-hybrid O(1/m⁴) and beyond-all-orders) is the only magnitude premise left, and it is still open [I]: the state-wise scalar O(1/m⁴) remainder costed to ν ≈ m exceeds tolerance at every sampled J ≥ 80 except 301.**
> * **Engine:** an exact leading-Regge EE plane-wave image engine with any polarisation and insertion (all three channels). It reproduces the heavy images (ratio 1.000000), the R∓ normalisation, and CHMRS ∂_p to 10⁻¹⁵.
> * **Computed constants:**
>   * Λ_B (block) = Λ_C2 = 1 exactly.
>   * The kernel's own O(r) term is ring class with κ_B = ±2, giving **Λ_B(+−, odd J) ≈ 2.0 (worse than the model)**; Λ_B(++) ≤ 0.09.
>   * C1(+−) ≈ 200× scalar (harmless); Λ_A ≈ 6–8 (a rescaling).
>   * Off-diagonal ++/+− at O(1/m²) = 0. The O(1/m) image after drift removal = 0, so P1′ is [DV] at leading-Regge kernel level and the (S) failure cannot act.
>   * **The monopole cancels between block and kernel (∓12.8 vs ±12.8), so PW-R1 holds for ⟨TTTT⟩ at O(1/m²).**
>   * Exact O(1/m²) cost/H: ≤ 0.41 (J = 26), ≤ 0.05 (J ≥ 64), ≤ 0.005 (J ≥ 200).
>   * κ_B(+−) ≈ 3.8 matches the (1100, 100) B fits [NE].
> * **Flat ν⁶ tower:** subtracted exactly at O(y): the block reproduces the flat partial waves E₂d⁰⁰ and E₂d⁴⁴, and with the kernel's O(r) term gives exactly the hybrid's stripped flat part. At O(y²) the ++ case matches up to an unexplained +(3/16)J₀y²; +− at this order is [D-sketch].
> * **P4 remainder:** scalar O(1/m⁴) minus the flat tower gives K₄ ≲ 3 for J ≤ 50 (exceptions J = 46: 26.6; J = 26, 34 marginal), 10–16 at J = 64, and 20–7·10³ at J ≥ 80. The cause is the AdS/flat kinematic mismatch when costed to ν ≈ Λ ≈ m, where the fixed-β series is invalid. **P4 needs a ν-window statement instead of the truncated full-range series.**
> * **Theorem W [P given P0, P1′, P4] + [DV]:** J\* = 26 (K₄ = 0), 30 (K₄ = 1), 38 (K₄ = 4), 50 (K₄ = 10), 100 (K₄ = 35), at θ = 1. With the state-wise computed scalar remainder there is no J\*: only J = 30, 33, 38, 42, 50, 64 and 301 pass.
> * **Open:**
>   * P4 (a ν-window bound, plus the spinning O(1/m⁴) symbols and the kernel O(r²) term);
>   * the (3/16)J₀y² residual;
>   * P1′ beyond the leading-Regge kernel;
>   * the J = 0 engine comparison;
>   * behaviour between samples.

> **Update (P4_WINDOW.md, unreviewed): P4 is reformulated as a ν-window statement and reduced to three premises (W, AG′, EDGE). It is NOT proved: the beyond-all-orders part remains an assumption [A], now with explicit constants. Do not write "proved for J ≥ 26 / 42".**
> * **The full-range P4 failure is a series artefact [DV].** On exact scalar heavy moments, the true remainder costs ≤ 1.5·10⁻³ of the budget (K₄-units 0.99 at J = 64, 8.1 at J = 100). The fixed-β series is asymptotic and worsens past x ≈ 0.3–0.4. In the Regge window ν < 3/β it converges (successive ratio ≤ 6.5·10⁻³ at all 47 states) and matches the exact cost to 0.3 %.
> * **Lemma P4W [P] (Cauchy).** cost(G) ≤ window sup-bound + boundary term + an edge term at ν = Λ (suppressed by κ⁻⁶) + the exact monopole. Premises:
>   * (W) a real-axis window bound;
>   * (AG′) exact − flat is a pair of analytic travelling waves on [3/β, Λ], with local decay 0.97β/√(1−x²);
>   * (EDGE) an amplitude bound at ν = Λ.
>
>   The 1/m series is never used outside the window. Checked on the exact scalar: bound 4.0–7.9 % of budget.
> * **Theorem W, 47 states [DV]:** computed part ≤ 0.49 H everywhere (k = 3 included); every state where the full-range series failed now passes. With the scalar-calibrated edge, J\* = 26 (edge multiple tolerated up to 114; spinning O(1/m⁴) coefficient tolerated up to 20.8).
>   * **If the +− edge carries the Ω^U helicity factor (1−x₁²)⁻²**, then (26, 1.026Λ), (30, 1.030Λ), (34, 1.033Λ) and (38, 1.037Λ) fail. J\* becomes 42 (edge multiple 1), 51 (×3), 100 (×10).
>   * The remaining inequality there is |ρ₊₋(x₁)| ≤ 0.15–0.63; the exact scalar has 0.27–0.69.
> * **Theorem L′ range (16 states, J = 223 to 10⁶):** P4 is discharged with ε_R ≤ 0.086 of the certified half-floor (0.28 at Λ_amp = 10). EDGE is not binding (tolerance ≥ 600× the scalar envelope).
> * **(3/16)J₀y² residual resolved:** it is a genuine, harmless beyond-hybrid O(1/m⁴) term.
> * **Remaining:**
>   * (AG′) beyond all orders [A]: local rate [NE] at x ≤ 0.86, [M] at x₁ > 0.94;
>   * the graviton edge amplitude and its +− factor [M];
>   * the spinning O(1/m⁴) window coefficients [M] (margin ×5–10⁴);
>   * the boundary factor 2 [D-sketch];
>   * sampled states only; float only.
>
>   Sensitivity: with threshold 10/κ the edge multiple tolerated drops from 114 to 18.7.

> **Update (EXT_FOLD.md, unreviewed): the external-machine points are folded in. J = 7 and J = 27 → [NE±]; J = 23 and J = 31 are [NE±] at their higher masses but stay [NE+M] at the twist edge; J = 29 still does not close sign-free (L 0.94). No state came out negative.**
> * **J = 7:** twist edge +0.2683, L = 6.7; 1.01–1.05Λ: L = 7.2–25. **[NE±]** (its own direct row; shadow coefficient computed and checked bit-for-bit).
> * **J = 27:** twist edge +0.976 ± 0.019, L 0.28 → **1.64** (signed 2.02); 1.03–1.10Λ: L = 1.93–3.21. **[NE±].**
> * **J = 23:** twist edge +0.923, L 0.45 → 0.85, stays [NE+M]. It is now limited by the degree-7 shape-fit representation term, which Monte Carlo noise amplifies; a noise-aware term would give L ≈ 1.52, but that is a protocol change. 1.03–1.10Λ: L = 1.16–2.33, **[NE±]**.
> * **J = 31:** it has no GK15 edge run, so the tail beyond x = 0.9 must be covered sign-free; twist edge to 1.05Λ: L = 0.31–0.83, stays [NE+M]. 1.075/1.10Λ: L = 1.53/1.91, **[NE±]**. A GK15 edge window on [0.8, 0.970] is needed.
> * **J = 29 sign-free:** not closed. The 16 B seeds give only 23 independent streams: in n3, replica r of seed s uses stream s + r, so consecutive seeds share 7 of 8 replicas; independent seeds must be spaced by ≥ 8. The streams are heavy-tailed (|z| ≈ 11), and the binding inputs are the n1 node at ν = 20 and the 10⁴ convention below ν = 10. Even with σ → 0 on ν = 10–17.5, L = 0.97. **SMALLNU's "≈16× n3 points would close it" is false.**
> * **Checks:**
>   * External points reproduce local recomputation to 1.5·10⁻¹² / 4.5·10⁻⁷.
>   * Direct-vs-K₂ differences are always favourable (median 58–89 c₂).
>   * The K₂ error is 0.6–8 % of the last Laplace step over seven states (HIGHSPIN's "4–5 %" is corrected).
>   * Job A's seed is irrelevant: n1 uses fixed Sobol seeds.
>   * Job C (J = 400) failed: a float SVD failure in the conserved-structure solver at large J.
> * **Code issues:** hs_odd.py's non-GK edge lacks a sign-free tail term, and its degree-7 representation probe is noise-amplified on direct rows.

> **Update (AG_PROOF.md, unreviewed): (AG′) is NOT proved, and in its literal local-rate form it is false near the gap edge. It is replaced by a strictly weaker pair of premises, (AG_U) + (EDGE-R), that suffices. With them Theorem W's sampled J\* = 26 (no longer 42 under the old premises), and Theorem L′ holds on the uniform rate alone.**
> * **Negatives:**
>   * **N1 [P + DV]:** the flat partial wave itself violates (AG′)'s pointwise bound near the edge (branch point at ν = m; ×4.1 excess at J = 26). The edge integral the lemma actually needs still holds (0.72–0.82 of (AG′)'s). P4_WINDOW's 15/κ-vs-10/κ sensitivity is an artefact.
>   * **G [NE]:** the exact scalar continued to complex ν (Arb, entire series) has an effective edge rate of 0.64–0.86 κ_loc at x₁ = 0.94–0.96, against the assumed 0.97. In the bulk it is 0.98–1.06. The uniform rate holds everywhere with margin (≥ 1.11κ₀; ≥ 2.3κ₀ at the edge).
>   * **N2 [P + DV]:** the plane-wave route cannot reach the near-gap edge. A positive average keeps only the worst rate, and at J ≤ 50, 14–35 % of the measure has |p| < Λ. It is harmless for J ≥ 100.
> * **Proved [P]:**
>   * **Lemma F:** an exact flat-wave split W₊ = (2/πi)C_J e^{i(J+1)θ}₂F₁(½, J+1; J+3/2; e^{2iθ}) with an explicit bound, so the uniform rate 0.97β holds for every explicit piece of the hybrid.
>   * **Lemma P4W-χ:** an entire erfc cutoff splits the complement cost into (EDGE-R), one real-axis integral E_χ[G] = ∫_{p_w}^1 ψ₂χG with no continuation and no rate, and (AG_U), a Cauchy part needing only the uniform rate 0.97Λβ. (AG′)+(EDGE) ⇒ (AG_U)+(EDGE-R).
> * **Constants [DV]:**
>   * The exact scalar edge integral is ≤ 0.30 S_E, 7–100× below (AG′)'s charge.
>   * Theorem W (47 states): everything except the edge integral uses ≤ 0.489 H. With the plain envelope the tolerated edge integral is ≥ 158 S_E, so **sampled J\* = 26**. With the Ω^U +− transfer the tolerance is ≥ 0.406 S_E at (26, 1.026Λ); J\* = 26 if the graviton +− edge is scalar-like, 38 if tolerance ≥ 1.2 is demanded.
>   * Theorem L′: the uniform rate alone suffices (ε_R ≤ 0.133; edge tolerated to 18.9× the scalar envelope).
> * **Remaining premises:**
>   * **(AG_U) beyond all orders [A]:** [D] to all orders; near the edge it only asks for the split at rate κ₀.
>   * **(EDGE-R) for the graviton [M]:** decided by real-axis graviton values near the edge at J ≤ 51.
>   * The exact scalar data stop at x₁ = 0.962.
>   * Unchanged: P0, P1′, W, the boundary factor 2 [D-sketch], the +− transfer convention [M]; sampled states only, float only.

> **Update (DSKETCH_UP.md, unreviewed):**
> * **c₁: universality and ν-independence → [D].**
>   * The regular second-order ASD solution is explicit [DV] (exact sympy at two kinematic points), and exists for all momenta [D].
>   * The remnant is Φ₃₄ = −b₊W₊ with b₊ ≠ 0, which agrees with the JJ21 check.
>   * J₀(E×W³) = 2 and polynomiality are [D]; the Mellin → ν map T ↦ (A²+ν²)/2A is exact by Barnes' first lemma [P+DV].
>   * So "(Pos) ⇒ (1)" is [D] when no light single trace couples to TT.
> * **Light single traces: NEGATIVE finding.** In CHMRS's exact Mellin representation, the per-ν action of a single light block grows like e^{cν} (c → π/2) above the block's mass. **The leading-Regge/hybrid "image" of a light block is not its AdS image: the UT3 §5.5 row "LR form = true AdS image" is refuted as stated.**
>   * Its light value is the arc at infinity: a polynomial in ν² of degree n − k, [D given (E′)]. (E′) is CHMRS's light-side formula carried to the spinning k = 2, 3 functionals, and is [D-sketch].
>   * For J = 0 the AdS effect is m² → μ² = (Δ_L−2)(Δ_L−1) + δ.
>   * Ω^EE on light states is robust: J = 0 positive iff μ² > −0.545Λ².
>   * LS6 (ε = 0.040) is positive iff μ² > −692.5 / −7748.7, against a natural |δ| ≈ 21 [M]. ε\* is expected to fail for Δ_L ≲ 6, so 0.040 is necessary.
> * **Caution (unchecked):** weights reaching ν above a near-gap state's mass (e.g. de-shift tails to Λ+100) may meet the same exponential growth. This needs a re-check.

> **Update (EDGE_R.md, unreviewed): (EDGE-R) for the graviton is NOT closed; it stays [M] at all 30 Theorem W states with J ≤ 51. No violation was found, but no graviton E_χ was computed at any W state. External compute is the route.**
> * **What decides it:** at the gap-mass states the value comes from x ∈ [0.94, x₁]. For the exact scalar the parts above and below ν = 0.94Δ cancel to a net 0.11–0.29 S_E2, so sign-free covers cost 4–10× the tolerance and the oscillating values are needed.
> * **Why n1 fails at ν/Δ ≳ 0.94:** the spin prefactor grows at complex angle and puts a second hump on the real contour (e^{14} above the saddle). Re-centring cannot fix it.
> * **n3:**
>   * scalar: works to the edge (within 1.6σ of the exact Arb scalar up to ν/Δ = 0.970);
>   * graviton k = 2 (EE): agrees with n1 where both work, and gives values beyond 0.94;
>   * **graviton k = 3 (C3H): broken in n3** (state values ×200–1000 too small, cause not found).
> * **Exact scalar at J = 26–38 (new):** 0.285 / 0.201 / 0.108 / 0.151 S_E2 at the states' own cutoffs. **Correction to AG_PROOF: "exact scalar ≤ 0.30 S_E over all cutoffs" is wrong.** At J = 26 the grid maximum is 0.405, essentially at the tolerance of 0.406. At the tabulated cutoff a scalar-like graviton passes with margin 1.42.
> * **Units correction:** the plain tolerance (≥ 158 S_E) and the transferred tolerance (≥ 0.406 S_E,t) are the same number. AG_PROOF's "plain margin ≥ 520" compares with the scalar edge, not the graviton.
> * **External route:** ext/jobs_batch3.txt (1586 n1 + 150 n3 jobs, about 90 CPU-h); ext/jobs_batch3_n1.txt holds the n1 part only. The 20 states with m ≥ 1.13Λ need only n1 and could reach [NE±]. The gap-mass states also need the n3 C3H fix and 2–4× less bulk MC noise.

> **Update (NU_TAIL_CHECK.md, unreviewed): the exponential growth above a state's mass affects NO computed number, but the published grid-weight DEFINITION at Δ_g = 1000 is defective, and the fix is free.**
> * **No effect on computed numbers [P, by inspection]:** every computed per-state number (crystal/D′, K₂ rows, direct n1/n3, GK edge nodes, windows, EDGE-R, hybrid floors) uses ν ≤ Λ ≤ m. The light-side constants are unaffected, and so is everything built from the continuum weights or the exact min-norm X: X and X^Q decay like e^{−(π/2)(ν−Λ)} beyond Λ, and tails beyond m are ≤ 4.4·10⁻⁷ (D) and ≤ 1.6·10⁻⁴ (Q) of a state's action.
> * **Problem (c):** the stored grids `prl/data/d_abel_X_1000.txt` and `prl/data/k2/xq_abel_X_1000.txt` run to Λ+312 and carry a non-decaying double-precision floor (|X| ≈ 10⁻¹⁵–10⁻¹⁶) beyond Λ. Read literally with θ(Λ+100, 10), that floor times the exact image growth (Γ ≈ 10¹⁶ at Λ+100) gives tails up to 650× (D) and 1300× (Q) a near-gap state's action. The sign is set by round-off, so (Pos) for that grid-defined object is untested at m − Λ ≲ 60–100. HK_CLOSE §5 and §6 claims are false for that object, though true for the computed values. Δ_g = 300 grids are unaffected (they end at Λ+29).
> * **Fix (changes no number):** define the functional by the continuum weights; or read the grids as zero beyond Λ+10; or use the cutoff (Λ+60, 6) for the k = 2 pieces (same leakage, tail 5·10⁻⁹ / 4.6·10⁻⁷).
> * **Growth numbers [DV]:** exact scalar Arb series. Γ(m) = 1.8·10³ and Γ(Λ+100) = 1.6·10¹⁶ at Δ_g = 1000, J = 0; the rate per unit ν is 0.13–0.19 at x = 1, rising to π/2 only as x → ∞.

## 0. Grades

| grade | meaning |
|---|---|
| [P] | proved |
| [certified] | computer-assisted proof: every inequality checked in exact rational or ball (interval) arithmetic, on top of written-out lemmas |
| [D] | derived (argument written out, not independently checked) |
| [D given X] | derived from the named inputs X; net grade = weakest of [D] and X |
| [D-sketch] | derivation sketched; at least one step asserted from a standard result without being written out or checked |
| [DV], [D + DV], [P + DV] | derived (or proved) and checked numerically/symbolically |
| [Q] | quoted from the literature |
| [NE±] | numerical evidence with an explicit error budget: every model input carries a stated safety factor and the conclusion survives it (PE §0) |
| [NE] | numerical evidence |
| [NE + M] | computed inside a model (near-gap chain / hybrid); the model's accuracy is itself a hypothesis |
| [M] | model |
| [A] | assumption |
| [I] | open |

**Weakest-link order:** [P] = [certified] > [DV] > [D] > [D-sketch] > [NE±] > [NE] > [NE + M] > [M] > [A] > [I].

**Three layers.** (1) Setting (A1′, A2). (2) Inputs of the implication (§4). (3) Per-state positivity (Pos), which the theorem
*assumes*; §5 lists the *evidence*. Evidence is never an input.

---

> **Update (INDEP_CHECK.md, independent reproduction from the statements with new code; no project code imported or read): the H_K EE reduction is reproduced, and Lemma Y is reproduced for the count and at the correlator level.**
> * **EE kernel = Σ_{b=−4..4} shifted scalar waves:** holds to 2.6·10⁻¹² at 81 cases; end to end to 10⁻⁹–5·10⁻⁸. **Corner exponents** (w2g −0.50, w2s ≈ −0.4, w4g −0.50, w4p +1.00) and the **thresholds** J₀ < 3/2 (k = 2, 3, 5) and J₀ < 5/2 (k = 4) are re-derived. The w4g exponent stays at −0.500 out to β = 13 in 45-digit arithmetic, so the c = ½ input is conservative. **f_pos(±3i/2) = 0** for Q's published k = 4 continuum weights (contour ≤ 10⁻⁶⁸; decay-rate test against a control with a pole).
> * **Lemma Y:** the count of 997,836 = 35,637 × 28 is reproduced, worst margin 0. The closed form σ = (g+e+f) − (z+e′+f′) + (k₁₃−k₂₄)/2 makes it [P] for any P·P content. Five explicit bulk spin-2 contacts (833 polynomial × D-function terms each) behave as A + B log y down to y = 1.2·10⁻⁶.
> * **Not independently reproduced:** the step from the correlator to the light-function bound for L6′ (it stays [D] on the original work); the non-dispersive counterexample slopes; HK_RESIDUAL §1's k = 4 kernel normalisation; K1′, K2ν and LRC-a.
> * **Minor corrections:** the failure of EESTEP §1.3's negative-b series is truncation, not round-off (not load-bearing); HK_CLOSE and RESULT_Q_STAGE3 differ by one factor of η in the raw k = 2 corner entry (both say "fails"); one displayed line in R4H §2.3 is garbled, and the closed form replaces it.

## 1. Setting, normalisations, couplings

| item | definition | source |
|---|---|---|
| theory `𝒯` | unitary CFT₃, parity-even or parity-violating, unique stress tensor `T` | UT2 §1 |
| `c_T` | free real scalar has `c_T = 1`; `8πG := 256/c_T` (AdS units), from `c_T = 32L²/(πG)`; independently confirmed | UT2 §1; REPRO §1 |
| **A1′ (higher-spin gap)** [A] | every single-trace primary in `T × T` with spin `J ≥ 4` has twist `Δ − J ≥ Δ_g`. Single traces of spin 0 and 2 (either parity) may have any dimension. Spins 1 and 3 do not occur in `T × T` (no `⟨TTO_J⟩` structure for identical conserved `T` at J = 1, 3 in d = 3). With A1′, `Λ` is the spin-4 gap (the CHLPS/CEMZ notion). | A1P §0, §2.1, §6 |
| A1 (old, full gap) | as A1′ but for every J including 0, 2. **No longer assumed.** A1 ⇒ A1′, so every statement below holds a fortiori under A1, where the light-state links of §4.7 are absent. | UT2 §1 |
| **A2** [A] | tree level, `O(1/c_T)`, `c_T ≫ Δ_g²`; bounds refer to `c_T → ∞` at fixed `Δ_g` | UT2 §1 |
| gap scale | `Λ := √((Δ_g−2)(Δ_g−1))`: `298.4996` (300), `998.4999` (1000) | UT2 §1 |
| heavy kinematics | `m² = (Δ−J−2)(Δ+J−1)`; `cosh β = 1 + (2J+1)²/(2m²)`; `b̂ = (2J+1)Λ/m` | UT2 §1 |
| cubic couplings | `⟨TTT⟩ = (c_T, t₄, α_T)`; `R := √(t₄²+α_T²) = 240|ĝ₃|` (CHLPS normalisation, AdS units) | UT2 §1 |
| `g₄` | `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴f(s,u)`, `f = 8πG/(stu) + … + g₄ + g₅t + …`; type II: `g₄M⁶/8πG = 2ζ(3)`, `M² = 4/α′`. **Under A1′, `g₄` is the R⁴ coefficient with the light fields kept in the amplitude (minimal light-exchange completion, CHLPS App. B convention)** [D] | UT2 §1; A1P §2.2, §6 |
| not bounded | all-plus `ĝ₄` (value 0 under the k = 4 functional [D]); `g₅` | F_NORM §2.4 |
| `g̃₄` | `g̃₄ := −½𝒜₄[light] = −½Ĉ₄^{ds}|_light` | F_NORM §1 |
| `c₁` | coefficient of the `O(t₄)` light term `16πGc₁t₄` of `Ĉ′_EE` (and `Ĉ₃`). One universal constant fixed by `⟨TT⟩` and `⟨TTT⟩` data (`c_T`, `t₄`) alone [D-sketch] (weakest link: the regular second-order ASD solution, L3d); uncomputed; **zero given premise F₂ [A]** (§4.3). Not a free parameter. | C1V §0, §4 |

---

## 2. The functionals (cut-off, defined by published weights)

### 2.1 Cutoff
`θ(ν) = ½[erf((ν+ν_c)/a) − erf((ν−ν_c)/a)]`; every de-shift weight `χ → χθ`.
* `Ω^EE_θ` and **both k = 2 pieces** (D's and Q's smeared one): `ν_c = Λ + 100`, `a = 10`. (CERT §3.2: `X^Q` extends to ν ≈ Λ, so Q's
  `ν_c = 150` must not be used for it. Correction to UT2 §2.1.)
* Q's C4, C5 pieces: `ν_c = 150`, `a = 6`.
* Relative effect on any computed number: `≲ e^{−180}` (Q C4/C5), `≲ 10⁻⁴¹` (k = 2) [D + NE] (HK §0). Required: a weight of compact
  ν-support regrows at the raw Regge rate (Paley–Wiener) [P] (HK §3).

### 2.2 `Ω^EE_θ` (Result D)
`Ω^EE_θ = ∫₀¹dp ψ₂Λ⁴Ĉ′_{EE,Λp} + ∫₀¹dp ψ₃½Λ⁶Ĉ_{3,Λp}`, NF02 weights (`nf02_psi.json`, exact rationals, `k(0) = 1`), k = 2 part with
dual structure `EE`, min-norm Abel de-shift `X`; k = 3 part `Ĉ₃` (`j′ = 1`, no de-shift). Full definition UT2 §2.2.
* **Weight realisations (CERT §3.2).** (a) the exact min-norm `X` (existence [D + NE]), tabulated on a grid in `d_abel_X_{300,1000}.txt`
  (identity holds only to `t ≈ 22–25`); (b) **new:** the explicit continuum weight `k2_continuum_D_{300,1000}.json` (`prl/data/k2/`),
  `χ = Q(ν)B(ν) + Σ_{j≤2} y_jν^{2j}(ν²+1)log(ν²+1)e^{−ν²/2}`, `Q = ν²/(ν²+¼)²` (double poles at `±i/2`, branch points at `±i`), with
  `|W − Φ| ≤ ε_W` **for all t ≥ 0**: `ε_W = 3.3·10⁻⁶ / 1.1·10⁻⁶` (`1.1·10⁻⁹ / 1.1·10⁻¹⁰` of `∫|ψ|`) [DV]. Equal to the grid weight to
  `≤ 4·10⁻⁷` of max ψ.
* Both give the same `K = k(0)(1+δk(0))`: `(1+dk₀)^{−½} = 1.00038354 / 1.00003415`, so `R₀` is unchanged (CERT §0) [DV].

### 2.3 `Ω^Q_up` (Result Q, upper end)
* **Flat target Q3\*** (definitional fix, CERT §2.1): the δ = 0.2 functional `Q3-full-delta0.2` **with its weights in factored form**
  `ψ₂ = p(1−p)³X(p)`, `ψ₃ = (1−p)³Y(p)`, `X, Y` the exact quotients of the stored doubles (dropped remainders `≤ 1.7·10⁻⁸`; fingerprint
  `sha256 80e0834a…33f2`). The stored binary weights carry `(1−p)³` only to rounding, and their flat action is **negative** at `+−`,
  `m = M`, `J ≥ 1990` [NE]; Q3\* removes this (min `+0.2000` there). `A = 2ψ₂′(0) = 17.240713`, `κ = 0.109908`, `a = 0.169405`, `c = ½`
  unchanged. Nominal LP margin δ = 0.2; **true flat minimum 0.1835** (`++`, `J ≈ 3600`, `m ≈ 100M`, `b̂ ≈ 72`) [NE], so δ = 0.2 is not
  met (by 8 %). (The Q3\* coefficients are not yet a separate published file.)
* **C4 value / derivative:** continuum weights `q_deshift_continuum_c4{v,d}.json`, `χ = Q·B` (poles at `±3i`, `±i/2`, `±5i/2`,
  `±i(3+2√2)/2` only). `ε_W ≤ 2.673·10⁻¹³` (value), `≤ 9.058·10⁻¹¹` (derivative) **[certified]** (HKR §3.3; DESHIFT's midpoint values
  were 1.2·10⁻¹³ / 4.2·10⁻¹¹).
* **C5:** `q_deshift_continuum_c5.json`, `ε_W = 2.4·10⁻⁷` [DV] (not certified; not load-bearing for H_K).
* **Smeared k = 2, 3 parts** with Q's `ψ₂, ψ₃`; the k = 2 de-shift weight `X^Q` is **now published** (`prl/data/k2/xq_abel_X_*.txt`,
  re-solved independently to `≤ 8·10⁻⁷` of max ψ₂^Q [DV]), and also as a continuum weight `k2_continuum_Q_*.json`
  (`ε_W = 1.9·10⁻⁷ / 7.4·10⁻⁸`) [DV] (CERT §3).

### 2.4 `Ω^Q_low`
The de-shifted k = 4 value piece alone (C4-value weight). Leading-Regge heavy image `P_{−½}(cosh β) > 0`.

### 2.5 Per-state positivity
For a primary `O` in the relevant set, `M_Ω(O)` is the Hermitian action on its blocks: 2×2 on (`++`, `+−`) for even J ≥ 4; 1×1 parity-odd
(OO) for parity-odd even J; **1×1 `+−` (parity-odd) for odd J ≥ 5**; 1×1 `++` and 1×1 OO for J = 0, 2 (no `+−` structure below J = 4)
(A1P §2.1; NG2 §1.1). **(Pos_Ω):** `M_Ω(O) ⪰ 0` for every heavy primary (`Δ − J ≥ Δ_g`) and, **under A1′, for every light J ≤ 2 single
trace in `T × T`** that enters the sum rule at that k (§4.7). Values: `λ/N_s` (D), `FF/N` (Q).

---

## 3. The theorem

> **Theorem (conditional).** Let `𝒯` satisfy A1′ and A2 at `Δ_g ∈ {300, 1000}`. Assume the inputs of §4 and (Pos) for `Ω^EE_θ`
> (item 1) and for `Ω^Q_up`, `Ω^Q_low` (item 2). Let `c₁` be the constant of §1. Write `R = √(t₄²+α_T²)`,
> `L_A = ∫₀¹ψ₂dp/(p²+9/(4Λ²))`, `K = k(0)(1+δ₃)`, `δ₃ = δk(0) + (H17 remainder)`.
>
> **(1) Cubic.** For every real `c₁`: `R ≤ R₀(x + √(1+x²))`, `R₀ = 240√(16L_A/K)/Λ⁴`, `x = 19941.0|c₁|/(Λ²√(16KL_A))`.
> At `K = 1`: `R₀ = 3.13206·10⁻⁶` (300), `2.89126·10⁻⁸` (1000); `x/|c₁| = 2.16006·10⁻³ / 1.67026·10⁻⁴`. Under A1′ the sum rule in
> addition bounds the light couplings: `K R²Λ⁸/57600 + Σ_light F_D^{(J_ℓ)}w_ℓ ≤ 16L_A + 16I_DΛ²c₁t₄`, every `F_D ≥ 0.448` (flat, §4.7).
> The constants are unchanged from A1.
>
> **(2) Quartic.** Add, if `(t₄, α_T) ≠ 0`, the signs `κ > 0` (upper) and `c_L > 0` (lower).
> * **(2a) Published constants.** If no light scalar with flat mass `m < 0.4602Λ` couples to `TT` (under A1 this is automatic):
>
>       −O(Λ⁻²) ≤ g̃₄Λ⁶/8πG ≤ 2L_AdS + 57600I_Q²c₁²/(κΛ⁴) = 72.041 + 6.8145·10⁻⁵c₁² (300),  92.814 + 5.4428·10⁻⁷c₁² (1000).
>
>   More generally the upper end is `2L_AdS + Σ_{light scalars} max(0, 0.430430 − 2.032073m_ℓ²/Λ²)|g₀(ℓ)|²Λ⁴/8πG + (c₁ term)` (A1P §4.3).
> * **(2b) Pure higher-spin gap (no coupling condition).** With `Ω^Q_up + εΩ^EE_θ`, **ε = 0.040**, and (Pos) also for `Ω^EE_θ`:
>
>       g̃₄Λ⁶/8πG ≤ 2L_AdS + 2εL_A + 57600(I_Q + εI_D)²c₁²/((κ + ε/8)Λ⁴) = 125.71 + 1.294·10⁻⁴c₁² (300),  164.51 + 1.033·10⁻⁶c₁² (1000).
>
>   (`c₁²` coefficients computed here with A1P's formula, flat κ and K = 1; at ε* = 0.038034 the same formula reproduces A1P's
>   1.2598·10⁻⁴ / 1.0062·10⁻⁶.) The lower end is unchanged in form.
> * If `t₄ = α_T = 0`, the `c₁` and κ terms are absent and no sign input is needed. For `R ≠ 0` the lower end reads
>   `g̃₄Λ⁶/8πG ≥ R²Λ⁸/(57600·4c_L) − O(Λ⁻²)`. **The lower end has no computed constant.**
>
> **(2′)** Given H_EFT(b) [A], `g̃₄ = g₄ + O(Λ⁻²)·(EFT data)`.

**Net grades (weakest link).**

| statement | net grade | weakest links |
|---|---|---|
| "(Pos) ⇒ (1)" | **[D-sketch]** | the `c₁` form (ν-independence, [D given `J₀(E×W³) = 2`], D-sketch); under A1′ the link "light value at k = 2, 3 = block action" [D-sketch]; H17 band [D-sketch + M] (only the ±0.24 % / ±0.021 % band). H_K is now [D]. |
| "(Pos) ⇒ (2), upper end" | **[D-sketch]**, plus **[A]** (H_EFT(b)) for (2′) | (P) for `g₅` [D-sketch + NE] (the upper end holds under either reading); under A1′ the k = 2, 3 light link [D-sketch] and, for (2b), the ε = 0.040 repair [D-sketch]; κ > 0 [D-sketch + M] if `R ≠ 0`. **No [I] link** (L6′ is now [D]). |
| "(Pos) ⇒ (2), lower end" | as upper, **no constant**; `c_L > 0` [D-sketch + M] if `R ≠ 0` | |
| (Pos) itself | **not established** | §5 |
| constants (`R₀`, `2L_AdS`, `x`, Λ-forms) | [DV], reproduced independently | REPRO §1; DG §1 |

**Dependence.** (1) uses only `Ω^EE_θ`. (2a) uses only `Ω^Q_up`, `Ω^Q_low` (the exact `c₁` form needs only κ > 0). **(2b) also needs
(Pos) for `Ω^EE_θ`**, so in the pure-higher-spin-gap form item (2) is no longer independent of item (1)'s positivity hypothesis (A1P §4.3).

**Sensitivity to `c₁`** (C1V §6; (2b) computed here). The bound moves by `< 1 %` iff:

| | 300 | 1000 |
|---|---|---|
| (1) | `|c₁| < 4.61` | `|c₁| < 59.6` |
| (2a) | `|c₁| < 103` | `|c₁| < 1306` |
| (2b) | `|c₁| < 99` | `|c₁| < 1262` |

Exact sample values of (1) and (2a) for `|c₁| = 0 … 1000`: UT2 §3 (unchanged). (The earlier natural-size estimate `10⁻³–10⁻¹` [M] is withdrawn: C1V shows c₁ is not a curvature remnant of the seagull, so the estimate had no basis.) Under F₂, `c₁ = 0` and
the headlines are exact as printed.

---

## 4. Inputs of the implication: full grade table

"D" = item (1), "Q" = item (2). Each row: statement, net grade, source. Changes since UT2 in **bold**.

### 4.1 Regge input and the interchange H_K

| # | statement | D | Q | grade | source |
|---|---|---|---|---|---|
| R1 | KKSDZ `J₀ ≤ 1` for `⟨TTTT⟩` (1904.05905 §4.2.4, eq. 4.55), caveat §4.2.6 | ✓ | ✓ | [Q] | UT2 §4 |
| R1′ (K1′) | uniform Rindler bound on the Regge box `B_c`, `|dDisc 𝒢| ≤ 2B(c)`, edge exponent 0; the caveat lies outside `B_c` | ✓ | ✓ | **[D + DV]**: spinning check on the exact free-scalar `⟨TTTT⟩` (63 polarisation sets, box to 10⁻⁶, edge row with τ_min = 1); interacting spinning correlator not exercised | LEMMAS §2; **HKR §2** |
| R2 | K2 power counting; K2ν `|𝒦_ν| ≤ K̄` | ✓ | ✓ | K2 [Q]/[D]; K2ν [D] | LEMMAS §3.1 |
| R3 | per-block Rindler positivity ⇒ swap lemma | ✓ | ✓ | [D] | HK §1 |
| DS (D, k = 2) | min-norm `X` exists, `χR_j ∈ L¹`, strip-analytic | ✓ | ✓ (k = 2 piece) | [D + NE]; **alternatively ε-DS for the published continuum k = 2 weight [DV]** | HK §4.1; **CERT §3.2** |
| DS (Q, C4/C5) | **not an H_K input** for the explicit continuum weights: the analyticity HK derived from DS holds **by construction [P]** | — | ✓ | **[P] (analytic structure)**; ε-DS is an input of the *values*: C4 **[D + certified]**, C5 [DV] | **HKR §0, §3** |
| LRC-a | leading-Regge form of the per-ν kernel in the Regge corner | ✓ | ✓ | scalar, EE k = 2 [D + DV]; **EE k = 4 [D + DV]** (EES kernel = Stage-4 split kernel to 1.8·10⁻⁷; numerical rank-3 heavy image = `R_j^{EE3}` to 1.1·10⁻¹¹) | LEMMAS §3.2; EES §1.4; **HKR §1** |
| LRC\* | exact per-ν kernel = continuous-spin partial wave; edge corrections resum to `₂F₁`; de-shift cancellation survives to all orders | ✓ | ✓ | scalar [D + DV]; **EE [D + DV]**: the EE kernel is an exact finite sum of 9 b-shifted scalar continuous-spin partial waves with ν-independent coefficients ([P] reduction; end to end 1.6·10⁻⁷–1.5·10⁻⁶); exponents of (exact − LR) −0.50 (k = 2), −0.50 / +1.00 (k = 4 regular / pole at 3i/2) | LRC §0–§4; **EES §1–§2** |
| f_pos(±3i/2) | position-space edge coefficient of Q's C4 weights at `±3i/2` | — | ✓ | **= 0**: [P] for the continuum weights (no pole at ±3i/2 by construction); [D given DS] + [NE] for the Stage-7 grid weights; ν-normalisation [DV] | **EES §3; HKR §1.4** |
| **H_K** | for each cut-off piece, `Σ_O Tr[P_O Ĉ_θ[G_O]] = −Ĉ_θ[light]`, absolutely convergent | ✓ | ✓ | **[D]** for D and for Q (Ω^Q_up defined by the continuum weights). If Ω^Q_up is instead taken as the Stage-7 grid weights: [D given DS], DS [NE]. | **EES §0, §6.2; HKR §0** |
| H_K criteria | `J₀ < 3/2` (k = 2, 3, 5); **`J₀ < 5/2` (k = 4 with Q's weights)**; `J₀ < 2` for generic k = 4 weights | ✓ | ✓ | [D]; under R1 the margins are **½** and **3/2** | **EES §0; HKR §1.4** |

### 4.2 Light side of `Ω^EE_θ` (item 1)

| # | statement | D | Q | grade | source |
|---|---|---|---|---|---|
| L1 | Adm; L0 (double zeros, contacts drop); H7 pole `16πG/(ν²+9/4)`, `Ĉ₃|_E = 0`; H7 = the Regge arc at n = k = 2 | ✓ | (smeared) | [D] | UT2 §4 |
| L2 | `K = k(0)(1+δ₃)`, `δk(0) = −7.666·10⁻⁴ / −6.830·10⁻⁵` (reproduced by the continuum weight); H17 band `±212.5/Λ²` | ✓ | (κ) | δk(0) [D + DV]; H17 [D-sketch + M] | UT2 §4; CERT §3.2 |

### 4.3 The constant `c₁`

| # | statement | grade | source |
|---|---|---|---|
| L3a | the `O(t₄)` light term is `16πGc₁t₄`, ν- and Δ_g-independent | [D given `J₀(E×W³) = 2` (D-sketch) and polynomiality] | C1 §0 |
| L3b | `u`-channel (shock) part = 0; **in d = 3 the preserved elements of the CHP impact-parameter operators vanish identically on every H₂ eigenfunction (all ν, L), so the shock route gives exactly 0 and cannot see a ν-independent term at all** | [D]; identity **[P]**; no-go [D] | C1 §3; **C1V §2** |
| L3c | k = 3 part = 0 | [P given Q2, Q8, F] | C1 |
| L3d | **ASD closure: the whole `O(ĝ₃)` preserved⊗preserved tree (s-, t-channel `E×W³` exchanges + seagull) collapses to an energy-channel `T` exchange `⟨T⁺T⁺T⁺⟩_{W₊}·[⟨T⁺T⁻T⁻⟩_E − κCS]/⟨TT⟩`; no independent seagull; so `c₁` is fixed by `⟨TT⟩`, `⟨TTT⟩` data alone (universal)** | [D-sketch] by the weakest-link rule (one input, a regular second-order ASD solution, is [D-sketch]); boundary relation checked against Jain–John–Mehta–Nizami–Suresh 2106.00016 [Q + DV] | **C1V §3** |
| L3e | **`c₁ = 0` given F₂** (the k = 2 functional sees the s-channel `T` family only through light-ray one-point matrices: diagonal for Einstein, pure helicity flip at `O(t₄)`) | `c₁ = 0` [D-sketch given F₂]; **F₂ [A]** (k = 3 analogue [D]); unconditional value [I]; no rigorous bound | **C1V §4** |

### 4.4 Light side of `Ω^Q` (item 2)

| # | statement | grade | source |
|---|---|---|---|
| L4 | one functional: `Ω_k[X] = D_LR[X] + 𝒜_k[X]`; `𝒜_k = 0` for Regge spin `n_X < k`; at `n_X = k` fixed by the top Regge coefficient | [D given C′ + Mellin convergence]; [DV] on a Witten family; **borderline mechanism checked at k = 4 on a contact, a crossed and a direct exchange (block included) to 5·10⁻⁷ – 1.4·10⁻⁶ [DV]** | LIGHTSIDE §0–§1; **A1P §5.2** |
| L5 | Lemma C′; Lemma A (`D_LR = 0` on tree terms with `J₀ < 4`) | [D]; Einstein instance derived, not evaluated | LEMMAS §1 |
| L6 (F) | `−𝒜₄[R⁴] = −Ĉ₄^{ds}|_{R⁴} = 2g₄`; `ĝ₄ → 0` | **[D]** (inputs Lemma A, H_EFT(a) [now D], CLS-D J = 0 heavy dictionary) | F_NORM §2.3; R4H §5 |
| **L6′** | the borderline R⁴ light function is integrable against `C_ν` at the collinear set `a = b`: `|L_{R⁴}| ≤ CΔφ⁻²(1+|log Δφ|)`, i.e. p = 2 against the criterion p < 7/2, margin 3/2; same bound for every contact and every heavy-exchange double-trace part | **[I] → [D]**, given (R1, R2) the polynomial × D-function form of spinning contact diagrams [Q: Paulos 1107.1504; Costa–Gonçalves–Penedones 1404.5625] and the Symanzik/Mellin representation [Q: Penedones 1011.1485]; counting **Lemma Y [P]** (exhaustive over 997 836 structures, worst margin 0); premise **dispersiveness** (finitely many non-zero s-channel residues) [D] for integer Mellin shifts, CLS-D's statement [Q]; every holographic `⟨TTTT⟩` tree term qualifies. Residual: the AdS₄ R⁴ `⟨TTTT⟩` itself was not evaluated; the physical Δφ scan is not yet asymptotic (slopes 2.5–2.8 at Δφ = 0.04) | **R4H §0, §2** |
| L7(a) | H_EFT(a): term-wise light value of a heavy exchange's double traces = actual value | **[D] for dispersive terms, pointwise in angle (hence ν-paired)**; [DV] on the dispersive spin-2 sector (`J_h = 0, 2, 4`, `τ_h = 41–81`, k = 2, 4: 10⁻¹⁸–10⁻³⁷). An unexplained ≈ −3 % leading-order ratio for `J_h = 2` stays [NE] (route II needs only a scalar χ) | **R4H §3** |
| L7(b) | H_EFT(b): light double-trace data = tree-level bulk EFT, no conspiracy in heavy sum / contact tower | **[A]**, identical to CHMRS 2106.10274 §5.1's standing assumption and CLS-D §2.2 | R4H §4 |
| L8 | (P): R⁴ per-ν light value constant [D]; `g₅`: **exact cancellation excluded** — the `++` derivative-piece correction is m-independent from m = Λ to 20Λ at J = 0 (direct, `s_d^∞ = 0.5081 / 0.5095 / 0.5096 / 0.5095` at m/Λ = 1 / 3.35 / 10 / 20) and at J = 2 (analytic-radial, 0.5102 / 0.5085 / 0.5095 / 0.5096); the double-counting reading holds (lower slope 0.255); C5 piece not run [I]. The upper end holds under either reading. | R⁴ [D]; **`g₅` [D-sketch + NE]** | **EES §5; HKR §4** |
| L9 | κ > 0, c_L > 0 (only if `R ≠ 0`) | [D-sketch + M] | UT2 §4 |

### 4.5 Q's de-shift values
| # | statement | grade | source |
|---|---|---|---|
| ε-DS | `|W[u_p] − Φ_p| ≤ ε_p` for all t; per-state images and light moments = targets ± ε; every bound moves by ≤ 6·10⁻¹¹ (margins) and ≤ 10⁻¹² (C4) / < 10⁻⁶ (C5) (light side) | C4 **[D + certified]**; C5 [DV] | HKR §3.2–§3.3; DESHIFT §2 |

### 4.6 Positivity-evidence premise used everywhere the 1/Δ expansion is used

| # | statement | grade | source |
|---|---|---|---|
| P_cyc | "the physical cycle lies in the class of D′" | **an identity: CLS-D define the functional (2311.04271, Lorentzian shadow representation) as an absolutely convergent integral over real Lorentzian points in the causal order 4 > 5 > 3 > 0 > 1 > 2; in the project gauge this region is exactly D′** (region [P] + [DV], 0 mismatches in 2·10⁵ configurations; integrand [D] + [DV]) | **PCYC §0–§1** |
| RD on D′ (Prop. D′) | engine = `∫_{D′}` + `O(Δ^N e^{−0.336Δ})` | [D], one [D-sketch] sub-step (convergence radius of the x₀-dependent slow factors) | RD_CLOSE; PCYC §0 |
| **H_RS** | engine image = physical integral at Q's binding states | **[D] with one [D-sketch] sub-step** (net [D-sketch] by the weakest-link rule). Numerical: direct `∫_{D′}` = engine K = 5 to 3·10⁻⁷ [DV]; scalar end to end at Δ = 20, 30, 50 (where 3–82 % of `∫|Φ|` lies outside the convergence region of the 1/Δ series): ≤ 2·10⁻⁸ / 6·10⁻⁹ / 1·10⁻⁹ [DV] | **PCYC §0, §1.4** |
| crystal cycle | saddle-point (crystal) cycle = D′ | **only a numerical method now, not a premise**: [P] in the Regge-scaled system; [I] at finite Δ (Stokes question); direct EE test consistent (`|T|/σ ≤ 1.7`) but weak [NE, weak] | PCYC §2–§3 |

### 4.7 Light single traces under A1′ (new)

| # | statement | D | Q | grade | source |
|---|---|---|---|---|---|
| LS1 | Regge spins in the helicity-preserving (MHV) sector the functionals see: light scalar (either parity) n = 3; light massive spin 2 n = 4; its crossed-channel n = 2 pole is double-flip and annihilated by `EE` (≤ 5·10⁻⁴⁶) and absent from C4/C5 | ✓ | ✓ | flat [Q CHLPS] + [D]; AdS top-Mellin-degree count [D-sketch] | A1P §2 |
| LS2 | drop-out: light J = 0 at k = 4, 5; light J = 2 at k = 5 (Lemma A on the light exchange, block included) | — | ✓ | [D]; toy [DV] (≤ 1.3·10⁻⁵ of scale) | A1P §3, §5.1 |
| LS3 | light J = 2 at k = 4: borderline arc ∝ top coefficient; flat `12|g₂|²p²/m⁴`: 0 in the value piece, favourable `+12a|g₂|²/m⁴` in the derivative piece | — | ✓ | structure [D + DV]; flat value [Q]; AdS normalisation [D-sketch] (only the favourable sign is used) | A1P §3, §4.3 |
| LS4 | at k = 2, 3 a light J ≤ 2 state contributes its block action (no drop-out; no preserved u-channel pole for J < 4) | ✓ | ✓ | **[D-sketch]** (flat [P]) | A1P §3 |
| LS5 | flat signs: NF02 positive on every light J = 0 (`11.317 + 20.772m²`) and J = 2 (min 0.448 at m = 0.675Λ) at every `0 < m ≤ Λ`; Q3 positive on light J = 2 (min 1.97), **negative on light J = 0 for m < 0.4602Λ** (`−0.4304 + 2.0321m²`) | ✓ | ✓ | **[P]** (exact polynomials, root isolation) | A1P §4.2–§4.3 |
| LS6 | AdS light-state action = flat value × (1 + O(Δ_L²/Λ²)) for `Δ_L ≪ Λ` | ✓ | ✓ | [D-sketch] | A1P §4.2 |
| LS7 | ε-repair: ε* = 0.038034 leaves zero flat margin at m → 0 (sign for `Δ_L = O(1)` set by uncontrolled O(Λ⁻²) terms: **[I]**); **ε = 0.040** gives flat margin 0.022N at m → 0, ×100 above every computed curvature effect | — | ✓ (2b) | ε = 0.040: [D-sketch] | NG2 §0, row "Repaired Q" |
| LS8 | KK towers: finitely many light modes below Δ_g, each J ≤ 2; favourable sums bounded by monotonicity; H_K covers light blocks | ✓ | ✓ | [D]; accumulating spectra not covered [I] | A1P §6 |

---

## 5. Per-state positivity: evidence (not an input)

**Disclosures (unchanged from UT2 §5).** `EE` and the Abel gauge were fixed after heavy-state failures had been seen; `Ω^EE` fixed since
20 Nov, `Ω^Q` since 22 Nov; `J = 18, 20, 24, 28` at 300 are in-sample; in the canonical gauge J = 24 at 300 fails. The statement is about the
min-norm functional.

### 5.1 Flat limit
| item | statement | grade | source |
|---|---|---|---|
| NF02 (D) | **`𝔉 ≥ N/5` on every flat heavy state** (both channels, all J, all m ≥ M); the corner `+−`, J ≥ 400, M ≤ m < 1.633M closed by Lemma B; sharp: `0.2000000103` at `+−`, J = 5, m = M | **[certified]** (new corner code unreviewed) | CERT §1 |
| Q3\* (Q) | `FF ≥ (9/50)N` on: all J < 400; J ≥ 400 with b̂ ≤ 120 (contains the global minimum); `++`, J ≥ 400, m ≤ 2M; J ≥ 400, b̂ ≥ 120, 2M ≤ m ≤ 10M. Not certified: `+−`, J ≥ 400, m ≤ 2M (float min 0.997) and J ≥ 400, b̂ ≥ 120, m ≥ 10M (float min 0.304). Flat minimum 0.1835 | [certified] on the listed regions; [NE] elsewhere | CERT §2 |
| bulk-point map | `Ω^EE` → stripped NF02 (`EE ≡ q²` on CLS-D's plane [P]; full contraction → 1 only as m → ∞ [NE]) | [P given map] + [NE] | UT2 §5.2 |
| Lemma C (curvature) | `|P_{−½+iν}(cosh β) − J₀(νB)| ≤ (1 − sech(β/2)) + min(2, ν(B − β))`, `B = 2sinh(β/2)`; bounds the hybrid's AdS curvature with `νB = b̂p` exactly | **[P]** (given Mehler–Dirichlet [Q]) + [DV] (worst ratio 0.73) | PE §2.1 |

### 5.2 Result D at Δ_g = 1000

| region / state | value | grade | source |
|---|---|---|---|
| **binding state** J = 5, `+−` (1×1), twist edge 1.0049925Λ | own-mass K₂ row: `λ/N_s = +0.20088`; sign-free budget (×3 on model terms) 2.6·10⁻³ ⇒ **`≥ +0.198`, L = 78** (130 with the λ-level truncation); certified floor − Lemma C alone ≥ +0.19993 N_s; 9 direct RQMC nodes agree (λ moves < 10⁻⁶); parity-odd engine agrees on ν = 12.5–38 | **[NE±]** | NG2 §0 |
| near gap J = 0, 2 | ≈ +1.0 | [NE±] (L ≥ 6) | PE §3.1 |
| near gap J = 4, 6 | `+0.20413`, `+0.22953` (own rows) | [NE±], sign-free L ≥ 22 / 5.0 | NG2 §0 |
| near gap J = 8 | `+0.42973` | [NE±] on the λ-level budget only (L = 4.3–8.4); sign-free marginal (1.1–1.8) | NG2 §0 |
| near gap J = 12, 16 | J = 12: `+0.84968` (twist edge), `+0.68574` (1.05Λ); J = 16: `+0.93421`, `+0.76316`, `+0.74581` | [NE + M] (budget does not close, L ≈ 1–3) **⟦HIGHSPIN⟧** | NG2 §0 |
| J = 15 twist edge | **`+0.6900`** (own row; was ≈ +0.69 [M]) | [NE + M] | NG2 §0 |
| other computed near-gap J (7, 9–11, 13, 14, 17, 31, 32, 36, 40, even 18–28, odd 19–29, 48) | min of these +0.2684 (J = 7) | [NE + M] (H_M) **⟦HIGHSPIN⟧ for J ≥ 12** | RD-S §5; PE §3.2 |
| light-cone edge window | computed: J = 18 (−0.0066), 19 (+0.0041), 24 (+0.024), 25 (+0.032) [NE + M]; **J = 29 not resolved** (plain w2 chain −3.62 from one K₂ outlier; with direct nodes +1.39 ± 0.03) [M]; not run: even 20, 22, 26–28, 30–32, odd 21, 23, 27, J = 48 [M] | as stated **⟦HIGHSPIN⟧ (J = 29)** | NG2 §0 |
| even J = 26, 30, 34, 38, 42–46 near gap | interpolation only | [M] (inside H_LJ(1000)) **⟦HIGHSPIN⟧** | PE §3.2 |
| near gap even J ≥ 50, odd J ≥ 33, m ≤ 1.3Λ | κ-model ≥ +0.90 | **H_LJ(1000): [D-sketch + M]**, leverage 2.7–4.8 (refined window; sign-free form violated at (1100, 100): −2154 ± 43 vs T_B = 1409) **⟦HIGHSPIN⟧** | HLJ; PE §1.1 D1 |
| m ≥ 1.3Λ, b̂ < 5 (J ≤ 2 at 1.3Λ … J ≤ 49 at 20Λ) | positive | **[NE±]**, L ≥ 8.7 (∝ m); bulk-ν model tolerated to ×76 | PE §2.4 |
| m ≥ 1.3Λ, 5 ≤ b̂ < b̂_A(m) | ≥ 0.32 | [NE + M], L ≥ 3 (bulk-model dominated) | PE §2.4 |
| m ≥ 1.3Λ, b̂ ≥ b̂_A(m) | hybrid ≥ 0.93 | H_R [M]; hybrid floor [certified] on b̂ ∈ [91.2, 2Λ] | PE §2.4; ROUTEA |
| OO block J ≥ 12 | = `++` diagonal | [D given Γ_OO = Γ_++] (identity [P] at order 0, [NE] beyond) | PE D5 |

**No `m_*` exists above which all J are covered** (PE §2.4). Do not write "positivity proved for m ≥ m_*".

### 5.3 Result D at Δ_g = 300
| item | value | grade | source |
|---|---|---|---|
| binding J = 5 `+−`, 1.0165851Λ | own row `+0.20377` (published `+0.2038`); budget 5.8·10⁻², **`≥ +0.146`, L = 3.5**; no direct check | [NE±] | NG2 §0 |
| computed sets | even J ≤ 24, 28, 32; odd J ≤ 15; `(336, 28)` `+0.759` (0.751–0.772), `(336, 32)` `+0.837`; J = 0, 2 on a superseded small-ν treatment | [NE + M] | UT2 §5.1 |
| m ≥ 1.3Λ, b̂ < 5 | L ≥ 2.1 (≥ 4.4 at m ≥ 3Λ) | [NE + M] (not claimed [NE±]) | PE §2.4 |
| **H_LJ(300) Regge window, J ≥ 34** | — | **[I]** (why the headline is 1000) | HLJ §0 |

### 5.4 Result Q
| item | value | grade | source |
|---|---|---|---|
| binding pair `++`/OO J = 0 at the gap | FF/N `+0.1427 / +0.1428`; `s_d = +0.5096 / +0.5077` vs `s_fail = 1.776` | given H_RS ([D], one D-sketch sub-step) + truncation [NE], leverage ×85 / ×850 | HRS; PE Q1; PCYC |
| `s_d` m-independence | J = 0 (direct) and J = 2 (analytic-radial): constant to 0.3 % for m/Λ = 1–20 | [NE] | EES §5; HKR §4 |
| measured channels (`++` J = 0–8, OO J = 0, 2, `+−` J = 4–9, 11) | full-set minimum `+0.1426 ± 0.0017` | [NE + M given H_RS] | RQ-S §3.1 |
| unmeasured channels | `s_fail ≥ 11` (≥ 19× the largest measured `|s_d|`) | [NE + M] (H_fwd^Q) | RQ-S |
| 1000, 1.3Λ ≤ m ≤ 3Λ, J²/Δ ≤ 1 | L ≥ 3.95 against the certified 9/50 | **[NE±]** | PE §2.5 |
| 1000, m > 3Λ, J²/Δ ≤ 1 | L ≥ 1.65 (≥ 5.0 with the computed flat floor) | [NE + M] | PE §2.5 |
| J²/Δ > 1, m ≥ 1.3Λ | sign-free tolerance ≥ 2·10⁴ c₂ (1000) / 1.9·10³ (300) | [M] | PE §5.2 |
| 300, m ≥ 1.3Λ | not converted (L < 1 at m ≥ 3Λ, b̂ ≥ 5) | HYB [M] | PE §2.5 |

### 5.5 Light J = 0, 2 single traces (A1′)
| item | value | grade | source |
|---|---|---|---|
| `Ω^EE` flat | positive at every mass | [P] | A1P §4.2 |
| `Ω^EE`, leading-Regge model, Δ_L from 2.6 (J = 0) / 4.5 (J = 2) to Δ_g − 1 | `F/N ≥ 0.788 / 0.638`; curvature moves F/N by ≤ 1.8·10⁻⁴; the LR image is certified by Lemma C for Δ_L ≥ 4.5 (J = 0, 1000), ≥ 3.5 (J = 0, 300), every computed J = 2 | [P given the LR form] + [NE] | NG2 §0 |
| `Ω^EE` just below the gap, beyond LR | own rows at Δ_L = 990: L = 4362 (J = 0) / 913 (J = 2); Δ_L ≳ 0.7Δ_g: L ≥ 28 (1000) | 0.99Δ_g: [NE±]; ≳ 0.7Δ_g: [NE + M] (was [I]) | NG2 §0 |
| LR form = true AdS image at Δ_L = O(1) | — | [D-sketch] | NG2; A1P |
| `Ω^Q_up + εΩ^EE` on light J = 0 | ε = 0.040: flat margin 0.022N at m → 0; ε*: zero, sign [I] | ε = 0.040 [D-sketch] | NG2 |

### 5.6 Covering hypotheses (replaces UT2 §5.2)

| # | statement | for | grade | leverage | source |
|---|---|---|---|---|---|
| H_LJ(1000) | near-gap large spin | D | [D-sketch + M] **⟦HIGHSPIN⟧** | 2.7–4.8 | PE §5.3 |
| H_LJ(300) | same at 300; window J ≥ 34 | D | **[I]** | — | HLJ |
| H_R′ | beyond-Regge at m ≥ 1.3Λ, **only b̂ ≥ 5** (the strip b̂ < 5 is [NE±]) | D | [M] (b̂ < b̂_A: [NE + M]) | measured ≥ 15; scalar ≥ 150 | PE §5.3 |
| H_M | near-gap chain | D | J = 0, 2, 4, 5, 6: [NE±]; J = 8: [NE±] (λ-level); J ≥ 12: [NE + M] **⟦HIGHSPIN⟧**; other computed J: [M], validated [NE] | ≈ 30 → 78 at the binding state | NG2; PE |
| H_edge | edge window `x ∈ (0.8, Λ/m]` | D | J = 18, 19, 24, 25 [NE + M]; J = 29 [M] **⟦HIGHSPIN⟧**; rest [M] | ≥ 50 | NG2 |
| H_OO | OO = `++` beyond order 0 | D, Q | [D given identity] | ≥ 10³ | PE |
| H_cont | continuity in m; J = 15 now computed | D | [NE] | ≥ 10 | NG2 |
| H_bulk | bulk-ν (ν ≥ 20, J > 0) beyond-Regge envelope ≤ 3(450 + 100J) c₂ | D, Q | [M] | ×76 on S_D | PE §5.3 |
| H_mono | corrections grow at most like m at fixed (J, m/Λ) | corollary | J = 0: [D given RD] + [NE] (Γ = 21.0m (ν-flat) + g(ν) + O(1/m)); otherwise [M] | — | PE §2.6 |
| HYB (Q) | scalar skeleton + leading Regge elsewhere | Q | [M] minus [1.3, 3]Λ at 1000 | ≥ 300 (scalar) | PE |
| H_fwd^Q, H_smear^Q, H_LJ^Q | as UT2 | Q | [NE + M], [NE + M], [M] | 19; 8·10³; 1–10³ | PE |
| Pos_light | (Pos) on light J ≤ 2 states (A1′ only) | D, Q | §5.5 | — | NG2; A1P |
| σ_hU | double counting of the helicity-odd correction | — | **now the consistent reading** (EES §5), but still not used by the theorem | — | EES §5 |

**Net.** At Δ_g = 1000 every covering hypothesis is [M] or better; the weakest are H_LJ(1000) [D-sketch + M] and H_R′, HYB, H_bulk [M]
**⟦HIGHSPIN⟧**. At 300 one is [I]. As in the flat-space literature, (Pos) is the one hypothesis supported numerically rather than proved.

---

## 6. Derivations behind §4 (new since UT2; UT2 §6 otherwise unchanged)

* **EE step (EES §1).** On the light-cone circle `X₂₀, X₄₀` are affine in cos θ; the EE dressing is a degree-≤ 8 trigonometric polynomial
  over `(X₂₀X₄₀)⁴`; its sin-odd part integrates to zero; with `t = X₂₀/X₄₀` the even part is a Laurent polynomial `Σ_{b=−4}^{4} c_b t^b`
  with ν-independent `c_b` [P]; each `t^b` insertion is the scalar continuous-spin partial wave with `(a, b) = (0, ±b)` [D], checked
  end to end [DV]. The b-shifted blocks' edge factors have poles only at `Im ν = 3/2 + 2p`, where `c_l` vanishes [P], so LRC\* transfers.
* **f_pos(±3i/2) = 0 (EES §3; HKR §1.4).** `R₋₃^{EE3}(−3i/2) = −524288 ≠ 0`, so a pole of χ at `−3i/2` would contradict the
  strip-analyticity of `χR₋₃`; for the continuum weights there is no pole by construction [P]. Hence k = 4 criterion `J₀ < 5/2`.
  Scalar-type weights would have f_pos ≠ 0 (control: residue 9·10⁻⁷ vs ≤ 3·10⁻²⁷), so the generic `J₀ < 2` is real.
* **ε-DS out of H_K (HKR §3).** For explicit weights the corner count uses only the weight's analytic structure; ε enters the bounds only
  through values, `|I_W − I_Φ| ≤ ε_p P_{−½}(cosh β)` for all β (Mehler–Dirichlet).
* **L6′ (R4H §2).** In the og/CLS-D gauge `y = 2sin²(Δφ/2)` enters only through `X₁₃` and the 1–3 contractions; `f₂, f₄ ∝ y⁻¹`. Lemma Y
  [P]: every term of a spin-2 contact diagram is O(y⁰ log y) in the light-ray configuration. Hence `|L| ≤ CΔφ⁻²(1+|log Δφ|)`. In
  non-dispersive toy classes (half-integer free-field links) the collinear growth reaches the edge (slope 4.7 → 3.5), so dispersiveness is
  the crux; LEMMAS §4.1 item 4 and LRC §6.5 item 3 (η^{5/2}, p = 5/2) are superseded.
* **H_EFT(a) (R4H §3).** For dispersive terms only finitely many double-trace residues at bounded Mellin arguments enter, so the
  term-wise expansion converges pointwise in angle. F_NORM's "false pointwise in angle" is a property of the non-dispersive toy classes
  only.
* **c₁ (C1V).** See L3b, L3d, L3e.
* **Physical chain = D′ (PCYC §1).** CLS-D's functional on a block uses the Lorentzian shadow representation over the causal diamond
  `4 > 5 > 3`; with the gauge `x₃ = −x₁ = y`, `x₄ = −x₂ = x/(−x²)`, `x₅ = e`, the causal conditions are exactly
  `{1 − s + q > 0, u < 1}` for x and y: D′.
* **Light states (A1P §2–§4).** CHLPS's light-matter sum rules (validated against their printed `F_{g3}` exactly) applied to NF02 and Q3.

---

## 7. Gap dependence and the Δ_g ≥ 300 corollary

Unchanged from UT2 §7–§8 (Λ-forms [DV]; corollary [D-sketch + M], nothing below 300). Two updates: H_mono at J = 0 is now
[D given RD] + [NE] (PE §2.6); the near-gap evidence at 1000 now includes [NE±] states (§5.2), which does not change the corollary's
hypotheses (it is anchored at 300, where H_LJ(300)'s window is [I]).

---

## 8. ABJM: an in-scope consistency check (A1P §7)

* In the type-IIA 't Hooft regime (λ ≫ 1, `g_s → 0`) the higher-spin gap is the string scale: the lightest J ≥ 4 single trace is at the
  first massive level, `M² = 4/α′`, twist ≈ ML − 4, `ML = 2^{5/4}π^{1/2}λ^{1/4}` [D-sketch]. Λ ↔ ML is the spin-4 gap of CHLPS/CHMRS at
  leading order in 1/(ML).
* **`T × T` contains no light single traces at tree level** through the orders computed (supergravity + R⁴): only the stress-tensor
  multiplet is exchanged (CPY 1804.00949 §3 [Q]; BCP 1906.07195 [Q]); KK modes couple to two zero-mode gravitons at two derivatives only
  through `∫√g χ(y) = 0` (CHLPS §5 [Q]); the only stress-tensor-multiplet member in `T × T` is T itself (the other bosons are R-charged or
  have J = 1, which never enters `T × T`) [D-sketch]. Caveat: for 𝒩 = 6, `⟨SSSS⟩` fixes only 2 of 5 superconformal structures, so the
  `⟨TTTT⟩` statement rests on the bulk three-point argument [D-sketch]. **So ABJM-IIA satisfies A1′ and A2 at leading order in 1/λ**; the
  "light R-singlet in T × T" of UT2 §1/§9 and referee M12/M13 is not supported at leading order (the breathing-mode remark of A1P §7 is
  recalled, not re-derived).
* The flat limit of `⟨TTTT⟩` is the 10d IIA amplitude at zero internal momentum, so `g₄M⁶/8πG₄ = 2ζ(3)` [D-sketch + Q]; `ĝ₃ = 0`.
* Placement: 2ζ(3) at 3.34 % (300) / 2.59 % (1000) of the upper end (2a); `(Λ/ML)⁶` shifts these by ≈ −11 % / −3 % relative, so the
  placement is meaningful to one digit. Against CHMRS's 10d window `[0, 3.000]`, 80 %.
* **Status: an in-scope consistency check [D-sketch for "in scope"], not evidence for (Pos)** — (Pos) is not shown for string states in
  AdS, and at 3 % of the upper end the check could only fail grossly. M-theory regime not covered.

---

## 9. What is not claimed

* No unconditional bound; (Pos) is not established.
* No Δ_g other than 300, 1000 except through the corollary; nothing below 300.
* No loop-level statement; no bound on `ĝ₄` or `g₅`; no cubic closure from item (2); no ratio `|ĝ₃|²/g₄`.
* Not the δ = 0.05 functional.
* No explicit lower constant.
* No numerical evaluation of Lemma A on the Einstein correlator, or of the AdS₄ R⁴ `⟨TTTT⟩` (L6′ is derived, not evaluated).
* No value of `c₁` without F₂.
* Nothing about theories violating A1′ (light single traces of spin ≥ 4 below Δ_g), or with light spectra accumulating below Δ_g.
* Not "positivity proved for m ≥ m_*" (no such m_* exists for D); not "error-bounded" for states graded [NE + M] or [M].

---

## 10. HIGHSPIN placeholders (to be filled from `HIGHSPIN.md`)

| where | current content (kept) | to be replaced by |
|---|---|---|
| §5.2 rows "near gap J = 12, 16", "other computed near-gap J" | [NE + M] (NG2) | grade / leverage of near-gap J ≥ 12 at 1000 |
| §5.2 row "light-cone edge window", J = 29 | [M], unresolved | J = 29 edge value and grade |
| §5.2 rows "even J = 26–46", "near gap even J ≥ 50, odd J ≥ 33"; §5.6 H_LJ(1000) | [D-sketch + M], leverage 2.7–4.8 | H_LJ(1000) grade and leverage |
| §5.6 H_M, H_edge rows; net paragraph | as above | as above |

---

## 11. Inconsistencies: UT2 §11 resolved, and new ones

### 11.1 UT2's twenty items
| # | UT2 item | status now |
|---|---|---|
| 1 | box order | **resolved** (procedural): boxes are in time order; this document uses deliverable dates and states its reading order in the header. UT's own header still says "newest first" (not edited here). |
| 2 | DS for Q [NE] vs [DV] | **resolved, superseded**: DS is not an H_K input for the continuum weights (HKR §3); ε-DS is a value input, [D + certified] (C4), [DV] (C5). |
| 3 | k = 4 criterion (5/2 vs 2) | **resolved**: generic k = 4 weights `J₀ < 2` (LRC); Q's weights have f_pos(±3i/2) = 0, so **`J₀ < 5/2`** (EES §3; HKR §1.4: [P] for the continuum weights). UT2 and the old letter's "J₀ < 2" for Q are corrected. |
| 4 | "min 0.200 [certified]" | **resolved**: CERT §1 certifies `𝔉 ≥ N/5` everywhere; the statement is now true. |
| 5 | linear ε_{c₁} | resolved in UT2 (exact form); unchanged. |
| 6 | sign of (F) | resolved in UT2 (`𝒜₄[R⁴] = −2g₄`). |
| 7 | a = b integrability | **resolved**: L6′ [D] (R4H §2). EES §4's intermediate "subsumed by H_EFT(a); free-standing [I]" is superseded. |
| 8 | far-region number | resolved in UT2 (e^{−26} after volume factor); **now moot for H_RS**, which no longer goes through the thimbles (PCYC). The number applies only to crystal-cycle numerics. |
| 9 | RD fit "within 2 %" vs "3 %" | **still open in the record** (UT's RD box says 2 %, RD_CLOSE §0/§2.4 says 3 %); adopted 3 %. |
| 10 | Q binding at 1000 | resolved by adoption (+0.1428 binding; +0.1426 ± 0.0017 full-set minimum). |
| 11 | D's 300 minimum +0.2038 vs +0.2054 | **resolved**: NG2's own-mass row gives +0.20377 (published +0.20384); +0.2054 is HLJ's protocol. |
| 12 | most fragile large-J state at 300 | resolved by adoption ((336, 28) +0.759). |
| 13 | "17–30×" vs "≥ 19×" | resolved by adoption (≥ 19×). |
| 14 | `s_d` law | **resolved with a note**: REPRO's `0.5076 + 147/Δ²` describes the production chain with its finite reference (Δ_ref = 2·10⁴); with an infinitely heavy reference `s_d^∞ = 0.50957 − 128/D²` (ν < 12), m-independent to 0.3 % (EES §5). The Δ-drift of the production values is mostly the finite-reference artefact (−2.8·10⁻⁴ at D = 300 per EES; larger at 1000). No margin changes. |
| 15 | Λ⁴ → Δ_g⁴ | resolved by adoption (−1.99 %). |
| 16 | δk(0) at 100 | irrelevant (100 not claimed). |
| 17 | `X^Q` not published | **resolved**: `prl/data/k2/xq_abel_X_*.txt`, re-solved [DV] (CERT §3.1). |
| 18 | D's functional vs "published weights" | **resolved**: explicit continuum k = 2 weights for D and Q, ε-DS for all t [DV] (CERT §3.2). Either realisation gives the same R₀. |
| 19 | (P)/g₅ deciding test | **resolved**: run at J = 0 (direct) and J = 2 (analytic-radial); exact cancellation excluded; double counting holds; C5 piece not run [I] (EES §5; HKR §4). |
| 20 | CHMRS equation numbers | **still open**; the letter cites the section only. |

### 11.2 New inconsistencies found in this pass
| # | conflict | resolution |
|---|---|---|
| N1 | EES: H_K for Q [D given DS], DS [NE]; HKR: H_K for Q [D] | adopted HKR (later; uses the continuum weights' analytic structure). |
| N2 | R4H §5 (written in parallel): Q's upper end weakest link "H_K through ⟦EESTEP⟧" | superseded by EES/HKR (H_K [D]); upper end still [D-sketch] through g₅-(P), light-state links and κ. |
| N3 | A1P §6: Q's net implication "[I] in the a = b link" | superseded by R4H (L6′ [D]). |
| N4 | A1P: ε* = 0.038034, upper end 123.08 / 160.99; NG2: use ε = 0.040, 125.71 / 164.51 | adopted NG2; `c₁²` coefficient recomputed here: 1.294·10⁻⁴ / 1.033·10⁻⁶. |
| N5 | PE §3.1: binding state "not error-bounded", J = 5 treated as 2×2; NG2: [NE±], 1×1 block, tolerance 3.6·10⁴ c₂ (not 6.6·10³) | adopted NG2 (PE §3.1 slip). |
| N6 | A1P: light J ≤ 2 near the gap [I]; NG2: [NE + M] for Δ_L ≳ 0.7Δ_g, [NE±] at 0.99Δ_g | adopted NG2. |
| N7 | UT2/old letter: Q's flat margin "δ ≈ 0.2 (dense min 0.1966)" | CERT §2.1: true flat minimum 0.1835 (`b̂ ≈ 72`, J ≈ 3600); certified ≥ 9/50 on most regions; stored weights defective at `+−`, m = M, J ≥ 1990 ⇒ Q3\* definitional fix. |
| N8 | DESHIFT ε_W (C4) 1.2·10⁻¹³ / 4.2·10⁻¹¹ (midpoint sums) vs HKR certified 2.673·10⁻¹³ / 9.058·10⁻¹¹ | the certified numbers are upper bounds; quote them. |
| N9 | UT2 §2.1 cutoff (150, 6) for all Q pieces vs CERT §3.2 | Q's k = 2 piece uses (Λ + 100, 10). |
| N10 | RD_CLOSE census "no critical point at raw radius 0.1–39" | **withdrawn** (PCYC §4: 0–276 critical points at 0.3 ≤ r < 1.5); does not touch D′. |
| N11 | C1_COMPUTE §3.5 (independent seagull, Maldacena–Pimentel analogy) and §7.1 third bullet (shock-vertex factorisation) | **corrected** by C1V §8: no independent seagull; the factorisation that matters is through the energy channel (F₂). ADS4_LIGHT §0 item 3 / §2.3 (`t₄[V(ν) − V(−3i/2)]/(ν²+9/4)`) is identically zero in d = 3. |
| N12 | F_NORM B3/§3.4/§4 "H_EFT(a) false pointwise in angle" | restricted to non-dispersive toy classes (R4H §3). |
| N13 | PCYC grades H_RS "[D]" while stating one [D-sketch] sub-step | **not resolvable by adoption**: the weakest-link rule gives [D-sketch]; this document writes "[D], one [D-sketch] sub-step" and uses [D-sketch] in net grades. (H_RS is evidence, not an implication input, so no theorem grade changes.) |
| N14 | A1P §7 recalls the ABJM breathing mode at Δ = 6 without re-deriving it | open (does not affect the leading-order statement, which rests on CPY/BCP's exchange content). |
| N15 | UT2 / old letter: "ABJM violates A1 (light KK and R-singlet single traces)" | under A1′ and at tree level in `T × T`, withdrawn at leading order in 1/λ (A1P §7). Under the old A1 (all single traces) ABJM still fails. |

---

## 12. Source table (numbers new since UT2)

| number | value | source |
|---|---|---|
| NF02 certified floor | `𝔉 ≥ N/5`; sharp `0.2000000103` | CERT §1 |
| Q3\* certified floor / flat min | 9/50 / 0.1835 (J ≈ 3600, m ≈ 100M) | CERT §2 |
| C4 ε_W certified | 2.673·10⁻¹³ / 9.058·10⁻¹¹ | HKR §3.3 |
| k = 2 continuum ε_W | D 3.3·10⁻⁶ / 1.1·10⁻⁶; Q 1.9·10⁻⁷ / 7.4·10⁻⁸; ≤ 1.1·10⁻⁹ of ∫|ψ| | CERT §3.2 |
| `X^Q` re-solve | ≤ 8.3·10⁻⁷ of max ψ₂^Q (ν < 2) | CERT §3.1 |
| H_K criteria / margins | 3/2, 5/2 / ½, 3/2 | EES §0; HKR §0 |
| EE kernel end to end | 1.6·10⁻⁷–1.5·10⁻⁶ (9 b-shifted waves) | EES §1.4 |
| LRC-a EE k = 4 | ratio 1 ± 1.8·10⁻⁷; R_j to 1.1·10⁻¹¹ | HKR §1 |
| `R₋₃^{EE3}(−3i/2)` | −524288 | EES §3 |
| s_d^∞ (J = 0) | 0.5081 / 0.5095 / 0.5096 / 0.5095 at m/Λ = 1 / 3.35 / 10 / 20; fit 0.50957 − 128/D² | EES §5 |
| s_d^∞ (J = 2) | 0.5102 / 0.5085 / 0.5095 / 0.5096 | HKR §4 |
| L6′ | p = 2 vs 7/2; 997 836 structures | R4H §0 |
| H_EFT(a) check | 10⁻¹⁸–10⁻³⁷ | R4H §0 |
| c₁ tolerances | §3 table | C1V §6; this document |
| D binding 1000 | +0.20088; ≥ +0.198; L = 78 (130) | NG2 §0 |
| D binding 300 | +0.20377; ≥ +0.146; L = 3.5 | NG2 §0 |
| J = 4, 6, 8 (1000) | +0.20413, +0.22953, +0.42973; L ≥ 22, 5.0; 4.3–8.4 | NG2 §0 |
| J = 15 edge | +0.6900 | NG2 §0 |
| strip S_D | L ≥ 8.7; ×76 | PE §2.4 |
| Q large mass | L ≥ 3.95 (1.3–3Λ) | PE §2.5 |
| light flat actions | F_D⁰ ≥ 11.32; F_D² ≥ 0.448; F_Q² ≥ 1.97; F_Q⁰ < 0 for m < 0.4602Λ | A1P §4 |
| repair | ε = 0.040: 125.71 / 164.51; c₁² 1.294·10⁻⁴ / 1.033·10⁻⁶ | NG2; this document |
| light EE rows | F/N ≥ 0.788 / 0.638; L = 4362 / 913 at 0.99Δ_g | NG2 §0 |
| P_cyc region check | 0 mismatches / 2·10⁵ | PCYC §1.2 |
| scalar end to end | 2·10⁻⁸ (Δ = 20), 6·10⁻⁹ (30), 1·10⁻⁹ (50) | PCYC §0 |
| ABJM | 3.34 % / 2.59 %; (Λ/ML)⁶ shift −11 % / −3 %; λ ≈ 2.6·10⁷ / 3.2·10⁹ | A1P §7 |

All other numbers: UT2 §12.
