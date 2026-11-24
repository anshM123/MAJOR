# P4_WINDOW: P4 as a ν-window statement. Exact scalar test, the window/contour lemma, Theorem W and L′ with the remaining premises

*23 Nov 2026. Unreviewed, not committed. No existing file was edited.*
* *New code:* `ads4/code/p4win/`, reproduced by `run_all.sh <regge_moments_order8.m>` (≈ 0.7 CPU-h, at most 2 processes,
  `ulimit -v 3000000`, `OMP_NUM_THREADS=1`). Caches are in `code/p4win/cache/`, whose `.gitignore` ignores everything except itself.
  * `p4w_series.py`: CHMRS's scalar series to all its orders 1/m² … 1/m⁸, Legendre-reduced, with the flat-scaling part g_k of each order.
  * `p4w_common.py`: P_{−½+iν}(cosh β) and ∂_βP by Mehler–Dirichlet with the exact β-derivative (checked against mpmath to ≤ 3·10⁻¹²);
    series evaluation; NF02 ψ₂, ψ₃.
  * `p4w_exact.py`, `p4w_check.py`, `p4w_316.py`: the exact-scalar tests (§3).
  * `p4w_contour.py`: the constants of Lemma P4W (§2).
  * `p4w_states.py`: Theorem W's 47 sample states (§4).
  * `p4w_L.py`: Theorem L′'s range (§5).
* *Read-only inputs:*
  * `spinfunc/cache/sm_exact_1000.pkl`: the exact scalar (Δφ = 3) heavy moments on the full ν range at the gap states J = 40, 48, 64, 100
    (Arb Laplace/Barnes route of `exactskel/h13_exact`; SMOOTH_BOUND).
  * `p2con/cache/{run,ksh}_*.json`, `jstar_bh.json`, `p4b.json` and `pc_jstar`/`pc_run` (the computed O(1/m²) correction F₂, block + kernel O(r)).
  * `lsw/lsw_common` (published grid, ψ₂, floors), `largespin/ls_core` (NF02 weights).
  * `spinfunc/cache/lj_scan_1000.txt` (HLJ's travelling-wave decay rates, quoted).
  * LRINF_TIGHT §5 (F_cert, F_num, quoted).
  * The CHMRS ancillary `regge_moments_order8.m` (scratch copy).

Grades as in UNIFIED_THEOREM §0. Conventions as in P2_CONTRACT, PW_IMAGES and LARGESPIN_W:
* f is the de-shifted image per unit leading image (f₀ = P := P_{−½+iν}(cosh β));
* the k = 2 cost of a correction is `2∫₀¹ψ₂ f dp + 2D·F(i/2)`, with p = ν/Λ and D = −9.98·10⁻⁴ (PW4);
* x = ν/m, y = x², b̂ = (2J+1)/(m/Λ), the window edge is ν_w = c/β (c = 3 unless stated), and p_w = ν_w/Λ;
* **G_e := exact − hybrid** for the entry e of M_Ω (the full beyond-hybrid correction). The hybrid is flat + leading-Regge curvature, so
  per unit leading image `f_hyb = flat_e + (P − J₀(b̂p))`;
* `F₂^bh` is P2_CONTRACT's computed O(1/m²) correction minus the hybrid's own flat O(y) tower;
* **R := G − F₂^bh/m²** is the P4 remainder: beyond-hybrid, beyond O(1/m²), all orders.

---

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

## 0. Summary

### 0.1 Table

| # | claim | grade | constants | coverage |
|---|---|---|---|---|
| **E1** | **The full-range failure of P4 is an artefact of the truncated series.** On the exact scalar, the true P4 remainder R_ex costs **≤ 1.5·10⁻³ of the budget** (K₄-equivalent 0.26–8.1). The truncated fixed-β series costs up to 0.72 of the budget (K₄ up to 8.5·10³), and gets *worse* with more terms (asymptotic). | [DV] (exact Arb data, 4 states) | J = 64: K₄(exact) 0.99 against the series' 16.1 / 75 / 441 (truncated at 1/m⁴ / 1/m⁶ / 1/m⁸); J = 100: 8.1 against 176 / 991 / 8.5·10³ | scalar Δφ = 3, gap, J = 40, 48, 64, 100 |
| E2 | **In the Regge window ν < 3/β the fixed-β series converges and agrees with the exact image.** The window cost of G equals the series' to ≤ 0.3 %. Successive orders fall by q ≤ 6.5·10⁻³ on every one of the 47 states. | [DV] | ratio q ≤ 6.5·10⁻³ (J = 26) … 5.5·10⁻⁵ (J = 250) | scalar; all 47 states (series), 4 states (exact) |
| E3 | **The (3/16)J₀y² residual of P2_CONTRACT is real and is beyond-hybrid.** The scalar hybrid E₂P_J has flat O(y²) term g₂ + (3/16)J₀ exactly; CHMRS's (and the exact moment's) is g₂. It is now inside the P4 window term, and harmless there. | [DV] (to 5 digits at 4 b) + [NE] (exact fit −0.12…−0.27 against −0.1875) | — | — |
| **L1** | **Lemma P4W (window/contour).** \|cost(G_e)\| ≤ \|window cost\| + 4·[V(p_w)M_w + V(1)M₁ + Top] + 2\|D\|\|Q\|, given (AG′) on [p_w, 1]. V and Top are explicit integrals of \|ψ_k(p + iy)\|e^{−κy}. **No 1/m series is used on ν ≥ ν_w; no analyticity is used on ν < ν_w** (a real-axis sup bound suffices there). | **[P]** (Cauchy) given (AG′); constants [DV] | §2 | every state, any J, m |
| E4 | Lemma P4W holds on the exact scalar, with the scalar-calibrated envelope e_sc (\|G\| ≤ 0.65 e_sc on every node). The bound (k = 2 + 3) is 4.0–7.9 % of the budget, against the true \|cost\| ≲ 10⁻³; it overestimates the complement by ×6–8 (k = 2) and ×8–17 (k = 3). | [DV] | bound / budget: 0.040, 0.045, 0.056, 0.079 (J = 40, 48, 64, 100) | 4 exact states |
| **W4** | **P4 at Theorem W's 47 states, window form.** With (AG′) and the scalar-calibrated edge envelope (Λ_AG = 1): **d/H ≤ 0.49 at every sampled state**, with window + boundary ≤ 0.49 (computed) and complement ≤ 9·10⁻³. Tolerated: edge multiple Λ_AG up to 114 (plain envelope); spinning O(1/m⁴) window coefficient K₄w up to 20.8; window sup-norm ε_w up to ≥ 2·10³ × sup\|S₄/m⁴\|. | [P given P0, P1′, W, AG′, EDGE] + [DV] | per state, §4 | all 47 states (J = 26–344, m ≤ 1.3Λ): **J\* = 26** with the plain envelope |
| W4′ | **The same with the Ω^U +− helicity transfer of the edge envelope (×(1−x₁²)^{−2}).** Λ_AG = 1 fails only at the four states m ≈ 1.03Λ, J = 26–38; at Λ_AG = 3 it fails at m ≤ 1.05Λ, J ≤ 50. Equivalently, the edge premise is \|ρ_{+−}(x₁)\| ≤ ρ_tol, with ρ_tol ≥ 0.15 everywhere and ≥ 1 for J ≥ 46 (all m), against the exact scalar's 0.27–0.69. | [P given …] + [DV]; transfer [M] | J\* = **42** (Λ_AG = 1), 51 (3), 100 (10), 150 (30); with \|ρ\| ≤ 0.7: 42; with \|ρ\| ≤ 1: 46 | §4.2 |
| **L4** | **P4 on Theorem L′'s range, window form.** With the O(1/m²) window term at Λ_amp = 3 (scalar-shaped, L′'s P_L′(i)), (AG′) and the transferred edge envelope, P4 plus the window O(1/m²) use **≤ 8.6 % of the certified half-floor** (28 % at Λ_amp = 10). The edge tolerance is ≥ 600× the scalar envelope. | [P given P_L′(i), AG′, EDGE] + [DV] | ε_R = 0.086 (J = 223) → 3.8·10⁻⁴ (J = 10⁶) | 16 states, J = 223 … 10⁶, m = m_gap … 5Λ (sampled) |
| AG′ | **What is assumed beyond all orders.** exact − flat is, on [ν_w, Λ] × [0, 25/(0.97β)], a pair of travelling waves, analytic, with the local decay rate 0.97β/√(1−x²) for y ≤ 15/κ and 0.97β beyond. The amplitude enters only at the edge p = 1 and on the top side (EDGE). | structure: [D] to all orders of the saddle expansion (HLJ §3); local rate: [NE] (lj_scan); beyond all orders: **[A]** (H_AG, now quantitative) | rate threshold T1 = 15/κ; T1 = 10/κ lowers Λ_AG,tol to 18.7 (plain) / 0.046 (transferred) | §6 |
| S4 | The spinning O(1/m⁴) symbols and the kernel's O(r²) term are **not needed** in computed form. They enter only through the window (tolerance K₄w ≥ 20.8, against the model's ≤ 3.8) and the O(1/m⁴) monopole (tolerance \|Q₄\| ≥ 9.7·10⁸, in the units where the O(1/m²) block monopole is 12.8). **Not computed.** | [DV] (tolerances); values [M] | — | — |

### 0.2 Bottom line
* **P4 is reformulated and reduced, not proved.**
  * The ν-window version replaces "the O(1/m⁴) series costed to ν ≈ m is ≤ K₄" (false: §1) by three statements:
    * **(W)** a real-axis sup bound in the Regge window ν < 3/β, where the series converges (ratio ≤ 6.5·10⁻³) and matches the exact
      scalar;
    * **(AG′)** the travelling-wave analytic structure on ν ≥ 3/β, i.e. the beyond-all-orders hypothesis H_AG, now with explicit
      constants;
    * **(EDGE)** an amplitude bound at the single point p = 1 (ν = Λ).
  * The implication "W ∧ AG′ ∧ EDGE ⇒ P4 with explicit constants" is **[P]** (Lemma P4W). The constants are [DV] at every sampled state.
* **What is closed.**
  * At all 47 Theorem W states the computed part (window O(1/m²) + O(1/m⁴…⁸) series + boundary term) uses ≤ 0.49 of H.
  * With the scalar-sized complement, everything uses ≤ 0.49 of H. So **Theorem W holds on all sampled J ≥ 26 given P0, P1′, W, AG′, EDGE**,
    with explicit tolerances.
  * At all 16 sampled Theorem L′ states P4 + window O(1/m²) use ≤ 8.6 % of the certified floor.
* **Where it does not close.**
  * The edge premise (EDGE) for the **+− helicity at m ≲ 1.05Λ and J ≤ 50**.
  * If the graviton's +− beyond-hybrid correction at ν = Λ is as large, relative to its own flat envelope, as the exact scalar's
    (\|ρ\| ≈ 0.3–0.7), then (26, 1.026Λ), (30, 1.030Λ), (34, 1.033Λ) and (38, 1.037Λ) fail the sign-free bound. They are
    near-twist-edge states, 3 % above the gap mass.
  * There the floor is also [NE] (LRINF_TIGHT §4: no certified floor for m < 1.3Λ below b̂ = 245).
* **What remains, precisely (§7):**
  * (AG′) beyond all orders [A];
  * the graviton edge amplitude at x₁ = Λ/m [M];
  * the spinning O(1/m⁴) window coefficients [M], tolerance ≥ 20.8;
  * the boundary-amplitude estimate (the factor 2 for the travelling-wave split) [D-sketch];
  * between sampled states [NE].

---

## 1. Why the full-range P4 fails, and why that does not matter (exact scalar, `p4w_exact.py`)

**Object.** At the four exact gap states (Δ_g = 1000, Δφ = 3, J = 40, 48, 64, 100):
* R_ex := f_ex − f_hyb − [T₁ − ν²g₁(b̂p)]/m² is the exact P4 remainder of the scalar analogue: everything beyond the hybrid and beyond
  CHMRS's O(1/m²).
* R_ser^{(K)} := Σ_{k=2}^{K}[T_k − ν^{2k}g_k(b̂p)]/m^{2k} is its truncated fixed-β series. K = 2 is P2_CONTRACT §2.2's R₄. My code
  reproduces pc_p4b's 16.1 (J = 64) and 176 (J = 100) exactly.

**Full-range costs** (budget = H − d₂ for the W sample states J = 64, 100; recomputed H at J = 40, 48):

| J | cost(R_ex)/budget | K₄-equiv (exact) | series K = 2 (1/m⁴) | K = 3 (1/m⁶) | K = 4 (1/m⁸) |
|---|---|---|---|---|---|
| 40 | +1.5·10⁻³ | 0.40 | 0.41 | 7.3 | 18 |
| 48 | −3.1·10⁻⁴ | 0.26 | 5.8 | 11 | 56 |
| 64 | +4.8·10⁻⁴ | **0.99** | **16.1** | 75 | 441 |
| 100 | −6.8·10⁻⁴ | **8.1** | **176** | 991 | 8.5·10³ |

(K₄-equivalent = \|cost\|·m⁴/cost(η²ν⁴P), P2_CONTRACT's unit; the K = 3, 4 columns are K₄-equivalents of the truncated sums.)

**Pointwise.**
* \|f_ex − (P + Σ_{k≤K}T_k/m^{2k})\| improves with K up to x ≈ 0.3–0.4: at J = 64, ν = 271 it is 1.1·10⁻², 6.8·10⁻⁴, 1.3·10⁻⁴,
  7.8·10⁻⁵ for K = 1…4.
* Beyond x ≈ 0.4 it gets worse with K (x = 0.8: 0.87, 3.3, 9.3, 51). This is the asymptotic, non-convergent fixed-β series.
* The exact beyond-hybrid correction stays bounded: \|G_ex\| ≤ 4(νβ)²/m² for x ≤ 0.7 and ≤ 7.3(νβ)²/m² for x ≤ 0.8. \|R_ex\| reaches ≈ 30 only near
  x → 1, where ψ₂ ∝ (1−p)⁵.
* It costs ≤ 1.5·10⁻³ of the budget because it oscillates at frequency ≈ β against the smooth ψ₂. This is the AG mechanism.

**Reading [DV].** P2_CONTRACT's "K₄ = 20 … 7·10³ at J ≥ 80" measures the series, not the remainder. At J = 100 the exact remainder is 22×
below the K = 2 series and inside P2_CONTRACT's tolerance (K₄ tol = 44.7). A full-range bound must not use the fixed-β series beyond
x ≈ 0.2–0.3. That is the whole content of the window reformulation.

## 2. Lemma P4W

### 2.1 Statement
Fix a state (J, m), with β = 2 asinh((J+½)/m), and an entry e. Let G := G_e. Choose the window edge p_w = ν_w/Λ, ν_w = c/β. Assume:
* **(W) window, real axis only.** For ν < ν_w, G = F₂^bh/m² + S₄/m⁴ + R_w, with \|R_w(ν)\| ≤ ε_w.
  * F₂^bh is computed (P2_CONTRACT).
  * S₄ is the O(1/m⁴) beyond-hybrid term: the scalar seed is computed; the spinning part is K₄w × (window costs of η²ν⁴P,
    η sinh β ν⁴∂_βP).
* **(AG′) complement, analytic.** On the rectangle Q_± = {p_w ≤ Re p ≤ 1, 0 ≤ ±Im p ≤ T}, T = 25/(0.97Λβ):
  * G = G₊ + G₋ with G_± analytic on Q_± and continuous on its closure;
  * G_± = g_± − (P − J₀(b̂p))_±, where g := exact − flat;
  * **g_± obeys \|g_±(p ± iy)\| ≤ M(Re p)·e^{−κ(Re p) y}**, with κ(p) = 0.97Λβ/√(1−x²) for y ≤ 15/κ(p), then the rate 0.97Λβ beyond (a
    continuous majorant);
  * (P − J₀)_± are the known Harish-Chandra / Hankel waves of the hybrid's curvature (uniform rates Λβ, b̂ ≥ Λβ).
* **(EDGE)** M(1) ≤ M₁: the amplitude of g at p = 1 (x₁ = Λ/m). On the top side, M(p) ≤ M_top.
* **Boundary.** At p_w the waves of G are those of its window series (ring classes, rate β; flat waves, rate b̂), with amplitude
  M_w := 2·sup_{[0.7, 1.3]ν_w}\|F₂^bh/m² + S₄/m⁴\|.

**Then**

    |cost₂(G)| ≤ 2|∫₀^{p_w}ψ₂(F₂^bh/m² + S₄/m⁴)| + 2ε_w∫₀^{p_w}|ψ₂|
               + 4[ V(p_w; κ_w) M_w + V(1; κ₁)(M₁ + A) + V(1; 0.97Λβ) A + Top ]  + 2|D||G_F(i/2)|,          (2.1)

with:
* V(p₀; κ) = ∫₀^T\|ψ₂(p₀ + iy)\|e^{−κy}dy;
* A = the amplitude at p = 1 of (P − J₀)_±. Numerically, A = 1.2[max_{p∈[0.9,1]}\|P − J₀\| + (β²/2)max\|P\|]: the two waves nearly cancel
  at small β;
* Top = e^{−0.97ΛβT}∫\|ψ₂(p + iT)\|M_top dp.

The k = 3 piece is the same, with weight u and ψ₃ (no de-shift).

**Proof [P].**
1. Linearity of the cost (PW4: the cost of 𝒜_g[F] is ∫ψ₂F/c + D·F(i/2), exactly) splits the real-axis integral at p_w.
2. On [0, p_w], (W) is used pointwise.
3. On [p_w, 1], each ∫ψ₂G_± is moved by Cauchy's theorem to the three other sides of Q_±.
4. On the right side, \|G_±\| ≤ \|g_±\| + \|(P − J₀)_±\|.
5. Each side is bounded by its majorant.
6. The monopole term is separate and exact (PW4). ∎

**Remarks.**
* No 1/m series is used on [p_w, 1]. No analyticity is used on [0, p_w).
* The fixed-β series is used only where it converges (E2), and only for the window integral and the boundary amplitude M_w.
* ψ₂ = (1−p)⁵pX(p)/k₀ vanishes to fifth order at p = 1, so V(1; κ) ≈ \|c₅\|·5!/κ⁶ with \|c₅\| = 8.48·10³. The edge term is suppressed by
  (Λβ/√(1−x₁²))⁻⁶. This is why the amplitude premise is needed only at the single point ν = Λ.
* The boundary term at p_w is the price of a sharp window. It is computed (it involves only the series), and it is the largest
  computed piece (§4).
* c = 2 lowers the computed part to ≤ 0.34 H; c = 5 raises it to 1.05 H at J = 344. So **c = 3 is used**.

### 2.2 Inputs and their grades
| input | value | grade |
|---|---|---|
| CHMRS scalar series, orders 1…4 | `cache/series.txt`; g₁, g₂ reproduce P2_CONTRACT; g₃, g₄ new | [Q] (ancillary) + [P] (Legendre reduction) |
| series condition | \|r₀P\| + \|r₁P′\| ≤ 15·\|T₄\| on the grid (no catastrophic cancellation) | [DV] |
| (P, ∂_βP) | Mehler–Dirichlet with exact derivative; against mpmath ≤ 3·10⁻¹² abs (ν ≤ 900) | [DV] |
| F₂^bh (graviton O(1/m²)) | p2con caches, block + kernel O(r), minus pc_jstar.flat_c1 | [DV] (P2_CONTRACT) |
| D, cost identity | PW4 / INDEP_CHECK2 | [DV] / [D] |
| local rate | lj_scan_1000 (HLJ): mean rate/β along vertical lines at ν_r = 898.6, y ≈ 2/β: 1.97, 1.91, 1.80, 1.69, 1.62, 1.51, 1.48 (J = 50, 64, 100, 150, 200, 300, 344), against 1/√(1−x²) = 1.95, 1.88, 1.75, 1.63, 1.54, 1.42, 1.39; at ν_r = 499: 1.14–1.09 against 1.14–1.08 | [NE] |

## 3. Lemma P4W on the exact scalar (`p4w_check.py`, `p4w_316.py`)

### 3.1 The envelope (EDGE calibration) [DV]

    e_sc(p) = [10 + (νβ)²(4 + 5(1−x²)^{−3})]/m²     (ν = Λp, x = ν/m)

* It bounds \|G_ex\| at every node of all four exact states: **max \|G\|/e_sc = 0.65** (k = 2), 0.43–0.52 (k = 3).
* The bulk (x ≤ 0.8) is \|G\| ≲ 4(νβ)²/m², i.e. the O(y) AdS/flat kinematic mismatch. The edge grows like (1−x²)^{−3}.
* Relative to its own envelope E₂/2, the exact scalar's edge correction is **\|ρ₂\| = 0.27 (J = 40), 0.34 (48), 0.38 (64), 0.69 (100)** at
  x ≥ 0.95 (0.9 for J = 100).

### 3.2 Window and complement [DV]

c = 3, k = 2 in units of the budget (k = 3 in brackets):

| J | window: exact G | window: series K = 1 / 2 / 4 | complement, actual | contour bound [bd, edge, top] (local rate) | bound/budget, total (k = 2 + 3) |
|---|---|---|---|---|---|
| 40 | +4.39·10⁻³ | +4.38 / 4.39 / 4.39·10⁻³ | −3.95·10⁻³ (−3.2·10⁻⁴) | [2.6·10⁻², 3.5·10⁻³, 6·10⁻⁵] (k=3 [3.7·10⁻³, 1.6·10⁻³]) | 0.040 |
| 48 | +5.11·10⁻³ | +5.10 / 5.11 / 5.11·10⁻³ | −4.81·10⁻³ (−5.3·10⁻⁴) | [3.2·10⁻², 2.9·10⁻³, 3·10⁻⁵] | 0.045 |
| 64 | +6.74·10⁻³ | +6.74 / 6.74 / 6.74·10⁻³ | −7.02·10⁻³ (−7.9·10⁻⁴) | [4.1·10⁻², 2.1·10⁻³, 1·10⁻⁵] | 0.056 |
| 100 | +9.67·10⁻³ | +9.66 / 9.66 / 9.66·10⁻³ | −9.63·10⁻³ (−9.3·10⁻⁴) | [6.0·10⁻², 1.3·10⁻³, 8·10⁻⁶] | 0.079 |

**Reading.**
* The window cost of the exact G is reproduced by its series to ≤ 0.3 %, already at K = 1 or 2.
* The complement is ≈ −(window): the sharp cut produces equal and opposite boundary pieces, and their sum (the true cost) is ~10⁻³.
* Lemma P4W's bound holds at every state and is conservative by ×6–8 on the k = 2 complement and ×8–17 on the k = 3 one.
* With the uniform rate (no local-rate premise) the edge term would be 10 / 4.7 / 1.5 / 0.26 of the budget. **The local-rate part of
  (AG′) is load-bearing near the edge for m → Λ.**

### 3.3 The (3/16)J₀y² residual (P2_CONTRACT residual 2) [DV + NE]
* **Hybrid side.** The y-expansion at fixed b of the scalar hybrid's flat part E₂(y)P_J(1−2y)/2 has O(y²) coefficient **g₂ + 0.18750 J₀(b)**,
  to 5 digits at b = 0.8, 1.7, 2.7, 3.6 (mpmath, 40 digits).
* **Exact side.** A fit of (R_ex − R_ser^{(4)})m⁴/ν⁴ = aJ₀(b) + c₀ on the window nodes gives a = −0.26, −0.27, −0.24, −0.14 (ν ≤ 60) and
  −0.25, −0.25, −0.22, −0.12 (ν ≤ 150) for J = 40, 48, 64, 100. This is consistent with −3/16 = −0.1875 within the data noise, which is
  ≈ 10⁻⁹ absolute at ν ≲ 20.
* **Conclusion.** CHMRS's flat part g₂ is the exact moment's. The hybrid's flat partial wave carries an extra +(3/16)J₀y², so **the exact
  moment differs from the hybrid by −(3/16)ν⁴J₀/m⁴ at O(1/m⁴).** It is a genuine beyond-hybrid term, not a normalisation slip.
  * Its window cost is included in §4 (S₄ includes it).
  * Its full-range K₄-equivalent ≤ 0.3 (P2_CONTRACT) is harmless.

## 4. Theorem W's 47 sample states (`p4w_states.py`; `cache/states15.txt` is the primary run, `states10`/`states40` the rate-threshold sensitivity)

### 4.1 Per state

**Columns.** c = 3; θ = 1; the local rate up to y = 15/κ; budget H from jstar_bh.
* "win" = max over entries of the window cost of F₂^bh (k = 2, plus k = 3 sign-free with u\|ψ₃\|), plus the scalar O(1/m⁴, ⁶, ⁸) window
  costs (incl. the 3/16 term).
* "bd" = the boundary term (k = 2 + 3).
* "AG" = edge + top with the plain scalar envelope at Λ_AG = 1, including the known curvature-wave term A.
* "d/H" = total at Λ_AG = 1: plain / with the +− helicity transfer.
* "Λ_AG tol" = the tolerated multiple of e_sc at the edge: plain / transferred.
* "ρ_tol" = the tolerated \|g_{+−}(x₁)\| relative to the +− hybrid envelope E₂(x₁)(1−x₁²)^{−2}/2.
* "K₄w tol" = the tolerated spinning O(1/m⁴) coefficient in the window (4a + 4b units, including its effect on M_w).
* "ε_w/S₄" = the window sup-norm tolerance over sup\|S₄/m⁴\|.
* "P2C K₄" = P2_CONTRACT's full-range series K₄ against its tolerance (for comparison).

| J | m/Λ | win/H | bd/H | AG/H | d/H plain / tr | Λ_AG tol plain / tr | ρ_tol(+−) | K₄w tol | ε_w/S₄ | P2C K₄ (tol) |
|---|---|---|---|---|---|---|---|---|---|---|
| 26 | 1.026 | 7.6e-3 | 0.06 | 9.0e-3 | 0.08 / **3.30** | 114 / **0.29** | **0.15** | 20.8 | 2.0e3 | 0.04 (21) |
| 26 | 1.128 | 5.4e-3 | 0.05 | 8.3e-3 | 0.06 / 0.21 | 130 / 6.0 | 0.40 | 24 | 2.3e3 | 1.07 (1.1) |
| 26 | 1.231 | 4.0e-3 | 0.03 | 8.6e-3 | 0.05 / 0.10 | 128 / 15 | 0.41 | 28 | 2.7e3 | 0.16 (1.2) |
| 30 | 1.030 | 8.5e-3 | 0.07 | 7.3e-3 | 0.09 / **2.09** | 140 / **0.46** | 0.26 | 31 | 2.9e3 | 0.52 (1.7) |
| 30 | 1.133 | 6.2e-3 | 0.05 | 6.9e-3 | 0.07 / 0.19 | 155 / 7.5 | 0.64 | 36 | 3.4e3 | 0.58 (1.9) |
| 30 | 1.236 | 4.6e-3 | 0.04 | 7.0e-3 | 0.05 / 0.10 | 156 / 19 | 0.66 | 41 | 3.9e3 | 0.17 (2.1) |
| 34 | 1.033 | 9.6e-3 | 0.08 | 6.3e-3 | 0.10 / **1.47** | 160 / **0.66** | 0.41 | 44 | 4.0e3 | 2.0 (2.3) |
| 34 | 1.137 | 7.0e-3 | 0.06 | 5.8e-3 | 0.07 / 0.17 | 181 / 9.3 | 0.97 | 51 | 4.8e3 | 0.54 (3.4) |
| 34 | 1.240 | 5.2e-3 | 0.04 | 5.9e-3 | 0.06 / 0.09 | 185 / 23 | 1.0 | 58 | 5.5e3 | 0.17 (3.9) |
| 38 | 1.037 | 1.1e-2 | 0.09 | 5.5e-3 | 0.11 / **1.08** | 183 / **0.92** | 0.63 | 60 | 5.4e3 | 0.33 (21) |
| 38 | 1.141 | 7.7e-3 | 0.07 | 5.1e-3 | 0.08 / 0.16 | 207 / 11 | 1.4 | 70 | 6.4e3 | 0.49 (7.5) |
| 38 | 1.245 | 5.8e-3 | 0.05 | 5.1e-3 | 0.06 / 0.09 | 214 / 27 | 1.4 | 80 | 7.4e3 | 0.13 (8.3) |
| 42 | 1.041 | 1.2e-2 | 0.10 | 4.9e-3 | 0.12 / 0.83 | 203 / 1.23 | 0.91 | 79 | 6.9e3 | 1.9 (4.8) |
| 42 | 1.145 | 8.4e-3 | 0.07 | 4.5e-3 | 0.09 / 0.15 | 231 / 13 | 1.9 | 92 | 8.4e3 | 0.35 (28) |
| 42 | 1.249 | 6.3e-3 | 0.05 | 4.4e-3 | 0.07 / 0.09 | 243 / 31 | 2.0 | 105 | 9.7e3 | 0.03 (34) |
| 46 | 1.045 | 1.3e-2 | 0.11 | 4.4e-3 | 0.13 / 0.67 | 224 / 1.6 | 1.3 | 102 | 8.8e3 | **26.6 (4.8)** |
| 46 | 1.150 | 9.1e-3 | 0.08 | 4.0e-3 | 0.09 / 0.15 | 256 / 15 | 2.5 | 119 | 1.1e4 | 0.01 (39) |
| 46 | 1.254 | 6.8e-3 | 0.06 | 3.9e-3 | 0.07 / 0.09 | 273 / 36 | 2.7 | 137 | 1.2e4 | 0.42 (30) |
| 50 | 1.049 | 1.3e-2 | 0.12 | 4.0e-3 | 0.13 / 0.56 | 243 / 2.0 | 1.7 | 128 | 1.1e4 | 2.9 (13) |
| 50 | 1.154 | 9.8e-3 | 0.09 | 3.7e-3 | 0.10 / 0.15 | 281 / 17 | 3.3 | 151 | 1.3e4 | 0.74 (19) |
| 50 | 1.259 | 7.3e-3 | 0.06 | 3.6e-3 | 0.08 / 0.09 | 302 / 41 | 3.4 | 173 | 1.5e4 | 1.4 (15) |
| 64 | 1.062 | 1.7e-2 | 0.15 | 3.0e-3 | 0.17 / 0.37 | 311 / 4.0 | 4.2 | 256 | 2.0e4 | 16.1 (16.7) |
| 64 | 1.168 | 1.2e-2 | 0.11 | 2.8e-3 | 0.12 / 0.15 | 366 / 26 | 7.0 | 304 | 2.4e4 | 10.5 (16.7) |
| 64 | 1.275 | 9.0e-3 | 0.08 | 2.6e-3 | 0.09 / 0.11 | 402 / 60 | 7.4 | 350 | 2.9e4 | 14.0 (15.5) |
| 80 | 1.077 | 2.0e-2 | 0.18 | 2.4e-3 | 0.20 / 0.31 | 385 / 7.4 | 9.1 | 476 | 3.2e4 | **163 (22)** |
| 80 | 1.185 | 1.4e-2 | 0.13 | 2.1e-3 | 0.15 / 0.17 | 460 / 38 | 14 | 572 | 4.1e4 | **35 (29)** |
| 80 | 1.293 | 1.1e-2 | 0.10 | 2.0e-3 | 0.11 / 0.12 | 515 / 83 | 14 | 664 | 5.0e4 | 20 (37) |
| 100 | 1.096 | 2.3e-2 | 0.21 | 1.9e-3 | 0.24 / 0.29 | 471 / 13 | 19 | 882 | 5.1e4 | **176 (45)** |
| 100 | 1.205 | 1.7e-2 | 0.15 | 1.7e-3 | 0.17 / 0.19 | 574 / 56 | 27 | 1.1e3 | 6.8e4 | 15 (245) |
| 120 | 1.114 | 2.6e-2 | 0.24 | 1.5e-3 | 0.27 / 0.30 | 552 / 21 | 34 | 1.5e3 | 7.4e4 | **300 (78)** |
| 120 | 1.225 | 1.9e-2 | 0.18 | 1.4e-3 | 0.20 / 0.21 | 683 / 76 | 46 | 1.8e3 | 1.0e5 | **5.5e3 (71)** |
| 150 | 1.140 | 3.1e-2 | 0.28 | 1.2e-3 | 0.31 / 0.33 | 659 / 36 | 69 | 2.7e3 | 1.1e5 | **634 (136)** |
| 150 | 1.254 | 2.2e-2 | 0.21 | 1.1e-3 | 0.23 / 0.24 | 837 / 112 | 86 | 3.4e3 | 1.6e5 | **297 (288)** |
| 200 | 1.183 | 3.6e-2 | 0.34 | 9.4e-4 | 0.38 / 0.38 | 783 / 69 | 160 | 5.7e3 | 1.9e5 | **2.8e3 (337)** |
| 250 | 1.225 | 4.1e-2 | 0.38 | 8.1e-4 | 0.42 / 0.43 | 813 / 108 | 290 | 1.0e4 | 3.0e5 | **4.3e3 (847)** |
| 300 | 1.265 | 4.4e-2 | 0.42 | 7.8e-4 | 0.46 / 0.47 | 748 / 147 | 470 | 1.6e4 | 2.7e5 | **7.2e3 (749)** |
| 344 | 1.299 | 4.7e-2 | **0.44** | 8.1e-4 | **0.49** / 0.49 | 643 / 171 | 650 | 2.3e4 | 2.4e5 | **1.5e3 (1.1e3)** |
| 33 | 1.033 | 1.8e-3 | 0.04 | 3.2e-3 | 0.04 / 0.77 | 337 / 1.31 | 0.80 | 88 | 8.1e3 | 0.03 (17) |
| 33 | 1.136 | 1.3e-3 | 0.03 | 3.0e-3 | 0.03 / 0.08 | 362 / 18 | 1.8 | 98 | 9.2e3 | 0.30 (5.9) |
| 33 | 1.239 | 9.1e-4 | 0.02 | 3.0e-3 | 0.03 / 0.04 | 370 / 45 | 1.9 | 112 | 1.1e4 | 0.38 (6.4) |
| 51 | 1.050 | 2.5e-3 | 0.06 | 1.9e-3 | 0.06 / 0.26 | 539 / 4.7 | 4.0 | 294 | 2.5e4 | **17 (13)** |
| 51 | 1.155 | 1.7e-3 | 0.04 | 1.8e-3 | 0.05 / 0.07 | 610 / 38 | 7.4 | 339 | 2.9e4 | 9.9 (15) |
| 51 | 1.260 | 1.2e-3 | 0.03 | 1.7e-3 | 0.04 / 0.04 | 647 / 89 | 7.7 | 383 | 3.4e4 | 2.4 (19) |
| 101 | 1.096 | 3.9e-3 | 0.11 | 9.2e-4 | 0.11 / 0.14 | 1.1e3 / 32 | 46 | 2.1e3 | 1.2e5 | 56 (167) |
| 101 | 1.206 | 2.7e-3 | 0.08 | 8.3e-4 | 0.08 / 0.09 | 1.3e3 / 126 | 61 | 2.5e3 | 1.6e5 | **949 (90)** |
| 201 | 1.184 | 5.7e-3 | 0.17 | 4.7e-4 | 0.18 / 0.18 | 2.1e3 / 185 | 430 | 1.5e4 | 5.1e5 | **726 (566)** |
| 301 | 1.266 | 6.3e-3 | 0.21 | 3.9e-4 | 0.22 / 0.22 | 2.2e3 / 430 | 1.4e3 | 4.8e4 | 7.8e5 | 3.0e3 (6.2e3) |

Bold in the last column marks the states where P2_CONTRACT's full-range series failed. **Every one of them passes in window form.**

**Reading.**
1. **The computed part (win + bd) is ≤ 0.49 H everywhere** (0.34 at c = 2). The O(1/m⁴…⁸) scalar window costs are ≤ 5·10⁻⁵ H, with
   ratio q ≤ 6.5·10⁻³. The boundary term dominates and grows with J. It is an artefact of the sharp window and is computed, not assumed.
2. **The complement, with the plain scalar-calibrated edge envelope, is ≤ 9·10⁻³ H.** The tolerated multiple is Λ_AG ≥ 114.
3. **The window premises are weak:**
   * the spinning O(1/m⁴) coefficient is tolerated up to K₄w = 20.8 (J = 26, gap) … 4.8·10⁴ (J = 301), in Theorem W's 4a/4b units.
     P2_SPINNING's model gives ≤ 3.8, and P2_CONTRACT's O(1/m²) block/kernel constants are 1–2× scalar (Λ_B ≈ 2);
   * equivalently, any real-axis window remainder up to ε_w ≥ 2·10³ × the scalar's sup\|S₄/m⁴\| is tolerated.
4. **The binding premise is the +− edge amplitude near m → Λ.**
   * If the +− beyond-hybrid correction at ν = Λ carries the Ω^U helicity enhancement (1−x₁²)^{−2} relative to the scalar's, four
     states fail at Λ_AG = 1 (bold in the d/H column).
   * In intrinsic units the requirement is \|ρ_{+−}(x₁)\| ≤ ρ_tol, the beyond-hybrid correction relative to the +− flat envelope there.
     ρ_tol = 0.15 / 0.26 / 0.41 / 0.63 at (26, 30, 34, 38) × ≈ 1.03Λ, against the scalar's 0.27–0.69.

### 4.2 Sampled J\* (sign-free, θ = 1; P0, P1′ and W assumed; AG′ structure assumed)

| edge-amplitude premise | failing sampled states | **J\*** |
|---|---|---|
| plain scalar envelope × Λ_AG, Λ_AG ≤ 114 | none | **26** (all sampled states) |
| +− transferred envelope, Λ_AG = 1 | (26, 1.026), (30, 1.030), (34, 1.033), (38, 1.037) | **42** |
| transferred, Λ_AG = 3 | the above + (33, 1.033), (42, 1.041), (46, 1.045), (50, 1.049) | 51 |
| transferred, Λ_AG = 10 | all m ≈ 1.03–1.08Λ up to J = 80, plus (26, 30, 34, 1.13Λ) | 100 |
| transferred, Λ_AG = 30 | up to (120, 1.114) | 150 |
| \|ρ_{+−}(x₁)\| ≤ 0.7 (scalar-like) | J = 26 (all m), 30 (all m), (34, 1.033), (38, 1.037) | 42 |
| \|ρ_{+−}(x₁)\| ≤ 1 | the above + (33, 1.033), (34, 1.137), (42, 1.041); (34, 1.240) is marginal (ρ_tol = 1.0003) | 46 |
| transferred envelope, sampled masses m ≥ 1.13Λ only | none for Λ_AG ≤ 5.9 (first failure: (26, 1.128) at Λ_AG = 6.0) | 26 |

**Sensitivity to the rate threshold.**
* T1 = 40/κ: Λ_AG,tol min 127 / 0.32.
* T1 = 15/κ (primary): 114 / 0.29.
* T1 = 10/κ: 18.7 / 0.046. With this threshold the transferred Λ_AG = 1 row also fails (33, 42, 46, 50) × ≈ 1.03–1.05Λ.

The lj_scan rates support T1 ≈ 15/κ at x ≈ 0.86: the average rate over y ≤ 82 is ≥ 0.97 × the local value. They do not probe x₁ > 0.94.
The m ≈ 1.03Λ states have x₁ = 0.96–0.975, so there **(AG′)'s local rate is [M], not [NE].**

## 5. Theorem L′'s range (`p4w_L.py`)

**Setup.**
* No computed O(1/m²) exists at J > 344. The window O(1/m²) term is therefore L′'s premise P_L′(i) restricted to the window:
  Λ_amp × the scalar beyond-hybrid B₂, costed sign-free.
* Budget: the per-entry half-floor F/2, with F_cert [certified] and F_num quoted from LRINF_TIGHT §5.
* The edge envelope is the transferred (+−) one, at Λ_AG = 1.

| state | J | m/Λ | β | ν_w | F_cert/2 | window O(1/m²)/budget | O(1/m⁴) | bd | edge + top | **d/budget, Λ_amp = 3 / 10 (cert)** | d/budget (num) | K₄w tol | edge tol / scalar (c₂) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gap | 223 | 1.203 | 0.370 | 8.1 | 9.6e-6 | 1.2e-2 | 2.6e-7 | 7.3e-2 | 1.5e-3 | **0.086 / 0.28** | 0.048 | 7.9e5 | 1.6e11 / 2.5e8 |
| gap | 300 | 1.265 | 0.471 | 6.4 | 8.9e-6 | 7.0e-3 | 7e-8 | 4.4e-2 | 4.0e-4 | 0.052 / 0.17 | 0.039 | 3.8e6 | 4.2e11 / 1.6e8 |
| gap | 500 | 1.415 | 0.695 | 4.3 | 7.5e-6 | 3.2e-3 | 7e-9 | 2.0e-2 | 1.1e-4 | 0.023 / 0.078 | 0.021 | 4.5e7 | 2.1e12 / 8.4e7 |
| gap | 857 | 1.648 | 1.000 | 3.0 | 6.5e-6 | 1.5e-3 | 9e-10 | 8.8e-3 | 9.2e-5 | 0.010 / 0.034 | 0.009 | 5.0e8 | 1.1e13 / 6.0e7 |
| gap | 10⁴ | 4.585 | 3.05 | 0.98 | 2.8e-6 | 1.6e-4 | 1e-9 | 6.7e-4 | 1.3e-4 | 9.6e-4 / 2.9e-3 | 7.8e-4 | 8.9e11 | 1.0e16 / 1.0e8 |
| gap | 10⁵ | 14.19 | 5.31 | 0.57 | 8.4e-7 | 8.3e-5 | 8e-10 | 2.4e-4 | 1.3e-4 | 4.5e-4 / 1.2e-3 | 3.7e-4 | 2.3e13 | 7.0e17 / 2.6e8 |
| gap | 10⁶ | 44.76 | 7.60 | 0.39 | 2.1e-7 | 7.3e-5 | 7e-10 | 2.0e-4 | 1.1e-4 | 3.8e-4 / 1.0e-3 | 2.8e-4 | 1.2e14 | 1.3e19 / 5.2e8 |
| 1.1 m_gap | 250 | 1.347 | 0.370 | 8.1 | 9.6e-6 | 9.0e-3 | 2e-7 | 5.7e-2 | 5.7e-4 | 0.066 / 0.22 | 0.037 | 1.3e6 | 6.8e10 / 4.0e7 |
| 1.3 m_gap | 307 | 1.652 | 0.371 | 8.1 | 9.6e-6 | 5.9e-3 | 7e-8 | 3.7e-2 | 2.2e-4 | 0.043 / 0.14 | 0.024 | 3.0e6 | 3.7e10 / 8.1e6 |
| 1.5Λ | 279 | 1.5 | 0.371 | 8.1 | 9.6e-6 | 7.2e-3 | 1e-7 | 4.5e-2 | 3.2e-4 | 0.053 / 0.17 | 0.029 | 2.0e6 | 4.6e10 / 1.5e7 |
| 2Λ | 372 | 2.0 | 0.371 | 8.1 | 9.6e-6 | 3.9e-3 | 3e-8 | 2.5e-2 | 1.3e-4 | 0.029 / 0.096 | 0.016 | 6.5e6 | 3.1e10 / 3.9e6 |
| 5Λ | 929 | 5.0 | 0.370 | 8.1 | 9.6e-6 | 6.1e-4 | 8e-10 | 3.9e-3 | 4.7e-5 | 0.0045 / 0.015 | 0.0025 | 2.6e8 | 3.1e10 / 1.4e6 |

(J = 250, 344, 400, 700 at the gap are in `cache/L.txt`. They are monotone between the rows shown.)

**Reading [DV].**
* On L′'s whole sampled range (J\* = ⌈(0.37211m − 1)/2⌉ … 10⁶, m from the gap to 5Λ), the window-form P4 together with the window O(1/m²)
  term uses **ε_R ≤ 0.086 of the certified half-floor** at Λ_amp = 3. The value is 0.28 at Λ_amp = 10, and is largest at J\* = 223.
* The edge, even with the transferred envelope, is ≥ 600× below tolerance, because x₁ = Λ/m ≤ 0.83 there.
* At large J everything is tiny (≤ 10⁻³). The boundary term ∝ ψ₂(p_w)/κ, with p_w = 3/(βΛ) → 0.
* **Theorem L′'s P_L′(ii) (the O(1/m⁴) and beyond-all-orders remainder ≤ ε_R) therefore holds with ε_R ≤ 0.09, given (W) and (AG′).** The edge
  premise is not binding on L′'s range.

## 6. The beyond-all-orders piece (AG): what is assumed, and what the window removes

**Assumed (AG′), precisely.** For each entry e, g_e := exact_e − flat_e (de-shifted, per unit leading image) satisfies the following on
Q = [p_w, 1] × [0, ±25/(0.97Λβ)]:
* it is a sum of two functions g_±, analytic in Q_± and continuous on its closure;
* \|g_±(p ± iy)\| ≤ M(p)·exp(−0.97Λβ y/√(1−x²)) for y ≤ 15/κ(p), and the uniform rate 0.97Λβ beyond;
* M(1) ≤ M₁ (EDGE).

**Status.**
* **Structure, to all orders of the travelling-wave saddle expansion: [D].** HLJ_ANALYTIC §3: non-degenerate saddles on all of
  [3, Λ], decay ≥ 1.002β; the graviton has the same exponent and saddles.
* **Local rate: [NE]** at x ≤ 0.86 (lj_scan vertical lines, 9 states, §2.2). **[M]** at x₁ > 0.94, i.e. the m < 1.06Λ states.
* **Beyond all orders: [A].** This is H_AG. HLJ's randomised critical-point search found no other low-frequency critical point [NE + M].
* **Amplitude M₁:**
  * scalar: [DV] at 4 states (\|G\| ≤ 0.65 e_sc);
  * graviton: [M]. The helicity transfer (1−x²)^{−2} for +− is the Ω^U convention, not derived.

**What the window removes.**
1. **In the window, no AG is used at all.**
   * The window term needs only the real-axis sup bound (W). Its tolerance is ε_w ≥ 2·10³ × the scalar O(1/m⁴) size, i.e. ≥ 10⁻⁷
     relative at J = 344 and ≥ 10⁻³ at J = 26.
   * Any beyond-all-orders term in the window is covered by that sup bound.
   * The exact scalar shows the series and the exact image agree in the window to the data noise (≈ 10⁻⁹).
2. **The fixed-β series is never used outside the window.** The former premise "the O(1/m⁴) remainder costed to ν ≈ m is ≤ K₄" is gone.
   It was false as a statement about the series (§1).
3. **On the complement, P4 and AG merge into one hypothesis on the exact function.** There is no longer a separate "O(1/m⁴) remainder" and
   "beyond-all-orders" premise. Its only quantitative content is the edge amplitude M₁ at the single point ν = Λ (and the top side, which
   is suppressed by e^{−25}).
4. **The monopole is in the window.** μ = i/2 is a small-ν point.
   * At O(1/m²) it is 0 (computed: block −12.8, kernel +12.8).
   * At O(1/m⁴) it is tolerated up to \|Q₄\| ≤ θH m⁴/(2\|D\|), which is ≥ 9.7·10⁸ (J = 344) … 1.6·10¹² over the 47 states. For
     comparison, the O(1/m²) block monopole is 12.8.

## 7. Updated theorem statements

**Theorem W (window form).** Let O be a heavy primary in H_LJ(1000) with couplings λ ∈ ℂ² (even J) or ℂ¹ (odd J). Assume:
* **(P0)** the hybrid floor f_h (unchanged; leading-Regge part [P + certified] for b̂ ≥ 245 by LRINF_TIGHT, [NE] below);
* **(P1′)** the O(1/m) part of M is a congruence ([DV] at leading-Regge kernel level, [D-sketch] beyond);
* **(W)** in the Regge window ν < 3/β: \|G_e − F₂^bh/m² − S₄^{sc}/m⁴\| ≤ ε_w, with ε_w at most the tabulated tolerance.
  Equivalently: the spinning O(1/m⁴) window coefficient is ≤ K₄w, with K₄w ≥ 20.8 at every sampled state;
* **(AG′)** as in §6;
* **(EDGE)** the edge amplitude of each entry is ≤ Λ_AG,tol × e_sc(x₁) (table §4.1).

Then M_Ω(O) ⪰ 0 for every coupling vector, at every sampled state with J ≥ 26 (m ≤ 1.3Λ), with sign-free margin d/H ≤ 0.49 at Λ_AG = 1.

**Grade.** [P given P0, P1′, W, AG′, EDGE] (Lemma P4W + the 2×2 criterion) + [DV] (sampled constants). Between sampled states: [NE].

**Changes.**
* **P4 is replaced by (W) ∧ (AG′) ∧ (EDGE).**
* **The full-range d₂ is no longer used:** the O(1/m²) correction enters only through its window cost and the boundary amplitude.
* **The k = 3 piece is included** (INDEP_CHECK2 D2): in the window with the same structures, sign-free; in the complement through ψ₃.

**Theorem W (window form, Ω^U +− transfer).** The same, with (EDGE) read as \|ρ_{+−}(x₁)\| ≤ ρ_tol. Sampled **J\* = 42** if the graviton's +−
edge correction is ≤ the scalar envelope × (1−x₁²)^{−2} (Λ_AG = 1). J\* = 42 / 46 if \|ρ_{+−}\| ≤ 0.7 / 1. J\* = 51 at Λ_AG = 3.

**Theorem L′ (unchanged threshold).** J ≥ J\*(m) = ⌈(0.37211m − 1)/2⌉ (223 at the twist edge). Premise P_L′(ii) is now discharged with
ε_R ≤ 0.086 (Λ_amp = 3) at every sampled L′ state, given (W) and (AG′). The premises are now:
* **P_L′(i)** (O(1/m²) magnitudes, Λ_amp ≤ 181; window part only needed);
* **(AG′)**;
* **(EDGE)**, not binding: ≥ 600× margin.

**Grade:** [D given P_L′(i), AG′, EDGE]. The floor is [P + certified]; the costs are [DV] at 16 sampled states.

## 8. Honest residual

1. **(AG′) beyond all orders is assumed [A].** This is H_AG. It is now the only structural hypothesis on ν ≥ 3/β, and it is formulated on
   the exact function. Its local-rate part is [NE] at x ≤ 0.86 and [M] at the near-gap edge x₁ ≳ 0.95.
2. **The +− edge at m ≲ 1.05Λ, J ≤ 50 is where P4 does not close.** The remaining inequality, per state, is

       4[V(1; κ₁)(M₁ + A) + V(1; 0.97Λβ)A] + Top ≤ θH − (window + boundary),

   i.e. M₁ ≤ M_tol = 1.3·10¹⁰ c₂ at (26, 1.026Λ) (table §4.1, "edge tolerance"). If M₁ = ρ × (the +− flat envelope there), the
   requirement is ρ ≤ 0.15 (J = 26), 0.26 (30), 0.41 (34), 0.63 (38) at m ≈ 1.03Λ. The exact scalar has ρ ≈ 0.27–0.69 at its own edges.
   * Nothing is known about the graviton there.
   * The floor at those states is itself only [NE].
3. **Spinning O(1/m⁴) symbols and the kernel's O(r²) term: not computed.** They enter only the window (tolerance K₄w ≥ 20.8, against the
   model's 3.8 and the measured O(1/m²) ratio ≈ 2) and the O(1/m⁴) monopole (tolerance \|Q₄\| ≥ 9.7·10⁸). Grade of (W) for the graviton: [M] with
   margin ≥ ×5 (J = 26) to ×10⁴ (J ≥ 250).
4. **Boundary amplitude M_w** (factor 2 over the real-axis sup near ν_w). It covers the Harish-Chandra wave split at νβ ≈ 3 (≤ 1.2) and the
   polynomial growth of ν² classes along the vertical side (≤ 1.25). That estimate is [D-sketch]. The boundary term is ≤ 0.44 H, so a
   factor-2 error here would push J = 300, 344 to d/H ≈ 0.9 (still < 1). c = 2 gives ≤ 0.34 H.
5. **Everything is float [DV], not certified.** The exact scalar data are Arb, but the costs are GL sums. A certification would need Arb
   evaluation of V(p₀; κ) (elementary integrals of polynomial moduli) and of the window integrals. This is straightforward but was not done.
6. **Sampled states only.** 47 W states plus 16 L′ states. Between them [NE]: smooth trends in every column.
7. Not covered by either theorem:
   * J < 26;
   * near-gap non-sampled 64 < J < 223;
   * m > 1.3Λ for W.
8. **Scalar → graviton transfer of the envelope shape** (bulk ∝ (νβ)², edge ∝ (1−x²)^{−3}) is [M]. Only its value at p = 1 matters for the
   bound, and the bulk shape enters only the e^{−25}-suppressed top side.

## 9. Suggested inventory changes (not applied)
* **P2_CONTRACT §2.2 / §5 item 1:**
  * "P4 [I]: the series remainder exceeds tolerance at J ≥ 80" → "the full-range fixed-β series is asymptotic and must not be costed beyond
    x ≈ 0.3. The exact scalar remainder costs ≤ 1.5·10⁻³ of the budget (K₄-equivalent ≤ 8 against the series' 176 at J = 100). P4 is
    reformulated as (W) ∧ (AG′) ∧ (EDGE) (P4_WINDOW)";
  * residual 2 (3/16): resolved [DV + NE], a genuine beyond-hybrid O(1/m⁴) term of the exact moment.
* **LARGESPIN_W / UNIFIED_THEOREM H_LJ(1000):** [D given P4] → **[P given P0, P1′, W, AG′, EDGE] + [DV] on all sampled states J ≥ 26 with the
  plain edge envelope; J ≥ 42 with the Ω^U +− transfer.** **Do not write "proved for J ≥ 26 / 42"**: AG′ beyond all orders is [A], and
  the edge amplitude and the spinning O(1/m⁴) window coefficients are [M].
* **Theorem L′:** P_L′(ii) holds with ε_R ≤ 0.09 given (W) and (AG′).
* **HLJ_ANALYTIC Lemma C:** extended with the curvature-wave term A and the local-rate threshold (this note's (2.1)). Its "E_edge" uses the
  local rate, which is load-bearing at m → Λ.

## 10. Reproduce

    cd ads4/code/p4win && bash run_all.sh <path/to/regge_moments_order8.m>     # ≈ 0.7 CPU-h, ≤ 2 processes

| file | content | grade |
|---|---|---|
| `cache/series.txt` | g₁ … g₄ (flat parts of CHMRS's orders 1/m² … 1/m⁸) | [Q + P] |
| `cache/exact.txt` | §1: exact vs series, full-range and window costs, envelope bins | [DV] |
| `cache/check.txt` | §3.1–3.2: envelope check, Lemma P4W on the exact scalar | [DV] |
| `cache/p316.txt` | §3.3: the 3/16 residual | [DV + NE] |
| `cache/states15.txt` (10, 40) | §4: 47 states, rate-threshold sensitivity | [DV] |
| `cache/L.txt` | §5: Theorem L′ range | [DV] |
