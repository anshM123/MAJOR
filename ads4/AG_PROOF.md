# AG_PROOF: the premise (AG′) of Lemma P4W — what can be proved, a weaker premise that suffices, and why (AG′) itself cannot be proved by the suggested routes

*23–24 Nov 2026. Unreviewed, not committed. No existing file was edited.*
* *New code:* `ads4/code/agproof/` (`run_all.sh`; ≈ 1.6 CPU-h, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`). Caches in
  `code/agproof/cache/`, whose `.gitignore` ignores everything except itself.
  * `ag_flat.py`: Lemma F (the flat partial wave as an exact pair of analytic travelling waves) and Proposition N1 (AG′'s majorant against it).
  * `ag_unif.py`: the uniform-rate and amplitude bounds of Lemma F on Lemma P4W's rectangle.
  * `ag_exact.py`: the **exact** scalar (Δφ = 3) k = 2 heavy moment at **complex ν** (first computation; Arb, entire Laplace/Barnes series).
  * `ag_lock.py`: test of a "flat-locked" two-wave model for exact − flat near the edge (negative, §3.3).
  * `ag_edge.py`, `ag_edge2.py`: the exact scalar's **edge contribution**, measured on the real axis with an entire cutoff (§4, §5.3).
  * `ag_chi.py`: **Lemma P4W-χ** (§5) on Theorem W's 47 states.
  * `ag_p4w.py`: Lemma P4W's tolerances under weaker or re-shaped decay premises (§6), 47 W states and 16 L′ states.
  * `ag_pwmeasure.py`: the part of the heavy factor's plane-wave measure beyond the twist edge (§7.1).
* *Read-only inputs:* `p4win/` (p4w_common, p4w_contour, p4w_states, p4w_L), `p2con/` caches, `exactskel/` (h13_exact → Laplace3H, Mack
  polynomial disk cache), `spinfunc/cache/sm_exact_1000.pkl` (real-axis cross-check only).

Grades as in UNIFIED_THEOREM §0. Conventions as P4_WINDOW: p = ν/Λ, x = ν/m, x₁ = Λ/m, β = 2 asinh((J+½)/m), f = de-shifted image per unit
leading image; flat = E₂(x)P_J(1 − 2x²)/2; g := exact − flat; G := exact − hybrid = g − (P − J₀); κ₀ := 0.97Λβ and κ_loc := κ₀/√(1−x₁²)
(p-units); T₀ := 25/κ₀ is P4_WINDOW's rectangle height.

---

## 0. Verdict

### 0.1 What is and is not proved

| # | claim | grade | constants | coverage |
|---|---|---|---|---|
| **V0** | **(AG′) is not proved.** No route available in the record gives the travelling-wave split of the exact (scalar or graviton) correction beyond all orders. | — | — | — |
| **N1** | **(AG′)'s pointwise majorant is violated by the explicit flat wave itself at the near-gap edge.** On the right side p = 1, \|flat₊(1+iy)/flat₊(1)\| exceeds exp(−κ_AG′(y)) by up to **×4.1** (J = 26, m = 1.026Λ, at y = 15/κ_loc); phase-only by ×40. Cause: the branch point of θ(ν) = 2 asin(ν/m) at ν = m, at distance (m − Λ) ≈ 26–49 from the edge, bends the profile from linear to ∝ √y once y ≳ m − Λ. **The integrated consequence survives:** the edge integral with the true flat profile is **0.72–0.82** of AG′'s at every one of the 67 sampled states (47 W, 4 exact, 16 L′). | [P] (Lemma F) + [DV] | table §2 | x₁ ≥ 0.95: violated; x₁ ≤ 0.94: holds |
| **N2** | **The plane-wave (spin-Gamma) route cannot prove a near-edge statement.** A positive ν-independent average of AG-class functions is AG-class with the *worst* rate of its support (Lemma R, [P]), but for J ≤ 50 near the gap **14–35 %** of the moment-weighted positive plane-wave measure has \|p\| < Λ, i.e. lies *beyond its own twist edge* inside [p_w, 1]. Those terms have no travelling-wave form there; the state's edge behaviour comes from cancellations between terms (the N-L1 mechanism). | Lemma R [P]; fraction [DV] (MC, 4·10⁵ samples) | table §7.1 | J ≤ 64 near the gap; harmless for J ≥ 100 (≤ 0.6 %) |
| **F** | **Lemma F: the flat partial wave is an exact pair of analytic travelling waves.** flat = flat₊ + flat₋, flat_± = E₂W_±(θ(ν))/2, with W₊(θ) = (2/πi)C_J e^{i(J+1)θ}₂F₁(½, J+1; J+3/2; e^{2iθ}), and **\|W₊(θ)\| ≤ (2/π)C_J e^{−(J+1)Im θ} min{2/\|sin Re θ\|, (1−e^{−2Im θ})^{−½}}**. It satisfies the uniform-rate premise on the whole rectangle (ratio ≥ 1.013 at every state) and has the local rate (J+½)·2/(m√(1−x²)) at y → 0. | **[P]** (DLMF 14.13.1 + Abel summation); checks [DV] to 3·10⁻²⁹ | §1 | all J, m > Λ |
| **χ** | **Lemma P4W-χ: (AG′) ∧ (EDGE) can be replaced by (AG_U) ∧ (EDGE-R).** An entire cutoff χ(p) = ½erfc((p_e − p)/σ) splits the complement cost into a **real-axis edge integral** E_χ[G] = ∫_{p_w}^1ψ₂χG (premise EDGE-R) and a Cauchy part with weight ψ₂(1−χ), which needs only **the uniform rate κ₀ — no local rate, no threshold, anywhere**. The edge amplitude enters that part suppressed by ½erfc((1−p_e)/σ). | **[P]** (Cauchy + \|erfc(a+ib)\| ≤ e^{b²}erfc(a)); constants [DV] | §5 | every state |
| **E** | **The exact scalar's true edge cost is 6–100× below what (AG′) + e_sc charges for it.** First evaluation of the exact scalar moment at complex ν (Arb; real axis reproduces sm_exact exactly). The edge contribution, isolated on the real axis by an entire cutoff (method validated on the flat wave against Lemma F's exact contour integral, ≤ 0.3 %), is **0.11 / 0.09 / 0.15 / 0.010** of P4W's edge term at J = 40 / 48 / 64 / 100 (x₁ = 0.962–0.913), i.e. **0.17 / 0.14 / 0.22 / 0.014 S_E** (S_E = the cost of a flat-profile wave of amplitude e_sc(1)). With Lemma P4W-χ's own cutoffs: {{EDGE2}}. | [DV] (Arb values, float quadrature; stable to ±5 % over 5 cutoffs) | §4, §5.3 | scalar, x₁ = 0.91–0.96 — **the near-gap edge where (AG′) was [M]** |
| **W-χ** | **Theorem W with (AG_U) ∧ (EDGE-R) instead of (AG′) ∧ (EDGE):** plain envelope, everything except E_χ uses **≤ {{DCHI_P}} H** and the tolerated edge cost is **≥ {{TOLE_P}} S_E** at every sampled state (exact scalar ≤ 0.22 S_E): **sampled J\* = 26**. With the Ω^U +− transfer: {{WCHI_TR}}. | [P given P0, P1′, W, AG_U, EDGE-R] + [DV] | §5.2 | 47 sampled states |
| **L-U** | **Theorem L′ needs only the uniform rate.** {{LU}} | [D given P_L′(i), AG_U, EDGE] + [DV] | §6.2 | 16 sampled L′ states |
| **θ** | **Weakest rate that still closes Lemma P4W (AG′'s shape, rate θβ/√(1−x²) without threshold):** {{THETA}} | [DV] | §6.1 | 47 W states |

### 0.2 Bottom line
* **(AG′) is not proved, and in its stated pointwise form it is not even true for the function it was modelled on** (N1): the flat partial
  wave's decay along the right side bends below the assumed local rate at y ≳ 10/κ_loc for x₁ ≥ 0.95. Only its integrated consequence (the
  edge integral) is right, with 20–28 % to spare. The 15/κ-vs-10/κ "threshold sensitivity" of P4_WINDOW is an artefact of this wrong shape.
* **The near-gap edge is not reachable through the plane-wave representation** (N2): at J ≤ 50 the representation's terms do not share the
  state's twist edge. No positivity argument can produce the local rate there.
* **What replaces it is weaker and better founded** (Lemma P4W-χ, [P]):
  * **(AG_U)** the uniform rate κ₀ = 0.97Λβ on the complement. Every explicit piece satisfies it [P] (flat waves: Lemma F, ratio ≥ 1.013;
    curvature waves: Harish-Chandra / Hankel, rates Λβ and b̂). For the exact correction it is [D] to all orders of the saddle expansion
    (HLJ: ≥ 1.002β) and [A] beyond all orders. It no longer contains any statement about the near-gap edge.
  * **(EDGE-R)** one real-axis number per entry: the edge integral E_χ[G_e]. It replaces both the [M] local rate at x₁ > 0.94 and the
    [M] amplitude premise (EDGE). It is **directly computable from real-axis data** (no continuation, no rate), and it is [DV] for the
    scalar analogue with margin ≥ {{MARGIN}}.
* **Coverage after this note:**
  * Theorem W, plain envelope: sampled J\* = 26 under (P0, P1′, W, AG_U, EDGE-R).
  * Theorem W, Ω^U +− transfer: {{WTR_SHORT}}
  * Theorem L′ (J ≥ ⌈(0.37211m − 1)/2⌉): {{LSHORT}}
* **What remains:** (AG_U) beyond all orders [A]; (EDGE-R) for the graviton [M] (scalar [DV]); the rest of Theorem W/L′'s premises
  unchanged (P0, P1′, W, the boundary factor 2 [D-sketch], sampled states, float).

---

## 1. Lemma F: the flat partial wave as a pair of analytic travelling waves [P]

### 1.1 Statement
For an integer n ≥ 0 let
* a_k = (½)_k(n+1)_k/(k!(n+3/2)_k), so 0 < a_{k+1} < a_k ≤ (½)_k/k!;
* F_n(q) = Σ_k a_k q^k = ₂F₁(½, n+1; n+3/2; q) (|q| < 1, continuous to |q| = 1, q ≠ 1);
* C_n = √π Γ(n+1)/(2Γ(n+3/2)) = 2^{2n}(n!)²/(2n+1)!;
* W₊(θ) = (2/(πi)) C_n e^{i(n+1)θ} F_n(e^{2iθ}),  W₋(θ) = −(2/(πi)) C_n e^{−i(n+1)θ} F_n(e^{−2iθ}).

Then for 0 < Re θ < π:
1. **P_n(cos θ) = W₊(θ) + W₋(θ)**, W₊ analytic for Im θ > 0, W₋ for Im θ < 0, both continuous up to the real interval (0, π) and each
   continued across it by W_± = P_n − W_∓.
2. **|W₊(θ)| ≤ (2/π) C_n e^{−(n+1)Im θ} · min{ 2/|sin Re θ|, (1 − e^{−2 Im θ})^{−½} }** for Im θ ≥ 0.
3. With θ(ν) = 2 asin(ν/m) (analytic on |Re ν| < m, mapping the upper half-plane into Im θ > 0 near the segment), **flat_± := E₂(x)W_±(θ(ν))/2**
   (n = J) split the flat image into two analytic travelling waves on Lemma P4W's rectangles Q_± (whose right side p = 1 has Re ν = Λ < m).
4. **Profile.** Im θ(ν_r + iy) ≥ 2 asinh(y/m) for every |ν_r| < m, so |flat₊| decays at least like e^{−(J+1)·2asinh(y/m)}. On Q_± (y ≤ T₀)
   this is **≥ the uniform rate κ₀**: min over the 67 sampled states and y ≤ T₀ of (J+1)2asinh(y/m)/(0.97βy) is **1.013** (J = 26, m = 1.231Λ).
   At y → 0 the rate is the local one, (J+½)·2/(m√(1−x²)) ≈ β/√(1−x²).

### 1.2 Proof
1. DLMF 14.13.1 with μ = 0 gives the Fourier series P_n(cos θ) = (4/π)C_n Σ_k a_k sin((n+2k+1)θ) for 0 < θ < π, which is
   (2/(πi))C_n[e^{i(n+1)θ}F_n(e^{2iθ}) − e^{−i(n+1)θ}F_n(e^{−2iθ})]. The series F_n converges absolutely for |q| < 1. On |q| = 1, q ≠ 1, it
   converges by Dirichlet's test (a_k decreasing to 0) and Abel's theorem gives continuity. Both sides are analytic in θ, so the identity
   continues.
2. *Second entry of the min:* a_k ≤ (½)_k/k! because (n+1)_k ≤ (n+3/2)_k, so |F_n(q)| ≤ (1 − |q|)^{−½}. *First entry:* Abel summation
   with a_k decreasing gives |Σa_kq^k| ≤ a₀ sup_N|Σ_{k<N}q^k| ≤ 2/|1 − q|, and |1 − re^{iφ}| ≥ |sin(φ/2)| for 0 ≤ r ≤ 1 (minimise over r).
3. asin is analytic off (−∞, −1] ∪ [1, ∞). E₂ has its only pole at x = ±1. Re ν ≤ Λ < m on the rectangles.
4. Im asin(x + iη) = log(α + √(α² − 1)) with α = (|z + 1| + |z − 1|)/2, which increases with |x| at fixed η, and equals asinh(η) at x = 0.

∎ **[P].** Checks (`ag_flat.py`, `ag_unif.py`, mpmath 30 digits): the identity to 3.2·10⁻²⁹ (n = 0…344, real and complex θ); bound 2 is
attained to 0.974 (second entry) and 0.948 (min); the monotonicity in 4 on a grid.

### 1.3 What Lemma F does and does not give
* It gives (AG′)-type structure **for the flat part only**, which is explicit. It is the first place in the record where the travelling-wave
  split of any piece of the near-edge function is rigorous rather than to-all-orders.
* Together with the Harish-Chandra waves of P_{−½+iν}(cosh β) (rate exactly β) and the Hankel waves of J₀(b̂p) (rate b̂/Λ ≥ β), **every
  explicit piece of the hybrid satisfies (AG_U) [P].**
* It says nothing about exact − flat. That is the whole content of (AG′).

---

## 2. Proposition N1: (AG′)'s majorant against the flat wave (`ag_flat.py`, `cache/flat.txt`) [P + DV]

(AG′) (P4_WINDOW §2.1, §6) assumes, on the right side p = 1,
|g₊(1 + iy)| ≤ M₁·exp(−κ_AG′(y)), with κ_AG′(y) = κ_loc·y for y ≤ 15/κ_loc and 15 + κ₀(y − 15/κ_loc) beyond.

The flat wave is the model this was abstracted from ("the local frequency of the flat partial-wave phase", p4w_contour docstring). Its own
profile on p = 1 is exact by Lemma F. Columns at y = c/κ_loc: −log|flat₊(1+iy)/flat₊(1)| (amplitude included) | phase only,
(J+½)Im[θ(1+iy) − θ(1)] | AG′'s exponent. Last columns: the edge integrals V = ∫₀^{T₀}|ψ₂(1+iy)|·profile dy, true flat (amplitude
included) / phase only, relative to AG′'s; and the worst pointwise ratio.

| state (J, m/Λ) | x₁ | c = 5 | c = 10 | c = 15 | V_flat/V_AG′ | V_phase/V_AG′ | max_y ratio (at y = 15/κ) |
|---|---|---|---|---|---|---|---|
| 26, 1.026 | 0.975 | 5.4 \| 4.8 \| 5.0 | 10.0 \| 8.4 \| 10.0 | **13.6 \| 11.3 \| 15.0** | 0.72 | **2.55** | **4.1** (phase: 40) |
| 30, 1.030 | 0.971 | 5.3 \| 4.9 \| 5.0 | 10.1 \| 8.9 \| 10.0 | 13.9 \| 12.0 \| 15.0 | 0.74 | 1.72 | 2.9 |
| 34, 1.033 | 0.968 | 5.3 \| 5.0 \| 5.0 | 10.1 \| 9.2 \| 10.0 | 14.2 \| 12.6 \| 15.0 | 0.75 | 1.36 | 2.2 |
| 38, 1.037 | 0.964 | 5.3 \| 5.0 \| 5.0 | 10.2 \| 9.4 \| 10.0 | 14.4 \| 13.1 \| 15.0 | 0.76 | 1.18 | 1.8 |
| 42, 1.041 | 0.960 | 5.2 \| 5.0 \| 5.0 | 10.2 \| 9.6 \| 10.0 | 14.6 \| 13.5 \| 15.0 | 0.77 | 1.07 | 1.5 |
| 50, 1.049 | 0.953 | 5.2 \| 5.1 \| 5.0 | 10.3 \| 9.8 \| 10.0 | 14.9 \| 14.1 \| 15.0 | 0.79 | 0.96 | 1.1 |
| 64, 1.062 | 0.942 | 5.2 \| 5.1 \| 5.0 | 10.3 \| 10.1 \| 10.0 | 15.2 \| 14.7 \| 15.0 | 0.80 | 0.88 | 1.0 |
| 100, 1.096 | 0.913 | 5.2 \| 5.2 \| 5.0 | 10.3 \| 10.3 \| 10.0 | 15.4 \| 15.3 \| 15.0 | 0.81 | 0.83 | 1.0 |
| 344, 1.299 | 0.770 | 5.2 \| 5.2 \| 5.0 | 10.4 \| 10.4 \| 10.0 | 15.6 \| 15.6 \| 15.0 | 0.78 | 0.78 | 1.0 |
| L′: 223 … 10⁶ | ≤ 0.83 | ≥ AG′ | ≥ AG′ | ≥ AG′ | ≤ 0.80 | ≤ 0.80 | 1.0 |

Max over all 67 states: V_flat/V_AG′ = 0.821 (k = 2), 0.821 (k = 3); V_phase/V_AG′ = 2.55; pointwise ratio 4.1.

**Reading.**
* **(AG′) as a pointwise statement fails for the flat wave at x₁ ≳ 0.95** (the m ≲ 1.05Λ states, J ≤ 51), beyond y ≈ 10/κ_loc. The reason
  is analytic, not numerical. The phase (J+½)θ(ν) has its branch point at ν = m, a distance m − Λ = 26 (J = 26) … 49 (J = 50) from the
  corner p = 1. Once y exceeds that distance, Im θ grows like √y, not linearly. At J = 26: m − Λ = 26 while 15/κ_loc = 67.
* **The amplitude rescues the integral.** |E₂| decreases away from its pole at ν = m along the same side, so the true flat profile
  (amplitude included) gives V within 0.72–0.82 of AG′'s at every state. A function that decays like the flat *phase* with a *constant*
  amplitude would exceed AG′'s edge term by up to 2.55.
* **Consequence for P4_WINDOW's threshold sensitivity** (T1 = 15/κ → 10/κ lowers the tolerance 114 → 18.7): the piecewise majorant is the
  wrong shape near the edge, so its threshold is not a meaningful parameter. The family 'Fa' of §6 (the exact flat profile, no
  threshold) is the natural replacement and gives the same or better tolerances than T1 = 15.

---

## 3. The exact scalar at complex ν (`ag_exact.py`, `ag_lock.py`) [DV / NE]

### 3.1 Method
exactskel's Arb Laplace/Barnes form (`h13_laplace.Laplace3H`, used read-only through `h13_exact`) evaluates the k = 2 moment as
m⁴C_{2,ν} = pre·Σ_{m<M_x} g_m(ν)R_m, g_m(ν) = Π_{l<m}((A/2 + l)² + ν²/4)/((A+l)(l+1))/Γ(A), with R_m independent of ν. **Each g_m is a
polynomial in ν², so the same series is the entire continuation of the exact moment.** `ag_exact.py` evaluates it in acb arithmetic
(400 bits). Checks:
* real axis: identical (relative difference 0.0 in double) to `sm_exact_1000.pkl` at 13 nodes per state;
* truncation: doubling M changes nothing in double (J = 40, 100);
* Arb radii ≤ 3.6·10⁻⁴⁴ relative on all lines;
* cost: 10–20 s per state (J = 40, 48, 64, 100 at the gap; m/Λ = 1.039, 1.047, 1.062, 1.096; x₁ = 0.962, 0.955, 0.942, 0.913).

### 3.2 Growth on vertical lines
On ν = pΛ + iy the function is dominated, for y ≳ 2/κ, by its growing wave g₋. For a real-analytic two-wave function the growth of g₋
upward mirrors the decay of g₊, up to the amplitude factor |a(ν+iy)a(ν−iy)|. |g/flat| along the lines (`cache/exact_lines_*.pkl`):

| J | line | y·κ_loc ≈ 1 | 3 | 6 | 9 | 12 | 18 | top (T₀) |
|---|---|---|---|---|---|---|---|---|
| 40 | p = 1 | 0.93 | 0.36 | 0.20 | 0.14 | 0.11 | 0.08 | 0.056 |
| 40 | p = 0.8 | 0.025 | 0.028 | 0.031 | 0.033 | 0.033 | 0.034 | 0.033 |
| 64 | p = 1 | 3.5 | 2.2 | 1.4 | 0.95 | 0.70 | 0.45 | 0.31 |
| 100 | p = 1 | 9.3 | 10.4 | 9.5 | 7.7 | 6.0 | 3.5 | 1.9 |

* arg(g/flat) drifts by only 1–3 rad over the full height (against a flat phase of 40–55 rad), so the exact correction's growing wave
  is **phase-locked to the flat wave to a few %** in the bulk and at the edge.
* **At the edge, |g/flat| decreases upward** (×17 at J = 40). Two readings are consistent with this: (i) g's phase is slower than the
  flat's (which would violate the local rate), or (ii) the *relative amplitude* of g decreases away from the twist edge like
  |1 − x²|^{−(1.6…2.5)} (the edge enhancement of e_sc), which is harmless. The vertical-line data alone cannot separate them.

**Growth rates** (`ag_growth.py`, `cache/growth.txt`): the mean slope of log|g| (and of log|flat|) over y ∈ [5, 15]/κ_loc(ν_r), where the
growing wave dominates by ≥ e¹⁰. For a real-analytic two-wave function this is the decay rate of the + wave, up to the slowly varying
amplitude factor |a(ν+iy)a(ν−iy)| [NE reading]. Entries: growth(g)/κ₀ | growth(g)/κ_loc | growth(flat)/κ_loc.

| J (x₁) | p = 1 | p = 0.95 | p = 0.9 | p = 0.8 | p = 0.6 | p = 0.4 |
|---|---|---|---|---|---|---|
| 40 (0.962) | 2.33 \| **0.64** \| 0.73 | 2.01 \| 0.82 \| 0.87 | 1.80 \| 0.90 \| 0.92 | 1.54 \| 0.98 \| 0.96 | 1.26 \| 1.03 \| 1.00 | 1.11 \| 1.03 \| 1.01 |
| 48 (0.955) | 2.48 \| **0.73** \| 0.82 | 2.09 \| 0.88 \| 0.91 | 1.84 \| 0.94 \| 0.95 | 1.56 \| 1.01 \| 0.98 | 1.26 \| 1.03 \| 1.01 | 1.11 \| 1.03 \| 1.01 |
| 64 (0.942) | 2.54 \| **0.86** \| 0.92 | 2.12 \| 0.95 \| 0.97 | 1.86 \| 0.99 \| 0.99 | 1.57 \| 1.04 \| 1.01 | 1.26 \| 1.04 \| 1.02 | 1.11 \| 1.03 \| 1.02 |
| 100 (0.913) | 2.44 \| 0.99 \| 1.00 | 2.05 \| 1.02 \| 1.01 | 1.80 \| 1.03 \| 1.02 | 1.56 \| 1.06 \| 1.02 | 1.24 \| 1.04 \| 1.03 | 1.11 \| 1.03 \| 1.03 |

* **In the bulk (x ≤ 0.86) the exact correction has AG′'s local rate** (0.98–1.06 κ_loc ≥ 0.97): this confirms lj_scan independently,
  with the exact function instead of the saddle expansion.
* **At the near-gap edge it does not:** 0.64 / 0.73 / 0.86 κ_loc at x₁ = 0.962 / 0.955 / 0.942 (and 0.99 at 0.913). The flat wave itself
  has 0.73 / 0.82 / 0.92 there (N1). **(AG′)'s local rate at x₁ > 0.94 is [NE]-contradicted for the scalar analogue**, at least in this
  averaged sense.
* **The uniform rate holds with margin everywhere:** growth(g) ≥ 1.11κ₀ on every line, ≥ 2.3κ₀ at the edge. This is the evidence for
  (AG_U) at the near-gap edge (§5).

### 3.3 Negative: no constructive split
`ag_lock.py` fits g = c̄·flat₋ + c·flat₊ (c, c̄ polynomials in ν, degree 8–12) by least squares on 300–400 exact complex values in a disc
at the edge (radius 78–125, avoiding ν = m). **It does not fit:** the residual is 0.3–5 times the local envelope at every degree and disc
tried. Real-analyticity (d = conj b) comes out to 10⁻¹³, as it must, but the model is wrong. So:
* exact − flat near the edge is **not** a flat-locked two-wave function with a slowly varying ratio, even though its growing wave is
  phase-locked on vertical lines;
* a numerical split g = g₊ + g₋ with a controlled decay cannot be constructed this way, and the extrapolated "g₊" profiles from these
  fits (V/V_AG′ = 1.4–3.3) are **not** evidence either way.

This is what forced the real-axis formulation of §4–5.

---

## 4. The exact scalar's edge contribution, measured on the real axis (`ag_edge.py`, `cache/edge.txt`) [DV]

**Method.** For f = f₊ + f₋ with analytic travelling waves, the edge (p = 1) part of ∫ψ₂f is, by Cauchy,
2Re[−i∫₀^∞ψ₂(1+iy)f₊(1+iy)dy]. It is isolated on the real axis by an entire cutoff:
E_χ[f] = ∫ψ₂(p)χ(p)f(p)dp, χ = ½erfc((p_e − p)/σ), with κ(p_e)σ ≥ 12 so that the cutoff contributes ≲ e^{−(κσ)²/4}. No continuation and
no split are needed. Real-axis values are the exact Arb moments; quadrature is composite Gauss–Legendre at ≥ 6 nodes per local wavelength.

**Validation on the flat wave (Lemma F gives the exact contour integral):**

| J | E_χ[flat] over 4–5 cutoffs | Lemma F right-side integral | agreement |
|---|---|---|---|
| 40 | −1.447 … −1.609·10⁻⁸ | −1.6128·10⁻⁸ | 0.25 % (p_e = 0.80) … 10 % (p_e = 0.90, narrow) |
| 48 | −1.048 … −1.084·10⁻⁸ | −1.0842·10⁻⁸ | 0.03–3 % |
| 64 | +1.684 … +1.710·10⁻⁹ | +1.7101·10⁻⁹ | 0.01–1.5 % |
| 100 | −5.453 … −5.466·10⁻¹¹ | −5.4657·10⁻¹¹ | ≤ 0.2 % |

**Exact scalar** (G = exact − hybrid; P4W's edge term at Λ_AG = 1 = 2[V(1;κ_loc, T1 = 15)(e_sc(1) + A) + V(1;κ₀)A], k = 2, per wave pair,
without the ring weight; S_E2 = 2V_flat(1)e_sc(1)):

| J | x₁ | \|E_χ[G]\| (5 cutoffs) | P4W edge term | **ratio** | in S_E2 | max\|G\| near p = 1 vs e_sc(1) |
|---|---|---|---|---|---|---|
| 40 | 0.962 | 4.04–4.38·10⁻⁸ | 3.80·10⁻⁷ | **0.11** | 0.17 | 26.2 vs 69.2 |
| 48 | 0.955 | 1.65–1.78·10⁻⁸ | 1.84·10⁻⁷ | **0.09** | 0.14 | 23.9 vs 58.3 |
| 64 | 0.942 | 8.84–9.25·10⁻⁹ | 5.99·10⁻⁸ | **0.15** | 0.22 | 14.9 vs 44.8 |
| 100 | 0.913 | 0.93–1.05·10⁻¹⁰ | 1.09·10⁻⁸ | **0.010** | 0.014 | 14.2 vs 30.3 |

**Reading.**
* **At the near-gap edge (x₁ = 0.91–0.96), where (AG′)'s local rate was only [M], the exact scalar's actual edge cost is 7–100× below
  what Lemma P4W charges with (AG′) and the scalar envelope.** About ×2–3 of this is e_sc(1) overestimating the amplitude; the rest is
  the true edge integral being smaller than the majorant.
* The test is insensitive to how the split would be defined: it measures the quantity the contour lemma bounds, on the real axis.
* E_χ[G] is stable to ±5 % over the cutoffs (±10 % at the narrowest), i.e. it is the edge contribution, not a cutoff artefact.

---

## 5. Lemma P4W-χ: the edge on the real axis, the complement at the uniform rate [P]

### 5.1 Statement
Fix a state (J, m), an entry e, G := G_e, the window edge p_w = 3/(βΛ), and parameters p_w < p_e < 1, σ > 0, 0 < T ≤ T₀. Let
χ(p) = ½erfc((p_e − p)/σ) (entire). Assume:
* **(W)** as in Lemma P4W (real axis, ν < 3/β);
* **(AG_U)** on Q_± = [p_w, 1] × [0, ±T]: G = G₊ + G₋, G_± analytic inside, continuous on the closure, and on the three sides
  |G_±(p ± iy)| ≤ M(p)e^{−κ₀y}, with M(p_w) ≤ M_w (the window series' waves, computed), M(1) ≤ M₁ + A, M ≤ M_top on the top. **No local
  rate is assumed anywhere;**
* **(EDGE-R)** |∫_{p_w}^1ψ₂χG dp| ≤ ε_{E,2} and |∫_{p_w}^1ψ₃χG dp| ≤ ε_{E,3} (real axis).

**Then**

    |cost(G)| ≤ [window terms of Lemma P4W] + 2ε_{E,2} + u ε_{E,3}
              + 4[ V_L M_w + V_R (M₁ + 2A) + Top_χ ]   (k = 2; k = 3 with weight u, ψ₃)  + 2|D||G(i/2)|,

with
* V_L = ∫₀^T|ψ₂(p_w + iy)|(1 + ½e^{y²/σ²}erfc((p_e − p_w)/σ))e^{−κ₀y}dy;
* **V_R = ½erfc((1 − p_e)/σ)·∫₀^T|ψ₂(1 + iy)|e^{y²/σ² − κ₀y}dy** (the edge, suppressed by the cutoff);
* Top_χ = e^{−κ₀T}∫_{p_w}^1|ψ₂(p + iT)|w_T(p)M_top(p)dp, w_T = 1 + ½e^{T²/σ²}erfc((p_e − p)/σ) for p < p_e, ½e^{T²/σ²}erfc((p − p_e)/σ) for p ≥ p_e.

**Proof [P].**
1. Linearity: ∫_{p_w}^1ψ₂G = ∫_{p_w}^1ψ₂χG + ∫_{p_w}^1ψ₂(1−χ)G. The first term is (EDGE-R).
2. 1 − χ(z) = ½erfc((z − p_e)/σ) is entire, so ψ₂(1−χ)G_± is analytic on Q_±. Cauchy moves ∫_{p_w}^1 to the left, top and right sides.
3. For all real a, b: erfc(a + ib) = (2/√π)∫₀^∞e^{−(a+ib+s)²}ds, hence **|erfc(a + ib)| ≤ e^{b²}erfc(a)**. This gives the weights: on the right
   side and for p ≥ p_e directly; for p < p_e via 1 − χ = 1 − ½erfc((p_e − p − iy)/σ).
4. Two waves (×2), ring weight (×2); the monopole term is separate and exact (PW4). ∎

**Remarks.**
* (AG′) ∧ (EDGE) imply (AG_U) ∧ (EDGE-R) (the latter by the same Cauchy argument applied to ψ₂χG), so the new premise pair is **weaker**.
* (EDGE-R) is a **single real-axis integral of the exact correction against a fixed smooth weight**. It needs no analytic continuation and
  no rate, and it is directly computable from real-axis data of the graviton correction near the edge, if such data are ever produced.
* The local rate of (AG′) entered only to suppress the edge amplitude by κ_loc⁻⁶. Here the same job is done by erfc((1 − p_e)/σ) on the
  Cauchy side, while the true edge contribution is carried by (EDGE-R), where its cancellations are kept (§4).

### 5.2 Theorem W's 47 states (`ag_chi.py`, `cache/chi.txt`) [DV]
Per state, (p_e, σ = sT, T = tT₀) is chosen on a grid (p_e ∈ {0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75}, t ∈ {1, 0.8, 0.65, 0.5}, s ∈ {0.35,
0.45, 0.6}, subject to κ(p_e)σ ≥ 10) to minimise the total at Λ_AG = 1. d_χ = window + boundary + right side + top at Λ_AG = 1 (i.e.
everything except E_χ); ε_E tol = H − d_χ, in units of S_E = (2·2V_flat,2 + u·2V_flat,3)·e_sc(1) (plain) or its +−-transferred version.

{{CHITABLE}}

### 5.3 The exact scalar's (EDGE-R) quantity with these cutoffs (`ag_edge2.py`, `cache/edge2.txt`) [DV]
{{EDGE2TABLE}}

---

## 6. Weaker or re-shaped decay premises in Lemma P4W itself (`ag_p4w.py`, `cache/p4w.txt`) [DV]

Families replacing (AG′)'s edge and top majorants (boundary side and window unchanged; rectangle height T = 25/(θΛβ) so that the top
stays e^{−25}-suppressed):
* **Lt(θ):** AG′'s local rate θΛβ/√(1−x₁²) on the whole side, no threshold;
* **Fa(θ):** |flat₊(1+iy)/flat₊(1)|^θ (Lemma F's exact profile, amplitude included);
* **U(θ):** the uniform rate θΛβ.

### 6.1 Theorem W, 47 states
{{FAMTABLE}}

### 6.2 Theorem L′, 16 states
{{LTABLE}}

---

## 7. Why (AG′) cannot be proved by the suggested routes

### 7.1 The plane-wave / spin-Gamma representation (LARGESPIN_L L2, Corollary B)
**Lemma R [P].** Let f(ν) = ∫dμ(ω)f_ω(ν) with μ ≥ 0 independent of ν, and suppose each f_ω = f_{ω,+} + f_{ω,−} with f_{ω,±} analytic on
Q_±, jointly measurable, |f_{ω,±}(p ± iy)| ≤ M_ω(p)e^{−φ_ω(p,y)} and ∫dμ M_ω < ∞. Then f_± := ∫dμ f_{ω,±} are analytic on Q_± (Morera +
Fubini) and |f_±(p ± iy)| ≤ ∫dμ(ω) M_ω(p) e^{−φ_ω(p,y)}. ∎

So the AG class is closed under the positive smearing that produces every heavy image (J enters only through a Gamma(J+1) average). Lemma R
reduces (AG′) for all states to a CFT-data-free statement about the plane-wave images 𝒜_ν(p), **with the rate of the worst part of the
support and with no cancellation between terms**. Two quantitative obstructions:

1. **Beyond-edge mass (`ag_pwmeasure.py`, Monte Carlo of the exact measure p = sn + q, s ~ Gamma(J+1), q ~ Riesz(λ = (τ−1)/2)):**

   | J (gap) | 1 − x₁ | sd\|p\|/m | P(\|p\| < Λ) | P(\|p\| < 1.02Λ) | moment-weighted (w = \|p\|⁻¹⁴): P(\|p\| < Λ) | P(\|p\| < 1.02Λ) |
   |---|---|---|---|---|---|---|
   | 26 | 0.025 | 0.031 | 0.21 | 0.43 | **0.35** | 0.60 |
   | 30 | 0.029 | 0.031 | 0.18 | 0.38 | **0.31** | 0.55 |
   | 38 | 0.036 | 0.031 | 0.12 | 0.29 | **0.23** | 0.46 |
   | 42 | 0.040 | 0.031 | 0.10 | 0.25 | 0.20 | 0.41 |
   | 50 | 0.047 | 0.031 | 0.06 | 0.18 | 0.14 | 0.32 |
   | 64 | 0.059 | 0.031 | 0.03 | 0.10 | 0.07 | 0.19 |
   | 100 | 0.087 | 0.030 | 0.002 | 0.010 | 0.006 | 0.03 |
   | ≥ 150 | ≥ 0.12 | 0.030 | 0 | 0 | 0 | ≤ 0.001 |

   The |p|-width (3 %, set by τ: 1/√(2λ)) exceeds the distance to the twist edge (2.5–4.7 %) for every Theorem W state with m ≲ 1.05Λ. A
   plane wave with |p| < Λ reaches its own edge ν = |p| *inside* [p_w, 1]; beyond it the flat-type image is non-oscillating and grows like
   e^{J·acosh}. Lemma R's bound ∫dμ|f_{ω,+}| is then not small at all. **The exact image's finite, oscillating edge behaviour is produced by
   cancellation between plane waves**, the same mechanism as N-L1 (smearing 1–50 % vs correction 10⁻³). No positivity argument reaches
   x₁ > 0.94. For J ≥ 100 (x₁ ≤ 0.91) the obstruction is absent (≤ 0.6 %), and the route reduces (AG_U)/(AG′) there to the open
   statement PW-AG for 𝒜_ν(p).
2. **Rate dilution.** Even inside its edge, the measure's rapidity spread (sd β(p) = 0.048–0.059, LARGESPIN_L §2.3) is comparable with β
   itself for J ≤ 64 (β = 0.05–0.12). The worst-support rate in Lemma R is then ≪ β.

### 7.2 The exact Laplace/Barnes series
The exact scalar moment is an explicit entire function (§3.1), so analyticity is free. What is not available is the **split**. The
oscillation in ν arises from a ≤ 10⁹-fold cancellation in Σ_m g_m(ν)R_m (EXACT_SKELETON §1), where every g_m is a non-oscillating
polynomial in ν² and R_m carries the Mack polynomial (degree J in the Mellin variable). There is no term-wise wave decomposition, and the
natural candidate (flat-locked waves with a smooth ratio) does not fit (§3.3). A split would need a Mellin-space analogue of
P_J = W₊ + W₋ for the Mack polynomial — not known.

### 7.3 Steepest descent with rigorous remainder (HLJ's crystal)
The ν-dependence is exactly e^{(iν/2)log A(z)} on a fixed contour (HLJ §1.1 [P]), and every saddle order is AG [D]. A beyond-all-orders
statement needs a rigorous thimble decomposition of a ≈ 10-dimensional integral with a soft fibre direction and a fold at ν ≈ 11 (HLJ
§1.3, HWIN_THIMBLE). That is beyond any explicit-constant method available here. More importantly for the edge: the saddle analysis was
validated (lj_scan) only at x ≤ 0.86, and N1 shows that the true near-edge profile is controlled by the singularity at ν = m, i.e. by the
exact function's behaviour near its own twist edge. That region is where the representation's terms spread across the edge (§7.1), and
where §3.3 shows the exact correction is not a simple deformation of the flat wave.

### 7.4 What would close (AG_U) and (EDGE-R)
* **(EDGE-R) for the graviton:** real-axis values of the exact graviton correction G_e on p ∈ [p_e − 3σ, 1] (p_e ≈ 0.2–0.75) at the near-gap
  states (J ≤ 51, m ≤ 1.05Λ). One integral per entry decides it. No complex continuation is needed.
* **(AG_U):** the uniform rate only. It is a statement about the bulk saddle structure, with no edge content. It is [D] to all orders
  (HLJ: ≥ 1.002β > 0.97β). Beyond all orders it is the irreducible form of H_AG.

---

## 8. Updated theorem statements

**Theorem W (window form, χ-version).** Let O be a heavy primary in H_LJ(1000) with couplings λ ∈ ℂ² (even J) or ℂ¹ (odd J). Assume:
* **(P0)**, **(P1′)**, **(W)** exactly as in P4_WINDOW §7;
* **(AG_U)** for every entry: G_e = G_{e,+} + G_{e,−} on Q_± = [p_w, 1] × [0, ±T] (T ≤ T₀ = 25/(0.97Λβ)), analytic, with the **uniform**
  rate 0.97Λβ and amplitudes M_w (left side, the window series' waves), Λ_AG·e_sc(1) + A (right side), Λ_AG·e_sc(p) (top);
* **(EDGE-R)** for every entry: |cost of E_χ[G_e]| ≤ ε_E tol (table §5.2), with the per-state cutoff of §5.2.

Then M_Ω(O) ⪰ 0 for every coupling vector, at every sampled state with J ≥ 26 (m ≤ 1.3Λ), plain envelope Λ_AG = 1.
{{WSTATEMENT}}

**Grade.** [P given P0, P1′, W, AG_U, EDGE-R] (Lemma P4W-χ + the 2×2 criterion) + [DV] (sampled constants). Between sampled states: [NE].

**Changes relative to P4_WINDOW §7.**
* (AG′) → (AG_U): the local rate and its threshold are gone; nothing is assumed about the near-gap edge in complex directions.
* (EDGE) → (EDGE-R): an amplitude bound at one point becomes one real-axis integral, which keeps the cancellations. The exact scalar
  satisfies it with margin ≥ {{MARGIN}}.
* The old statement (with (AG′) ∧ (EDGE)) remains valid as a corollary, but its premise (AG′) is now known to be pointwise false for the
  flat wave at x₁ ≥ 0.95 (N1), so it should not be quoted in that form.

**Theorem L′ (unchanged threshold J ≥ ⌈(0.37211m − 1)/2⌉).** {{LSTATEMENT}}

---

## 9. Honest residual
1. **(AG_U) beyond all orders is assumed [A].** It is now the only analytic hypothesis, and it is uniform (no edge content). To all orders
   [D] (HLJ), with margin 1.002/0.97. Every explicit piece satisfies it [P] (Lemma F, Harish-Chandra, Hankel).
2. **(EDGE-R) for the graviton is [M].** For the scalar analogue it is [DV] with the margins of §5.3, at x₁ = 0.91–0.96 only (no exact
   scalar exists at J = 26–38, m ≈ 1.03Λ).
3. The +− helicity transfer (Ω^U convention) is still a convention [M]. {{RESID_TR}}
4. (W), the window O(1/m⁴) coefficients, the boundary amplitude factor 2 [D-sketch], P0's [NE] part, sampled states only, float
   quadrature: all unchanged from P4_WINDOW §8.
5. The (EDGE-R) quantity depends on the cutoff (p_e, σ). The premise is stated per state for the tabulated cutoff. The exact scalar's value
   is stable to ±5 % across cutoffs (§4).
6. §3.3's negative (no flat-locked split) means that **no numerical evidence for (AG′)'s pointwise form exists at the near-gap edge**. The
   evidence is for its integrated consequence (§4) and for the weaker pair (AG_U, EDGE-R).

## 10. Suggested inventory changes (not applied)
* **P4_WINDOW §6 / §7, LARGESPIN_W, UNIFIED_THEOREM H_LJ(1000):** "[P given P0, P1′, W, AG′, EDGE]" → "[P given P0, P1′, W, AG_U, EDGE-R]
  (AG_PROOF Lemma P4W-χ)". Add: "(AG′)'s pointwise local-rate majorant is violated by the flat wave itself at x₁ ≥ 0.95 (AG_PROOF N1);
  only its integrated form is used, and it is now replaced."
* **P4_WINDOW §4.2 'sensitivity to the rate threshold':** an artefact of the piecewise majorant (N1); the exact flat profile (family Fa)
  has no threshold.
* **HLJ_ANALYTIC Lemma C / P4_WINDOW (2.1):** add Lemma P4W-χ as the edge-free variant.
* **LARGESPIN_L §2.4 (R1 reduction):** add N2: at J ≤ 50 near the gap, 14–35 % of the moment-weighted plane-wave measure lies beyond the
  twist edge; the representation cannot transfer near-edge structure.

## 11. Reproduce

    cd ads4/code/agproof && bash run_all.sh        # ≈ 1.6 CPU-h, ≤ 2 processes

| file | content | grade |
|---|---|---|
| `cache/flat.txt` | Lemma F checks; flat profile vs AG′ at 67 states (§1–2) | [P] + [DV] |
| `cache/unif.txt` | uniform-rate and Abel amplitude bounds for the flat waves (§1.1) | [DV] |
| `cache/exact_lines_<J>.pkl` | exact scalar on vertical lines (§3.2) | [DV] |
| `cache/lock.txt` | flat-locked fit (negative, §3.3) | [NE] |
| `cache/edge.txt`, `cache/edge2.txt` | edge contribution on the real axis (§4, §5.3) | [DV] |
| `cache/chi.txt` | Lemma P4W-χ on the 47 W states (§5.2) | [DV] |
| `cache/p4w.txt` | rate families, W and L′ (§6) | [DV] |
| `cache/pwmeasure.txt` | beyond-edge plane-wave mass (§7.1) | [DV] (MC) |
