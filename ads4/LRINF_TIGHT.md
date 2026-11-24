# LRINF_TIGHT: Lemma LR∞ certified with the tight constant ∫|ψ| — new threshold β₀, new J\* for Theorem L, and coverage of the near-gap band

*23 Nov 2026. Unreviewed. No existing file was edited and nothing was committed. Follows up INDEP_CHECK2 target 2 (D1).*
* *New code:* `ads4/code/lrtight/` (`run_all.sh`, ≈ 10 CPU-min, one process at a time). Outputs go to `code/lrtight/cache/`, whose `.gitignore`
  ignores everything except itself.
  * `lrt_common.py`: exact ψ_k, and the bound functions in Arb balls.
  * `lrt1_constants.py`: certified ∫|ψ_k|.
  * `lrt2_m1.py`: certified lower bounds on M₁.
  * `lrt3_rcert.py`: certified r(β), r̃(β), F_cert and the β ≥ 64 tail.
  * `lrt4_jstar.py`: the J\*(m) tables.
  * `lrt5_tolerances.py`: Theorem L′ tolerances, sanity checks, hybrid vs F_LR.
* *Inputs.*
  * The certified part (lrt1–lrt4) reads only `prl/data/nf02_psi.json`, and uses python-flint 0.9 (Arb, 128 bits) and sympy (exact Sturm
    sequences).
  * The [DV] part (lrt5) imports, read-only, `pwimg/pw5_theoremL.py` (and through it `lsl/lsl_core.py`, `n2sub/n2_common.py`,
    `pwimg/pw4_theoremW.py`) and `largespin/ls_core.py`. No file of theirs was written, and no bytecode was written.
* *Method sources.*
  * The M₁ lower bounds re-implement `exactskel/h13_cert_m1.py` region (b).
  * The bound formula is LARGESPIN_L §3 as stated, with `R = min(R_∞, a₁/x² + ∫|ψ|)`.

Grades: [P] proved; [certified] = Arb ball arithmetic together with proved monotonicity / rigorous Riemann sums; [D] derived; [DV] float numerics;
[NE] numerical evidence (sampled); [M] model; [I] open. Notation as in LARGESPIN_L: `sinh(β/2) = (J+½)/m`, `b̂ = 2Λ sinh(β/2)`, `u = Λ²/m²`,
`Λ = √(998·999) = 998.49987`, `F_LR = 2Ω₂ + uσΩ₃`.

---

> **Update (GAP_LARGESPIN.md, unreviewed): the large-spin floor is now certified for EVERY Δ_g ≥ 1000 (no upper limit), and the large-spin theorems' margins grow with the gap.**
> * **Lemma LR∞ for all Δ_g ≥ 1000 [P + certified]:** a scaling lemma makes Λ⁴·(each bound piece) non-increasing in Λ (envelope R̂ = min(R_∞, V/s⁴), V₂ ≤ 684625, V₃ ≤ 789636), plus a 13-step monotone Λ-grid below 2604.28. Thus r < 1 on β ≥ 0.37 and F_LR > 0 on β ≥ 0.245; J\*(m) = ⌈(0.37211m−1)/2⌉, with J\*/Λ → 0.22387 as Δ_g → ∞. This closes the tail above Λ = 10⁶.
>   * Thresholds improve with the gap: β₀ = 0.2375 / 0.1675 / 0.1175 for Δ_g ≥ 2606 / 5000 / 10⁴.
>   * Gap-unit floor: Λ²F_LR ≥ f_U(β) for all Δ_g ≥ 1000 [certified].
> * **ROUTEA curvature window uniform in the gap [D + certified]:** curv > 0 on b̂ ∈ [114.06, 2Λ] for m ≥ Λ (new exact identity ∫₀^{s₁}g₁ = −s₁³r(s₁)/24). **Lemma C′ is monotone in Λ [P]**, so b̂ ≤ 119 carries to every gap. **For m ≥ 1.3Λ the hybrid floor is rigorous on every b̂ ≤ 2Λ at every gap**, which closes the old (119, 245) hole for m in [1.3, 1.5)Λ.
> * **Theorems L′ / W at Δ_g = 1000 / 2000 / 5000 [DV, premises carried with Δ_g-independent AdS-unit constants]:**
>   * Λ²D = −995.3 / −999.1 / −1001.3.
>   * L′ remainder: max d/(F/2) = 0.19 / 0.035 / 0.005; edge tolerance 12.6 / 228 / 9050.
>   * **Theorem W: d_χ/H = 0.489 / 0.246 / 0.100; the tolerated graviton edge integral at the binding state is 0.406 / 3.45 / 57 S_E,t (∝ Λ^3.1).** So the EDGE-R premise gets much easier with the gap.
>   * Caveat: the 5000 C2/vector tolerance entries are not numerically converged (finite-difference β-derivative).
> * **Remaining:** the premises (P_L′ magnitudes, AG_U, EDGE-R, P1′, W) carry over only under the hypothesis of Δ_g-independent AdS-unit constants. The trend between computed gaps is [D-sketch]. The small-J twist-edge seam and β₀ < 0.10 are not addressed.

## 0. Summary

| quantity | old (LARGESPIN_L / PW_IMAGES / INDEP_CHECK2) | new (this note) | grade |
|---|---|---|---|
| small-x constant in R | code: Φ_max = Σ\|c_n\|/(n+1) = 7.92·10⁵ (k=2), 1.88·10⁶ (k=3); text: ∫\|ψ\| (uncertified, INDEP2: 10.39) | **∫\|ψ₂\| = 10.386029907129, ∫\|ψ₃\| = 12.835622283960** (balls ±7·10⁻²⁰) | **[certified]** (exact Sturm count + exact bisection + Arb) |
| r(β = 1) (sup over β ≥ 1) | 0.577 | **0.0162** (INDEP2: 0.0162, float) | [certified] |
| sup r over β ≥ 0.5 | — (INDEP2: 0.28, float, uncertified) | **0.280** | [certified] |
| smallest β₀ with r < 1 on [β₀, ∞) | 1 | **β₀ = 0.37** (b̂₀ = 0.37211Λ = 371.56) | [certified] (cells of width 0.005; the cell [0.365, 0.37] has r ≤ 1.056) |
| smallest β₀ with **F_LR > 0** (r̃ < 1, flat tail kept) | — | **β₀ = 0.245** (b̂₀ = 0.24561Λ = 245.24) | [certified] (the cell [0.2425, 0.245] has r̃ ≤ 1.003) |
| β > 64 tail | [D] (analytic, not evaluated) | **r ≤ 2.03·10⁻⁷ for all β ≥ 64** (explicit majorant built from decreasing terms) | [P + certified] |
| M₁ input | M₁ > 0 [certified] (h13), float values on the grid | certified **lower bounds** on 910 β-cells covering [0.1, 24], plus the region-(c) bound (re-derived, [P]) beyond | [certified] / [P] |
| Hadamard step of the lemma (LARGESPIN_L §3, step 3) | [D + DV] | **[P]**: three-line proof by s = sinh(t/2) (§1.1) | [P] |
| Lemma LR∞ overall | [P + DV] on β ≥ 1; tail [D] | **[P + certified] on every β ≥ 0.37 in the stated form `F_LR ≥ (1−r)A₁M₁/Λ² + A₂′/(2b̂³)`; F_LR > 0 [P + certified] on β ≥ 0.245** | [P + certified] |
| J\*(m), stated form | ⌈(1.0421 m − 1)/2⌉ | **⌈(0.37211 m − 1)/2⌉** (r < 1); ⌈(0.24561 m − 1)/2⌉ (F_LR > 0 only) | [certified] (integer decisions in Arb) |
| J\* at the twist edge | 857 | **223** (r < 1); 139 (F_LR > 0, small certified margin) | [certified] |
| J\* at 1.1 / 1.3 m_gap, 1.5Λ, 2Λ, 5Λ | 987 / 1275 / 780 / 1041 / 2602 | **250 / 307 / 279 / 372 / 929** | [certified] |
| certified floor F_cert / numerical min F_LR at J\* | (floor not certified) | 0.57 (β = 0.37), 0.82 (0.5), 0.98 (1.0) | [certified] / [DV] |
| monopole tolerance Q_tol at J\* (gap) | 3.6·10⁻³ (J = 857, numerical floor) | **4.8·10⁻³ (certified floor) / 8.5·10⁻³ (numerical), J = 223** | [DV] |
| Λ_amp tolerated, monopole only, at J\* (gap) | 1003 (J = 857) | **1025 (certified floor) / 1819 (numerical)**; ≥ 901 on 223 ≤ J ≤ 857 (cert., sampled) | [DV] |
| Λ_amp tolerated, full scalar-shaped vector R1′, at J\* (gap) | 1240 (J = 857) | **181 (certified floor) / 322 (numerical)**; rises to 1114 / 1240 at J = 857; 32 / 43 at J = 10⁶ | [DV] |
| near-gap band 345 ≤ J < 857 | no theorem; floor [NE] | **inside Theorem L′** (J ≥ 223 at the gap, every J, not only sampled J); floor (leading-Regge part) [P + certified] | Theorem L′ [D given P_L′]; tolerances [DV] |
| Theorem W (J\* = 64 with k = 3, sampled to 344) vs Theorem L at the gap | separated by 345–857 | **overlap on 223 ≤ J ≤ 344** (W's sampled gap states J = 250, 300, 301, 344 are inside L's range) | — |
| floor coverage, m ∈ [1.3, 1.5)Λ | b̂ ∈ (119, 1040.6) [NE] | b̂ ∈ (119, 245.2) [NE] | — |
| floor coverage, m ∈ [Λ, 1.3Λ) | b̂ < 1040.6 [NE] | b̂ < 245.2 [NE] | — |
| hybrid − F_LR in the new band | +1·10⁻⁶ (J ≥ 1000) … +1.7·10⁻³ (J = 100) | positive at every sample; ≤ 8.1·10⁻⁴ relative for J ≥ 139 (gap and 1.3 m_gap, +− worst) | [NE] |

**Bottom line.**
* **Certified β₀ = 0.37 for Lemma LR∞ as stated (r < 1), and 0.245 for F_LR > 0.** This holds for all m ≥ Λ (u ≤ 1) and both helicities.
* **New J\*(m) = ⌈(0.37211 m − 1)/2⌉, which is 223 at the twist edge** (857 before).
* **The 345–857 band is now inside Theorem L′'s range.** Its leading-Regge floor is [P + certified]. The theorem itself remains conditional
  on P_L′ (magnitudes), as for J ≥ 857.
* **Theorems W and L now overlap at the gap on J ∈ [223, 344].**
* This is not an unconditional theorem. The premise status (P_L′ magnitudes [NE]/[M]) is unchanged.

---

## 1. Certification method

### 1.1 The lemma (LARGESPIN_L §3), and a proof of its step 3

`Ω_k(β) = a₁M₁(β)/Λ² − a₂/(2b̂³) + E_k`, `|E_k| ≤ (√2/π)[T_a + T_b + T_c₁ + T_c₂]`. The definitions are those of LARGESPIN_L §3, with
`R(x) = min(R_∞(x), a₁/x² + ∫|ψ|)`.

I re-read the proof. Steps 1, 2, 4 and 5 are elementary and correct. Two points the original left implicit are checked here:
* the ρ_b bound needs ε = (cosh t − 1)/D ≤ ½ on t ≤ β/2. In fact ε ≤ 1/(4cosh²(β/4)) ≤ ¼ for **every** β > 0;
* `0 ≤ ρ ≤ K` on [β/2, β].

So the inequality holds for every β > 0, not only β ≥ 1. Nothing in the proof used β ≥ 1; the threshold came only from evaluating r.

**Step 3 [P] (previously [D + DV]).**
* Write `c(t) := 2∫₀^∞ν cos(νt)/(e^{2πν}+1)dν = cosh(t/2)/(4sinh²(t/2)) − 1/t²`. This follows from
  `∫₀^∞sin(at)/(e^{2πx}+1)dx = 1/(2a) − 1/(4sinh(a/2))`, the difference of the two Bose-type integrals.
* Mehler–Dirichlet and Fubini (absolute) give `M₁ = (√2/π)∫₀^βK c dt`.
* Split `c = −1/t² + h`, with `h = cosh(t/2)/(4sinh²(t/2))`. The claim `(√2/π)·[−f.p.∫₀^βK/t²] = M₁` is then equivalent to
  `f.p.∫₀^βK h dt = 0`.
* Substitute `s = sinh(t/2)`, `S = sinh(β/2)`. Then `K h dt = ds/(2s²√(2(S²−s²)))`, whose antiderivative is `−√(S²−s²)/(2√2S²s)`.
* This vanishes at s = S. At s → 0 it equals `−K₀/ε + O(ε)` (K₀ = 1/(√2S), s = sinh(ε/2)), so the finite part is 0. ∎

### 1.2 ∫|ψ_k| (`lrt1_constants.py`) [certified]
1. ψ_k = p(1−p)⁵X_k(p)/k₀ with exact rationals (k₀ = 27.48509817 is an exact decimal).
2. sympy's exact Sturm count gives **one root of X_k in [0, 1]** (k = 2: p = 0.91682…; k = 3: p = 0.92328…). It is confirmed by the count of
   ψ_k on [10⁻¹², 1 − 10⁻¹²].
3. Exact rational bisection brackets that root to width < 2⁻¹²⁰.
4. ∫|ψ| = |Q(r) − Q(0)| + |Q(1) − Q(r)|, with Q the exact antiderivative evaluated in Arb on the bracket.
5. Cross-check: rigorous Riemann sums on 4096 cells, with the sup taken from the exact Taylor expansion at each cell centre, bracket the
   value (10.38603 ≤ ∫|ψ₂| ≤ 10.39120).

| k | ∫\|ψ_k\| | ∫ψ_k | crude Φ_max (old code) | ratio |
|---|---|---|---|---|
| 2 | 10.386029907129 | 10.385946507 | 7.921·10⁵ | 7.6·10⁴ |
| 3 | 12.835622283960 | 12.835567817 | 1.876·10⁶ | 1.5·10⁵ |

a₁, a₂ are exact rationals: a₁ = 188.89611 / 50.64297 and a₂ = −1023.5457 / 7.28·10⁻⁵. Hence A₁,min = 2a₁⁽²⁾ − a₁⁽³⁾ = 327.14924 and
A₂′ = −2a₂⁽²⁾ − |a₂⁽³⁾| = 2047.09132 (exact).

### 1.3 Lower bounds on M₁ (`lrt2_m1.py`) [certified], [P]
* **β ∈ [0.1, 24]: 910 cells.** Cell widths: 0.0025 on [0.1, 0.3], 0.005 to 1.5, 0.01 to 3, 0.025 to 6, 0.05 to 12, 0.1 to 24.
  * Method: h13 region (b). Set t = β(1 − v²), take β as an Arb ball covering the whole cell, and form lower Riemann sums in v (1600–25600
    cells) of the analytic integrand.
  * c(t) on t ≤ 2 comes from its series, with an explicit tail bound.
  * Every cell bound is valid for **every** β in the cell, so no interpolation is involved.
  * Loss against float M₁ (min over the cell): 0.3–1.2 % on [0.1, 2], 7 % at β = 3, 11 % at β = 6, 23 % at β = 12 and 31 % at
    β = 24. It is irrelevant where r ≪ 1.
* **β ≥ 12 (used for β > 24): h13 region (c),**
  `I/K₀ ≥ 1/β − 1/(2sinh(β/2)) − √2(e^{−β/2}/4.7 + 4/β²) − 2√(2ln2)/(β−ln2)²`.
  I re-derived it line by line and grade it [P], given h13's certified c > 0 on (0, 4.7]:
  * r ≤ ½ ⇒ (1−r)^{−½} − 1 ≤ √2r;
  * sinh(t/2)/sinh(β/2) ≤ e^{(t−β)/2};
  * t₁ ≥ β − ln 2, because sinh(x − ln2/2) ≤ sinh(x)/√2;
  * sinh t₁ ≥ sinh²(β/2).

  On the cells (24, 64] it is evaluated with β as a ball.

### 1.4 Between-grid control of the remainder (`lrt3_rcert.py`) [P + certified]
Every piece of the bound is monotone in β (proved below). On a cell [β_a, β_b], each piece is therefore bounded by its value at an endpoint,
evaluated in Arb.

| piece | monotonicity | bound used on the cell |
|---|---|---|
| T_a = (K₀/Λ)d₁(Λβ) | K₀ = (cosh β − 1)^{−½} and d₁ (positive coefficients × X^{−n}) both decrease | value at β_a |
| T_b = K₀ sumc(Λβ)/(4DΛ³) | D increases, sumc decreases | value at β_a |
| T_c₁ = ∫₀^{β/2}R(Λt)ρ_b(t; β)dt | ρ_b = K₀[t⁴cosh t/(48D) + (3/8)2^{5/2}((cosh t−1)/D)²] decreases in β at fixed t; the range [0, β/2] grows | `K₀(β_a)[I₁(β_b/2)/(48D_a) + (3/8)2^{5/2}I₂(β_b/2)/D_a²]` |
| T_c₂ = R_∞(Λβ/2)(π/√2)P_{−½}(cosh β) | R_∞ decreases; P_{−½}(cosh β) = ₂F₁(½,½;1;−sinh²(β/2)) decreases (Euler integral) | value at β_a (Arb hypergeometric) |
| flat term A₂′/(2b̂³) | decreasing | value at β_b |
| M₁ | (none needed) | certified cell lower bound (§1.3) |

The two cell-independent integrals are `I₁(T) = ∫₀^T R(Λt)t⁴cosh t dt` and `I₂(T) = ∫₀^T R(Λt)(cosh t − 1)²dt`. They are bounded by
**upper Riemann sums** on one geometric t-grid (t₁ = 10⁻⁶, ratio 1.0005, 34 604 points, up to t = 32.5):
* R(Λt) decreases (a min of decreasing functions), so it is taken at the left end of each step;
* t⁴cosh t and (cosh t − 1)² increase, so they are taken at the right end;
* on [0, t₁], R ≤ a₁/(Λt)² + ∫|ψ|, t⁴cosh t ≤ t⁴cosh t₁ and (cosh t − 1)² ≤ t⁴cosh²t₁/4, integrated exactly.

The cumulative sums are stored. For T = β_b/2 the first grid point ≥ T is used, which is valid because the integrand is ≥ 0.

### 1.5 Tail β ≥ 64 [P + certified]
For β ≥ 64:
* M₁ ≥ (√2/π)K₀·ℓ/β with ℓ = 64·L(64) = 0.874007. Here βL(β) is increasing, since each subtracted term decreases for β ≥ 64.
* r(β) ≤ Λ²/(A₁,min ℓ) · Σ_k w_k[f_a + f_b + f_c₁ + f_c₂ + f_d](β), with w₂ = 2 and w₃ = 1.

Each f is an explicit decreasing function:
* f_a = βd₁(Λβ)/Λ;
* f_b = β sumc(Λβ)/(4DΛ³);
* f_c₁ = βI₁(32)/(48D) + ρ₀β/(96 sinh(β/2)), with ρ₀ = R_∞(32Λ)·32⁴. This uses R_∞(Λt)t⁴ ≤ ρ₀ for t ≥ 32;
* f_c₂ = (3/8)2^{5/2}[βI₂(32)/D² + R_∞(32Λ)(β²/(16sinh⁴(β/2)) + β cosh(β/2)/(8 sinh³(β/2)))];
* f_d = βR_∞(Λβ/2)(C₁ + C₂β)/√2, with C₁ = 2√2e^{1/4} and C₂ = (2/(1−e^{−1}))^{½}. This comes from
  `∫₀^βK ≤ e^{−β/4}(C₁ + C₂β)/√(2sinh(β/2))`, using sinh((β+t)/2) ≥ sinh(β/2)e^{t/2}.

Evaluated at 64 in Arb: **r(β) ≤ 2.03·10⁻⁷ for every β ≥ 64**. LARGESPIN_L's [D] tail is now closed.

### 1.6 The positivity corollary
For u ≤ 1 and σ = ±1 (all m ≥ Λ, both helicities), with M₁ > 0:

    F_LR ≥ A₁,min M₁/Λ² + A₂′/(2b̂³) − (2|E₂| + |E₃|) =: F_cert,
    r  := (2|E₂| + |E₃|)/(A₁,min M₁/Λ²)                      (LARGESPIN_L's r; F_LR ≥ (1−r)A₁M₁/Λ² + A₂′/(2b̂³)),
    r̃ := (2|E₂| + |E₃|)/(A₁,min M₁/Λ² + A₂′/(2b̂³))           (F_LR > 0 ⇐ r̃ < 1).

The cell bounds give certified upper bounds on r and r̃ and a certified lower bound F_cert, each valid on the whole cell.

### 1.7 Checks
* **Reproduction of INDEP_CHECK2's float values.** r = 0.280, 0.131, 0.0694, 0.0401, 0.0249, 0.0200, 0.0162 at β = 0.5 … 1, against
  INDEP2's 0.28, 0.13, 0.069, 0.040, 0.025, 0.020, 0.0162. Every value agrees.
* **[DV] sanity check at β = 0.25 … 3** (`lrt5`, direct Mehler quadrature with closed-form Φ):
  * the true |E_k| (2·10⁻¹¹ … 4·10⁻⁹) lies below the certified bound everywhere;
  * F_cert ≤ the numerical min F_LR everywhere;
  * F_cert/F_num = 0.035 (β = 0.25), 0.33 (0.30), 0.57 (0.37), 0.82 (0.50), 0.98 (1.0), 0.97 (3.0).
* **J\* reproduction.** At β₀ = 1 the table gives 857, 987, 1275, 780, 1041, 2602 (LARGESPIN_L §3.3); at β₀ = 0.5 it gives INDEP2's 324.

---

## 2. r(β) table [certified]

The sup is taken over all cells with β_a ≥ β′, together with the β ≥ 64 tail. "Cell at β′" is the cell [β′, β′ + h].

| β′ | sup_{β ≥ β′} r | sup r̃ | cell at β′: r | r̃ | F_cert (cell) | k = 2 split T_a / T_b / T_c₁ / T_c₂ (× √2/π) |
|---|---|---|---|---|---|---|
| 0.245 | 5.88 | **0.984** | 5.88 | 0.984 | 1.30·10⁻⁶ | 1.3e-7 / 2.0e-7 / 7.6e-6 / 1.1e-5 |
| 0.30 | 2.45 | 0.668 | 2.45 | 0.668 | 1.64·10⁻⁵ | 5.9e-8 / 8.8e-8 / 2.8e-6 / 4.7e-6 |
| 0.345 | 1.34 | 0.486 | 1.34 | 0.486 | 1.90·10⁻⁵ | — |
| 0.365 | 1.056 | 0.424 | 1.056 | 0.424 | 1.92·10⁻⁵ | — |
| **0.37** | **0.997** | 0.410 | 0.997 | 0.410 | 1.92·10⁻⁵ | 2.6e-8 / 3.8e-8 / 9.8e-7 / 2.0e-6 |
| 0.40 | 0.716 | 0.335 | 0.716 | 0.335 | 1.90·10⁻⁵ | — |
| 0.44 | 0.479 | 0.258 | 0.479 | 0.258 | 1.83·10⁻⁵ | — |
| 0.50 | 0.280 | 0.177 | 0.280 | 0.177 | 1.73·10⁻⁵ | 7.6e-9 / 1.1e-8 / 2.2e-7 / 5.8e-7 |
| 0.60 | 0.131 | 0.098 | | | 1.58·10⁻⁵ | |
| 0.70 | 0.0694 | 0.0572 | | | 1.48·10⁻⁵ | |
| 0.80 | 0.0401 | 0.0351 | | | 1.40·10⁻⁵ | |
| 0.90 | 0.0248 | 0.0226 | | | 1.35·10⁻⁵ | |
| 1.0 | **0.0162** (stated 0.577) | 0.0151 | | | 1.30·10⁻⁵ | 4.6e-10 / 6.4e-10 / 6.6e-9 / 3.4e-8 |
| 1.2 | 7.8·10⁻³ (0.225) | 7.5·10⁻³ | | | | |
| 1.5 | 3.2·10⁻³ (0.070) | 3.2·10⁻³ | | | | |
| 2.0 | 1.07·10⁻³ (0.015) | 1.07·10⁻³ | | | | |
| 3.0 | 2.6·10⁻⁴ (1.6·10⁻³) | 2.6·10⁻⁴ | | | | |
| 5.0 | 5.0·10⁻⁵ (1.0·10⁻⁴) | 5.0·10⁻⁵ | | | | |
| 24 | 2.1·10⁻⁶ | 2.1·10⁻⁶ | | | | |
| ≥ 64 | 2.03·10⁻⁷ (tail, all β) | — | | | | |

Reading:
* **Certified β₀ = 0.37 for r < 1** (0.44 for r < ½). **Certified β₀ = 0.245 for r̃ < 1, i.e. F_LR > 0.** Both hold for all m ≥ Λ and
  both helicities.
* At small β the bound is dominated by T_c₂ and T_c₁. R_∞ is already sharp there: R_∞(x) is within 1–26 % of the true envelope
  sup_{y≥x}|r(y)| for x ≥ 75. The non-oscillating p = 0 terms dominate it (k = 2, x = 100: 2.3·10⁻⁵ against 3·10⁻⁶ for the oscillating
  p = 1 terms).
* So a lower threshold would need a structurally sharper lemma, not sharper constants. Two options, neither done [I]:
  * treat the ψ‴(0)/x⁴ term exactly, like the a₁ term (a three-term form);
  * a sharper ρ_b.
* The true r is ≤ 10⁻⁶ (INDEP2), so the bound, not F_LR, is what stops at 0.245.

---

## 3. New J\*(m) [certified]

`J*(m) = ⌈(s₀m − 1)/2⌉`, `s₀ = 2sinh(β₀/2)`, for fixed m.

At the twist edge, J\* is the first J with `(2J+1)² ≥ s₀²c²·998·(999+2J)`. The left side divided by (999 + 2J) is increasing in J, so every
J ≥ J\* also qualifies. The integer decisions are made in Arb (`lrt4_jstar.py`).

| threshold | β₀ | b̂₀/Λ (= s₀) | m_gap | 1.1 m_gap | 1.3 m_gap | 1.5Λ | 2Λ | 5Λ |
|---|---|---|---|---|---|---|---|---|
| LARGESPIN_L (β₀ = 1), reproduced | 1 | 1.04219 | 857 | 987 | 1275 | 780 | 1041 | 2602 |
| INDEP2 side result (β₀ = 0.5), reproduced | 0.5 | 0.50522 | 324 | 365 | 453 | 378 | 504 | 1261 |
| r < ½ (certified) | 0.44 | 0.44356 | 276 | 310 | 382 | 332 | 443 | 1107 |
| **r < 1 (certified; Lemma LR∞ as stated)** | **0.37** | **0.37211** | **223** | **250** | **307** | **279** | **372** | **929** |
| r̃ < ½ (certified) | 0.345 | 0.34671 | 206 | 230 | 281 | 260 | 346 | 865 |
| r̃ < 1 (certified; F_LR > 0 only) | 0.245 | 0.24561 | 139 | 154 | 187 | 184 | 245 | 613 |

* **Operative J\* for Theorem L: β₀ = 0.37, J\*(m) = ⌈(0.37211m − 1)/2⌉ (≈ 0.186m), 223 at the twist edge.**
  * At β₀ = 0.37 the certified floor is ≥ 0.57 of the numerical one, so the premise tolerances lose less than a factor of 2.
  * The r̃ < 1 row (J\* = 139) proves F_LR > 0, but near β = 0.25 the certified floor is only 3–30 % of the true one, and the tolerances
    shrink accordingly (§5). It is reported as the limit of the present lemma, not as the theorem's J\*.
* m_gap/Λ at the new twist-edge J\* = 223 is 1.203. Near-gap states with J ∈ [223, 857] have m ∈ [1.20, 1.65]Λ.

---

## 4. Hybrid floor, Lemma C′ and ROUTEA below β = 1

**Which floor is certified.**
* Lemma LR∞ certifies the leading-Regge action **F_LR**.
* Lemma C′ (b̂ ≤ 119, m/Λ ∈ [1.3, 20], [P + DV]) and ROUTEA (b̂ ∈ [91.2, 2Λ], m ≥ 1.5Λ, [certified]) certify the **hybrid** floor (flat
  NF02 partial wave + curvature).
* Their difference, hybrid − F_LR, is a flat-space partial-wave-vs-Regge term. It is not bounded rigorously anywhere; it is [NE] here, as
  in LARGESPIN_L L1′.

**hybrid − F_LR in the new band** (`lrt5` (iii), `ls_core.hyb_fast`, relative to F_LR) [NE]:

| J | gap ++ | gap +− | 1.3 m_gap ++ | 1.3 m_gap +− |
|---|---|---|---|---|
| 139 | +8.6·10⁻⁵ | +7.4·10⁻⁴ | +9.1·10⁻⁵ | +8.1·10⁻⁴ |
| 223 | +2.2·10⁻⁵ | +2.0·10⁻⁴ | +3.4·10⁻⁵ | +2.6·10⁻⁴ |
| 344 | +5.5·10⁻⁶ | +4.7·10⁻⁵ | +8.3·10⁻⁶ | +7.2·10⁻⁵ |
| 500 | +1.6·10⁻⁶ | +1.2·10⁻⁵ | +2.7·10⁻⁶ | +2.1·10⁻⁵ |
| 700 | +6.4·10⁻⁷ | +3.4·10⁻⁶ | +9.8·10⁻⁷ | +6.2·10⁻⁶ |

The difference is positive at every sample (hybrid ≥ F_LR) and ≤ 8.1·10⁻⁴ relative, far below the certified margin F_cert. So the
"rigorous floor" statements below mean F_LR on b̂ ≥ b̂₀ and the hybrid on the C′/ROUTEA windows. The identification of the two is [NE],
exactly as before.

**Coverage of the leading-Regge / hybrid floor (b̂₀ = 245.2 for F_LR > 0; 371.6 for the r < 1 form):**

| m | before (LARGESPIN_L L1′) | now |
|---|---|---|
| m ≥ 1.5Λ | every b̂: C′ (≤ 119) ∪ ROUTEA [91.2, 2Λ] ∪ LR∞ (≥ 1040.6) | every b̂, unchanged; LR∞ now starts at 245.2, so it overlaps ROUTEA on [245.2, 2Λ] |
| 1.3Λ ≤ m < 1.5Λ | C′ (≤ 119) ∪ LR∞ (≥ 1040.6); **(119, 1040.6) [NE]** | C′ ∪ LR∞ (≥ 245.2); **(119, 245.2) [NE]** (with the r < 1 form: (119, 371.6)) |
| Λ ≤ m < 1.3Λ | LR∞ (≥ 1040.6); below [NE] | LR∞ (≥ 245.2); below [NE] |

* **Is "rigorous at every b̂" still consistent?** Yes. It is still claimed only for m ≥ 1.5Λ.
* **The m ∈ [1.3, 1.5)Λ caveat.**
  * The hole shrinks from b̂ ∈ (119, 1040.6) to (119, 245.2). In spin: J from ≈ (119m/Λ − 1)/2 ≈ 77–89 up to just below J\*(m) = 159–184
    (r̃ form), or 241–279 (r < 1 form).
  * To close it, F_LR > 0 is needed down to b̂ = 119 (β ≈ 0.12). The present lemma reaches only β = 0.245 (r̃(0.2) ≤ 1.43,
    r̃(0.1) ≤ 5.9). It would need the sharper lemma of §2.
* **Which m the new J\* covers.** J\*(m) applies to **every m ≥ Λ** (u ≤ 1). Lemma LR∞ is uniform in u and σ, and A₁,min, A₂′ are the
  worst cases over u ≤ 1. For m < Λ (u > 1) it is not claimed.
* **Near-gap states sit in the weakest m-band.** Near-gap states with 223 ≤ J ≲ 340 have m_gap < 1.3Λ. Below their J\* no rigorous floor
  exists, since neither C′ nor ROUTEA applies at m < 1.3Λ.

---

## 5. Theorem L, updated

**Theorem L (PW_IMAGES form L′, J\* lowered).** Assume A1′, A2 and Δ_g = 1000.
* **Scope:** a heavy primary O in T × T with b̂ ≥ 0.37211Λ, i.e. **J ≥ J\*(m) = ⌈(0.37211m − 1)/2⌉** (m ≥ Λ). That is **J\* = 223 at the
  twist edge**; 250 / 307 at 1.1 / 1.3 m_gap; 279 / 372 / 929 at m = 1.5 / 2 / 5 Λ.
* **Setup:** write `M_Ω = diag(F_LR) + δM`, with the k = 2 part of each entry exactly 𝒜_g[F_ab] (PW2), whose cost is
  `2[∫ψ₂F_ab/c + D·F_ab(i/2)]`.
* **Premise P_L′:** the O(1/m²) magnitude bounds (|κ_s| ≤ Λ_amp|κ_s^sc|, |κ_sh| ≤ K_sh, radial-shift family), plus a remainder fraction
  ε_R.

Then `M_Ω ⪰ 0` for every coupling vector whenever `Λ_amp/Λ_tol + K_sh/K_tol^sh + ε_R ≤ 1`.

**Grade [D given P_L′].** The floor is now **[P + certified]** on the whole range: Lemma LR∞ with certified constants, including the β ≥ 64
tail. The costs and tolerances are [DV]. The proof is that of LARGESPIN_L §5 / PW_IMAGES §6, verbatim, with the threshold β ≥ 1 replaced
by β ≥ 0.37.

**Tolerances at and above the new J\*** (`lrt5` (ii); the pw5_theoremL definitions re-run at new (J, m)). Each cell gives two values:
* "num" uses the numerical min_σ F_LR at the state's u, as in PW_IMAGES;
* "cert" uses the certified F_cert, valid for all u ≤ 1 and both σ on the β-cell of the state. This is the value a rigorous-floor statement
  should quote.

The columns:
* Q_tol = F/(4|D|);
* Λ_amp (monopole only) = Q_tol/|Q|_{Λ_amp=1};
* Λ_amp (vector) = the full scalar-shaped O(1/m²) vector in the c-ratio form R1′, worst radial shift.

| state | J | m/Λ | β | min F_LR (num) | F_cert | Q_tol num / cert | Λ_amp tol, monopole num / cert | **Λ_amp tol, vector R1′ num / cert** | K_tol sh±1 (cert) | K_tol C2 (cert) |
|---|---|---|---|---|---|---|---|---|---|---|
| gap | 139 | 1.131 | 0.247 | 8.22e-5 | 1.30e-6 | 2.1e-2 / 3.3e-4 | 4023 / 63 | 285 / **5** | 357 | 4.3 |
| gap | 175 | 1.162 | 0.301 | 5.12e-5 | 1.64e-5 | 1.3e-2 / 4.1e-3 | 2612 / 835 | 299 / 95 | 6.1e3 | 170 |
| gap | **223** | 1.203 | 0.370 | 3.41e-5 | 1.92e-5 | **8.5e-3 / 4.8e-3** | **1819 / 1025** | **322 / 181** | 9.4e3 | 221 |
| gap | 250 | 1.225 | 0.407 | 2.91e-5 | 1.89e-5 | 7.3e-3 / 4.7e-3 | 1590 / 1032 | 351 / 228 | 1.0e4 | 236 |
| gap | 300 | 1.265 | 0.471 | 2.37e-5 | 1.78e-5 | 5.9e-3 / 4.5e-3 | 1344 / 1010 | 445 / 335 | 1.1e4 | 333 |
| gap | 344 | 1.299 | 0.525 | 2.10e-5 | 1.69e-5 | 5.3e-3 / 4.2e-3 | 1226 / 990 | 542 / 437 | 1.1e4 | 399 |
| gap | 400 | 1.342 | 0.589 | 1.89e-5 | 1.60e-5 | 4.7e-3 / 4.0e-3 | 1140 / 965 | 702 / 595 | 1.1e4 | 607 |
| gap | 500 | 1.415 | 0.695 | 1.69e-5 | 1.49e-5 | 4.2e-3 / 3.7e-3 | 1067 / 939 | 929 / 818 | 1.1e4 | 1.1e3 |
| gap | 700 | 1.550 | 0.877 | 1.51e-5 | 1.36e-5 | 3.8e-3 / 3.4e-3 | 1018 / 914 | 1184 / 1063 | 1.0e4 | 3.5e3 |
| gap | 857 | 1.648 | 1.000 | 1.44e-5 | 1.30e-5 | 3.6e-3 / 3.3e-3 | 1003 / 901 | 1240 / 1114 | 1.0e4 | 1.2e4 |
| gap | 10⁴ | 4.585 | 3.046 | 6.86e-6 | 5.53e-6 | 1.7e-3 / 1.4e-3 | 540 / 435 | 555 / 448 | 5.4e3 | 5.0e3 |
| gap | 10⁵ | 14.19 | 5.305 | 2.08e-6 | 1.68e-6 | 5.2e-4 / 4.2e-4 | 164 / 132 | 166 / 133 | 1.8e3 | 1.9e4 |
| gap | 10⁶ | 44.76 | 7.603 | 5.51e-7 | 4.13e-7 | 1.4e-4 / 1.0e-4 | 43 / 33 | 43 / **32** | 449 | 1.0e3 |
| 1.1 m_gap | 250 | 1.347 | 0.370 | 3.43e-5 | 1.92e-5 | 8.6e-3 / 4.8e-3 | 2301 / 1287 | 407 / 228 | 1.2e4 | 277 |
| 1.3 m_gap | 307 | 1.652 | 0.371 | 3.46e-5 | 1.92e-5 | 8.7e-3 / 4.8e-3 | 3488 / 1934 | 616 / 342 | 1.8e4 | 402 |
| 1.5Λ | 279 | 1.5 | 0.371 | 3.44e-5 | 1.92e-5 | 8.6e-3 / 4.8e-3 | 2858 / 1594 | 502 / 280 | 1.5e4 | 316 |
| 2Λ | 372 | 2 | 0.371 | 3.48e-5 | 1.92e-5 | 8.7e-3 / 4.8e-3 | 5144 / 2835 | 907 / 500 | 2.6e4 | 576 |
| 5Λ | 929 | 5 | 0.370 | 3.54e-5 | 1.92e-5 | 8.9e-3 / 4.8e-3 | 32668 / 17722 | 5783 / 3137 | 1.6e5 | 3.8e3 |

The J = 857 "num" row reproduces PW_IMAGES (1003 / 1240 / 1.1·10⁴ / 1.3·10⁴).

Reading:
* **Premise tolerances at the new J\* (twist edge, J = 223):**
  * monopole tolerance Q_tol = 4.8·10⁻³ (certified floor; 8.5·10⁻³ numerical floor). This is *larger* than at J = 857, because F_LR
    grows as β decreases while the monopole cost does not depend on J;
  * tolerated Λ_amp = 1025 (monopole only) and **181 for the full scalar-shaped vector** (certified floor; 322 numerical).
* **The binding term changes with J.** Below J ≈ 600 the full vector is bound by the ring-class costs, mainly C2 (K_tol C2 ≈ 220 at
  J = 223), not by the monopole. The tolerated vector Λ_amp therefore *decreases* towards small J: 1114 (857) → 437 (344) → 181 (223).
* **Minimum over the theorem's range** (gap, sampled J, certified floor): **Λ_amp ≤ 181 on 223 ≤ J ≤ 10⁴, ≤ 133 to 10⁵, ≤ 32 to 10⁶**.
  PW_IMAGES's numerical-floor value was 43 at 10⁶.
  * Against the J = 100 fits (Λ_amp ≈ 2–3, [NE]) the margin is ×60–90 at J\* and ×10 at J = 10⁶.
  * Index shifts: K_tol sh±1 ≥ 449 against W's K_sh = 10.
* **Why not J\* = 139.** F_LR > 0 is certified there, but at J = 139 the certified floor is 1.6 % of the numerical one. The tolerated
  vector Λ_amp collapses to 5 (285 numerical). Between J = 175 and 223 it is 95–181. So J\* = 223 (r < 1) is the sensible theorem
  threshold. Theorem L with J\* = 175 is also valid, at Λ_amp ≤ 95.
* **Old (Loc) premise numbers.** LARGESPIN_L's original sufficient condition ε_out ≤ tol(1)/4 splits at t = 1 > β for β < 1. It is
  superseded by L′'s exact per-structure costs and was not recomputed.
* **ν-space window form.** At J ≥ 223 the binding window tolerance is §4.2's minimum ε_ds(ν<10) ≈ 6.7·10⁻⁴ (J ≈ 344–500, numerical floor),
  i.e. ≈ 5·10⁻⁴ with the certified floor. That gives a sufficient window size ≈ 2.5–3.3·10⁻⁴, against 4·10⁻⁴ before ([DV], sampled /
  interpolated at J = 223). Against the measured graviton non-ring remainder 2–5·10⁻⁴ (LARGESPIN_L §7.2) this old form is marginal.
  L′ is the operative form.

---

## 6. Coverage of the near-gap band and Theorem W / Theorem L

**What remains of 345 ≤ J < 857 (twist edge):**
* **The leading-Regge floor is now [P + certified] at every J ≥ 223 at the gap**, not only at sampled J. LARGESPIN_L had [NE] there.
* **Theorem L′ now applies to every state in the band** ([D given P_L′]), with tolerated Λ_amp ≥ 437 (certified floor, sampled J ≥ 344).
* **So the band no longer lacks a theorem.** It has exactly the conditional status that J ≥ 857 had. What remains is common to the whole
  large-spin range:
  1. P_L′ magnitudes for the graviton at J ≥ 150 [M];
  2. the beyond-all-orders remainder ε_R;
  3. the tolerances are [DV] at sampled J (smooth in J; 14 gap samples between 139 and 10⁶ here);
  4. the hybrid-vs-F_LR identification [NE] (positive, ≤ 2.6·10⁻⁴ relative at J ≥ 223).

**Theorem W vs Theorem L at the gap:**
* Theorem W, with the k = 3 piece included (INDEP_CHECK2 §4.4), has J\* = 64 at (Λ_amp, K_sh, K₄) = (3, 10, 10), θ = 1. It holds on the
  sampled H_LJ states up to J = 344, at masses m ≤ 1.3Λ.
* The new Theorem L starts at J\*(m_gap) = 223. For every fixed m ≤ 1.3Λ, J\*(m) ≤ ⌈(0.37211·1.3Λ − 1)/2⌉ = 241.
* **The two now overlap at the gap on 223 ≤ J ≤ 344**, which contains W's sampled gap states J = 250, 300, 301, 344. At Λ_amp = 3 both
  premises are met there with large margin: L′ tolerates ≥ 181.
* At the gap their union covers **J ≥ 64**:
  * W on its samples (J = 64, 80, 100, 101, 120, 150, 200, 201, 250, 300, 301, 344 near the gap);
  * L′ at every J ≥ 223.
* **Still not covered by any theorem:**
  * non-sampled J in 64 < J < 223 (W is sampled-only);
  * J < 64 (W's J\* with the k = 3 piece);
  * the premises of each theorem.
* **Side benefit for Theorem W.** Its premise P0 (the leading-Regge floor) is now certified at its sampled gap states with J ≥ 139
  (F_LR > 0) and J ≥ 223 (r < 1 form). W's floor f_h is the hybrid; the hybrid-vs-F_LR step is [NE], as above.

---

## 7. Caveats and reproduction

1. **What is certified.**
   * Certified: the constants, the bound evaluations, the cell-wise r, r̃ and F_cert, the β ≥ 64 tail and the integer J\*. The lemma's
     inequalities are [P] (LARGESPIN_L §3, re-read; step 3 proved in §1.1).
   * Not certified: the tolerances of §5 ([DV]: floating-point quadrature on the published X-grid, sampled states) and the hybrid
     comparison ([NE]).
2. **"Monotone in β"** is proved for each bound piece (§1.4). r itself is not monotone across cell-width changes (4 small upticks, at
   β = 2.99, 5.98, 11.98, 23.98), which is why sups are taken over cells.
3. **M₁ at β ≳ 3.** The cell lower bounds lose 7–31 % there (β = 3–24). The F_cert/F_num of 0.75–0.81 at J = 10⁴–10⁶ in §5 mixes this loss
   with the u-worst-case of F_cert. Finer cells would recover most of it; this was not needed for the conclusions.
4. **Scope.** Δ_g = 1000 only. The exact min-norm weights ψ_k are used. The weight-realisation caveat of LARGESPIN_L §3.3 (published grid
   weight valid to J ≲ 10¹²) is unchanged.
5. **Small-β limit (β < 0.245).** F_LR > 0 is not certified there, although the true F_LR is positive and the true r ≤ 10⁻⁴. The obstacle
   is the lemma's structure (§2), not its constants.

```
cd ads4/code/lrtight && bash run_all.sh       # ≈ 10 CPU-min, one process at a time
```

| file | content | grade |
|---|---|---|
| `lrt_common.py` | exact ψ_k; R_∞, R, d₁, sumc, P_{−½} in Arb | — |
| `lrt1_constants.py` | ∫\|ψ_k\| (Sturm + bisection + Arb), Riemann cross-check, A₁,min, A₂′ | [certified] |
| `lrt2_m1.py` | M₁ lower bounds on 910 β-cells [0.1, 24]; region-(c) bound | [certified] / [P] |
| `lrt3_rcert.py` | cell-wise r, r̃, F_cert on [0.1, 64]; β ≥ 64 tail majorant; β₀ | [certified] |
| `lrt4_jstar.py` | J\*(m) tables for each β₀ (reproduces 857 … 2602 and 324) | [certified] |
| `lrt5_tolerances.py` | sanity vs direct Ω_k; Theorem L′ tolerances at new J\*; hybrid vs F_LR | [DV] / [NE] |
