# GAP_MORSE: a certified holomorphic Morse chart and the Laplace remainder for the near-gap crystal integrals

*24 Nov 2026, gap-Morse agent. Unreviewed. I made no commit (an intermediate snapshot of `code/gapmorse/` was swept into commit
840898a by another session). No existing file was edited; `gapc3/`, `gapnear/`, `spinfunc/` were only read or imported.*
* *New code:* `ads4/code/gapmorse/` (`run_all.sh`; ≈ 6.5 CPU-h used in this session, ≈ 5.5 for run_all.sh; at most 2 processes; `ulimit -v 3000000`;
  `OMP_NUM_THREADS=1`). Every script installs `gm_guard.py`, which refuses writes outside `gapmorse/cache/`; that directory's
  `.gitignore` ignores everything but itself.
* *Read-only imports:* `gapc3/gc_node` (float Morse matrix, exact crystal ball), `gapc3/gc_scalar` (validation only), and through
  `gm_ball` the same fresh in-memory copies of `spinfunc/l1_crystal, l5_spin, l6_spin_action, w0_ee_common, l9_points` as GAP_C3.
* *Numbers taken from other documents:* the tolerances `G_cert(J)` (GAP_UNIFORM §3.2) and the targets rel*_J, Π = 3, μ = 1.3 of
  GAP_C3 §4.1; GAP_C3's Δ₁ table (`gapc3/cache/gc_table.json`) for comparison.

Grades as in GAP_C3: **[P]** proved; **[certified]** computer-assisted with Arb ball arithmetic; **[D]** derived; **[D-sketch]**;
**[DV]** derived + float check; **[NE]** numerical evidence; **[M]** model/convention; **[I]** open.
Notation as GAP_C3: ε = Δ^{−1/2}, ξ = ν/Δ, Φ_ξ the Δ-independent phase, v the Morse coordinates of Φ_ξ at the crystal, h = g/g_env the
envelope-normalised amplitude, K₂ = the Laplace expansion to relative O(ε⁴), rel = the K₂ truncation error relative to
(2π/Δ)^{7/2}e^{ΔΦ*}g_env.

---

## 0. Verdict

**A holomorphic Morse chart is now constructed and certified at every node, and the K₂ remainder it gives lowers the certified Δ₁ by
×16–×47 against GAP_C3 at every node. Certified: Δ₁ = 2.0·10⁴ – 1.4·10⁵ on ξ ∈ [0.3, 0.85] and 0.9·10⁵ – 4.2·10⁵ at ξ = 0.1 (scalar
k = 2 integrand with Q_J, J ∈ {0, 2, 4, 5, 6, 7, 8}); 4.5·10⁴ – 1.0·10⁵ (ξ ∈ [0.3, 0.85]) and 2.2–3.0·10⁵ (ξ = 0.1) for the
graviton EE/C3H amplitudes at J = 0, 2 (×16–×31 below GAP_C3 where it has values). The method does not
reach Δ_g = 1000–5000. The range 1000 ≤ Δ_g < Δ₁\*(J) ≈ 0.9·10⁵ (J = 0) … 4.2·10⁵ (J = 8) is still covered only by premise B_F,
with [NE±] spot checks at 1000, 2000 (5000).** GAP_C3 estimated Δ₁ ≈ 3·10³ – 1.3·10⁴ for a chart of radius 0.1. The radius certified
here is 0.023–0.052 for the scalar amplitude and 0.014–0.032 for the graviton prefactors, and Δ₁ scales as r_c⁻².

| # | statement | grade | coverage |
|---|---|---|---|
| **TM-B** | **Hermitian-ball Taylor models.** GAP_C3's Arb Taylor models with the polydisc degree bound Σ\|c_α\|ρ^k replaced by the Bombieri norm [p_k]ρ^k. This bound is submultiplicative and bounds sup\|p_k\| on the ball \|v\|₂ ≤ ρ, so every enclosure becomes one on the Hermitian ball. At ξ = 0.5 the binding elementary factor has q = 5.8ρ (polydisc: 21ρ), and the cubic norm is 2.2 (polydisc S₃ = 10.8). New amplitude forms: (a) one exp of a sum of logs; (b) closed forms `Z·Z = −m12m15m25`, `W·W = −m34m35m45`, `Z·W = −½(m13m25m45 − m14m25m35 − m23m15m45 + m24m15m35)` for the spin contraction; (c) the X-powers of the graviton prefactor as exps of logs. Together they cut amplitude remainders by 10³–10⁵ at equal radius. | **[certified]** (enclosures); [DV] validation ≤ 5·10⁻¹⁴ (points), ≤ 2·10⁻²⁰ (exact Taylor coefficients vs the project builders) | all kinds |
| **CH** | **Certified holomorphic Morse chart.** Let P be the degree-D Taylor polynomial of Φ_ξ − Φ*, and set `W(v) = Q(v)^{1/2}v`, `Q = −2∫₀¹(1−t)∇²P(tv)dt`. Then `P = −½W·W` exactly. Explicit bounds follow from Bombieri norms: ‖Q−1‖ ≤ η, ‖W′−1‖ ≤ L, \|log det W′\| ≤ λ_J. If η, L < 1 on B_{ρ₂}, the inverse V exists and is biholomorphic B_r → V(B_r) ⊂ B_{ρ₂}, r = √(1−η)ρ₂ (contraction mapping). | **[P]** (§2.2) | any node |
| **LM** | **Lemma M (remainder).** The Taylor-model remainder R stays in the exponent (κ = 0.01–0.06), H = h(V)·det V′, and Cauchy in ε on \|ε′\| = r_c/\|u\|. With a Stokes deformation from W(real ball) to the real ball and an explicit annulus bound, `rel(Δ) ≤ t₁ + t₂ + t₃ + t₄` holds for the K₂ truncation of the integral over the real Morse ball \|v\| ≤ s_out. Each term is non-increasing in Δ. | **[P]** (§2.3) | any node |
| **S** | **Scalar k = 2 integrand with Q_J, certified Δ₁** at J ∈ {0, 2, 4, 5, 6, 7, 8} × ξ ∈ {0.1, 0.3, 0.5, 0.7, 0.85}: all 35 nodes (table below; ξ = 0.1 binds). | **[certified]** + [P] | 35 nodes |
| **G** | **Graviton EE (k = 2) and C3H (k = 3), J = 0, 2**, same chart (Φ_ξ is kind-independent), amplitude = log-factored G × the project prefactor: Δ₁ = 4.5·10⁴ – 1.0·10⁵ at ξ ∈ {0.3, 0.5, 0.7, 0.85}, 2.2–3.0·10⁵ at ξ = 0.1 (§4.2). | **[certified]** for the frozen-⟨TTO⟩ amplitude | J = 0, 2, ξ-nodes of §4.2 |
| **I** | **ξ-intervals.** δξ as an 8th Taylor-model variable (crystal and phase coefficient exact in δξ, bi-graded norms, the remainder split into v-homogeneous parts). One certificate covers a whole interval. J = 0: Δ₁ = 8.9·10⁴ on [0.4975, 0.5025], 1.1·10⁵ on [0.495, 0.505], 5.0·10⁵ on [0.4875, 0.5125], 2.7·10⁵ on [0.295, 0.305], 1.0·10⁵ on [0.695, 0.705]. J = 7: 2.5·10⁵ on [0.495, 0.505]. | **[certified]** + [P] | 6 intervals (§4.3) |
| **E″** | **Explicit eventual positivity**, GAP_C3's E′ with the new Δ₁: Δ_g\* ≤ Δ₁\*(J) = 0.9·10⁵ (J = 0), 1.1·10⁵ (2), 1.8·10⁵ (4), 1.3·10⁵ (5), 2.2·10⁵ (6), 2.5·10⁵ (7), 4.2·10⁵ (8). Premises: TAIL beyond s_out; Π ≤ 3 [M]; sup\|Γ_K₂\| bounded [NE]; ξ-continuity between nodes [D-sketch]; graviton kinds only at J = 0, 2. | [P given premises] | certified nodes |
| **F** | **A defect of the GAP_C3 evaluator, fixed here.** gc_ball's log series uses float coefficients (1/3, 1/5, …), and the rounding (≲ 10⁻¹⁷ relative) is not in the enclosure. It is numerically irrelevant for GAP_C3, but strictly outside its certificates. | [D] (found by an exact-crystal test, §3.3) | GAP_C3 records |

**Certified Δ₁(J, ξ), scalar k = 2 integrand with Q_J** (smallest Δ with certified K₂ error ≤ rel*_J; valid for every Δ ≥ Δ₁):

| J | rel*_J | ξ = 0.1 | ξ = 0.3 | ξ = 0.5 | ξ = 0.7 | ξ = 0.85 | **Δ₁\*(J) = max** | GAP_C3 max |
|---|---|---|---|---|---|---|---|---|
| 0 | 3.5·10⁻² | 8.9·10⁴ | 3.2·10⁴ | 2.2·10⁴ | 2.2·10⁴ | 3.5·10⁴ | **8.9·10⁴** | 2.8·10⁶ |
| 2 | 1.1·10⁻² | 1.1·10⁵ | 3.3·10⁴ | 2.0·10⁴ | 2.4·10⁴ | 3.5·10⁴ | **1.1·10⁵** | 3.3·10⁶ |
| 4 | 2.0·10⁻³ | 1.8·10⁵ | 4.2·10⁴ | 3.0·10⁴ | 3.2·10⁴ | 4.7·10⁴ | **1.8·10⁵** | 4.0·10⁶ |
| 5 | 5.4·10⁻³ | 1.3·10⁵ | 3.8·10⁴ | 2.4·10⁴ | 2.2·10⁴ | 4.0·10⁴ | **1.3·10⁵** | 3.8·10⁶ |
| 6 | 6.4·10⁻⁴ | 2.2·10⁵ | 6.7·10⁴ | 4.5·10⁴ | 4.5·10⁴ | 6.7·10⁴ | **2.2·10⁵** | 4.7·10⁶ |
| 7 | 4.2·10⁻⁴ | 2.5·10⁵ | 7.1·10⁴ | 4.5·10⁴ | 4.2·10⁴ | 7.1·10⁴ | **2.5·10⁵** | 5.3·10⁶ |
| 8 | 5.6·10⁻⁵ | 4.2·10⁵ | 1.1·10⁵ | 8.4·10⁴ | 9.4·10⁴ | 1.4·10⁵ | **4.2·10⁵** | 1.0·10⁷ |

**Graviton amplitudes** (frozen ⟨TTO⟩ at Δ = 1002; §4.2):

| kind | J | ξ | rel*_J | Δ₁ (certified) | ρ_h (deg 6) | r_c | M_H | GAP_C3 Δ₁ |
|---|---|---|---|---|---|---|---|---|
| C3H | 0 | 0.1 | 3.5·10⁻² | 2.2·10⁵ | 0.015 | 0.014 | 1.85 | — |
| C3H | 0 | 0.3 | 3.5·10⁻² | 6.0·10⁴ | 0.035 | 0.0315 | 5.95 | — |
| C3H | 0 | 0.5 | 3.5·10⁻² | 4.5·10⁴ | 0.035 | 0.032 | 2.38 | 1.4·10⁶ |
| C3H | 0 | 0.7 | 3.5·10⁻² | 4.7·10⁴ | 0.035 | 0.0317 | 2.41 | — |
| C3H | 0 | 0.85 | 3.5·10⁻² | 7.5·10⁴ | 0.03 | 0.0267 | 3.62 | — |
| C3H | 2 | 0.1 | 1.1·10⁻² | 3.0·10⁵ | 0.015 | 0.014 | 1.92 | — |
| C3H | 2 | 0.3 | 1.1·10⁻² | 7.9·10⁴ | 0.03 | 0.0275 | 2.38 | — |
| C3H | 2 | 0.5 | 1.1·10⁻² | 5.0·10⁴ | 0.035 | 0.032 | 1.36 | 1.5·10⁶ |
| C3H | 2 | 0.7 | 1.1·10⁻² | 5.6·10⁴ | 0.035 | 0.0317 | 1.85 | 1.8·10⁶ |
| C3H | 2 | 0.85 | 1.1·10⁻² | 8.9·10⁴ | 0.025 | 0.0227 | 0.75 | — |
| EE | 0 | 0.1 | 3.5·10⁻² | 2.2·10⁵ | 0.015 | 0.014 | 1.76 | — |
| EE | 0 | 0.3 | 3.5·10⁻² | 6.3·10⁴ | 0.03 | 0.0275 | 2.81 | — |
| EE | 0 | 0.5 | 3.5·10⁻² | 5.6·10⁴ | 0.03 | 0.0278 | 1.66 | 1.0·10⁶ |
| EE | 0 | 0.7 | 3.5·10⁻² | 5.6·10⁴ | 0.03 | 0.0276 | 1.67 | — |
| EE | 0 | 0.85 | 3.5·10⁻² | 8.4·10⁴ | 0.03 | 0.0267 | 6.22 | — |
| EE | 2 | 0.1 | 1.1·10⁻² | 2.8·10⁵ | 0.015 | 0.014 | 1.73 | — |
| EE | 2 | 0.3 | 1.1·10⁻² | 7.9·10⁴ | 0.03 | 0.0275 | 2.25 | 1.7·10⁶ |
| EE | 2 | 0.5 | 1.1·10⁻² | 5.6·10⁴ | 0.03 | 0.0278 | 0.50 | 1.0·10⁶ |
| EE | 2 | 0.7 | 1.1·10⁻² | 6.3·10⁴ | 0.03 | 0.0276 | 1.03 | 1.2·10⁶ |
| EE | 2 | 0.85 | 1.1·10⁻² | 1.0·10⁵ | 0.025 | 0.0227 | 1.26 | — |

**Coverage (details §5).**
* *Certified, given the premises of E″:* Δ_g ≥ Δ₁\*(J) for the scalar data. For the graviton kinds this holds only at J = 0, 2
  (Δ_g ≥ 2.2·10⁵ for J = 0 and ≥ 3.0·10⁵ for J = 2, again set by ξ = 0.1).
* *Covered only by premise B_F (GAP_UNIFORM U7):* 1000 ≤ Δ_g < Δ₁\*(J). Inside that range the project has only point evaluations,
  all [NE±] or [DV]: K₂ rows at 1000 and 2000 (5000 for J = 5) and the exact scalar scan to 10⁴. GAP_C3 left [1000, 2.8·10⁶–10⁷).
* No (J, μ) state is newly covered at Δ_g ≤ 5000. No window statement (ν < 50, x > 0.9, x < 0.1) changes. J ≥ 4 graviton
  amplitudes are not certified.

**What remains, in order of value.**
1. **The gap range 1000 ≤ Δ_g < Δ₁\*.** At Δ₁ Lemma M has r_c/ε ≈ 7.5 (J = 0) to 16 (J = 8): the Cauchy radius must exceed many
   Gaussian widths, for the ε⁶ term and for the Stokes/annulus tails. So Δ₁ ≈ 56/r_c² (J = 0) to 260/r_c² (J = 8). With the certified
   radius r_c ≈ 0.035–0.05 (ξ ≥ 0.3) this is ≈ 2·10⁴–1.4·10⁵. Δ = 1000 would need r_c ≈ 0.24. That is the radius on which the elementary factors themselves change by O(1)
   (1/|∇ log Nn|₂ ≈ 0.21 at ξ = 0.5 and 0.10 at ξ = 0.1), so Taylor models centred at the crystal cannot reach it. Options:
   (a) higher Laplace orders computed exactly in ball arithmetic (Taylor coefficients of H beyond degree 4). This moves the Cauchy
   order up towards the optimal-truncation floor e^{−Δr_c²/2}; by a rough estimate it would halve Δ₁; (b) a C/Arb Taylor-model implementation at degree 12–16 (the remainders, not q, limit the radius now);
   (c) for the range below: a different certified method (e.g. rigorous quadrature of the thimble integral), or K₂ rows on a gap ladder.
2. **ξ-continuity.** The 8-variable certificate works (§4.3) but costs ≈ 95 s per interval of width 0.01 and loses ×4–8 in Δ₁; covering
   x ∈ [0.1, 0.85] needs ~75 intervals per spin (~2 CPU-h per spin), feasible but not run.
3. **Graviton kinds.** EE/C3H J = 0, 2 are certified with frozen ⟨TTO⟩ vectors; J ≥ 4 (two structures, helicity rotation) and odd J
   are not; the ⟨TTO⟩ prefactor's Taylor-model remainder (sums over null-circle samples) limits the radius to ≈ 0.02–0.03.
4. The other premises are unchanged: TAIL beyond the real Morse ball of radius s_out (now s_out ≈ 0.03–0.06, i.e. ≥ 9 Gaussian widths
   at Δ₁, so the premise concerns a region of Gaussian weight ≲ e^{−40}), Π ≤ 3 [M], sup_m|Γ_{K₂,m}| bounded [NE].



## 1. What is certified, and how it differs from GAP_C3

* **Object.** As in GAP_C3 §1: at fixed ξ the crystal integral is `c(Δ,ν)∫e^{ΔΦ_ξ(z)}g(z)dz` with Φ_ξ and (scalar) g
  Δ-independent (GAP_UNIFORM Prop. S). In the Morse coordinates v (z = z* + L_Φv, the same float L_Φ as GAP_C3), the certified
  integral is over the **real Morse ball** `|v|₂ ≤ s_out`, s_out = the certified radius of the node (0.025–0.06; GAP_C3 used the box
  ‖v‖_∞ ≤ τρ/2 ≈ 0.005–0.015). Premise TAIL (the rest of the cycle is negligible) is now needed only beyond s_out, i.e. beyond
  s_out·√Δ₁ ≈ 9–20 Gaussian widths.
* **Statement (Lemma M).** For every Δ ≥ Δ₁(J, ξ): `|I_ball(Δ)/((2π/Δ)^{7/2}e^{ΔΦ*}g_env) − (c₀ + c₂ε² + c₄ε⁴)| ≤ rel*_J`, with c_j
  the exact Laplace coefficients (the code's K₂ drops the Morse defect E ≤ 10⁻¹⁴, as in GAP_C3 §1) and rel*_J GAP_C3's Lemma-T
  target `0.9·G_cert(J)/(Πμ²Λ₀²)`.
* **Why the chart helps.** GAP_C3's Lemma L3 kept the cubic f₃ in the exponent; the complex Cauchy step then multiplies by
  `e^{κ|u|²}`, κ = sup|f₃|/(τρ)², which caps τρ at ≈ 1/(2S₃) ≈ 0.04 and forces ε ≪ 0.04. With the exact chart Φ_P = −½w·w the cubic
  is gone; only the Taylor-model remainder R (|R| ≤ 10⁻⁵–10⁻⁴, O(v⁹)) stays in the exponent, with κ = 0.01–0.06. The Cauchy radius
  becomes the radius r_c on which the chart and the amplitude are certified.

## 2. Lemma M: the certified Morse chart and the remainder [P]

### 2.1 Hermitian balls instead of polydiscs (Bombieri norms) [P]
* For a homogeneous polynomial p_k of degree k in v ∈ C⁷ with coefficients c_α, the Bombieri norm is
  `[p_k]² = Σ_{|α|=k} |c_α|² α!/k!`. It equals the Frobenius norm of the symmetric coefficient tensor T (p_k(v) = T(v,…,v)).
* (B1) `|p_k(v)| ≤ [p_k] |v|^k` for all complex v (Cauchy–Schwarz in the tensor). (B2) `[p q] ≤ [p][q]` (the product's tensor is
  the symmetrisation of T_p ⊗ T_q; symmetrisation is an orthogonal projection). (B3) for degree-(k−2) contractions,
  `‖T(v^{k−2},·,·)‖_F ≤ [p_k]|v|^{k−2}`, `(Σ_i ‖T(v^{k−3},e_i,·,·)‖_F²)^{1/2} ≤ [p_k]|v|^{k−3}`.
* So replacing GAP_C3's degree bound `B_k = Σ|c_α|ρ^k` (polydisc) by `B_k = [p_k]ρ^k` turns every Taylor model of gc_ball into a
  certified enclosure on the Hermitian ball `B_ρ = {|v|₂ ≤ ρ}` (products, compositions and the remainder formulas only use
  sup-bounds and B2). `gm_ball.py` is gc_ball with this one change (plus the exact-coefficient fix of §3.3).
* Why: the Laplace weight is Euclidean (|u|₂), and the elementary factors' log-gradients are spread over the 7 directions, so a ball
  of radius ρ costs much less than the polydisc of radius ρ: q ≈ 5.8ρ (ξ = 0.5) instead of ≈ 21ρ; the cubic Bombieri norm is
  [Φ₃] = 2.23 against the polydisc S₃ = 10.8.

### 2.2 The chart
Notation: Φ − Φ* = P + R on B_ρ (Taylor model): P = Σ_{k=2}^{D} Φ_k the exact Taylor polynomial (Arb-ball coefficients; Φ₁ = 0 by
the exact crystal, hypothesis (b) of GAP_C3), Φ₂ = −½vᵀ(1+E′)v, and |R| ≤ rem_Φ on B_ρ with R = O(|v|^{D+1}).

* **Chart of P.** `Q(v) = −2∫₀¹(1−t)∇²P(tv)dt = (1+E′) + Σ_{k≥3} Q_k(v)`, `Q_k(v) = −2T_k(v^{k−2},·,·)` (complex symmetric), so
  `P(v) = −½vᵀQ(v)v`. Put `S = Q^{1/2} = Σ_n binom(½,n)(Q−1)^n` (symmetric) and `W(v) = S(v)v`. Then **P = −½ W·W exactly**.
* **Bounds at |v| ≤ t** (all from Bombieri norms of the Φ_k, B1–B3; operator and Frobenius norms alike):
  * `η(t) = ‖E′‖_F + Σ_{k≥3} 2[Φ_k]t^{k−2} ≥ ‖Q−1‖`; `c = 1 − √(1−η) ≥ ‖S−1‖`, so `N := W − id` has `|N(v)| ≤ c|v|`.
  * `dQ(t) = Σ_k 2(k−2)[Φ_k]t^{k−3} ≥ ‖D_xQ‖` (|x| = 1); `L = c + dQ·t/(2√(1−η)) ≥ ‖DN‖` and `≥ ‖DN‖_F`.
  * Trace: `tr DN = ½trE′ + p_tr(v) + (rest)`, with the exact polynomial `p_tr = −Σ_k Lap(Φ_k)/k` (Lap = Σ_i∂_i², the holomorphic
    Laplacian; derived from Euler's identity), and |rest| ≤ (c − η/2) + dQ·t(√7 g₁ + g₂ − g₁) (g₁, g₂ explicit series in η).
  * `|log det W′| ≤ λ_J := T + (−log(1−L) − L)`, using `|tr Aⁿ| ≤ ‖A‖_F²‖A‖^{n−2}`.
* **Inverse (contraction).** If η(ρ₂), L(ρ₂) < 1 and `r = √(1−η(ρ₂))·ρ₂`, then for |w| ≤ r the map v ↦ w − N(v) sends B_{ρ₂} into
  itself with Lipschitz constant L(ρ₂) < 1. Its fixed point V(w) is holomorphic, `V = W⁻¹` on B_r, W is injective on B_{ρ₂}, and
  `|V(w)| ≤ |w|/(1 − c(ρ₂))`. `det V′ = 1/det W′(V)`.
* In w: `e^{Δ(Φ−Φ*)} h dv = e^{−Δw·w/2} e^{ΔR(V(w))} H(w) dw`, `H = h(V)·det V′`, `sup_{B_{r_c}}|H| ≤ M_H(r_c) = M_h(r_v)e^{λ_J(r_v)}`,
  `r_v = r_c/(1−c(ρ₂))`.

### 2.3 The remainder (Lemma M)
Real contour: the real Morse ball |v| ≤ s_out (= the certified radius). Statement: for every Δ with the set-up's monotonicity
condition `Δa ≥ 7/2`,
`|I_{s_out}(Δ)/((2π/Δ)^{7/2}e^{ΔΦ*}g_env) − (c₀ + c₂ε² + c₄ε⁴)| ≤ rel(Δ) = t₁ + t₂ + t₃ + t₄`, ε = Δ^{−1/2}, c_j the exact Laplace
coefficients (= E_γ[H_{2j}]), where, with r_c = τr, s < r_c, s(1 + c(s)) ≤ r, x = ε/r_c and κ = R_b(r_v)/r_c²
(R_b(t) = rem_Φ(t/ρ)^{D+1}):
* `t₁ = min_K M_H [Σ_{k=6..K even} x^k (1−2κ)^{−(7+k)/2} E|u|^k + x^{K+2}(1−2κ)^{−(9+K)/2}E|u|^{K+2}/(1−(s/r_c)²)]`
  (Cauchy in ε on |ε′| = r_c/|u| for `F(ε,u) = H(εu)e^{R(V(εu))/ε²}`, |F| ≤ M_H e^{κ|u|²}; odd orders vanish on the symmetric ball;
  R enters at ε^{D−1}, so F₀, F₂, F₄ = H₀, H₂, H₄; E|u|⁶ = 693 for u ~ N(0, 1₇)).
* `t₂ = M_H Σ_{k=0,2,4} x^k E[|u|^k; |u| > s/ε]` (K₂ polynomial outside the ball).
* `t₃ = M_H(r)|S⁶|c(s)s⁷(1+L(s))⁶e^{−Δa}(Δ/2π)^{7/2}`, `a = s²β/2 − R_b(r_{v,s})`, `β = 1 − 2c(s) − c(s)²`: Stokes. The real ball
  and its image W(real ball) are joined by G(v,t) = v + tN(v) inside B_r; the holomorphic 7-form is closed, so the two integrals differ
  by the integral over G(∂ball × [0,1]), where Re(w·w) ≥ s²β, and Hadamard bounds the Jacobian.
* `t₄ = Σ_shells M_h(s_{j+1})(1−θ_{j+1})^{−7/2} P(χ²₇ > Δ(1−θ_{j+1})s_j²)`: the real annulus s < |v| ≤ s_out, where
  `Re(Φ − Φ*) ≤ −(1−θ(|v|))|v|²/2`, `θ(t) = ‖E′‖_F + Σ_k 2[Φ_k]t^{k−2} + 2R_b(t)/t²`.
* Every term is non-increasing in Δ (t₃ once Δa ≥ 7/2), so `rel(Δ₁) ≤ rel*` implies it for all Δ ≥ Δ₁.

### 2.4 What fixes Δ₁ (structure of the bound) [D]
* At Δ₁ the interior term t₁ ≈ 693·M_H·(1−2κ)^{−6.5}(ε/r_c)⁶ dominates. The Stokes term t₃ and the annulus t₄ need s ≳ 6ε with
  s < r_c. So **Δ₁ ∝ r_c⁻²·(M_H/rel*)^{1/3}**. At Δ₁, r_c/ε ≈ 7.5 for J = 0 (rel* = 3.5·10⁻²) and ≈ 16 for J = 8
  (rel* = 5.6·10⁻⁵).
* r_c is limited by (i) the amplitude remainder (the Taylor model of h is accurate only up to ρ ≈ 0.04–0.06 at ξ ≥ 0.3 and 0.03 at
  ξ = 0.1), (ii) κ (the phase remainder; degree-10 phase records reduce it ×7–10), (iii) η(ρ₂) (the chart shrinks the radius by
  √(1−η) ≈ 0.85).
* Two factors make up M_H: the amplitude's own growth `e^{|∇h|r}` (|∇h|₂ = 7.2 at ξ = 0.5 and 18 at ξ = 0.1, for J = 0) and the
  Jacobian `e^{λ_J}` ≈ 1.3–1.7. Together M_H ≈ 2.1–2.7 for J = 0 and 0.3–2.5 for J ≥ 2 (envelope-normalised).

## 3. Evaluator: what changed relative to GAP_C3 [certified]

### 3.1 Bombieri degree bounds (Hermitian ball)
`gm_ball.py`: one function changed (`degbounds`), §2.1. Polydisc mode kept (`GM_NORM=poly`). Phase factor ranges on B_ρ at ξ = 0.5:
q = 0.161 / 0.223 / 0.290 / 0.362 / 0.440 at ρ = 0.03 … 0.07 (binding factor Nn; GAP_C3's polydisc: q = 0.43 at ρ = 0.02).

### 3.2 Amplitude in exp–log form, closed-form spin contraction
* `gm_scalar.py` (copy of gc_scalar, amplitude rewritten; phase unchanged):
  * `g/g* = exp(Σ_k p_k Log(f_k/f_k*)) · num/num*`: one exponential of a sum of the same elementary logs as the phase, the m-factors
    of f_sub moved into the exponent, only the subtraction numerator `num = m14m23 − m12m34` kept as a polynomial factor. (GAP_C3
    multiplied power series of (f/f*)^p; on the Hermitian ball with p up to −6 that construction gives an amplitude remainder ≈ 0.44
    at ρ = 0.03, ξ = 0.5; the exp–log form gives 3·10⁻⁶.)
  * **Spin contraction Q_J.** New closed forms (from dot(x_ab, x_ab) = −m_ab):
    `Z·Z = −m12m15m25`, `W·W = −m34m35m45`, `Z·W = −½(m13m25m45 − m14m25m35 − m23m15m45 + m24m15m35)`.
    Hence `Q_J = (ZZ·WW)^{J/2} P_J(cos θ)`, with `(ZZ·WW)^{J/2}` joining the exponent and `cos θ` = a sum of four exp-of-log-sums;
    P_J is expanded exactly around cos θ* = 1 − 2ξ² (exact Legendre coefficients). The direct Legendre sum of GAP_C3 / l1_crystal
    loses the cancellations between ZW^p ZZWW^{(J−p)/2} and between the lightcone components of Z, W in every remainder:
    at J = 6, ξ = 0.3, ρ = 0.04 the amplitude remainder drops from 2.0 to 9·10⁻⁴ (same exact Taylor polynomial).
* Validation (point backend, 50 random complex points |v| ~ 0.03): |h_new/h_old − 1| ≤ 5·10⁻¹⁴ against gapc3/gc_scalar (itself
  validated against l1_crystal) at (ξ, J) = (0.5,2), (0.3,4), (0.1,7), (0.85,8), (0.5,5), (0.3,6) [DV].

### 3.3 A small defect of the GAP_C3 evaluator (found and fixed here)
`gc_ball.log` builds the series Log(1+t) = Σ(−1)^{n+1}tⁿ/n with **float** coefficients 1/3, 1/5, 1/6, 1/7: they are rounded, and
the Taylor-model remainder does not account for the rounding. The enclosures of GAP_C3 are therefore off by a perturbation of relative
size ≲ 10⁻¹⁷ (vanishing to third order at the crystal). Numerically irrelevant for GAP_C3's conclusions, but outside the enclosure.
It became visible here: with ξ as a Taylor-model variable the crystal gradient's δξ² coefficient came out 2·10⁻¹⁵ ≠ 0. `gm_ball`
uses exact Arb coefficients for log and for the binomial series of powers. No existing file was edited.

### 3.4 Graviton amplitudes
`gm_spin.py`: g = exp(log-factored G_EE or G_C3H) × pref. Here G uses the same elementary logs as gm_scalar, with l9/w0's 3pt
exponents, −4 Log v042 and (EE) e^{lf} = num/(m13 m24). pref is the project's prefactor (w0 `_contract` of E⊗E, or q·o_H with l9's
`odd_structure` 'H' and the ⟨TTO⟩ samples of l5_spin), run on gm_ball. The one change is that the powers of X20 = −m20, X40 = −m40,
X24 = −m24 in it (1/(X20X40) and (X24X20X40³)^{−1/2}, …) are exps of the already-computed logs, not recip/powr of products. For C3H
that cuts the prefactor's range q at ρ = 0.02 from 0.36 to 0.10, and the remainder from 8.8 to 0.006. Validation against the project
builders (`w0.make_build_even_EE`, `l9.build_spin`) on the same backend: the exact Taylor coefficients of g/g(z*) agree to ≤ 2·10⁻²⁰
relative (EE and C3H, J = 0, 2, ξ = 0.3, 0.5). The ⟨TTO⟩ vectors are the float null vectors at Δ = 1002 (frozen, as in GAP_C3 §4.3).

## 4. Results

### 4.1 Scalar k = 2 integrand with Q_J: certified Δ₁ (`gm_run` → `cache/gm_cert_J*.json`, `gm_cert phase` → `cache/gm_phase.json`, `gm_table` → `cache/gm_table.{json,txt}`)
Nodes as GAP_C3: J ∈ {0, 2, 4, 5, 6, 7, 8} × ξ ∈ {0.1, 0.3, 0.5, 0.7, 0.85}, **all 35 certified** (every record used has bad = 0, an
enclosure of |∇Φ(z*)| ≤ 4·10⁻²⁹ containing 0, ‖E′‖_F ≤ 9·10⁻¹⁵). Taylor models of degree 8 (amplitude and phase) on Hermitian balls
ρ ∈ {0.03, …, 0.07} (ξ = 0.1: {0.02, …, 0.04}), plus J-independent degree-10 phase records. For each node, every pair (phase record
at ρ_Φ ≥ amplitude record at ρ_h) and a grid of (ρ₂, τ, s/r_c) is tried; Δ₁ is the smallest Δ on a 40-per-decade grid with
certified rel ≤ rel*_J (and every term non-increasing from there on).

| J | rel*_J | ξ = 0.1 | ξ = 0.3 | ξ = 0.5 | ξ = 0.7 | ξ = 0.85 | **Δ₁\*(J) = max** | GAP_C3 max |
|---|---|---|---|---|---|---|---|---|
| 0 | 3.5·10⁻² | 8.9·10⁴ | 3.2·10⁴ | 2.2·10⁴ | 2.2·10⁴ | 3.5·10⁴ | **8.9·10⁴** | 2.8·10⁶ |
| 2 | 1.1·10⁻² | 1.1·10⁵ | 3.3·10⁴ | 2.0·10⁴ | 2.4·10⁴ | 3.5·10⁴ | **1.1·10⁵** | 3.3·10⁶ |
| 4 | 2.0·10⁻³ | 1.8·10⁵ | 4.2·10⁴ | 3.0·10⁴ | 3.2·10⁴ | 4.7·10⁴ | **1.8·10⁵** | 4.0·10⁶ |
| 5 | 5.4·10⁻³ | 1.3·10⁵ | 3.8·10⁴ | 2.4·10⁴ | 2.2·10⁴ | 4.0·10⁴ | **1.3·10⁵** | 3.8·10⁶ |
| 6 | 6.4·10⁻⁴ | 2.2·10⁵ | 6.7·10⁴ | 4.5·10⁴ | 4.5·10⁴ | 6.7·10⁴ | **2.2·10⁵** | 4.7·10⁶ |
| 7 | 4.2·10⁻⁴ | 2.5·10⁵ | 7.1·10⁴ | 4.5·10⁴ | 4.2·10⁴ | 7.1·10⁴ | **2.5·10⁵** | 5.3·10⁶ |
| 8 | 5.6·10⁻⁵ | 4.2·10⁵ | 1.1·10⁵ | 8.4·10⁴ | 9.4·10⁴ | 1.4·10⁵ | **4.2·10⁵** | 1.0·10⁷ |

(In brackets in `cache/gm_report.txt`: GAP_C3's value at the same node; the improvement is ×16–×47 at every node.)

Binding set-ups (all records and set-ups in `cache/gm_table.json`; selection):

| J | ξ | Δ₁ | ρ_Φ (deg) | ρ_h | ρ₂ | r_c | η(ρ₂) | λ_J(ρ₂) | M_H | κ | t₁ | t₃ | t₄ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.1 | 8.9·10⁴ | 0.03 (10) | 0.03 | 0.03 | 0.0254 | 0.28 | 0.45 | 2.69 | 0.031 | 2.0e-02 | 8.0e-03 | 5.7e-04 |
| 0 | 0.5 | 2.2·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 2.49 | 0.018 | 1.6e-02 | 9.9e-03 | 7.6e-04 |
| 2 | 0.1 | 1.1·10⁵ | 0.03 (10) | 0.03 | 0.03 | 0.0254 | 0.28 | 0.45 | 2.54 | 0.031 | 8.8e-03 | 4.1e-04 | 1.6e-05 |
| 2 | 0.5 | 2.0·10⁴ | 0.07 (10) | 0.07 | 0.063 | 0.0522 | 0.31 | 0.53 | 0.56 | 0.035 | 5.3e-03 | 5.6e-03 | 2.6e-04 |
| 4 | 0.1 | 1.8·10⁵ | 0.03 (10) | 0.03 | 0.027 | 0.0234 | 0.25 | 0.38 | 1.87 | 0.012 | 1.9e-03 | 2.7e-05 | 2.5e-06 |
| 4 | 0.5 | 3.0·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 0.62 | 0.018 | 1.5e-03 | 6.2e-05 | 2.4e-06 |
| 5 | 0.1 | 1.3·10⁵ | 0.03 (10) | 0.03 | 0.027 | 0.0234 | 0.25 | 0.38 | 1.65 | 0.012 | 5.1e-03 | 2.2e-04 | 1.7e-05 |
| 5 | 0.5 | 2.4·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 0.51 | 0.018 | 2.7e-03 | 1.1e-03 | 6.1e-05 |
| 6 | 0.1 | 2.2·10⁵ | 0.03 (10) | 0.03 | 0.03 | 0.0254 | 0.28 | 0.45 | 1.61 | 0.031 | 6.2e-04 | 3.0e-06 | 2.5e-07 |
| 6 | 0.5 | 4.5·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 0.86 | 0.018 | 6.0e-04 | 2.2e-06 | 1.0e-07 |
| 7 | 0.1 | 2.5·10⁵ | 0.03 (10) | 0.03 | 0.027 | 0.0234 | 0.25 | 0.38 | 1.12 | 0.012 | 3.9e-04 | 3.2e-06 | 4.1e-07 |
| 7 | 0.5 | 4.5·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 0.54 | 0.018 | 3.7e-04 | 1.3e-06 | 6.2e-08 |
| 8 | 0.1 | 4.2·10⁵ | 0.03 (10) | 0.03 | 0.03 | 0.0254 | 0.28 | 0.45 | 0.98 | 0.031 | 5.3e-05 | 1.3e-09 | 1.3e-10 |
| 8 | 0.5 | 8.4·10⁴ | 0.06 (10) | 0.06 | 0.06 | 0.0503 | 0.30 | 0.49 | 0.54 | 0.018 | 5.2e-05 | 2.9e-07 | 2.2e-08 |

* The interior term t₁ binds everywhere. At J ≤ 2 the Stokes term t₃ is 0.05–1 × t₁, which forces s ≈ 0.7–0.8 r_c; at J ≥ 4 it is
  negligible.
* **ξ = 0.1 binds for every J** (r_c = 0.023–0.025 against 0.042–0.050 at ξ = 0.3–0.5): there the amplitude varies fast
  (|∇h|₂ = 18, [h₂] = 145 against 7.2 / 22 at ξ = 0.5) and the phase's Taylor coefficients grow ≈ 4× per degree (≈ 1.8× at ξ = 0.5).
  Nodes with ξ < 0.1 (x ∈ [0.05, 0.1)) were not computed and would be worse.
* Real contour certified: the real Morse ball of radius s_out = ρ_h (0.03 at ξ = 0.1, 0.045–0.06 otherwise), i.e.
  s_out√Δ₁ ≈ 9–20 Gaussian widths.

### 4.2 Graviton amplitudes EE (k = 2) and C3H (k = 3), J = 0, 2 (`gm_spin` → `cache/gm_spin.json`; `gm_report` → `cache/gm_spin_table.json`)
Amplitude Taylor models of degree 6 (degree 8 costs ≈ 6 min per node), ρ_h ∈ {0.015, 0.02, 0.025, 0.03} (C3H also 0.035; ξ = 0.1:
{0.01, 0.015}). The phase comes from the scalar records (Φ_ξ is kind-independent). ⟨TTO⟩ frozen at Δ = 1002 (GAP_C3 §4.3 caveat
unchanged).

| kind | J | ξ | rel*_J | Δ₁ (certified) | ρ_h (deg 6) | r_c | M_H | GAP_C3 Δ₁ |
|---|---|---|---|---|---|---|---|---|
| C3H | 0 | 0.1 | 3.5·10⁻² | 2.2·10⁵ | 0.015 | 0.014 | 1.85 | — |
| C3H | 0 | 0.3 | 3.5·10⁻² | 6.0·10⁴ | 0.035 | 0.0315 | 5.95 | — |
| C3H | 0 | 0.5 | 3.5·10⁻² | 4.5·10⁴ | 0.035 | 0.032 | 2.38 | 1.4·10⁶ |
| C3H | 0 | 0.7 | 3.5·10⁻² | 4.7·10⁴ | 0.035 | 0.0317 | 2.41 | — |
| C3H | 0 | 0.85 | 3.5·10⁻² | 7.5·10⁴ | 0.03 | 0.0267 | 3.62 | — |
| C3H | 2 | 0.1 | 1.1·10⁻² | 3.0·10⁵ | 0.015 | 0.014 | 1.92 | — |
| C3H | 2 | 0.3 | 1.1·10⁻² | 7.9·10⁴ | 0.03 | 0.0275 | 2.38 | — |
| C3H | 2 | 0.5 | 1.1·10⁻² | 5.0·10⁴ | 0.035 | 0.032 | 1.36 | 1.5·10⁶ |
| C3H | 2 | 0.7 | 1.1·10⁻² | 5.6·10⁴ | 0.035 | 0.0317 | 1.85 | 1.8·10⁶ |
| C3H | 2 | 0.85 | 1.1·10⁻² | 8.9·10⁴ | 0.025 | 0.0227 | 0.75 | — |
| EE | 0 | 0.1 | 3.5·10⁻² | 2.2·10⁵ | 0.015 | 0.014 | 1.76 | — |
| EE | 0 | 0.3 | 3.5·10⁻² | 6.3·10⁴ | 0.03 | 0.0275 | 2.81 | — |
| EE | 0 | 0.5 | 3.5·10⁻² | 5.6·10⁴ | 0.03 | 0.0278 | 1.66 | 1.0·10⁶ |
| EE | 0 | 0.7 | 3.5·10⁻² | 5.6·10⁴ | 0.03 | 0.0276 | 1.67 | — |
| EE | 0 | 0.85 | 3.5·10⁻² | 8.4·10⁴ | 0.03 | 0.0267 | 6.22 | — |
| EE | 2 | 0.1 | 1.1·10⁻² | 2.8·10⁵ | 0.015 | 0.014 | 1.73 | — |
| EE | 2 | 0.3 | 1.1·10⁻² | 7.9·10⁴ | 0.03 | 0.0275 | 2.25 | 1.7·10⁶ |
| EE | 2 | 0.5 | 1.1·10⁻² | 5.6·10⁴ | 0.03 | 0.0278 | 0.50 | 1.0·10⁶ |
| EE | 2 | 0.7 | 1.1·10⁻² | 6.3·10⁴ | 0.03 | 0.0276 | 1.03 | 1.2·10⁶ |
| EE | 2 | 0.85 | 1.1·10⁻² | 1.0·10⁵ | 0.025 | 0.0227 | 1.26 | — |

* The ⟨TTO⟩ prefactor limits the radius. It is a sum over 5×5×(2J+1) null-circle samples of high-degree polynomials with strong
  cancellation (sample sup ≈ 17, contracted value ≈ 0.1–1), so its Taylor-model remainder grows fast. The envelope-normalised
  remainder at ρ = 0.03 is 0.12–0.3 at ξ = 0.5–0.7 and 0.6–3.7 at ξ = 0.3, 0.85; at 0.035 (C3H) it is 0.5–7. That puts r_c at
  ≈ 0.023–0.032, against 0.04–0.05 for the scalar, and Δ₁ ×2–3 above the scalar value at the same node.
* ξ = 0.1 again binds. The amplitude varies fastest there (|∇h|₂ ≈ 18 for EE, 25 for C3H), and only ρ_h ≤ 0.015 is usable.

### 4.3 ξ-intervals: δξ as an 8th Taylor-model variable (`gm_cert8`, `gm_ball8` → `cache/gm_int_J*.json`; `gm_table int`)
* **Construction [certified].** The crystal z*(ξ) (a = ξ/(1+√(1−ξ²))) and the phase coefficient iξ/2 are Taylor models in δ = ξ − ξ₀
  on the disc |δ| ≤ r_x; degree bounds are bi-graded (Bombieri in v, plain in δ). The Morse matrix stays the float one at ξ₀, so the
  Hessian is −(1 + E′(δ)), E′ = O(δ), carried in ‖E′‖_F. For each fixed δ the remainder R(v, δ) (sup ≤ rem, vanishing to *total*
  order D+1) is split into v-homogeneous parts: R₀, R₁ drop exactly (P₁ + R₁ ≡ 0 because the crystal is exact for every ξ: the
  certified v-linear coefficients enclose 0 for every power of δ); R₂…R₆ join the chart polynomial with
  `[R_k]_B ≤ √C(k+6,6)·rem/ρ^k` (Bombieri² = dim H_k × the sphere mean of |p|²); R₇₊ stays in the exponent (|R₇₊| ≤ 8·rem·(|v|/ρ)⁷).
  The amplitude bound is divided by λ_min = min_δ env(ξ)/env(ξ₀) (certified from the δ-part of the envelope factor). Lemma M then
  holds verbatim **for every ξ in [ξ₀ − r_x, ξ₀ + r_x]**.
* Results:

| J | ξ-interval | ρ | ‖E′‖_F (uniform) | λ_min | Δ₁ (whole interval) | node Δ₁ at ξ₀ |
|---|---|---|---|---|---|---|
| 0 | [0.2950, 0.3050] | 0.02 | 0.047 | 0.862 | 2.7·10⁵ | 3.2·10⁴ |
| 0 | [0.4975, 0.5025] | 0.025 | 0.016 | 0.962 | 8.9·10⁴ | 2.2·10⁴ |
| 0 | [0.4950, 0.5050] | 0.025 | 0.032 | 0.923 | 1.1·10⁵ | 2.2·10⁴ |
| 0 | [0.4875, 0.5125] | 0.02 | 0.095 | 0.798 | 5.0·10⁵ | 2.2·10⁴ |
| 0 | [0.6950, 0.7050] | 0.025 | 0.032 | 0.950 | 1.0·10⁵ | 2.2·10⁴ |
| 7 | [0.4950, 0.5050] | 0.025 | 0.032 | 0.887 | 2.5·10⁵ | 4.5·10⁴ |

* Cost ≈ 95 s per interval (degree 8, 8 variables, 12 870 coefficients). Covering ξ ∈ [0.1, 0.85] at width 0.01 is ≈ 75 intervals,
  ≈ 2 CPU-h per spin; not run. The interval Δ₁ exceeds the node Δ₁ by ×4–8 (width ≤ 0.01) and ×23 (width 0.025): the R₂…R₆
  corrections grow like rem/ρ^k and the Hessian drift ‖E′(δ)‖ like r_x, which keeps ρ ≲ 0.025.

## 5. Coverage (details)

* **What Lemma M + the tables give [P + certified].** At each computed node (J, ξ) and for every Δ ≥ Δ₁(J, ξ), the K₂ truncation error
  of the crystal integral over the real Morse ball is ≤ rel*_J, the size GAP_C3 §4.1 derived from Lemma T (GAP_UNIFORM §3.1) for
  Λ ≫ Λ₀. At Δ ≥ Δ₁ ≥ 2·10⁴ one has (Λ/Λ₀)² ≥ 400, so the 0.9 in rel*_J (room for sup|Γ_K₂|) holds for every J here, including
  J = 8 (F = 0.85 at Λ₀, GAP_UNIFORM §3.4), as long as sup_m|Γ_{K₂,m}| stays bounded.
* **Explicit eventual positivity (GAP_C3's E′, with Δ₁ ×16–×47 smaller) [P given (i)–(v)].** Given the premises below, the near-gap
  states of spin J (at the computed x-nodes) are positive for every **Δ_g ≥ Δ₁\*(J) = max_ξ Δ₁(J, ξ)**: ≈ 0.9·10⁵ (J = 0), 1.1·10⁵ (2), 1.8·10⁵ (4), 1.3·10⁵ (5), 2.2·10⁵ (6), 2.5·10⁵ (7),
  4.2·10⁵ (8), set by ξ = 0.1 (on ξ ∈ {0.3, …, 0.85} alone: ≤ 3.5·10⁴, 3.5·10⁴, 4.7·10⁴, 4.0·10⁴, 6.7·10⁴, 7.1·10⁴, 1.4·10⁵).
  Premises, unchanged from GAP_C3 §0 E′: (i) TAIL beyond the real Morse ball of radius s_out [D-sketch]; (ii) Π ≤ 3 [M];
  (iii) sup_m|Γ_{K₂,m}| bounded [NE]; (iv) continuity in ξ between nodes (certified on the intervals of §4.3 only; [D-sketch]
  elsewhere); (v) **for the graviton kinds that enter Ω^EE, only EE/C3H J = 0, 2 are certified** (all five ξ-nodes;
  frozen ⟨TTO⟩; Δ_g ≥ 2.2·10⁵ / 3.0·10⁵); for J ≥ 4 the scalar k = 2 integrand stands in for the phase data, as in GAP_C3 §0 C3-S.
* **Combined with the existing numerics.** The project has, for these spins, error-budgeted K₂ rows at Δ_g = 1000 (NEARGAP2;
  direct rows for J ≥ 7) and 2000 (GAP_NEARGAP; also 5000 for J = 5), and the exact (Arb) scalar chart-A scan at
  Δ_g = 1000 … 10⁴ (GAP_UNIFORM §5.2). Those are point evaluations, and their K₂ truncation is itself [NE±], not certified. So:
  * **certified (given (i)–(v)):** Δ_g ≥ Δ₁\*(J), i.e. Δ_g ≳ 10⁵ (J ≤ 5) and ≳ 2.5·10⁵ (J = 6, 7), 4.2·10⁵ (J = 8);
  * **covered only by premise B_F (GAP_UNIFORM U7):** 1000 ≤ Δ_g < Δ₁\*(J), with [NE±] spot checks at 1000, 2000 (5000 for J = 5) and the
    [DV] exact scalar to 10⁴. GAP_C3 left [1000, 2.8·10⁶–10⁷); the B_F-only range is now [1000, 0.9–4.2·10⁵): its upper end is lower by
    ×24–×31, but it still spans 2–2.6 decades. For the graviton kinds at J = 0, 2 it is [1000, 2.2–3.0·10⁵).
* **Not reachable by this method at Δ_g = 1000.** At Δ₁ Lemma M has r_c/ε ≈ 7.5 (J = 0) to 16 (J = 8). At Δ = 1000 (ε = 0.032) that
  means r_c ≳ 0.24, while the integrand's elementary factors already vary by O(1) on Hermitian balls of radius
  1/|∇ log Nn|₂ ≈ 0.2 (ξ = 0.5) and 0.1 (ξ = 0.1). GAP_C3's estimate "Δ₁ ≈ 3·10³ (J = 0) to 1.3·10⁴ (J = 7) with a chart of radius 0.1" is consistent with what is found
  here at radius 0.05 (2.2·10⁴ / 4.5·10⁴ at ξ = 0.5, scaling as r_c⁻²), but radius 0.1 is not certifiable with the present Taylor
  models (degree 8 in Python; the remainders, not the branch checks, are what stop at ρ ≈ 0.06).

## 6. Reproduce

```
cd ads4/code/gapmorse && sh run_all.sh        # ~5.5 CPU-h, <= 2 processes; outputs in cache/
```

| file | content | grade |
|---|---|---|
| `gm_guard.py` | write guard (only `gapmorse/cache`) | — |
| `gm_ball.py` | gc_ball with Bombieri degree bounds (Hermitian ball) and exact series coefficients | [certified] (enclosures) |
| `gm_ball8.py` | the same with δξ as 8th variable, bi-graded bounds | [certified] |
| `gm_scalar.py` | log-factored integrand; exp–log amplitude; closed-form cos θ spin contraction; phase-only mode | [DV] (validation) |
| `gm_cert.py` | per-node Taylor-model data: q, ∇Φ, ‖E′‖_F, [Φ_k], p_tr norms, rem_Φ, [h_k], rem_h (`phase` mode: degree-10 phase records) | [certified] |
| `gm_cert8.py` | ξ-interval data (§4.3) | [certified] |
| `gm_chart.py` | Lemma M in Arb: chart bounds, inverse, M_H, t₁–t₄, Δ₁ | [P] |
| `gm_run.py` | scalar batch → `cache/gm_cert_J*.json` | [certified] |
| `gm_spin.py` | EE / C3H amplitudes (log-factored G × project prefactor); `check` validates against the project builders | [certified] (frozen ⟨TTO⟩) |
| `gm_table.py` | Δ₁ tables (`int`: ξ-intervals) → `cache/gm_table.*`, `cache/gm_int_table.json` | [certified] + [M] (targets) |
| `gm_report.py` | report, comparison with GAP_C3, graviton table → `cache/gm_report.txt`, `cache/gm_spin_table.json` | — |
| `gm_laneS.sh`, `gm_lane_int.sh` | graviton-amplitude lane (EE or C3H) and ξ-interval lane, called by run_all.sh | — |
| `gm_pause.sh` | runs a short analysis with one compute lane SIGSTOPped (≤ 2 running processes) | — |

*Cache provenance.* `cache/floatlog/`, `cache/abhorner/`, `cache/spin_v1/` hold superseded records made before the fixes of §3.3
(float log coefficients), the closed-form spin contraction and the graviton-prefactor powers. They are not used in any table, and
run_all.sh does not regenerate them.
