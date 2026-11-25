# LOWNU_SCOPE: a low-ν cubic functional with provable heavy positivity — scoping

*25 Nov 2026, low-ν scoping agent. Unreviewed, not committed. No existing file was edited; `T0_CERT.md`, `code/t0cert/` were not touched.*
* *New code:* `ads4/code/lownu/` (`run_lp.sh`, `ln1`–`ln7`); caches in `code/lownu/cache/`, whose `.gitignore` ignores everything but itself.
  Every script writes only through `ln_common.guard` (refuses writes outside `lownu/cache/`).
* *Read-only imports:* `robust/rcommon.py`, `robust/lp.py`, `robust/nf_lp.py` (no writes at import), `robust/functionals_noforward.json`.
* *Compute:* ≈ 0.4 CPU-h, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.
* *Inputs from the coordinator/reviewer* (dilation identity, analytic de-shift weight, K(L) formula, Gamma factorisation, ANGULAR_ROUTE)
  were re-checked where stated below; what was not re-checked is marked.

Grades: **[P]** proved; **[certified]**; **[D]** derived; **[D-sketch]**; **[DV]** derived + float check; **[NE]** numerical evidence;
**[M]** model / heuristic; **[I]** open. Notation: `a = p_max = ν_max/Λ` (support edge of the target weights), `Λ′ = aΛ`,
`ξ = ν/Δ`, `τ = √(t₄²+α_T²)`, `u′ = Λ′²/m²`. Flat units in the LP are `M′ = Λ′` (p′ = ν/Λ′ ∈ [0,1]).

---

## 0. Verdict

**Partially viable. The flat side (LF1, LF2) passes easily and gives better numbers than plain compression. The AdS remainder
(LF3), via the contour route that motivated the idea, is blocked at Δ_g = 1000. Two exact-representation routes, the plane-wave
route and the authors' angular route, avoid the blocking steps. Each still has one decisive open inequality, and neither can be
closed by Thursday with confidence. A rigorous 10⁻⁶ bound at Δ_g = 1000 is not in reach on any route: it would need a provable
ξ_max ≈ 0.18–0.2, against 0.088–0.107 available.**

| gate | status | key numbers |
|---|---|---|
| **LF1** sensitivity vs p_max | **passed** [P + DV] | **Dilation identity [P]:** compressing NF02 is NF02 at effective gap Λ′ = aΛ, so τ ∝ a⁻⁴. **New LP**, using the fact that only u′ ≤ a² must be positive: a **uniform ×19 gain** over dilation at every a [NE]. τ at Δ_g = 1000 (K(L) included): **2.2·10⁻³ / 1.56·10⁻⁴ / 3.3·10⁻⁵ / 1.11·10⁻⁵** for a = 0.025 / 0.05 / 0.075 / 0.1 (dilated NF02: 3.8·10⁻² / 2.9·10⁻³ / 6.3·10⁻⁴ / 2.1·10⁻⁴) |
| **LF2** flat positivity, all masses/spins/helicities | **passed** for dilated NF02 **[certified]**; **plausible** for the LP functionals [NE] | Dilated NF02: the existing certificate covers every needed state for a ≤ 0.612. Its only uncertified NF02 region (+−, J ≥ 400, m < 1.633M′) lies below the gap. Margin 3/20 certified, 0.2 dense. LP functionals: dense min 𝔉/N = 0.1991–0.2000 at every a (u′ ∈ (0, a²], J ≤ 6000, scaling limit to b̂′ = 2000, 1/b̂³ coefficient > 0). Certifying them reuses `nf_cert` (≈ 3–6 ah) |
| **LF3** relative AdS remainder, uniform in states and polarisations | **blocked** at Δ_g = 1000 on the contour route; **plausible but open** on the plane-wave and angular routes | Prop. D′ controls only the **far-region** (non-perturbative) piece. It does so for J = 0 kinds, at ξ < 0.107 (and below 0.088 once the Δ^N prefactor at Δ = 1000 is counted). The **near-region** integral, i.e. E itself, has no certified evaluation at Δ = 1000. Low ν removes CT (O1) and the fold (O4), but not the near-region resummation, J-uniformity (O3) or large-b̂ cancellation (O2); O2 starts earlier (b̂ ~ 1/a) |
| **LF4** physical bound | **estimated** | Light side: the exact-de-shift K(L) is reproduced (0.99674 / 0.98686 at L = 0.1Λ / 0.05Λ for NF02). The k = 2 analytic de-shift weight is reproduced at L = 0.05Λ and 0.1Λ. Its e^{πν}-weighted tail forces **L ≲ 0.06Λ** at Δ_g = 1000 if it must sit inside the Prop.-D′ core. **Best defensible estimate: τ ≈ 0.8–1.6·10⁻⁴ (LP, a = 0.05–0.06) at Δ_g = 1000**; 10⁻⁶ at a = 0.06 needs Δ_g ≈ 3000 |

**Comparison with the angular route (ANGULAR_ROUTE.md).**
* **Which is closer by Thursday.** The angular route. It never expands the block: per-harmonic positivity comes from the orthogonality
  of continuous dual Hahn polynomials (SP2), so it is uniform in spin and needs no near-region resummation, no contour theorem and no
  1/Δ remainder. Those are the items that block the contour route (§3). The plane-wave route (§3.4) is also exact and uniform in
  spin, but its decisive input is a pointwise inequality for the functional's plane-wave image, 𝓜(p) ⪰ c·f(β(p)) for M ≥ 800. My
  heuristic (§3.4) says this holds only for a ≲ 0.03 at Λ = 1000, which would push τ to ~10⁻³.
* **What the angular route needs.** Its items 1–3.
  * **Item 1 (heavy-family transition matrices).** My analysis says nothing about them. Every computation in the record and here
    treats one heavy conformal family at a time (a block's own action). The first retarded-time moment between distinct families is
    a new object. The selection rule |λ−h| ≥ 2 kills the harmonics that carry the light-side cubic term only if it holds for the full
    transition matrix, not just for momentum eigenstates.
  * **Item 3 (cubic sensitivity).** The flat LP is informative here. In every optimum, cubic sensitivity k = −2μ₂₆ + 2μ₃₄ comes from a
    **ψ₃ enhanced by 1/a²**: ψ₃″(0) = 730 / 2045 / 6383 / 1.5·10⁵ at a = 0.1 / 0.075 / 0.05 / 0.025. Its heavy action is
    **channel-odd** (σ_h) and is tolerated only because B₃ carries u′ ≤ a². So the cubic-sensitive part is **not positive term by term**.
    A functional that is positive harmonic by harmonic must get its sensitivity from harmonics that heavy states do not populate
    (|λ−h| ≤ 1). Whether a completely monotone smearing can put a nonzero, correctly signed weight there is the first thing to check
    (½ day, exact). If heavy states populate every harmonic that carries t₄, α_T, the bound is empty.
* **What both routes share.** Both reduce LF3 to an exact, Δ-free inequality family plus a measure-concentration tail. The Gamma tail
  bounds are reproduced (§3.4). The compressed light side, i.e. this note's LF1/LF4 numbers, is common to both.

---

## 1. Method: a general low-frequency certification scheme (what a functional must satisfy, and why it suffices)

A cubic functional ω = ∫χ₂Ĉ′_EE + ∫χ₃Ĉ₃ (physical, de-shifted weights χ_k) is admissible for the scheme if:

* **(H1) Low-frequency target.** The target weights ψ_k (what the de-shifted leading-Regge images reproduce) have support in
  p = ν/Λ ≤ a. They are C⁴ at the edge ((1−p/a)^r, r ≥ 3, so edge oscillations at large b̂ are O(b̂^{−r−1}), below the b̂⁻³ floor).
  Also ψ₃(0) = 0, and −2ψ₂″(0) > (1+δ)a²|ψ₃″(0)| (the large-b̂ coefficient).
* **(H2) Flat margin.** 𝔉 ≥ δ(|B₂|+|B₃|) on every state with u′ ≤ a², certified (interval or exact on (u′, J) cells, plus the
  scaling limit and Lemma-M asymptotics), with δ ≥ 0.15. The **dilation identity** (§2.1) makes any certified flat functional at gap
  M′ certified for this problem.
* **(H3) Analytic physical weights with a certified e^{πν}-tail.** χ_k entire on |Im ν| ≤ 1 (apart from singularities cancelled by
  zeros of R_∓). The de-shift identity holds exactly (or with an explicit ε, carried into the LR image), and
  `T(ν_t) := ∫_{ν_t}^∞ |χ|(|R₋|+|R₊|)e^{πν}dν ≤ η_tail`, with core cut **ν_t ≤ ξ_core·Δ_g**. Compact support is excluded by Paley–Wiener
  (HK_CLOSE Lemma PW [P]), so this is the only legitimate form.
* **(H4) Leading-Regge floor.** The **realised** LR image (of χ, not of ψ) satisfies M_main ⪰ c·N_s on every heavy state and every
  coupling vector, certified (T0/LR∞/ROUTEA machinery). The null directions of the 2×2 helicity matrix must be treated explicitly:
  the floor is on N_s, the component norm.
* **(H5) Core remainder.** For ν ≤ ν_t, a **relative** bound `|v†E_core v| ≤ η_core·N_s`, uniform in (Δ, J, v), on whichever exact
  representation is used (§3).
* **(H6) Tail remainder.** `|v†E_tail v| ≤ η_tail·N_s`, where E_tail is the part of the action from ν > ν_t, bounded by (H3) times an
  absolute bound on |Ĉ_ν| in heavy normalisation, **uniform in J**.

Then **M_full ⪰ (c − η_core − η_tail)·N_s ≻ 0** if η_core + η_tail < c. That is per-state positivity, and the bound follows from the
light side of the same ω (LF4).

**Why it suffices [D].** The functional's value on a heavy state is linear in ω's weights. Split into main + core + tail. Each piece is
bounded against the same N_s, and N_s > 0 for every non-zero coupling vector. Every assumption above is a statement about ω or about
one block. None is a statement about the CFT, apart from A1/A2 (spectrum and Regge boundedness, already in the theorem) and the
exchange of sum and integral, H_K. H_K holds for (H3)-type weights by HK_CLOSE Prop. K_c, given K1′ and K2ν.

**What the scheme does not supply.** (H5) and (H6) are where every route either fails or remains open (§3). (H1)–(H4) are routine.

---

## 2. Flat side (LF1, LF2) [P + DV + NE]

### 2.1 Dilation identity [P] (reviewer's identity, re-derived)
* Set ψ_{2,a}(p) = a⁻¹ψ₂(p/a) and ψ_{3,a}(p) = aψ₃(p/a) on [0, a]. Substituting p = aq in
  `B₂ = ∫ψ₂(2−up²)D_h(1−2up²)` and `B₃ = σ_h u∫ψ₃D_h` gives:
  * B₂^{(a)}(u, J) = B₂(a²u, J) and B₃^{(a)}(u, J) = B₃(a²u, J);
  * k_a = a⁶k and L_a(Λ) = a⁻²L(aΛ).
* Hence **τ_a(Λ) = τ(aΛ)**: the compressed NF02 is NF02 at gap Λ′ = aΛ, tested only on m ≥ Λ = Λ′/a.
* **Certification transfers.** NF02's certificate covers all m ≥ M′ except (+−, J ≥ 400, M′ ≤ m < 1.633M′). For a ≤ 1/1.633 = 0.612
  that region is below Λ. So the **dilated NF02 is flat-positive with certified margin 3/20 on every heavy state** [certified + P].
* **Numbers** (`ln1_dilation_tau.py`, exact rationals + 50-digit L_A; K = 1): τ = 2.89126·10⁻⁸ at a = 1 (reproduces the headline), and:

  | Δ_g | a = 0.025 | 0.05 | 0.075 | 0.1 | 0.107 | 0.118 | 0.163 | 0.35 |
  |---|---|---|---|---|---|---|---|---|
  | 1000 | 3.84e-2 | 2.90e-3 | 6.26e-4 | 2.10e-4 | 1.62e-4 | 1.12e-4 | 3.24e-5 | 1.70e-6 |
  | 2000 | 2.89e-3 | 2.09e-4 | 4.44e-5 | 1.47e-5 | 1.13e-5 | 7.78e-6 | 2.24e-6 | 1.16e-7 |
  | 4000 | 2.09e-4 | 1.47e-5 | 3.08e-6 | 1.01e-6 | 7.78e-7 | 5.33e-7 | 1.52e-7 | 7.75e-9 |

  These reproduce the reviewer's 2.10e-4 / 6.26e-4 / 2.90e-3 / 3.84e-2.

### 2.2 Re-optimised low-ν LP [NE] (`ln2_lp.py`, `ln3_check.py`, `run_lp.sh`)
* **Set-up.** The same class and constraints as `nf_lp.py` (ψ_k = (1−p)⁵p·Legendre₈, cutting planes for 𝔉 ≥ δN, scaling limit,
  1/b̂³ coefficient), with three changes:
  * states only for u′ ∈ (0, a²];
  * objective = the **exact** 16L_A(Λ′), not the asymptotic `A log + B` (that form is exploitable at Λ′ = 25–100 and gave an
    unbounded LP);
  * optional matter rows for m_ℓ ≤ Λ. They change nothing.
* **Results** (δ = 0.2; τ from the exact L_A with the 9/(4Λ′²) regulator, then K(L) from §4):

  | a | 16L_A(Λ′) | K(L) | **τ (Δ_g = 1000)** | dilated NF02 | gain | δ = 0.1: τ | dense min 𝔉/N | ψ₃″(0) |
  |---|---|---|---|---|---|---|---|---|
  | 0.025 | 11.0 | 0.885 | **2.18e-3** | 3.84e-2 | 18 | 1.45e-3 | 0.2000 | 1.5e5 |
  | 0.05 | 15.9 | 0.973 | **1.56e-4** | 2.90e-3 | 19 | 1.09e-4 | 0.1995 | 6.4e3 |
  | 0.075 | 18.8 | 0.990 | **3.33e-5** | 6.26e-4 | 19 | 2.34e-5 | 0.2000 | 2.0e3 |
  | 0.1 | 21.0 | 0.995 | **1.11e-5** | 2.10e-4 | 19 | 8.2e-6 | 0.1991 | 730 |
  | 1 (control) | 10616 | — | 2.49e-8 | 2.89e-8 | 1.2 | — | 0.2000 | — |

* **Mechanism.** B₃ carries an explicit u′ ≤ a². So ψ₃ can grow like a⁻² before its channel-odd heavy action competes with B₂, and
  it raises k. The scaling τ ∝ (aΛ)⁻⁴ is unchanged; the constant improves ×19. **Rule of thumb: τ_LP(a, Λ) ≈ 1.1·10⁻⁵ (0.1/a)⁴
  (998.5/Λ)⁴ × (log factor).**
* **Dense check** (off-grid u′, J ≤ 6000, both helicities). Worst cases: 0.1991 at a = 0.1 (u′ = 0.995a², ++, J = 74) and 0.1995 at
  a = 0.05 (+−, J = 209). Scaling-limit B₂ > 0 to b̂′ = 2000; large-b̂ coefficient 3.1–46 > 0. **Not a certificate.**
* **LF2 margin.** Certified 3/20 (dilated NF02); 0.199 [NE] (LP). Optimising the margin first: δ = 0.2 costs ×1.35 in τ against δ = 0.1.

### 2.3 Large impact parameter in the flat/LR layer [P]
* Under dilation, b̂ ↦ ab̂. The cancellation regime (value ∝ ψ″(0)/b̂³ from oscillatory cancellation) starts at **b̂ ~ 1/a**
  (≈ 10–40) instead of b̂ ~ 1.
* The LR layer is still fine: LR∞-type asymptotics hold with the 1/b̂³ coefficient positive (H1). For a sign-free per-ν error budget,
  the tolerance at given b̂ shrinks by (ab̂)⁻³ relative to NF02 at b̂. The low-ν functional is **worse**, not better, for O2 (§3.3).

---

## 3. AdS remainder (LF3)

### 3.1 What Prop. D′ actually controls (Question 1) [read from RD_CLOSE §2.3, GAP_UNIFORM §4.6, PCYC_CLOSE]
* **Statement.**
  `∫_{D′}Φ_κ = Σ_{k≤K} a^k I_k^{engine} + O(a^{K+1}) + O(Δ^{N_κ}e^{−cΔ})`, with c = 2 artanh(r/2) and split radius r ≤ 1/3
  (c(1/3) = ln(7/5) = 0.3365).
  * **[D]** for κ ∈ {scalar, EE}, **J = 0 kinds** (GAP_UNIFORM §4.5: "[D] for the J = 0 kinds").
  * One **[D-sketch]** sub-step: the convergence radius of the x₀-dependent slow factors on the polydisc (M4(a)).
  * D′ = the physical chain: [P] + [DV] (PCYC_CLOSE).
* **It controls only the far region** (raw radius > r) as an **absolute** bound, |B|^{Δ/2} ≤ e^{−2Δ[artanh(s_x/2)+artanh(s_y/2)]} [P].
  * The **O(a^{K+1})** truncation of the near region has **no explicit constant**. A certified version is M4(b), never done.
  * The engine series is a **fixed-ν** series. At ν ∝ Δ it must be resummed in ξ (GAP_UNIFORM §4.6).
  * **So Prop. D′ does not bound E.** E lives in the near region. Prop. D′ only says that the near region is all there is, up to
    O(Δ^N e^{−cΔ}).
* **The ξ limit.** The signed value is e^{−πν} below the absolute scale, so the far error relative to the value is
  Δ^N e^{−(c−πξ)Δ}:
  * ξ_crit = c/π = **0.107** (r = 1/3);
  * the polydisc analyticity used for the slow factors (|1 ± s + q| ≥ 1 − 2r − 2r² > 0) allows r < (√3−1)/2 = 0.366, so
    ξ_crit ≤ **0.118** [D];
  * the measured c(½) = 0.511 would give 0.163, but it is not rigorous (the r = ½ polydisc is not analytic);
  * r → 1 makes the "far" region empty, but then the whole problem is the near region.
  * **True limit of the D′ route: ξ_max ≈ 0.107–0.118.** The reviewer's "ξ_crit → ∞ as r → 2" does not apply: s < 1 on D′, and the
    analyticity radius stops at 0.366.
* **At Δ = 1000 the prefactor matters.** N_κ ≈ p_EE ≈ 7–8 (RD_CLOSE §2.4 fits) gives (0.3365 − πξ)Δ ≥ 8 ln Δ + ln(1/η), i.e.
  **ξ ≤ 0.088** at η = 0.1 [D, given the fitted N].
* **The prefactor depends on J [D-sketch].** For spin J the slow factor contains P_J(t), with boost growth e^{J|η_x−η_y|}
  (BELOW_FOLD §3.2, [NE]). In the far region one light-cone component is O(r) and its partner O(a), so P_J ~ (rΔ)^{J/2}:
  N_κ(J) ≈ 8 + J/2.
  * At Δ = 1000 and core ξ ≤ 0.08 the far factor e^{−85} covers **J ≲ 8–10** only.
  * At ξ ≤ 0.025 it covers J ≲ 58.
  * Uniformity in J (O3) fails structurally on D′. That is BELOW_FOLD's failure mechanism for J > 0 on real sections.

### 3.2 Contour route: what low ν removes and what it does not
| obstruction (UNCOND_INVENTORY) | with support ξ ≤ a ≤ 0.1 |
|---|---|
| O1 contour theorem (CT) | **removed** for the far region, J = 0 kinds (Prop. D′); J > 0 still needs the [D-sketch] boost estimate made [D] |
| O4 light-cone fold | **removed**: x = ν/m ≤ a + tail. Tail weights ≤ e^{−πν}, so the fold region sees only the absolute D′ bound (which is J-dependent, §3.1) |
| near-region evaluation of E at Δ = 1000 | **not removed, and worse.** GAP_MORSE: ξ = 0.1 already binds every chart-A certificate (r_c = 0.023–0.025 against Gaussian width 0.032 at Δ = 1000), and ξ < 0.1 "would be worse" (the atom divisors sit at raw distance ≈ ξ/2) |
| O3 J ~ √Δ uniformity | **not removed** (§3.1: N(J) ≈ 8 + J/2) |
| O2 large-b̂ cancellation | **not removed; it starts earlier** (b̂ ~ 1/a, §2.3) |

* **The only new near-region idea [D-sketch, not attempted].** In the Regge-scaled chart (M4(c)) the atom factors are exact
  Gamma/Beta integrals. Then:
  * the large parameter becomes ν/2 (Gamma saddle at X ≈ iν/2, relative width (ν/2)^{−½});
  * the residual coupling exponent is ≈ (ξ²Δ/4)(X/X\*)², relative size ξ/2;
  * the slow-factor radius in X/X\* units is ≈ (1/3)/(ξ/2) ≈ 7.
* That suggests an optimal-truncation floor of order e^{−cν} rather than e^{−Δρ²/2}. It would make the core certifiable for
  **ν ≳ 20–30** at Δ = 1000. Below that, the fixed-ν engine with certified truncation (M4(b)) is needed. Both are 10–20 ah, risk
  M–H, and only for J ≲ 10. This is a T1-type programme, not a Thursday item.

### 3.3 The relative target and where it fails
* Target: |v†Ev| ≤ η·v†M_main v with M_main ⪰ c·N_s (H4–H6).
* **Small J, b̂ ≲ 1/a.** The low-ν functional has little cancellation: κ = Σ|contributions|/value = O(1). So a per-ν relative error
  ≲ 0.1 would suffice. Contrast NF02's rel\* = 5.6·10⁻⁵ at J = 8. **This is the genuine advantage of low ν**: tolerances become O(1).
  What is still missing is a certified evaluation (§3.2).
* **b̂ ≳ 1/a.** The value ∝ (ab̂)⁻³ comes from cancellation. A per-ν error bar degrades the tolerance by (ab̂)³. No sign-free bound can
  work. This is O2, unchanged; only signed (complex-ν or exact-representation) methods help.
* **Null directions.** The 2×2 helicity matrix of M_main is positive definite with the N_s floor (flat certificate on |B₂|+|B₃|),
  so there are no exact null directions in the main term. The relative bound must hold per entry against N_s (Gershgorin), as in
  Theorem W. Its factor 1.5–3 (INDEP_CHECK2) eats into η.

### 3.4 Exact-representation routes (reviewer input, checked where marked)
* **Analytic de-shift weight (k = 2, EE)** (`ln4_series.py`, `ln5_deshift.py`).
  * R₋(−i) = 18/25 and R₋′/R₋(−i) = 5i/9 [verified, 40 digits].
  * **37/36 is exactly the ν²·log ν cancellation** [verified]. E(0) = 0, E′(0) = κ, and the log terms start at ν⁴ [verified].
  * The Gaussian part is even in ν (Hermitian symmetry), so **E′(0) = κ exactly**, whatever the fit.
  * Both prefactors (double and single pole at ±i/2) are admissible: R₋ has a **double** zero at −i/2, and the single-pole version
    is what the scalar kernel R₋/c(z+i) needs.
  * **Strip analyticity:** w = z²+1 is never real-negative in −1 ≤ Im z ≤ 0, except at w = 0 (z = −i). That is the integrable
    endpoint of the shifted contour (w log w → 0), so the contour move is legitimate [D].
  * **Fit** (ridge, s = 0.9, spacing 0.27):
    * L = 0.05Λ: 208 Gaussians, ∫|E−T|/∫|T| = 2.2·10⁻⁴.
    * L = 0.1Λ: 393 Gaussians, 6.7·10⁻⁵.
    * S-transform error over the LR floor 2κe^{−β/2}/(πβ): 5·10⁻⁵ / 3·10⁻⁴ / 7·10⁻⁴ / 2·10⁻³ / 7·10⁻³ at β = 0.1 / 0.5 / 1 / 2 / 4,
      matching the reviewer.
    * It is **0.21 / 0.046 / 0.018 at β = 8 / 16 / 32**, against the reviewer's 2.7 % / 6.6 % / 1.5 %. My unconstrained fit is
      worse at β = 8. The difference is fit regularisation (r̂ at t ~ 1/s² ≈ 1.2–8), not quadrature (ny = 400 and 1600 agree).
    * Reviewer's tail-lemma constants for my fit: A = ∫|r| = 0.11 / 0.07 and B = |r‴(0)| + ∫|r⁗| = 795 / 454. B is large because
      of the narrow Gaussians. A smoother, derivative-penalised fit is needed before the β⁻² lemma is useful.
  * **e^{πν}-weighted tail:**
    * L = 0.05Λ: ∫_{80}^∞|χ|(|R₋|+|R₊|)e^{πν} = 10^{−47.6} (reviewer: < 2·10⁻³³; consistent).
    * **L = 0.1Λ fails**: the tail from 130 is 10^{+15}. Gaussians centred up to L+6 need ν_t − μ_max ≳ s√(2πν_t) ≈ 26.
  * **Rule:** ν_t ≈ L + 6 + 0.9√(2πν_t). With ν_t ≤ ξ_core Δ_g, ξ_core ≤ 0.08–0.088 gives **L ≤ 0.05–0.06Λ** at Δ_g = 1000 [DV].
* **Plane-wave (Bernstein / spin-Gamma) route.** The heavy factor is a positive Riesz average over plane waves (LARGESPIN_L L2: [P]
  scalar, **[D-sketch] for ⟨TTO⟩⟨OTT⟩**, which needs J-independent polynomial insertions with index shifts r ≤ 4).
  * **Gamma factorisation and Markov bounds re-checked in exact rationals** (`ln6_measure.py`, λ = (τ−1)/2, h = 7), τ = 1000, J = 0,
    r = 4: P(M < 600) ≤ 1.3·10⁻⁴⁴, P(M < 700) ≤ 2.1·10⁻²², **P(M < 800) ≤ 5.7·10⁻⁹**, P(M < 900) ≤ 3·10⁻². These are smaller at larger
    J (7·10⁻⁷² at J = 489) and at τ = 2000 (10⁻²⁶⁸). They agree with the reviewer's 2e-44 / 3e-22 / 6e-9.
  * **Synergy [D].** AG_PROOF's N2 obstruction (14–35 % of the plane-wave measure lies beyond its own edge, |p| < ν) disappears for
    low-ν weights, since ν ≤ aΛ ≪ 800 ≤ M. That is what makes the plane-wave route usable at all.
  * **Decisive open inequality.** 𝓜(p) ⪰ c·f(β(p)) for M ≥ 800, with C/c < 3·10⁶ below. 𝓜 is the functional's plane-wave image
    (P2_CONTRACT engine), a Δ-free function of (M, β) per helicity entry.
  * **Heuristic estimate [M].** N-L1 says the kernel corrections in 𝓜 cancel the Bernstein smearing, whose size is
    ≈ 0.55ν²/M (fit to LARGESPIN_L §2.6: 50 % at ν = 30, τ = 1000). So 𝓜 deviates from the LR image by relative O(ν²/M). At
    ν ~ aΛ, M ≈ 1000: ≈ 5.5 / 1.4 / 0.34 for a = 0.1 / 0.05 / 0.025, against a 0.2 margin.
  * **Pointwise positivity is therefore plausible only for a ≲ 0.02–0.03 at Λ = 1000** (τ ≈ 2–4·10⁻³), and for a ∝ Λ^{−1/2} in
    general, i.e. τ ∝ Λ⁻². The first test is cheap (the image engine exists): compute 𝓜(p) for the LP/dilated weights on an
    (M, β) grid.
  * **Matrix structure.** LARGESPIN_L/P2_SPINNING give a positive scalar measure times **signed** polynomial insertions. No
    matrix-valued positive average is established [I]. That is the second obligation.
* **Angular route.** See §0.

---

## 4. Physical bound (LF4) [DV + estimate]
* **K(L) re-checked** (`ln7_K.py`, exact rationals):
  * NF02: 0.99996752310 at L = 998.5 (reviewer 0.99996752314; the last digits differ by k₀ rounding), 0.99674323 (0.1Λ),
    0.98686374 (0.05Λ), 0.9457 (0.025Λ).
  * LP: 0.995 / 0.990 / 0.973 / 0.885 at a = 0.1 / 0.075 / 0.05 / 0.025. Included in §2.2's τ.
* **Light side needs** for (H3) weights:
  * L_A with the realised weight. The pole term is exact; the ε-de-shift error enters as ε·16L_A.
  * The O(c₁Λ⁻²) term becomes O(c₁Λ′⁻²), i.e. ×a⁻² relative. It is still ≲ 10⁻² of the bound at Λ′ = 50, but it must be recomputed
    (x/|c₁| scales as Λ′⁻²).
  * The K-correction (H17 remainder, δ₃) must be redone at Λ′. The K_CERT machinery is reusable.
  * None of these changes the order of magnitude [M].
* **τ as a function of the provable ξ_max** (LP constant; support edge a = ξ_max − tail overhead ≈ ξ_max − 0.026 at Λ = 1000):

  | provable ξ_max | source | usable a | τ at Δ_g = 1000 | Δ_g for τ = 10⁻⁶ (same a) |
  |---|---|---|---|---|
  | 0.088 | Prop. D′, r = 1/3, Δ^8 prefactor, J = 0 | ≈ 0.06 | **≈ 8·10⁻⁵** | ≈ 3000 |
  | 0.107 | Prop. D′ without prefactor | ≈ 0.08 | ≈ 2.7·10⁻⁵ | ≈ 2300 |
  | 0.118 | polydisc limit r = 0.366 | ≈ 0.09 | ≈ 1.7·10⁻⁵ | ≈ 2000 |
  | 0.163 | measured c(½) (not rigorous) | ≈ 0.14 | ≈ 3·10⁻⁶ | ≈ 1300 |
  | 0.2 | needed for 10⁻⁶ at 1000 | ≈ 0.18 | ≈ 1·10⁻⁶ | 1000 |
  | — | plane-wave route heuristic, a ∝ Λ^{−½} | 0.025 | ≈ 2·10⁻³ | ≈ 5·10⁴ (τ ∝ Λ⁻²) |

  Dilated NF02 (certified flat) values are ×19 larger (§2.1).

---

## 5. Wed–Thu plan (if pursued)
1. **½ day, decisive for the angular route.** The cubic-sensitivity check (ANGULAR_ROUTE item 3): can a completely monotone smearing
   that is positive on every harmonic |λ−h| ≥ 2 keep k > 0? Exact, using the CDH formulas. If not, drop the angular route.
2. **½ day, decisive for the plane-wave route.** Compute 𝓜(p) (P2_CONTRACT engine) for the a = 0.025 / 0.05 LP weights on an
   (M ∈ [800, 4000], β ∈ [0, 8]) grid, scalar insertion first. This tests the O(ν²/M) heuristic of §3.4.
3. **In parallel, cheap and banked either way.**
   * Certify the a = 0.05 LP functional flat (`nf_cert` adaptation, 3–6 ah), or adopt the dilated NF02 (already certified).
   * Refit the analytic de-shift weight with derivative penalties to bring B down and make the β ≥ 8 error ≤ 3 %.
   * Recompute the light side at Λ′ (K, c₁ term).
4. **Not recommended:** the contour route at Δ_g = 1000 (§3.2). It needs M4(a)+(b)+(c), and still covers only J ≲ 10.

**Realistic Thursday outcome.** A certified flat functional with τ ≈ 1.6·10⁻⁴ (LP) or 2.9·10⁻³ (dilated NF02) at Δ_g = 1000. Its
leading-Regge positivity is plausible with the T0 machinery. **Full-action positivity stays conditional** unless item 1 or 2 succeeds
and its matrix/transition obligations close. I put that below 25 % by Thursday.

---

## 6. Reproduce (`ads4/code/lownu/`; `ulimit -v 3000000; OMP_NUM_THREADS=1`)
```
python3 ln1_dilation_tau.py                 # 5 s   -> cache/ln1_tau.txt (dilated NF02, exact moments)
./run_lp.sh                                 # 12 min -> cache/run_lp.txt, run_check.txt (LP a in {.025,.05,.075,.1}, delta .2/.1, matter off/on; dense checks)
python3 ln2_lp.py 1.0 5 8 0.2 0             # control at a = 1
python3 ln4_series.py                       # 2 s   -> cache/ln4_series.txt (R_-(-i), 37/36, E(0), E'(0), log terms)
python3 ln5_deshift.py 49.925 2 0.27 0.9 1e-8   # 2 s; also 99.85 -> cache/ln5_*.json (analytic de-shift fit, S errors, tails, A/B)
python3 ln6_measure.py                      # 1 s   -> cache/ln6_measure.txt (Markov bounds, exact rationals)
python3 ln7_K.py                            # 1 s   -> cache/ln7_K.txt (K(L) for NF02 and the LP functionals)
```
