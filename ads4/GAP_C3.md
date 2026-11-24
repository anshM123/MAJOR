# GAP_C3: explicit Laplace remainder constants C₃(ξ, J) for the chart-A crystal integrals

*24 Nov 2026, gap-C₃ agent. Unreviewed, not committed. No existing file was edited; `gapnear/`, `gapuni/`, `spinfunc/` were only read
or imported.*
* *New code:* `ads4/code/gapc3/` (`run_all.sh`, ≈ 1.5 CPU-h in total; at most 2 processes; `ulimit -v 3000000`; `OMP_NUM_THREADS=1`).
  Every script installs `gc_guard.py`, which refuses writes outside `gapc3/cache/`. That directory's `.gitignore` ignores everything
  but itself.
* *Read-only imports:* `spinfunc/l1_crystal, l1_tps, l5_spin, l5_shadow_matrix, l6_spin_action, l9_points, w0_ee_common,
  p3_odd_action, n1_pw` (fresh in-memory copies on a new backend; the only source substitutions are `np.zeros(..., complex)` →
  `object` and one `assert T.DEG in (6, 8)`).
* *Numbers taken from other documents:* the tolerances `G_cert(J)` of GAP_UNIFORM §3.2 and its Lemma T.

Grades as in GAP_UNIFORM: **[P]** proved; **[certified]** computer-assisted with interval (Arb ball) arithmetic; **[D]** derived;
**[D-sketch]**; **[DV]** derived + float check; **[NE]** numerical evidence; **[M]** model/convention; **[A]** assumption; **[I]** open.
Notation: `ε = Δ^{-1/2}`, `ξ = ν/Δ`, `Φ_ξ = (iξ/2) log A + ½ log B` (GAP_UNIFORM Prop. S), `v` = Morse coordinates of Φ_ξ at the
crystal (`z = z* + L_Φ v`, `L_Φ = √Δ · L`, L the l1_crystal Morse matrix; Δ-independent at fixed ξ), `P_ρ = {v ∈ C⁷ : |v_i| ≤ ρ}`.

---

## 0. Verdict

**C₃ is now explicit and certified, but it is ~10⁸–10¹² too large to be used at Δ_g = 1000. What it buys is an explicit,
certified "eventually" (Δ ≳ 10⁶) for the nodes computed, not positivity from 1000 on. The obstruction is structural: it is the
sup-norm Laplace method itself, not the evaluator. Even ideal (exact) sup norms cannot go below Δ ≈ 3.5·10⁵ with this lemma.**

| # | statement | grade | coverage |
|---|---|---|---|
| **TM** | **Pointwise complex evaluator with interval arithmetic.** `gc_ball.py` is an Arb **Taylor-model** backend with the `l1_tps` interface: polynomial in the 7 Morse variables (degree D = 3–8, ball coefficients = the exact Taylor coefficients at the crystal), plus a rigorous remainder ball valid on all of `P_ρ`. Logs and powers are continued from the crystal (the code's branch convention), and every log/power/reciprocal checks `sup_P |f/f* − 1| < 1`, which certifies analyticity on P_ρ. The project's builders run on it **unchanged**: scalar (l1_crystal), EE (w0 even/odd), C2/C3H (l9), odd J (p3). `gc_scalar.py` is a log-factored rewrite of the scalar integrand; it agrees with l1_crystal to 3·10⁻¹⁶ (phase) and 2·10⁻¹⁴ (amplitude) at random complex points [DV]. | **[certified]** (enclosures) | all kinds; certifies polydiscs up to ρ ≈ 0.02–0.03 |
| **L3** | **Explicit Laplace remainder (Lemma L3).** Take an analytic phase on P_ρ with an exact critical point at v = 0 and quadratic part `−½vᵀ(1+E)v`. Bound the cubic remainder by Taylor-model data, `|f₃(v)| ≤ Σ_{k=3}^{D} B_k(‖v‖_∞/ρ)^k + rem_Φ (‖v‖_∞/ρ)^{D+1}`, and the amplitude by `|h| ≤ M_h(τ)` on `P_{τρ}`. Then, for every Δ with `2ε < τρ`, the K₂ truncation (Laplace orders ε⁰, ε², ε⁴) of the integral over the real Morse box `‖v‖_∞ ≤ τρ/2` has relative-to-envelope error at most `rel(Δ) = (2π)^{-7/2}[2M_h(ε/τρ)⁶·7·15/(8b³)(π/b)^{7/2} + 6M_h e^{−bU²/2}U^{−5}·7Γ(3)c^{−3}(π/c)³]`, with `b = ½(1−7E) − κ(τ)`, `κ = (Σ_k B_kτ^{k−2} + rem_Φτ^{D−1})/ρ²`, `U = τρ/(2ε)`, `c = b/2`. The bound decreases in Δ, so `C₃ := rel(Δ₁)Δ₁³` is valid for all Δ ≥ Δ₁. | **[P]** (proof §2) | any node with certified data |
| **C3-S** | **Certified C₃ for the scalar k = 2 integrand with spin contraction Q_J** (this integrand shares Φ_ξ, and hence all the phase data, with the EE/C3H kinds). Nodes: J ∈ {0, 2, 4, 5, 6, 7, 8} × ξ ∈ {0.1, 0.3, 0.5, 0.7, 0.85}, 35 nodes, all certified. Results: **C₃ = 10¹⁴–10¹⁸**, and the Lemma-T accuracy `rel ≤ rel*_J` holds for every **Δ ≥ Δ₁ = 4.7·10⁵ – 10⁷** (table §4.1). The needed size at Δ₀ = 1000 is `C₃ ≲ rel*_J·10⁹ = 5·10⁴ (J = 8) … 3.5·10⁷ (J = 0)`. | **[certified]** (data) + [P] (L3) | the 35 (J, ξ) nodes, the box integral |
| **C3-G** | **Graviton amplitudes, certified sup norms** (EE k = 2 via the w0 builder; C3H k = 3 via the l9 builder, run on the TM backend). Done for single-structure J = 0, 2 at ξ = 0.3, 0.5 and 0.7, with ⟨TTO⟩ coefficient vectors frozen at Δ = 1002. EE J = 2: Δ₁ = 1.0·10⁶ (ξ = 0.5), 1.2·10⁶ (0.7), 1.7·10⁶ (0.3); EE J = 0: 1.0·10⁶ (0.5); C3H J = 2: 1.5·10⁶ (0.5), 1.8·10⁶ (0.7); C₃ = 10¹⁶ (§4.3). | **[certified]** for the frozen-coefficient amplitude | J = 0, 2 even; J ≥ 4 and odd J not run |
| **B** | **The sup-norm barrier.** Rerun Lemma L3 on *sampled* (ideal, uncertified) sup norms. It gives Δ₁^ideal = 3.5·10⁵ (J = 0), 4.7·10⁵ (J = 7), 4.2·10⁵ (EE J = 2) at ξ = 0.5, and enlarging ρ from 0.03 to 0.1 makes this *worse*. The reason: the genuine cubic size `S₃ = sup|f₃|/ρ³ ≈ 10–12` in Morse units (certified: 11.0 at D = 8 against 9.8 sampled) caps the complex Cauchy radius at `τρ ≲ 0.4/S₃ ≈ 0.035`. A remainder of this type needs `ε ≪ 0.035`. **So no sup-norm Laplace bound with a complex-ε Cauchy step can certify the K₂ truncation at Δ = 1000.** The actual next Laplace term (l9 K3 records) is ~2·10⁻⁵ relative at Δ = 1000: the bound is pessimistic by ~10¹⁰. | **[D]** (mechanism) + **[NE]** (sampled sup norms) | — |
| **E′** | **Explicit eventual positivity.** At each certified node, GAP_UNIFORM's non-effective Λ\* (Theorem E) becomes explicit, `Δ_g* ≤ Δ₁(J, ξ) ≈ 0.5–5·10⁶` (10⁷ at J = 8, ξ = 0.1). This holds given (i) premise TAIL (the cycle outside the Morse box is `O(e^{−δΔ})`: GAP_UNIFORM Lemma R [D-sketch]); (ii) the chain propagation factor Π ≤ 3 [M] (§4.2); (iii) `sup_m |Γ_{K₂,m}|` bounded (B_F restricted to the explicit K₂ function, [NE]); (iv) continuity in ξ between the certified nodes [D-sketch]; (v) for the graviton kinds, (C3-G) only at J = 0, 2. | [P given (i)–(v)] | certified nodes |

**What this does not change.**
* **U7 of GAP_UNIFORM (J ≤ 7 positive at every Δ_g ≥ 1000) still rests on premise B_F on the whole range 1000 ≤ Δ_g ≲ 10⁶.** The
  certificate only replaces B_F by a certified statement for Δ_g ≳ 10⁶ (at the computed nodes).
* No new (J, μ) is covered at Δ_g = 1000. No window statement (ν < 50, x > 0.9) is changed.

**What remains, in order of value.**
1. **A sharper remainder lemma.** The cubic term in the exponent is what kills Lemma L3. The next step is a *certified holomorphic
   Morse chart*: Φ = Φ* − ½y·y exactly, by Taylor-model inversion with a Newton–Kantorovich domain. The cubic then moves into the
   amplitude (as a Jacobian) and the Cauchy radius becomes the analyticity radius ρ_an ≈ 0.1–0.15 [NE, §5], instead of 0.4/S₃.
   Rough estimate with ideal sup norms: Δ₁ ≈ 3·10³ (J = 0) to 1.3·10⁴ (J = 7) [NE, §5]. **Even then Δ = 1000 is out of reach**: for an
   asymptotic series the attainable accuracy at Δ is ~`e^{−Δρ_an²/2}` (optimal truncation), i.e. 7·10⁻³–2·10⁻⁵ at Δ = 1000. That is
   below rel*₇ = 4.2·10⁻⁴ only if ρ_an ≳ 0.13.
2. **The range 1000 ≤ Δ_g ≤ Δ₁ by other means.** Two options: K₂ rows at a geometric ladder of gaps, with a rigorous interpolation
   in 1/m (the K₂ step is m-independent in c₂, GAP_UNIFORM §4.1); or certified quadrature of the thimble integral (not attempted:
   7 real dimensions).
3. **ξ-intervals instead of ξ-nodes.** A direct ξ-ball in the TM inputs fails, because cancellations are lost and E grows to ~10³
   for width 0.025 (tested). The fix is to make δξ an 8th Taylor-model variable; the lemma's quantities are then read off per ξ
   (§4.4).
4. **Graviton kinds at J ≥ 4 and odd J** (two structures / parity-odd 1×1): the evaluator runs (`gc_node.py`). Still needed: the
   ⟨TTO⟩ coefficients as functions of Δ (rational; only float null vectors exist), the shadow matrix, and the helicity basis
   (§4.3).

---

## 1. What is integrated, and what "K₂" means

At fixed ξ the crystal integral is (GAP_UNIFORM Prop. S [D])
```
I(Δ) = c(Δ, ν) ∫ e^{Δ Φ_ξ(z)} g(z) d⁷z,   g = e^{G(z)} × prefactor(z)
```
* **Prefactor by kind.** Scalar: Q_J. EE: the w0 contraction `E⊗E · ⟨TTO⟩·β⁻¹·⟨ÕTT⟩`. C3H: `q·o_H`.
* **Morse coordinates.** In v the integrand is `e^{ΔΦ*} e^{−Δ v·v/2} e^{Δ f₃(v)} g(v)`, with f₃ the Taylor remainder of Φ_ξ from
  order 3.
* **Scaling.** With `v = εu`, the u-integrand is `e^{−u·u/2} F(ε, u)`, where `F(ε,u) = e^{f₃(εu)/ε²} g(εu)` is analytic in ε and
  `F(−ε, −u) = F(ε, u)`. So the ε-Taylor coefficients `F_k(u)` have parity `(−1)^k` in u, and the odd ones integrate to zero.
* **K₂.** The code's K₂ (l1_tps.expand_pref: S0, S2, S4; l9/w1 rows: M0, M2, M4) is the Gaussian expectation of F₀, F₂, F₄. That is,
  the Laplace expansion to relative O(ε⁴) = O(Δ⁻²). The truncation certified here is therefore **K₂, remainder O(Δ⁻³)**.
* **Not a Laplace error.** The code drops the Morse defect E (|E| ≤ 3–4·10⁻¹⁵, certified, §3). That changes each coefficient by
  O(E) relative, i.e. ≤ 10⁻¹⁰ c₂ at m ≤ 10⁴, and it is a Δ-independent relative error, which Lemma T's tolerance (also relative)
  absorbs at every gap.
* **Envelope normalisation.** Γ is measured in "c₂ × LO shape": `Γ = m²[Re G/(a·env) − σD_h]` (l9_analysis). So the relevant
  remainder is relative to the envelope `g_env = g(z*)/D_h(ξ)`, not to g(z*), which vanishes at zeros of P_J. For the scalar
  `g_env = g(z*)/P_J(cos θ*) = e^{G*}·|ZZWW*|^{J/2}`; for EE `g_env = |g(z*)|/|P_J(1−2ξ²)|` (§4.3). With `h = g/g_env`, a relative
  error `rel` of the box integral moves Γ by `m²·rel`, up to the propagation factor Π of §4.2.

## 2. Lemma L3 (explicit remainder) [P]

**Hypotheses** (all about Δ-independent objects at fixed ξ):
* (a) Φ_ξ and h are analytic on P_ρ. This is certified by the TM checks `|f_k/f_k* − 1| < 1` for every elementary factor
  f_k that enters a log or a power.
* (b) `∇Φ_ξ(z*) = 0` exactly. This is the CLS-D crystal [Q]; the certified enclosure is `|∇Φ| ≤ 4·10⁻²⁹` at every node.
* (c) Hessian `−(1+E′)` with `|E′_ij| ≤ E`.
* (d) Taylor-model data: degree-k sums `B_k = Σ_{|α|=k}|Φ_α|ρ^k` (k = 3..D, ball upper bounds) and `rem_Φ ≥ sup_P|Φ − P_D|`; the
  same for h.

**Step 1 (cubic bound).** `Φ − P_D` vanishes to order D+1 at 0, so the Schwarz lemma on the complex lines λ ↦ λw, w ∈ ∂P_ρ, gives
`|f₃(v)| ≤ Σ_{k=3}^{D} B_k s^k + rem_Φ s^{D+1}` with s = ‖v‖_∞/ρ. Each homogeneous degree-k part obeys `|Φ_k(v)| ≤ B_k s^k`.

**Step 2 (Cauchy in ε).** Fix real u ≠ 0 and τ ∈ (0, 1], and put `R = τρ/‖u‖_∞`.
* For `|ε′| = R` we have `ε′u ∈ P_{τρ}`, so `|f₃(ε′u)|/R² ≤ κ(τ)‖u‖_∞²` and `|F(ε′,u)| ≤ M_h(τ) e^{κ‖u‖_∞²}`.
* By Cauchy, `|F(ε,u) − Σ_{k≤5} ε^k F_k(u)| ≤ M(u)(ε/R)⁶/(1 − ε/R)` and `|F_k(u)| ≤ M(u)R^{−k}`.

**Step 3 (integration).**
* *Inside the box* `‖u‖_∞ ≤ U = τρ/(2ε)` (so ε/R ≤ ½): use `|e^{−½uᵀ(1+E′)u}| ≤ e^{−½(1−7E)|u|²}` and `‖u‖_∞² ≤ |u|²`. The
  moment `∫e^{−b|u|²}‖u‖_∞⁶ ≤ 7·(5!!)(2b)^{−3}(π/b)^{7/2}` gives the first term of rel.
* *Outside the box*, the truncated series is bounded by `6M_h e^{κ|u|²}(‖u‖_∞/U)⁵`. Its Gaussian integral over `‖u‖_∞ > U` is at most
  `e^{−bU²/2}U^{−5}·7·Γ(3)c^{−3}(π/c)³` with c = b/2.
* *Odd orders.* They vanish on R⁷: the weight is even and F_k has parity (−1)^k.
* *Normalisation.* Dividing by the Gaussian normalisation (2π)^{7/2} gives the bound on `|I_box/(Gauss·g_env) − Σ_{k≤2}c_kε^{2k}|`,
  where `c_k = (2π)^{−7/2}∫e^{−½uᵀ(1+E′)u}F_{2k}`. ∎

*Remarks.*
* Monotonicity: for fixed τ, both terms decrease in Δ (ε⁶; and `e^{−bU²/2}U^{−5}` with U ∝ √Δ). So `rel(Δ₁) ≤ target` implies the
  same for every Δ ≥ Δ₁: this is the "uniformity for free" of GAP_UNIFORM §4.3.
* The order-2 variant (remainder after the leading term) is the same proof with ε² and moment 7·1/(2b).
* *Where it is weak.* The factor `e^{κ‖u‖²}` forces `κ(τ) < ½`, i.e. τρ ≲ 1/(2S₃) at small τ. The complex rotation of ε turns the
  bounded oscillation `e^{iIm f₃}` into exponential growth `e^{|f₃|}`. This is the barrier B of §0.

## 3. The evaluator (TM backend) [certified]

* **Representation.** A Taylor model of degree D in the 7 Morse variables on P_ρ is (polynomial with Arb ball coefficients) + (one
  remainder ball).
* **Product.** `(P₁+I₁)(P₂+I₂)`: the truncated product, plus the high part `Σ_{k+l>D}B_k(P₁)B_l(P₂)`, plus `B(P₁)|I₂| + B(P₂)|I₁| +
  |I₁||I₂|`.
* **Composition** (`1/a`, `log a`, `a^p`, `e^a`) is around the crystal value a₀, with `t = a/a₀ − 1` and `q = sup_P|t| < 1`
  (else rejected). Remainder bounds:
  * `q^{D+1}/(1−q)` (reciprocal);
  * `q^{D+1}/((D+1)(1−q))` (log);
  * `|binom(p,D+1)|q^{D+1}(1−q)^{p−D−1}` (powers, integral form, real p ≤ D+1);
  * `e^{|d|}|d|^{D+1}/(D+1)!` (exp).
* **Linear operations** act on the remainder ball automatically: balls centred at 0 add their radii. So the builders' linear algebra
  (tensordot, einsum, exact DFT on the null circles with Arb roots of unity, exact κ_l) is unchanged.
* **Branches.** Log and powers are `log a₀` (principal at the crystal, as l1_tps and n1_pw) plus `Log(1+t)`: the continuation from
  the crystal, analytic on P_ρ when q < 1.
* **Why a log-factored rewrite for the scalar phase.** The project's builders compute composite quantities: A has a denominator
  ≈ −775 at ξ = 0.5, and the subtraction factor is a difference ≈ −0.048 of O(1) terms. Their TM ranges are inflated, so the
  original l1_crystal certifies only ρ ≲ 0.01 (q = 0.43 at ρ = 0.02 for the denominator of A, and log f_sub fails). `gc_scalar.py`
  writes `Φ − Φ*` as `(iξ/2)Σa_k Log(f_k/f_k*) + ½Σb_k Log(f_k/f_k*)` over 23 elementary factors (m_ij, c13, c10, c30, Nn, …). It
  writes `g/g*` as a product of `(f_k/f_k*)^{p_k}`, with f_sub as an explicit rational function, so no log is taken of something
  that can vanish without being a branch point.
* **Validation.** Against n1_pw/l1_crystal at 50 random complex points (|v| ~ 0.03, J = 4): 3·10⁻¹⁶ (Φ) and 2·10⁻¹⁴ (g) [DV].
* **What it reaches.** The binding factor has `|∇log f_k|₁ ≈ 19–22` in Morse units. Individual factors vary fast even though Φ does
  not: their gradients cancel at the crystal. That gives q ≈ 20ρ (ξ = 0.5): q = 0.43 at ρ = 0.02, and the hard limit q = 1 is at
  ρ ≈ 0.047. Usable ρ: 0.015–0.025 (D = 6), 0.025–0.03 (D = 8).
* **Per-node data** at ξ = 0.5, J = 0:
  * D = 6, ρ = 0.01 / 0.015 / 0.02: S₃ = 10.8 / 11.3 / 12.6; M_h = 1.15 / 1.26 / 1.64; rem_Φ = 9·10⁻⁸ / 1.6·10⁻⁶ / 1.3·10⁻⁵.
  * D = 8, ρ = 0.02 / 0.025 / 0.03: S₃ = 11.0 / 11.3 / 11.8; M_h = 1.36 / 1.82 / 6.2.
  * Every node: E ≤ 9·10⁻¹⁵. |∇Φ(z*)|: the phase coefficient iξ/2 is now built from the exact ball ν/Δ, the same number as the
    crystal's, so the enclosure is ~10⁻²⁹ and contains 0. An earlier version used ξ rounded to double, giving enclosures up to
    5·10⁻¹⁷: a genuine linear term γ·v, |γ| ≤ 5·10⁻¹⁷, whose effect `e^{√Δ γ·u}` is ≤ 1 + 10⁻⁶ for Δ ≤ 10²⁰.
* **Speed.** 5 s per scalar node at D = 6, 50 s at D = 8; 20 s per EE node at D = 5.

## 4. Results

### 4.1 Certified scalar nodes (`gc_run.py` → `cache/gc_cert_J*.json`; `gc_table.py` → `cache/gc_table.txt`)

* **Target.** For Λ ≫ Λ₀ Lemma T needs `m²Π·rel + sup|Γ_K₂| ≤ G_cert(Λ/Λ₀)²`, i.e. `rel ≤ rel*_J = 0.9·G_cert(J)/(Π μ² Λ₀²)`,
  with Π = 3 and μ = 1.3 (the largest tabulated m/Λ). This gives rel* = 3.5·10⁻² (J = 0), 1.1·10⁻² (2), 2.0·10⁻³ (4), 5.4·10⁻³
  (5), 6.4·10⁻⁴ (6), 4.2·10⁻⁴ (7), 5.6·10⁻⁵ (8).
* **Columns.** Δ₁ is the smallest Δ (40 points per decade) with `rel ≤ rel*`; it holds for all Δ ≥ Δ₁. C₃ = rel(Δ₁)Δ₁³, best ρ and
  degree. Δ₁(o2) is the order-2 variant (no K₂ terms used).

| J | ξ = 0.1 | 0.3 | 0.5 | 0.7 | 0.85 | C₃ range | Δ₁(o2) at ξ = 0.5 |
|---|---|---|---|---|---|---|---|
| 0 | 2.8e6 | 8.9e5 | **4.7e5** (D = 8) | 7.5e5 | 1.7e6 | 3.5e15 – 7.3e17 | 1.3e7 |
| 2 | 3.4e6 | 1.0e6 | 6.0e5 | 7.9e5 | 1.6e6 | 2.3e15 – 3.0e17 | 6.7e6 |
| 4 | 4.0e6 | 1.1e6 | 7.9e5 | 9.4e5 | 2.0e6 | 9.1e14 – 1.1e17 | 6.7e7 |
| 5 | 3.8e6 | 1.1e6 | 6.7e5 | 7.5e5 | 1.8e6 | 1.6e15 – 1.8e17 | 1.3e7 |
| 6 | 4.7e6 | 1.4e6 | 1.0e6 | 1.1e6 | 2.4e6 | 5.7e14 – 6.1e16 | 2.5e8 |
| 7 | 5.3e6 | 1.6e6 | 1.1e6 | 1.1e6 | 2.2e6 | 2.9e14 – 6.0e16 | 2.5e8 |
| 8 | 1.0e7 | 2.4e6 | 1.3e6 | 1.6e6 | 3.6e6 | 1.0e14 – 5.2e16 | 1.1e9 |

**Needed C₃ at Δ₀ = 1000** (`rel*_J·10⁹`): 3.5·10⁷ (J = 0), 1.1·10⁷ (2), 2.0·10⁶ (4), 5.4·10⁶ (5), 6.4·10⁵ (6), 4.2·10⁵ (7),
5.6·10⁴ (8). The certified constants exceed these by 10⁸ (J = 0, ξ = 0.5) to 10¹² (J = 8, ξ = 0.1).

*Data provenance.* The table merges two sets of records:
* the first batch (`cache/prefix/`, phase coefficient ξ rounded to double, |∇Φ(z*)| ≤ 5·10⁻¹⁷);
* the re-run with the exact-ξ phase (`cache/`), which overrides them.

The re-run was stopped at the time box after J = 0 and J = 8 (28 nodes). On every node it finished, q, S₃, M_h and rem_h
agree with the first batch to 4·10⁻¹⁰ relative. `run_all.sh` regenerates everything with the exact-ξ phase.

### 4.2 Propagation to Γ [M]
* `Γ_{k,h}(x) = m²[Re G_{k,h}/(a_h env) − σD_h]`, where G is the state/reference cumulant ratio `pred(state, ref, K)` (l9_analysis).
  A relative-to-envelope error enters three ways:
  * δ_s, in the state, moves Γ by m²δ_s·(1 + O(Γ/m²));
  * δ_r, in the reference (same ν, Δ_ref ≥ Δ, ξ_ref ≤ ξ), moves Γ by m²δ_r;
  * through the per-state constant a_h, fitted on x ≤ 0.25, ν ≥ 15, it moves Γ by m²|δa/a|·|D_h|.
* Π = 3 encodes `|δ_r|, |δa/a| ≤ |δ_s|`. That is [M]: it needs the certificate at the reference's ξ_ref and at the fit nodes,
  which are small-ξ nodes where Δ₁ is largest (column ξ = 0.1).
* The de-shift and the helicity rotation (J ≥ 4) are linear in the same crystal values and are not analysed here.

### 4.3 Graviton amplitudes (`gc_spin.py` → `cache/gc_spin.json`, `gc_report.py`)
* The EE (w0 even) and C3H (l9) builders run unchanged on the TM backend; the amplitude's TM is taken from the builder. The phase
  data come from the log-factored integrand at the same (ξ, ρ, D), since Φ_ξ is kind-independent.
* The builder TMs stay inside q < 1 up to ρ ≈ 0.015 (EE J = 2: q_amp = 0.19 / 0.31 at ρ = 0.01 / 0.015; ρ = 0.02 is rejected at
  ξ = 0.3).
* Results (D = 5, envelope `|g(z*)|/|P_J(1−2ξ²)|`):

| kind | J | ξ | ρ | M_h (env) | D_h(ξ) | Δ₁ | C₃ |
|---|---|---|---|---|---|---|---|
| EE | 0 | 0.5 | 0.015 | 1.46 | 1 | 1.0e6 | 2.9e16 |
| EE | 2 | 0.3 | 0.01 | 0.63 | 0.51 | 1.7e6 | 3.9e16 |
| EE | 2 | 0.5 | 0.015 | 0.50 | 0.125 | 1.0e6 | 1.0e16 |
| EE | 2 | 0.7 | 0.015 | 1.07 | 0.50 | 1.2e6 | 1.4e16 |
| C3H | 2 | 0.5 | 0.01 | 1.15 | 0.125 | 1.5e6 | 3.5e16 |
| C3H | 2 | 0.7 | 0.01 | 2.92 | 0.50 | 1.8e6 | 5.0e16 |
| C3H | 0 | 0.5 | 0.01 | 2.17 | 1 | 1.4e6 | 9.3e16 |

  (C3H J = 0 was stopped at the time box after ξ = 0.5. These spin records predate the exact-ξ fix of §3: their phase enclosure has
  |∇Φ(z*)| ≤ 5·10⁻¹⁷, from ξ rounded to double precision in the phase coefficient; see §3 for its effect.)
* **Caveat (Δ-dependence of g).** The ⟨TTO⟩ null vectors `sol(Δ)`, `solt(3−Δ)` (l5_spin.conserved) are float numerical null vectors,
  frozen here at Δ = 1002. The certificate is for that fixed amplitude. Uniformity in Δ needs `sol(Δ)` as rational functions with a
  ball enclosure over `1/Δ ∈ [0, 1/Δ₀]`; the TM backend accepts ball coefficients, so this is input work only. The shadow matrix
  β⁻¹ and the helicity basis (J ≥ 4) are outside the certificate.

### 4.4 ξ-continuity
* The certificate is per node. Chart-A nodes `x_n = p_n/μ` and `ξ = x(1 − 3/(2m) + …)` form a continuum.
* A ξ-interval passed as a ball fails: at width 0.025 the enclosure gives E ≈ 1.5·10³ and |∇Φ| ≈ 37, because the cancellations
  between factors are lost.
* The remedy is δξ as an 8th TM variable. The lemma's per-ξ quantities follow from Cauchy estimates in v of the remainder (bounded
  by `rem·C(k+6,6)` for the v-degree ≤ 2 parts). Not done. Between nodes, the table varies smoothly (Δ₁ within ×2 over
  ξ ∈ [0.3, 0.7]) [NE].

## 5. Why Δ = 1000 is out of reach for sup-norm bounds [D + NE] (`gc_ideal.py` → `cache/gc_ideal.json`)

* **Ideal sup norms.** `M₃(ρ) = sup|f₃|` and `M_h(ρ)` are estimated by sampling the torus |v_i| = ρ (3000 points plus local ascent)
  through n1_pw, for the scalar J = 0, 7 and EE J = 2 at ξ = 0.5. They are fed into Lemma L3 in Schwarz form:

| kind | J | ρ | M₃/ρ³ | M_h (env) | Δ₁^ideal | rel bound at Δ = 10⁴ |
|---|---|---|---|---|---|---|
| S | 0 | 0.03 / 0.05 / 0.10 | 9.8 / 10.4 / 11.7 | 1.5 / 2.0 / 3.7 | 3.5e5 / 4.2e5 / 5.6e5 | 3.9e5 / 6.9e5 / 3.4e6 |
| S | 7 | 0.03 / 0.05 / 0.10 | 9.8 / 10.4 / 11.7 | 0.28 / 0.41 / 1.4 | 4.7e5 / 5.6e5 / 9.4e5 | 7e4 / 1.4e5 / 1.3e6 |
| EE | 2 | 0.03 / 0.05 / 0.10 | 9.8 / 10.4 / 11.7 | 0.29 / 0.44 / 1.0 | 3.4e5 / 4.0e5 / 5.6e5 | 7.5e4 / 1.5e5 / 9.3e5 |

  (Lemma L3 does not apply at Δ = 1000: `2ε = 0.063 > τρ`.)
* **The certified evaluator is close to ideal.** The certified S₃ at D = 8 (11.0–11.8) is within 20 % of the sampled one, and the
  certified Δ₁ = 4.7·10⁵ (J = 0, ξ = 0.5) is within ×1.4 of Δ₁^ideal. **The barrier is the lemma, not the evaluator.**
* **Mechanism [D].** κ < ½ forces `τρ ≲ 1/(2S₃) ≈ 0.04` (S₃ ≈ 10 is a genuine property of Φ_ξ: the sampled sup norm agrees). The
  Cauchy step needs `2ε‖u‖_∞ < τρ` over the Gaussian bulk (‖u‖_∞ ≈ 3–4). Hence `Δ ≳ (8/0.04)² ≈ 4·10⁴` before any constant, and
  with the 7-dimensional moment factor (≈ 10⁴), Δ₁ ≈ 4·10⁵.
* **Analyticity radius [NE].** Float probes of the integrand on the torus show winding of the log arguments reaching π (a branch
  cut is crossed) at ρ ≈ 0.2, and `sup|G − G*| = 0.29 / 0.76 / 1.65 / 4.4` at ρ = 0.02 / 0.05 / 0.1 / 0.2. So ρ_an ≈ 0.15.
* **Size of the true remainder [NE].** The l9 DEG-8 records give the next term directly. At Δ = 100, x = 0.5: |S6/S0| = 2·10⁻² for
  the scalar and |M6/M0| = 5·10⁻² for C2 J = 10, i.e. K₃ coefficients ~2·10⁴ and ~2·10⁻⁵ relative at Δ = 1000. That is ≈ 20 c₂,
  far inside every F_J. So the K₂ truncation is accurate; the obstacle is only to *certify* it.
* **What a better lemma could reach (§0 item 1) [NE, rough estimate].**
  * With an exact holomorphic Morse chart there is no `e^{κ|u|²}`, and the bound is `rel ≈ 210 M_h (ε/ρ_y)⁶` with ρ_y ≈ ρ_an.
  * At ρ_y = 0.1: Δ₁ ≈ 2.8·10³ (J = 0) and 1.3·10⁴ (J = 7).
  * At Δ = 1000 the optimal-truncation floor `e^{−Δρ_an²/2}` is 7·10⁻³ (ρ_an = 0.1) to 1.3·10⁻⁵ (0.15), against rel*₇ = 4.2·10⁻⁴.
    Certifying at Δ_g = 1000 for J = 7 would need a certified Morse chart on a polydisc of radius ≳ 0.13, plus high-order
    certified Laplace coefficients.

## 6. Reproduce

```
cd ads4/code/gapc3 && sh run_all.sh        # ~1.5 CPU-h, <= 2 processes; outputs in cache/
```

| file | content | grade |
|---|---|---|
| `gc_guard.py` | write guard (only `gapc3/cache`) | — |
| `gc_ball.py` | Arb Taylor-model backend with the l1_tps interface; loads ball copies of l1_crystal, l5_spin, l6, l9, w0, p3 | [certified] (enclosures) |
| `gc_scalar.py` | log-factored scalar integrand (point and TM backends) | [DV] (validation) |
| `gc_node.py` | Morse matrix, exact crystal ball, project builders on the TM backend (`build`, `node`) | [certified] |
| `gc_cert.py` | per-node certified data: q, E, ∇Φ, B_k, rem_Φ, S₃, B_k(h), rem_h, M_h | [certified] |
| `gc_lemma.py` | Lemma L3 in Arb: `rel_bound`, `best` (over τ), `delta1` | [P] |
| `gc_run.py` | scalar batch → `cache/gc_cert_J*.json` | [certified] |
| `gc_table.py` | Δ₁, C₃ table → `cache/gc_table.txt/.json` | [certified] + [M] (targets) |
| `gc_spin.py` | EE / C3H amplitude sup norms via the project builders → `cache/gc_spin.json` | [certified] (frozen ⟨TTO⟩) |
| `gc_ideal.py` | sampled (ideal) sup norms → `cache/gc_ideal.json` | [NE] |
| `gc_report.py` | spinning and ideal Δ₁ tables | [certified] / [NE] |
