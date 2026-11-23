# P_cyc: the physical chain is D′ by definition; a direct EE test of the crystal cycle; an extended census

*23 Nov 2026, P_cyc agent. Unrefereed, not committed. No existing document or code file was edited. New code in `ads4/code/pcyc/`:
`pcyc_region.py`, `pcyc_dlib.py`, `pcyc_dquad.py`, `pcyc_crys.py`, `pcyc_cv.py`, `pcyc_compare.py`, `pcyc_exact.py`, `pcyc_farw.py`,
`pcyc_census.py`, helpers `dqsum.py`, `crsum.py`, `census_table.py`, `mk_table.py`, and the queue scripts `q_conv.sh`, `dimscan.sh`,
`dimscanE.sh`, `dq_queue.sh`, `dq_prod.sh`, `cr_queue.sh`, `cr_seed2.sh`, `d50_queue.sh`, `exact_more.sh`, `farw_q.sh`,
`census_queue.sh`, `census_queue2.sh`, `census_big.sh` (`q_conv.sh`, `dq_queue.sh`, `census_queue.sh` were superseded mid-run). Caches in `ads4/code/pcyc/cache/` (git-ignored by its own `.gitignore`). Read-only imports:
`rd/rd_core`, `rd/rd_saddle`, `rd/rd_census` (RD_CLOSE), `spinfunc/n4_thimble`, `n3_hybrid`, `n3_flow`, `n3_thimble`, `n1_pw`,
`l1_crystal`, `n4_exact`, `n3_exact` (and through them `exactskel/skelexact`), `repro_engine/integrand` (through `rd_core.Qratio`).
No project cache was written (the RD saddle cache was redirected to `pcyc/cache`). The CLS-D source (arXiv:2311.04271, `draft.tex`) was
read from the session scratchpad. Nothing in `A1PRIME.md` or `code/a1p/` was read or touched. Compute ≈ 5 CPU-h, at most 2
processes, each under `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]** numerical evidence;
**[M]** model; **[A]** assumption; **[I]** open. Notation as in RD_CLOSE (raw `z = (x̃₀, u_x, v_x, u_y, v_y, x⊥, y⊥)`, `a = 1/Δ`,
`s = u + v`, `q = uv − x⊥²`, `B_x = (1 − s_x + q_x)/(1 + s_x + q_x)`, `Q = Φ_EE/(aΦ_scal)`).

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **(A1) The physical integration chain** (§1) | CLS-D define the functional as an integral over real Lorentzian points in the causal order `4 > 5 > 3 > 0 > 1 > 2` (Lorentzian shadow representation of the block). In the gauge `x₃ = −x₁ = y`, `x₄ = −x₂ = x/(−x²)`, `x₅ = e` this region **is** `D′`, the real section cut at the light cones of `x₅` (= the zeros of `B`). 2·10⁵ random configurations: 0 mismatches. On `D′` all 25 bases of the integrand are positive, `A, B > 0`, `Q` real, the integral converges absolutely: physical integrand = (ν-dependent constant) × the project's integrand. | **[P]** region + **[D]** integrand, **[DV]** |
| **End to end, scalar, including low Δ** (§1.4) | `∫_{D′}` of the project's scalar integrand, **exact in Δ with the cut**, against the exact CHMRS Mellin skeleton, ratio to Δ = 100: `|ε| ≤ 2·10⁻⁸` (Δ = 20), `≤ 6·10⁻⁹` (Δ = 30) for ν ≤ 5, at the quadrature error up to ν = 10–12. At these Δ the part of `D′` where the `1/Δ` series diverges carries O(1) weight, so this checks the chain itself, cut included. Also Δ = 50: ≤ 1·10⁻⁹ to ν = 10. | **[DV]** |
| **`P_cyc` as needed by `H_RS`** | "the physical cycle is in the class of `D′`" is an identity: the physical cycle *is* `D′`. | **[D]** |
| **`H_RS`** (engine = physical integral, Q's binding state) | = RD_CLOSE Prop. D′ (engine = `∫_{D′}` + `O(Δ^N e^{−0.336Δ})`) + §1. Numerically: REPRO_ENGINE's direct `∫_{D′}` (all nodes inside `D′`) = engine `K = 5` to 3·10⁻⁷. **`P_cyc` (crystal) is no longer a premise.** | **[D]** (one [D-sketch] sub-step inherited from Prop. D′), [DV] |
| (A2) crystal cycle vs `D′`, Regge-scaled system (§2.1) | `D` = the principal product thimble; intersection numbers `δ_{k0}` with every sheet; rank 1 in the populated sector. | **[P]** |
| (A2) scalar → EE transfer (§2.2) | `Q` rational, poles only on the shared divisors: classes transfer with the same constant (not periods). | **[D]** |
| (A2) crystal cycle vs `D′`, finite Δ (§2.3) | gap: no Stokes jump of the crystal cycle along `a ∈ (0, 1/Δ]`, `n± = 1` through the fold. | **[I]** |
| **(B) direct EE test, above the fold** (§3, n3, Δ = 50, 100, ν = 10–16) | `T = ⟨Q⟩_crys/⟨Q⟩_{D′} − 1` (normalisation-free; `P_cyc` ⇒ `T = 0` exactly): consistent at every point (`|T|/σ ≤ 1.05`); σ = 3·10⁻³ – 4·10⁻² (seed scatter and `D′` grid error included); best `(|T| + 2σ)/thr` = 4.8 at (Δ, ν) = (100, 16). At Δ = 30 the `D′` EE integral is not resolved for ν ≥ 10. | **[NE, weak]** |
| **(B) direct EE test, below the fold** (n4, Δ = 30, 100, ν = 3–10) | consistent (`|T|/σ ≤ 1.7`; 2.5σ at (100, 4), a 4-replica run), but σ = 0.1–1.4, i.e. 16–2000× the threshold. | **[I]** (uninformative) |
| (B) ν < 3 | not run: the scalar thimble is not EE-admissible there (RD_CLOSE §3.4). | [I] |
| **Census** (§4) | 12 × 4000 starts (Δ = 30–1000, ν = 2, 4, 8): critical points **do** exist in the far region `0.3 ≤ r < 1.5` (0–4 at Δ ≥ 300, 6–276 at Δ ≤ 100), many with `|B| ≈ 1` and heights +30…+94 above the crystal. RD_CLOSE's pilot "none at raw radius 0.1–39" is withdrawn. Not saturated (> 10³ points per run). | **[NE]** |

**Bottom line.**
* **The physical integration cycle is not the crystal cycle; it is `D′`.** CLS-D define the functional as an integral over real
  Lorentzian points in the causal order `4 > 5 > 3 > 0 > 1 > 2`, and in the project's gauge that region is exactly the real section cut
  at the light cones of `x₅`. The saddle-point (crystal) cycle is CLS-D's steepest-descent *representation* of this chain, introduced
  only to evaluate it at large Δ, ν.
* So `P_cyc`, in the form `H_RS` needs it, is an identity. With RD_CLOSE's Prop. D′ (engine = `∫_{D′}` up to `O(Δ^N e^{−0.336Δ})`), **`H_RS`
  is [D]**, and REPRO_ENGINE's direct `∫_{D′}` at the binding state is its numerical check [DV].
* End to end for the scalar: `∫_{D′}` exact in Δ equals the exact Mellin skeleton to 10⁻⁹–10⁻⁸ at **Δ = 20, 30, 50**, where 3–80 % of
  the absolute weight sits beyond raw radius 0.3, outside the convergence region of the `1/Δ` expansion. This checks the chain, its cut
  and its phases non-perturbatively [DV].
* The crystal-cycle question survives only as a question about the saddle-point *numerics*. In the Regge-scaled system the crystal is
  `D` [P]; at finite Δ it is open (a Stokes-wall question). The direct EE test of those numerics is consistent everywhere but weak
  (best: within a factor ≈ 5 of the failure threshold above the fold; uninformative below it).
* The extended census withdraws RD_CLOSE's "no critical point in the far region": there are many, with `|B| ≈ 1` and heights e^{30–90}
  above the crystal. This does not touch `D′` (its scalar integral is the exact one) but shows that the Picard–Lefschetz route to the
  crystal-cycle identity at finite Δ is far from closable.

---

## 1. Route (A), first half: the physical integration chain *is* D′

### 1.1 The definition (CLS-D §3.1–3.2, read from the source)
* The functional is a spacetime integral (eq. `scalarfunctional_plus_final`):
  `Ψ⁺_{k,ν}[𝒢] = 4A_{k,ν} ∫_{4>3>0>1>2} d^dx₁…d^dx₄ d^dx₀ /vol(SO(d,2)) ⟨…⟩⁻¹ ⟨celestial⟩ δ(x²_{1⁺0})δ(x²_{30})θθ ⟨Ω|[φ₄,φ₃][φ₁⁺,φ₂⁺]|Ω⟩`,
  over **real Lorentzian points** with the causal order `4 > 3 > 0 > 1 > 2`.
* On a block it uses the **Lorentzian shadow representation** (eq. `block_shadowrep0`), which integrates `x₅` over the causal diamond
  `4 > 5 > 3` with absolute-valued 3-pt structures `|⟨φφO(x₅)⟩||⟨Õ(x₅)φφ⟩|`. This represents the block itself (not the partial wave).
* So `Ψ⁺[G_{Δ,J}]` is one absolutely convergent integral over the real region `4 > 5 > 3 > 0 > 1 > 2`. After the gauge fixing
  (eq. `gaugefixing`: `x₃ = −x₁ = y`, `x₄ = −x₂ = x/(−x²)`, `x₅ = e`, `x, y > 0`; FP determinant eq. `FP_Det`; `z₀` fixed; the two
  delta functions put `x₀` on a circle, `S^{d−3} = {n₀ = ±1}` in d = 3) the remaining integral is over `(x, y, x₀)` in the image of that
  region, with exactly the measure and factors of `l1_crystal.build` (`meas`, `(−x²)^{−6}`, FP, dual, celestial, 3-pt, subtraction).
* The crystal saddle enters CLS-D only afterwards (§3.3: "deform the integration contour into steepest-descent flows and see which
  saddle locus the contour goes through"), as an evaluation method for `Δ, ν → ∞`. **The physical cycle is the real region; the
  crystal cycle is a deformation of it.**

### 1.2 The region in the gauge is exactly D′ [P + DV]
* `y > 0` and `3 > 0 > 1`: `y` future timelike, `x₀` on the circle (built into the coordinates).
* `5 > 3`: `e − y` future timelike, i.e. `(1 − u_y)(1 − v_y) > y⊥²`, `u_y < 1`.
* `4 > 5`: with `x₄ = x/(−x²)`, `(x₄ − e)² = −(e − x)²/x²` and `x² < 0`, so `x₄ − e` is timelike iff `e − x` is; future orientation by
  continuity from `x = e/2` (`x₄ = 2e`). So `4 > 5` ⇔ `(1 − u_x)(1 − v_x) > x⊥²`, `u_x < 1`.
* `4 > 2`, `1 > 2`, `4 > 3` follow by transitivity.
* Hence the physical region is `{1 − s_x + q_x > 0, u_x < 1} × {1 − s_y + q_y > 0, u_y < 1} × circle` — **exactly RD_CLOSE's `D′`**
  (the real section cut at the light cones of `x₅ = e`, which are the zeros of `B`, `m₄₅ = 0` and `m₃₅ = 0`).
* Check (`pcyc_region.py` → `cache/pcyc_region.txt`): 2·10⁵ random real configurations (u, v log-uniform in [10⁻³, 30] and [10⁻³, 1]),
  the causal predicates evaluated by an independent Minkowski implementation: **116 540 in both, 0 in only one**. On the 83 460
  configurations outside, `4 > 5` fails 57 122 times, `5 > 3` 57 360 times, and `B ≤ 0` on 53 000 of them.

### 1.3 The integrand on D′ is the physical one up to one constant [D + DV]
* CLS-D's integrand is built from absolute-valued structures and complex powers of real invariants. The project's integrand is the
  analytic function with the same factors.
* On `D′` **every base `f_k` of a non-integer or ν-dependent power is strictly positive** (`Nn, c₁₃, m₂₄, m₄₀, m₂₀, c₁₀, c₃₀, m₃₅,
  m₄₅, m₃₄, uv, m_A, m_B`; also the integer-power ones `m₁₂, m₁₅, m₂₅, −x², FP₁₋₃, v′−u′, m₁₃, m₁₄, m₂₃`), `A > 0`, `B > 0`, and `Q`
  is real (all 116 540 points; `min A = 2.7·10⁻⁷`, `min B = 3.5·10⁻⁷`). RD_CLOSE §2.2 found the same with 2·10⁵ other samples.
* `D′` is connected (the two sheets `n₀ = ±1` join at `y₀ = 0`, where the integrand is smooth in the angle `ψ`).
* Therefore physical integrand = `c(ν, Δ)` × analytic integrand on all of `D′`, with `c` independent of the point. Its only
  Δ-dependence is through `B^{Δ/2}` (no phase: `B > 0`) and the constants `β_{Δ,J}`, `sin²`; its ν-dependence is a phase convention for
  `A^{iν/2}`-type factors. Scalar and EE share it up to a sign (their exponents differ by integers). It cancels in `ρ = I_st/I_ref`
  after the ν = 0.02 normalisation (HRS_CHECK §1.1).
* Absolute convergence on `D̄′` (RD_CLOSE §2.2 [D + DV]): the atom divisors have exponents > −1 for scalar and EE; the `m₁₂ = 0`
  corner is integrable for Δ ≥ 27; the `B` zeros carry exponent `(Δ − 7)/2`. So the physical integral is an ordinary Lebesgue integral
  over `D′`, and the gauge fixing (Fubini) is legitimate.

### 1.4 End-to-end check for the scalar, now at Δ = 20, 30, 50 [DV]
* Before: the engine (= `∫_{D′}` up to `e^{−cΔ}`) equals the exact CHMRS Mellin skeleton to ≤ 5·10⁻⁸ (Δ = 100), 6·10⁻⁹ (300),
  1·10⁻⁹ (1000) (RD_CLOSE §4.1; `spinfunc/cache/bg_validate.txt`). At those Δ the part of `D′` beyond raw radius 0.3 weighs ≤ e^{−15}
  (Δ = 100) of `∫|Φ|`, so that check sees `D′` only through its Regge region.
* **New, Δ = 20, 30, 50** (`pcyc_exact.py`, `exact_more.sh` → `cache/pcyc_exact_D{20,30,50}_D100.txt`). The ratio
  `R_{D′}(ν) = I_S(Δ, ν)/I_S(100, ν)` of the **exact-in-Δ** integrals over `D′` (the product rule of §3.2, cut at the light cones of `x₅`
  included) against the exact skeleton `X(ν; Δ)/X(ν; 100)` (CHMRS Mellin; `n4_exact.record` / `n3_exact.X`, computed here, seconds per
  state), `ε = log(R_{D′}/R_ex) − const` (const fitted on ν ≤ 2: the dropped Δ-only factors). `|ε|` with the fine grids, and below it
  the coarse-minus-fine change of `ε` (a conservative estimate of the fine grid's quadrature error):

  | Δ | ν = 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 14 | 16 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | 20 | 7.2e-9 | 7.2e-9 | 7.1e-9 | 2.0e-8 | 6.5e-8 | 1.6e-6 | 2.0e-6 | 2.4e-5 | 1.0e-4 | 3e-2 | 0.6 |
  | (quad.) | 9e-9 | 9e-9 | 7e-8 | 5e-7 | 4e-6 | 4e-5 | 9e-5 | 3e-3 | 7e-2 | 0.5 | 1 |
  | 30 | 4.6e-9 | 4.6e-9 | 4.9e-9 | 5.2e-9 | 5.7e-9 | 7.9e-8 | 2.1e-8 | 7.2e-8 | 6.9e-7 | 1.8e-5 | 6.2e-4 |
  | (quad.) | 5e-10 | 5e-10 | 5e-9 | 3e-8 | 2e-7 | 5e-6 | 5e-7 | 1e-5 | 3e-4 | 3e-3 | 2e-2 |
  | 50 | 6.8e-10 | 6.8e-10 | 9.8e-10 | 9.9e-10 | 9.9e-10 | 5.3e-9 | 1.1e-9 | 1.1e-9 | 3.8e-8 | 1.8e-6 | 1.1e-4 |
  | (quad.) | 3e-11 | 3e-11 | 3e-10 | 2e-9 | 9e-9 | 9e-7 | 1e-8 | 1e-7 | 4e-7 | 3e-6 | 1e-4 |

  Wherever the quadrature resolves the integral, the residual is at the 10⁻⁹–10⁻⁸ level: **the exact functional equals `∫_{D′}` of the
  project's integrand at Δ = 20, 30, 50, cut included.** (ν = 6 sits next to the zero of
  `g + g̃` at 5.98, where every comparison is ill-conditioned.)
* Why this is a non-perturbative check (`pcyc_farw.py` → `cache/pcyc_farw.txt`): the share of `∫_{D′}|Φ_scal|` (ν-independent on
  `D′`) in the region `max(s_x, s_y) ≥ r` (raw light-cone sums; the `a`-series of `B` converges only for `r < 1/3`, RD_CLOSE §2.3) is

  | Δ | r ≥ 0.2 | r ≥ 0.3 | r ≥ 0.5 |
  |---|---|---|---|
  | 20 | 0.98 | 0.82 | 0.17 |
  | 30 | 0.81 | 0.47 | 0.027 |
  | 50 | 0.16 | 0.033 | 1.3·10⁻⁵ |
  | 100 | 3·10⁻⁴ | 3·10⁻⁷ | 10⁻¹⁴ |

  (node mask on the product rule at ν = 0; a rough estimate). At Δ = 20 and 30 most of the absolute weight lies where the `1/Δ`
  expansion does not converge and up to the light cones of `x₅`. A wrong chain (the uncut real section, a cut elsewhere) or a wrong
  phase on part of `D′` would show up at O(10⁻²–1) there; the agreement is 10⁻⁹–10⁻⁸. **[DV]**
* This tests the identification "physical = `∫_{D′}`" (region, integrand, phases, normalisation up to Δ-only constants) for the scalar.
  EE uses the same region and the same factors times the real rational `Q` (§1.3).

### 1.5 Consequence for H_RS
* `H_RS` (engine = physical integral) needs two statements:
  1. physical integral = `c · ∫_{D′}Φ` — §1.2–1.4, **[D]**, end-to-end **[DV]** for the scalar;
  2. `∫_{D′}Φ_EE` = engine + `O(a^{K+1})` + `O(Δ^N e^{−cΔ})`, `c = 2 artanh(r/2) ≥ 0.336` — RD_CLOSE Prop. D′, **[D]** (its sub-step "the
     `a`-series of the `x₀`-dependent slow factors converges on the polydisc of raw radius 1/3" is [D-sketch], supported by their
     sign-definiteness on `D′` [DV] and HRS_CLOSE B6), with HRS_CLOSE Prop. RS [D] for the order-by-order content.
* Numerically: REPRO_ENGINE's direct finite-Δ product rule (all nodes inside `D′` at Δ_g = 300 and 1000: its largest radial nodes,
  `X ≈ 34` (10-point) and `Y ≈ 52` (16-point), give `s_x ≤ 0.23`, `s_y ≤ 0.35` at 300) equals the engine's `K = 5` ratio to 3.4·10⁻⁷ (300) / 3·10⁻⁸ (1000) on ν < 4. That is a direct
  `∫_{D′}` evaluation of the EE integrand at the binding state, **[DV]**.
* **`P_cyc` is therefore not a premise of `H_RS`.** The crystal (saddle-point) cycle is CLS-D's steepest-descent *representation* of
  `D′`; in this project it enters only the crystal-cycle numerics (n1/n3/n4, the `K₂` Laplace rows, HRS_CHECK). Whether those numerics
  compute the same object for EE is the question of §2–3.
* What the definitional route does *not* cover: the identification of the project's EE integrand with CLS-D's spinning functional
  (shared definition, cross-checked point-wise by REPRO_ENGINE's independent implementation), and everything outside H_RS.

## 2. Route (A), second half: the crystal cycle versus D′ (Picard–Lefschetz)

### 2.1 The Regge-scaled local system: D is the principal product thimble [P]
* At `a → 0` (HRS_CLOSE §2.1) the integrand in the D-coordinates is `Π_j w_j(x_j)` times a Laurent polynomial, with one-variable
  twisted factors: four radial `X^{c_j + iνe_j} e^{−X}` on `(0, ∞)`; `t = s²`: `t^{−1/2}(1 − t)^{g + iνe}`; `cos^{h+iν}φ`; `sin^{k−iν}α`
  (exponents measured in `pcyc_dlib`, e.g. scalar `X₁^{2.25}`, `w^{1.75}`, `cos^{4.5}φ`, `sin^{−1/2}α`; EE differs by integers).
* In log coordinates (`X = e^ξ`, logistic for `t`, …) each factor is `e^{f_j(ξ)}` with `Re f_j → −∞` at both ends of the real line
  (every real exponent is > −1). Its critical points are `ξ_{j,k} = ξ_j* + 2πik`: one per sheet of the logarithm, i.e. exactly the
  "windings of `arg A`" that label the sheets in RD_CLOSE's census.
* Intersection numbers with the steepest-ascent cycles: `⟨ℝ, K_{j,k}⟩ = δ_{k0}` whenever the exponent has positive real part (the
  `Γ`/Beta-function statement `∫_ℝ e^{λξ − e^ξ}dξ = ∫_{J_0} = Γ(λ)`, `Re λ > 0`). The Morse function is a sum, so the decomposition of
  the product cycle is the product: **`D = Π_j J_{j,0}` = the thimble of the principal crystal of the Regge-scaled system, with
  intersection number 1, and every other critical point (every other sheet) has coefficient 0** [P].
* Each one-variable group is of rank 1 in the parity sector the integrand populates (radial: `ℂ*` with an exponential end; `t`:
  `ℂ \ {0, 1}`; the `φ`- and `α`-odd classes are exact or unpopulated, HRS_CLOSE §2.3, §3.3). So every closed cycle of the Regge-scaled
  system is `c_Γ(ν)·D`. This is HRS_CLOSE's Prop. RS in Picard–Lefschetz form.
* Boundary pieces: in the Regge-scaled system the `B` zero has moved to `X = ∞` and is the exponential end `e^{−X}` of both `D` and the
  thimble. At finite Δ it is the divisor `m₄₅ = 0` / `m₃₅ = 0`: the boundary of `D′` by definition, and the measured terminal divisor
  of 35–60 % of the crystal flow lines (RD_CLOSE §3.3 [NE]). The other terminal divisors of the crystal lines (`s → ±1`, `X_i → 0`, the
  angular corner) are atom divisors, where `D′` also ends.

### 2.2 Transfer from the scalar to EE [D]
* `Φ_EE = aQΦ_scal` with `Q` single valued and rational, all its poles on divisors of the scalar local system (RD_CLOSE §1.2 [P + DV]).
  So `ω_EE` and `ω_scal` are sections of the same twisted de Rham complex on the same complement `U`.
* Hence any homology relation `[Γ₁] = c[Γ₂]` in `H₇(U, L)` gives `∫_{Γ₁}ω_EE = c∫_{Γ₂}ω_EE` with the **same** `c`, provided both chains
  are admissible for `ω_EE` (integrable at their boundary divisors, so that the regularised and unregularised integrals agree) [D].
  * `D′` is admissible for EE (RD_CLOSE §2.2).
  * The crystal thimbles are admissible for EE except at the angular corner `y² = 0, α → ±i∞` for ν ≲ 2 (RD_CLOSE §3.4). There the
    class transfers but the untruncated scalar thimble must be replaced by its Hankel regularisation for EE.
* The transfer works for **classes**, not for **periods**: a scalar period check alone cannot exclude a class `F` with
  `∫_F ω_scal = 0 ≠ ∫_F ω_EE` (RD_CLOSE §4.2). That is why the scalar's [DV] does not by itself make `P_cyc` [D] for EE.

### 2.3 Finite Δ: what a proof needs, and the gap [D-sketch / I]
* At finite Δ the local system has further divisors at raw radius ~ 1 (`FP₁ = 0`, `FP₂ = 0`, the poles of `B`, the `x₀`-divisors
  far out), and the rank of the twisted homology grows. `D′` and the crystal cycle are both Regge-region chains with radial ends on the
  `B` zeros. Their difference `F` has all periods vanishing to every order in `1/Δ` (Prop. RS). A proof that `F = 0` needs one of:
  1. **continuity in `a`**: the class of the crystal thimble is locally constant in `a ∈ (0, 1/Δ]` away from Stokes walls; at `a → 0`
     it is `D` (§2.1). Missing: the absence of Stokes jumps of the crystal cycle along `a ∈ (0, 1/Δ_g]` at fixed ν, and across the fold
     ν* ≈ 10–11 where `S` splits into `A±` (`n± = 1` is calibrated against n3 at ν_c = 10, not derived);
  2. **a complete census** of the critical points of the flow function at the given (Δ, ν) with their intersection numbers
     `⟨D′, K_σ⟩`. A complete census of this seven-dimensional complex system is out of reach here; §4 extends the pilot.
* **Precise gap** [I]: *no Stokes jump of the crystal cycle along `a ∈ (0, 1/Δ]` at fixed ν < 10, and the `A±` continuation from ν_c with
  `n± = 1` stays in the class of the principal product thimble.* It is supported by the scalar periods (engine = Mellin skeleton to
  ≤ 10⁻⁹ at Δ = 300/1000; §1.4 at Δ = 20–50), the n3/n4 scalar validation (HRS_CHECK §2.2), §3 and §4.
* **(A) verdict.** The question "is the crystal cycle a divisor-free deformation of the physical cycle" is answered for the
  Regge-scaled system [P], and transferred to EE at the level of classes [D]; at finite Δ it stays open with the gap above. It is no
  longer load-bearing, because the physical cycle is `D′` itself (§1).

## 3. Route (B): direct test of the crystal cycle for the EE form at Δ = 30, 50 and 100

### 3.1 Test statistic [D]
* Because `Φ_EE = aQΦ_scal` point-wise, the EE crystal integral is `a∫_crys QΦ_scal`. Define
  ```
  ⟨Q⟩_Γ = ∫_Γ Q Φ_scal / ∫_Γ Φ_scal,        T(ν, Δ) = ⟨Q⟩_crys / ⟨Q⟩_{D′} − 1.
  ```
  `P_cyc` (crystal = `c(ν)·D′` as twisted cycles, §2.2) implies **`T = 0` exactly**, at every Δ and ν, all orders and non-perturbatively.
  No normalisation, branch constant, `a`-power, `β_{Δ,J}` or `c(ν)` enters: each cancels inside each ratio. No EE branch bookkeeping
  enters either: `Q` is single valued and is evaluated at the contour point with the `y₀` sheet replayed from the scalar integrand.
* `T` measures `δ_F(Δ, ν) = ∫_F ω_EE / ∫_{D′} ω_EE` of a foreign class `F` (with `∫_F ω_scal ≈ 0`).
* **Failure threshold, scaled to Δ = 30 / 50 / 100.** A worst-sign relative error of `ρ = I_st/I_ref` of 2.1·10⁻⁴ (Δ_g = 300) / 1.9·10⁻⁵
  (1000) on ν < 4 moves `s_d` to `s_fail` (REPRO_ENGINE §5). A Δ-independent `δ_F` cancels in `ρ` (state and reference at the same ν).
  The RD_CLOSE §4.2 scenario, the same in `c₂` at both Δ_g, is `δ_F(Δ) − δ_F(∞) = κ/Δ²` with `κ ≥ 2.1·10⁻⁴·300² = 19` (and 1.9·10⁻⁵·1000²
  = 19 from 1000: the same). At the test points this means
  * `|δ_F(30) − δ_F(∞)| ≥ 2.1·10⁻²`, `|δ_F(50) − δ_F(∞)| ≥ 7.6·10⁻³`, `|δ_F(100) − δ_F(∞)| ≥ 1.9·10⁻³`;
  * model-independently in the Δ-dependence of `κ/Δ²`: `|δ_F(30) − δ_F(100)| ≥ κ(1/30² − 1/100²) = 1.9·10⁻²`.
  An exponentially decaying `δ_F ∝ e^{−cΔ}` is larger at low Δ, so a low-Δ bound is conservative for it.

### 3.2 The D′ side: exact-in-Δ product rule [DV]
* **Coordinates** (`pcyc_dlib.py`). `D′` becomes a product domain: per light-cone pair `(u, v) = σ(r, 1 − r)` with
  `B_x = exp(−2t/Δ)` solved exactly for `σ(t, k)` (`k = r(1−r)(1−s²)`), so that `t ∈ (0, ∞)` ends exactly on the light cone of `x₅`
  (the `B` zero) and `B^{Δ/2} = e^{−t}`; `s, φ, α` as in BELOW_FOLD (`ψ = α − φ`, both sheets).
* **Rule** (`pcyc_dquad.py`, own code): 7-dimensional tensor-product product integration. Generalised Gauss–Laguerre in `t_x, t_y` (weight
  `t^p e^{−λt}`, `p = 6/8`, `λ = 1 − 1/Δ` scalar, `1 − 5/Δ` EE: the measured decay rates toward the `B` zero), Gauss–Jacobi in `r_x, r_y`,
  `t = s²`, `φ`, `α` with the measured end exponents (scalar `r^{2.25}(1−r)^{1.75}`, `w^{1.75}`, `cos^{4.5}φ`, `sin^{−1/2}α`; EE the same
  except `r_y^{0.25}(1−r_y)^{−0.25}`, `cos^{0.5}φ`), the ν-atoms `r^{∓iν/2}(1−r)^{±iν/2}`, `(1−t)^{−iν/2}`, `cos^{iν}φ`, `sin^{−iν}α`
  integrated exactly in the 1-D weights (fine tanh–sinh; rule self-test ≤ 3·10⁻¹¹ against Beta functions), and the remainder
  `exp(iν ε/2)` of `A^{iν/2}` evaluated exactly at the nodes (`|ε| ≤ 2.0` at Δ = 30, 0.07 at Δ = 100; no expansion). Both signs of `s`
  are evaluated, which makes the remainder analytic in all seven variables (transverse parity).
* **Convergence** (`cache/dimscan_30_S.txt`, `cache/log_dq_*`): per-dimension scan for the scalar at Δ = 30 (each dimension
  10 → 16 nodes changes `I_S` by ≤ 2·10⁻⁷ at ν ≤ 5, ≤ 2·10⁻⁴ at ν ≤ 10). Production pairs (coarse → fine):

  | Δ | kind | grid (fine) | rel. change coarse → fine at ν = 1 / 2 / 3 / 4 / 5 / 8 / 10 / 12 / 14 / 16 |
  |---|---|---|---|
  | 30 | S | 14,16,8,14,16,16,16 | 2e-10 / 9e-10 / 5e-9 / 3e-8 / 2e-7 / 5e-7 / 1e-5 / 3e-4 / 3e-3 / 2e-2 |
  | 30 | E | 12,14,8,12,14,24,16 | 3e-4 / 3e-5 / 2e-5 / 5e-6 / 8e-5 / 3e-3 / 4e-2 / 0.4 / 2 / 1 |
  | 100 | S | 14,16,8,14,16,16,16 | 8e-13 / 4e-12 / 1e-11 / 7e-11 / 3e-10 / 8e-10 / 4e-9 / 7e-8 / 2e-6 / 2e-4 |
  | 100 | E | 12,14,8,12,14,24,16 | 4e-4 / 4e-5 / 3e-5 / 1e-5 / 5e-5 / 3e-5 / 3e-5 / 1e-4 / 4e-4 / 4e-4 |

  The EE integrand at Δ = 30 and ν ≥ 10 is not resolved (cancellation `∫|Φ_EE|/|∫Φ_EE|` grows like `e^{1.6ν}` and the EE weight is broad
  at low Δ); those points are not used. ν = 1 lies next to the zero of the EE image (ν ≈ 1.07, HRS_CHECK), ν = 6 next to that of
  the scalar (5.98): relative errors are larger there.
* **`⟨Q⟩_{D′}`** (fine grids): Δ = 100: ν = 4: 2.316628·10⁻⁶; 8: 9.88642·10⁻⁷; 10: −6.67481·10⁻⁶; 12: −5.75613·10⁻⁶; 14: 8.2582·10⁻⁶;
  16: 4.1490·10⁻⁵. Δ = 50: ν = 12: −9.633·10⁻⁵; 14: 1.2472·10⁻⁴; 16: 6.3330·10⁻⁴. Δ = 30: ν = 3: 4.67187·10⁻⁴; 4: 3.26982·10⁻⁴;
  8: 1.961·10⁻⁴.

### 3.3 The crystal side: the project's n4 / n3 contours with Q evaluated point-wise
* **Contours** (`pcyc_crys.py`): the project's crystal cycle, unchanged (imported read-only, no cache written):
  * below the fold, the n4 Lefschetz thimbles `A₊ + A₋` (`n± = 1`), saddles tracked from ν_c = 10, contour `n3_hybrid C(ρ₁ = 0.8, ν_f = ν)`;
  * above it, the n3 hybrid contour `C(3.5, ν_f)` of the symmetric saddle `S` (`ν_f = 13` for ν < 13, `max(ν, 15)` above);
  * as an internal check, ν = 10 is run with both (n4 `A₊ + A₋` and n3 with `ν_f = 13`).
* **Points and sheet.** At each contour point the raw `z` and the sheet of `y₀ = ±√y₀²` are read off by replaying the scalar integrand's own
  branch references (the same log sequence as `pint`). Check: `[Re P_int − Re log J_w − Σ Re(e_k) log|f_k|]/(−ν/2) − Arg A ∈ 2πℤ` to
  ≤ 3·10⁻¹² on every weighted point of every run (a wrong `y₀` sheet changes `|f_k|` and `A`). `Q` is evaluated there with `rd_core.Qratio`;
  its double-precision condition number (`Q_condition`) is ≤ 1·10⁻¹⁵ on the weighted points.
* **Estimator.** Scrambled Sobol in the contour parameter, antithetic, Gaussian proposal moment-matched to `|F|(1 + |Q|/⟨|Q|⟩)`; R = 8
  replicas, jackknife errors; `A₊`, `A₋` combined with their `e^{logN}`.
* **Control variates** (`pcyc_cv.py`). On any closed twisted cycle `∫ e^Φ CV = 0` for `CV = ∂_{w_j}ψ + ψ ∂_{w_j}Φ` (Stokes; nothing about
  `P_cyc` is assumed). With ψ the monomials of degree ≤ 2 in (`z_j/c_j − 1`, `x̃₀`, `x⊥/a`, `y⊥/a`) and the exact gradient of `rd_core`,
  `Q − q` and `1` are regressed on the 252 CVs (weighted complex least squares on the *other* replicas, so no fit bias). The CV sums
  are consistent with zero on every run (median |Σ|/√Σ² ≈ 0.4–1, the noise level), which also checks that the contours are closed.
  * On the n3 contour (Δ = 100, ν = 16, 2¹¹·2·8 points) the CVs cut the error of `⟨Q⟩` from 7.3·10⁻³ to 8.3·10⁻⁴ (degree 2; degree 3:
    6.5·10⁻⁴).
  * On the n4 thimbles (ν ≤ 10) they do not help (Δ = 100, ν = 4: 0.22 → 0.40): the variance there comes from the flow Jacobians
    ("fingers") and from 1–9 % of cut flow lines, not from the smooth variation of `Q`.
* **Cost.** n3: ≈ 0.4 ms per point. n4: ≈ 20 ms per point (the finite-difference flows), with a relative error of `⟨Q⟩` of 0.1–1 at
  2 × 1.6·10⁴ points. Reaching 10⁻³ below the fold would need ≈ 10⁸ points per (Δ, ν), i.e. ≈ 500 CPU-h: out of reach with this machinery.

### 3.4 Results against the failure threshold
`cache/pcyc_compare.txt` (per run and combined over independent seeds). `σ` = max(jackknife, seed scatter)/√n_runs ⊕ the `D′`
quadrature error. The threshold column is the `c₂`-scaled failure threshold of §3.1 (`2.1·10⁻⁴(300/Δ)²`; it is the ν < 4 threshold,
conservative for ν ≥ 10, whose own tolerance is 3–10× larger, HRS_CLOSE §5).

| Δ | ν | contour | runs | `⟨Q⟩` rel. error per run | `D′` error | Re T | Im T | σ | `|T|/σ` | thr. | `(|T| + 2σ)/thr` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 3 | n4 | 1 | 6e-01 | 2e-05 | +5.8e-01 | +3.8e-01 | 6.5e-01 | 1.06 | 2.1e-02 | 95 |
| 30 | 4 | n4 | 1 | 2e-01 | 5e-06 | -1.4e-01 | +3.3e-02 | 2.0e-01 | 0.73 | 2.1e-02 | 26 |
| 30 | 8 | n4 | 1 | 7e-01 | 3e-03 | +1.6e-01 | +1.2e+00 | 7.5e-01 | 1.67 | 2.1e-02 | 1.3e+02 |
| 30 | 10 | n3 | 1 | 9e-02 | 4e-02 | -9.6e-02 | +2.7e-02 | 9.6e-02 | 1.05 | 2.1e-02 | 14 |
| 30 | 10 | n4 | 1 | 1e-01 | 4e-02 | -2.4e-02 | -1.3e-01 | 1.1e-01 | 1.29 | 2.1e-02 | 16 |
| 30 | 12 | n3 | 1 | 6e-03 | 4e-01 | +4.9e-02 | +1.5e-02 | 4.1e-01 | 0.13 | 2.1e-02 | 41 |
| 30 | 14 | n3 | 1 | 1e-02 | 2e+00 | -4.1e-01 | -7.9e-02 | 1.5e+00 | 0.27 | 2.1e-02 | 1.7e+02 |
| 30 | 16 | n3 | 1 | 1e-03 | 1e+00 | -5.0e-01 | -8.2e-02 | 9.8e-01 | 0.52 | 2.1e-02 | 1.2e+02 |
| 50 | 12 | n3 | 2 | 3e-02/2e-02 | 9e-03 | -1.2e-02 | -1.3e-03 | 1.9e-02 | 0.65 | 7.6e-03 | 6.8 |
| 50 | 14 | n3 | 2 | 6e-03/2e-02 | 4e-02 | +7.0e-03 | +5.2e-04 | 4.1e-02 | 0.17 | 7.6e-03 | 12 |
| 50 | 16 | n3 | 2 | 8e-04/5e-03 | 3e-02 | -1.1e-03 | -1.3e-03 | 3.4e-02 | 0.05 | 7.6e-03 | 9.3 |
| 100 | 4 | n4 | 1 | 1e-01 | 1e-05 | -3.5e-01 | -4.5e-02 | 1.4e-01 | 2.46 | 1.9e-03 | 3.4e+02 |
| 100 | 8 | n4 | 1 | 1e+00 | 3e-05 | +6.4e-01 | +1.2e+00 | 1.4e+00 | 0.98 | 1.9e-03 | 2.1e+03 |
| 100 | 10 | n3 | 2 | 4e-02/1e-02 | 3e-05 | +2.2e-03 | +5.1e-03 | 1.9e-02 | 0.30 | 1.9e-03 | 23 |
| 100 | 12 | n3 | 2 | 3e-02/1e-02 | 1e-04 | -8.6e-03 | +1.0e-02 | 1.4e-02 | 0.95 | 1.9e-03 | 22 |
| 100 | 14 | n3 | 2 | 7e-03/2e-02 | 4e-04 | +1.0e-02 | -4.5e-03 | 1.4e-02 | 0.79 | 1.9e-03 | 21 |
| 100 | 16 | n3 | 2 | 8e-04/4e-03 | 4e-04 | +2.0e-03 | -1.1e-03 | 3.4e-03 | 0.67 | 1.9e-03 | 4.8 |

* **Above the fold (n3).** Every point is consistent with `T = 0` (`|T|/σ ≤ 1.05`). The jackknife over 8 replicas underestimates the
  error, as HRS_CHECK found for `k = 4`: at (100, 16) the two seeds give `T = −1.4·10⁻³ ± 0.9·10⁻³` and `+5.3·10⁻³ ± 3.7·10⁻³`, i.e. a
  scatter ≈ 4× the smaller jackknife error; the combined σ uses the sample scatter. The best sensitivity, `|T| + 2σ`, is 4.8× the
  threshold, at (Δ, ν) = (100, 16) (σ = 3.4·10⁻³). At Δ = 50, ν = 16 the two seeds agree with `D′` to −4·10⁻³ and +2·10⁻³, which on
  the crystal side alone would give ≈ 1× the threshold; but the conservative `D′` grid error there (coarse − fine, 3·10⁻²) dominates the
  quoted σ, and a finer `D′` grid was not run.
* **Below the fold (n4).** Consistent, but the `A₊ + A₋` thimble estimates of `⟨Q⟩` carry 10–100 % errors (flow-Jacobian
  "fingers", 1–9 % of cut flow lines, partial cancellation between `A₊` and `A₋`); the CVs do not reduce them, so the plain estimator
  is used. The 2.5σ at (100, 4) comes from a 4-replica run; with HRS_CHECK's ×3 rule for the `k = 4` jackknife it is 0.8σ. These
  points are 16–2000 times above the threshold and carry no weight for `P_cyc`.
* **Internal check at ν = 10, Δ = 30:** n3 (`ν_f = 13`) and n4 (`A₊ + A₋`) give `⟨Q⟩ = −7.8·10⁻⁴ (±10 %)` and `−8.4·10⁻⁴ (±10 %)`;
  the `D′` value `−8.6·10⁻⁴` is itself only good to 4 % there.
* **ν < 3** was not run: the untruncated scalar thimble is not EE-admissible at ν ≲ 2 (RD_CLOSE §3.4), so a direct EE thimble number
  there would need a new contour (radial thimble × real angles).
* **(B) verdict.** No discrepancy between the project's crystal-cycle numerics and `D′` for the EE form at any tested point
  (Δ = 30, 50, 100; ν = 3–16). The test is weak: it reaches the failure threshold within a factor ≈ 5 above the fold, and not
  at all below it, where the decisive weight of `s_d` lies. **[NE, weak]** above the fold, **[I]** (uninformative) below. Because of
  §1, this no longer affects `H_RS`.

## 4. Extended critical-point census [NE]
* `pcyc_census.py` runs RD_CLOSE's `rd_census.run` unchanged (damped Newton on the exact gradient of the w-chart flow function `Φ_w`,
  random far starts, de-duplication in raw coordinates); only the seed, the number of starts and the saddle-cache location differ.
  One run of 4000 starts takes ≈ 75 s.
* Counts per raw-radius class (raw radius `r = max |u, v|`), 4000 starts, seed 23; in brackets the largest `dH` (the sheet-independent
  part of `Re Φ_w` above the crystal; each further sheet of `A^{iν/2}` lowers it by `πν`) (`cache/census_table.txt`,
  `cache/census_D*_s23.txt`):

| Δ | ν | distinct | Regge `r < 0.1` | mid `0.1–0.3` | **FAR `0.3–1.5`** | beyond `1.5–10` | raw ∞ `≥ 10` |
|---|---|---|---|---|---|---|---|
| 30 | 2 | 1462 | 39 (+23.06) | 151 (+42.55) | 76 (+48.29) | 300 (+41.92) | 896 (+29.11) |
| 30 | 4 | 1513 | 16 (+19.82) | 201 (+41.69) | 119 (+49.59) | 344 (+40.73) | 833 (+35.31) |
| 30 | 8 | 1556 | 6 (+18.05) | 126 (+45.48) | 276 (+45.98) | 418 (+49.64) | 730 (+35.13) |
| 100 | 2 | 1129 | 180 (+64.77) | 32 (+68.24) | 6 (+29.75) | 22 (+53.84) | 889 (+44.08) |
| 100 | 4 | 1201 | 227 (+63.89) | 43 (+39.12) | 6 (+35.78) | 30 (+36.18) | 895 (+43.15) |
| 100 | 8 | 1199 | 154 (+61.38) | 136 (+65.46) | 15 (+94.07) | 45 (+46.90) | 849 (+37.26) |
| 300 | 2 | 1017 | 204 (+84.24) | 3 (+37.50) | 1 (+34.10) | 5 (+46.41) | 804 (+64.38) |
| 300 | 4 | 1073 | 263 (+83.34) | 1 (+50.76) | 1 (+33.38) | 6 (+46.07) | 802 (+53.73) |
| 300 | 8 | 1067 | 286 (+81.83) | 1 (+49.45) | 4 (+43.60) | 14 (+48.48) | 762 (+60.88) |
| 1000 | 2 | 980 | 191 (+58.77) | 0  | 0  | 0  | 789 (+75.85) |
| 1000 | 4 | 1007 | 227 (+105.07) | 0  | 2 (+61.11) | 2 (+31.76) | 776 (+78.51) |
| 1000 | 8 | 991 | 269 (+58.80) | 0  | 2 (+58.35) | 4 (+44.71) | 716 (+63.47) |

* **RD_CLOSE's pilot statement "no critical point at raw radius 0.1–39" (Δ = 300, ν = 4, 400 starts) is withdrawn.** With more starts the
  far region is populated at every (Δ, ν): in `0.3 ≤ r < 1.5` (where the `B` zeros and `FP₁ = 0` live) 1–4 points at Δ = 300, 0–2 at
  1000, 6–15 at 100, 76–276 at 30; and more beyond (up to `r ≈ 10`). Already the 600-start seed-17 run at (300, 4) finds two at `r = 6.49`.
* Many of them have `(Δ/2) log|B|` between −4 and −1, i.e. `|B(σ)| ≈ 1`: exactly the Δ_g-independent candidates of RD_CLOSE §4.2. Their
  sheet-independent heights are +30 to +94 above the crystal; some lie within 10⁻²–10⁻³ of `m₄₅ = 0` with `log(|Q|/Q_eff)` up to +71.
* The number of distinct points is not saturated (≈ 250–330 with 600 starts, ≈ 1000–1550 with 4000). A complete census with
  intersection numbers (§2.3, route 2) is out of reach.
* The Regge-class points (`r < 0.1`, heights up to +105) are critical points of the w-chart flow function, which is not a sum of
  one-variable functions; §2.1 uses the D-chart Morse function, for which the Regge-scaled decomposition is clean. Picard–Lefschetz
  decompositions depend on the Morse function; the homology classes do not.
* **What this means.** On the principal sheets these points sit e^{30}–e^{90} above the crystal. A non-zero Picard–Lefschetz coefficient
  of any of them in a cycle would change that cycle's scalar integral by such factors, unless pushed down by `k ≳ dH/(πν)` windings. For
  `D′` this is moot: its scalar integral *is* the exact functional (§1.4). For the crystal cycle such coefficients are excluded only by
  the scalar validation of the crystal numerics (HRS_CHECK §2.2), not by geometry. **[NE]**

## 5. Grades
| item | before | now | evidence |
|---|---|---|---|
| physical integration domain of CLS-D, in the gauge, = `D′` | implicit ("physical saddle-point cycle") | **[P]** (causal-relation algebra) + **[DV]** (2·10⁵ configurations, 0 mismatches) | §1.2 |
| physical integrand = `c(ν)` × the project's analytic integrand on `D′` | — | **[D] + [DV]** (all 25 bases positive, `A, B > 0`, `Q` real) | §1.3 |
| `∫_{D′}Φ_scal` = exact Mellin skeleton (end to end) | [DV] at Δ ≥ 100 via the engine | **[DV]** at Δ = 20, 30, 50 (≤ 2·10⁻⁸ for ν ≤ 5; at the quadrature error up to ν = 10–12), where 3–80 % of `∫|Φ|` lies beyond raw radius 0.3 | §1.4 |
| **`P_cyc` as needed by `H_RS`** (physical cycle in the class of `D′`) | [DV] scalar periods; not proved | **[D]: an identity** (the physical chain *is* `D′`) | §1 |
| `P_cyc` as a statement about the crystal (saddle-point) cycle, Regge-scaled system | (Prop. RS, order by order) | **[P]** (product thimble, intersection numbers δ_{k0}) | §2.1 |
| same, scalar → EE transfer of classes | (Lemma F, same local system) | **[D]** (classes; not periods) | §2.2 |
| same, at finite Δ | [I] | **[I]**, precise gap: no Stokes jump along `a ∈ (0, 1/Δ]`, `n± = 1` through the fold | §2.3 |
| direct EE test of the crystal cycle, ν ≥ 10 (n3) | none | **[NE, weak]** | §3.4 |
| direct EE test of the crystal cycle, 3 ≤ ν < 10 (n4) | none | **[I]** (uninformative) | §3.4 |
| ν < 3 | none | **[I]**, not run (angular corner, RD_CLOSE §3.4) | §3.4 |
| critical-point census | pilot (1 point, 400 starts) | **[NE]**: far region populated (0–276 points in `0.3 ≤ r < 1.5` per run), `|B| ≈ 1` candidates exist; pilot's "none" withdrawn; not saturated | §4 |
| **`H_RS`** (engine = physical integral at Q's binding state) | [D given P_cyc] | **[D]** (Prop. D′ [D] with one [D-sketch] sub-step + §1), **[DV]** numerically (REPRO_ENGINE's direct `∫_{D′}` = engine `K = 5`, 3·10⁻⁷) | §1.5 |

**Not claimed:**
* that the crystal-cycle numerics (n3/n4/`K₂` rows, HRS_CHECK) compute `∫_{D′}` for EE below ν ≈ 10 at the needed precision;
* `P_cyc` for the crystal cycle at finite Δ;
* anything about the identification of the project's EE integrand with CLS-D's spinning functional beyond what REPRO_ENGINE checked
  (shared definition, independent point-wise implementation);
* anything for J ≥ 4, the `+−` channels, or the Result D window.

## 6. Wording for the letter and the supplement
**Letter, main text** (replacing "That this equals the integral over the physical saddle-point cycle (H_RS) is derived given one premise
P_cyc … [NE + D-sketch]."):
> `s_d` comes from an analytic `1/Δ` expansion integrated over the real section. By its definition the functional's action on a heavy
> block is an integral over the causal region `4 > 5 > 3 > 0 > 1 > 2` of real Lorentzian points; in the gauge used this region is
> exactly the real section cut at the light cones of `x₅` (`D′`). On `D′` the expansion holds to all orders with a remainder
> `O(Δ^N e^{−0.336Δ})`, so the identity of the expansion with the functional (`H_RS`) is derived (one step at the level of a sketch), and a
> direct finite-Δ integration over `D′` reproduces it at the binding state. The saddle-point cycle is used only as a numerical method;
> its homology class is no longer a premise.

**Letter, table row** `H_RS`: "Q's corrections, derived on the physical chain `D′` [D]" (instead of "[D given P_cyc], preliminary").

**Letter, outlook**: drop "a proof of `P_cyc`"; if a cycle item is kept: "an independent precision check of the saddle-point numerics
below ν ≈ 10".

**Supplement, `H_RS` paragraph, item (iii)** (replacing "What remains: P_cyc."):
> (iii) *Which cycle.* CLS-D define the functional as an integral over real points in the causal order `4 > 5 > 3 > 0 > 1 > 2`
> (Lorentzian shadow representation of the block). In the gauge `x₃ = −x₁ = y`, `x₄ = −x₂ = x/(−x²)`, `x₅ = e` this region is
> `D′` [P; 2·10⁵ random configurations, no mismatch], every non-integer-power factor is positive on it, and the integral converges
> absolutely. So the physical integral is `∫_{D′}` up to a ν-dependent constant that cancels in `ρ` [D]; for the scalar this is checked
> end to end against the exact Mellin skeleton (≤ 10⁻⁹ at Δ = 300, 1000 through the expansion; ≤ 2·10⁻⁸ directly at Δ = 20, 30, 50
> for ν ≤ 5, where 3–80 % of the weight of `D′` lies outside the convergence region of the expansion) [DV]. With (ii), `H_RS` is [D]. The
> saddle-point cycle is a steepest-descent representation of `D′`; its class identity with `D′` (`P_cyc`) is [P] in the Regge-scaled
> system, open at finite Δ (a Stokes-wall question), and is used only by the saddle-point numerics. A direct EE test of those numerics
> is consistent at every tested point (Δ = 30–100, ν = 3–16), with a sensitivity within a factor ≈ 5 of the failure threshold above the
> fold and none below it.

**Supplement, hypothesis table**: `H_RS` → "[D] on the physical chain `D′` (was [D given P_cyc])"; `P_cyc` → "not a premise of `H_RS`
(the physical chain is `D′` by definition); as a statement about the saddle-point numerics: [P] Regge-scaled, [I] at finite Δ; direct EE test [NE, weak] above the fold, [I] below".

**Supplement, reproduction paragraph**: "Both routes integrate over the real section, so P_cyc is not tested" → "Both routes integrate
over the real section, which is the physical chain `D′`".

**UNIFIED_THEOREM_v2 §5.2, row H_RS** (suggested): "**[D]** on the physical chain `D′` (CLS-D's causal domain, PCYC_CLOSE §1); numerically
[DV] at the binding state (REPRO_ENGINE: direct `∫_{D′}` = engine `K = 5` to 3·10⁻⁷). `P_cyc` is not a premise; the saddle-point
numerics agree with `D′` for EE wherever tested (Δ = 30–100, ν = 3–16; sensitive only above the fold, within a factor ≈ 5 of
the failure threshold)."

## 7. What remains
1. **Nothing for `H_RS`'s cycle question.** What `H_RS` still inherits is RD_CLOSE's [D-sketch] sub-step (convergence radius of the
   `a`-series of the `x₀`-dependent slow factors on the polydisc of raw radius 1/3) and the definition-level identification of the
   project's EE integrand with CLS-D's spinning functional (REPRO_ENGINE: shared definition, independent implementation).
2. **`P_cyc` for the crystal cycle at finite Δ** (§2.3): a Stokes-wall analysis of the flow function in `(a, ν)`, or a complete census
   with intersection numbers. §4 shows the latter is far away (≥ 10³ distinct critical points per (Δ, ν), not saturated).
3. **Precision below the fold.** A crystal-cycle EE number at 10⁻³ below ν ≈ 10 needs either a much cheaper thimble integrator
   (exact-gradient flows à la `rd_core`, deterministic quadrature on the thimble) or a different representative of the crystal class.
   It is a check of the saddle-point numerics only.
4. **Δ = 30, ν ≥ 10 for EE on `D′`**: the product rule does not resolve the cancellation there; the precise n3 crystal values at
   Δ = 30 (ν = 12, 16: 0.3–0.6 %) wait for a finer `D′` grid (≈ 2× nodes in `r_x, r_y, φ, α`; ≈ 5 CPU-h).
5. **Record corrections** suggested to the other notes (none edited): RD_CLOSE §0/§4.3 "no critical point at raw radius 0.1–39" →
   withdrawn (§4); RD_CLOSE §6.1 "identify the CLS-D functional's defining contour" → done (§1); the letter/supplement/UT wording of §6.

## 8. Reproduce (`ads4/code/pcyc/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)
```
python3 pcyc_region.py 200000              # 1 s     -> cache/pcyc_region.txt   (physical domain = D'; signs on D')
./dimscan.sh 30 S                          # 4 min   -> cache/dimscan_30_S.txt  (per-dimension convergence, scalar, Delta = 30)
./dq_prod.sh 30; ./dq_prod.sh 100          # 2 x 30 min -> cache/dq_D{30,100}_{S,E}_<grid>.npz  (D' product rule, nu = 1..16)
python3 pcyc_exact.py 30 100               # seconds -> cache/pcyc_exact_D30_D100.txt  (D' vs exact skeleton, end to end)
./exact_more.sh                            # Delta = 20, 50 likewise -> cache/pcyc_exact_D{20,50}_D100.txt
./cr_queue.sh; ./cr_seed2.sh; ./d50_queue.sh   # crystal runs: n3 ~1-3 min, n4 ~12 min each (+ pcyc_cv.py, ~1-2 min) -> cache/cr_*.pkl
python3 pcyc_farw.py                       # 5 min   -> cache/pcyc_farw.txt (far share of |Phi| on D')
python3 pcyc_crys.py 100 16 11 8           # single run;  env PC_SEED, PC_TAG, PC_NUSPLIT (9.9 forces n3 at nu = 10), PC_STORE
python3 pcyc_cv.py cache/cr_D100_nu16.pkl 2   # control variates -> appended to cache/pcyc_cv.txt
python3 pcyc_compare.py; python3 mk_table.py   # -> cache/pcyc_compare.txt, cache/compare_table.md (the T table of section 3.4)
./census_queue2.sh; ./census_big.sh        # -> cache/census_D<Delta>_nu<nu>_s<seed>.txt (+ .pkl); ~75 s per 4000 starts
python3 census_table.py                    # -> the table of section 4
python3 crsum.py cache/cr_D100_nu4_cvtest.pkl # per-thimble diagnostics
```
