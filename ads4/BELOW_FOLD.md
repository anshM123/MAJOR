# Below the fold: a saddle-free representation of the crystal integral, what it reaches, and why `K` is still not measured

*21 Nov 2026, below-fold agent. Unrefereed, not committed. No existing document or code file was edited; nothing of the `oo_*` agent was
touched or run. New code: `ads4/code/spinfunc/bf_core.py`, `bf_exact.py`, `bf_q.py`, `bf_qan.py`, `bf_validate.py`, `bf_qrun.sh`. Caches
`spinfunc/cache/bf_*` (git-ignored). Read-only imports: `n1_pw` (point-wise `l1_crystal`, `l5_spin`, `l6_spin_action`, `l9_points`),
`w0_ee_common`, `n3_exact`, `n4_exact` (record builder, CHMRS/Mellin `ExactState`), `resultq/stage4_hj`, Stage-6/7 weight caches.
Compute ≈ 1.5 CPU-h, ≤ 2 cores, < 1.5 GB.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]** open.

---

## 0. Answer in brief

1. **Structure of the integrand at small ν (new) [DV].** In the Regge scaling `u, v, x⊥, y⊥ = a·(O(1))`, `a = 1/Δ`:
   * every one of the seven directions is `ν`-soft (Hessian of `λF` in scaled units `= (iν/2)·Hess log A + O(ν²/Δ)`, singular values
     `(3.24, 2, 1.24, 1, 1, 1, 1)·ν/2`). The "two directions `∝ Δ`" of H0 §24.2 are the raw units of `x⊥, y⊥`. **At small ν there is no hard
     direction to Laplace-expand**; the only large parameter is Δ, and it enters as a *regular* expansion;
   * `m₁₂ = m₃₄` identically, so `B = B_x B_y` factorises exactly [P]:
     `B^{Δ/2} = [(1−s_x+q_x)/(1+s_x+q_x)]^{Δ/2} [(1−s_y+q_y)/(1+s_y+q_y)]^{Δ/2} → e^{−(X₁+X₂+Y₁+Y₂)}(1 + O(1/Δ²))` (`s = u+v`, `q = uv − p²`);
   * `A` is homogeneous of degree 0 at leading order, `A₀ = 4 (X₂²/q_x)(Y₁/Y₂) cos²φ / sin²(ψ + n₀φ)` (code agrees to `10⁻⁷` at Δ = 300);
   * with `y⊥ = aρ sin φ`, `x̃₀ = cos φ cos ψ`, `y₀ = n₀aρ sin ψ` (`ρ² = Y₁Y₂`) the delta-function Jacobian times `dy⊥ dx̃₀` is exactly `dφ dψ`
     (the celestial sphere, both `n₀`); HLJ's fibre Beta function `B(½, ¼ − iν/2)` is the `ψ`-integral, which decouples at leading order.
   * The scalar spin-J dependence is `P_J(t)` with `t ≥ 1` on the real section (`t − 1 = O(1/Δ²)`) [DV, ratio to the code's `Q_J` form `1 ± 10⁻¹⁴`].
2. **Method (task 1): the real future-cone section `D` (`bf_core.py`) [D + NE].** `X₁, X₂, Y₁, Y₂ > 0`, `|x⊥| < √(u_x v_x)`, `(φ, ψ)` on
   the sphere: all integrand factors are real and `A > 0`, so `ν` enters only as the phase `A^{iν/2}`; the integral is absolutely
   convergent for the scalar and **has no saddle, hence no fold**. Randomised QMC with the exact separable leading-Regge density as
   importance density; state and `Δ_ref = 2·10⁴` on common points in the scaled variables (their integrands differ by `1 + O(ν³/Δ²)`).
3. **Validation (task 2) — scalar: passed; spinning and graviton: failed [NE].**
   * **Scalar J = 0** (exact CHMRS/Mellin skeleton, state/ref): `|m²ε| ≤ 0.2–2` c₂ on `ν ≤ 2` and within errors (`±8 … ±30`) to `ν = 3` at
     Δ = 100 and 300; consistent within `±300–1000` to `ν = 25` (Δ = 300; at Δ = 100 two nodes `ν ≥ 19` at 2.6–3.5σ); Δ = 1000: within errors.
     **This is the first evaluation of the crystal integral at `ν < 5`, down to `ν = 0.02`, to 0.1–1 c₂.**
   * **Precision wall: the `e^{−πν/2}` cancellation of the radial phases `X^{±iν/2}`.** QMC errors (16 × 2¹⁷ points, Δ = 300):
     0.3 / 1 / 5 / 29 / 340 / 800 c₂ at `ν = 0 / 1 / 2 / 3 / 4 / 5`. Rotating the radial rays towards the crystal does not help
     (§2.3). So `D` is a small-`ν` method (`ν ≲ 3.5`), complementary to n3 (`ν ≳ 10`); the window `3.5 ≲ ν ≲ 8` stays open.
   * **Spin `J > 0`: correct at `ν ≲ 2.6`, wrong beyond.** J = 4 (Δ = 100) and J = 24 (Δ = 336) reproduce the exact state/ref ratio at
     `ν ≤ 2.6` (`≤ 1σ`, `±14` / `±330`) and then fail by 6–65σ, growing like `J(J+1)ν⁴`; J = 2 drifts at `ν ≥ 14` (3–9σ). Cause [D + NE]:
     the boost growth `P_J(t) ~ e^{J|η_x − η_y|}` makes the rapidity integrals heavy-tailed and reach the non-Regge far region where the
     real section is not equivalent to the crystal cycle. Validations (a) (n3/n4/K₂ at `ν ≥ 5`) and (c) (large J → leading Regge) are
     therefore **not achievable with `D`**.
   * **Graviton EE: insufficient [NE].** The Regge-limit EE `k = 4` integral on `D` changes sign at `ν = 1.05 ± 0.05`, the first zero
     1.0685 of Stage 4's `N̂₄^{EE}` (weak check of the spin-2 structures). But EE has the same boost growth (EE/scalar `∝ (Y₁/Y₂)^{≈1.7}`) and
     a `cos^{−4}φ` light-cone enhancement; after fixing the latter the state/ref ratio is still only known to `0.1–0.9 %`, i.e.
     **`±150 … ±2700` c₂ per harmonic, against the 3 c₂ needed.**
4. **Result Q (task 3): `K` not measured; δ verdict unchanged [I].** No EE `k = 4` image beyond the model was obtained at the required
   precision at any `ν`. Stage 7 stands: δ = 0.2 positive for `|K| ≤ 1` (min +0.068) in every gauge, δ = 0.05 not (−0.051); quote δ = 0.2.
5. **Result D (task 4): not attempted beyond the diagnosis [I].** `HYPOTHESES_CLOSE2.md` does not exist. The `ν < 8.3` window at
   `(336, 28/32)` and the point `(1150, 150)` need `J ≥ 28` at `ν ≳ 3`, exactly where `D` fails (item 3).
6. **What would close it (§5) [I].** The same structure gives a route without the two obstructions: expand in `1/Δ` at fixed scaled
   variables (each order a polynomial in the radial variables times the leading power law and `e^{−S}`), do the four radial integrals
   analytically (Gamma functions: they absorb both the `e^{−πν/2}` cancellation and the boost growth by continuation, exactly as the
   crystal contour does), and the three angular ones `(s, φ, ψ)` numerically (power-law, no cancellation). The graviton needs its
   polarisation prefactor as monomials in the radial variables.

---

## 1. Structure of the crystal integrand in the Regge scaling

### 1.1 Scaling and the absence of a hard direction [DV]
`l1_crystal` variables `z = (x̃₀, u_x, v_x, u_y, v_y, x⊥, y⊥)`, crystal `(0, −ia, ia, ia, −ia, 0, 0)`, `a = (Δ − √(Δ² − ν²))/ν ≈ ν/2Δ`.
Writing `u_x = aX₁` etc. (all four light-cone components and both transverse ones scaled by `a`):
* the Hessian of `log A` at the crystal is an `O(1)` matrix (singular values above); the Hessian of `log B` is `O(a²)`. So in natural units
  every direction has curvature `∝ ν`, and the Laplace series is an expansion in `1/ν` in all seven directions — this is why `K₂/K₃`
  degrade together and the saddle folds at a `Δ`- and `J`-independent `ν* ≈ 11.3` (EE_DG300_SMALLNU §1, HLJ §1.3);
* the raw-coordinate eigenvalues `~ 2Δ²/ν` and `±2Δ` quoted as "hard" come from `x⊥` and from the `(x̃₀, y⊥)` mixing `−iν/a = −2iΔ`:
  units, not hardness.
* Consequence for Stage 7's suggestion "hard-direction Laplace × exact soft block": at small `ν` the soft block is the whole integral.

### 1.2 Exact pieces [P / DV] (`bf_core`, scratch symbolic checks)
| piece | form | check |
|---|---|---|
| `B` | `m₁₂ = m₃₄` exactly, `B = (m₄₅/m₂₅)(m₃₅/m₁₅)`, `m₄₅/m₂₅ = (1 − s_x + q_x)/(1 + s_x + q_x)` (raw `s = u+v`, `q = uv − p²`), same for `y` | sympy; `(Δ/2) log B = −(X₁+X₂+Y₁+Y₂) + O(Δ⁻²)` at `a = 1/Δ` to `10⁻⁶` |
| `A` (leading) | `4 (X₂²/q_x)(Y₁/Y₂) cos²φ/sin²(ψ + n₀φ)`, `q_x = X₁X₂ − X⊥²` | code `/` formula `= 1 ± 10⁻⁷` at Δ = 300; exact `A = 4` at the crystal |
| measure | `(uv)^{½}(uv − y⊥²)^{−½}(uv(1−x̃₀²) − y⊥²)^{−½} dy⊥ dx̃₀ = dφ dψ` | algebra; `y₀/(aρ sin ψ) = 1` to round-off |
| leading density (scalar, J = 0) | `q_x^{7/4−iν/2} X₂^{−½+iν} (X₁Y₂ + X₂Y₁ − 2X⊥Y⊥) Y₁^{9/4+iν/2} Y₂^{7/4−iν/2} cos^{9/2+iν}φ |sin(ψ+n₀φ)|^{−½−iν} e^{−S}` | matches the code at `q = 0` to `10⁻⁵`; the `x–y` coupling `−2x·y/a² > 0` on `D` |
| spin J | `exp(G_J − G₀) = P_J(t)`, `t = Z·W/√(Z²W²) ≥ 1` on `D` | `1 ± 10⁻¹⁴` (J = 4, Δ = 100, 2·10⁴) |
| first-order corrections | `O(a)` only from the 3-point factors `m₁₅, m₂₅, m₃₅, m₄₅`; `FP`, `m₂₀`, `m₄₀`, `m₁₂` start at `O(a²)`; `A` has no `O(a), O(a²)` term in the expansion computed | sympy factor-by-factor |

### 1.3 Why the real section cannot fold
On `D` every factor is real and `A > 0`, so `|integrand|` does not depend on `ν`: changing `ν` rotates phases only. There is no saddle
whose existence depends on `ν`, hence no fold, no Stokes switch, no walls (diagnostic: no factor ratio winds; `Im G = 0` on all points).
At `ν → 0` the scalar integrand is positive. [D + NE]

## 2. Method (`bf_core.py`)

### 2.1 Section and sampler
* `D`: `X₁, X₂ ~ Γ(4), Γ(3)`; `x⊥ = s√(X₁X₂)` with `s ~ (1−s²)^{7/4}`; `Y₁, Y₂ ~ Γ(4), Γ(3.5)`; `φ ~ Beta(c+1, c+1)` (`c = 4.5` scalar,
  `0.5` graviton, §3.3); `α = ψ + n₀φ ~ π·Beta(½, ½)`, `ψ = α − n₀φ mod π`; both sheets `n₀ = ±1` (equal on `D` to `10⁻⁵`, as on the
  crystal). Scrambled Sobol, 8–16 replicas; error = replica spread.
* The integrand is the project's (`n1_pw` point-wise `l1_crystal.build`; graviton `w0_ee_common` EE builder; `k = 4` via Lemma B
  `R = ((v′−u′)² − 1)/(4u′v′)`); logs continued from a fixed interior point (all factors real, positive ratios).
* State and reference at the same scaled points (`a = 1/Δ` for each): the state/ref ratio is a CRN estimate.
* Cost: ≈ 10 s per 2¹⁷ points for scalar + EE, both Δ, all 150 `ν` nodes at once (the `ν`-dependence is `A^{iν/2}` on stored `log A`).

### 2.2 Variants tried and rejected [NE]
* **Rotated radial rays** (`X₁ → X₁e^{−iθ}`, `X₂ → X₂e^{iθ}`, `Y` alike; `θ → π/2` is the crystal direction), meant to remove the radial
  phase cancellation: for the scalar J = 0, `θ = 0.3` agrees at `ν ≤ 3` but shows 3–5σ systematics at `ν ≥ 9`; `θ ≥ 0.6` is worse at every
  `ν`; for J = 4 it moves the small-`ν` values by `+40ν²`. Heavy tails make the QMC errors unreliable off the real section. Kept `θ = 0`.
* **Rapidity sampler** (`X = re^{±η}`, logistic `η`): reaches the far region (`aX ~ 1`) where the real section meets zeros/branch loci of
  the exact `B` and the 3-point factors; garbage. Not used.

## 3. Validation (`bf_validate.py` → `cache/bf_validate.txt`; `bf_q.py scal` on the production points)

### 3.1 Scalar, exact Mellin skeleton (state / `Δ_ref = 2·10⁴`, same J; `m²ε = m²[log(R_D/R_exact) − const]`, const fitted on `ν < 1.4`)
| state | `ν ≤ 1` | `ν ≈ 2` | `ν ≈ 2.6–3` | `ν = 4–14` | `ν = 19–25` | grade |
|---|---|---|---|---|---|---|
| Δ = 100, J = 0 (2¹⁵·8) | `−0.2 … +0.2 ± 0.1–0.3` | `−2.9 ± 2.4` | `+11 ± 8` | `|·| ≤ 670`, all `≤ 1.1σ` | `−802 ± 228`, `−908 ± 344` (3.5σ, 2.6σ) | [NE] |
| Δ = 300, J = 0 (2¹⁷·16, production) | `≤ 0.9 ± 0.5` | `−6 ± 4` | `−16 ± 29` (3.0) | `≤ 1σ` (errors 100–1700) | `1.2–1.7σ` | [NE] |
| Δ = 1000, J = 0 (production) | `−6 ± 4`, `+12 ± 6`, `−15 ± 12` (one 2.2σ; `10⁻⁵` relative) | `−25 ± 29` | `+41 ± 194` | `≤ 1σ` (errors 700–16 000) | — | [NE] |
| Δ = 315, J = 2 / Δ = 1050, J = 2 (production) | `≤ 1.4 ± 1.5` / `+14 ± 4` at 0.62 (3.2σ), else `≤ 1.8σ` | `−3 ± 4` / `−27 ± 25` | `−21 ± 31` / `+41 ± 206` | `≤ 1σ` | — | [NE] |
| Δ = 100, J = 2 | `≤ 0.4 ± 0.3` | `−2.4 ± 2.2` | `+9 ± 9` | `≤ 2σ` | `+1384 ± 252`, `+3258 ± 353` (5.5σ, 9.2σ) | [NE] to `ν ≈ 12` |
| Δ = 100, J = 4 | `≤ 0.3 ± 0.4` | `−1.3 ± 2.5` | `+2 ± 14`; `+39 ± 44` (3.5) | `+3020 ± 477` (5.45), `+3683 ± 368` (13.8) | `+7983 ± 319`, `+25 317 ± 388` | [NE] `ν ≤ 2.6`; **fails** beyond |
| Δ = 336, J = 24 | `≤ 7 ± 11` | `+40 ± 66` | `−263 ± 334` | `+31 000 … +142 000` (7–35σ) | `+2·10⁵` (51σ) | [NE] `ν ≤ 2.6`; **fails** beyond |

(The exact ratio itself varies by `10³–10⁴` c₂ over these ranges, e.g. 1700 c₂ between `ν = 0` and 1.3 at (336, 24).)
* Near `ν ≈ 5.5–6.6` the comparison is ill-conditioned for every method: both `X(ν)` have the zero of `g + g̃` (`ν₀ = 5.979`).
* **Reading.** For non-spinning (J = 0) blocks `D` is equivalent to the crystal cycle to all orders tested (the Δ-dependence of the
  ratio, which is all of the beyond-LR content, is right at the `10⁻⁵–10⁻⁶` relative level at small `ν`). For J > 0 the agreement at
  `ν ≲ 2.6` (the `ν²` term, i.e. the leading `β²`-structure) is real, but higher orders in `ν` are not reproduced.

### 3.2 Why J > 0 fails [D + NE]
* `P_J(t)` with `t = cosh(η_x − η_y)`-like (`t ≥ 1` on `D`) grows like `e^{J|Δη|}`; the scalar density decays like `e^{−5.5|η|}` per pair.
  For J = 4 the rapidity integrand decays only like `e^{−1.5|η|}`; the Gamma sampler's tails give infinite variance, and the large-`|η|`
  region reaches `aX ~ 1`, where the real section is not the crystal cycle (zeros of the exact `B` and of `m₂₅`, branch loci).
* On the crystal cycle the rapidity is imaginary (`X₁/X₂ = −1`), `P_J` is bounded, and the integral is the analytic continuation. A
  real section cannot reproduce that for growing spin. [D-sketch for the far-region mechanism]

### 3.3 Graviton EE [NE]
* **Light-cone edge.** EE/scalar `∝ cos^{−4}φ` at the `y` light cone (median ratio ×5 at `cos φ ≈ 0.3`, ×400 at 0.15): the EE density is
  `∝ cos^{1/2}φ`; sampler changed to `c = 0.5` (`bf_q.CPHI_G`). This removed the dominant spikes.
* **Boost growth.** EE/scalar grows like `(Y₁/Y₂)^{≈1.7}` (median ×4.5 at `Y₂/Y₁ ≈ 0.1`, ×40 below 0.05): spin-2 polarisations carry boost
  weight, the same disease as §3.2.
* **Zero of `N̂₄^{EE}`.** The Regge-limit (Δ = 2·10⁴) EE `k = 4` integral on `D` changes sign between `ν = 1.0` and `1.1` (Stage 4:
  1.0685); the EE/scalar ratio follows the sign pattern of `N̂₄^{EE}/N̂₄^{sc}` but not its magnitude (factor ≈ 2 at `ν = 2–3`, errors
  20–40 %) [NE, weak].
* **Precision** (`cache/bf_ee_precision.txt`, 16 × 2¹⁷ points): `I_ref` itself relative error 0.24–1.1; CRN ratio `ρ = I_st/I_ref`
  `0.16–0.87 %` (Δ ≈ 300) and `0.04–0.25 %` (Δ ≈ 1000), i.e. `m²·δρ = 145–835` (300) and `465–2660` c₂ (1000), per harmonic, for
  `++` J = 0 and J = 2 alike. The scalar on the same points: `0–31` (300), `3–206` (1000).

## 4. Result Q (`bf_q.py`, `bf_qan.py`) [I]
* Computed: `++` J = 0 gap states (Δ = 300, 1000) and `++` J = 2 at `1.05Λ` (Δ = 314.93, 1049.93), EE `k = 4` and `k = 2` and the scalar,
  `ν = 0.02 … 5.98` on the `dν = 0.04` de-shift grid (`cache/bf_q_D*.pkl`).
* `bf_qan.py` implements the full chain (normalisation `N(ν) = I_ref/T_ref` through the zero of `N̂₄`, `Γ_K = Γ_exact − 0.75 Σ_j R_j(ν−ij)²`,
  `a_h` class removed, projections on H0, H2, H3, H4, SB, XP). With the §3.3 precision the projections are `s_d = O(10²–10³) ± O(10²–10³)`
  against a tolerance `|s_d| ≲ 0.5` (`3` c₂ × `A[1] ≈ 0.17`): **no information**. The `+−` J = 4, 5 states need the 2×2 helicity bases
  and are J > 0; not run.
* **Verdict: unchanged from Stage 7.** `K` [I]. δ = 0.2 recommended (min `+0.068` at `|K| = 1`, every gauge); δ = 0.05 not established.
  No `FF/N` value is changed by this note.

## 5. What would close the gap [I]
1. **Analytic-radial Regge expansion (recommended).** Write the exact integrand at `a = 1/Δ` as
   `H_ν(w) e^{−S}[1 + a f₁(w) + a² f₂(w) + …]` (§1.2; the `f_k` are polynomials in the radial variables `X₁, X₂, Y₁, Y₂` with angular
   coefficients, plus `P_J(t)` = polynomial in `a²(t−1)/a²`). Then for each angular point `(s, φ, ψ)` the radial integrals are sums of
   `Γ(p₁+1)Γ(p₂+1)Γ(p₃+1)Γ(p₄+1)` with `ν`-dependent complex `p_i`:
   * the `e^{−πν/2}` cancellation is then exact Gamma-function arithmetic (no QMC loss at any `ν`);
   * boost growth only shifts the `p_i` (continuation, as on the crystal), so J > 0 and the graviton are covered;
   * the angular integrals are power-law (`cos^{…}φ`, `|sin α|^{−½−iν}`, `(1−s²)^{7/4}`) and smooth: Gauss–Jacobi rules;
   * the scalar version can be validated against the exact skeleton at every `ν` (Δ-dependence = the `f₂` moment), the graviton
     against n3 at `ν ≥ 12` and Stage 4's `R_j` at leading order;
   * needed ingredient: the EE `⟨TTO⟩`/polarisation prefactor as monomials in the radial variables (symbolic `l5_spin.tto_samples`, or a
     numerical monomial fit per angular point).
   Truncation: relative `O(ν⁶/Δ⁴)` after the `a⁴` term, i.e. `≤ 1` c₂ for `ν ≤ 8` at Δ = 300 (order `a⁶` for `ν ≤ 18`).
2. For the scalar alone, `D` already gives the exact all-orders image at `ν ≲ 3` to `≤ 1–30` c₂ (Δ = 300) with ≈ 3 min of compute.

## 6. Graded claims
| claim | grade |
|---|---|
| all seven directions `ν`-soft in the Regge scaling; no hard direction at small `ν`; H0's "2 directions ∝ Δ" are units | [DV] |
| `m₁₂ = m₃₄`, `B = B_xB_y`, `B^{Δ/2} → e^{−S}(1 + O(Δ⁻²))` | [P] (sympy) + [DV] |
| leading `A₀`, flat celestial measure, leading scalar density, `G_J − G₀ = log P_J(t)`, `t ≥ 1` on `D` | [DV] |
| on `D` all factors real, `A > 0`: no saddle, no fold, `ν` only a phase | [D] + [NE] |
| `D` reproduces the exact scalar J = 0 state/ref ratio (Δ = 100, 300, 1000) to `≤ 0.2–2` c₂ at `ν ≤ 2`, within errors to `ν = 3` and (`±300–1000`) to `ν ≈ 17` | [NE] |
| `D` precision falls like the `e^{−πν/2}` radial cancellation: `ν ≲ 3.5` usable | [NE] |
| J > 0: correct to `ν ≈ 2.6`, fails beyond (`∝ J(J+1)ν⁴`); mechanism boost growth + far region | [NE] (failure) / [D-sketch] (mechanism) |
| EE on `D`: zero of `N̂₄^{EE}` at 1.05 ± 0.05 | [NE, weak] |
| EE `k = 4` beyond-model image of `++` J = 0, 2 to `≤ 3` c₂ | **not achieved** (errors `10²–10³` c₂) → `K` [I] |
| Result Q δ verdict | unchanged (Stage 7: δ = 0.2 [M given K ≲ 1]; δ = 0.05 [I]) |
| Result D below-fold points `(336, 28/32)`, `(1150, 150)` | not computed [I] |
| analytic-radial Regge expansion closes both obstructions | [I] (design) |

## 7. Files (`ads4/code/spinfunc/`)
| file | content | run |
|---|---|---|
| `bf_core.py` | section `D`, samplers (`Sample`, `SampleR` rejected), scalar/graviton point-wise parts, `t_ZW`, cross-ratios and Lemma-B `R`, `combine` (all `ν` at once) | library |
| `bf_exact.py` | exact Mellin records (via `n4_exact.record`) for new states → `cache/bf_exact.pkl` | seconds |
| `bf_validate.py` | scalar validation table §3.1 (Δ = 100 J = 0, 2, 4; Δ = 300 J = 0; (336, 24)) → `cache/bf_validate.txt` | 25 s |
| `bf_q.py`, `bf_qrun.sh` | production EE `k = 4/2` + scalar sums, state/ref, 16 × 2¹⁷ → `cache/bf_q_D*_J*.pkl`; `bf_q.py scal D J` = scalar check on those points | 3 min/state |
| `bf_qan.py` | `Γ_K` chain and gauge projections (H0, H2, H3, H4, SB, XP) → `cache/bf_qan.txt` (noise-dominated, §4) | 1 min |
| `cache/bf_ee_precision.txt` | EE vs scalar precision table §3.3 | — |
