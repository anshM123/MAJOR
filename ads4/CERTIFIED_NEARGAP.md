# Certified near-gap positivity of the D-safe heavy action (scalar-skeleton model), Δ_g = 100

> **Withdrawn functional (18 Nov 2026, later).** Theorems C and C′ below certify the **B′-uplift D-safe functional `Ω^{NF,SC}`
> (k = 2…7)**. That functional is withdrawn (`H0_SPINNING_FUNCTIONAL.md` §9–§11).
> * The route-A′ certificate is in **`CERTIFIED_NEARGAP_A.md`**. Route A′ is `∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`, with the stripped density and the
>   `η`-inversion.
> * Everything here about the certified moments (`cert_moments.py`) and the quadrature method (`cert_quad.py`) is reused there unchanged.

Status: 18 Nov 2026, derivation/numerics agent. **Unrefereed.** Code: `ads4/code/exactskel/cert_moments.py`, `cert_action.py`,
`cert_summary.py` (test: `cert_t1.py`). Caches and logs: `exactskel/cache/cert_100_Df*.{pkl,log}` and `cache/cert_100_summary.txt`.
The other files of `exactskel/` are used read-only; `omegau/` is used read-only.

Grades:
* **[certified]** a computer-assisted statement. Every floating-point step is an Arb ball (python-flint 0.9.0, 192 bits). Every
  truncation is replaced by a proven bound, and the inequalities used are checked in exact rational arithmetic.
* **[D]** analytic input, derived or quoted, not machine-checked.
* **[DV]** derived and numerically validated.
* **[NE]** numerical evidence.
* **[M]** a modelling definition.
* **[I]** an assumption.

## 0. What is certified

**Theorem C (computer-assisted).** Take `Δ_g = 100`, `Λ = √(98·99)`, and identical-scalar skeleton `Δφ ∈ {2, 3}`. Consider every state
in the table of §3. Its `m/Λ` values lie in a grid `G ⊂ [1.01, 2]` (`Δφ = 2`) or `G ⊂ [1.01, 3]` (`Δφ = 3`), with `J ∈ {0, 2, 4, …, 10}`
for `++` and `J ∈ {4, …, 10}` for `+−`, and `τ ≥ Δ_g`. For each such state the quantity `F_h^{(64)}(m, J)` of §1 satisfies
`F_h^{(64)}(m, J) > 0`. The certified lower bounds are listed in §3.

**What `F_h^{(64)}` is.** It is the D-safe (NF02 / `Ω^{NF,SC}`) heavy action of `OMEGA_NF.md`, with three specifications:
* exact scalar-skeleton moments (§1);
* the 64-node Gauss–Legendre `p`-rule of OMEGA_NF taken as the definition of the `p`-integral;
* odd-`k` moments defined by quartic interpolation in `k` from `k = 2, 4, 6, 8, 10`.

**What is certified about what.**

| object | status |
|---|---|
| the scalar-skeleton moments `m^{2k}C_{k,ν}[Δ, J]` at the nodes (even `k`), given the Mellin representation | [certified] enclosures |
| the D-safe action of the model as defined above | [certified] |
| positivity of the actual graviton (⟨TTTT⟩) moments | **not** addressed: skeleton → graviton transfer (H4), reading B (H2), odd-`k` two-pair moments (H3) stay [I] |

**Theorem C′ (computer-assisted, second pass: the exact `p`-integral).** Let `F_h(m, J) = Σ_k u^{k−2} ∫_0^1 φ_k(p) Mod_k^h(p) dp` be the
**exact** `p`-integral, with the same model as above. Then `F_h(m, J) > 0` for:
* every state of Theorem C except `++, J = 0, m = 1.01Λ`: 204 / 205 cases for `Δφ = 3` and 178 / 179 for `Δφ = 2`;
* the high-`J` states `J ∈ {11, 12, 13, 14, 16}`, `m/Λ ∈ {1.12, 1.15, 1.2, 1.3}` (`Δφ = 3`; §7.3).

The certified Gauss-error bound is ≤ 3·10⁻⁷ `N_s` for `m ≥ 1.03Λ`, and ≤ 0.43 `N_s` at `1.02Λ` (`J = 0`, where `F/N_s ≈ 1`). §7.1 has the
method; `cache/certq_100_summary.txt` has the numbers.

**Not achieved: positivity on continuous `m`-intervals.** It stays [NE]; §7.2 explains why and what it would cost.

## 1. The certified object

**State variables.** For a state with `Δ` (the double used by the scan; `α = Δ − J − 2`, `ᾱ = Δ + J − 1`, `m² = αᾱ`, `cosh β = (α²+ᾱ²)/(2m²)`):
* `u = Λ²/m²`;
* `p_i, w_i`: the 64 Gauss nodes and weights on `[0, 1]`, as doubles;
* `ν_i = Λ p_i`, `x_i = ν_i/m`, `y = x²`.

**The action.**

`F_h^{(64)} = Σ_{k=2}^{7} u^{k−2} Σ_i w_i φ_k(p_i) Mod_k^h(p_i)`.

`φ_k` are the exact rational NF02 weights.

| `k` | `Mod_k^h` |
|---|---|
| even | `ex_k + E_k (D_h − P_J)` |
| odd | `σ_h [ω_k (D_h − P_J) + Pc + ω_k (ex_k − 2Pc)/E_k]` |

with
* `E_k = (2−y)/(1−y)^{1+k/2}` and `ω_k = (1−y)^{−(k+1)/2}`;
* `Pc = P_{−½+iν}(cosh β) − J₀((2J+1)x)`;
* `D_{++} = P_J(1−2y)` and `D_{+−} = d^J_{44}`;
* `σ_{++} = 1`, `σ_{+−} = −1`.

This is algebraically OMEGA_NF's model, `e_k D_h + 2ε_kσ^k Pc + σ^k env_k ρ_k` with `ρ_k = (ex_k − E_kP_J − 2Pc)/E_k`.

**The moments.**
* Even `k`: `ex_k = m^{2k}C_{k,ν}[Δ, J]`, the exact heavy density of the CHMRS functional on the `d = 3` block.
* Odd `k` [M]: `ex_k := E_k·Σ_l L_l(k) ex_l/E_l`, a Lagrange interpolant over `l = 2…10`. At `Δ_g = 100` it differs from the
  Gamma-ratio `k`-continuation by ≤ 2·10⁻³ `E_k` near the gap and ≤ 3·10⁻⁵ `E_k` for `m ≥ 1.3Λ` [NE, `t_kinterp.py`].

**Evaluation.** Every function is evaluated in Arb:
* Legendre and Jacobi polynomials by their recurrences;
* `P_{−½+iν}` by `acb.legendre_p`;
* `J₀` by `arb.bessel_j`.

## 2. Method for the moments (`cert_moments.py`)

**Representation [D].** The derivation and the floating-point checks are in `EXACT_SKELETON_NEARGAP.md` §7.
* CHMRS's Mellin formula (eqs. CfromBMellin and Bkt action of arXiv:2106.10274).
* Mack polynomials from Costa–Gonçalves–Penedones (arXiv:1209.4355, app. A), mapped to CHMRS conventions as in `mack.py`. They are
  **[DV]**: they agree with the CHMRS order-1/m⁸ series to 1e−10 – 1e−12 at small `ν/m`, for `J ≤ 24`.
* Barnes' first lemma: `∫dT/(4πi) E(T)e^{Tu}/γ_a² = ₂F₁(c, c̄; A; 1−e^{−2u})/Γ(A)`.
* Laplace representations of the Pochhammer and bracket factors.

Together these give
`I(ν) = Σ_{j=0}^{J} Σ_{m≥0} h^{(j)}_m(ν) S_{j,m}`, with `S_{j,m} = Σ_{n≥0} q_{n,j} V_{n,m}`.

**Certified evaluation.**
* **Direct part.** For `n < N` and `m < M`, every quantity is an Arb ball. The Mack polynomials are exact rationals.
  * Defaults: `N = 0.6Δ²/(Δφ−1)`, `M ≈ 3Λ²/(2Δ) + 6√Λ + 80`.
  * One `Δφ = 2` run used `N` halved at low `m`. This is still rigorous, just wider.
* **`n`-tail (`n ≥ N`), `m < M`.** The Laplace kernels of `Γ(z)/Γ(z+κ)` and `Γ(z)/(zΓ(z+κ))` give
  `|V_{n,m}| ≤ |W_n|/Γ(κ)·[|c₂| β(x_n, m+κ) + β(x_n, m+κ+1)(x_n+m+κ)/(2κ(x_n−1))]`. The second term uses Gauss's ₂F₁(1) sum. The bound
  is then completed by:
  * `|W_n| ≤ |W_N|(x_N/x_n)^P`. **Certified:** the polynomial `(x+b₁)(x+b₂)x^P − (x+a₁)(x+a₂)(x+1)^P` has non-negative Taylor coefficients
    at `x_N` (exact rationals).
  * `|q_{n,j}| ≤ Q̄_j (x_n/x_N)^{d_j}`. Here `Q̄_j = |q_j(x_N)| Π(x_N+|r|)/(x_N−|r|)` over rigorous root enclosures (flint), with a
    `Σ|coeff|` fallback.
  * `β(x_n, ·) ≤ β(x_N, ·)(x_N/x_n)^{(·)}`.
  * The sum over `n` is bounded by an integral.
* **`m`-tail (`m ≥ M`), all `n`.**
  * `|S_{j,m}| ≤ K_j β(x₀, m+κ)(1 + (x₀+m+κ)/(x₀−1))/Γ(κ)`, with `x₀ = τ/2`. `K_j` is a direct Arb sum plus its own `n`-tail bound.
  * It is summed explicitly against `|h^{(j)}_m(ν)|` for `M ≤ m < M₂`.
  * Beyond `M₂` it uses `|h^{(j)}_m| ≤ (4(m+J))^J g_m(Λ)`, where `g_m` is decreasing for `m ≥ M₂`. **Certified:** the ratio inequality
    `t_{m+1}/t_m ≤ (m/(m+1))²` holds via a non-negative-Taylor-coefficient check, so `Σ_{m≥M₂} t_m ≤ t_{M₂}(1+M₂)`.
  * `M₂` is doubled until the far tail is ≤ 1e−9 of the absolute `m`-sum.
* **Normalisation.** `log q_{Δ,J}` and `a_{Δφ}` are computed in Arb.

**Checks.**
* The balls contain the earlier floating-point values (`e1`, `skelexact3`) at every tested node (`cert_t1.py`).
* Typical radii relative to `E_k`: 1e−4 (`k = 2`, dominated by the `n`-tail) and 1e−8 (`k = 6`).
* The certified `F/N_s` agree with OMEGA_NF / `e7` (exact skeleton, floating point): for example `+−, J = 4, 1.04Λ`: `[0.24122, 0.24138]`
  against 0.24133.

## 3. Results

**Every case of the grid is certified: 384 of 384 (state, channel) cases, `F_h^{(64)} > 0`** (`cache/cert_100_summary.txt`).

**Grid.** `m/Λ ∈ {1.01, 1.02, …, 1.10, 1.12, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2.0}`, plus `{2.5, 3.0}` for `Δφ = 3` only.
`J ∈ {0, 2, 4, …, 10}` (`++`) and `{4, …, 10}` (`+−`), restricted to physical states (`τ ≥ Δ_g`).

| skeleton | states / cases | `m/Λ ∈ [1.01,1.05)` | `[1.05,1.1)` | `[1.1,1.2)` | `[1.2,1.5)` | `[1.5,2]` | `[2.5,3]` |
|---|---|---|---|---|---|---|---|
| `Δφ = 3` | 145 / 205 | ≥ +0.241 (`+−` J=4, 1.04) | ≥ +0.227 (`+−` J=5, 1.05) | ≥ +0.296 | ≥ +0.307 | ≥ +0.494 | ≥ +0.792 |
| `Δφ = 2` | 127 / 179 | ≥ +0.241 | ≥ +0.072 (`++` J=6, 1.06) | ≥ +0.100 | ≥ +0.301 | ≥ +0.489 | — |

Entries are the certified lower bounds of `F/N_s`. `N_s` is the floating-point component norm, used for normalisation only; the certified
statement is `F > 0`.

**Ball radii.**
* `Δφ = 3`: ≤ 4·10⁻² `N_s`; typically 10⁻³.
* `Δφ = 2`: up to 0.25 `N_s` at `1.05–1.2Λ`. Those states were run with a halved direct `n`-range; the bound is wider but still certified.

The lowest certified values (+0.07 to +0.10 `N_s`, `Δφ = 2`, `++, J = 6`) are wide balls around midpoints of about +0.3 `N_s`.

## 4. What is not certified

1. **The `p`-integral: now certified** (Theorem C′, §7.1), except at `++, J = 0, m = 1.01Λ`. There the ellipse is too thin, because the
   singularity sits at `p = m/Λ = 1.01`.
2. **Continuity in `m` between grid points [NE].** Adjacent grid points (spacing 0.01Λ near the gap) differ by ≤ ~0.05 `N_s`, against
   margins ≥ 0.07–0.2 `N_s`. There is no certified Lipschitz bound.
3. **Outside the grid region.**
   * `m < 1.01Λ`: the `++, J = 0` state at `m = Λ` has `F/N_s ≈ +0.997` [NE].
   * `J > 10` near the gap: certified only on the sample of §7.3. `J = 18, 20` at `1.2Λ` give balls too wide to certify (midpoints
  ≈ +1.0 `N_s`); all others are [NE] (`e7` sample).
   * `m > 3Λ` (`Δφ = 3`) or `m > 2Λ` (`Δφ = 2`): OMEGA_NF / OMEGA_U [M].
   * `Δ_g = 300` and `1000`: [NE] in `EXACT_SKELETON_NEARGAP.md`. The direct `N ∝ Δ²` makes these runs about 10× and 100× slower.
4. **Analytic inputs.** The Mellin representation, the Mack-polynomial formula and its convention mapping are [D]/[DV], not formalised.
5. **Odd-`k` definition [M] and the graviton-level hypotheses** H2, H3, H4 [I].
6. **Large-`b̂` region (H13).** Untouched beyond the scanned `J ≤ 10`.

## 5. Grades for `RESULT_D`

H13′ is upgraded, for `Δ_g = 100` and the finite state set of §3 plus the high-`J` sample, from [NE] to **[certified]**. This now covers
the **exact `p`-integral** of the scalar-skeleton model (Theorem C′), except for one state.

It remains [NE] for:
* `m` between grid points;
* `J > 10`;
* `m` beyond the grid;
* `Δ_g ≥ 300`.

It remains [I] for the graviton dictionary H2–H4.

## 6. Run record

**Runs.**
* `cert_action.py 100 3` (all of `G`) and `cert_action.py 100 2 <G ∩ [1.01, 2]>`, run in parallel on 2 cores.
* About 25–120 s per state; about 2.5 h of wall time in total.
* Restarts:
  * after improving the `m`-far-tail choice (adaptive `M₂`);
  * after introducing the root-based `Q̄_j` (`Δφ = 2`).
* Cases certified earlier with the looser bounds stay valid, since those bounds were rigorous too.
* Uncertified cases were recomputed with the final code. All cases in the final pickles are certified.

**Continuity diagnostic [NE].**
* The largest change of `F/N_s` between adjacent grid points, relative to the smaller certified lower bound of the pair, is 0.76
  (`Δφ = 3`) and 0.85 (`Δφ = 2`), in both cases `+−, J = 5`, `1.15 → 1.2Λ`.
* The largest absolute jump is 0.24 `N_s` (`+−, J = 8`, `1.08 → 1.09Λ`), where the margins are ≥ 0.6.
* So linear interpolation between grid points stays positive everywhere. That is not a proof.

**`Δ_g = 300` attempt.** A spot check (`cert_300_Df3.log`, `m = 1.05Λ`, `J = 4`) gave useless balls (`±10⁶ N_s`). At `ν ≈ 300` the `n`- and
`m`-tail bounds are too loose with the default `N` and `M`. Tightening them (larger `N`, an `n`-tail bound that uses the decay of
`h^{(j)}_m` in `m`) is the next step. It is not done.

## 7. Second pass: exact `p`-integral, continuity, high `J`

### 7.1 Certified Gauss error (`cert_quad.py`) [certified]

**The bound.** Trefethen, ATAP Thm 19.3, mapped to `[0, 1]`: if `f` is analytic in the Bernstein ellipse `E_ρ` and `|f| ≤ M` there, then
`|∫₀¹ f − Σ_{i=1}^{64} w_i f(p_i)| ≤ ½·(64/15)·M ρ^{−128}/(ρ²−1)`.

**Analyticity of the integrand.** `f(p) = Σ_k u^{k−2} φ_k(p) Mod_k(p)`.
* The ex-moments are entire in `ν`: they are absolutely convergent series of polynomials in `ν²`.
* `P_{−½+iν}` and `J₀` are entire in `ν`.
* The only singularities are at `p = ±m/Λ` (`y = 1`) and on the branch cut `|p| ≥ m/Λ` of the half-integer powers.
* So `E_ρ` is chosen with its right vertex at `1 + f·(m/Λ − 1)`, `f ∈ {0.4, 0.6, 0.8, 0.9}` (`f = 0.3` for `m ≥ 1.3Λ`), and the best `ρ` is kept.

**Bound on `M`.** By the maximum-modulus principle it is taken on the boundary of `E_ρ`, covered by 720 complex balls:
* `|ex_k(ν)| ≤ |pre|·[Σ_j Σ_m ĥ^{(j)}_m(|ν|_max)(|S_{j,m}| + U_{j,m}) + window + far]`. Here `ĥ` is the majorant from
  `|(A/2+i)² + ν²/4| ≤ (A/2+i)² + |ν|²/4`, and `M₂` is re-certified for the larger `|ν|`.
* Odd `k`: `|ex_k| ≤ Σ_l |L_l| |ex_l| E_k^{up}/E_l^{lo}`.
* `E_k` and `(1−y)^{−(k+1)/2}`: bounded through real lower and upper bounds of `|1−y|` and `|2−y|`. This step is essential: Arb's complex
  power and inverse give enclosures containing 0 for such balls.
* `|P_{−½+iν}(cosh β)| ≤ cosh(|Im ν|β) P_{−½}(cosh β)` (Mehler–Dirichlet), and `|J₀(z)| ≤ cosh|Im z|`.
* `φ_k` by Horner in `acb`; `P_J` and `d^J_{44}` by `acb` recurrences.

**Results** (`cache/certq_100_summary.txt`). The certified quadrature error per `N_s`:

| `m/Λ` | error bound / `N_s` |
|---|---|
| 1.01, `J = 0` | ∞ (no ellipse works) |
| 1.02 | 0.43 (`J = 0`, `F/N_s ≈ 1`) |
| 1.03 – 1.1 | ≤ 3·10⁻⁷ |
| 1.1 – 1.2 | ≤ 3·10⁻⁹ |
| ≥ 1.2 | ≤ 10⁻¹⁸ |

The lowest certified exact-integral values, `F/N_s ≥`:
* `Δφ = 3`: +0.2266 (`+−, J = 5, 1.05Λ`);
* `Δφ = 2`: +0.0719 (`++, J = 6, 1.06Λ`; a wide ball, midpoint ≈ +0.3).

### 7.2 Continuity in `m`: not achieved [NE] — the obstruction and what it would take

**The obstruction.** A certificate on `m ∈ [m_a, m_b]` needs either `F` on Δ-intervals or a bound on `|∂_Δ F|`.

**Approaches tried or costed.**
1. **Interval-valued `Δ` in the moment series.** The per-`n` contributions `X_n(Δ)` have `|∂_Δ log X_n| = O(1–5)`, because the `n`-profile
   shifts with `Δ`, while the sum varies slowly. Ball arithmetic loses this cancellation, so the relative width is about `5δ` for a cell
   of half-width `δ`. A useful enclosure needs `δ ≲ 2·10⁻³` in `Δ`, i.e. about 10⁵ cells for `m ∈ [1.01, 3]Λ` and `J ≤ 10`.
   * A mean-value form of the combined log-weights `log(pre·W_n)`, with digamma enclosures, removes the large `Γ`-ratio
     dependency, but not the `O(1)` per-`n` drift.
   * The Mack coefficients would also have to be built as exact polynomials in `Δ` (not done).
2. **Cauchy / convexity bound.** Use `F ≥ min(F_a, F_b) − (h²/8) sup|F''|`, with `|F''| ≤ 2M_r/r²` and a bound `M_r` on `|F|` over complex
   `Δ` in a stadium of radius `r`. This needs moduli bounds of the moments at complex `Δ`.
   * With `r ≈ 1` and `M_r ≈ 10²|F|`, which is plausible since the per-`n` phases drift like `e^{iO(1)Im Δ}`, it requires `h ≲ 0.1` in `Δ`.
   * That is about 10× the present grid density, roughly 1.3·10³ cells for `Δφ = 3`.
   * Each cell costs one complex-`Δ` evaluation of the `N×M` sums, about 1–2 min. That is about 30 h on 2 cores, plus the
     `Δ`-polynomial Mack construction.
3. **Taylor models** (jets in `Δ`, order `K`, with a Cauchy remainder). These trade cell count for a cost of order `K²` per cell. With
   `K = 10`, `ε/r ≈ 0.025` and `r ≈ 5`, cells of width 0.25 in `Δ` would suffice, at an estimated 50 min per cell.

None of these fits the session budget.

**Numerical continuity evidence [NE].** Adjacent-grid variation relative to the certified margins is ≤ 0.85 (§6). So linear
interpolation stays positive, but this is not a proof.

### 7.3 High `J` near the gap (`Δφ = 3`) [certified at the listed points]

**Certified points.** `J ∈ {11, 12, 13, 14, 16}` at `m/Λ ∈ {1.12, 1.15, 1.2, 1.3}`, restricted to physical states:
* `F^{(64)} > 0`, with lower bounds `F/N_s ≥ +0.49` (`J = 16`, `1.2Λ`);
* the exact-integral certificate follows via §7.1.

**Not certified.** `J = 18` and `J = 20` (at `1.2Λ`) give balls `[−0.11, 2.1]` and `[−1.4, 3.4]` `N_s`. The `n`-tail factor `Q̄_j` and the
far-tail factor `(4(M₂+J))^J` grow with `J`, so these need a larger direct `N`. The `J = 20` build already takes 17 min.

### 7.4 `m ∈ [Λ, 1.01Λ]` and `Δ_g ≥ 300`

**`m ∈ [Λ, 1.01Λ]`.** Not done: the Gauss-error bound fails as `m → Λ`, and positivity on the interval would need §7.2 anyway.

**`Δ_g ≥ 300`.** Not attempted in this pass. A cheap tightening of the tail bounds does not exist: the `n`-tail factor at `ν ≈ 300` comes from
`h^{(j)}_m(ν)` growing like `e^{ν²/(2τ)}`.
