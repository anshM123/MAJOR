# Exact spinning Regge moments at x = ν/m ~ 0.1–1: feasibility study (no prototype)

Status: 18 Nov 2026, numerics agent. **Unrefereed. Not committed.** This is a feasibility study. **Outcome: no exact spinning moment was
computed.** None of the candidate routes reaches a validated prototype (≤ 10⁻³ against the existing exact position-space values) within
the ~5 h budget. Below are the reasons, the one route that is viable in principle, what its exact ingredients are, a staged plan with
validation targets, and the cost.

Marks: **[D]** derived; **[DV]** derived and checked symbolically or numerically; **[NE]** numerical evidence; **[I]** estimate or
assumption.

## 0. Summary

* **Why position space stops at `x ≈ 0.12–0.17` [NE].**
  * The real-section Regge moment of a heavy block is an exponentially small (`∼ e^{−πν}`) result of `O(1)` integrands.
  * `directgrav/c1_calib.py` measured float64 breakdown on the exact scalar skeleton at `ν ≈ 12` (`k = 2`) and `≈ 17` (`k ≥ 3`) for
    `m ≈ 100`.
  * Reaching `x = 1` (`ν ≈ 100`) in position space therefore needs **≈ 145 significant digits** throughout the spinning wave, plus an
    `O(ν)`-times denser `ξ`-grid.
  * For the Taylor-mode ⟨TTTT⟩ waves (25 seed blocks, 196 word pairs, derivative order ≤ 8, ~3000 grid points, 2.3 h per state in
    float64) that is out of reach: an estimated `10³–10⁴` core-hours per state [I].
* **The viable route is the Mellin/Barnes route with the CPPR differential basis** (§2). The spinning wave is already, in the project's own
  code, a finite sum of differential operators acting on scalar seed blocks with shifted external dimensions (`fourpt/spinblock2`, `tttt`).
  Three facts make the moment a finite sum of exact scalar-type Mellin moments:
  * **[DV] The frame of the moments is Mellin-friendly** (symbolic check, §2.2):
    * every `x_ij²` is a monomial in `r = √(ρρ̄)`, `(1−ρ)(1−ρ̄)` and `(1+ρ)(1+ρ̄)`;
    * so `u, v` and every power of `s_ij` in the structures are such monomials, which in Mellin space are shifts of `(s, t)` and of `k`;
    * the polarisation products are polynomials in `ρ, ρ̄, cos θ, sin θ` over powers of `r`.
  * **The derivatives are Euler derivatives** in `P = (log u + log v)/2` and `Q = (log v − log u)/2` (`spinblock2`, jets `G_ab`). In Mellin
    space these are exact polynomial multipliers.
  * **After the circle average, the extra `ρ^aρ̄^b` factors** become `r`-powers (shifts of `k`) times `cosh/sinh(nξ)`. Against `P_{−½+iν}`
    these give finite combinations at shifted `ν → ν ± i·m` (Legendre recursions).
* **Cost [I].** The route needs four new exact ingredients (§2.3). Two are derivations not available in the literature in usable form.
  The estimate is **1–2 weeks**, with a validation target at each stage. Within this budget it was not started beyond the structural checks.
* **Not viable or not exact:**
  * complex-rapidity contour shifts in position space (§1.2): the endpoint at `η = 1` restores the cancellation; the fix is the dispersive
    form, which is the Mellin route again;
  * eikonal / impact-parameter factorisation (§1.3): a `1/m` expansion, which exactskel showed fails at `x > 0.6`;
  * weight-shifting crossing (6j) to scalar exchanges (§1.4): does not remove the derivatives.

## 1. Routes assessed

### 1.1 Extended-precision position space — not tractable [NE + I]

* **Origin of the loss.** The moment `Π̂_{k,ν} = ∫[dη]P_{−½+iν}(η)∫dr r^{k−2}G` is a Mehler–Fock transform of a function analytic in a strip
  `|Im ξ| < π`. The result is `∼ e^{−πν}` relative to the integrand.
  * float64 calibration on the exact skeleton (`directgrav/cache/c1_100_1.040_4_Df3.txt`): `|R_pos − R_ex| ≲ 10⁻³` up to `ν ≈ 10.6`
    (`k = 2`; `2·10⁻³` at `13.8`) and `15.5` (`k ≥ 3`), then `O(1)`–`10¹⁰` errors one or two nodes later. This is consistent with `e^{−πν} ≈ 10⁻¹⁶`.
  * Digits needed: `≈ πν/ln10 + 6`, i.e. **20 / 47 / 88 / 143** digits at `x = 0.1 / 0.3 / 0.6 / 1` (`m ≈ 100`).
* **Second loss for `+−` at low `J`.** It needs a further `10²–10⁵` (`directgrav` §7.3: `|M|/Σσ`).
* **Grid.** The `u`-grid (72 nodes) must resolve `e^{iνξ}`, so about `ν/2` times more nodes at `ν ~ 100`.
* **Cost.**
  * The Taylor stage is about 6000 s per state in vectorised float64. Arb at 500 bits is `30–100×` slower per operation, on a `5–10×` larger
    grid.
  * The seed radial sums (`cblock12.c`, `N = 8000`, double only) would need a multiprecision rewrite.
  * Estimate: **`10³–10⁴` core-hours per state** [I].

### 1.2 Complex rapidity (contour rotation in `ξ`) — not a shortcut [D-sketch]

* **What it would do.** Shifting `ξ → ξ + iθ` moves the loss to `e^{−(π−θ)ν}`. It needs the seeds at complex cross-ratios. At fixed `r`,
  complex `ξ` makes `ρ, ρ̄` complex, inside the disk of convergence.
* **Why it fails.** The `η`-integral starts at the regular point `η = 1`, and the integrand is odd in `ξ`: `(ρ − ρ̄)` with `[dη] = sinh ξ dξ`.
  So the contour cannot be closed into the strip.
* **Splitting does not help.** Write `P_{−½+iν} ∝ Q_{−½+iν} − Q_{−½−iν}` and rotate each piece. Both pieces are log-singular at `η = 1`, and
  their endpoint contributions are `O(1)` and cancel to `e^{−πν}`.
* **What does work** is the dispersive form: `Π̂` as an integral of the discontinuity of `h(η)` on `η < −1` (the other operator ordering),
  with an explicitly `e^{−πν}` kernel. But this is the Froissart–Gribov / Mellin object again, with the extra difficulty that the
  `r`-domain (`r < e^{−ξ}`) couples to `η`. **No gain over §2.**

### 1.3 Eikonal / impact-parameter factorisation — not exact

* The spin dependence factorises into polarisation tensors times scalar `Π_{Δ−1}(b)` derivatives only at leading order in `1/m`, with
  `O(1/m²)` corrections that are not resummed.
* `EXACT_SKELETON_NEARGAP.md` §3 shows that the scalar `1/m` series is already useless at `x > 0.6`, which is the region needed.

### 1.4 Weight-shifting crossing (6j symbols, KKSD 1706.07813) — no simplification

* Crossing a weight-shifting operator through the block trades external-point derivatives for blocks with shifted internal `(Δ, J)`. It
  still acts with spinning differential operators.
* For the Regge moment, the derivatives are anyway exact multipliers in Mellin space (§2). So 6j crossing adds bookkeeping without removing
  any hard step.
* It could help to organise the 14×14 word pairs into fewer independent Mellin structures [I].

## 2. The Mellin/Barnes route with the CPPR differential basis

### 2.1 Structure

* **Spinning wave as differential operators on seeds.** For `l = 2` (`tttt/tmoments_l2.TaylorWaves`):
  * there are 14 words per side, giving 196 word pairs;
  * the external-dimension shifts are `(s₁, s₂) ∈ {(0,4), (1,3), (2,2), (3,1), (4,0)}` per side, giving **25 seed blocks**
    `g^{a,b}_{Δ,J}` with external dimensions `(3+s₁, 3+s₂, 3+s₃, 3+s₄)`, i.e. `Δ₁₂, Δ₃₄ ∈ {0, ±2, ±4}`;
  * the total derivative order is ≤ 8 (`word_budget`), in Euler variables `(P, Q)`.
* **Wave as a finite sum:** `W[X_hL, X_hR] = Σ_{wL,wR} c_{wL}c_{wR} Σ_{a+b≤8} Π^{(wL,wR)}_{ab}(s_ij, z·x) ∂_P^a∂_Q^b g_{seed(wL,wR)}(u, v)`.
  * For `l = 1` this is available symbolically (`fourpt/cache/gen4_l1.py`, jets to order 4).
  * For `l = 2` it is only implicit in the Taylor-mode engine. The symbolic expressions were too large, which is why Taylor mode was
    introduced.
* **Moment under a circle kernel.** Under `K` (`K*` or `𝔄`; both use only `(q, −q)` circle modes, so the harmonic function is `F^ν_0` for
  every term):
  * `Π̂_K = Σ_{seed} Σ_{a,b} ∫ measure · F^ν_0(ξ) · Π^K_{seed,ab}(ρ, ρ̄) · ∂_P^a∂_Q^b g_seed(u, v)`;
  * `Π^K` is the kernel-contracted, circle-averaged structure coefficient, a function of the frame only.

### 2.2 Frame structure [DV, sympy, frame of `e2e/moments.py`]

With `ρ = s₁²`, `ρ̄ = s₂²`, `r = s₁s₂`:

| invariant | value |
|---|---|
| `x₁₃²` | `−4r` |
| `x₂₄²` | `−4/r` |
| `x₁₂² = x₃₄²` | `−(1−ρ)(1−ρ̄)/r` |
| `x₁₄² = x₂₃²` | `−(1+ρ)(1+ρ̄)/r` |
| `u` | `(1−ρ)²(1−ρ̄)²/(16ρρ̄)` |
| `v` | `(1+ρ)²(1+ρ̄)²/(16ρρ̄)` |

* Example polarisation products:
  * `z₁·x₁₃ = 2√r`;
  * `z₁·x₁₂ = [2ρρ̄ + ρ(cos θ − 1) − ρ̄(cos θ + 1)]/(2r^{3/2})`;
  * `z₂·x₂₄ = −2/√r`;
  * `z₂·x₂₁ = [ρ(1 + cos θ) + ρ̄(1 − cos θ) − 2]/(2√r)`.
* **Consequences:**
  1. **Monomials in `s_ij`** are monomials in `u, v` and `r`. They are Mellin shifts of `(s, t)` and shifts of `k`.
  2. **Polarisation products** are polynomials in `ρ, ρ̄, cos θ, sin θ` over `r`-powers. After the circle average `Π^K` is a polynomial in
     `ρ, ρ̄` (times the monomials of item 1).
  3. **Each `ρ^aρ̄^b`** equals `r^{a+b}e^{(a−b)ξ}`.
     * The symmetric part `cosh((a−b)ξ)` is a polynomial in `η`. Against `P_{−½+iν}` it gives finite combinations of `P_{−½+i(ν∓ij)}`
       (recursion `(2μ+1)ηP_μ = (μ+1)P_{μ+1} + μP_{μ−1}`).
     * The odd part `sinh((a−b)ξ)` gives `P¹`-type kernels. Whether it survives the `M = 0` circle average is to be checked.
* So, given the ingredients of §2.3, **the spinning moment is an exact finite sum of scalar-type Mellin moments** with shifted `(Δ₁₂, Δ₃₄)`,
  `k`, `ν` (complex shifts) and polynomial Mellin multipliers.
* **What would be exact:**
  * all Mellin integrals (Barnes-lemma form, Arb);
  * the Mack polynomials (exact rationals);
  * the structure coefficients (exact rational polynomials).
* **Only truncations:** the `n`-sum tail and the polynomial degree, both controlled as in `skelexact3.py`. No `1/m` expansion and no
  position-space cancellation.

### 2.3 Missing ingredients, validation targets, and cost [I]

| # | ingredient | status in the repo | validation target | estimate |
|---|---|---|---|---|
| M1 | **Frame-Mellin kernel for unequal external dimensions with extra `r`-powers and `η`-polynomials**: the CHMRS (2106.10274) passage from the position-space moment to `∫dT Γ(T/2)²γ_{a−T}(ν)B̂_{k,T}` (or the Barnes form of `skelexact3`), generalised to `⟨φ₁φ₂φ₃φ₄⟩` stripping, shifted `k` and complex `ν` | exactskel covers identical `Δφ`, real `ν` only | `e2e/moments.ScalarMoments(Δ, J, D1 ≠ D2)` at `ν ≤ 10`, `≤ 10⁻⁸`; identical-`Δφ` limit = exactskel | 1–2 days (derivation) |
| M2 | **Mack polynomials with `Δ₁₂, Δ₃₄ ≠ 0`** (Costa–Gonçalves–Penedones 1209.4355, general form) and their pole structure in the `B̂` sum | `mack.py` has `Δ₁₂ = Δ₃₄ = 0` only | radial-series blocks `e2e/blocks.CBlock(a, b ≠ 0)` via the Mellin integral at sample points | 0.5–1 day |
| M3 | **Euler derivatives `∂_P^a∂_Q^b` and the `u^αv^β r^c` monomials as Mellin multipliers and shifts** | straightforward | inside M5 | 0.5 day |
| M4 | **Kernel-contracted structure coefficients `Π^K_{seed,ab}(ρ, ρ̄)` for `l = 2`**: run the Taylor engine on unit jets (`G_ab = δ`) at sample `(ρ, ρ̄)`, fit exact rational polynomials (bounded degree, by item 2), check odd `ξ` parts | not available; `l = 1` symbolic (`gen4_l1.py`) | reproduce `tttt` waves pointwise to `10⁻¹²` | 1–2 days |
| M5 | **Assembly; validation of the full spinning moment** | — | `⟨JJJJ⟩` (`l = 1`) and ⟨TTTT⟩ caches (`β = 1`: `m = 39, 70, 101`; directgrav A, C) at `x ≤ 0.05`, `++` to `≤ 10⁻³`, using the cleaned (t19) values as reference; `+−` where resolved (`𝔄`, `k = 5, 7`) | 1–2 days |
| M6 | **Production**: `x ∈ {0.1, 0.3, 0.6, 1.0}`, `k = 2…7`, `K*`/`𝔄`, A and C, both helicities; H4 (skeleton × `c*`) and H1 | — | — | hours per state (Arb; cf. skelexact3, 4–10 s per state per scalar moment set) |

**Total: about 1–2 weeks of focused work.**

**Risks:**
* M1 is the critical step: the domain `r < e^{−ξ}` and the stripping with unequal dimensions must reproduce a Barnes-type closed form.
* In M4 the polynomial degrees may be large (derivative order 8, two vertices), which makes the number of Mellin terms large but finite.
* The `+−` low-`J` cancellation (`10²–10⁵`) is harmless in Arb, but it requires all 196 word pairs exactly.

**Cheaper first milestone.** Do `l = 1` with one spinning pair (`⟨JφJφ⟩`, `SPINNING_E2E.md` caches) first. It needs M1–M3, M4 is already
symbolic (`gen4_l1.py`), and it has the exact references needed for validation.

## 3. What was done in this session

* The routes of §1 were assessed. The precision and cost statements of §1.1 use the `directgrav` calibration (`c1`) and the
  signal-to-noise data (`d3`).
* The frame invariants and polarisation products were factorised symbolically (sympy; the script is inline in this session and the
  results are quoted in §2.2).
* The seed, shift and derivative inventory of the `l = 2` Taylor-mode engine was taken from `structures2.words`, `embed2.shifts` and
  `taylor4.word_budget` (14 words per side, 25 seeds, order ≤ 8).
* The existing exact scalar machinery was checked for unequal external dimensions: `mack.py`, `exactmom.py` and `skelexact3.py` assume
  `Δ₁₂ = Δ₃₄ = 0` and real `ν`.
* **No exact spinning moment was computed. H4 and H1 at `x ≳ 0.1` remain where `DIRECT_GRAVITON_NEARGAP.md` left them ([I]).**

## 4. Recommendation

* Start M1–M3 with the `l = 1` milestone, validated against `e2e` unequal-dimension scalar moments and the `⟨JφJφ⟩` caches.
* Then do M4 for `l = 2` by unit-jet sampling.
* Until then, the near-gap graviton moments at `x ≳ 0.1` have no exact source in this programme. Position-space extended precision is not a
  realistic fallback (§1.1).
