# Certified near-gap positivity of the route-A′ D-safe action (scalar-skeleton model), Δ_g = 100

Status: 18 Nov 2026, H13′/certification agent. **Unrefereed. Not committed.**

**What changed.** This replaces `CERTIFIED_NEARGAP.md` (Theorems C, C′), which certified the withdrawn B′-uplift `Ω^{NF,SC}`. The
functional is now route A′, `Ω^{A′} = ∫ψ̃₂Λ⁴Ĉ₂ + ∫ψ₃½Λ⁶Ĉ₃` (`H0_SPINNING_FUNCTIONAL.md` §R, §10, §13).

**Code.**
* New: `ads4/code/exactskel/cert_actionA.py` and `certA_inv.py`.
* Reused read-only: `cert_moments.py`, `cert_quad.py`, `h13_common.py`, `h13a_common.py`, and the `spinfunc/h1`/`h2` pickle and helpers.
* Caches: `exactskel/cache/certA_100_Df{3,2}.pkl`, `.log`, `certA_inv_100.txt` (git-ignored).

**Grades.**
* **[certified]**: Arb (192 bits); every truncation replaced by a proven bound; inequalities checked in Arb or exact rationals.
* **[D]**: derived.
* **[DV]**: derived and numerically validated (floating point).
* **[NE]**: numerical evidence.
* **[M]**: model.
* **[I]**: assumption.

## 0. Results

**Theorem A (computer-assisted).** Setting:
* `Δ_g = 100`, `Λ = √(98·99)`;
* identical-scalar skeleton `Δφ ∈ {3, 2}`;
* the grid of `CERTIFIED_NEARGAP.md` §3: `m/Λ ∈ {1.01, …, 1.1, 1.12, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2.0}` (and `{2.5, 3.0}` for
  `Δφ = 3`), `J ∈ {0, 2, 4, …, 10}` for `++` and `{4, …, 10}` for `+−`, with `τ ≥ Δ_g`;
* plus the near-gap high-`J` sample `J ∈ {11, 12, 13, 14, 16}` at `m/Λ ∈ {1.12, 1.15, 1.2, 1.3}` (`Δφ = 3`).

For every such state and channel, the **exact `p`-integral** of the route-A′ action of §1 without the `η`-inversion,
`F^{A′}_0 = ∫₀¹[ψ₂Mod₂ + uσψ₃Mod₃]dp`, is **positive**. The one exception is `++, J = 0, m = 1.01Λ`: there only the 64-node sum is
certified, because no Bernstein ellipse fits below the singularity at `p = m/Λ`.

| skeleton | cases | `F₀ > 0` exact integral | `F^{(64)} > 0` | min certified lower bound `F₀/N_s` | min midpoint | max quad. error / `N_s` (`m ≥ 1.02Λ`) |
|---|---|---|---|---|---|---|
| `Δφ = 3` | 230 (205 with `J ≤ 10`, 25 with `J = 11–16`) | **229** | 230 | **+0.2251** (`+−`, `J = 5`, `1.05Λ`), `J ≤ 10`; +0.047 (`+−`, `J = 16`, `1.2Λ`; wide ball, midpoint ≈ 1) | +0.2256 | `7.6·10⁻⁷` |
| `Δφ = 2` | 179 | **178** | 179 | **+0.2240** (`+−`, `J = 5`, `1.05Λ`) | +0.2253 | `4.9·10⁻⁶` |

**By `m/Λ` band** (certified `F₀/N_s ≥`, `J ≤ 10`, `Δφ = 3` / `Δφ = 2`):

| `[1.01,1.05)` | `[1.05,1.1)` | `[1.1,1.2)` | `[1.2,1.5)` | `[1.5,2]` | `[2.5,3]` |
|---|---|---|---|---|---|
| 0.241 / 0.240 | **0.225 / 0.224** | 0.296 / 0.296 | 0.305 / 0.301 | 0.493 / 0.490 | 0.792 / — |

* The certified bounds sit within `5·10⁻⁴ N_s` of the floating-point model value (f6: min `0.2256`). The balls are much tighter than for the
  withdrawn functional: the old `Δφ = 2` bound was 0.07, from wide balls.
* **Exact `η`-inversion [DV, not certified].** `F^{A′} = F₀ + δF`, with `δF = ∫(ψ̃₂ − ψ₂)Mod₂` taken from the exact inversion of
  `spinfunc/h1` and `h2`.
  * On every certified case, `F₀,lo + δF − 0.05|δF| > 0`.
  * Minimum: `+0.2226` (`Δφ = 3`) and `+0.2214` (`Δφ = 2`) `N_s`, both at `+−`, `J = 5`, `1.05Λ`, for `J ≤ 10`. At `J = 16` it is
    `+0.012`, on the wide ball.
  * `max |δF|/N_s = 0.036`.
  * So **route A′ with the exact inversion is certified up to a floating-point inversion shift**: [certified] for `F₀` plus [DV] for `δF`.
* **`Δ_g = 300`, `1000`:** not certified. See §4.
* **Large-`b̂` window, secondary task.** With a signed region-I bound and a sharper oscillatory bound (`h13a_cert_curv.py`,
  `H13A_SIGNED=1`), the certified window of `ROUTEA_LARGE_B.md` §4.2 moves from `b̂₀ = 115.5` to **`91.2`** at `Δ_g = 100`, and from
  `162.6` to **__B300__** at `Δ_g = 300` (`1000`: __B1000__).

## 1. The certified object

**Model.** Route A′ (`H0_SPINNING_FUNCTIONAL.md` §R item 6), as in `spinfunc/f6_routeA2.py` and `ROUTEA_LARGE_B.md`, per state and node.

**Definitions.**
* `u = Λ²/m²`, `ν = Λp`, `x = p√u`, `y = x²`.
* `D̃ = P_J(1−2y)` for `++`; `P^{(0,8)}_{J−4}(1−2y) = d^J_{44}/(1−y)⁴` for `+−`.
* `fac = (1−y)²` for `++`; `(1−y)^{−2}` for `+−`.
* `E_k = (2−y)/(1−y)^{1+k/2}`.
* `Pc = P_{−½+iν}(cosh β) − J₀((2J+1)x)`.

**Integrand.**
* `Mod₂ = (2−y)D̃ + 2Pc + fac·(ex₂ − E₂P_J − 2Pc)`.
* `Mod₃ = D̃ + Pc + (1−y)^{−2}fac·(Σ_l L_l(3) ex_l/E_l − P_J − 2Pc/E₃)`.
* `F₀ = ∫₀¹[ψ₂Mod₂ + uσψ₃Mod₃]dp`.

These equal `f6`'s flat + curvature + remainder terms with the stripped envelope, rewritten per node. The certified midpoint reproduces
the floating-point model to `≤ 3·10⁻⁴ N_s`, and the float value lies inside every ball checked.

**Inputs.**
* `ex_l = m^{2l}C_{l,ν}`, `l = 2, 4, …, 10`: certified enclosures from `cert_moments.py`, unchanged (rigorous `n`- and `m`-tails;
  `CERTIFIED_NEARGAP.md` §2).
* `ex₃` by quartic interpolation in `k` [M].
* `ψ₂`, `ψ₃`: exact rationals (NF02).
* `N_s = |B₂+curv₂| + |B₃+curv₃|`, floating point, used only for reporting.

**What is certified about what.** As in `CERTIFIED_NEARGAP.md`, this is a statement about the scalar-skeleton **model**. It is not about
graviton moments: H4 (the stripped-envelope transfer), H1/H2 for `k = 2, 3`, and H12′ stay [I].

## 2. Method (`cert_actionA.py`)

1. **64-node sum.** `F^{(64)}`: every quantity is an Arb ball (192 bits). Jacobi and Legendre recurrences in Arb; `P_{−½+iν}` by
   `acb.legendre_p`; `J₀` by `arb.bessel_j`; `ψ_k` by Horner on exact rationals.
2. **Gauss error.** Trefethen ATAP Thm 19.3 on `[0, 1]`: `|∫f − Σw_if(p_i)| ≤ ½(64/15)Mρ^{−128}/(ρ²−1)`.
   * `M` bounds `|f|` on the boundary of a Bernstein ellipse whose right vertex is `1 + frac·(m/Λ − 1)`.
   * The boundary is covered by 720 complex balls. On each ball:
     * `|ex_l| ≤` the `cert_quad.ex_abs_bound` majorant;
     * real-modulus bounds for the powers of `|1−y|` and `|2−y|`;
     * `|P_{−½+iν}| ≤ cosh(|Im ν|β)P_{−½}(cosh β)` and `|J₀(z)| ≤ cosh|Im z|`.
   * Analyticity: the only singularities are at `y = 1` (and the branch cut of `(1−y)^{½}` beyond it), plus the pole of `1/E_l` at
     `y = 2`, which lies farther out.
3. **Positivity.** `F₀ ≥ F^{(64)}_lo − err > 0` is checked for each case.
4. **Run.** 2 cores, about 2.5 h; 20–180 s per state.

**`m = 1.02Λ`.** The route-A′ integrand behaves better there than the withdrawn one: the quadrature error is `6·10⁻⁸ N_s`, against
`0.43 N_s` before. Only `1.01Λ`, `J = 0` remains without an ellipse.

## 3. The exact `η`-inversion (`certA_inv.py`) [DV]

**Definition.** `ψ̃₂ = 𝔍[ψ₂]` is the exact Abel/ODE inversion (`H0_SPINNING_FUNCTIONAL.md` §R item 5, §13). Its numerics (checks V1–V3)
are in `spinfunc/h1_etainv_exact.py`, read through `h1_100.pkl`.

**Computation of `δF`.** As in `h2_neargap_exact.py`:
* the flat + curvature part on `h1`'s fine `ν`-grid;
* the remainder part on the 64 nodes, with the double-precision `ρ₂`;
* tails `p > 1` excluded [I].

**Why it is not certified.**
* `ψ̃₂` is analytic only in the strip `|Im ν| < ½`. A Bernstein-ellipse quadrature bound is therefore useless for `δF` (`ρ ≈ 1 + 1/Λ`).
* A sup-norm bound of `δψ₂` times `∫|Mod₂|` is about 100× too crude, because `δψ₂` is concentrated at `ν = O(1)`.
* Certifying it would need a rigorous enclosure of `ψ̃₂(ν)` on a fine grid: nested Arb integrals of the ODE solution, plus a
  derivative bound. That is not done.

**Result.** With a 5 % allowance for the `h1`/`h2` numerics, `F₀,lo + δF − 0.05|δF| > 0` on every certified case.

| | min, `J ≤ 10` | max `|δF|/N_s` |
|---|---|---|
| `Δφ = 3` | `+0.2226 N_s` (`+−`, `J = 5`, `1.05Λ`) | 0.036 |
| `Δφ = 2` | `+0.2214 N_s` (same state) | 0.029 |

`cache/certA_inv_100.txt`.

## 4. `Δ_g = 300`, `1000`

**Not certified.** A spot check at `Δ_g = 300` (`m = 1.3Λ`, `J = 4`; `cache/certA_300_Df3.log`) __S300__.

**Why.** As in `CERTIFIED_NEARGAP.md` §6 and §7.4, the certified `n`- and `m`-tail bounds of `cert_moments.py` grow like
`e^{ν²/(2τ)}` at `ν ≈ Λ`. Tightening them, or switching to the Arb tail quadrature of `h13_laplace.py` with a proven Euler–Maclaurin
remainder, is the step needed; it is not done. The floating-point model at 300 and 1000 (`f6`/`h2`: min `0.2052`, `0.2021`) stays [NE + M].

## 5. Not certified

* **Continuity in `m` between grid points** [NE]. As in `CERTIFIED_NEARGAP.md` §7.2, this needs a `Δ`-interval method.
* **The cases left out:** `++, J = 0` at `1.01Λ` (exact integral) and `m ∈ [Λ, 1.01Λ)`.
* **The `η`-inversion shift** [DV] (§3).
* **`Δ_g ≥ 300`.**
* **The model inputs:** the odd-`k` definition [M] and the graviton-level hypotheses [I].

## 6. Files

| file | role |
|---|---|
| `cert_actionA.py` | route-A′ action: 64-node Arb sum and certified Gauss error → `cache/certA_<Dg>_Df<Df>.pkl`, `.log` |
| `certA_inv.py` | exact-inversion shift `δF` (from `spinfunc/h1` and `h2`, read-only) → `cache/certA_inv_100.{txt,pkl}` |
| `certA_summary.py` | tables → `cache/certA_100_summary.txt` |
| `cert_moments.py`, `cert_quad.py` | certified moments and majorants (unchanged, reused) |
| `h13a_cert_curv.py` (`H13A_SIGNED=1`) | sharper large-`b̂` window → `cache/h13a_cert_curv_signed.txt` |

