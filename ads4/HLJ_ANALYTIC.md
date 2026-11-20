# H_LJ: analytic control of the non-oscillating beyond-Regge part at large `J` near the gap

*20 Nov 2026, H_LJ-analytic agent. Unrefereed, not committed. No existing document or code file was edited. New code in
`ads4/code/spinfunc/`: `lj_saddle.py` (travelling-wave saddles), `lj_scan.py` (scan over the H_LJ region, real and complex `ν`),
`lj_crit.py`, `lj_critmulti.py` (search for other critical points), `lj_phasecheck.py` (validation against the exact scalar), `lj_contour.py`
(constants of the contour lemma), `lj_scalwin.py` (exact scalar in the Bessel window), `lj_gravphase.py` (graviton phase test,
inconclusive, §3.5). Caches `spinfunc/cache/lj_scan_{1000,300}.txt` (`_v1` = first run, superseded: its complex-`ν` lines had a
branch error, §3.3), `cache/lj_contour.txt`. One core, < 0.3 GB. Project modules imported read-only.*

Grades: [P] proved; [P given X]; [D] derived, not independently checked; [D-sketch]; [NE] numerical evidence; [M] model; [I] open.

## 0. Summary (read this first)

**What is new.** H_LJ was "the non-oscillating (smooth) part of the beyond-leading-Regge correction stays above `−1 … −30 c₂`
(uniform) / `−27 … −3580 c₂` (fold window)", supported by a [D-sketch] ("every saddle order oscillates; smooth terms need the fold or
the edge") and untestable numerically for `J ≳ 100–150`. This note replaces the uniform/fold/edge tolerances by an analytic
structure statement plus a sign-free magnitude bound on a small Regge ("Bessel") window, and derives the structure statement.

1. **Exact structure of the CLS-D crystal integrand [P, given the gauge-fixed integrand as implemented in `l1_crystal.build`] (§1).**
   * `ν` enters only through `A(z)^{iν/2}` (and `z`-independent factors). Hence on any fixed contour
     `|integrand(ν + iy)| = |integrand(ν)|·|A(z)|^{−y/2}` exactly: the local `ν`-frequency `Re κ(z) = Re ½log A(z)` is also the local
     decay rate in `Im ν`.
   * `J` enters the scalar integrand **only through `P_J(t(z))`**, `t = Z·W/√(Z²W²)` (the 3-point powers cancel identically against
     `(Z²W²)^{J/2}`, because `Z² = −m₁₅m₂₅m₁₂`, `W² = −m₃₅m₄₅m₃₄`). This explains the observed `J`-independence of `Φ_ν` (EE_DG300_SMALLNU §1).
   * **Celestial fibration.** `B` (the `Δ`-exponent) and `t` do not depend on the light-ray parameter `x̃₀` (the celestial point); only
     `A` and the slow factors do. The soft direction of the crystal Hessian, where the fold at `ν* ≈ 11.3` lives, is this fibre
     direction [NE]. So **the fold is a property of a `J`- and `Δ`-independent fibre integral; it multiplies the spin-`J` waves and
     cannot create a non-oscillating term by itself** [D]. (Explains EE_DG300_SMALLNU's "fold independent of `J` and `Δ`".)
   * Graviton: same `A`, `B`, same fibration; the spin-`J` contraction of the `⟨TTO_J⟩` structures depends on the base through the same
     angle `t(z′)` (Jacobi polynomials of bounded index), times bounded-degree polarisation factors [D].
2. **Contour/endpoint lemma with explicit constants [P given AG; constants D, `lj_contour.py`] (§2).** Split `ν ∈ [0, Λ]` at
   `ν_B = c·m/(2J+1)` (`c = 3–5`: the Bessel window, `ℓθ ≲ c`). If on `[ν_B, Λ]` the beyond-model correction is a sum of two travelling
   waves analytic in the upper/lower half-strips with decay `e^{−β_*|Im ν|}` (class **AG**), then
   `|δF|/N_s ≤ W_B·max_{ν<ν_B}|Γ| + (E_bd + E_top + E_edge)·M`. With `β_* = 0.97β`: `E_top ≤ 3·10⁻¹³`, `E_edge ≤ 2·10⁻⁶` (`≤ 3·10⁻⁸`
   at `J ≥ 48`, 1000) per `c₂` of wave amplitude; the **window tolerance** (sign-free, worse helicity, C2) is

   | `Δ_g` | `J` | 48 | 64 | 100 | 150 | 200 | 300 | 344 |
   |---|---|---|---|---|---|---|---|---|
   | 1000 | `T_B(c=3)` / `(c=5)` (`c₂`) | 2456 / 1314 | 1913 / 1020 | 1402 / 734 | 1184 / 623 | 1114 / 633 | 1057 / 789 | 1095 / 881 |
   | 1000 | SMOOTH_BOUND uniform / `ν<30` | 29 / 8069 | 14 / 3728 | 5 / 1215 | 2 / 510 | 1.5 / 318 | 1.1 / 213 | 1.0 / 199 |

   | `Δ_g` | `J` | 34 | 40 | 50 | 64 | 80 | 100 | 103 |
   |---|---|---|---|---|---|---|---|---|
   | 300 | `T_B(c=3)` / `(c=5)` | 387 / 214 | 343 / 190 | 298 / 166 | 247 / 150 | 205 / 144 | 179 / 145 | 178 / 145 |
   | 300 | SMOOTH_BOUND uniform / `ν<30` | 11 / 273 | 8 / 187 | 5 / 114 | 3 / 70 | 2.4 / 48 | 1.9 / 36 | 1.8 / 35 |

   (`T_B` includes the non-analytic part of the min-norm EE weight, `|X − ψ₂|` on `ν ∈ [ν_B, 25]`, as a real-axis term.) C3H (`ψ₃`)
   tolerances are 7–15× looser. Under AG the uniform (`1–30 c₂`), edge-shaped and fold-above-`ν_B` items are **not needed**.
3. **AG holds to all orders of the travelling-wave saddle expansion, on the whole H_LJ region [D + NE] (§3).** Splitting
   `P_J = W₊ + W₋` (travelling waves `e^{±i(J+½)θ}`), each wave has the exponent `Φ_± = (iν/2)log A + (Δ/2)log B + G₀ ± i(J+½)θ(z)`.
   Computed (`lj_scan.py`; gap states and `m = 1.3·m_gap`; `J = 50–344` at 1000, `34–103` at 300):
   * the saddle `z_±(ν)` exists and is non-degenerate on the whole window `ν ∈ [3, Λ]` (Newton from the `J = 0` crystal; no fold; smallest
     base eigenvalue `69` (1000) / `32` (300), at `ν = Λ`; the only soft direction is the fibre);
   * **`min_ν Re κ̃_± / (±β) = 1.0022–1.0115`** on `ν ≥ 2m/(2J+1)` (`κ̃ = ½log A(z_±) − log 2`, the local frequency = decay rate);
   * along vertical lines `ν_r + iy` (`ν_r = ν_B … 0.9Λ`, `y ≤ 30/β`, branches continued) the wave decays at **`≥ 1.002 β`**;
   * validation against the **exact scalar** (Arb, independent method): the saddle phase predicts all 23–56 zeros of the exact action
     at `J = 40, 64, 100` (1000) and `64, 100` (300) with circular spread **0.004–0.06 rad**, 2–6× better than the flat phase `ℓθ`;
   * a randomised search finds **no other critical point with `|Re κ̃| < β`** in the bulk (`ν = 60, 300`); near the edge (`ν = 900`)
     only critical points with `|κ̃| ≥ 2.5β` (oscillating); the only low-frequency critical points are artefacts of the asymptotic
     wave form at `ℓθ ≈ 1–2`, inside the Bessel window.
   Each order of the Laplace expansion around a non-degenerate saddle analytic in `ν` is analytic in `ν` [standard]: so **to all orders
   the beyond-model correction on `[ν_B, Λ]` is in class AG with `β_* ≥ 0.97β`** [D]. Graviton: same exponent, same saddles [D].
4. **Checks [NE].** Exact scalar in the Bessel window: `max|Γ_s − κ| = 6.5–9.9 c₂` for every `J` (40–100, both `Δ_g`), **160–380× (1000)
   and 25–80× (300) below `T_B(c=3)`**; in the bulk the exact scalar correction is `10³–10⁵ c₂` pointwise yet integrates to `≤ 10⁻³ c₂`-equivalent — the AG
   mechanism in action. Graviton window data: `(1064, 64)` C2 `|Γ − κ| ≤ 350` at `ν = 13–33` (one noisy node `−2653 ± 1619` at
   `ν = 10`) vs `T_B = 1900`; C3H `≤ 1.8·10³ (±1.2·10³)` vs `1.5·10⁴`. `Δ_g = 300`, `J ≥ 34`: the window `ν < 6–14` is the
   thimble region, **no graviton data** there; `J = 24, 28` (300) have `Γ(ν→0) ≈ 430–480` (fold fits), above `T_B(c=5)` and near
   `T_B(c=3)` of their neighbours.
5. **Resulting grade of H_LJ (§5).** `H_LJ` ⇐ **`H_AG`** (bulk: [D] to all orders; beyond all orders only [NE + M] from a critical-point search)
   **+ `H_win`** (`|Γ − κ| ≤ T_B` on `ν < c·m/(2J+1)`) **+ `H_edge-amp`** (oscillating edge amplitude `≲ 10⁶–10⁷ c₂` at 1000, `≲ 10⁵–10⁶` at 300; margin over the exact-scalar envelope
   ×3–20 at 1000, ×1.6 at `J = 34` (300), rising with `J`).
   * **`Δ_g = 1000`:** `H_win` is [NE] at `J = 64` (C2 inside at 1σ; C3H inside), consistent at `J = 100`, [M] beyond (scalar window
     correction `J`-independent; graviton `10²–10³` measured at `J ≤ 100` vs `T_B ≥ 1050` (c = 3) for all `J` up to 344). **H_LJ(1000):
     [I] → [D + NE + M]**: no item now needs a cancellation of the smooth part to `10⁻²–10⁻⁴`; the sign-free window needs `|Γ| ≲ 10³`.
   * **`Δ_g = 300`:** `H_AG` [D] as at 1000. `H_win` needs `|Γ − κ| ≲ 180–390 c₂` (c = 3) on `ν < 6–20`, unmeasured for `J ≥ 34`, and the only
     graviton window data (`J = 24, 28`, `Γ(0) ≈ 450`) are of that size. It holds with a large margin **if** the window correction is
     leading-Regge shaped (H₂ harmonics at the state's `β`, the c-class of SUBLEADING_REGGE §3, tolerance `~10⁵–10⁶ c₂`) — a
     [D-sketch] (§5.2) supported by the `J = 24, 28` fold fits ("dominantly `∝ D_h`"). **H_LJ(300): [I] → [D] (bulk) + [D-sketch / M]
     (window).**
6. **What remains** (§6): a beyond-all-orders statement (thimble decomposition of each wave contains only the crystal saddle); the
   Regge-window structure at subleading order (the spinning `O(r)` Regge kernel, i.e. the `H_R` item, at near-gap `β ≤ 0.53`);
   graviton window data at `Δ_g = 300`, `J ≥ 34`; graviton edge amplitude (`x > 0.8`); the O(`1/Δ`) prefactors of the saddle expansion
   (needed for `c₂`-level comparisons, §3.5).

---

## 1. Exact structure of the crystal integrand

### 1.1 `ν`-dependence [P]
`l1_crystal.build`: the gauge-fixed CLS-D integrand is `exp(λF + G)` (times `z`-independent normalisations `N(ν)`, `1/β_{Δ,J}`), with
`λF = (iν/2)log A + (Δ/2)log B` and `G` independent of `ν` (measure, `(−x²)^{−2d}`, FP determinant, dual, celestial, subtraction,
three-point remainders). Therefore, for every contour `𝒞` (fixed, `ν`-independent):
`I(ν) = ∫_𝒞 e^{(iν/2)log A} w(z) dz` is a Fourier–Laplace transform in `L = ½log A`, and
**`|e^{λF+G}(z; ν + iy)| = |e^{λF+G}(z; ν)|·|A(z)|^{−y/2}`** pointwise. The normalisation `N(ν)` removes the constant
`½log A = log 2` of the `J = 0` crystal (checked: `κ_c = log 2` exactly at `ν = 12, 20, 200`), so the relevant local frequency is
`κ̃ = ½log A − log 2`.

### 1.2 `J`-dependence of the scalar [P]
`Z = x₁₅|x₂₅²| − x₂₅|x₁₅²|` gives `Z² = x₁₅²x₂₅²x₁₂²` identically, i.e. `Z·Z = −m₁₅m₂₅m₁₂` in the code's `msq = −x²` convention; likewise
`W² = −m₃₅m₄₅m₃₄`. The three-point slow term contains `−(J/2)[log m₁₂ + log m₁₅ + log m₂₅ + log m₃₄ + log m₃₅ + log m₄₅]` and
`Q_J = (Z²W²)^{J/2}P_J(t)`. Their product is `(±1)^{J/2}P_J(t)`. **The whole `J`-dependence at fixed `Δ` is `P_J(t(z))/β_{Δ,J}`.**
This is the mechanism behind "`Φ_ν` came out `J`-independent along the symmetric saddle family (mechanism not analysed)"
(EE_DG300_SMALLNU §1).

### 1.3 Celestial fibration [P structure + NE]
* `x̃₀` enters only the celestial point `x₀` (through `y₀² = [uv·uv(1−x̃₀²) − y_p²]/(uv − y_p²)`). `B = m₃₅m₄₅m₁₂/(m₁₅m₂₅m₃₄)` and
  `Z, W` are built from `x₁…x₄, e` only, so **`B` and `t` are exactly independent of `x̃₀`** [P]. Checked at random complex base points:
  `Δlog B = 0`, `Δt = 0` to round-off, while `Δlog A = O(1)` (scratch, `lj_saddle` functions).
* At the crystal base point, `A(x̃₀)/A(0) = 1/(1 − x̃₀²)` exactly on real `x̃₀` (to 12 digits; `κ̃ = log cosh w`, `x̃₀ = tanh w`); not
  at general base points.
* The crystal Hessian (with or without the spin term) has exactly one soft eigenvalue, `≈ 8.8 / 51 / 73 / 200` at `ν = 4 / 10 / 12 / 20`
  (`Δ = 1100`), all others `≥ 10⁴` (`ν ≤ 30`); its eigenvector is the `x̃₀` direction (fibre weight `≥ 0.99` for `ν ≤ 20` at
  `Δ = 1100, 1300`, `≥ 0.96` at `Δ = 364`; at `ν = 100` it mixes, weight `0.70–0.93`, but is no longer soft: eigenvalue `≥ 550`). Adding the spin term
  changes it by `< 0.3 %` at `J = 100` (`72.97` vs `73.13` at `ν = 12`): **`θ` is constant along the fibre** (first order exactly,
  second order `d²t ≈ 10⁻⁶` per unit soft coordinate, i.e. a phase `J·δθ ≲ 10⁻²`).
* Hence the integrand is (base factor carrying all `Δ`- and `J`-dependence) × (fibre factor carrying `x̃₀`), and the fold at
  `ν* ≈ 11.3` (where the soft fibre saddle structure changes) belongs to a `J`- and `Δ`-independent fibre integral
  `g_ν(z′) = ∫dx̃₀ A^{iν/2}e^{G₂}`. At large `ν` it is dominated by its own saddle; at `ν ≲ 11` it is an `O(1)`-width integral (at the
  crystal base point the real-`x̃₀` integrand is `(1 − x̃₀²)^{−3/4 − iν/2}` numerically, so `g` is the Beta function `B(½, ¼ − iν/2)`: poles
  at `Im ν = −(½ + 2n)`, power-law bounded elsewhere [NE]). **It multiplies the
  travelling waves `e^{±iℓθ(z′)}` of the base and cannot generate a non-oscillating term**; its analytic singularities sit at
  `|Im ν| = O(1)`, i.e. in the Regge (Bessel) window [D].
* The fold location depends on which slow terms are put in the exponent (`n1_direct` with the `w`-Jacobian: 11.3; here, `z`-variables:
  the symmetric saddle continues to `ν ≈ 9` and then jumps): a representation artefact, as SMOOTH_BOUND §1.3 suspected [NE].

### 1.4 Graviton [D]
`l6_spin_action.build_spin`: the same `λF` (same `A`, `B`) and the same geometry; `Gs` differs from `G` by external-dimension and
dual-structure powers. The three-point powers of `m₁₂, m₃₄` carry `+J/2` (the `⟨TTO_J⟩` structures contain `V_{5,12}^{J−a}` with
`V ∝ 1/x₁₂²`), so the spin-`J` dependence combines into `|Z|^J|W|^J × (polynomial in t of degree J)/(m…)^{J/2}`, i.e. Jacobi polynomials
`P^{(a,b)}_{J−n}(t)` with `|a|, |b| ≤ 8`, `n ≤ 4`, times polarisation factors of bounded (J-independent) degree. The light-ray
polarisations `z₁ = x₁₀`, `z₃ = x₃₀` and the dual structure depend on `x̃₀`: they are fibre factors of bounded degree. So the graviton
integrand is a finite sum of (fibre factor) × (base factor with the scalar's large exponent `(Δ/2)log B ± i(J+½)θ`) [D].

## 2. Contour/endpoint lemma

**Setting.** `δF = ∫₀¹ ψ₂(p) e_h(p) Γ(p) dp` (+ the `ψ₃` analogue), `e_h = env_h(x)/m²`, `x = Λp/m`, `Γ` the beyond-model correction in
`c₂` units (rows convention, κ-image subtracted). `ψ₂ = (1−p)⁵pX(p)/k₀` is a polynomial; `e_h` is analytic for `|x| < 1`.

**Lemma C [P].** Let `0 < p_B < 1`, `T > 0`, `R_± = {p_B ≤ Re p ≤ 1, 0 ≤ ±Im p ≤ T}`. Suppose that on `[p_B, 1]`, `Γ = Γ₊ + Γ₋` with `Γ_±`
analytic in `R_±`, continuous on its closure, and `|Γ_±(p)| ≤ M e^{−Λβ_*|Im p|}` there (**class AG**). Then by Cauchy
`∫_{p_B}^1 ψeΓ_± = ∫_{p_B}^{p_B ± iT} + ∫_{p_B ± iT}^{1 ± iT} − ∫_1^{1 ± iT}` and
`|δF| ≤ |∫₀^{p_B}ψeΓ| + M[2∫₀^T|ψe|(p_B+iy)e^{−Λβ_*y}dy + 2e^{−Λβ_*T}∫₀¹|ψe|(p+iT)dp + 2∫₀^T|ψe|(1+iy)e^{−Λβ_*y}dy]`.

**Constants [D; `lj_contour.py`].** Normalised by `N_s` (hybrid state), with `ψ₂(p_B + iy) ≈ ψ₂(p_B)(1 + y/p_B)` (`p_B ≤ 0.06`) and
`ψ₂ ≈ X(1)(1−p)⁵` at the edge: `W_B = Σ_{ν<ν_B} p₂ env/(m²N_s)`; `E_bd = 2|ψ₂(p_B)|env(x_B)(1 + 1/(ν_Bβ_*))/(m²Λβ_*N_s)`;
`E_top = 2e^{−β_*y_T}max|ψ₂ env|/(m²N_s)` at `y_T = ΛT = 30/β_*`; `E_edge = 2·5!|X(1)|env(x_max)/(m²φ′⁶N_s)`, `φ′ = Λβ_*/√(1−x_max²)`
(SMOOTH_BOUND's osc-edge, doubled). `β_* = 0.97β` (§3.3 gives `≥ 1.002β`).
The EE weight is the min-norm Abel solution `X(ν)`, not `ψ₂`: `|X − ψ₂|/ψ₂ = 0.3–0.9` at `ν ≤ 10`, 2–4 % at 12–18, and beyond `ν ≈ 20`
the analytic local symbol (`1.0 % / 0.24 % / 0.03 %` at `ν = 20 / 30 / 100`, 1000; the non-local null-direction tails are `≤ 10⁻⁴` of
`max ψ₂`, EE_FUNCTIONAL §3). So the deformed weight is the analytic local form, and `W_X = Σ_{ν∈[ν_B,25]}|X − ψ₂| env/(m²N_s)` is added
as a real-axis term (`u3_chi_EE` cache). `ψ₃` (`Ĉ₃`, no de-shift) is local. Results (`cache/lj_contour.txt`, worse helicity `+−`):

| `Δ_g`, `J` | `W_B` (c=3) | `E_bd` (c=3) | `E_top` | `E_edge` | `T_B` C2 without `W_X` (c=3 / 5 / 8) | `T_B` C3H (c=3) | both keys (c=3) |
|---|---|---|---|---|---|---|---|
| 1000, 48 | 1.3e-4 | 2.3e-4 | 5e-15 | 3.4e-8 | 2456 / 1314 / 691 | 18198 | 2164 |
| 1000, 64 | 1.7e-4 | 3.0e-4 | 1e-14 | 1.1e-8 | 1913 / 1020 / 526 | 14989 | 1696 |
| 1000, 100 | 2.2e-4 | 4.1e-4 | 5e-14 | 2.1e-9 | 1409 / 734 / 373 | 12066 | 1261 |
| 1000, 150 | 2.5e-4 | 4.6e-4 | 1e-13 | 4.2e-10 | 1211 / 624 / 315 | 11402 | 1095 |
| 1000, 200 | 2.4e-4 | 4.4e-4 | 2e-13 | 1.3e-10 | 1232 / 641 / 323 | 12581 | 1122 |
| 1000, 300 | 1.8e-4 | 3.3e-4 | 2e-13 | 2.0e-11 | 1582 / 824 / 410 | 18594 | 1458 |
| 1000, 344 | 1.5e-4 | 2.8e-4 | 3e-13 | 1.1e-11 | 1847 / 969 / 482 | 22945 | 1710 |
| 300, 24 | 6.2e-4 | 1.1e-3 | 2e-15 | 2.0e-6 | 522 / 291 / 160 | 3527 | 454 |
| 300, 34 | 8.2e-4 | 1.5e-3 | 1e-14 | 5.6e-7 | 392 / 214 / 115 | 3052 | 347 |
| 300, 64 | 1.1e-3 | 2.0e-3 | 6e-14 | 5.5e-8 | 286 / 153 / 79 | 2806 | 259 |
| 300, 100 | 1.0e-3 | 1.8e-3 | 1e-13 | 9.6e-9 | 297 / 159 / 81 | 3505 | 274 |

(`W_B`, `E_*` per `c₂`, in `N_s`.) With `W_X` (C2, c = 3 / 5): `W_X = 0` at 1000 `J ≤ 64` (`ν_B ≥ 25`), `3·10⁻⁶ … 3·10⁻⁴` at
`J = 100–344`; `4·10⁻⁶ … 2·10⁻³` at 300; the resulting `T_B` (c = 3 / 5) are those of the §0 table, e.g. 1000 `J = 300`: `1057 / 789`;
300 `J = 100`: `179 / 145` (`cache/lj_contour.txt`; `lj_contour_v1.txt` is the run without `W_X`). Reading:
* **The bulk (`ν ≥ ν_B`) contributes `≤ (E_bd + E_top + E_edge)·M`.** `E_top` is negligible for any conceivable `M`; `E_edge·M` stays
  `≤ 10⁻²` up to edge amplitudes `M ≈ 10⁶ c₂` at 1000 (`J ≥ 64`) and `≈ 10⁴–10⁵` at 300 (exact-scalar edge envelope `(200–700)(1−x²)⁻³`,
  i.e. `≤ 5·10⁵` at `J = 64`, 1000 and `≤ 3·10⁴` at `J = 64`, 300).
* `E_bd·M` is the price of the artificial boundary at `p_B`; it is of the same order as the window term, so the window tolerance is
  effectively a tolerance on `max(|Γ| on the window, envelope at ν_B)`.
* `T_B` is nearly `J`-independent (`ν_B ∝ m/(2J+1)` shrinks while `N_s ∝ b̂⁻³` falls).
* Not used: the cancellation between the two `E_bd` terms (they are complex conjugate for a real `Γ`), the `ψ₃`/`ψ₂` relative signs, any
  sign information. `T_B` is a crude, sign-free bound.

## 3. Travelling-wave saddles at large `J`

### 3.1 The exponent (`lj_saddle.py`)
`P_J = W₊ + W₋`, `W_± ≃ e^{±i((J+½)θ − π/4)}/√(2π(J+½)sin θ)`, `θ = arccos t`. Wave exponent
`Φ_±(z; ν) = (iν/2)log A + (Δ/2)log B + G₀ ± i(J+½)θ(z) − ½log sin θ(z)` (`G₀` = the `J = 0` slow exponent, exact by §1.2). Saddle by
Newton (TPS degree 2), continued in `J` from the `J = 0` crystal at `ν = Λ`, then in `ν`. Envelope theorem:
`dΦ_±(z_±(ν), ν)/dν = (i/2)log A(z_±) = i(κ̃_± + log 2)`. Along complex paths the branches of `(Z²W²)^{−1/2}` and of `θ` are continued
(`QREF`, `SREF`); the first run (`cache/lj_scan_*_v1.txt`) lacked this and showed spurious sign flips of the decay rate (§3.3).

### 3.2 Results on the real axis [NE] (`lj_scan.py`, `cache/lj_scan_{1000,300}.txt`)
| `Δ_g` | states | Newton | `min Re κ̃₊/β` (`ν ≥ 2m/(2J+1)`) | `Re κ̃₋` | min base eigenvalue | soft (fibre) eigenvalue at `ν = 3` |
|---|---|---|---|---|---|---|
| 1000 | gap, `J = 50, 64, 100, 150, 200, 250, 300, 344` | all converge, `ν ∈ [3, Λ]` | 1.0028–1.0036 (at `ν ≈ 57–69`) | `= −Re κ̃₊` | 69–71 (at `ν = Λ`) | 4.5–6.1 |
| 1000 | `m = 1.3 m_gap`, `J = 50, 100, 200, 300` | all | 1.0022–1.0027 | mirror | 1520–1820 | 4.8–5.9 |
| 300 | gap, `J = 34, 50, 64, 100` | all | 1.0096–1.0115 (at `ν ≈ 33–39`) | mirror | 32 (at `ν = Λ`) | 4.5–6.0 |
| 300 | `1.3 m_gap`, `J = 34, 64, 100` | all | 1.0074–1.0087 | mirror | 478–553 | 4.7–6.0 |

`Re κ̃₊` rises above `β` towards the edge (`0.95` vs `β = 0.47` at `ν = Λ`, `J = 300`; it equals the flat frequency `(2J+1)/√(m²−ν²)`
to `≤ 1.5 %` for `ν ≤ 0.7Λ`; near `ν = Λ` it exceeds it, by 5 % at `0.9Λ` and 22 % at `Λ` for `J = 300`) and in the Bessel window. No fold anywhere in `[3, Λ]`: the travelling-wave problem has a single
non-degenerate saddle per wave on the whole window, uniformly in `J` up to the twist cap of the near-gap region.

### 3.3 Complex `ν` [NE]
Ratio to the `J = 0` action at the same `Δ` (cancels `N(ν)` and common power laws), along `ν_r + iy`, `y ≤ 30/β`, branches continued:
decay rate `−d Re(Φ₊ − Φ₀)/dy ≥ 1.002β` at `ν_r = ν_B(c=5), 0.1Λ, 0.5Λ, 0.9Λ` for every state (1000: gap `J = 50, 64, 100, 150, 200,
300, 344` and `1.3m_gap`, `J = 100, 300`; 300: gap `J = 34, 50, 64, 100, 103` and `1.3m_gap`, `J = 34, 100`); `≥ 1.16β` at `0.9Λ`. (The
real-axis entries for `1.3m_gap`, `J = 50, 200` (1000) and `J = 64` (300) are from the first run, whose real-axis part was unaffected.) Spot checks with `ν_r = 6–30` and `y` up to `6ν_r` (`J = 64, 100, 300` at 1000; `34, 100` at
300): `≥ 0.97β`. **So AG holds at the saddle level in the full strip `y ≤ 30/β` over `ν ≥ ν_B`**, which is what Lemma C uses
(`E_top ∝ e^{−30}`). `β_* = 0.97β` is used in §2 with this margin.
*Pitfall (recorded):* with principal branches, `√(Z²W²)` flips sign where `t` crosses 0 and `θ → π − θ`; this produced apparent rate
flips `+β → −β` in the first run. Not physical.

### 3.4 Other critical points (`lj_crit.py`, `lj_critmulti.py`) [NE]
80 random starts per (state, `ν`, wave) (base perturbations `0.05–2·a` complex, `x̃₀ ∈ (−0.9, 0.9)(1 + 0.3i·N)`):
* bulk (`(1100,100)` at `ν = 60`; `(1300,300)` at `ν = 300`, both waves): only the crystal-type saddle (35–46 converged starts, 1 distinct);
* near the edge (`(1100,100)`, `ν = 900`; `(364,64)`, `ν = 200`): 2–3 further critical points, all with `|Re κ̃| ≥ 2.5β` (oscillating),
  some with larger `Re Φ` (+3.5 … +6.1) — off the integration cycle, since the computed actions agree with the crystal saddle;
* at `ν = 10` (`(1100,100)`, `(1300,300)`, `(364,64)`): one extra point with `|Re κ̃| ≈ 0.09–0.14β` and `ΔRe Φ = +0.2 … +0.8`, at
  `ℓθ ≈ 0.8–1.7`. It is created by the `−½log sin θ` term of the asymptotic wave (which diverges at `θ → 0`); the exact wave has only a
  logarithmic singularity there. It lies in the Bessel window (`ν < ν_B`), which Lemma C does not deform. For `(1300,300)`, `ν = 10`
  is `4.6 m/(2J+1)`, inside `ν_B(c=5) = 10.5`.
Grade: [NE + M] (a search, not a classification of the thimble decomposition).

### 3.5 Validation against the exact scalar, and a graviton phase test
* `lj_phasecheck.py` [NE]: zeros of the exact `ex₂(ν)` (sm_exact, Arb) vs the saddle phase `φ = Im(Φ₊ − Φ₀)`:

  | state | zeros | spread of `φ` mod π (saddle) | flat `ℓθ` | leading Regge `νβ` |
  |---|---|---|---|---|
  | `(1040, 40)` | 23 | 0.004 rad | 0.007 | 0.77 |
  | `(1064, 64)` | 37 | 0.008 | 0.020 | 0.81 |
  | `(1100, 100)` | 56 | 0.011 | 0.058 | 0.85 |
  | `(364, 64)` | 31 | 0.034 | 0.073 | 0.76 |
  | `(400, 100)` | 45 | 0.061 | 0.234 | 0.78 |

  The travelling-wave saddle **is** the oscillation of the exact action, including its `J²/m` phase shift.
* `lj_gravphase.py` [inconclusive]: on the graviton `K₂`/direct rows `J = 32, 36, 40` (1000) the scalar-saddle phase differs from the
  κ-phase by only `0.006–0.016` rad (`x = 0.2–0.6`), yet `m²·0.01 rad` is `~10³ c₂`, so the image built from the leading-order saddle
  phase (no `O(1/Δ)` prefactor phase) is worse than κ (rms 1100–1350 vs 106–323 on `x ∈ [0.1, 0.45]`). A `c₂`-level phase needs the
  prefactor corrections; not done, not used.

### 3.6 Derivation of AG to all orders [D]
For each wave and `ν` in `S = {ν_B ≤ Re ν ≤ Λ, |Im ν| ≤ 30/β}`: (i) by §1.3–1.4 the large parameters `Δ, J` sit in the base exponent
(`(Δ/2)log B ± i(J+½)θ`) and the `ν`-dependence in `A^{iν/2}`; the fibre is integrated exactly. (ii) The base saddle `z_±(ν)` is
non-degenerate on `S` (§3.2–3.3), so the Laplace expansion `H_± ≃ e^{Φ_±(z_±(ν),ν)}Σ_k c_k(ν)Δ^{−k}` has coefficients analytic in `ν`
on `S` (polynomials in derivatives of the exponent and prefactors at `z_±`). (iii) `|e^{Φ_±}|` decays like `e^{∓Re κ̃·Im ν}` with
`Re κ̃₊ ≥ 1.002β` on `S`. Hence every partial sum is in class AG with `β_* ≥ 0.97β`, and so is the κ-model (explicit: `D_h(θκ)` is a
Jacobi polynomial in `cos θκ`, `P_{−½+iν}` and `J₀` have the standard travelling-wave splits) — so the beyond-model correction is AG to all
orders. The same holds for the graviton (§1.4: same exponent, prefactors analytic near `z_±`). What is **not** derived: that no other
critical point enters the thimble decomposition of the original cycle with comparable weight and low frequency (beyond all orders);
§3.4 is evidence only.

## 4. Checks against data

### 4.1 Exact scalar (`lj_scalwin.py`) [NE, exact]
`Γ_s = m²(ex₂/E₂ − P_J)`, minus the κ-image (C2 `++` convention):

| state | `max|Γ_s − κ|`, `ν < ν_B(3)` | `ν < ν_B(5)` | bulk (`ν ≥ ν_B(5)`, `x < 0.8`): max / rms | integrated remainder (SMOOTH_BOUND) |
|---|---|---|---|---|
| `(1040, 40)` | 6.5 | 9.6 | 4923 / 1125 | `C_eq = 1.3·10⁻³ c₂` |
| `(1064, 64)` | 6.6 | 9.6 | 17326 / 3539 | `4.5·10⁻⁴` |
| `(1100, 100)` | 6.6 | 9.6 | 58207 / 10875 | `8·10⁻⁴` |
| `(324, 24)` | 6.6 | 9.5 | 1897 / 452 | `0.078` |
| `(340, 40)` | 6.7 | 9.6 | 6745 / 1415 | `1.5·10⁻³` |
| `(364, 64)` | 6.9 | 9.7 | 30322 / 5656 | `1.3·10⁻³` |
| `(400, 100)` | 7.3 | 9.9 | 99172 / 24016 | `3·10⁻⁴` |

* Window: `J`-independent `≈ 6.5–10 c₂`, i.e. `T_B(c=3)/Γ_win ≈ 160–380` (1000), `25–80` (300).
* Bulk: pointwise corrections of `10³–10⁵ c₂` (growing like `J²`) integrate to `≤ 10⁻³ c₂`-equivalent: the AG/Lemma C mechanism, now with a
  derivation of why (§3) instead of an observation.

### 4.2 Graviton (existing data, FOLD_WINDOW_LARGEJ §2, EE_DG300_SMALLNU, SMOOTH_BOUND §1.3) [NE]
| state | window `ν < ν_B` (c = 3 / 5) | data in the window | `T_B` C2 (with `W_X`) / C3H (c = 3) | verdict |
|---|---|---|---|---|
| `(1064, 64)` | 24.7 / 41.1 | `ν = 10–33.4`: C2 `Γ − κ = −111 … −344, +112` (± 80–385), `ν = 10`: `−2653 ± 1619`; C3H `−838 … +1845` (± 180–1240), `ν = 10`: `−2727 ± 1366` | 1913 / 14989 | C2 inside except the noisy `ν = 10` node (1σ consistent); C3H inside |
| `(1100, 100)` | 16.3 / 27.2 | `ν = 13.3–33.4` (node values not tabulated there); smooth-part fits C2 `r = −378 … +586` (± 125–250), C3H `±600–1500` | 1402 / 12066 | consistent [NE, weak] |
| `(1200, 200)` and beyond | 8.8 / 14.7 | none (one-node probe, 6–8 % errors) | 1114 / 12581 | [I] |
| `(336, 24/28)` (300; not H_LJ) | 19.7 / 32.8 | fold fit `Γ ≈ c₀D_h + r`, `c₀ = 427 ± 83 / 483 ± 122`, `r = +3 … +86` | 520 / 3527 (J = 24) | `|Γ(0)| ≈ T_B(3)`: the sign-free bound is marginal; the correction is `D_h`-shaped |
| `(364, 64)` | 8.3 / 13.8 | `ν ≥ 11.7` only (10 nodes): per-node errors `10–370`; window `ν < 8.3` unmeasured | 247 / 2806 | [I] |
| `(400, 100)` | 5.8 / 9.6 | `ν ≥ 10` (noisy `±1100–2400` at 10) | 179 / 3505 | [I] |

## 5. Resulting grade of H_LJ

### 5.1 The new decomposition
`H_LJ(Δ_g)` ⇐ **`H_AG`** ∧ **`H_win`** ∧ **`H_edge-amp`** [P given these, Lemma C], where
* **`H_AG`:** on `ν ∈ [ν_B, Λ]` the beyond-model correction (both helicities, C2 and C3H) is in class AG with `β_* = 0.97β`,
  `y_T = 30/β`, amplitude `M`. Grade: **[D] to all orders** (§3.6; exact structure [P], saddle data [NE]); beyond all orders [NE + M]
  (§3.4); graviton by the shared exponent [D].
* **`H_win`:** `max_{ν < ν_B}|Γ − κ| ≤ T_B` (§0 table; `c = 3`: ≥ 1050 C2 / ≥ 11000 C3H at 1000; ≥ 178 / ≥ 2800 at 300), the
  envelope of `Γ` at `ν_B` likewise, and the same bound on `ν ∈ [ν_B, 25]` for the part weighted by the non-local `X − ψ₂`.
* **`H_edge-amp`:** the oscillating edge amplitude obeys `E_edge·M ≲ 0.1`: `M ≲ 3·10⁶ / 9·10⁶` (1000, `J = 48 / 64`, growing with `J`),
  `M ≲ 1.8·10⁵ / 1.8·10⁶` (300, `J = 34 / 64`). Exact-scalar envelope `(200–700)(1−x_max²)⁻³` at `x_max`: `1.0·10⁶ / 4.7·10⁵` (1000),
  `1.1·10⁵ / 2.6·10⁴` (300): margin `×3` at `J = 48` (1000) and only `×1.6` at `J = 34` (300) with this crude bound (no `k = 2`/`k = 3`
  cancellation, no cancellation between the two waves' edge terms), `≥ 19` at `J ≥ 64`. Graviton `x > 0.8`: unmeasured (the S2 `H_edge` item for the computed states
  is a different, non-asymptotic question).

What disappeared: the uniform smooth tolerance (`1–30 c₂`, the item that needed a `10⁻²–10⁻⁴` cancellation and was untestable at
`J ≳ 100`), the smooth edge-shaped tolerance `A_E`, and the fold window above `ν_B`.

### 5.2 `H_win` and the Regge structure [D-sketch]
In the window `ν ≲ c·m/(2J+1) ≤ 54` the state is in the Regge regime (`ν ≪ m`, fixed `β`). There the leading-Regge action is exact up to
relative `O(1/m)` (H7 §18; `Re` removes the imaginary `O(1/m)` term), and every term of the split-representation derivation of H7 §18.1
— an integral over the celestial point of `(−z₀·x̂)^{−λ}(−z₀·ŷ)^{−λ′}` times polynomials — gives Legendre functions `P_{l+n}(cosh β)` of
the same `cosh β` (the equivariance argument of §18.3). If this persists at subleading order in `1/m` (the spinning `O(r)` Regge kernel,
`H_R`), the window correction is leading-Regge shaped, i.e. in SUBLEADING_REGGE's c-class, whose sensitivity is `(2–9)·10⁻⁷ N_s` per
`c₂` of `Γ(ν→0)` (`6·10⁻³ … 8·10⁻² N_s` per unit `c`, divided by `(2J+1)²`; computed at 1000 only): tolerance `≳ 10⁶ c₂`. Evidence: the exact scalar window
correction falls on the scale of `D_h`'s first lobe (SMOOTH_BOUND §1.3) and is `J`-independent (§4.1); the graviton fold fits at
`J = 24, 28` are dominantly `∝ D_h` (SMOOTH_BOUND §1.3). Grade [D-sketch + NE].

### 5.3 Grades
* **`Δ_g = 1000`:** `H_AG` [D] (+ [NE + M] beyond all orders). `H_win` [NE] at `J = 64` (C2, C3H) and [NE, weak] at `J = 100`; for
  `J ≥ 150` [M]: it needs `|Γ − κ| ≲ 1050–1200 c₂` (C2, `c = 3`, with `W_X`) on `ν < 9–25`, against measured `10²–10³` at
  `J = 64–100` and a `J`-independent `≈ 10 c₂` in the exact scalar; and [D-sketch] makes it harmless by `10³`. `H_edge-amp`: [M] with margin `≥ 3` (`J = 48`), `≥ 19`
  (`J ≥ 64`) against the scalar envelope.
  **Overall `H_LJ(1000)`: [I] → [D + NE + M]** — the remaining assumptions are magnitude statements (`|Γ| ≲ 10³`) in a Regge window
  and at the edge, plus the beyond-all-orders clause; none requires a cancellation of a smooth part.
* **`Δ_g = 300`:** `H_AG` [D] (+ [NE + M]). `H_win`: needs `|Γ − κ| ≲ 180–390 c₂` (C2, `c = 3`) on `ν < 6–20` (and `< 25` for the
  non-local weight), unmeasured for `J ≥ 34`;
  the only graviton window data (`J = 24, 28`) have `|Γ(0)| ≈ 450`, so the sign-free bound is **not** sufficient by itself; it is
  sufficient with the [D-sketch] Regge structure (margin `~10³`). **Overall `H_LJ(300)`: [I] → [D] (bulk) + [D-sketch + M] (window).**
* Odd `J` and the parity-odd helicity structures: same exponent, same lemma [D]; their windows are unmeasured at 300 (`J ≥ 17`) [M].

## 6. What remains (value order)
1. **Regge-window structure at subleading order** (turns `H_win` into [D] at both `Δ_g`): show that the `O(1/m)`, `O(1/m²)` terms of the
   heavy action at fixed `(ν, β)` are combinations of `P_{l+n}(cosh β)` (same `β`) — i.e. the spinning `O(r)` Regge kernel (`H_R`), at
   `β ≤ 0.53`. Alternatively, graviton window data: `ν < 8.3` at `(364, 64)`, `(400, 100)` with the travelling-wave thimble (below).
2. **Beyond all orders:** show that each wave's integration cycle decomposes into the crystal-type thimble only (Picard–Lefschetz on
   the base, fibre exact). §3.4 found no competitor with low frequency outside the Bessel window.
3. **A method without the sign problem** (now concrete): integrate each wave `W_±` separately on its own thimble through `z_±(ν)`
   (`lj_saddle` provides `z_±`, the Hessian and branch continuation). This gives `Γ` at `J²/m ≫ 1` (all `J` of H_LJ) and also tests 2.
   The exact `W_±` (Legendre `Q`-functions; Hilb/Hankel form near `θ → 0`) should replace the asymptotic form in the Bessel window.
4. Graviton edge amplitude at `x > 0.8` (shared with `H_edge`).
5. `O(1/Δ)` prefactor terms of the saddle expansion (phase accurate to `10⁻³` rad), for `c₂`-level checks (§3.5).
6. In the Bessel window itself the two waves merge (`ℓθ ≲ c`); a sign-free representation there is Hilb's `P_J(cos θ) ≈ (θ/sin θ)^{½}J₀(ℓθ)`
   with `J₀(w) = π⁻¹∫₀^π e^{iw cos φ}dφ`: the exponent `(iν/2)log A + (Δ/2)log B + iℓθ(z)cos φ` has a saddle in `(z, φ)` for every `ν`.
   Not tried.

## Hand-over
* **Nothing is running.** All `lj_scan` runs finished (`cache/lj_scan_{1000,300}.txt`; `_v1` files: first run, complex-`ν` lines invalid).
* **Reproduce** (from `ads4/code/spinfunc/`, `OMP_NUM_THREADS=1`):

  | command | output | time |
  |---|---|---|
  | `python3 lj_saddle.py 1100 100 200,150,100,60,40,30,20,15,12,10,8,6,4 1` | saddle track, `κ̃` vs flat / `β` | 1 s |
  | `python3 lj_scan.py 1000 50,64,100,150,200,300,344` (and `… 1.3`; `300 34,50,64,100,103`) | `cache/lj_scan_<Dg>.txt` | ~1 min/state |
  | `python3 lj_crit.py 1100 100 60 60 0.5`; `python3 lj_critmulti.py` | other critical points (§3.4) | 10 s; 1 min |
  | `python3 lj_phasecheck.py 1000 100` (`300 64 n0t`) | exact-scalar zero test | 5 s |
  | `python3 lj_contour.py 1000 48,64,…` / `300 24,…` | `cache/lj_contour.txt` (appends) | 10 s |
  | `python3 lj_scalwin.py` | §4.1 | 5 s |
  | `W2_DMIN=1000 python3 lj_gravphase.py 40:1040` | §3.5 (inconclusive) | 2 min |
* `lj_saddle` globals `SREF`, `QREF`, `QLAST`: set them along any complex path (see `lj_scan.vertical`); leave `None` on the real axis.
