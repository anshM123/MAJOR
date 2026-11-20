# The non-oscillating (smooth) beyond-Regge part of `Ω^EE` at large `J` near the gap: where it can come from, exact-scalar test, the `J = 48` offset, the curvature double count

*20 Nov 2026, smooth-bound agent. Unrefereed, not committed. No existing document or code file was edited. New code in
`ads4/code/spinfunc/`: `sm_exact.py`, `sm_exact_an.py`, `sm_nsum.py`, `sm_nancheck.py`, `sm_bound.py`, `sm_offset.py`,
`sm_curvcheck.py`, `sm_curvfix.py`, `sm_chainfix.py`, `sm_alias.py`, `sm_foldfit.py`, and the scripts `sm_queue.sh`, `sm_queue2.sh`,
`sm_queue3.sh`, `sm_chainfix_all.sh`.
Caches are in `spinfunc/cache/sm_*`. Project modules are imported read-only, and every patch is applied inside the process only.*

Grades: **[P]** proved, **[D]** derived, **[D-sketch]**, **[NE]** numerical evidence, **[M]** model, **[I]** open.
"Tolerance" means SUBLEADING_REGGE §0.1: the uniform non-oscillating `Γ`, in `c₂` units, that would take `0.9 λ` off the gap state.
It is reproduced exactly by `sm_bound.py` (`tol_all` column).

## 0. Summary

1. **The exact scalar skeleton at large `J` near the gap is now computed** [NE, exact; §2].
   * The `exactmom` NaN is a double-precision overflow. The Chebyshev Mack coefficients `c_i R_nⁱ` reach `10²⁶²` and then `inf` on the
     Euler–Maclaurin tail nodes (`n ~ 10⁹`). Even without the overflow, a double-precision Mellin n-sum cannot hold the large-`J`
     cancellation (H13).
   * The fix is to use the Arb Laplace/Barnes route of `exactskel/h13_exact.py` (`Laplace3H`), which is already validated to
     `J = 150`. It takes 1.5–6 min per state.
   * At `Δ_g = 300` near the gap (`J ≥ 40`), h13's default n-sum configurations (`N₀ = 200–320`) diverge at `ν/Λ ≳ 0.85`.
     `N₀ = 20–60` agree to `10⁻¹⁰`, independent of `nq` and `L`, and give n-sum spreads `≤ 6·10⁻⁹ N_s` (`sm_nsum.py`).
   * Checks:
     * the default and small-`N₀` configurations agree exactly at 1000 (`J = 64, 100`);
     * twice the p-nodes changes nothing at `(1040, 40)`.
   * Results, with `C_eq` the uniform-`c₂` equivalent of the whole remainder:

     | `Δ_g` | `J` (gap state) | `|rem/N_s|` `++` / `+−` | `C_eq` (`c₂`), worse helicity | tolerance | margin |
     |---|---|---|---|---|---|
     | 1000 | 40 | `3.8·10⁻⁶` / `2.5·10⁻⁵` | `1.3·10⁻³` | 46 | `> 3·10⁴` |
     | 1000 | 48 | `4.1·10⁻⁶` / `3.9·10⁻⁴` | `0.012` | 29 | `2·10³` |
     | 1000 | 64 | `3.5·10⁻⁶` / `2.9·10⁻⁵` | `4.5·10⁻⁴` | 14 | `3·10⁴` |
     | 1000 | 100 | `4.4·10⁻⁶` / `1.5·10⁻⁴` | `8·10⁻⁴` | 5 | `6·10³` |
     | 300 | 24 | `1.1·10⁻⁶` / `2.8·10⁻³` | `0.078` | 25 | 320 |
     | 300 | 40 | `3.6·10⁻⁵` / `1.7·10⁻⁴` | `1.5·10⁻³` | 8 | `5·10³` |
     | 300 | 64 | `2.6·10⁻⁵` / `3.6·10⁻⁴` | `1.3·10⁻³` | 3 | `2·10³` |
     | 300 | 100 | `4.0·10⁻⁵` / `1.5·10⁻⁴` | `3·10⁻⁴` | 2 | `6·10³` |

   * The scalar remainder is **not small pointwise**. `m²ρ₂` reaches `10⁵` near the light-cone edge. It is purely oscillating,
     and its `ψ`-integral is `10⁴–10⁶` times smaller than its band-wise pieces.
   * In the fold region, `m²ρ₂(ν → 0) = 6.5–6.9` `c₂` at every `J` (24–100) and both `Δ_g`. It then decreases smoothly on the scale
     of `D_h`'s first lobe (`(2J+1)ν/m ≈ 2.4`). There is no Airy feature at `ν* ≈ 11.3`.
2. **Analytic: where smooth terms can come from, and their weights** [D-sketch + D for the weights; §1].
   * Every order of the crystal-saddle expansion is linear in the spin-`J` prefactor and so oscillates. A smooth term needs a
     different critical structure: (a) the fold at `ν* ≈ 11.3` or (b) the light-cone edge.
   * (a) The fold is `J`- and `Δ`-independent. Its soft block has no large parameter, so there is no shrinking Airy layer. It is an
     `O(1)` function of `ν` confined to `ν ≲ ν_F ≈ 25–30`.
     * Its weight is `W_F ≈ X₀ν_F²/(Λ²m²N_s)`, with `X₀ = ψ₂′(0) = 188.9` [D; numerics agree to 5 % at 1000 and 15 % at 300].
     * That equals `(X₀/1024) ν_F² b̂³/(Λ²m²) ∝ ν_F²(2J+1)³Λ/m⁵`.
     * Tolerance on a fold-localised smooth term (`ν < 30`, worse helicity):
       * 1000: `12800 / 7750 / 3580 / 1170 / 490 / 306 / 205` `c₂` (`J = 40 / 48 / 64 / 100 / 150 / 200 / 300`);
       * 300: `653 / 440 / 186 / 70 / 36 / 27` (`J = 24 / 28 / 40 / 64 / 100 / 150`).
     * These are 30–280 times looser than the uniform table.
   * (b) The edge sits at `x_max = Λ/m`, where `ψ₂ ≃ X(1)(1−p)⁵` with `X(1) = −8648`.
     * An oscillating edge term is suppressed by `5!|X(1)|/φ′⁶`, `φ′ = b̂/√(1−x_max²)`: `≤ 3·10⁻⁸ N_s` per `c₂` of amplitude at
       1000, and `≤ 8·10⁻⁷` at 300.
     * A smooth edge-shaped term `A_E((1−x_max²)/(1−x²))³` has tolerance:
       * 1000: `A_E < 3·10⁴ / 1.2·10⁴ / 2900 / 390 / 80 / 32 / 12` `c₂` (`J = 40 → 300`);
       * 300: `3200 / 1600 / 360 / 63 / 18 / 8` (`J = 24 → 150`).
   * **Bound on the integrated smooth part** (§1.5):
     `|dF_sm/N_s| ≤ W_F(J,m)·max_{ν<ν_F}|Γ_sm| + W_E(J,m)·|A_E| + W_bulk·|r_bulk|`. Here `r_bulk` is a smooth part outside the fold
     and edge windows, which the saddle argument excludes [D-sketch].
     * Measured on the graviton rows:
       * `|r_bulk| ≤ 6–12` `c₂` (`J = 32–40` at 1000, `x ∈ [0.1, 0.45]`, K₂ and direct);
       * `|r_bulk| ≤ 20–36` (`J = 24–32` at 300);
       * fold-window smooth part `r = +3…+86` `c₂` (±22–37, depending on the fit class), from direct all-orders nodes at
         `ν = 10–33`, `Δ = 336`, `J = 24 / 28` (`sm_foldfit.py`). It is positive or zero, so `r ≥ −52` `c₂` at 2σ.
     * This is inside the fold tolerances (only the negative side matters) for every `J ≤ 300` at 1000 and `J ≤ 64` at 300. It is **not** inside the uniform table
       for `J ≥ 64` at 1000, and the graviton `r_bulk` is measured only up to `J = 40`.
     * Grade: [D-sketch + NE]. It is not a theorem: the graviton magnitudes are measured at `J ≤ 40`, and `ν < 10` and `x > 0.8` are
       unmeasured for the graviton.
3. **`J = 48` smooth offset: an artefact of the `J = 48` rows, not a fold/edge term and not the curvature double count** [NE; §3].
   * On the same nodes, after the κ-curvature image is subtracted, the fitted constant is:
     * exact scalar `(1048, 48)`: `r = +0.7 ± 3`;
     * graviton `J = 32 / 36 / 40` (1000) and `J = 24 / 28 / 32` (300, with `J²/m` up to 3.1): `r = 0 ± (4–36)`;
     * graviton `J = 48`: `r = +280…+650` (C2 `++`, `+−`, C3H `++`) and `−630…−690` (C3H `+−`), the same in K₂ and direct rows.
   * Residuals are 30–80 % at `J = 48`, against 1–10 % elsewhere.
   * It is not a `J²/m` effect, not the de-shift (C3H shows it too), not aliasing (the scalar control on the same nodes is clean), and
     not a `1/D_ref` reference-zero factor (adding that basis leaves `r`).
   * What singles out `J = 48`: its helicity bases (state 1048 and the separately computed `refb` reference) fit the rank-1 structure
     `10⁴–10⁶×` worse than every other `J`. Third singular value `7·10⁻⁵ / 9·10⁻⁶`, against `≤ 2·10⁻⁹`. Imaginary residue
     `0.057 / 9·10⁻⁴`, against `≤ 4·10⁻⁶`.
   * K₂ and direct rows share these bases, the reference, and the state normalisation constants (`logC`, `logK`). That is where the
     common offset must sit. A record inconsistency of `10⁻⁴` relative is `≈ 10²` `c₂` at `m² ≈ 1.1·10⁶`. The precise mechanism is
     **[I]**.
4. **Curvature double count fixed** (`sm_curvfix.py` + `sm_chainfix.py`; new outputs only) [NE; §4].
   * The rows' `Γ` contain the leading-Regge curvature in **partial-wave phase**. The κ-image
     `m²{A D_h(x_e) − D_h(x) + 2[P − A J₀(νβ)]/env}` matches with slope `1.0 ± 0.1` up to `x = 0.6` (`J = 40`). The hybrid's own
     Bessel-phase image `2m²(P − J₀)/env` drops out of phase at `x ≳ 0.35`.
   * The fix subtracts the κ-image from every row, after the Hampel filter so that the node set is unchanged, and from z4's
     non-local rows. Its integrated content equals `cv₂ + cv₃` to `< 10⁻³ N_s`.
   * **Corrected minima** (λ/N_s):
     * **binding states unchanged**: 1000 `+0.2021` (`J = 5`, `+−`, `1.01Λ`), 300 `+0.2052` (`J = 5`, `1.02Λ`);
     * 1000, large `J` (ee1000 shape chain): `J = 32` `+0.9036` (was 0.9038), `J = 36` `+0.9669`, `J = 40` `+0.9794` (K₂ and direct),
       `J = 17` `+0.8389` (was 0.8393). The w2 spline chain gives `J = 32` `+0.9245` (was 0.9535);
     * 300, z4 main chain (direct small-ν nodes, min-norm gauge): `J = 24` `+0.8204` (was 0.8446), `J = 28` `+0.8133` (was 0.8865),
       both at 1.5Λ; `n = 1`.
   * **Side finding** [NE]: the chains' representation of an oscillating function of the curvature's size is itself uncertain by
     more than the double count (`sm_alias.py`):
     * w2/q3 spline at 300: `+0.08 / +0.28 / +0.05–0.16 N_s` at `J = 24 / 28 / 32`, near the gap;
     * ee1000 shape fit with the Bessel image: up to `−0.9 N_s` (`J = 40`) and `−1.9` (`J = 48`); with the κ-image only
       `−0.005…−0.019`.
     * So the κ-image subtraction also removes most of the representation error. The spline chain's values at 300, `J ≥ 24`, carry a
       `±0.1–0.3 N_s` systematic.

**Grade of the headline.** The smooth part in the **exact scalar** is `≤ 0.08` `c₂`-equivalent at every computed large-`J` gap state,
`10²–10⁴` below tolerance [NE, exact]. For the **graviton**, the smooth part is measured to be consistent with zero:
* `≤ 6–36` `c₂` in the bulk at `J ≤ 40`;
* `+3…+86` (`≥ −52` at 2σ) in the fold window at `J = 24, 28` (300).

With the [D-sketch] localisation (fold/edge only), that suffices for `J ≤ 300` at 1000 and `J ≤ 64` at 300, **conditional** on:
* the [D-sketch] holding (no bulk smooth term at `J > 40`);
* the graviton fold-window smooth value staying above `−200` (1000) / `−30` (300, `J = 100–150`) `c₂` as `J` grows. A positive
  value is harmless;
* `ν < 10` and the edge `x > 0.8`, where there are no graviton data.

`n = 1` at large `J` near the gap therefore stays **[NE + M], conditional, with [I] narrowed** to these three items.

## 1. Analytic: sources and weights of smooth terms

### 1.1 Classification [D, from SUBLEADING_REGGE §3]
* Write `Γ(x) = (2J+1)²[c(x)D_h + s(x)xD_h′/(2J+1)] + r(x)`, with smooth `c, s, r`.
* After the `ψ`-integration the first two terms contribute only through the `p → 0` (Regge) endpoint: `6·10⁻³–8·10⁻² N_s` per unit
  `c`, `≤ 5·10⁻⁴` per unit `s`.
* `r` contributes with the full weight `W_all = Σψ₂ env r/(m²N_s)`, i.e. the tolerance table.
* The exact-scalar `++` totals (`≈ 4·10⁻⁶ N_s` at every `J`) are this endpoint term for `c₀ = 6.5`
  (`6·10⁻³·6.5/(2J+1)² ≈ 6·10⁻⁶` at `J = 40`).

### 1.2 Crystal saddle: every order oscillates [D-sketch]
* The heavy action is a 7-variable crystal integral `∫ e^{Φ} 𝒫_J`. Here `𝒫_J` is the spin-`J` prefactor: `Q_J = |Z|^J|W|^J P_J(cos)`
  for the scalar, and the polarisation contractions for the graviton.
* Every term of the saddle expansion around the crystal saddle is linear in `𝒫_J` and its derivatives at the saddle. These are
  `e^{±iℓθ}` times polynomials in `ℓ`.
* So every order, and the `J²/m` resummation, lies in the `c, s` class with `θ → θ_eff`.
* The `x`-dependence of the saddle's position and Hessian enters `c, s` smoothly. It cannot produce a term without `e^{±iℓθ}`.
* Non-oscillating terms therefore need a point where the saddle expansion fails:
  * (a) the fold, where the crystal saddle coalesces with the σ-conjugate pair;
  * (b) the light-cone edge, where the saddle reaches the boundary of the integration domain at `x → 1`.
* What is not shown: this argument does not prove that no third structure (e.g. a non-crystal saddle that is subdominant but
  non-oscillating) contributes at `x ∈ [0.1, 0.8]`. The exact scalar finds none at `≤ 0.1` `c₂` [NE]. The graviton rows find none
  at `≤ 6–36` `c₂` for `J ≤ 40` [NE].

### 1.3 (a) The fold [D-sketch + NE]
* `ν* ≈ 11.3–11.5` is independent of `J` and `Δ` (EE_DG300_SMALLNU §1). The soft 4×4 block of `Φ_ν` is `Δ`-independent to 0.2 %
  and `J`-independent [NE].
* Hence the soft integral near the fold has **no large parameter**. The uniform (Airy) asymptotics has width `δν ~ O(ν*)`, not
  `m^{-2/3}`. The fold contributes an `O(1)` function `g(ν)` on `ν ≲ ν_F ≈ 2ν*`, with no `m`- or `J`-scaling of its width.
* The function itself is analytic in `ν`: the fold is a property of the representation. The exact scalar confirms this.
  `m²ρ₂(ν → 0) = 6.5–6.9` `c₂` at `J = 24, 40, 48, 64, 100` and both `Δ_g`. It falls smoothly through zero where `D_h` does,
  i.e. at `ν ≈ 12–25` at 1000 and `ν ≈ 6–12` at 300 for `J = 40–64`. It shows no structure at `ν*` (probes at `ν = 9, 10, 11, 12, 13`).
* Weight of the window `ν < ν_F`, with `ψ₂ ≃ X₀p` and `env(0) = 2`:
  `W_F = Σ_{ν<ν_F} ψ₂ env/(m²N_s) ≃ X₀ν_F²/(Λ²m²N_s)`. With `N_s ≈ 1024/b̂³` this is `(X₀/1024)ν_F²(2J+1)³Λ/m⁵`.
  At the gap (`m ≈ Λ√(1 + 2J/Δ_g)`) it grows like `J³/Λ⁴` for `J ≪ Δ_g` [D].
* Numbers: `sm_bound.py` → `cache/sm_bound.txt`. The analytic `W_F` agrees with the quadrature to 4–5 % (1000) and 12–17 % (300).
* Inside the first lobe (`(2J+1)ν/m < 2.405`) a fold term cannot be told apart from `c₀D_h`. The question is only its size times
  `W_F`.
* Measured on the graviton:
  * Direct all-orders values at `Δ = 336`, `ν = 10–33`, minus the κ-curvature image, fit `c₀D_h + r` with `c₀ = 427 ± 83 / 483 ± 122`
    and `r = +3 ± 28 / +44 ± 22` `c₂` (`J = 24 / 28`, `++`). `+−` gives `+13 ± 29 / +49 ± 30`. With a `c₂x²D_h` term added,
    `J = 24` gives `r = +67 ± 36 (++) / +86 ± 33 (+−)` (`sm_foldfit.py` → `cache/sm_foldfit.txt`).
  * So the fold-region correction is dominantly `∝ D_h` (positive `c₀`). Its smooth part is `≲ 10²` `c₂`, positive or zero
    (`≥ −52` at 2σ).
  * A positive smooth term raises `λ`, so only the negative side matters for `n`.
  * `ν < 10` is measured only by thimbles, at `±5·10³` (EE_DG300_SMALLNU).

### 1.4 (b) The light-cone edge [D + D-sketch + NE]
* The `p`-support ends at `p = 1`, i.e. `x_max = Λ/m`. At the gap,
  `1 − x_max² = (m² − Λ²)/m² = 2J(Δ_g − 2)/m²`, i.e. `0.075` (`J = 40`, 1000) to `0.38` (`J = 300`), and `0.14` to `0.5` at 300.
* The light cone itself (`x = 1`, `ν = m`) is outside the support by `Δν = m − Λ ≈ J` [D].
* `ψ₂ = (1−p)⁵pX(p)` with `X(1) = −8648`, and `ψ₃` likewise vanishes to fifth order at `p = 1` [D, from nf_common].
* Exact scalar near the edge [NE]: the remainder oscillates, with envelope `m²ρ_k ≈ (200–700)(1−x²)^{-3}` for both `k = 2, 3`
  (`J = 40, 48, 64`).
  * The `k = 2` and `k = 3` edge contributions are each `0.2 N_s` in `+−`. They cancel in the `σ = −1` combination to `≤ 4·10⁻⁴`.
  * In `++` (`fac = (1−x²)²`) each is `≤ 10⁻³`.
  * This cancellation is a property of the NF02 weights acting on the scalar. For the graviton it would hold if its `k = 2, 3` edge
    envelopes stay equal [I].
  * Even without it, `0.2 N_s` is below the margin `0.9`.
* An oscillating edge term with `ψ ~ (1−p)⁵` contributes `≈ 5!|X(1)|·env(x_max)·A/(m²φ′⁶N_s)`, `φ′ = b̂/√(1−x_max²)`
  (`sm_bound` `osc-edge`). That is `≤ 3·10⁻⁸ N_s` per `c₂` of amplitude at 1000, and `≤ 8·10⁻⁷` at 300 [D].
* A smooth edge term would need the crystal saddle to reach the boundary. That happens at `x = 1`, outside the support, and its
  leakage to `x_max` is not derived [I]. The exact scalar sees none above `10⁻⁴ N_s` [NE].
* Tolerance on a smooth edge-shaped term `A_E((1−x_max²)/(1−x²))³` (amplitude at `x_max`): see §0.2(b).

### 1.5 The bound
`|dF_sm/N_s| ≤ W_F(ν_F)·|Γ_F| + W_E·|A_E| + W_all·|r_bulk|`, with the [D-sketch] setting `r_bulk = 0` and the weights in `sm_bound.txt`.

Against the data:

| `Δ_g`, `J` | fold tol (`ν < 30`) vs graviton `r_F ≥ −52` (2σ) | uniform tol vs graviton `|r_bulk|` (measured `J`) | scalar `C_eq` |
|---|---|---|---|
| 1000, 40 | 12800 ✓ | 46 vs `≤ 12` ✓ (`J = 40`) | `10⁻³` |
| 1000, 64–100 | 3580–1170 ✓ | 14–5 vs — (not measured) | `≤ 10⁻³` |
| 1000, 150–300 | 490–205 ✓ | 2–1 vs — | (not computed) |
| 300, 24–40 | 653–186 ✓ | 25–8 vs `≤ 36` (✓ at 24–32 with 1σ) | `≤ 0.08` |
| 300, 64 | 70 ✓ | 3 vs — | `1.3·10⁻³` |
| 300, 100–150 | 36–27: not demonstrated (2σ bound −52) | 2–1 vs — | `3·10⁻⁴` (`J = 100`) |

## 2. Exact scalar skeleton at large `J` near the gap [NE, exact]

* **Why `sr_scalcmp` failed** (`sm_nancheck.py 1040 40`): `Qpoly_cheb` returns `max|c_i| = 4·10⁷¹` (`n = 0`), `2.5·10¹⁰²` (`n = 399`),
  `2.6·10²⁶²` and then `inf` on the tail nodes (`n ≈ 7·10⁶–5·10⁹`). That is overflow in the conversion to doubles, and in the
  conditioning test `V @ cf`.
  * Rescaling per node would remove the overflow. It would not remove the large-`J` n-sum cancellation, which needs Arb weights
    (h13_laplace header).
  * So the Arb Laplace route is used: `sm_exact.py` → `cache/sm_exact_<Dg>.pkl`; analysis `sm_exact_an.py`.
* **`Δ_g = 300` n-sum** (`sm_nsum.py` → `cache/sm_nsum_364_64.txt`, `cache/sm_nsum_400_100.txt`):
  * At `(364, 64)` and `ν/Λ ≤ 0.6` every configuration agrees to `10⁻⁹`.
  * At `0.85, 0.95` the values depend on `N₀`: 200 gives `−44.8`, 320 gives `+3·10⁴`, 600 gives `2·10¹⁰`. This is unchanged by
    precision (800 bits), the w-series length (×3) and the Euler–Maclaurin order.
  * At `(400, 100)`, `N₀ = 20 / 40 / 60` agree to 10 digits at `ν/Λ = 0.95`, also with `nq = 1920` and `L/4`. `N₀ = 80` deviates at
    `4·10⁻⁴` and `120` gives `−104`.
  * So the explicitly summed terms at `n ≳ 60–80` are inaccurate at large `ν`, while the Euler–Maclaurin tail started early is not.
    The cause is not identified [I].
  * Production at 300 uses `SM_NCONF=20:960:1,60:960:1,40:960:1` (tag `n0t`). The quoted n-sum error is the spread of these three.
* **Records** (`m²ρ₂` in `c₂` units). `rem` uses the route-A′ weights (`k = 2, 3`; `k = 3` by quartic interpolation in `k` [M],
  odd-`k` sensitivity `≤ 10⁻⁴ N_s`):

  | state | `++`: rem/N_s (k2, k3) | `+−`: rem/N_s (k2, k3) | n-sum err | fold plateau `m²ρ₂(ν ≤ 10)` |
  |---|---|---|---|---|
  | (1040, 40) | `+3.8e-6` (`+4.0e-6`, `−2.3e-7`) | `−2.5e-5` (`+6.3e-3`, `−6.4e-3`) | 0 | 6.5 |
  | (1048, 48) | `+4.1e-6` | `+3.9e-4` (`+3.4e-3`, `−3.0e-3`) | 0 | 6.5 |
  | (1064, 64) | `+3.5e-6` | `+2.9e-5` (`−1.7e-3`, `+1.7e-3`) | 0 | 6.5–6.4 |
  | (1100, 100) | `+4.4e-6` | `−1.5e-4` (`−1.3e-4`, `−1.9e-5`) | 0 | 6.6–5.4 |
  | (324, 24) | `+1.1e-6` | `+2.8e-3` (`−5.5e-3`, `+8.3e-3`) | `10⁻¹¹` | 6.5–6.0 |
  | (340, 40) | `+3.6e-5` | `−1.7e-4` (`−2.0e-3`, `+1.8e-3`) | `10⁻¹⁵` | 6.6–3.1 |
  | (364, 64) | `+2.6e-5` | `−3.6e-4` (`−9.0e-4`, `+5.4e-4`) | `2·10⁻¹²` | 6.9–(−5.9) |
  | (400, 100) | `+4.0e-5` | `+1.5e-4` (`+1.2e-4`, `+3.0e-5`) | `6·10⁻⁹` | 7.2–(−10.7) |

* Band decomposition (`sm_exact_an`): the pieces on `x ∈ [0.1, 0.4]`, `[0.4, 0.7]`, … are `±0.1–5 N_s` and cancel to the totals above.
  * So the remainder is large and oscillating, and integrates out.
  * A "smooth projection" by least squares is unstable. It gives `±0.005–0.03 N_s`, against true totals `≤ 4·10⁻⁴`, so it is
    not used (as SUBLEADING_REGGE found).
* p-rule check: `(1040, 40)` with 2× p-nodes (219 instead of 141; tag `npf2`) reproduces every total to all printed digits.
* Output: `cache/sm_exact_an.txt`. The keys `n0s` (`N₀ = 80–160`) are superseded at 300 `J ≥ 64` by `n0t`.

## 3. The `J = 48` offset (`sm_offset.py` → `cache/sm_offset.txt`) [NE]

Fit on `x ∈ [0.1, 0.45]`: `Γ − κ-image = (c₀ + c₂x²)D_h + (s₀ + s₂x²)xD_h′/(2J+1) + r`.

| rows | `r` (C2 `++`, C2 `+−`, C3H `++`, C3H `+−`) | residual / rms |
|---|---|---|
| 1000: `J = 32` (1040), 36 (1036), 40 (1040 K₂ and direct) | all within `0 ± 6` | 2–10 % |
| 300: `J = 24, 28, 32` (336; `J²/m` = 1.7 / 2.4 / 3.1), `J = 24, 28` (500) | `0 ± (2–36)` | 2–40 % |
| exact scalar `(1048, 48)`, same nodes | `+0.7 ± 3` | 7 % |
| **`J = 48` (1048) K₂** | **`+279, +620, +275, −633`** (±120–230) | 70–80 % |
| **`J = 48` (1048) direct** | **`+363, +624, +361, −689`** (±90–230) | 55–75 % |

* Concentration: `Γ − κ` is `+400…+2400` at `x = 0.05–0.15` (`ν = 50–160`) and `+200…+400` at `0.175–0.43`. On `x ≥ 0.45` it swings by
  `±2000–4000`: the κ phase model fails there, but in the oscillating class.
* Excluded explanations:
  * the curvature double count: `r` is identical for `Γ`, `Γ − κ` and `Γ − Bessel` image, because the images lie in the oscillating
    class;
  * the fold: `ν ≥ 50 ≫ ν*`;
  * the edge: `x ≤ 0.45`;
  * aliasing or the fit: the exact scalar on the same nodes is clean;
  * the de-shift: C3H shows it too;
  * `J²/m`: the 300 rows up to 3.1 are clean;
  * the reference LO zero (`x₀ = 0.474`, node 0.475 despiked): a `D_h/D_ref` basis leaves `r` unchanged.
* What is specific to `J = 48` (`L.basis` quality, as l8: third singular value of the rank-1 helicity fit / imaginary residue of the
  basis vectors):

  | `J` | state row | reference `Δ = 2·10⁴` |
  |---|---|---|
  | 36 | `5·10⁻¹¹ / 3.6·10⁻⁶` | `8·10⁻¹² / 10⁻⁹` |
  | 40 | `2·10⁻⁹ / 4·10⁻⁶` | `5·10⁻¹⁰ / 3·10⁻⁸` |
  | **48** | **`7·10⁻⁵ / 0.057`** | **`8.5·10⁻⁶ / 8.8·10⁻⁴`** |

* The `J = 48` reference basis is built from the `J = 48` reference records, which include the separately computed `refb` set
  (ee1000_rows). The state basis is built from the `J = 48` K₂ records. K₂ and direct rows share both.
* Order of magnitude: in `Γ` units a relative inconsistency `δ` of the records is `m²δ ≈ 1.1·10⁶ δ` `c₂`. The `J = 48` rank-1 misfit
  (`7·10⁻⁵` state, `9·10⁻⁶` reference, relative) and the basis-phase residue (`0.057`) correspond to `10²–10³` `c₂`, the size of the
  offset. At `J ≤ 40` they correspond to `≤ 10⁻²` `c₂` [D-sketch].
* The mixing route alone is too small: `Im Γ_off` at `J = 48` has rms `0.05/m`, so `ε·m Im Γ_off·m ≈ 6` `c₂`. The precise mechanism
  stays **[I]**.
* **Verdict: artefact of the `J = 48` rows** [NE for the exclusions and the singled-out bases; the mechanism is D-sketch]. It is not
  evidence for a physical smooth term. The `J = 48` rows should be recomputed (records and `refb`) before any `J = 48` verdict.

## 4. Curvature double count: corrected chains [NE]

* `sm_curvcheck.py` → `cache/sm_curvcheck.txt`: slope of row `Γ` on the image, per band `[0, 0.2)`, `[0.2, 0.4)`, `[0.4, 0.6)`,
  `[0.6, 0.8)`, at `J = 40` (1040), C2 `++`:
  * κ-image: `1.01 / 1.03 / 0.88 / −0.48`, residual `89 / 182 / 421 / 3019`;
  * Bessel image: `1.01 / 0.69 / −0.31 / 0.51`.

  C3H and `+−` are alike. The rows contain the κ (partial-wave-phase) curvature.
* `sm_curvfix.py`:
  * `SM_CURV=row` subtracts the image (`SM_CURVIMG=kappa` default, or `bessel`) from every diagonal row `Γ`, even (`l9.grid` after
    w2's despike) and odd (`P4.grid`), and from z4's NL rows.
  * `SM_CURV=none` reproduces the old outputs exactly: w2 1000 and 300, z4 main `+0.8446 / +0.8865`, ee1000 f1/fd40/d48. The w2 300
    `J = 32` line differs only because the 300 agent's new `J = 32` 336 row is now on disk.
  * `SM_CURV=only` replaces `Γ` by the image, which gives the chains' representation error (`sm_alias.py`).
* Corrected runs (`sm_chainfix_all.sh` → `cache/sm_chainfix_<mode>_<tag>.txt/.pkl`):

  | chain | state | old | corrected (κ, row) | Bessel image (row) |
  |---|---|---|---|---|
  | w2 rerun, 1000 | binding `J = 5` | +0.2021 | **+0.2021** | +0.2021 |
  | w2 rerun, 1000 | `J = 32` (1.05Λ) | +0.9535 | +0.9245 | +0.911 |
  | ee1000 shape (f1), 1000 | `J = 32 / 36 / 40 / 17` | +0.9038 / 0.9668 / 0.9795 / 0.8393 | **+0.9036 / 0.9669 / 0.9794 / 0.8389** | +0.924 / 0.976 / 0.982 / 0.839 |
  | ee1000 direct (fd40) | `J = 40` | +0.9794 | **+0.9794** | +0.9824 |
  | z4 main, 300 | binding `J = 5` | +0.2052 | **+0.2052** | +0.2052 |
  | z4 main, 300 | `J = 24 / 28` | +0.8446 / +0.8865 | **+0.8204 / +0.8133** (1.5Λ) | +0.8195 / +0.8201 |
  | z4 main, 300 | 1.1Λ, `J = 24 / 28` | +1.056 / +1.446 | +1.013 / +1.570 | |

* The ee1000 Jacobi fix (`SM_JACFIX=1`) changes z4 `J = 28` at 1.1Λ by `−0.04` and the minima not at all.
* Representation errors (`sm_alias`, `(model − cv)/N_h` at the gap state):
  * w2 spline: 1000 `J = 32` `+0.05`; 300 `J = 24 / 28 / 32` `+0.08 / +0.28 / +0.05–0.16`, against `cv/N_h = 0.017 / 0.026 / 0.033`;
  * ee1000 shape: κ-image `−0.005 / −0.007 / −0.010 / −0.019` (`J = 32 / 36 / 40 / 48`), Bessel image `−0.3 / −0.36 / −0.92 / −1.9`.

## Hand-over (kept current)

* **Nothing is running.** All `sm_*` jobs have finished. The two `ee1000_queueR1.sh` queues of the previous agent were left alone.
* **Reproduce:**

  | command | output |
  |---|---|
  | `python3 sm_exact.py 1000 40,48,64,100` | exact scalar remainder (5 min/state) |
  | `SM_NCONF=20:960:1,60:960:1,40:960:1 SM_KEYTAG=n0t python3 sm_exact.py 300 24,40,64,100` | same at 300 (`n0t` key) |
  | `python3 sm_exact_an.py 1000 300` | `cache/sm_exact_an.txt` |
  | `python3 sm_foldfit.py 10 40` | `cache/sm_foldfit.txt` |
  | `python3 sm_bound.py 1000 40,48,64,100,150,200,300` and `python3 sm_bound.py 300 24,28,40,64,100,150` | `cache/sm_bound.txt` |
  | `W2_DMIN=1000 [EE1000_DSTATE=1 EE1000_JEVEN=40,48] python3 sm_offset.py "40:1040,48:1048"` | §3 |
  | `./sm_chainfix_all.sh` (≈ 25 min) | §4 corrected chains |
  | `SM_CURV=none ...` | reproduction of the old chains |

* **Next steps (value order):**
  1. Graviton fold window at `J ≥ 64` (300) and `J ≥ 150` (1000). It needs `|r_F| ≲ 30` / `200` `c₂` on `ν < 30`. Direct n3 nodes at
     `ν = 10–33` on one large-`J` row would decide whether `r_F` stays `≈ 0` as at `J = 24, 28`.
  2. Recompute the `J = 48` rows (records and `refb`). Check the basis quality (`sv3`, `imag`) first.
  3. Find the cause of the 300 n-sum `N₀`-dependence at `ν/Λ ≳ 0.85` (`sm_nsum.py`). Small `N₀` works, but it is not understood.
     This matters for H13′ at 300 too: its `N₀ = 200/320` records at `J ≳ 100`, `b̂ ≥ 150`, e.g. `(1.5, 112)`, show
     `S/N_s = −476 ± 500`.
  4. A proof of §1.2, with the edge leakage of §1.4.
