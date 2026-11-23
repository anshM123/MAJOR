# k = 2 de-shift weights (ads4/prl/data/k2)

Written 23 Nov 2026 by `ads4/code/cert/k2_export.py` (see `ads4/CERT_CLOSE.md` §3). New files only; nothing in `prl/data/` was
changed. Checksums: `SHA256SUMS` (full SHA-256, `sha256sum -c SHA256SUMS`). Conventions as in `../README.md` §4: `Λ = √((Δ_g−2)(Δ_g−1))`,
`p = ν/Λ`, kernels `R_∓` = `../kernels.json` → `D_k2_EE`, χ-normalisation `X = χ(R₋ + R₊)`.

## 1. `X^Q`: Result Q's k = 2 Abel weight on its grid (previously unpublished)

| file | content |
|---|---|
| `xq_abel_X_300.txt`, `xq_abel_X_1000.txt` | columns `ν_k`, `Re X`, `Im X`, `Re Y`, `Im Y`, `ψ₂^Q(ν_k)`; `ν_k = (k+½)dν`, `dν = 0.02` (to 327.68) / `0.04` (to 1310.72) |

* **Definition.** Min-norm solution of `C[aY] − i tanh t·S[bY] = C[ψ₂^Q]/cosh t`, `X = aY` (same equation, kernels, `n` and gauge as Result
  D's `../d_abel_X_*.txt`), with `ψ₂^Q(ν) := ψ₂(ν/Λ)` from `../q3_full_delta0.2.json` (no `1/Λ`; the smeared piece is
  `(1/Λ)Re Σ_k dν X_k Ψ̂_{EE,ν_k}/(h₋+h₊)(ν_k)`).
* **Source.** `ads4/code/resultq/cache/rv_k2eng_Y_{300,1000}.pkl` (RESULT_Q_REVISION §1.3; `u3_chi_EE.SolverAB`, LSQR residual
  `1.0·10⁻¹⁵` / `3.9·10⁻¹⁶`), exported unchanged.
* **Independent re-solve** (`code/cert/k2_xq_check.py`: DESHIFT's `d_abel` solver, own grid `dν·5/4`, FFT, CGLS) [DV]:
  `max|X_new − X_stored|/max ψ₂^Q ≤ 2.1·10⁻⁷ / 8.3·10⁻⁷` (ν < 2), `≤ 3·10⁻¹² / 1.3·10⁻¹¹` (25 ≤ ν ≤ Λ); graviton-pole identity V3 exact to
  `1.3·10⁻¹⁴ / 2.1·10⁻¹⁴` (stored); `dk₀^Q = −1.630146574·10⁻⁴ / −1.453689613·10⁻⁵` (agreement `3·10⁻¹³ / 8·10⁻¹¹`).
* Like D's grid weight, the stored double-precision grid satisfies the de-shift identity only to `t ≈ 18–22` (off-grid residual
  `~10⁻⁴` at `t = 22`). For all-`t` control use the continuum weights of §2.

## 2. Continuum k = 2 weights (Result D and Result Q), for all t

| file | weight |
|---|---|
| `k2_continuum_D_300.json`, `k2_continuum_D_1000.json` | `Ω^EE` (Result D): target `ψ = ψ₂^{NF02}(ν/Λ)/k₀`, `k₀ = 27.48509817` |
| `k2_continuum_Q_300.json`, `k2_continuum_Q_1000.json` | `Ω^Q_up`'s smeared k = 2 piece: target `ψ = ψ₂^Q(ν/Λ)` |
| `k2_continuum_*_table.txt` | `ν`, `Re χ`, `Im χ`, `Re X`, `Im X`, `ψ(ν)` on the stored grids (convenience only; the JSON formula is the definition) |

* **Form.** `χ = χ_re + iχ_im`, each part
  `Q(ν) Σ_k x_k [g_k(ν−c_k) + g_k(ν+c_k)] + Σ_{j=0}^{2} y_j ν^{2j} l(ν) e^{−ν²/2}`, with
  `Q(ν) = ν²/(ν² + ¼)²`, `g_k(z) = exp(−z²/2σ_k²)` (centres `c_k`, widths `σ_k` listed: `0.15/0.45` spacing/width to ν = 40, `0.3/0.9`
  beyond, to `Λ + 30`), `l(ν) = (ν²+1) log(ν²+1)` (principal branch). Coefficients are IEEE doubles (`repr`); double-precision
  evaluation agrees with 50-digit evaluation to `≤ 1.1·10⁻¹²` of `max|χ|`.
* **Singularity structure.** Double poles at `±i/2` (the double zeros of `R₋` at `−i/2` and of `R₊` at `+i/2`, the only common zeros the
  k = 2 shifts cross); double zero at 0 (cancels the `1/ν` poles of `R_∓`); branch points at `±i` (they produce the `|μ|` kink of the
  target, `ψ′(0) ≠ 0`, which no weight analytic on the closed strip can); analytic in `|Im ν| < 1` otherwise; super-Gaussian decay.
  Without the `±i/2` poles the fit does not reproduce the stored weight (`max|Δχ|/max χ ≈ 1`); with them it does (below).
* **Identity and its bound.** `W[χ](t) = ∫₀^∞ χ[R₋cos((ν+i)t) + R₊cos((ν−i)t)]dν = ½∫_ℝ E(μ)cos(μt)dμ` for every real `t`, with
  `E(μ) = 2Re[χ(μ−i)R₋(μ−i)]` (shifted form; no residue terms). Hence, for **all** `t ≥ 0`,
  `|W[χ_re](t) − Φ(t)| ≤ ε_W := ∫₀^∞|E − ψ|` and `|W[χ_im](t)| ≤ ε_W^{im}` (`χ_im` is a null vector), and
  `|image − target| ≤ ε_W P_{−½}(cosh β)` for all β.

  | weight | `ε_W` (χ_re) | `ε_W/∫|ψ|` | `ε_W^{im}` | `max|X − X_stored|/max ψ`, ν < 2 / 25–Λ | `dk₀` (continuum / stored) | V3 rel. error |
  |---|---|---|---|---|---|---|
  | D, 300 | `3.3·10⁻⁶` | `1.1·10⁻⁹` | `4.1·10⁻⁷` | `3.8·10⁻⁷ / 1.7·10⁻¹⁰` | `−7.666324·10⁻⁴` / same (`4·10⁻⁸`) | `6·10⁻¹¹` |
  | D, 1000 | `1.1·10⁻⁶` | `1.1·10⁻¹⁰` | `1.3·10⁻⁷` | `4.3·10⁻⁷ / 5.6·10⁻¹¹` | `−6.830150·10⁻⁵` / same (`4.5·10⁻⁷`) | `8·10⁻¹⁰` |
  | Q, 300 | `1.9·10⁻⁷` | `6.2·10⁻¹⁰` | `1.8·10⁻⁸` | `2.0·10⁻⁷ / 1.4·10⁻¹⁰` | `−1.630147·10⁻⁴` / same (`1·10⁻⁸`) | `8·10⁻⁹` |
  | Q, 1000 | `7.4·10⁻⁸` | `7.3·10⁻¹¹` | `5.6·10⁻⁹` | `2.2·10⁻⁷ / 2.1·10⁻¹¹` | `−1.453690·10⁻⁵` / same (`2·10⁻⁸`) | `3.5·10⁻⁷` |

  (V3: graviton-pole identity `∫χ[R₋D(ν+i) + R₊D(ν−i)] = ∫ψD`, `D(z) = 1/(z²+9/4)`; its error is bounded by `(4/9)ε_W`.)
* **Grade.** `ε_W` is a double-precision quadrature of `|E − ψ|` on `4.2·10⁵` points to `Λ + 60` (beyond, the Gaussian tail is
  `< 10⁻¹⁰⁰`): **[DV]**, not interval-verified (as DESHIFT's C4/C5 `ε_W`).
* **Relation to the stored weights.** The real part equals the stored min-norm grid weight to its own discretisation accuracy
  (`≤ 4·10⁻⁷` of `max ψ`, below the `1.5·10⁻⁶` spread between the two independent grid solvers of DESHIFT §4.1); the imaginary part (a
  null vector of the identity, 8 % of `max|Re χ|` at small ν) is reproduced to `≤ 4·10⁻⁹`. So these are the same functionals as the
  grid weights, now with uniform-in-`t` control. The cubic bound is unchanged: `(1 + dk₀)^{−½} = 1.00038354 / 1.00003415` for both
  weights, `R₀(K = 1) = 3.13206·10⁻⁶ / 2.89126·10⁻⁸` (weight-independent).
* **Cutoff.** Multiply by the erf plateau `θ` of UNIFIED_THEOREM_v2 §2.1 with **D's** parameters (`ν_c = Λ + 100`, `a = 10`) for both D's
  and Q's k = 2 piece (these weights extend to `ν ≈ Λ`; Q's `ν_c = 150` is meant for the C4/C5 weights and would truncate `X^Q`). The
  weights are `< 10⁻¹⁴` of their maximum beyond `Λ + 25` and decay super-Gaussianly beyond `Λ + 30`, so `θ` changes nothing computed.
