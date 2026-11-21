# Result Q, Stage 5: the spinning `Ĉ₅` — exact raw image, the 3-shift de-shift, validity, and large-`J` positivity

*21 Nov 2026, Stage-5 agent. Unrefereed, not committed by this agent. No existing document or code file was edited. New code:
`ads4/code/resultq/stage5_hj.py`, `stage5_check.py`, `stage5_deshift.py`, `stage5_positivity.py`, `stage5_gauge.py`; kept functional file
`ads4/code/resultq/functionals_q5_deshift.json`; outputs in `resultq/cache/stage5_*` (git-ignored). Read-only imports:
`h7i14/v3_regge_kernel.py`, `v4_hclosed.py`, `spinfunc/q1_c3_regge_kernel.py` (the existing `j′ = 1` light-ray side and `o_H q` dual side,
which Stage 4 had overlooked), `resultq/stage4_hj.py`, `stage12_forward_regge.py`. Compute: about 15 CPU-minutes on ≤ 2 cores.*

Grades: **[P]** proved (incl. exact computer algebra), **[D]** derived with the stated assumption, **[DV]** derived and checked numerically,
**[M]** model / numerical, **[NE]** numerical evidence, **[I]** open. Flat = D = 4 amplitudes; AdS = the CFT₃ functionals.
`𝒫(μ) := P_{−1/2+iμ}(cosh β)`, so `𝒫(ν − ij) = P_{l+j}`, `l = −½ + iν`. Classes: `++ := (h₁,h₂) = (+,−)`, `+− := (+,+)`.

---

## 0. Answer in brief

1. **Exact raw image [P + DV].** The leading-Regge heavy image of the raw spinning `Ψ_{5,ν}[o_H q]` (`j′ = 1`, `Ĉ₃` structures) is
   `N(ν) Σ_{j=0,±2} R_j(ν) 𝒫(ν − ij)` with `R_j` exact rationals (§1.3). Only `j = 0, ±2` occur.
   * **Carries `σ_h` exactly:** `R_j(+−) = −R_j(++)` for every `j`, and both classes are invariant under the total flip `(h₁,h₂) → (−h₁,−h₂)` [P].
   * **Flip-free:** a helicity flip on either line makes that side's rank-0 and rank-2 tensor vanish identically, so all off-diagonal
     terms are 0 [P]. Numerical check: `≤ 2·10⁻¹³` [DV].
   * The `n`-odd (ε-carrying) rank-2 components vanish on both sides [P].
   * Checks:
     * rank 0 reproduces q1's numerical `c₃` exactly (`10⁻¹⁵`);
     * the rank-2 scalar (`j′ = 0`) reproduces Stage 4's proxy `a⁽²⁾_j E₁₀(ν−ij)` to `3·10⁻¹⁶`;
     * an independent finite-difference evaluation in q1's engine at real β agrees to `≤ 1.3·10⁻⁶` (FD-limited).
2. **New structural facts.**
   * The coefficients are real-analytic (`conj[R_j/R₀](ν) = [R_j/R₀](−ν)`) but **not shadow-symmetric**. So no ν-phase makes them
     "real", unlike `k = 2, 4` and the scalar.
   * The forward normalisation `N̂₅ = Σ_j R_j` is complex and has **no real zero** on `ν > 0`: `|N̂₅|/Σ|R_j| ≥ 0.057` for `ν > 0.5`, and it
     → 0 only as `ν → 0⁺`, where the `1/ν` poles of `R_{±2}` cancel.
   * The raw `ν = 0` image is **exactly real** (normalised to 1 at β = 0):
     **`H₅ = κ₁ ∂_μP_μ|_{3/2} + κ₂ P_{3/2} + κ₃ P_{−1/2}`**, with `κ₁ = −192491/1028120`, `κ₂ = 1252043/1028120`,
     `κ₃ = −223923/1028120` [P].
     * `κ₁ < 0`, the **opposite sign to the scalar proxy** (`+5415/14984`).
     * So `H₅ ~ κ₁ βP_{3/2}` (`η^{3/2}log η`) is positive up to β ≈ 6.9 and **negative beyond**. The raw functional therefore fails in
       the **`++`** channel (proxy: `+−`).
3. **De-shift [D + DV].**
   * The restricted gauge (real χ on `Re[Ψ₅/N̂₅]`) is **not solvable**: the weights blow up (`∫|χ| = 4·10⁴ → 6·10⁹` under refinement).
     Im-only and two phase-twisted real gauges fail as well [NE].
   * With the general real sum rule `∫Re[c(ν)Ψ_{5,ν}/N̂₅]dν`, `c` complex, the Abel equation
     `C[a₀] + cosh2t C[A₂] − sinh2t S[B₂] = 1` is solved:
     * residual `7·10⁻¹⁵`;
     * `W(t ≤ 8) = 1 ± 5·10⁻¹⁰`;
     * mpmath direct image `≤ 3·10⁻¹²` (β ≤ 6);
     * light-side moments `(1, 2·10⁻¹³, 7·10⁻¹⁰)`;
     * `∫|c| = 1.1495` (small amplification);
     * converged: `ν_max` 41 → 61 changes `c` by `8·10⁻⁶`, `dν` 0.04 → 0.02 by `4·10⁻³`.
   * Null space: no spectral gap; about `N − 25` directions per grid, i.e. infinite-dimensional in the continuum (two real unknown
     functions, one equation). All null directions have zero leading-Regge image. Gauge used: min-norm.
4. **Lemmas A, B, C carry over [D].**
   * `R⁴`: `−Ĉ₅^{ds}|_{R⁴} = 0`; `g₅`: `−Ĉ₅^{ds} = g₅`; `|ĝ₃|²`: 0. All unchanged, via the moment identity checked above
     [P given (P), (F)].
5. **Validity [D].** With the same counting as D (raw `Ĉ₂` fails, `Ĉ₃` passes), the de-shifted `Ĉ₅` behaves as `κ = η^{−½}log η`.
   * It needs **`J₀ < 7/2`**, against raw `Ĉ₅`'s `3/2`: margin 2½ over KKSDZ's `J₀ ≤ 1`, the same as `Ĉ₄^{ds}`.
6. **Heavy positivity, both forward pieces de-shifted [M].** Leading-Regge model, β ≤ 30, i.e. `J` up to `10⁸–10⁹` at `m ≈ Λ`, `u ∈ [10⁻⁶, 1]`, all
   physical states with `b̂ ≥ 20`:
   * Minima:
     * δ = 0.05: **+0.206 / +0.065 / +0.042**;
     * δ = 0.2: **+0.287 / +0.221 / +0.212** (`Λ = 100 / 300 / 1000`).
   * All minima sit at the model edge `b̂ ≈ 20.7`, `m ≈ 20–26Λ`, `J ≈ 200–270`, where `Ĉ₅` (weight `u³ ~ 10⁻⁸`) plays no role.
   * At large `J`: β ≥ 2 gives min ≥ **+0.94**, β ≥ 6 gives ≥ **+0.987**.
   * With the raw spinning `Ĉ₅` both functionals go negative in **`++`** at β = 8.3 / 9.3 / 10.45, `m ≈ Λ`,
     **`J ≈ 3.2·10³ / 1.6·10⁴ / 9.3·10⁴`**, and reach `FF/Nrm → −1`.
   * The near-gap Stage-3 minima (`+0.0478 / +0.0499 / +0.0500`; `+0.1985 / +0.2003 / +0.2001`) are unchanged.
   * **`A`, `B`, `κ` unchanged.**
7. **Gauge [NE].** Natural gauges (min-norm, ν-penalised, as-real-as-possible) differ by **17–29 %** in `∫|c−c₀|/∫|c₀|`. The imaginary part
   cannot go below ≈ 25 % of `∫|c|`. Beyond-leading-Regge effects are therefore gauge-dependent at `O(1)` relative level. No canonical
   gauge is fixed [I].
   * Task 5b (the graviton near-gap `ν = 0` correction at `J = 4–8`) was not attempted [I].

---

## 1. The exact raw spinning image (`stage5_hj.py` → `cache/stage5_hj.{txt,pkl}`)

### 1.1 Kernel and method [P]
* **Kernel.** Lemma B gives `f⁻₅ = f⁻₃R` with `R → 16η²r²` at leading Regge. The q1 kernel of `Ĉ₃` is `(|x||y|)^{a−2}K₀`, per side
  `|x|^{b−3}` with `b = 8`. So the `k = 5` kernel is `16(x·y)²×` that, i.e. per side `|x|^{b−3}x_μx_ν` with the **same `b = 8`**, and
  `L₅ = +16⟨∂²Â, ∂²B̂⟩`. The relative `k = 3`/`k = 5` sign is kept: q1 drops the common sign of `f⁻`.
* **Side monomials.**
  * Dual side: q1's `A3_monomials('C3')` (`o_H q`, with the `z₀`-dependent `w_e = η·(z₀ × e)`).
  * Light-ray side: `B_monomials_j(1)` (the `j′ = 1` celestial moments).
  * Both with sympy coefficients.
* **Transform.** Each monomial `|x|^c(−z₀·x)^{−λ}(e·x)^i(f·x)^j(w·x)^k` maps to `K(c,λ)(e·∂)^i(f·∂)^j(w·∂)^kF_{c,λ}`.
  * Within each side all `K` differ by integer Γ-shifts from one reference. For `o_H q` the reference must be `k = 4 − δ`, not
    `1 − δ`, because of the odd `1/X` in `o_H`.
* **Rank-2 equivariance.**
  * Decomposition basis: `{pp, pẑ, ẑẑ, η}` and `n`-odd `{pn, ẑn}`, with `n = ε(p, ẑ, ·)`.
  * Contraction in the Gram basis `{p, p̄, ẑ_A}`, with `n_A = G⁻¹(0, V, 0)`, `n_B = G⁻¹(−rV, 0, 0)` and `V² = −det G`. Odd powers of `V`
    integrate to 0 (reflection fixing `p`, `p̄`).
  * Then `Σ C^k r^s → Σ_j R_j P_{l+j}` (stage4 `to_harmonics`).
* **Common factor.** `N(ν) = 2πK(c_A^ref, λ_A^ref)K(c_B^ref, λ_B^ref) × ⟨|s|^{1−δ}⟩` (dropped). Every forward-normalised quantity below
  is independent of `N`.

### 1.2 Checks
| check | result | grade |
|---|---|---|
| rank 0 (`k = 3`) vs q1's numerical `c₃ = h₀(++)` (8 ν values) | ratio `1 ± 5·10⁻¹⁵` | [DV] |
| rank 0: only `j = 0`; `R₀(+−) = −R₀(++)` (q1's `σ_h`) | exact | [P] |
| rank-2 scalar (`j′ = 0`) vs Stage 4's proxy `a⁽²⁾_j(l)E₁₀(ν−ij)` | `≤ 3·10⁻¹⁶` (ν = 0.3, 1.7, 5.2) | [DV]: the Stage-4 proxy with `b = 10` was the exact scalar |
| independent: q1 engine at real β, `x_μx_ν` by Richardson finite differences, 256-point `z₀` circle | `≤ 1.3·10⁻⁶`, both classes, ν = 0.4–7.5, β = 0.3–1.8 (FD-limited) | [DV] |
| free 5-harmonic fit (`j = −2…2`) of the numerical image | `|h_{±1}|/|h₀| ≤ 10⁻³` (conditioning-limited) | [NE] |
| `n`-odd rank-2 coefficients, `V`-odd contraction terms | identically 0 | [P] |

### 1.3 Result [P] (`++`; `+−` = −(`++`); `N` dropped)
```
R_{−2} = −i(2ν−5i)(2ν+5i)(2ν+9i)(2ν+13i)(8ν³+52iν²−22ν+153i)(32ν⁵+1072iν⁴−10000ν³−29240iν²+20922ν−1809i)
         / [4ν(ν+i)(2ν−13i)(2ν−9i)(2ν+17i)]
R_0    = −i(2ν−9i)(2ν−5i)(2ν+5i)(2ν+9i)(256ν⁸+6656iν⁷−35200ν⁶+51328iν⁵−502656ν⁴−282912iν³−971944ν²−481640iν−671769)
         / [2(ν−i)(ν+i)(2ν−13i)(2ν+13i)(2ν+17i)]
R_{+2} = −i(2ν−13i)(2ν−9i)(2ν−5i)(2ν−3i)(2ν+i)(2ν+5i)(4ν²+60iν+67)(8ν³−52iν²−22ν−153i) / [4ν(ν−i)(2ν+13i)(2ν+17i)]
k = 3:  R_0^{(3)} = −4i(2ν−5i)(2ν−3i)(2ν+i)(2ν+3i)(2ν+5i)(4ν²+60iν−53) / [(2ν−13i)(2ν−9i)(2ν+13i)(2ν+17i)]
```

### 1.4 Structure
* **Helicity [P].** `R_j(+−) + R_j(++) = 0` for `j = 0, ±2`, while `R_j(+−) − R_j(++) ≠ 0`. So the image **carries `σ_h`**; it is not
  helicity-blind. The total flip leaves both classes invariant.
* **Flips [P + DV].**
  * Dual-line flip (`z₂′ = z₄′`): side A vanishes identically, at rank 0 and rank 2, because `o_H` is `2↔4`-antisymmetric (as
    OFFDIAG §1.2).
  * Light-ray-line flip (`z₁ = z₃`): side B vanishes identically at rank 0 and rank 2. Numerically `|L_flip|/|L| ≤ 2·10⁻¹³`.
  * So **the leading-Regge helicity off-diagonal of `Ĉ₅` is exactly 0.**
* **Reality.**
  * Real-analytic: `conj[R_j/R₀](ν) = [R_j/R₀](−ν)` [P, numerical to 0].
  * Not shadow-symmetric: `R_{−2}(ν)/R₂(−ν) ≠ R₀(ν)/R₀(−ν)`, whereas the scalar gives exactly 1.
  * Consequence: no ν-only phase makes `h₀` real with `h₋₂ = conj h₂`. Both `Re` and `Im` of `Ψ_{5,ν}` carry leading-Regge image.
* **Factorisation.** `R_j^{sp}/R_j^{sc} = F(ν − ij)` fails (as at `k = 4`) [P].
* **Forward normalisation.** `N̂₅ = ΣR_j` is complex (arg winds by −3.49π on `(0, 400)`).
  * No real zero on `ν > 0`. On the grid `ν ≤ 60`, `|N̂₅|/Σ|R_j| ≥ 0.057` for `ν > 0.5`, with local dips at `ν ≈ 6–12`.
  * It → 0 only at `ν → 0⁺`, where the `1/ν` poles of `R_{±2}` (residues `±192491i/2056240` after normalisation) cancel.
  * The real zeros of Stage 4's scalar-proxy `N̂₅` (4.55, 9.39 / 6.17, 11.77) were artefacts of a real-structured normalisation; in the
    proxy, `|N̂₅^{sc}|` dips to `2·10⁻⁴` at 11.77.
* **Raw `ν = 0` image [P].** Expanding `ĥ_j = R_j/N̂₅` at `ν → 0`, `ĥ_{±2} = ±c/ν + d_± + O(ν)`, and
  `P_{l+2} − P_{l−2} = 2iν∂_μP_μ|_{3/2} + O(ν³)`. This gives the exact, real image
  **`H₅(β) = κ₁∂_μP_μ(cosh β)|_{μ=3/2} + κ₂P_{3/2} + κ₃P_{−1/2}`**, `H₅(0) = 1`:

  | | `κ₁ = 2ic` | `κ₂ = d₊+d₋` | `κ₃ = ĥ₀(0)` |
  |---|---|---|---|
  | spinning `o_H q`, `j′ = 1` | **−192491/1028120 = −0.18723** | 1252043/1028120 | −223923/1028120 |
  | scalar (= proxy `b = 10`) | +5415/14984 = +0.36139 | 11609/14984 | 3375/14984 |

  `H₅/(βP_{3/2})` = `+2.01, +0.44, +0.13, +0.025, −0.006, −0.028, −0.081, −0.124, −0.145` at β = `0.5, 2, 4, 6, 7, 8, 12, 20, 30`. It
  tends to `κ₁` (as `∂_μP_μ ≈ (β + const)P_μ`). The sign change sits at β ≈ 6.9.

## 2. The 3-shift de-shift (`stage5_deshift.py` → `cache/stage5_deshift.{txt,log}`, `stage5_deshift_{C5_pp,sc2}.pkl`)

### 2.1 Formulation [D]
* **Real sum rules available.** `Re[c(ν)Ψ_{5,ν}]` and `Im[…]` are separately valid sum rules for real correlators. So the general real
  functional is `∫₀^∞Re[c(ν)Ψ_{5,ν}/N̂₅(ν)]dν`, `c` complex. Its image is `Σ_ja_j𝒫(ν−ij)` with
  `a_j = ½(cĥ_j + conj(cĥ_{−j}))` and `ĥ_j = R_j/N̂₅`.
* **Abel form.** Mehler–Dirichlet with complex index (Stage 4 §2.1) gives
  `W[c](t) = C[Re cĥ₀] + cosh2t C[Re c(ĥ₂+ĥ₋₂)] − sinh2t S[Im c(ĥ₂−ĥ₋₂)] = Φ(t)`.
  * Abel injectivity: `W = 1` on `t ≥ 0` ⇔ image `P_{−1/2}(cosh β)`, i.e. "`Ĉ₅(0)`" pure.
  * No analyticity of `c` is used.
* **Normalisation (EE way).** `Y = c·n`, `n = |a₀| + 2|a₂|` (`~1/ν` at 0); rows divided by `cosh 2t`.
* **Discretisation.** Midpoint ν grid; t grid oversampled ×2 (`t_k = kπ/(2Ndν)`), which makes the system square for the two real
  unknowns; dense SVD; min-norm.
* **Why complex `c` is needed [NE + D-sketch].**
  * Restricted gauge `c = χ` real (`Re[Ψ/N̂₅]`, forward normalisation identically 1): smallest s.v. `6·10⁻⁶ → 4·10⁻¹¹` (`dν` 0.04 → 0.02),
    `∫|χ| = 4·10⁴ → 6·10⁹`, `W(6) = 1.10 → 930`. No bounded solution.
  * `c = −iχ` (Im only) has residual 0.38. `c = χ(ν−i)/(ν+i)` and `χ(ν−2i)/(ν+2i)` are also not solvable.
  * Mechanism (sketch): the large-`t` condition is a Riemann–Hilbert problem with symbol `a₂/ā₂`. `arg a₂` runs `+π/2 → 0` (spinning)
    versus `−π/2 → 0` (scalar). The sign is that of the `ν → 0` residue, i.e. of `κ₁`, so the index has the wrong sign for a
    real-weight solution.
  * The scalar control solves in the restricted gauge with 1–2 null directions, as Stage 4.

### 2.2 Numerics [DV]
| | spinning `C5_pp` (general `c`) | scalar `sc2` (general `c`) |
|---|---|---|
| residual (`dν = .04/.02/.02`, `N = 1024/2048/3072`) | `5·10⁻¹⁵ / 7·10⁻¹⁵ / 7·10⁻¹⁵` | `3–5·10⁻¹⁵` |
| `W(t)` direct, `t = 0…8` | `1 ± 5·10⁻¹⁰` | `1` |
| `∫Re c` (forward normalisation) | 1.0000000000 | 1.0000000000 |
| `∫|c|` | **1.1495** (1.14957 / 1.14951 / 1.14950) | 1.0014 |
| light-side moments `n = 0, 1, 2` | `1, 2.4·10⁻¹³, 6.8·10⁻¹⁰` | `1, −1.6·10⁻¹², −2.5·10⁻⁹` |
| refinement `ν_max` 41 → 61; `dν` .04 → .02 | `7.7·10⁻⁶`; `4.2·10⁻³` | `2.8·10⁻³`; `6·10⁻⁴` |
| support (`|c| > 10⁻⁴ max`) | `ν ≲ 23` (tail `≲ 2·10⁻⁷` at 30–40) | `ν ≲ 23` |
| rank deficiency (s.v. `< 10⁻¹¹`) | `N − 22 / N − 25 / N − 25`, **no gap** (`1.4·10⁻¹¹ │ 5·10⁻¹²`) | `N + 1–2`, clean gap (`1.6·10⁻¹¹ │ 2·10⁻¹⁴`); the `N` are `Im`-directions with zero image |
| `Im c` share of `∫|c|` | 0.39 | 0 |

* **Independent image check (V2).** Direct `Re Σ_j∫cĥ_jP_{−1/2+i(ν−ij)}(cosh β)dν`, using mpmath complex-degree Legendre functions on the
  full grid:

  | β | 0 | 0.5 | 1 | 2 | 4 | 6 |
  |---|---|---|---|---|---|---|
  | rel. diff | `7·10⁻¹⁶` | `7·10⁻¹⁶` | `3·10⁻¹⁵` | `2·10⁻¹⁵` | `1·10⁻¹⁵` | `2.8·10⁻¹²` |

  These are against cancellations of order `e^{2β}/ν`.
* **Profile.** `c` is concentrated on `ν ∈ [0.3, 4]`: `Re c ≈ 0.13–0.34`, `Im c ≈ 0.05–0.22`. It is `O(10⁻³)` at 6–15 and `< 10⁻⁶` beyond 30.

### 2.3 Lemmas A, B, C and the light side
* **Lemma A [D, unchanged].** `Ψ_{5,ν}[𝒢] = 0` pointwise in `ν` on Lemma A's tree correlators. Hence `Re` and `Im` of it vanish, and so
  does `∫Re[cΨ/N̂₅]` (Einstein, `E×W³`, MFT).
* **Lemmas B, C [P, unchanged].** They concern `f⁻₅ = f⁻₃R` only and are `ν`- and structure-independent.
* **Light side [P given (P), (F) for the pure-moment family].**
  * `W = Φ` gives the Taylor identities `Re Σ_j∫cĥ_j(ν−ij)^{2n} = (−1)ⁿΦ^{(2n)}(0)`, checked numerically `(1, 0, 0)` to `7·10⁻¹⁰`.
  * So a polynomial pure-family light side `p(μ²)` maps to `p(0)`:
    * `−Ĉ₅^{ds}|_{R⁴} = 0`, `|_{g₅} = g₅`, `|_{|ĝ₃|²} = 0`;
    * Einstein and `E×W³`: 0;
    * **unchanged from Stage 12**, including the exact `g₅` cancellation in `∂Ĉ₄ − Ĉ₅`.
  * The R⁴ normalisation `−Ĉ₄^{ds}|_{R⁴} = 2g₄` (Stage 4) is untouched.

## 3. Validity: Regge power counting [D]
Stage 3 §1.5 counting: kernel `r^aκ(η)`, corner `r ≲ 1/η`, K1′ `|dDisc| ≤ Br^{1−J₀}` ⇒ converges iff `J₀ < a + 1 − c` for `κ ~ η^c`.

| functional | `a` | `κ` | condition | verdict | reproduces |
|---|---|---|---|---|---|
| raw `Ĉ₂` | 0 | `η^{3/2}` | `J₀ < −½` | fails | D ✓ |
| `Ĉ₃` | 0 | `η^{−½}log` | `J₀ < 3/2` | passes | D ✓ |
| `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` | 2 | `η^{−½}log`, `η^{−½}log³` | `J₀ < 7/2` | passes | Stage 4 |
| raw spinning `Ĉ₅(0)` | 2 | `κ₁η^{3/2}log η` (exact, §1.4; `κ₁ < 0`) | `J₀ < 3/2` | converges, but heavy positivity fails (§4) | Stage 3/4 (now exact) |
| **de-shifted spinning `Ĉ₅^{ds}`** | 2 | **`η^{−½}log η`** (= `P_{−1/2}`, exact image identity of the heavy-normalised total kernel) | **`J₀ < 7/2`** | **passes** (margin 2½ over KKSDZ `J₀ ≤ 1`) | new |

* CLS-D's own condition for odd `j′` (`J₀ < k − 1 = 4`) is weaker than this.
* The corner structure is Stage 4's: `Ψ₅,ν = Ψ₃,ν[R·]`, with `R → 4ρ²/(1−ρ²)²` bounded and ν-independent, so `R` factors out of the
  `ν`-integral.
* Grade [D], at the level of D's I14(b): K1′, K2, K3 and KKSDZ's loophole are inherited.
* Same caveat as D and Stage 4: the swap must use the heavy-normalised total kernel. Here `c` has effective support `ν ≲ 25` (and the
  grid stops at 41/61).

## 4. Heavy positivity with both forward pieces de-shifted (`stage5_positivity.py` → `cache/stage5_positivity.txt`) [M]

**Model.** Stage 12's leading-Regge model:
* `B₄|₀ → 2u²P₀`;
* `∂B₄|₀ → −2u²Λ²D + c₄u³P₀`;
* `B₅|₀ → σ_hu³P₀` (de-shifted: exact by §2) or `σ_hu³H₅` (raw spinning, exact §1.4);
* smeared `ψ₂`, `ψ₃` pieces by Mehler with **exact polynomial cos transforms**.

The exact transforms were needed: Stage 12's interpolated transform is off by about 1 % in `I₂₀` at β = 5, Λ = 100, which matters at
`Λβ ~ 3·10⁴`. Checks: closed form vs Gauss–Legendre `1.6·10⁻¹³`; Mehler `H₅` vs mpmath `≤ 6·10⁻¹²`.

**Grid.** β from `5/Λ` to 30 (382 values, plus a dense 400×400 scan on `b̂ ∈ [20, 200]`), `u ∈ [10⁻⁶, 1]`, physical states
`2J+1 = Λ√(2(coshβ−1)/u) ≥ 2J_min+1`. Normalisation as Stage 12/4.

| `Λ` | functional | **min, `b̂ ≥ 20` (dense)** — location | min β ≥ 2 | min β ≥ 6 | raw spinning `Ĉ₅`: first negative |
|---|---|---|---|---|---|
| 100 | δ = 0.05 | **+0.206** (`++`, J = 226, m = 21.8Λ, b̂ = 20.7) | +0.952 | +0.990 | β = 8.30, **`++`**, m ≈ Λ, **J ≈ 3.2·10³** |
| | δ = 0.2 | **+0.287** (`+−`, J = 201, m = 19.3Λ, b̂ = 20.8) | +0.939 | +0.987 | same |
| 300 | δ = 0.05 | **+0.065** (`++`, J = 263, m = 25.5Λ, b̂ = 20.7) | +0.983 | +0.996 | β = 9.30, `++`, **J ≈ 1.6·10⁴** |
| | δ = 0.2 | **+0.221** (`+−`, J = 219, m = 21.1Λ) | +0.978 | +0.995 | same |
| 1000 | δ = 0.05 | **+0.042** (`++`, J = 268, m = 25.9Λ, b̂ = 20.7) | +0.995 | +0.999 | β = 10.45, `++`, **J ≈ 9.3·10⁴** |
| | δ = 0.2 | **+0.212** (`+−`, J = 219, m = 21.1Λ) | +0.993 | +0.999 | same |

* **De-shifted: positive everywhere scanned with `b̂ ≥ 20`, up to β = 30**, i.e. `J ≈ 2·10⁸ / 5·10⁸ / 2·10⁹` at `m ≈ Λ`. That is far
  beyond the old proxy failure points (`J ≈ 1.5·10³ / 8.6·10³ / 5.7·10⁴`) and the new raw-spinning ones. At the old failure points the
  value is `≥ +0.99`.
* The `b̂ ≥ 20` minima sit at the model's validity edge (`b̂ ≈ 20.7`), at `m ≈ 20–26Λ`, `J ≈ 200–270`. Two notes:
  * These are flat-regime states where the forward pieces enter with `u² ≲ 3·10⁻³` and `Ĉ₅` with `u³ ≲ 10⁻⁴`. **They are not caused by
    `Ĉ₅` and do not move with the de-shift.**
  * They are **lower than Stage 4's quoted `+0.80 / +0.74 / +0.37`.** That scan used β-spacing 0.125 and `u ≥ 10⁻⁴`, so it skipped
    `b̂ ∈ (20, ~55)`; Stage 12 itself had reported `+0.27`.
  * The δ = 0.05 value at Λ = 1000 (+0.042, in this model's normalisation `Nrm`) approaches the flat LP margin. For δ = 0.2 the margin is
    ≥ 0.21.
* **Raw spinning `Ĉ₅`:**
  * The failure is now in **`++`**, because `κ₁ < 0`; the proxy failed in `+−`.
  * It occurs at slightly larger `J` than the proxy predicted (`H₅` turns negative only beyond β ≈ 6.9).
  * `FF/Nrm → −1` at β ≈ 27–30.
  * No LP re-weighting can fix it: the term grows like `e^{3β/2}β`. **The de-shift is necessary and sufficient at leading Regge** [M].
* **`b̂ < 20`** (outside the model) is Stage 3's near-gap regime. The hybrid-exact pure-moment model there is unchanged, because it already
  used pure moments (Stage 4 §4.1) and the spinning de-shift realises the pure moment. Minima stay `+0.0478 / +0.0499 / +0.0500`
  (δ = 0.05) and `+0.1985 / +0.2003 / +0.2001` (δ = 0.2) [D identification + M].
  * The leading-Regge model alone gives `−0.14` (δ = 0.05) and `−0.09` (δ = 0.2) at `b̂ ≈ 7–8`, `J = 4`, `m ≈ 1.1–1.2Λ`. This is
    expected, since that is outside its validity.
* **Constants.** The light sides map exactly (§2.3), so:
  * δ = 0.05: `A = 11.448`, `B_eff = −12.80 / −12.94 / −12.97`, `κ = 0.0686`;
  * δ = 0.2: `A = 17.241`, `B_eff = −19.04 / −19.22 / −19.27`, `κ = 0.110`;
  * **all unchanged**. ABJM sits at 6.8 / 5.0 / 3.9 % (δ = 0.05) and 4.5 / 3.3 / 2.6 % (δ = 0.2) of the bound, also unchanged.

## 5. Gauge (`stage5_gauge.py` → `cache/stage5_gauge.txt`) [NE]
Weighted min-norm gauges on `dν = 0.04`, `N = 1024`. All have residual `≤ 4·10⁻¹⁴`, `W(0…6) = 1`, and moments `(1, ≤ 2·10⁻¹², ≤ 6·10⁻⁹)`:

| gauge | `∫|c|` | `Im` share | support (`10⁻⁴`) | `∫|c−c₀|/∫|c₀|` |
|---|---|---|---|---|
| G0 min-norm (used) | 1.150 | 0.39 | ν ≤ 23 | 0 |
| G1 penalise `(1+(ν/5)²)` | 1.112 | 0.30 | ν ≤ 19 | 0.17 |
| G2 `Im` ×10 (as real as possible) | 1.139 | 0.26 | ν ≤ 22 | 0.29 |
| G3 `Re` ×10 | 3.401 | 0.93 | ν ≤ 41 | 2.57 |

* The leading-Regge image, light side and validity are gauge-independent.
* Beyond leading Regge (`O(1/m)`, `O(J²/m)`, `x₀ ≠ 0` dual-structure terms), the forward `Ĉ₅` contribution near the gap depends on the
  gauge at `O(20–30 %)` of itself between reasonable gauges.
* Unlike the `Ĉ₂′` (H0 §R, 3-dim null space) and `Ĉ₄` (3 + 1) cases, the null space here is infinite-dimensional, because `Re` and `Im`
  both carry image. A canonical definition therefore needs an extra principle, e.g. a CHMRS-class analyticity in `ν` or a shadow-covariant
  pairing of `Ψ_ν` and `Ψ_{−ν}`. **Not fixed [I].**
* Since `Ĉ₅` enters the near-gap states with weight `a·u³` (`u ≈ 0.7` at the binding `+−`, `J = 5` state), this gauge ambiguity is part
  of the open graviton-specific near-gap correction (Stage 4 §5).

## 6. Remaining open items [I]
1. A canonical gauge for the 3-shift (complex, infinite-dimensional null space) and for the 4-shift weights. The heavy side beyond
   leading Regge is gauge-dependent (§5).
2. The graviton-specific near-gap (`ν = 0`, `J = 4–8`) corrections of both de-shifted forward pieces, which decide δ = 0.05 vs 0.2: not
   computed (Stage 4 §5's crude estimate stands). δ = 0.2 covers all estimates.
3. The leading-Regge model at `b̂ ≈ 20` (minima +0.042–0.065 for δ = 0.05 at Λ ≥ 300) has no finite-`b̂` correction. The matching to
   Stage 3's near-gap model across `b̂ ∈ [10, 30]` at large `m/Λ` is not done.
4. The validity argument inherits D's K1′, K2, K3 and KKSDZ's loophole (counting at the level of the heavy-normalised kernel).
5. The index explanation of the restricted gauge's failure is a sketch.
6. Inherited from Stages 0–4: `b₀`, `b₅`, `c₁`, spinning (F) and (P), A1, A2, microcausality of tree correlators, the odd-`k` CHMRS
   continuation used in Stage 3's near-gap model.

## 7. Graded claims
| claim | grade |
|---|---|
| exact `R_j` (`j = 0, ±2`) of raw spinning `Ψ_{5,ν}[o_H q]`; only `j = 0, ±2`; rank 0 = q1 exactly; rank-2 scalar = Stage-4 proxy (`b = 10`) | [P] (sympy) + [DV] |
| independent finite-difference check (q1 engine) | [DV] (`≤ 1.3·10⁻⁶`) |
| image carries `σ_h` exactly; total-flip invariant; both line flips vanish identically (off-diagonal 0); `n`-odd parts 0 | [P] + [DV] |
| real-analytic, not shadow-symmetric; `N̂₅` complex, no real zero on `ν > 0` | [P] / [NE] (grid `ν ≤ 60`) |
| raw `ν = 0` image `H₅ = κ₁∂P_{3/2} + κ₂P_{3/2} + κ₃P_{−1/2}`, `κ₁ = −192491/1028120` (opposite to scalar proxy) | [P] |
| restricted real-weight de-shift not solvable; mechanism | [NE]; [D-sketch] |
| general-`c` de-shift: residual `10⁻¹⁵`, image check `≤ 3·10⁻¹²`, moments `≤ 7·10⁻¹⁰`, `∫|c| = 1.15`, converged; null space infinite-dimensional, zero image | [DV] |
| Lemma A, B, C carry; light side `R⁴: 0`, `g₅: g₅`, `|ĝ₃|²: 0`, unchanged | [D] / [P] / [P given (P), (F)] |
| de-shifted `Ĉ₅` needs `J₀ < 7/2` | [D] (D-level counting) |
| leading-Regge positivity with both pieces de-shifted, `b̂ ≥ 20`, β ≤ 30: minima above; raw spinning fails in `++` at `J ≈ 3.2·10³ / 1.6·10⁴ / 9.3·10⁴` | [M] |
| near-gap minima, `A`, `B`, `κ`, ABJM placement unchanged | [D] (identification) + [M] |
| gauge dependence 17–29 % between natural gauges | [NE] |

## 8. Files (`ads4/code/resultq/`)
| file | content | run time |
|---|---|---|
| `stage5_hj.py` | rank-0/rank-2 exact `R_j` for `o_H q` / `j′ = 1` (4 helicity classes, 3 flip configurations, scalar control), checks → `cache/stage5_hj.{txt,pkl}` | 30 s |
| `stage5_check.py` | finite-difference q1-engine check, 5-harmonic fit, flip check → `cache/stage5_check.txt` | 30 s |
| `stage5_deshift.py` | restricted vs general-gauge 3-shift solver (`Deshift3S`), refinement, V2 (mpmath), V3, null space → `cache/stage5_deshift.{txt,log}`, `stage5_deshift_{C5_pp,sc2}.pkl` | 4 min |
| `stage5_positivity.py` | leading-Regge model, exact cos transforms, β ≤ 30, pure vs raw spinning `Ĉ₅` → `cache/stage5_positivity.txt` (the dense `b̂ ∈ [20, 200]` scan was run inline with `scan()` of this module, 400 β × 400 u per case) | 1 min (+ 5 min dense) |
| `stage5_gauge.py` | gauge comparison; writes **`functionals_q5_deshift.json`** (exact `R_j`, `R₀^{(3)}`, `κ₁…κ₃`, min-norm `c(ν)` (Re, Im), `dν = 0.02`, `N = 2048`) | 15 s |
