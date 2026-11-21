# Result Q, Stage 4: the 4-shift de-shifted forward `Ĉ₄`, its validity, and what it changes on the heavy side

*21 Nov 2026, Stage-4 agent. Unrefereed, not committed by this agent. No existing document or code file was edited. New code:
`ads4/code/resultq/stage4_hj.py`, `stage4_deshift.py`, `stage4_c5.py`, `stage4_save.py`; kept functional file
`ads4/code/resultq/functionals_q4_deshift.json`; outputs in `resultq/cache/stage4_*` (git-ignored). Read-only imports: `h7i14/v3_regge_kernel.py`,
`v4_hclosed.py`, `spinfunc/u2_hclosed_struct.py`, `s1_offdiag_kernel.py`, `resultq/stage12_forward_regge.py`. Compute: about 25 CPU-minutes.
(Note: `stage4_hj.py` appears in commit `628f94c`, which this agent did not make.)*

Grades: **[P]** proved (incl. exact computer algebra), **[D]** derived with the stated assumption, **[DV]** derived and checked numerically,
**[M]** model / numerical, **[I]** open. Flat = D = 4 amplitudes; AdS = the CFT₃ functionals. `𝒫(μ) := P_{−1/2+iμ}(cosh β)`.

---

## 0. Answer in brief

1. **Exact raw kernel [P].** The leading-Regge heavy image of the raw `Ψ_{4,ν}[EE]` is computed exactly with rank-3 equivariant tensors:
   `Σ_{j=±1,±3} h_j(ν)𝒫(ν − ij)`, with `h_j = N_com R_j` and `R_j` rational (§1).
   * Checks: the same code at rank 1 reproduces v4_hclosed (`q²`) and u2 (`EE`) exactly; at rank 3 the scalar reproduces Stage 12's
     proxy `a⁽³⁾_j(l)E₉(ν−ij)` to `10⁻¹⁶`.
   * Helicity-blind and flip-free: side `A` is 16× the scalar's, so flips carry `(ε·ε)² = 0`.
   * **New:** Result D's `k = 2` rule ("EE = scalar × c(ν′) per harmonic") **fails at `k = 4`**.
2. **Real zeros [P].** No single `h_j` vanishes on `ν > 0`. The forward normalisation `N̂₄ = Σ_j h_j` has **four real zeros**:
   `ν = 1.068482, 4.545932, 8.169433, 13.094015` (the scalar skeleton has three: 3.718, 7.536, 12.556).
   * They are handled as in EE: χ-normalisation `Y = χΣ|h_j|`, a bounded operator, and `χN̂` vanishes linearly at each zero.
3. **The de-shift [D + DV].** The 4-shift Abel equation
   `W[χ](t) := Σ_{j=1,3}[cosh jt C[χA_j] − sinh jt S[χB_j]] = Φ(t)`, `t ≥ 0`,
   is solved for `Φ = 1` (pure image `P_{−1/2}(cosh β)`, i.e. `Ĉ₄(0)`) and for `Φ = −t²/2` (pure image `−D(β)`, i.e. `∂_{ν²}Ĉ₄(0)`).
   * It realises the exact forward (`ν = 0`) pure moment. No `ν`-smearing is needed: the δ-target comes from the strip-edge poles, like
     D's `ν = 0` term.
   * Residual `5·10⁻¹⁵`.
   * Independent image check (direct `ν`-quadrature, complex-degree Legendre functions from mpmath): `3·10⁻¹²` (β = 4) and `2·10⁻¹⁰` (β = 6).
   * Null space: **3** exact directions (EE and scalar alike), plus one ill-conditioned EE direction (s.v. `10⁻⁴ → 10⁻⁷` under refinement).
     All null directions have zero leading-Regge image.
   * Weight size: `‖χN̂‖₁ = 1.13` (value), `0.071` (derivative). Small amplification.
4. **Validity [D].** Under the counting that reproduces D's verdicts (raw `Ĉ₂` fails, `Ĉ₃` passes), the de-shifted kernels behave as
   `κ = η^{−½}log η` (value) and `η^{−½}log³η` (derivative). The large-η corner then needs **`J₀ < 7/2`**: margin 2½ over KKSDZ's `J₀ ≤ 1`.
   Lemma B (`f₄ = f₂R`, `R` bounded in the corner) lets `R` factor out of the de-shift.
5. **Lemmas A, B, C carry over unchanged [D].** A: pointwise zero in `ν`. B, C: `ν`-independent.
   * **R⁴ normalisation unchanged, exactly:** `−Ĉ₄^{ds}|_{R⁴} = 2g₄·W(0) = 2g₄`, numerically `∫χN̂ = 1.0000000000`.
   * Every polynomial light side maps exactly: `Σ_j∫χh_j(ν−ij)^{2n} = (−1)ⁿΦ^{(2n)}(0)` (checked: `1, 0, 0` and `0, 1, 0`). So Stage 12's
     table (`2g₄ + b₀G|ĝ₃|²`, `4πG|ĝ₃|² + g₅`) holds for the de-shifted pieces [P given (P), (F) for the pure-moment family].
6. **Heavy side.**
   * Stage 3's near-gap model was already the pure-moment (de-shifted) one: CHMRS `C_{k,ν}`, `ex₄ ≈ 2P₀` to `10⁻³`. So **the images,
     minima and constants are unchanged**:
     * δ = 0.05: min `+0.0478 / +0.0499 / +0.0500`, `A = 11.448`, `B_eff = −12.80 / −12.94 / −12.97`, `κ = 0.0686`;
     * δ = 0.2: min `+0.1985 / +0.2003 / +0.2001`, `A = 17.241`, `B_eff = −19.04 / −19.22 / −19.27`, `κ = 0.110`;
     * no re-optimisation needed [M].
   * **New caveat on `Ĉ₅`.** The **raw** `Ĉ₅(0)` that Stage 3 kept has a leading-Regge image `∝ βP_{3/2}` (η^{3/2} log η) with sign
     `σ_h`. Stage 3/12's model silently used the pure moment.
     * With raw `Ĉ₅` both LP functionals go negative in the `+−` channel at `J ≈ 1.5·10³ / 8.6·10³ / 5.7·10⁴` (`m ≈ Λ`, `Λ = 100 / 300 / 1000`) [M, proxy].
     * No LP re-weighting can fix this. **A 3-shift de-shifted `Ĉ₅` is required.** It is solvable in the scalar proxy (§3), with real zeros
       of its forward normalisation, null dimension 1–2, residual `4·10⁻¹⁵`.
     * The spinning (`j′ = 1`) coefficients are [I].
7. **Task 3 (graviton-specific `ν = 0` corrections at J = 4–8) is not computed [I].** An extrapolation of D's measured EE near-gap
   corrections gives `≈ 0.007–0.10 N` at `Δ_g = 100` and `≲ 0.03` at `300` [M, crude]. So δ = 0.05 is undecided at 100 and plausible at
   ≥ 300; δ = 0.2 covers all estimates.
8. **ABJM** (`2ζ(3) = 2.404`): **6.8 / 5.0 / 3.9 %** of the δ = 0.05 bound, **4.5 / 3.3 / 2.6 %** of the δ = 0.2 bound, at
   `Δ_g = 100 / 300 / 1000`. Unchanged; inside, far from both boundaries [D given Q; bound M].

---

## 1. The exact raw kernel (`stage4_hj.py` → `cache/stage4_hj.txt`)

### 1.1 Method [P]
* By Lemma B, `f⁺₄ = f⁺₂R`, and `R → 16η²r²` at leading Regge. So the `k = 4` kernel is `16r²η³K₀ = −16(x·y)³/(|x||y|)K₀`.
  * Per side this is `|x|^{b−3}x^μx^νx^ρ` with the **same `b = 6`** as `k = 2`.
  * Hence `L₄ = −16⟨∂³Â, ∂³B̂⟩`, with `Â`, `B̂` the v3 plane-wave transforms (same monomials, same `K(c,λ)` ratios as v4_hclosed).
* **Rank-3 equivariance.**
  * The derivative patterns at one celestial point come from a truncated Taylor expansion of `S^{−λ}R^κ`.
  * Each side's symmetric 3-tensor is decomposed on `{p³, p²ẑ, pẑ², ẑ³, ηp, ηẑ}` plus 4 `n`-odd tensors. The `n`-odd ones vanish (asserted).
  * Contraction uses Gram algebra in `{p, p̄, ẑ_A}`. The `z₀` integral maps `C^k r^s` to `2πC^kP_{l+s}` (v4's `I(a)`), and the Legendre
    recurrence then collects the `P_{l+j}`.

### 1.2 Checks
* **(C1)** Rank 1 vs v4_hclosed / u2: differences `0` for scalar, `q²` and `EE` [P].
* **(C2)** Rank-3 scalar `R_j/R_{−1}` vs Stage 12's `a⁽³⁾_j(l)E₉(ν−ij)`: `≤ 3·10⁻¹⁶` at `ν = 0.3, 1.7, 5.2` [DV]. Stage 12's scalar proxy was
  therefore exact.
* **(H)** `R_j^{EE}(h = −1) = R_j^{EE}(h = +1)`: helicity-blind [P].
* Side `A` of `EE` is exactly 16× the scalar's on all six rank-3 tensors (x1 (E) extended). Helicity flips therefore vanish identically
  [P].

### 1.3 Result [P]
```
R_{−3} = −(2ν−i)(2ν+i)(2ν+5i)(2ν+9i)(2ν+11i)²(2ν+13i)(2ν+15i)(2ν+19i)(4ν²+12iν−1) / [8ν(ν+i)(ν+2i)(2ν−7i)(2ν−3i)]
R_{−1} = −3(2ν−11i)(2ν−5i)(2ν−3i)(2ν−i)(2ν+i)(2ν+5i)²(2ν+9i)(2ν+11i)(2ν+15i) / [8ν(ν−i)(ν+2i)(2ν−7i)]
R_{+j}(ν) = R_{−j}(−ν)
```
* `A_j := R_j + R_{−j}` and `B_j := −i(R_j − R_{−j})` are real on `ν > 0` (to `5·10⁻¹³`).
* The `1/ν` poles cancel in `N̂₄ = A₁ + A₃`, which is finite at 0.
* **The `k = 2` factorisation fails:** `R_j^{EE}/R_j^{sc}` is not `F(ν − ij)` for a single `F`. The differences are explicit non-zero rationals
  for `j ≠ −1` (printed in the cache file).

### 1.4 Zeros of the forward normalisation [P]
* Real zeros of `N̂₄^{EE}` on `ν > 0`: **1.068482, 4.545932, 8.169433, 13.094015**.
* Zeros in `|Im ν| ≤ 3.5` off the real axis: `±i/2, ±i, ±2i`.
* No `R_j` has a real zero on `ν > 0`.
* Scalar skeleton: real zeros 3.718, 7.536, 12.556.
* `N̂₄^{EE}(0⁺) < 0` in this convention, `N̂₄^{sc}(0⁺) > 0`. The overall sign is fixed by `W(0) = Φ(0)` (§2), so it is immaterial.

## 2. The 4-shift de-shift (`stage4_deshift.py` → `cache/stage4_deshift.txt`, `stage4_deshift_v2full.txt`)

### 2.1 Formulation [D]
* Mehler–Dirichlet holds for complex index: `cos((ν−ij)t) = cos νt cosh jt + i sin νt sinh jt`.
* So `∫₀^∞χΨ_{4,ν}dν` has image `(√2/π)∫₀^β W(t)(cosh β − cosh t)^{−½}dt`, with
  `W[χ](t) := Σ_{j=1,3}[cosh jt C[χA_j] − sinh jt S[χB_j]]`.
* Since the Abel transform is injective, `W = Φ := C[ψ]` on `t ≥ 0` is equivalent to the image being `∫ψ𝒫`. No analyticity of χ is used.
* **Forward targets:**
  * `ψ = δ(ν − 0⁺)` ⇔ `Φ = 1` (image `P_{−½}`);
  * `ψ = ∂_{ν²}δ` ⇔ `Φ = −t²/2` (image `−D(β)`).
* χ-normalisation (as u3): `Y = χn`, `n = Σ|R_j|`, rows divided by `cosh 3t`. The real zeros of `N̂₄` are never divided by.

### 2.2 Numerics [DV]
Midpoint grid with the paired DCT `t`-grid; dense SVD; min-norm solution.

| | EE value | EE deriv | scalar value | scalar deriv |
|---|---|---|---|---|
| residual (`dν = .02`, `N = 2048`) | `4.8·10⁻¹⁵` | `4.8·10⁻¹⁵` | `3.8·10⁻¹⁵` | `3.6·10⁻¹⁵` |
| `∫χN̂ = W(0)` (target) | `1.0000000000` (1) | `0.0000000000` (0) | `1.0000000000` | `0.0000000000` |
| `W(t)` by direct quadrature, `t ≤ 8` | `1 ± 6·10⁻⁶` | `−t²/2 (1 ± 2·10⁻⁸)` | `1` | `−t²/2` |
| `‖χN̂‖₁` | 1.130 | 0.071 | 1.008 | 0.082 |
| support (`|Y| > 10⁻⁶ max`) | `ν ≤ 27` | `ν ≤ 26` | `ν ≤ 27` | `ν ≤ 27` |
| null directions (s.v. `< 10⁻⁹`) | 3 (+1 at `1.4·10⁻⁷`) | same | 3 | same |

* **Independent image check.** Direct `Σ_j∫χR_j P_{−½+i(ν−ij)}(cosh β)dν` with mpmath Legendre functions on the full grid:

  | β | value | derivative |
  |---|---|---|
  | 4 | `2.9·10⁻¹²` | `5.5·10⁻¹⁰` |
  | 6 | `2.0·10⁻¹⁰` | `3.9·10⁻¹⁰` |

  This is against cancellations of order `e^{3β}`. (A stride-2 sub-sampled quadrature degrades to `10⁻²–1` at β ≥ 4; that is quadrature
  error, not a failure of the identity.)
* **Refinement.** `dν 0.04 → 0.02`: `Δχ/χ ≤ 6·10⁻⁴`. `ν_max 41 → 61`: `1.4 %` / `0.35 %`. This is a change of the min-norm gauge only:
  `W` is unchanged.
* **Profile.** `χN̂` is concentrated at `ν ∈ [0.5, 3]`: `+0.39` at 1.5, `+0.70` at 2. It vanishes linearly through the real zeros
  (`9.5·10⁻⁴` at 1.07, `3·10⁻⁶` at 4.55).
* **Gauge.** The leading-Regge image is gauge-independent: null vectors have `W ≡ 0`. A canonical (CHMRS-class) gauge was **not**
  attempted [I].

### 2.3 Validity: Regge power counting [D]
Stage 3 §1.5's counting is: kernel `r^aκ(η)`, corner `r ≲ 1/η`, K1′ `|dDisc| ≤ Br^{1−J₀}`; it converges iff `J₀ < a + 1 − c` for `κ ~ η^c`.

| functional | `a` | `κ` | condition | verdict | D check |
|---|---|---|---|---|---|
| raw `Ĉ₂` | 0 | `η^{3/2}` | `J₀ < −½` | fails | ✓ (D) |
| `Ĉ₃` | 0 | `η^{−½}log` | `J₀ < 3/2` | passes | ✓ (D) |
| raw `Ĉ₄(0)`, `∂Ĉ₄(0)` | 2 | `η^{5/2}` | `J₀ < ½` | fails | Stage 3 |
| **de-shifted `Ĉ₄^{ds}` (value)** | 2 | **`η^{−½}log η`** (= `P_{−½}`) | **`J₀ < 7/2`** | **passes** | new |
| **de-shifted `∂Ĉ₄^{ds}`** | 2 | **`η^{−½}log³η`** (= `D(β)`) | **`J₀ < 7/2`** | **passes** | new |
| raw `Ĉ₅(0)` | 2 | `η^{3/2}log η` (§3) | `J₀ < 3/2` | passes (convergence only) | Stage 3 |

* The `κ` of the de-shifted pieces is now an exact identity of the heavy-normalised total kernel (§2.2), not an estimate.
* **Corner structure.** By Lemma B, `Ψ₄,ν = Ψ₂,ν[R·]` with `ν`-independent `R`. In the large-η corner `R → 4ρ²/(1−ρ²)²` is bounded, so it
  factors out of the `ν`-integral; the de-shift acts on the `ν`-dependent `k = 2`-type kernel only.
* **Grade [D]**, at the level of D's I14(b) (Proposition K for the heavy-normalised total kernel; K1′, K2, K3 inherited).
* **Caveat (as for D).** In the un-normalised position kernel, the `ν`-integral is not absolutely convergent: `1/E_b(ν) ~ e^{π|ν|}` against
  a weight decaying like `e^{−ν}`. The swap must use the heavy-normalised total kernel, exactly as D's `ψ₂`-smeared `Ĉ′_EE` (supported to
  `ν = Λ`) already must.

### 2.4 Lemmas A, B, C and the light side
* **Lemma A [D, unchanged].** `Ψ_{4,ν}[𝒢] = 0` for every `ν` on Lemma A's tree correlators, so `∫χΨ_{4,ν}[𝒢]dν = 0`. No swap is needed.
* **Lemmas B, C [P, unchanged].** They concern `f₄` only and are `ν`-independent.
* **Light side [P given (P), (F) for the pure-moment family].**
  * Taylor coefficients of `W = Φ` at `t = 0` give `Σ_j∫χR_j(ν−ij)^{2n}dν = (−1)ⁿΦ^{(2n)}(0)`.
    * Value: `(1, 0, 0)`; derivative: `(0, 1, 0)` for `n = 0, 1, 2` (numerics: `1, 0, 4·10⁻⁹` and `0, 1, 0`).
    * So a polynomial light side `p(μ²)` of the pure moments maps to `p(0)` (value) and `p′(0)` (derivative).
  * Hence:
    * `−Ĉ₄^{ds}|_{R⁴} = 2g₄` exactly (**unchanged**);
    * `−∂Ĉ₄^{ds}|_{R⁴} = 0`;
    * `|ĝ₃|²`: `4πG|ĝ₃|²` in the derivative, `b₀G|ĝ₃|²` in the value;
    * `g₅` entries as in Stage 12.
  * Note: the naive `∫χN̂ν²` (`= 5.57` / `0.096`) is **not** the relevant moment. The forward-normalised raw `Ĉ₄(ν)` does not have a
    polynomial light side. (P) refers to the pure-moment family, which the de-shift realises.

## 3. The `k = 5` piece (`stage4_c5.py` → `cache/stage4_c5.txt`)

**(a) Raw image [D + DV, scalar proxy].**
* The `k = 5` kernel is `−16η²r² ×` (`k = 3` kernel), and the `k = 3` image is a pure `P_l`. So the raw image is `Σ_{j=0,±2}h_j𝒫(ν−ij)`.
* At `ν → 0` the `1/(2l+1)` poles cancel into a `∂_μ` term. Proxy `h_j = a⁽²⁾_j E_b(ν−ij)`:
  `H₅(β)/(βP_{3/2}) = 2.08, 0.72, 0.53, 0.44, 0.41` at `β = 0.5, 2, 4, 8, 12`, for both `b = 8` and `b = 10`.
* So the growth is `η^{3/2}log η`, insensitive to the proxy's radial weight.

**(b) Consequence [M].** In Stage 12's leading-Regge model (`b̂ ≥ 20`, β ≤ 20):

| `Λ` | functional | min with de-shifted (pure) `Ĉ₅` | raw `Ĉ₅`: first negative |
|---|---|---|---|
| 100 | δ = 0.05 / δ = 0.2 | +0.80 / +0.77 | β = 6.8, `+−`, `m = Λ`, `J ≈ 1.5·10³` |
| 300 | | +0.74 / +0.68 | β = 8.1, `J ≈ 8.6·10³` |
| 1000 | | +0.37 / +0.40 | β = 9.5, `J ≈ 5.7·10⁴` |

* The failing term is `aσ_hu³H₅`, growing like `e^{3β/2}`. It beats `2aΛ²u²D ~ β³e^{−β/2}` once `J ≳ m^{3/2}`.
* No choice of the LP weights can compensate it, so **positivity is lost with raw `Ĉ₅` and cannot be restored by re-optimisation.**
* Stage 3's verdict "raw `Ĉ₅` fine" holds for convergence only.

**(c) 3-shift de-shift [M (proxy) / I (spinning)].**
* `W = C[χh₀] + cosh 2t C[χA₂] − sinh 2t S[χB₂] = 1` is solvable.

  | `b` | residual | `W(t ≤ 6)` | null dim | `‖χN̂‖₁` | real zeros of `N̂₅` |
  |---|---|---|---|---|---|
  | 8 | `5·10⁻¹⁵` | 1 | 2 | 1.009 | 4.55, 9.39 |
  | 10 | `4·10⁻¹⁵` | 1 | 1 | 1.001 | 6.17, 11.77 |

* With it the `k = 5` image is `σ_hu³P₀`: exactly Stage 3/12's model.
* The exact spinning `j′ = 1` (`o_H·q`) coefficients need a `j′ = 1` version of v3's side-B monomials, which does not exist here [I]. By the
  `k = 4` experience (§1.3), they will differ from the proxy by `O(1)` rationals and may move the real zeros, but the same Abel form handles
  that.

## 4. Heavy images, constants, ABJM (task 2)

### 4.1 Why the Stage-3 numbers stand
* Stage 3's forward images were CHMRS's scalar `C_{k,ν}` (Mellin-defined). Its leading-Regge image is the pure `P_{−½+iν}`, and Stage 3
  measured `|ex₄(0)/(2P₀) − 1| ≤ 10⁻³`.
* So Stage 3 evaluated the **de-shifted** functionals all along, never the raw ones (whose image grows like `P_{5/2}`).
* Model consistency, in the scalar skeleton: the EE de-shift maps `Σ_jh_j^{EE}(ν)C^{sc}_{4,ν−ij}` to `C^{sc}_{4,0}`. This is the same
  `ex₄(0)` and `∂_{ν²}ex₄(0)` [D].

### 4.2 Table
Same functionals, same model; `B_eff = 2L_AdS − A log(2Λ/3)` with the exact pole.

| functional | `Δ_g` | min `FF/N` (location) | `A` | `B_eff` | `κ` | `2L_AdS` | ABJM |
|---|---|---|---|---|---|---|---|
| `S_min+B₃`, δ = 0.05 | 100 | +0.0478 (`+−` J = 5, 1.2Λ) | 11.448 | −12.804 | 0.0686 | 35.10 | 6.8 % |
| | 300 | +0.0499 (`+−` J = 5, 1.2Λ) | | −12.938 | | 47.66 | 5.0 % |
| | 1000 | +0.0500 (`++` J = 0, gap) | | −12.970 | | 61.45 | 3.9 % |
| re-opt. δ = 0.2 | 100 / 300 / 1000 | +0.1985 / +0.2003 / +0.2001 | 17.241 | −19.038 / −19.220 / −19.265 | 0.1099 | 53.11 / 72.04 / 92.81 | 4.5 / 3.3 / 2.6 % |

* `b̂ ≥ 20`: Stage 12's leading-Regge model with both forward pieces de-shifted gives min `+0.80 / +0.74 / +0.37` (δ = 0.05) and
  `+0.77 / +0.68 / +0.40` (δ = 0.2) (§3 (b), "pure").
* `κ`: Stage-0 value; the `∂Ĉ₄` coefficient `4πG|ĝ₃|²` is exact by §2.4.
* Lower bound: slope `≥ 0.2519` unchanged [M], with `b₀` [I].
* **No positivity lost** (given a de-shifted `Ĉ₅`), so **no re-optimisation is needed** [M].

### 4.3 ABJM
`g₄Λ⁶/(8πG) = 2ζ(3)`, `ĝ₃ = 0`. It sits inside, far from both boundaries, at the fractions in the table [D given Q; bound M]. The test
remains weak (log-dominated margin).

## 5. Task 3: graviton-specific corrections at `ν = 0` near the binding states [I]
* **Not computed.**
  * The de-shift removes the entire leading-Regge EE-vs-scalar difference exactly (§2).
  * What remains is beyond leading Regge: `O(1/m)`, `O(J²/m)`, and the `x₀ ≠ 0` dual-structure terms. No engine evaluates these at small `ν`.
  * The de-shift does not make this harder: `‖χN̂‖₁ = 1.13`. Raw-functional errors smooth in `ν` are mapped to their `ν⁰` (value) and
    `ν²` (derivative) Taylor coefficients.
* **Crude extrapolation [M].** EE_FUNCTIONAL §5 measured diagonal corrections of about `0.5 N_s` (J = 18, 20, `m ≈ 1.06Λ`, `Δ_g = 300`)
  and `0.02–0.13 N_s` at `Δ_g = 1000`. Scaling these to the binding state (`+−`, J = 5, `m ≈ 1.2Λ`), linearly in `J²/m` or steeper (the
  300 → 1000 data suggest `∝ m^{−2.7}`), gives:

  | `Δ_g` | estimate | δ = 0.05 margin |
  |---|---|---|
  | 100 | `0.007–0.10 N` | 0.048 |
  | 300 | `4·10⁻⁴–0.03` | 0.05 |
  | 1000 | `≤ 0.01` | 0.05 |

* The estimate mixes normalisations (`N_s` vs `N`) and pieces (smeared `k = 2` vs forward `k = 4`).
* **Verdict:**
  * δ = 0.05 is not established at `Δ_g = 100`, and plausible at ≥ 300;
  * δ ≈ 0.2 covers every estimate;
  * the decision itself remains [I].

## 6. Remaining open items [I]
1. The spinning (`j′ = 1`) 3-shift de-shift of `Ĉ₅`; here only the proxy (§3 (c)). It is **required**: the raw `Ĉ₅` breaks positivity at large `J`.
2. Graviton-specific beyond-leading-Regge corrections of both forward pieces near the gap (§5). These decide δ = 0.05 vs 0.2.
3. Canonical (CHMRS-class) gauge for the 4-shift weight (3 null directions + 1 ill-conditioned); the heavy side beyond leading Regge is
   gauge-dependent.
4. The validity argument inherits D's K1′, K2, K3 and KKSDZ's loophole. The counting is at the level of the heavy-normalised kernel.
5. Inherited from Stages 0–3: `b₀`, `b₅`, `c₁`, spinning (F) and (P), A1, A2, microcausality of tree correlators.

## 7. Graded claims

| claim | grade |
|---|---|
| exact `R_j` (`j = ±1, ±3`) of `Ψ_{4,ν}[EE]`; rank-1 reproduction of v4/u2; scalar = Stage-12 proxy | [P] (sympy) + [DV] |
| helicity-blind; flips vanish; side A = 16 × scalar at rank 3; `k = 2` per-harmonic factorisation fails at `k = 4` | [P] |
| real zeros of `N̂₄^{EE}`: 1.0685, 4.5459, 8.1694, 13.0940; `h_j` zero-free on `ν > 0` | [P] |
| Abel form (Mehler, complex index) ⇔ pure-moment image | [D] |
| solution for `Φ = 1`, `−t²/2`: residual `10⁻¹⁵`, direct image check `≤ 4·10⁻¹⁰`, 3 null directions | [DV] / [NE] |
| `−Ĉ₄^{ds}|_{R⁴} = 2g₄`; polynomial light sides map to `p(0)`, `p′(0)` | [P given (P), (F)] + [DV] |
| Lemma A, B, C unchanged | [D] / [P] / [P] |
| de-shifted `Ĉ₄`, `∂Ĉ₄` need `J₀ < 7/2` | [D] (D-level counting; heavy-normalised kernel) |
| Stage-3 heavy images = de-shifted images; minima, `A`, `B`, `κ` unchanged | [D] (identification) + [M] (numbers) |
| raw `Ĉ₅(0)` image `~ βP_{3/2}`, positivity lost at `J ≳ m^{3/2}` | [D] + [M] (proxy) |
| 3-shift de-shift of `Ĉ₅` solvable (proxy) | [M]; spinning [I] |
| ABJM at 6.8 / 5.0 / 3.9 % (δ = 0.05), 4.5 / 3.3 / 2.6 % (δ = 0.2) | [D given Q]; bound [M] |
| graviton-specific near-gap correction estimate | [M] (extrapolation); decision [I] |

## 8. Files (`ads4/code/resultq/`)

| file | content | run time |
|---|---|---|
| `stage4_hj.py` | rank-1/rank-3 exact `R_j` (scalar, `q²`, `EE`), checks C1, C2, H, factorisation test → `cache/stage4_hj.{txt,pkl}` | 40 s |
| `stage4_deshift.py` | 4-shift Abel solver (`Deshift4`), targets value/derivative, SVD null space, refinement, direct image check → `cache/stage4_deshift.{txt,log}`, `stage4_deshift_{EE3,sc3}.pkl` | 70 s; the full-grid check (`image_direct(S, chi, beta, stride=1)` at β = 4, 6) was run inline, output `cache/stage4_deshift_v2full.txt` (~2 min) |
| `stage4_c5.py` | raw `Ĉ₅` image (proxy), leading-Regge model raw vs de-shifted, 3-shift de-shift → `cache/stage4_c5.{txt,log}` | 18 min |
| `stage4_save.py` | writes `functionals_q4_deshift.json` (exact `R_j`, weights `χ_value`, `χ_deriv` for EE and scalar, `dν = 0.02`, `N = 2048`) | s |
