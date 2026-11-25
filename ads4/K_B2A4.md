# K_B2A4: the `|ĝ₃|²` light functions of `Ĉ₃` and `Ĉ′_EE` computed in closed form; `b₂ = −35/2`, `a₄ = −581/4`, `Λ²H17 = +35.7`

*25 Nov 2026, research agent. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/kb2a4/`
(`kb_celestial.py`, `kb_onepoint.py`, `kb_k2ratio.py`, `kb_h17.py`, `kb_run_all.sh`); outputs in `ads4/code/kb2a4/cache/`, which its own
`.gitignore` ignores. Read-only imports: `ads4/code/kcert/kc_common.py`, `kc_dk0.py`, `kc_grid.py` (run with `python3 -B`, so nothing is
written there) and the published data they read. Compute: about 8 CPU-min, at most 2 processes, one thread, `ulimit -v 3000000`.*

**Grades** as in UNIFIED_THEOREM_v3 / K_CERT: [P], [certified] (Arb ball arithmetic), [D], [D-sketch], [DV], [Q], [NE], [M], [I].
A chain has the grade of its weakest link.

---

## 0. Verdict

**Both coefficients are computed, exactly, not modelled.** H17 is no longer [D-sketch + M]. It is a definite number,
**`Λ²H17 = +35.7`**, with grade **[D]**: a derivation in which every step specific to this note is [P] or [DV], resting on inherited
[D] links. The arithmetic from the light functions to `K` and `τ₀` is [certified]. The model band `±424.97/Λ²` is replaced.

| quantity | result | grade |
|---|---|---|
| **`k = 3` light function** (`Ĉ₃`, `Π̂₅`-type admixtures included) | `−Ĉ₃|_{ĝ₃} = −2πG|ĝ₃|²(ν/Λ)⁴ [1 − (35/2)ν⁻² + (569/16)ν⁻⁴]`, **an exact polynomial** | [D]; shape [DV] to `1.5·10⁻¹³` |
| **`b₂`** (79 % of the old band) | **`−35/2 = −17.5`** (model: `|b₂| = 72`) ; `b₀ = 569/16` | [D] + [DV] |
| **`k = 2` light function** (`Ĉ′_EE`, h1 reading, i.e. what multiplies `Re X_θ`) | `2πG|ĝ₃|²(ν/Λ)⁶ g₂(ν)`, `g₂ = N₁₀(ν)/[64ν⁶(16ν⁴ − 728ν² + 969)]`, `N₁₀` below; **rational, not a polynomial** | [D]; algebra [P] |
| **`a₄`** | **`−581/4 = −145.25`** (model: `|a₄| = 108`); next terms `+21875/16 ν⁻⁴`, `−1129543/64 ν⁻⁶`, … | [D] |
| `k = 2` checks | (i) the same reduction reproduces the **exact** Lemma-R Einstein function of the EE functional, `∝ (16ν⁴ − 440ν² + 273)`, real zeros included; (ii) the `|ĝ₃|²`/Einstein ratio has top coefficient **exactly** `1/(8·240²)`, i.e. the flat `2πG|ĝ₃|²p⁶ : 16πG/p²` with the project's `√(t₄²+α_T²) = 240|ĝ₃|` | [P] (algebra) |
| **`Λ²H17`**: `k = 2` / `k = 3` / total | **300**: `+117.9268 / −82.0072 / +35.9196`; **1000**: `+117.7337 / −82.0165 / +35.7172`; **2000**: `+117.7195 / −82.0172 / +35.7024`; **5000**: `+117.7156 / −82.0174 / +35.6982` | [certified] given the light functions |

**`K = K₀ + δK₀ + H17` and the bound** (`K₀`, `δK₀`, `E` from K_CERT; continuum weights at 300/1000, grid weights at 2000/5000):

| `Δ_g` | **`K`** (replaces the band; ball radius ≤ `3·10⁻¹¹`) | `K − 1` (`Λ²(K−1)`) | **`τ₀(K)`** | certified upper end | K_CERT `τ₀(K_c)` (H17 = 0) |
|---|---|---|---|---|---|
| 300 | **`0.9996364968`** | `−32.389/Λ²` | **`3.132624·10⁻⁶`** | `≤ 3.1327·10⁻⁶` | `3.133256·10⁻⁶` |
| 1000 | **`0.999967523128`** | `−32.380/Λ²` | **`2.891310·10⁻⁸`** | `≤ 2.8914·10⁻⁸` | `2.891361·10⁻⁸` |
| 2000 | **`0.99999189318`** | `−32.379/Λ²` | **`1.928324·10⁻⁹`** | `≤ 1.9284·10⁻⁹` | `1.928332·10⁻⁹` |
| 5000 | **`0.99999870413`** | `−32.378/Λ²` | **`5.326180·10⁻¹¹`** | `≤ 5.3262·10⁻¹¹` | `5.326183·10⁻¹¹` |

* **`K_min`.** There is no band left: `K_min = K` as tabulated. The numerical uncertainty is the ball radius (≤ `3·10⁻¹¹`). The
  remaining uncertainty is structural (grade [D], §7), not a number.
* **Headline.** At `Δ_g = 1000`, **`τ₀ ≤ 2.90·10⁻⁸` holds with margin `5966/Λ²` in `K`**. A failure would need `H17 ≤ −5931/Λ²`, while
  the computed value is `+35.7/Λ²`. `|ĝ₃|Λ⁴ ≤ 120` holds as well. Even `τ₀ ≤ 2.892·10⁻⁸` holds (margin `476/Λ²`). At 300,
  `τ₀ ≤ 3.135·10⁻⁶` (margin `135/Λ²`).
* **The sign of H17 is positive.** The `k = 2` and `k = 3` corrections partly cancel (`+117.7 − 82.0`), and H17 partly offsets δK₀
  (`−68.1`). Net, `K − 1 = −32.38/Λ²` at every gap, to 0.03 %.
* **Size against the model.** The computed H17 is **8.4 %** of the old model band (`424.97/Λ²`). The leading linear form
  `−0.8106a₄ + 4.6867b₂ = +117.74 − 82.02 = +35.72` agrees with the exact `+35.717` at 1000, to `9·10⁻⁵`.
* **K_CERT §2.3 is settled.** K_CERT asked whether the `|ĝ₃|²` term enters shifted or unshifted, i.e. whether `δK₀` is a convention.
  The reduction gives the per-ν raw functional directly. Its Einstein value has the shifted (Lemma-R) form, both real zeros included,
  and its `|ĝ₃|²` value is a fixed rational multiple of that. So the "`Ĉ′`-primary" end (`δK₀ = 0`) is excluded: `δK₀` is physical,
  and `K = K₀ + δK₀ + H17` without a convention range.
* **Weighted norm (K_CERT §4.4).** `k = 3`: depends on `t₄² + α_T²` only, as before [P given Q2, Q8, F]. `k = 2`: the `W×W` value is
  invariant under the `α_T` phase if the parity-odd first-moment matrix obeys the same `M̂ = ½{Γ,Ĥ}` as the parity-even one
  [D-sketch, as in K_CERT].

**What is not claimed.**
* The grade is [D], not [P]: inherited links (F/F₂ factorisation, DZ with its parity-odd caveat, LIGHTSIDE_DIST's `Ω_OPE = Lemma R`)
  and this note's reduction lemma (§2), which is [D] with strong checks for `k = 2`.
* For `k = 3`, two new steps have no independent end-to-end check (§7): the `j′ = 1` sign structure, and the identification of
  `N̂₃` with the `J = 0` heavy action.

---

## 1. What is computed

**The light value.** By LIGHTSIDE_DIST (`Ω_OPE = Lemma R` on tree terms [D]) and DZ (the `[TT]` double zeros, [D] with the H0 §15.1
parity-odd caveat), the `O(|ĝ₃|²)` light value of each piece is the action of the per-ν CLS-D functional on the `W×W` part of the single
`s`-channel `T` block.

**Its factorisation.** By C1_FINAL's F₂ (`k = 2`) and H7's F (`k = 3`) [D], that block factorises through the light-ray one-point
matrices:

    T-block term = ⟨Ψ₄| A₃ P_T A₁ |Ψ₂⟩,   A ∈ {E, E⁽¹⁾},   E(n) = L[T](x₀ = ∞, (1,n)),   E⁽¹⁾ = ∫du u T = ½{E, −i n̄·D}.

**The task.** This was known pointwise. The missing step (K_CERT §4.2) was to turn it into an exact function of ν. It needs the
celestial and dual kernels and the `P`-measure, which this note supplies:
* the `k = 3` piece in §3;
* the `k = 2` piece, with first moments, in §4;
* `K` and `τ₀` in §5.

---

## 2. The reduction lemma: the functional at fixed ν is a rest-frame celestial integral [D]

**Setting.** CLS-D's functional (`spinningfunctional_definition`) integrates
`(dual)⁻¹ × ⟨P̃P̃P_δ⟩ × ⟨O₄[L₃,L₁]O₂⟩` over `x₂, x₄, x₀, z₀, z₁, z₃`, modulo `SO(3,2)`. Fix `x₀ = ∞` (i⁰).

**Step 1: Parseval.** The light rays are the energy detectors `E(n)`, which commute with `P̂`, and `f` is at most affine in the
retarded times. Translations then reduce `∫d³x₂d³x₄` to `∫d³x₂₄`. Parseval turns this into
`∫_{V⁺}d³P D̃_f(P)·G̃(P)`, with:
* `G̃(P)` the one-point-matrix product at momentum `P`;
* `D̃_f` the Fourier transform of (dual structure × `f`-weights).

**Step 2: gauge fixing.** Lorentz transformations and dilatation fix `P = (1,0,0)`. The integral over the remaining two variables
(three celestial points on `S¹`, modulo `SO(2)`) is exact. No Regge or large-ν expansion enters.

**Ingredients in this frame.**
* **Dual factor.** It is `(z₀·P)^{δ−1}`, constant. Covariance makes the dual's helicity entries `d_{h₄h₂}e^{i(h₄−h₂)φ₀}`.
* **Celestial structure (`δ₁ = δ₃ = 2`, `δ = ½+iν`).** It is
  `(2e₁₃)^{1+δ/2}(2e₁₀)^{−δ/2}(2e₃₀)^{−δ/2}`, with `e_{ij} = 1 − cos φ_{ij}` and `2e = chord²`.
  * For `j′ = 1` it is multiplied by the orientation sign of `(z₀,z₁,z₃)`, as in H0 §R and S1_S4 §1.1.
* **One-point matrices.** In the rest frame they are `⟨T,h′|E(n)|O,λ⟩ ∝ c_{h′λ} e^{i(λ−h′)φ_n}` (rotation covariance). The `T` block
  therefore enters through celestial harmonics of `φ₃₁`:
  * `m = 0`: Einstein (diagonal);
  * `m = ∓4`: `W×W` (flip × flip, `λ = −h`);
  * `m = ∓2`: a heavy `J = 0` block.

**Three facts make the ratio computable.**
* For `f = −1` (`k = 3`), each block's action is a block-independent `U(ν)` times celestial Fourier coefficients.
* `U(ν)` contains the dual's Fourier constant and all ν-only normalisations. It cancels in every ratio used below.
* For `k = 2`, `f` adds `P`-derivatives and `n̄·x` weights (§4). These are still computable in closed form.

**The Fourier coefficients** (`kb_celestial.py`). With `φ = φ₃ − φ₁`:

    k^{(j′)}(m) = ∫dφ e^{imφ}K^{(j′)}(φ) = (2π)² Σ_n b_n(−δ/2)² b_{n−m}(1+δ/2),
    b_n(s) = e^{−iπn}Γ(2s+1)/(Γ(1+s−n)Γ(1+s+n)),   n ∈ ℤ (j′ = 0),  n ∈ ℤ+½ (j′ = 1).

* The `j′ = 1` signed chord powers are antiperiodic, hence the half-integer modes.
* They agree with direct double quadrature to `6·10⁻¹²`, `9·10⁻¹⁰`, `5·10⁻¹²` and `3·10⁻¹⁰` (four cases, §6) [DV].
* For `j′ = 0` and harmonic 0 every monomial integral is Dixon's identity:

      ∫∫∫|x₁−x₃|^{2A}|x₁−x₀|^{2B}|x₃−x₀|^{2C} (on S¹)
          = (2π)³ Γ(1+A+B+C)Γ(1+2A)Γ(1+2B)Γ(1+2C) / [Γ(1+A)Γ(1+B)Γ(1+C)Γ(1+A+B)Γ(1+B+C)Γ(1+C+A)]   [P],

  checked against the Fourier sums to `8·10⁻¹⁴` and `2·10⁻¹²`.

**Grade: [D].** The steps are covariance, Parseval and the gauge quotient. They rest on C1_FINAL/H7's factorisation and on the
`x₀ = ∞` forms of the structures. The `k = 2` checks of §4.4 test the whole chain.

---

## 3. `k = 3`: `b₂` (`kb_celestial.py`)

**The structures.**
* **Weight and ordering.** `Ĉ₃ = Re[Ψ_{3,ν}[o_H q]/N̂₃]`, with `j′ = 1`, `f = −1` and no de-shift. Here `Ψ = 2Ψ⁺`: the `1↔3` image
  doubles, and the kernel is odd under `1↔3`.
* **The `W×W` term.** By H7 §4 (F), it is `Tr_d([Ĥ₃,Ĥ₁])` with `[Ĥ_W(φ₃),Ĥ_W(φ₁)] = −2i|b|² sin(4φ₃₁)σ_z`. The dual `o_H q` is
  parity-odd, so its diagonal satisfies `d₋₋ = −d₊₊`. Its flip entries cannot contribute, because the `W×W` commutator is diagonal.
  The result is `U₃d₊₊|b|²·(−4)k^{(1)}(4)`.
* **The normalisation `N̂₃`.** It is fixed by the heavy action, "real positive at `++`" (S1_S4 §1.4). Take a heavy parity-definite
  `J = 0` block:
  * `⟨T,h|E(n)|O⟩ = c_h e^{−ihφ}`, with `|c₊| = |c₋|`;
  * its preserving part gives `U₃d₊₊|c|²·(−4)k^{(1)}(2)`;
  * its flip part, `∝ e^{ih(φ₁+φ₃)}`, is `1↔3`-symmetric and is killed by the odd kernel;
  * the action factorises exactly into (ν) × (Δ), so the `m → ∞` limit at fixed ν is ν-independent up to a constant;
  * hence `N̂₃(ν) ∝ U₃(ν)d₊₊(ν)k^{(1)}(2)`.

**Consequence.**

    Ĉ₃|_{W×W}(ν) = (const > 0) · G₃(ν),    G₃ := k^{(1)}(4)/k^{(1)}(2).

`G₃` is real. Its sign at large ν is that of a heavy `++` state, so `−Ĉ₃|_{ĝ₃} ∝ −ν⁴` (§0) ✓.

**Result [DV].** At `ν = 0.1, 0.5, 1, 2, 3, 5.5, 7`, with `N = 3000` modes:

    G₃(ν) = (ν⁴ − (35/2)ν² + 569/16)/400      to ≤ 1.5·10⁻¹³ (truncation-limited; ≤ 3.5·10⁻¹² at N = 1500).

At `ν = 12` and `20` the truncation tail is larger: `1.7·10⁻¹⁰` and `1.1·10⁻⁵` at `N = 3000`. Both fall by about 20 per doubling of
`N`, as a truncation error should, and the imaginary parts shrink with them.

The values are exact binary fractions, for example `G₃(1) = 61/1280` and `G₃(2) = −59/1280`.

**Structural cross-check.** The `j′ = 0` ratio `k^{(0)}(4)/k^{(0)}(2)` is the same polynomial, and
`k^{(0)}(2)/k^{(0)}(0) = (4ν²+1)(4ν²+25)/1152`.

**Normalisation.** The flat limit `−Ĉ₃|_{ĝ₃} → −2πG|ĝ₃|²(ν/Λ)⁴` [Q: CHMRS `BknuFlatSpaceLimit`; CLS-D dictionary, H0 §15.2] fixes it.
Then:

    −Ĉ₃|_{ĝ₃}(ν) = −2πG|ĝ₃|²(ν/Λ)⁴ [1 − (35/2)ν⁻² + (569/16)ν⁻⁴],   b₂ = −35/2,  b₀ = 569/16,

**exactly**, with no `ν⁻⁶` or further term.
* So the polynomial ansatz of H0 §18 / K_CERT §4.3 is exact for `k = 3`.
* The `Π̂₅` admixture is not needed separately: the reduction evaluates the full functional.
* **`H17_{k=3}`** follows as an exact rational: `H17_{k=3} = (2/Λ²)[−(35/2)N₁ + (569/16)N₀/Λ²]`, i.e. `Λ²H17_{k=3} = −82.0072`
  (300) and `−82.0165` (1000).

---

## 4. `k = 2`: `a₄` (`kb_onepoint.py`, `kb_k2ratio.py`)

### 4.1 The subtraction factor at `x₀ = ∞`

With `x = x₂₄` and `X = (x₂+x₄)/2`, C1_FINAL's exact `f` (§2.1) is

    f = (2/(e₁₃x²)) [ (n̄₃·x)(n̄₁·X − u₁) − (n̄₁·x)(n̄₃·X − u₃) ]   (translation invariant).

So `f L₃L₁ → (2/(e₁₃x²))[−(n̄₃·x)E₃E⁽¹⁾_{X,1} + (n̄₁·x)E⁽¹⁾_{X,3}E₁]`, with the first moments taken about the midpoint `X`.

**Momentum space** (states `T(x,z)|0⟩ ↔ e^{iP·x}w(P)Π_TT[zz]`, `w² = (P²)^{3/2}`). Write `E⁽¹⁾ = −i[Ĥ n̄·∇ + ½(n̄·∇Ĥ)]` (F₂).
* The EE `(2,4)` contraction is the trace over TT tensors (C1_FINAL §2.4). Projector-derivative terms drop by `Π(∂Π)Π = 0`.
* The midpoint subtraction leaves `±(n̄·x)/2`.

The `T`-block integrand is then

    (2/e₁₃) w² { T₃₁ (n̄₃·∂)(n̄₁·∂)𝔉 + (n̄₃·∂𝔉)[T₃₁λ₁ + τ₁/2] + (n̄₁·∂𝔉)[T₃₁λ₃ + τ₃/2] },

with:
* `T₃₁ = Tr(Ĥ₃Ĥ₁)`, `τ₁ = Tr(Ĥ₃ n̄₁·∂Ĥ₁)`, `τ₃ = Tr((n̄₃·∂Ĥ₃)Ĥ₁)`, `λ = n̄·∂ log w = 3/2`;
* `𝔉(P) = ∫d³x e^{iPx}D(x)/x² ∝ (z₀·P)^{α}(P²)^{κ}` (covariance), with `α = δ−1` and `κ = −11/4 − iν/2`.
  * `D(x) ∝ (z₂·z₄)²(z₀·x)^{δ−1}(x²)^{11/4−iν/2}` is EE times the `Δ̃ = 0` scalar dual, of `x`-degree 5 as conformal weight requires.
  * `x`-weights become `P`-derivatives.

At rest:

    ∂𝔉/𝔉 = αe + 2κ,   ∂₃∂₁𝔉/𝔉 = (αe₁+2κ)(αe₃+2κ) − αe₁e₃ + 2κ(e₁₃ − 2).

### 4.2 One-point traces [DV] (`kb_onepoint.py`)

**The matrices** (covariant HM, `d = 3`):
* `Ĥ_E = ωΠ_TT`;
* `Ĥ_W = ω(|q⟩⟨q| − Π/4)`, with `q = Π_TT[NN]`, `N` the rest-frame direction of `n̄`, `ω = (P²)^{3/2}/(n̄·P)²`;
* `Ĥ = Ĥ_E + t₄Ĥ_W`.

**Calibration.** Flips are `t₄/4` (so `|flip|/diag = 1` at `|t₄| = 4`), with phase `arg = −1.48` at `φ = 0.37`. Both agree with
C1_FINAL §2.3(e).

**Traces at 40 digits** (finite-difference derivatives, DFT in `φ₃₁`):

| | `T₃₁` | `τ₁ = τ₃` |
|---|---|---|
| Einstein | `2` | `6` |
| `W×W` | `cos 4φ/8` | `3 cos 4φ/8` |

All other Fourier coefficients are `≤ 10⁻⁴²`. So the first-moment structure of `Ĥ_W` adds nothing beyond its `ω` factor, and

    B_W = (cos 4φ₃₁/16) · B_E   pointwise,
    B_E = 2Y,   Y = α(α−1)e₁e₃ + (2ακ+3α)(e₁+e₃) + 2κe₁₃ + 4κ² + 8κ.

### 4.3 Closed form [P given §4.1–4.2]

**Method.** Write `cos 4φ = T₄(1 − e₁₃)`. Every term is then a Dixon integral with integer-shifted exponents, so every ratio is
rational in ν (`kb_k2ratio.py`, sympy, exact).

**Einstein:**

    Ψ_ν[E×E] = U(ν) I₀(ν) (16ν⁴ − 440ν² + 273)/64.

**`W×W` over Einstein** (`t₄²` stripped):

    Θ(ν) = (4ν²+1)·N₁₀(ν) / [117964800 (16ν⁴ − 440ν² + 273)],
    N₁₀ = 1024ν¹⁰ − 195328ν⁸ + 8229504ν⁶ − 90780512ν⁴ + 251603732ν² − 114141375.

**Independent numerical check.** Fourier sums with the harmonic `m = ±4` inserted directly, with no Dixon formula and no `T₄`
expansion, agree to truncation level: `2.5·10⁻⁵` at `ν = 1.3` and `4·10⁻³` at `ν = 5`. These sums converge only like `N^{−3/2}`.

### 4.4 Two exact checks (the reason to trust §2 and §4.1–4.2)

1. **The Einstein function is reproduced, zeros included.**
   * The Lemma-R value of the EE functional uses H7 §5 with the EE leading-Regge coefficients `R_∓ = R^{sc}_∓ϱ(ν±i)` (the same as
     `prl/data/kernels.json`):

         R₋D(ν+i) + R₊D(ν−i) = −32(16ν⁴ − 440ν² + 273)/((4ν²+9)(4ν²+49)).

   * Its ratio to the reduction's result is `−2048/((4ν²+9)(4ν²+49))`. That has no real zeros or poles and is absorbed by the Γ-type
     `U(ν)I₀(ν)/N_com`.
   * The quartic `16ν⁴ − 440ν² + 273`, with real zeros at `ν = 0.797` and `5.183`, is produced by the interplay of the first-moment
     terms, the `f`-weights, the dual and the celestial structure. A sign or factor error in any of them moves those zeros.
   * Two independent routes are compared here: Lemma R (Regge residue, plane waves) and the energy-channel factorisation.
2. **The top coefficient equals the flat value exactly.**
   * `Θ = ν⁸/460800 + …`, and `460800 = 8·240²`.
   * The flat ratio `2πG|ĝ₃|²p⁶ : 16πG/p²` is `|ĝ₃|²ν⁸/8` (AdS units). With `t₄² + α_T² = 240²|ĝ₃|²` (C1_FINAL's normalisation list)
     this is exactly `t₄²Θ`.
   * The 240 comes from the flat three-graviton dictionary, independently of everything here. **The cross-check against `K₀`
     requested in the task passes exactly.**

### 4.5 The `k = 2` light function and `a₄`

**Per-ν value.** The per-ν Einstein value multiplying `X_θ = θχ_re(R₋+R₊)` is `16E_EE/(R₋+R₊)` (K_CERT §3, kc_grid). The `|ĝ₃|²` value
multiplying `X_θ` is therefore that times `57600|ĝ₃|²Θ`. In the h1 units of K_CERT §4.3:

    −Ĉ′_EE|_{ĝ₃} ≙ 2πG|ĝ₃|²(ν/Λ)⁶ g₂(ν),   g₂ = N₁₀(ν)/[64ν⁶(16ν⁴ − 728ν² + 969)]
                                           = 1 − (581/4)ν⁻² + (21875/16)ν⁻⁴ − (1129543/64)ν⁻⁶ + …

**So `a₄ = −581/4`.**

**Two remarks.**
* **The cancellation.** `(4ν²+1)` and the Einstein quartic cancel between `Θ` and `E_EE`. What remains in the denominator is
  `16ν⁴ − 728ν² + 969`, the factor of the forward normalisation `R₋+R₊`. It has real zeros (`ν = 1.17, 6.64`), so in the h1 reading
  `g₂` has poles there.
* **No poles in what enters `K`.** `X_θ ∝ (R₋+R₊)` cancels the poles. The product
  `X_θ·(light) = θχ_re·P(ν)/Λ⁶`, with `P = −(4ν²+1)N₁₀/(8(4ν²+9)(4ν²+49))`, is regular.
* **Consequence.** The polynomial ansatz `a₄, a₂, a₀` of K_CERT §4.3 is only asymptotic for `k = 2`. `H17` is therefore computed from
  `P` directly (§5), not from `a₂, a₀`.

---

## 5. Certified H17, `K`, `τ₀` (`kb_h17.py` → `cache/kb_h17.txt`)

**Formulas.**

    H17_{k=2} = −(2/Λ⁷) ∫₀^∞ θχ_re [P(ν) − (R₋+R₊)ν⁶] dν
              = −(2/Λ⁷) Σ_k dν Re[X_k/(R₋+R₊)(ν_k)] [P − (R₋+R₊)ν⁶](ν_k)     (grids)
    H17_{k=3} = (2/Λ²)[−(35/2)N₁ + (569/16)N₀/Λ²]                               (exact rational)

where `P − (R₋+R₊)ν⁶ = (4ν²+1)(148736ν⁸ − 8167488ν⁶ + 90780512ν⁴ − 251603732ν² + 114141375)/(8(4ν²+9)(4ν²+49))`.

**Method.** The same machinery as K_CERT:
* `acb.integral` on `[0, Λ+60]` in steps of 10;
* the analyticity checks of `Chi`;
* a rigorous tail beyond `Λ+60`. For `ν ≥ 300`, `|P − (R₋+R₊)ν⁶| ≤ 5000ν⁶`, since the numerator is `≤ 148737ν⁸` and the denominator
  `≥ 32ν⁸`. The tail is `< 10⁻²¹⁷`.

**Results.** `Λ²H17` (k=2 / k=3 / total):

| `Δ_g` | weight | `k = 2` | `k = 3` | **total** | `K` | `τ₀` |
|---|---|---|---|---|---|---|
| 300 | continuum | `117.926780 ± 8·10⁻⁷` | `−82.007162` | **`35.919618`** | `0.9996364968` | `3.1326244·10⁻⁶` |
| 300 | grid | `117.926780142` | `−82.007162` | `35.91961808` | `0.99963649682` | `3.1326277·10⁻⁶` |
| 1000 | continuum | `117.733686 ± 7·10⁻⁷` | `−82.016492` | **`35.717194`** | `0.999967523128` | `2.8913095·10⁻⁸` |
| 1000 | grid | `117.733686071` | `−82.016492` | `35.71719379` | `0.99996752316` | `2.8913186·10⁻⁸` |
| 2000 | grid | `117.719538` | `−82.017179` | **`35.702359`** | `0.99999189318` | `1.9283237·10⁻⁹` |
| 5000 | grid | `117.715586` | `−82.017371` | **`35.698215`** | `0.99999870413` | `5.3261796·10⁻¹¹` |

* Grid and continuum agree to `< 10⁻⁸`.
* `τ₀` uses `E` from K_CERT, the realised Einstein light value of each weight.

**Margins** `(K − K*(T))Λ²` for `τ₀ ≤ T`:

| `Δ_g` | `T` | margin |
|---|---|---|
| 300 | `3.135·10⁻⁶` | `+134.9` |
| 300 | `3.14·10⁻⁶` | `+417.9` |
| 1000 | `2.892·10⁻⁸` | `+476.0` |
| 1000 | **`2.90·10⁻⁸`** | **`+5966.3`** |
| 2000 | `1.93·10⁻⁹` | `+6934.9` |
| 5000 | `5.327·10⁻¹¹` | `+7695.5` |

**Quartic.** Unchanged: `K` enters (2b) only through `κ + εK/8 > 0` (K_CERT §5).

---

## 6. Checks

| check | result | grade |
|---|---|---|
| Fourier representation of the celestial kernels vs direct double quadrature (`ν = 0.7, j′=1, m=2`; `1.3, j′=1, m=4`; `1.3, j′=0, m=0`; `1.3, j′=0, m=4`) | `6·10⁻¹²`, `9·10⁻¹⁰`, `5·10⁻¹²`, `3·10⁻¹⁰` | [DV] |
| Dixon closed form vs Fourier (`j′ = 0`, `m = 0`; `ν = 0.9, 2.5`) | `8·10⁻¹⁴`, `2·10⁻¹²` | [P] + [DV] |
| `k^{(1)}(4)/k^{(1)}(2) = (ν⁴ − 35ν²/2 + 569/16)/400`, 7 values `ν ≤ 7`; imaginary part | `≤ 1.5·10⁻¹³`; `≤ 10⁻¹³` relative (at `ν = 12, 20`: truncation tail `1.7·10⁻¹⁰`, `1.1·10⁻⁵`, falling ×20 per doubling of `N`) | [DV] |
| `j′ = 0`: `k(2)/k(0) = (4ν²+1)(4ν²+25)/1152`, `k(4)/k(2)` = the same `G₃` polynomial (`ν = 0.5, 2, 4`) | `≤ 1.5·10⁻¹¹`, `≤ 1.3·10⁻¹⁰` | [DV] |
| HM normalisation and phase of `Ĥ_W` (flip `t₄/4`, `arg = −1.48` at `φ = 0.37`) | exact; matches C1_FINAL §2.3(e) | [DV] |
| first-moment traces `τ_W = 3cos4φ/8`, no extra term | other harmonics `≤ 10⁻⁴²` | [DV] |
| **`k = 2` Einstein: reduction ∝ Lemma-R function `(16ν⁴−440ν²+273)`** | identical; ratio `−2048/((4ν²+9)(4ν²+49))` | [P] (algebra) |
| **`k = 2` top coefficient `1/460800 = 1/(8·240²)`** | exact | [P] (algebra) |
| `Θ` by Fourier sums without Dixon or `T₄` | `2.5·10⁻⁵` (`ν = 1.3`), `4·10⁻³` (`ν = 5`): agreement at truncation level | [DV] |
| leading linear form `−0.8106a₄ + 4.6867b₂` vs exact `Λ²H17` | `35.720` vs `35.717` (1000); `35.956` vs `35.920` (300) | consistency |
| continuum vs grid H17 (300, 1000) | `< 10⁻⁸` | [certified] |

---

## 7. Where it could fail

1. **Inherited links, shared with the Einstein term and `c₁ = 0`.**
   * F/F₂ (factorisation through `P_T`, `E⁽¹⁾ = ½{E, −in̄·D}`) [D + DV].
   * DZ (`[TT]` double zeros; parity-odd odd-twist caveat [D-sketch], H0 §15.1).
   * `Ω_OPE = Lemma R` [D].
2. **The reduction lemma (§2) [D].**
   * It uses the `x₀ = ∞` forms of the structures and Parseval on the future cone.
   * For `k = 2` both exact checks of §4.4 test it end to end. Only the overall constant `U(ν)` is untested, and it cancels.
3. **`k = 3`-specific steps without an end-to-end check.**
   * (a) **The `j′ = 1` structure.** It is "orientation sign × `j′ = 0`", the project's own definition (S1_S4 §1.1), checked there only
     at Regge and the crystal.
   * (b) **The normalisation.** `N̂₃ ∝` the `J = 0` heavy action. This uses exact (ν)×(Δ) factorisation of that action and the
     equality of the forward and leading-Regge normalisations (S1_S4: the image is exactly `σ_hc₃P_l`).
   * (c) **The dual's helicity entries.** Only `d₋₋ = −d₊₊` for `o_H q` is used, and only its diagonal enters.
   * **Supporting evidence.** (i) The `j′ = 0` analogue gives the same `m = 4/m = 2` polynomial. (ii) The leading power `ν⁴` and the
     sign come out as in the flat sum rule.
   * **Tolerance.** If `b₂` were wrong by a factor 10 (`−175`), the headline `2.90·10⁻⁸` would still hold, with margin `5966 − 738`.
4. **Parity-odd `⟨TTT⟩` in `k = 2`.** The α_T² part needs `M̂_{α_T} = ½{Γ,Ĥ_{α_T}}`. That is known only on the diagonal [D-sketch],
   as in K_CERT §4.4. The `t₄²` part is [D].
5. **Not done.**
   * A Witten-diagram (Regge-limit) evaluation of the `R³` exchange. The task's "subleading Regge data" are bypassed, not computed.
   * A heavy-block check of `N̂₃` against S1_S4's numerical `c₃(ν)` (different conventions).

---

## 8. Proposed replacement text (not applied)

**K_CERT §0 / UNIFIED_THEOREM box.**
> "H17 computed (K_B2A4): `Λ²H17 = +35.92 / +35.72 / +35.70 / +35.70` (300/1000/2000/5000) [D; arithmetic certified]. `b₂ = −35/2`
> (`Ĉ₃`'s light function is exactly `−2πG|ĝ₃|²p⁴(1 − 35/(2ν²) + 569/(16ν⁴))`), `a₄ = −581/4` (`Ĉ′_EE`'s is rational, `g₂ =
> N₁₀/(64ν⁶(16ν⁴−728ν²+969))`). `K = 0.99963650 / 0.99996752`, `τ₀ = 3.13262·10⁻⁶ / 2.89131·10⁻⁸`. The model band `±424.97/Λ²` [M] is
> withdrawn."

**Supplement, Table S2, row L2.** Replace "H17 … modelled `|H17| ≤ 424.97/Λ²` [D-sketch + M]" by

```latex
H17 $=+35.92/\Lambda^2$ ($300$), $+35.72/\Lambda^2$ ($1000$): $\hat C_3$'s $|\hat g_3|^2$ light function is
$-2\pi G|\hat g_3|^2p^4(1-\tfrac{35}{2}\nu^{-2}+\tfrac{569}{16}\nu^{-4})$ exactly, $\hat C'_{EE}$'s is rational with
$a_4=-581/4$ (energy-channel $T$ block through light-ray one-point matrices, rest-frame celestial integrals) & D &
H17 [D] (was D-sketch\,+\,M); $K=0.99963650$ / $0.99996752$ [certified given D]\\
```

**Item (1) and the Letter.**
* The `τ₀` values become `3.1326·10⁻⁶` (300) and `2.8913·10⁻⁸` (1000).
* The Letter's "`2.89·10⁻⁸`" and "`|ĝ₃|Λ⁴ ≤ 120`" are unchanged.
* "`K = 1 + O(Λ⁻²)`" may be sharpened to "`K = 1 − 32.4/Λ²`".

---

## 9. Reproduce

    cd ads4/code/kb2a4 && ./kb_run_all.sh
        # cache/kb_onepoint.txt, kb_k2ratio.txt, kb_celestial.txt (full mode adds the direct quadratures), kb_h17.txt

The scratch exploration (not needed to reproduce) was done in the session scratchpad.
