# K_CERT: the `|ĝ₃|²` light coefficient `K`, with δK₀ certified and H17 reduced to a finite set of numbers that are still open

*25 Nov 2026, research agent. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/kcert/`
(`kc_common.py`, `kc_dk0.py`, `kc_grid.py`, `kc_bounds.py`, `kc_selftest.py`, `kc_run_all.sh`). Outputs go to `ads4/code/kcert/cache/`,
which its own `.gitignore` ignores. Inputs, all read-only: `prl/data/nf02_psi.json`, `prl/data/kernels.json` (`D_k2_EE`),
`prl/data/k2/k2_continuum_D_{300,1000}.json`, `prl/data/d_abel_X_{300,1000}_trunc.txt` and `prl/data/gap/d_abel_X_{2000,5000}_trunc.txt.gz`.
Compute: 45 CPU-s for everything, one process, one thread, `ulimit -v 3000000`.*

**Grades.** [certified]: Arb ball arithmetic (python-flint 0.9) with rigorous quadrature and tail bounds. The other grades are as in
UNIFIED_THEOREM_v3: [P], [D], [D-sketch], [DV], [NE], [M], [I].

---

> **Update (K_B2A4.md, unreviewed): H17 computed in closed form, so the model band is no longer needed. b₂ = −35/2, a₄ = −581/4 (leading), Λ²H17 = +35.7 at every gap [D]; the arithmetic to K and τ₀ is [certified]. K = 1 − 32.38/Λ².**
> * **Method:** at x₀ = ∞ the light-ray factorisation (C1_FINAL F₂, H7 F) plus a Fourier transform to the s-channel rest frame makes each piece's T-block light value an exact celestial-circle integral, evaluated in closed form. No Regge expansion and no Witten diagram are used.
>   * k = 3: −2πG|ĝ₃|²p⁴[1 − (35/2)ν⁻² + (569/16)ν⁻⁴] exactly.
>   * k = 2: rational, with poles at the zeros of R₋+R₊ that X_θ cancels; H17 is integrated from the exact function.
> * **Checks:** the same reduction reproduces exactly the Lemma-R Einstein function, including 16ν⁴ − 440ν² + 273. The flat top coefficient is exactly 1/(8·240²) (= K₀). Kernels match direct quadrature to 1e-10–1e-12. The value agrees with K_CERT's linear form to 3e-3. The physical shifted Einstein value settles δK₀ as physical.
> * **Certified K and headline:**
>   * K = 0.99963650 (300), 0.99996752 (1000), 0.99999189 (2000), 0.99999870 (5000).
>   * **τ₀ = 2.89131·10⁻⁸ (≤ 2.8914·10⁻⁸) at 1000** and 3.13262·10⁻⁶ at 300.
>   * τ₀ ≤ 2.90·10⁻⁸ holds with margin +5966/Λ² in K, and |ĝ₃|Λ⁴ ≤ 120 holds.
> * **Caveats:**
>   * Inherited [D] links: light-ray factorisation, the [TT] double zeros (parity-odd [D-sketch]), and Ω_OPE = Lemma R.
>   * k = 3 has two steps without end-to-end checks: the j′ = 1 structure sign, and N̂₃ = the J = 0 heavy action.
>   * The k = 2 α_T² part is [D-sketch].
>   * A factor-10 error in b₂ would still keep the headline.

## 0. Verdict

**Not achieved: H17 was not proved.** No proven bound on H17 exists, not even a weak one. Both the existing band `±212.5/Λ²` on `τ₀`
(`±424.97/Λ²` on `K`) and its `O(Λ⁻²)` scaling are still [D-sketch + M]. Proving any bound needs one missing computation: the
subleading Regge data of the AdS₄ `R³` stress-tensor exchange and the `Π̂₄`/`Π̂₅` admixtures of `Ĉ′_EE`, `Ĉ₃` (§4). That
computation was not attempted here.

**What was achieved.**
* **δK₀ is certified.** The interval contains the defined number with error `≤ 10⁻¹⁹`.
* **The Einstein light value of the actual weight is certified.**
* **The bound `τ₀(K)` is certified as a function of `K`.**
* **H17 reduces to a few coefficients.** It is now an exact linear form in a handful of unknown light-function coefficients, with
  certified weights.
* **The headline needs only a weak lower bound on H17.** "`τ₀ ≤ 2.90·10⁻⁸`" at `Δ_g = 1000` holds for every
  `H17 ≥ −5.9·10⁻³`, which is **14 times** the modelled band.

| item | `Δ_g = 300` | `Δ_g = 1000` | grade |
|---|---|---|---|
| `K₀ = −2∫ψ₂p⁶ + 2∫ψ₃p⁴` (NF02, `k₀ = 27.48509817` exact) | `1 + 3.1631·10⁻¹¹` (exact rational) | same | [P] |
| **`δK₀ = −2∫(Re X_θ − ψ₂)p⁶dp`** (continuum weight × θ, as defined) | **`−7.666324482972·10⁻⁴ ± 10⁻¹⁹`** | **`−6.830149983909·10⁻⁵ ± 10⁻¹⁹`** | **[certified]** (was [D + DV]) |
| `Λ²δK₀` | `−68.30848441` | `−68.09673194` | [certified] |
| `K_c := K₀ + δK₀` (H17 = 0) | `0.999233367583` | `0.999931698532` | [certified] |
| `E` = Einstein light value of the actual weight (replaces `16L_A`) | `10734.6054682539` (`16L_A·(1 − 5.8·10⁻¹¹)`) | `14339.5890757315` (`16L_A·(1 + 8.3·10⁻¹⁰)`) | [certified] |
| `τ₀(K = 1)` (published) | `3.132055·10⁻⁶` | `2.891263·10⁻⁸` | [certified] |
| **`τ₀(K_c)`** | **`3.133256·10⁻⁶`** | **`2.891361·10⁻⁸`** | [certified] |
| `τ₀` with `K ∈ K_c ± 424.97/Λ²` (also covering the `Ĉ′` convention, §2.3) | `[3.1246, 3.1408]·10⁻⁶` | `[2.8906, 2.8920]·10⁻⁸` | arithmetic [certified]; band [M] |
| **rounded headline, and the H17 it needs** | `τ₀ ≤ 3.15·10⁻⁶` iff `H17 ≥ −1.06·10⁻²` (2.2 × model band); `≤ 3.14·10⁻⁶` iff `H17 ≥ −4.3·10⁻³` (0.90 ×, **not robust**); `≤ 3.2·10⁻⁶` iff `H17 ≥ −4.1·10⁻²` (8.6 ×) | **`τ₀ ≤ 2.90·10⁻⁸` iff `H17 ≥ −5.95·10⁻³`** (14 × model band); `|ĝ₃|Λ⁴ ≤ 120` iff `H17 ≥ −4.13·10⁻³` (9.7 ×); `≤ 2.892·10⁻⁸` iff `H17 ≥ −4.4·10⁻⁴` (1.04 ×) | arithmetic [certified]; the H17 condition is open |
| H17 | not proved; reduced to `Λ²H17 = −0.8106a₄ + 4.6867b₂ + O(Λ⁻²)` (§4.3) | same, `−0.8106a₄ + 4.6867b₂` | **[D-sketch + M], unchanged** |

The `Δ_g = 2000` and `5000` columns use the grid weights, which are the weights used there (README §8):

| | `Δ_g = 2000` | `Δ_g = 5000` |
|---|---|---|
| `δK₀` | `−1.70458454586·10⁻⁵` (`Λ²δK₀ = −68.0811`) | `−2.72469183254·10⁻⁶` (`−68.0764`) |
| `τ₀(K_c)` | `1.928332·10⁻⁹` | `5.326183·10⁻¹¹` |
| model band | `[1.9282, 1.9285]·10⁻⁹` | `[5.3261, 5.3263]·10⁻¹¹` |
| `τ₀ ≤ 1.93·10⁻⁹` / `5.33·10⁻¹¹` | iff `H17 ≥ −1.7·10⁻³` (16 × band) | iff `H17 ≥ −1.4·10⁻³` (84 × band) |

**`K_min`.**
* **No certified `K_min` exists** while H17 is open.
* The certified statement is conditional: `K ≥ K_c − h` whenever `H17 ≥ −h`.
* On the model band: `K_min = 0.994432` (300) and `0.999505` (1000). These come from `K_c − 424.97/Λ²`, together with the
  `Ĉ′`-convention end `K₀ − 424.97/Λ²`.

**Quartic.** `K` does not enter the quartic upper end.
* (2a) uses only `Ω^Q_up`.
* In (2b), `Ω^Q_up + εΩ^EE`, with `c₁ = 0`, `K` appears only in the sign condition `κ + εK/8 > 0`. That condition holds for every
  `K > 0` once `κ > 0`.
* The numbers `72.0406 / 92.8135` and `125.71 / 164.51` are unchanged.
* The Einstein value `2εL_A` moves by `≤ 10⁻⁹` relative.

**The Letter.** Its numbers do not change:
* `2.89·10⁻⁸` and `|ĝ₃|Λ⁴ ≤ 120` hold at `K_c` and on the whole model band.
* Only the phrase "at `K = 1`" should be replaced (§6).

**Net grade of item (1).** Unchanged: [D] without light single traces, plus the H17 band [D-sketch + M]. What changes:
* the δK₀ part of row L2 is now [certified];
* the conclusion "`≤ 2.90·10⁻⁸`" is robust to an H17 error 14 times the modelled size.

---

## 1. What `K` is

**The sum rule.** Item (1) is the root of

    K τ² Λ⁸/57600 ≤ E,    E := the Einstein light value of Ω^EE (= 16L_A up to the realisation of the weight, §3),
    τ ≤ τ₀(K) = (240/Λ⁴) √(E/K).

**The definition of `K`.** `K` is the `|ĝ₃|²` light value of `Ω^EE`, normalised by its flat value:

    K = [light value of the smeared k = 2 piece Re∫X_θ Ψ̂_EE,ν/(h₋+h₊) on the |ĝ₃|² term]
      + [light value of ∫ψ₃ ½Λ⁶Ĉ_{3,Λp} on the |ĝ₃|² term]   (normalised so that the flat value is K₀).

**The flat limit.** The light functions are flat:
* `−Ĉ′_EE|_{ĝ₃} ∝ +(ν/Λ)⁶`;
* `−Ĉ₃|_{ĝ₃} ∝ −(ν/Λ)⁴` (H0 §15.2, §18).

This gives the flat value `K₀ = −2∫ψ₂p⁶ + 2∫ψ₃p⁴`.

| piece | flat value |
|---|---|
| `k = 2` piece | `−0.283854` |
| `k = 3` piece | `+1.283854` |
| **sum `K₀`** | `1 + 3.16·10⁻¹¹` |

The `k = 3` piece dominates.

**The decomposition `K = K₀ + δK₀ + H17`.**
* **δK₀.** In the published ("h1") convention the `k = 2` weight `Re X_θ` multiplies the flat `p⁶` image. This gives the shift
  `δK₀ = −2∫(Re X_θ − ψ₂)p⁶`. Only the `k = 2` piece is de-shifted; `ψ₃` is used directly.
* **H17.** Everything else is H17: the AdS (`O(ν⁻²)`-relative) parts of the two light functions.

---

## 2. δK₀ certified (`kc_dk0.py`, `kc_grid.py`)

### 2.1 Definition and method
* **Definition.** The weight is the published continuum weight times the cutoff,
  `X_θ = θ·χ_re·(R₋+R₊)`, with `θ = θ(Λ+100, 10)` (Supplement §S2). Then

      δK₀ = −(2/Λ⁷) ∫₀^∞ (θχ_re(R₋+R₊) − ψ)ν⁶ dν,    ψ(ν) = ψ₂(ν/Λ) on [0, Λ], 0 beyond.

  This is exactly `code/cert/k2_cont.py::checks()`, which printed `−7.666324482·10⁻⁴ / −6.830149984·10⁻⁵`.
* **`χ_im` does not enter.** `R₋+R₊ = 2Re R₋` is real on ℝ and the functional takes `Re`, so `χ_im` drops out of every real light value.
* **Data, read exactly.**
  * `χ_re`: the Gaussian-RBF coefficients and centres are IEEE doubles, and `y_j` likewise.
  * `R_∓`: exact Gaussian-rational functions.
  * `ψ₂`, `ψ₃`: exact rationals, divided by the exact decimal `k₀`.
  * `Λ = √((Δ_g−2)(Δ_g−1))`, computed as a ball.
* **Quadrature.**
  * Arb's `acb.integral`: Gauss–Legendre with Bernstein-ellipse error bounds. The integrand is checked analytic on every complex ball.
  * The integrand is written as `θ·χ·ν⁵·[ν(R₋+R₊)]`, where `ν(R₋+R₊)` is a polynomial ratio with no pole at 0.
  * Poles at `±i/2` (in `Q`), `±3i/2`, `±7i/2` and the `log` branch points at `±i` make the ball non-finite automatically. The
    `log(ν²+1)` cut is tested explicitly.
  * Pieces: `[0, Λ]` in steps of 10, then `[Λ, Λ+60]`. The breakpoint at `Λ` is needed because `ψ` has compact support.
* **The Gaussian sum on a ball.**
  * Only centres within 10 of the ball's real range are summed.
  * All omitted centres together are bounded by `4·max|x|·e^{−(D²−Y²)/(2σ_max²)}/(1 − e^{−Dh/σ_max²})`, with `D = 10 − 10⁻⁶`,
    `σ_max = 0.9`, `h = 0.15`, `Y = sup|Im z|`. This bound is added as a ball.
* **Tail beyond `Λ + 60`.** Rigorous elementary bound: `|R_∓| ≤ 18ν²`, `|Q| ≤ 1/ν²`, and Gaussian/log-derivative estimates. It is
  `< 10⁻²⁰⁰` and is added as a ball.
* **Precision.** 128-bit.
* **Self-tests** (`kc_selftest.py`, all pass):
  * windowed vs full sum `≤ 3·10⁻²⁴`;
  * vs the published table `≤ 1.3·10⁻¹³` of `max χ`;
  * small-`ν` series branch vs direct evaluation;
  * a closed-form Gaussian integral to `10⁻¹⁹`;
  * agreement with CERT_CLOSE's double-precision δK₀ to `1·10⁻¹¹` relative.

### 2.2 Results

| `Δ_g` | weight | **δK₀** [certified] | `Λ²δK₀` | `K₀ + δK₀` | `(K₀+δK₀)^{−½}` |
|---|---|---|---|---|---|
| 300 | continuum (definition) | `−7.66632448297211·10⁻⁴ ± 7·10⁻²⁰` | `−68.30848441` | `0.999233367583334` | `1.00038354` |
| 1000 | continuum (definition) | `−6.83014998390935·10⁻⁵ ± 3·10⁻²⁰` | `−68.09673194` | `0.999931698531792` | `1.00003415` |
| 300 | grid `d_abel_X_300_trunc` | `−7.66632416821·10⁻⁴` | `−68.308482` | | |
| 1000 | grid `d_abel_X_1000_trunc` | `−6.83014691062·10⁻⁵` | `−68.096701` | | |
| 2000 | grid (weight used at 2000) | `−1.70458454586·10⁻⁵` | `−68.081141` | `0.999982954186` | `1.00000852` |
| 5000 | grid (weight used at 5000) | `−2.72469183254·10⁻⁶` | `−68.076431` | `0.999997275340` | `1.00000136` |

* **Grid values.** The grid value is the exact light value of the grid functional `Re Σ_k dν X_k Ψ̂/(h₋+h₊)`. It is a finite sum of
  exact decimals, with no quadrature error.
* **Grid vs continuum.** The two agree to `4·10⁻⁸` relative (300) and `4.5·10⁻⁷` relative (1000).
* **As a function of Λ.** `Λ²δK₀` decreases monotonically in magnitude, from `−68.31` (300) to `−68.08` (5000). This is consistent
  with the first-order value `−68.076`.
* **No certified statement "for all `Δ_g ≥ 1000`" is possible.** No functional is defined at gaps other than 300, 1000, 2000 and
  5000: there is no published weight elsewhere (README §8). The supplement's "`δK₀ = −68.08/Λ²`" is accurate to 0.3 % at 300 and to
  0.02 % at 1000.

### 2.3 δK₀ is a convention; its physical content belongs to H17 [certified computation; D for the identity]

**The alternative convention.** Suppose instead that `Ĉ′_EE` is taken as the primary object, with a flat `|ĝ₃|²` image, so that
`Ω = ∫ψ₂Re Ĉ′`. This is H7_I14 §23.2's remark: "this shift is then part of H17".

**What the de-shift identity gives in that convention.**
* Take the shifted form of the identity (CERT_CLOSE §3.2), `W(t) = ½∫E(μ)cos μt`, and its even `t`-moments at `t = 0`.
* The `|ĝ₃|²` value is then the same as for the Einstein pole (V3): `∫θχ[R₋(ν+i)⁶ + R₊(ν−i)⁶]`.
* Certified results:

| `Δ_g` | `Λ²δK₀^{Ĉ′}` | weight |
|---|---|---|
| 300 | `−5.93·10⁻⁸` | continuum |
| 1000 | `+2.77·10⁻⁸` | continuum |
| 2000 | `+6.9·10⁻⁵` | grid |
| 5000 | `+4.4·10⁻⁴` | grid; the grid satisfies the identity only to `t ≈ 18` |

  So in that convention `δK₀ = 0` exactly, up to the realisation of the weight.

**Which reading holds is part of H17.** The question is whether the `|ĝ₃|²` light term enters the de-shifted functional
* unshifted, `X(ν)F(ν)`, or
* shifted, `χ[R₋F(ν+i) + R₊F(ν−i)]`.

The argument that settles this for the Einstein pole is H7 §5.2 E4: the `ν`-factor `γ_a²M²` is common to heavy and light states when
`k = j`. That argument does not apply to the `|ĝ₃|²` term, where `k = 2 ≠ j = 4`. So the reading is part of H17.

**Consequence.** The certified range of `K` at H17 = 0 is the interval between `K₀ + δK₀` and `K₀`. The table in §0 covers both ends.

---

## 3. The Einstein light value of the actual weight [certified]

**Definition.**

    E = 16Λ·V3,   V3 = ∫₀^∞ θχ_re[R₋D(ν+i) + R₊D(ν−i)]dν,   D(z) = 1/(z²+9/4).

**Comparison target.** The exact target `L_A/Λ = ∫ψD` is evaluated in closed form: an exact rational polynomial division by
`p² + 9/(4Λ²)`, plus `log` and `atan`, at 512 bits.

| `Δ_g` | `16L_A` (exact) | `E` | `E/16L_A − 1` |
|---|---|---|---|
| 300 | `10734.6054688713` | `10734.6054682539` | `−5.75·10⁻¹¹` |
| 1000 | `14339.5890638318` | `14339.5890757315` | `+8.30·10⁻¹⁰` |
| 2000 (grid) | `16427.1900682` | `16427.3300184` | `+8.52·10⁻⁶` (midpoint rule at the `ν = 0` kink) |
| 5000 (grid) | `19192.1199281` | `19192.2598783` | `+7.29·10⁻⁶` |

* The continuum values are within the a-priori bound `(4/9)ε_W`.
* The certified `τ₀` uses `E`, which is the actual light value. It therefore no longer depends on `ε_W`, which is only [DV].

---

## 4. H17

### 4.1 What it is
H17 is `K − K₀ − δK₀`: the AdS part of the two `|ĝ₃|²` light functions, beyond their flat monomials (H0 §15.2, §18; REVIEW_EE m4).
* **Its `k = 2` part.** `−Ĉ′_EE,ν|_{ĝ₃}` has Regge spin `n = 4 > k + 1`. So Lemma R (H7 §3) gives it as
  `(2πi/4)([𝒦^{R³}]_{r⁻¹} + c₄^{(2)} ⊛ [𝒦^{R³}]_{r⁻³})`. This needs:
  * the subleading Laurent coefficients `r⁻¹`, `r⁻³` of the AdS₄ `t₄`/`α_T` stress-tensor exchange, contracted with `EE`;
  * the `Π̂₄` admixture `c₄^{(2)}` of `Ĉ′_EE`. That admixture inherits the poles of `c(ν)⁻¹` at `±i/2`, `±3i/2` (REVIEW_EE m4).
* **Its `k = 3` part.** The same, with `r⁻²`, `r⁻⁴` and the `Π̂₅` admixture of `Ĉ₃`.
* **The `δK₀` convention.** §2.3: the choice of shift prescription for the `k = 2` term.
* **A possible split.** Separate `t₄²` and `α_T²` coefficients ("weighted norm") are a possible further part of H17.

### 4.2 Why it is not proved here
Proving any bound needs the exact light functions, and none of the project's tools gives them.
* **The shock/CHP route is blind.** C1_VALUE proved that the helicity-preserving `(2,4)` element of every `d = 3` CHP operator vanishes
  on `H₂` eigenfunctions [P]. So the `|ĝ₃|²` light term is carried by the energy-channel `T`-family term `I[W,W]`, not by a shock.
* **The factorisation route lacks two ingredients.** C1_FINAL's F₂ factorisation [D] reduces `I[W,W]` to products of the exact
  flip one-point matrices: `Ĥ_W ∝ e^{∓4iφ}` in the rest frame of `P`, and `M̂_W = 0`. Turning this into a function of `ν` needs:
  * the explicit celestial and `(2,4)` kernels of the CLS-D functional;
  * the `P`-measure, evaluated in full (not pointwise).

  These were never written down for the Einstein term either. The Einstein term was obtained from Lemma R plus H7's heavy/light
  argument, which is valid only at `k = j`.
* **No a-priori inequality is available.** The sum rule constrains `R` given `K`, not `K`. Positivity of the heavy side carries no
  information on `K`.

**Hence: no proven bound on H17, not even a weak one.** The `O(Λ⁻²)` scaling itself rests on the polynomial (Regge-counting) form,
which is [D-sketch] (H0 §15.2, Lemma R's R4 [D-sketch]).

### 4.3 What is proved: H17 is an exact linear form with certified coefficients
Suppose the light functions are the polynomials of H0 §18 (degree 3 and 2 in `ν²`):

    −Ĉ′_EE|_{ĝ₃} ∝ (ν/Λ)⁶[1 + a₄ν⁻² + a₂ν⁻⁴ + a₀ν⁻⁶],    −Ĉ₃|_{ĝ₃} ∝ −(ν/Λ)⁴[1 + b₂ν⁻² + b₀ν⁻⁴]    (h1 reading),

Then, exactly:

    H17 = −(2/Λ²)[a₄M₂ + a₂M₁/Λ² + a₀M₀/Λ⁴] + (2/Λ²)[b₂N₁ + b₀N₀/Λ²],

with `M_j = ∫Re X_θ p^{2j}dν/Λ` [certified] and `N_j = ∫ψ₃p^{2j}` [P, exact]:

| | `M₀` | `M₁` | `M₂` | `M₃` | `N₀` | `N₁` | `N₂` |
|---|---|---|---|---|---|---|---|
| 300 | `10.3859465067` | `1.52135733768` | `0.406104720474` | `0.142310277367` | `12.8355678175` | `2.34335451314` | `0.641926961159` |
| 1000 | `10.3859465067` | `1.51879115380` | `0.405293807716` | `0.141961111893` | same | same | same |

* **Leading order.** `Λ²H17 = −0.8106a₄ + 4.6867b₂ + O(Λ⁻²)`.
* **The model.** The heuristic `|a₄| = 3·36` and `|b₂| = 2·36` (`h3_H17_model.py`) reproduces `424.969398/Λ²` exactly, using `ψ₂`
  moments in place of `M₂`.
* **The weight of the two pieces.** The `Ĉ₃` coefficient `b₂` carries 79 % of the model band (`337.4` of `425.0`). Computing `b₂`
  alone would remove most of the uncertainty.
* **Tolerances for "`≤ 2.90·10⁻⁸`"** (leading order; the `a₂`, `a₀`, `b₀` terms are `O(Λ⁻⁴)`). The headline holds for
  * `|b₂| ≤ 1247` with `a₄` at its model size `108`, or
  * `|a₄| ≤ 6900` with `b₂` at its model size `72`.

  With the other coefficient zero, the limits are `|b₂| ≤ 1265` and `|a₄| ≤ 7316`.

### 4.4 The weighted-norm caveat
* **For `Ĉ₃` it is excluded [P given Q2, Q8, F].** H7 §4 route A gives the `T`-block light side of `Ĉ₃` as `∝ [H(n₁), H(n₃)]`, with
  `H = a𝟙 + b(e^{i(4φ+θ)}σ₊ + h.c.)`. This equals `2ib² sin(4φ₁₃)σ_z`. It depends on `b² ∝ t₄² + α_T²` and not on the phase `θ`.
  So the dominant (79 %) piece has no separate `t₄²` and `α_T²` coefficients.
* **For `k = 2` it is only [D-sketch].** Under F₂ the `W×W` trace is `2b²cos 4φ₃₁`, again without `θ`. But the first-moment matrix `M̂`
  of the parity-odd structure is known only on the diagonal (C1_FINAL §2.3(d)). So the statement is [D-sketch].
* **Status.** This is not a bound on H17. It narrows only the weighted-norm reading.

---

## 5. Combined result (`kc_bounds.py` → `cache/kc_bounds.txt`)

All entries are ball arithmetic; upper ends are rounded up.

| `Δ_g` | `τ₀(K=1)` | `τ₀(K_c)` | model band | robust rounding (H17 tolerance) |
|---|---|---|---|---|
| 300 | `3.13205501·10⁻⁶` | `3.13325626·10⁻⁶` | `[3.1246, 3.1408]·10⁻⁶` | `≤ 3.15·10⁻⁶` (`H17 ≥ −944/Λ²`, 2.2 × band); `≤ 3.2·10⁻⁶` (8.6 ×) |
| 1000 | `2.89126258·10⁻⁸` | `2.89136133·10⁻⁸` | `[2.8906, 2.8920]·10⁻⁸` | **`≤ 2.90·10⁻⁸`** (`H17 ≥ −5931/Λ²`, 14 ×); `|ĝ₃|Λ⁴ ≤ 120` (9.7 ×) |
| 2000 | `1.92831588·10⁻⁹` | `1.92833231·10⁻⁹` | `[1.9282, 1.9285]·10⁻⁹` | `≤ 1.93·10⁻⁹` (16 ×) |
| 5000 | `5.32617611·10⁻¹¹` | `5.32618337·10⁻¹¹` | `[5.3261, 5.3263]·10⁻¹¹` | `≤ 5.33·10⁻¹¹` (84 ×) |

**Quartic upper end.**
* Unchanged: `72.0406 / 92.8135` (2a) and `125.71 / 164.51` (2b).
* In (2b), `K` enters only through the sign condition `κ + εK/8 > 0` (A1PRIME §4, with `c₁ = 0`).

---

## 6. Draft replacement text (not applied)

**Supplement, Table S2, row L2.**

```latex
L2 & $|\hat g_3|^2$ light coefficient $K=K_0+\delta K_0+{\rm H17}$; $K_0=1+3.2\times10^{-11}$ (exact);
$\delta K_0:=-2\!\int(\mathrm{Re}X_\theta-\psi_2)p^6=-7.66632448\times10^{-4}$ / $-6.83014998\times10^{-5}$
($\Lambda^2\delta K_0=-68.308$ / $-68.097$; $0$ in the $\hat C'$-primary convention, the difference belonging to H17);
H17: AdS part of the $\hat C'_{EE}$, $\hat C_3$ light functions (subleading Regge data of the $R^3$ $T$ exchange,
$\hat\Pi_4$/$\hat\Pi_5$ admixtures), modelled $|{\rm H17}|\le424.97/\Lambda^2$
& D & $K_0$ [P]; $\delta K_0$ [certified]; H17 [D-sketch\,+\,M]; Eq.~\eqref{eq:taubound} needs only
${\rm H17}\ge-5.9\times10^{-3}$ at $1000$ for $\tau_0\le2.90\times10^{-8}$\\
```

**Supplement, Result statement, item (1)** (replaces "`K = 1 + O(Λ⁻²)` … At `K = 1`, `τ₀ = 3.13206·10⁻⁶` (300) and `2.89126·10⁻⁸` (1000)").

```latex
$K=K_0+\delta K_0+{\rm H17}$ is the $|\hat g_3|^2$ light coefficient, $K_0=1$ (to $3\times10^{-11}$),
$\delta K_0=-7.666\times10^{-4}$ / $-6.830\times10^{-5}$ [certified], H17 $=O(\Lambda^{-2})$ [D-sketch\,+\,M].
With ${\rm H17}=0$, $\tau_0=3.13326\times10^{-6}$ ($\Delta_g=300$) and $2.89136\times10^{-8}$ ($1000$); over the modelled
H17 band $\tau_0\in[3.1246,3.1408]\times10^{-6}$ / $[2.8906,2.8920]\times10^{-8}$; and
$\tau_0\le3.15\times10^{-6}$ / $2.90\times10^{-8}$ for every ${\rm H17}\ge-1.06\times10^{-2}$ / $-5.9\times10^{-3}$,
$2.2$ / $14$ times the modelled band (arithmetic certified, \texttt{code/kcert}).
```

**Supplement, the paragraph after Eq. `taubound`.**

```latex
($3.13326\times10^{-6}$ / $2.89136\times10^{-8}$ at $K=K_0+\delta K_0$, i.e.\ H17 $=0$; $3.13206\times10^{-6}$ /
$2.89126\times10^{-8}$ at $K=1$). Here $16L_A$ stands for the Einstein light value of the realised $k=2$ weight,
$16\Lambda\int\theta\chi[R_-D(\nu+i)+R_+D(\nu-i)]$, which equals $16L_A$ to $6\times10^{-11}$ / $8\times10^{-10}$ [certified].
```

**Supplement, the same paragraph, the H17 sentence** (replaces "modelled as a `±0.24 %` / `±0.021 %` band [M]").

```latex
(H17, modelled as $\pm424.97/\Lambda^2$ on $K$, i.e.\ $\pm0.24\%$ / $\pm0.021\%$ on $\tau_0$ [M]; $\Lambda^2{\rm H17}
=-0.811a_4+4.687b_2+O(\Lambda^{-2})$ in terms of the uncomputed subleading coefficients, the $\hat C_3$ one carrying 79\%)
```

**Supplement, I3.** Replace "`δK₀ = −68.08/Λ²` [D; DV]" by:

```latex
$\Lambda^2\delta K_0=-68.31$ / $-68.10$ ($300$ / $1000$; $-68.08$ / $-68.08$ at $2000$ / $5000$) [certified]
```

**Supplement, Table `tab:dg`.**
* Optionally add a column `τ₀(K_c)`: `3.1333·10⁻⁶`, `2.8914·10⁻⁸`, `1.9283·10⁻⁹`, `5.3262·10⁻¹¹`.
* The H17-band column stays.
* No `K_c` exists for 500; leave that entry blank or quote `K = 1`.

**Letter.** The numbers are unchanged. Replace "At `Δ_g = 1000` and `K = 1`, `τ₀ = 2.89·10⁻⁸`" by:

```latex
At $\Delta_g=1000$, $\tau_0=2.89\times10^{-8}$ ($K=1-6.8\times10^{-5}$ up to an uncomputed $O(\Lambda^{-2})$
remainder; $\tau_0\le2.90\times10^{-8}$ unless that remainder is below $-0.6\%$), that is $|\hat g_3|\Lambda^4\le120$~\cite{SM}.
```

The earlier sentence "`K = 1 + O(Λ⁻²)`" may stay.

---

## 7. What remains

1. **Compute H17.** This is the one open item. It needs:
   * the `r⁻¹…r⁻⁴` Regge Laurent coefficients of the AdS₄ `R³` `T` exchange contracted with `EE` / the odd structure;
   * the `Π̂₄`/`Π̂₅` admixtures of `Ĉ′_EE`, `Ĉ₃`.

   Equivalently, evaluate `I[W,W]` through F/F₂ with the explicit CLS-D kernels. **Start with `b₂` (`Ĉ₃`).**
   * It carries 79 % of the band.
   * Its factorisation F is the simplest: `f = −1`, no first moments, no de-shift.
   * It is already free of the weighted-norm issue (§4.4).

   Any computed `|a₄|, |b₂|` turns §4.3 into a certified `K_min` with no further work.
2. **Weights at other gaps.** A certified δK₀ at gaps other than 300, 1000, 2000 and 5000 needs weights there, which are not published.
3. **`ε_W` is still [DV].** It is no longer used for `K` or `E`: both are computed directly.

---

## 8. Reproduce

    cd ads4/code/kcert && ./kc_run_all.sh      # ~45 s; cache/kc_dk0_{300,1000}.txt, kc_grid.txt, kc_bounds.txt, kc_selftest.log (PASS)
