# LRC\* and LC\*: the uniform corner remainder of the de-shifted kernels, and what H_K now needs

*23 Nov 2026. Unreviewed, not committed. Answers the open item of LEMMAS_CLOSE §3.3–3.4 ("LRC uniform remainder [I]") and the
matching gap in H7 §14's corollary for `Ĉ′_EE`. No existing file was edited.*
* *New code:* `ads4/code/lrc/lrc_blocks.py` (module), `lrc1_identify.py`, `lrc2_deshift.py`, `lrc3_ee.py`.
* *Outputs:* `ads4/code/lrc/cache/lrc{1,2,3}_*.txt` (git-ignored).
* *Read-only imports and reuse:* `lemmas/lem3_k2nu_lrc.py` (part (b), re-run inside `lrc1`); the (x₀, z₀)-dependent factors and the EE
  structure are read off the header of `repro_engine/integrand.py`; the radial Casimir recursion is that of `reggemoment/blocks3d.py`
  (generalised to complex spin, and checked against it at integer spin to `1·10⁻¹⁵`).
* *Literature read from source (scratchpad copies):* KKSDZ 1904.05905 §4.2–4.4; CLS-D 2311.04271 §2–3.
* *Compute:* about 1 CPU-minute in total; one process at a time, one thread, < 1 GB.

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]**; **[M]**; **[Q]** quoted;
**[A]** assumption; **[I]** open. A chain carries the grade of its weakest link.

---

## 0. Verdict

1. **LRC\* holds for `k = 2`, and holds partially for `k = 4`.** The de-shift cancellation survives the O(ρ) corrections to all orders in ρ.
   * **What the exact kernel is [D + DV].** The exact (all orders in ρ) scalar `j′ = 0` per-ν kernel of CLS-D is
     `N(ν) · r² · Ĝ_{1,l}(ρ, ρ̄)`, with `l = −½ + iν`.
     * `Ĝ_{1,l}` is the `d = 3` Casimir eigenfunction of continuous spin `l` and dimension 1 in the (13)(42) channel, normalised to
       `P_l(η)` at `r → 0`.
     * Checked against a direct 2-d quadrature of the CLS-D integrand to `3·10⁻⁸` (§1). The leading-Regge form alone misses by up to 70 %;
       dimension 2 misses by up to 30 %.
   * **The ν-structure of the corrections [D].**
     * At each order `rⁿ` the correction is a finite sum of shifted harmonics `P_{l+m}(η)`, `|m| ≤ n`. In heavy language these are
       shifted copies `E_b(ν − ij′)`. The coefficients are **rational** in ν, not polynomial.
     * On the Regge box the dangerous (light-cone-edge) piece resums to all orders into one multiplier,
       **`F(ρ; ν) = ₂F₁(½, h; h+½; ρ²)`, with `h = (½ + iν)/2`** (§2). Its first term is `ρ²(1+2iν)/(2(3+2iν))`. There is no O(ρ¹) term:
       the O(ρ) shifts of the individual bases (LEMMAS §3.3) cancel in the kernel.
     * `F` is analytic for `0 ≤ Im ν < 3/2` and bounded uniformly for `ρ ≤ c < 1`. Its poles sit at `ν = i(3/2 + 2p)`, where the
       leading-Regge edge coefficient `c_l ∝ 1/Γ(l+1)` vanishes.
   * **`k = 2` (`Ĉ′_EE`).** The de-shift strip has depth 1 < 3/2, so the contour shift of HK Lemma S2 goes through unchanged. The exact
     cut-off kernel obeys `|K_c| ≤ C r^0 η^{−½} log η` on the whole box. **Criterion: `J₀ < 3/2`** (HK's full margin).
   * **`k = 4` (`Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`).** The strip has depth 3, so the pole of `F` at `ν = 3i/2` is crossed.
     * It leaves a residue `∝ f_pos(3i/2) · ρ² · η^{+1}`, where `f_pos` is the position-space edge coefficient of the weight.
     * The residue vanishes iff `f_pos(3i/2) = 0`. That happens automatically when the position weight is regular at `3i/2`.
     * Residual growth `c = 1`, hence **criterion `J₀ < 2` in general** (HK claimed 5/2), and `J₀ < 5/2` if `f_pos(±3i/2) = 0`.
   * **Numbers [DV]** (§3), exact kernels at `ρ = 0.2, 0.5, 0.85` and `β ≤ 16`:
     * the local growth exponent of (exact − leading Regge) is **−0.50** for the entire `k = 2` weight and for the regular `k = 4`
       weight, and tends to −½ with the expected `β η^{−½}` form (−0.38 at `β ≈ 13`) for the `k = 2` weight with target poles at `±i`;
     * it is **+1.00** for the `k = 4` weight with a pole at `3i/2`, and there the residue formula matches to `10⁻⁶`.
   * **EE structure [D-sketch]** (§4). The EE dressing of both `Ĉ′_EE` and `Ĉ₄^{ds}` is ν-independent and contains no `z₀` [P + DV].
     * So the (13) side, which carries the light-ray pair, is the scalar one, and the EE kernel is a (13) Casimir eigenfunction of the same
       representation `(1, l)`.
     * Its edge functions are collinear blocks with `a = 0` (1↔3 symmetry) and integer-shifted `b`. These have the **same pole set**
       `Im ν = 3/2 + 2p` [P + DV].
     * The step "spinning partial wave = scalar partial waves with shifted (24) labels" is standard (CPPR; weight-shifting) but was not
       checked numerically here.
2. **LC\* is not needed, and is not implied by standard results.**
   * At fixed ρ, KKSDZ's *asymptotic light-cone expansion* gives `κ = τ_*/2 ≥ ½`. That expansion is an unproved conjecture of
     Hartman–Jain–Kundu, argued only for identical scalar pairs.
   * Uniformity in the Regge corner is in neither KKSDZ nor CLS-D. KKSDZ bound the corner only by Rindler/Cauchy–Schwarz, with edge
     exponent 0: that is K1′.
   * Grade of LC\*: **[A]** (§5).
3. **H_K.** **[D-sketch]**, now with **no [I] link, no `J₀ < ½` and no LC\*.**
   * Explicitly: **[D given K1′, K2ν, K3, DS, LRC-a, LRC\*]**, with:
     * LRC\* [D + DV] for the scalar-type kernel and [D-sketch] for the EE structure;
     * K1′ and K2ν [D];
     * DS [D + NE] (D) and [NE] (Q).
   * Criteria actually achieved:

     | piece | criterion | margin under `J₀ ≤ 1` |
     |---|---|---|
     | `Ĉ′_EE` (`k = 2`) | `J₀ < 3/2` | ½ |
     | `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` (`k = 4`) | `J₀ < 2` (5/2 if `f_pos(±3i/2) = 0`) | 1 |
     | `Ĉ₅^{ds}` (raw), `Ĉ₃` | `J₀ < 3/2`, unchanged | ½ |
4. **H7 §14's corollary (Result D, `Ĉ′_EE`).** It holds for the **exact** kernel of the cut-off functional, with margin ½ (`J₀ < 3/2`).
   * Grade: **[D-sketch]**, through the EE step.
   * For the scalar-type kernel the same statement is **[D + DV]**.
   * As in HK, it holds for `Ω_c`, not for compactly supported grid weights.
5. **Letter sentence:** see §6.3.

---

## 1. The exact per-ν kernel is a continuous-spin Casimir eigenfunction (`lrc1_identify.py`)

### 1.1 Statement and derivation

**The integrand.** In CLS-D's integrand (eq. `scalarfunctional_definition` / `spinningfunctional_definition`; the project's
`repro_engine/integrand.py`), ν enters **only** through `A^{iν/2}`, where `A = Nn²c₁₃/(m₂₄m₂₀m₄₀c₁₀c₃₀)`. In the H7 §18.1 frame:
* `x₁ = −y`, `x₃ = y` and `x₂ = −x₄ = x/x²`;
* `x₀ = |y|n̂_θ` is exact: it is the circle cut out by the two light-ray δ-functions;
* `z₀` runs over the celestial circle.

So the f-stripped scalar kernel is the 2-d integral

    K_ν(r, β) = ∫dθ ∫dφ  A^{iν/2} |Nn|^{−½} m₂₄^{¾} (m₂₀m₄₀)^{−¼} c₁₃^{5/4} (c₁₀c₃₀)^{−¼}.

**The dual side.** With `Nn ∝ V_{0,42}X₂₄`, the (2,4) factor is exactly `V_{0,42}^{l}(X₂₄/(X₂₀X₄₀))^{(1+l)/2}`, with `l = −½ + iν`.
* This is the three-point structure of an operator of dimension 1 and spin `l` with two weight-0 partners.
* This agrees with CLS-D's dual structure `⟨0|φ₄L[O]φ₂|0⟩^{−1}` with `O = (d/2 + iν, J = −1)`. Indeed, `L[O]` carries `(2, −½−iν)`, and
  its dual under the Lorentzian pairing carries the shadow `(1, −½+iν)`.
* The (1,3) celestial factor is its partner.

**Consequence [D].** The Casimir of the (13) pair acting on the light-ray side returns `C = 1·(1−3) + l(l+1) = −9/4 − ν²`. So `K_ν` is a
(13)(42)-channel Casimir eigenfunction.
* Its leading Regge form is `r²P_l(η)`.
* Hence `K_ν = N(ν) · r^p · Ĝ_{Δ′,l}`, where `Ĝ` is generated from `P_l` by the radial Casimir recursion and `Δ′ ∈ {1, 2}`.

### 1.2 Check [DV] (`cache/lrc1_identify.txt`)

**The test.** It is prefactor-free: `Q(x) = [K_ν(x)/K_{ν₀}(x)] / [Ĝ_ν(x)/Ĝ_{ν₀}(x)]` must be independent of `x = (ρ, β)`.
* Grid: 9 points, `ρ ∈ {0.25, 0.55, 0.85}`, `β ∈ {0.6, 1.4, 2.2}`.
* `ν₀ = 0.8`.
* Quadrature: `θ` by trapezoid; `φ` by a logistic map at the two celestial singularities. Its own accuracy is `5·10⁻¹⁰`.

| ν | spread of `Q`, `Δ′ = 1` | `Δ′ = 2` | leading Regge only (`Ĝ → P_l`) |
|---|---|---|---|
| 0.4 | **2.6·10⁻⁸** | 6.8·10⁻² | 0.13 |
| 1.3 | **1.3·10⁻⁷** | 0.30 | 0.67 |
| 2.7 | **1.5·10⁻⁶** (quadrature) | 0.25 | 0.74 |

**The prefactor.** `K_{ν₀}/Ĝ_{ν₀}` over the 9 points:
* it has constant phase, to `5·10⁻⁹`;
* `log|·| = 2.000000 log r + const`, with residual `3·10⁻⁸`.

So the exact f-stripped kernel is `N(ν) r² Ĝ_{1,l}`, with nothing else depending on `x`.

**Reproduction.** `lrc1` first re-runs LEMMAS' part (b): the individual bases `X₂₀`, `z₀·x₂₀` deviate by exactly `ρ`. The deviation is
real, but it organises into a Casimir eigenfunction.

---

## 2. The ν-structure of the corrections (`lrc_blocks.py`, `lrc2` (a)–(b))

**Radial expansion [D].**
* With `K(E, j) = E(E−3) + j(j+1) − C`, the recursion of `blocks3d.py` holds verbatim for complex `j`. Every Legendre re-expansion
  coefficient is rational in `j`, with no truncation.
* It gives

      Ĝ_{1,l}(r, η) = Σ_{n even} rⁿ Σ_{m=−n, step 2}^{n} B_{n,m}(ν) P_{l+m}(η),     B_{0,0} = 1.

* So **at each order `rⁿ` the correction is a finite sum of shifted harmonics with rational coefficients.** In heavy language it is
  `Σ_{j′}(rational)·E_b(ν − ij′)P_{−½+i(ν−ij′)}`, with `|j′| ≤ n + n_f`. This answers the question as posed: shifted copies, yes, but
  with rational and not polynomial coefficients.

**Poles of `B_{n,m}` [P].**
* They are the zeros of `K(1+n, l+m) = n(n−1) + m(m + 2iν)`, i.e. `ν = i[n(n−1) + m²]/(2m)`.
* The smallest `|Im ν|` are `1.5, 3.5, 4, 5.5, …`. The top shell `m = n` gives `ν = i(n − ½)`.
* The Legendre denominators `(2j+1)(2j±1)` add poles at `ν ∈ iℤ`. These cancel between `±m` (`Ĝ` is real on the real axis: `|Im Ĝ|/|Ĝ| ≤ 4·10⁻¹⁴`).

**Organisation in ρ: the edge shells [D].** On the box, `rⁿP_{l+m} ∼ ρ^{(n+m)/2}ρ̄^{(n−m)/2}`. So the expansion reorganises into edge
shells `ρ̄^{h̄+s}g_s(ρ)`:
* **The top shell `s = 0` (`m = n`) carries the raw growth.**
* At leading order in `z̄` the Casimir equation reduces to the SL(2) one, `D_z g₀ = h(h−1)g₀`, with `h = (Δ′+l)/2`. The `r → 0`
  boundary condition selects
  `g₀ = c_l z^h₂F₁(h, h; 2h; z) = c_l(4ρ)^h ₂F₁(½, h; h+½; ρ²)`.
* Hence

      exact edge piece = (leading-Regge edge piece c_l e^{lβ}) × F(ρ; ν),     F = ₂F₁(½, h; h+½; ρ²),

  to all orders in ρ, with corrections `O(ρ̄/ρ) = O(e^{−2β})`.
* **Check** (`lrc2` (a); ν = 0.7, 3.1; ρ = 0.3, 0.85): `|Ĝ − 2Re[c_l e^{lβ}F]|/|Ĝ|` falls as `e^{−2β}`. It equals
  `1.0·10⁻² … 3.6·10⁻⁸` at `β = 2 … 8`, while `|F − 1|` ranges over 0.025–0.67.

**Analytic structure of `F` [P].**
* Its poles are at `h + ½ = −p`, i.e. **`ν = i(3/2 + 2p)`**.
* There `l = −2 − 2p`, so `c_l = Γ(l+½)/(√πΓ(l+1)) = 0`: the pole of `F` sits on a zero of the leading-Regge edge coefficient.
* On `0 ≤ Im ν ≤ 1`, `Re(h + ½) ≥ ¼`, so `|F| ≤ C_c` uniformly in ν for `ρ ≤ c < 1`. At large `|ν|`, `F → (1−ρ²)^{−½}`.

**Subleading shells `s ≥ 1` [D-sketch + DV].** They carry an extra `(ρ̄/ρ)^s = e^{−2sβ}` relative to the top shell. Their coefficient poles
(`m < n`) lie at `|Im ν| ≥ 4` and are irrelevant to raw counting.

---

## 3. Survival of the de-shift cancellation (`lrc2` (c)–(d))

### 3.1 Argument [D for the scalar-type kernel]

**Setting.**
* Let `f_pos(ν)` be the position-space leading-Regge edge coefficient of the cut-off weight: `w_c` times the top-harmonic factor.
* By HK Lemma S2/S4, transferred to position space as in LEMMAS §3.2 (`E_b s` has no zeros or poles in the strip), `f_pos` is analytic:
  * in `0 < Im ν < 1` for `k = 2`, with its target singularity on `Im ν = 1`;
  * to depth 2 for `k = 4` [D], depth 3 given HK's NE.
* The exact edge piece of the smeared kernel is `ηⁿ ∫_ℝ f_pos(ν)F(ρ; ν)e^{lβ}dν`.

**The contour shift.**
* **`k = 2`.** Shift to `Im ν = 1`. `F` is analytic and bounded there, so the result is the leading-Regge target multiplied by the bounded
  factor `F(ρ; ν_target)`, and is `O(η^{−½}log η)` uniformly in `ρ ≤ c`. The shells `s ≥ 1` are `O(η^{1−½−2})` raw.
* **`k = 4`.** The shift passes `ν = 3i/2`.
  * The residue is `2πi f_pos(3i/2)·Res F·η³e^{−2β} ∝ ρ²η^{+1}`. It is zero iff `f_pos(3i/2) = 0`.
  * Beyond it, the shift continues to depth 2 (or 3), where the next pole of `F` lies at `7i/2`.
  * The shell `s = 1` is raw `O(η^{½})`.
* **The cutoff.** `c(ν)` is entire and bounded on strips, so Lemma C applies unchanged.
* **The measure and f-factors** are ν-independent and bounded relative to their leading-Regge forms [P, LEMMAS §3.3 (c)]. The kernel
  prefactor is exactly `r²` [DV].

**Uniform corner bounds (Lemma LRC\*).** On the box `B_c`:
* `|K_c^{k=2}| ≤ C_c η^{−½}log η`;
* `|K_c^{k=4}| ≤ C_c r²(ρ²η·|f_pos(3i/2)| + η^{½})`.

**Corner count (HK §1: `c < a + 1 − J₀`).**
* `k = 2` (`a = 0`, `c = −½`): **`J₀ < 3/2`**.
* `k = 4` (`a = 2`, `c = 1`): **`J₀ < 2`**; with `f_pos(±3i/2) = 0`, `c = ½` and `J₀ < 5/2`.
* HK's `7/2` (NE) would need in addition the shell `s = 1` de-shifted. It is not claimed.

### 3.2 Numbers [DV] (`cache/lrc2_deshift.txt`)

**Set-up.**
* `S = ηⁿ∫₀^∞ w Ĝ_{1,l}dν`, against `S_LR = ηⁿ∫w P_l`.
* Exact kernel: radial order 260 (converged at `ρ ≤ 0.85`).
* Quadrature: 700 Gauss–Legendre nodes on `[0, 40]`.

The table gives the local exponent `d log|S − S_LR|/dβ` over `β ≈ 8.5–13.5`, at `ρ = 0.2 / 0.5 / 0.85`.

| weight | type | exponent of (exact − LR) | LR target |
|---|---|---|---|
| `w2g = ν²e^{−ν²/9}` | `k = 2`, entire | **−0.50 / −0.50 / −0.50** | −0.50 |
| `w2s = ν² sech(πν/2)` | `k = 2`, poles at `±i` (like D's weights) | −0.33 → −0.38 (≈ `βη^{−½}`) | same form |
| `w4g = ν²(ν²+1)(ν²+4)e^{−ν²/9}` | `k = 4`, regular at `3i/2` | **−0.50** (β ≤ 9; double-precision floor beyond) | −0.50 |
| `w4p = w4g/(ν² + 9/4)` | `k = 4`, pole at `±3i/2` (`f_pos(3i/2) ≠ 0`) | **+1.00 / +1.00 / +1.00** | −0.50 |

* **Raw rates** are +0.5 (`k = 2`) and +2.5 (`k = 4`). The exact corrections never show them.
* **`w4p` residue formula.** `η³Re[2πi Res_{3i/2} w c_l e^{lβ}F]` reproduces `S − S_LR` to `1.000003` (`ρ = 0.2`), `1.000002` (0.5) and
  `1.000000` (0.85) at `β = 10`. The same residue for `w4g` is `≤ 2·10⁻²⁷`.
* **Uniformity.** The amplitude of the exact correction grows with ρ: about `ρ²`, then the full `F`. Its rate is independent of ρ up to 0.85.

---

## 4. The EE kernels of `Ĉ′_EE` and `Ĉ₄^{ds}` (`lrc3_ee.py`)

**The dressing [P + DV].** Both functionals use `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²` together with `V_{0,42}^{−4}`. So the EE
kernel is the scalar `(θ, z₀)` integrand times

    D = (H₂₄ − 2V_{2,04}V_{4,02})²/(X₂₀X₄₀)².

* `D` is ν-independent and **contains no `z₀`**: its variation over `z₀` is `≤ 1.1·10⁻¹⁵`, for random null polarisations normalised in the
  frame of `x̂`.
* Its `x₀`-dependence is **bounded but not small**. `max|D − D_LR|/max|D_LR|` is up to 12 (ρ = 0.3), 255 (ρ = 0.5) and `1.2·10⁴` (ρ = 0.8),
  consistent with `(1−ρ)^{−8}`.

**Argument [D-sketch].**
1. The (1,3) side, with the light-ray pair and the celestial `j′ = 0` structure, is the scalar one (H0: "the celestial factor is the scalar
   one"). So the (13) Casimir still returns `C = −9/4 − ν²` [D].
2. The 1↔3 symmetry of the light-ray pair gives `a = 0` in every collinear block [D].
3. The (24) spin enters through ν-independent covariant factors. By CPPR's reduction of spinning partial waves (or weight-shifting), these
   amount to scalar partial waves of the same representation `(1, l)` with integer-shifted (24) labels `b` [Q-standard, not checked here].
4. The edge coefficients are then `(h)_p(h+b)_p/((2h)_p p!)`. With `a = 0` the factor `(h)_p` cancels every even-`2h` pole. **So the only
   poles in `0 < Im ν ≤ 3` sit at `Im ν = 3/2`, for every `b`** [P + DV].

| `(a, b)` | poles in `0 < Im ν ≤ 3` |
|---|---|
| `(0, 0)`, `(0, ±1)`, `(0, ±2)` | 3/2 |
| `(0, ½)` | none |
| control `(½, ½)` | ½, 5/2 |
| control `(−1, 1)` | ½, 3/2 |

**Consequences.** The `k = 2` EE kernel keeps the full de-shift cancellation (`J₀ < 3/2`). The `k = 4` EE kernel has the same single residue
at `3i/2` (`J₀ < 2`). EE3's own leading-Regge `R_j` have poles at `±3i/2` (LEMMAS §3.2). If they coincide with this residue, the bound
acquires at most a factor `log η`, which does not change the strict criterion.

**Precise remaining gap (EE only).** An end-to-end check that the EE per-ν kernel equals the predicted spinning partial wave (the analogue
of §1.2) was not done. It needs `d = 3` spinning blocks with continuous spin, or the 2-d quadrature with the EE dressing at `β ≳ 6`
(features of width `e^{−β}`). This is a check of a standard reduction, not an open physics question.

---

## 5. LC\*: what the literature gives

**KKSDZ (1904.05905).**
* §4.2: the Rindler/Cauchy–Schwarz argument bounds the second sheet by first-sheet correlators, uniformly for `r < r₀` and all `σ`. That
  is the uniform corner bound with edge exponent 0, i.e. K1′.
* For the light-cone limits they invoke the **asymptotic light-cone expansion** (their eq. `lightconeasymptotic`): the OPE is valid
  asymptotically on the second sheet. They write: "a natural expectation, … we don't have a general proof that it holds".
  * It is due to Hartman–Jain–Kundu (1509.00014), who argued it only for identical scalar pairs with real coordinates, using positivity of
    the t-channel sum.
* With it, the double commutator decays in the light-cone limits at the rate of the first twist accumulation point `τ_*`, because finite
  sums of s-channel blocks have zero dDisc. In our variables this gives `κ = τ_*/2 ≥ (d−2)/2 = ½ > ¼`.
* For the Regge corner: "For the Regge limit, we don't have an analogue of [the asymptotic expansion]". The corner is controlled only by
  `J₀^{dDisc} ≤ 1` from Rindler positivity.

**CLS-D (2311.04271)** quotes only `J₀ ≤ 1` (citing HJK, Caron-Huot and KKSDZ) and fixes the Regge convergence by subtraction factors
(`J₀ < k`). It makes no light-cone uniformity statement.

**Verdict.**
* LC\* at fixed ρ is KKSDZ/HJK's asymptotic light-cone assumption: **[A]**, standard but unproved for spinning operators.
* LC\* uniform over the Regge box ("Regge and light-cone limits commute with a power bound") is in neither paper and does not follow from
  them: the asymptotic expansion is pointwise in the other variable. It is **[A]**, beyond the literature.
* With LRC\* in hand it is **not needed**. If it were used instead, H_K would be [D given A].

---

## 6. Final statements

### 6.1 H_K (replaces LEMMAS_CLOSE §3.5's H_K row)

> **H_K.** Take each cut-off de-shifted piece `Ĉ_c ∈ {Ĉ′_EE, Ĉ₄^{ds}, ∂Ĉ₄^{ds}, Ĉ₅^{ds}}` (HK §4: `χ → χ·c`), and `Ĉ₃`. Assume:
> * K1′ (uniform Rindler bound on the Regge box, at `J₀ = 1`) [D];
> * K2ν [D];
> * K3 / Lemma S [D];
> * DS: [D + NE] for D, [NE] for Q.
>
> Then `Σ_O Tr[P_O Ĉ_c[G_O]] = −Ĉ_c[light]`, with absolute convergence.
> * The exact corner kernels satisfy `|K_c| ≤ C r^a η^{c}` with `c = −½ + 0` for `k = 2` and `c = 1` for `k = 4` (Lemma LRC\*, §3.1).
> * The heavy sum therefore converges for `J₀ < 3/2` (`k = 2`, `k = 5`, `Ĉ₃`) and `J₀ < 2` (`k = 4`; `5/2` if the weight's
>   position-space edge coefficient vanishes at `ν = ±3i/2`).
> * This holds in particular under `J₀ ≤ 1`.
>
> **Grade: [D-sketch].**
> * The weakest links are LRC\* for the EE structure (§4) and LRC-a for EE at `k = 4` [D-sketch, LEMMAS]. DS for Q is [NE].
> * For the scalar-type kernel the chain is [D + DV].
> * **No [I] link remains. `J₀ < ½` (H_ip) is not needed, and neither is LC\*.**

### 6.2 H7 §14 corollary, Result D

> For `Ĉ′_EE` with the analytic cutoff, the large-η corner integral `∫|K_c||dDisc 𝒢|` of the **exact** kernel is finite under K1′ with
> margin `η^{−½}` (`J₀ < 3/2`).
> * The leading-Regge statement of H7 §14 extends to all orders in ρ, because the corrections multiply the edge coefficient by `F(ρ; ν)`,
>   which is analytic in the de-shift strip.
> * Grade: **[D-sketch]** (EE step); **[D + DV]** for the scalar-type kernel.
> * The corollary applies to `Ω_c`, not to the stored grid weights (HK §3).

### 6.3 Sentence for the letter

> "For the analytically cut-off de-shifted functionals the heavy OPE sum may be interchanged with the Regge-corner and ν integrations
> by Tonelli's theorem, using only the Rindler bound `J₀ ≤ 1` made uniform over the Regge box. The exact per-ν kernels are
> continuous-spin conformal partial waves. At the light-cone edge, their corrections beyond leading Regge order resum into a
> hypergeometric factor with no singularity inside the de-shift strip, so the de-shift cancellation persists to all orders. The heavy
> sum then converges absolutely for `J₀ < 3/2` (cubic functional) and `J₀ < 2` (quartic functionals). We derive this for scalar-type
> kernels and verify it numerically; for the spin-2 structures it rests on the standard reduction of spinning partial waves to scalar
> ones."

### 6.4 Proposed table rows (UNIFIED_THEOREM §3.1, LEMMAS_CLOSE §5; not edited here)

| # | before | proposed |
|---|---|---|
| LRC uniform remainder (LEMMAS §3.3) | [I] | **LRC\*: [D + DV]** for the scalar-type kernel; **[D-sketch]** for EE. `k = 2` full; `k = 4` one residue at `3i/2` (`c = 1`) |
| exact per-ν kernel | — | `N(ν) r² Ĝ_{1,−½+iν}` (continuous-spin Casimir eigenfunction) **[D + DV, 3·10⁻⁸]** |
| LC\* | [I] | **[A]**; not needed. KKSDZ/HJK asymptotic light-cone expansion (fixed ρ) [A, standard]; uniform version beyond the literature |
| H_K | [D-sketch] with [I] link | **[D-sketch]**, no [I] link; criteria `J₀ < 3/2` (`k = 2`), `J₀ < 2` (`k = 4`) |
| HK §0 / §4.5 "`J₀ < 5/2` for `k = 4`" | [D] | **correct to `J₀ < 2`** for generic weights (5/2 only if `f_pos(±3i/2) = 0`) |
| H7 §14 corollary (`Ĉ′_EE`) | [P given K1′] (leading Regge only) | exact kernel: **[D-sketch]** (EE), [D + DV] (scalar-type), margin ½ |

### 6.5 What remains

1. **EE end-to-end check** (§4): the identity of §1.2 with the EE dressing, or spinning continuous-spin blocks. This would raise LRC\* for
   EE from [D-sketch] to [D + DV].
2. **`f_pos(±3i/2)` for Q's EE3 weights.** If it vanishes, the `k = 4` criterion returns to `J₀ < 5/2`. No other change follows; this is
   not needed under `J₀ ≤ 1`.
3. **Unchanged from LEMMAS §4:** the borderline R⁴ arc at `a = b` (a light-side question in the same corner) stays [I]. This note controls
   the kernel side only.

## 7. Reproduce

    cd ads4/code/lrc
    ulimit -v 3000000
    OMP_NUM_THREADS=1 python3 lrc1_identify.py   # ~10 s -> cache/lrc1_identify.txt
    OMP_NUM_THREADS=1 python3 lrc2_deshift.py    # ~45 s -> cache/lrc2_deshift.txt
    OMP_NUM_THREADS=1 python3 lrc3_ee.py         # ~2 s  -> cache/lrc3_ee.txt

| file | content |
|---|---|
| `lrc_blocks.py` | continuous-spin `d = 3` radial Casimir recursion; scaled Legendre functions of complex degree; edge factor `F`; `c_μ` |
| `lrc1_identify.py` | LEMMAS §3.3 (b) re-run; exact CLS-D scalar kernel by 2-d quadrature vs `Ĝ_{Δ′,l}` (`Δ′ = 1, 2`) and vs leading Regge; prefactor `r²` |
| `lrc2_deshift.py` | edge structure; pole table; smeared exact vs leading-Regge kernels for `k = 2, 4` weights; growth exponents; residue prediction at `3i/2` |
| `lrc3_ee.py` | EE dressing is ν- and `z₀`-free and bounded; pole content of collinear coefficients for `a = 0` vs controls |
