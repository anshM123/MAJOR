# Independent check 2: the de-shift identity, the monopole cost, Lemma LR∞, the order of limits, Theorem W, and the ⟨TTO⟩ counts

*23 Nov 2026, independent checker. Asked for by REFEREE_REPORT_3 §M6. No existing file was edited and nothing was committed.*
* *Code:* `ads4/code/indep2/` (new; `run_all.sh` reproduces everything in ≈ 1.5 CPU-h). Outputs are in `ads4/code/indep2/cache/`; its
  `.gitignore` ignores everything except itself.
* *Independence.*
  * All numerical code is new. Nothing from `ads4/code/` outside `indep2/` was imported.
  * The (0, b) block recursion re-implements the operator written in the docstring of the previous independent check
    (`code/indep/cblocks.py`, INDEP_CHECK.md). It was re-coded so that it works over any field: sympy with l and b symbolic, or mpmath
    with complex l.
  * Inputs:
    * `prl/data/` only: `kernels.json`, `d_abel_X_{1000,300}.txt`, `nf02_psi.json`;
    * the statements in PW_IMAGES, LARGESPIN_L, LARGESPIN_W, REGGE_TREE and P2_SPINNING.
  * **Two disclosed exceptions.**
    1. Target 4 needs the Route W floors f_h (certified flat NF02 action + curvature). These are not in `prl/data/`, so they were read as
       data from `code/lsw/cache/costs_1000.json`, together with Route W's cached costs (used for comparison only). Every cost entering
       my J\* was recomputed with my own code.
    2. Target 2: my evaluation of the stated remainder bound came out 36× tighter than the stated r(β) table. To diagnose this I read
       `code/lsl/lsl_lrinf.py` (the bound routine, lines 38–107) *after* my result was fixed (§2.3).
* *Literature.* arXiv sources and PDFs of 2106.10274 (CHMRS), 2311.04271 v3 (CLS-D), 1904.05905 (KKSDZ) and 1708.05718 (DKKPS) were
  downloaded and read at the cited passages. Equation numbers were read off the PDFs.
* *Compute.* ≈ 1.5 CPU-h; at most 2 processes; `ulimit -v 3000000`; `OMP_NUM_THREADS=1`.

Grades as in the project: [P], [D], [D-sketch], [DV], [NE], [M], [A], [Q], [I].

---

## 0. Summary

| # | claim (source) | independent result | agrees? | grade after check |
|---|---|---|---|---|
| **1a** | **Magic-block identity** (PW_IMAGES (2.1)): for each b-shifted block, 2η(1+r²)Ĝ¹_b(ν) = (1+i/2ν)Ĝ⁰_b(ν+i) + (1−i/2ν)Ĝ⁰_b(ν−i), with b-independent coefficients. Claimed [P] to O(r⁸) for b = 0…4, [DV] to r⁴⁰ | **Exact (sympy) with l and b both symbolic to O(r⁶)**: 0 nonzero residuals out of 35 coefficients. So the identity is a polynomial identity in b at each order, not just a fact about b = 0…4. **Exact to O(r⁸), l symbolic, every integer b = −4…4**: 0 nonzero out of 30 (b = 0) / 54 (b ≠ 0), the same counts as the original. **50 digits, n ≤ 40**, four ν (two complex), b = −4…4, 0.37, 7, 10: residual ≤ 4·10⁻⁴⁸. Pointwise with directly evaluated P_{l+m}(η) of complex degree: the residual is pure truncation, reaching 10⁻³⁰ at N = 60. Controls (swapped coefficients, Δ shifted on either side) fail at O(1) | **yes** | **[P] to O(r⁶) for all b, [P] to O(r⁸) for b ∈ {−4…4}, [DV] to r⁴⁰; independently reproduced.** All orders: [D] (unchanged) |
| 1a′ | PW2 side claims: f = v′−u′ = η(1+r²)/(2r) in the H7 frame; Ĝ⁰_b even in μ; analytic at μ = ±i/2, ±3i/2 | Frame identity **exact** (sympy, standard cross-ratios). Evenness to **9·10⁻⁶⁰** (b = −4…4, N ≤ 14). Removability: max\|B\| is the same at ε = 10⁻¹², 10⁻²⁴, 10⁻³⁶ from the would-be poles (a pole would scale like 1/ε). The algebra (2.1) ⇒ (2.2) checks | **yes** | [P] (frame), [DV] (evenness, removability) |
| **1b** | **Monopole cost** (PW_IMAGES (4.1)): cost(𝒜_g[F]) = ∫ψ₂F/c + D·F(i/2), D = −9.98332·10⁻⁴ (Δ_g = 1000), −1.09754·10⁻² (300); D independent of t | **D(1000) = −9.9833211·10⁻⁴, D(300) = −1.0975445·10⁻²** (D·Λ² = −995.3 / −977.9), identical in every printed digit at t = 0 … 6. **New exact mechanism:** 1/c = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9), so **D = (15/2)·defect[1/(4ν²+1)]**, and t-independence follows because (P_t − 1)/(4ν²+1) is strip-analytic. Measured defect of that function: ≤ 10⁻¹⁸. Formula (4.1) holds on 8 random mixtures with ν² factors: 2·10⁻¹⁵ … 5·10⁻¹⁰. The Ω̃(t) table is reproduced entry by entry | **yes**, digit for digit | **[DV], independently reproduced**; t-independence and the form of (4.1) upgraded to **[D]** (partial fractions + de-shift identity); value of D [DV] |
| **2** | **Lemma LR∞** (LARGESPIN_L §3): Ω_k = a₁M₁/Λ² − a₂/(2b̂³) + E_k; for β ≥ 1, u ≤ 1, σ = ±1: F_LR ≥ (1−r)A₁M₁/Λ² + 2047.1/(2b̂³) > 0, with the r(β) table | The two-term form is accurate to **\|E₂\|/main ≤ 3·10⁻⁷** at every β ∈ [1, 32]: high-precision Mehler–Dirichlet with exact Φ, cross-checked against a direct p-quadrature to 12 digits. **The true r(β) is ≤ 8·10⁻⁷**, against the stated bound of 0.577 … 1.6·10⁻⁷. F_LR,min = 2Ω₂ − \|Ω₃\| exceeds the stated bound at all 19 β checked (1 ≤ β ≤ 32). The Hadamard identity of step 3 holds to 15 digits; a₁, a₂, A₁ = 327.149, −A₂ = 2047.09 and the J\*(m) table are all reproduced. **The r-table is not reproduced from the formula as written:** with R = min(R_∞, a₁/x² + ∫\|ψ\|) the stated bound gives r(1) = 0.016, not 0.577. The table was computed with the cruder Φ_max = Σ\|c_n\|/(n+1) ≥ ∫\|ψ\| (diagnosed in the code). Both are valid bounds, so this is a documentation slip in the conservative direction | **yes** (lemma); table: conservative, its formula mis-stated | **[P + DV], independently reproduced.** Side result: with ∫\|ψ\| the bound stays < 1 down to β = 0.5 (r = 0.28), i.e. J\* = 324 at the gap instead of 857 (needs certification) |
| **3** | **Order of limits, route N** (REGGE_TREE §0, §4) and its literature quotes | **Logic sound.** It is CHMRS §1.1.3, eqs. (1.12)–(1.14), verbatim: the exact sum rule at finite ε, heavy side used only through its sign, only the light side expanded. Every quote in REGGE_TREE §2 was checked against the source and is verbatim or faithful. Equation numbers (1.12)–(1.14), (2.8), (2.39), (4.45), (4.55), (4.70), (4.73), (4.79)–(4.83) are correct. Minor inaccuracies: 3 section labels; CLS-D is cited for a "logic" it only restates; f⁻ = −1, not 1. The burden moves to the uniformity of the light side (A2 [A]); REGGE_TREE says so | **yes** | [D given A2 restated], as REGGE_TREE says; literature [Q] **verified** |
| **4** | **Theorem W's implication** (LARGESPIN_W §4.0) and the J\* table (§4.3; P2_SPINNING §7.5) | 2×2 criterion d ≤ H = f₁f₂/(f₁+f₂): correct and sharp (0/200 000 random violations; indefinite at d = 1.001H). The conclusion Q ≥ (1−θ)min f\|λ\|²/2 holds with a factor 2 to spare. Gershgorin d ≤ min f/3: correct. "Lowers tolerances by ≤ 1.5": true only when f₊₊ ≈ f₊₋ and f_OO ≥ min(f₊₊, f₊₋) (sampled floors give 1.50–1.53; it can reach 3 otherwise). **Own Lemma-I costs** agree with Route W's to ≤ 1·10⁻³ (A, B, C1, C3), 2·10⁻² (C2), 1·10⁻² (sh); stripped costs agree to ≤ 10⁻³. **All 6 R1 rows, the R2 row and the 4 P2 §7.5 entries of the J\* tables are reproduced exactly** (30/38, 50/100, 64/none(300, 344); R2: none(300, 344); model: 38/50) | **yes**, with two qualifications (§4.4) | implication [P] confirmed; J\* [DV] independently reproduced. **Not covered by d: the k = 3 piece of the corrections.** Costed the same way, it moves J\* 30 → 33, **50 → 64**, 64 → 80 (θ = 1) |
| **5** | ⟨TTO_{Δ,J}⟩ conserved structures in d = 3: generic even J 2⁺ + 1⁻, odd J 1⁻; integer J = 0…9: 1+1, 0, 1+1, 0, 2+1, 0+1, 2+1, 0+1, 2+1, 0+1 (P2_SPINNING §1) | DKKPS `TTTT.tex` table: O₀ 1⁺+1⁻, O₂ 1⁺+1⁻, O_{2n≥4} 2⁺+1⁻, O_{2n+1≥5} 1⁻, plus "a unique parity-odd tensor structure … for odd ℓ ≥ 5". Exact match. An independent helicity count (h = ±2, \|h₁−h₂\| ≤ ℓ, parity, Bose (−1)^ℓ) reproduces all totals for ℓ = 0…11. It needs one convention datum for which parity the (2, −2) orbit carries | **yes** | [Q] verified + [D] (counting) |

**Net.**
* **No load-bearing claim was refuted.**
  * The de-shift identity is reproduced and strengthened to an identity in symbolic b.
  * The monopole constant D is reproduced to 8 digits, and its t-independence now has a two-line proof.
  * Lemma LR∞ is reproduced, with a true remainder six orders of magnitude below the stated r.
  * The order-of-limits argument is sound and its quotes are accurate.
  * Theorem W's algebra is correct, and its J\* entries are reproduced with my own costs.
* **Discrepancies (none overturns a stated result; details in §6):**
  1. LARGESPIN_L's formula for R uses ∫\|ψ\|, but the r(β) table was computed with a cruder constant. The table is conservative by
     ×36 at β = 1.
  2. Theorem W's cost d omits the k = 3 (Ĉ₃) part of the corrections (LARGESPIN_W caveat 1 only). Included with the same structures,
     it moves the headline J\* = 50 to 64.
  3. Route W's cached `4a`/`4b` costs are stored with 1/m², not the stated 1/m⁴. The J\* tables are nevertheless consistent with 1/m⁴.
  4. Several small citation labels in REGGE_TREE are wrong.
  5. "Gershgorin lowers tolerances by ≤ 1.5" needs an unstated assumption on f_OO.

---

## 1. Target 1: the EE kernel's de-shift identity and the monopole

### 1.1 Blocks (`blocks2.py`)
* **Operator.** The d = 3 Dolan–Osborn ½-Casimir with (a, b) = (0, b) in radial variables, multiplied by 2(1−r²)Den. It is the operator
  recorded in INDEP_CHECK §1.1 and checked there against Dolan–Osborn to 10⁻³¹.
* **Expansion.** G = r^Δ Σ_N r^N Σ_m B[N, m] P_{l+m}(η), with B[0, 0] = 1, so Ĝ → P_l at r → 0. This matches the (4r)^{Dp} normalisation of
  PW_IMAGES: the r⁰ term of (2.1) is the Legendre recursion.
* **Recursion.** Generic in the field: it runs with sympy expressions (l and b symbolic, exact cancellation) or with mpmath numbers.
* **Coefficients.** 1+i/2ν = 2l/(2l+1) and 1−i/2ν = 2(l+1)/(2l+1), since ν → ν ± i means l → l ∓ 1.
* **Basis.** Both sides of (2.1) are expanded in the common basis r^N P_{l+m}: a Dp = 0 block at l ∓ 1 contributes to P_{l+(m∓1)}.
* **Parity rule.** m ≡ N (mod 2). So the Dp = 0 denominator zero at (N, m) = (3, 0), which vanishes for every l, is never reached. Its
  odd orders vanish identically at b = 0.

### 1.2 The identity (`t1_magic.py`, `t1d_pointwise.py` → `cache/t1_magic.txt`, `cache/t1d_pointwise.txt`)

| test | result |
|---|---|
| (a) exact, **l and b symbolic**, N ≤ 4 / 5 / 6 | 0 / 0 / 0 nonzero residual coefficients (of 20 / 27 / 35) |
| (a2) exact, l symbolic, b = −4, …, +4, N ≤ 8 | 0 nonzero residuals, of 30 (b = 0) and 54 (each b ≠ 0) coefficients |
| (b) 50 digits, N ≤ 40, ν ∈ {0.7, 5, 2.3+0.2i, 2.078−0.279i}; b = −4…4, 0.37, 7, 10 | max \|residual\|/max\|LHS\| at the same order ≤ **3.7·10⁻⁴⁸** (typically 10⁻⁴⁹) |
| (c) controls at b = 0, 2, −3 | swapped coefficients: 0.42–0.64; Dp = 1 on the right: 0.34–0.69; (Dp_L, Dp_R) = (2, 1): 0.52–0.80. The blocks themselves depend on b (\|B^{b=2} − B^{b=0}\| ~ 2·10³) |
| (d) pointwise, mpmath `legenp` of complex degree, 4 points (r, β), b = −4…4, two ν | residual at N = 30 / 45 / 60 falls geometrically to the 30-digit floor (e.g. b = 0, (0.2, 0.6): 7·10⁻¹⁷ / 4·10⁻²³ / 3·10⁻³⁰). Negative b converges more slowly (b = −4, (0.3, 0.4): 0.04 / 6·10⁻⁶ / 3·10⁻¹⁰), as INDEP_CHECK §1 found for the direct negative-b series |

**Reading.**
* The identity holds as a polynomial identity in b at each order checked. That is stronger than the original (b = 0…4, 7, 10 sampled).
* It also covers both sign conventions for the DO b (INDEP_CHECK: b_DO = −b), and the full EE range b = −4…4.

### 1.3 PW2 side claims (`t1c_even_analytic.py` → `cache/t1c_even_analytic.txt`)
* **Frame.** With x₁ = −rŷ, x₃ = rŷ, x₂ = −x̂, x₄ = x̂, ŷ² = x̂² = −1 and −ŷ·x̂ = η, and the standard cross-ratios u = x₁₂²x₃₄²/(x₁₃²x₂₄²),
  v = x₁₄²x₂₃²/(x₁₃²x₂₄²): **v − u = η(1+r²)/(2r) exactly** (sympy). This assumes the project's u′, v′ are these cross-ratios.
* **Evenness.** B(l)[N, m] = B(−1−l)[N, −m] to 8.9·10⁻⁶⁰. Since P_{−1−l+m} = P_{l−m}, this is Ĝ⁰_b(μ) = Ĝ⁰_b(−μ).
* **Removability** at μ = ±i/2 (l = 0, −1) and ±3i/2 (l = 1, −2), b = −4…4, N ≤ 14:
  * max\|B\| = 3.62762·10⁹ (respectively 5.67948·10⁶) at ε = 10⁻¹², 10⁻²⁴ and 10⁻³⁶;
  * \|B(10⁻³⁶) − B(10⁻²⁴)\|/max\|B\| ≤ 7·10⁻²²;
  * so the would-be poles are removable.
* **(2.1) ⇒ (2.2).** Multiply by N r² C_b v^{b/2}/(4r) and sum over b. The C_b, the v^{b/2} and f are ν-independent, so this is linear
  algebra. It is correct.

### 1.4 The monopole (`t1b_monopole.py` → `cache/t1b_monopole.txt`)
* **Definitions** (from the README of `prl/data` and the Abel equation):
  * 𝒜[f] = (R₋f(ν+i) + R₊f(ν−i))/(R₋+R₊), with the EE R∓ of `kernels.json`;
  * cost = Re Σ_k dν X_k 𝒜[f](ν_k)/Λ;
  * target Σψ₂ f dν/Λ.
* **Convention check.** The orientation (R₋ with f(ν+i)) was confirmed by the de-shift identity itself: defect[P_t] ≤ 3·10⁻¹⁵ at
  t ≤ 6.
* **Legendre functions.** P_{−½+iμ}(cosh t) for complex μ comes from a vectorised Mehler–Dirichlet quadrature. It agrees with mpmath to
  ≤ 5·10⁻¹³, also at μ = 1000 + i.
* **Exact mechanism (new).** Partial fractions give `1/c(μ) = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9)`.
  * The 1/32 and 1/(4μ²+9) pieces are analytic in |Im μ| ≤ 1, so their defects vanish: measured −2·10⁻¹⁵ and −8·10⁻¹⁹.
  * For even F in the Mehler class, (F − F(i/2))/(4μ²+1) is strip-analytic. Its measured defect is ≤ 1·10⁻¹⁸ at every t.
  * Hence

        defect[F/c] = (15/2)·F(i/2)·δ₁,   δ₁ := defect[1/(4ν²+1)] = −1.331109479·10⁻⁴ (1000), −1.463392666·10⁻³ (300).

    This is exactly (4.1), with **D = (15/2)δ₁**. The t-independence is a consequence, not a numerical coincidence: [D].

| Δ_g | D from (15/2)δ₁ | D direct (defect of 1/c) | D·Λ² | PW_IMAGES |
|---|---|---|---|---|
| 1000 | −9.9833211·10⁻⁴ | −9.9833211·10⁻⁴ | −995.3 | −9.983321·10⁻⁴ |
| 300 | −1.0975445·10⁻² | −1.0975445·10⁻² | −977.9 | −1.097545·10⁻² |

The defect of P_t/c at t = 0, 0.01, 0.03, 0.1, 0.3, 1, 2, 3, 6 equals D in every printed digit (8 significant). The table of Ω₂, Ω_c and
Ω̃ = Ω_c + D is reproduced entry by entry, for example:

| t | Ω₂ | Ω_c | Ω̃ |
|---|---|---|---|
| 0 | 10.386 | 0.32562 | 0.32463 |
| 0.03 | 1.947·10⁻² | 1.4422·10⁻³ | 4.439·10⁻⁴ |
| 0.1 | 5.214·10⁻⁴ | 6.852·10⁻⁴ | −3.131·10⁻⁴ |
| 1 | 7.582·10⁻⁶ | 2.678·10⁻⁴ | −7.306·10⁻⁴ |
| 6 | 7.034·10⁻⁷ | 6.314·10⁻⁶ | −9.920·10⁻⁴ |

* **The claim \|Ω̃\| ≤ 1.0·10⁻³ for t ≥ 0.03** holds at the sampled t. Since Ω_c ≥ 0 decreases from 1.44·10⁻³, Ω̃ ∈ (D, 4.4·10⁻⁴] on
  the samples.
* **Mixtures** F = Σ w_k(1 + a_kν²)P_{t_k}, with F(i/2) = Σ w_k(1 − a_k/4): (4.1) holds to 1.6·10⁻¹⁵ … 5.4·10⁻¹⁰ (1000) and
  7·10⁻¹⁴ … 7·10⁻¹³ (300).

## 2. Target 2: Lemma LR∞

### 2.1 Re-derivation of the large-β form
* **Step 1.** Mehler–Dirichlet gives Ω_k(β) = (√2/π)∫₀^β Φ(Λt)(cosh β − cosh t)^{−½}dt, with Φ(x) = ∫₀¹ψ cos(px)dp.
* **Step 2.** For a polynomial ψ with ψ(0) = 0, repeated integration by parts gives
  Φ(x) = Σ_{j odd}(−1)^{(j+1)/2}ψ^{(j)}(0)/x^{j+1} + Re Σ_j(−1)^jψ^{(j)}(1)e^{ix}/(ix)^{j+1}.
  * The leading term is −a₁/x².
  * The p = 1 terms start at j = 5, because (1−p)⁵ divides ψ.
* **Step 3, the −a₁/x² term.** Its Hadamard finite part against K is (π/√2)·M₁. Equivalently ∫₀^∞ν P_{−½+iν}dν (f.p.) = M₁: write
  ν = ν tanh πν + 2ν/(e^{2πν}+1) and use the Mehler–Fock orthogonality of the first piece for β > 0.
* **Step 4, the −a₂/(2b̂³) term.** It is the flat-space (Bessel) moment, ∫₀^∞p²J₀(b̂p)dp = −1/b̂³ (regulated). In the Mehler form it is
  the K₀t²/(4D) term.
* **Numerical confirmation.**
  * The identity (√2/π)[K₀/β − ∫₀^β(K−K₀)t⁻²dt] = M₁(β) holds to 15 digits at β = 0.1, 1, 4, 10.
  * The (K − K₀)/t² integrand had to be evaluated in cancellation-free form.

### 2.2 Numbers (`t2_lrinf.py` → `cache/t2_lrinf.txt`, `cache/t2_lrinf_extra.txt`)
**Method.**
* mpmath at 40 digits.
* Gauss–Legendre panels of phase ≤ π/2 while Λt ≤ 3000; beyond, only the smooth part of Φ (the oscillatory remainder is < 10⁻²⁵ there).
* t = β − u² near the endpoint.
* Exact NF02 polynomials.
* Control: a direct p-quadrature with `legenp` agrees to 12 digits at β = 0.3 and 1.

**Constants.**
* a₁ = 188.8961 / 50.64297 and a₂ = −1023.5457 / 7.3·10⁻⁵ (k = 2 / 3).
* A₁ = 327.1492, −A₂ = 2047.091, ∫ψ₂ = 10.385947.

| β | Ω₂ | E₂/two-term | Ω₃ | \|E₃\|/(a₁⁽³⁾M₁/Λ²) | r_actual | stated r | F_LR,min = 2Ω₂ − \|Ω₃\| | stated bound | ratio |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 7.570137·10⁻⁶ | 3.1·10⁻⁷ | 1.907797·10⁻⁶ | 1.2·10⁻⁶ | 5.7·10⁻⁷ | 0.577 | 1.32325·10⁻⁵ | 6.12142·10⁻⁶ | 2.16 |
| 1.5 | 6.404700·10⁻⁶ | 2.8·10⁻⁷ | 1.686119·10⁻⁶ | 2.9·10⁻⁶ | 7.8·10⁻⁷ | 0.070 | 1.11233·10⁻⁵ | 1.03608·10⁻⁵ | 1.074 |
| 2.0 | 5.382527·10⁻⁶ | 2.5·10⁻⁷ | 1.432442·10⁻⁶ | 2.7·10⁻⁶ | 7.0·10⁻⁷ | 0.015 | 9.33261·10⁻⁶ | 9.19381·10⁻⁶ | 1.015 |
| 3.0 | 3.527467·10⁻⁶ | 1.6·10⁻⁷ | 9.43929·10⁻⁷ | 1.8·10⁻⁶ | 4.7·10⁻⁷ | 1.6·10⁻³ | 6.11100·10⁻⁶ | 6.10125·10⁻⁶ | 1.0016 |
| 5.0 | 1.238826·10⁻⁶ | 4·10⁻⁸ | 3.32051·10⁻⁷ | 4.4·10⁻⁷ | 1.1·10⁻⁷ | 1.0·10⁻⁴ | 2.14560·10⁻⁶ | 2.14539·10⁻⁶ | 1.0001 |
| 10 | 6.838588·10⁻⁸ | −2·10⁻⁸ | 1.833418·10⁻⁸ | 2.2·10⁻⁷ | 5.7·10⁻⁸ | 1.0·10⁻⁴ | 1.18438·10⁻⁷ | 1.18426·10⁻⁷ | 1.0001 |
| 32 | 4.051824·10⁻¹³ | −4·10⁻⁹ | 1.086292·10⁻¹³ | 3.9·10⁻⁸ | 1.0·10⁻⁸ | 1.0·10⁻⁴ | 7.01736·10⁻¹³ | 7.01665·10⁻¹³ | 1.0001 |

* All 19 grid points β ∈ {1, 1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32} satisfy F_LR,min > 0,
  F_LR,min ≥ the stated bound and r_actual ≤ the stated r. F_LR,min = 2Ω₂ − \|Ω₃\| is the minimum over u ∈ [0, 1] and σ = ±1, i.e. both
  helicities and all m ≥ Λ.
* **Below β = 1** (outside the lemma): r_actual = 1.35·10⁻⁴ (β = 0.3), 1.1·10⁻⁵ (0.5), 2·10⁻⁶ (0.7), 6·10⁻⁷ (0.85). This is consistent
  with LARGESPIN_L's remark that at β = 0.3 the k = 3 two-term form misses ~10⁻³.
* **J\*(m) table** (`t2b_jstar_m.py`), with m² = (Δ−J−2)(Δ+J−1) (CHMRS d = 3, the same as Route W's cache):
  * J\* = 857, 987, 1275 at m = 1, 1.1, 1.3 m_gap;
  * J\* = 780, 1041, 2602 at m = 1.5, 2, 5 Λ;
  * the β ≥ 1.5 and β ≥ 2 rows are also reproduced exactly;
  * the rounded constant 1.0421 in "⌈(1.0421m−1)/2⌉" gives 2601 instead of 2602 at 5Λ; the table used 2 sinh ½ = 1.0421906.

### 2.3 The r(β) table (`t2c_rbound.py` → `cache/t2c_rbound.txt`)
I implemented the stated bound literally: T_a, T_b (closed form), T_c (quadrature), with R = min(R_∞, a₁/x² + ∫\|ψ\|). The grid was the
stated one (0.01 / 0.05 / 0.5).

| β ≥ | 1.0 | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 64 |
|---|---|---|---|---|---|---|---|
| sup r, formula as written (this check) | **0.0162** | 0.0078 | 0.0032 | 1.06·10⁻³ | 2.4·10⁻⁴ | 4.7·10⁻⁵ | 1.4·10⁻⁷ |
| stated | 0.577 | 0.225 | 0.070 | 0.015 | 1.6·10⁻³ | 1.0·10⁻⁴ | 1.6·10⁻⁷ |

* **Diagnosis** (I read `lsl_lrinf.py` lines 38–107 after obtaining the above). The code uses Φ_max = Σ\|c_n\|/(n+1), commented
  "≥ ∫\|ψ\| (crude, rigorous)", in place of ∫\|ψ\| = 10.39. With polynomial coefficients up to ~10¹⁰ this inflates the small-x branch
  of R, and hence T_c₁.
* Both choices are rigorous upper bounds. The stated table is therefore valid but conservative, and LARGESPIN_L §3 should say
  "Φ_max = Σ\|c_n\|/(n+1)".
* **Side result** (not claimed by the project): with the tighter ∫\|ψ\|, r_bound = 0.28, 0.13, 0.069, 0.040, 0.025, 0.020 at
  β = 0.5, 0.6, 0.7, 0.8, 0.9, 0.95. Lemma LR∞'s threshold could therefore be lowered to β ≥ 0.5, i.e. **J\* = 324 at the gap instead of
  857**. This would need:
  * a certified ∫\|ψ₃\| (here computed numerically);
  * a denser grid on [0.5, 1].

## 3. Target 3: the order of limits (REGGE_TREE route N), read against the sources

### 3.1 Is the logic sound?
**Yes.** Route N is:
1. the exact sum rule at finite c_T, where J₀ ≤ 1 (KKSDZ (4.55), (4.70), (4.73));
2. (Pos) for every state with Δ − J ≥ Δ_g, a c_T-independent statement about (Δ, J, parity);
3. hence Ω[light](c_T) ≤ 0 for every c_T;
4. Ω[light] = L₁/c_T + o(1/c_T) by A2;
5. so L₁ ≤ 0.

This is word for word the CHMRS procedure: §1.1.3, eqs. (1.12) ω|light, ω|heavy; (1.13) "−ω|light = ω|heavy"; (1.14) "−ω|light ≥ 0",
with "It is easy to see that ω|light = O(ε)". The heavy side enters only through its sign. No heavy-side 1/c_T expansion, convergence of
an O(1/c_T) heavy sum, or planar intercept is used.

**Two points a referee may still press, both acknowledged in REGGE_TREE §4 and §6.**
* (i) The interchange of the 1/c_T expansion with the infinite light double-twist sums is now an assumption about the **light** side
  (A2, uniformity in spin [A]). That is where a KKSDZ §4.5-type non-commutation could re-enter. It concerns the lightcone (fixed-twist,
  large-J) regime, not the Regge limit, and CHMRS make the same assumption implicitly ("It is easy to see …").
* (ii) Route N needs (Pos) for all states above the gap, multi-traces included. That is a (Δ, J)-generic statement, so nothing new is
  needed, but the Letter's wording "heavy single traces" must change, as REGGE_TREE §6 says.

### 3.2 Are the quotes accurate? (sources read: `sec_introduction.tex`, `sec_light.tex`, `draft.tex`, `superconvergence.tex`, PDFs)

| REGGE_TREE quote | source | verdict |
|---|---|---|
| CHMRS "In weakly coupled AdS EFTs, where tree-level interactions are suppressed by ε (which we assume to be much smaller than any power of 1/Δ_gap)" | §1 (before §1.1), p. 3 | verbatim. **Label:** given as "§1.1.1"; the parenthetical "(§1, p. 3)" is right |
| "we make no assumptions about the primaries with τ > Δgap … can be strongly coupled" | §1.1.1 | verbatim |
| (1.12)–(1.14), "It is easy to see … −ω\|light ≥ 0" | §1.1.3 | verbatim; equation numbers correct |
| "fully rigorous … (with Regge intercept ≤ 1) … non-perturbative level" | §1.1.3 | verbatim |
| §4.1 light data = tree-level EFT | §4.1 | faithful |
| §4.4 (4.45) "a universal prediction for the k = 2 moment … in any holographic CFT" | §4.4, eq. `universalreggegravity` = (4.45) | verbatim, number correct |
| §4.4 eikonal: "Naively expanding … divergent integral ∼∫dΔ Δ … nice consistency check … the eikonal model satisfies the sum rules" | §4.4 | verbatim. CHMRS add at once "this success cannot by any means be used to justify the eikonal approximation", which does not affect the use made of it |
| (2.8) lim M/s² = 0, flat-space | §2.1, (2.8) | correct |
| §5.7 "slightly stronger, nonperturbative result … twice-subtracted" | §5.7 | verbatim |
| footnote 9 "(G − 1) vanishes … second sheet … never been rigorously proven" | footnote 9 (§2.4, printed at the top of §2.5) | verbatim |
| "§5.3 (conclusion of the eq. (5.3) discussion) is the technical assumption" | the technical assumption is **eq. (5.3) in the §5 preamble**, not §5.3 | label wrong (immaterial) |
| CLS-D §3.1 "J₁+J₃−1 > J₀" | §3.1 | verbatim |
| CLS-D §3.4 "we only know that J₀ ≤ 1"; "becomes J₀ < k, which is indeed true for any k ≥ 2" | §3.4 | verbatim |
| CLS-D §3.5.2 "k > J₀ for even j′ and k − 1 > J₀ for odd j′ … k ≥ 2 … k ≥ 3" | §3.5.2 | verbatim. For J_e = 2: f⁺_{2;2} = v′−u′ and **f⁻_{3;2} = −1** (REGGE_TREE writes "f⁻ = 1"; sign immaterial) |
| CLS-D §2.2.2 "We leave an exploration of the Regge regime … to future work" | §2.2.2 | verbatim |
| CLS-D §1 footnote 1 "One must additionally check other regimes" | §1, footnote | verbatim |
| "This is exactly the logic of … CLS-D (§2.1, §3.4)" | CLS-D restates the CHMRS light/heavy logic in **§2.2.2** (ω\|light = heavy sum, positivity ⇒ −ω\|light ≥ 0); §2.1 is the light-transform review | **label wrong.** Also, CLS-D never discusses the c_T expansion or the order of limits. "CLS-D restate CHMRS's logic and use J₀ ≤ 1" is accurate; "the logic of CLS-D" overstates |
| KKSDZ (4.55), (4.70), (4.73); §4.2.6 "we don't have a general proof" | §4.2–4.4 | verbatim, numbers correct |
| KKSDZ §4.4 "3 > J₀ = J₀^dDisc" | eq. (4.79), §4.4 | correct |
| KKSDZ §4.5 (4.80) "J₀^planar ≤ 2 … saturated in holographic theories … tree-level graviton exchange" | §4.5, (4.80) | verbatim |
| KKSDZ (4.81)–(4.83): "does not commute with taking the Regge limit … term-wise integration is reliable for the orders at which it converges" | §4.5 | verbatim |
| KKSDZ §5.3 planar sum rules term-wise (single traces) | §5.3 "Structure of the general sum rule" | faithful |

**Conclusion.** The literature attribution in REGGE_TREE is accurate in substance. The corrections are three section labels and one
overstatement (CLS-D). The verdict "no paper states a tree-level heavy-sum Regge assumption for k = 2" agrees with what I read.

## 4. Target 4: Theorem W (`t4_theoremW.py` → `cache/t4_theoremW.txt`)

### 4.1 The implication
* **2×2.** Diagonal ∈ [f_h − d, f_h + d], \|offdiagonal\| ≤ d. Then (f₁ − d)(f₂ − d) ≥ d² ⇔ d ≤ H := f₁f₂/(f₁+f₂), and this is sufficient
  for M ⪰ 0.
  * 2·10⁵ random Hermitian matrices satisfying the premises: 0 indefinite.
  * The worst case at d = 1.001H is indefinite in 1000/1000 cases, so the criterion is sharp.
  * LARGESPIN_W writes "iff". That is correct for the worst case, and the step needs only "if".
* **Hermiticity.** Bounding \|M₁₂\| by \|κ\|·\|cost\| needs the off-diagonal image coefficients to be Hermitian (κ^{h′h} = κ̄^{hh′}).
  Then M₁₂ = κ₁₂·Re(∫X𝒜[S]) = κ₁₂·cost. This is implicit in P2 and should be stated.
* **Conclusion bound.** Q(λ) ≥ (1−θ)min f\|λ\|²/2 holds. Numerically, λ_min ≥ 2 × that bound over 2·10⁵ random (f, θ).
* **Gershgorin (3×3).** d ≤ min f/3 suffices: 0/10⁵ violations. It is exact for equal floors.
  * "Lowers tolerances by ≤ 1.5" is H/(min f/3) = 1.5 only when f₊₊ = f₊₋ and f_OO ≥ min(f₊₊, f₊₋).
  * Examples: floors (1, 10, 10) give 2.7 and (1, 3, 1) gives 2.25 (exact worst-case 3×3 tolerances 0.83 and 0.43).
  * On Route W's sampled even-J floors f₊₊/f₊₋ ∈ [0.965, 1.000], giving 1.50–1.53. f_OO is not tabulated.

### 4.2 Own costs
**Definitions.**
* Lemma-I (de-shift, R1) cost of structure S: 2Ω₂[S]/m², with Ω₂[S] = ∫₀¹ψ₂(p)S(Λp; β)dp, exact NF02 ψ₂.
* S ∈ {ηP, ην²P, P′(η), (η²−1)ν²P′(η), (η²−1)P′(η), Re P_{ν∓i}, Re P_{ν∓2i}}, with P′ from (η²−1)P′_μ = μ(ηP_μ − P_{μ−1}).
* The O(1/m⁴) shapes {η²ν⁴P, η sinh β ν⁴∂_βP} are divided by m⁴.
* Stripped (R2) cost: 2 Re Σ X_k S(ν_k)dν/Λ/m², with the published X.

**Comparison with Route W's cache** (47 states):
* R1: ≤ 1.1·10⁻³ (A), 1.3·10⁻⁴ (B), 4·10⁻⁵ (C1, C3), 2.3·10⁻² (C2), 1.1·10⁻² (sh±1), 5·10⁻⁴ (sh±2).
* R2: ≤ 10⁻⁶ except C2 (2·10⁻⁴) and 4a/4b (9·10⁻⁴).
* **Units.** The cache stores `4a`/`4b` with 1/m², i.e. ×m² relative to P4's stated 1/m⁴ units (the ratio to my costs is m² to 10⁻³).
  The J\* tables are reproduced only with 1/m⁴, so the J\* code must divide. This is a documentation inconsistency in
  `costs_1000.{json,txt}`, not an error in J\*.

### 4.3 J\* recomputed (floors from Route W's cache; every cost my own)
**Convention that reproduces the tables.** Ring classes in R1 form, index-shifted harmonics in **stripped** form (LARGESPIN_W §3.4:
"tolerated at K_sh ≥ 52 (stripped)"), k = 2 only, d = max over helicity entries.

| (Λ_amp, K_sh, K₄), θ | claimed | own costs |
|---|---|---|
| (1, 10, 1), 1 | 30 | **30** (fails 26) |
| (1, 10, 1), ½ | 38 | **38** |
| **(3, 10, 10), 1** | **50** | **50** (fails 26–46) |
| (3, 10, 10), ½ | 100 | **100** |
| (10, 30, 10), 1 | 64 | **64** |
| (10, 30, 10), ½ | none (fails 300, 344) | **none (fails 300, 344, and < 100)** |
| (3, 10, 10), R2 everywhere, 1 | none (fails 300, 344) | **none (fails 300, 344, and < 50)** |
| P2 §7.5 model (Λ_A, Λ_B, Λ_C2, Λ_C3) = (3.17, 1, 1, 2.03), K₄ = 4, K_sh = 10, θ = 1 / ½ | 38 / 50 | **38 / 50** (Λ_C1 = 1 or 3.17 alike) |

**Note.** If the index-shifted terms are also given their R1 cost, (10, 30, 10) at θ = ½ gives J\* = 100 with no failure at 300/344. The
failures there are entirely the stripped index-shift cost.

### 4.4 Qualifications
* **The k = 3 piece.** d contains only the k = 2 action of the corrections (LARGESPIN_W caveat 1), while the floors contain the k = 3
  term uσΩ₃.
  * Costing the same structures with the same coefficients in the full hybrid functional, (2Ω₂[S] ± uΩ₃[S])/m², worse helicity, changes
    the costs by up to 95 % (A–C3) and ×8 (C2).
  * **J\* becomes 33 (1, 10, 1), 64 (3, 10, 10) and 80 (10, 30, 10) at θ = 1.** At (3, 10, 10) the states J = 50 and 51 near the gap
    fail at d/H = 1.24 / 1.23 (k = 2 only: 0.85).
  * This is optimistic for k = 3: PW_IMAGES caveat 5 says Ĉ₃ corrections are counted raw with \|ψ₃\|.
  * **Theorem W's premises should list "the k = 3 part of the O(1/m²) corrections is negligible" explicitly**, or include it. The
    θ = ½ values are unchanged.
* **Gershgorin factor** "≤ 1.5": see §4.1.

## 5. Target 5: ⟨TTO_{Δ,J}⟩ counts (`t5_counts.py` → `cache/t5_counts.txt`)
* **DKKPS 1708.05718, `TTTT.tex` §2 table** (Δ strictly above unitarity): O₀ 1⁺+1⁻; O₂ 1⁺+1⁻; T 2⁺+1⁻; O_{2n≥2} 2⁺+1⁻; O_{2n+1≥5} 1⁻.
  The appendix adds "There exists a unique parity-odd tensor structure for ⟨TTO_ℓ⟩ for odd ℓ ≥ 5" and two parity-even structures for
  even ℓ ≥ 4. Odd ℓ = 1, 3 are absent, i.e. zero.
* **Comparison.** This is exactly P2_SPINNING's quotation, its generic counts (even 2⁺+1⁻, odd 1⁻) and its integer table J = 0…9.
* **Independent count.** Conserved spin-2 operators carry helicities ±2, and a generic spin-ℓ operator admits (h₁, h₂) iff
  \|h₁ − h₂\| ≤ ℓ (Kravchuk–Simmons-Duffin counting). Parity pairs (h₁, h₂) ↔ (−h₁, −h₂), and Bose symmetry acts with (−1)^ℓ.
  * The (2, 2) orbit gives 1⁺+1⁻ for even ℓ and nothing for odd ℓ.
  * The (2, −2) orbit exists for ℓ ≥ 4 and gives one structure of parity (−1)^ℓ. Its absolute parity is fixed by one datum, e.g. ℓ = 4.
  * Totals: ℓ = 0…11 → 1+1, 0, 1+1, 0, 2+1, 0+1, 2+1, 0+1, 2+1, 0+1, 2+1, 0+1. All agree.
* **Remarks.**
  * The single odd-J structure lies in the (2, −2) = "+−" helicity class and is parity-odd. That agrees with P2's correction of
    LARGESPIN_W §1.
  * Not re-done: the explicit Todorov-operator solve over ℚ(Δ, J). The counts are verified; the explicit bases are not.

## 6. Discrepancies (all non-fatal)

| # | where | issue | effect |
|---|---|---|---|
| D1 | LARGESPIN_L §3 statement vs r(β) table | The text defines R with ∫\|ψ\|; the table was computed with Φ_max = Σ\|c_n\|/(n+1). Evaluated as written, r(1) = 0.016, not 0.577 | Table valid but conservative (×36). The text should name the constant used. Opportunity: β-threshold 0.5, J\* ≈ 324 (uncertified) |
| D2 | LARGESPIN_W §4, PW_IMAGES §5 | d omits the k = 3 part of the corrections (caveat 1 only) | With the same structures costed in the full hybrid functional, J\* = 50 → 64 (θ = 1, (3, 10, 10)); 30 → 33; 64 → 80. Should be a named premise |
| D3 | `lsw/cache/costs_1000.{json,txt}` | `4a`/`4b` stored with 1/m², not the stated 1/m⁴; the `costs_1000.txt` ratios (e.g. "4a: 1.1e+04") are in those units | Documentation only. J\* is consistent with 1/m⁴ |
| D4 | LARGESPIN_W §4.0 Gershgorin remark | "lowers tolerances by ≤ 1.5" needs f₊₊ ≈ f₊₋ and f_OO ≥ min f | Factor up to 3 otherwise; f_OO not tabulated |
| D5 | REGGE_TREE §0, §2 | CHMRS ε-quote is in §1, not §1.1.1; the technical assumption is eq. (5.3) in the §5 preamble, not §5.3; CLS-D's restatement of the logic is §2.2.2, not §2.1, and CLS-D does not discuss the order of limits; f⁻_{3;2} = −1 | Citation hygiene only |
| D6 | PW_IMAGES §4.1 | "D is [DV]": the t-independence and the form of (4.1) follow exactly from 1/c's partial fractions; only the value of δ₁ is numerical | Upgrade, not a defect |
| D7 | PW_IMAGES §2.1 | "[P] to O(r⁸) for b = 0…4" undersells it | The identity holds with b symbolic (exact to O(r⁶)) and for b = −4…4 to O(r⁸) |

## 7. Reproduce

    cd ads4/code/indep2 && bash run_all.sh          # ≈ 1.5 CPU-h, at most one process at a time

| file | content |
|---|---|
| `blocks2.py` | (0, b) block recursion over any field (sympy / mpmath), magic-identity residual in the r^N P_{l+m} basis |
| `t1_magic.py` | (a) exact, l and b symbolic; (a2) b = −4…4 exact to r⁸; (b) 50 digits to r⁴⁰; (c) controls; (d) pointwise |
| `t1d_pointwise.py` | pointwise identity vs truncation order (N = 30/45/60), per b |
| `t1c_even_analytic.py` | evenness, removability at ±i/2, ±3i/2, frame identity f = η(1+r²)/(2r) |
| `t1b_monopole.py` | D by partial fractions and directly, t-independence, controls, Ω̃ table, mixtures (Δ_g = 1000, 300) |
| `t2_lrinf.py` | Ω₂, Ω₃ at 40 digits, two-term form, F_LR vs the stated bound, Hadamard identity (β = 1…32; argv for extra β) |
| `t2b_jstar_m.py` | J\*(m) table |
| `t2c_rbound.py` | the stated remainder bound evaluated as written, sup r on the stated grid (argv for extra β) |
| `t4_theoremW.py` | 2×2 / Gershgorin algebra; own R1/R2 costs; J\* tables (floors read from `code/lsw/cache/costs_1000.json`) |
| `t5_counts.py` | helicity count of conserved ⟨TTO_ℓ⟩ structures |
