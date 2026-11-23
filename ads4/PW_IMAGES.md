# PW_IMAGES: the EE functional's plane-wave images beyond leading Regge. An exact scalar de-shift identity, the Regge-spin-0 monopole, and what they do to Theorems W and L

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. `P2_SPINNING.md` and `code/p2spin/` were not touched.*
* *New code:* `ads4/code/pwimg/` (`pw_common.py`, `pw1_magic.py`, `pw1b_e2e.py`, `pw2_radial.py`, `pw3_cost.py`, `pw4_theoremW.py`,
  `pw5_theoremL.py`, `run_all.sh`). Caches in `code/pwimg/cache/`; its `.gitignore` ignores everything except itself.
* *Read-only imports:* `eestep/ee_blocks.py`, `eestep/ee_quad.py` (the EE kernel, its 9-block reduction and the (0, b) block
  recursion), `h7i14/v3_regge_kernel.py`, `h7i14/v4_hclosed.py`, `spinfunc/u2_hclosed_struct.py`, `spinfunc/s1_offdiag_kernel.py`
  (exact leading-Regge heavy images), `n2sub/n2_common.py` (published X, ψ₂, R∓ of EE), `lsw/lsw_common.py` + `lsw/cache/costs_1000.json`
  (Route W's states, floors and structures), `lsl/lsl_core.py` (F_LR, Ω_k). No bytecode was written.
* *Compute:* ≈ 1.5 CPU-h, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as UNIFIED_THEOREM §0. Notation as LARGESPIN_L/W: `l = −½ + iν` (Regge spin), `P_μ := P_{−½+iμ}(cosh β)`, H7 frame
`x₁ = −rŷ, x₃ = rŷ, x₂ = −x̂, x₄ = x̂` (|x| = 1, η = cosh B), `f := v′ − u′`. `Ĝ^{Dp}_b(r, η; μ) = Σ_n r^n Σ_m B^{Dp}_b[n, m](μ) P_{l+m}(η)`
is the (0, b) continuous-spin block of EESTEP §1.3 divided by (4r)^{Dp}; `E^{Dp}_μ := Σ_{b=0}^{4} C_b(x; pols) v^{b/2} Ĝ^{Dp}_b` with the
ν-independent EE dressing coefficients `C_b` of EESTEP §1.2. The raw EE per-ν kernel is `K^EE_ν = N(ν) r² E¹_ν` (EESTEP §1.4).
`c(μ) = 32(4μ²+1)(4μ²+9)/((4μ²+25)(4μ²+81))` (EE_FUNCTIONAL §2). "n2 units": `Σ ψ₂ f dν/Λ = ∫ψ₂ f dp` on the published grid; Theorem L's
`F_LR` and Route W's floors are in these units.

---

> **Update (INDEP_CHECK2.md, independent reproduction with new code): all five targets reproduced, nothing load-bearing refuted. Two substantive qualifications: Theorem W's J\* = 50 → 64 once the k = 3 piece is included, and Lemma LR∞'s r(β) table uses a cruder constant than the stated formula.**
> * **Magic-block identity:** exact with l and b symbolic to O(r⁶), so it is an identity in b. It is exact to O(r⁸) for every b = −4…4, and holds to 4·10⁻⁴⁸ at 50 digits to r⁴⁰ (including b = 0.37, 7, 10 and complex ν); controls fail at O(1). **The monopole constant D is reproduced to 8 digits, and its t-independence is now derived:** 1/c = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9), so D = (15/2)·defect[1/(4ν²+1)].
> * **Lemma LR∞:** reproduced; the true r(β) is ≤ 8·10⁻⁷ on β ∈ [1, 32]. **Discrepancy:** the stated formula (with ∫|ψ|) gives r(1) = 0.016, while the table (0.577) used Σ|c_n|/(n+1). Both are rigorous; the text must name the constant used. With the tighter constant the bound holds down to β = 0.5 (r = 0.28), which would give J\* = 324 at the gap instead of 857 (needs a certified ∫|ψ₃| and a denser grid).
> * **Order of limits (route N):** sound, and the quotes are verbatim or faithful. Minor fixes:
>   * the CHMRS ε quote is in §1, not §1.1.1;
>   * the "technical assumption" is eq. (5.3), not §5.3;
>   * CLS-D restates the logic in §2.2.2 but never discusses the order of limits, so do not say "the logic of CLS-D";
>   * f⁻ = −1.
> * **Theorem W:** the 2×2 criterion d ≤ H is correct and sharp (200 000 random tests). All J\* table entries are reproduced. The Gershgorin factor is ≤ 1.5 only if f₊₊ ≈ f₊₋ (otherwise up to 3; 1.50–1.53 on the samples). **Main qualification: d omits the k = 3 piece of the corrections. Including it moves J\* 50 → 64 at (3, 10, 10), 30 → 33, and 64 → 80**; the J = 50, 51 states fail at d/H ≈ 1.24. This must be stated as a named premise or used as the J\*. The cache stores the 4a/4b shapes with 1/m² instead of 1/m⁴ (documentation slip only).
> * **⟨TTO⟩ counts:** verified against DKKPS's TTTT.tex table and by an independent helicity count (ℓ ≤ 11).

## 0. Summary

### 0.1 Table

| # | claim | result | grade | effect on J\* |
|---|---|---|---|---|
| **PW1** | **Magic-block identity for every b-shifted block of the EE kernel:** `2η(1+r²) Ĝ¹_b(ν) = (1 + i/2ν) Ĝ⁰_b(ν+i) + (1 − i/2ν) Ĝ⁰_b(ν−i)` for b = 0…4, with **b-independent** coefficients (the CHMRS coefficients of b = 0) | holds identically | **[P]** to O(r⁸) (exact rational-function identity, sympy, l symbolic, all b = 0…4); **[DV]** to O(r⁴⁰) at 50 digits (residual ≤ 8·10⁻⁴⁸; also b = 7, 10); all orders [D] | — |
| **PW2** | **Kernel-level scalar de-shift of the raw EE functional.** Since `f = η(1+r²)/(2r)` in the H7 frame: `f K^EE_ν = N(ν)(r/4)[(1+i/2ν) E⁰_{ν+i} + (1−i/2ν) E⁰_{ν−i}]` for **every polarisation set, all orders in r**. Hence the raw image of **every** input (plane wave in any gauge, any polarisation/polynomial insertion, any heavy block) is exactly `𝒜_g[F]` (CHMRS's scalar de-shift), F even in μ and analytic in \|Im μ\| ≤ 3/2 | exact | **[P + DV]** given EESTEP's representation (end to end against the direct light-ray quadrature: 2·10⁻⁷–1.1·10⁻⁶, quadrature-limited); analyticity [DV] (would-be poles at ±i/2, ±3i/2 removable) + [D] | **Route W's stripped (R2) worst case cannot occur** |
| **PW3** | **PW-R1 as stated** (de-shift form w.r.t. the EE operator `𝒜 = 𝒜_g∘c` with f = F/c strip-analytic) **⇔ F(±i/2) = 0**, i.e. the Regge-spin-0 (l = 0, −1) value of the scalar-de-shifted image vanishes. It holds at leading Regge (F₀ ∝ c). **It is not protected:** the zeros of the EE/scalar ratio at ±i/2 are a property of the leading radial weights. Under the radial-weight shifts that O(1/m²) corrections produce they move: `|y|²` (side 13): zeros → ±3i/2, ±7i/2; `|x||y|`: → ±1.085i, ±2.515i; only `|x|²` (side 24) keeps them | **PW-R1 is refuted as a structural statement**; generically false at O(1/m²) | exact ratios **[P]** (sympy, v4_hclosed method); "generic" [D]; the actual O(1/m²) value of F(±i/2) for ⟨TTTT⟩ **[I]** (needs Route W's P2 coefficients) | — |
| **PW4** | **Exact cost of every raw correction:** `cost(𝒜_g[F]) = ∫ψ₂ F/c + D·F(i/2)`, **D = −9.98332·10⁻⁴** (Δ_g = 1000), −1.09754·10⁻² (300); D ≈ −990/Λ². D is **independent of the impact parameter** (spread 6·10⁻¹⁵ over t ∈ [0, 12]), because `P_t(±i/2) ≡ 1`. Equivalently: Lemma I holds for the scalar-de-shifted density σ_F with Ω₂ replaced by `Ω̃(t) = ∫ψ₂P_t/c + D`, **bounded by 1.0·10⁻³ for all t ≥ 0.03** and → D as t → ∞ (no e^{t/2} growth) | the non-de-shift part is a single **monopole** term | **[DV]** (formula to 10⁻¹²–4·10⁻⁸ on test mixtures; controls F = cP_t and F = P_t/(ν²+9/4) give defect 0) | replaces P3/R1 by a magnitude bound |
| **PW5** | **PW-Loc as stated** (no density at t < 1) | **false at O(1/m²)**: an index-shifted harmonic has the exact density `σ[Re P_{ν∓i}] = cosh β·δ_β − ½ d(cosh s)·1_{s<β}` (mass 0.271 per unit coefficient on t < 1); `σ[Re P_{ν∓2i}] = cosh 2β·δ_β − ½(3cosh s + cosh β) d(cosh s)·1_{s<β}`; a monopole (F(i/2) ≠ 0) carries the massless H₂ Green's-function density (infinite total variation). **But (Loc) is no longer needed:** costs are computed exactly per structure (PW4) | densities **[P]** (Abel inversion, closed form); costs [DV] | — |
| **TW′** | **Theorem W without P3.** With R1′ (exact c-ratio form, monopole included, sign-free, worst radial shift) in place of the stripped R2 | the scalar-shaped O(1/m²) vector is tolerated at ×24–29 for every sampled J ≥ 150 (stripped: ×1.8–7); **J\* = 50** at (Λ_amp, K_sh, K₄) = (3, 10, 10), **100 at margin ×2**, the same as with P3 assumed; W §9.3's window [50, 250] closes | [P given P0, P1, P2, P2′, P4] + [DV]; P2′ (multipliers in the radial-shift family) [D-sketch] | **J\* = 50 without R1** (was: no J\*, window [50, 250]) |
| **TL′** | **Theorem L with P_L reduced to magnitudes.** (R1) becomes an identity plus an explicit monopole term; (Loc) is replaced by exact per-structure costs. Monopole tolerance `Q_tol = min F_LR/(4\|D\|)` = 3.6·10⁻³ (J = 857), 1.7·10⁻³ (10⁴), 1.4·10⁻⁴ (10⁶). A scalar-shaped O(1/m²) correction (CHMRS (3.20) coefficients × Λ_amp, worst radial multiplier) is tolerated up to **Λ_amp = 1240 (J = 857), 555 (10⁴), 166 (10⁵), 43 (10⁶)** at the gap | J\* = 857 at the gap (unchanged; set by Lemma LR∞) | **[D given P_L′]** (magnitudes + remainder only) | none |

### 0.2 Bottom line

* **The plane-wave images were computed structurally to all orders.** The EE functional is, exactly and for every polarisation, a
  two-term finite difference in ν of a single de-shifted kernel family (PW1–PW2). So its plane-wave images beyond leading Regge are
  **always** of de-shift form with respect to the scalar operator 𝒜_g, at every order and in every gauge. The project's EE weight is
  built for `𝒜 = 𝒜_g∘c`. The two agree exactly when the de-shifted image vanishes at Regge spin l = 0.
* **Refutation (reported prominently).** PW-R1 in its literal form, "the images enter through 𝒜 exactly as the leading term does", is
  **not protected**. It needs `F(±i/2) = 0`. That holds at leading order only because of the leading radial weights.
  * Concretely [P], the plane-wave image of the leading EE kernel with an O(1/m²) insertion `|y|²` or `r = |x||y|` is **not** in EE
    de-shift form: its EE/scalar ratio is nonzero at μ = i/2 (0.914 and 0.500).
  * Whether the O(1/m²) terms of a full ⟨TTTT⟩ block cancel in the sum is open [I]; no mechanism for a cancellation was found.
  * In this language, Route W's "two extra channels" correspond to the O(1/m²) terms that do not carry the factor c(μ).
  * PW-Loc in its literal form is also false at O(1/m²): index-shifted harmonics have explicit density down to t = 0.
* **Why the refutation does not hurt.** The failure has an exact and benign form. It is a single number per state, the monopole
  `Q = F₁(i/2)`, with cost `2D·Q` independent of the impact parameter. It does not grow with β, unlike the stripped (R2) form that made R1
  load-bearing in Routes W and L (N-L3: tolerance ∝ 1/J). Every finite-order structure is then costed exactly.
* **Effect.**
  * **Theorem W:** the R1 premise P3 is no longer needed. **J\* = 50** at (Λ_amp, K_sh, K₄) = (3, 10, 10) (100 at margin ×2), the same
    as with P3 assumed. The window [50, 250] closes (§5).
  * **Theorem L:** J\* = 857 unchanged. The premise is reduced to magnitude bounds: tolerated Λ_amp = 1240 (J = 857), 555 (10⁴), 166
    (10⁵), **43 (10⁶)** at the gap (§6). Both theorems still rest on the unknown spinning O(1/m²)
  coefficients (Route W's P2) and on the beyond-all-orders premise. **Neither becomes unconditional.**

---

## 1. The question in this note's language

* **Inputs.** Route L (Proposition G, Corollary B) writes every heavy block as a positive average of plane waves `e^{p·w}`, with finitely
  many polynomial and polarisation insertions for ⟨TTO⟩⟨OTT⟩ (§2.5 there, [D-sketch]). Route W (W0) writes the same objects as
  weight-shifted scalar seeds. Either way the image of a heavy state is a linear superposition of the EE kernel's images of
  elementary inputs.
* **De-shift form (N2's R1, Route L's convention).** A raw correction is of de-shift form if it equals
  `𝒜[f] = (R₋^EE f(ν+i) + R₊^EE f(ν−i))/(R₋^EE + R₊^EE)` with f strip-analytic. Then the de-shift identity
  `Re∫X 𝒜[f] = ∫ψ₂ f` makes its k = 2 cost the ψ₂ cost of f, which Lemma I turns into a density integral.
* **The two open statements (Route L §5).**
  * (PW-R1): every plane-wave image beyond leading Regge is of this form.
  * (PW-Loc): its Lemma-I density is negligible at t < 1.
* **Known at leading Regge.**
  * `R₋^EE = R₋^sc c(ν+i)` and `R₊^EE = R₊^sc c(ν−i)` (EE_FUNCTIONAL §2, [P given H7 §18.1]). So `𝒜 = 𝒜_g ∘ c`, where 𝒜_g is CHMRS's
    scalar de-shift (same R^sc for every structure).
  * The leading raw image is `𝒜[P_β]`.

## 2. PW1–PW2: an exact, all-orders, polarisation-independent de-shift of the EE kernel

### 2.1 The identity (`pw1_magic.py` → `cache/pw1_magic.txt`)
EESTEP §1.2 reduces the exact EE kernel to `K^EE_ν = N(ν) r² Σ_{b=0}^{4} C_b(x; pols) v^{b/2} Ĝ¹_b(x; ν)`. Each piece is a **Dp = 1**
continuous-spin (0, b) block, and the `C_b` are ν-independent and carry all of the polarisation dependence. The raw integrand
contains the ν-independent factor `f = v′ − u′` (repro_engine header), and `f = η(1+r²)/(2r)` exactly in the H7 frame (checked:
6·10⁻¹⁵). For b = 0 CHMRS's magic-block identity (`app_superconvergence.tex`, eq. `magicblockidentity`, [Q]) states
`(v′−u′)G̃_{−1,ν} = …G̃_{−2,ν+i} … G̃_{−2,ν−i}`. Its leading Regge term is the Gegenbauer recursion for ηP_l.

**Claim (M_b).** For b = 0, …, 4 and all complex ν,

    2η(1+r²) Ĝ¹_b(r, η; ν) = (1 + i/2ν) Ĝ⁰_b(r, η; ν + i) + (1 − i/2ν) Ĝ⁰_b(r, η; ν − i),                          (2.1)

with **the same coefficients for every b**: `1 + i/2ν = 2l/(2l+1)` and `1 − i/2ν = 2(l+1)/(2l+1)`.

**Evidence.** Both sides are expanded in `r^n P_{l+m}(η)`, using the recursion of `ee_blocks.py` re-implemented with generic arithmetic.
* (a) **Exact.** With sympy, l symbolic, b = 0…4 and n ≤ 8, every coefficient of the residual is the zero rational function. That is
  30 coefficients for b = 0 and 54 for each b ≥ 1. **[P] to O(r⁸).**
* (b) At 50 digits, n ≤ 40, two random complex ν: max |residual|/|coefficient| ≤ 7.9·10⁻⁴⁸ for b = 0…4. For b = 7 and b = 10 it is
  ≤ 2·10⁻⁴⁷: the b-independence is not special to b ≤ 4. **[DV] to O(r⁴⁰).**
* (c) A least-squares fit that admits extra terms (Ĝ⁰_b(ν), rĜ¹_b(ν)) returns zero coefficients for them at ν = 0.7, 2.3 + 0.2i and 5
  (scratch test, reproduced by (a)).

All orders is [D]. The coefficients at each order are polynomial in b (EESTEP §1.3), the identity is CHMRS's at b = 0, and the
coefficients are fixed at r⁰.

### 2.2 Consequence: the raw EE kernel is a scalar de-shift, exactly
Multiplying (2.1) by `N r² C_b v^{b/2}/(4r)` and summing over b (the C_b are ν-independent):

    f K^EE_ν(x; pols) = N(ν) (r/4) [ (1 + i/2ν) E⁰_{ν+i}(x; pols) + (1 − i/2ν) E⁰_{ν−i}(x; pols) ].                    (2.2)

Since `(1+i/2ν)(γ₇(ν+i)/γ₇(ν))² = g(ν)` (CHMRS; v4_hclosed), the coefficients are CHMRS's g, g̃ up to a ν-only normalisation (N ∝ γ₇⁻²).
Pair (2.2) with any input I (plane wave in Route L's or the symmetric gauge, any polynomial or polarisation insertion, a full heavy
block). Then

    raw image(ν) = R₋^sc(ν) F(ν + i) + R₊^sc(ν) F(ν − i) =: 𝒜_g[F](ν),     F(μ) := ⟨(r/4) E⁰_μ, I⟩/κ_sc(μ),                  (2.3)

with κ_sc ∝ γ₇(μ)² the scalar leading normalisation, which is analytic and zero-free for |Im μ| < 13/2.

* **F is even in μ.** `Ĝ⁰_b(μ) = Ĝ⁰_b(−μ)` to 8·10⁻¹⁵ (pw1 (d)).
* **F is analytic in |Im μ| ≤ 3/2.**
  * The only zeros of the Dp = 0 Casimir denominators in that strip are at μ = ±i/2 (from (n, m) = (1, ±1), (2, ±2), (3, ±1)) and
    ±3i/2 (from (3, ±3)).
  * They are removable: |Ĝ⁰_b(μ₀ + ε)| is constant to 10⁻⁶ as ε goes 10⁻² → 10⁻⁶, for b = 0…4 (pw1 (d)) [DV].
  * Higher orders only raise the radial weights of the Laplace transforms, so no Γ-pole enters the strip [D].

**This proves R1 with respect to 𝒜_g** for every plane-wave image, to all orders and for every polarisation channel. That includes
the two channels that are invisible at leading order (Route W §9.1, rank 3). The R1 failure Route W's §3.4 test was looking for (a
non-𝒜 combination of index shifts) cannot occur. The de-shift is a property of the blocks, not of the polarisation tensors.

### 2.3 End-to-end check against the project's kernel code (`pw1` (c), `pw1b_e2e.py` → `cache/pw1b_e2e.txt`)
The test compares the directly integrated light-ray kernel `K^EE_ν` (`ee_quad`, 2-d quadrature) with the right side of (2.2) divided
by f, with E⁰ at complex ν ± i from the Dp = 0 block series. N(ν) is fixed at one point; there are four points
(ρ, B) ∈ {(0.25, 0.6), (0.55, 1.4), (0.85, 2.2), (0.4, 0.3)} and two random polarisation sets.

| ν | pol. set | spread of K/pred (nth = 128, ds = 0.02, N = 140) | refined (256, 0.01, 200) | identity on block sums, \|r²E¹/pred − 1\| |
|---|---|---|---|---|
| 2.5 | 0 | 5.3·10⁻⁶ | **2.9·10⁻⁷** | 1.6·10⁻¹⁰ |
| 6.0 | 0 | 4.3·10⁻³ | **7.0·10⁻⁷** | 3.9·10⁻¹⁰ |
| 2.5 | 1 | 1.7·10⁻⁴ | **2.2·10⁻⁷** | 2.9·10⁻¹⁰ |
| 6.0 | 1 | 7.6·10⁻³ | **1.1·10⁻⁶** | 3.0·10⁻¹⁰ |

The residual is the quadrature/series error of EESTEP's representation (its own end-to-end figure is 1.6·10⁻⁷–1.5·10⁻⁶). It falls ×10³
under refinement, while the identity itself holds on the block sums to 10⁻¹⁰.

## 3. PW3: what is left of PW-R1. The Regge-spin-0 condition

### 3.1 Reduction [P]
Write `𝒜_g[F] = 𝒜[F/c]` (because R^EE = R^sc c(ν ± i)). The raw image is therefore of **EE** de-shift form with strip-analytic
f = F/c **iff F/c is analytic in |Im μ| ≤ 1**.
* c has zeros at μ = ±i/2 (inside the strip) and ±3i/2 (outside |Im μ| ≤ 1). Poles at ±3i/2 do not matter: the control
  F = P_t/(ν² + 9/4) has defect ≤ 10⁻¹⁷ (pw3 (b)).
* So, F being even, **PW-R1 ⇔ F(i/2) = 0.**
* `P_{−½+iμ}(cosh t) = P₀ = 1` at μ = ±i/2 (l = 0, −1). So F(i/2) is the total mass of F's H₂ density, the **Regge-spin-0 (monopole)
  value** of the scalar-de-shifted image.
* At leading Regge `F₀ = c κ P_β` [P given H7 §18.1] and F₀(i/2) = 0.

### 3.2 The zero is not protected (`pw2_radial.py` → `cache/pw2_radial.txt`) [P]
A double-divergence (∇^a∇^b T_ab on H₂) structure would protect `c ∝ C(C+2)`, C = μ² + ¼, at every order. To test for it, the exact
leading-Regge heavy image (v4_hclosed/u2 method: equivariant gradient at one celestial point, all K-ratios rising factorials, no fits)
was recomputed with the radial weight of each side shifted, `|x|^c → |x|^{c+s_A}` (pair 24) and `|y|^c → |y|^{c+s_B}` (pair 13). These are
the shifts that O(1/m²) terms produce: block-side `|x|², |y|², x·y` insertions (Route W §3.1) and the kernel's own r¹ term (r = |x||y|).
For each shift, `c_s(μ) := R₋^EE/R₋^sc(μ − i)`, and the (HN) consistency `R₊^EE/R₊^sc = c_s(ν − i)` holds exactly in every case:

| shift (s_A, s_B) | c_s(μ) | zeros | c_s(i/2) |
|---|---|---|---|
| (0, 0) leading | 32(4μ²+1)(4μ²+9)/((4μ²+25)(4μ²+81)) | ±i/2, ±3i/2 | 0 |
| (2, 0) `|x|²` | same as (0, 0) | ±i/2, ±3i/2 | 0 |
| **(0, 2) `|y|²`** | 32(4μ²+9)(4μ²+49)/((4μ²+81)(4μ²+169)) | ±3i/2, ±7i/2 | **0.914** |
| **(1, 1) `r = |x||y|`** | 32(16μ⁴+120μ²+119)/((4μ²+49)(4μ²+121)) | ±1.0845i, ±2.5147i | **0.500** |
| (2, 2) | as (0, 2) | ±3i/2, ±7i/2 | 0.914 |
| (3, 3) | 32(16μ⁴+376μ²+1143)/((4μ²+121)(4μ²+225)) | ±1.894i, ±4.462i | 1.250 |

* **Reading.** The zeros of c come from the (13) light-ray side: the ratio of the spinning to the scalar light-ray Fourier–Laplace
  factors at the leading |y| weight. They are not a tensor-structural (double-divergence) property.
* So **PW-R1 in its literal form fails at O(1/m²) unless the O(1/m²) coefficients of the various terms conspire** to cancel
  `Σ_s c_s(i/2)·[term](i/2)`. No mechanism for such a cancellation was found.
* Route W's two extra channels (§9.1, rank 3) are, in this language, the O(1/m²) pieces whose EE/scalar ratio is c_s ≠ c.
* **Grade.** "Not protected": [P] (exact ratios). "Generically false for ⟨TTTT⟩": [D]. The actual monopole of a ⟨TTTT⟩ heavy block at
  O(1/m²) needs the spinning O(1/m²) coefficients (Route W's P2, §9.2 items 1–3): **[I]**.
* Polarisation tilts (non-weight-0 tensors at p̂) are a second O(1/m²) mechanism. They were not tested and could only add further
  c_s ≠ c terms.

## 4. PW4–PW5: the exact cost of a non-de-shift correction, and the densities

### 4.1 The monopole formula (`pw3_cost.py` → `cache/pw3_cost.txt`) [DV]
With the published X, ψ₂, R∓ (n2_common), compute for `F = P_t`:
`defect(t) := Re Σ X 𝒜[P_t/c] − Σ ψ₂ P_t/c` (n2 units).

| t | Ω₂(t) = Σψ₂P_t | Ω_c(t) = Σψ₂P_t/c | **Ω̃(t) = raw cost of 𝒜_g[P_t]** | defect(t) |
|---|---|---|---|---|
| 0 | 10.39 | 0.3256 | 0.3246 | −9.98332·10⁻⁴ |
| 0.01 | 0.443 | 1.48·10⁻² | 1.38·10⁻² | −9.98332·10⁻⁴ |
| 0.03 | 1.95·10⁻² | 1.44·10⁻³ | 4.44·10⁻⁴ | −9.98332·10⁻⁴ |
| 0.1 | 5.21·10⁻⁴ | 6.85·10⁻⁴ | −3.13·10⁻⁴ | −9.98332·10⁻⁴ |
| 0.3 | 2.67·10⁻⁵ | 5.02·10⁻⁴ | −4.96·10⁻⁴ | −9.98332·10⁻⁴ |
| 1 | 7.58·10⁻⁶ | 2.68·10⁻⁴ | −7.31·10⁻⁴ | −9.98332·10⁻⁴ |
| 2 | 5.39·10⁻⁶ | 1.19·10⁻⁴ | −8.80·10⁻⁴ | −9.98332·10⁻⁴ |
| 3 | 3.54·10⁻⁶ | 5.38·10⁻⁵ | −9.45·10⁻⁴ | −9.98332·10⁻⁴ |
| 6 | 7.03·10⁻⁷ | 6.31·10⁻⁶ | −9.92·10⁻⁴ | −9.98332·10⁻⁴ |
| 12 | (−2.7·10⁻⁷, grid-limited) | 1.54·10⁻⁷ | −9.98·10⁻⁴ | −9.98332·10⁻⁴ |

* **D(1000) = −9.983321·10⁻⁴.**
  * Spread over 18 values of t ∈ [0, 12]: 6·10⁻¹⁵.
  * Halving the grid changes D by 1·10⁻⁴ relative.
* **D(300) = −1.097545·10⁻²** (spread 2·10⁻¹²). Hence `D·Λ² = −995 / −978`: D ≈ −990/Λ², the scale of the floor's `A₁M₁/Λ²`.
* Controls:
  * F = cP_t (EE de-shift form) has defect ≤ 6·10⁻¹² (10⁻¹⁴–10⁻⁵ relative, at the floor of the large-t quadrature);
  * F = P_t/(ν² + 9/4) has defect ≤ 3·10⁻¹⁷ (1000) / 6·10⁻¹⁶ (300).
* Four random mixtures `F = Σ w_k(1 + a_kν²)P_{t_k}` satisfy the formula below to 10⁻⁹–4·10⁻⁸ relative.

**Statement.** For every even F with an H₂ density σ_F (`F = ∫σ_F(dt)P_t`, total variation finite):

    k = 2 cost of the raw correction 𝒜_g[F]  =  ∫ψ₂ F/c  +  D·F(i/2)  =  ∫σ_F(dt) Ω̃(t),      F(i/2) = ∫σ_F.              (4.1)

The k = 2 entry of M_Ω carries ring weight 2 (Theorem L / W conventions). So a monopole Q := F(i/2), in units where the leading image is
`F₀ = c P_β` (f₀ = P_β), costs **2D·Q = −2.0·10⁻³·Q**.

**Why it is benign.** The stripped (R2) form costs grow like e^{t/2} (raw `P_{ν±i}(cosh t)` growth; Route L N-L3; N2 §1.2 table).
Ω̃(t) instead is **bounded**: `|Ω̃(t)| ≤ 1.0·10⁻³` for all t ≥ 0.03, and Ω̃ → D as t → ∞. For t < 0.03 it rises like Ω₂(t)/32 (to
0.325 at t = 0), i.e. 1/32 of the EE-de-shift penalty of the same density (`c → 32` at large ν). The non-de-shift part therefore costs
an impact-parameter-independent amount per unit monopole. It is not penalised more at large β, which is exactly where R1 was
load-bearing.

### 4.2 PW-Loc: the densities, in closed form [P]
Abel inversion (`P_μ(cosh s) = (√2/π)∫₀^s cos(μt)(cosh s − cosh t)^{−½}dt`; `Re P_{ν∓ik}(cosh β)` has Abel profile `cosh(kt)K(β,t)`):

    σ[Re P_{ν∓i}(cosh β)]  = cosh β · δ_β  − ½ d(cosh s) · 1_{s<β},
    σ[Re P_{ν∓2i}(cosh β)] = cosh 2β · δ_β − ½(3 cosh s + cosh β) d(cosh s) · 1_{s<β}.                                       (4.2)

Checks: the total masses equal F(i/2) = (1 + η)/2 and (η + P₂(η))/2, as they must.
* **Ring classes** (A: ηP, B: ην²P, C1–C3: β-derivatives) have no density at t < β: ε_in = 0, as before.
* **Index-shifted harmonics** have a smooth negative density on [0, β]. On t < 1 it has mass ½(cosh 1 − 1) = **0.271** per unit
  coefficient for the ν∓i shift, and ½(3/2·sinh²1 + cosh β(cosh 1 − 1)) for ν∓2i.
  * With an O(1/m²) coefficient κ_sh (Route W's units) this is 0.271·κ_sh/m²: 9.5·10⁻⁸κ_sh at J = 857 and 1.4·10⁻¹⁰κ_sh at J = 10⁶.
  * Theorem L's old sufficient condition ε_in ≤ tol(0⁺)/4 ≈ 10⁻⁷ … 4·10⁻⁹ is then met only for κ_sh ≲ 1 … 30.
  * **PW-Loc in its literal form is false at O(1/m²)** whenever index shifts are present. The kernel's own r-expansion produces them
    (EESTEP §1.3; Route W W2). [D]
* **Monopole part.** If Q ≠ 0, then f = F/c has a pole at ±i/2. Its Lemma-I density contains the massless H₂ Green's function
  (−Δ)^{-1}: a density at every t, with infinite total variation. Lemma I in f-form is not available for it, and (4.1) replaces it.
* **Conclusion.** (Loc) is the wrong premise. Every finite-order structure has an exact cost by (4.1), in n2 units, with no density
  bound needed (§5, §6). Only the beyond-all-orders remainder still needs a density-type premise, as before (RD/H_AG).

## 5. Theorem W′: Theorem W without P3 (`pw4_theoremW.py` → `cache/pw4_theoremW.txt`)

**What changes.**
* Route W's P3 (de-shift form of the O(1/m²) ring terms) was [I], and its worst case R2 ("stripped": the raw correction equals S
  itself) set the J\* window [50, 250] (W §9.3). **PW2 excludes R2.**
* The worst case allowed by the exact identity is R1′: the term enters as `𝒜_g[ρ_s S]` with an EE/scalar multiplier ρ_s ≠ c. Its exact
  cost (4.1) is

      R1′(S, s) = 2[ Σψ₂ (c_s/c) S dν/Λ + D c_s(i/2) S(i/2) ] / m²,                                                            (5.1)

  with `S(i/2)` = η (A), −η/4 (B), 0 (C1, C2, C3, 4b), η²/16 (4a), (1+η)/2 (sh±1), (η + P₂(η))/2 (sh±2).
* Below, "R1′" is the sign-free maximum over the O(1/m²) shifts s ∈ {(1, 1), (0, 2), (3, 3)} of pw2; (3, 3) is included as a margin.

**Premise (P2′) [D-sketch].** The multipliers of the O(1/m²) terms lie in the radial-shift family {c, c_(1,1), c_(0,2)}. This is what
block-side `|x|², |y|², x·y` insertions and the kernel's r¹ term produce on the leading split-representation kernel. Polarisation tilts
are not covered (§8, caveat 3).

**Class-wise tolerances** K_tol = H/cost (minimum over the sampled masses; the same states, floors and structures as Route W §4.2, read
from `lsw/cache/costs_1000.json`). Entries are R1 (P3 assumed) / R2 (stripped; now excluded) / **R1′ (P3 dropped, exact)**:

| J | B | C2 | index-shifted (sum) | scalar-shaped vector ×1 |
|---|---|---|---|---|
| 26 | 18.8 / 18.6 / **18.8** | 0.70 / 0.70 / **0.70** | 2.8e4 / 5.8e4 / **2.0e4** | 3.26 / 3.26 / **3.27** |
| 50 | 96 / 71 / **104** | 5.5 / 5.6 / **5.5** | 4.4e4 / 1.1e4 / **3.5e3** | 23.2 / 22.0 / **23.5** |
| 100 | 406 / 48 / **114** | 69 / 38 / **30** | 1.8e4 / 3.1e3 / **1.1e3** | 140 / 23 / **40** |
| 150 | 989 / 61 / **75** | 69 / 6.5 / **42** | 4.6e3 / 1.1e3 / **625** | 225 / 17 / **31** |
| 200 | 1630 / 22 / **71** | 150 / 4.1 / **28** | 2.1e3 / 197 / **411** | 419 / 7.1 / **27** |
| 250 | 2480 / 7.6 / **77** | 292 / 7.5 / **24** | 1.3e3 / 89 / **280** | 682 / 3.5 / **27** |
| 300 | 3570 / 5.1 / **88** | 305 / 3.9 / **25** | 843 / 59 / **196** | 812 / 2.3 / **28** |
| 344 | 4500 / 5.0 / **101** | 2930 / 1.4 / **25** | 632 / 52 / **146** | 1520 / 1.8 / **29** |
| 301 (odd) | 7250 / 10 / **176** | 2160 / 7.5 / **47** | 1.7e3 / 118 / **389** | 2400 / 4.6 / **56** |

A (≥ 697), C1 (≥ 101) and C3 (≥ 2.1·10³ for J ≥ 100) are harmless in every form (`cache/pw4_theoremW.txt`).

**Reading.**
* The exact worst case costs **10–40× less** than the stripped (R2) form at J ≥ 200. It costs 2–60× more than full de-shift (R1) at
  J ≥ 100.
* **The scalar-shaped vector is now tolerated at ×24–29 for every sampled J ≥ 150**, against ×1.8–7 stripped. That was the obstruction
  of W §9.3.
* For the ring classes A and B the R1′ cost is dominated by the monopole term 2D·c_s(i/2)S(i/2)/m². For the index-shifted harmonics the
  monopole is the whole R1′ excess (the "monopole part" column of `pw4_theoremW.txt`).

**J\* (Theorem W′).** The criterion is Route W's (`d ≤ θH`), with R1′ for every class:

| (Λ_amp, K_sh, K₄) | θ | J\*: P3 assumed (W §4.3) | P3 dropped, stripped (W §9.3) | **P3 dropped, exact (this note)** |
|---|---|---|---|---|
| (1, 10, 1) | 1 | 30 | 30 | **30** |
| (1, 10, 1) | ½ | 38 | none (fails 300, 344) | **38** |
| **(3, 10, 10)** | **1** | **50** | none: window [50, 250] (fails 300, 344) | **50** |
| (3, 10, 10) | ½ | 100 | none | **100** |
| (10, 30, 10) | 1 | 64 | none | **64** |
| (10, 30, 10) | ½ | 100 | none | none (fails 300, 344) |
| (30, 100, 30) | 1 | 100 | none | none (ν² classes, J ≥ 120) |

**Theorem W′.** Theorem W (LARGESPIN_W §4.0) holds **with (P3) deleted**. Its place is taken by PW2 (exact) and (P2′) [D-sketch]. **J\* = 50**
at (Λ_amp, K_sh, K₄) = (3, 10, 10) with θ = 1, and **J\* = 100 at margin ×2**: the same values as with P3 assumed. The window
[50, 250] of W §9.3 closes. **Grade:** [P given P0, P1, P2, P2′, P4] + [DV] (sampled costs). P2 magnitudes stay [NE] (J ≤ 100) / [M]
(J ≥ 150).

## 6. Theorem L′ (`pw5_theoremL.py` → `cache/pw5_theoremL.txt`)

**Theorem L′.** Assume A1′, A2 and Δ_g = 1000. Let O be a heavy primary with **b̂ ≥ 1.0421Λ**, i.e. J ≥ J\*(m) = ⌈(1.0421m − 1)/2⌉
(**J\* = 857 at the twist edge**, unchanged). Write `M_Ω = diag(F_LR(β, u, σ_h)) + δM`.
* By PW2, the k = 2 part of every entry's beyond-LR raw correction is **exactly** `𝒜_g[F_ab]`, F_ab even and strip-analytic. There is no
  R1 premise.
* By (4.1) its cost is `2[∫ψ₂F_ab/c + D·F_ab(i/2)]`.

Assume premise **P_L′** (magnitudes only):
* (i) at O(1/m²), `F_ab = Σ_s ρ_s κ_s^{ab} S_s/m² + Σ_k ρ_k κ_{sh,k}^{ab} Re P_{ν∓ik}·(η/m²)^{k−1}/m²`, with the multipliers ρ in the
  radial-shift family {c, c_(1,1), c_(0,2)} (P2′ of §5), |κ_s| ≤ Λ_amp|κ_s^sc| and |κ_sh| ≤ K_sh. The index-shift coefficients are
  normalised as Route L's lightcone uniformity (L4) requires;
* (ii) the O(1/m⁴) and beyond-all-orders remainder costs at most a fraction ε_R of ½ min_h F_LR, e.g. via its scalar-de-shifted density
  and the bounded kernel Ω̃ (|Ω̃| ≤ 1.0·10⁻³ on t ≥ 0.03).

Then `M_Ω ⪰ 0` for every coupling vector whenever `Λ_amp/Λ_tol + K_sh/K_tol^{sh} + ε_R ≤ 1`, with the tolerances below (2×2 criterion of
Theorem L, half of the budget per entry). **[D given P_L′]**; the floor is [P + DV] (Lemma LR∞); the costs are [DV].

**Tolerances.**
* Gap states (and 1.3 m_gap), worse helicity, half of Theorem L's per-entry budget ½ min F_LR. All entries are maxima of the premise
  constants that the criterion tolerates.
* The monopole of a scalar-shaped O(1/m²) correction is bounded by `|Q| ≤ c*(|κ_A| + |κ_B|/4)η/m²`, with c\* = 0.914.
* The C-classes and every β-derivative have zero monopole (`∂_β P_{±i/2} = 0`).

| J | m/m_gap | β | min F_LR | Q_tol = minF_LR/(4\|D\|) | η/m² | **Λ_amp tolerated (monopole only)** | **Λ_amp tolerated, scalar vector, R1′** (R1) | K_tol sh±1, R1′ | K_tol C2 |
|---|---|---|---|---|---|---|---|---|---|
| 857 | 1 | 1.000 | 1.44·10⁻⁵ | 3.6·10⁻³ | 5.7·10⁻⁷ | 1003 | **1240** (6.5·10⁴) | 1.1·10⁴ | 1.3·10⁴ |
| 1000 | 1 | 1.100 | 1.40·10⁻⁵ | 3.5·10⁻³ | 5.6·10⁻⁷ | 992 | **1220** (1.1·10⁵) | 1.1·10⁴ | 1.8·10⁴ |
| 2000 | 1 | 1.611 | 1.20·10⁻⁵ | 3.0·10⁻³ | 5.2·10⁻⁷ | 912 | **1050** (1.3·10⁵) | 1.0·10⁴ | 3.5·10⁵ |
| 3000 | 1 | 1.948 | 1.08·10⁻⁵ | 2.7·10⁻³ | 5.1·10⁻⁷ | 832 | **925** (1.3·10⁵) | 9.5·10³ | 2.9·10⁵ |
| 10⁴ | 1 | 3.046 | 6.86·10⁻⁶ | 1.7·10⁻³ | 5.0·10⁻⁷ | 540 | **555** (2.9·10⁴) | 6.7·10³ | 6.1·10³ |
| 10⁵ | 1 | 5.305 | 2.08·10⁻⁶ | 5.2·10⁻⁴ | 5.0·10⁻⁷ | 164 | **166** (6.7·10⁴) | 2.2·10³ | 2.3·10⁴ |
| **10⁶** | 1 | 7.603 | 5.51·10⁻⁷ | 1.4·10⁻⁴ | 5.0·10⁻⁷ | 43 | **43** (7.6·10³) | 599 | 1.35·10³ |
| 10⁴ | 1.3 | 2.582 | 8.48·10⁻⁶ | 2.1·10⁻³ | 1.9·10⁻⁷ | 1790 | 1830 | 2.2·10⁴ | 9.4·10⁴ |
| 10⁶ | 1.3 | 7.079 | 7.48·10⁻⁷ | 1.9·10⁻⁴ | 1.8·10⁻⁷ | 168 | 169 | 2.3·10³ | 1.75·10⁴ |

Index shifts by 2 with the lightcone-consistent coefficient κη/m⁴ (L4) are tolerated at K ≥ 8·10⁵ (J = 10⁶) to 10¹⁰ (J = 857). With a
bare κ/m⁴·m² normalisation (Route W's units) they would be tolerated only at K ≈ 15 (10⁵) and 0.4 (10⁶). That normalisation would
violate L4's J-uniformity (relative size η²/m² unbounded), so it is not a physical worst case.

**Reading.**
* The monopole is the binding piece of R1′ at large J. Its cost 2|D||Q| does not decay with J, while the floor decays like M₁(β).
  η/m² → 1/(2(τ−2)²) is J-independent.
* The tolerated Λ_amp therefore falls slowly: ≈ 1000 at J = 857–3000, 555 at 10⁴, 166 at 10⁵, **43 at 10⁶** (gap; 168 at 1.3 m_gap).
* Against the J = 100 fits (Λ_amp ≈ 2–3, Route W C2) the margin is ×15–20 at J = 10⁶ and ×300–500 near J\*.
* Extrapolating `F_LR ∝ M₁(β) ~ βe^{−β/2}`, a J-independent monopole of scalar size × Λ_amp = 3 would exhaust the budget near
  J ~ 10⁹ (β ≈ 14). That is beyond the theorem's sampled range but inside the weight-realisation range (J ≲ 10¹²; LARGESPIN_L §3.3).
  **Theorem L′ as stated covers J\* ≤ J ≤ 10⁶ (sampled) for Λ_amp ≤ 43.** Beyond 10⁶ it needs Λ_amp ≲ 43·(F_LR(β)/F_LR(7.6)).

**Status of P_L′ for the graviton.**

| ingredient | before (Route L) | now |
|---|---|---|
| (R1) | [I]: PW-R1 | **identity** (PW2, [P to O(r⁸)] + [DV]) plus an explicit monopole term (PW4, [DV]); PW-R1 in literal form refuted as a structural property (PW3) |
| (Loc), finite orders | [I]: PW-Loc | **not needed**: exact per-structure costs (4.1); literal PW-Loc false at O(1/m²) (closed-form densities (4.2)) |
| magnitudes at O(1/m²) | [NE] (J ≤ 100), [M] | unchanged: [NE]/[M]. Requirement Λ_amp ≤ 43 (J ≤ 10⁶, gap) |
| beyond all orders | RD for EE [D given P_cyc]; not checked at J ≥ 150 | unchanged; may now be stated with the bounded kernel Ω̃ |

## 6.1 Suggested inventory changes (not applied)
* LARGESPIN_W §9.3 / UT inventory, "R1 at O(1/m²) [I]":
  * replace by "R1 w.r.t. 𝒜_g: identity (PW_IMAGES PW2). EE de-shift form fails at most through a Regge-spin-0 monopole with exact
    cost 2D·Q, D = −9.98·10⁻⁴ (PW4)";
  * Theorem W's J\* = 50 then holds without P3 (window closed).
* LARGESPIN_L Theorem L, premise P_L:
  * "(R1) [I], (Loc) [I]" → "(R1) identity + monopole bound Q ≤ Q_tol (§6); (Loc) replaced by exact per-structure costs at finite
    orders";
  * grade [D given P_L′] (magnitudes only); J\* = 857 unchanged; tolerated Λ_amp ≥ 43 on 857 ≤ J ≤ 10⁶.
* LARGESPIN_L N-L3 ("raw-form tolerance collapses ∝ 1/J"): this applies only to the stripped form, which PW2 excludes. The exact
  non-de-shift cost per unit monopole is J-independent.
* The EESTEP letter sentence can add: "the EE kernel is exactly a finite difference in ν of a single de-shifted kernel family, at all
  orders and for every polarisation."
* **Remark (not tested).** The monopole exists only because the EE weight is built for 𝒜_g∘c. A weight solving the Abel equation for
  𝒜_g with target 32ψ₂/c, if one exists, would have no defect. The (constant-in-t) solvability condition is the obstruction [I].

## 7. Numerical verification (against the project's own EE kernel code)

| # | check | code | result |
|---|---|---|---|
| V1 | f = v′ − u′ = η(1+r²)/(2r) in the H7 frame | pw1 (0) | 6·10⁻¹⁵ |
| V2 | (M_b) exact, b = 0…4, n ≤ 8, l symbolic | pw1 (a) | 0 nonzero residuals (30 + 4×54 coefficients) |
| V3 | (M_b) 50 digits, n ≤ 40, complex ν; b = 0…4, 7, 10 | pw1 (b) | ≤ 2·10⁻⁴⁷ relative |
| V4 | Dp = 0 blocks analytic at μ = ±i/2, ±3i/2; even in μ | pw1 (d) | removable (constant to 10⁻⁶ for ε = 10⁻²…10⁻⁶); evenness 8·10⁻¹⁵ |
| V5 | (2.2) end to end vs `ee_quad` direct light-ray quadrature, 2 polarisation sets, 4 points, ν = 0.7–6 | pw1 (c), pw1b | 2·10⁻⁷–1.1·10⁻⁶ (refined; quadrature-limited); identity on block sums 10⁻¹⁰ |
| V6 | leading c(μ) and (HN) consistency reproduced (shift (0,0)) | pw2 | exact (sympy), equals EE_FUNCTIONAL §2 |
| V7 | defect of F = P_t/c is t-independent; controls F = cP_t, P_t/(ν²+9/4) | pw3 (a), (b) | spread 6·10⁻¹⁵; controls ≤ 6·10⁻¹² |
| V8 | (4.1) on random mixtures | pw3 (c) | 10⁻⁹–4·10⁻⁸ |
| V9 | n2-grid Ω₂(β) vs `lsl_core.Omega`/Lemma LR∞ form | pw5 | 0.1–0.4 % (trapezoid on the 0.04 grid) |
| V10 | densities (4.2) | pw6 | see `cache/pw6_density.txt` |
| V11 | W's R1/R2 costs reproduced from `lsw/cache/costs_1000.json` (same states, floors, structures) | pw4 | by construction (read-only reuse) |

## 8. Caveats

1. **(M_b) for all orders is [D].** It is proved exactly to O(r⁸), which covers the r⁰ and first subleading coefficients of the
   de-shifted kernel even though the C_b have r⁻⁸ poles. It is checked to O(r⁴⁰). A structural all-orders proof (a weight-shifting
   or Casimir argument for (0, b) blocks) was not written down. CHMRS's b = 0 identity is [Q].
2. **PW2 rests on EESTEP's representation K^EE = N r² E¹** ([P] reduction, [D + DV] identification of the Ψ^{(b)} as (0, b) blocks).
3. **"Generic failure" of PW-R1 (PW3) is shown on the radial-shift family only.** The physical O(1/m²) monopole for ⟨TTTT⟩ was not
   computed; it is a definite linear combination of the c_s(i/2) with Route W's unknown P2 coefficients. A cancellation cannot be
   excluded, but no mechanism for one was found. Polarisation-tilt terms were not tested.
4. **D is computed on the published Abel weight** (0.04 grid): grid halving moves it by 10⁻⁴ relative, and the formula (4.1) holds to
   10⁻⁸. The underlying contour-shift mechanism (residues of 1/c(ν ± i) at ν = ∓i/2 picking `P_t(±i/2) = 1`) is [D]; the value of D is
   [DV].
5. **The k = 3 piece (Ĉ₃) is not an EE functional** and is untouched. As in Route L it is counted raw with |ψ₃|.
6. **Theorems W′ and L′ still need the magnitudes** (Λ_amp, K_sh, K₄, and the remainder). This note removes a *structural* premise
   (R1, and (Loc) at finite orders); it does not compute the spinning coefficients.
7. **Δ_g = 1000** for the theorem numbers; D(300) is quoted for scaling only.

## 9. Reproduce

    cd ads4/code/pwimg && bash run_all.sh       # ≈ 1.3 CPU-h, one process at a time

| file | content |
|---|---|
| `pw_common.py` | (0, b) blocks at complex μ, E^{Dp}_μ, c(μ); read-only imports of eestep |
| `pw1_magic.py` | (M_b): exact (sympy) and 50-digit checks, analyticity/evenness, end-to-end vs `ee_quad` |
| `pw1b_e2e.py` | refined end-to-end check |
| `pw2_radial.py` | exact EE/scalar ratio c_s under radial-weight shifts (v4_hclosed/u2 method) |
| `pw3_cost.py` | Ω₂, Ω_c, Ω̃, the monopole constant D (Δ_g = 1000, 300), formula (4.1) |
| `pw4_theoremW.py` | R1′ costs on Route W's states, K_tol, J\* without P3 |
| `pw5_theoremL.py` | monopole tolerance, exact per-structure costs and tolerated Λ_amp, J = 857 … 10⁶ |
| `pw6_density.py` | check of the closed-form densities (4.2) |
