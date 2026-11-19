# Foundations of `Ω^{EE} = ∫ψ₂Ĉ′_EE + ∫ψ₃Ĉ₃`: admissibility of `EE`, light side and swappability, Lemma N, the canonical gauge, `Ĉ₃`

*19 Nov 2026, EE-foundations agent. Unrefereed, not committed. No existing document was edited (`EE_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, `H0_SPINNING_FUNCTIONAL.md`, `RESULT_D.md` were read only). New code: `ads4/code/spinfunc/x1_admissibility.py`,
`x2_gauge_EE.py`. Outputs: `spinfunc/cache/x1_admissibility.txt`, `x2_gauge_EE.txt` (+ `.log`), `x2_chi_EE_can_{100,300,1000}.pkl`.
Light numerics only (one core, under 3 min in total). Project modules (`u2`, `u3`, `v4_chi`, `v4_hclosed`, `g4_images`) are imported read-only.*

*Source: CLS-D = Chang–Landau–Simmons-Duffin, arXiv:2311.04271, LaTeX `draft.tex` from the session scratchpad
(`…/scratchpad/dict/src/2311.04271/`), read as data. Equation labels are CLS-D's own.*

Grades: **[P]** proved (here this includes exact algebra); **[P given X]**; **[Q]** quoted; **[D]** derived here, not independently
checked; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]** open.

---

## 0. Verdict

| item | before (EE_FUNCTIONAL) | now | what remains |
|---|---|---|---|
| **`EE` admissible as a `j′ = 0` dual (2,4) structure** | [D-sketch] | **[D]**. Weights, transversality and membership in CLS-D's general structure space are [P] (§1.2–§1.3). CLS-D's printed `j′ = 0` list is a basis of the bulk-point (saddle) quotient only (§1.3). Vanishing and convergence of `Ψ_EE` do not depend on the dual structure (§1.4). | The §1.4 reduction is [D]. Convergence inherits K1 (KKSDZ, with its loophole), as for `q²`. |
| H7 for `Ĉ′_EE` (graviton pole exact) | [P given Lemma R, Q3, Q4, D2, §18.1 (D)] | **same grade, with a smaller input.** (HN) now holds in its original form: the scalar CHMRS magic-block kernel up to the ν-only factor `c(ν)`. The (2,4) side of the kernel is exactly 16 × the scalar's. Lemma V's polarisation referral is no longer needed. | Lemma R's R4; D2; two of the three [D] inputs of H7 §18.1 |
| I14(a) for `EE` | (assumed as for `q²`) | **[P given K1, K2]**, same as `q²`. The K2 rows are independent of the dual structure (§1.4). | K1 loophole; K2 endpoints [D] |
| I14(b), K2 at `η → ∞` | (assumed) | **[P given K1, K1′, K2, K3]**. Proposition K applies verbatim, because the de-shifted heavy kernel is `∫ψ₂P_l`. Spinning components [D], simpler than for `q²` (§2.4). | K1′ uniformity; K3 [D] |
| Lemma N for `EE`'s 2-dim null space | [P given NE regularity] (asserted) | **[P given NE regularity] + [DV]**: annihilation `≤ 10⁻¹⁵`, residue formula to `10⁻¹⁶`, `X_h = O(ν²)` (§2.5) | — |
| **Canonical gauge** | [I] | **The CHMRS class selects a unique weight [NE + D].** The conditions are `F(−i/2) = F′(−i/2) = 0` (2 complex conditions on 2 complex parameters). The 2×2 matrix is non-singular (singular values `0.223`, `1.45·10⁻⁴`, the same at every `Δ_g`). The solution is grid-stable. **Adopted; the gauge-choice caveat is removed at the level of the definition.** **But the canonical weight is large at `ν ≲ 20`:** `X_can/ψ₂ ≈ −185` at `ν = 5` and `+65` at `ν = 10`. The `L¹` norm of `X_can − X_min` is `9.4 / 1.0 / 0.09 × ‖ψ₂‖₁` at `Δ_g = 100 / 300 / 1000`. | The heavy side beyond leading Regge at `ν ≲ 20` has to be re-evaluated in this gauge [I] (§3.4) |
| `Ĉ₃` unchanged | — | **confirmed** (§4). Its definition, normalisation, light side, convergence and flat image use only its own structure `o_H·q`. The only contacts with `Ĉ₂` are the flat-limit relative normalisation (unchanged, since `EE ≡ q²` at the crystal with ratio exactly 1) and `δk(0)`. | — |

---

## 1. Task 1: is `EE` a legitimate dual (2,4) structure?

`EE := [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]² = E²` and `q² = [(2H₂₀H₄₀ + V_{0,42}²H₂₄)/(X₂₀X₄₀)]²`, with CLS-D's building blocks
(eq. `CFT_building_blocks_VH`).

### 1.1 What "admissible" has to mean

CLS-D define the spinning functional (eq. `spinningfunctional_definition`) by integrating the light-ray commutator
`⟨O₄[L[O₃](x₀,z₃), L[O₁](x₀,z₁)]O₂⟩` against two ingredients:
* a dual structure `⟨O₄ L[O] O₂⟩^{(a)−1}`, where `O` has `Δ = 3/2 + iν`, `J = J₁+J₃−1 = 3` and transverse rep `j′`;
* a celestial structure.

The dual structure enters three things:
* (i) conformal invariance of the integrand (weights and transversality);
* (ii) which linear combination of the coefficient functions `C_{ab}(ν)` the functional measures (CLS-D, line after
  eq. `lambdaprime_selectionrule`);
* (iii) possibly, convergence.

The vanishing on physical correlators is pointwise: the integrand contains the commutator, which vanishes. So it holds for any dual
kernel once the integral converges. Admissibility therefore reduces to (i) and convergence. Item (ii) only says what the sum rule
measures.

### 1.2 Weights and transversality [P] (`x1` (A))

* The homogeneity degrees of `q²` and `EE` in each embedding variable agree exactly. They are `(X₀, Z₀, X₂, Z₂, X₄, Z₄) = (0, 4, 0, 2, 0, 2)`
  relative to the common scalar dual prefactor (round-off `10⁻¹⁵`, in `d = 3` Lorentzian and `d = 4` Euclidean).
* Both are invariant under `Z_i → Z_i + αX_i` to `≤ 4·10⁻¹⁵`.
* This is exact algebra: `H_ij` and `V_{i,jk}` are transverse; `V_{2,04}V_{4,02}` and `H₂₄` have the same degrees; and `V₀²` matches `H₂₀H₄₀`.
* **So `EE` has exactly the weights of `q²`.** Whatever scalar dual prefactor multiplies `q²/V₀⁴` (H7 §18.1, `l6`) also makes
  `EE/V₀⁴` an invariant integrand.

### 1.3 Span: CLS-D's general family vs their printed list [P + Q]

* **CLS-D's general family.** The labels `(a)` of the dual structures range over **all** continuous-spin structures of
  `⟨O₂O₄O⟩`, counted by eq. `dual_spinning_counting`. They are defined as the dual basis, under the Lorentzian pairing, of a basis of
  `⟨O₄L[O]O₂⟩` structures (eq. `dual_struct_definition_spinning`).
  * In `d = 3` with two spin-2 operators the count is `25 = 13 + 12`: 13 structures for `j′ = 0` (even under the transverse `O(1)`, i.e.
    parity-even) and 12 for `j′ = 1`.
  * **Check [DV, x1 (B)].** The 14 parity-even monomials `V₂^aV₄^bV₀^cH₂₄^{n₁}H₂₀^{n₂}H₄₀^{n₃}` (with `a+n₁+n₂ = b+n₁+n₃ = 2` and a common
    prefactor) have numerical rank **13 in `d = 3`** (one Gram relation, printed in `x1_admissibility.txt`) and 14 in `d = 4`.
  * `EE = V₀^J(H₂₄² − 4H₂₄V₂V₄ + 4V₂²V₄²)/(…)` is a combination of three of these monomials, so it lies in the 13-dimensional `j′ = 0` space
    [P].
* **The printed list `H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ` is a basis choice at the saddle, not a classification [Q + DV].**
  * CLS-D derive it from eq. `kernel_identity_saddle`, `−V_{2,40}V_{0,24} = H₂₀` at the scattering crystal for polarisations with zero time
    component (large-Δ conservation). In their words: *"This identity allows us to just remove all `V_{2,40}, V_{4,02}` structures when
    writing down the basis for the spinning dual structures."*
  * Reproduced (`x1` (C)): the identity holds to `2·10⁻¹⁶` at `x = 0.2, 0.5, 0.8`, both `S⁰` points, for 3 polarisations, for both 2 and 4.
  * Their count "`min(J₂,J₄)+1`" checks the list against eq. `spinning_dual_counting`, which is the large-ν count **for conserved
    externals**. Away from the saddle (Regge, `x₀ = 0`) the `V`-structures are independent; there `EE` and `q²` differ.
  * In `d = 3` that count is `2`, not `3`: the `SO(2)` spin-2 rep restricted to `O(1)` splits as `1 + 1`, so `j′ = 0` has 2. The printed
    3-element list is therefore **overcomplete at the `d = 3` saddle by one**: one combination is null there. This is H0 §12's null
    direction (`a = 0` members of `aP² + bPH + cH²`) and EE_FUNCTIONAL's null `EH` [D, consistent with both].
* **Rank tests [P/DV].**
  * `q²` lies in the printed span with coefficients exactly `(4, 4, 1)` (residual `5·10⁻¹⁶`).
  * `EE` does not (residual `0.9–1.0`), but lies in the full span (`10⁻¹⁵`).
  * At the saddle `EE ≡ q²` (EE_FUNCTIONAL §1.2, [P]). So **on the saddle quotient `EE` is the printed element `(4, 4, 1)`**.
* **Geometric form [DV, x1 (D)].** `e₂₄ := (H₂₄ − 2V_{2,04}V_{4,02})/X₂₄ = (I(x₀₂)z₂)·(I(x₀₄)z₄)` to `3·10⁻¹³` (`d = 3` both signatures,
  `d = 4`): both polarisations are inversion-transported to the light-ray point `x₀` and contracted there. So
  `E = V₀²X₂₄e₂₄/(X₂₀X₄₀)` is the "`z₂·z₄` at `x₀`" structure, of which EE_FUNCTIONAL's `x₀ = 0` reduction is the special case.

### 1.4 Validity of the sum rule and convergence do not depend on the dual structure [D]

* **Complete gauge fixing in `d = 3`.** Use conformal symmetry to fix `(x₂, x₄, x₀, [z₀])` in the region "`2 > 4`, `0` spacelike to
  `2, 4`". The parameter count is `3 + 3 + 3 + 1 = 10 = dim SO(3,2)`, and for generic `[z₀]` the residual stabiliser is finite.
* In this gauge the dual structure is a **constant polarisation tensor** `S(z₂, z₄)` (and the Faddeev–Popov factor is constant). The
  remaining integral is the celestially smeared light-ray product at the fixed point `x₀`:
  `𝒲 := ∫Dz₁Dz₃⟨P̃₁P̃₃P_δ⟩ ∫dα₁dα₃ f(u′,v′)⟨O₄[O₃,O₁]O₂⟩`.
  Equivalently, `Ψ_S = ⟨S, 𝒲⟩` is a finite-dimensional pairing with the 3-point object `𝒲`. Since `f` is conformally invariant, `𝒲` is a
  covariant structure (CLS-D: `Ψ ≡ C_{ab}(ν)`).
* **Consequences.**
  * (a) `Ψ_S[G] = 0` for every `S` of the right weights once `𝒲` converges.
  * (b) `|S·G| ≤ ‖S‖ Σ_components |G_comp|`, with `O₂`, `O₄` at fixed points. Absolute convergence for any `S` therefore follows from
    absolute convergence of the light-ray integrals of every polarisation component. These are exactly the KKSDZ conditions used in H7 §7.1:
    * the Regge corner `J₀ < 2` (CLS-D §3.5.2: set by `f⁺_{2;2}` alone);
    * the `1 ≈ 3` endpoint;
    * the light-cone edges.
  * (c) Absolute integrability of the contracted integrand is invariant under change of variables. It is therefore the same in the CHMRS
    frame that H7 uses.
* **Hence `EE` needs exactly the same convergence conditions as `q²`: `J₀ < 2`, with the proven `J₀, J₀^{dDisc} ≤ 1` (KKSDZ) [Q].**
* **The `[TT]` double-trace zeros** come from the block's commutator factor `2 sin²(π(τ−10)/2)` (CLS-D eq.
  `heavyaction_definition_generalrep`). They are independent of `S`. So the tree-level light side is again the `T` block alone, with the
  same parity-odd caveat as `q²` (H0 §15.1).
* **A caveat on the CHMRS-frame bookkeeping [D + DV].** Contracted absolute convergence is frame-independent; component-wise bounds are
  not.
  * Near `x₀₂² = 0`, `q ~ X₀₂⁻¹` but `E ~ X₀₂⁻²` (x1 (D): fitted `−0.9994`, `−1.9991`). The singular locus is the same; the order is higher.
  * In the CHMRS frame (`x₁ = −y`, `x₃ = y`, `x₂ = −x₄ = x/x²`, `x₀ = |y|n̂_θ`) one has
    `x₀₂²|x|² = r² + 2r sinh β cos θ − 1 ≤ (ρ−1)(ρ̄+1) < 0` for `ρ < 1`. So this locus is reached only at the light-cone edge `ρ → 1`
    (`x₁₂² → 0`), at `θ = 0`. It never touches the Regge corner.
  * There, component-wise bounds on the `EE` kernel are worse than for `q²` by two powers of `x₀₂²`. Bounds of the type used in Lemma S
    should therefore be applied in the fixed-`(2,4,0)` gauge, or to the contracted integrand. Lemma R uses only the contracted `𝒦`, so it is
    unaffected.

### 1.5 Grade

| statement | grade |
|---|---|
| `EE` has the weights of `q²` and is transverse | [P] (x1 (A)) |
| `EE` lies in CLS-D's general `j′ = 0` structure space (13-dim in `d = 3`); not in the printed 3-element list | [P] (x1 (B)) |
| the printed list is a saddle-quotient basis for conserved externals (CLS-D `kernel_identity_saddle`); overcomplete by 1 in `d = 3` | [Q + DV] / [D] |
| `Ψ_EE[G] = 0` on the non-perturbative `⟨TTTT⟩`, with the same convergence conditions as `q²` (`J₀ < 2`) | **[P given K1, K2 and the §1.4 reduction (D)]** |
| **overall: `EE` is an admissible dual structure** | **[D]** (was [D-sketch]) |

---

## 2. Task 2: H7, I14 and Lemma N for `EE`

### 2.1 The leading-Regge kernel: (2,4) side = scalar (`x1` (E)) [P given H7 §18.1]

* In the exact plane-wave transform (v4_hclosed/u2 method), the (2,4)-side covector coefficients satisfy
  **`c₁^{EE}/c₁^{sc} = c₂^{EE}/c₂^{sc} = 16`** (sympy, exact).
* The reason: at `x₀ = 0`, `16·EE/(z₀·x̂)⁴ = (z₂′·z₄′)²`, which depends on neither `x̂` nor `z₀`.
* Hence `h₋^{EE}/h₋^{sc} = 16 d₁ᵀ/d₁^{sc} = c(ν+i)` and `h₊^{EE}/h₊^{sc} = c(ν−i)` exactly (both differences `= 0` symbolically), with
  `c(ν) = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))`. So **`c` is entirely a light-ray (1,3)-side quantity**, and that side is common to all
  `j′ = 0` structures.
* `(z₂′·z₄′)²` contracted with a heavy tensor `τ₂₄` is an `O(2,1)`-invariant trace. It is therefore the same under any common Lorentz
  referral of the two polarisations. **Lemma V's inversion referral, one of the three [D] inputs of H7 §18.1, is not needed for `EE`.**
* For the same reason the (2,4)-side flip vanishes manifestly: `(ε_h·ε_h)² = 0` for null `ε`.

### 2.2 H7 (graviton pole exact)

| ingredient | `q²` (H7 §12–§20) | `EE` |
|---|---|---|
| Lemma R | [P given R1–R5] | same. R1, R2, R4, R5 do not involve `S`. **R3 [D]**: `EE` is rational in `x_ij²`, `z·x` times the same scalar-dual powers, with its singularities on the **same locus** (`x₀₂² = 0`, `x₀₄² = 0`) as `q²` (only the order differs, §1.4); dDisc commutes with it for the same reason. Lemma R uses only the contracted `𝒦`. |
| D3 (heavy vs light plane-wave bookkeeping) | [P] | [P], structure-independent; it needs only `a = 7`, which holds for `EE` (identical weights, §1.2; `s(∞) = 23/2`, `c(∞) = 32`) |
| `h₀ ≡ 0`, helicity blindness | [P given §18.1] | [P given §18.1] (the §18.3 proof holds for any (2,4) structure) |
| **(HN)** | only via the spinning de-shift, by construction; fails in the scalar-`g` form (§19) | **holds in its original form [P given §18.1]**: `h∓ = (g, g̃)·c(ν±i)`, so `Ĉ′_EE = 𝒜_h⁻¹[Ψ_EE] = c(ν)⁻¹𝒜_g⁻¹[Ψ_EE]`. This is CHMRS's own de-shift with the ν-only normalisation `N̂ = c`, i.e. H0 §24.1's "same combination for every correlator" is exact for `EE` with `N_T = c` (the special `N_T*` of H7 §19.1 is the explicit rational `c`) |
| §18.1 kernel [D] inputs | Ψ⁺ vs −(1↔3) reading; ν-only constants; Lemma V referral | the first two (the second is irrelevant: everything is normalisation-free); **the referral is not needed** (§2.1) |
| smeared pole identity V3 | `≤ 8·10⁻¹²` (v4) | `≤ 1.1·10⁻¹⁴` (u3, min-norm); **`≤ 3·10⁻¹³` in the canonical gauge** (x2) |

**Grade of H7 for `Ĉ′_EE`: [P given Lemma R, Q3, Q4, D2, §18.1 (D; two inputs)]**, so `−Ĉ′_EE|_E = 16πG/(ν²+9/4)` exactly. This is the
same grade as `q²` with one [D] input fewer. (HN) is satisfied in the scalar-CHMRS form, which is the form H7 §5.2 E2 originally asked for.
**`A = 3022.34`, `B_AdS = −6552.66`: unchanged.**

### 2.3 I14(a)

**[P given K1, K2]**, the same as `q²`. By §1.4 every K2 row of H7 §7.1 is a light-ray (KKSDZ) condition independent of `S`:
* the Regge corner, via the CLS-D count with `f⁺_{2;2}`;
* the `1 ≈ 3` endpoint;
* the light-cone edges.

### 2.4 I14(b) and K2 at `η → ∞`

* **Lemma S.** Same proof (per-family Rindler positivity, Cauchy–Schwarz), applied in the fixed-`(2,4,0)` gauge, where `|ℋ|` is the
  constant tensor `S`. Grade **[P given K1, K2, K3]**, with K3 [D].
* **K2 at `η → ∞` (Proposition K, H7 §14).**
  * The de-shifted, heavy-normalised total kernel of `∫ψ₂Ĉ′_EE` is `∫X𝒦^{EE}_ν(η)dν = ∫ψ₂P_l(η)dν`. This is exact by the Abel identity
    (u3 residual `≤ 3·10⁻¹⁵`; x2: `≤ 2.5·10⁻¹³` in the canonical gauge). It is gauge-independent, because the null vectors annihilate
    every `Â_h`-image (§2.5).
  * This is literally the same function `𝒦_w`, `w = ψ₂`, as for `q²`. **Proposition K applies verbatim [P]**, with decay
    `η^{−½}/log η`; hence **[P given K1′]** for the Fubini corollary.
  * The dDisc-space weight changes only `M₀, M₂, c₄`, not the rate.
* **Spinning components [D].** As H7 §14, but simpler. At leading Regge the (2,4) side of the `EE` kernel is the scalar's times a
  **constant** tensor (§2.1), so its only tensor growth comes from the light-ray (1,3) side, which is identical to `q²`'s.
* **Extra caution for `EE`.** The forward-normalised local blocks `½[(1±u_EE)P_{l∓1}]` are singular at the real zeros `1.17`, `6.64` of
  `h₋+h₊`, not just growing like `η^{3/2}`. The swap must therefore use the total kernel in the χ-normalisation, as u3 does.

### 2.5 Lemma N for `EE` (`x2` (L)) [P given NE regularity] + [DV]

* **Null vectors.** u3's two vectors (`dν = 0.025`), re-projected:
  * `|AY|/|Y| ≤ 3·10⁻¹⁵`;
  * `|X|` peaks at `ν = 0.56` and `4.46`;
  * **`X_h = O(ν²)`** (`X/ν² = 10.66, 1.90` at the first two nodes);
  * `|X_h(30)| ≤ 6·10⁻⁹`, `|X_h(50)| ≤ 4·10⁻¹⁶`.
* **(ii) Annihilation of `Â_h[g]`** for even, strip-analytic `g`:
  * `1/(ν²+a²)`, `a = 1.5, 3`; `e^{−ν²/4}`, `e^{−ν²/25}`; `cos(0.5ν)`, `cos(2ν)`; the graviton pole `1/(ν²+9/4)`: `|I|/scale ≤ 1.1·10⁻¹⁵`;
  * `a = 1.1` (closest to the strip edge): `6·10⁻¹²`;
  * Mehler–Fock kernel `P_l(η)`, `η = 1.5, 3`: `≤ 3·10⁻⁸` (truncated at `ν ≤ 25`).
* **(iii) Residue formula** `I_h = −(π/a)F_h(−i(1−a))` for `a = 0.3, 0.6, 0.9`: agreement to **`≤ 8·10⁻¹⁶`** relative (`F_h` from the
  representation of §3.1).
* **Generic functions are not annihilated**: `ν²e^{−ν²/25}` and `1/(ν²+9/4)` paired directly with `X_h` give `0.4–0.7`.
* So the two null directions are invisible to:
  * the leading-Regge action;
  * the graviton pole;
  * the `|ĝ₃|²` term (given its C-level strip analyticity);
  * the flat limit.

  They are visible only through the part of `Ψ_EE[state]` that is not `Â_h` of a strip-analytic density, which is H7 §24.0 item 2's [I].

---

## 3. Task 3: the canonical (CHMRS-class) gauge

### 3.1 The class and the conditions [D]

* `R₊^{EE}(ν) = R₋^{EE}(−ν)` and `N_com` is even, so `h₊(ν) = h₋(−ν)` and `ρ(ν)ρ(−ν) = 1`. Hence `χ := f/h₋` (with `f = χh₋` on `ν > 0`
  and `χh₊` on `ν < 0`) is single-valued and even on `ℝ`.
* `χ` is analytic on the closed strip `|Im ν| ≤ 1` iff:
  * `F/h₋` is analytic in `−1 ≤ Im ν < 0` (`F` is the continuation of `f`, analytic there for every Abel solution);
  * the upper half follows by evenness;
  * the exception is the `ν = 0` kink that `ψ₂ ∝ |ν|` forces in every class, the "allowed pole term", as for the scalar.
* In `−1 ≤ Im ν < 0`, `h₋^{EE}` has **exactly one zero, a double zero at `−i/2`** (`g`'s zero times `c(ν+i)`'s zero), and no pole.
* **The class is therefore `F(−i/2) = F′(−i/2) = 0`: 2 complex conditions on the 2 complex null coefficients.**
* (`N_com` is zero-free in the strip, so using `R₋` instead of `h₋` changes nothing.)

**Evaluating `F` without exponential error amplification [DV].** For `t > 0` the Abel equation gives
`f̂(t) = 2e^{−t}Φ(t) − e^{−2t}f̂(−t)`. For `−1 < Im z < 0` this yields the exact real-axis representation

`F(z) = (1/2π)[∫_ℝ f(ν)(i/(ν−z) − 1/(2−iz−iν))dν + ∫₀^∞ψ₂(ν)(1/(1−iz−iν) + 1/(1−iz+iν))dν]`,

and its `z`-derivative; the `ψ₂` term is absent for null vectors. Checks:
* synthetic `f`: `F`, `F′` to `≤ 9·10⁻¹⁶`;
* it reproduces H7 §24.2's `|F_h(−i/2)| = 3.1·10⁻⁶ / 5.1·10⁻³ / 0.108` and `|F_h(−0.2928i)| = 3.5·10⁻⁶ / 2.8·10⁻³ / 0.150` (q² null vectors);
* it reproduces §24.6's scalar `|F_h(−i/2)|/∫|f_h| = 1.141326·10⁻³` (H7: `1.141327·10⁻³`) at `dν = 0.05` and `0.025`.

### 3.2 Result [NE] (`x2` (G); null vectors normalised to `max|X_h| = max ψ₂`, the H7 §24 convention)

| `Δ_g` | `sv(M)/max ψ` | `c₁` | `c₂` | `|F|, |F′|` after | grid `dν → dν/2`: `|Δc|`, `max|ΔX_can|/max ψ` |
|---|---|---|---|---|---|
| 100 | `0.2228, 1.454·10⁻⁴` | `−13.78 + 2.89i` | `24.94 + 76.09i` | `2·10⁻¹⁶, 1·10⁻¹⁵` | `2·10⁻⁴, 2·10⁻³`; `5·10⁻⁵` |
| 300 | `0.2229, 1.454·10⁻⁴` | `−4.596 + 0.964i` | `8.32 + 25.41i` | `3·10⁻¹⁶` | `8·10⁻⁴, 3·10⁻³`; `2·10⁻⁴` |
| 1000 | `0.2232, 1.454·10⁻⁴` | `−1.378 + 0.289i` | `2.50 + 7.63i` | `≤ 10⁻¹⁶` | `1·10⁻³, 4·10⁻³`; `7·10⁻⁴` |

* **Existence and uniqueness [NE + D].**
  * `M` depends only on `h_∓`; its singular values are the same at all `Δ_g` and both grids. It is non-singular, but ill-conditioned in
    one direction (`1.45·10⁻⁴`): the two null vectors have nearly parallel `(F, F′)(−i/2)`, with `F′/F ≈ −1.59i` for both.
  * The solution is unique, and stable to `10⁻⁴` relative under grid halving.
  * `c ∝ 1/Λ` (ratios `3.00`, `3.33` against `Λ` ratios `3.03`, `3.35`), because `F₀ ∝ ψ₂′(0) ∝ 1/Λ`.
* **`χ` is analytic in the canonical gauge.** On circles `|z + i/2| = 0.2 / 0.1 / 0.05`:
  * min-norm gauge: `max|χ| = 0.096 / 0.34 / 1.25` (`∝ r⁻²`, the double pole);
  * canonical gauge: `4.3·10⁻³ / 3.4·10⁻³ / 2.9·10⁻³` (bounded).

  The same pattern holds at 300 and 1000.
* **Gauge-independent quantities are unchanged.**
  * V3 (graviton pole): `≤ 3·10⁻¹³`.
  * Abel residual: `≤ 2.5·10⁻¹³`.
  * `δk(0)` changes only at `Δ_g = 100`: `−7.2452·10⁻³ → −7.4103·10⁻³`. The bound factor `(1+δk)^{−½}` goes `1.00364 → 1.00372`, which is
    negligible. At 300 the change is `2·10⁻⁸`; at 1000 none.
* **The canonical weight is large at small `ν` [NE].** The profile in units of the local `ψ₂` is `Δ_g`-independent:

  | `ν` | 3 | 5 | 10 | 15 | 18 | 20 | 25 |
  |---|---|---|---|---|---|---|---|
  | `Re X_can/ψ₂` (100 / 300 / 1000) | `−158 / −151 / −150` | `−184 / −169 / −163` | `+66 / +54 / +51` | `+3.9 / +3.1 / +2.9` | `1.24 / 1.16 / 1.14` | `1.06 / 1.03 / 1.03` | `1.02 / 1.01 / 1.00` |
  | `Re X_min/ψ₂` (u3) | `2.6–2.8` | `1.2` | `0.71–0.73` | `1.04–1.07` | `1.02–1.04` | `1.01–1.03` | `1.00–1.02` |

  * `‖X_can − X_min‖₁/‖ψ₂‖₁ = 9.4 / 1.03 / 0.093` at `Δ_g = 100 / 300 / 1000`.
  * The difference is `< 10⁻² max ψ` beyond `ν ≈ 21 / 20 / 19`.
  * At `ν = 15` the canonical weight is still about `3ψ₂`, against about `1.05ψ₂` in the min-norm gauge.

### 3.3 Verdict on the gauge

**Yes: the CHMRS class (`χ` even and analytic in `|Im ν| ≤ 1`, apart from the `ν = 0` term) selects a unique weight for `Ĉ′_EE`.**
* It is the analogue of the scalar case (H7 §24.6: one condition, one parameter), now with two conditions and two parameters. The class
  exists because `ρ(ν)ρ(−ν) = 1`; for `q²` it does not exist (H7 §24.4).
* **It is adopted here as the definition of `Ω^{EE}`'s weight: `X := X_can`.** The weight is then canonical, and the "gauge must be fixed by
  convention" caveat of H7 §24.4 / H0 §R is removed for `Ω^{EE}`.
* Grade: [NE + D]. The conditions are [D]; existence, uniqueness and stability are [NE] (2×2 non-singular, grid-converged). Analyticity
  near `−i/2` is checked [NE].

### 3.4 What adopting it costs (stated, not hidden)

* By Lemma N, nothing gauge-independent moves: `A`, `B_AdS`, the bounds, the leading-Regge per-state positivity (EE_FUNCTIONAL §4,
  [P] + [M]), and the flat limit / NF02.
* The heavy side **beyond leading Regge at `ν ≲ 20`** is now weighted by `X_can`. That weight oscillates with amplitude up to about
  `180ψ₂` (`ν ≈ 5`), `65ψ₂` (`ν ≈ 10`) and `3ψ₂` (`ν ≈ 15`).
* Any non-range (not `Â_h[strip-analytic]`) beyond-Regge part of `Ψ_EE[state]` there is amplified by those factors. **No existing chain
  was computed in this gauge.**
  * EE_FUNCTIONAL §5's `J = 18, 20` chain at `Δ_g = 300` used the local de-shift at nodes with `ν ≥ 15`. That is the large-`ν` symbol, and
    the canonical weight is about 3× it at `ν = 15`.
* **Rough scale, using H7 §24.3's per-unit sensitivities [M].** These were computed for `q²`'s null vectors and are only indicative here:
  * `Δ_g = 1000`: `≤ 3·10⁻⁵ × |c| ≈ 2·10⁻⁴ N_s`, negligible;
  * `Δ_g = 300`: `≤ 1·10⁻³ × 27 ≈ 0.03 N_s` (real part), against EE's computed margins of `0.45 / 0.48` at `J = 18, 20`;
  * `Δ_g = 100`: `0.3–0.5 × 80`, not a usable estimate. **`n` at `Δ_g = 100` in the canonical gauge is [I]**; it was already [I] for the
    `x ≲ 0.1` region.
* **Alternative, if the cost proves prohibitive.** Keep the min-norm gauge as an explicit convention (as H7 §24.4 recommends for `q²`).
  The canonical choice is then available as a cross-check. Either way the definition is now unambiguous; the choice between them is a
  design decision with a precisely located consequence.

---

## 4. Task 4: `Ĉ₃` is unchanged

`Ĉ₃` is CLS-D's `j′ = 1`, `f = −1` functional with its own dual structure `o_H·q/V⁴` (`l9_points` kind `C3H`). The `q` there is a factor of
**`Ĉ₃`'s own definition**; replacing `q²` by `EE` in `Ĉ₂` does not touch it.

| `Ĉ₃` item | where | depends on `Ĉ₂`'s (2,4) structure? |
|---|---|---|
| normalisation `N̂₃` (`ε₃m⁶Ĉ₃ → 2P`) | H0 §R.3 | no: its own heavy action |
| `Ĉ₃|_E = 0`, route A (bare commutator, F, Q8) | H7 §4 | no |
| route B, "no `Π₂`" (H7 §15) | `Ĉ₃/Ĉ₂[q²] ∝ x²` scaling | `Ĉ₂[q²]` serves only as an `m⁻⁴` reference; the conclusion concerns `Ĉ₃` alone |
| convergence (`k − 1 > J₀`), no inversion | H7 §7, H0 §14 | no (§1.4: structure-independent anyway) |
| crystal image `O⊗Q`, stripped `B₃` | H0 §12 | no |
| leading-Regge off-diagonal `≡ 0` | OFFDIAG §1.2 | no |
| `|ĝ₃|²` light term, H17 | H0 §18 | no |
| near-gap rows (`C3H`, `D = bm(q, o)`) | l9/p8/n3 | no (EE_FUNCTIONAL §5 reused them unchanged) |
| **flat-limit relative normalisation of `ψ₂Ĉ₂` vs `ψ₃Ĉ₃` (NF02)** | H0 §R.7 | **yes, through `Ĉ₂`'s crystal image**: unchanged, because `EE ≡ q²` at the crystal with normalisation ratio `1.00000000` (EE_FUNCTIONAL §1.2 [P]) and both `Ĉ₂`'s are heavy-normalised to `P_l` |
| `k(0) = −2∫ψ₂p⁶ + 2∫ψ₃p⁴` | H0 §R.4 | only through `δk(0)` of `Ĉ₂` (§3.2) |

**Confirmed [D]: nothing in `Ĉ₃`'s analysis used `Ĉ₂`'s (2,4) structure.** The two contact points (flat normalisation, `k(0)`) are
unchanged for `EE`, up to `δk(0)`. (Separately: the parity-odd sector chains use `Ĉ₂`'s `+−` entries computed with `q` in `p3_odd_action`.
That is a `Ĉ₂` item, EE_FUNCTIONAL §5 item 3, not a `Ĉ₃` one.)

---

## 5. Regraded table (for whoever edits `RESULT_D.md` / H0 §R; not edited here)

| statement | grade |
|---|---|
| `EE` has the weights of `q²`, is transverse, and lies in CLS-D's general 13-dim `j′ = 0` space (`d = 3`); not in the printed saddle basis | [P] |
| the printed list is a saddle-quotient basis (CLS-D `kernel_identity_saddle`, reproduced to `2·10⁻¹⁶`); overcomplete by 1 in `d = 3` | [Q + DV] / [D] |
| `e₂₄ = (I(x₀₂)z₂)·(I(x₀₄)z₄)` | [DV] (`3·10⁻¹³`) |
| vanishing and convergence of `Ψ_S` are independent of the dual structure `S` (complete gauge fixing; finite-dim pairing) | [D] |
| **`EE` admissible** | **[D]** (was [D-sketch]) |
| (2,4) side of `EE`'s leading-Regge kernel = 16 × scalar; `c(ν±i) = 16d₁ᵀ/d₁^{sc}` (light-ray side) | [P given H7 §18.1] |
| (HN) in the scalar-CHMRS form, `N̂ = c` | [P given H7 §18.1] |
| H7: `−Ĉ′_EE|_E = 16πG/(ν²+9/4)` | [P given Lemma R, Q3, Q4, D2, §18.1 (D; two inputs, referral not needed)] |
| Lemma R for `EE` | [P given R1–R5], R3 [D] (same singular locus) |
| I14(a) / I14(b) | [P given K1, K2] / [P given K1, K1′, K2, K3]; Prop K verbatim [P]; spinning components [D] |
| Lemma N for `EE` | [P given NE regularity] + [DV] (`10⁻¹⁵–10⁻¹⁶`) |
| canonical (CHMRS-class) gauge exists and is unique | [NE + D] |
| canonical weight: `|X_can/ψ₂|` up to about 180 at `ν ≲ 10`; `‖X_can−X_min‖₁/‖ψ₂‖₁ = 9.4 / 1.0 / 0.09` | [NE] |
| `A`, `B_AdS`, bounds, V3, leading-Regge positivity: gauge-independent | [P given Lemma N] + [DV] |
| `n` beyond leading Regge in the canonical gauge (`ν ≲ 20`) | [I]; rough scale `≲ 0.03 N_s` at 300 and `≲ 2·10⁻⁴` at 1000 [M, q² sensitivities]; 100 open |
| `Ĉ₃` independent of `Ĉ₂`'s structure | [D] |

## 6. Hand-over

**Files.**
* `ads4/code/spinfunc/x1_admissibility.py` → `cache/x1_admissibility.txt` (8 s). Contents:
  * (A) weights and transversality;
  * (B) monomial rank (13 in `d = 3`, 14 in `d = 4`), printed-span tests, and the `d = 3` Gram relation;
  * (C) CLS-D's saddle identity;
  * (D) singular orders at `x₀₂² → 0`, and the `e₂₄` geometric identity;
  * (E) exact (2,4)-side ratio `16` and `c(ν±i)` from the light-ray side.
* `ads4/code/spinfunc/x2_gauge_EE.py` → `cache/x2_gauge_EE.txt` (about 2 min). Contents:
  * (S) synthetic validation;
  * (V) reproduction of H7 §24.2 / §24.6 numbers;
  * (L) Lemma N for `EE`;
  * (G) canonical gauge.

  It writes `cache/x2_chi_EE_can_<Δg>.pkl`, with keys `nu`, `X` (canonical), `Y`, `X_minnorm`, `c`, `Lam` on the `C.SETTINGS` grids.
  `Fz(z, nu, dnu, fp, fm, phi, deriv)` is the strip-continuation routine; it is reusable for any Abel solution.

**Next steps (value order).**
1. **Heavy side in the canonical gauge.**
   * Redo EE_FUNCTIONAL §5's `J = 18, 20` chain at `Δ_g = 300` with the node weights multiplied by `X_can/X_local`. In the window `ν ∈ [15, 20]`
     that ratio is about `3 → 1`.
   * Then do the even `J ≤ 16` and odd rows. If the change is `≲ 0.05 N_s`, keep the canonical gauge. Otherwise state min-norm as the
     convention and quote the canonical value as the cross-check (§3.4).
2. **`Δ_g = 100`.** The canonical weight reaches `−180ψ₂` at `ν ≈ 5`, so any `x = 0.05` node there must use contour data and the exact
   non-local `X_can`. A local de-shift is meaningless at that weight.
3. **Exact spinning magic-block identity beyond leading Regge** (H7 §24.8 item 2), now for `EE`. If the `O(r)` term of the split
   representation lies in the range of `Â_h`, the gauge question becomes moot for all states, and the cost in §3.4 disappears.
   * For `EE`'s (2,4) side, `x₀ ≠ 0` gives `e₂₄ = (I(x₀₂)z₂)·(I(x₀₄)z₄)` in closed form (§1.3), which should make this computation simpler
     than for `q²`.
4. For RESULT_D/H0 §R:
   * replace "the covariant `EE` lies outside CLS-D's printed list [D-sketch]" by §1.5;
   * record (HN) for `EE` in the scalar form with `N̂ = c`;
   * record the gauge decision of §3.3–§3.4.
