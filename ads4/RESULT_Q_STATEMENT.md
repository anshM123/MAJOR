# Result Q: two-sided bound on the R⁴ coupling `g₄` (consolidated statement, 22 Nov 2026)

*Consolidation of 22 Nov 2026. This is a writing pass with no new computation; every number is quoted from the source named next to
it. This file replaces the version of commit `17dc133`, whose eight stacked update boxes superseded its own body. That version is in
git history; Appendix A summarises the sequence.*

**Status.** This is a conditional result.
* The implication "per-state positivity of the two fixed functionals ⇒ the bounds" is [D], given the analytic inputs of §4.
* Per-state positivity is computed [NE + M] on the binding and near-binding channels, and assumed elsewhere.
* The decisive below-fold input, `H_RS` for the binding state at `ν < 6`, is **[I]**.

**Sources and precedence** (later wins):
1. the previous `RESULT_Q_STATEMENT.md` (body and boxes);
2. `RESULT_Q_STAGE0–7`, `RESULT_Q_REVISION`, `RESULT_Q_V4PRIME/V4DPRIME`;
3. `RESULT_Q_V4FINAL` (V4F);
4. `REVIEW_RESULTQ_2026-11`;
5. `RESULT_Q_FWD_CLOSE` (FWD);
6. `REVIEW_UNIFIED_2026-11` (RU; authoritative on grades);
7. `HRS_CHECK` (HRS);
8. `RESULT_Q_OG_CHECK` (OG);
9. `RESULT_Q_FERMION_HEFT` (FH).

Note: HRS, OG and FH are later than RU and refine three of its items. Each such refinement is marked.

**Grades.**
* **[P]**: proved. **[D]**: derived. **[D-sketch]**: derivation sketched. **[DV]**: derived and checked numerically.
* **[Q]**: quoted. **[NE]**: numerical evidence. **[NE + M]**: computed in the chain/HYB model. **[M]**: model.
* **[A]**: assumption. **[I]**: open.

---

> **Update (HRS_CLOSE.md, unreviewed): H_RS for Q's binding state [I] → [D-sketch], i.e. [D] given premise RD.** The argument is analytic and covers every ν window, including ν < 1.
> - At every order of the 1/Δ expansion, the engine's integrand is one fixed multivalued function times a polynomial. Integration by parts reduces every term to a rational multiple of one base integral, the same on any cycle (s-/φ-odd terms vanish; there are no α-odd terms).
> - So the crystal cycle = c_Γ(ν) × engine order by order, and c_Γ cancels in the state/reference ratio.
> - The hypotheses of the argument are checked on the engine's pieces to a⁴–a⁵ [DV]: lattice-shifted exponents shared by scalar/EE2/EE4, no interior poles, no integer exponents, no boundary layers.
> - **Remaining premise RD** (the crystal integral is Regge-dominated, with the far region e^{−O(Δ)}): [DV] for the scalar, [NE] for EE at ν ≥ 10, [D-sketch] for EE below 10.
> - **New K = 5 images:**
>   - truncation is ≤ 0.01 c₂ per harmonic (300), with |Δs_d| ≤ 0.015 (300) / 0.0015 (1000);
>   - the binding FF/N is unchanged (+0.1427 / +0.1428 at δ = 0.2), with the margin ≥ 85× the bound.
> - **Covered:** ++ J = 0, 2. OO J = 0 is covered via Γ_OO = Γ_++ (direct application [D-sketch]). **Not covered:** Result D's (336, 28/32) window, where the obstacle is the J²/Δ resummation, not the cycle.

> **Update (RESULT_Q_EINSTEIN_CHECK.md, unreviewed): Lemma A + D_LR stays [D]; the Einstein `⟨TTTT⟩` itself is still not evaluated [I].**
> * **New evidence [DV].** The mechanism was checked on interacting spinning correlators: `T′ = σ ⊗ T_free`, where `G_σ` is an exact sum of AdS₄ Witten diagrams.
>   * The diagrams: contacts with `J₀ = 0, 2, 4`, and a massive spin-2 exchange (`τ = 2`, `J₀ = 2`) in the s, t and u channels.
>   * The pointwise zero holds at `Ψ₄` and `Ψ₅` to 10⁻¹⁴–2·10⁻¹², and in the collinear approach down to Δφ = 0.025.
>   * Controls fail at O(1).
> * **Not exercised:** the graviton's `τ = 1`, its own tensor structures, and the exact Einstein coefficients at external Δ = 3.
> * **Definition of the light side (clarified, [D-sketch]).** The light side is **D_LR on the `J₀ < 4` terms plus term-wise (F)+(P) on the `J₀ ≥ 4` terms**. It is not D_LR applied to the full correlator.
>   * In the test family the swapped-form zero holds whenever `J₀ ≤ k`, because the joint Regge corner then has an integrable exponent `k − 1 − J₀`.
>   * So the R⁴ value `−Ĉ₄ = 2g₄` cannot sit in the fixed-angle integrand. Applied to the full correlator, D_LR would give it 0.
>   * This is how `g̃₄` is already defined (§1 table); the Einstein zero is unaffected.
> * **Weight control qualified.** The light-transform-weight control is one-sided: lower weights also give 0. It therefore tests convergence of the light-ray integral, not light-transform covariance (RESULT_Q_OG_CHECK §3.4).

> **Update (REPRO_ENGINE.md, independent): the engine ratio ρ(ν) and Q's binding `s_d` are reproduced by a different route, with no discrepancy beyond stated precision.**
> * **Route.** A direct 7-D quadrature of the exact pointwise integrand on D at finite Δ. It uses no 1/Δ expansion, and the conserved ⟨TTO⟩ coefficients were re-derived symbolically.
> * **ρ agreement:** `≤ 3·10⁻⁷` (Δ_g = 300) and `≤ 3·10⁻⁸` (Δ_g = 1000) relative to the engine's K = 5 ratio on ν < 4.
> * **`s_d`:** `+0.50967` / `+0.50776`, against `+0.5096` / `+0.5077`.
> * **FF/N:** `+0.14265` / `+0.14282`.
> * **Sensitivity.** Reaching the fail point `s_d = 1.776` needs a worst-sign shape error of `2.1·10⁻⁴` / `1.9·10⁻⁵` in ρ, i.e. ≈ 2.5 % of its beyond-leading-Regge deviation. That is ~600× the measured direct-vs-engine difference.
> * **New grade.** The engine's algebra and truncation for the binding state are [DV] on D.
> * **Not tested:** premise RD (crystal cycle ≡ real section D), because both routes integrate over D. H_RS stays [D-sketch].
> * **Record note.** Production K = 4 is missing some a⁴ pieces (HRS_CLOSE §3.5). Its effect on `s_d` is +0.0004, inside ±0.015.

> **Update (RD_CLOSE.md, unreviewed): premise RD at ν < 10 is reduced to one cycle-class premise P_cyc, which EE shares with the scalar.**
> * **EE = a·Q·scalar exactly [P + DV].** Q is single-valued rational, every pole of Q is already singular for the scalar, and Q is regular at m₁₄ = 0 and at s → ±1.
> * **Far-region suppression on D′** (the real section cut at the zeros of B): `O(Δ^N e^{−cΔ})`, with `c = 2 artanh(r/2)` (≥ 0.336) [D]. The integrated check gives fitted exponents within 2 % of analytic for EE and scalar [DV].
> * **On the actual crystal thimbles** the far weight at r ≥ 0.3 is pointwise `≤ e^{−66}` (Δ_g = 300) and `≤ e^{−277}` (1000). After the thimble volume factor (≲ e^{40} at 300) this is `≲ e^{−26}`, against failure thresholds `e^{−8.5}` / `e^{−10.8}` [NE + D-sketch]. The r ≥ 0.2 figure (e^{−37}) does not survive the volume factor and must not be quoted (REFEREE_REPORT M10).
> * **Net:** RD for EE at ν < 10 is [D given P_cyc], and H_RS for Q's binding state is [D given P_cyc].
> * **What P_cyc says:** the physical cycle lies in the class of D′. For the scalar it is supported through its periods (engine = exact skeleton).
> * **Side finding.** At ν ≲ 2 the scalar thimble has an angular end on which EE does not decay. This plausibly explains HRS_CHECK's signal loss below ν ≈ 2 [D-sketch].

> **Update (C1_COMPUTE.md, unreviewed): c₁ is one uncomputed universal constant, not a free parameter. Exact (non-linearised) bounds hold for every c₁.**
> * **What vanishes, in H7's prescription:**
>   * the k = 3 part [P given Q2, Q8, F];
>   * the u-channel E×R³ T-exchange Regge residue [D]. At O(t₄) it is helicity-flip, and EE kills every flip configuration; checked across all 16 polarisation assignments.
> * **What c₁ is.** The remaining contact-like piece is ν- and Δ_g-independent [D, given J₀(E×W³) = 2 (D-sketch) and polynomiality].
>   * Its value is [I].
>   * Its natural size is ~10⁻³–10⁻¹ [M].
>   * No bound on it exists [I].
> * **Cubic bound, exact:** `√(t₄²+α_T²) ≤ R₀(x + √(1+x²))`, with `x = 2.160·10⁻³|c₁|` (Δ_g = 300) or `1.670·10⁻⁴|c₁|` (1000) [D].
>   * |c₁| = 100 gives 3.88·10⁻⁶ / 2.94·10⁻⁸.
>   * |c₁| = 1000 gives 1.42·10⁻⁵ / 3.41·10⁻⁸.
> * **Quartic upper end, exact (given κ > 0):** `2L_AdS + 6.81·10⁻⁵c₁²` (300) / `+ 5.44·10⁻⁷c₁²` (1000), taking the supremum over unknown t₄, α_T [D].
>   * It no longer depends on Result D, which removes RU issue 14.
>   * **The `0.57|c₁|` / `0.059|c₁|` figures used above are superseded.**

> **Update (LIGHTSIDE_DIST.md, unreviewed): M3 resolved structurally. Result Q's light side is ONE functional, and it is the one the heavy blocks use [D].**
> * **The functional.** `Ω_k` is CLS-D's term-wise swapped double-commutator form on the s-channel OPE, with no separate collinear or ν-boundary terms.
> * **Decomposition, by Cauchy in Mellin space:** `Ω_k = D_LR + 𝒜_k`, where `𝒜_k` is the Regge arc [D].
>   * The arc is 0 for Regge spin `n < k`; that is Lemma A, now read as "no arc below the borderline" [DV].
>   * At `n = k` the arc is fixed by the top Regge coefficient alone (P) [D + DV].
> * **Checks on EIN's Witten family [DV]:**
>   * term-wise values = whole-function values to 1e-11–3e-5;
>   * the borderline contact gives a finite light value, e.g. −8.7642e-3 at P1, equal to its top-Regge term to ≤ 7e-7;
>   * a borderline exchange (the H7 analogue) matches its top-Regge term to 8e-8.
> * **Consequences:**
>   * The piecewise "D_LR on J₀ < 4 + (F)+(P) on J₀ ≥ 4" light side is now the *derived* value of Ω₄.
>   * `g̃₄ = −½·𝒜₄[light] = −(heavy sum)` given H_K [D].
>   * Result D's H7 pole is the arc of the same functional at `n = k = 2` [D].
>   * No number changes.
> * **Corrections to the Einstein-check box above:**
>   * Its claim that J₀ = k values sit at a = b or in the ν integral is **withdrawn**.
>   * Its "J₀ = 4 contact under Ψ₄" is subcritical in that product family (effective n = J₀ − 2), so it is not an R⁴ control.
> * **Remaining:**
>   * (F), the normalisation `𝒜₄[R⁴] = 2g₄` [D-sketch];
>   * the primaries-vs-Mellin-residue rearrangement [D-sketch];
>   * a = b integrability against C_ν [I].

## 0. Headline

**Theorem-form statement (conditional).** Assume A1 and A2 at `Δ_g ∈ {300, 1000}`, the analytic inputs of §4, and per-state positivity of the upper
functional `Ω^Q_up` and the lower functional `Ω^Q_low` on every twist-allowed heavy block. Then the CFT-defined R⁴ coefficient `g̃₄`
(§1; equal to the bulk `g₄` given (F)+(P)) satisfies

    −O(Λ⁻²) ≤ g̃₄ Λ⁶/(8πG) ≤ 2L_AdS − κ|ĝ₃|²Λ⁸ + ε_{c₁},
    2L_AdS = 17.241·log(2Λ/3) + B_eff,   B_eff = −19.220 (Δ_g = 300) / −19.265 (Δ_g = 1000),
    i.e.  2L_AdS = 72.04 (300) / 92.81 (1000).

* **Cubic terms.**
  * **Case `t₄ = α_T = 0`** (CFT-checkable; then `ĝ₃ = 0`): the `κ` term drops, and so does `ε_{c₁}` (∝ `t₄`).
  * **Case general `ĝ₃`:** the upper bound needs **`κ ≥ 0`** [D-sketch + M], with `κ = 0.110`. Then `ĝ₃ ≠ 0` only tightens it.
  * If `κ` could be negative, Result D's `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴` would allow a shift of order `10³`, about 15× `2L_AdS` (RU issue 7).
  * The lower bound for `ĝ₃ ≠ 0` reads `g̃₄Λ⁶/8πG ≥ |ĝ₃|²Λ⁸/(4c_L) − O(Λ⁻²)` with `c_L > 0`. It uses the sign of the W³W³ light term
    of the lower functional (§2.4).
* **`c₁`** (free) is the coefficient of the `E × W³` light term (`k = 2` and `k = 3` parts).
  * `|ε_{c₁}| = 2|c₁t₄|Λ²∫ψ₂`, with `∫ψ₂ = 1.016` for the δ = 0.2 weights.
  * With Result D's `t₄` bound: `|ε_{c₁}| ≤ 0.57|c₁|` (0.79 % of `2L_AdS`, 300) and `≤ 0.059|c₁|` (0.06 %, 1000).
  * Importing D's bound makes Q inherit D's `H_LJ`, `H_M`, `H_R`, `H_edge` and `H_RS`. Carrying `c₁t₄` explicitly avoids that.
* **Heavy side.**
  * Full-set minimum of FF/N for the δ = 0.2 functional: **`+0.1426 ± 0.0017`** (SB gauge `+0.1409`).
  * Binding pair: `++` J = 0 and OO J = 0 at the gap, at both `Δ_g` [NE + M given `H_RS`].
  * The binding margin in the derivative piece is `s_fail/s_d = 1.78/0.51 ≈ 3.5`.
* **Not claimed.**
  * `Δ_g = 100` or any other `Δ_g`.
  * Loop level.
  * The δ = 0.05 functional as a headline (its margin is `+0.0060`, 10 % of δ, and depends on `σ_hU`).
  * **Any cubic closure `|ĝ₃|²Λ⁸ ≤ …` from Q** (§5.3).
* **ABJM:** `g₄Λ⁶/8πG = 2ζ(3)` lies at 3.3 % (300) / 2.6 % (1000) of the bound. It violates A1, so this is a consistency check only (§6).

---

## 1. Setting, normalisations, functionals

* **A1, A2, `c_T`, `G`, `Λ`, `m`** as in `RESULT_D_STATEMENT.md` §1:
  * free real scalar `c_T = 1`;
  * `8πG = 256/c_T`;
  * `Λ = √((Δ_g−2)(Δ_g−1)) = 298.50 / 998.50`;
  * `m² = (Δ−J−2)(Δ+J−1)`;
  * gap state `m_gap² = (Δ_g−2)(Δ_g−1+2J)`.

  These conventions are identical in `resultq/stage6_model.py` and Result D (RU §1).
* **Couplings.** `g₄` is the flat R⁴ Wilson coefficient in Stage 0's normalisation: the type II Virasoro–Shapiro amplitude has
  `g₄/(8πG) = 2ζ(3)/M⁶`, `M² = 4/α′`. `ĝ₃` is the CHLPS `R³` normalisation, with `√(t₄²+α_T²) = 240|ĝ₃|` (AdS units).
* **`g̃₄` (CFT-defined).** After Lemma A removes every tree term with Regge intercept `J₀ < 4`, `g̃₄` is defined by the light side of the
  de-shifted `k = 4` value functional, `−Ĉ₄^{ds}|_light = 2g̃₄`. It equals the bulk `g₄` given (F)+(P) [D-sketch].
* **Upper functional `Ω^Q_up`.**
  * The flat functional is δ = 0.2, `Q3-full-delta0.2` of `code/resultq/functionals_q3.json`: `A = 2ψ₂′(0) = 17.241`, `κ = 0.1099`.
  * It is realised in AdS by the de-shifted `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`, with the C4/C5 weights of Stage 7
    (`functionals_q7_xprec.json`), plus the smeared `k = 2, 3` parts `∫ψ₂Ĉ′_EE`, `∫ψ₃Ĉ₃` (Result D's objects with Q's weights).
  * `s_d` and FF/N are quoted in gauge H0, with H2/H3/H4/SB/XP as variants (spread `≤ 0.016`).
* **Lower functional `Ω^Q_low`.** The value piece: heavy image `P_{−½}(cosh β) > 0` at leading Regge, flat
  `2(|c⁺⁺|²+|c⁺⁻|²)/m⁸ ≥ 0`, light side `2g̃₄` (RU issue 17; `RESULT_Q_REVISION` §3.3).
* **Fixed as of 22 Nov 2026.**

## 2. Light side

### 2.1 Lemma A + D_LR (tree terms with `J₀ < 4` contribute zero)
* **Pointwise Lemma A [D]** (V4F §2.2). For every tree term with `J₀ < 4` (Einstein, `E × W³`, MFT), the integrand integrated over
  `t₁, t₃` and symmetrised under 1↔3 vanishes identically at `φ₁ ≠ φ₃`.
  * Inputs: microcausality, CLS-D endpoint convergence (`J₀ < k`), Lemmas B and C, `J₀(E×W³) ≤ 3`.
  * So this light side is the zero function in every order of the ν-, angle- and `t`-integrations.
  * The earlier ν-inside/outside circularity is resolved.
* **D_LR [D] (consequence [P]).** CLS-D's domain `x₁₃² > 0`, with the operators on their null contours, fixes the cone-preserving
  prescription. On the light-transform manifold, `W_u·f₄ ≡ 0`. The timelike-tube regularisation is not a regularisation of CLS-D's
  integral.
* **Checks at `O(1/C_T)` [DV]** (OG, FH).
  * **Correlators.** `Ψ₄` and `Ψ₅` were evaluated in swapped form on the exact connected `⟨TTTT⟩` of:
    * the free conformally coupled scalar: all 25 (2,4) components, `k = 2–5`, residual `≤ 8·10⁻¹⁵`; per ν (ν = 1, 3) `≤ 5·10⁻¹³`;
    * the free Dirac/Majorana fermion: `≤ 3.8·10⁻¹⁵`; per ν `≤ 1.4·10⁻¹⁵`.

    Hence they vanish on every boson + fermion mixture, including the `t₄ = 0` mixture that has Einstein's `⟨TTT⟩`.
  * **Generic exponents and logs.** GFF `[φφ]₀,₂` and `[ψ̄ψ]₀,₂` families and their Δ-derivatives: `≤ 3·10⁻¹¹` / `1·10⁻¹⁰`.
  * **Controls fail at O(1).**
    * Scalar: non-primary `T` 0.5–1.7; wrong light-transform weight 0.4–1.8; no T-shift 0.03–0.04; wrong symmetrisation 2.00.
    * Fermion: weight, primarity and frame controls 0.04–1.8. The no-T-shift control is degenerate there; the branch rule is pinned by
      the 2- and 3-point validation instead (`κ = i` gives `−1`).
  * **Corner.** The collinear-corner law `Δφ^{−1/2}` (scalar) confirms V4DPRIME's twist-1 homogeneity lemma [DV].
  * **Limitation.** The fermion is blind to the (13) link and the collinear corner in the CLS-D gauge (two Wick classes vanish
    identically [P]). The (13)-channel/corner part of D_LR is therefore checked on the scalar only.
* **Grade.** Lemma A + D_LR [D]; mechanism [DV] at `O(1/C_T)` on two free theories. **The Einstein instance is unchecked [I].** No
  light-ray-usable Einstein `⟨TTTT⟩` exists in the repository. Its leverage on `g₄` is `10¹⁰–10²²` (RU issue 9).

### 2.2 Contacts and heavy-exchange double traces (`J₀ ≥ 4`): (F)+(P) and `H_EFT`
* **(F)+(P) [D-sketch].** The spinning flat limit and polynomial light side: `−Ĉ₄^{ds}|_{R⁴} = 2g₄`, and the `g₅` cancellation. The R⁴
  normalisation control in position space is not reachable: no `J₀ = 4` spinning correlator exists (OG §5.1, FH §6).
* **`H_EFT` splits in two (FH §4.5, refining RU issue 9 and OG §5.2).**
  * **(a)** The term-wise (P) value of a heavy exchange's double-trace part equals its actual (OPE-side) double-trace value.
    * **[DV] in the scalar-Mellin toy** (`J_h = 0, 2, 4, 6`; `Δ_h = 20, 50, 100`). The error at the asymptotic optimum or numerical
      floor is `≤ 1·10⁻⁹` (`Δ_h = 20`) to `≤ 1.5·10⁻⁶` (`Δ_h = 100`) of the heavy image, `≥ 10⁵×` below the bound's sensitivity.
    * Spinning and position-space versions: [D-sketch].
  * **(b)** The heavy-sum / contact-tower interchange ("no conspiracy"): **[A]**, untested. The absolute moments grow super-geometrically,
    so the term-wise series is at best asymptotic.
  * Caution: "double traces = −heavy block" holds only for `J_h ≤ 3`, where it is Lemma A, checked `D = −H` to the numerical floor. For
    `J_h ≥ 4`, `Arc/H ≈ −30 … −35` (`J_h = 4`) and `≈ −2.5·10³ … −3.4·10³` (`J_h = 6`), tending to the flat values `−35`, `−3423`.
    Result Q does not use the `J_h ≤ 3` shortcut.

### 2.3 The `|ĝ₃|²` sector and `c₁`
* **W³W³.** The u-channel diagram has Regge spin 2 and is annihilated by `k = 4`. The s- and t-channel diagrams act like contacts:
  `−Λ²∂_{ν²}Ĉ₄(0)|_{W³W³} = 4πG|ĝ₃|²` [DV algebra in the scalar-Mellin form; spinning D-sketch] (`RESULT_Q_REVISION` §0.3).
* **`κ = 0.110`.** Flat moments of the δ = 0.2 functional, plus H17-type `O(Λ⁻²)` terms. **`κ ≥ 0`: [D-sketch + M]** (RU §4). The earlier
  table grade [I] is superseded.
* **`c₁`:** free (§0). It includes the `Ĉ₃` (`k = 3`) `E × W³` term as well as the `k = 2` one (RU issue 14).

### 2.4 Lower bound
* At `ĝ₃ = 0`: `g̃₄Λ⁶/8πG ≥ −O(Λ⁻²)` **[D-sketch]**, limited by (F)+(P) and `H_K` (RU issue 17). It uses only the value piece. It is
  unaffected by `σ_hU` and by `s_d`. "≥ 0" really means "≥ −O(Λ⁻²) × (EFT data)", from `b₀`, `s_v` and the contact tower.
* For `ĝ₃ ≠ 0`: `≥ |ĝ₃|²Λ⁸/(4c_L)`, with `c_L = 0.509` (measured `s_d`; slope `0.49`) or `0.98` (`s_d` double counted; slope `0.255`).
  * The slope is positive in both readings [D-sketch + M]: W³W³ sign as in §2.3, plus HYB.
  * This item is not in RU §5's wording. It is added here because the lower functional's `|ĝ₃|²` term enters once `ĝ₃ ≠ 0`.

---

## 3. Heavy side: per-state positivity

### 3.1 Which channels can bind (FWD §1, tolerances `s_fail` = the `s_d` at which a channel's FF reaches 0)

| channel | `s_fail` (δ = 0.2) | `s_fail` (δ = 0.05) | measured `s_d` (gap; 300 / 1000) | sign | source |
|---|---|---|---|---|---|
| **`++` J = 0** | **+1.78** | +0.58 | `+0.509 ± 0.016` / `+0.508 ± 0.002` | unfavourable | BF2 |
| **OO J = 0** | +1.78 | +0.58 | `+0.509 ± 0.016` / `+0.510 ± 0.002` | unfavourable (= `++`) | FWD §2 |
| `++` / OO J = 2 | +2.81 | +0.83 | `++` `+0.511 ± 0.010` / `+0.508 ± 0.001`; OO `+0.511 ± 0.010` / `+0.510 ± 0.001` (1.05Λ) | unfavourable | BF2; FWD §2 |
| `+−` J = 4 | +2.38 / +2.10 | +1.04 / +0.75 | `−0.516 ± 0.039` / `−0.510 ± 0.004` | favourable | BF3 |
| `+−` J = 5 (parity-odd) | +5.72 | +1.78 | `−0.539 ± 0.033` / `−0.512 ± 0.004` | favourable | PO §2 |
| `+−` J = 6 | +9.5 / +9.4 | +1.99 | `−0.564 ± 0.012` / `−0.515 ± 0.003` | favourable | FWD §3 |
| `+−` J = 7, 8, 9, 11 | ≥ 20 | ≥ 8.9 | 1000 (used): `−0.48 … −0.52`; 300: `−0.43 … −0.59` | favourable | FWD §3 |
| `++` J = 6, 8 | +36 / +66 | +28 / +52 | `+0.511 ± 0.008` / `+0.54 ± 0.41` | unfavourable | FWD §3 |

* The correction is helicity-odd and universal: `Γ_K^h ≈ σ_h U(ν)`, with `σ_{++} = σ_OO = +1`, `σ_{+−} = −1`, and `|s_d| ≈ 0.51 ± 0.08`.
* Unmeasured channels (`++` J ≥ 10, `+−` J = 10 and ≥ 12, OO J ≥ 4, masses `> 1.475Λ`) have `s_fail ≥ 11`, or are `(Λ/m)⁶`-suppressed.
* **Full-set minimum.**
  * δ = 0.2: `+0.1426 ± 0.0017` (SB `+0.1409`). The 300 state alone is `+0.1427 ± 0.0017`; the 1000 state `+0.1426 ± 0.0002`.
  * δ = 0.05: `+0.0060 ± 0.0013` (SB `+0.0047`).
  * Both `Δ_g`, all gauges [NE + M] (FWD §0; unchanged by HRS).

### 3.2 `H_RS` for the binding state (HRS; refines RU issue 3)

* **Weight.** `s_d = 0.04 Σ_ν χ_d(ν) Re Γ_K(ν)`.
  * Per window, engine `s_d` (H0) is `+0.10` [0,1), `+0.57` [1,2), `−0.05` [2,3), `−0.10` [3,4), `−0.01` [4,8), `≈ 0` above.
  * `ν < 2` carries 68 % of the weight, and `ν < 4` 89 %.
* **Crystal-cycle test** (direct n3/n4 nodes, `J = 0` at the gap at 300, 1000, `Δ = 100`; `J = 2` at 1.05Λ). Crystal = engine beyond
  leading order:
  * **[NE] at `ν ≥ 10`**: `Γ_exact` up to `−209` per harmonic reproduced to 0.3–1 %; e.g. `ν = 20.02`: `−0.3 ± 1.5` against
    `Γ_exact = −73.4`;
  * **[NE, weak] at `6 ≤ ν < 10`**;
  * **[I] at `2 ≤ ν < 6`**: consistent, but the errors (`±39 … ±1800` per harmonic, with ×3 seed inflation) are 3–30× the δ = 0.2
    thresholds;
  * **untested at `ν < 2`**.
* **Propagated.** On `ν ∈ [5, 25]` at 300: `Δs_d = −0.10 ± 0.35` and `ΔFF/N = +0.011 ± 0.039` (×3 errors), against the margin
  `s_fail − s_d = 1.27`. On `[2, 25]`: no information.
* **To break δ = 0.2.** A uniform per-harmonic error of 30 (ν ∈ [0,1)) or 17 ([1,2)) would break it.
* **`H_RS` for Q's binding state stays [I].** `s_d = +0.51` is [NE given `H_RS`].

### 3.3 Other heavy-side items

| item | content | grade | source |
|---|---|---|---|
| `H_OO^Q` | OO J = 0, 2 measured, equal to `++` (helicity lemma at `k = 4`: order 0 to `10⁻⁹`; beyond, `≤ 1.3·10⁻³ c₂` per harmonic; `|Δs_d| ≤ 5·10⁻⁶`) | [NE + M], i.e. **[D given `H_RS`]**: OO and `++` use the same engine (RU issue 11) | FWD §2; PO §1 |
| `H_fwd^Q` | unmeasured forward channels; the `Δ_g = 300` J-drift of `s_d` (`−0.07` at J = 7, 2.4× the truncation error, cause [I]), J ≥ 7 quoted from 1000 (a substitution; tolerance ×17 at δ = 0.05) | [NE + M] given `H_RS`; drift [I], harmless at δ = 0.2 | FWD §3; RU issue 12 |
| `H_smear^Q` | graviton corrections to the smeared `k = 2, 3` parts at the binding states: `≤ 3·10⁻⁵` of B₂ against a tolerance `≈ 24 %` | [NE + M]; edge `x > 0.8` and `k = 3` below the fold bounded by tolerances | REVISION §1 |
| `H_LJ^Q` | large-J near-gap corrections stay below Q's tolerance: `≥ 1.1·10⁴ c₂` (`+−` J = 8, 300), `≥ 1.6·10⁵` (J ≥ 16), against measured `10²–10⁴`. The forward term `a(2J(J+1) − c₄)u³` dominates. | [M], margin ×10–10³ | RU issue 13 |
| HYB | scalar skeleton + leading Regge accuracy away from measured states; `β > 30` | [M] (was [NE]) | RU §4 |

---

## 4. Hypotheses (single list)

| # | hypothesis | role | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|
| 1 | **A1**, **A2** | setting | [A] | — | — |
| 2 | **K1** (KKSDZ 1904.05905, with its loophole), **K1′**, K2/K3 for the de-shifted functionals (validity `J₀ < 7/2`) | Regge | [Q] + [D-sketch] | V4F §5 | closing the loophole; writing out K1′ |
| 3 | **`H_K`**: heavy OPE sum ↔ (corner, ν) interchange for the heavy-normalised `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}` (same statement as Result D's) | the only load-bearing Regge input | **[D-sketch]** ([D given `H_ip`], raw corner `J₀ < ½`) | V4F §2.4 | Proposition K for the heavy-normalised kernel |
| 4 | **Lemma A + D_LR** | `J₀ < 4` light side `= 0` | **[D]**; mechanism [DV] at `O(1/C_T)` (free scalar + free fermion); **Einstein instance [I], unchecked** | V4F §2.2; OG; FH | one swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` (needs O(G) `[TT]` data or a light-ray-usable representation) |
| 5 | **(F)+(P)** spinning | EFT light values (`2g₄`, `g₅` cancellation) | [D-sketch] | V4F §5 | a `J₀ = 4` spinning position-space correlator (AdS contact diagram) |
| 6a | **`H_EFT`(a)**: term-wise (P) of heavy-exchange double traces = actual value | `J₀ ≥ 4` light side | [DV] scalar-Mellin toy; [D-sketch] spinning | FH §4 | spinning toy |
| 6b | **`H_EFT`(b)**: no conspiracy in the heavy sum / contact tower | same | **[A]** | — | — |
| 7 | **De-shift** existence: continuum C4/C5 solutions for all β | heavy images beyond `t = 30` | [NE] to `t = 30`; C5 spinning coefficients [I]. The stored sc3 weights have growing high moments (`m₄ = −1.4·10⁻²`, `m₅ = −27`, ratio `≈ 1.9·10³`, a ν-grid cutoff effect), harmless for `Λ ≥ 300` (ratio `≤ 0.02` per order). | Stage 7; FH §4.4 | — |
| 8 | Adm, H7 (inherited from D) | `k = 2` smeared part | [D] | Result D §4 | — |
| 9 | **`κ ≥ 0`** | only if `ĝ₃ ≠ 0` (upper) | [D-sketch + M] | REVISION §3.3 | spinning W³W³ |
| 9′ | `c_L > 0` (lower slope) | only if `ĝ₃ ≠ 0` (lower) | [D-sketch + M] | REVISION §3.3 | same |
| 10 | **`c₁`** | `E × W³` light term | free parameter, carried | RU issue 14 | Regge residue of `E × R³` |
| 11 | **`H_RS`** | every `s_d`, incl. the binding state | **[I]** for `ν < 6`; [NE] `ν ≥ 10`; [NE, weak] 6–10 | HRS | an EE `k = 4` crystal value at `ν < 4` to ~10–20 per harmonic (lower Δ, variance control for the `R` factor, or an analytic continuation argument) |
| 12 | `H_OO^Q` | OO J = 0, 2 | [D given `H_RS`] | FWD §2 | — |
| 13 | `H_fwd^Q` | unmeasured channels, 300 J-drift | [NE + M]; drift [I] | FWD §3 | — |
| 14 | `H_smear^Q` | smeared parts | [NE + M] | REVISION §1 | — |
| 15 | `H_LJ^Q` | large-J near-gap | [M], margin ×10–10³ | RU issue 13 | — |
| 16 | HYB accuracy; `β > 30` | rest of the heavy action | [M] | RU §4 | — |
| — | `σ_hU` double counting | **not needed** for the δ = 0.2 headline. It decides δ = 0.05 and the (withdrawn) closure. | [I]; the flat-limit argument favours double counting | RU issue 8 | — |
| — | `H_νout`, H_LT, V₄/V₄′/V₄″/V₄‴ | **deleted/replaced** by item 4 | — | V4F | — |

---

## 5. Caveats

1. **The weakest links.**
   * `H_RS` at `ν < 6` for the binding state [I]. An engine error of ×3.5 in the derivative piece at J = 0 would break δ = 0.2, and
     nothing independent excludes it.
   * Then `H_EFT`(b) [A], the unchecked Einstein instance of Lemma A, `H_K` and (F)+(P) [D-sketch].
2. **`g̃₄` vs `g₄`.** Only `g̃₄` is CFT-defined. Its identification with the bulk coupling is (F)+(P).
3. **Cubic closure withdrawn** (RU issue 8).
   * The values "`|ĝ₃|²Λ⁸ ≤ 120–198 / 154–257`" are not claimed.
   * The `Δ_g`-independent lower slope 0.49 (with the measured `σ_hU` correction) is 1.9× CHLPS's flat optimum `0.26`, which an AdS
     bound whose sum rules reduce to the flat ones cannot exceed in the flat limit. This favours the double-counting reading.
   * At most a remark: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257` [D-sketch + M]. No ratio `|ĝ₃|²/g₄` is claimed.
   * Result D's `1.07·10⁴ / 1.43·10⁴` stands under fewer hypotheses.
4. **δ = 0.05.** Its minimum `+0.0060 ± 0.0013` is positive but depends on `σ_hU`. It is not the headline.
5. **Rounding.** `B_eff(1000) = −19.2647` (REVIEW_RESULTQ §1) has been quoted as "−19.27" in earlier documents; to two decimals it is
   −19.26. The effect on `2L_AdS = 92.81` is `< 0.01`.

## 6. ABJM consistency check

* **Dictionary** (Stage 0; arithmetic re-checked by RU §1):
  * `M² = 4/α′`, `L²/α′ = π√(2λ)`, so `Λ ↔ ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.22 λ^{1/4}`. The identification `Λ ↔ ML` holds only at leading
    order in `1/(ML)`.
  * The type IIA `R⁴` term `ζ(3)α′³/32` gives `g₄M⁶/8πG = 2ζ(3) = 2.404`, independent of `λ`.
* **Placement.** `2ζ(3)/2L_AdS = 3.3 %` (300) and `2.6 %` (1000), recomputed: `0.0334`, `0.0259`.
* **`ĝ₃ = 0`** because 10d type IIA supersymmetry forbids an `R³` term. It does *not* follow from "maximal SUSY": ABJM at `k > 2` has
  `N = 6`.
* **`c_T` does not enter the placement.** Do not confuse flavour-current coefficients with `c_T`:
  * for `N = 2` theories, `∂²_m F` gives a flavour-current two-point coefficient, not `c_T`;
  * the identity `c_T = (64/π²)∂²_m F` (free real scalar `c_T = 1`) relies on the mass coupling to a current in the stress-tensor
    multiplet. It is used only in `ABJM_THEOREM`'s separate `c_T` result, consistently with `8πG = 256/c_T`.
* **Outside A1** (light KK modes); only the flat image (F) is used; log-dominated. **Consistency check, not evidence.**

## 7. Literature

* **CHMRS 2106.10274** (Caron-Huot–Mazáč–Rastelli–Simmons-Duffin) is the **methodological precedent**: two-sided `log(Δ_gap)` bounds
  on contact couplings in AdS₄ with dynamical gravity, from uplifted flat sum rules, for scalar correlators. Q's light-side formula and
  `H_EFT` are of CHMRS's (5.3) type.
* **CLS-D 2311.04271** (Chang–Landau–Simmons-Duffin; JHEP 04 (2025) 016) build the spinning light-ray sum rules and the crystal
  dictionary, and leave the uplift of flat bounds to AdS as future work. **Result Q is a first explicit-constant, conditional execution
  of that programme for `g₄`**. Method priority: CLS-D and CHMRS.
* **CHLPS 2201.06602:** the flat D = 4 two-sided graviton bounds (log IR). Stage 0 reproduces the flat ¼ from the same forward rules.
* **To be discussed (not only cited):** Kundu 2104.11238 and the `N = 6` bootstrap (Binder–Chester–Jerdee–Pufu 2011.05728).
* **Novelty [I] at search level.** No `⟨TTTT⟩` sum-rule bound on `R⁴` was found (abstract-level search plus one extended search, RU §1).

## 8. Open problems (value order)

1. `H_RS` at `ν < 4` for the binding state.
2. One swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` (routes: O(G) `[TT]` data in every helicity structure, days; momentum space, impractical).
   Also the (13)-channel/corner part of D_LR on a correlator with a spin-1/spin-2 single-trace exchange.
3. `H_EFT`(b); spinning `H_EFT`(a).
4. `H_K`.
5. (F)+(P) in position space; the R⁴ normalisation control.
6. Spinning W³W³ (would firm up `κ ≥ 0` and `c_L > 0`); `σ_hU`.
7. The 300 J-drift of `s_d`.
8. `c₁`.

---

## Appendix A. History (superseded content, in order)

1. **Stages 0–7** built the forward functionals.
   * δ = 0.05 (`A = 11.45`) was replaced by δ = 0.2 (`A = 17.241`) as the quoted bound.
   * The `Δ_g = 100` column was dropped.
2. **REVIEW_RESULTQ: major revision.**
   * `κ`, the ¼ slope and `4πG|ĝ₃|²` were graded [I].
   * V₄ (swapped form, `[TT]₀` tail) was graded [D-sketch].
3. **REVISION.**
   * `H_smear^Q` → [NE + M].
   * The W³W³ residue in `∂Ĉ₄(0)` was found to be zero (Mellin).
   * Cubic terms were restored at [D-sketch + M], and a closure `|ĝ₃|²Λ⁸ ≤ 120–198 / 154–257` was proposed. **Later withdrawn (RU).**
4. **The V₄ sequence.**
   * V₄′ was false at O(1) for the raw `Ψ₄,ν`; the tube contact term grows like `ε^{−1/4}`.
   * V₄″ was false [DV], because of `e^{πν/4}` growth.
   * V₄‴ worked per ν with H_LT and H_νout, but left a circularity.
   * **V4F replaced all of them** by the pointwise Lemma A [D] and D_LR [D], deleted H_νout, and isolated `H_K`.
5. **PARITY_ODD_CLOSE / FWD.**
   * `+−` J = 5, OO J = 0, 2 and `+−` J = 6–11 were measured.
   * `H_OO^Q` [I] → [NE + M]; `H_fwd^Q` coverage [I] → [NE + M].
   * The minimum was unchanged.
6. **RU.**
   * Sound with corrections.
   * New or renamed: `H_EFT`, `H_RS` [I], `H_LJ^Q`; `H_K` [D-sketch].
   * Cubic closure withdrawn; `κ ≥ 0` required for `ĝ₃ ≠ 0`; `c₁t₄` carried; HYB [M].
7. **HRS.** `H_RS` regraded by window; the binding state is still [I]. The "512 c₂" item is the resolved `k = 2` comparison artefact
   (BF3 §2), not an open `k = 4` discrepancy.
8. **OG, FH.**
   * Lemma A mechanism [DV] at `O(1/C_T)` on free scalar and fermion; the Einstein instance is still [I].
   * `H_EFT` split into (a) [DV toy] and (b) [A].
   * The OG reading "`H_EFT` is an assumption only for `J_h ≥ 4`" is refined by FH: (b) is assumed for all `J_h`.

**Numbers unchanged throughout:** `2L_AdS = 72.04 / 92.81`; full-set min `+0.1426 / +0.1427` (δ = 0.2).

**Superseded statements of the previous body:**
* hypotheses 3 ("V₄ incl. the `[TT]₀` tail") and 9 ("`H_OO^Q` assumed");
* the regrading table's "`κ` [I]" and "Lemma A [D-sketch]";
* "Next steps" 1–3 (all done: smeared `k = 2` correction, swapped-form `Ψ₄` on MFT, W³W³ residue).

**Unresolved between documents.** The absolute de-shift moments of the sc3 weights are `2.2, 83, 1.4·10⁴, 3.6·10⁶` in V4F but
`1.49, 33, 3.5·10³, 6.6·10⁵` in FH ("different normalisation or key"). No bound depends on them.
