# Result D: bound on the cubic couplings `t₄`, `α_T` from the EE light-ray functional (consolidated statement, 22 Nov 2026)

*Consolidation of 22 Nov 2026. This is a writing pass with no new computation; every number below is quoted from the source named next
to it. This file replaces the version of commit `c473d9e` (20 Nov, with its post-review box). That version is kept in git history,
and Appendix A summarises what changed.*

**Status.** This is a conditional result and **not a theorem at any `Δ_g`**. The implication "per-state positivity ⇒ bound" is [D]
and its constants are [DV]. Per-state positivity is computed on a finite set of states [NE + M] and is assumed elsewhere, through the
named hypotheses of §4.

**Sources and precedence.** Where sources disagree, the later item in this list wins:
1. the previous `RESULT_D_STATEMENT.md` and its post-review box;
2. `RESULT_D.md` (all "Update"/"Correction" bullets);
3. `REVIEW_EE_2026-11.md`;
4. `HYPOTHESES_CLOSE.md` (HC), `HYPOTHESES_CLOSE2.md` (HC2);
5. `PARITY_ODD_CLOSE.md` (PO);
6. `HLJ_ANALYTIC.md`, `HWIN_THIMBLE.md`, `FOLD_WINDOW_LARGEJ.md`, `SMOOTH_BOUND.md`;
7. `BELOW_FOLD*.md` (BF, BF2, BF3, BF_SPINJ), with BF_SPINJ as regraded by the unified review;
8. `HRS_CHECK.md`;
9. `RESULT_D_CLEANUP.md` (CL);
10. `REVIEW_UNIFIED_2026-11.md` (RU). Its issues 1–6 and its §4 table are authoritative on grades.

**Grades.**
* **[certified]**: computer-assisted proof.
* **[P]**: proved. **[P given X]**: proved from the named inputs.
* **[D]**: derived, not independently checked. **[D-sketch]**: derivation sketched.
* **[DV]**: derived and checked numerically.
* **[Q]**: quoted from the literature.
* **[NE]**: numerical evidence. **[NE + M]**: direct crystal values inside the near-gap chain model.
* **[M]**: model, i.e. an approximation whose accuracy is itself a hypothesis.
* **[A]**: standing assumption of the setting.
* **[I]**: open.

A chain of inputs carries the grade of its weakest link.

---

> **Update (HLJ_CLOSE.md, unreviewed): H_LJ(300) is not closed. Its Regge window for J ≥ 34 stays [I], and H_LJ(1000) stays [D-sketch + M].**
> * **Diagnosis.** The large a¹ term is a J-independent EE normalisation, of size 16–21/Δ relative.
>   * The genuinely J-dependent part is the helicity structure. It is small in the window (25–60 c₂ per order), but its 1/Δ series at fixed J/Δ is asymptotic, not convergent.
>   * No resummation was achieved [I].
> * **Analytic normalisation `ρ_an`.** It has no fitted constant, agrees with the old fits within their errors, and Prop. RS bookkeeping holds at J = 28, 32 [DV].
> * **Window states at Δ_g = 300, `ν < 8.3`: [M] → [NE + M].**
>   * (336, 28): `λ/N_s = +0.759`; stresses give 0.751–0.772.
>   * (336, 32): `+0.837`.
>   * The binding margin is the structure error, at 3–6× the engine's structure terms. Every measured discrepancy has the safe sign.
>   * (336, 28) at +0.759 is now the lowest computed large-J value at 300, below the previously quoted J = 32 +0.783.
> * **High J (J = 64, 100).** The engine's last term is 20–55 × T_B, so it cannot decide H_win either way.
> * **Unexplained:** a negative structure offset of 20–330 c₂ at ν = 5–10 [I].
> * **Net for Result D:** the only remaining [I] is the H_LJ(300) Regge window for J ≥ 34. The headline at Δ_g = 1000 stands (PRL draft switch = open).

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

> **Update (HK_CLOSE.md, unreviewed): H_K no longer needs J₀ < ½. With an analytic cutoff on the de-shift weights it follows under the known bound J₀ ≤ 1, with margin ½.**
> * **No-go [P].** Any de-shift weight supported on a bounded ν-range, including every stored grid weight, regrows at the raw Regge rate by Paley–Wiener. The interchange then needs J₀ < ½, the old H_ip.
> * **The fix: Proposition K_c [D given K1′, K2ν, LRC, K3, DS].**
>   * Each de-shift weight is multiplied by an entire erf cutoff c(ν): ν_c = 150, a = 6 for Q; ν_c = Λ+100, a = 10 for D.
>   * The de-shift identity forces strip-analyticity [D], so the cutoff does not reintroduce growth.
>   * Tonelli then gives an absolutely convergent heavy sum for J₀ < 3/2 (k = 2) and J₀ < 5/2 (k = 4).
>   * Checks: exactly solvable toys [DV]; real Q weights [NE].
> * **Every computed number is unchanged.** ε_c ≲ e^{−180} (Q) and ≲ 10⁻⁴¹ (D).
> * **Ĉ₅^{ds} and Ĉ₃ need nothing new** (plain Fubini).
> * **Net grade.** H_K stays [D-sketch] by weakest link, but it is no longer an independent hypothesis. It reduces to K1′ (existing), K2ν (new: per-ν kernel growth ≤ e^{γν} away from the corner, [D-sketch]), LRC and DS.
> * **Required change.** The cutoff must be added to the functional definitions (§1 above; RESULT_D_STATEMENT §2; RESULT_Q_STATEMENT).
> * **Side note.** The C5 de-shift is poor beyond t ≈ 15, which affects per-state positivity at β ≳ 15 (the HYB region).

> **Update (LEMMAS_CLOSE.md, unreviewed).**
> * **Lemma C′: → [D]**, with a closed-form check [DV] at ≤ 2·10⁻⁹.
>   * On Ω = {Re t₁, Re t₃ > 0} ∖ (0, ½] the swapped integrand is one single-valued analytic function, by an exact identity keeping z, z̄ off [1, ∞) [P].
>   * Lemma A is [D given C′]. The Einstein instance is [D] (not numerically evaluated), which answers M7.
> * **K1′: → [D]**, uniform on the Regge box at J₀ = 1, with no edge dependence.
> * **K2ν: → [D].** `|𝒦_ν| ≤ K̄` for all real ν.
> * **LRC leading form:** [D + DV] for scalar k = 2, 4 and EE k = 2; [D-sketch] for EE at k = 4.
> * **Correction to the H_K box above: the LRC uniform remainder is [I].**
>   * Over the large-η part of the Regge box the corrections are O(ρ = re^β), not O(r). Prop K_c's "extra factor r" step therefore fails uniformly, and with K1′ alone the remainder again needs J₀ < ½.
>   * The same gap sits in H7 §14's corollary for Ĉ′_EE.
>   * It closes given either of two inputs, both [I]:
>     * LRC\*: the de-shift cancellation survives the O(ρ) corrections;
>     * LC\*: uniform light-cone suppression of dDisc with exponent κ > ¼.
>   * **H_K stays [D-sketch].**
> * **a = b integrability:**
>   * the kernel bound is [D] (integrable iff growth p < 7/2);
>   * J₀ < 4 terms carry nothing at a = b [D];
>   * the borderline R⁴ arc in that corner is [I], tied to LRC\*/LC\*.

## 0. Headline

1. **Implication [D]; constants [DV].** Let `𝒯` satisfy A1 and A2 and the analytic inputs of §4. Suppose the fixed functional `Ω^EE`
   (§2) acts positively on every twist-allowed heavy conformal block, in the sense of a positive semidefinite per-state action matrix.
   Then
   **`√(t₄² + α_T²) ≤ 3.13·10⁻⁶` at `Δ_g = 300` and `≤ 2.89·10⁻⁸` at `Δ_g = 1000`** (tree level, AdS units).
   * The coefficient `c₁` of a possible `E × W³` light term is **not computed and is carried as a free parameter**. It moves the bound by
     `∓0.22 %` (300) or `∓0.017 %` (1000) per unit `c₁`.
   * H17 adds `±0.24 %` / `±0.021 %`.
   * There is no argument for `|c₁| ≤ 1`, and none is assumed.
2. **Positivity is not established.**
   * **`Δ_g = 1000`.** The EE block is computed [NE + M] on:
     * even `J ≤ 24`, `28`, `32`, `36`, `40`, `48`;
     * odd `J ≤ 31`;
     * masses from the twist edge to `5Λ`.

     The minimum is `+0.2009` (`J = 5`, `+−`, hybrid state at the twist edge `1.00499Λ`). The large-spin near-gap region is the
     hypothesis `H_LJ(1000)`, graded [D-sketch + M].
   * **`Δ_g = 300`.** The EE block is computed [NE + M] on even `J ≤ 24`, `28`, `32` and odd `J ≤ 15`, with direct small-`ν` nodes on
     every computed `J ≥ 4` (`J = 0, 2` are still on the superseded small-`ν` treatment, at margin `≈ 1`). The minimum is `+0.2038` (`J = 5`, `+−`, twist edge `1.0166Λ`). `H_LJ(300)` is [D-sketch/M] in the bulk
     and **[I] in the Regge window for `J ≥ 34`**. The `ν < 8.3` window of `J = 28, 32` is [M].
   * **The parity-odd even-`J` (OO) block** reduces to the EE `++` diagonal through an identity that is [NE] (`H_OO`, §4).
   * **`Δ_g = 100`** is not established (§5.4). Nothing is claimed at any other `Δ_g`.
3. **The inputs are mostly computational.** The only physics-level inputs are A1 (twist gap for all single traces), A2 (tree level) and
   the KKSDZ Regge bound with its stated loophole. `H_M`, `H_R`, `H_edge`, `H_OO`, `H_RS` and `H_LJ` are statements about explicit
   conformal-block integrals. `H_K` (swap of the heavy sum with the corner/`ν` integration) is an analytic statement, graded [D-sketch].
4. **Literature (§6.4).** Rewritten in CHLPS form, the bound is `|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21` [D].
   * It is about 10× weaker than the transplanted CHLPS optimum. That optimum is asserted by CHLPS, not derived in AdS.
   * Its parametric form differs from AKT's, and the constants are not comparable.
   * No published `⟨TTTT⟩` sum-rule bound on R³ couplings was found. This rests on literature searches of 19 and 22 Nov, so novelty is
     [I] at search level.
   * The method is CLS-D's (2311.04271) programme of uplifting flat bounds to AdS, which CLS-D left as future work. CHMRS (2106.10274)
     is the methodological precedent, for scalars.

---

## 1. Setting and conventions

* **A1 (spectrum) [A].**
  * `𝒯` is a unitary CFT₃, parity-even or parity-violating, with a unique stress tensor `T`.
  * Every exchanged single-trace operator other than `T` has twist `Δ − J ≥ Δ_g`, for all `J` (including `J = 0, 1, 2`). There are no
    other light single traces.
  * This is the convention of every scan (`omegau/common.state`). The older "`Δ − J ≥ Δ_g − 2` for `J > 2`" (`RESULT_D.md` §1) is wrong.
  * A1 excludes all known holographic CFT₃s, ABJM included, because they have light Kaluza–Klein towers.
* **A2 (tree level) [A].** `O(1/c_T)` truncation with `c_T ≫ Δ_g²`. Loops are dropped.
* **Normalisations.**
  * `c_T` is normalised so that a free real scalar has `c_T = 1`.
  * `8πG := 256/c_T` in AdS units. This follows from `c_T = 32L²/(πG)`, reproduced by RU §1.
* **Heavy states `(Δ, J)`.**
  * `α = Δ − J − 2 ≥ Δ_g − 2` and `m² = α(α + 2J + 1)`.
  * `Λ := √((Δ_g − 2)(Δ_g − 1))`, which is `98.499 / 298.50 / 998.50` at `Δ_g = 100 / 300 / 1000`.
  * The gap state of spin `J` has `m_gap² = (Δ_g − 2)(Δ_g − 1 + 2J)`.
  * Further variables: `b̂ = (2J+1)Λ/m`, `cosh β = 1 + (2J+1)²/(2m²)`, `x = ν/m`, `p = ν/Λ`.
  * "Near the gap" means `m ≤ 1.3Λ`.
* **Couplings.**
  * `⟨TTT⟩` is written in the basis `(c_T, t₄, α_T)`.
  * In the bulk, `t₄² + α_T² = 57600|ĝ₃|²/L⁸`, so `√(t₄² + α_T²) = 240|ĝ₃|` in AdS units (`DICTIONARY_T4.md`).
  * The `ĝ₃` normalisation is CHLPS's: the project reproduces their (3.7), `1598 log − 3211`, as `1596.7, −3208.0` [DV].
* **Units.** `c₂` = `m² ×` relative correction (the graviton correction `Γ` of a block). `N_s` is the per-state normalisation of the
  action. "`λ/N_s`" is the smallest eigenvalue of the per-state action matrix in units of `N_s`.

## 2. The functional `Ω^EE` (fixed; not to be changed)

`Ω^EE := ∫₀¹ dp ψ₂(p) Λ⁴ Ĉ′_{EE,Λp} + ∫₀¹ dp ψ₃(p) ½Λ⁶ Ĉ_{3,Λp}`
(sources: `EE_FUNCTIONAL` §1–§3, `EE_FOUNDATIONS` §1–§3, `H0` §R, `S1_S4_FIXES` for `Ĉ₃`).

| ingredient | definition | grade |
|---|---|---|
| building block | CLS-D's light-ray-commutator functional `Ψ_{k,ν,j′}[f G]` on `⟨TTTT⟩` (arXiv:2311.04271, eq. `spinningfunctional_definition`) | [Q] |
| `k = 2` part | `j′ = 0`, `f = v′ − u′` (unsubtracted), dual (2,4) structure `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²` | definition; admissible [D] |
| de-shift | `Ĉ′_EE := c(ν)⁻¹𝒜_g⁻¹[Ψ_EE]`, `c(ν) = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))` (EE's leading-Regge kernel `h₋ = g·c(ν+i)`, `h₊ = g̃·c(ν−i)`, CHMRS's scalar `g`) | [P given the H7 §18.1 kernel (D)] |
| smeared form | `Re∫dν X(ν)Ψ̂_{EE,ν}`, where `X` solves the generalised Abel equation `C[aY] − i tanh t·S[bY] = C[ψ₂]/cosh t`, `X = aY`. `X` vanishes at the real zeros `ν = 1.171513, 6.642858` of `h₋ + h₊`. | [D] + [NE] (residual `≤ 3·10⁻¹⁵`) |
| **gauge** | the Abel solution has a 2-dimensional null space (index `w = 5`, supported on `ν ≲ 18`). **Min-norm gauge** (correction orthogonal to the null space). | definition |
| `k = 3` part | `Ĉ₃ := Re[Ψ_{3,ν}[o_H·q]/N̂₃(ν)]`, `j′ = 1`, `f = −1` (unsubtracted), no de-shift (its spinning leading-Regge image is exactly `σ_h c₃(ν)P_l`) | [P given kernel reading (D)] + [DV `4·10⁻¹³`] |
| reality | `Re` of the complex sum rule; phase convention `logK → logK − iπ(Δ+J)` on every record | definition; legitimacy [D] |
| weights | NF02: `ψ_k(p) = (1−p)⁵ p X_k(p)/k₀`, degree-8 polynomials (`code/robust/functionals_noforward.json`), `k₀ = 27.48509817`, `k(0) := −2∫ψ₂p⁶ + 2∫ψ₃p⁴ = 1` | definition |

**Properties.**
* **Flat limit [P given H0 §17's crystal + polarisation map] + [NE].**
  * On CLS-D's physical polarisation plane at the `D = 4` crystal, `EE ≡ q²` holds exactly. So the bulk-point limit of `Ω^EE` is the
    stripped NF02 functional, whose flat positivity certificate applies: minimum `𝔉/N = 0.200` [certified].
  * One corner is float-checked only: `+−`, `J ≥ 400`, `1.5Λ ≤ m < 1.633Λ`.
  * In the engine's full null-circle contraction, `EE/q² → 1` only as `m → ∞` at fixed `J` [NE] (e.g. `J = 18`, `m ≈ 1600`:
    `0.89 − 0.19i`; REVIEW_EE m7). The previous grade "Flat [P]" was too high.
* **Leading-Regge heavy action.** It is helicity-diagonal, with `P_{−½+iν}(cosh β)` on both helicities. The leading-Regge helicity
  off-diagonal vanishes identically [P]; this was the obstruction of the `q²` functional.
* **What the chains compute.**
  * Even `J` at `Δ_g = 300`: the non-local min-norm `X` below `ν = 25`, on direct nodes (z4 protocol).
  * Elsewhere, and for odd `J`: the local first-order de-shift `(1 ± u_EE)/2` [M].
  * The local and non-local versions differ by `0.02–0.10 N_s` at 300, in the favourable direction (REVIEW_EE S4). The previous
    "`≤ 3·10⁻³ N_s`" was wrong. The difference is negligible at 1000.
* **Chronology of the choice (REVIEW_EE S3).**
  * Both the structure (`q² → EE`, motivated by the leading-Regge off-diagonal `≡ 0` [P]) and the gauge (canonical → min-norm) were
    fixed after heavy-state failures had been seen.
  * Any fixed solution of the Abel equation is a valid functional, so the bound is unaffected.
  * It does make `J = 18, 20, 24, 28` at `Δ_g = 300` **in-sample** evidence; they are labelled "in-sample" in §5.
  * **`Ω^EE` (min-norm, NF02) has been fixed since 20 Nov.** Everything computed after it was fixed is out-of-sample.
  * With the canonical gauge the same data give `n = 1.55` (spline), or `−22 … −40` (direct, CV) at `J = 24` (300). So the statement is
    about the min-norm functional specifically.

## 3. Light side and the bound

**Light side** (tree level; `H7_I14_LIGHTSIDE`, `EE_FOUNDATIONS` §2, `H0` §R.8).
* Both functionals are unsubtracted, so every `[TT]` double trace sits at a double zero:
  * parity-even: `2 sin²(π(τ−10)/2)`;
  * parity-odd: `2 cos²(πτ/2)` at odd `τ`.

  Contact terms drop out (L0, [D]). The tree-level light side is the identity plus the single `T` block.
* `−Ĉ′_{EE,ν}|_E = 16πG/(ν² + 9/4)` holds exactly (H7, [D]), and `Ĉ₃|_E = 0` [P given Q2, Q8, F].
* The `|ĝ₃|²` term enters with coefficient `k(0)(1 + δ)`. Here `δk(0)_EE = −7.245·10⁻³ / −7.666·10⁻⁴ / −6.830·10⁻⁵`
  (`Δ_g = 100 / 300 / 1000`), plus the H17 remainder.
* **`c₁` (free).** A possible ν-independent `E × W³` light term `16πG c₁t₄`. `c₁` collects both the `k = 2` (`Ĉ′_EE`) and the `k = 3`
  (`Ĉ₃`) parts (RU issue 14). It is not computed [I].
  * Remark (REVIEW_EE m5, unproved): `c₁ = 0` is plausible at leading Regge, because EE's kernel is helicity-diagonal while the linear
    `ĝ₃` term lives in the helicity-flip amplitude.

**Sum rule.** With heavy `≥ 0`: `16L_A − k(0)(1+δ)|ĝ₃|²Λ⁸ ≥ 0`, where `L_A = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`. `G` cancels. Hence

**`√(t₄² + α_T²) ≤ 240·√(16L_A(Λ)/(k(0)(1+δ)))/Λ⁴ · (1 + O(c₁Λ⁻²))`**,

with `16L_A = A log Λ + B_AdS + O(1/Λ)`, `A = 16ψ₂′(0) = 3022.34`, `B_AdS = −6552.66`.

| `Δ_g` | `Λ` | `16L_A` (exact) | `A log Λ + B_AdS` | **bound** (exact `p`-integral, `k = 1`) | `δk(0)` factor | H17 band | per unit `c₁` | cubic closure `|ĝ₃|²Λ⁸ ≤` | status of positivity |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 98.499 | 7508.51 | 7320.0 | 2.21·10⁻⁴ | ×1.0036 | ±2 % | ∓2.37 % (corrected from ≈ ∓1.7 %, DG_SCAN) | 7.51·10³ | **not established** (§5.4) |
| 300 | 298.50 | 10734.61 | 10670.9 | **3.13·10⁻⁶** (3.1321·10⁻⁶) | ×1.0004 | ±0.24 % | ∓0.22 % | 1.07·10⁴ | conditional (§4, §5.3) |
| 1000 | 998.50 | 14339.59 | 14320.4 | **2.89·10⁻⁸** (2.8913·10⁻⁸) | ×1.00003 | ±0.02 % | ∓0.017 % | 1.43·10⁴ | conditional (§4, §5.2) |

Cross-checks:
* `A`, `B_AdS` and the bounds were re-derived independently by REVIEW_EE §0. They were recomputed again for this note from `16L_A` and
  `Λ` (`3.13206·10⁻⁶`, `2.89126·10⁻⁸`).
* They are gauge-independent (Lemma N). The log form underestimates the exact integral by 0.6 % at 300; the table uses the exact integral.

---

## 4. Hypothesis table (one list; grades per RU §4 and later documents)

| # | hypothesis | grade | evidence | what would upgrade it |
|---|---|---|---|---|
| **A1** | twist gap `Δ − J ≥ Δ_g` for every single trace `≠ T`, all `J`; unique `T` | [A] | setting | — |
| **A2** | tree level, `c_T ≫ Δ_g²` | [A] | setting | — |
| **K1** | KKSDZ Regge bound (1904.05905) `J₀ ≤ 1` for spinning correlators; needed `J₀ < 2` | [Q], **with KKSDZ's stated loophole** (first/second-sheet boundary, their §4.2.6) | `H7_I14` §2, §7.1 | closing the loophole |
| K1′ | uniformity of K1 towards the light-cone edges | [D-sketch] | review V | write it out |
| K2, K3 | kernel power counting at every boundary; per-family Rindler positivity ⇒ swap lemma | [Q]/[D] | `H7_I14` §7, §14 | — |
| **H_K** | interchange of the heavy OPE sum with the (corner, `ν`) integration for the heavy-normalised de-shifted kernel `Ĉ′_EE` (the old "V(b)/I14(b)"; `η → ∞` decay via Prop. K). Same statement as Result Q's `H_K`. | **[D-sketch]** ([D given `H_ip`]: raw per-`ν` heavy-sum convergence, `J₀ < ½`, not implied by KKSDZ) | RU issue 6; `RESULT_Q_V4FINAL` §2.4 | a Proposition K for the heavy-normalised kernel, or `H_ip` |
| V(a) | `Ω^EE` annihilates the non-perturbative `⟨TTTT⟩` | [P given K1, K2] | `H7_I14` §7 | — |
| **Adm** | `EE` is an admissible `j′ = 0` dual structure | [D] | `EE_FOUNDATIONS` §1 | — |
| **H7** | `−Ĉ′_EE|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` | [D] | `H7_I14` §3–§5, §18–§23 | one independent small-`ν` check of `h₋ + h₊` (a defect would move `B_AdS`, not `A`) |
| **L0** | `[TT]` double zeros; contacts drop out | [D] | `H0` §15; `PARITY_ODD_ROUTEA` §12 | — |
| Flat | flat limit = stripped NF02; certificate applies | [P given H0 §17 map] + [NE]; NF02 [certified] (one float corner) | `EE_FUNCTIONAL` §1, §7; `NOFORWARD_FUNCTIONAL` | — |
| Re | `Re` prescription; strip `logK − iπ(Δ+J)` | definition; legitimacy [D]. The strip **is load-bearing for odd `J ≥ 17` at 1000**: state `Im logK ≡ 0`, reference `≡ π`. The code strips and asserts this (`q3_analysis._strip`), and the numbers are right (REVIEW_EE m3). | `PARITY_ODD_ROUTEA` §12 | — |
| Gauge / Lemma N | min-norm gauge; null directions invisible to `A`, `B_AdS`, flat limit, graviton pole, leading-Regge action | definition; Lemma N [P given NE-observed regularity] + [DV `10⁻¹⁵`] | `EE_FOUNDATIONS` §2.5, §3 | — |
| **H17** | `|ĝ₃|²` light coefficient `k(0)(1+δ)` | **[D-sketch + M], carried over from `q²`**. The `Π̂₄` admixture of `Ĉ′_EE` was not derived. Separate `t₄²`, `α_T²` coefficients at `O(Λ⁻²)` are not excluded; if present, the bound is on a weighted norm. Effect on the bound: `±0.24 % / ±0.021 %`. | `H0` §18; REVIEW_EE m4 | spinning `R³` exchange to `O(r²)` |
| **c₁** | the `E × W³` light term (`k = 2` and `k = 3`) | **free parameter, carried explicitly** (`∓0.22 % / ∓0.017 %` per unit) [I] | REVIEW_EE m5; RU issue 14 | Regge residue of the `E × R³` `T` exchange (Lemma R with `n = 2`) |
| **H_RS** | the analytic-radial engine's EE image beyond leading order equals the crystal-cycle integral | **[I]** where D uses it. Scalar [DV]; EE at LO [DV]; EE `k = 2, 4` at `J = 0, 2`, `ν ≥ 10` [NE] and `6 ≤ ν < 10` [NE, weak] (HRS_CHECK); `ν < 6` [I]. At `J = 28` the entry-resolved test fails at 3–16σ (RU issue 2). | HRS_CHECK; RU issues 2–3 | an independent crystal check below the fold at large `J` |
| **H_M** | the near-gap chain (model M) is accurate on the computed states | [M], validated [NE] on every row with direct data | §5; HC §5; CL §0.3 | interval versions of the chains |
| — `K₂` beyond criterion | the `K₂` Laplace rows at large `J²/m` | `K₂` **fails** (×2–5 in amplitude and in phase from `x ≈ 0.2`) at `(336, 28)` (`J²/m = 2.35`), `(336, 32)` (3.08) and `(1048, 48)` (2.2). All three states are now on complete direct sub-rows. Validated up to `J²/m = 1.73` (`(336, 24)`, state values `≤ 1 %`) and 1.54 (`(1040, 40)`, `≤ 0.05 N_s`). The band 1.6–2.2 is not systematically tested [I]; e.g. the secondary `(500, 28)` row (1.58) is still `K₂`. | HC §5; CL §0.3, §5 item 4 | direct sub-rows in the 1.6–2.2 band |
| **H_edge** | the light-cone-edge window `x ∈ (0.8, Λ/m]` | **[NE + M]** where computed (direct `n1` nodes `x = 0.8–0.9`). Computed on every row with a near-gap state at 300 for `J = 4–24`, and at 1000 for `J = 4–14`, `16`, `17`, `31`, `36`, `40`. Largest unfavourable shift of a near-gap state `−0.016 N_s`; of a per-`J` minimum `−0.006`. **[M]** on the rest: 300 `J = 28, 32` and `J ≥ 26`; 1000 even `J = 18–32`, odd `J = 19–29`, `J = 48`. | HC §2; CL §0.2 | edge rows at larger `Δ` for `(336, 28/32)`; the 1000 rows listed |
| **H_OO** | positivity of the parity-odd even-`J` block (CHL's `X^O_same`) | **[D given `Γ_OO = Γ_++`]**. The identity holds exactly at order 0 [P given the helicity lemma (D, DV)] and to `≤ 0.05 / 0.12 c₂` beyond, on all `K₂` rows `J = 4–22` and on direct nodes at `ν ≥ 8.3` [NE]. Beyond order 0 its support uses the engine (`H_RS`). Untested below the fold. Directly computed OO minima: `J = 4–10` at both `Δ_g`, `J = 12–22` at 300. | HC2 §1; PO §1; RU issue 11 | a derivation of `X = 0` beyond order 0 |
| **H_R** | Regge/large-`b̂` region: beyond-leading-Regge corrections at `m ≳ 1.3Λ` stay below the leading-Regge + curvature margin | [M]. Leading-Regge EE band model positive (min `0.2056`). Computed corrections `≤ 0.005 N_s` (1000, `J ≤ 32`), `≤ 0.02` (`J ≤ 40`), `−0.04 … −0.49` with every state passing (300, `J = 12–28`). Scalar remainder `|S|/F_hyb ≤ 122/m²` [NE]; but the 300 records behind it used `N₀` configurations with an unexplained pathology at `ν/Λ ≳ 0.85` (REVIEW_EE m9). | `EE_NEARGAP` §6; `ROUTEA_LARGE_B` | spinning `O(r)` Regge kernel |
| **H_LJ(Δ_g)** | near-gap, large-spin states outside the computed set (§5.2/§5.3): `H_LJ ⇐ H_AG ∧ H_win ∧ H_edge-amp` [P given these] | **1000: [D-sketch + M]. 300: [D-sketch/M] in the bulk, [I] in the Regge window for `J ≥ 34`. `(336, 28/32)` window `ν < 8.3`: [M].** | HLJ_ANALYTIC §5; HWIN_THIMBLE; RU issue 4 | the spinning `O(r)` Regge kernel at `β ≤ 0.53` (makes `H_win` [D]); graviton window data at 300, `J ≥ 34` |
| H_cont | continuity in `m` between grid masses | [NE]. Twist-edge masses are now sampled for every computed `J` at 300 (CL dense grid `J = 4–22`; HC: `J = 24` at `1.078Λ`, `J = 28` at `1.091Λ`, `J = 32` at `1.102Λ`) and at 1000 for `J = 4–16` (CL) and `J ≥ 17` (gap masses in the chains). Exception: `J = 15` at 1000 is not in CL's table. | CL §0.4; REVIEW_EE m2 | — |
| H3, H9a/b, H10, H14, H15 | odd-`k` moments; forward/smeared contacts; EFT power counting; object identity | **not needed** (only unsubtracted `k = 2, 3`; contacts drop out) | REVIEW_ROUTEA §3b | — |

---

## 5. Per-state positivity: coverage, one value per state

### 5.1 Conventions

* Per state, the test is that the EE 2×2 (`++`, `+−`) action matrix is PSD. The OO 1×1 block is listed separately where it was
  computed. Values are min `λ/N_s` over the masses of the chain, from the twist edge to `5Λ`, at the quoted protocol.
* **Protocol labels (REVIEW_EE p1):**
  * **direct** = direct small-`ν` nodes (z4 protocol: n3 nodes `ν ≥ 8.3`, non-local min-norm `X` below `ν = 25` for even `J`, value
    held below the lowest node), with `K₂` rows inside their validated range elsewhere;
  * **direct + sub-row** = complete direct `n1` sub-row replacing `K₂`;
  * **+ edge** = direct edge nodes `x = 0.8–0.9`;
  * **conv. A** = the superseded small-`ν` treatment (`K₂` + spline below `x = 0.05`, local weights);
  * **model** = leading-Regge + curvature (κ-uniform / hybrid / band).
* **In-sample** marks states that were seen before the functional was fixed (REVIEW_EE §0, S3). Everything else is out-of-sample.
* "Passes" means "passes in the adopted protocol". Every grade is [NE + M] unless stated.
* **The minimum over all computed states is always the parity-odd `+−`, `J = 5` state at the twist edge.**

### 5.2 `Δ_g = 1000`

| `J` | value (min `λ/N_s` [`m/Λ`]) | protocol | source | grade |
|---|---|---|---|---|
| **5 (binding)** | **`+0.2009`** (hybrid, twist edge `1.00499Λ`); `+0.2021` at `1.01Λ` on the exact-skeleton grid | w2 chain + edge | `EE_DG1000_LARGEJ` §4; REVIEW_EE §2; HC §2 (edge, unchanged) | [NE + M] |
| 4, 6, 7, 8, 9, 10 | `+0.2041` [1.004] / `+0.2296` / `+0.2684` / `+0.4309` / `+0.4542` / `+0.5886` | twist-edge grid + edge | CL §0.2, §0.4 | [NE + M] |
| 11, 12, 13, 14, 16 | `+0.3482` / `+0.6759` / `+0.4076` / `+0.6982` / `+0.7191` | twist-edge grid + edge | CL §0.2 | [NE + M] |
| 15 | `+0.869` (grid from `1.01Λ`; no twist-edge mass, no edge nodes) | w2 chain | previous statement §9.2 | [NE + M]; `H_edge` [M] |
| even 0, 2, 18–24, 28 | positive, above the binding state; diagonal corrections `≤ 0.09`, off-diagonal `≤ 0.001 N_s` | EE_NEARGAP chain (grid from `1.01Λ`) | `EE_NEARGAP` §2; `SMOOTH_BOUND` §4 | [NE + M]; `H_edge` [M] for `J = 18–28` |
| 17 | `+0.840` | rows + edge | HC §2 | [NE + M] |
| odd 19, 21, 23, 25, 27, 29 | `+0.916 / +0.942 / +0.887 / +0.918 / +0.954 / +0.979` | rows `Δ = 1000 + J`, odd-row despike | previous statement §9.1 | [NE + M]; `H_edge` [M] |
| 31 | `+0.988` | rows + edge | HC §2 | [NE + M] |
| 32 / 36 / 40 | `+0.904` / `+0.957` / `+0.9795` | own-mass rows (+ edge at 36, 40); `J = 40` validated by a 25-node direct sub-row (`≤ 0.05 N_s`) | `EE_DG1000_LARGEJ` §0, §4; HC §2 | [NE + M] |
| **48** | **`+0.915 ± 0.08`** (jackknife) at `1.15Λ`; variants `+0.67 … +0.99` | direct + sub-row (27 nodes, `x = 0.125–0.8`) | CL §0.3 | [NE + M]; `H_edge` [M] |
| even 26, 30, 34, 38, 42–46 | no rows; neighbours positive with margin `≥ 0.9` | model (interpolation in `J`) | — | inside `H_LJ(1000)` |
| even `J ≥ 50`, odd `J ≥ 33`, near the gap, up to the twist cap (`J ≈ 344` at `m ≤ 1.3Λ`) | κ-uniform model `≥ +0.90` (`J = 300`: `+0.902`) | model | `SUBLEADING_REGGE` §0.3 | **`H_LJ(1000)`** |
| `m ≥ 1.3Λ`, `J ≤ 40` | corrections `≤ 0.02 N_s` (`≤ 0.005` for `J ≤ 32`) | rows interpolated/held in `1/m` | `EE_NEARGAP` §6 | [NE + M] |
| `m ≥ 1.3Λ`, larger `J`; large `b̂`; wedge | band model min `0.2056`; `F_hyb > 0` [certified] on `b̂ ∈ [274.1, 1997]` and `[91.2, 2Λ]` (scalar hybrid) | model | `EE_FUNCTIONAL` §4; `ROUTEA_LARGE_B` | `H_R` [M] |
| OO block, `J = 4 / 6 / 8 / 10` | `+0.215 / +0.309 / +0.503 / +0.605` | `K₂` OO rows | HC §3 | [NE + M] |
| OO block, `J ≥ 12` | not computed; equals the EE `++` diagonal via `Γ_OO = Γ_++` | — | HC2 §1; PO §1 | `H_OO` |

**`H_LJ(1000)`, precise form (HLJ_ANALYTIC §5.1, refined by HWIN_THIMBLE).** The region is: every twist-allowed state with
`m ≤ 1.3Λ` and even `J ≥ 50`, even `J ∈ {26, 30, 34, 38, 42, 44, 46}`, or odd `J ≥ 33`. `H_LJ ⇐ H_AG ∧ H_win ∧ H_edge-amp`:
* **`H_AG`**: on `ν ∈ [ν_B, Λ]`, `ν_B = c·m/(2J+1)`, the beyond-model correction is a sum of two travelling waves in class AG
  (decay `≥ 0.97β` in `|Im ν|`).
  * [D] to all orders of the saddle expansion: non-degenerate saddles for `J = 50–344`, and the saddle phase hits all 23–56 zeros of the
    exact scalar action.
  * [NE + M] beyond all orders.
* **`H_win`, refined form**: `|Γ − κ − c₀D_h| ≤ T_B` on `ν < ν_B`, with the `D_h` (leading-Regge-shaped) class sent to its own much
  larger tolerance.
  * That step is HLJ §5.2 [D-sketch]; its class tolerance is computed at 1000 only.
  * **The sign-free form `|Γ − κ| ≤ T_B` is violated at `(1100, 100)`**: C2 `+−` = `−2154 ± 43` against `T_B = 1409` (fw: `−2123 ± 172`).
    The excess is `D_h`-shaped (`c₀ ≈ 3–6.5·10³`), and the remainder after removing it is inside `T_B`.
  * At `(1064, 64)`: `|Γ − κ| ≤ 209 ± 216` on `ν = 13–25`, against `T_B = 1913` [NE].
  * `ν < 10` is not measurable at any `J`. For `J ≥ 200` the whole window lies below `ν ≈ 10`.
* **`H_edge-amp`**: oscillating edge amplitude `M ≲ 3·10⁶ / 9·10⁶` at `J = 48 / 64`, against the exact-scalar envelope. Margin `×3` at
  `J = 48`, `≥ 19` at `J ≥ 64`. Graviton `x > 0.8` is unmeasured at large `J` [M].

| `J` (1000) | 48 | 64 | 100 | 150 | 200 | 300 | 344 |
|---|---|---|---|---|---|---|---|
| `T_B(c = 3)` / `(c = 5)` (`c₂`, C2, worse helicity) | 2456 / 1314 | 1913 / 1020 | 1402 / 734 | 1184 / 623 | 1114 / 633 | 1057 / 789 | 1095 / 881 |

**Grade: `H_LJ(1000)` [D-sketch + M]** (RU issue 4). `RESULT_D.md`'s "[D + NE + M]" is superseded.

### 5.3 `Δ_g = 300`

| `J` | EE block: value (min `λ/N_s` [`m/Λ`]) | OO block | protocol | source | grade / label |
|---|---|---|---|---|---|
| 0 / 2 | `+0.999` [1.0] / `+1.000` [1.7] | — | conv. A (small-`ν` irrelevant at this margin) | HC §4.4 | [NE + M] (superseded small-`ν` treatment) |
| 4 | `+0.2137` [1.013, gap] | `+0.230` | direct spot nodes + edge (shift 0) | HC §4.4; CL §0.2 | [NE + M] |
| **5 (binding)** | **`+0.2038`** [twist edge `1.0166Λ`, hybrid] | — | direct (rows 200 + 400) + edge (all four rows) | CL §0.1, §0.2 | [NE + M] |
| 6 | `+0.3075` [1.25] | `+0.3078` | direct | PO; CL §0.1 | [NE + M] |
| 7 | `+0.7033` [1.095] | — | direct (rows 200 + 400) + edge | CL | [NE + M] |
| 8 | `+0.4268` [1.055] | `+0.4787` [1.65] | direct + edge | CL | [NE + M] |
| 9 | `+0.4506` [1.03] | — | direct (rows 200 + 400) + edge | CL | [NE + M] |
| 10 | `+0.6033` [2.0] (near gap `+0.6698`) | `+0.6038` [2.0] | direct + edge | CL | [NE + M] |
| 11 | `+0.6687` [1.036] | — | direct + edge | CL | [NE + M] |
| 12 | `+0.6706` [1.055] | `+0.886` | direct + edge | CL §0.2; HC2 §1 | [NE + M] |
| 13 | `+0.7746` [1.10] | — | direct + edge | CL | [NE + M] |
| 14 | `+0.7487` [1.046] | `+0.998` | direct + edge | CL §0.2; HC2 §1 | [NE + M] |
| 15 | `+0.7853` [1.07] | — | direct (rows 336/400) + edge | CL | [NE + M] |
| 16 | `+0.884 ± 0.022` (edge shift `−0.004`, measured on the conv. A chain) | `+1.007` | direct | HC2 §0; HC §2.3 | [NE + M] |
| 18 | `+0.872` [4.0] (near gap `+0.922` with edge) | `+0.994` | direct + edge | HC §4.4; CL §0.2; HC2 §1 | [NE + M]; **in-sample** |
| 20 | near gap `+0.905` [1.065] with edge; all `m` `≥ +0.863` | `+0.936` | direct + edge | HC §4.4; CL §0.2; HC2 §1 | [NE + M]; **in-sample** |
| 22 | `+0.844` [1.071] | `+0.842` | direct + edge | CL §0.2; HC2 §1 | [NE + M] |
| 24 | `+0.820` [1.5] (gap `1.078Λ`: `+1.011`) | — | direct + edge | HC §2.3; `SMOOTH_BOUND` §4 | [NE + M]; **in-sample** |
| **28** | **`+0.799 ± 0.011`** [1.5] (gap `1.091Λ`: `+0.927 ± 0.235`) | — | direct + sub-row (16 nodes) + pooled validated thimble node at `ν = 5.0` (s2 + oo3, scalar `≤ 1.7σ`), held below `ν = 5` | HC2 §3.2–3.3 | bulk and `ν ≥ 5`: [NE + M]; `ν < 5`: [M]; **in-sample** |
| **32** | **`+0.783`** [1.15] (jk `±0.016`; twist edge `1.102Λ`: `+0.826`) | — | direct + sub-row (16 nodes), held below `ν = 8.3` | HC §5 | bulk [NE + M]; window `ν < 8.3` [M] (no validated direct datum) |
| even 26, 30 | no rows; neighbours positive | — | model | — | inside `H_LJ(300)` |
| even `J ≥ 34` (largest near-gap `J ≈ 103`), odd `J ≥ 17`, near the gap | κ-uniform model `≥ +0.92` for `J ≤ 150` | — | model | `SUBLEADING_REGGE` §0.3 | **`H_LJ(300)`** |
| `m ≥ 1.3Λ`, `J = 12–28` | corrections `−0.04 … −0.49 N_s`, every state passes | — | rows | `EE_NEARGAP` §6 | [NE + M] |
| `m ≥ 1.3Λ`, larger `J`; large `b̂`; wedge | band min `0.2056`; `F_hyb > 0` [certified] on `b̂ ∈ [162.6, 597]` and `[91.2, 2Λ]` (scalar hybrid); wedge `J > 20` `≥ 0.27` | — | model | `EE_FUNCTIONAL` §4; `ROUTEA_LARGE_B` | `H_R` [M] |

Notes.
* Odd `J` use the local de-shift, because z4's non-local term is implemented for even `J` only [M] (CL §5 item 2).
* Below the lowest direct node (`ν ≈ 8.3`) the value is held [M]. This is irrelevant at `J ≤ 22`, where failure needs a uniform
  `Γ_c ≲ −10⁴ c₂` there. It is decisive at `J ≥ 28`.
* **The `(336, 28/32)` window `ν < 8.3`.** Failure thresholds depend on the baseline (RU issue 18); all are uniform `Γ_c`.

  | `J` | baseline | failure threshold |
  |---|---|---|
  | 28 | no thimbles | `≈ −5400 c₂` on `ν < 8.3` |
  | 28 | validated pooled thimble node | `≈ −7.3·10³ c₂` on `ν < 5` |
  | 32 | no thimbles | `≈ −4250 c₂` on `ν < 8.3` |
  | 32 | pooled thimble seeds that fail the scalar `Im` check at `11σ` | `≈ −1.4·10³ c₂` on `ν < 5` |
  | 28, 32 | BF_SPINJ engine values | an extra `≈ −4.3·10³ c₂` |

  * Validated data at `ν = 8.3` and `ν = 5.0` are `O(10³)`.
  * The BF_SPINJ engine gives `Γ_c = O(10²–10³)` with both signs below the fold, and `+0.787` (`J = 28`, order a²) / `+0.834`
    (`J = 32`, order a¹). This is **[M] evidence only, given `H_RS`**:
    * its one-constant calibration does not validate the engine;
    * out of sample, the entry-resolved residuals miss by 3–16σ (a¹) and 4–90σ (a²) (RU issue 2).
  * RU's reading is "plausibly positive, margin ×2–4". **`J = 32` is the most fragile state at 300.**

**`H_LJ(300)`, precise form.** The region is: every twist-allowed state with `m ≤ 1.3Λ` and even `J ≥ 34`, even `J = 26, 30`, or odd
`J ≥ 17`, plus the `ν < 8.3` window at `J = 28, 32`. Same decomposition as at 1000.
* `H_AG` is [D].
* `H_win` needs `|Γ − κ| ≲ 180–390 c₂` (C2, `c = 3`) on `ν < 6–20`. **No graviton window data exist for `J ≥ 34`.**
* The only window data (`J = 24, 28`, `Γ(ν → 0) ≈ 430–480`) are of the size of `T_B`. The sign-free bound is therefore insufficient.
  `H_win` holds with margin only if the window correction keeps the leading-Regge (`D_h`) shape at subleading order (HLJ §5.2,
  [D-sketch]).
* At `(364, 64)` and `(400, 100)` the X-range exceeds `T_B` in the `D_h` class. The remainder (`≈ 150–250`) is at the level of
  `T_B ≈ 250`, not clearly inside.
* `H_edge-amp` has margin only `×1.6` at `J = 34`.
* Odd `J ≥ 17` have no rows, and the odd correction was still growing at `J = 15` in the old treatment.

| `J` (300) | 34 | 40 | 50 | 64 | 80 | 100 | 103 |
|---|---|---|---|---|---|---|---|
| `T_B(c = 3)` / `(c = 5)` (`c₂`) | 387 / 214 | 343 / 190 | 298 / 166 | 247 / 150 | 205 / 144 | 179 / 145 | 178 / 145 |

**Grade: `H_LJ(300)` [D-sketch/M] in the bulk, [I] in the window for `J ≥ 34`; `(336, 28/32)` window [M]** (RU issue 4, §4).
`H_LJ(300)` is materially weaker-supported than `H_LJ(1000)`.

### 5.4 `Δ_g = 100` (status only)

* Main protocol: `n = 1` on even `J ≤ 10` and odd `J = 5, 7, 9`; binding odd `J = 5` at 1.05Λ, `0.2286` [NE + M].
* It is not robust:
  * the odd extra-row state (`J = 9`, 1.1Λ) gives `+0.03`, with MC `0.005 ± 0.044` [I];
  * near-gap `J ≥ 16` has no rows, and both transfer models fail at `J ≥ 18–24` [M].
* **Not established.** The nominal `2.21·10⁻⁴` is not claimed (`EE_DG100.md`).

---

## 6. Statements

### 6.1 The implication [D]

*Let `𝒯` satisfy A1 and A2. Assume K1 (with KKSDZ's loophole), K1′, K2, K3, `H_K`, Adm, H7, L0 and H17, and let `c₁` be the
(uncomputed) `E × W³` coefficient. If every twist-allowed heavy block has a positive semidefinite per-state action matrix under `Ω^EE`,
then*

`√(t₄² + α_T²) ≤ 240·√(16L_A(Λ)/(k(0)(1+δ)))/Λ⁴ · (1 + O(c₁Λ⁻²))`,

*i.e. `≤ 3.13·10⁻⁶` (`Δ_g = 300`) and `≤ 2.89·10⁻⁸` (`Δ_g = 1000`), `∓0.22 % / ∓0.017 %` per unit `c₁`, `±0.24 % / ±0.021 %` from H17.*

Proof: heavy `≥ 0` ⇒ `−Ω_light ≥ 0` ⇒ `k(0)(1+δ)|ĝ₃|²Λ⁸ ≤ 16L_A`. Constants [DV] (REVIEW_EE §0).

### 6.2 Strongest supported statement

> **Result D (conditional).** Assume A1 and A2 [A], and the analytic inputs: K1 [Q, loophole], K1′ [D-sketch], `H_K` [D-sketch],
> Adm, H7, L0 [D], H17 [D-sketch + M], with `c₁` free.
> * **At `Δ_g = 1000`**, `√(t₄²+α_T²) ≤ 2.89·10⁻⁸` holds given `H_M`, `H_edge`, `H_OO`, `H_R` and `H_LJ(1000)` [D-sketch + M].
>   Per-state positivity of the EE block is computed [NE + M] on even `J ≤ 24`, `28`, `32`, `36`, `40`, `48` and odd `J ≤ 31`, from the
>   twist edge to `5Λ`, with minimum `+0.2009`. It holds in the model [M] elsewhere.
> * **At `Δ_g = 300`**, `√(t₄²+α_T²) ≤ 3.13·10⁻⁶` holds given the same plus `H_LJ(300)` ([I] in the window for `J ≥ 34`) and the
>   `(336, 28/32)` window [M].
>   * Computed [NE + M]: every even `J ≤ 24` and `28`, `32`; odd `J ≤ 15`; minimum `+0.2038`.
>   * `J = 18, 20, 24, 28` are in-sample.
> * **At `Δ_g = 100`** it is not established. Nothing is claimed at any other `Δ_g`.

**What would make it stronger, in order of value:**
1. Graviton window data (`ν < 8.3`) at 300 for `J ≥ 34` and at `(336, 32)`. A validated below-fold method is needed; it does not exist.
2. The spinning `O(r)` Regge kernel at near-gap `β ≤ 0.53`. This makes `H_win` [D] at both `Δ_g` and removes most of `H_R`.
3. `H_K`, i.e. Proposition K for the heavy-normalised kernel, or `H_ip`.
4. `c₁`.
5. `H_edge` on the remaining rows. `K₂` in the `J²/m = 1.6–2.2` band.
6. Interval versions of the chains (`H_M`).
7. The KKSDZ loophole and K1′.

### 6.3 Conditional fallback

Let `U(Δ_g)` be the `H_LJ` region. Assume positivity on its complement, and **`H_U(ε)`**: the states in `U` contribute to `Ω^EE` at
least `−ε·16L_A`. Then `√(t₄²+α_T²) ≤ 240√((1+ε)16L_A)/Λ⁴` [D given the hypotheses]:
* `ε = 0.1`: `3.29·10⁻⁶ / 3.03·10⁻⁸`;
* `ε = 1`: `4.43·10⁻⁶ / 4.09·10⁻⁸`.

`H_U` is implied by `H_LJ` and is a property of explicit block integrals. It is not a spectral condition, unlike the `q²` route's
`H_band` (`RESULT_D_CONDITIONAL`).

### 6.4 Literature

* **CHLPS form [D].** The AdS graviton pole `16πG/(ν² + 9/4)` makes `m_IR = 3/(2R_AdS)` exact at log accuracy. The result is
  `|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21` (`+ O(1/Λ)`).
* **Comparison with the CHLPS optimum.**
  * The log coefficient 3022.34 compares with CHLPS's SDPB optimum 24.9 (2201.06602 (4.4)). The ratio is 121.4, i.e. ×11.0 in the bound.
  * Decomposition [NE]:
    * ×1.05 from CHLPS's extra rules;
    * ×33.6 from the absence of forward terms (forced in AdS by uncontrolled contact content);
    * ×2.58 from the margin δ = 0.2;
    * ×1.34 from the large-`b̂` proxy.
  * Transplanting CHLPS's optimum with the same `m_IR` gives `3.09·10⁻⁷ / 2.80·10⁻⁹` (`|ĝ₃|²Λ⁸ ≤ 104.5 / 134.5`). This is asserted by
    CHLPS, not derived in AdS: it would need forward and helicity-flip sum rules in AdS.
* **Parametric form.** `√(t₄²+α_T²) ≲ 240√A·√(ln Δ_g)/Δ_g⁴`. AKT (1805.07393) give `t₄ ≲ ln Δ_g/Δ_g⁴` (`d = 3`, unknown constant). The
  parametric forms differ, and the **constants are not comparable**, because AKT's normalisation was not converted. The earlier claim
  "stronger by `√ln Δ_g` than AKT" is withdrawn.
* **Positioning (RU issue 16).**
  * CLS-D (2311.04271; JHEP 04 (2025) 016) build the spinning light-ray sum rules and the crystal dictionary to flat-space graviton sum
    rules, and list "the remaining steps for uplifting flat space bounds to AdS" as future work. Result D is **a first explicit-constant,
    conditional execution of that uplift** for the cubic couplings. The priority on method is CLS-D's.
  * CHMRS (2106.10274) is the methodological precedent: two-sided `log Δ_gap` AdS₄ bounds on scalar contact couplings with dynamical
    gravity, from uplifted flat sum rules.
  * The curvature-induced helicity off-diagonal, and its removal by `EE`, is not in the literature as far as the 19 Nov search found.
* **Novelty [I] at search level.** No published paper derives `R³` (or `t₄`, `α_T`) bounds from `⟨TTTT⟩` sum rules. This rests on an
  abstract-level search (19 Nov) plus one extended search (22 Nov).

---

## 7. Open items (value order)

1. `ν < 8.3` at 300 for `J ≥ 28`, with `J = 32` the most fragile state. No validated below-fold method reaches `ν ≈ 2–5` at the needed
   `|Γ_c| ≲ 2000 c₂`. Thimble sums fail the scalar `Im` check at `J = 32`.
2. `H_LJ(300)` window for `J ≥ 34`, and the Regge-window structure at subleading order (spinning `O(r)` kernel).
3. `H_RS` below the fold at large `J`: entry- or helicity-resolved engine tests, or an analytic normalisation (RU issue 2).
4. `H_K`; `c₁`; H17 for `EE` (`Π̂₄` admixture).
5. `H_edge` at `(336, 28/32)`, at 1000 even `J = 18–32`, odd `J = 19–29`, `J = 48`, and `J = 15` (1000).
6. `K₂` adequacy at `J²/m ≈ 1.6–2.2`.
7. A derivation of `Γ_OO ≡ Γ_++` beyond order 0. Its `O(Δ⁻⁴)` residual mechanism is [I].
8. The 300 `N₀` pathology behind `H_R`'s scalar remainder (REVIEW_EE m9). The numerical `x₀₂² = 0` edge-scaling check for `EE`
   (REVIEW_EE m10).
9. `Δ_g = 100`.

---

## Appendix A. History and superseded content

**A.1 What this rewrite changed (RU issue 1 items, all fixed here).**
* Added: `H_edge` row; `K₂`-beyond-criterion row; `H_RS` row; `H_K` replacing V(b), graded [D-sketch]; `c₁` as a free parameter.
* Regraded:
  * `H_LJ` regraded per RU issue 4;
  * "Flat [P]" → [P given H0 §17] + [NE];
  * the "Re" row now says the strip is load-bearing on odd rows at 1000;
  * H17 marked as carried over;
  * BF_SPINJ window evidence regraded to [M];
  * "first evidence of real-section ≡ crystal beyond LO" retracted (RU issue 2).
* Restated or withdrawn:
  * "Every computed state passes" restated with the EE-block qualifier and protocol labels;
  * the local vs min-norm difference at 300 corrected;
  * the AKT sentence softened;
  * "`|c₁| ≤ 1`" removed;
  * in-sample states labelled.
* `HYPOTHESES_CLOSE(2)`, `PARITY_ODD_CLOSE`, `HLJ_ANALYTIC`, `HWIN_THIMBLE`, `BELOW_FOLD_SPINJ`, `HRS_CHECK` and `RESULT_D_CLEANUP` are
  now incorporated.

**A.2 One value per state: superseded values and why.**

| state | quoted values in older documents | adopted value | reason |
|---|---|---|---|
| **`(300, J = 28)`** | `+1.45 ± 0.58` (EE_DG300_SMALLNU, 1.1Λ, pre-curvature fix; post-fix central `+1.570`); `+0.81` (curvature-fixed, `K₂` bulk); `+0.69` (REVIEW_EE S5, repeat thimble seed, `K₂` bulk, 1.1Λ); `+0.835` / `+0.794` (HC §5, direct bulk, main / repeat seed); **`+0.799 ± 0.011`** (HC2 §3, direct bulk + pooled validated `ν = 5` node); `+0.787` / `+0.739` (BF_SPINJ engine, a² / a¹) | **`+0.799 ± 0.011`** (HC2) | latest value with every input validated. `K₂` bulk is wrong at `J²/m = 2.35`; the main thimble seed is a `2⁹` outlier; BF_SPINJ is [M] (RU issue 2) and agrees within its calibration (±0.05) |
| **`(300, J = 32)`** | `+0.657` (previous statement §9.3, `K₂` bulk); **`+0.783`** (HC §5, direct bulk); `+0.865 ± 0.12` (HC2, pooled thimble seeds that fail the scalar check); `+0.834` (BF_SPINJ a¹, [M]); `−1.63` (HC2, failed thimble node, not a measurement) | **`+0.783`** (HC) | latest with validated inputs; window `ν < 8.3` [M] |
| binding state, 300 | `+0.2052` (`1.02Λ` grid point); `+0.2043` (hybrid at the twist edge, uncorrected); `+0.2038` | **`+0.2038`** (twist edge `1.0166Λ`, corrected; HC, CL) | the `+0.2052` is the grid point `1.02Λ`, not the minimum; `+0.2043` is the uncorrected value |
| binding state, 1000 | `+0.2021` (`1.01Λ` grid); `+0.2009` (hybrid, `1.00499Λ`) | **`+0.2009`** | the twist-edge state is the minimum; both reproduced by REVIEW_EE |
| `(300, J = 22)` | `+0.15` / `+0.123` (`t* = 1.1`) / `+0.046` (twist edge), all conv. A | **`+0.844`** (CL, direct + edge) | conv. A artefact (`K₂` at `x = 0.05`) |
| `(300, J = 24)` | 1.1Λ `+1.06 ± 0.14` (pre-fix) / `+1.013` (post-fix); min `+0.82` | **`+0.820`** (HC, with edge) | — |
| `(300, J = 18)` OO | `+0.222` (HC, conv. A) | **`+0.994`** (HC2) | spline aliasing + conv. A artefact |
| `(300, J = 15)` | `+0.707` / `+0.630` (conv. A) / `+0.829` (PO, sparse grid) | **`+0.7853`** (CL, dense grid + edge) | PO's grid missed the minimum |
| `(300, J = 9, 11, 13)` | `+0.438` / `+0.645` / `+0.728` (conv. A) | `+0.4506` / `+0.6687` / `+0.7746` (CL) | conv. A pessimistic |
| `(1000, J = 48)` | "smooth offset", [I]; `+0.932` (fixed-engine `K₂` row) | **`+0.915 ± 0.08`** (CL, direct sub-row) | `K₂` wrong at `J²/m = 2.2` |
| `(1000, J = 4, 6, 7, 11, 13, 16)` | `+0.2103` / `+0.2508` / `+0.3114` / `+0.4099` / `+0.6046` / `+0.7477` (grid from `1.01Λ`) | values of §5.2 (CL) | twist-edge masses not sampled before |
| `(1000, J = 29)` | `−1.20` (raw chain) | `+0.979` | isolated `K₂` outlier in a reference node; odd-row despike |
| `H_LJ` tolerances | uniform smooth `1–30 c₂`, fold window, `A_E` (SMOOTH_BOUND) | `H_AG ∧ H_win ∧ H_edge-amp` with `T_B` | HLJ_ANALYTIC §5.1 |
| `H_LJ(1000)` grade | [I] → [D + NE + M] (`RESULT_D.md`) | [D-sketch + M] | RU issue 4 |
| `H_OO` | leading order only, [M] | [D given `Γ_OO = Γ_++`] | HC2, PO |

**A.3 Withdrawn functionals and claims (unchanged from the previous statement §7, condensed).**
* The B′-uplift functional `Ω^{NF,SC}` and its "`n = 1`", and `Ω^U` on M02, are withdrawn. CFT light-ray sum rules map to the *stripped*
  flat rules, not to B′.
* The literal CLS-D `Ĉ₂` is withdrawn: its heavy normalisation vanishes at `ν₀ = 5.979`. The scalar-`g` `Ĉ₂` violates (HN).
* Route A′ with `Ĉ₂′` (`q²`) is withdrawn because of the leading-Regge helicity off-diagonal: near-gap failures at 300 for `J ≥ 16–18`,
  and modelled failures at 1000 from `J ≈ 32`.
* "`Δ_g ≥ 300`" and "certified window `b̂ ∈ [91.2, 2Λ]`" as graviton statements are withdrawn; the latter concerns the scalar skeleton
  model.
* Bugs and artefacts that were fixed:
  * scipy `eval_jacobi` for complex arguments;
  * the curvature double count (λ was too high by `≤ 0.0135` at 1000 and up to `0.07` at 300);
  * the spin-`J` engine precision loss `ε²4^J`;
  * undespiked odd rows;
  * the canonical gauge reversed to min-norm.

**A.4 Contradictions between documents (status).**
* **Resolved here.**
  * Binding-state values (A.2).
  * `J = 28` and `J = 32` values (A.2).
  * `H_LJ` grade (RU wins).
  * The "512 c₂ `K₂` discrepancy" of RU issue 3 is the `k = 2` (not `k = 4`) comparison-script artefact `log(ΣM_k)`, resolved in BF3 §2:
    `K₂` = engine to 2 c₂. It is positive evidence for `H_RS` at `ν = 20–40` (HRS_CHECK §5). RU's description of it is superseded; its
    grade for `H_RS` ([I] where used) stands.
  * `K₂` convergence criterion: "`J²/m ≲ 1.5`" (previous statement) vs "`≲ 2.2`" (`EE_DG1000_LARGEJ`, as quoted by CL). CL shows 2.2
    fails. Adopted: validated to 1.73; 1.6–2.2 untested.
* **Unresolved.**
  * CL §0.2 says `H_edge` is computed at 1000 "for `J ≤ 17`", but its table has no row for `J = 15` (nor `J = 0, 2`), and its open-item
    list does not mention them. Here `J = 15` (1000) is treated as edge not computed, on the `1.01Λ` grid (`+0.869`).

**A.5 Older reproduction notes.** The reproduction recipes of the previous statement §9 (odd `J = 19–31` at 1000; `J = 32` at 300; odd-row
despike check) are in git (`c473d9e`). The recipes for the current values are in HC §6, HC2 §5, PO, CL §6 and HRS_CHECK §8.
