# Einstein gravity through quartic order in AdS₄/CFT₃ with a large gap: a conditional theorem

**Ansh Mishra and Aryan Senthilkumar**

*22 Nov 2026. Unrefereed. This document unifies Result D (`RESULT_D_STATEMENT.md`) and Result Q (`RESULT_Q_STATEMENT.md`). Its starting
point is `REVIEW_UNIFIED_2026-11.md` (RU) §5's recommended wording, updated with `HRS_CHECK.md`, `RESULT_Q_OG_CHECK.md`,
`RESULT_Q_FERMION_HEFT.md` and `RESULT_D_CLEANUP.md`. No new computation was made for it, and every number is quoted from the two
statements, which name its source.*

**This is a conditional theorem.**
* The implication "per-state positivity of three fixed functionals ⇒ the bounds" is [D], given a list of analytic inputs, themselves
  graded [Q], [D], [D-sketch] or [A].
* **Per-state positivity is not established in general.** It is computed [NE + M] on finite sets of heavy states and assumed elsewhere,
  through named hypotheses, two of which are open [I].

**Grade legend.**
* **[P]**: proved. **[certified]**: computer-assisted proof.
* **[D]**: derived. **[D-sketch]**: derivation sketched. **[DV]**: derived and checked numerically.
* **[Q]**: quoted from the literature.
* **[NE]**: numerical evidence. **[NE + M]**: computed inside the near-gap chain/HYB model.
* **[M]**: model. **[A]**: assumption. **[I]**: open.

A chain carries the grade of its weakest link.

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

> **Update (DG_SCAN.md, unreviewed): Δ_g dependence.**
> * **Bound values as explicit functions of Λ [DV].** The functionals are held fixed, and the forms reproduce every quoted digit at 100, 300 and 1000:
>   * cubic: `√(t₄²+α_T²) ≤ 1.3194·10⁴·√(ln(Λ/8.7415) + 6.384/Λ + O(lnΛ/Λ²))/Λ⁴`;
>   * quartic upper end: `17.2407·ln(2Λ/3) − 19.2782 + 10.314/Λ + (860.3·ln(2Λ/3) − 2508.6)/Λ² + …`.
> * **Where the constant drift comes from.** The shift −19.220 → −19.265 is the finite-Λ tail of the AdS pole, and the constant tends to the flat value −19.2782.
> * **Error terms in Λ:**
>   * c₁ (Result D): `∓19941/(Λ²√(16L_A))` per unit c₁;
>   * H17: `±212.5/Λ²` [M].
> * **Positivity threshold not computed.** The model bracket is `140 ≲ Δ_g* ≲ 300` [M]. The sign-free large-J window closes only from Δ_g ≈ 275, with margin ×1.1 at 300. So 300 is the lowest gap where the current method closes.
> * **Proposed corollary.** Under H_mono (corrections fall monotonically at fixed (J, m/Λ)) and H_cont (continuity in m/Λ), the bounds hold for all `Δ_g ≥ 300` [D-sketch + M]. Nothing is claimed below 300.
> * **Q's binding state is Δ_g-stable.** The law `s_d ≈ 0.50 + 1440/Δ²` was mis-fitted (REPRO D1). The cached data give `s_d ≈ 0.5076 + 147/Δ²` [NE/M], which is more stable than stated.
> * **Record correction.** RESULT_D_STATEMENT's c₁ entry at Δ_g = 100 should read `∓2.37 %`, not `≈ ∓1.7 %`, by the same formula as the 300/1000 entries.

> **Update (REPRO.md, independent reproduction with fresh code): every headline number reproduces within its stated precision.**
> * **Bounds and constants (fully independent):**
>   * cubic bound `3.13206·10⁻⁶` / `2.89126·10⁻⁸`;
>   * `2L_AdS` `72.0406` / `92.8135`;
>   * `B_eff` and the DG_SCAN Λ-forms.
> * **Error terms (fully independent):** c₁, H17 and δk(0).
> * **Binding states:**
>   * Q binding FF/N `+0.14270` / `+0.14282`. Partly independent: the engine ratio ρ(ν), the de-shift weights χ and `R_j` were taken as given.
>   * D binding states: hybrid `+0.20427` (300) and `+0.20090` (1000).
> * **ABJM:** placement `3.34 %` / `2.59 %`, with `8πG = 256/c_T` confirmed.
> * **Not independently reproduced:**
>   * the analytic-radial engine ρ(ν), which carries all of `s_d`;
>   * the Abel de-shift weights;
>   * D's crystal corrections;
>   * the NF02 flat certificate (grid check only);
>   * the full 14000-state Q scan.
> * **Record fixes applied:**
>   * H17 band `±0.24 %` (was `±0.2 %`);
>   * twist edge `1.00499Λ` (was `1.0048Λ`);
>   * the `s_d` vs Δ_g law corrected.
> * **Remaining minor items (none binding):** D4, D5, D7 and D9 in REPRO §2.

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

## 1. Setup and normalisations

| item | definition |
|---|---|
| theory | `𝒯`: a unitary CFT₃ (parity-even or parity-violating) with a unique stress tensor `T` |
| `c_T` | normalised so that a free real scalar has `c_T = 1` |
| Newton constant | `8πG := 256/c_T` (AdS units `L = 1`; from `c_T = 32L²/(πG)`) |
| **A1** (gap) | every single-trace operator other than `T` has twist `Δ − J ≥ Δ_g`, for all `J` (including `J = 0, 1, 2`) |
| **A2** (tree level) | `O(1/c_T)` truncation, `c_T ≫ Δ_g²`; loops dropped |
| gap scale | `Λ := √((Δ_g − 2)(Δ_g − 1))`: `298.50` (`Δ_g = 300`), `998.50` (`Δ_g = 1000`) |
| heavy kinematics | `m² = (Δ − J − 2)(Δ + J − 1)`; gap state `m_gap² = (Δ_g − 2)(Δ_g − 1 + 2J)` |
| cubic couplings | `⟨TTT⟩` in the basis `(c_T, t₄, α_T)`; bulk dictionary `√(t₄² + α_T²) = 240|ĝ₃|` (AdS units; `ĝ₃` in CHLPS normalisation) |
| quartic coupling | `g₄`: flat R⁴ coefficient in Stage 0's normalisation (type II Virasoro–Shapiro: `g₄/(8πG) = 2ζ(3)/M⁶`, `M² = 4/α′`). `g̃₄`: the CFT-defined light-side coefficient `−Ĉ₄^{ds}|_light = 2g̃₄` after Lemma A; `g̃₄ = g₄` given (F)+(P). |
| `c₁` | coefficient of a possible `ν`-independent `E × W³` light term `16πG c₁t₄`, collecting the `k = 2` (`Ĉ′_EE`) and `k = 3` (`Ĉ₃`) parts. **Not computed; a free real parameter.** |

RU §1 checked that the two results use the same conventions: gap, `Λ`, `G`, twist convention, and the `k = 2` object `Ĉ′_EE` in the
min-norm gauge.

**The fixed functionals** (fixed as of 22 Nov 2026; not to be changed).
* **`Ω^EE`** (Result D). It is built from CLS-D's light-ray functional (2311.04271):
  * `Ω^EE = ∫₀¹dp ψ₂ Λ⁴Ĉ′_{EE,Λp} + ∫₀¹dp ψ₃ ½Λ⁶Ĉ_{3,Λp}`;
  * the `k = 2` part uses the dual structure `EE`, de-shifted by `c(ν)⁻¹𝒜_g⁻¹`, in the **min-norm gauge**;
  * the `k = 3` part `Ĉ₃` (`o_H·q`, `j′ = 1`) has no de-shift;
  * NF02 weights `ψ₂, ψ₃` from `code/robust/functionals_noforward.json`, with `k₀ = 27.48509817` and `k(0) = 1`.

  Full definition: `RESULT_D_STATEMENT.md` §2.
* **`Ω^Q_up`** (Result Q, upper).
  * The δ = 0.2 functional `Q3-full-delta0.2` of `code/resultq/functionals_q3.json` (`A = 17.241`, `κ = 0.1099`).
  * It is realised by the de-shifted `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}` with the Stage 7 weights (`functionals_q7_xprec.json`), plus the
    smeared `k = 2, 3` parts with Q's weights.
* **`Ω^Q_low`** (Result Q, lower): the de-shifted `k = 4` value piece (`RESULT_Q_STATEMENT.md` §1, §2.4).

---

## 2. The theorem

> **Theorem (conditional; Einstein gravity through quartic order at `Δ_g ∈ {300, 1000}`).**
> Let `𝒯` satisfy A1 and A2 (§1). Assume the analytic inputs of §3:
> * the KKSDZ Regge bound with its loophole, and K1′;
> * `H_K`;
> * Adm, H7, L0;
> * Lemma A with D_LR;
> * (F)+(P);
> * `H_EFT`;
> * existence of the de-shift;
> * H17.
>
> Assume also **per-state positivity** of `Ω^EE`, `Ω^Q_up` and `Ω^Q_low` on every twist-allowed heavy block. Positivity is computed
> [NE + M] on the finite sets of §3.2. Elsewhere it is assumed through `H_M`, `H_R`, `H_edge`, `H_OO`, `H_RS`, `H_LJ(Δ_g)`, `H_LJ^Q`,
> `H_fwd^Q`, `H_smear^Q` and HYB. Let `c₁` be the (uncomputed, free) coefficient of the `E × W³` light term. Then:
>
> 1. **(cubic)** `√(t₄² + α_T²) ≤ 240·√(16L_A/(k(0)(1+δ)))/Λ⁴ · (1 + O(c₁Λ⁻²))`, with `L_A = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`. In AdS
>    units:
>    * **`≤ 3.13·10⁻⁶` at `Δ_g = 300`** and **`≤ 2.89·10⁻⁸` at `Δ_g = 1000`**;
>    * `∓0.22 %` / `∓0.017 %` per unit `c₁`, and `±0.24 %` / `±0.021 %` from H17.
>
>    Equivalently, `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴ / 1.43·10⁴`.
> 2. **(quartic)** Add one of two hypotheses:
>    * **`κ ≥ 0`** [D-sketch + M] for the upper end and **`c_L > 0`** [D-sketch + M] for the lower end;
>    * or simply **`t₄ = α_T = 0`**.
>
>    Then
>
>    `−O(Λ⁻²) ≤ g̃₄Λ⁶/(8πG) ≤ 17.241·log(2Λ/3) − 19.220 (Δ_g = 300) / − 19.265 (Δ_g = 1000) + ε_{c₁}`,
>
>    i.e. the upper end is `72.04 / 92.81 + ε_{c₁}`, with `|ε_{c₁}| = 2.03·|c₁t₄|Λ²`.
>    * Using item 1 to bound `t₄`: `|ε_{c₁}| ≤ 0.57|c₁|` (300) and `≤ 0.059|c₁|` (1000).
>    * For `ĝ₃ ≠ 0` and `κ ≥ 0`, the term `−κ|ĝ₃|²Λ⁸` (`κ = 0.110`) only tightens the upper end.
>    * The CFT-defined `g̃₄` equals the bulk `g₄` given (F)+(P).
>
> **Grade.**
> * The implication "positivity ⇒ (1), (2)" is **[D]** given the analytic inputs; the constants are [DV].
> * **Positivity itself is not established.** It holds at the grades of §3: two items are [I] (`H_RS` for Q's binding state; the Regge
>   window of `H_LJ(300)` for `J ≥ 34`), and the analytic inputs `H_K` and (F)+(P) are [D-sketch].
> * The Einstein instance of Lemma A has not been checked.
>
> Nothing is claimed at any other `Δ_g`, at `Δ_g = 100`, at loop level, or for the ratio `|ĝ₃|²/g₄` (no cubic closure from (2)).

**Remarks.**
1. **Dependence structure.**
   * Item 1 uses only `Ω^EE`, and so does not need Lemma A, (F)+(P), `H_EFT` or anything forward.
   * Item 2 needs `Ω^Q_up` and `Ω^Q_low`. It needs item 1 only through `ε_{c₁}`.
   * If item 1 is used to bound `ε_{c₁}`, item 2 at `Δ_g = 300` inherits `H_LJ(300)`, `H_M`, `H_R`, `H_edge` and `H_RS`. Carrying `c₁t₄`
     explicitly avoids that.
2. **Two deviations from RU §5's wording.**
   * The lower end for `ĝ₃ ≠ 0` additionally needs the sign `c_L > 0` of the lower functional's `|ĝ₃|²` slope (measured 0.255–0.49,
     `RESULT_Q_REVISION` §3.3). RU's wording covered only `κ`.
   * `B_eff(1000)` is written `−19.265`. The value `−19.2647` (REVIEW_RESULTQ §1) has been rounded to "−19.27" in earlier documents; to two
     decimals it is −19.26. The effect is `< 0.01` on `92.81`.
3. **Where the margins are.**
   * D's binding states: `+0.2038` (300) and `+0.2009` (1000) in units of `N_s`, against a functional built with δ = 0.2.
   * Q's binding pair (`++`/OO J = 0 at the gap): FF/N `+0.1426 ± 0.0017`, with derivative-piece margin `s_fail/s_d = 1.78/0.51 ≈ 3.5`.

---

## 3. Consolidated hypothesis table

"D"/"Q" mark which result needs the item. Grades follow RU §4, updated by HRS_CHECK (HRS), RESULT_Q_OG_CHECK (OG),
RESULT_Q_FERMION_HEFT (FH) and RESULT_D_CLEANUP (CL).

### 3.1 Setting and analytic inputs

| # | hypothesis | D | Q | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|---|
| A1 | twist gap `Δ − J ≥ Δ_g`, all `J`; unique `T` | ✓ | ✓ | [A] | setting | — |
| A2 | tree level, `c_T ≫ Δ_g²` | ✓ | ✓ | [A] | setting | — |
| K1 | KKSDZ Regge bound (1904.05905), `J₀ ≤ 1` for spinning correlators | ✓ | ✓ | [Q], with KKSDZ's loophole (first/second-sheet boundary, their §4.2.6) | `H7_I14_LIGHTSIDE` §2 | closing the loophole |
| K1′ | uniformity of K1 towards the light-cone edges | ✓ | ✓ | [D-sketch] | review V | write it out |
| K2, K3 | kernel power counting; per-family Rindler positivity ⇒ swap lemma; validity `J₀ < 7/2` of the de-shifted `k = 4, 5` functionals | ✓ | ✓ | [Q]/[D]; `J₀ < 7/2` [D-sketch] | `H7_I14` §7; V4FINAL §5 | — |
| **H_K** | heavy OPE sum ↔ (corner, `ν`) interchange for the heavy-normalised de-shifted kernels (`Ĉ′_EE`; `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`) | ✓ | ✓ | **[D-sketch]** ([D given `H_ip`]: raw per-`ν` heavy-sum convergence, `J₀ < ½`, not implied by KKSDZ) | RU issue 6; V4FINAL §2.4 | Proposition K for the heavy-normalised kernel, or `H_ip` |
| Adm | `EE` is an admissible `j′ = 0` dual structure | ✓ | ✓ | [D] | `EE_FOUNDATIONS` §1 | — |
| H7 | `−Ĉ′_EE|_E = 16πG/(ν² + 9/4)`, `Ĉ₃|_E = 0` | ✓ | ✓ | [D] | `H7_I14` §3–§5, §18–§23 | an independent small-`ν` check of `h₋ + h₊` |
| L0 | `[TT]` double zeros at `k = 2, 3`; contacts drop out | ✓ | — | [D] | `H0` §15; `PARITY_ODD_ROUTEA` §12 | — |
| **Lemma A + D_LR** | tree terms with `J₀ < 4` (Einstein, `E × W³`, MFT) are annihilated pointwise by `Ψ₄`, `Ψ₅`; the collinear set carries nothing | — | ✓ | **[D]**. MFT instance [P]. Mechanism **[DV] at `O(1/C_T)`** on the exact free-scalar and free-fermion `⟨TTTT⟩` (residual `≤ 8·10⁻¹⁵`; controls fail at O(1)). Also **[DV] on spinning Witten-diagram correlators** with a `J₀ = 2` spin-2 exchange (EIN). **Einstein `⟨TTTT⟩` itself not evaluated [I].** | V4FINAL §2.2; OG; FH; EIN | one swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` at one `ν` |
| **(F)+(P)** | spinning flat limit and polynomial light side (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`, `g₅` cancellation) | — | ✓ | [D-sketch] | V4FINAL §5 | a `J₀ = 4` spinning position-space correlator |
| **H_EFT** | (a) term-wise (P) of heavy-exchange double traces = their actual value; (b) no conspiracy in the heavy sum / contact tower | — | ✓ | (a) [DV] in the scalar-Mellin toy (`J_h = 0–6`), [D-sketch] spinning; **(b) [A]** | FH §4.5 (refines RU issue 9 and OG §5.2) | (a) a spinning toy; (b) — |
| De-shift | continuum C4/C5 solutions exist for all `β` | (C2: [D]) | ✓ | [NE] to `t = 30`; C5 spinning `j′ = 1` coefficients [I] | Stage 7; V4FINAL §5 | — |
| H17 | `|ĝ₃|²` light coefficient `k(0)(1+δ)` | ✓ | (κ) | [D-sketch + M], carried over from `q²` | `H0` §18; REVIEW_EE m4 | spinning `R³` exchange to `O(r²)` |
| **κ ≥ 0** | sign of Q's upper `|ĝ₃|²` coefficient (`κ = 0.110`) | — | ✓ if `ĝ₃ ≠ 0` | [D-sketch + M] | REVISION §3.3; RU issue 7 | spinning W³W³ |
| c_L > 0 | sign of the lower functional's `|ĝ₃|²` slope (0.255–0.49) | — | ✓ if `ĝ₃ ≠ 0` | [D-sketch + M] | REVISION §3.3 | same |
| **c₁** | `E × W³` light term (`k = 2` and `k = 3`) | ✓ | ✓ | **free parameter, carried explicitly** | RU issue 14 | Regge residue of the `E × R³` `T` exchange |

### 3.2 Per-state positivity: what is computed, and the hypotheses that cover the rest

| # | hypothesis | D | Q | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|---|
| — | **computed set, D** | ✓ | — | [NE + M]. 1000: even `J ≤ 24`, 28, 32, 36, 40, 48; odd `J ≤ 31`; min `+0.2009`. 300: even `J ≤ 24`, 28, 32; odd `J ≤ 15`; min `+0.2038`. One value per state: `RESULT_D_STATEMENT` §5. | CL; HC; HC2; `EE_DG1000_LARGEJ` | — |
| — | **computed set, Q** | — | ✓ | [NE + M] given `H_RS`. Measured: `++` J = 0, 2, 4, 6, 8; OO J = 0, 2; `+−` J = 4–9, 11. Every other channel needs `|s_d| ≥ 17–30×` the measured `0.51` to fail. | BF2; BF3; PO; FWD | — |
| **H_RS** | analytic-radial engine image = crystal-cycle integral for EE beyond leading order | (window 28/32; OO beyond order 0) | ✓ (every `s_d`) | **[I]** where it matters: `ν < 6` for Q's binding state (`ν < 4` carries 89 % of the `s_d` weight; `ν < 2` is untested). [NE] at `ν ≥ 10`; [NE, weak] at 6–10; scalar and EE LO [DV]. At `J = 28` the entry-resolved test fails at 3–16σ. | HRS (refines RU issue 3); RU issue 2 | an EE `k = 4` crystal value at `J = 0`, `ν < 4`, to ~10–20 per harmonic |
| H_M | the near-gap chain is accurate on the computed states (incl. `K₂` validity) | ✓ | — | [M], validated [NE] on rows with direct data. `K₂` fails at `J²/m ≥ 2.2`; all such states are now on direct sub-rows (incl. `(1048, 48)`: `+0.915 ± 0.08`). The 1.6–2.2 band is untested. | HC §5; CL §0.3 | direct sub-rows in the band; interval chains |
| H_edge | light-cone-edge window `x ∈ (0.8, Λ/m]` | ✓ | (tolerance `≥ 760 c₂`) | [NE + M] where computed (300: `J = 4–24`; 1000: `J = 4–14`, 16, 17, 31, 36, 40); [M] elsewhere | HC §2; CL §0.2 | the remaining rows |
| H_OO | parity-odd even-`J` block | ✓ | ✓ (J = 0, 2) | [D given `Γ_OO = Γ_++`]. The identity is [P given the helicity lemma] at order 0 and [NE] beyond (`≤ 0.05 / 0.12 c₂`; `|Δs_d| ≤ 5·10⁻⁶`). It reduces to EE `++` and `H_RS`. | HC2; PO; FWD; RU issue 11 | a derivation beyond order 0 |
| H_R | beyond-leading-Regge at `m ≳ 1.3Λ`, large `b̂` | ✓ | (HYB) | [M] | `EE_NEARGAP` §6; `ROUTEA_LARGE_B` | spinning `O(r)` Regge kernel |
| **H_LJ(1000)** | near-gap large spin: `H_AG ∧ H_win (refined) ∧ H_edge-amp` | ✓ | — | **[D-sketch + M]**. The sign-free `H_win` is violated at `(1100, 100)`; the refined form rests on HLJ §5.2 [D-sketch]. | HLJ_ANALYTIC; HWIN_THIMBLE; RU issue 4 | spinning `O(r)` kernel at `β ≤ 0.53` |
| **H_LJ(300)** | same at 300, plus the `(336, 28/32)` window `ν < 8.3` | ✓ | (only if `t₄` is imported) | **[D-sketch/M] bulk; [I] Regge window for `J ≥ 34`** (no graviton data); `(336, 28/32)` window [M] (no validated direct datum below `ν = 8.3` at `J = 32`) | RU issue 4; HC2 §3; RU issue 2 | graviton window data at 300; a validated below-fold method |
| H_LJ^Q | large-`J` near-gap corrections below Q's tolerance (`≥ 1.1·10⁴ c₂`) | — | ✓ | [M], margin ×10–10³ | RU issue 13 | `s_d` at larger `J` |
| H_fwd^Q | unmeasured forward channels; the 300 J-drift (J ≥ 7 quoted from 1000) | — | ✓ | [NE + M] given `H_RS`; drift [I], harmless at δ = 0.2 | FWD; RU issue 12 | — |
| H_smear^Q | smeared graviton corrections at the binding states (`≤ 3·10⁻⁵` of B₂ vs tolerance `≈ 24 %`) | — | ✓ | [NE + M] | REVISION §1 | — |
| HYB | scalar skeleton + leading Regge elsewhere; `β > 30` | (hybrid) | ✓ | [M] | RU §4 | — |
| σ_hU | double counting of the helicity-odd correction | — | δ = 0.05 and cubic closure only | [I]; **not used by the theorem** | RU issue 8 | — |

---

## 4. What is proved, derived, numerical, assumed

* **Proved / certified / derived and checked.**
  * The constants `A = 3022.34`, `B_AdS = −6552.66`, `16L_A = 10734.61 / 14339.59` and the cubic bounds [DV]. They were reproduced
    independently by REVIEW_EE, and recomputed for this note from `16L_A` and `Λ`.
  * `2L_AdS = 72.041 / 92.814` [DV] (RU §1).
  * NF02's flat positivity certificate [certified] (one float-checked corner).
  * The crystal identity `EE ≡ q²` on CLS-D's physical polarisation plane [P]. Vanishing of `Ω^EE`'s leading-Regge helicity
    off-diagonal [P].
  * On the light-transform manifold `W_u·f₄ ≡ 0` [P]. The MFT instance of Lemma A [P].
  * The Lemma A mechanism on the free scalar and the free fermion at `O(1/C_T)` [DV].
  * The crystal saddle sits on a helicity eigen-polarisation [P], so `Γ_OO = Γ_++` at order 0 [P given the helicity lemma (D, DV)].
  * The ABJM arithmetic and the `c_T ↔ G` normalisation [DV].
* **Derived [D].**
  * The two implications (theorem items 1 and 2) given their inputs.
  * Adm, H7, L0; Lemma A + D_LR (Einstein instance unchecked); V(a).
* **Derived in sketch [D-sketch].**
  * `H_K`; K1′; validity `J₀ < 7/2`; (F)+(P); spinning `H_EFT`(a); H17, `κ ≥ 0` and `c_L > 0` ([D-sketch + M]).
  * The lower end `g̃₄ ≥ −O(Λ⁻²)`.
  * The refined `H_win` and the Regge-window structure behind `H_LJ`.
* **Numerical [NE], [NE + M].**
  * Per-state positivity on the computed sets (§3.2).
  * `H_edge` where computed; the `Γ_OO = Γ_++` identity beyond order 0.
  * `H_RS` at `ν ≥ 10`; `H_smear^Q`; `H_fwd^Q`; de-shift existence to `t = 30`; `H_EFT`(a) in the toy.
* **Model [M].** `H_M`, `H_R`, HYB, `H_LJ^Q`, `H_edge` on uncomputed rows, the `(336, 28/32)` window, and the large-spin and
  large-`b̂` regions.
* **Quoted [Q].** The KKSDZ Regge bound, with its loophole. CLS-D's definition of the light-ray functional.
* **Assumed [A].** A1, A2, `H_EFT`(b).
* **Open [I].**
  * `H_RS` below `ν ≈ 6` for Q's binding state.
  * The `H_LJ(300)` Regge window for `J ≥ 34`.
  * The Einstein instance of Lemma A (derived, never checked).
  * The C5 spinning de-shift coefficients.
  * `c₁` (carried as a free parameter rather than assumed).

## 5. Weakest links, ranked

1. **`H_RS` at `ν < 6` for Q's binding state (`++`/OO J = 0 at the gap)** [I]. Result Q.
   * `ν < 4` carries 89 % of the `s_d` weight; `ν < 2` is untested; no available crystal-cycle method reaches the needed precision.
   * A ×3.5 engine error in the derivative piece would break the δ = 0.2 functional at the gap.
   * HRS_CHECK confirmed the engine at `ν ≥ 10` and found nothing against it, but cannot test where it matters.
2. **The `H_LJ(300)` Regge window for `J ≥ 34`, and the `(336, 32)` window `ν < 8.3`** [I] / [M]. Result D at 300, and Q at 300 if
   `t₄` is imported.
   * There are no graviton window data at `J ≥ 34`.
   * The sign-free tolerance (`180–390 c₂`) is of the size of the only measured window values (`≈ 450` at `J = 24, 28`). It holds with
     margin only if the window correction stays leading-Regge shaped [D-sketch].
   * `J = 32` has no validated datum below `ν = 8.3`.
3. **The Einstein instance of Lemma A** [D, unchecked]. Result Q.
   * Its leverage on `g₄` is `10¹⁰–10²²`.
   * The mechanism is [DV] on two free theories, and the controls show what the zero depends on: primarity of `T`, the light-transform
     weight, the T-shift, microcausality. Einstein satisfies all of these [D].
   * But Einstein's (13) channel (`T` + `[TT]` with anomalous dimensions) is probed only by the free scalar, because the fermion is
     blind to it in this gauge.
4. **`H_EFT`(b), no conspiracy in the heavy sum / contact tower** [A]. Result Q. Untested; the term-wise series is at best asymptotic.
5. **`H_K`** [D-sketch]. Both results. The main analytic hypothesis; it has no criterion short of `H_ip`.
6. **(F)+(P)** [D-sketch]. Result Q. It also carries the identification `g̃₄ = g₄`. There is no position-space R⁴ normalisation control.
7. **`H_LJ(1000)`** [D-sketch + M]. Result D at 1000. The sign-free `H_win` is violated at `(1100, 100)`; the refined form's class
   tolerance is computed at 1000 only; `ν < 10` is unmeasurable.
8. **KKSDZ's loophole and K1′** [Q] + [D-sketch]. Both.
9. **`κ ≥ 0` and `c_L > 0`** [D-sketch + M]. Q, only for `ĝ₃ ≠ 0`. Not needed if `t₄ = α_T = 0`.
10. **Model-level items** [M]: `H_M` (the `K₂` band `J²/m = 1.6–2.2`), `H_R` (the `N₀` pathology at 300), `H_edge` on uncomputed rows,
    HYB, `H_LJ^Q`, and the C5 spinning de-shift coefficients [I]. RU does not list the last among the weakest, but its grade is [I].

## 6. What is NOT claimed

* **No other `Δ_g`.** Only 300 and 1000 were computed; "`Δ_g ≥ 300`" is not claimed.
* **Not `Δ_g = 100`.** Positivity is contradicted at model level at near-gap `J ≳ 18`.
* **No loop level** (A2).
* **No cubic closure** `|ĝ₃|²Λ⁸ ≤ …` from item 2, and no ratio `|ĝ₃|²/g₄`.
  * The earlier "`120–198 / 154–257`" is withdrawn. Its lower slope 0.49 is 1.9× CHLPS's flat optimum 0.26, inconsistent with the flat
    limit (RU issue 8).
  * At most a remark: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257` [D-sketch + M].
* **The Einstein instance of Lemma A is unchecked.** Only free-field proxies at `O(1/C_T)` were evaluated.
* **Per-state positivity is not established in general**, and neither result is a theorem without it.
* **Not the δ = 0.05 functional.** Its margin `+0.0060` depends on `σ_hU`.
* **Nothing about theories that violate A1**, e.g. ABJM (§7).
* **Novelty only at search level** [I].

## 7. ABJM consistency check

* **Dictionary** (Stage 0; re-checked by RU §1):
  * `M² = 4/α′` and `L²/α′ = π√(2λ)`, so `Λ ↔ ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.22 λ^{1/4}`. The identification holds only at leading order
    in `1/(ML)`.
  * The type IIA `R⁴` term `ζ(3)α′³/32` gives `g₄M⁶/(8πG) = 2ζ(3) = 2.404`, independent of `λ`.
* **Placement.** `2ζ(3)` sits at **3.3 %** (`Δ_g = 300`) and **2.6 %** (`Δ_g = 1000`) of the upper end of item 2 (`2.404/72.04`,
  `2.404/92.81`).
* **`ĝ₃ = 0`.** This follows from **10d type IIA supersymmetry** (no `R³`), not from "maximal SUSY": ABJM at `k > 2` has `N = 6`.
* **No `c_T` confusion.** `c_T` does not enter the placement, because `2ζ(3)` is `λ`-independent and `Λ` depends only on `λ`.
  * The mass-derivative identity `c_T = (64/π²)∂²_m F` relies on the mass coupling to a current in the stress-tensor multiplet. It is used
    only in `ABJM_THEOREM`'s separate `c_T` result, with the same `c_T` (free real scalar `= 1`, `8πG = 256/c_T`).
  * In an `N = 2` theory `∂²_m F` gives a **flavour-current** coefficient, not `c_T`. The two must not be confused.
* **Status.** ABJM violates A1 (light Kaluza–Klein modes). The check uses only the flat image (F), and the bound is log-dominated.
  **A consistency check, not evidence for the theorem.**

## 8. Literature positioning

* **CLS-D (Chang–Landau–Simmons-Duffin, arXiv:2311.04271, v3 Feb 2025; JHEP 04 (2025) 016).**
  * They construct the spinning light-ray dispersive sum rules on `⟨TTTT⟩` and the dictionary to flat-space graviton sum rules in the
    bulk-point limit.
  * They list "the remaining steps for uplifting flat space bounds to AdS" as **future work**.
  * Results D and Q are precisely that programme. They are presented as **a first explicit-constant, conditional execution of CLS-D's
    uplift** for the cubic and quartic graviton couplings.
  * Priority on the method is CLS-D's. The new ingredients are:
    * the `EE` structure, which removes the curvature-induced helicity off-diagonal;
    * the de-shifted `k = 4, 5` functionals;
    * the pointwise Lemma A;
    * the near-gap computations.
* **CHMRS (Caron-Huot–Mazáč–Rastelli–Simmons-Duffin, arXiv:2106.10274)** is the **methodological precedent**: two-sided
  `log(Δ_gap)` bounds on contact couplings in AdS₄ with dynamical gravity, from uplifted flat sum rules, for scalar correlators. Q's
  light-side formula and `H_EFT` are of CHMRS's (5.3) type, and D's de-shift uses CHMRS's scalar kernel `g`. Credit is due to both
  CLS-D and CHMRS.
* **CHLPS (arXiv:2201.06602)** gives the flat D = 4 graviton bounds.
  * D's CHLPS-form rewriting, `|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21`, is about 10× weaker than the transplanted CHLPS optimum.
  * That optimum would need forward and helicity-flip sum rules in AdS, which are absent; it is asserted by CHLPS, not derived in AdS.
* **Parametric `Δ_g⁻⁴` bounds** on `⟨TTT⟩` deviations: CEMZ (1407.5597), Meltzer–Perlmutter (1712.04861), Costa–Hansen–Penedones
  (1707.07689), AKT (1805.07393; `t₄ ≲ ln Δ_g/Δ_g⁴`, unknown constant). Item 1 has a different parametric form, `√(ln Δ_g)/Δ_g⁴`, and
  explicit constants; the constants are not comparable.
* **To be discussed:** Kundu (2104.11238); the `N = 6` bootstrap (Binder–Chester–Jerdee–Pufu, 2011.05728).
* **Novelty [I] at search level.** No published bound on `R³`/`R⁴` couplings from `⟨TTTT⟩` sum rules was found (abstract-level search
  of 19 Nov plus one extended search of 22 Nov).

## 9. Open problems

1. **`H_RS` below the fold.** An EE `k = 4` crystal value at `J = 0`, `ν < 4`, to ~10–20 per harmonic. Routes: lower `Δ` with a
   higher-order engine, variance control for the Lemma-B factor `R`, or an analytic proof that the Γ-continuation of the `R`-weighted
   monomials equals the crystal-cycle integral.
2. **A validated below-fold method at large `J`.** This would give graviton window data at `Δ_g = 300`, `J ≥ 34`, and `(336, 32)` at
   `ν < 8.3`. The thimble sums currently fail the scalar `Im` check.
3. **The spinning `O(r)` Regge kernel** at near-gap `β ≤ 0.53`. It would make `H_win` [D] at both `Δ_g` and remove most of `H_R`.
4. **One swapped-form `Ψ₄` on Einstein `⟨TTTT⟩`** (O(G) `[TT]` data in every helicity structure). Also the (13)-channel/corner part of
   D_LR on a correlator with a spin-1 or spin-2 single-trace exchange (the free fermion is blind to it).
5. **`H_K`** (Proposition K for the heavy-normalised kernel) and the KKSDZ loophole.
6. **(F)+(P) in position space**, including an R⁴ normalisation control from a `J₀ = 4` spinning AdS contact diagram. Spinning
   `H_EFT`(a); any handle on `H_EFT`(b).
7. **`c₁`:** the Regge residue of the `E × R³` `T` exchange (`k = 2` and `k = 3`). It may vanish at leading Regge for `EE`.
8. **Spinning W³W³** (firms up `κ ≥ 0`, `c_L > 0`); `σ_hU`; the 300 J-drift of `s_d`.
9. **H17 for `EE`** (the `Π̂₄` admixture of `Ĉ′_EE`).
10. **`H_edge` on the remaining rows; `K₂` in the `J²/m = 1.6–2.2` band; interval versions of the near-gap chains.**
11. **Other `Δ_g`**, and `Δ_g = 100` in particular.
12. **Release checklist of RU §5:**
    * statements rewritten: done (22 Nov);
    * `H_RS` check: done, inconclusive at `ν < 6`;
    * spinning O(G) Einstein `Ψ₄`: not done (free-field proxies only);
    * direct sub-row for `(1048, 48)`: done (`+0.915 ± 0.08`);
    * BF_SPINJ entry-by-entry residuals: computed by RU itself, regrade adopted. The BF_SPINJ authors have not re-reported.

---

### Pointers

* Result D: `RESULT_D_STATEMENT.md` (one value per state, history of superseded values).
* Result Q: `RESULT_Q_STATEMENT.md`.
* Review: `REVIEW_UNIFIED_2026-11.md`.
* Checks: `HRS_CHECK.md`, `RESULT_Q_OG_CHECK.md`, `RESULT_Q_FERMION_HEFT.md`, `RESULT_D_CLEANUP.md`.
