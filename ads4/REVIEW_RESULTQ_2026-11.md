# Hostile referee report on Result Q (AdS₄/CFT₃ two-sided bound on the R⁴ coupling g₄), part 1

*21 Nov 2026, referee agent. Scope: RESULT_Q_STAGE0/12/3/4/5/6/7, BELOW_FOLD/2/3, PARITY_ODD_CLOSE, HYPOTHESES_CLOSE2
(context: RESULT_D_STATEMENT top box, EE_FUNCTIONAL). Result D's high-spin below-fold sector is out of scope. No existing file
was edited; nothing was committed; no `bj_*` file was touched. The spot checks below were run in memory and wrote no files.*

Grades as in the project: [P] [P given X] [D] [D-sketch] [Q] [DV] [NE] [M] [I].

---

## 0. Verdict

**Major revision.** Result Q can be kept only as a **conditional statement**, with a much longer hypothesis list (§6) and the
regradings in §5. I found no error that kills the headline upper bound at `ĝ₃ = 0`. But:

* "heavy positivity over the full set, with measured graviton corrections" is overclaimed;
* the `|ĝ₃|²`-dependence (`¼`, `κ`) is graded far too high;
* several "exact" light-side statements rest on an assumption, (P), that does not cover the objects it is applied to;
* the `Δ_g = 100` column has no support.

The numbers reproduce (§1). The chain of logic has the gaps listed in §2–§4.

---

## 1. Spot checks (re-run, in memory only)

| check | claimed | re-computed |
|---|---|---|
| `2L_AdS = 2∫ψ₂/(p²+9/(4Λ²))`, δ = 0.2 functional, `Λ = √((Δ_g−2)(Δ_g−1))` | 53.11 / 72.04 / 92.81 | 53.107 / 72.041 / 92.814 |
| `B_eff = 2L_AdS − A log(2Λ/3)`, `A = 2ψ₂′(0) = 17.2407` | −19.038 / −19.220 / −19.265 | −19.0378 / −19.2197 / −19.2647 |
| ABJM `2ζ(3)/2L_AdS` | 4.5 / 3.3 / 2.6 % | 4.53 / 3.34 / 2.59 % |
| `Λ = 2^{5/4}π^{1/2}λ^{1/4}`; `g₄M⁶/(8πG) = ζ(3)α′³/32 · 64/α′³` | 4.2156; 2ζ(3) | 4.21563; 2ζ(3) ✓ (`Λ⁶`, not `Λ⁸`, is the right power: `g₄/G ~ L⁶`) |
| Lemma B: `f⁺₄ = f⁺₂R`, `f⁻₅ = f⁻₃R` | [P] | sympy: both differences ≡ 0 ✓ |
| `++` J = 0 gap state, flat parts of the δ = 0.2 functional | FF/N = 0.200 | `B₂ = 1.805`, `B₃ = −0.203`, forward pair `a(2𝒥 − c₄ + σ) = 0`, `cB₄ = 1`, `N = 3.008`, FF/N = 0.2000 ✓ |
| effect of the measured `s_d = 0.509` there | 0.200 → 0.1427 | `2a s_d/N = 0.1126 × 0.509 = 0.0573` → 0.1427 ✓ |

So the constants and the bookkeeping are right. The binding state's positivity is carried entirely by the **smeared** parts:
the forward pair is identically zero at `++` J = 0, and `B₂` is 60 % of `N`. This matters for §3.1.

---

## 2. Light side: Lemmas A/B/C, V₄, (F), (P)

### 2.1 Lemma A (Einstein, `E×W³`, MFT annihilated by `Ψ₄`, `Ψ₅`) — correct in structure, over-graded
* The commutator argument is the standard one (CHMRS: "graviton exchange contributes only to `k = 2`"; CLS-D = Chang–Landau–Simmons-Duffin,
  arXiv:2311.04271). I see **no surviving `O(G)` term**:
  * the forward `G`-pole lives only in `k = 2`;
  * the identity and `T` terms are inside `Ψ[1 + T + [TT]]`, which Lemma A sets to zero as a whole.

  So the "light side exactly zero" claim is plausible.
* **But the grade [D] is not earned. It should be [D-sketch].**
  * The vanishing is proved in the **commutator** form. The sum rule is used in the **swapped (s-channel)** form. Their equality, with the
    `1/(u′v′)` poles, is V₄, and V₄ has **no independent test**.
  * The MFT "direct test" is vacuous. The light-ray-polarised `⟨T₃T₁⟩` numerator vanishes for **any** `f`, so the test cannot see the
    subtraction poles (Stage 3 §1.6 admits this).
  * The swapped MFT test and the swapped Einstein test were not done.
* **Factual error in Lemma C / §1.4 (Stage 3).** The text says the edge produces "only the finitely many light-family terms (identity, T,
  leading `[TT]`)". The leading family `[TT]_{0,J}` is **infinite**:
  * it has anomalous dimensions `γ_{0,J} ~ 1/J` from `T` exchange in `d = 3`;
  * with only a single zero left by `f`, its whole large-`J` tail contributes.
  * Per-block Gel'fand–Shilov continuation followed by summation over `J` need not equal the continuation of the sum. That tail sums to
    the crossed-channel light-cone singularity, which is exactly where the poles of `f` sit.

  This is the unverified step of V₄. It is also the place where an `O(G)` remnant would hide, with the `10¹⁰–10²²` leverage of the
  Stage 12 table.
* **Required:** one swapped-form evaluation of `Ψ₄` on MFT, or on Einstein `⟨TTTT⟩`, at one `ν`. Stage 12 already listed it as the
  "cheap independent check"; it was never done.

### 2.2 "`−Ĉ₄^{ds}|_{R⁴} = 2g₄` exactly", "`−Ĉ₅^{ds}|_{g₅} = g₅` exactly" — not [P given (F),(P)]
* (P), the polynomial light side, is invoked for the "pure-moment family". That family exists only as a **leading-Regge** construct:
  `raw Ψ_{4,ν} = Σ_j h_j(ν) × pure(ν − ij)` is derived from `R → 16η²r²`.
  * The actual functional is `∫χ(ν) Ψ_{4,ν}^{raw} dν`.
  * Its light side on a contact is `∫χ L_raw(ν)`. `L_raw` is never computed.
  * Contact light sides are not Regge-dominated, so the exact `R = 16η²r²(1 + O(r))` corrections enter.
  * Stage 4 §2.4 even concedes that the raw `Ĉ₄(ν)` "does not have a polynomial light side".
* Consequences:
  * the coefficient of `g₄` in the bound (hence `g̃₄ = g₄` and the ABJM evaluation) is [D-sketch];
  * the `g₅` cancellation between `∂Ĉ₄^{ds}` and the complex-gauge `Ĉ₅^{ds}` (two different de-shifts, one with an
    infinite-dimensional null space) is [D-sketch].
* Damage if the `g₅` cancellation is imperfect, by `ε g₅`: the shift is `≲ a·ε·g₄ = 0.17ε g₄`, using the flat `|g₅|M² ≤ g₄`, which
  itself needs its own AdS functionals. That is not fatal, but it is not "exact".

### 2.3 The `|ĝ₃|²` coefficients (`4πG|ĝ₃|²` "exact", `κ = 0.110`, lower slope `¼`) — **[I], not [P given (F),(P)]**
* `W³×W³` is a graviton **exchange** between two `W³` vertices, with a t-channel pole. It is not a contact, so (P) does not apply.
  Stage 12's own table lists it with `J₀ = 4`, i.e. not annihilated, and with a pole.
* Its AdS light side is `G|ĝ₃|²F(ν)` with `F(ν) = 4πν² + b₀ + (decaying residues, e.g. c/(ν²+9/4)ⁿ)`, the same residue structure as
  CHMRS's Lemma R. (F) fixes only the `ν²` growth.
* In the **derivative** piece, which reads `Λ²∂_{ν²}` at `ν = 0`, a decaying residue contributes `O(1)·G|ĝ₃|²`, **not** suppressed by
  `Λ⁻²`. The flat bookkeeping (`G_AdS|ĝ₃|²_AdS = G_f|ĝ₃|²_f/Λ¹⁰`) gives exactly that.
* The Stage 12 leverage table includes `b₀`, which is `Λ⁻²`-suppressed in the value piece, but **omits this term**.
* So:
  * `κ` and the lower-bound slope `¼` are uncertain at `O(1)` relative;
  * the lower bound "`¼|ĝ₃|²Λ⁸ ≤ g₄Λ⁶/8πG`" in AdS is [I];
  * the bound at `ĝ₃ = 0` is unaffected.
* The project computed exactly this boundary-layer effect on the **heavy** side, where it is not small. The measured `s_d ≈ 0.51` is
  `Δ_g`-independent (0.509 at 300, 0.508 at 1000), i.e. an `O(1/m²)` correction at `ν = O(1)` that survives the flat limit in the
  derivative piece. Assuming the light-side analogue is zero is not consistent with that.
* **Unremarked consequence.** Taken at face value, the two-sided statement implies a closure
  `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ+¼) ≈ 150 / 200 / 260` (Δ_g = 100 / 300 / 1000).
  * That is 30–50× stronger than Result D's headline (`≈ 7.5–14·10³`).
  * Either Result D's headline is superseded, which no document says, or Q's `|ĝ₃|²` sector is not trustworthy at that level.
  * I hold the latter: until the `W³W³` residue in `∂Ĉ₄` is computed, no `|ĝ₃|` statement should be drawn from Q.
* Minor: Stage 3 §2.3 writes `[1/(4c_L) − b₀/(16π)]`. That drops the `Λ⁻²` the Stage 12 leverage numbers use.

### 2.4 Contact tower and `c₁t₄`
* **Hidden assumption: the contact tower is controlled.** The light side is "Einstein + contacts". Every contact of degree `n ≥ 4`
  enters the de-shifted pieces through its `ν = 0` curvature constant. The statement "higher contacts `O(Λ^{−2j})`" is the CHMRS (5.3)
  assumption: infinite towers cannot conspire to give contributions enhanced at large `Δ_gap`. It must be stated as a hypothesis. It is
  not on any list.
* **`c₁t₄`.** The `k = 2` piece carries `16πG c₁ t₄` with `c₁` free.
  * Controlling it needs **Result D's `t₄` bound**, with all of D's per-state hypotheses (`H_M`, `H_R`, `H_LJ`, `H_OO`, `H_edge`).
  * Leverage is 0.35–2.9 per unit `c₁` at `Λ = 100`.
  * Result D at `Δ_g = 100` is "not established", so the `Δ_g = 100` column of Q has **no control** of this term.
  * Either carry `+c₁t₄` explicitly or drop `Δ_g = 100`.

---

## 3. Heavy side

### 3.1 Graviton corrections of the smeared parts were never computed (serious)
* Every "measured graviton correction" (BELOW_FOLD2/3, PARITY_ODD_CLOSE) is a projection `s_d` onto the **forward derivative** weight
  `χ_d` of the `k = 4` image. `bg_ffscan` / `bh_ffscan` / `po_ffscan` apply only `FF → FF − a·2u³s_d`.
* The smeared `∫ψ₂Ĉ′_EE` and `∫ψ₃Ĉ₃` parts are modelled throughout by "flat-exact + leading-Regge curvature" (HYB), with the **scalar**
  skeleton as the only beyond-LR check (Stage 3: ±0.005; the `k = 3` remainder is not computed at all).
* Ω^Q is a **different functional** from Ω^EE:
  * its ψ₂ is supported to `ν = Λ`;
  * it includes the fold region, the `ν < 8.3` windows and the `x → 1` edge window.
* Result D needed `H_M`, `H_LJ`, `H_edge`, `H_OO` and below-fold windows precisely for the graviton-specific corrections of these
  structures. None of that work carries over to Ω^Q, and none was redone.
* At the binding `++` J = 0 gap state, `B₂` is 60 % of `N`. A relative graviton correction `ε` of the smeared `k = 2` image shifts
  FF/N by `≈ 0.6ε`, so the δ = 0.2 margin (0.143) survives `|ε| ≲ 24 %`.
* By D's typical sizes (`|Γ| ~ 10²–10³ c₂`) this is likely safe at Δ_g ≥ 300 (`ε ~ 0.1–1 %`). It is not safe at Δ_g = 100
  (`ε ~ 1–10 %`), and it is unknown at near-gap moderate J.
* **Regrade** "min FF/N = +0.143 over the full set, with measured graviton corrections" to:
  * [M] (scalar skeleton + LR) for the smeared parts, whose graviton corrections are [I];
  * [NE given the engine equivalence] for the forward-derivative piece, measured in 5 `(h, J)` channels at one or two masses each, at
    Δ_g = 300 and 1000 only.

### 3.2 Coverage of the physical region
* Covered at the model level:
  * all J;
  * masses from the gap to `10⁴Λ`;
  * `++`/`+−`;
  * parity-odd even-J taken equal to `++` (PARITY_ODD_CLOSE lemma, order 0 [P given (A)]; beyond order 0 [NE], and only for the `k = 2, 3`
    functionals).
* Not covered:
  * **OO at `J = 0, 2`** for the `k = 4` de-shifted derivative piece. This is the **binding** channel, with the same margin as `++` J = 0.
    It is assumed equal to `++` by extending a lemma checked only on Ĉ′_EE and Ĉ₃ at J ≥ 4. Lemma B makes `R` polarisation-blind, so it is
    plausible, but it is [D-sketch].
  * odd-J `+−` beyond J = 5 [M];
  * the `s_d` m-dependence within a channel [M];
  * β > 30;
  * the 2×2 helicity mixing beyond the measured `m²κ ≈ 1.3`.
* **The realised functional is not the positive functional.**
  * Stored double weights fail from β = 15.6.
  * The XP C4 weights hold to t ≈ 30.
  * C5 is realised only to t ≈ 18 (no spectral gap). Stage 7 extrapolates C5 failure at β ≈ 31–32.
  * Positivity at larger J holds only for the **continuum** Abel solution. Its existence and uniqueness of image are not proved for the
    complex 3-shift (infinite-dimensional null space): [D-sketch].
* The flat LP grid stopped at J = 2500. The out-of-grid flat dips (+0.185 at δ = 0.2) set the Δ_g → ∞ floor. They are positive. Fine,
  but the "LP certified" language should not be used.

### 3.3 Double counting (BELOW_FOLD3 §4.2)
* The measured correction is `σ_h U(ν)`: universal, helicity-odd and **Δ_g-independent**. Per order, `s_d` = `+1.66 −3.95 +2.52 +0.27`
  at Δ_g = 300 and `+1.66 −4.21 +2.96 +0.09` at 1000. The `a⁰` and `a¹` terms are identical across Δ_g and come from the Δ-dependence of
  the conserved `⟨TTO⟩` coefficients (BELOW_FOLD2 §3.2).
* A correction that survives `Δ_g → ∞` must be a flat-space effect. The flat normal form HYB uses, `(c₄_h − 2𝒥)u³` with a helicity-odd
  constant ∓20, is exactly where such a term lives.
* So one of two things holds:
  * (a) it is double-counted against HYB's `γ_h`; or
  * (b) the flat image of the de-shifted `∂Ĉ₄` at `p = 0` is **not** CHLPS's `∂_{p²}B₄|₀`. The boundary layer at `ν = O(1)` is invisible
    to the pointwise flat limit (F).
* Under (b), the light side must carry the matching boundary layer, which reinforces §2.3.
* Either way δ = 0.2 survives:
  * under (a), `++` J = 0 returns to ≈ 0.198 and `+−` to its K = 0 value ≈ 0.20;
  * under (b), the quoted +0.143 stands.
* δ = 0.05 is undecided. The documents already say so; the point is that (b) is a statement against "flat image [P]" for the derivative
  piece and should be recorded as such.

### 3.4 Real section ↔ crystal cycle for EE beyond LO
* Assumed, not shown. The real-section QMC itself **fails** for J > 0 beyond ν ≈ 2.6 (boost growth, far region).
* The analytic-radial engine silently replaces divergent radial and angular integrals by Γ/Beta continuations. For EE these include the
  `cos^{−4}φ` light-cone enhancement and the `(Y₁/Y₂)^{1.7}` growth.
* Validated:
  * scalar, all J ≤ 32 [DV];
  * EE at a⁰ (zeros of `N̂₄`);
  * EE `k = 2` against the K₂ chain at ν = 20–40 to 2 c₂ (BELOW_FOLD3 §2), a genuine beyond-LO EE check, but at large ν and for `k = 2`;
  * one odd-J n3 comparison at 1.5σ.
* No EE check exists at ν < 10, where ≈ 96 % of `χ_d`'s weight lies.
* Grade: [NE given an assumed equivalence], as the documents say. Do not upgrade it.

### 3.5 Validity of the de-shifted functionals
* The counting "`J₀ < 7/2`" is done on the **heavy-normalised total kernel**. Stage 4 §2.3 concedes the position-space `ν`-integral is
  not absolutely convergent (`1/E_b ~ e^{π|ν|}`). With `ν_max = 123` (XP), the position kernel is dominated by truncation tails.
* So the de-shifted functional is defined in heavy-image space, while each raw `Ψ_{4,ν}` fails KKSDZ (needs `J₀ < ½`). The swap argument is
  as circular as D's I14(b).
* **Regrade:** validity [D] → [D-sketch] (inherited).

---

## 4. Normalisations, ABJM, novelty

* **R⁴ normalisation:** `−Ĉ₄(0)|_{R⁴} = 2g₄` follows from the flat `B₄|₀` low side, given (F) and (P). Both are D-sketch for spinning
  functionals. The ratio to the `16πG` pole is convention-free only if both are heavy-normalised in the same way; that is asserted, not
  checked.
* **Dictionary:**
  * The BCP/Penedones chain reproduces `B⁴₄/B¹_SG = 105ζ(3)/(64√2π³)` [DV].
  * `Λ = 2^{5/4}π^{1/2}λ^{1/4}` ✓; `2ζ(3)` ✓; `Λ⁶` ✓.
  * `c_T ↔ 8πG = 256/c_T` is quoted from ADS4_LIGHT [Q + D]. Flavour-current coefficients play no role, as expected for `⟨TTTT⟩`.
  * `b₀`, `b₅`, `c₁` are uncomputed. `b₀` and `b₅` are harmless (Λ⁻²) in the value piece. `c₁` is not harmless at Δ_g = 100 (§2.4).
  * The `W³W³` derivative residue is missing from the budget (§2.3).
* **ABJM placement:**
  * Grade it [D-sketch given Q]: it needs (F) + (P) spinning, and the identification `g̃₄ = g₄` (§2.2).
  * ABJM is outside A1 by the project's own statement (light KK modes), so this is a **consistency check of a flat image, not a test of
    the bound**.
  * The placement at 2.6–4.5 % stands numerically, but it is log-dominated and essentially uninformative: any theory with
    `Λ`-independent `g₄M⁶/G` of `O(1)` passes.
  * Do not present it as evidence for Result Q.
* **Novelty:**
  * I searched for a published AdS₄ bound on the `⟨TTTT⟩` R⁴ contact from dispersive CFT sum rules and did not find one. CLS-D
    (2311.04271) constructs the spinning sum rules and the crystal dictionary but derives no R³/R⁴ bound ("remaining steps for uplifting"
    are left open). CHLPS 2201.06602 only *interprets* its flat bounds in AdS₄.
  * The novelty claim (first explicit-constant dispersive CFT₃ `⟨TTTT⟩` R⁴ bound, conditional) therefore **plausibly stands**.
  * But the literature table must add the **closest real precedent**: CHMRS 2106.10274. It gives two-sided, `log(Δ_gap)`-type AdS₄ bounds
    on scalar contact couplings with dynamical gravity, through the same subtracted `C_{k,ν}`. The table cites CHMRS only as a source of
    exact moments.
  * Kundu (2104.11238, gap-scaling of higher-derivative couplings) and the N = 6 numerical bootstrap (Binder–Chester–Jerdee–Pufu
    2011.05728) should be discussed, not merely listed.
  * My search was abstract-level, as was the project's. Grade [I].

---

## 5. Required regradings

| claim (where) | as graded | should be |
|---|---|---|
| Lemma A (Stage 12 §2.1, Stage 3 §1.7, V₄) | [D] | [D-sketch]; "finitely many light-family terms" is false (infinite `[TT]_{0,J}` tail) |
| MFT test "exactly 0" (Stage 3 §1.6) | [P]+[DV] | correct but vacuous for `f`; not evidence for V₄ |
| `−Ĉ₄^{ds}|_{R⁴} = 2g₄`, `−Ĉ₅^{ds}|_{g₅} = g₅` (Stages 4, 5) | [P given (F),(P)] | [D-sketch]: (P) is applied to an undefined pure-moment family; the raw light side is uncomputed |
| `−Λ²∂Ĉ₄(0)|_{W³W³} = 4πG|ĝ₃|²` exact; `κ = 0.110`; lower slope `¼` in AdS | [P given (F),(P)] / [M] | [I]: exchange, not contact; `O(1)` residues in the derivative piece; missing from the leverage table |
| full-set min FF/N = +0.143 "with measured graviton corrections" | [NE + M] | [M] (smeared parts: graviton corrections [I]) + [NE given D ≡ crystal] (forward piece, 5 channels, Δ_g = 300/1000) |
| OO (parity-odd) at J = 0, 2 for the `k = 4` piece | implicit = `++` | [D-sketch]; list as H_OO^Q |
| validity of de-shifted `Ĉ₄`, `Ĉ₅` (`J₀ < 7/2`) | [D] | [D-sketch] (heavy-normalised kernel only; position kernel not absolutely convergent; continuum C5 solution unproved) |
| `B_eff` at Δ_g = 100 in the bottom line | quoted | unsupported (no measured correction at 100; Result D's `t₄` bound absent at 100, so `c₁t₄` is uncontrolled) |
| ABJM 2.6–4.5 % | [D given Q] | [D-sketch given Q], outside A1: a consistency check, not a test |
| "`¼` lower bound survives [M]" | [M] | not re-evaluated with the measured `s_d`; depends on §2.3: [I] |
| HYB "covers every physical state" (Stage 6) | [M] | [M] for the scalar skeleton + LR only; β > 30 and the realised-weight limits (t ≈ 18 for C5) excluded |

---

## 6. The hypothesis list Result Q actually needs

1. **A1**: twist gap Δ_g for all heavy states in T×T, all J (excludes ABJM).
2. **A2**: tree level, i.e. O(1/c_T) heavy sum.
3. **EFT light sector**: the `[TT]` data equal Einstein + `E×W³` + `W³W³` + local contacts. The contact tower cannot conspire at large
   Δ_g (the CHMRS (5.3)-type assumption).
4. **Regge/validity**: KKSDZ `J₀ ≤ 1` with its loophole; K1′, K2, K3 for the de-shifted `k = 4, 5` functionals, swapped via the
   heavy-normalised kernel.
5. **V₄**: swapped-form edge/endpoint bookkeeping with the `1/(u′v′)` poles, including the infinite leading `[TT]₀` tail. Also
   microcausality of tree correlators and `J₀(E×W³) ≤ 3`.
6. **(F) + (P) spinning**: for the `g₄` coefficient, and for the `g₅` cancellation between two different de-shifts.
7. **`W³W³` derivative-piece residue** (new). Needed only for the `|ĝ₃|²` terms (`κ`, `¼`, the implied ĝ₃ closure).
8. **Existence of the continuum de-shift solutions** (C4, C5) with the claimed images for all β. The realised weights hold only to
   t ≈ 30 / 18.
9. **Result D's `k = 2` inputs**: H7 (exact pole), EE admissibility, `c₁`. Plus either D's `t₄` bound, with all of D's per-state
   hypotheses, or an explicit `c₁t₄` term in Q.
10. **H_smear^Q**: graviton-specific (non-skeleton) corrections of `∫ψ₂Ĉ′_EE` and `∫ψ₃Ĉ₃` in Ω^Q stay below the remaining margin
    (≈ 24 % relative of `B₂` at `++` J = 0, gap). Never computed.
11. **H_fwd^Q**: the forward-derivative corrections at the unmeasured (h, J, m) equal the measured universal `σ_h U` or are smaller;
    m-independence within a channel. Includes the real-section ≡ crystal equivalence for EE beyond LO.
12. **H_OO^Q**: the parity-odd even-J channels (J = 0, 2 binding) equal `++` for the `k = 4` piece.
13. **HYB accuracy** (scalar skeleton + LR) away from the measured states; β > 30, `m > 10⁴Λ`.
14. *(δ = 0.05 only)* resolution of the `σ_h U` double-counting question.
15. *(`|ĝ₄| ≤ g₄` only, not part of Q's headline)* the all-plus non-shock rule `B₄⁽⁴⁾`.

**Honest headline:** given 1–6 and 8–13, at Δ_g = 300 and 1000,
`0 ≤ g̃₄Λ⁶/(8πG) ≤ 17.24 log(2Λ/3) − 19.22 / −19.27 (+ c₁t₄ term)`, at `ĝ₃ = 0`. Here `g̃₄` is the CFT-defined light-side coefficient,
equal to the bulk `g₄` given (F) + (P). The `|ĝ₃|²` dependence (`¼`, `κ`) requires item 7 as well.

---

## 7. Cheap next steps that would most change the grades

1. A swapped-form `Ψ₄` on MFT `⟨TTTT⟩` at one ν (tests V₄ and the `[TT]₀` tail).
2. The `W³W³` residue in `∂Ĉ₄(0)`: one AdS exchange with W³ vertices in the Regge/Lemma-R form.
3. The smeared `k = 2` EE correction at `++` J = 0, 2 at the gap. The analytic-radial engine already produces the EE `k = 2` images
   (`bg_eeq`). Project them on ψ₂ at ν ≲ 40, and bound the rest by D's sizes.
4. Drop Δ_g = 100, or state the `c₁t₄` dependence explicitly.
5. Record the implied `|ĝ₃|` closure, or explicitly disclaim it.

Sources consulted for §4: [CLS-D 2311.04271](https://arxiv.org/abs/2311.04271);
[CHMRS 2106.10274](https://arxiv.org/pdf/2106.10274); [CHLPS 2201.06602](https://web3.arxiv.org/pdf/2201.06602);
[Knop–Mazáč 2203.11170](https://arxiv.org/pdf/2203.11170); [CCTZ 2412.00249](https://arxiv.org/pdf/2412.00249).
