# R³ and one R⁴ coupling of holographic CFT₃ from ⟨TTTT⟩ light-ray sum rules: the conditional theorem, consolidated (v2)

**Ansh Mishra and Aryan Senthilkumar**

*23 Nov 2026. Writing pass: no new computation. It replaces `UNIFIED_THEOREM.md` (UT) as the source of truth for the letter. UT's
body and its fourteen Update boxes were read in full; where they conflict, the later box wins (the boxes were appended in time order:
HRS_CLOSE → EINSTEIN_CHECK → DG_SCAN → REPRO → REPRO_ENGINE → RD_CLOSE → HLJ_CLOSE → C1_COMPUTE → LIGHTSIDE_DIST → HK_CLOSE → F_NORM →
LEMMAS_CLOSE → LRC_CLOSE → DESHIFT_CLOSE). Every number below carries a source tag, resolved in §12. Unless marked "independent",
every source document is unreviewed.*

**One item is still being computed elsewhere:** the EE (spin-2) step of LRC\* (`EESTEP_CLOSE.md`, not yet written). It is marked
**⟦EESTEP⟧** wherever it enters. Until it lands, H_K is graded [D-sketch] through that step.

---

## 0. Grades and how to read this document

| grade | meaning |
|---|---|
| [P] | proved |
| [certified] | computer-assisted proof (interval arithmetic) |
| [D] | derived (argument written out, not independently checked) |
| [D given X] | derived from the named inputs X; its net grade is the weakest of [D] and the grades of X |
| [D-sketch] | derivation sketched; at least one step is asserted from a standard result without being written out or checked |
| [DV], [D + DV], [P + DV] | derived (or proved) and checked numerically/symbolically |
| [Q] | quoted from the literature |
| [NE] | numerical evidence |
| [NE + M] | computed inside a model (the near-gap chain / hybrid model, §5); the model's accuracy is itself a hypothesis |
| [NE + D-sketch] | numerical evidence whose interpretation needs a sketched step |
| [M] | model: an approximation whose accuracy is a hypothesis |
| [A] | assumption |
| [I] | open |

**Weakest-link rule.** A chain carries the grade of its weakest link, ordered [P] = [certified] > [DV] > [D] > [D-sketch] > [NE] >
[M] > [A] > [I].

**Three layers, kept apart throughout.**
1. **Setting** (A1, A2): what theories the statement is about.
2. **Inputs of the implication** (§4): analytic statements about the functionals and the light side. These are what "positivity ⇒
   bound" needs.
3. **Per-state positivity** (§5): the one hypothesis about the heavy spectrum. The theorem *assumes* it; §5 lists the *evidence* for
   it, which is numerical and model-based on finite sets. Evidence is never an input.

---

## 1. Setting, normalisations, couplings

| item | definition | source |
|---|---|---|
| theory `𝒯` | a unitary CFT₃, parity-even or parity-violating, with a unique stress tensor `T` | UT §1 |
| `c_T` | normalised so that a free real scalar has `c_T = 1` | UT §1 |
| Newton constant | `8πG := 256/c_T` (AdS units `L = 1`), from `c_T = 32L²/(πG)`; independently confirmed | UT §1; REPRO §1 |
| **A1** (gap) [A] | every single-trace primary other than `T` that appears in the `T × T` OPE has twist `Δ − J ≥ Δ_g`, for all `J` including `J = 0, 1, 2`. (At `O(1/c_T)` only single traces in `T × T` and the `[TT]` double traces enter `⟨TTTT⟩`, so this is all the sum rules see. Every scan used the stronger "all single traces".) A1 excludes every known holographic CFT₃; in particular ABJM has light Kaluza–Klein and R-singlet single traces. | UT §1; RD-S §1; REFEREE M12 |
| **A2** (tree level) [A] | `O(1/c_T)` truncation with `c_T ≫ Δ_g²`; loops dropped. All bounds are statements about the `c_T → ∞` limit at fixed `Δ_g`. | UT §1 |
| gap scale | `Λ := √((Δ_g − 2)(Δ_g − 1)) = Δ_g − 3/2 − 1/(8Δ_g) + O(Δ_g⁻²)`: `298.4996` (300), `998.4999` (1000) | UT §1; DG §0 |
| heavy kinematics | `m² = (Δ − J − 2)(Δ + J − 1)`; gap state `m_gap² = (Δ_g − 2)(Δ_g − 1 + 2J)`; `cosh β = 1 + (2J+1)²/(2m²)` | UT §1; RD-S §1 |
| cubic couplings | `⟨TTT⟩ = (c_T, t₄, α_T)`; `t₄` is the d = 3 Hofman–Maldacena parameter, `α_T` the coefficient of the parity-odd structure (exists only in d = 3). `t₄`, `α_T` are dimensionless CFT data. Bulk dictionary `√(t₄² + α_T²) = 240|ĝ₃|` with `ĝ₃` the `R³` coupling in CHLPS normalisation in AdS units. | UT §1; C1 §1.1; DICTIONARY_T4 |
| **quartic coupling `g₄`** | defined by the flat four-graviton amplitude `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u)`, `f = 8πG/(stu) + … + g₄ + g₅t + …` (CHLPS (2.1), (2.7a)). Normalisation: type II gives `g₄M⁶/8πG = 2ζ(3)`, `M² = 4/α′` (CHMRS 2102.08951 §3.8). | F_NORM §1 |
| quartic structures **not** bounded | In D = 4 there are three real eight-derivative four-graviton couplings: `g₄` (MHV) and the complex all-plus coefficient `ĝ₄` (two real). The k = 4 functional assigns `ĝ₄` the value 0 [D], so `ĝ₄` is not bounded. `g₅` (ten derivatives) is not bounded either. | F_NORM §0, §2.4 |
| `g̃₄` (CFT-defined) | `g̃₄ := −½ 𝒜₄[light] = −½ Ĉ₄^{ds}|_light` (= minus the heavy sum, given H_K) | F_NORM §1; LIGHTSIDE §0 |
| `c₁` | the coefficient of the `O(t₄)` light term of `Ĉ′_EE` (and `Ĉ₃`): `16πG c₁t₄`. One universal pure number, uncomputed (§6.5). | C1 §0 |

---

## 2. The functionals (cut-off, defined by the published weights)

All functionals are built from CLS-D's per-ν light-ray functional `Ψ_{k,ν,j′}` on `⟨TTTT⟩` (2311.04271, eq.
`spinningfunctional_definition`), in heavy normalisation `Ψ̂_ν = Ψ_{k,ν}/N_com(ν)`. The data files are in `ads4/prl/data/`
(README, SHA-256 prefixes in `MANIFEST.txt`).

### 2.1 The analytic cutoff
`θ(ν) := ½[erf((ν + ν_c)/a) − erf((ν − ν_c)/a)]` (even, entire, bounded on strips, super-Gaussian decay). Every de-shifted piece is
defined with weight `χ_θ := χ·θ`:
* Result Q pieces: `ν_c = 150`, `a = 6`;
* Result D pieces: `ν_c = Λ + 100`, `a = 10`.

The cutoff changes no computed number: relative effect `≲ e^{−180}` (Q), `≲ 10⁻⁴¹` (D) [D + NE] (HK §0, §5). It is required: any
weight of compact ν-support (every grid weight) regrows at the raw Regge rate by Paley–Wiener, so the heavy-sum interchange would need
`J₀ < ½` [P] (HK §3). (Notation: UT/HK call the cutoff `c(ν)`; renamed here because `c(ν)` is already Result D's EE kernel factor.)

### 2.2 `Ω^EE_θ` (Result D, item 1)
`Ω^EE_θ := ∫₀¹dp ψ₂(p) Λ⁴ Ĉ′_{EE,Λp} + ∫₀¹dp ψ₃(p) ½Λ⁶ Ĉ_{3,Λp}`, with:
* **weights** NF02: `ψ_k(p) = (1−p)⁵ p X_k(p)/k₀`, exact rational coefficients in `nf02_psi.json`, `k₀ = 27.48509817`, normalisation
  `k(0) := −2∫ψ₂p⁶ + 2∫ψ₃p⁴ = 1`;
* **k = 2 part:** `j′ = 0`, unsubtracted `f = v′ − u′`, dual (2,4) structure
  `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²`, de-shifted: the smeared form is `Re∫dν X_θ(ν)Ψ̂_{EE,ν}`, `X_θ = X·θ`, where
  `X` is the **min-norm** solution of the generalised Abel equation `C[aY] − i tanh t·S[bY] = C[ψ₂]/cosh t`, `X = aY` (kernels
  `R_∓`, `a`, `b`, `n` in `kernels.json`; 2-dimensional null space; min-norm = `Y − ψ₂` of minimal `L²(dν)` norm);
* **k = 3 part:** `Ĉ₃ := Re[Ψ_{3,ν}[o_H·q]/N̂₃(ν)]`, `j′ = 1`, unsubtracted, no de-shift.

**Status of the definition.** Exact existence of the min-norm `X` is [D + NE] (Abel operator bounded and solvable; residual
`≤ 3·10⁻¹⁵`). The published `d_abel_X_{300,1000}.txt` tabulate `X` on the grids `dν = 0.02` (to 327.68) / `0.04` (to 1310.72); they
were reproduced by an independent solver to `≤ 1.5·10⁻⁶` (ν < 2) and `≤ 3·10⁻¹¹` (ν > 25) of max ψ₂ [DV] (DESHIFT §4.1). Unlike Q's
weights, **no continuum k = 2 weight has been built** (DESHIFT §7 item 2): `Ω^EE_θ` is defined by the exact `X` times `θ`, and the
published grid is its numerical realisation. The stored double-precision grid satisfies the de-shift identity only to `t ≈ 22–25`
(HK §0, DESHIFT §7); the analytic control beyond is HK's Lemma S2 [D].

### 2.3 `Ω^Q_up` (Result Q, upper end)
The flat target is `Q3-full-delta0.2` (`q3_full_delta0.2.json`): on a flat heavy state
`FF = ∫ψ₂(2 − y)D_h + s u∫ψ₃D_h + a(2J(J+1) − c₄ + s)u³ − 2cu²`, `A = 2ψ₂′(0) = 17.240713`, `κ = 0.109908`, `a = 0.169405`,
`c = ½`, margin δ = 0.2 on the LP grid (dense-grid flat minimum 0.1966, i.e. δ ≈ 0.2 is nominal at the 2 % level; REPRO D9). It is
realised in AdS by:
* **C4 value / C4 derivative** (`Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`): continuum weights `q_deshift_continuum_c4v.json`, `…_c4d.json`,
  `χ = Q(ν)B(ν)` with the explicit rational `Q` and a 401-term Gaussian basis `B` (coefficients to 60 digits; evaluate in ≥ 300-bit
  arithmetic); de-shift identity `|W − Φ| ≤ ε_W` **for all t ≥ 0**, `ε_W = 1.2·10⁻¹³` (value), `4.2·10⁻¹¹` (derivative) [DV]; equal to
  the stored Stage-7 XP weights to `1.9·10⁻¹⁷` / `1.1·10⁻¹⁶` of max on their grid (DESHIFT §2.3);
* **C5** (`Ĉ₅^{ds}`): `q_deshift_continuum_c5.json`, `ε_W = 2.4·10⁻⁷` for all `t ≥ 0` (with the pole `a₅` defined as the exact root
  of the quintic given in the README); a different gauge from the stored XP weight, same leading-Regge image (DESHIFT §2.3);
* **smeared k = 2, 3 parts** `∫ψ₂^Q Ĉ′_EE` and `∫ψ₃^Q Ĉ₃` with Q's `ψ₂, ψ₃`: Result D's objects, the k = 2 part de-shifted by the
  min-norm Abel weight `X^Q` of Q's `ψ₂` (same equation as §2.2; LSQR residual `10⁻¹⁵`; RESULT_Q_REVISION §1.3). **`X^Q` is not among
  the published files**; it is fixed by the definition plus the published `ψ₂^Q` and kernels, but publishing it is still owed.

All continuum weights decay super-Gaussianly beyond ν ≈ 60, so they are already cut-off weights; multiplying by `θ` (150, 6) changes
them by `< e^{−300}` relative (DESHIFT §2.4).

### 2.4 `Ω^Q_low` (Result Q, lower end)
The de-shifted k = 4 value piece alone, i.e. `Ĉ₄^{ds}` with the C4-value continuum weight. Leading-Regge heavy image `P_{−½}(cosh β) > 0`;
flat image `2(|c⁺⁺|² + |c⁺⁻|²)/m⁸ ≥ 0`; light side `2g̃₄` (RQ-S §1).

### 2.5 What "per-state positivity" means
For a heavy primary `O = (Δ, J)` with `Δ − J ≥ Δ_g` (any parity), the action of a functional `Ω` on the conformal blocks of `O` in
`⟨TTTT⟩` is a Hermitian matrix `M_Ω(O)` in the helicity/parity basis of `⟨TTO⟩`: a 2×2 block on (`++`, `+−`) and, for parity-odd
even-`J` exchange, a 1×1 block (OO). By unitarity the OPE-coefficient matrix `P_O` is PSD, and the sum rule involves `Tr[P_O M_Ω(O)]`.

**(Pos_Ω)** `M_Ω(O) ⪰ 0` for every twist-allowed heavy primary `O`.

Numerical values are quoted as `λ/N_s` (Ω^EE) or `FF/N` (Ω^Q): the smallest eigenvalue of `M_Ω(O)` in units of the state's flat
normalisation, the unit in which the flat functionals have minimum `≈ 0.200`. Unit for beyond-Regge corrections:
`c₂ := m² × (relative correction)`. (RD-S §1, §5.1; RQ-S §3.)

---

## 3. The theorem

> **Theorem (conditional).** Let `𝒯` satisfy A1 and A2 at `Δ_g ∈ {300, 1000}`. Assume the inputs of §4 and (Pos_Ω) for
> `Ω^EE_θ` (item 1) and for `Ω^Q_up`, `Ω^Q_low` (item 2). Let `c₁` be the universal constant of §6.5. Write
> `R := √(t₄² + α_T²)`, `L_A := ∫₀¹ψ₂(p)dp/(p² + 9/(4Λ²))`, `K := k(0)(1 + δ₃)` with `δ₃ = δk(0) + (H17 remainder)`.
>
> **(1) Cubic.** For every real `c₁`,
>
>     R ≤ R₀ (x + √(1 + x²)),   R₀ = 240 √(16L_A/K) / Λ⁴,   x = 1920 I_D |c₁| / (Λ² √(K·16L_A)) = 19941.0 |c₁| / (Λ² √(K·16L_A)),
>
> with `I_D = ∫ψ₂ = 10.385947`. At `K = 1`: `R₀ = 3.13206·10⁻⁶` (`Δ_g = 300`) and `2.89126·10⁻⁸` (`Δ_g = 1000`);
> `x/|c₁| = 2.16006·10⁻³` / `1.67026·10⁻⁴`. The `t₄` interval is asymmetric: for `c₁ > 0`, `−R₀(√(1+x²) − x) ≤ t₄ ≤ R₀(√(1+x²) + x)`.
> Equivalently (at `c₁ = 0`) `|ĝ₃|²Λ⁸ ≤ 16L_A/K = 1.07·10⁴ / 1.43·10⁴`.
>
> **(2) Quartic.** Add, if `(t₄, α_T) ≠ 0`, the sign inputs `κ > 0` (upper) and `c_L > 0` (lower) of §4. Then
>
>     −O(Λ⁻²) ≤ g̃₄Λ⁶/8πG ≤ 2L_AdS − κ R²Λ⁸/57600 + 2 I_Q c₁ t₄ Λ²,      2L_AdS = 17.240713·ln(2Λ/3) + B_eff,
>
> with `I_Q = ∫ψ₂^Q = 1.0160366`, `B_eff = −19.2197` (300) / `−19.2647` (1000), i.e. `2L_AdS = 72.041 / 92.814`. Taking the supremum
> over the unknown `(t₄, α_T)`, for every real `c₁`:
>
>     g̃₄Λ⁶/8πG ≤ 2L_AdS + 57600 I_Q² c₁² / (κΛ⁴) = 72.041 + 6.8145·10⁻⁵ c₁²  (300),   92.814 + 5.4428·10⁻⁷ c₁²  (1000).
>
> If `t₄ = α_T = 0` the `c₁` and `κ` terms are absent and no sign input is needed. For `(t₄, α_T) ≠ 0` the lower end reads
> `g̃₄Λ⁶/8πG ≥ R²Λ⁸/(57600·4c_L) − O(Λ⁻²)`.
> **The lower end has no computed constant:** "`−O(Λ⁻²)`" means `≥ −O(Λ⁻²) × (EFT data)`, from the curvature constants of the light
> terms (`b₀`, `s_v`, the contact tower).
>
> **(2′) Bulk coupling.** Given H_EFT(b) [A] (§4, row L7), `g̃₄ = g₄ + O(Λ⁻²)·(EFT data)`, so (2) bounds the bulk Wilson coefficient `g₄`
> at the same order.

**Net grades (weakest link; §4 gives every link).**

| statement | net grade | weakest links |
|---|---|---|
| implication "(Pos) ⇒ (1)" | **[D-sketch]** | H_K through ⟦EESTEP⟧ (and LRC-a for EE at k = 4, Q only); `J₀(E×W³) = 2` behind the universal form of `c₁` [D-sketch]; the `±0.24 % / ±0.021 %` H17 band is [M] |
| implication "(Pos) ⇒ (2), upper end" | **[I] in one technical link; otherwise [A]** | a = b integrability of the borderline R⁴ arc [I]; H_EFT(b) [A] (only for (2′)); H_K [D-sketch]; (P) for `g₅` [D-sketch] (upper end survives either reading, §6.3); `κ > 0` [D-sketch + M] if `(t₄, α_T) ≠ 0` |
| implication "(Pos) ⇒ (2), lower end" | same, **and no constant** | as above; `c_L > 0` [D-sketch + M] if `(t₄, α_T) ≠ 0` |
| (Pos) itself | **not established** | computed [NE + M] on finite sets; one covering hypothesis [I] at 300 (§5) |
| bound constants (`R₀`, `2L_AdS`, `B_eff`, `x`, Λ-forms) | **[DV]**, reproduced independently | REPRO §1; DG §1 |

**Dependence.** Item (1) uses only `Ω^EE_θ`: it needs neither Lemma A, (F), (P), H_EFT nor any k = 4, 5 input. Item (2) needs nothing
from item (1): the exact `c₁` form uses only `κ > 0` (C1 §4.2), so the old route "bound `t₄` by item (1)" and its inheritance of
Result D's positivity hypotheses are gone. The old linear `|ε_{c₁}| ≤ 0.57|c₁| / 0.059|c₁|` is **superseded** (it dropped the
`−κR²` term).

**Sample values of the exact forms** (C1 §4, re-checked arithmetic here):

| `|c₁|` | 0 | 1 | 10 | 100 | 1000 |
|---|---|---|---|---|---|
| (1) at 300 | **3.1321·10⁻⁶** | 3.1388·10⁻⁶ | 3.2004·10⁻⁶ | 3.8808·10⁻⁶ | 1.4221·10⁻⁵ |
| (1) at 1000 | **2.8913·10⁻⁸** | 2.8917·10⁻⁸ | 2.8961·10⁻⁸ | 2.9400·10⁻⁸ | 3.4142·10⁻⁸ |
| (2) upper at 300 | **72.041** | 72.041 | 72.047 | 72.72 | 140.2 |
| (2) upper at 1000 | **92.814** | 92.814 | 92.814 | 92.819 | 93.36 |

The natural size of `c₁` is `10⁻³–10⁻¹` [M] (not a bound). The linearisation `R₀(1 + x)` is good to 1 % for `|c₁| ≤ 70.5` (300) and
`≤ 911` (1000).

---

## 4. Inputs of the implication (precise statements)

"D" = needed for item (1), "Q" = for item (2). Every row states a mathematical claim; the grade is the net grade of that row.

| # | statement | D | Q | grade | source |
|---|---|---|---|---|---|
| **R1** (Regge) | KKSDZ: `J₀ ≤ 1` for `⟨TTTT⟩` non-perturbatively (1904.05905 §4.2.4, eq. 4.55; generalises Caron-Huot 1703.00278), with the caveat of their §4.2.6 (first/second-sheet boundary) | ✓ | ✓ | [Q] | UT §3.1; LIT §1 |
| **R1′** (uniform Regge) | On the Regge box `B_c = {0 < ρ ≤ c, 0 < 1/ρ̄ ≤ c}`, `c < 1`, for all normalised polarisations: `|dDisc 𝒢| ≤ 2B(c)`, `B(c) < ∞` independent of η (edge exponent 0), from Rindler positivity + Cauchy–Schwarz on first-sheet comparison correlators. KKSDZ's loophole lies outside `B_c` (`ρ̄ ≥ 1/c > 1`) [P]. | ✓ | ✓ | [D]; scalar examples [DV] | LEMMAS §2 |
| **R2** (kernels) | (K2) endpoint/edge power counting of every per-ν kernel at every boundary [Q]/[D]; (K2ν) `|𝒦_ν(x)| ≤ K̄(x)` for all real ν (`K̄` = modulus kernel; `γ = 0` in position normalisation, π in heavy normalisation); `K̄` integrable against `|dDisc 𝒢|` away from the corner | ✓ | ✓ | K2 [Q]/[D]; K2ν [D], celestial factor [DV] | UT §3.1; LEMMAS §3.1 |
| **R3** (swap lemma) | (K3 / Lemma S) per-block Rindler positivity ⇒ `Σ_O∫|K||dDisc G_O| ≤ ∫|K| dDisc 𝒢` | ✓ | ✓ | [D] | HK §1; LRC §6.1 |
| **DS** (de-shift) | D: the min-norm Abel solution `X` of §2.2 exists, `χR_j ∈ L¹`, exponentially decaying, `O(ν²)` at 0 [D + NE]. Q: the published continuum weights satisfy `|W − Φ| ≤ ε_W` for all `t ≥ 0` (ε-DS) [DV]; exact C4 existence [NE + D-sketch] is not needed. | ✓ | ✓ | D [D + NE]; Q [DV] | HK §4.1; DESHIFT §0, §2 |
| **LRC-a** | leading form of the per-ν position kernel in the Regge corner: `r^a c_k(ν)ηⁿ P_{−½+iν}(η)(1 + O(r))`, `(a,n) = (0,1)` (k = 2), `(2,3)` (k = 4); heavy map `E_b(ν − ij)` | ✓ | ✓ | [D + DV] scalar k = 2, 4 and EE k = 2; **[D-sketch] EE k = 4** | LEMMAS §3.2 |
| **LRC\*** | the exact scalar `j′ = 0` per-ν kernel is `N(ν)r²Ĝ_{1,−½+iν}(ρ, ρ̄)`, a d = 3 continuous-spin Casimir eigenfunction (dimension 1) in the (13)(42) channel [D + DV, 3·10⁻⁸]; its light-cone-edge piece is the leading-Regge one times `F = ₂F₁(½, h; h+½; ρ²)`, `h = (½+iν)/2`, analytic for `0 ≤ Im ν < 3/2` and bounded for `ρ ≤ c` [P]. Hence `|K_θ| ≤ C η^{−½}log η` (k = 2) and `≤ C r²(ρ²η|f_pos(3i/2)| + η^{½})` (k = 4) on `B_c`. **EE structure:** the EE dressing is ν-independent and `z₀`-free [P + DV]; that the EE kernel equals the spinning partial wave with integer-shifted (24) labels (the **EE step**) is standard but not checked end to end. ⟦EESTEP: grade of the EE step from `EESTEP_CLOSE.md`, currently [D-sketch]⟧ | ✓ | ✓ | scalar [D + DV]; **EE [D-sketch] ⟦EESTEP⟧** | LRC §0, §1–§4 |
| **H_K** (interchange) | For each cut-off piece `Ĉ_θ ∈ {Ĉ′_EE, Ĉ₄^{ds}, ∂Ĉ₄^{ds}, Ĉ₅^{ds}}` and `Ĉ₃`: `Σ_O Tr[P_O Ĉ_θ[G_O]] = −Ĉ_θ[light]` with `Σ_O |Tr[P_O Ĉ_θ[G_O]]| < ∞`. Follows by Tonelli from R1′, R2, R3, DS, LRC-a, LRC\* whenever `J₀ < 3/2` (k = 2, 3, 5) and `J₀ < 2` (k = 4; `5/2` if `f_pos(±3i/2) = 0`). Under `J₀ ≤ 1`: margins **½** and **1**. No `J₀ < ½` and no light-cone hypothesis (LC\*, now graded [A], not used). | ✓ | ✓ | **[D-sketch]** (⟦EESTEP⟧; LRC-a EE k = 4 for Q) | HK §4; LEMMAS §3; LRC §6.1 |
| **L1** (D light side) | Adm: `EE` is an admissible `j′ = 0` dual structure [D]. L0: unsubtracted k = 2, 3 put every `[TT]` at a double zero and drop contacts [D]. H7: `−Ĉ′_{EE,ν}|_E = 16πG/(ν² + 9/4)`, `Ĉ₃|_E = 0` [D]; H7's value is the Regge arc of the same term-wise functional at `n = k = 2` [D] (LIGHTSIDE §5.2). | ✓ | (smeared part) | [D] | UT §3.1; LIGHTSIDE §0 |
| **L2** (`|ĝ₃|²` coefficient) | the `R²` light term of `Ω^EE_θ` is `K = k(0)(1 + δ₃)`: `δk(0) = −7.666·10⁻⁴ / −6.830·10⁻⁵` [D + DV]; the H17 remainder (`Π̂₄` admixture) is modelled: band `±212.5/Λ²` on the bound. Separate `t₄²`, `α_T²` coefficients at `O(Λ⁻²)` are not excluded; then (1) bounds a weighted norm. | ✓ | (κ) | δk(0) [D + DV]; **H17 [D-sketch + M]** | RD-S §3–4; DG §1 |
| **L3** (`c₁` form) | the `O(t₄)` light term of `Ĉ′_EE` and `Ĉ₃` is exactly `16πG c₁t₄`, ν- and `Δ_g`-independent, `c₁` one universal number; its `u`-channel T-exchange part and its k = 3 part vanish | ✓ | ✓ | form [D given `J₀(E×W³) = 2` (D-sketch) and polynomiality]; `c₁^{(U)} = 0` [D]; `c₁^{(3)} = 0` [P given Q2, Q8, F]; **value [I], carried exactly** | C1 §0, §6 |
| **L4** (one functional) | The k = 4 light side is CLS-D's term-wise s-channel functional `Ω_k = Σ_O p_O ω_k[O]` evaluated on the light OPE content, the same object that acts on the heavy blocks. On every tree term `X`: `Ω_k[X] = D_LR[X] + 𝒜_k[X]`, `𝒜_k` the u-channel Regge arc (Mellin contour at infinity), with `𝒜_k[X] = 0` if X's Regge spin `n_X < k` and, at `n_X = k`, a linear functional of X's top Regge coefficient only. There are no separate collinear or ν-boundary terms. | — | ✓ | [D given C′ + Mellin convergence]; arc structure [DV] on a Witten family | LIGHTSIDE §0–§1; LEMMAS §1.4 |
| **L5** (Lemma C′, Lemma A) | C′: on `Ω = {Re t₁, Re t₃ > 0; t_i ∉ (0, ½]}` the swapped-form integrand `t₁ʷt₃ʷf_k G` is one single-valued analytic function, whose four one-sided boundary values on the double-commutator segment are the four Wightman orderings. Lemma A: every tree term with Regge intercept `J₀ < 4` (Einstein, `E × W³`, MFT) has `D_LR = 0` pointwise in angle for `Ψ₄`, `Ψ₅`; the collinear set `a = b` carries nothing for them. Einstein instance: derived, **not numerically evaluated**. | — | ✓ | C′ [D] (closed form [DV] 2·10⁻⁹); Lemma A [D given C′] = [D]; MFT [P]; mechanism [DV] on free scalar, free fermion, spinning Witten diagrams | LEMMAS §1, §4; EIN §0 |
| **L6** (F: R⁴ normalisation) | `−𝒜₄[R⁴] = 2g₄`, i.e. `−Ĉ₄^{ds}|_{R⁴} = 2g₄` (the R⁴ value is ν-independent, so the de-shift weight enters only through `W(0) = 1`); the all-plus `ĝ₄` gets 0 | — | ✓ | **[D given Lemma A (term-wise, J₀ = 3 exchange), H_EFT(a) ν-paired at leading order, CLS-D J = 0 heavy dictionary]**; scalar analogue of each step [DV] | F_NORM §0, §2.3–2.5, §4 |
| **L6′** (a = b) | the borderline R⁴ arc's light function is integrable against `C_ν` at the collinear set (growth slower than `Δφ^{−7/2}`); the kernel bound and criterion are [D], the `J₀ < 4` part is [D] | — | ✓ | **[I]** (located: same light-cone edge of the Regge corner as LRC\*; LRC_CLOSE controls only the kernel side) | LEMMAS §4; LRC §6.5 |
| **L7** (H_EFT) | (a) the term-wise (P) value of a heavy exchange's double-trace part equals its actual OPE-side value **in ν-paired form** (false pointwise in angle [NE]); (b) the light double-trace data are those of a tree-level bulk EFT with a truncated derivative expansion, with no conspiracy in the heavy sum / contact tower (CHMRS 2106.10274 §5.1's standing assumption; also CLS-D §2.2) | — | ✓ | (a) [DV] scalar-Mellin toy, [D-sketch] spinning; **(b) [A]** | F_NORM §0, §3.4; RQ-S §2.2 |
| **L8** (P) | (P) for R⁴: the per-ν R⁴ light value is constant [D]. (P) for `g₅`: the `D²R⁴` light value of `Λ²(∂Ĉ₄ − Ĉ₅)^{ds}` cancels iff the `++` `J ≤ 3` heavy-image correction vanishes as `m → ∞`. The measured `s_d ≈ 0.51` is `m`-independent, so exact cancellation and the measured correction cannot both hold; the consistent reading is the double counting of RU issue 8. **The δ = 0.2 upper end holds under either reading** (`s_fail > 0` in every channel). | — | ✓ | R⁴ [D]; **g₅ [D-sketch], with tension** | F_NORM §2.6 |
| **L9** (signs) | `κ > 0`: sign of `Ω^Q_up`'s `R²` light coefficient (flat value 0.109908); `c_L > 0`: sign of `Ω^Q_low`'s `R²` slope (0.255 or 0.49 depending on the L8 reading). Needed only if `(t₄, α_T) ≠ 0`. | — | ✓ | [D-sketch + M] | UT §3.1; RQ-S §2.3–2.4 |

**What changed relative to UT's body.**
* **H_K:** from "[D-sketch], [D given H_ip (J₀ < ½)], no criterion short of H_ip" to [D-sketch] **with no [I] link**, under the known
  `J₀ ≤ 1`, margins ½ / 1, provided the functionals carry the cutoff `θ` (HK → LEMMAS → LRC; HK's `J₀ < 5/2` for k = 4 is corrected
  to `J₀ < 2`).
* **K1′:** [D-sketch] → [D]. **K2ν:** new, [D]. **LC\*:** not needed.
* **Light side of Q:** from "D_LR on J₀ < 4 plus term-wise (F)+(P) on J₀ ≥ 4, *defined* piecewise [D-sketch]" to one functional,
  `Ω₄ = D_LR + 𝒜₄` [D]; EIN's "J₀ = k values sit at a = b or in the ν integral" is withdrawn.
* **(F):** [D-sketch] → [D] (calibration), with the sign `𝒜₄[R⁴] = −2g₄`.
* **Lemma C:** withdrawn; **C′** [D]; Lemma A [D given C′]; Einstein instance [D], unevaluated (was "[I]").
* **De-shift:** "[NE] to t = 30; C5 spinning coefficients [I]" → ε-DS [DV] for all t; the C5 "[I]" was a stale label (coefficients
  [P + DV] since Stage 5); the upper end is exactly independent of the C5 weight [D].
* **`c₁`:** "free parameter" → one uncomputed universal constant, with exact bounds for every value.

---

## 5. Per-state positivity: evidence (not an input)

**Disclosures.**
* **Chronology.** The EE structure (`q² → EE`) and the Abel gauge (canonical → min-norm) were fixed *after* heavy-state failures had
  been seen. `Ω^EE` has been fixed since 20 Nov; `J = 18, 20, 24, 28` at `Δ_g = 300` are in-sample. In the canonical gauge the
  `J = 24` (300) state fails (`n = 1.55` spline; `−22 … −40` direct). The statement is about the min-norm functional specifically.
  `Ω^Q` was fixed on 22 Nov. (RD-S §2.)
* **Finite sets.** Positivity is computed on finite sets of states, inside the near-gap chain / hybrid model [NE + M]. The heavy
  spectrum is infinite in `J` and `m`; the rest is covered by the hypotheses below, several [M] and one [I].
* **Not certified.** The flat NF02 certificate proves `𝔉 ≥ (3/20)N` (not 0.200) on all flat states **except** `+−`,
  `J ≥ 400`, `M ≤ m < 1.633M`, which is float-checked only (min `𝔉/N = 0.765` to `J = 20000`) [NE]; the dense float minimum is
  `0.2000` [NE] (NOFORWARD §0, §4). UT and RD-S's "min 0.200 [certified]" overstate this.

### 5.1 Computed sets

| item | content | grade | source |
|---|---|---|---|
| D, `Δ_g = 1000` | even `J ≤ 24`, 28, 32, 36, 40, 48; odd `J ≤ 31`; masses twist edge → 5Λ. Min `+0.2009` (`J = 5`, `+−`, twist edge `1.00499Λ`); hybrid baseline `+0.20090` reproduced independently. `J = 48`: `+0.915 ± 0.08` (direct sub-row). `J = 15`: quoted `+0.869` is a grid minimum from `1.01Λ`; hybrid at its twist edge `≈ +0.69` [M]. | [NE + M] | RD-S §5.2; REPRO D3, D4 |
| D, `Δ_g = 300` | even `J ≤ 24`, 28, 32; odd `J ≤ 15`. Min `+0.2038` (`J = 5`, `+−`, `1.0166Λ`; hybrid `+0.2043`, AdS correction `−0.0005` not independently reproduced). Large J: `(336, 28)` `+0.759` (gap state; `0.751 … 0.772` under stresses), `(336, 32)` `+0.837` — the lowest computed large-J value is now `J = 28`. `J = 0, 2` still on the superseded small-ν treatment (margin ≈ 1). `J = 16`: quoted `+0.884` is near-gap only; plausibly `≤ 0.80` at large mass. HLJ_CLOSE's chain protocol gives `+0.2054` for the J = 5 minimum; the adopted `+0.2038` is RD-S's protocol. | [NE + M] | RD-S §5.3; HLJ §4; REPRO D5 |
| Q, both `Δ_g` | measured: `++` J = 0, 2, 4, 6, 8; OO J = 0, 2; `+−` J = 4–9, 11. Binding pair `++`/OO J = 0 at the gap: FF/N `+0.1427 / +0.1428` (δ = 0.2), full-set minimum `+0.1426 ± 0.0017`. Controlled by the derivative-piece coefficient `s_d = +0.5096 / +0.5077` against `s_fail = 1.776` (margin 3.48). Unmeasured channels have `s_fail ≥ 11`, i.e. ≥ 19× the largest measured `|s_d|` (0.59). | [NE + M] given H_RS | RQ-S §3.1; HRS §6 |

### 5.2 Hypotheses covering the rest (evidence structure; not inputs of the implication)

| # | statement | for | grade | source |
|---|---|---|---|---|
| H_RS | the analytic `1/Δ` expansion ("radial engine") of the EE k = 4 heavy image equals the physical crystal-cycle integral, at Q's binding state | Q | **[D given P_cyc]** (preliminary; was [I]). `P_cyc` (the physical cycle lies in the twisted-homology class of the cut real cycle `D′`) is [DV] for the scalar through its periods, **not proved**, and not tested by either reproduction (both integrate over the real section). Far region on the thimbles: `≤ e^{−66}` (300) pointwise at `r ≥ 0.3`, `≲ e^{−26}` after the volume factor `≲ e^{40}`, against `e^{−8.5}`; the `r ≥ 0.2` figure `e^{−37}` must not be quoted. Covers `++` J = 0, 2 and OO J = 0; other measured Q channels use the same engine without the bookkeeping, at larger margins. | HRS §0; RD §0, §5; UT RD box |
| H_M | the near-gap chain is accurate on the computed states | D | [M], validated [NE] on rows with direct data; `K₂` band `J²/m = 1.6–2.2` untested | RD-S §4 |
| H_edge | light-cone-edge window `x ∈ (0.8, Λ/m]` | D | [NE + M] where computed; [M] elsewhere | RD-S §4 |
| H_OO | parity-odd even-J block via `Γ_OO = Γ_++` | D, Q | [D given the identity]; identity [P given helicity lemma] at order 0, [NE] beyond | RD-S §4 |
| H_R | beyond-leading-Regge at `m ≳ 1.3Λ`, large `b̂` | D, (Q) | [M] | RD-S §4 |
| H_LJ(1000) | near-gap large spin at 1000: `H_AG ∧ H_win(refined) ∧ H_edge-amp`. The sign-free `H_win` is **violated at (1100, 100)** (`−2154 ± 43` vs `T_B = 1409`); the refined form (leading-Regge-shaped `D_h` class removed) rests on a [D-sketch] step and on direct nodes at `ν ≥ 13`; the engine cannot decide it at `J²/Δ ≳ 9` | D | **[D-sketch + M]** | RD-S §4; HLJ §5–6; REFEREE M9 |
| H_LJ(300) | the same at 300: bulk [D-sketch/M]; **Regge window `J ≥ 34`: [I]** (no graviton data; engine series asymptotic); `(336, 28/32)` windows `ν < 8.3`: [NE + M] (was [M]) | D | **[I]** (window) | HLJ §0, §6 |
| H_LJ^Q | large-J near-gap corrections below Q's tolerance (`≥ 1.1·10⁴ c₂`) | Q | [M], margin ×10–10³ | RQ-S §3.3 |
| H_fwd^Q | unmeasured forward channels; the 300 `J`-drift | Q | [NE + M] given H_RS; drift [I], harmless at δ = 0.2 | RQ-S §3.3 |
| H_smear^Q | smeared graviton corrections at the binding states (`≤ 3·10⁻⁵` of B₂ vs ≈ 24 %) | Q | [NE + M] | RQ-S §3.3 |
| HYB | scalar skeleton + leading Regge elsewhere; `β > 30` (re-scanned to β = 45 with the new C5 weight: min `+0.9866` at Λ = 100) | Q | [M] (images [DV]) | DESHIFT §2.5 |
| H_cont | continuity in `m` between sampled masses (exception: `J = 15` at 1000) | D | [NE] | RD-S §4 |
| Flat | bulk-point limit of `Ω^EE` = stripped NF02 (`EE ≡ q²` on CLS-D's plane [P]; `EE/q² → 1` only as `m → ∞` in the full contraction [NE]); NF02 certificate as stated above | D | [P given map] + [NE]; certificate as in §5 | RD-S §2 |
| `σ_hU` | double counting of the helicity-odd correction | — | [I]; **not used** (decides only δ = 0.05 and the withdrawn closure) | UT §3.2 |

**Net:** at `Δ_g = 1000` every covering hypothesis is [M] or better (weakest: H_LJ(1000) [D-sketch + M], H_M, H_R, HYB [M]). At
`Δ_g = 300` one is [I] (H_LJ(300) window, `J ≥ 34`), which is why the headline is `Δ_g = 1000`. H_LJ(300) is needed only by item (1)
at 300; item (2) at 300 no longer needs it (§3, Dependence).

---

## 6. What is now derived (summary of the derivations behind §4)

### 6.1 The light side is one functional [D] (LIGHTSIDE §1.2)
In Mellin form on the swapped contour `C × C`, the kernel `K(δ₁₂, δ₁₃)` is entire in `δ₁₂` (`|U/V| ≤ 1` on `C`). Three evaluations of
the same object: the whole-function fixed-angle integral (`D_LR`, "t outside"), the straight Mellin line, and the term-wise
s-channel sum (close `δ₁₂` to the left: double-trace double poles at `δ₁₂ = −n`, single-trace poles). Cauchy at fixed `δ₁₃`:
line − term-wise = −(arc at infinity). After 1↔3 (anti)symmetrisation the `D_LR` value of any microcausal term vanishes as an integrand
identity, so `Ω_k[X] = 𝒜_k[X]`. Along the arc `M ~` Regge growth and `K ~` corner behaviour: `𝒜 = 0` for `n_X < k`; at `n_X = k`,
`𝒜 ∝` top Regge coefficient × a universal kernel — CHMRS's `eftformula` in position space. Checks: term-wise = whole-function to
`10⁻¹¹–3·10⁻⁵`; a borderline contact gives a finite light value (`−8.7642·10⁻³` at P1) equal to its top-Regge term to `≤ 7·10⁻⁷`; a
borderline exchange (H7 analogue) to `8·10⁻⁸` [DV].

### 6.2 Lemma C′ [D] (LEMMAS §1.2)
Lemma W [P]: for `Re s₁, Re s₃ ≥ 0` (`s = 1/t`) and `y ∈ (0, 2]`, neither root of `λ² − (1+U−V)λ + U` lies in `[1, ∞)`, by the exact
identity `8y[λ² − (1+U−V)λ + U] = (s₁+4λ−2)(s₃+4λ−2) − 8(2−y)λ(λ−1)` and `|s+4λ−2| ≥ 4λ−2` for `λ ≥ 1`. So `|ρ(z)|, |ρ(z̄)| < 1`
on `Ω`, the s-channel OPE converges there [Q], `log U` continues as `L = Log(2−s₁) + Log(2−s₃) − log 8y`, and the integrand is one
analytic function on `Ω` (Weierstrass). Tree-level terms: Witten diagrams are analytic on the s-channel sheet away from `{0,1,∞}`; the
bulk-point singularity is on the sheet reached by crossing `[1, ∞)`, which `Ω` never reaches.

### 6.3 (F) by calibration [D given inputs] (F_NORM §2.3)
Auxiliary tree: a heavy scalar χ coupled through `cχ(W⁺)² + c̄χ(W⁻)²`. Its exchange has Regge spin 3 < 4, so Lemma A (term-wise) gives
`p_χĈ₄^{ds}[G_χ] + Ĉ₄^{ds}[DT(W_χ)] = 0`. The heavy block's value is `2u²|c^{++}|²(1 + O(m⁻²)) = 2g₄^χ(1 + O(m⁻²))` by the CLS-D
J = 0 heavy dictionary and the flat `B₄⁽¹⁾` sum rule. Integrating out χ gives an AdS R⁴ contact with couplings `g₄^χ ∝ |c|²`,
`ĝ₄^χ ∝ c²`, plus curvature terms of Mellin degree ≤ 3 (no arc). By H_EFT(a) (ν-paired, leading order),
`Ĉ₄^{ds}[DT(W_χ)] = 𝒜₄[R⁴(g₄^χ, ĝ₄^χ)] + O(g₄^χ/m²)`. The arc is m-independent at fixed coupling, so `𝒜₄[R⁴] = −2g₄` exactly;
rotating the phase of `c` shows the `ĝ₄` coefficient is 0. Checks [DV]: CHMRS's contact functional reproduced exactly; Penedones'
top coefficient vs the heavy-density coupling, extrapolated ratio 0.9999990; calibration identity to `5·10⁻¹⁰`; Lemma A on heavy
σ-exchanges `≤ 1.7·10⁻⁴`; (P) polynomial ratios to `10⁻⁷`. Negative [NE]: at fixed celestial angle the calibration fails at O(1), so
H_EFT(a) holds only ν-paired.

### 6.4 H_K by Tonelli with cut-off weights [D-sketch] (HK §1, §3–§5; LEMMAS §3; LRC §1–§3, §6.1)
1. Corner count: for a kernel `|K| ≲ r^a η^c` on `B_c`, R1′ gives absolute convergence iff `c < a + 1 − J₀`; R3 turns this into
   `Σ_O∫|K||dDisc G_O| < ∞`. Thresholds `a = 0` (k = 2, 3), `a = 2` (k = 4, 5).
2. No-go [P]: a compactly supported ν-weight has raw growth `c = ½` (k = 2) / `5/2` (k = 4) by Paley–Wiener, so `J₀ < ½`.
3. Strip analyticity [D]: the de-shift identity itself forces `χR_j` to continue into the strips `−1 < Im ν < 0` (k = 2) and to depth
   2 (k = 4) with uniformly `L¹` lines (Lemmas S2, S4); the entire cutoff then spoils the cancellation only by `ε_θ`:
   `|I_θ − I| ≤ ε_θ P_{−½}(cosh β)` (k = 2) (Lemma C).
4. Beyond leading Regge (LRC\*): the exact per-ν kernel is a continuous-spin Casimir eigenfunction; its light-cone-edge correction is
   the factor `F(ρ; ν)`, analytic below `Im ν = 3/2`, so the contour shift survives to all orders in ρ. k = 2: `c = −½`, `J₀ < 3/2`.
   k = 4: the shift crosses the pole of `F` at `3i/2`, residue `∝ f_pos(3i/2)ρ²η`, `c = 1`, `J₀ < 2`. Checks [DV]: exponent of
   (exact − leading Regge) `−0.50` for entire k = 2 and regular k = 4 weights, `+1.00` for a k = 4 weight with a pole at `3i/2`
   (residue formula to `10⁻⁶`).
5. EE: ⟦EESTEP⟧.
6. `Ĉ₅^{ds}` (raw `c = 3/2 < 2`) and `Ĉ₃` (no de-shift): plain Fubini under `J₀ < 3/2`.

### 6.5 `c₁` and the exact bounds [D] (C1 §3–§4)
At `O(t₄)` the d = 3 one-point matrix is pure helicity flip; `EE` annihilates every flip configuration of the plane-wave kernel ([P]
dual line; [DV] `≤ 4·10⁻¹⁴` light-ray line). Hence the u-channel T-exchange part of `c₁` vanishes [D], and the k = 3 part vanishes
[P given Q2, Q8, F]. The remaining contact-like part is a polynomial in ν²; flat power counting forbids `ν^{2n}`, `n ≥ 2`, and the
`ν²` coefficient is the flat `O(ĝ₃)` MHV term, which vanishes: a single constant [D given `J₀(E×W³) = 2`]. Its value is the `1/L²` remnant of a vanishing flat quantity:
uncomputed [I], natural size `10⁻³–10⁻¹` [M], no bound [I].
* Item (1): the sum rule `K R²Λ⁸/57600 − 16I_DΛ²c₁t₄ − 16L_A ≤ 0` with `t₄, α_T` free is a quadratic inequality; its root gives
  `R ≤ R₀(x + √(1+x²))`.
* Item (2): `U(t₄, α_T) = 2L_AdS − κR²Λ⁸/57600 + 2c₁t₄Λ²I_Q`; maximising over `(t₄, α_T)` at `t₄* = 57600 I_Q c₁/(κΛ⁶)` gives
  `2L_AdS + 57600 I_Q²c₁²/(κΛ⁴)`; `57600 I_Q²/κ = 541017`.

### 6.6 De-shift weights explicit and published [DV] (DESHIFT §2)
Shifted form of the identity: shifting each `j`-term to `μ = ν − ij` gives `W[u](t) = Σ(residue terms) + ½∫E(μ)cos(μt)dμ` for every
real `t`; imposing "residue part = Φ" exactly leaves `|W − Φ| ≤ ε_W := ∫₀^∞|E_even|` for all `t ≥ 0`, hence
`|image − target| ≤ ε_W P_{−½}(cosh β)` for all β. Independent reproduction: D's NF02 Abel weights (own discretisation, FFT, CGLS) and
Q's C4 (different method) to `10⁻¹⁶`.

---

## 7. Gap dependence (functionals held fixed) [D]; constants [DV] (DG §1)

`a = 3/(2Λ)`; for a polynomial weight `L(ψ; a) = ∫₀¹ψ/(p² + a²)` is exact term by term (`I₀ = arctan(1/a)/a`, `I₁ = ½ln(1 + a⁻²)`,
`I_n = 1/(n−1) − a²I_{n−2}`).
* `16L_A = A ln(2Λ/3) + B_∞ + b₁/Λ + (b₂L ln(2Λ/3) + b₂)/Λ² + O(Λ⁻³ln Λ)`, `A = 3022.3377`, `B_∞ = −5327.2096` (`= B_AdS + A ln(3/2)`,
  `B_AdS = −6552.662`), `b₁ = 19293.4`, `b₂L = −13549`, `b₂ = −14705`.
* `R₀ (K = 1) = 1.3194·10⁴ √(ln(Λ/8.7415) + 6.384/Λ + O(ln Λ/Λ²))/Λ⁴` (relative error of the bracketed form `≲ 5·10⁻⁵` at 300, 1000).
* `2L_AdS = 17.2407 ln(2Λ/3) − 19.2782 + 10.314/Λ + (860.3 ln(2Λ/3) − 2508.6)/Λ² + O(Λ⁻³ ln Λ)` (absolute error `≤ 10⁻³` with
  all terms; the truncation at `1/Λ` gives 72.017 at 300, on the unsafe side, and must not be quoted as a bound).
* `B_eff` drift `−19.2197 → −19.2647`: the finite-Λ tail of the AdS pole, converging to the flat `−19.2782`.
* Error terms: `x = 19941|c₁|/(Λ²√(16L_A))`; H17 `±212.5/Λ²`; `δk(0) = −68.08/Λ²(1 + O(Λ⁻¹))`; Q's `c₁` term `541017c₁²/(κ_flat Λ⁴)`.
* `Λ⁴ → Δ_g⁴` is not innocuous: `−1.99 %` at 300 (REPRO D6; DG's −2.2 % is wrong).

These are values of the implication; they are bounds at a given `Δ_g` only where (Pos) holds.

---

## 8. The `Δ_g ≥ 300` corollary (stated only with its hypotheses)

> **Corollary (conditional, [D-sketch + M]).** Assume (i) A1, A2 and the inputs of §4 at every `Δ_g ≥ 300`; (ii) (Pos) at
> `Δ_g = 300` for the relevant functional — which at 300 includes the open window of H_LJ(300) [I]; (iii) H_mono: at fixed `(J, m/Λ)`
> the beyond-leading-Regge correction (in `c₂`) grows at most like `m`, and in the Regge window it is bounded uniformly in `Δ_g` [M];
> (iv) H_cont for states between the 300 and `Δ_g` twist edges [NE]; (v) H_LJ(Δ_g) [D-sketch + M]. Then (Pos) holds at every
> `Δ_g ≥ 300`, and §7's forms with their error terms are bounds there.

* It is not a proof: corrections are measured only where direct data exist. Proxies reuse computed weights and corrections
  (validated on the 100/300/1000 verdicts): correction-driven margins improve monotonically; the binding baseline decreases to the
  flat 0.200; the window tolerance grows like `Λ^{1.0–1.2}` (DG §2).
* The sign-free large-spin window closes only from `Δ_g ≈ 275` **with the [D-sketch] refinement** (margin ×1.1 at 300), and from
  `≈ 440` without it. "300 is the lowest gap where the method closes" holds only with the refinement.
* Anchoring at 1000 instead would avoid H_LJ(300) but DG does not state that instance; not used.
* Nothing is claimed below 300. At `Δ_g = 100` positivity fails at model level for near-gap `J ≳ 18`.

---

## 9. ABJM: a normalisation check, outside A1

* `g₄` normalisation: type II gives `g₄M⁶/8πG = 2ζ(3) = 2.404` (CHMRS 2102.08951; Stage 0 Virasoro–Shapiro check [DV]).
* For ABJ(M) the same coefficient is *determined* in the CFT from localisation plus the stress-tensor-multiplet four-point function:
  IIA regime Binder–Chester–Pufu 1906.07195 (key `BCPW`); M-theory regime Chester–Pufu–Yin 1804.00949. Our `g̃₄` is a different CFT
  object (a 4d light-ray functional on `⟨TTTT⟩`); identifying the two needs (F), 4d restriction of the 10d amplitude, and SUSY relations.
* With `Λ ↔ ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.2156 λ^{1/4}` (leading order in `1/(ML)`), 2ζ(3) lies at 3.34 % (300) / 2.59 % (1000) of
  the upper end. This small fraction mostly reflects the 4d IR logarithm `17.24 ln(2Λ/3)`; against CHMRS's ten-dimensional
  maximal-supergravity window `0 ≤ g₀M⁶/8πG ≤ 3.000`, 2ζ(3) sits at 80 %.
* `ĝ₃ = 0` by ten-dimensional type IIA supersymmetry.
* ABJM violates A1 (light KK and R-singlet single traces) and its flat limit is ten-dimensional, so this checks sign and
  normalisation only. It is not evidence for the theorem.

---

## 10. What is not claimed

* No unconditional bound; (Pos) is not established.
* No `Δ_g` other than 300 and 1000, except through §8's corollary; nothing below 300; not `Δ_g = 100`.
* No loop-level statement.
* No bound on `ĝ₄` (all-plus R⁴) or on `g₅`; no cubic closure `|ĝ₃|²Λ⁸ ≤ …` from item (2) and no ratio `|ĝ₃|²/g₄` (the earlier
  closure was withdrawn: its slope 0.49 is 1.9× CHLPS's flat optimum 0.26). At most a remark: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257`
  [D-sketch + M].
* Not the δ = 0.05 functional (margin `+0.0060`, depends on `σ_hU`).
* No explicit lower constant.
* The Einstein instance of Lemma A is derived, not evaluated numerically.
* Nothing about theories violating A1.

---

## 11. Inconsistencies found among the sources (and what this document does)

1. **Box order.** UT's header says boxes are "newest first", but they were appended oldest first (LEMMAS corrects HK; LRC corrects
   both; DESHIFT is last). This document treats the last box as newest.
2. **DS for Q.** LRC's box lists "DS for Q [NE]" among H_K's weakest links; the later DESHIFT box replaces DS by ε-DS [DV] for the
   published continuum weights. Adopted: [DV].
3. **H_K criterion for k = 4.** HK: `J₀ < 5/2` (`7/2` given NE); DESHIFT §2.4: `7/2` "without the NE input"; LRC: `J₀ < 2` for generic
   weights beyond leading Regge. DESHIFT's 7/2 is a leading-Regge statement; adopted: LRC's `J₀ < 2`.
4. **"Minimum 0.200 [certified]"** (UT §4, RD-S §2, old letter) vs NOFORWARD: the certificate is `𝔉 ≥ (3/20)N`, with 0.200 a dense float
   minimum [NE], and one corner not certified. Adopted: NOFORWARD's statement.
5. **`ε_{c₁}`.** UT §2, DG, RQ-S, old letter: `|ε_{c₁}| ≤ 0.57|c₁| / 0.059|c₁|`; C1: superseded (drops `−κR²`), exact form
   `6.81·10⁻⁵c₁² / 5.44·10⁻⁷c₁²`. Adopted: C1.
6. **Sign of (F).** LIGHTSIDE §4 and its UT box write `𝒜₄[R⁴] = 2g₄`; F_NORM: `𝒜₄[R⁴] = −2g₄` (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`). Adopted: F_NORM.
7. **a = b integrability.** LEMMAS and LRC keep the borderline R⁴ arc's a = b integrability [I]; F_NORM grades (F) [D] and LIGHTSIDE
   grades the light side [D] without discussing it; UT's boxes do not carry it into the net grade of item (2). This document lists it
   as L6′ [I] and makes item (2)'s net grade [I] in that one link.
8. **Far-region number for H_RS.** RD box/letter: `≤ e^{−37}` at `r ≥ 0.2`; RD §3.2 and REFEREE M10: that does not survive the volume
   factor `≲ e^{40}`. Adopted: `e^{−66}` at `r ≥ 0.3`, `≲ e^{−26}` after the volume factor.
9. **RD fit accuracy.** UT RD box "within 2 %"; RD §0/§2.4 "within 3 %". Adopted: 3 %.
10. **Q binding value at 1000.** RQ-S §3.1 "`+0.1426 ± 0.0002`"; HRS, REPRO, REPRO_ENGINE `+0.1428`. Adopted: `+0.1428` for the binding
    state, `+0.1426 ± 0.0017` for the full-set minimum.
11. **D's 300 minimum.** RD-S `+0.2038`; HLJ_CLOSE's chain protocol `+0.2054`. Adopted: `+0.2038` (RD-S protocol), with the other
    recorded.
12. **Most fragile large-J state at 300.** RD-S/old letter: `J = 32`, `+0.783`; HLJ_CLOSE: `(336, 28)` `+0.759`, `(336, 32)` `+0.837`.
    Adopted: HLJ_CLOSE.
13. **"Every unmeasured channel needs `|s_d| ≥ 17–30×`"** (UT §3.2) vs RQ-S §3.1 `s_fail ≥ 11`: `11/0.59 = 18.6`, `11/0.51 = 21.6`.
    Adopted: "≥ 19× the largest measured `|s_d|`".
14. **`s_d` law.** DG `0.500 + 1440/Δ²` misses its anchors; REPRO: `0.5076 + 147/Δ²`. Adopted: REPRO.
15. **Λ⁴ → Δ_g⁴.** DG `−2.2 %`, REPRO `−1.99 %`. Adopted: REPRO.
16. **`δk(0)` at 100.** RD-S vs DG differ (REPRO D7); irrelevant (100 not claimed).
17. **`X^Q` not published.** Ω^Q_up's smeared k = 2 de-shift weight is computed (RESULT_Q_REVISION §1.3) but absent from `prl/data/`,
    while DESHIFT §0 says "Q3-full-δ0.2" and "the kernels" are exported. The definition is still fixed (min-norm solution of a stated
    equation), but the referee's M6 request is only partly met.
18. **D's functional and "published weights".** For Q the published continuum weights *are* the definition; for D they are a grid
    realisation of the exact min-norm `X` (no continuum k = 2 weight exists; the stored grid satisfies the identity only to `t ≈ 22–25`).
19. **(P)/g₅.** F_NORM finds exact g₅ cancellation incompatible with the measured `m`-independent `s_d`; the deciding test (engine
    `++` J = 0 derivative-piece correction at `m = 3–10Λ`) has not been run. The δ = 0.2 upper end survives either reading.
20. **CHMRS swamp equation numbers.** LITERATURE cites eqs. 3.35, 3.43 for the normalisation; F_NORM could not match "(3.43)" (finds
    the lines near (3.31)/(3.39)). The letter cites the section only.

---

## 12. Source table

Abbreviations: UT = UNIFIED_THEOREM.md; RD-S = RESULT_D_STATEMENT.md; RQ-S = RESULT_Q_STATEMENT.md; DG = DG_SCAN.md; C1 =
C1_COMPUTE.md; HK = HK_CLOSE.md; LEMMAS = LEMMAS_CLOSE.md; LRC = LRC_CLOSE.md; LIGHTSIDE = LIGHTSIDE_DIST.md; F_NORM; DESHIFT =
DESHIFT_CLOSE.md; HRS = HRS_CLOSE.md; RD = RD_CLOSE.md; HLJ = HLJ_CLOSE.md; EIN = RESULT_Q_EINSTEIN_CHECK.md; REPRO, REPRO_ENGINE
(independent); NOFORWARD = NOFORWARD_FUNCTIONAL.md; LIT = prl/LITERATURE.md; DATA = prl/data/README.md.

| number / statement | value | source |
|---|---|---|
| `8πG = 256/c_T` | — | UT §1; REPRO §1 |
| Λ at 300 / 1000 | 298.50 / 998.50 | UT §1; DG §0 |
| `16L_A` | 10734.61 / 14339.59 | DG §1.2; RD-S §3; REPRO §1 |
| `R₀` (K = 1) | 3.13206·10⁻⁶ / 2.89126·10⁻⁸ | REPRO §1; C1 §4.1 |
| `|ĝ₃|²Λ⁸ ≤` | 1.07·10⁴ / 1.43·10⁴ | RD-S §3 |
| `I_D`, `1920I_D` | 10.385947, 19941.0 | C1 §1.1, §4.1 |
| `x/|c₁|` | 2.16006·10⁻³ / 1.67026·10⁻⁴ | C1 §4.1 |
| exact (1) at `|c₁| = 1, 10, 100, 1000` | table §3 | C1 §4.1 (re-checked) |
| linearisation good to 1 % | `|c₁| ≤ 70.5 / 911` | C1 §4.1 |
| `I_Q`, `κ`, `57600I_Q²/κ` | 1.0160366, 0.109908, 541017 | C1 §1.1, §4.2 |
| `c₁²` coefficients | 6.8145·10⁻⁵ / 5.4428·10⁻⁷ | C1 §4.2 |
| exact (2) at `|c₁| = 100, 1000` | 72.72, 140.2 / 92.819, 93.36 | C1 §4.2 |
| `c₁` natural size | 10⁻³–10⁻¹ [M] | C1 §0 |
| `2L_AdS` | 72.041 / 92.814 (72.0406 / 92.8135) | DG §1.2; REPRO §1 |
| `B_eff` | −19.2197 / −19.2647 | DG §1.2 |
| `A_Q`, `B_∞^Q`, `b₁`, `b₂L`, `b₂`, `Λ₀^Q` | 17.240713, −19.278176, 10.3137, 860.35, −2508.6, 4.5889 | DG §1.1 |
| `A`, `B_∞`, `B_AdS`, `b₁`, `b₂L`, `b₂`, `Λ₀^D` | 3022.3377, −5327.2096, −6552.662, 19293.4, −13549, −14705, 8.7415 | DG §1.1 |
| cubic Λ-form | 1.3194·10⁴, 8.7415, 6.384 | DG §0, §1.3 |
| quartic Λ-form | 17.2407, −19.2782, 10.314, 860.3, −2508.6 | DG §0 |
| truncation 72.017 at 300 | — | REFEREE minor 7; DG §1.3 |
| H17 band | ±212.5/Λ² = ±0.24 % / ±0.021 % | DG §0; REPRO D2 |
| `δk(0)` | −7.666·10⁻⁴ / −6.830·10⁻⁵; −68.08/Λ² | RD-S §3; DG §0 |
| `Λ⁴ → Δ_g⁴` | −1.99 % | REPRO D6 |
| cutoff parameters, `ε_θ` | Q (150, 6), `≲ e^{−180}`; D (Λ+100, 10), `≲ 10⁻⁴¹` | HK §0, §4.4 |
| H_K criteria | `J₀ < 3/2` (k = 2, 3, 5), `J₀ < 2` (k = 4) | LRC §0, §6.1 |
| exact kernel check | 3·10⁻⁸ (spread 2.6·10⁻⁸ – 1.5·10⁻⁶) | LRC §1.2 |
| LRC\* exponents | −0.50, +1.00; residue to 10⁻⁶ | LRC §3.2 |
| `ε_W` | 1.2·10⁻¹³ / 4.2·10⁻¹¹ / 2.4·10⁻⁷ | DESHIFT §2.3; DATA §3 |
| continuum vs XP | 1.9·10⁻¹⁷ / 1.1·10⁻¹⁶ | DESHIFT §2.3 |
| XP / Stage-4 / XP-C5 fail at | t ≈ 41 / 12 / 20 | DESHIFT §2.3 |
| D Abel reproduction | ≤ 1.5·10⁻⁶ (ν < 2), ≤ 3·10⁻¹¹ (ν > 25) | DESHIFT §4.1 |
| C5 tolerance `G*` | 1.18·10⁵ c₂ (300, new weight) | DESHIFT §1.3 |
| HYB re-scan | +0.9645 → +0.9866 | DESHIFT §0, §2.5 |
| Q3 δ nominal | dense min 0.1966 | REPRO D9 |
| LIGHTSIDE checks | 10⁻¹¹–3·10⁻⁵; −8.7642·10⁻³; ≤ 7·10⁻⁷; 8·10⁻⁸ | LIGHTSIDE §0 |
| C′ closed-form check | ≤ 2·10⁻⁹ | UT LEMMAS box; LEMMAS §1.3 |
| F_NORM checks | 0.9999990; 5·10⁻¹⁰; 1.7·10⁻⁴; 10⁻⁷ | F_NORM §0 |
| Lemma A mechanism | ≤ 8·10⁻¹⁵ (scalar), ≤ 3.8·10⁻¹⁵ (fermion), 10⁻¹⁴–2·10⁻¹² (Witten) | RQ-S §2.1; EIN §0 |
| `K̄` | 2 log(1/Δφ) + 8.69 | LEMMAS §3.1 |
| D minima | +0.2038 [1.0166Λ] / +0.2009 [1.00499Λ] | RD-S §5; REPRO D3 |
| D hybrid | +0.2043 / +0.2009 (+0.20427 / +0.20090) | REPRO §1 |
| (336, 28), (336, 32) | +0.759 (0.751–0.772), +0.837 | HLJ §0, §4.2 |
| (1048, 48) | +0.915 ± 0.08 | RD-S §5.2 |
| J = 15 (1000), J = 16 (300) caveats | ≈ +0.69; ≤ 0.80 | REPRO D4, D5 |
| (1100, 100) violation | −2154 ± 43 vs T_B = 1409 | REFEREE M9; HLJ §5 |
| Q binding FF/N | +0.1427 / +0.1428 (+0.14270 / +0.14282; direct ρ: +0.14265 / +0.14282) | HRS §6; REPRO; REPRO_ENGINE §0 |
| Q full-set min | +0.1426 ± 0.0017 | RQ-S §3.1 |
| `s_d` | +0.5096 / +0.5077 (direct: +0.50967 / +0.50776) | HRS §0; REPRO_ENGINE §0 |
| `s_fail`, margin | 1.776; 3.48 | REPRO §2; REPRO_ENGINE §0 |
| truncation `|Δs_d|` | ≤ 0.015 / 0.0015 (×85 / ×850) | HRS §6 |
| ρ agreement | ≤ 3.4·10⁻⁷ / 3·10⁻⁸ | REPRO_ENGINE §0–§1 |
| sensitivity | 2.1·10⁻⁴ / 1.9·10⁻⁵; ~600× | REPRO_ENGINE §0, §5 |
| unmeasured channels | `s_fail ≥ 11` | RQ-S §3.1 |
| far region | e^{−66} (r ≥ 0.3), volume ≲ e^{40}, ≲ e^{−26}; thresholds e^{−8.5} / e^{−10.8} | RD §0, §3.2; UT RD box |
| far exponent on `D′` | `c = 2 artanh(r/2)` ≥ 0.336; fit within 3 % | RD §0, §2.3–2.4 |
| scalar periods | ≤ 6·10⁻⁹ (300), ≤ 1·10⁻⁹ (1000) | RD §0 |
| NF02 certificate | `𝔉 ≥ (3/20)N` except `+−`, J ≥ 400, M ≤ m < 1.633M (float 0.765) | NOFORWARD §0, §4 |
| `s_d` law | 0.5076 + 147/Δ² | REPRO D1 |
| window closes from | ≈ 275 (refined) / ≈ 440 (unrefined); ×1.1 at 300 | DG §0, §2.4 |
| Δ_g = 100 | fails at model level, J ≳ 18 | UT §6; DG §2.1 |
| ABJM | 2ζ(3) = 2.404; 3.34 % / 2.59 %; ML = 4.2156λ^{1/4}; 80 % of 3.000 | UT §7; REPRO §1; LIT §3 |
| withdrawn closure | 0.49 vs 0.26; 200 / 257 | RQ-S §5 |
| δ = 0.05 margin | +0.0060 | RQ-S §0 |
