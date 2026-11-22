# `c₁`: the `E × W³` light term of `Ω^EE` and `Ω^Q_up`. What it is, which part vanishes, and the exact (non-linearised) bounds

*22 Nov 2026, derivation agent, responding to `prl/REFEREE_REPORT.md` §M5. This note is unreviewed and nothing has been committed. No
existing file was edited. New code is in `ads4/code/c1/` and its outputs are in `ads4/code/c1/cache/`. `LIGHTSIDE_DIST.md`, written by
another agent, did not exist when this note was finished; this note is consistent with it as described in §1.2.*

**Grades.**
* [P] proved; [P given X] proved from the named inputs X.
* [D] derived; [D-sketch] derivation sketched; [DV] derived and checked numerically or symbolically.
* [Q] quoted from the literature.
* [NE] numerical evidence; [M] model or natural-size estimate; [A] assumption; [I] open.

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **Which prescription** | H7's: the light term is the Lemma-R (CHMRS App. E) residue of the full tree correlator, i.e. its `r^{1−k}` Laurent coefficient at the `u`-Regge point. It is mapped to `Ĉ′_EE` by D3 (the same plane-wave map for heavy and light) and normalised by the flat limit (E5). This prescription produces `−Ĉ′_EE|_E = 16πG/(ν²+9/4)`. | definition (§1) |
| **`k = 3` part of `c₁`** | **`= 0`**. Two independent routes: H7 route A (`[H₁,H₃]` has no term linear in the flip) and route B (Lemma R with no `Π₂` moment). | **[P given Q2, Q8, F]**; route B [P given Lemma R, `n = 2`, no-`Π₂` (D + NE)] |
| **`k = 2` part from the `u`-channel `T` exchange** (the "Regge residue of the `E × R³` `T` exchange" named by the referee and REVIEW_EE m5) | **`= 0`**, both its pole part and its off-shell vertex part. At `O(t₄)` the exchange is "flip on one Regge line, preserved on the other". `EE` annihilates every flip configuration of the plane-wave kernel exactly (dual line, [P]) and to `≤ 4·10⁻¹⁴` (light-ray line, [DV]). The `q²` control has an `O(1)` dual-line flip. | **[D]** (inputs at H7's level) |
| **What `c₁` actually is** | All of `c₁` comes from the *contact-like* part of the `O(t₄)` correlator in the `u`-Regge limit: the `s`- and `t`-channel `E × W³` exchanges plus the `O(ĝ₃)` quartic contact, in the preserved ⊗ preserved helicity sector. This part is a polynomial in `ν²`. Flat-limit power counting forbids `ν^{2n}` with `n ≥ 2`. The `ν²` coefficient is the flat `O(ĝ₃)` MHV-type term, which vanishes. **So the light term is exactly `16πG c₁t₄`, independent of `ν`, with `c₁` one universal pure number.** It does not depend on `Δ_g`, on the theory beyond `t₄`, or on `ν`. This confirms the project's ansatz and its `Λ⁻²` scaling, which REVIEW_EE m5 called "asserted, not derived". | form **[D given J₀(E×W³) = 2 (D-sketch) and polynomiality (E)]**; flat part `= 0` [P] |
| **Value of `c₁`** | **Not computed [I].** It is the AdS curvature (`1/L²`) remnant, at zero impact parameter, of a flat quantity that vanishes. No argument in hand forces it to zero (§3.5). Natural size: `|c₁| ~ 10⁻³ – 10⁻¹`. | value **[I]**; size **[M]** (not a bound) |
| **Rigorous bound on `|c₁|`** | **None** (§3.6). | [I] |
| **Exact `D` bound, every `c₁`** | `√(t₄²+α_T²) ≤ R₀(x + √(1+x²))`, with `x = 19941.0·|c₁|/(Λ²√(K·16L_A))`. The `t₄` interval is asymmetric (§4.1). | **[D]** (same grade as the `c₁ = 0` implication) |
| **Exact `Q` bound, every `c₁`** | `g̃₄Λ⁶/8πG ≤ 2L_AdS + 57600 I_Q² c₁²/(κΛ⁴) = 2L_AdS + 6.81·10⁻⁵c₁² (300) / 5.44·10⁻⁷c₁² (1000)`. **It needs no input from Result D**, only `κ > 0`, which Q already assumes. The old linear `0.57|c₁|` / `0.059|c₁|` dropped the `−κ|ĝ₃|²` term. | **[D]** given `κ > 0` ([D-sketch + M], as before) |

**Consequences for the headline (§5).**
* The referee's "Regge residue of the `E×R³` `T` exchange" is zero, and so is the whole `k = 3` part.
* What survives is one universal ν-independent constant. Its natural size moves the cubic bound by `≲ 0.02 % / 0.002 %` [M].
* The exact quadratic bound stays within 1 % of the linearised one up to `|c₁| = 70` (300) and `|c₁| = 911` (1000). At `|c₁| = 1000` it is 4.54× (300) / 1.18× (1000) the `c₁ = 0` value.
* **Result Q is decoupled from Result D entirely.** At `Δ_g = 300` the `c₁` term is below 1 % of `2L_AdS` for `|c₁| ≤ 100`.

---

## 1. Definitions and prescription

### 1.1 The object
`c₁` is the coefficient of a possible `O(t₄)` light term `16πG c₁t₄` in `−Ĉ′_{EE,ν}`, and collects the matching term in `Ĉ₃`
(RESULT_D_STATEMENT §3; UNIFIED_THEOREM §1; RU issue 14).
* Result D's sum rule with this term is
  `16L_A + 16 I_D Λ² c₁t₄ − K|ĝ₃|²Λ⁸ ≥ 0`, with `I_D = ∫ψ₂^{NF} = 10.385947` and `K = k(0)(1+δ)`.
* Result Q's upper end carries `+2 I_Q Λ² c₁t₄`, with `I_Q = ∫ψ₂^Q = 1.0160366` (`dgscan/bounds.py`).
* `Ω^Q_up` uses the same `Ĉ′_EE` and `Ĉ₃` objects (RESULT_Q_STATEMENT §1), so `c₁` is one number for both results.
* Normalisations:
  * `8πG = 256/c_T`, with `c_T` normalised to a free real scalar;
  * `√(t₄²+α_T²) = 240|ĝ₃|` (DICTIONARY_T4);
  * `t₄` is the d = 3 Hofman–Maldacena parameter, with `⟨E(n)⟩ ∝ 1 + t₄(|εnn|²/ε*ε − ¼)`. The d = 3 `t₂` structure vanishes identically (§3.1);
  * `α_T` is the parity-odd structure.

  No flavour-current coefficient enters.

### 1.2 Prescription (task requirement)

**This note's prescription = H7's.**
* **Lemma R.** The light side is `−ω|_light = (2πi/4)[𝒦_+]_{r^{−1}}` for `k = 2` and Regge spin `n = 2`. Here `𝒦` is the
  kernel-contracted full tree correlator: exchange in all channels plus contacts (H7 §3: "the residue of the full `𝒢^E`").
* **D3/E4.** That residue is a superposition of plane-wave pairs on `H₂`, mapped by the *same* leading-Regge plane-wave kernel as the
  heavy states (EE_FOUNDATIONS §2.2: D3 [P]).
* **E5.** Normalisation by the flat limit.

This prescription gives `−Ĉ′_EE|_E = 16πG/(ν²+9/4)`.

**Two other prescriptions now circulating in the project.**
1. **The pointwise / fixed-angle prescription** (`D_LR`, RESULT_Q_EINSTEIN_CHECK §3.3–3.5).
   * Every `J₀ ≤ k` term is annihilated pointwise, so `c₁ = 0` trivially.
   * But the Einstein pole is then also 0, so this is **not** Result D's prescription. Its values "live at the collinear set / ν
     integration", which is exactly where Lemma R's residue sits.
2. **The term-wise `s`-channel OPE prescription** (`Ω_OPE` in `code/resultq/lsd_mellin.py`; CLS-D's `ω|_light`).
   * L0 removes every contact and every `[TT]`. The light side is then the `s`-channel `T` block.
   * In the light-transform helicity basis, Einstein `⟨TTT⟩` is purely preserving and `t₄` purely flip (Q8, §3.1). So its `E × t₄`
     element is the helicity-off-diagonal action of `Ĉ′_EE` on the state `T`, and its leading-Regge image vanishes for `EE` (§3.2).
   * Whether the *full* action of `Ĉ′_EE` on the `T` block's off-diagonal structure vanishes is not established [I].
   * If `Ω_OPE` and Lemma R agree, which is what LIGHTSIDE_DIST is testing, that full action equals the `c₁` of §3.5.

All statements below are in H7's prescription unless marked otherwise.

---

## 2. Decomposition of the `O(t₄)` tree correlator at the `u`-Regge point

Write `𝒢 ⊃ t₄𝒢₁` (LIGHT_INPUTS §3). Seen from the `u`-Regge point, `𝒢₁` has two parts.
* **(U)** The `u`-channel `E × W³` graviton exchange, between the Regge lines `(1,3)` and `(2,4)`. It has a genuine `j = 2` Regge pole:
  pole part plus off-shell vertex part.
* **(C)** The `s`- and `t`-channel `E × W³` exchanges and the `O(ĝ₃)` quartic contact (the quartic expansion of `W³`). These are
  contact-like there: no `u`-channel pole, and polynomial in `ν²` (ADS4_LIGHT §0 item 4 uses the same reading for `W³⊗W³`).

**Regge spin.** `J₀(E×W³) = 2` [D-sketch].
* **Flat 4D.** The `O(ĝ₃)` four-graviton amplitudes are `M(−+++) ~ s²/t` and `M(++++) ∝ stu`. Both have Regge spin 2, and there is
  no `O(ĝ₃)` MHV amplitude [P: little group, dimension, factorisation].
* **AdS.** The Mellin-degree-3 part of the AdS amplitude is the flat amplitude. Every lower-degree part has Regge spin ≤ 2.
* Stage 12's "`≤ 3`" is the conservative count. Only `n = 2` is used below, as in H7.

**Helicity language.**
* "Preserved/flip on a line" means `a₃ = −a₁` (resp. `a₃ = a₁`) on the light-ray line `(1,3)`, and `a₄ = −a₂` (resp. `a₄ = a₂`) on the
  dual line `(2,4)`. These are the polarisations of the plane-wave kernel `L(a₁,a₃;a₂,a₄)(ν,β)` (OFFDIAG_KERNEL §1.2).
* In the bulk, a preserved line has one self-dual (`Ψ`) and one anti-self-dual (`Ψ̄`) linearised Weyl tensor. A flipped line has two
  of the same type.

---

## 3. Results

### 3.1 The vertex at `O(t₄)` is pure flip, on-shell and off-shell [Q + DV; D]

**(a) On-shell (pole).** The d = 3 HM/CDP one-point matrix in the helicity basis is `H = 𝟙 + (t₄/4)(e^{−4iφ}σ₊ + e^{4iφ}σ₋)`, with
`α_T` entering only off-diagonally (Q8; LIGHT_INPUTS `t2`).
* Re-checked here (`c1_symbolic.py (A)`):
  * `⟨±|t₄-structure|±⟩ = 0` and `⟨∓|t₄-structure|±⟩ = e^{∓4iφ}/4`;
  * the `t₂` structure vanishes identically in d = 3.
* The diagonal is the Ward-identity normalisation and does not depend on `t₄`. **So the residue of the `j = 2` pole has no preserved
  ⊗ preserved `O(t₄)` part.**

**(b) Off-shell (vertex polynomial).** ADS4_LIGHT §2.3 introduced `c₁` precisely as `t₄[V(ν) − V(−3i/2)]/(ν²+9/4)`. Two independent
arguments show that this vanishes in the preserved sector.

* **(b1) Self-dual / anti-self-dual decomposition [D].**
  * In 4D, `Tr W³ = Tr W₊³ + Tr W₋³` holds identically, for any metric.
  * The cubic `W³` vertex therefore needs both legs on a line to carry the same type of Weyl tensor (`Ψ,Ψ` or `Ψ̄,Ψ̄`). That is the flip
    configuration, whatever the off-shell exchanged leg is.
  * Linearised gravitons on AdS₄ (vanishing background Weyl) split into self-dual and anti-self-dual parts. In the eikonal limit this
    split coincides with the transverse helicity along the line [D-sketch].
  * On-shell this reproduces (a), and flat `W³` 3-point amplitudes are all-plus/all-minus.
* **(b2) Weyl factorisation [D + DV].**
  * `W³` couples to the exchanged (shock-type) field `h_uu = f(x)g(u)` only through that field's linearised Weyl tensor.
  * For the AdS₄ shock, `W_{uiuj} = −(g(u)/2z)·D_ij f` exactly, with `D_ij f = ∇_i∇_j f − ½g_ij□f` the traceless Hessian on `H₂` and
    `f = F/z`. Every other Weyl component vanishes, and `R_uu + 3g_uu = −½(g/z)(□−2)f` (`c1_symbolic.py (B)`, sympy, arbitrary `F`, `g`).
  * The preserved part of the vertex is a scalar, `SO(2,1)`-invariant differential operator on `H₂`. Such operators are polynomials
    `P(□)` (Harish-Chandra) [P].
  * `P(□)` must annihilate `ker D ⊇ {1, X^A}`, which have `□ = 0` and `2`. Both are checked in `(C2)`. Hence **`P(0) = P(2) = 0`**, and
    by `(C1)` the minimal form is `∇^i∇^jD_ij f = ½□(□−2)f`.
  * Then `α_pres(ν) = P(□)/(2−□)` evaluated at `□ = −(ν²+¼)` (eigenvalue checked in `(C3)`). For degree 2 this is
    `a(ν² + ¼)` (`(D)`): the constant is locked to the `ν²` coefficient.
  * The `ν²` coefficient is a flat-space quantity: the preserved `u`-channel `W³` vertex. It is zero by (b1) in flat space [P], so
    `a = 0`. Higher degrees are excluded by flat-limit finiteness.

  For comparison, a scalar `(∂φ)⁴`-type contact corresponds to `P = (□−2)·const` with `P(0) ≠ 0`. That is CHMRS's ν-independent
  `g₂`. It is possible only because the scalar vertex is not built from the Weyl tensor.

### 3.2 `EE` annihilates every flip plane-wave configuration [P / DV]
`c1_flipcheck.py` evaluates the exact leading-Regge plane-wave image `L(a₁,a₃;a₂,a₄)(ν,β)` of the de-shifted `k = 2` functional. It
uses the `EE` dual structure through `s1_offdiag_kernel.img_vec` and `A_monos_struct('EE')`, imported read-only. All 16 polarisation
assignments are run, at `ν ∈ {0.3, 1.5, 4, 12, 40}` and `β ∈ {0.05, 0.2, 0.5, 1, 1.8}`. The table gives `max|L|/|L_pres|`:

| structure | (1,3) flip × (2,4) preserved | (1,3) preserved × (2,4) flip | double flip |
|---|---|---|---|
| **`EE`** | `≤ 3.9·10⁻¹⁴` | `≤ 4.9·10⁻³²` | `≤ 4.8·10⁻⁴⁶` |
| `q²` (control) | `≤ 1.9·10⁻¹⁴` | **`0.17 – 3.9`** | `≤ 5·10⁻¹⁵` |
| `Ĉ₃` (`o_H q`) | `≤ 9·10⁻¹⁶` | `≤ 10⁻¹⁵` | `≤ 3·10⁻³²` |

* The dual-line zero for `EE` is exact: `(z₂′·z₄′)² = (ε_h·ε_h)² = 0` (EE_FOUNDATIONS §2.1) [P].
* The light-ray-line zero is [DV] for every `j′ = 0` structure; the `(1,3)` side is common to all of them.
* **The `q²` control shows that REVIEW_EE m5's reasoning would have failed for the `q²` functional** (its `u`-channel `c₁` part would
  be `O(1)`). It holds for `EE`.

### 3.3 `k = 2`, `u`-channel part: `c₁^{(U)} = 0` [D]
By Lemma R (`n = 2`), D3/E4 and the eikonal factorisation of the `u`-channel exchange into shock vertices (CHP, Q3 [Q]), the `O(t₄)` part
of (U) is a superposition of plane-wave pairs.
* Its polarisation content is the vertex product `[𝟙 + t₄·flip]_{(1,3)} ⊗ [𝟙 + t₄·flip]_{(2,4)}` at `O(t₄)`.
* By §3.1 this is flip ⊗ preserved plus preserved ⊗ flip, for the pole part (a) and the off-shell part (b) alike.
* By §3.2 the `EE` kernel maps each of these configurations to zero pointwise in `(ν, β)`.

**Grade [D].** The inputs are H7's: Lemma R [P given R1–R5], D3 [P], Q3/Q8 [Q], plus §3.1 [D/DV] and §3.2 [P/DV]. This makes
REVIEW_EE m5's one-line remark precise, and extends it to the off-shell vertex, which m5 did not address.

### 3.4 `k = 3`: `c₁^{(3)} = 0` [P given Q2, Q8, F]
* **Route A** (H7 §4). The whole `O(1/c_T)` light side of `Ĉ₃` is `∝ [H(n₁), H(n₃)]`. With `H = 𝟙 + c·flip`, this is
  `−2i|c|² sin4φ₁₃ σ_z` (LIGHT_INPUTS `t2` [DV]): it does not depend on the identity part, so no term is linear in `t₄` or `α_T`.
* **Route B.** Lemma R with `n = 2` plus "`Ĉ₃` has no `Π₂` moment" (H7 §4, §15 [D + NE]). This kills *every* Regge-spin-2 `O(t₄)`
  term, including the contact-like part (C), in H7's prescription.
* Consistently, the `Ĉ₃` plane-wave image of every flip configuration vanishes (§3.2).

### 3.5 What is left: the contact-like part (C) is exactly `ν`-independent, and its coefficient is `c₁` [D]; value [I]
* **Polynomial.** (C) has no `u`-channel pole. A pole `∝ t₄/(ν²+9/4)` is excluded by §3.1(a). Other non-polynomial residues are
  removed by the `γ²` normalisation (H7 §5.1 step 3), up to ADS4_LIGHT's input (E). So
  `−Ĉ′_EE|_{(C)} = 16πG t₄ Σ_n c_{1,n} ν^{2n}` in the preserved ⊗ preserved sector; flip sectors are killed by §3.2.
* **Flat-limit power counting.** `t₄ ∝ ĝ₃/L⁴`. A term `t₄ν^{2n}` relative to the Einstein `ν⁻²` scales as `ĝ₃ p^{2n+2}L^{2n−2}`.
  * `n ≥ 2` diverges in the flat limit, so it is excluded.
  * `n = 1` is the flat `O(ĝ₃)` preserved ⊗ preserved term of the flip-free flat functional (`EE ≡ q²` at the crystal, flip-free
    [P given H0 §17 map]). It vanishes because no `O(ĝ₃)` four-graviton MHV amplitude exists in 4D [P].
  * `n = 0` survives with relative size `Λ⁻²`.
* **Hence `−Ĉ′_EE ⊃ 16πG c₁t₄` exactly, with `c₁ := c_{1,0}`** [D given J₀ = 2 (D-sketch) and (E)]. It is one pure number:
  * universal (it depends on nothing but the GR + `W³` tree structure);
  * independent of `Δ_g` and `ν`;
  * it is the `L = 0` (zero impact parameter) remnant, at order `1/L²`, of a flat quantity that vanishes.

  The same argument, with `α_T` in place of `t₄`, gives a possible `α_T` constant. It is parity-odd, so it is absent from the
  parity-even `Re` functional for parity-invariant theories (ADS4_LIGHT §0 item 3).
* **Why §3.1 does not force it to zero.**
  * In the preserved ⊗ preserved sector, (C) is fed by the `s`- and `t`-channel exchanges. For example, `W₊(1,2,int)` with lines
    `(1,3)` and `(2,4)` both preserved has legs `1, 2` of type `Ψ` and `3, 4` of type `Ψ̄`.
  * It is also fed by the quartic seagull `W₊^{(1)}[1]W₊^{(1)}[2]W₊^{(2)}[3,4]`.
  * In flat space these sum to zero (no MHV at `O(ĝ₃)`), but nothing shown here forbids a `1/L²` remnant at Regge spin 2. In 3d CFT
    language, this is the analogue of the non-flat `+++` components of Einstein `⟨TTT⟩` (Maldacena–Pimentel).
  * The premise that would make it vanish ("the `O(ĝ₃)` leading-Regge residue factorises through shock vertices; no seagull at
    `L = 0`") is **not derived**. I do not claim it.
* **Natural size [M].**
  * The `|ĝ₃|²` light term `2πG|ĝ₃|²ν⁶ = 16πG(c_vĝ₃ν⁴)²/ν²` fixes the `W³` vertex normalisation at `c_v = 1/√8 ≈ 0.35` (ADS4_LIGHT §0
    item 4).
  * A curvature remnant of one vertex gives `c₁t₄ ~ c_vĝ₃·O(1)`, i.e. `c₁ ~ O(1)·c_v/240 ≈ 1.5·10⁻³·O(1)`.
  * With `O(1)` read generously as up to 30–60, this gives `|c₁| ≲ 0.1`.
  * This is an estimate, not a bound.

### 3.6 Bound on `|c₁|` (task item c)
* **No rigorous or derived bound is available [I].** `c₁` enters every sum rule in hand only as the product `c₁t₄`. Every CFT with
  `t₄ = 0` satisfies all of them for any `c₁`.
* What is derived is weaker but useful:
  * the *form*: a universal constant, no `ν`- or `Δ_g`-dependence [D];
  * its flat counterpart vanishes [P];
  * the `u`-channel and `k = 3` contributions are zero [D] / [P given …].
* What closes it is a computation, not a bound: the `u`-Regge `r^{−1}` coefficient, preserved ⊗ preserved and constant part, of the
  AdS₄ `s`- and `t`-channel `E × W³` exchange diagrams plus the `O(ĝ₃)` contact. Equivalently, in `Ω_OPE`, the full action of `Ĉ′_EE`
  on the `T` block's `E × t₄` structure (§7).

---

## 4. Exact quadratic-root bounds for every `c₁` (`code/c1/c1_quadratic.py` → `cache/c1_quadratic.txt`)

### 4.1 Result D [D]
Let `R = √(t₄²+α_T²)` and `K = k(0)(1+δ)` (headline `K = 1`). The sum rule
`K R²Λ⁸/57600 − 16I_DΛ²c₁t₄ − 16L_A ≤ 0`, with `t₄` and `α_T` free, gives:

* **`R ≤ (57600/(KΛ⁸))·[8I_DΛ²|c₁| + √(64I_D²Λ⁴c₁² + KΛ⁸·16L_A/57600)] = R₀·(x + √(1+x²))`**, where
  * `R₀ = 240√(16L_A/K)/Λ⁴`;
  * `x = 1920 I_D|c₁|/(Λ²√(K·16L_A)) = 19941.0·|c₁|/(Λ²√(K·16L_A))`;
  * `x/|c₁| = 2.16006·10⁻³` (300) and `1.67026·10⁻⁴` (1000).
* **`t₄` interval (asymmetric):** for `c₁ > 0`, `−R₀(√(1+x²) − x) ≤ t₄ ≤ R₀(√(1+x²) + x)`; for `c₁ < 0`, mirror it.
* **Cubic closure:** `|ĝ₃|²Λ⁸ ≤ (16L_A/K)(x + √(1+x²))²`.
* The linearised form `1 ± x` is the first order of this. The relative error of the linearisation is `√(1+x²) − 1 ≈ x²/2`, which
  reaches 1 % at `x = 0.152` and 10 % at `x = 0.585`.

| `Δ_g` | `|c₁|` | `x` | **`√(t₄²+α_T²) ≤`** (exact) | ratio to `c₁ = 0` | linearised `(1+x)R₀` | `t₄` range (`c₁ > 0`) | `|ĝ₃|²Λ⁸ ≤` |
|---|---|---|---|---|---|---|---|
| 300 | 0 | 0 | **3.1321·10⁻⁶** | 1 | 3.1321·10⁻⁶ | `±3.132·10⁻⁶` | 10734.6 |
| 300 | 1 | 0.00216 | 3.1388·10⁻⁶ | 1.00216 | 3.1388·10⁻⁶ | `[−3.125, +3.139]·10⁻⁶` | 10781.1 |
| 300 | 10 | 0.0216 | 3.2004·10⁻⁶ | 1.0218 | 3.1997·10⁻⁶ | `[−3.065, +3.200]·10⁻⁶` | 11208.5 |
| 300 | 100 | 0.216 | 3.8808·10⁻⁶ | 1.2391 | 3.8086·10⁻⁶ | `[−2.528, +3.881]·10⁻⁶` | 16480.8 |
| 300 | 1000 | 2.16 | 1.4221·10⁻⁵ | 4.5404 | 9.8975·10⁻⁶ | `[−0.690, +14.22]·10⁻⁶` | 221294 |
| 1000 | 0 | 0 | **2.8913·10⁻⁸** | 1 | 2.8913·10⁻⁸ | `±2.891·10⁻⁸` | 14339.6 |
| 1000 | 1 | 1.67·10⁻⁴ | 2.8917·10⁻⁸ | 1.00017 | 2.8917·10⁻⁸ | `[−2.891, +2.892]·10⁻⁸` | 14344.4 |
| 1000 | 10 | 0.00167 | 2.8961·10⁻⁸ | 1.00167 | 2.8961·10⁻⁸ | `[−2.886, +2.896]·10⁻⁸` | 14387.6 |
| 1000 | 100 | 0.0167 | 2.9400·10⁻⁸ | 1.0168 | 2.9396·10⁻⁸ | `[−2.843, +2.940]·10⁻⁸` | 14826.7 |
| 1000 | 1000 | 0.167 | 3.4142·10⁻⁸ | 1.1809 | 3.3742·10⁻⁸ | `[−2.448, +3.414]·10⁻⁸` | 19996.2 |

Notes on the table:
* The bound on `R` is even in `c₁`. Negative `c₁` mirrors the `t₄` interval.
* The linearisation is good to 1 % for `|c₁| ≤ 70.5` (300) and `|c₁| ≤ 911` (1000), and to 10 % for `|c₁| ≤ 271` / `3502`.
* `x = 1` (bound × `1+√2`) is reached at `|c₁| = 463` / `5987`. This matches the referee's "≈ 460".
* With `K = 1 + δk(0)`, `R₀ = 3.13326·10⁻⁶` / `2.89136·10⁻⁸`.

### 4.2 Result Q [D given `κ > 0`]
The upper end is `U(t₄, α_T) = 2L_AdS − κ(t₄²+α_T²)Λ⁸/57600 + 2c₁t₄Λ²I_Q`, with `κ = 0.109908` (flat). Taking the supremum over the
unknown `(t₄, α_T)`:

**`g̃₄Λ⁶/8πG ≤ 2L_AdS + 57600 I_Q² c₁²/(κΛ⁴)`**, attained at `t₄* = 57600 I_Q c₁/(κΛ⁶)` and `α_T = 0`.

* The coefficient `57600 I_Q²/κ = 541017`. Dividing by `Λ⁴` gives **`6.8145·10⁻⁵`** (300) and **`5.4428·10⁻⁷`** (1000).
* Imposing Result D's exact `t₄` interval changes nothing for any `|c₁| ≤ 1000`. For example, at `c₁ = 1000`, Δ_g = 300,
  `t₄* = 7.5·10⁻⁷` lies inside D's `1.4·10⁻⁵`.

| `|c₁|` | `ε_{c₁}` (exact, 300) | upper end (300) | old linear (300) | `ε_{c₁}` (exact, 1000) | upper end (1000) | old linear (1000) |
|---|---|---|---|---|---|---|
| 0 | 0 | **72.041** | 0 | 0 | **92.814** | 0 |
| 1 | 6.8·10⁻⁵ | 72.041 | 0.567 | 5.4·10⁻⁷ | 92.814 | 0.0586 |
| 10 | 0.0068 | 72.047 | 5.67 | 5.4·10⁻⁵ | 92.814 | 0.586 |
| 100 | 0.68 | 72.72 (+0.9 %) | 56.7 | 0.0054 | 92.819 | 5.86 |
| 1000 | 68.1 | 140.2 (×1.95) | 567 | 0.544 | 93.36 (+0.6 %) | 58.6 |

**Record correction.** UNIFIED_THEOREM §2, DG_SCAN, RESULT_Q_STATEMENT §0 and `prl/letter.tex` quote `|ε_{c₁}| ≤ 0.57|c₁|` /
`0.059|c₁|`.
* That figure bounds `2|c₁t₄|Λ²I_Q` alone by importing D's `t₄`. It drops the `−κ|ĝ₃|²Λ⁸` term, which is quadratic and dominates.
* **The exact `c₁`-dependence of Q is `+6.8·10⁻⁵c₁²` / `+5.4·10⁻⁷c₁²`, and it needs no input from Result D.** So RU issue 14's
  "importing D's bound makes Q inherit `H_LJ`, `H_M`, `H_R`, `H_edge`, `H_RS`" disappears.
* The only input is `κ > 0`, which Q already needs for general `ĝ₃` ([D-sketch + M]; κ < 0 breaks Q anyway, RU issue 7).
* If `κ` were smaller than its flat value, the coefficient scales as `1/κ`. For example, `κ = 0.01` gives `7.5·10⁻⁴c₁²` at 300.

---

## 5. Consequences for the headline bounds

1. **Cubic bound (Result D, item 1 of UNIFIED_THEOREM).** Replace "`(1 + O(c₁Λ⁻²))`, `c₁` free" with:
   * "`(x + √(1+x²))`, `x = 19941|c₁|/(Λ²√(16L_A))`, with `c₁` a single universal constant whose `u`-channel (`T`-exchange) and `k = 3`
     parts vanish [D]; its remaining contact-like part is uncomputed [I], natural size `≲ 0.1` [M]".
   * At the natural size the effect is `≤ 0.02 %` (300) and `≤ 0.002 %` (1000). That is 10× below the H17 band.
   * The bound is finite for every `c₁`. Even `|c₁| = 100` (≈ 10³ × natural) moves it by 24 % (300) / 1.7 % (1000).
2. **Quartic bound (Result Q, item 2).** Replace `|ε_{c₁}| = 2.03|c₁t₄|Λ² ≤ 0.57|c₁|` / `0.059|c₁|` with
   `ε_{c₁} ≤ 6.8·10⁻⁵c₁²` / `5.4·10⁻⁷c₁²` [D given `κ > 0`].
   * **Item 2 then no longer depends on item 1 at all.** This removes UNIFIED_THEOREM §2 Remark 1's only link, and RU issue 14.
   * For the natural size the effect is `< 10⁻⁶` of `2L_AdS`.
3. **Referee M5.**
   * Problem 1 ("`c₁` is a definite number, calling it free hides that it is uncomputed"): it is now *one* definite universal number,
     and the parts the referee named are computed to be zero. The uncomputed remainder is localised to one well-defined AdS₄
     Witten-diagram coefficient.
   * Problem 2 ("unbounded, linearisation valid only for `|c₁| ≪ 460`"): replaced by the exact forms of §4, valid for all `c₁`.
   * Suggested wording: "uncomputed universal constant", not "free parameter", in the theorem and in Table I.

---

## 6. Grade table

| statement | grade |
|---|---|
| prescription = H7's Lemma-R residue + D3 + E5 | definition |
| d = 3 one-point matrix: `t₄`, `α_T` purely flip; `t₂` absent; diagonal `t₄`-free | [Q] + [DV] (sympy, `c1_symbolic (A)`) |
| `W³` cubic vertex is flip-only off-shell (self-dual/anti-self-dual decomposition) | [D]; SD/ASD ↔ eikonal helicity [D-sketch] |
| shock Weyl = traceless `H₂` Hessian; `R_uu ∝ (□−2)f`; `∇^i∇^jD_ij = ½□(□−2)`; `ker D ⊇ {1, X^A}` | [DV] (sympy, exact) |
| preserved vertex `P(□)` with `P(0) = P(2) = 0` ⇒ `α_pres = a(ν²+¼)`; flat ⇒ `a = 0` | [P] (algebra) + [D] (flat input) |
| `EE` annihilates dual-line flips | [P] (+ DV `≤ 5·10⁻³²`) |
| `EE` annihilates light-ray-line flips and double flips | [DV] (`≤ 4·10⁻¹⁴`, `≤ 5·10⁻⁴⁶`; `q²` control `O(1)`) |
| `c₁^{(U)}` (`u`-channel `T` exchange, `k = 2`) `= 0` | **[D]** |
| `c₁^{(3)}` (`k = 3`) `= 0` | **[P given Q2, Q8, F]**; route B [P given Lemma R, `n = 2`, no-`Π₂`] |
| `J₀(E×W³) = 2` | [D-sketch] (flat [P]; Mellin-degree argument) |
| `O(t₄)` light term of `Ĉ′_EE` `= 16πG c₁t₄` exactly (ν-independent, universal) | [D given J₀ = 2, (E)] |
| value of `c₁` | **[I]** |
| `|c₁| ~ 10⁻³–10⁻¹` | [M] (natural size; not a bound) |
| rigorous bound on `|c₁|` | none [I] |
| exact D bound `R₀(x+√(1+x²))`, `t₄` interval, cubic closure | [D] (given the `c₁ = 0` implication's own inputs) |
| exact Q correction `57600I_Q²c₁²/(κΛ⁴)`, D-free | [D given `κ > 0`] |
| Q's old `0.57|c₁|` / `0.059|c₁|` | superseded (overestimate; drops the `−κ|ĝ₃|²` term) |

---

## 7. What remains

1. **The value of `c₁` [I].** It is the constant term of the preserved ⊗ preserved `u`-Regge `r^{−1}` coefficient of the AdS₄
   `s`- and `t`-channel `E × W³` exchange diagrams plus the `O(ĝ₃)` quartic contact, normalised as in H7 (§3.5). Two routes:
   * a spinning Mellin / conformal-Regge evaluation of those diagrams. The Regge-spin-2 part only is needed, and its flat part is
     known to vanish, so only the `1/L²` remnant must be computed;
   * the full action of `Ĉ′_EE` on the `T` block's `E × t₄` structure (`Ω_OPE` route; ties to LIGHTSIDE_DIST).
   * A proof that the `O(ĝ₃)` leading-Regge residue factorises through shock vertices (no zero-impact-parameter seagull) would give
     `c₁ = 0` at once, by §3.1(b).
2. **Prescription consistency.** H7's Lemma R and `Ω_OPE` assign Regge-spin-`k` contact pieces differently: CHMRS's scalar `g₂`
   enters `C₂` in Lemma R, while L0 drops contacts. `c₁` is exactly such a piece. It vanishes in `D_LR`. In Lemma R it equals the `(C)`
   constant; in `Ω_OPE` it equals the `T`-block off-diagonal action. LIGHTSIDE_DIST's outcome decides which one Result D needs.
3. Inherited inputs: (E) (no non-polynomial residues); `J₀(E×W³) = 2` written out for spinning Mellin amplitudes; H7's D2, D3 and
   §18.1 kernel; the SD/ASD ↔ eikonal-helicity identification.

---

## 8. Reproduce

Every command runs on one core with `ulimit -v 3000000; export OMP_NUM_THREADS=1`. Each writes to `ads4/code/c1/cache/`.

```
cd ads4/code/c1
python3 c1_flipcheck.py    # §3.2: EE / q2 / C3 plane-wave flip images, 16 configs  -> cache/c1_flipcheck.txt   (~1 min)
python3 c1_symbolic.py     # §3.1: HM d=3 matrix, AdS4 shock Weyl, H2 identities     -> cache/c1_symbolic.txt    (~1 min)
python3 c1_quadratic.py    # §4:   exact D and Q bounds vs c1 at Delta_g = 300, 1000 -> cache/c1_quadratic.txt   (seconds)
```

The scripts import `spinfunc/s1_offdiag_kernel.py`, `h7i14/v3_regge_kernel.py`, `q1_c3_regge_kernel.py`, `dgscan/bounds.py` and
`omegau/nf_common.py` read-only. They write no bytecode, and no other file is touched.
