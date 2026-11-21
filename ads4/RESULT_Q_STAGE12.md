# Result Q, Stages 1–2: the CFT₃ lift of the forward pieces, and the exact light-side gate

*21 Nov 2026, Stage-1/2 agent. Unrefereed, not committed. No existing document or code was edited. New code:
`ads4/code/resultq/stage12_*.py` (outputs in `resultq/cache/stage12_*.txt`, `.log`). Downloads (data only): scratchpad
`…/scratchpad/resultq12/` (arXiv abstract pages of 1201.6452, 1712.02800, 1804.00949). The CLS-D, CHMRS and BCP LaTeX sources already in
the scratchpad (`dict/src/2311.04271`, `chmrs_src`, `resultq/bcp`) were read as data. Total compute: about 15 CPU-minutes on 2 cores.*

Grades as in `RESULT_D_STATEMENT.md`: **[certified] [P] [P given X] [D] [D-sketch] [Q] [DV] [NE] [M] [I]**. A chain has the grade of its
weakest link. "[Q, abstract]" means only the abstract was read. "Unverified" marks literature quoted from memory.

Conventions: Stage 0 (`RESULT_Q_STAGE0.md`). Flat units `M = 1`, `u = M²/m²`, `𝒥 = J(J+1)`. AdS units `R = 1`, `Λ = MR`, `ν = Λp`.
Dictionary for the scalings: `G_AdS = G_f/Λ²`, `|ĝ₃|²_AdS = |ĝ₃|²_f/Λ⁸`, `g₄,AdS = g₄,f/Λ⁸`, `x² := |ĝ₃|²_fM⁸`. A "flat-normalised" AdS
functional `Ĉ_k` is normalised so that `Λ^{2k}Ĉ_k(Λp) → B_k(p)` on heavy blocks.

---

## 0. Answer in brief

1. **Stage 1, the lift: possible; no obstruction found.** The forward pieces lift to CLS-D subtracted light-ray functionals on `⟨TTTT⟩`.
   * `B₄⁽¹⁾|₀`, `∂_{p²}B₄⁽¹⁾|₀` ← `Ĉ₄(ν)`: `j′ = 0`, CLS-D's `f⁺_{4;2} = ((v′−u′)³ − (v′−u′))/(4u′v′)`, Result D's dual structure `EE`.
     Take `ν = 0` and `Λ²∂_{ν²}` at `ν = 0`.
   * `B₅⁽¹⁾|₀` ← `Ĉ₅(0)`: `j′ = 1`, `f⁻_{5;2} = −((u′−v′)² − 1)/(4u′v′)`, the structure `o_H·q` of Result D's `Ĉ₃`.
   * Smeared `B₂` ← Result D's `Ĉ′_EE`; optional `B₃` ← `Ĉ₃`.

   Properties:
   * **Flat images.** At the scattering crystal, `f_4/f_2 = f_5/f_3 = u²/(s(s+u))` exactly [P, sympy]. These are the two extra flat
     subtractions of `B₄` vs `B₂`. The structures are unchanged, so `EE ≡ q²` (flip-free, `q²` stripped density) carries over
     [P given CLS-D crystal (Q) + H0 §17 map].
   * **Regge convergence.** Both need `J₀ < 4` (CLS-D: `k > J₀` for even `j′`, `k−1 > J₀` for odd `j′`) [Q]. The margin over KKSDZ's
     `J₀ ≤ 1` is 3 units, against 1 unit for `k = 2`.
   * **Regge kernels.** `f⁺₄ = f⁺₂·16η²r²(1+O(r))` and `f⁻₅ = f⁻₃·16η²r²(1+O(r))` [P]. So the leading moments are `Π₄`, `Π₅`, and the
     leading-Regge helicity off-diagonal of the `EE` functional stays zero [D].
   * **Not pure harmonics.** The leading kernels carry `η³` (`k = 4`) and `η²` (`k = 5`). The heavy images mix `P_{l±1}, P_{l±3}` (resp.
     `P_l, P_{l±2}`), so a pure-moment ("de-shifted") forward functional needs a 4-shift (3-shift) Abel inversion. That inversion is not
     constructed [I].
     * Scalar proxy: the raw `ν = 0` image is positive and grows like `P_{5/2}(cosh β) ~ (2J/m)⁵` [DV for the scalar proxy; `EE`: I].
   * **The `2𝒥` growth of `∂B₄|₀` survives curvature, with the right sign** in the pure-moment model.
     * `−∂_{ν²}P_{−½+iν}(cosh β)|₀ = (√2/2π)∫₀^β t²(cosh β − cosh t)^{−½}dt > 0` for every `β` [P]. Its ratio to `P_{−½}` equals `𝒥u/Λ²`
       for `β ≪ 1`, and keeps growing like `≈ (0.75–1)β²/4` (`log²J`) for `β ≳ 1` [DV].
     * The Stage-0 functionals are positive in this leading-Regge model on every physical state with `b̂ ≥ 20`, at `Λ = 100, 300, 1000`
       (min `+0.27`) [M].
2. **Stage 2, the exact tree-level light side: computed (analytically), with the dangerous terms exactly zero.**
   * **Lemma A** (§2.1): `Ψ₄`, `Ψ₅` annihilate, for every `ν` (so also `∂_ν`, any smearing and any de-shift):
     * the tree-level Einstein correlator (`T` exchange in all channels, plus the EH contact; `J₀ = 2`);
     * the linear `E×W³` correlators (`t₄`, `α_T`; `J₀ ≤ 3`);
     * the `O(c_T⁰)` MFT correlator.

     **Their light sides are identically zero.** This is a pointwise commutator identity, not a polynomial-structure argument, so it also
     excludes the decaying residues `G(ν²+9/4)^{−n}` that carried Stage 0's leverage.
     Grade: [P given V₄ (CLS-D validity at `J₀ < 4`, applied to tree correlators), microcausality of tree correlators (D), and the `J₀`
     counts (Q/D)] → **[D]**.
   * Exact light sides (tree level, flat-normalised; `b₀, b₅` are `O(1)` AdS numbers that are not computed):

     | bulk term | `−Ĉ′_EE` (smeared `k = 2`) | `−Ĉ₄(0)` | `−Λ²∂_{ν²}Ĉ₄(0)` | `−Ĉ₅(0)` |
     |---|---|---|---|---|
     | Einstein (`T` in all channels + EH contact) | `16πG/(ν²+9/4)` | **0** | **0** | **0** |
     | `E×W³` (linear `t₄`, `α_T`) | `16πGc₁t₄` (`c₁` [I], Result D) | **0** | **0** | **0** |
     | `W³×W³` (`∝ |ĝ₃|²`) | `2πG|ĝ₃|²(ν⁶+…)` | `b₀G|ĝ₃|²` [I] | `4πG|ĝ₃|²` | 0 |
     | `R⁴` MHV (`g₄`) | 0 | **`2g₄`** (exact, `ν`-independent) | **0** | **0** |
     | `R⁴` all-plus (`ĝ₄`) | 0 | 0 | 0 | 0 |
     | `D²R⁴` (`g₅`) | 0 | `b₅g₅` (curvature, `Λ⁻²`) | `g₅` | `g₅` |
     | higher contacts | 0 | `O(Λ^{−2j})` | `O(Λ^{−2j})` | `O(Λ^{−2j})` |

     Zeros: Lemma A. The `R⁴`, `|ĝ₃|²`-top and `g₅` entries: [P given (F) + (P)] (§2.2). The `k = 2` column is Result D.
   * Not done: no explicit `⟨TTTT⟩` Witten-diagram integral was evaluated; none exists in a usable form here (Raju 1201.6452 is
     momentum space) [I]. The values follow from Lemma A plus flat-limit identities. The cross-checks are the scalar CHMRS results (read
     from source), Raju's abstract, and an exact check of the BCP/Penedones chain (§3).
3. **Gate: residual leverage ≤ 10² ⇒ GO for Stage 3**, conditional on Lemma A's inputs (§2.4).
   * **Before:** an unknown unit `O(G)` residue would shift the bound on `g₄Λ⁶/(8πG)` (signal 35–80) by
     `3.9·10¹⁰ / 2.8·10¹³ / 3.9·10¹⁶` (`Ĉ₄(0)`) and `1.4·10¹⁴ / 9·10¹⁷ / 1.4·10²²` (`∂Ĉ₄`), at `Λ = 100 / 300 / 1000`.
   * **After:** those coefficients are exactly 0. What is left, per unit of the uncomputed `O(1)` numbers:
     * `b₀`: `2·10⁻⁴` (`x² ≤ 110`) to `1.5·10⁻²` (Result D's `x` bound) at `Λ = 100`;
     * `b₅`: `5·10⁻⁵`;
     * Result D's open `c₁`: `0.35` to `2.9`.
   * **Residual leverage ≤ 3 ≪ 10².**
   * **Bonus:** the exact pole makes the AdS IR log exact. `2L_AdS = 2∫ψ₂/(p² + 9/(4Λ²))` replaces Stage 0's [M] `A log(2Λ/3) + B`, and
     differs from it by `+0.17 / +0.04 / +0.009` [P given the pole].
   * **The no-go branch of the gate is not triggered.** The gate is passed analytically, not by a numerical Witten-diagram evaluation.
     If Lemma A's input V₄ failed for the subtracted functionals, the leverage would be back to `10¹⁰–10²²`. V₄ is the single
     load-bearing input (§4).
4. **ABJM normalisation (§3).** Exact chain [DV given Q]: `g₄Λ⁶/(8πG) = 2ζ(3)` at leading order in `1/λ`, with
   `Λ = 2^{5/4}π^{1/2}λ^{1/4}`.
   * Penedones' formula maps a degree-`n` Mellin term with `(L²/2)ⁿ/Γ(n+½)` (`h = ½` for `S`). This reproduces BCP's
     `B⁴₄/B¹_SG = 105ζ(3)/(64√2π³)` exactly.
   * **CFT recipe:** with `ρ_M := λ^{−3/2}B⁴₄/B¹_SG`, the ratio of `R⁴` to SUGRA top-degree Mellin coefficients in `⟨SSSS⟩`,
     **`g₄/(8πG) = ρ_M/105`** in AdS units.
     * General component: `g₄/(8πG) = ρ_M(½)³Γ(h+n_SG)/Γ(h+n_SG+3)`.
   * In our functionals no Mellin convention enters. `−½Ĉ₄(0)|_{R⁴}` and the pole coefficient `16πG` are both heavy-normalised, so
     `−Ĉ₄(0)|_{R⁴}/(16πG) = g₄/(8πG)` exactly [P given (F), (P)].
   * Placement: ABJM sits at **6.8 % / 5.0 % / 3.9 %** of the exact-pole bound (robust `S_min + B₃`) at `Λ = 100 / 300 / 1000`
     (`λ ≈ 3.2·10⁵ / 2.6·10⁷ / 3.2·10⁹`).

---

## 1. Stage 1: the lift

### 1.1 Definitions
CLS-D's spinning functional (arXiv:2311.04271, eq. `spinningfunctional_definition`) is
`Ψ^{(a),(c)}_{k,ν,j′,f}[𝒢] = Ψ⁺[f·𝒢] − (1↔3)`, with subtraction factors `f^±_{k;J_e}` (eqs. `subtraction_factor_evenjprime/oddjprime`, read
from source). For `J_e = 2`:

| `k` | `j′` | `f` | type | used for |
|---|---|---|---|---|
| 2 | 0 | `v′ − u′` | unsubtracted | Result D's `Ĉ′_EE` (smeared `B₂`) |
| 3 | 1 | `−1` | unsubtracted | Result D's `Ĉ₃` (smeared `B₃`, optional) |
| **4** | **0** | **`((v′−u′)³ − (v′−u′))/(4u′v′)`** | subtracted | **`Ĉ₄`: `B₄|₀`, `∂_{p²}B₄|₀`** |
| **5** | **1** | **`−((u′−v′)² − 1)/(4u′v′)`** | subtracted | **`Ĉ₅`: `B₅|₀`** |

Choices:
* `Ĉ₄` uses `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²` and the unique `j′ = 0` celestial structure.
* `Ĉ₅` uses the `(1,3)`/`(2,4)` structures of `Ĉ₃` (`o_H·q`).
* Both are flat-normalised by their forward heavy normalisation (`N̂₄(ν)`, `N̂₅(ν)`, as `N̂₃` for `Ĉ₃`).
* Forward pieces: `Ĉ₄(0)`, `Λ²∂_{ν²}Ĉ₄|₀`, `Ĉ₅(0)`.
* The AdS version of Stage 0's family is
  `Ω^Q := ∫ψ₂Λ⁴Ĉ′_EE + ∫ψ₃½Λ⁶Ĉ₃ − a(Λ¹⁰∂_{ν²}Ĉ₄|₀ − Λ¹⁰Ĉ₅(0)) − cΛ⁸Ĉ₄(0)`, with Stage 0's weights.

### 1.2 Regge structure and convergence (`stage12_kernel.py`, `cache/stage12_kernel.txt`)
* **Frame.** CHMRS frame `x₁ = −y`, `x₃ = y`, `x₂ = −x₄ = x/x²`, with `r = |x||y|`, `η = −x·y/(|x||y|)`. Exactly:
  `u′ = (1 − 2rη + r²)²/(16r²)`, `v′ = (1 + 2rη + r²)²/(16r²)`. Hence `v′ − u′ = η(1+r²)/(2r)` and `u′v′` is a perfect square
  [P; numerical check on 2000 random configurations: `5·10⁻¹³`].
* **Laurent leading terms** [P, sympy]:
  * `f⁺₂ = η/(2r)`, `f⁻₃ = −1`;
  * **`f⁺₄ = 8η³r(1+O(r))`, `f⁻₅ = −16η²r²(1+O(r))`**;
  * `f⁺₆ = 128ηr³`, `f⁻₇ = −256r⁴`, `f⁺₈ = 2048η³r⁵`.

  So `f₄/f₂ = f₅/f₃ = 16η²r²(1+O(r))`. The `k = 4, 5` kernels are the `k = 2, 3` kernels with two more powers of `r`, i.e. leading Regge
  moments `Π₄`, `Π₅` (given `Π₂`, `Π₃` for `k = 2, 3`; for `k = 3` see H7 §4's route-B caveat, which does not affect the gate, §2.2).
  The `η`-degrees `1, 0, 3, 2, 1, 0` for `k = 2…7` agree with `SUPERCONVERGENT_UPLIFT.md` item 1.
* **Convergence.**
  * CLS-D: `f⁺` raises the Regge spin of the kernel by `k − 2J_e + 1`, `f⁻` by `k − 2J_e`. The integrals converge for `J₀ < k` (even
    `j′`) and `J₀ < k − 1` (odd `j′`) [Q]. Here that is `J₀ < 4` for both.
  * Non-perturbative input: KKSDZ `J₀, J₀^{dDisc} ≤ 1` [Q, with their loophole]. **Margin: 3 units** (Result D's `k = 2` had 1).
  * The new feature compared with `k = 2, 3` is the factor `1/(u′v′)`. It has poles on `x₁₂² = 0`, `x₃₄² = 0`, `x₁₄² = 0`, `x₂₃² = 0`,
    where the commutator support begins. These are the "finitely many non-vanishing double traces" of subtracted sum rules. Validity of
    the sum rule with such factors is CLS-D's statement (as for CHMRS's scalar `C_{k≥4}`) [Q]. The endpoint and edge bookkeeping
    (Result D's K2) is not redone for `k = 4, 5` [D-sketch, as K2].
* **Helicity structure.**
  * `f` depends only on the points, not on the polarisations. So the `k = 4` (`k = 5`) kernel has exactly the tensor structure of the
    `k = 2` (`k = 3`) kernel.
  * `EE`'s `(2,4)` side at leading Regge is 16 × the scalar one (EE_FOUNDATIONS §2.1).
  * Hence the **leading-Regge helicity off-diagonal of `Ĉ₄[EE]` vanishes**, as for `Ĉ′_EE` [D]. `Ĉ₅` inherits `Ĉ₃`'s `σ_h`-diagonal
    form [D].

### 1.3 Flat images (scattering crystal)
* CLS-D's crystal gives `u′ = m⁴/ν⁴`, `v′ = (m²−ν²)²/ν⁴`. The values of `f_k` for `k = 2…7` agree with CLS-D's crystal formulas
  [P, sympy]. (CLS-D's displayed "`f⁻_{3;2}|crystal = 1`" has the opposite sign to their general formula and to the definition. This is
  a convention-level sign.)
* **`f₄/f₂ = f₅/f₃ = ν⁴/(m²(m²−ν²)) = u²/(s(s+u))`** at `s = m²`, `u = −ν²`: exactly the two extra subtractions of CHLPS's kernel
  `(2s+u)/(s(s+u))^{k/2+1}`.
* The structures are unchanged. So `EE`'s crystal identity `EE ≡ q²` (EE_FUNCTIONAL §1.2) carries over: `Ĉ₄` is helicity-flip-free at
  the bulk point, with `q²`'s stripped density.
* Its flat image is CHLPS's `B₄⁽¹⁾(p)` (MHV, `s↔t` symmetric). `Ĉ₅`'s is `B₅⁽¹⁾(p)` (as `Ĉ₃ → B₃⁽¹⁾`) [P given CLS-D crystal (Q) +
  H0 §17 polarisation map].
* The forward heavy actions are those of Stage 0: `2u²`; `(c₄ − 2𝒥)u³` with `c₄ = 1 (++)`, `41 (+−)`; and `σ_hu³`.
* The all-plus rule `B₄⁽⁴⁾` (for `|ĝ₄| ≤ g₄`) is still not obtained. It is not a shock rule (LIT_HELICITY_MIXING §1) [I]. It is not
  needed for ABJM (`ĝ₄ = 0`).

### 1.4 Leading-Regge heavy image: not pure harmonics (`stage12_rawkernel.py`)
* **Multiplication by `η`.** It acts on `H₂` harmonics by the Legendre recurrence. The `k = 4` leading kernel is `η³ ×` (f-stripped
  kernel), so its heavy image is `Σ_{j=±1,±3} a_j(l)E_b(ν−ij)P_{l+j}(cosh β)`, with:
  * `l = −½ + iν`;
  * `E_b(ν) = 2^{2b−2}π Γ((b−½−iν)/2)²Γ((b−½+iν)/2)²` from v3_regge_kernel's `K(c,λ)`;
  * `b = 9` for `k = 4`, `b = 7` for `k = 2`.
* **Validation.** For `η¹`, `b = 7` this reproduces CHMRS's magic-block ratio `g/g̃` exactly (`10⁻⁴¹` at `ν = 0.3, 2.7, 9.1`) [DV].
* **Scalar proxy at `ν → 0`.**
  * Forward weights (normalised): `0.35 (j = ±3)`, `0.15 (j = ±1)`.
  * The raw normalised image `H(β)` equals 1 at `β → 0` (flat). It grows like `P_{5/2}(cosh β)`: `H/P_{5/2} = 1.00, 1.25, 1.81, 2.98` at
    `β = 0.01, 2, 4, 8`, and is positive throughout [DV, scalar proxy].
* **For `EE`.**
  * The `k = 2` rule "EE = scalar × `c(ν′)` per harmonic" gives `c(±i) = −0.297`, `c(±3i) = −61.1`: same sign, so if the rule held at
    `b = 9` the `j = ±3` term would dominate even more.
  * Whether it holds at `b = 9` is not computed [I].
* **Consequences.**
  * A **pure-moment (de-shifted) `Ĉ₄`** requires inverting a 4-shift difference operator (3-shift for `Ĉ₅`: `P_l`, `P_{l±2}`). This is
    the `k = 4` analogue of Result D's Abel problem, with possibly more null directions and real zeros of the forward normalisation.
    **Not constructed [I].** A forward (`ν = 0`) pure-moment functional is then an `O(1)`-width smearing in `ν` around 0 (flat width
    `1/Λ` in `p`), which has the correct flat limit.
  * The **raw** forward functional is a legitimate sum rule (the de-shift is only a re-weighting in `ν`). If the scalar-proxy sign holds
    for `EE`, its large-`J` action is positive and grows like `(J/m)⁵`. That would help positivity but makes the Regge region
    forward-dominated. Swappability then needs K2-type control at `η → ∞`, where the raw kernel grows like `η³` [I].
  * None of this touches the light side (§2), which is linear in `Ψ_ν` and vanishes on Einstein for every `ν`.

### 1.5 Does the `2𝒥` growth of `∂B₄|₀` survive curvature? (`stage12_forward_regge.py`)
* **Pure-moment model.** A pure Regge moment, flat-normalised, acts as `(flat normal form) × P_{−½+iν}(cosh β)`,
  `cosh β = 1 + (2J+1)²/(2m²)` (CHMRS; `EE`: same in both channels). At `ν = 0`:
  * `B₄|₀ → 2u²P₀(β)`;
  * `∂_{p²}B₄|₀ → −2u²Λ²D(β) + c₄u³P₀`;
  * `B₅|₀ → σ_hu³P₀`;

  with `P₀ = P_{−½}(cosh β)` and `D := −∂_{ν²}P_{−½+iν}(cosh β)|₀`.
* **Positivity [P].** From Mehler–Dirichlet (standard; checked against mpmath to `3·10⁻¹¹`),
  **`D(β) = (√2/2π)∫₀^β t²(cosh β − cosh t)^{−½}dt > 0`**. So the growth term has the flat sign in **both** channels for every `β`.
* **Rate [DV].** Let `R := D/P₀`.

  | `β` | 0.01 | 0.3 | 1 | 2 | 4 | 8 | 12 | 20 |
  |---|---|---|---|---|---|---|---|---|
  | `R/(β²/4)` | 1.0000 | 0.9991 | 0.990 | 0.966 | 0.909 | 0.831 | 0.789 | 0.747 |

  * Flat region: `Λ²R → Λ²β²/4 = u(J+½)² ≈ 𝒥u` (Stage 0's `2𝒥u³`).
  * For `β ≳ 1` (`J ≳ m`) the growth continues as `≈ (0.75–1)Λ²β²/4`, i.e. like `log²(J/m)`, against the `O(1)` constants `c₄`, `c/a`.
  * All forward terms share the factor `P₀ ~ βe^{−β/2}`. The smeared terms are suppressed relative to them by `~1/(Λβ)²` at large `β`.
* **Answer to Stage 0's [I]: the N3 mechanism survives AdS curvature, with the right sign, in the pure-moment leading-Regge model.**
* **Full Stage-0 functionals in this model** [M] (`c = ½`; physical states `b̂ = Λβ ≥ (2J_min+1)√u`; min of `FF_R/Nrm`):

  | `Λ` | functional | all physical | `b̂ ≥ 20` | `β ≥ 0.1` | `β ≥ 1` |
  |---|---|---|---|---|---|
  | 100 | `S_min + B₃` robust | −0.010 (`b̂ = 2.9`) | **+0.62** | +0.17 | +0.86 |
  | 100 | `S_min` robust | +0.025 | +0.70 | +0.38 | +0.98 |
  | 300 | `S_min + B₃` | −0.14 (`b̂ = 7.8`) | **+0.75** | +0.75 | +0.96 |
  | 1000 | `S_min + B₃` | −0.14 (`b̂ = 8`) | **+0.27** | +0.76 | +0.97 |
  | 1000 | `r = 5`, `ψ₃(0) = 0` | −0.15 | +0.43 | +0.99 | +0.99 |

  * The negatives sit at `b̂ < 10` (`J ≲ 4`). There the large-`J` model is not valid, and the exact flat LP (true `d`-functions) is
    positive.
  * Every state with `b̂ ≥ 20` is positive, and the curvature region (`β ≥ 1`) has margin `≥ 0.86`.
  * **This is not Stage 3.** Beyond-leading-Regge and near-gap corrections (Result D's `H_M`, `H_LJ` analogues) are not touched.

### 1.6 Verdict on the lift
**Possible; no obstruction found.**
* The forward pieces are CLS-D subtracted light-ray sum rules with:
  * the right flat images [P given Q];
  * more Regge margin than Result D's `k = 2` [Q];
  * no helicity flips at the crystal [P given] or at leading Regge [D];
  * the `∂B₄` growth intact [P/DV in the pure-moment model].
* Open [I], all heavy-side:
  * the pure-moment construction (4-/3-shift de-shift) or, alternatively, control of the raw kernel at `η → ∞`;
  * the K2-type edge bookkeeping with the `1/(u′v′)` poles;
  * all per-state positivity beyond leading Regge.

---

## 2. Stage 2: the exact light side (the gate)

### 2.1 Lemma A (annihilation of the dangerous terms)
**Statement.** Let `k ∈ {4, 5}` and let `𝒢` be an `O(1/c_T)` tree-level `⟨TTTT⟩` correlator that is:
* (i) microcausal, i.e. the `O(1/c_T)` coefficient of a family of microcausal correlators;
* (ii) of Regge intercept `J₀(𝒢) < 4`;
* (iii) inside the domain of CLS-D's validity statement at the light-cone edges and endpoints (V₄).

Then `Ψ_{k,ν}[𝒢] = 0` for every `ν`. Hence `∂_ν^nΨ_{k,ν}[𝒢]|₀ = 0`, and every smearing or de-shift `∫X(ν)Ψ_{k,ν}[𝒢]dν = 0`.

**Proof.** The integrand contains the commutator of two light-ray operators on the same null plane, which vanishes pointwise for a
microcausal `𝒢` once the integral converges (CLS-D §2, §3.5.2; KKSDZ §4) [Q]. Convergence is CLS-D's count `k > J₀` (`k − 1 > J₀` for odd
`j′`) plus V₄. ∎

**Applications** (the light side `Σ_{light} = Ψ[1 + T + [TT]]` restricted to each bulk term equals `Ψ` of that term's correlator, by
swapping its own s-channel OPE, with margin `4 − J₀ ≥ 1` [D]):

| correlator | `J₀` | source | Lemma A |
|---|---|---|---|
| Einstein (`T` in s, t, u + EH quartic) | 2 | Costa–Hansen–Penedones 1707.07689 (spin-2 exchange: `j(ν) = 2`) [Q]; flat `G s²/t` | **annihilated** |
| `E×W³` (one `W³` vertex: `t₄`, `α_T`) | ≤ 3 | u-channel exchange `j = 2` [Q]; s-, t-channel contact-like with Mellin degree ≤ 3 (dimension: `Gĝ₃p⁶`) [D] | **annihilated** |
| MFT (`O(c_T⁰)`, identity + `[TT]`) | ≤ 1 | KKSDZ (unitary GFF) [Q] | **annihilated** (s-channel identity + low-`n` `[TT]` sum to 0) |
| `W³×W³` (`|ĝ₃|²`) | 4 | dimension `G|ĝ₃|²p¹⁰/p²` | not annihilated (§2.2) |
| `R⁴` contacts (`g₄`, `ĝ₄`) | 4 | Mellin degree 4 | not annihilated (§2.2) |
| `g₅` and higher contacts | ≥ 5 | — | not annihilated (§2.2) |

**Grade.** [P given V₄ + (i) + the `J₀` counts]. (i) is [D] (order-by-order vanishing of commutators of the `1/c_T` expansion). The
`E×W³` count is [D], the Einstein one [Q]. Overall **[D]**. This is the same class of input as Result D's I14(a) [P given K1, K2], with
the subtraction factor added.

**What it replaces.** `ADS4_LIGHT.md`'s input (E) and Stage 0's (E)-type risk concerned decaying `O(G)` residues
`16πG(ν²+9/4)^{−n}` in `C⁺₄`, `C⁻₅` that the flat limit (F) cannot exclude. Lemma A excludes them at every `ν` without any Regge-moment
or Laurent-coefficient bookkeeping.
* The CHMRS residue formula (Lemma R) gives the same answer: Regge spin 2 has no `r^{−3}`, `r^{−4}` Laurent coefficient.
* Lemma A is stronger: it needs neither R4 nor R5, nor the moment expansion.

**Consistency checks.**
* Flat: for `f = 8πG/(stu)`, the `B_{k≥4}` low side is a sum of residues of `dz/z^{k/2+2}` (`z = s(s−q)`), hence zero [P].
* Scalar AdS: CHMRS state "The graviton exchange grows as spin two … and therefore contributes only to the `k = 2` moments", and
  `−C_{4,ν}|_{g₄} = 4` exactly (`sec_light.tex`, eqs. `C2C4Light`, `MGravityRegge`, read from source) [Q].
* `SUPERCONVERGENT_UPLIFT.md` §4 made the same annihilation argument for its non-meromorphic `Ψ_k` (proven reach there `J₀ < 3`).
  Here the meromorphic CLS-D factors give `J₀ < 4` [Q].

### 2.2 The `J₀ = 4` terms: flat-limit identities (F) + polynomial structure (P)
* **(F)** At fixed `p`, `Λ^{2k}Ĉ_k(Λp)|_light → B_k(p)|_low`. The light-side analogue of the crystal dictionary: [Q for scalars, CHMRS
  `BknuFlatSpaceLimit`; D-sketch for spinning].
* **A consistency argument for (F) [D].** A heavy exchange Witten diagram of spin `< 4` is annihilated by `Ψ₄` (Lemma A). So its
  induced `[TT]` ("contact") light side equals minus its heavy-block action, which is flat-normalised. Contacts generated by heavy states
  therefore inherit the heavy normalisation.
* **(P)** The light side of a contact-like term is a polynomial in `ν²`. [Q for scalars, CHMRS `BknuContactsLeading`: "`C_{k,ν}[M_{a,b}]` is
  a polynomial in `ν²` of degree at most `2a+3b−k`", `−C_{2a,ν}[M_{a,0}] = 2^a`; D-sketch for spinning (Lemma R's R4).]
* **The flat-unit argument.** In AdS units the light side of a term with coupling `g` is `g × (a pure function of ν)`; `Λ` enters only
  through the dictionary. Under (P) that function is a polynomial whose top coefficient is fixed by (F). Every other coefficient
  multiplies a positive power of `Λ` in flat units, so it must vanish if the flat coefficient is zero. Hence [P given (F), (P)]:

  | term | `Ĉ₄(ν)` | `∂_{ν²}Ĉ₄(0)` | `Ĉ₅(ν)` |
  |---|---|---|---|
  | `g₄`: `Λ⁸g₄,AdS = g₄,f` | `−2g₄` (constant; any `ν²`-term would scale as `Λ²`) | 0 | 0 (would scale as `Λ²`) |
  | `ĝ₄`: flat `B₄⁽¹⁾ ∌ ĝ₄` | 0 | 0 | 0 |
  | `W³W³`: `Λ⁸G|ĝ₃|²_AdS = G_f|ĝ₃|²_fΛ⁻²` | `−4πG|ĝ₃|²ν² − b₀G|ĝ₃|²` | `−4πG|ĝ₃|²` (exact) | 0 |
  | `g₅`: `Λ⁸g₅,AdS = g₅,fΛ⁻²` | `−g₅(ν² + b₅)` | `−g₅` | `−g₅` (exact) |

  * The `g₅` cancellation in `∂Ĉ₄ − Ĉ₅` survives exactly (ADS4_LIGHT item 5).
  * `b₀`, `b₅` are AdS numbers of the bulk Lagrangian's curvature completion. They are not computed [I].
* **Without (P).** Non-polynomial contact or `|ĝ₃|²` terms (e.g. `g₄/(ν²+a)`) would survive at `ν = 0` with `O(1)` coefficients.
  * For `g₄` they would make `g̃₄ := −½Ĉ₄(0)|_light` differ from the bulk `g₄` at relative `O(1)`. The bound would then still be a
    correct bound on the CFT quantity `g̃₄`, but its identification with the bulk coupling, and the ABJM evaluation, would need (P).
  * For `|ĝ₃|²` the leverage is tiny anyway (§2.4).
  * **(P) is not load-bearing for the gate; it is load-bearing for the interpretation and for ABJM.**

### 2.3 What is and is not "exact"
* **Exact (given the stated inputs):**
  * every Einstein and `E×W³` entry of the table in §0 (zero);
  * the `g₄` constant `2g₄`;
  * the `∂Ĉ₄` coefficients `4πG|ĝ₃|²` and `g₅`;
  * the `Ĉ₅` coefficient `g₅`;
  * Result D's pole `16πG/(ν²+9/4)` in the `k = 2` piece (H7, [D]).
* **Not computed [I]:**
  * `b₀` (AdS constant of `W³W³` in `Ĉ₄(0)`);
  * `b₅` and the higher-contact curvature constants;
  * Result D's `c₁`.

  §2.4 shows they are harmless at the gate level.
* **Explicit Witten diagrams: not evaluated.** A direct numerical evaluation of `Ψ₄`, `Ψ₅` on Einstein `⟨TTTT⟩` Witten diagrams (≥ 12
  digits, as the task asked) was **not possible here**.
  * No position- or Mellin-space tree-level `⟨TTTT⟩` is available in a form the light-ray integrals can use.
  * Raju (1201.6452) is momentum space. Its abstract: "We write down an explicit formula for the four point correlator with arbitrary
    external helicities and momenta. We verify that … the MHV amplitude of gravitons … appears as the coefficient of a specified
    singularity" [Q, abstract]. This supports (F) for Einstein, not Lemma A.
  * Lemma A makes such an evaluation unnecessary for the zeros. It would be a valuable independent check of V₄ (§4).
* **Cross-checks performed:**
  * the CHMRS scalar statements (§2.1) [Q, source];
  * the flat `z`-residue zero [P];
  * the BCP/Penedones chain for the `R⁴` normalisation (§3) [DV];
  * the `k = 2` magic-block reproduction of the kernel machinery (§1.4) [DV].

### 2.4 The leverage gate (`stage12_leverage.py`, `cache/stage12_leverage.txt`)
**Setup.** With `c = ½`: `g₄/(8πG) ≤ 2L_AdS − κx² − a(δ_{∂B₄} − δ_{B₅})/(8πG) − δ_{B₄}/(16πG)`, where `δ_X` is any extra light-side term.
"Leverage" means the shift of the bound on `g₄Λ⁶/(8πG)` per unit of the input. The signal is `2L_AdS = 35 / 48 / 61` (robust full set).

| input | `Λ = 100` | 300 | 1000 | status after §2.1–2.2 |
|---|---|---|---|---|
| (i) unit residue `16πG(ν²+9/4)^{−4}` in `−Ĉ₄(0)` | `3.9·10¹⁰` | `2.8·10¹³` | `3.9·10¹⁶` | **× 0** (Lemma A) |
| (ii) its `∂_{ν²}` in `−∂Ĉ₄` (× `a`) | `1.4·10¹⁴` | `9.1·10¹⁷` | `1.4·10²²` | **× 0** |
| (iii) unit residue `(ν²+9/4)^{−5}` in `−Ĉ₅(0)` | `3.5·10¹³` | `2.3·10¹⁷` | `3.5·10²¹` | **× 0** |
| (iv) unit linear `16πGt₄` in `−Ĉ₄(0)` (`x = 10`) | `2.4·10⁷` | `2.2·10⁸` | `2.4·10⁹` | **× 0** (`J₀(E×W³) ≤ 3`) |
| (v) `b₀G|ĝ₃|²` in `−Ĉ₄(0)`, per unit `b₀`: `x² ≤ 110` / Result D's bound | `2.2·10⁻⁴` / `1.5·10⁻²` | `2.4·10⁻⁵` / `2.4·10⁻³` | `2.2·10⁻⁶` / `2.9·10⁻⁴` | remains [I], harmless for `|b₀| ≲ 10³` |
| (vi) `b₅g₅` (relative, `|g₅| ≤ g₄`), per unit `b₅` | `5·10⁻⁵` | `6·10⁻⁶` | `5·10⁻⁷` | remains [I], harmless |
| (vii) Result D's `16πGc₁t₄` in the `k = 2` piece, per unit `c₁` (`x = 10.5` / D's bound) | `0.35` / `2.9` | `0.04` / `0.38` | `0.003` / `0.04` | remains [I]; `O(1)` at `Λ = 100` |

Notes on the table:
* `x² ≤ 110` is Stage 0's own `|ĝ₃|²` closure with matter. Result D's bound is `x² ≤ 16L_A = 7.5·10³ / 1.07·10⁴ / 1.43·10⁴`.
* (Stage 0 quoted `~Λ⁸` for (i). In `g₄`'s own normalisation the recount gives `(4/9)⁴Λ⁶` for (i) and `Λ⁸` for (ii), (iii). The
  conclusion is the same.)

**Gate verdict: residual leverage ≤ 3 ≪ 10² ⇒ GO for Stage 3** (heavy positivity).
* **Conditions:** V₄ and (i) of Lemma A, at grade [D] overall. They are not decided by numbers here.
* **If V₄ failed** (some boundary contribution of the `1/(u′v′)`-weighted light-ray integral on a `J₀ = 2` correlator), the zeros
  would become unknown residues, and rows (i)–(iv) would apply: **then the precise no-go is** "a forward `g₄` bound in AdS₄ needs the
  `ν = 0` light side to be known to `10⁻¹¹ – 10⁻²²` relative to `G`; only an exact vanishing theorem can supply that".

**Exact IR log [P given the pole (D)].** The `k = 2` pole makes Stage 0's [M] translation unnecessary:
`g₄Λ⁶/(8πG) ≤ 2L_AdS(Λ) − κx² − …`, with `L_AdS = ∫₀¹ψ₂(p)dp/(p² + 9/(4Λ²))`.

| `Λ` | 100 | 300 | 1000 |
|---|---|---|---|
| `2L_AdS`, robust `S_min + B₃` | 35.272 | 47.719 | 61.471 |
| Stage 0 `A log(2Λ/3) + B` | 35.101 | 47.678 | 61.462 |
| `2L_AdS`, robust `S_min` | 47.159 | 63.552 | 81.598 |

---

## 3. ABJM: the normalisation of the `R⁴` contact (`stage12_abjm.py`, `cache/stage12_abjm.txt`)

**Chain** (exact, sympy; inputs [Q] from the BCP 1906.07195 source in the scratchpad):
1. BCP eq. (flatSSSS) maps the MHV function `f₁` to the six `⟨SSSS⟩` flat amplitudes.
   * `f₁ = 1/(stu)` gives `(tu/s, su/t, st/u, s/2, t/2, u/2)`, which is BCP's `M_SG` large-`s,t` form (SGLimit, `g₁ = 1`).
   * `f₁ = 1` gives `(t²u², …)`, BCP's `M₄` top form (largeA4).
   * Ratio `stu` in every component [P].
2. Penedones' formula as used by BCP (FlatLimit): `𝒜 = 𝒩√π∫dα/(2πi) e^α α^{−½}M(L²s/(2α), L²t/(2α))`. A degree-`n` term maps with
   `(L²/2)ⁿ/Γ(n+½)` (`h = ΣΔ/2 − d/2 = ½`). Degree 4 vs degree 1: `L⁶/105` [P].
3. With `𝒜 = 𝒜_SG[1 + ζ(3)ℓ_s⁶stu/32]` (BCP IIAAMP) and `L⁸/ℓ_s⁸ = 4π⁴λ²` (BCP `class`): **`B⁴₄/B¹_SG = 105ζ(3)/(64√2π³)` is reproduced
   exactly** (BCP IIARes1) [DV].
4. Result-Q variables.
   * `f = 8πG/(stu)·(1 + (g₄/8πG)stu + …)`, so `g₄/(8πG) = ζ(3)α′³/32` (`ℓ_s² = α′`).
   * `M² = 4/α′`, the first massive level, which contains VS's `J = 4` `(+−)` state (Stage 0 §5).
   * `Λ = ML = 2^{5/4}π^{1/2}λ^{1/4} = 4.2156λ^{1/4}`.
   * **Result: `g₄Λ⁶/(8πG) = 2ζ(3)`** [DV given Q].
5. **CFT-side recipe.** Let `ρ_M := λ^{−3/2}B⁴₄/B¹_SG` be the ratio of top-degree `R⁴` and SUGRA Mellin coefficients in `⟨SSSS⟩`.
   * In AdS units, **`g₄/(8πG) = ρ_M Γ(3/2)/(8Γ(9/2)) = ρ_M/105`**. With BCP's value this gives `2ζ(3)/Λ⁶` exactly.
   * General rule: for a component correlator with `h = ΣΔ_i/2 − d/2`, SUGRA top degree `n_SG` and `R⁴` top degree `n_SG + 3` in the
     same flat helicity structure, `g₄/(8πG) = ρ_M(½)³Γ(h+n_SG)/Γ(h+n_SG+3)`.
   * For a direct `⟨TTTT⟩` Mellin amplitude, `h` and `n_SG` depend on the spinning Mellin convention [D-sketch]; the rule is stated, not
     applied.
6. **In the Result-Q functionals no Mellin convention enters.** `−Ĉ′_EE`'s pole coefficient `16πG` and `−Ĉ₄(0)|_{R⁴} = 2g₄` are both
   heavy-normalised (§2.2). So `−Ĉ₄(0)|_{R⁴}/(16πG) = g₄/(8πG)`, with `8πG = 256/c_T` (ADS4_LIGHT, [Q + D]). The ABJM tree correlator
   gives `g̃₄Λ⁶/(8πG) = 2ζ(3)(1 + O(λ^{−1/2}))`.
   * The `O(λ^{−1/2})` (`≈ 18/Λ²`) comes from degree-4 Mellin pieces at higher order in `1/λ` (BCP's `B⁵₄M₄` at `λ^{−2}`, coefficient
     not used here) and from the `ν = 0` constants of `D⁴R⁴, …` [I, sizes `O(10⁻³)` at `Λ = 100` unless the coefficients are `≫ 1`].
   * `λ → λ − 1/24` changes `Λ` by `< 10⁻⁵` relative for `λ ≥ 10³`.
   * Genus one (finite `μ`, outside A2): `2ζ(3) → 2ζ(3) + (2π²/3)g_s²` [Q].
7. **Placement** (exact-pole bound, robust sets):

   | `λ` | `Λ` | `S_min + B₃` bound | ABJM fraction | `S_min` bound | ABJM fraction |
   |---|---|---|---|---|---|
   | `3.2·10⁵` | 100.3 | 35.3 | 6.8 % | 47.2 | 5.1 % |
   | `2.6·10⁷` | 301.0 | 47.7 | 5.0 % | 63.6 | 3.8 % |
   | `3.2·10⁹` | 1002.7 | 61.5 | 3.9 % | 81.6 | 2.9 % |

   **ABJM (IIA) is inside**, at 3–7 % [D given Q; the AdS bound itself awaits Stage 3]. As Stage 0 noted, this is a weak test (log-dominated
   margin).

---

## 4. What would make the gate a theorem, and the remaining risks

| item | role | grade | how to settle |
|---|---|---|---|
| **V₄**: CLS-D validity (vanishing + absolute convergence) of `Ψ₄[f⁺₄, EE]`, `Ψ₅[f⁻₅, o_H q]` on tree correlators with `J₀ < 4`, incl. the `1/(u′v′)` light-cone poles and `EE`'s `x₀₂²`-edge | **load-bearing** for the gate | [Q (CLS-D statement) + D (application)] | write out K2 for `k = 4, 5`; or one independent numerical evaluation of `Ψ₄` on a known AdS₄ spinning tree correlator (e.g. Einstein `⟨TTTT⟩` from Raju's momentum-space formula, transformed), or on the MFT `⟨TTTT⟩` |
| (i) microcausality of the `O(1/c_T)` tree correlators | load-bearing | [D] | standard; write out |
| `J₀(E×W³) ≤ 3` | load-bearing for the `t₄`-linear zero | [D] (Mellin degree / dimension) | Mellin degree of the `W³`-vertex exchange (Sleight–Taronna-type formulas, unverified) |
| (F) light-side flat dictionary for spinning functionals | normalisation of `g₄`, ABJM | [D-sketch] (heavy-normalisation consistency argument §2.2) | CLS-D's stated open problem ("derive a similar dictionary … low energy EFT contribution") |
| (P) polynomial contact light sides | identification `g̃₄ = g₄`, ABJM | [Q scalar / D-sketch spinning] | spinning Lemma R (R4) |
| `b₀`, `b₅`, `c₁` | residual leverage ≤ 3 | [I] | not needed at the gate |
| de-shift for `k = 4, 5` (or raw-kernel `η → ∞` control) | heavy side (Stage 3) | [I] | `k = 4` analogue of EE_FUNCTIONAL §3 |
| `B₄⁽⁴⁾` (all-plus) | `|ĝ₄| ≤ g₄` | [I] | non-shock rule; not needed for ABJM |

---

## 5. Graded claims

| claim | grade |
|---|---|
| frame cross-ratios; Laurent expansions of `f^±_{k;2}`; `f₄/f₂ = f₅/f₃ = 16η²r²(1+O(r))` | [P] (sympy) + [DV] |
| crystal values of `f_k`, `k = 2…7`; `f₄/f₂ = u²/(s(s+u))` | [P] (sympy; CLS-D crystal [Q]) |
| flat images `B₄⁽¹⁾`, `B₅⁽¹⁾`, flip-free, `q²` density | [P given CLS-D crystal (Q) + H0 §17 map] |
| convergence `J₀ < 4` for `k = 4, 5` | [Q] (CLS-D) |
| leading-Regge off-diagonal zero for `Ĉ₄[EE]` | [D] |
| `η³`/`η²` mixing; scalar-proxy raw image grows like `P_{5/2}`, positive | [P] (recurrence) + [DV] (scalar proxy); `EE` [I] |
| `D(β) > 0` for all `β`; table of `R(β)` | [P] + [DV] |
| Stage-0 functionals positive on `b̂ ≥ 20` in the pure-moment leading-Regge model | [M] |
| **Lemma A**: Einstein, `E×W³`, MFT light sides of `Ψ₄`, `Ψ₅` vanish identically (all `ν`) | **[D]** ([P given V₄, (i), `J₀` counts]) |
| `R⁴`: `−Ĉ₄(0) = 2g₄`, `∂Ĉ₄`, `Ĉ₅` zero; `|ĝ₃|²` top `4πG`; `g₅` exact; `Ĉ₅ ∌ |ĝ₃|²` | [P given (F), (P)] |
| residual leverage ≤ 3 (per unit of the `O(1)` unknowns) | [D] (arithmetic [DV]) |
| exact IR log `2L_AdS` | [P given the pole (D)] + [DV] |
| BCP coefficient from Penedones; `g₄Λ⁶/(8πG) = 2ζ(3)`; recipe `ρ_M/105` | [DV given Q] |
| ABJM at 3.9–6.8 % of the bound | [D given Q]; AdS bound [I] until Stage 3 |
| **gate: GO** | [D] (conditional on V₄, (i)) |

---

## 6. Files (`ads4/code/resultq/`)

| file | content | run time |
|---|---|---|
| `stage12_kernel.py` | frame cross-ratios, Laurent expansions of `f^±_{k;2}` (`k = 2…8`), crystal values, flat ratios → `cache/stage12_kernel.txt` | ~1 min |
| `stage12_rawkernel.py` | `η`-recurrence heavy images; `k = 2` check against CHMRS `g/g̃`; `k = 4` scalar proxy at `ν → 0` → `cache/stage12_rawkernel.txt` | seconds |
| `stage12_forward_regge.py` | Mehler–Dirichlet `D(β)`, `R(β)`; pure-moment leading-Regge model of the Stage-0 functionals → `cache/stage12_forward_regge.{txt,log}` | ~5 min (2 cores) |
| `stage12_leverage.py` | leverage table, exact IR log, `c₁` item → `cache/stage12_leverage.txt` | seconds |
| `stage12_abjm.py` | BCP/Penedones chain, `2ζ(3)`, placement → `cache/stage12_abjm.txt` | seconds |

---

## 9. Hand-over

**State.** Stages 1 and 2 are complete at the level above.
* Lift: possible.
* Gate: **GO**, conditional on V₄ and (i) of Lemma A (grade [D]).
* ABJM normalisation: fixed (`2ζ(3)`, recipe `ρ_M/105`).

**Next steps (in order of value):**
1. **Make V₄ explicit for `k = 4, 5`.** K2-type power counting at every boundary with `f⁺₄`, `f⁻₅`, `EE`. This is the only load-bearing
   input of the gate.
   * A cheap independent check: evaluate `Ψ₄` (raw) numerically on the MFT `⟨TTTT⟩` (known in closed form; must give 0). This tests
     the `1/(u′v′)` edges.
2. **Stage 3 prerequisites.**
   * Construct the pure-moment `Ĉ₄` (4-shift Abel problem: exact `h_j` via the v4_hclosed/u2 machinery extended to `η³`, i.e. rank-3
     equivariant tensors), or decide to use the raw forward functional and control `η → ∞`.
   * Compute the `EE` factor at `b = 9`, `11` (`stage12_rawkernel.py` gives the scalar part).
3. Stage 3 proper: per-state positivity of `Ω^Q` beyond leading Regge (near-gap rows, Regge corrections), as Result D did for `Ω^EE`.
4. Optional: `b₀` (one AdS `W³W³` exchange Witten diagram in the `u`-Regge limit), the spinning (P), and a direct `⟨TTTT⟩` Mellin
   check of the `R⁴` normalisation.

**Not done:** explicit Witten-diagram evaluations; the de-shift; `B₄⁽⁴⁾`; anything beyond leading Regge.
