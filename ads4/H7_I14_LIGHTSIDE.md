# H7 and I14 for route A′: the graviton-pole light side of `Ĉ₂`, `Ĉ₃`, and swappability

*Written 19 Nov 2026 by a derivation agent. Not refereed, not committed. It touches no other file. Code: `ads4/code/h7i14/`.*

**Grades used here.** [P] proved; [P given X] proved from the named inputs X; [D] derivation written out here, not checked
independently; [D-sketch]; [Q] quoted (source and section named); [DV] derived and checked numerically or symbolically; [NE]; [I] open.

**Sources.** The LaTeX sources were read in the session scratchpad.
* CHMRS: arXiv:2106.10274, `sec_light.tex`, `app_omega_light.tex` (App. E = `app:omegalight`), `sec_sumrules.tex`, `sec_heavy.tex`.
* CLS-D: arXiv:2311.04271, §2.1, §3.5.1–3.5.2.
* KKSDZ: arXiv:1904.05905, "Shocks, superconvergence…", §2 (flat shocks; §2.3 GR), §3.3 (AdS shocks, eq. `adsshock`) and §4 (light-ray products: §4.1 existence vs commutativity, §4.2.3–4.2.6 Rindler positivity, Regge and light-cone limits, §4.3 double commutator, §4.4 summary).
* CHP: Costa–Hansen–Penedones, arXiv:1707.07689, §2.3–2.5.
* KPZ: Kulaxizi–Parnachev–Zhiboedov, arXiv:1705.02934, §3.
* CMRS-D: arXiv:2008.04931, `sec_subtractions.tex` and `sec_MFTbounds.tex`.

Kravchuk–Simmons-Duffin arXiv:1805.00098 and Caron-Huot arXiv:1703.00278 were downloaded but are cited only through KKSDZ, who
extend their Regge-boundedness argument to spinning operators.

---

## 0. Verdict

| item | before | after this note | what remains |
|---|---|---|---|
| **H7**, `Ĉ₃` part: `Ĉ₃|_E = 0` exactly | [D-sketch] | **[P given Q2, Q8, F]**, plus an independent route [P given Lemma R + no Π₂ moment] | F: factorisation of the T-block term of the bare commutator through one-point shock matrices ([D], LIGHT_INPUTS `t2` [DV]) |
| **H7**, `Ĉ₂` part: `−Ĉ₂|_E = 16πG/(ν²+9/4)` exactly | [D-sketch] | **[D]: complete derivation, given Lemma R, (HN), Q3–Q5** | (HN) is the heavy normalisation premise, already load-bearing on the heavy side. D2 (CHP's conserved-operator impact-parameter form at both vertices) and D3 (plane-wave bookkeeping) are this note's [D]. Not [P]. |
| Lemma R (CHMRS residue formula for spinning `Ĉ₂`, `Ĉ₃`) | [D-sketch] | **[P given R1–R5]** | R4 (tree-level `⟨TTTT⟩` has a single-valued integer Laurent expansion at `r = 0`) is [D-sketch]; R5 is the decay assumption of CHMRS App. E |
| **I14(a)**: `Ĉ₂`, `Ĉ₃` annihilate the non-perturbative `⟨TTTT⟩` (convergence + microcausality) | [I] | **[P given K1, K2]** | K1: KKSDZ's Rindler-positivity Regge bound, with their stated loophole. K2: kernel power counting at every boundary, [Q] at the Regge corner and [D] at the endpoints |
| **I14(b)**: heavy sum and functional commute | [I] | **[P given K1, K2, K3]** | K3: per-family Rindler positivity, [D] here. K2 for the total de-shifted `Ĉ₂` kernel as `η → ∞` is **[D-sketch]** |
| Is `j′ = 0` beyond the proven Regge bound? | — | **No.** Both functionals need `J₀ < 2` (CLS-D §3.5.2). Proven: `J₀, J₀^{dDisc} ≤ 1` (KKSDZ §4.2.4 `sec:reggelimit`, §4.3 `sec:doublecommutatorconvergence`). Margin: one unit | — |
| H16 sub-items P1(`j′=0`), `j_k(ν)`, "Einstein absent from `Ψ_k`" | [I] | **Not needed under route A′** (§8) | — |

**Net effect.** The light-side constants `A`, `B_AdS` rest on:
* (HN);
* the literature inputs Q1–Q5, Q8;
* two [D] steps (D2, D3);
* I14, now [P given K1–K3], with one [D-sketch] at `η → ∞` for the de-shifted `Ĉ₂`.

The derivation is **not** [P]. The honest overall grade is: **H7 [D] (was [D-sketch]); I14 [P given K1–K3], with K2 partly [D-sketch] (was [I])**.

---

## 1. Statements

Notation: `Ĉ_{2,ν}`, `Ĉ_{3,ν}` are as in `H0_SPINNING_FUNCTIONAL.md` §R, with the §24.1 definition.
* `Ĉ₂` is the CHMRS-type (inverse magic-block, "de-shifted") combination of the CLS-D `Ψ_{2,ν′}[q²]`.
* `Ĉ₃` is the CLS-D `j′ = 1`, `f = −1` functional with the `o_H·q` structure, combined in the same way.
* Both are normalised by `N̂_k(ν)` (§R.3).
* `𝒢^E` is the tree-level (`O(1/c_T)`) `⟨TTTT⟩` of pure Einstein gravity in AdS₄ (`t₄ = α_T = 0`). "`|_E`" is the part of the light
  side ∝ `λ²_{TTT,E}`. Cross terms `E × t₄` and the `|ĝ₃|²` terms are H6/H17, not H7.

* **(L2)** `−Ĉ_{2,ν}|_E = 16πG/(ν² + 9/4)` for all real `ν`.
* **(L3)** `Ĉ_{3,ν}|_E = 0` for all real `ν`.
* **(S)** For the full non-perturbative `⟨TTTT⟩`, `Ĉ_{k,ν}[𝒢] = 0` (k = 2, 3), and `Ω^{A′}[𝒢] = Σ_O p_O Ω^{A′}[G_O]`, with the sum absolutely
  convergent.

---

## 2. Literature inputs (each a statement in the cited source)

| id | statement | where | grade |
|---|---|---|---|
| **Q1** | **Residue formula.** Take a dispersive functional `ω[G] = ∫∫ dDisc_s[ℋ G]` over `w, w̄ ∈ (−∞, 0]`, with `ℋ` commuting with `dDisc_s` up to terms supported away from the u-channel Regge point. Let it be `s↔t` antisymmetric, with moment expansion `ω ∼ Π_{k,η} + O(Π_{n+1})`, and let `𝒢^{EFT}` be a tree-level EFT correlator of Regge spin `n`. Assume (i) `dDisc_s 𝒢^{EFT}` vanishes near the u-Regge limit, (ii) `𝒢^{EFT}_+` is single-valued with integer Laurent powers at `r = 0`, and (iii) the u-channel Euclidean terms decay. Then `−ω|_light = (2πi/4)[𝒢^{EFT}_+]_{r^{1−k}}`, and in general (`moregeneralformula`) `−ω|_light = ¼Σ_{k ≤ n}∫dη a_k(η)∮dr r^{k−2}𝒢_+`. | CHMRS `ssec:lightResidue` eq. `eftformula`, `moregeneralformula`; proof `app:omegalight` | [Q] |
| **Q2** | The CLS-D functional is `Ψ = Ψ⁺ − (1↔3)`, an integral of the double commutator `⟨[𝒪₄, L𝒪₃][L𝒪₁, 𝒪₂]⟩` (manifestly `1↔3` antisymmetric, the `s↔t` of Q1). It vanishes on physical correlators if the light-ray integrals converge. `f^±` change the Regge spin; convergence needs `k > J₀` (even `j′`) and `k − 1 > J₀` (odd `j′`). For `J_e = 2`, `k = 2, 3` are unsubtracted (zero on all double traces). | CLS-D §2.1 (eq. after `light_transform_vacuum`), §3.5.1–3.5.2 | [Q] |
| **Q3** | **Impact-parameter representation** of a Regge limit, `A = ∫dp dp̄ e^{−2ip·x−2ip̄·x̄} B`, with `ℬ = ∫dν S^{j(ν)−1}Σ_kβ_k(ν)𝒟̂_kΩ_{iν}(L)`, and `α(ν) = γ(ν)γ(−ν)β(ν)`; `β` has poles only at bulk fields. For conserved external `T`, `p·∂_z ℬ = 0`: all longitudinal structures vanish (`β₄ = … = β₁₀ = 0`), so `ℬ` lives on transverse tensors of `H_{d−1}` at fixed `E`; `𝒟̂₁ = 𝟙` (minimal coupling). | CHP §2.3–2.5, eqs. `AtoB`, `BRegge`, `AlphaBeta`, `conservation_B`, `conservation_cond_beta`, `DH3` | [Q] (for `⟨TT𝒪𝒪⟩`) |
| **Q4** | **Einstein shock.** The AdS shock `h ∝ z^d/(z²+y²)^{d−1}` is an exact solution of Einstein's equations, even with higher-derivative corrections (Horowitz–Itzhaki). In Einstein gravity the one-point (transfer) matrix is the identity on graviton polarisations: in `d = 3`, `H = [[1, (t₄/4)w^{−4}], [(t₄/4)w⁴, 1]]`, identity at `t₄ = 0` (LIGHT_INPUTS `t2` [DV]). The Witten exchange diagram's impact-parameter transform is the shock phase shift. | KKSDZ §3.3 eqs. `adsshock`, `shockprofile` (and §2.3 flat GR shocks); KPZ §3 | [Q] |
| **Q5** | **Heavy Regge moment.** `Π̂_{k,ν}[G_{Δ,J}] → γ_{2Δφ+k−1}(ν)² P(η_AdS)/m^{2k}` via the Fourier–Laplace identity applied to plane waves `e^{x·p}`. Spinning analogue: transverse-helicity plane waves give `γ_a²M(ν;a)²P(cosh β)`, the same `M_𝔖` for `++` and `+−`. | CHMRS `sec:reggemomentsheavy2` eqs. `Fourier identity`, `finalreggemomentofablock`; `REGGE_MOMENT_TTTT.md` Prop. §5 [D + DV pieces] | [Q] / [D + DV] |
| **Q6** | **Regge boundedness.** `J₀ ≤ 1` and `J₀^{dDisc} ≤ 1` for spinning Wightman functions and double commutators in unitary CFTs, by Rindler positivity and Cauchy–Schwarz (normalised by two-point functions, generic polarisations). The light-ray-product convergence conditions are listed in §4.4. KKSDZ flag "a small loophole" at the first/second-sheet boundary (§4.2.6 `sec:asymptoticLC`). | KKSDZ §4.2.3–4.2.4, §4.3, §4.4 | [Q] ("we prove", with loophole) |
| **Q7** | **Swapping for scalar dispersive functionals.** If the kernel gives a finite integral with `|𝒢|` replaced by a constant at `|w| → ∞` (spin-2 decay), the functional swaps with the OPE. For CHMRS's `C_{k,ν}` (and the `Π_{k,η}` inversions) swappability is "confirmed directly in examples" only. | CMRS-D `sec_subtractions.tex` (after eq. `cHNu01Decomposition`), `sec_MFTbounds.tex` ("Swappability" item); CHMRS `sec:reggemomentsheavy` ("Swappability … confirmed directly in examples") | [Q] |
| **Q8** | One-point structure `⟨T|E(n)|T⟩` in `d = 3`: Einstein part ∝ `𝟙` in helicity space; flips ∝ `t₄`, `α_T` (phase). | Hofman–Maldacena; KKSDZ §2.3–2.4, §3.3; LIGHT_INPUTS `t2` [DV] | [Q] |

---

## 3. Lemma R: the CHMRS residue formula for `Ĉ₂`, `Ĉ₃`

**Lemma R.** For `ω ∈ {Ĉ_{2,ν}, Ĉ_{3,ν}}` and `𝒢^{EFT}` the tree-level `⟨TTTT⟩` of a bulk EFT with Regge spin `n`, Q1 holds with
`𝒢^{EFT}_+` replaced by the kernel-contracted scalar function `𝒦_+(r, η) = Σ_a ℋ_a 𝒢^a_+`. **Grade: [P given R1–R5].**

**Proof.** CHMRS `app:omegalight` uses five things:
* that `𝒦(w, w̄) := ℋ𝒢` is a function of two variables, analytic off the cuts `(−∞, 0]`, `[1, ∞)`;
* that `I_left = −I_right` (antisymmetry);
* that the Euclidean (opposite-`iε`) terms decay at large `|w|`;
* that `dDisc_s 𝒦 → 0` near the u-Regge point;
* the single-valued Laurent expansion in `r` at fixed `η`.

It never uses that `𝒢` is a scalar correlator. A spinning functional reduces, by conformal invariance, to such an integral of the
contracted function `𝒦 = Σ_a ℋ_a 𝒢^a`. The tensor structures and the dual/celestial contractions are absorbed into `ℋ_a`. So the proof
carries over once:

* **R1 [Q, Q2].** `Ĉ` is a `dDisc_s` integral of `𝒦` over the CHMRS region: CLS-D's double-commutator form. CLS-D note that the
  scalar `Ψ_{k,ν}` is CHMRS's `C̃_{k,ν}`.
* **R2 [Q, Q2].** `1↔3` antisymmetry: `Ψ = Ψ⁺ − (1↔3)`, with the signature `(−1)^{j′}` of the celestial structure compensated by
  `f^{(−1)^{j′}}`. For `Ĉ₂`, `f = v′−u′` is odd and `j′ = 0`. For `Ĉ₃`, `f = −1` is even and `j′ = 1`.
* **R3 [D].** `dDisc_s` commutes with `ℋ`:
  * `f₂ = v′−u′` and `f₃ = −1` are polynomials;
  * the dual and celestial structures are products of powers of `x_ij²`, `z·x`, with `ν`-dependent exponents only through the
    SO(1,1)×SO(2,1)-covariant `η`-dependence (`P_{−½+iν}(η)` in the scalar case);
  * they are analytic near the left cut away from the Regge point, exactly as CHMRS's `C_{k,ν}` kernel.
* **R4 [D-sketch].** `𝒦^{EFT}_+` has an integer-power Laurent expansion at `r = 0`, at fixed `η`.
  * CHMRS assert it for tree-level scalar EFTs.
  * For `⟨TTTT⟩` it follows if the tree-level spinning Mellin amplitudes are polynomially bounded with Laurent expansions at
    large `S` (true for finite sums of exchange and contact Witten diagrams). The `d = 3` parity-odd building blocks contribute
    factors `∝ ρ − ρ̄ ∝ r√(η²−1)`, also single-valued in `r`.
  * Not written out for spinning diagrams.
* **R5 [as in CHMRS].** `dDisc_s 𝒢^{EFT}` decays along the negative axes and vanishes near the u-Regge point, and the u-channel
  Euclidean terms decay. True for finite sums of Witten diagrams (CHMRS `app:omegalight`). The spinning case is assumed, not checked.

**Consequence: the counting rule.**
* Only moments `Π_k` with `k ≤ n` contribute: `[𝒢_+]_{r^{1−k}} = 0` for `k > n`.
* Einstein: `n = 2` (Q3 with `j = 2`: `ℬ ∝ S^{1}`; tree GR has Regge spin 2 in every helicity component [D, flat]).
* So:
  * only the **leading `Π₂` moment** of `Ĉ₂` sees `𝒢^E`;
  * `Ĉ₂`'s `Π₃, Π₄, …` admixtures are irrelevant for H7. In particular the question whether `Ĉ₂`'s Regge expansion has only even
    steps (§15.2's "[I-minor]") **does not matter for the Einstein term**; it matters only for `|ĝ₃|²` (H17).
* Check (d) of `v1` [DV]: `(1/4)∮dr r^{k−2}Σ_{n≥−1}c_nr^n = (iπ/2)c_{−1}` for `k = 2`, and `0` for `k = 3, 4, 5`.

**What Lemma R does *not* need.**
* That the parity-odd `[TT]` double traces have double zeros (§15.1 caveat).
* That contacts drop out.

The light side is defined as the full tree-level light contribution, `∫∫dDisc_s[ℋ𝒢^{EFT}]`, which includes whatever the double
traces do. The residue formula evaluates it. For the Einstein part this is the residue of the full `𝒢^E` (exchange in all channels
plus the quartic contact).

---

## 4. (L3): `Ĉ₃|_E = 0`

Two independent routes.

**Route A (bare commutator; [P given Q2, Q8, F]).**
* `Ĉ₃` has `f = −1`, so `Ψ_{3,ν}` is the bare light-ray commutator `⟨T₄[L T₃, L T₁]T₂⟩`, contracted with fixed dual and celestial
  structures.
* Its `O(1/c_T)` light side is `λ²_{TTT}Ψ₃[G_T]`: double traces drop at `O(1/c_T)` by Q2 (unsubtracted), and the identity drops
  because `L[T]|Ω⟩ = 0`.
* **F [D; DV in LIGHT_INPUTS `t2`].** The s-channel T-family term of `⟨T₄|L₃L₁|T₂⟩` factorises through `⟨T|L[T]|T⟩`, i.e. the
  one-point shock matrices `H(n)`. So `Ψ₃[G_T] ∝ [H(n₁), H(n₃)]` (contracted).
* With `H = a𝟙 + b(e^{i(4φ+θ)}σ₊ + h.c.)` (Q8), the commutator is `2ib² sin(4φ₁₃)σ_z`, independent of `a`.
  * Check (e) of `v1` [DV]: `∂_a[H₁, H₃] = 0`.
  * So the Einstein term (`a²`) and the linear terms (`a·b`) both vanish. Only `|flip|² ∝ t₄² + α_T²` survives.
* The CHMRS-type de-shift is linear in `Ψ_{3,ν′}`, so it preserves the zero.
* This route needs neither Lemma R nor any Regge counting.

**Route B (Regge counting; [P given Lemma R + "no Π₂ moment"]).**
* The leading heavy action of `Ĉ₃` scales as `m^{−6}` at fixed `β` (§R.3, `ε₃m⁶Ĉ₃ → 2P`, [DV at the crystal]).
* By CHMRS `sec:reggemomentsheavy`, a `Π₂` component `a₂(η)` would contribute `γ²â₂(ν)P(cosh β)/m⁴` at fixed `β`. The heavy map `a₂ ↦ â₂(ν)P(cosh β)` is
  injective (harmonic analysis on `H₂`). So `a₂ ≡ 0`, and `Ĉ₃ ∼ Π₃ + …`.
* With `n = 2`, Lemma R then gives `Ĉ₃|_E = 0`.
* **Caveat.** CLS-D's convergence condition for odd `j′` (`k − 1 > J₀`) would suggest spin-2 decay, i.e. a `Π₂` term. If that
  condition is tight, route B fails and route A carries (L3) alone. The two statements are consistent if CLS-D's odd-`j′` count is
  only sufficient. This was not resolved here; it is noted as a consistency question for a referee.

**Grade of (L3): [P given Q2, Q8, F]**, with F [D] / [DV]. Route B is supporting.

---

## 5. (L2): `−Ĉ₂|_E = 16πG/(ν²+9/4)`

### 5.1 The scalar case, rederived (consistency of the mechanism)

The following steps reproduce CHMRS `BknuGravity`, `−C_{2,ν}|_grav = 8πG/(ν² + d²/4)`. CHMRS only checked it numerically; their
footnote suggests a partial-wave route.

1. Q1 with `k = 2`: `−C_{2,ν}|_E = (a_{Δφ}/γ_{2Δφ+1}(ν)²) · (2πi/4) · ∫[dη]P_{−½+iν}(η)[𝒢^E_+]_{r^{−1}}(η)`.
   * The `a/γ²` is CHMRS's own normalisation: `C_{2,ν} = (a/γ²)Π̂_{2,ν} + O(Π̂₄)`.
2. Q3: the harmonic transform of the leading Regge coefficient is `α(ν) = γ_{CHP}(ν)γ_{CHP}(−ν)β(ν)`, with `j = 2`.
3. **Identity (check (a) of `v1`, exact, any `d`):** `γ_{2Δφ+k−1}(ν)² = γ_{CHP}(ν)γ_{CHP}(−ν)` iff `k = j`.
   * At `k = 2`, `j = 3` the ratio is `ν`-dependent.
   * So the `Δφ`-dependent `Γ`'s cancel, leaving `−C_{2,ν}|_E ∝ β(ν)`, the impact-parameter (phase-shift) partial amplitude.
4. Q4: `β(ν)` is the harmonic transform of the shock profile, i.e. of the `H_{d−1}` propagator of dimension `d−1`.
   * Check (b) of `v1`: `(□_{H_{d−1}} − (d−1))K_{d−1} = 0` for `d = 3, 4, 5`, and KKSDZ's profile `= z·K_{d−1}`.
   * Its harmonic transform is `1/(ν² + (d/2)²)`. Check (c) of `v1`: `∫₁^∞Q_{Δ₀−1}P_{−½+iν} = 1/(ν²+(Δ₀−½)²)` to `≤ 10⁻¹¹`,
     `Δ₀ = 1, 2, 3`, five `ν` each.
   * `Δ₀ = 2` gives `9/4`. A double-pole profile gives `(ν²+9/4)^{−2}` (control, also checked).
5. Normalisation: by the flat limit, `−𝒞_{2,u} = 8πG/(−u)` (CHMRS).
6. End to end, check (f) of `v1` [DV, independent of CHMRS's `R_m` constants]: `−C₂/(8πG) = 1/(ν²+9/4)` to `2·10⁻¹³` at
   `Δφ = 1.2, 2.6` and `ν = 0, 3`.
   * Steps 1–5 explain why the result is `Δφ`-independent.

So in the scalar case the exactness of the pole is **[P given Q1, Q3, Q4]**, and the numerics confirm it.

### 5.2 The spinning case

* **E1 [P given Lemma R].** Only the `Π₂` moment of `Ĉ₂` matters (§3):
  `−Ĉ_{2,ν}|_E = (2πi/4)·𝒫₂,ν[𝒦^E_+|_{r^{−1}}]`, where `𝒫₂,ν` is the leading-moment pairing (harmonic projection in `η`, with the
  `q²` contraction).
* **E2: the premise (HN).** The leading-Regge heavy action of `Ĉ₂` on every heavy `⟨TTTT⟩` block is
  `N̂₂(ν)·2P_{−½+iν}(cosh β)/m⁴ × (ρ⁺⁺ + ρ⁺⁻)`, with `N̂₂` depending on `ν` only.
  * This is §R.3, the definition of the heavy normalisation, and it is used on the whole heavy side (Theorem A, the `n` scans).
  * Support: `REGGE_MOMENT_TTTT.md` Prop. §5 [D + DV pieces]; §24.1 (the de-shift reduces to it iff the spinning `Ψ_ν` has the
    scalar magic-block kernel up to a `ν`-only factor) [I].
  * **(HN) is not new to H7.**
* **E3 [Q3, with D2].** The leading Regge term of `𝒢^E` is a superposition of plane-wave pairs `(p, p̄)` with amplitude `ℬ_E`.
  * Conservation at **both** vertices makes `ℬ_E` a transverse tensor at each vertex.
  * **D2 [D].** CHP prove the transversality for `⟨TT𝒪𝒪⟩`. Their argument is vertex-local (`p·∂_{z_i}ℬ = 0` for each conserved
    operator), so it applies to the `(T₁T₂)` and `(T₃T₄)` vertices separately.
  * Einstein = minimal coupling: `ℬ_E = c·S·Π_{Δ₀=2}(L)·𝟙⊗𝟙` (Q4: shock profile, identity transfer matrix; KPZ §3 for the
    exchange diagram).
  * Non-minimal structures `𝒟̂₂, 𝒟̂₃` carry `t₄`, `α_T` and are excluded from "`|_E`".
* **E4: the same map for heavy and light (D3) [D].** `𝒫₂,ν` is linear. Acting on one plane-wave pair with transverse-helicity
  polarisations, it gives `c_k·γ_a(ν)²M_𝔖(ν;a)²·(|p||p̄|)^{−a}·P_{−½+iν}(cosh L(p̂, p̄̂))·δ_{helicity-preserving}`.
  * Here `a = 2Δ_e + k − 1` is a property of the **functional**, not of the state (Q5; `REGGE_MOMENT_TTTT` Stage 2b [DV] per
    plane wave).
  * Heavy blocks are plane waves at a point of `H₂` (`L = β`) with `|p| ~ m`. Their image is (HN).
  * The light term is a superposition over `H₂` with weight `Π₂(L)` and energy dependence `S^{j−1}`.
  * The energy integral is scale-free. Continued in `k`, it has a pole at `k = j = 2` whose residue is the `r^{1−k}` coefficient of
    Lemma R (LIGHT_INPUTS §4.1). It contributes a `ν`-independent factor (phases from `q = −2ip`, `k`-dependent only).
  * The `ν`-dependent factor `γ_a²M_𝔖²` is **the same** for heavy and light, since `k = j` (scalar analogue: check (a)), and
    `M_𝔖` is helicity-blind on the preserving sector (`REGGE_MOMENT_TTTT` §0 item 4).
  * Dividing by `N̂₂`: `−Ĉ_{2,ν}|_E = C·∫_{H₂}Π₂(L)P_{−½+iν}(cosh L) = C/(ν² + 9/4)`, with `C` independent of `ν`.
* **E5 [Q, CLS-D crystal + flat `B₂`].** `C = 16πG` from the flat limit, `−B₂ ∋ 16πG/p²` (`ν = Λp`).
  * The factor 2 relative to the scalar `8πG` is the two helicity classes of `𝟙⊗𝟙`, each `8πG·Π₂` (LIGHT_INPUTS §4.2).

**Grade of (L2): [D]**, given Lemma R ([P given R1–R5]), (HN), Q3–Q5, and the [D] steps D2, D3. It was [D-sketch].

**Why the grade is not [P].**
1. D2 and D3 are this note's derivations. The tensor bookkeeping of D3 (`q²` against `𝟙⊗𝟙` plane waves, phases) is argued, not
   computed.
2. Q3 and Q4 are leading-Regge (eikonal) statements. Exchanging the plane-wave integral with the Regge moment is assumed, as in
   CHMRS.
3. (HN) is itself [D + DV pieces] / [I] (§24.1).

**What the derivation excludes.**
* **A double pole `(ν²+9/4)^{−2}`.** It would need a non-simple bulk pole in `β`. Q3: `β` has only the simple bulk poles, and Q4:
  Einstein = single exchange with an exact shock.
* **A rational mismatch factor `r(ν) = 1 + O(1/ν²)`** (the LIGHT_INPUTS §4.3 worry). It would need longitudinal components, which
  conservation forbids (Q3, D2), or a heavy/light difference in `γ_a²M²`, which D3 rules out for `k = j`.
* **Additive polynomial terms in `ν²`.** Excluded by the flat limit: Einstein has no `s²·t^n` terms.
* **What would survive.** A failure of (HN) itself (shifted harmonics `w_ε P_{ν−iε}`). It would act on heavy and light alike:
  * light: `16πG Σ_ε w_ε/((ν−iε)² + 9/4)`;
  * heavy: `Σ_ε w_ε P_{ν−iε}`.

  H7 is exact **to the same extent as the heavy normalisation**.

---

## 6. Numerical and symbolic checks (`ads4/code/h7i14/v1_lightside_checks.py` → `cache/v1_lightside_checks.txt`; ≈ 1 min, 2 cores)

| check | result |
|---|---|
| (a) `γ_{2Δφ+k−1}(ν)² / [γ_CHP(ν)γ_CHP(−ν)]` at `k = j` (`2, 3, 4`), general `d`, `Δφ` | `= 1` exactly (sympy); `ν`-dependent at `k ≠ j` |
| (b) shock profile `= z·K_{d−1}`, `□K_{d−1} = (d−1)K_{d−1}` | exact, `d = 3, 4, 5` |
| (c) `∫Q_{Δ₀−1}P_{−½+iν} = 1/(ν²+(Δ₀−½)²)`; double-pole control | rel. err `≤ 10⁻¹¹` (`Δ₀ = 1`), `≤ 2·10⁻¹⁶` (`Δ₀ = 2, 3`, control) |
| (d) residue counting `k = 2…5` against Regge spin 2 | `iπc_{−1}/2, 0, 0, 0` |
| (e) `[H(n₁), H(n₃)]`, `d = 3` helicity basis | `2ib² sin(4φ₁₃)σ_z`; `∂_a = 0` |
| (f) scalar end to end (independent route) | `1/(ν²+9/4)` to `2·10⁻¹³`, `Δφ = 1.2, 2.6`, `ν = 0, 3` |
| (g) NF02 `A`, `B_AdS`, shape sensitivity | `A = 3022.3377`, `B_AdS = −6552.6620`. At `L = 98.5`, `c → c·e^{0.1}` shifts `16L_A` by `−283` (asymptotic `−A·0.1 = −302`) |

**Minor finding.** `B_AdS` is `−6552.662` exactly (`NOFORWARD_FUNCTIONAL.md` §53 and check (g)). `RESULT_D.md` and
`H0_SPINNING_FUNCTIONAL.md` quote `−6552.5` (ten places). The difference, `2.5·10⁻⁵` relative, is not load-bearing, but the quoted
value should read `−6552.7`.

---

## 7. I14: convergence and swappability for `j′ = 0` (`Ĉ₂`) and `j′ = 1` (`Ĉ₃`)

### 7.1 Validity of the sum rules (I14a)

`Ĉ_{k,ν}[𝒢] = 0` for the non-perturbative `⟨TTTT⟩` needs absolute convergence of the light-ray / dDisc integral. KKSDZ §4.1: absolute
convergence of the Wightman-function integral gives existence and commutativity.

| region | condition | status |
|---|---|---|
| Regge corner (`r → 0`, fixed angle) | `Ĉ₂`: `k = 2 > J₀`; `Ĉ₃`: `k − 1 = 2 > J₀` (CLS-D, sufficient) | **K1:** `J₀ ≤ 1`, `J₀^{dDisc} ≤ 1` [Q6]. **Margin one unit for both.** `j′ = 0` is **not** beyond the proven bound. |
| `1 ≈ 3` endpoint (`α₁, α₃ → ∞`, `x₁₃² → 0`) | `f₂ = v′−u′` grows like `α₁` (one power). The KKSDZ condition `J₁+J₃ > 2 − Δ₀′` becomes `3 > 2 − Δ₀′` | [D], satisfied (`Δ₀′ ≥ ½`) |
| light-cone edges, first and second sheet | `J₁+J₃ − |Δ₁−Δ₃| > 2 − τ₀′` (Wightman) or `> 2 − τ_*` (dDisc), shifted by one for `f₂` | [D], satisfied (`3 > 2 − τ`, `τ ≥ ½`) |
| uniformity in the angle towards the edges | KKSDZ bound "at fixed `σ`", with the edges treated in §4.2.5 | [Q], plus KKSDZ's loophole |

**Grade: [P given K1 (Q6 incl. loophole) and K2 (the table)].**
* The Regge row is [Q] (CLS-D counting + KKSDZ).
* The endpoint rows are this note's power counting [D].

### 7.2 Swapping (I14b)

**Lemma S.** Let `ω` have kernel `ℋ` in the double-commutator form (R1). Suppose
`∫|ℋ|·sup_{pol}[dc^{RP}]^{1/2}[dc^{RP}]^{1/2} < ∞`, where `dc^{RP}` is the double commutator in the Rindler-reflected configuration and
polarisations. Then `ω[Σ_O p_O G_O] = Σ_O p_O ω[G_O]`, with absolute convergence.

**Proof [D, standard ingredients].**
* Write the double commutator as `⟨Ω|J(B)·A|Ω⟩` with `A = [𝒪₁, 𝒪₂]`, `B` its partner, and `J` the CRT/Rindler reflection.
* Tomita–Takesaki gives `⟨Ω|J(B)P_O A|Ω⟩ = ⟨Δ^{1/4}P_O B*Ω|Δ^{1/4}P_O AΩ⟩`. Here `P_O` (the projector onto the conformal family of
  `O`) commutes with the modular flow `Δ^{it}`, a boost in the conformal group.
* So each family gives a positive-semidefinite form. Cauchy–Schwarz, twice (per family, then over families), gives
  `Σ_O|⟨…P_O…⟩| ≤ dc^{RP}(A,A)^{1/2}dc^{RP}(B,B)^{1/2}`.
* The two diagonal terms are bounded by KKSDZ's Rindler bounds K1. Fubini then gives the swap.
* This is KKSDZ's eq. `finalCSbound` with a projector inserted. For scalars it reduces to the positivity of `dDisc` per block used by
  Caron-Huot and CMRS-D (Q7).
* **K3 [D]:** the step "`P_O` commutes with `Δ^{it}`". It is elementary but not quoted from a source.

**What K2 must cover for swapping (absolute kernel bounds).**
1. **Regge corner.** `|kernel| ~ r^{k−2}` against `r^{1−J₀^{dDisc}}`: integrable iff `J₀^{dDisc} < k` (`Ĉ₂`, `k = 2`), resp. `< 2`
   for `Ĉ₃` on CLS-D's count. Satisfied, with one unit to spare [Q + K1].
2. **Endpoints and edges.** As in §7.1 [D].
3. **`ν`-smearing and the de-shift.**
   * `Ω^{A′}` integrates `Ĉ_{2,ν}` against `ψ₂(ν/Λ)` on `[0, Λ]`, and `Ĉ₂` is the CHMRS-type combination.
   * In the shifted-contour representation (`ν′ = ν ± i`), the individual kernels grow like `η^{+½}` as `η → ∞` (the light-cone edges
     inside the Regge region). The swap must therefore be done with the **total** kernel, not term by term.
   * In the §R.5 / §24.1 Abel–ODE form, the total kernel is the `ψ`-smeared harmonic kernel divided by `η` (the inverse of
     multiplication by `η`). It decays like `η^{−3/2}/(log η)²` (endpoint `ψ₂ ∝ p`).
   * **[D-sketch]:** the decay was not bounded rigorously, and the exact kernel `gP_{l−1} + g̃P_{l+1}` (§24.1) was not checked
     against it.
   * CHMRS themselves only check swappability of `C_{k,ν}` in examples (Q7). This is the one technical gap left in I14.
4. **`Ĉ₃`** needs no inversion (`j = 0`, §R.5). Item 3 does not arise for it at the level of §R.5; the §24.1 combination for `Ĉ₃` is
   [I] there.

**Grade of I14(b): [P given K1, K2, K3]**, with:
* K1 [Q];
* K2 [Q] at the Regge corner, [D] at the endpoints, [D-sketch] at `η → ∞` for de-shifted `Ĉ₂`;
* K3 [D].

**Is anything a new assumption, as opposed to a literature theorem?**
* **Theorems, modulo the stated loophole:**
  * `J₀ ≤ 1` and `J₀^{dDisc} ≤ 1` for spinning correlators (KKSDZ §4.2.4 `sec:reggelimit`, §4.3 `sec:doublecommutatorconvergence`);
  * the convergence conditions (KKSDZ §4.4);
  * microcausality ⇒ vanishing of the convergent light-ray commutator (KKSDZ §4.1);
  * the residue formula (CHMRS `app:omegalight`, under its hypotheses).
* **Derived here [D]:**
  * endpoint power counting with `f₂`;
  * per-family Rindler positivity (K3);
  * Lemma S.
* **Not established:**
  * the `η → ∞` bound for the de-shifted `Ĉ₂` kernel [D-sketch];
  * KKSDZ's own loophole.

  Neither is a new *physical* assumption: both are technical.

### 7.3 Tree level versus non-perturbative

* The Einstein correlator alone is marginal for `Ĉ₂` (`J₀ = 2`). That is the graviton pole; it is not a convergence failure.
* The sum rule is applied to the full CFT correlator (K1), and `𝒢^{EFT}` enters only through Lemma R, where growth is allowed.
* The tree-level truncation (loops excluded, `O(1/c_T²)`) is the standing assumption of RESULT_D, unchanged.

---

## 8. H16 under route A′

| H16 sub-item | status under route A′ |
|---|---|
| I14 | §7: [P given K1–K3], with K2 [D-sketch] at `η → ∞` |
| P1(`j′ = 0`) (light-ray mode map for the project's `Ψ_k`) | **not used**: the route-A′ heavy side is CLS-D's crystal formula itself (`H0` §R, §21–§24) |
| `j_k(ν)` (untracked Jacobian) | **absorbed**: `N̂_k(ν)` is defined by the heavy action (§R.3), and the light side is computed for the same normalised `Ĉ_k` (E2–E4). Any `ν`-only Jacobian cancels |
| Einstein term absent from `Ψ_k` (k ≥ 3) | reduces to (L3) `Ĉ₃|_E = 0`, §4 [P given Q2, Q8, F]. There is no `k ≥ 4` in route A′ |

---

## 9. Proposed rows (for whoever edits `RESULT_D.md`; not edited here)

| # | hypothesis | proposed grade |
|---|---|---|
| H7 | (R) graviton-pole residue exact: `−Ĉ₂|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` | **[D]** (was [D-sketch]): Lemma R [P given R1–R5, R4 D-sketch]; `Ĉ₃|_E = 0` [P given Q2, Q8, F]; `Ĉ₂` pole exact given (HN) + CHP/KKSDZ/KPZ [Q] + D2, D3 [D]; scalar mechanism [P given Q] and [DV]. `H7_I14_LIGHTSIDE.md` |
| H16 | I14 for route A′ (P1, `j_k`, Einstein-absence not needed / absorbed) | **[P given K1–K3]** (was [I]): K1 = KKSDZ Rindler-positivity Regge bound [Q, with their loophole]; K2 kernel bounds [Q/D; D-sketch at `η → ∞` for de-shifted `Ĉ₂`]; K3 per-family Rindler positivity [D]. `j′ = 0` is within the proven bound (needs `J₀ < 2`, have `≤ 1`) |

Also: `B_AdS = −6552.66`, not `−6552.5` (§6).

---

## 10. Files

* `ads4/H7_I14_LIGHTSIDE.md`: this note.
* `ads4/code/h7i14/v1_lightside_checks.py`: checks (a)–(g). Output in `ads4/code/h7i14/cache/v1_lightside_checks.txt`.
* Literature sources (scratchpad, not in the repo): CHMRS and CLS-D LaTeX already present. KKSDZ 1904.05905, CHP 1707.07689,
  KPZ 1705.02934, KS 1805.00098 and Caron-Huot 1703.00278 were downloaded to `…/scratchpad/h7i14_dl/`.

## 11. Hand-over

1. **Highest value: D3 as a computation.** For one transverse-helicity plane-wave pair, compute the `q²`-contracted
   `Π₂`-moment for both the Laplace form (heavy) and the Fourier form (light, `q = −2ip`). Confirm that the ratio is `ν`-independent.
   This is `REGGE_MOMENT_TTTT` Stage 2b with the continuation added. It would move (L2) from [D] to [P given (HN), Q3, Q4].
2. **(HN) for the de-shifted `Ĉ₂`.** Check that the spinning `Ψ_{2,ν}[q²]` has the scalar magic-block Regge kernel `gP_{l−1} + g̃P_{l+1}`
   up to a `ν`-only factor (§24.1 [I]). Extract the leading Regge kernel from the crystal machinery at two `ν` values.
3. **K2 at `η → ∞`.** Bound the total de-shifted `Ĉ₂` kernel (Abel/ODE form, exact `g`) at large `η` uniformly in `r`. Needs
   `O(η^{−1−ε})` against the light-cone `dDisc` growth.
4. **Resolve the `Ĉ₃` Regge-spin question** (§4 route B caveat): is CLS-D's odd-`j′` condition `k − 1 > J₀` tight, i.e. does `Ĉ₃` have
   a `Π₂` moment? (L3) does not depend on it (route A), but Lemma-R-based H17 for `Ĉ₃` does.
5. R4 for spinning tree diagrams (integer Laurent expansion) and KKSDZ's loophole remain as quoted caveats.
6. Fix the `B_AdS = −6552.5` transcription (→ `−6552.7`) when RESULT_D is next edited.
