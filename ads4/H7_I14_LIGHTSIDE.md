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

---

# Round 2 (19 Nov 2026, second derivation agent): D3, (HN), K2 at `η → ∞`, `Ĉ₃`'s `Π₂` moment

*Appended by a second agent. §0–§11 are unchanged; where they differ, §12–§17 supersede them. Not refereed, not committed.*
*New code: `ads4/code/h7i14/v2_d3_planewave.py`, `v2_k2_kernel.py`, `v2_hn_model.py`; output in `ads4/code/h7i14/cache/v2_*.txt`
(about 4 min in total on 2 cores; `v2_d3_planewave.py` takes about 2.5 min).*

## 12. D3: heavy (Laplace) versus light (Fourier) bookkeeping — now a lemma

**Set-up.** Use the CHMRS frame (3.23), `d = 3`, mostly-plus metric on `ℝ^{2,1}`. `V⁺` is the future cone and `H₂` its unit hyperboloid,
with `x̂ = x/|x|`.

* The leading-Regge term of the functional's kernel is homogeneous: `K(x,y) = (|x||y|)^{a−3} 𝒦(x̂,ŷ; z)`, with `a = 2Δ_e + k − 1`
  (`= 7` for `⟨TTTT⟩` at `k = 2`). Here `𝒦` carries the polarisation slots `z`.
* **Heavy side** (Lemma V, `REGGE_MOMENT_TTTT` §3): a heavy block is a plane-wave pair `e^{2p·x+2p̄·y}` times a constant transverse
  helicity tensor. Its leading action is the **Laplace** transform `L(p, p̄)` of `K`, contracted with that tensor.
* **Light side** (CHP eq. `AtoB`): the Regge term is `∫dp dp̄ e^{−2ip·x−2ip̄·y} ℬ(p, p̄)`, with `ℬ` transverse at each vertex.
  * `ℬ` is homogeneous in `p` and in `p̄`. At `k = j` its degree is complementary to that of `K`; this is v1 check (a), and the
    reason the scale integral is the marginal `∫dr/r` that Lemma R evaluates.

**Lemma D3.** **[P]** (elementary). Uses only homogeneity, Fubini and `SO(2,1)` invariance.

* **(i) Per plane wave.** Assume `∫_{H₂}|𝒦(ŷ;z)|(−ŷ·p̂)^{−Re a}dŷ < ∞`; this holds for `a = 7`, since tensor growth is at most `e^{4χ}`. Then,
  for `Re λ > 0`,

  `∫_{V⁺}d³y e^{λy·p}K(y) = Γ(a) λ^{−a}|p|^{−a} ∫_{H₂}𝒦(ŷ)(−ŷ·p̂)^{−a}dŷ = λ^{−a} L_K(p)`.

  * *Proof.* Write `y = Rŷ`, so `d³y = R²dR dŷ`. Fubini applies by absolute convergence. Use
    `∫₀^∞R^{a−1}e^{−λRc}dR = Γ(a)(λc)^{−a}` with `c = −ŷ·p > 0` (principal branch). ∎
  * The Fourier transform (`λ = −i`, with the Wightman `iε`) is the boundary value of this analytic function of `λ`. It is a
    tempered distribution supported by the cone, and for `p` in the open cone it equals `e^{iπa/2}L_K(p)`.
  * The angular integral is **literally the same** in both cases. So the heavy factors `γ_a(ν)²M_𝔖(ν;a)²` of `REGGE_MOMENT_TTTT`
    Lemma L, and the polarisation dependence, are identical for the Laplace and the Fourier plane wave. The only difference is the phase
    `e^{iπa/2}` per pair, which depends on `a` (i.e. `k`) only.
* **(ii) Parseval on the cone (two pairs).** Let `Ĝ(x̂,ŷ)` be the angular part of the light Regge term (the `r^{−1}` coefficient of
  Lemma R, before projection). Then the leading-moment angular pairing of Lemma R / E1 is

  `𝒫 := ∫_{H₂}dŷ ⟨𝒦(x̂₀,ŷ), Ĝ(x̂₀,ŷ)⟩ = e^{iπa} ∫_{H₂}dp̄̂ ⟨ℬ̂(p̂₀,p̄̂), L(p̂₀,p̄̂)⟩`,

  where `L(p̂,p̄̂) = Γ(a)²2^{−2a}∫∫𝒦(x̂,ŷ)(−p̂·x̂)^{−a}(−p̄̂·ŷ)^{−a}` is the **heavy** per-plane-wave action.
  * *Proof.* (i) applied in both pairs, then Fubini on `(H₂)³`. Absolute convergence: `ℬ̂_E ∝ Π₂ ~ e^{−2L}` and `(cosh d)^{−7}` against
    `P_l ~ e^{−d/2}` and tensor growth `≤ e^{4d}`.
  * The integrand is `SO(2,1)`-invariant, so fixing `x̂₀` on one side is equivalent to fixing `p̂₀` on the other (unimodular group,
    invariant measures). ∎
  * At `a = 7`: `e^{iπa} = e^{7iπ} = −1` for either sign of `iε`, so there is no branch ambiguity.

**Consequence for (L2).** Insert Einstein's `ℬ̂_E = c·Π₂(L)·𝟙_{p̄̂}⊗𝟙_{p̂}` (Q3, Q4, D2). `𝟙` contains only helicity-preserving
components. Let `𝔥^h_ν(β)` be the normalised leading-Regge heavy kernel of `Ĉ₂` on helicity class `h ∈ {++, +−}`: the preserving part of
`L`, divided by the heavy normalisation `N̂₂(ν)`. Then

**`−Ĉ_{2,ν}|_E = 8πG Σ_{h∈{++,+−}} T̂_ν[𝔥^h_ν]`, where `T̂_ν[f] := (1/2π)∫_{H₂}Q₁(cosh d) f(cosh d) dŷ`**,

with the overall constant fixed by the flat limit (E5).
* This holds without (HN). Its grade is **[P given Lemma R, Q3, Q4, D2]**.
* `N̂₂(ν)` cancels **identically**, because it is built from the same `L` (D3(i)). The phases `e^{iπa}` and `2πi/4` (Lemma R) and the
  dDisc phases are `ν`-independent.
* **With (HN)** (`𝔥^h_ν = P_{−½+iν}` for both classes): `T̂_ν[P_l] = 1/(ν²+9/4)` (v1 (c); v2 (D)), so `−Ĉ_{2,ν}|_E = 16πG/(ν²+9/4)`.

**What D3 removes from §5.2 "why not [P]".**
* Item 1 (tensor bookkeeping and phases): it is now a lemma, checked numerically below.
* Item 2 (exchange of the plane-wave integral with the moment): at the level used here (angular pairing at fixed `r`, scale handled by
  Lemma R) the exchange is Fubini with absolute convergence, as in (ii).
* "A rational mismatch factor `r(ν)`" can now come **only** from `𝔥 ≠ P_l`, i.e. from a failure of (HN). D3 rules out any heavy/light
  difference in `γ_a²M²`.

**Numerical confirmation (`cache/v2_d3_planewave.txt`).**

| check | content | result |
|---|---|---|
| (A) | `Γ(a)h_a(ν) = 2^{a−1}√π γ_a(ν)`, `h_a = 2π∫₁^∞x^{−a}P_l`. This is the angular core of (B.2), the "ν-only factor" itself. mpmath, 30 digits; `a = 7, 6.5, 9`; `ν = 0…6` | `≤ 2·10⁻³¹` |
| (B) | Spin-2 rows of Lemma L (`D₁², D₁D₂, D₂²`, preserving vs `M⁰P_l`, flip vs `M^{flip}P^{(4)}`), `a = 7`, `ν = 0…3`, `ξ = 0.4, 1.1`. Not integrated numerically before (`REGGE_MOMENT_TTTT` §4.2) | preserving `≤ 1.1·10⁻¹²`; flip ratio `1.0000000000` (imaginary part `≤ 3·10⁻¹²`) |
| (C) | Direct 3-d transform at `λ = e^{−iθ}`, radial integral done numerically; all spin-2 structures, both helicity sectors, `a = 7, 6.5`, `ν = 0.7, 3` | `I(λ)/(λ^{−a}I(1)) − 1 ≤ 4·10⁻¹⁴` (`θ = π/6`), `7·10⁻¹³` (`π/3`), `2·10⁻⁹` (`0.45π`), `3·10⁻⁷` (`0.47π`; quadrature-limited) |
| (D) | Two-pair Parseval with the shock: `Ĝ(t)` built from `ℬ̂ = Q₁(cosh L)` by two explicit `H₂` convolutions (`ν`-independent), then projected on `P_l` | `𝒫/[h_a(ν)²·2π/(ν²+9/4)] − 1 ≤ 4·10⁻¹³` (`a = 7, 6.5`; `ν = 0, 0.7, 1.5, 3`) |
| (E) | Heavy image of the kernel `ηP_l(η)`, direct `H₂` integral vs `[l h_7(ν+i)P_{l−1} + (l+1)h_7(ν−i)P_{l+1}]/(2l+1)` | `≤ 1.2·10⁻¹²` |
| (E′) | `g(ν)` of `l2_scalar_validate.gfac` vs `(2l/(2l+1))(γ₇(ν+i)/γ₇(ν))²`, and `g̃ = ḡ` | `≤ 3·10⁻¹⁷` |

(E) and (E′) show that **the CHMRS magic-block coefficient is exactly the heavy image of multiplication by `η`**: the subtraction factor
`f₂ = v′−u′ ∝ ηr` at leading Regge order. This is used in §13 and §14.

**Grade of the `Ĉ₂` part of H7.** **[P given (HN), Lemma R, Q3, Q4, D2]**; it was [D].
* D2 (transversality at both vertices) is CHP's vertex-local conservation argument applied twice. It stays [D], but it is a one-line
  extension of a quoted proof, not a new assumption.
* Lemma R is [P given R1–R5] (R4 [D-sketch], as before).
* The general formula `−Ĉ₂|_E = 8πG Σ_h T̂_ν[𝔥^h]` is **[P given Lemma R, Q3, Q4, D2]**, with no reference to (HN).

## 13. (HN): what it must deliver, how it could fail, and what was (not) checked

**Exactly what is needed.** By §12, H7's Einstein term is `8πG Σ_h T̂_ν[𝔥^h_ν]` whatever the heavy kernel is. So:
* (HN) in the form `𝔥^{++}_ν = 𝔥^{+−}_ν = P_{−½+iν}` is necessary and sufficient for `16πG/(ν²+9/4)`;
* any failure of (HN) changes heavy and light sides **coherently**: the light side is the `Π₂`-transform of the same kernel that the heavy
  side uses.

**Structure of the spinning kernel.**
* **[P] Magic block = `η`.** (§12 (E).) The literal CLS-D `Ψ_{2,ν}` has the leading-Regge kernel `f·(index-l harmonic)`, with `f ∝ η`.
  For scalars this gives the heavy image `½[gP_{l−1} + g̃P_{l+1}]`, with `g` exactly CHMRS's coefficient. So the de-shift **is** the
  inversion of `η`-multiplication through the heavy map, as §R.5 / §24.1 described.
* **[D] Spinning: three terms at most.** For the spinning `Ψ_{2,ν}[q²]` the kernel is `η × (index-l covariant tensor harmonic)`. This
  is [D]: it assumes that the light-ray structure of `O_{3/2+iν}` gives an `SO(2,1)`-covariant eigen-kernel of index `l` at leading
  Regge order; the rest of the argument is [P].
  * `η` is the matrix element of the 3-dimensional representation of `SO(2,1)`, so `η·(index l)` decomposes into indices
    `l−1, l, l+1`. The `l` term is absent for scalars and for derivative-free transverse structures.
  * Lemma L maps each preserving index-`l′` component to `M⁰(ν′)P_{l′}`.
  * Hence **`𝔥^h ∈ span{P_{l−1}, P_l, P_{l+1}}`**, and (HN) holds iff, after the `ν`-only normalisation, the `l` coefficient vanishes
    and the `l±1` coefficients are `g, g̃` for both helicity classes.
* **[P] `Re` removes conjugation-odd admixtures.**
  * For real `ν`, `conj(P_{l−1}) = P_{l+1}`. The de-shift `𝔇` (inverse of `Y ↦ gY(ν+i) + g̃Y(ν−i)`) maps real data to real data,
    because `g̃ = ḡ`.
  * So an admixture `δ(ν)(gP_{l−1} − g̃P_{l+1})` with real `δ` is purely imaginary before and after `𝔇`, and `Re` removes it.
  * This is the mechanism behind the helicity-odd imaginary `1/m` term of §22.2 (`m·Im D₊₊ = −m·Im D₊₋`, the same for
    `Δ = 400…6400`). In the `δ` language it is a real, helicity-odd, roughly `ν`-independent `δ_h`, which is exactly the harmless case.
  * The §23.5 prescription is therefore not only legitimate but **needed** for (HN).

**A model of failure (`v2_hn_model.py`, `cache/v2_hn_model.txt`) [M].**
* Suppose the kernel were `η × (𝔖⊗𝔖)P_l` with a fixed derivative-free transverse structure `𝔖`. `η` commutes with `𝔖`, so by Lemma L
  exactly `g → g·ρ`, with `ρ(ν) = [M_𝔖(ν+i)/M_𝔖(ν)]²`.
* This is conjugation-symmetric, so `Re` does **not** remove it.
* `|ρ − 1|`, `a = 7`:

  | `𝔖` | `ν = 0` | 1 | 1.5 | 3 | 10 | 25 | 50 |
  |---|---|---|---|---|---|---|---|
  | `D₁` | 0.063 | 0.13 | 0.17 | 0.21 | 0.066 | `6.9·10⁻³` | `9·10⁻⁴` |
  | `D₁²` | 0.27 | 0.46 | 0.51 | 0.46 | 0.12 | `1.3·10⁻²` | `1.9·10⁻³` |
  | `D₂²` | 1.7 | 3.9 | 9.8 | 1.3 | 0.28 | `1.6·10⁻²` | `1.9·10⁻³` |

  Large-`ν` behaviour: `ρ − 1 = O(ν⁻³)` for `D₁², D₂²`.
* **Reading.**
  1. (HN) is a genuine premise. A generic spinning kernel violates it by `O(1)` at `ν ≈ 1–3`, the region that carries the graviton pole
     `1/(ν²+9/4)`.
  2. The violation decays like `ν⁻³`. At `ν ≥ 25`, where all heavy-side crystal numerics live, it is `≲ 1–2 %` of the de-shift
     coefficient. The de-shift itself is an `O(1/m²)` effect at fixed `x`, so the crystal data **cannot** detect it.
  3. A violation decaying like `ν⁻³` changes `∫ψ₂(−Ĉ₂|_E)` by a convergent amount. It shifts **`B_AdS`, not `A`** (the coefficient of
     `log L` comes from `ν ∈ (n, Λ)` uniformly in `log ν`).
  4. In this model the size would be `|δB_AdS| ~ A × O(0.3)` (`~10³`), i.e. `~10–20 %` of `B_AdS`. This is model-dependent [M].

**Was it checked with the crystal code (`spinfunc/n1_pw`, `l*`)?** No. It could not be done in a meaningful way, and no new runs were made.
* The non-perturbative crystal integral fails for `ν ≲ 15` at every `Δ` (§24.2). The curved contour of §25 (n3, read-only for this
  agent) targets `ν = 5–20` at fixed `J ≤ 10`.
* A fixed-`β` Regge limit needs `J ∝ m`, e.g. `J ≈ 100` at `Δ = 400` for `β ≈ 0.5`. That is far outside the tested range of the
  conserved-structure and `P_J` machinery.
* At accessible points (`ν ≥ 20`, `J ≤ 10`, `β ≲ 0.2`) the model deviation is `≤ 10⁻²` of an `O(1/m²)` term, below the RQMC
  errors (§24.3–§24.4).
* The existing evidence is all consistent with (HN) at large `ν`, and silent at small `ν`:
  * `N̂₂` is helicity-independent to `2·10⁻⁶` at `Δ_ref = 2·10⁴` (§23.5(b)), which bounds a helicity-dependent `P_l`-admixture at `β → 0`;
  * the imaginary helicity-odd `1/m` term is removed by `Re` (above);
  * the scalar-`g` de-shift reproduces the exact skeleton at `ν ≥ 20` (§24.3).

**Status of (HN): [I] for `ν ≲ 15`; consistent at `ν ≳ 20` [NE, via §22–§24].**

**Recommended way to close it** (and in my view the right definition):
* Define `Ĉ₂` on `⟨TTTT⟩` with the **spinning** de-shift coefficients: invert the actual preserving heavy kernel
  `h₋P_{l−1} + h₀P_l + h₊P_{l+1}` (a three-term difference operator in `ν`) instead of using the scalar `g`.
* Then (HN) holds by construction, and §12 makes H7 hold exactly. The cost moves to computing `h_{−,0,+}(ν)` at small `ν`. That is the
  leading-Regge kernel of `Ψ_{2,ν}[q²]` in frame (3.23), a Lemma-L-type computation.
* The heavy side needs the same object at `ν ≲ 15` (§24.7 / §25). At `ν ≥ 20` the change is `≤ 1–2 %` of the shift term, i.e. invisible.
* **Caveat:** if `h^{++} ≠ h^{+−}` beyond `Re`-odd parts, no `ν`-only de-shift gives `P_l` for both helicity classes. The exact
  statement is then §12's helicity sum, and the heavy side must use `𝔥^h` too.

## 14. K2 at `η → ∞`: the total de-shifted kernel, with a bound

**Object.** At leading Regge order the total kernel of `∫₀^Λdν ψ₂(ν/Λ)Ĉ_{2,ν}` is
`𝒦_w(η) = ∫₀^Λdν w(ν)P_{−½+iν}(η)`, `η = cosh β`. Two normalisations:
* `w = ψ₂(ν/Λ)`: heavy / Regge-moment normalisation;
* `w = ψ₂(ν/Λ)·a/γ₇(ν)²`: the dDisc-space kernel, CHMRS's `C_{k,ν}`.

By §12 (E), in the CLS-D building-block form this is `η × [η-inverted superposition of P_{l′}]` = `𝒦_w` itself. So "the exact kernel
`gP_{l−1} + g̃P_{l+1}`" of §7.2 item 3 is **identical** to the inverted form: there is nothing separate to check.

**Proposition K [P].** Let `w ∈ C⁴[0,Λ]` with `w(0) = 0` and `w′(Λ) = w″(Λ) = 0`. Put `M₀ = ‖w‖₁`, `M₂ = |w′(0)| + ‖w″‖₁`,
`c₄ = |w‴(0)| + ‖w⁗‖₁`. Then for `β ≥ 2`

`|(π/2)e^{β/2}𝒦_w(β) − w′(0)/β| ≤ c₄/(3β³) + 2M₀(e^{−β/2} + βe^{−β}/2) + (8M₂/β²)(e^{−1} + βe^{−β}/2) + 2.6M₂/(β−1)²`.

**Proof.**
1. Mehler–Dirichlet gives `𝒦_w = (√2/π)∫₀^βΦ(t)(cosh β − cosh t)^{−½}dt`, with `Φ(t) = ∫₀^Λw cos(νt)dν`.
2. Two and four integrations by parts give `|Φ| ≤ M₂/t²` and `|Φ + w′(0)/t²| ≤ c₄/t⁴`. Cosine inversion with `w(0) = 0` gives
   `∫₀^∞Φ = 0`, so `∫₀^βΦ = w′(0)/β + O(c₄/β³)`.
3. Write `(cosh β − cosh t)^{−½} = √2e^{−β/2}κ`. Then `0 ≤ κ − 1 ≤ 2(e^{t−β} + e^{−t−β})` on `[0, β−1]`, and
   `∫_{β−1}^β κ ≤ 2.6`. ∎

**Consequences.**
* `|𝒦_w(η)| ≤ (2/π)η^{−½}[|w′(0)|/log η + O(1/log²η)]`, with `w′(0) = X(0)/(k₀Λ) ≠ 0` for NF02.
* **Correction to §7.2:** the decay is `η^{−½}/log η`, not `/log²η`. The extra `η⁻¹` of §7.2 is the `r`-range `r ≤ r_max = e^{−β}`; with
  it the large-`η` integrand is `η^{−3/2}/log η`.

**Numbers (`cache/v2_k2_kernel.txt`; NF02, `Δ_g = 100` and `300`).**
* The Mehler evaluation agrees with direct `ν`-integration (mpmath) to `≤ 4·10⁻¹⁰` at `β = 1, 3, 6`.
* `(π/2)e^{β/2}𝒦·β/w′(0) = 0.259, 0.792, 0.955, 0.9891, 0.9965` at `β = 2, 8, 32, 128, 400`, approaching 1. The remainder is
  `≈ −1.39 w′(0)/β²`.
* The bound holds at every tested `β`; it is loose by `10³` at `β ≤ 8` and by `≈ 10×` asymptotically.
* The per-`ν` (unsmeared) kernel behaves like `(2/π)e^{−β/2}(β + 2log2)` at `ν = 0` (no `log` gain, still `η^{−½}log η`).
* The individual shifted building blocks `ηP_{l±1}` grow like `η^{3/2}`.

**Corollary (Fubini at the large-`η` Regge corner) [P given K1′].** Take the region `{ρ = re^β ≤ c}`, `c < 1`.
* **Hypothesis K1′:** `|dDisc 𝒢(r,η)| ≤ B r^{1−J₀}` uniformly in `η` there. This is KKSDZ's Rindler bound "at fixed `σ`" made uniform
  in the angle, which KKSDZ treat in §4.2.5, with their loophole. Both radial coordinates are `≤ c`, so the reflection-positive comparison
  correlator is in its convergent `u`-channel OPE region.
* **Conclusion:** the large-`η` part of `∫|ℋ||dDisc 𝒢|` is `≤ 2B∫^∞dη |𝒦_w(η)|(ce^{−β})^{2−J₀}/(2−J₀) < ∞` for `J₀ ≤ 1`.
* **Margin:** a factor `η^{−½}`. It survives even an `η^{½−ε}` growth of `B`, and it holds per `ν` without the smearing.
* The swap must use the **total** kernel. The shifted building blocks are not absolutely integrable (`η^{3/2}·η^{−1}`).
* The dDisc-space weight `a/γ₇(ν)²` (`~e^{πΛ}`) only enlarges `M₀, M₂, c₄`, not the rate.

**Spinning components [D].**
* The same proof applies to each tensor component of the de-shifted spinning kernel, provided its leading-Regge `η`-dependence is a
  `ν`-smearing of index-`l` functions with real `ν`, with smooth coefficients, and polarisation factors normalised as in KKSDZ.
* Given (HN) in the spinning-de-shift form of §13, that is the case.
* The `ρ ∼ 1` light-cone edge (`ρ̄ → 0` at finite `ρ`) is the separate "light-cone edges" row of §7.1, still [D].

**K2 at `η → ∞`: [D-sketch] → [P given K1′] (scalar-type total kernel), [D] (spinning components).**

## 15. Does `Ĉ₃` have a spin-2 Regge moment? (matters for H17 only)

**Answer: no.** Grade [D + NE].

* **Heavy scaling [NE, §23.3 data].** A `Π₂` component of `Ĉ₃` would give a heavy action `∝ m⁻⁴` at fixed `ν` (Q5). `Π₃` gives `m⁻⁶`.
  * §23.3 found `Ĉ₃/Ĉ₂ ∝ x²` with ratio `4.000` (to `≤ 2·10⁻³`) between `Δ = 800` and `1600` at equal `ν` (`J = 4`, both odd structures).
  * Write `Ĉ₃/Ĉ₂ = A₂/B + (A₃/B)m⁻²`. The ratio `4(1 + ε/4)/(1 + ε)` then bounds the `Π₂` share of `Ĉ₃` at `Δ = 1600` by `ε ≤ 7·10⁻⁴`.
* **Kernel [Q, §R table].** The leading kernel is `−r¹P_l` (`f = −1`, no `η`, so no `r⁰` term).
* **Consistency with the flat limit.** `−B₃ ∋ −2πG|ĝ₃|²p⁴ ≠ 0` does not need a `Π₂` moment. The `|ĝ₃|²` exchange has residue-formula Regge
  spin `n = 4` (§18), so `Ĉ₃` sees it through `[𝒢]_{r^{−2}}` (and `r^{−4}`), as §18 already wrote.
* **Hence:**
  * CLS-D's odd-`j′` convergence count `k − 1 > J₀` is sufficient, not tight (§4 caveat resolved in that direction);
  * route B for (L3) becomes [P given Lemma R + "no `Π₂`" (D + NE)], supporting route A;
  * for H17, Lemma R for `Ĉ₃` needs only `k ≥ 3` moments.
* **Caveats.**
  * The `m`-scaling test was made only at the crystal nodes `ν` of §23.3 and `J = 4`. An analytic `â₂(ν)` that is `≤ 10⁻³` there is
    not excluded by the data alone.
  * The `r¹` kernel statement is quoted, not re-derived for spinning externals.

## 16. Regraded table (supersedes §0, §9 where they differ)

| item | §0–§11 | now | remaining |
|---|---|---|---|
| D3 (heavy vs light plane-wave bookkeeping, phases) | [D] | **[P]** (Lemma D3), [DV] to `10⁻¹²–10⁻³¹` | — |
| H7, `Ĉ₂` part, general form `−Ĉ₂|_E = 8πG Σ_h T̂_ν[𝔥^h_ν]` | — | **[P given Lemma R, Q3, Q4, D2]** | Lemma R's R4 [D-sketch]; D2 [D] (vertex-local CHP) |
| H7, `Ĉ₂` part, `= 16πG/(ν²+9/4)` | [D] | **[P given (HN), Lemma R, Q3, Q4, D2]** | (HN) |
| (HN) | [D + DV pieces]/[I] | **[I] at `ν ≲ 15`**; consistent at `ν ≳ 20` [NE]; a generic spinning kernel violates it by `O(1)` at `ν ∼ 1–3` [M] | Spinning leading-Regge kernel of `Ψ_{2,ν}[q²]`, or redefine `Ĉ₂` with the spinning de-shift (§13) |
| H7, `Ĉ₃` part | [P given Q2, Q8, F] | unchanged; route B now [P given Lemma R + no `Π₂` (D + NE)] | F [D] |
| I14(a) | [P given K1, K2] | unchanged | K1 loophole; K2 endpoints [D] |
| I14(b) | [P given K1–K3], K2 [D-sketch] at `η → ∞` | **[P given K1 (with K1′ uniformity), K2, K3]**; K2 at `η → ∞` now [P given K1′] (scalar-type kernel), [D] for spinning components; K2 endpoints/edges [D]; K3 [D] | K1′ = KKSDZ §4.2.5 uniformity + loophole |
| `Ĉ₃` `Π₂` moment | open | **absent** [D + NE] | analytic derivation |

**Proposed RESULT_D rows (not edited here).**
* H7: **[P given (HN) + Lemma R + CHP/KKSDZ/KPZ]** for `Ĉ₂`, with `(HN)` [I] at small `ν`. State the unconditional form
  `8πG Σ_h T̂_ν[𝔥^h]`. Any `(HN)` failure shifts `B_AdS` (model: `≲ 20 %`), **not `A`** (§13).
* H16/I14: [P given K1 (+K1′), K2 (endpoints [D]), K3 [D]].

## 17. Hand-over (supersedes §11)

1. **(HN) is now the only premise between the heavy side and an exact graviton pole.**
   * Either compute the spinning leading-Regge kernel of `Ψ_{2,ν}[q²]` in frame (3.23), i.e. the coefficients `h_{−,0,+}(ν)` of §13:
     expand the CLS-D dual `q²` and the celestial structure at `r → 0`, then apply Lemma L (the transverse table is now fully checked,
     §12 (B)).
   * Or adopt the spinning-de-shift definition of `Ĉ₂` (§13) and compute `h_{−,0,+}` for the heavy side at `ν ≲ 15`.
   * Watch for a helicity-even `P_l` component (`h₀`): it is not removed by `Re`.
2. If (HN) is replaced by computed `𝔥^h`, evaluate `T̂_ν[𝔥^h]` (one-dimensional `Q₁`-transform; `v2_d3_planewave.two_pair` does the
   two-pair version) and redo `L_A`, `B_AdS`. `A` is unchanged if `T̂_ν[𝔥]·(ν²+9/4) → 1` as `ν → ∞`.
3. K1′: write out the uniformity in the angle of KKSDZ's Regge bound on `{ρ, ρ̄ ≤ c}` (their §4.2.5) for spinning normalised correlators.
   That would make I14(b) [P given KKSDZ incl. loophole].
4. `Ĉ₃` `Π₂` moment: a derivation of the spinning `r¹` kernel would upgrade §15 from [D + NE] to [P].
5. Unchanged from §11: R4 (integer Laurent expansion for spinning tree diagrams); the `B_AdS = −6552.5 → −6552.7` transcription.

---

# Round 3 (19 Nov 2026, third derivation agent): the spinning Regge kernel computed; (HN) decided; `Ĉ₂′`

*Appended by a third agent. §0–§17 are unchanged; where they differ, §18–§22 supersede them. Not refereed, not committed.*
*New code: `ads4/code/h7i14/v3_regge_kernel.py` (kernel, heavy image, validations), `v3_hn_analysis.py` (tables, (HN) tests),
`v3_heavy_impact.py` (near-gap chains with the spinning de-shift; read-only imports of `spinfunc/l9_*`, `n1_analysis`). Output in
`ads4/code/h7i14/cache/v3_*.txt` (about 2 min for `v3_hn_analysis.py`, 1 core).*

**Summary.**
* The leading-Regge kernel of `Ψ_{2,ν}[q²]` is computed in closed form (split representation through the celestial point `z₀`, §18.1).
  Its heavy image is **exactly `h₋(ν)P_{l−1} + h₊(ν)P_{l+1}`: `h₀ ≡ 0`, and `++` and `+−` are identical** (§18.3).
* So the only possible (HN) defect is the **ratio `h₋/h₊`**, the spinning analogue of the CHMRS magic-block coefficients `g/g̃`. It differs
  from the scalar one: `O(1)` for `ν ≲ 8`, `≈ −180/ν³` in the local shift coefficient for `ν ≳ 20` (§18.4).
* **(HN) verdict (§19):** not exact for the current `Ĉ₂` (scalar-`g` de-shift). The defect is `O(1)` and, in the implemented forward
  normalisation, singular at `ν₀ = 5.98` for `ν ≲ 8`, and `≤ 0.5 %` of the graviton-pole term for `ν ≥ 15`.
* **Fix (§20):** `Ĉ₂′ := 𝒜_h⁻¹[Ψ_{2,ν}[q²]]`, the de-shift with the computed spinning coefficients. (HN) holds by construction, so H7 is
  exact and **`A = 3022.34`, `B_AdS = −6552.66`, `k₀`, NF02 and the bounds `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` are unchanged**. The heavy
  side changes only through the local shift coefficient (§20.3).

## 18. The leading-Regge kernel of `Ψ_{2,ν}[q²]` and its heavy image

### 18.1 Derivation: split representation through the celestial point [D]

Frame (3.23): `x₁ = −y`, `x₃ = y`, `x₂ = x/x²`, `x₄ = −x/x²` (as `l6_spin_action.geom`). Regge limit `x, y → 0`, `r = |x||y|`.

* **Light-ray side (operators 1, 3).** The two δ-functions of `Ψ⁺` put `x₀` on the circle `x₀ = |y|n̂_θ`, `n̂·y = 0`, `n̂² = 1`, with uniform
  measure `dθ` (rest frame of `y`). The light-ray polarisations are `z₁ = x_{1⁺0} = −|y|(ŷ + n̂)` and `z₃ = x₃₀ = |y|(ŷ − n̂)`: antipodal null
  directions. The celestial `j′ = 0` structure (`δ₁ = δ₃ = 2`, `δ = ½ + iν`) is `[(−z₀·(ŷ+n̂))(−z₀·(ŷ−n̂))]^{−δ/2}` up to a constant
  (`(ŷ+n̂)·(ŷ−n̂) = −2`) and a ν-only phase (all bases have a fixed sign on the region). With `m̂` the direction of `z₀` at `ŷ`:
  `B(ŷ, z₀)[τ₁₃] = (−z₀·ŷ)^{−δ} ∫₀^{2π}dθ |sin(θ − φ_m)|^{−δ} τ₁₃(ŷ + n̂_θ, ŷ − n̂_θ)`.
  The θ-moments are Beta functions; `n̂n̂`-products reduce to `m̂` and `η + ŷŷ` (exact, `v3_regge_kernel.B_monomials`).
* **Dual side (operators 2, 4).** `x₀ = O(|y|)` against `|x₂| = 1/|x|`, so at leading Regge order the dual structure is evaluated at
  `x₀ = 0` (corrections `O(r)`, i.e. `Π₃`, invisible at leading Regge order). With the polarisations of 2, 4 referred to the origin by the
  inversion (`Z₂ ≃ Inv·Z(x, z₂′)`, as in Lemma V), exactly at `x₀ = 0`, `|x| = 1`:
  `V_{0,42} = (z₀·x̂)/2`, `H₂₀ = z₂′·z₀`, `H₄₀ = z₄′·z₀`, `X₂₀ = X₄₀ = 1`, so
  **`q = 2(z₂′·z₀)(z₄′·z₀) − (z₀·x̂)²(z₂′·I(x̂)·z₄′)`**, `I(x̂) = η + 2x̂x̂`, and
  `A(x̂, z₀)[τ₂₄] = (−z₀·x̂)^{−(1−δ)} · 16q²/(z₀·x̂)⁴` (scalar dual with `Δ̃₂ = Δ̃₄ = 0` times `q²/V⁴`, as in `l6`).
* **Kernel.** `𝒦₀(x̂, ŷ) = ∫Dz₀ A(x̂, z₀)B(ŷ, z₀)` (`z₀` on the celestial circle; total `z₀`-degree `−1`). The full leading kernel is
  `(|x||y|)^{a−3} η 𝒦₀` with `f = v′ − u′ = η(1+r²)/(2r)` exactly in this frame.
* **Radial powers.** Relative to the scalar `Δφ = 3` functional, `⟨TTTT⟩` adds `|x|⁴` (`q²/V⁴`, `V ∝ 1/|x|`), `|y|⁴` (`τ₁₃(x₁₀, x₃₀)`,
  degree 4) and `r⁻⁴` (`f` without `1/(u′v′)`). These cancel pair by pair, so **`a = 7` for both pairs**, as assumed in §12. Independent
  confirmation below: the large-ν shift coefficient (§18.4) and the crystal ν-only factor (§18.2, check X).
* **Scalar control.** `q²/V⁴ → 1`, `τ₁₃ → 1` give `B ∝ (−z₀·ŷ)^{−δ}`, `A ∝ (−z₀·x̂)^{−(1−δ)}`, so `𝒦₀ ∝ P_{−½+iν}(η)` (split representation),
  i.e. CHMRS's `½ηP_l`.

**What is assumed (why [D], not [P]).** The leading Regge kernel is read off `Ψ⁺` (the `−(1↔3)` term is the same shape acting on the other
ordering, as for scalars); ν-only constants and phases are dropped (they cancel in every statement below); `Lemma V`'s inversion-referenced
polarisations for 2, 4. Each is the same input as in the scalar case, where the result is exact (check V1).

### 18.2 Heavy image: method and validation [DV]

* Plane-wave pair `e^{p·x + p̄·y}`, constant preserving tensors `τ₂₄ = ε_{h₂}²⊗ε_{−h₂}²` at `p̂`, `τ₁₃ = ε_{h₁}²⊗ε_{−h₁}²` at `p̄̂` (parallel frame
  of Lemma V). Then `L = −η^{μν}∂_{p^μ}∂_{p̄^ν}∫Dz₀ Â₆(p, z₀)B̂₆(p̄, z₀)` (the η-trick of §12 (E), `b = a − 1`).
* Every term of `A`, `B` is `|x|^c(−z₀·x)^{−λ}(ε·x)^i(ε′·x)^j`, whose transform is `(ε·∂_p)^i(ε′·∂_p)^j F_{c,λ}`,
  `F_{c,λ} = K(c,λ)(−z₀·p)^{−λ}(−p²)^{(2λ−c−3)/2}`, `K = 2^{μ−1}√πΓ(1+c/2)Γ((c+3)/2−λ)`, `μ = c+3−λ`. Derivatives by an exact
  truncated Taylor engine. The `z₀` integral is smooth and periodic (trapezoid, converged at 64 nodes to `10⁻¹⁵`).
* `L(β)` at 12–16 values of `β` is fitted to `{P_{l−1}, P_l, P_{l+1}}(cosh β)` (and the f-stripped image to `P_l`); the residual tests
  the structure.

| check | content | result |
|---|---|---|
| V0 | `K(c,λ)` vs `2π∫₁^∞P_{λ−1}t^{−μ}` (mpmath); full transform vs direct `H₂` quadrature | `≤ 3·10⁻¹⁶`; `8·10⁻¹¹` (`10⁻⁴` in a slowly convergent case, quadrature-limited) |
| V1 | scalar functional through the same pipeline: f-stripped image `= c₀P_l`; η-image vs CHMRS `½c₀[gP_{l−1} + g̃P_{l+1}]` | residual `≤ 1.4·10⁻¹⁵`; `h₀/h₋ ≤ 5·10⁻¹⁴`; `h₋/(c₀g/2) − 1 ≤ 4·10⁻¹⁴` (`ν = 0.3, 1, 2.5, 6`) |
| X | forward (`β = 0`) graviton/scalar ratio `(h₋ᵀ+h₊ᵀ)/(h₋ˢᶜ+h₊ˢᶜ)` at `ν = 40…1000`, fit `log R = p log ν + c + a₁/ν + a₂/ν² + a₃/ν³` | `p = 0.0000`, **`a₁ = 22.06i`**, `a₂ = 92 − 6i`. The crystal (H0 §22.1, independent code, bulk-point regime) found the graviton/scalar ν-only factor **`exp(22i/ν + 82/ν² + …)`** |

Check X is an end-to-end test of the spinning ingredients (`q` at `x₀ = 0`, the light-ray polarisations, `a = 7`) against `l6`'s exact
integrand in a different regime. The `a₂` difference (92 vs 82) is within the truncation of both fits (H0 §22.1 extracted it from `K₂`).

### 18.3 Structure of the spinning heavy image [P given §18.1] + [DV]

1. **The f-stripped kernel is `l`-isotypic.** Its preserving heavy image is `c₀ᵀ(ν)P_l` (fit residual `≤ 3·10⁻¹⁵`, all ν).
2. **`h₀ ≡ 0`.** Fit: `|h₀|/(|h₋|+|h₊|) ≤ 10⁻¹³` at every ν in `0.05…300`, both classes.
   * **Proof (given §18.1).** Write `η = −x̂·ŷ = −x̂^μŷ_μ`. Then `L = −∫Dz₀ V^μ(p̂, z₀)W_μ(p̄̂, z₀)`, with `V^μ` the transform of `x̂^μA` contracted
     with `τ₂₄`. A preserving `τ₂₄` has weight 0 at `p̂`, so it is fixed by the stabiliser `SO(2)_{p̂}` and `p̂ ↦ τ₂₄(p̂)` is
     `SO(2,1)`-equivariant. Hence `V^μ` is an equivariant vector function of `(p̂, z₀)`: `V = (−z₀·p̂)^{−λ}[c₁p̂ + c₂z₀/(−z₀·p̂)]`, and the
     same for `W` at `p̄̂` with `λ′ = 1 − λ`. The `z₀` integral of each of the four products is a split-representation integral with exponents
     `(λ, λ′)`, `(λ−1, λ′+1)` or `(λ+1, λ′−1)`. The first comes with `p̂·p̄̂ = −η_β` and gives `ηP_l ∈ span{P_{l±1}}`; the others give
     `P_{l∓1}`. No term gives `P_l`. ∎
   * The argument uses only the scalar celestial point (the celestial "sphere" is `S¹` in `d = 3`) and the weight-0 tensors. Accordingly,
     `h₀ = 0` was also found numerically for `(H₂₀H₄₀)²`, `(V²H₂₄)²`, the `2↔4`-asymmetric `H₂₀²V₄²` and `H₂₀²H₄₀V₄`, and for a light-ray side
     with the `1↔3` symmetry deliberately broken.
   * **So the §13 worry (a helicity-even `P_l` component not removed by `Re`) does not arise.**
3. **Helicity blindness.** `(h₁,h₂) = (+,−)` [`++`] and `(+,+)` [`+−`] (and their reflections) give identical images, to `3·10⁻¹⁶`.
   * Reason [P given §18.1]: `θ → θ+π` swaps `z₁ ↔ z₃` and leaves the weight invariant (the `1↔3` symmetry of `j′ = 0`). This maps `h₁ → −h₁`;
     combined with parity `(h₁,h₂) → (−h₁,−h₂)`, `++ ≡ +−`. It holds for any `(2,4)` structure.
   * Consistent with H0 §22.1 (helicity-independent ν-only factor to `3·10⁻⁸`).

### 18.4 The spinning magic-block coefficients `h_∓(ν)` [DV]

Normalisation-free statements use ratios at the same ν. Two parametrisations (`cache/v3_hn_analysis.txt`, full table there):
* `a_∓ := 2h_∓ᵀ/(c₀ˢᶜ·(g, g̃))`, i.e. `Ψᵀ` heavy `= ½c₀ˢᶜ[g a₋P_{l−1} + g̃ a₊P_{l+1}]` in the normalisation of the scalar `Ψ` (same dual and
  celestial prefactors as `l6`). The scalar has `a₋ = a₊ = 1`.
* `u_T := (h₋ − h₊)/(h₋ + h₊)` (forward-normalised local shift coefficient); scalar `u_sc = (g − g̃)/(g + g̃) = i Im g/Re g`. Under `Re` only
  `σ := Im u` matters (§19.2).

| `ν` | `a₋` | `a₊` | `|h₋/h₊|` | `Re u_T` | `σ_T` | `σ_sc` | `ν³Δσ` |
|---|---|---|---|---|---|---|---|
| 0.05 | −0.76+0.17i | −0.77−0.08i | 1.010 | 0.065 | +3.49 | +8.45 | 0.00 |
| 0.5 | +0.24+1.96i | −0.68−1.16i | 1.467 | 0.351 | −0.89 | +0.69 | −0.2 |
| 1 | +4.53+1.01i | +0.95−2.15i | 1.977 | 0.598 | +0.81 | +0.10 | +0.7 |
| 2 | −1.9−12.2i | +2.7+3.0i | 3.079 | 1.825 | +0.42 | −0.49 | +7.3 |
| 3 | −22.9−6.7i | −3.0+10.8i | 2.133 | 0.362 | −0.04 | −1.10 | +28.8 |
| 4 | −33.0+11.5i | −9.4+17.7i | 1.739 | 0.481 | −0.83 | −2.14 | +84 |
| 5 | −32.1+27.8i | −12.9+23.9i | 1.559 | 1.242 | −1.85 | −5.07 | +402 |
| 6 | −26.1+38.5i | −13.3+29.1i | 1.452 | 4.582 | −1.93 | **+265.9** (pole at `ν₀ = 5.98`) | — |
| 7 | −18.8+44.3i | −11.5+33.0i | 1.377 | 3.739 | +3.03 | +5.91 | −989 |
| 8 | −11.9+46.9i | −8.5+35.6i | 1.322 | 1.347 | +2.66 | +3.19 | −270 |
| 10 | −0.6+47.0i | −1.5+37.7i | 1.245 | 0.413 | +1.632 | +1.766 | −134 |
| 15 | +15.1+39.2i | +12.0+34.8i | 1.142 | 0.116 | +0.8685 | +0.9062 | −127 |
| 20 | +21.9+31.8i | +19.4+29.6i | 1.092 | 0.060 | +0.6119 | +0.6294 | −140 |
| 30 | +27.3+22.4i | +25.9+21.7i | 1.047 | 0.026 | +0.3928 | +0.3986 | −157 |
| 50 | +30.3+13.9i | +29.7+13.7i | 1.018 | 0.0095 | +0.23185 | +0.23321 | −171 |
| 100 | +31.6+7.0i | +31.4+7.0i | 1.005 | 0.0024 | +0.11522 | +0.11540 | −177 |
| 300 | +31.95+2.35i | +31.93+2.34i | 1.0005 | 0.0003 | +0.03834 | +0.03835 | −180 |

Large-ν fits (`ν = 40…300`):
* `νσ_T = 11.5025 − 0.75/ν + 264/ν²` and `νσ_sc = 11.5010 − 0.31/ν + 415/ν²`: the same `23/2`, as CHMRS's `23/4` for `Im g/(g+g̃)`
  (H0 §24.1). This confirms `a = 7` for the spinning kernel.
* `ν³Δσ = −180.3 + 104/ν + …`, with `Δσ := σ_T − σ_sc`.
* `ν²Re u_T = 24.0 − 8.4/ν`. This part is conjugation-odd and removed by `Re` (§19.2).
* `a_∓ → ≈ 32 e^{22i/ν}`: the spinning `Ψ` is, at large ν, `≈ 32×` the scalar one times the crystal's ν-only factor (check X).

## 19. (HN) for the current `Ĉ₂` (scalar inversion kernel `gP_{l−1} + g̃P_{l+1}`): verdict

### 19.1 What (HN) now reduces to [P given §18]

* By §18.3, `Ψᵀ_ν := Ψ_{2,ν}[q²]` has preserving leading-Regge heavy action `h₋(ν)P_{l−1} + h₊(ν)P_{l+1}` on every heavy `⟨TTTT⟩` state
  (times the positive state factors `ρ⁺⁺ + ρ⁺⁻`). There is no `P_l` term and no helicity dependence.
* The current `Ĉ₂` (H0 §24.1) is the scalar-`g` de-shift `𝒜_g⁻¹[Ψᵀ_ν/N_T(ν)]`, with `𝒜_g: Y ↦ gY(ν+i) + g̃Y(ν−i)`. It satisfies (HN)
  iff there is a real-analytic `c(ν)` with `h₋/N_T = g c(ν+i)` and `h₊/N_T = g̃ c(ν−i)`. Then the heavy action is `c(ν)P_l` and `N̂₂ = c`.
* **This depends on the ν-normalisation `N_T` of `Ψᵀ`, which the current definition does not fix.** H0 §24.1 says "the same combination
  for every correlator", but CLS-D's spinning `Ψ` has no canonical prefactor. Any `N_T(ν)` rescaling changes `𝒜_g⁻¹[Ψᵀ/N_T]`, because
  `𝒜_g` does not commute with ν-dependent rescaling.
* For one specific `N_T*` (solve `c(ν+i)/c(ν−i) = (h₋/h₊)(g̃/g)` for `c`, then `N_T* = h₋/(g c(ν+i))`), (HN) holds exactly. That
  functional is `Ĉ₂′` of §20. The question is therefore whether the normalisations actually used are `N_T*`.

### 19.2 The two normalisations in use, and the size of the defect [DV]

**(i) Same prefactors as the scalar `Ψ` (the `l6` integrand: scalar dual `× q²/V⁴`, scalar celestial).**
* Exact criterion: `a₋(μ − i) = a₊(μ + i)` for real μ (then `c(μ) = a₋(μ−i)`). Evaluated by analytic continuation of the closed forms
  (fit residuals `≤ 10⁻¹⁴`).
* `|a₋(μ−i)/a₊(μ+i) − 1| = 0.35, 0.90, 0.15, 0.036, 0.009` at `μ = 0.7, 2, 4, 9, 20`.
* **Fails**, by `O(1)` at `μ ≲ 4` and roughly as `μ^{−1.8}` beyond. In addition `a_∓` themselves range over a factor `≈ 40` and carry the
  phase `e^{22i/ν}` (§18.4), so (i) is not even a sensible normalisation of the flat limit.

**(ii) Forward / heavy-limit normalisation, as implemented** (H0 §22.1 Fix 2 and §23.5(b): `N̂ = ` the `m → ∞` limit at fixed ν and small
`J`, i.e. `β = 0`; then `shift()` of `l9_analysis` with the scalar `g`; then `Re`).
* Per ν: `Re[Ψᵀ/N̂] = ½[(1 + iσ_T)P_{l−1} + (1 − iσ_T)P_{l+1}]`. `Re` removes `Re u_T`, because `P_{l−1} − P_{l+1} = 2i Im P_{l−1}` is
  purely imaginary for real ν.
* The scalar-`g` de-shift inverts `½[(1 + iσ_sc)P_{l−1} + c.c.]`. So the heavy kernel is `𝔥 = P_l − 𝒜̂_g⁻¹[Δσ·Im P_{l−1}]`,
  `Δσ = σ_T − σ_sc`; locally, `𝔥 ≈ P_l − Δσ·Im P_{l−1}`.
* By §12, the coherent light side is
  `−Ĉ₂|_E = 16πG(1 + ε(ν))/(ν² + 9/4)`, with `ε(ν) = Δσ·2ν(ν² + 9/4)/((ν² + 5/4)² + 4ν²)` (first order).

| `ν` | 1 | 3 | 5 | 6 | 7 | 8 | 10 | 15 | 20 | 30 | 50 | 100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ε(ν)` | +0.51 | +0.51 | +1.1 | **−80** (pole of `σ_sc` at `ν₀ = 5.98`) | −0.76 | −0.12 | −0.026 | −4.9·10⁻³ | −1.7·10⁻³ | −3.9·10⁻⁴ | −5.4·10⁻⁵ | −3.5·10⁻⁶ |

* Large ν: `Δσ ≈ −180/ν³`, `ε ≈ −360/ν⁴`.
* For `ν ≲ 8` the defect is `O(1)`, and it is singular at `ν₀`: there `g + g̃ = 0`, while the graviton forward normalisation `h₋ + h₊`
  has no zero.
* This is H0 §24.1's ill-conditioning of the local de-shift, now shown to be a genuine (HN) failure for the graviton, not only a numerical
  problem.

**Verdict on (HN) for the current `Ĉ₂`.**
* `h₀`-type failure (helicity-even `P_l`): **none**, exactly [P given §18.1; DV].
* Helicity dependence: **none**, exactly [P given §18.1; DV].
* Magic-block ratio: **(HN) fails** for both normalisations in use. In the implemented one (ii) the per-ν defect of the graviton-pole term
  is `|ε| ≥ 0.1` at `ν ≲ 8` (singular at `5.98`), `≤ 2.6 %` at `ν ≥ 10`, `≤ 0.5 %` at `ν ≥ 15` and `∝ ν⁻⁴` beyond. **[DV]**, given the
  §18.1 kernel [D].
* Only the special normalisation `N_T*` makes the scalar-`g` de-shift exact, and it is equivalent to §20's `Ĉ₂′`.

**Correction to §13's failure model.** The model there (`g → gρ`, `ρ = [M_𝔖(ν+i)/M_𝔖(ν)]²`) is of the form `𝒜_g[M_𝔖²P_l]`, so with the
scalar normalisation it satisfies (HN) exactly, with `N̂₂ ∝ M_𝔖²`. It fails only in the forward normalisation (ii). The actual kernel is not
of that form (criterion (i) fails), but its size estimate was the right order: `O(1)` at `ν ∼ 1–8`, `ν⁻³` decay of the shift coefficient.

## 20. The corrected functional `Ĉ₂′` and its impact

### 20.1 Definition [definition; properties P given §18 + D]

**`Ĉ′_{2,ν} := 𝒜_h⁻¹[Ψ_{2,ν}[q²]]`**, with `𝒜_h: Y ↦ h₋(ν)Y(ν+i) + h₊(ν)Y(ν−i)` and `h_∓` the spinning coefficients of §18.4.
* There is no `h₀` term (§18.3), so the three-term difference operator of §13/§17 reduces to two terms.
* `Ω^{A′}` keeps `Re`, `ψ₂`, `ψ₃` and `Ĉ₃` unchanged.

Smeared form: `∫ψ₂(ν)Ĉ′_{2,ν}dν = ∫χ(ν)Ψᵀ_ν dν`, with `χ(ν−i)h₋(ν−i) + χ(ν+i)h₊(ν+i) = ψ₂(ν)`. This is the same type of equation as the
scalar `ψ̃`, with `(g, g̃) → (h₋, h₊)`. Pole terms arise if `h_∓` have poles or zeros in `|Im ν| ≤ 1`. The scalar `g` has the `1/ν` pole at
`ν = 0`; `h_∓` inherit it (the two terms degenerate at `l = −½`). The analytic structure of `h_∓` in the strip was not mapped [I].

Properties:
* **Normalisation-free [P]:** `𝒜_{φh}[Y] = φ𝒜_h[Y]`, so any ν-only rescaling of `Ψᵀ` drops out. This removes the ambiguity of §19.1.
  Equivalently, `Ĉ₂′ = 𝒜_g⁻¹[Ψᵀ/N_T*]/c` with the special normalisation `N_T*` of §19.1.
* **(HN) exact by construction** [P given §18]: the heavy leading-Regge action is `P_{−½+iν}(cosh β)·(ρ⁺⁺ + ρ⁺⁻)` with `N̂₂ ≡ 1`, for both
  helicity classes (§18.3). `Re` is then inert at leading Regge order; it still removes the subleading imaginary terms of H0 §22.2.
* **H7 exact:** by §12, `−Ĉ′₂|_E = 8πG Σ_h T̂_ν[P_l] = 16πG/(ν² + 9/4)`, **[P given Lemma R, Q3, Q4, D2 and the §18.1 kernel (D)]**.
  * At the smeared level the pole term is invariant, `∫χ Re[h₋D₊ + h₊D₋] = ∫ψ₂/(ν²+9/4)` with `D_± = 1/((ν±i)² + 9/4)`, by the same
    contour shift as for the scalar `ψ̃` (H0 §R.9). This holds given the analyticity of `χh_∓` in the strip, the [I] above.
* **Flat limit unchanged:** the de-shift is an `O(1/m²)` effect at fixed `x`, and `Ĉ₂′` and the current `Ĉ₂` differ only in its coefficient,
  by `Δσ = O(ν⁻³)`. The crystal / bulk-point action, NF02 and `k₀ = 27.48509817` are untouched [D].
  * `δk(0)` (the `ψ̃ − ψ₂` mass, `−1.5·10⁻³` at `Δ_g = 100`) changes by relative `Δσ/σ ≈ 1.5·10⁻³` at `ν ∼ Λ`, i.e. by `≈ 2·10⁻⁶`
    absolute [M estimate].
* **No shadow coefficient is used** anywhere in §18–§20, so the `PARITY_ODD_ROUTEA.md` §12 branch artefact does not arise.

### 20.2 Light side: `A`, `B_AdS` and the bounds [P given §20.1]

With `Ĉ₂′` the graviton-pole term is exactly `16πG/(ν²+9/4)`, helicity-blind, so `L_A` and every light-side constant are unchanged:
* **`A = 3022.34`, `B_AdS = −6552.66`** (v1 check (g)), divided by `k(0)_{A′}` as before;
* **bounds `√(t₄²+α_T²) ≤ 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` at `Δ_g = 100 / 300 / 1000`, unchanged** (`n` as before, §20.3).

For the record, the current forward-normalised `Ĉ₂` (§19.2 (ii)) with its coherent light side:
* above `ν = 10 / 15 / 20` the defect `ε(ν)` shifts `B_AdS` by `A∫ε ν dν/(ν² + 9/4) = −19 / −4.0 / −1.4` (`0.3 % / 0.06 % / 0.02 %` of
  `B_AdS`, not `A`);
* below `ν ≈ 8` it is `O(1)` and singular at `ν₀`, so that functional has **no well-defined light side**. The §13 model's "`~10–20 %` of
  `B_AdS`" was the right order for this region.
* `Ĉ₂′` removes the issue, rather than this being a correction to be applied to `B_AdS`.

### 20.3 Heavy side [NE + M] (`v3_heavy_impact.py`, `v3_heavy_split.py`; `cache/v3_heavy_*.txt`)

Method: the existing l9/n1 chains (Laplace rows; and at `Δ_g = 100` the direct rows on `x ≥ 0.2`, as H0 §24.5). The only change is the
`k = 2` local de-shift coefficient `σ_sc → σ_T`. `Ĉ₃` keeps the scalar `g` (its shift relation is a separate open item).
* The baseline reproduces H0 §24.5 exactly: direct rows `0.1814 / 0.2037 / 0.2742 / 0.4885` for `J = 4/6/8/10`.
* At `Δ_g = 300` the baseline gives `t* = (−10.5, +8.0)`, not the older `(−5.1, +5.0)`, because of the `J = 10`, `Δ = 200` row added
  since. `n` and `λ_min` are as before.

| run | `Δ_g = 100` (Laplace rows) | `Δ_g = 100` (direct rows, `x ≥ 0.2`) | `Δ_g = 300` | `Δ_g = 1000` |
|---|---|---|---|---|
| baseline (scalar `g`) | `n = 1`, min `λ/N_s = 0.1855` (`J = 4`, `1.04Λ`) | `n = 1`, `0.1814` | `n = 1`, `0.2054` (`+−`, `J = 5`, unaffected) | `n = 1`, `0.2021` (unaffected) |
| `Ĉ₂′` at all nodes | **`n = 1.09`**: `J = 8`, `1.08Λ` gives `−0.006`; `J = 10`: `0.069` | `n = 1`, min **`0.055`** (`J = 8`, `1.08Λ`) | `n = 1`, `0.2054`; `J = 8`: `0.4304 → 0.4332` | `n = 1`, `0.2021`; changes `≤ 10⁻⁴` |
| `Ĉ₂′` only at `ν ≥ 15` | `n = 1`, `0.1868` | `n = 1`, `0.1827` | — | — |
| `Ĉ₂′` only at `ν ≥ 20` | `n = 1`, `0.1861` | `n = 1`, `0.1820` | — | — |

Reading:
1. **On the controlled window (`ν ≥ 15–20`, where the crystal is computed non-perturbatively) the change is `≤ 0.002 N_s`.** This is the
   `ν⁻³` decay of `Δσ`. At `Δ_g = 300 / 1000` every node has `ν ≳ 15 / 50`, and nothing moves beyond `3·10⁻³ / 10⁻⁴`. **`n = 1` at
   `Δ_g = 300, 1000` stands** (same grade as before, [NE + M]).
2. **At `Δ_g = 100` the whole effect comes from the two nodes `x = 0.05, 0.1` (`ν ≈ 4.9, 9.8`).**
   * There `σ_sc = −4.52` (next to the pole at `5.98`) against `σ_T = −1.72`; at `9.8`, `1.84` against `1.69`.
   * The local (first-order) de-shift is not valid at these ν for either functional (H0 §24.1), so neither row is a computation of
     `Ĉ₂′` there.
   * What is established is that **the `x < 0.2` graviton correction at `Δ_g = 100` is sensitive at `O(1)` to the de-shift coefficient**:
     `J = 8` moves from `0.23` to `−0.006`. It must be redone with the exact non-local `𝒜_h⁻¹` and non-perturbative small-ν crystal data
     (the n3 contour).
   * This sharpens H0 §24.8, but does not change its status: `n` at `Δ_g = 100` remains **[I]**.
3. Theorem A (`CERTIFIED_NEARGAP_A.md`) is a statement about the scalar-skeleton model, whose leading Regge action is `2P_l` exactly. It
   does not involve the graviton de-shift and is unchanged. Its relevance to `⟨TTTT⟩` now rests on (HN), which `Ĉ₂′` provides exactly.

## 21. Regraded table (supersedes §16 where they differ)

| item | §16 | now | remaining |
|---|---|---|---|
| leading-Regge kernel of `Ψ_{2,ν}[q²]` | not computed | **closed form** (split representation, §18.1) **[D]**; validated by the scalar CHMRS `g` (`10⁻¹⁴`) and the crystal ν-only factor `exp(22i/ν)` (`a₁ = 22.06i`) **[DV]** | the `Ψ⁺` vs `−(1↔3)` reading and Lemma V's polarisation referral, both as in the scalar case |
| `h₀` (helicity-even `P_l`) | feared | **`≡ 0`** **[P given §18.1]**, [DV] `≤ 10⁻¹³` | — |
| helicity dependence of the Regge kernel | open | **none** **[P given §18.1]**, [DV] `3·10⁻¹⁶` | — |
| `h_∓(ν)` | — | tabulated `ν = 0.05…300` [DV]; `Δσ = −180/ν³`, `ν²Re u_T → 24` | analytic form; zeros and poles in `|Im ν| ≤ 1` [I] |
| (HN) for the current `Ĉ₂` (scalar `g`) | [I] at `ν ≲ 15` | **fails**: not well-posed without `N_T`; in both normalisations in use the defect is `O(1)` at `ν ≲ 8` (singular at `5.98`), `≤ 0.5 %` at `ν ≥ 15`, `∝ ν⁻⁴` in `ε` **[DV given the kernel]** | — |
| (HN) for `Ĉ₂′` | — | **exact by construction** [P given §18.1] | non-local `χ` (as for the scalar `ψ̃`) not computed |
| H7, `Ĉ₂′` | [P given (HN), …] | **[P given Lemma R, Q3, Q4, D2, §18.1 (D)]**: `16πG/(ν²+9/4)` exactly; **`A`, `B_AdS`, bounds unchanged** | R4; D2; the §18.1 [D] |
| `n = 1` at `Δ_g = 300, 1000` under `Ĉ₂′` | [NE + M] | unchanged (`Δλ ≤ 3·10⁻³`) [NE + M] | as before |
| `n` at `Δ_g = 100` | [I] (`x < 0.2`) | **[I], sharper**: with the spinning local coefficient on `x < 0.2`, `J = 8` reaches `−0.006` (Laplace rows) / `+0.055` (direct rows); on `ν ≥ 15` the change is `≤ 0.002` | exact `𝒜_h⁻¹` + n3 small-ν crystal |

## 22. Hand-over (supersedes §17)

1. **Adopt `Ĉ₂′` in H0 §R / §24.1** (not edited here; read-only for this agent).
   * The current text ("the same combination for every correlator") does not define `Ĉ₂` on `⟨TTTT⟩`: it depends on an unstated
     ν-normalisation of `Ψᵀ` (§19.1), and the implemented one violates (HN) at `ν ≲ 8`.
   * Proposed wording: "`Ĉ_{2,ν} := 𝒜_h⁻¹[Ψ_{2,ν}[q²]]` with the spinning magic-block coefficients `h_∓` of `H7_I14_LIGHTSIDE.md` §18.4; on scalars
     `h_∓ ∝ (g, g̃)` and this is CHMRS's `C_{2,ν}`."
   * With it, RESULT_D's H7 row can read [P given Lemma R, Q3, Q4, D2, H7 §18.1], and (HN) is no longer a premise.
2. **Heavy side at `Δ_g = 100`, `x < 0.2`.**
   * Replace `gfac` by `(1 ± u_T)/2` (`cache/v3_uT_table.pkl`, or `v3_heavy_impact.sigma_T`) in any future de-shift.
   * Better: do the de-shift non-locally, i.e. solve for `χ` on the small-ν crystal data of the n3 contour.
   * The `J = 8`, `m = 1.08Λ` state is the one to watch.
3. **Analytic `h_∓`.** The structure of §18.3 (proof of `h₀ = 0`) gives `h_∓` as finite sums of Γ-ratios times polynomials in ν. A closed
   form would make the strip analysis (zeros and poles of `h_∓`, pole terms of `χ`) [P]. Worth doing before computing `χ`.
4. **`Ĉ₃`**: the same split-representation computation with `o_H·q`, `f = −1` and the `j′ = 1` celestial sign would give its own shift
   relation (H0 §24.1's [I]). `v3_regge_kernel` needs only a new `A_monomials` (the `ε`-structure) and the sign.
5. Unchanged from §17: K1′ uniformity; R4; the `B_AdS = −6552.5 → −6552.7` transcription.

**Files (this round).**
* `ads4/code/h7i14/v3_regge_kernel.py`: kernel, closed-form transforms, fits, checks V0 and V1 → `cache/v3_regge_kernel.txt`.
* `ads4/code/h7i14/v3_hn_analysis.py`: helicity and `h₀` tests, `h_∓` table, asymptotics, exact (HN) test, crystal cross-check X →
  `cache/v3_hn_analysis.txt`, `cache/v3_uT_table.pkl`.
* `ads4/code/h7i14/v3_heavy_impact.py`, `v3_heavy_split.py`: near-gap chains with the spinning de-shift → `cache/v3_heavy_impact.txt`,
  `cache/v3_heavy_split.txt` (and `cache/v3_heavy_tmp/`, `cache/v3_sigmaT_memo.pkl`).

## 23. Consolidating `Ĉ₂′` (v4): exact `h_∓`, strip structure, the non-local weight `χ`, parity-odd sector

*Written 19 Nov 2026 by the `Ĉ₂′`-consolidation agent. Unrefereed. New code `ads4/code/h7i14/v4_*.py`, caches `ads4/code/h7i14/cache/v4_*`.
Nothing in `spinfunc/` (n3_*, p*, l9_*, …) was edited; it was only imported.*

### 23.1 Exact `h_∓(ν)` and their analytic structure (`v4_hclosed.py` → `cache/v4_hclosed.txt`, `.pkl`) [P given §18.1] + [DV]

**Method (computer algebra, no fits).** §18.3's equivariance argument is turned into a closed-form computation.
* The gradient of each side's plane-wave transform is `V_μ = (−z₀·p)^{−λ}[c₁p_μ + c₂z₀_μ/(−z₀·p)]`, with `λ = 1−δ` on the (2,4) side and
  `λ′ = δ` on the (1,3) side. Its component outside `span(p, z₀)` vanishes identically (sympy). Numerically it is `≤ 3·10⁻¹⁶`, and the
  coefficients do not depend on the celestial point.
* Then `∫dφ(−z₀·p)^{−a}(−z₀·p̄)^{a−1} = 2πP_{−a}(cosh β)` (checked to `10⁻¹⁶`). This gives
  `h₊ = 2π[c₁c₁′(l+1)/(2l+1) + c₁c₂′]` and `h₋ = 2π[c₁c₁′ l/(2l+1) + c₂c₁′]`.
* `c₁, c₂` come from a single celestial point, `p = (1,0,0)`, `z₀ = (1,0,1)`. Within one side, all monomial constants `K(c,λ)` differ by
  integer shifts of both Γ-arguments, so the ratios reduce to rising factorials.
* Hence **`h_∓ = 2πK_A^{ref}(δ)K_B^{ref}(δ)·R_∓(ν)` with `R_∓` rational, exactly.**
* Checks:
  * the scalar functional through the same code gives `R₋^{sc}/R₊^{sc} = g/g̃` identically (CHMRS's `g`) [P];
  * spinning vs the v3 fits: `a_∓`, `u_T` agree to `2·10⁻¹¹` at all 53 table nodes [DV];
  * the assembled `h_∓` agree with the v3 β-fits to `10⁻¹³`, including at complex ν (`1 − i`, `3 + 0.5i`) [DV].

**Closed forms** (`l = −½ + iν`). In the **CHMRS normalisation** (`Ψᵀ` with the scalar's prefactors, divided by the scalar's ν-only factor,
so that the scalar has kernel exactly `(g, g̃)`), `H_∓ = N_com(ν)·R^T_∓(ν)` with:
* `N_com = Γ(11/4 − iν/2)²Γ(11/4 + iν/2)²/(2Γ(13/4 − iν/2)²Γ(13/4 + iν/2)²)`: even in ν, positive on ℝ, and **free of zeros and poles for
  `|Im ν| < 11/2`**. The scalar is `(g, g̃) = N_com·(R₋^{sc}, R₊^{sc})`, with `R₋^{sc} = −(2ν+i)(2ν+11i)²/(16ν)` and `R₊^{sc}(ν) = R₋^{sc}(−ν)`.
* `R^T₋ = −2(2ν−i)(2ν+i)(2ν+3i)(2ν+5i)·Q₅(ν) / [ν(2ν−7i)²(2ν−3i)²(2ν+7i)²]`
* `R^T₊ = −2(2ν−5i)(2ν−3i)(2ν−i)²(2ν+i)·Q₄(ν) / [ν(2ν−7i)²(2ν+3i)²(2ν+7i)²]`
* `Q₅ = 32ν⁵ + 848iν⁴ − 6160ν³ − 14536iν² + 18042ν + 4185i`, `Q₄ = 16ν⁴ + 288iν³ − 824ν² + 408iν + 4185`.
* So `ρ := h₊/h₋` and `a_∓ = h_∓ᵀ/h_∓^{sc}` are rational, and `a_∓ → 32` at large ν (§18.4).
* `ρ(ν)ρ(−ν) ≠ 1`. Unlike the scalar, `Ψᵀ_{−ν}` is **not** proportional to `Ψᵀ_ν` [P given §18.1]. Only `ν ≥ 0` is used, so nothing depends
  on this, but the evenness of `Ĉ₂′` in ν is a convention, not inherited from `Ψ`.

**Zeros and poles in `|Im ν| ≤ 1`** (all exact roots of the polynomials above; CHMRS normalisation):

| function | poles | zeros |
|---|---|---|
| `h₋` | `ν = 0` (simple) | `±i/2`, `−0.29280i` |
| `h₊` | `ν = 0` (simple) | `+i/2` (double), `−i/2`, `±1.77877 + 0.64888i` |
| `h₋ + h₊` (forward normalisation) | none (the residues at 0 cancel) | `±i/2`, **`±0.72846 − 0.08264i`**, **`±1.90917 − 0.26418i`**, **`±6.34580 − 0.79227i`** |
| `ρ = h₊/h₋` | `−0.29280i` | `+i/2`, `±1.77877 + 0.64888i` |
| scalar `g` / `g̃` / `g + g̃` | `0` / `0` / none | `−i/2` / `+i/2` / **`±5.97913` (real)** |

Just outside the strip: `h₋` has a double pole at `+3i/2` and zeros at `±1.42630 − 1.30525i`, `−3i/2`; `h₊` has a double pole at `−3i/2`.

**Real axis, `ν ∈ (0, 300]`** [P + numerical scan at `10⁻⁴` spacing]:
* `h₋`, `h₊` have no zeros (`|h_∓|ν ≥ 5.81`). Their only singularity is the simple pole at `ν = 0` that `g` also has, with opposite residues;
  the kernel `h₋P_{l−1} + h₊P_{l+1}` is finite at `ν = 0`.
* `h₋ + h₊` has **no real zero**. The relative forward normalisation `|h₋+h₊|/(|h₋|+|h₊|)` has local minima 0.257 (`ν = 0.730`), 0.514
  (1.930) and **0.173 (6.439)**.
* `u_T = (h₋−h₊)/(h₋+h₊)` is bounded on `(0, ∞)` except for `u_T = i s₀/ν + O(1)` at `ν → 0`, with **`s₀ = 465/2606 = 0.1784`**. For
  comparison, `s₀ = 11/26` for the scalar `g` and `½` for `ηP_l`. At large ν, `s → 23/2`. Further, `max|u_T| = 5.77` on `[0.2, 300]`, at
  `ν = 6.44`.

**The old de-shift singularity at `ν₀`.** For the scalar, `g + g̃ ∝ R₋^{sc} + R₊^{sc} = −(4ν² − 143)/4`. So `ν₀ = √143/2 = 5.97913`
**exactly**, and this is the only zero for `|Im ν| < 11/2` [P].
* For `Ĉ₂′` this zero **moves off the real axis to `ν* = 6.34580 − 0.79227i`**. It stays inside the strip, 0.79 below the contour, and
  leaves a shallow minimum (0.173) of the forward normalisation at `ν = 6.44`.
* Two further complex zeros appear closer to the axis, at `0.72846 − 0.08264i` and `1.90917 − 0.26418i`.

**What this means for `Ĉ₂′` [D].**
1. *The functional itself.* `Ĉ₂′ = 𝒜_h⁻¹[Ψᵀ]` never divides by `h₋ + h₊`; it is normalisation-free (§20.1).
   * The forward-normalised objects the chains use, the local de-shift coefficients `(1 ± u_T)/2`, are finite everywhere on `ν > 0`.
   * **Nothing is singular on the real contour**, in contrast to the scalar-`g` de-shift (pole at `5.979`).
   * The near-zeros make the *local* de-shift large, though not singular, at `ν ≈ 0.6–0.9` and `5–7.5` (`|u_T| ≤ 5.8`). That is relevant only
     for the `x = 0.05` node at `Δ_g = 100` (`ν ≈ 4.9`).
2. *Contour-shift (difference-equation) form of the smeared functional*, `χ(ν−i)h₋(ν−i) + χ(ν+i)h₊(ν+i) = ψ₂`, which §20.1 used.
   * The strip contains the pole of `h_∓` at `ν = 0` (on the contour, as for `g`: CHMRS's pole term) and the zeros above. In particular
     `ρ` has a pole at `−0.29280i`, where the scalar has one at `−i/2`.
   * In that form, residue terms or a contour deformation are needed: **yes, the ν-space form needs them**.
   * The Abel form of §23.2 needs **no analyticity at all** and only the real-axis values of `u_T`. It is the form adopted below.
3. *Index.* The winding number of `ρ` along `(0, ∞)` (`ρ(0⁺) = −1 → ρ(∞) = +1`) is `Δarg ρ/π = 7` for `Ĉ₂′`, against 3 for the scalar `g`
   and 1 for `ηP_l`.
   * This sets the dimension of the null space of the Abel inversion: `(7 − 1)/2 = 3` (§23.2).

### 23.2 The non-local weight `χ` for `Ĉ₂′` (`v4_chi.py` → `cache/v4_chi_100_300_1000.txt`, `v4_chi_<Δ_g>.pkl`, `v4_chi_null.pkl`)

**Definition (generalising H0 §R.5 / §13).**
* Write `Ω^{A′}_{k=2} = Re∫₀^∞ψ₂Ĉ′₂ dν = Re∫₀^∞X(ν)Ψ̂_ν dν`, with `Ψ̂ := Ψᵀ/(h₋+h₊)` (the forward-normalised spinning functional, kernel
  `𝒦_ν = ½[(1+u_T)P_{l−1} + (1−u_T)P_{l+1}]`) and `X := χ(h₋+h₊)` (normalisation-free). `X` is the analogue of h1's `ψ̃₂`.
* Defining identity: `∫X𝒦_ν(η)dν = ∫ψ₂P_l(η)dν` for all `η ≥ 1`.
* Mehler–Dirichlet holds for complex index: `P_{−½+i(ν±i)}(cosh ξ) = (√2/π)∫₀^ξcos((ν±i)t)dt/√(cosh ξ − cosh t)`. The Abel transform is
  injective. So the identity is **equivalent** to
  **`E[X](t) := cosh t·C[X](t) − i sinh t·S[u_T X](t) = Φ(t) := C[ψ₂](t)` for all `t ≥ 0`** [D],
  where `C`, `S` are the cosine and sine transforms on `[0, ∞)`.
* No analyticity of `X`, `ψ₂` or `u_T` enters.
* For `u = i/(2ν)` (`ηP_l`) this is exactly h1's ODE. For `u = ic/ν` it is `H′ + c tanh t·H = Φ/cosh t` with `H = S[X/ν]`.
* Equivalently, with `f(ν) = χh₋(ν)` for `ν > 0` and `f(−ν) = ρ(ν)f(ν)`, it reads `E(t) = ½[e^t f̂(t) + e^{−t}f̂(−t)]` (a Riemann–Hilbert
  type problem on the half line; hence the index of §23.1).

**Numerics [NE].**
* Grids: midpoint ν-grid with the DCT-II/DST-II-paired t-grid; collocation of `E[X]/cosh t = Φ/cosh t`.
* Solver: LSQR (exact adjoint) for the minimum-norm correction `X − ψ₂`. GMRES with the exact ODE preconditioner is used for the
  `ηP` check.
* Grid sizes: `dν = 0.01 / 0.02 / 0.04`, `ν_max = 164 / 328 / 1311` at `Δ_g = 100 / 300 / 1000`.
* Refinement: halving `dν` changes `X` by `≤ 7·10⁻⁹` (relative to `max ψ₂`) on `ν ≥ 20`; doubling `ν_max` changes it by `≤ 4·10⁻¹¹`.

**Validation on the old kernel (`ηP_l`)** [DV]:
* h1's `ψ̃₂` is reproduced to `5.5·10⁻⁶ / 2.1·10⁻⁶ / 7.3·10⁻⁷` (rel. `max ψ₂`).
* `δk(0) = −1.507423·10⁻³ / −1.637649·10⁻⁴ / −1.463195·10⁻⁵`, identical to h1 in all printed digits.
* V3 holds to `≤ 3·10⁻¹³`. GMRES converges in 2–3 iterations; LSQR agrees to `10⁻¹²`.

**Uniqueness: a 3-dimensional null space [NE + D-sketch].**
* A dense SVD of the collocated operator (it depends on `u_T` only, not on `Λ` or `ψ₂`) has exactly three singular values that go to 0
  under refinement. At `dν = 0.05 → 0.025` they are `8·10⁻⁶ → 1.6·10⁻⁷`, `7·10⁻¹¹ → 6·10⁻¹⁴` and `6·10⁻¹⁶ → 1·10⁻¹⁶`; all others are
  `≥ 0.65`.
* The count is the same in the `χ`-normalisation (`n = |h₋|+|h₊|`), and is invariant under `ν → ν/λ`, `λ = 0.5…6` (it drops to 2 at
  `λ = 15` only because one vector then leaves the `ν`-box).
* Rule `(Δarg ρ/π − 1)/2`, checked on three kernels: `ηP_l` → 0, scalar `g` → **1**, `Ĉ₂′` → 3.
* The RHS is consistent: its overlap with the left null vectors goes to 0 with the singular values, so solutions exist.
* The null vectors `X_h` are converged in shape. They peak at `ν = 4.36`, `1.49` and `0.36`, and **all three are supported on `ν ≲ 18`**
  (`|X_h| < 10⁻³ max` beyond; `≤ 2·10⁻⁴` at `ν = 20`; `10⁻¹⁵` at `ν = 50`).
* What they are:
  * By construction, `∫X_h𝒦_ν(η)dν = 0`. So they are **invisible to every term with an η-density**: the leading-Regge heavy action of
    every state and the graviton-pole light term [P].
  * They are **gauge for the full functional** (`∫X_hΨ̂ ≡ 0`) if, on every state, `Ĉ′_ν` has a Mehler–Fock density decaying faster than
    `η^{−3/2}` (equivalently, `Ĉ′_ν` analytic in `|Im ν| ≤ 1` with mild growth). Then `∫X_hΨ̂ = ∫dη ρ ∫X_h𝒦 = 0` [D-sketch]. Whether this
    holds for the `O(1/m²)` heavy remainder at `ν ≲ 18` is **[I]**.
* **Consequence:** the leading-Regge identity fixes `χ` only up to `span{X_h}`. The ambiguity lives at `ν ≲ 18`:
  * `x ≲ 0.18` at `Δ_g = 100`, i.e. the region of the other agent;
  * `x ≲ 0.06 / 0.018` at `Δ_g = 300 / 1000`.
* Gauge used below: `X − ψ₂ ⊥ span{X_h}` in `l²(dν)`.
* **The scalar-`g` (§24.1) definition has the same issue (a 1-dimensional null space) [NE], which H0 §24.1 does not mention.**

**Results** (`Ĉ₂′`, min-norm gauge; `δk(0) := −2∫(X − ψ₂)p⁶dp` in h1's convention):

| `Δ_g` | V2: residual | **V3: `∫X𝒦[D]` vs `∫ψ₂D`** (without inversion) | `δk(0)` (`ηP`, h1) | `Λ²δk(0)` | `max|X−ψ₂|/max ψ₂`, `ν ≥ 20` (all `ν`, at) | tail `ν > Λ` |
|---|---|---|---|---|---|---|
| 100 | `2.5·10⁻¹¹` | **`2.3·10⁻¹⁵`** (`0.38`) | **`−7.2007·10⁻³`** `+2.3·10⁻⁴i` (`−1.5074·10⁻³`) | 69.9 | `2.7·10⁻²` (0.38 at `ν = 6.1`) | `3·10⁻⁹` |
| 300 | `1.5·10⁻¹⁰` | **`7.9·10⁻¹²`** (`0.29`) | **`−7.661·10⁻⁴`** (`−1.638·10⁻⁴`) | 68.3 | `6.1·10⁻³` (0.14 at 6.1) | `6·10⁻¹²` |
| 1000 | `2.3·10⁻¹⁴` | **`1.5·10⁻¹⁵`** (`0.22`) | **`−6.830·10⁻⁵`** (`−1.463·10⁻⁵`) | 68.1 | `1.4·10⁻³` (0.045 at 6.1) | `10⁻¹⁴` |

* **Exact η-inversion identity: holds** [D + NE].
  * The leading-Regge identity holds to the residuals above.
  * The graviton-pole identity V3, `∫X·½[(1+u_T)D(ν+i) + (1−u_T)D(ν−i)] = ∫ψ₂/(ν²+9/4)` (complex LHS; `Im ≤ 2·10⁻¹¹`), holds to
    `≤ 8·10⁻¹²`, against an `O(30 %)` violation without the inversion.
  * V3 is gauge-independent: a null vector normalised to `max ψ₂` moves it only by its own interpolation residual.
  * So **`L_A`, `A = 3022.3`, `B_AdS = −6552.7` are unchanged given I9** (the graviton term is an η-density, `Q₁(η)`).
  * The **pole term needs no contour deformation** in the Abel form.
* **Interior.** `X − ψ₂` agrees with the first-order symbol `½ψ₂″ − (σ̃ψ₂)′`, `σ̃ = −iu_T ≈ 23/(2ν)`, to 0.1–4 % at `p = 0.2…0.8`, e.g.
  `+1.311·10⁻²` vs `+1.286·10⁻²` at `p = 0.4`, `Δ_g = 100`. It is ≈ 23× the `ηP` symbol's `ψ/ν²`-part. This is why `δk(0)` is 4.7× h1's.
* **Flat limit unchanged** [D + NE].
  * `ψ₂`, `ψ₃`, NF02 and `k₀ = 27.48509817` are not touched (the de-shift acts on `Ĉ₂`, not on the weights).
  * `X → ψ₂` as `Λ⁻²`: `Λ²δk(0) = 69.9 / 68.3 / 68.1`, and the interior deviation scales the same way.
  * The bulk-point limit of `Ĉ₂′` is the stripped `B₂` (§20.1).
* **Size of `δk(0)` and what it does.** In h1's convention (`ψ̃` multiplies the forward-normalised `Ψ` whose `|ĝ₃|²` image is taken
  flat), the factors are:
  * `A, B_AdS × 1/(1+δk) = 1.0073 / 1.00077 / 1.00007`;
  * bounds `× (1+δk)^{−½} = 1.0036 / 1.0004 / 1.00003`.
  * With `Ĉ₂′` as the primary object (`Ω = ∫ψ₂ Re Ĉ′`), `k(0) = 1` exactly if `Ĉ′`'s `|ĝ₃|²` image is flat. This shift is then part of
    H17 (`k_AdS − k(0)`; model budget 4.4 % / 0.48 % / 0.043 %), which covers it.
  * At the quoted precision the bounds are unchanged at `Δ_g = 300, 1000`. At 100, `2.21·10⁻⁴` becomes `2.22·10⁻⁴` in h1's convention,
    inside the H17 budget [NE + M].
  * **RESULT_D's `δk(0) = −1.507·10⁻³ …` refers to the `ηP` kernel. For `Ĉ₂′` (and for the scalar-`g` §24.1 functional, which has the same
    `23/(2ν)`) the analogue is `−7.20·10⁻³ / −7.66·10⁻⁴ / −6.83·10⁻⁵`.**
* **Small ν.**
  * `X − ψ₂` is `O(1)` relative at `ν ≲ 10`: the de-shift is genuinely non-perturbative there, and `|u_T|` is large near the
    forward-normalisation minima.
  * `X` dips to `0.11 max ψ₂` at `ν = 6.36`, against `ψ₂ = 0.48`. This mirrors the minimum of `|h₋+h₊|` (`X = χ(h₋+h₊)`, while `χΨᵀ` stays
    smooth).
  * Part of this region is gauge-dependent (null space).
  * Feeding `χ` into the `x < 0.2` heavy side at `Δ_g = 100` was **not** done (the other agent's region).

### 23.3 Parity-odd sector at `Δ_g = 300, 1000` with `Ĉ₂′`, and combined `n` (`v4_odd_impact.py` → `cache/v4_odd_impact.txt`) [NE + M]

Protocol: PARITY_ODD_ROUTEA.md §6 exactly (p4 grid: `K₂` linear chain, `Δ_ref = 2·10⁴`, `a₊₋` from `Ĉ₂`; p6 near gap: `J = 5, 7, 9`
corrected, odd `J ≥ 11` base). The **only change** is in the `k = 2` local de-shift: `g/(g+g̃) → (1+u_T)/2`, with the **exact rational**
`u_T` of §23.1, as v3_heavy_impact does for the even sector. `Ĉ₃` keeps the scalar `g`. No shadow coefficient is touched: the p3
caches are on even-integer rows, and the §12 strip is the identity for state/reference ratios at equal `J`.
* The baseline reproduces `cache/p4_grid.pkl` exactly (`max|ΔΓ| = 0`), and min `F/N_s` `0.2072 / 0.2023`.
* `Γ₂₊₋` changes (max over `x`): `≤ 19 / 67 / 130` c₂ on the `Δ = 200` row (`J = 5/7/9`), `≤ 4 / 10 / 18` on `Δ = 400`, falling as `ν⁻³`.
  On the `Δ = 100` row (only relevant at `Δ_g = 100`) the change is `0.4–1.2·10³`.

| `Δ_g` | run | odd `n = 1`? | min `F/N_s` (state) | odd `t*` | per `J = 5 / 7 / 9` |
|---|---|---|---|---|---|
| 300 | baseline (scalar `g`) | yes | 0.2072 (`+−`, `J = 5`, `1.02Λ`) | (−56.9, `> 10⁴`) | 0.2072 / 0.7100 / 0.4572 |
| 300 | **`Ĉ₂′`** | **yes** | **0.2072** (`J = 5`, `1.02Λ`) | (−54.2, `> 10⁴`) | 0.2072 / 0.7103 / 0.4592 |
| 1000 | baseline | yes | 0.2023 (`J = 5`, `1.01Λ`) | (−210.5, `> 10⁴`) | 0.2023 / 0.3125 / 0.4554 |
| 1000 | **`Ĉ₂′`** | **yes** | **0.2023** (`J = 5`, `1.01Λ`) | (−210.2, `> 10⁴`) | unchanged to `10⁻⁴` |

(`Ĉ₂′` only at `ν ≥ 15`: identical to the above to `≤ 10⁻³`.)

**Combined even + odd** (even sector from `v3_heavy_impact`, recomputed here with the exact `u_T` → identical: `J = 4` `0.2205 / 0.2103`,
`t* = (−10.2, +8.3) / (−36.8, +39.1)`). The sectors are block-diagonal (PARITY_ODD_ROUTEA §3).
* **`n = 1` at `Δ_g = 300` and `1000`.**
* Binding state `+−`, `J = 5` (parity-odd): **min `F/N_s` 0.2072 / 0.2023**, as before.
* Combined `t*` `(−10.2, +8.3) / (−36.8, +39.1)` (set by the even sector).
* Grade unchanged **[NE + M]**; even `J ≥ 12` and odd `J ≥ 11` are uncorrected (base model).

*Addendum to §23.2.* The `δk(0)` integral restricted to `ν ∈ [10, Λ]` equals the full one: `−7.2008·10⁻³ / −7.6611·10⁻⁴ / −6.8297·10⁻⁵`. So
`δk(0)` is gauge-independent and insensitive to `ν < 10`. The first-order symbol gives `−7.063·10⁻³ / −7.645·10⁻⁴ / −6.828·10⁻⁵`
(agreement 2 % / 0.2 % / 0.02 %).

### 23.4 Parity-odd sector at `Δ_g = 100` with `Ĉ₂′`, contour points at `x = 0.1, 0.15`, and combined `n` (`v4_odd100.py` → `cache/v4_odd100.txt`, `v4_odd100_xrows.txt`) [NE + M + I]

(Coordinator's addition, after H0 §25.)

**Protocol.**
* p8_analysis's chain (PARITY_ODD_ROUTEA §11.4 main protocol: the §6 rows, `P8_XROWS = 0`, direct values on the `Δ = 100` row) with:
  * direct p8 values on `x ≥ 0.2`;
  * **the n3 Picard–Lefschetz points** (`spinfunc/cache/n3_odd_pts_{C2,C3H}_J{5,7,9}.pkl`, H0 §25.7) at `x = 0.1, 0.15` (state `Δ = 100` and
    reference `2·10⁴`), merged into p8's direct database and accepted by p8's rule (relative RQMC error ≤ 5 %; all were accepted);
  * `K₂` only at `x = 0.05` (and on the rows `Δ = 200 … 3200`).
* β phase: stripped as usual (`logK → logK − iπ(Δ+J)` on state and reference, which is p8's `P8_PHASE = 1`). Asserted on the entries used:
  `|Im logK − π(Δ+J)| = 0` exactly.
* De-shift: `k = 2` with the exact `u_T` (§23.1); `Ĉ₃` keeps the scalar `g`.
* Near gap: p6 conventions, all odd-`J` states (`J ≥ 11` base).
* `t*`: scale on the whole correction. "x < 0.075 window": the correction on Gauss–Legendre nodes `x < 0.075` scaled by `t` (as n3_analysis
  does for the even sector). Both are exact intervals, since the action is linear in `t`.
* The baseline reproduces §11.4 exactly (`0.2262`, `+−`, `J = 5`, `1.05Λ`).

**Direct vs `K₂` at the new nodes** (`Ĉ₂′`, `Δ = 100` row; `Γ` in c₂ units, `±` RQMC):

| `J` | `Ĉ₂`, `x = 0.10` / `0.15`: direct, `K₂`, ratio | `Ĉ₃`, `x = 0.10` / `0.15`: ratio |
|---|---|---|
| 5 | `−314 ± 83` / `−251 ± 13`; `K₂` `−477 / −319`; 0.66 / 0.79 | 0.94 / 0.94 |
| 7 | `−1061 ± 57` / `−468 ± 7`; `K₂` `−1437 / −575`; 0.74 / 0.81 | 0.98 / 1.00 |
| 9 | `−1491 ± 35` / `+93 ± 7`; `K₂` `−1942 / +49`; 0.77 / (near a zero of the LO shape) | 0.99 / 1.09 |

The pattern is the same as in the even sector (§25.6b): `Ĉ₃` is at 94–100 % of `K₂`, and `Ĉ₂` at 66–81 %, i.e. `K₂` overstates the small-`x`
`Ĉ₂` correction. With the scalar `g` the `Ĉ₂` ratios are 0.83–0.92.

**Results (odd sector, `Δ_g = 100`):**

| run | odd `n = 1`? | min `F/N_s` (state) | per `J = 5 / 7 / 9` | `t*` (whole) | `x < 0.075` window |
|---|---|---|---|---|---|
| scalar `g`, `K₂` below 0.2 (= §11.4) | yes | 0.2262 (`J = 5`, 1.05Λ) | 0.226 / 0.589 / 0.897 | (−153, +5.8) | (−∞, +3.45] |
| scalar `g` + n3 points | yes | 0.2304 (`J = 5`, 1.05Λ) | 0.230 / 0.628 / 0.933 | (−4.8, +9.1) | (−∞, +4.47] |
| `Ĉ₂′`, `K₂` below 0.2 | yes | 0.2349 (`J = 5`, 1.05Λ) | 0.235 / 0.646 / 0.816 | (−23.6, +5.5) | (−∞, +3.20] |
| **`Ĉ₂′` + n3 points (MAIN)** | **yes** | **0.2392** (`+−`, `J = 5`, 1.05Λ) | 0.239 / 0.693 / 0.947 | **(−6.9, +29.1)** | **(−∞, +4.94]** (edge: `J = 9`, 1.09Λ) |

* MAIN with the `x < 0.075` correction scaled by `t_low = −1 / 0 / 2 / 3`: min `0.2405 / 0.2399 / 0.2386 / 0.2379`, always `J = 5`.

**Combined even + odd at `Δ_g = 100`** (even sector, `Ĉ₂′`: coordinator's / H0 §25.6b number from `n3_analysis.py` with `N3_SPIN = 1`:
`n = 1`, min 0.2206 at `J = 4`, 1.04Λ; `x < 0.075` window `[−32.6, +3.62]`; `t* = (−1.1, +3.0)`). The sectors are block-diagonal.
* **Combined `n = 1`. The binding state is even, `J = 4` at 1.04Λ, 0.2206**; the best odd state is `+−`, `J = 5` at 1.05Λ, 0.2392.
* Combined `x < 0.075` window: **`[−32.6, +3.62]`**, set by the even sector (the odd sector alone allows `t ≤ 4.94`).
* Combined whole-correction `t*`: `(−1.1, +3.0)`, also even-set.
* Grade **[NE + M + I]**: controlled down to `x = 0.1`; `x = 0.05` from `K₂` with the local de-shift (now `Ĉ₂′`'s `u_T`, `|u_T| ≈ 2.05` at
  `ν = 4.9`).

**Robustness: the extra `m`-rows (`P8_XROWS = 1`, `m = 1.05/1.1/1.2/1.3Λ`; `cache/v4_odd100_xrows.txt`). Not robust.**
* With these rows, which have `K₂` on all of `x < 0.2` (n3 points exist only on the `Δ = 100` row), the odd sector gives **`n = 0`** in
  every variant:
  * scalar `g`: `J = 9` at 1.3Λ, `−0.532` (§11.6 found `−0.48`);
  * `Ĉ₂′`: `J = 9` at 1.2Λ, `−0.550`.
* The whole failure comes from the `x < 0.075` nodes: with that part set to 0 the minimum is `+0.2450`.
* On those rows the `x = 0.05` node sits at `ν = x·m ≈ 5.4–6.4`. That is the scalar's `ν₀ = 5.98` and the spinning near-zero
  `6.35 − 0.79i` of §23.1 (`|u_T|` up to 5.8), where both the `K₂` series and the local de-shift are unreliable.
* Its contribution to `J = 9` oscillates in `m`: `+0.85 / −1.36 / +0.29 N_s` at 1.1 / 1.2 / 1.3Λ.
* With these rows the `x < 0.075` window is `[−0.93, +0.60]` (`Ĉ₂′`), which excludes `t = 1`.
* So the odd `n = 1` at `Δ_g = 100` holds in the main protocol with a wide window. Like §11.6, it **depends on the unresolved `x = 0.05`
  node at `m = 1.1–1.3Λ`**, which the main protocol only reaches by interpolation in `1/m` between the `Δ = 100` and `Δ = 200` rows.
* Deciding this needs direct (contour) values at `x = 0.05`, `ν ≈ 5–6.5`, below the region where H0 §25.4 validated the method (`ν ≥ 9`),
  together with the exact non-local `χ` instead of the local de-shift there (§23.2). It is **[I]**.

### 23.5 Regraded table (supersedes §21 where they differ)

| item | §21 | now |
|---|---|---|
| `h_∓(ν)` | tabulated [DV]; strip structure [I] | **exact rational × common Γ-factor** (`h_∓ = N_com R_∓`) **[P given §18.1]**, agreeing with the fits to `2·10⁻¹¹` [DV] |
| zeros and poles in `|Im ν| ≤ 1` | [I] | **mapped [P given §18.1]** (§23.1 table); real axis `(0, 300]`: no zeros of `h₋`, `h₊`, `h₋+h₊`; one simple pole at 0, as for `g` |
| old `ν₀` | singular at 5.98 | scalar: `ν₀ = √143/2` exactly [P]; `Ĉ₂′`: the zero moves to `6.34580 − 0.79227i`, plus `0.72846 − 0.08264i`, `1.90917 − 0.26418i`; **no singularity on the contour** [P given §18.1]; local de-shift large (`|u_T| ≤ 5.8`) near them |
| non-local `χ` for `Ĉ₂′` | not computed | **Abel form [D]**; computed [NE]; graviton-pole identity to `≤ 8·10⁻¹²` (gauge-independent) → `A`, `B_AdS` unchanged given I9 |
| uniqueness of `χ` | — | **3-dim null space** at `ν ≲ 18` [NE; index rule (w−1)/2 checked on 3 kernels]; invisible to every η-density [P]; gauge for the full functional iff `Ĉ′_ν` is Mehler–Fock-representable in `|Im ν| ≤ 1` **[I]**. Scalar-`g` (§24.1) has a 1-dim one [NE] |
| `δk(0)` (h1 convention) | `−1.507·10⁻³ / −1.638·10⁻⁴ / −1.463·10⁻⁵` (`ηP`) | **`−7.20·10⁻³ / −7.66·10⁻⁴ / −6.83·10⁻⁵`** for `Ĉ₂′` (and `≈` scalar `g`); `→ 0` as `68/Λ²`; within the H17 budget [NE + M] |
| flat limit, NF02, `k₀` | unchanged [D] | unchanged [D]; `X → ψ₂` as `Λ⁻²` [NE] |
| odd `n` at `Δ_g = 300 / 1000` with `Ĉ₂′` | (not done) | **`n = 1`, min 0.2072 / 0.2023** (`+−`, `J = 5`), odd `t* ≥ −54 / −210` [NE + M] |
| combined even + odd at 300 / 1000 | `n = 1` | **`n = 1`**, binding `+−`, `J = 5` 0.2072 / 0.2023; `t*` (−10.2, +8.3) / (−36.8, +39.1) [NE + M] |
| odd `n` at `Δ_g = 100` with `Ĉ₂′` | [I] | main protocol + n3 points: **`n = 1`, min 0.2392** (`J = 5`, 1.05Λ), `x < 0.075` window `(−∞, 4.94]`; with the extra `m`-rows `n = 0` (`J = 9`, 1.2Λ, `−0.55`, entirely from the `K₂` `x = 0.05` node at `ν ≈ 5.9`) → **[I]** |
| combined at `Δ_g = 100` | [I] | main protocol: **`n = 1`, min 0.2206** (even `J = 4`), window `[−32.6, +3.62]` [NE + M + I]; not robust to the extra odd `m`-rows (above) |

### 23.6 Hand-over (supersedes §22 where they differ)

1. **For H0 §R / §24.1 (not edited here).** Adopt `Ĉ₂′` and replace "η-inversion with the exact kernel `gP_{l−1} + g̃P_{l+1}`" by the
   Abel-form equation of §23.2 with `u_T`. Add two notes:
   * The leading-Regge identity determines `χ` only modulo a 3-dimensional space (1-dimensional for scalar `g`) supported at `ν ≲ 18`. The
     weight is unique given the analyticity of `Ĉ′_ν` in `|Im ν| ≤ 1` [I].
   * RESULT_D's `δk(0)` values belong to the `ηP` kernel; for `Ĉ₂′` they are 4.7× larger (§23.2). The headline bounds are unaffected at
     `Δ_g ≥ 300`; at 100, `2.21 → 2.22·10⁻⁴` in h1's convention, inside H17.
2. **`Δ_g = 100`, odd sector.** The open item is the `x = 0.05` node at `m = 1.05–1.3Λ` (`ν ≈ 5–6.5`).
   * It needs contour values there, which is below H0 §25.4's validated range.
   * It also needs the non-local `χ` instead of the local de-shift. `cache/v4_chi_100.pkl` holds `X(ν)` in the min-norm gauge. Use
     `X/ψ₂ − 1` as the ν-dependent correction to the weight only after the gauge question (null space) is settled for the heavy remainder.
   * Watch state: `J = 9`, `m = 1.2Λ`.
3. **Gauge question (new, definitional).** Show, or refute, that `Ĉ′_ν` on a heavy state is analytic in `|Im ν| ≤ 1`. Then the null
   directions annihilate the full functional.
   * A test: evaluate `∫X_hΨ̂_ν[state]dν` with crystal data at `ν ≲ 18`. Doing that requires the small-ν crystal, so it is the other
     agent's region.
4. Unchanged from §22: `Ĉ₃`'s own shift relation (`v3_regge_kernel` + `o_H`); K1′; R4.

**Files (this round).**
* `ads4/code/h7i14/v4_hclosed.py`: exact `h_∓` (sympy), strip zeros and poles, `N_com`, `s₀` → `cache/v4_hclosed.txt`, `cache/v4_hclosed.pkl`
  (`load()`, `numfuncs()`).
* `ads4/code/h7i14/v4_chi.py`: Abel-form `χ` (LSQR; GMRES + ODE preconditioner), null space, winding, checks V0, V2, V3, R →
  `cache/v4_chi_100_300_1000.txt`, `cache/v4_chi_<Δ_g>.pkl`, `cache/v4_chi_null.pkl`, `cache/v4_chi.log`.
* `ads4/code/h7i14/v4_odd_impact.py`: odd sector at 300/1000 with `Ĉ₂′`, plus even re-run with the exact `u_T` → `cache/v4_odd_impact.txt`.
* `ads4/code/h7i14/v4_odd100.py`: odd sector at 100 with n3 contour points and `Ĉ₂′`, windows → `cache/v4_odd100.txt` (`P8_XROWS=1`:
  `cache/v4_odd100_xrows.txt`), `.pkl`.
