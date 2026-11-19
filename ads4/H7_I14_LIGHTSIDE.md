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
