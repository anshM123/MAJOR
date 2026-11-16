# Inputs (E) and (R) of the AdS4 light side: superconvergence, Regge counting, and what is left

Status: written 16 Nov 2026 by a derivation agent. **Revised after `REFEREE_LIGHT_INPUTS.md`** (M1–M3, m1–m6; list in §7).
**Not a referee-approved result.** Scope: inputs (E) and (R) of `ADS4_LIGHT.md` §2.1, §2.3, §4 (and `REFEREE_ADS4_LIGHT.md` M1, m2).
Code: `ads4/code/lightinputs/` (`bash run_all.sh`, ≈ 4 s, single core; logs in `lightinputs/cache/`, git-ignored). No file of another
agent was modified, except one marked reconciliation note in `ADS4_LIGHT.md` §3 (M1).

Marks: **[D]** derived here (from the items named on the line); **[D-sketch]** derivation outlined, steps not all written;
**[DV]** derived and checked by the code; **[Q]** quoted, with sources as follows:
* KKSDZ = Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov arXiv:1904.05905, LaTeX read;
* CLS-D = arXiv:2311.04271, text read; equation numbers below are arXiv (3.xx), with the HTML-text numbers in brackets;
* CHMRS = arXiv:2106.10274, LaTeX read: the appendices "omega light" and "superconvergence";
* MP = Meltzer–Perlmutter arXiv:1712.04861.

**[I]** inferred, not derived. Nothing here is fitted.

## 0. Results

1. **Counting rule [D from Q].**
   * **Statement.** Take a crossing-antisymmetric dispersive functional whose Regge expansion starts at the tensor moment `Π̂_k` (kernel
     `r^{k−2}` in the CHMRS frame, `SPINNING_E2E.md` §3.1). It gives exactly zero (all `ν`, all frame modes) on the tree-level light side of any EFT
     term with Regge intercept `J₀ < k`. The threshold is the same for both parities.
   * **CLS-D off-by-one.** CLS-D's text after (3.80) [(111)] says `k − 1 > J₀` for odd `j′`. That is a slip, and it is confirmed by the referee (§2.3).
     So a tree exchange of Regge spin `J` is annihilated iff `J < k`.
   * **Standard content.** This is the standard superconvergence mechanism (KKSDZ, CHMRS App. "superconvergence", CLS-D).
2. **Polarisations add no powers of `r` [DV, `t1`].**
   * In the CHMRS frame all four polarisations are `r`-independent.
   * Every weight-0 tensor building block is analytic at `r = 0` with a finite limit.
   * So the frame modes `𝒢_q` are **bounded** by `r^{1−J₀}`, with `J₀` equal to KKSDZ's up to an `r`-independent change of frame [D-sketch, §2.1].
3. **The `k = 3` bare commutator kills Einstein exactly [D given Q], but this is not our `C⁻₃` [DV, `t5`].**
   * **The bare functional.** The odd `k = 3` graviton functional of CLS-D is the bare light-ray commutator: `f⁻_{3;2}` is a constant (−1 from (3.76)/(3.77)
     [(107)/(108)]; (3.80) [(111)] prints +1; the sign is immaterial). Its `O(1/c_T)` light side is exactly `λ²_TTT·Ψ₃[G_T]`. In `d = 3` this is
     `∝ [H(n₁),H(n₃)] = −2i|flip|² sin 4φ₁₃ σ_z` [DV, `t2`; reproduced independently by the referee]. So it has no Einstein term and no linear `t₄`/`α_T`
     term; only `t₄² + α_T²` survives.
   * **Relation to KKSDZ.** This is the `d = 3` specialisation, with `α_T` included, of KKSDZ's own planar sum rules (M2, §3.2).
   * **The mode map (I1) fails.** The light-transform kinematics put the `(1,3)` polarisations at antipodal points of the frame circle. There the project's odd kernels
     `D₅ = ε(z₁,z₃,ŷ)` and `𝔄 = D₅D₁` vanish identically. So the `C⁻₃` defined in `ADS4_FUNCTIONAL.md` (pres–pres `D₅⊗D₅`/`𝔄⊗𝔄`) is **not** the leading
     moment of any bare-commutator functional (§3.1).
   * **Consequence.** The clean `C⁻₃` zero of the first version **does not apply to the project's `C⁻₃`**. It applies to the CLS-D functional, which is a
     different AdS uplift of flat `B₃`.
4. **(E) for the project's `C⁻₃`, `C⁺₄`, `C⁻₅` [D given RB].**
   * Via the counting rule, the Einstein light side of a functional with leading moment `k` is the `r^{1−k′}` Laurent coefficients (`k′ ≥ k`) of `𝒢^E_+`.
   * These vanish if `J₀(G^E) < k`, i.e. if every frame mode of the tree Einstein `⟨TTTT⟩` is `O(r^{2−k})`:
     * `C⁻₃` needs **RB₂** (`J₀(G^E) ≤ 2`);
     * `C⁺₄` needs **RB₃** (`J₀(G^E) ≤ 3`);
     * `C⁻₅` needs **RB₄** (`J₀(G^E) ≤ 4`).
   * Status: **[I]**. KKSDZ assert `J₀^planar = 2` for holographic theories, from tree graviton exchange [Q]. No AdS4 Einstein `⟨TTTT⟩` Regge limit was computed.
   * The forward terms (leverage `10¹²`–`10¹³`) need only RB₃, a full unit of margin. The `C⁻₃` term (leverage `~10⁴`) needs RB₂ with **no** margin.
   * Also required: that the completion of `C⁻₃` to a spin-3 dispersive functional exists. This is [I], because `C⁻₃` is not in the CLS-D family.
5. **(R) [reduced, not derived].**
   * At `k = J` the `γ_{2Δ+k−1}(ν)²` divided out by our `C^±` are the Laplace factors of each plane wave and cancel identically [D-sketch, §4.1].
   * Given Lemma IP [I] and the helicity-independent Einstein shock density [Q], `−C⁺₂|_E = 16πG/(ν²+9/4)`.
     Lemma IP: the `r^{−1}` coefficient of `𝒢^E_+` is a superposition of transverse-helicity plane waves of heavy type.
   * Exclusion of a `(ν²+9/4)^{−2}` admixture: **[I]**. It is [D] only given the MP impact-parameter form *and* that failures of IP stay within re-weighted
     transverse plane waves (`t4`); longitudinal pieces were never computed (M3).
6. **`⟨JJJJ⟩` warm-up.**
   * For currents `J₁+J₃−1 = 1`. All CLS-D current functionals are subtracted, including `k = 3` (`n = 2`).
   * Photon exchange (`J = 1`) is killed by all `k ≥ 2`, and graviton exchange by all `k ≥ 3` [D given RB for the Einstein–Maxwell tree, (I)].
   * KKSDZ's Einstein–Maxwell shock matrix `H_ij = δ_ij + a₂(n_in_j − δ_ij/(d−1))` governs `⟨J|EE|J⟩`, i.e. `⟨JTTJ⟩`, **not** `⟨JJJJ⟩` (m3).

Net effect on `ADS4_LIGHT.md` §4:
* (E) stays conditional for all three terms. It is now a precise Regge bound, RB₂ for `C⁻₃` and RB₃/RB₄ for the forward terms, all [I] and quoted from KKSDZ.
* There is a concrete alternative: replace `C⁻₃` in `Ω` by the CLS-D bare-commutator `k = 3`, `j′ = 1` functional. For that term (E) is [D] given only the
  planar chaos bound, and no `[TT]` family enters at `O(1/c_T)`. Its heavy action in the project normalisation was not computed.
* (R) is [D-sketch given IP].
* (P), the polynomiality of the `|ĝ₃|²` and contact light sides, is not addressed.

## 1. What is needed

`Ω` uses `C⁺₂` (all `ν`), `C⁻₃` (all `ν`, `ψ₃(0) = 162`), `∂_{ν²}C⁺₄|₀` and `C⁻₅|₀` (`ADS4_LIGHT.md` §0). At tree level the light side is `ω[G^EFT]`, with
`G^EFT = G^E + t₄G₁ + (t₄²+α_T²)G₂ + contacts`.
* **(E):** `G^E` gives exactly 0 in `C⁻₃`, `C⁺₄` and `C⁻₅`.
* **(R):** `G^E` gives exactly `16πG/(ν²+9/4)` in `C⁺₂`.

## 2. The general mechanism

**2.1 Frame lemma F [DV, `t1`; frame comparison D-sketch].**
* **The frame.** Take `(x₁,x₂,x₃,x₄) = (−y, x/x², y, −x/x²)` with `|x| = 1`, so `x₁,x₃ = ∓rŷ` and `x₂,x₄ = ∓x̂`. The polarisations are
  `z(θ) = ŷ + cos θ e_∥ + sin θ e_⊥` at points 1, 3 and `I(x̂)z′(θ)` at points 2, 4. All of these are `r`-independent.
* **Blocks are finite and analytic.** `t1` evaluates the 12 normalised `V̂_{i,jk}` and the 6 `Ĥ_ij` at `r = 10⁻²…10⁻⁵` and `η = 1.3, 2, 5`.
  * Every block has a finite limit: fitted exponent `|p| ≤ 0.013`, consistent with 0 given the `O(r)` corrections.
  * The corrections go as `B(r) − B(0) ∝ r` or `r²`.
  * `Ĥ₁₃ = −2(z₁·z₃+2)` is exactly `r`-independent (referee's analytic check).
* **Consequence.** The stripped tensor is `u^{−Δ}Σ_I τ_I f_I(u,v)` with `τ_I` polynomial in these blocks. So every frame mode `𝒢_q` is bounded by the growth of
  the scalar-like `u^{−Δ}f_I`. Cancellations can make a mode smaller, never larger.
* **Comparison with KKSDZ's frame.** KKSDZ (J0definition) use a different frame (`x₂ = −z₂/α₂`, …), in which they show the rescaled polarisations stay finite.
  * Both frames are one-parameter families reaching the same conformally invariant limit (`ρ, ρ̄ → 0` at fixed `ρ/ρ̄` on the second sheet).
  * The conformal map between the two configurations at equal cross-ratios depends on `r` only through a smooth function. It acts on polarisations by a
    matrix that is analytic and invertible at `r = 0`.
  * Hence a bound `O(r^{1−J₀})` in one frame implies the same bound in the other. Status [D-sketch]: the map was not written out.
* **What this settles.** Polarisation contractions add no powers of `r` in this frame. The Regge limit is a dilatation, which does not act on polarisations.

**2.2 Spinning carry-over of CHMRS (eftformula) [D given CHMRS App. E hypotheses].**
* **Hypotheses, as in CHMRS's appendix:**
  * (a) `ω` is a `dDisc_s` integral, antisymmetric under `s ↔ t` (CLS-D construct the spinning functionals so, (3.74)–(3.75) [(103)–(105)]);
  * (b) `dDisc` commutes with the kernel;
  * (c) `G^EFT_+` is single-valued near `r = 0`, with integer Laurent powers, `dDisc_s G^EFT → 0` near the `u`-Regge limit, and Euclidean decay.
* **Result.** With Lemma F the contour deformation is done mode by mode at fixed `θ`'s, which gives
  `−ω|_light = (2πi/4) Σ_{k′≥k} Σ_q ∫[dη] K^{(k′)}_q(η) [𝒢^EFT_{+,q}]_{r^{1−k′}}`.
* **Caveat (m5).** The tower `Σ_{k′≥k}` is the *full* functional's Regge expansion. The leading moment alone does not determine the light side of terms
  with `J₀ ≥ k`, such as contacts. Example: for a bare commutator the `T`-block route (§3) shows that no contact enters at `O(1/c_T)`, since contacts
  only shift double-trace data, which have double zeros. In the tower language this is a cancellation among the `k′ ≥ k` moments, which is invisible at
  leading moment. The **vanishing** statements (terms with `J₀ < k`) need no knowledge of the higher moments.

**2.3 Exact threshold [Q + D, `t3`; confirmed by the referee].**
* **KKSDZ.** `[L[O₁],L[O₃]] = 0` on a common null plane if `J₁+J₃ > J₀+1`, with `J₀` defined as in (J0definition). For `T` this means `J₀ < 3`.
* **CLS-D subtraction factors.** The crystal values (3.77) [(108)] give `f^± ~ m^{−2(k−2J_e+1)}`, i.e. a Regge shift of `k−2J_e+1` for both parities.
  Hence the condition is `J₀ < (2J_e−1) + (k−2J_e+1) = k`.
* **The cleanest proof is internal to CLS-D (referee).** They state that the bare commutator converges for `J₁+J₃−1 > J₀`, and `f⁻_{3;2}` is constant. Together
  these give `3 > J₀`, not `2 > J₀`.
* **Graviton functionals (`J_e = 2`):**
  * `k = 2`: `f⁺ = v′−u′`, unsubtracted; marginal for Einstein, which gives the graviton pole.
  * `k = 3`: `f⁻ =` const, the bare commutator.
  * `k ≥ 4`: subtracted.

## 3. The `k = 3` bare commutator, and why it is not our `C⁻₃`

**3.1 Mode map (I1) [D + DV, `t5`].**
* **Antipodal kinematics.** CLS-D's functional applies `L[O₁](x₀,z₁)`, `L[O₃](x₀,z₃)`, and the polarisation of a light-transformed operator equals its null
  direction. In the CHMRS frame, `x₁ = −y`, `x₃ = y` and `x₀ = |y|n̂` (`n̂ ⟂ y`) give `z₁ ∝ ŷ+n̂` and `z₃ ∝ ŷ−n̂`. These are antipodal: `θ₃ = θ₁+π`
  (as in `SPINNING_E2E.md` §3.3, re-derived).
* **What a bare commutator sees.** It sees the `(1,3)` tensor only through `Σ_q 𝒢_{q₁q₃}(−1)^{q₃}e^{i(q₁+q₃)φ}`. The weight is symmetric under `q₁ ↔ q₃` for
  even `M = q₁+q₃` and antisymmetric for odd `M`.
* **Our kernels.** The project's `C⁻` kernels sit in `M = 0`, the pres–pres harmonic. They are `D₅ = ε(z₁,z₃,ŷ)` for currents and `𝔄 = D₅D₁` for stress tensors.
  Both are antisymmetric in `q₁ ↔ q₃` at `M = 0`, and both vanish identically at antipodal polarisations (`t5`: `max|K| ≤ 10⁻¹⁵`, whereas `D₁` and `D₁²` give 1).
* **Conclusion.** At the level of frame modes, our `C⁻₃` is not in the span of bare-commutator leading moments. Not checked: whether conservation identities
  among the modes of conserved correlators could make the two coincide on `⟨TTTT⟩` [I: no mechanism is known].
* **Resolution of M1.** The first version's (I1) is false.
  * `ADS4_LIGHT.md` §0 item 6 / §3 ("`⌈k/2⌉` families enter `C⁻₃` at `O(1/c_T)`", [I]) is **not contradicted** for the project's `C⁻₃`, which is not a
    bare commutator. It remains [I].
  * For the CLS-D bare commutator no family enters (§3.2).
  * A marked note was added to `ADS4_LIGHT.md` §3.

**3.2 The bare commutator itself [D given Q].**
* **Double traces.** In MFT `[T,T]` is a c-number, so `⟨0|[T₁,T₂]|[TT]_{n,ℓ}⟩ = O(γ)`. The double commutator in `Ψ⁺` is then `O(γ²) = O(1/c_T²)` for every
  family and every structure (referee's argument; it replaces the CLS-D quote).
* **Single block.** The `T` block grows like `r^{−1}`, so `Ψ₃[G_T]` converges. The sum rule itself requires `J₀ < 3` for the full planar correlator. That holds by
  the chaos bound `J₀^planar ≤ 2` (MSS; KKSDZ §"convergence in perturbation theory") [Q]. So no assumption on the Einstein *truncation* `G^E` is needed, but a
  Regge assumption on the planar correlator is (m2).
* **Block action.**
  * Put `x₀ = ∞`, where `L[T](∞,(1,n)) = ½E(n)`, and Fourier transform in `x₂, x₄`.
  * `E(n)` conserves momentum.
  * The `T` family at fixed timelike `q` is `{ε·T(q)}`.
  * So the `T` part of `⟨T|[E₁,E₃]|T⟩` is `(q⁰/Ω)² ε†[H₁,H₃]ε`. This holds pointwise in `x₂, x₄, z₁, z₃`, hence for every kernel.
* **`d = 3` [DV, `t2`].** KKSDZ's `H_{ij,kl}(n)` in the helicity basis is `H = [[1, (t₄/4)w^{−4}],[(t₄/4)w⁴, 1]]`, with `w = e^{iφ}`.
  * The `t₂` structure vanishes identically, as it must in `d = 3`.
  * The diagonal is exactly 1, which independently confirms `ADS4_LIGHT.md` §2.3.
  * The referee re-derived it without the KKSDZ tensor, from the HM one-point form, SO(2) covariance and the energy sum rule: `H = 𝟙 + flip` is forced, and
    `α_T` enters only as `Im c`, `c = (−t₄+iα_T)/4`.
  * Result: `[H₁,H₃] = −2i|c|² sin 4φ₁₃ σ_z`, `|c|² = (t₄²+α_T²)/16`.

**3.3 Relation to KKSDZ (M2).** KKSDZ already derive this structure in general `d`, (superconvergenceexamples) and their App.:
* their `{2,2|2|2,2}` sum rules have the `T` contribution `∝ ((4+d)t₂+4t₄)²` (in `d = 4`: `(t₄+2t₂)²`), with no Einstein and no linear term;
* they are "satisfied if all single traces have minimal couplings";
* double traces enter at `O(1/N⁴)`.

§3.2 is their `d = 3` specialisation. In `d = 3` the `t₂` structure is absent, and `α_T` enters through `|c|²`. **New here:**
* the CLS-D off-by-one (an erratum; harmless for CLS-D);
* the explicit `d = 3` two-line commutator including the parity-odd coupling;
* the frame lemma;
* the mode-map finding of §3.1, which is negative for the project;
* the bookkeeping.

## 4. (E) for the project's moments, and (R)

**4.0 Einstein term in `C⁻₃`, `C⁺₄`, `C⁻₅` [D given RB_{k−1}].**
* **Statement.** By §2.2 the Einstein light side of a moment-`k` functional is a combination of the `r^{1−k′}` (`k′ ≥ k`) coefficients of `𝒢^E_{+,q}`. It vanishes if
  `𝒢^E_{+,q} = O(r^{2−k})` with integer Laurent powers. Equivalently it vanishes if `ω[G^E] = 0`, which follows from microcausality of tree Witten diagrams [I] and
  `J₀(G^E) < k`.
* **Required bounds:**
  * `C⁻₃`: RB₂;
  * `C⁺₄`: RB₃;
  * `C⁻₅`: RB₄.
* **Status of "tree Einstein `⟨TTTT⟩` has `J₀ = 2`": [I] (m6).**
  * KKSDZ assert it for holographic theories [Q].
  * Conformal Regge theory: the leading pole of a spin-2 exchange is `σ^{1−j}` times `r`-independent structures (MP (aregdiff)) [Q].
  * In flat space, all GR helicity amplitudes have Regge spin 2 [D].
  * The explicit AdS4 result (Raju 1201.6452; AdS4 spinor-helicity results) was not checked.
* **Two extra caveats for `C⁻₃`:**
  * RB₂ has no margin. An `r^{−2}` term would correspond to Regge spin 3. That is excluded by the planar chaos bound for the *full* correlator, but not
    obviously for the Einstein truncation alone.
  * The existence of a spin-3 dispersive completion of `C⁻₃` is [I], since it is not in the CLS-D family.

**4.1 (R): the `Γ`-factors [D-sketch].**
* **Heavy side.** Each heavy plane wave `e^{2P·x}` with a transverse helicity tensor transforms into `c(a,½±iν)|2P|^{−a}M(ν;a)`, with `c = 2^{a−1}√πγ_a`
  (`SPINNING_E2E.md` §4.4).
* **Light side.** The light Regge term in MP's impact-parameter form uses the same master formula at `q = −2ip`, which only adds `k`-dependent phases. For a Regge term
  `r^{1−J}` the `|p|`-integral is a pure power. So `Π̂_k` has a pole at `k = J`, and its residue is the eftformula value, with `γ` taken at `a = 2Δ+J−1`.
* **Cancellation.** That is exactly what `C_J` divides out. Scalar check: CHMRS `8πG/(ν²+d²/4)`, `Δφ`-independent [DV in `ADS4_LIGHT.md` §1].
* **Not written out.** The phase bookkeeping, and the equality of the residue in `k` with the `r`-contour of §2.2 for spinning modes.

**4.2 Lemma IP and the Einstein density.**
* **Lemma IP [I].** The `r^{−1}` coefficient of `𝒢^E_+` is a superposition of transverse-helicity plane waves of heavy type, with no `z·p̂` components. Conservation
  makes this plausible (MP use `z·p̂ → 0`).
* **Einstein density [Q].** The AdS shock profile `∝ z^d/(z²+y²)^{d−1}` is `Π^{H₂}_{Δ₀=2}`, and the transfer matrix is `H = 𝟙` (KKSDZ (adsshock); CEMZ). So the density is
  helicity-independent: `8πG·Π_{Δ₀=2}` per helicity class.
* **Consequence.** Given IP, `−C⁺₂|_E = 16πG/(ν²+9/4)` with the constant fixed by (F), and `C⁻₂|_E = 0`.

**4.3 Possible form of an (R) violation [D given MP form + transverse class; otherwise I] (`t4`).**
* **Assumptions.** The light density is (polynomial differential operator) × `Π_{Δ₀=2}` (MP form), and failures of IP are re-weightings within transverse
  plane waves.
* **Where singularities can sit.** A mismatch then enters `C⁺₂` through ratios of rational `M`-table entries. Their singularities are at:
  * `ν² = −81/4` and `−25/4`;
  * the roots of `N₇`: `ν² = −413/4 ± √9049` (≈ `−8.1`, `−198.4`).

  None of these is at `−9/4`, so `(ν²+9/4)^{−2}` cannot arise.
* **Not covered.** Longitudinal (`z·p̂ ≠ 0`) pieces, whose tables were never computed, could behave differently. So the exclusion is **[I]** in general (M3).

## 5. Updated input table (proposed; replaces the (E), (R) rows of `ADS4_LIGHT.md` §4 if a referee accepts it)

| input | status |
|---|---|
| counting rule: moment-`k` functional annihilates EFT terms with `J₀ < k` (both parities) | [D] given CHMRS App. E hypotheses + Lemma F [DV; frame comparison D-sketch]; CLS-D off-by-one [D, refereed] |
| (E) project `C⁻₃` (`D₅⊗D₅`/`𝔄⊗𝔄`): no Einstein term | **[D given RB₂]**, RB₂ [I] (KKSDZ [Q]); existence of its spin-3 completion [I]; (I1) **false** [DV, `t5`] |
| (E) CLS-D bare commutator `k = 3`, `j′ = 1` (alternative uplift of flat `B₃`): light side `= λ²Ψ₃[G_T] ∝ t₄²+α_T²`, no `[TT]` at `O(1/c_T)` | [D] given planar `J₀ < 3` [Q] and the HM/CDP one-point matrix [Q]; `d = 3` algebra [DV]; this is KKSDZ's sum rule in `d = 3` |
| (E) `C⁺₄`, `C⁻₅`: no Einstein term | **[D given RB₃, RB₄]**, [I] (KKSDZ [Q]) |
| (R) `R ≡ 1` in `C⁺₂` | [D-sketch given IP [I]] + Einstein eikonal `= 𝟙·Π_{Δ₀=2}` [Q]; `(ν²+9/4)^{−2}` excluded only within MP form + transverse class, otherwise [I] |
| (P) polynomiality of `|ĝ₃|²`/contact light sides | not addressed ([I]) |

**Not done**, in decreasing order of value:
1. The heavy action of the CLS-D `k = 3`, `j′ = 1` functional in the project normalisation, so that it can replace `C⁻₃` in `Ω`.
2. The Regge limit of AdS4 Einstein `⟨TTTT⟩` (RB₂–RB₄).
3. The direct `T`-exchange `C⁺₂`.

## 6. Code (`ads4/code/lightinputs/`, `bash run_all.sh`)

| file | content | log |
|---|---|---|
| `t1_frame_regularity.py` | weight-0 tensor blocks in the CHMRS frame as `r → 0`: finite, analytic (Lemma F) | `cache/t1.log` |
| `t2_shock_commutator.py` | KKSDZ shock matrix in `d = 3`, helicity basis; `[H₁,H₃]`; `t₂` drops out; diagonal = 1 | `cache/t2.log` |
| `t3_clsd_counting.py` | Regge counting of the CLS-D subtraction factors; threshold `J₀ < k` for both parities | `cache/t3.log` |
| `t4_mtable_poles.py` | singular points of the `l = 2` `M`-tables (`k = 2…5`) vs the `T` pole | `cache/t4.log` |
| `t5_mode_map.py` | (I1): `D₅`, `𝔄` vanish at antipodal (light-transform) polarisations; `M = 0` modes antisymmetric | `cache/t5.log` |

## 7. Revisions after `REFEREE_LIGHT_INPUTS.md`

* **M1:** the mode map was carried out (§3.1, `t5`). (I1) is false for the project's `C⁻₃`. The `C⁻₃` zero is now [D given RB₂]. The bare-commutator result is
  kept for the CLS-D functional, with the referee's `O(γ)` double-trace argument. `ADS4_LIGHT.md` §3 got a marked note. Its `⌈k/2⌉` statement stands as [I]
  for the project's `C⁻₃`.
* **M2:** KKSDZ's general-`d` sum rules are cited and compared (§3.3), and novelty is stated.
* **M3:** the `(ν²+9/4)^{−2}` exclusion is relabelled [I] (it is [D] only given MP form + transverse class), and the `γ` cancellation is relabelled [D-sketch].
* **m1:** the CLS-D sign source is corrected.
* **m2:** the planar `J₀ < 3` requirement is stated.
* **m3:** the Einstein–Maxwell `H` now refers to `⟨JTTJ⟩`.
* **m4:** a frame comparison was added, and "bounded by" replaces "equal".
* **m5:** the tower vs `T`-block route is reconciled (§2.2).
* **m6:** `J₀(G^E) = 2` is kept as [I].

## 8. Follow-up: `Ω′` (see `ads4/OMEGA_PRIME.md`, revised after `REFEREE_OMEGA_PRIME.md`)

The alternative of §0 is worked out in `OMEGA_PRIME.md`: replace `C⁻₃` by `C̃⁻₃`. This is a bare light-ray commutator functional with the CLS-D odd celestial
kernel at `(1,3)` and the lift of `𝔄` at `(2,4)`. It is not shown to be a CLS-D family member; it is valid by the operator statement.
* **Heavy action [DV given the light-ray mode map, `k_r = 3`, and the two-pair heavy formula].** Its leading heavy action equals that of `C⁻₃`, with the exact
  `N(ν) = −(32/7)(2ν+i)(4ν²+9)/(4ν²+201)` (`d₀` convention). There is no flip interference at either pair (structural).
* **Flat limit.** It is inherited, so Theorem 1 applies.
* **Reality.** `C̃⁻₃` is real and even in `ν` [DV at leading-kernel level, `t7`; D-sketch for the full functional].
* **Light side.** It is `∝ t₄²+α_T²` only, with `k(0)` unchanged and no `O(Λ⁻¹)` shift.
* **Open.** The light-ray Jacobian `j(ν)` is not tracked ([I]).
