# Inputs (E) and (R) of the AdS4 light side: superconvergence, Regge counting, and what is left

Status: written 16 Nov 2026 by a derivation agent. **Not refereed.** Scope: inputs (E) and (R) of `ADS4_LIGHT.md` §2.1, §2.3, §4
(see also `REFEREE_ADS4_LIGHT.md` M1, m2). Code: `ads4/code/lightinputs/` (`bash run_all.sh`, < 1 min, single core; logs in
`lightinputs/cache/`, git-ignored). No file of another agent was modified.

Marks: **[D]** derived here (from the items named on the line); **[DV]** derived and checked by the code; **[Q]** quoted
(KKSDZ = Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov arXiv:1904.05905, LaTeX read; CLS-D = arXiv:2311.04271, text read;
CHMRS = arXiv:2106.10274, LaTeX read, App. "omega light" and App. "superconvergence"; MP = Meltzer–Perlmutter 1712.04861, LaTeX
read); **[I]** inferred, not derived; **[NE]** numerical evidence. Nothing here is fitted.

## 0. Results

1. **Exact counting rule [D from Q].** A dispersive functional whose Regge expansion starts at the tensor moment `Π̂_k`
   (`SPINNING_E2E.md` §3.1: kernel `r^{k−2}` in the CHMRS frame) annihilates **exactly** (all `ν`, all frame modes, not
   only at leading order in `ν`) the tree-level light side of any EFT term whose Regge intercept `J_0` (KKSDZ definition) is `< k`.
   The threshold is `J_0 < k` for **both** parities. CLS-D write `k − 1 > J_0` for odd `j′`. By their own crystal formulas
   (108), (111), that is an off-by-one slip (`t3`, §2.3). So a tree exchange of Regge spin `J` is annihilated iff `J < k`, not `J < k − 1`.
2. **Polarisations add no powers of `r` [DV, `t1`].** In the CHMRS frame, the polarisations of all four operators are
   `r`-independent. Every weight-0 tensor building block (`V̂_{i,jk}`, `Ĥ_ij`) is analytic at `r = 0` with a finite limit. Hence
   the Regge growth of every frame mode `𝒢_q` of the stripped tensor is `r^{1−J_0}`, with `J_0` exactly KKSDZ's (fixed generic
   polarisations). The `C^±` normalisation (`γ_a²`, `M`-tables, `|c|²`) is `r`-independent and cannot change this.
3. **(E) for `C⁻₃`: derived, with no Regge assumption on the Einstein correlator [D given (I1) and Q].** The `k = 3` odd
   graviton functional is the *bare* light-ray commutator (CLS-D (111): `f⁻_{3;2} = −1`). This matches flat `B₃`, which has no subtraction. It has
   double zeros on all `[TT]` families, so at `O(1/c_T)` its light side is **exactly** `λ²_TTT·Ψ₃[G_T]`. The `T`-family
   contribution to `⟨T|[E(n₁),E(n₃)]|T⟩` is `(q⁰/Ω)² ε†[H(n₁),H(n₃)]ε` [D + Q KKSDZ]. In `d = 3`, `H(n) = 𝟙 + flip`, with flip `∝ t₄ e^{±4iφ}`
   (and `α_T`, CDP) [DV, `t2`]. Therefore **`[H₁,H₃] = −2i|flip|² sin 4φ₁₃ σ_z`**. The Einstein term is **exactly zero**, the linear `t₄`, `α_T` term is
   **exactly zero**, and the light side of `C⁻₃` is purely `∝ t₄² + α_T² ∝ |ĝ₃|²`. This agrees with `ADS4_LIGHT.md` §2.2 and §2.4.
   (I1) [I]: `C⁻₃` is the leading Regge moment of *some* bare-commutator functional. The zero holds pointwise in `x₂, x₄, z₁, z₃`,
   so it holds for every kernel.
4. **(E) for `C⁺₄`, `C⁻₅` (forward terms) [D given RB₃ and the CHMRS App. E hypotheses].** These are subtracted functionals. Their
   Einstein light side is the `r^{−3}` (resp. `r^{−4}`, …) Laurent coefficient of `𝒢^{E}_+`, or equivalently `ω[G^E]`. It
   vanishes exactly if **RB₃: `J_0(G^E) ≤ 3`**, i.e. every frame mode of the tree-level Einstein `⟨TTTT⟩` is `O(r^{−2})` on the `+` sheet,
   with an integer Laurent expansion. RB₃ is much weaker than the expected `J_0 = 2`. RB₃ is **[Q]**: KKSDZ §"convergence in perturbation
   theory" state `J_0^planar = 2` from tree-level graviton exchange. It is **[I]** here: we computed no AdS4 `⟨TTTT⟩` Regge limit.
   **No component is expected to be nonzero.** If RB₃ held only in the weaker form `J_0 ≤ 4`, `C⁺₄` could pick up an
   Einstein term and `C⁻₅` could not.
5. **(R) [reduced, not derived].** At `k = J` the `γ_{2Δ+k−1}(ν)²` that our `C^±` divide out are *exactly* the Laplace factors
   `c(a,½±iν)` of each plane wave. They cancel identically, by a residue-in-`k` argument (§4.1) [D]. So `R ≡ 1` reduces to
   **Lemma IP**: the `r^{−1}` coefficient of `𝒢^E_+` is a superposition of the *same* transverse-helicity plane waves that a
   heavy state produces (`SPINNING_4PT.md` §4.1), with helicity-independent density `∝ Π^{H₂}_{Δ₀=2}`. The helicity independence
   is [Q] (Einstein AdS shocks: `H = 𝟙`, KKSDZ (adsshock); CEMZ). The transverse plane-wave form is [I]. Given IP:
   `−C⁺₂|_E = 16πG/(ν²+9/4)` exactly and `C⁻₂|_E = 0`. **[D]** what a failure of IP can look like: only a rational factor whose
   singularities sit at `ν² = −(a−5/2)², −(a−9/2)²` or at the zeros of the `M`-table numerators. None of these is at `ν² = −9/4`
   (`t4`). So the referee's `(ν²+9/4)^{−2}` admixture is impossible, because the exchange is a single `H₂` propagator.
   **No direct computation of the `T`-exchange `C⁺₂` was done.**
6. **`⟨JJJJ⟩` warm-up.** The bare current commutator has `J₁+J₃−1 = 1`, so there is no unsubtracted shortcut. Photon exchange (`J = 1`) is
   annihilated by every `k ≥ 2` and graviton exchange by every `k ≥ 3`. Both statements are [D given RB for the Einstein–Maxwell tree]. KKSDZ's Einstein–Maxwell shock
   matrix `H_ij = δ_ij + a₂(n_in_j − δ_ij/2)` gives the analogous statement at `k = 3`: the `J`-block commutator is `∝ a₂²`.

Net effect on `ADS4_LIGHT.md` §4: the row "(E) [I]" splits. **`C⁻₃`: [D given (I1)]**, and this was the item with leverage `~10⁴`. The forward
terms `C⁺₄`, `C⁻₅` are **[D given RB₃]**, where RB₃ is a quoted, unproved Regge bound with a large margin. **(R)** is **[D given IP]**.
Part (P) of the old (E) row is **not** addressed: (P) is the claim that the `|ĝ₃|²` and contact light sides are polynomial in `ν`.

## 1. What is needed

Ω uses `C⁺₂` (all `ν`), `C⁻₃` (all `ν`, with weight `ψ₃(0) = 162`), `∂_{ν²}C⁺₄|₀` and `C⁻₅|₀` (`ADS4_LIGHT.md` §0). The light side at
tree level is `ω[G^EFT]`, where `G^EFT = G^E + t₄G₁ + (t₄² + α_T²)G₂ + contacts`.

* **(E)**: the Einstein piece `G^E` gives exactly 0 in `C⁻₃`, `C⁺₄` and `C⁻₅`.
* **(R)**: the Einstein piece gives exactly `16πG/(ν²+9/4)` in `C⁺₂`, with no leftover `Γ`-ratio.

## 2. The general mechanism

**2.1 Frame lemma F [DV, `t1`].** In `(x₁,x₂,x₃,x₄) = (−y, x/x², y, −x/x²)` with `|x| = 1`:
* `x₁, x₃ = ∓r ŷ` and `x₂, x₄ = ∓x̂`.
* The polarisations are `z(θ) = ŷ + cos θ e_∥ + sin θ e_⊥` at points 1 and 3, and `I(x̂)z′(θ)` at points 2 and 4.
* All of these depend on `(x̂, ŷ, θ)` only, not on `r`.

`t1` evaluates all 12 normalised `V̂_{i,jk} = V_{i,jk}[(P_i·P_j)(P_i·P_k)/(P_j·P_k)]^{−1/2}` and all 6 `Ĥ_ij = H_ij/(P_i·P_j)` at
`r = 10⁻²…10⁻⁵`, `η = 1.3, 2, 5`, at random `θ`. Results:
* Every block has a finite nonzero limit (fitted exponent `|p| ≤ 0.013`, consistent with 0 given `O(r)` corrections).
* `B(r) − B(0) ∝ r¹` or `r²`. `Ĥ₁₃` and `Ĥ₂₄` are exactly `r`-independent.

The stripped tensor is `s₁₃^Δ s₂₄^Δ ⟨…⟩ = u^{−Δ} Σ_I τ_I f_I(u,v)`, with `τ_I` polynomials in these blocks. So every frame mode `𝒢_q(r,η)` grows
like the scalar-like functions `u^{−Δ}f_I`, i.e. like `r^{1−J_0}`. The Fourier projection onto `q` and the `C^±` normalisation are
`r`-independent. This answers the brief's (ii): polarisation contractions add no powers of `r` (of "`s`") in this frame. The flat-space growth of
helicity polarisations under boosts has no analogue here, because the Regge limit is a dilatation, which does not act on polarisations.

**2.2 Spinning carry-over of CHMRS (eftformula) [D given CHMRS App. E hypotheses].** CHMRS App. "omega light" uses the following:
* (a) `ω` is a `dDisc_s` integral, antisymmetric under `s ↔ t`. CLS-D build the spinning functionals antisymmetric under `1 ↔ 3`, (103)–(105).
* (b) `dDisc` commutes with the kernel.
* (c) `G^EFT_+` is single-valued near `r = 0` with integer Laurent powers, `dDisc_s G^EFT → 0` near the `u`-channel Regge limit, and the Euclidean
  terms decay. This holds for tree Witten diagrams.

By Lemma F the contour deformation in `(r,η)` can be done mode by mode at fixed `θ`'s. Then
`−ω|_light = (2πi/4) Σ_{k′≥k} Σ_q ∫[dη] K^{(k′)}_q(η) [𝒢^EFT_{+,q}(r,η)]_{r^{1−k′}}`, the same formula as CHMRS (eftformula) and (moregeneralformula).
Consequence: an EFT term with `𝒢_{+,q} = O(r^{1−n})` for all `q` contributes exactly zero to every functional with leading moment `k > n`.

**2.3 Superconvergence form and the exact threshold [Q + D, `t3`].** KKSDZ (their (J0definition) and the summary of convergence conditions)
define `J_0` by `normalised ⟨O₄O₁O₂O₃⟩ ∈ O(r^{1−J_0})` with fixed generic polarisations; the light-ray commutator `[L[O₁],L[O₃]]` on a common
null plane vanishes if `J₁ + J₃ > J_0 + 1` (and is nonzero iff a Regge pole sits at `J₁+J₃−1`, ref. AnecOPE). For `T`: `J_0 < 3`. Planar
holographic theories have `J_0 = 2` from tree-level graviton exchange (KKSDZ §convergence in perturbation theory) [Q].
CLS-D add subtraction factors `f^±_{k;J_e}`. From their crystal formulas, `f⁺ ~ f⁻ ~ m^{−2(k−2J_e+1)}`, so both raise the kernel's
Regge spin by `k − 2J_e + 1`. The convergence condition is then `J_0 < (2J_e − 1) + (k − 2J_e + 1) = k` for both parities. For `J_e = 2`:
* `k = 2`: `f⁺ = v′−u′` (unsubtracted, marginal for Einstein, which is the graviton pole).
* `k = 3`: `f⁻ = −1`, the bare commutator. It converges for `J_0 < 3`, which agrees with KKSDZ and contradicts CLS-D's text `k − 1 > J_0`.
* `k ≥ 4`: subtracted, `(u′v′)^{−1}`.

## 3. `C⁻₃`: the Einstein term vanishes exactly

* **Light side = one block [D from Q].** For an unsubtracted functional all double-twist families sit at double zeros. CLS-D: "the
  action on all double-traces will remain zero". So `[TT]_{n,ℓ}` enter at `O(γ²) = O(1/c_T²)`, and the `O(1/c_T)` light side is `λ²_TTT Ψ₃[G_T]`.
  The `T` block alone grows like `r^{−1}`, so `Ψ₃[G_T]` converges (`2 < 3`). This needs no statement about `G^E`. Contrast
  `k = 2`: there the single-block action diverges logarithmically and the light side must be taken from the arc, which gives the graviton pole.
* **The block action [D + Q].** By conformal invariance put the common point of the two light transforms at infinity, so that
  `L[T](∞,(1,n)) = ½E(n)`. Fourier transform in `x₂, x₄`.
  * `E(n)` commutes with momentum.
  * The `T` family at fixed timelike `q` consists of the states `ε·T(q)` (in `d = 3`, rest frame: two helicities).
  * Hence the `T`-family part of `⟨T(q)|E(n₁)E(n₃)|T(q)⟩` is `(q⁰/Ω)² Σ_{ε′}⟨ε|H(n₁)|ε′⟩⟨ε′|H(n₃)|ε⟩`.
  * The two orderings decompose in the `(32)` and `(12)` channels, which are the two copies of the `s`-channel `T` block in the antisymmetrised functional.
  * So the `T`-block contribution to the commutator is `ε†[H₁,H₃]ε` pointwise in `x₂, x₄, z₁, z₃`, and every `k = 3` bare-commutator functional (any kernel, any `ν`, any `j′`) gets exactly this.
* **`d = 3` [DV, `t2`].** KKSDZ's shock matrix `H_{ij,kl}(n)` (`t₂`, `t₄`), restricted to the transverse plane, in the helicity basis:
  * `H = [[1, (t₄/4)w^{−4}], [(t₄/4)w⁴, 1]]`, with `w = e^{iφ}`.
  * The `t₂` structure drops out identically in `d = 3`, consistent with `⟨TTT⟩` having two parity-even structures.
  * The diagonal is exactly 1 for all `t₄`. This independently confirms `ADS4_LIGHT.md` §2.3 "diagonal = 1" and the CDP off-diagonal `|t₄|/4`.
  * With `α_T` (CDP), the flip becomes `(−t₄ + iα_T)/4`.
  * Result: `[H₁,H₃] = −2i|flip|² sin(4φ₁₃) σ_z`, i.e. Einstein `0`, `∂_{t₄}|₀ = 0`, and `|flip|² = (t₄²+α_T²)/16`.
  * The helicity-odd `σ_z` is what an odd (`D₅⊗D₅`) moment can see.
* **Remaining step (I1) [I].** Our `C⁻₃` (`D₅⊗D₅` pres–pres, `k = 3`) must be the leading Regge moment of a bare-commutator functional. Support:
  * flat `B₃` is the unsubtracted odd sum rule (`FLAT_GRAVITON_BOUND.md` (2.26b) at `k = 3`);
  * CLS-D (111) identify the unique unsubtracted odd graviton functional with it;
  * `SPINNING_E2E.md` §3.3 maps the `j′ = 1` CLS-D functionals to odd `q₁+q₃` modes.

  The exact mode combination of our `C⁻₃` inside that family was not checked. Because the zero is pointwise, it suffices that `C⁻₃` lies *in* the family.
* What is **not** derived: the `ν`-dependence of the surviving `|ĝ₃|²` term (`(P)`: polynomial `ν⁴ + f₁ν² + f₀`). It is now a
  well-defined single-block action `λ²Ψ₃[G_T]|_{|flip|²}`, computable in principle.

## 4. `C⁺₄`, `C⁻₅`, and (R)

**4.0 Forward terms [D given RB₃].** By §2.2, the Einstein light side of `C⁺₄` (respectively `C⁻₅`) is the `r^{−3}` (`r^{−4}`) Laurent coefficient
of `𝒢^E_{+,q}`, plus higher `k′`. It is zero if `J_0(G^E) ≤ 3` (`≤ 4` for `C⁻₅`). Equivalently, by §2.3, `ω[G^E] = 0`, provided `G^E` obeys boundary
microcausality (as tree Witten diagrams do [I]) and `J_0(G^E) < k`. The required bound has a full unit of margin over the expected
`J_0 = 2`. Status:
* [Q] KKSDZ, and conformal Regge theory for the leading Regge pole of spinning exchanges (MP (aregdiff): `σ^{1−j(ν)}` times `r`-independent
  differential structures).
* Supporting [D]: in flat space every GR helicity amplitude has Regge spin 2 (`f = 8πG/(stu)`, prefactor `~ s⁴`).
* Not checked here: the explicit AdS4 Einstein `⟨TTTT⟩`.

**4.1 (R): why the `Γ`-factors cancel [D].** In the heavy derivation (`SPINNING_E2E.md` §4.4–4.5), each plane wave `e^{2P·x}` with a transverse
helicity tensor transforms through the master formula into `c(a,½±iν)|2P|^{−a}·M(ν;a)`, with `c = 2^{a−1}√π γ_a(ν)`. The
light Regge term in impact-parameter form (MP: `∫dp dp̄ e^{−2ip·x−2ip̄·x̄} ℬ/((−p²)^{h−Δ}…)`) uses the same master formula at `q = −2ip`. The
`ν`-dependence per plane wave is again `γ_a M`; continuation only adds `k`-dependent phases. For a Regge term `r^{1−J}` the `|p|`-integral is a pure power, so
`Π̂_k` has a pole at `k = J`. Its residue is the eftformula value (Mellin in `r`), with `γ` evaluated at `a = 2Δ+J−1`. That is exactly the factor that `C_k` divides
out at `k = J`. So **`γ`'s cancel identically** and `C^±_J|_light` equals the `H₂` transform of the impact-parameter density per helicity class
(scalar check: CHMRS `8πG/(ν²+d²/4)`, `Δφ`-independent, [DV] in `ADS4_LIGHT.md` §1).

**4.2 Lemma IP and the Einstein density.** IP [I]: the `r^{−1}` coefficient of `𝒢^E_+` is a superposition of transverse-helicity plane waves
(no `z·p̂` components), helicity-preserving at both pairs. Plausible from conservation (MP use transversality `z·p̂ → 0`) but not shown.
Einstein density [Q]: the AdS shock is an exact solution with profile `∝ z^d/(z²+y²)^{d−1}`, i.e. the `H_{d−1}` propagator of dimension
`d−1 = 2`, and the probe transfer matrix is `H = 𝟙` (KKSDZ (adsshock), (shockequationgravity); CEMZ). So the density is
`8πG·Π^{H₂}_{Δ₀=2}` per helicity class, helicity-independent. Then `−C⁺₂|_E = 2·8πG/(ν²+9/4)`, with the constant fixed by (F), and `C⁻₂|_E = 0` [D given IP].

**4.3 Possible form of an (R) violation [D, `t4`].** If IP fails, the mismatch enters `C⁺₂` as a ratio of rational `M`-table
entries. In pres–pres even, the `g₁S` parts cancel, and the normalisation `R_LR_R ∝ N_a²/(𝖽_{5/2}𝖽_{9/2})²` with `a = 7` contributes the possible singularities:
`ν² = −81/4`, `−25/4`, and the roots of `N_7` at `ν² = −413/4 ± √9049` (≈ `−8.1`, `−198.4`). None of these is at `−9/4`; `N_7(−9/4) = 18432`. The exchange contributes one simple
pole. Hence `r(ν) ∝ (ν²+9/4)^{−2}` cannot occur. A violation would be `16πG·Q(ν)/(ν²+9/4)` with `Q` rational, `Q → 1` (F). Its
effect on `B_AdS`, `C₁` is given by the formulas of `ADS4_LIGHT.md` §4 (m2).

## 5. Updated input table (replaces the (E), (R) rows of `ADS4_LIGHT.md` §4 if accepted by a referee)

| input | status |
|---|---|
| counting rule: functional with leading moment `k` annihilates EFT terms with `J_0 < k` (both parities) | [D] given CHMRS App. E hypotheses + Lemma F [DV]; CLS-D off-by-one noted [D] |
| (E) `C⁻₃`: no Einstein term, no linear `t₄`/`α_T` term; light side `= λ²Ψ₃[G_T] ∝ t₄²+α_T²` | **[D]** given (I1) [I] and KKSDZ/CDP one-point matrix [Q]; `d = 3` algebra [DV] |
| (E) `C⁺₄`, `C⁻₅`: no Einstein term | **[D] given RB₃** (`J_0(G^E) ≤ 3`): [Q] KKSDZ, [I] here |
| (R) `R ≡ 1` in `C⁺₂` | `γ` cancellation [D]; **[D] given IP [I]** and Einstein eikonal `= 𝟙·Π_{Δ₀=2}` [Q]; violation form restricted [D, `t4`]; `(ν²+9/4)^{−2}` excluded [D] |
| (P) polynomiality of `|ĝ₃|²`/contact light sides | **not addressed** (still [I]) |

**Not done** (would make everything [D]):
* the explicit Regge limit of AdS4 Einstein `⟨TTTT⟩` (momentum space: Raju 1201.6452; spinor-helicity AdS4 results);
* the action of the `k = 2` functional on the `T` exchange in our `C⁺` normalisation;
* the mode map of `C⁻₃` into the CLS-D `k = 3`, `j′ = 1` family.

## 6. Code (`ads4/code/lightinputs/`, `bash run_all.sh`)

| file | content | log |
|---|---|---|
| `t1_frame_regularity.py` | weight-0 tensor blocks in the CHMRS frame as `r → 0`: finite, analytic (Lemma F) | `cache/t1.log` |
| `t2_shock_commutator.py` | KKSDZ shock matrix in `d = 3`, helicity basis; `[H₁,H₃]`; `t₂` drops out; diag = 1 | `cache/t2.log` |
| `t3_clsd_counting.py` | Regge counting of CLS-D subtraction factors; threshold `J_0 < k` for both parities | `cache/t3.log` |
| `t4_mtable_poles.py` | singular points of the `l = 2` `M`-tables (`k = 2…5`) vs the `T` pole | `cache/t4.log` |
