# Unified-theorem scoping: logical map, candidate statements, next step

*21 Nov 2026, scoping agent. Analysis and planning only — no manuscript. Unrefereed, not committed. No existing document or code
was edited. Light computations only (two scripts in the session scratchpad, quoted in §6). Grades as in `RESULT_D_STATEMENT.md`:
**[certified] [P] [P given X] [D] [D-sketch] [Q] [DV] [NE] [M] [I]**; a chain carries the grade of its weakest link.*

**Foundational result (accepted, user decision): Result D, conditional** (`RESULT_D_STATEMENT.md` post-review box). Under A1 (twist
gap `Δ_g` for every single trace ≠ `T`), A2 (tree level, `c_T ≫ Δ_g²`), the analytic inputs V, Adm, H7, L0, H17, L3, per-state positivity
`P` of the fixed functional `Ω^EE` implies `√(t₄²+α_T²) ≤ 2.89·10⁻⁸` (`Δ_g = 1000`), `3.13·10⁻⁶` (`Δ_g = 300`) [D; constants DV],
i.e. the CHLPS form `|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21`. `P` is [NE + M] on the computed near-gap states and is a named hypothesis
(`H_M, H_OO, H_R, H_LJ`) elsewhere.

---

## 0. Answer in brief

1. **No single unified theorem exists today.** The repository holds one foundational conditional result (D) and five independent results.
   Only two of them enter D's proof at all, and neither enters the near-gap positivity that limits D (§1). Nothing in the repository follows from D except
   restatements (§2). The honest description is a **programme**: D at the centre, with independent supporting theorems.
2. **One logical link was not recorded before and is now checked [DV].** The "AdS-only" large-impact-parameter moment `M₁(β)` in D's
   heavy action (`NF02_LARGE_B.md` §4.1, certified positive separately by `h13_cert_m1.py`) is **identical** to Lemma Q's `Q₁(β)`
   (`POSITIVITY_LEMMA.md`). The next moment `M₂ = Q₂` holds in closed form. So Lemma Q gives an analytic proof that both leading terms of
   D's large-`b̂` hybrid action are positive. That region is part of the support for `H_R`, not the binding near-gap region.
3. **ABJM gives no test of D's `t₄` bound.** Every 3d N ≥ 1 SCFT has `t₄ = α_T = 0` exactly ([Q] Buchbinder–Kuzenko–Samsonov
   arXiv:1503.04961; Buchbinder–Stone arXiv:2102.04827; [D] identification below). So D's inequality holds trivially there.
   * D's proof needs a weaker spectral assumption than A1 as written. **A1′:** no single trace other than `T` in `T × T` at `O(c_T^{−1/2})` below twist `Δ_g`.
   * In the type IIA regime ABJM plausibly satisfies A1′ with `Δ_g ~ λ^{1/4}` [D-sketch + Q-physics]. So D *applies* to a known theory but *predicts nothing new* there.
   * The ABJM free-energy theorems are orthogonal to D. `F/c_T → π²/64` tests the 4-derivative/normalisation sector, which `⟨TTT⟩` does not see. `G_N` cancels in D.
4. **The AdS3 HRT theorem has no logical link to D.** Different dimension, opposite direction of implication, no shared input. The link is thematic only.
5. **Recommended next step (new gravitational physics):** extend D from the cubic to the **quartic graviton couplings** (`R⁴`-type
   contacts in `⟨TTTT⟩`). ABJM in its IIA regime supplies the first **non-trivial** test: its `R⁴` coefficient is fixed by localization
   (Binder–Chester–Pufu arXiv:1906.07195). The plan starts with a cheap flat-space go/no-go (§4). The main risk is exactly the forward-term issue
   that made D forward-free (`NOFORWARD_FUNCTIONAL.md` item 1; H9a).

---

## 1. Logical map

Classes: **(a)** used in the proof of D (step, grade); **(b)** follows from D (or D + X); **(c)** logically independent, same physical story;
**(d)** unrelated.

| result | class | precise relation | grade of the relation |
|---|---|---|---|
| **Lemma Q** (`Q₁ > 0`, `Q₂ < 0`, d = 3; `POSITIVITY_LEMMA.md`) | **(a), large-`b̂` step of `P` only** | D's hybrid heavy action (flat NF02 + leading AdS Regge curvature) has the large-`b̂` expansion `F_hyb ≃ A₁M₁(β)/Λ² + A₂M₂(β)/(2Λ³)`, with `A₁ ≥ 355`, `A₂ = −2047` (`ROUTEA_LARGE_B.md`, `NF02_LARGE_B.md` §4.1). **`M₁ ≡ Q₁` and `M₂ ≡ Q₂`.** So Lemma Q gives the sign of both terms for all `β > 0`. The repository certified `M₁ > 0` by a separate route (Arb + large-`β` asymptotics) without noting the identity. **Where it enters:** the certified windows `F_hyb > 0` on `b̂ ∈ [162.6, 597]` (300) and `[274.1, 1997]` (1000), and the [M] model support of `H_R`. **Not:** the near-gap states (`m ≲ 1.3Λ`), which bind `P` and carry `H_LJ`. | identity: `M₂ = Q₂` [P] (`cosh β − 1 = 2 sinh²(β/2)`); `M₁ = Q₁` [DV to 12 digits at β = 0.3, 1, 2.5, 6, three representations] + [D-sketch] (both are the Abel/analytic regularisation of `∫ν P dν`; the `tanh` part integrates to 0 for β > 0 by Mehler–Fock); use in D: [P given the derived envelope bound] for the windows, [M] elsewhere |
| **Lemma M** (monotone improvement ratio; same file) | (c) | scalar CHMRS functional only; not used by `Ω^EE` | [D] (by inspection of D's hypothesis table) |
| **Uniform positivity** (`UNIFORM_POSITIVITY.md`: `∫f(p)P_{−½+iMp}dp > 0`, `f = (1−p)³p(65+155p+47p²)`) | **(c) for D as accepted**; was (a) for the superseded B′-uplift chain (`RESULT_D.md` §3a) | `f` is the CHLPS `F_{g3}` `B₂` weight. `Ω^EE` uses the NF02 weights `ψ_k = (1−p)⁵pX_k(p)/k₀`, so the theorem does not apply as stated. **Transferable method:** the same three-region proof (Mehler kernel, small-`M` monotonicity, expansion at small `Mβ`) applied to D's leading-Regge action `∫(ψ₂,ψ₃)P dp` would upgrade the leading-Regge part of `H_R` ("positive on the whole band, min 0.2263", [M]) to [certified]. | [D] (weights compared); transfer [I] |
| **Flat graviton bound / NF02 certificate** (`FLAT_GRAVITON_BOUND.md`, `NOFORWARD_FUNCTIONAL.md`) | (a) | flat limit of `Ω^EE` = stripped NF02. Its certificate supplies the bulk-point positivity and `A = 3022.34` | [certified] (one float corner) + [P] (`EE ≡ q²` at the crystal) |
| **Dictionary** `ĝ₃ ↔ (t₄, α_T)` (`DICTIONARY_T4.md`) | **(a), final conversion step** | `√(t₄²+α_T²) = 240|ĝ₃|` (AdS units). The convention-free combination is all D uses; the sign/normalisation of `α_T` alone is not needed | [D] + [DV] (three published cases; reproduced by the foundation referee) |
| **Lightcone theorem** (`LIGHTCONE_THEOREM.md`) | **(c)**: neither light-side nor double-trace input | D's light side is the identity + the single `T` block (Lemma R residue, H7 [D]). `[TT]` double traces drop out through the **double zeros** of CLS-D's unsubtracted kernels (L0 [D given Q]). No anomalous-dimension sign is used: at tree level the factor `2sin²(π(τ−10)/2)` vanishes; beyond tree level it is `≥ 0` whatever the sign of `γ`. The lightcone theorem concerns identical external **scalars**, averaged `[φφ]` twists, any `c_T`, no gap. It cannot follow from D (different correlator, no large-`c_T` input). Shared theme: both are imprints of `T` exchange (attraction at large spin vs. Einstein vertex at the bulk point). | [D] (by inspection of L0 and H7) |
| **ABJM** `N^{3/2}` law, finite-`N` bounds, `16N ≤ c_T ≤ 16N²`, `F/c_T → π²/64` (`ABJM_THEOREM.md`) | **(c)** + a worked instance of D's setting; **no (b) beyond triviality** | detailed in §1.1 | see §1.1 |
| **d = 3 energy-correlator sum rule (★), `(S_k)`, kernels** (`EINSTEIN_PROGRAM.md` §3, `KERNELS_D3.md`) | **(c)**, shares one input; plus a **(b) corollary** and a comparison | §1.2 | §1.2 |
| **AdS3 HRT theorem** (repo root) | **(d) logically; (c) thematically** | §1.3 | [D] |

### 1.1 ABJM in detail

* **Does SUSY force `t₄ = α_T = 0`? Yes, already at N = 1.**
  * In a 3d N = 1 SCFT the supercurrent three-point function has a **single** functional form [Q, BKS arXiv:1503.04961 abstract]. Supersymmetry
    forbids parity-odd structures in supercurrent/flavour-current three-point functions [Q, Buchbinder–Stone arXiv:2102.04827, search
    excerpt; full text not read].
  * That single form must include the free N = 1 scalar multiplet (one real scalar + one Majorana). There `c = (n_s − n_f)/(n_s + n_f) = 0`
    (`EINSTEIN_PROGRAM.md` §3 checks; `DICTIONARY_T4.md` §0 `t₄ = 4(n_B − n_F)/(n_B + n_F)`). Hence **`t₄ = α_T = 0` in every 3d N ≥ 1 SCFT**
    [D given Q]. ABJM (N = 6, or 8) has this at every `N, k`, `λ`. The repository already states it (`EINSTEIN_PROGRAM.md` §3, citing 1503.04961).
* **Is ABJM inside D's hypotheses?** Literally no: A1 as written excludes it (light KK towers; `RESULT_D_STATEMENT.md` §1). But D's
  proof only sees operators exchanged in `⟨TTTT⟩`:
  * the heavy side is the set of blocks in `T × T`;
  * the light side is the identity, `T`, and `[TT]` (double zeros);
  * a light single trace `O ∉ T × T` enters only through `[OO]`, with `λ²_{TT[OO]} = O(c_T⁻²)`, i.e. one loop.

  So the proof uses **A1′: no single trace other than `T` appears in `T × T` at `O(c_T^{−1/2})` with twist `< Δ_g`** [D-sketch; to be
  confirmed by a referee that no step of the chain, e.g. a near-gap scan or the `H_LJ` region, uses A1 for operators outside `T × T`].
* **ABJM against A1′:**
  * In 11d/IIA supergravity on `AdS₄ × S⁷/Z_k` (resp. `CP³`), pure Einstein + `Λ` is a consistent truncation (minimal N = 2 gauged
    supergravity is a consistent truncation for any Sasaki–Einstein₇, Gauntlett–Varela arXiv:0707.2315 [Q]; setting the graviphoton to zero is
    consistent [D]). So no supergravity KK mode is linearly sourced by two gravitons. For N = 8 this agrees with representation theory:
    only the `[0,0,2,0]` multiplet contains SO(8) singlets [D-sketch].
  * Hence at tree level the single traces in `T × T` are `T` plus **string states**. In the IIA regime (`λ = N/k ≫ 1`) these start at
    `Δ_g ~ λ^{1/4}` [Q-physics; coefficient not needed].
  * In the M-theory regime (fixed `k`) there are no tree-level heavy states: the UV is 11d Planckian. A2's "heavy states at `O(1/c_T)`"
    fails, so **the M-theory regime is outside D's derivation** [I].
* **Does D + ABJM give a non-trivial prediction?** **No** for `t₄, α_T`: D gives `√(t₄²+α_T²) ≲ C√(ln Δ_g)/Δ_g⁴` and ABJM has exactly 0.
  What survives:
  * **(i)** A2's numerical precondition holds rigorously given localisation and the `N ≥ 6` Ward identity: `c_T ≥ 16N`, so
    `c_T/Δ_g² ≳ 16√(Nk) ≫ 1` with `Δ_g² ~ λ^{1/2}` [P given those inputs + Q-physics for `Δ_g`].
  * **(ii)** With `t₄ = 0` D's sum rule becomes an **equality budget** for heavy data: `Σ_heavy = 16L_A` exactly (§2, corollary C3). In
    the flat limit this is an identity for the Virasoro–Shapiro amplitude (no test). In curved `AdS₄ × CP³` it would test the string spectrum,
    which is not known.
* **Does the finite-`N` free-energy bound combine with D?** **No.**
  * D's bound is `G_N`-independent (the graviton pole and `|ĝ₃|²` term are both `∝ G`, which cancels; `RESULT_D_STATEMENT.md` §3). `G_N`
    enters only through A2.
  * `F/c_T → π²/64` checks that `L²/G_N` read from `F` and from `c_T` agree. Its deviation `F − (π²/64)c_T = 32π²c₁` (Bobev–Charles–Hristov–Reys,
    arXiv:2006.09390, [Q via `ABJM_THEOREM.md`]) measures the **4-derivative** sector. In AdS₄ that sector is Gauss–Bonnet (topological) plus
    Ricci terms (field-redefinable), and it does not affect `t₄`, which first arises from `W³` (`DICTIONARY_T4.md` §0, §3.4 [D]).
  * The two are **orthogonal Einstein diagnostics**: ABJM's checks the 2-point/partition-function normalisation in one theory; D checks
    the 3-vertex in a class. A Gauss–Bonnet coupling is invisible to every local correlator. So **no dispersive argument of D's type
    can bound `F/c_T − π²/64`** [D] — the combination cannot be built from D's machinery.

### 1.2 Energy-correlator sum rules (★)/(S_k)

* **Not used by D.** D uses CLS-D light-ray-commutator functionals on `⟨TTTT⟩`, not the energy–energy correlator in `|T, ±2⟩`. The one
  **shared input** is the conformal-collider one-point structure `⟨T|E(n)|T⟩ ∝ 𝟙 + flips ∝ (t₄, α_T)`: (E5) in Proposition S, Q8 in D's
  `Ĉ₃|_E = 0` (`H7_I14_LIGHTSIDE.md` §4) [Q].
* **Comparison: D is strictly stronger than anything (S_k) + energy-correlator positivity can give.**
  * Without a twist gap on scalars and spin 2, `KERNELS_D3.md` §7.3(ii) shows the system admits `|c|² ≃ 2700/Δ_gap⁴` [P given §5.1].
  * Under D's own A1 (twist gap for all `J`) I reran the LP (`lp_A1.py`, §6). Its feasible points are exact solutions of the spin-truncated
    system: `4|c| = √(t₄²+α_T²) = 1.12, 0.40, 0.146` at `Δ_g = 10, 20, 40` (`J_max = 40`; at `Δ_g = 20`, `J_max = 80` changes it by 1.4 %) [NE].
  * These are ≥ 10⁻¹-sized. D's bound, wherever established, is ≤ 3·10⁻⁶. So **the entire gain of D comes from the spin-2 sum rule with
    the graviton pole** (route S-A), as `KERNELS_D3.md` §7.3 anticipated. The (S_k)-consistent "spectra" with large `|c|` are excluded in
    tree-level A1 theories only through D [D given P].
* **(b) corollary:** with D, (★) gives `Σ_{J≥6}|â^O_{6,2}|² − Σ_{O≠T}|â^O_{−2,2}|² = |c|² ≤ 5.2·10⁻¹⁷` (`Δ_g = 1000`) [D given P + (E)]. This is
  a restatement.

### 1.3 AdS3 HRT theorem

* **No logical link.** The HRT theorem *assumes* a vacuum Einstein (locally AdS₃) bulk and derives entropy inequalities: Einstein ⇒
  entropy cone. D *derives* near-Einstein couplings from CFT data: CFT ⇒ Einstein. Composing them would need either:
  * an AdS₄ covariant entropy-cone theorem (does not exist); or
  * a d = 2 analogue of D. In d = 2 that analogue is **vacuous**: `⟨TTT⟩` is fixed by Virasoro (`c` only), and 3d gravity has no
    propagating gravitons and no `t₄`-type coupling [P, standard].
* **Thematic link only:** both say what Einstein gravity looks like from the boundary.
* **Side observation, AdS3 only, not part of any unified statement.** Theorem 1 of `MANUSCRIPT.md` already allows arbitrary positive weights
  `w_L, w_R`. The only parity-odd local modification of 3d gravity (gravitational Chern–Simons ↔ `c_L ≠ c_R`) would therefore be covered,
  **if** the covariant entropy prescription of topologically massive gravity on these quotients reduces to `min_c (c_L ℓ_L + c_R ℓ_R)/(…)`.
  That reduction is not checked [I].

---

## 2. What follows from Result D (class (b)), honestly

All of these carry D's grade: [D given P] plus the analytic inputs. None is new physics beyond D. They are its physical readings.

| # | corollary | content | value |
|---|---|---|---|
| C1 | **Conformal-collider isotropy** | energy flux in a graviton-created state, `⟨E(θ)⟩ = (E/2π)[1 + |c|cos(4θ−φ)]`, has `|c| = √(t₄²+α_T²)/4 ≤ 7.2·10⁻⁹` (`Δ_g = 1000`), `≤ 7.8·10⁻⁷` (300) | restatement; the cleanest observable form |
| C2 | **Effective boson–fermion balance** | in the free-field basis `t₄ = 4(n_B − n_F)/c_T`, so `|n_B − n_F|_eff/c_T ≤ 7.2·10⁻⁹` (1000). Large-gap CFT3s are "supersymmetric in `⟨TTT⟩`" to that precision; N ≥ 1 SCFTs (ABJM) sit exactly at 0 | restatement; consistent with §1.1 |
| C3 | **Heavy–graviton coupling budget** | every heavy contribution is ≥ 0 and the total is `16L_A − k(0)(1+δ)|ĝ₃|²Λ⁸ ≤ 16L_A` (§6.1 of the statement). So **each** heavy state `X` near the gap obeys `λ²_{TTX}·F_X ≤ 16L_A` (in the `G`-normalisation of the light side), with `F_X/N_s ≥ 0.20`: no massive state couples to two gravitons more strongly than `O(G log Δ_g)` in gap units. Conversion to a bulk coupling `g_{hhX}` not done [I]. In N ≥ 1 theories the budget is saturated (`Σ_heavy = 16L_A`) | mildly new reading (spinning analogue of CHMRS-type heavy-coupling bounds); not checked against any data |
| C4 | **Comparison with (S_k)** | §1.2 | clarifies why route S-B was abandoned |

---

## 3. Candidate unified statements

### U0 (the honest default): **there is no unified theorem; there is a programme**

* **Statement.** Result D is the foundational, conditional gravitational result. The AdS3 HRT theorem, the ABJM theorems, the lightcone theorem,
  Lemma Q/uniform positivity and the energy-correlator sum rules are independent, separately verifiable results. Of these, Lemma Q (large-`b̂`
  signs), the NF02 certificate and the dictionary enter D; the rest do not.
* **Grade.** Each component keeps its own grade. No combined statement is claimed.
* **Recommendation.** Keep this framing in any outward description until U2 or a theorem-grade D exists.

### U1: "Gap-suppressed Einstein universality of the AdS₄ graviton three-vertex"

* **Statement.** Let `𝒯` be a unitary CFT₃ with unique `T`, large `c_T`, satisfying A1′ (§1.1) and A2. Under the analytic inputs and
  per-state positivity `P` of `Ω^EE`, the two non-Einstein `⟨TTT⟩` parameters obey `√(t₄²+α_T²) ≤ 240√(16L_A(Λ)/k(0))/Λ⁴` (2.89·10⁻⁸ at
  `Δ_g = 1000`, 3.13·10⁻⁶ at 300). Equivalently C1/C2 hold. Every N ≥ 1 SCFT satisfies the conclusion exactly; ABJM in its IIA regime is
  an instance of the hypotheses.
* **Combines.** D; dictionary; collider one-point structure; SUSY literature; ABJM as instance.
* **New vs repackaging.** Repackaging of D, except two items: the A1 → A1′ widening (which makes D apply to the known string-regime
  holographic CFT₃s), and the observation that SUSY theories saturate trivially.
* **Missing for a theorem.** Everything D lacks: `H_LJ`, `H_R`, `H_M`, `H_OO`, the KKSDZ loophole, `c₁`. Also a referee check of A1′.
* **Grade.** [D given P] (D); A1′ [D-sketch]; ABJM instance [Q-physics + D-sketch].

### U2 (the real target on the AdS₄ side): "Einstein gravity through quartic order from CFT₃ data"

* **Statement (target).** Under A1′, A2 and positivity of explicit dispersive functionals, the low-energy graviton S-matrix of the bulk
  dual equals Einstein's through quartic order, up to explicit gap-suppressed corrections:
  * (i) the cubic vertex obeys D's bound;
  * (ii) the leading quartic contacts (the D = 4 `R⁴`-type helicity structures `g₄`, `g₅`, … of CHLPS's `f(s,u)` expansion) obey two-sided
    bounds `|g_i|Λ^{n_i} ≤ C_i log(2Λ/3) + D_i` (in units of `8πG`).
* **Test.** ABJM in the IIA regime must satisfy (ii) with its localization-determined `R⁴` coefficient (Binder–Chester–Pufu
  arXiv:1906.07195 [Q]). This would be the **first non-trivial check of the programme's dispersive bounds against an exactly known CFT₃ coefficient**.
  It is non-trivial because `R⁴ ≠ 0` there, whereas `t₄ = 0`.
* **Combines.** D (cubic); a new quartic result; ABJM (localization data + the rigorous `c_T(N,k)` behind A2); dictionary.
* **New physics.** (ii) and the ABJM test are new. No explicit-constant AdS₄/CFT₃ quartic graviton bound from spinning sum rules is known
  to the repository's literature notes; the novelty check is still to be done [I]. The N = 8 numerical bootstrap (Agmon–Chester–Pufu,
  "M-theory archipelago") bounds the `R⁴` coefficient against `c_T` in the M-theory regime. That is a different statement: no gap, SUSY input.
* **Missing.** All of (ii); see §4 for the obstruction.
* **Grade today.** (i) [D given P]; (ii) [I]; test [Q] data, comparison not done.

### U3 (rejected): "Einstein emergence at the 2- and 3-point level"

The idea is to combine D with `F/c_T → π²/64`. Rejected because the 4-derivative sector is topological/field-redefinable in AdS₄ and invisible to
local dispersive bounds (§1.1). The two results have disjoint hypotheses (a class vs. one matrix integral). Any combined sentence would be a
juxtaposition, not a theorem [D].

---

## 4. Highest-value next research step: the quartic graviton couplings ("Result Q"), with the ABJM test

**Why this step.**
* It is the next rung of the same physical question ("is the bulk Einstein?") and is new gravitational physics.
* It reuses D's machinery (CLS-D spinning functionals, CHMRS η-inversion, Lemma R light side, NF02-type certificates).
* It turns ABJM from a vacuous instance into a **test**.
* The foundation-hardening item (`H_LJ` analytic bound, `RESULT_D_STATEMENT.md` §6.2) stays the top priority *for D itself* and should run in
  parallel (owned by the numerics agent).

**The known obstruction (must be faced first).**
* In the stripped basis `B₂⁽¹⁾, B₃⁽¹⁾` see **no** contacts.
* Smeared `B₄…B_K` cannot isolate contacts: each `B_k` sees infinitely many at `O(1)` (`NOFORWARD_FUNCTIONAL.md` item 1 [D]).
* Only **forward terms** (`p = 0`: `B₄|₀`, `∂_{p²}B₄|₀`, `B₅|₀`, …) isolate finitely many contacts.
* In AdS these are the `ν = 0` forward terms whose double-trace/contact content was uncontrolled (H9a, leverage up to 10¹³). D avoided them by design.
* For a **contact bound** the contact content is the signal. The double-trace content of a subtracted CLS-D functional at tree level is
  computable in principle: it is the dDisc of the crossed-channel `T` exchange, the same object Lemma R evaluates for D's pole
  (`H0_SPINNING_FUNCTIONAL.md` comparison table: "finitely many subtracted double-trace families, (P) analysable").

**Plan (go/no-go gates; ≤ 1 core except Stage 3).**

| stage | task | output | cost | gate |
|---|---|---|---|---|
| 0 | **Flat D = 4 feasibility.** Use the AdS-liftable sum-rule set {smeared `B₂⁽¹⁾, B₃⁽¹⁾`} ∪ {the lowest forward terms of `B₄, B₅`}. With the existing NF02/CHLPS code and a linear program on a `(m, J)` grid, find two-sided bounds on `g₄` (and `g₅`) in units of `8πG/M⁶`·log. Check where type II Virasoro–Shapiro (`g₄ ∝ ζ(3)`) lies; compare with CHLPS's published region [Q, to read] | flat two-sided `g₄` bound with the minimal forward content; list of forward terms actually needed | 2–4 days | **No-go** if a two-sided bound needs forward terms beyond `k = 5` or derivatives `∂^{≥2}_{p²}` (AdS lift then unrealistic) |
| 1 | **Which CLS-D functionals lift those forward terms.** The d = 3 specialisation of CLS-D's *subtracted* `k = 4, 5` graviton sum rules (`j′`, `f`). Their crystal images; flip-freeness (the EE analysis of `EE_FOUNDATIONS.md`) | explicit AdS objects `Ĉ₄^{fwd}, Ĉ₅^{fwd}` with flat limit = Stage 0's forward terms | 1–2 weeks | no-go if no flip-free dual structure exists |
| 2 | **Light side exactly at tree level.** (a) `T`-exchange dDisc into the subtracted functionals: Lemma R with subtraction, plus the `[TT]` single-zero families. (b) Contact Witten diagrams of the 8-derivative graviton invariants → their images (polynomial in Mellin). (c) Cross-check (a)+(b) against the known tree-level Einstein `⟨TTTT⟩` in AdS₄ momentum space (Raju arXiv:1201.6452 [Q, to verify]) | `light = pole + Σ_i g_i K_i(ν) + double-trace(ν)` with every term computed | 2–4 weeks | the leverage of residual errors on `g₄` must be ≤ 10² (vs 10¹³ for D's old forward terms); else stop and report |
| 3 | **Heavy positivity** of the combined functional (NF02-part + forward part) near the gap: the same chains as `Ω^EE`; reuse `H_M`-type models | `n = 1` or not at `Δ_g = 300, 1000` | numerics agent; weeks | — |
| 4 | **ABJM test.** Express BCP's IIA `R⁴` coefficient (1906.07195) in the `g₄` normalisation (via `⟨SSSS⟩ → ⟨TTTT⟩` SUSY Ward identity at the flat-space/contact level) and compare with the Stage 3 bound at `Δ_g ≈` lowest IIA string state (`∝ λ^{1/4}`, coefficient from the literature) | pass/fail of a non-trivial consistency test | 1 week | — |

**Fallback if Stage 0 or 1 fails.** Record the negative result: "forward-free spinning sum rules cannot bound quartic contacts in AdS₄". This
is a genuine statement about which Einstein-emergence statements are reachable by this method. Then return to U1 + theorem-grading D (`H_LJ`).

**Cheap side tasks (hours each; foundation, not new physics):**
1. Replace `h13_cert_m1` by Lemma Q: write a short proof that `M₁ = Q₁`, i.e. that the two regularisations agree (§1 table); cite Lemma Q
   in `NF02_LARGE_B.md`/`ROUTEA_LARGE_B.md` via whoever edits those files.
2. Apply `UNIFORM_POSITIVITY.md`'s method to the NF02 leading-Regge action, to certify the leading-Regge part of `H_R`.
3. Referee check of A1 → A1′ (§1.1); if confirmed, change "A1 excludes ABJM" in `RESULT_D_STATEMENT.md` §1.

---

## 5. Grades of the new claims in this note

| claim | grade |
|---|---|
| `M₂ ≡ Q₂` | [P] |
| `M₁ ≡ Q₁` | [DV] (12 digits, 4 values of β, 3 representations) + [D-sketch] |
| Lemma Q gives both signs of D's large-`b̂` envelope | [P given the envelope expansion ([D] in `NF02_LARGE_B.md`)] |
| uniform positivity not applicable to `Ω^EE` as stated | [D] |
| lightcone theorem not used in D | [D] |
| `t₄ = α_T = 0` in all 3d N ≥ 1 SCFTs | [D given Q] |
| D needs only A1′ | [D-sketch], review needed |
| ABJM (IIA) satisfies A1′ with `Δ_g ~ λ^{1/4}`; M-theory regime outside A2 | [D-sketch + Q-physics]; [I] |
| no dispersive bound on `F/c_T − π²/64` | [D] |
| (S_k) under A1 admits `√(t₄²+α_T²) = 0.40` at `Δ_g = 20` | [NE] (LP feasible point, grid positivity) |
| C1–C3 | [D given P] |
| AdS3 weighted-entropy extension | [I] |

---

## 6. Hand-over

* **Files created:** this note. Scratchpad (not in repo; may be cleared):
  `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/unified/m1q1.py` (`M₁` vs `Q₁`, three
  representations, ~1 min) and `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/unified/lp_A1.py`
  (variant of `ads4/code/kernels/lp_balance.py` with A1: every operator at twist `≥ Δ_g`; usage `python3 lp_A1.py gap Jmax`, seconds).
  Outputs:
  * `M₁ = Q₁ = 0.0412709827709, 0.037558523835, 0.0231910297119, 0.00369640605981` at β = 0.3, 1, 2.5, 6;
  * LP max `|c|²` = 0.0779 (`Δ_g` = 10, `J_max` 40), 0.00978/0.00991 (20; `J_max` 40/80), 0.00132 (40, 40); unconverged 6.3·10⁻⁴ (80) and
    2.5·10⁻⁴ (160) at `J_max = 60 < Δ_g` (not used).
* **Literature consulted (abstract level):**
  * Buchbinder–Kuzenko–Samsonov, [arXiv:1503.04961](https://arxiv.org/abs/1503.04961) (single supercurrent structure);
  * Buchbinder–Stone, [arXiv:2102.04827](https://arxiv.org/pdf/2102.04827) and [arXiv:2108.01865](https://arxiv.org/abs/2108.01865) (no
    parity-odd supercurrent structures; tension only for superspin-2 currents), via search excerpts;
  * Binder–Chester–Pufu, [arXiv:1906.07195](https://arxiv.org/abs/1906.07195) (ABJM `R⁴` from localization, IIA regime);
  * CHLPS [arXiv:2201.06602](https://arxiv.org/abs/2201.06602) (two-sided Wilson-coefficient bounds, D = 4).

  Quoted from memory and **to be verified before use**: Gauntlett–Varela arXiv:0707.2315 (consistent truncation); Raju arXiv:1201.6452
  (AdS₄ `⟨TTTT⟩`); Agmon–Chester–Pufu (N = 8 numerical bootstrap of `R⁴`).
* **Open questions for the next agent:** (1) Stage 0 of §4 — the decisive cheap test. (2) A1′ review. (3) The three cheap side tasks.
  (4) Novelty search for an AdS₄ quartic graviton bound from spinning sum rules.
* **Not done:** no full-text reading of the SUSY papers; no check of CHLPS's quartic results; no check that the TMG entropy functional
  reduces to weighted chiral lengths.
