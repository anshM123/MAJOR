# Referee report: the foundation of Claim D (`RESULT_D.md`) and its new, unrefereed inputs

Referee: hostile foundation review, 18 Nov 2026. Scope: `RESULT_D.md`, `STATUS.md`, and the unrefereed inputs
`EXACT_SKELETON_NEARGAP.md`, `SPINNING_TTTT_EXT.md`, `PARITY_ODD_HEAVY.md`, `DIRECT_GRAVITON_NEARGAP.md`, `SPINNING_MELLIN.md`,
`OBJECT_IDENTITY.md`, `FUNCTIONAL_REDEFINITION.md`, `DICTIONARY_T4.md`. `CERTIFIED_NEARGAP.md` was not reviewed.
Referee scripts: `ads4/code/referee_foundation/` (`rf1_kernels.py`, `rf2_bound.py`, `rf3_skeleton.py`, `rf4_misc.py`; outputs in its
`cache/`, git-ignored). Nothing else was modified.

## 0. Verdict

**The numbers are right; the claim built on them is not the claim the numbers support.**

* Every load-bearing number I re-derived came out as stated:
  * the exact scalar-skeleton moments, reproduced by an independent implementation (my own `J = 0` Mack data, Barnes lemma,
    quadrature and `n`-sum) to `≤ 5·10⁻⁸` (`Δ_g = 300`), `≤ 5·10⁻⁷` (`Δ_g = 1000`) and `≤ 7·10⁻⁶` (`Δ_g = 100`, worst node `x = 0.99`);
  * the corrected circle weights `(1/36, 1/16, 1)`, from the `M`-coefficients recomputed from source;
  * the flip-freeness of `K*` and `𝔖_ff`;
  * the `ν`-dependence of `N₄^{K*}` and `N₆^{K*}`;
  * the bounds `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` (D-safe) and `4.24·10⁻⁵ / 5.28·10⁻⁷ / 4.13·10⁻⁹` (D-strong), with `A`, `B_AdS` and
    the factor 240;
  * the bulk-point density ratios of `OBJECT_IDENTITY.md`.
* **But the headline, "D-safe is positive down to the gap, `n = 1`", is not established for the functional that Claim D defines,
  nor for gravitons.**
  * The heavy action is computed in a model: CHMRS's scalar `C_{k,ν}` on scalar skeletons, with the graviton spin put in through
    flat B′ densities. The scalar remainder this model computes exactly is `≤ 0.02 N_s`. The remaining `≈ 98 %` of the near-gap action
    is the flat B′ dictionary at `x = ν/m → 1`. For gravitons that dictionary has been tested only at `x ≲ 0.1`, on nodes that carry
    `0.4 %` of the norm.
  * By the authors' own `OBJECT_IDENTITY.md`, the objects the functional actually uses (`Ψ_k/N_k`, `k ≥ 3`) have different bulk-point
    densities. So the flat limit of `Ω^{NF,SC}` is **not** NF02. The flat certificate, and the constants `A`, `B_AdS`, belong to a
    different functional.
  * The proposed fix, a spinning CHMRS-type `C_k`, is a functional that does not yet exist (H14′). **H14 is not even listed among the
    hypotheses of either version of Claim D.**
* Several grades in `RESULT_D.md` are wrong in the optimistic direction:
  * H9 "absent in D-safe" is false;
  * H13′ outside the scan quotes the wrong functional's envelope;
  * the H12 tolerance quoted is the wrong scope;
  * the "established" table contains conditional items.
* Three inputs the chain uses are missing from the hypothesis list:
  * `|ĝ₃|²`-freedom of `Ψ_{5,6,7}`;
  * the I14 swappability;
  * the existence of the spin-2 dispersive functional `Ĉ₂` that carries the whole light side.
* **Novelty.** There is no prior explicit-constant `Δ_g` bound on `t₄`, `α_T` in CFT₃ (I agree with `DICTIONARY_T4.md` §7). But:
  * the parametric statement and the `log(M R_AdS)` translation are in arXiv:2201.06602 §5;
  * the forward-free `B₂ + B₃` class is that paper's footnote functional;
  * the D-safe constant `A = 3022` is `121×` CLPS's optimum, so the explicit bound is far from what the method could give.

**Recommendation:** do not describe Claim D as "D-safe positive down to the gap". Restate it as conditional on a model, with H14/H14′
and the missing inputs added (§3, §4).

## 1. Findings, ranked

### FATAL (for the headline as worded)

**F1. The D-safe statement is not about a well-defined functional with the properties used.** Three different objects carry three
different parts of the argument.

1. **Near-gap heavy action.** It is computed for CHMRS's scalar `C_{k,ν}` (`OBJECT_IDENTITY.md` §1: "the skeleton model identifies every
   term, near the gap, with CHMRS's dispersive `C_{k,ν}`").
2. **The functional actually defined** (`OMEGA_NF.md` "Definition"): `Ĉ_{3..7} = Ψ_k/N_k`. These are defined only through their leading
   Regge kernel, with `j_k(ν)` untracked. Their near-gap action is "not defined by the project's own definitions" (`OBJECT_IDENTITY.md`
   §0.1).
   * Their bulk-point densities differ from `C`'s at `O(x⁴)`, and at `O(1)` as `x → 1`. Reproduced in `rf4`; e.g. `Ψ/C = 0.317, 0.004`
     at `x = 0.9, 0.99` for `k = 4`.
   * The flat (bulk-point) limit of `Ω^{NF,SC}` is therefore **not NF02**:
     * the NF02 certificate (`NOFORWARD_FUNCTIONAL.md`) does not apply to it;
     * its flat heavy positivity is a float scan only: `J ≤ 80`, min `0.178` for the actual `mix` assignment (`OBJECT_IDENTITY.md` §3.1);
     * its light side (`k(0)`, `A`, `B_AdS`, the contact cancellation) was never computed. `OBJECT_IDENTITY.md` §4 lists it as [I].
3. **`Ĉ₂`, which carries the graviton pole and the entire `|ĝ₃|²` light side.** For gravitons it is "a spin-2 dispersive functional whose
   leading Regge moment is `Π̂₂` … (I10-type assumption)" (`OBJECT_IDENTITY.md` §1 table). Its existence beyond the leading kernel is
   assumed.
4. **The redefinition** (`FUNCTIONAL_REDEFINITION.md`).
   * It removes H14 only for external **scalars**, i.e. for `⟨φφφφ⟩`, which does not contain `t₄` or `α_T`.
   * For `⟨TTTT⟩` it trades H14 for **H14′**: the existence of a spinning `C`-type functional that has not been built. That document says
     building it is research-level.
   * Odd `k` stays `Ψ_k`, whose density (`S_k` factors) is not the B′ odd density. So even after the redefinition the flat limit is not
     NF02.

**Evidence of the omission.**
* `RESULT_D.md` l.30–31: the hypothesis columns of both versions omit H14.
* l.73: "H14 removed at scalar level for k = 2, 4, 6".
* The table row still advertises "n = 1 … at all three `Δ_g`".

**Fix.**
* State Claim D for an explicitly named functional, `Ω^{NF,C}` with spinning `C`-type `Ĉ_k`.
* Add **H14′** (existence, Einstein-freeness, bulk-point density equal to B′, analysable (P)) and **H0** (existence of the spinning
  dispersive `Ĉ₂` with the stated leading kernel) to the hypothesis column of both versions.
* Or keep `Ψ_k` and recompute the flat limit, its certificate and `A`, `B_AdS` for the `Ψ` densities.
* Replace "positive down to the gap" by "positive on the scanned states **in the scalar-skeleton B′ model**".

### MAJOR

**M1. H13′ outside the scanned set quotes the wrong functional; the D-safe scan is much narrower than stated.**

*Where the wrong envelope is quoted.* `RESULT_D.md` l.72, and `EXACT_SKELETON_NEARGAP.md` §8 "Not done": "beyond, OMEGA_U's
`κ ≥ 1.3–2.6` [M]". Those `κ` are for `Ω^U` (M02).

*What the refereed NF02 analysis says* (`OMEGA_NF.md`, I15, after `REFEREE_OMEGA_NF.md` M1):
* `b̂ > 20` is **open, [I]**;
* positivity there needs a pointwise relative AdS correction `≲ 4·10⁻³` (`Δ_g = 100`) or `≲ 9·10⁻⁵` (`Δ_g = 1000`);
* the `C₁ = 1` envelope was withdrawn as uncalibrated.

*The actual exact-skeleton NF02 coverage* (`rf4_misc.py`, from `e7_<Δ_g>_recs.pkl`):

| `Δ_g` | `m/Λ` scanned | `J` | `b̂` max | records with `b̂ > 20`, `m ≥ 1.5Λ` |
|---|---|---|---|---|
| 100 | `[1, 3]` | ≤ 28 (≤ 9 at `m < 1.1Λ`) | 43.8 | 40 |
| 300 | `[1, 2]` (`J ≤ 10`, `b̂ ≤ 14` above `1.5Λ`) | ≤ 20 | 38.3 | **0** |
| 1000 | **`[1, 1.5]` only** | ≤ 20 | 40.2 | **0** |

* "`b̂ ≤ 40–60`" (`RESULT_D.md` l.72; `STATUS.md` l.27) overstates the scan. The maximum is 38–44.
* At `Δ_g = 1000` no NF02 state with `m > 1.5Λ` was evaluated with exact moments.
* The region where `OMEGA_NF`'s truncated series failed (`b̂ ∈ [10, 40]`, `m` up to `~3Λ`) is sampled only at `Δ_g = 100`. The statement
  "the series' negative `κ` was an artefact" is extrapolated from there to `Δ_g = 300, 1000`.

**Fix.** Grade H13′ for D-safe as **[NE] on the listed states, [I] elsewhere** (I15). Give the scan envelope per `Δ_g`. Drop "60".

**M2. H9 is not "absent in D-safe".**
* (P), the double-trace and contact content, is unanalysed for the smeared even replacements `Ψ₄`, `Ψ₆`. Their weight
  `sgn(v′−u′)r^{k−3}` is non-meromorphic, so no `dDisc` rewriting is available (`SUPERCONVERGENT_UPLIFT.md` item 4; `OMEGA_NF.md` I8).
* Its normalised test-residue leverage in `Ω^{NF,SC}` is **48 / 200 / 5 / 16 / 0.4** (`k = 3…7`, `Δ_g = 100`) and **2.6·10³ / 1.1·10⁴ / …**
  at `Δ_g = 1000` (`OMEGA_NF.md` item 3). It **grows with `Δ_g`**.
* `RESULT_D.md` l.68 says "absent in D-safe", and the D-safe column omits H9.

**Fix.** Split H9 into H9a (forward pieces, D-strong only) and H9b ((P) for smeared even `Ψ_k`, both versions, [I], leverage as quoted).
Add H9b to the D-safe column.

**M3. Inputs that the chain uses but `RESULT_D.md` does not list.**

| input | where it is stated | grade there |
|---|---|---|
| (a) `|ĝ₃|²` absent from the light side of `Ψ_{5,6,7}`; needs convergence for `J₀ ≤ 4`, which is shown only in the Regge region | `SUPERCONVERGENT_UPLIFT.md` item 4 (withdrawn to [I]); `OMEGA_NF.md` I10 | [I] |
| (b) I14: convergence/swappability of `C⁺₂` and the higher moments; it is used in the Einstein-absence argument that the "established" table calls [D] | `OMEGA_U.md` input list | [I] |
| (c) P1(`j′=0`) and the sign property for even `Ψ_k`; `j_k(ν)` untracked and possibly `k`-dependent | `OMEGA_NF.md` I17 | — |
| (d) H0, the existence of the spinning `Ĉ₂` dispersive functional | F1 | — |
| (e) `k_AdS − k(0)`, the `O(Λ⁻²)` correction to the `|ĝ₃|²` coefficient | `omegau/cache/nf5.log` | 4.2 % at `n = 1`, `Δ_g = 100` |

Input (a) bears directly on the coefficient `k` that sets the bound. Input (e) is small, but it belongs in the error budget.

**Fix.** List all five, with grades.

**M4. The D-safe binding state is a parity-odd exchange; H12′ is load-bearing, and the quoted tolerance is the wrong scope.**
* **Odd `J` gives a parity-odd structure.** For identical `T` the single odd-`J` `⟨TTO⟩` structure is parity-odd: `KERNELS_D3.md` §3.1
  table, "odd `J ≥ 5`: 0 even / 1 odd"; `DIRECT_GRAVITON_NEARGAP.md` §0. So every odd-`J` `+−` heavy state is a parity-odd exchange.
* **The global minimum of the D-safe action** at all three `Δ_g` is `+−, J = 5`: `+0.2267` (100), `+0.2057` (300), `+0.2022` (1000)
  (`rf4`). Over parity-even exchanges only, the minimum is `+−, J = 4` (`+0.2412 / 0.2203 / 0.2103`).
* **No test covers these states.** Every `⟨TTTT⟩` end-to-end test is at even `J` (20, 36, 52, 14, 4, 6) on parity-even structures, and
  so is reading B (H2).
* **The narrative and the tolerance are wrong.**
  * `PARITY_ODD_HEAVY.md` §2.3 and §3.3 describe parity-odd exchanges as "ordinary same-helicity contributions". That holds only at even
    `J`.
  * `RESULT_D.md` l.71 quotes the `++`-scope tolerance `|c₁| ≲ 8` (D-strong). The relevant "all" scope is `|c₁| ≤ 2.2 / 2.1` at
    `Δ_g = 100 / 300` (`δ = 0.03`, `PARITY_ODD_HEAVY.md` §2.5). That is **below** the mutation scale `4–10`. The D-safe tolerance
    (`≥ 19`) is unaffected.

**Fix.**
* State that H12′ covers all odd-`J` `+−` states, including the D-safe binding state.
* Quote the "all"-scope tolerance.
* Make the `⟨TTO⟩` odd-`J` end-to-end point the first item of §6.

**M5. What "exact near-gap moments" buys is mislabelled.**
* **The model.** The action is `F_exact = F_hyb + S_exact`. On every scanned D-safe state `|S_exact| ≤ 0.02 N_s`, and
  `F_exact/N_s ≈ F_hyb/N_s` to `≤ 0.006` (`e7_*.txt`). So the D-safe `n = 1` is, to `98 %`, the flat NF02 margin (`0.200`) transported
  by the B′ dictionary (I1 / H2) with the curvature term.
* **The scalar part is genuine and correct.** I reproduced it independently (§2a). But it is not where the uncertainty is.
* **Where the uncertainty is.** The load-bearing items are the graviton densities at `x ~ 1`: H2, H4 and H14. These are tested only at
  `x ≲ 0.1`, on `≤ 0.4 %` of the remainder norm (`DIRECT_GRAVITON_NEARGAP.md` §7.4). `+−` at `k ≤ 6` is noise-limited.
* **The mislabel.** The `RESULT_D.md` table header "`n` (exact skeleton moments, scanned states)" suggests that the threshold is computed
  from exact moments of the relevant objects.

**Fix.** Label `n` as "scalar-skeleton B′ model [M + NE]". State that the exact part is `≤ 0.02 N_s`.

**M6. H1 is overgraded in the "established" table, and the kernel switch is a tautology.**
* **What is tested.** The `⟨TTTT⟩` [DVE] covers `β = 1`, `ν ≤ 3`, three masses, even `J`, parity-even structures.
* **What is not.**
  * `k = 2` rests on an [I] cleaning of defective caches.
  * `¼` is untested end to end.
  * `𝔖_ff` is the kernel `Ω^{NF,SC}` is defined with (`OMEGA_NF.md`; `DIRECT_GRAVITON_NEARGAP.md` §1 table). It is **not resolved** at
    low `J` (`β ≈ 0.09`) for `k ≤ 4` even at same helicity (`DIRECT_GRAVITON_NEARGAP.md` §2, §4). It fails at `β = 0.4`, `k = 4`, opposite
    helicity (`0.196`), and at `k = 2`, opposite helicity everywhere.
* **The kernel switch.** `RESULT_D.md` §6 says "the even kernel is now `K*`", but no document redefines `Ω^{NF,SC}` with `K*`. The
  sentence "switching the even kernel to `K*` reproduces every threshold and bound exactly" (l.60) is true by construction: `e10` applies
  a response ratio that is asserted to be 1 on every node (`EXACT_SKELETON_NEARGAP.md` §9). It is no evidence.
* **The finite-`m` sensitivity** (`c₂ ≤ 268`) uses `c₂` measured at `ν ≤ 3` and applies it at `x ~ 1`. `c₂(ν)` is non-monotone
  (`−56, −12, −305`).

**Fix.**
* Move ⟨TTTT⟩ H1 out of the "established" table, or restrict the row to "`β = 1`, `ν ≤ 3`, even `J`, parity-even, `k ≥ 4`".
* Redefine `Ω^{NF,SC}` with `K*` explicitly.
* Mark the kernel-switch statement as model-tautological.

### MINOR

* **m1. Twist gap, not dimension gap.** `RESULT_D.md` §1 says "a gap `Δ_g` to all single-trace operators of spin `J > 2`". All scans
  impose `m ≥ Λ`, i.e. `τ = Δ − J ≥ Δ_g`, which is CHMRS's convention (their footnote: "`Δ_gap` is conceptually a twist gap"). With a
  dimension gap, states with `m < Λ` are allowed: `m/Λ = 0.978` at `J = 20` and `0.86` at `J = 50` (`Δ = 100`, `rf4`). Say "twist gap".
* **m2. Error term.** The `A log + B_AdS` form has an `O(n/Λ)` error, not `O(n²/Λ²)`, because `ψ₂″(0) ≠ 0`. Exact `L_AdS` gives `7509`
  against `7320` at `Δ_g = 100` (`rf2`). It is harmless, because the tables use the exact integral, but fix the formula.
* **m3. Conditional rows in the "established" table.** It includes the spinning graviton pole [D given F, R] and "Einstein term absent"
  [D given `J₀ ≤ 2`, microcausality, I14]. Move both to the hypotheses, or mark them conditional.
* **m4. The Arb error radius is not an error bound.** "Arb error radius ≤ 1e−96" (`EXACT_SKELETON_NEARGAP.md` §7) does not bound the
  `n`-tail quadrature (Gauss nodes from doubles), the `M`-truncation or the odd-`k` interpolation. My independent reproduction agrees to
  `10⁻⁹–7·10⁻⁶`, so nothing is wrong, but do not call it a certificate.
* **m5. Flip responses do not all vanish.** `PARITY_ODD_HEAVY.md` §2.3 says the flip responses "vanish identically" for all kernels of
  the functionals. `Ψ₄` and `Ψ₆` carry a `j′ = 0` celestial structure at `(1,3)` with a **nonzero** flip response
  (`SUPERCONVERGENT_UPLIFT.md` item 3; residual `≤ 1.7·10⁻⁶ N_s`). The `λ_Eλ_O` cancellation still follows from parity invariance, as
  §2.2 argues, so the conclusion stands. Correct the sentence.
* **m6. `STATUS.md` is stale.** Link 11 still says "H1 … NOT for opposite helicity under the even kernels (open, load-bearing)". It
  contradicts `RESULT_D.md`. Harmonise.
* **m7. Scope of H11.** "No other light single traces" excludes every known holographic CFT₃ (ABJM has light KK modes). CLPS treat
  light matter as entering with a helpful sign; NF02 was not checked for that. State the scope.
* **m8. Novelty and attribution.**
  * The scaling and the AdS₄ translation are in CLPS §5: "translates to a bound on three-point OPE coefficients of stress tensors in
    holographic 3d CFTs. In AdS₄, infrared logarithms are simply replaced by `log(MR_AdS)`".
  * The `B₂ + B₃` forward-free class is CLPS's footnote functional (`1598 log − 3211`).
  * CLS-D (2311.04271, §5) list Regge positivity, the non-shock sum rules and the light-side dictionary as the open steps. The present
    chain assumes rather than closes them (H6, H13′, F1).
  * The `−240` dictionary coefficient appears new as a printed formula, and is checked (ECG, QTG).
  * State the D-safe constant against CLPS's optimum: `√(3022/24.9) ≈ 11×` weaker.
* **m9. D-strong sensitivity to the object.** For D-strong the object choice moves `n` within `1.00–1.10`, i.e. the bound by up to
  `×1.4` (`rf2`: the bound doubles from `n = 1` to `1.2`). The quoted single value `4.2·10⁻⁵` hides this. The same applies to the
  `H12′`-"all" scope (`n = 1.15` at `c₁ = 16`).

## 2. Independent computations (scripts in `ads4/code/referee_foundation/`)

**(a) Exact scalar-skeleton moments (`rf3_skeleton.py`).**

*What I implemented independently of `mack.py`, `skelexact*.py` and `exactmom.py`:*
* the `J = 0` Mack data from the Dolan–Osborn expansion: `Q_n ∝ (b)_n²/(n!(Δ−½)_n)`, `b = 1−Δφ+Δ/2`;
* Barnes' first lemma, re-derived analytically (`|Γ(c)|⁴ ₂F₁(c, c̄; A; 1−e^{−2u})/Γ(A)`), with mpmath `₂F₁` on a Chebyshev grid;
* the `y`-integral by Gauss–Laguerre, and the `n`-sum directly to 1500 plus a mapped Euler–Maclaurin remainder.

*Normalisation.* One global constant is fixed at `k = 2`, `ν = 0.03`. Everything else is a test.

| state (`J = 0`, `m = Λ`, `Δφ = 3`) | method behind the cache | `k = 2` max rel. diff | `k = 4, 6` max rel. diff | nodes |
|---|---|---|---|---|
| `Δ_g = 100` | double contour | `7.2·10⁻⁶` (`x = 0.992`); `≤ 1.8·10⁻⁹` for `x ≤ 0.79` | `≤ 1.2·10⁻⁷` | 8 nodes, `x` up to `0.9997` |
| `Δ_g = 300` | Arb Laplace | `≤ 2.8·10⁻⁸` | `≤ 4.8·10⁻⁸` | 7 nodes, `x` up to `0.9997` |
| `Δ_g = 1000` | Arb Laplace | `≤ 5·10⁻⁷` | `≤ 2.6·10⁻⁸` | 7 nodes, `x` up to `0.9997` |

* The `10⁻⁶` residuals at `Δ_g = 100`, `x > 0.94` are the double-contour values. They match the authors' own Arb-vs-double check
  (`e2b`: max `|Δρ| = 1.1·10⁻⁵`).
* An earlier version of my code disagreed at `6.7·10⁻⁵` (`k = 6`, `Δ_g = 1000`). The cause was a cancellation in my evaluation of
  `G_κ(t)`. After an `expm1`-based fix the disagreement fell to `2·10⁻¹¹`, so the cache was right.

Verdict: the exact skeleton moments are **correct**. Note that at `x = 0.9997` the exact `m^{2k}C` is `10⁻⁹–10⁻⁸` of the bulk-point
`E_k`, which is why the near-gap remainder is `O(1)` in `E_k` units.

**(b) Even kernels (`rf1_kernels.py`).**
* The `M`-coefficients `l = 2` were recomputed from `laplace.M_coeff`, not from the pickled cache.
* The one-pair flip responses per circle mode are in the ratio `36 : 16 : 1` exactly, for `a = 7, 9, 11` and `ν = 0…10`.
* With `n_|q| = (1/36, 1/16, 1)` the Lemma-L flips are `1 : −1 : 1` to machine precision. The preserved ratios at `ν = 0`, `a = 7` are
  `−0.270947`, `+0.129489`, as quoted.
* `K*/𝔖_ff` preserved equals `(4ν²+(4a−7)²)/(2(4ν²+8(a−2)²+1))` at all test points.
* `K*` and `𝔖_ff` flips are `0` (`≤ 2·10⁻¹⁶`).
* The alternative `n = 1/C(4,2+q)` fails the Lemma-L flip test (`−1.143`, `1.714`). So the `1/C²` weights are the right ones.

**(c) Closed forms of `N^{K*}`.** Symbolically, `N₄^{K*} ∝ (16ν⁴+216ν²+263)/((4ν²+81)(4ν²+841))` and
`N₆^{K*} ∝ (16ν⁴+536ν²+2023)/((4ν²+169)(4ν²+1369))`.
* The `ν`-dependence matches exactly. The overall constant differs by 16 because my kernel is normalised `1/16` of the project's; it
  cancels in `Ψ_k/N_k`.
* Both are positive on ℝ and `g₁`-free.
* I also checked the celestial weights against the Fourier coefficients of `|sin ψ|^{−δ}`.

**(d) Bound formula (`rf2_bound.py`).**
* `A`, `B_AdS` were recovered from the exact `ψ₂`: NF02 `3022.33`, `−6552.5`; M02 `51.737`, `−86.297`. `B_AdS = B_flat + A log(2/3)`,
  i.e. `m_IR = 3/(2L)`.
* All six tabulated bounds are reproduced from the exact `L_AdS` integral.
* `240 = √57600 = 2(d−1)d(d+1)(d+2)|_{d=3}`.
* The ECG cross-check `t₄ = −240·(21/4)μ = −1260μ` and the QTG cross-check `3780μ` are reproduced from the stated actions.
* Sensitivity: `bound(n = 1.2)/bound(1) = 2.0`, `bound(1.5)/bound(1) = 4.7`.

**(e) Parity-odd interference.**
* The structural argument is correct: a parity-invariant functional gives `ω[W_EO] = −ω[W_EO] = 0`. Products `D₅ ⊗ D₅` and `j′`-matched
  `(1,3) ⊗ (2,4)` are parity-even.
* But see m5 (nonzero `Ψ₄`, `Ψ₆` flip) and M4 (odd-`J` `+−` states are parity-odd and bind D-safe).

**(f) Bulk-point densities (`rf4_misc.py`).**
* `Π̂/C = (1−r*²)^k/(1+r*²)` and `Ψ/C = (1−r*²)^{k+2}/(1+r*²)` are identical to the stated density formulas (to `10⁻¹²`).
* Table values: `k = 4`, `Π̂`: `0.997 / 0.940 / 0.443 / 0.023`; `Ψ`: `0.996 / 0.917 / 0.317 / 0.004`.
* The small-`x` slope `(k+1)/16` confirms `O(x⁴)`.
* The `Ψ` kernel itself remains [D-sketch]: the large-`ν` identity was not checked independently.

## 3. Re-grading of the hypotheses

| # | `RESULT_D.md` grade | referee grade | reason |
|---|---|---|---|
| H1 | [DVE] (`β = 1`, `k = 4–7`, both helicities, three kernels); `k = 2` [NE] | **[NE]** for ⟨TTTT⟩ as used. [DVE] only for `β = 1`, `ν ≤ 3`, even `J`, parity-even, `k ≥ 4`, under `K*`/`𝔄`. `k = 2` [NE + I] (cleaning). `𝔖_ff` (the defined kernel) unresolved at low `J`, `k ≤ 4` | M6 |
| H2 | [NE] | [NE] at `x ≤ 0.1`, one `b̂`, even `J`; **[I] at `x ~ 1`**, where it is load-bearing | M5 |
| H3 | [DVE] / [M] | [M] (odd `k` by continuation or interpolation); odd-`J` exchanges are H12′ | — |
| H4 | [I] | [I] (agree); the `++` evidence covers `≤ 0.4 %` of the norm | M5 |
| H5 | not needed | agree for the skeleton part | — |
| H6 | [I] | [I]; must also contain H0 (existence of spinning `Ĉ₂`) | F1, M3 |
| H7 | [D-sketch] | [D-sketch] | — |
| H8 | [Q] | [Q] | — |
| H9 | [I], "absent in D-safe" | **[I] in both**: H9a forward (D-strong); H9b (P) for smeared even `Ψ_k` (D-safe, leverage up to `1.1·10⁴` at `Δ_g = 1000`) | M2 |
| H10 | [A] | [A] | — |
| H11 | [A]/[I] | [A]/[I]; scope excludes all known holographic CFT₃s; twist gap | m1, m7 |
| H12 | [D] flat; interference-free; H12′ [I] | [D] flat; [D] interference-free (flip sentence corrected); **H12′ [I] and load-bearing for the D-safe binding state** (`+−`, `J = 5`); D-strong tolerance `|c₁| ≤ 2.2` ("all") | M4, m5 |
| H13′ | [NE] inside, [M] outside | [NE] inside (envelope of §1 M1); **[I] outside for D-safe** (I15); [M] outside for D-strong only | M1 |
| H14 | [I], bounded [M]; "removed at scalar level" | **[I], and missing from both hypothesis columns.** Redefinition ⇒ H14′ [I] (non-existent object). The scalar-level removal is irrelevant to ⟨TTTT⟩ | F1 |
| (new) H0 | — | [I]: existence of the spinning `Ĉ₂` dispersive functional | F1 |
| (new) `|ĝ₃|²`-free `Ψ_{5,6,7}` | — | [I] | M3 |
| (new) I14, P1(`j′=0`), `j_k(ν)` | — | [I] | M3 |

**"Established" rows:**

| row | `RESULT_D.md` | referee |
|---|---|---|
| Flat `F_{g3}` positivity, flat bound | [certified] | agree (refereed before) |
| JHEP `F_{g3}` not positive | [D] | agree |
| M01, M02 | [certified] | agree |
| Unstripped = stripped span | [D] | agree |
| No forward-free smeared `B_{k≥4}`; NF02 | [D]; certified except a corner | agree for NF02 **as a flat functional**. It is not the flat limit of `Ω^{NF,SC}` as defined (F1) |
| Heavy Regge formula, ⟨TTTT⟩ | [DVE] | restrict the scope (M6); not established for odd `J`, parity-odd or low `J` with `𝔖_ff` |
| Graviton pole in `C⁺₂` | scalar [DV]; spinning [D given F, R] | conditional; move to the hypotheses (m3) |
| Einstein term absent from `Ψ_k` | [D] given `J₀ ≤ 2` | conditional on RB₂, microcausality and I14 (m3, M3) |
| Uniform positivity, scalar-kernel `B₂` | [certified] | agree |
| Exact scalar-skeleton moments | [D] method; numbers | **agree, independently reproduced** (§2a). Scalar only |

## 4. Statements in `RESULT_D.md` that are overclaimed

1. l.30, the D-safe row: "`n` = **1** at all three `Δ_g`". This is a model threshold on scanned states (F1, M5). The hypotheses column
   omits H14, H9b, H0 and the M3 items.
2. l.30–31: neither hypothesis column lists H14.
3. l.54: the "exact scalar-skeleton moments" row reads as if near-gap positivity were exact. It is exact for the scalar remainder only
   (M5).
4. l.60: "switching the even kernel to `K*` reproduces every threshold and bound exactly". True by construction (M6). Also "`k = 2` …
   both helicities": this rests on the [I] cleaning.
5. l.50 and l.60: ⟨TTTT⟩ "[DVE] … both helicities" without the restrictions to even `J`, parity-even, `ν ≤ 3` and `β = 1` (M6, M4).
6. l.68: H9 "absent in D-safe" (M2).
7. l.71: H12 "tolerance … `≲ 8` (D-strong)". The relevant scope gives `2.2` (M4). "Flip responses vanish for `K*`, `𝔖_ff`, `𝔄`" is not
   true of `Ψ₄`, `Ψ₆` (m5).
8. l.72: H13′ "high-`J` samples (`J ≤ 28`, `b̂ ≤ 40–60`) … beyond, `OMEGA_U`'s `κ ≥ 1.3–2.6` [M]". The maximum `b̂` is 38–44. For D-safe
   the `OMEGA_U` envelope does not apply; it is [I] (M1).
9. l.73: H14 "**H14 removed at scalar level for `k = 2, 4, 6`**" and "bounds … unchanged". The removal concerns external scalars, which
   do not see `t₄`; for ⟨TTTT⟩ it becomes H14′ (F1).
10. l.16: "gap … to all single-trace operators of spin `J > 2`". It must be a twist gap (m1).
11. l.24: "`+ O(n²/Λ²)`". It is `O(n/Λ)` (m2).
12. l.37–39: "here with explicit constants and an explicit hypothesis list". The list is incomplete (M3, F1).
13. `STATUS.md` l.27: "D-safe … positive down to the gap on all scanned states incl. `b̂ ≤ 40–60`". Same as items 1 and 8.

## 5. What would close the gaps (in order of value)

1. **Fix the object.** Define the D-safe functional explicitly: the spinning `C`-type `Ĉ_k`, `K*` for even `k`.
   * Either construct it (H14′, H0),
   * or keep `Ψ_k` and recompute its flat limit, certify it, and recompute its light side (`k(0)`, `A`, `B_AdS`, contacts).
2. **H12′ at odd `J`.** One end-to-end `+−`, odd-`J` (parity-odd `⟨TTO⟩`) point. It decides the D-safe binding state.
3. **H13′ for NF02.** Exact-skeleton NF02 scans at `Δ_g = 300, 1000` for `m ∈ [1.5, 10]Λ` and `b̂ ∈ [20, 200]` (the Laplace method is
   fast: 4–10 s per state), plus a large-`b̂` asymptotic envelope.
4. **H9b.** (P) for the even `Ψ_k`, or replace them by meromorphic-weight, `η`-inverted functionals (`FUNCTIONAL_REDEFINITION.md` §3),
   for which (P) is analysable.
5. **H2/H4 at `x ~ 1` for gravitons.** This is the real near-gap content. Exact skeleton moments do not touch it.
