# Hostile referee report: Result D, route A′ with `Ĉ₂′` (state of 19 Nov 2026)

*Referee agent, 19 Nov 2026. Read: `RESULT_D.md` (update box and hypothesis table), `H0_SPINNING_FUNCTIONAL.md` §R, §22–§25,
`H7_I14_LIGHTSIDE.md` (all rounds, §0–§24), `PARITY_ODD_ROUTEA.md` (§0–§12), `CERTIFIED_NEARGAP_A.md` (Theorem A),
`ROUTEA_LARGE_B.md`, `NF02_LARGE_B.md` §1, and the code in `code/spinfunc`, `code/h7i14`, `code/exactskel`, `code/omegau`.
No existing document or code was edited. Nothing was committed. The cheap re-runs (one core) are listed in §4; their scripts are in the
session scratchpad, not in the repository. They only import project modules read-only, and each one is described in §4 well enough to
rebuild it.*

Grades used: **[certified]**, **[P]** proved, **[P given X]**, **[D]** derived (not independently checked), **[D-sketch]**, **[Q]**
quoted, **[DV]**, **[NE]** numerical evidence, **[M]** model, **[I]** assumed or open, **[A]** standing assumption.

---

## 0. Summary verdict

* **The implication "heavy positivity ⇒ bound" is sound, and the constants are right.** I recomputed `A = 16ψ₂′(0) = 3022.337`,
  `k(0) = 1`, `16L_A = 7508.5 / 10734.6 / 14339.6` and the bounds `2.209·10⁻⁴ / 3.132·10⁻⁶ / 2.891·10⁻⁸` from the NF02 data (exact
  `p`-integral). The sign chain is correct: heavy ≥ 0 ⇒ `−Ω_light ≥ 0` ⇒ `k|ĝ₃|²Λ⁸ ≤ 16L_A`. I found no fatal error in the logic.
* **What is not established is the premise "`n = 1` at `Δ_g = 300`" over all heavy states.**
  * Graviton-specific heavy actions exist only for even `J ≤ 10` and odd `J = 5, 7, 9` near the gap, in a model (§2, item M).
  * Every other state (near-gap `J ≥ 11`, the whole Regge / large-`b̂` region, the wedge `m < 1.5Λ` at large `J`) is covered only by
    the scalar-skeleton model, with no graviton-specific correction at all.
  * The computed graviton corrections **grow steeply with `J`**. A crude extrapolation reaches the base-model margin at `J ≈ 15–16`
    (`Δ_g = 300`). This is the most serious open point (finding S1).
* **The `Re` prescription is legitimate.** The β-phase strip does not touch any `Δ_g ≥ 300` number. The gauge (null space), Lemma N and
  `δk(0)` move nothing at 300/1000 beyond `10⁻³`.
* **A definitional inconsistency in `Ĉ₃` is visible in the numbers** (finding S4). The chains apply the scalar `k = 2` de-shift to `Ĉ₃`,
  which §R.5 says needs none. At `Δ_g = 300`, `J = 10`, this choice changes the largest diagonal correction from 0.021 to 0.072 `N_s`.
* **The defensible statement at `Δ_g ≥ 300` is a conditional one.** The chain "positivity ⇒ bound" is [D]. Positivity is [NE + M] on
  `J ≤ 10` near the gap and [I] for gravitons elsewhere. The bound holds at the two computed gaps, 300 and 1000; it is not shown "for
  all `Δ_g ≥ 300`" (§2).

---

## 1. Findings, ranked by severity

Severity scale: **fatal** (the stated conclusion is false or does not follow), **serious** (a load-bearing step is unestablished or
mis-graded), **minor** (a correct step is mis-stated, or the effect is small and quantified), **presentational**.

### Fatal

None found. The derivation from per-state positivity to the bound is valid. Every constant I checked is reproduced. No computed number
breaks `n = 1` at `Δ_g = 300` or `1000`.

### Serious

**S1. `n = 1` at `Δ_g = 300` is unestablished for `J ≥ 11` near the gap, and the computed corrections grow fast in `J`.**
* *Location.*
  * `RESULT_D.md` update box: "Combined even + odd … `n = 1` at `Δ_g = 300` and `1000`"; the Ĉ₂′ consolidation bullet (iii).
  * `PARITY_ODD_ROUTEA.md` §9 says the uncovered `J` have "margin ≥ 0.4 in the base model".
  * `H0` §23.4 (rows), `l9_neargap.py` (`JMAX = 10`).
* *What the evidence is.*
  * Graviton corrections exist only for even `J ≤ 10` (2×2) and odd `J = 5, 7, 9` (1×1).
  * The exact-skeleton state set at 300 has `J ≤ 20` (`e1_300_Df3_hiJ`: `J = 11–20` at `m/Λ ∈ [1.04, 1.3]`). `J = 11–20` enter
    uncorrected.
* *Why the "margin ≥ 0.4" argument fails.* I re-ran the `Ĉ₂′` chain (`C23H_K2`; §4 item 5) and read off the per-state corrections at
  `Δ_g = 300`. Units are `N_s`, worst `m`, for `J = 4 / 6 / 8 / 10`:

  | quantity at `Δ_g = 300` | `J = 4` | `J = 6` | `J = 8` | `J = 10` |
  |---|---|---|---|---|
  | off-diagonal `|ΔF₊₋|/√(N₊₊N₊₋)` (definition-independent) | 0.001 | 0.014 | 0.037 | 0.088 |
  | largest diagonal `|ΔF_hh|/N_h`, with the `k = 3` scalar-`g` shift (as implemented) | 0.004 | 0.009 | 0.024 | 0.072 |
  | largest diagonal `|ΔF_hh|/N_h`, without the `k = 3` shift | 0.000 | 0.002 | 0.005 | 0.021 |

  * The off-diagonal grows by ×2.4–2.6 per `ΔJ = 2`.
  * Extrapolated geometrically it reaches ≈ 0.2 at `J = 12`, ≈ 0.5 at `J = 14` and ≈ 1.2 at `J = 16`.
  * The uncorrected diagonals there are 0.7–1.0, so a 2×2 with `|off| ≳ 0.85` loses positivity at `J ≈ 15–16`.
  * This is a crude extrapolation through the regime `J²/m ≳ 0.5`, where H0 §22.3 already calls the Laplace series unusable (`J = 12` at
    `Δ < 400`). It does not show a failure. **It shows that `J ≥ 12` is load-bearing and unknown, not safe by margin.**
  * At `Δ_g = 1000` the corrections are ≈ 10× smaller (`∝ 1/m²`: `J = 10` off-diagonal 0.020), so the same crossing sits at
    `J ≳ 30` [extrapolation].
* *Why it matters.* The bound needs positivity on every state. As things stand, `n = 1` at 300 holds for gravitons only on `J ≤ 10`,
  with `J ≥ 11` resting on an assumption.
* *Fix.*
  1. Compute `Γ_{2,3}` (Ĉ₂′, Ĉ₃) for `J = 12, 14, 16, 20` (even, 2×2) and `11, 13, 15` (odd) on the `Δ = 400, 800, 1600` rows. The
     Laplace chain is usable there for `J²/m ≲ 0.5`. Use the n1/n3 direct integrator where it is not.
  2. Add `J = 21–40` near-gap states to the exact-skeleton set at `Δ_g = 300`.
  3. Until then, state `n = 1` at 300 as "[NE + M] for `J ≤ 10`; [I] for `J ≥ 11`".

**S2. The Regge / large-`b̂` region is shown only for the scalar-skeleton hybrid, never for gravitons. The near-gap wedge
`m ∈ (Λ, 1.5Λ)`, `J > 20` had no scan at all.**
* *Location.* `ROUTEA_LARGE_B.md` §0, §4.3 coverage map; `RESULT_D.md` box (H13′ bullet).
* *What is there.*
  * `F_hyb` (flat NF02 + leading AdS Regge term) and the scalar-skeleton remainder `S` on `m ≥ 1.5Λ`, `b̂ ≥ 20`.
  * The "[certified]" window `b̂ ∈ [91.2, 2Λ]` certifies `F_hyb`, a scalar-skeleton object.
  * For gravitons this is the right leading-Regge object (with `Ĉ₂′`, (HN) is exact, so the graviton leading-Regge action is
    `P_l(ρ⁺⁺ + ρ⁺⁻)`).
  * But at large `b̂` positivity is carried by the `O(Λ⁻²)` AdS term `A₁M₁(β)/Λ²`.
  * Any graviton-specific beyond-leading-Regge correction is compared against that `O(Λ⁻²)` term, not against `O(1)`. None has been
    computed. Its scalar analogue is calibrated at `|S|/F_hyb ≤ 122/m²`.
  * The steep `J`-growth in S1 shows that graviton corrections need not track the scalar ones.
* *The wedge.*
  * The twist gap allows `J` up to `≈ Λ(m²/Λ² − 1)/2`: 187 at `m = 1.5Λ`, `Δ_g = 300`, and 625 at 1000.
  * The near-gap sets stop at `J = 20` (`m ≤ 1.3Λ`); the large-`b̂` scans start at `m = 1.5Λ`.
  * I scanned `F_hyb` over the wedge (§4 item 2). Every state is positive. The minimum is at small `J`: 0.2056 at 300 (`+−`, `J = 5`,
    `1.02Λ`) and 0.2009 at 1000 (`+−`, `J = 5`, `1.005Λ`). States with `J > 20` have `F_hyb/N_s ≥ 0.27` (300) and `≥ 0.21` (1000).
  * So the model is fine there, but this was not in any document, and it is model-only.
* *Fix.*
  1. Compute the `O(r)` (subleading-Regge) term of the spinning kernel. The split-representation machinery of `v3_regge_kernel` with the
     dual structure at `x₀ ≠ 0` is the step H7_I14 §24.8 already proposes for Lemma N. It gives the graviton-specific large-`b̂`
     correction directly.
  2. Add the wedge to `h13a_hyb` (`m/Λ ∈ [1.005, 1.5]`).
  3. Grade the large-`b̂` graviton statement [M] (scalar-skeleton hybrid) + [I] (graviton beyond leading Regge), not [NE].

**S3. The model behind "[NE + M]" is not stated in `RESULT_D`, and its graviton part has no non-perturbative check at the masses that
matter for `Δ_g ≥ 300`.**
* *Location.* Every "[NE + M]" in the box; `l9_neargap.py`, `p6_odd_neargap.py`, `v3_heavy_impact.py`, `v4_odd_impact.py`.
* *What M is* (from the code):
  1. **Base.** Stripped flat NF02 densities, plus the exact leading-Regge curvature `P_l − J₀`, plus the exact CHMRS scalar-skeleton
     remainder `ρ_k` (`Δφ = 3` or `2`, both run). The remainder is transplanted onto `+−` with the stripped envelope `(1−y)^{−2}`. No
     `η`-inversion.
  2. **Graviton correction.** `Γ_{k,h}(x; m)`: the `K₂` Laplace expansion of the graviton/scalar crystal ratio in linear form.
     * The ν-only factor is removed by a same-`J` heavy reference at `Δ_ref = 2·10⁴`.
     * The local first-order de-shift is applied: `u_T` for `k = 2`, the scalar `g` for `k = 3` (see S4), then `Re`.
     * A per-state constant `a_h` is fitted at `x ≤ 0.25`, `ν ≥ 15`. It is a congruence once both `k` are corrected.
  3. **Interpolation.**
     * Cubic spline in `x` on the nodes `x = 0.05 … 0.8`, held beyond.
     * Linear in `1/m` between the rows `Δ = 100, (200), 400, 1600`, held outside.
     * For `J ≤ 8` at `Δ_g = 300`, a state at `m ≈ 1.04Λ` takes ≈ 10 % of its correction from the `Δ = 100` row. That row is
       unconverged at small `x` (H0 §23.4, §25.5: `t = 0.66–1.04` against direct values).
     * For `J = 10` it takes ≈ 29 % from the `Δ = 200` row.
  4. **Coverage.** Even `J ≤ 10`, odd `J = 5, 7, 9`; the parity-odd even-`J` block is uncorrected (S5); per-state positive
     semidefiniteness of EE ⊕ OO.
* *Soundness.*
  * It is a reasonable model. The scalar part is exact.
  * The graviton part was checked against a non-perturbative evaluation only on the `Δ = 100` row (`m ≈ 98`): `K₂` vs direct
    `|ΔΓ| ≤ 117` c₂ on `x ≥ 0.25`; `t = 0.66–1.04` at `x = 0.1, 0.15`.
  * At `Δ_g = 300` (`m ≈ 300–600`) the `K₂` expansion should be better, but no direct point exists on the `Δ = 400` row.
  * Robustness within the model is good:
    * the whole-correction window is `t* = (−10.2, +8.3)` at 300;
    * dropping the `Δ = 100, 200` rows keeps `n = 1` (min 0.2054, `t* = (−5.6, +5.1)`; §4 item 6).
  * `t*` scales the computed correction along one direction only. It does not bound errors of another shape, for example state-localised
    ones or ones growing in `J`.
* *Fix.*
  1. Write M into `RESULT_D` §3 as a named model.
  2. Add a handful of direct (n1) graviton points on the `Δ = 400` row (`J = 8, 10`; `x = 0.1, 0.3, 0.5`) to validate `K₂` at the
     masses that matter for 300.
  3. State that "[NE + M]" means "in model M, `J ≤ 10`".

**S4. `Ĉ₃` is defined three different ways, and the chains use the one that the definition section rules out.**
* *Location.*
  * `H0` §R.5: "`k = 3` needs no inversion (`j = 0`)".
  * `H0` §R definition update: "`Ĉ₃` is unchanged".
  * `H0` §24.1: "`Ĉ₃` likewise (CHMRS-type `C₃` …; its own shift relation is not derived — [I])".
  * Code: `l9_analysis.shift` and `p4_odd_analysis.shift` apply the **`k = 2` scalar `g` (`γ₇`)** de-shift to `k = 3`.
    `v3_heavy_impact` and `v4_*` keep this "on purpose".
* *Why it matters.*
  * `f₃ = −1` carries no `η`, so `Ĉ₃`'s leading kernel is `−r¹P_l` ([Q], §R table). On that reading the de-shift is spurious.
  * The spurious shift term is of the `2αJ(J+1)` type and grows like `J²`.
  * Re-running the `Ĉ₂′` chain without the `k = 3` shift (§4 item 5) gives:
    * `n` unchanged at 300 and 1000;
    * `J = 6 / 8 / 10` minima at 300 move by `+0.004 / −0.006 / +0.002`;
    * **but the largest `J = 10` diagonal correction at 300 drops from 0.072 to 0.021 `N_s`**.
  * So the ambiguity is not load-bearing at `J ≤ 10`. It is a sizeable part of the computed correction at `J = 10`, and it feeds the
    `J`-growth of S1.
* *Fix.*
  1. Fix `Ĉ₃` in §R as `Re[Ψ_{3,ν}[o_H q]/N̂₃(ν)]`, with **no** de-shift. That is justified if its spinning leading-Regge heavy image is
     `σ_h c₃(ν)P_l`. Check this with `v3_regge_kernel`'s split representation and an `o_H` dual; it is about a day's work.
  2. Remove the `k = 3` shift from the chains.
  3. If the spinning image is not `l`-isotypic, define `Ĉ₃′` by its own spinning de-shift, as for `Ĉ₂′`.

### Minor

**M1. The twist-gap convention in `RESULT_D` §1 disagrees with the code.**
* `RESULT_D` §1 says: twist `Δ − J ≥ Δ_g − 2` for `J > 2`, and "no other light single traces" for `J ≤ 2`.
* Every scan (`omegau/common.state`: "None if unphysical (`tau < Dg`)") uses `Δ − J ≥ Δ_g` for **all** `J`.
* The §1 convention admits states with `α ∈ [Δ_g − 4, Δ_g − 2)`. I checked them in the hybrid model (§4 item 4): min 0.2021 (300),
  0.2005 (1000), both at `+−`, `J = 5`, `α = Δ_g − 4`. They are positive, but uncovered by the graviton-corrected chain.
* *Fix.* State the hypothesis as the code uses it: "every exchanged single-trace operator other than `T` has `Δ − J ≥ Δ_g`". This also
  fixes the vague H11 for `J ≤ 2`.

**M2. The coverage map's "`b̂ < 20`: near-gap scans" is wrong above `m = 2Λ` (300) and `1.5Λ` (1000).**
* The exact-skeleton sets stop at `m/Λ = 2.0` (300) and `1.5` (1000).
* I scanned `F_hyb` with `b̂ < 20` at `m/Λ ∈ [1.6, 30]` (§4 item 3): min 0.482 (`++`, `J = 8`, `1.6Λ`), rising to ≥ 0.976 at `10Λ`.
* *Fix.* Add this scan to `ROUTEA_LARGE_B` and correct the map.

**M3. A term linear in `t₄` in the `Ĉ₂` light side is not excluded at `O(1/Λ²)`.**
* *Location.* `H0` §15.2; `ADS4_LIGHT.md` item 3: a ν-independent `16πG c₁ t₄`, with `c₁` not computed.
* It is not addressed in the H7 note ("Cross terms `E × t₄` … are H6/H17, not H7"), and it is not in the H17 budget.
* Size: at the bound it shifts `16L_A` by **`0.43 % × c₁` at 300 and `0.033 % × c₁` at 1000** (§4 item 1).
* It is harmless for `|c₁| ≲ 10`. It is [I].
* *Fix.* Compute `c₁` (the Regge residue of the `E × R³` `T`-exchange, Lemma R with `n = 2`), or bound it. List it explicitly (H6-lin
  below).

**M4. H17 is a heuristic model.**
* The model takes `|a_j| = D × 36` per power lowered. The sign is unknown.
* The budget at 300 is ±0.2 % on the bound. Even a 10× underestimate would give ±2 %, irrelevant to "≲".
* The `Ĉ₂′` `δk(0)` (`−7.66·10⁻⁴` / `−6.83·10⁻⁵`, h1 convention) belongs inside this budget. `RESULT_D` and `H0` §R.9 still quote the
  ηP values (`−1.638·10⁻⁴` / `−1.463·10⁻⁵`).

**M5. The light-side shape at `ν ≲ 8` rests on the §18.1 kernel, which is [D].**
* That kernel is a split-representation reading of `Ψ⁺`, with Lemma V polarisation referral and dropped ν-only constants. It is
  validated at large ν (check X, `a₁ = 22.06i` vs the crystal's `22i`) and in the scalar limit, but not at `ν ≲ 8`.
* `B_AdS` comes from `ν = O(1)`. An `O(1)` error in `h_∓` there would shift `B_AdS` by up to `~10³`, i.e. ≲ 5–10 % on the bound at 300.
  The heavy side would move coherently, with negligible weight at 300.
* Calling H7 "[P given …]" hides that the given includes this [D], plus D2 [D] and R4 [D-sketch]. The honest overall grade is **[D]**.
* *Fix.* One independent small-ν check of `h₋ + h₊` (forward normalisation). For example, the n3 contour on the `Δ_ref = 2·10⁴`
  reference at `ν ≈ 9–15`, where it is validated, compared with the `N_com R` closed form.

**M6. The leading Regge kernel of `Ĉ₃` for spinning externals ("(HN₃)") is quoted from CLS-D, not derived.**
* The large-`b̂` hybrid uses `σu∫ψ₃[P_l − J₀]` (the `σu` part of `A₁` is 50.6 of 377.8).
* `A₁ ≥ 355 > 0` survives any `O(1)` change in that part, so this is not load-bearing for the sign. It is needed for the definition (S4).

**M7. The functional as defined is not the functional as computed; the difference is negligible at `Δ_g ≥ 300`.**
* Defined: `Re∫XΨ̂` with a complex non-local Abel weight `X` in the min-norm gauge.
* Computed: the local first-order de-shift `(1 ± iσ_T)/2`, with `Im X·Im Ψ̂` dropped.
* The whole de-shift moves `λ_min` by `≤ 3·10⁻³` at 300 (H7_I14 §20.3). The gauge moves it by `< 10⁻⁴` (§24.5). The null space lives at
  `ν ≲ 18`, i.e. `x ≲ 0.06`, a negligible weight at 300.
* *Fix.* Say so in §R ("heavy chains use the local form; error ≤ the full de-shift effect").

**M8. Continuity in `m` and the gap edge.**
* The per-`J` minimum sits at the twist edge `α = Δ_g − 2`, which is generally not a grid point.
* Hybrid at the edge for `J = 5`, `+−`, 300: 0.2043 at `m = 1.0166Λ`, against the grid's 0.2054 at `1.02Λ` (§4 item 4).
* Not load-bearing, but the minimum should be quoted at the edge.

**M9. The parity-odd even-`J` block (`X^O_same`, couples like `++`) is never corrected.**
* *Location.* PARITY_ODD_ROUTEA §6 covers odd `J` only; H12′ "Not run: even-J parity-odd".
* At leading order it is the `++` density with `N = tan²(πτ/2) ≥ 0`, so its sign is safe.
* Its `O(1/m²)` graviton correction is uncomputed. It is probably as small as the EE `++` one (`≤ 0.004 N_s` at `J = 4`, 300), but that
  is [M].
* *Fix.* One `p6`-type run at even `J`.

**M10. "`Δ_g ≥ 300`" is an extrapolation from two gaps.**
* `n` was computed at `Δ_g = 300` and `1000` only.
* Monotonic behaviour (corrections `∝ 1/Λ²`, minima 0.2054 → 0.2021) is plausible, not shown.
* *Fix.* State the bound at the computed gaps. Or add one intermediate gap (500) and one larger (3000) in the base + `K₂` chain; the
  analysis is cheap once the rows exist.

### Presentational

* **P1.** The tables of `RESULT_D` §2–§4 describe the withdrawn B′-uplift functional. "D-strong" (`Ω^U`, M02, B′ objects) was never
  re-examined under the stripped dictionary and should be marked withdrawn. Replacements are in §3 below.
* **P2.** `B_AdS = −6552.5` still appears in `RESULT_D` §2 and `H0` §R.9, §15.3. The value is `−6552.66`.
* **P3.** `H0` §R still has the pre-`Ĉ₂′` items:
  * item 3: normalisation by `N̂₂`; `Ĉ₂′` is normalisation-free, `N̂₂ ≡ 1`;
  * item 5: the ηP `ψ̃₂` ODE; it should be the Abel equation with `u_T`;
  * item 9: ηP `δk(0)`;
  * item 10: still lists H12′/H4 in their old form.

  The update box above §R supersedes them, but a self-contained §R should be rewritten.
* **P4.** Several grades are written "[P given X]" where X contains [D]/[D-sketch]/[I] items: H7, I14 (K1′ is not written out; K2 for
  spinning components is [D]; K3 is [D]). The overall grade should be the weakest input: [D].
* **P5.** The `RESULT_D` box is 20 bullets long, and later bullets silently supersede earlier ones. For example, the "δk(0)" values are
  stated three times with two different kernels, and the "n = 1.5" (§22) result is withdrawn two bullets later. A referee cannot read
  off the current state without this report. Collapse the box to the current state.

### Checked and found sound (no finding)

* **(a) Bound formula and constants.**
  * `k(0)|ĝ₃|²Λ⁸ ≤ 16L_A` follows from heavy ≥ 0. It uses `−Ĉ₂′|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` and the flat images
    `−B₂ ∋ 2πG|ĝ₃|²p⁶`, `−B₃ ∋ −2πG|ĝ₃|²p⁴`. `G` cancels.
  * The `n`-placement in §R.9 is consistent.
  * `A`, `B_AdS`, `L_A` and the bounds are reproduced (§4 item 1).
  * The log form differs from the exact `p`-integral by 0.6 % at 300 (it is the `O(1/L)` term). The headline correctly uses the exact
    integral.
* **(c) `Re` prescription.**
  * The real part of a valid complex sum rule with real OPE-coefficient products is a valid sum rule. Conditions (a) and (b) of H0 §23.5
    hold for the EE and OO blocks: Lorentzian-real structures, real `β`, EO cancelled by the `S⁰` sum.
  * Condition (c) is automatic for `Ĉ₂′`: its heavy leading action is exactly `P_l`, and the light side `Â_h[D]` is real (V3).
  * The prescription is load-bearing, since `absmix` fails at 300. It is part of the definition, and stated.
* **(c) β-phase strip.**
  * Every `Δ_g ≥ 300` row has even-integer `Δ`, and every normalisation is a same-`J` state/reference ratio. So `e^{iπ(Δ+J)}`, which is
    `−1` for odd `J` on even rows, cancels identically.
  * `N = tan²(πτ/2)` is a same-state ratio, so it is phase-free and ≥ 0.
  * No 300/1000 number depends on the strip.
* **(f) Null space, gauge and Lemma N.** These are irrelevant at 300/1000:
  * worst case `≤ 10⁻³` per unit null coefficient (real part);
  * natural gauges `< 10⁻⁴`;
  * the light side is gauge-independent given Lemma N's [NE]-observed regularity.
* **(f) K1 loophole, K1′, R4, R5.** These are inherited technical caveats on validity and swappability (H16). None is specific to 300.

---

## 2. The statement the evidence supports at `Δ_g ≥ 300`

**Result D (route A′, `Ĉ₂′`) at `Δ_g ∈ {300, 1000}`: conditional.**

Let `𝒯` be a unitary CFT₃, parity-even or parity-violating, with stress tensor `T`. Let `Ω^{A′} = Re∫X Ψ_{2,ν}[q²]/(h₋+h₊) dν +
Re∫ψ₃ ½Λ⁶ Ĉ_{3,ν}` be the functional of `H0` §R with the `Ĉ₂′` definition:
* NF02 weights, `k₀ = 27.48509817`;
* `X` the min-norm-gauge Abel weight;
* `Ĉ₃` with the `o_H q` structure, fixed as in S4.

Assume the following.

| label | hypothesis | grade |
|---|---|---|
| **A1** (spectrum) | every exchanged single-trace operator other than `T` has `Δ − J ≥ Δ_g` (all `J`, incl. `0, 1, 2`); no other light single traces (excludes all known holographic CFT₃s) | [A] |
| **A2** (tree level) | `O(1/c_T)` truncation, `c_T ≫ Δ_g²` | [A] |
| **V** (validity + swap, H16) | `Ω^{A′}` annihilates the non-perturbative `⟨TTTT⟩` and commutes with the OPE sum | [D]: K1 = KKSDZ `J₀ ≤ 1` [Q, with their loophole]; K1′ (angular uniformity) not written out [D-sketch]; K2 Regge [Q], endpoints/edges [D], `η → ∞` [P given K1′] (scalar-type kernel) / [D] (spinning); K3 [D] |
| **L0** (light content) | only the identity and the `T` block survive: parity-even `[TT]` at double zeros [D given Q]; parity-odd `[TT]` at double zeros of `cos²` [P given Q + P-sketch]; contacts drop out (unsubtracted) [D given Q] | [D] |
| **L1** (H7) | `−Ĉ′_{2,ν}|_E = 16πG/(ν²+9/4)`, `Ĉ_{3,ν}|_E = 0` | [D]: Lemma R [P given R1–R5; R4 D-sketch, R5 assumed for spinning]; D3 [P]; D2 [D]; §18.1 kernel [D] (validated at large ν and in the scalar limit [DV]); `Ĉ₃` part [P given Q2, Q8, F (F [D/DV])] |
| **L2** (H17) | the `|ĝ₃|²` light coefficient is `k(0)(1+δ)`, `|δ| ≤ 0.5 %` (300) / `0.05 %` (1000), incl. `Ĉ₂′` `δk(0)` | [D-sketch + M] |
| **L3** (H6-lin) | any `E × t₄` light term is `16πG c₁t₄` with `|c₁| ≲ 10` | [I] (not computed; effect `0.43 % c₁` / `0.033 % c₁` on `16L_A`) |
| **P** (heavy positivity, `n = 1`) | every heavy state's per-state action matrix (EE 2×2 ⊕ OO 1×1) is positive semidefinite for all states allowed by A1 | see below |

**Evidence for P**, by region:

| region | evidence | grade |
|---|---|---|
| near gap, even `J ≤ 10`, odd `J = 5, 7, 9`; `m/Λ ∈ [1, 2]` (300), `[1, 1.5]` (1000) | model M (S3) with `K₂` graviton corrections. `n = 1`, min `λ/N_s` 0.2072 / 0.2023 (`+−`, `J = 5`, parity-odd); `t* = (−10.2, +8.3)` / `(−36.8, +39.1)`; insensitive to the `Ĉ₃` de-shift choice, the `Δ = 100, 200` rows, the gauge and the de-shift (all `≤ 6·10⁻³`) | **[NE + M]** |
| near gap, even `J ≥ 12`, odd `J ≥ 11` | base model only. Graviton corrections uncomputed, growing ×2.5 per `ΔJ = 2` | **[M] + [I]** (S1) |
| parity-odd even `J` (OO) | leading order `N = tan² ≥ 0` × `++`; corrections uncomputed | [NE] (LO) + [M] (S5/M9) |
| wedge `m ∈ (Λ, 1.5Λ)`, `J > 20` | scalar-skeleton hybrid scan (this report), `≥ 0.27` / `≥ 0.21` | [M] |
| `m ≥ 1.5Λ`, `b̂ ∈ [20, 91.2)` | hybrid [NE]; exact scalar `S` on a sample [NE + M] | [M] |
| `m ≥ 1.5Λ`, `b̂ ∈ [91.2, 2Λ]` | hybrid positive [certified]; `|S| < F` [I] (scalar calibration `122/m²`) | [certified] for the model; graviton beyond leading Regge [I] |
| `b̂ > 2Λ` | hybrid scan to `b̂ = 1000` at 5 masses [NE]; two-term asymptotic form [D + NE]; `M₁ > 0` [certified] | [M] |
| `m > 2Λ` (300) / `> 1.5Λ` (1000), `b̂ < 20` | hybrid scan (this report), `≥ 0.48`; flat NF02 limit [certified] | [M] |
| between grid points; at the twist edge | hybrid at the edge (this report) | [NE] |

**Conclusion.** Under A1, A2, V, L0–L3 and P:

`√(t₄² + α_T²) ≤ 240·√(16L_A(Λ)/(k(0)(1+δ))) / Λ⁴ (1 + O(c₁t₄Λ²))`
` = 3.13·10⁻⁶ (Δ_g = 300), 2.89·10⁻⁸ (Δ_g = 1000)`.

* `L_A(Λ) = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`, `Λ = √((Δ_g−2)(Δ_g−1))`, AdS units `L_AdS = 1`, `t₄`, `α_T` as in `DICTIONARY_T4.md`.
* Uncertainty: `±0.4 %` / `±0.04 %` from L2 and L3 with `|c₁| ≤ 1`.

**Overall grade at `Δ_g ≥ 300`.** The implication P ⇒ bound is **[D]**, and the constants are [DV]. P is **[NE + M] for `J ≤ 10` near
the gap, and [M] + [I] for gravitons everywhere else** (near-gap `J ≥ 11`; every state at `b̂ ≥ 20`; graviton beyond-leading-Regge
corrections).

The statement is therefore **conditional on P beyond `J = 10`**, not merely "[NE + M]". It holds at the two computed gaps. "For all
`Δ_g ≥ 300`" is an extrapolation (M10). Shown: `Δ_g = 300` and `1000`. Not shown: any other `Δ_g`, and `Δ_g = 100`, which stays [I].

---

## 3. Proposed replacement tables for `RESULT_D.md`

### 3a. Headline numbers (replaces the §2 table)

| `Δ_g` | `Λ` | `16L_A` (exact) | `A log Λ + B_AdS` | **`√(t₄²+α_T²)` bound** (exact `p`-integral, `k = 1`) | `Ĉ₂′` `δk(0)` → bound factor | H17 band (bound) | `c₁` sensitivity (bound, per unit `c₁`) | `n = 1`? (binding state, min `λ/N_s`, `t*`) |
|---|---|---|---|---|---|---|---|---|
| 100 | 98.499 | 7508.5 | 7320.0 | **2.21·10⁻⁴** | `−7.20·10⁻³` → ×1.0036 (2.22·10⁻⁴) | ±2 % | ≈ ∓1.7 % (scaled estimate) | **[I]**: main protocol `n = 1` (even `J = 4`, 0.2206); fails with extra odd `m`-rows (`J = 9`, `1.2Λ`, `−0.55`) |
| 300 | 298.50 | 10734.6 | 10670.9 | **3.13·10⁻⁶** | `−7.66·10⁻⁴` → ×1.0004 | ±0.2 % | ∓0.22 % | **[NE + M] for `J ≤ 10`**: `+−`, `J = 5` (parity-odd), 0.2072, `t* = (−10.2, +8.3)`; **[I] for `J ≥ 11` and `b̂ ≥ 20`** |
| 1000 | 998.50 | 14339.6 | 14320.4 | **2.89·10⁻⁸** | `−6.83·10⁻⁵` → ×1.00003 | ±0.02 % | ∓0.017 % | **[NE + M] for `J ≤ 10`**: `+−`, `J = 5`, 0.2023, `t* = (−36.8, +39.1)`; [I] for `J ≥ 11`, `b̂ ≥ 20` |

* `A = 3022.34`, `B_AdS = −6552.66` (flat NF02, `k(0) = 1`).
* The `Ω^U`/M02 "D-strong" row is **withdrawn**: B′ objects, not re-derived for the stripped dictionary.

### 3b. Hypotheses (replaces the §4 table; route A′ with `Ĉ₂′`)

| # | hypothesis | grade | evidence / what would settle it |
|---|---|---|---|
| **H0** | **Definition of `Ω^{A′}`**: `Ĉ₂′ := 𝒜_h⁻¹[Ψ_{2,ν}[q²]]` (smeared: `Re∫XΨ̂`, Abel weight, min-norm gauge); `Ĉ₃ := Re[Ψ_{3,ν}[o_H q]/N̂₃]`; NF02 `ψ₂, ψ₃`; `Re` | definition; ingredients: CLS-D functionals [Q]; `d = 3` specialisation, `q`, `o_H` [DV at crystal, D covariant]; `h_∓ = N_com R_∓` [P given §18.1 (D)]; Abel weight [D + NE]; gauge = convention | **`Ĉ₃` definition inconsistent** (S4): fix to no de-shift and remove the `k = 3` shift in the chains, after checking `Ĉ₃`'s spinning leading-Regge image |
| H1 | heavy action on `⟨TTTT⟩` blocks at leading order = stripped flat densities | [Q] (CLS-D crystal) + [DV] (`d = 3`, both helicities, both parities, `J = 4–10`; odd `J = 5, 7, 9`) | merged with H2 |
| H2 | bulk-point dictionary (stripped `B^{(1)}_{2,3}`, `(1−x²)^{−4}`) | [Q + DV] | — |
| H3 | odd-`k` moments | **not needed** (only `k = 2, 3`) | — |
| **H4** | graviton-specific `O(1/m²)` corrections near the gap (model M of S3) | **[NE + M] for even `J ≤ 10`, odd `J = 5, 7, 9`** at `Δ_g = 300, 1000` (`K₂`, rows `Δ = 100–1600`; direct checks only on `Δ = 100`); **[I] for `J ≥ 11`**, where the computed corrections grow ×2.5 per `ΔJ = 2` | `Γ` for `J = 11–20` on `Δ = 400–1600`; direct points on `Δ = 400` |
| H5 | near-gap tolerance convention | not needed | — |
| **H6** | (F) flat limit of the light side; (P) contacts drop out; (lin) `E × t₄` term | (F) [Q + D]; (P) [D given Q]; **(lin) [I]**, effect `0.43 % c₁` (300) | compute `c₁` |
| **H7** | `−Ĉ′₂|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` | **[D]** (Lemma R [P given R1–R5, R4 D-sketch]; D3 [P]; D2 [D]; §18.1 kernel [D], validated at large ν [DV]); `Ĉ₃` part [P given Q2, Q8, F] | independent small-ν check of `h₋ + h₊` (M5) |
| H8 | Einstein Regge spin 2 | [Q] | — |
| H9a/H9b | forward-term / smeared even-`k` contact content | **not needed** (unsubtracted `k = 2, 3`) | — |
| H10 | EFT power counting of contacts | **not needed at tree level** (contacts drop out) | — |
| **H11** | spectrum: all single traces `≠ T` have `Δ − J ≥ Δ_g` (all `J`); tree level, `c_T ≫ Δ_g²`; `[TT]` double zeros (even [D given Q], odd [P given Q + P-sketch]) | [A] / [D] | fix §1 wording (M1) |
| **H12′** | parity-odd sector | odd `J`: shape [NE], `N = tan²(πτ/2) > 0` [NE + P given Q], EO cancels [DV], `J = 5, 7, 9` corrected [NE + M]; odd `J ≥ 11` [I]; **even-`J` OO block uncorrected [M]** | p6 run at even `J` |
| **H13′** | positivity off the corrected set | near-gap grid + continuity [NE]; gap edge [NE] (this report); wedge `m < 1.5Λ`, large `J` [M] (this report); `b̂ ∈ [20, 91.2)` [M]; `b̂ ∈ [91.2, 2Λ]` hybrid [certified]; `b̂ > 2Λ` [M]; `b̂ < 20` at `m > 2Λ` [M] (this report); **graviton beyond leading Regge [I]**; `|S| < F` off sample [I] | `O(r)` spinning Regge kernel (S2) |
| H14 | object identity | **resolved by definition** (`Ĉ` is the CLS-D object; the base model is CHMRS, identical for scalars) | — |
| H15 | `|ĝ₃|²` absent from `Ψ₅, Ψ₆, Ψ₇` | **not needed** | — |
| **H16** | validity and swappability (I14) | **[D]**: K1 [Q, KKSDZ loophole], K1′ [D-sketch], K2 [Q/D; P given K1′ at `η → ∞`], K3 [D] | write out K1′ |
| **H17** | `k_AdS − k(0)` (incl. `Ĉ₂′` `δk(0)`) | [D-sketch + M]; ±0.2 % (300), ±0.02 % (1000) on the bound | spinning `R³` exchange to `O(r²)` |
| **H18** (new) | `Re` prescription and phase convention | definition; legitimacy [D]; conditions (a)–(c) hold [DV/NE] | — |
| **H19** (new) | null-space gauge; Lemma N | gauge = convention; Lemma N [P given NE-observed regularity]; effect at 300/1000 `< 10⁻⁴` | — |
| **H20** (new) | `Ĉ₃` leading-Regge heavy image `σ_h c₃(ν)P_l` (HN₃) | [Q] (scalar-type kernel) + [NE] (crystal: `σ_h` pattern, `m⁻⁶`) | split-representation computation with `o_H` (S4) |

### 3c. "What is established" for route A′ (replaces §3b–§3c, which refer to the B′ objects)

| ingredient | status | file |
|---|---|---|
| NF02 flat positivity (stripped), min 0.200 | [certified] (one `+−` corner float, empty at 300) | `NOFORWARD_FUNCTIONAL.md` |
| `A = 3022.34`, `B_AdS = −6552.66`, bounds | [DV] (reproduced here) | `H7_I14_LIGHTSIDE.md` v1 (g); this report |
| `Ĉ₂′` leading-Regge kernel `h₋P_{l−1} + h₊P_{l+1}`, `h₀ ≡ 0`, `++ = +−`; closed form | [P given §18.1 (D)] + [DV] | `H7_I14_LIGHTSIDE.md` §18, §23.1 |
| scalar-skeleton near-gap action at `Δ_g = 100` | [certified] (Theorem A, grid) | `CERTIFIED_NEARGAP_A.md` |
| scalar-skeleton hybrid at large `b̂`, `b̂ ∈ [91.2, 2Λ]`; `M₁ > 0` | [certified] | `ROUTEA_LARGE_B.md` |
| graviton near-gap corrections `J ≤ 10` at `Δ_g = 300, 1000` | [NE + M] | `H0` §23, `PARITY_ODD_ROUTEA.md` §6, §9; `H7_I14` §20.3, §23.3 |
| EO interference cancels; OO sign `N = tan² ≥ 0` | [DV]; [NE + P given Q] | `PARITY_ODD_ROUTEA.md` §3, §12 |

---

## 4. Re-runs made for this report (one core, read-only imports)

1. **Constants.** `omegau/nf_common.NF`, scipy `quad`.
   * Results: `k(0) = 1.0000`, `A = 16ψ₂′(0) = 3022.3369`, `16L_A = 7508.515 / 10734.605 / 14339.589`; bound (exact `p`) =
     `2.2094·10⁻⁴ / 3.1321·10⁻⁶ / 2.8913·10⁻⁸`; log form `2.181·10⁻⁴ / 3.123·10⁻⁶ / 2.889·10⁻⁸`.
   * `∫ψ₂ = 10.386`. Per unit `c₁`, a linear `t₄` term shifts `16L_A` by `4.3·10⁻³` (300) and `3.3·10⁻⁴` (1000).
2. **Wedge scan** (`h13a_common.hybridA`, `h13_common.flat_rows`, `omegau/common.state`).
   * Grid: `m/Λ ∈ {1.005, …, 1.5}` (15 values), every allowed `J ≥ 4` up to the twist cap (`J ≤ 186` at 300, `≤ 624` at 1000), both
     channels.
   * All positive. Min 0.2056 (300; `+−`, `J = 5`, `1.02Λ`), 0.2009 (1000; `+−`, `J = 5`, `1.005Λ`). The per-mass minima are always at
     `J = 4–6`.
3. **`b̂ < 20` above the exact-skeleton sets.** Same tools, `m/Λ ∈ {1.6, 2, 2.5, 3, 4, 6, 10, 30}`. Min 0.482 (`++`, `J = 8`, `1.6Λ`) at
   both gaps.
4. **Twist strip and gap edge.** `hybridA` at `α = Δ_g − 4, −3, −2`, `J = 3–40`. Min 0.2021 (300) / 0.2005 (1000), `+−`, `J = 5`,
   `α = Δ_g − 4`. At the code's edge (`α = Δ_g − 2`) `J = 5` gives 0.2043 at `1.0166Λ` (300).
5. **`Ĉ₃` de-shift sensitivity.** `v3_heavy_impact` (k = 2 → `σ_T`), with `l9_analysis.shift` for `k = 3` either kept or set to zero;
   `l9_neargap` `C23H_K2`, 300/1000.
   * `n = 1` either way.
   * 300, `J = 6 / 8 / 10` minima `0.3049 / 0.4332 / 0.6017` (kept) vs `0.3086 / 0.4275 / 0.6037` (none).
   * Max `|ΔF_hh|/N_s` at `J = 10`: `0.072` vs `0.021`.
   * `t* (−10.2, +8.3)` vs `(−10.2, +8.7)`.
   * The per-state off-diagonal table of S1 is from the same run.
6. **Row sensitivity.** Same chain with the `Δ = 100, 200` rows removed (corrections held at `Δ = 400` below `m(400)`): `n = 1`,
   min 0.2054 (300) / 0.2021 (1000), `t* (−5.6, +5.1)` / `(−36.8, +39.1)`, `J = 10` off-diagonal max 0.170 (300).
7. **Base margins for `J ≤ 20`** (`spf_common.load_states`, f6 base). At 300, per-helicity minima for `J = 11–20` lie in `0.68–1.00`.
   At 1000 they lie in `0.41–1.00`; the lowest is `J = 11` `+−` at 0.411.

Outputs went to the session scratchpad only (`…/scratchpad/ng/cache/`). No project cache was written.
