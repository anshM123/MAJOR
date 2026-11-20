# Hostile review of Result D with `Ω^EE` (`RESULT_D_STATEMENT.md`, state of 20 Nov 2026)

*20 Nov 2026, referee agent. Target: `ads4/RESULT_D_STATEMENT.md` (commit `50d2260`). Read in full: the statement, `EE_FUNCTIONAL.md`,
`EE_FOUNDATIONS.md`, `EE_NEARGAP.md`, `EE_DG300_SMALLNU.md`. Read in the relevant parts: `EE_DG1000_LARGEJ.md`, `SUBLEADING_REGGE.md`,
`SMOOTH_BOUND.md`, `FOLD_WINDOW_LARGEJ.md`, `ROUTEA_LARGE_B.md`, `REVIEW_ROUTEA_2026-11.md`, `H0_SPINNING_FUNCTIONAL.md` §15/§18,
`H7_I14_LIGHTSIDE.md` (Lemma R), `PARITY_ODD_ROUTEA.md` §0–§1, and the update box of `RESULT_D.md`. Code was read in
`spinfunc/{w2_chain, ee1000_chain, q3_analysis (_strip), l9_analysis (a_h fit), z3_nodes, sm_chainfix, st_odd_despike}.py`.
No project file other than this one was edited and nothing was committed. Every re-run wrote its output under a new `zzrev_*` tag,
and those files were moved to the session scratchpad (`…/scratchpad/rev/`), so the project cache is unchanged.*

Severity scale: **fatal** (the implication or the headline numbers are wrong) / **serious** (a grade or coverage claim in the
headline must change) / **minor** (a local error or an unflagged hole; no verdict changes) / **presentational**.

---

## 0. Verdict (read this first)

1. **No fatal error was found.** The step "per-state PSD of `Ω^EE` on every twist-allowed heavy block ⇒
   `√(t₄²+α_T²) ≤ 240√(16L_A/(k(0)(1+δ)))/Λ⁴`" is sound at the stated grade [D].
   * I re-derived the constants independently: `16L_A = 7508.51 / 10734.61 / 14339.59`, `A = 16ψ₂′(0) = 3022.337`, `k(0) = 1.0000`,
     bounds `2.2094·10⁻⁴ / 3.1321·10⁻⁶ / 2.8913·10⁻⁸`.
   * Nothing on the light side silently uses `q²`. The pole uses EE's own (HN), `k(0)` uses only the crystal images, and `δk(0)` uses EE's
     own `X`. The one structure-specific light input that was never recomputed for `EE` is H17's `Π̂₄` admixture. It is small (§3, m4).
2. **The functional is legitimate but was selected post hoc.**
   * Both the structure (`q² → EE`) and the gauge (canonical → min-norm) were changed after heavy-state failures had been seen. The
     chronology is in §1.
   * This cannot invalidate the bound: any fixed solution of the Abel equation is a valid functional, and positivity of any one of them
     suffices.
   * It does make a minority of the "computed" states in-sample evidence: `J = 18, 20, 24, 28` at `Δ_g = 300`.
   * Genuine out-of-sample evidence exists and passes: everything computed after 23:41 on 19 Nov (`J = 32–48` and odd `J = 17–31` at
     1000; the direct small-`ν` nodes at 300).
3. **"Every computed state passes" is true in the adopted protocol, but the statement over-grades several of those states.**
   * `J = 28, 32` at 300 and `J = 48` at 1000 rest on `K₂` beyond the project's own convergence criterion (S1).
   * Every computed state uses a pure, unvalidated extrapolation on the light-cone-edge window `x ∈ (0.8, Λ/m]`. That window can move
     the gap state by `O(1) N_s` (S2).
   * The `J ≤ 22` evidence at 300 comes from the protocol that was diagnosed as an artefact for `J ≥ 24`. Post-fix, `J = 22` passes
     there with `t* = 1.1` (S4).
   * An independent, unreported thimble repeat lowers `J = 28` at 1.1Λ (300) from `+1.57` to `+0.69` (S5).
4. **Honest claim** (§6): a conditional statement.
   * At **`Δ_g = 1000`** positivity is computed [NE + M] for the EE block on a near-gap set covering roughly the first 40 spins (with
     interpolation in `J` and an extrapolated edge window). Everything else is model-level, and the large-spin near-gap region `H_LJ`
     is open with a tolerance of `1–30 c₂` against typical corrections of `10²–10⁴ c₂`.
   * At **`Δ_g = 300`** the computed set is materially smaller: even `J ≤ 24` and odd `J ≤ 15`. `J = 28, 32` should be regraded to [M]/[I].
   * The statement's own framing ("not a theorem at any `Δ_g`"; 300 weaker than 1000; 100 not established) is correct and should stay.

---

## 1. Chronology of design decisions (from `git log` and the documents)

| time (UTC) | event | triggered by |
|---|---|---|
| 18–19 Nov | route A′ with `q²` fails near the gap at 300 (`J = 16–20`), and at 1000 from `J ≈ 32` (model) | S1_S4, LARGEJ_OFFDIAG |
| 19 Nov 16:58 | `EE` adopted (`beaa80b`). The alternatives tried were `EH` (null at the crystal) and `EE + λq²` (LP infeasible). It was tested immediately on the two failing states `J = 18, 20` at 300. | the `q²` failures; motivated by the leading-Regge kernel (off-diagonal `≡ 0`) — a genuine theoretical reason |
| 19 Nov 17:25 | `EE_FOUNDATIONS` §3.3 **adopts the canonical (CHMRS-class) gauge as the definition** | analytic (canonical class exists) |
| 19 Nov 23:41 | `EE_NEARGAP` §4: the canonical gauge gives `n = 1.55` at 300 (spline); the doc recommends min-norm "as the stated convention" | heavy-state signs |
| 20 Nov 04:52 | `EE_DG300_SMALLNU`: canonical `−37.5 / −33.1` at `J = 24 / 28` on the best data; min-norm passes | heavy-state signs |
| 20 Nov 07:27 | `eval_jacobi` fix | a failure (`J = 32` at `−3.63`) |
| 20 Nov 08:56 | curvature double count fixed | *not* failure-triggered; it moved λ **down** (good sign) |
| 20 Nov 12:57 / 22:02 | engine `ε²4^J` fix | an anomaly (`J = 48` offset) |
| 20 Nov 22:25 | odd-row Hampel added | a failure (`J = 29` at `−1.20`) |

**Assessment.**
* Min-norm was the original default (u3, before any heavy test). The reversal is a return to it, not a search over many gauges. The
  selection is between two candidates.
* The structure change has an independent derivation (`o_EE ≡ 0` at leading Regge [P]). Neither choice is a tuning parameter fitted to
  data. The look-elsewhere effect is small.
* Debugging was, however, asymmetric. Every bug hunt except the curvature fix started from a *negative* number, and no audit of
  suspiciously favourable numbers was made. Examples that should have been audited:
  * `J = 31` at 1000: the gap state is `+0.21` from two held `x ≥ 0.775` nodes of size `±2000–3600`;
  * `J = 48`, variant `x ≤ 0.6`: a held-extrapolation correction of `+4.1 N_s`.
* None of this affects the implication, because any fixed functional that is positive gives the bound. It affects only how far the
  computed states can be trusted to stand for the uncomputed ones.

---

## 2. Spot checks re-run for this review (one core; outputs in `…/scratchpad/rev/`)

| check | command (abridged) | result | agrees? |
|---|---|---|---|
| constants | `omegau/nf_common.NF`, `quad` | `16L_A`, `A`, `k(0)`, bounds as in §0.1 | yes |
| binding states | `sm_chainfix.py w2 300 1000` (`Q3_JEVEN=4 Q3_JODD=5`, `SM_CURV=row`, rows 200…1600) | 300: `+0.2052` (`+−`, `J = 5`, 1.02Λ); 1000: `+0.2009` (hybrid, 1.005Λ), `+0.2021` on the skeleton grid | yes |
| statement §9.1, raw | `sm_chainfix.py ee1000`, odd `J = 17–31` | `J = 29`: `−1.2031` at 1.0286Λ (shape-fit residual 0.24/0.32, against 0.03–0.08 for the other `J`) | yes |
| statement §9.1, despiked | `st_odd_despike.py ee1000` | `+0.839 / 0.916 / 0.942 / 0.887 / 0.918 / 0.954 / 0.979 / 0.995`; dropped nodes `J = 19`: 0.775, 0.8; `J = 23`, `29`: 0.8 | yes |
| 300, `J = 24, 28` (curvature- and Jacobi-fixed) | `sm_chainfix.py z4 300`, `Z4_EST=best`, rows 336/500/1600 | `+0.8204 / +0.8133` (1.5Λ) | yes (SMOOTH_BOUND §4) |
| **300, `J = 28` with the independent thimble repeat** | same, `Z3_PREFER=x,s2` | **min `+0.6905` at 1.1Λ** (was `+1.570` there). `Γ_c(ν = 5.0)`, `++`: `+10073+769i` (main) vs `−1108+3610i` (repeat, 2× points) | new (S5) |
| edge-window sensitivity, 1000 | `st_odd_despike.py ee1000`, `EE1000_XMAX = 0.8 / 0.7 / 0.6`, `J = 27, 31, 36, 40` | `J = 36`: `+0.967 / +0.922 / +0.534`; `J = 40`: `+0.979 / +0.979 / +0.650`; `J = 31`: `+0.995 / +0.979 / +0.763` | new (S2) |
| `K₂` order ratios on the 300 rows | `M₀, M₂, M₄` of the cached records (whole-matrix max, median over `x`) | `|M₄|/|M₂|`: `(336,18)` 0.44, `(336,22)` 0.51, `(336,24)` 0.47, **`(336,28)` 0.79, `(336,32)` 1.20**; `(1048,48)` 0.70, `(1064,64)` 0.77 | new (S1) |
| phase strip on odd rows | `ee1000_p3_EE_J19.pkl` | state `(1019, 19)`: `Im logK ≡ 0`; reference `(2·10⁴, 19)`: `Im logK ≡ π` | new (m3) |
| non-finite drops | grep of chain outputs | 0 everywhere | yes |
| reference-key collisions in the 300 direct caches | keys of `z2*_pts_*_J{24,28,32}.pkl` | one `Δ` per `J`; no collision | yes |

---

## 3. Findings

### Fatal

None. Specifically:
* The bound needs only the light side, which is functional-specific in exactly two places: the pole via EE's own (HN), and `δk(0)` via
  EE's own `X`. Both were recomputed for `EE`.
* "`EE ≡ q²` only on the physical polarisation plane" does not touch the implication. The flat limit enters only the *evidence* for
  positivity (as the base model), never the inequality.
* The 2-dim null space and the min-norm convention give a well-defined fixed functional, grid-converged to `3·10⁻⁵`.

### Serious

**S1. States counted as "computed" rest on `K₂` beyond the project's own convergence criterion.**
* *Location.* Statement §0.3, §5.1 row "even `J = 48`", §5.2 rows "`J = 28`" and "`J = 32`", §6.2. Its own H_M comment says "`K₂` usable
  only for `J²/m ≲ 1.5`" (EE_DG1000_LARGEJ §3). FOLD_WINDOW §3 says "the engine fix moves `J = 48` (2.2) into the usable range, not
  `J = 64` (3.9)".
* *Why it matters.*
  * At 300 the `Δ = 336` rows have `J²/m = 1.73 / 2.35 / 3.08` at `J = 24 / 28 / 32`. Direct data replace `K₂` only for `x ≤ 0.1`
    (`ν ≤ 33`); everything above is `K₂`.
  * `(336, 24)` was validated against direct state values (`≤ 1 %`, `x ≤ 0.4`). `(336, 28)` and `(336, 32)` were not validated at all.
  * The crude order ratio of §2 shows `|M₄| ≥ |M₂|` on `(336, 32)` (1.20), worse than `(1064, 64)` (0.77), which the project declares
    non-converged. `(336, 28)` (0.79) sits beside `(1048, 48)` (0.70).
  * `J = 48` at 1000 has no post-fix direct-vs-`K₂` comparison. Its `+0.932` rests on shape residuals and the vanished offset, and its
    variants span `+0.73 … +0.99`, including a `+4.1 N_s` held-extrapolation correction.
* *Fix.* Regrade `J = 28, 32` (300) to **[M]** and move them into `H_LJ(300)`, or compute direct nodes on `(336, 28/32)` for
  `x = 0.15–0.6` (≈ 15 min per node; the z2 machinery exists). State `J = 48` (1000) as "[NE + M] at the `K₂` convergence limit,
  unvalidated".

**S2. The light-cone-edge window `x ∈ (0.8, Λ/m]` is a pure extrapolation for every computed state and is not a named hypothesis.**
* *Location.* `ee1000_chain` (`EE1000_XMAX = 0.8`; the shape model and spline hold `Γ` beyond the last node up to `x_max = Λ/m ≈ 0.90–0.99`
  at the gap). The same holds for the w2/z4 chains (rows end at 0.8). Statement §5 lists "no graviton data at `x > 0.8`" only under
  `H_LJ` (`J > 48`).
* *Why it matters.*
  * The `+−` envelope `(2−y)/(1−y)⁴` makes this window heavy at the gap. One held node at `x = 0.8` (`J = 29`, `Γ = −1.8·10⁴`) moved the
    gap state by `−2.2 N_s`.
  * Truncating the data at `x = 0.6` or `0.7` and holding beyond (§2) moves `J = 36 / 40 / 31` by `−0.43 / −0.33 / −0.23 N_s`.
  * Neighbouring `x ≈ 0.8` node values are `±2000–3600` at `J ≈ 30`. The failure threshold of a uniform held value is roughly
    `−8000` at `J = 29` (scaling the `−1.8·10⁴ → −2.2 N_s` spike), presumably falling like `N_s ∝ b̂⁻³` with `J`. That is a margin of only ×2–4 on an uncomputed region.
  * `K₂` is "validated to 0.8" only.
* *Fix.* Add an explicit **`H_edge`** to `H_M` ("`Γ` on `x ∈ (0.8, Λ/m]` is within the held value ± T_edge(J)"), with the tolerance
  table, for every computed `J`. Add direct nodes at `x = 0.85, 0.9` on one row per `Δ_g` (`J = 24`/300, `J = 40`/1000).

**S3. Post hoc selection of the functional and protocol; asymmetric debugging (§1).**
* *Location.* Statement §2 ("gauge: definition"), §7 ("canonical gauge … reversed"), §8.1.
* *Why it matters.*
  * The statement presents min-norm as a convention and says "any gauge gives a valid functional". That is true and sufficient for the
    implication.
  * But the statement then cites `J = 18, 20, 24, 28` (300) as positivity evidence, and those states selected both the structure and
    the gauge.
  * The canonical-gauge functional is equally valid, and on the best available data it is *negative* at `J = 24, 28` (`−37.5 / −33.1`).
    The two gauges differ by a null-direction functional `N` with zero light side, so `Σ_heavy p·N = 0`. `N` is O(10) `N_s` on near-gap
    300 states and is driven entirely by `ν ≲ 6`, where the data have 5–100 % errors.
  * So positivity at 300 near the gap is a statement about how the functional weights a window the data do not resolve.
* *Fix.*
  * Freeze the definition (min-norm, NF02, `EE`) with a date.
  * Label each computed state in-sample (computed before 19 Nov 23:41) or out-of-sample.
  * Add one sentence: "the gauge was chosen after the canonical gauge failed on `J = 24, 28` (300) in the then-available data".
  * Audit the favourable outliers listed in §1 with the same tools used on the negative ones.

**S4. The `Δ_g = 300` evidence for `J ≤ 22` (and odd `J ≤ 15`) comes from the protocol that was diagnosed as an artefact for `J ≥ 24`,
and is quoted pre-fix.**
* *Location.* Statement §5.2 row 1 ("`[NE + M]` in every variant … worst even `J = 22` `+0.15` (convention A)"); §2 "What the chains
  compute … the difference is `≤ 3·10⁻³ N_s` at 300".
* *Why it matters.*
  * `J ≤ 22` use the `K₂` node at `x = 0.05` (`ν ≈ 17`) plus a natural-spline extrapolation below it, with local `ψ₂` weights
    (convention A). For `J = 24, 28, 32` this is exactly what EE_DG300_SMALLNU showed to be wrong in sign.
  * In the curvature- and Jacobi-fixed run (`cache/sm_chainfix_w2_row_kappa.txt`), convention A gives `J = 22`: **`+0.123`, `t*_J = 1.1`**;
    `J = 18 / 20`: `+0.472 / +0.473`. The statement quotes the pre-fix `+0.147`, `+0.476`, `+0.495`.
  * The observed `K₂` small-`ν` error had the unfavourable sign (`J = 18`, `x = 0.075`: `K₂` `−364` vs direct `−65`). The `J ≤ 22` values
    are therefore probably conservative, but that was not shown.
  * The local-vs-min-norm difference at 300 is not `≤ 3·10⁻³ N_s`. EE_NEARGAP §4 gives `0.023 / 0.032 / 0.066 / 0.096 / 0.046` at
    `J = 18 / 20 / 22 / 24 / 15 (odd)`, all in the favourable direction. The `3·10⁻³` is the `q²`-era review number (M7, `J ≤ 10`). So
    the reported `J ≤ 22` values are not values of the defined (min-norm) functional.
* *Fix.*
  * Run the z4 protocol (direct n3 nodes `ν = 8.3–33` + non-local min-norm `X` below `ν = 25`) for even `J = 16–22` and odd `J = 13, 15` at
    300: about 8 nodes × 15 min per `J`.
  * Until then, write "`J ≤ 22`: [NE + M] in convention A (superseded small-`ν` treatment), margin `t* = 1.1` at `J = 22`".
  * Correct the `≤ 3·10⁻³` sentence.

**S5. `J = 28` at `Δ_g = 300`: a finished independent repeat of the decisive thimble nodes was never analysed, and it moves the formerly
failing state by `−0.9 N_s`.**
* *Location.* Statement §5.2 row "`J = 28`" (`1.1Λ` `+1.45 ± 0.58`); EE_DG300_SMALLNU §10 (`z2_queueS2.sh`, "pending"). The repeat finished
  at 06:36 on 20 Nov (`cache/z2s2_pts_*_J28.pkl`, `cache/z2_queueS2.log`) and is cited nowhere.
* *Result (this review).*
  * With `Z3_PREFER=x,s2` the curvature- and Jacobi-fixed chain gives `J = 28` min **`+0.6905` at 1.1Λ** (was `+1.570`; main minimum
    `+0.8133` at 1.5Λ).
  * The `ν = 5.0` node goes from `Γ_c,++ = +10073+769i` to `−1108+3610i`, and C3H from `−36` to `+4717`. The repeat has 2× the points and
    scalar errors of 4–29 %.
  * Positive either way, but the statement's "`+1.45 ± 0.58`" is an upward fluctuation, and the scatter between two seeds (`≈ 10⁴ c₂`) is
    larger than the quoted error.
* *Fix.* Report both runs; grade `J = 28` (and `J = 32`, which has no `ν < 8.3` data at all) as [M] for the `ν < 8.3` window, as the
  statement already half-does.

**S6. The scale of `H_LJ` is understated by the headline.**
* *Location.* §0.3–§0.4, §6.2 ("holds given `H_M`, `H_OO`, `H_R` and `H_LJ`"; "decidable by computation in principle").
* *Why it matters.*
  * At 1000, `H_LJ` covers every near-gap even `J ≥ 50` up to the twist-allowed `J ≈ 344` (`m ≤ 1.3Λ`), every odd `J ≥ 33`, and seven
    unrowed even `J` below 48. At 300 it covers every near-gap even `J ≥ 34` (to `J ≈ 103`), `J = 26, 30`, every odd `J ≥ 17`, and the
    `ν < 8.3` window.
  * Its tolerance for a uniform smooth correction is `1–30 c₂` (1000) and `1–25 c₂` (300). Typical graviton corrections are
    `10²–10⁴ c₂`, so `H_LJ` asserts a cancellation of the non-oscillating part to `10⁻²–10⁻⁴`.
  * The supporting evidence is the exact scalar skeleton (`≤ 4·10⁻⁴`) and a [D-sketch] saddle argument. The documents say it is
    "undecidable with the present direct method" at `J ≥ 150` (1000) and `J ≥ 100` (300); the only `J = 64` datum at 300 (`−51 ± 43`
    against 70) is already at the edge.
  * "Decidable in principle" is true but needs "not with any method the project has".
* *Fix.*
  * In §0.3 and §6.2 add one clause: "`H_LJ` covers most of the near-gap spin range, with a sub-percent tolerance, and is beyond present
    numerical methods for `J ≳ 150` (1000) / `≳ 100` (300)".
  * Keep the [I] grade.

**S7. `H_OO` is a whole helicity sector with no data beyond leading order, and "every computed state" silently means "the EE block".**
* *Location.* §0.3 ("Every computed state passes"); §5 rows "parity-odd even-`J` (OO) block"; §4 `H_OO`.
* *Why it matters.*
  * `X^O_same` couples like `++`. At 300 near the gap, the `++` correction in `K₂` was up to `1.4 N_s` (convention A, `J = 28`), and it
    took direct small-`ν` nodes to make it benign.
  * The OO block has had neither a `K₂` row nor a direct node at any `J ≥ 2`. So for even `J` the per-state test is incomplete on
    *every* computed state.
  * The OO sector decouples, so this does not contaminate the EE numbers, but positivity requires it as well.
* *Fix.* Say "EE block" in the headline. Compute one OO `K₂` row (`J = 18` or `24`, `Δ = 336`, `p6`-type) before claiming any even-`J`
  state at 300 as computed.

### Minor

**m1. Despike rules: coverage and reproducibility.**
* The Hampel filter in `w2_chain._grid` acts on keys `(kind, '++'|'+−')` only. **The off-diagonal keys are never despiked.** They matter
  at 300 for `J ≥ 22` (`|dF₊₋|` up to 0.036).
* `bad[0] = False` exempts the first node, i.e. `x = 0.05`, which is exactly the node with the documented wrong-sign `K₂` (EE_DG300 §0).
* The odd-row despike exists only as an in-process wrapper (`spinfunc/st_odd_despike.py`; the scratchpad copy is identical). The
  statement's §9.1 results are not reproducible from `ee1000_chain.py` alone.
* *Fix.* Move the odd and off-diagonal despike into the chain modules, behind a flag that defaults on, and print the drop log in every
  output.

**m2. Unflagged mass holes at the twist edge (300).**
* The z4/w2 chain states start at 1.1Λ. The twist-allowed edge is `1.077Λ` (`J = 24`) and `1.090Λ` (`J = 28`), so `m ∈ [edge, 1.1)Λ` is
  unsampled there. Only `J = 32`'s `[1.102, 1.15)` is flagged.
* Below the lowest row (`m = 333.6 = 1.118Λ`), the correction is *held*. Convention A failed exactly at the lowest mass, and the
  correction grows as `m` falls.
* More generally, the hybrid base oscillates in `m` with period `≈ 0.03–0.05Λ` at `J ≳ 17` (e.g. `J = 17`, 1000: `0.992 / 0.839 / 0.970`
  at `1.017 / 1.05 / 1.075Λ`). The grid `{1.02, 1.05, 1.1, …}` does not resolve it.
* This is harmless at 1000: the base is floored by the flat certificate (0.200), and the corrections are `≤ 0.02`. At 300, where
  `dF ~ O(1) N_s` near the gap, it is not demonstrated.
* *Fix.* Add the edge masses to `Q3_EXTRA_M`, and a `0.01Λ` mass grid for `m ≤ 1.3Λ` at the `J` with large corrections.

**m3. The phase-strip claim in the hypothesis table is wrong for odd `J` (the numbers are right).**
* §4 row "Re" says "no `Δ_g ≥ 300` number depends on the strip (even-integer `Δ`, same-`J` ratios)".
* On the ee1000 odd rows (`Δ = 1000 + J`), the state has `Im logK ≡ 0` and the reference `(2·10⁴, J)` has `Im logK ≡ π` (checked,
  `J = 19`). Without the strip, every odd-`J ≥ 17` correction at 1000 would change sign.
* The code strips correctly and asserts it (`q3_analysis._strip`), so nothing changes. But the strip *is* load-bearing there, and the
  sentence should say so.

**m4. H17 and `δk(0)` for `EE`: the structure-specific part was not recomputed.**
* The `|ĝ₃|²` light term has Regge spin 4, so its AdS correction involves the `Π̂₄` admixture `c₄^{(2)}` of the specific functional
  (H0 §18).
* For `Ĉ′_EE = c(ν)⁻¹𝒜_g⁻¹Ψ_EE` that admixture was never derived. `c(ν)⁻¹` has poles at `±i/2, ±3i/2`, inside the strip used by
  Lemma N, so "the `|ĝ₃|²` term is invisible to the null directions (given its C-level strip analyticity)" is an assumption.
* The ±2 % / ±0.2 % / ±0.02 % band is the `q²`-era heuristic (`|a| = D × 36`), carried over unchanged. Separate `t₄²` and `α_T²`
  coefficients at `O(1/Λ²)` are not excluded either; the bound is then on a weighted norm.
* Even a ×10 error in the band moves the 300 bound by only about 2 %.
* *Fix.* State "H17 for `EE`: model carried over from `q²`; `Π̂₄` admixture of `Ĉ′_EE` not derived".

**m5. `c₁`.**
* The headline quotes "`|c₁| ≤ 1`" (±0.4 %). There is no argument for that normalisation, and the `Λ⁻²` scaling of the linear term is
  asserted, not derived.
* Plausibly `c₁ = 0` at leading Regge for `EE`: its kernel is helicity-diagonal, and the linear `ĝ₃` term lives in the flip amplitude.
  That would be worth one line.
* *Fix.* Quote the bound with `c₁` as an explicit parameter ([I]).

**m6. Fold-window evidence for `H_LJ` is not reproducible from the current caches.**
* `(1064, 64)` and `(400, 100)` lost their references to cache collisions (FOLD_WINDOW §2.6). The numbers were logged before the
  overwrite, and the statement cites them as evidence (iii) without that caveat.
* At `(1100, 100)` the statement quotes the class value `−378 ± 243`. The other class gives `+586 ± 179`. Both are inside the tolerance;
  quote the range.

**m7. "Flat [P]" is over-graded.**
* The crystal identity `EE ≡ q²` is exact on CLS-D's physical polarisation plane.
* In the engine's full null-circle contraction, `V_{2,04}` is *not* ∝ `H₂₀`, and `EE/q² → 1` only as `m → ∞` at fixed `J` [NE] (e.g.
  `J = 18`, `m ≈ 1600`: `0.89 − 0.19i`).
* EE_FUNCTIONAL §7 grades it correctly: [P given H0 §17's crystal + polarisation map]. Use that grade.
* Consequence for the evidence (not the bound): the `O(J²/m)` real EE−`q²` difference lives in the rows. Where only the hybrid/band
  model is used, it is part of `H_R`/`H_LJ`.

**m8. Pre-fix numbers quoted.** None changes a verdict:
* 300, `J ≤ 22` (S4);
* `J = 24` at 1.1Λ "`+1.06 ± 0.14`" (pre-curvature-fix; post-fix `+1.013`);
* `J = 28` at 1.1Λ "`+1.45 ± 0.58`" (post-fix central `+1.570`; repeat `+0.69`, S5).

The jackknife errors were never re-run after the curvature and Jacobi fixes.

**m9. The 300 n-sum `N₀` pathology (SMOOTH_BOUND hand-over item 3) is unexplained.**
* ROUTEA_LARGE_B excluded the records it flagged. But the "reliable" 300 records behind `H_R`'s `|S|/F_hyb ≤ 122/m²` were produced with
  the same default `N₀ = 200–320` configurations that diverge at `ν/Λ ≳ 0.85`.
* *Fix.* Re-run those records with the `n0t` configuration (`N₀ = 20/40/60`).

**m10. The `x₀₂² = 0` double pole of `E`.**
* The admissibility and convergence argument (EE_FOUNDATIONS §1.4, complete gauge fixing) is reasonable at [D].
* No numerical check of the light-cone-edge integrability with `EE` exists (the `q²` checks were single-pole).
* *Fix.* One numerical edge-scaling check in the CHMRS frame at `θ = 0`.

### Presentational

* **p1. Define "computed".** Use one protocol label per state: "direct small-`ν` + `K₂` (validated range)", "`K₂` outside validated
  range", "convention A", "model". Say "passes in the adopted protocol" rather than "passes".
* **p2. Use "in every variant", "main protocol" and "central values" consistently.**
  * For `J ≤ 22` at 300, "every variant" silently excludes the canonical gauge and the `q²` basis (`J = 18, 20`: `−0.63 / −0.25`).
  * For `J = 24` it excludes convention A (`−0.13`, post-fix) and the canonical gauge (`−37.5`).
  * The legitimate exclusions (a different functional; a diagnosed artefact) should be named.
  * I count at least nine independent analysis axes in the 300/1000 chains: estimator; extrapolation; gauge; basis; model/`K`;
    `XMAX`/`XMIN`; curvature image; Jacobi fix; despike; row set.
* **p3. Align the `H_LJ(300)` sentence with its own table.** The `J = 32` mass hole and the `J = 24, 28` edge masses (m2) are in neither.
* **p4. The header says "not committed".** The file is in commit `50d2260`.
* **p5. Literature (§6.4).**
  * "Stronger by `√(ln Δ_g)` than AKT" compares a conditional bound's functional form with a heuristic bound whose normalisation was not
    converted. Write "different parametric form; constants not comparable".
  * The novelty sentence rests on one literature search on 19 Nov; mark it as such.
  * The CHLPS comparison (×11 weaker, decomposition) is fair.

---

## 4. Answers to the six questions of the brief

1. **Logical chain with `EE` replacing `q²`.** Sound.
   * The light side is computed for `EE` itself: the pole by EE's (HN); `δk(0)` from EE's `X`.
   * `k(0)` comes from the crystal images, where `EE ≡ q²` on the physical plane suffices: at the bulk point only the physical
     polarisations survive, which the engine confirms [NE].
   * Two silent carry-overs remain: the H17 model (m4) and the `c₁` guess (m5). Both are small.
   * `EE ≠ q²` off the plane affects only the heavy *evidence* (finite `m`). The rows capture it there; elsewhere it is in `H_R`/`H_LJ`.
2. **Legitimacy of `EE` and of the gauge.**
   * Admissibility [D] is acceptable. Weights and span are [P]. The structure-independence of vanishing and convergence is [D].
   * The double pole: m10. The real zeros of `h₋+h₊` are harmless in the χ-normalised Abel form (`X` vanishes there, residual `10⁻¹⁵`).
   * Min-norm is a fixed, grid-converged functional. It was chosen *after* the canonical gauge failed (17:25 → 23:41 on 19 Nov), but it
     was the prior default, and the choice is between two candidates.
   * This does not affect validity. It does make `J = 18–28` (300) in-sample evidence (S3).
3. **Numerics.**
   * The binding states and the §9.1 odd-`J` results reproduce exactly.
   * Bug classes still open:
     * `K₂` used outside its convergence range (S1);
     * held extrapolation at the light-cone edge (S2);
     * despike gaps on off-diagonals and on the first node (m1);
     * pre-fix numbers quoted (S4, m8);
     * an unanalysed finished run (S5);
     * irreproducible fold numbers (m6).
   * Not found: non-finite drops (0), phase-strip violations (0), reference-key collisions in the 300 direct caches (none).
4. **Selection and protocol dependence.** See S3, S4 and p2. "Passes in every variant" is true only within a protocol that was revised
   after failures. Most revisions are justified on independent grounds, and the out-of-sample states pass.
5. **Scope of `H_M`, `H_OO`, `H_R`, `H_LJ`.**
   * `H_M` must explicitly include `H_edge` (S2) and the `K₂`-range exceptions (S1).
   * `H_OO` is a full sector (S7).
   * `H_LJ` is correctly graded [I] but undersold in size (S6).
   * Holes in the coverage map:
     * twist-edge masses at 300 (m2);
     * mass-grid resolution where `dF ~ O(1)` (m2);
     * `x > 0.8` for all `J` (S2);
     * `J = 26, 30` (300) and `26, 30, 34, 38, 42–46` (1000) without rows (listed).
   * Helicity sectors are complete apart from OO. Odd `J` is purely parity-odd, `J = 1, 3` are empty, and `J = 0, 2` are covered.
6. **Over- and under-claiming.**
   * Over: "every computed state passes" (S1, S2, S4, S7); "decidable in principle" (S6); the AKT comparison (p5); "Flat [P]" (m7);
     "`|c₁| ≤ 1`" (m5).
   * Under: nothing material. The statement is candid about 100 and about 300 being weaker.

---

## 5. Required changes to `RESULT_D_STATEMENT.md` (none made here)

1. §0.3 and §6.2: say "EE block"; regrade `J = 28, 32` (300) to [M]; add `H_edge`; add the S6 clause on the size of `H_LJ`.
2. §2: correct "difference `≤ 3·10⁻³ N_s` at 300" (S4); state the gauge chronology (S3).
3. §4: add the rows `H_edge` and "`K₂` beyond `J²/m ≈ 2`". Change the "Re" row (m3) and "Flat" to [P given H0 §17] + [NE] (m7). Note under
   H17 that it was carried over (m4).
4. §5.2: give post-fix numbers for `J ≤ 22` (convention A, `J = 22` `+0.123`, `t* = 1.1`) and the `J = 28` repeat (`+0.69`). Add the
   edge-mass holes.
5. §6.4: soften the AKT comparison and the novelty claim.

## 6. What Result D can honestly claim

> **Result D (conditional).** Assume A1, A2, V [D, K1 with KKSDZ's loophole], Adm [D], H7 [D], L0 [D], H17 [model carried over from
> `q²`] and `c₁` as a free parameter. Then per-state positivity of the fixed functional `Ω^EE` (min-norm gauge, NF02 weights)
> implies `√(t₄²+α_T²) ≤ 3.13·10⁻⁶` (`Δ_g = 300`) and `≤ 2.89·10⁻⁸` (`Δ_g = 1000`) [D; constants DV].
> * **Positivity at `Δ_g = 1000`** is computed [NE + M] on the EE block for the near-gap states:
>   * even `J ≤ 40` with rows (`J ≤ 24`, `28`, `32`, `36`, `40`);
>   * `J = 48` at the `K₂` convergence limit;
>   * odd `J ≤ 31`;
>
>   with the light-cone-edge window extrapolated. It is computed [M] elsewhere. The OO block is at leading order only. The large-spin
>   near-gap region (`J ≳ 50`) is an open hypothesis (`H_LJ`) with a tolerance of `1–30 c₂`.
> * **Positivity at `Δ_g = 300`** is computed for even `J ≤ 24` and odd `J ≤ 15`. Only `J = 24` uses the validated small-`ν` treatment;
>   `J ≤ 22` uses a superseded one, with a 10 % margin at `J = 22`. Everything at `J ≥ 26` is model-level or open.
> * At `Δ_g = 100` the bound is not established.
>
> None of the conditions is a spectral assumption on the CFT beyond A1. All of them are statements about explicit block integrals, and
> the large-spin ones are out of reach of present numerical methods.

This is a respectable, clearly conditional result. The 1000 claim rests on substantially more than the 300 claim. It should not be
described as "every computed state passes" without the protocol and EE-block qualifiers.
