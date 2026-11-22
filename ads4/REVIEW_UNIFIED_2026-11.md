# Hostile review before unification: Result D (`t₄, α_T`) + Result Q (`g₄`) → "Einstein gravity through quartic order"

*22 Nov 2026, referee agent. Unrefereed, not committed. No existing file was edited. Scope: the reading list of the brief (Result D:
`RESULT_D_STATEMENT`, `RESULT_D` update bullets, `REVIEW_EE_2026-11`, `HYPOTHESES_CLOSE(2)`, `PARITY_ODD_CLOSE`, `BELOW_FOLD(2,3,_SPINJ)`,
`HLJ_ANALYTIC`, `FOLD_WINDOW_LARGEJ`, `SMOOTH_BOUND`, and `HWIN_THIMBLE` where `RESULT_D.md` cites it; Result Q: `RESULT_Q_STATEMENT`,
`REVIEW_RESULTQ_2026-11`, `RESULT_Q_REVISION`, `RESULT_Q_V4PRIME/V4DPRIME/V4FINAL`, `RESULT_Q_STAGE0/12/3–7`, `RESULT_Q_FWD_CLOSE`,
`ABJM_THEOREM`). New code: `ads4/code/review/rv2_bj_entry.py`, `rv2_q_c1t4.py`, `rv2_q_largeJ_tol.py`; outputs in
`ads4/code/review/cache/rv2_*`. Everything ran with ≤ 1 process, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`, and read the `spinfunc/` and
`resultq/` caches read-only. The `cl_*` jobs were not touched.*

Severity: **fatal** (the stated implication or headline number is wrong) / **major** (a grade, a coverage claim or a hypothesis
list must change before unification) / **minor**.

---

## 0. Verdicts

| result | verdict | one line |
|---|---|---|
| **D** (`√(t₄²+α_T²) ≤ 3.13·10⁻⁶ / 2.89·10⁻⁸`) | **sound with corrections** | The implication "per-state positivity ⇒ bound" stands at [D] and the constants reproduce. But the authoritative statement is stale: six later documents are not in it, and several corrections from the prior review were adopted in the box but not in the body. The new below-fold engine evidence at `(336, 28/32)` fails an entry-resolved out-of-sample test (issue 2), so that window is [M], not [NE + M]. `H_LJ` is over-graded in `RESULT_D.md`. `J = 48` at 1000 still rests on `K₂` beyond its criterion. |
| **Q** (`0 ≤ g̃₄Λ⁶/8πG ≤ 17.24 log(2Λ/3) − 19.2`, at `ĝ₃ = 0`) | **sound with corrections** (upper and lower bound at `ĝ₃ = 0`); **flawed as stated** for the `|ĝ₃|²` sector and the implied cubic closure | The light-side circularity is resolved: the pointwise Lemma A [D] together with D_LR. Two items need to be named, though. One is an unnamed EFT-continuation hypothesis (`𝒢_rest` contains heavy-exchange double-trace parts that Lemma A does not cover). The other is that `H_K` carries the same grade as D's I14(b), so the two must be harmonised. Every `s_d`, binding states included, rests on an assumed real-section ≡ crystal equivalence for which no independent check exists at the relevant `ν`. The slope 0.475–0.49 (and the 120–159 closure) exceeds CHLPS's flat-space optimum 0.26 by ×1.9 with a `Δ_g`-independent correction. That is inconsistent with the flat limit unless the `σ_h U` term is double counted. |
| **Unified theorem** | **possible** with one hypothesis list (§4) and the wording of §5 | Gap, `Λ`, the `G` normalisation, the twist convention and the `k = 2` object (`Ĉ′_EE`, min-norm) are consistent. Three things must be explicit: the upper bound on `g₄` for `ĝ₃ ≠ 0` needs `κ ≥ 0` (issue 7); Q at `Δ_g = 300` inherits `H_LJ(300)` through `c₁t₄` unless `c₁t₄` is carried; and the cubic closure is not part of the theorem. |

No fatal error was found in either implication.

---

## 1. Independent checks run for this review

| check | script / output | result |
|---|---|---|
| **Entry-resolved version of BELOW_FOLD_SPINJ's calibration test.** Per entry of the 2×2 conserved-basis matrix: the state/ref double ratio direct/engine, normalised by the fitted `ρ`, with CV, jackknife | `code/review/rv2_bj_entry.py 28 2 5 1` and `… 28 3 4 2` → `cache/rv2_bj_entry_J28_d{2,3}.txt` | Scalar column reproduces `bj_calib_J28_d2/d3.txt` exactly. **Entries disagree out of sample** (table in issue 2): a¹, ν = 16.66: `−379 / +649 / +918 / +124 ± 195`; ν = 20.83: `−390 / +256 / +431 / −73 ± 27` (9–16σ). a², ν = 13.33: `+1251 … +133 ± 288`; 16.66: `−1641 … +595 ± 193`. |
| `bj_calib` files (existing) read | `spinfunc/cache/bj_calib_J28_{d2,d3,N4_K1_d2}.txt` | ν = 20.83 scalar residual: `+1.5 ± 27` (a¹, N = 5), **`+1437 ± 131` (a¹, N = 4)**, **`+1692 ± 27` (a², N = 4)**. ν = 24.996: `−1301`, `+14291`, `+15297`. |
| `c₁t₄` leverage in Q's δ = 0.2 functional, using D's `t₄` bound | `rv2_q_c1t4.py` | `2L_AdS = 72.041 / 92.814` (reproduced); `∫ψ₂ = 1.016`; shift `2c₁t₄Λ²∫ψ₂` ≤ **`0.567·c₁` (0.79 %) at 300, `0.059·c₁` (0.06 %) at 1000** |
| Does Q need an `H_LJ` analogue? Tolerance of Q's HYB action to a uniform or adversarial-sign graviton correction at near-gap `J = 8–100` | `rv2_q_largeJ_tol.py` → `cache/rv2_q_largeJ_tol.txt` | FF/N `≥ 0.29` (`+−` J = 8) and `≥ 0.97` for J ≥ 16. The adversarial tolerance is **`≥ 1.1·10⁴ c₂` (`+−` J = 8, 300), `≥ 1.6·10⁵` for J ≥ 16**, against D's `1–30 c₂`. The forward term `a(2J(J+1) − c₄)u³` dominates `N`. |
| Cubic closures compared | arithmetic | CHLPS optimum transplanted (D §6.4): `|ĝ₃|²Λ⁸ = 104.5 / 134.5`. Q with K = 0: `197 / 254`. Q with `s_d`: `120 / 154`. Result D: `1.07·10⁴ / 1.43·10⁴`. |
| ABJM dictionary | algebra | `M² = 4/α′`, `L²/α′ = π√(2λ)` ⇒ `ML = 2^{5/4}π^{1/2}λ^{1/4}` ✓. `g₄M⁶/8πG = (ζ(3)α′³/32)(64/α′³) = 2ζ(3)` ✓. `c_T = 32L²/(πG)` (free real scalar = 1) ⇒ `8πG = 256/c_T` ✓. |
| Gap and `Λ` conventions | `resultq/stage6_model.py` header vs `RESULT_D_STATEMENT` §1 | identical: `m² = (Δ−J−2)(Δ+J−1)`, `Λ² = (Δ_g−2)(Δ_g−1)`, `u_gap = (Δ_g−1)/(Δ_g+2J−1)` ⇒ `m_gap² = (Δ_g−2)(Δ_g−1+2J)` |
| Literature | arXiv abstract pages; one extended search | CLS-D 2311.04271 (v3, Feb 2025; JHEP04(2025)016) derives the dictionary and lists "remaining steps for uplifting flat space bounds to AdS" as open. CHMRS 2106.10274 uplifts two-sided flat bounds to AdS₄ for scalar contacts. No paper found that bounds R³/R⁴ from `⟨TTTT⟩` sum rules. |

---

## 2. Issues

### Major

**1. `RESULT_D_STATEMENT.md` is not authoritative any more; the prior review's corrections were boxed, not carried through.**
* *Evidence.* The statement cites none of HYPOTHESES_CLOSE, HYPOTHESES_CLOSE2, PARITY_ODD_CLOSE, HLJ_ANALYTIC, HWIN_THIMBLE or BELOW_FOLD_SPINJ (grep: 0 hits). Its
  `H_LJ` (§5.1/§5.2 tables: uniform `1–30 c₂`, fold, `A_E`) is the pre-HLJ_ANALYTIC definition. HLJ_ANALYTIC §5.1 replaces it by
  `H_AG ∧ H_win ∧ H_edge-amp`. Body items the box says are superseded but that still stand:
  * §0.1 "`|c₁| ≤ 1`" (review m5);
  * §0.3/§6.2 "every computed state passes", with no EE-block qualifier (S7);
  * §2 "`≤ 3·10⁻³ N_s` at 300" (S4);
  * §4 has no `H_edge` row and no "`K₂` beyond criterion" row; the "Re" row still claims no `Δ_g ≥ 300` number depends on the strip (m3); "Flat [P]" (m7);
  * §5.2 quotes `J = 22` `+0.15`, `J = 24` `+1.06 ± 0.14`, `J = 28` `+1.45 ± 0.58` (S4, m8, S5);
  * §6.4 says "stronger by `√ln Δ_g` than AKT" (p5);
  * the S3 in-/out-of-sample labelling was not done.

  The numbers now conflict across documents. `J = 28` at 300 is given as `+1.45 ± 0.58`, `+0.69`, `+0.835`, `+0.794`, `+0.799`, `+0.787` and `+0.739`. `J = 32` as `+0.657`, `+0.783`, `+0.865 ± 0.12` and `+0.834`. The binding state at 300 as `+0.2052` and `+0.2038`.
* *Fix.* Rewrite §0–§6 from the update bullets before unification. Give one protocol label and one number per state (min-norm EE, NF02, the
  HYPOTHESES_CLOSE §5 protocol), and import HLJ_ANALYTIC's decomposition with HWIN_THIMBLE's refined `H_win`.

**2. BELOW_FOLD_SPINJ: the calibration with one fitted constant per state does not validate the engine, and the "first evidence of
real-section ≡ crystal beyond LO" does not survive an entry-resolved test.**
* *What is done* (`bj_chain.calibrate`). There is one complex least-squares projection `κ = Σ conj(I)V/Σ|I|²` per node and replica, and
  `ρ = κ_st/κ_ref` is averaged over the **three noisiest nodes** (ν = 8.33 / 10.0 / 11.66, errors ±2756 / ±496 / ±254 c₂). The "shape
  check" then compares only this **scalar projection** at the other nodes. The quoted "χ² ≈ 3.4/6" includes the three calibration nodes.
  The `+1.5 ± 27` agreement at 20.8 is the only precise out-of-sample point.
* *Entry-resolved test* (`rv2_bj_entry.py`; c₂, ± jackknife; entries (0,0) (0,1) (1,0) (1,1) of the conserved basis):

  | ν | a¹ (depth 2, N = 5; used for J = 32 and the a¹ rows) | a² (depth 3, N = 4; the headline `+0.787`) |
  |---|---|---|
  | 13.33 | `+252 / −248 / −367 / +11 ± 287` | `+1251 / +353 / +133 / +810 ± 288` |
  | 16.66 | `−379 / +649 / +918 / +124 ± 195` | `−1641 / +121 / +595 / −776 ± 193` |
  | 20.83 | `−390 / +256 / +431 / −73 ± 27` | `+974 / +2162 / +2488 / +1558 ± 27` |

  The scalar `+1.5 ± 27` at 20.8 is a cancellation between entries that disagree at 9–16σ. It also becomes `+1437 ± 131` at N = 4 and
  `+1692 ± 27` at a² (the authors attribute this to N-non-convergence at 20.8; if so, the node cannot be evidence either way). Out of sample, the
  helicity structure (the off-identity part of the 2×2 matrix, which is what separates `++` from `+−`) is wrong by `0.5–2·10³ c₂` at
  ν = 13–21, in both orders.
* *Second problem.* The test region (ν ≥ 8.3, above the fold) is disjoint from the use region (ν < 8.3). Below the fold the crystal
  integral is a two-thimble sum (BELOW_FOLD §0, NEARGAP100), and BELOW_FOLD §3.2 shows the real section is *not* equivalent to the crystal
  cycle for J > 0 at ν ≳ 2.6 in the sampled form. The analytic-radial continuation may cure this, but it is not tested where it is used.
* *What survives.* The failure scale at `(336, 28/32)` is an extra uniform `Γ_c ≈ −4.3·10³`. Entry-level engine errors of `10³`, applied
  adversarially, cost roughly `0.2 N_s` (sensitivity ≈ `1.7·10⁻⁴ N_s/c₂`). Against `λ ≈ +0.74 … +0.83`, the window is plausibly positive,
  but **at [M], margin ×2–4**, not [NE + M]. The C3H window below 8.3 is not computed by the engine at all ("kept or dropped"). At J = 32
  there are no validated C3H data below 8.3 (HC2 §3.3: both thimble seeds fail the scalar Im check, 3.4σ / 19σ). "Closes HC2 §4 item 1
  at J = 28, 32" is therefore an overclaim.
* *Fix.*
  * Report entry- or helicity-resolved residuals.
  * Calibrate on the precise nodes and test on the noisy ones, or better, fix the normalisation analytically (BF_SPINJ §5 item 4).
  * Regrade the `(336, 28/32)` window to [M] with an engine systematic of `±2·10³ c₂` per entry.
  * Retract "first evidence of real-section ≡ crystal beyond LO" (replace it by "consistent in the scalar projection; inconsistent entrywise at 3–16σ").

**3. "Real-section ≡ crystal-cycle for EE beyond LO" (call it `H_RS`) is not flagged where it matters most.**
* *Where it is used.*
  * D: the `(336, 28/32)` window (issue 2); PARITY_ODD_CLOSE's beyond-order-0 support for `Γ_OO = Γ_++` (flagged there, §2.3 l.121, but
    not in the `RESULT_D.md` bullet).
  * Q: **every** `s_d`, including the binding `++`/OO J = 0, 2 gap states; `H_OO^Q`; the lower slope; the `k = 2` engine images of
    `RESULT_Q_REVISION` §1.3.
* *Where it is flagged.*
  * Flagged: BELOW_FOLD2/3, PARITY_ODD_CLOSE, RESULT_Q_FWD_CLOSE, RESULT_Q_STATEMENT box and table, V4FINAL item 8.
  * **Not flagged:** RESULT_D_STATEMENT (no mention), and the `RESULT_D.md` bullets for HYPOTHESES_CLOSE2, PARITY_ODD_CLOSE and BELOW_FOLD_SPINJ. The
    last calls the window "now computed [NE + M]".
* *Evidence base.*
  * Scalar: all J [DV].
  * EE: a⁰ (zeros of `N̂₄`); `k = 2` against `K₂` at ν = 20–40.
  * EE `k = 4`, J = 0, Δ = 400: an unresolved **512 c₂** discrepancy against `K₂` (BELOW_FOLD2 §3.3).
  * J = 28 (issue 2: fails entrywise).
  * **No independent EE beyond-LO check at ν < 10, where ≈ 96 % of `χ_d`'s weight lies.**
* *Why it matters for Q.* The binding margin at δ = 0.2 is `s_fail/s_d = 1.78/0.51 = 3.5`. An engine error of ×3.5 in the derivative-piece
  correction at J = 0 would break the δ = 0.2 functional at the gap. Nothing independent excludes that.
* *Fix.* Make `H_RS` a named hypothesis of both results, graded **[I]** (scalar [DV], EE LO [DV]). The highest-value single computation for
  Q is an independent EE `k = 4` crystal value at J = 0, `ν ≈ 2–8`, Δ ≈ 300 (n3/thimble, scalar-validated).

**4. `H_LJ(1000)` is over-graded in `RESULT_D.md` ([D + NE + M]).**
* HWIN_THIMBLE shows that the sign-free `H_win` is **violated** at `(1100, 100)`: C2 `+−` ≈ `2154 ± 43` against `T_B = 1409`, confirmed by fw's
  `−2123 ± 172`.
* The refined form (D_h class + remainder ≤ `T_B`) depends on HLJ §5.2, which is [D-sketch]. Its class tolerance is computed only at
  `Δ_g = 1000`.
* ν < 10 is unmeasurable at every J, and for J ≥ 200 (1000) / J ≥ 64 (300) the whole Bessel window lies below ν ≈ 10.
* `H_edge-amp` is [M] (margin ×3 at J = 48; ×1.6 at J = 34, 300).
* By the project's weakest-link rule (P4): **`H_LJ(1000)` = [D-sketch + M]** (data [NE] at J = 64, 100 on ν ≥ 13). **`H_LJ(300)` =
  [D-sketch / M] in the bulk, [I] in the window** for J ≥ 34 (no graviton window data). BELOW_FOLD_SPINJ adds nothing at high J (it says so).

**5. `J = 48` at `Δ_g = 1000` is still listed [NE + M] although it rests on `K₂` at `J²/m = 2.2`.**
* HYPOTHESES_CLOSE §5 showed that `K₂` at `J²/m = 2.35` (`(336, 28)`) is wrong by ×2–5 in amplitude and in phase.
* `(1048, 48)` has order ratio 0.70 (review S1), beside `(336, 28)`'s 0.79.
* No direct sub-row exists (HC2 §4 item 4 lists it as open).
* *Fix.* Regrade it to [M] and move J = 48 into `H_LJ(1000)`, or compute the `x = 0.125–0.8` direct sub-row (16 nodes, as for `(336, 28)`).

**6. The swap hypothesis is graded differently in D and Q although it is the same mathematical statement.**
* V4FINAL §2.4 shows that the de-shifted total kernel is not a position-space object: its weight grows like `e^{κν}`, `κ → π/2⁺`, and the
  raw per-ν heavy sum needs `J₀ < ½`. The heavy-sum ↔ (corner, ν) interchange is therefore a "Proposition K for the heavy-normalised kernel"
  (`H_K`, [D-sketch]), "the class of Result D's I14(b)". It says explicitly that D's `Ĉ′_EE` has the same `e^{π|ν|}` caveat.
* RESULT_D_STATEMENT §4 grades `V` (b) **[P given K1, K1′, K2, K3] → [D]**, and `K2`'s `η → ∞` decay as [P given K1′] via Prop. K. The two cannot both be right.
* *Fix.* One hypothesis `H_K` for both results, graded **[D-sketch]** (it becomes [D given `H_ip`]). D's `V` drops to [D-sketch].

**7. Q's upper bound is stated "at `ĝ₃ = 0`"; in a unified theorem `ĝ₃ ≠ 0` and the sign of `κ` becomes load-bearing.**
* The upper functional gives `g₄Λ⁶/8πG ≤ 2L_AdS − κ|ĝ₃|²Λ⁸`, with `κ = +0.110`. RESULT_Q_STATEMENT's own table grades `κ` [I], its box
  [D-sketch + M].
* With `κ > 0`, `ĝ₃ ≠ 0` only tightens the bound. If `κ` could be negative, Result D's `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴` would allow a shift of
  `|κ|·10⁴ ~ 10³`, i.e. 15× `2L_AdS`.
* "`ĝ₃ = 0`" is not a CFT-checkable condition. The CFT version is `t₄ = α_T = 0`.
* *Fix.* State the unified quartic bound either (a) for CFTs with `t₄ = α_T = 0`, or (b) for general `ĝ₃` with **`κ ≥ 0`** as an explicit
  hypothesis ([D-sketch + M]: the flat part of `κ` is exact, the `W³W³` AdS slope term is [D-sketch]).

**8. The implied cubic closure `|ĝ₃|²Λ⁸ ≤ 120–198 / 154–257` is overclaimed, and partly inconsistent with the flat limit.**
1. **Internal contradiction.** The RESULT_Q_STATEMENT box says "cubic terms restored [D-sketch + M]" and gives the closure. The headline and the
   regrading table, which come later in the same file, say the `|ĝ₃|²` terms are [I].
2. **Flat-limit inconsistency of the 0.475–0.49 slope.**
   * The measured `σ_h U` correction is `Δ_g`-independent (0.509 at 300 and 1000). The lower slope with it, `1/(4s_d) = 0.49`, therefore
     persists as `Δ_g → ∞`.
   * That is **1.9× CHLPS's flat-space numerical optimum `g₄/8πG ≥ 0.26|ĝ₃|²M²`** (Stage 0 reproduces the flat ¼ from the same forward rules).
   * An AdS bound whose flat limit beats the flat optimum is not credible if the AdS sum rules reduce to the flat ones. So either `s_d` is
     double counted against HYB's flat helicity-odd constant (review reading (a): slope 0.255, closure 197 / 254), or (b) holds and the
     light side owes a matching boundary layer that nobody has computed.
   * Either way the 120–159 numbers should go.
   * This is also the most concrete evidence on the open `σ_h U` question: it favours (a).
3. Even at K = 0, the closure needs Q's whole list, including the EFT/contact-tower hypothesis that D avoids by construction. It also needs
   the `W³W³` slope, which rests on the scalar-Mellin proxy (spinning [D-sketch]), and `H_RS`/`H_OO^Q` at the lower functional's binding
   state (`++`/OO J = 0).
4. Consistent with D (54–90× stronger, more hypotheses) and with CHLPS (197 / 254 vs the transplanted optimum 104 / 135, ×1.9). Nothing is contradicted there.
* *Fix.* At most a remark: "conditional corollary [D-sketch + M]: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257`". It should not be part of the
  theorem.

**9. Q's light side: Einstein is zero, but "zero + contacts" contains an unnamed hypothesis, and the highest-leverage lemma has no O(G) check.**
* *Not circular.* The pointwise Lemma A (V4FINAL §2.2) is a sound [D] argument: microcausality at fixed angles, CLS-D endpoint convergence,
  Lemmas B and C. D_LR is a faithful reading of CLS-D's domain. The ν-inside/outside circularity is genuinely gone.
* *However:*
  * **No O(G) check.** V4FINAL item 4 admits that no spinning O(G) Einstein evaluation was made. The only checks are MFT (exact) and the
    scalar-Mellin proxy. The Einstein term has leverage `10¹⁰–10²²` on `g₄`. Grade: **[D], unchecked**. The prior review's top request
    (one swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` at one ν) is still open.
  * **Heavy-exchange double traces.** `𝒢_rest` contains "the s-channel double-trace parts of heavy exchange diagrams" (REVISION §2.2). These
    have `J₀ = J_heavy ≥ 4`, so Lemma A does not cover them. Equating `Ψ̂[𝒢_rest]` with the sum of contact values from (F)+(P) requires
    the J₀ ≥ 4 part of `𝒢_rest` to be an EFT whose light-ray integrals, which diverge term by term (V4FINAL table: "position integral
    diverges"), may be replaced by the term-wise continued values, with the contact tower convergent. This is CHMRS's (5.3) plus their
    arc assumption.
* *Fix.* Name it **`H_EFT`** ([A]/[I]). It underlies both `g₄ ≤ …` and `g₄ ≥ 0`.

### Minor

**10. The Q statement's own lists are stale.**
* Hypotheses 3 ("V₄ incl. the [TT]₀ tail") and 9 (`H_OO^Q` assumed) predate V4FINAL and FWD_CLOSE.
* "Next steps" 1–3 are done.
* Rewrite the statement with V4FINAL §5's list plus the items of §4 here.

**11. `H_OO^Q` is not independent evidence.**
* The helicity lemma makes OO ≡ `++` to `10⁻⁹` at order 0 and to `5·10⁻⁶` in `s_d`, on the same engine.
* So `H_OO^Q` reduces to `H_RS` + `++`, i.e. [D given `H_RS`]. The same holds for D's `H_OO` (reduction [D] given an identity that is [NE] on
  `K₂` and on direct nodes at ν ≥ 8.3, untested below the fold).

**12. The `Δ_g = 300` J-drift of `s_d` is harmless, but it is a substitution.**
* The drift is `−0.07` at J = 7, 2.4× the truncation error, and of favourable sign for `+−`.
* The thresholds for J ≥ 7 are `≥ 8.9` (δ = 0.05) and `≥ 20` (δ = 0.2), against `|s_d| ≈ 0.5`. The binding J = 0, 2 have `J²/Δ ≤ 0.013`
  (300: 0.509 vs 1000: 0.510). The drift cannot matter at δ = 0.2.
* Quoting the J ≥ 7 channels "from `Δ_g = 1000`" for the 300 bound is still a substitution. Record it as part of `H_fwd^Q` [M], with
  the ×17 tolerance, and keep the cause [I].

**13. Q has no named `H_LJ` analogue. It needs only a mild one.**
* `rv2_q_largeJ_tol.py`: at near-gap J = 8–100 the δ = 0.2 functional's forward term `a(2J(J+1) − c₄)u³` dominates. A correction of either
  sign, adversarially shaped, must exceed `1.1·10⁴ c₂` (`+−` J = 8, 300) or `≥ 1.6·10⁵ c₂` (J ≥ 16) to matter, against D's measured raw
  corrections of `10²–10⁴`.
* State **`H_LJ^Q`** as [M] with margin ×10–10³. It rests on the HYB forward normal form at large J, where `s_d` is measured only to J = 11.
* This also locates D's `H_LJ` problem: it is the price of the forward-free functional.

**14. `c₁t₄`.**
* Using D's `t₄` bound, Q's bound moves by `0.57|c₁|` (0.79 %) at 300 and `0.059|c₁|` (0.06 %) at 1000 (`rv2_q_c1t4.py`).
* Importing D's bound makes Q inherit `H_LJ(Δ_g)`, `H_M`, `H_R`, `H_edge` and `H_RS`. Carrying `c₁t₄` explicitly avoids that.
* `c₁` must include the `k = 3` (`Ĉ₃`) `E × W³` term as well as `k = 2`. D's §3 mentions only `16πG c₁t₄` in `Ĉ′_EE`.

**15. ABJM.**
* The arithmetic is correct (§1), and the placement 3.3 % / 2.6 % reproduces.
* `c_T` does not enter the placement: `2ζ(3)` is λ-independent, and `Λ` depends on λ only. There is no flavour-current/`c_T` confusion in
  Result Q.
* The `c_T = (64/π²)∂²_m F` identity, which relies on the mass coupling to a flavour/R current in the stress-tensor multiplet, appears
  only in ABJM_THEOREM's independent `c_T` theorem. It is used correctly there, normalised to a free real scalar `c_T = 1`, consistent with
  `8πG = 256/c_T`.
* Wording: "ĝ₃ = 0 (maximal SUSY)" is wrong for k > 2 (N = 6). The reason is 10d IIA supersymmetry (no R³).
* It is outside A1 (light KK modes), it uses only the flat image (F), and it is log-dominated. It is a consistency check only. Also,
  `Λ ↔ Δ_g` is identified only at leading order in `1/(ML)`.

**16. Literature positioning.**
* CHMRS 2106.10274 is the methodological precedent: two-sided log(Δ_gap) AdS₄ bounds on contact couplings with dynamical gravity, from
  uplifted flat rules, for scalars. It must be cited as such, not as "a source of moments".
* CLS-D 2311.04271 (v3 2025; JHEP04(2025)016) builds the spinning sum rules and the crystal dictionary, and lists the uplift of flat
  bounds to AdS as the remaining step. Results D and Q are precisely that announced programme. Present them as "a first explicit-constant,
  conditional execution of CLS-D's uplift", with the CLS-D/CHMRS authors' priority on method stated.
* No published R³/R⁴ `⟨TTTT⟩` bound was found (abstract-level search plus one extended search). Novelty stays [I] at search level.
* D §6.4's AKT sentence still needs the softening of review p5.

**17. Lower bound `g₄ ≥ 0`: sound at the grade of its inputs.**
* It needs only the value piece: heavy image `P_{−½}(cosh β) > 0` at leading Regge; flat `2(|c⁺⁺|²+|c⁺⁻|²)/m⁸ ≥ 0`; the AdS correction `s_v`
  is `Λ⁻²`-suppressed.
* Its light side is `2g₄` by Lemma A + (F)+(P) + `H_EFT` + `H_K`.
* It is unaffected by `σ_h U` and by `s_d`. Grade **[D-sketch]**, limited by (F)+(P) spinning and `H_K`.
* Caveat: the `Λ⁻²` value-piece corrections (`b₀`, `s_v`, contact tower) mean "≥ 0" is really "≥ −O(Λ⁻²) × (EFT data)". This is harmless,
  but state it.

**18. Minor inconsistencies in D.**
* The binding state at 300 is `+0.2052` (statement) vs `+0.2038` (HC/HC2).
* HC2's J = 32 thresholds ("`−1.4·10³` on ν < 5", pooled unvalidated seeds) and BF_SPINJ's ("extra `−4.3·10³`") use different baselines.
  State both with the baseline.
* m2 (twist-edge holes) is partly closed by HC (edge masses added). Say so.

### Prior review (REVIEW_EE_2026-11): carried through?

| item | in box | in body / later docs |
|---|---|---|
| S1 (`K₂` beyond criterion) | yes | J = 28, 32 fixed by direct sub-rows (HC §5); **J = 48 (1000) not** (issue 5) |
| S2 `H_edge` | yes | computed on 8 rows (HC §2), [M] elsewhere; **no row in §4** |
| S3 chronology / in-sample labels | partly | labels never added |
| S4 300 `J ≤ 22` | yes | HC §4 recomputed J = 18–22; **§2 "≤ 3·10⁻³" and §5.2 "+0.15" remain** |
| S5 J = 28 repeat | yes | superseded by direct bulk; §5.2 still `+1.45 ± 0.58` |
| S6 size of `H_LJ` | yes | reworked by HLJ_ANALYTIC; then over-graded (issue 4) |
| S7 OO block | yes | HC/HC2 reduce OO to `++` (issue 11) |
| m3, m5, m7, p4, p5 | partly | §4 "Re" and "Flat [P]" rows and §6.4 unchanged; header says "not committed" |

---

## 3. Answers to the brief's questions

* **Is the calibration with one fitted constant per state in BELOW_FOLD_SPINJ legitimate?** As a modelling step ([M]), yes: the engine
  drops Δ-dependent constants, so a ν-independent normalisation is expected. As *validation*, no.
  * It is fitted on the three noisiest nodes and tested only through a scalar projection.
  * Entry by entry, out of sample, the a¹ engine misses by 3–16σ (`0.4–0.9·10³ c₂`) and the a² engine by 4–90σ (`10³–2.5·10³ c₂`).
  * The one precise agreement is a cancellation, and it is unstable under N and order.
  * It also absorbs a `c₀D_h`-shaped physical shift: "most of [the a² term's] effect on λ comes through the calibration constant".

  See issue 2.
* **Is the real-section ≡ crystal-cycle beyond LO assumption adequately flagged?**
  * Yes in BELOW_FOLD2/3, PARITY_ODD_CLOSE, FWD_CLOSE and the Q statement.
  * No in RESULT_D_STATEMENT (absent), in the `RESULT_D.md` bullets for HC2, PARITY_ODD_CLOSE and BF_SPINJ ("now computed [NE + M]"), and in
    REVISION §1.3 (engine `k = 2` images used without the caveat).
  * Nowhere is it said that the only large-J test is above the fold while the use is below it (issue 3).
* **Is the light side zero for tree-level gravity?** Yes at [D]: pointwise Lemma A + D_LR, with no circularity left. It is unchecked at O(G)
  in spinning form, and the leverage is enormous. The remaining non-trivial content is `H_K` (swap) and the unnamed `H_EFT` (issue 9).
* **Is `H_K` adequate?** It is honest as a label, but it is essentially "the de-shifted sum rule holds", with no criterion short of `H_ip`
  (`J₀ < ½` at the raw corner, not implied by KKSDZ). It is the main analytic hypothesis of **both** results. Grade [D-sketch] in both
  (issue 6).
* **Is `g₄ ≥ 0` sound?** Yes, at [D-sketch] (issue 17).
* **ABJM?** Correct arithmetic, uninformative, outside A1, and no `c_T`/flavour confusion in Q (issue 15).
* **`σ_h U` and the 300 J-drift?** `σ_h U` decides nothing at δ = 0.2 for the upper bound, but the flat-limit argument of issue 8.2 favours
  double counting and kills the 0.49 slope. The drift is harmless at δ = 0.2 (issue 12).
* **Is Q's implied `|ĝ₃|²Λ⁸ ≤ 120–257` consistent with D?** It is consistent (stronger, under more hypotheses), but overclaimed: the
  internal grade conflict, the flat-limit violation at 120–159, and dependence on `H_EFT`/`H_RS` (issue 8).

---

## 4. Regraded unified hypothesis table

Weakest-link grades. "D" / "Q" mark which result needs each item.

| # | hypothesis | D | Q | grade (this review) | was |
|---|---|---|---|---|---|
| A1 | twist gap `Δ − J ≥ Δ_g` for every single trace ≠ T, all J; unique T | ✓ | ✓ | [A] | same |
| A2 | tree level `O(1/c_T)`, `c_T ≫ Δ_g²` | ✓ | ✓ | [A] | same |
| K1 | KKSDZ Regge bound, with its loophole; K1′ uniformity | ✓ | ✓ | [Q] + [D-sketch] | same |
| **H_K** | heavy OPE sum ↔ (corner, ν) interchange for the heavy-normalised de-shifted kernels (`Ĉ′_EE`; `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`) | ✓ | ✓ | **[D-sketch]** ([D given `H_ip`]) | D: [D]; Q: [D-sketch] |
| Adm | EE admissible `j′ = 0` structure | ✓ | ✓ | [D] | same |
| H7 | `−Ĉ′_EE|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` | ✓ | ✓ | [D] | same |
| L0 | `[TT]` double zeros (k = 2, 3), contacts drop out | ✓ | — | [D] | same |
| Lemma A + D_LR | J₀ < 4 tree terms annihilated pointwise by `Ψ₄`, `Ψ₅` | — | ✓ | [D] (MFT [P]; **O(G) unchecked**) | [D] |
| **H_EFT** | the J₀ ≥ 4 part of `𝒢_rest` (contacts + heavy-exchange double traces) evaluates term-wise by (F)+(P); the contact tower does not conspire | — | ✓ | **[A]/[I]** (CHMRS (5.3)-type) | "local EFT" (unlisted for double traces) |
| (F)+(P) | spinning flat-limit and polynomial light side (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`, `g₅` cancellation) | — | ✓ | [D-sketch] | same |
| De-shift | continuum C4/C5 solutions exist, images for all β | (C2: [D]) | ✓ | [NE] to t = 30; C5 spinning [I] | same |
| H17 | `|ĝ₃|²` light coefficient `k(0)(1+δ)` | ✓ | (κ) | [D-sketch + M], carried over from `q²` | same |
| **κ ≥ 0** | sign of Q's `|ĝ₃|²` coefficient (needed for `ĝ₃ ≠ 0`) | — | ✓ | **[D-sketch + M]** | Q statement: [I]/[D-sketch + M] (conflict) |
| c₁ | `E × W³` light term (k = 2 and k = 3) | ✓ (±0.22 %/±0.017 % per unit) | ✓ (±0.79 %/±0.06 % per unit with D's `t₄`) | **free parameter, carried explicitly** | "|c₁| ≤ 1" in D §0 |
| **H_RS** | real-section ≡ crystal-cycle equivalence of the analytic-radial engine for EE beyond LO | (window 28/32; OO beyond order 0) | ✓ (every `s_d`) | **[I]** (scalar [DV]; EE LO [DV]; J = 28 entrywise fails at 3–16σ) | unlisted in D; [I] in Q |
| H_M | near-gap chain accurate on computed states (incl. `K₂` validity, `H_edge` on non-computed rows) | ✓ | — | [M], validated [NE] on rows with direct data; **J = 48 (1000): [M] only** | [M]/[NE + M] |
| H_edge | edge window `x ∈ (0.8, Λ/m]` | ✓ | (tolerance ≥ 760 c₂) | [NE + M] on 8 computed rows; [M] elsewhere | [I] → [NE + M] |
| H_OO | parity-odd even-J block | ✓ | ✓ (J = 0, 2) | [D given `Γ_OO = Γ_++`]: identity [NE] above the fold; reduces to EE `++`/`H_RS` | [NE + M] |
| H_R | beyond-leading-Regge at `m ≳ 1.3Λ`, large `b̂` | ✓ | (HYB) | [M] | same |
| **H_LJ(Δ_g)** | `H_AG ∧ H_win (refined) ∧ H_edge-amp` near the gap, large J; window ν < 8.3 at J = 28, 32 (300) | ✓ | — | **1000: [D-sketch + M]; 300: [D-sketch/M] bulk, [I] window J ≥ 34; (336, 28/32) window [M]** | `RESULT_D.md`: [D + NE + M] / [NE + M] |
| H_LJ^Q | large-J near-gap corrections below Q's tolerance (≥ 1.1·10⁴ c₂) | — | ✓ | [M], margin ×10–10³ (new) | unlisted |
| H_smear^Q | smeared graviton corrections at the binding states | — | ✓ | [NE + M] (`3·10⁻⁵` of B₂) | same |
| H_fwd^Q | unmeasured forward channels; J-drift at 300 | — | ✓ | [NE + M] given `H_RS`; drift [I], harmless at δ = 0.2 | same |
| σ_hU | double counting | — | δ = 0.05 and the cubic closure only | [I]; flat-limit argument favours double counting | open |
| HYB | accuracy elsewhere; β > 30 | (hybrid) | ✓ | [M] | [NE] |

---

## 5. Recommended wording of the unified conditional theorem

> **Theorem (conditional; "Einstein gravity through quartic order" at `Δ_g ∈ {300, 1000}`).** Let `𝒯` be a unitary CFT₃ with a unique
> stress tensor, `c_T` normalised so that a free real scalar has `c_T = 1`, and `8πG := 256/c_T`. Assume:
> * **A1:** every single-trace operator other than `T` has `Δ − J ≥ Δ_g` for all `J`;
> * **A2:** tree level (`O(1/c_T)`, `c_T ≫ Δ_g²`).
>
> Let `Λ := √((Δ_g − 2)(Δ_g − 1))`. Let `Ω^EE` be the min-norm-gauge NF02 EE functional (`RESULT_D_STATEMENT` §2) and `Ω^Q` the δ = 0.2
> functional (`functionals_q3.json`, `Q3-full-delta0.2`), both fixed as of 22 Nov 2026. Assume:
> * **analytic inputs:** KKSDZ's Regge bound with K1′, `H_K` [D-sketch], Adm, H7, L0 [D], Lemma A with D_LR [D], (F)+(P) [D-sketch],
>   `H_EFT` [A], de-shift existence [NE/I], H17 [D-sketch + M];
> * **per-state positivity** of `Ω^EE` and of `Ω^Q` on every twist-allowed heavy block. It is computed [NE + M] on the finite sets listed in
>   §4, and elsewhere assumed through `H_M`, `H_R`, `H_OO`, `H_edge`, `H_RS`, `H_LJ(Δ_g)`, `H_LJ^Q`, `H_fwd^Q`, `H_smear^Q`.
>
> Then, with `c₁` the (uncomputed, free) coefficient of the `E × W³` light term:
> 1. **(cubic)** `√(t₄² + α_T²) ≤ 240·√(16L_A/(k(0)(1+δ)))/Λ⁴ · (1 + O(c₁Λ⁻²))`, i.e. **`≤ 3.13·10⁻⁶` (`Δ_g = 300`) and
>    `≤ 2.89·10⁻⁸` (`Δ_g = 1000`)** in AdS units, `∓0.22 %` / `∓0.017 %` per unit `c₁`.
> 2. **(quartic)** If in addition `κ ≥ 0` [D-sketch + M] (or `t₄ = α_T = 0`), the CFT-defined R⁴ coefficient `g̃₄` (equal to the bulk `g₄`
>    given (F)+(P)) satisfies **`−O(Λ⁻²) ≤ g̃₄Λ⁶/(8πG) ≤ 17.24·log(2Λ/3) − 19.22` (`Δ_g = 300`) / `− 19.27` (`Δ_g = 1000`)** `+ ε_{c₁}`,
>    with `|ε_{c₁}| ≤ 0.57|c₁|` / `0.059|c₁|` when item 1 is used to bound `t₄`.
>
> The implication "positivity ⇒ (1), (2)" is [D] given the analytic inputs. Positivity itself is not established. It holds at the grades of
> §4, the weakest being `H_RS` and the large-spin window of `H_LJ(300)` ([I]), and `H_LJ(1000)`, `H_K`, (F)+(P) ([D-sketch]). Nothing is
> claimed at other `Δ_g`, at `Δ_g = 100`, at loop level, or for the ratio `|ĝ₃|²/g₄` (no cubic closure from (2)). ABJM satisfies (2)
> numerically (`2ζ(3)` at 3.3 % / 2.6 % of the bound), but it violates A1, so this is a consistency check only.

Before this can be released:
1. Rewrite `RESULT_D_STATEMENT` and `RESULT_Q_STATEMENT` (issues 1, 10).
2. Run one independent EE `k = 4` crystal check at J = 0, ν ≈ 2–8 (`H_RS`, issue 3).
3. Do one spinning O(G) Einstein evaluation of `Ψ₄` (issue 9).
4. Compute a direct sub-row for `(1048, 48)` (issue 5).
5. Report BF_SPINJ residuals entry by entry, and regrade (issue 2).

## Files

* `ads4/code/review/rv2_bj_entry.py` → `ads4/code/review/cache/rv2_bj_entry_J28_d2.txt`, `rv2_bj_entry_J28_d3.txt`
  (reproduce: `cd ads4/code/review; python3 rv2_bj_entry.py 28 2 5 1`; `… 28 3 4 2`; ≈ 1 min each).
* `ads4/code/review/rv2_q_c1t4.py` (stdout; seconds).
* `ads4/code/review/rv2_q_largeJ_tol.py` → `ads4/code/review/cache/rv2_q_largeJ_tol.txt` (≈ 2 min).
