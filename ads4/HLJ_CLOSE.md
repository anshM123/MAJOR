# Closing `H_LJ(300)`? The spin-J engine with an analytic normalisation, its convergence structure, and the `ν < 8.3` windows

*22 Nov 2026, H_LJ-close agent. Unrefereed, not committed. No existing document, code file or cache was edited or overwritten. New code:
`ads4/code/spinfunc/lj_rs_*.py`, `lj_rs_*.sh`. New caches: `spinfunc/cache/lj_rs_*`, plus the unchanged `hc_chain.py`'s own outputs
`cache/hc_chain_z4_lj_rs*` (new names). Read-only imports: `bj_ee`, `bj_eeq`, `bj_chain` (functions that write nothing), `bj_gamma`,
`fw_gamma`, `bf_core`, `z3_nodes`, `bg_*`, `hc_chain` (run unchanged, in-process wrappers). ≤ 2 worker processes, `OMP/OPENBLAS = 1`,
`ulimit -v 3000000` on every job (one exploratory scratch script without the limit was killed; it wrote nothing). Nothing in
`code/resultq` was touched.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model;
**[A]** assumption; **[I]** open. `c₂` = `m² ×` relative correction. "Engine" = the closed-form spin-J analytic-radial `1/Δ` engine
(BELOW_FOLD_SPINJ, `bj_*`), images `I_ac(ν)` in the 2×2 conserved basis; order `aᵏ` = relative geometric order `k` (depth `k+1`).
"Direct" = the independent n3/z2 crystal records (state 336 and reference 2·10⁴ on common random numbers, CV estimator).

---

> **Update (LARGESPIN_ANALYTIC.md, unreviewed): no large-spin hypothesis is replaced by a proof. H_R′ is upgraded [M] → [D-sketch + M]; H_LJ(1000) and H_bulk keep their grades. A possible convention problem (N2) affects H_LJ's leverage.**
> * **Common missing input:** the beyond-leading-Regge correction inside the Regge window (ν ≲ 3/β). It needs the ⟨TTTT⟩ heavy-block Regge formula at subleading order (the spinning analogue of CHMRS (3.17)–(3.20)), which is not derived.
> * **Lemma I [P]:** for an image ∫P_{−½+iν}(cosh t)σ(dt), the action is ∫Ω_LR(t)σ(dt). This gives explicit tolerance profiles for H_LJ (≥ 0.062 of the ring weight at t = β) and tolerated sizes of candidate correction shapes (ring amplitude a ≥ 80/m, ν² a ≥ 37/m², index-shifted a ≥ 10.8/m, …). The remaining model input R_win (non-ring window correction): exact scalar ≈ 10 c₂ (×38); graviton J = 100 ≈ 550 c₂ (leverage 1.1–1.5); no data at J ≥ 150.
> * **H_R′ [D-sketch + M] for b̂ ≥ b̂_A(m).** Lemma C′ [P + DV], a sign-free curvature bound using the decay of Φ_k, makes the hybrid floor positive for b̂ ≤ 119 at every m/Λ from 1.3 to 20 (the old Lemma C reached b̂ < 10). With ROUTEA, the leading-Regge hybrid is rigorous for all b̂ ≤ 2Λ when m ≥ 1.5Λ. The AG saddle test passes on 5 states. 5 ≤ b̂ < b̂_A stays [NE + M].
> * **H_bulk [M], restated:** the bulk correction may reach ≥ 3.9 % of the flat envelope (9.8 % at 1.3Λ) before any state fails. There is no analytic route.
> * **N2 (open, may weaken H_LJ):** Ω^EE's k = 2 part acts on raw images with the de-shift weight X, not ψ₂, and |X−ψ₂|/ψ₂ = 0.3–0.9 at ν ≤ 10. HLJ's T_B used ψ₂ inside the window. With |X|, the leverage at (1100, 100) drops to ≈ 1.1, and the sign-free ring-amplitude tolerance (3.2) is below the measured excess (2.7–5.9). Which convention the fw/pw Γ data use must be checked. **Until then quote the H_LJ(1000) leverage as 1.1–4.8 (not 2.7–4.8).**
> * **Corrections:** ROUTEA's certified window [91.2, 2Λ] holds only for m ≥ 1.5Λ (UT3 §5.2 and PE §2.4 quote it for m ≥ 1.3Λ). HLJ §5.2's [D-sketch] (subleading terms are P_{l+n}(cosh β) of the same β) does not by itself make window corrections harmless (N3). PE's "no sign-free argument can work" is true for δ but false for the curvature.

> **Update (N2_SUBREGGE.md, unreviewed): N2's premise is right (X, not ψ₂, multiplies the window data), but the claimed leverage collapse does not happen. A dir–cv estimator systematic at (1100, 100) is now the dominant uncertainty.**
> * **What the rows are [DV].** Window rows Γ are raw images / (R₋+R₊), and X multiplies them; Re ΣX𝒜[P_β] = Σψ₂P_β holds to 2·10⁻¹². Γ has poles at the real zeros of R₋+R₊ (ν = 1.171, 6.643), which explains HLJ_CLOSE's "local Γ singular at every J". H_win should be stated on the regular aΓ_c with weight |Y|.
> * **T_B with the exact weight:** (1100, 100): 1409 → 1439 (|X|) / 1152 (|Y|). Over all sample states the change is −28 % to +3 %. LARGESPIN's "1.1" and "3.2" were artefacts (ψ₂ + |X−ψ₂| overcounts; a different floor and structure were switched in; signed ring tolerance 31 vs measured 4.5–6.6). The D_h excess is an a_h normalisation, already in de-shift form (cross-checked against k = 3, which shares a_h). The published Γ − κ contains a stripped-κ artefact of ±120–460 c₂.
> * **Corrected H_LJ(1000) leverage at (1100, 100)** (cv estimator, raw-form refit): **3.0–4.2 central; ≥ 4.1 at 2σ with the full ring class.** (364, 64): 20–36. **With the dir estimator: 0.8–1.4.** This dir–cv systematic is unresolved and dominant.
> * **New explicit hypothesis R1 [I] for EE beyond leading Regge:** ring-type window corrections enter in de-shift form. It holds exactly for the a_h mechanism and for the scalar. In stripped form they would cost 27–228× the ring action at J = 150–344.
> * **Subleading Regge:** CHMRS (3.20) is checked against the exact scalar skeleton (9 states). All its terms sit on the state's own ring, with cost 6.5·10⁻⁴ (J = 100) to ≤ 6·10⁻⁵ (J ≥ 150) of the floor, so margins are ×17 to ×10⁴. The graviton window data fit the same ring classes, with non-ring remainder consistent with zero. The **spinning formula is not derived [I]** (obstacle: the EE = Σ_b identity decomposes the kernel, not the spinning block; it needs ⟨TTO⟩ structure ratios). At each finite 1/m² order the corrections are ring-located or index-shifted at the state's β [D-sketch].
> * **Corrections:** the "leverage 1.1" in the LARGESPIN box above is superseded. Quote "3.0–4.2 (cv), ≥ 4.1 at 2σ; 0.8–1.4 with dir". HWIN_THIMBLE: refit the D_h class in de-shift form, after which the (364, 64) remainder is 7–13 c₂.

## 0. Summary

**Bottom line.** `H_LJ(300)` is **not closed**: its Regge-window clause for `J ≥ 34` stays **[I]**. What changed: the two computed
large-J states (336, 28) and (336, 32) move from [M] to **[NE + M]**, with no fitted parameter, and the obstruction for `J ≥ 34` is now
identified precisely. `H_LJ(1000)` stays [D-sketch + M].

| item | result | grade |
|---|---|---|
| Diagnosis (§1) | The large `a¹` term (`≈ 6·10³ c₂` at Δ ≈ 300, `2·10⁴` at 1100) is a **J-independent** normalisation of the EE radial moments (coefficient `≈ 16–21/Δ`, present at J = 0, geometric with ratio `≈ 1/25` there through `a⁴`), common to all entries. It is not the `J²/Δ` obstruction. The `J`-enhanced part is the **helicity structure** (entry differences): `≤ 25–60 c₂` per order in the window at J = 28, 32, but its series in `1/Δ` at fixed `J/Δ` is **asymptotic** (ratio `≈ 1–2` at ν ≳ 5 for `J²/Δ = 2.3`), growing like `ν⁶–ν⁸` and with `J²/Δ`. | [DV + NE] |
| Resummation (§2.3) | **Not achieved** for the structure series [I]. Replaced by optimal truncation with error estimates validated against data. The common part needs no resummation. | [I] / [NE] |
| Analytic normalisation (§2.1) | The engine/direct constant is fixed point-wise on the real section, `ρ_an = (a_st/a_ref)⁶ c_ref/c_st`: constant to `6·10⁻⁵` (median), ν-independent to `10⁻¹⁰`, truncation-free to `≤ 2 c₂` (J = 28). It matches the fitted constants within their errors (`+203 ± 295`, `−137`, `−120 ± 540 c₂`). For J ≥ 64 it needs fw_frame's Fourier cure. | [D + DV] |
| Prop. RS at J = 28, 32 (§2.2) | B1–B4 hold on the full J = 28 and J = 32 engines (identical exponents to J = 0; no negative powers in 5.4·10⁶ monomials; parity violators `≤ 9·10⁻¹¹`). So HRS_CLOSE's cycle argument covers the large-J window below the fold too. | [DV]; H_RS there [D given RD] |
| Validation, entry-resolved, no fit (§3) | Common part: the `a¹` engine matches the direct nodes within 1–1.5σ on ν = 10–16.7 (J = 28) and 11.65 (J = 32). Structure part (errors `±1–14 c₂`): the "3–16σ entry misses" are the engine's own truncation; residuals track the omitted terms with one factor (`2.2–3.1`, all components, ν = 10–20.8), reproduced out of sample at J = 32. **Unexplained:** a negative structure offset of `20–120 c₂` (J = 28) and `30–330` (J = 32) at ν = 5–10. | [NE] |
| (336, 28) window (§4) | **`λ/N_s = +0.759`** (C3H kept) / **+0.765** (dropped), gap state; `0.751 … 0.772` under `±1`, `±3.3 ×` the last term; `+0.80` with the measured structure offset; `n = 1`. Failure needs a uniform common `Γ_c ≈ −4.6·10³` (margin `≥ 15`) or a uniform **positive** structure error `≈ +215 c₂`, which is `3–6×` the engine's structure terms and opposite in sign to every measured discrepancy. | **[NE + M]** (was [M]) |
| (336, 32) window (§4) | **`λ/N_s = +0.837`** (`a²`; min at `m/Λ = 1.3`), `0.835 … 0.843` under the stresses, `n = 1`; failure at a common `−4.4·10³` or a structure `≈ +490`. | **[NE + M]** (was [M]) |
| High J (§5) | (1100, 100) and (364, 64) through fw_gamma with `ρ_an`: the `a¹` engine's window values, *after* removing the `D_h` class, are `6·10³–1.5·10⁴ c₂` (`4–11 × T_B = 1409`, 1100) and `6–7·10³` (`25 × T_B = 255`, 364), and its last term after the same removal is `20–55 × T_B`. The engine **cannot decide** `H_win` at `J²/Δ ≳ 9`, in either form. J = 34, 48, 150, 200 were not computed (no fw/`K₂` rows; and the structure asymptotics would forbid it). | [I] (negative) |
| Sign-free bound (target 3) | Not obtainable from the engine for J ≥ 64 (300) or J = 100 (1000): the truncation estimate exceeds `T_B`. | [I] |


---

## 1. What the obstruction actually is (diagnosis from the cached images) [DV + NE]

### 1.1 The large `a¹` term is a J-independent EE normalisation, not a `J²/Δ` effect
Per entry, `m²·P₁/P₀` (state, local ratio of the `a¹` to the `a⁰` image), ν = 2–10, away from the entries' common zeros (ν ≈ 1.07, 5.5–7, 13.3):

| state | `a¹` | `a²` | `a³` | `a⁴` | `a¹` coefficient `f₁ = Δ·(a¹ term)/m²` |
|---|---|---|---|---|---|
| J = 0, Δ = 300 (`bg_eeq_D300_J0_K4`) | +6230 ± 40 | +270 ± 40 | +9 | +0.3 | 20.8 |
| J = 28, Δ = 336 (`bj_eeq_D336_J28_d3`) | +4800 … +6500 | +160 … +305 | — | — | 15–19 |
| J = 100, Δ = 1100 (`bj_eeq_D1100_J100_d2_p`) | +1.7 … 2.1·10⁴ | — | — | — | 16–19 |
| reference Δ = 2·10⁴, J = 28 / 100 | 4.2·10⁵ (`m²` of 2·10⁴) | ≤ 0.6 % of it | — | — | 20.9 / 20.9 |

* The `a¹` term is (to 0.3 %) **the same for all four entries** and nearly ν-independent; its coefficient `f₁ ≈ 16–21` is independent of `J`
  and of `J/Δ` (it is present at J = 0). It is the radial-moment normalisation of the EE density (`∝ X₁^{9.25}`, `⟨X₁⟩ ≈ 10`), and in
  `c₂` units it is `≈ f₁ m`: `≈ 6·10³` at Δ ≈ 300, `≈ 2·10⁴` at Δ ≈ 1100. At J = 0 the series is geometric with ratio `≈ 1/25` through
  `a⁴`. The "`a¹ ≈ 2–4·10³` (J = 28), `8·10³` (1100, 100)" of BELOW_FOLD_SPINJ §2 are this term after subtracting its `ν < 1` mean, i.e.
  mostly the shift of the entries' zero at ν ≈ 1.07; they are not `J²/Δ`-enhanced terms.
* A common factor of this kind is a normalisation. With a fitted constant (BELOW_FOLD_SPINJ) it is absorbed differently at each order
  (the fitted ρ moves by `−137 → +203 c₂` from `a¹` to `a²`, §2.1), which is why the window `λ` moved by `+0.05` per order there.
  With the analytic normalisation of §2.1 it no longer moves (§4).

### 1.2 The genuinely J-dependent part: the helicity structure
Split the 2×2 state/reference log-ratio into a **common part** (entry mean) and a **structure part** (entry `(a,c)` minus entry `(0,0)`;
normalisation- and calibration-free; it is what separates `++` from `+−`). Per-order terms at J = 28 (`c₂`, `m²·Re log`):

| ν | 2 | 3 | 4 | 5 | 7 | 8.33 | 10 | 11.66 | 13.33 | 16.66 | 20.83 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| common `a¹` | +6305 | +5643 | +5331 | +4877 | +6220 | +5371 | +4676 | +3565 | (+464, zero) | +12099 | +7381 |
| common `a²` | +282 | +233 | +163 | +18 | +430 | +138 | −95 | −399 | (−1049) | +452 | −770 |
| structure `a¹` (01, 10, 11) | 6, 2, 3 | 7, 4, 3 | 9, 6, 3 | 15, 13, 7 | −28, −39, −15 | −8, −15, −5 | 19, 18, 9 | 99, 116, 49 | 425, 522, 212 | −1289, −1627, −645 | −1143, −1452, −570 |
| structure `a²` | −9, −13, −5 | −4, −7, −2 | 3, 2, 2 | 21, 24, 11 | −45, −59, −22 | −5, −8, −2 | 35, 42, 18 | 116, 145, 59 | 392, 493, 197 | −746, −950, −370 | −472, −602, −233 |

(ν = 6, the entries' zero: common `+318 / −1512`, structure `123–142 / 239–300`.)
* **Common part:** geometric, ratio `0.03–0.07` on ν ≤ 8.3 (outside the zero at 5.5–7); the ratio grows with ν (`≈ 0.1` at 11.7).
* **Structure part:** small in the window (`|a¹|, |a²| ≤ 25 c₂` on ν ≤ 5 and at 8.33; `≤ 60` at 7) but **not convergent**: `a²/a¹ ≈ 1.2–2`
  at ν = 5–11.7, `≈ 0.9` at 13.3, and it grows like `ν⁶–ν⁸`. At the reference (`J²/Δ = 0.04`) the same ratio is `≈ 0.06`. This is the
  `J`-enhanced non-uniformity; at J = 28 (`J²/Δ = 2.3`) the structure series is asymptotic, not convergent, at ν ≳ 5.
* At J = 32 (order `a¹` only): structure `a¹` = 7–43 `c₂` on ν = 2–5, 14–85 at 7–8 (zero at 6: local blow-up). At (1100, 100) (`J²/Δ = 9.1`):
  structure `a¹` = 15–150 (entries 01, 10) and 48–590 (entry 11) on ν = 2–10, 2.8·10³ at 11.7.


## 2. Method: analytic normalisation, bookkeeping, optimal truncation (and what was not achieved)

### 2.1 Analytic normalisation `ρ_an` (`lj_rs_norm.py`) [D + DV]
* The engine integrand (`bg_geom.assemble` of the `bj` pieces, with the `B`-factor, exactly what `bj_ee.images` integrates) is the
  `a`-expansion of the **same holomorphic integrand** as the direct one (w0 EE integrand of the z2/n3 records), written in the scaled
  variables of the real section `D` with every overall `a`-power dropped (log factors, `v′ − u′`, measure, `pref`) and its own branch
  constants. So point-wise on `D` (where both are defined and the series converges absolutely at small scaled radii)
  `f_eng(p; ν)/f_dir(p; ν) = c(Δ, J)` is `p`-independent up to the engine's truncation, with `f_dir` = `bf_core.grav_parts` (n1_pw
  point-wise w0 integrand, sheet `n₀ = +1`, `ψ = α − φ ∈ (0, π)`).
* By HRS_CLOSE Prop. RS (one populated cohomology class, `c_Γ(ν)` common to state and reference) the crystal integrals then obey
  `V_st/V_ref = (a_st/a_ref)⁶ (c_ref/c_st) I_st/I_ref`, i.e. `ρ_an ≡ κ_st/κ_ref = (a_st/a_ref)⁶ c_ref/c_st` (`bf_core` drops `a⁶`).
  **No parameter is fitted.**
* Checks (10 points, 4 entries, 4 ν):

| state | radial scale ×0.3 / ×1 / ×3: max (median) `|r/c − 1|` | `ρ_an` (×0.3 / ×1 / ×3) | ν-dependence |
|---|---|---|---|
| (336, 28), depth 3 | 1.8·10⁻⁴ (1.6·10⁻⁵) / 6.0·10⁻⁴ (6.4·10⁻⁵) / 1.7·10⁻³ (2.6·10⁻⁴) | 8.798991 / 8.798973 / 8.798927 ·10³¹ | `≤ 6·10⁻¹¹` |
| (336, 32), depth 2 | 7.5·10⁻⁴ (4·10⁻⁵) / 3.6·10⁻³ (1.8·10⁻⁴) / 2·10⁻² (2·10⁻³) | 8.79916 / 8.80107 / 8.81918 ·10³¹ | `≤ 6·10⁻¹¹` |
| references (2·10⁴, 28 / 32) | `≤ 10⁻⁵` | — | — |

  The spread shrinks with the radial scale as a truncation effect must. The ×0.3 value is the truncation-free limit to `≲ 2·10⁻⁶`
  (J = 28) / `≲ 2·10⁻⁵` (J = 32, depth 2). Production: the ×1 value for the J = 28 runs (0.2 `c₂` from the ×0.3 one) and the J = 32
  `a¹` runs (**24 `c₂`** from it, i.e. `≈ 0.004` in `λ/N_s`); the ×0.3 value for every later run (J = 32 `a²`, structure stresses).
* **Against the fitted constant** (bj_chain.calibrate, ±295 / ±540 `c₂` jackknife): `m²(ρ_fit/ρ_an − 1)` = **+203 − 115i** (J = 28, `a²`),
  **−137 − 114i** (J = 28, `a¹`), **−120** (J = 32, `a¹`). The analytic constant lies inside the fit's error at both orders; the fit
  moves by 340 `c₂` between orders, the analytic one does not move (it is order-independent by construction).

### 2.2 Bookkeeping at J = 28 (`lj_rs_book.py`): Prop. RS applies [DV]
On the full J = 28 engine integrand (`a²`, depth 3, N = 4; 5.4·10⁶ monomials per entry, all four entries):
* **B1** `E₁` (the ν-part of the leading exponents) is identical to the J = 0 EE2 one: `max|ΔE₁| = 0`;
* **B2** `E₀ − E₀(J = 0)` = 0 exactly (same lattice point);
* **B3** no negative power of `s, sf, ca, sa` (0 of 5.4·10⁶ monomials);
* **B4** the parity rules hold up to round-off (largest violating `|c|/max|c|` = `9·10⁻¹¹`).

The same holds at **J = 32** (depth 3, 7.1–7.5·10⁶ monomials per entry: B1, B2 exact, B3 0 violations, B4 `≤ 3·10⁻¹¹`;
`cache/lj_rs_book_J32_D336_d3.txt`).
So the J-resummed contraction (`D̂`-series, closed-form `B_k`) does not leave HRS_CLOSE's local system: the J = 28 engine is, order by
order, `c_Γ(ν)⁻¹ ×` the crystal-cycle integral, with the same `c_Γ` as the scalar and the reference, **below the fold as well**. This
answers the review's "test region disjoint from use region" objection (RU issue 2, second problem) at the level of Prop. RS:
[D given RD]. Not re-checked at J = 28: B5 (non-resonance; inherited, the exponents are identical) and B6 (no boundary layer) for the
new J-dependent factors (`√(V²W²)⁻¹`, `D̂`).

### 2.3 Resummation: what was attempted and what was achieved
* **Common part**: needs no resummation. Its large `a¹` is the J-independent normalisation of §1.1; the series is geometric with ratio
  `≤ 0.07` in the window (`≈ 1/25` at J = 0 through `a⁴`). An exponentiated (cumulant) form changes the `a²` image by
  `≈ ½(P₁/P₀)²m² ≈ 5 c₂` at (336, 28): irrelevant.
* **Structure part**: the `J`-enhanced terms are here. At J = 28 their series is asymptotic (ratio `≈ 1–2` at ν ≳ 5); the brief's
  suggested resummations (exponentiating a phase drift, or expanding at fixed `J²/Δ`) were not carried out: the J-dependence of the
  geometric coefficients enters through the conserved coefficients and the closed-form kernel in a way that is not a single phase, and a
  fixed-`J²/Δ` reorganisation needs all geometric orders of the `J²`-leading terms, which the engine does not produce. **A resummed
  engine was not obtained [I].**
* **What is used instead: optimal truncation with an error estimate, validated parameter-free against the direct nodes (§3).**
  Production: `a²` for both parts (the structure terms are `≤ 25 c₂` in the window); error band: `±1` and `±3.3 ×` the last computed
  term (whole image), plus separate structure-only stresses at the size of the unexplained structure offsets of §3.


## 3. Validation against the direct nodes, entry-resolved, no fitted constant (`lj_rs_val.py`) [NE]

Per node and common-seed replica (CV estimator), `m²Re[(V_st/V_ref)_ac/(ρ_an I_st/I_ref)_ac − 1]`; jackknife over the 8 replicas.
Errors of the **structure** combinations (entry − entry (0,0)) are computed on the differences themselves: they are `±1–9 c₂`, 20–500×
smaller than the common errors, because the MC noise is almost entirely common to the four entries (the review quoted per-entry errors
only, RU issue 2). `cache/lj_rs_val_J28_d3_N4_K2.txt`, `cache/lj_rs_val_J28_d3_N4_K2_win.txt` (with the n4 nodes 5.0, 6.67), `cache/lj_rs_val_J32_d2_N5_K1.txt`, `cache/lj_rs_val_J32_d3_N4_K2.txt`.

### 3.1 (336, 28): residual direct − engine at each engine order
Common part (entry mean; `R_k` = residual with the engine truncated at `aᵏ`):

| ν | 8.33 | 10.00 | 11.66 | 13.33 | 16.66 | 20.83 |
|---|---|---|---|---|---|---|
| `R₁` (`a¹`) | −1512 ± 1466 | −373 ± 401 | **−75 ± 168** | **−208 ± 193** | **+229 ± 151** | (+1232 ± 18) |
| `R₂` (`a²`) | −1650 ± 1466 | −278 ± 401 | +324 ± 168 | +841 ± 193 | −223 ± 151 | (+2002 ± 18) |
| last term `a²` (common) | +138 | −95 | −399 | −1049 | +452 | −770 |

(20.83: the `D̂`-series is not converged at N = 4, `|N4 − N5| ≈ 10³`, BELOW_FOLD_SPINJ §2; the node tests nothing.) The engine with
**no free constant** reproduces the direct nodes on 10–16.7 within 1.5σ at its optimal truncation (`a¹`), and within the size of the
first omitted term at `a²` (`|R₂| ≤ 0.8 |a² term|` on 11.7–16.7). The common series is asymptotic above ν ≈ 11 (the `a²` term is
already larger than the agreement of `a¹`).

Structure part (entries (0,1), (1,0), (1,1) minus (0,0); `R₀`, `R₁`, `R₂`; first omitted terms `s₁`, `s₂` from §1.2):

| ν | `R₀` | `R₁` | `R₂` (± jk) | `s₁` | `s₂` |
|---|---|---|---|---|---|
| 8.33 | −45, −56, −23 | −37, −41, −18 | −32, −33, −16 (± 4, 5, 9) | −8, −15, −5 | −5, −8, −2 |
| 10.00 | −50, −63, −21 | −69, −81, −30 | −104, −123, −48 (± 3, 3, 6) | +19, +18, +9 | +35, +42, +18 |
| 11.66 | −60, −77, −23 | −159, −193, −72 | −275, −338, −131 (± 2, 2, 4) | +99, +116, +49 | +116, +145, +59 |
| 13.33 | −83, −105, −33 | −508, −627, −245 | −900, −1120, −442 (± 3, 4, 7) | +425, +522, +212 | +392, +493, +197 |
| 16.66 | −269, −337, −148 | +1020, +1290, +497 | +1766, +2240, +867 (± 4, 5, 8) | −1289, −1627, −645 | −746, −950, −370 |
| 20.83 | −425, −537, −217 | +718, +915, +353 | +1190, +1517, +586 (± 1, 1, 2) | −1143, −1452, −570 | −472, −602, −233 |

* The review's "3–16σ entry misses" are the structure part. They are **the engine's own truncation**: at every node ν ≥ 11.7 the
  residual of the order-`aᵏ` engine has the size and the sign of the omitted terms, the order-`a⁰` structure is the closest
  (`|R₀| < |s₁|` at 11.7–20.8), and `R₂/(a¹ − a²)` is the same number `2.2–2.5` for all three structure components and all of
  ν = 11.7–20.8 (a factor 20 in size, two sign changes): the residual tracks the engine's series term by term. That is what an
  asymptotic series at its optimal truncation looks like; it is not a cycle (H_RS) discrepancy, which would not follow the engine's
  own `a²` term.
* **Unexplained**: an offset of `−20 … −60 c₂` (J = 28) at ν = 8.33 and 10.0 that no engine order reproduces (`|R_k| > |s_k|` there),
  3–15σ on the precise structure errors; also at the n4 thimble node ν = 5.0 (`R₂ = −66, −77, −33 ± 5`, against `s₂ = 21, 24, 11`). Same sign and pattern at J = 32 (§3.2), larger: `−40 … −130` at 8.32, `−80 … −220` at 10.0.
  Candidates: a systematic of the below-fold z2 nodes (8.33, 10.0 are below ν* ≈ 11.3), or a beyond-all-orders piece of the structure
  series [I]. It is carried as a structure stress in §4.

### 3.2 (336, 32), engine order `a¹` (ρ_an, no fit)
| ν | 8.32 | 9.99 | 11.65 |
|---|---|---|---|
| common `R₁` | −2602 ± 1659 | −533 ± 598 | +14 ± 400 |
| structure `R₁` (± jk) | −94, −37, −40 (± 6, 8, 14) | −224, −83, −82 (± 5, 7, 12) | −731, −278, −276 (± 5, 7, 13) |
| structure `a¹` term `s₁` | −36, −14, −9 | +58, +24, +30 | +500, +197, +207 |

Common part consistent at every node (≤ 1.6σ). Structure: same pattern as J = 28, `|R₁| ≈ 1.5–4 |s₁|`.

**(336, 32) at order `a²`** (new depth-3 images, 4.9 h; `ρ_an` at radial scale ×0.3; `cache/lj_rs_val_J32_d3.log`):

| ν | 8.32 | 9.99 | 11.65 |
|---|---|---|---|
| common `R₂` / `R₁` | −2680 / −2573 ± 1658 | −336 / −511 ± 600 | +756 / +43 ± 402 |
| structure `R₂` (± jk) | −79, −30, −35 (± 6, 8, 14) | −332, −123, −125 (± 5, 7, 12) | −1237, −473, −472 (± 5, 8, 13) |
| `R₂ / (a¹ − a²)`, structure | −4.9, −4.3, −7 | 3.1, 3.1, 3.0 | 2.46, 2.46, 2.43 |

* **Out-of-sample check of §3.1's reading.** From J = 28 one predicts, for J = 32, that (i) the common part agrees best at `a¹` and
  `a²` overshoots above ν ≈ 11 — observed (`+43 ± 402` vs `+756`); (ii) the structure residuals follow the engine's last term with the
  same factor `≈ 2.3–3` on ν = 10–12 — observed `3.0–3.1` and `2.4–2.5` (J = 28: `2.7–3.0`, `2.2–2.4`); (iii) an unexplained negative
  structure offset at 8.3 — observed (`−79, −30, −35`). The engine's truncation behaviour is therefore a property of the series, not a
  fit to J = 28.


## 4. The `ν < 8.3` windows at `Δ_g = 300`: (336, 28) and (336, 32) through the unchanged chain (`lj_rs_chain.py`) [NE + M]

Protocol: exactly BELOW_FOLD_SPINJ §3.1 / `bj_chain.chain2 + chain_win` (two passes of the unchanged `hc_chain.py z4 300`; engine C2
records at every engine node `ν < 8.3`, the k = 2 NL density below 8.3 replaced by the engine's `Γ_NL`; everything above 8.3 is the
HYPOTHESES_CLOSE §5 protocol), **with `ρ := ρ_an` for every replica** (the per-replica reference constant still comes from the
reference records, which the reference engine matches to `10⁻⁷–10⁻⁵`). C3H below 8.3: `keep` (the n3/thimble nodes, as the main run)
or `drop`. Environments: `lj_rs_env28.sh`, `lj_rs_env32.sh` (identical to the `bj_chain` logs). `n = 1.0` in every run (the Δ_g = 300
minimum stays J = 5, `+0.2054`).

### 4.1 Engine values below the fold, (336, 28), `a²`, `ρ_an` (`Re Γ_c`, `c₂`; `cache/lj_rs_nl_J28_a2.pkl`)
| ν | 0.5 | 1.0 | 1.5 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|
| `++` | +808 | +854 | +689 | +705 | +661 | +586 | +472 | +204 | +892 | +478 |
| `+−` | −2020 | −5647 | +1944 | +401 | −216 | −562 | −1066 | −2656 | +3517 | +150 |
| `++`, `a¹` | +556 | +696 | +354 | +424 | +444 | +449 | +459 | +521 | +66 | +293 |

(ν ≈ 6–7 straddles the zero of `h₋ + h₊`, ν ≈ 1 the zero of `N̂`: the local `Γ_c` swings there, the NL density used by the chain is finite.)
BELOW_FOLD_SPINJ's fitted-constant `a²` values were `++ +903 … +1291`: the difference, `≈ 250–400`, is the fitted constant (+203 `c₂`
of ρ plus its jackknife).

### 4.2 `λ/N_s` per state (`hc_chain_z4_lj_rs*`, gap state is the minimum in every run)

| state | engine | C3H < 8.3 kept / dropped | at `m/Λ` |
|---|---|---|---|
| (336, 28) | **`a²`, ρ_an (production)** | **+0.7592 / +0.7652** | 1.091 (gap) |
| | `a¹`, ρ_an | +0.7568 / +0.7630 | gap |
| | `a² ± 1 ×` last term | +0.7616 / +0.7673 ; +0.7568 / +0.7630 | gap |
| | `a² ± 3.3 ×` last term | +0.7672 / +0.7723 ; +0.7512 / +0.7580 | gap |
| | uniform `Γ_c` shift −2000 / −4000 / −8000 (kept) | +0.431 / +0.103 / −0.553 → **failure at `≈ −4.6·10³`** | gap |
| | for comparison: BELOW_FOLD_SPINJ (fitted ρ), `a¹` / `a²` | +0.739 / +0.787 (kept) | gap / 1.3 |
| (336, 32) | `a¹`, ρ_an (×1) | +0.8405 / +0.8405 | 1.3 |
| | uniform shift −2000 / −4000 | +0.540 / +0.087 → failure at `≈ −4.4·10³` | gap (1.102) |
| | for comparison: BELOW_FOLD_SPINJ (fitted ρ) | +0.834 | 1.3 |
| (336, 32) | **`a²`, ρ_an (×0.3) (production)** | **+0.8370** (C3H: no node below 8.3 at J = 32, kept = dropped) | 1.3 |
| | `a² ± 1 ×` last term / `± 3.3 ×` | +0.8362 / +0.8390 ; +0.8350 / +0.8425 | 1.3 |

**Structure-part stress** (`LJ_SSTRUCT`: `I_ac → I_ac(1 + d_ac/m²)` on every engine node `ν < 8.3`, entry (0,0) fixed; C3H kept).
The chain is **20× more sensitive to the structure part than to the common part** (it splits `++` from `+−`):

| state | `(d₀₁, d₁₀, d₁₁)` [c₂] | `λ/N_s` (min, at `m/Λ`) |
|---|---|---|
| (336, 28) | (−50, −60, −25) ≈ the measured direct − engine offset (§3.1) | +0.799 (1.5) |
| | (−150, −180, −75) / (−500, −600, −250) | +0.813 (1.5) / +0.815 (1.5) |
| | (+66, +77, +33) (the measured offset at ν = 5.0, opposite sign) | +0.530 (gap) |
| | (+150, +180, +75) / (+250, +300, +125) | +0.229 / −0.125 (gap) → **failure at `d₀₁ ≈ +215`** |
| (336, 32) | (−224, −83, −82) (measured at ν = 10) / ×3 | +0.945 (1.5) / +0.938 (1.2) |
| | (+672, +249, +246) | −0.308 (gap) → failure at `d₀₁ ≈ +490` (linear) |

* Every measured structure discrepancy at ν = 5.0 (n4 thimble node), 8.33 and 10.0 is **negative** (J = 28: `−66/−77/−33 ± 5`,
  `−32/−33/−16 ± 6`, `−104/−123/−48 ± 4`; J = 32: `−94/−37/−40`, `−224/−83/−82`), i.e. in the direction that *raises* `λ`.
* The engine's own structure terms in the window are `≤ 36 c₂` on ν ≤ 5 and `≤ 60` at 7–8 (larger only in the narrow zero region
  ν ≈ 6, which the NL convention weights little: the `±3.3×` stress, which multiplies those `120–300 c₂` terms, moves `λ` by `≤ 0.013`).
* So failure needs a **uniform positive** structure error `≈ +215 c₂` (J = 28) / `+490` (J = 32): `3–6×` the engine's structure terms
  and of the opposite sign to every measured discrepancy. This is the binding margin of the window, not the common-part `−4.6·10³`.
* The thimble node at ν = 5.0 cannot test the common part: its scalar has a 23 % MC error (CV factor `0.78 + 0.15i`), so its common
  residual (`−2.2·10⁴ ± 3.9·10³` after CV) is noise; its structure combination is clean (`±5`).


**Reading.**
* With the normalisation fixed analytically the window result **no longer depends on the engine order** (J = 28: `a¹ → a²` moves
  `λ` by +0.002, against +0.048 with the fitted constant) and is stable under `±3.3×` the last term (`0.751 … 0.772`).
* Common part: failure needs an *additional* uniform `Γ_c ≈ −4.5·10³ c₂`, against an engine error `≲ 300 c₂` (last term) and
  `≤ 2 c₂` (analytic constant): margin `≥ 15`.
* Structure part: failure needs a uniform positive structure error `≈ +2·10² c₂` (J = 28), against engine structure terms `≲ 36–60`
  and measured discrepancies of the opposite sign: margin `≈ 3–6` in magnitude, and the sign is favourable.
* J = 32 at order `a²` behaves as J = 28: `a¹ → a²` moves `λ` by −0.002 (−0.0035 including the 24 `c₂` change of ρ_an's radial
  scale), the series stresses span `0.835 … 0.843`. (At order `a¹` alone the `±1 ×` stress is not an error estimate: the last term is
  then the full `a¹` normalisation term, `≈ 5·10³ c₂`, and `a¹ ± (a¹ − a⁰)` spans `+0.23 … +0.96`; `cache/lj_rs_chain32.log`.)




## 5. High J: `Γ − κ` in the H_win window from the engine, fw conventions, `ρ_an` (`lj_rs_gamma.py`) [NE: negative]

* Done for the two high-J rows that have fw `K₂` helicity bases and direct anchors: **(1100, 100)** (Δ_g = 1000, window `ν < 16.3`,
  `T_B(c=3) = 1409`) with the cached depth-2 engine, and **(364, 64)** (Δ_g = 300, window `ν < 8.3`, `T_B ≈ 255` with `W_X`) with a new
  depth-2 engine (fw patched basis, N = 4, 5; 3.7 h). Engine nodes injected below 13.0 / 9.9 into `fw_gamma` exactly as `bj_gamma`
  (a_h from the fw anchors `x = 0.1, 0.2`; fw `K₂` bases), but with `ρ := ρ_an`.
* `ρ_an` at high J needs `fw_frame`'s cure of the `ε²4^J` Fourier round-off in the point-wise direct integrand (without it the
  point-wise ratio scatters by `O(1)` at J = 100); with it the ratio is constant to `7·10⁻⁴` (max) / `8·10⁻⁵` (median), J = 100, and
  `ρ_an` is ν-independent to `4·10⁻¹²`. (`pw_grav`'s closed-form contraction was tried as the direct side: it agrees with w0 to `2·10⁻⁹`
  at (336, 28) but loses all precision at Δ = 2·10⁴, `t − 1 ~ 10⁻⁹`; not used.)
* The fitted constants that `bj_gamma` would use scatter by `±4·10⁴` (1100) and `−7·10³ … +8·10³` (364) between the direct nodes
  ν = 11.7–16.7: the engine is not converged there, so no fitted calibration is possible at high J. `ρ_an` removes that problem but not the
  next one.
* Results (C2; `max` over the window without the de-shift zeros ν ∈ (0.75, 1.5) ∪ (5.5, 7.5), where the local `Γ` of the fw convention
  is singular at every J):

| row | engine `max|Γ − κ|` (`a¹`) | after the `c₀D_h + s xD_h′` class (`a¹`) | last term `a¹ − a⁰`: max / after its own class | `T_B(c=3)` |
|---|---|---|---|---|
| (1100, 100) `++` / `+−` | 7.4·10³ / 6.7·10³ | 1.5·10⁴ / 6.2·10³ | 2.3·10⁴ / 6.9·10⁴ ; 2.3·10⁴ / 3.1·10⁴ | 1409 |
| (364, 64) `++` / `+−` | 3.8·10³ / 4.9·10³ | 6.9·10³ / 6.3·10³ | 7.1·10³ / 1.1·10⁴ ; 6.7·10³ / 1.4·10⁴ | 255 |

  (Excluding (5.0, 8.0) instead: remainders 1.6·10⁴ / 3.6·10³ and 4.7·10³ / 5.0·10³.) Per-node values: `cache/lj_rs_gamma.txt`,
  `cache/lj_rs_gamma_J*_K*.pkl`.
* **Reading.** At `J²/Δ = 9–11` the engine's window values after the D_h class are `4–25×` the H_win tolerance and its last computed
  term after the same removal is `20–55×` it, and the structure series (§1.2) is not convergent there (structure `a¹` alone is `10²–10³ c₂` in the window at (364, 64)).
  The engine at the available depth **cannot decide H_win at (1100, 100) or (364, 64)**, in either the sign-free or the refined form.
  It gives no number that could be compared with `T_B`. The large values in the table are engine truncation, not a measured violation.
* Not computed: J = 34, 48 (300) and J = 150, 200 (1000). There are no fw/`K₂` basis rows for (334, 34), (348, 48), (1150, 150),
  (1200, 200), so `Γ − κ` cannot be formed in the H_win conventions, and §1.2's structure asymptotics (worse with `J²/Δ`) says the
  engine would not be accurate enough at J ≥ 48 even with rows. At J = 34 (`J²/Δ = 3.5`, window ν < 14.4) the engine might reach the
  common part to `≲ 10²` on ν ≲ 8, but not the window above ν ≈ 10, where the structure terms are `10²–10³`.


## 6. Regrading

| item | before (UNIFIED_THEOREM §3.2 / RU) | now | basis |
|---|---|---|---|
| `(336, 28)` window `ν < 8.3` | [M] (fitted constant; entry misses 3–16σ) | **[NE + M]** | §2–§4: no fitted parameter (ρ_an, `≤ 2 c₂`); Prop. RS bookkeeping at J = 28 [DV] (cycle question below the fold reduced to RD); common part reproduced parameter-free at 1–1.5σ on 10–16.7; structure misses = engine truncation (asymptotic series), except an unexplained `≲ 60 c₂` offset at 8.3–10; `λ/N_s = +0.759 / +0.765` (C3H kept / dropped), `0.751 … 0.772` under the series stresses, `+0.80` with the measured structure offset applied; failure needs an extra uniform common `Γ_c ≈ −4.6·10³` or a uniform positive structure error `≈ +2·10²` (`3–6×` the engine's structure terms, opposite in sign to every measured discrepancy). M = the chain (HYPOTHESES_CLOSE §5), C3H below 8.3 not computed by the engine, RD, and the structure margin. |
| `(336, 32)` window | [M] (no validated datum below 8.3) | **[NE + M]** | same construction; `λ/N_s = +0.837` at `a²` (`0.835 … 0.843` under the series stresses; `+0.945` with the measured structure offset, `a¹` engine); failure needs a common shift `≈ −4.4·10³` or a uniform positive structure error `≈ +5·10²`; the J = 28 truncation pattern is reproduced out of sample (§3.2). No C3H node below 8.3 (HC2's thimble seeds failed). |
| `H_LJ(300)` Regge window, J ≥ 34 | [I] | **[I]** (unchanged) | §5: at J = 64 the `a¹` engine's window values after the D_h class are `≈ 25 × T_B` and its last term `45–55 × T_B`; no rows at J = 34, 48; J = 100, 150 not run. The engine cannot decide it. |
| `H_LJ(300)` bulk | [D-sketch/M] | unchanged | not addressed |
| `H_LJ(1000)` | [D-sketch + M] | **[D-sketch + M]** (unchanged) | §5: at (1100, 100) the `a¹` engine's window values after the D_h class are `4–11 × T_B` and its last term `20–50 × T_B`; it neither supports nor violates H_win. The refined H_win remains supported only by HWIN_THIMBLE's direct nodes at ν ≥ 13 [NE]. |
| H_RS for Result D's window (J = 28, 32) | "not covered" (HRS_CLOSE §7) | **[D given RD]** at J = 28, 32 for the `1/Δ` series (B1–B4 [DV]; B6 for the J-dependent factors not checked) | §2.2. Note: at fixed J/Δ the structure series is asymptotic (ratio ≥ 1 at ν ≳ 5), so "all orders" means the asymptotic series. |

**Net for the theorem.** The two computed large-J states at `Δ_g = 300` are now [NE + M] like the other computed states, with margin
`≈ 0.75 N_s`. **`H_LJ(300)` stays [I]** because its window clause for J ≥ 34 is not decided, and **`H_LJ(1000)` stays
[D-sketch + M]**. UNIFIED_THEOREM §2's "two items are [I]" becomes, for D, one item: the `H_LJ(300)` Regge window for `J ≥ 34`.
(This note changes no existing document; the wording is a recommendation.)

## 7. What remains

1. **The `J`-enhanced structure series.** At fixed `J/Δ` it is asymptotic, not convergent (ratio `≈ 1–2` at ν ≳ 5 for
   `J²/Δ = 2.3`; residuals of the truncated engine track the omitted terms, §3.1). Closing `H_LJ(300)` for J ≥ 34 needs a genuine
   resummation of it, e.g. a uniform-in-`J²/Δ` treatment of the helicity structure of the conserved `⟨TTO_J⟩` contraction, or the
   per-wave thimble with a working below-fold contour (HWIN_THIMBLE §7 item 1). Neither was achieved here [I].
2. **The structure offset at ν = 8.3–10** (`−20 … −60 c₂` at J = 28, `−40 … −220` at J = 32), 3–15σ on the precise structure errors,
   not reproduced by any engine order. Small against the window's failure scale, but unexplained: a below-fold systematic of the z2
   nodes, or a beyond-all-orders piece [I]. The cheapest test: the same structure combination on the n4 thimble nodes or on a second
   n3 seed at ν = 8.3, 10.
3. **C3H below 8.3** at (336, 28/32): still the thimble/n3 nodes or nothing (kept/dropped differ by `≤ 0.006`); the engine has no
   C3H integrand.
4. **RD** (Regge dominance of the crystal integral) for the EE at J = 28 below the fold: [D-sketch], as in HRS_CLOSE.
5. **B6** (no boundary layers) for the J-dependent factors (`√(V²W²)⁻¹`, `D̂`) at J = 28, 32.
6. High J: rows (`K₂` bases, anchors) for (334, 34), (348, 48), (1150, 150), (1200, 200) and depth ≥ 3 engines there; neither helps
   without item 1.


## 8. Reproduce (`ads4/code/spinfunc/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)

```
# analytic normalisation (needs the cached bj pieces; references: any depth)            ~25 min per J
python3 lj_rs_norm.py 28 336 3 4 2                    # -> cache/lj_rs_norm_J28_D336_d3_N4_K2.{txt,pkl}
BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 python3 lj_rs_eeq.py 20000 32 2 4     # reference pieces for J = 32 (3 min) -> cache/lj_rs_bj_*
python3 lj_rs_norm.py 32 336 2 5 1
BJ_PATCH=1 python3 lj_rs_norm.py 100 1100 2 4 1 p ; BJ_PATCH=1 python3 lj_rs_norm.py 64 364 2 4 1 p
# bookkeeping B1-B4 at J = 28 (12 min)
python3 lj_rs_book.py 28 336 3 4 2 ; python3 lj_rs_book.py 32 336 3 4 2   # -> cache/lj_rs_book_J{28,32}_D336_d3.txt
# validation, entry/structure resolved, no fit (2 min)
python3 lj_rs_val.py 28 3 4 2 ; python3 lj_rs_val.py 32 2 5 1 ; LJ_VMIN=4 LJ_VTAG=_win python3 lj_rs_val.py 28 3 4 2
# window chains (unchanged hc_chain z4; ~3 min per pass)
./lj_rs_queue28.sh                                    # J = 28: a2, a1, +-1, +-3.3 x last term, uniform shifts -> cache/lj_rs_chain28.log
./lj_rs_queueR.sh                                     # J = 32 reference pieces, norm, val, a1 chains     -> cache/lj_rs_chain32.log
./lj_rs_queue32d3.sh                                  # J = 32 depth 3 (a^2) images (~4 h), val, chains
./lj_rs_queueS.sh                                     # structure stresses                                 -> cache/lj_rs_chainS.log
# high J: (364, 64) depth-2 engine (fw basis), normalisation, fw_gamma-convention Gamma - kappa
./lj_rs_queueH.sh ; ./lj_rs_queueG.sh                 # -> cache/lj_rs_gamma.txt, cache/lj_rs_gamma_J*_K*.pkl
```
Single steps: `LJ_FEXT=f python3 lj_rs_chain.py nl J depth N Kr TAG` (engine `I_K + f(I_K − I_{K−1})`), `LJ_SSTRUCT=d01,d10,d11 …`
(structure stress), then `BJ_C3H=keep|drop [LJ_GSHIFT=g LJ_TAGSFX=_g] python3 lj_rs_chain.py win J TAG` with `source lj_rs_env<J>.sh`.
`LJ_KR=K BJ_NUCUT=ν FW_EST=dir python3 lj_rs_gamma.py J:Δ:depth:N:Nref lo hi`.
