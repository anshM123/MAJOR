# H_RS check: is the analytic-radial engine's EE `k = 4` image the crystal-cycle integral beyond leading order? (Result Q's binding `++` J = 0 gap state)

*22 Nov 2026, H_RS-check agent. Unrefereed, not committed. No existing document or code file was edited. New code:
`ads4/code/spinfunc/rs_points.py`, `rs_engine.py`, `rs_compare.py`, `rs_prop.py`, `rs_queue{A,B,B3,C,D,E,F,G}.sh`. Caches:
`spinfunc/cache/rs_*` (git-ignored). Read-only imports: `z2_points` (EE_DG300_SMALLNU's direct-node protocol, unchanged; only its
in-process patches are used), `n3_hybrid`, `n3_flow`, `n3_run`, `n4_thimble`, `n4_run`, `n1_pw`, `w0_ee_common`, `l5_spin`, `l1_crystal`,
`bg_eeq`, `bg_qan`, `bf_qan`, `bf_exact`, `n4_exact`, `n3_exact`, `exactskel/exactmom`, `resultq/stage6_*`, `common_q`. Existing caches
(`bg_eeq_*`, `bg_qan_T_*`, `bf_exact.pkl`, `n4_exact.pkl`, `z3_exact.pkl`, the Stage 6/7 weights) were read and never written. Compute:
≤ 2 worker processes, `ulimit -v 3 GB`, `OMP/OPENBLAS = 1`, ≤ 0.5 GB per process, ≈ 8 CPU-h. The short analysis runs were made with one
worker paused (SIGSTOP). No `cl_*` file or process was touched.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]** open.
Unit: **`c₂` per harmonic** = `Γ/Σ_j|R_j(ν)|` (BELOW_FOLD2 §4.2's unit; Stage 7's tolerances are quoted in it). `s_d = 0.04 Σ_ν χ_d(ν)
Re Γ_K(ν)` (Stage 7 §2). **Engine** = BELOW_FOLD2's analytic-radial `1/Δ` engine on the real section D (`bg_eeq` caches, K = 4).
**Crystal** = the crystal-cycle integral evaluated with the direct-node protocol of EE_DG300_SMALLNU / HYPOTHESES_CLOSE(2): the n3 hybrid
contour above the fold and the n4 `A₊ + A₋` Lefschetz thimbles below it.

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **Method** | The z2 direct-node protocol, unchanged, with one new integrand kind `EE4` = w0 EE (`k = 2`) × the Lemma-B factor `R = ((v′−u′)²−1)/(4u′v′)`. This is exactly the integrand of the engine's `EE4` image. Each node is a state/reference pair (`Δ_ref = 2·10⁴`) on common random numbers. The exact scalar is computed on the same samples. | — |
| **Scalar validation** (exact Mellin skeleton, same samples) | n3 nodes (`ν ≥ 10`): all within 1.8σ, `|Δc₂| ≤ 2`. n4 thimble nodes (`ν = 2–8`): within 2.4σ. At `ν = 6.02` the check is ill-conditioned (zero of `g + g̃` at 5.98). The n3 contour **below** the fold fails (`ν = 6.02`: `+2800 ± 95`, 30σ), so it is not used there. χ² (scalar): 26/28 (Δ = 1000), 26/32 (300), 22/20 (100), 17/10 (J = 2). | [NE] |
| **Above the fold, `ν = 16–25`** | Crystal ≡ engine **beyond LO**. The engine's whole `Γ_exact` there (−10 … −209 per harmonic) is reproduced to **±0.5–1.8** (Δ = 100, 300) and ±3.7–5.4 (1000). Example, ν = 20.02: `Γ_exact = −73.4`, crystal − engine `= −0.3 ± 1.5` (300), `−0.1 ± 0.6` (100), `+0.2 ± 1.7` (J = 2). | **[NE]** |
| **`ν = 10–14`** | First seed: a smooth **−3 … −28 per harmonic** deviation, nominally up to 4.9σ. **It is a sampling fluctuation, not a discrepancy:** (i) it scales like `Δ` in c₂ (−2.9 / −9.0 / −27.9 at Δ = 100 / 300 / 1000, ν = 12.02), i.e. like `1/Δ` in the log-ratio, which is the signature of a shared Monte Carlo fluctuation of the pointwise `O(1/Δ)` structure; (ii) independent seeds give `−1.6 ± 1.4`, `+5.1 ± 4.4`, `+2.7 ± 2.1` (4× points) at 300, and `−0.1 ± 0.7` (100), `−6.1 ± 4.7`, `+15 ± 14` (1000); (iii) it does not depend on the contour (`ν_f` = 12.4 / 13 / 14.5). The `k = 2` image on the same samples agrees to **±0.3** (the R factor adds heavy-tailed variance). Lesson: the 8-replica jackknife **underestimates the `k = 4` error by about ×3** (seed scatter χ² = 16/2). Every bound below uses ×3. | [NE] (consistent) |
| **Below the fold, `ν = 2–8`** (n4) | Consistent at every node: ≤ 1.9σ before the ×3 inflation, ≤ 0.6σ after. But the errors (×3) are **±39 (7.02) … ±520 (3.02) … ±1800 (2.02)** per harmonic at Δ = 300, about 3× worse at 1000 and about 3× better at 100. **The c₂ error grows like Δ**: Δ_g = 1000 is the *least* sensitive state, not the cleanest. All central values below ν = 6 are positive (correlated: shared seed). | consistent, uninformative below ν ≈ 6 |
| **`ν < 2`** | Not reached. n4 tracks down to 2.02, but the common-random-number correlation is lost there (Δ = 100: 30 % per-block, ±6000). | [I] |
| **Weight coverage** (H0 `0.04 Σ|χ_d|Σ|R_j|` = 0.171) | ν < 2: **68 %** of the weight (and 131 % of the engine's `s_d`: +0.67 of +0.51). 2–4: 21 %. 4–6: 4 %. ≥ 6: 7 %. The uniform per-harmonic error that breaks δ = 0.2 at the binding state is **17 (ν ∈ [1,2]) … 65 ([3,4]) … 200 ([4,6]) … 266 ([6,8]) … 344–980 (ν ≥ 8)**. For δ = 0.05 it is 0.9 … 3.6 … 11 … 15 … 19–55. | — |
| **Propagated `Δs_d`, `ΔFF/N`** (binding `++`/OO J = 0, gap) | On the window **`ν ∈ [5, 25]`**: `Δs_d` = **−0.10 ± 0.35** (Δ_g = 300, H0, ×3 errors; SB +0.02 ± 1.08), −0.52 ± 0.97 (1000). That is `ΔFF/N(δ = 0.2)` = **+0.011 ± 0.039** (300) / +0.058 ± 0.11 (1000), against the margin `s_fail − s_d = 1.27`. On `ν ∈ [2, 25]`: `Δs_d = +3.9 ± 10.7` (300), i.e. **no information**. `ν < 2` is not tested. | [NE] (ν ≥ 5) / [I] |
| **J = 2** (`++`, 1.05Λ, Δ = 314.935) | Same picture: ν = 10–20 crystal ≡ engine (`ν = 20.02`: `+0.2 ± 1.7` against `Γ_exact = −72.0`); 3.02 / 5.02: `+18 ± 166` / `+57 ± 44` (jackknife). | [NE] above the fold / [I] below |
| **"512 c₂"** | This is BELOW_FOLD2 §3.3's EE **`k = 2`** (not `k = 4`) J = 0, Δ = 400 comparison with the `K₂` Laplace rows. BELOW_FOLD3 §2 traced it to the plain-sum `log(Σ M_k)` resummation in the comparison script (`bg_k2check`). In Result D's chain form `K₂` = engine to 2 c₂. It is the same item as the "log(ΣM_k) artefact"; the review's issue 3 bullet lists it as an unresolved `k = 4` discrepancy, which mislabels it and misses the resolution. It is **positive** (Laplace-at-the-crystal) evidence for H_RS at ν = 20–40, `k = 2`, and says nothing below the fold. | [DV] (from BF3, re-read) |
| **H_RS regrade** | Split by ν. **ν ≥ 10: [NE]** (J = 0, 2; `k = 4` and `k = 2`; beyond LO). **6 ≤ ν < 10: [NE, weak]**: bounds sufficient for δ = 0.2 there, not for δ = 0.05. **2 ≤ ν < 6: [I]**: consistent, but errors are 3–30× the δ = 0.2 thresholds. **ν < 2: [I], untested**. For Result Q's binding state H_RS **stays [I]**: the decisive 89 % of the weight (ν < 4) is not constrained. | — |

**Bottom line.** This is the first independent crystal-cycle test of the engine's EE `k = 4` image beyond leading order.
* **Where it can be done** (ν ≥ 6, J = 0 and 2, Δ = 100–1000), the engine passes. Above the fold it passes to the 1 % level of its whole beyond-LO content.
* **Where Result Q needs it** (ν ≲ 4, which carries 89 % of the `s_d` weight and more than all of the measured `s_d = +0.51`), no crystal-cycle method available reaches the precision. The needed relative precision of the state/reference ratio there is 10⁻⁵–10⁻⁴; the n4 thimbles deliver 10⁻³–10⁻², and the n3 contour is invalid below the fold.
* So the δ = 0.2 binding margin (×3.5) is neither confirmed nor refuted by this check. **The review's issue 3 stands for ν < 4.**

---

## 1. What is compared

### 1.1 Observable [D]
* The engine and the crystal both give the state/reference ratio `ρ(ν) = I_st(ν)/I_ref(ν)` of the EE `k = 4` image (J = 0, 2: one conserved
  `⟨TTO⟩` structure, so the prefactor is 1×1).
* Result Q's chain (`bg_qan.gamma_K`) uses only `Γ_exact = m²(c T_ref ρ − T_st)`, with `c` fixed at ν = 0.02. Then
  `Γ_K = Γ_exact − 0.75 Σ_j R_j(ν − ij)² − const·N̂₄`.
* Per node we form `d(ν) = log ρ_crystal − log ρ_engine` (complex). If H_RS holds, `d` is a ν-independent constant κ: the engine drops
  Δ-only constants (the overall `a`-power, `β_{Δ,J}`). κ is fitted (weighted, complex); the fit is dominated by the ν = 16–25 nodes.
* The difference in the chain's unit is `dG(ν) = Re[m² T_st(ν)(d(ν) − κ)]/Σ_j|R_j(ν)|` (c₂ per harmonic).
* **κ cannot fake or hide an `s_d` effect** [DV]:
  * `T_st/N̂₄ = 1` to 10⁻³ on the whole grid, so a κ shift is a multiple of `N̂₄`;
  * the chain removes exactly that class, and `0.04 Σ χ_d N̂₄ = 0` to 10⁻¹³ in all six gauges (`rs_prop` prints it).
  So the unmeasured ν = 0.02 normalisation and κ drop out of `s_d`.
* **"Entry by entry".** At J ≤ 2 there is one entry per ν. Its two components, `Re d` and `Im d`, are both reported.
  `Re Γ_K`, hence `s_d`, uses essentially `Re d` (`T_st` is real here). `Im d` is a pure consistency test.

### 1.2 Precision requirement [D]
* `dG` error = `(m²|T_st|/Σ|R_j|) × δ(log ρ)`. The conversion factor is `0.9·10⁶` (Δ = 1000), `0.8·10⁵` (300), `0.9·10⁴` (100) at
  ν = 2–3. It is 10–100× smaller near the zeros 1.07, 4.55, 8.17, 13.09 of `N̂₄`.
* The δ = 0.2 thresholds of §4 (17–80 per harmonic at ν < 4) therefore need `δ(log ρ) ≈ 10⁻⁵` (1000), 10⁻⁴ (300) or 10⁻³ (100).
* The common-random-number (CRN) ratio error comes from the pointwise `O(1/Δ)` difference of the state and reference integrands, so
  `δ(log ρ) ∝ 1/Δ` at fixed sampling, and **the c₂ error grows like `m²/Δ ∝ Δ`**. Δ = 100 is the most sensitive test of the same
  integral identity. Δ = 1000 is the least sensitive, although its engine truncation is the cleanest.

## 2. Method and scalar validation

### 2.1 Crystal nodes (`rs_points.py`)
* `z2_points.n3_pair` / `n4_pair`, imported unchanged:
  * n3 hybrid contour `C(ρ₁ = 3.5, ν_f)` with `ν_f = max(ν, 15)` for ν ≥ 13 and `ν_f = 13` below;
  * n4 `A₊ + A₋` thimbles, saddles tracked from ν_c = 10.5;
  * proposal adapted on the state and re-used for the reference, with the same Sobol seeds (CRN);
  * 2¹⁰·2·8 points (n3) or 2⁹·2·8 (n4) per block; error = jackknife over the 8 scramblings.
* New kind `EE4` = w0 EE × `R` (computed from the contour point's raw `z`; `R` is rational, so no branch replay is needed).
* `bg_ee` builds the engine's `EE4` from the same w0 builder and the same `R` (pointwise checked in BELOW_FOLD2 §1.3).
* Nodes sit on the engine grid (`ν = k + 0.02`).
* **One file per state** (`rs_pts_EE4_J<J>_D<Δ><tag>.pkl`). In a first version, states sharing a file overwrote each other's reference
  records (the reference key `(2·10⁴, ν)` is state-specific under CRN). The affected Δ = 1000 n3 pairs were re-run; the Δ = 1000 n4
  pair at 8.02 was dropped. The collided files are kept, renamed `*.collided.pkl`, and are not used.

### 2.2 Scalar validation (same samples; `d_s = log ρ_MC − log ρ_exact`, `Δc₂ = m² Re d_s · ex₂/E₂` as in z3_nodes)

| Δ | n3, ν = 10–25 (10 nodes) | n4, ν = 2–8 | notes |
|---|---|---|---|
| 1000 | `|Δc₂| ≤ 1.6 ± 1.1`; `|Re d_s| ≤ 1.6·10⁻⁶` | 2.02: −99 ± 126; 3.02: +1 ± 22; 4.02: +10 ± 24; 5.02: −104 ± 105; 7.02: **−29 ± 12 (2.4σ)** | 6.02: +2717 ± 2403 (ill-conditioned) |
| 300 | `|Δc₂| ≤ 1.7 ± 1.0`; `|Re d_s| ≤ 2·10⁻⁵` | 2.02: −768 ± 694; 3.02: −23 ± 24; 4.02: −8 ± 40; 5.02: −16 ± 24; 7.02: −6 ± 10; 8.02: −51 ± 48 | 6.02: +1965 ± 4832 (ill-cond.) |
| 100 | `|Δc₂| ≤ 1.8 ± 1.0` | 2.02: +585 ± 820; 3.02: +108 ± 128; 5.02: +14.5 ± 7.3 (2.0σ) | 2.02: CRN lost |
| 314.9 (J = 2) | `≤ 1.5 ± 1.5` | 3.02: −8 ± 25; 5.02: −39 ± 21 (1.9σ) | |
| pilot, n3 below the fold (Δ = 1000, `ν_f = 13`) | — | 4.02: +3042 ± 1581; **6.02: +2801 ± 95 (30σ)** | n3 is **not** valid below ν ≈ 8 |

**The method passes the exact scalar where it is used** [NE]. The scalar's log-ratio residual, converted with the EE factor (column
"scal-sys" of `rs_compare`), is ≤ 20 per harmonic on every used node except 2.02 (130 / 1200) and the ill-conditioned 6.02.

## 3. Entry-by-entry comparison (crystal − engine)

Columns:
* engine `Γ_exact` and `Γ_K` (per harmonic; the a⁴ term is ≤ 0.3 except where noted);
* `dG ± σ_jk` (per harmonic; multiply σ by 3 for the error model of §3.5);
* `Re/Im(d − κ)` in units of 10⁻⁴.

Seed 4242 for all main nodes. Files: `cache/rs_compare_D*_J*_EE4_n3_n4.txt`.

### 3.1 Δ_g = 1000, `++` J = 0 at the gap (Δ = 1000)

| ν | meth. | `Γ_exact` | `Γ_K` | **dG ± σ_jk** | Re(d−κ) | Im(d−κ) |
|---|---|---|---|---|---|---|
| 2.02 | n4 | −7.64 | −10.41 | +278 ± 690 | +3.0 ± 7.5 | +11.6 ± 6.7 |
| 3.02 | n4 | +7.69 | −7.58 | +650 ± 666 | +8.2 ± 8.4 | −20.2 ± 9.0 |
| 4.02 | n4 | +14.13 | −2.14 | +510 ± 376 | +23 ± 17 | −10.5 ± 8.8 |
| 5.02 | n4 | +7.86 | +0.87 | +188 ± 183 | −16 ± 16 | +10 ± 14 |
| 6.02 | n4 | −1.10 | +1.45 | +18 ± 117 | −1.0 ± 6.8 | +2.2 ± 2.9 |
| 7.02 | n4 | −6.15 | +0.93 | +67 ± 42 | −6.6 ± 4.1 | −4.5 ± 5.6 |
| 10.02 | n3 | +1.51 | −0.58 | −56 ± 36 | −7.2 ± 4.7 | +4.1 ± 2.2 |
| 11.02 | n3 | +5.46 | −0.59 | −40 ± 15 | −5.5 ± 2.0 | +2.0 ± 1.7 |
| 12.02 | n3 | +7.70 | −0.42 | −27.9 ± 5.7 | −6.1 ± 1.3 | +1.1 ± 2.2 |
| 13.02 | n3 | +7.53 | −0.12 | −23.7 ± 11.0 | −68 ± 32 | +6 ± 30 |
| 14.02 | n3 | +4.57 | +0.25 | −12.4 ± 5.1 | +2.7 ± 1.1 | +0.3 ± 2.0 |
| 16.02 | n3 | −10.12 | +1.08 | −5.5 ± 5.4 | +0.35 ± 0.35 | +0.35 ± 0.44 |
| 20.02 | n3 | −73.18 | +2.70 | −1.3 ± 4.5 | +0.04 ± 0.13 | +0.09 ± 0.16 |
| 25.02 | n3 | −206.5 | +4.35 | +0.3 ± 3.7 | −0.01 ± 0.07 | −0.05 ± 0.09 |

χ² (Re, Im of `d − κ`; κ fitted on all nodes): n4 17.9/12, n3 51.7/16 (the seed-4242 fluctuation at 10–14, §3.5).

### 3.2 Δ_g = 300, `++` J = 0 at the gap (Δ = 300)

| ν | meth. | `Γ_exact` | `Γ_K` | **dG ± σ_jk** | Re(d−κ) | Im(d−κ) |
|---|---|---|---|---|---|---|
| 2.02 | n4 | −7.64 | −10.41 [a⁴ −0.18] | +414 ± 607 | +50 ± 74 | +10 ± 82 |
| 3.02 | n4 | +7.71 | −7.56 | +142 ± 173 | +20 ± 24 | **−127 ± 54** |
| 4.02 | n4 | +14.16 | −2.11 | +162 ± 87 | +83 ± 45 | −25 ± 31 |
| 5.02 | n4 | +7.88 | +0.89 | +67 ± 70 | −64 ± 67 | +45 ± 68 |
| 6.02 | n4 | −1.11 | +1.45 | −29 ± 43 | +19 ± 28 | +27 ± 19 |
| 7.02 | n4 | −6.16 | +0.91 | +18.9 ± 13.1 | −21 ± 14 | −18 ± 19 |
| 8.02 | n4 | −6.24 | +0.19 | +19.3 ± 20.2 | −186 ± 190 | +3 ± 78 |
| 10.02 | n3 | +1.51 | −0.58 | −18.5 ± 12.6 | −27 ± 18 | +15 ± 8 |
| 11.02 | n3 | +5.47 | −0.58 | −13.1 ± 5.2 | −20 ± 8 | +7 ± 6 |
| 12.02 | n3 | +7.72 | −0.40 | −9.0 ± 2.1 | −22 ± 5 | +4 ± 8 |
| 13.02 | n3 | +7.55 | −0.10 | −7.8 ± 2.9 | −251 ± 95 | +20 ± 91 |
| 14.02 | n3 | +4.59 | +0.26 | −3.7 ± 1.7 | +8.8 ± 4.0 | +2.2 ± 6.8 |
| 16.02 | n3 | −10.14 | +1.06 | −1.4 ± 1.8 | +1.0 ± 1.3 | +1.3 ± 1.4 |
| 18.02 | n3 | −36.40 | +1.83 | −0.7 ± 1.7 | +0.3 ± 0.7 | +0.7 ± 0.7 |
| 20.02 | n3 | −73.40 | +2.49 [a⁴ −0.27] | −0.3 ± 1.5 | +0.09 ± 0.48 | +0.30 ± 0.46 |
| 25.02 | n3 | −207.3 | +3.57 [a⁴ −1.03] | +0.3 ± 1.2 | −0.05 ± 0.25 | −0.24 ± 0.26 |

χ²: n4 18.5/14, n3 48.2/18. The 3.02 `Im` value (2.4σ) recurs at Δ = 1000 (2.2σ) and at J = 2 (2.6σ). All three use the same seed
and nearly the same reference thimbles, so this is one ≈ 2.5σ observation. It does not enter `Re Γ_K`.

### 3.3 Δ = 100, J = 0 (sensitivity state; engine built here, `rs_eng_D100_J0_K4.pkl`)

| ν | meth. | `Γ_exact` | `Γ_K` [a⁴] | **dG ± σ_jk** |
|---|---|---|---|---|
| 2.02 | n4 | −7.41 | −10.18 [−1.16] | −7000 ± 6100 (CRN lost: 30 % per block) |
| 3.02 | n4 | +7.80 | −7.47 [−0.42] | +87 ± 131 |
| 5.02 | n4 | +7.84 | +0.85 [+0.31] | +76 ± 89 |
| 10.02 / 11.02 / 12.02 / 14.02 | n3 | +1.5 / +5.5 / +7.7 / +4.6 | −0.6 / −0.6 / −0.4 / +0.3 | −7.0 ± 5.4 / −4.4 ± 2.4 / −2.9 ± 1.1 / −1.0 ± 0.7 |
| 16.02 / 20.02 / 25.02 | n3 | −10.1 / −73.6 / −209.1 | +1.1 / +2.3 [−2.0] / +1.8 [−7.6] | −0.36 ± 0.70 / −0.06 ± 0.61 / +0.06 ± 0.52 |

At Δ = 100 the engine's own a⁴ step is −2.0 (ν = 20) and −7.6 (ν = 25), yet crystal = K = 4 engine to ±0.6. The engine's truncation
estimate is conservative there.

### 3.4 J = 2, `++` at 1.05Λ (Δ_g = 300; Δ = 314.935, the BELOW_FOLD2 state)

| ν | meth. | `Γ_exact` | `Γ_K` [a⁴] | **dG ± σ_jk** |
|---|---|---|---|---|
| 3.02 | n4 | +7.78 | −7.49 [−0.86] | +18 ± 166 (`Im(d − κ) = −109 ± 42·10⁻⁴`, 2.6σ, see §3.2) |
| 5.02 | n4 | +7.84 | +0.85 [−0.20] | +57 ± 44 |
| 10.02 | n3 | +1.57 | −0.52 [−0.13] | −17.8 ± 12.1 |
| 16.02 | n3 | −10.27 | +0.93 [+0.53] | −1.3 ± 1.6 |
| 20.02 | n3 | −72.01 | +3.88 [+2.40] | +0.2 ± 1.7 |

χ²: n4 8.9/4, n3 4.6/6.

### 3.5 The ν ≈ 10–14 deviation: seed, contour, `k = 2` and Δ tests (ν = 12.02; `rs_queueD–G`)

| Δ | variant | dG ± σ_jk |
|---|---|---|
| 300 | seed 4242, `ν_f = 13` (main) / proposal adapted on EE4+EE | −9.0 ± 2.1 / −6.3 ± 1.3 |
| 300 | seed 4242, `ν_f = 12.4` / `14.5` | −10.6 ± 3.1 / −13.1 ± 4.8 (11.02: −25 ± 12) |
| 300 | seed 9191 / 5151 / 5151 with 2¹² (4×) points | **−1.6 ± 1.4 / +5.1 ± 4.4 / +2.7 ± 2.1** |
| 300 | **`k = 2` (EE, no R) on the same samples** (seed 4242 / 5151, 4×) | `Re(d − κ) = −0.7 ± 0.7 / +0.3 ± 0.4 ·10⁻⁴` (k = 4: −15.5 ± 3.2 / +6.6 ± 5.2) |
| 1000 | seed 4242 / 9191 / 5151 | −27.9 ± 5.7 / **−6.1 ± 4.7 / +15 ± 14** |
| 100 | seed 4242 / 9191 | −2.9 ± 1.1 / **−0.06 ± 0.71** |

* **Reading [NE].**
  * The seed-4242 points carry a correlated fluctuation of the `R`-weighted pointwise `O(1/Δ)` deviation. It is shared by all Δ, by the
    ν_f variants (same base points) and by the neighbouring ν, and it scales as 1/Δ in `log ρ` (−6.5 / −2.2 / −0.61 ·10⁻³ at
    100 / 300 / 1000).
  * A structural (Γ-level) error would scale as 1/Δ² in `log ρ`, i.e. as a Δ-independent `dG`. That is excluded at Δ = 100 (−0.06 ± 0.7
    against −9 at 300).
  * Independent seeds scatter around zero.
* The `k = 2` image on the same samples has 4–8× smaller errors. The `R` factor (poles at `u′v′ = 0`) produces heavy tails, which the
  8-replica jackknife underestimates. The seed scatter at 12.02 (Δ = 300) gives χ² = 16 for 2 dof, i.e. a factor ≈ 2.8.
* **Error model used for all bounds: σ = 3·σ_jk.** With it every node of §3.1–3.4 is ≤ 1.7σ.
* The means over all seeds at 12.02 (×3 errors) are −2.0 ± 3.0 (300), −0.9 ± 1.8 (100) and −13 ± 11 (1000).

## 4. Propagation into `s_d` and FF/N (`rs_prop.py`; Stage 7 §2 chain, `bg_qan.ffn` read-only)

`dG` is interpolated linearly between nodes. Outside the measured range it is *not tested* and set to 0. Binding state: `++`/OO J = 0
at the gap. Slope `∂(FF/N)/∂s_d = −0.1126` (δ = 0.2) and `−0.0863` (δ = 0.05). Margins: `s_fail − s_d = 1.78 − 0.51 = 1.27` (δ = 0.2),
`0.58 − 0.51 = 0.07` (δ = 0.05).

| window (`ν`) | Δ_g | `Δs_d` H0 (σ_jk / ×3) | `Δs_d` SB (×3) | `ΔFF/N` δ = 0.2 (×3) | `ΔFF/N` δ = 0.05 (×3) |
|---|---|---|---|---|---|
| [5.02, 25.02] | 300 | **−0.10 ± 0.12 / ± 0.35** | +0.02 ± 1.08 | **+0.011 ± 0.039** | +0.009 ± 0.030 |
| [5.02, 25.02] | 1000 | −0.52 ± 0.32 / ± 0.97 | −0.65 ± 2.85 | +0.058 ± 0.110 | +0.045 ± 0.084 |
| [2.02, 25.02] | 300 | +3.9 ± 3.6 / ± 10.7 | +2.3 ± 9.7 | −0.44 ± 1.21 | −0.34 ± 0.93 |
| [2.02, 25.02] | 1000 | +17 ± 13 / ± 40 | +12 ± 33 | −2.0 ± 4.5 | −1.5 ± 3.4 |
| [0, 2.02) | both | **not tested** | | | |

The full-set minimum of Result Q (+0.1426 at δ = 0.2, +0.0060 at δ = 0.05) is therefore **unchanged at the precision of this check**.
* The ν ≥ 5 part of `s_d` is confirmed at Δ_g = 300 to ±0.35 (×3 errors, 2σ = 0.7 < 1.27).
* Nothing is learned about the ν < 5 part, which is where `s_d` lives: engine `s_d` by window (H0) is +0.10 [0,1), +0.57 [1,2), −0.05
  [2,3), −0.10 [3,4), −0.01 [4,8), ≈ 0 above.

**Per-window tolerance vs achieved bound** (H0 weight `0.04 Σ|χ_d|Σ|R_j|`; breaking = the uniform `|dG|` that uses up the margin):

| ν | share of weight | breaks δ = 0.2 / 0.05 at | best achieved 1σ (×3), any Δ | decision-relevant? |
|---|---|---|---|---|
| [0, 1) | 25 % | 30 / 1.7 | — | not measured |
| [1, 2) | 44 % | 17 / 0.9 | — | not measured |
| [2, 3) | 9 % | 80 / 4.5 | ±1800 (300) | no (×20) |
| [3, 4) | 11 % | 65 / 3.6 | ±390 (100), ±520 (300) | no (×6–8) |
| [4, 6) | 4 % | 200 / 11 | ±130 (J = 2), ±210 (300) | marginal for δ = 0.2 |
| [6, 8) | 3 % | 266 / 15 | ±39 (300, 7.02) | **yes for δ = 0.2**, no for 0.05 |
| [8, 10) | 1 % | 980 / 55 | ±61 (300) | yes / no |
| [10, 18) | 4 % | 344–478 / 19–27 | ±2–6 (100, 300) | **yes for both** |
| [18, 41) | < 0.1 % | ≥ 9000 / 500 | ±2–5 | yes |

## 5. The "512 c₂ `K₂`" item

* **Source.** BELOW_FOLD2 §3.3. `bg_k2check.py` compares the engine's EE **`k = 2`** image (`kinds=('EE2',)`; rows
  `w1_pts_EE_J0_k2_b`, Δ = 400, J = 0, ν = 20 → 40) with the `K₂` Laplace rows. The rows were formed as `log(M0 + M2 + M4 + M6)`:
  +512 c₂, and −107 c₂ for the scalar against the exact skeleton.
* **Resolution** (BELOW_FOLD3 §2, `bh_k2forms.py`, `cache/bh_k2forms.txt`, re-read here).
  * The rows' `M2/M0` is a Δ-constant (+0.052) and `M4/M0` a ν-only piece shared by the reference. The plain sum therefore carries the
    spurious cross term `−(A2 − R2)A4 ≈ 600 c₂`.
  * In the cumulant / chain (`l9_analysis.pred`) form, `K₂` = engine to **−2.0 c₂** (J = 0, Δ = 400, ν = 20 → 40). The values at the
    other points are −15 c₂ (Δ = 1600, within the engine's a⁴ step of 16) and −22 c₂ (J = 2, a⁴ step 107).
  * The scalar is +2.5 c₂ against the exact skeleton.
* **So the "512 c₂ discrepancy" and the "log(ΣM_k) comparison-script artefact" are the same thing.** REVIEW_UNIFIED issue 3 ("EE
  `k = 4`, J = 0, Δ = 400: an unresolved 512 c₂ discrepancy against `K₂`") gets the order wrong (it was `k = 2`) and its status wrong
  (resolved in BF3 §2, which the review lists as read). It is not an independent problem for H_RS.
* **Bearing on H_RS.**
  * Once resolved, it is *positive* evidence. `K₂` is a Laplace expansion at the crystal saddle, i.e. a crystal-cycle evaluation to
    relative `O(Δ⁻²)`, which is all of the c₂-level content. It matches the engine beyond LO at `k = 2`, ν = 20–40.
  * This note adds direct all-orders crystal nodes at `k = 4` (ν = 10–25, J = 0, 2) and `k = 2` (ν = 12–20).
  * Neither reaches ν < 10.

## 6. Regrading H_RS

| part of H_RS | before (REVIEW_UNIFIED §4) | now | evidence |
|---|---|---|---|
| scalar, all J | [DV] | [DV] | unchanged (BF2) |
| EE at LO | [DV] | [DV] | unchanged |
| **EE `k = 4` (and `k = 2`), J = 0, 2, ν ≥ 10, beyond LO** | [I] (only the `k = 2` `K₂` check at ν = 20–40, mis-reported) | **[NE]** | §3: `Γ_exact` up to −209 reproduced to 0.3–1 % (ν = 16–25); ν = 10–14 consistent at the ×3 error model (independent seeds); `k = 2` to ±0.3 |
| EE `k = 4`, J = 0, 6 ≤ ν < 10 | [I] | **[NE, weak]** | consistent; bounds ±39–130 (×3) against δ = 0.2 thresholds 266–980 |
| EE `k = 4`, J = 0, 2, 2 ≤ ν < 6 | [I] | **[I]** (consistent, uninformative) | errors 3–30× the δ = 0.2 thresholds |
| EE `k = 4`, ν < 2 | [I] | **[I], untested** | 68 % of the weight; the engine's whole `s_d` |
| J = 28 entry-resolved (Result D window, `bj_*`) | fails at 3–16σ | not re-tested | outside the scope of this check |
| **H_RS as used by Result Q's binding state** | [I] | **[I]** | the decisive window ν < 4 is not constrained |

`s_d = +0.51` stays [NE given H_RS]. The full-set min FF/N (+0.1426 / +0.0060) stays [NE + M given H_RS]. Recommended wording for the
statements: "H_RS (EE `k = 4` real-section ≡ crystal beyond LO) is confirmed numerically above the fold (ν ≥ 10, [NE]) and for
6 ≤ ν < 10 at the δ = 0.2 level ([NE, weak]); for ν < 6, which carries ≈ 90 % of the `s_d` weight, it is untested [I]."

## 7. What remains [I]

1. **ν < 4 at the 10–20 c₂-per-harmonic level** (δ = 0.2) or ≲ 1 (δ = 0.05). Possible routes:
   * (a) lower Δ: the c₂ error ∝ Δ, so Δ = 30–50 gains ×6–10 over Δ = 300. The engine there needs K ≥ 6 or a check of its own
     convergence;
   * (b) more statistics with variance control of the `R` factor. A pointwise regression control variate would be possible, since `k = 2`
     on the same samples is 4–8× quieter;
   * (c) n4 below ν = 2: the saddles track to 2.02, but the CRN correlation collapses there;
   * (d) an analytic proof that the Γ-continuation of the `R`-weighted monomials equals the crystal-cycle integral. `R`'s poles at
     `u′v′ = 0` are the only new singular locus relative to the validated `k = 2` / scalar integrands. That would turn §5.4 of BELOW_FOLD2
     from [D-sketch] into [D].
2. Error model: the jackknife over 8 scramblings underestimates `k = 4` errors by ≈ ×3. Future direct `k = 4` nodes need independent
   seeds.
3. Not covered: `+−` J = 4–6, OO, J ≥ 4, other Δ_g states, and the J = 28/32 window of Result D.
4. The 2.5σ `Im` residual at ν = 3.02 (one seed, three states): repeat with an independent seed.

## 8. Reproduce (`ads4/code/spinfunc/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)

```
# crystal nodes (≈ 2–3 min per n3 pair, 12–25 min per n4 pair); queues as run:
./rs_queueA.sh; ./rs_queueB.sh    # (B was cut after its n4 job; its Delta = 100 part ran as rs_queueB3.sh)
./rs_queueB3.sh; ./rs_queueC.sh; ./rs_queueD.sh; ./rs_queueE.sh; ./rs_queueF.sh; ./rs_queueG.sh 1000; ./rs_queueG.sh 100
#   single node, e.g.: RS_TAG=_n4 RS_NUSPLIT=9.5 Z2_M=9 Z2_R=8 python3 rs_points.py 0 3.02 300
#   (env: RS_KINDS=EE4[,EE], Z2_SEED, N3_NUF, N3_NUF_LOW; Delta 1000 n3 pairs were re-run after the shared-file fix)
python3 rs_engine.py 100 0                      # engine image for Delta = 100 (2 min); 300/1000/314.935 read from bg_eeq caches
# comparison tables (§2.2, §3):
python3 rs_compare.py 0 1000 _n3,_n4;  python3 rs_compare.py 0 300 _n3,_n4;  python3 rs_compare.py 0 100 _n3,_n4
python3 rs_compare.py 2 314.935 _n3,_n4
python3 rs_compare.py 0 300 _n3k2 EE   # k = 2;  seeds / contours: tags _n3s2 _n3s3 _n3m12 _n3f124 _n3f145 _n3f17
python3 rs_compare.py 0 1000 _pilotn3,_pilot   # n3 below the fold fails the scalar (6.02: 30 sigma)
# propagation (§4):
RS_INFL=3 python3 rs_prop.py 300 _n3,_n4 dir [2.02,3.02,4.02]   # (and 1000; RS_INFL=1)
```

Outputs:
* `cache/rs_pts_EE4_J{0,2}_D*_{n3,n4,…}.pkl`, `rs_raw_D*_*.pkl`;
* `cache/rs_compare_D*_J*_{EE4,EE}_*.txt`, `rs_prop_D*_dir_x{1,3}[_ex…].txt`;
* `cache/rs_exact.pkl` (exact scalar for Δ = 314.935, J = 2), `rs_eng_D100_J0_K4.pkl`, `rs_T_D100_J0.pkl`;
* logs `cache/rs_queue*.log`.
