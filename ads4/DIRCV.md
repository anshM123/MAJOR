# DIRCV: the dir/cv estimator systematic at (Δ, J) = (1100, 100), and the corrected H_LJ(1000) leverage

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. `HIGHSPIN.md` and `code/highspin/` were not touched.*
* *New code:* `ads4/code/dircv/`. Caches are in `cache/`; its `.gitignore` ignores everything except itself.
  * `dc_guard.py` is a process-local write guard, imported first by every script. It refuses any open-for-write, replace, rename,
    remove or mkdir outside `dircv/cache`.
  * The project pipelines (`pw_gamma`, `fw_gamma`, `pw_points`, `z3_nodes`, `fw_scalar`, `n2_common`) are imported **read-only and
    unchanged**. Redirection is done by monkey-patching.
  * `run_all.sh` reproduces the light part in about 5 CPU-min. `run_new.sh` and `run_tail.sh` are the MC runs.
* *Compute:* ≈ 2.2 CPU-h, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.
  * New MC data: two independent seeds (1001, 2002) of the full production node (waves + P_J twin + reference, 24 replicas each,
    3× the original) at ν = 13.3 and 16.3. The `_x3`/`_x4` runs at 20.8 and 25 were stopped unused.
  * Per-point samples (8192 per mode) at ν = 13.3 for the tail diagnostics.

Grades as in UNIFIED_THEOREM_v3 §0. `r` is the constant (non-ring) coefficient of N2_SUBREGGE's raw-form fit, in c₂ units. "F1" is the
A + C ring class + r, "F1B" is A + B + C + r. T_B (c = 3) is N2_SUBREGGE §2.1: `++` 1158 (`|Y|`) / 1441 (`|X|`), `+−` 1152 / 1439.
Leverage `L = T_B/(|r| + rms)`, quoted from `|Y|` to `|X|`. "At 2σ" means `T_B^{|Y|}/(|r| + 2σ + rms)`.

---

## 0. Verdict

**Neither estimator is correct. Both are biased, in opposite directions, by the same per-wave systematic. The two "disagreeing"
remainders were also quoted with errors that were not errors.** [NE]

1. **The quoted errors were wrong.** The ±10 / ±108 (cv) and ±257 / ±131 (dir) on r are residual-scaled least-squares errors of a
   1-dof fit. Propagating the replica jackknife through the whole pipeline (a_h, bases, cv factor, fit) gives:
   * cv: r = −160 ± 213 / +327 ± 216;
   * dir: r = −1310 ± 526 / −810 ± 529;
   * paired same-seed difference: dir − cv = −1150 ± 479, i.e. **2.4σ, not ≫ 5σ**.
2. **The whole dir − cv gap is one node.** The linear response of r to each node's record predicts r_dir − r_cv = **−1140** from the
   scalar deviations; the measured value is −1150. Of that, **−1309 comes from ν = 13.3** (ℓθ₀ = 2.4), where
   ∂r/∂(relative record error) = 7.8·10⁴ c₂ (§2.2).
3. **The CV mean is exact.** The exact scalar ratio is correct to ≤ 10⁻¹⁰ (ρ₂ continuation), and ≤ 3·10⁻⁸ for the omitted reference
   ρ₂. So cv is not biased by a mis-specified mean (§2.1).
4. **The per-wave scalar bias is real, but seed 4242 overstated it at ν = 13.3.**
   * Four independent seeds give d = −0.0075 ± 0.0021 (3.6σ) at 13.3 and +0.0011 ± 0.0003 (3.3σ) at 16.3.
   * Seed 4242 alone had −0.0168 at 13.3, a 1.4σ upward fluctuation that cv imported into the graviton.
5. **The graviton per-wave bias is about half the scalar's.** Write λ = graviton relative wave bias / scalar relative wave bias (dir
   assumes λ = 0, cv assumes λ = 1). Three tests give:

   | test | λ |
   |---|---|
   | raw state-level test, 14 nodes, 4 states | 0.55 ± 0.29 |
   | paired comparison with the unbiased P_J twin on the same random numbers, pooled data | 0.58 ± 0.37 |
   | k = 3 cross-check, jackknife errors | 0.43 ± 0.50 |

   **Adopted: λ = 0.5 ± 0.25.**
   * dir carries the graviton's own wave bias.
   * cv replaces it by the scalar's, which is about twice as large, and so over-corrects by about the same amount.
6. **No heavy-tail pathology.** At ν = 13.3:
   * no point is dropped;
   * the Hill tail index is 2.6–3.2 (k = 50) and 3.6–5.6 (k = 200);
   * the top 10 points carry ≤ 3 % of the sum;
   * bootstrap/naive error = 1.00 ± 0.02.

   So the dir jackknife errors are not underestimated. The bias is a smooth systematic of the per-wave contours (HWIN's cut sweep),
   not missed tail mass.

**Corrected remainder and leverage at (1100, 100).** Pooled data, the λ-corrected per-wave estimator, error =
stat ⊕ 0.25·|r_cv − r_dir|. Worse helicity first, and the assumption-free P_J-cv twin in brackets:

| class (fit) | r (worse h) | L central | L at 2σ |
|---|---|---|---|
| A + C, 4 pw nodes | **−441 ± 180** `++` (P_J-cv −393 ± 174); `+−` +66 ± 181 | **2.3–2.8** (2.6–3.3) | **1.3** (1.5) |
| A + C, 6 nodes (+ fw 16.7, 33.4) | **+147 ± 77, rms 279** `++`; `+−` −179 ± 76 | **2.7–3.4** (3.1–3.8) | **2.0** (1.9) |
| A + B + C (full CHMRS ring class), 6 nodes | **−102 ± 122 / −100 ± 122, rms 58** | **7.2–9.0** (7.7–10.0) | **2.9** (2.7–2.8) |

* **Grade: [NE].**
* **H_LJ(1000) leverage at (1100, 100):**
  * **2.3–3.4 central** (A + C, worse helicity, depending on node set);
  * **1.3–2.0 at 2σ**;
  * **7–9 central, ≥ 2.7 at 2σ with the full ring class**.
* **This replaces both readings in N2_SUBREGGE:**
  * the cv headline: "3.0–4.2; ≥ 4.1 at 2σ";
  * the dir reading: "0.8–1.4".
* **A leverage ≤ 1 is excluded at ≥ 2σ in every fit:**
  * the worst case is A + C on 4 nodes, `++`: 1.3 (1.5 with P_J-cv);
  * dir's 0.8–1.4 came from the seed-4242 fluctuation at ν = 13.3 plus dir's own bias.
* **But the margin is modest.** With only the A + C class, the 2σ leverage is 1.3–2.0, not ≥ 3. The figure ≥ 2.7 needs the CHMRS
  B-class (ν²D) to be admitted as ring class (N2_SUBREGGE §3, R1).

---

## 1. What the two estimators compute (code)

* **Shared steps.** Both run the same pipeline: `pw_points` records → `pw_gamma.to_z2` → `fw_gamma.direct_records` →
  `z3_nodes.to_record` → `l9_analysis.grid` (q3/w2) → `fw_gamma.analyse` / `n2_refit` F1. They differ only in one per-node complex
  factor on the **state** record (`fw_gamma.py:83–89`).
* **State record.** At each window node the state record is the sum of two per-wave integrals, `W₊` and `W₋` (`pw_wave`: Legendre-Q
  waves, `W₊ + W₋ = P_J`). Each wave is integrated over the 7-dim crystal cycle, deformed to the n3_hybrid contour `C(ρ₁ = 3.5, ν_f)`
  of **its own** wave saddle (`pw_core.saddle`, homotopy in λ log W±).
  * Proposal: `n3_hybrid.adapt` on W₊ (seed 31, deterministic). It is reused for W₋, for the P_J twin and for the reference.
  * Sampling: 8 replicas × 2⁹ scrambled-Sobol points × 2 antithetic.
  * Normalisation: `logN = Pint(w_c) + log det L`, `− log β_{Δ,J}`.
  * Graviton: `pw_grav` (closed-form spin-J contraction, linear in f = W±) on the same contour points, with its own log replay and its
    own winding filter (`pw_core.grav_on_paths`).
* **Reference record.** Δ_ref = 2·10⁴, P_J on the n3 contour, same random numbers (CRN).
* **dir.** The plain CRN ratio, state/reference.
* **cv.** The graviton record is multiplied, as one complex scalar common to both helicities and both kinds (EE, C3H), by
  `exp[(lx_ex,st − lx_ex,ref) − log(S_MC,st/S_MC,ref)]` (`z3_nodes.scalar_check`). Here:
  * `S_MC` is the scalar integrand (the per-wave scalar for the state, P_J for the reference) on the same points;
  * `lx_ex` is the exact CLS-D scalar object (`n3_exact.X`: exact hybrid at ν ± i plus sm_exact's Arb ρ₂ continued by barycentric
    interpolation; heavy reference with ρ₂ := 0, `fw_scalar.install`).
  * It is a ratio control variate with coefficient fixed at 1.
* **Consequences.**
  * cv is unbiased iff the per-wave scalar is unbiased, or the graviton's relative wave bias equals the scalar's (λ = 1).
  * dir is unbiased iff the graviton per-wave estimate is unbiased (λ = 0).
  * The P_J twin (mode 0 on the n3 contour, `_wP`; P_J is entire, so no cut sweep) is unbiased under either estimator.
  * Its cv uses the P_J scalar, which is unbiased, so P_J-cv is a legitimate variance reduction.

Reproduction: `dc_runall.py` reproduces all six stored `pw_gamma.pkl` entries at (1100, 100) exactly, before any change.

## 2. Evidence

### 2.1 The CV mean is exactly known (`dc_cvmean.py` → `cache/dc_cvmean.txt`) [DV]
Measured as the relative change of the exact state/reference log ratio at the nodes 10–33.4:

| source of error in the CV mean | relative effect | in scalar Δc₂ |
|---|---|---|
| ρ₂ continuation (even/odd halves, NLOC 8/12) | ≤ 7·10⁻¹¹ | 0.00 |
| the whole ρ₂ term (set to 0) | ≤ 4·10⁻⁵ | ≤ 10.5 |
| reference ρ₂ := 0, bounded by the state ρ₂ scaled by (m_st/m_ref)² | ≤ 3·10⁻⁸ | ≤ 0.01 |

The scalar MC deviations being corrected are 10⁻³–10⁻². So a mis-specified mean is excluded.

### 2.2 Jackknife errors, and where dir − cv comes from (`dc_fit.py`, `dc_compare.py`, `dc_sens.py`) [NE]
* `dc_fit.fit` refits on each leave-one-replica-out record set. The a_h, the bases, the cv factors (`cv_jk`) and κ are recomputed
  each time (`dc_pipe`).
* Seed-4242 data, F1, 4 pw nodes:

  | estimator | r `++` | r `+−` | N2's quoted errors |
  |---|---|---|---|
  | cv | −160 ± 213 | +327 ± 216 | ±10 / ±108 |
  | dir | −1310 ± 526 | −810 ± 529 | ±257 / ±131 |
  | P_J-cv | −1193 ± 725 | −676 ± 727 | — |
  | P_J-dir | −597 ± 1084 | −86 ± 1093 | — |

* Paired (same replicas):

  | difference | `++` | `+−` |
  |---|---|---|
  | dir − cv | −1150 ± 479 | −1137 ± 475 |
  | cv − P_J-cv | +1033 ± 677 | — |
  | dir − P_J-cv | −117 ± 672 | — |

* **Linear response** (`dc_sens`). ∂r/∂ε_ν for a relative change (1 + ε) of the state record at one node, through the full
  pipeline:

  | node ν | 13.3 | 16.3 | 20.8 | 25 |
  |---|---|---|---|---|
  | ∂r/∂ε_ν (F1, 4 nodes) | +7.8·10⁴ | +4.4·10⁴ | +1.5·10⁵ | −2.7·10⁵ |

  Combined with the seed-4242 scalar deviations, the predicted dir − cv is −1140 (`++`) / −1128 (`+−`). Of this, **−1309 comes from
  ν = 13.3**; the other nodes give +88 (16.3), +5 (20.8), +76 (25).
* **The helicity independence is structural.** All four EE record components move together (`dc_raw`): the deviation is one common
  relative factor per node, so `++` and `+−` shift by the same amount.

### 2.3 The scalar per-wave bias, all seeds (`dc_scalseeds.py`) [NE]

| ν (ℓθ₀) | seeds (4242, `_s2`, 1001, 2002) | pooled d | χ²/dof | scalar Δc₂ |
|---|---|---|---|---|
| 13.3 (2.4) | −0.0168 ± 0.0067, −0.0114 ± 0.0057, −0.0093 ± 0.0041, −0.0038 ± 0.0029 | **−0.0075 ± 0.0021 (3.6σ)** | 4.2/3 | +144 ± 40 |
| 16.3 (3.0) | +0.0020 ± 0.0007, +0.0002 ± 0.0009, +0.0009 ± 0.0007, +0.0010 ± 0.0005 | **+0.0011 ± 0.0003 (3.3σ)** | 2.9/3 | −329 ± 101 |
| 25 (4.6) | two seeds | −0.00027 ± 0.00010 (2.7σ) | — | — |

HWIN's per-wave bias is confirmed. But at 13.3 it is less than half of what seed 4242 alone showed, and seed 4242 is the one the
graviton cv uses.

### 2.4 How much of the scalar bias the graviton shares: λ (`dc_mstest.py`, `dc_lambda.py`, `dc_lev.py`) [NE]
* **Raw state-level test.** The reference and the exact mean drop out.
  * `ms = log(G_w/S_w) − log(G_P/S_P)` (waves-cv vs P_J-cv): CV predicts E[ms] = 0.
  * `md = log(G_w/G_P) + [log(S_P/ref) − exact]` (waves-dir vs P_J-cv): DIR predicts E[md] = 0.
  * Data: 14 nodes at (1100, 100) (3 seeds), (1064, 64), (364, 64), (400, 100). The EE helicity trace gives χ²_CV = 6.0 and
    χ²_DIR = 6.4.
  * Fitting ms = (λ − 1)·d_w gives **λ = 0.55 ± 0.29** (all nodes), 0.57 ± 0.30 (window nodes ℓθ₀ ≤ 3.2) and 0.47 ± 0.34 (|d_w| > 2σ
    nodes).
  * C3H is too noisy to constrain λ (2.7 ± 1.3).
* **Γ-level paired test against the P_J-cv twin.** One-parameter family: the record times (exact/MC scalar)^λ. On the pooled data,
  `r(λ) − r(P_J-cv)` vanishes at λ = 0.59 ± 0.37 (F1, 4 nodes), 0.35 ± 1.4 (F1, 6 nodes) and 0.54 ± 0.72 (F1B).
* **k = 3 cross-check** (C3H shares a_h and has no de-shift, so its D_h coefficient must equal σ_h c0(C2)). With jackknife errors of
  the paired difference (`dc_k3.py`):
  * seed 4242, `+−`, 4 nodes: dir +110 ± 1543 (0.1σ), cv +5572 ± 2116 (2.6σ), P_J-cv −2826 ± 9214 (0.3σ), P_J-dir +252 ± 1832;
  * pooled: λ = 0.43 ± 0.50.
  * **N2's "P_J misses k = 3 by 1.5–2.1σ while dir/cv pass" was an artefact of the 1-dof lsq errors.** With jackknife errors P_J passes
    (0.1–0.3σ), and it is seed-4242 cv that misses (2.6σ on 4 nodes, 1.2σ on 6).
* **Adopted: λ = 0.5 ± 0.25** (rounded, conservative combination). λ = 1 (cv) and λ = 0 (dir) are each disfavoured at about 2σ.

### 2.5 Per-point diagnostics at ν = 13.3 (`dc_tail.py`, `dc_tailan.py`; 8192 points per mode, two seeds) [NE]

| mode | Hill α, scalar / EE (k = 50) | top 10 \|x\| share of Re Σ | bootstrap / naive error | graviton/scalar correlation | β_opt (CV) |
|---|---|---|---|---|---|
| W₊ | 3.2 / 2.7 | ≤ 1 % | 0.98–1.00 | 0.95 | 1.11 |
| W₋ | 2.8 / 2.6 | ≤ 3 % | 1.01–1.02 | 0.95 | 1.16 |
| P_J | 2.4 / 2.3 | ≤ 5 % | 1.00–1.03 | 0.95 | 1.00 |

* There are no dropped points, scalar or graviton (winding filter).
* Variance is finite and the errors are honest.
* The CV coefficient 1 is near-optimal **for variance**: about 10× variance reduction. This is why cv looked so precise.
* The point-wise graviton/scalar ratio varies by an IQR of 26–32 % on the wave contours, so a common relative bias is not expected.
  λ ≈ 0.5 is plausible for a cut-sweep contribution coming from a region where this ratio differs.

### 2.6 Pooled data (`dc_pool.py`: seeds 4242 + 1001 + 2002 at 13.3 and 16.3 → 8 super-replicas × 7; seed 4242 elsewhere)
F1, 4 nodes, `++` / `+−`:

| estimator | r `++` | r `+−` |
|---|---|---|
| dir | −693 ± 235 | −183 ± 236 |
| cv | −188 ± 63 | +316 ± 62 |
| λ = 0.5 | −441 ± 129 | +66 ± 131 |
| P_J-cv | −393 ± 174 | +106 ± 173 |

* The gap shrinks from 1150 to 505: most of the seed-4242 gap was the scalar fluctuation at 13.3.
* Paired against P_J-cv: dir −301 ± 177, cv +205 ± 196, λ = 0.5 −48 ± 148.

The full table, including the 6-node F1 and F1B fits, is in `cache/dc_lev_pool.txt` and `cache/dc_final_pool.txt`.

## 3. Error budget for the quoted r (pooled, λ = 0.5)

| item | size (F1, 4 nodes / F1, 6 nodes / F1B, 6 nodes) | treatment |
|---|---|---|
| statistics (jackknife, 8 super-replicas) | 129–131 / 71 / 102 | in σ. The jackknife σ itself is uncertain by ≈ ±25 % (8 groups). |
| λ = 0.5 ± 0.25 | 125 / 28 / 68 | in σ (0.25·\|r_cv − r_dir\|) |
| fit model (rms of residuals) | 12–68 / 106–279 / 58 | added linearly in L (N2 convention) |
| CV mean | ≤ 10 scalar Δc₂, i.e. ≲ 1 c₂ in r | negligible |
| heavy tails / point drops | none found | — |
| λ varying from node to node | covered by the assumption-free P_J-cv twin, which agrees within 0.3σ | cross-check column |
| per-row a_h | absorbed exactly by c0 in the raw-form fit (to ≲ 10 c₂, N2 §1) | — |

## 4. J ≥ 150 sample state: not computed
* No J ≥ 150 state at Δ_g = 1000 has K₂ basis rows or anchors (HWIN §5: "no basis rows for Δ = 1150").
* The c = 3 window (ν < 11.3 at J = 150, < 8.8 at 200) lies at or below the fold, where the per-node relative errors are ≥ 1 %. The
  sensitivity ∂r/∂ε of this section (≈ 10⁵ c₂) then gives σ_r ≳ 10³.
* Nothing "cheap" exists. H_LJ(1000) at J ≥ 150 stays [M].
* Add one caveat: the per-wave bias found here grows as ℓθ₀ falls (§2.3), and every J ≥ 150 window has ℓθ₀ ≲ 2.

## 5. Corrections to existing documents (not applied)
1. **N2_SUBREGGE.md §0** (rows "Corrected leverage", "caveat (not N2): estimator", "Task 2: graviton window data"), **§2.4** (table
   and bullets) and **§3.4**:
   * Replace "L = 3.0–4.2 (cv, A + C, central; 1.8 at 2σ); ≥ 4.1 at 2σ with the full class" and "dir: L ≈ 0.8–1.4; unresolved" by:
     "**L = 2.3–3.4 central (A + C, worse helicity), 1.3–2.0 at 2σ; A + B + C 7–9 central, ≥ 2.7 at 2σ** (λ-corrected per-wave and
     P_J-cv twin, pooled seeds; DIRCV)".
   * The r errors quoted there (±10, ±83, ±108, ±257, ±131) are 1-dof residual-scaled lsq errors. The jackknife errors are ±145–216
     (cv) and ±526 (dir) on the seed-4242 data.
   * N2-e's C3H numbers carry the same caveat.
2. **Update boxes** (UNIFIED_THEOREM.md:282–288, HWIN_THIMBLE.md:11–17, HLJ_CLOSE.md:25–31, LARGESPIN_ANALYTIC.md:18–24,
   POS_ERRORS.md:38–44):
   * Replace "3.0–4.2 (cv), ≥ 4.1 at 2σ; 0.8–1.4 with dir … unresolved and dominant" by "**2.3–3.4 central, 1.3–2.0 at 2σ (A + C);
     ≥ 2.7 at 2σ with the full ring class. The dir/cv systematic is resolved: both are biased, the graviton carries λ = 0.5 ± 0.25
     of the scalar per-wave bias, and the unbiased P_J-cv twin agrees [NE]**".
   * The dominant uncertainty is now statistics at ν ≈ 13 plus λ.
3. **HWIN_THIMBLE §0.2 / §3.1** ("the graviton uses the scalar control variate; |dir − CV| is carried as a systematic"):
   * The CV over-corrects. The graviton's per-wave bias is λ = 0.5 ± 0.25 times the scalar's.
   * The §4 window values (cv) carry a systematic of ≈ 0.5·|Γ_dir − Γ_CV|. The unbiased P_J-cv values are the reference.
   * §3.2 (1100, 100): pooled scalar bias −0.0075 ± 0.0021 (13.3) and +0.0011 ± 0.0003 (16.3), with four seeds each.
4. **POS_ERRORS / HLJ_CLOSE "quote 1.1–4.8"** → "2.3–3.4 (1.3–2.0 at 2σ), A + C".

## 6. What would tighten this (value order)
1. **More P_J-twin statistics at ν ≈ 13 (×4).** This brings σ_r to ≈ 90 with **no λ assumption**, because P_J-cv is bias-free. About
   1.5 CPU-h with `dc_points.py`.
2. **Put both waves on one contour inside the window (ℓθ₀ ≲ 3)**, which removes the cut sweep by construction (HWIN §7.5). Or measure
   λ per node.
3. **More window nodes (ν ≈ 11–15) to increase the dof of the F1/F1B fits.** At ν = 13.3, r currently responds with 7.8·10⁴ c₂ per
   unit relative error.

## Hand-over
* **Code:** `ads4/code/dircv/`.
  * Analysis: `dc_pipe.py` (Γ − κ with jackknife replicas), `dc_fit.py` (F1/F1B with jackknife), `dc_compare.py`, `dc_sens.py`,
    `dc_k3.py`, `dc_cvmean.py`, `dc_raw.py`, `dc_reps.py`, `dc_scalseeds.py`, `dc_mstest.py`, `dc_lambda.py`, `dc_lev.py`,
    `dc_final.py`.
  * MC: `dc_points.py` (in-process `pw_points`, output redirected), `dc_pool.py`, `dc_tail.py`, `dc_tailan.py`.
* **Caches:** `cache/pw_pts_J100_D1100_x{1,2}{,P}.pkl` (new seeds) and `_pool{,P}.pkl`; per-point samples in `tail_*.pkl`; the text
  results named above.
* **Nothing is running.**

---

## 7. Addendum: 4× P_J-twin statistics at ν = 13.3 (the λ-free route) [NE]

*Code:*
* `dc_ptwin.py` (P_J twin only: state + reference, same contour and proposal as production);
* `run_ptwin.sh`: seeds 5001 and 6001, 84 replicas each, for +168 replicas;
* `dc_pool.py` now takes separate P_J tags and an output tag; `DC_INTERLEAVE=1` sets the grouping;
* `dc_pjseed.py`.

*Compute:* +1.65 CPU-h, two processes. The total for DIRCV is ≈ 3.9 CPU-h.

**Data.** The P_J twin at ν = 13.3 now pools 224 replicas (8 + 24 + 24 + 84 + 84), 4× the 56 of §2.6. The waves are unchanged. The
pool is `_pool2`, contiguous groups of 28; `_pool2i` is the same pool grouped by interleaving.

**Result: the error did not drop to ≈ 90.** The P_J twin at ν = 13.3 has a heavy tail at the replica level.
* The seeds agree with each other: the state quantity log(G_P/S_P) has χ² = 3.6/4 over the 5 seeds (`dc_pjseed.txt`). So there is no
  sign of bias.
* But one replica of seed 6001 has a leave-one-out pseudo-value ≈ 8× the typical one. That run also had 6 dropped points and maximum
  winding 3.07. Seed 6001 alone has twice the error of seed 5001 at equal size.
* This matches the P_J per-point Hill index α ≈ 2.3–2.4 (§2.5): the variance is barely finite.
* The jackknife error of r therefore depends on the grouping (contiguous vs interleaved: 306 vs 380 on 4 nodes). The larger one
  (interleaved: seeds spread evenly, so the groups are closer to i.i.d.) is quoted below.
* **The earlier 56-replica P_J-cv values (§0 brackets, §2.6) were not robust.** For example, F1 4-node `++` moved from −393 ± 174 to
  +23 ± 380.

**λ-free (P_J-cv) remainder and leverage at (1100, 100).** Worse helicity in bold.

| class (fit) | r `++` (rms) | r `+−` (rms) | L central | L at 2σ |
|---|---|---|---|---|
| A + C, 4 pw nodes | +23 ± 380 (39) | **+534 ± 387 (96)** | **1.8–2.3** | **0.8** (`++`: 1.4) |
| A + C, 6 nodes | **+289 ± 124 (133)** | −35 ± 124 (250) | **2.7–3.4** | **1.7** (`+−`: 2.2) |
| A + B + C, 6 nodes | +181 ± 231 (60) | +184 ± 232 (59) | **4.7–6.0** | **1.6** |

(With the contiguous grouping the r errors are 306 / 313, 87 / 91 and 171 / 175. The 2σ leverages are then 1.7 / 0.9, 1.9 / 2.5 and
2.0 / 1.9.)

**Effect on λ.** The P_J-cv centre moved towards cv.
* Compared without pairing with the per-wave family r(λ), it now prefers λ ≈ 1.4–1.7 ± 0.7.
* The raw state-level test is unchanged (0.55 ± 0.29; the waves were not rerun), and k = 3 gives 0.43 ± 0.50.
* Combined: **λ ≈ 0.6 ± 0.25**. The λ-corrected table of §0 (λ = 0.5 ± 0.25) moves by less than its errors (e.g. A + C 4 nodes `++`:
  r ≈ −390 ± 200, L = 2.5, 1.35 at 2σ).

**Net (supersedes the bracketed P_J-cv numbers in §0):**
* **λ-corrected per-wave estimator** (§0, unchanged): A + C 2.3–3.4 central, 1.3–2.0 at 2σ; A + B + C 7–9 central, 2.9 at 2σ.
* **λ-free P_J-cv (this section):** A + C 1.8–3.4 central, **0.8–1.7 at 2σ**; A + B + C 4.7–6.0 central, **1.6 at 2σ**.
* **Correction to §0's last bullet.** "A leverage ≤ 1 is excluded at ≥ 2σ in every fit" holds only for the λ-corrected estimator.
  Without any λ assumption, the 4-node A + C fit reaches **0.8 at 2σ** (`+−`). The 6-node A + C fit and the full ring class stay at
  ≥ 1.6 at 2σ.
* **Central leverages are > 1.8 in every variant.**
* **Next step.** More P_J statistics at ν = 13.3 do not help efficiently, because of the heavy tail. Either:
  * a variance-bounded twin at 13.3 (both waves on one contour, or a tail-robust proposal); or
  * more window nodes (ν ≈ 11–15) to reduce the 4-node fit's leverage on ν = 13.3.
