# EXT_FOLD: the external direct points (jobs A, B) folded into the Δ_g = 1000 near-gap positivity budgets

*23 Nov 2026. Unreviewed, not committed by this agent. Another agent's commit 182cb14 swept in a snapshot of `code/extfold/`, identical to
the working tree. No existing file was edited (the external pickles, `code/highspin/`, `code/smallnu/`,
`code/neargap2/` and `code/spinfunc/` and their caches were only read); `git status` shows only this file.*
* *New code:* `ads4/code/extfold/`. `ef_guard.py` is the process-local write guard: every write outside `extfold/cache/` raises. It also
  registers itself as `ng_guard`, `hs_guard` and `sn_guard`. The scripts are:
  * `ef_core.py`: merged point reader, J = 7 shadow registration, output redirection;
  * `ef_merge.py`: merges the point sets and pools the J = 29 streams;
  * `ef_shadow.py`: the odd ⟨TTO⟩ shadow of (1007, 7);
  * `ef_points.py`: local re-computation of external points;
  * `ef_oddrun.py`: `highspin/hs_odd.py` run unchanged;
  * `ef_window.py`: `smallnu/sn_odd.py` run unchanged;
  * `ef_odddiag.py`: representation and J = 31 edge diagnostics;
  * `ef_j7.py`: `neargap2/ng_bind.run_state` with the own direct row, the `hs_even` protocol for odd J;
  * `ef_j7nodes.py`;
  * `ef_streams.py`: the J = 29 stream-by-stream analysis;
  * `ef_j29env.py`: the J = 29 envelope decomposition and what-ifs;
  * `ef_summary.py`.
* *Caches:* `ads4/code/extfold/cache/`. Its `.gitignore` ignores everything but itself. The merged point sets are in `cache/pts/{old,base,all}/`.
* *Compute:* ≈ 0.25 CPU-h (two local re-computations of external points ≈ 9 min; all analyses ≈ 4 min), at most two processes,
  `ulimit -v 3000000`, one thread, no bytecode written.

Grades follow UNIFIED_THEOREM_v2 §0, POS_ERRORS, NEARGAP2, HIGHSPIN and SMALLNU_WINDOW: [P], [certified], [D], [DV], [NE±], [NE],
[NE + M], [M], [I]. `L = value / budget`. A "×3 envelope" is SMALLNU's sign-free measured window envelope (node bound |Γ| + 3σ, ×3, and
10⁴ c₂ below the lowest node). "Signed" means the window is integrated from the direct values. The odd-J λ-level budget is HIGHSPIN's,
×3 on every model term and 3σ on MC.

---

## 0. Verdict

| state (Δ_g = 1000) | value `λ/N_s` (old → new) | budget (new) | `L` old → new | grade old → new |
|---|---|---|---|---|
| **J = 7, twist edge `1.00698Λ`** | published K₂ chain (rows 800/1600, 1/m-interpolated) `+0.26822` → **own direct row `+0.26832 ± 0.00004`** | `0.040` sign-free (edge x > 0.775 at 10⁴: 0.038; ν < 50 at 10⁴: 0.0024; truncation 0; MC 0.0001) | — (no budget existed) → **6.7** (λ-level 6.7) | [NE + M] → **[NE±]** |
| J = 7, 1.01 / 1.02 / 1.05Λ | `+0.3115 / +0.5832 / +1.0000` | `0.044 / 0.056 / 0.040` | → **7.2 / 10.3 / 25** | [NE + M] → **[NE±]** |
| **J = 27, twist edge `1.02667Λ`** | K₂ bulk + GK15 edge `+0.959` → **direct bulk + GK15 edge `+0.976 ± 0.019`** (λ_c `+0.998`) | `0.347` (MC 0.058, quadrature 0.064, missing GK 0.136, representation 0.089; truncation 0) + window ×3 `0.263` = `0.610` | `0.28` → **1.64** (signed 2.02) | [NE + M] → **[NE±]** |
| J = 27, 1.03 / 1.04 / 1.05 / 1.075 / 1.10Λ | `+0.977 / +0.973 / +0.944 / +0.980 / +0.957` | `0.518 / 0.398 / 0.364 / 0.349 / 0.303` | `— / 0.31 / 0.32 / 0.34 / 0.37` → **1.93 / 2.50 / 2.65 / 2.86 / 3.21** | [NE + M] → **[NE±]** |
| **J = 23, twist edge `1.02276Λ`** | `+0.883` → **`+0.923 ± 0.016`** (λ_c `+0.901`) | `0.925` (MC 0.047, quadrature 0.023, missing GK 0.303, **representation 0.553**) + window ×3 `0.131` = `1.056` | `0.45` → **0.85** (signed 0.90; noise-aware representation 1.52) | [NE + M] → **[NE + M]**. The blocker is no longer K₂ truncation: it is the degree-7 shape fit of the noisy direct row (§3.2) and the failed GK node x = 0.941 |
| J = 23, 1.03 / 1.05 / 1.075 / 1.10Λ | `+0.945 / +1.021 / +0.970 / +1.018` | `0.793 / 0.590 / 0.455 / 0.429` | `0.51 / 0.57 / 0.61 / 0.64` → **1.16 / 1.69 / 2.09 / 2.33** | [NE + M] → **[NE±]** |
| **J = 31, twist edge `1.03056Λ`** | hs_odd spline edge `+0.786` → `+0.805 ± 0.023`; with the unmeasured edge tail removed `+0.862` (λ_c `+0.881`); published shape continuation `+1.02` | `hs_odd` protocol `0.479 + 0.427` (L 0.91, **understated**, §3.3); corrected: MC 0.069, edge interpolation 0.147, **edge tail x > 0.9 (sign-free) 1.846**, representation 0.377, window 0.427 = `2.87` | `0.08` → **0.31** (corrected) | [NE + M] **unchanged**. J = 31 has **no GK15 edge run**, only HC's 8 equispaced nodes x = 0.8–0.9 |
| J = 31, 1.04 / 1.05Λ | `+0.843 / +0.897` (tail removed `+0.870 / +0.909`) | corrected `1.67 / 1.11` | `0.10 / 0.11` → `0.53 / 0.83` (signed 0.61 / 1.02) | [NE + M] unchanged |
| J = 31, 1.075 / 1.10Λ | `+0.926 / +0.964` | corrected `0.618 / 0.511` | `0.16 / —` → **1.53 / 1.91** | [NE + M] → **[NE±]** |
| **J = 29, twist edge `1.02862Λ`** | `+0.972` (unchanged: the bulk and edge are HIGHSPIN's) | signed `0.789 → 0.784` (window MC 0.041 → 0.035) | signed **1.23 → 1.24**; **sign-free ×3: 0.94 → 0.94** | [NE±] (thin, signed) unchanged; **the sign-free budget still does not close** (§4) |
| J = 29, 1.04 – 1.10Λ | unchanged | — | signed `2.24…5.70 → 2.27…5.86`; ×3 `1.47…2.81` (unchanged) | [NE±] unchanged |

(J = 31 at 1.03Λ lies below its twist edge `1.0306Λ`, so it is not an allowed state; it was computed only as a row value.)

**Negative findings (reported prominently).**
1. **No state is negative.** Every value computed is `≥ +0.268`; the lowest is J = 7 at its twist edge, as expected next to the
   J = 5 binding state, `+0.20`. The lowest at large J is J = 31's twist edge, `+0.805` with the spline edge.
2. **J = 29 does not close sign-free.** The ×3 envelope stays at `L = 0.94`. Three things block it:
   * **The 16 "independent seeds" contain only 23 independent RQMC streams, not 128.** `z2_points.n3_block` draws replica r of a job with
     `Z2_SEED = s` from stream `s + r` with a deterministic proposal, so consecutive seeds 5001…5016 share 7 of their 8 replicas. This was
     checked bit for bit on all 840 shared replicas.
   * **The n3 streams are heavy-tailed.** Single streams sit at `|z| = 11.4` (ν = 12.5, C3H) and `11.7` (ν = 15, C3H), so pooling all 31
     distinct streams does not shrink the jackknife σ there (C3H at ν = 12.5: `398 → 661`).
   * **The binding node is not in the B window.** Even with σ → 0 on all of ν = 10–17.5, the ×3 envelope gives only `L = 0.97`. The
     n1 node at ν = 20 (C3H `−761 ± 727`, bound 2943 c₂) and the 10⁴ convention below ν = 10 set it. Closing needs that node resolved:
     `L = 1.05` if it were known to its central value.

   SMALLNU's "≈ 16× the n3 points would close it" is therefore not right.
3. **`hs_odd`'s non-GK edge budget has a gap.** Without GK nodes, the spline edge holds the last node's value (J = 31: `Γ(0.9) ≈ 1.4·10⁴ c₂`)
   over `(0.9, Λ/m]`. Its "quadrature" term compares this only with another held continuation (`0.032`). Covering that tail sign-free
   with HIGHSPIN's own light-cone envelope costs `1.85 N_s` at the twist edge. J = 31's protocol `L = 0.91` thus becomes `0.31`.
   SMALLNU's J = 31 budget (9.9, not closing anyway) also omits this.
4. **E_repr, the degree-7 probe of the representation term, is MC-dominated on noisy direct rows.** At J = 23 its jackknife σ is 5× that
   of the degree-6 main fit, and degree 8 runs away (`+2.0`). It alone blocks J = 23's twist edge. SMALLNU's J = 21 (`+0.815` at degree 7)
   is the same effect. This is a protocol issue, not a negative state (§3.2).
5. **Two smaller items:**
   * **Job A's "seed 4242" is irrelevant.** The n1 contour always uses Sobol seeds `2000 + r`; `Z2_SEED` acts on n3 only. The external
     n1 points are deterministic and reproduce locally to `1.5·10⁻¹²`.
   * **Job C (J = 400) failed on the external machine** (float SVD in the conserved-structure solver, per the commit message). There is no
     data for the 345 ≤ J < 857 band.

---

## 1. Inputs

* **Job A** (`ext/unpacked/out/A_bulk_J{7,23,27,31}/`): 34 n1 points per state, x = 0.05 … 0.775 on the own row Δ = 1000 + J, state and
  Δ_ref = 2·10⁴, 8 replicas. They use the settings of HIGHSPIN's and SMALLNU's bulk rows (`HC_ADAPT = 1`, `HC_M = HC_MREF = 11`). No point
  failed: the replica error is < 20 % everywhere, and `hs_core.failed_x` is empty on the bulk.
* **Job B** (`ext/unpacked/out/B_n3_J29_s50xx/`): n3 at ν = 10, 12.5, 15, 17.5 on (1029, 29), 16 jobs `Z2_SEED = 5001 … 5016`, each with
  8 replicas, so streams `s … s + 7`.
* **Existing direct points merged** (read-only; `ef_merge`, `cache/ef_merge.txt`):
  * J = 23, 27: HIGHSPIN's GK15 edge nodes (11 each, `[0.8, Λ/m]` capped at 0.945; J = 23's node x = 0.941 failed) and SMALLNU's 5 n1
    nodes at ν = 20–47.
  * **J = 31: no GK15 edge run exists.** Its edge has only HC's 8 equispaced n1 nodes x = 0.8 … 0.9 (`spinfunc/cache`, identical copy
    in `highspin/cache`), plus SMALLNU's 5 small-ν n1 nodes.
  * J = 7: no direct point existed, and **no K₂ row at Δ = 1007 exists either**: the chain's J = 7 rows are 200/400/800/1600.
  * J = 29: HIGHSPIN's 54 (bulk + GK15 + 0.8–0.9), NEARGAP2's 5 edge nodes (identical copies) and SMALLNU's 10 window nodes.
* **Collisions:** between old sources only, all bit-identical. The external points share no (Δ, ν) key with any existing direct point.

## 2. Methods (`code/extfold/`)

* **Merged reader.** `ef_core` replaces `hc_nodes.pts` (and `sn_core._pts`, `sn_nus`, `rel_err`) with a reader of the merged sets:
  * `all`: every direct point;
  * `base`: highspin/cache only, which is sn_core's baseline normalisation (`BASE[0]`);
  * `old`: without the external points.

  Everything downstream is the unchanged project code: `hs_core` → `neargap2/ng_val` → spinfunc chain, `hs_odd`, `sn_odd`, `sn_window`,
  `ng_bind`.
* **Reproduction of the old numbers.** With set `old` the wrapper reproduces SMALLNU's `sn_odd` results for J = 23, 27, 29, 31 to
  every printed digit: J = 29 `L` ×3 / signed `0.94 / 1.23`; J = 23 / 27 / 31 `0.451 / 0.282 / 0.081`. `hs_odd` on the merged J = 23,
  27, 31 rows reproduces HIGHSPIN's K₂-bulk values (`K2/gk +0.8830 / +0.9585`, `K2/espl +0.7858` vs `+0.8831 / +0.9587 / +0.7855`) to
  ≤ 3·10⁻⁴. The residual comes from the small-ν nodes moving the row normalisation slightly.
* **Odd J = 23, 27, 31** (`ef_oddrun` = SMALLNU's `sn_oddrun` for J = 21): `hs_odd.py` unchanged. With 34/34 direct bulk nodes the main
  bulk is `dir` (shape fit degree 6). The edge is GK15 quadrature of the direct GK nodes (J = 23, 27) or the spline through the HC nodes
  (J = 31). The truncation term is 0 because the bulk is direct. The window comes from `ef_window` (`sn_odd.main` unchanged), which uses
  the **FULL normalisation** (the grid with all direct points) exactly as SMALLNU did for its J = 21 direct-bulk row (SMALLNU §4.2: the
  window must be quoted in the normalisation of the value it is combined with).
* **J = 7** (`ef_j7`): HIGHSPIN's even-J own-row protocol (`hs_even` → `ng_bind.run_state`, unchanged), applied to the odd 1×1 state.
  * The direct row 1007 is injected as an own-mass row into the chain's 1/m model.
  * The budget is NEARGAP2 §1.2 sign-free: truncation, interpolation, `ν < 50` at 10⁴ c₂ with the gauge factor, and the edge `x > 0.775`
    at `G_edge = max(10⁴, 10 max|Γ|)`. On top of that, `E_mc` = 3 × the jackknife σ.
  * Turning direct points into chain records needs the parity-odd ⟨TTO⟩ shadow of (1007, 7), which did not exist. `ef_shadow` computed it
    with `p2_shadow_odd.shadow_all` at its defaults (ξ = 0.5, 8 triples, seed 7), and first reproduced the cached shadows of (800, 7) and
    (1023, 23) **bit for bit**.
* **J = 29 pooling** (`ef_merge.pool_n3`): every distinct stream (4242–4249 from SMALLNU, 5001–5023 from the external jobs; 31 in all, with
  identical proposals and `logN` agreeing to `≤ 10⁻⁹`) is put into one record. The record has 8 super-replicas (stream groups 4,4,4,4,4,4,4,3,
  state and reference grouped identically), so the chain's 8-replica jackknife runs unchanged. `ef_streams` also evaluates Γ stream by
  stream through `sn_core.small_nodes`.
* **Diagnostics** (`ef_odddiag`):
  * the degree-7 fit's own jackknife;
  * a model-free bulk (natural spline through the direct nodes);
  * a "noise-aware" representation term: ×3 the largest deviation of degrees 4 and 5, the two half-samples and the spline bulk,
    without the degree-7 fit;
  * for J = 31, the edge with its unmeasured tail removed and covered sign-free.

## 3. Results, odd J

### 3.1 J = 27: closes (`cache/hs_odd_J27.txt`, `sn_odd_J27.txt`)

* **Direct vs K₂:** in the bulk, `|Δ|` median 64 c₂, max 462 (x = 0.75), median 1.7σ. The K₂ row overestimates the oscillation amplitude
  at x ≳ 0.6, as at J = 23 and 31. The twist edge moves `+0.959 → +0.976` (+0.018, favourable).
* **Representation is stable:** degree 7 gives `+0.977` (0.0σ from degree 6), and the spline bulk `+0.993`.
* **Window:** measured in the direct-bulk normalisation, C2 `+120 … +198`, C3H `−165 … −215` on ν = 25–47; ν = 20: `−554 ± 597` (C3H).
* **Result:** twist edge `L = 1.64` (×3), `2.02` (signed), `1.09` (×10). The 10⁴ convention would give `0.94`. All masses to 1.10Λ have
  `L ≥ 1.64`. **[NE±].**

### 3.2 J = 23: masses ≥ 1.03Λ close; the twist edge is representation-limited (`cache/hs_odd_J23.txt`, `ef_odddiag_J23.txt`)

* **Direct vs K₂:** median 58 c₂, max 385, median 1.6σ. The twist edge moves `+0.883 → +0.923` (+0.040, favourable).
* **Budget at the twist edge** (λ-level): MC 0.047, GK15 − G7 0.023. The missing-node term 0.303 comes from x = 0.941 (failed) and the
  nodes above the 0.945 cap. The **representation term is 0.553**: the degree-7 shape fit gives `+1.107` against `+0.923`.
* **Diagnosis of the representation term.** The degree-7 fit's jackknife σ is `0.077` (degree 6: `0.016`), and its deviation is 2.4σ of
  that. Degree 8 gives `+2.00`. Degrees 4 and 5, both half-samples and the model-free spline bulk give `+0.893 … +0.925`. So the degree-7
  probe mostly amplifies the MC noise of the 34 direct nodes (direct σ ≈ 36 c₂, against a fit with 24 parameters). It is not clean
  evidence of a representation ambiguity, but 2.4σ is not nothing.
* **Result:** by the protocol, `L = 0.85` (×3) / `0.90` (signed) at the twist edge, so [NE + M]. With the noise-aware representation
  term (`0.089`), `L = 1.52`, but that is a changed protocol, so it is graded [NE + M], not [NE±].
* **Masses 1.03–1.10Λ:** `L = 1.16 / 1.69 / 2.09 / 2.33` by the protocol. **[NE±].**

### 3.3 J = 31: no GK15 edge; the edge tail blocks the low masses (`cache/hs_odd_J31.txt`, `ef_odddiag_J31.txt`)

* **Reference zero and bulk.** The reference zero sits at x = 0.7477, and `hs_core.refbad` drops the K₂ node 0.75. The direct node at
  x = 0.75 (`|x − x₀| = 0.0023`) has normal errors (±65). Direct vs K₂: median 89 c₂, max 962 (x = 0.775, next to the reference zero),
  median 2.0σ.
* **Spline edge value.** With HC's equispaced edge nodes (spacing 0.0125, against an edge oscillation period ≈ 0.05, with `Γ` growing
  from 1.7·10³ at 0.8 to 1.4·10⁴ at 0.9), `hs_odd` takes the spline edge: `+0.805 ± 0.023` (K₂ bulk `+0.786`).
* **The spline edge holds `Γ(0.9) = ±1.4·10⁴` constant out to Λ/m = 0.970**, and `hs_odd`'s E_quad for non-GK edges (`3|espl − held|`)
  compares it only with another held continuation. With that unmeasured tail set to 0 the twist edge is `+0.862`.
* **Corrected budget.** The tail is covered sign-free with HIGHSPIN's E_miss envelope `3A((1 − x_l)/(1 − x))²` (A = 1.48·10⁴ c₂,
  x_l = 0.9), and the [0.8, 0.9] spline is checked against linear interpolation and the two every-other-node subsamples (×3):

  | m/Λ | tail cost | edge interpolation | corrected `L` (×3 / signed) |
  |---|---|---|---|
  | 1.0306 (twist) | 1.85 | 0.15 | 0.31 / 0.33 |
  | 1.04 | 0.80 | 0.09 | 0.53 / 0.61 |
  | 1.05 | 0.32 | 0.06 | 0.83 / 1.02 |
  | 1.075 | 0.014 | 0.012 | **1.53 / 2.15** |
  | 1.10 | 0 | 0.001 | **1.91 / 2.83** |
* **Representation:** 0.12–0.38. The half-samples give `+0.93` against `+0.805` at the twist edge; the edge representation and the bulk
  fit interact here.
* **Result:** the twist edge to 1.05Λ stays **[NE + M]**; 1.075 and 1.10Λ are **[NE±]**. **What would close it:** the GK15 edge window
  `[0.8, 0.9703]`, 15 n1 nodes as HIGHSPIN ran for J = 23, 27, 29. Expect ≈ 4 nodes above x ≈ 0.94 to fail; HIGHSPIN's E_miss covers
  those at a cost of ≈ 0.1–0.3.

### 3.4 J = 7 (`cache/ef_j7.txt`, `ef_j7nodes.txt`)

* **Direct row:** `|Γ| ≤ 212 c₂` everywhere, with σ = 4–25 c₂. Against the published model of the same row (1/m interpolation of the
  K₂ rows 800 and 1600, spline in x): `|Δ|` median 14, max 53 c₂ (median 1.2σ, max 5.5σ, C3H at small x). **H_M is confirmed at J = 7.**
* **Value:** the state moves `+0.26822 → +0.26832` (twist edge) and `+0.58304 → +0.58318` (1.02Λ).
* **Budget:** sign-free `0.040` at the twist edge, almost entirely the edge `x > 0.775` held at 10⁴ c₂ (×47 above the largest measured
  `|Γ|`). The truncation term is 0 on the own row, interpolation is `9·10⁻⁶`, MC `1.2·10⁻⁴`. **`L = 6.7`** sign-free (λ-level 6.7);
  `L = 7.2 / 10.3 / 25` at 1.01 / 1.02 / 1.05Λ.
* **Certified-floor route:** `(N_f/5 − Lemma C)/N_s = +0.1988`, so `λ/N_s ≥ +0.159` with no hybrid-quadrature input. **[NE±].**
* For comparison, the published rows alone with NEARGAP2's budget already give `L = 4.3`. J = 7 was [NE + M] only because no own-mass row
  had been built and the 1/m interpolation itself was unchecked.

### 3.5 K₂ calibration, extended (HIGHSPIN §4)

The table compares the actual K₂ error, `λ_direct − λ_K2`, with the last Laplace step, `λ_K2 − λ_K1`, at the twist edge, both with the
same edge treatment:

| J | 12 | 16 | 21 | 23 | 27 | 29 | 31 |
|---|---|---|---|---|---|---|---|
| actual / step | 0.003 / 0.069 | 0.008 / 0.160 | 0.018 / 0.40 | 0.040 / 0.493 | 0.018 / 0.984 | 0.046 / 1.28 | 0.020 / 3.16 |
| ratio | 4 % | 5 % | 4.5 % | **8 %** | 1.8 % | 3.6 % | 0.6 % |

The ratio is **0.6–8 % of the step over 7 states, always with direct > K₂ (favourable)**. SMALLNU's "calibrated" truncation (5 % of the
step, ×3 = 15 %) covers every case, including J = 23's 8 %. HIGHSPIN's "4–5 %" should read "0.6–8 %".

## 4. J = 29: the pooled sign-free window (`cache/ef_streams.txt`, `ef_j29env.txt`, `sn_odd_J29.txt`)

**Streams.** There are 31 distinct n3 streams per ν. Γ − κ in c₂, SMALLNU's baseline normalisation:

| ν | kind | SMALLNU (8 streams) | pooled (31, super-replica jackknife) | 31-stream mean ± sd/√31 | max single-stream `|z|` | block χ² (O, A, B, C; 3 dof) |
|---|---|---|---|---|---|---|
| 10 | C2 / C3H | `−681 ± 634` / `+1113 ± 528` | `−437 ± 344` / `+274 ± 343` | `−482 ± 462` / `+261 ± 327` | 1.9 / 4.1 | 1.19 (p 0.76) / 4.27 (p 0.23) |
| 12.5 | C2 / C3H | `−173 ± 338` / `−87 ± 398` | `−69 ± 155` / `−265 ± 661` | `−86 ± 208` / `−255 ± 593` | 3.5 / **11.4** | 0.51 (0.92) / 1.41 (0.70) |
| 15 | C2 / C3H | `+50 ± 361` / `−235 ± 381` | `−73 ± 120` / `+209 ± 431` | `−80 ± 168` / `+439 ± 644` | 1.6 / **11.7** | 0.66 (0.88) / 2.48 (0.48) |
| 17.5 | C2 / C3H | `+66 ± 412` / `−393 ± 403` | `−64 ± 132` / `+12 ± 272` | `−74 ± 161` / `−4 ± 246` | 1.6 / 2.9 | 1.57 (0.67) / 2.49 (0.48) |

The independent blocks are O = 4242–4249 (SMALLNU), A = 5001–5008, B = 5009–5016 and C = 5017–5023. **Every block χ² is consistent.**
**The 16 jobs as submitted are not independent:** their per-job means drift smoothly from job to job (e.g. ν = 12.5 C3H
`+794 … +1652` for jobs 5001–5006, then `−682 … −2628`), because neighbouring jobs share 7 of 8 streams. A χ² over the 16 jobs is
therefore meaningless, and was not formed. Every value is within 1.3σ of 0. C2 shrinks by ×1.8–3.1 as expected. **C3H at ν = 12.5 and 15 does
not shrink**: one stream each (5014, 5017) sits at `z ≈ 11`, the heavy tail already seen at (1100, 100) (DIRCV §7).

**Envelope** at the twist edge. The window may cost at most `λ_c − 0.697 = 0.268` for L = 1:

| scenario | ×3 envelope | `L` |
|---|---|---|
| S0 SMALLNU | 0.327 | 0.94 |
| S1 pooled (main) | 0.330 | **0.94** |
| S2 pooled, σ = sd/√31 | 0.333 | 0.94 |
| S3 ν = 10–17.5 at `|Γ|` (σ → 0) | 0.293 | 0.97 |
| S4 = S3 + the ν = 20 n1 node at `|Γ|` | 0.223 | 1.05 |
| S5 pooled, ν = 20 node dropped | 0.261 | 1.01 |

* The per-segment cost of S1 is 0.014 / 0.017 / 0.012 on 10–17.5, **0.029 + 0.069 on 17.5–25 (the ν = 20 node)**, 0.136 on 25–51, and
  0.052 below ν = 10.
* A uniform node bound `|Γ| + 3σ ≤ 1.4·10³ c₂` on [10, 25) would close it; the ν = 20 node is at 2.9·10³.
* The **signed window** improves slightly: its MC/interpolation error goes `0.041 → 0.035`, and `L` goes **1.23 → 1.24** (1.04–1.10Λ:
  `2.27 … 5.86`).

**Verdict J = 29.** **The twist-edge budget does not close sign-free:** `L = 0.94 → 0.94`. J = 29 stays [NE±] on the signed window
(thin, `L = 1.24`), and sign-free at 1.04–1.10Λ (`L = 1.47 … 2.78`), as before.
**What would close it:**
* n3 at ν = 20 and 22.5, to replace the near-fold n1 node;
* ≥ 64 **distinct** streams at ν = 12.5–20, with `Z2_SEED` spaced by ≥ 8 (e.g. 6000, 6008, …) and a heavy-tail-robust error;
* the ν < 10 convention stays (0.052).

Even with all of that the margin is thin (S4: 1.05).

## 5. Consistency checks

1. **External vs local machine.** One job-A point (J = 23, ν = 51.06) and one job-B point (J = 29, ν = 17.5, s5001) were recomputed here
   with `ef_points` (the same code path as `ext/point.py`). The replicas agree to `1.5·10⁻¹²` (n1) and `≤ 4.5·10⁻⁷` relative (n3: the
   reference record), which is `≤ 0.5 c₂` in Γ. The n3 `logN` agrees to `3·10⁻¹⁰`.
2. **Overlaps.** No external point shares a key with an existing direct point. Because n1 is deterministic, an overlap would only have
   reproduced it (check 1). Old sources that collide (HC/HIGHSPIN J = 31, NEARGAP2/HIGHSPIN J = 29) are bit-identical.
3. **Joins.** The external x = 0.05 node against SMALLNU's x ≈ 0.046 n1 node: J = 23 `+0.5σ / −0.4σ`, J = 27 `−1.4σ / +1.6σ`, J = 31
   `−0.1σ / −0.0σ` (C2 / C3H). The degree-6 shape fit of the external bulk rows has rms/rms 0.05–0.06, at the MC level.
4. **Direct vs K₂** (bulk, both kinds): median `|Δ|` 58 / 64 / 89 c₂ and 1.6 / 1.7 / 2.0σ (J = 23 / 27 / 31). This is the same picture as
   HIGHSPIN's J = 29 (50 / 66): K₂ has real errors of `10²` c₂ that oscillate with x and partly cancel in λ (net +0.02–0.04, favourable).
5. **J = 29 streams.** Shared streams are bit-identical (840 checks); the independent blocks are χ²-consistent (§4); the pooled record
   (super-replicas) and the 31-stream mean of per-stream values agree within 0.4σ. The largest gap is ν = 15 C3H, `+209` vs `+439 ± 644`,
   a ratio-of-means vs mean-of-ratios difference driven by the z = 11.7 stream.
6. **Normalisation transform.** On the old nodes it holds to `≤ 0.13 c₂` (J = 29, baseline), and to 5·10⁻¹¹ for J = 23–31, where the FULL
   normalisation needs no transform.

## 6. Revised positivity inventory, Result D at Δ_g = 1000

| region / state | before (HIGHSPIN §6 + SMALLNU) | now | basis |
|---|---|---|---|
| flat limit | [certified] | [certified] | CERT_CLOSE |
| `m ≥ 1.3Λ`, `b̂ < 5` | [NE±] (`L ≥ 8.7`) | unchanged | POS_ERRORS §2.4 |
| J = 5 `+−` (binding) | [NE±] (`L = 78`) | unchanged | NEARGAP2 §1.3 |
| J = 0, 2, 4, 6 | [NE±] | unchanged | NEARGAP2 |
| **J = 7** | [NE + M] | **[NE±]**, `+0.2683`, `L ≥ 6.7` sign-free (own direct row) | §3.4 |
| J = 8 | [NE±] (λ-level, `L ≥ 4.3`) | unchanged | NEARGAP2 |
| J = 9–11, 13–15 | [NE + M] (J = 15 own K₂ row) | unchanged | NEARGAP2 |
| J = 12 | [NE±], `L ≥ 3.6` | unchanged | SMALLNU |
| J = 16 | [NE±], `L ≥ 1.7` | unchanged | SMALLNU |
| J = 17–20, 22, 24–26, 28 | [NE + M] (K₂ rows; edges computed at 18, 19, 20, 22, 24, 25, 28; J = 26 no row [M]) | unchanged | HIGHSPIN, NEARGAP2 |
| J = 21 | [NE±] (thin), `L ≥ 1.20` | unchanged | SMALLNU |
| **J = 23** | [NE + M] (K₂ truncation) | **1.03–1.10Λ [NE±]** (`L ≥ 1.16`); **twist edge [NE + M]** (`+0.923`, `L = 0.85`; representation-limited, 1.52 noise-aware) | §3.2 |
| **J = 27** | [NE + M] (K₂ truncation) | **[NE±]**, `+0.976` at the twist edge, `L ≥ 1.64` | §3.1 |
| J = 29 | [NE±] (thin; signed `L = 1.23`; ×3 0.94) | [NE±] unchanged: signed `1.24`, **×3 still 0.94** | §4 |
| **J = 31** | [NE + M] (K₂ truncation, `L = 0.08`) | **1.075–1.10Λ [NE±]** (`L ≥ 1.53`); **twist edge – 1.05Λ [NE + M]** (`+0.81 … +0.90`; blocked by the unmeasured edge tail, no GK15 run) | §3.3 |
| J = 32 | [NE + M] (window measured, no budget) | unchanged | SMALLNU |
| J = 36, 40, 48 | [NE + M]; J = 48 edge [M] | unchanged | — |
| even J = 30, 34–46 (no rows), even ≥ 50, odd ≥ 33 | [M] / [D-sketch + M] (`H_LJ(1000)`) | unchanged; the exploratory J = 400 job failed, no data | HIGHSPIN §5 |
| `m ≥ 1.3Λ`, `b̂ ≥ 5` | [NE + M] / [M] (`H_R`) | unchanged | POS_ERRORS |

**Bottom line at Δ_g = 1000.** Error-bounded ([NE±] or better):
* the flat limit, and `b̂ < 5` at `m ≥ 1.3Λ`;
* near the gap, J = 0, 2, 4, 5, 6, **7**, 8, 12, 16, 21, **27**, 29 (signed);
* J = **23 above 1.03Λ** and J = **31 at 1.075–1.10Λ**.

J = 23's and J = 31's twist edges stay [NE + M], with the missing pieces named: a noise-robust representation term (J = 23) and a GK15
edge window (J = 31). Nothing computed is negative.

## 7. Suggested corrections to existing documents (not applied)

1. **`ext/README_EXT.md` / `jobs.txt`:**
   * **B's "16 independent seeds" share streams.** Use `Z2_SEED` spaced by ≥ `R = 8`.
   * **A's seed is irrelevant.** n1 always uses Sobol seeds 2000 + r.
   * C (J = 400) failed.
2. **SMALLNU §3.1 / §0 item 2:** "to close a sign-free ×3 envelope … ≈ 16× the n3 points" → false. With σ → 0 on ν = 10–17.5 it gives
   `L = 0.97`. The n1 node at ν = 20 (C3H σ = 727) and heavy-tailed n3 streams bind (§4 here).
3. **HIGHSPIN §0 / §6, SMALLNU §0, and the update boxes in NEARGAP2 / POS_ERRORS / UNIFIED_THEOREM:**
   * "J = 23 / 27 / 31 stay [NE + M], limited by the K₂ truncation" → J = 27 [NE±] (`L ≥ 1.64`); J = 23 [NE±] at ≥ 1.03Λ, twist edge
     [NE + M] (representation); J = 31 [NE±] at ≥ 1.075Λ, twist edge to 1.05Λ [NE + M] (edge tail).
   * "J = 7 … [NE + M]" → J = 7 [NE±], `L ≥ 6.7`.
4. **HIGHSPIN §4 calibration:** "actual K₂ error 4–5 % of the step (three states)" → "0.6–8 % over seven states, always favourable".
5. **`hs_odd.py` (not changed):**
   * (a) Without GK nodes the spline edge holds the last node over `(x_l, Λ/m]` with no budget term for it. Add HIGHSPIN's E_miss
     envelope for that tail (J = 31: 1.85 N_s at the twist edge). This affects SMALLNU's J = 31 numbers, which did not close anyway.
   * (b) E_repr's degree-(KB+1) probe on a noisy direct row is MC-amplified (J = 23: its jackknife σ is 5× the main fit's). Consider
     comparing within jackknife errors, or using lower degrees and a model-free bulk.
6. **HIGHSPIN §3 table / task wording "21/23/27 had edge runs":** correct; J = 31 has **no GK15 run**, only HC's 8 equispaced nodes
   x = 0.8–0.9.

## 8. Reproduce
```bash
cd ads4/code/extfold          # every command under: ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1
python3 ../../ext/unpack.py   # only if ext/unpacked/ is missing (writes ext/unpacked/, git-ignored)
python3 ef_merge.py                                   # -> cache/pts/{old,base,all}/, cache/ef_streams_J29.pkl, cache/ef_merge.txt
python3 ef_shadow.py                                  # (1007, 7) odd shadow -> cache/ef_p2shadow.pkl
python3 ef_points.py n1 23 1023 51.061483 4242 A_J23; python3 ef_points.py n3 29 1029 17.5 5001 B_J29_s5001   # §5.1
python3 ef_oddrun.py 23 1.03,1.05,1.075,1.1; python3 ef_oddrun.py 27 1.03,1.04,1.05,1.075,1.1; python3 ef_oddrun.py 31 1.03,1.04,1.05,1.075,1.1
python3 ef_window.py 23; python3 ef_window.py 27; python3 ef_window.py 31; python3 ef_window.py 29
python3 ef_window.py 29 old old_repro                 # SMALLNU reproduction (likewise 23, 27, 31)
python3 ef_odddiag.py 23 1.03,1.05,1.075,1.1; python3 ef_odddiag.py 27 1.03,1.04,1.05,1.075,1.1; python3 ef_odddiag.py 31 1.03,1.04,1.05,1.075,1.1
python3 ef_j7.py; python3 ef_j7nodes.py
python3 ef_streams.py; python3 ef_j29env.py
python3 ef_summary.py                                 # -> cache/ef_summary.txt / .json
```

## 9. Caveats

1. The direct n1/n3 integrals are the project's reference standard [NE], not certified.
2. **J = 7:**
   * Its own-row records rest on a shadow computed here; the procedure reproduces two cached shadows bit for bit.
   * Its edge (x > 0.775) and window (ν < 50) are the 10⁴ conventions: ×47 above the measured row maximum, but not measured there.
3. **J = 23:** the noise-aware `L = 1.52` uses a representation term that differs from the protocol. It is quoted as [NE + M] and not used
   for the grade.
4. **J = 31's corrected budget** uses HIGHSPIN's light-cone envelope `3A((1 − x_l)/(1 − x))²` for the tail. This is a model, ×3 above
   the growth measured at J = 29.
5. Below the lowest measured node, `|Γ| ≤ 10⁴ c₂` remains a model input everywhere: ν < 20 for J = 23–31, ν < 10 for J = 29, ν < 50 for
   J = 7.
6. **J = 29:** the heavy-tailed streams make the jackknife σ itself uncertain. The χ² tests use block σ and are not sensitive to a single
   outlier stream at the level that matters here.
7. Only the masses listed were evaluated; between them, `H_cont` as before.
