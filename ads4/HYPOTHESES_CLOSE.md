# Closing the computational hypotheses of Result D (`H_edge`, `H_OO`, `Δ_g = 300` small-`ν`, `K₂` beyond its criterion)

*20–21 Nov 2026, hypothesis-closing agent. Unrefereed, not committed. No existing document or code file was edited. New code
`ads4/code/spinfunc/hc_*.py` (+ `hc_queue*.sh`), caches `spinfunc/cache/hc_*` (one file per state; no shared reference keys).
The functional `Ω^EE` (min-norm gauge, NF02 weights) and the analysis protocol (sm_chainfix chains: curvature double count removed,
`SM_CURV=row`, κ image; Jacobi fix; Hampel despike, odd rows included as in `st_odd_despike.py`) are FIXED; nothing was tuned on
the signs found here.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

**STATUS: in progress** (see hand-over at the end).

## 1. Methods and validation

### 1.1 Direct points (`hc_points.py`, `hc_nodes.py`)
* **Edge and bulk nodes (`HC_METHOD=n1`).** The n1 real re-centred contour (`n1_run.point` integrand for even `J`; `p8_odd_pw`'s
  point-wise `p3_odd_action` for odd `J`, kind `EE` = `w0.make_build_odd_EE`), with
  * the `Q_J` overflow rescaling of z2 (`set_qsc` at the n1 saddle; without it every heavy-reference point at `ν ≲ 650` is masked,
    which is why `ee1000_direct` had to use `K₂` references),
  * the scalar `Q_J` by the Legendre recurrence and the spin-`J` Fourier floor `fw_frame.patch_S` (FOLD_WINDOW §1),
  * a moment-matched Gaussian proposal (`n1_direct.adapt`, pilot seed 5) and RQMC seeds `2000 + r`, `r = 0…7`, **shared by the state and
    its reference** (`Δ_ref = 2·10⁴`, same `ν`). At `(1040, 40)`, `x = 0.9` the adaptive proposal lowers the element errors ×8 at equal
    cost (fixed `s = 1.2`: 28–40 %; adaptive: 4–5 % at `2⁹·2·8` points; `0.9–1.4 %` on the helicity diagonals at `2¹³·2·8`).
  * `fw_frame.patch_conserved_clean` is **not** applied (needed only at `J ≥ 150`): it rotates the SVD basis of the conserved `⟨TTO⟩`
    null space, and the points must be in the structure basis of the `K₂` rows they are injected into.
* **Small-`ν` nodes (`HC_METHOD=n3`).** z2's n3 hybrid contour (`ν_f = 13` below 13), proposal adapted on the state and re-used for the
  reference (common random numbers), `Z2_NUSPLIT=5` (n3 also at `ν = 8.3, 10`; no thimbles), `2¹⁰·2·8` points.
* **Records.** `M0 = E b⁻ᵀ` (even) / `E/b_OO` (odd), `M2 = M4 = 0`, `logC = logN`, `logK` of the shadow; the shadow `b` is the one whose
  stored fit residual equals the `shadow_resid` of the row's `K₂` records (fw_fix's rule), so direct and `K₂` records share the
  structure/shadow basis. One cache file per state (`cache/hc_pts_<kind>_J<J>_D<Δ>.pkl`), holding the state and its own references.
* **Validation (state values).** `(1040, 40)`, `x = 0.8`, EE: the adaptive point agrees element by element with the independent
  `ee1000_n1` point (fixed proposal, `2¹⁴` points): ratios `1.0079 / 0.9977 / 0.9962 / 0.9993` against errors `0.64 % / 0.18–0.5 %`
  (≤ 1.2σ), identical `logN` [NE]. So the hc points are in the same structure basis and normalisation as the project's direct points.

### 1.2 Chains (`hc_chain.py`)
* `sm_chainfix`'s modes re-run in-process (`ee1000`, `w2`, `z4`), curvature double count removed (`SM_CURV=row`, κ image), Jacobi fix,
  odd-row Hampel despike as `st_odd_despike.py`. The hc records are injected into the record loaders of the named rows (both kinds,
  state and reference); a `K₂` node at the same `ν` is replaced; the helicity basis is computed **without** the injected records.
* **Reproduction:** with nothing injected, `ee1000` (`J = 36, 40, 17, 31`, `XMAX = 0.8`) gives `+0.9669 / +0.9794 / +0.8389 / +0.9945`,
  the statement's numbers exactly.

## Hand-over (kept current)

* Running: `hc_queueA.sh` (edge nodes), `hc_queueOO.sh` (OO K₂ rows), then `hc_queueB.sh`, `hc_queueA2.sh` (wait on the first two).
  All resumable: after a container restart re-launch the same scripts from `ads4/code/spinfunc/` with
  `setsid nohup ./hc_queueX.sh > cache/hc_queueX.log 2>&1 &` (finished points are skipped).
