# Graviton fold window and bulk at large `J` near the gap: precision of the spin-`J` engine, direct values, `J = 48`

*20 Nov 2026, fold-window agent. Unrefereed, not committed. No existing document or code file was edited. New code in
`ads4/code/spinfunc/`: `fw_frame.py` (engine fixes), `fw_shadow.py` (robust shadow matrix), `fw_points.py` (direct points),
`fw_rows.py` (K₂ / basis records), `fw_gamma.py` (Γ and smooth-part extraction), `fw_scalar.py` (exact-scalar validation),
`fw_chain48.py` (`J = 48` chain), queues `fw_queueA.sh` (superseded), `fw_queueA2.sh`, `fw_queueB.sh`, `fw_after.sh`.
Caches `spinfunc/cache/fw_*`. Project modules are imported read-only; every patch is in-process.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

**STATUS: in progress** (see the hand-over at the end; this file is kept current).

## 0. Summary (current)

1. **The spin-`J` crystal engine loses precision like `ε²4^J`** [NE + D for the mechanism; §1]. This affects every graviton record
   at large `J`: K₂ rows, helicity bases, direct points, and the shadow matrix.
   * The spin-`J` index contraction `Σ_n T_n S_{−n}/κ_n` (l5_spin) divides by `κ_{±J} = 2^{−J}`. At the edge modes `|n| ≈ J`,
     both sides' true Fourier coefficients are far below the FFT round-off floor. Two round-off coefficients divided by `κ_J` give a
     relative error `≈ ε²4^J√J`.
   * Measured by rotating the sampling phase of the null circle (an exact identity), the relative change of the order-0 `M0` is:

     | point | unpatched engine | fixed engine |
     |---|---|---|
     | `(1048, 48)`, `x = 0.3` | `2.5·10⁻⁴` (`5.7·10⁻³` at another phase) | `≤ 10⁻¹⁴` |
     | `(2·10⁴, 48)`, `x_ref = 0.3` | (unpatched vs fixed: `1.4·10⁻³`) | `≤ 10⁻¹³` |
     | `(1040, 40)`, `x = 0.3` | (unpatched vs fixed: `3.8·10⁻⁸`) | `≤ 2·10⁻¹⁴` |

     Extrapolated: `O(1)` at `J ≈ 64`, `≫ 1` at `J ≥ 100`.
   * **Fix** (`fw_frame.patch_S`): Fourier coefficients below `10³ε√(2J+1)·max|samples|` (per TPS coefficient, or per point in the
     point-wise backend) are set to zero before the contraction. Dropped true terms are `~|F_n G_{−n}|/κ_n` with both coefficients at the
     floor. These are negligible because both spectra decay: at `J = 48`, `|F_n| ≈ 3·10⁻⁵` at `n = 24` and `4·10⁻¹²` at `n = 32`.
   * Three further large-`J` defects were found and fixed (§1.2–1.4):
     * the scalar `Q_J` was evaluated in the power basis of `P_J`, with cancellation `~2^J` (`6.6·10⁻⁷` at `J = 24`, `10⁻¹` at
       `J = 100`). It is now evaluated by the Legendre recurrence;
     * **the shadow matrix `b` is wrong at `J ≥ 56`** (unscaled least squares over random polarisation triples). The cached
       `(1064, 64)` value has `det b = 8·10⁻¹⁸`, and seed 11 gives `det 0.57`. The equilibrated fit is seed- and `ξ`-independent to
       `10⁻⁴` (`fw_shadow.py`). The `J = 48` shadows were fine (`det` agrees to 3 digits);
     * the conserved `⟨TTO⟩` null space fails (`SVD did not converge`) at `J ≥ 150`, from over/underflowing random sample points.
       Those points are now skipped; the null space is unchanged at `J = 48` (`1 − cos ≤ 3·10⁻¹⁶`).
   * **Consequence:** graviton records with `J ≳ 48` from the unpatched engine carry errors of `m²·10⁻⁴ ≈ 10²` `c₂` (`J = 48`,
     `Δ ≈ 10³`) and are meaningless at `J ≥ 56` (shadow and contraction). This covers the `ee1000` `J = 48` rows, refb `J = 56, 64`,
     `(1064, 64)`, and r1's `J = 72`. EE_DG1000_LARGEJ's "K₂ divergent at `J = 64`" diagnosis was made with such records and should
     be re-checked [I].
2. **Direct fold-window values at large `J`** (§2): in progress. First finding: the direct method's **statistical** precision at
   `J²/m ≫ 1` near the gap is far from the tolerance [NE].
   * Relative error of the state/reference helicity ratio, per `2¹⁰·2·8` points: `≈ 3 %` (EE, dir) at `(400, 100)`.
   * Projected per-node `Γ` errors and tolerances:
     * `(400, 100)`: `δΓ ≈ m²|D_h|δ ≈ 10³` `c₂` (dir), `≈ 400` (CV), against a fold tolerance of 36;
     * `(1064, 64)`: `≈ 400` (CV), against 3580.
   * The cause is a sign problem: the spin-`J` prefactor oscillates across the Gaussian proposal on the spin-independent thimble.
3. (pending) Bulk test.
4. (pending) `J = 48`.

## 1. The spin-`J` precision problem and its fix [NE + D]

### 1.1 Contraction edge modes
* l5_spin contracts spin-`J` tensors from their values on the null circle `w(φ) = (1 − cos φ, 1 + cos φ, sin φ)`. The formula is
  `⟨T, S⟩ = Σ_{n=−J}^{J} T_n S_{−n}/κ_n`, with `κ_n = (−1)^{J+n}2^{−J}C(2J, J+n)`.
* `κ_0 ≈ 2^J/√(πJ)` is large and `κ_{±J} = 2^{−J}` is tiny.
* At `(1048, 48)`, `x = 0.3` (scratch diagnostics, reproduced by `fw_frame` docstring numbers):
  * the `F`-side spectrum is `|F_n| = 3, 5, 0.3, 3·10⁻⁵, 4·10⁻¹², 6·10⁻¹⁵, 3·10⁻¹⁵` at `n = 0, 8, …, 48` (the last three at round-off);
  * the `G` side is alike;
  * the true contraction is `~10⁻¹¹`.
* The noise sum `Σ ε max|F|·|G_{−n}|/κ_n + (F ↔ G)` is `1.5·10⁻¹⁴`, i.e. relative `≈ 2·10⁻³`, dominated by `n = ±48`.
* Rotating the sampling phase by `0.37` (exact invariance) moves `M0` by `2.5·10⁻⁴`.
* Tried and rejected: a complex Lorentz frame adapted to the spin-`J` covectors (`FW_FRAME=1`, kept as an option). It makes the
  error worse (`1.7·10⁻³`) because the noise is at the edge modes, not in the bulk of the spectrum.
* With the floor (`FW_THR = 10³`): phase invariance `≤ 10⁻¹³` at `J = 40, 48`, for the state and the reference, for orders `M0, M2, M4`.

### 1.2 Scalar `Q_J`
* `l1_crystal.build` uses `Q_J = Σ_p c_p (Z·W)^p (Z²W²)^{(J−p)/2}`, with the power-basis Legendre coefficients `c_p`.
* Replacement (even `J`): `Q_J = (Z²W²)^{J/2} P_J(t)`, with `t = Z·W/√(Z²W²)` and the recurrence. It composes with z2's overflow
  rescaling `QSC`.
* Agreement: `≤ 6·10⁻¹⁴` at `J ≤ 8`. The difference is `3·10⁻¹¹` at `J = 16` and `1.4·10⁻⁶` at `J = 24`; that is the old
  evaluation's error.

### 1.3 Shadow matrix
* l5_shadow_matrix fits `S_c(t) = Σ_b b_{cb} r_b(t)` over 6 random null triples by an unscaled `lstsq`. At large `J` one triple
  dominates.
* `fw_shadow.robust` uses row and column equilibration and 12 triples.
* Seed/`ξ` dependence of the robust fit:

  | point | seed/`ξ` dependence |
  |---|---|
  | `(1040, 40)` | `1.2·10⁻⁵` |
  | `(1048, 48)` | `2·10⁻⁵` |
  | `(2·10⁴, 48)` | `3·10⁻⁹` |
  | `(1064, 64)` | `≤ 10⁻⁴` (5 seeds/`ξ`) |

* The old `(1064, 64)` (seed 7) is singular. The robust `b` has `det 0.573`; that is the seed-11 value.
* The fit residual (Laplace truncation of the shadow integral at order 2) is `≤ 3·10⁻⁵` at `Δ ≥ 1000`, but `9·10⁻⁴`
  (`(364, 64)`), `1.7·10⁻³` (`(400, 100)`) and `3·10⁻³` (`(450, 150)`) at `Δ_g = 300`.
* So at 300 the structure basis itself carries a `10⁻³` systematic. In `Γ` it enters as `≈ m²·10⁻³ ≈ 150` `c₂` times the other
  helicity's shape, which is in the oscillating class [D-sketch].

### 1.4 Conserved structures at `J ≥ 150`
* Random sample points with unbounded null `z` overflow. The fix normalises every `z` (homogeneity) and skips non-finite points.

### 1.5 Helicity basis at `J²/Δ ≫ 1`: generalised (non-symmetric) rank-1 pair [NE]
* The shadow matrix is itself a Laplace series (around `x₅ = e`). Change of `b` from order 1 to order 2 (`fw_shadow` equilibrated
  fit, scratch `shord.py`):

  | state | `J²/Δ` | change of `b`, order 1 → 2 |
  |---|---|---|
  | `(1040, 40)` | 1.5 | `7·10⁻⁴` |
  | `(1048, 48)` | 2.2 | `1.0·10⁻³` |
  | `(1064, 64)` | 3.8 | `1.8·10⁻³` |
  | `(1100, 100)` | 9.1 | `4·10⁻³` |
  | `(364, 64)` | 11.3 | `1.5·10⁻²` |
  | `(400, 100)` | 25 | `3·10⁻²` (not converged) |

* Consequence at `(400, 100)`: the order-0 matrices `M0 = S b^{−T}` are not symmetric (`|M01 − M10|/|M| ≈ 0.7`). l8's
  symmetric rank-1 basis fails there: `sv3 = 6·10⁻⁵ ≈ sv2/2.5`, `imag 0.37`. This is **not** an engine error: phase invariance holds to
  `2·10⁻¹⁴`, and the result is independent of the floor threshold (`10¹–10⁶`) to `10⁻⁷`.
* A wrong `b` is a fixed linear map on one index, so the x-independent two-helicity structure survives as `M0(x) = A diag(d(x)) Bᵀ`.
  `fw_gamma.gbasis` finds the two rank-1 matrices `a_h b_hᵀ` in the 2-dimensional span of `{M0(x)}` and uses `Ti = A⁻¹` on the
  left and `B⁻¹` on the right. It fixes a constant phase per helicity and identifies `++` by flatness, as l8 does.
* Quality of the generalised basis:

  | point | `x` window | `sv3/sv1` | off-diagonal residue | `imag` |
  |---|---|---|---|---|
  | `(400, 100)` | `[0.15, 0.4]` | `6·10⁻⁸` | `6·10⁻⁸` | `2·10⁻¹³` |
  | `(2·10⁴, 48)` | `[0.3, 0.6]` | `3·10⁻⁷` | | |

  On `x ∈ [0.3, 0.6]` the `(400, 100)` span acquires a third direction at `7·10⁻⁴`, and on `[0.15, 0.8]` at `0.06`, so the window is
  restricted to `x ≤ 0.4` for states and `x_ref ∈ [0.3, 0.6]` for references. At `x_ref = 0.8` the fixed-engine `(2·10⁴, 48)`
  order-0 action is asymmetric at `3·10⁻²`, and `0.7` at `6·10⁻⁵`.
* Check at `(1040, 40)`: the generalised and l8 helicity diagonals agree up to one real constant per helicity (absorbed in `a_h`) to
  `10⁻⁷`, at every `x`.

## 2. Direct fold-window values (in progress)

Method: as EE_DG300_SMALLNU (n3 contour, `ν_f = 13` below 13, state and `Δ_ref = 2·10⁴` on common random numbers, proposal adapted on
the state), at the gap states `Δ = Δ_g + J` (`τ = Δ − J = Δ_g`, the code's twist convention `α = Δ − J − 2 ≥ Δ_g − 2`). Nodes
`ν = 10, 13.3, 16.7, 20.8, 25, 33.4`, `2¹⁰·2·8` points per state.

Precision tests (`ν = 16.7`, relative jackknife error of the state/reference ratio of an element; fixed engine):

| state | `J²/m` | points | EE dir | EE CV | C3H dir | C3H CV |
|---|---|---|---|---|---|---|
| `(1064, 64)` | 3.9 | `2⁹·2·4` | 1.2 % | 0.35 % | 1.9 % | 1.0 % |
| `(364, 64)` | 11.5 | `2⁹·2·4` | 7.4 % | 1.1–1.5 % | 5.4 % | 3.1 % |
| `(400, 100)` | 26 | `2⁸·2·4` | 8.8 % | 3.2 % | 48 % | 41 % |
| `(400, 100)`, unpatched engine | 26 | `2¹⁰·2·4` | 24 % | 34 % | 30 % | 37 % |

### 2.1 Structure-basis consistency (`fw_fix.py`) [D, code]
* l5_spin.conserved returns the null space in an SVD-arbitrary basis that depends on its random sample points.
  `fw_frame.patch_conserved` (needed at `J ≥ 150`) changes those points, so it rotates the basis: `P = pinv(sol_pat)·sol_unp` is a
  non-trivial `2×2`, with residual `10⁻¹⁴`.
* So records made before the patch was wired in (10:10 UTC) and after it are in different structure bases. K₂ records also store
  `M = S b^{−T}` with whatever shadow `b` was current.
* `fw_fix` brings every record to one basis at load time (exact linear algebra):
  * K₂ records: the `b` actually used is identified by its stored `shadow_resid`, and `M ← M b_usedᵀ b_cur^{−T}`;
  * direct points before the patch (`(400, 100)` and its reference at `ν = 16.7`): `E ← P^{−T} E Pt^{−1}`.
* The `J = 48` refb records from before the patch were removed from this agent's cache and recomputed (`fw_rmrefb48.py`).
* The fw shadows live in `cache/fw_shadowk2.pkl`, first run in `cache/fw_shadowrun1.pkl`. These names were chosen not to match the
  `*shadow_k2*` globs of other agents' loaders. A transient copy `cache/fw_shadow_k2.pkl` is used only by an fw_rows process started
  before the rename and is to be deleted once it ends.

### 2.2 `Δ_g = 300`, `J = 100`, gap state `(400, 100)` (`m = 385.6`, `J²/m = 26`) [NE, weak]
Setup:
* Nodes `ν = 13.3, 16.7, 20.8, 25.0, 33.4`, n3 contour, `2¹⁰·2·8` points.
* Bases: generalised, state on `x ∈ [0.15, 0.4]` (`sv3 = 6·10⁻⁸`), reference on `x_ref ∈ [0.3, 0.5]` (`sv3 = 2·10⁻⁶`).
* `a_h` is fitted on the 4 nodes with `ν ≥ 15`.
* Fit `y = Γ − κ-image = c D_h + s xD_h′/(2J+1) + r` (`cs`), plus variants with `x²D_h` (`cs+x2`) or a linear term (`cs+lin`).
* Errors are `lsq` (residual-scaled, i.e. the fit-model error) and `jk` (replica jackknife).

Exact-scalar validation:
* Scalar on the same samples (`fw_scalar.py`, against sm_exact `n0t`): `Re(d − d_ref)` is within `0.3–1.6σ` at all 5 nodes.
  `Δc₂ = −105 … +405`, with errors `32–306`.
* Smooth fit applied to the exact scalar at the same nodes (`fw_scalfit.py`): `r_s = −2.1 ± 8.1` (`cs`), `+1.4 ± 2.2` (`cs+x2`),
  `−1.0 ± 10.9` (`cs+lin`). So the extraction has no bias at the `≈ 10` `c₂` level on a known object.

Graviton, smooth part `r` (`c₂`), CV estimator (scalar control variate; dir in brackets):

| key | `cs` | `cs+x2` | `cs+lin` |
|---|---|---|---|
| C2 (EE) `++` | `−96 ± 214 ± 24` (`−124 ± 293 ± 90`) | `−5 ± 52 ± 14` | `−68 ± 290 ± 20` |
| C2 `+−` | `+156 ± 172 ± 27` (`+128 ± 130 ± 92`) | `+104 ± 184 ± 16` | `+191 ± 204 ± 23` |
| C3H `++` | `−184 ± 202 ± 125` (`−137 ± 276 ± 110`) | `−97 ± 4 ± 119` | `−169 ± 288 ± 141` |
| C3H `+−` | `+31 ± 84 ± 131` (`−16 ± 59 ± 113`) | `+65 ± 53 ± 125` | `+24 ± 120 ± 147` |

* Per-node `Γ − κ` errors (jk): C2 `20–140` (CV), `150–440` (dir); C3H `170–540` (both estimators).
* A preliminary 4-node fit (without `ν = 20.8`) gave C2 `r ≈ +4000 ± 800`. That shows a 2-parameter-short fit is unstable; it is
  superseded.

**Reading.**
* The fold-window smooth part at `(400, 100)` is consistent with zero, `|r| ≲ 200` `c₂` at 1σ, in every key and estimator.
* It is **not** shown to be above the fold tolerance (`−36` `c₂`; "≳ −30" in the task). The uncertainty is dominated by the fit
  model with 5 nodes (`lsq` errors `50–290`) and, for C3H, by statistics (`≈ 125`).
* Grade: [NE, weak], and the `J = 100`, `Δ_g = 300` fold item stays **[I]**.
* Closing it needs `≳ 10` nodes with C3H statistics `×10` (`≈ 100×` the points), or a method without the state-side sign problem.

## Hand-over (kept current)

* Running (detached, resumable; after a restart re-run the same commands, finished points are skipped):
  * `cd ads4/code/spinfunc && setsid nohup ./fw_queueA2.sh > cache/fw_queueA2.log 2>&1 &` (direct nodes);
  * `setsid nohup ./fw_queueB.sh > cache/fw_queueB2.log 2>&1 &` (basis records, `J = 48` rows).
* Engine fixes: `fw_points.setup()` (direct) / `fw_rows.setup()` (K₂). They are in-process only.
