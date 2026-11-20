# Graviton fold window and bulk at large `J` near the gap: precision of the spin-`J` engine, direct values, `J = 48`

*20 Nov 2026, fold-window agent. Unrefereed, not committed. No existing document or code file was edited. New code in
`ads4/code/spinfunc/`: `fw_frame.py` (engine fixes), `fw_shadow.py` (robust shadow matrix), `fw_points.py` (direct points),
`fw_rows.py` (K₂ / basis records), `fw_fix.py` (load-time basis consistency), `fw_gamma.py` (Γ, generalised helicity basis,
smooth-part extraction), `fw_scalar.py` / `fw_scalfit.py` (exact-scalar validation), `fw_probe.py` (one-node precision),
`fw_chain48.py` (`J = 48`/`64` chains), `fw_rmrefb48.py`, queues `fw_queue{A,A2,A3,A4,A6,A7,B,B3,B4,B5,B6}.sh`, `fw_after.sh`.
Caches `spinfunc/cache/fw_*`. Project modules are imported read-only; every patch is in-process.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

**STATUS: final for this session (hand-over 22:00 UTC).** One job is still running: the reference block of the `(1200, 200)` probe (see the hand-over).

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
     `(1064, 64)`, and r1's `J = 72`. EE_DG1000_LARGEJ's "K₂ divergent at `J = 64`" diagnosis was made with such records. Re-checked
     here with the fixed engine (§3): K₂ is still not converged at `(1064, 64)`, so its conclusion stands. Its numbers there were
     engine noise.
2. **Fold-window smooth part (`ν ≈ 10–33`, direct n3 nodes, shared-seed reference, generalised bases, κ-image subtracted,
   oscillating class + constant fitted)** [NE; validated on the exact scalar at every state]:

   | `Δ_g`, `J` (state) | fold tolerance | smooth part `r` (`c₂`), C2 / C3H, CV | verdict |
   |---|---|---|---|
   | 300, 64 `(364, 64)` | 70 | 10 nodes: C2 `−51 ± 43` (classes `−72 … −3`), C3H `±50 ± 60` | **at the tolerance edge**: inside at 1σ, not at 2σ [NE, weak] → [I] |
   | 300, 100 `(400, 100)` | 36 | C2 `−104 … +116`, C3H `−133 … +51` (± `150–460`) | **not decided** [I] |
   | 300, 150 `(450, 150)` | 27 | one-node probe only: errors as at `J = 100` | **not decidable** with this method [NE + M] |
   | 1000, 64 `(1064, 64)` | 3580 | C2 `−200 ± 210`, C3H `−365 ± 700 / +440 ± 720` (weighted) | **inside tolerance** at 2σ [NE] |
   | 1000, 100 `(1100, 100)` | 1170 | C2 `+586 ± 179 / −378 ± 243` (weighted `cs`; classes `−378 … +586`); C3H `+592 ± 700 / +406 ± 590` (classes `−1588 … +1608`, errors `±600–1500`) | C2 **inside** at 2σ [NE]; C3H inside at about 1σ only [NE, weak] |
   | 1000, 200 `(1200, 200)` | 306 | one-node probe: state-block element errors 6.0 % (EE) / 8.5 % (C3H) per `2¹⁰·2·8` points, about `10×` those at `(1100, 100)` (0.6 % / 0.9 %); reference block still running | **not decidable** with this method [NE + M] |

   * Exact-scalar control on the same nodes: `|r_s| ≤ 3` `c₂` at every state. The scalar MC is consistent with the exact
     state/reference ratio at `≤ 2.1σ` on every node (`≤ 1.8σ` on all fold nodes).
   * The limitation at `Δ_g = 300` is statistical and structural:
     * the state-side sign problem (`J²/m = 11–53`): per-node `Γ` errors `10²–10³`;
     * C3H has `20–40 %` element errors at `J ≥ 100`;
     * with 5–6 nodes the constant is not separable from the oscillating class to better than `±10²`.
3. **Bulk smooth term (task 2)** [NE, negative; §3]:
   * Direct bulk points at `(1064, 64)`, `x = 0.2 / 0.4 / 0.6`, give `Γ − κ = +56 / −742 / −793` (C2 `++`, CV, ± `40–90`). That is the
     known κ-model failure at `x ≳ 0.4`.
   * The fixed-engine K₂ row at `(1064, 64)` disagrees with direct by `10³–3·10³` at `x ≥ 0.4`: K₂ is not converged at
     `J²/m = 3.9`.
   * So a bulk smooth term at `J > 48` is **not testable at the uniform tolerance** (`≈ 14` `c₂` at `J = 64`) with these data. It
     stays [I].
4. **`J = 48` (task 3)** [NE + M; §4]:
   * Rows recomputed with the fixed engine. The helicity bases are normal: `sv3 = 1.2·10⁻⁷ / 1.1·10⁻⁷` (was `7·10⁻⁵ / 8.5·10⁻⁶`).
   * The "smooth offset" is gone: `r = 0 ± 5` in all four keys (was `+280 … +650` / `−630 … −690`). **It was the engine's precision
     loss.**
   * The ee1000 shape chain with sm_curvfix gives min `λ/N_s = +0.932` (variants `+0.73 … +0.99`). **`n = 1` at `J = 48`, `Δ_g = 1000`.**
5. **Status of `n = 1` at large `J` near the gap** (on top of SMOOTH_BOUND's [NE + M, conditional]):
   * `Δ_g = 1000`:
     * `J = 48` is now **[NE + M]** (was [I]).
     * Fold window: `J = 64` inside tolerance [NE]; `J = 100` inside for C2 [NE], for C3H only at about 1σ [NE, weak]; `J = 200`
       has a one-node probe only.
     * The conditional "fold value `≥ −200` at `J ≳ 150–300`" is **not tested**. At `(1200, 200)` the state-block element errors
       are `6–8.5 %`, about `10×` those at `(1100, 100)`, which implies `r` errors `≳ 10³` against 306 [NE + M].
     * Bulk smooth term at `J > 48`: **[I]** (K₂ not converged at `J = 64`; 3 direct bulk points cannot separate the classes).
   * `Δ_g = 300`:
     * `J = 64`: fold smooth part `≈ −50 ± 45` `c₂` (C2), at the tolerance edge (70) [NE, weak].
     * `J = 100`: `|r| ≲ 200` with errors `≈ 200–500` against 36. **[I]**.
     * `J = 150`: not decidable with this method. **[I]**.
     * The "≳ −30 `c₂` at `J = 100–150`" condition is **not demonstrated**. Nothing measured points to a violation: every central
       value is `≥ −170` and consistent with zero.
   * **Overall: `n = 1` at large `J` near the gap stays [NE + M], conditional, with [I] narrowed.**
     * `J = 48` (1000) is closed.
     * The `J = 64` (1000) fold item is closed.
     * The `J = 100` (1000) C2 fold item is closed.
     * Still open: `Δ_g = 300` `J ≥ 64` (fold), C3H at 1000 `J ≥ 100`, `J ≥ 150` at 1000, the bulk term at `J > 48`, and the
       `ν < 10` / edge items of SMOOTH_BOUND (not addressed here).

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

**Final (6 nodes, `ν = 10.0` added; reference basis on `x_ref ∈ [0.15, 0.45]`, 7 records, `sv3 = 1.6·10⁻⁷`)** (`c₂`; `lsq ± jk`):

| key | CV `cs` | CV `cs+x2` | CV `cs+lin` | dir `cs` |
|---|---|---|---|---|
| C2 `++` | `−95 ± 165 ± 129` | `−68 ± 103 ± 156` | `−84 ± 201 ± 64` | `−15 ± 292 ± 306` |
| C2 `+−` | `−104 ± 460 ± 132` | `−165 ± 433 ± 161` | `+26 ± 280 ± 65` | `−22 ± 273 ± 311` |
| C3H `++` | `−27 ± 308 ± 158` | `+25 ± 190 ± 181` | `−91 ± 297 ± 128` | `+35 ± 359 ± 190` |
| C3H `+−` | `−38 ± 135 ± 157` | `−40 ± 166 ± 180` | `−10 ± 129 ± 130` | `−101 ± 150 ± 186` |

* The `ν = 10` node is noisy: jk `1100–2400`. Scalar check there: `0.2σ`.
* Exact-scalar fit on the 6 nodes: `r_s = −0.7 ± 6.7` / `+0.6 ± 1.8` / `−0.8 ± 8.4`.
* Central values are `−165 … +98` `c₂` and combined errors `≈ 200–500`.
* **The `(400, 100)` fold-window smooth part is consistent with zero at the `±200–500` level. It is not shown to be `≥ −36`.**
  Stays [I]. [NE, weak]: no evidence of a large negative term (none `< −170` centrally).

**Reading.**
* The fold-window smooth part at `(400, 100)` is consistent with zero, `|r| ≲ 200` `c₂` at 1σ, in every key and estimator.
* It is **not** shown to be above the fold tolerance (`−36` `c₂`; "≳ −30" in the task). The uncertainty is dominated by the fit
  model with 5 nodes (`lsq` errors `50–290`) and, for C3H, by statistics (`≈ 125`).
* Grade: [NE, weak], and the `J = 100`, `Δ_g = 300` fold item stays **[I]**.
* Closing it needs `≳ 10` nodes with C3H statistics `×10` (`≈ 100×` the points), or a method without the state-side sign problem.

### 2.3 `Δ_g = 1000`, `J = 64`, gap state `(1064, 64)` (`m = 1060.5`, `J²/m = 3.9`) [NE]

Setup:
* Nodes `ν = 13.3, 16.7, 20.8, 25.0, 33.4`, same method. Bases: state `sv3 = 6·10⁻¹⁰`, reference `2·10⁻¹²`.
* `(2J+1)x = 1.6 … 4.1`: the window lies in the first one or two lobes of `D_h`. So only the 3-parameter class `cs` is determined;
  `cs+x2` / `cs+lin` give `r = ±10³–10⁴` with errors of the same size.

Validation:
* Scalar on the same samples: `0.1–1.1σ` at all nodes.
* Exact-scalar smooth fit on the same nodes: `r_s = −3.1 ± 0.2` (`cs`), `+0.7`, `−5.5`. The exact scalar's own fold plateau is
  `m²ρ₂ ≈ 6 → −9` `c₂` here.

Graviton (`c₂`, `cs` class, `lsq ± jk`):

| key | CV | dir | `Γ − κ` at the nodes (CV, ± jk) |
|---|---|---|---|
| C2 `++` | `−253 ± 108 ± 164` | `−449 ± 390 ± 277` | `−111 … −344`, `+112` (± `80–385`) |
| C2 `+−` | `−242 ± 110 ± 159` | `−450 ± 392 ± 277` | similar |
| C3H `++` | `−162 ± 863 ± 547` | `−180 ± 337 ± 619` | `−838 … +1845` (± `180–1240`) |
| C3H `+−` | `+171 ± 834 ± 536` | `+198 ± 369 ± 616` | `−1802 … +790` (± `180–1170`) |

**Update with the 6th node (`ν = 10.0`, n3 `ν_f = 13`).** At `ν = 10` the node is noisy: `Γ − κ = −2653 ± 1619` (C2 `++`, CV) and
`−2727 ± 1366` (C3H `++`). Fits including it:

| fit | C2 `++` | C2 `+−` | C3H `++` | C3H `+−` |
|---|---|---|---|---|
| unweighted (`lsq ± jk`) | `−828 ± 455 ± 523` | `−821 ± 459 ± 522` | `−1157 ± 1002 ± 501` | `+1213 ± 1025 ± 471` |
| jk-weighted (`FW_WFIT=1`; error `√(χ²/dof)`-scaled) | **`−202 ± 183 ± 109`** | **`−185 ± 190 ± 98`** | **`−365 ± 568 ± 414`** | **`+441 ± 598 ± 400`** |

The `(400, 100)` weighted fits for comparison (CV): C2 `+116 ± 176` / `−71 ± 254`; C3H `−133 ± 201 ± 118` / `+51 ± 176 ± 120`.
The reference-basis window (`x_ref ∈ [0.3, 0.6]` vs `[0.15, 0.45]`) changes nothing at the printed precision.

**Reading.**
* The fold tolerance at `J = 64` (1000) is `3580` `c₂` on a `ν < 30` smooth term. Weighted fits: C2 `r ≥ −620`, C3H `r ≥ −1765` at
  2σ, both inside tolerance [NE]. Unweighted fits, dominated by the noisy `ν = 10` node: `r ≥ −2200` (C2) and `≥ −3400` (C3H) at 2σ,
  still inside.
* Earlier 5-node estimate: For C2, `r ≥ −250 − 2·200 ≈ −650` at 2σ: **inside
  tolerance** [NE].
* For C3H the statistical error is `≈ 10³`, and the extraction degenerates in the first lobe. Every node value is `≥ −1802 ± 1174`, so the
  whole window correction, oscillating or not, is above `−3580` at `≈ 1.5σ` [NE, weak].
* SMOOTH_BOUND's note applies: inside the first lobe a fold term cannot be told apart from `c₀D_h`. The statement is therefore about the
  size of the window correction, which is `|Γ − κ| ≲ 2·10³` `c₂`.

### 2.3b `Δ_g = 300`, `J = 64`, gap state `(364, 64)` (`m = 356.7`, `J²/m = 11.5`) [NE, weak]

Setup:
* 5 nodes (`ν = 13.3 … 33.4`). Bases `sv3 = 7·10⁻⁹` (state), `2·10⁻¹²` (reference).
* Scalar check `0.6–1.5σ`. Exact-scalar smooth fit `r_s = +1.2 ± 3.0` / `−0.4 ± 1.8` / `+1.3 ± 3.4`.
* Per-node `Γ − κ` errors (CV, jk): C2 `10–140`, C3H `110–370`.

Smooth part (`lsq ± jk`):

| key | CV `cs` | CV `cs+x2` | CV `cs+lin` | CV weighted `cs` | dir weighted `cs` |
|---|---|---|---|---|---|
| C2 `++` | `−129 ± 163 ± 22` | `−168 ± 227 ± 24` | `−116 ± 46 ± 20` | `+60 ± 81 ± 5` | `−7 ± 103` |
| C2 `+−` | `−130 ± 170 ± 23` | `−168 ± 242 ± 24` | `−117 ± 50 ± 21` | `−1 ± 132 ± 6` | `−10 ± 107` |
| C3H `++` | `+41 ± 32 ± 44` | `+23 ± 9 ± 48` | `+41 ± 45 ± 42` | `+54 ± 106 ± 52` | `+15 ± 52` |
| C3H `+−` | `−45 ± 35 ± 44` | `−25 ± 4 ± 48` | `−45 ± 49 ± 42` | `−57 ± 104 ± 48` | `−15 ± 52` |

**Reading.**
* The fold tolerance at `J = 64` (300) is `70` `c₂`.
* C3H is within `±60` in every class.
* C2 depends on the class, from `−155` (`cs+lin`, weighted) to `+60` (`cs`, weighted), with errors `35–240`. One class is below
  the tolerance at `1–2σ`; the others are consistent with zero.
* **Not decided** [I]. There is no evidence of a negative smooth term beyond the fit-class spread.

**With 10 nodes** (`ν = 11.7, 13.3, 15, 16.7, 18.5, 20.8, 23, 25, 29, 33.4`; `fw_queueA6.sh`):
* Checks:
  * scalar check `0.1–0.7σ` on the new nodes;
  * exact-scalar smooth fit `r_s = +0.3 ± 2.5` / `−0.2 ± 0.7` / `+0.3 ± 2.6`.
* The fit classes now agree.

| key | CV `cs` | CV `cs+x2` | CV `cs+lin` | CV weighted `cs` / `cs+lin` | dir weighted `cs` |
|---|---|---|---|---|---|
| C2 `++` | `−51 ± 32 ± 29` | `−50 ± 34 ± 26` | `−53 ± 14 ± 30` | `−3 ± 23 ± 4` / `−66 ± 9 ± 14` | `−31 ± 51 ± 22` |
| C2 `+−` | `−48 ± 33 ± 28` | `−47 ± 36 ± 26` | `−50 ± 16 ± 29` | `−19 ± 23 ± 6` / `−72 ± 12 ± 14` | `−29 ± 50 ± 21` |
| C3H `++` | `+50 ± 38 ± 52` | `+45 ± 30 ± 48` | `+51 ± 33 ± 54` | `−6 ± 45 ± 16` / `+27 ± 60 ± 42` | `+17 ± 42 ± 9` |
| C3H `+−` | `−54 ± 38 ± 51` | `−48 ± 27 ± 48` | `−55 ± 35 ± 53` | `+2 ± 44 ± 16` / `−29 ± 60 ± 44` | `−17 ± 43 ± 7` |

**Reading (10 nodes).**
* Every central value is `|r| ≤ 72` `c₂`, and most are within `−55 … +50`. Errors are `25–60` (combined).
* Against the fold tolerance `70`:
  * C2 and C3H, `cs` class: inside at 1σ but not at 2σ (C2 `−51 ± 43` → `−137` at 2σ);
  * C2 weighted `cs+lin`: `−66 … −72 ± 17`, i.e. at the edge.
* **`(364, 64)`: smooth part `≈ −50 ± 45` `c₂` (C2), at the tolerance edge** [NE, weak]. Not demonstrated; not excluded.
* The C2 and C3H smooth parts of `+−` are both negative in the CV `cs` class (`−48`, `−54`). If they enter the same helicity with
  comparable weights, they add. sm_bound's tolerance is per `Γ` (`ψ₂`-weighted), and C3H enters with `ψ₃u` [I, not combined here].

### 2.6 Two pitfalls found late: reference-cache collision, and `a_h` at first-lobe-only rows [D, code + NE]

**1. Reference-cache collision.**
* z2_points keys points by `(Δ, ν)`. Two states of the same `J` that share a node `ν` therefore share the `(2·10⁴, ν)` reference
  key, and the later run overwrites the earlier reference.
* Affected pairs and times:
  * `(364, 64)` overwrote `(1064, 64)`'s references at `ν = 13.3, 16.7, 20.8, 25, 33.4`, from about 16:35 UTC;
  * `(1100, 100)` overwrote `(400, 100)`'s from about 18:30.
* The analyses reported in §2.2 and §2.3 were run **before** the overwrites, on intact common-random-number pairs. Evidence: an
  unweighted rerun of `(1064, 64)` at 16:40 reproduced the earlier numbers exactly.
* Re-running `fw_gamma` on those two rows now gives different (invalid) numbers. For `ν = 16.7` at `(400, 100)` it also applies the
  pre-patch basis conversion to a post-patch reference.
* The reported numbers are recorded in `cache/fw_gamma.txt` (appended log).
* `fw_points` now writes every run to a per-state file `fw_pts_<kind>_J<J>_D<Δ>.pkl`, and `fw_gamma.fwpts` prefers it. Reproducing
  §2.2/§2.3 from caches needs the overwritten reference blocks recomputed: about 15 min each, using the stored state proposal in
  `fw_raw.pkl` and the same seeds. Not done.

**2. Per-row normalisation `a_h` at `(1100, 100)`.**
* All 5 fold nodes lie in the first lobes, `(2J+1)x = 2.4–6.1` (`x ≤ 0.03`). l9's `a_h` fit (`D_h, x², x⁴` on `x ≤ 0.25`, `ν ≥ 15`)
  is ill-conditioned there: it returned `a_{+−}/a_{++} = 6·10⁻³` instead of `O(1)`, and `Γ ~ 10⁵–10⁷`.
* A `D_h`-only `a_h` fit (`FW_AFIT=D`) also leaves `|Γ| ~ 3·10⁵`.
* The large-`J` local de-shift is `O(m²[cosh((2J+1)/m) − 1]) ~ 10⁴` `c₂` per unit `a_h` error, so a well-determined `a_h` is
  essential.
* At `(1064, 64)` the bulk node `x = 0.2` anchors `a_h`. At `(364, 64)` both `a_h` variants give the same smooth part (`FW_AFIT=D`,
  weighted `cs`: C2 `−7 ± 16 / −16 ± 14`, C3H `−6 ± 44 / +2 ± 44`).
* For `(1100, 100)` the nodes `x = 0.2, 0.1` were queued (`fw_queueA7.sh`) to fix `a_h`.

**`(1100, 100)` with the `x = 0.2` anchor** (`m = 1093.9`, `J²/m = 9.1`; 5 fold nodes `ν = 13.3–33.4` plus `ν = 218.8`, which
enters `a_h` only):
* A third bug was found here. `fw_fix` also converted the `(2·10⁴, 16.7)` reference as "unpatched", but after the overwrite of item 1
  it was the `(1100, 100)` run's patched reference. That one conversion had produced `|Γ| ~ 10⁵–10⁷` and, through `a_h`, the
  ill-conditioned values above. It is fixed: UNPATCHED now lists only the `(400, 100)` state point.
* The `a_h`-conditioning point above therefore stands only as a warning. With the anchor, `a_{++} = 2.4·10⁻⁵`, `a_{+−} = 0.9·10⁻⁵`.
* Scalar check at `x = 0.2`: `2.1σ`; fold nodes `0.1–1.3σ`. Exact-scalar smooth fit `r_s = +2.3 ± 1.2 / +0.4 / +1.9`.

| key | CV weighted `cs` | CV `cs` | CV weighted `cs+x2` / `cs+lin` |
|---|---|---|---|
| C2 `++` | `+661 ± 100 ± 133` | `+385 ± 180 ± 234` | `+360 ± 282 ± 307` / `+483 ± 190 ± 233` |
| C2 `+−` | `−347 ± 199 ± 135` | `−72 ± 147 ± 240` | `+369 ± 283 ± 307` / `+100 ± 191 ± 235` |
| C3H `++` | `+592 ± 639 ± 294` | `−518 ± 738 ± 817` | `−1166 ± 1153 ± 1049` / `−511 ± 831 ± 750` |
| C3H `+−` | `+407 ± 510 ± 301` | `+982 ± 457 ± 842` | `+1163 ± 1168 ± 1062` / `+875 ± 843 ± 763` |

**Reading.**
* The fold tolerance at `J = 100` (1000) is `1170`.
* C2: worst class `−347 ± 241`, i.e. `≥ −830` at 2σ: **inside** [NE].
* C3H: errors `600–1500`; the worst central value is `−1166 ± 1540`. Inside at about 1σ, not at 2σ [NE, weak].
* **With both anchors (`x = 0.1` and `0.2`; final).** `a_{++} = 2.41·10⁻⁵`, `a_{+−} = 0.92·10⁻⁵`.

  | key | CV weighted `cs` | CV `cs` | classes (weighted + unweighted) |
  |---|---|---|---|
  | C2 `++` | `+586 ± 124 ± 129` | `+250 ± 209 ± 232` | `−52 … +586` |
  | C2 `+−` | `−378 ± 204 ± 132` | `−91 ± 143 ± 240` | `−378 … +348` |
  | C3H `++` | `+592 ± 636 ± 292` | `−517 ± 736 ± 814` | `−1588 … +592` (± `700–1500`) |
  | C3H `+−` | `+406 ± 509 ± 300` | `+982 ± 457 ± 842` | `+406 … +1608` |

  * C2: worst `−378 ± 243` → `≥ −865` at 2σ, inside `1170` [NE].
  * C3H `++`: worst class `−1588 ± 1525` (`cs+x2`, unweighted), `cs` class `−517 ± 1100`. Inside at about 1σ only [NE, weak].
  * The anchors themselves have `Γ − κ ≈ ±(0.7–1.0)·10⁴` at `x = 0.1` and `±(2.5–3.7)·10³` at `x = 0.2`. That is the κ-model's
    failure at `(2J+1)x ≳ 20`, in the oscillating class, and gives no smooth-term information.

### 2.4 Precision probes (`fw_probe.py`; one node at `ν = 16.7`) [NE]

Relative jackknife error of the state/reference helicity-diagonal ratio, `2¹⁰·2·8` points per state:

| state | `J²/m` | EE dir | EE CV | C3H dir | C3H CV | reference basis `sv3` |
|---|---|---|---|---|---|---|
| `(1064, 64)` | 3.9 | 0.37 % | 0.16 % | 0.70 % | 0.49 % | `2·10⁻¹²` |
| `(400, 100)` | 26 | 2.9 % | 0.9 % | 31 % | 30 % | `1.6·10⁻⁷` |
| `(450, 150)` | 53 | 3.3 % | — (no exact scalar computed) | 22–24 % | — | `7·10⁻⁵` |
| `(1200, 200)` | 34 | state block: 6.0 % per element | — | state block: 8.5 % | — | `5·10⁻⁸` (x_ref ≤ 0.45; refs incomplete) |

* The rough proxy `m²|D_h|·δ` reproduces the fw_gamma jackknife errors only to within a factor `3–7`.
* At `(450, 150)` the per-node errors are like those at `(400, 100)`, i.e. `r` errors of a few hundred `c₂`, against a fold tolerance
  of `27` (`J = 150`, 300).
* So **`J = 150` at `Δ_g = 300` is not decidable with this method either** [NE for the precision, M for the projection]. No full row
  was computed.

### 2.5 Bulk points at `(1064, 64)`: `x = 0.2, 0.4, 0.6` (`ν = 212, 424, 636`; n3, `ν_f = ν`) [NE]
* Scalar check: `0.3σ / 1.8σ / 1.0σ`.
* `Γ − κ-image`, CV (± jk):

  | key | `x = 0.2` | `0.4` | `0.6` |
  |---|---|---|---|
  | C2 `++` | `+56 ± 93` | `−742 ± 38` | `−793 ± 68` |
  | C2 `+−` | `+79 ± 103` | `−671 ± 83` | `+816 ± 223` |
  | C3H `++` | `+284 ± 94` | `−752 ± 114` | `−175 ± 198` |
  | C3H `+−` | `−274 ± 100` | `+535 ± 124` | `−1317 ± 261` |

* So the κ-model leaves `≈ ±800` `c₂` at `x ≥ 0.4`, as SUBLEADING_REGGE found at `J = 40` ("κ fails at `x ≳ 0.45`, in the oscillating
  class").
* Three points cannot separate a smooth part from the oscillating class. The bulk smooth test uses the fixed-engine K₂ row
  `(1064, 64)` (`fw_queueB5.sh`), validated against these points (§3).

## 3. Bulk smooth term at `J = 64` (`Δ_g = 1000`): K₂ row with the fixed engine, against direct [NE, negative]

Setup:
* Fixed-engine K₂ row `(1064, 64)`, X35, EE + C3H, state + reference (`fw_rows.py row EE,C3H 64 1064 X35`).
* Generalised bases: state `sv3 = 2.5·10⁻⁶` on `x ∈ [0.15, 0.6]`, reference `2·10⁻¹²`.
* Run by `FW_J=64 FW_NODIRECT=1 FW48_VAR=k2 python3 fw_chain48.py` → `cache/fw_chain64_k2.txt`.

Smooth-offset fit (`Γ − κ`, oscillating class + `r`):

| key | `x ∈ [0.1, 0.45]` | `x ∈ [0.1, 0.7]` | residual on `[0.1, 0.45]` |
|---|---|---|---|
| C2 `++` | `−55 ± 47` | `+57 ± 107` | 25 % |
| C2 `+−` | `−50 ± 42` | `−306 ± 214` | 25 % |
| C3H `++` | `+572 ± 502` | `+410 ± 319` | 70 % |
| C3H `+−` | `−352 ± 312` | `+128 ± 324` | 70 % |

Shape-model residuals: `5 %` (C2 `++`), `28–41 %` (the other keys; max `|res| ≈ 6·10³`).

K₂ against direct at the bulk points (same pipeline, dir estimator, `Γ − κ`, C2 `++`):

| `x` | K₂ | direct |
|---|---|---|
| 0.2 | `+430` | `+297` |
| 0.4 | `−1398` | `−246` |
| 0.6 | `+3348` | `+589` |

C3H shows the same pattern: `−1407` vs `−367` at `x = 0.4`, `+3714` vs `+1197` at 0.6.

**Reading.**
* **K₂ is not converged at `(1064, 64)`, `J²/m = 3.9`, even with the fixed engine.** The disagreement with direct is
  `10³–3·10³` `c₂` at `x ≥ 0.4`, and the K₂ chain gives a meaningless min `λ/N_s = −4.87` (`max |dF_hh|/N_h = 8.4`). So
  EE_DG1000_LARGEJ's K₂ criterion (`J²/m ≲ 1.5–2.2`) stands. The engine fix moves `J = 48` (`2.2`) into the usable range, not `J = 64`.
* The bulk smooth term at `J = 64` cannot be tested at the uniform tolerance (`14` `c₂`) with these data.
  * The K₂ C2 fit (`−55 ± 47`) is not trustworthy, because K₂ is off by `≈ 10³` at `x = 0.4`.
  * The 3 direct points (`Γ − κ` within `±800` CV, `±600` dir at `x = 0.2–0.6`) cannot separate a constant from the oscillating
    class.
  * **Bulk smooth term at `J > 48`: still [I].** A direct row with `≳ 15` nodes on `x ∈ [0.1, 0.5]` would be needed: about 5 core-hours
    per row at `J = 64`, with per-node errors `40–200` (CV).

## 4. `J = 48` at `Δ_g = 1000`: recomputed rows, normal bases, no offset [NE + M]

Rows (`fw_rows.py row EE,C3H 48 1048 X35`, refb recomputed in the patched basis):
* All 35 nodes (EE and C3H), state and reference, recomputed with the fixed engine. Run by `FW48_VAR=k2 python3 fw_chain48.py` →
  `cache/fw_chain48_k2.txt`.
* Old versus new records (same keys) differ by a fixed `2×2` map, from the different structure basis and shadow, plus the engine
  noise. `fw_fix` reconciles 23 records with an older shadow.

Helicity bases (generalised):

| point | `x` window | `sv3/sv1` | `imag` | old l8 basis (`sv3` / `imag`) |
|---|---|---|---|---|
| `(1048, 48)` state | `[0.15, 0.6]`, 19 records | `1.2·10⁻⁷` | `9·10⁻¹¹` | `7·10⁻⁵` / `0.057` |
| `(2·10⁴, 48)` reference | `[0.3, 0.6]` | `1.1·10⁻⁷` | `2·10⁻¹¹` | `8.5·10⁻⁶` / `9·10⁻⁴` |

These are now as normal as `J ≤ 40` [NE].

Smooth-offset fit of SMOOTH_BOUND §3 on the new rows (`Γ − κ`, oscillating class + constant):

| key | `x ∈ [0.1, 0.45]` | `x ∈ [0.1, 0.7]` |
|---|---|---|
| C2 `++` | `+0 ± 5` | `+20 ± 34` |
| C2 `+−` | `+1 ± 5` | `+11 ± 18` |
| C3H `++` | `−1 ± 3` | `+19 ± 34` |
| C3H `+−` | `+0 ± 5` | `−8 ± 14` |

* The old rows gave `+280 … +650` / `−630 … −690`. **The `J = 48` offset was the engine's precision loss** (the contraction noise
  `≈ m²·10⁻⁴ ≈ 10²` `c₂`, plus the defective basis it produced). It was not a physical smooth term [NE].
* Shape-model residuals: `++` `1.7 %` / `0.8 %` (C2 / C3H); `+−` `23 %` / `25 %` (K = 4), `13 %` (K = 5), `2 %` (K = 6). Before
  the fix they were 20–33 % for every key.

Chain (ee1000 shape chain, hybrid states `m_gap … 5Λ`, sm_curvfix κ, rows 1048):
* Main (`K = 4`, `x ∈ [0.05, 0.8]`): **min `λ/N_s = +0.932`** (`1.2Λ`; uncorrected `+0.989`). Max `|dF_hh|/N_h = 0.065`,
  max `|dF_{+−}|/√(NN) = 0.010`.
* Variants:

  | variant | min `λ/N_s` | max `|dF_hh|/N_h` |
  |---|---|---|
  | `K = 5` | `+0.935` | |
  | `K = 6` | `+0.990` | `0.26` |
  | `SM_CURV = none` | `+0.930` | |
  | state basis `x ≤ 0.4` | `+0.932` | |
  | `x ≤ 0.6` | `+0.990` | 4.1 (held extrapolation) |
  | spline | `+0.730` | 2.7 (aliasing, as in EE_DG1000_LARGEJ §4) |

* **`n = 1` at `J = 48`, `Δ_g = 1000`: [NE + M]**, min `+0.73 … +0.99` over all variants. Old grade: [I].

## Hand-over (kept current)

**Running** (detached; resumable — after a container restart re-run the same commands from `ads4/code/spinfunc/`; finished points are
skipped):
* `setsid nohup ./fw_queueA7.sh > cache/fw_queueA7.log 2>&1 &`: `(1100, 100)` nodes at `x = 0.2, 0.1`, which anchor `a_h`.
* `setsid nohup ./fw_queueB6.sh > cache/fw_queueB6.log 2>&1 &`: its last item, the `(1200, 200)` probe at `ν = 16.7`, was still in its
  reference block at hand-over (22:00 UTC). It writes `cache/fw_pts_*_J200_D1200.pkl`; then run `python3 fw_probe.py 200:1200 16.7`.
  fw_probe sets `CUR[J]['D']`, so it reads the per-state file.
* Then analyse:
  * `FW_EST=cv FW_WFIT=1 python3 fw_gamma.py 100:1100 8 34` (fold fit; the `x = 0.1, 0.2` nodes enter `a_h` and are printed but are
    outside the fit window);
  * `python3 fw_probe.py 200:1200 16.7`.
* Superseded queues (do not re-run): `fw_queueA.sh`, `A2`, `A3`, `A4`, `A6` (done), `B`, `B3`, `B4`, `B5`.
* The two `ee1000_queueR1.sh` queues of the earlier agent were left alone.
* **Do not re-run** fw_gamma on `(400, 100)` or `(1064, 64)` and expect §2.2/§2.3: their shared-`ν` references were overwritten
  (§2.6).

**Analysis commands** (each ≤ a few minutes):

| command | output |
|---|---|
| `FW_EST=cv [FW_WFIT=1] python3 fw_gamma.py J:Delta 8 34` | fold-window `Γ`, `Γ − κ`, smooth fits (also `FW_EST=dir`); appends `cache/fw_gamma.txt` |
| `python3 fw_scalar.py J:Delta` | scalar check against sm_exact |
| `python3 fw_scalfit.py J:Delta nu1,...` | exact-scalar smooth-fit control |
| `python3 fw_probe.py J:Delta nu` | one-node precision |
| `FW48_VAR=k2 python3 fw_chain48.py` | `J = 48` chain, `cache/fw_chain48_k2.txt`; variants via `EE1000_K`, `EE1000_MODEL`, `EE1000_XMAX`, `SM_CURV`, `FW_XB1` |
| `FW_J=64 FW_NODIRECT=1 FW48_VAR=k2 python3 fw_chain48.py` | the `(1064, 64)` K₂ chain (not converged) |

**Caches** (`spinfunc/cache/`):
* `fw_pts_{EE,C3H}_J{64,100,150,200}.pkl`: direct points, z2/n1 format.
* `fw_k2_{EE,C3H}_J{48,64,100,150,200}.pkl`: K₂ records.
* `fw_shadowk2.pkl`: robust shadows in the patched basis; `fw_shadowrun1.pkl`: superseded first run, kept for `fw_fix`.
* `fw_sols.pkl`: structure-basis maps.
* `fw_raw.pkl`.
* Logs: `fw_queue*.log`.

**Code** (`fw_*.py`):
* `fw_frame`: Fourier floor, `Q_J` recurrence, conserved-points fix.
* `fw_shadow`: equilibrated shadow fit.
* `fw_fix`: load-time basis consistency.
* `fw_points` / `fw_rows`: direct points and K₂ records.
* `fw_gamma`: `Γ`, generalised basis, fits.
* `fw_scalar`, `fw_scalfit`, `fw_probe`, `fw_chain48`, `fw_rmrefb48`.

**For other agents' code**:
* Any graviton record at `J ≳ 44` made by l5_spin / l9_points / n1_pw without `fw_frame.patch_S` carries the `ε²4^J` error. This
  covers `ee1000_pts_*_J{48,56,64}`, `ee1000_n1_*_J{48,64}`, r1's `J = 72`, the `w1` 336-row `J = 32` (small, `~10⁻¹¹`), etc.
* Shadows at `J ≥ 56` from `l5_shadow_matrix` (unscaled fit) are unreliable.
* The conserved basis is SVD-arbitrary: never mix records computed with different `conserved()` sampling.

**Next steps (value order):**
1. `(1100, 100)`: finish and fit (queued).
2. `(364, 64)`: refit with 10 nodes (queued).
3. A method without the state-side sign problem for `J²/m ≫ 1`, needed for `Δ_g = 300`, `J ≥ 64`. For example, split `P_J` into its
   two travelling waves and give each its own (`J`-dependent) thimble.
4. A direct bulk row at `(1064, 64)` with `≳ 15` nodes (about 5 core-hours) for the bulk smooth term.
5. Re-check EE_DG1000_LARGEJ's `J = 56–64` conclusions with the fixed engine. K₂ at `(1064, 64)` is still not converged, but its
   numbers there were engine noise.
