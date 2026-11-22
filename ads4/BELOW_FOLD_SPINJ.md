# Below the fold at large spin: the closed-form spin-J EE engine, the `ν < 8.3` window at `(336, 28/32)`, and the high-J window

*21–22 Nov 2026, below-fold-spinJ agent. Unrefereed, not committed. No existing document or code file was edited; nothing of the
`oo_*` agent was touched or run. New code `ads4/code/spinfunc/bj_*.py` (+ `bj_queue*.sh`), caches `spinfunc/cache/bj_*`, `z2bj_*`,
`hc_chain_z4_bj*` (git-ignored; the last are outputs of the unchanged `hc_chain.py`). ≤ 2 processes, `OMP_NUM_THREADS=1`,
`ulimit -v 3 GB` on every job (one early unconstrained ring test was OOM-killed at ~13 GB; all later jobs ran under the limit).
Read-only imports: `bg_ser`, `bg_geom`, `bg_ee`, `bh_spinJ`, `bg_scalarJ`, `l5_spin`, `fw_frame`, `z3_nodes`, `hc_chain`/`z4_chain`/
`w2_chain` (run unchanged, in-process wrappers only), `l9_analysis`, `fw_gamma`.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model;
**[I]** open. `c₂ = m² ×` relative correction. "Engine" = the analytic-radial 1/Δ expansion on the real section D (BELOW_FOLD2); its
EE beyond-LO equivalence with the crystal cycle is assumed there — §2.2 below is the first test of it at large J.

## 0. Answer in brief

1. **Construction [P + DV].** The spin-J O-index contraction is now closed form inside the engine's series algebra. Every conserved
   `⟨TTO_J⟩` structure is `(V·z)^{J−4} × quartic`, so the contraction is a 9×9 kernel `B_{nn′}` between the quartics' null-circle Fourier
   modes, `B = N_J R^J Σ_{k≤8} P_J^{(k)}(t) B_k` with `B_k` from the BELOW_FOLD3 polarisation formula written in manifestly small
   variables (§1). `R^J` cancels the J-parts of the 3-pt exponents identically; `P_J^{(k)}` is resummed in J/Δ exactly as bg_scalarJ
   (D̂-series, all residual orders). Cost is J-independent apart from the D̂-series length; J = 28–100 run in 1–1.5 h per state at
   geometric depth 2 (relative order a¹), 3–5 h at depth 3.
2. **Validation (§2).** ring vs 60-digit sampled pairing ≤ 1·10⁻¹³ (J = 4–32) [DV]; ASer `B` vs exact numeric at the same point
   8·10⁻⁴ / 1·10⁻⁵ / 3·10⁻⁷ (K = 1/2/3, J = 28) [DV]; scalar J = 24–32 through the same P-series path = bg_scalarJ (itself = exact
   skeleton) to ≤ 0.01 c₂ (K = 4) [DV]; **EE at J = 28 against the independent z2 n3 direct nodes (ν = 8.3–20.8): engine − direct
   = −1374 ± 2756, −230 ± 496, +72 ± 254, −50 ± 286, +242 ± 195, +1.5 ± 27 c₂** (after one fitted constant; first test of "D ≡ crystal
   cycle" for EE beyond LO at large J) [NE]; reference (2·10⁴, J = 28) vs direct: entrywise 2·10⁻⁷ … 5·10⁻⁵ (ν = 8.3–33) [NE];
   J = 4 vs the sampled engine: agreement at the level of the sampled engine's own truncation (§2.1) [NE, not decisive].
   **Not achieved:** a convergent comparison with the K₂ rows above the fold at J ≥ 4 (the engine's D̂-series at N = 4–5 is not
   converged at ν ≳ 20 for J = 28) [I]. **Main limitation found [NE]:** the geometric 1/Δ series converges slowly at large `J²/Δ`
   (J = 28, Δ = 336: per-entry a¹ term 2–4·10³ c₂, a² term ≤ 250 c₂ on ν ≤ 5 but 370–730 at ν = 10–12; at (1100, 100) the a¹ term
   is ≈ 8·10³ c₂) — the J/Δ resummation of P_J does not resum these J-enhanced geometric terms.
3. **Result D, `ν < 8.3` at `Δ_g = 300` (§3) [NE + M].** Engine `Γ_c` below the fold (C2, complex, chain convention):
   `(336, 28)`: `++` `Re Γ_c ≈ +210 … +430` on ν = 1.5–6 (`+317 / +367` at ν = 3 / 5), `+−` `≈ −580` mean (`−700 / −800`);
   `(336, 32)`: `++` `+380 … +560`, `+−` `≈ −490` mean. **None is near the failure scale (uniform ≲ −4·10³ … −5·10³).**
   Fed into the unchanged hc_chain z4 protocol (engine C2 NL density below 8.3, everything else as HYPOTHESES_CLOSE §5):
   **min λ/N_s at the gap: J = 28 +0.787 (relative order a², depth 3; C3H thimble nodes kept) / +0.793 (C3H below 8.3 dropped);
   order a¹: +0.739 / +0.746; order a⁰: +0.590. J = 32 (order a¹ only): +0.834 at 1.3Λ (twist edge 1.102Λ: +0.974).** Calibration
   (one fitted constant per state) ±295 / ±540 c₂ → ±0.05 / ±0.06. Failure needs an additional uniform `Γ_c ≲ −4.3·10³ c₂` on top
   of the engine values (both J). Grade [NE + M]: the window is now computed, positive with margin ≈ 0.7, but the engine's own
   uncertainty there is O(10²–10³ c₂) (slow geometric convergence, calibration), not O(10).
4. **High J / H_LJ / H_win (§4) [I].** The engine runs at (1100, 100) (pieces + images, fw basis) and (364, 64) (pieces only, run
   stopped), but there it is not accurate enough for the H_win tolerances (T_B ≈ 1.4·10³ / 2.5·10²): at (1100, 100) the D̂-series is
   converged (N = 4 vs 5 ≤ 6 c₂) only for ν ≤ 10 (560 c₂ at 13.3, where the fw direct nodes start, so no calibration), and the
   geometric a¹ term alone is ≈ 8·10³ c₂ (J²/Δ = 9). No λ/N_s or Γ − κ is quoted for high J.
   **Verdict: the fold region at J ≥ 64 (300) / J ≥ 150 (1000) still requires H_win/H_LJ** [I]; what is missing is technical
   (a resummation of the J²/Δ-enhanced geometric terms, and N ≳ 8 within 3 GB, §5).

## 1. Construction
### 1.1 Quartic remainder and the 9×9 pairing kernel [P + DV]
* Every parity-even conserved `⟨T T O_J⟩` monomial has `a₅ = J − b₁₅ − b₂₅ ≥ J − 4`, so in the O polarisation
  `f_a(z₅) = (V·z₅)^{J−4} φ_a(z₅)` with a quartic `φ_a` (`V = V_{5,12}`; likewise `g_c = (W·z₅)^{J−4} γ_c`, `W = V_{5,34}`) [P].
  The quartic is exactly bg_ee.tto evaluated with "J = 4 powers" `V₅^{4−b₁₅−b₂₅}` and the conserved coefficients of the TRUE `(Δ, J)`
  (the monomial list of l4 is J-independent for J ≥ 4).
* On the null cone a quartic is fixed by its 9 null-circle Fourier modes, and `e_n(z) = (ζ_{sgn n}·z)^{|n|}(u·z)^{4−|n|}` with
  `ζ_± = (−1, 1, ±i)`, `u = (−1, −1, 0)` restricts to `e^{inφ}` (`ζ_±·z = e^{±iφ}`, `u·z = 1`). Hence
  `⟨f_a, g_c⟩ = Σ_{n,n′} φ_{a,n} γ_{c,n′} B_{nn′}`, `B_{nn′} = ⟨(V·z)^{J−4}e_n, (W·z)^{J−4}e_{n′}⟩` [P], and at J = 4 `B = δ_{n+n′,0}/κ_n`
  (the old contraction).
* `B` from BELOW_FOLD3 §3 (Legendre identity + polarisation): with `R = √(V²W²)`, `t = V·W/R`,
  `B_{nn′} = N_J R^J Σ_{k=0}^{8} P_J^{(k)}(t) B_k[n,n′]`, `B_k` the Taylor coefficient `[ε_s^{|n|}ε_u^{4−|n|}η_{s′}^{|n′|}η_u^{4−|n′|}]` of
  `(1+s_a)^J(1+ρ_a)^{J/2}(1+s_b)^J(1+ρ_b)^{J/2} τ^k/k!` with `s_a = v·δa`, `ρ_a = (δa² − s_a²)/(1+s_a)²`,
  `τ = t[(1+ρ_a)^{−½}(1+ρ_b)^{−½} − 1] + (q_a + q_b + δa·δb − t s_a s_b)/[(1+s_a)(1+s_b)√((1+ρ_a)(1+ρ_b))]`, `q_a = w_⊥·δa`,
  `w_⊥ = (W − (V·W/V²)V)/|W|` (manifestly O(1/Δ)) [D]. Checked against the 60-digit sampled pairing for J = 4–32, random complex
  kinematics: **≤ 1·10⁻¹³ relative** (`cache/bj_ring_test.txt`) [DV].
* `R^J = (V²W²)^{J/2}` cancels the J-dependent parts of the 3-pt exponents identically (`V² = −m₁₅m₂₅/m₁₂`, `W² = −m₃₅m₄₅/m₃₄`;
  asserted in code) [P], so the logs are the J = 0 ones; `N_J` is J-only and drops out of every state/reference ratio.

### 1.2 a-expansion bookkeeping [D + DV]
* First attempt (BELOW_FOLD3 §3.2's suggestion, J counted as `a⁻¹`) is **wrong at J = O(10–30)**: the cross terms `δa·δb ~ 1/J²`
  are then counted as `a²` and truncated, but numerically they are O(1/J²) ≫ a² (at J = 6 the `B_{4,−4}` entry, 14× the `B_{00}` one,
  was dropped). Kept: **every J-dependent number exact**; only geometric a-orders are truncated, relative to O(1) leading objects. For
  that the generators are rescaled by 1/J (`s_a → S/J`, `δa², δa·δb → /J²`) and the manifestly small `q_a = a·q̂_a` enters as `q̂_a`
  (d = 0) times the number `a/J`. Then every term `P_J^{(k)}τ^k ~ J^{2k}(a/J)^f J^{−2s} = (Ja)^f` is O(1) at fixed J/Δ with a uniformly
  truncated geometric factor [D]. Check: ASer `B` vs the exact numeric closed form at the same geometric point (V, W evaluated),
  J = 28, Δ = 336: `max|ΔB|/max|B| = 8·10⁻⁴ / 1.0·10⁻⁵ / 3.4·10⁻⁷` at K = 1 / 2 / 3 [DV].
* J/Δ resummation exactly as bg_scalarJ: `t = 1 + a²D̂`, `P_J^{(k)}(t) = Σ_n c_{n+k}(J)(n+k)!/n! a^{2n} D̂^n` (exact numbers),
  `pref_J/(N_J R^J) = Σ_{n ≤ N} D̂^n Σ_k γ_{k,n}X_k`, full residual order for every n (a triangular truncation `r_max(n) = K − n` was
  tested and is NOT adequate: it costs 10–10³ c₂, the same lesson as BELOW_FOLD2's R = 2).
* Contraction `X_k[a][c] = Σ_j κ₂(j)⁻¹ Σ_{n,n′} Em_a[j][n] B_k[n,n′] G_c[4−j][n′]` (`Em = Σ_i D[i][j]F_a[4−i][n]/κ₂(i)`, `D = E²`).
  Entries `(k, n, n′)` whose contribution `|P^{(k)}(t)B_k F_n G_{n′}|` stays below `10⁻⁶` of the largest at sample points (radial values
  ×1 and ×2.5) are dropped (J = 28: 157 of 402 entries kept at 10⁻⁶, 402 at 10⁻¹⁰).

### 1.3 Two precision/perf items found on the way
* **`bg_ser.Poly.prune` precision bug** [DV]: it drops coefficients below 10⁻¹⁵ of the largest coefficient. With high powers of D̂
  (degree ≈ 2n in radial atoms whose values in the integration region are O(10)) such coefficients carry O(1) values: `D̂⁵·10⁻¹⁶` added to
  an O(1) polynomial changed its value by 10⁻⁴. In every bj_* process the prune is value-weighted (a term is dropped only if
  `|c·monomial| < 10⁻¹⁶ … 10⁻¹²` of the largest at either of two reference points, radial atoms 1 and 4). The old engine (J ≤ 4) does not
  reach such powers; no existing result is affected as far as checked (bg_scalarJ forms moments per n, no mixing of scales).
* Ring arithmetic on bg_ser.ASer costs ~2 ms per product of small objects (86k products ≈ 200 s at K = 2 per B); `bj_xs` (one
  table per ring element, ε-graded, vectorised outer products, chunked combine) does the same in 1–7 min at K = 2–3. The contraction
  uses a batched multiply-accumulate (one combine per output order; dense-box bincount when the exponent box is ≤ 1.2·10⁷ cells).

## 2. Validation

### 2.1 Table
| check | result | grade |
|---|---|---|
| closed-form 9×9 kernel vs 60-digit sampled pairing, J = 4, 6, 8, 12, 24, 32, random complex V, W, 6 (n, n′) each | ≤ 1·10⁻¹³ relative (exact zeros at J = 4 to 10⁻¹⁵ absolute) | [DV] |
| ASer/XS `B` (bookkeeping §1.2) vs exact numeric closed form at the same point, J = 28, Δ = 336 | 8·10⁻⁴ / 1.0·10⁻⁵ / 3.4·10⁻⁷ at K = 1 / 2 / 3 (also 5.8·10⁻⁴ / 2.6·10⁻⁶ at J = 6) | [DV] |
| scalar through the bj P-series path (pref ≔ 1, t from the EE vectors) vs bg_scalarJ N = 12, R = 4 | (336, 24/28/32) Kr = 4: ≤ 0.000 / 0.000 / 0.008 c₂ (ν ≤ 18); Kr = 3: ≤ 0.1 / 0.4 / 3.8; Kr = 2: 3 / 9 / 80 | [DV] |
| J = 4, Δ = 300, pref pointwise vs bg_ee (sampled), depth 2 | rel. diff 2·10⁻⁴ … 1.8·10⁻³; the sampled engine itself moves by 6·10⁻⁴ … 3.4·10⁻³ from depth 2 to 3 at the same points → consistent, not decisive (equal-depth-3 check not run) | [NE, weak] |
| EE J = 28: engine (depth 2, N = 5) vs z2 n3 direct nodes, state/ref ratio, CV estimator, after one constant | §0.2 numbers; χ² ≈ 3.4 / 6 (8.3–20.8) | [NE] |
| EE J = 28 reference Δ = 2·10⁴ vs direct | entrywise 1.6·10⁻⁷ (8.3) … 5·10⁻⁵ (33.3) | [NE] |
| leading-Regge large-J limit | the a⁰ image at Δ = 2·10⁴ (leading Regge) reproduces the direct reference to ≤ 5·10⁻⁵ at J = 28 (above); no separate kernel comparison | [NE] |
| K₂ rows above the fold (l9 pred form), J = 28 | not conclusive: engine N = 5 unconverged at ν ≥ 20 (N = 4 vs 5: 10³–10⁴ c₂ at 25) | [I] |
| internal: D̂-series N = 4 vs 5 (J = 28 / 32) | ≤ 0.2 / 1.2 c₂ on ν ≤ 10; 12 / 5 at 13.3 (J = 32 has a zero of entry (0,0) there) | [NE] |
| internal: geometric order a¹ (Kr = 0 → 1) | 10³–4·10³ c₂ per entry; at the chain level the a¹ term moves `Γ_c` by 10³–2·10³ and λ(gap) by +0.15 (J = 28) | [NE] |

### 2.2 What the EE check at J = 28 means
The overlap nodes are independent crystal integrals (n3 contour, CRN reference). After one complex constant per state the engine's
state/reference ratio follows them within errors from ν = 8.3 to 20.8, including the most precise node (20.8: `+1.5 ± 27 c₂`). This is
the first evidence that the real-section engine equals the crystal cycle for the EE beyond leading order at large J (BELOW_FOLD2 §5.4
assumed it) [NE]. At ν = 25 it fails by `−1300 ± 27` — the D̂-series (N = 5) is not converged there (N = 4 → 5 changes it by 10⁴),
not a contour discrepancy.

## 3. Result D: the `ν < 8.3` window at `(336, 28)`, `(336, 32)`, `Δ_g = 300`

### 3.1 Protocol (`bj_chain.py`) [D]
* Engine state (336, J) and reference (2·10⁴, J), depth 2 (relative a¹), D̂-series N = 5 (state) / 4 (reference), kind EE2 = C2,
  conserved basis of `l5_spin.conserved` (= the z2 records' basis).
* Calibration: per common-seed replica, `κ_st`, `κ_ref` = least-squares constants of the engine 2×2 matrices against the z2 records
  (`reps·e^{logN}`) at ν = 8.33, 10.0, 11.66 (CV estimator: the replica's MC scalar ratio replaced by the exact one); `ρ = κ_st/κ_ref`
  jackknifed: ±295 c₂ (J = 28), ±540 c₂ (J = 32). The constant is a normalisation, i.e. a `c₀ D_h`-shaped Γ shift.
* Records at ν′ = x·m (x = round(ν/m, 6), 0.25 grid), replica-wise calibration carried. Pass 1: these records + C3H copies (only to
  keep the nodes in q3's rows; C2 depends on the EE kind alone, `a_h` on ν ≥ 15, bases without injected records) → `z4_chain.NLROWS`
  `Γ_c`, `Γ_NL` per helicity. Pass 2: the unchanged chain with its C3H as in the main run (or C3H below 8.3 dropped) and the k = 2 NL
  density below 8.3 replaced by the engine's `Γ_NL` (cubic spline in ν). Everything above 8.3 (direct n3 nodes, hc n1 bulk, K₂ rows
  500/1600 or 1040) is the HYPOTHESES_CLOSE §5 protocol.

### 3.2 Engine values below the fold (`Re Γ_c`, c₂; complex values in `cache/bj_nl_J*_bj*d2.pkl`)
| ν | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 | 5.0 | 6.0 | 7.0 | 8.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| J = 28 `++` | +418 | +564 | +213 | +288 | +317 | +335 | +367 | +493 | −206 | +159 |
| J = 28 `+−` | −2603 | −5961 | +1096 | −283 | −700 | −770 | −800 | −964 | −100 | −543 |
| J = 32 `++` | +629 | +827 | +378 | +467 | +491 | +496 | +508 | +557 | −95 | +264 |
| J = 32 `+−` | −2767 | −6324 | +1328 | −153 | −590 | −657 | −696 | −974 | +251 | −381 |
(ν ≈ 6.3–7.0: the local `Γ_c` passes the zero of `h₋ + h₊` at 6.64; the chain uses `Γ_NL = a_χ Γ_c`, finite there. ν ≲ 1.3: the `+−`
channel swings through 10³–10⁴ where its forward normalisation is small; the window weight there is small (HYPOTHESES_CLOSE2: 98 % at
ν < 6, 58 % at 3–5).) For comparison: the pooled thimble node of HYPOTHESES_CLOSE2 at ν = 5.0, `Re Γ_c = −1300 ± 980`, vs the engine
`++ +367 / +− −800` at J = 28 — within 1.7σ / 0.5σ.

### 3.3 λ/N_s (hc_chain z4, `cache/hc_chain_z4_bj*.txt`)
| J | variant | gap state | min over m [m/Λ] |
|---|---|---|---|
| 28 | baseline (main thimble run, held below 8.3) | +1.728 (1.091) | +0.835 [1.5] |
| 28 | **engine below 8.3 (order a², depth 3, N = 4), C3H kept** | **+0.787** (1.091) | **+0.785** [1.3] |
| 28 | same, C3H below 8.3 dropped | +0.793 | +0.787 [1.3] |
| 28 | engine order a¹ (depth 2, N = 5), C3H kept | +0.739 (1.091) | +0.739 [gap] |
| 28 | engine order a¹, C3H below 8.3 dropped | +0.746 | +0.746 [gap] |
| 28 | engine also replacing the 8.33 node (cut 8.4), C3H kept / dropped | +0.741 / +0.806 | +0.741 / +0.790 [1.3] |
| 28 | order a¹, calibration ±1σ (uniform ±295) | +0.789 / +0.690 | +0.786 / +0.690 |
| 28 | geometric order a⁰ only (Kr = 0) — for scale | +0.590 | +0.590 |
| 32 | baseline (held, oowin x_oo2) | +0.884 (1.102) | +0.825 [1.15] |
| 32 | **engine below 8.3** | **+0.974** (1.102) | **+0.834** [1.3] |
| 32 | calibration ±1σ (uniform ±540) | +1.097 / +0.852 | +0.877 / +0.771 |
| 28 / 32 | stress: engine + uniform −2000 / −4000 / −8000 | | +0.411 / +0.084 / −0.572 ; +0.521 / +0.068 / −0.838 |

**Reading.** With direct values everywhere below the fold the window is no longer a held extrapolation: both states pass with
λ/N_s ≈ +0.74 (J = 28) and +0.83 (J = 32), the failure scale is an *additional* uniform `Γ_c ≈ −4.3·10³ c₂`, and the measured engine
values are O(10²–10³) with both signs. Grade **[NE + M]** for the engine values (M = depth 2 (a¹) truncation and the calibration
constant), conditional on the C3H window (below 8.3 still the thimble/n3 C3H nodes or nothing, as in the main run; the engine has no
C3H integrand). **Order by order (J = 28):** λ(gap) = +0.590 (a⁰) → +0.739 (a¹) → +0.787 (a²). Per entry the a² term is ≤ 250 c₂
on ν ≤ 5 (where 98 % of the window weight sits) and 370–730 c₂ at the calibration nodes ν = 10–11.7, so most of its effect on λ comes
through the calibration constant (a c₀D_h-type shift). The a² Γ_c table: `++` +1016 / +1052 / +903 / +911 / +854 / +760 / +612 / +247 /
+1291 / +680, `+−` −1817 / −5461 / +2160 / +607 / −25 / −390 / −927 / −2615 / +3917 / +351 at the ν of §3.2. Against the direct nodes
above the fold the a² engine is not closer than the a¹ one (13.3: +705 ± 288; 16.7: −581 ± 194; 20.8: +1692 ± 27, the last
N-unconverged) — the geometric series is converging slowly here, and the engine's window values are uncertain at the 10²–10³ c₂
level [NE]. J = 32 was run at order a¹ only.

## 4. High J (H_LJ / H_win)
* **What was run.** `(1100, 100)`, Δ_g = 1000 (J = 100, Δ = 1100, `m = 1093.9`, H_win window ν < 16.3 (c = 3)): pieces at depth 2
  (95 min), images N = 4, 5; reference (2·10⁴, 100); fw (patched) conserved basis, so that the records sit in the fw rows' basis.
  `(364, 64)`, Δ_g = 300: pieces at depth 2 (stopped at k = 5 to end the session); reference (2·10⁴, 64) done. `(1150, 150)`,
  `(1300, 300)`, `(334, 34)`, `(348, 48)`: not run (no time; also no basis rows for Δ = 1150, 1300, 334, 348 — the chain would need them).
* **Convergence at (1100, 100)** (`cache/bj_conv.txt`): D̂-series N = 4 vs 5: ≤ 1.5 c₂ on ν ≤ 8, 6 at 10, 26 at 11.7, **560 at 13.3 and
  16.7**; geometric order a⁰ → a¹: 8–1.8·10⁴ c₂ per entry over the whole window (J²/Δ = 9.1, vs 2.3 at (336, 28) where a¹ ≈ 2–4·10³).
* **fw_gamma conventions** (`bj_gamma.py`, `cache/bj_gamma_1100.log`): engine nodes ν ≤ 13 injected into the fw (1100, 100) row with
  the fw direct nodes ≥ 13.3 (C3H copies only to keep the nodes; C2 only), calibrated on the 13.3 node: the engine's `Γ − κ` runs from
  −7·10⁴ (ν = 1) to −2·10⁴ (ν = 12) and jumps to +1262 at the direct 13.3 node — i.e. the calibration node is where the engine is not
  converged, and the engine values are not usable in Γ units. (Also: the fw local de-shift `Γ` is singular at the zero of `h₋ + h₊`,
  ν ≈ 6.64; the z4 NL convention would be needed in the window.) The CV estimator cannot be used for engine records (no MC scalar);
  `FW_EST=dir` was used.
* **Verdict for H_LJ / H_win.** The closed-form engine removes the 2J+1-sample bottleneck and computes the window below the fold at
  any J, but at J²/Δ ≳ 2 its 1/Δ geometric series converges too slowly to resolve `|Γ − κ| ≲ T_B` (1.4·10³ at (1100, 100), 2.5·10² at
  (364, 64)): the a¹ term alone is 2·10³–2·10⁴ c₂. **The fold region still requires H_LJ / H_win at J ≥ 64 (300) and J ≳ 100–150
  (1000)** [I]. The one large-J point where it does decide something is (336, 28/32) (§3), because there the failure scale (≈ −4·10³)
  is far above the engine's uncertainty.

## 5. Remaining open items [I]
1. **J²/Δ-enhanced geometric terms.** The a-expansion at fixed J/Δ is uniform only for the P_J factor; the EE structure (the q, s
   attachments of §1.2 and the conserved coefficients) carries terms growing like J²/Δ per order (a¹ ≈ 2·10³ c₂ at J = 28 vs ≈ 40 at
   J = 0). A resummation in J²/Δ (e.g. expanding around the exact frame of V, W rather than around the forward point) is the next step;
   without it the engine is a 10²–10³ c₂ tool at J ≈ 30 and a 10⁴ c₂ tool at J ≈ 100.
2. **N ≳ 8 at depth ≥ 2 within 3 GB**: pref_sum + assembly of D̂^n·Y_n grow to 10⁶ terms per order; per-n integration (linear) instead
   of summing the polynomial first would bound the memory. Needed at (1100, 100) above ν ≈ 11 and at (364, 64) above ν ≈ 10.
3. **C3H window** at (336, 28/32): the engine has no C3H (o_H q) integrand; below 8.3 the chain keeps the n3/thimble C3H nodes or
   drops them (both variants quoted, they agree to 0.01–0.05).
4. **Calibration constant**: one fitted complex constant per state (±295 / ±540 c₂); an analytic normalisation (the Δ-dependent constants
   of the n3/w0 integrand vs the engine's dropped powers) would remove it.
5. **J = 4 equal-depth check** against the sampled engine (depth 3 both) and the image-level comparison with `bh_eeq_D304.03_J4_K4`
   (depth 5) were not run; the K₂-row check at J ≥ 4 above the fold needs N-converged engine values at ν ≥ 20.
6. **J = 32 at order a²**, (364, 64) images, and the high-J states without rows — not run.
7. Inherited: the real-section ≡ crystal-cycle equivalence for the EE beyond LO is now supported at J = 28 on ν = 8–21 (§2.2) [NE], not
   proved; the chain model [M] (HYPOTHESES_CLOSE).

Hand-over: nothing is running. To resume: `BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 python3 bj_eeq.py 336 32 3 4 win` (≈ 4 h; pieces are
cached once computed), then `python3 bj_chain.py calib 32 4 2 3`, `HC_TAG=bj32d3 python3 bj_chain.py nl 32 4 2 3`,
`HC_TAG=bj32d3_win_keep python3 bj_chain.py win 32 bj32d3` with the J = 32 hc/z4 environment of §3 (`HC_INJ=32:336 Z4_JS=32
Q3_JEVEN=32 W2_ROWS=336,1040 Q3_EXTRA_M=1.102,1.1,…,5 HC_IXMIN=0.12 HC_IXMAX=0.81 HC_DESPIKE_INJ=0 Z4_EST=best Q3_JODD=5
SM_CURV=row`).

## Graded claims
| claim | grade |
|---|---|
| conserved ⟨TTO_J⟩ = (V·z)^{J−4}×quartic; contraction = 9×9 kernel in the quartic's Fourier modes; R^J cancels the J-parts of the 3-pt exponents | [P] |
| closed-form kernel (polarised Legendre identity, manifest-small variables) = sampled pairing to 10⁻¹³, J = 4–32 | [DV] |
| bookkeeping: exact J numbers, geometric truncation relative to O(1) objects; ASer B → exact as K grows (8·10⁻⁴ → 3·10⁻⁷) | [D] + [DV] |
| "J ~ a⁻¹" bookkeeping (BELOW_FOLD3 §3.2) is wrong at J = 6–30 (drops O(1/J²) cross terms) | [DV] |
| bg_ser prune drops O(1)-valued terms of high-degree polynomials (fixed process-locally in bj_*) | [DV] |
| scalar through the bj P-series path = bg_scalarJ to ≤ 0.01 c₂ (J = 24–32) | [DV] |
| EE engine at J = 28 agrees with independent n3 direct nodes on ν = 8.3–20.8 within errors (one constant) | [NE] |
| (336, 28/32) window: engine Γ_c = O(10²–10³) c₂ with both signs; λ/N_s = +0.787 (28, a²) / +0.834 (32, a¹); failure needs ≈ −4.3·10³ extra | [NE + M] |
| engine accuracy at (336, 28) in the window: O(10²–10³) c₂ (slow geometric convergence + calibration) | [NE] |
| high-J (J²/Δ ≳ 2) window not resolvable at the T_B level with the present engine; H_LJ/H_win still needed | [NE] / [I] |

## 6. Files (`ads4/code/spinfunc/`; caches in `cache/`)
| file | content |
|---|---|
| `bj_ring.py` | closed-form 9×9 pairing kernel: numeric version + test vs `bh_spinJ` 60-digit sampling (`cache/bj_ring_test.txt`); production `bk_matrices_xs` (ASer, bookkeeping §1.2), `_test_xs` |
| `bj_xs.py` | ε-graded vectorised series tables (ring arithmetic), batched multiply-accumulate `macc` / `macc_dense`, low-memory combine |
| `bj_ee.py` | J-generic EE pieces (quartic F, G; D; X_k), value-weighted prune + low-memory `Poly.__mul__` (process-local patches), P-series `pref_sum`, `images`, `value` |
| `bj_eeq.py` | production: pieces → `cache/bj_pieces_D*_J*_d*.pkl`, images for several (N, Kr) → `cache/bj_eeq_D*_J*_d*[_p].pkl` (`_p`: fw patched basis) |
| `bj_valscal.py` | scalar check vs bg_scalarJ (`cache/bj_valscal_keep.txt`) |
| `bj_conv.py` | internal N / Kr convergence (`cache/bj_conv.txt`) |
| `bj_chain.py` | calibration vs z2 records, K₂ check, records, two-pass z4 chain with engine `Γ_NL` below 8.3 (`cache/bj_calib_J*`, `bj_nl_J*`, `hc_chain_z4_bj*`) |
| `bj_gamma.py` | engine nodes in fw_gamma conventions (high J) |
| `bj_sum.py` | per-state λ/N_s tables from the chain pickles |
| `bj_queue.sh`, `bj_queue2.sh` | the run queues |

