# EDGE_R: the graviton (EDGE-R) quantity of Lemma P4W-χ — not closed; what blocks it, a method that half-works, the exact scalar at the binding states, and the external job list

*24 Nov 2026. Unreviewed, not committed. No existing file was edited. One new file outside the code directory: `ads4/ext/jobs_batch3.txt`.*
* *New code:* `ads4/code/edger/` (`run_all.sh`, ≈ 1.9 CPU-h, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`).
  * `er_guard.py`: process-local write guard. Every write must go under `edger/cache/`. It registers itself as `ng_guard` and `hs_guard`.
  * `er_core.py`: E_χ from real-axis Γ data, the units H, S_E and S_E,t, and the cutoffs.
  * `er_diag.py`, `er_diag2.py`, `er_cut.py`: why n1 fails (§2.1).
  * `er_n3test.py`, `er_sval.py`: the n3 contour on the scalar (§2.2–2.3).
  * `er_points.py`: `spinfunc/hc_points.py`, unchanged, under the guard (§2.4).
  * `er_n3edge.py`, `er_row.py`, `er_existing.py`: the J = 29 row (§3).
  * `er_scalar.py`: the exact scalar (§4).
  * `er_jobs.py`: the job list (§6).
  * `er_existing.py` is the first version of `er_row.py` (n1 data only). It runs for J = 29 only; for J = 21, 23, 27 the direct bulk rows are in the EXT pickles, not in `highspin/cache`.
* *Caches:* `edger/cache/`, whose `.gitignore` ignores everything except itself. The merged point files are copies in `cache/merged/`.
* *Read-only inputs:* `highspin/` (hs_core, hs_odd, hs_points.parse_x, direct points), `agproof/` (ag_exact, ag_flat, chi.json),
  `p4win/`, `lsw/` (floors), `spinfunc/` (n1, n3, z2, hc_points), `exactskel/`.
* *Interruption:* a container restart killed the first runs of the n3 pilot and of `er_sval` midway. Finished points survived. The
  remaining x values were rerun (`cache/sval_J29_part{1,2}.out`).

Grades as in AG_PROOF / UNIFIED_THEOREM §0. Units as in AG_PROOF §5:
* S_E is the k = 2 + 3 edge cost of a flat-profile wave pair of amplitude e_sc(1).
* S_E,t = S_E·(1 − x₁²)⁻² is its +− transferred version.
* S_E2 is the k = 2 part, used for the scalar.
* Γ is in c₂ units, the EE convention (κ subtracted).

---

## 0. Verdict

**(EDGE-R) for the graviton is NOT closed. Its grade stays [M] at every Theorem W state with J ≤ 51. Theorem W's sampled J\* is unchanged: J\* = 26, still conditional on (EDGE-R).** No state was found to violate it; no state was shown to satisfy it.

| # | finding | grade |
|---|---|---|
| **1** | **At the gap-mass states (m ≈ 1.03–1.05Λ, x₁ = 0.953–0.975), (EDGE-R) is decided on x ∈ [0.94, x₁].** For the exact scalar at J = 26–48 the part of E_χ from p > p_c (ν > 0.94Δ) is **0.4–1.9 S_E2**. It cancels the part from p ≤ p_c to a net 0.11–0.29 S_E2. Sign-free covers of that window cost 4–10× the transferred tolerance, so values with their oscillation are needed there. | [DV] (exact scalar, §4) |
| **2** | **Why the direct n1 method fails at ν/Δ ≳ 0.94.** The spin-J prefactor (P_J, and the graviton's spin-J Fourier projection) grows like e^{J·Im γ} at complex angle. On the real ξ-contour this produces a second hump about 11σ out along the soft Hessian direction. The hump is e^{14} above the saddle at x = 0.96 (J = 29), and the moment-matched proposal jumps to it (\|μ\| = 5.9, 94 % error). Plain proposals have heavy tails (14–345 %). Radial cuts have no plateau: already at x = 0.9 the soft direction contributes out to \|ξ\| ≈ 5. **Re-centring or cutting the n1 proposal cannot fix this.** | [DV] (§2.1) |
| **3** | **The n3 hybrid contour (Φ-thimble flow beyond ρ₁ = 3.5) works for the scalar up to the edge.** It is validated against the exact scalar (Arb): the double ratio T(ν) = [Y(1029)/Y(1100)]/[ex₂(1029)/ex₂(1100)] is constant to within 1.6σ over x = 0.88–0.972 (ν/Δ up to 0.970), with per-node errors of 1–9 %. | [DV] (§2.3) |
| **4** | **n3 for the graviton: the EE (k = 2) kind works; the C3H (k = 3) kind fails at the state for both parities and every large ν tested.** EE agrees with n1 where both work: J = 26, x = 0.93, all four entries within 1 % (0.4 % on the dominant one); J = 29, x = 0.9365, Γ_C2 = 3.6·10⁴ vs 3.1·10⁴ c₂, within the ≈ 10⁴ c₂ MC error. At the n1-failed node x = 0.9605 it gives Γ_C2 = −7.8·10⁴ c₂ (3.2 % RQMC), consistent with the (1−x)⁻² growth of the valid nodes. **C3H state values come out ×200–1000 too small, with 30–600 % errors**, already at x = 0.93 where n1 is fine. The Δ_ref = 2·10⁴ reference values agree with n1. This is a defect of n3's C3H evaluator at large ν, not an edge effect. | [DV]; the C3H defect is [I] (§2.4) |
| **5** | **J = 29 twist edge (odd, not a W state), existing n1 data.** The measured part (x ≤ 0.937) is cost_E = **0.27 S_E,t**. Error budget: MC 3σ **0.52**, bulk representation 0.34, reference-zero GK node (x = 0.8007) sign-free 1.24, four failed edge nodes sign-free 1.72. Upper bound **6.0 S_E,t**, against this state's own tolerance ≈ 1.1 S_E,t (odd-J H). **Not closable on existing data.** The MC term alone exceeds the even-J gap tolerances of 0.41–0.63. | [NE] value; budget [DV] (§3) |
| **6** | **Exact scalar at the binding W states (new; AG_PROOF had none at J = 26–38).** At the W state's own cutoff, \|E_χ[G]\| = **0.285 / 0.201 / 0.108 / 0.151 S_E2** at J = 26 / 30 / 34 / 38. The maxima over AG_PROOF's cutoff grid are 0.405 / 0.210 / 0.298 / 0.158. **This corrects AG_PROOF's "exact scalar ≤ 0.30 S_E for every cutoff"**: the grid maximum at J = 26 is 0.405, at the tolerance 0.406. At the tabulated cutoff the scalar-like reading passes with margin 0.406/0.285 = 1.42 (AG_PROOF: 1.35). J = 40 and 48 reproduce AG_PROOF (0.166, 0.143). | [DV] (§4) |
| **7** | **Extrapolating from x ≤ 0.94 is not a controlled route.** A shape model (c D + s xD′/(2J+1) + r, polynomials of degree K) fitted on exact scalar values at x ≤ x_c and continued to x₁ misses E_χ by 0.000 / 0.012 / 0.25 / 0.44 / 0.23 / 0.10 S_E2 (K = 3; J = 26 … 48), i.e. up to 3× the value itself. K = 1, 2 are worse (0.23–1.12). | [DV] (§4) |
| **8** | **External route (main route): `ext/jobs_batch3.txt`.** 1586 n1 jobs and 150 n3 jobs. Measured cost: n1 146 s per node (J = 26, M = 11), n3 490–630 s per node (M = 10). Total ≈ 64 + 23 ≈ **90 CPU-h ≈ 4–5 h on a 24-core machine (22 workers)**. The P3 group (the 20 states with m ≥ 1.13Λ, x₁ ≤ 0.887) is n1-only, needs no edge method, and has tolerances ≥ 7.4 S_E,t. Its ≈ 1220 jobs (≈ 50 CPU-h) can move those states to [NE±]. **For the 9 gap-mass states the k = 3 edge stays blocked by finding 4** until n3's C3H evaluator is fixed. | — (§6) |

### Per-state status (Theorem W, J ≤ 51)

The tolerance is AG_PROOF's, in S_E,t (transferred) units; the plain tolerance is ≈ the same absolute number (§1). "Scalar-like" means the exact scalar's own-cutoff value in its own units (k = 2).

| state (J, m/Λ) | x₁ | tolerance (tr.) | exact scalar, scalar-like (own cutoff / grid max) | graviton E_χ | grade of (EDGE-R) |
|---|---|---|---|---|---|
| **(26, 1.026)** | 0.975 | **0.406** | 0.285 / **0.405** | not computed (no data; x > 0.94 needs n3, C3H broken) | **[M]** |
| **(30, 1.030)** | 0.971 | **0.629** | 0.201 / 0.210 | not computed | **[M]** |
| **(34, 1.033)** | 0.968 | **0.886** | 0.108 / 0.298 | not computed | **[M]** |
| (38, 1.037) | 0.964 | 1.22 | 0.151 / 0.158 | not computed | [M] |
| (42, 1.041), (46, 1.045), (50, 1.049) | 0.960–0.953 | 1.59 / 2.05 / 2.57 | (J = 40: 0.166; J = 48: 0.129) | not computed | [M] |
| (33, 1.033), (51, 1.050) (odd) | 0.968 / 0.952 | 1.77 / 5.89 | — | not computed | [M] |
| all 20 states with m ≥ 1.13Λ | ≤ 0.887 | ≥ 7.45 | — | not computed; **n1 alone suffices**: P3 jobs | [M] → [NE±] expected after batch 3 (not yet) |
| (J = 29, 1.029) — not a W state | 0.972 | ≈ 1.1 (own odd H) | — | 0.27 measured (x ≤ 0.937) ± 0.52 (3σ MC); upper 6.0 with sign-free edge | — |

**Theorem W's J\*:** unchanged.
* Plain envelope: sampled J\* = 26 given (P0, P1′, W, AG_U, EDGE-R). (EDGE-R) is [M] at every J ≤ 51 state.
* +− transfer: J\* = 26 if the graviton is scalar-like. The exact scalar now has margin 1.42 at (26, 1.026) at the tabulated cutoff, but it touches the tolerance (0.405 vs 0.406) at another grid cutoff.
* No graviton evidence was added at the near-gap edge. External compute is needed.

---

## 1. Units and identification (`er_core.py`) [D]

The EE pipeline's correction of a state is ∫ w_k(x)Γ_k(x)dx, with w_k = ψ_k(p)·env_{h,k}(x)[u]/m²·dp/dx (`hs_core.state_weights`). P4W's k = 2 cost is 2∫ψ₂G. So for the entry h:
* k = 2: G_e = env_{h,2}Γ_{C2,h}/(2m²);
* k = 3: G_{e,3} = env_{h,3}Γ_{C3H,h}/m² (cost u∫ψ₃G₃);
* cost_E = 2\|E_χ,2\| + u\|E_χ,3\| = \|∫χw₂Γ₂\| + \|∫χw₃Γ₃\|.

Both pipelines use the same ψ_k (checked to round-off) and the same leading-image normalisation. P4W's H = f₊₊f₊₋/(f₊₊ + f₊₋) for even J and f₊₋ for odd J, with f = N_flat/5 + min(0, curv) (lsw floors). The EE base equals N_flat to 10⁻³, and H ≈ N_s/10 (even) or N_s/5 (odd). So cost_E compares directly with AG_PROOF's tolerance (1 − d_χ)H.

A consequence: **the plain tolerance (≥ 158 S_E) and the transferred tolerance (≥ 0.406 S_E,t) are the same absolute number** to within d_χ. For example, at (26, 1.026): 166 S_E = 0.41 S_E,t, since S_E,t/S_E = 408. "Plain envelope, margin ≥ 520 over the scalar" compares the graviton tolerance with the *scalar's* edge. It says nothing about the graviton.

Off-diagonal entries: the direct even-J records give Γ_off ≈ 0 in the min-norm gauge (J = 12, 16), as in HIGHSPIN. They are not treated separately here.

## 2. Methods beyond ν/Δ = 0.94

### 2.1 n1: diagnosis (`er_diag*.py`, `er_cut.py`; J = 29, Δ = 1029) [DV]
* The full saddle exists and is smooth through the edge. The smallest Hessian eigenvalue of −Φ falls from 139 (x = 0.96) to 123 (x = 0.97); the largest is 8·10³–10⁴.
* Scalar crystal (Q_J prefactor), moment-matched proposal:
  * x = 0.93: \|μ\| = 0.17, 1.0 % error.
  * x = 0.96: \|μ\| = 5.85, covariance eigenvalues from 0.01 to 3.8, 94 % error, and 7836 masked points.
* Along the line from the saddle to μ, log\|f\| goes −1.7 (saddle) → −3.6 (valley, \|ξ\| ≈ 5) → +12.5 (hump, \|ξ\| ≈ 11). The smooth exponent Φ stays within e^{4.5} of its saddle value. **The growth is in the prefactor:** P_J at complex angle, e^{J·Im γ}, and the graviton's spin-J null-circle projection likewise. Φ, which excludes P_J, does not see it.
* Plain proposals at x = 0.96: 14 % (s = 1.0), 19 % (s = 0.8), 345 % (s = 1.2) error. Radial cuts at \|ξ\| = 3, 4, 5, 6 give 0.93, 1.43, 1.94, 1.97: no plateau. Even at x = 0.90, cuts at 4 and 5 differ by 4 %.

### 2.2 n3 on the scalar (`er_n3test.py`)
n3_hybrid with ρ₁ = 3.5 and the adaptive proposal, M = 9:
* x = 0.93: Es = 0.1783 − 0.0285i ± 1.4 % (n1: 0.1761 − 0.0271i ± 1.0 %).
* x = 0.96: Es = 0.2950 − 0.0610i ± **4.0 %**, with no masked or bad points.

The flow beyond ρ₁ takes the outer points off the hump direction.

### 2.3 n3 against the exact scalar (`er_sval.py`, `cache/sval_J29_part{1,2}.out`) [DV]
T(ν) = [Y_n3(1029, ν)/Y_n3(1100, ν)] / [ex₂(1029, ν)/ex₂(1100, ν)] must be constant in ν. The ν-only factors of the crystal normalisation cancel in the Δ-ratio (n1_scalar_validate's a_ν + b_J form).

| x (ν/Δ) | 0.88 (0.878) | 0.90 | 0.92 | 0.935 | 0.95 (0.948) | 0.96 (0.958) | 0.967 | 0.972 (0.970) |
|---|---|---|---|---|---|---|---|---|
| Re T | 2.565 | 2.557 | 2.568 | 2.711 | 2.618 | 2.585 | 2.650 | 2.435 |
| rel. sd | 0.15 % | 0.8 % | 1.0 % | 3.5 % | 1.3 % | 1.1 % | 3.9 % | 7.5 % |
| T / mean(0.88–0.92) | 1.000 | 0.997 | 1.001 | 1.058 | 1.021 | 1.008 | 1.033 | 0.950 |

\|Im T\| ≤ 0.08. Every deviation is within 1.6σ of its replica error. **n3 reproduces the exact scalar to the light-cone edge, at 1–9 % per node (M = 10).** The error grows toward x₁, so more samples are needed there.

### 2.4 n3 on the graviton (`er_points.py`, `er_n3edge.py`) [DV]
State and reference are on common random numbers, M = 10, R = 8. Each node takes 490–630 s on one core.

| node | EE (k = 2) state | C3H (k = 3) state | references |
|---|---|---|---|
| J = 26, x = 0.93 (n1 works) | n3 −1.218·10⁻⁷ ± 2 % vs n1 −1.223·10⁻⁷ ± 1.5 % (all 4 entries within 1 %) | **n3 ≈ 10⁻⁹ (31–111 %) vs n1 2.36·10⁻⁷ i ± 1.7 %** | EE, C3H agree to 0.3 % |
| J = 26, x = 0.96 | 1.9 % | **99–600 %**, ≈ 3·10⁻⁹ | 0.35 % |
| J = 29, x = 0.9365 | n3/n1 ratio differs by 2 %; Γ_C2 3.58·10⁴ vs 3.12·10⁴ c₂ | **n3 −3.3·10⁻⁹ (78 %) vs n1 4.65·10⁻⁶ i (0.9 %)** | agree to 0.4 % |
| J = 29, x = 0.9605 (n1 failed) | 3.2 % → **Γ_C2 = −7.8·10⁴ c₂** (trend (1−x)⁻² from 3.1·10⁴ at 0.9365 predicts ≈ 8·10⁴) | 256 % | 0.5 % |

* **EE: usable beyond 0.94** [DV, two nodes].
* **C3H: broken at the state** for both parities (p8 odd / o_H q even integrand on the scalar thimble's flow lines), already at x = 0.93. Γ is m² × a small relative difference, so a 1 % error in the state/reference ratio gives ≈ 10⁴ c₂. That is why n3's EE noise must also come down near the edge.
* The C3H failure is the reason n3 has been "poor at large ν". Its cause was not found here [I]. Candidates:
  * a branch or sign continuation of the C3H prefactor along the flow lines (it cancels to ~10⁻³);
  * the graviton-adapted proposal (adapt_g) ignoring C3H's dominant region.

## 3. The J = 29 twist-edge row (`er_row.py`, `cache/row_J29.txt`) [NE / DV]

Data:
* HIGHSPIN's direct n1 row: 34 bulk nodes, GK15 edge on [0.8, 0.9722].
* The four GK nodes with x ≥ 0.95 are failed.
* The GK node 0.8007 sits on the Δ_ref reference zero (HIGHSPIN §2).

Cutoff: the nearest W state's own (J = 30: p_e = 0.4, σ = 0.157).

* cost_2 = +2.8·10⁻⁵ (bulk +2.3·10⁻⁵, edge +0.5·10⁻⁵); cost_3 = −2.1·10⁻⁴. cost_E = **0.266 S_E,t** = 0.22 H.
* Budget, in S_E,t:
  * MC 3σ: 0.52. It is bulk-dominated: 0.36 per kind from the bulk, 0.13 from the edge, at M = 11.
  * quadrature: 1.91. This number is meaningless while nodes are missing.
  * representation: 0.34.
  * reference-zero node: 1.24.
  * missing nodes: 1.72.
* Upper bound **6.0 S_E,t** (= 4.9 H). The 38-cutoff grid gives the same: median 5.9.
* Against the odd-J tolerance ≈ (1 − 0.08)H = 1.13 S_E,t: **fails as a bound; not a violation.**
* With n3-EE nodes the k = 2 edge could be filled. **The k = 3 part (here 8× the k = 2 part) cannot.**

## 4. The exact scalar at the binding states (`er_scalar.py`, `cache/scalar.txt`) [DV]

k = 2 throughout, G = ex/2 − E₂P_J/2 − (P − J₀) in P4W units, and S_E2 = 2V_Fa2(1)e_sc(1).
* x_c is where ν = 0.94Δ (the n1 limit).
* The extrapolation error is that of a shape model fitted on exact values at x ≤ x_c (0.0125 grid plus the GK15 nodes of [0.8, x₁] below x_c) and continued to x₁.

| J (x₁) | own cutoff | grid max | E(p ≤ p_c) | E(p > p_c) | sign-free \|edge part\| | ×3 growth envelope | extrap. error K = 1 / 2 / 3 |
|---|---|---|---|---|---|---|---|
| 26 (0.975) | **0.285** | **0.405** | −1.33 | +1.04 | 1.14 | 4.55 | 0.76 / 0.68 / 0.000 |
| 30 (0.971) | 0.201 | 0.210 | +1.98 | −1.78 | 1.87 | 4.33 | 1.12 / 0.74 / 0.012 |
| 34 (0.968) | 0.108 | 0.298 | −2.00 | +1.89 | 2.09 | 4.55 | 1.08 / 0.76 / 0.25 |
| 38 (0.964) | 0.151 | 0.158 | +1.45 | −1.30 | 1.61 | 4.14 | 0.83 / 0.60 / 0.44 |
| 40 (0.962) | 0.166 | 0.166 | +1.35 | −1.18 | 1.19 | 3.62 | 0.58 / 0.34 / 0.23 |
| 48 (0.955) | 0.129 | 0.143 | +0.52 | −0.40 | 0.40 | 0.75 | 0.23 / 0.15 / 0.10 |

**Reading.**
* The real-axis edge integral is a near-cancellation (×5–20) between the parts below and above x_c. So (EDGE-R) at the gap states is decided by G on [0.94, x₁], with its phase.
* At J = 26 the exact scalar's own-cutoff value is 0.285, and over AG_PROOF's grid it reaches 0.405 = the transferred tolerance 0.406. The scalar-like margin quoted in AG_PROOF (1.35) holds at the tabulated cutoff (1.42), **but the grid choice matters at the 1 % level.**
* No extrapolation protocol tried is controlled. K = 3 is excellent at J = 26, 30 and fails by 2–3× the value at J = 34, 38.

## 5. What it would take to close (EDGE-R) at the gap states
1. **Fix n3's C3H state evaluator at large ν.** Validate it as in §2.4: the n1/n3 overlap at x = 0.90–0.93, plus the Δ-ratio test of §2.3 if a k = 3 exact reference is wanted (sm_exact's k = 3 is interpolated, [M]). Alternatively, a travelling-wave split of the spin-J projection: integrate the null-circle angle as an 8th variable on saddle-point arcs (the analogue of Lemma F's W±), so that each wave has its own re-centred contour. Not implemented.
2. **Lower the MC noise ×2–4.** The J = 29 MC budget (0.52 S_E,t at 3σ) is above the even-J gap tolerances (0.41–0.63). It is bulk-dominated: run M = 12–13 on the bulk nodes, with two independent n3 seeds near the edge.
3. **Avoid the reference zero.** Place nodes ≥ 0.004 away from x₀ (done in the job list).
4. Then cost_E with a budget (MC + K15/G7 + representation) per state, against the tolerances of §0.

## 6. External job list (`ext/jobs_batch3.txt`, generated by `er_jobs.py`)

**Format:** as `ext/jobs.txt` (TAG METHOD J DELTA SEED ν-list); run with `run_ext.py N jobs_batch3.txt`.

**Per state:**
* the row Δ is set by α(α + 2J + 1) = m²;
* χ's support: x ≥ max(p_w, p_e − 3σ)/mfl;
* bulk: the 0.025 grid (n1);
* edge: two GK15 windows on [0.8, x₁] (n1 for ν/Δ < 0.93, n3 above);
* the two highest nodes with ν/Δ ∈ [0.89, 0.93) are also run with n3 (overlap check);
* every n3 node is run twice with independent seeds (6000 + 8k and +4096; README_EXT's spacing rule);
* nodes are kept ≥ 0.004 from the reference zero (the GK window start is shifted if needed).

**Groups:**

| group | states | jobs | note |
|---|---|---|---|
| P1 | (26, 1.026), (30, 1.030), (34, 1.033) | 161 n1 + 56 n3 | thinnest tolerances |
| P2 | (33, 38, 42, 46, 50, 51) at the gap | 331 n1 + 94 n3 | |
| P3 | 20 states with m ≥ 1.13Λ | 1094 n1, no n3 | x₁ ≤ 0.887, so n1 covers the whole support |

**Cost (measured here):** n1 2.4 min per node (J = 26, M = 11, state + reference, both kinds), n3 8–10 min per node (M = 10). In total ≈ 64 CPU-h (n1) + 23 CPU-h (n3) ≈ **90 CPU-h, ≈ 4–5 h with 22 workers**. Memory ≤ 0.3 GB per job here (README_EXT: 0.7 GB).

**Recommended order:**
1. **P3 and all `*_n1` jobs now.** P3 closes (or fails) 20 states without any new method. The n1 parts of P1/P2 are needed in any case.
2. **The `*_n3` jobs give the k = 2 (EE) edge only** until the C3H fix (§5.1). Their C3H output is unusable (the header of the job file says so).

For the MC target of §5.2, set `HC_M=12` in the environment for the n1 pass; it applies to every job of that run.

**Analysis after the run:** `er_row.py` is written for the odd-J own-mass rows (HIGHSPIN layout). The even-J / general-mass version (bulk shape fit + GK15 edge, both helicities, same budget) is the remaining analysis code. It takes a few lines of `hs_even`-style injection; it is not written yet because no data exist.

## 7. Honest residual
1. **(EDGE-R) for the graviton is [M] at all 30 W states with J ≤ 51.** Nothing here computes a graviton E_χ at a W state.
2. The n3 validation is two EE nodes and eight scalar nodes at one or two J. C3H is broken.
3. The scalar comparisons are k = 2 only (no exact k = 3).
4. The J = 29 row is not a W state. Its tolerance was formed with the neighbouring state's d_χ.
5. AG_PROOF's "exact scalar ≤ 0.30 S_E over all cutoffs" is **corrected** to ≤ 0.405 (J = 26, grid maximum). At the tabulated cutoffs, ≤ 0.285.
6. CPU used ≈ 1.9 h, within the 6 h budget. Two processes at most. A restart interrupted the runs once; the finished points were kept and the rest rerun.

## 8. Suggested inventory changes (not applied)
* **AG_PROOF §5.3 / §0 row E / W-χ:** "exact scalar ≤ 0.30 S_E over every cutoff" → "≤ 0.405 S_E2 over the grid (J = 26, the binding state; EDGE_R §4), ≤ 0.285 at the tabulated cutoff; scalar-like +− margin 1.42 there, ≈ 1.0 over the grid".
* **AG_PROOF §7.4 / §9 item 2:** "real-axis graviton values near the edge at J ≤ 51 would decide it" → add: "on x ∈ [0.94, x₁], where n1 fails (a P_J-prefactor hump, EDGE_R §2.1). n3 works for the scalar and for EE there, but n3's C3H is broken at large ν (EDGE_R §2.4). Real-axis extrapolation from x ≤ 0.94 is not controlled (EDGE_R §4)."
* **HIGHSPIN §1 / EXT_FOLD (method notes):** "n1 fails for ν/Δ ≳ 0.94" → add the cause (EDGE_R §2.1) and "n3 C3H fails at the state for ν ≳ 950 even where n1 works".
