# Premise RD for the EE integrand below ν = 10: exact reduction to the scalar, an explicit far-region exponent, and what is left

*22 Nov 2026, RD-close agent. Unrefereed, not committed. No existing document or code file was edited. New code in
`ads4/code/rd/`: `rd_core.py`, `rd_check_core.py`, `rd_divisors.py`, `rd_saddle.py`, `rd_flow.py`, `rd_queue_flow.sh`,
`rd_summary.py`, `rd_corner.py`, `rd_dlib.py`, `rd_dprime.py`, `rd_dprime_int.py`, `rd_census.py`. Caches in `ads4/code/rd/cache/`
(git-ignored by its own `.gitignore`). Read-only imports: `spinfunc/l1_tps` (complex-symmetric Cholesky), `spinfunc/n1_direct` and
`spinfunc/l1_crystal` (only in the validation `rd_check_core`), `repro_engine/integrand` (the EE polarisation contraction).
Nothing in `HLJ_CLOSE.md` or `spinfunc/lj_*` was read or touched. Compute ≈ 1.2 CPU-h, at most 2 processes, each under
`ulimit -v 3000000`, `OMP_NUM_THREADS=1`.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]** derivation sketched; **[DV]** derived and checked numerically; **[NE]** numerical
evidence; **[M]** model; **[A]** assumption; **[I]** open.
Notation as in HRS_CLOSE: `a = 1/Δ`; raw variables `z = (x̃₀, u_x, v_x, u_y, v_y, x⊥, y⊥)` of `l1_crystal`; scaled `X = u/a` etc.;
`B = B_x B_y`, `B_x = (1 − s_x + q_x)/(1 + s_x + q_x)` with `s = u + v`, `q = uv − p²` (raw); `λF = (iν/2) log A + (Δ/2) log B`.
"Raw radius" `r` = `max(|u_x|, |v_x|, |u_y|, |v_y|)` (thimble tables) or `max(s_x, s_y)` (real-section tables); RD's far region
`|raw x| ~ 1` is `r ≳ 0.3–0.5` in either unit.

---

> **Update (PCYC_CLOSE.md, unreviewed): P_cyc, as H_RS needs it, is an identity. The physical integration chain IS D′. H_RS [D given P_cyc] → [D].**
> * **Definition [P] + [DV].** CLS-D (arXiv:2311.04271) define the functional as an absolutely convergent integral over real Lorentzian points in the causal order 4>5>3>0>1>2. In the gauge x₃=−x₁=y, x₄=−x₂=x/(−x²), x₅=e this region is exactly D′, the real section cut at the light cones of x₅ (the zeros of B). Checked with an independent Minkowski code: 0 mismatches in 2·10⁵ configurations.
> * **Integrand [D] + [DV].** On D′ all 25 factor bases are positive, A, B > 0 and Q is real. So the physical integrand = point-independent constant × the project's integrand, and the constant cancels in ρ.
> * **End to end, scalar [DV].** ∫_{D′}, exact in Δ with the cut, against the exact CHMRS Mellin skeleton: |ε| ≤ 2·10⁻⁸ (Δ = 20), ≤ 6·10⁻⁹ (Δ = 30) for ν ≤ 5, and ≤ 1·10⁻⁹ to ν = 10 at Δ = 50. At these Δ, 82 % / 47 % / 3 % of ∫|Φ| lies beyond raw radius 0.3, outside the convergence region of the 1/Δ series. So this checks the chain itself, non-perturbatively.
> * **H_RS [D]** = RD_CLOSE Prop. D′ (engine = ∫_{D′} + O(Δ^N e^{−0.336Δ}); one [D-sketch] sub-step: the convergence radius of the x₀-dependent slow factors) + the identification above. REPRO_ENGINE's direct ∫_{D′} = engine K = 5 to 3·10⁻⁷ [DV].
> * **The crystal (saddle-point) cycle is only a numerical method; its class is no longer a premise.** Crystal = D: [P] in the Regge-scaled system; [I] at finite Δ (gap: no Stokes jump for a ∈ (0, 1/Δ], n± = 1 through the fold). A direct EE test T = ⟨Q⟩_crys/⟨Q⟩_{D′} − 1 finds no discrepancy (|T|/σ ≤ 1.7), but it is weak: (|T|+2σ)/threshold ≥ 4.8 above the fold and 16–2000 below it [NE, weak].
> * **Retraction.** RD_CLOSE's census statement "no critical point at raw radius 0.1–39" is withdrawn. There are 0–276 critical points per run at 0.3 ≤ r < 1.5, many with |B| ≈ 1, at heights 30–94 above the crystal (12 × 4000 starts, not saturated) [NE]. This does not touch D′.

## 0. Summary

| item | result | grade |
|---|---|---|
| **Reduction to the scalar** (§1) | `Φ_EE = a·Q·Φ_scal` **pointwise and exactly**. `Q = m₁₂⁻²m₃₄⁻²V₀₄₂⁻⁴·pref·((v′−u′)² − 1)/4` is a single-valued rational function of the raw point (and of `y₀`). Checked against the reproducer's independent EE integrand: log-difference constant to ≤ 9·10⁻⁹ while `log|Q|` varies by 14–16 units, at Δ = 30, 300, 1000; the constant is `log a` to 10⁻¹¹. | [P] (algebra) + [DV] |
| **Divisor table of `Q`** (§1.2, 60-digit arithmetic) | `Q` has poles only on `FP1 = 0` (`m₁₂ = 0`, order 4), `m₂₀ = 0`, `m₄₀ = 0` (order 4), `m₁₅, m₂₅, m₃₅, m₄₅ = 0` (order 2), `c₁₃ = 0` and `y² = 0` (order 2). **Each of these is a singular divisor of the scalar** with exponent `−11/2`, `−1/4 − iν/2`, `∓Δ/2 (+3/2)`, `+11/4`. `Q` is **regular** on the scalar's only other integer-order divisor `FP2 = 0` (`m₁₄ = 0`, simple pole), on `x² = 0` (`s → ±1`), and on `Nn`, `c₁₀`, `c₃₀`, `y₀ = 0`, `v′ = u′`. So **EE and scalar live on the same complement `U` with the same local system**. | [P] + [DV] |
| **The real cycle `D′`** (§2) | `D′` = real section cut at the zeros of `B` (`m₄₅ = 0`, `m₃₅ = 0`: the light cones of `x₅ = e` from `x₄`, `x₃`). On `D′` every non-integer-power factor keeps its sign, `A > 0`, `FP1 ≥ ½(B_x^num + B_y^num) > 0`, `FP2 ≥ 1`, `|m₁₅| ≥ 1`, `|m₂₅| ≥ 4`. So `D′` is an admissible twisted cycle for scalar and EE. The **untruncated** section `D` is **not**: beyond the `B` zero it crosses `FP1 = 0` (exponent −11/2 scalar, −19/2 EE). That is the concrete content of BELOW_FOLD's "far region of D". | [D] + [DV] (2·10⁵ samples per Δ) |
| **Analytic far-region bound** (§2.3) | On `D′`: `|B|^{Δ/2} ≤ exp(−2Δ[artanh(s_x/2) + artanh(s_y/2)])` [P; checked, max excess −2·10⁻⁴]. On the polydisc of raw radius ≤ 1/3 every slow factor and `Q` is analytic and non-vanishing off the atom divisors, so the `a`-series converges there. Hence **`∫_{D′}Φ_κ = (engine series, all orders) + O(Δ^{N_κ} e^{−cΔ})`** with `c = 2 artanh(r/2)` for the split radius `r`: 0.30 at `r = 0.3`, 0.336 at the guaranteed `r = 1/3`; the measured weight follows the same law up to `r = 1/2` (0.51). This holds for κ = scalar and EE alike. EE differs only through `N_κ`. | [D] |
| **Check on `D′`, integrated far weight** (§2.4) | `W(r, Δ) = ∫_{D′, far}|Φ|J / ∫_{D′}|Φ|J`, importance-sampled at Δ = 30, 60, 100, 200 and fitted to `e^{−cΔ}Δ^p`. **Fitted exponents:** scalar `c = 0.194 / 0.300 / 0.513` and **EE `c = 0.200 / 0.309 / 0.522`**, against the analytic **0.201 / 0.302 / 0.511** (r = 0.2 / 0.3 / 0.5). EE/scalar costs only a power, `p_EE − p_scal ≈ 2–2.5`. | [DV] |
| **The crystal descent cycle itself** (§3) | Downward flow lines of n4's scalar flow function `Φ_w` from the crystal saddles: `A±` at ν = 1, 2, 3, 4, 6, 8; `S` at 12, 16. Integrated in the action parametrisation, so the exact weight is `e^{−s}`, with ≈ 3000 lines per saddle at Δ = 30, 100, 300, 1000. **Far-region exponent `c(r) = min s/Δ` at raw radius ≥ r:** the fit `c_∞ Δ − p log Δ + b` over Δ = 100–1000 gives **0.102, 0.200, 0.319, 0.562, 0.881** at r = 0.1, 0.2, 0.3, 0.5, 0.7. The analytic value is `artanh r` = 0.100, 0.203, 0.310, 0.549, 0.867. EE (`e^{−s}|Q|/Q_eff`): 0.102, 0.206, 0.316, 0.535, 0.879. | [NE] + [D-sketch] (analytic leading order) |
| **Loci the descent cycle must avoid** (§3.3) | Over all ≈ 10⁵ lines: `|m₁₅| ≥ 0.90`, `|m₂₅| ≥ 4.2` (the pole of `B` is never approached: `h → +∞` there), `|FP1| ≥ 0.125` (`m₁₂ = 0` never approached). Lines terminate only on `m₄₅ = 0`, `m₃₅ = 0` (the `B` zeros, i.e. **the boundary of `D′`**), on `s → ±1`, and on a corner `y² = 0, α → ±i∞`. At `s → ±1` `Q` is regular (60 digits), and the EE weight decreases on every line. | [NE] |
| **Far-region EE weight at the binding Δ_g** (§3.2) | Pointwise maximum of `e^{−s}|Q|/Q_eff` on the crystal thimbles (worst ν < 10), r ≥ 0.2 / 0.3 / 0.5: **e^{−37} / e^{−66} / e^{−135}** (Δ = 300), **e^{−171} / e^{−277} / e^{−507}** (Δ = 1000). The failure thresholds are **e^{−8.5}** (2·10⁻⁴) and **e^{−10.8}** (2·10⁻⁵). | [NE] |
| **New finding: an angular corner** (§3.4) | On the scalar thimbles at **ν ≤ 2**, lines ending on the corner `y² = 0, α → ±i∞` (`cos φ → 0` with `x₀ → ∞`) carry an EE weight that **rises** toward the corner: on 13–75 % of those lines at ν = 2, 97 % at ν = 1, none at ν ≥ 3. Its pointwise size beyond n4's cut `s = 20` is up to `e^{−0.2}` (Δ = 100) / `e^{−7.6}` (300) / `e^{−7.3}` (1000) of the near region. **The untruncated scalar thimble is not an admissible cycle for the EE integrand at ν ≲ 2.** The effect is Δ-independent and lies in the Regge-scaled angular sector, so it is **not an RD (far-region) effect**. It plausibly explains HRS_CHECK's loss of the EE signal below ν ≈ 2. The real-angle cycle `D′` is admissible. | [NE] + [D-sketch] |
| **Verdict on RD for EE, ν < 10** (§4) | **The far-region part of RD is closed.** On `D′` it is [D] with `c = 2 artanh(r/2)` ([DV], integrated). On the crystal thimbles it is [NE], with `c(r) → artanh r` matching the analytic leading order. **The EE-specific part is closed.** `Φ_EE = aQΦ_scal` with every pole of `Q` on the scalar's own singular divisor [P+DV]; EE far exponents equal the scalar ones to ≤ 5 % (thimbles) and ≤ 3 % (`D′`), at a polynomial cost `Δ^{2–3}`. **What remains** is the form-independent premise `P_cyc`: the physical cycle lies in the twisted-homology class of `D′` (equivalently, of the regularised crystal thimbles). `P_cyc` is the same premise for scalar and EE, because by Lemma F the homology is literally the same, and it is [DV] for the scalar through its periods. | far part **[D]/[DV]**; RD_EE net **[D given P_cyc]**, `P_cyc` [DV] (scalar periods) |
| **H_RS for Q's binding state** | From "[D given RD], net [D-sketch]" (HRS_CLOSE) to **[D given P_cyc]**. Nothing EE-specific is left unproved. The residual failure mode is sharp: a cycle class outside the Regge region whose scalar period vanishes by cancellation, to < 5·10⁻⁸ (Δ = 100) / 6·10⁻⁹ (300) / 1·10⁻⁹ (1000) relative at every ν ≤ 10 (J = 0; J = 2, 4 checked at Δ = 100 to ≲ 10⁻⁶), while its EE period is ≥ 2·10⁻⁴ / 2·10⁻⁵ and **the same at Δ_g = 300 and 1000**. A pilot census found **no critical point** of the flow function at raw radius 0.1–39 (§4.3). | see §5 |

**Bottom line.**
* The "far region" of RD (raw radius ≳ 0.3, where `m₁₂ = 0`, `m₁₄ = 0` and the pole of `B` live) carries weight `e^{−cΔ}` with
  `c ≈ 0.3` (split at r = 0.3; 0.336 guaranteed at r = 1/3) to `0.5` (r = 1/2, measured). That holds on the real cycle `D′` (derived, and checked as an integrated weight,
  fitted exponent within 3 % of `2 artanh(r/2)`) and on the crystal descent cycle itself (measured along its flow lines at
  Δ = 30–1000, fitted exponent within 1–3 % of `artanh r`).
* The EE integrand is the scalar one times a rational function with no new divisor. It has the same exponent and costs only a power
  of Δ.
* At the binding states the far EE weight is ≤ `e^{−37}` (300) and ≤ `e^{−171}` (1000) already at r ≥ 0.2, against required
  `e^{−8.5}` / `e^{−10.8}`.
* What is not proved is the form-independent cycle identity `P_cyc`, which the validated scalar shares.

---

## 1. Reduction to the scalar: the exact factorisation (`rd_core.Qratio`, `rd_check_core.py`, `rd_divisors.py`)

### 1.1 The identity [P + DV]
* The EE `k = 4` integrand (REPRO_ENGINE §3.1, built from `bg_geom`/`bg_ee`) and the scalar `k = 2`, J = 0 integrand
  (`l1_crystal.build`) share every factor raised to a non-integer or ν-dependent power: `A^{iν/2}`, `B^{Δ/2}`, measure, dual and
  celestial factors, `(−x²)^{−6}`, FP.
* They differ by `m₁₂⁻² m₃₄⁻² V₀₄₂⁻⁴ (v′−u′) R·pref / f`, with `f = (v′−u′)/(u′v′)` and `R = ((v′−u′)²−1)/(4u′v′)`:
  ```
  Φ_EE = a · Q · Φ_scal,      Q = prefĥ · ((v′−u′)² − 1) / (4 m₁₂² m₃₄²),      prefĥ = pref/V₀₄₂⁴ = F·Ê·Ê·G,  Ê = E/v₀₄₂²
  ```
  `F`, `G` are the conserved ⟨TTO⟩ quadratic forms and `Ê = (H₂₄ − 2V₂,₀₄⊗V₄,₀₂)/(X₂₀X₄₀)`. All of these are polynomial in the raw point
  and in `y₀`, divided by products of `m_ij`, so `Q` is single valued on the double cover (`y₀ = ±√y₀²`).
* Check (`cache/rd_check_core.txt`):
  * 400 random points of `D` per Δ, with `log|Q|` spanning 14–16 units;
  * `lG_EE(repro_engine) − [Σ_k Re e_k log|f_k| + log J_D] − log|Q|` is constant to **4·10⁻⁹ / 9·10⁻⁹ / 3·10⁻⁹** (Δ = 30 / 300 / 1000);
  * the constant equals `log a` to ≤ 3·10⁻¹¹;
  * `Q` is real on `D`.
  * The scalar factor list of `rd_core` reproduces the TPS gradient of `n1_direct.phi_w` to ≤ 10⁻¹² at the crystal and at random
    complex points.

### 1.2 Divisor table (`cache/rd_divisors.txt`; Δ = 300; three generic complex points per divisor; `Q` in 60-digit arithmetic)

| divisor | meaning | scalar exponent (Re) | `Q` exponent | EE exponent | status |
|---|---|---|---|---|---|
| `FP1 = 0` | `m₁₂ = m₃₄ = 0` (`m₁₂ = −FP1/x²`) | −5.50 | −4 | −9.5 | scalar branch divisor (half-integer) |
| `FP2 = 0` | `m₁₄ = m₂₃ = 0` | −1.00 | 0 | −1 | integer pole, **same for both** |
| `FP3 = 0` | `1 − x²y²` | +1 | 0 | +1 | zero |
| `m₁₅ = 0`, `m₂₅ = 0` | pole of `B` | −150.0 (= −Δ/2) | −2 | −152 | in `S` (pole of order Δ/2 at even Δ, branch otherwise) |
| `m₃₅ = 0`, `m₄₅ = 0` | zero of `B` | +148.5 (= (Δ−3)/2) | −2 | +146.5 | branch, strongly vanishing |
| `m₂₀ = 0`, `m₄₀ = 0` | `x₂ − x₀`, `x₄ − x₀` null | −0.25 (`−¼ − iν/2`) | −4 | −4.25 | branch |
| `c₁₃ = 0` | `A` divisor | +2.75 | −2 | +0.75 | branch, EE still vanishing |
| `y² = 0` | `u_y v_y − y⊥² = 0` (`cos φ = 0`) | +2.75 | −2 | +0.75 | branch, EE still vanishing |
| `x² = 0` | `s → ±1` | +1.75 | **0** | +1.75 | `Q` regular |
| `Nn`, `c₁₀`, `c₃₀` | `A` divisors | −½, −¼, −¼ | 0 | same | `Q` regular |
| `y₀ = 0`, `v′ = u′` | measure, zero of `f` | −½, +1 | 0 | same | `Q` regular |

* The denominators of `Q` are read off its construction: the `V_{i,jk}` give `m₁₅, m₂₅, m₃₅, m₄₅, X₀₂, X₀₄, X₂₄`; `Ê` gives
  `X₂₀X₄₀`; `v′−u′` gives `m₁₃ m₂₄`; and the explicit `m₁₂m₃₄` appears. So the list is complete.
* `m₁₃ = 4q_y` is the `y² = 0` row. `m₂₄ = −4/x²` vanishes only at `x² = ∞`.
* `Q` changes sign on `D′` (it has zeros). That is harmless: its exponents are integers.

### 1.3 What this does and does not give [D]
* **Same complement, same local system.** Every pole component of `Q` is a component of the scalar's singular divisor `S`. The only
  integer-exponent divisor of the scalar besides `B`'s poles is `FP2 = 0`, and there `Q` is regular. Hence:
  * `U_EE = U_scal`, and `L_EE = L_scal` (B1, B2 of HRS_CLOSE, now exact rather than order by order);
  * every chain whose interior avoids `S` (every chain on which the scalar integrand is holomorphic) is one on which `Φ_EE` is
    holomorphic.
* **EE inherits the scalar's cycles** except possibly at a cycle's **boundary divisors**. There `Q`'s pole lowers the exponent by an
  integer, so the integrability exponent can differ:
  * `D′` (§2.2): EE exponent ≥ scalar exponent − 2 and still positive on every boundary divisor;
  * the crystal thimbles (§3.4): not true at one angular corner for ν ≲ 2.
* **The factorisation holds as an identity, not as an inequality.** So "RD_EE ⇐ RD_scal + polynomial bound" is exactly true for
  any representative cycle on which `|Φ_scal|` itself is `e^{−cΔ}`-small in the far region, provided `|Q|/Q_eff` is polynomially
  bounded there (measured: §2.4, §3.2).
* It **fails as a reduction** only for a far cycle class whose scalar period is small by cancellation rather than by smallness of
  `|Φ_scal|`. That is the residual gap of §4.

## 2. The real cycle `D′` and the analytic far-region bound (`rd_dprime.py`, `rd_dprime_int.py`)

### 2.1 Why `D` must be cut, and where [D + DV]
* On the untruncated real section, `B_x` passes through zero at `aX ≈ 1` (the light cone `m₄₅ = 0`). Beyond it `|B_x| → 1` as raw
  `u_x → ∞`.
* `x·y ≈ −½` is reached at `aX₁·aY₂ ≈ 1`, i.e. `FP1 = 0`, `m₁₂ = 0`. There the scalar has exponent −11/2 and the EE −19/2, which is
  not integrable.
* Along a ray in `X₁` at Δ = 100 the scalar integrand returns to `e^{+10}` (vs `e^{−30}` in the near region) at `aX₁ ≈ 50`
  [NE, scratch check, not cached]. The `FP1` row of §1.2 gives the exponents.
* So `D` is not a cycle. The natural real cycle is
  ```
  D′ = {X_i > 0, |s| < 1, |φ| < π/2, α ∈ (0, π) [both sheets] : 1 − s_x + q_x > 0, 1 − s_y + q_y > 0, component of the origin},
  ```
  the causal region of `x₃`, `x₄` bounded by the light cones of `x₅ = e`.
* The engine (Regge series, exponentially decaying radial weights) and every QMC on `D` with Gamma proposals effectively integrate
  `D′`.

### 2.2 Admissibility of `D′` [D + DV] (`cache/rd_dprime.log` (A), 1.96·10⁵ uniform samples of `D′` per Δ, Δ = 30 and 300)
* No factor with a non-integer exponent changes sign or turns complex; `A > 0` everywhere (so `|A^{iν/2}| = 1`).
* `FP1/(B_x^num + B_y^num) ≥ 0.51`, so `m₁₂ = 0` touches `D̄′` only at the corner `B_x^num = B_y^num = 0`. There
  `|Φ_EE| ≲ (B_x^num B_y^num)^{(Δ−7)/2}(B_x^num + B_y^num)^{−19/2}`, which is integrable for Δ ≥ 27.
* `FP2 ≥ 1`, `FP3 ≥ 0.22`, `|m₁₅| ≥ 1.0`, `|m₂₅| ≥ 4.0`.
* **Boundary of `D′`:**
  * the atom divisors (`X_i = 0`, `w = 0`, `cos φ = 0`, `α = 0, π`), where EE exponents are ≥ the scalar's (HRS B2, §3.2) and `Q` is
    regular or has at most the poles of §1.2 with EE exponent still > −1;
  * `m₃₅ = 0`, `m₄₅ = 0` (EE exponent (Δ − 7)/2).

  Both integrands vanish or are integrable on all of them. So `D′` is a loaded twisted cycle for both.

### 2.3 The bound [P / D]
* On `D′`, `q ≤ s²/4` (because `q ≤ uv ≤ s²/4` for `w ∈ [0,1]`), and `B_x = 1 − 2s/(1 + s + q)` increases with `q`. Hence
  ```
  |B|^{Δ/2} ≤ ((1 − s_x/2)/(1 + s_x/2))^Δ ((1 − s_y/2)/(1 + s_y/2))^Δ = exp{−2Δ[artanh(s_x/2) + artanh(s_y/2)]}     [P]
  ```
  Numerically the maximal excess is −2·10⁻⁴ (Δ = 30) and −1·10⁻³ (300).
* **Convergence region of the Regge expansion.** On the polydisc of raw radius ≤ 1/3 (light-cone and transverse components,
  angles real):
  * `|FP1 − 1|, |FP2 − 1| ≤ 2|x·y| + |x²y²| ≤ 0.49`;
  * `|1 ± s + q| ≥ 1/9`;
  * the `x₀`-dependent factors are FP-type in `x·x₀` with `|x₀|² = −y²`.

  So no divisor of `S` other than the atoms is met. The `a`-expansion of every slow factor, `Q` included, converges uniformly there,
  with ratio ≤ 3r. This is [D] for the `FP` and `B` factors. It is [D-sketch] for the `x₀`-dependent ones (`m₂₀`, `m₄₀`, `Nn`,
  `c₁₀`, `c₃₀`, `c₁₃`), supported by their sign-definiteness on `D′` [DV] and HRS_CLOSE B6.
* **Proposition D′ [D].** For κ ∈ {scalar, EE} and any `K`,
  `∫_{D′}Φ_κ = Σ_{k≤K} a^k I_{κ,k}^{engine} + O(a^{K+1}) + O(Δ^{N_κ} e^{−cΔ})`, with `c = min(r, 2 artanh(r/2)) = 2 artanh(r/2)` for
  split radius `r ≤ 1/3` (`c(1/3) = 0.336`).
  * The near part is the convergent series termwise.
  * The engine's termwise integrals over all of `D` differ from those over the near part by Regge tails `≤ e^{−rΔ}·poly`.
  * The far part of `D′` is bounded by the inequality above times `sup|slow|`. The slow factors are polynomial in the raw
    coordinates away from `S`; near the boundary divisors they are compensated as in §2.2.
  * EE and scalar differ only in `N_κ`, because `Q` adds a polynomial of bounded degree.
* **Measured** (§2.4): the far weight follows `e^{−cΔ}` with exactly this `c` even for split radii up to 1/2.

### 2.4 Numerical check on `D′` [DV] (`rd_dprime_int.py` → `cache/rd_dprime_int.txt`; 2¹⁸ importance samples per entry)
`W(r, Δ)` = integrated `|Φ|J` over `D′ ∩ {max(s_x, s_y) ≥ r}`, divided by the integral over `D′`. `|Φ|` on `D′` does not depend on ν.

| r | analytic `2 artanh(r/2)` | scalar: log W at Δ = 30 / 60 / 100 / 200 | fitted `c` (p) | EE: log W | fitted `c` (p) |
|---|---|---|---|---|---|
| 0.2 | 0.2007 | +0.3 / −2.3 / −7.5 / −23.6 | **0.194** (+4.9) | +0.6 / −0.7 / −5.0 / −20.3 | **0.200** (+6.9) |
| 0.3 | 0.3023 | −0.7 / −6.1 / −15.3 / −41.6 | **0.300** (+5.4) | +0.1 / −3.8 / −12.1 / −37.6 | **0.309** (+7.8) |
| 0.5 | 0.5108 | −4.1 / −15.7 / −33.4 / −80.8 | **0.513** (+5.5) | −2.5 / −12.5 / −29.4 / −76.0 | **0.522** (+8.0) |

* At Δ = 30 the near sampler under-represents the total (W > 1 there). The fit nonetheless returns the analytic exponent within
  3 % for both kinds.
* **The EE/scalar factor costs a power only**: `p_EE − p_scal ≈ 2.0–2.5`.
* Extrapolated EE far weight at r = 0.3 / 0.5: `e^{−65}` / `e^{−125}` (Δ = 300) and `e^{−272}` / `e^{−480}` (1000), from the fitted `(c, p)` and the Δ = 200 point.
* Converting to the signed integral costs at most the near-region cancellation factor `∫|Φ|/|∫Φ|`. That factor is
  `≲ e^{1.6ν} ≤ e^{16}` at ν < 10 (BELOW_FOLD §0 item 3: the QMC error growth `e^{−πν/2}`-type). Irrelevant here.
* Pointwise maxima on `D′` (`rd_dprime.py` (C), scalar, Nelder–Mead; `cache/rd_dprime.log`):
  * `c(0.5) = 0.334 / 0.434 / 0.482` at Δ = 100 / 300 / 1000, approaching 0.511 from below with `O(log Δ/Δ)`.
  * The EE pointwise maximisation was too slow and non-monotone in `r`, a sign of optimiser failures. It was stopped and is not
    used; the integrated EE check above replaces it.

## 3. The crystal descent cycle itself (`rd_saddle.py`, `rd_flow.py`, `rd_summary.py`, `rd_corner.py`)

### 3.1 Method
* **Saddles.** The n4 construction, re-implemented with an exact (dual-number) gradient. The flow function `Φ_w = λF + G_scal + log J_w`
  (n0 = +1, `n1_direct.zmap` coordinates) is continued from ν = 20 down to ν_c = 10 with σ-symmetrisation. Below the fold it is
  kicked to `A±` and continued in steps ≤ 0.25; `S` is used at ν = 12 and 16.
  * `|∇Φ| ≤ 3·10⁻¹³`.
  * The A± values at Δ = 300 reproduce the n4 pattern: raw `X₁ ≈ 4.35 − 2.36i` at ν = 4; `asym` = 1.14 at ν = 4, 0.50 at ν = 2.
  * At ν = 1 the tracked A± coincide (`asym` = 0).
* **Flow lines.** Directions `n ∈ S⁶` in the saddle's normal coordinates (`L^T(−H)L = 1`): 2048 random directions plus 2 rounds of
  local refinement around the best lines at every radius, i.e. ≈ 3000 lines per saddle. Each line `dw/ds = −L conj(Lᵀ∇Φ)/|Lᵀ∇Φ|²` is
  integrated by RK4.
  * A step is accepted only if the **measured** change `Σ_k e_k Log(f_k′/f_k) + Δlog J` satisfies `|ΔΦ + ds| ≤ 2·10⁻³ ds` and
    `|Im ΔΦ| ≤ 2·10⁻³ ds`. So `s` is the exact drop of `Re Φ` and `Im Φ` stays constant (accumulated `|Im ΔΦ|` ≤ 0.23 over a whole line, with `s` up to 2Δ).
  * The scalar weight of a thimble point relative to the saddle is `e^{−s}`. The EE weight is `e^{−s}|Q|/Q_eff`, with `Q_eff` the
    Gaussian-weighted `|Q|` of the near thimble.
  * `Q` is frozen at its last trusted value where `|1 − x⊥²/(u_x v_x)| < 0.03`: the double-precision contraction is unreliable there
    (REPRO E3), and §1.2 shows `Q` regular at `s → ±1`.
* **Grid.** Δ = 30, 100, 300, 1000; ν = 2, 4, 6, 8 (A±), 12, 16 (S); plus ν = 1, 3 at Δ = 300. 26 runs, ≈ 20 s per saddle.

### 3.2 Far-region exponent on the thimbles [NE] (`cache/rd_summary.txt`; worst saddle per (Δ, ν))

| Δ | ν | c(0.1) | c(0.2) | c(0.3) | c(0.5) | c(0.7) | cEE(0.2) | cEE(0.3) | cEE(0.5) |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 2 / 4 / 8 | 0.021 / 0.017 / 0.012 | 0.091 / 0.091 / 0.079 | 0.172 / 0.180 / 0.164 | 0.394 / 0.394 / 0.375 | 0.67 / 0.71 / 0.67 | 0.108 / 0.057 / 0.045 | 0.214 / 0.131 / 0.118 | 0.402 / 0.328 / 0.308 |
| 300 | 2 / 4 / 8 | 0.057 / 0.057 / 0.055 | 0.145 / 0.152 / 0.145 | 0.248 / 0.260 / 0.248 | 0.467 / 0.490 / 0.490 | 0.80 | 0.146 / 0.128 / 0.122 | 0.237 / 0.231 / 0.219 | 0.462 / 0.455 / 0.451 |
| 1000 | 2 / 4 / 8 | 0.082 / 0.086 / 0.082 | 0.179 / 0.188 / 0.188 | 0.291 | 0.523 | 0.85 | 0.171 / 0.175 / 0.174 | 0.279 / 0.277 / 0.276 | 0.508 / 0.514 / 0.507 |
| analytic `artanh r` | | 0.100 | 0.203 | 0.310 | 0.549 | 0.867 | | | |

ν = 6, 12, 16 and Δ = 30 are in the cache file. The ν-dependence is weak; ν = 1 and 3 at Δ = 300 give 0.260 / 0.248 at r = 0.3.

* **Δ-scaling fit** (drop `= c_∞Δ − p log Δ + b`, worst ν < 10, Δ = 100, 300, 1000):

  | r | scalar `c_∞` | EE `c_∞` | `artanh r` |
  |---|---|---|---|
  | 0.1 | 0.1016 | 0.1017 | 0.1003 |
  | 0.2 | 0.2004 | 0.2058 | 0.2027 |
  | 0.3 | 0.3190 | 0.3155 | 0.3095 |
  | 0.5 | 0.5624 | 0.5345 | 0.5493 |
  | 0.7 | 0.8806 | 0.8792 | 0.8673 |

  Plain linear fits give 0.091 … 0.872. **The exponential Δ-scaling is clean and its exponent is the analytic one.**
* **Analytic leading order [D-sketch].**
  * In the Regge-scaled system only the four radial one-variable factors `X^c e^{−X}` let a descent line leave a compact region.
    The `s`, `t = sin φ`, `z = e^{iα}` factors end their descent on their divisors.
  * Asymptotically along such a line `dX/dt = |X|²`, i.e. `X` becomes real positive.
  * With one light-cone component `u` large and the rest `O(a)`, the exact `(Δ/2) log B_x = −Δ artanh(u) + O(1)`, and the other
    factors contribute `O(log Δ)`. So the first passage at raw radius `r` costs `Δ artanh r − O(log Δ)`.
  * Two large components cost more: `2Δ artanh r` at `u = v = r`.
* **Far EE weight vs threshold** (pointwise maximum of `e^{−s}|Q|/Q_eff` beyond r; worst ν < 10):

  | | r ≥ 0.1 | r ≥ 0.2 | r ≥ 0.3 | r ≥ 0.5 | threshold |
  |---|---|---|---|---|---|
  | Δ = 300 | e^{−11.5} | **e^{−36.7}** | e^{−65.7} | e^{−135} | e^{−8.5} |
  | Δ = 1000 | e^{−73} | **e^{−171}** | e^{−277} | e^{−507} | e^{−10.8} |

  * r = 0.1 (`X ≈ 30–100`) is still inside the Regge-convergent region, which the engine integrates. It is listed only to show the
    trend.
  * The thimble volume factor is at most polynomial, `≲ (rΔ)^7 ≤ e^{40}` at Δ = 300. It does not change the conclusion for r ≥ 0.3.
* `max log(|Q|/Q_eff)` beyond r reaches +30 … +69. These values occur only near the `B` zeros (`Q ∝ m₄₅⁻²`) and at the corner of
  §3.4, where `s` is larger still. That is why cEE tracks c to ≤ 3 % at Δ ≥ 300.

### 3.3 Loci avoided and terminal divisors [NE] (all runs; `cache/rd_flow_*.txt`)
* **Never approached:**
  * the pole of `B`: `min|m₁₅| = 0.90`, `min|m₂₅| = 4.2` (raw units; on `D′` the minima are 1.0 and 4.0);
  * `m₁₂ = 0`: `min|FP1| = 0.125`.

  This is the HRS_CLOSE §2.6 claim, now measured: `h → +∞` at the pole of `B`, so descent lines cannot reach it.
* **Terminal divisors** (smallest normalised distance at the last point):
  * `m₄₅ = 0` and `m₃₅ = 0` (the zeros of `B`): 35–60 % of lines. **The radial ends of the crystal cycle coincide with the boundary
    of `D′`.**
  * `s → ±1`: ≈ 25 %. `Q` is regular there; on every such line the EE weight decreases monotonically (`rd_corner`: 0 rising lines).
  * The corner `y² = 0, α → ±i∞`: ≈ 25 %. See §3.4.
  * A handful of `X_i → 0`.

### 3.4 New finding: the angular corner at ν ≲ 2 [NE + D-sketch] (`cache/rd_corner.txt`)
* About a quarter of the lines end where `u_y v_y − y⊥² → 0` while `x̃₀` stays finite. Then `y₀ → ∞`; in `D` coordinates this is
  `cos φ → 0` together with `ψ, α → ±i∞`, i.e. `z = e^{iα} → 0` or `∞`.
* `Q` grows there; this was confirmed in 60-digit arithmetic along a traced line: `log|Q| ≈ 25 → 60` while `|y²/(u_y v_y)|`
  falls from 5·10⁻⁴ to 4·10⁻⁷.
* EE weight along those lines, beyond n4's cut `s = 20`:

  | ν (Δ = 300) | 1 | 2 | 3 | 4–16 |
  |---|---|---|---|---|
  | corner lines whose EE weight still rises | 97 % | 13–75 % | 0 | 0 |
  | max EE log-weight beyond `s = 20` | −3.7 | −7.6 | −8.1 … −9.8 | ≤ −10.6 |

  At ν = 2, Δ = 30 / 100 / 1000 the rising fraction is 26–34 % / 37–75 % / 12–66 % and the maxima are −5.1 / −0.2 / −7.3.
* **Reading [D-sketch].**
  * In the Regge-scaled system the scalar's α-factor `u^{λ_u}v^{λ_v}` decays like `e^{−|Im α|/2}` as `α → ±i∞`.
  * The EE Laurent polynomial carries positive powers of `cos α`, `sin α`, i.e. `e^{d|Im α|}`, and the ν-dependent phases shift the
    balance.
  * Below ν ≈ 2–3 the EE integrand is therefore **not integrable at the end of the scalar thimble** in the α direction. The scalar
    thimble (untruncated) is not an admissible cycle for EE there, while the real-angle cycle (`D′`, the engine's cycle) is.
* **Consequences.**
  1. As a twisted-homology class the crystal cycle is still well defined for EE: regularise the end with a Hankel loop, as allowed
     in HRS_CLOSE Prop. RS, since the exponent is non-integer. But **numerical EE integrals over the scalar thimble become
     cutoff-dependent below ν ≈ 2–3**. This is a plausible mechanism for HRS_CHECK's loss of the common-random-number signal and
     the ±6000 errors at ν = 2.02, and for the heavy tails of the `k = 4` nodes [D-sketch; not re-run with HRS_CHECK's code].
  2. It is a **near-region, Δ-independent** effect in the angular sector, so it is not part of RD. It does mean that any future
     direct crystal check of EE at ν < 3 should keep the angles real (e.g. a radial thimble × real sphere), or use the EE flow
     function.

## 4. RD for EE at ν < 10: the argument assembled

### 4.1 Statement and status
**RD (EE, ν < 10).** `I_phys,EE(Δ, ν) = c_Γ(ν) · [engine EE series] + O(Δ^N e^{−cΔ})`.

Decompose into:

* **(a) Far-region suppression on a representative cycle.**
  * On `D′`: [D] (Prop. D′) with `c = 2 artanh(r/2)` (0.336 at the guaranteed split r = 1/3); [DV] integrated, fitted 0.200 / 0.309 /
    0.522 vs 0.201 / 0.302 / 0.511.
  * On the crystal thimbles: [NE], `c(r) → artanh r` (fit 0.200 / 0.319 / 0.562), with analytic leading order [D-sketch].
* **(b) EE-specific inputs.** Exact factorisation and divisor table [P + DV]; EE exponent = scalar exponent, power cost `Δ^{2–2.5}`
  [DV]; `D′` admissible for EE [D + DV]; descent cycles avoid `m₁₅ = m₂₅ = 0` and `m₁₂ = 0` [NE]; `Q` regular at `s → ±1` [DV,
  60 digits].
* **(c) Near region.** HRS_CLOSE Prop. RS [D]: any closed (regularised) twisted cycle of the Regge-scaled system gives
  `c_Γ(ν) ×` engine, with one constant for scalar and EE.
* **(d) Class identity `P_cyc`.** The physical cycle is homologous, in the twisted homology of `(U, L)`, to `D′` (equivalently to the
  regularised crystal thimbles; both have their radial ends on the `B` zeros, §3.3).
  * By §1.3 `(U, L)` is the same for scalar and EE, so `P_cyc` is one statement about one topological object. It is not
    EE-specific.
  * For the scalar it is [DV] at the level of periods. With Prop. D′, `∫_{D′}Φ_scal` = engine = exact Mellin skeleton to:
    * ≤ 5·10⁻⁸ relative at Δ = 100, J = 0, ν ≤ 3 (≤ 3·10⁻⁷ at ν = 7–10), from `spinfunc/cache/bg_validate.txt`; J = 2 / 4 there:
      ≤ 3·10⁻⁷ / 3·10⁻⁶ on ν ≤ 3;
    * ≤ 6·10⁻⁹ at Δ = 300 and ≤ 1·10⁻⁹ at Δ = 1000 (HRS_CLOSE §2.6).

(a)–(c) give RD_EE for any physical cycle satisfying (d). **RD_EE (ν < 10) is therefore [D given P_cyc], with `P_cyc` [DV] through
the scalar.**

### 4.2 Why `P_cyc` cannot fail "for EE only" except in one sharp way [D]
Suppose `Γ_phys = D′ + F` with a non-zero class `F`. Split `F = F_R + F_far` at raw radius r ≤ 1/3.

* **`F_R` (Regge region).** By Prop. RS, `∫_{F_R}Φ_EE = c_F · engine_EE` and `∫_{F_R}Φ_scal = c_F · engine_scal`, with the same
  `c_F`. Lemma 2 of HRS_CLOSE (one populated class) applies; the α-odd class is not populated.
* **`F_far` with a representative on which `|Φ_scal| ≤ e^{−cΔ}`** (descent cycles, real cycles; §2, §3). Its EE period is
  `≤ e^{−cEE·Δ}`, `cEE ≈ c`.
* The scalar check then forces `|c_F| ≤ ε_scal + e^{−cΔ}`, and the EE error is `≤ ε_scal + e^{−cEE·Δ}`: 6·10⁻⁹ (300) and 1·10⁻⁹
  (1000), against the needed 2·10⁻⁴ / 2·10⁻⁵.
* **The only remaining failure mode.** `F_far` has no representative on which the scalar is small, so its scalar period is small by
  **cancellation**:
  * < 5·10⁻⁸ (Δ = 100), 6·10⁻⁹ (300), 1·10⁻⁹ (1000) at every ν ≤ 10 for J = 0 (and ≲ 10⁻⁶ for J = 2, 4 at Δ = 100);
  * while its EE period is ≥ 2·10⁻⁴ (300) and ≥ 2·10⁻⁵ (1000);
  * and, for the failure scenario of HRS_CLOSE, **the same at Δ_g = 300 and 1000** in c₂.

  By Lefschetz decomposition such a class is a sum of thimbles of non-crystal critical points σ, with contribution
  `∝ e^{h_σ}·amp(σ)`. Δ_g-independence needs `|B(σ)| ≈ 1`. The scalar's silence then needs `amp_scal(σ) = 0` identically in ν while
  `amp_EE(σ) = Q(σ)·amp_scal(σ) + O(1/Δ) ≠ 0`. That is possible only if σ sits on a zero divisor of the scalar slow factor at which
  `Q` has a compensating pole. §1.2 shows `Q` has **no pole on any zero divisor of the scalar slow factor**: `FP3`, `v′ = u′`,
  `x² = 0` are regular; the pole divisors `c₁₃`, `y²` keep EE exponent > 0. So the failure would need a degenerate (non-Morse)
  coincidence [D-sketch].

### 4.3 Pilot census of critical points [NE, pilot] (`rd_census.py 300 4 400` → `cache/rd_census_D300_nu4.txt`)
* Damped Newton from 400 random far starts gave 172 distinct critical points of `Φ_w` (Δ = 300, ν = 4).
* 71 lie in the Regge region (raw radius < 0.1); these are the other critical points of the product of one-variable systems, covered
  by Prop. RS.
* 101 lie at raw radius ≥ 39, beyond the `B` zeros and outside both `D′` and the crystal thimbles. There `(Δ/2) log|B| ∈ [−17, −2]`
  (`|B| → 1` at raw infinity).
* **None lie at raw radius 0.1–39**, i.e. nowhere in the RD far region `|raw x| ~ 1`.
* The raw-infinity points have sheet-independent heights `dH ∈ [−26, +47]` above the crystal, up to `−(ν/2)(arg A + 2πk)`. On the
  principal sheet several sit far above it. Any non-zero coefficient of such a thimble in the physical cycle would dominate the
  scalar by `e^{+40}` and is excluded by the scalar check. Only sheets shifted by `k ≳ 5` windings of `arg A` could hide one.
* This census is a pilot: one (Δ, ν), 400 starts, completeness not established. It supports §4.2 and does not replace it.

## 5. Grades

| item | before (HRS_CLOSE §2.6, §7) | now | evidence |
|---|---|---|---|
| EE = scalar × rational `Q`, exactly | [D-sketch] (code reading) | **[P] + [DV]** (≤ 9·10⁻⁹) | §1.1 |
| `Q`'s poles ⊂ scalar singular divisor; `Q` regular on `FP2 = 0`, `x² = 0` | partly stated | **[P] + [DV]** (60-digit exponents) | §1.2 |
| `D` vs `D′`: `D` crosses `FP1 = 0`; `D′` admissible for scalar and EE | [I] ("far region of D") | **[D] + [DV]** | §2.1–2.2 |
| far-region suppression on `D′`, `c = 2 artanh(r/2)` | [D-sketch] | **[D]**, checked integrated **[DV]** (fit within 3 %, EE and scalar) | §2.3–2.4 |
| far-region suppression on the crystal thimbles, `c(r) → artanh r` | [D-sketch] | **[NE]** (Δ = 30–1000, ν = 1–16) + **[D-sketch]** analytic leading order | §3.2 |
| descent cycles avoid `m₁₅ = m₂₅ = 0` and `m₁₂ = 0` | [D-sketch] | **[NE]** (min distances 0.90, 4.2, 0.125) | §3.3 |
| EE/scalar polynomially bounded where it matters, incl. `s → ±1` | [D-sketch] | **[DV]**: `Q` regular at `s → ±1`; power cost `Δ^{2–2.5}` (integrated); cEE within 3 % of c | §1.2, §2.4, §3.2 |
| angular corner: scalar thimble not EE-admissible at ν ≲ 2 | — | **[NE]** + [D-sketch] mechanism; not RD | §3.4 |
| class identity `P_cyc` (physical cycle ≃ `D′` ≃ regularised crystal) | implicit | **[DV] for the scalar** (periods); form-independent by §1.3, so the same premise for EE; residual failure mode §4.2 | §4 |
| **RD for EE, ν < 10** | [D-sketch] | far part **[D]** (`D′`) / **[NE]** (thimbles); **net [D given P_cyc]**, `P_cyc` [DV] via the scalar | §4.1 |
| **H_RS for Q's binding state** (`++`/OO J = 0, Δ_g = 300, 1000; `++` J = 2) | [D given RD], net [D-sketch] | **[D given P_cyc]**. No EE-specific premise left. Error budget unchanged: truncation `|Δs_d| ≤ 0.015 / 0.0015`; non-perturbative `≤ e^{−37}` / `e^{−171}` relative already at r ≥ 0.2 | §4 |

**Not claimed:**
* [P] for RD;
* any statement for J ≥ 4 graviton states, OO beyond what Γ_OO = Γ_++ gives, or the Result D window (J²/Δ non-uniformity,
  untouched);
* completeness of the critical-point census;
* that HRS_CHECK's ν < 3 EE thimble numbers are wrong: §3.4 says only that they are cutoff-sensitive.

## 6. What remains
1. **`P_cyc` as a theorem.** Identify the CLS-D functional's defining contour and show it is `D′` (the causal region bounded by the
   light cones of `x₅` from `x₃, x₄`) or homologous to it. The `B` zeros as the radial ends of both `D′` and the crystal thimbles
   (§3.3) is strong circumstantial evidence. A complete census of critical points with their intersection numbers, at one
   (Δ, ν) per fold region, would turn §4.2 into [D].
2. **§3.4 follow-up.** Rerun one HRS_CHECK ν = 2.02 EE node with the angles kept real (radial thimble × real sphere) to test the
   cutoff explanation of the ν < 2 CRN loss.
3. **EE pointwise sup on `D′`** (`rd_dprime.py` (C), EE part): the optimiser was too slow and unreliable and was stopped. The
   integrated check (§2.4) supersedes it.

## 7. Reproduce (`ads4/code/rd/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)
```
python3 rd_check_core.py            # 2 s    -> cache/rd_check_core.txt   (gradient vs n1_direct TPS; EE = a Q scalar on D)
python3 rd_divisors.py              # 3 min  -> cache/rd_divisors.txt     (local exponents, 60-digit Q)
python3 rd_saddle.py 300 16,12,10,8,6,4,3,2,1   # seconds; also 30, 100, 1000 (cache/rd_saddles_D*.pkl)
./rd_queue_flow.sh                  # ~15 min, 2 chains: rd_flow.py Delta nu 2048 2.0 for Delta in 30 100 300 1000, nu in 2 4 6 8 12 16
RD_NSTEP=4000 python3 rd_flow.py 300 3 2048 2.0; RD_NSTEP=4000 python3 rd_flow.py 300 1 2048 2.0
python3 rd_summary.py               # -> cache/rd_summary.txt  (c(r), cEE(r), Delta fits, far EE weights)
python3 rd_corner.py                # -> cache/rd_corner.txt   (EE weight toward each terminal divisor)
python3 rd_dprime_int.py            # 3 min  -> cache/rd_dprime_int.txt (integrated far weight on D', fits)
python3 rd_dprime.py                # (A),(B) in 1 min; (C) scalar sup ~4 min (EE part slow; stopped, see section 6.3) -> cache/rd_dprime.log
python3 rd_census.py 300 4 400      # 10 min -> cache/rd_census_D300_nu4.txt (pilot)
```
