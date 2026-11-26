# LEMU_CERT: a rigorous Lemma U certificate pipeline for the sliver at small/middle L (Δ_g = 1000)

*27 Nov 2026. Research agent. Unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*New code: `ads4/code/lemU_cert/` (outputs in `lemU_cert/cache/`, gitignored except `.gitignore`). Read-only imports:
`epswcert/ew_common.py` (realised weight), `fullpos3/fp3_htdata.py` (EPSW residual cells, Gaussian tail, ∫|ψ|e^{κΛp}, Rint),
`lrtight/cache/lrt3_rcert.json` (F_cert), `prl/data/nf02_psi.json`, `epswcert/cache/k2_continuum_D_1000_rep.json`.*
*Compute: ≤ 3 processes on a 4-core box shared with two other agents (load ≈ 9); all timings below are CPU seconds unless marked wall.*

Grades: **[certified]** Arb balls / exact rationals with every truncation carried as a ball radius; **[P]** proved here on paper;
**[DV]** numerical; **[retained]** an input proved elsewhere (author or other agent) and not re-derived here.

## 0. Verdict (status report; the session ended before the first certified pass)

**Nothing is certified yet.** The pipeline is built and its pieces are validated, but the end-to-end single-c runs (L = 2 and
L = 0.5, both b, realised weight) **FAILED: the final balls were uninformative** (rad P̂₀ ≈ 1e439 at L = 2 and ≈ 1e116 at L = 0.5,
with Ŷ_0 = A₀ correct to 0.7 % (6.1e-6 ± 7.6e-8) and the residual radius as budgeted). The blow-up is in the ball contraction, and its cause is
**not yet located** (§2.5 lists the hypotheses that were tested and excluded). No margin is claimed.

What does stand:
1. **Design measurements [DV]** (§2.3): the contraction cancellation κ is e^{238} for the integrated G at L = 2 and e^{70} at L = 0.5.
   Per frequency, κ ≈ 1 for ν ≤ 300 and 1e110–1e140 near ν ≈ 1000. So the EPSW residual can be bounded absolutely only below
   ν ≈ 350, and above that the realised weight has to be integrated exactly.
2. **Rigorous ingredients [P + Arb], implemented and checked:**
   * Lemma J, a complex-frequency jet envelope via Laplace's integral, 2–47 nats above the true jets;
   * GL errors on Bernstein ellipses, with rigorous nodes;
   * exact evaluation of the realised weight f_θ and its disc sups (reusing EPSW's code);
   * the residual-cell bound and the Gaussian tail;
   * the precision rule for ball jets (a measured loss of ~3400 bits at j = 3003).
3. **c-direction [P, validated DV]** (§3):
   * Exact identity R_{m,j} = x^m d^{2m−j}R'_{m,j}(v), v = d/c, where R' is a polynomial of degree ≤ m whose coefficients in −v are ≥ 0.
   * On top of it, a rigorous Taylor model on real cells with every remainder carried as a ball radius.
   * Validation: it reproduces direct evaluation to all 15 printed digits (`lu_cell_check.py`).
   * Rigorous status: the cell certificate is rigorous as coded (no Lemma C). Its only crude inputs are the Q-majorant tail, which
     enters × a Chernoff binomial tail, and the Cauchy remainder of the jet shift.
   * Not yet run end to end: the test cell was stopped behind the point failure.
4. **Launcher** (`run_lemU_scan.py`): resumable, one process per cell, `--summary`. Its parameter table is provisional until one
   cell passes.
5. **Low L (coordinator's update)**: the global-coefficient route beats per-c jets for L ≲ 0.82 (§4a, DV evidence).

---

## 1. Target and reduction

Sliver target: H_ang,b > 0 for J ≥ 1001, b ∈ {0, 4}, a ≥ 998, β_c = 2 atanh√c ∈ [0.486, 5.1].
Following ANGULAR_ROUTE §80 / FULLPOS4 §8, at fixed c every state with a ≥ 998 has h = x d² ∈ (0, h_*(c)] and

  H_M = P(c, h)/D⁴,  P = Σ_m h^m [𝔄_m + q_b 𝔄_{m−1} + 𝔅_{m−1}],  D = 1 + q_b h,

with 𝔄_m = Σ_j R'_{m,j}(v) Y'^A_j (§3.1 for the R' form; at a single c this is exactly FULLPOS4's A_m y^m with y = h/h_*).
Positivity of P on the twist interval implies H_ang,b > 0 once |H_ang − H_M| ≤ e^{−84} (author §97, **retained**, see §6).
The test is the signed Chebyshev test of §80: lower(P̂₀) − upper(Σ_{ℓ≥1}|P̂_ℓ|) > D_max⁴ e^{−84}, all in ball arithmetic.

## 2. Single impact parameter: the certified Lemma U check (`lu_cert1.py`)

### 2.1 Decomposition of the integrated jets (all errors become ball radii of Y_j)
Y_j = d³(cd)^j G^{(j)}(c)/j!. G_A uses the REALISED repaired weight f_θ = 2Re[θ(μ−i)χ(μ−i)R₋(μ−i)] (θ: ν_c = Λ+100, a = 10),
G_A = (2/Λ)∫₀^∞ f_θ F_{3,μ} dμ; G_D = Λ²∫₀¹ψ₃F_{4,Λp}dp (ideal, as instructed). Split at a frequency ν_cut chosen per c:

| piece | method | error control |
|---|---|---|
| (1) ψ₂ on [0, ν_cut] | Gauss–Legendre in p, rigorous nodes (`arb.legendre_p_root`) | ATAP Thm 19.3: (len/2)(64/15)ρ^{−2N}/(ρ²−1)·sup_{E_ρ}\|ψ y_j\|, sup by Lemma J (§2.2) |
| (2) residual r_θ = f_θ − ψ on [0, ν_cut] | not integrated | \|Y^r_j\| ≤ (2/Λ)Σ_cells mass·sup_cell\|y_j\| (EPSW certified panel masses via `fp3_htdata.eps_cells`, Lemma J at real ν), plus the θ-correction 2\|θ(μ−i)−1\|\|g\| |
| (3) f_θ on [ν_cut, NU_E] | GL panels in μ (width ≤ 4), f_θ evaluated exactly (all RBF terms within 60 of μ, the rest as a radius) | same ellipse formula; sup\|f_θ\| on the disc from `ew_common.f_disc_sup` and a ring-box sup of \|θ\| |
| (4) μ > NU_E | Gaussian tail \|f_θ\| ≤ 33·2S e^{1/(2σ_min²)}e^{−(μ−c_max)²/(2σ_max²)} (\|QR₋\| ≤ 16.3, \|θ\| ≤ 1.01 for μ ≥ 1000) | integrated against the Lemma-J envelope |
| G_D | GL in p on [0,1] | as (1) |

ν_cut is the largest EPSW checkpoint at which the cumulative residual bound (2), measured through the column sums ω_j = Σ_m|Ω_{m,j}|,
is still ≤ 5% of F_cert. N, ρ (and panel orders) are chosen so that Σ_j ω_j·err_j ≤ 10⁻¹⁵ (a choice; the radii themselves are
rigorous whatever ω is). Then Y^B, the R recursion, P and the Chebyshev coefficients are computed in ball arithmetic.

### 2.2 Lemma J (jet envelope at complex frequency) [P + Arb] (`lu_bounds.py`)
For real c, 0 < r < 1/c, p ∈ {3,4} and ν = ν_r + iν_i:
  |y_j(ν)| ≤ r^{−j}(d/min_{|t|=r}|1−z|)^p exp((½+|ν_i|)Lw(r) + |ν_r|Θ(r)),
with Lw ≥ max|log|w||, Θ ≥ max|arg w| over z = c + cd·r e^{iα} and w(φ) = (1+z+2√z cos φ)/(1−z).
*Proof.* Cauchy in t; Laplace's first integral P_μ(X) = π⁻¹∫₀^π w^μ dφ (Re X > 0 ⇔ |z| < 1); |w^μ| = |w|^{−½−ν_i}e^{−ν_r arg w}.
w(φ) runs over the segment [e^{−ξ}, e^{ξ}]; if 0 is not on it, arg w is monotone along it (so max|arg w| = |Im ξ| at an end point)
and min|w| is the distance to the segment, attained at the clamped projection (convexity). Arb boxes in α only, adaptive. []
*Check [DV]:* against true jets at L = 2, ν = 0.1/5/300/990, j ≤ 3003: the bound exceeds |y_j| by 2–47 nats (median ≈ 10–25).
(A first version with φ-boxes overestimated arg w near φ = π by up to 330 nats; replaced by the exact segment argument.)

### 2.3 Measured cancellation (why the design is what it is) [DV, `lu_explore.py`]
κ := Σ_{m,j}|R_{m,j}||y_j| / Σ_ℓ|P̂_ℓ| (absolute contraction over the signed one):

| | integrated G | ν = 0.5 | 5 | 50 | 300 | 700 | 990 | 1050 |
|---|---|---|---|---|---|---|---|---|
| L = 2, log₁₀κ | 103.5 | 0.0 | 0.0 | 0.0 | 0.7 | 12.3 | 110 | 139 |
| L = 0.5, log₁₀κ | 30.2 | 0.0 | 0.0 | 0.0 | 0.5 | 12.1 | 37.5 | 41.3 |
| per-frequency sup_y\|P_ν\| (L = 2) | | 0.26 | 0.096 | 1.0e-3 | 1.8e-2 | 3.6e-2 | 0.16 | 3.8e3 |

So below ν ≈ 350 the residual can be bounded absolutely (κ ≈ 1); above it the realised weight must be integrated exactly.
The integrated-G cancellation (e^{238} at L = 2) is what drives every precision and every c-direction choice.

### 2.4 Precision
Ball jets wrap like the majorant of the jet recurrence, ≈ e^{0.76 j} over the true c^j: ≈ 3400 bits lost at j = 3003 at L = 2
(measured at ν = 5…1100: 2630–3400 bits; a-priori estimate `lu_core.wrap_bits` 3459). Jets therefore run at prec + wrap
(L = 2: 4959 bits; L = 0.5: 5937 bits); R recursion and contraction at 1500 (L = 2) / 3000 (L = 0.5) bits.
At low L the wrap grows: 6600 bits at L = 0.2, 12100 bits at L = 0.06 (cost of batch 2, §5).

### 2.5 Results
| run | nodes (A / D / E) | jet prec | CPU (jets) | outcome |
|---|---|---|---|---|
| L = 2, prec 1500 | 362 / 1806 / 9762, ν_cut = 348.85 | 4959 | 5600 s | **FAIL**: P̂₀ = [± 2.4e439] (b = 0), [± 1.8e439] (b = 4) |
| L = 0.5, prec 3000 | 181 / 847 / 6380, ν_cut = 348.85 | 5937 | 4700 s | **FAIL**: P̂₀ = [± 6.7e115] (b = 0), [± 4.7e116] (b = 4) |

Common to both runs:
* A₀ = 6.1e-6 ± 7.6e-8 (L = 2) and 1.29e-5 ± 3.4e-8 (L = 0.5). The radius is the budgeted residual: 2.7 % and 0.4 % of F_cert.
* The budgeted quadrature errors are ≤ 1e-15 in P.
* An earlier L = 2 run at jet precision 1500 failed for a different reason: ball jets lose ~3400 bits. That was fixed (§2.4).

**Excluded causes of the blow-up.** The debug scripts are in the package.
1. *Budgeted radii alone* (`LU_NOJETS=1`): contracted radii ≤ 1e-3 in B₀ and ≤ 4e-8 in A₀ (tiny in P).
2. *Rounding in the contraction*: with DV midpoints plus the budgeted radii, the result has radius 3e-7 at both 1500 and 4959 bits.
3. *Per-node jet precision*: the windowed loss is 2690–3400 bits at both 4959 and 9000 bits.
4. *Seeds*: they lose ≤ 126 bits.
5. *Node-sum rounding (a priori)*: it needs ≤ 631 bits, against 1559 available.

The L = 2 log reports "max rel radius of jets A = 10⁰", so **some summed jet Y_j has lost all precision**, and the mechanism is not
one of 1–5. The next step is `lu_debug_radius.py 2.0 1500 1e5`, which prints the first m where |R|rad(Y) or rad(R)|Y| explodes and
the j responsible. It needs ~1.3 CPU-h. The run here was stopped at 17 CPU-min to report.
`lu_cert1.run` now caches the jets (`cache/point_*.pkl`) and contracts at the jet precision, so the next attempt can iterate on the
contraction without recomputing jets.

## 3. The c-direction: an exact Taylor model on real cells (`lu_cell.py`)

### 3.1 Exact structure [P; checked]
* **R' form.** R_{m,j}(c, x) = x^m d^{2m−j} R'_{m,j}(v), v = d/c, R'_0 = e₀,
  R'_{m,j+1} = j(j+s)Bv[j] + j(2j+s−1)Bv[j−1] + j(j−1)Bv[j−2], Bv[j] = R'_{m−1,j+1} − vR'_{m−1,j}.
  (Substitute R = x^m d^{2m−j}R' into the FULLPOS4 recursion: every term carries d^{2m−1−j}; checked to 4e-178 at m = 6.)
  So R'_{m,j} is a polynomial of degree ≤ m in v, and R'_{m,j}(v) = Q_{m,j}(−v) with Q having **nonnegative** coefficients
  (induction: Bv = Q[j+1](−v) + (−v)Q[j](−v)). With h = xd² and Y'_j = d^{−j}Y_j = d³c^jG^{(j)}(c)/j!, P = Σ_m h^m Σ_j R'_{m,j}Y'_j(+dressing).
* **Twist over a cell.** h_*(c) = x_*(c)d² is decreasing in c, so on a cell every state has h ∈ (0, h̄], h̄ = h_*(c_lo).
* **Taylor shift of the jets (exact).** Y'_j(c) = (d/d₀)³(c/(c₀d₀))^j Σ_{i≥j} C(i,j) Y_i(c₀) τ^{i−j}, τ = (c−c₀)/(c₀d₀).

### 3.2 The rigorous Taylor model in w = v − v₀, |w| ≤ δ (real)
* Y'_j: exact degree-D Taylor polynomial from the anchor jets (needs Y_i, i ≤ 3M+3+D) + remainder
  crude_j(ρ)(δ/ρ)^{D+1}/(1−δ/ρ), crude_j(ρ) ≥ sup_{|w|≤ρ, complex}|Y'_j| from |Y_i(c₀)| (balls, i ≤ K_j) and the integrated Lemma-J envelope
  beyond (with ∫|ψ|e^{κΛp} and Rint from fp3_htdata).
* R'_{m,j}(v₀ + w): exact Taylor coefficients to order I₁ (the truncated recursion is exact); tail beyond I₁ bounded by
  Q_{m,j}(v₀+δ)·P(Bin(m, δ/(v₀+δ)) > I₁) [P: |[w^i]Q(−v₀−w)| ≤ [w^i]Q(v₀+w) and Σ_{i>I}[w^i]Q(v₀+w)δ^i = Σ_l a_l(v₀+δ)^l P(Bin(l,p) > I)],
  Chernoff for the binomial tail; orders I₀ < i ≤ I₁ enclosed by Horner on the ball [−δ, δ]; all as constant-term radii.
* Dressing (c(w), d(w), 1/(2cd) as exact truncated series + Cauchy remainders), products as ball polynomials, truncations of
  products enclosed by Horner on [−δ, δ].
* Test: lower(P̂_{0,0}) − Σ_{ℓ≥1}|P̂_{0,ℓ}| − Σ_{k≥1}δ^k S_k − D_max⁴e^{−84} > 0, S_k the Chebyshev (k ≤ 40) or monomial abs-sums.
* **Validation [DV]** (`lu_cell_check.py`, single frequency ν = 5 as G, M = 20): the cell polynomial evaluated at w = 0, ±δ equals
  the direct contraction at c(v₀ + w) to all 15 printed digits (−0.0961071532676617, −0.0961811054695036, −0.0960307606398316).

### 3.3 Why no cheaper rigorous c-direction was found [DV numbers]
Any rigorous continuum argument needs an upper bound of |P| on a complex neighbourhood of the cell (or all Taylor coefficients).
Every computable bound is termwise in j and/or uses the R-majorant, and loses
* the contraction cancellation κ (e^{238} at L = 2, e^{70} at L = 0.5), and
* the sign cancellation of the R recursion: Q_{m,j}(v₀) against |R'_{m,j}(v₀)| differs by up to e^{204} already at m = 300 (L = 2;
  median e^{127}), growing roughly linearly in m.
A plain u-Taylor with remainder from such a crude bound S needs order K with K·ln(ρ/δ) ≳ ln(S/F) ≈ 10³ — infeasible; the author's
proven strip bounds (e^{10500}F, e^{1650}) are far worse. The design of §3.2 avoids both losses where it can: the R-majorant enters
only multiplied by a binomial tail beyond order I₁, and κ only multiplied by (δ/ρ)^{D+1} or by the exact high-order Taylor
coefficients of R' (no majorant). The price is I₁ ≈ 200–300, I₀ ≈ 80–150, D ≈ 90–140 and polynomial (not scalar) contractions.
Lemma C (|H_M| ≤ K·F on a complex strip, K ≈ 1.3 measured) would reduce this to a ~10th-order Taylor step; it is not used.

### 3.4 Test cell
Not completed. The anchor-jet stage of a b = 0 test cell (u ∈ [3.3144, 3.3154], I₁ = 200, I₀ = 80, D = 90, ρ = 20δ) was stopped:
the same contraction problem as §2.5 would hit it, and CPU was needed for the point runs.

Measured component costs (wall, box at load ≈ 9; CPU ≈ 0.4 × wall), per channel and per cell:
* R' recursion: ~280–440 µs per (m, j) entry;
* contraction: ~0.5–1.4 ms per entry, at (I₀, D) = (80, 90) … (150, 135);
* total ⇒ 1400–2700 s per channel.

## 4. Laptop launcher (`run_lemU_scan.py`)
Command, to be run only after one cell passes (§2.5):
```
cd ads4/code/lemU_cert
python3 run_lemU_scan.py --plan                      # writes scan_plan.json
python3 run_lemU_scan.py 20 --batch 1 --mem 1500     # beta in [1.47, 5.1]
python3 run_lemU_scan.py --summary
```
How it runs:
* Each cell is one single-threaded process with its own JSONL/log in `cache/scan/`.
* Cells that already have both channel records are skipped.
* `--split-failed` halves failed cells.

Provisional sizing:
* the defaults are width 0.004 in u and (I₁, I₀, D) = (260, 60, 60);
* §5 says D and I₀ must be larger (D ≈ 150–200, I₀ ≈ 100–160);
* edit `params()` once the test cell has passed.

## 5. Cost model and core-hour estimate
Estimates, not measurements of a passing run. Per anchor at L = 2:
* rigorous jets: ~1 CPU-h at ~5000–5500 bits, 12–13k nodes;
* Y' polynomials: minutes;
* per channel: R' recursion + contraction at 1.5e6 entries, 0.3–1 CPU-h each at L = 2.

Sizing for (I₁, I₀, D) (§3.3):
* I₁ ln(I₁/(mp)) − I₁ ≳ 960, with mp ≈ 1000·coth(u/2)·δ_u;
* I₀ from (mp)^{I₀}/I₀! ≲ e^{−260};
* D·ln(ρ/δ) ≳ ln κ + Λρ_u + 20.

With these, the cost per unit of u is roughly constant in δ_u ∈ [0.002, 0.02] at ≈ 150–300 core-h (L = 2, κ = e^{238}).

Totals:
* **batch 1 (β ∈ [1.47, 5.1]): ~600–1100 core-h**, i.e. ~1.5–2 days on 22 cores. It could drop if κ falls at larger L, which is not
  measured.
* **Batch 2 (β ∈ [0.486, 1.47])**: jets need 6600–12100 extra bits, which is prohibitive per c. Use the global-coefficient route
  (§4a) instead.
* Per-process memory: ~400–700 MB at I₁ ≈ 260–360. That exceeds the 500 MB target at the upper end; the R' row is now streamed
  (one row in memory).

## 4a. Low L: global coefficients instead of per-c jets (coordinator / author sliver_next) [DV]

`lu_hybrid_dv.py 0.5 3249 3000 2000`:
* A_k = [z^k]G_A (ideal) come from the author's triangular moment recurrence: 38 s in python, midpoint.
* The jets at c follow by an exact Taylor shift, Y_j(c) = d³(cd)^jΣ_k C(k,j)A_kc^{k−j}: 27 s.
* The DV per-c quadrature jets took 535 s; a certified per-c run takes ~1.3 CPU-h.
* The two agree to 1e-10 at j = 0, ≥ 7 digits at j = 100, and 11 digits at j = 500 and 1000. The j = 1–10 differences of
  1e-8–1e-4 are attributed to the coarse DV quadrature (Np = 2000).
* **They disagree completely at j = 2000 and 3000. Reason: K = 3249 global coefficients are too few for the high jets.** The shift
  tail Σ_{k>K}C(k,j)A_kc^{k−j} converges only for k ≳ j/(1 − c/0.99). With j = 3003, c ≤ 0.56, one needs K ≈ 7000–8000
  (≈ 2.3·3003).

Recommendation for L ≲ 0.82:
* Use global A_k, D_k (and realised A_k from the author's sliver_L23 tools; not duplicated here) with K ≈ 8000.
* Take the jets at each c, and the cell polynomials Y'_j(w), from the exact shift. Y'_j(w) is then an exact polynomial in c, so the
  D-remainder (the κ-driven crude piece of §3.2) disappears.
* Keep the Lemma U contraction and the cell model of §3.
* This removes both the per-c quadrature and the 6600–12100-bit jet precision. The remaining crude input is the R'-tail.
* The author's bivariate prefix numerator is an alternative. Its h-degree is K ≈ 1200–3250, against M + 1 = 1001 for Lemma U.

## 6. What is rigorous, what is retained, what is DV
| item | status |
|---|---|
| Lemma J, GL ellipse errors, Taylor-model remainders, R' identity and Q ≥ 0, binomial tail | [P], implemented in Arb |
| Realised weight f_θ: formula, repaired coefficients | read from EPSW (`ew_common`, repaired JSON) [certified there] |
| EPSW residual masses ∫\|r\| on cells (log-parsed, +2e-10 per cell), Gaussian-tail parameters, ∫\|ψ\|e^{κΛp}, Rint | [retained] (EPSW_CERT, FULLPOS3 `fp3_htdata`) |
| Bound \|QR₋(μ−i)\| ≤ 16.3, \|θ\| ≤ 1.01 for μ ≥ 1000 | [P] (stated in code comment; elementary) |
| Action identity H/w³ = ℒ(G_A + τ𝒟_bG_A + τσ_bG_D); the R recursion and the x_* parametrisation (FULLPOS4/author) | [retained] |
| Truncation \|H_ang − H_M\| < e^{−84} (M = 1000) | **[retained, author §97, directed floating point, not Arb]**; not reproduced here |
| End-to-end single-c certificate | **not achieved** (§2.5) |
| c-direction cells | implemented and validated (DV); not yet run end to end |
| Ideal G_D (not realised) | as instructed |

## 7. Files and reproduction (`ads4/code/lemU_cert/`; `ulimit -v 3000000; export OMP_NUM_THREADS=1`)
| file | content |
|---|---|
| `lu_core.py` | jets, integrated DV jets, dressing, R recursion/contraction, y-polynomial, Chebyshev, wrap estimate |
| `lu_bounds.py` | Lemma J envelope (Arb boxes) |
| `lu_weight.py` | realised repaired weight f_θ (exact evaluation, disc sups) via ew_common |
| `lu_cert1.py` | single-c certificate (`build` = rigorous integrated jets, `run` = both channels) |
| `lu_cell.py` | cell certificate (exact Taylor model in w = d/c − v₀); `LU_CHANNELS=0` or `4` selects a channel; anchor-jet cache |
| `lu_cell_check.py` | validation of the cell model against direct evaluation |
| `lu_explore.py` | κ measurements (§2.3) |
| `run_lemU_scan.py` | resumable parallel launcher (`--plan`, `--summary`, `--batch`, `--split-failed`, `--mem`) |

```
python3 lu_cert1.py 2.0 1500            # single c, L = 2 (both b)
python3 lu_cert1.py 0.5 3000            # single c, L = 0.5
python3 lu_cell_check.py 3.3149 0.001 5.0 20
LU_CHANNELS=0 python3 lu_cell.py 3.3144089 3.3154089 1500 200 80 90 cache/cell_test.jsonl 20 900
python3 run_lemU_scan.py --plan; python3 run_lemU_scan.py 20; python3 run_lemU_scan.py --summary
```
