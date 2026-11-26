# GAP_LEMMAS: lemmas to close the sliver gap 0.33 ≲ L < 1.5 (shared file: designer ↔ prover)

Protocol (coordinator-defined):
* The DESIGNER agent appends lemma entries below. Each entry has: id (G1, G2, …), precise statement (all quantities defined, all hypotheses),
  intended use (which route, which c/L range), feasibility evidence (numerical tests with script paths and results), expected cost of a
  rigorous certificate, and a status line: `STATUS: PROPOSED` → `STATUS: READY FOR PROOF` (designer is confident it is true and useful).
* The PROVER agent takes entries marked READY FOR PROOF (lowest id first), writes a proof (as a subsection under the entry), implements the
  rigorous check (Arb / exact integer / directed rounding), runs it, and sets `STATUS: PROVED (+ check PASS)` or `STATUS: FALSE (counterexample …)`
  or `STATUS: GAP (…)`. The prover never edits the designer's statement; if it needs a change it appends a `PROVER NOTE` and sets `STATUS: NEEDS REVISION`.
* Both agents append only (no rewriting of others' text). Code goes to code/gap_design/ (designer) and code/gap_prove/ (prover).

Context: ANGULAR_ROUTE.md §§104–123. Certified: L ≤ 0.305 (direct method; author claims also [0.312, 0.329]), L ≥ 2.3; laptop launcher
is running L ∈ [1.5, 2.3]. Open: ≈ 0.305 < L < 1.5 (c = 1−e^{−L} ∈ (0.263, 0.777)).

## Seed lemmas (from the coordinator; prover may start on these immediately)

### G1 (physical twist cap for the direct method) — STATUS: READY FOR PROOF
For every physical sliver state (J ≥ 1001, b ∈ {0,4}, s = 2b+1, real a ≥ 998, n = J − b, N = n + a, D = N(N+s), h = 1/D, c = n(n+s)h):
h ≤ min( c/((1001−b)(1002+b)), (1−√c)²/998² ). Hence on a cell c ∈ [c_lo, c_hi] the rectangle h ∈ [0, min(c_hi/((1001−b)(1002+b)), (1−√c_lo)²/998²)]
contains all physical states. Use: replace the h-range in the author's direct_sign_triangle/tight drivers (code/sliver_direct_cont/) by this cap
(as NEW driver files; do not modify the author's files), so cells above c ≈ 0.25 contain no non-physical points (§122). Needed for c up to 0.56
(coefficient table K = 3299 in code/sliver_direct_ext/tables/). Also needed: a smaller h-radius option when 1.5·h·K² ≥ 1 (above c ≈ 0.51).

### G2 (seed-error budget for c ≤ 0.9) — STATUS: PROPOSED
The coherent seed / truncation action allowance used by the direct checkers (1e-18) is proved only for c ≤ 0.56. For c ∈ [0.56, 0.9] the
existing bounds give pole correction < 1.8e-15 and log-seed < 1e-20 (pole_bound_full / raw_seed_budget_full); establish the allowance to use
there (and that it is ≪ the expected lower bounds ~1e-6–1e-5).

### G3 (conformal route) — STATUS: PROPOSED
Author's conformal prefix theorem (code/sliver_direct_cont/EXTENSION.md): with w = (1−√(1−z))/(1+√(1−z)), 5334 transformed coefficient pairs
suffice for |H − H^{(w)}_{5333}| < 2^{−110} on c ≤ 0.9. Missing: (a) the transformed seed-error budget; (b) a continuum-certified, stable way to
evaluate the transformed prefix (five-term recurrence C0..C4 for V_j = T_h[w^j]) on (c, h) cells; (c) generation of the transformed table.

---
### G1 — PROVER: proof and check (prover agent, 26 Nov)

**Proof.** Fix b ∈ {0,4}, s = 2b+1 ≥ 1, an integer J ≥ 1001, n = J − b ≥ 1001 − b, and a real a ≥ 998. Put N = n + a,
D = N(N+s), h = 1/D, c = n(n+s)h. All of n, N, h, c are positive and c < 1.

(i) Spin cap. h = c/(n(n+s)), and n ↦ n(n+s) is increasing on n > 0, so n(n+s) ≥ (1001−b)(1001−b+s) = (1001−b)(1002+b).
Hence h ≤ c/((1001−b)(1002+b)).

(ii) Twist cap. Let f(x) = √(x(x+s)) on x > 0. Then f′(x) = (2x+s)/(2√(x(x+s))), and (2x+s)² − 4x(x+s) = s² ≥ 0, so f′(x) ≥ 1.
By the mean value theorem, f(N) − f(n) ≥ N − n = a ≥ 998. Divide by f(N) = √D = h^{−1/2}: 1 − √c = 1 − f(n)/f(N) ≥ 998·√h. Both sides are
nonnegative, so (1 − √c)² ≥ 998²·h, i.e. h ≤ (1−√c)²/998². This uses only that a is real and ≥ 998. No integrality of a is needed.

(iii) Cell form. On a cell c ∈ [c_lo, c_hi], c ↦ (1−√c)² is decreasing on [0,1]. So every physical state with c in the cell has
0 < h ≤ min(c_hi/((1001−b)(1002+b)), (1−√c_lo)²/998²). The rectangle [c_lo, c_hi] × [0, hcap], with hcap ≥ that minimum, therefore
contains every physical state of the cell. ∎

**Rounding used in the driver** (code/gap_prove/gap_triangle.py, `g1_cap`):
* s_lo = isqrt(⌊c_lo·10^40⌋)/10^20, with s_lo² ≤ c_lo asserted. Then s_lo ≤ √c_lo < 1, so (1−s_lo)² ≥ (1−√c_lo)².
* Both caps are rounded upward to denominator 10^15.
* hcap = min of the two rounded caps. Each rounded cap is a rigorous upper bound, so the minimum is one too.

**Downstream facts the direct checker needs**, all still valid under the smaller cap:
* hcap ≤ 998^{−2} (prefix tail lemma, README §3).
* 1 − hcap·(K−1)(K−1+s) > 0 (asserted).
* normalisation (1−c)³/(1+qh)⁴ ≥ (1−c_hi)³/(1+q·h_hi)⁴ on each h-piece.

**Check** (code/gap_prove/check_g1.py): PASS. It verifies:
* the identity (2x+s)² − 4x(x+s) = s² symbolically;
* both inequalities on 180 exact rational states (n from 1001−b to 10^8, a from 998 to 10^9, including non-integer a), using the squared
  integer-free form 4·998²·n(n+s) ≤ (N(N+s) − n(n+s) − 998²)², with the right-hand base positive;
* that the driver's rounded cap2 dominates (1−√c_lo)²/998² exactly at 20 cell endpoints.

STATUS: PROVED (+ check PASS). Execution (sign cells using the cap) is recorded below under "G1 — cells".

---

## Designer entries (designer agent, 26–27 Nov)

Designer code: code/gap_design/ (non-rigorous numerics use python-flint arb at 3400–4800 bits on table MIDPOINTS; they are
feasibility evidence only). Tables used: K = 3299 (code/sliver_direct_ext/tables, decompressed to scratch) and the author's
corrected 17,201-pair vector coeff_long_17200.txt (full-band package, scratch; not copied into the repo).

### G4 (designer): Route A on the whole band — physical-cap cells, overshoot rule, h-pieces, precision guard — STATUS: READY FOR PROOF

**Statement (computational lemma; conditional on the retained inputs of §§115–120 and SPECTRAL_DEGREE.md).**
Fix b ∈ {0,4}, s = 2b+1, q = b(b+1). Cover c ∈ [0.268, 0.777] by cells [c_lo, c_hi] whose *twist overshoot*
ω := (1−√c_lo)²/(1−√c_hi)² satisfies ω ≤ 1.08 (b = 0) and ω ≤ 1.06 (b = 4) (see evidence: ω ≈ 1.10 still passes for b = 0 at c ≈ 0.4
but b = 4 fails at ω ≈ 1.135). Let H = G1 cap = min(⌈c_hi/((1001−b)(1002+b))⌉, ⌈(1−√c_lo)²/998²⌉), let K be the least integer with
tail(c_hi, K) < 10^{−7}, where tail is the author's exact ordinary-tail bound with envelope 2^{2337} (tight_spectral_tail.py of the
full-band package; any K with tail < 2^{−110}, e.g. need_k, is also fine but ≈ 10 % larger). Split [0, H] into m equal h-pieces, m minimal
such that the complex polydisk with h-radius ratio 2 is pole-free: 1 − (h_hi + dh)(K−1)(K−1+s) > 0 on the top piece (m = 1 up to
c ≈ 0.5, m = 2 to c ≈ 0.65, m = 2–3 above). Run the triangle certificate (author's direct_sign_triangle.cpp, or the unscaled variant
below) on each (cell, piece) with coefficient table coeff_long_17200.txt and working precision D ≥ 100 + log10 G + 30, where
G = max_k Π_{j<k} g_j is the a-priori growth of the C++ l1 rounding recurrence, g_j = (|a|+dc+2r)/(d−r) at step j.
**Claim:** every (cell, piece, b) returns lower = (centre − L1 − roundoff − analytic tail)·(1−c_hi)³/(1+qH)⁴ − allowance > 0, with
allowance = 10^{−7} + 10^{−18} for c_hi ≤ 0.56 and 10^{−7} + 10^{−14} for c_hi ≤ 0.9 (SPECTRAL_DEGREE.md §1; this is how G2 is discharged,
provided that budget is accepted as a retained input). Predicted lower bounds ≈ 0.95·(normalised action) ≈ 1e−5 (see table).
Together with G1 (every physical state of the cell lies in [c_lo,c_hi]×[0,H]) this certifies H_ang,b > 0 on c ∈ [0.268, 0.777].

**Why the overshoot rule (new numerics).** On the physical range y := h/hcap_tw(c) ∈ [0,1] the prefix is affine in y to relative
1e−8 (code/gap_design/scan_h2.py), but just beyond the twist cap it explodes (onset.py; relative deviation from the affine interpolant):

| c | b | y = 1.04 | 1.06 | 1.08 | 1.10 | 1.12 |
|---|---|---|---|---|---|---|
| 0.30 | 0 | −1.6e−6 | 3.8e−5 | 3.1e−4 | −3.5e−2 | 1.2 |
| 0.50 | 0 | −1.7e−6 | −5.3e−5 | 1.4e−3 | 6.6e−3 | −1.5 |
| 0.77 | 0 | −1.1e−7 | −9.1e−5 | 7.4e−4 | 5.2e−2 | −2.4 |
| 0.30 | 4 | −2.3e−7 | 5.4e−4 | −9.4e−3 | −9.3e−2 | 12 |
| 0.50 | 4 | −1.1e−5 | −1.1e−4 | 7.2e−3 | −9.6e−2 | −2.8 |
| 0.77 | 4 | −1.4e−6 | −1.5e−4 | 2.0e−3 | 6.8e−2 | −4.2 |

The onset (y ≈ 1.1) is the same at all c, so the overshoot ω, not the width, is the right cell-size parameter; b = 4 is ~5× more
sensitive. (The author's [0.268, 0.35] FAIL and the c = 0.35, h = 2.33e−7 witness are at y ≈ 1.4.)

**Feasibility evidence.**
* Prediction tool: code/gap_design/box_cheb_h.py samples the prefix at 16×12 tensor Chebyshev nodes of the rectangle and returns the
  centre and L1 of the non-constant Chebyshev coefficients — exactly the quantities the triangle checker certifies. Agreement with real
  runs is to 4–5 digits (below).
* Real certificate runs (author's checker, my driver code/gap_design/sign_hsplit.py = sign_next.py + h-pieces + precision guard):

  | cell | b | piece | K | D | centre | L1 | lower (normalised) | time |
  |---|---|---|---|---|---|---|---|---|
  | [0.3, 0.3097] (ω=1.04) | 0 | [0,1]H | 1403 | 100 | 4.2235e−5 | 1.5339e−6 | **1.3288e−5 PASS** | 48 s |
  | [0.6, 0.6068] (ω=1.04) | 0 | [0,½]H | 3316 | 556 | 1.6857e−4 | 4.3994e−6 | **9.880e−6 PASS** | 220 s |
  | [0.6, 0.6068] | 0 | [½,1]H | 3316 | 100 | — | — | FAIL: roundoff 1e+256 (D too small, see below) | 118 s |

  sampled predictions: [0.3,0.3097] b=0: 4.2235e−5 / 1.5337e−6; [0.6,0.6068] b=0 lower piece: 1.6857e−4 / 4.3994e−6.
* Predicted L1/centre at ω = 1.04 (box_cheb_h.py, logs/arect.log): 0.032 (c = 0.5), 0.026 (0.6, both pieces), 0.029 (0.7),
  0.028 (0.77), both b. I.e. every cell is expected to pass with ~97 % of the centre left.
* The prover's chain (gap_prove/run_chain.py, cells of width 0.02–0.022, ω ≈ 1.08–1.10) is already passing c ∈ [0.28, 0.418] (b = 0) and
  [0.29, 0.396] (b = 4) at 50–185 s per cell, consistent with this lemma.

**Precision guard (important above c = 0.5).** direct_sign_triangle.cpp stores A_k, G_k × 2^{−k}, which doubles the Horner ratio;
its l1 rounding error is multiplied by g_j = (|a|+dc+2r)/(d−r) ≈ 2(|c−h₀p|+dc+dh·p)/(1−(h₀+dh)p) per step. For c > 1/2 the product
reaches 10^{427} (c = 0.6, lower h-piece, K = 3316): with D = 100 the run returns garbage (roundoff 1e256, centre 1e68) — not a sign
failure. Either set D from the a-priori growth (done in sign_hsplit.py and in the prover's run_chain.py) or use the **unscaled variant**
code/gap_design/direct_sign_triangle_u.cpp + sign_hsplit_u.py (coefficients stored without 2^{−k}; ratio not doubled; three lines
changed; reproduces the [0.3, 0.3097] cell bit-for-bit in centre/L1/lower, but is 2× slower there (97 s vs 48 s) because the integers are
larger; timing comparison at c = 0.6 pending).

**Expected certificate cost** (cost_model.py: exact replica of the driver's majorant/Q/terms selection; time ≈ 1.3e−7 s × K × terms,
calibrated on 5 real runs): per channel, ω = 1.04: 1.7 (c = 0.3), 3.0 (0.4), 6.8 (0.5), 16 (0.6), 45 (0.7), 107 (0.77) CPU-h per
unit c, i.e. ≈ 11 CPU-h per channel for [0.268, 0.777] at ω = 1.04; ω ≈ 1.08 roughly halves this (the prover's 0.022-wide cells run at
≈ 1.3 CPU-h per unit c near c = 0.4). Above c = 0.56 the K = 3299 table is insufficient: use coeff_long_17200.txt (K ≤ 6544 at
c = 0.777 with the 1e−7 tail).

STATUS: READY FOR PROOF (on the critical path: it is the cheapest complete route found; nothing new beyond G1 + driver changes).

### G4 — designer addendum (high-c runs, h-pieces, unscaled variant timing)

* **[0.6, 0.6068], b = 0, both h-pieces PASS** (17,201-table, K = 3316, tail < 1e−7, allowance 1e−7 + 1e−14), my runs (non-official;
  the prover should rerun under its own driver): piece [0,½]H lower 9.8801e−6 (centre 1.6857e−4, L1 4.3994e−6), piece [½,1]H lower
  9.9645e−6 (centre 1.7000e−4, L1 4.4441e−6). Unscaled variant (D = 130): 159 s + 181 s; scaled variant with precision guard (D = 556):
  220 s for the lower piece. ⇒ above c = 0.5 the unscaled variant is ~30 % faster and needs no D growth; below c = 0.5 keep the scaled
  checker (2× faster). Results are identical in centre/L1/lower between the two variants.
* **h-piece choice at high c** (cost_model.py, piece_opt.py, c = 0.77, ω = 1.08, b = 0, est. seconds per cell): m = 2: 2785;
  m = 3 uniform: **1928**; m = 4 uniform: 2088; non-uniform [0,.5,.75,.9,1]: 2212. The *bottom* piece is the bad one when it is wide
  (|c − h₀p| does not cancel there), so prefer uniform m = 3 for c ≳ 0.7, m = 2 for 0.55 ≲ c ≲ 0.7.
* **Projected total cost (ω = 1.08, b = 0; b = 4 similar per cell but ω ≤ 1.06):** per channel ≈ 1.5 CPU-h on [0.268, 0.56] and
  ≈ 5.3 CPU-h on [0.56, 0.777] (per-unit-c costs 6.5, 9.8, 15, 29, 43, 66 CPU-h at c = 0.56, 0.6, 0.65, 0.7, 0.74, 0.77).
  **Whole band, both channels: ≈ 15–20 CPU-h** (≈ 5 h wall on 4 cores).

### G5 (designer): twist-direction structure on the physical domain; consequences for the author's degree-600 twist polynomial — STATUS: PROPOSED

**Statement.** For c ∈ [0.268, 0.788], b ∈ {0,4}, let P(c, y) := S₀(c, h = y·(1−√c)²/998²) (unnormalised prefix numerator, K ≥ tail rule),
y ∈ [0,1] (y = x/x_phys(c) with x = h/c; [0,1] is exactly the twist-allowed range, the spin cap is weaker for c ≥ 0.26).
Let t_m(c) be the Chebyshev coefficients of y ↦ P(c,y) on [0,1]. Then (i) |t₁| ≤ 0.03·t₀; (ii) Σ_{m≥2} |t_m| ≤ 10^{−8}·t₀; (iii) the
monomial Taylor coefficients P_m(c) := [y^m]P(c,y) satisfy max_m |P_m(c)| ≥ 10^{40}·t₀(c) (they are NOT small).

**Intended use.** Answers the coordinator's question about certifying the author's degree-600 twist polynomial:
* A dominant-term argument in the *monomial* x-basis (x⁰ term minus Σ_{m≥1}|p_m| x_max^m) is **hopeless**: by (iii) the monomial
  coefficients carry 40–90 orders of hidden cancellation inside the physical range (y-Taylor at c = 0.3: |P_m| reaches 2e37 at m ≈ 250;
  at c = 0.6: 9e85 at m ≈ 200; values 4e−5, 1.6e−4).
* The right basis is Chebyshev in y = x/x_phys(c) on [0,1]: by (i)–(ii) the dominant-term test t₀ − |t₁| − Σ_{m≥2}|t_m| > 0 has ≥ 97 %
  margin, and degree ≈ 100 in y already gives 1e−26 accuracy (the degree-600 representation is far more than needed). If the author's
  polynomial is supplied as exact coefficients p_m(c) (as polynomials/Taylor models in c on rectangles), the certificate is: per c-subcell,
  convert to y-Chebyshev on [0, x_phys(c_lo)] by exact rational arithmetic and apply the dominant-term test, with the c-dependence handled
  by a 1-D Taylor model in c — but the c-direction keeps the full 10^{430–630} cancellation of the ordinary series, so this is **not
  cheaper** than G4 (which already handles c and h jointly at ≈ 15–20 CPU-h for the whole band). Recommendation: do not build the
  twist-polynomial sign checker; use G4.
* The x⁰ term is the infinite-spin action S₀(c,0) = Σ A_k c^k; it is ≈ 98 % of the action on the physical domain (t₁/t₀ ≈ 1.3–2 %),
  so positivity is driven by it, but a rigorous bound on the twist correction needs the same 2-D machinery.

**Feasibility evidence.** code/gap_design/ycheb.py (c = 0.3, b = 0, 400 Chebyshev nodes): t₀ = 4.152e−5, t₁ = 7.73e−7, |t₂…t₅₀| ≤ 2e−14,
|t₇₅| = 4e−18, |t₁₀₀| = 1e−26, |t₂₀₀| = 1e−88 (super-geometric decay after m ≈ 60). box_cheb_y.py on the author's 21-rectangle c-grid
(10 rows from 0.268 to 0.7879, both b; 24×8 nodes; logs/ycells.log): L1/centre of the full 2-D Chebyshev expansion 0.12–0.25 (mostly
the c-linear term), y-linear term 1.3–2 % of centre, y-degree ≥ 4 tail ≤ 1.5e−11 absolute. xseries.py / xseries2.py: monomial y-Taylor
coefficients quoted above. scan_h2.py: deviation from affine in y is ≤ 3e−13 (relative 2e−8) on [0,1] at c = 0.3 and 0.5.

**Expected certificate cost.** Not recommended (see above); if used, same order as G4.

STATUS: PROPOSED (diagnostic; not on the critical path).

### DESIGNER NOTE (routes R2, R3 assessed; not proposed as lemmas)

* **R2, conformal variable w (G3).** Exact transformed coefficients Ã_j = Σ_k 4^k(−1)^{j−k}C(j+k−1, j−k)A_k from the 17,201-table
  midpoints (code/gap_design/conformal_sizes.py; the midpoint error after amplification is < 10^{−480} for j ≤ 1600, so the sizes are
  reliable): log10|Ã_j| = 336 (j = 100), 895 (400), 1130 (600), 1292 (800), 1382 (1000), then a plateau ≈ 1390 (j = 1100–1600)
  (the rigorous envelope is 2^{4813} ≈ 10^{1449}). The cancellation depth max_j |Ã_j| w(c)^j is ≈ 10^{683} at c = 0.5 (ordinary
  basis: 10^{533}) and ≈ 10^{947} at c = 0.777 (ordinary: 10^{626}). So the conformal prefix is shorter (J ≈ 3.1k vs K ≈ 6.5k at
  c = 0.777) but ~1.5× deeper in cancellation, i.e. Q (Taylor degree of the certificate) ~1.5× larger; cost ∝ J·Q² comes out the
  same as the ordinary basis, with extra work (transformed seed budget, five-term recurrence evaluation, a new table). **Not
  recommended.** G3 can stay PROPOSED/parked.
* **R3, partial fractions / small-h expansion (§123).** The twist dependence is NOT a small correction in any monomial/partial-fraction
  sense: on the physical range the prefix is affine in h to 1e−8 relative, but its h-Taylor coefficients (in units of the physical cap)
  reach 2e37 (c = 0.3) and 9e85 (c = 0.6) against values 4e−5 and 1.6e−4 (xseries.py, xseries2.py), and the function explodes
  just past the cap (G4 table). A first-order-in-h expansion with a Cauchy/Lagrange remainder therefore needs a bound on |S| on a complex
  h-disk, which the only available (termwise) majorant puts at 10^{500–700}: the remainder bound fails by > 500 orders. Partial
  fractions in h at the poles h = 1/(j(j+s)) have residues Π_{i≠j}(c − p_i/p_j)/(1 − p_i/p_j) of factorial size with alternating signs —
  the same cancellation in another guise. **Not recommended.** (The author's termination-preserving x = h/c coordinate is fine as a
  coordinate, but it does not remove the c-direction cancellation either; cost model route B ≈ route A per unit c, plus finite-spin checks.)
* **Summary ranking:** G4 (route A, physical-cap cells, ω ≤ 1.08, m h-pieces, 17,201-table above 0.56, unscaled checker above 0.5)
  ≈ 15–20 CPU-h total ≫ everything else. Route B (termination-preserving x-box + 10,214 finite-spin 1-D checks) is a valid fallback
  with similar cost. R2, R3, and the monomial twist-polynomial dominance are not competitive.

### G2 — designer addendum: the allowance to use for 0.56 < c ≤ 0.777 (recommended READY FOR PROOF as stated here)

**Statement.** For the author's corrected 17,201-pair vector (coeff_long_17200.txt, seed-substituted, spectral-degree M = 6000), on
physical states with c ≤ 0.9 and any prefix length K for which the exact ordinary tail tail(c_hi, K) (envelope 2^{2337},
tight_spectral_tail.py) is < 10^{−7}: |H_ang,b − normalised prefix| < 10^{−7} + 10^{−14}. Components (all in the package, exact
rational checks): ordinary tail < 10^{−7} (tight_spectral_tail.tail); old rational-pole seed correction < 1.8e−15 and log-seed
< 1e−20 (pole_bound_full.py / raw_seed_budget_full.py, c ≤ 0.9); new Gaussian M₀ seed < 1e−45 (verify_new_seeds.py); spectral-degree
truncation < 2^{−1910} (spectral_cut.py). These are exactly the numbers sign_next.py subtracts (allowance = 1e−7 + 1e−14 when
c_hi > 0.56). **Use:** G4 cells above c = 0.56. **Margin:** predicted lower bounds ≈ 1e−5 (G4), i.e. the allowance is 100× smaller.
**Cost:** seconds (rerun the four scripts; check they print PASS). Retained inputs: SPECTRAL_DEGREE.md §§1–3 and §118's seed proofs.

---
### G2 — PROVER: proof and check of the designer-addendum form (prover agent, 26 Nov)

I took the designer addendum ("recommended READY FOR PROOF as stated here") as the statement.

**Statement proved.** Fix b ∈ {0,4}. Take the 17,201-pair seed-substituted vector coeff_long_17200.txt, a prefix length K, and a cell with
c_hi ≤ 0.9. On physical states (h ≤ 998^{−2}, U_k ≤ c^k), |H_ang,b − normalised prefix_K| ≤ tail(c_hi, K) + 10^{−14}. Here tail is the exact
ordinary-tail bound of tight_spectral_tail.py, with envelope 2^{2337} on |z| = 199/200.

**Proof (assembly of retained, separately proved bounds).** The true action differs from the normalised prefix by six terms.
Each term majorises absolute per-k contributions with U_k ≤ c^k, h ≤ 998^{−2} and normalisation ≤ 1, so it holds for every prefix length.

(f) Ordinary truncation beyond K: ≤ tail(c_hi, K). This is the same estimate as EXTENSION.md "tight ordinary-power tail", with envelope
2^{2337} (SPECTRAL_DEGREE.md §5).

(e) Spectral-degree truncation M = 6000: < 2^{−1910} (SPECTRAL_DEGREE.md §§2–3, spectral_cut.py, c ≤ 0.9).

(d) New Gaussian M₀ seed substitution: < 10^{−45} (SPECTRAL_DEGREE.md §1, verify_new_seeds.py).

(a) Old rational-pole seeds. Each seed error is an imaginary atom at ν² = −a/4. Its action contribution is
|r_a|·(mass)·(exact polynomial action bound at c = 0.9)/998, using pole_bound_full_results.json. The masses are 7.1e−24/a from the cutoff
discrepancy plus 10^{−40} numerical width. The sum is < 1.7e−15 in each channel. The exact values are 1.6679e−15 for b = 0 and
1.6670e−15 for b = 4.

(b) Correlated log seed: < 10^{−21}, exact value 8.18e−22 (raw_seed_budget_full.py; its geometric sums use r = 0.9/0.99, i.e. c ≤ 0.9).

(c) Dropping θ on the log component: < 10^{−400}.

(a)+(b)+(c)+(d)+(e) < 1.7e−15 + 10^{−21} + 10^{−45} + 2^{−1910} + 10^{−400} < 10^{−14}. ∎

**Remark.** Below c = 0.56 my cells keep the K = 3299 table (old vector), need_k and the original 10^{−18} allowance (README §2). Above
0.56 they switch to the new vector, with allowance tail(c_hi, K) + 10^{−14} and K the least integer with tail(c_hi, K) < 10^{−7}. The two
vectors differ coefficientwise by coherent seed substitutions only (|ΔA_k| ≤ 10^{−107}, |ΔG_k| ≤ 10^{−1050} for k ≤ 3298;
code/gap_prove/compare_tables.py is the diagnostic). Each vector carries its own action-level budget, and the certificate uses one vector per
cell.

**Check** (code/gap_prove/check_g2.py PKGDIR, run on an unmodified copy of the package scripts): PASS. It runs tight_spectral_tail.py
(E = 2337), pole_bound_full.py, raw_seed_budget_full.py, verify_new_seeds.py and spectral_cut.py, and re-asserts (a)–(e) with exact
rationals. Example tails: tail(0.6, 3316) = 4.2e−24 and tail(0.7, 4700) = 5.4e−13.

**Retained inputs:** SPECTRAL_DEGREE.md §§1–3 (the envelope and the spectral-cut proof) and the §118 seed proofs.

STATUS: PROVED (+ check PASS), conditional on the retained inputs named above, which are the same as for the c ≤ 0.56 cells.

### G4 — designer addendum 2 (hardest cell verified; ω rule for b = 4; ready-made plan for [0.56, 0.777])

* **Top-of-band cell [0.77, 0.7781] (ω = 1.08), b = 4, m = 3 uniform h-pieces, K = 6677 (tail < 1e−7), 17,201-table, unscaled checker
  (D = 130): all three pieces PASS** — lowers 7.0366e−6, 7.0158e−6, 6.9521e−6 (centres 6.85e−4, L1 3.2–3.6e−5), 398 + 452 + 528 s.
  This overlaps the laptop range (c ≥ 1 − e^{−1.5} = 0.77687). Together with [0.6, 0.6068] b = 0 (addendum 1) and [0.3, 0.3097] b = 0,
  G4 is now verified by real certificate runs at the bottom, middle and top of the band, matching the sampling predictions.
* **ω rule for b = 4 refined** (omega_scan.py, logs/omega4.log; predicted L1/centre): ω = 1.06: 0.046 / 0.043 / 0.037 at c = 0.3/0.5/0.7;
  ω = 1.08: 0.064 / 0.058 / 0.055; ω = 1.10: 0.45 / 0.22 / 0.23 (under-resolved, near failure). ⇒ **ω ≤ 1.08 for both channels**
  (the G4 statement's 1.06 for b = 4 is conservative). This explains the prover's b = 4 lows at 0.022-wide cells (9.2e−7 at
  [0.44, 0.462], ω ≈ 1.11; FAIL at [0.462, 0.484], ω ≈ 1.12): choose widths by ω, not by a fixed width.
* **Ready-made plan:** code/gap_design/plan_highc_056_0777.json — 21 cells on the 1/2000 grid covering [0.56, 0.777] with ω ≤ 1.08,
  K from the exact 1e−7 tail, m = 2 (c < 0.7) / 3, estimated ≈ 10 CPU-h for both channels (unscaled checker ~30 % less). Driver:
  code/gap_design/sign_hsplit_u.py (unscaled; needs direct_sign_triangle_u.cpp, the author's tight_spectral_tail.py/base_bounds.py/
  resid_budget.py/check_lowfreq.py/band_masses.csv/residual_masses.csv in the working dir) — or the prover's own driver with the same
  parameters. Allowance: G2 addendum (1e−7 + 1e−14).
* Table consistency: prefix values from the K = 3299 table and from the 17,201 vector agree to > 15 digits at c = 0.3, 0.5, 0.55
  (table_xcheck.py), so cells below 0.56 certified with either table are interchangeable.

---
### G1 — cells: execution ledger, interim (prover agent, 26 Nov, 19:10)

**Code** (code/gap_prove/; the author's files are unmodified; the C++ binary is the author's direct_sign_triangle.cpp, compiled unchanged):
* `gap_triangle.py`: the author's triangle driver plus the G1 cap, equal h-pieces (m) and input precision D.
* `gap_stair.py`: **staircase** variant. [0, hcap] is split into mb equal pieces below h_u ≈ (1−√c_hi)²/998², plus mt pieces above it.
  The piece starting at h_i uses only c ∈ [c_lo, c_i], with c_i = ⌈(1−998·r_i)²⌉₁₀⁻⁶, where r_i ≤ √h_i. This is valid by G1:
  h ≥ h_i forces c ≤ (1−998√h_i)².
  * Effect: it removes nearly all of the non-physical overshoot (designer's ω), and L1 drops 5–10×.
  * Example: b = 4, [0.462, 0.484] FAILS as a plain m = 2 rectangle (−6.4e−5), but PASSES as a 2+2 staircase (1.04e−5).
* `gap_stair_long.py`: c ≥ 0.56. It uses coeff_long_17200.txt, K = least K with tail(c_i, K) < 1e−7 (envelope 2^2337), and allowance
  1e−7 + 1e−14 (G2).
* `run_chain*.py`: adaptive chains.
* D is set from the a-priori growth of the C++ l1 rounding recurrence. The reported roundoff is ≤ 6e−74 in every cell.
* `verify_gap_cells.py`: modelled on verify_extension.py. It checks exact transfer of every coefficient ball from the table, the prefix
  condition (need_k, or the tail < 1e−7), the G1 caps recomputed independently, exact tiling of [0, hcap] by the h-pieces, the staircase
  inequality, real and complex pole gaps, the Cauchy majorant and triangle tail, LOWER = CENTRE − L1 − ERROR, and the final rational
  inequality. It writes `ledger.json`.
* `cells.json`: the cell list.
* Decompressed tables and the .in/.out files are in the session scratchpad, not in the repo.

**Verified PASS cells so far** (verify_gap_cells.py: PASS, 27 channel cells). The lower bound is the normalised action lower bound, minimum
over pieces. K is the largest K over the pieces.

| c-cell | b | kind | K | pieces | lower | CPU s |
|---|---|---|---|---|---|---|
| [0.2800, 0.2900] | 0 | tri | 1512 | 1 | 1.3618e-05 | 49 |
| [0.2900, 0.3100] | 0 | tri | 1598 | 1 | 1.2840e-05 | 72 |
| [0.3100, 0.3300] | 0 | tri | 1688 | 1 | 1.2572e-05 | 87 |
| [0.3300, 0.3520] | 0 | tri | 1793 | 1 | 1.1939e-05 | 89 |
| [0.3520, 0.3740] | 0 | tri | 1903 | 1 | 1.1558e-05 | 99 |
| [0.3740, 0.3960] | 0 | tri | 2021 | 1 | 1.0968e-05 | 133 |
| [0.3960, 0.4180] | 0 | tri | 2146 | 1 | 1.0509e-05 | 184 |
| [0.4180, 0.4400] | 0 | tri | 2280 | 1 | 9.5265e-06 | 257 |
| [0.4400, 0.4620] | 0 | tri | 2425 | 2 | 8.5258e-06 | 176 |
| [0.4620, 0.4840] | 0 | tri | 2580 | 2 | 6.3341e-06 | 197 |
| [0.4840, 0.5060] | 0 | stair | 2749 | 2+2 | 1.0419e-05 | 285 |
| [0.5060, 0.5335] | 0 | stair | 2981 | 2+2 | 9.7684e-06 | 424 |
| [0.5335, 0.5600] | 0 | stair | 3230 | 2+2 | 9.4623e-06 | 639 |
| [0.5600, 0.5850] | 0 | long (G2) | 3087 | 2+2 | 9.0956e-06 | 520 |
| [0.2800, 0.2900] | 4 | tri | 1512 | 1 | 1.3075e-05 | 52 |
| [0.2900, 0.3100] | 4 | tri | 1598 | 1 | 1.1682e-05 | 73 |
| [0.3100, 0.3300] | 4 | tri | 1688 | 1 | 1.1599e-05 | 89 |
| [0.3300, 0.3520] | 4 | tri | 1793 | 1 | 8.5239e-06 | 91 |
| [0.3520, 0.3740] | 4 | tri | 1903 | 1 | 8.2361e-06 | 95 |
| [0.3740, 0.3960] | 4 | tri | 2021 | 1 | 7.4508e-06 | 133 |
| [0.3960, 0.4180] | 4 | tri | 2146 | 1 | 4.7262e-06 | 191 |
| [0.4180, 0.4400] | 4 | tri | 2280 | 1 | 3.6399e-06 | 259 |
| [0.4400, 0.4620] | 4 | tri | 2425 | 2 | 9.2279e-07 | 206 |
| [0.4620, 0.4840] | 4 | stair | 2580 | 2+2 | 1.0447e-05 | 228 |
| [0.4840, 0.5115] | 4 | stair | 2793 | 2+2 | 9.8117e-06 | 356 |
| [0.5115, 0.5415] | 4 | stair | 3053 | 2+2 | 9.2085e-06 | 558 |
| [0.5415, 0.5600] | 4 | stair | 3230 | 2+2 | 9.7207e-06 | 496 |

**Contiguous coverage: c ∈ [0.28, 0.56] in both channels (L ∈ [0.3285, 0.8210]); b = 0 continues to 0.585.** Each result is conditional on
the retained inputs of the direct method. Chains continue upward with long-table staircase cells.

### G4 — PROVER NOTE (prover agent)
I am executing G4's conclusion (positivity on the band, both b) with the staircase cell shape above, not with the literal "ω ≤ 1.08,
uniform m" layout. The staircase cover contains every physical state by G1, so it certifies the same conclusion. It is cheaper and has more
margin: L1/centre is about 6 % against 30–60 % for plain rectangles at ω ≈ 1.10.

Two deviations from G4 as stated:
* below c = 0.56 I use need_k (2^{−110}-type prefix) with the 1e−18 allowance and the K = 3299 table, rather than the 1e−7 tail;
* D is chosen from the same a-priori growth rule.

G4 status will be set when the chains finish. It will read PROVED over the covered range, or GAP for whatever range remains uncovered.

### G4 — designer addendum 3 (plan pre-validated)

All 102 (cell, h-piece, b) rectangles of code/gap_design/plan_highc_056_0777.json were pre-screened by tensor-Chebyshev sampling
(validate_plan.py → logs/validate_highc.log): predicted L1/centre ∈ [0.004, 0.053], i.e. no cell of the plan is expected to fail.
The prover's staircase driver (bottom pieces + thin top strips, seen in its run directory) is an improvement over the uniform
h-pieces here and is consistent with this lemma; either can be used. Designer processes: none left running.
