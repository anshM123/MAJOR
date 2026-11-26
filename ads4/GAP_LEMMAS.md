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
