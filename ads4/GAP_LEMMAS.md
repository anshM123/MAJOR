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
