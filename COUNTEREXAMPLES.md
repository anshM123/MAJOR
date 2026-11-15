# COUNTEREXAMPLES, FAILED CONJECTURES, ADVERSARIAL TESTS

## 1. A natural lemma that is false as stated

**False lemma.** "Let `γ_1,…,γ_L` be cut classes and `f` a contraction map. The homotopy class of the
cut-and-paste output `∂U_r` is the same for every minimal-position arrangement of representatives of the
`γ_l`."

Why it fails: a triangle move replaces a triangular cell with bit string `x` by one with `x ⊕ 111` (on the
three curves involved). Take the coordinate function `m = f(·)_r` equal to 1 on `000, 100, 010, 001` and 0 on
`110, 011, 101, 111`. Before the move `U_r` is connected inside the disc (the triangle joins the three
weight-1 sectors) and `∂U_r` consists of three arcs, each going around a weight-2 sector. After the move the
three weight-1 sectors are disjoint and `∂U_r` consists of three arcs, each going around a weight-1 sector.
The six endpoints on the boundary of the disc are matched differently. A single coordinate of a contraction
map can have this pattern (the contraction property constrains only the weighted sum over `r`).

Consequence: the "obvious" proof (cut and paste separately in `X_L` and `X_R` and compare) does not work.
The proof in PROOFS.md never needs this lemma.

## 2. Functionals that are not of Crofton type do violate the inequalities

`F_max(A) = min_c max(ℓ_L(c), ℓ_R(c))` violates monogamy of mutual information in 923/3000 random cylinder
configurations and in 299/600 four-holed-sphere configurations, and the 5-party cyclic inequality in 559/3000
(EXPERIMENTS.md). So "any reasonable combination of the two chiral lengths with a common class" is not
enough; linearity in the measure is used.

## 3. Regularisation artefact

Comparing finite parts (cut-off term dropped) for an inequality that is not balanced at some endpoint gives
meaningless signs: strong subadditivity "failed" in 1830/3000 single-metric cylinder trials. With the cut-off
kept there are none. Theorem 1 is accordingly stated with a cut-off and an `o(1)`.

## 4. Invalid test set-ups caught (not counterexamples)

* Random Schottky groups for the four-holed sphere: essential curves entered a minimal cut in ≤ 7% of trials,
  so "zero violations" was nearly vacuous.
* Wrong-side Fenchel–Nielsen gluing: non-discrete representation (trace 1.195 element); the single-metric
  control then showed 325/500 spurious violations.

## 5. First-phase results (quantum information; archived in archive_qi/)

* Ordinary-Petz recovery bound for conditional mutual information: violated. Independently rediscovered here
  (ratio `I/(−log F)` down to ≈ 0.69 at second order for qubits, verified with an independent 80-digit
  implementation), but **already known**: certified counterexamples were recorded in September 2026, and the
  general-channel form was disproved in arXiv:2309.02074. No novelty claimed.
* Zhang–Yeung inequality for von Neumann entropies: no violation found by three different searches. Open.
* Transpose-degradable versus degradable channels: no separation found by SDP in dimensions 2×2, 2×4, 4×2,
  3×3. Suggests, but does not prove, that transpose-degradable implies degradable.
