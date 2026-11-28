# FULLCONE: the full RT cone for arbitrary parties (intervals with endpoints + whole boundaries)

Working notes, appended in the order the work was done. Nothing outside `fullcone/` was edited.
Notation is that of PROOFS.md §§1–5 (`X = H²/Γ` with funnel ends, marked points `P`, pieces, cut classes
`𝒞(A)`, test classes, `i(·,·)`, horoball cut-offs `B_p(ε)` normalised as in PROOFS §3.2, `e_p(c) ∈ {0,1}` the
number of arc-ends of `c` at `p`, `W = Σ_j w_j`).

Grades: **[proved]** = complete argument here, using only PROOFS §§1–4 and standard facts named on the spot;
**[proved mod X]** = complete given the cited statement X; **[conjecture]**; **[false]**.

---

## 0. Summary of what was found (filled in at the end; see §9)

---

## 1. The key new identity: renormalised Crofton formula (replaces the `κ`-asymptotics of Lemma B)

### 1.1 The infinite-mass classes

For `p ∈ P` let `τ_p` be the test class of the short arc that runs from the piece on one side of `p` to the
piece on the other side, around `p` (it cuts off a half-disc of `X̄` containing `p` and no other marked point,
no limit point). In the chart of PROOFS §3.2 (`p̃ = ∞`, `B_p̃(ε) = {y > 1/ε}`) let

`q_- := min(Λ ∪ Γ·P̃ \ {∞})`, `q_+ := max(Λ ∪ Γ·P̃ \ {∞})`, `d_p := q_+ − q_- > 0`

(the set is closed and bounded, PROOFS §3.2; `q_±` are the lifts of the marked points adjacent to `p` along its
boundary circle, possibly translates of `p` itself). We assume `q_- < q_+`, which fails only for `X = H²` with `|P| ≤ 2`
(a single interval in pure AdS3, `N = 1`, trivial); then also `τ_p ≠ τ_{p'}` for `p ≠ p'`. The two lifted pieces adjacent to `∞` are `(q_+, ∞)` and
`(−∞, q_-)`, and the geodesics of class `τ_p` are, in this chart and up to `Γ`, exactly

`τ̃_p = {g = (a, b) : a < q_-, b > q_+}`   (stabiliser of `p̃` is trivial, so this is a fundamental set).

**Lemma 1.1 [proved].** For every cut class `c` (geodesic representative `c^X`) and every `ĝ ∈ τ_p` with ends
off `Γ·P̃`: `#(ĝ ∩ c^X) = e_p(c)`, and the unique crossing (if `e_p(c)=1`) is with the lift `x = x_0` of the
arc-end of `c` at `p` (`x_0 ∈ [q_-, q_+]`), at Euclidean height `√((x_0 − a)(b − x_0))`.

*Proof.* Lifts of components of `c^X` have endpoints in `Λ ∪ Γ·P̃`. Those not ending at `∞` have both
endpoints in `[q_-, q_+] ⊂ (a, b)` and are not linked with `g`. The lifts ending at `∞` are in bijection with
the arc-ends of `c` at `p` (trivial stabiliser of `p̃`), so there are `e_p(c) ≤ 1` of them (condition 1 of
PROOFS §1.2); such a lift is a vertical line `x = x_0`, `x_0 ∈ Γ·P̃ ⊂ [q_-, q_+]`, linked with `g` once. ∎

**Lemma 1.2 [proved].** Every test class `τ` other than the `τ_p` and the classes homotopic into
`∂_∞X \ P` ("trivial") has finite `m̄`-mass `m̄(τ) > 0`. Closed test classes have mass `0`.

*Proof.* A closed class is one point of `𝒢/Γ`. Let `τ` be an arc class with lifted end-pieces `I ∋ a`,
`J ∋ b`. (i) If `I, J` lie in different components of `Ω̃ = S¹ \ Λ`, a lift crosses the axis bounding the
component of `I`; with `β` the corresponding boundary geodesic, `m̄(τ) i(τ, β) ≤ ℓ(β) < ∞` (PROOFS §5.3). (ii)
If `I, J` lie in one component `(u, v)` of `Ω̃`, the side of the chord `ab` inside `(u,v)` contains no limit
point and contains `k` lifted marked points: `k = 0` is trivial, `k = 1` is `τ_p`; for `k ≥ 2` the closures of
`I` and `J` are disjoint (they are separated by the `≥ 2` marked points on that side, and on the other side by
`u, v`, or, for trivial `Γ`, by further marked points unless the other side has exactly one marked point `q`,
in which case `τ = τ_q`), so `m(I × J) = ∫∫_{I×J} da db/(b−a)² < ∞`; the stabiliser of the pair is trivial.
Positivity: `I × J` is open. ∎

Let `T := ⋃_p τ_p` (as a subset of `𝒢/Γ`) and `m̄' := m̄|_{(𝒢/Γ)\T}`. For a cut class `c` put
`i(m̄', c) := ∫ i(τ(ĝ), c) dm̄'(ĝ) = Σ_τ m̄'(τ) i(τ, c)` (a priori in `[0, ∞]`).

### 1.2 The identity

**Proposition 1 (renormalised Crofton) [proved].** For every cut class `c` and every `ε ≤ ε_0`,

`ℓ^ε(c) = i(m̄', c) + Σ_{p∈P} e_p(c) · log(1/(ε d_p))`,   (R)

and `i(m̄', c) = ℓ^{ε_0}(c) + Σ_p e_p(c) log(ε_0 d_p) < ∞`.

*Proof.* Let `c^ε = c^X \ ⋃_p B_p(ε)`, a finite union of compact geodesic segments and closed geodesics. By
Crofton (PROOFS §3.1), `ℓ^ε(c) = ∫_{𝒢/Γ} #(ĝ ∩ c^ε) dm̄ = I_1(ε) + Σ_p I_p(ε)` with
`I_1(ε) = ∫_{(𝒢/Γ)\T} #(ĝ ∩ c^ε) dm̄` and `I_p(ε) = ∫_{τ_p} #(ĝ ∩ c^ε) dm̄`.

*`I_1`.* As `ε ↓ 0`, `c^ε ↑ c^X`, so `#(ĝ ∩ c^ε) ↑ #(ĝ ∩ c^X) = i(τ(ĝ), c)` for `m̄`-a.e. `ĝ` (B0), and by
monotone convergence `I_1(ε) ↑ i(m̄', c) ∈ [0, ∞]`.

*`I_p`.* If `e_p(c) = 0`, `I_p = 0` by Lemma 1.1. Otherwise unfold to the fundamental set `τ̃_p`. By Lemma 1.1
the only crossing is with `x = x_0`, at height `h(g) = √((x_0−a)(b−x_0))`. That line is a lift of `c^X`, which
meets the horoball family only in its two terminal rays (Fact 3.1): the ray to `∞` (removed iff `h > 1/ε`)
and the ray to the marked point `x_0` (a horoball tangent at `x_0` of Euclidean diameter `η_ε → 0`). So
`I_p(ε) = m{g ∈ τ̃_p : h(g) < 1/ε} − m{g ∈ τ̃_p : h(g) < η_ε}`. The second term tends to `0` (decreasing sets
of finite measure with empty intersection, since `h > 0` on `τ̃_p`). For the first, with `u = x_0 − a > u_0
:= x_0 − q_-`, `v = b − x_0 > v_0 := q_+ − x_0`, `Y = 1/ε`, `u_0 + v_0 = d_p`:

`m{uv < Y²} = ∫_{u_0}^{Y²/v_0} [1/(u+v_0) − u/(u²+Y²)] du`
`= log((Y²/v_0 + v_0)/d_p) − ½ log((Y⁴/v_0² + Y²)/(u_0² + Y²)) = log Y − log d_p + o(1)`

(for `v_0 = 0` the same limit, by the same computation with the upper limit `∞`). **The limit does not depend
on `x_0`**, i.e. on `c`. Hence `I_p(ε) = e_p(c)[log(1/ε) − log d_p] + o(1)`.

*Exactness.* By (3.1), `ℓ^ε(c) − e(c) log(1/ε) = ℓ^{ε_0}(c) + e(c) log ε_0` is independent of `ε ≤ ε_0`,
and `e(c) = Σ_p e_p(c)`. Therefore `I_1(ε) = ℓ^{ε_0}(c) + e(c) log ε_0 + Σ_p e_p(c) log d_p + o(1)` is
bounded, so its monotone limit `i(m̄', c)` is finite and equals `ℓ^{ε_0}(c) + Σ_p e_p(c) log(ε_0 d_p)`. Then
`ℓ^ε(c) = ℓ^{ε_0}(c) + e(c) log(ε_0/ε) = i(m̄', c) + Σ_p e_p(c) log(1/(ε d_p))` exactly. ∎

*Consistency with Lemma B.* The `τ_p`-geodesics that lie in `𝒢_ε` (radius `< 1/ε`) have mass
`∫_{d_p}^{2Y} (s − d_p)/s² ds = log(2Y/d_p) − 1 + o(1)`, while the `τ_p`-crossings of `c^ε` have mass
`log(Y/d_p) + o(1)`; the difference is `1 − log 2 = κ`, which is exactly Lemma B's end constant. (Checked
numerically, §6.)

*Remarks.* (a) The identity holds for each class separately; no uniformity in `c` is used. (b) Only the
truncated part `I_1` converges to an intersection pairing; the divergence is carried entirely by the single
class `τ_p`, whose intersection number with **every** cut for `A` is `[p ∈ ∂A]`, a constant on `𝒞(A)`. That is
the whole mechanism of the extension to `P ≠ ∅`. (c) Changing the horoball normalisation at `p` by a factor
changes `d_p` by the same factor; (R) is covariant.

---

## 2. Consequence for the multi-metric functional

Setting of PROOFS §4: `(X_j, P_j)`, `φ_j`, weights `w_j > 0`, common `ε_0`. Transport each `m̄'_j` to test
classes of `X_1` and put

`ν := Σ_j w_j m̄'_j` (a countable positive combination `Σ_τ ν(τ) τ` over arc classes `τ ∉ {τ_p} ∪ trivial`,
with `0 < ν(τ) < ∞`),
`t_p(ε) := Σ_j w_j log(1/(ε d_p^{(j)}))`,  `G(A) := min_{c ∈ 𝒞(A)} i(ν, c)`.

(`τ_p` and "trivial" are topological notions, preserved by `φ_j`.)

**Corollary 2 [proved].** For every region `A` and `ε ≤ ε_0`:

`F^ε(A) = G(A) + Σ_{p∈∂A} t_p(ε)`,  i.e.  `F^ε = G + Σ_p t_p(ε) b_p`,   (D)

where `b_p(A) := [p ∈ ∂A]`. The minimum defining `G(A)` is attained, at the (finite, `ε`-independent) set of
minimisers of `F^ε(A)`.

*Proof.* By (R) for each `j`, `Σ_j w_j ℓ^ε_{X_j}(c) = i(ν, c) + Σ_p e_p(c) t_p(ε)`, and `e_p(c) = [p ∈ ∂A]`
for every `c ∈ 𝒞(A)` (PROOFS §1.2, condition 1). So on `𝒞(A)` the two functionals differ by a constant. ∎

`b_p` is the entropy vector of one Bell pair shared by the two parties adjacent at `p` (if one of them is
the purifier, `b_p(A) = [i ∈ A]`; if both pieces have the same party, `b_p = 0`). It is the min-cut function
of the graph with one edge of weight 1 between the two terminals: a graph-model vector.

---

## 3. Theorem 3 for arbitrary parties — exact graph realisation (route R2, made rigorous)

**Theorem 3* [proved, given PROOFS §§1–4 incl. the surgery step of Lemma A].** `G` is exactly the min-cut
vector of a finite weighted graph with terminals labelled by parties. Hence:

1. for `ε ≤ ε_1 := min_p Π_j (d_p^{(j)})^{−w_j/W}` (so that all `t_p(ε) ≥ 0`), `F^ε` itself is the min-cut
   vector of a finite graph (the graph for `G` plus, for each `p`, an edge of weight `t_p(ε)` joining the
   terminals of the two parties adjacent at `p`), so **`F^ε` lies in the RT holographic entropy cone `C_N`**;
2. for every `ε ≤ ε_0` and every linear inequality `Q(S) := Σ_A q_A S(A) ≥ 0` valid on `C_N`,

   `Q(F^ε) = Q(G) + Σ_p t_p(ε) Δ_p(Q)`, with `Q(G) ≥ 0` and `Δ_p(Q) := Q(b_p) ≥ 0`;

   so `Q(F^ε) ≥ 0` exactly (and `= Q(G)`, independent of `ε`) if `Q` is balanced at every marked point, and
   `Q(F^ε) − Q(G) = W·(Σ_p Δ_p) log(1/ε) + const → +∞` otherwise.

No contraction map is used: this holds for every valid RT inequality, whether or not it admits one.

*Construction.* Fix for every region `A ≠ ∅` a minimiser `c_A` of `G(A)` (Corollary 2), in `X_1`-geodesic
position, and let `𝒜 := ⋃_A c_A^{X_1}`, a finite union of closed geodesics and complete geodesic arcs ending
at marked points. As in Lemma A, vertices are crossing points of distinct geodesics of `𝒜` (finitely many,
PROOFS §1.3(c)); *edges* are the segments of `𝒜` between consecutive vertices, plus the terminal rays from
the last vertex to a marked point, plus vertex-free closed geodesics; *cells* are the components of
`X_1 \ 𝒜` (finitely many: `𝒜 ∪ ∂_∞X` is a finite graph in the compact surface `X̄_1`).

*Pieces and terminals.* Every elementary piece lies in the closure of exactly one cell (closed geodesics stay
in the convex core; arcs approach `∂_∞X` only at marked points). Two pieces of different parties `i ≠ j` lie
in different cells: one of them, say `i`, is not the purifier, and `c_{{i}}` separates them (bit strings as in
Lemma A). So each cell is either *internal* (contains no piece) or a *terminal* of a unique party.

*The measure on geodesics.* For each arc class `τ` in the support of `ν` let `π_τ := m̄_1|_τ / m̄_1(τ)`, a
probability measure on the `X_1`-geodesics of class `τ` (Lemma 1.2: `0 < m̄_1(τ) < ∞`), and
`ν̂ := Σ_τ ν(τ) π_τ`, a measure on `𝒢/Γ` of `X_1`. It is absolutely continuous w.r.t. `m̄_1`, so `ν̂`-a.e. `ĝ`
has ends off `Λ ∪ Γ·P̃`, is not a geodesic of `𝒜`, misses all vertices, and is transverse to `𝒜`. By (B0)
(valid for **every** geodesic with ends off `Λ ∪ Γ·P̃`, PROOFS (1.1)),

`∫ #(ĝ ∩ c^{X_1}) dν̂ = Σ_τ ν(τ) i(τ, c) = i(ν, c)` for every cut class `c`.   (3.1)

No boundary map between `X_1` and `X_j` is needed: the metric dependence sits in the numbers `ν(τ)`.

*Weights.* `w(e) := ν̂{ĝ : ĝ crosses e}` (with multiplicity). Each edge lies on some `c_A`, so
`w(e) ≤ i(ν, c_A) < ∞` by (3.1). The graph `Γ_𝒜` has the cells as vertices and, for each edge `e`, an edge of
weight `w(e)` between the two cells on its sides (distinct cells, since `e` lies on some `c_A` whose two
sides differ in bit `A`). Zero-weight edges may be deleted.

*Claim: `mincut_A(Γ_𝒜) = G(A)` for every `A`.*

(≤) Let `S` be the set of cells contained in the homology region `W_A` of `c_A`. It contains the terminals of
the parties in `A` and no others. An edge separates `S` from its complement iff it lies on `c_A`. So the cut
weight is `Σ_{e ⊂ c_A} w(e) = ∫ #(ĝ ∩ c_A) dν̂ = i(ν, c_A) = G(A)` by (3.1).

(≥) Let `S` be admissible for `A` (terminals of parties in `A` in `S`, other terminals not). Let `U` be the
closure of the union of the cells in `S`. Its frontier `∂U` is the union of the edges across which membership
changes. Apply the surgery of Lemma A (PROOFS §2, "Surgery", with `U_r := U`, `B_r := A`): resolve at vertices,
pair rays inside small horoballs at marked points. Since the outermost sectors at `p` contain the two
adjacent pieces, the result `∂U'` is a cut for `A`, of some class `δ ∈ 𝒞(A)`. For `ν̂`-a.e. `ĝ` (transverse,
missing vertices) choose the discs and horoballs of the surgery so small that `ĝ` misses them; then
`i(τ(ĝ), δ) ≤ #(ĝ ∩ ∂U') = #(ĝ ∩ ∂U) = Σ_{e ⊂ ∂U} #(ĝ ∩ e)`. (The class `δ` does not depend on the sizes.)
Integrating, `G(A) ≤ i(ν, δ) ≤ Σ_{e ⊂ ∂U} w(e) = weight(S, S^c)`. ∎ (Claim)

*Proof of 1 and 2.* Merge the terminals of each party (identifying same-label terminals does not change any
min-cut, since they are always on the same side). Adding, for each `p`, an edge of weight `t_p(ε) ≥ 0`
between the terminals of the two parties adjacent at `p` adds exactly `t_p(ε)[p ∈ ∂A]` to every cut
separating `A` from its complement, hence to every min-cut. So `F^ε = G + Σ t_p b_p` is a graph-model vector
for `ε ≤ ε_1`. Part 2: `Q(G) ≥ 0` and `Q(b_p) ≥ 0` because `G`, `b_p` are graph-model vectors; linearity and
(D). The divergence: `t_p(ε) = W log(1/ε) + const`. ∎

*Remarks.*
* The `ε_1` threshold is convention-dependent (horoball normalisation), the statement "for all sufficiently
  small `ε`" is not. With the "adapted" normalisation `B_p(ε)` rescaled by `d_p^{(j)}` one has `ε_1 = 1`.
* This is literally the argument of Bao et al. (arXiv:1505.07839, Lemma 3: cells of the union of all RT
  surfaces) run in the **virtual static geometry** `(X_1, ν̂)`: a curve `σ` transverse to `𝒜` gets the
  "length" `ν̂{ĝ crossing σ}`, the `X_1`-geodesic representative of a class is length-minimising in its class
  (by (1.1)), and min over cuts of this length is `G`. It is the Burger–Iozzi–Parreau–Pozzetti /
  Martínez-Granado–Thurston pseudo-metric of the current `ν̂`, `d(x,y) = ½(ν̂(𝒢[x,y)) + ν̂(𝒢(x,y]))`
  (MG–T arXiv:2605.04031 eq. (8.4), Corollary B, Remark 8.9 "straight"), with the infinite-mass classes `τ_p`
  removed and accounted for by the Bell-pair edges. See §5 for what is and is not taken from MG–T.
* The same proof with `P = ∅` reproves Theorem 3 of PROOFS §5.3 **exactly** (no limit, no Lemma F, no
  closedness of the cone), and with Bao et al.'s "graph models obey contraction-map inequalities" reproves
  Theorem 1. The surgery step of Lemma A is the only topological input.

---

## 4. Route R1: the "six missing steps" of PROOFS §5.4, filled (independent second proof of `G ∈ C_N`)

PROOFS §5.4 pushed forward `Σ_j w_j m̄_j|_{𝒢_ε}` (an `ε`-dependent measure). With Proposition 1 one pushes
forward the **`ε`-independent** measure `ν` of §2 instead, and all six items become short.

Enumerate the support of `ν` as `τ_1, τ_2, …`, starting with the filling family `T_0` of Lemma F (if `X̄` is a
disc, `T_0 = ∅`). Put `μ_n := Σ_{k≤n} ν(τ_k) τ_k` and `G_n(A) := inf_{c∈𝒞(A)} i(μ_n, c)`.

1. *Finite masses* [proved]: Lemma 1.2 (each `τ_k` has `0 < ν(τ_k) < ∞`; the infinite-mass classes are
   exactly the trivial ones and the `τ_p`, and the latter are not in the support of `ν`).
2. *Existence of minimisers* [proved]: the proof of Lemma G shows `inf = min-cut` and the cut `∂W` built from
   an optimal vertex set attains it (its `i(μ_n,·)` is `≤` the min-cut weight). Lemma G applies verbatim with
   `P ≠ ∅`; its hypothesis "not homotopic into `∂_∞X`" is only used to exclude classes with `i ≡ 0`, and the
   classes enclosing `≥ 2` marked points are handled by the same proof (their geodesic representatives meet
   `∂_∞X` only at their two ends, so the thin collar of `A` can be chosen to miss them).
3. *Confinement, single stage* [proved]: a minimiser `c_n` of `G_n(A)` has
   `Σ_{τ∈T_0} ν(τ) i(τ, c_n) ≤ G_n(A) ≤ G(A)`, so by Lemma F and `min_{T_0} ν > 0` the `c_n` range over a
   finite set `𝒟_A` independent of `n` (no `ε` is involved any more).
4. *Uniformity of the end estimate*: **not needed**. Proposition 1 is an identity class by class.
5. *Positive mass of `T_0`* [proved]: each class of `T_0` contains an open set of geodesics of every `X_j`,
   so `ν(τ) > 0`; there is no `𝒢_{ε_0}` restriction to worry about.
6. *Explicit bound replacing `O(1)`*: **not needed**; (D) is exact.

Then `G_n(A) = min_{c∈𝒟_A} i(μ_n, c) ↑ min_{c∈𝒟_A} i(ν, c) ≥ G(A) ≥ G_n(A)`, so `G_n → G` pointwise; each
`G_n ∈ C_N` by Lemma G; `C_N` is closed (polyhedral, Bao et al. Prop. 7 / Thm 5). So `G ∈ C_N`
[proved mod closedness of `C_N`, cited], and the conclusions of Theorem 3* follow from (D) exactly as in §3.
This route needs Lemma F (normal coordinates, standard) and the closedness of `C_N`; route R2 (§3) needs
neither and gives an exact finite graph.

---

## 5. Route R2 and Martínez-Granado–Thurston (arXiv:2605.04031): what is used, what is not

Read: abstract, §1 (Theorem A, Corollary B, Corollary C), §8.4 (proof of Corollary B, eq. (8.4), Remarks
8.9–8.11, Question 8.12), §11.1 (open surfaces). Findings:

* MG–T work on **closed** orientable hyperbolic surfaces. §11.1 states that extending to open surfaces
  (boundary, punctures) is open: they formulate Conjecture 11.1 (a continuous intersection pairing
  `Curr_cusp × Curr_open`), Question 11.2 and Conjecture 11.3, and give an example showing that additivity +
  smoothing are **not** sufficient for a dual current in the open setting. So Theorem A (functional ⇒
  current) is **not available** for funnels, arcs to marked points, or cut-offs, and **we do not use it**.
* We do not need it: our functional comes with its current. `ν̂` is explicit (§3). The direction we use is
  current ⇒ pseudo-metric, which is elementary: MG–T eq. (8.4), due to BIPP21 §4, `d(x,y) =
  ½(μ(𝒢[x,y)) + μ(𝒢(x,y]))`, "straight" along geodesics (Remark 8.9). MG–T Remark 8.11 warns `d` need not be
  continuous or positive in general (it is continuous iff `μ` has no atoms, DRMG25); our `ν̂ ≪ m̄_1` has no
  atoms, but we never need continuity: only crossing measures of the finitely many edges of `𝒜`.
* The static-RT-to-graph step (Bao et al. Lemma 3) needs, of the "metric", only: (a) additivity of length over
  the edges of an arrangement, (b) minimality of the chosen representatives within their class, (c) that
  resolving an arrangement at vertices does not increase length. For `ν̂`-length these are (3.1), (1.1) and
  the transversality argument in the (≥) step. So the "virtual static metric" route works on surfaces with
  funnel ends, arcs to marked points and cut-offs, **provided the infinite-mass classes `τ_p` are split off
  first** (Proposition 1). Without that split the terminal rays have infinite `ν̂`-length and the graph is
  not defined; that, not anything about MG–T, was the obstacle in PROOFS §5.4.
* Verdict on R2: **works** [proved], in the form of §3. It is also the cleanest presentation of Theorem 1:
  "`F^ε` is the RT functional of the pseudo-metric of the current `Σ_j w_j L_{X_j}` (with the endpoint classes
  split off as Bell pairs), hence a graph model". The paper should cite MG–T for Corollary B / eq. (8.4)
  (and BIPP21 §4 for the construction), and say explicitly that MG–T's Theorem A is for closed surfaces and is
  not used; MG–T do not discuss entropy, cuts, homology regions or arcs, so they do not scoop Theorem 3*.
  (BIPP21 itself was not read; cited via MG–T.)

---

## 6. Numerical checks, part 1: facet lists and the simply connected case

All scripts: `python3 -I fullcone/<script>.py` (they add their own directory to `sys.path`).

**Facet lists (`cones.py`, `check_facets.py`).** `C_3`: SA (6) + MMI (1) = 7 facets; `C_4`: SA (10) + MMI
(10) = 20; `C_5`: SA 15, MMI 20, MMI-with-a-pair 45, cyclic 72, and the four further orbits 90, 10, 60, 60:
**372 facets**, matching the orbit sizes stated by Avis–Hernández-Cuenca (arXiv:2102.07535, §6 and Table 2).
Checks: every orbit representative except one admits a weighted contraction map (two independent searches,
`fc_common.contraction_exists` and `contraction_mrv.py`); 200 random graph models per `N` satisfy all
facets (min value `−1.5·10⁻¹⁵` relative); Bao et al.'s non-holographic ray (`S=1` for `|I|∈{1,5}`, `2`
otherwise) violates the facets (value `−2`).

*Side finding (literature-consistent, relevant to G4).* The fifth new `C_5` orbit (Avis–HC Table 2 no. 5 =
Bao et al.'s last inequality) has **no weighted contraction map** in the sense of PROOFS §1.5 (both searches
terminate with "none" after constraint propagation). Avis–HC (their eq. (4.19) and text) report the same and
prove it by splitting the RHS terms with coefficient ≥ 2 into unit-coefficient copies. Consequence for the
draft: Theorem 2 as stated ("admits a contraction map", weighted Hamming norms) does **not** literally cover
this facet of `C_5`. Lemma A works verbatim for the split form (repeated `B_r` are allowed), so the draft's
Theorem 2 should say "admits a contraction map after splitting repeated terms" — or, better, be replaced by
Theorem 3*, which needs no contraction maps at all.

**Simply connected case, exact (`disc_exact.py`).** `X = H²`, `k ∈ {1,2,3}` metrics = `k` configurations
of the same cyclically ordered marked points (pure AdS3 is `k = 2`, `w = ½`), random horoball normalisations,
random party assignments with `n ∈ [N+2, N+4]` pieces. Results (seed 7):

| N | trials | max `|G − Σ_τ ν(τ) i(τ,·)|` (identity (R)+(D)) | max `|G − mincut(cell graph)|` (Thm 3*) | min facet `F^ε` | min facet `G` |
|---|---|---|---|---|---|
| 3 | 300 | 4.7e-14 | 4.8e-14 | −4e-16 | −7e-15 |
| 4 | 200 | 5.3e-14 | 5.3e-14 | −5e-16 | −6e-15 |
| 5 | 60 | 5.0e-14 | 4.8e-14 | −1e-15 | −1.5e-14 |

Here `G` is computed as `F^ε − Σ_{p∈∂A} t_p`, independently as `min_c Σ_τ ν(τ) i(τ,c)` with `ν(τ)` the
Liouville masses (log cross-ratios) of pairs of non-adjacent pieces, and independently as the min-cut of the
explicit cell graph of §3 (arrangement of the minimisers in the Klein model of metric 1, `ν̂` = Dirac masses
on one generic chord per class). All three agree to rounding; all 372/20/7 facets hold for `F^ε` and for `G`
(note `G` must also satisfy the **unbalanced** facets SA/Araki–Lieb, which is a non-trivial check since
`G` has no divergent terms).

*Multiply connected checks of Proposition 1 itself (`cylinder_identity.py`).* Cylinder
`H²/⟨z ↦ e^ℓ z⟩`, 2–5 marked points on one boundary, random horoball normalisations; for single arcs
`(x_i, λ^m x_k)` winding `m ∈ {0,1,2}` times, `i(m̄', c)` computed as the explicit sum over test classes
(pairs of lifted pieces, log cross-ratios, plus the classes running to the other funnel), truncated at
40 periods (largest truncated mass `< 10⁻¹⁶`): 181 cases, max `|ℓ^ε(c) − i(m̄',c) − Σ_p e_p log(1/(ε d_p))|
= 5.7·10⁻¹⁴`. Core geodesic: `Σ_τ m(τ) i(τ, core) = ℓ` to `4·10⁻¹⁶` (Crofton sanity check).

---

## 7. Gravity: Theorem 3'' (HRT, arbitrary parties) and the cone equality

### 7.1 Theorem 3'' [proved, given Theorem H as in PROOFS §6 (acausal-slice definition)]

Let `M` be a vacuum 2+1 spacetime of the class of PROOFS §6.1 (any topology with funnel ends; rotation
allowed), parties finite unions of intervals of one acausal boundary Cauchy slice and of whole boundaries,
`S^ε` the HRT entropy of PROOFS §6.3. By Theorem H, `4G_N S^ε = F^ε` with `k = 2`, `w_L = w_R = ½`. Then:

* `4G_N S^ε = G + Σ_p t_p(ε) b_p` exactly, with `G` the min-cut vector of a finite graph and
  `t_p(ε) = log(1/ε) − ½ log(d_p^L d_p^R)`: **the regulated HRT vector is a graph-model vector plus one Bell
  pair across each entangling point, of strength `(1/4G_N) t_p(ε)`** (`= (c/6) log(1/ε) + O(1)`).
* For every `ε` below an explicit threshold, `S^ε ∈ C_N` (the RT holographic entropy cone).
* Every linear inequality valid on `C_N` holds for `S^ε` at every `ε ≤ ε_0`: with equality of the slack to
  the `ε`-independent number `Q(G)/4G_N ≥ 0` if it is balanced at every endpoint, and with slack
  `→ +∞` otherwise. **No contraction map, and hence no completeness claim about contraction maps
  (Bao–Furuya–Naskar), is used.** This removes OPEN_GAPS G4 for vacuum 2+1, and item 19 of G6a.

The caveats of Theorem H are inherited unchanged (acausal, possibly non-smooth, asymptotically lightlike
homology slices when the regions have endpoints; smooth slices only for whole-boundary parties).

### 7.2 Cone equality [proved mod Bao et al. arXiv:1505.07839 §3.1 (entropy-preserving graph moves), Lemma 4, Theorem 5, Proposition 7]

Let `𝓗_N` be the set of vectors `4G_N S` over all spacetimes of the class above (any topology, any rotation),
all party assignments with `N` parties, and (for intervals) all cut-offs `ε ≤ ε_1(M)`. Then

`cl(R_{>0} · 𝓗_N) = C_N`, and already `cl(R_{>0} · 𝓗_N^{static, whole boundaries}) = C_N`.

*⊆:* Theorem 3'' (and `C_N` is a closed cone). *⊇:* Bao et al. Lemma 4 (read in full, §3.2): for a graph
model with positive weights, trivalent bulk vertices and each boundary colour on a single degree-one vertex
(any graph can be brought to this form by their §3.1 moves without changing entropies), glue hyperbolic
pairs of pants and half-collars with gluing geodesics of length `4G_N w(e)`; for curvature sufficiently
negative (equivalently, at curvature `−1`, after scaling all weights by a small `λ > 0`) the RT entropies of
the boundary components equal the min-cuts. The surface (with the half-collars continued to funnels) is
`X = H²/Γ`, Γ convex cocompact, and the time-symmetric vacuum spacetime with `X_L = X_R = X` has, by
Theorem H (or directly, since it is static), `4G_N S = λ·(graph vector)`. Graph vectors of connected graphs
are dense in `C_N` (`C_N` is the set of graph vectors, Bao et al. Theorem 5/Prop. 7; a disconnected graph is
the limit of connected ones joined by edges of weight `→ 0`). Hence `C_N ⊆ cl(R_{>0} 𝓗_N^{static})`.

Remarks. (i) The equality is of *closed cones of rays*; with `G_N` fixed only small entropy vectors of a ray
are realised by the Bao et al. construction (lengths `4G_N λ w(e)` must be short). (ii) The extreme rays of
`C_N` are in particular realised (up to scale) by static multiboundary wormholes, as Bao et al. state
(their §2.2 for `N ≤ 4`, Lemma 4 in general). (iii) Corollary of Theorem 3* + Lemma 4: **for whole-boundary
parties the HRT entropy vector of any rotating vacuum wormhole is, up to a positive factor, the RT vector of
some static multiboundary wormhole with `N+1` boundaries** [proved mod Bao et al. Lemma 4] — rotation does
not enlarge the set of entropy rays.

---

## 8. What changes in the draft (recommendations; no existing file was edited)

1. PROOFS §5.4: replace "NOT proved" by Proposition 1 + Corollary 2 + Theorem 3* (or the R1 fill of §4).
   Theorem 3 and Theorem 3' then hold for arbitrary parties, in the strong form "`F^ε ∈ C_N` for small `ε`".
2. PROOFS §3, Lemma B: keep (it is used nowhere else once Theorem 3* is adopted), or replace by
   Proposition 1, which is exact and simpler; Lemma B's `κ = 1 − log 2` is recovered as a corollary.
3. Theorem 1/Theorem 2: can be presented as corollaries of Theorem 3* ("virtual static metric" form, §5).
   The contraction-map statement should in any case say "admits a contraction map after splitting
   repeated RHS terms", because one `C_5` facet has no weighted contraction map (§6; Avis–HC eq. (4.19)).
4. OPEN_GAPS G4: resolved for this class (no reliance on contraction-map completeness). G6a item 19:
   resolved. MANUSCRIPT §5/§10: "Theorem 3 for regions with endpoints: not proved" → proved.
5. Cite MG–T arXiv:2605.04031 (Corollary B, eq. (8.4); BIPP21 §4) for the current ⇒ pseudo-metric step and
   state that their Theorem A is for closed surfaces and is not used. Cite Bao et al. Lemma 4 for the cone
   equality, Avis–Hernández-Cuenca for the `C_5` facet list.
6. Do not claim more than is proved: the HRT statements inherit Theorem H's slice caveats.

---

## 9. Route R0 (contraction-map completeness shortcut), checked at the coordinator's request

### 9.1 What Bao–Furuya–Naskar actually state (arXiv:2506.18086 = JHEP 12 (2025) 140; full text read: §1, §2, §4)

* **Statement (their Theorem 4.1).** For a candidate inequality `Σ_i α_i S_{P_i} ≥ Σ_j β_j S_{Q_j}` with
  `α_i, β_j ∈ Z_+`, `P_i, Q_j` non-empty subsets of the `n` parties (purifier `O` separate), it is a valid
  HEI iff there is a contraction map `f : {0,1}^M → {0,1}^N` with `M = Σα_i`, `N = Σβ_j`, i.e. for the
  **unit-split form** (each term repeated with multiplicity, unweighted Hamming distances; this is "case
  (1b)" of Grimaldi–Headrick–Hubeny–Shteyner arXiv:2601.09987 App. A), with the boundary conditions
  `f(x_{A_k}) = y_{A_k}` for all `n+1` colours. "Valid" = holds for every RT geometry (equivalently every
  graph model, Bao et al.). Nothing restricts it to balanced inequalities: as stated it covers SA/Araki–Lieb
  and non-superbalanced facets too. Rational coefficients reduce to integers by scaling.
* **Their own disclaimers (§1):** (i) reliance on geometry ↔ graph duality; (ii) only linear inequalities
  with rational coefficients; (iii) "we take some of those insights, for instance, certain adjacency
  conditions, as physical requirements to construct our proof. We justify their existence and provide
  supporting evidence".
* **Why it is not a proof (my reading of §4).** Lemma 4.2 asserts: "for a non-contraction map `f̃`
  satisfying the homology conditions, there exists at least one geometry that does not satisfy the HEI".
  As stated this is false for every valid HEI that has a contraction map (MMI has many maps that are not
  contractions, yet no geometry violates MMI; BFN's own §3 calls such maps ones that "fail a valid HEI").
  The proof builds a graph `G̃` containing "virtual edges" `0̃_Ẽ` of **infinite** weight (eq. (58): the
  resulting `S*` "is divergent by definition", and "any geometry realized by `G̃` cannot give the proper
  holographic entanglement entropy"), and evaluates the RHS as the weight of the cuts **assembled through
  `f̃`** (eq. (60)), not as the RHS min-cuts. A violation of `Σ (assembled cut weights)` by an
  infinite-weight graph says nothing about the min-cut inequality on legitimate graphs. The step from
  "every non-contraction map is defective" to "some legitimate graph violates the inequality" — the whole
  content of completeness — is missing. **Grade of BFN Theorem 4.1: not established** (it may well be true;
  it is consistent with all known `N ≤ 6` data, but only in the split form).
* **Others' hesitation, confirmed.** GHHS 2601.09987 v1 §2.4: "an argument was recently given for the
  converse [BFN]. Throughout this paper, we will assume that this is the case"; v2 (comments field read):
  "clarified that our results apply to contraction-provable HEIs", i.e. they withdrew the assumption.
  GHHS App. A also records that the **weighted** form (case (2b), the one in PROOFS §1.5) fails for one
  `N = 5` facet (Avis–HC eq. (4.19)); our §6 search reproduces this. So "every valid HEI has a contraction
  map" is false in the weighted form and unproved in the split form.

### 9.2 What R0 would give, and its grade

(b) *Rescaling and splitting are harmless for Theorem 2* **[proved]**. Multiplying an inequality by a
positive integer multiplies both sides of (A1) and of Theorem 1's conclusion by the same factor. Splitting
`α_l S(A_l)` into `α_l` copies is allowed in Lemma A (the `A_l` need not be distinct: take the same minimiser
`γ_l` for repeated terms; the `δ_r` for repeated `B_r` may differ, and `F(B_r) ≤ Σ_j w_j ℓ_j(δ_r)` holds
for each copy); `Δ_p` is computed with multiplicity and is the same number. So Theorem 2 holds for every
inequality with a split-form contraction map, which covers the `C_5` facet with no weighted map.

(c) *Endpoint regulators.* Given (b), R0 would read: if BFN's Theorem 4.1 holds, every valid `Q` has a
split contraction map, so by Theorem 1, `Q(F^ε) ≥ 0` for all `ε ≤ ε_0` when `Q` is balanced at every marked
point, and `Q(F^ε) − W Δ_Q log(ε_0/ε) = const` with `Δ_Q = Σ_p Δ_p ≥ 0` otherwise (PROOFS Lemma A (A2)
gives `Δ_p ≥ 0` for every contraction-provable `Q`; for an arbitrary valid `Q`, `Δ_p(Q) = Q(b_{ij}) ≥ 0`
because the Bell pair is a graph). Since `C_N` has finitely many facets, `F^ε ∈ C_N` for `ε` below some
(non-explicit) threshold. **Grade: [conditional on BFN Theorem 4.1, which is not established].**

On Ju–Czech (arXiv:2611.21726: not every facet is superbalanced): no effect on either route. The per-endpoint
slack is `t_p(ε) Q(b_{ij})` with `(i,j)` the parties adjacent at `p`; `Q(b_{ij}) ≥ 0` for **every** valid `Q`
(Bell pairs are graph models), whether or not `Q` is superbalanced. A non-superbalanced facet with
`Q(b_{ij}) > 0` has divergent positive slack exactly when `i, j` are adjacent at some endpoint, and finite
`ε`-independent slack `Q(G)` otherwise. The `Δ ≥ 0` step of PROOFS Lemma A (A2)/Theorem 1(2) uses only the
contraction property at the two adjacent colours, so it too is valid for every contraction-provable facet,
balanced or not; nothing there assumes superbalance.

**Strongest correct statement for interval regions (this is Theorem 3*, unconditional):**
`F^ε = G + Σ_p t_p(ε) b_{i_p j_p}` exactly for all `ε ≤ ε_0`, with `G` the min-cut vector of an explicit
finite graph and `t_p(ε) = W log(1/ε) − Σ_j w_j log d_p^{(j)}`; hence (1) `F^ε ∈ C_N` for all
`ε ≤ ε_1 = min_p Π_j (d_p^{(j)})^{−w_j/W}` (explicit); (2) for every valid `Q` and every `ε ≤ ε_0`,
`Q(F^ε) = Q(G) + Σ_p t_p(ε) Q(b_{i_p j_p})`, both terms `≥ 0` once `t_p ≥ 0`, the first independent of `ε`;
(3) scheme independence: changing the horoball normalisation at `p` shifts `t_p` by a constant and `G` by a
multiple of `b_{i_p j_p}`, so "`F^ε ∈ C_N` for all sufficiently small `ε`" and the `ε`-independent slack of
balanced inequalities are scheme-independent, while the statement "`G ∈ C_N`" holds in the adapted scheme.

### 9.3 Decision

R0 does not work as a rigorous route: its key input is an unproved claim (and false in the weighted form
the draft uses). **R2/Theorem 3* remains the main route**; it needs no contraction maps, gives an explicit
finite graph and an explicit cut-off threshold, and is strictly stronger than what R0 would give. R0 is
recorded only as a conditional consistency check: both routes predict the same slack structure.

### 9.4 (d) Reverse inclusion `C_N ⊆ cl(R_+ 𝓗_N)`: literature check

* Bao et al. arXiv:1505.07839 (read §§2.2, 3.2, 3.3, App. A): **Lemma 4** (graph → hyperbolic surface;
  pairs of pants per trivalent bulk vertex, half-collars per boundary vertex, gluing geodesics of length
  `4G_N w(e)`, curvature sufficiently negative; proof by the collar lemma) and **Theorem 5** ("(1) can be
  restricted to two-dimensional manifolds of constant negative curvature"). The surfaces are `t = 0` slices
  of static multiboundary wormholes (their §5; Skenderis–van Rees, Brill). Parties are **whole boundaries,
  one boundary per colour**. Exact realisation, up to the overall scale fixed by the curvature/`G_N`.
  **No higher dimensions are needed** in 2+1; §3.3 only remarks that the construction generalises.
* Every extreme ray, indeed every graph vector, is therefore realised up to scale by a static vacuum
  2+1 spacetime with whole-boundary parties; with our §7.2 this gives `cl(R_+ 𝓗_N) = C_N` **[proved mod
  Bao et al. Lemma 4, §3.1 moves, Thm 5/Prop. 7]**.
* **Not claimed / not checked:** realisation of all extreme rays with *interval* parties on a fixed
  topology (e.g. intervals on the single boundary of pure AdS3 or BTZ). Bao et al. App. A only sketches a
  discretisation ("we believe that the method can … be made mathematically rigorous"), which goes
  geometry → graph, not graph → interval geometry. The cone equality does not need it, because whole
  boundaries are among the allowed parties.
