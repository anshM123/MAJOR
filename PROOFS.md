# PROOFS (revision 5: Hypothesis H removed)

Status legend: **[proved]** = complete argument written here, relying only on cited standard facts;
**[standard]** = textbook fact, cited; **[assumed]** = hypothesis, stated explicitly where used.

§§1–5 are mathematics about hyperbolic surfaces and contain no physics. Gravity enters only in §6. Since
revision 5 the link to the HRT formula is a theorem (Theorem H, §6.7): the HRT entropy is *defined* in §6.3 by
the usual prescription (minimal extremal curve subject to the spacelike homology condition) and *proved* to
equal the minimum over topological cut classes.

**Attribution.** The mechanism — perform the cut-and-paste of the contraction-map proof in an auxiliary
geometry, keep only crossing counts (which are topological), and weight them by a positive measure on
geodesics — is due to Czech and Dong (arXiv:1905.03787, §§4, 7), who used it for simply connected 2+1 bulks.
Their extension to multiply connected bulks via the universal cover is incorrect (Grado-White, Grimaldi,
Headrick, Hubeny, arXiv:2407.07165: "an infinite number of surfaces enter on both sides of the inequality").
What is done here: the argument is run on the quotient surface, with test classes for all homotopy classes
and the measure on `𝒢/Γ` (Lemma A, Lemma B), and positivity is supplied by the two chiral Liouville measures
of a vacuum solution. This repairs the multiply connected step for vacuum solutions only.

---

## 1. Setting

### 1.1 Surfaces, marked ideal points, regions

* `Γ ⊂ PSL(2,R)` is finitely generated, torsion-free, with every non-trivial element hyperbolic, and of the
  second kind (limit set `Λ ≠ S¹`). `X = H²/Γ` is a complete orientable hyperbolic surface of finite topology
  all of whose ends are funnels. The trivial group (`X = H²`) and cyclic groups (hyperbolic cylinder) are
  included. `Ω = S¹ \ Λ`, and `∂_∞X = Ω/Γ` is a finite union of circles, one per funnel (for trivial `Γ`, the
  whole circle `S¹`). `X̄ = X ∪ ∂_∞X`.
  **[standard]** `Γ` is convex cocompact, the Hausdorff dimension of `Λ` is `< 1`, and in particular `Λ` has
  Lebesgue measure zero (Patterson–Sullivan theory; see e.g. Patterson, *Acta Math.* 136 (1976), and Nicholls, *The Ergodic Theory of
  Discrete Groups*. These citations were not re-checked against the sources.)
* `P ⊂ ∂_∞X` is a finite set of *marked ideal points*; the components of `∂_∞X \ P` are the *elementary
  pieces*.
* A *party assignment* maps each elementary piece to one of `1,…,N,O` (`O` = purifier). A *region*
  `A ⊂ {1,…,N}` is identified with the closure of the union of the pieces of parties in `A`. `∂A ⊂ P` is the
  set of marked points at which exactly one of the two adjacent pieces lies in `A`.

### 1.2 Cut classes

A *cut* for `A` is a finite disjoint union `c` of smoothly embedded curves in `X`, each a closed curve or a
proper arc whose two ends converge to marked points, such that

1. exactly one arc-end converges to each `p ∈ ∂A` and none to any other marked point;
2. there is a closed set `W ⊂ X`, equal to the closure of its interior, whose topological boundary in `X` is
   `c`, and whose closure in `X̄` meets `∂_∞X` exactly in `A` (the *homology region*).

`𝒞(A)` is the set of cuts for `A` modulo isotopies of `X̄` fixing `P`, deletion of null-homotopic closed
components, and deletion of pairs of isotopic closed components. Both deletions preserve condition 2 (swap
the colour of the open disc, resp. open annulus, removed; its closure contains no point of `∂_∞X`). Every
class has a representative without null-homotopic components and without repeated components.
`e(c) := |∂A|` is the number of arc-ends, the same for every `c ∈ 𝒞(A)`.

### 1.3 Geodesic representatives **[standard]**

(a) Every non-null-homotopic closed curve is freely homotopic to a unique closed geodesic; every proper arc
between marked points `p ≠ p'` is homotopic, rel. its ideal endpoints, to a unique complete geodesic.

(b) Two geodesics of `H²` meet iff their endpoint pairs strictly link in `S¹`, and then once, transversally.

(c) For curves `σ, σ'`, each a closed geodesic or a complete geodesic with ideal endpoints, with `σ ≠ σ'`,
the number of intersection points counted with multiplicity (pairs of parameters) is the minimum over all
pairs of transverse curves homotopic to them rel. ideal endpoints. (Lift: by (b) the count equals the number
of `Γ`-orbits of linked pairs of lifts, and linked lifts of any representatives must cross.)

(d) A cut without null-homotopic or repeated components is isotopic (fixing `P`) to the union of the
geodesic representatives of its components, which is again a cut for `A`. (Epstein, *Acta Math.* 115 (1966);
Farb–Margalit, *A Primer on Mapping Class Groups*, §1.2; for ideal arcs the same argument applies in `X̄`.
Two geodesic arcs ending at the same ideal point are disjoint by (b).)

### 1.4 Test classes and intersection numbers

A *test class* `τ` is the free homotopy class of a closed curve in `X`, or the proper homotopy class of a
proper arc whose ends converge to points of `∂_∞X \ P`, homotopies being allowed to move the ends within
`∂_∞X \ P`. For a test class `τ` and a multi-curve class `c` as in §1.2,

`i(τ, c) := min #(τ' ∩ c')`,

the minimum over transverse representatives, intersections counted with multiplicity; `i(τ,c) = 0` if `τ` is
null-homotopic. If `τ^X` is a geodesic representative of `τ` (any ideal endpoints in the prescribed elementary
pieces) and `c^X` the geodesic representative of `c`, then

`i(τ, c) = #(τ^X ∩ c^X)` (if a closed component `κ` of `c^X` equals `τ^X` or a geodesic that `τ^X` covers, replace
`τ^X` by a parallel copy before counting).   (1.1)

Indeed lifts of components of `c^X` have endpoints in `Λ ∪ Γ·P̃`; moving an ideal endpoint of a lift of `τ^X`
inside a component of `Ω \ Γ·P̃` changes no linking, and (c) applies. For a closed geodesic `κ` that is simple,
`i(κ,κ) = 0`.

### 1.5 Contraction maps

`Σ_{l=1}^{L} α_l S(A_l) ≥ Σ_{r=1}^{R} β_r S(B_r)` (`α_l, β_r > 0`) *admits a contraction map* if there is
`f: {0,1}^L → {0,1}^R` with `‖f(x) − f(x')‖_β ≤ ‖x − x'‖_α` for all `x, x'` (weighted Hamming norms) and
`f(x^{(i)}) = y^{(i)}` for all `i ∈ {1,…,N,O}`, where `x^{(i)}_l = [i ∈ A_l]`, `y^{(i)}_r = [i ∈ B_r]`
(Bao–Nezami–Ooguri–Stoica–Sully–Walter, arXiv:1505.07839).

---

## 2. Lemma A (topological contraction lemma) **[proved]**

**Lemma A.** Let `(X, P)` carry a party assignment, let `Σ α_l S(A_l) ≥ Σ β_r S(B_r)` admit a contraction map
`f`, and let `γ_l ∈ 𝒞(A_l)` be arbitrary. Then there are `δ_r ∈ 𝒞(B_r)` such that for every test class `τ`

`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)`,   (A1)

and for every marked point `p`

`Δ_p := Σ_l α_l [p ∈ ∂A_l] − Σ_r β_r [p ∈ ∂B_r] ≥ 0`.   (A2)

The `δ_r` depend on the `γ_l`, on `f`, and on the hyperbolic structure `X` used in the construction, not on `τ`.

*Proof.* (A2). Let `a, b` be the pieces adjacent to `p`, of parties `i_a, i_b`. Then
`[p ∈ ∂A_l] = |x^{(i_a)}_l − x^{(i_b)}_l|`, `[p ∈ ∂B_r] = |y^{(i_a)}_r − y^{(i_b)}_r|`; apply the contraction
property to `x^{(i_a)}, x^{(i_b)}`.

*Construction.* Replace each `γ_l` by its geodesic representative (§1.3(d)), with homology region `W_l`. Let
`𝒜` be the union of all these geodesics, a finite union of closed geodesics and complete geodesics. A
*vertex* is a point where two distinct geodesics of `𝒜` cross; there are finitely many. For a geodesic
`κ ⊂ 𝒜` let `Λ(κ) = {l : κ ⊂ γ_l}`.

The components of `X \ 𝒜` are *cells*. For a cell `C` put `x(C)_l = [C ⊂ W_l]`; this is well defined because
`W_l` is the closure of a union of components of `X \ γ_l`. If `l ∈ Λ(κ)`, the two sides of `κ` near a
non-vertex point lie one in `W_l` and one outside (as `W_l` is the closure of its interior with boundary
`γ_l`); if `l ∉ Λ(κ)` both lie on the same side. So the cells on the two sides of `κ` there are distinct and
their bit strings differ exactly on `Λ(κ)`. A cell whose closure in `X̄` contains an elementary piece of party
`i` has bit string `x^{(i)}`.

Let `U_r` be the closure of the union of cells `C` with `f(x(C))_r = 1`. By the boundary conditions its
closure in `X̄` contains each piece of a party in `B_r` and no other piece. `∂U_r` is the union of the closed
edges of `𝒜` across which `f(x(·))_r` changes.

*Surgery.* Near a vertex `v`, the sectors around `v` alternate in membership of `U_r` across the half-edges
of `∂U_r` at `v`, so there is an even number of such half-edges; inside a small disc around `v` replace `U_r`
by a set whose boundary consists of disjoint arcs joining these half-edges in cyclically adjacent pairs. Near
a marked point `p`, the geodesics of `𝒜` ending at `p` are pairwise disjoint or equal; they divide a
neighbourhood of `p` into sectors, and those belonging to `∂U_r` are rays `ρ_1,…,ρ_k` in angular order, with
`k ≡ [p ∈ ∂B_r] (mod 2)` since the outermost sectors contain the two pieces adjacent to `p`. Inside a horoball
`H_p` at `p` modify `U_r` by joining `ρ_1` to `ρ_2`, `ρ_3` to `ρ_4`, … along `∂H_p` and discarding the parts of
those rays inside `H_p`; a ray remains iff `p ∈ ∂B_r`. Call the result `U'_r`. It is the closure of its
interior, its closure in `X̄` meets `∂_∞X` exactly in `B_r`, and its boundary is an embedded multi-curve with
exactly one arc-end at each point of `∂B_r` and none elsewhere. So `∂U'_r` is a cut for `B_r`; let `δ_r` be
its class. The class does not depend on the sizes of the discs and horoballs.

*Proof of (A1).* If `τ` is null-homotopic both sides vanish. Otherwise let `τ^X` be a geodesic
representative. If `τ^X` is a closed geodesic contained in `𝒜`, replace it by a parallel copy on one side.
Perturb near vertices so as to miss them. The resulting `τ'` represents `τ`, is transverse to `𝒜`, misses the
vertices, and satisfies `#(τ' ∩ κ) = i(τ, κ)` for every geodesic `κ ⊂ 𝒜`: by (1.1) when `κ ≠ τ^X` (the
perturbations slide crossings along `κ` without creating or destroying any), and because both sides are `0`
when `κ = τ^X` (a component of a cut is simple). As `τ'` is a fixed curve whose ends do not converge to
marked points, it avoids all sufficiently small `H_p` and vertex discs, so `τ' ∩ ∂U'_r = τ' ∩ ∂U_r`.

Each `q ∈ τ' ∩ 𝒜` (with multiplicity) lies on a unique geodesic `κ_q` and separates cells `C_q, C'_q` with bit
strings differing exactly on `Λ(κ_q)`. Therefore

`Σ_l α_l i(τ,γ_l) = Σ_l α_l Σ_{κ ⊂ γ_l} #(τ' ∩ κ) = Σ_q Σ_{l ∈ Λ(κ_q)} α_l = Σ_q ‖x(C_q) − x(C'_q)‖_α`
`≥ Σ_q ‖f(x(C_q)) − f(x(C'_q))‖_β = Σ_r β_r #(τ' ∩ ∂U_r) = Σ_r β_r #(τ' ∩ ∂U'_r) ≥ Σ_r β_r i(τ, δ_r)`.

The last inequality holds because `τ'` and `∂U'_r` are transverse representatives of `τ` and of a cut whose
class, after deleting null-homotopic and repeated components, is `δ_r`. ∎

---

## 3. Lemma B (Crofton formula with horoball cut-offs) **[proved]**

### 3.1 Normalisations

`𝒢` = unoriented geodesics of `H²`; `m` = the invariant measure with

`m{g : g ∩ s ≠ ∅} = length(s)` for every geodesic segment `s`   **[standard]** (Santaló, *Integral Geometry
and Geometric Probability*, Part III).

In the upper half-plane, with Euclidean centre `c` and radius `ρ`, `dm = dc dρ/(2ρ²)` (equivalently
`da db/(b−a)²` in the endpoints `a < b`). Check: geodesics meeting `{x = 0, y_1 < y < y_2}` are
`{|c| < ρ, y_1² < ρ² − c² < y_2²}`; with `c = ρ sin φ` the measure is `½∫cos φ dφ ∫dρ/ρ = log(y_2/y_1)`.

`m̄` = induced measure on `𝒢/Γ`. For a compact geodesic segment `s ⊂ X`, `∫ #(ĝ ∩ s) dm̄ = length(s)`.

Since `Λ` is null and `Γ·P̃` countable, `m̄`-a.e. `ĝ` has both endpoints in `Ω \ Γ·P̃`; such `ĝ` is a proper arc
with ends in `∂_∞X \ P`, hence defines a test class `τ(ĝ)`, and by (1.1)

`#(ĝ ∩ c^X) = i(τ(ĝ), c)` for `m̄`-a.e. `ĝ`, for every class `c`.   (B0)

### 3.2 Choice of cut-offs

For each `p ∈ P` fix a lift `p̃ ∈ Ω` and an upper half-plane coordinate with `p̃ = ∞`. The set
`Λ ∪ (Γ·P̃ \ {∞})` is bounded, say contained in `[−Q_p, Q_p]` (its accumulation points lie in `Λ`, and `∞` is an
interior point of `Ω`). Put `B_p̃(ε) = {y > 1/ε}`, extended `Γ`-equivariantly. Fix `ε_0` with
`1/ε_0 > max_p Q_p` and so small that the horoballs `B_p(ε_0) ⊂ X` are embedded and pairwise disjoint.

**Fact 3.1.** A geodesic of `H²` with both endpoints in `Λ ∪ Γ·P̃` and not ending at `∞` has Euclidean radius
`≤ Q_p`, hence is disjoint from `B_p̃(ε_0)`. Consequently, for every class `c` and every `ε ≤ ε_0`, `c^X` meets
`⋃_p B_p(ε)` exactly in the terminal rays of its arc-ends, and

`ℓ^ε(c) := length(c^X \ ⋃_p B_p(ε)) = ℓ^{ε_0}(c) + e(c)·log(ε_0/ε)`.   (3.1)

### 3.3 The lemma

Let `𝒢_ε ⊂ 𝒢/Γ` be the set of geodesics of `X` disjoint from `⋃_p B_p(ε)`.

**Fact 3.2.** If `g ∈ 𝒢` has both endpoints outside `Γ·P̃`, then `g` misses all horoballs of the family at
scale `ε` for `ε` small. (The endpoint pairs of the translates `γg` accumulate only on `Λ`, so only finitely
many translates have an endpoint in a given compact subset of the component of `Ω` containing `p̃`, and
for each of these the Euclidean radius is finite; hence `Γg` misses `{y > 1/ε}` for small `ε`. Do this for
each `p`.) Thus the sets `(𝒢/Γ) \ 𝒢_ε` decrease, as `ε ↓ 0`, to an `m̄`-null set.

**Lemma B.** For every class `c`, as `ε → 0`,

`ℓ^ε(c) = ∫_{𝒢_ε} #(ĝ ∩ c^X) dm̄(ĝ) + κ·e(c) + o(1)`,  `κ = 1 − log 2`.

Moreover `Def_ε(c) := ℓ^ε(c) − ∫_{𝒢_ε} #(ĝ ∩ c^X) dm̄ ≥ 0` for all `ε ≤ ε_0`.

*Proof.* For `ĝ ∈ 𝒢_ε`, `ĝ ∩ c^X = ĝ ∩ c^ε` where `c^ε = c^X \ ⋃B_p(ε)`. By Crofton,
`Def_ε = ∫_{(𝒢/Γ)\𝒢_ε} #(ĝ ∩ c^ε) dm̄ ≥ 0`. Split `c^ε = K ∪ ⋃_{ends} σ_ε`, `K = c^{ε_0}` compact, and for each
arc-end (at `p`) `σ_ε` the segment of the terminal ray between `∂B_p(ε_0)` and `∂B_p(ε)`.

(i) *Core.* `∫ #(ĝ ∩ K) dm̄ = length(K) < ∞`; by Fact 3.2 and dominated convergence
`∫_{(𝒢/Γ)\𝒢_ε} #(ĝ ∩ K) dm̄ → 0`.

(ii) *One end.* Lift so that the terminal ray is `{x = x_0}` going to `p̃ = ∞`, `|x_0| ≤ Q_p`; translate to
`x_0 = 0` (then the other relevant boundary points lie in `[−2Q_p, 2Q_p]`). `σ̃_ε = {x = 0, Y_0 < y < Y}`,
`Y_0 = 1/ε_0`, `Y = 1/ε`. Unfolding, the contribution of this end is `m(E_1) + m(E_2)`,

`E_1 = {g : g ∩ σ̃_ε ≠ ∅, g ∩ {y > Y} ≠ ∅}`,
`E_2 = {g : g ∩ σ̃_ε ≠ ∅, g ∩ {y > Y} = ∅, g meets another horoball of the family at scale ε}`.

`g` meets `{y > Y}` iff `ρ > Y`, and meets `σ̃_ε` iff `|c| < ρ`, `Y_0 < ρ cos φ < Y` (`c = ρ sin φ`). Hence

`m(E_1) = ½∫_{−π/2}^{π/2} cos φ dφ ∫[ρ > Y, Y_0 < ρ cos φ < Y] dρ/ρ = ∫_0^{π/2} cos φ · log(1/max(cos φ, Y_0/Y)) dφ`,

which increases to `∫_0^{π/2} cos φ log(1/cos φ) dφ = 1 − log 2` as `Y → ∞`.

For `E_2`: the other horoballs are based in `[−2Q_p, 2Q_p]` and, being disjoint from `{y > Y_0}`, have Euclidean
diameter `≤ Y_0` at scale `ε_0`; they lie in the box `R = {|x| ≤ 2Q_p + Y_0/2, 0 < y ≤ Y_0}`. Let `E_2'(ε) ⊃ E_2`
be the set of geodesics meeting the ray `{x = 0, y > Y_0}` and some other horoball at scale `ε`. A geodesic
through `(0,h)`, `h > Y_0`, has endpoints `a < 0 < b` with `|a| b = h²`. If both endpoints have modulus
`> Q' := 2Q_p + 3Y_0/2`, the semicircle has radius `> Q'` and its height above every `x` with
`|x| ≤ 2Q_p + Y_0/2` exceeds `Y_0`, so it misses `R`. Hence a member of `E_2'` has an endpoint of modulus
`≤ Q'`; by symmetry take `|a| ≤ Q'`; then
`b > Y_0²/|a|`, and

`m(E_2'(ε_0)) ≤ 2∫_0^{Q'} d|a| ∫_{Y_0²/|a|}^{∞} db/(b + |a|)² = 2∫_0^{Q'} |a| d|a| /(Y_0² + a²) < ∞`.

The sets `E_2'(ε)` decrease as `ε ↓ 0` to a set of geodesics with an endpoint in `Γ·P̃` (Fact 3.2), which is
null. So `m(E_2) → 0`. (All these sets, for every terminal ray `x_0 ∈ [−Q_p, Q_p]`, are contained in the single
finite-measure set `{a ∈ [−Q', Q'] , |b − a| ≥ 2Y_0}` up to translation by `x_0`; intersected with "meets another horoball at scale `ε`" this is one
`x_0`-independent family of sets decreasing to a null set, so the convergence is uniform in the ray.)

Summing over the `e(c)` ends, `Def_ε → κ e(c)`. ∎

---

## 4. Theorem 1 (multi-metric Ryu–Takayanagi inequalities) **[proved]**

**Setting.** A reference surface `(X_1, P_1)` with party assignment, and for `j = 2,…,k` surfaces `(X_j, P_j)`
as in §1 with homeomorphisms `φ_j : X̄_1 → X̄_j`, `φ_j(P_1) = P_j`, carrying the party assignment. Cut classes
and test classes on `X_j` are identified with those on `X_1` via `φ_j`; intersection numbers are defined by a
minimum over representatives (§1.4) and are therefore preserved. Fix cut-offs as in §3.2 on each `X_j` with a
common `ε_0`, and weights `w_j > 0`. For a region `A`,

`F^ε(A) := min_{c ∈ 𝒞(A)} Σ_j w_j ℓ^ε_{X_j}(c)`,  `ε ≤ ε_0`.

By (3.1), `Σ_j w_j ℓ^ε_{X_j}(c) = Σ_j w_j ℓ^{ε_0}_{X_j}(c) + (Σ_j w_j)|∂A| log(ε_0/ε)` for every `c ∈ 𝒞(A)`. For
each bound `M` only finitely many classes have `ℓ^{ε_0}_{X_1}(c) ≤ M` (their geodesic representatives lift to
finitely many `Γ`-orbits of geodesics meeting a compact set with bounded length, by discreteness). So the
minimum exists, the set of minimisers is finite and independent of `ε ≤ ε_0`, and

`F^ε(A) = F^{ε_0}(A) + (Σ_j w_j)|∂A| log(ε_0/ε)`.   (4.1)

**Theorem 1.** Let `Σ_l α_l S(A_l) ≥ Σ_r β_r S(B_r)` admit a contraction map, and let
`Δ = Σ_p Δ_p = Σ_l α_l|∂A_l| − Σ_r β_r|∂B_r| ≥ 0` (Lemma A (A2)).

1. If `Δ = 0` (the inequality is balanced at every marked point; in particular if `P = ∅`), then for every
   `ε ≤ ε_0`: `Σ_l α_l F^ε(A_l) ≥ Σ_r β_r F^ε(B_r)`.
2. If `Δ > 0`, then `Σ_l α_l F^ε(A_l) − Σ_r β_r F^ε(B_r) = C + (Σ_j w_j) Δ log(ε_0/ε) → +∞` as `ε → 0`, with `C`
   independent of `ε`.

*Proof.* Fix minimisers `γ_l` of `F^ε(A_l)` (independent of `ε`). Lemma A, applied on `X_1`, gives
`δ_r ∈ 𝒞(B_r)` with (A1), independent of `ε`. For each `j`, integrate (A1) with `τ = τ(ĝ)` over
`ĝ ∈ 𝒢_ε(X_j)` against `m̄_j`, using (B0) on `X_j` (the class of an `X_j`-geodesic with ends off `P_j`
corresponds under `φ_j^{-1}` to a test class on `X_1`):

`Σ_r β_r ∫_{𝒢_ε(X_j)} #(ĝ ∩ δ_r^{X_j}) dm̄_j ≤ Σ_l α_l ∫_{𝒢_ε(X_j)} #(ĝ ∩ γ_l^{X_j}) dm̄_j`,

all integrals finite. By Lemma B, for the finitely many classes involved,

`Σ_r β_r ℓ^ε_{X_j}(δ_r) ≤ Σ_l α_l ℓ^ε_{X_j}(γ_l) − κΔ + o(1)`.

Multiply by `w_j`, sum, and use `F^ε(B_r) ≤ Σ_j w_j ℓ^ε_{X_j}(δ_r)`:

`Σ_l α_l F^ε(A_l) − Σ_r β_r F^ε(B_r) ≥ κ(Σ_j w_j)Δ − o(1)`  (`ε → 0`).   (4.2)

By (4.1) the left side equals `C + (Σ_j w_j)Δ log(ε_0/ε)`. If `Δ = 0` it is the constant `C`, and (4.2) gives
`C ≥ −o(1)`, i.e. `C ≥ 0`. If `Δ > 0` statement 2 is (4.1). ∎

**Remarks.** (0) In the balanced case the proof gives more: `Σ_r β_r ℓ^{ε_0}_{X_j}(δ_r) ≤ Σ_l α_l ℓ^{ε_0}_{X_j}(γ_l)` for
each `j` separately. (1) For `k = 1` this is the Ryu–Takayanagi statement on a hyperbolic surface. (2) Only two
properties of length are used: it is the integral of intersection numbers against a positive measure on test
classes, and its divergence at each end is universal.

---

## 5. Theorem 3 (membership in the full Ryu–Takayanagi cone) **[proved for `P = ∅`; NOT proved for `P ≠ ∅`]**

This section removes the restriction to inequalities with contraction maps, for parties made of whole
boundary circles. The observation that the functional is a graph min-cut for finite families of test curves
is due to the referee, who also verified Lemma G numerically on the disc (5464/5464 cases).

### 5.1 Lemma G (graph realisation) **[proved]**

Let `μ = Σ_{k=1}^{m} w_k τ_k` be a finite positive combination of test classes, none null-homotopic and none
homotopic into `∂_∞X`. Without loss of generality the `τ_k` are distinct (merge equal classes, adding
weights) and primitive (replace `σ^n` by `σ` with weight `n w`; `i(σ^n, c) = n·i(σ, c)`). Put
`F_μ(A) := min_{c ∈ 𝒞(A)} i(μ, c)`, `i(μ,c) = Σ_k w_k i(τ_k, c)`.

**Lemma G.** There is a finite weighted graph `G`, with some vertices labelled by parties, such that for every
region `A`, `F_μ(A)` is the minimum weight of an edge cut of `G` separating the vertices labelled by parties
in `A` from those labelled by parties not in `A` (purifier included). Hence `(F_μ(A))_A` is a graph-model
entropy vector in the sense of Bao et al. (several terminals with the same label may be merged through an
auxiliary vertex joined to them by edges of very large weight; a party with no terminal is an isolated
terminal), and lies in the holographic entropy cone.

*Proof.* Represent the `τ_k` by geodesics `τ_k^X` (closed, or with ideal endpoints in the prescribed pieces,
chosen pairwise distinct), and let `T` be their union. No perturbation is made. The vertices of `G` are: all
points of `T` through which more than one branch passes (any multiplicity); the ideal ends of the arcs, each
labelled by the party of the piece containing it; and one extra vertex on each closed `τ_k^X` carrying no
other vertex. The edges are the segments of `T` between consecutive vertices; an edge on `τ_k^X` has weight
`w_k`.

(≥) Let `c ∈ 𝒞(A)` with geodesic representative `c^X` and homology region `W`. If a component of `c^X` coincides
with some closed `τ_k^X`, push it to a parallel copy (its count with `τ_k` becomes `0 = i(τ_k, ·)` for that
component, counts with the other branches are unchanged). Slide `c^X` off the vertices of `T`, keeping one
crossing with each branch through the vertex. Then `#(c ∩ τ_k^X) = i(τ_k, c)` for all `k`. Let `S` be the set
of vertices in `W`; an ideal end lies in `S` iff its piece lies in `A`. An edge with exactly one endpoint in
`S` is a segment from a point of `W` to a point outside `W`, so it meets `c`. Hence
`weight(S, S^c) ≤ i(μ, c)`.

(≤) Let `S` be a vertex set containing the ideal ends lying in pieces of `A` and no other ideal end. Let `W` be
the union of a small closed regular neighbourhood of the vertices of `S` and of the edges with both endpoints
in `S`, and a thin closed collar in `X̄` of the closed set `A ⊂ ∂_∞X`, so thin that it meets `T` only near ideal
ends lying in `A` (which are in `S`). `W` is the closure of its interior, its closure meets `∂_∞X` exactly in
`A`, and `∂W` is an embedded multi-curve whose arcs end exactly at the points of `∂A`: a cut for `A`. It crosses
each edge between `S` and `S^c` once and no other edge, so `F_μ(A) ≤ weight(S, S^c)`. ∎

### 5.2 Lemma F (filling) **[proved, using the standard theory of normal coordinates]**

**Lemma F.** There is a finite set `T_0` of test classes of proper arcs, none homotopic into `∂_∞X`, each
containing an open set of geodesics, such that for every `M` the set
`{c ∈ ⋃_A 𝒞(A) : Σ_{τ ∈ T_0} i(τ, c) ≤ M}` is finite.

*Proof.* If `X̄` is a disc, `T_0 = ∅` and the claim is the finiteness of non-crossing pairings of `P`. Otherwise
choose a system of pairwise disjoint simple essential geodesic arcs with ideal endpoints in `∂_∞X \ P`
cutting `X̄` into closed discs (one arc for an annulus; `2g + b − 1` in general); every boundary circle then
carries an arc end. Moving the ideal endpoints slightly keeps the class, so each class contains an open set
of geodesics. A cut in minimal position with respect to the system meets each disc in arcs joining sides of
the disc or marked points on its ideal sides, without bigons; by the standard normal-coordinate argument its
class is determined by the numbers of arcs of each of the finitely many combinatorial types in each disc,
and these numbers are bounded in terms of the intersection numbers with the system. ∎

### 5.3 The case `P = ∅` (parties are unions of whole boundary circles) **[proved]**

**Theorem 3.** With the notation of Theorem 1 and `P = ∅`, the vector `(F(A))_A`,
`F(A) = min_c Σ_j w_j ℓ_{X_j}(c)`, lies in the holographic entropy cone: it satisfies every linear inequality
valid for Ryu–Takayanagi entropies (equivalently, for all graph models).

*Proof.* By Crofton and (B0), `ℓ_{X_j}(c) = ∫ i(τ(ĝ), c) dm̄_j(ĝ)`. Geodesics whose class is homotopic into
`∂_∞X` (both ends in one boundary circle, cutting off a half-plane) have `i ≡ 0` with every cut; discard them
(their total measure is infinite). Each remaining test class `τ` has finite `m̄_j`-measure `μ_τ^{(j)}`: let `β` be the
boundary geodesic of a funnel containing an end of `τ`. A lift of an arc not homotopic into `∂_∞X` has its
endpoints in two different components of `Ω`, so it crosses the axis bounding each of them; hence
`i(τ, β) ≥ 1`. Every geodesic of class `τ` meets `β` exactly `i(τ, β)` times, so
`μ_τ^{(j)} · i(τ, β) ≤ ℓ_{X_j}(β) < ∞` by Crofton. Test classes form a countable set, so the push-forward of `Σ_j w_j m̄_j`
(transported to `X_1`) is a countable sum `μ = Σ_τ μ_τ τ` with finite coefficients and
`i(μ, c) = Σ_j w_j ℓ_{X_j}(c) < ∞` for every `c`. Let `μ_n` be finite partial sums by whole classes, increasing
to `μ`, containing `T_0` from the start. By Lemma G, `(F_{μ_n}(A))_A` lies in the cone. `F_{μ_n}(A)` is
non-decreasing in `n` and `≤ F(A)`. A minimiser `c_n` of `F_{μ_n}(A)` (it exists: by Lemma F only finitely many
classes have `i(μ_n, ·) ≤ F(A)`) satisfies `Σ_{τ∈T_0} μ_τ i(τ,c_n) ≤ F(A)`, so the `c_n` range over a finite
set `𝒟`; then `lim_n F_{μ_n}(A) = lim_n min_{c∈𝒟} i(μ_n,c) = min_{c ∈ 𝒟} i(μ,c) ≥ F(A)`. Hence
`F_{μ_n}(A) → F(A)` for every `A`. The cone is closed (it is polyhedral, Bao et al.). ∎

### 5.4 The case `P ≠ ∅` **[NOT proved]**

**Expected statement.** For every inequality `Q` valid on the holographic entropy cone: if `Q` is balanced at
every marked point then `Q(F^ε) ≥ 0` for all `ε ≤ ε_0`; otherwise `Q(F^ε) → +∞`.

`Δ_p ≥ 0` for such `Q` is proved (apply `Q` to the graph with one unit edge joining the two parties adjacent
at `p`). The remainder is a sketch: push `Σ_j w_j m̄_j|_{𝒢_ε}` forward to a countable measure `ν_ε`, apply Lemma
G to finite partial sums, and control the minimisers by splitting `ν_ε` into the part from `𝒢_{ε_0}` and the
rest. The second referee judged the skeleton sound and completable, and listed what is missing:

1. discard the infinite-mass trivial classes and show each remaining class has finite `ν_ε`-mass;
2. existence of the minimisers `δ_r^n`;
3. a two-stage argument: first, for fixed `ε`, the `δ_r^n` lie in a finite set (Lemma F and
   `i(μ_n, δ_r^n) ≤ i(ν_ε, c_ref)`), and only then is the tail `i(ν_ε − μ_n, δ_r^n)` small for `n ≥ n(ε)`;
4. uniformity in the class of the end estimate (available from the uniform domination noted in Lemma B);
5. `T_0` must have positive mass in `𝒢_{ε_0}`: choose `T_0` first, then shrink `ε_0` (Fact 3.2);
6. the explicit bound replacing "O(1)".

Until these are written out, **Theorem 3 is claimed only for `P = ∅`.**

---

## 6. Application to vacuum 2+1 gravity

### 6.1 The spacetimes **[standard]**

`AdS_3` = `PSL(2,R)` with `ds² = −det(dx)` (computations are done with matrix representatives in `SL(2,R)`,
defined up to sign); `PSL(2,R) × PSL(2,R)` acts by `x ↦ g_L x g_R^{-1}`.

Let `S` be an oriented surface of finite type with non-empty boundary and negative Euler characteristic, or a
disc or an annulus, and `ρ_L, ρ_R : π_1(S) → SL(2,R)` two representations whose images are groups as in §1.1
and which are holonomies of hyperbolic structures `X_L, X_R` on the interior of `S` **in the same component
of Teichmüller space** (i.e. `X_R` is obtained from `X_L` by an orientation-preserving homeomorphism isotopic
to the identity relative to the marking; boundary elements go to boundary elements). Then
`Γ = {(ρ_L(γ), ρ_R(γ))}` acts properly discontinuously on the *invisible domain* `D̃ ⊂ AdS_3` of its limit set
(the open set of points not causally related to any point of the limit set `Λ_Γ ⊂ ∂AdS_3`, which is
the graph of the boundary map over `Λ_L`; a Cantor set for these groups); `M = D̃/Γ` is a
locally-`AdS_3` spacetime with one asymptotically-AdS region per end of `S` (Mess 1990; Barbot, *Causal
properties of AdS-isometry groups I, II*; Bonsante–Krasnov–Schlenker, arXiv:math/0610429). Physically:
rotating BTZ (`π_1 = Z`), rotating multi-boundary wormholes, wormholes with handles
(Aminneborg–Bengtsson–Holst 1998). `ρ_L = ρ_R` is the time-symmetric case; for `ρ_L ≠ ρ_R` the spacetime has
no time-reflection symmetry, and for non-abelian `π_1` it is not stationary.

**Scope.** Only these quotients are covered: orientable, all holonomies hyperbolic (no extremal horizons,
no cusps, no conical defects).

The `i`-th asymptotic boundary is `(I^i_L × I^i_R)/⟨h_i⟩`, `I^i_{L,R}` being the components of `Ω_{L,R}`
stabilised by the peripheral element `h_i`. A *boundary Cauchy slice* is, in each asymptotic boundary, an
acausal Cauchy curve of the boundary cylinder; its lift to `I^i_L × I^i_R` is the graph of an increasing homeomorphism
`ψ_i : I^i_L → I^i_R` commuting with `h_i`.

**Marking.** Let `φ̃ : S¹ → S¹` be the `ρ_R ∘ ρ_L^{-1}`-equivariant homeomorphism equal to the boundary map of
the marking on `Λ_L` and to the translates of the `ψ_i` on `Ω_L`. It is a monotone bijection of the circle, hence a homeomorphism; it descends to
`φ : ∂_∞X_L → ∂_∞X_R` and extends equivariantly to the disc (e.g. by the Douady–Earle extension), giving a
homeomorphism `X̄_L → X̄_R` in the isotopy class of the marking. A
finite set of points of the boundary slice gives `P_L` and `P_R = φ(P_L)`; regions of the slice are regions
(§1.1) for both surfaces; and for boundary points `p, q` of the slice with lifts `(ũ_p, ṽ_p)`, `(ũ_q, ṽ_q)`,
`ṽ = φ̃(ũ)`, the arc class on `X_L` of the geodesic from `ũ_p` to `ρ_L(γ)ũ_q` corresponds under `φ` to the arc
class on `X_R` of the geodesic from `ṽ_p` to `ρ_R(γ)ṽ_q`. This fixes the identification of cut classes on the
two surfaces, including boundary twists.

### 6.2 Proposition (length formula) **[proved]**

Boundary points of `AdS_3` are rank-one matrices up to nonzero scale, `p = a b^T`; put `u = [a]`,
`v = [Jb]` with `J = [[0,1],[−1,0]]`, so that `g_L p g_R^{-1}` corresponds to `(g_L u, g_R v)`.
`ω(a,c) := a^T J c`.

(a) For `γ ≠ 1`, `(ρ_L(γ), ρ_R(γ))` preserves a spacelike geodesic of `AdS_3` on which it acts by translation
of length `½(ℓ_L(γ) + ℓ_R(γ))`. This geodesic joins the attracting fixed point `(u^+, v^+)` to the repelling
fixed point `(u^−, v^−)`, both on the limit set of `Γ`, and lies in `D̃` (which contains the convex hull of the
limit set). If `ℓ_L(γ) ≠ ℓ_R(γ)` there is exactly one other invariant spacelike geodesic, with translation
length `½|ℓ_L(γ) − ℓ_R(γ)|` (the inner horizon in BTZ); it joins `(u^+, v^−)` to `(u^−, v^+)` and is **not** in
`D̃`. If `ℓ_L(γ) = ℓ_R(γ)` the second geodesic is pointwise fixed and there are further invariant null
geodesics. In all cases the first geodesic is the **only** spacelike geodesic of `D̃` that is invariant under
`(ρ_L(γ), ρ_R(γ))` and translated non-trivially by it; this is what Theorem H uses.

*Proof.* Conjugate to `ρ_L(γ) = diag(e^{a/2}, e^{−a/2})`, `ρ_R(γ) = diag(e^{b/2}, e^{−b/2})`, `a = ℓ_L`,
`b = ℓ_R > 0`. Then `x ↦ g_L x g_R^{-1}` multiplies the entries `x_{11}, x_{12}, x_{21}, x_{22}` by
`e^{(a−b)/2}, e^{(a+b)/2}, e^{−(a+b)/2}, e^{(b−a)/2}`. The curve `s ↦ [[0, e^{s}],[−e^{−s}, 0]]` is a unit-speed
spacelike geodesic (`−det(x') = 1`, `x'' = x`), translated by `(a+b)/2`; its endpoints are the matrix units
`E_{12}, E_{21}`, the most expanded and most contracted directions, i.e. the attracting and repelling fixed
points. The curve `s ↦ diag(e^{s}, e^{−s})` is a unit-speed spacelike geodesic translated by `(a−b)/2`. A
geodesic is the intersection of the quadric with a 2-plane through the origin; an invariant geodesic spans
an invariant 2-plane, and for `a ≠ b` the four eigenvalues are distinct, so the invariant 2-planes are the
six coordinate planes, of which only `{x_{11} = x_{22} = 0}` and `{x_{12} = x_{21} = 0}` meet `det = 1` in a
spacelike curve. Finally, for `x = diag(e^{s}, e^{−s})`, `tr(adj(x) E_{12}) = 0`: every point of the second
geodesic is lightlike-related to the limit point `E_{12}`, hence lies outside the open domain `D̃` of points
not causally related to the limit set (argument supplied by the second referee). For `a = b`: a geodesic
translated non-trivially has as endpoints two rank-one eigen-matrices with different eigenvalues; the
rank-one eigen-matrices are `E_{12}` (eigenvalue `e^{a}`), `E_{21}` (`e^{−a}`) and `E_{11}, E_{22}` (eigenvalue
`1`), and `N(E_{12}, E_{11}) = N(E_{12}, E_{22}) = N(E_{21}, E_{11}) = N(E_{21}, E_{22}) = 0`, so the only pair
joined by a spacelike geodesic is `E_{12}, E_{21}`.

*All causal types, all `a, b > 0` (first referee of revision 5).* Geodesics of `AdS_3` are the intersections
with 2-planes of `M_2(R)`. A plane invariant under `X ↦ ρ_L(h) X ρ_R(h)^{-1}` is the sum of its intersections
with the eigenspaces (`RE_{12}`, `RE_{21}`, and `RE_{11}`, `RE_{22}` — the last two forming one eigenspace when
`a = b`). If it is not `span{E_{12}, E_{21}}` it does not contain both `E_{12}` and `E_{21}`, hence lies in
`span{E_{11}, E_{22}, E_{12}}` or in `span{E_{11}, E_{22}, E_{21}}`. On the first `N(·, E_{12}) ≡ 0`, on the second
`N(·, E_{21}) ≡ 0`. Since `E_{12}, E_{21}` represent limit points and `D̃` is defined by `N(·, q) > 0` for all limit
points `q` (§6.3), such a plane does not meet `D̃`. So `hull{E_{12}, E_{21}}` is the only `h`-invariant geodesic
in `D̃`, of any causal type. ∎

(b) Let `p = a b^T`, `q = c d^T` with `N := ω(a,c)·ω(b,d) > 0`. Then `x(t) = (e^{t} p + e^{−t} q)/√N`,
`t ∈ R`, is the unit-speed geodesic of `AdS_3` from `q` to `p`, and the length of its part where the
coefficient of `p` and that of `q` are both `≤ 1/ε` is `log N + 2 log(1/ε)`.

*Proof.* `det(s·ab^T + t·cd^T) = st·ω(a,c)ω(b,d)` (expand the determinant of a sum of two rank-one matrices).
So `det x(t) = 1`, `−det(x'(t)) = 1`, and `x'' = x` (geodesic of the quadric). The truncation is
`e^{|t|}/√N ≤ 1/ε`, of length `2 log(√N/ε)`. ∎

In `PSL(2,R)` the sign of `N` is a choice of representatives. For two distinct points of the graph of a
monotone map (points of the lifted boundary slice, or of the limit set) `ω(a,c) ≠ 0` and `ω(b,d) ≠ 0`, and
the vectors can be chosen along the graph so that `N > 0` for every pair.

In `H²`, with horocycles at `[a]`, `[c]` determined by the vectors `a, c`, the signed distance between the
horocycles is `2 log|ω(a,c)|` (Penner, *Comm. Math. Phys.* 113 (1987)); the square root of its exponential is
the lambda length). Hence, for every class `c` of cuts and compatible cut-offs — on `X_L` the horocycle at
`u_p` determined by `a_p`, on `X_R` that at `v_p` determined by `b_p`, where the spacetime cut-off is the scale
of `p = a_p b_p^T` —

`L^ε(c) = ½( ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c) )`,   (6.1)

where `L^ε(c)` is the total regulated length of the spacetime geodesics of class `c`: for a closed component in
the conjugacy class of `γ`, the geodesic of (a) in `D̃`; for an arc from `p` to `q` labelled by `γ`, the
geodesic of (b) from `(ũ_p,ṽ_p)` to `(ρ_L(γ)ũ_q, ρ_R(γ)ṽ_q)`. By Fact 3.1 horoball removal and horocycle
distance agree once both chiral horoballs are below `ε_0`; rescaling `p → p/ε` shrinks both. For inequalities
balanced at every endpoint the conclusion does not depend on the cut-off at all. Equation (6.1) is the length
formula of Appendix A of arXiv:2407.07165 in hyperbolic-surface language.

### 6.3 Definition of the HRT entropy in these spacetimes

Up to here nothing has been said about which extremal curves the HRT prescription admits. This subsection
fixes that as a definition; §§6.4–6.7 then *prove* that the resulting entropy is the minimum over cut classes.
No hypothesis remains.

**Positive lifts.** An *acausal circle* is the graph `C' ⊂ RP¹ × RP¹ = ∂AdS_3` of an orientation-preserving
homeomorphism. It has a continuous lift `p : C' →` (rank-one matrices), unique up to multiplication by a
positive function and a global sign, with `N(p,q) := tr(adj(p) q) > 0` for all distinct `p, q ∈ C'` (take
`a = (cos θ, sin θ)`, `b = (cos φ(θ), sin φ(θ))`, `θ ∈ [0, π)`, `φ` increasing; then
`N = sin(θ' − θ) sin(φ' − φ)`). `N` is twice the polarisation of `det`: `det(X + Y) = det X + det Y + N(X,Y)`.
For a closed set `Δ ⊂ C'` let `K(Δ)` be the set of finite positive combinations `X = Σ λ_i p_i`, `p_i ∈ Δ`,
`λ_i > 0`, involving at least two distinct points; then `det X = Σ_{i<j} λ_iλ_j N(p_i,p_j) > 0`, and
`[X] := X/√det X` is a point of `AdS_3`. `hull(Δ) := {[X] : X ∈ K(Δ)}`. For two points `hull` is the geodesic
between them (§6.2(b)); for three it is the totally geodesic *ideal triangle* they span.

**Spacelike separation.** `x, y ∈ AdS_3 = PSL(2,R)` are *spacelike separated* if `½|tr(x^{-1}y)| > 1`, i.e. if
they are joined by a spacelike geodesic. For `X, Y` in cones `K(·)` over the same acausal circle,
`tr([X]^{-1}[Y]) = N(X,Y)/√(det X det Y) > 0`.

*Caution (first referee of revision 5).* In `D̃` the criterion must be read with signs. `D̃` (defined next) has
a continuous *positive lift* to `SL(2,R)`, the one with `N(X,q) > 0` for the limit points `q`. Two points of
`D̃` with `tr(X^{-1}Y) < −2` in the positive lift are joined by a spacelike geodesic of `PSL(2,R)` that leaves
`D̃`, and can be causally related inside `D̃` (explicit example in rotating BTZ found by the referee). The
correct condition is `tr(X^{-1}Y) > 2` in the positive lift; Lemma 6.3a shows that slices satisfy it.

**The spacetime.** `D̃ := {[X] ∈ AdS_3 : N(X, q) > 0 for all q ∈ Λ_Γ}` (for one of the two lifts `±X`), the
invisible domain of the limit set; `M = D̃/Γ` (§6.1). For every acausal circle `C' ⊃ Λ_Γ`,
`hull(C') ⊂ D̃`: `N(Σλ_i p_i, q) = Σ λ_i N(p_i,q) > 0` since at most one `p_i` equals `q`.

**Boundary slices through `∂A`.** Fix the boundary Cauchy slice and marked points as in §6.1, giving the
`Γ`-invariant acausal circle `C` (the lifted slice together with `Λ_Γ`). For a region `A`, a *slice through
`∂A`* is any `Γ`-invariant acausal circle `C' ⊃ Λ_Γ` containing the lifts of the points of `∂A`; it carries the
region `A'` lying over the same intervals of the left ideal boundary as `A` (the domain of dependence of `A'`
in the boundary is that of `A`).

**Admissible bulk slices.** A bulk slice anchored on `C'` is a set `Σ̃ ⊂ D̃` such that
(S1) `Σ̃ ∪ C'` is a closed topological disc embedded in `AdS_3 ∪ ∂AdS_3`, with interior `Σ̃` and boundary `C'`;
(S2) `Σ̃` is `Γ`-invariant;
(S3) any two distinct points of `Σ̃` are spacelike separated.
`Σ := Σ̃/Γ ⊂ M`.

**Lemma 6.3a ((S1) and (S3) give `tr > 2` in the positive lift; acausality).** Let `Σ̃ ⊂ D̃` satisfy (S1) and
(S3), and let `X, Y` denote positive lifts. Then `tr(X^{-1}Y) > 2` for all distinct `x, y ∈ Σ̃`; the spacelike
geodesic segment `{[sX + tY] : s, t ≥ 0}` joins them inside `D̃`; and `x, y` are not joined by any causal curve
in `D̃`.

*Proof.* `Σ̃` is an open disc, so the set of ordered pairs of distinct points of `Σ̃` is connected; on it
`τ(x,y) = tr(X^{-1}Y)` is continuous, `|τ| > 2` by (S3), and `τ → 2` as `y → x`. Hence `τ > 2` throughout.
Then `det(sX + tY) = s² + t² + stτ > 0` and `N(sX + tY, q) = sN(X,q) + tN(Y,q) > 0` for limit points `q`, so
the segment lies in `D̃`; it is the geodesic `cosh`-combination of `X` and `Y`, spacelike because `τ > 2`.
The last claim is proved in §6.3b (Theorem C, Step 1, with Lemma K): the two points are not causally related
in the conformal model `E`, a fortiori in `D̃`. No cited fact is used. ∎

**HRT candidates.** A *candidate* for `A` is a quadruple `(C', Σ̃, γ, R)`: `C'` a slice through `∂A`; `Σ̃` a bulk
slice anchored on `C'`; `γ ⊂ Σ` a finite disjoint union of complete geodesics of `M`, each *simple* (no self-intersection), each closed or having
both ends at points of `∂A`, with exactly one end at each point of `∂A`; and `R ⊂ Σ` a closed set, equal to the
closure of its interior, with frontier `γ`, such that the closure of its preimage in `Σ̃ ∪ C'` meets
`C' \ Λ_Γ` exactly in the lift of `A'`. (`γ` is an extremal curve anchored on `∂A`, and `R` is a spacelike
region with `∂R = γ ∪ A'`: the homology condition of Headrick–Takayanagi and of Wall's maximin.)

**Definition.** `S^ε(A) := (1/4G_N) · inf { L^ε(γ) : (C', Σ̃, γ, R) a candidate for A }`, with `L^ε` the regulated
length of §6.2.

*Remarks on the definition.* (i) (S3) makes the slice acausal (Lemma 6.3a). It does **not** make it smooth
or uniformly spacelike; see (ii). Theorem H below holds verbatim for any class of slices between "(S1), (S2)
only" and "(S1)–(S3) and piecewise totally geodesic with spacelike faces", because its proof of `≥` uses only
(S1), (S2) and its proof of `≤` produces slices of the latter kind. (ii) The slices produced in §6.6 are acausal and have spacelike faces, but they are bent along geodesics,
and at a closed component of `γ` with `ℓ_L ≠ ℓ_R` their faces accumulate and become asymptotically lightlike
(Remark after Proposition P). If the definition required the homology region to be *smooth* and *uniformly*
spacelike, an additional argument would be needed: rounding a crease is elementary where the bending is
locally finite, but near such closed components the slice has to be replaced by a boost-twisted smooth
strip. That is **not carried out here**; Theorem H is proved for the definition as stated, with (S3).
(iii) Components of `γ` are required to be simple, as for an HRT surface that is an embedded submanifold; the
proof of `≥` in Theorem H uses this (the first referee of revision 5 exhibited, for static BTZ, a
self-crossing arc plus the horizon with a checkerboard region, which the earlier wording admitted). For
self-intersecting collections with a mod-2 homology condition, as in arXiv:2407.07165, Lemma 6.8 shows that
the corresponding collections of curves *on the surface* never beat the minimum; the step from a spacetime
configuration that lies on no single slice to such a collection on the surface is **not written out**, so no
claim is made here for definitions of HRT that do not place the homology region in a slice.
(iv) *Comparison with the definitions in the literature* (wording checked against the sources, 15 Nov 2026).
Headrick–Hubeny–Lawrence–Rangamani (arXiv:1408.6300, §2.2 and footnotes 4, 10) require the homology region to
be "a region of a bulk Cauchy slice" of the spacetime with its boundary, and "require all Cauchy slices to be
acausal"; they call this the "spacelike homology" condition. Wall (arXiv:1211.3494, Definitions 7–8 and
footnote 9) uses "complete achronal" slices, "a Cauchy surface [being] permitted to be null". Neither demands
smoothness or a uniformly spacelike region. §6.3b proves that the sets satisfying (S1)–(S3) are exactly the
acausal Cauchy slices of the spacetime with its conformal boundary (Theorems C and E), so the definition
above *is* the first of these. Wall's larger class (achronal, possibly null slices) is **not** treated.

### 6.3b The physical formulation: acausal Cauchy slices of the spacetime with its boundary **[proved]**

This subsection proves that the slices of §6.3 are *exactly* the acausal Cauchy slices of the spacetime with
its conformal boundary, so that the definition of §6.3 is the HRT prescription with the spacelike homology
condition of Headrick–Hubeny–Lawrence–Rangamani (region of an acausal bulk Cauchy slice). All causal
statements are made with the signed criterion; nothing uses `½|tr| > 1`.

**The conformal model.** Let `D̄` be the closed hemisphere `{(χ,φ) : 0 ≤ χ ≤ π/2}` with the round metric and
distance `d`, and `E := D̄ × R ∋ (z, θ)` with the metric `−dθ² + dχ² + sin²χ dφ²`. With
`R(θ) = [[cos θ, sin θ],[−sin θ, cos θ]]`, `S(φ) = [[cos φ, sin φ],[sin φ, −cos φ]]` and `sinh ρ = tan χ`, the map
`(z,θ) ↦ X = cosh ρ R(θ) + sinh ρ S(φ)` is a covering of `SL(2,R)` by the interior `E°`, and
`−det(dX) = cos^{−2}χ (−dθ² + dχ² + sin²χ dφ²)` (`det` is the quadratic form with `R, R(θ+π/2)` of norm `1`,
`S, S(φ+π/2)` of norm `−1`, all four orthogonal). A boundary point `(φ, θ) ∈ ∂E` corresponds to the ray of the
rank-one matrix `b = R(θ) + S(φ)`. `PSL(2,R) ∪ ∂` is the quotient by `(θ,φ) ∼ (θ+π, φ+π)`. A *causal curve*
in `E` is a non-constant curve along which `θ(s') − θ(s) ≥ d(z(s), z(s'))` for `s < s'`; `J^±` are defined
accordingly. Since the metric is a product and `D̄` is geodesically convex, `y ∈ J^+(x)` iff
`θ_y − θ_x ≥ d(z_x, z_y)`.

**Lemma K (signed causal criterion).** Two points of `E` are *not* causally related iff `|Δθ| < d(z,z')`.
In terms of the matrices attached to the points: for two interior points iff `|Δθ| < π` and `N(X,Y) > 2`; for
an interior point and a boundary point iff `|Δθ| < π` and `N(X,b) > 0`; for two boundary points iff
`|Δθ| < π` and `N(b,b') > 0`.

*Proof.* `N(R(θ),R(θ')) = 2cos Δθ`, `N(S(φ),S(φ')) = −2cos Δφ`, `N(R,S) = 0`. Hence
`N(X,Y) = 2cosh ρ cosh ρ' cos Δθ − 2 sinh ρ sinh ρ' cos Δφ`, while
`cos d(z,z') = cos χ cos χ' + sin χ sin χ' cos Δφ = (1 + sinh ρ sinh ρ' cos Δφ)/(cosh ρ cosh ρ')`. So
`N(X,Y) > 2 ⇔ cos Δθ > cos d`. As `d ∈ [0,π]`: if `|Δθ| ≤ π`, `cos Δθ > cos d ⇔ |Δθ| < d`; if `|Δθ| > π ≥ d` the
points are causally related. The other two cases are the limits `ρ' → ∞`, and `ρ, ρ' → ∞`, after dividing by
`cosh ρ'` (and `cosh ρ`). ∎

**Lemma L (lifting acausal circles and the domain).** (i) An acausal circle `C'` with positive lift `p` lifts to
a closed curve `C̃' ⊂ ∂E`, the graph `θ = c(φ)` of a function with `|c(φ) − c(φ')| < d(φ,φ')` for `φ ≠ φ'`, with
`p ∈ R_+·(R(c(φ)) + S(φ))`; any two of its points are not causally related.
(ii) Let `Λ̃ ⊂ C̃` be the lift of `Λ_Γ` (`Γ` non-trivial, so `Λ_Γ` has at least two points) and
`𝔼_Γ := {x ∈ E : x is not causally related to any λ ∈ Λ̃}`. Then `𝔼_Γ` is open in `E` and causally convex (a
causal curve of `E` between two of its points lies in it), and `(z,θ) ↦ [X]` is a homeomorphism of `𝔼_Γ ∩ E°`
onto `D̃` whose inverse is given by the positive lift.

*Proof.* (i) Write `p(u) ∈ R_+(R(θ(u)) + S(φ(u)))`. `N(p,p') > 0` reads `cos Δθ > cos Δφ`, so `φ(u) ≠ φ(u')` for
`u ≠ u'`: `u ↦ φ(u)` is a continuous injection of a circle into a circle, hence a homeomorphism, and `θ` is
locally a function `c` of `φ` with `|Δc| < |Δφ|` for nearby points. A continuous lift `c : R → R` is therefore
1-Lipschitz with `c(φ + 2π) = c(φ) + 2πk`, `|k| ≤ 1`; if `|k| = 1` then `|c(φ+π) − c(φ)| = π`, contradicting
`cos Δc > cos π`. So `k = 0`, `|Δc| ≤ d` and, since `cos Δc > cos d`, `|Δc| < d`. Lemma K gives the last claim.
(ii) Open: `Λ̃` is compact and `{(x,λ) : |θ_x − θ_λ| ≥ d}` is closed. Causally convex: if `w` on a causal curve
from `x` to `y` were in `J^+(λ)` then so would `y` be; if in `J^−(λ)`, so would `x`. If `x ∈ 𝔼_Γ ∩ E°`, Lemma K
gives `N(X, λ) > 0` for all `λ`, so `[X] ∈ D̃` and `X` is its positive lift. Injectivity: another preimage of
`[X]` is `(θ + kπ, φ + kπ)`, `k ≠ 0`; for odd `k` its distance to `λ` is `π − d` and `|Δθ + kπ| ≥ π − |Δθ| > π − d`;
for even `k`, `|Δθ + kπ| > π ≥ d`; so it is causally related to `λ`. Surjectivity: let `X` be a positive lift of
a point of `D̃` and `λ_0 ∈ Λ̃`. `N(X,λ_0) > 0` gives `cos(θ − θ_{λ_0}) > −1`, so the preimage `x` of `X` can be
chosen with `|θ_x − θ_{λ_0}| < π`, and then `|θ_x − θ_{λ_0}| < d_0 := d(z_x, z_{λ_0})`. For another `λ`, with
`d_λ := d(z_x, z_λ)`: `N(X,λ) > 0` says that `θ_x − θ_λ` is within `d_λ` of a multiple `2πm` of `2π`. If `m ≠ 0`
then `|θ_λ − θ_{λ_0}| > 2π − d_λ − d_0`, whereas `|θ_λ − θ_{λ_0}| < d(z_λ, z_{λ_0}) ≤ 2π − d_0 − d_λ` (path through
the antipode of `z_x`). So `m = 0` and `x ∈ 𝔼_Γ`. ∎

**The spacetime with boundary.** `M̄ := 𝔼_Γ/Γ`. Its interior is `M = D̃/Γ`; its boundary `∂M̄ = (𝔼_Γ ∩ ∂E)/Γ`
consists of the boundary points not causally related to the limit set, i.e. (§6.1) the domains of
dependence of the components of `C \ Λ_Γ` — one cylinder per asymptotic region, on which the boundary theory
lives. This is the bulk-with-boundary of arXiv:1408.6300. A *Cauchy slice of `M̄`* is a subset met exactly
once by every inextensible causal curve of `M̄`; it is then acausal (a causal curve between two of its points
extends to an inextensible one). Causal curves of `M̄` are the projections of those of `𝔼_Γ`.

**Theorem C (slices are Cauchy).** Let `C'` be a `Γ`-invariant acausal circle containing `Λ_Γ` and let `Σ̃`
satisfy (S1)–(S3) with boundary `C'`. Then `Σ̄ := (Σ̃ ∪ (C' \ Λ_Γ))/Γ` is a Cauchy slice of `M̄`, and
`Σ̄ ∩ ∂M̄ = (C' \ Λ_Γ)/Γ` is a Cauchy slice of `∂M̄`.

*Proof.* Identify `Σ̃` with its preimage in `𝔼_Γ` (Lemma L) and let `G` be its closure in `E`, i.e.
`Σ̃ ∪ C̃'` (S1).
*Step 1: no two points of `G` are causally related.* For `x ≠ y` in `Σ̃`, `N(X,Y) > 2` in the positive lift
(Lemma 6.3a, first part). `N(X,Y) > 2` implies `cos Δθ > −1`; `Δθ` is continuous on the connected set of pairs
of distinct points of the disc `Σ̃` and tends to `0` at the diagonal, so `|Δθ| < π`; Lemma K. Two points of
`C̃'`: Lemma L(i). For `x ∈ Σ̃`, `b ∈ C̃'`: `b` is a limit of points of `Σ̃`, so `|θ_b − θ_x| ≤ d(z_x, z_b)`. The
same closure argument shows that `G` contains no two points on a vertical line and is the graph of a
1-Lipschitz function over its projection. If `|θ_b − θ_x| = d(z_x,z_b)`, the Lipschitz function has slope
exactly one along the minimising segment from `z_x` to `z_b`, whose points other than `z_b` are interior; two
of them give points of `Σ̃` with `|Δθ| = d`, contradicting the first case. So `|θ_b − θ_x| < d`.
*Step 2: `G` is the graph of a continuous function `f : D̄ → R`.* The projection `G → D̄` is continuous,
injective (Step 1: points on a vertical line are causally related), maps the boundary circle `C̃'`
homeomorphically onto `∂D̄` (Lemma L(i)) and the open disc `Σ̃` into the open disc; by invariance of domain
the image of `Σ̃` is open, and it is closed in the open disc because `G` is compact; so it is everything.
*Step 3: every inextensible causal curve `κ` of `𝔼_Γ` meets `G`.* Otherwise `κ` lies in `{θ < f}` or in
`{θ > f}`; by time reversal assume the former, and parametrise `κ` towards the future. `θ` increases along
`κ` and is bounded by `max f`, and `d(z(s),z(s')) ≤ θ(s') − θ(s)`, so `κ` converges to a point `q ∈ E` with
`θ_q ≤ f(z_q)`. `q ∉ 𝔼_Γ` because `κ` is inextensible. So `q` is causally related to some `λ ∈ Λ̃`. If
`q ∈ J^−(λ)`, then `κ ⊂ J^−(q) ⊂ J^−(λ)`, contradicting `κ ⊂ 𝔼_Γ`. If `q ∈ J^+(λ)` and `q ≠ λ`, the point
`q' = (z_q, f(z_q)) ∈ G` lies in `J^+(q) ⊂ J^+(λ)` and `q' ≠ λ` (else `q ∈ J^+(λ) ∩ J^−(λ) = {λ}`); but `q'` and
`λ` are both in `G`, contradicting Step 1.
*Step 4.* A point of `κ ∩ G` lies in `𝔼_Γ ∩ G = Σ̃ ∪ (C̃' \ Λ̃)`. By Step 1 there is exactly one. An
inextensible causal curve of `M̄` lifts to one of `𝔼_Γ`, and its intersections with `Σ̄` correspond to the
intersections of the lift with `Σ̃ ∪ (C̃' \ Λ̃)`, which is `Γ`-invariant: exactly one. The boundary statement
is the same argument in `∂E` (or: `C̃' \ Λ̃` is a graph over each gap of `∂D̄ \ z(Λ̃)`). ∎

**Theorem E (converse: every acausal Cauchy slice is of this form).** Let `Σ̄` be a Cauchy slice of `M̄` and
`Σ̂ ⊂ 𝔼_Γ` its preimage. Then the closure of `Σ̂` in `E` is `Σ̂ ∪ Λ̃`, `C̃' := (Σ̂ ∩ ∂E) ∪ Λ̃` is the lift of a
`Γ`-invariant acausal circle `C'`, and `Σ̃ := Σ̂ ∩ E°` satisfies (S1)–(S3) with boundary `C'`. Hence
`Σ̃ ↦ Σ̄` is a bijection between the slices of §6.3 and the Cauchy slices of `M̄`.

*Proof.* Two points of `Σ̂` are not joined by a causal curve of `𝔼_Γ` (its projection would meet `Σ̄` twice,
or be a closed causal curve through `Σ̄`, which extended periodically meets it infinitely often); by causal
convexity they are not causally related in `E`. So `Σ̂` is the graph of a function `g` over a subset `U ⊂ D̄`
with `|g(z) − g(z')| < d(z,z')`. For `z ∈ D̄ \ z(Λ̃)` the set `{θ : (z,θ) ∈ 𝔼_Γ} = ⋂_λ (θ_λ − d(z,z_λ), θ_λ + d(z,z_λ))`
is an open interval, non-empty (it contains the point of the graph of Theorem C for the pleated slice of
any cut class, e.g. the empty one with auxiliary points), and the vertical segment over it is an
inextensible causal curve of `𝔼_Γ`; it meets `Σ̂` exactly once. So `U = D̄ \ z(Λ̃)`. `g` is 1-Lipschitz, hence
extends continuously to `D̄`; since `|g(z) − θ_λ| < d(z, z_λ)`, the extension takes the value `θ_λ` at `z_λ`.
Therefore the closure of `Σ̂` is the graph of the extension `ḡ` — a closed disc — and equals `Σ̂ ∪ Λ̃`; its
boundary circle `C̃'` is the graph of `ḡ` over `∂D̄`, with `|Δḡ| < d` for all pairs (pairs in `Σ̂`: shown;
pairs involving `λ`: definition of `𝔼_Γ`, and Lemma L(i) for `λ, λ'`), so by Lemma K it projects to a circle
with a positive lift, i.e. an acausal circle. This is (S1). (S2): `Σ̂` is `Γ`-invariant. (S3) in the signed
form: Lemma K gives `N(X,Y) > 2` for distinct points of `Σ̃`, and the matrices are the positive lifts
(Lemma L(ii)). Bijectivity: Theorem C shows the map is well defined, and `Σ̃` is recovered as the interior
part of the preimage. ∎

**Corollary (the definition of §6.3 is the HRT prescription with spacelike homology).** By Theorems C and
E, a candidate `(C', Σ̃, γ, R)` in the sense of §6.3 is the same thing as: a Cauchy slice `Σ̄` of the spacetime
with boundary `M̄` (necessarily acausal) whose boundary passes through `∂A`; a finite disjoint union `γ ⊂ Σ̄`
of simple complete geodesics of `M`, closed or anchored at `∂A`; and a region `R ⊂ Σ̄` with `∂R = γ ∪ A'`. For
`C' ∋ A` this is verbatim the condition "`ℛ_A` is a region of a bulk Cauchy slice with `∂ℛ_A = A ∪ ℰ_A`" with
Cauchy slices "required to be acausal" (arXiv:1408.6300, §2.2, footnotes 4, 10). Restricting to slices that
contain `A` itself does not change the infimum: the lower bound of Theorem H holds for all slices through
`∂A`, and the candidate constructed for the upper bound has `C' = C ⊃ A`.

### 6.4 Separation Lemma **[proved]**

**Lemma 6.4.** Let `C'` be an acausal circle, `r ≠ r'` two of its points, and `α, β` the two closed arcs of
`C'` with endpoints `r, r'`. For all `X ∈ K(α)`, `Y ∈ K(β)`:

`N(X, Y) ≥ 2 √(det X · det Y)`,

with equality if and only if `[X] = [Y]` lies on the geodesic `hull{r, r'}`. Hence every point of `hull(α)` is
spacelike separated from every point of `hull(β)` other than itself.

*Proof.* Apply an element of `PSL(2,R) × PSL(2,R)` (an isometry, preserving everything in the statement)
taking `r` to `(∞, ∞)` and `r'` to `(0, 0)`. Then `C'` is the graph of an increasing homeomorphism `φ` of
`R ∪ {∞}` fixing `0` and `∞`; `α = {u ∈ [0, ∞]}`, `β = {u ∈ [−∞, 0]}`; and a positive lift is
`p(u) = (u, 1)^T (φ(u), 1) = [[uφ(u), u], [φ(u), 1]]` for finite `u`, `p(∞) = E_{11}`
(`N(p(u), p(u')) = (u − u')(φ(u) − φ(u')) > 0`, `N(E_{11}, p(u)) = 1`). Therefore

`X = [[x_{11}, x_{12}], [x_{21}, x_{22}]]`, `Y = [[y_{11}, −y_{12}], [−y_{21}, y_{22}]]`, with all `x_{ab}, y_{ab} ≥ 0`,

because on `α` both `u` and `φ(u)` are `≥ 0` and on `β` both are `≤ 0`. Then
`N(X,Y) = tr(adj(X) Y) = x_{22}y_{11} + x_{11}y_{22} + x_{12}y_{21} + x_{21}y_{12}`, and

`x_{22}y_{11} + x_{11}y_{22} ≥ 2√(x_{11}x_{22} · y_{11}y_{22}) ≥ 2√(det X · det Y)`,

since `x_{11}x_{22} = det X + x_{12}x_{21} ≥ det X` and likewise for `Y`. Equality forces `x_{12}x_{21} = 0`,
`y_{12}y_{21} = 0`, `x_{12}y_{21} = x_{21}y_{12} = 0` and `x_{22}y_{11} = x_{11}y_{22}`. Now `x_{12} = Σ λ_i u_i`
vanishes iff every `p_i` occurring in `X` is `p(0)` or `p(∞)`, and the same for `x_{21} = Σλ_i φ(u_i)`; so
`x_{12}x_{21} = 0` means `X` is a combination of the two cut points only, i.e. `[X] ∈ hull{r,r'}`, and then
`X = diag(x_{11}, x_{22})`. Likewise `Y = diag(y_{11}, y_{22})`, and `x_{22}y_{11} = x_{11}y_{22}` says
`X ∝ Y`. ∎

**Lemma 6.4a (one triangle).** For three distinct points `p_1, p_2, p_3` of an acausal circle, `det` restricted
to their span has signature `(1,2)`, the ideal triangle `hull{p_1,p_2,p_3}` lies in the totally geodesic
spacelike plane `{det = 1} ∩ span`, and any two distinct points of it are spacelike separated.

*Proof.* The Gram matrix of `N` on `p_1,p_2,p_3` has zero diagonal and positive off-diagonal entries; its
trace is `0` and its determinant `2N_{12}N_{13}N_{23} > 0`, so its eigenvalues have signs `(+,−,−)`. On a space
of signature `(1,2)`, vectors with `det > 0` in the same component of `{det > 0}` satisfy the reversed
Cauchy–Schwarz inequality `N(X,Y) ≥ 2√(det X det Y)` with equality iff proportional; `K{p_1,p_2,p_3}` is
connected and contained in `{det > 0}`. ∎

**Prior work (found by the first referee of revision 5; statements checked against the introduction of the
paper, Section 3 only through a summary).** Mazzoli–Viaggi, arXiv:2206.06946, prove for maximal
representations of surface groups into `SO_0(2,n+1)` (Theorem A): "For every maximal lamination `λ`, the pleated
set `Ŝ_λ ⊂ H^{2,n}` is an embedded Lipschitz disk which is also acausal, that is, every pair of points
`x, y ∈ Ŝ_λ` is joined by a spacelike geodesic"; their Lemma 3.6 is, in substance, the two-geodesic case of
Lemma 6.4 (geodesics with unlinked endpoints on an acausal curve have acausal union), and their
Propositions 3.7–3.8 treat pleated sets over an acausal curve. With `n = 1` this is the same phenomenon as
Lemma 6.4 and Proposition P, **which are therefore not new in substance**. What differs here: the boundary
curve is a boundary Cauchy slice strictly larger than the limit set, the group is a free group acting with
funnels, the family contains arcs ending at boundary points and half-plane fillings, and the proofs (sign
pattern + arithmetic–geometric mean; piecewise-affine chart) are self-contained. Their theorem also covers
laminations with closed leaves and leaves spiralling onto them, which is independent support for the
closed-geodesic case of Lemma T / Proposition P.

Lemma 6.3' of the previous revision (two geodesics with unlinked endpoints) is the special case of Lemma 6.4
in which `X` and `Y` are combinations of two points each.

### 6.5 Lemma T (invariant spun families of ideal triangles containing a given cut) **[proved, given one cited topological fact]**

Let `c ∈ 𝒞(A)` with geodesic representative `c^X` on `X_L` (its components are pairwise distinct simple
geodesics, none null-homotopic); let `𝓕` be the set of lifts of its components to `H²` (chords of `S¹`): axes
of hyperbolic elements for the closed components, chords between lifted marked points for the arcs.

**An obstruction, and what the lemma therefore does not claim.** Let `A ∈ 𝓕` be an axis, with stabiliser
`⟨h⟩ ⊂ Γ`. In *no* `Γ`-invariant family of ideal triangles with pairwise disjoint interiors is `A` an edge of
a triangle: if `A` were an edge of `τ`, then `hτ` would be a triangle of the family with the edge `A` on the
same side (`h` preserves the two sides of its axis), so `τ` and `hτ` would overlap near `A` unless `hτ = τ`;
but `h` fixes only the two endpoints of `A`, so it preserves no ideal triangle. Consequently a cut containing
a closed geodesic **cannot** be extended to a locally finite invariant ideal triangulation having that
geodesic as an edge. The family constructed below is of a different kind: an axis of `𝓕` is not an edge and
is not adjacent to any triangle; infinitely many triangles accumulate on it from each side, and `h` permutes
them without fixing any. In the quotient the triangles spiral around the closed geodesic.

**Lemma T.** There is a `Γ`-invariant set `𝒯` of ideal triangles of `H²`, with pairwise disjoint interiors,
such that
(T1) the edges of the triangles are pairwise non-crossing and none crosses a chord of `𝓕`;
(T2) every arc-chord of `𝓕` is an edge of two triangles; every axis `A ∈ 𝓕` is disjoint (in `H²`) from all
closed triangles, and neither endpoint pair of an edge equals the endpoint pair of an axis of `𝓕`;
(T3) the union of the closed triangles and of the axes in `𝓕` is `H²`;
(T4) `𝒯` is locally finite in `H²` minus the axes in `𝓕`, and every sequence of distinct triangles meeting a
fixed compact set subconverges, in the Hausdorff sense in the closed disc, to an axis of `𝓕`;
(T5) all vertices of the triangles lie in the left coordinates of points of `C`: they are lifted marked
points (original or auxiliary) or endpoints of axes of `𝓕`.

*Proof.* **Step 0 (auxiliary points, boundary chords).** Add auxiliary marked points so that every ideal
circle carries at least one (at least three in all if `Γ` is trivial), and let `V_Ω ⊂ Ω` be the set of lifts of all marked points; stabilisers of points
of `V_Ω` are trivial. On `X̄_L` the geodesics joining consecutive marked points of each ideal circle
("boundary chords"; a loop around the funnel if the circle carries one point) cut off regions whose lifts are
half-planes `H_e` over the elementary arcs `e = [v, v']`, with trivial stabiliser. An arc of `c^X` may coincide
with a boundary chord.

**Step 1 (half-planes).** For each `Γ`-orbit of elementary arcs choose a homeomorphism `e ≅ [0,1]` and
transport it by `Γ`; take the triangles with vertices `k2^{-n}, (2k+1)2^{-n-1}, (k+1)2^{-n}`. They have
disjoint interiors, are locally finite in `H_e`, and their union is `H_e`: a point of `H_e` in no triangle would
lie in the regions cut off by a nested sequence of chords over dyadic intervals shrinking to a point of `e`,
which is impossible.

**Step 2 (pieces and corners).** The lifts of `c^X` and of the boundary chords form a locally finite
`Γ`-invariant family of pairwise disjoint chords. Let `Ỹ` be the closure in `H²` of a complementary component
that is not a half-plane `H_e`, `Γ_Y` its stabiliser, `Y = Ỹ/Γ_Y` (a compact surface once its finitely many
ideal vertices are added). `Ỹ` is closed and convex; its frontier consists of *axis sides* (axes of `𝓕`) and
*chord sides* (arc-chords of `𝓕`, boundary chords), and the closure `D_Y` of `Ỹ` in the closed disc is a
closed disc. A *corner* of `Ỹ` is an endpoint of a chord side (a point of `V_Ω`) or an axis side. Corners
correspond to pairwise disjoint closed subsets of the circle `∂D_Y` (a point, or the closure of an axis side),
so they are cyclically ordered. For every axis side `A` of `Ỹ` choose one of its two endpoints,
`e(A, Ỹ) ∈ {a^−, a^+}`, equivariantly: `e(gA, gỸ) = g·e(A, Ỹ)` for `g ∈ Γ`. This is possible because the
stabiliser of `(A, Ỹ)` is generated by a hyperbolic element with axis `A`, which fixes both endpoints. (The
two sides of an axis may receive different choices.) Let `ε` map a corner in `V_Ω` to itself and an axis side
`A` to `e(A, Ỹ)`. Distinct axes of a discrete group have no common endpoint and `V_Ω ∩ Λ_Γ = ∅`; hence **`ε` is
injective and preserves the cyclic order**, and the endpoint of `A` other than `e(A,Ỹ)` is not in the image
of `ε`.

**Step 3 (a topological decomposition; the cited fact).** An *edge* in `Y` is an embedded arc whose ends are
ideal vertices or points of the closed geodesics in `∂Y` ("smooth circles"), with interior in the interior of
`Y`. **Standard fact** (existence of triangulations of a compact surface with prescribed vertices in the
interior and on the boundary: Hatcher, *Topology Appl.* 40 (1991) 189–194; equivalently, a maximal family of
pairwise disjoint, pairwise non-homotopic essential arcs on the surface `Y^♯` obtained by collapsing each
smooth circle to a point cuts it into triangles): there is a finite family of pairwise disjoint edges in
`Y`, containing the chord sides, such that, in the universal cover `Ỹ`,
(D1) every lifted edge joins two *distinct* corners, and distinct lifted edges join distinct pairs of
corners;
(D2) every complementary region of the lifted edges has exactly three corners, pairwise joined by the three
lifted edges in its frontier (the rest of its frontier consists of segments of axis sides at its corners);
every lifted edge that is not a chord side is in the frontier of exactly two regions, a chord side of one.
The hypotheses of the fact are that `Y^♯` has at least one vertex and is not a sphere with one or two
vertices, nor a disc with one or two vertices all on the boundary; these cases would require a component of
`c^X` or a boundary chord to be null-homotopic, or two distinct ones to be homotopic, which does not happen
for distinct geodesics. (D1) is the statement that the arcs are essential and pairwise non-homotopic,
together with: an element `g ≠ 1` of `Γ_Y` cannot fix or exchange two distinct corners — corners in `V_Ω` have
trivial stabiliser, two distinct axes have stabilisers meeting trivially, and `Γ` has no torsion.
Nothing in this step involves the closed geodesics as edges: on `Y^♯` they are vertices.

**Step 4 (straightening; proved here).** For a lifted edge `ẽ` joining corners `κ, κ'` let `σ(ẽ)` be the chord
of `H²` joining `ε(κ), ε(κ')`; for a region `t̃` with corners `κ_1, κ_2, κ_3` let `σ(t̃)` be the ideal triangle
with vertices `ε(κ_i)`. Then `σ(gẽ) = gσ(ẽ)`, `σ(g t̃) = gσ(t̃)` for `g ∈ Γ_Y`.

(a) *`σ(ẽ)` is a chord contained in `Ỹ`, is not an axis side, and has trivial stabiliser in `Γ`; `σ` of a chord
side is that chord side; `σ(t̃)` is a non-degenerate ideal triangle contained in `Ỹ` and disjoint in `H²` from
the axis sides.* Injectivity of `ε` gives non-degeneracy; convexity of `Ỹ` gives containment. An axis side has
only one endpoint in the image of `ε`, so no `σ(ẽ)` is an axis side, and a closed ideal triangle inside the
closed half-plane of `A` meets `A` in `H²` only if `A` is one of its edges. If `g ≠ 1` in `Γ` preserved
`σ(ẽ)`, then `σ(ẽ)` would be the axis of the hyperbolic element `g` (or `g²`), with both endpoints in `Λ_Γ`;
both are then of the form `e(A,Ỹ)`, `e(A',Ỹ)` with `A ≠ A'`, and the axis of `g` would share an endpoint with
the distinct axis `A` — impossible in a discrete group.

(b) *The chords `σ(ẽ)` are pairwise distinct and pairwise non-crossing.* Two distinct lifted edges are
disjoint arcs in the closed disc `D_Y` with ends on its boundary circle, so their end pairs are unlinked
there (Jordan curve theorem); as `ε` preserves the cyclic order of corners, the endpoint pairs of the two
chords are unlinked or share points, i.e. the chords do not cross. They are distinct by (D1).

(c) *Triangles have disjoint interiors.* Let `t̃ ≠ t̃'` be regions. `t̃'` lies in a component of `D_Y` minus
the closure of `t̃`, which is cut off by one lifted edge `ẽ` of `t̃`; so the corners of `t̃'` lie in the closed
arc of `∂D_Y` beyond `ẽ`, and the third corner of `t̃` in the complementary one. Applying `ε`, `σ(t̃)` and
`σ(t̃')` lie in the two opposite closed half-planes of `σ(ẽ)`.

(d) *Local finiteness off the axes.* Let `K ⊂ H²` be compact and disjoint from all axes of `𝓕`; it is then at
distance `≥ δ > 0` from them (they form a locally finite family). Suppose infinitely many distinct chords
`g_nσ(ẽ)` of one `Γ_Y`-orbit meet `K`, at points `g_n y_n`, `y_n ∈ σ(ẽ)`. The image of `y_n` in `X_L` lies in the
image of `K`, which is compact and at distance `≥ δ` from the closed geodesics of `c^X`. But each end of
`σ(ẽ)` either ends at a point of `V_Ω ⊂ Ω`, in which case it leaves every compact subset of `X_L`, or ends at
an endpoint of an axis side `A`, in which case it is asymptotic to `A` and its image comes within `δ` of the
closed geodesic. So the `y_n` stay in a compact segment `J` of `σ(ẽ)`, and `g_n J ∩ K ≠ ∅` for infinitely many
distinct `g_n` (distinct because the stabiliser of `σ(ẽ)` is trivial), contradicting proper discontinuity.
There are finitely many orbits, so `K` meets finitely many chords `σ(ẽ)`.

(e) *The triangles cover `Ỹ` minus its axis sides.* Let `z ∈ Ỹ`, not on an axis side, and suppose `z` lies in
no closed triangle `σ(t̃)`. Fix a region `t̃_0` and an interior point `w_0` of `σ(t̃_0)`. The complement of a
closed ideal triangle consists of three disjoint open half-planes, one beyond each edge; so `z` lies strictly
beyond exactly one edge `σ(ẽ_0)` of `σ(t̃_0)`. `ẽ_0` is not a chord side, since nothing of `Ỹ` lies strictly
beyond a chord side. Let `t̃_1` be the other region with `ẽ_0` in its frontier (D2); by (c) `σ(t̃_1)` lies in
the closed half-plane beyond `σ(ẽ_0)`. Repeat: `z` lies strictly beyond exactly one edge `σ(ẽ_1)` of
`σ(t̃_1)`, and `ẽ_1 ≠ ẽ_0` because `z` and `σ(t̃_1)` are on the same side of `σ(ẽ_0)`. Inductively one gets
edges `ẽ_i` and open half-planes `H_0 ⊋ H_1 ⊋ …` bounded by the chords `σ(ẽ_i)`, with `z ∈ H_i` and
`w_0 ∉ H_i` for all `i`. Hence the infinitely many distinct chords `σ(ẽ_i)` all meet the segment `[w_0, z]`.
This segment is compact, lies in `Ỹ`, and meets the frontier of `Ỹ` at most at `z`, which is not on an axis
side; so it is disjoint from all axes of `𝓕`, contradicting (d).

By (a), (c), (e) the triangles `σ(t̃)` have disjoint interiors and their union is `Ỹ` minus its axis sides.
(*Description, not used below:* the regions with a corner at an axis side `A` with stabiliser `⟨h⟩` form an
`h`-periodic sequence along `A`; their images are the triangles `(e(A), v_k, v_{k+1})`, `k ∈ Z`,
`v_{k+m} = h v_k`, which accumulate on `A` as `k → −∞` and shrink to the ideal point `e(A)` as `k → +∞` when
`e(A)` is the attracting point of `h`.)

**Step 5 (assembly).** Carry out Steps 2–4 for one piece in each `Γ`-orbit and transport by `Γ`; let `𝒯` be
the set of all triangles `σ(t̃)` together with the half-plane triangles of Step 1. `𝒯` is `Γ`-invariant.
Interiors are disjoint (different pieces and half-planes have disjoint interiors). (T1): edges lie in
closed pieces or half-planes, whose interiors no chord of `𝓕` or boundary chord meets; within a piece, (b).
(T2): an arc-chord is a chord side of two pieces, or of a piece and a half-plane, and is an edge of one
triangle in each; axes: (a). (T3): Step 1 and (e). (T5): by construction. (T4): local finiteness off the
axes holds in each half-plane (Step 1) and each piece (d), and the pieces and half-planes are locally finite.
If distinct triangles `τ_n` meet a compact set, pass to a subsequence whose vertices converge; the limit is
not a non-degenerate triangle (two members of the sequence would then both contain a fixed interior point of
the limit), and it is not a point of the circle; so it is a chord `λ` through a point of the compact set,
and it is a limit of infinitely many distinct edges. By local finiteness off the axes, `λ` contains no point
off the axes of `𝓕`, i.e. `λ` is an axis of `𝓕`. ∎

### 6.6 Proposition P (pleated slices) **[proved, given Lemma T]**

Let `c`, `𝓕`, `𝒯` be as in Lemma T and `C` the acausal circle of the boundary slice, `p : C →` positive lift,
`u ↦ p(u)` its parametrisation by the left coordinate. For a triangle `τ = (u_1,u_2,u_3) ∈ 𝒯` put
`T(τ) = hull{p(u_1),p(u_2),p(u_3)}`, for an axis `A = (a^−,a^+) ∈ 𝓕` put `g(A) = hull{p(a^−),p(a^+)}`, and

`Z := ⋃_{τ ∈ 𝒯} T(τ) ∪ ⋃_{A} g(A) ⊂ AdS_3`.

**Proposition P.** `Z` is a bulk slice anchored on `C` in the sense of §6.3 — (S1), (S2), (S3) hold and
`Z ⊂ D̃` — it is a union of totally geodesic spacelike ideal triangles, and it contains the spacetime geodesic
`hull{p(u), p(u')}` of every chord `(u,u') ∈ 𝓕`.

*Proof.* `Z ⊂ hull(C) ⊂ D̃` (§6.3), and `Z` is `Γ`-invariant because `𝒯`, `𝓕` and `p` (up to positive scale) are.
Every chord of `𝓕` is an edge of a triangle or an axis, so its hull is in `Z`.

(S3). Let `x ≠ y` in `Z`, in pieces `T, T'` (triangles or axes). If `T = T'`, Lemma 6.4a (or §6.2(b) for an
axis). Otherwise the vertex sets of `T` and `T'` are non-crossing (T1), so all vertices of `T'` lie in one
closed arc of `C` between two vertices `r, r'` of `T` that are consecutive among the vertices of `T`, and the
vertices of `T` lie in the complementary closed arc. Lemma 6.4 with cut points `r, r'` gives spacelike
separation, since `x = y` is excluded by assumption.

(S1). Use the Klein model, in which `H²` is the open unit disc `D`, chords are straight segments, and the
triangles of `𝒯` are Euclidean triangles inscribed in the unit circle. Choose a linear functional `ℓ` on
`M_2(R)` positive on `p(C)` (e.g. `ℓ = N(·, p_1 + p_2 + p_3)` for three points of `C`) and normalise the lift
so that `ℓ ∘ p ≡ 1`; work in the affine chart `{ℓ = 1}` of projective space, in which `AdS_3 ∪ ∂AdS_3` near
`hull(C)` is a subset of `R³`, geodesics are segments and totally geodesic planes are affine planes. Define
`Ψ : D̄ → {ℓ = 1}` by `Ψ(u) = p(u)` on the circle, by the affine map with `Ψ(u_i) = p(u_i)` on each triangle
`(u_1,u_2,u_3) ∈ 𝒯`, and by the affine map onto the segment `[p(a^−), p(a^+)]` on each axis. On a common edge
of two triangles the two definitions agree (an affine map of a segment is determined by its endpoints). By
(T3) `Ψ` is defined on all of `D̄`, and `Ψ(D) = Z`.

*Continuity.* At a point of `D` off the axes, by local finiteness (T4). At `z` on an axis `A = (a^−, a^+)`,
`z = s a^+ + (1−s) a^−`: if `z_n → z` with `z_n` in triangles `τ_n`, then either `τ_n` is eventually constant
or, by (T4), the vertices of `τ_n` tend to `{a^−, a^+}`; writing `z_n` in Euclidean barycentric coordinates of
`τ_n`, the total weight on vertices tending to `a^+` tends to `s`, so
`Ψ(z_n) → s p(a^+) + (1−s) p(a^−) = Ψ(z)` by continuity of `p`. At `u` on the unit circle: if
`z_n = Σ_i s_i^{(n)} v_i^{(n)} → u`, then since `|z_n|² = 1 − Σ_{i<j} s_i s_j |v_i − v_j|²` and `|u| = 1`,
`Σ_i s_i^{(n)} |v_i^{(n)} − u|² → 0`, hence `Ψ(z_n) = Σ_i s_i^{(n)} p(v_i^{(n)}) → p(u)` by uniform continuity
of `p`.

*Injectivity.* `Ψ` is injective on the circle, and on each closed triangle and axis (an affine map between
non-degenerate simplices: three distinct points `p_1, p_2, p_3` of `p(C)` are linearly independent, because
`det(λp_1 + μp_2) = λμ N(p_1,p_2)` vanishes on the plane they span only on the two lines `Rp_1`, `Rp_2`, so
the rank-one matrix `p_3` is not in that plane). For two distinct pieces `T, T'` (triangles or axes) with cut points `r, r'` as in the proof
of (S3), the equality case of Lemma 6.4 says: if `Ψ(x) = Ψ(y)` with `x ∈ T`, `y ∈ T'`, then, written as positive
combinations of the vertices of their pieces, both have weight only on vertices in `{r, r'}`. So the two
images meet only in the image of the common face spanned by the common vertices among `r, r'`: a common
edge (both common, where the two affine definitions of `Ψ` agree), a common ideal vertex (a point of `p(C)`,
not in `Z`), or nothing. In particular, for an axis `A = (a^−, a^+)` and a triangle `τ`: by (T2) `τ` does not
have both `a^−` and `a^+` as vertices, so `T(τ) ∩ g(A) = ∅`; this is the case of the triangles spiralling
towards `A`, which share with `A` at most the ideal point `e(A)`. Two distinct axes have no common endpoint,
so their images are disjoint. Images of points of `D` do not lie on `p(C)` (they have `det > 0`).

A continuous injective map from the compact disc `D̄` to a Hausdorff space is an embedding. So `Z ∪ C = Ψ(D̄)`
is an embedded closed disc with boundary `C` and interior `Z`. ∎

**Remark (the pleated slice is asymptotically lightlike at rotating closed geodesics).** Let `A` be an axis of
`𝓕` with `(ρ_L(h), ρ_R(h)) = (diag(e^{a/2}, e^{−a/2}), diag(e^{b/2}, e^{−b/2}))` in the coordinates of §6.2(a),
so that `g(A)` joins `E_{21}` and `E_{12}`. The triangles of `Z` adjacent to `g(A)` on one side are
`hull{E_{12}, ρ(h)^k p(v), ρ(h)^{k+1} p(v)}`, `k → −∞`. Writing `p(v)` in the eigenbasis, modulo `E_{12}` and up
to scale `ρ(h)^k p(v) = E_{21} + e^{−|k|a} c_{11} E_{11} + e^{−|k|b} c_{22} E_{22}`, so the planes of these
triangles converge to `span{E_{12}, E_{21}, E_{11}}` if `a < b`, to `span{E_{12}, E_{21}, E_{22}}` if `a > b`,
and are all equal to one spacelike plane if `a = b`. On `span{E_{12}, E_{21}, E_{ii}}` the form `det` is
degenerate: these are the two lightlike planes through `g(A)`. Thus, when `ℓ_L(h) ≠ ℓ_R(h)` — a closed geodesic
carrying angular momentum — `Z` is acausal (S3) and each of its faces is spacelike, but the faces accumulating
at `g(A)` tilt towards a lightlike plane: `Z` is **not uniformly spacelike** near `g(A)`. This is forced for
any invariant surface that is totally geodesic on one side of `g(A)`: `h` acts on the normal plane of `g(A)`
by a boost of rapidity `(a−b)/2`, and the only boost-invariant normal directions are null. An invariant
surface that is smooth and spacelike along `g(A)` must have a tangent plane that is boosted steadily along
the geodesic; such strips exist locally, but fitting one into a global invariant slice is not done here. This
is the content of the caveat about smooth slices in §6.3.

### 6.6a Proposition S (smooth spacelike slices for cuts made of closed geodesics) **[proved, given cited theorems whose hypotheses are listed]**

The pleated slice is not smooth and not uniformly spacelike at rotating closed geodesics. For cuts without
arcs — the case of regions that are unions of whole boundaries — a smooth spacelike slice exists. (Route
suggested by the first referee of revision 5.)

**Proposition S.** Let `c ∈ 𝒞(A)` consist of closed geodesics only, with axes `𝓕`. Then there is a set `S̃`
satisfying (S1)–(S3) which is a smooth spacelike surface and contains `g(A')` for every axis `A' ∈ 𝓕`.

*Cited inputs.*
(C1) For an acausal circle `C`, the invisible domain `Ω(C) := {[X] : N(X,p) > 0 for all p ∈ C} ⊂ D̃` is the
domain of dependence of every properly embedded acausal disc with boundary `C`; in particular it is globally
hyperbolic and such a disc is a Cauchy surface of it. (C2) Every Cauchy surface of `Ω(C)` is a disc properly
embedded in `AdS_3` whose closure in `AdS_3 ∪ ∂AdS_3` is a closed disc with boundary `C`. (Mess; Barbot;
Bonsante–Seppi, arXiv:2004.14414, §4. **Not re-checked against the sources.**)
(C3) Bernal–Sánchez, arXiv:gr-qc/0512095, Theorem 1.1 (wording checked): "Let `(M,g)` be a globally
hyperbolic spacetime, and let `H ⊂ M` be a spacelike and acausal compact `m`-submanifold with boundary. Then,
there exists a spacelike Cauchy hypersurface `S` such that `H ⊂ S`." It is used with `m = 1` and `∂H = ∅`; that
the empty boundary is allowed is the usual convention and consistent with their Remark 4.14 on
non-connected `H`, but was **not verified in their proof**.

*Proof.* Let `Z` be the pleated slice of Proposition P for `c`. `Z ⊂ Ω(C)`: for `X = Σλ_i p(u_i)` and `p ∈ C`,
`N(X,p) > 0`. By (C1) `Z` is a Cauchy surface of `Ω(C)`. `Γ` preserves `C`, hence `Ω(C)`, and acts on it freely
and properly discontinuously (it does on `D̃`). `N := Ω(C)/Γ` is globally hyperbolic with Cauchy surface `Z/Γ`:
an inextensible causal curve of `N` lifts to one of `Ω(C)`, which meets `Z` exactly once, so the curve meets
`Z/Γ` exactly once.
Let `H ⊂ N` be the image of `⋃_{A' ∈ 𝓕} g(A')`: finitely many disjoint simple closed spacelike geodesics, a
compact 1-submanifold. `H` is acausal in `N`: a causal curve of `N` between two points of `H` would lift to a
causal curve of `Ω(C) ⊂ D̃` between two distinct points of `⋃ g(A') ⊂ Z`, contradicting Lemma 6.3a applied to
`Z`. By (C3) there is a smooth spacelike Cauchy hypersurface `S ⊃ H` of `N`. Its preimage `S̃ ⊂ Ω(C)` is smooth,
spacelike, `Γ`-invariant (S2), contains the `g(A')`, and is a Cauchy surface of `Ω(C)`: an inextensible causal
curve of `Ω(C)` projects to an inextensible causal curve of `N`, which meets `S` at exactly one parameter
value, so the curve meets `S̃` exactly once. By (C2) `S̃ ∪ C` is a closed disc with boundary `C` (S1). (S3): for
distinct `x, y ∈ S̃` with positive lifts and `τ = tr(X^{-1}Y)`: if `−2 < τ ≤ 2`, the segment `{[sX + tY]}` lies
in `Ω(C)` (`N(sX+tY, p) > 0`, `det = s² + t² + stτ > 0`) and is a timelike or lightlike geodesic from `x` to
`y`, contradicting acausality; so `τ ∉ (−2, 2]`; `τ → 2` near the diagonal and the set of pairs of distinct
points of the disc `S̃` is connected, so `τ > 2` for all pairs. ∎

**Corollary (Theorem H with smooth spacelike homology regions, `∂A = ∅`).** If `A` is a union of whole
boundary components, the infimum in the definition of `S^ε(A)` is unchanged when slices are in addition
required to be smooth spacelike surfaces: it is `≥` the unrestricted infimum, and Proposition S applied to a
minimising class (which has no arcs), with the homology region defined by the parity rule of §6.7, gives a
candidate with a smooth slice. Given (C1)–(C3). For regions with endpoints the smooth version remains
**open**: the arcs are not compact and (C3) does not apply.

### 6.7 Theorem H (the HRT entropy is the minimum over cut classes) **[proved, given Lemma T]**

**Theorem H.** For the spacetimes of §6.1 and every region `A`,

`S^ε(A) = (1/4G_N) · min_{c ∈ 𝒞(A)} ½( ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c) )`,

and the infimum in the definition of `S^ε(A)` is attained.

*Proof.* `≥` (every candidate is a cut class). Let `(C', Σ̃, γ, R)` be a candidate. A closed component of `γ` in
the free homotopy class of `g ∈ π_1(S)` lifts to a spacelike geodesic of `D̃` invariant under
`(ρ_L(g), ρ_R(g))` and translated non-trivially by it; by §6.2(a) this is `hull{p(a^−), p(a^+)}`, the
geodesic joining the repelling and attracting fixed points, which lie in `Λ_Γ ⊂ C'`. An arc component lifts
to the geodesic joining two lifted points of `∂A`, also on `C'`. Since `γ` is a finite union of
disjoint simple complete geodesics, each closed or running out to the boundary at both ends, it is a
properly embedded 1-manifold in `Σ`, and its preimage is one in `Σ̃`. So the preimage `γ̃ ⊂ Σ̃` of `γ` is a
`Γ`-invariant, locally finite family of pairwise disjoint arcs, each closed in the closed disc `Σ̃ ∪ C'` with
its two endpoints on `C'`. By the Jordan curve theorem in the closed disc, each such arc separates it, and
two disjoint ones have unlinked endpoint pairs. The left coordinate identifies `C'` with `S¹` equivariantly
and preserves cyclic order, so the chords of `H²` with the corresponding endpoints are pairwise
non-crossing; they are the lifts of the `X_L`-geodesics of finitely many classes (the same group elements
and the same lifted endpoints), hence form a `Γ`-invariant locally finite family, whose quotient is an
embedded geodesic multi-curve `c^X` on `X_L` with one arc-end at each point of `∂A`.

*Homology region.* Fix an elementary piece `e_0 ⊂ C' \ Λ_Γ`. For a point `z` of `Σ̃ \ γ̃` (or of `H²` off the
chords) let `n(z)` be the number of arcs of `γ̃` (chords) separating `z` from `e_0`; it is finite by local
finiteness. Since `R̃` is the closure of its interior with frontier `γ̃`, membership in `R̃` changes exactly
when an arc is crossed, so `z ∈ R̃ ⇔ n(z) + [e_0 ⊂ Ã'] ≡ 1 (mod 2)`. The same parity, computed with the chords
in `H²`, defines a set `W̃ ⊂ H²`, the closure of its interior, with frontier the chords. Whether an arc
(chord) separates an elementary piece from `e_0` depends only on the cyclic order of endpoints, so `W̃` and
`R̃` meet the ideal boundary in corresponding pieces: exactly the lifts of `A`. The parity function is
`Γ`-invariant: `n(gz) − n(z)` is, mod 2, the number of arcs separating `e_0` from `g e_0`, which is
`[g e_0 ⊂ Ã'] − [e_0 ⊂ Ã'] = 0`. Hence `W = W̃/Γ` is a homology region and `c := [c^X] ∈ 𝒞(A)`. Distinct
components of `c^X` are distinct geodesics, and no component is null-homotopic, so `c` is a reduced
representative. By (6.1), `L^ε(γ) = ½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))`.

`≤` (every cut class is a candidate). Let `c ∈ 𝒞(A)`, `𝓕` its chords, `𝒯` as in Lemma T, `Z` as in Proposition
P; take `C' = C`, `Σ̃ = Z`. Let `γ̃ = ⋃_{(u,u') ∈ 𝓕} hull{p(u),p(u')} ⊂ Z`, and `γ = γ̃/Γ`: a finite disjoint
union of complete geodesics of `M`, closed for the axes (§6.2(a)) and joining points of `∂A` for the arcs,
one end at each point of `∂A`; they are pairwise disjoint because `Ψ` is injective and the chords are. Let
`W̃ ⊂ H²` be the lift of the homology region of `c^X` and `R̃ := Ψ(W̃) ⊂ Z`. Since `Ψ` is a homeomorphism of
`D̄` onto `Z ∪ C` restricting to `u ↦ p(u)` on the circle and mapping each chord of `𝓕` onto its hull, `R̃` is
closed in `Z`, is the closure of its interior, has frontier `γ̃`, and its closure meets `C \ Λ_Γ` in the lifts
of `A`. It is `Γ`-invariant as a set, because `W̃` is a union of triangles of `𝒯` and axes (its frontier
consists of edges and axes), `𝒯` is `Γ`-invariant, and `Ψ` maps the triangle with vertices `u_i` onto the
triangle with vertices `p(u_i)` for every triangle. So `(C, Z, γ, R̃/Γ)` is a candidate, of length
`½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))` by (6.1).

The minimum over `𝒞(A)` exists (§4), so the infimum is attained. ∎

**What Theorem H replaces.** In earlier revisions "HRT = minimum over classes" was Hypothesis H, with parts H1
(candidates are embedded cut classes) and H2 (every embedded class bounds a spacelike region). H1 is the `≥`
part of the proof; H2 is Proposition P. The remaining inputs are: the definition of §6.3; the description of
the spacetimes in §6.1 (cited); and one topological fact in Lemma T (existence of triangulations of a compact
surface with prescribed vertices; cited). The spiralling of triangles around closed geodesics is proved in
Lemma T, Step 4.

### 6.8 Lemma (non-embedded collections do not win) **[proved]**

Let `c'` be a finite collection of closed curves and arcs between marked points, transverse, whose union has
the mod-2 intersection property of a cut for `A` (i.e. `X \ c'` admits a 2-colouring changing across each
branch of `c'`, with the pieces of `A` on colour 1 and exactly one arc-end at each point of `∂A`). Then there
is `c ∈ 𝒞(A)` with `i(τ, c) ≤ i(τ, c')` for all test classes `τ`, hence `ℓ^ε_{X}(c) ≤ ℓ^ε_X(c') + o(1)` for
every hyperbolic structure `X`; with (3.1), `ℓ^ε_X(c) ≤ ℓ^ε_X(c')` for all `ε ≤ ε_0` when `c'` has the same
ends.

*Proof.* Replace each component of `c'` by its geodesic representative, deleting null-homotopic components.
The existence of a 2-colouring with the stated boundary behaviour is equivalent to: every closed curve, and
every proper arc between two elementary pieces, meets `c'` in a number of points whose parity is prescribed
by `A` (even for closed curves; for arcs, the parity of the number of its ends lying in `A`). These parities
are homotopy invariants of `c'`, so the geodesic collection again admits such a colouring; geodesics that
coincide are reduced mod 2 (an even number of copies is deleted; a multiply covered closed component counts
with its multiplicity). Let `W` be the closure of the colour-1 set.
At each crossing all half-edges are colour changes; resolve as in Lemma A. Nothing needs resolving at marked
points (one end each). The result is a cut `c` for `A` made of pieces of `c'`. For a test class `τ` take a
geodesic representative, pushed off coincident closed geodesics to a parallel copy and perturbed off the vertices as in Lemma A; it
realises `i(τ, c')`
and meets `c` in at most that many points. Integrate as in Theorem 1. ∎


### 6.9 Theorems 2 and 3' (covariant entropy inequalities in vacuum 2+1 gravity) **[proved]**

Let `M`, the boundary slice and the parties be as in §6.1, and `S^ε(A)` the HRT entropy of §6.3. By Theorem H
and (6.1), `4G_N S^ε(A) = F^ε(A)` for `k = 2`, `w_L = w_R = ½`.

**Theorem 2.** Every inequality `Σ α_l S(A_l) ≥ Σ β_r S(B_r)` admitting a contraction map holds for the HRT
entropies `S^ε`: exactly, for every `ε ≤ ε_0`, if it is balanced at every endpoint of the regions (in
particular when all parties are unions of whole boundary components); and with the difference of the two
sides tending to `+∞` as `ε → 0` otherwise.

**Theorem 3'.** If all parties are unions of whole boundary components, the vector of HRT entropies
`(S(A))_A` lies in the Ryu–Takayanagi holographic entropy cone (no cut-off is involved). For parties
containing intervals the corresponding statement is **not proved** (§5.4).

*Proof.* Theorem H, then Theorems 1 and 3. ∎

**Hypotheses of Theorems 2 and 3', in full.**
* *Geometry:* `M = D̃/Γ` with `Γ = (ρ_L, ρ_R)(π_1 S)`, `S` an oriented surface of finite type with non-empty
  boundary; `ρ_L, ρ_R` convex cocompact Fuchsian of the second kind, holonomies of two hyperbolic structures
  on the interior of `S` with the same marking (§6.1). No time-symmetry or stationarity is assumed.
  Not allowed: extremal horizons, cusps, conical defects, non-orientable quotients, matter.
* *Regions:* finite unions of intervals of one acausal boundary Cauchy slice and of whole boundary
  components.
* *Regulator:* horoball (equivalently horocycle, equivalently boundary-point rescaling) cut-offs at the
  endpoints, below the threshold `ε_0` of §3.2.
* *Gravitational input:* none beyond the definition of the HRT entropy in §6.3 (minimal extremal curve
  anchored on `∂A`, subject to the spacelike homology condition on a bulk slice satisfying (S1)–(S3)). The
  null curvature condition, maximin and minimax are not used.
* *Cited mathematics:* the structure of the spacetimes (Mess; Barbot; Bonsante–Krasnov–Schlenker); measure
  zero of the limit set; the Crofton formula; existence of triangulations of marked surfaces and spinning
  (Lemma T); the graph characterisation and closedness of the cone (Bao et al.) for Theorem 3'.

## 7. Status table (§§1–5 after three referee passes; §§6.3–6.9: revision 6)

| Statement | Status |
|---|---|
| Lemma A (topological contraction lemma) | proved; mechanism due to Czech–Dong; three referee passes, no error found; its existence statement (A1) tested on the cylinder (the construction of `δ_r` is not executed in code) |
| Facts 3.1, 3.2; Lemma B (`κ = 1 − log 2`); (3.1), (4.1) | proved; constant and normalisation re-derived numerically by the referee |
| Theorem 1 (multi-metric RT, contraction-map inequalities, any regions) | proved |
| Lemma G (graph realisation) | proved; idea due to the referee; verified numerically on the disc by the referee |
| Lemma F (filling) | proved modulo the standard normal-coordinate argument |
| Theorem 3, `P = ∅` (full cone, whole-boundary parties) | proved |
| Theorem 3, `P ≠ ∅` (full cone, intervals) | **not proved**; missing steps listed in §5.4 |
| Length formula 6.2(b), (6.1), marking | proved |
| 6.2(a): outer geodesic in `D̃`, inner one not | proved (argument for the second part supplied by the referee) |
| Lemma 6.8 (non-embedded collections do not win) | proved |
| Lemma 6.4 (Separation Lemma), Lemma 6.4a | proved; tested numerically (160,000 cases, control fails) |
| Lemma T (invariant spun family of ideal triangles containing a cut; closed geodesics are not edges, triangles spiral onto them) | proved given one cited topological fact (triangulations of a compact surface with prescribed vertices); the spiralling step is proved in the text |
| Proposition P (pleated slice: embedded, invariant, pairwise spacelike separated, contains the geodesics) | proved, given Lemma T |
| §6.3b: Lemma K (signed causal criterion), Lemma L, Theorem C (slices of §6.3 are acausal Cauchy slices of the spacetime with boundary), Theorem E (converse) | written in revision 7; **not yet reviewed by any referee** |
| Theorem H (HRT entropy of §6.3 = minimum over cut classes) | proved, given Lemma T and the standard fact in Lemma 6.3a, **for the definition of §6.3** (homology region on an acausal, possibly non-smooth and asymptotically lightlike slice; simple components). One independent automated audit of revision 5: no error in the statement, gaps G2–G4 repaired in revision 6. Closed-geodesic case re-audited separately: see OPEN_GAPS.md G6 |
| Theorem 2 (HRT, contraction-map inequalities, arbitrary regions) | proved for the HRT entropy of §6.3 (acausal slices); with smooth spacelike slices only when no region has endpoints |
| Theorem 3' (HRT, full cone, whole-boundary parties) | proved for the HRT entropy of §6.3; also with smooth spacelike homology regions, given (C1)–(C3) of Proposition S |
| Lemma 6.4 / Proposition P: novelty | **not new in substance**: Mazzoli–Viaggi, arXiv:2206.06946 (Lemma 3.6, Theorem A); independent proofs here in a different setting |
| Proposition S (smooth spacelike slice through a cut made of closed geodesics) | proved given cited (C1)–(C3); (C1), (C2) not re-checked against sources; use of (C3) with empty boundary not verified in its proof |
| Variant of the definition with *smooth* spacelike slices, regions with endpoints | **not proved** |
| Matter; higher dimensions; cusps; conical defects; non-orientable quotients; endpoints not on one slice | not covered |

**Standing assumptions** (collected at the referee's request): every `X_j` is orientable with funnel ends
only; the identifications `φ_j` preserve peripheral structure, marked points and the party assignment;
weights `w_j > 0` (positivity is essential — the `max` functional violates the inequalities); the same class
is used in all metrics; regions are unions of elementary pieces on one boundary Cauchy slice; `ε_0` depends
on the geometry and on the horoball normalisation; coefficients `α_l, β_r > 0`; `ρ_L, ρ_R` are holonomies of
the same marked oriented surface; the HRT entropy is the quantity defined in §6.3 (the leading classical
term; no quantum corrections; slices satisfying (S1)–(S3)); graph models are
allowed several terminals per party; citations in §1.1 and §6.1 were not re-checked against the sources.
