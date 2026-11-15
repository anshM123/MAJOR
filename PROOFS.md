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
(the open convex set of points not causally related to any point of the limit set `Λ_Γ ⊂ ∂AdS_3`, which is
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
`D̃`. (If `ℓ_L(γ) = ℓ_R(γ)` the second geodesic is pointwise fixed and there are further invariant null
geodesics; none of this is used.)

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
not causally related to the limit set (argument supplied by the second referee). ∎

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
balanced at every endpoint the conclusion does not depend on the cut-off at all. Equation (6.1) is equations (58),
(60) of arXiv:2407.07165 in hyperbolic-surface language.

### 6.3 Definition of the HRT entropy in these spacetimes

Up to here nothing has been said about which extremal curves the HRT prescription admits. This subsection
fixes that as a definition; §§6.4–6.7 then *prove* that the resulting entropy is the minimum over cut classes.
No hypothesis remains.

**Positive lifts.** An *acausal circle* is the graph `C' ⊂ RP¹ × RP¹ = ∂AdS_3` of an orientation-preserving
homeomorphism. It has a continuous lift `p : C' →` (rank-one matrices), unique up to multiplication by a
positive function, with `N(p,q) := tr(adj(p) q) > 0` for all distinct `p, q ∈ C'` (take
`a = (cos θ, sin θ)`, `b = (cos φ(θ), sin φ(θ))`, `θ ∈ [0, π)`, `φ` increasing; then
`N = sin(θ' − θ) sin(φ' − φ)`). `N` is twice the polarisation of `det`: `det(X + Y) = det X + det Y + N(X,Y)`.
For a closed set `Δ ⊂ C'` let `K(Δ)` be the set of finite positive combinations `X = Σ λ_i p_i`, `p_i ∈ Δ`,
`λ_i > 0`, involving at least two distinct points; then `det X = Σ_{i<j} λ_iλ_j N(p_i,p_j) > 0`, and
`[X] := X/√det X` is a point of `AdS_3`. `hull(Δ) := {[X] : X ∈ K(Δ)}`. For two points `hull` is the geodesic
between them (§6.2(b)); for three it is the totally geodesic *ideal triangle* they span.

**Spacelike separation.** `x, y ∈ AdS_3 = PSL(2,R)` are *spacelike separated* if `½|tr(x^{-1}y)| > 1`, i.e. if
they are joined by a spacelike geodesic. For `X, Y` in cones `K(·)` over the same acausal circle,
`tr([X]^{-1}[Y]) = N(X,Y)/√(det X det Y) > 0`.

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

**HRT candidates.** A *candidate* for `A` is a quadruple `(C', Σ̃, γ, R)`: `C'` a slice through `∂A`; `Σ̃` a bulk
slice anchored on `C'`; `γ ⊂ Σ` a finite disjoint union of complete geodesics of `M`, each closed or having
both ends at points of `∂A`, with exactly one end at each point of `∂A`; and `R ⊂ Σ` a closed set, equal to the
closure of its interior, with frontier `γ`, such that the closure of its preimage in `Σ̃ ∪ C'` meets
`C' \ Λ_Γ` exactly in the lift of `A'`. (`γ` is an extremal curve anchored on `∂A`, and `R` is a spacelike
region with `∂R = γ ∪ A'`: the homology condition of Headrick–Takayanagi and of Wall's maximin.)

**Definition.** `S^ε(A) := (1/4G_N) · inf { L^ε(γ) : (C', Σ̃, γ, R) a candidate for A }`, with `L^ε` the regulated
length of §6.2.

*Remarks on the definition.* (i) (S3) is the strongest reasonable reading of "spacelike": it implies
achronality and acausality in `D̃`. Theorem H below holds verbatim for any class of slices between "(S1), (S2)
only" and "(S1)–(S3) and piecewise totally geodesic with spacelike faces", because its proof of `≥` uses only
(S1), (S2) and its proof of `≤` produces slices of the latter kind. (ii) If one insists on *smooth* spacelike
slices, the slices produced in §6.6 must be smoothed along their bending lines while keeping `γ` fixed. This
is elementary where the bending is locally finite (the tangent planes of a rounded crease interpolate
between two spacelike planes through a common spacelike line, with normals in the same time cone) but is
**not carried out here**, in particular not near closed components of `γ`, where bending lines accumulate.
(iii) Allowing non-embedded configurations with a mod-2 homology condition, as in arXiv:2407.07165, does not
change the infimum (Lemma 6.8).

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

Lemma 6.3' of the previous revision (two geodesics with unlinked endpoints) is the special case of Lemma 6.4
in which `X` and `Y` are combinations of two points each.

### 6.5 Lemma T (invariant ideal triangulations containing a given cut) **[proved modulo two cited standard facts]**

Let `c ∈ 𝒞(A)` with geodesic representative `c^X` on `X_L`; let `𝓕` be the set of lifts of its components to
`H²` (chords of `S¹`): axes of hyperbolic elements for the closed components, chords between lifted marked
points for the arcs.

**Lemma T.** There is a `Γ`-invariant set `𝒯` of ideal triangles of `H²`, with pairwise disjoint interiors,
such that
(T1) the edges of the triangles are pairwise non-crossing and none crosses a chord of `𝓕`;
(T2) every arc-chord of `𝓕` is an edge of two triangles; every axis `A ∈ 𝓕` is disjoint from the interiors of
all triangles;
(T3) the union of the triangles and of the axes in `𝓕` is `H²`;
(T4) `𝒯` is locally finite in `H²` minus the axes in `𝓕`, and every sequence of distinct triangles meeting a
fixed compact set converges, in the Hausdorff sense in the closed disc, to an axis of `𝓕`.

*Proof.* Add auxiliary marked points so that every ideal circle carries at least one, and let `V_Ω ⊂ Ω` be the
set of lifts of all marked points. On `X̄_L` the geodesics joining consecutive marked points of each ideal
circle ("boundary chords") cut off regions whose lifts are half-planes `H_e` over the elementary arcs
`e = [v, v']`; these have trivial stabiliser.

*Half-planes.* For each `Γ`-orbit of elementary arcs choose a homeomorphism `e ≅ [0,1]` and transport it by
`Γ`; take the triangles with vertices `k2^{-n}, (2k+1)2^{-n-1}, (k+1)2^{-n}`. They have disjoint interiors,
are locally finite in `H_e`, and their union is `H_e`: a point of `H_e` in no triangle would lie in the regions
cut off by a nested sequence of chords over dyadic intervals shrinking to a point of `e`, which is impossible.

*The compact part.* Cut `X̄_L` along `c^X` and the boundary chords, and discard the half-plane regions. Each
remaining piece `Y` is a compact surface whose boundary consists of closed geodesics of `c` ("smooth circles")
and of circles made of geodesic arcs meeting at ideal vertices. Collapse each smooth circle to a point (a
*puncture*) to get a compact surface `Y^♯` with marked points on its boundary and punctures. No `Y^♯` is a
sphere with at most two punctures, or a disc with at most two marked points and no puncture (two distinct
components of `c^X` or boundary chords are not homotopic, and none is null-homotopic).
**Standard fact 1** (Hatcher, *Topology Appl.* 40 (1991); Fomin–Shapiro–Thurston, *Acta Math.* 201 (2008),
§2): such a marked surface has a triangulation whose vertex set is the set of marked points and punctures and
whose edges are pairwise non-homotopic essential simple arcs, disjoint except at vertices, containing the
boundary arcs.
**Standard fact 2** ("spinning"; Thurston; Bonahon, *Ann. Fac. Sci. Toulouse* 5 (1996), §§1–3): choose, for each
smooth circle and each of its two sides, one of the two directions along it; replacing each edge by the
geodesic of `X_L` in its homotopy class which, at every end lying at a puncture, spirals towards the
corresponding closed geodesic in the chosen direction, gives pairwise disjoint simple geodesics whose
complementary regions in `Y` are the interiors of ideal triangles. In the universal cover: an edge ending at a
puncture lifts to a chord ending at the endpoint `e(A, side) ∈ {a^−, a^+}` of the relevant axis `A` selected
by the chosen direction; chords from one axis-side all share that endpoint; the triangles with a vertex at
the puncture lift to a `⟨h⟩`-orbit `(e(A), v_k, v_{k+1})`, `k ∈ Z`, of ideal triangles filling the region
between `A` and the chain of chords `(v_k, v_{k+1})`, with `v_k` tending to the other endpoint of `A` as
`k → −∞` and to `e(A)` as `k → +∞`.

The lifts of all these triangles, together with the half-plane triangles, form `𝒯`. (T1)–(T3) hold by
construction. (T4): away from the axes the family is the lift of a finite triangulation together with the
locally finite half-plane families; distinct triangles have disjoint interiors and equal hyperbolic area, so
a sequence of distinct triangles meeting a compact set cannot converge to a non-degenerate triangle, and its
Hausdorff limits are chords not crossing any edge, through a point of the compact set, which by (T3) are axes
of `𝓕`. ∎

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
non-degenerate simplices). For two distinct pieces, the equality case of Lemma 6.4 shows that their images
meet only in points of `hull{r,r'}` that are combinations of the cut points `r, r'` *in both pieces*, i.e. in
a common edge, where `Ψ` is consistently defined. Images of interior points do not lie on `p(C)` (they have
`det > 0`).

A continuous injective map from the compact disc `D̄` to a Hausdorff space is an embedding. So `Z ∪ C = Ψ(D̄)`
is an embedded closed disc with boundary `C` and interior `Z`. ∎

### 6.7 Theorem H (the HRT entropy is the minimum over cut classes) **[proved, given Lemma T]**

**Theorem H.** For the spacetimes of §6.1 and every region `A`,

`S^ε(A) = (1/4G_N) · min_{c ∈ 𝒞(A)} ½( ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c) )`,

and the infimum in the definition of `S^ε(A)` is attained.

*Proof.* `≥` (every candidate is a cut class). Let `(C', Σ̃, γ, R)` be a candidate. A closed component of `γ` in
the free homotopy class of `g ∈ π_1(S)` lifts to a spacelike geodesic of `D̃` invariant under
`(ρ_L(g), ρ_R(g))` and translated non-trivially by it; by §6.2(a) this is `hull{p(a^−), p(a^+)}`, the
geodesic joining the repelling and attracting fixed points, which lie in `Λ_Γ ⊂ C'`. An arc component lifts
to the geodesic joining two lifted points of `∂A`, also on `C'`. So the preimage `γ̃ ⊂ Σ̃` of `γ` is a
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
the spacetimes in §6.1 (cited); and the two standard facts in Lemma T (cited).

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

## 7. Status table (§§1–5 after three referee passes; §§6.3–6.9 rewritten in revision 5)

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
| Lemma T (invariant ideal triangulation containing a cut) | proved modulo two cited standard facts (triangulations of marked surfaces; spinning) |
| Proposition P (pleated slice: embedded, invariant, pairwise spacelike separated, contains the geodesics) | proved, given Lemma T |
| Theorem H (HRT entropy of §6.3 = minimum over cut classes) | proved, given Lemma T; **not yet refereed** when first written (see OPEN_GAPS.md for the review outcome) |
| Theorem 2 (HRT, contraction-map inequalities, arbitrary regions) | proved; no hypothesis H |
| Theorem 3' (HRT, full cone, whole-boundary parties) | proved; no hypothesis H |
| Variant of the definition with *smooth* spacelike slices | **not proved** (smoothing of the pleated slice not carried out) |
| Matter; higher dimensions; cusps; conical defects; non-orientable quotients; endpoints not on one slice | not covered |

**Standing assumptions** (collected at the referee's request): every `X_j` is orientable with funnel ends
only; the identifications `φ_j` preserve peripheral structure, marked points and the party assignment;
weights `w_j > 0` (positivity is essential — the `max` functional violates the inequalities); the same class
is used in all metrics; regions are unions of elementary pieces on one boundary Cauchy slice; `ε_0` depends
on the geometry and on the horoball normalisation; coefficients `α_l, β_r > 0`; `ρ_L, ρ_R` are holonomies of
the same marked oriented surface; the HRT entropy is the quantity defined in §6.3 (the leading classical
term; no quantum corrections; slices satisfying (S1)–(S3)); graph models are
allowed several terminals per party; citations in §1.1 and §6.1 were not re-checked against the sources.
