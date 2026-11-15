# PROOFS (revision 2, after first referee round)

Status legend: **[proved]** = complete argument written here, relying only on cited standard facts;
**[standard]** = textbook fact, cited; **[assumed]** = hypothesis, stated explicitly where used.

§§1–5 are mathematics about hyperbolic surfaces and contain no physics. Gravity enters only in §6, and the
link to the HRT formula goes through Hypothesis H (§6.3), which is **assumed**, not proved.

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
  Lebesgue measure zero (Beardon, *Amer. J. Math.* 93 (1971); Patterson, *Acta Math.* 136 (1976)).
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

`i(τ, c) = #(τ^X ∩ c^X)` (components of `c^X` equal to `τ^X` omitted).   (1.1)

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
diameter `≤ Y_0` at scale `ε_0`. Let `E_2'(ε) ⊃ E_2` be the set of geodesics meeting the ray `{x = 0, y > Y_0}`
and some other horoball at scale `ε`. A geodesic through `(0,h)`, `h > Y_0`, has endpoints `a < 0 < b` with
`|a| b = h²`; if it enters the region `{|x| ≤ 3Q_p, y ≤ Y_0}` below height... it must have an endpoint of
modulus `≤ 3Q_p + Y_0 =: Q'` (a semicircle whose both endpoints have modulus `> Q'` and which separates them
passes above every point of `{|x| ≤ 3Q_p, y ≤ Y_0}` it lies over). By symmetry take `|a| ≤ Q'`; then
`b > Y_0²/|a|`, and

`m(E_2'(ε_0)) ≤ 2∫_0^{Q'} d|a| ∫_{Y_0²/|a|}^{∞} db/(b + |a|)² = 2∫_0^{Q'} |a| d|a| /(Y_0² + a²) < ∞`.

The sets `E_2'(ε)` decrease as `ε ↓ 0` to a set of geodesics with an endpoint in `Γ·P̃` (Fact 3.2), which is
null. So `m(E_2) → 0`.

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

**Remarks.** (1) For `k = 1` this is the Ryu–Takayanagi statement on a hyperbolic surface. (2) Only two
properties of length are used: it is the integral of intersection numbers against a positive measure on test
classes, and its divergence at each end is universal.

---

## 5. Theorem 3 (membership in the full Ryu–Takayanagi cone) **[proved for `P = ∅`; argument for `P ≠ ∅` in §5.4, less detailed]**

This section removes the restriction to inequalities with contraction maps. The observation that the
functional is a graph min-cut for finite families of test curves is due to the first referee of this file.

### 5.1 Lemma G (graph realisation) **[proved]**

Let `μ = Σ_{k=1}^{m} w_k τ_k` be a finite positive combination of non-null-homotopic test classes, and
`F_μ(A) := min_{c ∈ 𝒞(A)} i(μ, c)`, `i(μ,c) = Σ_k w_k i(τ_k, c)`.

**Lemma G.** There is a finite weighted graph `G` with a labelling of some vertices by parties such that, for
every region `A`, `F_μ(A)` equals the minimum weight of an edge cut of `G` separating the vertices labelled by
parties in `A` from those labelled by parties not in `A` (purifier included). Hence `(F_μ(A))_A` is a
graph-model entropy vector and lies in the holographic entropy cone of Bao et al.

*Proof.* Represent the `τ_k` by geodesics (closed, or with ideal endpoints in the prescribed pieces, chosen
distinct), perturbed so that they are in general position (only transverse double points) and still realise
all pairwise and self intersection numbers minimally; let `T` be their union. `G` has as vertices the double
points of `T`, the ideal ends of the arcs (each labelled by the party of the piece containing it), and one
extra vertex on each closed `τ_k` without double points; its edges are the segments of `T` between
consecutive vertices, an edge on `τ_k` having weight `w_k`.

(≥) Let `c ∈ 𝒞(A)`; choose a representative transverse to `T`, missing the vertices, with
`#(c ∩ τ_k) = i(τ_k, c)` for all `k` (geodesic representative of `c`, slightly perturbed), with homology
region `W`. Let `S` be the set of vertices in `W`; an ideal end belongs to `S` iff its piece lies in `A`. Every
edge with exactly one endpoint in `S` is crossed by `c`, so the cut `(S, S^c)` has weight `≤ i(μ, c)`.

(≤) Let `S` be a vertex set containing the ideal ends in pieces of `A` and no other ideal end. Let `W` be the
union of a small closed regular neighbourhood of the vertices of `S` and of the edges with both endpoints in
`S`, together with a thin closed collar, in `X̄`, of every elementary piece of `A`, chosen so thin that it
meets `T` only near ideal ends lying in those pieces (which are in `S`). `W` is the closure of its interior,
its closure meets `∂_∞X` exactly in `A`, and `∂W` is an embedded multi-curve whose arcs end exactly at the
points of `∂A`: a cut for `A`. It crosses each edge between `S` and `S^c` once and no other edge, so
`F_μ(A) ≤ weight(S, S^c)`. ∎

### 5.2 Filling

**Lemma F.** There is a finite set `T_0` of test classes of proper arcs, each containing an open set of
geodesics (hence of positive `m̄`-measure and, for small `ε`, of positive `m̄|_{𝒢_ε}`-measure), such that for
every `M` the set `{c ∈ ⋃_A 𝒞(A) : Σ_{τ ∈ T_0} i(τ, c) ≤ M}` is finite.

*Proof.* Choose finitely many pairwise disjoint simple proper geodesic arcs with ideal endpoints in
`∂_∞X \ P` cutting `X̄` into closed discs, each of which meets `∂_∞X` in a union of intervals (possible since
every end of `X` is a funnel: take an ideal triangulation of the surface with boundary). Moving the ideal
endpoints slightly keeps the class, so each class contains an open set of geodesics. A cut in minimal
position with respect to these arcs meets each disc in arcs joining sides of the disc or marked points on
its ideal sides, without bigons; its class is determined by the number of arcs of each of the finitely many
combinatorial types in each disc, and these numbers are bounded in terms of the intersection numbers. ∎

### 5.3 The case `P = ∅` (parties are unions of whole boundary circles) **[proved]**

**Theorem 3 (`P = ∅`).** With the notation of Theorem 1 and `P = ∅`, the vector `(F(A))_A`,
`F(A) = min_c Σ_j w_j ℓ_{X_j}(c)`, lies in the holographic entropy cone; i.e. it satisfies every linear
inequality valid for Ryu–Takayanagi entropies (equivalently, for all graph models).

*Proof.* By Crofton and (B0), `ℓ_{X_j}(c) = ∫ i(τ(ĝ), c) dm̄_j(ĝ)`. Test classes form a countable set, so the
push-forward of `Σ_j w_j m̄_j` (transported to `X_1`) under `ĝ ↦ τ(ĝ)` is a countable sum `μ = Σ_τ μ_τ τ` with
`i(μ, c) = Σ_j w_j ℓ_{X_j}(c) < ∞` for every `c`. Let `μ_n` be finite partial sums increasing to `μ`, containing
`T_0` from the start. By Lemma G, `(F_{μ_n}(A))_A` lies in the cone. `F_{μ_n}(A)` is non-decreasing in `n` and
bounded by `F(A)`. A minimiser `c_n` of `F_{μ_n}(A)` satisfies `Σ_{τ∈T_0} μ_τ i(τ,c_n) ≤ F(A)`, so by Lemma F
the `c_n` range over a finite set `𝒟`; then
`lim_n F_{μ_n}(A) = lim_n min_{c∈𝒟} i(μ_n,c) = min_{c ∈ 𝒟} i(μ,c) ≥ F(A)`. Hence `F_{μ_n}(A) → F(A)` for every
`A`, and the cone is closed. ∎

### 5.4 The case `P ≠ ∅` **[argument given; the bookkeeping is less detailed than §§2–4]**

**Theorem 3 (`P ≠ ∅`).** For every linear inequality `Q(S) = Σ_l α_l S(A_l) − Σ_r β_r S(B_r) ≥ 0` valid on the
holographic entropy cone: if `Q` is balanced at every marked point then `Q(F^ε) ≥ 0` for all `ε ≤ ε_0`;
otherwise `Q(F^ε) → +∞`.

*Argument.* Validity on the cone gives `Δ_p ≥ 0` at each `p` (apply `Q` to the graph with one unit edge
joining the two parties adjacent at `p`, all other parties isolated). By (4.1) it remains to show
`Q(F^ε) ≥ −o(1)`. Fix `ε`, let `γ_l` be the (`ε`-independent) minimisers, `ν_ε` the push-forward of
`Σ_j w_j m̄_j|_{𝒢_ε(X_j)}` to test classes, and `μ_n ↑ ν_ε` finite partial sums containing `T_0`. By Lemma G,
`Σ_r β_r F_{μ_n}(B_r) ≤ Σ_l α_l F_{μ_n}(A_l) ≤ Σ_l α_l i(ν_ε, γ_l)`. Let `δ_r^{n}` be minimisers of
`F_{μ_n}(B_r)`. Split `ν_ε = ν^{far} + ν^{near}`, `ν^{far}` = the part coming from `𝒢_{ε_0}`. Two uniform facts:

* for every class `c` with an arc-end at `p`, the `ν^{near}`-mass crossing its terminal ray between
  `∂B_p(ε_0)` and `∂B_p(ε)` is at least `(Σ_j w_j)(log(ε_0/ε) − κ − θ(ε))` with `θ(ε) → 0` independent of `c`
  (this is `log(ε_0/ε) − m(E_1) − m(E_2)` of Lemma B; the bounds there depend only on `Q_p` and `ε_0`);
* `i(ν^{near}, c_ref) ≤ (Σ_j w_j)|∂B| log(ε_0/ε) + O(1)` for a fixed reference class `c_ref ∈ 𝒞(B)`.

Comparing `δ_r^n` with `c_ref` (both have `|∂B_r|` ends) these give
`i(μ_n^{far}, δ_r^n) ≤ O(1) + (ν_ε − μ_n)`-tail terms, uniformly in `ε`; for `n` large the tail is `≤ 1`. By
Lemma F the `δ_r^n` then lie in a finite set independent of `n` and `ε`. Passing to a subsequence with
constant `δ_r` and letting `n → ∞`, `Σ_r β_r i(ν_ε, δ_r) ≤ Σ_l α_l i(ν_ε, γ_l)`, and Lemma B, now applied to
finitely many classes, gives `Q(F^ε) ≥ κ(Σ_j w_j)Δ − o(1)` as in Theorem 1. ∎

---

## 6. Application to vacuum 2+1 gravity

### 6.1 The spacetimes **[standard]**

`AdS_3` = `SL(2,R)` with `ds² = −det(dx)`; `SL(2,R) × SL(2,R)` acts by `x ↦ g_L x g_R^{-1}`.

Let `S` be an oriented surface of finite type with non-empty boundary and negative Euler characteristic, or a
disc or an annulus, and `ρ_L, ρ_R : π_1(S) → SL(2,R)` two representations whose images are groups as in §1.1
and which are holonomies of hyperbolic structures `X_L, X_R` on the interior of `S` **in the same component
of Teichmüller space** (i.e. `X_R` is obtained from `X_L` by an orientation-preserving homeomorphism isotopic
to the identity relative to the marking; boundary elements go to boundary elements). Then
`Γ = {(ρ_L(γ), ρ_R(γ))}` acts properly discontinuously on a convex domain `D̃ ⊂ AdS_3`; `M = D̃/Γ` is a
locally-`AdS_3` spacetime with one asymptotically-AdS region per end of `S` (Mess 1990; Barbot, *Causal
properties of AdS-isometry groups I, II*; Bonsante–Krasnov–Schlenker, arXiv:math/0610429). Physically:
rotating BTZ (`π_1 = Z`), rotating multi-boundary wormholes, wormholes with handles
(Aminneborg–Bengtsson–Holst 1998). `ρ_L = ρ_R` is the time-symmetric case; for `ρ_L ≠ ρ_R` the spacetime has
no time-reflection symmetry, and for non-abelian `π_1` it is not stationary.

**Scope.** Only these quotients are covered: orientable, all holonomies hyperbolic (no extremal horizons,
no cusps, no conical defects).

The `i`-th asymptotic boundary is `(I^i_L × I^i_R)/⟨h_i⟩`, `I^i_{L,R}` being the components of `Ω_{L,R}`
stabilised by the peripheral element `h_i`. A *boundary Cauchy slice* is, in each asymptotic boundary, a
closed acausal curve; its lift to `I^i_L × I^i_R` is the graph of an increasing homeomorphism
`ψ_i : I^i_L → I^i_R` commuting with `h_i`.

**Marking.** Let `φ̃ : S¹ → S¹` be the `ρ_R ∘ ρ_L^{-1}`-equivariant homeomorphism equal to the boundary map of
the marking on `Λ_L` and to the translates of the `ψ_i` on `Ω_L`. It descends to
`φ : ∂_∞X_L → ∂_∞X_R` and extends to a homeomorphism `X̄_L → X̄_R` in the isotopy class of the marking. A
finite set of points of the boundary slice gives `P_L` and `P_R = φ(P_L)`; regions of the slice are regions
(§1.1) for both surfaces; and for boundary points `p, q` of the slice with lifts `(ũ_p, ṽ_p)`, `(ũ_q, ṽ_q)`,
`ṽ = φ̃(ũ)`, the arc class on `X_L` of the geodesic from `ũ_p` to `ρ_L(γ)ũ_q` corresponds under `φ` to the arc
class on `X_R` of the geodesic from `ṽ_p` to `ρ_R(γ)ṽ_q`. This fixes the identification of cut classes on the
two surfaces, including boundary twists.

### 6.2 Proposition (length formula) **[proved]**

Boundary points of `AdS_3` are rank-one matrices up to positive scale, `p = a b^T`; put `u = [a]`,
`v = [Jb]` with `J = [[0,1],[−1,0]]`, so that `g_L p g_R^{-1}` corresponds to `(g_L u, g_R v)`.
`ω(a,c) := a^T J c`.

(a) For `γ ≠ 1`, `(ρ_L(γ), ρ_R(γ))` preserves exactly two geodesics of `AdS_3`, both spacelike, with
translation lengths `½(ℓ_L(γ) + ℓ_R(γ))` and `½|ℓ_L(γ) − ℓ_R(γ)|`. The first joins the pair of attracting
fixed points to the pair of repelling fixed points, both of which lie on the limit curve of `Γ`, and is the
one contained in `D̃` (the convex hull of the limit curve lies in `D̃`); it is the outer horizon in BTZ, the
second the inner horizon.

*Proof.* Conjugate to `ρ_L(γ) = ±diag(e^{a/2}, e^{−a/2})`, `ρ_R(γ) = ±diag(e^{b/2}, e^{−b/2})`, `a = ℓ_L`,
`b = ℓ_R > 0`. Then `x ↦ g_L x g_R^{-1}` multiplies the entries `x_{11}, x_{12}, x_{21}, x_{22}` by
`e^{(a−b)/2}, e^{(a+b)/2}, e^{−(a+b)/2}, e^{(b−a)/2}` (up to a common sign). The curves
`s ↦ [[0, e^{s}],[−e^{−s}, 0]]` and `s ↦ [[e^{s}, 0],[0, e^{−s}]]` are unit-speed spacelike geodesics
(`−det(x') = 1`), translated by `(a+b)/2` and `(a−b)/2`. An invariant geodesic has invariant endpoints on
the boundary, i.e. rank-one eigen-matrices, which are the four matrix units; the two geodesics above are the
only ones joining pairs of them inside `AdS_3`. ∎

(b) Let `p = a b^T`, `q = c d^T` with `N := ω(a,c)·ω(b,d) > 0`. Then `x(t) = (e^{t} p + e^{−t} q)/√N`,
`t ∈ R`, is the unit-speed geodesic of `AdS_3` from `q` to `p`, and the length of its part where the
coefficient of `p` and that of `q` are both `≤ 1/ε` is `log N + 2 log(1/ε)`.

*Proof.* `det(s·ab^T + t·cd^T) = st·ω(a,c)ω(b,d)` (expand the determinant of a sum of two rank-one matrices).
So `det x(t) = 1`, `−det(x'(t)) = 1`, and `x'' = x` (geodesic of the quadric). The truncation is
`e^{|t|}/√N ≤ 1/ε`, of length `2 log(√N/ε)`. ∎

For `N > 0` the points are spacelike related in the sense needed; for points of a boundary Cauchy slice and
their `Γ`-translates this holds after fixing the signs of the vectors along the lifted slice (the slice
together with the limit curve is an achronal circle).

In `H²`, with horocycles at `[a]`, `[c]` determined by the vectors `a, c`, the signed distance between the
horocycles is `2 log|ω(a,c)|` (Penner, *Comm. Math. Phys.* 113 (1987)); the square root of its exponential is
the lambda length). Hence, for every class `c` of cuts and compatible cut-offs — on `X_L` the horocycle at
`u_p` determined by `a_p`, on `X_R` that at `v_p` determined by `b_p`, where the spacetime cut-off is the scale
of `p = a_p b_p^T` —

`L^ε(c) = ½( ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c) )`,   (6.1)

where `L^ε(c)` is the total regulated length of the spacetime geodesics of class `c`: for a closed component in
the conjugacy class of `γ`, the geodesic of (a) in `D̃`; for an arc from `p` to `q` labelled by `γ`, the
geodesic of (b) from `(ũ_p,ṽ_p)` to `(ρ_L(γ)ũ_q, ρ_R(γ)ṽ_q)`. By Fact 3.1 horoball removal and horocycle
distance agree for `ε ≤ ε_0`. Rescaling `p → p/ε` shrinks both horoballs. Equation (6.1) is equations (58),
(60) of arXiv:2407.07165 in hyperbolic-surface language.

### 6.3 Hypothesis H **[assumed]**

> **H.** In the spacetime `M` of §6.1, for a boundary region `A` made of intervals of a boundary Cauchy slice
> and of whole boundary components, the HRT entropy is
> `S^ε(A) = (1/4G_N)·min_{c ∈ 𝒞(A)} L^ε(c)`.

What H contains, item by item:

* **H1 (candidates are embedded).** Every extremal surface that is spacelike-homologous to `A` is the
  geodesic realisation of some `c ∈ 𝒞(A)`. Reason: its homology region is a spacelike surface with boundary,
  contained in a Cauchy surface homeomorphic to `S`, so the surface is an embedded cut. *Not written out as a
  proof here.*
* **H2 (every class is admissible).** For every `c ∈ 𝒞(A)` the spacetime geodesics of class `c` exist in `M`,
  are pairwise spacelike separated, and bound a spacelike region together with `A`. Existence follows from
  (6.2)(a),(b) and convexity of `D̃`. Pairwise spacelike separation and the existence of the spacelike region
  are **not proved here**. (The referee reports a computation indicating that two boundary-anchored geodesics
  with unlinked endpoints on an acausal curve are spacelike separated; I have not verified it and do not rely
  on it.)
* **H3 (the literature's prescription).** arXiv:2407.07165 minimises over collections of group elements with
  a homology condition checked on abelianised words, which allows non-embedded collections. Lemma 6.4 below
  shows these never beat embedded ones, so the two minima agree.

If H2 failed for some classes, Theorems 2 and 3' would still be true statements about the function
`min_c L^ε(c)`, but that function could be smaller than the HRT entropy for some regions, and the conclusion
about HRT would not follow.

### 6.4 Lemma (non-embedded collections do not win) **[proved]**

Let `c'` be a finite collection of closed curves and arcs between marked points, transverse, whose union has
the mod-2 intersection property of a cut for `A` (i.e. `X \ c'` admits a 2-colouring changing across each
branch of `c'`, with the pieces of `A` on colour 1 and exactly one arc-end at each point of `∂A`). Then there
is `c ∈ 𝒞(A)` with `i(τ, c) ≤ i(τ, c')` for all test classes `τ`, hence `ℓ^ε_{X}(c) ≤ ℓ^ε_X(c') + o(1)` for
every hyperbolic structure `X`; with (3.1), `ℓ^ε_X(c) ≤ ℓ^ε_X(c')` for all `ε ≤ ε_0` when `c'` has the same
ends.

*Proof.* Put `c'` in geodesic position; let `W` be the closure of the colour-1 set and resolve the crossings
as in Lemma A (vertices and marked points). The result is a cut for `A` made of pieces of `c'`, so any
`τ'` realising `i(τ, c')` meets it in at most `#(τ' ∩ c')` points. Integrate as in Theorem 1. ∎

### 6.5 Theorems 2 and 3' (covariant entropy inequalities in vacuum 2+1 gravity) **[proved given H]**

Let `M`, the boundary slice and the parties be as in §6.1, and `S^ε` as in H. By (6.1),
`4G_N S^ε(A) = F^ε(A)` for `k = 2`, `w_L = w_R = ½`.

**Theorem 2.** Every inequality `Σ α_l S(A_l) ≥ Σ β_r S(B_r)` admitting a contraction map holds for `S^ε`:
exactly, for every `ε ≤ ε_0`, if it is balanced at every endpoint of the regions (in particular when all
parties are unions of whole boundary components); and with the difference of the two sides tending to `+∞`
otherwise.

**Theorem 3'.** The same for every inequality valid on the holographic entropy cone (proved for whole-boundary
parties, §5.3; argument of §5.4 for intervals). Consequently, given H, the HRT entropy vectors of these
vacuum spacetimes lie in the Ryu–Takayanagi holographic entropy cone.

*Proof.* Theorems 1 and 3. ∎

---

## 7. Status table

| Statement | Status |
|---|---|
| Lemma A (topological contraction lemma) | proved; mechanism due to Czech–Dong |
| Lemma B (Crofton with horoballs, `κ = 1 − log 2`) | proved; constant and normalisation re-derived numerically by the referee |
| Theorem 1 (multi-metric RT, contraction-map inequalities) | proved |
| Lemma G (graph realisation), Lemma F (filling) | proved (Lemma G suggested by the referee) |
| Theorem 3, `P = ∅` (full cone, whole-boundary parties) | proved |
| Theorem 3, `P ≠ ∅` (full cone, intervals) | argument given, bookkeeping not at the level of §§2–4 |
| Length formula (6.1) and marking | proved |
| Lemma 6.4 | proved |
| Hypothesis H: H1 | reason given, not written as a proof |
| Hypothesis H: H2 (admissibility of every embedded class) | **assumed** |
| Theorems 2, 3' (statements about HRT) | proved **given H** |
| Matter; higher dimensions; cusps; conical defects | not covered |
