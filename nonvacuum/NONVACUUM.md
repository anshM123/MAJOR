# NONVACUUM: contraction-map inequalities for HRT in 2+1 with matter (NEC)

Working notes. Nothing outside this directory is edited. Builds on the vacuum proof (PROOFS.md), which is
taken as correct; only its purely topological / measure-theoretic parts are reused: Lemma A (§2), the
integration step of Theorem 1 (§4), and Lemma G / Lemma F / Theorem 3 (§5).

Grades: **[proved]** complete argument here; **[proved mod X]** complete except for the cited/assumed X;
**[conj + evidence]**; **[false]**.

## 0. Progress log
* Report 1 (start). Plan fixed; the route differs from OPEN_GAPS G2 in one important way (see §1): the
  analogue of Theorem H is obtained from maximin + Wall's null projection, *without* pleated slices and
  without needing a common slice for the different component geodesics.

## 1. Setting, definitions, hypotheses

**Spacetime.** `M` is a smooth, time- and space-orientable, asymptotically AdS_3 spacetime whose conformal
compactification `M̄ = M ∪ ∂M` is globally hyperbolic as a manifold with timelike boundary, with Cauchy
slices `Σ̄` that are compact orientable surfaces with `b ≥ 1` boundary circles (one per asymptotic region).
Any genus, any `b`, rotating, time dependent. `M` satisfies the null curvature condition (NEC)
`R_{μν}k^μk^ν ≥ 0` for null `k` (Einstein's equations with `T_kk ≥ 0`). `Γ := π_1(Σ̄)` (free; trivial for a
disc, `Z` for an annulus). Excluded: interior ends (cusps / extremal horizons), conical defects,
non-orientable slices.

**Boundary slice, admissible slices.** Fix a boundary Cauchy slice `C ⊂ ∂M` (an acausal circle in each
boundary cylinder). `𝒮(C)` = acausal Cauchy slices `σ̄` of `M̄` with `σ̄ ∩ ∂M = C` (the
Headrick–Hubeny–Lawrence–Rangamani convention; the same one used in PROOFS §6.3b). Assume `𝒮(C) ≠ ∅`. All
`σ̄ ∈ 𝒮(C)` are identified with one another (and with an abstract `Σ̄`) by the flow of a timelike vector
field tangent to `∂M`, which is the identity on `C`; the inclusion `σ̄ ⊂ M̄` is a homotopy equivalence. Hence
homotopy classes do not depend on the slice:
* a *closed class* is a free homotopy class of loops in `M̄` (↔ conjugacy class in `Γ`);
* an *arc class* `κ` between `p, q ∈ C` is a homotopy class rel endpoints of paths in `M̄` from `p` to `q`.

**Lengths.** Fix a defining function `z` (Fefferman–Graham type) near `∂M`. For a curve with `e` ends on `∂M`,
`ℓ^ren := lim_{ε→0} [ℓ(curve ∩ {z > ε}) − e·log(1/ε)]`. For `σ ∈ 𝒮(C)` and a class `κ`,

`ℓ_σ(κ) := inf { ℓ^ren(η) : η ⊂ σ̄ a rectifiable curve in κ }`,   `L(κ) := sup_{σ ∈ 𝒮(C)} ℓ_σ(κ)`.

`L(κ)` is the *class-restricted maximin length*. For a multi-class `c = {c_1,…,c_m}` (e.g. a cut class) put
`L(c) := Σ_k L(c_k)`.

**Cut classes, regions, test classes.** As in PROOFS §1 (topological notions on `Σ̄`, with marked points
`P ⊂ C`, regions = unions of elementary pieces of `C`, `𝒞(A)` = cut classes for `A`). To use PROOFS §§1–5
we fix an *auxiliary* hyperbolic structure `X_0` with funnel ends on `int Σ̄` (for the disc: `H²`; annulus:
hyperbolic cylinder). It plays no geometric role: it only supplies an equivariant homeomorphism of the
ideal circle `∂_∞Σ̃ := C̃ ⊔ ∂Γ` (lifted boundary lines plus Gromov boundary, with its cyclic order) onto
`S¹ = ∂H²`, carrying `C̃` onto `Ω_0 = S¹ ∖ Λ_0`, and geodesic representatives for Lemma A. `Λ_0` is nowhere
dense, so `C̃` is dense in `S¹`.

**HRT entropy.** A *candidate* for `A` is `(σ, γ, R)` with `σ ∈ 𝒮(C)`, `γ ⊂ σ` a finite disjoint union of
simple complete spacelike geodesics of `M` (extremal curves), each closed or with both ends in `∂A`, one end
at each point of `∂A`, and `R ⊂ σ̄` closed, the closure of its interior, with frontier `γ` and `R ∩ C = A`.
`S(A) := (1/4G_N)·inf ℓ^ren(γ)` (below, `4G_N = 1`). This is PROOFS §6.3 with the matter spacetime in place
of `D̃/Γ` and the boundary slice held fixed.

**Maximin quantities.** `min_σ(A) := inf` of renormalized perimeters of closed sets `R ⊂ σ̄` (finite
perimeter, closure of interior) with `R ∩ C = A`; equivalently inf of `ℓ^ren(∂R)` over regions bounded by
rectifiable curves.

### Hypotheses (all standard in the covariant holographic literature; none is proved here)

* **(MM) Maximin.** (a) For every region `A`, there is a candidate `(σ_*, γ_*, R_*)` with
  `ℓ^ren(γ_*) = min_{σ_*}(A) = sup_{σ∈𝒮(C)} min_σ(A)`. (b) For every single class `κ` (arc or closed), the
  sup defining `L(κ)` is attained by a slice `σ_κ` and an `ℓ_{σ_κ}`-minimiser `m_κ` that is a simple extremal
  geodesic; in particular `L(κ) < ∞`.
  (a) is Wall's maximin theorem (arXiv:1211.3494; existence is proved there for spacetimes without horizons
  and for Kasner-like singularities, and extended by Marolf–Wall–Wang, arXiv:1901.03879, "restricted
  maximin" anchored on a boundary Cauchy slice — exactly our `𝒮(C)` setting). (b) is the same statement with
  the minimisation restricted to one homotopy class; Wall's deformation argument is local near the
  minimiser and is insensitive to the restriction. **Only the consequence `S(A) ≤ sup_σ min_σ(A)` of (a) is
  used.**
* **(MM2) Joint maximin for a nested pair.** For two arc classes `κ_1, κ_2` with lifts `ã→d̃` and `b̃→c̃`,
  `a<b<c<d` cyclically on `∂_∞Σ̃` (the second nested inside the first),
  `sup_σ [ℓ_σ(κ_1) + ℓ_σ(κ_2)] = L(κ_1) + L(κ_2)`.
  This is the class-restricted form of the common-slice property for nested regions that drives Wall's
  proof of strong subadditivity (and MWW's restricted maximin). Sketch why it should hold: maximise
  `F(σ) = ℓ_σ(κ_1)+ℓ_σ(κ_2)`; at a maximiser the two minimisers are `σ`-geodesics in different classes, so they
  meet in isolated transverse points; Wall's local deformation at a non-extremal point of one of them, away
  from the other, would increase `F`; so both are extremal, and each extremal achronal curve has length
  `≥ L` of its class by Lemma P below. The gap is the same as in Wall's single-surface argument
  (non-unique minimisers / stability), no worse.
* **(R) Regularity / asymptotics.** (R1) the null-hypersurface regularity needed for Wall's projection
  (area theorem for achronal null boundaries, Chruściel–Delay–Galloway–Howard), used in Lemma P;
  (R2) bounded families of extremal geodesics with endpoints in a compact part of `C` are precompact, with
  convergence of renormalized lengths; (R3) every slice that is relevant (one with `ℓ_σ(κ) > −∞` for arc
  classes) is asymptotically hyperbolic near `C` with constants depending on `σ` only, so that two curves on
  `σ` ending at the same boundary point can be joined at the cut-off at bounded cost.

## 2. The HRT entropy is the minimum over cut classes of class-maximin lengths (Theorem H')

### Lemma P (componentwise null projection) **[proved mod Wall's projection lemma + (R1)]**

Let `g` be a complete spacelike geodesic of `M`, closed or with both ends on `C`, lying on some
`σ_0 ∈ 𝒮(C)`, and let `κ = [g]`. Under NEC, `ℓ_σ(κ) ≤ ℓ^ren(g)` for every `σ ∈ 𝒮(C)`. Hence `L(κ) ≤ ℓ^ren(g)`.

*Proof.* Let `H = 1` if `g` is an arc and `H = ⟨h⟩` if `g` is closed in class `h`; let `M̂ = M̃/H` (a cover of `M`,
globally hyperbolic with timelike boundary), `ĝ` a lift of `g`, `σ̂_0, σ̂` the lifts of the slices. `ĝ` is a
properly embedded line (arc case: from `â` to `b̂` on the lifted boundary) or an embedded closed curve
(closed case, `σ̂_0` an annulus) separating `σ̄̂_0` into two sides. Let `r ⊂ σ̄̂_0` be the closure of one side
("left"); its boundary-at-infinity part is the left arc `C_L` of the lifted boundary slice, and its edge is
`ĝ`. Because `ĝ` is extremal both null expansions of `ĝ` vanish; by NEC and Raychaudhuri
(`dθ/dλ ≤ −θ²` for a one-dimensional cross-section, no shear in 2+1) the outward null congruences from `ĝ`
have `θ ≤ 0` until they leave the achronal boundary. This is exactly the hypothesis of Wall's projection
lemma (the "representative" construction of arXiv:1211.3494 used to prove that the maximin surface is the
minimal extremal surface): the set `q := σ̄̂ ∩ (J^+(D(r)) ∪ J^-(D(r)))` is closed in `σ̄̂`, and its frontier in
`σ̂` lies on the outward null sheets of `ĝ`, has `H^1`-measure `≤ ℓ(ĝ)` (lengths measured on the null sheet,
whose induced metric is `ψ² ds²` with `ψ ≤ 1`).

Two topological facts finish the proof. (i) `q ∩ C̃ = C_L`: if `c ∈ C_R` were in `J^±(D(r))`, extending the
causal curve inextendibly through `D(r)` would produce a causal curve between `c ∈ σ̄̂_0` and a point of
`r ⊂ σ̄̂_0`, contradicting acausality of `σ̄_0`. (ii) The frontier separates `C_L` from `C_R` in `σ̄̂`: it is the
frontier of a closed set containing `C_L` and disjoint from `C_R`. A closed set of finite `H^1`-measure
separating the two boundary arcs of a closed disc (resp. the two ends of an annulus) contains a simple arc
from `â` to `b̂` (resp. a simple closed curve generating `π_1`) [plane topology: a separating closed set
contains a separating continuum (Janiszewski); a continuum of finite length is arcwise connected]. That
curve lies in class `κ` after projection to `M` and has length `≤ H^1(frontier) ≤ ℓ(ĝ)`. Renormalization:
near the ends the outward null generators from points of `ĝ` at `z ≈ ε` move `z` only by `O(z²)` (true in
exact AdS_3, where null normals to a radial geodesic are tangent to `z = const`; in general by the
asymptotic expansion), so the cut-off contributes `o(1)` — this is part of (R). ∎

*Why the cover.* Wall's lemma is applied not to the homology region of a whole HRT configuration but to one
side of a single component in the cover where that component is the full edge. This is what makes the
statement *componentwise* and *class by class*, which is what the contraction argument needs.

### Lemma Q (mod-2 resolution on one slice) **[proved]**

For `σ ∈ 𝒮(C)`, a region `B` and any `δ = {δ_1,…,δ_m} ∈ 𝒞(B)`: `min_σ(B) ≤ Σ_k ℓ_σ(δ_k)`.

*Proof.* Take near-minimising curves `η_k ⊂ σ̄` in the classes `δ_k` (they may cross each other and
themselves). Let `δ^emb` be an embedded representative on `σ̄ ≅ Σ̄` with homology region `W`. A homotopy
(rel endpoints, pushed off `C` except at the endpoints) from `δ^emb_k` to `η_k` gives a mod-2 2-chain `H_k`
with `∂H_k = η_k + δ^emb_k` (mod 2). The mod-2 flat 2-chain `R := W + Σ_k H_k` is (the indicator of) a set of
finite perimeter with `∂R = Σ_k η_k + (B-part of C)` mod 2, so its perimeter is at most the mass of `Σ η_k`,
i.e. `≤ Σ_k ℓ(η_k)`, and `R ∩ C = B`. (This is PROOFS Lemma 6.8 on a Riemannian, possibly Lipschitz,
surface, with resolution of crossings replaced by mod-2 chains.) ∎

### Theorem H' **[proved mod (MM)(a), Lemma P]**

For every region `A`: `S(A) = min_{c ∈ 𝒞(A)} L(c)`, the minimum is attained by the class `c_A` of any HRT
configuration, and for that class `L(c_A) = Σ_k ℓ^ren(g_k)` (each HRT component has the class-maximin length
of its class). Equivalently, with `ℓ_σ(c) := Σ_k ℓ_σ(c_k)`,

`sup_σ min_{c} ℓ_σ(c)  =  S(A)  =  min_c Σ_k sup_σ ℓ_σ(c_k)`   (a minimax theorem; the saddle is the HRT
configuration on its maximin slice).

*Proof.* `≥`: let `(σ_A, γ, R)` be any candidate, components `g_k`. On `σ̄_A ≅ Σ̄`, `γ` is an embedded
multi-curve bounding `R` with `R ∩ C = A`, so its class `c` lies in `𝒞(A)` (delete null-homotopic and
repeated closed components; their lengths are `≥ 0`). Each `g_k` lies on `σ_A`, so Lemma P gives
`L(c_k) ≤ ℓ^ren(g_k)`. Summing, `ℓ^ren(γ) ≥ L(c) ≥ min_c L(c)`; take the inf over candidates.
`≤`: for any `c ∈ 𝒞(A)`, by (MM)(a) and Lemma Q,
`S(A) ≤ sup_σ min_σ(A) ≤ sup_σ Σ_k ℓ_σ(c_k) ≤ Σ_k L(c_k)`. The HRT configuration exists by (MM)(a); its class
attains the min and the chain of inequalities is then an equality. ∎

**Remarks.** (1) This replaces the whole pleated-slice machinery of PROOFS §§6.4–6.7. The common-slice
problem for the *right-hand side* classes `δ_r` of a contraction argument disappears: `S(B) ≤ L(δ)` holds
for every class `δ`, even when the extremal geodesics of its components cross each other or lie on no
common slice, because maximin lets each slice be tested separately. (2) Non-uniqueness of extremal
geodesics in a class, and extremal geodesics in a class that are not minimal, are irrelevant: `L(κ)` is
the class maximin, which by Lemma P and (MM)(b) is the minimum length of an extremal geodesic of class `κ`
lying on an acausal slice anchored at `C`. (3) Horizons and entanglement shadows play no role: classes
are topological, and closed extremal geodesics behind horizons are just closed classes.
(4) Vacuum check: for `M = D̃/Γ`, every class has a unique extremal geodesic, so `L(c) = ½(ℓ_L(c)+ℓ_R(c))` and
Theorem H' reproduces PROOFS Theorem H (with fixed boundary slice) without pleated slices — but at the
price of (MM), which in vacuum the pleated slice proves by hand.

## 3. The function D on the ideal circle and the crossing inequality

For distinct `x, y ∈ C̃` let `κ_{x,y}` be the arc class in `M` of the projection of a path in `σ̄̃` from `x` to `y`
(well defined: `σ̄̃` is simply connected, and slices are identified by the flow), and

`D(x, y) := L(κ_{x,y})`.

So `D(x,y)` is the class-maximin renormalized length; by Lemma P and (MM)(b) it is the least length of an
extremal geodesic of `M̃` from `x` to `y` lying on (the lift of) an acausal slice anchored at `C`. This is the
`D` of OPEN_GAPS G2, made well defined between points on *different* lifted boundary lines: it is a
quantity on the quotient (one homotopy class), never an HRT quantity of the non-compact cover.
`D(γx, γy) = D(x,y)` for `γ ∈ Γ`.

### Lemma X (crossing inequality) **[proved mod (MM2)]**

For `a < b < c < d` cyclically ordered in `C̃`: `D(a,c) + D(b,d) ≥ D(a,d) + D(b,c)`.

*Proof.* Fix `σ`. Let `η_{ac}, η_{bd}` be curves in `σ̄̃` with these endpoints and near-minimal length. Their
endpoints are linked, so they meet at a point `x` (Jordan, in the closed disc `σ̄̃ ∪ ∂_∞`). Exchanging the tails
at `x` gives curves from `a` to `d` and from `b` to `c` with the same total length. Projecting to `σ̄`:
`ℓ_σ(κ_{ad}) + ℓ_σ(κ_{bc}) ≤ ℓ_σ(κ_{ac}) + ℓ_σ(κ_{bd}) ≤ L(κ_{ac}) + L(κ_{bd})`. Take `sup_σ` of the left side and
use (MM2) for the nested pair `(ad, bc)`. ∎

*Comment.* This is strong subadditivity on the universal cover, `S(AB)+S(BC) ≥ S(B)+S(ABC)` for the
"intervals" `[a,c], [b,d], [b,c], [a,d]` of `∂_∞Σ̃`, but formulated on the quotient so that no maximin theory on
the non-compact cover (infinitely many boundaries, Cantor set of ends) is needed. NEC enters through (MM2)
(the extremality of the joint maximiser and Lemma P). Without NEC it is false in general: see §7 (Vaidya
test, NEC-violating sign).

### Lemma C (continuity) **[proved mod (MM)(b), (R)]**

`D` is continuous on `{(x,y) ∈ C̃² : x ≠ y}`.

*Proof.* Lower semicontinuity: `D = sup_σ ℓ_σ`, and for each relevant `σ` (R3) `ℓ_σ(κ_{x,y})` is continuous in
the endpoints. Upper semicontinuity: if `x_n → x`, let `m_n` be the class-maximin extremal geodesics of
(MM)(b), `ℓ^ren(m_n) = D(x_n, y)`. By (R2) a subsequence converges to an extremal geodesic `m` of class
`κ_{x,y}` with `ℓ^ren(m_n) → ℓ^ren(m)`; `m` lies on a limit of the slices `σ_{κ_n}`, an achronal slice anchored
at `C`; Lemma P (whose proof only uses achronality of the slice carrying `g`, via acausality of `C`) gives
`D(x,y) ≤ ℓ^ren(m) = lim D(x_n,y)`. ∎

### The measure ν **[proved, given Lemmas X, C]**

For `x_1 < x_2 < y_1 < y_2` cyclically in `C̃` put

`ν([x_1,x_2) × [y_1,y_2)) := ½ [D(x_1,y_1) + D(x_2,y_2) − D(x_1,y_2) − D(x_2,y_1)]` (diagonals minus sides).

By Lemma X this is `≥ 0`; it is finitely additive on such rectangles (telescoping), `Γ`-invariant, and
symmetric under exchanging the two factors. Since `C̃` is dense in `S¹` and `D` is continuous (Lemma C), the
usual Lebesgue–Stieltjes construction gives a unique `Γ`-invariant Radon measure `ν` on the space `𝒫` of
unordered pairs of distinct points of `S¹` (= unoriented geodesics of `H²` for `X_0`), and

**(N1)** `ν({x} × J) = 0` for `x ∈ C̃` and `J` compact away from `x` (continuity of `D` at `x`).

Normalisation check (vacuum / pure hyperbolic): if `D(x,y) = 2 log|x−y| + f(x) + f(y)` (horocyclic
renormalisation), `ν` is the Liouville measure `dx dy/(x−y)²` with `ν(geodesics crossing a segment) = length`.
In the vacuum spacetime `D = ½(D_L + D_R)` and `ν = ½(m_L + m_R)`, the measure used in PROOFS §4 with
`w_L = w_R = ½`.

## 4. Lengths are ν-integrals of intersection numbers (Bonahon-type identities)

### Lemma N3 (closed classes; G2 condition (ii)) **[proved mod (MM)(b), (R3); uses Lemma X]**

Let `1 ≠ h ∈ Γ` with attracting/repelling fixed points `h^±`, and `x, y ∈ C̃` on the two sides of the axis
`(h^−, h^+)`, oriented so that `h` moves both sides from `h^−` towards `h^+`. Put
`d_n := D(x, h^n y) − D(x, h^{n−1} y)`. Then `d_n` is non-decreasing and `lim_{n→∞} d_n = L(h)`. Consequently,
with `F_h := [x, hx) × (h^+, h^−)_y` (a fundamental domain for `⟨h⟩` acting on pairs linking the axis),

`ν(F_h) = L(h)`.

*Proof.* Monotone: for `n ≥ 2` the points `x < hx < h^n y < h^{n−1} y` are cyclically ordered, and by
invariance `ν([x,hx) × [h^n y, h^{n−1} y)) = ½[D(x,h^ny) + D(x,h^{n−2}y) − 2D(x,h^{n−1}y)] = ½(d_n − d_{n−1}) ≥ 0`.
Hence `lim d_n = lim_n D(x, h^n y)/n` (Cesàro), possibly `+∞`.
Upper bound: fix `σ`; let `η` be a curve in `σ̄̃` from `x` to `y`, `β̃` an `h`-invariant lift of a closed curve of
class `h` on `σ`. `η` meets `β̃` at some `p` (linked ends). The path `η|_{[x,p]} · β̃|_{[p, h^n p]} · h^n(η|_{[p,y]})` runs
from `x` to `h^n y` and has length `ℓ(η) + n ℓ(β)`. So `ℓ_σ(κ_{x,h^ny}) ≤ ℓ_σ(κ_{x,y}) + n ℓ_σ(h) ≤ D(x,y) + nL(h)`;
take `sup_σ`: `D(x,h^ny) ≤ D(x,y) + nL(h)`, so `lim d_n ≤ L(h)` (no maximin needed).
Lower bound: fix `σ` and work in the annulus `A_σ = σ̄̃/⟨h⟩`, where `x, y` project to opposite ends. An arc
from `x` to `h^n y`, closed up by a fixed arc from `y` to `x` (joined at the cut-off at bounded cost, (R3)), is a
closed curve of winding number `n` in `A_σ`. *In a Riemannian annulus a closed curve of winding `w` has length
`≥ |w|·ℓ_1`*, `ℓ_1` the least length of a winding-one curve (induct on self-intersections: a simple closed
curve has `|w| ≤ 1`; splitting at a double point gives `w = w_1 + w_2` and `L = L_1 + L_2 ≥ (|w_1|+|w_2|)ℓ_1`;
approximate rectifiable curves by generic ones). And `ℓ_1 = ℓ_σ(h)`. So `ℓ_σ(κ_{x,h^ny}) ≥ nℓ_σ(h) − C_σ`,
hence `lim d_n ≥ ℓ_σ(h)` for every `σ`, i.e. `≥ L(h)`.
Measure: `F_h` is the increasing union over `N, M → ∞` of `[x,hx) × [h^N y, h^{−M} y)`, whose measure is, by
invariance, `½[(D(x,h^Ny) − D(x,h^{N−1}y)) + (D(x,h^{−M−1}y) − D(x,h^{−M}y))] = ½(d_N + d'_{M+1})`, where `d'`
is the same sequence for `h^{−1}`. Both tend to `L(h) = L(h^{−1})`. ∎

So (ii) of OPEN_GAPS G2 holds, with the translation length *defined* as the class-maximin length `L(h)` —
which by Lemma P / (MM)(b) is the length of the *minimal* (achronal) closed extremal geodesic in the class;
non-minimal closed extremal geodesics in the same class are invisible to `D`.

### Lemma N2 (no mass on the limit set) **[proved, given N1, N3, (MM)(b)]**

`ν` gives measure zero to pairs with an endpoint in `Λ_0` (the Gromov boundary). Hence `ν`-a.e. pair has both
endpoints in `C̃ ∖ Γ·P̃`, i.e. defines a test class (PROOFS §1.4), and PROOFS (1.1)/(B0) hold `ν`-a.e.

*Proof.* (a) Fixed points. Let `u = h^+`, `J_0` a fundamental interval for `h` on one side. If
`ν({u} × J_0) = c > 0`, then `ν({u} × h^n J_0) = c` for all `n`. Take an arc chord `(x, y_0)`, `x, y_0 ∈ C̃` on opposite
sides of the axis; the pairs `(u, v)` with `v` in the translates `h^nJ_0` lying beyond `y_0` (infinitely many,
accumulating at `h^−`) all link the chord and avoid small windows at `x, y_0`; these pairs lie in the
rectangle `R(x_±, y_{0,±})` of Lemma N4, which has finite measure. So `c = 0`.
(b) Other points. Let `{κ_j}` be a finite filling family of closed classes (none needed for the disc; the
annulus has `Λ_0 = {h^±}`, done by (a)). The `X_0`-geodesic ray to a point of `Λ_0` that is not a fixed point
of a peripheral element crosses `∪κ_j` infinitely often (otherwise it ends in a complementary disc —
impossible — or in a peripheral annulus, and then it spirals onto the boundary of the convex core and ends at
a peripheral fixed point). Distinct crossings are with distinct lifts. Unfolding,
`∫_{𝒫/Γ} #(lifts of κ_j linked) dν̄ = ν(F_{h_j}) = L(κ_j) < ∞` (N3, (MM)(b)), so the set where the integrand is
infinite is null. ∎

### Lemma N4 (arcs, with windows) **[proved, given N1, Lemma C]**

For `p̃ ≠ q̃ ∈ C̃` and windows `W_p = (p_−, p_+) ∋ p̃`, `W_q = (q_−, q_+) ∋ q̃` in `C̃`, the pairs linking the chord
`(p̃, q̃)` with no endpoint in `W_p ∪ W_q` form the rectangle `R = (p_+, q_−) × (q_+, p_−)`, and exactly

`ν(R) = ½[D(p_+,q_+) + D(p_−,q_−) − D(p_+,p_−) − D(q_+,q_−)]`, so
`D(p̃,q̃) = ν(R) + λ(W_p) + λ(W_q) + η(W)`, `λ(W_p) := ½D(p_+,p_−)`, `η → 0` as the windows shrink (Lemma C).

`λ(W_p)` depends only on the window (not on `q̃` or on the class), and is the same for all `Γ`-translates.

### Corollary I (length = intersection with ν) **[proved, given N1–N4]**

Fix windows `W_p` around one lift of each marked point `p ∈ P`, extended equivariantly, and let
`𝒫_W ⊂ 𝒫` be the pairs with no endpoint in any window. For a cut class `c` (closed components `κ`, arcs
`p → q`) put `i_W(ν, c) := ∫_{𝒫_W/Γ} i(τ(ĝ), c) dν̄(ĝ)`. Then

`L(c) = i_W(ν, c) + Σ_{ends of c at p} λ(W_p) + err_W(c)`, `err_W(c) → 0` as the windows shrink.

*Proof.* Unfold: a closed component contributes `ν(F_h ∩ 𝒫_W) ↑ ν(F_h) = L(h)` (N3; the removed pairs
decrease to pairs with an endpoint in `Γ·P̃`, null by N1). An arc component contributes `ν(R ∩ 𝒫_W)`, and
`ν(R) − ν(R ∩ 𝒫_W)` is the measure of pairs in `R` with an endpoint in another window; these sets are
dominated by a fixed finite-measure rectangle and decrease to a null set (N1), so this tends to `0`; then
use N4. By N2 the integrand is `i(τ, c)` a.e. (PROOFS (1.1)). ∎

## 5. Main theorem

**Theorem M.** Let `M` be as in §1 (NEC, any topology, rotating, time dependent), assume (MM), (MM2), (R),
and let the parties be unions of elementary pieces of one boundary Cauchy slice `C` (intervals and whole
boundaries). Let `Σ_l α_l S(A_l) ≥ Σ_r β_r S(B_r)` admit a contraction map.
1. If it is balanced at every endpoint (`Δ_p = 0` ∀p; automatic when all parties are whole boundaries),
   it holds for the renormalized HRT entropies.
2. Otherwise `Δ_p ≥ 0` (PROOFS Lemma A (A2)) and, for cut-off entropies `S^ε = S + |∂·| log(1/ε) + o(1)`, the
   difference of the two sides is `Σ_p Δ_p log(1/ε) + O(1) → +∞`.

**Status: [proved mod (MM), (MM2), (R)].**

*Proof of 1.* By Theorem H' choose the HRT classes `γ_l = c_{A_l}`, so `S(A_l) = L(γ_l)`. Apply PROOFS Lemma A
on `(X_0, P)` — a purely topological statement — to get `δ_r ∈ 𝒞(B_r)` with
`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)` for every test class `τ`. Integrate over `𝒫_W/Γ` against `ν̄ ≥ 0` (finite
integrals by Corollary I): `Σ_r β_r i_W(ν, δ_r) ≤ Σ_l α_l i_W(ν, γ_l)`. By Corollary I and `Δ_p = 0` (the window
terms `λ(W_p)` cancel point by point),
`Σ_l α_l L(γ_l) − Σ_r β_r L(δ_r) ≥ err_W → 0`, and the left side does not depend on the windows, so it is
`≥ 0`. Finally `S(B_r) ≤ L(δ_r)` (Theorem H', `≤` half). Hence
`Σ_l α_l S(A_l) = Σ_l α_l L(γ_l) ≥ Σ_r β_r L(δ_r) ≥ Σ_r β_r S(B_r)`. ∎
*Proof of 2.* Each `S^ε(A)` differs from `S(A)` by `|∂A| log(1/ε) + o(1)` (each HRT end contributes the
universal divergence). ∎

**What the proof uses, in one line.** Positivity (`ν ≥ 0`) comes from SSA on the universal cover in
class-restricted form (Lemma X); the "Bonahon identity" `L = i(ν, ·)` is elementary given positivity and
the two-sided estimate of Lemma N3; the HRT identification needs only maximin and the componentwise null
projection (Lemma P). Pleated slices, a common slice for different HRT surfaces, and uniqueness of extremal
geodesics are never needed.

## 6. Full RT cone for whole-boundary parties (Theorem 3'-NEC) **[proved mod (MM), (MM2), (R)]**

If all parties are unions of whole boundary circles (`P = ∅`), `(S(A))_A` lies in the RT holographic entropy
cone. *Proof.* PROOFS §5.3 verbatim with `Σ_j w_j m̄_j` replaced by `ν̄`: by N2, `ν̄`-a.e. pair is a test class;
discard the classes homotopic into `∂_∞X` (`i ≡ 0`); for any other test class `τ` and a peripheral class `β` at
one of its ends, `ν̄(τ)·i(τ,β) ≤ i(ν,β) = L(β) < ∞` (N3, (MM)(b)), so the push-forward is a countable positive
combination `μ = Σ μ_τ τ` with `i(μ, c) = L(c)` for closed multi-classes (Corollary I without windows). Then
Lemma G, Lemma F and closedness of the cone, exactly as there; and `S = min_c i(μ, c)` by Theorem H'. ∎
For intervals the full-cone statement is open here for the same reason as in vacuum (PROOFS §5.4).
