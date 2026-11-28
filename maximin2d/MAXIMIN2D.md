# MAXIMIN2D: removing the maximin hypotheses from the 2+1 non-vacuum theorem

Working notes, appended as I go. Nothing outside `maximin2d/` is edited. Builds on `nonvacuum/NONVACUUM.md`
(cited as NV), `PROOFS.md`. Grades: **[proved]**, **[proved mod X]** (X an exact cited result or a named
hypothesis), **[conj]**, **[false]** (with counterexample).

## 0. Log
* Report 1. Read NV fully, PROOFS §6.3b–6.4. Plan: (i) restructure so each hypothesis is used as little as
  possible; (ii) prove what compactness + 2D gives for free; (iii) literature check on what is proved
  about maximin (Wall, MWW, AEPU, GMW, HHLR); (iv) try to avoid maximin for crossing positivity.
* Report 2. Literature read (§1). Restructuring (§2), the 2+1 facts that are free (§3), a 2D resolution of
  Wall's non-uniqueness gap given regularity of the maximin slice (§4), and a maximin-free *local* proof of
  crossing positivity from NEC via a cooperative Jacobi system, with numerics that also show NEC is
  necessary (§5). Status table and verdict (§§6–7).

## 1. Literature: what is actually proved (marks: [R] text read via ar5iv/arXiv HTML, [A] abstract only)

| Source | What it proves | What it assumes |
|---|---|---|
| Wall, arXiv:1211.3494v4 [R] | Thm 9: min(A,Σ) exists on a fixed achronal slice (compactness + lsc of area). Thm 10: maximin slice exists on *horizonless* spacetimes (Arzelà–Ascoli on slices, area usc). Thm 11: also with Kasner-like singularities, because area `∝ t^{1−p_{D−1}} → 0`. Thm 14 (NEC): maximin surface has no null-separated points. Thm 15: a *stable* maximin surface is extremal. Thm 17: for `A ⊃ B`, maximise `c₁Area[M₁]+c₂Area[M₂]` on one slice; minimisers non-crossing; ⇒ m(A), m(B) minimal on a common slice. Thm 18 (SSA), 19 (MMI) from 17. | Def. 8b + §3.5: "we will proceed under the assumption that stable maximin surfaces can always be found"; Thm 12: unique ⇒ stable; several ⇒ one stable "at least if Σ is smooth and spacelike". Discussion: results "only follow if one either assumes or proves the existence of the needed maximin surfaces". |
| Marolf–Wall–Wang, arXiv:1901.03879 [R] | Restricted maximin (slices with `Σ∩∂M = C∂`, exactly our `𝒮(C)`): extremal, equals HRT when in a smooth region; existence when anchored slices reach singularities only at Kasner-like ones; extends to generic charged/spinning holes. | NEC, generic condition, AdS-hyperbolicity; stability "proven only when the slice is both spacelike and smooth … it remains an assumption more generally". Homology only; no homotopy classes, nothing 2+1-specific. |
| Akers–Engelhardt–Penington–Usatyuk, arXiv:1912.02799 [A] | Quantum maximin; "first general proof" of EWN and SSA for the QES prescription. | Existence argued, not proved; perturbative `G_Nħ`. |
| Grado-White–Marolf–Weinberg, arXiv:2008.07022 [A] | With a radial cutoff, *HRT* can violate SSA without time-reflection symmetry; *restricted maximin* obeys SSA, EWN, MMI (NEC). | Maximin existence/stability as in MWW. Relevant warning: HRT and maximin can differ once the setting is non-standard. |
| Headrick–Hubeny–Lawrence–Rangamani, arXiv:1408.6300 [R, first 100k chars] | NEC ⇒ any spacelike-homologous extremal surface lies in the causal shadow (boundary causality). Footnote 11: "all our theorems apply to any spacelike-homologous extremal surface" (no minimality, no maximin). | **Does not prove SSA** (the only occurrence of SSA is about static RT). The premise "HHLR proved SSA for HRT" in the task is not borne out by the text read. |
| Headrick–Hubeny, arXiv:2208.10507 [A] | Covariant bit threads / minimax reformulations of HRT. | Starts from HRT. |
| Grado-White–Grimaldi–Headrick–Hubeny, arXiv:2502.09894 [R, first 100k chars] | Minimax time-sheets; "we will assume that the sup and inf are achieved"; assumes existence of maximin/HRT surfaces; "cooperating" configurations proved for two crossing regions only, and "the cooperating conjecture … may be false". | — |

**Searches for Theorem H'** (`S_HRT(A) = min_c Σ_k sup_σ ℓ_σ(c_k)`, i.e. HRT = min over cut classes of
class-restricted maximin, `sup_σ min_c = min_c sup_σ`): nothing found [S: two web searches; Balasubramanian et
al. 1406.2663 and Maxfield 1412.0687 minimise over geodesics in homotopy classes only in *vacuum*, where the
geodesic in each class is unique; Wall's Cor. 17h — several disjoint regions minimal on one slice — is the
closest covariant statement]. Verdict: not in the literature as stated; any HRT expert would regard it as an
"expected" consequence of maximin, i.e. folklore-adjacent, *but it is exactly as conditional as maximin*.
No rigorous (GMT-level) proof of maximin existence + stability in any generality was found (search not
exhaustive).

## 2. Restructuring: three functionals, and what each identification costs

Notation of NV §1. For a region `A` put
* `S_mm(A) := sup_σ min_σ(A)` (restricted maximin value),
* `S_L(A) := inf_{c∈𝒞(A)} L(c)`, `L(c) = Σ_k L(c_k)`, `L(κ) = sup_σ ℓ_σ(κ)` (class maximin),
* `S_HRT(A)` (NV §1).

**Proposition 2.1 (sandwich) [proved mod Lemma P].** `S_mm(A) ≤ S_L(A) ≤ S_HRT(A)`.
*Proof.* On each slice, `min_σ(A) ≤ ℓ_σ(c)` for every `c ∈ 𝒞(A)` (NV Lemma Q, proved); take `sup_σ` then
`inf_c`, using `sup_σ inf_c ≤ inf_c sup_σ` and `sup_σ Σ_k ≤ Σ_k sup_σ`. Second inequality: NV Theorem H', `≥`
half (each HRT component is an extremal geodesic on an acausal slice; Lemma P). ∎
(The second inequality needs **no** maximin hypothesis; the first needs none either.)

**Proposition 2.2 (the contraction argument only sees `S_L`) [proved, given NV's ν-machinery].** If `ν ≥ 0`
and the Bonahon identities (NV Cor. I) hold for `L`, then `Σ_l α_l S_L(A_l) ≥ Σ_r β_r S_L(B_r)` for every
balanced contraction-map inequality.
*Proof.* NV Theorem M, proof of 1, with `γ_l` a near-optimal class for `S_L(A_l)` (inf, so up to `ε`) and
the final step `S_L(B_r) ≤ L(δ_r)`, which is now the definition. ∎

So the hypotheses enter in exactly two places:
* **positivity/identities for `L`** (NV Lemmas X, S, C, N3): used (MM)(b), (MM2), (R2), (R3);
* **identification `S_HRT = S_L`**: by 2.1 only `S_HRT ≤ S_L` is missing, i.e. *for every region `B` and every
  cut class `δ` there is an HRT candidate for `B` of length `≤ L(δ)`*. This is an **existence** statement
  for extremal curves and cannot be obtained from focusing/projection arguments, which only produce upper
  bounds on `ℓ_σ` (lower bounds on extremal lengths). Some variational existence theory (maximin) is
  unavoidable here. It follows from the standard "HRT = maximin" (Wall/MWW) since `S_mm ≤ S_L` is free:
  `S_HRT = S_mm ≤ S_L`.

**Proposition 2.3 (one hypothesis implies all three) [proved].** Let (JM) be: *for every embedded
multi-class `κ = {κ_1,…,κ_m}` (pairwise disjoint simple representatives) in `M` or in a finite cover, the
sup of `F(σ) = Σ_i ℓ_σ(κ_i)` over 𝒮(C) is attained on a slice carrying pairwise disjoint minimisers `m_i`
that are complete extremal geodesics and whose union with `C` is acausal.* Then (MM)(a) in the form
`S_HRT ≤ S_L`, (MM)(b) and (MM2) hold.
*Proof.* (MM)(b): `m = 1`. (MM2): `Σ ℓ(m_i) = max F ≤ Σ L(κ_i) ≤ Σ ℓ(m_i)` (Lemma P for each `m_i`, after
extending `∪m_i ∪ C` to an acausal Cauchy slice — standard extension result, Bernal–Sánchez CMP 2006 for the
smooth compact case; the version with timelike boundary is Aké–Flores–Sánchez arXiv:1808.04412-type and was
not checked line by line). `S_HRT(B) ≤ L(δ)`: apply (JM) to `δ ∈ 𝒞(B)`; the `m_i` are disjoint simple
curves on `σ* ≅ Σ̄` componentwise homotopic to the components of an embedded multicurve bounding a region
with trace `B`; by Epstein's theorem (Acta Math. 115, 1966: componentwise-homotopic disjoint essential simple
multicurves/arcs are ambient isotopic, rel boundary for arcs) the region is carried to one bounded by `∪m_i`,
so `(σ*, ∪m_i, R)` is a candidate of length `Σ ℓ(m_i) = Σ L(δ_k) = L(δ)`. (Inessential/peripheral
components: delete, as in NV.) ∎

So the four hypotheses of the task collapse to **(JM) + (R)**; (JM) is Wall's Theorems 10/11 + 15 + 17
with the minimisation restricted to homotopy classes.

## 3. What 2+1 gives for free (no stability, no extremality)

**Set-up.** Fix a Cauchy temporal function on `M̄` adapted to the timelike boundary, `M̄ ≅ ℝ × Σ̄` with
`C = {0} × ∂Σ̄` (Aké–Flores–Sánchez, arXiv:1808.04412, cited, not re-proved), metric
`g = −N²dt² + h_t(dx + βdt, dx + βdt)`. Every achronal set met once by each `t`-line is a graph `t = f(x)`.
Let `𝒮̄_*(C)` be the *achronal* slices anchored at `C` satisfying the boundary-orthogonality bound
`|f(z,φ) − t_C(φ)| ≤ C_0 z^{1+δ_0}` in Fefferman–Graham coordinates (fixed `C_0, δ_0 > 0`). Define `L`, `D`, `Λ`
with sups over this class (or its lifts). This loses nothing: a slice that is boosted at an endpoint has
`ℓ^ren = −∞` there, and once (JM) holds `L` over any class between `𝒮_*` and all achronal slices is the same
number `ℓ(m)` (Lemma P is valid for any achronal target slice).

**(K) Compactness hypothesis.** The closure of the Wheeler–DeWitt set `W(C) := M̄ ∖ (I^+(C) ∪ I^−(C))` in
`M̄` is compact. (Every slice of 𝒮(C) lies in `W(C)`.)

### Proposition U (upper semicontinuity in the slice) **[proved mod Ioffe's theorem and FG asymptotics]**
If `σ_n → σ` uniformly (graphs `f_n → f`) in `𝒮̄_*(C)`, then `limsup ℓ_{σ_n}(κ) ≤ ℓ_σ(κ)` for every class.
The same holds for `F = Σ_i ℓ_σ(κ_i)` and, locally uniformly, for slices of covers.

*Proof.* Let `η ⊂ σ` be a competitor with spatial path `x(s)` and project it vertically:
`η_n(s) = (f_n(x(s)), x(s)) ∈ σ_n`. With `w = (f∘x)'`, `g(η',η') = Q(s,t,w) := h_t(x'+βw, x'+βw) − N²w²`, a
quadratic polynomial in `w` with leading coefficient `h(β,β) − N² < 0` (`∂_t` timelike): **concave in `w`**,
and `≥ 0` a.e. on achronal slices. Hence `Ψ(s,t,w) = √Q` (−∞ where `Q < 0`) is concave in `w`, usc in `(t,w)`,
bounded above by `C|x'|`. `f_n∘x → f∘x` uniformly with `(f_n∘x)' ⇀* (f∘x)'` in `L^∞` (uniform Lipschitz bound on
compact sets), so by Ioffe's lower-semicontinuity theorem for normal convex integrands (SIAM J. Control
Optim. 15 (1977) 521) applied to `−Ψ`, `limsup ∫Ψ(s,f_n∘x,(f_n∘x)') ≤ ∫Ψ(s,f∘x,(f∘x)')` on `{z ≥ δ}`. On
`{z < δ}` the orthogonality bound and FG asymptotics give tails `log(δ/ε) + O(δ^{δ_0})` uniformly in `n` for
radial competitors, and on `σ` a competitor can be made radial near its ends at cost `O(δ^{δ_0})`. So
`limsup ℓ_{σ_n}(κ) ≤ ℓ^ren(η) + O(δ^{δ_0})`; take `η` near-optimal and `δ → 0`. ∎

**The opposite semicontinuity is [false].** Minkowski `ℝ^{2,1}`, `σ = {t = 0}`, `σ_n = {t = f_n(x)}`, `f_n` a
zigzag of slope `±0.99` and amplitude `1/n`. `σ_n → σ` uniformly, but the curve `y = const` on `σ_n` has length
`√(1−0.99²) ≈ 0.14` per unit `x`, so the distance on `σ_n` between the points over `(0,0)` and `(1,0)` is
`≤ 0.14`, while on `σ` it is `1` (the same happens with boundary-anchored classes). Length is
therefore **not** lower semicontinuous under convergence of slices; minimisers on `σ_n` need not converge to
minimisers on `σ`. This is why the maximisation step is easy (usc is the right direction) and the
extremality step is hard (§4). The task's premise "length is lower semicontinuous" is right only on a
*fixed* spacelike slice.

### Proposition E (existence of maximin slices) **[proved mod (K), Prop U]**
Under (K), `𝒮̄_*(C)` is compact in `C^0` (Arzelà–Ascoli: uniformly Lipschitz graphs over compact `Σ̄`;
uniform limits of achronal graphs are achronal since `I^±` are open; the orthogonality bound is closed), so
every usc functional — `ℓ_σ(κ)`, `Σ_i ℓ_σ(κ_i)`, `min_σ(A)` — attains its sup. Moreover `L(κ) < ∞` for every
class, with no maximin input: a curve on an achronal slice is no longer than its spatial path measured in
`h_t + (β-terms)`, which is bounded on the compact set. (Finiteness of `D`, `Λ` thus no longer needs Scott's
theorem or (MM)(b).)

**(K) is restrictive.** Already for eternal BTZ with `t_L = t_R = 0`, `W(C)` is the Penrose diamond
`|y| ≤ 1−|x|`, whose top and bottom vertices lie on the `r = 0` circle: the closure is not compact (marginally).
With infalling NEC matter the singularity moves towards the past and, I expect, `W(C)` meets it in an open
set [conj; not computed]. So (K) holds for horizonless spacetimes (stars, solitons) but essentially *never*
for the multiply connected NEC spacetimes of interest (wormholes always have horizons by topological
censorship). Wall's Kasner barrier (Thm 11) needs the minimal *area* to vanish near the singularity; for
curves in 2+1 that requires both Kasner exponents positive. For BTZ-type degenerations (exponents (1,0), only
the `φ`-circle shrinks) it fails for arc classes. **Existence (hypothesis (1)) is therefore proved here only
under (K); in the black-hole case it is exactly as open as in Wall/MWW, and slightly worse for arcs.**

### Lemma C revisited (continuity of `D`) **[proved mod (K), FG; (MM)(b) and (R2) are no longer used]**
`D(x,y) = sup_σ̃ ℓ_σ̃(x→y)` is lsc as a sup of functions continuous in the endpoints (orthogonality bound).
usc: for `x_n → x` take `σ̃_n` with `ℓ_{σ̃_n}(x_n→y) ≥ D(x_n,y) − 1/n`; a subsequence converges locally
uniformly (local compactness of the lifted family under (K)); joining at the cut-off,
`ℓ_{σ̃_n}(x_n→y) ≤ ℓ_{σ̃_n}(x→y) + o(1)`, and Prop U gives `limsup ≤ ℓ_σ̃(x→y) ≤ D(x,y)`. ∎

### Minimisers on one slice **[proved for uniformly spacelike Lipschitz slices; not needed otherwise]**
On a uniformly spacelike Lipschitz slice the induced length structure is a compact length space
bi-Lipschitz to the background; Arzelà–Ascoli and lsc of the length-space length give minimisers in every
class (renormalisation via the orthogonality bound). On achronal slices with null pieces the induced distance
is only a pseudo-metric (null generators have length 0), and minimising sequences need not have bounded
background length; existence of curve minimisers there is not proved here. Nothing below uses it.

## 4. Extremality: a 2D resolution of the non-uniqueness gap, and what is left (regularity)

**(Reg) Regularity hypothesis** for an embedded multi-class `κ = {κ_1,…,κ_m}` (`m = 1` allowed): some
maximiser `σ*` of `F(σ) = Σ_i ℓ_σ(κ_i)` over `𝒮̄_*(C)` is uniformly spacelike and Lipschitz, and is `C²` on an
open neighbourhood `U` of `K := ∪_i ∪𝓜_i`, where `𝓜_i` is the set of `ℓ_{σ*}`-minimisers in `κ_i`.

### Proposition W2 **[proved mod (K), (Reg), Sion's minimax theorem, FG asymptotics, and the acausal-extension theorem]**
Under (K) and (Reg), for each `i` some `m_i ∈ 𝓜_i` is a complete extremal geodesic of `M`; the `m_i` are
simple and pairwise disjoint, and `∪m_i ∪ C` is acausal. Hence **(JM) holds for `κ`** (Prop 2.3).

*Proof.* (1) `σ*` uniformly spacelike Lipschitz ⇒ `(σ*, d_{σ*})` is a compact length space (after the
cut-off), minimisers exist and near-minimising sequences subconverge uniformly, with convergence of length,
to minimisers. In `U` minimisers are `C²` geodesics of the `C²` metric of `σ*`; `𝓜_i` is compact in `C¹`
there (geodesic ODE).
(2) *No transversal crossings.* Two minimisers, or a minimiser and a deck translate of a minimiser of the
same or another `κ_j`, cannot cross transversally: lift to the universal cover; the lifts have non-linking
endpoints (simple, non-linking classes), so a transversal crossing comes with a second crossing; exchanging
the arcs between them gives two curves with the *same endpoints*, hence the same classes, the same total
length, and corners; cutting a corner shortens one of them below the minimum. (Same endpoints with one
crossing: exchange at that point.) So two minimisers through `p ∈ U` are tangent at `p`, hence coincide near
`p` (ODE uniqueness on a `C²` surface). Therefore `h(p) := ⟨H_M(p), n(p)⟩` — the spacetime curvature vector of
any minimiser `M ∋ p`, paired with the unit normal `n` of `σ*` — is a well-defined continuous function on `K`.
(Its component tangent to `σ*` vanishes: `M` is a geodesic of `σ*`.)
(3) *Danskin.* For `φ ∈ C²_c(U ∖ C)` let `σ_t` be `σ*` pushed along `tφn` (still in `𝒮̄_*` for small `|t|`; it
equals `σ*` off `supp φ`). With `A_φ(M) := −∫_M hφ ds` (first variation of length):
`d/dt|_{0+} ℓ_{σ_t}(κ_i) = min_{M∈𝓜_i} A_φ(M)`. Upper bound: push each `M`. Lower bound: pull back by the flow,
`ℓ_{σ_t}(Ψ_t c) = ℓ_{σ*}(c) + tA_φ(c) + O(t²)` uniformly for curves of bounded length (the change is supported in
`supp φ`, where everything is `C²`); take near-minimisers `c_t`, which converge (subsequence, by (1)) to some
`M ∈ 𝓜_i` with `c_t' → M'` in `L¹` on `supp φ` (uniform convergence + convergence of length in a Riemannian
region), so `A_φ(c_t) → A_φ(M)`.
(4) Maximality of `σ*` gives `Σ_i min_{M_i∈𝓜_i} A_φ(M_i) ≤ 0` for every `φ`.
(5) *Minimax.* `(φ, (μ_i)) ↦ Σ_i ∫A_φ dμ_i` is bilinear, continuous in `μ = (μ_i) ∈ Π_i P(𝓜_i)` (compact
convex, weak-*) and in `φ ∈ V := C²_c(U∖C)` (convex). Sion's theorem:
`min_μ sup_φ Σ_i∫A_φ dμ_i = sup_φ min_μ Σ_i∫A_φdμ_i = sup_φ Σ_i min_{M_i} A_φ(M_i) ≤ 0`. Linearity in `φ` makes
the left sup `0` or `+∞`; so some `μ` has `Σ_i ∫A_φ dμ_i = 0` for **all** `φ`.
(6) *The measure argument (this is where 2D is used).* Let `λ := Σ_i ∫ H¹⌊M dμ_i(M)`, a Radon measure on `U`
carried by `K`. By (2), `Σ_i∫A_φ dμ_i = −∫hφ dλ`, so `hλ = 0`, i.e. `h = 0` λ-a.e.; hence for every `i` and
`μ_i`-a.e. `M`, `h = 0` `H¹`-a.e. on `M`, so `h ≡ 0` on `M` (continuity): `M` is a spacetime geodesic. Each `μ_i`
is a probability measure, so such an `m_i` exists.
(7) Distinct classes give distinct geodesics of `σ*`; by (2) they cannot touch without coinciding: disjoint.
Simplicity: a self-crossing of a minimiser in a simple class is excluded by the same exchange.
(8) *Acausal.* If `p, q ∈ ∪m_i` were causally related, achronality of `σ*` forces `q ∈ J^+(p) ∖ I^+(p)`, and every
point of the connecting null geodesic lies in `J^+(σ*) ∖ I^+(σ*)` and `J^−(σ*) ∖ I^−(σ*)`, hence in `σ*`; but a null
geodesic cannot lie in `σ*` near `p`, where `σ*` is spacelike. Pairs involving `C`: the orthogonality bound
makes the compactified slice `C¹` and spacelike (conformal metric) at `C`. Extension of the acausal set
`∪m_i ∪ C` to an acausal Cauchy slice: Bernal–Sánchez (Commun. Math. Phys. 257 (2005) 43; Lett. Math. Phys. 77
(2006) 183) in the boundaryless smooth case; with timelike boundary the analogous result is expected from
Aké–Flores–Sánchez — **cited, not checked**. ∎

**Relation to the literature.** Wall's Thm 12 and MWW already assert that stability holds when the
maximin slice is smooth and spacelike; W2 is a complete, self-contained version of that statement for
curves, which in addition (i) handles a *continuum* of minimisers via the measure `λ` (the step that is
genuinely 2D: minimisers are curves on a surface and cannot cross, so their curvatures do not cancel on
average), and (ii) handles several classes at once (Wall's Thm 17 structure). It does **not** remove the
real gap; it locates it:

**What is left: (Reg).** Nothing proved here or in the literature guarantees a maximiser that is smooth and
spacelike near its minimisers. Why it is hard: maximisers form a large non-compact-in-`C¹` family (any
modification away from the minimisers that creates no short-cut is still a maximiser), uniform limits of
smooth slices can acquire null pieces exactly where lengths are being maximised (cf. the zigzag example,
§3), and on a null piece the length functional is degenerate, so the first-variation step (3) has no meaning
there. In vacuum the pleated slices of PROOFS §6.6 are maximisers but are bent; for closed geodesics PROOFS
§6.6a produces smooth spacelike slices, and for arcs I expect (not checked) a pleated slice bent only away
from the arc, which would satisfy (Reg) [conj]. A proof of (Reg) in general would require a regularity
theory for a Lorentzian sup–inf problem; I know of none.

**Hypothesis (2) of the task, sharpened.**
* "the maximin curve is an extremal geodesic": **[proved mod (K), (Reg)]** (W2 with `m = 1`).
* "of minimal length among extremals in κ": restricted to extremals lying on acausal slices of `M̄` anchored at
  `C`, **[proved mod Lemma P]** (`L(κ) ≤ ℓ(g)` for every such `g`, and the maximin one attains `L`).
  Unrestricted it is **[false]**: in the maximal analytic extension of rotating BTZ the inner bifurcation
  circle `r = r_−` is a closed extremal geodesic in the horizon class with length `2πr_− < 2πr_+ = L(h)`; it lies
  beyond the Cauchy horizon, on no slice of `𝒮(C)`.
* For classes without simple representatives the statement is **[false]** already in vacuum (NV §1).

## 5. Avoiding maximin

### 5a. Null projection / focusing alone (route (a)) **[not viable as a replacement; explanation]**
Every focusing argument (Wall's representative, Lemma P, HHLR's causal-shadow theorem) takes an extremal
surface and produces, on *another* slice, a curve that is *no longer*. It therefore gives **upper** bounds
on `ℓ_σ` (`sup_σ ℓ_σ(κ) ≤ ℓ(g)`) and **lower** bounds on extremal lengths. The two places where the theorem
needs more are of the opposite type: (i) `S_HRT(B) ≤ L(δ)` asks for an extremal configuration *no longer*
than a sup — an existence statement; (ii) crossing positivity asks for an *upper* bound on `D(a,d) + D(b,c)`,
which by Lemma P would need extremal geodesics in the classes `ad`, `bc` with total length
`≤ D(a,c)+D(b,d)`; focusing applied to `m_ac, m_bd` produces no extremal curves at all. HHLR (1408.6300) prove
causality (extremal surfaces in the causal shadow) for any spacelike-homologous extremal surface without
maximin, but **not** SSA (text checked, §1). The one known general proof of HRT-SSA is Wall's, through joint
maximin (Thm 17–18). So (a) cannot remove (1) or (3).

### 5b. Proposition J (local crossing positivity from NEC via a cooperative Jacobi system) **[proved mod (ND) and FG asymptotics]**

Let `m` be a complete spacelike extremal geodesic, unit tangent `τ`, arclength `s ∈ ℝ`, and `k, l` a
*parallel* frame of future null normals, `⟨k,l⟩ = −1` (exists: the normal bundle of a geodesic is preserved
by Levi-Civita transport), `k` pointing spatially to side O of `m`, `l` to side I. Write a normal Jacobi field
(`J'' = −R(J,τ)τ`) as `J = αk + βl`. Using the 3D identities (checked numerically, `check_jacobi_identities.py`)

`⟨R(l,τ)τ,l⟩ = Ric(l,l)`, `⟨R(k,τ)τ,k⟩ = Ric(k,k)`, `ρ := ⟨R(k,τ)τ,l⟩ = 1 − T(k,l)` (`Λ = −1`, `8πG = 1`),

the Jacobi equation is the 2×2 system

`α'' = ρα + Ric(l,l)β`,  `β'' = Ric(k,k)α + ρβ`.

With `p := −α`, `q := β`: `𝓛(p,q) := (−p'' + ρp − Ric(l,l)q, −q'' + ρq − Ric(k,k)p) = 0`. **Under NEC the
off-diagonal couplings `−Ric ≤ 0`: the system is cooperative.** (Vacuum: decoupled, `ρ = 1`.) The invariant
cone `{p, q ≥ 0}` is `{−pk + ql}` = the *inward spacelike* cone of normals.

**(ND)** `𝓛` admits a strict positive supersolution on `ℝ` (`φ_i ≥ c > 0`, `𝓛φ ≥ cφ`). Sufficient:
`T(k,l) + max(T(k,k), T(l,l)) ≤ 1 − c` along `m` (take `φ = (1,1)`); holds in vacuum and for weak matter.

**Statement.** Let `(x,y) ↦ m_{xy}` be a `C²` family of extremal geodesics on an open set `U` of pairs of
boundary points with `D = ℓ^ren(m_{xy})` on `U`, each satisfying (ND). For boundary null directions at the
endpoints labelled by the frame (`k`-type, `l`-type; first index at `x`, second at `y`):

`D_{ll} ≤ 0,  D_{kk} ≤ 0,  D_{lk} ≥ 0 (> 0),  D_{kl} ≥ 0 (> 0)`,  with `D_{ll} = D_{kk} = 0` in vacuum.

Hence along any spacelike boundary curves through `x` and `y` (tangent `a l − b k` at `x`, `a'k − b'l` at `y`,
`a,b,a',b' > 0`) the crossing density `∂_s∂_t D = aa'D_{lk} − ab'D_{ll} − ba'D_{kk} + bb'D_{kl} > 0`: **`ν|_U > 0`**.

*Proof.* Near the ends (FG falloff) `ρ → 1`, `Ric(k,k), Ric(l,l) → 0`, so solutions behave as `c_± e^{±s}`.
Moving `x` along an `l`-type (resp. `k`-type) boundary null vector gives the Jacobi field `J_x` with growing
mode at `−∞` equal to a positive multiple of `l` (resp. `k`), and no growing mode at `+∞`; it exists and is
unique by (ND). The mixed second derivative of the renormalised length is the Wronskian
`W(J_x,J_y) = ⟨J_x',J_y⟩ − ⟨J_x,J_y'⟩` (constant along `m` because `R(·,τ)τ` is self-adjoint), and evaluating at
`+∞` pairs the decaying coefficients `(c_p, c_q)` of `J_x` with the growing coefficient of `J_y`:
`D_{·k} = κ c_q`, `D_{·l} = −κ c_p` with a universal `κ > 0` (fixed by the vacuum value
`D = log Δu + log Δv`). For `l`-type data `(p,q)(−∞) = (0,+)`: the maximum principle for cooperative systems
(de Figueiredo–Mitidieri, Comm. PDE 1990; Sweers, Math. Z. 1992: it holds iff a strict positive
supersolution exists), applied on `[−S,S]` and letting `S → ∞`, gives `p, q ≥ 0`, hence `c_p, c_q ≥ 0`
(strong maximum principle: `c_q > 0`), i.e. `D_{ll} ≤ 0 ≤ D_{lk}`. `k`-type data flip the signs of `(p,q)`,
giving `D_{kk} ≤ 0 ≤ D_{kl}`. In vacuum `p ≡ 0`. ∎

*Meaning.* `D_{ll} ≤ 0` is exactly boundary SSA for null-deformed intervals
(`S(x',y) + S(x,y') ≥ S(x,y) + S(x',y')`, `x' = x + δl`, `y' = y + δl`), and the cooperative coupling is the
linearisation of entanglement-wedge nesting. The first-order-in-matter version (entanglement density ∝
integrated NEC) is in Bhattacharya–Hubeny–Rangamani–Takayanagi, arXiv:1412.5472 [S]; I did not find the
non-perturbative cooperative statement [S, one search].

### 5c. Numerical test (`jacobi_sign_vaidya.py`, `nearnull_slice.py`; logs alongside) **[numerical evidence]**
Smooth planar Vaidya–AdS3 `ds² = z^{−2}(−(1−m(v)z²)dv² − 2dvdz + dx²)`, `m = m_0 + Δm(1+tanh((v−v_0)/w))/2`;
NEC ⇔ `m' ≥ 0` (`check_nec_sign.py`: `R_ab + 2g_ab = (z m'/2) δ^v_a δ^v_b`). Exact geodesics by Hamiltonian
shooting (DOP853, rtol 1e−12, cut-off `2·10^{−4}`, mass continuation), 4-point mixed differences, `h = 0.02`
(step check `h = 0.04, 0.02, 0.01`: `D_{+−} = −0.24989, −0.24860, −0.24828`). `+` = future-right, `−` =
future-left, `x_1 < x_2`; so `D_{+−} = D_{ll}`, `D_{−+} = D_{kk}`, `D_{++} = D_{lk}`, `D_{−−} = D_{kl}`.

| geometry | `D_{++}` | `D_{−−}` | `D_{+−}` | `D_{−+}` | ν-density, slice slopes ±0.5/∓0.3 |
|---|---|---|---|---|---|
| vacuum (planar BTZ), 4 configs | 0.64 … 10.9 | 0.72 … 45.8 | `≤ 6·10^{−7}` (= 0) | `≤ 6·10^{−7}` | > 0 |
| NEC `Δm = +1.5`, 4 configs | > 0 | > 0 | −0.11 … −0.40 | −0.005 … −0.06 | > 0 |
| NEC `Δm = +4`, 4 configs | > 0 | > 0 | −0.06 … −0.45 | −0.008 … −0.22 | > 0 |
| NEC-violating `Δm = −0.8` | > 0 | > 0 | **+0.06 … +0.23** | **+0.003 … +0.04** | > 0 |
| NEC-violating `Δm = −0.95` | > 0 | > 0 | **+0.013 … +0.11** | **+0.001 … +0.06** | > 0 |

All 16 NEC sign predictions hold; vacuum cross terms vanish to solver accuracy; **every** NEC-violating
case flips the sign of both cross terms. Consequence for the crossing inequality on a legitimate
(spacelike) boundary slice that is nearly null at the endpoints (slope `+S` at `P = (0,0)`, `−S` at
`Q = (0.4,1.5)`):

| | `S = 0.9` | `S = 0.99` |
|---|---|---|
| vacuum | +0.182 | +0.019 |
| NEC `Δm = +1.5` | +0.400 | +0.265 |
| NEC-violating `Δm = −0.8` | +0.061 | **−0.117** |

So **without NEC the crossing inequality (SSA on the cover) is violated** on an admissible boundary
Cauchy slice — the necessity that NV §7 could not exhibit with the thin-shell solver (consistent with
Callan–He–Headrick, arXiv:1204.2309 [S]: SSA in AdS3-Vaidya holds iff NEC).

### 5d. What Prop J does and does not buy globally
* **[proved mod (ND)]** On any open set of pairs where `D` is the length of a `C²` nondegenerate family of
  extremal geodesics, `ν ≥ 0` there, with **no joint-maximin hypothesis**.
* **Ties.** Where two branches compete, `D = min(D_1, D_2)`, `ν` acquires the singular part
  `−∂_x g ∂_y g · δ_Γ` on `Γ = {g := D_1 − D_2 = 0}`. It is `≥ 0` iff `∂_x g ∂_y g ≤ 0` on `Γ` (**(Tie)**: the
  dominant branch switches in an SSA-compatible direction). In static RT, (Tie) follows from cut-and-paste;
  covariantly it is implied by joint maximin and I see no local proof. It would follow from nesting of the
  two tied geodesics (one in the inner wedge of the other), which is Wall's Cor. 17g — again stability.
* `D = ℓ(m_{xy})` itself needs single-class maximin ((MM)(b)) or must be taken as a definition.
* Net effect: (MM2) can be replaced by **(MM)(b) + (ND) + (Tie)** — conditions on extremal geodesics that are
  checkable in examples, but not more standard than (MM2). Prop J's real value is (i) a maximin-free proof of
  the local positivity, showing exactly where NEC enters (the two couplings `Ric(k,k), Ric(l,l)`), and (ii)
  with 5c, a demonstration that NEC is necessary.

### 5e. (K) fails in the NEC wormholes of NV §7, *because of* NEC **[proved for thin-shell Vaidya–BTZ; script `wdw_singularity.py`]**
Shell from the right boundary at `v_0 < 0` (i.e. `C = {t_L = t_R = 0}` lies after the shell, the NV §7 regime),
`M_1 → M_2`. The top vertex of `W(C)` is where the right ray `v = 0` meets the left ray from `t_L = 0`, which in
the future interior obeys `dr/dv = (r² − M(v))/2`. For `M ≡ M_1` it reaches `r = 0` exactly at `v = 0` (square
Penrose diagram, `r = tanh(−v/2)`). For `v > v_0` the right-hand side is pointwise smaller iff `M_2 > M_1`, so by
ODE comparison **under NEC the left ray hits `r = 0` at some `v* < 0`, and `W(C)` meets the singularity along an
interval**; under NEC violation the vertex is in the interior. Numbers (`wdw_singularity.log`):
`v_0 = −1`: `M_2 = 1.5 → v* = −0.352`, `1.1 → −0.098`; `0.9 → +0.122`, `0.5 → +1.21`.
Together with topological censorship (Galloway–Schleich–Witt–Woolgar, PRD 60 (1999) 104039: under NEC, a
horizonless asymptotically AdS spacetime has simply connected domain of outer communication, so in 2+1 the
slice is a disc), **(K) is essentially incompatible with the multiply connected NEC case** — the case the
theorem is about. In that case existence of (class-)maximin slices requires a barrier argument near the
singularity which, for arcs in 2+1, Wall's Kasner argument does not supply.

## 6. Status of the four hypotheses (this note)

| Hypothesis / claim | Grade |
|---|---|
| (1) maximin existence for regions: attainment of `sup_σ min_σ(A)` | **[proved mod (K)]** (Prop E: usc + Arzelà–Ascoli). Without (K): **open**; (K) fails for NEC wormholes (§5e). In the form actually needed (`S_HRT ≤ S_L`) it is implied by (JM) (Prop 2.3) or by Wall/MWW "HRT = restricted maximin". Cannot be obtained from focusing (§5a). |
| (1′) the maximin value is computed by an extremal configuration (HRT = maximin) | [proved mod (K) + (Reg) for the optimal cut class] (W2 + 2.3); in general **[standard hypothesis: Wall 1211.3494 §3.5, MWW 1901.03879]** — stability is assumed in the literature except for smooth spacelike maximin slices |
| (2) class maximin attained | [proved mod (K)] |
| (2) maximin curve extremal | [proved mod (K), (Reg)] (W2) — the non-uniqueness part of Wall's gap is closed in 2D; regularity of the maximising slice remains |
| (2) minimal among extremals in κ | [proved mod Lemma P] if restricted to extremals on slices of `𝒮(C)`; **[false]** unrestricted (inner circle of rotating BTZ); **[false]** for non-simple classes (NV) |
| (3) joint maximin, non-linking simple pairs | [proved mod (K), (Reg)] (W2 with `m = 2`); **[false]** for linking pairs (NV) |
| (3′) local crossing positivity | **[proved mod (ND), FG] without any maximin** (Prop J); numerically confirmed; NEC shown necessary (§5c) |
| (3″) crossing across branch transitions | equivalent to (Tie); [conj]; implied by (JM) |
| (4) near-boundary regularity | (R2) **no longer needed** (Lemma C via usc); (R3) built into 𝒮_* at no cost; (R1) = Chruściel–Delay–Galloway–Howard area theorem [cited]; FG falloff of matter assumed |
| `L`, `D`, `Λ` finite | [proved mod (K)] (no Scott / (MM)(b) needed) |
| Lemma C (continuity of `D`) | [proved mod (K), FG] |
| Sandwich `S_mm ≤ S_L ≤ S_HRT` | [proved mod Lemma P] |
| Contraction inequalities for `S_L` | [proved mod positivity of ν + NV Lemmas S, N3] — these need (MM)(b), (MM2) |
| Theorem H′ (sup–min = min–sup) | [proved mod (JM)]; not found in the literature [S]; folklore-adjacent |

**Net reduction.** Theorem M holds under **(JM) + (R)**, one hypothesis replacing four; and (JM) holds under
**(K) + (Reg)**. Since (K) essentially excludes the multiply connected case, the honest form is:
**Theorem M is conditional on (JM)**, i.e. on homotopy-restricted, multi-class *stable restricted maximin*.

## 7. Verdict

* **Unconditional? No.** Two genuinely open inputs remain, and both are the open problems of the covariant
  holographic literature itself, not artefacts of this proof: (i) *existence* of maximising slices when the
  Wheeler–DeWitt set reaches a singularity — which, by §5e, is the generic situation for NEC wormholes; (ii)
  *regularity/stability* of the maximiser near its minimal curves. Nothing in 2+1 removes either: curves make
  the minimisation and the non-uniqueness problem easy (Prop W2), but the maximisation over Lorentzian slices
  is a sup of a non-convex, non-lsc functional with no regularity theory.
* **What is new and unconditional-ish.** (a) The reduction to a single hypothesis (JM) (Prop 2.3), and the
  proof of (JM) from regularity alone in 2+1 (W2: the measure argument resolves Wall's non-uniqueness gap for
  curves). (b) A maximin-free local proof of crossing positivity from NEC (Prop J), with the NEC entering as the
  sign of the couplings of a cooperative Jacobi system, and numerics showing the sign flips — and the crossing
  inequality fails on an admissible boundary slice — when NEC is violated. (c) Lemma C, finiteness of `L, D`,
  existence of maximisers, without (MM)(b)/(R2), under (K).
* **What a PRL can state.**
  > *Theorem.* Let `M` be an asymptotically AdS₃ spacetime with the null energy condition, globally hyperbolic
  > with timelike boundary, of arbitrary spatial topology. Assume the restricted-maximin property of Wall and
  > Marolf–Wall–Wang for curves restricted to homotopy classes: for every finite set of disjoint simple homotopy
  > classes (in `M` and its finite covers), the supremum over anchored achronal slices of the total minimal
  > length is attained by disjoint extremal geodesics on an acausal slice. Then the HRT entropies of unions of
  > boundary intervals and boundaries satisfy every balanced inequality admitting a contraction map; for
  > whole-boundary parties they lie in the RT entropy cone.

  plus the sentence: "the assumption holds whenever the maximising slice is smooth and spacelike near its
  minimal curves and the Wheeler–DeWitt patch is precompact; in vacuum it is proved by pleated slices." It
  must **not** be called unconditional, and "HRT = maximin" must not be cited as a theorem for black-hole
  spacetimes with non-Kasner singularities — MWW themselves say stability "remains an assumption".
* **How an HRT-community referee would read it.** As a result "at the level of rigour of Wall's proof of SSA
  and MMI for HRT": accepted practice in the field, which routinely assumes maximin existence and stability.
  Points the referee would press: (1) the homotopy-class restriction of maximin is not in the literature —
  they would want the adaptation of Wall's Thms 10/11/15/17 spelled out (W2 + Prop U do this for curves,
  modulo (K)/(Reg)); (2) finite covers: that maximin is assumed in covers as well (legitimate, covers are
  spacetimes of the same class, but it should be said); (3) the reliance on non-invariant slices of the
  universal cover in `D`; (4) the boundary-slice dependence (fixed `C`). A mathematically minded referee would
  note that the conditional part is precisely the analytically hard part (existence and regularity of a
  Lorentzian sup–inf), so the theorem is a *reduction*: "all contraction-map inequalities for HRT in 2+1
  follow from class-restricted maximin", which is still a strong and clean statement (to my knowledge — searches
  not exhaustive — no prior result covers the multiply connected case of arXiv:2407.07165 with matter), but it is not a proof of the covariant
  cone for 2+1 from the Einstein equations alone.
* **Single most useful addition to the PRL**: Prop J + §5c (local positivity directly from NEC, and the
  NEC-violating counterexample to crossing), as evidence that the conditional input is not hiding a failure.
  None of this has been checked by anyone but its author.

Scripts (all `python3 -I`, in `maximin2d/`): `check_nec_sign.py` (NEC ⇔ `m' ≥ 0`), `check_jacobi_identities.py`
(3D curvature identities of Prop J), `jacobi_sign_vaidya.py` (5 min; log `jacobi_sign_vaidya.log`),
`nearnull_slice.py` (25 s; log `nearnull_slice.log`), `wdw_singularity.py` (log `wdw_singularity.log`).
