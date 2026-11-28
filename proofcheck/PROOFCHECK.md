# Independent proof-check of `prl_hrt/supplement.tex` (SM) and `prl_hrt/letter.tex`

Scope: mathematical correctness only (wrong inequalities, missing quantifiers, unverified hypotheses, circularity,
non-sequiturs, Letter statements stronger than the SM proves). Significance and style are not judged. No existing file
was edited. Numbering below is the SM's printed numbering (checked against `supplement.pdf`); line numbers are lines of
`supplement.tex` / `letter.tex`.

Scripts written for this check (all in this directory, run with `python3 -I`):

| script | what it checks | result |
|---|---|---|
| `check_kk_identity.py` | Prop. S18.2 / Letter Eq. (KKric), exact rational 2-jet computation, m=1 and m=2, generic fibre metric and gauge fields; tensor form: Δ := Ric_ĝ − [Ric_D^hor + du⊗du + ¼Tr(P P) + ½G F·F] is pure trace | residual **exactly 0** for m=1,2; controls (F-coefficient 1, P-coefficient ½) fail. Output `kk_output.txt` |
| `check_finite_graph.py` | Theorem S18.9, telescoping inverse (S-T9) and (S-T8), on random nonnegative weights (some zero), n=3..8, every union of atoms; independent chord-arrangement + max-flow code | T8∘T9 = id to 1e-16; min-cut = NC-matching minimum to 4e-15 in all 24 trials; boundary-only Bell graph fails (control). Output `finite_graph_output.txt` |
| `check_misc.py` | (1) endpoint integral of Prop. S4.4 and the end constant of Remark S4.5(d); (2) the explicit strip at a rotating horizon (Sec. S10); (3) Lemma S14.4 bookkeeping in vacuum; (4) the 3d Jacobi coefficients ρ = 1 − T(k,l), ⟨R(k,τ)τ,k⟩ = T(k,k) of Sec. S17 | all identities hold (the equivariance residual printed is identically zero after expanding exponents; checked by hand). Output `misc_output.txt` |

---------------------------------------------------------------------------------------------------------------------

## 1. Proposition S15.2 (Bell-pair decomposition with matter), lines 1685–1730

**Verdict: the derivation is correct *as a deduction from* Lemmas S14.1–S14.6 / Cor. S14.7, but it inherits three gaps
(one major) from Sec. S14, and one of them (G2 below) bites exactly at the step where S15.2 reuses Theorem S5.2.**

Line-by-line:

* *Coefficients* (1702–1706). Non-trivial arc class with end pieces in different components of Ω₀: every geodesic of the
  class crosses the bounding peripheral axis β at least once, so ν̄(τ) ≤ ν(F_β) = Θ(β) = L(β) < ∞. Correct (uses S14.4,
  S14.1). Same component: Lemma S4.3(ii) is topological, the closures of I, J are disjoint and their endpoints are lifted
  marked points (pieces on a lifted line with marked points are bounded by marked points; a line without marked points is
  a single piece and then k=0, trivial), so ν(I×J) is given by Eq. (S-nu) and is finite. Correct.
  **But** "ν̄ is carried by test classes and closed classes are null" is quoted from Lemma S14.6, which does not cover the
  axis pair of a peripheral element (gap G2). If ν had an atom on a peripheral axis of a boundary carrying marked points,
  that closed class is non-trivial, ν′ is *not* a combination of arc classes, and the construction of Theorem S5.2
  (π_τ = m̄₁|_τ / m̄₁(τ)) is undefined for it (m̄₁(τ)=0).
* *Endpoint classes* (1707–1716). τ̃_p = (a₋, p̃)×(p̃, a₊) is the correct fundamental set (it is the matter analogue of
  {a<q₋, b>q₊}; a_± are lifted marked points because every lifted boundary line with a marked point carries infinitely many
  of them, and for the disc they are the neighbouring marked points; |P|≤2 on the disc is correctly excluded).
  τ̃_p ∩ G_W = [a₋⁺, p₋]×[p₊, a₊⁻] is right (no other point of Γ·P̃ lies between a₋ and a₊). The algebra
  ν(τ̃_p∩G_W) + ½D(p₊,p₋) = ½[D(a₋⁺,p₊)+D(p₋,a₊⁻)−D(a₋⁺,a₊⁻)] is correct (checked by hand from (S-nu) with
  (x₁,x₂,y₁,y₂) = (a₋⁺,p₋,p₊,a₊⁻)). The limit "→ κ_p by continuity of D at the three pairs" needs **joint**
  continuity of D (both arguments move); Lemma S14.3 asserts separate continuity only and its upper-semicontinuity proof
  is wrong (G1).
* *Identity* (1717–1722). Correct given Cor. S14.7: L(c) fixed, err_W → 0, bracket → κ_p, hence i′_W(c) converges, and by
  monotone convergence (G_W increases to pairs with no endpoint in Γ·P̃, complement null by (N1) and countable additivity)
  the limit is i(ν′,c). Finiteness follows. The τ_p-part equals e_p(c)·ν(τ̃_p∩G_W) by Lemma S4.2 (topological) and
  Eq. (S-B0pre) a.e. (needs S14.6).
* *Decomposition and cone* (1723–1729). 4G_N S = G_ν + Σ κ_p b_p follows from Theorem S13.4 and e_p(c) = b_p(A) on Cut(A)
  (Def. S2.1(i)). Running Theorem S5.2 with η-minimisers: lower bound (surgery) unchanged, upper bound G_ν+η; limit of graph
  vectors, and C_N is closed and equals the set of graph vectors (polyhedral), so G_ν is itself a graph-model vector — the
  Letter's wording "graph-model vector G" (letter.tex 307) is therefore fine. The Q-statements follow by linearity and (R4).
* *Finite covers* (inherited, G3). κ_p uses D(a₋,p̃), D(p̃,a₊), D(a₋,a₊). If p is the only marked point on its boundary
  circle, a_± = h^{±1}p̃, and these are arcs from p to p in M; their finiteness and the identity D = L̂ (Lemma S14.1)
  require a finite cover in which the two endpoints become *distinct* points. The covers constructed in S14.2 step 2 only
  avoid elements producing *linking* translates; they must also avoid the finitely many elements identifying endpoints
  (here h, h²). Fixable in one line, but as written the case is not covered.

Statement-strength remark: Letter 306–308 says interval parties "obey every valid inequality in the sense of
Theorem 2(ii)". In the matter case the slack is ε-independent only for the *renormalised* entropies; regulated ones carry
o(1) corrections (SM Remark S15.3(2) says so correctly). The Letter should say "for the renormalised entropies".

---------------------------------------------------------------------------------------------------------------------

## 2. Renormalised Crofton formula, decomposition, Theorem S5.2, Theorem 2 of the Letter

**Verdict: correct.** Answers to the specific questions:

* **Lemma S4.2 (= FULLCONE Lemma 1.1), lines 328–338.** Lifts not ending at p̃=∞ have both ends in
  Λ ∪ Γ·P̃∖{∞} ⊂ [q₋,q₊] ⊂ (a,b), hence are unlinked with g=(a,b) ∈ τ̃_p. Lifts ending at ∞ correspond bijectively to arc
  ends of c at p because Stab(p̃) is trivial. **e_p(c) ≤ 1 is forced for every cut, not only minimisers:** Def. S2.1(i)
  demands exactly one arc end at each p ∈ ∂A and none elsewhere, and an arc cannot have both ends at p for the same
  reason. The HRT candidates of Defs. S7.3/S12.1 have the same property, so the decomposition is exact on all of Cut(A),
  including every non-minimal competitor in the min. After reduction (Fact S2.2(d)) no two components share a geodesic
  representative (distinct closed classes; arcs are determined by their ends), so the count is not spoiled by repeated
  components.
* **The measure computation (lines 386–394).** ∫_{v₀}^{Y²/u} dv/(u+v)² = 1/(u+v₀) − u/(u²+Y²); integrating over
  u ∈ (u₀, Y²/v₀) gives the displayed closed form, whose limit is log(Y/d_p) (u₀+v₀ = d_p). The degenerate cases u₀=0 or
  v₀=0 (x₀ = q_±, possible since x₀ ∈ Γ·P̃) also give log(Y/d_p). Verified numerically (`check_misc.py`), as is the end
  constant 1 − log 2 of Remark S4.5(d) (m{s<2Y} = log(2Y/d) − 1 + d/2Y). The second term m{h<η_ε} → 0 (empty for small ε
  if u₀v₀>0; finite and decreasing otherwise). "Exactness" (396–399) is a correct bootstrap: I₁ is monotone and equals a
  bounded expression + o(1), so its limit is finite and the identity holds exactly for every ε ≤ ε₀ by Fact S4.1.
* **Lemma S4.3.** Correct (closed classes are single points; arcs crossing a peripheral axis are bounded by Crofton; arcs in
  one component of Ω with k≥2 marked points on the inner side have disjoint end closures).
* **Corollary S5.1.** Correct: on Cut(A) the functionals differ by the constant Σ_p t_p b_p(A); finiteness of
  minimisers from discreteness.
* **Theorem S5.2 (graph).** The arrangement is the union of chosen minimisers c_A for **all** 2^N − 1 non-empty regions
  A ⊂ {1,…,N} (line 459), so for every region both bounds are available: "≤" uses that c_A ⊂ 𝒜 (cells in W_A; an
  edge separates iff it lies on c_A; cut weight = ∫#(ĝ∩c_A)dν̂ = i(μ′,c_A)); "≥" uses only the surgery of Lemma S3.1 on
  an arbitrary admissible cell set, which produces a genuine cut for A with one arc end per point of ∂A (the paired rays
  at p are joined along ∂H_p, so for ν̂-a.e. ĝ missing the small horoballs #(ĝ∩∂U′) = #(ĝ∩∂U)), whence
  G(A) ≤ i(μ′,δ) ≤ weight(S,Sᶜ). Weights are finite (≤ i(μ′,c_A)); ν̂ ≪ m̄₁ gives transversality and vertex-avoidance a.e.;
  the class-by-class normalisation makes ∫#(ĝ∩c)dν̂ = i(μ′,c) without any boundary map. Terminals: every piece lies in one
  cell (in the p̃-chart the region x>q₊ near the piece meets no lift of 𝒜), and pieces of different parties are separated
  by c_{i}. The Bell edges add exactly t_p b_p(A) to every admissible cut. ε₁ is the correct threshold for t_p ≥ 0 and
  matches the Letter with w = ½,½.
* **Corollaries S5.3, S5.4, Lemma S3.1.** Correct (A2 by the contraction property applied to x^{(i_a)}, x^{(i_b)}; A1 by
  counting crossings of a minimal-position test representative; parity at marked points from the outermost sectors).
* **Letter Theorem 2 / Prop. 1** are exactly SM Theorem S11.1 / Prop. S4.4; no over-statement. (Wording only: Letter 200
  "all other ideal points" means limit points and lifted marked points, as in the SM.)

---------------------------------------------------------------------------------------------------------------------

## 3. Proposition S10.1 (smooth slices via the conformal double), lines 1047–1131

**Verdict: correct; only cosmetic remarks.**

* (1) Doubling: gluing two copies of E_Γ along E_Γ∩∂E is open in S²×ℝ; ι commutes with the holonomy in the ambient model
  Q(X,w) = −det X + w²; free/proper action follows by folding a compact set (F commutes with Γ). The ι-invariant
  representative g_d + ι*g_d is pointwise conformal. Fine.
* (2) Folding: with an ι-invariant metric F is an isometry on each half, so a future causal Lipschitz curve folds to a
  future causal curve; inextendibility is preserved (fibre of F over a seam point is one point, over an interior point two
  separated points). F⁻¹(Σ̄₀) is met by every inextendible causal curve and is acausal (a causal curve between mirror points
  folds to a closed causal curve, impossible in M̄ since Σ̄₀ is Cauchy). Geroch applies. Fine.
* (3) Endpoint parametrisation verified: [(X(t),1)] = [(p+u²q, u√N)] (multiply by u = e^{−t}); the point is Q-null, the
  tangent (2uq, √N) is Q-orthogonal to it and has Q-norm N > 0 (det q = 0), smooth through u=0 with reflection u ↦ −u.
  Each arc doubles to a smooth spacelike circle; distinct arcs have distinct endpoints, so K is an embedded compact
  1-manifold, acausal as a subset of the doubled pleated slice.
* (4) Ribbon: the codimension-one ribbon is exactly what makes Bernal–Sánchez Theorem 1.1 applicable (compact, acausal,
  spacelike submanifold *with boundary* of a globally hyperbolic spacetime). The acausality argument (closedness of J⁺,
  acausality of K, then local temporal function in a causally convex neighbourhood) is correct.
* (5) Restriction: transversality to the timelike seam; inextendible causal curves of M̄ are inextendible in M̂; so Σ̄ is
  an acausal Cauchy slice and its trace a Cauchy curve of each cylinder; Theorem S7.7 gives (S1)–(S3).
* (6) Compactness via the peripheral relation: π₁(S̄) = Γ (contractible lift), the boundary relation gives a non-zero
  element of ker(H₁(∂S̄)→H₁(S̄)), while for non-compact connected S̄, H₂(S̄,∂S̄) ≅ H⁰_c(S̄) = 0 would force injectivity.
  Correct (Lefschetz duality for non-compact manifolds with compact boundary). It applies to every connected oriented
  acausal Cauchy slice, in particular to the pleated Σ̄₀, which (7) needs.
* (7) Homology transport: the hitting map Σ̄₀ → Σ̄ along a complete timelike field tangent to the boundary is a continuous
  bijection of compact Hausdorff spaces, fixes γ̄_c pointwise and preserves boundary intervals. Correct.
* Corollary S10.2 and Letter Theorem 1 ("minimum attained in both classes") follow.

Remark (not an error): the cited "Remark 4.14" of Ref. [BernalSanchez2006] for disconnected K could not be checked here;
Theorem 1.1 itself does not require connectedness of K in the formulations known to me.

---------------------------------------------------------------------------------------------------------------------

## 4. Matter section (S12–S17)

| item | verdict |
|---|---|
| (JM) (1286–1290) | Used correctly in: Thm S13.4 "≤" (S12.5(c)); Lemma S14.1 "≥" (single class in a cover); Lemma S14.2 step 3 (joint, nested pair in a cover); Lemma S14.4 upper bound (via S14.1). Gap: (JM) does not say the minimisers m_i are *simple*, but S12.5(c) needs simple components for a candidate (line 1362 asserts it). Add "simple" to (JM) (or prove it by the exchange argument of S16.3(7)). |
| Prop. S12.4 ((JM) in vacuum) | Correct *given Lemma S13.1* (hence (R1), (R4)), as the SM says. The Letter (283–285) says "It holds in vacuum" without this proviso. (R1) is essentially the CDGH area theorem, so this is minor. |
| Prop. S12.5 | (a),(b) correct (inequalities close up). (c) correct modulo the simplicity point above; Epstein's theorem applies to disjoint, pairwise non-homotopic essential curves/arcs. |
| Remark S12.6 (the corrected statements) | Correct. (i) linked lifts must cross on any slice, so disjoint minimisers of both classes cannot exist on any admissible slice; (JM) is only assumed for embedded multi-classes. (ii) correctly hedged ("need not"); D, Θ are defined on covers where every chord/core class is simple. |
| Lemma S13.1 (null projection) | Correct in structure (Wall's lemma applied to one side of one component in the cover M̃/H where it is the whole edge; (i) and (ii) are right). Gaps (minor, fixable): (a) for arcs ℓ(ĝ)=∞, so "a closed set of finite H¹-measure … contains a simple arc from one end of ĝ to the other" (1409–1411) must be run in {z ≥ ε} and closed up near the anchors with (R4)/(R3); (b) the statement should assume g simple (its lift separates σ̂₀) — true in every application; (c) Lemma S14.1 applies it to non-invariant slices of the universal/annular cover, which the statement ("σ anchored at C′" for M′ finite cover) does not literally include, though the proof covers it. |
| Lemma S13.2, Prop. S13.3, Thm S13.4 (Theorem H′) | Correct. Sandwich uses sup inf ≤ inf sup and sup Σ ≤ Σ sup in the right direction; deletion of null-homotopic/repeated closed extremal components lowers length. |
| Lemma S14.1 | Correct (the "≤/≥" labels are swapped relative to the displayed equality but both directions are proved: invariant slices lift, so L̂ ≤ D; Lemma S13.1 in the cover gives D ≤ ℓ^ren(m) = L̂). The finiteness claim "every class becomes simple in some finite cover (arcs: residual finiteness)" needs G3. |
| Lemma S14.2 (crossing = SSA on the cover) | Correct: step 1 exchange of tails gives ℓ(a→d)+ℓ(b→c) ≤ D(a,c)+D(b,d) on every slice; step 3 gives sup ≥ D(a,d)+D(b,c) via (JM) for the pair in a finite normal cover. Gap G3: step 2 must also exclude the finitely many γ with γ{a,d} ∩ {a,b,c,d} ≠ ∅ or γ{b,c} ∩ {a,b,c,d} ≠ ∅ (otherwise in M̂ the two arcs share an endpoint, or an arc becomes a loop at one boundary point, and the pair is not an embedded multi-class). |
| **Lemma S14.3 (continuity)** | **Error in the proof (G1).** The "upper semicontinuity" paragraph (1541–1545) proves D(x,y) ≤ ℓ^ren(m) = lim D(x_n,y) — that is again *lower* semicontinuity. Lemma S13.1 only ever gives D ≤ (length of an extremal curve on an admissible slice); to get limsup D(x_n,y) ≤ D(x,y) one needs the limit curve m to be an ℓ_σ-*minimiser* on its slice, or uniform control of joining costs over the maximising slices, neither of which is in (R2)/(R3) ((R3) constants "depend on σ only"). The proven properties (Γ-invariance, crossing inequality, lsc, Lemma S14.4) are all compatible with discontinuous D, e.g. D₀ − c·1[x∈S] − c·1[y∈S] for an invariant discrete S. Moreover only separate continuity is addressed, while S14.5 (η(W)→0) and S15.2 need joint continuity. |
| Construction of ν (1549–1558) | Needs continuity (countable additivity of the Lebesgue–Stieltjes premeasure, (N1)); so it inherits G1. Positivity and finite additivity are correct (diagonals ≥ sides is Lemma S14.2 with (a,b,c,d) = (x₁,x₂,y₁,y₂)). |
| Lemma S14.4 (d_n limits) | Correct. Cyclic order x < hx < hⁿy < h^{n−1}y and ν([x,hx)×[hⁿy,h^{n−1}y)) = ½(d_n − d_{n−1}) checked; upper bound by concatenation on a Γ′-invariant slice with h ∈ Γ′ (Hall), powers of h give parallel (unlinked) translates; lower bound by the winding inequality on annular-cover slices; ν(F_h) = ½(d_N + d′_{M+1}) → Θ(h). Vacuum sanity check in `check_misc.py` (d_n ↑ log λ). Needs G3 for the cover. Countable additivity (increasing union) needs G1. |
| Lemma S14.5 | Algebra correct; η(W) → 0 requires joint continuity (G1). |
| **Lemma S14.6 (no mass on Λ₀)** | **Gap (G2).** (a) shows ν({h⁺}×J₀) = 0 for fundamental intervals J₀ on either side, which covers {h⁺}×(S¹∖{h⁺,h⁻}) but **not the single pair (h⁺,h⁻)**. For non-peripheral h this pair is covered by (b) (the axis crosses a filling family infinitely often). For **peripheral h** (and the core of the annulus, where (a) is said to "suffice") neither (a) nor (b) excludes an atom ν({(h⁺,h⁻)}) > 0. Nothing in (JM)/(R) visibly forbids it. Consequences: harmless for Thm S15.1(iii) (P=∅: peripheral classes are trivial), harmless for Thm S15.1(i) (integration of Lemma S3.1's inequality works for closed test classes too, B0pre holds with the parallel-copy caveat), but it breaks the claim in S15.2 that ν′ is a combination of *arc* classes and the π_τ-normalisation in the reuse of Theorem S5.2. Fix: either prove the absence of such atoms, or (simpler) extend Theorem S5.2 to closed classes with positive weight by using a fixed transverse representative τ′ in minimal position with 𝒜 (parallel copy if τ ⊂ 𝒜), w(e) += μ(τ)#(τ′∩e); both bounds go through verbatim. |
| Cor. S14.7 | Correct given G1, G2 (dominated/monotone convergence of the window deficits is right). |
| Thm S15.1 (Theorem 4 of the Letter) | The logic of (i)–(iii) is correct (window terms cancel point by point when Δ_p = 0; (ii) by (R4); (iii) via the S5.2 graph with η-minimisers and closedness of C_N). All three parts depend on ν being a measure with the identities of S14.4–S14.7, hence on G1. |
| Prop. S16.1–S16.3, (K), (Reg) | Not needed for Theorem 4. S16.3 step (8) acausality argument is correct. But Letter 284–285 paraphrases (Reg) as "a maximising slice is C² and spacelike near its minimal curves"; (Reg) also requires the maximiser to be uniformly spacelike and Lipschitz *everywhere* (used in steps (1) and (8)), and S16.3 proves (JM) only for M, while (JM) is required in all finite covers ((K) lifts to finite covers; (Reg) would have to be assumed there). |
| Prop. S17.1 (Jacobi system) | Coefficients verified symbolically from the 3d Schouten form with Λ=−1, 8πG=1: ρ = ⟨R(k,τ)τ,l⟩ = 1 − T(k,l), ⟨R(k,τ)τ,k⟩ = T(k,k) = R_kk; the (p,q) system and its cooperativity under the NEC are right; the elementary maximum principle step is correct. (Lighter check of the Wronskian sign bookkeeping.) |

How (JM) enters, summarised: realisation of every cut class by a candidate (Theorem H′ "≤"); D = L for simple classes (S14.1);
SSA on the cover (S14.2 step 3); closed lengths (S14.4 upper bound); and, in the intended but currently incorrect argument,
continuity (S14.3). The NEC enters only through Lemma S13.1. No circularity was found (S13.1 → S12.5 → S13.4; S12.5 → S14.1 →
S14.2 → S14.3 → S14.4/5/6 → S14.7 → S15.1/S15.2).

---------------------------------------------------------------------------------------------------------------------

## 5. Strips (S18)

**Verdict: correct as a conditional theorem under (H1)–(H4).**

* **KK / effective NEC identity (Prop. S18.2, Letter Eq. KKric):** verified independently and exactly (`check_kk_identity.py`):
  for m=1 and m=2 with generic G_ij(x) and A^i(x), Ric_ĝ − Ric_D^hor − du⊗du − ¼Tr(P_aP_b) − ½G_ijF^i_acF^j_b^c is pure
  trace, so the identity holds on all null k; perturbing either coefficient breaks it. The isotropic value m(m+1)(kφ)² and
  the positivity argument (PSD of g⁻¹ on the annihilator of a null vector; trace of a product of PSD matrices) are correct.
* **Averaging/coarea (Lemma S18.3):** correct (convexity and isometry invariance of total variation; for invariant f,
  |Df|(Σ) = V∫λ|df|; weighted coarea; a.e. level set has trace A since the trace of f is χ_A).
* **Relative-cutoff lemma (S18.6):** correct (BV gluing in a layer, L¹ convergence of slice traces, averaging, and gluing to
  X's own collar so that spatial minimality applies). The edge-band estimate H^{d−1}(E_ε)/V ≲ ε^{−(d−1)}δ_ε has the right
  power; same-anchor asymptotics supply δ_ε = o(ε^{d−1}).
* **Matching formula (S-T7):** correct under (H3), (H4) (lower bound: invariant HRT surface → disjoint arcs, each an
  admissible single-strip extremal candidate; upper bound: laminar joint maximin, assumed).
* **Finite graph theorem (S18.9) and (S-T9):** telescoping verified by hand and numerically; min-cut = matching minimum
  verified on random instances including zero weights (`check_finite_graph.py`); the "Bell graph fails for disconnected
  unions" remark and the 2 min(w₀₂, w₁₃) gap are correct.
* **UV shift a_p (S-T10):** correct: D* = d + a_a + a_b (a≠b, D*_aa=0) changes only the adjacent weight w_{p−1,p}, by +a_p;
  non-adjacent weights involve each endpoint once with each sign, so local counterterms cancel and SSA gives w⁰ ≥ 0 there;
  F_{D*} = F_d + Σ_{∂A} a_p because each matching uses each endpoint once.
* Minor: the statement "for a<b<c<d the required inequalities are SSA and weak monotonicity" (2187) is redundant — only the
  "diagonals ≥ nested sides" combination enters w_ij for non-adjacent atoms (and on a circle WM is SSA for the rotated
  quadruple). Harmless.
* Lemma S18.8: the compactness step needs a finite cover of 𝓜 by neighbourhoods on which b ≠ 0 in the interior to fix one
  cutoff χ; implicit but standard.

---------------------------------------------------------------------------------------------------------------------

## 6. Vacuum core (Lemma S3.1, Theorem S9.1, Theorem S7.6/S7.7, Theorem S11.2) — lighter pass

**Verdict: correct.** Separation Lemma S8.1 (sign pattern of the positive lift, AM–GM, equality case) and Lemma S8.2
checked; Prop. S6.1(b) truncation length log N + 2 log(1/ε) checked; Theorem S7.6 Step 3 (limit point argument) checked;
Theorem S9.1 "≥" (unlinked endpoints via Jordan, uniqueness of the invariant geodesic, parity reconstruction of W) and "≤"
(Ψ-image of W̃) checked; the explicit horizon strip of Sec. S10 checked symbolically (`check_misc.py`). Theorem S11.2 uses the
collar-lemma realisation of graphs (cited to Bao et al.), fine.

Terminology error (Letter 257, 262; SM 1198, 1207, 1214): "static multiboundary wormholes". For non-abelian π₁ the
time-symmetric spacetimes with X_L = X_R admit **no** Killing field at all (a Killing field of M lifts to an
Ad(Γ)-invariant element of sl₂⊕sl₂, and a non-elementary Fuchsian group fixes no non-zero vector of sl₂), so they are not
static. Every use in the proofs only needs time-reflection symmetry (RT on the fixed slice). Replace "static" by
"time-symmetric (non-rotating)". The Letter sentence (105) "for X_L ≠ X_R and non-abelian π₁ it is not even stationary" is
true but implicitly suggests the X_L = X_R case is stationary, which it is not.

---------------------------------------------------------------------------------------------------------------------

## Ranked list of issues to fix before submission

1. **(Major, matter) Lemma S14.3: the upper-semicontinuity argument is wrong (it re-proves lsc), and joint continuity is
   never addressed.** Continuity of D is used for the construction of ν as a measure, (N1), Lemma S14.5 (η(W)→0),
   Cor. S14.7, Thm S15.1 (i)–(iii) and Prop. S15.2 (κ_p limit at three pairs). *Fix:* strengthen the hypotheses, e.g.
   require in (JM) that the maximising slice can be taken in S̄*(C′) with one fixed orthogonality bound (C₀, δ₀) for all
   multi-classes and covers (then D(x_n,y_n) = ℓ_{σ_n}(x_n→y_n) ≤ ℓ_{σ_n}(x→y) + o(1) ≤ D(x,y) + o(1) by joining at the
   cutoff with uniformly small cost, giving joint usc; lsc is already proved), or add continuity of D to (R), or require
   that limits in (R2) of the (JM) minimisers be minimisers on their limit slices. The Letter's hypothesis list must be
   updated accordingly.
2. **(Moderate, matter) Lemma S14.6 does not exclude an atom of ν on the axis pair (h⁺,h⁻) of a peripheral element (in
   particular the annulus core).** Affects the claim in Prop. S15.2 that ν′ is a combination of arc classes and its reuse of
   Theorem S5.2. *Fix:* prove absence, or extend Theorem S5.2 to positively weighted closed classes via fixed transverse
   representatives (straightforward); also restate S14.6 accordingly.
3. **(Minor–moderate, matter) Finite-cover choices** in Lemma S14.2 step 2 (and hence S14.1, S14.3, S14.4) must also avoid
   the finitely many elements that identify endpoints of the chords (γx ∈ {endpoints}); otherwise the classes in the cover
   are loops at one boundary point or share endpoints and are not an embedded multi-class. Needed concretely in S15.2
   when p is the only marked point on its circle (a_± = h^{±1}p̃).
4. **(Letter, statement strength) Theorem 4's hypotheses** are described as "(JM) and the asymptotic regularity stated in the
   SM" (letter 289), but (R) also contains (R1) (Wall's projection/area theorem) and (R2) (precompactness of extremal
   families with convergence of renormalised lengths), and after fix 1 a further uniformity assumption. State them as
   regularity *and compactness* hypotheses.
5. **(Letter, statement strength) Sufficient conditions for (JM)** (letter 283–285): (Reg) also requires the maximiser to be
   uniformly spacelike and Lipschitz globally, and Prop. S16.3 proves (JM) only in M, not in its finite covers; "It holds in
   vacuum" is proved given Lemma S13.1, i.e. given (R1), (R4).
6. **(Minor, matter) Lemma S13.1:** run the separating-continuum/finite-H¹ argument in {z ≥ ε} for arcs (length is infinite),
   assume g simple, and state it for non-invariant slices of the universal and annular covers as used in S14.1.
7. **(Minor, matter) (JM) should state that the minimisers are simple** (used in Prop. S12.5(c)).
8. **(Terminology, vacuum) "static multiboundary wormholes"** → "time-symmetric"; they have no Killing field for
   non-abelian π₁ (Letter Thm 3 and SM S11.2–S11.3).
9. **(Minor, Letter 306–308)** "every valid inequality in the sense of Theorem 2(ii)" with matter holds with ε-independent
   slack only for the renormalised entropies (regulated ones have o(1) corrections, SM Remark S15.3(2)).
10. **(Cosmetic)** Lemma S14.1's "≤/≥" labels are reversed relative to the displayed identity; S18 "weak monotonicity"
    is redundant; Letter 200 "all other ideal points" should read "all other limit points and lifted marked points".

## Overall verdict

* **Theorem 1 (Theorem H, both slice classes), Proposition 1 (renormalised Crofton), Theorem 2 (G an exact finite graph
  min-cut, S^ε ∈ C_N for ε ≤ min(ε₀,ε₁), every valid inequality with the stated slack behaviour) and Theorem 3 (cone
  equality):** proved at the strength stated, for the stated class, given the cited inputs. The only correction is the
  terminology "static" (item 8), which does not affect any argument.
* **Theorem 4 (matter):** **not proved at the strength stated.** The reduction (Theorem H′), the crossing inequality and the
  closed-length identity are sound, but the construction of ν and all window/endpoint limits rest on continuity of D,
  whose proof in Lemma S14.3 is incorrect (item 1). With an added hypothesis of the kind described in item 1 (and the
  small repairs 3, 6, 7) Theorem 4(i)–(iii) would follow as written. The interval-party extension stated in the Letter
  (306–308, SM Prop. S15.2) additionally needs item 2.
* **Theorem 5 (strips):** proved conditional on (H1)–(H4) as stated in the Letter; the KK identity and the finite graph
  theorem were independently verified.
