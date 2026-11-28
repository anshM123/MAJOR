# CFT2 reformulation of the vacuum covariant-cone theorem

Working notes, 28 Nov 2026. Appended as the work proceeds. Scripts in this directory run with `python3 -I`.
Grades: **[proved]**, **[proved mod cited X]**, **[conjecture]**, **[false]**. Sources: **[read]** (text read in this
session), **[abstract]** (abstract/listing only), **[memory]** (not re-checked this session).

Repo inputs used (not re-proved here): Theorem H, Prop. S4.3 (renormalised Crofton), Cor. S5.1 / Thm. S5.2 (G is a
graph min-cut), Sec. S6 length formula (`prl_hrt/supplement.tex`). Referee minor issue 3 (`referee_hrt/REPORT.md`):
G is independent of the horoball normalisation.

## 0. Scripts

| script | what it checks |
|---|---|
| `banados_map_check.py` | the Roberts map pulls Poincaré AdS3 back to the Bañados metric with L = -(1/2){f,x}, Lb = -(1/2){fb,xb} (numerical residual ~1e-13 at random points) |
| `fg_vs_horocycle.py` | in a dressed (Bañados) geometry, FG-cutoff geodesic lengths = frame-horocyclic lengths + c(p,q) eps^2 + O(eps^4), with c(p,q) **not** endpoint-local |
| `g_boundary.py` | canonical G vs the CFT-defined G_CFT = S - sum_p (1/2) I(adjacent) b_p: disc, rotating BTZ; frame (Weyl) independence; cone facets; free fermion |
| `ising_i3.py` | lattice I3 of three spin blocks (ED) in the critical Ising and XX chains, and of fermionic blocks (correlation matrix) |
| `hill_kernel.py` | t_p from ⟨T⟩ alone: chiral kernel K(x,y) = Hill solution = uniformiser formula (1e-12) |
| `adjacent_subtraction_graphs.py` | Bell-subtraction lemma on random graph models |

## 1. Boundary gravitons: already inside the repo's class (only the regulator changes)

Units: lengths l = 4G_N S; c = 3/(2G_N), so S = (c/6) l. Chiral lengths l_L, l_R; spacetime length (l_L+l_R)/2.

### 1.1 Bañados geometries are the repo's spacetimes with moved boundary slices [proved mod cited classification]

Every asymptotically AdS3 solution of pure gravity with flat boundary cylinder (x^± = φ ± t, φ ~ φ+2π) is,
in FG gauge, a Bañados metric
`ds² = dz²/z² + L(x⁺)dx⁺² + L̄(x⁻)dx⁻² + (1/z² + z² L L̄) dx⁺dx⁻`
(FG expansion terminates in 3D; Skenderis–Solodukhin hep-th/9910023 [abstract]). The map
`w = f − 2z²f'²f̄''/(4f'f̄'+z²f''f̄'')`, `w̄ = f̄ − 2z²f̄'²f''/(…)`, `u = 4z(f'f̄')^{3/2}/(…)`
pulls Poincaré `(du²+dw dw̄)/u²` back to it with `L = −½{f,x⁺}`, `L̄ = −½{f̄,x⁻}` (Roberts 1204.1982 [abstract];
**checked** in `banados_map_check.py`). Equivalently f = ψ₁/ψ₂ with ψ'' = Lψ (Hill), and the holonomy around the
circle is the Hill monodromy M₊ ∈ SL(2,R) (resp. M₋).

Classification of the orbits (Balog–Fehér–Palla hep-th/9703045 [abstract]; Garbarz–Leston 1403.3367 [abstract];
energy statements [memory]):

| (M₊, M₋) orbit type | geometry | repo status |
|---|---|---|
| both hyperbolic, winding 0 (orbit of a constant L₀>0) | rotating BTZ + boundary gravitons | **covered** (π₁ = Z) |
| both "vacuum" orbit (monodromy −1, f ∈ Diff⁺(S¹)) | global AdS3 + boundary gravitons | covered by the Γ = 1 case (Sec. 1.4) |
| elliptic (−1/4 < L₀ < 0) | conical defect = point particle (matter) | excluded; angles 2π/n via `particles/FINITE_ORBIFOLD_THEOREM.md` |
| parabolic | massless/extremal BTZ (cusp) | excluded; partial results in `particles/CUSPS_EXTREMAL_BTZ.md` |
| exotic orbits (no constant representative, winding ≥ 1, L₀ = −n²/4 with n ≥ 2) | no single-valued regular quotient / energy unbounded below [memory] | excluded: not physical states of a unitary CFT with positive energy |

**Lemma 1 (dressing = moving the boundary slice).** Let B be a Bañados geometry with hyperbolic winding-0
monodromies on each asymptotic boundary j, with uniformisers f^{(j)}_± : R → I^j_{L,R} equivariant under the
peripheral holonomy (h_j: ρ_L(h_j) = M^{(j)}_+, ρ_R(h_j) = M^{(j)}_−). Then B is isometric to (a neighbourhood of
the boundary in) the repo's spacetime M = D̃/Γ for (ρ_L, ρ_R), and the boundary map (x⁺,x⁻) ↦ (f₊(x⁺), f₋(x⁻)) is a
conformal diffeomorphism of the boundary cylinder onto (I_L × I_R)/⟨h⟩. Every boundary Cauchy curve t = τ(φ),
|τ'| < 1, goes to the graph of the increasing, h-equivariant map ψ = f₋∘(…)∘f₊⁻¹, i.e. an acausal boundary Cauchy
slice in the sense of Sec. S6 of the supplement; boundary intervals go to boundary intervals; whole boundaries to
whole boundaries. Proof: f± increasing and equivariant; acausality is preserved by increasing reparametrisations
of the two null coordinates separately. HRT (Def. S7.2) is diffeomorphism invariant. [proved, given that the bulk of
B relevant to HRT is the maximal globally hyperbolic development of its boundary data, which for these holonomies is
the Mess–Barbot quotient M — **proved mod cited** Mess 2007 / Barbot 2008 / Bonsante–Seppi as already used in the
repo (Sec. S6 "Cited input")]. Different boundaries, and the two sides of a TFD, can be dressed independently and
evolved to independent times.

So the repo's class "two convex-cocompact Fuchsian reps + arbitrary acausal boundary slice" **already contains every
boundary-graviton dressing at every time**. Nothing in the bulk changes; the only new ingredient is which cutoff the
physical boundary frame induces.

### 1.2 The induced cutoff

In the uniformised (Poincaré-type) chart the Bañados FG cutoff z = ε is, by the map above,
`u = ε √(f₊'(x⁺) f₋'(x⁻)) (1 + O(ε²))`, `w = f₊ + O(ε²)`.
In the repo's language (Sec. S6, lambda lengths): the spacetime cutoff at p is the rank-one matrix
`p = a_p b_pᵀ/ε` with the **frame lifts**
`a_p = (f₊(x⁺_p), 1)/√f₊'(x⁺_p)`, `b_p = (f₋(x⁻_p), 1)/√f₋'(x⁻_p)` (up to the repo's J-convention).
The corresponding (exact) horocyclic regulated length of a chord is
`l_frame(p,q) = log[ (f₊(p)−f₊(q))(f₋(p)−f₋(q)) / (ε² √(f₊'(p)f₊'(q)f₋'(p)f₋'(q))) ]`.

**Lemma 2 (FG vs horocycles).** For each chord, `l_FG(p,q) = l_frame(p,q) + c(p,q) ε² + O(ε⁴)`, both for the
point-to-point FG length and for the boundary-anchored geodesic truncated at z = ε; the coefficients c(p,q) are
**not** endpoint-local (`fg_vs_horocycle.py`: c₀₂+c₁₃−c₀₃−c₁₂ = ∓3.85 ≠ 0 for a random dressing). Hence
`4G S_FG = G + Σ_p t_p(ε) b_p + R(ε)`, `R = O(ε²)` not of Bell-pair form. [proved: the leading identity is exact
algebra; the O(ε²) form is the elementary estimate in H² (crossing height of a horocycle of diameter D is
D(1−O(D²/r²))) and is confirmed numerically]. This reproduces the referee's point (M2c).

**Lemma 3 (twist-field regulator = horocyclic regulator, exactly).** In CFT, S_A = lim_{n→1} (1−n)⁻¹ log of
`ε^{2Δ_n·|∂A|} ⟨Π σ_n⟩`, Δ_n = (c/12)(n − 1/n): the cutoff enters only as one normalisation factor per twist field,
so it shifts S by exactly `(c/6) log(1/ε)` per endpoint; under a conformal map the twist fields (primaries)
transform with `(f₊'f₋')^{Δ_n/2}` per point. Hence the CFT's standard conformal regularisation **is** the frame
horocyclic cutoff, and in it `4G S = G + Σ t_p b_p` holds **exactly**, with no O(ε²) remainder. [proved; standard
conformal covariance]. The O(ε²) terms of Lemma 2 are scheme artefacts of the same order as lattice corrections,
which are non-universal in a CFT; they carry no physical content.

### 1.3 What the dressing does to G and t_p [proved]

By Prop. S4.3 and the referee's minor issue 3, G depends only on the two marked hyperbolic surfaces
(X_L, P_L), (X_R, P_R), not on the horocycles. With Lemma 1:
`G_dressed(A) = G_M(A; P_L = f₊(x⁺_∂A), P_R = f₋(x⁻_∂A))`, and only t_p sees the frame:
`t_p = log(1/ε) − ½ log(d_p^L d_p^R)`, d_p computed with the frame lifts. Boundary gravitons act on the finite part
**only by moving the entangling points on the two ideal circles**. In CFT language this is the exact statement
`S_{UρU†}(A) = S_ρ(f⁻¹A) + (c/12) Σ_{p∈∂A} log(f₊'(p) f₋'(p))` for U implementing (f₊, f₋) — true in every CFT2,
holographic or not, because twist fields are primaries. So the dressing extension is CFT-elementary once the
undressed result holds on arbitrary acausal slices; the repo's arbitrary-slice theorem is exactly what is needed.

### 1.4 The simply connected case (global AdS3 + gravitons = vacuum Virasoro module)

Γ = 1, X_L = X_R = H², marked points u_p, v_p on two circles in the same cyclic order. Sections S2–S5 of the
supplement are stated for groups including the trivial one (Lemma S4.2 treats "trivial Γ"; d_p > 0 needs |P| ≥ 3).
Theorem H for Γ = 1 is the standard HRT statement in pure AdS3 (cuts = non-crossing matchings; Czech–Dong
1905.03787, Sec. 4) — the repo's pleated-slice proof is local to the disc and goes through, but I did **not**
re-check it line by line for Γ = 1. The graph theorem is checked to 1e-14 in `fullcone/disc_exact.py` and
`g_boundary.py` (disc). [proved mod repo Thm S5.2 with Γ = 1 and mod Theorem H for Γ = 1 (standard)].

### 1.5 Goal theorem (bulk form)

**Theorem B.** Let M be a globally hyperbolic, locally AdS3 spacetime (any vacuum solution of 3D gravity with
Λ < 0) whose conformal boundary is a finite union of timelike cylinders with arbitrary boundary metrics in the
conformal class (arbitrary boundary-graviton dressing), with orientable Cauchy surface of finite type, and holonomy
either trivial or a pair of convex-cocompact Fuchsian representations (all non-trivial holonomies hyperbolic).
Fix any acausal boundary Cauchy slice, N parties made of intervals and whole boundaries, endpoint set P (the points
where the party changes; |P| ≥ 3 if Γ = 1). Then, with the CFT's twist-field (= frame-horocyclic) regulator,
`4G S(A) = G(A) + Σ_p t_p(ε) b_p(A)` **exactly** for all ε ≤ ε₀, where G is the min-cut function of a finite graph,
independent of ε and of the boundary metric (Weyl/conformal frame), and b_p is a Bell pair between the parties
adjacent at p. With an FG cutoff the same holds up to O(ε²) non-Bell corrections.
Consequences: S ∈ C_N for small ε; every superbalanced valid inequality holds as an identity between finite
quantities `Q(S) = Q(G) ≥ 0` (all regulators; Q(b_p) = 0); every other valid inequality has slack → +∞.
**Grade: [proved mod repo Thms H, S4.3, S5.2 and the cited Mess–Barbot classification]**; Γ = 1 additionally
mod Theorem H for the disc (standard). Excluded and genuinely open within pure gravity: parabolic holonomy
(extremal/massless BTZ), and "exotic" Virasoro orbits (not physical). Conical defects are matter, not vacuum.

Regulator issue settled as follows: the **regulator-free** statement is about G (or about lim_{ε→0} of
superbalanced combinations); the exact-at-finite-ε statement holds in the CFT's own conformal regularisation;
FG at O(ε²) is non-universal. A sharper claim — that with any per-geodesic regulator exact membership survives,
because a saturated inequality Q(G) = 0 forces equal multisets of minimising geodesics on both sides (a rigidity
statement for multi-arcs from intersection numbers with all test arcs) — is a **[conjecture]**; it would follow
from inclusion–exclusion over lifted boxes on the disc, not checked in general.

## 2. CFT2 statement

### 2.1 Chiral data of a state

At leading order in c, a "vacuum-sector" state ω of a CFT2 on k circles is characterised by two marked hyperbolic
structures: the holomorphic and antiholomorphic developing maps f^{(j)}_± of each circle (solutions of the Hill
equations ψ'' = (−6/c)⟨T_{±±}⟩ψ, up to the sign convention of L), and the holonomy representations ρ_L, ρ_R of
π₁(S). On a single circle ρ is fixed by the monodromy of ⟨T⟩, ⟨T̄⟩; for several circles (wormhole states) the
relative holonomies are extra data, not visible in one-point functions.

**Hypothesis (V) (vacuum-channel / HRT dominance).** At leading order in c, for every union A of intervals and
whole circles, `S_ω(A) = (c/6) min_{c∈Cut(A)} ½[l_L(c) + l_R(c)]`, with chiral lengths from the developing maps
(twist normalisation). Equivalently: the large-c entropy is HRT in M_ω. In CFT terms this is "one OPE channel
(cut class) dominates, the same for both chiralities, and the vacuum Virasoro block of that channel gives its value".
The common channel for both chiralities is the CFT shadow of the repo's common-class constraint (Theorem H).

### 2.2 Which states satisfy (V), as far as the literature supports

| state / evolution | status of (V) | source |
|---|---|---|
| vacuum, N intervals on one time slice (line or circle) | derived to all orders in a series expansion assuming vacuum-block dominance (Hartman); from handlebody dominance (Faulkner) | 1303.6955, 1303.7221 [abstract] |
| vacuum, intervals on an arbitrary acausal slice (different times) | analytic continuation z → u, z̄ → v of the Euclidean answer (spacelike twists commute); standard but not separately derived for N intervals | [memory] |
| Virasoro coherent states U_{f₊}U_{f₋}|0⟩; any evolution by H(t) = ∫[v₊(t,x)T + v₋(t,x)T̄] (inhomogeneous / Möbius / SSD quenches, Floquet drives) | **exactly equivalent** to the previous row by conformal covariance of twist primaries (Sec. 1.3) | [proved] |
| thermal (rotating, β_L ≠ β_R) state on a circle above Hawking–Page, and its TFD purification, at any times on both sides | single interval incl. plateau at large c (Chen–Wu and others); TFD time evolution matched (Hartman–Maldacena); **multi-interval on the torus not derived from CFT** | 1412.0763, 1303.1080 [abstract]; HKS 1405.5137 for dominance of BTZ in Z [memory] |
| thermal states + Virasoro dressing and inhomogeneous quenches of thermal states (Goto–Nozaki–Ryu–Tamaoka–Tan "black-hole-like excitation") | equivalent to the thermal row by conformal covariance | 2112.14388 [abstract]; [proved] given thermal row |
| multi-boundary states from Euclidean path integrals on Riemann surfaces (wormhole phase) | dual identified in a region of moduli space; time-symmetric (ρ_R = ρ̄_L); rotating ones are bulk states without a standard path-integral preparation | 1406.2663, 0912.2090 [abstract] |
| general Bañados states / states with finitely many heavy insertions, single interval | monodromy method + Wilson lines, vacuum block dominance assumed | Chen–Wu 1605.06753 [abstract + text summary] |
| heavy primary eigenstates (h > c/24), local operator quenches | single interval only; conjectured dominance; dual = BTZ microstate / boosted black hole | ABGH 1410.1392 [abstract] |
| global quenches (Calabrese–Cardy boundary states; collapse by many light operators) | **not covered**: duals contain end-of-world branes or a matter shell (Vaidya) | Hartman–Maldacena; Anous–Hartman–Rovai–Sonner 1603.04856 [abstract] |
| conical-defect primaries (h < c/24) | **not covered** (matter); 2π/n orbifolds via the repo's particles theorem | `particles/` |

Pure heavy states need a modified prescription (no horizon homology): S'(A) = min(F(A), F(A∪O)) where O is the
second boundary of BTZ. This is the min-cut of the same graph with the terminal O turned into a bulk vertex, so it
is again in C_N with the same Bell-pair decomposition [proved for the functional]; that it is the large-c answer for
multiple intervals in a heavy primary is a **[conjecture]**.

### 2.3 The CFT theorem

**Theorem C (vacuum-sector CFT2 states).** Let ω be a state of a CFT2 satisfying (V) whose chiral data are trivial
(vacuum module) or a pair of convex-cocompact Fuchsian representations (thermal, rotating, TFD, wormhole states),
dressed by arbitrary Virasoro transformations on each circle and evolved by arbitrary stress-tensor Hamiltonians.
Take N parties built from intervals on any acausal Cauchy slice of each circle and from whole circles. Then at
leading order in c, in twist normalisation:
1. `S(A) = (c/6)[G(A) + Σ_{p∈∂A} t_p(ε)]` exactly, with `t_p = log(1/ε) + O(1)`.
2. G is the min-cut function of a finite graph; it is UV finite, independent of the regulator and of the conformal
   frame; under a Virasoro dressing it changes only through the uniformised positions of the endpoints.
3. Every superbalanced valid inequality (by He–Hubeny–Rangamani 2002.04558 [abstract], every non-redundant HEC
   inequality other than SA) holds as `Q(S) = Q(G) ≥ 0`, a statement about UV-finite, scheme-independent numbers.
   SA instances hold, with divergent slack when the two regions share endpoints.
4. The closure of the set of rays {G} equals C_N: dynamics generated by the stress tensor, rotation, and wormhole
   topology produce **no** finite entanglement pattern that a time-symmetric (static RT) state does not already
   produce.

**Grade:** [proved mod (V) and the repo's theorems]. (V) itself is a theorem of CFT nowhere; it is the standard
large-c dictionary, derived under vacuum-block-dominance assumptions only for the vacuum (any N), single intervals in
heavy/thermal/Bañados states, and TFD evolution.

**The cleanest formulation (headline form).** *In a holographic CFT2 at leading order in c, the entanglement of
any number of intervals and whole circles in any state of the universal Virasoro sector — vacuum, thermal and
rotating states, thermofield doubles, multi-boundary wormhole states, and all their evolutions under
stress-tensor Hamiltonians (inhomogeneous quenches, Floquet drives), at all times — is a sum of Bell pairs, one per
entangling point with weight (c/6)log(1/ε), plus a UV-finite, frame-independent min-cut vector of a finite graph.
Hence every holographic entropy inequality holds at all times in this sector.*

## 3. The canonical G as a boundary observable

### 3.1 The Bell-pair weight t_p in CFT language [proved]

With the frame lifts of Sec. 1.2, `t_p = ½(t_p^L + t_p^R)`, where for each chirality
`t_p^{L} = log[ K(p₋,p) K(p,p₊) / (ε K(p₋,p₊)) ]`, `K(x,y) = (f(y)−f(x))/√(f'(x)f'(y)) = ψ(y)` with
`ψ'' = Lψ, ψ(x) = 0, ψ'(x) = 1`, L ∝ ⟨T⟩ (`hill_kernel.py`: agreement 1e-12). Here p₋, p₊ are the neighbouring
marked points of p **on the universal cover** of its circle (for a circle with only p and q: q and its translate
around the circle). So t_p is a functional of the stress-tensor one-point function along the arc p₋ → p₊ alone.
Examples (S-units, `S_p = (c/6) t_p`):
- vacuum on the line: `(c/6) log[ a b / (ε (a+b)) ]`, a, b the distances to the neighbours;
- thermal (β; rotating: average over β_±): `(c/6) log[ (β/πε) sinh(πa/β) sinh(πb/β) / sinh(π(a+b)/β) ]`.
Both are exactly **half the mutual information of the two adjacent elementary intervals**, computed with the
connected (boundary-homotopic) saddle. Hence:

**Prop. 3.1.** `t_p = ½·4G·I(α_p⁻ : α_p⁺)` whenever α_p⁻, α_p⁺ and α_p⁻∪α_p⁺ are all in their connected
(non-plateau) phase and share only p (≥ 3 marked points on the circle). In the vacuum module this always holds.
[proved; numerics in `g_boundary.py` (disc): t_CFT = t_can at every point, 1e-14].

### 3.2 Boundary formula for G, simply connected case [proved]

On the disc (vacuum module, any dressing, any acausal slice, |P| ≥ 3):
`G(A) = 4G[ S(A) − ½ Σ_{p∈∂A} I(α_p⁻ : α_p⁺) ]`.
Kinematic-space form: for pieces α = (p_i,p_{i+1}), β = (p_j,p_{j+1}) the Liouville mass of the test class
joining them is `w_{αβ} = ½·4G·I(α:β|γ)` with γ the pieces between them on either side (Czech–Dong's weights,
CLMS 1505.05515); the adjacent boxes w_{α,α±1} are the infinite classes τ_p that the Bell pairs remove. So
`G(interval A) = ½ Σ_{α⊂A, β⊄A, α,β non-adjacent} 4G·I(α:β|γ_{αβ})`, a finite sum of conditional mutual
informations, and for unions `G(A) = min over non-crossing matchings of ∂A of Σ G(chords)` (each endpoint is used once,
so subtracting t_p commutes with the minimum). Numerically G_can ≡ G_CFT on the disc to 7e-15 (`g_boundary.py`).

### 3.3 BTZ and wormholes [proved mod repo]

G(A) = min_c Σ_τ μ'(τ) i(τ,c), τ ranging over test arc classes between lifted pieces; μ'(τ) = ½·4G·I_univ(α̃:β̃|γ̃)
is the CMI of the lifted pieces **in the universal-cover state** (for thermal BTZ in its own frame: the thermal state
on the line). Winding classes contribute ~ e^{−k·2π·(2π/β)} per extra winding. The horizon length is recovered as
the μ'-mass of the arcs crossing it (closed classes have zero mass). This is the multiply connected analogue of
kinematic space run **on the quotient**, the step that fails on the universal cover (GGHH).

### 3.4 A CFT definition valid in every state, and its relation to the canonical G

For any state of any CFT2 and any configuration with ≥ 3 marked points on each circle carrying intervals, define
`G_MI(A) := S(A) − ½ Σ_{p∈∂A} I(α_p⁻ : α_p⁺) b_p(A)` (S-units).
- UV finite and independent of the regulator and of the conformal frame: the endpoint divergence (c/6)log(1/ε)
  and the non-universal endpoint constant enter S(A) once per endpoint, and ½I(α_p⁻:α_p⁺) contains exactly one copy of
  each at p and none at p± (they cancel between S(α_p^±) and S(α_p⁻α_p⁺)). [proved mod locality of the UV structure
  of entanglement entropy in CFT2 (standard)]. This needs no holography and no ⟨T⟩.
- Holographic states: `G_MI = G_can − Σ_p δ_p b_p` with `δ_p = ½I(α_p⁻:α_p⁺) − t_p ≥ 0`.
  Proof of δ_p ≥ 0: G_can ∈ C_N (repo Thm 2(i)), so SA for G_can on the pieces α_p⁻, α_p⁺ gives
  `S(α⁻)+S(α⁺)−S(α⁻α⁺) − 2t_p ≥ 0`. [proved]. δ_p > 0 exactly when α_p⁻∪α_p⁺ (or a piece) is in its plateau
  (horizon-wrapping) phase: in 13 582 sampled BTZ points (3 000 random configurations), δ_p = 0 to 1e-14 in every all-connected case and
  δ_p ∈ [1e-3, 6.5] in all 1 338 plateau cases.
- **Bell-subtraction lemma.** If S ∈ C_n for the fine-grained parties (pieces, whole boundaries, purifier) and
  0 ≤ w_ij ≤ ½ I_S(i:j), then S − Σ w_ij b_ij ∈ C_N for every coarse-graining. Proof: non-SA facets are
  superbalanced (He–Hubeny–Rangamani 2002.04558 [abstract]) hence blind to every b_ij; for SA,
  `S'(X)+S'(Y)−S'(XY) = I_S(X:Y) − 2Σ_{i∈X,j∈Y} w_ij ≥ I_S(X:Y) − Σ_{i∈X,j∈Y} I_S(i:j) ≥ 0`, the last step by
  iterated MMI (superadditivity of MI); pairs not split between X and Y contribute 0; positivity is SA with the
  purifier. [proved mod HHR]. Hence **G_MI ∈ C_N for every state whose fine-grained entropy vector is in the cone**,
  in particular for all states of Theorem C. Checks: 3 000 (N=3) + 800 (N=4) BTZ configurations with δ_p > 0 at
  1 732 points, no violation; 600 random graph models (`adjacent_subtraction_graphs.py`), no violation.

So there are two regulator-free finite parts: G_MI (defined in any CFT state, from entropies only) and G_can
(defined from entropies and ⟨T⟩; canonical geometrically). They coincide in the vacuum module and for "small"
adjacent intervals; both are graph-model vectors in the holographic vacuum sector.

**Is G new?** Honest answer: G is a linear re-coordinatisation of the entropy vector (S minus a Bell-pair vector
whose weights are themselves entropies or ⟨T⟩-functionals); it carries no information beyond S (and ⟨T⟩). Its value is
that it is the UV-finite, scheme- and frame-independent representative, so the cone statement becomes a statement
about well-defined numbers, and "S ∈ C_N for small ε" becomes "G ∈ C_N". Whether such an "adjacent-MI-subtracted
entropy" has been named before: mutual-information regularisations of entropies are standard (Casini–Huerta and
others [memory]); I found no prior use of this specific Bell-pair-subtracted vector as a cone object. [not found ≠
absent].

### 3.5 Contrast: when is G not a graph?

- **Bell-pair-type states (I3 ≡ 0 on all unions of pieces)**: free massless Dirac fermion (fermionic subalgebras),
  quasiparticle-picture quenches. Extensivity gives S(A) = ½ Σ_{α∈A,β∉A} I(α:β), so
  `G_MI(A) = ½ Σ_{α∈A, β∉A, non-adjacent} I(α:β)` with I ≥ 0: a terminal-only graph, always in the cone. [proved].
  Checked for the Dirac fermion on the circle: I3 = 0 to 2e-14, G_MI equals the Bell-pair-sum formula to 2.5e-14,
  all SA/MMI instances satisfied (`g_boundary.py`).
- **Generic CFT2 vacua: G is NOT a graph.** MMI is balanced (superbalanced), so Q_MMI(G) = Q_MMI(S) = −I3. For a
  CFT whose lightest primary is a scalar O with small C_OOO, Agón–Bueno–Casini 2109.09179 eq. (3)/(34) give I3 > 0 at
  long distances (they state the formula is "completely general" for a scalar lightest operator; [read]). Example:
  the **critical Ising CFT**, lightest primary σ (Δ = 1/8), C_σσσ = 0 by Z2: I3 ∼ +r^{−3/4} > 0, so G violates MMI
  and is not in C_3. [proved mod ABC applied at d = 2]. Exact diagonalisation (`ising_i3.py`, periodic chains
  L = 18 and 22): I3 = +0.08 … +0.19 for three blocks of 1–2 spins at every separation.
- **Same CFT, different algebra.** The XX chain (Dirac fermion CFT): for **spin** blocks (bosonic algebra; lightest
  operator σ^± ~ charged vertex operator, C = 0 by U(1)) ED gives I3 = +0.03 … +0.10 > 0; for **fermionic** blocks
  (correlation matrix, infinite chain) I3 → 0 (−4e-2 at gap 1 site, −3e-9 at gap 128). So "G is a graph" is a
  property of the state *and* the subsystem algebra, not of the CFT.
- In holography at leading order, G fails to be a graph only with matter that violates the NEC (repo: Vaidya shell
  violates SSA).

## 4. Literature: is "holographic entropy inequalities hold at all times for vacuum-block-dominated CFT2 states" stated anywhere?

Not found. Details:

| source | mark | what it says relevant here |
|---|---|---|
| Czech–Dong 1905.03787 | [read: full text via extracted HTML, searched] | bulk theorem: every contraction-map inequality holds for HRT in 2+1 bulks with the null curvature condition; Sec. 6 extends to multiboundary via the universal cover (later shown incorrect by GGHH). **No discussion of CFT2 state classes** (no vacuum module, Virasoro, thermal/BTZ, quench); "pure states of CFTs on connected manifolds", mixed via purification. Discussion: SSA "selects the channel in which the identity conformal block dominates". In CFT terms their valid (simply connected) case already gives the contraction-map inequalities for the whole vacuum module at all times. |
| GGHH 2407.07165v2 | [read: pdftotext, searched] | no occurrence of CFT2 / central charge / Virasoro / Bañados; leaves the multiply connected case open. |
| Czech–Lamprou–McCandlish–Sully 1505.05515 | [abstract] | kinematic space: metric from CMI; static slice of AdS3; lengths = kinematic volumes. Not an inequality statement. |
| Bao et al. 1505.07839 | [memory] | static RT cone = graph models. |
| He–Hubeny–Rangamani 2002.04558 | [abstract] | all non-redundant HEC inequalities other than SA are superbalanced (used in Sec. 3.4). |
| Balasubramanian–Kang–Cummings–Murdia–Ross, PRL 136, 031602 (2026), 2509.03621 | [abstract] | inequality on reflected entropy + multipartite signal; "time-symmetric holographic states can never have purely GHZ-like entanglement". Not an entropy-cone statement, restricted to time symmetry; the present theorem does **not** extend it (reflected entropy is not a function of the entropy vector). |
| Hartman 1303.6955; Faulkner 1303.7221 | [abstract] | vacuum, N intervals, large c, vacuum-block / handlebody dominance ⇒ RT. |
| Asplund–Bernamonti–Galli–Hartman 1410.1392 | [abstract] | heavy states and local quenches, single interval, conjectured dominance. |
| Anous–Hartman–Rovai–Sonner 1603.04856 | [abstract] | collapse (Vaidya) from CFT, single interval; matter — outside the vacuum class. |
| Chen–Wu 1605.06753 | [abstract; text summary via tool] | single interval in general states with heavy insertions (Bañados), monodromy + Wilson lines; Hill equation ψ'' + (6/c)Tψ = 0. |
| Roberts 1204.1982 | [abstract] | Poincaré → Bañados diffeomorphism (used in Sec. 1). |
| Goto–Nozaki–Ryu–Tamaoka–Tan 2112.14388 | [abstract] | Möbius/SSD quenches of thermal states, holographic dual = deformed BTZ horizon (in Theorem C's class). |
| Jiang 2611.07010 (thesis); 2511.16729 | [abstract] | SSD as change of foliation; membrane theory of multiparty signals. No cone statement in the abstracts. |
| Agón–Bueno–Casini 2109.09179 | [read] | long-distance I3 > 0 for scalar lightest primary with small OPE coefficient (Sec. 3.5). |
| Casini–Huerta 0812.1773 | [abstract] | extensive MI (I3 = 0) holds for 2D massless free fermions only. |
| Balog–Fehér–Palla hep-th/9703045; Garbarz–Leston 1403.3367 | [abstract] | Virasoro coadjoint orbits = boundary gravitons; energy bounds; exotic orbits. |
| Skenderis–van Rees 0912.2090; Balasubramanian–Hayden–Maloney–Marolf–Ross 1406.2663 | [abstract] | CFT states dual to multiboundary wormholes. |

Conclusion: the CFT2 statement of Theorem C appears to be new as a stated result. Its simply connected,
contraction-map part is implicit in Czech–Dong + Hartman/Faulkner + conformal covariance, but nobody seems to have
written it down. The thermal/rotating/TFD/wormhole part, the non-contraction facets and the finite-part (G)
formulation are not in the literature I could find.

## 5. Grade table

| claim | grade |
|---|---|
| Bañados geometries with hyperbolic winding-0 monodromies (BTZ + gravitons, wormholes + gravitons) are the repo's spacetimes with a moved acausal boundary slice; all times included | [proved mod cited Mess/Barbot classification] |
| Boundary gravitons require extending the theorem | [false]: already included; only the regulator changes |
| FG cutoff z = ε ↔ frame horocycles up to O(ε²), with non-endpoint-local O(ε²) | [proved] (+ numerics) |
| CFT twist-field normalisation = frame horocyclic cutoff, decomposition exact at finite ε | [proved] |
| SM remark "changing horoball normalisation shifts G by a multiple of b_p" | [false] (referee minor 3; G is horocycle-independent) |
| Dressing acts on G only by moving marked points on the ideal circles | [proved] |
| Theorem B (all vacuum solutions with trivial or convex-cocompact holonomy, any dressing, any time) | [proved mod repo + Mess/Barbot; Γ = 1 mod standard disc HRT] |
| Exact cone membership for every per-geodesic regulator (FG) at small ε | [conjecture] (multi-arc rigidity) |
| Virasoro coherent states / stress-tensor evolutions ≡ acausal configurations of the base state | [proved] (any CFT2) |
| Theorem C (CFT statement) | [proved mod hypothesis (V) and the repo]; (V) is derived in CFT only for vacuum (N intervals), single intervals in heavy/thermal/Bañados states, TFD evolution |
| Heavy pure states: functional with O as a bulk vertex is in C_N | [proved] for the functional; its validity as the large-c answer [conjecture] |
| t_p = ½·4G·I(adjacent) | [proved] in the vacuum module and in connected phases; [false] in general (BTZ plateau phases, two-point circles) |
| t_p is a functional of ⟨T⟩ along p₋ → p₊ (Hill kernel) | [proved] |
| Disc: G = S − ½Σ I(adjacent) b_p = CMI sum (kinematic space without adjacent boxes) | [proved] |
| G_MI UV finite, regulator- and frame-independent in any CFT2 state | [proved mod standard locality of EE divergences] |
| δ_p = ½I − t_p ≥ 0, > 0 only in plateau phases | ≥ 0 [proved]; "only in plateau phases" [numerical, 13 582 points] |
| Bell-subtraction lemma ⇒ G_MI ∈ C_N whenever fine-grained S ∈ C_n | [proved mod He–Hubeny–Rangamani superbalance] |
| Bell-pair-type (I3 ≡ 0) states: G is always a graph | [proved] |
| Ising CFT2 vacuum: G not a graph (MMI violated at long distance) | [proved mod ABC formula at d = 2]; ED evidence L = 18, 22 |
| Dirac CFT: fermionic blocks graph, spin blocks not | fermionic [proved]; spin blocks [numerical, ED] |
| "Inequalities hold at all times for vacuum-sector CFT2 states" stated in the literature | [not found] |

## 6. Verdict

**Strongest correct CFT2 headline.** *In any holographic CFT2 at leading order in c, for every state in the
universal Virasoro sector — the vacuum, thermal and rotating states, thermofield doubles and multi-boundary
wormhole states, together with all their evolutions by stress-tensor Hamiltonians (Möbius/SSD/inhomogeneous quenches,
Floquet drives) and at all times — the entanglement entropies of arbitrarily many intervals and whole circles are
exactly one Bell pair of weight (c/6)log(1/ε) per entangling point plus a UV-finite, scheme- and frame-independent
vector that is the min-cut function of a finite graph. Consequently every holographic entropy inequality holds,
and these dynamics generate no entanglement pattern beyond those of static (time-symmetric) states.* Contrast:
the Ising CFT vacuum already violates the graph property (I3 > 0), while free fermions satisfy it trivially.
Conditional on the standard large-c dictionary (V) for multi-interval entropies in these states.

**Does it answer the referee's importance objection (M1)?** Partly, not decisively.
- Gains: a statement in CFT language about a named, physically central class (thermal/rotating states, TFD
  evolution, wormhole states, inhomogeneous quenches) at all times; a regulator-free formulation (G, or superbalanced
  Q's as finite numbers) that disposes of the FG objection (M2c) cleanly via the twist normalisation; a sharp
  diagnostic contrast (graph vs Ising) that is of interest beyond holography.
- Limits a hostile referee will press: (i) the boundary-graviton extension is CFT-elementary (conformal covariance of
  twist primaries) and was already in the bulk theorem, so it adds breadth of wording, not new physics; (ii) the
  sector is exactly the one with no local bulk degrees of freedom — the objection is restated in CFT terms ("states
  chirally uniformisable to static ones"), not removed; the genuinely non-static content is ℓ_L ≠ ℓ_R and the
  common-channel constraint; (iii) the physically most interesting non-equilibrium protocols — global quenches,
  collapse, multi-interval local quenches, heavy microstates — are outside (branes/matter) or conjectural;
  (iv) hypothesis (V) for multi-interval thermal and wormhole states is the holographic dictionary itself, not a CFT
  derivation.
- My estimate: the CFT reformulation raises the Letter's case from "special case" to "universal sector of large-c
  CFT2 dynamics", which is the right framing to lead with, but the answer remains the expected one. PRL acceptance
  stays doubtful (I would put it at ~25–35%); the reformulation plus the G/Ising contrast is a clear win for a
  JHEP/CMP-level paper and for the Letter's introduction.
