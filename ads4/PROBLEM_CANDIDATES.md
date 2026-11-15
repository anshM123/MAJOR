# AdS4/CFT3 programme — open-problem candidates and recommendation (15 November 2026)

Source marks as in `ads4/LITERATURE.md`: **[R]** text read / searched this session, quotations verbatim;
**[A]** abstract only; **[M]** memory, unverified. "Derivation" = physics-level; "theorem" = mathematical
rigour. Scores are my judgement, not measurements.

Bottom line first.

* **Recommended target: P2** — a rigorous theorem that stress-tensor exchange produces a universal *negative*
  shift, of the predicted size, in the mean twist of large-spin two-particle operators of any unitary CFT3
  (long-range gravitational attraction in AdS4 with strength `∝ 1/c_T`).
* **Fallback: P6** — rigorous large-N asymptotics of the ABJM three-sphere partition function
  (first the `N^{3/2}` law with an error bound, then the Airy form).
* **The strip-region entropy-cone idea (P3) is not in the literature I searched, but adds essentially nothing
  beyond Czech–Dong plus Wall plus the GGHH `π_1 = Z` argument.** It ranks ninth of ten on the five criteria
  set by the coordinator. See §1.4.
* The candidates that would matter most for four-dimensional quantum gravity (P1, P8, P4) are not tractable
  with the stated resources.

---

## 1. Requested checks on specific papers

### 1.1 Bousso–Kaya

Two papers match "Bousso and Kaya on holographic entropy inequalities"; neither is dated 2025 *and* about HRT.

**(i) arXiv:2404.11644, "Geometric Quantum States Beyond AdS/CFT" (April 2024) [R].**

* Theorem 4 (static case, verbatim): "Let a, b and c be open subsets of Σ such that
  `a ∩ E(bc) = b ∩ E(ca) = c ∩ E(ab) = ∅`. Then
  `A[E(a)] + A[E(b)] + A[E(c)] − A[E(ab)] − A[E(bc)] − A[E(ca)] + A[E(abc)] ≤ 0`."
  Here Σ is a time-reflection-symmetric Cauchy slice and `E` the generalized entanglement wedge.
* Theorem 29 (time-dependent, verbatim): "(Monogamy of Mutual Information): With x, y, and z as in
  Convention 26, `A[e(xy)] + A[e(yz)] + A[e(zx)] ≥ A[e(x)] + A[e(y)] + A[e(z)] + A[e(xyz)]`."
* Hypotheses: independence of the input wedges (Convention 26), matter entropy neglected, normality and
  weak-focusing properties of generalized entanglement wedges taken from Bousso–Penington
  (their refs. Bousso:2022hlz, Bousso:2023sya).
* What it does **not** prove (their footnote 2, verbatim): "aside from MMI, the holographic entropy cone for
  the time-dependent case (HRT) has not been established even in AdS/CFT."
* Status: theorem-style GR argument for **MMI only**, in the time-dependent setting, for generalized wedges.

**(ii) arXiv:2502.03516, "Holographic Entropy Cone Beyond AdS/CFT" (February 2025) [R].**

* Theorem 12 ("Generalized Static Holographic Entropy Cone", verbatim core): "Let `a_0 = ∅`; and let
  `a_1, …, a_n` be open subsets of the static Cauchy slice Σ such that `a_i ⊂ [E(⋓_{j≠i} a_j)]′` for all i
  … If there exists an α,β contraction map `f : {0,1}^L → {0,1}^R` with `f(x^i) = y^i` for all
  `i ∈ {0,…,n}`, then `Σ_l α_l A[E(v_l)] ≥ Σ_r β_r A[E(w_r)]`."
* Hypotheses: **static** (time-reflection-symmetric slice, their Definition 1, which says the terminology is
  chosen "for forward compatibility with the time-dependent case not considered in this paper");
  mutual independence of the bulk input regions; classical limit.
* Their own statement of the state of the art: "With the exceptions of SSA and MMI, the HEC inequalities have
  so far been proven only in static settings".
* Status: theorem, static only. **Says nothing about HRT.**

### 1.2 Akers–Hernández-Cuenca–Rath, arXiv:2108.07280 (JHEP 11 (2021) 177) [R]

* Setting: entropies given by the QES formula "at a moment of time symmetry". Verbatim: "A generalization to
  nontrivial time-dependent settings will be no less subtle than that of understanding the analogue of the HEC
  for HRT, a problem which remains barely understood except for specific situations such as the
  low-dimensional case studied in Ref. [Czech–Dong]."
* Theorem 2 (verbatim): "If `S_bulk` obeys all HEC inequalities, then `S_QES` obeys all HEC inequalities too."
* Proposition 2 (verbatim): "Bulk MMI implies boundary MMI", via the seven-party bulk inequality
  `S(ABDG) + S(ACEG) + S(BCFG) ≥ S(A) + S(B) + S(C) + S(ABCDEFG)`.
* General fact stated: an n-party inequality holds for `S_QES` if the bulk state obeys an n′-party inequality
  with n′ "generally doubly-exponential in n".
* Explicit conjectures (verbatim, strongest first): "If a given n-party inequality no stronger than HEC ones is
  satisfied in the bulk, then the same n-party inequality is satisfied on the boundary." / "If the complete set
  of n-party HEC inequalities is satisfied in the bulk, then every n-party HEC inequality is satisfied on the
  boundary." / a weaker `n_k` version. "Such a situation may be hard to prove but easy to falsify."
* Status: theorems (contraction-map combinatorics) for a **static** QES cone; conjectures open as far as the
  38 INSPIRE citations show by title.

### 1.3 Bao–Furuya–Naskar

Three papers by these authors (plus one by Bao–Naskar) concern holographic entropy inequalities; all are about
the **RT** (static) cone.

* **arXiv:2506.18086** (JHEP 12 (2025) 140), "On the completeness of contraction map proof method for
  holographic entropy inequalities" [R]. Theorem 4.1 (verbatim): "Given a candidate entropy inequality
  `Σ α_i S_{P_i} ≥ Σ β_j S_{Q_j}`, `α_i, β_j ∈ Z_+`, it is a valid HEI if and only if a contraction map
  `f : {0,1}^M → {0,1}^N` exists (where `M = Σ α_i` and `N = Σ β_j`) satisfying homology conditions on the
  boundary." Definition 4.1: "An entropy inequality is a valid HEI if and only if every holographic geometry
  obeying the RT formula satisfies the inequality." Hypothesis flagged by the authors: the proof is "upto
  certain disclaimers", the first being that it relies "on the existence of a graph model corresponding to a
  smooth bulk geometry (and vice versa)". The words HRT, covariant, time-dependent do not occur in the text.
  Acceptance is not universal: arXiv:2601.09987 v2 restricts its results to "contraction-provable HEIs".
* **arXiv:2409.17317** (JHEP 03 (2025) 117), "Towards a complete classification of holographic entropy
  inequalities" [A; text searched]: algorithm generating all contraction-map inequalities via partial cubes;
  completeness "argued". RT only.
* **arXiv:2408.04741** (JHEP 10 (2024) 251), generalized toric inequalities [A]: conjectures and partial
  proofs. RT only.
* **arXiv:2403.13283** (Bao–Naskar, JHEP 06 (2024) 39) [A]: deterministic construction of contraction maps;
  "an argument for the completeness". RT only.

None of the Bao–Furuya–Naskar papers contains a time-dependence result.

### 1.4 Does a "geodesic current ⇒ contraction-map inequalities" theorem, applied to strips in planar-symmetric time-dependent AdS4, give anything new?

**Set-up.** Planar-symmetric asymptotically AdS4 spacetime, regions = finite unions of parallel strips whose
edges lie on a common boundary Cauchy slice. Translation-invariant extremal surfaces are curves in the 2+1
quotient, extremal for the area functional (length weighted by the transverse warp factor). The reduced
spatial slice has `π_1` trivial (one-sided collapse, or two-sided brane with non-compact strip direction) or
`Z` (strip direction compactified, two-sided).

**What the existing results give.**

1. *Simply connected.* Czech–Dong's proof [R] is the identity (their eq. (15))
   `Σ α_l S(I_l) = Σ β_r S(J_r) + ½ Σ_{i≠j} (h_LHS − h_RHS)(Z_i,Z_j) · CMI(i,j)`
   together with a combinatorial lemma that a contraction map gives `h_LHS ≥ h_RHS`. The only physical inputs
   are (a) the entropy of a union of intervals is a minimum over pairings of sums of *single-interval*
   entropies, and (b) strong subadditivity of single-interval entropies, "which the RT and HRT proposals are
   known to obey". Both inputs are available for strips: (a) by symmetry reduction, given a symmetric HRT
   surface (GGHH Corollary B.3 covers finite symmetry groups; the continuous-translation, per-unit-area
   version needs a sentence); (b) from Wall's maximin SSA applied in four dimensions to adjacent strips.
   So the strip statement follows from the published argument essentially verbatim. Czech–Dong's remark that
   the proof "is unlikely to extend to higher-dimensional cases" refers to general regions, for which there is
   no ordering; strips have one.
2. *`π_1 = Z`.* GGHH §3.1 [R] is **not** restricted to vacuum. Hypotheses, verbatim: "Let M be a 2+1
   dimensional spacetime obeying standard classical holographic assumptions (asymptotically locally AdS,
   globally hyperbolic, obeying the Einstein equation and null energy condition), with spatial topology of a
   cylinder." And: "Here we are not restricting ourselves to vacuum solutions". The paper presents this as
   "an argument", not a theorem: it uses `S(nA) = n S(A)` on n-fold covers, Czech–Dong on the universal
   cover, and the estimate `lim_{n→∞} [S(ñA) − S(nA)]/n = 0` (their eq. (36)), about which they write:
   "(This statement seems intuitively clear to us, but we have not attempted to prove it formally.)"
   For other topologies "the estimate (36) does not obviously carry over". Their *numerics* are vacuum only.
   The same argument transfers to strips on a compactified planar-symmetric AdS4 spacetime with the same gap.
3. *The three papers of §1.1–1.3.* None implies or bears on the strip statement: Bousso–Kaya is static for
   the full cone (time-dependent only for MMI, and for generalized wedges); AHR is a static QES cone;
   Bao–Furuya–Naskar is a statement about which inequalities are RT-valid.

**Answer.** A positive-measure ("geodesic current") theorem applied to AdS4 strips would give:

* for trivial `π_1`: **nothing new in substance** — a restatement of Czech–Dong's identity (the positive
  measure is the collection of conditional mutual informations) with the observation that it applies to
  strips. Possible small additions: an explicit weighted planar graph realising the entropies (hence all
  RT-valid inequalities, without appeal to arXiv:2506.18086); the spherically symmetric analogue with
  latitude bands, where the rotation axis acts as one extra boundary point of the reduced disc and
  cap-type surfaces must be included (my observation, unpublished as far as I found, and easy);
* for `π_1 = Z`: a **rigorous replacement for GGHH's unproved estimate (36)**, i.e. a proof where the
  literature has an argument. This is a tidy-up, not a new physical statement;
* in both cases the hypothesis "lengths are integrals against a positive measure" must itself be derived
  from SSA on the (universal cover of the) reduced spacetime, which is exactly the input Czech–Dong and GGHH
  already use. The measure formulation adds real content only for reduced topologies beyond `Z`, which do
  not arise for strips in planar-symmetric black-brane geometries.

**Pure AdS4.** Whether the RT inequalities beyond SSA/MMI hold for HRT in vacuum AdS4 with regions on an
arbitrary non-time-symmetric boundary Cauchy slice is **open**. Evidence for the open status:
Grimaldi–Headrick (arXiv:2610.33905, Oct 2026) [R]: "The validity of higher HEIs for the HRT formula remains
an open problem"; Bousso–Kaya footnote quoted above; GGHH list the known cases as Wall (SSA, MMI),
Bao–Mezei (membrane regime), Czech–Dong (2+1, simply connected). Known sub-cases in pure AdS4: slices that
are images of a time-symmetric slice under a bulk isometry (trivially RT); neighbourhoods of light-cone
configurations to the order tested by the majorization theorems (arXiv:2508.21823, 2601.09987, 2601.09989);
Zhao's conditional theorem (arXiv:2602.04888) if its exposed-region hypothesis holds, which has been verified
for no class. I found no result special to vacuum AdS4.

---

## 2. Candidates

Each entry: statement; why it matters; what is known; nearest proved result and the gap; strategy;
computation; risk of prior solution (with the search done); tractability T (1–5, for a small team with
pen-and-paper and python/mpmath/sympy); significance S (1–5).

### P1. HPPS for the stress tensor in d = 3 (bootstrap / dispersive)

* **Statement.** Let `{C_n}` be unitary CFT3s with unique stress tensor, `c_T → ∞` with large-N
  factorisation, in which every single-trace primary of spin > 2 has dimension `≥ Δ_gap`. Then the two
  non-Einstein parameters `(t_4, α_T)` of ⟨TTT⟩ obey `|t_4| + |α_T| ≤ K (log Δ_gap)^a / Δ_gap^b` with
  explicit universal `K, a, b` (dimensional analysis suggests `b = 4`, from Riemann³ terms in 4D [M]).
* **Why it matters.** This *is* the emergence of Einstein gravity in AdS4 from CFT3 data.
* **Known.** Disc bound on `(t_4, α_T)` (arXiv:1707.03007, derivation); d = 4 analogues (1610.09378,
  1904.05892, 1904.05905, derivations); 4D flat-space graviton bounds with an IR logarithm (2201.06602);
  scalar-correlator AdS4 bound in which AdS regulates the logarithm (2106.10274).
* **Nearest result and gap.** Chang–Landau–Simmons-Duffin (2311.04271) supply the bulk-point dictionary for
  spinning sum rules and state it is "a crucial step for establishing the HPPS conjecture for stress tensor
  correlators"; missing by their own account: the Regge regime for spinning correlators, light-state
  contributions, sum rules for the remaining polarisations.
* **Strategy.** Uplift the 4D functionals of 2201.06602 with the 2311.04271 dictionary; prove positivity in
  the Regge regime using d = 3 spinning blocks.
* **Computation.** Spinning conformal blocks and semidefinite optimisation at a scale beyond the stated
  resources.
* **Risk already solved.** Low (INSPIRE citations of 2311.04271 checked by title, to Nov 2026).
* **T = 1, S = 5.**

### P2. Rigorous universal large-spin attraction from stress-tensor exchange (rigorous lightcone bootstrap)

* **Setting.** d = 3 (the statement makes sense for d ≥ 3). `G(z, z̄) = Σ p_{τ,ℓ} g_{τ,ℓ}(z, z̄)`, the
  four-point function of identical scalars of dimension `Δ`, with `p ≥ 0`, discrete spectrum, crossing
  `((1−z)(1−z̄))^Δ G(z, z̄) = (z z̄)^Δ G(1−z, 1−z̄)`, blocks normalised so that
  `g_{τ,ℓ}(z, z̄) = z^{τ/2}[k_{τ+2ℓ}(z̄) + O(z)]`, `k_{2h}(x) = x^h ₂F₁(h, h; 2h; x)`.
  Hypothesis (H): apart from the identity, the operators of twist below `τ_* + δ` (some `δ > 0`) are finitely
  many, all of twist exactly `τ_*`; put `b_* := Σ_{τ_O = τ_*} p_O Γ(τ_* + 2ℓ_O) / Γ(τ_*/2 + ℓ_O)² > 0`.
  For `0 < z < 1` define the positive measure on twists
  `μ_z := z^Δ Σ_{τ,ℓ} p_{τ,ℓ} k_{τ+2ℓ}(1−z) δ_τ` (it weights spins `ℓ ≲ z^{-1/2}`).
* **Statement.** Under (H), and for `Δ` and `δ` in a range to be determined (my heuristic estimates suggest
  `Δ > τ_*/(2 − τ_*)` suffices for a first theorem), for every fixed `0 < η < 1`:
  `∫_{|τ − 2Δ| < η} (τ − 2Δ) dμ_z(τ) = −2 b_* z^{τ_*/2} (1 + o(1))` as `z → 0`.
  In particular the leading double-twist family contains operators with `τ < 2Δ` at arbitrarily large spin,
  and its spin-averaged anomalous dimension is `γ_0 ℓ^{−τ_*}` with the coefficient of
  Fitzpatrick–Kaplan–Poland–Simmons-Duffin, eq. (43) [R]. For a CFT3 in which the stress tensor is the unique
  operator of minimal twist `τ_* = 1` in `φ × φ`, `b_* = p_T · 128/(3π)` with `p_T ∝ Δ²/c_T`
  (normalisation of `p_T` not verified here): the binding energy of two particles in AdS4 at large separation
  is negative, universal, and proportional to `G_N ∝ 1/c_T`.
* **Why it matters.** It is the theorem-level form of "gravity is attractive and universal at long distance in
  the emergent AdS4": a necessary property of any gravitational bulk, derived from CFT3 axioms alone, with the
  Newton constant identified with `1/c_T`.
* **Known.** Derivations: arXiv:1212.3616, 1212.4103, 1403.6829; Caron-Huot's inversion formula
  (1703.00278) gives large-spin expansions "with controlled errors" but is not written at theorem level.
  Theorems: twist accumulation (Pal–Qiao–Rychkov, 2212.04893); spectrum and OPE coefficients tend to GFF
  for all towers (van Rees, 2412.06907).
* **Nearest result and gap.** van Rees proves the `z → 0` limit of `μ_z` (and its analogues for higher
  towers). Missing: the first correction. His outlook: "It will be interesting to understand the size of the
  corrections to the results presented here."
* **Strategy.** (i) Crossing at small `z` gives, for `w = log(1 − z̄) < 0`,
  `∫ e^{(τ/2−Δ)w} dμ_z = (1−e^w)^{−Δ}[1 + z^{τ_*/2} Σ_O p_O k_{τ_*+2ℓ_O}(1−e^w) + …]`, and
  `k_{2h}(1−y) = −(Γ(2h)/Γ(h)²)(log y + 2ψ(h) − 2ψ(1)) + O(y log y)`: the term linear in `w` at the bottom
  tower is `−b_* z^{τ_*/2} w` (my computation; consistent with eq. (43) of 1212.3616 after converting the
  `z`-weight to spin via `∫_0^∞ t^{2a−1} K_0(2t) dt = Γ(a)²/4`). (ii) A positive measure whose Laplace
  transform is `1 + ε(a + bw) + o(ε)` on an interval has first moment `bε + o(ε)` and variance `o(ε)`
  (second differences in `w` are non-negative). (iii) The difficulty: the higher towers `τ ≈ 2Δ + 2n`
  contribute `e^{nw}` at order one, so isolating the bottom tower to accuracy `o(z^{τ_*/2})` requires
  `1 − z̄ ≪ z^{τ_*/2}`, a double-lightcone window in which the direct-channel remainder must be bounded
  uniformly, in the manner of Pal–Qiao–Rychkov's "DLC" estimates. My crude bound closes the window only
  under extra gap assumptions; removing them is the research content.
* **Computation.** mpmath checks on generalized free fields, on mean-field theory plus one exchange, and on
  published large-N or numerical-bootstrap data; no cluster.
* **Risk already solved.** Medium. Not found in the 23 papers citing 2412.06907 nor the 51 citing 2212.04893
  (INSPIRE, titles; abstracts of the closest three read). van Rees states an intention to continue:
  "we plan to use this result as a stepping stone towards proving several of the claims". Physicists regard
  the *result* as known; the novelty is the proof.
* **T = 3, S = 4.**

### P3. Covariant entropy cone for symmetric region families in AdS4 (entropy cone / HRT)

* **Statement.** Let `(M, g)` be asymptotically AdS4, globally hyperbolic in the AdS sense, with a
  hypersurface-orthogonal spacelike translation symmetry (planar case) or SO(2) rotation symmetry with
  axis (spherical case), obeying the null curvature condition, and let regions be finite unions of parallel
  strips (resp. latitude bands) with edges on one boundary Cauchy slice. If the reduced spatial slice is a
  disc, the HRT entropy vector (per unit transverse length in the planar case) equals the min-cut vector of
  an explicit weighted planar graph, hence satisfies every RT inequality.
* **Why it matters.** Would be the first statement of the full cone for genuinely time-dependent states in
  bulk dimension four.
* **Known / gap.** See §1.4: for strips this is Czech–Dong's identity plus Wall's SSA; for `π_1 = Z` it is
  GGHH's argument with one unproved estimate. Numerics for SSA in a charged Vaidya-type background:
  1304.3398 [A] (that it uses strips in higher dimension is from memory [M]).
  Missing from the literature: the written statement; the band case with caps; a proof of GGHH's (36).
* **Strategy.** Chord-arrangement graph from the non-negative second differences of single-strip entropies;
  for bands add one boundary vertex for the axis and check the crossing inequalities involving it
  (they are SSA and weak monotonicity of actual HRT entropies).
* **Computation.** Optional Vaidya-AdS4 shooting for illustration.
* **Risk already solved.** Not found (full-text searches listed in LITERATURE.md §(f); 43 citing papers of
  Czech–Dong by title). Risk that referees call it a corollary: high, and justified.
* **T = 5, S = 2** (S = 1 for emergent Einstein gravity specifically).

### P4. Covariant entropy cone in pure AdS4 on arbitrary boundary slices (entropy cone / gravitational constraints)

* **Statement.** For vacuum AdS4 and any regions on any spacelike boundary Cauchy slice, HRT entropies
  satisfy every contraction-map inequality. First test case: the five-region cyclic inequality for unions
  of boosted discs.
* **Why it matters.** Simplest setting in which the higher-dimensional HRT cone is unknown; a proof would be a
  new Lorentzian integral-geometry theorem, a counterexample a major surprise.
* **Known / gap.** Open (§1.4). No ordering structure, no chiral decomposition, no known positive Crofton
  formula for codimension-two spacelike surfaces in AdS4.
* **Strategy.** None with a credible path: minimax "cooperating time-sheets" (2502.09894) or verification of
  Zhao's exposed-region condition (2602.04888) in vacuum.
* **Computation.** Numerical extremal surfaces for a falsification search are feasible but not a proof.
* **Risk already solved.** Low.
* **T = 1–2, S = 3.**

### P5. Bulk-to-boundary entropy-cone conjecture for quantum extremal surfaces (QEC / QES)

* **Statement (AHR conjecture, verbatim).** "If the complete set of n-party HEC inequalities is satisfied in
  the bulk, then every n-party HEC inequality is satisfied on the boundary." First open case `n = 5`; the
  stronger inequality-by-inequality form is also open.
* **Why it matters.** Characterises which boundary entropy vectors are compatible with a semiclassical bulk
  carrying quantum matter: a necessary condition one order beyond classical gravity.
* **Known.** Theorem: bulk HEC for all party numbers ⇒ boundary HEC; Proposition: bulk MMI ⇒ boundary MMI
  (arXiv:2108.07280). Static setting.
* **Gap.** The n-party-to-n-party statement for `n ≥ 5`.
* **Strategy.** Enumerate QES configurations as bitstring labellings; for each five-party inequality search
  for a certificate (bulk five-party inequalities + SA/SSA implying the required many-party bulk inequality)
  by linear programming, or for a counterexample entropy vector. The authors call it "easy to falsify".
* **Computation.** Linear programming over entropy vectors; modest.
* **Risk already solved.** Low–medium (38 citing papers, none by title).
* **T = 4, S = 2.**

### P6. Rigorous large-N asymptotics of the ABJM sphere partition function (mathematical QFT / nonperturbative data)

* **Statement.** For integer `k ≥ 1` let
  `Z(N,k) = (1/N!) Σ_{σ∈S_N} (−1)^σ ∫ dᴺx Π_i ρ_k(x_i, x_{σ(i)})` with `ρ_k` as in Mariño–Putrov eq. (3.3).
  Stage 1: `log Z(N,k) = −(π√(2k)/3) N^{3/2} + O(N)`. Stage 2: error `O(N^{1/2})` with the coefficient
  fixed by `B = k/24 + 1/(3k)`, and the `−¼ log N` term. Stage 3:
  `Z(N,k) = e^{A(k)} C^{−1/3} Ai[C^{−1/3}(N − B)] (1 + O(e^{−c√N}))`, `C = 2/(π²k)`.
* **Why it matters.** `Z_{S³}` is the one explicit nonperturbative observable of an AdS4/CFT3 pair; the
  `N^{3/2}` coefficient is `L²/G_N`, and the logarithm is the one-loop quantum-gravity test of 1210.6057.
* **Known.** All three stages are derivations or conjectures (1007.3837, 1106.4631, 1110.4066, 1306.1734);
  Boldis–Korchemsky–Testa (2602.10196): "their derivation remains largely conjectural."
* **Nearest result and gap.** Weyl law for related functional-difference operators
  (Laptev–Schimmer–Takhtajan, 1510.00045, theorem). No theorem for the ABJM operator or for `Z(N,k)` found.
* **Strategy.** `Ξ(z) = det(1 + zρ_k) = Π(1 + zλ_n)` with `ρ_k` positive trace class. (i) `Z(N)` is the
  N-th elementary symmetric function of the `λ_n`, so `z^N Z(N)/Ξ(z)` is the law of a sum of independent
  Bernoulli variables; a local limit theorem gives
  `log Z(N) = J(μ_*) − μ_* N − ½ log(2π J″(μ_*)) + o(1)` with `J = log Ξ(e^μ)`, `J′(μ_*) = N`.
  (ii) Leading `J(μ) = Cμ³/3 + O(μ²)` from trace inequalities: `Tr ρ^s` versus the exact classical value of
  `Tr e^{−sU/2}e^{−sT}e^{−sU/2}` (Araki–Lieb–Thirring), log-majorisation, and eigenvalue comparison with
  `(|x|+|p|)/2`. (iii) Stage 2–3 need `Tr ρ^s` near `s = 0` to higher accuracy (spectral zeta function,
  cf. 1503.07883): hard.
* **Computation.** Exact `Z(N,k)` at small N are published (1207.4283, 1207.5066) for checks; mpmath.
* **Risk already solved.** Medium: possibly in the spectral-theory literature under another name
  (arXiv searches for ABJM with "rigorous", "Weyl law", "mathematically" found nothing).
* **T = 4 (stage 1), 2 (stage 3); S = 3.**

### P7. Sharp distal split property for the N = ∞ single-trace net in AdS4 (operator algebras)

* **Statement.** For a countable family of free fields on AdS4 with boundary dimensions `Δ_i` (the `N = ∞`
  single-trace sector), the inclusion of causal-wedge algebras of nested boundary diamonds is split iff the
  separation exceeds an explicit function of `β_H = inf{β : Σ_i e^{−βΔ_i} < ∞}`; in particular it is split
  for all separations iff the spectrum is sub-Hagedorn. Extension: the black-hole conjecture of
  arXiv:2510.01556.
* **Why it matters.** A candidate algebraic criterion for bulk locality (Faulkner–Li conjecture:
  "split property for 𝒞 ⟺ local quantum field theory").
* **Known.** Chiral theorem (math-ph/0603083); free-string model in Minkowski and AdS (2510.01556,
  derivation).
* **Gap.** Theorem-level statement with sharp threshold in AdS4; the black-hole case.
* **Strategy.** Modular nuclearity for second-quantised nets; conformal covariance of generalized free fields.
* **Risk already solved.** High for the qualitative claim (2510.01556), medium for a sharp theorem.
* **T = 3, S = 2.**

### P8. Penrose inequality for asymptotically AdS4 initial data (gravitational constraints)

* **Statement.** For asymptotically hyperbolic time-symmetric data with scalar curvature `≥ −6` and outermost
  minimal boundary of area `A`, the mass is bounded below by the Schwarzschild–AdS value at area `A`.
* **Why it matters.** Implied by AdS/CFT (Engelhardt–Horowitz, arXiv:1903.00555: "a holographic argument in
  favor of the AdS Penrose inequality"); a failure would contradict the duality.
* **Known.** Special classes by arXiv number and title only: graphs (1201.3321, 1304.7887, 1309.6248),
  perturbations of Schwarzschild–AdS (1402.4317, 2209.01234), locally hyperbolic with non-positive mass
  (1310.3002), cohomogeneity one (2404.13247), suboptimal constant in the hyperboloidal regime (2504.10641),
  planar/spherical with NEC (2204.12239) [all A, titles and first lines].
* **Gap.** The general case.
* **Risk.** Heavily worked by geometric analysts; the repo's own LITERATURE.md records AI-generated
  manuscripts on AdS Penrose inequalities (not re-verified).
* **T = 1, S = 3.**

### P9. Positivity lemma behind the AdS4 scalar locality bound (bootstrap, lemma level)

* **Statement.** For d = 3 let `Q_n(β) = ∫_0^∞ dν ν^n P_{−1/2+iν}(cosh β)` (conical Legendre function,
  regularised). Prove `Q_1(β) > 0` for all `β > 0` (and classify the `(n, d)` with `Q_n > 0`).
* **Why it matters.** arXiv:2106.10274 uses it for the AdS4 bound (236) and states it on numerical grounds:
  "By checking numerical examples, we find that Q_n(β) is everywhere positive … for 0 ≤ n ≤ 1 in d = 3".
* **Known (this audit, numerics only).** Using the Mehler–Dirichlet representation I get
  `Q_0 = 1/(2 sinh(β/2))` and
  `Q_1 = (√2/π)[g(0)/β − ∫_0^β (g(t) − g(0)) t^{−2} dt]`, `g(t) = (cosh β − cosh t)^{−1/2}`;
  mpmath gives `Q_1 > 0` on `β ∈ [0.02, 20]`, `Q_1(0⁺) ≈ 0.041667` (consistent with 1/24), and the stated
  large-β behaviour. Script (session scratchpad, may not persist):
  `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/tools/q1.py`.
  This is evidence, not a proof, and the representation itself should be re-derived before use.
* **Gap.** An analytic proof; then the remaining asymptotic steps of 2106.10274 would still be physics-level.
* **T = 5, S = 2.**

### P10. Covariant entropy cone for 2+1 spacetimes with matter and arbitrary topology (repo item G2)

* **Statement.** As in `OPEN_GAPS.md` G2. Listed for completeness: it is a three-dimensional statement and
  its four-dimensional corollaries are those of P3 only.
* **Status.** Open by GGHH ("remains an open question") beyond `π_1 = Z`. **T = 3, S = 1 for the 4D
  objective** (higher for the 2+1 programme).

---

## 3. Ranking on the five prescribed criteria

Criteria, in priority order: **C1** fundamental significance for four-dimensional quantum gravity;
**C2** genuine mathematical novelty; **C3** direct implications for emergent Einstein gravity or
nonperturbative consistency; **C4** availability of a precise, independently verifiable theorem;
**C5** likelihood of an important physical conclusion rather than a technical extension.
Scores 1–5; T = tractability with the stated resources (not one of the five criteria, used only to choose
the target).

| Rank | Candidate | C1 | C2 | C3 | C4 | C5 | T |
|---|---|---|---|---|---|---|---|
| 1 | P1 HPPS for stress tensors, d = 3 | 5 | 5 | 5 | 3 | 5 | 1 |
| 2 | **P2 universal large-spin attraction (rigorous)** | 4 | 3 | 4 | 4 | 3 | 3 |
| 3 | **P6 ABJM partition function asymptotics** | 3.5 | 3 (stage 3: 4) | 4 | 5 | 3 | 4 / 2 |
| 4 | P8 AdS4 Penrose inequality | 3 | 5 | 3 | 5 | 3 | 1 |
| 5 | P4 pure-AdS4 covariant cone | 3 | 5 | 2 | 4 | 3 | 1–2 |
| 6 | P5 AHR quantum-extremal-surface cone conjecture | 2 | 3 | 2 | 5 | 2 | 4 |
| 7 | P9 positivity lemma for the AdS4 functional | 2 | 2 | 3 | 5 | 1 | 5 |
| 8 | P7 sharp distal split in AdS4 | 2 | 2 | 3 | 3 | 2 | 3 |
| 9 | P3 strips/bands covariant cone in AdS4 | 2 | 1 | 1 | 5 | 1 | 5 |
| 10 | P10 2+1 with matter, all topologies | 1 | 4 | 1 | 4 | 1 | 3 |

Remarks on the ranking.

* P3 is placed ninth on the merits: it is a corollary-level statement (C2 = 1), its proof uses Einstein
  dynamics only through Wall's SSA (C3 = 1), and it yields no new physical conclusion (C5 = 1).
  Its resemblance to the completed AdS3 work has been given no weight.
* No tractable candidate scores above 3 on C5. P2 and P6 turn believed physics into theorems; they do not
  change what is believed. The candidates that could (P1, P4, P8) have T ≤ 2.

### Recommendation

**Target: P2.** Highest-ranked candidate with T ≥ 3. Exact statement as in §2 P2 ("Statement"). Reasons:
it concerns the gravitational interaction itself, holds for every unitary CFT3 satisfying (H) (no large-N
assumption), has a sharp, numerically checkable constant, builds on two recent theorem-level papers, and is
explicitly left open by the more recent of them. Honest weaknesses: physicists consider the result known;
the first provable version may need restrictive hypotheses on `Δ` and on the twist gap; the group that
proved the leading-order theorem has announced follow-up work.

**Fallback: P6, stage 1** (then stage 2). Exact statement: for each integer `k ≥ 1`,
`log Z(N,k) = −(π√(2k)/3) N^{3/2} + O(N)` for the integral `Z(N,k)` defined in §2 P6. Reasons: fully
explicit object, a proof strategy with standard tools, exact small-N data for verification, and direct
contact with the bulk Newton constant. Weakness: stage 1 confirms a number nobody doubts; the physically
informative part (Airy function, `A(k)`, instantons) is stage 3 with T = 2.

Not recommended as targets: P3 (for the reasons in §1.4), P9 (worth doing as a by-product: one page if the
inequality yields), P5 (tractable, but about the static cone with quantum corrections rather than about
four-dimensional gravity).

---

## 4. What I could not verify

1. **Completeness of searches.** arXiv abstract queries, INSPIRE citation lists read by title, full-text
   pattern searches of about 25 papers. A statement of P3, a proof of P2 or of P6 stage 1 could exist in a
   paper these searches missed (in particular in the mathematics literature for P6).
2. **My own sketches.** The first-moment mechanism and the window condition `Δ > τ_*/(2 − τ_*)` in P2; the
   Bernoulli/local-limit and trace-inequality route in P6; the axis-vertex construction for latitude bands in
   P3; the formula for `Q_1` in P9. None has been written out or checked by anyone else.
3. **Normalisations.** `p_T` in terms of `Δ` and `c_T`; the exponent `b` in P1; `A(k)` in P6.
4. **Items quoted from memory [M].** The identification of `(t_4, α_T)` with Riemann³ couplings; the classical
   split/nuclearity theorems for free fields in Minkowski space (Doplicher–Longo, Buchholz–Wichmann);
   attribution of authors to the Penrose-inequality papers in P8 (only arXiv numbers and titles were seen).
5. **Level of rigour of "proofs" in physics papers** (Maldacena–Zhiboedov; Hofman et al.; Kologlu et al.;
   Furuya et al.; Engelhardt–Liu 2509.05413; Zhao 2509.23119): classified from abstracts, not audited.
6. **Papers known by abstract or title only**: all entries marked [A] in LITERATURE.md; in particular
   Wall's precise hypotheses (taken from the abstract and from GGHH's description), Zhao 2602.04888 beyond
   pattern searches, and the 2025–2026 follow-ups 2512.06815, 2604.22612, 2610.00670.
7. **AI-assisted items** listed in `/home/user/MAJOR/LITERATURE.md` §1 (non-arXiv identifiers, vendor
   releases): not re-checked.
8. **GGHH Corollary B.3** is stated for finite symmetry groups; its use for continuous translations per unit
   transverse area (needed for strips) is my extrapolation.
9. **The file `ads4/EXPLORATORY_MEASURE.md`** was present in the directory when I wrote these files; I did not
   create, read or modify it.
