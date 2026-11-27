# math_transfer: existing theorems for the open mathematical steps (P1–P4)

27 Nov 2026, research agent. Nothing outside `ads4/code/math_transfer/` was edited, and nothing was committed.
Context read: `prl/sec_transport.tex` (Remarks 1–10), `code/cft_etas/THEOREM_C.md` (the ledger the task calls
`cft_etas/README.md`; there is no README.md in that folder), and `code/m1_remainder/README.md`.

**Grades.** [T] a published theorem whose statement I checked against the source (the source is named).
[T?] a published theorem stated from standard knowledge; I did not see the source text, so the theorem number or
exact form is unverified. [P] proved here, with the short proof given. [D] derived here. [NE] numerical evidence from
the scripts in this folder. I give no theorem numbers I could not check.

---

## Ranked list of transfers (most useful first)

1. **P1: Theorem C is classical, and no method can beat it** [T + P]. In the Mellin variable, Theorem C is optimal
   recovery in H^∞ of a half plane. Its exact-data factor β_N = ∏ n/(n+2γ) is *exactly* the Müntz–Szász product
   (checked to 60 digits, `out_p1_muntz.txt`). The Euclidean data E(n) are Hausdorff moments of the relaxation
   function, so P1 is the classical "extrapolate the moment problem to order −1" problem. By Smolyak's lemma
   (optimal recovery theory), the sup in Prop. C′ is the smallest worst-case error of *any* algorithm, linear or
   not. With noisy data this holds up to a factor of at most 2. So the log rate is final for the class (S)+(E). The paper can cite the classical pieces rather than present
   Theorem C as new: Müntz 1914 / Szász 1916, Blaschke products, Micchelli–Rivlin 1977, Baym–Mermin 1961,
   Cuniberti–De Micheli–Viano 2001.
2. **P1: a ladder of hypotheses with sharp rates** [T / P / D / NE]. Let L = log(M/ε). In P1 only the shape of the
   analyticity domain *at s = ∞*, where the Matsubara points accumulate, sets the rate. Analyticity near 0 does not:
   * half plane: L^{−2γ} (Theorem C);
   * sector of opening 2α > π: exp(−Θ(L^{1−π/2α})) [D, NE];
   * completely monotone relaxation (a Stieltjes class with a sum-rule TV bound): exp(−Θ(√L)) [NE];
   * one-signed deviation, data on a continuum, or analyticity at s = ∞: Hölder ε^θ [P/T].

   A power law in Δ_g needs one of the Hölder rows. Euclidean CFT data supply none of them.
3. **P1: positivity of the relaxation function gives an M-free *lower* bound on η/s from Euclidean data** [P + NE].
   If Φ ≥ 0, then E(n) is a Hausdorff moment sequence, and the truncated Hausdorff problem (an LP) gives
   4πη/s ≥ LB_N(E(1..N)). For example LB₂ = E(1)²/E(2). On the Einstein values this gives 4πη/s ≥ 0.715 (N = 2),
   0.897 (N = 4) and 0.982 (N = 8), and stays at 0.968 with ε = 10⁻⁶ data noise. The caveat: Einstein's own Φ_E
   changes sign at late times (complex QNMs), so the hypothesis must be relaxed to "Φ ≥ −φ₋ with a known φ₋".
   This is the only transfer that produces a usable KSS-direction inequality from Euclidean data.
4. **P2: Chebyshev–Markov–Krein / Isii duality settles which moment bounds exist** [T?/P]. A bound on Σcᵢf(mᵢ) from
   bounds on Σcᵢu_k(mᵢ) exists iff f ≤ Σa_k u_k pointwise on [m_gap, ∞) for some a_k ≥ 0. This is pointwise
   domination; f need not be a positive combination. Σcᵢ/mᵢ² is therefore **not** bounded by any finite set of
   higher inverse moments (m⁻⁴, m⁻⁶, …), whatever the gap [P]. Explicit optimal formulas are given below and checked
   by LP (`out_p2_cmk.txt`).
5. **P3: the support lemma is the equality case of a sharp, classical Cauchy–Schwarz statement, and no cone
   restriction helps unless it bounds the coupling ratio** [T + P]. A linear interference term is bounded iff the
   generating vectors are confined, |g₃/g₂| ≤ κ on the support. Remark 4's scalar family is a unitary single exchange,
   so it lies in *every* crossing-compatible cone. Without a six-point diagonal or a bound on the XTTT/XTT ratio the
   step is impossible. Precedents for confinement: null constraints (Caron-Huot–Van Duong 2021), spectrahedral cones
   (Li–Xu–Yang–Zhang–Zhou 2021), interference collider bounds (Córdova–Maldacena–Turiaci 2017).
6. **P4: no cancellation identity exists, and Mehler–Dirichlet proves it** [T]. The back-to-back point z = −1 is
   X = 0 for P_{−1/2+iν}(X). Mehler–Dirichlet (DLMF 14.20.9) writes this kernel as a *positive* superposition of
   cosh(νφ), φ ∈ [0, π/2]. H_b(−1) is therefore a positive average of the impact-parameter transform at imaginary
   impact parameter iφ, and its size e^{πΛ/2}/√(2πΛ) is fixed by Paley–Wiener growth. Any cancellation has to come
   from crossing (the s↔u image of the same block), not from a special-function identity.

---

## P1. Stability of continuation from Matsubara points

### P1.0 Dictionary (Theorem C in classical language) [P]

Use units with 2πT = 1. THEOREM_C.md defines E(s) = c∫₀^∞ Φ(t)e^{−st}dt. Substituting x = e^{−t} gives

  E(s) = c ∫₀¹ x^{s−1} Φ(−log x) dx = ∫₀¹ x^s dν(x),   dν = cΦ(−log x) dx/x.

So the Euclidean data E(1), E(2), … are the **Hausdorff moments** ∫x^n dν, and 4πη/s = E(0) is the zeroth
moment. The data start at n = 1, so mass near x = 0 (late times) is invisible to them. Analyticity and boundedness in
Re s > −γ is control of ∫x^{−γ}d|ν|. That makes P1 three classical problems at once:

* (a) optimal recovery of f(0) from f(1..N) in H^∞ of a half plane;
* (b) extrapolation of a Hausdorff moment sequence to an index outside the data;
* (c) a quantitative Müntz–Szász theorem.

In the Mellin picture the Hardy space H²(Re s > −γ) is L²(x^{2γ−1}dx) on (0,1), evaluation at s is the inner
product with x^s, and the worst |F(0)| with F(1..N) = 0 and ‖F‖ ≤ 1 is dist(1, span{x¹..x^N}). The Cauchy/Gram
determinant evaluates this distance as (2γ)^{−1/2}∏ₙ n/(n+2γ) = (2γ)^{−1/2}β_N. That is Theorem C's β_N exactly
(`p1_muntz_check.py`: agreement to 10⁻³³–10⁻⁶¹).

### P1.1 The classical theorems (a)

* **Blaschke-product extremal property / Schwarz–Pick** [T?]. Take f in the unit ball of H^∞(D) with f(z_j) = 0,
  j ≤ N. Then f = B·g with ‖g‖ ≤ 1, so sup|f(z₀)| = |B(z₀)|, attained by f = B. This is the standard factorization
  (e.g. Garnett, *Bounded Analytic Functions*, Ch. II; Rudin, *Real and Complex Analysis*, Ch. 15). In the half plane
  this is THEOREM_C step 4 and Prop. C′ (k = 0).
* **Uniqueness (N → ∞, ε = 0)** [T?]. A bounded analytic function on a half plane vanishing on {z_n} with
  Σ Re z_n/(1+|z_n|²) = ∞ vanishes identically (Blaschke condition, Garnett Ch. II). For z_n = n + γ the sum
  diverges like log N. The rate of divergence is exactly what makes β_N ~ N^{−2γ}: only logarithmic. The same
  statement is **Müntz–Szász** (Müntz 1914, *H. A. Schwarz Festschrift* 303–312; Szász 1916, *Math. Ann.* 77, 482–496)
  [T; the Rudin chapter containing it is confirmed, but I did not verify its number], and the physics version is
  **Baym & Mermin, "Determination of thermodynamic Green's functions", J. Math. Phys. 2, 232 (1961)** [T, abstract:
  the transform is uniquely determined by its Matsubara values although they have no interior limit point]. That the
  proof uses Carlson's theorem I did not verify.
* **Constructive Matsubara-to-real-time interpolation**: G. Cuniberti, E. De Micheli, G. A. Viano, "Reconstructing
  the thermal Green functions at real times from those at imaginary times", *Commun. Math. Phys.* 216, 59–83 (2001),
  arXiv:cond-mat/0109175 [T, abstract]. They identify G_R as the unique "Carlsonian" interpolant of the Matsubara
  coefficients and build it with Pollaczek polynomials plus regularization. This is the closest published analogue
  of THEOREM_C, but I did not see a sharp stability rate there.
* **Optimal-recovery principle** [T?]. C. A. Micchelli & T. J. Rivlin, "A survey of optimal recovery", in *Optimal
  Estimation in Approximation Theory* (Plenum, 1977), pp. 1–54, following Smolyak (1965) and Bakhvalov (1971). For a
  convex, centrally symmetric class K and a linear functional ℓ, the minimal worst-case error of *any* algorithm
  using data |I f| ≤ ε equals, up to a factor at most 2, sup{|ℓf| : f ∈ K, |If| ≤ ε}. With exact data it equals it
  exactly. Consequence: Prop. C′'s lower bound is a bound on every conceivable reconstruction, including Nevanlinna
  and MaxEnt.
  * Related: Micchelli, Rivlin & Richter-Dyn, "Blaschke products and optimal recovery in H^∞", *Calcolo* (1987)
    [T, bibliographic record only]. Osipenko's papers on optimal recovery in H^∞ with δ-errors [T, abstract: the
    δ-error equals sup|f(z₀)| over the class with |f(z_j)| ≤ δ].
* **Hausdorff moment ill-posedness (b)**: G. Talenti, "Recovering a function from a finite number of moments",
  *Inverse Problems* 3, 501 (1987) [T, bibliographic; it contains a stability theorem whose rate I did not see].
  D. D. Ang, R. Gorenflo, V. K. Le, D. D. Trong, *Moment Theory and Some Inverse Problems in Potential Theory and Heat
  Conduction*, Lecture Notes in Math. 1792 (Springer, 2002) [T, TOC: a Hausdorff-moment chapter with error
  estimates; statements unverified].
* **General ill-posed continuation (Hölder inside, log at the boundary)**:
  * F. John, *Comm. Pure Appl. Math.* 13 (1960) 551–585 [T, review: Hölder dependence under an a priori bound];
  * M. M. Lavrent'ev, *Some Improperly Posed Problems of Mathematical Physics* (Springer, 1967) [T?];
  * K. Miller, *SIAM J. Math. Anal.* 1 (1970) 52–74 [T, via Trefethen's references];
  * S. Vessella, "A continuous dependence result in the analytic continuation problem", *Forum Math.* 11 (1999)
    695–703 [T, bibliographic, via Trefethen; statement not seen];
  * **L. N. Trefethen, "Quantifying the ill-conditioning of analytic continuation", BIT 60 (2020) 901–915,
    arXiv:1908.11097** [T, text read]. Thm 5.1: data on a continuum E with no boundedness give *no* bound off E.
    Thm 5.2: if also ‖g‖_Ω ≤ 1, then |g(z)| ≤ ε^{α(z)} with α ∈ (0,1) (Hölder). Thm 4.1: half-strip with data on
    the end segment, α(x) = (4/π)e^{−πx/2} − (4/π − 1)e^{−πx}. A matching example has exponent (4/π)e^{−πx/2}, so
    digits are lost exponentially.

**Is THEOREM_C itself in the literature?** Its ingredients are all classical: the Blaschke bound, Lagrange
interpolation with Blaschke weights, the Müntz product and optimal recovery. I did not find the combined statement
"half-plane H^∞ + integer samples ⇒ (log M/ε)^{−2γ}, sharp". It is a short exercise for experts. The paper should
cite the pieces and call the result a standard consequence.

### P1.2 What upgrades log to Hölder: the rate ladder

L = log(M/ε). "Data" means |F(n)| ≤ ε for n = 1..N (any N, or all n).

| a-priori class | worst \|F(0)\| | status | physics reading |
|---|---|---|---|
| (S) \|F\| ≤ M on Re s > −γ (THEOREM_C) | ≍ M L^{−2γ} | [P] THEOREM_C + Prop C′; [NE] `out_p1_rates.txt` class A: W·L^{2γ} ≈ 1.0–1.25 (γ = ½), ≈ 2.8–3.3 (γ = 1) | relaxation gap only |
| (S) plus any regularity/analyticity of δη near ω = 0 | same | [P] F_sh is meromorphic with poles only at Re s ≤ −1−2γ and at −a, so it is analytic and smooth at 0 | *does not help* |
| (S) plus a bound on a bigger domain D whose boundary near ∞ is still a vertical line (e.g. H_γ ∪ a disk around 0) | M L^{−2Re φ(0)} | [D] map D → right half plane with φ(s) = s + c + O(1/s) at ∞; the Blaschke product at 0 is ∏\|φ(n)−φ(0)\|/\|φ(n)+φ(0)*\| ~ N^{−2Re φ(0)}. The exponent changes, the log type does not | *does not help* qualitatively |
| sector: \|F\| ≤ M on \|arg(s+γ)\| < α, π/2 < α ≤ π | exp(−Θ(L^{1−p})), p = π/2α | [D] THEOREM_C proof in ζ = (s+γ)^p; β_N = exp(−(2γ^p/(1−p))N^{1−p}), Λ_N = e^{c_p N} with c_p = 1.75 (p = 1), 2.00 (¾), 2.21 (0.6), 2.38 (½) (`out_p1_sector.txt`); lower bound by the Prop C′ construction in ζ [D-sketch] | needs \|F\| ≤ M on part of the second sheet. Einstein QNMs sit at ≈ 30–35° from the negative s axis, so α ≲ 145–150° and p ≳ 0.6: rate exp(−cL^{0.4}) |
| **CM relaxation** (both theory and reference Φ are positive mixtures of e^{−λt}, λ ≥ γ; total masses Φ(0) bounded by the sum rule) | exp(−Θ(√(γL))) | [NE] dual certificates `out_p1_stieltjes_dual.txt`: exact-data error e_N ≈ exp(−3.3√(γN)) (→ 4√(γN), the Blaschke product in √(s+γ)); coefficients ~ e^{1.9N}. LP primal agrees (`out_p1_rates.txt` class B). Faster than every power of L, slower than every power of ε. Not proved asymptotically | **fails for holography**: Einstein's E_E is *not* Stieltjes; the LP residual is 1.7×10⁻³ already with N = 3 points (`out_p1_einstein_positivity.txt`), as its lowest QNM is complex |
| **one-signed deviation**: δΦ ≥ 0 (or ≤ 0) with ∫x^{−γ}dν ≤ M | M^{1/(1+γ)} ε^{γ/(1+γ)} exactly; one data point suffices | [P] Lyapunov (below); [NE] LP class C reproduces it to 4 digits | Hölder, but needs a sign of the deviation |
| data on a continuum E ⊂ H_γ (e.g. all real s ≥ 1) | M^{1−ω}ε^{ω} | [T] two-constants theorem; Trefethen Thm 5.2. For E = [1, ∞): ω = (2/π)arctan(γ/√(1+2γ)) (1/3 at γ = 1) [P] | the Kubo–Mori/form data of Remark 10 (every ν > 0) give this. Euclidean data do not |
| F analytic and bounded at s = ∞ (convergent large-frequency expansion) | Hölder | [D] the data now accumulate at an *interior* point. Lagrange on K clustered nodes gives εC₁^K + Mq^K with q < 1 | false in holography (QNMs accumulate at ∞), so the thermal OPE in ω is only asymptotic |

**Lyapunov proof [P].** s ↦ log∫x^s dν is convex for ν ≥ 0 (Hölder's inequality). Interpolating between s = −γ and
s = 1 at s = 0 gives ∫dν ≤ (∫x^{−γ}dν)^{1/(1+γ)}(∫x dν)^{γ/(1+γ)} ≤ M^{1/(1+γ)}ε^{γ/(1+γ)}. Equality holds for a single
atom. ∎

**Two-constants for E = [1, ∞) [P].** Map Re u > 0, u = s+γ, by u ↦ u² ↦ ψ = u²/((1+γ)²−u²) ↦ √ψ. The slit becomes
{|Im| ≥ 1} on the imaginary axis of the right half plane, and s = 0 becomes r = γ/√(1+2γ). The Poisson integral gives
ω = (2/π)arctan r. ∎

**Bottom line for P1 (d).**
* Theorem C cannot be improved within (S)+(E): this is proved, and the classical theory confirms it.
* Positivity and complete monotonicity improve the rate only if they bound the *total variation of the
  deviation*, or if the deviation is one-signed:
  * CM with a sum-rule mass bound gives exp(−c√L);
  * one-signed deviation gives Hölder.
* No positivity statement available for a generic holographic CFT turns L^{−2γ} into Δ_g^{−p}.
* A power law needs data off the integers (Remark 10's form bounds or real-time norms) or a sign of the deviation.
  This answers the question "is (E) at precision Δ_g^{−4} enough?": no, by optimal-recovery theory, for every method.

### P1.3 Complete monotonicity: what Bernstein/Hausdorff/Widder/Stieltjes actually give

* **Bernstein** [T?]: E is completely monotone on (0, ∞) iff E is the Laplace transform of a positive measure. In
  our setting that is iff Φ ≥ 0 (Widder, *The Laplace Transform*, 1941, Ch. IV). **Hausdorff** [T?] (Math. Z. 9,
  1921): (E(n))_{n≥0} is a Hausdorff moment sequence iff it is a completely monotone sequence, (−Δ)^k E(n) ≥ 0.
* **Φ completely monotone** (Φ = ∫e^{−λt}dμ, μ ≥ 0) **⇔ E is a Stieltjes function** E(s) = ∫dμ(λ)/(λ+s) (Widder
  Ch. VIII; Schilling–Song–Vondraček, *Bernstein Functions*, Ch. 2) [T?]. The two-sided bounds on E(0) from E at
  positive points are the multipoint Padé bounds: M. Barnsley, "The bounding properties of the multipoint Padé
  approximant to a series of Stieltjes on the real line", J. Math. Phys. 14, 299–313 (1973) [T, abstract].
  Convergence conditions on the nodes: G. López Lagomasino, Math. USSR Sb. 35 (1979) 363–376 [T, abstract only].
* **Hölder?** Only for a one-signed deviation (Lyapunov above). For a two-signed deviation between two CM systems the
  rate is exp(−Θ(√L)) [NE]. It is better than THEOREM_C because the TV bound forbids the exponentially large
  residues of F_sh, but it is still not Hölder.
* **New usable statement: an M-free lower bound** [P + NE, `p1_einstein_positivity.py`]. If Φ ≥ 0, then 4πη/s =
  ν([0,1]) ≥ LB_N := min{ν([0,1]) : ν ≥ 0, ∫x^n dν = E(n), n ≤ N}. This is an LP, and also the lower principal
  representation of the truncated Hausdorff problem. LB_N increases to the true value, because a Hausdorff measure is
  determined on (0,1] by its moments n ≥ 1. On the Einstein numbers (`out_einstein_ref.txt`, 4πη/s = 1):

  | N | 1 | 2 | 4 | 6 | 8 |
  |---|---|---|---|---|---|
  | LB_N, exact data | 0.413 | 0.715 (= E₁²/E₂) | 0.897 | 0.958 | 0.982 |
  | LB_N, ε = 10⁻⁴ | | 0.715 | 0.888 | 0.908 | 0.912 |
  | LB_N, ε = 10⁻⁶ | | 0.715 | 0.897 | 0.957 | 0.968 |

  * With the M-bound added (γ = 1, M = 5, N = 8) the bounds are two-sided: 0.982 ≤ 4πη/s ≤ 1.021. M = 1.5 at γ = ½
    is already *inconsistent* with the data.
  * **Caveat:** Φ_E is not ≥ 0. Its late-time tail is 2Re(Re^{−iω₁t}) with complex ω₁, so it oscillates. The first
    eight Matsubara values cannot see this (the LP is feasible). A physical use needs Φ ≥ −φ₋ with ∫φ₋ known, which
    lowers LB_N by at most c∫φ₋.

### P1.4 Nevanlinna/Pick continuation in condensed matter and lattice QCD

* **Fei, Yeh, Gull, "Nevanlinna Analytical Continuation", PRL 126, 056402 (2021), arXiv:2010.04572** [T, abstract].
  * It uses the Pick criterion (an interpolant exists iff the Pick matrix is PSD) and the Schur algorithm to
    parametrize *all* causal interpolants.
  * The matrix-valued version is Fei, Yeh, Zgid, Gull, arXiv:2107.00788.
  * The bosonic version is Nogaki & Shinaoka, J. Phys. Soc. Jpn. 92, 035001 (2023).
* **Rigorous uncertainty bounds**: T. Bergamaschi, W. I. Jay, P. R. Oare, arXiv:2305.16190 (2023) [T, abstract].
  * At each point of the *open* upper half plane the set of values consistent with the data is a disk (Nevanlinna's
    Wertevorrat; R. Nevanlinna, *Ann. Acad. Sci. Fenn.* 1919, as cited by Rogosinski 1936). This gives rigorous
    bounds on *smeared* spectral densities.
  * Follow-ups: Abbott, Jay, Oare, arXiv:2508.01377 (Hankel/moment version), arXiv:2602.11260 (review, with a proof
    that the Pick set is convex); Fields & Christ, arXiv:2510.12136 (uncertain data).
  * Lexing Ying, J. Comput. Phys. 469, 111549 (2022), arXiv:2202.09719 is algorithmic and states the problem is
    ill-conditioned. I found no rate theorem there.
* **What they give at ω → 0.** [P] In the disk picture the Wertevorrat radius at z is at most |B_N(z)|: two
  interpolants differ by B_N·g with |g| ≤ 2. [D] As z = iy → 0, |B_N(iy)| = ∏|(ω_n−y)/(ω_n+y)| → 1, so the upper
  bound on the uncertainty of G_R(iy) does not shrink. After division by y it gives no control of η. Remark 5's
  explicit deformations show that η really is unconstrained.
  * Every published rigorous bound is at *finite smearing* y > 0. For exact data its radius shrinks like
    N^{−y/πT} (the Blaschke product).
  * With an a-priori strip bound it is THEOREM_C again. **Nothing in this literature reaches ω = 0 with a rate**,
    and by P1.1 nothing can within the class.

**P1 verdict (d):**
* *Partially solves*: THEOREM_C is classical and sharp for every algorithm.
* *Proves impossible*: a power law from Matsubara data plus analyticity, positivity or analyticity of δη near 0.
* *New input that would work*:
  * a sign of the deviation (Hölder exponent γ/(1+γ));
  * data on a real-s continuum (Remark 10's form bounds; Hölder exponent (2/π)arctan(γ/√(1+2γ)));
  * a sum-rule TV bound with CM (exp(−c√L) only).

---

## P2. Inclusive moment bounds

Setting: weights cᵢ ≥ 0 at masses mᵢ ≥ m_gap. Put x = (m_gap/m)² ∈ (0,1], so inverse moments Σcᵢmᵢ^{−2a} =
m_gap^{−2a}∫x^a dρ. Known: ∫u_k dρ ≤ B_k (or = B_k). Target: ∫f dρ.

**(a) Theorems.**
1. **Duality** [T?] (K. Isii, "The extrema of probability determined by generalized moments (I)", Ann. Inst.
   Statist. Math. 12 (1960) 119–134, confirmed bibliographically; S. Karlin & W. J. Studden, *Tchebycheff Systems*,
   Interscience 1966, Ch. XII):
   sup{∫f dρ : ρ ≥ 0, ∫u_k dρ ≤ B_k} = inf{Σa_kB_k : a_k ≥ 0, Σa_k u_k ≥ f on the support},
   under a Slater-type condition. With equality data the a_k are free in sign. The sup is attained by a measure with
   at most K atoms (Richter–Rogosinski–Carathéodory).
2. **Chebyshev–Markov–Krein** [T?] (Karlin–Studden Ch. II–IV; M. G. Krein & A. A. Nudel'man, *The Markov Moment
   Problem and Extremal Problems*, AMS Transl. Math. Monographs 50, 1977). Suppose (u₀..u_n) and (u₀..u_n, f) are
   T-systems on [a,b] and the data lie in the interior of the moment space. Then min and max of ∫f dρ are attained
   uniquely at the *lower* and *upper principal representations*. For power moments these are the Gauss and
   Gauss–Radau quadratures of the data (nodes are zeros of the orthogonal/quasi-orthogonal polynomials, one Radau node
   at an endpoint). Equivalently, if f^{(2n)} ≥ 0 and the known moments are 1, x, …, x^{2n−1}, then the Gauss rule is
   the exact minimum (Gauss remainder f^{(2n)}(ξ)/(2n)!·∫π_n²dρ ≥ 0).
3. **Physics precedents**: the arcs and Wilson coefficients of forward amplitudes are Hausdorff moments.
   * Bellazzini, Miró, Rattazzi, Riembau, Riva, "Positive moments for scattering amplitudes", PRD 104, 036006 (2021),
     arXiv:2011.00037 [T, abstract].
   * Arkani-Hamed, Huang, Huang, "The EFT-Hedron", JHEP 05 (2021) 259 [T?].

**(b, c) Explicit optimal bounds** [P, checked by LP, `p2_cmk_check.py` → `out_p2_cmk.txt`].
Known ∫x^{a_k}dρ ≤ B_k with a₁ < a₂ < …; target ∫x^b dρ.
* **b < a₁** (a *lower* inverse power than every known one, e.g. Σc/m² from Σc/m⁴, Σc/m⁶): **sup = ∞**. Put one
  state of weight B₁(m/m_gap)^{2a₁} at m → ∞. Then the target grows like (m/m_gap)^{2(a₁−b)}. No gap helps. **This
  settles the example in the task negatively.**
* **b ≥ max a_k**: sup = min_k B_k, attained by one state at m = m_gap.
* **a₁ ≤ b ≤ a₂, two data**: sup = B₁^θ·min(B₁,B₂)^{1−θ}, θ = (a₂−b)/(a₂−a₁), attained by one atom at
  x* = (B₂/B₁)^{1/(a₂−a₁)}. Example: Σc/m⁶ ≤ √(B₄B₈) when B₈ ≤ B₄ in m_gap units. LP agreement is to 6 digits.
* **Equality data, lower bounds**: inf ∫x^b dρ with ∫x^{a₁}, ∫x^{a₂} fixed and a₁ < b < a₂ equals the given value
  of ∫x^{a₂}dρ, since x^b ≥ x^{a₂} on (0,1]. It is not attained: an atom at x = 1 plus mass escaping to x → 0. More data give the canonical-representation formulas of
  (a2).

**What is needed for a moment that is not a positive combination.** The right condition is not "f is a positive
combination of the u_k" but **pointwise domination f ≤ Σa_k u_k on [m_gap, ∞), including m → ∞**. For inverse
moments this decides by the slowest-decaying known weight. To bound Σc/m² one needs a known quantity decaying no
faster than m⁻², for example:
* Σcᵢ (a less-subtracted sum rule, which needs better Regge behaviour);
* Σcᵢ/mᵢ² itself;
* a UV cutoff on the spectrum, m ≤ m_max, which makes x ≥ x_min > 0 and restores CMK bounds.

(d) *Insta-solves* the bookkeeping (optimal formulas, LP certificates). *Proves impossible* any lower-power moment
from higher-power ones.

---

## P3. Positive functionals, Schur complements, cones smaller than PSD

**(a) Theorems.**
* **Block positivity** [T, operator version stated in arXiv:1904.08680, Thm 5.11, which cites R. Bhatia, *Positive
  Definite Matrices* (Princeton 2007), Prop. 1.3.2 for a related proof; I did not see Bhatia's text]: [[A, X],[X*, D]] ⪰ 0 iff A, D ⪰ 0 and X = A^{1/2}KD^{1/2} for a
  contraction K. Hence |tr W(BC† + CB†)| ≤ 2√(tr(WBB†)·tr(WCC†)) for W ⪰ 0. The support lemma (Remark 7′) is the
  equality case tr(WBB†) = 0. The bound is sharp: the interference ranges over the full disk allowed by the two
  diagonals.
* **Kadison–Schwarz** [T?]: for a unital 2-positive map Φ, Φ(a)*Φ(a) ≤ Φ(a*a). (R. V. Kadison, Ann. Math. 56 (1952)
  494–503; M.-D. Choi, Illinois J. Math. 18 (1974) 565–574.) Again, a linear term is bounded only by a quadratic one.
* **PSD completion** [T?] (Grone, Johnson, Sá, Wolkowicz, Linear Algebra Appl. 58 (1984) 109–124): a partial matrix
  with unspecified diagonal entries places no constraint on the specified off-diagonal entries.
* **Constrained cones / Positivstellensatz** [T?]: Schmüdgen, Math. Ann. 289 (1991) 203–206; Putinar, Indiana Univ.
  Math. J. 42 (1993) 969–984. Measures supported on a semialgebraic set K are characterized by PSD moment *and
  localizing* matrices. Symmetry reduction: Gatermann & Parrilo, J. Pure Appl. Algebra 192 (2004) 95–128.

**(b) Hypotheses in our setting and (c) what they give** [P]. Write the heavy-state contribution as a positive
measure ρ over coupling vectors v = (g₂, g₃)/m. Then (a, b, c) = ∫(g₂²/2, g₂g₃, g₃²/2)dρ/m², and

  **sup{|b| : a ≤ A, ρ ≥ 0} < ∞ ⟺ κ := sup_{supp ρ}|g₃/g₂| < ∞, and then sup|b| = 2κA.**

*Proof.* |g₂g₃| ≤ κg₂² pointwise gives the bound; single atoms with |g₃/g₂| → κ saturate it. ∎

This is the P2 domination principle again. A smaller cone helps exactly when it confines the generating vectors. That
is what crossing does in the known cases:
* **Null constraints** (S. Caron-Huot & V. Van Duong, "Extremal Effective Field Theories", JHEP 05 (2021) 280,
  arXiv:2011.02957) [T, text read]:
  * g₂ > 0 and g₃ ≤ 3g₂/M² follow from positivity alone (their eq. 2.22).
  * The *lower* bound on g₃ needs the crossing null constraint n₄ (eq. 2.23). It gives g̃₃ ≥ −10.61249 in d = 4,
    improving to −10.346 with more constraints.
* **Spectrahedral cones from crossing** (X. Li, H. Xu, C. Yang, C. Zhang, S.-Y. Zhou, "Positivity in Multifield
  Effective Field Theories", PRL 127, 121601 (2021), arXiv:2101.01191) [T, abstract]: the bounds are the dual of a
  spectrahedron built from crossing, with extremal rays given by irreducible UV states.
* **Interference collider bounds**: C. Córdova, J. Maldacena, G. J. Turiaci, "Bounds on OPE coefficients from
  interference effects in the conformal collider", JHEP 11 (2017) 032, arXiv:1710.03199 [T, abstract].
  * ANEC in superposition states bounds Σ_O f(Δ)λ²_{TTO} ≤ N_B, i.e. the ∫ΦW² coupling (our g₂).
  * It is a Cauchy–Schwarz with both diagonals known, the T state and the O state.

**Application to Remark 4 / 7′.**
* The scalar family (g₂, g₃ arbitrary) is a single unitary exchange. It lies in every cone generated by physical
  states, crossing-symmetric or not, so κ = ∞ and no functional bounds b by a.
* The CMT-type bound controls g₂, not g₃. I know of no theorem bounding the XTTT/XTT ratio (or ⟨TTTX⟩ structures) over
  heavy X.

(d) **Proves impossible without new input.** Needed: either the six-point diagonal c (then |b| ≤ 2√(ac), sharp), or
a confinement theorem |g₃| ≤ κ|g₂| for all heavy states, e.g. from a five- or six-point collider/ANEC argument that
does not exist yet.

---

## P4. Spin sums at the back-to-back point (low priority)

**(a) Identities** [T, DLMF text read; checked to 12+ digits in `p4_conical_check.py`].
* In m1_remainder, Ψ_ν(z) = zP_{−½+iν}(X)/(1−z)² with X = (1+z)/(1−z). The back-to-back point z = −1 is X = 0, and
  z ∈ (−1, 0) is X = cos θ ∈ (0, 1).
* DLMF 14.5.1 at μ = 0 (with 14.3.1): P_{−½+iτ}(0) = √π/|Γ(¾ + iτ/2)|² exactly.
* DLMF 14.20.15: P_{−½+iτ}(cos θ) ~ (θ/sin θ)^{1/2}I₀(τθ), hence P_{−½+iτ}(0) ≈ e^{πτ/2}/√(2πτ). At τ = 300 this
  agrees to 7 digits. At τ = Λ = 1000 it is 10^{682}, the size found in m1.
* **DLMF 14.20.9 (Mehler–Dirichlet)**: P_{−½+iτ}(cos θ) = (2/π)∫₀^θ cosh(τφ)/√(2(cos φ − cos θ)) dφ > 0.
* Signature/crossing [T?]: the even-J sum is Σ_{J even}h_J z^J = ½[H(z) + H(−z)], and P_J(−x) = (−1)^J P_J(x).

**(b, c) What follows** [P].
* For the undifferentiated kernel, ∫G(ν)P_{−½+iν}(0)dν = (2/π)∫₀^{π/2}dφ (2cos φ)^{−1/2} Ĝ_c(φ), where
  Ĝ_c(φ) = ∫G(ν)cosh(νφ)dν is the ν-transform of the weight continued to imaginary impact parameter iφ.
  * The kernel is positive, so cancellation can come only from G's sign structure, never from the special function.
  * H_b also contains θ-derivatives in z, i.e. X-derivatives of P. Those kernels grow like e^{πν/2} times powers
    of ν (differentiate 14.20.15), but I do not claim they are positive.
* G is supported in ν ≤ Λ, so Ĝ is entire of exponential type Λ: |Ĝ_c(φ)| ≤ ‖G‖₁e^{Λφ} (Paley–Wiener). The observed
  e^{πΛ/2} is the generic saturation at φ = π/2. This is the m1 statement "heavy image at imaginary impact parameter
  iπ/2", now derived from a classical identity.
* For *physical* partial waves the back-to-back value ½[A(x) + A(−x)] at x = −1 is the crossed (u-channel) amplitude,
  which is finite. For functional actions ĥ_b(J) it has no such meaning unless the functional is s↔u symmetric.

(d) **No identity forces the cancellation.** The only possible mechanism is crossing bookkeeping: if the dispersive
functional is s↔u symmetric, the u-channel image of the same light block must appear elsewhere in the sum rule and
could cancel the X < 1 part. Check this next. It fits the m1 conclusion that the block-wise light side, not the
size estimate, is what fails. Not a special-function question.

---

## Files

| script | what | output |
|---|---|---|
| `p1_muntz_check.py` | Gram/Cauchy-determinant distance = (2γ)^{−1/2}β_N (Theorem C = Müntz) | `out_p1_muntz.txt` |
| `p1_rates_lp.py` | LP worst cases: class A (half-plane/TV, L^{−2γ}), B (CM/Stieltjes), C (one-signed, Lyapunov) | `out_p1_rates.txt` |
| `p1_stieltjes_dual.py` | 120-digit dual certificates for class B (rigorous up to grid sup), ε down to 10⁻⁴⁰ | `out_p1_stieltjes_dual.txt` |
| `p1_sector_rate.py` | sector class: β_N, Λ_N growth, optimized bound | `out_p1_sector.txt` |
| `p1_einstein_positivity.py` | Einstein data: Φ ≥ 0 lower bounds on 4πη/s, two-sided with M, Stieltjes infeasibility | `out_p1_einstein_positivity.txt` |
| `p2_cmk_check.py` | optimal moment bounds, closed form vs LP | `out_p2_cmk.txt` |
| `p4_conical_check.py` | P_{−½+iτ}(0): Γ formula, asymptotics, Mehler–Dirichlet | `out_p4_conical.txt` |

LPs use scipy/HiGHS on grids (fine log grids; tolerances 10⁻¹⁰). They are numerical evidence, not certificates,
except where marked.

## Sources consulted (verification level)

* Text read: Trefethen BIT 2020 (Thms 4.1, 5.1, 5.2 and reference list); Caron-Huot–Van Duong arXiv:2011.02957
  (eqs. 2.22–2.23, g̃₃ bounds); DLMF §§14.12, 14.20.
* Abstracts or records seen:
  * Cuniberti–De Micheli–Viano CMP 2001; Baym–Mermin JMP 1961; Fei–Yeh–Gull PRL 2021; Bergamaschi–Jay–Oare 2023;
    Abbott–Jay–Oare 2025/2026;
  * Bellazzini et al. PRD 2021; Li et al. PRL 2021; Córdova–Maldacena–Turiaci JHEP 2017;
  * Demanet–Townsend FoCM 2019 (Hölder extrapolation from equispaced samples on an interval, exponent
    −log r(x)/log ρ, under oversampling M* ≤ ½√N);
  * Barnsley JMP 1973; López Lagomasino Sb. Math 1979; Isii 1960; Talenti 1987; Ang et al. LNM 1792; John CPAM 1960;
    Lexing Ying JCP 2022.
* Standard results stated without checking the source text (marked [T?]): Blaschke/Garnett, Müntz–Szász (the Rudin
  chapter exists; its number is unverified), Bernstein/Hausdorff/Widder, Karlin–Studden and Krein–Nudel'man
  principal representations, Kadison–Schwarz, Grone et al., Schmüdgen/Putinar, Gatermann–Parrilo, Micchelli–Rivlin
  1977, EFT-hedron.
