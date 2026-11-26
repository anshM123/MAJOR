# TAIL_CERT: the large-impact-parameter tail L ≥ L_tail of full AdS positivity at Δ_g = 1000

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*Code: `ads4/code/tailcert/` (outputs in `cache/`, git-ignored). Compute: ≤ 2 processes, `OMP_NUM_THREADS=1`, `ulimit -v 4000000`.*

Grades: **[certified]** Arb balls / exact rationals with every remainder bounded; **[P]** proved on paper here (proof given or
sketched with every inequality stated); **[DV]** sampled numerics; **[author]** claim not reproduced.

---

## 0. Verdict

**The tail closes for L ≥ 20, i.e. β_c ≥ 21.4, every twist a ≥ 998, every spin, both channels, including all connection terms
j ≥ 1, the first-moment dressing and the k = 3 piece** [certified arithmetic on top of the [P] lemmas of §§1–5; the CONTOUR_BOUNDS
Arb tables are used as certified inputs]:

  **H_ang,b ≥ (1 − ρ)·𝓜(L),  ρ ≤ 0.066 for all L ≥ 20**,  𝓜(L) ≥ (2/Λ)√w ∫₀^{3/4}W(y)e^{−Ly}dy > 0  (w = 1 − r = e^{−L}).

| piece (ratio to the certified main-term lower bound) | L = 20 | L = 24 | L = 30 | grade |
|---|---|---|---|---|
| ℰ₂-type: detour + leg (j = 0, incl. the "1" parts and the dressing) | 8.2·10⁻⁵ | 1.2·10⁻⁶ | 1.5·10⁻⁸ | certified (given §2 lemmas) |
| segment D (local multiplier on W ≥ 0) | 3.3·10⁻⁶ | same | same | certified |
| connection terms 1 ≤ j < 200 (all pieces) | 5.6·10⁻¹⁰ | 5.6·10⁻¹⁰ | 5.5·10⁻¹⁰ | certified (given §3 recurrence) |
| remainder j ≥ 200 (localisation) | < 10⁻⁴³⁵ | smaller | smaller | [P] + certified arithmetic |
| k = 3 piece (absolute, real axis, its own j-sum) | **0.065** | 1.4·10⁻³ | 4.2·10⁻⁶ | certified (given §2) |
| k = 3 remainder | < 10⁻⁶⁵² | | | [P] + certified arithmetic |

* **Threshold.** L_tail = 20 (β_c ≥ 21.4) is far below the cell certificate's reach (β ≤ 64 ⇔ L ≲ 62.6): wide overlap.
  The binding piece is k = 3 (∝ e^{−L}). **L ≥ 18 also passes: ρ ≤ 0.68** (k = 2: 4.8·10⁻³, k = 3: 0.673, remainders < 10⁻²⁶⁴), so
  H_ang,b ≥ 0.32·𝓜 for every L ≥ 18 (β_c ≥ 19.4) (§6).
* **What changed against §§38/40/43.** The three author inputs that CONTOUR_BOUNDS had to *assume* (\|D\| ≤ 100/a², 300/a² on the
  detour, 5·10¹² on the leg) are now **proved** (§2), and much sharper: \|D\| ≤ 3.3·10⁻⁶ on the segment/detour; on the leg
  \|w^βS(β)\| ≤ 1.8–20 for t ≤ 45, where all of ∫\|fB\| sits. This alone moves the ℰ₂ threshold from 30.16 to below 20.
* **Item 3 (realised-weight residual).** Not an issue in this route (§R.3): the main term uses the realised f exactly. The coordinator's
  numerical request and the author's Lemma R are answered in §R.
* **Independent of the author's §82 Lemma R′.** I reached the same structure by a different exact representation (Theorem P below:
  a q-integral with an exact Beta(N+1, n+s) average, instead of §82's kernel K(q)). §82's beta representation and the sampled
  K(q)-bound were spot-checked [DV] (§2.4); §82's own constants were not needed.
* **Caveats** (§7): the per-cell f-integrals are CONTOUR's certified tables (not recomputed); a handful of e^{−1000}-size terms in
  Theorem P and Lemma Loc are bounded on paper, not in balls; the whole document is unreviewed.

---

## 1. Exact framework [exact / DV] (`tc_exact.py`, `cache/tc_exact.txt`)

State (a, n, b): a = Δ − J − 2 ≥ 998, n = J − b, s = 2b+1, N = n + a, T = N(N+s) + b(b+1), C = J(J+1), r = C/T, **w = 1 − r = e^{−L}**,
T − C = a(N+n+s).
* (E1) U_k = Π_{j<k}(C−u_j)/(T−u_j) = (−n)_k(n+s)_k/((−N)_k(N+s)_k), so **S(β) := ℒ((1−z)^{−β}) = ₃F₂(β, −n, n+s; −N, N+s; 1)**.
* (E2) 𝒱_b = ℒ(· + 𝒟_b·/T) on powers: 𝒱_b((1−z)^{−β}) = S(β) + T⁻¹[A₁S(β+1) + A₀S(β)], A₁ = β(β+1)/2 + (4−2b)β,
  A₀ = −β(β+1)/2 − (7/2−2b)β + 9/2 − b(b−1) (= §36's form).
* (E3) **Exact three-term recurrence** (from the annihilator Σ_kU_k[(T−u_{k−1})F_{k−1} − (C−u_k)F_k] = 0 of the ladder):
  [N(N+s) − β(β−s)]S(β−1) = [a(N+n+s) + β(s−2β−1)]S(β) + β(β+1)S(β+1).
* (E4) Branch formula: (1−z)^{−3}P_{−½+iν}((1+z)/(1−z)) = 2Re[B(ν)Σ_j c_j(ν)(1−z)^{j−5/2−iν}], B = Γ(2iν)/Γ(½+iν)²,
  c_j = ((½−iν)_j)²/((1−2iν)_j j!) (|1−z| < 1).
* (E5) Generating function: Σ_kP_kz^k = (2/Λ)(1−z)^{−3}∫₀^∞E_θ(μ)P_{−½+iμ}(X)dμ with the shifted density (§32) — checked [DV] on a toy
  weight with the same analytic structure (`scratchpad`, 15 digits, k ≤ 4).
* E1–E3 hold in exact rationals on 30 random states and rational β; E4 to 15 digits at three (ν, z).

Hence (FULLPOS3 §4.1(a)) **H_ang,b = w³[𝒱_b(G_A) + σ_bT⁻¹ℒ(G_D)]**, G_A = (2/Λ)∫₀^∞f(ν)F_{3,ν}dν (f = CONTOUR's realised density),
G_D = Λ²∫₀¹ψ₃(p)F_{4,Λp}dp. For every J₀,
  F_{3,ν} = Σ_± B(±ν)Σ_{j<J₀}c_j(±ν)(1−z)^{j−5/2∓iν} + R_{J₀}(z;ν),
where each truncated term is analytic in |z| < 1, so 𝒱_b may be applied termwise (finite sums); R_{J₀} is handled in §4.

**Contour (per j).** I_j := 2Re∫₀^∞fBc_j𝒱_b(E_{α−j})dν (α = 5/2+iν, E_γ = (1−z)^{−γ}) is deformed separately for every j to: quarter
circle at 0 (gives f(0)X_j(0)/2), segment ν = iy, 0 ≤ y ≤ ¾ (integrand W(y)X_j(iy), W = −2B(iy)Im f(iy) ≥ 0 [CONTOUR, certified]),
short detour t + 3i/4, vertical ¼ + iy (¾ ≤ y ≤ 9/4), leg t + 9i/4. Here X_j := w^{α−j}𝒱_b(E_{α−j}), and
\|w³fBc_j𝒱_b(E_{α−j})\| = √w e^{−L Im ν}\|fB\|\|c_j\|w^j\|X_j\|. (Analyticity: c_j has poles only in Im ν < 0; B and f are singular only on
the imaginary axis, avoided by the detour; Gaussian decay of f at Re ν → ∞ against at most e^{π\|ν\|} growth of X_j.)

**Main term.** 𝓜(L) := (2/Λ)√w[f(0)/2 + ∫₀^{3/4}W e^{−Ly}dy], f(0) = r(0) = 9.69·10⁻¹⁴ > 0 (EPSW). The lower bound used is
(2/Λ)√w∫ minW·e^{−Ly}(1 − δ_seg) over CONTOUR's 7500 certified W-cells.

---

## 2. The multiplier of a pure power: proved bounds (`tc_power.py`)

Ŝ(β) := w^βS(β). Needed: (i) \|Ŝ − 1\| small for β near the real axis (segment, detour); (ii) \|Ŝ\| moderate on the leg β = ¼ + it
(and β+1 = 5/4 + it), t ≤ 1146; (iii) σ = 7/2, 9/2 on the real axis for k = 3.

### 2.1 Lemma NB (near-real) [P + certified constants]
Ŝ(β) = (w/d)^β Σ_kπ_k(β)E_k, d = 1 − c, c = U₁, π_k = (β)_kc^kd^β/k! (Σπ_k = 1), E_k = U_k/c^k ∈ [0,1] (FULLPOS3 §4.1(b)).
* For real β > 0: π_k ≥ 0, so Ŝ − (w/d)^β ∈ [−Σπ_k(1−E_k), 0]. With 1 − E_k ≤ −log E_k ≤ (1+ε)(d/c)hS₁(k), S₁(k) ≤ k³/3 + sk²/2 for
  k ≤ n/2 (hu_j/c ≤ ¼), the exact NB factorial moments E[(K)_m] = (β)_m(c/d)^m, h/d² = N(N+s)/(a²(N+n+s)²) ≤ 0.2501/a², and Markov at
  20th order for k > n/2: **\|Ŝ(β) − 1\| ≤ 1.0001·0.2501[(β)₃/3 + O(d)]/a² + 10⁻¹²**, i.e. ≈ (β)₃/(12a²) (the [DV] fixed-twist value is
  exactly −(β)₃/(12a²), `scratchpad/nlim.py`). At β = 5/2: **3.30·10⁻⁶** (author §35/§38: 94/a² ≈ 9.4·10⁻⁵, 100/a²; §82: 12/a²).
* Complex β with \|τ\| ≤ ¼: \|(β)_k\| ≤ (σ)_k exp(τ²(σ⁻² + σ⁻¹)/2) =: C_τ(σ)_k, so \|Σπ_k(β)(E_k−1)\| ≤ C_τ·(the real bound at σ).
  (w/d)^β = (1 − b(b+1)/T)^β: \|·−1\| ≤ 21\|β\|/T. Dressing (E2): T⁻¹/w ≤ w/(3.9a²), T⁻¹ ≤ w²/(3.9a²).

### 2.2 Theorem P (oscillatory) [P + certified evaluation]
**Exact representation** (Re β > 0), verified [DV] in `scratchpad/chk_rep.py` (10 digits) and its ingredients exactly:
  S(β) = κ(β)∫₀^∞q^{β−1}g(q)dq,  κ = Γ(N+s)/(Γ(β)Γ(N+s−β)),  g(q) = Σ_{K≤n}Bin(M,p)(K)·C(N−K,a)/C(N,a),
  M = N+n+s, p = q/(1+q). (Beta integral with T ~ Beta(β, N+s−β): E T^k = (β)_k/(N+s)_k, then Pfaff for the terminating
  ₂F₁(−n, n+s; −N; t).) On the real axis g = A − T with
  **A(q) = E_V[((1 − λ_Vq)/(1+q))^a], V ~ Beta(N+1, n+s), λ_V = 1/V − 1** (exact: E V^{−j} = (M)↓_j/(N)↓_j; checked to 12 digits),
  and T(q) = Σ_{K>N}(…) is ≤ e^{−M·KL(N/M‖p)}·(1.01)^a, astronomically small for p ≤ 0.479 (Chernoff).
**Proof outline.**
1. Split ∫₀^∞ = ∫₀^{q₁}(A − T) + ∫_{q₁}^∞g. Fubini over V on the first piece; for each V deform [0, q₁] to the path
   q = iy (0 ≤ y ≤ ρ), arc \|q\| = ρ down to q₁ = ρ (integrand analytic in the quarter disc, integrable branch at 0).
2. V-concentration (Beta–binomial duality P(V ≤ v) = P(Bin(M,v) ≥ N+1), Chernoff + Pinsker): outside \|V − N/M\| ≤ η,
   η² = (0.75a + πτ/2 + 60)/M, the contribution is ≤ 10⁻³⁰ (uses E V^{−2a} ≤ 2.0001^{2a}); inside, **\|λ − 1\| ≤ w + 6.6√w ≤ 7√w**
   (a/N ≤ 0.5003w, M ≥ a/(1.0001w)).
3. Prefactor: \|w^βκ\| ≤ (w(N+s))^σ e^{ε_Γ}/\|Γ(β)\|, w(N+s) ≤ 2a, \|Γ(x)/Γ(x−β)\| ≤ x^σ exp(2(\|β\|²+\|β\|)/x) (Binet; x ≥ 998e^{L}).
   In x = 2aq units the path integral is (2a)^{−σ}∫\|x\|^{σ−1}e^{−τ arg x}\|h̃\|\|dx\|, h̃ = ((1−λx/2a)/(1+x/2a))^a.
4. Imaginary leg: \|h̃(iy)\| ≤ exp((λ²−1)₊y²/(8a)) → G·Y₁^σ/σ·e^{c₁}, G := sup_cell\|1/Γ(β)\|e^{−πτ/2} (Arb boxes).
5. Arc, exponent a·f(φ), f = u(π/2 − φ) + log((\|1−ρe^{iφ}\| + δ)/\|1+ρe^{iφ}\|), u = τ₂/a, δ ≥ \|λ−1\|ρ:
   * **regime A** (u ≤ 0.52, ρ = u, Y₁ = 2τ₂): log((1−x)/(1+x))/2 ≤ −x gives f ≤ u(π/2−φ) − 2ρcosφ/(1+ρ²) + δ/(1−ρ), and
     (π/2−φ) ≤ (π/2)cos φ·… : 2ρ/(1+ρ²) ≥ (π/2)u ⟺ u ≤ 0.5227. So the arc costs ≤ (π/2)e^{c₂}, c₂ = δτ₂/0.48, **uniformly in a**.
   * **regime B** (u > 0.52, so 998 ≤ a ≤ τ₂/0.52 ≤ 2204; ρ = 0.92): max_φ f certified on a two-level φ-grid (Arb), tabulated in u
     (f is nondecreasing in u), cost e^{a·max(f,0)}; at u ≤ 0.8 the max is ≤ 2.4·10⁻⁴, at u = 1: 5.5·10⁻³, u = 1.148: 0.060.
6. Real tail ∫_{q₁}^∞: on [q₁, 0.92], A ≤ e^{−1.99aq}(δ-perturbed), giving \|1/Γ(β)\|Γ(σ, 0.995Y₁)/0.995^σ; beyond 0.92,
   g ≤ (1 − 0.99(M/N)p₁)^a + Chernoff ≤ 0.0516^a + e^{−10⁶}, and g ≤ 2^{−a−s}((1+q)/4)^{−(N+s)} for q ≥ 3: ≤ 10⁻³⁰ for τ ≤ 1146.
Result (all tail states, L ≥ 20): **\|Ŝ(¼+it)\| ≤ 2.3 (t = ¼), 3.1 (1), 8.9 (10), 20 (45), 55 (300), 440 (t ≈ 600), 600 (t ≈ 800),
3.2·10⁴ (t ≈ 1000, regime B at a ≈ 998)**, against the [DV] fixed-twist truth 1.0–5.6 for t ≤ 998 at a = 998 (`scratchpad/nleg.py`). Valid for t ≤ 1146.
For t ≥ 1146: \|Ŝ(σ+it)\| ≤ Γ(σ)\|1/Γ(σ+it)\|·Ŝ(σ) ≤ Γ(σ)Γ(1−σ)cosh(πt)/π (reflection; U_k ≥ 0, \|Γ(x+iy)\| ≤ Γ(x)).
Sanity [DV] (`tc_power_dv.py`): exact S at n = 3·10⁵, a = 998 (L = 5.0): \|Ŝ(¼+10i)\| = 1.00003, \|Ŝ(¼+100i)\| = 1.0031, below the bound.

### 2.3 Connection terms by the recurrence [P]
From E3, with \|T/D − 1\|, \|β(s−2β−1)\|/(w\|D\|) ≤ 10⁻⁶ and w²\|D\| ≥ 3.9a² (w(N+s) ≥ 1.999a):
  **\|Ŝ(β−1)\| ≤ (1+10⁻⁶)\|Ŝ(β)\| + \|β(β+1)\|/(3.9a²)·\|Ŝ(β+1)\|**  (\|β\| ≤ 1500).
So every X_j (j ≥ 1) is bounded from \|Ŝ(α)\|, \|Ŝ(α+1)\| alone; the growth is ≤ (½ + ½√(1 + (\|α\|+j)²/(0.975a²)))^j ≤ 1.5^j for
j ≤ 200, against w^j ≤ e^{−20j}. This replaces §40's integration-by-parts/derivative-kernel argument and §82's Mellin–Barnes sum.

### 2.4 The author's §82 inputs, spot checks [DV] (`tc_author_dv.py`)
* §82 beta representation S_b(α) = Γ(N+1+α)/(Γ(N+1)Γ(α))∫q^{α−1}(1−q)^aK(q)dq: agrees with the direct sum to 12 digits (real α) and
  7 digits (α = ¼ + 3i, quadrature-limited) at (n,a,s) = (12,4,1), (15,6,9).
* \|K(q)\| ≤ 2\|1+q\|^{−a}: sampled max of \|K(q)\|\|1+q\|^a = 1.0006 at n = 1024a², a = 4 (15 points of the quarter disc).
* Not reproduced in Arb: §82's own constants (12/a², 4·10⁴e^{G}, MB remainder e^{−2800}); not needed, since §§2.1–2.3 and §4 give
  independent, sharper certified bounds for the same quantities.

---

## 3. Assembly of the k = 2 contour (`tc_assemble.py`, `cache/tc_assemble_L20.0.*`) [certified given §§1–2]
* Inputs: CONTOUR's certified per-cell upper bounds of ∫\|fB\| (cb_short 25 cells, cb_vert 24 000, cb_horiz_2.25 14 568 cells to
  T_END = 1054; row values include the Cauchy/far errors) and W-cells (cb_w). Leg tail [1054, 1146) recomputed with CONTOUR's tail
  formula times Theorem P per unit cell; [1146, 6000) with the reflection bound; beyond 6000 the Gaussian factor is e^{−1.6·10⁷}.
* Every cell carries \|X₀\| + Σ_{1≤j<200}\|c_j\|w^j\|X_j\| (c_j enclosed on the cell in acb balls; X_j by §2.3, dressing included).
* ρ(L) ≤ Err(L_k)/𝓜_lb(L_{k+1}) on cells of width 0.05 over [20, 80] (all errors decrease in L, 𝓜_lb decreases); beyond L = 80 every
  off-segment error carries e^{−3L/4}, the segment j ≥ 1 part is a fixed 5.3·10⁻¹⁰ fraction, so the ratio keeps decreasing [P].
* **Result: ρ_{k=2}(L) ≤ 8.23·10⁻⁵ on [20, 80]** (L = 20: leg 5.6·10⁻⁵, short 1.6·10⁻⁵, vertical 1.0·10⁻⁵, j ≥ 1 on the segment
  5.6·10⁻¹⁰, residue j ≥ 1 10⁻²⁰); 1.5·10⁻⁸ at L = 30.
* Leg: Σ\|fB\|(\|X₀\| + j-sum) = 2.56·10¹² against ∫\|fB\| = 5.30·10¹¹, i.e. an effective multiplier 4.8 instead of §38's 5·10¹².

## 4. The remainder R_{J₀}, J₀ = 200 (`tc_remainder.py`) [P + certified arithmetic]
**Lemma Loc.** If R is analytic in \|ζ\| < 1 and \|R(ξ)\| ≤ K\|1−ξ\|^m on \|ξ\| ≤ ρ₂ = (1+ρ)/2, ρ = n/(n+a−m), 2 ≤ m ≤ a−3, then
  \|𝒱_b(R)\| ≤ 2.53K(1.0006w)^{m−1}Q_m(1 + 40·1.5^m),  Q_m = ((999+m)/(997−m))^{m+1}.
*Proof.* ℒ(F) = E_B[[z^n](1−z)^{−a−1}F(Bz)]/C(N,a), B ~ Beta(n+s, a) (§34, checked); Cauchy on \|z\| = ρ; \|1−Bz\| ≤ \|1−z\| + (1−B);
(1/2π)∫\|1−ρe^{iθ}\|^{−γ}dθ ≤ (1−ρ)^{1−γ}√πΓ((γ−1)/2)/(4√ρΓ(γ/2)) (sin x ≥ 2x/π); with ρ = n/(n+γ′), γ′ = a−m, the entropy bound
for C(n+γ′, γ′) and Gautschi's inequality give I₁ ≤ 2.53(a/N)^m exactly (the e^{a²/n} factors cancel), I₂ ≤ 2.53((a+m)/(a−m))^m(a/N)^m
(E(1−B)^m = (a)_m/(N+s)_m); a/N ≤ 0.5003w. Dressing by Cauchy on circles of radius ½min(\|1−ζ\|, 1−ρ) (T⁻¹/(1−ρ)² ≤ (a−m)^{−2}). ∎
**K for R_{J₀}** (m = 196): near \|1−ξ\| ≤ ½ the branch-series tail (\|arg(1−ξ)\| < π/2 gives e^{πν/2}); far: FULLPOS3 Lemma B
(\|P_{−½+iν}(cosh ξ)\| ≤ (π/2)e^{ν\|Im ξ\|}) plus the truncated sum; \|ν\| ≤ 0.4 by the maximum principle in ν (R is analytic there:
the ν = 0 poles of B(±ν) cancel) with the Pfaff-series bound \|F_{p,ν}(ξ)\| ≤ \|1−ξ\|^{−p}2^{0.9}e^{0.2π}/(1−\|ξ\|).
\|B(ν)\|² = coth(πν)/(4πν) exactly for real ν. Weight: \|f_θ\| ≤ sup\|ψ\| + sup\|r\| + sup\|f_θ − f\| ≤ 21.34 on [0, 1300] (EPSW; sup\|f_θ − f\| ≤ 10⁻²¹ by `tc_cutoff.py`,
within the 5·10⁻¹⁰ slack of the code), CONTOUR's Gaussian tail formula (θ included exactly) beyond. **Result: log₁₀(remainder/𝓜) ≤ −435 for L ≥ 20.**

## 5. The k = 3 piece (`tc_k3.py`) [certified given §2]
H^{(3)} = w³σ_bT⁻¹ℒ(G_D) bounded in absolute value on the real axis (no contour needed): T⁻¹/w = 1/(T−C) ≤ 0.2515w/a², and
\|ℒ(F_{4,ν})\| ≤ 2\|B\|w^{−7/2}[\|Ŝ(7/2+iν)\| + Σ_j\|c_j\|w^j\|Ŝ(7/2+iν−j)\|] + \|ℒ(R′)\|. ψ₃ is enclosed per p-cell by an exact
fmpq Taylor shift (a naive ball evaluation of the degree-14 polynomial overestimates by 10³, an earlier run's bug, fixed).
∫₀¹\|ψ₃\|\|B(Λp)\|X₃ dp ≤ 329 (∫\|ψ₃B\| ≈ 0.2 [DV]; X₃ ≤ 330 for ν ≤ 510, ≤ 2·10⁴ above, regime B).
**Result: \|H^{(3)}\|/𝓜 ≤ 0.065 at L = 20, decreasing like e^{−L}L**; k = 3 remainder ≤ 10⁻⁶⁵².

## 6. Threshold (`cache/tc_assemble_L18.0.*`, `tc_k3_L18.0_J200.txt`, `tc_remainder_k2_L18.0_J200.txt`)
Theorem P's constants depend on L₀ only through \|λ − 1\| ≤ 7e^{−L₀/2} (the arc table `ehat_d_table_L<L₀>.json` is rebuilt per L₀).

| L₀ | ρ_{k=2} (max over [L₀, 80]) | ρ_{k=3} | k=2 remainder | total ρ | conclusion |
|---|---|---|---|---|---|
| 20 | 8.23·10⁻⁵ | 0.0650 | < 10⁻⁴³⁵ | **≤ 0.066** | H_ang,b ≥ 0.93𝓜 > 0 for L ≥ 20 |
| 18 | 4.85·10⁻³ | 0.673 | < 10⁻²⁶⁴ | **≤ 0.68** | H_ang,b ≥ 0.32𝓜 > 0 for L ≥ 18 |

Going lower is limited by k = 3 (×e per unit of L) and the leg (×e^{2.25}). Nothing more is needed for the overlap: L ≥ 18
(β_c ≥ 19.4; for b = 4, d = 1 − c ≥ w, so L_d = 2 log cosh(β_c/2) ≤ L_w) is deep inside the cell certificate's β ≤ 64.
**Note:** the cell certificate covers a ≥ a*(β) only; the tail covers every a ≥ 998.

## 7. Caveats, gaps, and what was not done
1. **Inputs taken from other certificates:** CONTOUR's per-cell ∫\|fB\| tables and W-cells (Arb, CONTOUR_BOUNDS §§2–4); EPSW's r(0),
   ε_W, B, sup\|r\|; LRINF's F_cert cells (Lemma R only). Not recomputed.
2. **[P] steps not ball-evaluated:** the 10⁻³⁰-size pieces of Theorem P (T-part, V-concentration complement, q ≥ 0.92 real tail),
   the e^{−11000} cutoff piece on [0, 1], the Binet bound for Γ(x)/Γ(x−β), the analyticity/decay justifications of the contour shifts,
   and Lemma Loc's inequalities (all stated with explicit constants above). Elementary, but unreviewed.
3. **Theorem P validity:** τ ≤ 1146 (the reflection bound is used beyond). Regime-B arc maxima are certified on a u-grid using
   monotonicity in u.
4. **E5** (generating function with the shifted realised density) is checked on a toy weight [DV]; it is the identity FULLPOS3 §4.1
   and ANGULAR_ROUTE §32 use, not re-derived here.
5. The DV fixed-twist profiles (`scratchpad/nlim.py`, `nleg.py`, not in the repo) are sanity checks only.
6. Not done: an Arb reproduction of the author's own §82 constants (12/a², 4·10⁴e^{G}, MB e^{−2800}). Unnecessary for the tail, since
   §§2–4 give independent certified bounds for the same objects.

## 8. Reproduce (`ads4/code/tailcert/`; `ulimit -v 4000000; export OMP_NUM_THREADS=1`)
```
python3 tc_exact.py                 # E1-E4                                  (seconds)
python3 tc_power.py                 # Theorem P / Lemma NB spot values       (seconds; builds cache/ehat_d_table_L20.json)
python3 tc_assemble.py 20 200       # k = 2 contour, j < 200                 (~3 min)
python3 tc_remainder.py 20 200      # R_200 by Lemma Loc                     (~5 min)
python3 tc_k3.py 20 200             # k = 3                                  (~3 min)
TC_L0TAB=18 python3 tc_assemble.py 18 200; TC_L0TAB=18 python3 tc_k3.py 18 200; python3 tc_remainder.py 18 200
python3 tc_cutoff.py; python3 tc_lemmaR.py; python3 tc_lownu_dv.py      # cutoff audit, Lemma R, item-3 numbers
python3 tc_power_dv.py 300000; python3 tc_author_dv.py                  # [DV] sanity checks
```

---

## R. Item 3 and the author's Lemma R (coordinator request, done first)

### R.1 Signed low-frequency residual, measured [DV] (`tc_lownu_dv.py`, `cache/tc_lownu_dv.txt`)
G_r^low(β) := (2/Λ)∫₀^{ν₁} r(ν)P_{−½+iν}(cosh β)dν, r = E − ψ₂(ν/Λ)/k₀, E = 2Re[χR₋](ν−i), repaired weight (θ = 1 to 10⁻²⁰⁰ here).
This is the low-ν part of (1−c)³G^r(c) (FULLPOS3 §4.1(c)), the term whose *absolute* bound grows like β² in FULLPOS5 §3.
P from the two-branch formula (checked against mpmath `legenp` to 14 digits at β = 30).

| β | F_cert | signed G_r^low/F_cert, ν₁ = 0.31 | ν₁ = 1.0 | absolute bound (2/Λ)∫₀^{ν₁}\|r\|·P_{−½}/F_cert, ν₁ = 0.31 |
|---|---|---|---|---|
| 30 | 1.176·10⁻¹² (LRINF cell [30, 30.5]) | **+0.0021** | +0.0021 | 0.011 |
| 64 | 3.61·10⁻²⁰ (fcert_tail) | **+0.0031** | +0.0031 | 0.032 |
| 100 | 3.52·10⁻²⁸ | **+0.0036** | +0.0037 | 0.077 |
| 200 | 3.40·10⁻⁵⁰ | **+0.0040** | +0.0041 | 0.306 |

The signed piece is **positive** and saturates at a few ×10⁻³ of the floor, while the absolute bound grows like β² (0.31 at β = 200),
exactly as FULLPOS5 §3 found. So the β² growth is an artefact of taking |r|.

### R.2 Lemma R (ANGULAR_ROUTE §81) reproduced in Arb [certified, given EPSW's certified inputs]
(`tc_lemmaR.py`, `cache/tc_lemmaR.txt`)
* **Inputs.** From `epswcert/cache/ew_cert_1000_re_rep.json` (EPSW_CERT, certified): r(0) = 9.6925·10⁻¹⁴, ε_W ≥ ∫\|r\| = 2.68136·10⁻⁶,
  B_EPSW ≥ \|r′(0⁺)\| + ∫\|r″\| = 1.89971·10⁻³. For r̃ = r_θ − r(0)sech(πν):
  * E_in = ε_W + r(0)/2 + δE_θ ≤ **2.68236·10⁻⁶** (author: < 2.683·10⁻⁶ ✓);
  * B_in = B_EPSW + π²r(0) + δB_θ ≤ **1.89971·10⁻³** (author: < 1.901·10⁻³ ✓).
  * Here ∫sech(πν) = ½, ∫\|(sech πν)″\| ≤ π² (sech″ = sech − 2sech³). **Cutoff correction, corrected form (audit §85)** [certified,
    `tc_cutoff.py`, `cache/tc_cutoff.txt`]: f_θ = 2Re[θ(ν−i)g], g = χR₋(ν−i), so \|f_θ − f\| ≤ 2\|θ(ν−i) − 1\|\|g\| (not via \|Re g\|).
    Enclosed on 0.05-cells over [1, Λ+80] (θ at the complex argument, RBF disc bounds, Cauchy for the second derivative) plus CONTOUR's
    Gaussian tail beyond: **δE_θ = ∫\|f_θ − f\| ≤ 1.7·10⁻²³, δB_θ = ∫\|(f_θ − f)″\| ≤ 4.4·10⁻²², ∫e^{0.02ν}\|f_θ − f\| ≤ 1.3·10⁻¹⁴**
    (author §85: < 10⁻⁷); on [0, 1] \|θ(ν−i) − 1\| ≤ e^{−12000} and the 1/ν parts of (g)″ cancel in the real part ([P]).
    *My first version of this item argued with "1 − θ ≤ erfc(7) times \|E\|", i.e. the \|Re g\| form the audit flags; it is replaced by
    the computation above. The 10⁻⁹ allowances used in `tc_lemmaR.py` remain valid upper bounds, so no number changes.*
* **Mehler step, re-derived [P].** W(t) = ∫r̃ cos νt, \|W\| ≤ min(E_in, B_in/t²) (two integrations by parts), ∫₀^∞W = (π/2)r̃(0) = 0.
  With k(t) = (cosh β − cosh t)^{−½} on [0, β) and K₀ = k(0):
  I_{r̃} = (√2/π)∫₀^∞W(t)[k(t)1_{t<β} − K₀]dt, and
  **\|I_{r̃}\| ≤ (√2/π)K₀[B/β + 8B log2/β² + E·ε₁(β)]**, ε₁ = (sinh(β/2) − β/2)/(2(cosh β − cosh(β/2))) ≈ e^{−β/2}/2.
  * [0, β/2]: √(1+x) ≤ 1 + x/2 gives k − K₀ ≤ K₀(cosh t − 1)/(2(cosh β − cosh t)).
  * [β/2, β]: \|W\| ≤ 4B/β², and ∫₀^β(k/K₀ − 1)dt ≤ ∫₀¹((1−x)^{−½} − 1)dx/x = 2 log 2 exactly (substitute x = (cosh t−1)/(cosh β−1);
    dt/dx ≤ 1/x). This makes the author's 8B log2/β² term exact, for every β.
  * [β, ∞): K₀∫_β^∞B/t² = K₀B/β.
* **Result [certified arithmetic].**
  * β ∈ [30, 64], against LRINF's certified cells: (2/Λ)\|I_{r̃}\| ≤ **0.0249·F_cert** (worst cell [30, 30.5]); author: < 0.027F_*.
  * β ≥ 64, against `fcert_tail`: the ratio 2Λ[B + 8B log2/β + Eβε₁]/((1−2.1·10⁻⁷)A₁ℓ) decreases in β, so ≤ its β = 64 value
    **0.01443** (author: < 0.015); limit 0.01328. A₁ = 2ψ₂′(0) − ψ₃′(0) = 327.149, ℓ = 64L(64) = 0.87401.
  * The extracted piece r(0)sech(πν) contributes +(2/Λ)r(0)/(2cosh(β/2)) > 0 (and, at finite mass, r(0)/Λ·w³𝒱_b((1−z)^{−5/2}) ≥ 0,
    all coefficients nonnegative — consistent with the author).
* **Consequence.** In the FULLPOS3/5 h-expansion route the m = 0 realised-residual term G^r is ≥ −0.025F_cert for every β ≥ 30, uniformly:
  the β² obstruction of FULLPOS5 §3.2 is removed **for the m = 0 term**. What Lemma R does not cover (author agrees): the finite-mass
  (ladder) transfer of r̃ and the high-frequency residual pieces in the local terms h^mP_m(ϑ)G^r, which in FULLPOS5 grow like β.

### R.3 How Lemma R sits in the contour (tail) route
In the contour route (§§1–3 below) the realised f is never split into ψ + r: the main term
𝓜(L) = √w[f(0)/Λ + (2/Λ)∫₀^{3/4}W e^{−Ly}dy] uses the realised f exactly (W ≥ 0 certified by CONTOUR, f(0) = r(0) > 0 by EPSW), and the
finite-mass transfer is the factor (1 + D) with D the ladder deviation of a pure power, applied to the whole f. So the contour route
needs neither Lemma R nor any bound on r: the "higher-order finite-mass transfer of r̃" that Lemma R leaves open is, in this route,
the D-bound of §2, which acts on f as a whole.

---
