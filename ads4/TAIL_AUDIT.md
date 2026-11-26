# TAIL_AUDIT: independent audit of TAIL_CERT.md (the large-impact-parameter tail L ≥ 20 at Δ_g = 1000)

*26 Nov 2026. Independent audit agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
* *Code:* `ads4/code/tailaudit/`, outputs in `cache/` (git-ignored: `*`, `!.gitignore`).
* *Re-runs of the tailcert scripts:* done in a private mirror, `tailaudit/cache/mirror/` (`ta_rerun.sh`). Nothing in
  `code/tailcert/cache/` was overwritten.
* *Compute:* ≤ 3 processes, `OMP_NUM_THREADS=1`, `ulimit -v 4000000`, about 1 CPU-hour in total.

Grades used here:
* **correct**: every step checked, and no change is needed.
* **correct with fix**: the conclusion and the numbers stand, but a proof line is wrong or missing. The fix is given.
* **wrong**: the claim fails.
* **inherited / unverifiable here**: taken from another document and not re-derived. Where possible it is spot-checked numerically.

Severity of problems:
* **high**: invalidates the tail claim.
* **medium**: the claim needs a real repair.
* **low**: a proof-wording gap with a straightforward fix; no number changes.
* **info**: cosmetic, or a statement of scope.

---

## 0. Verdict

**"Tail closed for L ≥ 20" is confirmed: H_ang,b ≥ (1 − ρ)𝓜(L) with ρ ≤ 0.066, for every heavy state with L ≥ 20, every a ≥ 998
(integer or not, up to a → ∞), every spin, b ∈ {0, 4}.** This is conditional only on the inherited inputs listed in §0.3.

What the audit found:
* No high-severity or medium-severity problem.
* Four low-severity proof-wording gaps, each with a fix that changes no number: problems 1–4 in §0.2.
* The certified arithmetic reproduces digit for digit:
  * k = 2 part: ρ ≤ 8.2251·10⁻⁵ (L ≥ 20) and 4.8480·10⁻³ (L ≥ 18);
  * k = 3 part: 0.06501 (L ≥ 20) and 0.6728 (L ≥ 18);
  * remainders: 10^−435.5 / 10^−264.0 (k = 2) and 10^−652.2 / 10^−480.7 (k = 3);
  * cutoff, Lemma R, and the E1–E4 identities.
* Every lemma's key inequality was evaluated directly at adversarial states and holds:
  * the states: a = 998, 998.5, 1234.567, 10⁴, 10⁶; b = 0 and 4; L = 18, 20, 25, 40;
  * complex β on the detour, the vertical piece and the leg, and the k = 3 frequencies.
* **The binding piece is very conservative.** With the actual multiplier in place of Theorem P's bound, the same absolute-value
  k = 3 term is **4.8·10⁻⁵ 𝓜** at L = 20, against the certified 0.065 (§5.6).
  * The one sharp constant is in Lemma NB: the bound sits only 0.05 % above the true value at a = 998.
  * Lemma NB enters ρ only at the 10⁻⁶ level, so this sharpness has no effect on the result.

### 0.1 Per-lemma verdicts

| item (TAIL_CERT §) | verdict | numerical spot check (this audit) |
|---|---|---|
| §1 E1–E3 (ladder, first-moment operator on powers, three-term recurrence) | **correct** (exact; also at non-integer a) | E2, E3 at 4 non-integer-a states × 4 complex β: rel. dev. ≤ 2.5·10⁻⁶⁰ |
| §1 E4 (branch formula) | **correct** | reproduced (15 digits, 3 points) |
| §1 E5 / H/w³ = ℒ(G_A + τ𝒟_bG_A + τσ_bG_D) | **inherited** (FULLPOS3 §4.1(a), ANGULAR_ROUTE §32) | not re-derived |
| §1 per-j contour deformation, arc term f(0)X_j(0)/2, segment W = −2B(iy)Im f(iy) | **correct** | real-axis integral = arc + segment + horizontal line to 8·10⁻⁹ with the realised f (§5.4); W inside every sampled cb_w cell |
| §2.1 Lemma NB | **correct with fix** (problem 1) | holds at all 120 tested (state, β); ratio to bound 0.9995 at a = 998 (sharp), 0.78 at β = ¼ + i/4, ≤ 0.01 for a ≥ 10⁴ |
| §2.2 Theorem P (representation, V-window, regimes A/B, real tail, reflection beyond 1146) | **correct with fix** (problem 2) | representation exact to quadrature precision at non-integer a; \|Ŝ\| ≤ bound at every tested point; worst ratio 0.948 (σ = 9/2, t = ¼), typically 10⁻²–10⁻⁴ |
| §2.3 recurrence bound for the connection terms | **correct** (doc's \|β\| ≤ 1500 is narrower than the use, problem 3) | constants hold up to \|β\| = 6200: max 2.9·10⁻⁷ ≤ 10⁻⁶, min w²\|D\|/(3.9a²) = 1.026 |
| §3 assembly (k = 2 contour, j < 200) | **correct** (reproduced: 8.2251·10⁻⁵) | CONTOUR rows recomputed independently (§5.5) |
| §4 Lemma Loc + K for R₂₀₀ | **correct** (constants have large slack) | I₁, I₂ ≤ 0.40 × the stated bounds; Lemma B and the \|ν\| = 0.4 Pfaff bound ≤ 0.52 and 0.22 of the stated bounds |
| §4 weight bound sup\|f_θ\| ≤ 21.34 on [0, 1300] | **correct with fix** (problem 4) | — |
| §5 k = 3 piece | **correct** (reproduced: 0.06501; σ_b handled by absolute value) | true-multiplier version: 4.8·10⁻⁵ |
| §6 threshold table (L ≥ 20: 0.066; L ≥ 18: 0.68) | **correct** (reproduced) | — |
| §R.2 Lemma R (side item, not in the tail chain) | **correct** (reproduced: 0.0249, 0.01443) | Mehler-step inequalities re-derived |
| §R.2 cutoff correction, \|f_θ − f\| ≤ 2\|θ(ν−i) − 1\|\|g\| | **correct** (the right form is used everywhere) | sup over [1, Λ+80] of the cell enclosure = 3.6·10⁻²⁴ (doc: ≤ 10⁻²¹) |

### 0.2 Problems found

| # | severity | where | problem | fix / effect |
|---|---|---|---|---|
| 1 | low | §2.1 Lemma NB | The bound −log E_k ≤ (1+ε)(d/c)hS₁(k) is justified by "k ≤ n/2 (hu_j/c ≤ ¼)". That only gives 1 + ε = (1 − u_j/X)⁻¹ ≤ 4/3, not 1.0001. The stated constant is sharp: the computed \|Ŝ − 1\| equals 0.9995 × the bound at a = 998. So the proof must really deliver ε < 5·10⁻⁴. | Use −log E_k ≤ (d/c)hΣ_{j<k}u_j/(1 − u_j/X). Then the relative excess is E[S₁(K)u_K/X]/E[S₁(K)] ≈ (β+3)(β+4)/(4a²) ≤ 2·10⁻⁵ (β ≤ 9/2, a ≥ 998). Markov at 20th order handles K > n/2. Even with the literal 4/3, ρ changes by ≈ 10⁻⁶. |
| 2 | low | §2.2 Theorem P, steps 1–2 | Two statements are integer-a statements: "g = A − T, T = Σ_{K>N}(…)" and the Beta–binomial duality P(V ≤ v) = P(Bin(M,v) ≥ N+1). For generic real a, two things go wrong: E_V[((1−λ_Vq)/(1+q))^a] is not defined for V < p, and the q-deformation crosses the branch point q = 1/λ_V when V is far outside the window. | The exact statement for every real a > 0 is **g(q) = E_V[(1 − p/V)₊^a]**, V ~ Beta(N+1, n+s) (verified here to 10⁻⁶⁰ at non-integer a). The repair has three parts. (i) Deform only the in-window V, where λ_V ∈ 1 ± 7√w (measured max: 0.46 × 7√w). (ii) Keep the out-of-window V on the real axis, where 0 ≤ (1 − p/V)₊^a ≤ 1. (iii) Take the window tail from Beta sub-Gaussianity, P(\|V − EV\| > η) ≤ 2e^{−2(M+2)η²} (Marchal–Arbel), which replaces the duality. With η² = (0.75a + πτ/2 + 60)/M the contribution is ≤ 2e^{−1.5a−πτ−120}·(2a)^σ\|1/Γ(β)\|q₁^σ/σ ≤ e^{−1.5a−πτ/2−110}(2a)^σ ≪ 10⁻³⁰. No number changes. |
| 3 | low | §2.3 vs `tc_assemble.py` | The text restricts the recurrence bound to \|β\| ≤ 1500. The code applies it on the leg tail [1146, 6000), with j ≤ 199, i.e. \|β\| ≤ 6200. | The constants hold there: \|β(s−2β−1)\|/(w\|D\|) ≤ 2.94·10⁻⁷ and w²\|D\| ≥ 1.026 × 3.9a², at L ≥ 18 and a ∈ {998, 998.5, 10⁴, 10⁶, 10¹²}. Change the stated range. |
| 4 | low | §4, `tc_remainder.py` | sup\|f\| ≤ sup\|ψ\| + sup\|r\| (+ sup\|f_θ − f\|) is used on [0, 1300]. EPSW's certified sup\|r\| covers only the panels on [0, ⌈Λ⌉+61] = [0, 1060]. | On (1060, 1300], f is ≤ e^{−500} by CONTOUR's Gaussian tail formula (all RBF centres ≤ 1029, widths ≤ 0.9), so the bound is true. Cite the tail formula. |
| 5 | info | §0, §6 | "L ≥ 20 ⟺ β_c ≥ 21.4". In fact L_d = 2 log cosh(β_c/2), so L = 20 corresponds to β_c = 21.386. | The implication actually used, β_c ≥ 21.4 ⇒ L_w ≥ L_d ≥ 20.014, is in the safe direction. Cosmetic. |
| 6 | info | ANGULAR_ROUTE §86 | It quotes Lemma NB as 1.0003(β)₃/(12a²). TAIL_CERT's constants give 1.0001·4·0.2501 = 1.0005; the code uses 1.0005. | Cosmetic. |
| 7 | info | §2.2 | The quoted "[DV] fixed-twist truth 1.0–5.6 for t ≤ 998" is not reproduced exactly. Here: \|Ŝ(¼+it)\| ∈ [1.000, 4.13] for t ≤ 998 at a = 998; it is 7.6 at σ = 5/4, 29.6 at σ = 7/2 and 53.7 at σ = 9/2 (t = 998). | Immaterial: the certified bounds exceed the computed values by ≥ 10³ there. |
| 8 | info | §5 | k = 3 is binding only because Theorem P's regime-B arc (a ≈ 998, ν > 519) is lossy. The certified X₃-weighted integral is 329; with the actual multiplier it is 0.245. | ρ_{k=3}(20) ≈ 4.8·10⁻⁵ with the true multiplier, against 0.065 certified. The certificate is valid but far from optimal, and L_tail could be lowered substantially. |

### 0.3 Inherited inputs (not re-derived here)
1. **The physical action.**
   * The formula: H/w³ = ℒ(G_A + τ𝒟_bG_A + τσ_bG_D), with G_A built from the realised f (θ-cutoff, χ_re) and G_D from the ideal ψ₃.
   * The generating-function identity E5.
   * Sources: FULLPOS3 §4.1(a) [P] and ANGULAR_ROUTE §32. TAIL_CERT also only cites these.
2. **CONTOUR_BOUNDS Arb tables.** These are cb_w, cb_short, cb_vert, cb_horiz_2.25, and the tail formula.
   * They are used as certified.
   * Spot-checked here in double precision (§5.5) and consistent.
   * Analyticity of f in Re ν > 0 rests on CONTOUR's singularity list, which is consistent with the exact kernel polynomials.
3. **EPSW_CERT constants.** r(0), sup\|r\|, ε_W and B, taken from `ew_cert_1000_re_rep.json`. The tail chain uses only r(0) and sup\|r\|.

---

## 1. Reproduction of the certified numbers (`ta_rerun.sh`, mirror in `cache/mirror/`)

Every tailcert script was re-run from a copy of the code. The read-only inputs (contour, epswcert, lrtight, prl) were symlinked.

| output | TAIL_CERT | re-run | match |
|---|---|---|---|
| E1–E3 exact, E4 | True; 3 points to 15 digits | identical | ✓ |
| Theorem P / NB spot values (`tc_power.py`) | 2.266, 3.078, 8.891, 20.17, 55.41, 431, 3.09·10⁴; NB 3.2961·10⁻⁶ | identical | ✓ |
| ehat_d tables (L₀ = 20, 18) | — | byte-identical (`cmp`) | ✓ |
| k = 2 assembly, L ≥ 20 / 18 | 8.2251·10⁻⁵ / 4.8480·10⁻³ | identical (only timings differ) | ✓ |
| leg Σ∫\|fB\|(\|X₀\| + j-sum) | 2.5639·10¹² (L₀ = 20) | identical | ✓ |
| k = 2 remainder | 10^−435.5 / 10^−264.0 | identical | ✓ |
| k = 3 | 3.2891·10², 0.06501 / 5.0410·10², 0.6728; remainder 10^−652.2 / 10^−480.7 | identical | ✓ |
| cutoff | δE ≤ 1.659·10⁻²³, δB ≤ 4.387·10⁻²², ∫e^{0.02ν}\|f_θ−f\| ≤ 1.333·10⁻¹⁴ | identical | ✓ |
| Lemma R | 0.024895, 0.014427 | identical | ✓ |

Total: ρ(L ≥ 20) ≤ 8.23·10⁻⁵ + 0.06501 + 10^−435 + 10^−652 = 0.06509 ≤ 0.066.

The segment deviation δ_seg ≤ 3.3·10⁻⁶ is folded into 𝓜_lb. So H ≥ (1 − δ_seg)𝓜 − ρ_err𝓜_lb ≥ (1 − 0.0651)𝓜, which is the
stated inequality with the full 𝓜, including the f(0)/2 term.

The [DV] scripts were not re-run: `tc_lownu_dv.py`, `tc_author_dv.py` and `tc_power_dv.py`. They are not in the certified chain.
Their role is covered by §5.1.

---

## 2. Line-by-line checks of the [P] lemmas

### 2.1 Framework (§1)

**Identities.** E2, E3, the Theorem P representation and Lemma Loc's B-representation of ℒ were re-checked at non-integer a
(`ta_identities.py`, 60 digits). The states were (n, a, b) = (7, 3.37, 0), (11, 5.5, 4), (20, 12.25, 0), (15, 17.9, 4), with β
complex, including Re β < 0 for E2 and E3. `tc_exact.py` only tested integer a. Maximum deviations:

| identity | max deviation |
|---|---|
| E2 | 1.4·10⁻⁶⁰ |
| E3 | 2.5·10⁻⁶⁰ |
| Theorem P representation | 4·10⁻¹⁵ (quadrature) |
| g = E_V[(1−p/V)₊^a] | 2·10⁻⁶⁰ |
| ℒ(F) = E_B[[zⁿ](1−z)^{−a−1}F(Bz)]/C(N,a) | 3·10⁻⁶¹ |

**Per-j contour.** Each piece of the argument was checked.
* *Singularities in the swept region.* There are none in {Re ν > 0, 0 < Im ν ≤ 9/4}:
  * B has poles at im and zeros at i(m+½);
  * c_j has poles at −i(m+1)/2;
  * f has poles at ±3i/2, ±5i/2, 9i/2 and branch points at 0, ±2i, all on the imaginary axis;
  * Log(ν′²+1) = Log ν + Log(ν−2i) holds in Re ν > 0, and the conjugate part uses Log ν + Log(ν+2i);
  * X_j is entire, since S(β) is a polynomial in β.
* *Closing at Re ν → ∞.* On 0 ≤ Im ν ≤ 9/4 the Gaussian RBF decay e^{−(Re ν − 1029)²/1.62} beats the e^{π\|ν\|}·poly growth of X_j,
  which follows from the reflection bound together with the recurrence.
* *Integrability at ν = 0.* f = F₁ + F₂ log ν with F₂(0) = 0 (EPSW §1). The j-th integrand on the real axis is integrable at 0,
  because 2Re[h(0)/(2πiν)] = 0 for real h(0) = c_j(0)w^{1/2+j}X_j(0).
* *Arc term.* It is (iπ/2)·f(0)c_j(0)w^{1/2+j}X_j(0)/(2πi) = f(0)c_j(0)w^{1/2+j}X_j(0)/4, giving the stated
  (2/Λ)·f(0)c_j(0)w^{1/2+j}X_j(0)/2 after 2Re.
* *Segment.* c_j(iy), B(iy) and X_j(iy) are real, so the integrand is W(y)c_j w^jX_j√w e^{−Ly}.

**Numerical check of the bookkeeping** (`ta_contour_bookkeeping.py`). This uses an independent numpy evaluation of the realised f
(repaired JSON, θ included) with X ≡ 1. The real-axis value 2Re∫₀^∞fBw^{1/2−iν}dν equals arc + ∫₀^{0.3}W√w e^{−Ly}dy + the
horizontal line at Im ν = 0.3, to 7.3·10⁻⁹ at L = 3 and 8.1·10⁻⁹ at L = 4. W(y) computed this way lies inside CONTOUR's certified
cb_w cell enclosure at all 100 sampled cells. The arc term (f(0) ≈ 10⁻¹³) is below this test's resolution; its factor ½ is the
standard quarter-residue.

### 2.2 Lemma NB (§2.1): correct with fix (problem 1)

The following parts were checked and are correct:
* the decomposition Ŝ = (w/d)^βΣπ_k(β)E_k with Σπ_k = 1 for complex β;
* E_k ∈ [0,1]: each factor (X − u)/(Y − u) ≤ X/Y;
* h/d² = N(N+s)/(a²(N+n+s)²) ≤ 0.2501/a² in the tail (a/N ≈ w/2);
* S₁(k) = (k)₃/3 + (s+1)(k)₂/2, consistent with FULLPOS3 §4.1(b); the code's O(d) coefficient dominates (s+1)(β)₂d/(2c);
* C_τ = exp(τ²(σ⁻² + σ⁻¹)/2) from Π\|1 + iτ/(σ+j)\|;
* (w/d)^β = (1 − 20/T)^β for b = 4, and w = d exactly for b = 0.

The one gap is the (1+ε) step (problem 1). The per-cell use of max(NB(σ₀), NB(σ₁)) on vertical cells is valid because
NB(σ, ¼) is increasing on [¼, 7/2]: d log NB/dσ = −(2σ⁻³ + σ⁻²)/32 + Σ1/(σ+i) > 0 already at σ = ¼.

Numerics (`ta_power_check.py`, `ta_power_nonint.py`). The independent evaluation of Ŝ is described in §5.1.

| state | β | \|Ŝ − 1\| | bound | ratio |
|---|---|---|---|---|
| a = 998, b = 0/4, L = 18, 20, 25 | 5/2 | 3.2944·10⁻⁶ | 3.2961·10⁻⁶ | 0.9995 |
| a = 998, L = 40 | 5/2 | 3.2910·10⁻⁶ | 3.2961·10⁻⁶ | 0.9985 |
| a = 998 | 7/4 + i/4 | 1.5350·10⁻⁶ | 1.5537·10⁻⁶ | 0.988 |
| a = 998 | ¼ + i/4 | 8.54·10⁻⁸ | 1.10·10⁻⁷ | 0.776 |
| a = 998 | 7/2 | 7.2476·10⁻⁶ | 7.2514·10⁻⁶ | 0.9995 |
| a = 998.5, 1234.567 (non-integer) | 5/2, 7/4 + i/4, ¼ + i/4 | ≤ bound | — | ≤ 0.999 |
| a = 10⁴, 10⁶ | all | ≈ (σ)₃/(12a²) | (a = 998 constant) | ≤ 0.01 |

The fixed-twist value −(β)₃/(12a²) quoted in the document is confirmed.

### 2.3 Theorem P (§2.2): correct with fix (problem 2)

The following were checked line by line:
* **Representation.** It is correct as a Beta integral followed by Pfaff (terminating in n). It also holds for non-integer a, numerically.
* **Prefactor.**
  * w(N+s) ≤ 2a: w(N+s) = a(N+n+s)(N+s)/T ≤ a(2N − a + s)/N ≤ 2a.
  * The Binet factor uses x = N+s ≥ 998e^L and |β|² + |β| ≤ (t₂+5)².
* **Window.** |λ − 1| ≤ a/N + 4.0001η ≤ w/2 + 4√(1.0001w(0.75 + 1860/998)) = w/2 + 6.47√w. Measured at the window edges:
  ≤ 0.46 × 7√w.
* **Imaginary leg.** |(1−λiy)/(1+iy)|^a ≤ exp(a(λ²−1)₊y²/2), and λ² − 1 ≤ 3|λ − 1|.
* **Regime-A arc.**
  * −atanh(x) ≤ −x.
  * Jordan's inequality π/2 − φ ≤ (π/2)cos φ.
  * Together these give f ≤ cos φ[uπ/2 − 2u/(1+u²)] ≤ 0 ⟺ u² ≤ 4/π − 1 ⟺ u ≤ 0.5227.
  * The δ-term is aδ/(1−ρ) = |λ−1|τ/(1−ρ) ≤ |λ−1|τ₂/0.48.
* **Regime B.**
  * Ehat_d is a rigorous Arb enclosure on φ-balls that cover [0, π/2].
  * Monotone in u, since ∂_u f = π/2 − φ ≥ 0.
  * The a-cells use the conservative ends:
    * the larger u = t₂/a₁;
    * Y₁ with a₂ and c₁ with a₁;
    * the loss term a₂·E (for E > 0);
    * Γ(σ, ·) at a₁.
* **Real tail.**
  * (1−λq)/(1+q) ≤ e^{−1.99q} on [q₁, 0.92] needs |λ−1| ≤ (0.01 + 2q²/3)(1−q), i.e. |λ − 1| ≤ 0.0099, which holds at L ≥ 18.
  * g ≤ 0.0516^a beyond 0.92 kills e^{πτ/2} ≤ e^{1800}.
  * Γ(σ, X) ≤ X^{σ−1}e^{−X} holds for σ ≤ 1.
* **Reflection bound beyond 1146.** |1/Γ(¼+it)| ≤ Γ(¾)cosh(πt)/π. For σ = 5/4 the extra factor is ¼/|¼+it| ≤ 1. Both are as coded.

The gap is problem 2: the integer-a form of g = A − T and the Beta–binomial duality. Its repair is given in §0.2 and is
rigour-neutral.

**Numerics.**
* *Method.* An independent quadrature of the exact representation, at finite n, with an exact Beta average over V. It is
  validated against the direct ladder sum at n = 3·10⁵ to 3·10⁻¹².
* *States.* a ∈ {998, 10⁴, 10⁶} (and 998.5, 1234.567), b ∈ {0, 4}, L ∈ {20, 25, 40} (and 18).

\|Ŝ(σ+it)\| against the certified K_osc:

| σ \ t | 0.25 | 10 | 45 | 300 | 500 | 800 (regime B) | 998 (regime B) |
|---|---|---|---|---|---|---|---|
| ¼: max \|Ŝ\| / bound | 1.000 / 1.80 | 1.000 / 8.5 | 1.0006 / 18.4 | 1.030 / 51.6 | 1.094 / 74.8 | 1.368 / 478 | 4.13 / 2.8·10⁴ |
| 5/4 | 1.000 / 1.59 | 1.0001 / 7.28 | 1.0011 / 15.9 | 1.054 / 46.6 | 1.173 / 69.0 | 1.712 / 505 | 7.61 / 5.2·10⁴ |
| 7/2 (k = 3) | 1.000 / 1.08 | 1.0001 / 26.1 | 1.0023 / 60.0 | 1.111 / 182 | 1.372 / 261 | 2.84 / 8772 | 29.6 / 2.2·10⁵ |
| 9/2 (k = 3) | 1.000 / 1.054 | 1.0001 / 48.1 | 1.0028 / 115 | 1.138 / 353 | — | 3.55 / 3.1·10⁴ | 53.7 / 4.1·10⁵ |

The regime-B values were computed on the ρ = 0.92 path and agree with the regime-A path to all printed digits (path independence).
Finite-n and fixed-twist values agree to 4 digits, so the finite-n correction is ≲ 10⁻⁴ here.

### 2.4 Recurrence (§2.3): correct (scope note, problem 3)

From E3:

  Ŝ(β−1) = (T/D)Ŝ(β) + β(s−2β−1)/(wD)·Ŝ(β) + β(β+1)/(w²D)·Ŝ(β+1),  D = N(N+s) − β(β−s).

Constants checked at L ∈ {18, 20, 30}, a ∈ {998, 998.5, 10⁴, 10⁶, 10¹²}, b ∈ {0, 4}, |β| ≤ 6200:

| quantity | measured | required |
|---|---|---|
| \|T/D − 1\| | ≤ 2.2·10⁻¹⁵ | ≤ 10⁻⁶ |
| \|β(s−2β−1)\|/(w\|D\|) | ≤ 2.9·10⁻⁷ | ≤ 10⁻⁶ |
| min w²\|D\|/(3.9a²) | 1.026 | ≥ 1 |

The code iterates the triangle inequality exactly, with |β| ≤ |α| + j − 1. It does not use the 1.5^j envelope. The early break in
`jsum` is safe because the term ratio is ≤ w\|ν\| ≪ 1.

**Dressing constants.** The required ratios are all ≤ 1:

| ratio | max measured |
|---|---|
| (T⁻¹/w)/(w/(3.9a²)) | 0.975 |
| T⁻¹/(w²/(3.9a²)) | 0.975 |
| (T⁻¹/w)/(0.2515w/a²) (k = 3) | 0.994 |

The |A₁| and |A₀| majorants in `jsum` dominate both channels.

### 2.5 Assembly (§3): correct

* **Main term.** `main_lb` uses minW·(1 − δ_seg)·∫e^{−Ly} with δ_seg including the dressing. The arc term f(0)X₀(0)/2 > 0 is dropped.
* **Errors.** The segment j ≥ 1 term uses maxW and has no extra factor 2, since W already contains it. The off-segment terms carry
  2·e^{−L Im ν}·∫|fB|·(|X₀| + j-sum), with Im ν at the lower cell end on the vertical piece.
* **Leg tail.** [1054, 1146) uses CONTOUR's ∫_T^∞ per unit cell × K_osc at [T, T+1]. [1146, 6000) uses the reflection bound
  e^{π(T+1)}; it underflows to 0 because of e^{−(T−1029)²/1.62}.
* **Monotonicity in L.** All errors decrease in L and 𝓜_lb decreases, so Err(L_k)/𝓜_lb(L_{k+1}) is valid.
* **Beyond L = 80.** The segment j ≥ 1 and residue ratios printed by the code are computed at w(L₀). They are therefore flat or
  slightly rising (res: 10⁻²⁰ → 3·10⁻²⁰). The true terms carry w(L)^j, so the [P] argument beyond 80 is right.
* **Hand recomputation at L = 20.**
  * Leg: 2e⁻⁴⁵·2.56·10¹²/𝓜_lb(2.62·10⁻³) = 5.6·10⁻⁵, matching the code.
  * Short detour: ~2·10⁻⁵, against the code's 1.6·10⁻⁵ (the code integrates per cell).

### 2.6 Lemma Loc and the remainder (§4): correct (weight-range note, problem 4)

**Proof.** Every step was re-derived:
* the B-representation (verified exactly);
* Cauchy on |z| = ρ, with |Bz| ≤ ρ < ρ₂;
* |1 − Bz|^m ≤ 2^{m−1}(|1−z|^m + (1−B)^m);
* (1/2π)∫|1 − ρe^{iθ}|^{−γ} ≤ (1−ρ)^{1−γ}√πΓ((γ−1)/2)/(4√ρΓ(γ/2)), obtained by the substitution θ = π(1−ρ)u/(2√ρ);
* the dressing, by Cauchy on circles of radius ½min(|1−ζ|, 1−ρ). These stay in |ξ| ≤ ρ₂, and the factor ≤ 1.5^m loses one power
  of |1−ζ|. A direct count gives a constant ≤ ~30·1.5^m T⁻¹, well inside 40·1.5^m.

**Numerics.** I₁ and I₂ were computed exactly via ₂F₁(γ/2, γ/2; 1; ρ²):
* states: L ∈ {8, 20}, a ∈ {998, 1234.5, 10⁴}, m ∈ {2, 50, 196};
* I₁ ≤ 0.396 × the stated 2.53(a/N)^m;
* I₂ ≤ 0.396 × the stated bound;
* 2a/N = 0.9994 × (1.0006w).

**K(ν).**
* Near |1−ξ| ≤ ½: the branch-series tail with |arg(1−ξ)| < π/2.
* Far: FULLPOS3 Lemma B, |P_{−½+iν}(cosh ζ)| ≤ (π/2)e^{ν|Im ζ|}. It holds on 1600 random samples (|ξ| up to 0.999, ν up to 60),
  with max ratio 0.51.
* |ν| ≤ 0.4: maximum principle, since R is analytic in |ν| < ½ (the ν = 0 poles of B(±ν) cancel pairwise). The Pfaff bound
  2^{0.9}e^{0.2π}/(1−|ξ|) holds with max ratio 0.22, and 1/(1−ρ₂) ≤ 4.0002/((1−m/998)w) gives the e^{L} factor.

**Code.** Checked:
* both branches (×2);
* the max over near and far;
* ∫_T^∞|f|·K(T+10) summed over steps of 10;
* the w^{5/2} prefactor;
* the "decreasing in L" claim, since e^{−195L} beats e^{L}.

The weight bound is problem 4.

### 2.7 k = 3 (§5): correct

* **Formula.** |H^{(3)}|/𝓜 ≤ (T⁻¹/w)Λ³∫|ψ₃||B(Λp)|X₃dp/I_W, with T⁻¹/w ≤ 0.2515w/a² (measured ≤ 0.994 of it).
* **Weight.** |ψ₃B| = |ψ₃/p|·|νB(ν)|/Λ, with |νB|² = ν coth(πν)/(4π) exactly. ψ₃/p is enclosed by an exact fmpq Taylor shift.
* **Multiplier.** X₃ comes from K_osc(7/2), K_osc(9/2) and the recurrence. It includes the dressing, which G_D does not have, so
  this is conservative.
* **Sign.** σ_b is irrelevant, since the term is bounded in absolute value.
* **Monotonicity in L.** The ratio decreases for L ≥ 20, because d log(e^{−L}/I_W)/dL = −1 + ⟨y⟩ ≤ −¼.
* **Reproduced value.** 0.06501.
* **Remainder.** 10^−652. The extra factor Λ and the sup|ψ₃| ≤ Σ|coeff| in the code are conservative.

### 2.8 Lemma R and the cutoff (§R.2): correct

**Lemma R.**
* *Inputs.* E_in = ε_W + r(0)/2 + δE and B_in = B + π²r(0) + δB. Here ∫|(sech πν)″| ≤ π², with sech″ = sech − 2sech³.
* *Mehler step.* All three pieces were re-derived:
  * on [0, β/2], (1−x)^{−½} − 1 ≤ x/(2(1−x));
  * on [β/2, β], dt/dx ≤ 1/x and ∫₀¹((1−x)^{−½} − 1)dx/x = 2 log 2;
  * on [β, ∞), K₀B/β.
* *Reproduced values.* 0.024895 and 0.014427.
* Lemma R is not part of the tail chain.

**The cutoff inequality (task item 4).**
* *Where it is used.* f_θ − f appears in TAIL_CERT in only three places, and each uses the correct form
  |f_θ − f| ≤ 2|θ(ν−i) − 1||g|:
  1. Lemma R's δE_θ and δB_θ;
  2. the remainder's weight bound;
  3. the main term's f(0), where f_θ(0) − f(0) = O(e^{−12000}).
* *How `tc_cutoff.py` encloses it.*
  * θ − 1 = −½erfc((ζ+ν_c)/10) − ½erfc((ν_c−ζ)/10), with |erfc(x+iy)| ≤ e^{y²}erfc(x). The proof shifts the path to t + iy.
  * |g| is bounded through |QR₋|·(RBF disc bound) + |YlR₋|.
  * The second derivative is bounded by Cauchy on radius ρ − h.
* *No |Re g| bound anywhere.* The contour route never uses such a bound, because CONTOUR's f is f_θ itself:
  G = θ(ν′)χ(ν′)R₋(ν′).
* *Sup.* The certified sup over [1, Λ+80] of the cell enclosure is 3.65·10⁻²⁴ (at ν ≈ 1025.5), so the quoted ≤ 10⁻²¹ holds.

---

## 3. Inputs from CONTOUR_BOUNDS and EPSW_CERT (task item 4)

| input | used for | how used | check here |
|---|---|---|---|
| cb_w (7500 cells, minW, maxW; W ≥ 0) | main term 𝓜_lb, segment j ≥ 1 | lower bound with minW(1−δ_seg), errors with maxW | 100 cells: numpy W inside the enclosure |
| cb_short / cb_vert / cb_horiz_2.25 rows (∫\|fB\| per cell) | detour, vertical, leg | ∫\|fB\|·sup-cell \|X\| | 25 / 212 / 164 rows recomputed: max computed/certified = 0.961 / 0.99981 / 0.984 |
| cb_lines.tail (Gaussian tail formula) | leg t > 1054, remainder weight beyond 1300, cutoff tail | ∫_T^∞\|f\| (θ ≤ e^{y²/a²}), sup\|B\| | formula read: conservative (e.g. (2R₀−7) twice in qbar) |
| EPSW r(0) = 9.6925·10⁻¹⁴ | residue j ≥ 1 term (as \|f(0)\|), positivity of the dropped arc term | upper bound only | repaired weight, r(0) > 0 ✓ |
| EPSW sup\|r\| = 5.556·10⁻⁶ | remainder weight on [0, 1300] | sup\|ψ\| + sup\|r\|·1.0001 | valid on [0, 1060] only (problem 4) |
| EPSW ε_W, B | Lemma R only | — | — |

The density is consistent across the documents:
* CONTOUR's f, EPSW's E with θ, and the §32 "E_eff,θ" are the same realised, θ-cut, χ_re density.
* The main term uses the realised f exactly; the contour route never splits f = ψ + r.

---

## 4. Scope: all a ≥ 998 at L ≥ 20, and the overlap (task item 5)

**Every bound is uniform in a ≥ 998.**

| bound | why it is uniform |
|---|---|
| NB constant 0.2501/a² | evaluated at a = 998 |
| Dressing and recurrence | use A_MIN = 998 |
| Theorem P regime A | uniform in a, via \|λ − 1\| ≤ 7√w, w(N+s) ≤ 2a, the Binet bound with x ≥ 998e^{L} |
| Theorem P regime B | only needed for 998 ≤ a ≤ t/0.52 ≤ 2204, and covered by geometric a-cells |
| Lemma Loc | Q_m decreases in a; (a/N)^m ≤ (0.5003w)^m |
| k = 3 | T⁻¹/w ≤ 0.2515w/a² is worst at a = 998 |

* *Numerical support.* Checks at a up to 10⁶ (Ŝ), up to 10¹² (constants), and at non-integer a all pass.
* *Non-integer a.* All exact identities hold. The integer-specific steps of Theorem P are repaired in problem 2.
* *Channels and dressing.* Both channels use the worse constants (s = 9 in NB; |A₁|, |A₀| majorants). The first-moment dressing is
  in every X_j.

**L ≥ 20 covers every such state.** For any J with w = 1 − C/T ≤ e^{−20}, n ≈ 2a e^{L} is forced. Nothing else is assumed.

**Overlap.**
* The cell certificate covers β_c ≤ 64 only at a ≥ a*(β), with a* ≈ 3200–3600 (§§73, 77, 84). For β ∈ [27.5, 64] the laptop scan
  was still running at §84.
* The tail does not need the overlap: for L ≥ 20 it alone covers every a ≥ 998. TAIL_CERT §6 says this ("the tail covers every
  a ≥ 998").
* The claimed "wide overlap" is real only for a ≥ a*(β), and only on completed cells.
* The open region is unchanged by this audit: the sliver 998 ≤ a < a*(L), 0.04 ≲ L < 20, as in §86.

**Coordinates.**
* L_d = −log d = 2 log cosh(β_c/2) ≤ L_w = −log w for b = 4, and they are equal for b = 0.
* So β_c ≥ 21.4 ⇒ L_w ≥ 20.014, which is inside the tail.

---

## 5. Numerical spot checks: method and files

### 5.1 Independent evaluation of Ŝ
`ta_power_check.py` computes

  Ŝ(β) = w^βκ(2a)^{−β}∫_path x^{β−1}E_V[((1 − λ_V x/2a)/(1 + x/2a))^a]dx,  V ~ Beta(N+1, n+s).

* E_V is taken by 24-point Gauss–Hermite on the exact Beta mean and variance.
* The path runs i[0, Y₁] → arc → real axis. Y₁ = 2τ in regime A and 1.84a in regime B.
* The prefactor is evaluated at 90 digits, and the quadrature in double precision.
* There is no cancellation in regime A. In regime B the arc amplification is ≤ e^{6} for t ≤ 998.

This uses the exact representation, but none of Theorem P's estimates. Validation against the direct ladder sum (n = 3·10⁵,
a = 998, L = 5.02):

| β | discrepancy |
|---|---|
| ¼ + 10i | 3.1·10⁻¹² |
| 5/2 | 2.5·10⁻¹⁴ |
| 5/4 + 3i | 8.7·10⁻¹⁴ |

### 5.2 Lemma NB, Theorem P
See §§2.2–2.3.
* 120 (state, β) NB checks and ~600 Theorem-P evaluations; all pass, with assertions in the code.
* Non-integer a and L = 18 (against the L₀ = 18 table): `ta_power_nonint.py`, all assertions pass.

### 5.3 Constants
`ta_lemmas_check.py` covers four groups:
* the recurrence, the dressing, and the 0.2515 constant;
* the λ-window;
* Lemma Loc's I₁ and I₂;
* the Legendre bounds used for K(ν).

### 5.4 Contour bookkeeping
`ta_contour_bookkeeping.py` (§2.1).

### 5.5 CONTOUR rows
`ta_contour_rows.py`: per-cell ∫|fB| recomputed by 30-point Gauss–Legendre with the independent f. All rows are below the
certified values, which are 0.1–4 % above the computed ones.

### 5.6 How conservative k = 3 is
`ta_k3_truth.py` takes the same formula as `tc_k3.py`, with X₃ replaced by the actual fixed-twist |Ŝ(7/2+iΛp)| at a = 998. Its
maximum is 28.2, reached near ν ≈ a.

| quantity | value |
|---|---|
| ∫\|ψ₃B\| | 0.198 (doc: ≈ 0.2) |
| X₃-weighted integral, actual multiplier | 0.245 (certified: 329) |
| k = 3 ratio at L = 20, actual multiplier (not a bound) | 4.8·10⁻⁵ (certified: 0.065) |

---

## 6. Files

`ads4/code/tailaudit/`:
* `ta_rerun.sh`: re-runs every tailcert script in a mirror (`cache/mirror/ads4/code/tailcert/cache/*`);
* `ta_identities.py` → `cache/ta_identities.txt`;
* `ta_power_check.py` → `cache/ta_power_check.txt`;
* `ta_power_nonint.py` → `cache/ta_power_nonint.txt`;
* `ta_lemmas_check.py` → `cache/ta_lemmas_check.txt`;
* `ta_contour_bookkeeping.py` → `cache/ta_contour_bookkeeping.txt`;
* `ta_contour_rows.py` → `cache/ta_contour_rows.txt`;
* `ta_cutoff_sup.py` → `cache/ta_cutoff_sup.txt` (Arb);
* `ta_k3_truth.py` → `cache/ta_k3_truth.txt`.

Reproduce (`ulimit -v 4000000; export OMP_NUM_THREADS=1`):

```
./ta_rerun.sh                     # ~10 min, 3 processes
python3 ta_identities.py          # ~1 min
python3 ta_power_check.py         # ~15 min
python3 ta_power_nonint.py; python3 ta_lemmas_check.py; python3 ta_contour_bookkeeping.py; python3 ta_contour_rows.py
python3 ta_cutoff_sup.py; python3 ta_k3_truth.py
```
