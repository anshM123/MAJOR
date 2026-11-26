# ANALYTIC_SLIVER: an analytic route to H_ang,b > 0 on the sliver (Δ_g = 1000)

*26 Nov 2026. Research agent. Running document; unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*New code: `ads4/code/analytic/` (outputs in `cache/`, git-ignored except its `.gitignore`).*
*Grades: **[P]** proof written here; **[certified]** Arb/exact scalar constants; **[DV]** numerics (high-precision, sampled); **[conj]**.*

The sliver: 998 ≤ a < a*(L) (a* ≈ 3100–3600), J ≥ 1001, b ∈ {0, 4}, 0.06 ≲ L < 20 (u = β_c ≥ 0.486).

---

## 0. Verdict (running)

**PARTIALLY CLOSED.** The proof is analytic and certified (modulo review of the new lemmas listed in §9) for **every sliver state with L = −log d ≥ 5**:
* all a ≥ 998;
* every spin;
* b ∈ {0, 4}.

Coverage by method:
* [5, 8] by Route T′ with Lemma Loc′ (§7.3–7.4);
* [8, 16] by Route T′ (§7);
* L ≥ 15.5 by TAIL+E′ (§8);
* L ≥ 20 by TAIL_CERT.

**Open: 0.06 ≲ L < 5.** The obstruction is the turning-point residual cell ν ≈ 1031 (termwise hump × Lemma X). Below L ≈ 4 it is also
Lemma Loc′ and the top line (Lemma T). See §9.

---

## 1. Direction 1: how far down can the TAIL_CERT route go? (quantified) [DV + reading of TAIL_CERT's certified logs]

TAIL_CERT's per-piece budgets (cache/tc_assemble_L18.0.txt, L20.0.txt) scale as follows.

| piece | L = 20 | L = 18 | scaling in L | origin of the size |
|---|---|---|---|---|
| leg (Im ν = 9/4) | 5.6e-5 | 4.6e-3 | e^{−2.25L}·L | ∫\|fB\| on the leg = **5.3e11** (realised f has off-axis Gaussian growth e^{y²/2σ²}, σ ≈ 0.45–0.9); multiplier only ×4.8 |
| vertical Re ν = 1/4 | 1.0e-5 | 1.6e-4 | e^{−0.75…2.25L} | same off-axis growth of f |
| short detour | 1.6e-5 | 6.5e-5 | e^{−0.75L} | |
| k = 3 (real axis, absolute) | 0.065 | 0.673 | e^{−L}L | absolute real-axis bound; true value 4.8e-5𝓜 (TAIL_AUDIT) |
| remainder R₂₀₀ (Lemma Loc) | 1e-435 | 1e-264 | ~e^{−196(L−L*)} | (1.0006w)^{m−1}Q_m with m = 196 |

* **Sharpening Theorem P / k = 3 alone** (k = 3 by the same contour as k = 2, or with |X − 1| instead of |X|) removes the k = 3 obstruction,
  but then the leg binds: ratio 5.6e-5·e^{2.25(20−L)} reaches 0.5 at **L ≈ 16**. With difference multipliers |X − 1| ~ |β|³/(12a²) on the
  leg (t ≲ 45, where all of ∫|fB| sits) the leg gains ~10³–10⁴ and the floor moves to **L ≈ 12–13**. This is the ceiling of the route as built:
  **the obstruction is the realised density's off-axis growth (5.3e11), not the multiplier.**
* **Ideal-weight variant (new, §4).** Splitting f = ψ₂/k₀ + r and moving only the polynomial ψ₂ part onto the contour removes the 5.3e11:
  ∫\|ψ₂B\| over a leg of height Y is ~10³ and grows only polynomially in Y. A termwise j-sum scan (`t9.py`, multiplier envelope
  5·∏(1 + \|β\|²/(3.9a²))) gives leg/F_cert = 2e-6 (L = 3, Y = 10), 1e-16 (L = 4, Y = 10), 2e-22 (L = 5, Y = 10). The leg then no longer binds.
  What binds instead:
  1. **Lemma Loc (j-truncation remainder).** With TAIL's K ≈ e^{1600} and the optimal m: log(bound) = −40 at L = 4.5 (m = 671), −384 at L = 5
     (m = 705), +285 at L = 4. **Floor L ≈ 4.6** (a = 998; better at larger a).
  2. **Theorem P at moderate L.** Its arc bound uses \|λ_V − 1\| ≤ 7√w; at L = 5 the regime-A cost c₂ = δτ/0.48 is ~e^{300} at τ ~ 10³. The leg
     can absorb this only by raising Y (e^{−LY}); below L ≈ 8–10 Theorem P has to be re-derived with the V-average kept signed.
  3. **The j-series hump.** Termwise Σ_j\|c_j(ν)\|d^j ~ e^{\|ν\|d/2}: harmless for L ≳ 4, e^{67} at L = 2, e^{470} at small L. **Intrinsic to the
     (1−z)-connection expansion; no choice of leg height cures it below L ≈ 2.**

**Conclusion for direction 1.** The TAIL route, even rebuilt with the ideal-weight split, cannot go below L ≈ 4.5 (Lemma Loc, j-series hump),
and with the existing Theorem P not below L ≈ 8–10. It cannot reach L ≈ 0.5.

---

## 2. A finite-mass Mehler kernel in t-space does not exist in tame form (dead end, recorded) [DV]

Mehler: F_{3,ν}(z) = (2/π)coth(πν)∫_{β(z)}^∞ sin(νt) q_t(z) dt with q_t(z) = (1−z)^{−5/2}(1 − z/z_t)^{−1/2}/(2 sinh(t/2)), z_t = tanh²(t/2).
q_t is analytic at z = 0 for every t > 0, so k_fm(t) := ℒ(q_t) is well defined (finite sum of Taylor coefficients). `t6.py` at (J, a, b) = (300, 998, 0),
β_c = 0.471: k_fm(t)·2sinh(t/2) = 1.2045 / 1.3724 / 3.4907 at t = 2 / 0.9 / 0.5 (infinite twist: 1.2045 / 1.3724 / 3.5050), but **4.9·10⁶⁴ at t = 0.3**.
Reason: for t < β_c the singularity z_t lies inside \|z\| < c and Σ_k U_k z_t^{−k} ~ exp((2/3)√(c/(hd))·log^{3/2}(c/z_t)). So the "smeared
Mehler kernel" picture is false: the finite-mass kernel is exponentially large in the classically forbidden range t < β_c and the t-space
pairing with W(t) (large near t = 0) is not absolutely controllable.

---

## 3. The canonical finite-mass outgoing kernel [DV, exact identity checked]

For real ν, K(ν) := ℒ(F_{3,ν}) = Σ_{k≤n} U_k g_k(ν) (a polynomial of degree 2n in ν). With α = 5/2 + iν, a₀ = ½ − iν, b₀ = 1 − 2iν, and
B(ν) = Γ(2iν)/Γ(½ + iν)², the infinite-twist outgoing function is (TAIL E4, checked here including complex ν against −(i/π)coth(πν)Q_{−½−iν}):

  κ_∞(ν) = B(ν)Σ_j c_j(ν) d^{j−α} = (1−c)^{−3}P⁺_ν(β_c),  c_j = (a₀)_j²/((b₀)_j j!),  K_∞ = 2Re κ_∞.

**(a) The j-series is only asymptotic at finite mass.** κ_J(ν) := B Σ_{j<J} c_j S(α−j), S(β) = ℒ((1−z)^{−β}) = ₃F₂(β, −n, n+s; −N, N+s; 1).
Since S(−m) ~ C(m, n)U_n for m ≫ n, the terms grow like j^{n−1}: **the series diverges for every state** (`t8.py`: at (n, a) = (5, 4.5),
2Re κ_J = 1.2749, 1.2823, 1.4208, 2.19, 9.19 at J = 10, 20, 30, 40, 59 against the exact 1.27727). Where d is not close to 1 the minimal
term is astronomically small (e^{−600}…e^{−1400}), and the optimally truncated series reproduces the ladder to 1e-61 (`kplus.py`, `t4.py`).
At small c (d → 1) the hump e^{\|ν\|d/2} meets the onset of divergence (j ~ (3/x)^{1/3}) and the series is useless (e.g. (1001, 3200, 0), ν = 990).

**(b) Exact convergent form (Mellin–₃F₂).** With TAIL's Theorem-P representation S(β) = κ(β)∫₀^∞ q^{β−1}g(q)dq, κ(β) = Γ(N+s)/(Γ(β)Γ(N+s−β)),
g(q) = E_V[(1 − p/V)₊^a] (p = q/(1+q), V ~ Beta(N+1, n+s)), and κ(α−j) = κ(α)(−1)^j(1−α)_j/(N+s−α)_j, the j-sum resums inside the q-integral:

  **κ^M(ν) := B(ν)κ(α)∫₀^∞ q^{α−1} g(q) ₃F₂(a₀, a₀, 1−α; b₀, N+s−α; −1/q) dq.**

The integral converges for every ν (integrand ~ q^{α−1+a₀} = q² and ~const at q → 0; g decays at ∞). **Check (`t8.py`, 30 digits, non-integer
a = 4.5, n = 5, ν = 1.3): 2Re κ^M = 1.27726500574794198760375971779 = K(ν) exactly.** So κ^M is a canonical outgoing/incoming split of the exact
finite-mass kernel, K(ν) = κ^M(ν) + κ^M(−ν), with no truncation remainder. (Proof of the identity in general: not written; [DV].)

**(c) Stieltjes form (for Im ν > 3/2).** Using Euler's integral for the ₃F₂ and x = u/q:
  κ^M(ν) = B(ν)(sin πα/π)∫₀^∞ x^{−α−1} ₂F₁(a₀, a₀; b₀; −x) G_N(x) dx,  G_N(x) = E_{u,V}[((x − λ_V u)/(x + u))₊^a],
u ~ Beta(1, N+s−1), λ_V = 1/V − 1; at infinite twist G_∞(x) = x/(x + d). This exhibits the finite-mass operator as a profile G_N replacing a
simple pole at x = −d by order-a poles at x = −u smeared over u. For bounds it has an intrinsic e^{π Re ν} cancellation (|sin πα|) unless
rotated onto the negative axis, where G_N's singular support sits: the same forward/backward-heat cancellation as everywhere else.

---

## 4. The new structure: Regge contour in ν for the polynomial weight [exact identity + DV]

**Identity.** Let ψ = ψ₂ (degree 14, ψ(0) = 0, (1−p)⁵ | ψ). For 0 < δ < 1 and Y > 0, by Cauchy's theorem for κ = κ^M (analytic in
{Re ν > 0}: B has its poles at ν = ik on Re ν = 0, the c_j-poles are in Im ν < 0, ψ(ν/Λ) is a polynomial):

  ℒ(G_A^ψ) = (2/Λ){ ∫₀^δ ψK dν + 2Re[ ∫_{δ}^{δ+iY} ψκ dν − ∫_{Λ}^{Λ+iY} ψκ dν + ∫_{δ+iY}^{Λ+iY} ψκ dν ] },

and the same with κ_∞ gives G_A^ψ(c). Pieces: *bottom* (real, small ν), *left* (Re ν = δ: β = α − j = 5/2 − y − j + iδ, near-real exponents),
*right* (ν = Λ + iy: weight ψ(1 + iy/Λ) = O((y/Λ)⁵), the Lemma-A endpoint), *top* (Im ν = Y: outgoing decay e^{−β_cY}).

**DV (`ct_check.py`, `t5.py`, `realint.py`; exact ladder, 3000-bit balls; δ = 1/4).** Units: the value ℒ(G^ψ_A).

| state (J, a, b), c, β_c | Y | bottom: \|inf\|, \|fm−inf\|, max\|m−1\| | left | right | top | contour total vs real-axis exact |
|---|---|---|---|---|---|---|
| (2000, 998, 0), 0.445, 1.611 | 12 | 1.6e-4, 2.5e-10, 1.5e-6 | 5.0e-4, 3.5e-10, 6.8e-5 | 1.1e-14, 1.3e-14, 9.4 | 1.6e-8, 1.3e-8, 1.77 | 7.22710203e-5 vs 7.2271020312e-5 |
| (1500, 1500, 4), 0.250, 1.099 | 14 | 8.8e-5, 1.3e-10, 1.5e-6 | 4.7e-4, 8.7e-10, 5.0e-5 | 6.2e-14, 4.3e-14, 1.77 | 4.5e-7, 3.1e-7, 1.88 | 3.46288204e-5 vs 3.46288206e-5 |

(Columns: (2/Λ)∫\|ψκ_∞\|, (2/Λ)∫\|ψ(κ_fm − κ_∞)\|, max\|κ_fm/κ_∞ − 1\| along the piece; value of ℒ(G^ψ) ≈ 7.2e-5 resp. 3.5e-5.)

What this shows:
1. **The real-axis integral (∫\|·\| / value ≈ 3.5·10⁴) becomes a contour integral with ∫\|·\| / value ≈ 5–15.** The cancellation that
   per-frequency methods could not capture is converted into the outgoing decay e^{−β_c Im ν}.
2. **On bottom and left the finite-mass multiplier is 1 + O(10⁻⁶ … 10⁻⁴)** (near-real exponents, Lemma-NB regime).
3. **On the top line the multiplier is O(1)** (max\|m − 1\| ≈ 1.8–1.9; at (20000, 998, 0) and (5000, 3000, 4) the sampled top-line values are
   \|m−1\| ≈ 1.00–1.7, i.e. κ_fm is not larger than κ_∞ there and often much smaller). Its contribution is suppressed by e^{−β_cY}.
4. **The right edge carries the turning point** (\|m\| up to ~10 for a ≈ Λ) but has weight 10⁻¹⁴.
5. The infinite-twist part of each piece is explicit; the identity is exact (contour totals agree with the real-axis ladder integral to
   the quadrature accuracy, 1e-9 and 8e-9 relative).

**Consequence (reduction) [P, given the identity].** For the ideal part,
  \|ℒ(G^ψ_A) − G^ψ_A(c)\| ≤ (2/Λ)[∫₀^δ\|ψ\|\|K−K_∞\| + 2∫_left\|ψ\|\|κ−κ_∞\| + 2∫_right\|ψ\|(\|κ\|+\|κ_∞\|) + 2∫_top\|ψ\|(\|κ\|+\|κ_∞\|)],
so positivity on the sliver follows from per-frequency bounds on κ_fm **on three edges only**, plus the residual r = f − ψ₂/k₀ on the real
axis (mass 2.68e-6, needs only \|K(ν)\| ≤ C per frequency) and the same construction for the dressing and for ψ₃ (k = 3).
The binding analytic input is:

> **Lemma T (open).** For every sliver state and t ∈ [δ, Λ], \|κ_fm(t + iY)\| ≤ M·e^{−β′Y}/√t with explicit M, β′ > 0, for one
> Y = Y(β_c) with (2/Λ)∫_top\|ψ\|·M e^{−β′Y}/√t ≤ 0.1·F_cert(β_c).

DV says M ≈ 3 with β′ = β_c works. A termwise (j-series) proof works only for L ≳ 4 (§1); at small L it needs a resummed (non-termwise)
handle on κ_fm in the upper half-plane.

---
## 5. Theorem E / E′: uniform-in-spin per-power multipliers from the spin recurrence [P + DV]

### 5.1 The recurrence is the dual-Hahn difference equation [P]
S_n(α) := ℒ_{n,N}((1−z)^{−α}) = ₃F₂(α, −n, n+s; −N, N+s; 1) is the dual Hahn polynomial R_{−α}(λ(n); γ = N+s−1, δ = −N, N) of
"degree" −α in λ(n) = n(n+s) (Koekoek–Swarttouw (1.6.5); a formal ₃F₂ contiguity identity, valid for real N because the series
terminates through −n). Its difference equation in n reads

  αS_n = B_nS_{n+1} − (B_n + D_n)S_n + D_nS_{n−1},  B_n = (N+n+s)(n+s)(N−n)/((2n+s)(2n+s+1)),  D_n = −n(N+n+s)(N−n)/((2n+s−1)(2n+s)),

with B_n + D_n = (s−1)(N+n+s)(N−n)/((2n+s−1)(2n+s+1)). **For b = 0 (s = 1), B = −D and S_{n+1} = S_{n−1} + ε_nαS_n,
ε_n = 2(2n+1)/((N−n)(N+n+1))** — the author's §92 recurrence, now with a proof. (`dualhahn_check.py`: 48 cases, s = 1 and 9,
N = 57.3 and 1234.75, complex α: max relative defect 9e-83 at 300 bits.)

### 5.2 Theorem E′ (b = 0) [P]
Let α = σ + iτ with σ ≥ 0, 0 ≤ τ < a := N − n, κ := τ/a, and d_n = 1 − n(n+1)/(N(N+1)). Then, for every integer n ≥ 1 and real N,

  **|d_n^α S_n(α)| ≤ 2c₀ · e^{σ/((a+1)√(1−κ²))} · (2/(1 + √(1−κ²)))^σ · e^{q·arcsin κ} · (1−κ²)^{−q²/2},  q := |α|/τ,**
  c₀ = max(|1−λ₋|, |λ₊−1|)/(2|r|) at n = 0, = 1 + O(|α|/N²) (≤ 1 + 2|α|/a²).
  The growth as κ → 1 is a **power**, (1−κ²)^{−q²/2}, not an exponential: e.g. ≤ 2·31.6·4.7·2^σ at κ = 0.9995.

*Proof.* Write the recurrence as v_{n+1} = M_nv_n, v_n = (S_n, S_{n−1}), M_n = [[2w_n, 1],[1, 0]], w_n = ε_nα/2 = u_n + iy_n. The frozen
eigenvalues are λ± = w ± r, r = √(1+w²), λ₊λ₋ = −1, and with P_n = [[λ₊, λ₋],[1, 1]] the coordinates c_n = P_n^{−1}v_n satisfy
c_{n+1} = (P_{n+1}^{−1}P_n)Λ_nc_n, Λ_n = diag(λ₊, λ₋). Hence, in the max norm,
  |c_{n+1}| ≤ |λ₊(n)|·(1 + (|Δλ₊| + |Δλ₋|)/(2|r_{n+1}|))·|c_n|   (P_{n+1}^{−1}(P_n − P_{n+1}) = (2r_{n+1})^{−1}[[−Δλ₊, −Δλ₋],[Δλ₊, Δλ₋]]),
since |λ₋| = 1/|λ₊| ≤ 1 (Re asinh w ≥ 0 for Re w ≥ 0). Three elementary estimates:
1. log|λ₊| = Re asinh(w) ≤ u/√(1−y²) (integrate (1+t²)^{−1/2} from iy to iy + u; |1 + (iy+x)²| ≥ 1 − y²).
2. σ log(d_k/d_{k+1}) ≥ σ·2(k+1)/P_k ≥ σε_k/2, P_k = (N−k)(N+k+1) (d_{k+1}/d_k = 1 − 2(k+1)/P_k). So the d^σ normalisation cancels the
   leading growth exactly, leaving Σ_k(σε_k/2)((1−y_k²)^{−1/2} − 1).
3. ε_k ≤ 2/m and y_k ≤ τ/m with m = N − k ≥ a + 1; the summands are increasing in (ε, y) and decreasing in m, so
   Σ_k(σε_k/2)((1−y_k²)^{−1/2} − 1) ≤ σ∫_a^∞ m^{−1}(m(m²−τ²)^{−1/2} − 1)dm = σ log(2/(1+√(1−κ²))), and, with |dλ±/dw| ≤ 1 + |w|/|r|,
   |r| ≥ √(1−y²), |w| ≤ qy: Σ_{k≤n−2}(|Δλ₊|+|Δλ₋|)/(2|r_{k+1}|) ≤ Σ q f(y_{k+1})Δy_k, f(y) = (1−y²)^{−1/2} + qy/(1−y²). Since ε_k is convex
   in k (Δy_k ≤ Δy_{k+1}) and f increases, f(y_{k+1})Δy_k ≤ ∫_{y_{k+1}}^{y_{k+2}}f, so the sum is ≤ q∫₀^{y_n}f ≤ q[arcsin κ + (q/2)log(1/(1−κ²))],
   using y_n = τ(2n+1)/(a(2n+a+1)) < κ (the formula for ε at index n, where N − n = a).
The final step gives |S_n| ≤ (|λ₊| + |λ₋|)|λ₊||c_{n−1}|-type factors with |λ₊| + |λ₋| ≤ 2e^{σ/((a+1)√(1−κ²))}; with c₀ = |c_0| the claim follows. ∎
(`thmE.py`: 40 random states, n ∈ [1001, 12000], a ∈ [998, 3600], σ ∈ {¼, 5/4, 5/2, 7/2, 9/2}, τ ≤ 0.95a: max exact/bound = 0.497.)

Sample values of the closed form (σ = 5/2, a = 998): 3.0 (κ = 0.3), 6.2 (κ = 0.6), 32 (κ = 0.9), 246 (κ = 0.99). The bound does not depend on n (spin) or on L.

### 5.3 Numerical comparison [DV] (`energy_dv.py`, `energy_gen.py`, `energy2.py`; exact values from 6000–8000-bit recurrences, checked against
direct ladder sums)
| (n, a), α | exact \|Ŝ\| | Theorem E (energy form) | E′ (eigen-coordinates, stepwise) |
|---|---|---|---|
| (1001, 998), 2.5 + 900i | 1.199 | 5.76 | 4.28 |
| (3000, 998), 3.5 + 900i | 1.903 | 90.9 | 8.39 |
| (10⁵, 998), 2.5 + 900i | 3.37 | 3483 | 19.0 |
| (10⁵, 998), 2.5 + 990i | 9.30 | 4.0e6 | 88.8 |
| (10⁶, 998), 3.5 + 990i | — | — | 200.8 |
| (5000, 3000), 3.5 + 900i | 1.033 | 2.42 | 2.57 |
For comparison TAIL_CERT's Theorem P (valid only for L ≥ L₀) gives log\|Ŝ\| ≤ 10.4 / 123 / 325 at τ ≈ 1000 for L₀ = 20 / 10 / 8
(`tp/kosc_probe.py`, `logs/kosc*.log`). **Theorem E′ is uniform in spin and in L, and is sharp to a factor 2–20 below κ ≈ 0.95.**

**b = 4 (s = 9).** The recurrence is S_{n+1} = (h_n + αe_n)S_n + g_nS_{n−1} with e_n = 1/B_n, g_n = −D_n/B_n = n(2n+s+1)/((n+s)(2n+s−1)) ∈ [0, 1),
h_n = 1 − g_n. The energy form E_n = \|S_n\|² + g_n\|S_{n−1}\|² − 2y_n Im(S_nS̄_{n−1}) satisfies E′_{n+1} = g_nE_n + 4u_n²\|S_n\|² + 4g_nu_nRe(S_nS̄_{n−1})
(u = (h + σe)/2, y = τe/2), and g + √g(1 − g) ≤ 1; the stepwise bound (`energy_gen.py`) gives 6.6–7300 against exact 1.0–2.7 on the same
test grid (the factor ~6.6 is the early-n transient where g_n is small). Closed form for s = 9: not written.

**Limits.** τ < a is essential (turning point): for a ∈ [998, 998.5) the ideal weight reaches beyond the turning point, and for
a ∈ [998, ~1040] the realised residual band ν ∈ [a, 1040] does too. An Airy-zone version is needed there (author's §92 claims
21000e^{G_a(τ)}; not reproduced).

---
## 6. The reduction theorem (what an analytic proof of the sliver now needs) [P]

**Splitting lemma (exact, [P]).** Let κ be any function analytic on the closed rectangle R = [δ, Λ] × [0, Y] (0 < δ < 1) with
K(ν) = 2Re κ(ν) + ρ(ν) for real ν ∈ [δ, Λ]. Then for the ideal weight (ψ(ν/Λ) is a polynomial)

  ∫₀^Λ ψK dν = ∫₀^δ ψK dν + ∫_δ^Λ ψρ dν + 2Re[ ∫_{δ}^{δ+iY} − ∫_{Λ}^{Λ+iY} + ∫_{δ+iY}^{Λ+iY} ] ψκ dν.

Admissible choices: (i) κ = κ^M (§3b; ρ = 0, identity checked [DV]); (ii) κ = κ_J (truncated j-series; ρ = ℒ(R_J), Lemma Loc);
(iii) κ = κ^h := Σ_{m≤M} h^mP_m(ϑ)[(1−z)^{−3}P⁺_ν(β(z))]_{z=c} (FULLPOS3's local h-expansion applied to the *outgoing* function; ρ = FULLPOS3's
far remainder ℛ_M(ν) of the full kernel, exactly, since F = F⁺ + F⁻ and P_m(ϑ) is linear). Choice (iii) is exact by construction and
uses only objects already certified elsewhere (P_m, ℛ_M); the price is that Lemma T must then be proved for κ^h.

**Theorem R (reduction).** For a sliver state, write d³ℒ(G) − M ≥ −(E_ψ₂ + E_ψ₃ + E_r + E_ρ) with M ≥ F_cert(β_c) (FULLPOS3 §4.1(c)), where
* E_ψ₂, E_ψ₃: the four contour pieces of the splitting lemma for ψ₂/F₃ (with the dressing 𝒱_b) and ψ₃/F₄ (prefactor τΛ), bounded by
  \|κ − κ_∞\| on bottom and left (near-real exponents: Lemma NB for Re β ≥ ¼, and the difference recurrence below for Re β < ¼),
  and by \|κ\| + \|κ_∞\| on the top line and the right edge;
* E_r = (2/Λ)d³∫\|r\|\|K\| (realised residual, ∫\|r\| ≤ 2.682e-6, on the real axis; needs a per-frequency bound on K);
* E_ρ = the real-axis remainder of the chosen splitting.
Then **H_ang,b > 0 whenever E_ψ₂ + E_ψ₃ + E_r + E_ρ < F_cert(β_c).** The DV of §4 shows the sizes: bottom/left ≲ 10⁻⁵ of the value, top
≲ 10⁻³ at Y ≈ 25/β_c, right ≲ 10⁻⁹, i.e. a margin of 10²–10⁴ against F_cert at every tested state.

**Lemma V (difference recurrence for near-real exponents) [P].** With Δ(β) := d^βS(β) − 1 and D = N(N+s), the exact three-term recurrence
(TAIL E3) gives [D − β(β−s)]Δ(β−1) = [D + β(s−2β−1)/d]Δ(β) + (β(β+1)/d²)Δ(β+1) + R(β), R(β) = (βc/d²)(βc + 1 + sd), and with
d²D ≥ 3.9a² (TAIL §2.3): |Δ(β−1)| ≤ (1 + (3\|β\|² + (2s+1)\|β\|)/(3.9a²))|Δ(β)| + (\|β\|² + \|β\|)/(3.9a²)·|Δ(β+1)| + \|β\|(\|β\| + 1 + s)/(3.9a²)·(1 + 10⁻⁶),
for \|β\|² ≤ 0.1·3.9a². Started from Lemma NB at Re β ∈ {¼, 5/4} (\|Δ\| ≤ 10⁻⁶), this gives at a = 998, s = 9 (worst case):
\|Δ\| ≤ 3.2e-5, 1.8e-4, 1.1e-3, 3.3e-3, 7.3e-3 at Re β = −4.75, −9.75, −19.75, −29.75, −39.75, uniformly in spin (scalar evaluation in §6
of the log); the left-line weight beyond Re β ≈ −10 is ≤ e^{−10β_c} of the total.
(DV, §4: max\|m − 1\| on the left line ≤ 7.8e-4 up to Y = 24, i.e. Re β ≈ −21.)

**Remaining analytic inputs, precisely:**
1. **Lemma T (top line), small L.** A per-frequency bound on the finite-mass outgoing kernel at Im ν = Y with constant ≪ e^{β_cY}. For
   L ≳ 4 it follows from Theorem E′ (σ = ¼, 5/4) + the β-recurrence + termwise j-sums. **For L ≲ 4 every termwise route fails by the
   hump e^{\|ν\|d/2} of the (1−z)-connection series**, so a resummed (full-kernel) estimate is required.
2. **Turning point.** For a ∈ [998, ~1160] the frequencies ν ∈ [0.95a, 1100] (ideal weight near Λ; realised residual up to 1100) are at or
   beyond the turning point ν = a, where Theorem E′ does not apply. The weights there are small (ψ ~ (1−p)⁵; residual mass ≤ 4e-10 in
   [1000, 1040], ≤ 1e-22 beyond) but need a bound of the turning-point type e^{G_a(ν)}, G_a ≈ (2/3)√(2/a)(ν−a)^{3/2} (author §92, not
   reproduced). Theorem P at moderate L gives e^{100…325} there, which is useless.
3. **E_ρ.** [Resolved for L ≳ 4 by Lemma Loc′, §7.3.] Lemma Loc (choice (ii)) works only for L ≳ 7 (K ≈ e^{2640} from \|F_ν\| ≤ (π/2)e^{πν/2} on the disc); choice (i) has no
   remainder but its κ^M has no bounds yet; choice (iii) has the certified FULLPOS3 far remainder (fails at a ≈ 998 for β ≈ 2 in FULLPOS3's
   implementation, e^{−74} claimed by the author's §63, not reproduced).
4. **b = 4 version of Theorem E′** (closed form): only the stepwise energy bound exists (§5.3).

---
## 5.4 b = 4 by the helicity identity [P + DV]
S^{(9)}_{n,N}(β) = Σ_{i=0}^{4} C(4,i)(−1)^i S^{(1)}_{n+4,N+4}(β − i) / U^{(1)}_4(n+4, N+4) (same twist a) — the §72 identity at the level of a
single power (`dualhahn_check.py`-type check: relative defect 1e-80…1e-86 at three states). With d = 1 − n(n+9)/(N(N+9)) and
d₀ the shifted s = 1 value, d/d₀ = (D+20)/D, U₄ ≥ (1−d₀)⁴(1 − 10⁻⁶), so for σ ≥ 0
  |Ŝ^{(9)}(β)| ≤ (1 + 10⁻⁵)[|Ŝ^{(1)}(β)| + Σ_{i≥1} C(4,i) d₀^i 2^i max(|Ŝ^{(1)}(σ)|, |Ŝ^{(1)}(σ+1)|)] / ((1−d₀)⁴(1 − 10⁻⁶))
(downward β-recurrence growth ≤ 2 per step). For L ≥ 7 the i ≥ 1 terms are ≤ 10⁻²; at small L this loses ((1+d)/c)⁴ (the §72 corner).

## 5.5 Lemma X: crossing the turning point (b = 0) [P; stated with explicit constants in `crossing.py` / `rt_pp.cross`]
For α = σ + iτ, σ ≥ 0, 900 ≤ τ ≤ 1150 and every state with a ≥ 998 (uniformly): split the spin steps by y_k = ε_kτ/2 ∈ [0, 1 − δ]
(phase I: eigen-coordinates as in E′), (1 − δ, 1 + δ) (window: T_k = i^{−k}S_k obeys T_{k+1} = 2μ_kT_k − T_{k−1}, μ = y − iu; the norm
|T| + |D|/√η, D = T_k − T_{k−1}, grows by ≤ 1 + √η + η per step, η = 2(δ + u); entry/exit conversions use |λ± − i| = O(√δ)), and
[1 + δ, τ/(a+1)] (phase III: eigen-coordinates with |λ₊| = e^{acosh y}; Σ acosh y_k ≤ G_a(τ) = τ arccos(a/τ) − a·acosh(τ/a)),
δ = cτ^{−2/3}. Result: **|Ŝ(σ+iτ)| ≤ e^{G_998(τ) + C(σ, τ)}** with C ≈ 8.4 (σ = ¼) … 10.4 (σ = 5/2) at τ ≈ 1000 (`cache/crossing_dv.txt`:
holds at all 70 test points against exact multipliers, e.g. (n, a) = (20000, 998), τ = 1080: exact log\|Ŝ\| = 12.6, G = 21.7, bound 32.9).
Two exact inequalities make the phase-I part sharp: **ε_{k+1} − ε_k ≥ ε_k²/2** (so Δy_k ≥ y_k²/τ, turning Σ(σε_k/2)((1−y_k²)^{−1/2} − 1) into
σ∫((1−y²)^{−1/2} − 1)dy/y), and the convexity of ε_k in k (adiabatic sums by left Riemann sums).

---

## 7. Route T′: certified closure of L ∈ [8, 16], extended to [5, 8] in §7.3–7.4 (all a ≥ 998, all spins, b ∈ {0, 4}) [P + certified arithmetic]

**Statement (Theorem T′).** For every state with a ≥ 998, b ∈ {0, 4}, n = J − b ≥ 1001 and L = −log d ∈ [L₁, 16]:
  d³ℒ(G) ≥ M − (E₂ + E₃ + E_r + E_ρ) > 0,  M ≥ F_cert(β_c),
with the four error terms of Theorem R (§6) bounded as follows (`rt_cert.py`, `rt_pp.py`, `rt_rmass2.py`, `rt_rem.py`):
* **E₂, E₃ (Regge contour of ψ₂ with the dressing, and of ψ₃):** contour 0 → ¼ → ¼ + iY → Λ + iY → Λ with Re α = ¼ − 2 on the top line
  (Y = 4.25 for ψ₂, 5.25 for ψ₃), truncated j-series κ_{J₀} − κ_∞ bounded termwise (|c_j| in Arb on cells, 59 explicit terms plus a
  rigorous geometric tail with ratio sup_{j≥59}|c_{j+1}/c_j|·2d < ½), |B| on cells by |B(mid)|e^{r·sup|(log B)′|} (digamma via the
  shifted Stirling bound), |ψ| by exact Taylor shifts; multipliers: Lemma NB (re-derived with explicit d: g[(σ)₃/3 + (1+s)(σ)₂d/2] +
  (4/15)g²(σ)₅, g = h/d² ≤ ¼(1+s/N)/(1−d/2)²/a²) and Lemma V on bottom/left; Theorem E′ (+ helicity for b = 4) at σ ∈ {¼, 5/4} for
  τ ≤ 0.95·998, Lemma X above, then the d-normalised β-recurrence [D − β(β−s)]Ŝ(β−1) = [D + β(s−2β−1)/d]Ŝ(β) + (β(β+1)/d²)Ŝ(β+1)
  (d²D ≥ 3.9a² for N ≥ 40a) on top line and right edge.
* **E_r (realised residual on the real axis):** |2Re κ_{J₀}(ν)| per frequency (same per-power bounds, σ = 5/2, 7/2 for the dressing),
  against certified residual masses: EPSW ranges for ν < 898, and for ν ≥ 898 **EPSW's own panel bounds regrouped into unit cells**
  (`rt_rmass2.py`, total 3.3e-9 on [898, ∞)); [0, 0.05] by the maximum principle on |ν| = 0.05 (g = κ(ν) + κ(−ν) is analytic in |ν| < ½).
* **E_ρ (j-truncation remainder):** TAIL_CERT's Lemma Loc with J₀ ∈ {200, 400, 600, 750} (constants lI₂ = 2777.6, 3315.5, 3804.5, 4157.3;
  `cache/rt_rem.json`); this is what limits L₁ (≈ 7.9: at L = 8, E_ρ/F ≈ 1.5e-26 with J₀ = 750; at L = 7.5 it is e^{+300}).
* **M:** F_cert(β_c) from LRINF_TIGHT's certified cells (min over the β-range of each L-cell).

**Final scan, L ∈ [8, 16], cells of width ½ (`cache/rt_cert_final_8_16.txt`; errors at the lower end of each cell, F_cert minimised
over the cell; all PASS).**
| L-cell | F_cert ≥ | ψ₂ bottom / left / top / right | ψ₃ (all) | residual (worst piece) | E_ρ (J₀ = 750) | **ρ** |
|---|---|---|---|---|---|---|
| [8, 8.5] | 1.11e-7 | 1.7e-5 / 1.2e-5 / 5e-12 / 5e-9 | ≤ 5e-10 | 2.6e-2 (cell [1031, 1032]: 1.4e-2) | 1.5e-26 | **0.026** |
| [9, 9.5] | 6.22e-8 | 1.9e-5 / 1.3e-5 / ... | ≤ 1e-10 | 2.6e-2 | < 1e-300 | **0.026** |
| [10.5, 11] | 2.30e-8 | 2.3e-5 / 1.4e-5 | ≤ 4e-11 | 3.2e-2 | 0 | **0.032** |
| [12, 12.5] | 9.75e-9 | 2.6e-5 / 1.5e-5 | ≤ 1e-11 | 3.6e-2 | 0 | **0.036** |
| [14, 14.5] | 3.14e-9 | 3.0e-5 / 1.5e-5 | ≤ 1e-12 | 4.2e-2 | 0 | **0.042** |
| [15.5, 16] | 1.35e-9 | 3.2e-5 / 1.6e-5 | ≤ 1e-13 | 4.6e-2 | 0 | **0.046** |
(max over all 16 cells: ρ ≤ 0.0464.) The contour pieces are negligible (≤ 5e-5); the realised residual near the turning point
(ν ≈ 1031, a ≈ 998, Lemma X) and at ν < 0.4 dominates. **The floor L₁ = 8 is set by Lemma Loc** (E_ρ/F ≈ e^{5909 − 748L + 16}):
at L = 7.5 it is e^{+300}, although every other piece is ≤ 0.035 there (§7, cell [7, 7.5]).

**Rigor notes.** Arb for every special-function and polynomial bound; positive upper bounds are then accumulated in floats with a
final 1 + 10⁻⁹ factor (≤ 10⁶ operations). Inputs taken as certified from elsewhere: LRINF_TIGHT F_cert cells, EPSW_CERT panel machinery
and range totals, CONTOUR's Gaussian tail formula (ν ≥ 1150, contributes < e^{−700}), TAIL_CERT's Lemma Loc and K-functions, FULLPOS3
§4.1 (main term ≥ F_cert), TAIL_CERT E1–E5. New and unreviewed: Theorem E′, Lemma X, the helicity bound §5.4, the splitting lemma §6,
Lemma NB with explicit d, Lemma V.

### 7.3 Lemma Loc′: a localised remainder bound (replaces TAIL's Lemma Loc below L = 8) [P, new, unreviewed + certified arithmetic]
TAIL's Lemma Loc bounds R_J by one global constant, K\|1−ξ\|^m with K ≈ e^{1600–4200}, over the whole disc. That K is what set the floor
L₁ ≈ 7.9. Lemma Loc′ keeps the Cauchy/Beta representation but splits the integration domain, so that each piece of R_J is
bounded only where it lives (`rt_loc3.py`).

**Setting.** R = R_J(ξ;ν) = F_{p,ν}(ξ) − T_J, T_J = Σ_± B(±ν)Σ_{j<J}c_j(±ν)(1−ξ)^{j−e₀∓iν}, e₀ = p − ½. On \|1−ξ\| < 1 it is also the tail
Σ_± B Σ_{j≥J}. The representation is ℒ(R) = E_B[(2πi)⁻¹∮_{\|z\|=ρ}(1−z)^{−a−1}R(Bz)z^{−n−1}dz]/C(N,a), with B ~ Beta(n+s, a) and ρ = n/N.

**Lemma Loc′.** Let W(θ) = \|1−ρe^{iθ}\|^{−a−1}, K₀ = ρ^{−n}(2π)⁻¹∫W/C(N,a), and ξ = Bρe^{iθ}. Then
  \|𝒱_b ℒ(R)\| ≤ 1.6·41·[ sup_{G₁}\|R\| + f(θ₀)sup_{G₂}\|R\| + Σ_k (1.02d/b_k)^{998} sup_{T_k}\|R\| + f(θ₁)sup_{Th}\|R\| ],
where the regions and their bounds are:
* G₁ = {B ≥ 1−b₁, \|θ\| ≤ θ₀} and G₂ = {B ≥ 1−b₁, θ₀ < \|θ\| ≤ θ₁}: \|1−ξ\| ≤ 0.51d + θ + b₁ + δ (tail series);
* T_k = {1−b_{k+1} ≤ B < 1−b_k, \|θ\| ≤ θ₁}: b_k − δ ≤ \|1−ξ\| ≤ 0.51d + θ₁ + b_{k+1} + δ. The b_k are geometric (ratio 1.6) up to 0.95, and the
  last tier is B < 0.05 with \|1−ξ\| ≤ 2 + δ. On the tiers \|R\| ≤ \|F\| + \|T_J\|;
* Th = {\|θ\| > θ₁}: 0.49d − δ ≤ \|1−ξ\| ≤ 2 + δ;
* f(θ) = π√a·e^{0.51}(1−ρ)⁻¹((1−ρ)/q(θ))^{a+1} with q² = (1−ρ)² + 4ρsin²(θ/2), and δ = 0.255d.

*Proof.*
1. The θ-weights f(θ) bound the W-mass of \|θ\| > θ (from I_W ≥ e^{−0.51}(1−ρ)^{−a}/(π√a)).
2. P(B < 1−b) ≤ E(1−B)^{998}/b^{998} = (a)_{998}/((N+s)_{998}b^{998}) ≤ (2a/(Nb))^{998} ≤ (1.02d/b)^{998}.
3. K₀ ≤ 1.6 (Stirling/Robbins for C(N,a), Gautschi).
4. The dressing is handled as in TAIL_CERT §4 (Cauchy on circles of radius ½(1−ρ) ≤ δ, constant 40, used with m = 0); this enlarges every
   region by δ and gives the factor 41.
5. On each region: \|F_{p,ν}\| ≤ (π/2)e^{π\|ν\|/2}\|1−ξ\|^{−p} (FULLPOS3 Lemma B), \|(1−ξ)^{∓iν}\| ≤ e^{π\|ν\|/2} since \|arg(1−ξ)\| < π/2, and
   \|c_j\| ≤ ∏_{i<j}((i+½)² + ν_hi²)/(√((i+1)² + 4ν_lo²)(i+1)), all in Arb.
6. The tail ratio is sup_{j≥J}\|c_{j+1}/c_j\| ≤ (1 + ν/(J+1))².
7. **Small ν ∈ [0, 0.4]** (pole of B at 0). Write B = E/(2iν) with E = Γ(1+2iν)/Γ(½+iν)². Then
   Σ_±B(±ν)Φ(±ν) = (2iν)⁻¹∫_{−ν}^{ν}(EΦ)′, which is ≤ (1/r′)sup\|EΦ\| on {\|Re w\| ≤ 0.45, \|Im w\| ≤ r′ = 0.05}. Here sup\|E\| is in Arb,
   and the c_j bound for complex w is \|a₀+i\|² ≤ (i+½+\|y\|)² + x² and \|b₀+i\| ≥ i+1−2\|y\|. ∎

State relations used: 0.49d ≤ a/N ≤ 0.51d and Td² ≥ 3.9a². Exact algebra gives a/N ≈ (d/2)(1 + d/4); the worst sampled values over all states with
L ∈ [4.45, 8] are a/N ≤ 0.5015d and Td²/a² ≥ 3.977.
[DV] The exact values at a ∈ {998, 1500, 3000} and L ∈ {4.5, 6, 8} are:
* K₀ = 1.0003 (bound 1.6);
* the θ-fraction at θ = 0.02, a = 998, L = 4.5: e^{−1314} exactly vs bound e^{−1301}.

**Numbers (J = 1500, a = 998 worst).**
* log\|𝒱_bℒ(R_J)\| ≤ −2753, −1881, −898, −873 at ν = 1, 500, 1029, 1040 (L = 5).
* At L = 4.5: −2504 … −382.
* At L = 3.5 it is +597 at ν = 1029, so **Lemma Loc′ itself holds down to L ≈ 4**.
* ν-integrated against \|f\| (sup\|f\| ≤ 21.34 on [0, 1030]; CONTOUR's tail(T) beyond, with tail(1300) = e^{−45450}):
  log∫\|f\|\|𝒱_bℒ(R)\| = −403 on [4.5, 5], −894 on [5, 5.5], −2100 on [7.5, 8]. The k = 3 analogue is −491 … −2186.
* **J₀ = 1500 is admissible:** the j-sum growth bound (M′_j ×2 per step, REC) needs (A+B)/den ≤ 2 at
  \|β\| = hypot(J₀+3, 1151) + 2. Its value is 1.946 at L = 4.5; it is asserted in `rt_cert.remainder_terms_loc3`.

### 7.4 Route T′ with Lemma Loc′ on [4.5, 8] (`python3 rt_cert.py L1 L2 0.5 2 1500`; logs/rtL3_*.log, cache/rt_cert_*_loc3J1500.txt)
| L-cell | F_cert ≥ | ψ₂ bottom / left / top / right | residual (cell ν ∈ [1031,1032]) | E_ρ (Loc′) | **ρ** |
|---|---|---|---|---|---|
| [4.5, 5] | 8.56e-7 | 1.4e-5 / 1.3e-5 / 2.1e-4 / 2.0e-5 | 15.9 (16) | 1e-177 | **15.9 FAIL** |
| [5, 5.5] | 6.40e-7 | 1.4e-5 / 1.3e-5 / 6.5e-6 / 9.9e-7 | 0.947 (0.94) | < e^{−745} | **0.947 PASS** |
| [5.5, 6] | 4.78e-7 | 1.4e-5 / 1.3e-5 / 4.0e-7 / 1.5e-7 | 0.179 (0.17) | < e^{−745} | **0.179 PASS** |
| [6, 6.5] | 3.57e-7 | 1.5e-5 / 1.2e-5 / 3.5e-8 / 4.2e-8 | 0.070 (0.060) | < e^{−745} | **0.070 PASS** |
| [6.5, 7] | 2.66e-7 | 1.5e-5 / 1.2e-5 / 3.5e-9 / 1.8e-8 | 0.042 (0.031) | < e^{−745} | **0.042 PASS** |
| [7, 7.5] | 1.99e-7 | 1.6e-5 / 1.2e-5 / 3.9e-10 / 1.0e-8 | 0.032 (0.021) | < e^{−745} | **0.032 PASS** |
| [7.5, 8] | 1.49e-7 | 1.6e-5 / 1.2e-5 / 4.5e-11 / 6.7e-9 | 0.028 (0.016) | < e^{−745} | **0.028 PASS** |

(ψ₃ pieces are ≤ 6.4e-8 everywhere. "< e^{−745}" means the float exp underflowed: the certified log-bounds are −894 and below, see
cache/rt_loc3.json.)

**What now limits L₁ = 5: the realised residual at the turning point**, not the remainder. The binding piece is the unit cell
ν ∈ [1031, 1032] (a = 998, Lemma X). Its size is roughly e^{G_a(ν)} × the termwise hump Σ\|c_j\|d^j\|Ŝ_j\| ~ e^{1.3νd/2}. It is 0.94 at L = 5 and
16 at L = 4.5. Refining ν-cells does not help, because the residual masses are certified per unit cell (EPSW panels).

## 8. TAIL route with Theorem E′ (TAIL+E′): the tail now closes at L ≥ 15.5 [certified arithmetic on TAIL_CERT's lemmas]
`tp/tcE_common.py` replaces Theorem P by min(Theorem P, Theorem E′(+helicity)) in copies of TAIL_CERT's assembly and k = 3 scripts
(`tp/tcE_assemble.py`, `tp/tcE_k3.py`, `tp/tcE_remainder.py`; tailcert/ untouched, its Theorem-P table rebuilt in analytic/cache).
| L₀ | k = 2 ρ (leg / vertical / short) | k = 3 | remainders | total ρ |
|---|---|---|---|---|
| 16 | 0.163 (0.16 / 6.6e-3 / 2.6e-4) | 7.4e-3 (was 0.673·e² at 16) | 1e-92.5, 1e-309 | **≤ 0.171** |
| 15.5 | 0.486 (0.47 / 1.8e-2 / 3.7e-4) | 1.35e-2 | 1e-49.6, 1e-266 | **≤ 0.50** |
The k = 3 integral ∫|ψ₃||B|X₃ drops from 504 to 0.83 (the multiplier on the real axis was the whole looseness, as TAIL_AUDIT
suspected). Below 15.5 the realised-f leg (∫|fB| = 5.3e11) binds; Route T′ takes over there.

---

## 9. Final summary (end of session)

**Theorem (sliver, L ≥ 5).**
* Hypotheses: Δ_g = 1000, a ≥ 998, J ≥ 1001 (n = J − b), b ∈ {0, 4}, and L = −log d ≥ 5.
* Conclusion: d³ℒ(G) ≥ M − (E₂ + E₃ + E_r + E_ρ) with Σ ≤ ρ·F_cert(β_c) ≤ ρ·M, ρ < 1, hence **H_ang,b > 0**.
* Covered region: the whole sliver with L ≥ 5, uniformly in a, which is more than a < a*(L).

| L | method | decisive constant |
|---|---|---|
| [5, 8] | Route T′ + Lemma Loc′ (J₀ = 1500) | ρ ≤ 0.947 on [5, 5.5]; ≤ 0.179, 0.070, 0.042, 0.032, 0.028 on the next cells |
| [8, 16] | Route T′ + TAIL's Lemma Loc (J₀ ≤ 750) | ρ ≤ 0.0464 |
| ≥ 15.5 | TAIL+E′ | ρ ≤ 0.50 (15.5), 0.171 (16) |
| ≥ 20 | TAIL_CERT | ρ ≤ 0.066 |

**Grades.** Everything above is [P] + [certified] arithmetic (Arb upper bounds, float accumulation × (1 + 10⁻⁹)), modulo the items below.
* **New and unreviewed:**
  - Theorem E′ (§5.2);
  - Lemma X (§5.5);
  - the helicity bound (§5.4);
  - the splitting lemma and Theorem R (§6);
  - Lemma NB with explicit d;
  - Lemma V;
  - **Lemma Loc′ (§7.3)**;
  - the state relations at L ≥ 4.5, which are exact algebra, only spot-checked on a grid.
* **Taken as certified from elsewhere:**
  - LRINF_TIGHT F_cert cells;
  - FULLPOS3 §4.1(c) (main term ≥ F_cert) and Lemma B (\|F_{p,ν}\| bound);
  - EPSW_CERT panel bounds and range masses;
  - CONTOUR's tail formula;
  - TAIL_CERT E1–E5, its Lemma Loc and its dressing constant 40;
  - sup\|f\| ≤ 21.34 on [0, 1300].

**Obstructions below L = 5 (exact).**
1. **Turning-point residual.** The realised-residual unit cell ν ∈ [1031, 1032] at a = 998 contributes 0.94·F_cert at L = 5 and
   16·F_cert at L = 4.5. Its bound is Lemma X's e^{G_a(ν)} times the termwise j-sum Σ\|c_j\|d^{j}\|Ŝ_j\| ≈ e^{1.3νd/2}. This is the
   (1−z)-connection hump at ν ≈ a; the residual masses are only certified per unit cell. Possible cures:
   - a resummed (non-termwise) bound of K_J near the turning point;
   - finer certified residual masses near ν ≈ 1031 together with a pointwise K bound that decays across the cell.
2. **Lemma Loc′** itself fails near L ≈ 4: log bound +597 at L = 3.5, ν = 1029; the tier/Th terms grow like e^{πν/2}Σ\|c_j\|.
3. **Top line (Lemma T).** It is 2.1e-4·F at L = 4.5, and from the §1 analysis termwise top-line bounds die near L ≈ 2–4 (hump e^{\|ν\|d/2}).
   A full-kernel outgoing bound in the upper half plane is needed for small L.
4. **j-sum growth.** J₀ ≤ ~1500 is required by the growth constant (A+B)/den ≤ 2 at \|β\| ≈ 1900. It is not binding.

**Files (all new).**
* `ads4/code/analytic/`:
  - Route T′: `rt_cert.py` (certificate; 5th argument J > 0 selects Lemma Loc′), `rt_pp.py` (per-power bounds), `rt_loc3.py` (Lemma Loc′;
    `rt_loc2.py` is the superseded probe), `rt_rmass2.py`, `rt_rem.py`;
  - Lemma X: `crossing.py`;
  - Theorem E/E′ DV: `thmE.py`, `energy*.py`;
  - TAIL+E′: `tp/`;
  - contour DV: `ct_check.py` and `t*.py`.
* Outputs:
  - logs in `logs/rtF_*.log` and `logs/rtL3_*.log`;
  - certificates in `cache/rt_cert_final_8_16.txt` and `cache/rt_cert_*_loc3J1500.txt`;
  - `cache/rt_loc3.json`;
  - `cache/` is git-ignored except its `.gitignore`.
