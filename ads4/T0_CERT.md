# T0_CERT: the leading-Regge (hybrid) action of Ω^EE is ≥ 0.199·N_s on every heavy state, for every Δ_g ≥ 300 [certified]

*25 Nov 2026. Unreviewed, not committed. No existing file was edited. Target T0 of `UNCOND_INVENTORY.md` (items M1, M20).*
* *New code:* `ads4/code/t0cert/`: `t0_common.py`, `t1_curv.py` (Parts K, T), `t2_lrinf.py` (Part L), `t3_flat.py` (Part F),
  `t0_checks.py` ([DV] sanity checks), `run_all.sh`. Logs and JSON in `t0cert/cache/`, whose `.gitignore` ignores everything but itself.
* *Read-only imports:* `lrtight/lrt_common.py`, `gaplarge/gl_common.py` (and through it `lrtight/lrt2_m1.py` and the caches
  `lrtight/cache/lrt1_constants.json`, `lrt2_m1.json`), `robust/{nf_cert,rcert}.py` (and `flatbound/common.py`), `largespin/ls_core.py`
  (checks only). No file of theirs was written and no bytecode was written.
* *Compute:* ≈ 0.6 CPU-h in total (including exploration), ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as in UNIFIED_THEOREM_v3 §0. **[certified]** means every inequality is either proved here (or in a cited [P]/[certified]
document) or decided in python-flint/Arb ball arithmetic or exact rational arithmetic. Continua (b̂, u, J, Λ) are covered by
cells with proved between-grid bounds; nothing is sampled.

---

## 0. Verdict

### 0.1 Theorem T0 (as proved)

Let N_s and H be the hybrid normalisation and action of a heavy state (§1.1). For the heavy-state kernels
* `h = ++`: every even J ≥ 0, with Legendre kernel P_J;
* `h = +−`: every J ≥ 4, even and odd, with Jacobi kernel P^{(0,8)}_{J−4};

every mass m ≥ Λ (so u = Λ²/m² ∈ (0, 1]; this contains every state at or above the twist edge, and also the unphysical
strip Λ ≤ m < m_gap), every impact parameter b̂ = (2J+1)Λ/m ∈ (0, ∞), and every gap Δ_g ≥ 300 (Λ = √((Δ_g−2)(Δ_g−1)), real Δ_g):

    H_h(J, m; Δ_g)  ≥  c · N_s,h(J, m; Δ_g),        c = 199/1000.                         [certified]

At **Δ_g ≥ 1000**, the same holds with **c = 1999/10000** [certified].

* **c is essentially sharp.** As Λ → ∞ the curvature vanishes and H/N_s → 𝔉/N, whose infimum is 0.2000000103 (NF02 at `+−`,
  J = 5, m = M; CERT_CLOSE §1.1). So no c > 0.2000000103 can hold at all gaps. A float scan (`t0_checks.py` [DV]) gives
  min H/N_s = 0.19999969 at Δ_g = 300 (`+−`, J = 4, m = Λ).
* **Matrix form.** The hybrid off-diagonals vanish (POS_ERRORS §2.1). So for even J the 2×2 hybrid matrix
  diag(H₊₊, H₊₋), normalised by N_s, has λ_min ≥ c for every coupling vector.
* **Parities and the OO block.** Odd J has only the parity-odd `+−` structure, whose hybrid kernel is the `+−` kernel above.
  The OO block of even J is covered through the identity "OO hybrid = κ × `++` hybrid, κ > 0" (order-0 identity of
  PARITY_ODD_CLOSE (C)). That identity is **[P given helicity lemma (A) [D]]**, so for OO the grade of T0 is
  [certified given (A)]. Every other sector is [certified] with no premise.
* **No physics input.** A1, A2 and the CFT data are not used: T0 is a statement about blocks and the functional's weights.

### 0.2 What "the leading-Regge action of Ω^EE" means here (the M20 question)
T0 is a statement about the **de-shifted (target-weight) leading-Regge action**: the k = 2 weight enters only through
ψ₂(ν/Λ). That object is defined at every Δ_g, because ψ_k is Δ_g-independent by definition (GAP_SCALING §1.1). By the de-shift
identity it equals the leading-Regge action of Ω^EE whenever the realised k = 2 weight X_Λ satisfies that identity exactly. Exact
existence of X_Λ is [NE + D-sketch] (DESHIFT_CLOSE §2.6).

**M20 is not closed.** A continuum X_Λ with a Λ-uniform ε_W was **not** built. For the **realised** functional (published
continuum weights) T0 holds as a corollary at **Δ_g = 1000 and Δ_g = 300** only (§4.2). It uses ε_W [DV] and needs a cap
β ≲ 200 on the impact parameter (J + ½ ≲ 10⁴³·m). Both caveats come from the weight, not from the hybrid. So, as the brief
allowed: **every-gap T0 for the de-shifted leading-Regge action, and Δ_g ∈ {300, 1000} for the realised weight.**

### 0.3 Coverage and how the inventory's M1 sub-items are disposed of
The proof does **not** go through Lemma C′, ROUTEA-U or the m/Λ bands. It splits only on b̂ (Λ-uniformly), so the four M1
sub-items disappear:

| inventory item | status here |
|---|---|
| hole m < 1.3Λ, b̂ ≲ 114 | covered: Part F for b̂ ≤ 16 (every J, u); Part K for 16 ≤ b̂ ≤ 2Λ (every u ≤ 1 at once) |
| Lemma C′ on μ-intervals (sampled masses) | not needed: the curvature bound depends on (b̂, u, σ, Λ) only; u ∈ [0, 1] is covered by affinity in u |
| m/Λ > 20 | covered (u → 0 is the endpoint u = 0 of the affine check; Part F's region S covers J ≥ 100 at every u) |
| b̂ > 2Λ (hybrid vs F_LR) | covered: Part L (LR∞ with a Λ_min anchor, plus a new flat-Regge remainder lemma) |
| ROUTEA's [D] derivation | not used; replaced by Lemma 2 [P] (§2.2) |

**What remains.** Nothing inside T0's statement. Outside it:
1. **M20.** The realised weight at gaps other than 300 and 1000, and a certified ε_W (now [DV]) — §4.
2. **The OO identification** (helicity lemma (A) [D]).
3. **CERT_CLOSE's corner code is marked "unreviewed" there.** T0 inherits that flat certificate (𝔉 ≥ N/5) as an input.
4. Everything **beyond leading Regge** (δ), which is the actual content of (Pos), is untouched (UNCOND_INVENTORY §1.3).

---

## 1. Setting

### 1.1 The object
For a heavy state (J, m), helicity kernel h, σ₊₊ = +1, σ₊₋ = −1, u = Λ²/m², b̂ = (2J+1)√u and sinh(β/2) = (J+½)/m
(so b̂ = 2Λ sinh(β/2)). With ψ_k the NF02 weights (`prl/data/nf02_psi.json`), the definitions are those of
`spinfunc/ee1000_common.hyb` / `largespin/ls_core.hyb_fast`:

    B₂ = ∫₀¹(2 − up²)ψ₂ D_h dp,   B₃ = σu∫₀¹ψ₃ D_h dp,   D_h = D_h(1 − 2up²)   (P_J, or P^{(0,8)}_{J−4} for +−),
    cv₂ = 2C₂,  cv₃ = uσC₃,   C_k(b̂, Λ) = ∫₀¹ψ_k(p)[P_{−½+iΛp}(cosh β) − J₀(b̂p)]dp,
    H = B₂ + B₃ + cv₂ + cv₃,    N_s = |B₂ + cv₂| + |B₃ + cv₃|.

The flat part (B₂, B₃) depends only on (h, J, u). The curvature depends only on (b̂, Λ, u, σ) and is affine in u.

### 1.2 Lemma 0 (sign patterns) [P]
`H − cN_s = min over (s₂, s₃) ∈ {±1}² of [FF_g(c) + (1 − cs₂)cv₂ + (1 − cs₃)cv₃]`, where
`FF_g(c) := (1 − cs₂)B₂ + (1 − cs₃)B₃` is exactly the robust-margin action of ROBUST/NOFORWARD (`rcert.gens` with δ = c,
a = 0). *Proof:* |x| + |y| = max over signs of (s₂x + s₃y). ∎

### 1.3 Lemma 1 (flat) [certified, CERT_CLOSE §1]
𝔉 = B₂ + B₃ ≥ N/5 on every flat state, with N = |B₂| + |B₃|. Hence for c ≤ 1/5 and every g:
FF_g(c) ≥ 𝔉 − cN ≥ (1/5 − c)N ≥ 0.

So T0 follows from:
* **(K/L)** on b̂ ≥ b_pos: the curvature part of every pattern is ≥ 0;
* **(T/F)** on b̂ ≤ b_pos: FF_g(c) exceeds the largest possible negative curvature part.

Here b_pos = 16.

---

## 2. The curvature in gap units

### 2.1 Lemma 2 (gap-unit form of the curvature) [P]
For b̂ ≤ 2Λ and ε := Λ⁻²:

    Λ² C_k(b̂, Λ) = a₁ₖ M₁(β) + (1/24) ∫₀¹ χ_k''(p) J₀(b̂p) dp + ε (2/π)∫₀^{π/2} h_k(b̂ sin φ) dφ,

where
* a₁ₖ = ψ_k′(0) (188.896 / 50.643 for k = 2 / 3);
* χ = 2ψ − pψ′, an exact polynomial; χ'' has the factor (1−p)², and χ″(0) = 0;
* |h_k(x)| ≤ h̄_k(x) := (3x⁵/640 + x⁵/192)|r′(x)| + (3x⁴/128)|r(x)| + (x⁶/1152)·sup_{[5x/6, x]}|r″|, with r = Φ + a₁/x² and
  Φ(x) = ∫₀¹ψ cos(px) dp;
* M₁ is exactly LRINF_TIGHT's M₁ (the function with certified cell lower bounds in `lrt2_m1`).

*Proof.*
1. **Mehler–Dirichlet** in Lemma C's variables (POS_ERRORS §2.1, steps 1–3, with s = B sin φ) gives
   C_k = (2/π)∫₀^{π/2}[wΦ(τ) − Φ(x)]dφ. Here x = b̂ sin φ, τ = 2Λ asinh(x/2Λ) and w = (1 + x²/4Λ²)^{−½}; Fubini is absolute
   since |P_{−½+iν}| ≤ 1 and |J₀| ≤ 1.
2. **Split** Φ = r − a₁/x². The a₁ part is (2/π)∫a₁[1/x² − w/τ²]dφ = a₁M₁(β)/Λ², with
   M₁(β) = (2/π)∫₀^{π/2}[σ⁻² − (t² cosh(t/2))⁻¹]dφ, where σ = 2 sinh(β/2) sin φ and 2 sinh(t/2) = σ.
   * The change of variables φ → t sends (2/π)dφ/cosh(t/2) to (√2/π)K dt (Mehler) and (2/π)dφ/σ² to
     (√2/π)K·cosh(t/2)/(4 sinh²(t/2)) dt.
   * Hence M₁ = (√2/π)∫₀^βK(t)c(t)dt with c = cosh(t/2)/(4 sinh²(t/2)) − 1/t². This is LRINF_TIGHT §1.1's definition, and
     the integral is absolutely convergent.
3. **Taylor expansion with alternating-series bounds.** On y = x/2Λ ∈ [0, 1]:
   * 0 ≤ asinh y − (y − y³/6) ≤ 3y⁵/40;
   * 0 ≤ (1+y²)^{−½} − (1 − y²/2) ≤ 3y⁴/8.

   Both series are alternating with terms decreasing in magnitude. Hence:
   * τ = x − εx³/24 + δτ with 0 ≤ δτ ≤ 3ε²x⁵/640;
   * w = 1 − εx²/8 + δw with 0 ≤ δw ≤ 3ε²x⁴/128;
   * 0 ≤ x − τ ≤ εx³/24 ≤ x/6.

   Writing w r(τ) − r(x) = (w−1)r + w r′(τ−x) + wR₂, with |R₂| ≤ ½ sup_{[τ,x]}|r″|(τ−x)², gives
   w r(τ) − r(x) = ε g₁(x) + ε² h(x). Here g₁ = −(x²/8)r − (x³/24)r′, and the remainder comes from three terms:
   δw·r, [w δτ + (1−w)εx³/24]·r′ with 1 − w ≤ εx²/8, and w R₂. This gives the bound h̄ above.
4. **The g₁ identity.** g₁ = −(x³r/24)′ = (1/24)∫₀¹χ″(p)cos(px)dp. Proof: 3x²Φ + x³Φ′ = x²∫χ cos(px)dp by one
   integration by parts (ψ(0) = ψ(1) = 0). Two more integrations by parts give x²∫χ cos(px) = −a₁ − ∫χ″cos(px). The
   boundary data are χ(1) = χ′(1) = 0 and χ′(0) = a₁, all asserted exactly in `t0_common`.
5. Finally (2/π)∫₀^{π/2}cos(pb̂ sin φ)dφ = J₀(pb̂). ∎

Checks [DV], `t0_checks.py`:
* the g₁ identity agrees to 10⁻¹⁰ (Arb, 8 points);
* Lemma 2 against direct Mehler–Dirichlet quadrature of Λ²C_k agrees to ≤ 10⁻⁵ relative at Δ_g = 10⁴, and to the expected
  O(ε) at 300.
* The one mismatch (Δ_g = 10⁴, b̂ = 1.8·10⁴) is outside `ls_core.OmFlat`'s documented quadrature range (x ≲ 4000), not a
  failure of the lemma.

### 2.2 Lemma 3 (enclosure of F_k(b) := ∫χ″J₀(bp)dp) [P + certified]
* **At a point:** the J₀ power series with exact moments μ₂ₘ = ∫p^{2m}χ″. The tail beyond M = 2b + 60 is bounded by
  (4/3)·(next term)·∫|χ″|.
* **On a cell** [b_c − h, b_c + h]: F(b_c) + F′(b_c)[−h, h] + [−1, 1]·(h²/2)∫|χ″|p². This uses
  |J₀″| = |J₀ − J₂|/2 ≤ 1, since |J_n| ≤ 1 from Bessel's integral.
* **Tail:** χ″ = p f with f(1) = f′(1) = 0. The Bessel equation gives ∫f J₀(bp)p dp = −b⁻²∫(pf″ + f′)J₀(bp)dp, hence
  |F_k(b)| ≤ CT_k/b², with CT₂ ≤ 2.0225·10⁵ and CT₃ ≤ 2.7944·10⁵.

Constants (Arb): ∫|χ″| ≤ 743.85 / 1643.77 and ∫|χ″|p² ≤ 341.52 / 616.34 (k = 2 / 3). The |P| bounds use exact Taylor shifts on
2000 p-cells.

### 2.3 Lemma 4 (H_k := sup_x h̄_k(x)) [certified]
|r^{(n)}| (n ≤ 2) is bounded by different majorants on different x-ranges:

| x-range | majorant |
|---|---|
| x ≤ 10⁻³ | crude ∫\|ψ\| + c_n a₁/x^{2+n}, closed form |
| [10⁻³, 0.5] | min of the crude and the exact finite integration-by-parts majorant (`gl7`'s `RB`, copied) |
| [0.5, 80] | point Arb values of r, …, r⁗ (cosine series with exact moments) plus a third-order Taylor remainder, on cells of width 0.01 |
| [80, 4000] | the RB majorants, using the Taylor cells for sup r″ wherever [5x/6, x] reaches below 80 |
| x ≥ 4000 | every term of the asymptotic majorant is non-increasing (x^{3−j}, j ≥ 3; x^{6−a}, a ≥ 6), so its value at 4000 bounds it |

Result: **H₂ ≤ 12124 and H₃ ≤ 12273**. The asymptotic value is 3059 / 3508. So ε₀(1+c)(2H₂ + H₃) ≤ 0.044 at Δ_g = 1000 and
0.49 at 300, in Λ⁻² units.

### 2.4 Lemma 5 (flat-Regge remainder, for b̂ > 2Λ) [P]
Write ψ = p f, with f(1) = f′(1) = 0 and f₁ = f′(0) = ψ″(0)/2. Then Ω_flat(b) := ∫ψJ₀(bp)dp = −f₁/b³ + E_fl(b), with

    |E_fl(b)| ≤ b⁻³[ |f₁|(1 + 2/b) + |g(1)| + (|g′(1)| + ∫|pg″ + g′|)/b ],      g = Σ_{j≥2} j² f_j p^{j−2}.

*Proof.*
1. Two applications of the identity ∫ f p y dp = −b⁻²{[boundary] + ∫(pf′)′y dp}, with y = J₀(bp) and (py′)′ = −b²py.
2. For f the boundary terms vanish, and pf″ + f′ = f₁ + pg.
3. For g the boundary terms are g(1)bJ₁(b) + g′(1)J₀(b).
4. ∫₀¹J₀(bp)dp = b⁻¹(1 − ∫_b^∞J₀) uses ∫₀^∞J₀ = 1 (DLMF 10.22.43). Also
   ∫_b^∞J₀ = −J₁(b) + J₀(b)/b − ∫_b^∞J₀/t², so |·| ≤ 1 + 2/b.
5. Only |J_n| ≤ 1 is used. ∎

The −f₁/b³ term is exactly LR∞'s flat tail −a₂/(2b̂³), so it cancels in C_k.

---

## 3. The certificate

Fix c and Λ_min = Λ(DGMIN). Every bound below is valid for **all** Λ ≥ Λ_min.

**Part K: 16 ≤ b̂ ≤ 2Λ (`t1_curv.py`).** For each pattern g, σ = ±1 and u ∈ {0, 1}:

    K_g = S_g·M₁low(β_max) + (1−cs₂)·2·F₂/24 + (1−cs₃)uσ·F₃/24 − ε₀[(1−cs₂)·2H₂ + (1−cs₃)uH₃],
    S_g = (1−cs₂)·2a₁₂ + (1−cs₃)uσa₁₃ > 0.

* u ∈ {0, 1} suffices because K_g is affine in u.
* M₁low is the certified minimum of M₁ over β ≤ β_max = 2 asinh(b̂_hi/2Λ_min): the lrt2 cells, shaved by 2⁻⁴⁰ as in
  `gl_common`, and on (0, 0.1] the [P] bound of GAP_LARGESPIN §3.2.
* b̂ is covered by 2840 cells of width 0.1 on [16, 300]. The tail [300, 2Λ] uses Lemma 3's CT/b² and M₁ ≥ 0.030442 on
  β ≤ 2 asinh 1.

By Lemma 2, Λ²·(curvature part of pattern g) ≥ K_g > 0.

**Part L: b̂ > 2Λ, i.e. β > 2 asinh 1 (`t2_lrinf.py`).**
* C_k = a₁M₁/Λ² + E_k − E_fl,k (LR∞'s form, with Lemma 5).
* LR∞ uses the envelope R̂ anchored at Λ_a = Λ_min. Then Λ²(2|E₂| + |E₃|) ≤ Λ_a²N̂(β; Λ_a) for every Λ ≥ Λ_a
  (GAP_LARGESPIN §1.2, [P + certified]). It is evaluated on the 564 lrt2/gl cells in (1.76, 64], plus `tail_majorant` for β ≥ 64.
* Lemma 5 gives Λ²|E_fl| ≤ e(b̂_min)/(Λ_a s³), with s = 2 sinh(β/2), decreasing in Λ.
* The check is S_min M₁ > (1+c)[LR∞ remainders + flat remainders], with S_min = (1−c)·2a₁₂ − (1+c)a₁₃ = 241.89.

**Part T: threshold on b̂ ≤ 16 (`t1_curv.py`).** Drop the M₁ term, which is ≥ 0. Then
T_gap := max over cells, g, σ and u of −K_g (without M₁), and T_abs := T_gap/Λ_min². Then (curvature part) ≥ −T_abs.

**Part F: b̂ ≤ 16 (`t3_flat.py`).** FF_g(c) ≥ T_abs on every flat state with b̂ ≤ 16. The computation uses unnormalised NF02,
so the shift is k₀·T_abs:
* **(a) J < 100:** the exact Lemma-P polynomial in u (`rcert.poly_u`) minus the shift, on u ∈ [0, min(1, (16/(2J+1))²)].
  Positivity by exact Descartes/bisection. 584 polynomials (`++` even J, `+−` J ≥ 4).
* **(b) J ≥ 100:** CERT_CLOSE's region-S machinery, copied with the shift subtracted from the Taylor model's constant term.
  It uses (b̂, u) Taylor/Cauchy cells with u ≤ ((b̂ + h)/201)², and the small-b̂ series for b̂ ≤ 1. 672 cells, 4 patterns × 2 kernels.

Then Lemmas 0 and 1 give T0. Results:

| run | c | Part K: cells / min K_g (Λ⁻² units) / tail | Part L: max (1+c)·rem/(S_min M₁) | Part T: T_gap, T_abs | Part F: shift used | verdict |
|---|---|---|---|---|---|---|
| Δ_g ≥ 1000 | 0.199 | 2840 / ≥ 4.92 (b̂ = 19.6, `++`, s = (+,−), u = 1) / ≥ 6.94 | 0.0457 (β = 1.76); β ≥ 64: 2·10⁻⁷ | 28.99; 2.9075·10⁻⁵ | 3·10⁻⁵ (PASS, 584 polys + 672 cells) | **PASS** |
| Δ_g ≥ 300 | 0.199 | 2840 / ≥ 4.58 / ≥ 6.49 | 0.1905 | 29.44; 3.3036·10⁻⁴ | 3.31·10⁻⁴ (PASS) | **PASS** |
| Δ_g ≥ 1000 | 0.1999 | 2840 / ≥ 4.92 / PASS | 0.0458 | 29.01; 2.9097·10⁻⁵ | 2.9097·10⁻⁵ (PASS) | **PASS** |

(Logs: `cache/t{1,2,3}_{1000,300,1000_c1999}.log`; JSON/TXT `cache/t*_<DGMIN>_<c>.*`. The Part-F shifts in the first two rows
are deliberately above T_abs: the slack, 9·10⁻⁷ and 6·10⁻⁷, absorbs the realised-weight term of §4.2.)

**Why it works, in one line.** In gap units the curvature is Λ⁻²·[a₁M₁ + F(b̂)/24] + O(Λ⁻⁴). It is negative only for
b̂ ≲ 13 (§5), where the flat action is O(1), and it is ≥ +4.6/Λ² beyond. So positivity is limited only by the sharp flat 1/5.

---

## 4. Realised weights (M20) and the two coordinator inputs

### 4.1 M20: status
* The every-gap theorem (§0.1) is about the de-shifted leading-Regge action, which is defined at every Δ_g through ψ_k(ν/Λ).
* For the **realised** Ω^EE, the k = 2 weight X_Λ must exist at each gap. Published continuum weights exist only at 300 and
  1000, with ε_W = 3.3·10⁻⁶ / 1.1·10⁻⁶ [DV] (CERT_CLOSE §3.2).
* A Λ-uniform continuum family (ψ₂ + S/Λ² + O(Λ⁻⁴) with an explicit remainder, or certified Abel solves on Λ-balls) was
  **not** attempted: it is a separate numerical-analysis task, and its ε_W would again be [DV] unless the 3·10³-term RBF
  residual were interval-integrated.
* **M20 open.**

### 4.2 Corollary (realised weights at Δ_g = 1000 and 300) [certified + DV(ε_W)]
CERT_CLOSE gives |image − target| ≤ ε_W P_{−½}(cosh β) for the realised k = 2 leading-Regge image. In hybrid units this
changes H − cN_s by at most 2(1+c)ε̂_W P_{−½}(cosh β), with ε̂_W = (ε_W/∫|ψ₂|_ν)·∫|ψ₂| ≈ 1.1·10⁻⁹ (1000) / 1.1·10⁻⁸ (300) and
P_{−½} ≤ 1. Against the certified margins:
* **Part F:** the slack between the shift used and T_abs (§3) is ≥ 6·10⁻⁷, ≫ 2.7·10⁻⁸.
* **Part K:** the margin is ≥ 4.58/Λ² ≥ 5·10⁻⁵ (300), ≫.
* **Part L:** the margin is ≥ 0.81·S_min M₁/Λ². With M₁ ≥ (√2/π)K₀ℓ/β (β ≥ 64) and P_{−½} = (√2/π)∫K ≤
  (√2/π)e^{−β/4}(C₁ + C₂β)/√(2 sinh(β/2)) (LRINF_TIGHT §1.5), the ratio is P_{−½}/M₁ ≤ (C₁ + C₂β)β/(√2ℓ) ≈ 1.44β². The
  corollary therefore holds for **β ≤ β_W ≈ 220** (both gaps), i.e. J + ½ ≤ m sinh(110).

So at Δ_g ∈ {300, 1000}, T0 with c = 0.199 holds for the realised functional on every heavy state with β ≤ 200, given ε_W [DV].
The cap is an artefact of the crude |W − Φ| ≤ ε_W, which does not decay in t.

### 4.3 The reviewer's relative de-shift lemma (coordinator input 1)
The lemma: if the residual r(ν) = (effective weight) − (target) has r(0) = r′(0) = 0 and ∫|r⁗| < ∞, then
|r̂(t)| ≤ min(A, B/t⁴). With ∫₀^∞r̂ = 0, this gives |S[r](β)| = O(K₀β⁻³).
* The four integrations by parts are correct: boundary terms r(0)sin/t, r′(0)/t² and r″ sin/t³ vanish, leaving r‴(0)/t⁴.

On the three questions:
* **(a) Floor scaling.** Confirmed for our weights. ψ_k(0) = 0 exactly (ψ_k = p(1−p)⁵X_k), and the leading-Regge floor is
  A₁M₁(β)/Λ² with M₁ ≥ (√2/π)K₀·0.874/β (LRINF_TIGHT §1.5), i.e. ∝ e^{−β/2}β⁻¹. Part L confirms the floor dominates: the LR∞
  and flat remainders are ≤ 4.6 % (Δ_g ≥ 1000) / 19 % (≥ 300) of it on every β-cell, and ≤ 2·10⁻⁷ for β ≥ 64. A weight with
  ψ(0) ≠ 0 would add a ψ(0)Λ⁻¹∫P_{−½+iν}dν-type term ∝ e^{−β/2}, as the reviewer suspected; not our case.
* **(b) The published continuum weights do not satisfy the hypotheses.** In CERT_CLOSE's shifted form r(μ) = E(μ) − ψ₂(|μ|/Λ):
  * r(0) = E(0) = 2Re[χ(−i)R₋(−i)] is not constrained to vanish (it is small only at the fit level);
  * r′(0⁺) = 0 would need the fitted log-basis coefficients to reproduce the kink κ = ψ₂′(0) exactly; least squares gives
    this only approximately;
  * the kink-carrying basis functions (ν²+1)log(ν²+1)e^{−ν²/2}, evaluated at ν = μ ∓ i, contain μ² log|μ|. So r‴ ∼ 1/μ and
    ∫|r⁗| = ∞, i.e. **B = ∞**.

  Meeting the hypotheses needs a new weight whose residual is C⁴ on [0, ∞) with r(0) = r′(0) = 0: an exact match of the kink
  and log structure, which is essentially the exact X_Λ. That is M20-type work, not done here. A and B were therefore not
  computed. A certified A alone would need interval quadrature of the RBF fit, and A alone does not help.
* **(c) Effect on T0.** For T0 as proved (the de-shifted action) there is no residual, and the large-impact-parameter part is
  already closed analytically and certified (Part L: all β > 1.76, with an analytic β ≥ 64 tail). The lemma would matter only
  for the realised-weight corollary §4.2, where it would remove the β ≲ 200 cap, **if** (b) could be met.

### 4.4 The monopole-free weight (coordinator input 2)
T0 is **not** limited by Q_tol. The monopole D·F(i/2) (D ≈ −990/Λ²) is a beyond-leading-Regge O(1/m²) cost of Theorem L′.
T0 uses only the exact leading-Regge image, on which the de-shift identity holds whatever the weight's pole structure. A
monopole-free weight (single pole at z = ±i/2, D = 0) therefore changes nothing in T0. It is relevant to Theorem L′'s
tolerances, and, as a new weight, to M20 and §4.2. The algebra quoted (R₋(z)/c(z+i) = −(2z+i)(2z+11i)²/(16z)) was not
re-checked here.

---

## 5. Diagnostics [DV]
* **Where the curvature is negative.** Gap-unit float scan, b̂ ≤ 200, all patterns: K_g < 0 only for b̂ ≤ 12.5 (c = 0),
  13.0 (c = 0.2). The certified T_gap is attained at b̂ ≈ 2.5, pattern (−, −), σ = +1, u = 1.
* **Flat margin on b̂ ≤ 14** (float, J < 400): min(𝔉 − cN) = 0.0025 (`++`, J = 4, u = 1) / 0.0091 (`+−`, J = 6, u = 1) at
  c = 0.199. Against T_abs this is ×85 (1000) / ×7.5 (300) at c = 0.199, which is why c = 0.1999 passes at 1000 but would not
  at 300.
* **Ratio scan** (`t0_checks.py`, 3 gaps × 18 (h, J) × 7 masses): min H/N_s = 0.19999969, consistent with c ≤ 0.2000000103.
* **Not pushed below 300.** At Δ_g = 100 the LR∞ remainders scale ∝ Λ⁻² (×~100 vs 1000, so Part L's 0.046 → ~4), and T_abs
  grows ×100. A lower gap would need a different large-b̂ argument.

## 6. Reproduce
    cd ads4/code/t0cert && bash run_all.sh        # ≈ 20 CPU-min, ≤ 2 processes

| file | content | grade |
|---|---|---|
| `t0_common.py` | exact χ, F_k enclosures, h̄/H_k, M₁ lower bounds, Lemma-5 constants | — |
| `t1_curv.py` | Parts K and T | [certified] |
| `t2_lrinf.py` | Part L | [certified] (LR∞ [P + certified] as cited) |
| `t3_flat.py` | Part F (exact Descartes + region-S cells with shift) | [certified] |
| `t0_checks.py` | identity and ratio checks | [DV] |

## 7. Caveats
1. **Cited inputs.** Lemma 1 (CERT_CLOSE §1; its corner code is marked unreviewed there), LR∞ with the R̂ anchor
   (LRINF_TIGHT, GAP_LARGESPIN §1.2), the lrt2 M₁ cells, and the small-β M₁ bound (GAP_LARGESPIN §3.2). All are graded
   [P]/[certified] in their documents and were not re-derived, except that M₁'s definition was matched (§2.1, step 2).
2. **Mehler–Dirichlet** and ∫₀^∞J₀ = 1 are classical [Q]. Everything else in §2 is proved here.
3. The T0 statement covers u ≤ 1 for every J and therefore also the unphysical strip Λ ≤ m < m_gap(J). This is stronger than
   needed, and it is what makes the result automatically uniform in Δ_g: the curvature bounds use only Λ ≥ Λ_min.
4. T0 says nothing about δ (beyond leading Regge). It does not upgrade any (Pos) grade of the full action.
