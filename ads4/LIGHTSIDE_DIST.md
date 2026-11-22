# The light side of Result Q as one functional: term-wise s-channel action, the Regge arc, and a check on EIN's Witten-diagram family

*22 Nov 2026, light-side agent. Unrefereed, not committed. Answers PRL referee report §M3 (and the M7 "Einstein instance" wording) and
RESULT_Q_EINSTEIN_CHECK (EIN) §5. No existing file was edited. New code: `ads4/code/resultq/lsd_mellin.py`, `lsd_core.py`, `lsd_pos.py`,
`lsd_run.py`, `lsd_hairpin.py`; outputs `resultq/cache/lsd_*` (git-ignored). Read-only imports: `ein_dbar`, `ein_run`, `og_free`. Literature
read from source (arXiv LaTeX): CLS-D 2311.04271 (§2.2, eq. `omega_on_light`), CHMRS 2106.10274 (§`sec:lightAction`, App.
`app:omegalight`, eq. `eftformula`). Compute ≈ 1.5 CPU-h, ≤ 2 processes (one brief overlap of 3 for ~1 min), 1 thread, < 1 GB each.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked; **[NE]** numerical evidence; **[M]** model;
**[Q]** quoted; **[A]** assumption; **[I]** open.

---

## 0. Summary

**Outcome: (A), in a corrected form.** There is one functional, it is the one the heavy side already uses, and on tree terms it splits
into "Lemma A value + Regge arc". The arc is zero below the borderline and is the polynomial (P) value at the borderline. Two parts of
the current framing are wrong and must change.
* The J₀ = k value does **not** live at the collinear set a = b, nor in the ν/de-shift integration.
* EIN's "J₀ = 4 contact annihilated by Ψ₄" is **not** an R⁴ analogue. In the product T′ = σ ⊗ T_free that contact is two units below
  the borderline.

1. **The functional [Q + D].** Ω_k is CLS-D's own sum-rule object (2311.04271 eq. `omega_on_light`; CHMRS App. `app:omegalight`):
   * the swapped (double-commutator) functional applied **term by term to the s-channel (12)(34) OPE**, with the 1↔3 image;
   * at every fixed celestial pair (a, b) with a ≠ b, before the C_ν/χ pairing.

   Heavy blocks enter it block by block, so (i) holds by construction. Σ_O Ω[O] = 0 is the full-correlator statement plus the
   OPE/(corner, ν) interchange H_K, which is unchanged.
2. **Its value on a tree term X [D, exact Mellin algebra].** In the Mellin representation on the swapped contour C × C:
   * the kernel K(δ₁₂, δ₁₃) is entire in δ₁₂;
   * term-wise s-channel = sum of the left δ₁₂-residues (double-trace double poles δ₁₂ = −n, plus single-trace poles);
   * the whole-function fixed-angle integral (EIN's E⁺, D_LR) = the straight Mellin line.

   Hence **Ω_k[X] = E⁺_{D_LR}[X] + 𝒜_k[X]**, where 𝒜_k is the δ₁₂ contour at infinity, i.e. the u-channel Regge arc:
   * 𝒜 = 0 when X's effective Regge spin n_X < k;
   * at n_X = k it depends **only on X's top Regge coefficient**. That is CHMRS's `eftformula` and the content of (P).
3. **Why D_LR cannot see it [DV].** On C × C the 1↔3 (anti)symmetrised whole-function integrand cancels as an *integrand identity*
   (EIN's pair identities). So the D_LR value of any microcausal term is 0 at every k:
   * at the borderline the fixed-angle integral itself diverges: joint-corner exponent −3, e.g. EIN's c4 at k = 2;
   * the light value is carried entirely by the arc.
4. **Check on EIN's engine [DV]** (`T′ = σ ⊗ T_free`, three angle/polarisation sets P1–P3, k = 2, 3, 4, 5). Term-wise values are fine-grid
   residue sums to n = 150 with power-law tails; "sym" = E(a,b) + (−1)^k E(b,a), summed over the six free classes.
   * **Validation.** The Mellin form reproduces EIN's D-function engine:
     * pointwise on C × C to 1·10⁻¹⁴;
     * on whole-function families for every subcritical piece: c0–c3, xs, xt to 10⁻⁷ – 10⁻¹¹ of S; c4 at k = 4 and xu to ≤ 3·10⁻⁵,
       limited by the residue tail and quadrature (k = 2, 3, 4, 5).
   * **Below the borderline the term-wise value is zero.** Pieces: c0, c2, c3, xu, xs + xt at k = 2–5; c4 and xu4 at k = 3, 4, 5.
     * The sym partial sums fall like N^{−(1.2 – 4)} and cancel against the fitted tail to ≤ 10⁻⁶ of S.
     * There is no collinear, boundary or corner remainder: term-wise = D_LR, family by family.
   * **The borderline in this family is k = 2 for c4** (effective Regge spin n = J₀(σ) − 2 = 2). The term-wise light value is
     **non-zero, finite and convergent**:

     | | c4 | r4 (its top Regge term only) |
     |---|---|---|
     | P1 (Δφ = 1.2) | −8.7642·10⁻³ | −8.7642·10⁻³ (rel. 4.8·10⁻⁷) |
     | P2 (Δφ = 4.1) | −2.6187·10⁻⁴ | −2.6187·10⁻⁴ (6.7·10⁻⁷) |
     | P3 (Δφ = 0.3) | −5.7015 | −5.7015 (3.2·10⁻⁷) |

     * This holds although c4 and r4 differ by up to 12× in each single ordering.
     * The D_LR (whole-function) value is 0 by the integrand identity, and divergent before symmetrisation.
   * **Borderline exchange (an H7 analogue).** xu4 is a u-channel exchange with a degree-4 polynomial residue (Regge spin 2). Its value
     at k = 2 equals that of its top Regge term r4u = (δ₁₂ − δ₁₄)⁴/(16(δ₁₃ − 1)) to 6–8·10⁻⁸ at P1–P3 (−3.4123·10⁻³ at P1).
5. **Consequences.**
   * **Result Q.** The light side "D_LR on J₀ < 4 terms + (F)+(P) on J₀ ≥ 4 terms" is the value of the single functional Ω₄ on the
     light OPE content.
     * The J₀ < 4 part is Lemma A plus "no arc" [D; DV toy].
     * The J₀ = 4 part is the arc, a linear functional of the top Regge coefficient only [D; DV toy].
     * The bounds do not move numerically.
     * The only part still at [D-sketch] is the **normalisation** arc(R⁴) = 2g₄, i.e. (F): the spinning flat limit of the arc.
       Unchanged, not tested here.
   * **Result D.** The Einstein exchange has n = 2 = k under Ψ₂, so H7's value is the arc of the same term-wise functional. This is the
     one sentence M3 asks for; CHMRS compute their graviton light term exactly this way.

---

## 1. Definition

### 1.1 The functional [Q + D]
For fixed celestial angles a ≠ b, polarisations, and k:
* the single-ordering swapped form is `E⁺[G](a,b) = ∮_C∮_C t₁^w t₃^w f_k G` (og_free / EIN conventions);
* the block action is `ω_k[O](a,b) := E⁺[G^s_O](a,b) + (−1)^k E⁺[G^s_O](b,a)`, one s-channel conformal block at a time.

**Ω_k[G](a,b) := Σ_{O ∈ (12) OPE of G} p_O ω_k[O](a,b).** The de-shifted functional is `Ĉ^{ds} = ∫dν χ(ν) ∫dφ₁dφ₃ C_ν Ω_k / N_com`.

* **This is the definition CLS-D and CHMRS use** [Q]. The sum rule is "Σ over the s-channel conformal block decomposition of ω[G^s]",
  split into heavy and light (CLS-D eqs. `omega_on_light`, `heavyaction_scalar_definition`; CHMRS `ω|_light = Σ_{τ≤Δgap} p ω[G^s] +
  subtractions`).
* **The heavy side is this object by construction** [D]. Result Q's per-state images are ω_k on individual heavy blocks, ν outside,
  with the same χ.
* **The sum rule.** Σ_O p_O ω_k[O] = 0 follows from:
  * the whole-correlator zero (Lemma A on 𝒢_full, n < k);
  * plus the interchange of the full OPE sum with the (corner, ν) integrations, i.e. **H_K**.

  Nothing new is assumed.
* **There are no separate collinear or ν-boundary terms.**
  * Ω_k is defined pointwise at a ≠ b.
  * The a = b set enters only through the C_ν pairing, as for any heavy block.
  * Integrability of the light function near a = b against C_ν is the open item of M7 problem 3. It is not addressed here [I].
* **The light side** is −Σ_heavy = Σ_light p_O ω_k[O]: T, and the [TT]_{n,J} with their full O(G) data.
  * By linearity this equals Σ over tree terms X (Einstein, contacts, heavy-exchange double-trace parts) of **Ω_k[X] := the term-wise
    s-channel value of X**.
  * This is the H_EFT split of FH §4.5 and is unchanged.

### 1.2 Term-wise vs whole function: Mellin form on the swapped contour [D]
For a correlator with a Mellin representation (EIN's σ-part: `I_Δ = ½·(2πi)⁻²∫[dδ]ΠΓ(δ_ij)X_ij^{−δ_ij}`; the ½ checked in §3.1):
* **The product of powers.**
  * `ΠX^{−δ} = (X₁₄X₂₃)^{−2}(U/V)^{−δ₁₂}V^{δ₁₃}`.
  * On C every X_ij has its principal branch: 2t − 1 = −e^{iθ}, arg t ∈ (−π/2, π/2), Re(2t+1) > 0.
  * The MB representation evaluated with these branches *is* EIN's physical-route continuation [DV, §3.1].
* **The kernel.** `K[c](δ₁₂, δ₁₃) := ∮_C∮_C H[c]·(U/V)^{−δ₁₂}V^{δ₁₃}`, with `H = t₁⁵t₃⁵ f_k G_free[c](X₁₄X₂₃)^{−2}/2`.
  * It is **entire in δ₁₂**, because |U/V| ≤ 1 on C and |X₁₂| = |X₃₄| = 1.
  * It is separable in (t₁, t₃).
* **Three evaluations of the same object:**
  * **whole function (EIN, D_LR):** t outside;
  * **line:** Mellin outside (Fubini partner);
  * **term-wise s-channel:** close δ₁₂ to the left. The poles are Γ(δ₁₂)² at δ₁₂ = −n (the s-channel double-trace families, including
    their O(λ) "anomalous-dimension" logs), plus the single-trace pole of an s-exchange. Each residue is a convergent C × C integral of
    integer powers, so it needs no regularisation and is the same on the hairpin as on C.
* **Cauchy [P].** At fixed δ₁₃:

      line − term-wise = −(arc of Γ(δ₁₂)²Γ(δ₁₄)²M K on the left semicircle at infinity)

  After 1↔3 (anti)symmetrisation and the δ₁₃ integral, the D_LR whole-function value is 0 (§1.3). Hence
  **Ω_k[X] = 𝒜_k[X] := (symmetrised) contour at infinity**.
* **Power counting [D].** Along the arc, M ~ (Regge growth) and K ~ (corner behaviour of H):
  * if X's Regge spin n_X < k, the arc integrand is o(|δ|⁻¹) and 𝒜 = 0;
  * at n_X = k it is ~ c/δ₁₂, and 𝒜 ∝ (the coefficient of the top Regge monomial of M) × (a universal Regge kernel of Ψ_k);
  * lower monomials, which have n < k, drop out.

  This is the position-space statement of CHMRS `eftformula` ("−ω|_light = ¼∮dr r^{k−2}𝒢⁺_EFT", the r^{1−k} Laurent coefficient)
  [Q], and it is exactly the structure (P) asserts.
* **(F)** is the separate statement that, for the spinning R⁴ ⟨TTTT⟩ at n = 4 = k, the universal kernel's flat limit returns 2g₄ in
  Stage 0's normalisation. That statement is not touched here [D-sketch].

### 1.3 Why the fixed-angle integrand never shows the value [DV]
* **The identity.** EIN's pair identities `E[p×c](a,b) + (−1)^k E[π(p)×π(c)](b,a) = 0` hold node by node on the symmetric C × C grid,
  given single-valuedness on C × C (EIN §2.5, order independence 10⁻¹¹). So the symmetrised D_LR value is 0 for every microcausal term,
  for every k, with any 1↔3-symmetric regularisation of a divergent corner.
* **Example at the borderline.** EIN's own engine at c4, k = 2 gives E(a,b) = −1.73429·10⁻³, E(b,a) = +1.73435·10⁻³, with sym
  6·10⁻⁸ (n = 60).
  * The integral is divergent: EIN §3.5 joint exponent −3.
  * Those single values are quadrature artefacts.
* **Consequence.** "Lemma A holds pointwise for J₀ ≤ k" (EIN §3.5) is true but carries no information about light values. The term-wise
  functional differs from D_LR exactly by 𝒜.

---

## 2. What changes in the reading of EIN

* **The relevant spin is the product's.** EIN's J₀ is the σ-part's Mellin growth at fixed δ₁₃. The Regge spin that matters for the sum
  rule is that of the product σ ⊗ T_free in light-ray polarisation.
* **Where the borderline sits [DV, §3.3].** Term-wise values vanish at k = 3, 4, 5 for c4 and are non-zero at k = 2. The J₀(σ) = 2
  pieces vanish at k = 2. Hence **n_eff = J₀(σ) − 2**: c4 → 2, c2/c3/xu → 0, xs/xt → −1, c0 → −2.
* **The joint-corner exponent rule.**
  * With n_eff, EIN's measured exponents (§3.5: c0 +3, xs/xt +2, c2/c3/xu +1, c4 −1 at k = 4; −1/−3 at k = 2) all satisfy
    **e = k − 3 − n_eff** [NE, empirical in this family].
  * So the borderline is a *power-divergent* (e = −3) fixed-angle integral, not the integrable e = −1 case.
* **Consequences for EIN's text.**
  * "The J₀ = 4 contact under Ψ₄" is n = 2 < 4: two units subcritical. Its zero is correct, and its (P) value is also 0.
  * The "J₀ = 2 exchange at k = 2" is n = 0: subcritical.
  * EIN §5's implication "light-side values of J₀ = k terms live at the collinear set or in the ν integration" is **withdrawn**. They live
    in the Regge arc of the s-channel OPE sum, at every angle pair (§3.3: non-zero at Δφ = 0.3, 1.2, 4.1).
* **For the actual R⁴ ⟨TTTT⟩ under Ψ₄**, n = 4 = k is fixed by the flat dictionary [Q + D]:
  * CLS-D §6.5: the k = 4 graviton sum rules are the lowest-subtracted `−∮ds(s−t)/(−st)^{(k−2)/2}M^{(i)}`;
  * R⁴ has constant M^{(i)}, hence a non-zero arc.

  So it *is* at the borderline, and its value is its arc. The Einstein term under Ψ₄, Ψ₅ (n = 2) is subcritical, which is Lemma A.

---

## 3. Checks

### 3.1 Representation (`lsd_mellin.py`; `scratchpad` diagnostics reproduced by `lsd_run`)
* **Normalisation.** At a Euclidean point, MB/EIN = 2.0000000000 for c0, c2, c3, c4 (|dev| ≤ 8·10⁻¹¹). So I_Δ = MB/2.
* **Pointwise on C × C, Δφ = 1.** Six node pairs (θ₁, θ₃) ∈ {(1,2), (0.3,0.5), (0.3,6), (5.5,5.8), (3,0.05)} against EIN's
  physical-route D-function values (mp jets):
  * c0, c4, xs (contour Re δ₁₂ = 1.3), xu (Re δ₁₃ = 1.3) agree to **≤ 3·10⁻¹⁴**;
  * at (0.08, 0.1), near the joint corner, the MB trapezoid is truncation-limited (10⁻³), as expected on the Regge sheet.
* **Whole-function families, k = 4, coarse grid (n ≤ 60).** Line and term-wise vs EIN, per piece × class:
  * 10⁻⁹ – 10⁻¹⁶ of S for classes 1234/s/t/u;
  * 1243/1324: c0 ≤ 6·10⁻⁹; c2, c3, xs, xt, xu ≤ 4·10⁻⁵; c4 ≤ 8·10⁻² (line) / 1·10⁻² (term-wise). These are coarse-grid residue noise at
    n ≳ 40, removed on the fine grid below.
* **Fine-grid production.** Graded C, s₀ = 10⁻⁹, 16-point panels, ratio 1.8, 1248 nodes. Alternative grid: s₀ = 10⁻⁸, 14 points, ratio
  2.2.

### 3.2 Term-wise vs whole function, subcritical pieces (`lsd_run report`, P1; class sums)

| k | c0 | c2 | c3 | c4 | xu | scale |
|---|---|---|---|---|---|---|
| 5 | 4·10⁻¹² | 6·10⁻¹² | 2·10⁻¹³ | 5·10⁻¹² | 1.5·10⁻⁹ ¹ | 10⁻⁴ |
| 4 | 7·10⁻¹² | 4·10⁻¹² | 3·10⁻¹³ | 6.9·10⁻⁹ (tail est. 2.6·10⁻⁸) | 1.2·10⁻⁸ ¹ | 2·10⁻⁴ |
| 3 | 4·10⁻¹² | 7·10⁻¹² | 7·10⁻¹³ | 7·10⁻⁸ ² | 2.4·10⁻⁹ ¹ | 10⁻⁴ |
| 2 | 4·10⁻¹² | 9·10⁻¹⁰ | 5·10⁻¹⁰ | — (borderline, §3.3) | 1.3·10⁻⁸ ¹ | 10⁻⁴ |

¹ The δ₁₃ line Re = 1.3 lies 0.3 from the pole δ₁₃ = 1. With h = 0.2 the trapezoid error is ~10⁻⁴ relative. The xs/xt/xu re-run at
h = 0.05 is in §3.4.

² The same sign and size in both orderings, so it cancels in the antisymmetric k = 3 combination. EIN's n = 60 quadrature at e = −2 is
slow.

* **Symmetrised term-wise values, which are the light values, for every subcritical piece and every k = 2–5.**
  * The partial sums at N = 50 / 100 / 150 decrease like N^{−(p−1)}, with fitted residue decay p = 3.1 – 8.3.
  * They cancel against the fitted tail to ≤ 10⁻⁶ of the scale. Example: c2 at k = 2, P1: 1.72·10⁻⁹ against tail −1.63·10⁻⁹, scale
    3.2·10⁻⁵.
  * c4 and xu4 at k = 3, 4, 5 behave the same way: at k = 3, c4 gives −4.1·10⁻⁸ → −9.2·10⁻⁹ → −4.0·10⁻⁹ (tail −1.1·10⁻⁹) at scale
    2.2·10⁻⁵.
  * P2 at k = 4: c4 3.3·10⁻¹⁰ against tail −2.8·10⁻¹⁰ (scale 2.5·10⁻⁴).
  * **Lemma A' (term-wise = D_LR = 0 below the borderline) [DV].**
* **Exchanges xs, xt (with the s-channel single-trace residue taken term-wise):** see §3.4.

### 3.3 The borderline (k = 2): non-zero light values, determined by the top Regge coefficient

| tag (Δφ) | c4: sym at N = 50 / 100 / 150; tail | r4 = (δ₁₂−δ₁₄)⁴/4 | xu4 (u-exchange, residue Q_u²) | r4u = (δ₁₂−δ₁₄)⁴/(16(δ₁₃−1)) |
|---|---|---|---|---|
| P1 (1.2) | −8.7500 / −8.7609 / −8.7642 ·10⁻³; −5.4·10⁻⁶ | −8.7641689·10⁻³ (c4 −8.7641647·10⁻³) | −3.4122821·10⁻³ | −3.4122823·10⁻³ |
| P2 (4.1) | −2.6146 / −2.6177 / −2.6187 ·10⁻⁴; −1.7·10⁻⁷ | rel. diff 6.7·10⁻⁷ | −7.55698·10⁻⁵ | rel. 8·10⁻⁸ |
| P3 (0.3) | −5.6924 / −5.6996 / −5.7015; −4.9·10⁻³ | rel. diff 3.2·10⁻⁷ | −8.04989 | rel. 6·10⁻⁸ |

* **Convergence.** Residue decay p ≈ 1.9 – 2.2, so the sums converge absolutely. Extrapolated: P1 −8.7696·10⁻³, P3 −5.706.
* **Grid robustness.** The alternative grid gives −8.7585·10⁻³ at N = 100, against −8.7609·10⁻³ (2.7·10⁻⁴). Its N = 150 value is
  noise-limited.
* **The top Regge coefficient decides the value.**
  * In single orderings c4 and r4 differ by 0.7 % (P1) to 12× (P2).
  * c4 − r4 = g(δ₁₃)(δ₁₂−δ₁₄)² + g² has n = 0 < 2 and drops out of the light value to 3–7·10⁻⁷.
  * The same holds for xu4 − r4u.

  This is (P), checked on a position-space Witten-diagram family [DV].
* **Where the value sits.** By ordering and class it lies in the 1243 class: P1, (a,b) −5.977·10⁻³ and (b,a) −2.787·10⁻³. Class 1324
  contributes ≤ 5·10⁻⁹. This is the class whose (13) link carries T_free's Regge growth.
* **The D_LR whole-function value of the same term** is power-divergent. Its symmetric quadrature gives 0 by §1.3.

### 3.4 Exchanges and single-trace (block-like) terms (`lsd_run exch`)
* **Setup.** δ₁₃ step h = 0.05, |Im δ₁₃| ≤ 10, P1. The s-exchange's single-trace pole residue (the "block-like" term) is included
  term-wise.

| k | piece | term-wise (ab) | EIN whole (ab) | rel. diff (ab / ba) | single-trace residue (ab) |
|---|---|---|---|---|---|
| 4 | xs | +1.13501963·10⁻⁵ | +1.13501963·10⁻⁵ | 1.1·10⁻¹¹ / 1.0·10⁻¹¹ | −2.961·10⁻⁵ |
| 4 | xt | +1.24609949·10⁻⁵ | +1.24609949·10⁻⁵ | 1.0·10⁻¹¹ / 1.1·10⁻¹¹ | — |
| 4 | xu | −4.33640483·10⁻⁵ | −4.33640484·10⁻⁵ | 2.4·10⁻⁹ / 2.3·10⁻⁹ | — |
| 2 | xs | −2.36749928·10⁻⁶ | −2.36749856·10⁻⁶ | 3.1·10⁻⁷ / 1.1·10⁻⁷ | +2.880·10⁻⁵ |
| 2 | xt | −6.03166171·10⁻⁶ | −6.03166246·10⁻⁶ | 1.2·10⁻⁷ / 3.0·10⁻⁷ | — |
| 2 | xu | +5.67022377·10⁻⁵ | +5.67013817·10⁻⁵ | 1.5·10⁻⁵ / 1.4·10⁻⁵ (residue tail, p ≈ 3.2) | — |

* **What the table shows.** The single-trace residue is larger than the total (e.g. −2.96·10⁻⁵ against +1.14·10⁻⁵). So the agreement
  checks the split "block + double traces" term by term, not a dominant piece.
* **Symmetrised.** xs + xt gives 9·10⁻¹⁸ (k = 4) and 4·10⁻¹⁴ (k = 2), at scale 2.4·10⁻⁵ / 8.4·10⁻⁶.
* **Grade [DV].** For subcritical exchanges, including their single-trace (block) part, the term-wise functional equals the D_LR
  whole-function value. There is no collinear or boundary remainder.
* **The h = 0.2 runs of §3.2 (xu column, ¹)** differ from this only by the δ₁₃ quadrature (3·10⁻³ for xs/xt at h = 0.2).

### 3.5 Position-space hairpin attempt (`lsd_hairpin.py`) — not used
* **Contour.** A contour family C_τ that leaves t = 0 tangentially to the cut (τ = 1) and reduces to EIN's circle at τ = 0.
* **τ = 0** reproduces EIN to 10⁻⁸ (n = 40).
* **τ = 0.5** is already off by 1 % at n = 40, and τ = 1 at n = 60 by O(10²) on subcritical pieces. That is quadrature-limited, as in
  EIN's κ = 0.5 note.
* **Verdict.** A direct position-space hairpin evaluation of the borderline light value would need n ≳ 200 with mp jets per node (days).
  Not attempted further [I]. The Mellin term-wise evaluation needs no such contour.

---

## 4. Regrading

| claim | before | after |
|---|---|---|
| Light side of Q "D_LR on J₀ < 4 + (F)+(P) on J₀ ≥ 4" | definition [D-sketch] (UT, Q §2, EIN §5); referee M3: patch | **derived**: it is Ω₄ = term-wise s-channel functional (the heavy-side object) on the light OPE content. Its value on a tree term is D_LR + 𝒜₄, with 𝒜₄ = 0 for n < 4 and 𝒜₄ = arc for n = 4 **[D]**, exact Mellin–Cauchy plus Regge power counting; **[DV]** on the σ ⊗ T_free Witten family, k = 2–5, 3 angle sets |
| Lemma A + D_LR for n < k tree terms | [D] | **[D]**, with its role restated: D_LR is the term-wise value because the arc vanishes. Term-wise = D_LR checked [DV] (§3.2) |
| (P): polynomial light value fixed by the top coefficient | [D-sketch] | structure **[D + DV]** (c4 = r4, xu4 = r4u to ≤ 7·10⁻⁷). The per-ν form Σ_j R_j p((ν−ij)²) still needs (F) |
| (F): flat normalisation arc(R⁴) = 2g₄ (`−Ĉ₄^{ds}` = 2g₄) | [D-sketch] | **[D-sketch], unchanged**. The R⁴ ⟨TTTT⟩ itself was not evaluated; the toy has no flat normalisation |
| "J₀ = k light values sit at a = b or in the ν integration" (EIN §5) | [D-sketch] | **withdrawn**. They sit in the Regge arc of the s-channel OPE sum at every a ≠ b [DV] |
| "Pointwise zero for J₀ ≤ k; the J₀ = 4 contact is an R⁴ control" (EIN §3.5) | [DV] | the zero is an integrand identity on C × C valid for all k [DV]. In the product the Regge spin is J₀(σ) − 2, so c4 under Ψ₄ is subcritical; the R⁴ analogue is c4 under Ψ₂ **[DV/NE]** |
| R⁴ is at the borderline of Ψ₄ (n = 4 = k) | implicit | **[Q + D]** (CLS-D flat dictionary) |
| H7 (D's k = 2 Einstein pole) is the value of the same prescription | asserted | **[D]**: Einstein exchange has n = 2 = k at Ψ₂, so its light value is its arc (CHMRS compute the graviton term this way [Q]). Toy analogue xu4 = r4u at k = 2 **[DV]** |
| g̃₄ | "minus the heavy sum", identified by (F)+(P) | the arc of Ω₄ on the light content, equal to minus the heavy sum given H_K [D]; = 2g₄ requires (F) [D-sketch] and H_EFT(b) [A] |

---

## 5. Consequences

### 5.1 Result Q
* **No number changes.** The value assigned to every light term is the same as before. Lemma A gives 0 below the borderline, and the
  arc gives (P) at it.
* **The referee's "possibly W" is removed.** The functional the sum rule uses does not annihilate R⁴: D_LR's zero on R⁴ (if it occurred)
  would be the integrand identity of §1.3, not the light value.
* **Recommended wording for the letter / UT §1.** "The light side is Σ_{light O} p_O ω₄[O], the same block functional as the heavy
  side. On each tree term it equals the fixed-angle (D_LR) value plus the Regge arc (CHMRS eq. 5.3/`eftformula`). The arc vanishes for
  Regge spin < 4 (Lemma A) and, for the R⁴ contact (spin 4), depends only on its top Regge coefficient; its flat normalisation 2g₄ is
  (F)."
* **Hypothesis list.**
  * Row 4 (Lemma A + D_LR) keeps [D] and gains "arc = 0 below the borderline [D; DV]".
  * Row 5 ((F)+(P)) becomes "(P) structure [D + DV], (F) normalisation [D-sketch]".
  * H_K and H_EFT are unchanged.
* **The lower bound** uses the same light value 2g̃₄ and is unaffected.

### 5.2 Result D
* The H7 value is the k = 2 arc of the Einstein exchange under the same term-wise functional, since Einstein has n = 2 = k [D].
* The sentence owed in M3: "the k = 2 Einstein pole is the Regge-arc value of the same term-wise s-channel functional that is applied to
  the heavy blocks; the fixed-angle integrand does not contain it".
* Whether D's specific H7 normalisation equals this arc is D's existing derivation (CHMRS-type graviton Regge limit). It is not
  re-derived here.

### 5.3 What the referee's fix list now has
* "Define Ω^Q as a distribution with its collinear and ν-boundary terms": done (§1.1). It has none; its extra piece is the Regge arc,
  defined at every a ≠ b.
* "Show on one explicit spinning J₀ = 4 position-space contact that it returns 2g₄": the mechanism and (P) are shown on a borderline
  position-space Witten term (k = 2 in the product family). The literal "2g₄" needs the R⁴ ⟨TTTT⟩ itself and (F): still open.
* "Show that the same definition is used for heavy blocks": yes, by definition (§1.1).

---

## 6. What remains

1. **The R⁴ ⟨TTTT⟩ arc in Stage 0 normalisation** [D-sketch → DV].
   * Needed: the spinning Mellin/Regge form of the AdS₄ R⁴ contact and the k = 4 kernel's Regge coefficient.
   * Equivalently, the spinning analogue of CHMRS's `C_{k,ν}[M_{a,b}]` computation.
   * The (F) flat limit then gives 2g₄ + O(Λ⁻²).
2. **The calibration e = k − 3 − n** is empirical in the σ ⊗ T_free family [NE]. A derivation for Δ = 3 T with light-ray polarisations
   would fix where Einstein and R⁴ sit without invoking the flat dictionary.
3. **Primaries vs Mellin poles.**
   * The term-wise sum here is over Mellin residues (twist levels with descendants). Equality with the sum over s-channel primaries is a
     rearrangement.
   * It is legitimate if the double series converges absolutely. The residues decay with p ≥ 1.9 [D-sketch].
4. **Integrability of the light function near a = b** against C_ν (it grows as Δφ → 0: −5.7 at Δφ = 0.3) [I], as in M7 problem 3.
5. **H_K, H_EFT(b)** are unchanged.

---

## 7. Reproduce

```bash
cd ads4/code/resultq
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_pos.py P1 60)          # EIN whole-function reference, k = 4,2,5 (~6 min) -> cache/lsd_pos_P1_n60.pkl
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_pos.py P1 60 3)        # same, k = 3                                   -> cache/lsd_pos_P1_n60_k3.pkl
for k in 2 3 4 5; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py profile $k P1); done   # term-wise residues, ~3 min each
for t in P2 P3; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py profile 2 $t); done
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py profile 4 P2)
(ulimit -v 3000000; OMP_NUM_THREADS=1 LSD_GRID=alt python3 lsd_run.py profile 2 P1)                 # grid robustness
for k in 2 3 4 5; do python3 lsd_run.py report $k P1; done; python3 lsd_run.py report 2 P2; python3 lsd_run.py report 2 P3
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py exch 4 P1); (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_run.py exch 2 P1)
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_hairpin.py P1 2,4 0.0 40)                       # tau = 0 reproduces EIN (§3.5)
```

| file | content |
|---|---|
| `lsd_mellin.py` | Mellin amplitudes of EIN's pieces + r4, r4u, xu4; principal-branch logs on C; MB point evaluator; graded contour |
| `lsd_core.py` | `Setup` (H matrices, separable kernel K), `omega_line`, `omega_ope`, `residues_d12`; diagnostics `residues_right`, `closed_fast` (right closure diverges, not used) |
| `lsd_pos.py` | EIN whole-function reference values (calls `ein_run.Etable`) |
| `lsd_run.py` | `profile` (term-wise residues), `report` (sums, tails, sym, vs EIN), `exch` |
| `lsd_hairpin.py` | C_τ contour family on EIN's engine (§3.5; τ = 0 check only) |
| `cache/lsd_report_k*_P*.txt`, `lsd_exch_k*_P1.txt` | the numbers quoted in §3 |

Conventions: angle pairs/polarisations from `lsd_pos.config` (P1 (0.7, −0.5), P2 (−2.3, 1.8), P3 (0.25, −0.05)); weight t⁵; f_k of
og_free; Mellin contours (Re δ₁₂, Re δ₁₃) = (2/3, 2/3) contacts, (1.3, 0.35) xs, (0.35, 0.35) xt, (0.35, 1.3) xu/xu4/r4u.
