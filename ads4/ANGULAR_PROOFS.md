# ANGULAR_PROOFS: symbolic-J collinear identities, absolute kernels, and the two-scale crossing lemma

*26 Nov 2026. Research agent, about 4 h. Unreviewed and not committed. No existing file was edited. BATCH4_FOLD.md,
FREE_NULLPLANE.md, `code/b4fold/`, `code/freenp/` and `prl/` were not touched.*

* **Purpose.** This file writes out the two items that `ANGULAR_AUDIT.md` §0 asked for:
  * (A) the symbolic-J derivation of the three collinear identities and of the absolute detector kernels, with their envelopes;
  * (B) the two-scale crossing lemma, stated with its operator ordering, and the corollary chain down to 𝒟 = 0.
* **Code.** It is in `ads4/code/angproofs/` (map in §8). Outputs go to `angproofs/cache/`, which ignores everything but its own
  `.gitignore`. `code/kernels/{structures,pipeline,kernels_lib,canonical,closed_forms}.py` are imported read-only, with
  `KERNELS_CACHE` redirected to `angproofs/cache`.
* **Compute.** About 15 CPU-min in total, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

**Grades.**

| grade | meaning |
|---|---|
| [P] | proved: hand proof, or exact computer algebra with **symbolic** J and Δ |
| [P\|X] | proved, given the named input X |
| [DV] | verified exactly at finitely many J |
| [D] | derived, not machine-checked |
| [S] | taken from the literature |

---

## 0. Verdict

**Both items are now written as proofs. 𝒟(c_T) = 0 at finite c_T moves from [D] (audit) to [P | standard axioms + four named inputs].**

### What is now proved for all J

| item | before (audit) | now |
|---|---|---|
| canonical conserved structures S^{++}, S^{+−}, S^{odd} (even J ≥ 4) and S^{oddJ} (odd J ≥ 5), as explicit polynomial coefficient vectors in (s, J) | [C] for each J ≤ 30 separately | **[P]** conservation and Bose symmetry for symbolic J (§1) |
| collinear identities F₋₂,₋₂, F₋₁,₋₁, F₋₂,₋₁, including the J-independent factor 1/16 (∓i/16 odd) | [DV], J ≤ 30 | **[P \| (St)]** for all J (§2) |
| extreme detector amplitudes M_{∓J} of the four canonical structures, and the vanishing of the absent sector | — | **[P]** symbolic J, by Chu–Vandermonde (§3) |
| absolute kernels \|K^{++}₀\|², \|K^{odd}₀\|², \|K^{+−}₄\|², \|K^{oddJ}₄\|² (KERNELS_D3 §5.3) | [DV], J ≤ 30 | **[P \| (L)]** for all J (§3.4) |
| envelopes \|K_even\|² ≤ c·4^Δ2^{−J}√(J+1)√(sh)Δ^{−15/2}, \|K_odd\|² ≤ c·4^Δ2^{−J}√(J+1)Δ^{−15/2}/√(sh), for all J ≥ 4 and h ≥ 5 | [DV] numerics | **[P]**, explicit constants, Wendel's inequality only (§4) |
| block lower bound k^{(δ)}_{2h}(1−S⁻²) ≥ c(δ)4^hS^{1/2+2δ} on h ∈ [S,2S] | [DV] numerics | **[P]**, explicit c(δ), from Euler's integral (§5.3) |
| two-scale crossing lemma Σ4^Δ2^{−J}\|F_{m₁m₂}\|² ≤ 𝒦S^{11/2+M}U^{11/2−M}, in the mirror ordering ⟨A B B†A†⟩ | [D] | **[P \| axioms (U), (OPE), (T)]** (§5) |
| shell bound, ΣΔ²\|a₀\|² < ∞, μ((L,∞)) = O(L⁻³), O(L⁻⁴) in the band, B_N = O(1/N), 𝒟 = 0 | [P] given the above | **[P \| all of the above]** (§6) |
| low-spin sectors J = 0, 2 (F grows like Δ⁰, Δ (odd J=0), Δ³ (even J=2), Δ² (odd J=2)) | [DV] | **[P \| KERNELS_D3 J = 0, 2 data]**, exact in Δ (§6.5) |

### The named inputs (nothing else is assumed)

* **(St) Structure count [S].** For even J ≥ 4 the conserved, Bose-symmetric ⟨TTO⟩ structures form a space of dimension
  2 (even) + 1 (odd); for odd J ≥ 5 the dimension is 0 + 1 (1708.05718 §2.1; Kravchuk–Simmons-Duffin counting).
  * Here I prove that the explicit vectors are conserved and symmetric for every J. I do not prove that nothing else is.
  * The repo computes the nullspace dimension exactly for every J ≤ 30 [C].
* **(L) The helicity ladder (LADDER_ALLJ).** The normalised ratios W⁰_J = |a_{−J}/a₀|² and W⁴_{J−4} = |a_J/a₄|². They are
  [P] given the light-transform conservation identity D_z·[K, E] = 0 ([D], authors), together with the all-m norm formula
  N_{O,m} of KERNELS_D3 §3.3.
  * That norm formula is symbolic only for J ≤ 5, and [DV] for J ≤ 30 (audit, independent `norm_poly`).
  * Its extreme case m = ∓J is proved here for all J (§3.3).
* **(U), (OPE), (T)** are the standard CFT axioms, stated in §5.1: unitarity and reflection positivity; OPE convergence
  in radial quantisation with s–t crossing; a unique stress tensor with finite c_T.
* **§§13–16 of ANGULAR_ROUTE** (endpoint rule, comparable-twist ladder bound, R-tail criterion), as graded [P] there.
  * They rest on the retained physical inputs listed in ANGULAR_AUDIT §2: the physical ladder and first-moment matching,
    the contact coefficient (Δ_φ > 1), fixed-window planar regularity, and KKSZ (4.66).
  * Those inputs are unchanged by this file.

### What remains [DV] or open

* The structure count (St), and the all-m norm formula inside (L), for J > 30. Both are standard representation theory,
  but neither is re-derived here.
* B_N = O(1/N) is sharp for this input (audit §4). §18's O(N⁻²) and "absolute angular swapping" are still **not** delivered.
  They are not needed.
* The Casimir/K⁰ step of §21 item 5 remains [D] and unnecessary (audit).

---

## 1. Canonical conserved structures for symbolic J  [P]

**Conventions.** These are those of KERNELS_D3 §3.1:
* embedding building blocks V₁ = V_{1,34}, V₃ = V_{3,41}, V₄ = V_{4,13}, H_{ij}, X_{ij}, ε_{ij};
* even monomials V₁^{m₁}V₃^{m₃}V₄^{J−k}H₁₃^{n₁₃}H₁₄^{n₁₄}H₃₄^{n₃₄}, over X₁₃^{5−s}X₁₄^{s}X₃₄^{s};
* odd monomials ε_{ij}X_{ij}·(…), over the half-shifted prefactor;
* s = (Δ+J)/2, h = (Δ−J)/2, p(h) = (h−1)(h−2)(h−3).

The repo's basis uses the 13 even monomials other than H₁₃²V₄^J (one d = 3 Gram relation), and 12 odd monomials, of types
ε₁₃ and ε₁₄.

**Candidates (heuristic step, `ap1_fit.py`).** I read the canonical coefficient vectors of `code/kernels` for
J = 4,…,12, 14, 16, 18, 20, 24, 30 (even) and J = 5, 7, 9, 11, 15, 21 (odd), and interpolated every coefficient in J.
* After the repo's normalisation, every coefficient is a **polynomial in (s, J)**: total degree ≤ 4 even, ≤ 3 odd.
* All 50 coefficients reproduce the two held-out J values exactly.
* The full tables are in `cache/ap1_out.txt`.

The entries that matter below are:

| structure | V₁²V₃²V₄^J | V₁V₃H₁₃V₄^J | normalising entry (KERNELS_D3 §5.3) |
|---|---|---|---|
| S^{++} | 2p(h)(7s+J−4) | 8(s−1)p(h) | V₁V₃V₄^JH₁₃ has leading term 8s⁴ ✓ |
| S^{+−} | p(h)(9s−J−4) | 4s p(h) | H₁₄²H₃₄²V₄^{J−4} = (128s⁴ + …)/16, leading term 8s⁴ ✓ |

| structure | ε₁₃X₁₃V₁V₃V₄^J (c_VV) | ε₁₃X₁₃H₁₃V₄^J (c_H) | normalising entry |
|---|---|---|---|
| S^{odd} (even J) | 10J³−20J²s+56J²−28Js+54J+8s³−28s²+14s+15 | (2s−J−5)(4s²−4Js−J²−12s+7J+9) | ε₁₄X₁₄V₁V₃V₄^{J−2}H₃₄ = J(J−1)(J−2) ✓ |
| S^{oddJ} (odd J) | 10J³−20J²s+80J²−140Js+190J+8s³+36s²−194s+135 | (2s−J−1)(4s²−4Js−J²−4s−9J−15) | ε₁₃X₁₃V₁V₃V₄^J has leading term 8s³ ✓ |

**Theorem 1 [P].** For every J (and symbolic s), each of the four vectors defines a ⟨O T T⟩ structure that is:
* conserved at the first T;
* symmetric under (x₁,z₁) ↔ (x₃,z₃).

*Proof (`ap2_conservation.py`).*
1. **Reduction to one configuration.** Δ_T = 3 = J_T + d − 2, so the divergence ∂_{x₁}·D_{z₁}S of a conformally covariant
   structure is again covariant (Todorov operator D, d = 3). The same holds for S(1,3) − S(3,1).
   * The conformal group acts on triples of mutually spacelike points with a one-dimensional stabiliser, so the orbit of
     one triple is open in the 9-dimensional configuration space.
   * Both objects are real-analytic, so each vanishes identically if and only if it vanishes at one configuration for all
     polarisations.
2. **The configuration.** I take x₄ = 0, x₃ = (0,0,1), x₁ = (0,0,−2).
   * x₁ is kept symbolic until after differentiation.
   * z₁ is fully symbolic and off the cone; the Todorov operator is interior.
   * z₃ = ζ(a₃,1) and z₄ = ζ(a₄,1), with a₃, a₄ symbolic. These are generic null polarisations up to scale.
3. **The symbolic exponents.** The prefactor V₄^J X₁₃^{−e₁₃}X₁₄^{−e₁₄}X₃₄^{−e₃₄} has symbolic exponents. It is differentiated
   logarithmically (V₄ depends on x₁ but not on z₁); every remaining exponent is an integer.
4. **The check.** The three z₁-components of the divergence, and the symmetry difference, are rational functions of
   (a₃, a₄, s, J) [and a₁]. Their numerators expand to **0** identically for all four structures (10–13 s each).
5. **Negative controls.** The test fails, as it should, for:
   * S^{++} with one coefficient shifted by 1;
   * S^{++} with J → J+2 in its coefficients;
   * the single monomial V₁²V₃²V₄^J.

   (`cache/ap2_controls.txt`.) ∎

**Identification with the canonical structures [P | (St), (L)].**
* With (St), the even conserved space is spanned by the two even vectors, and the odd space by the odd vector.
* §3 proves, for symbolic J, that M_{+J}[S^{++}] = 0 while M_{+J}[S^{+−}] ≠ 0, and that M_{−J}[S^{+−}] = 0.
* By the ladder (L), M_{±J} = 0 is equivalent to the vanishing of the whole b = 4 (resp. b = 0) sector. So the two vectors
  are exactly the KERNELS_D3 structures defined by "P_{m≥4} = 0" and "P_{m≤0} = 0".
* The normalisation entries agree with KERNELS_D3 §5.3 (tables above).

---

## 2. The collinear identities for symbolic J  [P | (St)]

**Definition** (ANGULAR_AUDIT 1a; `ca1_collinear.py`).
* Put the three points collinear in the Euclidean plane t = 0: x₁ = (0,−1,0), x₃ = (0,1,0), x₄ = 0.
* O is polarised by e₊ = (0,1,i). This is its top planar-spin component, a 2d quasi-primary of weights (h, s) in (w, w̄).
* The T polarisations are z_i = ζ(a_i,b_i). The coefficient of a^{2+m}b^{2−m} is the planar-spin-m component T_m.
* F_{m₁m₃} := 2^{6−Δ} × [a₁^{2+m₁}b₁^{2−m₁}a₃^{2+m₃}b₃^{2−m₃}]⟨O T T⟩.
* The 2d three-point structure has modulus 2^{Δ−6} here (|w₁₃| = 2, |w₁₄| = |w₃₄| = 1), so F is the 2d OPE coefficient up
  to a phase.

**Lemma 2.1 (building blocks at the collinear point) [P].** Direct evaluation (`ap3_collinear.py`) gives:

| block | value |
|---|---|
| V₁ | −(a₁²+b₁²) |
| V₃ | −(a₃²+b₃²) |
| **V₄** | **−½** |
| H₁₃ | −2(a₁a₃+b₁b₃)² |
| H₁₄ | −a₁² |
| H₃₄ | −a₃² |
| ε₁₃ | i(a₁a₃+b₁b₃)(a₁b₃−a₃b₁) |
| ε₁₄ | −2i a₁b₁ |
| ε₃₄ | −2i a₃b₃ |
| X₁₃ | 4 |
| X₁₄, X₃₄ | 1 |

Hence the prefactor 2^{6−Δ}V₄^J X₁₃^{−e₁₃}X₁₄^{−e₁₄}X₃₄^{−e₃₄} equals 2^{−4−2η}(2V₄)^J = 2^{−4−2η}(−1)^J, with η = 0 (even) or ½ (odd).
The J-dependence enters **only** through the polynomial coefficients and the sign (−1)^J. This is why the normalisation
factor is J-independent.

**Collapse of the monomial sums.**
* For a general component, F_{m₁m₃} is a finite sum over the monomial labels. Each term contains a binomial expansion of
  (a₁²+b₁²)^{m₁}(a₃²+b₃²)^{m₃}(a₁a₃+b₁b₃)^{2n₁₃}, times a₁^{2n₁₄}a₃^{2n₃₄}.
* At the three components used here, every term with n₁₄ > 0 or n₃₄ > 0 drops out (H₁₄ ∝ a₁², H₃₄ ∝ a₃² carry too many a's),
  and so do the ε₁₄ and ε₃₄ terms (∝ a₁b₁, a₃b₃).
* Only the V₄^J monomials survive. So the binomial sums collapse to a few terms, and no summation identity is needed in
  this frame. (The binomial sums do appear in the detector amplitudes, §3.)
* By counting a/b parity, even structures contribute only to components with m₁ + m₃ even, and odd structures only to
  those with m₁ + m₃ odd.

Reading off the coefficients:

* **b₁⁴b₃⁴.** It comes from V₁²V₃² (coefficient +1) and from V₁V₃H₁₃ (coefficient −2). So, for even J,
  F₋₂,₋₂ = (1/16)[c(V₁²V₃²V₄^J) − 2c(V₁V₃H₁₃V₄^J)].
  * S^{++}: 2p(h)(7s+J−4) − 16(s−1)p(h) = −2p(h)(h−4).
  * S^{+−}: p(h)(9s−J−4) − 8sp(h) = p(h)(h−4).
* **a₁b₁³a₃b₃³.** It comes only from V₁V₃H₁₃, with coefficient −4. So F₋₁,₋₁ = −(1/4)c(V₁V₃H₁₃V₄^J).
  * S^{++}: −2(s−1)p(h).
  * S^{+−}: −s p(h).
* **b₁⁴a₃b₃³ (odd).** Only ε₁₃X₁₃V₁V₃V₄^J (weight −i/8) and ε₁₃X₁₃H₁₃V₄^J (weight +i/4) contribute. So
  F₋₂,₋₁ = (−1)^J(i/8)(2c_H − c_VV).
  * In both odd families, 2c_H − c_VV = (2h−7)(2h−5)(2h−3), identically in (s, J). This is a polynomial identity, checked
    by expansion.

**Theorem 2 [P | (St)].** For all even J ≥ 4, all odd J ≥ 5 and all s (with h = s − J):

* 16·F₋₂,₋₂ = p(h)(h−4)(λ₊₋ − 2λ₊₊)
* 16·F₋₁,₋₁ = −16p(h)[2(s−1)λ₊₊ + sλ₊₋]
* (16i)·F₋₂,₋₁ = −2(2h−7)(2h−5)(2h−3)λ_odd (even J), and the same with (−16i) for odd J (λ_odd ≡ λ_oddJ).

Moreover:
* no odd structure enters F₋₂,₋₂ or F₋₁,₋₁;
* no even structure enters F₋₂,₋₁.

So ANGULAR_ROUTE §21 item 1 holds exactly, with the common factor 1/16 (∓i/16 odd) found by the audit, for every J.
(`cache/ap3_out.txt` shows the symbolic ratios 1/16, 1/16, −i/16 and i/16 = (−1)^J(−i/16).)

---

## 3. Extreme detector amplitudes and absolute kernels for symbolic J

### 3.1 Method  [P]  (`ap4_extreme.py`)

**Chain.** I use the chain (F1)–(F3) of KERNELS_D3 §3.2:
* source T(e₊), detector direction θ = 0, rest frame p = (1,0,0);
* M_m(θ) = ∫d³x e^{−ip·x}⟨0|O_m(x)†E(n)T(0)|0⟩.

**Extreme components.** The component m = −J (resp. +J) is the coefficient of t^{2J} (resp. t⁰) for the bra ζ(t,1), i.e.
the structure evaluated at bra polarisation ζ = e₊/2 (resp. e₋/2). In the detector limit:
* V₄ = ((ζ·x)P₄ − (z·ζ)x²)/P₃, with z·ζ = ½;
* the other blocks are as in KERNELS_D3 (F1).

For a monomial with V₄^n, n = J − k, the binomial expansion gives

  V₄^n = Σ_{j=0}^n C(n,j)(−½)^{n−j} P₃^{−n}P₄^j (x²)^{n−j}(ζ·x)^j.

Each term then goes through the α-integral (F2) and the Fourier transform (F3), with (ζ·x)^{n₂+j} → (iζ·∂_p)^{n₂+j}.

**Three facts make every term a hypergeometric term in j.**
1. **Null direction.** ζ is null and ζ·p = 0 at rest, so ζ·∂_p acts only on (−4z·p)^{−A}. The bounded number of other
   derivative directions (z₃, an ε-form) contribute a finite generating-function sum, (1−X)^g(1−M)^{−A} with
   X = O(τ_small). The j-dependence therefore enters through (A)_{N−r+R} with N = n₂ + j and A = A₀ − j, i.e. through
   Γ(A₀ + N₀ + …)/Γ(A₀ − j).
2. **No j-dependent ratio.** The powers (−z·ζ)^{−j} e^{iπj/2} 4^{−j} i^j 4^j (L_u)^j combine to (L_u/z·ζ)^j = 1, since
   L_u = z·ζ.
3. **Cancellation of Γ(A).** Γ(a+b−1)/Γ(A) = (n−j+c−l)_l is a polynomial in j, with c = esum−1−i_T−j_T.

So the summand is C(n,j)·Q(j)/(Γ(x−j)Γ(y+j)), with x = e₁₄ − j_T, y = e₁₄ − k_T − n and Q a polynomial. Expanding Q in
falling factorials and applying **Chu–Vandermonde** in Γ-form,

  Σ_j C(n,j) j^{(r)}/(Γ(x−j)Γ(y+j)) = n^{(r)}Γ(x+y+n−r−1)/(Γ(x+y−1)Γ(x−r)Γ(y+n)),

gives a closed form for every monomial. (The identity is a terminating Gauss 2F1(−m,·;·;1); it is proved by Vandermonde's
convolution of generalised binomials.)

**Checks of the implementation.** It reproduces the read-only pipeline (independent flint code, exact) at J = 4 and 14
(15 for odd J), at two non-integer s each, to all printed digits (`cache/ap4_out.txt`, `cache/ap4_pp_out.txt`).

### 3.2 Result  [P]  (`ap5_kernels.py`)

**Reduction.** The sums are reduced exactly for J = 2L (resp. 2L+1), with L symbolic:
* all Γ's are reduced to four base Γ's by integer shifts and the Legendre duplication formula;
* powers in L are reduced to 2^L and (−1)^L;
* the remainder is cancelled as a rational function of (s, L).

The results (`cache/ap5_out.txt`) are:

| structure | M_{−J} | M_{+J} |
|---|---|---|
| S^{++} | 4i(−1)^{J/2}e^{iπs}π³2^{−J}(J+3)! / (Γ(s)Γ(Δ+3)Γ(1−h)) | **0** |
| S^{odd} | 8(−1)^{J/2}e^{iπs}π³2^{−J}(J+3)! / (Γ(s+½)Γ(Δ+3)Γ(3/2−h)) | **0** |
| S^{+−} | **0** | −2i(−1)^{J/2}e^{iπs}π³2^{−J}(J−1)! / (Γ(s)Γ(Δ−1)Γ(1−h)) |
| S^{oddJ} | **0** | −8i(−1)^{(J−5)/2}e^{iπs}π³2^{−J}(J−1)! / (Γ(s+½)Γ(Δ−1)Γ(3/2−h)) |

* "(amplitude)/(claimed closed form)" simplifies to exactly **1** in each nonzero case.
* "(amplitude)/(first term)" simplifies to **0** in each absent sector.

The nonzero entries are the KERNELS_D3 §5.3 closed forms for M₀ and M₄, times the ladder factors (4)_J/(Δ−1−J)_J and
(4)_{J−4}/(Δ−1−J)_{J−4}.

### 3.3 Extreme norm  [P]  (`ap9_norm_kernel.py`)

Take the bra e₊/2 and the ket e₋/2 in ⟨O(x,ζ′)O(0,ζ)⟩ = (ζ′·ζ x² − 2ζ′·xζ·x)^J/(x²)^{Δ+J}.
* Both are null with ζ·p = 0, so the only derivative contractions are ζ′·ζ = ½.
* The Fourier transform gives π^{5/2}2^{4−2Δ−J}Γ(Δ−½)^{−1}Σ_k C(J,k)(−1)^k k!/Γ(Δ+k).
* By Gauss, this sum is 2F1(−J,1;Δ;1)/Γ(Δ) = (Δ−1)_J/Γ(Δ+J).

So

  N_{O,∓J} = π^{5/2}2^{4−2Δ−J}(Δ−1)_J/(Γ(Δ+J)Γ(Δ−½)),

which is the KERNELS_D3 §3.3 formula at m = ∓J.

### 3.4 Absolute kernels  [P | (L)]

The extreme kernels are |K_{∓J}|² = 12|M_{∓J}|²/N_{∓J}. The ladder (L) gives the normalised ratios
* W⁰_J = C(J+3,3)²(1)_{2J}/(Δ−1−J)_{2J} = |a_{−J}/a₀|², and
* W⁴_{J−4} = C(J−1,3)²(1)_{J−4}(J+5)_{J−4}/((Δ−1−J)_{J−4}(Δ+3)_{J−4}) = |a_J/a₄|².

Then |K₀|² = |K_{−J}|²/W⁰_J and |K₄|² = |K_J|²/W⁴_{J−4}. With symbolic J and Δ, these reproduce **exactly** the four
KERNELS_D3 §5.3 closed forms (`gammasimp` ratio = 1). After the duplication formula they take the factorised form

* |K^{++}₀|² = 864π^{3/2} sin²(πh) 4^Δ [2^J/C(2J,J)] [Γ(s+½)/Γ(s)] [Γ(h)/Γ(h−½)] Γ(Δ−1)Γ(Δ−½)/Γ(Δ+3)²
* |K^{odd}₀|² = 3456π^{3/2} cos²(πh) 4^Δ [2^J/C(2J,J)] [Γ(s)/Γ(s+½)] [Γ(h−½)/Γ(h)] Γ(Δ−1)Γ(Δ−½)/Γ(Δ+3)²
* |K^{+−}₄|² = 216π^{3/2} sin²(πh) 4^Δ [2^J/C(2J,J+4)] [Γ(s+½)/Γ(s)] [Γ(h)/Γ(h−½)] Γ(Δ−5)Γ(Δ−½)/(Γ(Δ+3)Γ(Δ−1))
* |K^{oddJ}₄|² = 3456π^{3/2} cos²(πh) 4^Δ [2^J/C(2J,J+4)] [Γ(s)/Γ(s+½)] [Γ(h−½)/Γ(h)] Γ(Δ−5)Γ(Δ−½)/(Γ(Δ+3)Γ(Δ−1))

The factorised forms agree with `closed_forms.py` to 10⁻³⁸ (`ap6_envelopes.py`).
* sin²(πh) comes from Γ(1−h)⁻² = sin²(πh)Γ(h)²/π², and cos²(πh) from Γ(3/2−h)⁻².
* The double zeros of K² at h = 1, 2, 3, 4 (even) and h = 3/2, 5/2, 7/2 (odd) are the zeros of the collinear F of §2.

---

## 4. Envelopes with explicit constants  [P]  (`ap6_envelopes.py`)

**Wendel's inequality (1948).** For x > 0 and 0 < a < 1, (x/(x+a))^{1−a} ≤ Γ(x+a)/(x^aΓ(x)) ≤ 1. With a = ½ this gives:
* (W1) Γ(x+½)/Γ(x) ≤ √x;
* (W2) Γ(x+½)/Γ(x) ≥ x/√(x+½).

These give, for every J ≥ 0 and h > ½ (no asymptotics; this replaces Stirling with remainder):
1. 2^J/C(2J,J) = 2^{−J}√π Γ(J+1)/Γ(J+½) ≤ 2^{−J}√π √(J+½), by W1 at x = J+½.
2. Γ(s+½)/Γ(s) ≤ √s and Γ(h)/Γ(h−½) ≤ √h, by W1.
3. Γ(s)/Γ(s+½) ≤ s^{−1/2}(1+1/(2s))^{1/2} and Γ(h−½)/Γ(h) ≤ h^{−1/2}·h/(h−½), by W2.
4. Γ(Δ−½)/Γ(Δ+3) ≤ √(Δ+3)/((Δ+5/2)(Δ−½)(Δ+½)(Δ+3/2)), by W2 at x = Δ+5/2. Then √(Δ(Δ+3)) ≤ Δ+5/2, and every factor
   Δ/(Δ+k) with k ≥ 0 is at most 1. Hence:
   * Δ^{15/2}Γ(Δ−1)Γ(Δ−½)/Γ(Δ+3)² ≤ g₀(Δ) := (Δ/(Δ−1))(Δ/(Δ−½));
   * Δ^{15/2}Γ(Δ−5)Γ(Δ−½)/(Γ(Δ+3)Γ(Δ−1)) ≤ g₄(Δ) := Δ⁴/((Δ−5)(Δ−4)(Δ−3)(Δ−2))·Δ/(Δ−½).

   Both are decreasing in Δ.
5. sin², cos² ≤ 1, and C(2J,J)/C(2J,J+4) = ρ_J := (J+1)(J+2)(J+3)(J+4)/((J−3)(J−2)(J−1)J), which decreases from ρ₄ = 70 to 1.

**Theorem 4 [P].** For all J ≥ 4 and h ≥ 5 (so Δ ≥ 14 and s ≥ 9):
* |K^{++}₀|² ≤ c₊₊ · 4^Δ2^{−J}√(J+1)√(sh)Δ^{−15/2},
* |K^{odd}₀|² ≤ c_odd · 4^Δ2^{−J}√(J+1)Δ^{−15/2}/√(sh),
* |K^{+−}₄|² ≤ c₊₋(J) · 4^Δ2^{−J}√(J+1)√(sh)Δ^{−15/2},
* |K^{oddJ}₄|² ≤ c_oddJ(J) · 4^Δ2^{−J}√(J+1)Δ^{−15/2}/√(sh).

The constants are:

| constant | formula | value |
|---|---|---|
| c₊₊ | 864π²g₀(14), with g₀(14) = (14/13)(28/27) = 1.1168 | 9.53·10³ |
| c_odd | 3456π²(19/18)^{1/2}(10/9)g₀(14) | 4.35·10⁴ |
| c₊₋(J) | 216π²ρ_J g₄(14), with g₄(14) = 3.3534 | ≤ 5.01·10⁵ (J = 4), → 7.1·10³ as J → ∞ |
| c_oddJ(J) | 3456π²(19/18)^{1/2}(10/9)ρ_J g₄(14) | ≤ 3.3·10⁶ (J = 5), → 1.3·10⁵ |

*Proof.* Insert items 1–5 into the factorised forms of §3.4, using √(J+½) ≤ √(J+1) and π^{3/2}·√π = π². ∎

**Scan.** On J ≤ 10⁴ and h ∈ [5, 10⁴], the largest ratios |K|²/envelope are 7.9·10³, 3.7·10⁴, 3.2·10⁵ and 2.2·10⁶
(`cache/ap6_out.txt`), all below the proven constants. With sin² = ½ the audit's ranges are reproduced.

---

## 5. The two-scale crossing lemma  (`ap8_blocks.py`)

### 5.1 Setting and assumptions

**Theory.** Euclidean 3d CFT, the plane Π = {t = 0}, with complex coordinate w = x¹ + ix².
* The global conformal algebra of Π, sl(2)×sl(2) (L₀, L_{±1} built from D, M₁₂, P_w, K_w and their bars), is a
  subalgebra of so(4,1).
* In radial quantisation about a point of Π, L_{−1}† = L₁.
* Every 3d primary component restricted to Π is annihilated by L₁ and L̄₁, so it is a 2d quasi-primary.
* O(e₊) has weights (h_p, h̄_p) = (h, s) = ((Δ−J)/2, (Δ+J)/2). T_m has weights ((3−m)/2, (3+m)/2): T₋₂ = (5/2, ½), T₋₁ = (2, 1).

**Assumptions.**
* **(U) Unitarity / reflection positivity.** The radial-quantisation Hilbert space is positive.
  * Consequence: every state other than the vacuum has h, h̄ ≥ ¼. Write E = h + h̄ and m = h − h̄ for energy and planar spin.
    Unitarity of so(3,2) gives E − j ≥ ½ for every state of spin j ≥ |m| (scalars Δ ≥ ½, spinors Δ ≥ 1, J ≥ 1: Δ ≥ J+1,
    and descendants raise E at least as fast as j). Hence h, h̄ = (E ± m)/2 ≥ ¼.
* **(OPE) Convergence and crossing.** The OPEs A×B, B×B† and A×A† converge in radial quantisation at every Euclidean
  configuration. The s- and t-channel expansions of a four-point function agree where both converge (associativity).
* **(T) Stress tensor.** T is the unique spin-2, Δ = 3 primary, and 0 < c_T < ∞.
  * So the components T_m have finite, nonzero two-point norms N_m.
  * ⟨O T T⟩ = Σ_a λ_a S_a with the canonical structures (with (St)).
  * The 2d OPE coefficient of T_{m₁}×T_{m₂} → O(e₊) is F_{m₁m₂} up to a phase and a J-independent constant κ₀ (§2).
  * O(e₊) has norm (e₊*·e₊)^J = 2^J.

No assumption is made on the spectrum beyond (U).

### 5.2 Positivity at independent real z, z̄  [P | (U), (OPE)]

**Lemma B1.** Let A = T_{m₁}, B = T_{m₂}, with weights (h_A, h̄_A) and (h_B, h̄_B). In the frame (0, z, 1, ∞), take the
**mirror ordering**

  G(z,z̄) := lim_{w₄→∞}|w₄|^{…}⟨A†(w₄) B†(1) B(z) A(0)⟩,   i.e. ⟨A B B†A†⟩.

For independent real z, z̄ ∈ (0,1):

  G(z,z̄) = z^{−h_A−h_B}z̄^{−h̄_A−h̄_B} Σ_p (|C_{ABp}|²/N_p) k^{(δ)}_{2h_p}(z) k^{(δ̄)}_{2h̄_p}(z̄),

with:
* k^{(δ)}_{2h}(z) = z^h 2F1(h+δ, h+δ; 2h; z);
* δ = h_B − h_A, δ̄ = h̄_B − h̄_A;
* p running over an orthonormal basis of sl(2)×sl(2) quasi-primaries;
* every term ≥ 0 and the series absolutely convergent.

*Proof.*
1. **The expansion.** Insert the sl(2)×sl(2) decomposition of the Hilbert space into ⟨Ψ'|Ψ⟩, with Ψ = B(z)A(0)|0⟩. The
   coefficient of the p-block is C_{ABp}C_{B†A†p†}/N_p.
2. **Positivity of the coefficients.** Reflection positivity (U) identifies C_{B†A†p†} = conj(C_{ABp}).
3. **Positivity of the blocks.** In standard notation the 2F1 parameters are (h_p − h₁₂, h_p + h₃₄), with
   h₁₂ = h_A − h_B and h₃₄ = h_B − h_A. Both equal h_p + δ, so the Taylor coefficients ((h_p+δ)_n)²/((2h_p)_n n!) are ≥ 0, and
   every block is positive on (0,1).
4. **Convergence.** The expansion is a double series in z^{h_p+n}z̄^{h̄_p+n'} with nonnegative coefficients. At z = z̄ = r < 1
   it is a convergent Euclidean correlator (OPE). For max(z, z̄) ≤ r it is dominated termwise.
5. **Equality with G.** The series converges absolutely in the cut bidisc, so it defines a function holomorphic in (z, z̄)
   there. It coincides with G on the Euclidean section z̄ = z*, which is totally real of maximal dimension, and hence on the
   connected domain containing (0,1)². ∎

*Remark (ordering).* With ⟨A B A†B†⟩ the parameters become (h_p+δ, h_p−δ) and the coefficients C_{ABp}C_{A†B†p}. Neither is
sign-definite, so the argument fails, as the audit noted. The mirror ordering is essential.

### 5.3 The block lower bound  [P]

**Lemma B2.** For |δ| ≤ 1, S ≥ 4 and h ∈ [S, 2S]:

  k^{(δ)}_{2h}(1−S⁻²) ≥ c(δ) 4^h S^{1/2+2δ},

  c(δ) = (1−2/S₀)(1−1/(S₀−1))(1+1/(2S₀))^{−1/2}(2√π)^{−1}(1−3/S₀)^{2S₀} 2^{−δ−1}3^{−max(δ,0)}, at S₀ = 4

(each factor increases with S₀). Numerically:
* c(0) = 6.76·10⁻⁷, c(½) = 2.76·10⁻⁷, c(−½) = 9.57·10⁻⁷;
* at S₀ = 8 or 16, c(0) = 4.8·10⁻⁵ or 1.5·10⁻⁴.

*Proof.* Let z = 1 − ε with ε = S⁻². Since 2h > h + δ > 0, Euler's integral applies:

  2F1(h+δ, h+δ; 2h; z) = Γ(2h)/(Γ(h+δ)Γ(h−δ)) ∫₀¹ t^{h+δ−1}(1−t)^{h−δ−1}(1−zt)^{−h−δ} dt,

with a positive integrand. Restrict to u := 1 − t ∈ [1/S, 2/S]. There:
1. 1 − zt = u + εt ∈ [1/S, 3/S], and t ≥ 1 − 2/S.
2. tu/(1−zt) ≥ (1−2/S)·(1/S)/(1/S + 1/S²) ≥ 1 − 3/S, because u/(u+ε) increases in u. So [tu/(1−zt)]^h ≥ (1−3/S)^{2S}.
3. t^{δ−1} ≥ 1 (δ ≤ 1).
4. u^{−δ−1} ≥ (2/S)^{−δ−1} (δ ≥ −1).
5. (1−zt)^{−δ} ≥ 3^{−max(δ,0)}S^{δ}.

Integrating over an interval of length 1/S gives ∫ ≥ (1−3/S)^{2S}2^{−δ−1}3^{−max(δ,0)}S^{2δ}. For the prefactor:
* Γ(2h)/Γ(h)² = 4^hΓ(h+½)/(2√πΓ(h)) ≥ 4^h√h(1+1/(2h))^{−1/2}/(2√π), by W2.
* By Gauss's product, Γ(h)²/(Γ(h+δ)Γ(h−δ)) = Π_{n≥0}(1 − δ²/(h+n)²) ≥ 1 − δ²Σ_{n≥0}(h+n)^{−2} ≥ 1 − 1/(h−1).
* Finally z^h ≥ (1−S⁻²)^{2S} ≥ 1 − 2/S, and √h ≥ √S. ∎

Numerically (`cache/ap8_out.txt`), the true minimum of the ratio over S ∈ [4, 4096] is ≈ 0.009–0.014 for
δ ∈ {0, ±½, ±1}, consistent with the bound and with the audit's [0.009, 1.24].

### 5.4 The crossed channel  [P | (U), (OPE)]

**Lemma B3.** For u := 1−z ∈ (0, ¼] and ū := 1−z̄ ∈ (0, ¼]:

  0 < G(z,z̄) ≤ 𝒞_{AB} u^{−2h_B} ū^{−2h̄_B},   𝒞_{AB} = N_AN_B + (M_A M_B)^{1/2},

where:
* M_B := ¼^{−2h_B}¼^{−2h̄_B}G_{B B† B B†}(¼,¼) − N_B²;
* M_A := ¼^{−2h_A}¼^{−2h̄_A}G_{A†A A†A}(¼,¼) − N_A²;
* the four-point functions are taken in mirror ordering at the Euclidean point ¼.

*Proof.*
1. **The t-channel.** Crossing (x ↦ 1−x, a rotation by π in Π, which gives a phase only) gives
   |G(z,z̄)| = |Σ_q C_{B†Bq}C_{AA†q}N_q^{−1} u^{h_q−2h_B}ū^{h̄_q−2h̄_B} b_q(u)b̄_q(ū)|, with b_q(u) = 2F1(h_q,h_q;2h_q;u) > 0.
   The 2F1 parameters are equal because the external pairs (B,B†) and (A,A†) have equal weights.
2. **The identity term.** It gives N_AN_B u^{−2h_B}ū^{−2h̄_B}.
3. **Cauchy–Schwarz.** For the remaining terms, Σ|C_{B†Bq}||C_{AA†q}|w_q ≤ (Σ|C_{B†Bq}|²w_q)^{1/2}(Σ|C_{AA†q}|²w_q)^{1/2},
   with w_q = N_q^{−1}u^{h_q}ū^{h̄_q}b_q b̄_q ≥ 0.
4. **The two positive sums.** They are, by Lemma B1 applied with δ = 0, the s-channel series of the mirror-ordered
   ⟨B B† B B†⟩ and ⟨A† A A† A⟩, with the identity removed and the prefactor stripped.
5. **Monotonicity.** Each is a double series with nonnegative coefficients and exponents ≥ ¼ (U). It therefore increases in
   u and in ū separately, and on (0,¼]² it is at most its value at the Euclidean point u = ū = ¼, which is finite by (OPE).
   The two variables are independent here; this is the "two-scale" step. ∎

### 5.5 The lemma

**Theorem B (two-scale crossing) [P | (U), (OPE), (T)].** Let (m₁, m₂) ∈ {(−2,−2), (−1,−1), (−2,−1)} and
M = m₁ + m₂. For all S, U ≥ 4:

  Σ_{O: s∈[S,2S], h∈[U,2U]} 4^{Δ} 2^{−J} |F^O_{m₁m₂}|² ≤ 𝒦_{m₁m₂} S^{11/2+M} U^{11/2−M},

  𝒦_{m₁m₂} = 𝒞_{AB}/(κ₀ c(δ)c(δ̄)),   δ = −δ̄ = (m₁−m₂)/2.

* The sum runs over an orthonormal basis of 3d primaries, so degeneracies are included.
* 𝒦 depends on the CFT only through N_{T_m}, the two values M_{T_m} (finite-c_T data) and κ₀.
* F^O is defined from the actual ⟨O T T⟩, so the theorem does not use (St). (St) enters only when F is rewritten through λ.

*Proof.*
1. **Lower bound on G.** Choose z = 1 − U⁻² and z̄ = 1 − S⁻². In Lemma B1, keep only the quasi-primaries p = O(e₊) with
   (h_p, h̄_p) = (h, s) in the shell. All dropped terms are ≥ 0, and z^{−h_A−h_B}z̄^{…} ≥ 1.
2. **Coefficients.** By (T), |C_{ABp}|²/N_p = κ₀2^{−J}|F^O_{m₁m₂}|².
3. **Blocks.** By Lemma B2, k_{2h}(z) ≥ c(δ)4^hU^{1/2+2δ} and k_{2s}(z̄) ≥ c(δ̄)4^sS^{1/2+2δ̄}.
4. **Upper bound on G.** By Lemma B3, G ≤ 𝒞_{AB}U^{4h_B}S^{4h̄_B}.
5. **Combination.** Dividing, the exponents are 4h_B − ½ − 2(h_B−h_A) = 2(h_A+h_B) − ½ = 11/2 − M (U) and
   2(h̄_A+h̄_B) − ½ = 11/2 + M (S), since h_{T_m} + h̄_{T_m} = 3 and h_{T_m} − h̄_{T_m} = −m. ∎

**Inverting the identities [P].** For h ≥ 5:
* p(h)(h−4) ≥ (4·3·2·1/5⁴)h⁴ = 0.0384h⁴;
* p(h) ≥ 0.192h³;
* (2h−7)(2h−5)(2h−3) ≥ 0.84h³.

With X₁ = λ₊₋ − 2λ₊₊ and X₂ = 2(1−1/s)λ₊₊ + λ₊₋, one has λ₊₊ = (X₂−X₁)/(4−2/s) and λ₊₋ = X₁ + 2λ₊₊, so
|λ₊₊|² + |λ₊₋|² ≤ 2.8(|X₁|² + |X₂|²) for s ≥ 9. Theorem B with M = −4, −2, −3 then gives, on each shell (S, U ≥ 4, h ≥ 5):

  Σ 4^Δ2^{−J}(|λ₊₊|²+|λ₊₋|²) ≤ 𝒦′ S^{3/2}U^{3/2},   Σ 4^Δ2^{−J}|λ_odd|² ≤ 𝒦′ S^{5/2}U^{5/2}.

Here M = −4 gives U^{19/2}/U⁸; M = −2 gives S^{7/2}U^{15/2}/(U⁶S²); M = −3 gives S^{5/2}U^{17/2}/U⁶.

---

## 6. Corollary chain  [P | Theorems 2, 4, B, inputs]

Notation (ANGULAR_ROUTE §§13–16):
* |a_{O,0,0}|² = (λ₊₊²|K^{++}₀|² + λ_odd²|K^{odd}₀|²)/C_T and |a_{O,4,0}|² = λ₊₋²|K^{+−}₄|²/C_T (odd J: λ²|K^{oddJ}₄|²/C_T);
* C = J(J+1), T = (Δ−1)(Δ−2), T − C = (2h−2)(2s−1);
* R = T/(T−C), and μ = Σ_{O,b} C|a_{O,b,0}|²δ_R.

**6.1 Shell bound.** For S ≥ U ≥ 4 and h ≥ 5: Σ_{shell}Σ_b|a_{O,b,0}|² ≤ 𝒦″S^{−5}U².

*Proof.* Combine Theorem 4 with the inverted Theorem B. On the shell:
* J+1 ≤ 5S, √(sh) ≤ 2√(SU), and Δ ≥ s ≥ S.
* Even: S^{3/2}U^{3/2}·√(5S)·2√(SU)·S^{−15/2} ∝ S^{−5}U².
* Odd: S^{5/2}U^{5/2}·√(5S)/√(SU)·S^{−15/2} ∝ S^{−5}U². ∎

**6.2 Endpoint budget.** Σ_{τ≥10, J≥4}Δ²|a₀|² ≤ Σ_{S dyadic}16S²·𝒦″S^{−5}Σ_{U≤2S}U² ≲ Σ_S S^{−1} < ∞.

**6.3 Tail.** On a shell with U ≥ 4: R ≤ (4S)²/((2U−2)(2S−1)) ≤ 6.1S/U, and C ≤ 4.5S². So R > L forces U < 6.1S/L, and

  μ((L,∞)) ≤ Σ_S 4.5S²𝒦″S^{−5}·(4/3)(6.1S/L)² ≲ Σ_{S ≳ L}S^{−1}L^{−2} ≲ L^{−3}.

In the intermediate band (τ ≥ κ√J and τ/J → 0, where J ≳ R²) only S ≳ L² contributes, which gives O(L⁻⁴). The sector
τ < 10 is bounded twist and contributes nothing (§14).

**6.4 Boundary term and defect.** By §15, B_N ≤ ((N+2)(N+3)/6)∫e^{−N/R}dμ(R) for T > C. Moreover
∫e^{−N/R}dμ = ∫₀^∞e^{−u}μ((N/u,∞))du ≤ ∫e^{−u}min(μ_tot, K u³N^{−3})du ≤ 6K N^{−3}. Hence:
* B_N ≤ K(N+2)(N+3)N^{−3} = O(1/N);
* the bounded-twist sector is o(1) by §14;
* J = 0, 2 do not enter for N > 2 (a_N = 0 for N > J − b).

So 𝒟 = lim B_N = 0, and equivalently μ((L,∞)) = o(L⁻²) (§15 R-tail criterion). B_N = O(1/N) is sharp for this input (audit §4).

**6.5 Low-spin sectors [P | KERNELS_D3 J = 0, 2].** Theorem B applies with S ≍ U ≍ Δ₀/2, so the budget on a dyadic Δ-shell
is Σ4^Δ|F|² ≲ Δ₀^{11}. The collinear components, exact in Δ (`ap7_lowspin.py`; the same evaluation as §2, on the
KERNELS_D3 J = 0, 2 structures), are:

| sector | F per unit λ | growth Δ^g | \|K\|²/4^Δ (KERNELS_D3 §5.2–5.3) | Σ_shell Δ²\|a\|² |
|---|---|---|---|---|
| J = 0 even | (Δ−8)(Δ−6)/(16Δ(Δ+2)) | g = 0 | Δ^{−29/2} | Δ₀^{2+11−29/2} = Δ₀^{−3/2} |
| J = 0 odd | i(Δ−7)/8 | g = 1 | Δ^{−25/2} | Δ₀^{2+11−2−25/2} = Δ₀^{−3/2} |
| J = 2 even | −(Δ−10)(Δ−8)(Δ−6)/64 | g = 3 | Δ^{−17/2} | Δ₀^{2+11−6−17/2} = Δ₀^{−3/2} |
| J = 2 odd | i(Δ−9)(Δ−7)/8 | g = 2 | Δ^{−21/2} | Δ₀^{2+11−4−21/2} = Δ₀^{−3/2} |

* **Zeros.** The zeros of F (Δ = 6, 8; 7; 6, 8, 10; 7, 9) are among the double zeros of K² (sin²(πΔ/2) or cos²(πΔ/2)), so
  |K|²/|F_coef|² is finite and continuous, and |a|² = |λ|²|K|² = |F|²(|K|²/|F_coef|²).
* **Exponents.** They are elementary consequences of the closed forms (Stirling for fixed J). The scan
  (`cache/ap7_out.txt`) shows the per-shell sup of Δ^{13}|K|²/(4^Δ|F_coef|²) falling by 0.105–0.116 per factor 4 in Δ₀
  (Δ₀^{−3/2} gives 0.125).
* **Conclusion.** Σ_{low spin}Δ²|a|² < ∞.

---

## 7. Assumption ledger for 𝒟(c_T) = 0

| # | input | status |
|---|---|---|
| 1 | (U) unitarity/reflection positivity; (OPE) convergence and s–t crossing; (T) unique T, finite c_T | standard axioms |
| 2 | (St) dimension of the conserved ⟨TTO⟩ space (2+1 even J, 1 odd J) | [S]; [C] for J ≤ 30 |
| 3 | (L) helicity ladder ratios W⁰_J, W⁴_{J−4} | [P \| D_z·[K,E] = 0 ([D]) and the all-m norm formula ([DV], J ≤ 30)] |
| 4 | ANGULAR_ROUTE §§13–16 (endpoint rule, comparable-twist bound, R-tail criterion) | [P] there, given the retained physical inputs (ANGULAR_AUDIT §2) |
| — | everything in §§1–6 above | [P] relative to 1–4 |

---

## 8. Reproduction

```
cd ads4/code/angproofs                       # KERNELS_CACHE is redirected by common.py
python3 ap1_fit.py                 # candidate symbolic-J vectors (heuristic; out-of-sample checks)    ~40 s
python3 ap2_conservation.py        # [P] conservation + Bose symmetry, symbolic J                      ~1 min
python3 ap2_conservation.py controls   # negative controls
python3 ap3_collinear.py           # [P] the three collinear identities, symbolic J                    ~20 s
python3 ap4_extreme.py pp; python3 ap4_extreme.py pm odd oddJ   # extreme amplitudes, vs pipeline       ~3 min
python3 ap5_kernels.py             # [P] closed forms / vanishing, symbolic J                          ~30 s
python3 ap9_norm_kernel.py         # [P] N_{-J}; assembly of |K_0|^2, |K_4|^2 = KERNELS_D3 closed forms
python3 ap6_envelopes.py           # [P] factorised kernels, explicit envelope constants, scan
python3 ap8_blocks.py              # [P] block lower-bound constants c(delta), numerical check
python3 ap7_lowspin.py             # J = 0, 2 collinear F and the Delta^{-3/2} shell scaling
```

`ap4` needs `cache/kern_J{4,14,15}_*.pkl` (copied from `angaudit/cache`) only for its numerical cross-check against the
pipeline; the symbolic computation does not use them.

---

## 9. Supplement text

To paste as a subsection of `\section{An argument without positivity: the angular route}` (house style of
`prl/supplement.tex`; about one page).

```latex
\subsection{The finite-$c_T$ spectral estimate: collinear identities and two-scale crossing}
\emph{Collinear identities} [P, all $J$]. Place $x_1=(0,-1,0)$, $x_4=0$, $x_3=(0,1,0)$ on a line of the plane $t=0$,
polarise $\mathcal O$ by $e_+=(0,1,i)$ (a 2d quasi-primary of weights $(h,s)=(\frac{\Delta-J}2,\frac{\Delta+J}2)$) and let
$F_{m_1m_3}$ be the coefficient of the 2d three-point structure in the planar-spin components $T_{m_1}$, $T_{m_3}$.
There $V_4=-\frac12$, $H_{14}\propto a_1^2$, $H_{34}\propto a_3^2$, $\epsilon_{14}\propto a_1b_1$, so only the $V_4^J$
monomials reach $m_1,m_3\le-1$. The canonical structures are explicit polynomial vectors in $(s,J)$, conserved and Bose
symmetric for symbolic $J$ (checked as rational-function identities at one conformal frame, which suffices by
covariance), and their extreme detector amplitudes $M_{\mp J}$ follow in closed form by expanding $V_4^{J-k}$
binomially and summing with Chu--Vandermonde; $S^{++}$, $S^{\rm odd}$ have $M_{+J}=0$ and $S^{+-}$, $S^{\rm odd}_J$ have
$M_{-J}=0$. Hence, for all $J$ and in the canonical normalisation of the structures,
$16F_{-2,-2}=p(h)(h-4)(\lambda_{+-}-2\lambda_{++})$, $16F_{-1,-1}=-16p(h)[2(s-1)\lambda_{++}+s\lambda_{+-}]$,
$\pm16iF_{-2,-1}=-2(2h-7)(2h-5)(2h-3)\lambda_{\rm odd}$, $p(h)=(h-1)(h-2)(h-3)$. With the ladder ratios, the absolute
kernels are, e.g., $|K^{++}_0|^2=864\pi^{3/2}\sin^2(\pi h)\,4^\Delta\frac{2^J}{\binom{2J}J}\frac{\Gamma(s+\frac12)}{\Gamma(s)}
\frac{\Gamma(h)}{\Gamma(h-\frac12)}\frac{\Gamma(\Delta-1)\Gamma(\Delta-\frac12)}{\Gamma(\Delta+3)^2}$, and Wendel's inequality gives,
for all $J\ge4$, $h\ge5$, $|K_{\rm even}|^2\le c\,4^\Delta2^{-J}\sqrt{J+1}\sqrt{sh}\,\Delta^{-15/2}$ and
$|K_{\rm odd}|^2\le c\,4^\Delta2^{-J}\sqrt{J+1}\,\Delta^{-15/2}/\sqrt{sh}$ with explicit $c$ ($c_{++}<10^4$).

\emph{Two-scale crossing} [P given unitarity, OPE convergence with $s$--$t$ crossing, a unique $T$ with finite $c_T$, and
the count of conserved structures]. Let $A=T_{m_1}$, $B=T_{m_2}$ and take the mirror ordering
$G(z,\bar z)=\langle A^\dagger(\infty)B^\dagger(1)B(z,\bar z)A(0)\rangle$. Decomposing the 3d Hilbert space into
$sl(2)\times sl(2)$ representations of the plane, $G=z^{-h_A-h_B}\bar z^{-\bar h_A-\bar h_B}\sum_p\frac{|C_{ABp}|^2}{N_p}
k^{(\delta)}_{2h_p}(z)k^{(\bar\delta)}_{2\bar h_p}(\bar z)$ with $k^{(\delta)}_{2h}=z^h{}_2F_1(h+\delta,h+\delta;2h;z)$,
$\delta=h_B-h_A$: every term is $\ge0$ at independent real $z,\bar z\in(0,1)$ (the ordering $\langle ABA^\dagger B^\dagger\rangle$
would give $(h+\delta,h-\delta)$ and indefinite signs). In the crossed channel the identity gives
$N_AN_B(1-z)^{-2h_B}(1-\bar z)^{-2\bar h_B}$, and Cauchy--Schwarz bounds the rest by the identity-subtracted mirror
correlators $\langle BB^\dagger BB^\dagger\rangle$, $\langle A^\dagger AA^\dagger A\rangle$, monotone in each variable
(all non-vacuum weights are $\ge\frac14$) and finite at the Euclidean point $\frac14$. Euler's integral gives
$k^{(\delta)}_{2h}(1-S^{-2})\ge c(\delta)4^hS^{1/2+2\delta}$ for $h\in[S,2S]$, $S\ge4$, $|\delta|\le1$, with explicit
$c(\delta)$. Taking $z=1-U^{-2}$, $\bar z=1-S^{-2}$ and $|C_{AB\,\mathcal O(e_+)}|^2/N=\kappa_02^{-J}|F_{m_1m_2}|^2$:
\[
\sum_{s\in[S,2S],\,h\in[U,2U]}4^\Delta2^{-J}|F_{m_1m_2}|^2\le\mathcal K\,S^{11/2+M}U^{11/2-M},\qquad M=m_1+m_2 ,
\]
with $\mathcal K$ depending on the CFT only through $T$-component norms and two Euclidean four-point values.

\emph{Consequences} [P given the above and the ladder]. Inverting the identities ($p(h)(h-4)\ge0.038h^4$ for $h\ge5$)
gives $\sum4^\Delta2^{-J}(|\lambda_{++}|^2+|\lambda_{+-}|^2)\lesssim S^{3/2}U^{3/2}$ and
$\sum4^\Delta2^{-J}|\lambda_{\rm odd}|^2\lesssim S^{5/2}U^{5/2}$ per shell, hence $\sum_{\rm shell}\sum_b|a_{\mathcal O,b,0}|^2
\lesssim S^{-5}U^2$, $\sum_{\tau\ge10}\Delta^2|a_0|^2<\infty$, $\mu((L,\infty))=O(L^{-3})$ ($O(L^{-4})$ in the band),
$B_N=O(1/N)$ (sharp for this input) and $\mathcal D=0$. The sectors $J=0,2$ close because $F$ grows like $\Delta^0,\Delta,
\Delta^3,\Delta^2$ (even/odd $J=0$, even/odd $J=2$), giving $\Delta_0^{-3/2}$ per dyadic shell. The ladder ratios rest on
light-transform conservation and the all-$m$ helicity norms (checked for $J\le30$); everything else is exact for symbolic
$J$. Derivations and scripts: \texttt{ANGULAR\_PROOFS.md}, \texttt{code/angproofs/}.
```
