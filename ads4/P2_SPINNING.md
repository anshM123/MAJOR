# P2_SPINNING: conserved ⟨TTO_{Δ,J}⟩ in d = 3, heavy symbols to O(1/m²), and P2 of Theorem W (LARGESPIN_W §9.2)

*23 Nov 2026. Unreviewed, not committed by this agent. No existing file was edited. `LARGESPIN_L.md` and `code/lsl/` were not
touched.*
* *New code:* `ads4/code/p2spin/`. Caches are in `cache/`, whose `.gitignore` ignores everything except itself. `run_all.sh` reproduces
  everything in about 0.5 CPU-h, one process.
  * `p2_core.py`: CPPR building blocks (V, H, ε) and exact dual numbers.
  * `p2_series.py`: exact truncated Laurent series in ε = 1/m.
  * Step 1: `s1_structures.py`, `s1b_canonical.py`, `s1c_validate.py`.
  * Step 2: `s2_symbols.py`, `s2_main.py`.
  * Step 3: `s3_mech.py`, `s3_spin.py`, `s3b_leakage.py`.
  * J\*: `s4_jstar.py`.
* *Read-only inputs:*
  * `repro_engine/tto_conserved.py` and `repro_engine/integrand.py` (imported for the J = 0 check);
  * `lsw/cache/costs_1000.json` (Lemma-I cost table).
* *Literature, fetched from arXiv source:*
  * Dymarsky–Kos–Kravchuk–Poland–Simmons-Duffin, *The 3d stress-tensor bootstrap*, 1708.05718 (`TTTT.tex`: structure table and
    appendix);
  * CHMRS 2106.10274: `sec_heavy.tex`, `app_heavy.tex`, and the ancillary `regge_moments_order8.m` ("series Π").

  Note: 1705.04278, as cited in the request, is Dymarsky–Penedones–Trevisani–Vichi (continuous global symmetry). The 3d stress-tensor
  bootstrap is 1708.05718.
* *Compute:* ≈ 0.5 CPU-h, mostly exact sympy. At most 2 processes, `ulimit -v 3000000`.

Grades are as in UNIFIED_THEOREM_v3 §0. **[M]** marks the one modelling step of §3.2. It is flagged wherever it is used.

---

## 0. Summary

| step | result | grade | constants |
|---|---|---|---|
| **1. Conserved ⟨TTO_{Δ,J}⟩, general (Δ, J)** | **Solved exactly over ℚ(Δ, J).** Conservation is imposed with the Todorov operator. The kinematic relations are quotiented out: 14 even monomials with 1 relation, 21 ε-monomials with 9 relations, so 13 + 12 = 25 = 5×5. **Generic J even: 2 parity-even + 1 parity-odd. Generic J odd: 0 parity-even + 1 parity-odd.** Integer J = 0…9 gives 1+1, 0, 1+1, 0, 2+1, 0+1, 2+1, 0+1, 2+1, 0+1. | **[DV] + [Q]**: identical to DKKPS 1708.05718's table | explicit bases rational in (Δ, J): `cache/s1b_canonical.txt` |
| 1. validations | (a) J = 0 is identical to `tto_conserved.py` and `integrand.coef` (exact difference 0). (b) The integer-J solution spaces lie inside the generic family at J = 0, 2, 4–9. (c) An independent mpmath check (30 digits, numerical derivatives, non-integer Δ and J) gives a Todorov divergence at **both** T's ≤ 4·10⁻²⁷ relative. | [DV] | — |
| 1. correction to LARGESPIN_W §1 | The single odd-J structure is an **ε-structure** (parity-odd relative to O). It is the project's "+−". | [DV] | — |
| **2. Heavy symbols to O(1/m²)** | Exact Laurent coefficients R⁽⁰⁾, R⁽¹⁾, R⁽²⁾ of every conserved structure in CHMRS's symmetric frame at fixed β, as polynomials of degree ≤ k in (X, Y). They are obtained by exact interpolation, with zero residual on extra points. (i) R⁽⁰⁾ is **X, Y-independent and pure helicity** in the frame z₁ → B(e^{−β/2}), z₂′ → B(e^{+β/2}) (z₂′ the inversion-referred polarisation): only the modes (±2, ±2) appear. The exact amplitudes are rational in q = e^{β/2}. (ii) The shadow side follows from the Weyl reflection: **R_R(ε) = R_L(−ε) holds exactly** for all four structures (W1a extended to spin). (iii) O(1/m): a constant (basis drift) plus a **purely imaginary transverse linear term**, a Wigner phase. | **[D] + [DV]** (exact rational arithmetic, β = 0.19, 0.81, 1.39) | helicity basis in closed form (§2.2) |
| 2. P1 (parity) | P1 reduces to one group-theoretic input (S): the shadow coefficients are diagonal in the exact collinear helicity basis, with positive leading entries. Given (S), Klein invariance of EE, and the [P] fact that the leading off-diagonal vanishes, **the O(1/m) part of M is diagonal, hence a congruence**. | **[D given (S)]**; (S) [D-sketch] | — |
| **3a. CHMRS ∂_p mechanism, validated** | An exact symbolic map from polynomial insertion to {P, P′} at the state's η. On the scalar seed (Δφ = 3, k = 2) it reproduces CHMRS's ancillary "series Π" at O(1/m²) **exactly** at three β: the ν²P coefficient (11/6)η, the ν²P′ coefficient (η²−1)/6, and the whole ν⁰P′ coefficient. The ν⁰P coefficient differs by **exactly 8η**: the dropped normalisation. | **[DV]** (exact) | α = 2Δφ + k − 1 = 7 |
| 3a. by-product | N2's unidentified scalar residual "−8.0…−9.1·P" (N2 §3.2) equals **−8η·P** at all 8 N2 states, to ≤ 0.05 (e.g. (300, 100): 8η = 9.09 vs 9.06; (1000, 100): 8.13 vs 8.13). | [NE] (functional form); origin [I] | — |
| **3b. P2 magnitudes, helicity-diagonal kernel model** | The model keeps one kernel channel per helicity class. The spinning symbol's O(1/m²) insertion then contributes **exactly zero to the ν² classes B and C2**, for both helicities, both pairing conventions and all three β. It contributes ν⁰ terms only: **+14ηP + 2(η²−1)P′**, helicity-independent (direct pairing). | **[M]** (model) on exact [DV] symbols | **Λ_amp(B) = Λ_amp(C2) = 1 uniformly in the coupling ratios** (2×2 ∝ 𝟙); Λ_A ≤ 3.2, Λ_C3 ≤ 2.0 (direct); conjugate-pairing variant Λ_A ≈ 12, Λ_C1 ≈ 200, Λ_C3 ≈ 6.4 |
| 3c. helicity leakage (R1/P3 channels) | Modes outside the leading classes appear at O(1/m) (≤ 8 relative, plus pure-constant drift into the other ±2 class, removable by congruence) and at O(1/m²) (≲ 10² intrinsic; ≤ 1.8·10³ including drift products). These are what the two non-de-shift kernel channels can see. | [DV] sizes; their kernel projection [I] | §3.3 table |
| 3. not done | No explicit extraction of the three kernel tensors and contraction at O(1/m²); no K_sh; no K₄; no J = 0 engine O(a²) comparison. | — | K_sh, K₄ stay **[M]** |
| **J\* (Theorem W)** | With the model constants (Λ_B = Λ_C2 = 1, Λ_A = 3.17, Λ_C3 = 2.03): **J\* = 50 at margin ×1 even in stripped form (R2)** for (K_sh, K₄) = (10, 10). J\* = 30 for K₄ = 1. At margin ×2, R2 fails at J = 300 and 344, through **K_sh**. | [P given P0–P4] with P2 **[M]** | §4 |

**Bottom line.**
* **Step 1 is done** [DV + Q]. **Step 2 is done at the symbol level** [D + DV], including the reflection identity and the helicity frame.
* **Step 3 is done only within an explicitly labelled model [M].** In that model the spinning block-side O(1/m²) terms contain no ν²
  ring content at all: the decisive classes B and C2 keep their scalar CHMRS values, so Λ_amp = 1 there.
* If the model holds, R1 (P3) is **no longer load-bearing at margin ×1**: J\* = 50 in either form. What decides J ≥ 300 is then K_sh
  (the index-shifted / non-de-shift channels), not the ν² ring classes.
* **Not proved:**
  * that the true three-channel EE kernel acts helicity-diagonally on the O(1/m²) symbol;
  * K_sh and K₄.
* **H_LJ(1000) stays [D-sketch + M].**

---

## 1. Step 1: conserved ⟨T T O_{Δ,J}⟩ in d = 3

### 1.1 Set-up (`p2_core.py`, `s1_structures.py`)
* **Conventions.** These are exactly `tto_conserved.py`'s: light-cone components, embedding P = (1, x², x), Z = (0, 2x·z, z), and
  CPPR's V_{i,jk}, H_ij.
* **Even monomials.** `V₁^{2−a−b}V₂^{2−a−c}V₅^{J−b−c}H₁₂^aH₁₅^bH₂₅^c / (X₁₂^{(10−Δ−J)/2}X₁₅^{(Δ+J)/2}X₂₅^{(Δ+J)/2})`, with
  a + b ≤ 2 and a + c ≤ 2: 14 monomials.
* **Odd monomials.**
  * ε₁₂ = ε(Z₁, Z₂, P₁, P₂, P₅)X₁₂, ε₁₅ = ε(Z₁, Z₅, P₁, P₂, P₅)X₁₅ and ε₂₅ = ε(Z₂, Z₅, P₁, P₂, P₅)X₂₅, times even monomials of
    degree (1, 1, J), (1, 2, J−1) and (2, 1, J−1): 21 monomials.
  * They share the common factor (X₁₂X₁₅X₂₅)^{−1/2}, which is moved into the X-prefactor, so that every row is rational.
* **Constraints.** Each is imposed as a function identity at random rational configurations, in exact Fractions:
  1. **(C1) Conservation at x₁.** Apply the d = 3 Todorov operator to the z₁-quadratic form, then the x₁-divergence. First x₁-derivatives
     are taken exactly with dual numbers, and the factor V₅^J X^{−α} is carried as e^L. Each row is A + ΔB + JC.
  2. **(S) Bose symmetry.** f(1↔2) = f, with V₅^J → σV₅^J and σ = (−1)^J as a sector label.
* **Solve.** The null space is computed over the fraction field ℚ(Δ, J) (sympy DomainMatrix). The kinematic relations are the null space of
  the monomial-evaluation matrix; they solve everything trivially and are quotiented out.

### 1.2 Result (`cache/s1_structures.txt`, `cache/s1b_canonical.txt`)

| sector | generic J even | generic J odd | relations |
|---|---|---|---|
| parity-even (14 monomials) | **2** | **0** | 1 (the d = 3 Gram relation, coefficients 4,4,1,4,2,1,4,2,2,2,1) |
| parity-odd (21 monomials) | **1** | **1** | 9 |

* **Explicit bases.**
  * **Even.** Two vectors rational in (Δ, J), normalised on H₁₂² and V₁²V₂².
  * **Odd.** One vector each, normalised on ε₁₂V₁V₂.
  * The denominators are gauge artefacts of the relation quotient: for example the cubic
    Δ³ − Δ²J − 37Δ² − 27ΔJ + 116Δ − 8J² + 52J − 96 in the even basis.
* **Helicity basis.** The physical (helicity) combinations are fixed in §2.2.

### 1.3 Validations (`cache/s1c_validate.txt`)
1. **Integer J.** Exact solves on the legal monomials only (V₅ power ≥ 0), with σ = (−1)^J:

   | J | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | even + odd | 1+1 | 0 | 1+1 | 0 | 2+1 | 0+1 | 2+1 | 0+1 | 2+1 | 0+1 |

   * **Literature [Q].** DKKPS 1708.05718 (`TTTT.tex`, §2 table and appendix) gives O₀: 1⁺+1⁻; O₂: 1⁺+1⁻; O_{2n≥4}: 2⁺+1⁻;
     O_{2n+1≥5}: 1⁻. It also states "a unique parity-odd structure for odd ℓ ≥ 5". This is an exact match.
   * **Containment.** Every integer-J space lies inside the generic family specialised at that J (rank test on function values).
2. **J = 0.** `(V₁²V₂², V₁V₂H₁₂, H₁₂²) ∝ (Δ(Δ+2)/(2(Δ²−6Δ+6)), 2Δ(Δ−4)/(Δ²−6Δ+6), 1)` is identical to `tto_conserved.py` (difference 0)
   and to `integrand.coef` at Δ = 300 and 2.5.
3. **Independent check.** mpmath at 30 digits, with numerical x- and z-derivatives, at (Δ, J) = (7.3, 4.6) and (12.1, 9.3). The Todorov
   divergence of each canonical structure at x₁ **and at x₂** is ≤ 4·10⁻²⁷ of the single-monomial scale.

**Remark (correction of record).** LARGESPIN_W §1 lists "odd J ≥ 5: one, +−". This is right, but the structure is an ε-structure (§2.2
confirms that its helicity class is the project's "+−").

## 2. Step 2: heavy symbols to O(1/m²) (`s2_symbols.py`, `s2_main.py` → `cache/s2_main.txt`)

### 2.1 Definition
* **Frame.** CHMRS's symmetric frame: x₁ = −y, x₂ = −x/(−x²), x₃ = y, x₄ = x/(−x²), x₅ = e. Rescale x, y = εX, εY with ε = 1/m, and
  set Δ + J − 1 = q/ε, Δ − J − 2 = 1/(qε), q = e^{β/2}. So J + ½ = m sinh(β/2), as in LARGESPIN_W.
* **Polarisations.** e = (1, 0, 0) and n = (1, 1, 0) in Cartesian coordinates. The far polarisation is inversion-referred,
  z₂ = I(x)z₂′, which is ε-independent.
* **Left symbol.** `R_L = Σ_μ c_μ(Δ, J) N_μ / X₁₂²`. This is the spinning ⟨TTO⟩ divided by CHMRS's scalar seed with Δφ = 3. For the
  ε-structures there is an extra exact factor (|X₁₂X₁₅X₂₅|)^{−1/2}.
* **Right symbol.** R_R is built at (x₃, x₄) with the Weyl-reflected coefficients c_μ(3−Δ, −1−J).
* **Expansion.** Exact Laurent series in ε. The order-k coefficients are degree-≤ k polynomials in (X, Y), obtained by exact rational
  interpolation from 32 points, with 4 extra points having residual exactly 0.

### 2.2 Leading order: pure helicity [D + DV]
In the frame z₁ = B(q⁻¹)(1, cos t₁, sin t₁), z₂′ = B(q)(1, cos t₂, sin t₂), with B(s) the boost along n by rapidity ln s, the leading
symbols contain **only** e^{i(±2t₁ ± 2t₂)}. Exact amplitudes:

| structure | modes (2,2) = (−2,−2) | modes (2,−2) = (−2,2) |
|---|---|---|
| even basis 0 | −q²(2q⁶+4q⁴−6q²+1)/(2(q²+1)²) | −q²(2q⁶+4q⁴−6q²−1)/(2(q²+1)²) |
| even basis 1 | q⁴/4 | q⁴/4 |
| odd, J even (OO) | 0 | ∓2iq³/(q⁶−5q⁴−5q²+5) (antisymmetric) |
| odd, J odd | ±2iq³/(q⁶−5q⁴−5q²+5) (antisymmetric) | 0 |

* **The two classes.** The leading symbol is X, Y-independent (checked at two points, exactly). The helicity classes are:
  * **"same"** (modes (±2, ±2)) = the project's **+−**. The inversion reverses orientation, and the odd-J structure lives here.
  * **"opp"** (modes (±2, ∓2)) = the project's **++**. OO is its antisymmetric partner.
* **Helicity basis for even J** (unit leading amplitude):
  * same = −(q²+1)²/q²·b₀ − 2(2q⁶+4q⁴−6q²−1)/q⁴·b₁;
  * opp = (q²+1)²/q²·b₀ + 2(2q⁶+4q⁴−6q²+1)/q⁴·b₁.

### 2.3 Orders 1 and 2 [DV] (`cache/s2_main.txt`; β = 0.1906, 0.8109, 1.3863)
* **Reflection identity.** R_R(ε) = R_L(−ε) holds exactly for all four structures. This is W1a for spinning structures: the frame points
  coincide, c(3−Δ, −1−J) = c(Δ, J)|_{m→−m}, and the referral is the same.
* **O(ε), leading mode.** It is a constant plus a purely imaginary transverse-linear term. At q = 11/10:
  * same: 44.16 − (221/55)i(X_y − Y_y);
  * opp: −3.74 + (21/55)i(X_y + Y_y).

  The X, Y part is a Wigner phase, and it is identical for same and odd, and for opp and OO. The constants are basis drift.
* **O(ε²), leading mode.** Constant plus a real and imaginary quadratic polynomial. The real quadratic part contains (X_z + Y_z)²-type
  pieces, i.e. momentum-shift insertions. The X, Y-dependent parts are identical across structures of one class, apart from one
  linear coefficient.
* **Leakage.** Other modes appear from O(ε). See §3.3.

### 2.4 P1, sharpened [D given (S)]
Write M = g·B, where:
* g_{hh′} = Ω[∫⟨TTO⟩_h⟨ÕTT⟩^{refl}_{h′}] in the Weyl-reflected exact helicity basis;
* B is the matrix of shadow coefficients.

Four facts then fix the O(1/m) part of M:
1. **Reflection.** By the reflection identity and the Klein (13)(24) invariance of the EE kernel (LARGESPIN_W §2.3), g(−m) = g(m)ᵀ.
   Hence g₁ is antisymmetric, with zero diagonal.
2. **(S), the one input.** The stabiliser SO(2) × parity of the collinear frame (0, e, ∞) separates the three structures, because they have
   distinct (|h₁ + h₂|, parity) labels. The Lorentzian shadow transform acts on O's leg only. So B is diagonal in the exact helicity
   basis, B = diag(b_h), and its leading entries are positive.
3. **Symmetry.** M is symmetric, by Klein invariance.
4. **Leading order.** g₀ is diagonal ([P], EE_FUNCTIONAL).

**Conclusion.** (g₁)_{hh′}(b_h + b_{h′}) = 0, so g₁ = 0. M₁ = g₀B₁ is diagonal: a normalisation, i.e. a congruence (W1c).

**Grade.** [D given (S)]. (S) is [D-sketch]: group theory, not checked against KKSD's explicit 3d shadow coefficients. **P1: [D-sketch] →
[D given (S)].**

## 3. Step 3: the ∂_p mechanism, the P2 magnitudes in a model, and the leakage

### 3.1 The mechanism, validated exactly (`s3_mech.py`)
* **Definition.** An insertion X^αY^β becomes (½∂_P)^α(½∂_P̄)^β acting on |P|^{−a}|P̄|^{−a}𝒫_ν(−P·P̄/|P||P̄|). The Fourier–Laplace
  identity is used twice, with independent P and P̄. Then 𝒫″, 𝒫‴, … are reduced by Legendre's equation.
* **Seed series.** CHMRS's scalar seed T(x, y) (Δφ = 3, d = 3) is expanded exactly:
  * s₁ ≡ 0 (CHMRS evenness);
  * the leading exponent gives |P| = |P̄| = 1 and −P·P̄ = cosh β, as it should.
* **Result** (scalar seed at O(1/m²), a = 7):

  | q | ν²P | ν²P′ | ν⁰P′ vs "series Π" (α = 7) | ν⁰P: mechanism − series Π |
  |---|---|---|---|---|
  | 11/10 | (11/6)η | (1/6)(η²−1) | identical (difference 0) | 8.1458 = **8η** |
  | 3/2 | (11/6)η | (1/6)(η²−1) | identical | 97/9 = **8η** |
  | 2 | (11/6)η | (1/6)(η²−1) | identical | 17 = **8η** |

* **The 8η residue.** The mechanism omits the moment's normalisation factors (1/b_{Δ,J}, m^{2k}, γ's, constant logs), which only multiply P.
  The difference is therefore exactly ν-flat. Its value is 8ηP.
* **Link to N2 [NE].** N2 §3.2 found, for the project's exact scalar skeleton against CHMRS (3.20), a residual c·P with
  c = −8.01/−8.03/−8.05/−8.13 (1000; J = 40/48/64/100) and −8.07/−8.18/−8.50/−9.06 (300; J = 24/40/64/100). Using
  η = 1 + 2(J+½)²/m_gap²:
  * 8η = 8.02/8.03/8.06/8.13 at Δ_g = 1000;
  * 8η = 8.09/8.23/8.52/9.09 at Δ_g = 300.

  So **c = −8η** to ≤ 0.05 (≤ 0.6 %) at all 8 states [NE]. That is the functional form of N2's unexplained residual: a ν-flat ηP (class A) term. Its
  origin (which normalisation convention) is not identified [I].

### 3.2 P2 magnitudes in the helicity-diagonal kernel model [M] (`s3_spin.py` → `cache/s3_spin.txt`)
* **Model [M].** The EE kernel acts on each helicity class t ∈ {same, opp} as the scalar Regge-moment kernel times a positive class
  normalisation, and annihilates every other mode.
  * At leading order this reproduces the exact facts: pure-helicity leading symbols (§2.2) and the helicity-diagonal leading image
    (EE_FUNCTIONAL, [P]).
  * Beyond leading order it is **not derived**. The kernel has three polarisation channels (LARGESPIN_W §9.1), and the model keeps one per
    class.
* **Computation.** In the model, M_{hh′} ∝ Σ_t n_t L_{h,t}(ε)L_{h′,t}(−ε)·(1 + ε²s₂). Here L is the class-t coefficient of structure h,
  normalised to the class leading amplitude. The ε² polynomial is pushed through §3.1. Two pairing conventions of the conjugate modes are
  run ("direct", "conj").
* **The X, Y-independent part.** It is basis and normalisation drift, a congruence class under W1c, and is listed separately.

**Results (exact rational, three β):**

| class | seed (= scalar CHMRS moment) | spinning symbol, direct pairing | spinning symbol, conj pairing |
|---|---|---|---|
| **ν²P (B)** | (11/6)η | **0** (both helicities, all β) | **0** |
| **ν²P′ (C2)** | (1/6)(η²−1) | **0** | **0** |
| ν⁰P, X, Y-dependent (A) | — | **+14η** exactly (helicity-independent) | helicity-dependent, e.g. −145.3 (same) / −33.3 (opp) at q = 2 |
| ν⁰P′ (C1/C3) | CHMRS value | **+2(η²−1)** exactly (C3) | fit ≈ −34 (C1) − 10.6(η²−1) (C3), same; smaller for opp |
| O(ε) part of M | — | **0** on the diagonal; off-diagonal antisymmetric (± 30.1, ± 102.5 at q = 3/2) | same |
| O(ε²) off-diagonal (same ↔ opp) | — | 0 in every ring class | 0 |

* **Why the ν² content vanishes.** The real quadratic insertions of the symbol are momentum-shift combinations, (∂_{P_z} + ∂_{P̄_z})²-like,
  for which ∂η = 0 at the state. So no 𝒫″ (no ν²) is generated. Products of the imaginary Wigner-phase terms are likewise ν²-free.
  This is exact within the model. It agrees with LARGESPIN_W §3.5's expectation that ν² coefficients are flat-space-universal.
* **Model premise constants**, against CHMRS's scalar κ (κ_A = 6.46, κ_B = 11/6, κ_C1 = κ_C2 = 1/6, κ_C3 = −11.75/6), sign-free:

  | pairing | Λ_B | Λ_C2 | Λ_C1 | Λ_A | Λ_C3 |
  |---|---|---|---|---|---|
  | direct | **1** | **1** | 1 | (6.46 + 14)/6.46 = 3.17 | (1.96 + 2)/1.96 = 2.03 |
  | conj | **1** | **1** | ≈ 205 | ≈ 11.5 | ≈ 6.4 |

  These hold uniformly in the coupling ratios: the 2×2 is ∝ 𝟙 for B, C2 and the direct-pairing ν⁰ classes.

**Validation against the graviton-window data [NE]: inconclusive.**
* The model gives B = 11/6 = 1.83 per unit.
* The fitted B coefficients are 3.7–5.0 at (1100, 100) (DIRCV §8 / N2 §3.4) and 0.11 at (364, 64)@300. These bracket 1.83, at ×2–2.7 and
  ×0.06.
* The fits have 1–2 dof, so this neither confirms nor refutes Λ_amp(B) = 1.
* The J = 0 engine's O(a²) images were **not** compared: the k = 4 engine's normalisation and light-ray x₀ integral were not set up here.

### 3.3 Helicity leakage: what the two non-de-shift channels can see (`s3b_leakage.py` → `cache/s3b_leakage.txt`)
Maximum coefficient (constant, X, Y polynomial), relative to the structure's own leading amplitude:

| structure | O(ε): other ±2 class (pure constant) | O(ε): \|k\| < 2 modes | O(ε²): other ±2 class | O(ε²): \|k\| < 2 modes |
|---|---|---|---|---|
| same (+−), β = 0.19 / 0.81 / 1.39 | 45 / 102 / 225 | ≤ 8 / ≤ 4 | ≤ 2.5·10⁴ (drift products) | ≤ 1.8·10³ / ≤ 9·10² |
| opp (++) | 3.8 / 30 / 111 | ≤ 8 / ≤ 4 | ≤ 1.2·10⁴ | ≤ 9·10² / ≤ 4.5·10² |
| OO, odd-J (single pure structures) | 0 | ≤ 8 / ≤ 4 | 0 | ≤ 73 / ≤ 42 |

**Reading.**
* **O(ε) leakage into the other ±2 class is X, Y-independent.** It is basis drift of the fixed (q-only) helicity basis and is removable by
  an O(1/m) congruence (W1c). By §2.4 its O(1/m) effect is antisymmetric anyway.
* **The intrinsic leakage is ≲ 10 at O(ε) and ≲ 10² at O(ε²).** This is the |k| < 2 content, best measured on the single pure structures
  OO and odd-J. Contracted with the non-de-shift channels, whose singular values are 0.24–0.37 and 3·10⁻⁴–4·10⁻³ of the leading one
  (LARGESPIN_W §9.1), this suggests κ ≲ 30 in c₂ units [M].
* **Functional class.**
  * The channel functions are rational-in-ν combinations of P_{l+s}: [D-sketch], from EESTEP §1.3's B^{(b)}[n, m] structure.
  * The ∂_p mechanism turns any polynomial insertion into ring derivatives at the state's β.
  * So these terms are **index-shifted harmonics with rational coefficients, possibly carrying ν² factors** [D-sketch].
  * K_sh tolerates κ ≲ 52–59 at J = 300–344 (stripped), so the estimate sits at the margin there and well inside at J ≤ 250. The ν²-weighted
    index-shifted classes have no tabulated costs yet.
* **R1/P3 status.** It is localised further, to "do the |k| < 2 and drift-free leakage polynomials project onto ℓ₂, ℓ₃ with κ ≲ 50?".
  **[I]**, with the explicit inputs above.

## 4. Theorem W's J\* with the computed / model constants (`s4_jstar.py` → `cache/s4_jstar.txt`)
* **Method.** Same criterion and cost table as `lsw_jstar.py`. Its rows are reproduced exactly: (1, 10, 1) → 30/38;
  (3, 10, 10) → 50/100; R2 fails at 300, 344.
* **Premise constants.** Class-wise Λ_s from §3.2.

| premise constants | form | θ = 1 | θ = ½ |
|---|---|---|---|
| model-direct, K_sh = 10, K₄ = 1 | R1 / **R2** | **30** / **30** | 38 / none (300, 344) |
| **model-direct, K_sh = 10, K₄ = 10** | R1 / **R2** | **50** / **50** | 80 / none (300, 344) |
| model-direct, K_sh = 30, K₄ = 10 | R1 / R2 | 50 / none (300, 344) | none (300, 344) |
| model-direct, K_sh = 50, K₄ = 10 | R1 / R2 | 50 / none (300, 344) | none (250, 300, 344) |
| model-conj, K_sh = 10, K₄ = 10 | R1 / R2 | 50 / none (200, 300, 344) | 100 / none |
| *previous (LARGESPIN_W): Λ_amp = 3, (10, 10)* | R1 / R2 | 50 / none (300, 344) | 100 / none |

**Reading.**
* With Λ_B = Λ_C2 = 1, the ν² ring classes stop being binding at large J even in stripped form. **R1 is not needed for J\* = 50 at
  margin ×1.**
* At margin ×2, J\* improves 100 → 80 (R1).
* The binding premises are now **K₄** (J < 50, the full-range ν⁴ artefact, unchanged) and **K_sh** (J ≥ 300). Both are uncomputed [M].
* The conj variant is worse at J = 200 in R2, through the large C1 coefficient.

## 5. Exactly what remains
1. **The true three-channel kernel contraction at O(1/m²)** (replaces the [M] of §3.2).
   * Extract 𝒯₁,₂,₃ as explicit tensors from ee6's L_ν[pols]. It is linear in the dressing, so one can run ee_quad's Laurent decomposition on
     the tensor components of x₁₀^⊗2 ⊗ x₃₀^⊗2 ⊗ M ⊗ M (1296 components).
   * Contract with the exact symbols of §2: R⁽⁰…²⁾, in `cache/s2_main.pkl`.
   * Push the result through §3.1. This requires the harmonic decomposition of ℓ_i(ν, B) with rational coefficients, or a 4-d angular
     quadrature on H² × H² after the analytic radial integrals.
   * Decides: Λ_amp beyond the model, and the projection of the leakage (§3.3) onto ℓ₂, ℓ₃ (R1/P3).
2. **K_sh.** The kernel's own O(r) terms, EESTEP's B^{(b)}[1, m], contracted with R⁽⁰⁾. Also tabulate costs for "ν² × P_{l±s}" classes in
   `lsw_costs.py`'s format.
3. **K₄.** CHMRS's order-1/m⁴ ancillary series applies to the seed. The symbol expansion extends to ORDER = 4 by changing one constant in
   `p2_series.py`. This is cheap.
4. **(S) for P1.** Check diagonality and the sign of the 3d shadow coefficients in the helicity basis (KKSD 1809.05111).
5. **Validation targets not reached.**
   * The J = 0 engine's O(a²) images: needs the k = 4 light-ray normalisation.
   * The (1100, 100) B fit: the model's 1.83 sits inside the spread of the data, which is too noisy to decide.
6. **Origin of N2's −8η·P** (normalisation convention): [I], harmless (ν-flat, class A).

## 6. Effect on the Δ_g = 1000 inventory (suggested; not applied)
* **H_LJ(1000).** It stays **[D-sketch + M]**.
* **What changes:**
  * P1 becomes [D given (S)].
  * P2's ν² part is Λ_amp = 1 **in a stated model** [M], not [NE]. P2's structure is now explicit: exact symbols, a validated mechanism,
    and the helicity basis.
  * R1 is no longer needed at margin ×1 if the model holds. Its residual is localised in the helicity leakage (§3.3) and K_sh.
* **Do not write** "Λ_amp = 1 is derived", or "proved for J ≥ 50".

---

## 7. Follow-ups (coordinator request): K₄, the (S) check, the J = 0 engine, the −8η factor, and J\*

*Same day. The three-channel kernel contraction was not attempted: it is covered by `PW_IMAGES.md` / `code/pwimg/` (another agent;
not touched). New code: `s5_order4.py` (O(1/m⁴) seed + symbols), `s5b_dphi.py` (Δφ-dependence of CHMRS's O(1/m⁴) top-ν terms),
`s6_nu6cost.py` (Lemma-I cost of the ν⁶ shape). Compute ≈ 1.5 CPU-h.*

### 7.1 (a) O(1/m⁴): the seed, validated, and a gap in P4's representatives

**Seed at O(1/m⁴) [DV], exact.**
* The mechanism of §3.1 is extended with memoised derivatives and Legendre reduction to order 10.
* It is compared with CHMRS's ancillary "series Π" at order 1/m2², parsed from `regge_moments_order8.m` and Legendre-reduced. At
  β = 0.81 and 1.39 (q = 3/2, 2):
  * the **ν⁶ and ν⁴ coefficients of P and P′ are identical**;
  * the ν² and ν⁰ differences equal **exactly** 8η × (my first-order coefficients);
  * so CHMRS's series = (1 − 8η/m² + n₂/m⁴) × (mechanism), i.e. the same ν-flat normalisation as at first order (§3.1).
* The O(1/m⁴) mechanism is therefore validated exactly (`cache/s5_order4.txt`).

**The O(1/m⁴) seed contains a ν⁶ term.** By `s5b_dphi.py` (Δφ = 1, 2, 3, 5, 8; k = 2; d = 3):

| coefficient at O(1/m⁴) (per unit leading P) | value | Δφ-dependence |
|---|---|---|
| ν⁶ P | **−(η²−1)/72** | none (flat-universal) |
| ν⁴ P′ | (11/30) η(η²−1) | none (flat-universal) |
| ν⁴ P | (5189η² − 1595)/1440 at Δφ = 3 | yes (Δφ = 5: (389η² + 3205)/1440) |

* **K₄ units.** Theorem W's P4 bounds the O(1/m⁴) remainder by K₄ in units of {η²ν⁴P, η sinhβ ν⁴∂_βP} only (LARGESPIN_W §4.0, caveat 4).
  * In those units the scalar seed has K₄ = |5189η² − 1595|/(1440η²) ≤ 3.6 (ν⁴P) and 0.367 (ν⁴P′). That is within K₄ = 10.
  * **But the ν⁶ shape is not among the representatives.** Its Lemma-I cost per unit, relative to the 4a unit
    (`cache/s6_nu6cost.txt`, R2), is 2.5·10³ (J = 26) … 2.3·10⁵ (J = 120). With the 1/72 coefficient this is **equivalent to K₄ ≈ 35
    (J = 26), 135 (J = 50), 410 (J = 100), 2·10³–3·10³ (J = 120–200), 560 (J = 344)**.
* **Reading [D-sketch].**
  * The ν⁶ term is Δφ-independent, i.e. part of CHMRS's flat-space top-ν tower. These are the terms the hybrid floor resums (flat partial
    wave, LARGESPIN_W §3.5 / §4.3 note).
  * They must therefore be **excluded from P4, not costed**. P4 should be restated as "the beyond-hybrid O(1/m⁴) part ≤ K₄".
  * If they are costed as Regge-series terms over the full range ν ≤ Λ, Theorem W fails at every sampled J. This is the full-range artefact
    that LARGESPIN_W already flags for K₄ at J < 50, now visible at all J.
  * **This is a gap in Theorem W's P4 as stated**, independent of spin. It needs either the explicit hybrid subtraction at O(1/m⁴) or a
    restatement. Graded [I] for the subtraction; the identification of the ν⁶ term is [DV].

**Spinning symbols at O(1/m⁴), helicity-diagonal model [M].** See §7.1b; the numbers are filled in from `cache/s5_order4.txt`.

### 7.2 (b) The (S) check against Karateev–Kravchuk–Simmons-Duffin 1809.05111: (S) is NOT what KKSD establish

KKSD's source was fetched; the relevant parts are §§2.3, 4.1–4.3 and the 3d sections "Implementation in 3d", "Shadow transforms and
pairings with conserved operators" and "OPE coefficients of 3d stress-tensor MFT".

**What KKSD say [Q].**
1. **Shadow transforms are matrices.** The shadow transform of a three-point structure is a matrix S(…)^d_e on the structure space. It is
   in general **non-diagonal**: their stress-tensor MFT formula carries S([T̃₁]T̃₂𝒪)^d_e S(T₁[T̃₂]𝒪)^e_a.
2. **Diagonality holds only in momentum space.** The shadow transform is diagonal only in the **momentum-space SO(d−1) basis of the
   shadowed leg**, with eigenvalues 𝒜_{λ_{d−1}}. The Fourier transform is diagonal in the position-space SO(d) basis.
   * The algorithm needs two **non-diagonal** changes of basis between them ("the coefficients ℱ_{Δ,J} and 𝒜_{λ_{d−1}} for different J and
     λ_{d−1} will mix").
   * In 3d, the position-space structures are labelled by SO(2) weights [q₁, q₂, q₃] of the collinear frame. Only the **pairing** is
     diagonal in that label ([q] ↔ [−q]), not the shadow transform.
3. **Conserved T adds gauge redundancy.** For conserved T the shadow structures ⟨T̃T̃𝒪⟩ carry the gauge redundancy T̃ ~ T̃ + ∂ξ. This does
   not make S diagonal.

**Consequence.**
* §2.4's assumption (S) ("the shadow matrix is diagonal in the exact helicity basis of the T legs") is **not supported**. The two
  parity-even structures (++, +−) can mix under the shadow transform. (The parity-odd one cannot.)
* My §2.4 argument that "SO(2) × parity separates them" was wrong. Equivariance holds for SO(d) on the T legs at (0, ∞), whose labels are
  the SO(3) channels of 2 ⊗ 2. It does not hold for the collinear SO(2) helicities.
* **P1 therefore stays [D-sketch]. The upgrade of §2.4 is withdrawn.**

**What survives.**
1. **[D, from known facts] B₀ is diagonal and positive.** The leading heavy block is helicity-diagonal ([P], EE_FUNCTIONAL), and the
   leading g₀ is diagonal with positive entries (§2.2). Since M₀ = g₀B₀, the leading shadow matrix B₀ is diagonal and positive in the
   helicity basis.
2. **[DV] + [M] The model-level O(1/m) off-diagonal is constant.** At O(1/m) the off-diagonal symbol content between the two parity-even
   classes is X, Y-**independent** (exact, §3.3: the "other ±2 class" leakage at O(ε) has zero X, Y-coefficients). So in the
   helicity-diagonal model, the O(1/m) off-diagonal of M is a constant multiple of the leading image: a congruence, whatever B₁ is.
3. **Model-independently, still [D-sketch].** The X, Y-dependent O(1/m) content sits only in |k| < 2 modes (≤ 8 relative). Its kernel
   projection is the three-channel question (PW_IMAGES).
4. **Next step if wanted.** KKSD's 3d algorithm (their eqs. for 𝒜_{1,m} and the SO(3) ↔ SO(2) Clebsch–Gordan change of basis) gives
   S(T₁[T̃₂]𝒪_{Δ,J}) explicitly. In the heavy limit one would check that its off-diagonal (++ ↔ +−) part is O(1/m) and odd.
   Not done here.

### 7.3 (c) The J = 0 engine O(a²) comparison: not done

* **What the engine is.** `repro_engine/integrand.py` is the **k = 4** EE integrand at **J = 0** in the a = 1/Δ expansion. It includes the
  light-ray x₀ integral and its own state/reference normalisation at ν = 0.02. At J = 0 the heavy state sits at β = O(1/m), outside the
  fixed-β expansion of §§2–3.
* **What the comparison would need.**
  * The full EE kernel contraction (§3.2's model is not enough at J = 0, where everything is ++ and the test is of the kernel, not the
    symbol). That is PW_IMAGES' object.
  * The engine's normalisation conventions.
* Neither is cheap. **Step 3's validations remain:**
  * the exact reproduction of CHMRS's scalar series at O(1/m²) and O(1/m⁴) (§§3.1, 7.1);
  * the J = 0 structure identity (§1.3).

### 7.4 The −8η factor
* **Normalisation-side behaviour [DV].** The relation "CHMRS series Π = (1 − 8η/m² + n₂/m⁴) × (bare mechanism)" holds exactly at
  O(1/m²) and O(1/m⁴) (§7.1). So −8η/m² is the first-order term of the moment's **ν-independent normalisation**. In the bare mechanism
  that normalisation is set to |P| = |P̄| = 1 at the rescaled momenta, and 1/b_{Δ,J}, m^{2k}, γ² and the dropped constants are omitted.
* **Not pinned to a single factor.** Which factor carries it (1/b_{Δ,J}'s Stirling tail, or |p|² = Δ² − J² against
  m² = (Δ−J−2)(Δ+J−1)) was not determined; that needs CHMRS's exact normalisation chain.
* **Sign relative to N2.** N2's residual (exact skeleton − (3.20) = −8ηP) has the sign of "CHMRS × (1 − 8η/m²)". The project's skeleton
  therefore carries this normalisation once more than CHMRS's (3.20). This is consistent with a double-counted or differently defined m²
  normalisation, but that is not verified [I].
* **Harmless either way.** It is ν-flat and multiplies the leading image (W1c).

### 7.1b Spinning symbols at O(1/m⁴), helicity-diagonal model [M] (`s5_order4.py` → `cache/s5_order4.txt`)

**Method.**
* Exact series to ε⁴ for the even-J structures at q = 3/2 and 2.
* Polynomials of degree ≤ k in (X, Y), from exact interpolation on 210 monomials with 6 check points (residual 0). The series construction
  is [DV]; the kernel action is the §3.2 model.
* The normalisation drift is factored out as a congruence (W1c). That removes the X, Y-constant parts, the constant at each order, and
  constant × s₂.
* The O(1/m⁴) insertion σ₄ + s₂σ₂ − c₂(s₂ + σ₂) is pushed through the mechanism.

**Results:**

| β | ν⁴P | ν⁴P′ | ν²P (×1/m⁴) | ν²P′ (×1/m⁴) |
|---|---|---|---|---|
| 0.81 | −0.1698 = **−(5/24)(η²−1)** (both helicities) | **0** | 175.5 (+−) / 65.8 (++) | 19.4 / 4.9 |
| 1.39 | −0.7324 = **−(5/24)(η²−1)** (both helicities) | **0** | 426.2 / 217.2 | 72.6 / 36.5 |

**Reading.**
* In the model, the spinning symbol adds to the ν⁴ classes only the helicity-independent −(5/24)(η²−1)ν⁴P. In K₄ units that is
  ≤ (5/24)(1 − η⁻²) < **0.21** at every β.
* **K₄ (model, beyond the flat ν⁶ term)** = the seed's Δφ-dependent |5189η² − 1595|/(1440η²) ≤ 3.60, plus 0.21. So **K₄ ≤ 3.8**,
  uniformly in the coupling ratios and in β. Grade: [M], on [DV] series.
* The new ν²/m⁴ terms (≤ 430/m⁴) are a ≤ 4·10⁻⁴ correction to κ_B at m ≥ Λ. They are negligible.
* **Model-independently**, K₄ is not bounded here. The three-channel contraction (PW_IMAGES) and the hybrid subtraction of the ν⁶ tower
  (§7.1) are both needed.

### 7.5 Theorem W's J\* with the follow-up constants (`s4_jstar.py` → `cache/s4_jstar.txt`)

| premise constants | P2 grade | form | θ = 1 | θ = ½ |
|---|---|---|---|---|
| **model: Λ_B = Λ_C2 = 1, Λ_A = 3.17, Λ_C3 = 2.03, K₄ = 4, K_sh = 10** | [M] | R1 / **R2** | **38** / **38** | 50 / none (300, 344) |
| same, K_sh = 30 | [M] | R1 / R2 | 38 / none (300, 344) | none (300, 344) |
| model, K₄ = 10 (earlier §4), K_sh = 10 | [M] | R1 / R2 | 50 / 50 | 80 / none (300, 344) |
| LARGESPIN_W baseline (Λ_amp = 3, 10, 10) | [NE] / [M] | R1 / R2 | 50 / none | 100 / none |

**Caveats on this table.**
1. **The flat ν⁶ term.** "K₄ = 4" holds only if the flat-universal ν⁶ tower of §7.1 is assigned to the hybrid floor. This is [D-sketch]
   and is now a named premise. If it were costed as a Regge term, every sampled J would fail.
2. **K_sh is the binding premise at J ≥ 300.** It is uncomputed [M] and is PW_IMAGES' object.
3. **P1** stays [D-sketch] (§7.2).

### 7.6 Net of the follow-ups
* **(a) K₄.**
  * The seed's O(1/m⁴) is validated exactly against CHMRS's ancillary series [DV].
  * The model gives **K₄ ≤ 3.8** [M], and J\* improves 50 → **38** (θ = 1, either form).
  * **New gap in P4:** CHMRS's O(1/m⁴) contains a Δφ-independent ν⁶(η²−1)/72 term outside P4's representatives. Its full-range
    Lemma-I cost is equivalent to K₄ ≈ 35–3000. Theorem W needs the hybrid subtraction of this flat tower to be explicit [I].
* **(b) (S).** Not supported by KKSD 1809.05111. Shadow matrices are diagonal only in the momentum-space SO(d−1) basis, and ++ and +− can
  mix. **P1 stays [D-sketch]**; the §2.4 upgrade is withdrawn. What remains [D]: B₀ is diagonal and positive, and the model-level O(1/m)
  off-diagonal is a constant (a congruence).
* **(c) J = 0 engine O(a²).** Not done. It needs the full kernel contraction (PW_IMAGES) and the engine's normalisation.
* **−8η.** It is the first-order term of the moment's ν-flat normalisation [DV, at O(1/m²) and O(1/m⁴)]. The single responsible factor is
  not pinned down. Harmless.
