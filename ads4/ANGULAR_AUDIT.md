# ANGULAR_AUDIT: hostile audit of ANGULAR_ROUTE §21 (finite-c_T spectral defect via crossing + Casimir)

*26 Nov 2026. Independent audit agent, about 2.5 h. Unreviewed and not committed. No existing file was edited. BATCH4_FOLD.md,
FREE_NULLPLANE.md, `code/b4fold/` and `code/freenp/` were not touched.*

* **Code.** New code is in `ads4/code/angaudit/` (outputs in `angaudit/cache/`, which ignores everything but its `.gitignore`).
  `code/kernels/{structures,pipeline,kernels_lib,closed_forms,canonical}.py` are imported read-only, with `KERNELS_CACHE`
  redirected to `angaudit/cache`.
* **Compute.** About 25 CPU-min, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

**Grades.** [P] proved, or exact algebra valid in general; [DV] verified exactly/numerically in finite cases;
[D] derived, plausible; [I] unsupported; WRONG.

---

## 0. Verdict: **confirmed with fixes**

No error was found that invalidates the chain items 1–4 ⇒ 𝒟(c_T) = 0. Every number I could test came out as claimed:

* the three collinear identities hold **exactly**, in the repo's conventions, for every J tested (even J = 4–30, odd J = 5–21);
* the weight 4^Δ2^{−J} is the correct one;
* the two-scale crossing estimate is legitimate;
* the kernel envelopes hold;
* the tail arithmetic is right.

Four fixes are needed before §21 can be called a proof.

1. **All-J status.** The identities (item 1) and the absolute kernel closed forms (item 3) are verified exactly only for J ≤ 30.
   They are [DV], not [P]. The author's "symbolic J" derivation should be written out.
2. **The crossing lemma (item 2) must be written.** It is correct only with the **mirror ordering** ⟨A B B†A†⟩ (A = T_{m₁},
   B = T_{m₂}), with a 2d global-block decomposition at independent real z, z̄, and with a Cauchy–Schwarz bound on the crossed
   channel. With the ordering ⟨A B A†B†⟩ the coefficients carry signs and the argument fails. Grade [D].
3. **Rate.** The chain gives B_N = O(1/N), and this rate is sharp for the crossing input (§4). It does **not** give §18's
   O(N⁻²). It therefore does not give §18's "absolute angular swapping" (that sum is log-divergent at this rate). 𝒟 = 0 is
   unaffected, because the R-tail criterion (§15) needs only μ((L,∞)) = o(L⁻²).
4. **Low spin.** The stated one-line argument ("kernel O(4^ΔΔ^{−29/2}) against a budget Δ^{11}") covers only J = 0 even.
   J = 0 odd and J = 2 have larger kernels (4^ΔΔ^{−25/2}, 4^ΔΔ^{−17/2}, 4^ΔΔ^{−21/2}). They close only because their
   collinear F grows like Δ, Δ³, Δ². With that growth, every low-spin sector gives Δ^{−3/2} per dyadic shell [DV].

**Item 5 (Casimir) is not load-bearing.** 𝒟 = 0 follows from items 1–4 together with §§14–16 alone. Item 5 proves only that
I_η stays finite. It does not prove the null-plane **identity**, and it is no longer needed: §18 assumption 2 is retired, not
proved.

---

## 1. Per-link table

| # | link | grade | finding |
|---|---|---|---|
| 1a | Meaning of F_{m₁m₂} | [P] (convention) | 3 points collinear in the Euclidean plane t = 0; O polarised by exactly e₊ = (0,1,i), i.e. its top planar-spin component, a 2d quasi-primary of weights (h, s) in (w, w̄); T components = coefficients of a^{2+m}b^{2−m} in T(ζ(a,b)); F = coefficient of the 2d 3-pt structure. With O polarised e₋ you get the mirror identities (s ↔ h roles exchanged) |
| 1b | F₋₂,₋₂ = p(h)(h−4)(λ₊₋−2λ₊₊), F₋₁,₋₁ = −16p(h)[2(s−1)λ₊₊+sλ₊₋], F₋₂,₋₁ = −2(2h−7)(2h−5)(2h−3)λ_odd | **[DV]** exact in s | Hold identically in s for even J = 4–20, 24, 30 and odd J = 5–15, 21 (single structure, F₋₂,₋₁ ∝ λ_oddJ with the same polynomial). The common factor (1/16; ∓i/16 for odd) is **J-independent**, so the λ normalisation is the KERNELS_D3 canonical one. No even structure enters F₋₂,₋₁ and no odd structure enters F₋₂,₋₂ or F₋₁,₋₁ |
| 1c | Structure counting | [P] | Even J ≥ 4: 3 structures, 3 identities, invertible. Odd J ≥ 5: 1 and 1. J = 2 and J = 0: 1 even + 1 odd, F computed (§3). J = 1, 3: none. Degenerate operators are harmless, since every bound is a sum of squares |
| 2a | Two-scale bound from reflection positivity + crossed identity | **[D]** | Valid. The 3d Hilbert space decomposes into unitary so(2,2) reps. O(e₊) is a 2d quasi-primary with norm (e₊*·e₊)^J = 2^J, orthogonal to the other families. At independent real z, z̄ ∈ (0,1), ⟨A B B†A†⟩ = Σ_p \|c_p\|² k_{2h_p}(z)k_{2h̄_p}(z̄) with every term ≥ 0, so the shell (s ~ S, h ~ U) is selected by z = 1−U⁻², z̄ = 1−S⁻². The crossed channel is bounded by Cauchy–Schwarz against the positive sums for ⟨B B† B B†⟩ and ⟨A†A A†A⟩ at small arguments, and every non-identity p has h_p, h̄_p ≥ 1/4. **Must use the mirror ordering** |
| 2b | Exponents S^{11/2+M}U^{11/2−M} | **[DV]** | The w-weights of T_{−2}, T_{−1} are 5/2, 2 and their w̄-weights 1/2, 1; O(e₊) has (h, s). The crossed identity gives U^{2(h_A+h_B)−1/2} S^{2(h̄_A+h̄_B)−1/2} = U^{11/2−M}S^{11/2+M}. Calibrated on the 1d GFF ⟨ABBA⟩: Σ_{[U,2U]}c·4^h / U^{2(h_A+h_B)−1/2} → const for six weight pairs |
| 2c | k^{(δ)}_{2h}(1−S⁻²) ≳ 4^h S^{1/2+2δ} on h ∈ [S,2S] | **[DV]** | For δ ∈ {−1, ±½, 0, 1, 2} and S = 8–4096, the ratio lies in [0.009, 1.24], stable in S. Here δ = h_B − h_A, the standard 2F1(h−h₁₂, h+h₃₄; 2h) convention |
| 2d | Weight 4^Δ2^{−J} | **[P]** | 4^Δ = 4^h·4^s comes from the two blocks; 2^{−J} = 1/⟨O(e₊)\|O(e₊)⟩. It is consistent only because F's normalisation is J-independent (1b) |
| 2e | Off-diagonal positivity (m₁ ≠ m₂) | [P] with the mirror ordering | In ⟨A B B†A†⟩ the coefficient is \|F_{m₁m₂}\|²2^{−J}/(N_A N_B) ≥ 0. Note \|F₋₂,₋₁\| = \|F₋₁,₋₂\| |
| 3a | \|K_even\|² ≲ 4^Δ2^{−J}√(J+1)√(sh)/Δ^{15/2}, \|K_odd\|² ≲ 4^Δ2^{−J}√(J+1)/(√(sh)Δ^{15/2}) | **[DV]**; [P] from the closed forms | Pipeline K² (independent `norm_poly` norms) = closed forms to 30 digits for J ≤ 30, h = 5¼–1000¼. Envelope ratio: pp ∈ [2.7e3, 4.3e3], pm ∈ [9.9e2, 2.9e4], odd ∈ [1.3e4, 1.8e4], oddJ ∈ [1.7e4, 1.1e6] over J ≤ 2·10⁴, h ≤ 2·10⁴. Stirling on the closed forms gives the claim (the sin²πh factor is ≤ 1). The closed forms themselves are [DV] (J ≤ 30) |
| 3b | Conditioning of F → λ | **[P]** | After row scaling by p(h)(h−4) ≥ 0.038h⁴ and 16p(h)s ≥ 3.07h³s (h ≥ 5), the 2×2 matrix is [[−2, 1], [2(1−1/s), 1]], with det = −4 + 2/s. Uniformly well conditioned |
| 4a | Σ_shell\|a₀\|² ≲ S^{−5}U²; Σ_{τ≥10}Δ²\|a₀\|² < ∞ | **[P]** given 1–3; [DV] model | Arithmetic correct, including the odd parts and the b = 4 (pm, oddJ) ladders. In the saturating model (§4), shell sums/(S^{−5}U²) ≤ 540 and ΣΔ²\|a₀\|² converges |
| 4b | Low-spin sectors | **[DV] with fix** | See fix 4: Δ^{−3/2} per dyadic shell in all four sectors (J = 0, 2; even, odd). Low spins never enter B_N for N ≥ 4 (a_N = 0 for N > J − b) |
| 4c | μ((L,∞)) = O(L⁻³), O(L⁻⁴) in the band | **[P]** given 1–3 | Uses R = (Δ−1)(Δ−2)/((2h−2)(2s−1)) ≈ (s+h)²/(4sh). Model: L³μ ∈ [2.0, 2.9] for L = 2–64; L⁴μ_band bounded. "Overall" should read "τ ≥ 10"; τ < 10 is handled by §14 |
| 4d | ΣΔ²\|a_N\|² = O(N³), B_N = O(1/N), 𝒟 = 0 | **[P]** | ∫e^{−N/R}dν ≲ N⁻³ via the ladder (C/T)^N ≤ e^{−N/R}; bounded twist O(N²) by §14. Model: N·B_N → 2.64, N²B_N grows ⇒ **O(1/N) is sharp** for this input. See fix 3 |
| 5a | ‖(K⁰)^{1/2}w‖² ≲ ‖∇_p w‖² + (1+Δ²+J²)‖w‖² on a compact timelike patch | [D] | It is an upper form bound. The second-order part has an indefinite but bounded coefficient (~\|p\|). The first-order, spin and measure-derivative terms (∝ Δ, J, Δ/\|p\|) are absorbed by Cauchy–Schwarz. There are no boundary terms for compact support in the open cone |
| 5b | Uniform decay of the boosted-profile Fourier coefficients | [D] | True for a fixed smooth celestial profile F: Möbius aberration with bounded rapidity gives joint smoothness, hence decay faster than any power, uniformly with p-derivatives. **Not** uniform in the harmonic label m of f_m; that uniformity is no longer needed |
| 5c | H vs K⁰ form-comparability | [P] in the needed direction | K^μ lies in the closed forward cone in every unitary positive-energy rep, so every null component satisfies K± ≤ 2K⁰. **Caveat:** the chart is not written. In the Hofman–Maldacena chart, translations along the null plane map to a null **momentum** P₋, and then A(f)ψ ∈ Dom(H^{1/2}) is elementary for momentum wavepackets, since E(F) preserves momentum support |
| 5d | lim I_η < ∞ | [D] | Dominated convergence from Aψ ∈ Dom(H^{1/2}) (needs ψ ∈ Dom(H)). It gives finiteness only, not [A,[H,A]] = ½A((f′)²), and it is downstream of items 1–4 (it uses ΣΔ²\|a_k\|² = O(k³)). **Not needed for 𝒟 = 0** |
| 6 | Scope | — | See §2 |

---

## 2. What is proved, and under which inputs

**Statement reached.** Take a unitary 3d CFT at finite c_T with the retained inputs below. Then 𝒟(c_T) = 0, i.e. the endpoint
rule holds with right-hand side exactly 3. Hence √(t₄²+α_T²) ≤ 4√(3/K(g) + 𝓛_g(c_T)). Taking c_T → ∞ afterwards, with 𝓛_g → 0,
gives limsup √(t₄²+α_T²) ≤ 0.017814 at g = 1000.

**Grade of 𝒟 = 0:** [D], with these parts:
* [DV]: the all-J identities and closed forms, J ≤ 30;
* [D]: the 2d crossing lemma;
* [P]: the arithmetic and §§14–16.

**Retained inputs (unchanged by §21):**
1. physical matching: the all-J ladder and first-moment operator, from [z·K, E] = 6iE⁽¹⁾ and D_z·[K, E] = 0;
2. contact matching (the exact BMS/KZ coefficient);
3. fixed-window planar regularity (𝓛_g → 0);
4. Δ_φ > 1 for scalars in T×T, with conserved higher-spin sectors excluded;
5. KKSZ (4.66), i.e. separated-direction vanishing of [E, E⁽¹⁾] at the light-cone edge, together with existence of the physical
   mixed commutator in momentum eigenstates.

**New inputs introduced by §21:**
* only standard CFT axioms: reflection positivity, and OPE convergence at independent real z, z̄ in both channels;
* the two [DV] all-J extrapolations.

The CFT-dependent crossing constants are harmless, because 𝒟(c_T) = 0 is exact at each c_T before the planar limit.

**What is not delivered:**
* §18's O(N⁻²) rate and absolute angular swapping (only O(1/N));
* the regulated null-plane **identity** (only I_η finiteness, item 5).

Neither is needed for the bound. §18 assumption 2 should be **removed** from the assumption list, not marked proved.

---

## 3. Low-spin data (ca5)

Here F is per unit canonical coefficient, with O polarised e₊.

| sector | F used | \|K\|²/4^Δ | budget per shell | Σ_shell Δ²\|a\|² (ratio between shells 2048 and 512) |
|---|---|---|---|---|
| J = 0 even | F₋₂,₋₂ = (s−4)(s−3)/(16s(s+1)) → 1/16 | Δ^{−29/2} | Δ^{11} | 0.119 (Δ^{−3/2}: 0.125) |
| J = 0 odd | F₋₂,₋₁ = i(2s−7)/8 | Δ^{−25/2} | Δ^{11} | 0.121 |
| J = 2 even | F₋₂,₋₂ = −(s−6)(s−5)(s−4)/8 | Δ^{−17/2} | Δ^{11} | 0.117 |
| J = 2 odd | F₋₂,₋₁ = i(2s−11)(2s−9)/8 | Δ^{−21/2} | Δ^{11} | 0.119 |

The zeros of F coincide with the double zeros of K² at the same Δ:
* J = 0: Δ = 6, 8 (even) and Δ = 7 (odd);
* J = 2: Δ = 6, 8, 10 (even) and Δ = 7, 9 (odd);
* general J: h = 1–4 (even) and h = 3/2, 5/2, 7/2 (odd).

So \|K\|²/\|F\|² stays finite there.

## 4. Saturating-model check of items 3–4 (ca4)

**Model.** One parity-even operator per unit cell of the lattice h = n + ¼ (n ≥ 5), even J = 4–6000, with
4^Δ2^{−J}\|λ₊₊\|² = √(sh). This saturates the claimed shell budget S^{3/2}U^{3/2}, and the closed-form K is used.

**Results.**
* ΣΔ²\|a₀\|² converges (3.95, 4.05, 4.10, 4.11 for cut-offs Δ < 500, 1000, 2000, 3000).
* L³μ((L,∞)) ∈ [2.0, 2.9] for L ∈ [2, 64], and L²μ → 0.
* N·B_N: 0.68, 1.10, 1.59, 2.06, 2.42, 2.55, 2.61, 2.64 at N = 1, 2, 4, 8, 16, 24, 32, 40.

So B_N → 0 exactly at the rate 1/N, and no faster.

## 5. Reproduction

```
cd ads4/code/kernels && KERNELS_CACHE=../angaudit/cache python3 compute_kernels.py J    # J = 14,15,16,18,20,21,24,30 (~1 min each)
cd ads4/code/angaudit
python3 ca1_collinear.py            # item 1, J = 4..12 from data/*.json (~10 min)
python3 ca2_blocks.py               # item 2: block lower bound, GFF calibration (~1 min)
python3 ca3_kernels.py 4 6 ... 30   # items 1+3: pipeline K vs closed forms vs envelope; identities for J >= 14 (~2 min/J)
python3 ca4_tail.py                 # item 4: saturating model (~1 min)
python3 ca5_lowspin.py              # low-spin sectors (~1 min)
```
