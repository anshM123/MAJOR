# Result Q: free-fermion O(1/C_T) check of Lemma A + D_LR, and a scalar-Mellin toy check of H_EFT

*22 Nov 2026, fermion/H_EFT agent. Unrefereed, not committed. No existing document or code file was edited (the `og_*` files included);
no `cl_*` / `spinfunc` file was touched. New code: `ads4/code/resultq/of_fermion.py` (engine), `of_run.py` (task 1), `of_heft.py`
(task 2); outputs `resultq/cache/of_check.txt`, `of_run.log`, `of_heft.txt`, `of_heft.log`, `of_heft_moments.txt`, `of_heft_diag_nc40/64.txt`. Read-only imports:
`og_free` (contour, `f_k`, null basis), `v4c_core`, `v4d_core`, `v4c_wu.Kfun`, `exactskel/exactmom.py` (+ `mack.py`). Sources re-read:
CLS-D 2311.04271 `draft.tex` (eq. `light_transform_definition` and its embedding form `∫dα 𝒪(Z − αX, −X)`), CHMRS 2106.10274
`sec_heavy.tex` (eqs. `BMellinDefinition`, `Bkt action`, `CfromBMellin`) and `sec_light.tex` (eq. `eftformula`) (scratchpad copies).
Compute: task 1 4 min, task 2 10 min (+ 2 × 4 min diagnostics), 1 core each, < 0.5 GB.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and numerically verified; **[NE]** numerical evidence; **[M]** model; **[I]** open.

---

## 0. Summary

**Task 1 (free Dirac fermion, d = 3, exact connected `⟨TTTT⟩` at O(1/C_T)).**

1. **The spinor T-shift rule is derived, not tuned [D], and validated [DV] (§2).**
   * Derivation. The patch-1 field `ψ(x⁺)` obeys the free Dirac equation and patch Lorentz covariance, so every link is a constant times
     the flat formula `x̸(σx²)^{−3/2}`. The constants are then fixed by CLS-D's embedding definition `T(x⁺) = T_emb(−X)`, applied to the
     tensor 2-pt and 3-pt functions. The result is the same rule as the previous agent's scalar rule:
     * mixed links: positive base `(−x²)^{−p}`, flat numerator `x̸`, coefficient +1;
     * same-type links (12), (34): the patch formula with the Wightman iε.
   * Validation:
     * `⟨T₃T₁⁺⟩` equals the CLS-D continuation `H²/|X₁₃|^{Δ+2}` to 13 digits;
     * `⟨T₁⁺T₂⁺T₃⟩` and `⟨T₁⁺T₃T₄⟩` equal the CLS-D continuation of the conformal 3-pt function to ≤ 10⁻¹¹, also for generic Δψ;
     * the alternative branch κ = i ("no T-shift") gives exactly −1 there;
     * the MFT s- and t-terms reproduce v4c's exact edge terms (6·10⁻¹⁵ – 3·10⁻¹⁰).
2. **Lemma A holds pointwise to float64 precision [DV] (§3).**
   * **All 25 (2,4) tensor components, k = 2, 3, 4, 5, three angle pairs:** residual ≤ 3.8·10⁻¹⁵ of the largest family. The wrong
     symmetrisation gives 2.00.
   * **Per ν (ν outside, j′ = 0 kernel, ν = 1, 3):** ≤ 1.4·10⁻¹⁵ of a single family.
   * **Generic exponents.** The generalised-free-fermion `[ψ̄ψ]₀,₂` family (Δψ = 0.8–2.5) gives ≤ 1·10⁻¹¹, and its Δψ-derivative (logs in
     every link) gives 1·10⁻¹⁰.
3. **A structural difference from the scalar, which limits what the fermion adds (§3.1).**
   * In the light-ray gauge, the Wick classes 1243 and 1324 vanish **identically** for the fermion. Both contain the (13) propagator
     sandwiched between the two light-ray vertices, and `z̸₃(αz̸₁ + βz̸₃)z̸₁ = 0`.
   * So the whole fermion test sits in class 1234 and cancels as `1234(a,b) ± 1234(b,a) = 0`. It is non-trivial: families are O(S),
     residual 10⁻¹⁵.
   * Consequences:
     * every rotation-invariant structure (EE included) vanishes class by class, so the EE `Ψ₄` light side of the fermion is 0 at every
       angle;
     * the fermion does **not** exercise the (13) link, the twist-1 (13) channel or the collinear corner. Its families go like `Δφ⁻¹`,
       which becomes `Δφ^{+3/2}` with the kernel, so there is no corner singularity at all;
     * for the t₄ = 0 boson + Majorana mixture (the Einstein `⟨TTT⟩`), the (13)-channel content in this gauge is the scalar's alone.
4. **Controls (§3.5).** Each of these fails at O(1):
   * light-transform weight `t⁴`, `t⁵`: 0.04–1.3;
   * non-primary vertex: 0.13–1.8;
   * reflected spinor frame at the ⁺ points: 0.57–0.97;
   * wrong symmetrisation: 2.00.

   **The "no T-shift" control is degenerate for the fermion.** It is a constant phase per mixed link, which becomes an overall −1 on the
   single surviving class, so it gives 2·10⁻¹⁵. The branch rule is therefore pinned by the 2-pt/3-pt validation (κ = i gives −1 there),
   not by the Lemma-A test.
5. **Regrading (§5).**
   * Lemma A + D_LR stays **[D]**.
   * Its mechanism is now **[DV] on two independent exact spinning O(1/C_T) correlators** (free scalar and free fermion), and hence on
     every boson + fermion mixture, including the t₄ = 0 one, whose `⟨TTT⟩` is Einstein's.
   * **The Einstein instance is still [I].** The (13)/corner part of D_LR still rests on the scalar check alone.

**Task 2 (H_EFT, scalar-Mellin toy, single heavy exchange, Df = 3, k = 4, pure-moment C₄,ν at ν = 0 = the de-shifted value functional
of the scalar proxy).**

6. **The contour identity, and Lemma A in Mellin form [DV, toy] (§4.3).**
   * At fixed T, the residues of `M_h × kernel` split exactly into the heavy block `H` (s- and t-poles), the **actual double-trace part**
     `D` (kernel poles at double-twist locations, OPE side) and the arc at infinity. `H + D = Arc` is the exchange's actual total, by
     CHMRS's contour definition.
   * For J_h = 0, 2: `D = −H` (Arc = 0) at Δ_h = 20, 50, 100 to the numerical floor (10⁻¹⁰ … 1.3·10⁻⁴ at NC = 40, ν = 0, 2, 5; 4·10⁻¹² …
     1.4·10⁻⁶ at NC = 64, ν = 0), curvature corrections included. This is the J_h ≤ 3 control the OG check specified.
   * For J_h = 4, 6: `D ≠ −H`. `Arc/H = −30.3, −34.5, −34.9` (J_h = 4) and `−2477, −3308, −3396` (J_h = 6) at Δ_h = 20, 50, 100. These
     tend to the flat forward values −35 and −3423.
7. **Term-wise (P) of the double-trace part reproduces the actual double-trace value for every J_h, including J_h ≥ 4 [DV, toy] (§4.3).**
   * The order-R low-energy expansions `TP_R` converge to D, with best agreement 5·10⁻¹³ – 1.5·10⁻⁶ of the heavy image. Beyond about 10⁻⁸
     this is the numerical floor.
   * At Δ_h = 20, J_h = 6 the series is visibly asymptotic, with its optimum at R ≈ 40 ≈ πτ_h.
   * The mismatch is ≥ 10⁵× below the bound's sensitivity (δ = 0.2 margin +0.143; H_smear tolerance 24 % of B₂).
8. **Reading (§4.5).**
   * "Heavy block + term-wise double traces = actual total" holds for all J_h (Cauchy). But for J_h ≥ 4 that total is a non-zero arc, 35×
     (J_h = 4) to 3.4·10³× (J_h = 6) the heavy image. So the J_h ≤ 3 shortcut "double traces = −heavy block" must not be extended.
   * Result Q does not use it. It uses the full sum rule and evaluates the double traces term-wise, which the toy confirms exchange by
     exchange.
   * **H_EFT therefore splits** into:
     * (a) term-wise (P) of heavy-exchange double traces = their OPE-side value: **[DV] in the scalar-Mellin toy, any J_h**; spinning or
       position-space **[D-sketch]**;
     * (b) the heavy-sum / contact-tower interchange ("no conspiracy"): **[A], unchanged**.
   * The OG check's statement "H_EFT is a genuine assumption only for J_h ≥ 4" is refined accordingly.
9. **De-shift weights [NE] (§4.4).** The stored sc3 χ has exact moments `m₀…₃ = δ_{n0}` but `m₄ = −1.4·10⁻²`, `m₅ = −27`, … (ratio
   ≈ 1.9·10³ ≈ ν_max², a grid-cutoff effect). This is harmless for EFT terms suppressed by `Λ ≥ 300` (ratio ≤ 0.02), and should be stated
   with hypothesis 5.

---

## 1. What was and was not evaluated

| item | evaluated? | how |
|---|---|---|
| exact connected `⟨TTTT⟩` of the free massless Dirac fermion, d = 3 (Majorana = ½ Dirac) | **yes** | `of_fermion.py`: six oriented 4-cycles (classes 1234, 1243, 1324), `−tr[…]` per loop, exact |
| spinor T-shift and branch rule from CLS-D's embedding definition | **derived** (§2.2) | covariance + Dirac equation in patch 1, constants fixed by the tensor 2-pt and 3-pt functions |
| its validation on 2-pt, 3-pt and fermion-MFT pieces, incl. generic Δψ | **yes** | §2.3 |
| swapped form of `Ψ⁺_k`, k = 2–5, all 25 (2,4) tensor components, D_LR domain | **yes** | contour C of og_free, pointwise in angle |
| EE (and the two other rotation-invariant structures) | **yes** | 0 class by class (§3.2) |
| per-ν light side at ν = 1, 3, ν outside, j′ = 0 kernel | **yes** (generic component; EE is trivially 0) | §3.4 |
| controls: wrong weight, non-primary vertex, reflected frame, no T-shift, wrong symmetry | **yes** | §3.5 |
| generic exponents + logs (GFF `[ψ̄ψ]₀,₂`, Δψ = 0.8–2.5, ∂/∂Δψ) | **yes** | §3.6 |
| t₄ = 0 boson + Majorana mixture | **yes** (linear; numbers only) | §3.7 |
| Einstein `⟨TTTT⟩` | **no** | unchanged from RESULT_Q_OG_CHECK §5.1 |
| H_EFT scalar-Mellin toy, J_h = 0, 2 (control) and 4, 6, Δ_h = 20, 50, 100 | **yes** | `of_heft.py`, §4 |
| spinning H_EFT; position-space (F)+(P) continuation for J₀ ≥ 4 | **no** | the toy is Mellin-space (CHMRS's definition of the functional on M) |

---

## 2. Task 1: correlator, T-shift rule, validation

### 2.1 Correlator [D, exact]
* **Conventions.** Mostly-plus metric, real 2×2 gammas `γ⁰ = iσ²`, `γ¹ = σ¹`, `γ² = σ³`.
* **Propagator.** `S(x) = ⟨ψ(x)ψ̄(0)⟩ = x̸(x²)^{−p}`, with p = Δψ + ½ = 3/2. This is the analytic function; the Wightman functions are its
  boundary values. Its constant drops out because everything is normalised by `⟨TT⟩`.
* **Stress tensor.** `T(x,z) = ψ̄ z̸ (z·∂ψ) − (z·∂ψ̄) z̸ ψ`, z null. It is primary as it stands:
  `⟨TT⟩/[((Iz)·w)²/(x²)³] = 12.000000000000` at 6 configurations, spread 3·10⁻¹⁴. Hence `C2 = 144`, and every number below is in v4c
  `C_T = 1` MFT units.
* **Wick expansion.**
  * Connected part: the six oriented 4-cycles, `−tr[Γ S Γ S Γ S Γ S]` each, with derivative flags from the two vertex terms.
  * MFT: products of two 2-loops.
* **Links.** `d_i^m d_j^n [x̸ w^{−p}]` for m, n ≤ 1, in closed form (`link_mats`). They are checked against central finite differences to 10⁻⁷, the FD accuracy (development check, not in `of_run.py`).
* **Gauge.** Exactly v4c / og_free: `x₂ = −e`, `x₄ = e`, `x₁ = −t₁z₁` (patch point of 1⁺), `x₃ = t₃z₃`. The light-ray polarisations are
  `z₁`, `z₃`, the weight is `t₁³t₃³`, and the subtraction factor is `f_k`.
* **Exact (13) invariants.** `w₁₃ = 2y t₁t₃` and `x̸₁₃ = −t₁z̸₁ − t₃z̸₃`.

### 2.2 The spinor T-shift and branch rule [D]
CLS-D define the light transform in embedding space as `∫dα 𝒪(Z − αX, −X)`. Its endpoint operator is the **tensor** `T(x⁺) = T_emb(−X)`.
For J = 2 the sign of the polarisation is irrelevant.

**Why the spinor case needs care.** For spinors no canonical lift exists:
* in d = 3, `X → −X` is not in SO(3,2) (det = −1), and `x ↦ x⁺` reverses spatial orientation on R × S²;
* a twistor `S` with `S∧S = X` goes to `λS` with `λ² = −1`, or to a real det −1 frame change.

The rule is therefore fixed through T only, in three steps.

1. **Form of the links.**
   * Translations, dilatations and patch Lorentz transformations commute with `x ↦ x⁺`.
   * The free fermion in patch 1, in a frame where the patch Lorentz group acts in the standard way, obeys the same Dirac equation.
   * So each link is `c·x̸_ij |x_ij²|^{−p}` on its domain, with a constant c. Schur's lemma excludes any other matrix.
   * The constants are:
     * `c` for `⟨ψ(x⁺)ψ̄(x)⟩`;
     * `c′` for `⟨ψ(x)ψ̄(x⁺)⟩`;
     * `d` for the ⁺⁺ link (12);
     * `η₊` for the normalisation of the ⁺-vertex.
   * The invariants under field rescaling are `η₊cc′` and `η₊d`.
   * A reflected frame instead of the standard one would turn `T(x⁺,z)` into `T(x⁺,Pz)`. That is not covariant; it is control (f).
2. **2-pt function.**
   * CLS-D: `⟨T(x₃)T(x₁⁺)⟩ = H₁₃²/X₁₃^{Δ+2}` evaluated at `X₁ → −X₁`. `H₁₃` changes sign, `H²` does not, and the base is positive. So the
     value is `N_flat/|x₁₃²|^{Δ+2}`.
   * Wick terms have base power `P = 2p + a + b`, where a, b count derivatives hitting a base. Write each term as `Q(x²)^{Δ+2−P}/(x²)^{Δ+2}`.
   * The CLS-D value is then `Σ(−1)^{a+b}Q/|x²|^P`. This equals the positive-base Wick sum, because each base derivative flips one sign.
   * Hence **`η₊cc′ = +1`**.
3. **3-pt function.**
   * A conformal 3-pt structure is (polynomial in `V_i`, `H_ij`) over `Π X_ij^{τ/2}`. Under `X₁, X₂ → −X₁, −X₂` the numerator picks up
     `(−1)^{J₁+J₂−2#H₁₂} = +1`.
   * So the CLS-D value is `N_flat / [(x₁₂²)^{τ/2}|x₁₃²|^{τ/2}|x₂₃²|^{τ/2}]`.
   * The same counting (`τ/2 − P_ij = 1 − a_ij` for every link) gives the positive-base Wick sum times `η₊²dcc′`. Hence **`η₊d = +1`**.
   * The same holds for one ⁺ point and for any Δψ.

**Result.** The fermion rule is the scalar rule of og_free with spinor numerators:
* mixed links `x̸(−x²)^{−p}`, coefficient +1;
* (12), (34) by the patch formula with iε;
* everything differentiated in patch coordinates.

Only `η₊cc′` and `η₊d` enter. In a 4-cycle the ⁺ set is entered and left equally often, so the connected classes carry
`(η₊d)(η₊cc′)` (1234, 1243) and `(η₊cc′)²` (1324).

### 2.3 Validation (`of_check.txt` §1)

| check | result |
|---|---|
| Dirac `⟨TT⟩` / v4c structure | 12.000000000000, spread 3·10⁻¹⁴ (primary; C2 = 144) |
| GFF `[ψ̄ψ]₀,₂` primary: `ψ̄z̸∂ψ + a ∂ψ̄z̸ψ` | conformal iff a = −1 at Δψ = 0.8, 1.2, 1.5, 2.0, 2.5 (spread ≤ 2·10⁻¹⁴; a = −0.9 / −0.5 / 0: 0.01 / 0.4–0.8 / 4.4) |
| `⟨T₃T₁⁺⟩` (x₁₃ past-timelike) with the rule vs CLS-D `H²/\|X₁₃\|^{Δ+2}` | 1.0000000000000 (Δψ = 1, 0.8, 1.7; 5 configurations each) |
| `⟨T₁⁺T₂⁺T₃⟩`, `⟨T₁⁺T₃T₄⟩` vs CLS-D continuation of the 3-pt function (generic null polarisations) | 1.000000000000 (≤ 1·10⁻¹¹; Δψ = 1, 0.8, 1.6; 3 configurations each). **Control κ = i: −1.000** |
| MFT s, t (EE, k = 4, n = 100) vs v4d exact `∓RtEE` | s: 1·10⁻¹³ – 3·10⁻¹⁰; t: 6·10⁻¹⁵ – 8·10⁻¹³ |
| MFT u (exact 0, see §3.1: it contains the (13) link) | `\|u\|/\|t\|` ≤ 1·10⁻⁵ (roundoff amplified by the (13) base at small y) |
| MFT s, t, random null (2,4) polarisations, vs v4c `edge_residue_full` | t ≤ 3·10⁻¹³; s ≤ 5·10⁻⁶ (v4c's FFT accuracy) |
| Cauchy: circle radius 0.4 / 0.75 vs 0.5; n = 60 vs 100 (class 1234) | ≤ 1.3·10⁻¹⁵; 1.3·10⁻¹⁴ |

The 3-pt check must use generic polarisations at 1 and 3. With light-ray polarisations, the triangles that contain (13) vanish
identically (§3.1), and the comparison is 0/0.

---

## 3. Task 1 results

### 3.1 A kinematic identity: classes 1243 and 1324 vanish identically [P]
* In the light-ray gauge `x₁₃ = −t₁z₁ − t₃z₃ ∈ span(z₁, z₃)`.
* Every derivative of the (13) link is a combination of `z̸₁`, `z̸₃` and `x̸₁₃`. In both classes that link sits between `z̸₃` and `z̸₁`,
  and `z̸₃(αz̸₁ + βz̸₃)z̸₁ = 0` because `z̸² = z² = 0`.
* So 1243 = 1324 = 0 for all complex t₁, t₃ and all pol₂, pol₄. Numerically, at 40 random complex points `|1243|/|1234|` ≤ 6·10⁻¹¹
  (roundoff).
* The same identity kills the MFT u-term and every triangle that contains (13) with light-ray polarisations.
* **Consequences for the test.**
  * The fermion's Lemma-A content is `1234(a,b) ± 1234(b,a) = 0`. The links there are (12) ⁺⁺, (23), (41) mixed, (34).
  * The fermion does not probe the (13) link, its exact-invariant bookkeeping, the twist-1 (13) channel, or the collinear corner.

### 3.2 Rotation-invariant structures, k = 4

EE = `T₄·T₂`, `T₄₀₀T₂₀₀` and `T₄₀μT₂₀^μ` all satisfy `E(a,b) = E(b,a)`, so Lemma A predicts that every class vanishes. The reference
scale S is the largest `|1234|` over the 25 basis components at the same Δφ.

| Δφ | max over the three structures of \|1234\|/S | \|1243\|, \|1324\| | S (v4c units) |
|---|---|---|---|
| 10⁻³ | 3.0·10⁻¹³ | ≤ 1·10⁻⁵ S (roundoff; identically 0) | 1.5·10² |
| 10⁻² | 3.9·10⁻¹⁴ | ≤ 7·10⁻¹¹ S | 15 |
| 0.1 | 1.9·10⁻¹⁵ | ≤ 7·10⁻¹⁴ S | 1.5 |
| 0.5 – 3 | ≤ 6·10⁻¹⁶ | ≤ 3·10⁻¹⁵ S | 0.04 – 0.3 |

So the fermion's EE `Ψ₄` integrand is 0 at every angle. In the scalar the zero came from 1243 = −1324; here every class is zero on its own.

### 3.3 All (2,4) tensor components, k = 2–5 (25 null-basis pairs, n = 80)

Residual = `max |E(a,b) + (−1)^k E(b,a)| / S`, S = largest family.

| (φ₁, φ₃) | k = 2 | k = 3 | k = 4 | k = 5 | wrong symmetry (median) |
|---|---|---|---|---|---|
| (0.7, −1.1) | 1.9·10⁻¹⁵ | 1.2·10⁻¹⁵ | 1.8·10⁻¹⁵ | 2.4·10⁻¹⁵ | 2.00 |
| (2.6, 0.4) | 1.4·10⁻¹⁵ | 7.5·10⁻¹⁶ | 1.5·10⁻¹⁵ | 1.5·10⁻¹⁵ | 2.00 |
| (−2.2, 1.9) | 1.9·10⁻¹⁵ | 2.7·10⁻¹⁵ | 3.1·10⁻¹⁵ | 3.8·10⁻¹⁵ | 2.00 |

* For even k, 5 of the 25 components vanish identically (the same count as the scalar's k = 2 selection rule).
* The per-class residual of 1234 is ≤ 3·10⁻¹⁵.
* The gauge is complete (EE_FOUNDATIONS §1.4), so this covers every dual structure, including the j′ = 1 (`Ĉ₃`) structure of `Ψ₅`.

### 3.4 Collinear corner and per-ν light side
* **Corner [DV].** `Δφ·1234 → 3.25906·10⁻²` and `2.70222·10⁻³` (two components) at Δφ = 10⁻⁴…10⁻², stable to 6 digits. So each family goes
  like `Δφ⁻¹` (from `f₄`; 1234 has no (13) link). With the kernel `Δφ^{5/2}` the corner integrand is `Δφ^{+3/2}`.
  * The Lemma-A residual is still 2·10⁻¹¹ at Δφ = 10⁻⁴.
  * D_LR is therefore trivial for the fermion at O(1/C_T): no corner contribution of any prescription.
* **Per ν (ν outside, j′ = 0 kernel `2A₁₃^{(2+δ)/2}K_ν`, Φ = 0.3, Δφ ∈ [10⁻⁴, π], 64 nodes, generic component).** EE is identically 0
  (§3.2), so this uses a generic component.

| component | ν | I[1234(a,b)] | I[1234(b,a)] | \|I[sum]\|/\|I[1234(a,b)]\| |
|---|---|---|---|---|
| (0,2) | 1 | −0.536936 − 0.492452i | +0.536936 + 0.492452i | 1.3·10⁻¹⁵ |
| (0,2) | 3 | −0.270763 − 0.270720i | +0.270763 + 0.270720i | 8.8·10⁻¹⁶ |
| (3,1) | 1 | −0.846850 − 0.776690i | +0.846850 + 0.776690i | 1.4·10⁻¹⁵ |
| (3,1) | 3 | −0.254145 − 0.254104i | +0.254145 + 0.254104i | 1.0·10⁻¹⁵ |

### 3.5 Controls (random null (2,4) polarisations, two celestial pairs; residual relative to the largest class)

| control | what it breaks | k = 4 | k = 5 |
|---|---|---|---|
| none | — | 9·10⁻¹⁶, 3·10⁻¹⁴ | 3·10⁻¹⁵, 2·10⁻¹⁴ |
| (a) weight `t⁴` | light-transform covariance (Δ + J − 2 = 3) | 0.097, 0.93 | 1.1, 0.81 |
| (a) weight `t⁵` | same | 1.0, 1.3 | 0.043, 1.3 |
| (b) vertex a = 0 (`ψ̄z̸∂ψ`, = ½(T + ∂J)) | conformal primarity | 0.86, 0.33 | 0.91, 0.16 |
| (b) vertex a = −½ | same | 1.1, 1.8 | 1.4, 0.13 |
| (f) reflected frame at 1⁺, 2⁺ (`T(x⁺, Pz)`) | the spinor frame of §2.2 | 0.76, 0.97 | 0.59, 0.57 |
| (e) no T-shift (phase i per mixed link) | microcausality of Step 1 | 2·10⁻¹⁵, 4·10⁻¹⁴ — **degenerate** | 7·10⁻¹⁶, 5·10⁻¹⁵ — **degenerate** |
| wrong symmetry | — | 2.00 | 2.00 |

**(e) is uninformative for the fermion.** A constant phase per mixed link reweights the classes by (−1, −1, +1), and only 1234 is
non-zero. It is caught instead by the 2-pt/3-pt validation, where κ = i gives exactly −1 (§2.3). The "danger" flagged in
RESULT_Q_OG_CHECK §5.1 (a relative phase between mixed links failing silently) is excluded there, not by Lemma A.

### 3.6 Generic exponents and logs (GFF Dirac `[ψ̄ψ]₀,₂`, `Δ_O = 2Δψ + 1`, weight `t^{2Δψ+1}`, two random tensor/angle sets)

| Δψ | k = 4 residual | k = 5 residual | \|E\| range |
|---|---|---|---|
| 0.8 | 2·10⁻¹⁴, 3·10⁻¹⁵ | 1·10⁻¹⁴, 5·10⁻¹⁵ | 1.2–100 |
| 1.2 | 1·10⁻¹³, 1·10⁻¹⁴ | 1·10⁻¹⁴, 5·10⁻¹⁵ | 0.7–56 |
| 1.5 | 3·10⁻¹³, 1·10⁻¹⁵ | 2·10⁻¹³, 4·10⁻¹⁵ | 0.5–33 |
| 2.0 | 1·10⁻¹², 4·10⁻¹⁵ | 2·10⁻¹², 4·10⁻¹⁵ | 0.2–12 |
| 2.5 | 1·10⁻¹¹, 2·10⁻¹⁴ | 6·10⁻¹², 1·10⁻¹⁴ | 0.07–4.2 |

**∂/∂Δψ at Δψ = 1.5 (a log `x_ij²` in every link):**
* k = 4: `dE(a,b) = −60.909`, `dE(b,a) = +60.909`;
* k = 5: `+7.086`, `+7.086`;
* symmetrised residual 1·10⁻¹⁰ for both.

### 3.7 t₄ = 0 mixture
* In d = 3 the Majorana and real-scalar C_T are equal, and t₄ = 0 ⇔ n_s = n_M.
* Per field in C_T units, the scalar and Majorana connected families are of the same order, e.g. k = 5: −1.62 and −0.89. Their sum obeys
  Lemma A to ≤ 1.4·10⁻¹³.
* This follows from linearity. It is recorded because it is the correlator with Einstein's `⟨TTT⟩`.
* In the light-ray gauge its (13) channel is purely the scalar's (§3.1).

---

## 4. Task 2: H_EFT, scalar-Mellin toy check (`of_heft.py`)

### 4.1 Setting and definitions [D]
* **External operators and functional.** Identical scalars Df = 3 (= Δ_T), d = 3, k = 4.
  * Functional: CHMRS's Mellin sum rule `B̂_{4,T}[M] = ∮_{|S′|=∞} dS′/(2πi) M(S′,T′) K(S′)`, with
    `K = 16/[(2Df−S′)²(2Df+2−S′)(S′−T)(S′−T+2)]` and `T′ = T + 2Df − S′`.
  * Its pure moment: `C_{4,ν}[X] ∝ ∫_{Re T=2} dT Γ(T/2)² γ_{A−T}(ν) X(T)`, A = 8.5.
  * At ν = 0 this is the de-shifted value functional of the scalar proxy (STAGE4 §0 item 6: heavy image `P_{−1/2}`, polynomial light value
    `p(0)`).
* **Heavy exchange (three channels, residues = the block's Mack residues, mack.py conventions):**
  `M_h(S,T) = Σ_n [ρ_n(T)/(S − s_n) + ρ_n(S)/(T − s_n) + ρ_n(T)/(U − s_n)]`, with `s_n = τ_h + 2n` and `ρ_n = K_n Qpoly^n`.
* **Cauchy at fixed T.** The residues of `M_h K` split into:
  * `H(T)`: the s- and t-channel poles of M_h, i.e. `B̂_{4,T}`[heavy block];
  * `D(T)`: the kernel poles `S′ = 2Df` (double), `2Df + 2`, `T`, `T − 2`. This is the **double-trace part, OPE side, exact**.
  * Their sum `Arc := H + D` is the contour at infinity, i.e. the exchange's total contribution as CHMRS define `B̂` on M.
  * For J_h < 4, `M_h K ~ S′^{J_h−5}` gives `Arc = 0`. This is the Mellin form of Lemma A for a J₀ = J_h tree term.
* **Term-wise (P).**
  * `M^{(R)}` is the low-energy expansion of `M_h` to order R in `1/s_n`: `ρ_n(T)(1 − (S/s_n)^{R+1})/(S − s_n) + (t) + (u)`. It is a
    crossing-symmetric polynomial.
  * Its (P) light value is its kernel-pole sum (= its residue at ∞, CHMRS eq. `eftformula`).
  * `TP_R := C₄,ν[D[M^{(R)}]]`. "Term-wise continued light value of the double-trace part" = `lim TP_R`.
* **The caller's test**, "heavy block + term-wise double traces = actual total", reads `H + lim TP_R = Arc`. Two separate statements are
  checked:
  * **(i)** `lim TP_R = D`: term-wise (P) gives the actual OPE-side double-trace value;
  * **(ii)** whether `Arc = 0`.

### 4.2 Validation (`of_heft.txt` §0)

| check | result |
|---|---|
| `H(T)` from my residue sum vs `exactmom.Bhat_{4,T}` (block), 4 complex T | ratio −2.000000000000 (the −½ density factor) |
| n-sum convergence of `M_h`, `H` (N0 = 200 / 400 / 800) | 13 digits |
| kernel-circle radius 0.5 vs 0.8 | 13 digits |
| `C₄,₀[H]` on Re T = 2 vs exactmom's exact heavy image (ν → 0 fit), states (20,0), (20,2), (50,4), (30,6) | normalised ratio −1.6260083599 for all four (10 digits; state-independent) |
| fixed T: `Arc/H` for J_h = 0, 2 (Δ_h = 20, 50) | 3·10⁻¹³ – 7·10⁻⁹ (Arc = 0) |

* **Floor at J_h ≤ 3.** The residual of `(H + D)/H` comes from float64 cancellation among the four kernel-pole residues, each ≫ H ∝ m⁻⁸.
  * It is J-independent and falls with the circle resolution.
  * `of_heft_diag_nc40/64.txt`, NC = 40 → 64: 1.1·10⁻¹⁰ → 4·10⁻¹² (Δ_h = 20); 2.5·10⁻⁷ → 6·10⁻⁹ (50); 7.4·10⁻⁵ → 1.4·10⁻⁶ (100).
  * The T-integral itself cancels only mildly (κ ≈ 7–12).
  * The production run uses NC = 40.

### 4.3 Results (`of_heft.txt` §1, §2)

**Fixed T (no integral).**
* `Arc(T)/H(T)` for J_h = 0, 2 is ≤ 7·10⁻⁹ at T = 2 + 3i and 1 + 7i (Δ_h = 20, 50).
* For J_h = 4 and 6 it is O(10) and O(10²–10³): e.g. Δ_h = 50, J_h = 4 gives `|Arc/H|` = 26 and 22 at the two T.

**De-shifted value functional `C₄,₀` (pure moment, ν = 0), Re T = 2.**

`TP_R` is the term-wise (P) value of the order-R low-energy expansion. "best" is the minimum over R ∈ {0,…,60} of `|TP_R − D|/|H|`.

| Δ_h | J_h | D/H | Arc/H = (H + D)/H | \|TP_R − D\|/\|H\| at R = 3, 5, 8, 12, 20 | best (at R) |
|---|---|---|---|---|---|
| 20 | 0 | −1.0000000001 | −1.1·10⁻¹⁰ | 1.0, 0.32, 2.2·10⁻², 4.4·10⁻⁴, 1.9·10⁻⁷ | 1.0·10⁻¹² (40) |
| 20 | 2 | −1.0000000001 | −1.1·10⁻¹⁰ | 1.5, 2.8·10⁻², 1.9·10⁻², 8.5·10⁻⁴, 1.1·10⁻⁶ | 1.7·10⁻¹² (60) |
| 20 | 4 | −31.32685 | **−30.33** | 1.0, 0.72, 6.9·10⁻³, 2.2·10⁻³, 1.2·10⁻⁵ | 5.0·10⁻¹³ (40) |
| 20 | 6 | −2478.094 | **−2477.1** | 16, 1.2, 0.12, 2.6·10⁻³, 1.4·10⁻⁴ | 1.1·10⁻⁹ (40); 7.6·10⁻⁸ at R = 60 (asymptotic) |
| 50 | 0 | −1.00000025 | −2.5·10⁻⁷ (floor) | 1.0, 1.8·10⁻², 1.4·10⁻⁵, 1.6·10⁻⁸, 1.8·10⁻⁸ | 1.6·10⁻⁸ (12) |
| 50 | 2 | −1.00000026 | −2.6·10⁻⁷ (floor) | 0.34, 1.5·10⁻², 1.5·10⁻⁵, 1.7·10⁻⁸, 1.7·10⁻⁸ | 1.7·10⁻⁸ (12) |
| 50 | 4 | −35.54104 | **−34.54** | 0.98, 8.1·10⁻³, 1.7·10⁻⁵, 1.2·10⁻⁸, 7.2·10⁻⁹ | 7.2·10⁻⁹ (16) |
| 50 | 6 | −3309.437 | **−3308.4** | 2.4, 4.3·10⁻³, 2.0·10⁻⁵, 4.3·10⁻⁹, 7.3·10⁻⁹ | 4.3·10⁻⁹ (12) |
| 100 | 0 | −1.000074 | −7.4·10⁻⁵ (floor) | 1.0, 1.3·10⁻³, 8.9·10⁻⁸, 1.3·10⁻⁸, 1.3·10⁻⁸ | 1.3·10⁻⁸ (12) |
| 100 | 2 | −1.000073 | −7.3·10⁻⁵ (floor) | 0.83, 1.2·10⁻³, 4.1·10⁻⁷, 4.2·10⁻⁷, 4.2·10⁻⁷ | 4.1·10⁻⁷ (8) |
| 100 | 4 | −35.89424 | **−34.89** | 0.43, 1.1·10⁻³, 1.6·10⁻⁶, 1.5·10⁻⁶, 1.5·10⁻⁶ | 1.5·10⁻⁶ (12) |
| 100 | 6 | −3396.497 | **−3395.5** | 0.15, 9.3·10⁻⁴, 6.6·10⁻⁷, 5.5·10⁻⁷, 5.5·10⁻⁷ | 5.5·10⁻⁷ (12) |
| flat (Δ_h → ∞), forward | 4 / 6 | — | **−35 / −3423** | (`C₄,u=0`: arc = s⁴ coefficient of the u-channel `P_J(1 − 2s/m²)/m²` + s-, t-channel quotients; H = −2P_J(−1)/m¹⁰) | |

**Notes on the table.**
* For R ≤ 3 and J_h ≤ 2, `TP_R = 0` exactly (error 1.0). The k = 4 functional annihilates polynomials of degree < 4, so the light value
  comes entirely from the s⁴-and-higher terms of the expansion.
* Below the optimum the convergence is geometric. The fitted ratio is about 1/2.6 per order at τ_h = 20 and 1/10 at τ_h = 50, i.e.
  roughly `(6–8)/τ_h`. For Δ_h = 20, J_h = 6 the minimum at
  R ≈ 40 ≈ πτ_h followed by growth is the asymptotic character of the T-integrated series.
* At Δ_h ≥ 50 the "best" values are the numerical floor of §4.2, not the series.

**Other ν (`Arc/H`; tests the identity, not the term-wise sum).**
* ν = 2:
  * J_h = 0, 2: ≤ 8·10⁻⁵ (floor);
  * J_h = 4: −37.7, −35.5, −35.1 (Δ_h = 20, 50, 100);
  * J_h = 6: −4537, −3555, −3455.
* ν = 5:
  * J_h = 0, 2: ≤ 1.3·10⁻⁴;
  * J_h = 4: +680 (Δ_h = 20, where H(ν) is near a zero of `P_J(1 − 2ν²/m²)`-type), −41.7, −36.4;
  * J_h = 6: +5766, −5480, −3788.


### 4.4 The stored sc3 de-shift weights on polynomials (`of_heft_moments.txt`) [NE]
* **What a moment does.** Under the realised de-shifted value functional (`χ_value` of `functionals_q4_deshift.json`, sc3), a light-side
  polynomial `p(ν²) = Σ p_n ν^{2n}` is mapped to `Σ p_n m_n`, with `m_n = Σ_j ∫χR_j(ν − ij)^{2n}`.
* **Low moments.** `m₀ … m₃ = 1, −2·10⁻¹², −3·10⁻⁹, −7·10⁻⁶` (the exact δ_{n0} of V4FINAL §3).
* **High moments are not zero:** `m₄ = −1.4·10⁻²`, `m₅ = −27`, `m₆ = −5.2·10⁴`, …, ratio ≈ 1.9·10³ ≈ 1.16 ν_max² per order.
  * This is the grid cutoff (ν_max = 40.95), not a property of the continuum de-shift.
* **Effect.** It adds `Σ_{n≥4} p_n m_n` to the term-wise light value of a degree-2n EFT term.
  * For terms suppressed by `Λ^{−2n}` the series is geometric with ratio ≈ 1.9·10³/Λ²: 0.02 at Λ ≈ Δ_g = 300, 0.002 at 1000.
  * The first correction is `|m₄|/Λ⁸`, negligible there.
  * In the toy at Δ_h = 20–50 (m² ≲ 2.5·10³) the realised sc3 weights would not converge term-wise. This is why §4.3 uses the exact pure
    moment `C₄,₀`.
* **Side note.** The absolute moments for sc3 as stored are 1.49, 33, 3.5·10³, 6.6·10⁵ (n = 0–3). V4FINAL's 2.15, 82.6, … refer to a
  different normalisation or key.

### 4.5 Reading [D, toy-level DV]
1. **Term-wise (P) of an exchange's double-trace part equals its actual (OPE-side) double-trace contribution, for every J_h, including
   J_h ≥ 4.**
   * It converges to the numerical floor by order R ≈ 8–40.
   * The series is asymptotic in `1/τ_h`, with optimum `R* ≈ πτ_h` (visible at Δ_h = 20, J_h = 6).
   * This is the content of H_EFT that concerns heavy-exchange double traces. In the toy it holds with an error at the asymptotic
     optimum or the numerical floor: ≤ 1·10⁻⁹ of the heavy image at Δ_h = 20, ≤ 1.5·10⁻⁶ at Δ_h = 100.
   * The bound tolerates O(10⁻¹) (δ = 0.2 positivity margin +0.143, H_smear tolerance 24 % of B₂). The mismatch is ≥ 10⁵× below the
     bound's sensitivity.
2. **"Heavy block + term-wise double traces" reproduces the exchange's actual total contribution exactly.** This is Cauchy: Arc is
   *defined* as H + D. The non-trivial numerical content is item 1 (TP_R → D) and Arc = 0 for J_h ≤ 3.
   * For J_h ≤ 3 that total is 0 (Lemma A, Mellin form): `D = −H` to the floor, curvature corrections included, at Δ_h = 20–100. This is
     the control of RESULT_Q_OG_CHECK §5.2, now run.
   * For J_h ≥ 4 the total is the arc, `Arc/H ≈ −30 … −35` (J_h = 4) and `−2.5·10³ … −3.4·10³` (J_h = 6). It tends to the flat forward
     values −35 and −3423, with `O(1/Δ_h²)` corrections.
3. **What this means for H_EFT.**
   * The J_h ≤ 3 statement "double-trace value = −(heavy-block value)" must **not** be extended to J_h ≥ 4. It would be wrong by 35× and
     3.4·10³× the heavy image, far outside any tolerance.
   * Result Q does not use it. The bound uses the *full* correlator's sum rule (heavy side = all heavy blocks) and evaluates the
     double-trace parts term-wise. Item 1 shows that this evaluation is correct exchange by exchange.
   * The individual arcs never enter.
     * Formally, for a Regge-bounded completion they cancel against each other and against the explicit contacts:
       `Σ_h Arc_h + D[contacts] = 0`, if the heavy sum commutes with the contour at infinity. That commutation is doubtful, since the arcs
       grow fast with J_h.
     * The bound needs no such cancellation.
   * So within the scalar-Mellin toy, J_h ≥ 4 exchanges require nothing beyond what J_h ≤ 3 exchanges require.
   * What the toy does not test is the interchange of the sum over infinitely many heavy exchanges and the contact tower with the
     functional ("no conspiracy"). That stays hypothesis 1.
4. **Scope.**
   * Scalar externals only, Mellin-space definition of (P) (CHMRS), pure moments `C₄,ν` at ν = 0 (with ν = 2, 5 for the Arc identity)
     instead of the realised sc3 χ (§4.4).
   * Untested:
     * the spinning version;
     * the position-space (F)+(P) continuation that Result Q actually uses for J₀ ≥ 4 pieces;
     * the large-J_h regime, where the per-exchange expansion parameter is `J_h²/m²` rather than `1/m²`.

---

## 5. Regradings

### 5.1 Lemma A + D_LR

| claim | before (RESULT_Q_OG_CHECK) | after |
|---|---|---|
| pointwise Lemma A for J₀ < 4 tree terms (as used by Result Q) | [D] | **[D]** (unchanged) |
| proof mechanism (swapped form = 0 after 1↔3 symmetrisation; Lemmas B, C, T-shift, edges, t → 0 endpoint) | [DV] at O(1/C_T), free scalar + GFF scalar | **[DV] at O(1/C_T) on two independent exact spinning correlators**: free scalar and free Dirac/Majorana fermion (+ GFF fermion with generic exponents and logs). This covers every boson + fermion mixture, incl. the t₄ = 0 one with Einstein's `⟨TTT⟩`. All 25 components, k = 2–5, pointwise ≤ 4·10⁻¹⁵; per ν ≤ 1.4·10⁻¹⁵ |
| spinor T-shift / branch rule | not derived | **[D]** from CLS-D's embedding definition (via T's 2-pt and 3-pt functions); **[DV]**: 2-pt 13 digits, 3-pt ≤ 10⁻¹¹, κ = i gives −1 |
| D_LR, collinear corner at O(1/C_T) | [DV] (scalar: `Δφ^{−1/2}`) | unchanged; the fermion adds nothing: its corner integrand is `Δφ^{+3/2}`, and the (13) channel is identically 0 in this gauge |
| sensitivity of the zero | controls fail at O(1) (scalar) | fermion: weight, primarity, frame and symmetry controls fail at O(1); **the "no T-shift" control is degenerate** (only 1234 survives), and the branch is pinned by 2-/3-pt validation instead |
| Lemma A on Einstein `⟨TTTT⟩` (and E×W³) | [D], unchecked [I] | **unchanged [I]** |

**Recommended grade line:** "Lemma A + D_LR: [D]; mechanism [DV] at O(1/C_T) on the free scalar and the free fermion (RESULT_Q_OG_CHECK,
RESULT_Q_FERMION_HEFT), hence on the t₄ = 0 free mixture; (13)-channel/corner part checked on the scalar only; Einstein instance [I]."

**Why "closest proxy to Einstein" is weaker than hoped.**
* The fermion's t₄ = 0 property matters only through `⟨TTT⟩`.
* In the CLS-D light-ray gauge the fermion's T exchange in the (13) channel is invisible: the triangle with two light-ray vertices from a
  common x₀ vanishes.
* So for the mixture the (13) content is the scalar's, and Einstein's (13) channel (T + `[TT]` with anomalous dimensions) is still
  probed only by proxies.

### 5.2 H_EFT

| claim | before | after |
|---|---|---|
| H_EFT as one hypothesis (review issue 9) | [A]/[I] | **split into (a) and (b)** |
| (a) the term-wise (P) value of a heavy exchange's double-trace part equals its actual double-trace contribution | part of [A] | **[DV] in the scalar-Mellin toy** for J_h = 0, 2, 4, 6, Δ_h = 20–100, pure-moment ν = 0 functional: asymptotic, error ≤ 1.5·10⁻⁶ of the heavy image (floor). Holds equally for J_h ≥ 4. Spinning / position-space version: **[D-sketch]**, untested |
| (b) heavy-sum and contact-tower interchange ("no conspiracy") | [A] | **[A]** (unchanged; not tested) |
| "J_h ≤ 3 exchanges are fixed by Lemma A" (OG check §5.2) | [D] | **[DV] in the toy**: D = −H to ≤ 1.3·10⁻⁴ (NC = 40, ν = 0–5) / 1.4·10⁻⁶ (NC = 64, ν = 0), curvature corrections included |
| "H_EFT is a genuine assumption only for J_h ≥ 4" (OG check) | [D] | **refined.** J_h ≥ 4 exchanges are not fixed by Lemma A (`Arc/H ≈ −35, −3.4·10³`), but their double traces are still evaluated correctly term-wise. What remains assumed is (b), for all J_h |

---

## 6. What remains

1. **Einstein `⟨TTTT⟩` itself [I].** Unchanged; routes and costs as RESULT_Q_OG_CHECK §5.1 items 2–4. The fermion route (item 1 there) is
   now done.
2. **(13)-channel / collinear-corner part of D_LR on anything other than the free scalar.** The free fermion is blind to it in the CLS-D
   gauge. A GFF with a spin-1 or spin-2 single-trace exchange, or a 1/N vector model at its first non-trivial order, would be next.
3. **H_EFT (b)** (no conspiracy in the contact tower and the heavy sum) [A]. Also (a) for spinning externals and for the position-space
   (F)+(P) continuation.
4. **The realised sc3/EE3 de-shift weights** have `m_n ≠ 0` for n ≥ 4 (grid cutoff). This is harmless for Λ ≥ 300 (ratio ≤ 0.02 per order)
   but should be stated next to hypothesis 5.
5. **R⁴ normalisation control (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`)** in position space: still no J₀ = 4 spinning correlator. In the scalar-Mellin toy
   the polynomial light values are exact by construction.

---

## 7. Reproduce

```bash
cd ads4/code/resultq
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 of_run.py)            # task 1, all sections, ~4 min -> cache/of_check.txt
OMP_NUM_THREADS=1 python3 of_run.py val kin                           # single sections -> cache/of_check_val_kin.txt
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 of_heft.py)           # task 2, ~10 min -> cache/of_heft.txt
OMP_NUM_THREADS=1 python3 of_heft.py mom                              # sc3 moments -> cache/of_heft_moments.txt
OMP_NUM_THREADS=1 OF_NC=64 python3 of_heft.py diag                    # floor diagnostic -> cache/of_heft_diag_nc64.txt
```

| file | content |
|---|---|
| `of_fermion.py` | engine. Gammas, `slash`, closed-form link derivatives (`link_mats`), Wick loops (`loop`, `correlator`: classes 1234/1243/1324, MFT s/t/u), gauge links with the derived branch rule and exact (13) invariants (`links_gauge`), `Eplus(k, φ₁, φ₃, pol₂, pol₄, p, vert, q)`, 2-pt helpers |
| `of_run.py` | task 1 sections: `val` (§2.3), `kin` (§3.1), `ee` (§3.2), `tens` (§3.3), `corner`, `pair` (§3.4), `ctrl` (§3.5), `gff` (§3.6), `mix` (§3.7) |
| `of_heft.py` | task 2: `Exchange` (three-channel Mellin amplitude from exactmom/mack residues), `Hfast`, `Lkernel` / `Lchunk` (kernel-pole residues by circle trapezoid, `OF_NC` nodes), truncations `M^(R)`, ν-weights `tgrid`; sections (0)–(2), `mom` (§4.4), `diag` |
| `cache/of_check.txt`, `of_run.log` | recorded task-1 output (232 s) |
| `cache/of_heft.txt`, `of_heft.log`, `of_heft_moments.txt`, `of_heft_diag_nc40.txt`, `of_heft_diag_nc64.txt` | recorded task-2 output |

**Conventions to keep.**
* Divide fermion values by `C2 = 144` (Dirac, this engine's normalisation) to compare with v4c `C_T = 1` MFT numbers.
* The Majorana value is ½ Dirac.
* In the fermion gauge only class 1234 is non-zero. Residuals are relative to the largest family.
* In of_heft every ratio is relative to the heavy image H of the same exchange. The normalisation of M_h is arbitrary but common to
  H, D and TP.
