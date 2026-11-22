# Result Q, O(G)-order check of Lemma A + D_LR: the swapped-form `Ψ₄`, `Ψ₅` on an exact spinning `O(1/C_T)` `⟨TTTT⟩`

*22 Nov 2026, O(G)-check agent. Unrefereed, not committed. No existing document or code file was edited; no `cl_*` / `spinfunc` file was
touched. New code: `ads4/code/resultq/og_free.py` (engine), `og_gff.py` (generic-exponent family), `og_run.py` (production); output
`resultq/cache/og_check.txt` (+ `og_run.log`). Read-only imports: `v4c_core`, `v4d_core`, `v4c_wu.Kfun`. Source re-read: CLS-D
(2311.04271) `draft.tex`, eqs. `light_transform_definition`, `scalarfunctional_plus_final`, `spinningfunctional_definition` (scratchpad
copy). Compute: 2 min per full run, 1 core, < 1 GB.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and numerically verified; **[NE]** numerical evidence; **[M]** model; **[I]** open.

---

## 0. Summary

1. **No Einstein `⟨TTTT⟩` was evaluated.** The repo has no representation of the d = 3 Einstein correlator that the light-ray integrals can
   use (STAGE12 §2: "No position- or Mellin-space tree-level `⟨TTTT⟩` is available in a form the light-ray integrals can use").
   * A search of `ads4/code` for Einstein / graviton-exchange code finds only:
     * flat-space low sides (`omegau/u5_light.py`);
     * the KKSDZ shock matrix (`lightinputs/t2_shock_commutator.py`);
     * scalar CHMRS Mellin `k = 2` (`adslight/t1_scalar_mellin.py`);
     * three-point dictionaries;
     * heavy-block engines (`fourpt/tttt`, `e2e`). Writing one from scratch (momentum-space
   spinor-helicity results Fourier-transformed, or an O(G) `[TT]` data set in every spin structure) is not a few-hour job. **The
   Einstein-specific instance of Lemma A therefore stays unchecked [I].**
2. **Instead, the strongest available spinning check was done at the same order `O(1/C_T)`.** The correlator is the exact connected
   `⟨TTTT⟩` of the conformally coupled free scalar in d = 3. It is a unitary, microcausal, spinning `O(1/C_T)` correlator. It contains
   T exchange in all three channels (free-boson `⟨TTT⟩`), the twist-1 tower and half-integer branch cuts at every edge.
   * **Ψ₄ (k = 4) and Ψ₅ (k = 5)** are evaluated in the **swapped (double-commutator) form**, in the project's complete gauge, with
     light-ray polarisations and the CLS-D T-shift (D_LR domain).
   * **Lemma A holds pointwise to float64 precision [DV].** The test is `E⁺(a,b) ± E⁺(b,a) = 0` at fixed celestial angles, after the
     null-line integrals.
     * **EE structure:** the whole swapped-form integrand vanishes at every angle, `|E⁺|/|single family| ≈ 1·10⁻¹³` for Δφ from 10⁻⁴ to
       3. The individual families range from 3·10⁻⁶ up to 3·10⁹.
     * **All 25 (2,4) tensor components, k = 2, 3, 4, 5:** residual ≤ 8·10⁻¹⁵ of the largest family.
     * **Per-ν light side (ν outside, j′ = 0 kernel, EE):** at ν = 1 and ν = 3 the total is 1.1·10⁻¹³ and 5.1·10⁻¹³ of one family.
   * **Cancellation structure.** It is the `O(1/C_T)` analogue of MFT's `R_s(a,b) = −R_t(b,a)`:
     * Wick class `1234` vanishes on its own;
     * classes `1243` and `1324`, which 1↔3 exchanges, cancel each other.
3. **The zero is not automatic: the controls fail at O(1).** Each control keeps microcausality and Wightman analyticity:
   * a non-primary "T" (improvement ξ = 0 or ¼ instead of ⅛) gives residual **0.5–1.7**;
   * the wrong light-transform weight (`t⁴`, `t⁵` instead of `t³`) gives **0.4–1.8**;
   * dropping the T-shift (1 timelike to 3) gives **0.03–0.04**;
   * the wrong symmetry combination gives **2.00**.

   So the check tests the light-transform/vacuum-annihilation step of Lemma A (conformal covariance), not only microcausality.
4. **Generic exponents [DV].** The same holds (residual ≤ 3·10⁻¹¹) for the spin-2 double-trace primary `[φφ]₀,₂` of a generalised free
   field with Δφ = 0.6, 0.9, 1.3, 2.0. All edges then have generic, non-half-integer powers. The Δφ-derivative of the identity, which
   puts a `log x_ij²` into every link (the structure anomalous dimensions produce), vanishes to 10⁻¹¹.
5. **D_LR at O(1/C_T) [DV + D].**
   * Each family grows like `Δφ⁻³` at the collinear corner (coefficient 3.106·10⁻³, stable to 5 digits). With the kernel that gives
     `Δφ^{−1/2}`, exactly the twist-1 exponent of V4DPRIME's homogeneity lemma (`τ − 3/2`, τ = 1). It is integrable, so no contact term
     sits at φ₁ = φ₃.
   * The only prescription-sensitive piece remains the O(1) u-channel identity (already treated in V4FINAL §1).
6. **Regrading (§4).** Lemma A + D_LR stays **[D]**. Its proof chain is now **[DV] at O(1/C_T)** on exact spinning correlators, with
   sensitivity shown. **The Einstein instance is still [I].** The 10¹⁰–10²² leverage is not discharged; it is now only exposed to an
   Einstein-specific failure of an otherwise verified mechanism.
7. **R⁴ normalisation control (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`): not reachable.** No position-space J₀ ≥ 4 spinning correlator exists in the
   repo. The Regge probe below the free intercept (k = 1, 0) did not give a clean failure (§3.4).
8. **H_EFT (§5.2).**
   * For heavy exchanges with J_h ≤ 3, Lemma A fixes the double-trace light-side value at −(heavy-block value). So the term-wise (F)+(P)
     rule can be tested there.
   * A cheap scalar-Mellin toy check with existing ingredients is specified. It was not run.

---

## 1. What was and was not evaluated

| item | evaluated? | how |
|---|---|---|
| Einstein tree-level `⟨TTTT⟩` (any representation) | **no** | no usable representation exists in the repo (§5.1) |
| s-channel T block + `[TT]` O(G) data of Einstein | **no** | `γ_J`, `δp_J` of the spinning `[TT]_{n,J}` families in d = 3 are not available |
| exact `O(1/C_T)` spinning correlator: free-scalar connected `⟨TTTT⟩` | **yes** | `og_free.py`, all Wick cycles, exact |
| swapped (double-commutator) form of `Ψ⁺_k` at fixed angles, D_LR domain | **yes** | double contour integral = edge distributions (§2.3) |
| k = 4 (`Ψ₄`, j′ = 0) and k = 5 (`Ψ₅`, j′ = 1); also k = 2, 3 | **yes** | `f_k` of v4c_core / Stage 3 (Lemma B) |
| EE and every constant (2,4) tensor (25-component basis) | **yes** | null-basis decomposition |
| per-ν pairing with the j′ = 0 celestial kernel at ν = 1, 3 (ν outside) | **yes** (EE) | `v4c_wu.Kfun` quadrature |
| full de-shifted `Ĉ₄^{ds}` (χ-weighted ν-integral) | **implied, not computed** | the integrand is 0 pointwise (relative 10⁻¹³), so every ν-weight gives 0 |
| j′ = 1 celestial pairing for `Ψ₅` | implied | pointwise antisymmetric part = 0 at every angle |
| generic exponents / log structure | **yes** | GFF `[φφ]₀,₂` family, Δφ = 0.6 … 2.0, and ∂/∂Δφ |
| controls: non-primary T, wrong weight, no T-shift, wrong symmetry | **yes** | §3.4 |
| R⁴ contact (`2g₄`) | **no** | no position-space J₀ ≥ 4 spinning correlator |
| timelike-tube vs D_LR on the connected part | **not computed** | absolute corner integrability (§3.6) makes them equal by dominated convergence [D] |
| free fermion (would give Einstein's `⟨TTT⟩` as the t₄ = 0 boson + fermion mix) | **no** | spinor T-shift frame not derived (§5.1) |

**Why the free scalar is the right proxy.**
* Lemma A's proof (V4FINAL §2.2) uses only four inputs, and the free correlator has all of them:
  * microcausality at fixed angles;
  * Wightman analyticity, with f carrying the same iε (Lemma C);
  * conformal covariance of the light transform (vacuum annihilation);
  * endpoint convergence (J₀ < k).
* The free correlator is a genuine correlator of a unitary CFT with J₀ ≤ 1.
* Its connected part is the `O(1/N) = O(1/C_T)` term, the same order as Einstein's O(G).
* What it lacks relative to Einstein:
  * anomalous dimensions (logs at fixed external Δ = 3);
  * J₀ = 2 (it has J₀ ≤ 1);
  * the Einstein `⟨TTT⟩` mix (it has the pure-boson structure).

---

## 2. Representation and its validation

### 2.1 Correlator [D, exact]
* **Field and stress tensor.** Conformally coupled scalar, `⟨φφ⟩ = (x²)^{−1/2}`. The stress tensor (null z) is
  `T(x,z) = ¾(z·∂φ)² − ¼ φ(z·∂)²φ`, i.e. ξ = ⅛.
* **Connected 4-point function.** Sum over the 48 Wick matchings that form one 4-cycle, in three classes `1234`, `1243`, `1324` (16
  matchings each).
* **Disconnected part.** The 12 matchings with two 2-cycles (s, t, u). These give the MFT check.
* **Derivatives.** Directional derivatives `(a·∂)^m(−b·∂)^n (σX²)^{−p}` (m, n ≤ 2, null a, b) come from sympy closed forms.
* **Normalisation.** `⟨TT⟩_free = (3/2)×` the v4c `C_T = 1` structure at 5 random configurations, to 12 digits. So T is primary, and
  results are quoted divided by `(3/2)² = 2.25`.

### 2.2 Gauge and T-shift branches [D]
* **Gauge (v4c_core).**
  * `x₀ = 0`, `x₂ = −e` (patch point of 2⁺), `x₄ = e`;
  * `x₃ = t₃z₃`, `x₁ = −t₁z₁` (patch point of 1⁺);
  * polarisations z₁, z₃ (light-ray); measure `t₁³t₃³`; `f_k(u′,v′)` with patch `s_ij`.
* **Correlator.** `⟨T₄(x₄)T₃(x₃)T₁(x₁⁺)T₂(x₂⁺)⟩`, which is CLS-D's `scalarfunctional_plus_final` configuration.
* **Branches.** CLS-D define the operator on the light ray as `O(Z − αX₀, −X₀)` in embedding space (`light_transform_definition`), so
  `O(x⁺) = O_emb(−X_x)`.
  * Since `−X_x = Ω(cos(τ+π), sin(τ+π), −n)` with the same Ω > 0, `O(x⁺) = Ω^{−Δ}O_cyl(τ+π, −n)` carries **no phase**.
  * Every pair that links a ⁺ point to an unshifted point — (13), (14), (23), (24) — is spacelike on the cylinder.
  * Its Wightman function is therefore the positive real `(−x_ij²)^{−p}` (continued within the spacelike region from equal cylinder
    time).
  * The pairs (12) and (34) carry the ordinary patch formula `(x_ij²)^{−p}` with Wightman iε. "1 left of 2" is `t₁ + i0`; "4 left of 3"
    is `t₃ + i0`.
* **Consistency with the project.** For MFT this is v4c_core's rule ("−1 per mixed T link, naive patch formulas, sign +1 per Wick term")
  and reproduces its residues (§2.4).
* **The alternative assignment fails.** Giving each mixed φ-link a phase ±i (operator 1 at the patch point, not 1⁺) breaks Lemma A
  (§3.4 (e)).

### 2.3 Swapped form as a double contour integral [D]
* **One analytic function.** For `t₁, t₃ > ½` the pairs (12) and (34) are spacelike, so all four orderings of `⟨[T₄,T₃][T₁,T₂]⟩`
  coincide. They are the boundary values `t₁ ± i0`, `t₃ ± i0` of one analytic function `h = t₁³t₃³f_k G` (Lemma C: f's poles take the
  same iε).
* **The contour formula.** Hence

      E⁺(a,b) := ∫dt₁dt₃ t₁³t₃³ f_k ⟨[T₄,T₃][T₁,T₂]⟩ = ∮_C dt₁ ∮_C dt₃ h(t₁,t₃),   C: |t − ½| = ½, counter-clockwise from t = 0.

  This is the hairpin around the double-commutator support `(0, ½]` in each variable.
  * **Edges.** The edges t = ½ become the analytic continuation (Stage 3 §1.4, Gel'fand–Shilov).
  * **Rational h.** For rational h the formula reduces to v4c_core's `(−2πi)² Res Res`.
  * **Endpoint.** The only endpoint is t = 0 (`x₁⁺ → x₀`, the Regge end of the light transform). The integrand there is `O(t^{−1/2})`
    for k ≥ 2 (§3.4 (d)) and is made smooth by `θ = π(1 − cos ψ)`.
  * **Branch cuts.** On C, `(2t−1)^{−p}` is principal (cut = support). The (13) link is `(2y)^{−p}t₁^{−p}t₃^{−p}` with
    `y = 2sin²(Δφ/2)`. All other bases have positive real part.
* **Exact (13) invariants are required.** `x₁₃² = −2t₁t₃y`, `z₁·x₁₃ = t₃y`, `z₃·x₁₃ = t₁y` must be used in closed form.
  * Computing them from O(1) vectors loses `Δφ²` in float64.
  * A first run that did so showed spurious failures at small Δφ and a spurious k ≤ 1 "Regge failure". **Both disappeared with the exact
    invariants.**

### 2.4 Validation (`og_check.txt` §1, §2)

| check | result |
|---|---|
| `⟨TT⟩_free` / conformal structure, 5 random configurations | 1.5 exactly (12 digits): primary |
| MFT s-term (contour) vs exact `−RtEE(c)`, Δφ = 0.3, 1.3, 2.2, 3.0 | rel. 1·10⁻¹²–6·10⁻¹⁰ |
| MFT t-term (contour, two mixed links) vs exact `RtEE(c)` | rel. 2·10⁻¹⁵–2·10⁻¹³ |
| MFT u-term (light-ray `N ≡ 0`) | `|u|/|t|` ≤ 2·10⁻¹⁰ |
| random null (2,4) tensors: s, t vs v4c_core `edge_residue_full` (independent Cauchy-FFT code) | t: ≤ 8·10⁻⁹; s: ≤ 3·10⁻⁶ (v4c's own FFT accuracy for the 6th-order pole) |
| Cauchy independence: circle radius 0.4 / 0.5 / 0.75 (Δφ = 1) | family value stable to 12 digits; total ≤ 5·10⁻¹² |
| quadrature n = 60 vs 100 | family drift ≤ 5·10⁻¹³ (1·10⁻¹⁰ at Δφ = 3, where the families are smallest) |
| GFF `[φφ]₀,₂` primary: `a₀/a₁ = −Δφ/(Δφ+1)` (−⅓ at Δφ = ½, the improved T) | 2-pt conformality spread ≤ 1·10⁻¹² |

---

## 3. Results

### 3.1 EE structure, k = 4: the swapped-form integrand at fixed angle (units: v4c `C_T = 1` MFT)

Lemma A for EE (E⁺(a,b) = E⁺(b,a) by rotation + reflection) predicts **E⁺ ≡ 0**.

| Δφ | class 1234 | class 1243 | class 1324 | total | \|total\|/\|1243\| |
|---|---|---|---|---|---|
| 1·10⁻⁴ | 7·10⁻⁹ | +3.106088091681·10⁹ | −3.106088091681·10⁹ | 3·10⁻⁴ | 9·10⁻¹⁴ |
| 1·10⁻³ | 2·10⁻¹⁰ | +3.106085680580·10⁶ | −3.106085680580·10⁶ | 5·10⁻⁷ | 2·10⁻¹³ |
| 1·10⁻² | 2·10⁻¹² | +3.105844571251·10³ | −3.105844571250·10³ | 4·10⁻¹⁰ | 1·10⁻¹³ |
| 0.1 | 3·10⁻¹⁵ | +3.081741332083 | −3.081741332082 | 4·10⁻¹³ | 1·10⁻¹³ |
| 0.5 | 4·10⁻¹⁶ | +2.003541423826·10⁻² | −2.003541423826·10⁻² | 2·10⁻¹⁵ | 1·10⁻¹³ |
| 1.0 | 1·10⁻¹⁶ | +8.751676416542·10⁻⁴ | −8.751676416543·10⁻⁴ | 3·10⁻¹⁶ | 3·10⁻¹³ |
| 2.0 | 3·10⁻¹⁷ | −2.093323266576·10⁻⁴ | +2.093323266576·10⁻⁴ | 2·10⁻¹⁷ | 9·10⁻¹⁴ |
| 3.0 | 1·10⁻¹⁶ | −3.048441213115·10⁻⁶ | +3.048441213141·10⁻⁶ | 4·10⁻¹⁶ | 2·10⁻¹⁰ (quadrature drift 1·10⁻¹⁰) |

Size relative to MFT. A single connected family is `O(1)` times the MFT edge terms at the same angle (e.g. Δφ = 1.3: MFT `R_s = 1.9·10⁻⁴`;
connected 1243 `≈ 1.3·10⁻⁴` per free field). With N fields it is `1/N` of MFT, i.e. `O(1/C_T)`.

### 3.2 Every (2,4) tensor component, k = 2, 3, 4, 5 (connected part, 25 null-basis pairs, n = 80)

Residual `= max over components |E(a,b) ± E(b,a)| / S`, where S = the largest family value over all components. The sign is + for
even k and − for odd k.

| (φ₁, φ₃) | k = 2 | k = 3 | k = 4 | k = 5 | wrong-symmetry combination (median) |
|---|---|---|---|---|---|
| (0.7, −1.1) | 2.4·10⁻¹⁵ | 8.0·10⁻¹⁶ | 7.9·10⁻¹⁵ | 1.6·10⁻¹⁵ | 2.00 |
| (2.6, 0.4) | 1.5·10⁻¹⁵ | 3.8·10⁻¹⁵ | 8.3·10⁻¹⁵ | 3.1·10⁻¹⁵ | 2.00 |

* **The cancellation holds family by family.** `1234(a,b) ± 1234(b,a) = 0` and `1243(a,b) ± 1324(b,a) = 0`, each ≤ 6·10⁻¹⁵ S.
* **Some components vanish identically.** For k = 2, 5 of the 25 components vanish identically (a selection rule).
* **The gauge is complete** (EE_FOUNDATIONS §1.4), so this covers every dual structure, including the Ĉ₃ (j′ = 1) structure of `Ψ₅`.

### 3.3 Per-ν light side, swapped form, ν outside (EE, j′ = 0 kernel `A₁₃^{(2+δ)/2}K_ν`, Δφ ∈ [10⁻⁴, π], 64 nodes)

| ν | I[1234] | I[1243] | I[1324] | I[total] | \|I[total]\|/\|I[1243]\| | \|I[total]\|/∫\|ker\|\|1243\| |
|---|---|---|---|---|---|---|
| 1 | 1·10⁻¹⁵ | −1.339122·10⁻² − 1.228179·10⁻²i | +1.339122·10⁻² + 1.228179·10⁻²i | 2·10⁻¹⁵ | **1.1·10⁻¹³** | 7.5·10⁻¹⁴ |
| 3 | 7·10⁻¹⁶ | +1.996420·10⁻³ + 1.996098·10⁻³i | −1.996420·10⁻³ − 1.996098·10⁻³i | 1·10⁻¹⁵ | **5.1·10⁻¹³** | 9.3·10⁻¹⁴ |

* **Tail below 10⁻⁴.** The single-family tail below Δφ = 10⁻⁴ (`Δφ^{−1/2}` law, ~1 % of a family) is left out of both columns. The total's
  tail is zero, like its integrand.
* **De-shifted `Ĉ₄^{ds}`.** The integrand of the total is zero pointwise, so `∫χ(ν)N₀/N(ν) I_ν[total] dν = 0` for any weight. This holds in
  either ν-order and needs no corner convergence and no ν-decay, exactly as V4FINAL §2.2 asserts.
* **A single family need not have a finite de-shift.** For MFT, V4FINAL §2.5 found that such a de-shift diverges. It was not computed here.

### 3.4 Controls (connected part, random null (2,4) polarisations; residual relative to the largest family)

> **Qualification (RESULT_Q_EINSTEIN_CHECK §3.4).** The weight control is one-sided. Lowering the weight (q = −1, −2, −3) also gives 0 (10⁻⁶–10⁻¹⁵). Rows (a) therefore test the decay of the light-ray integral (weight ≤ Δ + J − 2), not exact light-transform covariance.

| control | what it breaks | k = 4 | k = 5 |
|---|---|---|---|
| none (ξ = ⅛, weight `t³`) | — | 1.4·10⁻¹⁴ | 9.2·10⁻¹⁵ |
| (a) weight `t⁴` | light-transform covariance (Δ + J − 2 = 3) | 0.39 | 8.6·10⁻³ |
| (a) weight `t⁵` | same | 1.8 | 1.8 |
| (b) ξ = 0 (non-primary T, still microcausal) | conformal covariance of T | 0.51 | 1.7 |
| (b) ξ = ¼ | same | 0.57 | 0.72 |
| (e) no T-shift (1 at x₁, timelike to 3) | microcausality of (13) in Step 1 | 0.039 | 0.025 |
| wrong symmetry `E(a,b) ∓ E(b,a)` | — | 2.00 | 2.00 |

**(c) Regge probe (neutral, not a control).** Fewer subtraction powers: `f₁ = (v′−u′)²` ("k = 1" = J₀ of free theory) and
`f₀ = (v′−u′)³` ("k = 0").
* k = 1: still zero (7·10⁻¹³, converged in n).
* k = 0: zero only to the quadrature noise (2·10⁻⁷–2·10⁻⁴, not converging in n).
* **(d) Endpoint exponents.** At t → 0 the exponents are `t^{−1/2}` (k ≥ 2), `t^{−3/2}` (k = 1) and `t^{−5/2}` (k = 0). For k ≤ 1 the
  iterated contour integral is therefore only conditionally convergent.
* **Conclusion.** The free correlator does not show a sharp failure at k = J₀. This probe cannot stand in for the R⁴ (J₀ = 4 = k) control.

### 3.5 Generic exponents and the log structure (GFF `[φφ]₀,₂`, Δ_O = 2Δφ + 2, weight `t^{2Δφ+2}`, two random tensor/angle sets)

| Δφ | k = 4 residual | k = 5 residual | \|E\| range |
|---|---|---|---|
| 0.6 | 2·10⁻¹³, 5·10⁻¹⁵ | 2·10⁻¹⁴, 2·10⁻¹³ | 0.06–18 |
| 0.9 | 3·10⁻¹⁴, 2·10⁻¹⁴ | 8·10⁻¹⁴, 5·10⁻¹⁴ | 1.8–66 |
| 1.3 | 1·10⁻¹², 2·10⁻¹⁴ | 9·10⁻¹³, 2·10⁻¹⁴ | 9–150 |
| 2.0 | 3·10⁻¹¹, 1·10⁻¹⁵ | 5·10⁻¹¹, 2·10⁻¹⁴ | 13–190 |

**∂/∂Δφ at Δφ = 0.9 (logs in every link).**

| k | dE(a,b)/dΔφ | dE(b,a)/dΔφ | relative residual of the symmetrised derivative |
|---|---|---|---|
| 4 | +202.17 | −202.17 | 1·10⁻¹¹ |
| 5 | −12.48 | −12.48 | 7·10⁻¹² |

The edge mechanism therefore works for generic powers `(x₁₂²)^{−p}` and for `log x₁₂²` terms, the forms Einstein's `[TT]` anomalous
dimensions produce. One caveat: here the external dimension moves with Δφ, whereas in Einstein the logs sit in internal data at
fixed Δ = 3.

### 3.6 Collinear corner and D_LR [DV + D]
* **Corner law [DV].** `E⁺_{1243}(Δφ)Δφ³ → 3.106088·10⁻³`, `3.106086·10⁻³`, `3.105845·10⁻³` at Δφ = 10⁻⁴, 10⁻³, 10⁻². So each family is
  `∝ Δφ⁻³`, and with the kernel `Δφ^{5/2}` the corner integrand is `∝ Δφ^{−1/2}`.
* **Agreement with V4DPRIME.** This is the homogeneity-lemma exponent `τ − 3/2` for a twist-1 (13) exchange: the free (13) channel
  has T and the twist-1 currents. **First exact spinning confirmation of that lemma.**
* **Integrability.** It is integrable with margin ½, per family.
* **Consequences for D_LR [D].**
  * At O(1/C_T) the collinear set carries nothing.
  * The total is 0 at every Δφ ≠ 0, so D_LR's value is 0.
  * Any regularisation that converges to the pointwise integrand (the tube, by dominated convergence) gives the same. The tube was not
    run numerically on the connected part.
* **Where prescription dependence lives.** It is confined to the O(1) u-channel identity, as in V4FINAL §1.4.

---

## 4. Regrading of Lemma A + D_LR

| claim | before | after |
|---|---|---|
| pointwise Lemma A for J₀ < 4 tree terms (statement, as used by Result Q) | [D] (MFT [P]); O(G) unchecked | **[D]**; no change for Einstein proper |
| the proof mechanism, swapped form = 0 after symmetrisation: Lemmas B, C, T-shift bookkeeping, edge continuation, t → 0 endpoint, relabelling | [D] | **[DV] at O(1/C_T)** on exact spinning correlators (free-scalar `⟨TTTT⟩_conn`; GFF `[φφ]₀,₂` family with generic exponents and logs). All (2,4) structures, k = 2–5, pointwise ≤ 10⁻¹⁴, per-ν ≤ 5·10⁻¹³. Sensitivity shown (§3.4). |
| Lemma A on the Einstein `⟨TTTT⟩` (and E×W³) | [D], unchecked | **[D], unchecked [I]**: no Einstein evaluation was made |
| D_LR: collinear set carries nothing at O(1/C_T) | [D] + [DV] (twist-1 model structures) | **[DV]**: exact spinning `Δφ^{−1/2}` corner law, total 0 pointwise |
| twist-1 homogeneity lemma (V4DPRIME) | [D] | **[DV]** on an exact spinning correlator |

**Reading for the unified statement.**
* The review's top request ("one swapped-form Ψ₄ on Einstein `⟨TTTT⟩` at one ν") is **not** met.
* What is met: one swapped-form `Ψ₄` (and `Ψ₅`) on an exact spinning `O(1/C_T)` `⟨TTTT⟩` at ν = 1 and ν = 3, plus pointwise in angle.
  It is zero to 10⁻¹³, and the controls fail at O(1).
* **Residual risk.** For the 10¹⁰–10²² leverage to be triggered, Einstein would have to violate one of the inputs that the controls
  show the zero depends on: conformal primary T, microcausality, J₀ < k. The swapped-form machinery itself would not be the cause.
  Tree-level Einstein satisfies all of these [D].
* **Recommended grade line.** "Lemma A + D_LR: [D]; mechanism [DV] at O(1/C_T) (RESULT_Q_OG_CHECK); Einstein instance [I]."

---

## 5. What remains

### 5.1 Toward the Einstein instance (in order of cost)
1. **Free fermion `⟨TTTT⟩_conn`.**
   * **Gain.** Linear combination with the boson at `t₄ = 0` makes the T-exchange content in every channel equal to Einstein's (up to
     normalisation).
   * **Needed.** The spinor version of §2.2 (frame of `ψ(x⁺)`).
   * **Danger.** A relative phase between mixed links would silently fail §3.4 (e). Validate on the fermion MFT edge terms first.
   * **Cost.** Hours.
2. **Einstein in the s-channel light-family form.** The ingredients are: the T block (Einstein `⟨TTT⟩`, light-ray polarised) plus
   `[TT]_{n,J}` with O(G) `γ`, `δp` for every spin structure in d = 3. The engine of §2.3 accepts any analytic h(t₁, t₃). The data do not
   exist (Lorentzian inversion of the crossed-channel T exchange in every helicity structure). **Cost: days.**
3. **Einstein from momentum space** (Raju 1201.6452 spinor-helicity `⟨TTTT⟩`). Fourier transforming to the light-ray configuration with
   complex t is not practical at 10⁻¹⁰ precision.
4. **R⁴ normalisation control.** It needs a J₀ = 4 spinning position-space correlator (AdS contact diagram). Not available.

### 5.2 H_EFT (review issue 9): concrete comment and a cheap check
* **Where Lemma A reaches.** A heavy exchange Witten diagram of spin J_h has J₀ = J_h.
  * **J_h ≤ 3.** The whole diagram is a J₀ < 4 tree correlator. Lemma A covers it (and §3 now verifies the mechanism at this order for
    generic exponents). So **its s-channel double-trace part has light-side value exactly −(its heavy-block value)**, and no new
    hypothesis is needed for it.
  * **J_h ≥ 4.** Neither the diagram nor its double-trace part has a convergent position-space functional. The term-wise (F)+(P) value
    is then a definition, so **H_EFT is a genuine assumption [A] only for J_h ≥ 4 exchanges and the contact tower.**
* **Cheap check for J_h ≤ 3 (not run).** It uses the scalar analogue that exists in the repo: the `sc3` entry of
  `functionals_q4_deshift.json` and the exact CHMRS Mellin moments in `exactskel/exactmom.py`. For a scalar (J_h = 0) and a spin-2
  (J_h = 2) heavy exchange with Δ_h ∈ {20, 50, 100}, compare:
  1. the exact heavy-block image under the de-shifted value functional (finite Δ_h, all orders in 1/m);
  2. minus the term-wise (P) value of the exchange's low-energy Mellin expansion. With V4FINAL §3's exact moments (value weight returns
     `p(0)` for degree ≤ 3), this is the forward low-energy coefficient: a closed-form sum over the exchange's Mellin poles.
* **What agreement would mean.** At J_h ≤ 3 the answer is fixed by Lemma A, so agreement validates the term-wise rule and the (F)+(P)
  normalisation at finite Δ_h (curvature corrections included). Repeating at J_h = 4, 6 then shows whether the term-wise value equals the
  analytic continuation in J_h of (1). That is the only sense in which an "actual contribution" exists for J₀ ≥ k.
* **Limits of the check.** In the flat limit, (1) = −(2) is the tautology "Taylor coefficient = −pole residue". The content is in the
  1/Δ_h² terms.
* **Cost.** A few hours (both ingredients exist).
* **What the check would not cover.** The convergence of the contact tower ("no conspiracy"). V4FINAL's absolute moments grow
  super-geometrically: 2.2, 83, 1.4·10⁴, 3.6·10⁶ for n = 0–3. So the term-wise series is at best asymptotic once
  `n ≳ log(Δ_g²)/log(ratio)`. That is a separate statement of hypothesis 1.

---

## 6. Reproduce

```bash
cd ads4/code/resultq
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 og_run.py)          # all sections, ~2 min, -> cache/og_check.txt
OMP_NUM_THREADS=1 python3 og_run.py val ee                          # single sections -> cache/og_check_val_ee.txt
```

| file | content |
|---|---|
| `og_free.py` | engine. Free-scalar T vertex; Wick matchings (classes 1234/1243/1324, MFT s/t/u); directional derivatives of `(σX²)^{−1/2}`; links with the T-shift branch rules and exact (13) invariants; `f_k`; contour C (cosine-substituted Gauss, optional graded panels); `Eplus(k, φ₁, φ₃, pol₂, pol₄)`; EE via 5-null-vector basis |
| `og_gff.py` | GFF `[φφ]₀,₂` family: general power p = Δφ, primary coefficient `a₀/a₁ = −Δφ/(Δφ+1)` with conformality check, weight `t^{2Δφ+2}` |
| `og_run.py` | sections `val` (§2.4), `ee` (§3.1), `tens` (§3.2), `ctrl` (§3.4), `pair` (§3.3), `gff` (§3.5) |
| `cache/og_check.txt`, `cache/og_run.log` | recorded output of the full run quoted here |

**Conventions to keep.**
* Divide by `C2 = 2.25` to compare with v4c `C_T = 1` MFT numbers.
* "Residual" is always relative to the largest single Wick-class value.
* Angles are `(φ₁, φ₃)` with `x₁ = −t₁z(φ₁)` and `x₃ = t₃z(φ₃)`.
